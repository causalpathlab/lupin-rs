//! A choice popup: a question, numbered options with a line explaining the
//! selected one, chosen with the arrows and Enter or the option's key.

use ratatui::crossterm::event::KeyCode;

/// What an option does once chosen.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Action {
    /// The prior orders nothing: pick a type to start from.
    PickRoot,
    /// Start from this type (`--root`).
    Root(String),
    /// Close the menu and state the order in the order table.
    StateOrder,
    /// Pick a precedence file (`--prior`).
    PriorFile,
    /// Pick a `label<TAB>CL:id` file (`--label-cl`).
    LabelCl,
    /// Annotate the run, then run the trajectory on its labels.
    Annotate,
    /// Pick a `cell<TAB>type` labels file.
    LabelsFile,
    /// After a pass: run the trajectory on the new labels.
    RunTrajectory,
    /// After a pass: leave the trajectory for later.
    NotNow,
    /// List the run's family: the run, its rounds, its trajectories (`g`).
    ListRuns,
    /// Open this member of the run's family.
    OpenRun(std::path::PathBuf),
    /// Open it, dropping the unsaved edits.
    OpenRunDropping(std::path::PathBuf),
    /// Keep what is on screen.
    Stay,
}

pub struct Choice {
    pub key: char,
    pub label: String,
    /// One line on what choosing it does.
    pub detail: String,
    pub action: Action,
}

pub struct Menu {
    pub question: String,
    pub choices: Vec<Choice>,
    pub sel: usize,
    /// The status line when the menu is closed without a choice; the
    /// trajectory's "nothing started" when unset.
    pub on_cancel: Option<String>,
    /// Drawn across the screen, not over the clusters column.
    pub wide: bool,
}

/// What a key did to the menu.
#[derive(Debug, PartialEq, Eq)]
pub enum Outcome {
    /// Still open (moved, or a key it does not take).
    Open,
    Cancel,
    Chosen(Action),
}

impl Menu {
    pub fn new(question: impl Into<String>, choices: Vec<Choice>) -> Self {
        Self {
            question: question.into(),
            choices,
            sel: 0,
            on_cancel: None,
            wide: false,
        }
    }

    /// ↑↓ (and the other list keys) move, Enter or an option's key takes it,
    /// esc cancels.
    pub fn key(&mut self, code: KeyCode) -> Outcome {
        match code {
            KeyCode::Esc => Outcome::Cancel,
            KeyCode::Enter => self.take(self.sel),
            KeyCode::Char(c) => match self.choices.iter().position(|o| o.key == c) {
                Some(i) => self.take(i),
                None => Outcome::Open,
            },
            _ => {
                super::app::step(&mut self.sel, self.choices.len(), code);
                Outcome::Open
            }
        }
    }

    fn take(&mut self, i: usize) -> Outcome {
        match self.choices.get(i) {
            Some(o) => {
                self.sel = i;
                Outcome::Chosen(o.action.clone())
            }
            None => Outcome::Open,
        }
    }
}

/// `label`, keyed `1`, `2`, … in order (then `a`, `b`, … past nine).
pub fn numbered(items: Vec<(String, String, Action)>) -> Vec<Choice> {
    items
        .into_iter()
        .enumerate()
        .map(|(i, (label, detail, action))| Choice {
            key: if i < 9 {
                char::from(b'1' + i as u8)
            } else {
                char::from(b'a' + (i - 9) as u8 % 26)
            },
            label,
            detail,
            action,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn menu() -> Menu {
        Menu::new(
            "which?",
            numbered(vec![
                ("one".into(), "the first".into(), Action::StateOrder),
                ("two".into(), "the second".into(), Action::PriorFile),
            ]),
        )
    }

    #[test]
    fn arrows_move_enter_or_a_key_takes_esc_cancels() {
        let mut m = menu();
        assert_eq!(m.key(KeyCode::Down), Outcome::Open);
        assert_eq!(m.sel, 1);
        assert_eq!(m.key(KeyCode::Enter), Outcome::Chosen(Action::PriorFile));
        assert_eq!(
            m.key(KeyCode::Char('1')),
            Outcome::Chosen(Action::StateOrder)
        );
        assert_eq!(m.key(KeyCode::Char('9')), Outcome::Open, "no such option");
        assert_eq!(m.key(KeyCode::Esc), Outcome::Cancel);
    }

    #[test]
    fn options_past_nine_are_keyed_by_letters() {
        let items = (0..11)
            .map(|i| (format!("CT{i}"), String::new(), Action::StateOrder))
            .collect();
        let keys: String = numbered(items).iter().map(|o| o.key).collect();
        assert_eq!(keys, "123456789ab");
    }
}
