//! A choice popup: a question, numbered options with a line explaining the
//! selected one, chosen with the arrows and Enter or the option's key.

use crate::annotate::rounds::ClusterId;
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
    /// Show this member of the run's family; with `drop_edits`, even over
    /// unsaved edits.
    OpenRun {
        pick: crate::manifest::family::Pick,
        drop_edits: bool,
    },
    /// Label cluster `cluster` with this candidate.
    Label { cluster: ClusterId, label: String },
    /// Label the selected cluster on the panel tree.
    LabelOnTree,
    /// Label cluster `cluster` with a name typed by hand.
    LabelByName { cluster: ClusterId },
}

pub struct Choice {
    /// The key that takes it; none for options past the keys there are.
    pub key: Option<char>,
    pub label: String,
    /// One line on what choosing it does.
    pub detail: String,
    pub action: Action,
}

pub struct Menu {
    pub question: String,
    pub choices: Vec<Choice>,
    pub sel: usize,
    /// The status line when the menu is closed without a choice.
    pub on_cancel: String,
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
    pub fn new(
        question: impl Into<String>,
        choices: Vec<Choice>,
        on_cancel: impl Into<String>,
    ) -> Self {
        Self {
            question: question.into(),
            choices,
            sel: 0,
            on_cancel: on_cancel.into(),
            wide: false,
        }
    }

    /// ↑↓ (and the other list keys) move, Enter or an option's key takes it,
    /// esc or q cancels, as they close the form and the file browser.
    pub fn key(&mut self, code: KeyCode) -> Outcome {
        match code {
            KeyCode::Esc => Outcome::Cancel,
            KeyCode::Char(c) if CLOSES.contains(c) => Outcome::Cancel,
            KeyCode::Enter => self.take(self.sel),
            KeyCode::Char(c) => match self.choices.iter().position(|o| o.key == Some(c)) {
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

/// The keys besides esc that close a menu, which no option is keyed by.
const CLOSES: &str = "q";

/// `label`, keyed `1`, `2`, … in order (then `a`, `b`, … past nine, and
/// `A`, `B`, … past `z`).
pub fn numbered(items: Vec<(String, String, Action)>) -> Vec<Choice> {
    numbered_but(items, "")
}

/// As [`numbered`], skipping the keys in `taken`, which the menu's own
/// options have.
pub fn numbered_but(items: Vec<(String, String, Action)>, taken: &str) -> Vec<Choice> {
    let mut keys = ('1'..='9')
        .chain('a'..='z')
        .chain('A'..='Z')
        .filter(|c| !CLOSES.contains(*c) && !taken.contains(*c));
    items
        .into_iter()
        .map(|(label, detail, action)| Choice {
            key: keys.next(),
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
            "closed",
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
        let keys: String = numbered(items).iter().filter_map(|o| o.key).collect();
        assert_eq!(keys, "123456789ab");
    }

    #[test]
    fn q_cancels_and_no_option_is_keyed_by_a_closing_key() {
        let mut m = menu();
        assert_eq!(m.key(KeyCode::Char('q')), Outcome::Cancel);
        let items = |n| {
            (0..n)
                .map(|i| (format!("CT{i}"), String::new(), Action::StateOrder))
                .collect()
        };
        let keys: Vec<Option<char>> = numbered(items(70)).iter().map(|o| o.key).collect();
        assert!(!keys.contains(&Some('q')), "q closes the menu");
        let keyed: std::collections::BTreeSet<char> = keys.iter().flatten().copied().collect();
        assert_eq!(keyed.len(), keys.iter().flatten().count(), "no key twice");
        assert_eq!(keys[9..11], [Some('a'), Some('b')]);
        assert_eq!(keys[25], Some('r'), "after p comes r");
        assert_eq!(keys[34], Some('A'), "past z, capitals");
        assert_eq!(keys[69], None, "past the keys, only the arrows reach it");
        // Keys an option of the menu's own takes are skipped.
        let keys: String = numbered_but(items(30), "tn")
            .iter()
            .filter_map(|o| o.key)
            .collect();
        assert!(!keys.contains(['t', 'n', 'q']), "{keys}");
    }
}
