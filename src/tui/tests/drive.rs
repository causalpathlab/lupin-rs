//! Driving an [`App`] in the TUI's tests: keys pressed, jobs ticked to their
//! end.

use super::app::App;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};

pub(super) fn press(app: &mut App, code: KeyCode) {
    app.key(KeyEvent::new(code, KeyModifiers::NONE));
}

/// Tick until the job is taken in.
pub(super) fn finish(app: &mut App) {
    for _ in 0..200 {
        app.tick();
        if app.child.is_none() {
            return;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    panic!("the job never ended");
}
