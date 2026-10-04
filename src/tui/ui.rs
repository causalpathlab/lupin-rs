//! Drawing the screen from an [`App`].

use super::app::{App, Focus, GeneView, TreeMode, CONTESTED, SETTINGS};
use super::ontology::{OntologyView, Role};
use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::celltype_tree::TreeSource;
use crate::annotate::markers::label_key;
use enrichment::UNASSIGNED_LABEL;
use ratatui::layout::{Constraint, Flex, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{
    Block, BorderType, Borders, Cell, Clear, Paragraph, Row, Table, TableState, Wrap,
};
use ratatui::Frame;

const BAR: usize = 8;

pub fn draw(f: &mut Frame, app: &App) {
    let [header, body, summary, log, keys] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(12),
        Constraint::Length(2),
        Constraint::Length(6),
        Constraint::Length(1),
    ])
    .areas(f.area());
    // Left: the clusters and their genes; middle: the ontology and the cell
    // types competing for the selected cluster; right, when the round scored
    // them: the cluster's GO terms.
    let (left, right, go) = if app.has_go() {
        let [l, m, g] = Layout::horizontal([
            Constraint::Percentage(30),
            Constraint::Percentage(35),
            Constraint::Percentage(35),
        ])
        .areas(body);
        (l, m, Some(g))
    } else {
        let [l, r] = Layout::horizontal([Constraint::Percentage(45), Constraint::Percentage(55)])
            .areas(body);
        (l, r, None)
    };
    let [clusters, genes] =
        Layout::vertical([Constraint::Min(8), Constraint::Length(16)]).areas(left);
    let [tree, candidates] =
        Layout::vertical([Constraint::Min(8), Constraint::Length(12)]).areas(right);

    draw_header(f, header, app);
    draw_clusters(f, clusters, app);
    draw_candidates(f, candidates, app);
    draw_genes(f, genes, app);
    draw_tree(f, tree, app);
    if let Some(go) = go {
        draw_go(f, go, app);
    }
    draw_summary(f, summary, app);
    draw_log(f, log, app);
    f.render_widget(Paragraph::new(help(app)).dim(), keys);
    if app.settings_open {
        draw_settings(f, app);
    }
    if app.help_open {
        draw_guide(f);
    }
    if app.prompt.is_some() {
        draw_prompt(f, app);
    }
}

/// The bottom line: where the key guide is, and the pane's main keys.
fn help(app: &App) -> &'static str {
    match app.focus {
        Focus::Clusters => {
            " ? keys · ↑↓ cluster · 1-6 take · k keep · ] next flagged · tab pane · s save · q quit"
        }
        Focus::Genes => {
            " ? keys · ↑↓ gene · a add · A add to a type · d drop · x hide · m view · tab pane"
        }
        Focus::Tree if matches!(app.tree_mode, TreeMode::Order(_)) => {
            if app.figures.as_ref().is_some_and(|v| v.shown) {
                " ? keys · v next figure · V table · t labels · c colour · m layout · , . components · p export · f exports · r run"
            } else {
                " ? keys · ↑↓ type · space mark · > precedes · - unrelated · r run · v figures · t tree"
            }
        }
        Focus::Tree => {
            " ? keys · ↑↓ node · enter label · space mark · + mixed label · o ontology · / search · t order"
        }
        Focus::Go => " ? keys · ↑↓ term · ← → read a long name · tab pane · esc clusters",
    }
}

/// Every key, by pane, in a popup (`?`).
const GUIDE: &[(&str, &[(&str, &str)])] = &[
    (
        "anywhere",
        &[
            ("?", "this guide (any key closes it)"),
            (
                "tab / shift-tab",
                "next / previous pane: clusters → genes → tree → GO terms (when scored)",
            ),
            ("pgup/dn home/end", "a page / to either end of the pane's list"),
            ("r", "cluster & run: Leiden and pass settings, enter runs"),
            (
                "x",
                "while a pass or save runs: stop it (asks again; in genes, x hides)",
            ),
            ("s", "save the edits as the next round, and export"),
            ("e", "export the open round"),
            ("q / ctrl-c", "quit (asks again with unsaved edits)"),
        ],
    ),
    (
        "clusters",
        &[
            ("↑↓", "select a cluster"),
            ("1-6", "label it with that candidate"),
            ("k", "keep its label (✓ = decided)"),
            ("u", "unassign it"),
            ("enter", "find its label in the tree"),
            ("⌫", "undo its edit"),
            ("]", "next flagged (?) cluster not yet decided"),
        ],
    ),
    (
        "tree",
        &[
            ("↑↓", "select a node"),
            ("enter", "label the cluster with it"),
            ("← →", "fold / unfold (ontology: up / into a term)"),
            ("o", "panel tree ↔ the full Cell Ontology"),
            (
                "/",
                "search the Cell Ontology (names, synonyms, abbreviations)",
            ),
            (
                "space / +",
                "mark nodes / give the cluster their mixed label (A+B)",
            ),
        ],
    ),
    (
        "genes",
        &[
            ("↑↓ / space", "select / mark genes"),
            ("m", "specific genes ↔ the label's markers"),
            ("a", "add as markers of the cluster's label"),
            (
                "A",
                "add as markers of any cell type (a new name makes a new type), then offer to label the cluster with it",
            ),
            ("d", "drop from the label's markers"),
            ("⌫", "undo the last marker edit"),
            ("x / X", "hide the genes / hide by pattern (MT-*)"),
            ("H", "show hidden genes"),
        ],
    ),
    (
        "GO terms",
        &[
            ("↑↓", "the selected cluster's top terms, by effect"),
            ("← →", "slide a long name back / on, a word at a time"),
            ("r", "GO terms in the settings: score them in the next pass"),
        ],
    ),
    (
        "order view (t in the tree pane)",
        &[
            ("↑↓", "a cell type on screen, with its cells"),
            ("space", "mark a type; mark two, first the earlier one (in the table, not the figures)"),
            (">", "the first marked type precedes the second (asks why)"),
            ("-", "the two marked types are unrelated (asks why)"),
            ("r", "run lupin trajectory on the labels on screen (asks for the output prefix)"),
            ("t / esc", "back to the tree / to the clusters"),
            ("", "statements go to the project's precedence.tsv; edges show source and, after a trajectory run, verdict"),
        ],
    ),
    (
        "figures (v in the order view, after a trajectory run)",
        &[
            ("v / V", "next figure: scatter, order, connectivity / back to the table"),
            ("t", "type labels on the scatter: small, medium, large, largest, off"),
            ("c", "colour the cells by pseudotime, cell type, lineage or component"),
            ("m", "the scatter's layout: each layout the run has (PHATE first), then the diffusion map"),
            (", .", "another pair of diffusion components"),
            ("p", "export the figure as {run}.trajectory.{figure}.svg + .pdf and log it"),
            ("f", "the exports strip: ↑↓ select, enter list an unlisted file, M move, d unlist, D D delete"),
            ("R", "re-read the log and check every export against its files"),
        ],
    ),
    (
        "marks",
        &[
            ("?  ✓", "cluster flagged (weak or no call) / decided"),
            ("new", "a specific gene no panel type lists"),
            (
                "weak",
                "a marker not enriched in this cluster (log2 FC ≤ 0)",
            ),
            ("◀", "the cluster's current label in the tree"),
        ],
    ),
];

fn draw_guide(f: &mut Frame) {
    let mut lines: Vec<Line> = Vec::new();
    for (section, keys) in GUIDE {
        if !lines.is_empty() {
            lines.push(Line::default());
        }
        lines.push(Line::from(format!(" {section}")).bold().cyan());
        for (key, what) in *keys {
            lines.push(Line::from(vec![
                Span::from(format!("   {key:<16}")).bold(),
                Span::from(*what),
            ]));
        }
    }
    let height = (lines.len() as u16 + 2).min(f.area().height);
    let [area] = Layout::vertical([Constraint::Length(height)])
        .flex(Flex::Center)
        .areas(f.area());
    let [area] = Layout::horizontal([Constraint::Length(84)])
        .flex(Flex::Center)
        .areas(area);
    f.render_widget(Clear, area);
    f.render_widget(
        Paragraph::new(lines)
            .style(POPUP)
            .block(popup(" keys ".into(), " any key closes ")),
        area,
    );
}

/// A popup's text: light on a dark ground, apart from the panes under it
/// whatever the terminal's theme.
const POPUP: Style = Style::new().fg(Color::Indexed(255)).bg(Color::Indexed(236));

/// A popup's frame: a bright border and title over [`POPUP`].
fn popup(title: String, hint: &'static str) -> Block<'static> {
    Block::default()
        .borders(Borders::ALL)
        .border_type(BorderType::Thick)
        .border_style(Style::new().fg(Color::LightYellow).bg(Color::Indexed(236)))
        .title(Line::from(title).bold())
        .title_bottom(Line::from(hint).fg(Color::Indexed(250)))
        .style(POPUP)
}

fn pane(title: String, focused: bool) -> Block<'static> {
    let b = Block::default().borders(Borders::ALL).title(title);
    if focused {
        b.border_style(Style::default().fg(Color::Cyan))
    } else {
        b
    }
}

fn highlight(focused: bool) -> Style {
    if focused {
        Style::default().add_modifier(Modifier::REVERSED)
    } else {
        Style::default().add_modifier(Modifier::BOLD)
    }
}

fn bar(share: f32) -> String {
    "█".repeat(((share.clamp(0.0, 1.0) * BAR as f32).round() as usize).max(1))
}

fn draw_header(f: &mut Frame, area: Rect, app: &App) {
    // The status first, where it is seen; the run's file last.
    let mut spans = vec![
        Span::from(" lupin annotate ").bold().reversed(),
        Span::from(format!(" {} ", app.status)).bold(),
    ];
    if !app.edits.is_empty() {
        spans.push(Span::from(format!("· {} unsaved ", app.edits.len())).yellow());
    }
    if app.stale {
        spans.push(Span::from("· settings changed, r to re-run ").yellow());
    }
    spans.push(Span::from(format!("· {}", app.source.display())).dim());
    f.render_widget(Paragraph::new(Line::from(spans)), area);
}

fn draw_clusters(f: &mut Frame, area: Rect, app: &App) {
    let focused = app.focus == Focus::Clusters;
    let Some(r) = &app.round else {
        f.render_widget(
            Paragraph::new(" no pass yet: r to run one").block(pane(" clusters ".into(), focused)),
            area,
        );
        return;
    };
    let rows: Vec<Row> = r
        .clusters
        .iter()
        .map(|c| {
            let label = r.label_of(c.id, &app.edits);
            let shown = label.clone().unwrap_or_else(|| UNASSIGNED_LABEL.into());
            let flag = if app.decided(c.id) {
                "✓"
            } else if c.flagged() {
                "?"
            } else {
                ""
            };
            let row = Row::new([
                format!("K{}", c.id),
                c.cells.to_string(),
                shown,
                format!("{:.2}", c.top_share()),
                flag.to_string(),
            ]);
            if app.edited(c.id) {
                row.yellow()
            } else if label.is_none() {
                row.dim()
            } else {
                row
            }
        })
        .collect();
    let header = Row::new(["", "cells", "label", "share", ""]).bold();
    let done = r.clusters.iter().filter(|c| app.edited(c.id)).count();
    let mut state = TableState::default().with_selected(Some(app.cluster_sel));
    let t = Table::new(
        rows,
        [
            Constraint::Length(5),
            Constraint::Length(6),
            Constraint::Min(20),
            Constraint::Length(5),
            Constraint::Length(1),
        ],
    )
    .header(header)
    .block(pane(
        format!(" clusters ({}, {done} edited) ", r.clusters.len()),
        focused,
    ))
    .row_highlight_style(highlight(focused));
    f.render_stateful_widget(t, area, &mut state);
}

/// The clustering and pass settings, as a popup over the screen.
fn draw_settings(f: &mut Frame, app: &App) {
    let [area] = Layout::vertical([Constraint::Length(SETTINGS.len() as u16 + 4)])
        .flex(Flex::Center)
        .areas(f.area());
    let [area] = Layout::horizontal([Constraint::Length(64)])
        .flex(Flex::Center)
        .areas(area);
    let rows: Vec<Row> = SETTINGS
        .iter()
        .map(|s| Row::new([s.name().to_string(), s.value(&app.args)]))
        .collect();
    let note = if app.fixed_clusters {
        "clusters from the manifest"
    } else {
        "Leiden on the cell embedding"
    };
    let mut state = TableState::default().with_selected(Some(app.setting));
    let t = Table::new(rows, [Constraint::Length(14), Constraint::Min(8)])
        .style(POPUP)
        .block(popup(
            format!(" cluster & run · {note} "),
            " ↑↓ setting · ←→ change · enter run (on output: edit) · esc close ",
        ))
        .row_highlight_style(highlight(true));
    f.render_widget(Clear, area);
    f.render_stateful_widget(t, area, &mut state);
}

/// `v` to `digits` places, `—` when unknown.
fn num(v: Option<f32>, digits: usize) -> String {
    v.map_or_else(|| "—".into(), |v| format!("{v:.digits$}"))
}

/// A p- or q-value: three places, scientific below 1e-3.
fn pval(v: Option<f32>) -> String {
    match v {
        Some(v) if v < 1e-3 => format!("{v:.0e}"),
        v => num(v, 3),
    }
}

/// The GO column's effect, p and q widths, before the term name.
const GO_STAT_WIDTHS: [u16; 3] = [6, 6, 6];

fn draw_candidates(f: &mut Frame, area: Rect, app: &App) {
    let Some(c) = app.selected() else {
        f.render_widget(
            Block::default().borders(Borders::ALL).title(" cluster "),
            area,
        );
        return;
    };
    let now = app.current_label();
    let mut rows: Vec<Row> = c
        .candidates
        .iter()
        .enumerate()
        .map(|(k, cand)| {
            let row = Row::new(vec![
                Cell::from(format!("{}", k + 1)).dim(),
                Cell::from(cand.label.clone()),
                Cell::from(format!("{:.2}", cand.share)),
                Cell::from(bar(cand.share)).fg(if cand.share >= CONTESTED {
                    Color::Green
                } else {
                    Color::Yellow
                }),
                Cell::from(num(cand.nes, 2)),
                Cell::from(pval(cand.p)),
                Cell::from(pval(cand.q)),
            ]);
            if now.as_deref() == Some(cand.label.as_str()) {
                row.bold()
            } else {
                row
            }
        })
        .collect();
    // A label picked off the candidates (a new type, or one that did not
    // pass FDR) is shown too: it is scored once the edits are saved.
    if let Some(l) = now.as_deref().filter(|l| {
        !c.candidates
            .iter()
            .any(|cand| label_key(&cand.label) == label_key(l))
    }) {
        let share = c
            .shares
            .iter()
            .find(|(t, _)| label_key(t) == label_key(l))
            .map(|(_, s)| *s);
        rows.insert(
            0,
            Row::new(vec![
                Cell::from("→"),
                Cell::from(Line::from(vec![
                    Span::from(l.to_string()).bold(),
                    Span::from(pinned_note(app, l, share)).dim(),
                ])),
                Cell::from(share.map_or_else(|| "—".into(), |s| format!("{s:.2}"))),
                Cell::from(share.map(bar).unwrap_or_default()),
            ]),
        );
    }
    if rows.is_empty() {
        rows.push(Row::new(["", "no candidate passed FDR"]).dim());
    }
    let path = now.as_deref().and_then(|l| app.tree.node_of(l)).map(|i| {
        app.tree
            .path(i)
            .iter()
            .map(|&n| app.tree.nodes[n].name.clone())
            .collect::<Vec<_>>()
            .join(" › ")
    });
    let scored = if app.rescoring.is_some() {
        " · rescoring…"
    } else if app.recorded.is_some() {
        " · rescored with your marker edits"
    } else {
        ""
    };
    let title = format!(
        " K{} · {} cells · {}{scored} ",
        c.id,
        c.cells,
        now.as_deref().unwrap_or(UNASSIGNED_LABEL)
    );
    let t = Table::new(
        rows,
        [
            Constraint::Length(1),
            Constraint::Min(14),
            Constraint::Length(5),
            Constraint::Length(BAR as u16),
            Constraint::Length(6),
            Constraint::Length(6),
            Constraint::Length(6),
        ],
    )
    .header(Row::new(["", "type", "share", "", "NES", "p", "q"]).bold())
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(title)
            .title_bottom(Line::from(format!(" {} ", path.unwrap_or_default())).cyan()),
    );
    f.render_widget(t, area);
}

/// Why a cluster's label is not among its candidates.
fn pinned_note(app: &App, label: &str, share: Option<f32>) -> &'static str {
    let edited = app.edits.iter().any(|e| {
        matches!(e, super::round::Edit::Markers { label: l, .. } if label_key(l) == label_key(label))
    });
    match share {
        Some(s) if s > 0.0 => " · below the top candidates",
        _ if app.rescoring.is_some() => " · rescoring…",
        _ if edited && app.recorded.is_none() => " · scored on save",
        _ => " · not called (did not pass FDR)",
    }
}

fn draw_genes(f: &mut Frame, area: Rect, app: &App) {
    let focused = app.focus == Focus::Genes;
    let (Some(c), Some(r)) = (app.selected(), &app.round) else {
        f.render_widget(pane(" genes ".into(), focused), area);
        return;
    };
    let now = app.current_label();
    let key = now.as_deref().map(label_key);
    let mark = |g: &str| {
        if app.marked.iter().any(|m| m == g) {
            "●"
        } else {
            " "
        }
    };
    // Only the rows that fit are built: the lists run to thousands of genes.
    let height = area.height.saturating_sub(2) as usize;
    let start = app.gene_sel.saturating_sub(height.saturating_sub(1));
    let window = start..start + height;
    let (rows, title): (Vec<Row>, String) = match app.gene_view {
        GeneView::Specific => {
            let fc: std::collections::HashMap<&str, f32> =
                c.genes.iter().map(|(g, s)| (g.as_str(), *s)).collect();
            let listed = app.listed_genes();
            let rows = listed[window.start.min(listed.len())..window.end.min(listed.len())]
                .iter()
                .map(|g| {
                    let of = r.marker_of(g);
                    let hidden = app.hidden.hides(g);
                    // A gene no panel type lists: the panel missed it.
                    let tag = if hidden {
                        "hidden".to_string()
                    } else if of.is_empty() {
                        "new".to_string()
                    } else {
                        of.join(",")
                    };
                    let row = Row::new([
                        mark(g).to_string(),
                        g.clone(),
                        format!("{:+.2}", fc.get(g.as_str()).copied().unwrap_or(f32::NAN)),
                        tag,
                    ]);
                    if hidden {
                        row.dim()
                    } else if key.as_deref().is_some_and(|l| of.contains(&l)) {
                        row.green()
                    } else if of.is_empty() {
                        row.yellow()
                    } else {
                        row.dim()
                    }
                })
                .collect();
            let shown = if app.show_hidden {
                "H: hide them"
            } else {
                "H: show hidden"
            };
            (
                rows,
                format!(
                    " {} specific genes · log2 FC · new = in no panel · {shown} ",
                    listed.len()
                ),
            )
        }
        GeneView::Markers => {
            let Some(label) = now.as_deref() else {
                f.render_widget(
                    Paragraph::new(" unassigned: no markers")
                        .block(pane(" markers ".into(), focused)),
                    area,
                );
                return;
            };
            let listed = app.label_markers();
            let rows = listed[window.start.min(listed.len())..window.end.min(listed.len())]
                .iter()
                .map(|(g, fc, added)| {
                    // A marker not enriched in this cluster: missed here.
                    let weak = fc.is_none_or(|v| v <= 0.0);
                    let tag = match (added, weak) {
                        (true, _) => "added",
                        (false, true) => "weak",
                        (false, false) => "",
                    };
                    let row = Row::new([
                        mark(g).to_string(),
                        g.clone(),
                        fc.map_or_else(|| "—".into(), |v| format!("{v:+.2}")),
                        tag.to_string(),
                    ]);
                    match (added, weak) {
                        (true, _) => row.yellow(),
                        (false, false) => row.green(),
                        (false, true) => row.dim(),
                    }
                })
                .collect();
            let weak = listed
                .iter()
                .filter(|(_, fc, _)| fc.is_none_or(|v| v <= 0.0))
                .count();
            (
                rows,
                format!(
                    " {} markers of {label} · {weak} weak (FC ≤ 0 in K{}) ",
                    listed.len(),
                    c.id
                ),
            )
        }
    };
    let mut state = TableState::default().with_selected(focused.then_some(app.gene_sel - start));
    let widths = [
        Constraint::Length(1),
        Constraint::Length(14),
        Constraint::Length(7),
        Constraint::Min(8),
    ];
    let t = Table::new(rows, widths)
        .block(pane(title, focused))
        .row_highlight_style(highlight(focused));
    f.render_stateful_widget(t, area, &mut state);
}

/// The selected cluster's top GO terms: each term's effect (its genes' mean
/// in the cluster against the other clusters) and where it came from.
fn draw_go(f: &mut Frame, area: Rect, app: &App) {
    let focused = app.focus == Focus::Go;
    let Some(c) = app.selected() else {
        f.render_widget(pane(" GO terms ".into(), focused), area);
        return;
    };
    // A term name is a phrase: wrapped, up to three lines, to stay readable
    // in a narrow column. The list is short (each cluster's top terms), so
    // the table scrolls it whole.
    // Borders, the stat columns and a space after each.
    let fixed: u16 = 2 + GO_STAT_WIDTHS.iter().map(|w| w + 1).sum::<u16>();
    let width = area.width.saturating_sub(fixed).max(8) as usize;
    let rows: Vec<Row> = c
        .terms
        .iter()
        .enumerate()
        .map(|(i, t)| {
            // The selected name starts `go_shift` words in (← →).
            let text = match app.go_shift {
                n if n > 0 && i == app.go_sel => {
                    let rest: Vec<&str> = t.term.split_whitespace().skip(n).collect();
                    format!("… {}", rest.join(" "))
                }
                _ => t.term.clone(),
            };
            let lines = wrap(&text, width, 3);
            let height = lines.len() as u16;
            Row::new([
                format!("{:+.2}", t.effect),
                pval(t.p),
                pval(t.q),
                lines.join("\n"),
            ])
            .height(height)
        })
        .collect();
    let title = if c.terms.is_empty() {
        format!(" GO terms · none for K{} ", c.id)
    } else {
        format!(" GO terms of K{} · effect · p · q ", c.id)
    };
    let mut state = TableState::default().with_selected(focused.then_some(app.go_sel));
    let widths = GO_STAT_WIDTHS
        .map(Constraint::Length)
        .into_iter()
        .chain([Constraint::Min(8)]);
    let t = Table::new(rows, widths)
        .block(pane(title, focused))
        .row_highlight_style(highlight(focused));
    f.render_stateful_widget(t, area, &mut state);
}

/// `text` in lines of at most `width` characters, broken between words, at
/// most `max` lines (the last ends in `…` when cut).
fn wrap(text: &str, width: usize, max: usize) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    for word in text.split_whitespace() {
        match lines.last_mut() {
            Some(l) if l.chars().count() + 1 + word.chars().count() <= width => {
                l.push(' ');
                l.push_str(word);
            }
            _ => lines.push(word.chars().take(width).collect()),
        }
    }
    if lines.len() > max {
        lines.truncate(max);
        let last = &mut lines[max - 1];
        let keep: String = last.chars().take(width.saturating_sub(1)).collect();
        *last = format!("{keep}…");
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    lines
}

fn draw_tree(f: &mut Frame, area: Rect, app: &App) {
    if let TreeMode::Order(v) = &app.tree_mode {
        if let Some(fig) = app.figures.as_ref().filter(|fig| fig.shown) {
            return draw_figures(f, area, app, fig);
        }
        return draw_order(f, area, app, v);
    }
    if let (TreeMode::Ontology(v), Some(cl)) = (&app.tree_mode, &app.cl) {
        return draw_ontology(f, area, app, v, cl);
    }
    let focused = app.focus == Focus::Tree;
    let c = app.selected();
    let now = app.current_label().map(|l| label_key(&l));
    let rows: Vec<Row> = app
        .tree
        .visible()
        .into_iter()
        .map(|i| {
            let n = &app.tree.nodes[i];
            let glyph = match (n.children.is_empty(), app.tree.is_folded(i)) {
                (true, _) => "·",
                (false, true) => "▸",
                (false, false) => "▾",
            };
            let under = app.tree.labels_under(i);
            let share: f32 = c.map_or(0.0, |c| {
                c.candidates
                    .iter()
                    .filter(|c| under.contains(&label_key(&c.label).as_str()))
                    .map(|c| c.share)
                    .sum()
            });
            let here = now.as_deref() == Some(app.tree.label(i));
            let marked = app.tree_marked.iter().any(|m| m == app.tree.label(i));
            let name = format!(
                "{}{}{glyph} {}{}",
                if marked { "●" } else { " " },
                "  ".repeat(n.depth),
                n.name,
                if here { "  ◀" } else { "" }
            );
            let row = Row::new([
                name,
                if share > 0.0 {
                    format!("{share:.2}")
                } else {
                    String::new()
                },
                n.cl_id.clone().unwrap_or_default(),
            ]);
            if here {
                row.bold().cyan()
            } else if share > 0.0 {
                row
            } else {
                row.dim()
            }
        })
        .collect();
    let title = match (&app.tree.source, &app.tree.release) {
        (TreeSource::CellOntology, r) => format!(
            " cell ontology · {} · evidence of the selected cluster ",
            r.as_deref().unwrap_or("release unknown")
        ),
        (TreeSource::MarkerSharing, _) => " cell types, grouped by shared markers ".into(),
    };
    let mut state = TableState::default().with_selected(Some(app.tree_sel));
    let t = Table::new(
        rows,
        [
            Constraint::Min(30),
            Constraint::Length(6),
            Constraint::Length(12),
        ],
    )
    .block(pane(title, focused))
    .row_highlight_style(highlight(focused));
    f.render_stateful_widget(t, area, &mut state);
}

/// The Cell Ontology around a term: its parents, itself, its children (or a
/// search's hits), each with the selected cluster's evidence under it and
/// how many panel types it covers.
fn draw_ontology(f: &mut Frame, area: Rect, app: &App, v: &OntologyView, cl: &ClTerms) {
    let focused = app.focus == Focus::Tree;
    // The selected cluster's candidates on their terms, with every term above.
    let evidence: Vec<(f32, std::collections::BTreeSet<String>)> = app
        .selected()
        .map(|c| {
            c.candidates
                .iter()
                .filter_map(|cand| {
                    let t = super::ontology::term_of(cl, &app.tree, &cand.label)?;
                    Some((cand.share, cl.ancestors_or_self(&t)))
                })
                .collect()
        })
        .unwrap_or_default();
    let now = app
        .current_label()
        .and_then(|l| super::ontology::term_of(cl, &app.tree, &l));
    let rows: Vec<Row> = v
        .rows
        .iter()
        .map(|r| {
            let glyph = match r.role {
                Role::Parent => "↑",
                Role::Focus => "●",
                Role::Hit => "⌕",
                Role::Child if cl.children(&r.id).is_empty() => "·",
                Role::Child => "▸",
            };
            let indent = if r.role == Role::Child { "  " } else { "" };
            let share: f32 = evidence
                .iter()
                .filter(|(_, above)| above.contains(&r.id))
                .map(|(s, _)| s)
                .sum();
            let types = app
                .panel_ancestry
                .iter()
                .filter(|(_, above)| above.contains(&r.id))
                .count();
            let here = now.as_deref() == Some(r.id.as_str());
            let marked = !app.tree_marked.is_empty()
                && app
                    .tree_marked
                    .contains(&super::ontology::term_label(cl, &app.tree, &r.id));
            let name = format!(
                "{}{indent}{glyph} {}{}",
                if marked { "●" } else { " " },
                cl.name(&r.id).unwrap_or(&r.id),
                if here { "  ◀" } else { "" }
            );
            let row = Row::new([
                name,
                if share > 0.0 {
                    format!("{share:.2}")
                } else {
                    String::new()
                },
                if types > 0 {
                    format!("{types} type{}", if types == 1 { "" } else { "s" })
                } else {
                    String::new()
                },
                r.id.clone(),
            ]);
            if here {
                row.bold().cyan()
            } else if r.role == Role::Focus {
                row.bold()
            } else if share > 0.0 || types > 0 {
                row
            } else {
                row.dim()
            }
        })
        .collect();
    let title = match (&v.typing, &v.query) {
        (Some(t), _) => format!(" search the Cell Ontology: {t}▏ "),
        (None, Some(q)) => format!(" terms matching {q:?} · ← back "),
        (None, None) => {
            let path: Vec<String> = cl
                .lineage(&v.focus)
                .iter()
                .map(|id| cl.name(id).unwrap_or(id).to_string())
                .collect();
            let tail = path.len().saturating_sub(4);
            let prefix = if tail > 0 { "… › " } else { "" };
            format!(" {prefix}{} ", path[tail..].join(" › "))
        }
    };
    let mut state = TableState::default().with_selected(Some(v.sel));
    let t = Table::new(
        rows,
        [
            Constraint::Min(30),
            Constraint::Length(6),
            Constraint::Length(9),
            Constraint::Length(12),
        ],
    )
    .block(pane(title, focused).title_bottom(Line::from(" cell ontology · o: panel tree ").dim()))
    .row_highlight_style(highlight(focused));
    f.render_stateful_widget(t, area, &mut state);
}

fn draw_summary(f: &mut Frame, area: Rect, app: &App) {
    let Some(r) = &app.round else { return };
    let s = r.summary(&app.edits);
    let total = s.iter().map(|(_, n)| n).sum::<usize>().max(1);
    let spans: Vec<Span> = s
        .iter()
        .flat_map(|(l, n)| {
            let pct = 100.0 * *n as f64 / total as f64;
            [
                Span::from(format!(" {l} ")).bold(),
                Span::from(format!("{pct:.1}% ·")).dim(),
            ]
        })
        .collect();
    f.render_widget(
        Paragraph::new(Line::from(spans)).wrap(Wrap { trim: true }),
        area,
    );
}

fn draw_log(f: &mut Frame, area: Rect, app: &App) {
    let h = area.height.saturating_sub(2) as usize;
    let lines: Vec<Line> = app.log[app.log.len().saturating_sub(h)..]
        .iter()
        .map(|l| Line::from(l.as_str()))
        .collect();
    f.render_widget(
        Paragraph::new(lines).block(Block::default().borders(Borders::ALL).title(" log ")),
        area,
    );
}

fn draw_prompt(f: &mut Frame, app: &App) {
    let Some(p) = &app.prompt else { return };
    // Wider for longer questions (a long gene list), then as tall as the
    // question and the answer need once wrapped.
    let screen = f.area();
    let question = p.title.trim();
    // A choice takes keys, not text: no answer line.
    let choice = matches!(p.pending, super::app::Pending::TrajectoryLabels);
    let answer = if choice {
        String::new()
    } else {
        format!("{}▏", p.text)
    };
    let want = question.chars().count().max(answer.chars().count()) as u16 + 4;
    // At least 60 columns (or the screen), at most 90% of the screen.
    let lo = screen.width.min(60);
    let width = want.clamp(lo, (screen.width * 9 / 10).max(lo));
    let inner = width.saturating_sub(2).max(1) as usize;
    // Word wrapping can take a row more than the characters alone.
    let rows = |t: &str| {
        let n = t.chars().count().div_ceil(inner).max(1) as u16;
        n + u16::from(n > 1)
    };
    let height = (rows(question) + 1 + rows(&answer) + 2).min(screen.height);
    let [area] = Layout::vertical([Constraint::Length(height)])
        .flex(Flex::Center)
        .areas(screen);
    let [area] = Layout::horizontal([Constraint::Length(width)])
        .flex(Flex::Center)
        .areas(area);
    let lines = vec![
        Line::from(question.to_string()).bold(),
        Line::default(),
        Line::from(answer),
    ];
    f.render_widget(Clear, area);
    f.render_widget(
        Paragraph::new(lines)
            .style(POPUP)
            .wrap(Wrap { trim: false })
            .block(popup(
                String::new(),
                if choice {
                    " enter/l: pick a file · a: annotate · esc: cancel "
                } else {
                    " enter: ok · esc: cancel "
                },
            )),
        area,
    );
}

/// The order view: the types on screen, then the prior's direct edges with
/// their source and, when a trajectory run recorded them, their verdicts.
fn draw_order(f: &mut Frame, area: Rect, app: &App, v: &super::order::OrderView) {
    let focused = app.focus == Focus::Tree;
    let mark = |t: &str| {
        app.tree_marked
            .iter()
            .position(|m| m == t)
            .map_or("  ".to_string(), |i| format!("{} ", i + 1))
    };
    let mut rows: Vec<Row> = v
        .types
        .iter()
        .map(|(t, n)| {
            Row::new(vec![
                Cell::from(format!("{}{}", mark(t), t)),
                Cell::from(format!("{n}")),
                Cell::from(""),
                Cell::from(""),
            ])
        })
        .collect();
    rows.push(Row::new(vec![Cell::from("")]));
    let edges = v.edges();
    rows.push(
        Row::new(vec![
            Cell::from(format!("{} direct edge(s)", edges.len())),
            Cell::from(""),
            Cell::from("source"),
            Cell::from("data"),
        ])
        .style(Style::default().add_modifier(Modifier::DIM)),
    );
    for e in &edges {
        rows.push(Row::new(vec![
            Cell::from(format!("  {} → {}", e.from, e.to)),
            Cell::from(""),
            Cell::from(e.source.clone()),
            Cell::from(e.data.map_or(String::new(), |d| d.summary())),
        ]));
    }
    for s in v.unrelated() {
        rows.push(Row::new(vec![
            Cell::from(format!("  {} ∥ {}", s.from, s.to)),
            Cell::from(""),
            Cell::from(s.source.as_str()),
            Cell::from("unrelated"),
        ]));
    }
    let candidates = v.candidates();
    if !candidates.is_empty() {
        rows.push(
            Row::new(vec![
                Cell::from(format!(
                    "{} pair(s) the data connects but the prior does not order",
                    candidates.len()
                )),
                Cell::from(""),
                Cell::from(""),
                Cell::from("data"),
            ])
            .style(Style::default().add_modifier(Modifier::DIM)),
        );
        for d in candidates.iter().take(8) {
            rows.push(Row::new(vec![
                Cell::from(format!("  {} — {}", d.a, d.b)),
                Cell::from(""),
                Cell::from(""),
                Cell::from(d.summary()),
            ]));
        }
    }
    let red = Style::default().fg(Color::Red);
    if let Some(e) = &v.error {
        rows.push(Row::new(vec![
            Cell::from(format!("  not a DAG: {e}")).style(red)
        ]));
    }
    for p in &v.problems {
        rows.push(Row::new(vec![
            Cell::from(format!("  unreadable: {p}")).style(red)
        ]));
    }
    let file = v
        .file
        .as_ref()
        .map_or("nowhere to write".to_string(), |p| p.display().to_string());
    let table = Table::new(
        rows,
        [
            Constraint::Min(24),
            Constraint::Length(7),
            Constraint::Length(14),
            Constraint::Min(16),
        ],
    )
    .header(
        Row::new(vec!["type", "cells", "", ""]).style(Style::default().add_modifier(Modifier::DIM)),
    )
    .block(pane(format!(" order · {file} "), focused))
    .row_highlight_style(highlight(focused));
    let mut state = TableState::default().with_selected(if focused { Some(v.sel) } else { None });
    f.render_stateful_widget(table, area, &mut state);
}

/// The figure pane: the current figure, and the exports strip when open.
fn draw_figures(f: &mut Frame, area: Rect, app: &App, v: &super::figure_pane::FigurePane) {
    let focused = app.focus == Focus::Tree;
    let (figure, strip) = if v.exports.open {
        let [a, b] = Layout::vertical([Constraint::Min(6), Constraint::Length(8)]).areas(area);
        (a, Some(b))
    } else {
        (area, None)
    };
    let title = format!(
        " {} ({}/{}) ",
        v.data.title(v.current()),
        v.sel + 1,
        v.panels.len()
    );
    let block = pane(title, focused);
    let inner = block.inner(figure);
    f.render_widget(block, figure);
    match &*v.protocol(inner) {
        Ok(p) => f.render_widget(ratatui_image::Image::new(p), inner),
        Err(e) => f.render_widget(
            Paragraph::new(e.as_str()).wrap(Wrap { trim: true }).dim(),
            inner,
        ),
    }
    if let Some(strip) = strip {
        draw_exports(f, strip, &v.exports, focused);
    }
}

/// The exports strip: the gallery's entries with their check, then the
/// figure files beside the run that it does not list.
fn draw_exports(f: &mut Frame, area: Rect, x: &super::figure_pane::Exports, focused: bool) {
    let now = super::gallery::now();
    let mut rows: Vec<Row> = x
        .gallery
        .entries
        .iter()
        .zip(&x.statuses)
        .map(|(e, s)| {
            let row = Row::new(vec![
                Cell::from(super::app::file_name(&e.path)),
                Cell::from(super::gallery::ago(e.when, now)),
                Cell::from(s.as_str()),
            ]);
            match s.color() {
                Some(c) => row.style(Style::default().fg(c)),
                None => row,
            }
        })
        .collect();
    rows.extend(x.not_listed.iter().map(|p| {
        Row::new(vec![
            Cell::from(super::app::file_name(p)),
            Cell::from(""),
            Cell::from("not listed"),
        ])
        .style(Style::default().add_modifier(Modifier::DIM))
    }));
    let note = x.gallery.set_aside.as_ref().map_or(String::new(), |p| {
        format!(" · unreadable log set aside as {}", p.display())
    });
    let table = Table::new(
        rows,
        [
            Constraint::Min(24),
            Constraint::Length(12),
            Constraint::Length(20),
        ],
    )
    .block(pane(
        format!(
            " exports ({}){note} · enter list · M move · d unlist · D D delete ",
            x.gallery.entries.len()
        ),
        focused,
    ))
    .row_highlight_style(highlight(focused));
    let mut state = TableState::default().with_selected(Some(x.sel));
    f.render_stateful_widget(table, area, &mut state);
}

#[cfg(test)]
mod wrap_tests {
    use super::wrap;

    #[test]
    fn a_long_name_wraps_between_words_and_is_cut_at_the_last_line() {
        assert_eq!(wrap("one two three", 7, 3), ["one two", "three"]);
        assert_eq!(wrap("aa bb cc dd", 5, 2), ["aa bb", "cc dd"]);
        assert_eq!(wrap("aa bb cc dd ee", 5, 2), ["aa bb", "cc d…"]);
        assert_eq!(wrap("abcdefgh", 4, 3), ["abcd"]);
        assert_eq!(wrap("", 4, 3), [""]);
    }
}
