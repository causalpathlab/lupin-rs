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
use ratatui::widgets::{Block, Borders, Cell, Clear, Paragraph, Row, Table, TableState, Wrap};
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
    let [left, right] =
        Layout::horizontal([Constraint::Length(56), Constraint::Min(50)]).areas(body);
    let [clusters, settings] = Layout::vertical([
        Constraint::Min(6),
        Constraint::Length(SETTINGS.len() as u16 + 2),
    ])
    .areas(left);
    let [tree, detail] =
        Layout::vertical([Constraint::Min(8), Constraint::Length(14)]).areas(right);
    let [candidates, genes] =
        Layout::horizontal([Constraint::Percentage(58), Constraint::Percentage(42)]).areas(detail);

    draw_header(f, header, app);
    draw_clusters(f, clusters, app);
    draw_settings(f, settings, app);
    draw_candidates(f, candidates, app);
    draw_genes(f, genes, app);
    draw_tree(f, tree, app);
    draw_summary(f, summary, app);
    draw_log(f, log, app);
    f.render_widget(Paragraph::new(help(app)).dim(), keys);
    if app.prompt.is_some() {
        draw_prompt(f, app);
    }
}

fn help(app: &App) -> &'static str {
    match app.focus {
        Focus::Clusters => {
            " ↑↓ cluster · 1-6 take candidate · enter pick in tree · u unassign · ⌫ undo · ] next flagged · tab pane · s save · e export · r run · q quit"
        }
        Focus::Genes => {
            " ↑↓ gene · space mark · a add as markers of the cluster's label · d drop them · m specific genes / markers · esc back · s save"
        }
        Focus::Tree => match app.tree_mode {
            TreeMode::Panel => {
                " ↑↓ node · enter give the cluster this label · ←→ fold · o Cell Ontology · esc back · tab pane · s save"
            }
            TreeMode::Ontology(_) => {
                " ↑↓ term · enter give the cluster this term · → into · ← up · / search · o panel tree · esc back · s save"
            }
        },
        Focus::Settings => " ↑↓ setting · ←→ change · r run a pass with them · tab pane",
    }
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
    let mut spans = vec![
        Span::from(" lupin annotate ").bold().reversed(),
        Span::from(format!(" {} ", app.source.display())),
    ];
    if !app.edits.is_empty() {
        spans.push(Span::from(format!("· {} unsaved ", app.edits.len())).yellow());
    }
    if app.stale {
        spans.push(Span::from("· settings changed, r to re-run ").yellow());
    }
    spans.push(Span::from(format!("· {}", app.status)));
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
            let flag = if c.flagged() { "?" } else { "" };
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

fn draw_settings(f: &mut Frame, area: Rect, app: &App) {
    let focused = app.focus == Focus::Settings;
    let rows: Vec<Row> = SETTINGS
        .iter()
        .map(|s| Row::new([s.name().to_string(), s.value(&app.args)]))
        .collect();
    let note = if app.fixed_clusters {
        "clusters from the manifest"
    } else {
        "Leiden on the cell embedding"
    };
    let mut state = TableState::default().with_selected(focused.then_some(app.setting));
    let t = Table::new(rows, [Constraint::Length(14), Constraint::Min(8)])
        .block(pane(format!(" pass settings · {note} "), focused))
        .row_highlight_style(highlight(focused));
    f.render_stateful_widget(t, area, &mut state);
}

fn draw_candidates(f: &mut Frame, area: Rect, app: &App) {
    let Some(c) = app.selected() else {
        f.render_widget(
            Block::default().borders(Borders::ALL).title(" cluster "),
            area,
        );
        return;
    };
    let now = app.current_label();
    let num =
        |v: Option<f32>, digits: usize| v.map_or_else(|| "—".into(), |v| format!("{v:.digits$}"));
    let pval = |v: Option<f32>| match v {
        Some(v) if v < 1e-3 => format!("{v:.0e}"),
        v => num(v, 3),
    };
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
    let title = format!(
        " K{} · {} cells · {} ",
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
    let (rows, title): (Vec<Row>, String) = match app.gene_view {
        GeneView::Specific => {
            let rows = c
                .genes
                .iter()
                .map(|(g, s)| {
                    let of = r.marker_of(g);
                    let row = Row::new([
                        mark(g).to_string(),
                        g.clone(),
                        format!("{s:+.2}"),
                        of.join(","),
                    ]);
                    if key.as_deref().is_some_and(|l| of.contains(&l)) {
                        row.green()
                    } else if of.is_empty() {
                        row
                    } else {
                        row.dim()
                    }
                })
                .collect();
            (
                rows,
                " specific genes (log2 FC) · green: marker of its label · m: its markers ".into(),
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
            let rows = listed
                .iter()
                .map(|(g, fc, added)| {
                    let fc = *fc;
                    let row = Row::new([
                        mark(g).to_string(),
                        g.clone(),
                        fc.map_or_else(|| "—".into(), |v| format!("{v:+.2}")),
                        if *added {
                            "added".into()
                        } else {
                            String::new()
                        },
                    ]);
                    match fc {
                        _ if *added => row.yellow(),
                        Some(v) if v > 0.0 => row.green(),
                        Some(_) => row,
                        None => row.dim(),
                    }
                })
                .collect();
            (
                rows,
                format!(
                    " {} markers of {label} (log2 FC in K{}) · m: specific genes ",
                    listed.len(),
                    c.id
                ),
            )
        }
    };
    let mut state = TableState::default().with_selected(focused.then_some(app.gene_sel));
    let t = Table::new(
        rows,
        [
            Constraint::Length(1),
            Constraint::Length(12),
            Constraint::Length(7),
            Constraint::Min(8),
        ],
    )
    .block(pane(title, focused))
    .row_highlight_style(highlight(focused));
    f.render_stateful_widget(t, area, &mut state);
}

fn draw_tree(f: &mut Frame, area: Rect, app: &App) {
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
            let name = format!(
                "{}{glyph} {}{}",
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
            let name = format!(
                "{indent}{glyph} {}{}",
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
    let [area] = Layout::vertical([Constraint::Length(3)])
        .flex(Flex::Center)
        .areas(f.area());
    let [area] = Layout::horizontal([Constraint::Percentage(70)])
        .flex(Flex::Center)
        .areas(area);
    f.render_widget(Clear, area);
    f.render_widget(
        Paragraph::new(format!("{}▏", p.text)).block(
            Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(Color::Yellow))
                .title(p.title.clone())
                .title_bottom(Line::from(" enter: keep · esc: cancel ").dim()),
        ),
        area,
    );
}
