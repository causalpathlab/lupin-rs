//! Browsing the Cell Ontology around a term, to label a cluster with any
//! term, not only the panel's: the term's parents, the term, its children,
//! or the terms a search finds.

use crate::annotate::celltype_tree::ClTerms;
use crate::annotate::markers::label_key;
use crate::annotate::panel_tree::PanelTree;
use ratatui::crossterm::event::KeyCode;
use std::collections::{BTreeMap, BTreeSet};

/// Search hits listed.
const HITS: usize = 200;

/// Where a row stands relative to the focused term.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    Parent,
    Focus,
    Child,
    Hit,
    /// A term of the tree of the terms in the data.
    Tree,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    pub id: String,
    pub role: Role,
    /// How deep a [`Role::Tree`] row sits.
    pub depth: usize,
    /// Terms above a [`Role::Tree`] row folded into its line: each not in
    /// the data and with only this way down, top first.
    pub chain: Vec<String>,
}

impl Row {
    fn new(id: &str, role: Role) -> Self {
        Self {
            id: id.to_string(),
            role,
            depth: 0,
            chain: Vec::new(),
        }
    }
}

/// Which terms the ontology view lists.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Scope {
    /// The terms the data's types sit on, under their lowest common
    /// ancestor, with chains of terms that only lead on folded into one line.
    #[default]
    Data,
    /// The whole ontology, around a term.
    All,
}

/// The ontology pane's state.
/// What [`OntologyView::key`] did with a key.
#[derive(Debug, PartialEq, Eq)]
pub enum ViewKey {
    Taken,
    /// `d`, with no type of the data on a term: nothing to narrow to.
    NoData,
    /// Not one of its keys.
    Other,
}

pub struct OntologyView {
    pub focus: String,
    pub rows: Vec<Row>,
    pub sel: usize,
    /// A search being typed.
    pub typing: Option<String>,
    /// The query whose hits are listed, when they are.
    pub query: Option<String>,
    pub scope: Scope,
    /// The data's types by the term they sit on, with their cell counts.
    pub data: BTreeMap<String, Vec<(String, usize)>>,
    /// Rows whose folded chain is shown term by term.
    expanded: BTreeSet<String>,
    /// Terms of the data's tree listed on all their children, not only
    /// those leading to the data.
    opened: BTreeSet<String>,
    /// Terms of the data's tree showing only the data's terms below them.
    hidden: BTreeSet<String>,
    /// The data's tree's top when climbed above the common ancestor.
    top: Option<String>,
    /// Terms keeping a row of their own, never folded into a chain: the
    /// tops climbed from.
    own: BTreeSet<String>,
}

impl OntologyView {
    /// Focused on `focus`, the cursor on it.
    #[must_use]
    pub fn at(cl: &ClTerms, focus: &str) -> Self {
        let mut v = Self {
            focus: focus.to_string(),
            rows: Vec::new(),
            sel: 0,
            typing: None,
            query: None,
            scope: Scope::All,
            data: BTreeMap::new(),
            expanded: BTreeSet::new(),
            opened: BTreeSet::new(),
            hidden: BTreeSet::new(),
            top: None,
            own: BTreeSet::new(),
        };
        v.refocus(cl, focus);
        v
    }

    /// The terms of `data` (term → its types with their cell counts), the
    /// cursor on `focus` when it is listed; the whole ontology around
    /// `focus` when there is no data.
    #[must_use]
    pub fn in_data(
        cl: &ClTerms,
        data: BTreeMap<String, Vec<(String, usize)>>,
        focus: &str,
    ) -> Self {
        let mut v = Self::at(cl, focus);
        v.data = data;
        if !v.data.is_empty() {
            v.scope = Scope::Data;
            v.show_data(cl, focus);
        }
        v
    }

    /// List the data's tree, the cursor on `id` (or the row folding it in).
    fn show_data(&mut self, cl: &ClTerms, id: &str) {
        self.query = None;
        self.rows = data_rows(cl, self);
        self.sel = self.row_of(id).unwrap_or(0);
    }

    /// The row showing `id`, on its own or folded into its chain.
    fn row_of(&self, id: &str) -> Option<usize> {
        self.rows
            .iter()
            .position(|r| r.id == id || r.chain.iter().any(|c| c == id))
    }

    /// Between the terms in the data and the whole ontology, the cursor on
    /// the same term; `false` when there is no data to show.
    /// The keys every ontology view shares: ↑↓ and pages move, → into a
    /// term (or open a chain), ← up, `/` search, `d` in the data ↔ all, `h`
    /// hide / show the terms not in the data below the selected one.
    pub fn key(&mut self, cl: &ClTerms, code: KeyCode) -> ViewKey {
        if super::app::step(&mut self.sel, self.rows.len(), code) {
            return ViewKey::Taken;
        }
        match code {
            KeyCode::Right => self.enter(cl),
            KeyCode::Left => self.up(cl),
            KeyCode::Char('/') => self.typing = Some(String::new()),
            KeyCode::Char('d') if !self.toggle_scope(cl) => return ViewKey::NoData,
            KeyCode::Char('d') => {}
            KeyCode::Char('h') if self.scope == Scope::Data && self.query.is_none() => {
                if let Some(id) = self.selected().map(|r| r.id.clone()) {
                    if !self.hidden.remove(&id) {
                        self.hidden.insert(id.clone());
                    }
                    self.show_data(cl, &id);
                }
            }
            _ => return ViewKey::Other,
        }
        ViewKey::Taken
    }

    pub fn toggle_scope(&mut self, cl: &ClTerms) -> bool {
        let id = self
            .selected()
            .map_or_else(|| self.focus.clone(), |r| r.id.clone());
        match self.scope {
            Scope::Data => {
                self.scope = Scope::All;
                self.refocus(cl, &id);
            }
            Scope::All if self.data.is_empty() => return false,
            Scope::All => {
                self.scope = Scope::Data;
                self.show_data(cl, &id);
            }
        }
        true
    }

    /// Focus on `id`: its parents, itself and its children, the cursor on it.
    pub fn refocus(&mut self, cl: &ClTerms, id: &str) {
        self.focus = id.to_string();
        self.query = None;
        let mut rows: Vec<Row> = cl
            .parents(id)
            .iter()
            .map(|p| Row::new(p, Role::Parent))
            .collect();
        self.sel = rows.len();
        rows.push(Row::new(id, Role::Focus));
        rows.extend(cl.children(id).iter().map(|c| Row::new(c, Role::Child)));
        self.rows = rows;
    }

    /// List the terms `query` finds: in the data's tree when it lists any of
    /// them, else in the whole ontology (`true`: the view went to it).
    pub fn search(&mut self, cl: &ClTerms, query: &str) -> bool {
        let mut hits = cl.search(query, HITS);
        let mut left_data = false;
        if self.scope == Scope::Data {
            let shown: BTreeSet<&str> = self
                .rows
                .iter()
                .flat_map(|r| std::iter::once(&r.id).chain(&r.chain))
                .map(String::as_str)
                .collect();
            let inside: Vec<String> = hits
                .iter()
                .filter(|h| shown.contains(h.as_str()))
                .cloned()
                .collect();
            if inside.is_empty() && !hits.is_empty() {
                self.scope = Scope::All;
                left_data = true;
            } else {
                hits = inside;
            }
        }
        self.rows = hits.iter().map(|id| Row::new(id, Role::Hit)).collect();
        self.sel = 0;
        self.query = Some(query.to_string());
        left_data
    }

    /// The cursor on `id`, opened on all its children (an opened term is
    /// never folded): in the data's tree when it lists `id` or `id` is above
    /// its top (the tree then starts at it), else the whole ontology around it.
    pub fn open_at(&mut self, cl: &ClTerms, id: &str) {
        if self.scope == Scope::Data && self.row_of(id).is_none() {
            let above = self
                .rows
                .first()
                .is_some_and(|t| cl.ancestors_or_self(&t.id).contains(id));
            if above {
                self.top = Some(id.to_string());
            } else {
                self.scope = Scope::All;
            }
        }
        if self.scope == Scope::Data {
            self.opened.insert(id.to_string());
            self.show_data(cl, id);
        } else {
            self.refocus(cl, id);
        }
    }

    /// Whether `id` shows only the data's terms below it.
    #[must_use]
    pub fn is_hiding(&self, id: &str) -> bool {
        self.hidden.contains(id)
    }

    /// Whether `id` is opened on all its children in the data's tree.
    #[must_use]
    pub fn is_open(&self, id: &str) -> bool {
        self.opened.contains(id)
    }

    #[must_use]
    pub fn selected(&self) -> Option<&Row> {
        self.rows.get(self.sel)
    }

    /// Descend into the selected row's term; in the data's tree, show a
    /// folded chain term by term, else open the term in place on all its
    /// children.
    pub fn enter(&mut self, cl: &ClTerms) {
        let Some(r) = self.selected().cloned() else {
            return;
        };
        if self.scope == Scope::Data {
            // A search hit goes back to the tree, on it; a folded chain
            // unfolds first.
            if self.query.is_none() && (r.chain.is_empty() || !self.expanded.insert(r.id.clone())) {
                self.opened.insert(r.id.clone());
            }
            return self.show_data(cl, &r.id);
        }
        self.refocus(cl, &r.id);
    }

    /// Up to the focus's first parent; out of a search, back to the focus.
    /// In the data's tree: close an opened term, else fold an unfolded
    /// chain again, else up to the row's parent; on the top row, climb to
    /// the term above it, opened on all its children.
    pub fn up(&mut self, cl: &ClTerms) {
        if self.scope == Scope::Data {
            if self.query.is_some() {
                let focus = self.focus.clone();
                return self.show_data(cl, &focus);
            }
            let Some(r) = self.selected().cloned() else {
                return;
            };
            if self.opened.remove(&r.id) || self.expanded.remove(&r.id) {
                return self.show_data(cl, &r.id);
            }
            if r.depth == 0 {
                if let Some(p) = cl.parents(&r.id).first().cloned() {
                    self.opened.insert(p.clone());
                    self.top = Some(p);
                    self.own.insert(r.id.clone());
                    self.show_data(cl, &r.id);
                }
                return;
            }
            if let Some(p) = self.rows[..self.sel]
                .iter()
                .rposition(|p| p.depth < r.depth)
            {
                self.sel = p;
            }
            return;
        }
        if self.query.is_some() {
            let focus = self.focus.clone();
            self.refocus(cl, &focus);
        } else if let Some(p) = cl.parents(&self.focus).first().cloned() {
            let from = self.focus.clone();
            self.refocus(cl, &p);
            // Keep the cursor on the term we came from.
            if let Some(i) = self.rows.iter().position(|r| r.id == from) {
                self.sel = i;
            }
        }
    }
}

/// The tree of `data`'s terms: under their lowest common ancestor, each term
/// below the deepest of its parents in the tree, children by name. A term
/// not in the data with only one way down is folded into the line of the
/// term it leads to, unless that term is `expanded`. An `opened` term lists
/// all its children, in the data or not; neither it nor an `own` term is
/// folded. A `top` above the common ancestor starts the tree higher; below a
/// `hidden` term only the data's terms are listed.
fn data_rows(cl: &ClTerms, v: &OntologyView) -> Vec<Row> {
    let data: BTreeSet<String> = v.data.keys().cloned().collect();
    let (expanded, opened) = (&v.expanded, &v.opened);
    let rows = tree_rows(cl, &data, expanded, opened, &v.own, v.top.as_deref());
    hide_unmatched(rows, &data, &v.hidden)
}

fn tree_rows(
    cl: &ClTerms,
    data: &BTreeSet<String>,
    expanded: &BTreeSet<String>,
    opened: &BTreeSet<String>,
    own: &BTreeSet<String>,
    top: Option<&str>,
) -> Vec<Row> {
    let above: Vec<BTreeSet<String>> = data.iter().map(|d| cl.ancestors_or_self(d)).collect();
    let Some(first) = above.first() else {
        return Vec::new();
    };
    // Every term above the data with its own ancestors, walked once: a
    // term's depth is their number.
    let ancestors: BTreeMap<&str, BTreeSet<String>> = above
        .iter()
        .flatten()
        .map(|t| (t.as_str(), cl.ancestors_or_self(t)))
        .collect();
    let depth_of = |id: &str| ancestors.get(id).map_or(0, BTreeSet::len);
    // The lowest common ancestor: the deepest term above every data term.
    let top = top.map(String::from).or_else(|| {
        first
            .iter()
            .filter(|t| above.iter().all(|a| a.contains(*t)))
            .max_by_key(|t| (depth_of(t), std::cmp::Reverse((*t).clone())))
            .cloned()
    });
    let nodes: BTreeSet<String> = ancestors
        .iter()
        .filter(|(_, a)| top.as_ref().is_none_or(|top| a.contains(top)))
        .map(|(t, _)| (*t).to_string())
        .collect();
    let mut children: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    let mut roots = Vec::new();
    for n in &nodes {
        let parent = cl
            .parents(n)
            .iter()
            .filter(|p| nodes.contains(*p))
            .max_by_key(|p| (depth_of(p), std::cmp::Reverse((*p).clone())));
        match parent {
            Some(p) if top.as_deref() != Some(n.as_str()) => {
                children.entry(p.as_str()).or_default().push(n);
            }
            _ => roots.push(n.as_str()),
        }
    }
    for o in opened {
        let kids = children.entry(o.as_str()).or_default();
        for c in cl.children(o) {
            if !kids.contains(&c.as_str()) {
                kids.push(c);
            }
        }
    }
    let name = |id: &str| cl.name(id).unwrap_or(id).to_lowercase();
    for kids in children.values_mut() {
        kids.sort_by_cached_key(|k| name(k));
    }
    roots.sort_by_key(|k| name(k));
    let mut rows = Vec::new();
    let mut stack: Vec<(&str, usize)> = roots.iter().rev().map(|r| (*r, 0)).collect();
    while let Some((n, d)) = stack.pop() {
        // Fold the terms that only lead on, unless the user unfolded them.
        let (mut chain, mut end) = (Vec::new(), n);
        if d > 0 {
            while !data.contains(end) && !opened.contains(end) && !own.contains(end) {
                match children.get(end).map(Vec::as_slice) {
                    Some([only]) => {
                        chain.push(end.to_string());
                        end = only;
                    }
                    _ => break,
                }
            }
        }
        let mut d = d;
        if expanded.contains(end) {
            for c in chain.drain(..) {
                rows.push(Row {
                    depth: d,
                    ..Row::new(&c, Role::Tree)
                });
                d += 1;
            }
        }
        rows.push(Row {
            depth: d,
            chain,
            ..Row::new(end, Role::Tree)
        });
        for k in children.get(end).into_iter().flatten().rev() {
            stack.push((k, d + 1));
        }
    }
    rows
}

/// Below each row in `hidden`, only the rows of `data`'s terms (their folded
/// chains dropped), each under the nearest of them above it.
fn hide_unmatched(rows: Vec<Row>, data: &BTreeSet<String>, hidden: &BTreeSet<String>) -> Vec<Row> {
    if hidden.is_empty() {
        return rows;
    }
    let mut out = Vec::with_capacity(rows.len());
    let mut rows = rows.into_iter().peekable();
    while let Some(row) = rows.next() {
        if !hidden.contains(&row.id) {
            out.push(row);
            continue;
        }
        let top = row.depth;
        out.push(row);
        // The kept rows above the current one, by their depth in `rows`.
        let mut above: Vec<usize> = Vec::new();
        while let Some(r) = rows.next_if(|r| r.depth > top) {
            if !data.contains(&r.id) {
                continue;
            }
            while above.last().is_some_and(|&d| d >= r.depth) {
                above.pop();
            }
            out.push(Row {
                depth: top + 1 + above.len(),
                chain: Vec::new(),
                ..r
            });
            above.push(r.depth);
        }
    }
    out
}

/// The label a cluster gets from CL term `id`: the panel type sitting on
/// exactly that term when there is one, so a pick agrees with the panel,
/// else the term's name.
#[must_use]
pub fn term_label(cl: &ClTerms, tree: &PanelTree, id: &str) -> String {
    // The first by name when several panel types sit on one term, so a pick
    // names the same label every session.
    tree.typed_terms()
        .filter(|(_, t)| *t == id)
        .map(|(l, _)| l)
        .min()
        .map_or_else(|| label_key(cl.name(id).unwrap_or(id)), String::from)
}

/// The CL term a label names: a panel type's, else the term it matches
/// by name, synonym or abbreviation.
#[must_use]
pub fn term_of(cl: &ClTerms, tree: &PanelTree, label: &str) -> Option<String> {
    if let Some(id) = tree.cl_of(label) {
        return Some(id.to_string());
    }
    let (mapped, _) = cl.map_labels([label]);
    mapped.into_values().next()
}

/// Each panel type on a CL term with every term above it: which types a
/// term covers.
#[must_use]
pub fn type_ancestry(cl: &ClTerms, tree: &PanelTree) -> Vec<(String, BTreeSet<String>)> {
    tree.typed_terms()
        .map(|(l, id)| (l.to_string(), cl.ancestors_or_self(id)))
        .collect()
}

/// Mixed labels: a label that stands for several cell types a cluster's
/// evidence cannot separate, one per line of `mixed_labels.tsv`
/// (`name<TAB>type<TAB>type…`), read from the user's config and the project's
/// `lupin/` folder and added to from the TUI. Only a listed label is a mix:
/// a `+` inside a panel type's name (`CD14+ Monocytes`) means nothing.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mixed {
    by_name: std::collections::BTreeMap<String, Vec<String>>,
}

/// The file mixed labels are kept in.
pub const MIXED: &str = "mixed_labels.tsv";

/// How mixed-label names compare: as lupin labels, ignoring case.
fn name_key(name: &str) -> String {
    label_key(&name.trim().to_lowercase())
}

impl Mixed {
    fn extend(&mut self, text: &str) {
        for line in text.lines().map(str::trim) {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split('\t').map(str::trim);
            let name = cols.next().unwrap_or_default();
            let parts: Vec<String> = cols.filter(|p| !p.is_empty()).map(label_key).collect();
            if !name.is_empty() && parts.len() > 1 {
                self.by_name.insert(name_key(name), parts);
            }
        }
    }

    /// The user's and the project's mixed labels.
    pub fn load(search: &crate::manifest::data_files::SearchPath) -> anyhow::Result<Self> {
        let mut m = Self::default();
        for text in search.user_and_project(MIXED)? {
            m.extend(&text);
        }
        Ok(m)
    }

    /// Name `parts` `name` from now on, and append it to `file` (once).
    pub fn add(
        &mut self,
        name: &str,
        parts: &[String],
        file: &std::path::Path,
    ) -> anyhow::Result<()> {
        if self
            .by_name
            .get(&name_key(name))
            .is_some_and(|p| p == parts)
        {
            return Ok(());
        }
        crate::manifest::data_files::append_line(
            file,
            "Mixed cell-type labels: the name, then the types it stands for, tab-separated.",
            &std::iter::once(name)
                .chain(parts.iter().map(String::as_str))
                .collect::<Vec<_>>()
                .join("\t"),
        )?;
        self.by_name.insert(name_key(name), parts.to_vec());
        Ok(())
    }

    /// The cell types `label` stands for, when it is a listed mix.
    #[must_use]
    pub fn parts(&self, label: &str) -> Option<Vec<String>> {
        self.by_name.get(&name_key(label)).cloned()
    }
}

#[cfg(test)]
#[path = "tests/ontology.rs"]
mod tests;
