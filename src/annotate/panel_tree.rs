//! A marker panel's cell types on a tree, to pick a cluster's label from.
//!
//! On the Cell Ontology, each type hangs under the analysis classes (see
//! [`ClTerms`]) it shares with another panel type, nearest first, so the
//! tree offers every level worth calling a cluster at: a fine type, or a
//! class over several. Without an ontology, the groups of types that share
//! markers ([`TypeTree::from_markers`]) stand over their types. Folding a
//! node only hides its subtree; it changes no label.
//!
//! Pure, like [`crate::annotate::celltype_tree`].

use crate::annotate::celltype_tree::{panel_types, ClTerms, TreeSource, TypeTree};
use crate::annotate::markers::label_key;
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// One node: a CL term, a marker-sharing group, or a panel type the
/// ontology has no term for.
#[derive(Debug, Clone)]
pub struct PanelNode {
    /// The CL term, for an ontology node.
    pub cl_id: Option<String>,
    pub name: String,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// Panel types (in [`label_key`] form) placed on this node itself.
    pub labels: Vec<String>,
    /// 0 at a root.
    pub depth: usize,
}

impl PanelNode {
    /// A node on its own: no parent, children or panel types yet.
    fn new(cl_id: Option<String>, name: String) -> Self {
        Self {
            cl_id,
            name,
            parent: None,
            children: Vec::new(),
            labels: Vec::new(),
            depth: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PanelTree {
    pub nodes: Vec<PanelNode>,
    pub roots: Vec<usize>,
    pub source: TreeSource,
    pub release: Option<String>,
    /// Each panel type's ([`label_key`]) node.
    of_label: HashMap<String, usize>,
    folded: BTreeSet<usize>,
}

impl PanelTree {
    /// The panel's types on the Cell Ontology. `None` when fewer than two
    /// types match a term, as for [`TypeTree::from_ontology`].
    #[must_use]
    pub fn from_ontology(terms: &ClTerms, panel: &[(String, String)]) -> Option<Self> {
        let types = panel_types(panel);
        let (mapped, unmapped) = terms.map_labels(types.iter().map(String::as_str));
        if mapped.len() < 2 {
            return None;
        }
        let ancestry: BTreeMap<&str, BTreeSet<String>> = mapped
            .values()
            .map(|id| (id.as_str(), terms.ancestors_or_self(id)))
            .collect();
        // The nodes: every matched term, and every analysis class above two
        // or more of them.
        let mut sharers: BTreeMap<&str, usize> = BTreeMap::new();
        for a in ancestry.values() {
            for t in a.iter().filter(|t| terms.is_class(t)) {
                *sharers.entry(t.as_str()).or_default() += 1;
            }
        }
        let ids: BTreeSet<&str> = ancestry
            .keys()
            .copied()
            .chain(sharers.iter().filter(|(_, &n)| n >= 2).map(|(t, _)| *t))
            .collect();
        // How far below the ontology's top each node is, to order ancestors.
        let rank: BTreeMap<&str, usize> = ids.iter().map(|id| (*id, terms.depth(id))).collect();

        let mut tree = Self::empty(TreeSource::CellOntology, terms.release.clone());
        let index: BTreeMap<&str, usize> = ids.iter().enumerate().map(|(i, id)| (*id, i)).collect();
        for id in &ids {
            tree.nodes.push(PanelNode::new(
                Some((*id).to_string()),
                label_key(terms.name(id).unwrap_or(id)),
            ));
        }
        // Each node hangs under its nearest ancestor among the nodes.
        for (i, id) in ids.iter().enumerate() {
            let above = ancestry
                .get(id)
                .cloned()
                .unwrap_or_else(|| terms.ancestors_or_self(id));
            tree.nodes[i].parent = above
                .iter()
                .filter(|t| t.as_str() != *id)
                .filter_map(|t| index.get(t.as_str()).map(|&j| (t, j)))
                .max_by(|(x, _), (y, _)| {
                    rank[x.as_str()]
                        .cmp(&rank[y.as_str()])
                        .then_with(|| y.cmp(x))
                })
                .map(|(_, j)| j);
        }
        for (label, id) in &mapped {
            let i = index[id.as_str()];
            let key = label_key(label);
            tree.nodes[i].labels.push(key.clone());
            tree.of_label.insert(key, i);
        }
        // A type the ontology does not have goes where the first round's grouping puts it
        // ([`TypeTree::from_ontology`]): under the class whose types share most of its markers,
        // with the types it shares them with, or on its own when it shares none. Standing alone
        // is a family of its own at the top, which costs every test there.
        let unmapped: BTreeSet<String> = unmapped.iter().map(|l| label_key(l)).collect();
        let groups = TypeTree::from_ontology(terms, panel).map_or_else(Vec::new, |t| t.groups);
        for g in groups {
            let loose: Vec<&String> = g
                .members
                .iter()
                .filter(|m| unmapped.contains(&label_key(m)))
                .collect();
            if loose.is_empty() {
                continue;
            }
            let parent = match g.cl_id.as_deref().and_then(|id| index.get(id)) {
                Some(&i) => Some(i),
                None if g.members.len() > 1 => {
                    tree.nodes.push(PanelNode::new(None, label_key(&g.name)));
                    Some(tree.nodes.len() - 1)
                }
                None => None,
            };
            for m in loose {
                tree.push_leaf(m, parent);
            }
        }
        tree.link();
        Some(tree)
    }

    /// The groups of a marker-sharing [`TypeTree`], each over its types.
    #[must_use]
    pub fn from_type_tree(types: &TypeTree) -> Self {
        let mut tree = Self::empty(types.source.clone(), types.release.clone());
        for g in &types.groups {
            if let [only] = &g.members[..] {
                tree.push_leaf(only, None);
                continue;
            }
            tree.nodes
                .push(PanelNode::new(g.cl_id.clone(), g.name.clone()));
            let group = tree.nodes.len() - 1;
            for m in &g.members {
                tree.push_leaf(m, Some(group));
            }
        }
        tree.link();
        tree
    }

    fn empty(source: TreeSource, release: Option<String>) -> Self {
        Self {
            nodes: Vec::new(),
            roots: Vec::new(),
            source,
            release,
            of_label: HashMap::new(),
            folded: BTreeSet::new(),
        }
    }

    fn push_leaf(&mut self, label: &str, parent: Option<usize>) {
        let key = label_key(label);
        self.nodes.push(PanelNode {
            parent,
            labels: vec![key.clone()],
            ..PanelNode::new(None, key.clone())
        });
        self.of_label.insert(key, self.nodes.len() - 1);
    }

    /// Children, roots and depths from the parents; siblings by name, the
    /// tree before types standing alone.
    fn link(&mut self) {
        for i in 0..self.nodes.len() {
            match self.nodes[i].parent {
                Some(p) => self.nodes[p].children.push(i),
                None => self.roots.push(i),
            }
        }
        let by_name = |nodes: &[PanelNode], v: &mut Vec<usize>| {
            v.sort_by(|&a, &b| nodes[a].name.cmp(&nodes[b].name));
        };
        let mut roots = std::mem::take(&mut self.roots);
        by_name(&self.nodes, &mut roots);
        roots.sort_by_key(|&r| self.nodes[r].children.is_empty());
        let mut stack: Vec<(usize, usize)> = roots.iter().map(|&r| (r, 0)).collect();
        while let Some((i, d)) = stack.pop() {
            self.nodes[i].depth = d;
            let mut kids = std::mem::take(&mut self.nodes[i].children);
            by_name(&self.nodes, &mut kids);
            stack.extend(kids.iter().map(|&c| (c, d + 1)));
            self.nodes[i].children = kids;
        }
        self.roots = roots;
    }

    /// The label a cluster gets from node `i`: a node that is a panel type
    /// keeps the panel's label (`T_cells` for CL `T cell`), else its name.
    #[must_use]
    pub fn label(&self, i: usize) -> &str {
        let n = &self.nodes[i];
        n.labels.first().map_or(n.name.as_str(), String::as_str)
    }

    /// The CL term panel type `label` sits on, if the ontology has it.
    #[must_use]
    pub fn cl_of(&self, label: &str) -> Option<&str> {
        let i = *self.of_label.get(&label_key(label))?;
        self.nodes[i].cl_id.as_deref()
    }

    /// The panel types (label keys) placed on a CL term, with their terms.
    pub fn typed_terms(&self) -> impl Iterator<Item = (&str, &str)> {
        self.of_label.iter().filter_map(|(l, &i)| {
            self.nodes[i]
                .cl_id
                .as_deref()
                .filter(|_| self.nodes[i].labels.contains(l))
                .map(|id| (l.as_str(), id))
        })
    }

    /// The node a label names: a panel type, or a node called by it.
    #[must_use]
    pub fn node_of(&self, label: &str) -> Option<usize> {
        let key = label_key(label);
        self.of_label
            .get(&key)
            .copied()
            .or_else(|| (0..self.nodes.len()).find(|&i| self.label(i) == key))
    }

    /// `i` and its ancestors, root first.
    #[must_use]
    pub fn path(&self, i: usize) -> Vec<usize> {
        let mut out = vec![i];
        while let Some(p) = self.nodes[*out.last().unwrap_or(&i)].parent {
            out.push(p);
        }
        out.reverse();
        out
    }

    #[must_use]
    pub fn is_folded(&self, i: usize) -> bool {
        self.folded.contains(&i)
    }

    /// `←` on node `i`, as in a file tree: fold an open branch, else the
    /// parent to move to.
    pub fn left(&mut self, i: usize) -> Option<usize> {
        if self.nodes[i].children.is_empty() || self.is_folded(i) {
            self.nodes[i].parent
        } else {
            self.fold(i, true);
            None
        }
    }

    /// Hide (`fold`) or show a branching node's subtree.
    pub fn fold(&mut self, i: usize, fold: bool) {
        if fold && !self.nodes[i].children.is_empty() {
            self.folded.insert(i);
        } else {
            self.folded.remove(&i);
        }
    }

    /// Unfold whatever hides `i`.
    pub fn reveal(&mut self, i: usize) {
        for a in self.path(i) {
            if a != i {
                self.folded.remove(&a);
            }
        }
    }

    /// Visible rows, depth-first; a folded node's subtree is hidden.
    #[must_use]
    pub fn visible(&self) -> Vec<usize> {
        let mut out = Vec::new();
        let mut stack: Vec<usize> = self.roots.iter().rev().copied().collect();
        while let Some(i) = stack.pop() {
            out.push(i);
            if !self.folded.contains(&i) {
                stack.extend(self.nodes[i].children.iter().rev().copied());
            }
        }
        out
    }

    /// This tree as TreeBH's hypothesis tree over `celltypes` (the columns of an enrichment
    /// pass): a root over the tree's roots, then every node; each cell type on a leaf of its own,
    /// which is its node when that node is childless and holds only it, else a new leaf under the
    /// node, so a coarse type (`T_cells` over `T_helper_cells`) is tested as its own hypothesis
    /// inside its class's family. Types not on the tree are left to hang under the root.
    #[must_use]
    pub fn treebh(&self, celltypes: &[Box<str>]) -> enrichment::treebh::TypeTree {
        let n = self.nodes.len();
        let mut children: Vec<Vec<usize>> = vec![self.roots.iter().map(|&r| r + 1).collect()];
        children.extend(
            self.nodes
                .iter()
                .map(|nd| nd.children.iter().map(|&c| c + 1).collect()),
        );
        debug_assert_eq!(children.len(), n + 1);
        let leaf = celltypes
            .iter()
            .map(|t| {
                let i = *self.of_label.get(&label_key(t))?;
                let node = &self.nodes[i];
                if node.children.is_empty() && node.labels.len() == 1 {
                    return Some(i + 1);
                }
                children.push(Vec::new());
                let l = children.len() - 1;
                children[i + 1].push(l);
                Some(l)
            })
            .collect();
        enrichment::treebh::TypeTree {
            children,
            root: 0,
            leaf,
        }
    }

    /// The panel types (label keys) at or under `i`.
    #[must_use]
    pub fn labels_under(&self, i: usize) -> Vec<&str> {
        let mut out = Vec::new();
        let mut stack = vec![i];
        while let Some(n) = stack.pop() {
            out.extend(self.nodes[n].labels.iter().map(String::as_str));
            stack.extend(self.nodes[n].children.iter().copied());
        }
        out
    }
}

#[cfg(test)]
#[path = "tests/panel_tree.rs"]
mod tests;
