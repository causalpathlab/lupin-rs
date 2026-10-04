pub mod scatter;
pub mod strand;
pub mod topic;

pub use legume_plot::{hull, palette, rasterize, svg_emit};

use rasterize::{DataBounds, Extent};

/// Data → pixel with y pointing up (larger data-y → higher on screen), the
/// convention every layout figure shares. Hull vertices, label anchors and
/// raster layers all go through here so they align.
#[must_use]
pub(crate) fn to_pixel(p: (f32, f32), bounds: &DataBounds, ext: Extent) -> (f32, f32) {
    let (x, y) = bounds.to_pixel(p, ext);
    (x, ext.h as f32 - y)
}

/// The bounding box of the points `idx` of `(x, y)`.
pub(crate) fn bounds_of(x: &[f32], y: &[f32], idx: &[usize]) -> DataBounds {
    let (mut x0, mut x1, mut y0, mut y1) = (f32::MAX, f32::MIN, f32::MAX, f32::MIN);
    for &i in idx {
        x0 = x0.min(x[i]);
        x1 = x1.max(x[i]);
        y0 = y0.min(y[i]);
        y1 = y1.max(y[i]);
    }
    DataBounds::from_minmax(x0, x1, y0, y1)
}

/// Map a label (cell type, batch, …) to a filesystem-safe basename:
/// keep ASCII alphanumerics and `-_.`, replace everything else with `_`.
pub(crate) fn sanitize_filename(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '_' | '-' | '.' => c,
            _ => '_',
        })
        .collect()
}

/// Inverse of [`legume_numeric::matrix::dense_mat_io::axis_id_names`]: the
/// integer id of a `{prefix}{c}` column, or of a legacy bare-integer name.
#[must_use]
pub(crate) fn parse_axis_id(name: &str, prefix: &str) -> Option<i64> {
    name.strip_prefix(prefix)
        .and_then(|rest| rest.parse::<i64>().ok())
        .or_else(|| name.parse::<i64>().ok())
}

/// Every column's axis id, or `None` if any column carries none.
#[must_use]
pub(crate) fn try_parse_axis_ids(cols: &[Box<str>], prefix: &str) -> Option<Vec<i64>> {
    cols.iter().map(|c| parse_axis_id(c, prefix)).collect()
}

/// [`try_parse_axis_ids`], numbering columns `0..n` when they carry no ids.
#[must_use]
pub(crate) fn axis_ids_or_positions(cols: &[Box<str>], prefix: &str) -> Vec<i64> {
    try_parse_axis_ids(cols, prefix).unwrap_or_else(|| (0..cols.len() as i64).collect())
}
