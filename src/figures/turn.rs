//! A figure set on its page turned a quarter: drawn along the page's other
//! side, the PDF page turned (`/Rotate`), so its text, set over the raster
//! as `legume_tui::export::write_pdf` sets it, turns with it and stays text.
//! On paper it is turned left, as a journal sets a landscape figure on a
//! portrait page (the reader turns the page to read it); on a slide, which
//! no one turns, right, so no label reads upside down. `legume_tui`'s writer
//! has no turned page; this one writes one page as it does, turned.

use anyhow::{Context, Result};
use legume_tui::export::font::{helvetica_width, win_ansi};
use legume_tui::export::Page;
use legume_tui::map::colour::{to_srgb_senna, Lin};
use pdf_writer::types::{LineJoinStyle, TextRenderingMode};
use pdf_writer::{Content, Filter, Finish, Name, Pdf, Rect, Ref, Str};
use std::path::Path;

/// Point size of text whose raster line is `line` pixels tall, and where its
/// baseline sits below the top of that line, as fractions of it (as
/// `legume_tui::export` sets them).
const EM_PER_LINE: f32 = 0.74;
const BASELINE: f32 = 0.78;

/// Which way a figure is turned onto its page, a quarter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Turn {
    Left,
    Right,
}

impl Turn {
    /// Degrees clockwise a PDF reader shows the page turned.
    pub(crate) fn degrees(self) -> i32 {
        match self {
            Turn::Left => 270,
            Turn::Right => 90,
        }
    }

    /// `img` turned so.
    pub(crate) fn image(self, img: &image::RgbaImage) -> image::RgbaImage {
        match self {
            Turn::Left => image::imageops::rotate270(img),
            Turn::Right => image::imageops::rotate90(img),
        }
    }
}

/// `page` as a one-page PDF at `path`, `dpi` raster pixels to the inch,
/// shown turned `turn`.
pub(crate) fn write_turned_pdf(page: &Page, dpi: f32, turn: Turn, path: &Path) -> Result<()> {
    anyhow::ensure!(
        dpi.is_finite() && dpi > 0.0,
        "{dpi} dpi is not a resolution"
    );
    let (w, h) = (page.img.width(), page.img.height());
    anyhow::ensure!(w > 0 && h > 0, "the page is empty");
    let (catalog, tree, regular, bold) = (Ref::new(1), Ref::new(2), Ref::new(3), Ref::new(4));
    let (page_id, content_id, image_id) = (Ref::new(5), Ref::new(6), Ref::new(7));
    let mut pdf = Pdf::new();
    pdf.catalog(catalog).pages(tree);
    pdf.pages(tree).kids([page_id]).count(1);
    for (id, name) in [(regular, &b"Helvetica"[..]), (bold, b"Helvetica-Bold")] {
        pdf.type1_font(id)
            .base_font(Name(name))
            .encoding_predefined(Name(b"WinAnsiEncoding"));
    }
    let pt = 72.0 / dpi;
    let (wpt, hpt) = (w as f32 * pt, h as f32 * pt);
    let mut p = pdf.page(page_id);
    p.media_box(Rect::new(0.0, 0.0, wpt, hpt))
        .rotate(turn.degrees())
        .parent(tree)
        .contents(content_id);
    let mut res = p.resources();
    res.x_objects().pair(Name(b"Im"), image_id);
    res.fonts()
        .pair(Name(b"F1"), regular)
        .pair(Name(b"F2"), bold);
    res.finish();
    p.finish();

    let mut rgb = Vec::with_capacity(page.img.as_raw().len() / 4 * 3);
    rgb.extend(page.img.pixels().flat_map(|p| [p[0], p[1], p[2]]));
    let packed = miniz_oxide::deflate::compress_to_vec_zlib(&rgb, 3);
    let mut im = pdf.image_xobject(image_id, &packed);
    im.filter(Filter::FlateDecode);
    im.width(w as i32).height(h as i32).bits_per_component(8);
    im.color_space().device_rgb();
    im.finish();

    let mut c = Content::new();
    c.save_state()
        .transform([wpt, 0.0, 0.0, hpt, 0.0, 0.0])
        .x_object(Name(b"Im"))
        .restore_state();
    c.begin_text().set_line_join(LineJoinStyle::RoundJoin);
    for run in &page.text {
        let size = run.line * EM_PER_LINE * pt;
        let bytes = win_ansi(&run.text);
        let width = helvetica_width(&bytes, run.bold) * size;
        let mut x = run.x * pt;
        if run.centred {
            x += 0.5 * (run.w * pt - width);
        }
        // Turned a quarter left, the baseline is a vertical line
        // `BASELINE` of a line in from the column's left.
        let matrix = if run.vertical {
            [
                0.0,
                1.0,
                -1.0,
                0.0,
                (run.x + BASELINE * run.line) * pt,
                hpt - run.y * pt,
            ]
        } else {
            let y = hpt - (run.y + BASELINE * run.line) * pt;
            [1.0, 0.0, 0.0, 1.0, x, y]
        };
        let font = Name(if run.bold { b"F2" } else { b"F1" });
        c.set_font(font, size);
        // The halo: a stroke in the page colour under the letters.
        let [r, g, b] = srgb(run.bg);
        c.set_stroke_rgb(r, g, b)
            .set_line_width(0.22 * size)
            .set_text_rendering_mode(TextRenderingMode::Stroke)
            .set_text_matrix(matrix)
            .show(Str(&bytes));
        let [r, g, b] = srgb(run.ink);
        c.set_fill_rgb(r, g, b)
            .set_text_rendering_mode(TextRenderingMode::Fill)
            .set_text_matrix(matrix)
            .show(Str(&bytes));
    }
    c.end_text();
    pdf.stream(content_id, &c.finish());
    std::fs::write(path, pdf.finish()).with_context(|| format!("writing {}", path.display()))
}

/// A linear colour as the PDF's sRGB fractions, by senna's table.
fn srgb(c: Lin) -> [f32; 3] {
    to_srgb_senna(c).map(|v| f32::from(v) / 255.0)
}
