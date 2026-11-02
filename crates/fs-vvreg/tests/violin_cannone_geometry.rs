//! M8-W3c violin geometry specimen (IRP bead `frankensim-irp-v5-instrument-realism-pexfp.199`): Paganini's
//! violin "il Cannone" (Guarneri del Gesu, 1743) from the ESRF synchrotron CT data collection
//! 10.15151/ESRF-DC-2413577471 (CC BY 4.0).
//!
//! A data-form and consistency check, never a physical claim. It checks that the five specimen tables parse
//! and agree with each other, that every derived scalar equals the data table it was computed from, that the
//! data tables are well formed, and that the map points lie inside the plate outlines. The specimen is NOT the
//! UIowa violin of the recordings, and the boundary stays recorded in `unknowns.tsv`.

use fs_vvreg::specimen::{
    Dimensions, KeyedTable, SPECIMEN_GEOMETRY_ROOT, Sigma, Specimen, SpecimenDescription,
    SpecimenTable,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

const SPEC: &str = "violin-cannone-guarneri-1743-esrf2026";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

fn read(file: &str) -> String {
    let path = root().join(SPECIMEN_GEOMETRY_ROOT).join(SPEC).join(file);
    fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn load() -> SpecimenDescription {
    let text = |table: SpecimenTable| read(table.file_name());
    SpecimenDescription::try_new(
        Specimen::parse(&text(SpecimenTable::Specimen)).expect("specimen.tsv"),
        KeyedTable::parse(SpecimenTable::Frames, &text(SpecimenTable::Frames)).expect("frames"),
        Dimensions::parse(&text(SpecimenTable::Dimensions)).expect("dimensions"),
        KeyedTable::parse(SpecimenTable::Parts, &text(SpecimenTable::Parts)).expect("parts"),
        KeyedTable::parse(SpecimenTable::Unknowns, &text(SpecimenTable::Unknowns))
            .expect("unknowns"),
        None,
    )
    .expect("tables agree with each other")
}

/// Rows of a data table: schema line and header checked, fields split, every row as wide as the header.
fn table(file: &str, header: &str) -> Vec<Vec<String>> {
    let text = read(file);
    let mut lines = text.lines();
    let schema = lines.next().expect("schema line");
    assert!(
        schema.starts_with("frankensim.specimen-") && schema.ends_with(".v1"),
        "{file} schema"
    );
    assert_eq!(lines.next(), Some(header), "{file} header");
    let width = header.split('\t').count();
    lines
        .map(|l| {
            let row: Vec<String> = l.split('\t').map(str::to_owned).collect();
            assert_eq!(row.len(), width, "{file}: {l}");
            row
        })
        .collect()
}

fn num(s: &str) -> f64 {
    s.parse::<f64>().unwrap_or_else(|_| panic!("number {s}"))
}

fn value(d: &Dimensions, part: &str, feature: &str, quantity: &str) -> f64 {
    d.get(part, feature, quantity)
        .unwrap_or_else(|| panic!("missing {part}/{feature}/{quantity}"))
        .value_si
}

#[test]
fn cannone_tables_parse_and_agree() {
    let d = load();
    assert_eq!(d.specimen.id(), SPEC);
    assert_eq!(d.specimen.source_ids(), ["esrf-il-cannone-ct-2026"]);
    assert!(d.dimensions.rows().count() >= 100);
    for part in d.dimensions.parts() {
        assert!(d.frames.row(part).is_some(), "frame for {part}");
        assert!(d.parts.row(part).is_some(), "parts row for {part}");
    }
    let relation = d.specimen.get("target_relation").expect("target_relation");
    assert!(relation.contains("NOT the UIowa violin"));
    assert_eq!(d.specimen.get("doi"), Some("10.15151/ESRF-DC-2413577471"));
    for id in [
        "u-not-iowa-violin",
        "u-tension-state",
        "u-wood-material",
        "u-segmentation-error",
        "u-repairs-and-modifications",
    ] {
        assert!(d.unknowns.row(id).is_some(), "unknown {id} stays recorded");
    }
}

#[test]
fn derived_scalars_equal_the_tables_they_came_from() {
    let dims = load().dimensions;
    // Soundpost: equivalent diameter from the mean section area.
    let post = table(
        "soundpost-sections.tsv",
        "y_mm\tarea_mm2\tcentroid_x_mm\tcentroid_z_mm\tx_min_mm\tx_max_mm\tz_min_mm\tz_max_mm",
    );
    let mean_area = post.iter().map(|r| num(&r[1])).sum::<f64>() / post.len() as f64;
    let diameter = 2.0 * (mean_area / std::f64::consts::PI).sqrt() * 1e-3;
    let stated = value(&dims, "soundpost", "section", "equivalent-diameter-mean");
    assert!(
        (diameter - stated).abs() < 1e-8,
        "{diameter} against {stated}"
    );
    // Bass bar: largest section area.
    let bar = table(
        "bass-bar-sections.tsv",
        "z_mm\tarea_mm2\tcentroid_x_mm\tcentroid_y_mm\tx_min_mm\tx_max_mm\ty_min_mm\ty_max_mm",
    );
    let largest = bar.iter().map(|r| num(&r[1])).fold(0.0, f64::max) * 1e-6;
    assert!((largest - value(&dims, "bass-bar", "section", "area-max")).abs() < 1e-10);
    // Top plate: widest outline row.
    let outline = table("top-plate-outline.tsv", "z_mm\tx_left_mm\tx_right_mm");
    let widest = outline
        .iter()
        .map(|r| num(&r[2]) - num(&r[1]))
        .fold(0.0, f64::max)
        * 1e-3;
    assert!((widest - value(&dims, "top-plate", "footprint", "width-max")).abs() < 1e-8);
    // First f-hole: raster area from its runs (0.5 mm rows).
    let holes = table("top-plate-fholes.tsv", "hole\tz_mm\tx_start_mm\tx_end_mm");
    let area: f64 = holes
        .iter()
        .filter(|r| r[0] == "1")
        .map(|r| (num(&r[3]) - num(&r[2])) * 0.5)
        .sum::<f64>()
        * 1e-6;
    let stated = value(&dims, "top-plate", "f-hole-1", "area");
    assert!(
        ((area - stated) / stated).abs() < 1e-5,
        "{area} against {stated}"
    );
    // Vibrating length: mean nut contact of D, A, E minus mean bridge contact of the four strings.
    let contacts = table(
        "string-contacts.tsv",
        "string\tcontact\tz_mm\tx_mm\ty_mm\tslope_dy_dz_before\tslope_dy_dz_after",
    );
    let mean = |contact: &str, strings: &[&str]| {
        let z: Vec<f64> = contacts
            .iter()
            .filter(|r| r[1] == contact && strings.contains(&r[0].as_str()))
            .map(|r| num(&r[2]))
            .collect();
        assert_eq!(z.len(), strings.len());
        z.iter().sum::<f64>() / z.len() as f64
    };
    let length = (mean("nut", &["d", "a", "e"]) - mean("bridge", &["g", "d", "a", "e"])) * 1e-3;
    let stated = value(
        &dims,
        "fingerboard",
        "nut-to-bridge",
        "string-vibrating-length",
    );
    assert!((length - stated).abs() < 1e-6, "{length} against {stated}");
}

#[test]
fn bar_post_and_strings_sit_on_the_sides_the_frame_text_says() {
    let dims = load().dimensions;
    // x is positive toward the treble side: bass bar at negative x, soundpost at positive x.
    assert!(value(&dims, "bass-bar", "position", "centroid-x-at-bridge-z") < 0.0);
    assert!(value(&dims, "soundpost", "position", "centroid-x-at-mid-length") > 0.0);
    // The strings are ordered G, D, A, E from negative to positive x at the bridge and thin out the same way.
    let x = |s: &str| value(&dims, s, "bridge-contact", "x");
    assert!(
        x("string-g") < x("string-d")
            && x("string-d") < x("string-a")
            && x("string-a") < x("string-e")
    );
    let diameter = |s: &str| value(&dims, s, "section", "equivalent-diameter-median");
    assert!(
        diameter("string-g") >= diameter("string-d") && diameter("string-d") > diameter("string-a")
    );
    assert!(diameter("string-a") > diameter("string-e"));
    // The back plate is thicker than the top plate at the median, and both are thinner than the ribs are tall.
    let top = value(&dims, "top-plate", "thickness", "normal-median-interior");
    let back = value(&dims, "back-plate", "thickness", "normal-median-interior");
    assert!(back > top && back < value(&dims, "ribs", "section", "height-y-median"));
}

#[test]
fn plate_maps_are_well_formed_and_inside_their_outlines() {
    let header = "x_mm\tz_mm\touter_y_mm\tinner_y_mm\tthickness_normal_mm\tdist_to_outline_mm\tvalid_cells_of_36";
    for (plate, top) in [("top-plate", true), ("back-plate", false)] {
        let maps = table(&format!("{plate}-maps.tsv"), header);
        let outline = table(
            &format!("{plate}-outline.tsv"),
            "z_mm\tx_left_mm\tx_right_mm",
        );
        assert!(maps.len() > 3000, "{plate}: {} cells", maps.len());
        for row in &maps {
            let [x, z, outer, inner, t, dist, valid] = [0, 1, 2, 3, 4, 5, 6].map(|i| num(&row[i]));
            assert!(t > 0.0 && t <= 8.0, "{plate} thickness {t}");
            assert!(dist >= 0.0 && (18.0..=36.0).contains(&valid));
            // The outer surface is on the string side for the top plate and on the far side for the back.
            assert!(
                if top { outer < inner } else { outer > inner },
                "{plate} {x} {z}"
            );
            let near = outline
                .iter()
                .find(|r| (num(&r[0]) - z).abs() <= 0.75)
                .unwrap_or_else(|| panic!("{plate}: no outline row at z {z}"));
            assert!(
                x >= num(&near[1]) - 3.0 && x <= num(&near[2]) + 3.0,
                "{plate}: x {x} outside at z {z}"
            );
        }
    }
}

#[test]
fn section_tables_are_well_formed() {
    let sections =
        "z_mm\tarea_mm2\tcentroid_x_mm\tcentroid_y_mm\tx_min_mm\tx_max_mm\ty_min_mm\ty_max_mm";
    for file in [
        "bass-bar-sections.tsv",
        "neck-sections.tsv",
        "fingerboard-sections.tsv",
        "tail-block-sections.tsv",
        "string-g-sections.tsv",
        "string-d-sections.tsv",
        "string-a-sections.tsv",
        "string-e-sections.tsv",
    ] {
        let rows = table(file, sections);
        assert!(!rows.is_empty(), "{file}");
        for row in &rows {
            assert!(num(&row[1]) > 0.0, "{file} area");
            assert!(
                num(&row[4]) <= num(&row[5]) && num(&row[6]) <= num(&row[7]),
                "{file} extents"
            );
        }
    }
    let rings = table(
        "ribs-sections.tsv",
        "theta_deg\tcentre_x_mm\tcentre_z_mm\tarea_mm2\tr_min_mm\tr_max_mm\ty_min_mm\ty_max_mm\tmean_width_mm",
    );
    assert!(rings.len() >= 170);
}

/// Open cells of each hole of a raster-runs table, on a grid of `cell_mm` with cell `(i, j)` covering
/// `[i, i + 1) x [j, j + 1)` in cells; a value off the grid fails.
fn raster_cells(file: &str, cell_mm: f64) -> BTreeMap<u32, BTreeSet<(i64, i64)>> {
    let grid = |v: f64| {
        let k = (v / cell_mm).round();
        assert!((k * cell_mm - v).abs() < 1e-6, "{file}: {v} off the {cell_mm} mm grid");
        k as i64
    };
    let mut holes: BTreeMap<u32, BTreeSet<(i64, i64)>> = BTreeMap::new();
    for row in table(file, "hole\tz_mm\tx_start_mm\tx_end_mm") {
        let j = grid(num(&row[1]) - 0.5 * cell_mm);
        let cells = holes.entry(row[0].parse().unwrap()).or_default();
        for i in grid(num(&row[2]))..grid(num(&row[3])) {
            assert!(cells.insert((i, j)), "{file}: cell listed twice");
        }
    }
    holes
}

/// The 0.5 mm cells with at least `min_open` of their 25 cells of 0.1 mm open.
fn coarsened(fine: &BTreeSet<(i64, i64)>, min_open: usize) -> BTreeSet<(i64, i64)> {
    let mut count: BTreeMap<(i64, i64), usize> = BTreeMap::new();
    for &(i, j) in fine {
        *count.entry((i.div_euclid(5), j.div_euclid(5))).or_default() += 1;
    }
    count.into_iter().filter(|&(_, n)| n >= min_open).map(|(c, _)| c).collect()
}

/// The wall-table value of `(hole, quantity)`.
fn wall(hole: &str, quantity: &str) -> f64 {
    let rows = table("top-plate-fhole-walls.tsv", "hole\tquantity\tvalue\tunit\tnote");
    let row = rows
        .iter()
        .find(|r| r[0] == hole && r[1] == quantity)
        .unwrap_or_else(|| panic!("wall row {hole} {quantity}"));
    num(&row[2])
}

#[test]
fn f_hole_raster_is_the_half_open_coarsening_of_the_0p1_mm_openings() {
    // top-plate-fholes.tsv is fhole_openings.py's 0.5 mm raster: a cell is open when at least 13 of its 25 cells of
    // 0.1 mm are open through the plate (top-plate-fholes-0p1mm.tsv). Planted negative: the rule it replaced, a cell
    // open only when all 25 are ("edge as wood"), gives other cells and falls short of the stated area by more than its
    // sigma.
    let dims = load().dimensions;
    let fine = raster_cells("top-plate-fholes-0p1mm.tsv", 0.1);
    let coarse = raster_cells("top-plate-fholes.tsv", 0.5);
    assert_eq!(fine.keys().collect::<Vec<_>>(), [&1, &2]);
    assert_eq!(coarse.keys().collect::<Vec<_>>(), [&1, &2]);
    for hole in [1u32, 2] {
        assert_eq!(coarsened(&fine[&hole], 13), coarse[&hole], "hole {hole}: not the half-open coarsening");
        let row = dims
            .get("top-plate", &format!("f-hole-{hole}"), "area")
            .expect("area row");
        let Sigma::Value(sigma) = row.sigma else { panic!("hole {hole}: area sigma unstated") };
        let area = coarse[&hole].len() as f64 * 0.25e-6;
        assert!(((area - row.value_si) / row.value_si).abs() < 1e-5, "hole {hole}: {area} against {}", row.value_si);
        let edge_as_wood = coarsened(&fine[&hole], 25);
        assert_ne!(edge_as_wood, coarse[&hole], "hole {hole}: the planted rule must differ");
        let short = row.value_si - edge_as_wood.len() as f64 * 0.25e-6;
        assert!(short > sigma, "hole {hole}: edge-as-wood falls short by {short} m2, sigma {sigma}");
        // sigma is the rule's error against the 0.1 mm area, written by the extraction
        assert!((sigma - wall(&hole.to_string(), "area-0p5mm-rule-error") * 1e-6).abs() < 1e-11);
        assert!((wall(&hole.to_string(), "area-0p1mm") * 1e-6 - fine[&hole].len() as f64 * 0.01e-6).abs() < 1e-11);
    }
}

#[test]
fn f_hole_bores_equal_the_wall_table() {
    // dimensions.tsv's bores are the wall table's: each hole's measured through-resistance times its raster area, and
    // the one bore of both that keeps their parallel bore conductance.
    let dims = load().dimensions;
    let coarse = raster_cells("top-plate-fholes.tsv", 0.5);
    let mut parallel = 0.0;
    for hole in [1u32, 2] {
        let h = hole.to_string();
        let bore = value(&dims, "top-plate", &format!("f-hole-{hole}"), "bore");
        assert!((bore - wall(&h, "bore") * 1e-3).abs() < 1e-9, "hole {hole} bore {bore}");
        let area_mm2 = coarse[&hole].len() as f64 * 0.25;
        assert!((wall(&h, "bore") - wall(&h, "through-resistance") * area_mm2).abs() < 1e-4);
        parallel += 1.0 / wall(&h, "through-resistance");
    }
    let total_mm2 = (coarse[&1].len() + coarse[&2].len()) as f64 * 0.25;
    let pair = value(&dims, "top-plate", "f-holes", "bore-pair");
    assert!((pair - wall("0", "bore-pair") * 1e-3).abs() < 1e-9);
    assert!((wall("0", "bore-pair") - total_mm2 / parallel).abs() < 1e-4);
}
