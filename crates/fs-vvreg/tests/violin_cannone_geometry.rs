//! M8-W3c violin geometry specimen (IRP bead `frankensim-irp-v5-instrument-realism-pexfp.199`): Paganini's
//! violin "il Cannone" (Guarneri del Gesu, 1743) from the ESRF synchrotron CT data collection
//! 10.15151/ESRF-DC-2413577471 (CC BY 4.0).
//!
//! A data-form and consistency check, never a physical claim. It checks that the five specimen tables parse
//! and agree with each other, that every derived scalar equals the data table it was computed from, that the
//! data tables are well formed, and that the map points lie inside the plate outlines. The specimen is NOT the
//! UIowa violin of the recordings, and the boundary stays recorded in `unknowns.tsv`.

use fs_vvreg::specimen::{
    Dimensions, KeyedTable, SPECIMEN_GEOMETRY_ROOT, Specimen, SpecimenDescription, SpecimenTable,
};
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
