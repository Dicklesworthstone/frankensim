//! Potential-flow conductance of plate apertures (`fs_bem::plate_aperture`).
//!
//! Oracles: Rayleigh's thin disc `K = 2a` and thin ellipse `K = 2a AGM(1,
//! b/a)`, and the flanged long duct, whose end correction tends to 0.8216 a
//! at each end (Norris and Sheng 1989). Symmetry: rotations, a mirror and a
//! translation of the raster. Measurement: the il Cannone f-holes from the
//! CT raster, thin (pexfp.366) and through their measured bores, where the
//! specimen's conductance table must come out digit for digit (pexfp.377).

use std::{collections::BTreeMap, time::Instant};

use fs_bem::plate_aperture::{
    ApertureRaster, PlateApertureConductance, PlateApertureOptions, plate_aperture_conductance,
};

const SPECIMEN: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/"
);

fn options(thickness: f64, subdivision: usize) -> PlateApertureOptions {
    PlateApertureOptions {
        thickness,
        subdivision,
        tolerance: 1e-10,
        max_cells: 1 << 22,
        max_fft_points: 1 << 25,
        max_voxels: 1 << 25,
        workers: 16,
        ..PlateApertureOptions::default()
    }
}

fn solve(raster: &ApertureRaster, o: &PlateApertureOptions) -> (PlateApertureConductance, f64) {
    let start = Instant::now();
    let result = plate_aperture_conductance(raster, o).unwrap();
    let flux = (result.top_flux - result.conductance).abs() / result.conductance;
    assert!(flux <= 1e-8, "top and bottom flux differ by {flux:e}");
    (result, start.elapsed().as_secs_f64())
}

/// Cells whose centre lies inside the disc of radius `rpx` cells.
fn disc(rpx: i64, cell: f64) -> ApertureRaster {
    let mut open = Vec::new();
    for i in -rpx - 1..=rpx {
        for j in -rpx - 1..=rpx {
            let (x, y) = (i as f64 + 0.5, j as f64 + 0.5);
            if x * x + y * y <= (rpx * rpx) as f64 {
                open.push((i, j));
            }
        }
    }
    ApertureRaster { cell, open }
}

fn agm(a: f64, b: f64) -> f64 {
    let (mut g1, mut g2) = (a, b);
    for _ in 0..40 {
        (g1, g2) = (0.5 * (g1 + g2), (g1 * g2).sqrt());
    }
    g1
}

/// `log2` of the ratio of successive differences, and the Richardson limit
/// from the last three rungs (refinement ratio 2).
fn richardson(k: &[f64]) -> (f64, f64) {
    let n = k.len();
    let (d1, d2) = (k[n - 2] - k[n - 3], k[n - 1] - k[n - 2]);
    let p = (d1 / d2).log2();
    (p, k[n - 1] + d2 / (2f64.powf(p) - 1.0))
}

#[test]
fn rayleigh_thin_disc_ladder() {
    // A thin circular hole of radius a has K = 2a. The raster is the
    // discretization: the disc of radius 1 m in cells of 1/rpx.
    let mut errors = Vec::new();
    for rpx in [4, 8, 16, 32, 64] {
        let raster = disc(rpx, 1.0 / rpx as f64);
        let (r, seconds) = solve(&raster, &options(0.0, 1));
        let a_eq = (raster.open.len() as f64 / std::f64::consts::PI).sqrt() / rpx as f64;
        let error = r.conductance / (2.0 * a_eq) - 1.0;
        println!(
            "{{\"test\":\"rayleigh_disc\",\"rpx\":{rpx},\"cells\":{},\"K_over_2a_eq\":{:.8},\"K_over_2R\":{:.8},\"iterations\":{},\"residual\":{:.2e},\"seconds\":{seconds:.2}}}",
            r.cells,
            r.conductance / (2.0 * a_eq),
            r.conductance / 2.0,
            r.iterations,
            r.residual
        );
        errors.push(error.abs());
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order >= 0.8,
            "disc ladder order {order:.3} below first order: {errors:?}"
        );
    }
    assert!(*errors.last().unwrap() <= 2.5e-3, "{errors:?}");
}

#[test]
fn rayleigh_thin_ellipse_at_the_ct_pixel() {
    // An ellipse as slender as an f-hole's equivalent (42 x 3.1 mm), rastered
    // at the CT's 0.5 mm. The subdivision ladder converges to the raster's
    // own conductance; the gate is how far that sits from the smooth ellipse.
    let (a, b, cell) = (42e-3, 3.1e-3, 0.5e-3);
    let exact = 2.0 * a * agm(1.0, b / a);
    let (na, nb) = ((a / cell) as i64 + 2, (b / cell) as i64 + 2);
    let mut open = Vec::new();
    for i in -na..na {
        for j in -nb..nb {
            let (x, y) = ((i as f64 + 0.5) * cell / a, (j as f64 + 0.5) * cell / b);
            if x * x + y * y <= 1.0 {
                open.push((i, j));
            }
        }
    }
    let raster = ApertureRaster { cell, open };
    let mut k = Vec::new();
    for m in [1, 2, 4, 8] {
        let (r, seconds) = solve(&raster, &options(0.0, m));
        println!(
            "{{\"test\":\"rayleigh_ellipse\",\"m\":{m},\"cells\":{},\"K_over_exact\":{:.8},\"iterations\":{},\"residual\":{:.2e},\"seconds\":{seconds:.2}}}",
            r.cells,
            r.conductance / exact,
            r.iterations,
            r.residual
        );
        k.push(r.conductance / exact);
    }
    let (p_first, _) = richardson(&k[..3]);
    let (p, limit) = richardson(&k);
    println!(
        "{{\"test\":\"rayleigh_ellipse\",\"order\":[{p_first:.3},{p:.3}],\"limit_over_exact\":{limit:.6}}}"
    );
    assert!(p_first >= 0.8 && p >= 0.8, "orders {p_first} {p}: {k:?}");
    assert!(
        (limit - 1.0).abs() <= 5e-3,
        "the 0.5 mm raster sits {limit} from the ellipse"
    );
    assert!((k[3] - 1.0).abs() <= 5e-3, "{k:?}");
}

#[test]
fn flanged_long_duct_end_correction() {
    // A circular duct four radii long: K = pi a^2 / (t + 2 delta a), and the
    // ends no longer interact (the first evanescent mode has decayed by
    // e^-15), so delta is the flanged end correction 0.8216.
    let rpx = 8;
    let raster = disc(rpx, 1.0 / rpx as f64);
    let a_eq = (raster.open.len() as f64 / std::f64::consts::PI).sqrt() / rpx as f64;
    let t = 4.0;
    let mut delta = Vec::new();
    for m in [1, 2, 4] {
        let (r, seconds) = solve(&raster, &options(t, m));
        let d = (std::f64::consts::PI * a_eq * a_eq / r.conductance - t) / (2.0 * a_eq);
        println!(
            "{{\"test\":\"flanged_duct\",\"m\":{m},\"cells\":{},\"layers\":{},\"end_correction\":{d:.6},\"iterations\":{},\"inner_iterations\":{},\"residual\":{:.2e},\"seconds\":{seconds:.2}}}",
            r.cells, r.layers, r.iterations, r.inner_iterations, r.residual
        );
        delta.push(d);
    }
    let (p, limit) = richardson(&delta);
    println!("{{\"test\":\"flanged_duct\",\"order\":{p:.3},\"limit\":{limit:.6}}}");
    assert!(p >= 0.8, "end-correction ladder order {p}: {delta:?}");
    assert!((limit - 0.8216).abs() <= 0.01 * 0.8216, "{limit} {delta:?}");
}

#[test]
fn rotations_mirror_and_translation_keep_the_conductance() {
    // An asymmetric polyomino through a plate 1.5 cells thick.
    let base: Vec<(i64, i64)> = vec![
        (0, 0),
        (1, 0),
        (2, 0),
        (2, 1),
        (2, 2),
        (3, 2),
        (4, 2),
        (4, 3),
        (0, 1),
        (0, 2),
        (7, 1),
    ];
    let transforms: [fn((i64, i64)) -> (i64, i64); 5] = [
        |(i, j)| (i, j),
        |(i, j)| (j, -i),
        |(i, j)| (-i, -j),
        |(i, j)| (-i, j),
        |(i, j)| (i + 1000, j - 77),
    ];
    for thickness in [0.0, 1.5e-3] {
        let k: Vec<f64> = transforms
            .iter()
            .map(|f| {
                let raster = ApertureRaster {
                    cell: 1e-3,
                    open: base.iter().map(|&c| f(c)).collect(),
                };
                let mut o = options(thickness, 2);
                o.tolerance = 1e-12;
                solve(&raster, &o).0.conductance
            })
            .collect();
        for v in &k[1..] {
            assert!(((v - k[0]) / k[0]).abs() <= 1e-9, "t {thickness}: {k:?}");
        }
    }
}

/// The open cells of each f-hole: rows of 0.5 mm cells, `z_mm` the row
/// centre, `[x_start_mm, x_end_mm)` the run.
fn read_f_holes() -> [Vec<(i64, i64)>; 2] {
    let path = format!("{SPECIMEN}top-plate-fholes.tsv");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("frankensim.specimen-raster-runs.v1"));
    assert_eq!(lines.next(), Some("hole\tz_mm\tx_start_mm\tx_end_mm"));
    let mut holes = [Vec::new(), Vec::new()];
    let grid = |v: f64| {
        let k = (v / 0.5).round();
        assert!((k * 0.5 - v).abs() < 1e-9, "{v} is off the 0.5 mm grid");
        k as i64
    };
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let f: Vec<f64> = line.split('\t').map(|v| v.parse().unwrap()).collect();
        let j = grid(f[1] - 0.25);
        for i in grid(f[2])..grid(f[3]) {
            holes[f[0] as usize - 1].push((i, j));
        }
    }
    holes
}

/// A specimen table after its schema line: the header's columns and the rows,
/// each as wide as the header.
fn read_specimen_table(name: &str, schema: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let text = std::fs::read_to_string(format!("{SPECIMEN}{name}")).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(schema), "{name} schema");
    let header: Vec<String> = lines.next().unwrap().split('\t').map(str::to_owned).collect();
    let rows: Vec<Vec<String>> = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    assert!(rows.iter().all(|r| r.len() == header.len()), "{name} row width");
    (header, rows)
}

/// The conductance table, keyed by configuration, each row by column name.
fn read_conductance_table() -> BTreeMap<String, BTreeMap<String, String>> {
    let (header, rows) = read_specimen_table(
        "top-plate-fhole-conductance.tsv",
        "frankensim.specimen-fhole-conductance.v1",
    );
    rows.into_iter()
        .map(|r| (r[0].clone(), header.iter().cloned().zip(r).collect()))
        .collect()
}

/// A bore of the wall table [mm], as written: `("0", "bore-pair")` or
/// `("1", "bore")`.
fn wall_bore(hole: &str, quantity: &str) -> String {
    let (_, rows) = read_specimen_table(
        "top-plate-fhole-walls.tsv",
        "frankensim.specimen-fhole-walls.v1",
    );
    let found: Vec<&Vec<String>> = rows
        .iter()
        .filter(|r| r[0] == hole && r[1] == quantity)
        .collect();
    assert_eq!(found.len(), 1, "wall row {hole} {quantity}");
    assert_eq!(found[0][3], "mm");
    found[0][2].clone()
}

/// The h8 body table's `cavity` neck radius and length [mm].
fn read_h8_neck() -> (f64, f64) {
    let text = std::fs::read_to_string(format!("{SPECIMEN}shell-assembly-h8.tsv")).unwrap();
    let row: Vec<f64> = text
        .lines()
        .find(|l| l.starts_with("cavity\t"))
        .unwrap()
        .split('\t')
        .skip(1)
        .map(|v| v.parse().unwrap())
        .collect();
    (row[1] * 1e3, row[2] * 1e3)
}

/// One ladder at subdivisions 1, 2, 4, 8 with its JSON lines: the conductances
/// [m], the observed orders of the first and last triples, the Richardson limit
/// and its extrapolation step [m], in the arithmetic of the example that writes
/// the conductance table; a ladder that is not monotone or has an order below
/// 0.5 fails.
fn cannone_ladder(
    name: &str,
    cells: &[(i64, i64)],
    thickness_mm: f64,
) -> (Vec<f64>, f64, f64, f64, f64) {
    let raster = ApertureRaster {
        cell: 0.5e-3,
        open: cells.to_vec(),
    };
    let mut k = Vec::new();
    for m in [1, 2, 4, 8] {
        let (r, seconds) = solve(&raster, &options(thickness_mm * 1e-3, m));
        println!(
            "{{\"test\":\"cannone_f_holes\",\"holes\":\"{name}\",\"thickness_mm\":{thickness_mm},\"m\":{m},\"h_mm\":{:.4},\"cells\":{},\"layers\":{},\"fft_points\":{},\"K_mm\":{:.6},\"iterations\":{},\"inner_iterations\":{},\"residual\":{:.2e},\"seconds\":{seconds:.1}}}",
            r.cell_side * 1e3,
            r.cells,
            r.layers,
            r.fft_points,
            r.conductance * 1e3,
            r.iterations,
            r.inner_iterations,
            r.residual
        );
        k.push(r.conductance);
    }
    let (p_first, _) = richardson(&k[..3]);
    let (p, limit) = richardson(&k);
    let step = (k[3] - k[2]) / (2f64.powf(p) - 1.0);
    println!(
        "{{\"test\":\"cannone_f_holes\",\"holes\":\"{name}\",\"thickness_mm\":{thickness_mm},\"order\":[{p_first:.3},{p:.3}],\"K_limit_mm\":{:.4},\"K_limit_m\":{limit:.6e}}}",
        limit * 1e3
    );
    assert!(
        k.windows(2).all(|w| w[1] > w[0]),
        "{name} t {thickness_mm}: not monotone {k:?}"
    );
    assert!(
        p_first >= 0.5 && p >= 0.5,
        "{name} t {thickness_mm}: orders {p_first} {p}: {k:?}"
    );
    (k, p_first, p, limit, step)
}

/// The equivalent ellipse of the body example on the raster cells: the area,
/// the semi-axis from the principal second moment, Rayleigh's K0 [mm].
fn equivalent_ellipse_mm(cells: &[(i64, i64)]) -> (f64, f64) {
    let h = 0.5;
    let n = cells.len() as f64;
    let centre = |c: &(i64, i64)| ((c.0 as f64 + 0.5) * h, (c.1 as f64 + 0.5) * h);
    let (mx, mz) = cells
        .iter()
        .map(centre)
        .fold((0.0, 0.0), |s, p| (s.0 + p.0 / n, s.1 + p.1 / n));
    let (mut sxx, mut szz, mut sxz) = (h * h / 12.0, h * h / 12.0, 0.0);
    for c in cells {
        let (x, z) = centre(c);
        sxx += (x - mx) * (x - mx) / n;
        szz += (z - mz) * (z - mz) / n;
        sxz += (x - mx) * (z - mz) / n;
    }
    let lambda = 0.5 * (sxx + szz) + (0.25 * (sxx - szz) * (sxx - szz) + sxz * sxz).sqrt();
    let area = n * h * h;
    let a = 2.0 * lambda.sqrt();
    let b = area / (std::f64::consts::PI * a);
    (area, 2.0 * a * agm(1.0, b / a))
}

#[test]
fn cannone_f_hole_conductance_ladder() {
    // Measurement, not a band check: the conductance of the CT f-holes, each
    // hole alone and both together, on a subdivision ladder of the 0.5 mm
    // raster, as a thin plate and against the equivalent ellipse.
    let holes = read_f_holes();
    let both: Vec<(i64, i64)> = holes[0].iter().chain(&holes[1]).copied().collect();
    let configs: [(&str, &Vec<(i64, i64)>); 3] =
        [("hole1", &holes[0]), ("hole2", &holes[1]), ("both", &both)];
    let mut k0s = Vec::new();
    for (name, cells) in &configs[..2] {
        let (area, k0) = equivalent_ellipse_mm(cells);
        println!(
            "{{\"test\":\"cannone_f_holes\",\"holes\":\"{name}\",\"cells\":{},\"area_mm2\":{area:.2},\"ellipse_K0_mm\":{k0:.4}}}",
            cells.len()
        );
        k0s.push(k0);
    }
    let limits: Vec<f64> = configs
        .iter()
        .map(|(name, cells)| cannone_ladder(name, cells, 0.0).3 * 1e3)
        .collect();
    let reference: f64 = k0s.iter().sum();
    println!(
        "{{\"test\":\"cannone_f_holes\",\"thickness_mm\":0,\"K_both_mm\":{:.4},\"K_hole1_plus_hole2_mm\":{:.4},\"mutual_factor\":{:.5},\"ellipse_sum_mm\":{reference:.4},\"ratio_to_ellipse_sum\":{:.5}}}",
        limits[2],
        limits[0] + limits[1],
        limits[2] / (limits[0] + limits[1]),
        limits[2] / reference
    );
}

#[test]
fn cannone_f_hole_conductance_table_reproduces() {
    // top-plate-fhole-conductance.tsv (crates/fs-bem/examples/fhole_conductance.rs)
    // is the solver's ladder on the raster through the wall table's bores,
    // digit for digit, and the h8 body table's neck carries its pair row: the
    // flanged circular neck of that conductance through the pair bore.
    let holes = read_f_holes();
    let both: Vec<(i64, i64)> = holes[0].iter().chain(&holes[1]).copied().collect();
    let table = read_conductance_table();
    let (pair_bore, own) = (
        wall_bore("0", "bore-pair"),
        [wall_bore("1", "bore"), wall_bore("2", "bore")],
    );
    let configs: [(&str, &Vec<(i64, i64)>, &str); 5] = [
        ("pair", &both, &pair_bore),
        ("hole1-at-pair-bore", &holes[0], &pair_bore),
        ("hole2-at-pair-bore", &holes[1], &pair_bore),
        ("hole1-at-own-bore", &holes[0], &own[0]),
        ("hole2-at-own-bore", &holes[1], &own[1]),
    ];
    assert_eq!(
        table.keys().collect::<Vec<_>>(),
        {
            let mut names: Vec<&str> = configs.iter().map(|c| c.0).collect();
            names.sort_unstable();
            names
        },
        "the table's configurations"
    );
    for (name, cells, bore) in configs {
        let row = &table[name];
        assert_eq!(row["cells"], cells.len().to_string(), "{name} cells");
        assert_eq!(row["thickness_mm"], bore, "{name} bore");
        let (k, p_first, p, limit, step) = cannone_ladder(name, cells, bore.parse().unwrap());
        let computed = [
            format!("{:.6}", k[0] * 1e3),
            format!("{:.6}", k[1] * 1e3),
            format!("{:.6}", k[2] * 1e3),
            format!("{:.6}", k[3] * 1e3),
            format!("{p_first:.4}"),
            format!("{p:.4}"),
            format!("{:.6}", limit * 1e3),
            format!("{:.6}", step * 1e3),
        ];
        let columns = [
            "K_m1_mm", "K_m2_mm", "K_m4_mm", "K_m8_mm", "order_first", "order_last", "K_mm", "step_mm",
        ];
        for (column, value) in columns.iter().zip(&computed) {
            assert_eq!(&row[*column], value, "{name} {column}");
        }
    }
    let k_pair: f64 = table["pair"]["K_mm"].parse().unwrap();
    let (radius, length) = read_h8_neck();
    let pair_bore_mm: f64 = pair_bore.parse().unwrap();
    assert!((length - pair_bore_mm).abs() < 1e-9, "h8 neck length {length} mm");
    let k_neck = std::f64::consts::PI * radius * radius / (length + 2.0 * 0.8216 * radius);
    assert!(
        ((k_neck - k_pair) / k_pair).abs() < 1e-6,
        "h8 neck K {k_neck} mm against the table's {k_pair} mm"
    );
}
