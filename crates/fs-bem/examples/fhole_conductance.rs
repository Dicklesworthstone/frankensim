//! The il Cannone f-holes' potential-flow conductance for the body tables' neck (`fs_bem::plate_aperture`), written
//! as a `frankensim.specimen-fhole-conductance.v1` table.
//!
//!   fhole_conductance RASTER.tsv WALLS.tsv > top-plate-fhole-conductance.tsv
//!
//! RASTER is the 0.5 mm f-hole raster (`frankensim.specimen-raster-runs.v1`, top-plate-fholes.tsv) and WALLS its
//! wall table (`frankensim.specimen-fhole-walls.v1`, top-plate-fhole-walls.tsv): the `bore` row of each hole is the
//! straight bore with that hole's measured through-resistance, and the `bore-pair` row the one bore of both holes
//! that keeps their parallel bore conductance (tools/geometry/violin-cannone-esrf/fhole_openings.py).
//!
//! Five ladders, each at subdivisions 1, 2, 4 and 8 of the raster cell with the Richardson limit of the last three;
//! a ladder that is not monotone, an observed order below 0.5, or top and bottom fluxes that differ by more than 1e-8
//! refuse:
//! - `pair`: both holes in one raster through the pair bore (the neck of the shell tables);
//! - `hole1-at-pair-bore`, `hole2-at-pair-bore`: each hole alone through the pair bore (`pair` over their sum is the
//!   holes' mutual factor);
//! - `hole1-at-own-bore`, `hole2-at-own-bore`: each hole alone through its own bore.
//!
//! The plate is flat, rigid and infinite, with a half-space of air on each side (the model of
//! `fs_bem::plate_aperture`): the arching, the finite body and the cavity side of each hole are not represented.
//! Nothing is fitted.
use std::{collections::BTreeMap, error::Error, fs};

use fs_bem::plate_aperture::{ApertureRaster, PlateApertureOptions, plate_aperture_conductance};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

const CELL_MM: f64 = 0.5;

/// Open cells of each hole of a raster-runs table, keyed by hole number. A value off the 0.5 mm grid refuses.
fn read_raster(text: &str) -> Result<BTreeMap<u32, Vec<(i64, i64)>>> {
    let mut lines = text.lines();
    if lines.next() != Some("frankensim.specimen-raster-runs.v1")
        || lines.next() != Some("hole\tz_mm\tx_start_mm\tx_end_mm")
    {
        return Err("not a frankensim.specimen-raster-runs.v1 hole table".into());
    }
    let grid = |v: f64| -> Result<i64> {
        let k = (v / CELL_MM).round();
        if (k * CELL_MM - v).abs() > 1e-9 {
            return Err(format!("{v} mm is off the 0.5 mm raster").into());
        }
        Ok(k as i64)
    };
    let mut holes: BTreeMap<u32, Vec<(i64, i64)>> = BTreeMap::new();
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 4 {
            return Err(format!("raster run needs 4 fields: {line}").into());
        }
        let hole: u32 = f[0].parse()?;
        let (z, a, b): (f64, f64, f64) = (f[1].parse()?, f[2].parse()?, f[3].parse()?);
        let row = grid(z - 0.5 * CELL_MM)?;
        let cells = holes.entry(hole).or_default();
        for column in grid(a)?..grid(b)? {
            cells.push((column, row));
        }
    }
    if holes.keys().copied().collect::<Vec<_>>() != [1, 2] {
        return Err("the raster must hold holes 1 and 2".into());
    }
    Ok(holes)
}

/// The value [mm] of a wall-table row, as written.
fn wall_value(text: &str, hole: &str, quantity: &str) -> Result<String> {
    let mut lines = text.lines();
    if lines.next() != Some("frankensim.specimen-fhole-walls.v1")
        || lines.next() != Some("hole\tquantity\tvalue\tunit\tnote")
    {
        return Err("not a frankensim.specimen-fhole-walls.v1 table".into());
    }
    let rows: Vec<Vec<&str>> = lines
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f.len() >= 4 && f[0] == hole && f[1] == quantity)
        .collect();
    match rows.as_slice() {
        [row] if row[3] == "mm" => Ok(row[2].to_owned()),
        [_] => Err(format!("wall row {hole} {quantity} is not in mm").into()),
        _ => Err(format!("wall table needs exactly one row {hole} {quantity}").into()),
    }
}

/// One ladder: conductances [m] at subdivisions 1, 2, 4, 8, the observed orders of the first and last triples, the
/// Richardson limit and its extrapolation step [m].
fn ladder(raster: &ApertureRaster, thickness_m: f64) -> Result<([f64; 4], f64, f64, f64, f64)> {
    let mut k = [0.0; 4];
    for (slot, subdivision) in k.iter_mut().zip([1, 2, 4, 8]) {
        let options = PlateApertureOptions {
            thickness: thickness_m,
            subdivision,
            tolerance: 1e-10,
            max_cells: 1 << 22,
            max_fft_points: 1 << 25,
            max_voxels: 1 << 25,
            workers: 16,
            ..PlateApertureOptions::default()
        };
        let started = std::time::Instant::now();
        let r = plate_aperture_conductance(raster, &options)
            .map_err(|e| format!("refused at subdivision {subdivision}: {e:?}"))?;
        let flux = (r.top_flux - r.conductance).abs() / r.conductance;
        if flux > 1e-8 {
            return Err(format!("top and bottom fluxes differ by {flux:e}").into());
        }
        eprintln!(
            "  subdivision {subdivision}: {} cells, {} layers, K {:.6} mm, residual {:.2e}, {:.1} s",
            r.cells,
            r.layers,
            r.conductance * 1e3,
            r.residual,
            started.elapsed().as_secs_f64()
        );
        *slot = r.conductance;
    }
    let order = |a: f64, b: f64, c: f64| ((b - a) / (c - b)).log2();
    let (p_first, p) = (order(k[0], k[1], k[2]), order(k[1], k[2], k[3]));
    if !k.windows(2).all(|w| w[1] > w[0]) || !(p_first >= 0.5 && p >= 0.5) {
        return Err(format!("ladder not converging: K {k:?}, orders {p_first} {p}").into());
    }
    let step = (k[3] - k[2]) / (2f64.powf(p) - 1.0);
    Ok((k, p_first, p, k[3] + step, step))
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() != 2 {
        return Err("usage: fhole_conductance RASTER.tsv WALLS.tsv > top-plate-fhole-conductance.tsv".into());
    }
    let holes = read_raster(&fs::read_to_string(&args[0])?)?;
    let walls = fs::read_to_string(&args[1])?;
    let pair_bore = wall_value(&walls, "0", "bore-pair")?;
    let own_bore = [wall_value(&walls, "1", "bore")?, wall_value(&walls, "2", "bore")?];
    let both: Vec<(i64, i64)> = holes.values().flatten().copied().collect();
    let configs: [(&str, Vec<(i64, i64)>, &str); 5] = [
        ("pair", both, &pair_bore),
        ("hole1-at-pair-bore", holes[&1].clone(), &pair_bore),
        ("hole2-at-pair-bore", holes[&2].clone(), &pair_bore),
        ("hole1-at-own-bore", holes[&1].clone(), &own_bore[0]),
        ("hole2-at-own-bore", holes[&2].clone(), &own_bore[1]),
    ];
    println!("frankensim.specimen-fhole-conductance.v1");
    println!("config\tcells\tthickness_mm\tK_m1_mm\tK_m2_mm\tK_m4_mm\tK_m8_mm\torder_first\torder_last\tK_mm\tstep_mm");
    for (name, open, bore) in configs {
        let thickness_mm: f64 = bore.parse()?;
        eprintln!("{name}: {} cells through {bore} mm", open.len());
        let cells = open.len();
        let raster = ApertureRaster {
            cell: CELL_MM * 1e-3,
            open,
        };
        let (k, p_first, p, limit, step) = ladder(&raster, thickness_mm * 1e-3)?;
        println!(
            "{name}\t{cells}\t{bore}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{p_first:.4}\t{p:.4}\t{:.6}\t{:.6}",
            k[0] * 1e3,
            k[1] * 1e3,
            k[2] * 1e3,
            k[3] * 1e3,
            limit * 1e3,
            step * 1e3
        );
    }
    Ok(())
}
