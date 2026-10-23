//! Generate the soundboard shipped with the live engine (run NATIVELY).
//!
//! It is the offline reference's own computation: `steinway_d::build_resolved(16, 0.03, true)` ->
//! `board_geometry::BoardGeometry::prepare_reduced` (the binary's `--board-geometry
//! board-rt16.fsb --board-band-hz BAND --board-reduction SPEC`). The result is written in
//! the layout documented at `Description::board_bytes` and embedded by the
//! library with `include_bytes!`. Re-run it whenever the reference board
//! modules change; the output path must not exist.
//!
//! usage: gen_board OUT.bin [source_band_hz=2400] [spec=pod,32,16,16,8000] [mesh=16,0.03] [damping=authored]
//! `spec` is the reference's --board-reduction. The live boards take N/2 exact low modes and the
//! POD of every key's static and 16 harmonic bridge responses from the first unprotected mode to
//! 8 kHz (`pod,N,N/2,16,8000`). Above 2 kHz they follow the complete slice about 3x more closely
//! in radiated pressure and 7-17x in bridge conductance than the pivoted
//! `N,N/2,200,400,800,1600,2300` boards did, and improve steadily with N
//! (examples/board_admittance.rs; grand_piano/BOARD_REDUCTION.md).
//! `mesh` is DIVISIONS,RIB_NORMAL_STEP_M of `steinway_d::build_resolved` (the board the graded
//! offline renders use: piano_geometry --divisions 16 --rib-normal-step-m 0.03), or `preset` for
//! the unrefined `steinway_d::build(8)`.
//! `damping` is `authored` (the shipped `damping,0.015`) or `material`: the board declares
//! `damping,material` and every mode takes its energy-weighted loss from the loss packs
//! (pexfp.373; claims measured at other frequencies taken by declared extrapolation, sugar pine
//! loss unknown and omitted as a printed lower bound). A material board's loss receipt is printed
//! per source mode and per band against the authored 0.015.
use fs_live_piano::refmods::{board_geometry, geometry, steinway_d};
use fs_live_piano::Description;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(out) = args.get(1) else {
        eprintln!("usage: gen_board OUT.bin [source_band_hz] [spec] [mesh]");
        std::process::exit(2);
    };
    let band: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2400.0);
    let spec = args.get(3).cloned().unwrap_or_else(|| "pod,32,16,16,8000".to_owned());
    let mesh = args.get(4).map_or("16,0.03", String::as_str);
    let material = match args.get(5).map_or("authored", String::as_str) {
        "authored" => false,
        "material" => true,
        other => { eprintln!("damping must be authored or material, not {other}"); std::process::exit(2); }
    };
    let policy = steinway_d::LossPolicy { extrapolate_frequency: true, sugar_pine_loss_unknown: true };
    let materials = steinway_d::reference_materials().expect("board materials");
    let preset = match (mesh, material) {
        ("preset", false) => steinway_d::build(8),
        ("preset", true) => steinway_d::build_with_material_loss(8, true, &materials, policy),
        (mesh, _) => {
            let (divisions, step) = mesh.split_once(',').expect("mesh is DIVISIONS,STEP or preset");
            let (divisions, step) = (divisions.parse().expect("divisions"), step.parse().expect("rib-normal step"));
            if material { steinway_d::build_resolved_with_material_loss(divisions, step, true, &materials, policy) }
            else { steinway_d::build_resolved(divisions, step, true) }
        }
    }.expect("board geometry");
    let courses = geometry::read_scale(fs_live_piano::kernel::SCALE_CSV).expect("embedded scale");
    let start = std::time::Instant::now();
    let desc = Description::with_reduction_on(&preset, courses.clone(), band, &spec).expect("reduced board");
    println!("prepared in {:.3} s: {} coordinates ({} exact low modes), {} surface points, {:.6} kg, {:.6} m^2",
        start.elapsed().as_secs_f64(), desc.board.len(), desc.exact_low_modes, desc.surface.len(), desc.board_mass_kg, desc.board_area_m2);
    println!("{}", desc.board_provenance);
    for (i, b) in desc.board.iter().enumerate() {
        println!("coordinate {i:2}: {:10.4} Hz{}", b.frequency_hz, if i < desc.exact_low_modes { " (exact FE mode)" } else { " (projected-pencil Ritz)" });
    }
    if material {
        let options = board_geometry::ritz::RitzOptions::parse(&spec).expect("reduction spec");
        let keys: Vec<u8> = courses.iter().map(|c| c.midi).collect();
        let (_, receipt) = board_geometry::BoardGeometry::read(&preset.geometry).expect("board geometry")
            .prepare_reduced_with_loss_receipt(&keys, band, false, false, false, false, &options).expect("material loss");
        let receipt = receipt.expect("a damping,material board returns its loss receipt");
        print!("{}", receipt.tsv());
        for b in receipt.band_summary(&[0.0, 200.0, 400.0, 800.0, 1600.0, band]) {
            println!("band {:.0}-{:.0} Hz: {} source modes, damping ratio mean {:.5} (min {:.5}, max {:.5}) = {:.3} x the authored 0.015",
                b.band_hz.0, b.band_hz.1, b.modes, b.mean_ratio, b.min_ratio, b.max_ratio, b.mean_ratio / 0.015);
        }
        for caveat in &receipt.caveats { println!("caveat: {caveat}"); }
        for (i, b) in desc.board.iter().enumerate() {
            println!("coordinate {i:2}: damping ratio {:.5} = {:.3} x the authored 0.015", b.damping_ratio, b.damping_ratio / 0.015);
        }
    }
    let bytes = desc.board_bytes();
    // Round trip before publishing: the loader must reproduce what was written.
    let back = Description::from_board_bytes(courses, &bytes).expect("round trip");
    assert_eq!(back.board.len(), desc.board.len());
    for (a, b) in back.board.iter().zip(&desc.board) {
        assert_eq!(a.frequency_hz.to_bits(), b.frequency_hz.to_bits());
        assert_eq!(a.bridge, b.bridge);
    }
    assert_eq!(back.damping, desc.damping);
    let worst = back.surface.iter().zip(&desc.surface).flat_map(|(a, b)| a.mode_shape.iter().zip(&b.mode_shape))
        .map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
    let scale = desc.surface.iter().flat_map(|p| p.mode_shape.iter()).fold(0.0_f64, |m, x| m.max(x.abs()));
    println!("single-precision surface shapes: max |error| {worst:.3e} against max |shape| {scale:.3e}");
    std::fs::OpenOptions::new().write(true).create_new(true).open(out)
        .and_then(|mut f| f.write_all(&bytes)).unwrap_or_else(|e| { eprintln!("{out}: {e}"); std::process::exit(1); });
    println!("wrote {} bytes to {out}", bytes.len());
}
