//! Generate the soundboard shipped with the live engine (run NATIVELY).
//!
//! It is the offline reference's own computation: `steinway_d::build_resolved(16, 0.03, true)` ->
//! `board_geometry::BoardGeometry::prepare_reduced` (the binary's `--board-geometry
//! board-rt16.fsb --board-band-hz BAND --board-reduction SPEC`). The result is written in
//! the layout documented at `Description::board_bytes` and embedded by the
//! library with `include_bytes!`. Re-run it whenever the reference board
//! modules change; the output path must not exist.
//!
//! usage: gen_board OUT.bin [source_band_hz=2400] [spec=32,16,200,400,800,1600,2300] [mesh=16,0.03]
//! `mesh` is DIVISIONS,RIB_NORMAL_STEP_M of `steinway_d::build_resolved` (the board the graded
//! offline renders use: piano_geometry --divisions 16 --rib-normal-step-m 0.03), or `preset` for
//! the unrefined `steinway_d::build(8)`.
use fs_live_piano::refmods::{geometry, steinway_d};
use fs_live_piano::Description;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(out) = args.get(1) else {
        eprintln!("usage: gen_board OUT.bin [source_band_hz] [spec] [mesh]");
        std::process::exit(2);
    };
    let band: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2400.0);
    let spec = args.get(3).cloned().unwrap_or_else(|| "32,16,200,400,800,1600,2300".to_owned());
    let mesh = args.get(4).map_or("16,0.03", String::as_str);
    let preset = if mesh == "preset" {
        steinway_d::build(8)
    } else {
        let (divisions, step) = mesh.split_once(',').expect("mesh is DIVISIONS,STEP or preset");
        steinway_d::build_resolved(divisions.parse().expect("divisions"), step.parse().expect("rib-normal step"), true)
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
