//! Geometry-first grand-piano composition of existing FrankenSim owners.
//! Run `cargo run -p fs-couple --example grand_piano -- --dump-scale scale.csv`.
//! The default scale is ESTIMATED, not a measured Steinway digital twin.
mod geometry;

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let scale = match args.as_slice() {
        [flag, path] if flag == "--scale" => geometry::read_scale(&std::fs::read_to_string(path).map_err(|e| e.to_string())?)?,
        [flag, path] if flag == "--dump-scale" => {
            let scale = geometry::demonstration_scale()?;
            let text = format!("# ESTIMATED demonstration scale; not measured Steinway geometry\n{}", geometry::write_scale(&scale));
            std::fs::write(path, text).map_err(|e| e.to_string())?;
            scale
        }
        [] => geometry::demonstration_scale()?,
        _ => return Err("usage: grand_piano [--scale measurements.csv | --dump-scale estimates.csv]".into()),
    };
    println!("Model D published envelope: {} x {} m; board {} -> {} m (center -> edge).",
        geometry::D_LENGTH_M, geometry::D_WIDTH_M, geometry::D_BOARD_CENTER_M, geometry::D_BOARD_EDGE_M);
    println!("{} courses; {} speaking strings. Source authority belongs to the input, not the model name.", scale.len(), scale.iter().map(|c| c.unison).sum::<usize>());
    for c in scale { println!("key {}: L={:.4} m, mu={:.6} kg/m, T={:.2} N, f1={:.3} Hz", c.midi, c.length_m, c.linear_density_kg_m, c.tension_n, c.partial_hz(1, c.tension_n)); }
    Ok(())
}

fn main() {
    if let Err(error) = run() { eprintln!("grand_piano: {error}"); std::process::exit(1); }
}
