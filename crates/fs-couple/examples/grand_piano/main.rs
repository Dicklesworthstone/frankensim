//! Geometry-first grand-piano composition of existing FrankenSim owners.
//! `cargo run --release -p fs-couple --example grand_piano -- --render piano.wav`
//! Defaults are ESTIMATED, not a measured Steinway digital twin.
mod geometry;
mod linear;
mod board;
mod felt;
mod engine;

fn render(path:&str)->Result<(),String>{
    let scale=geometry::demonstration_scale()?;
    let mut piano=engine::Instrument::new(scale,&board::demonstration(),48_000,4,24,true)?;
    piano.set_sustain(1.0).map_err(|e|e.to_string())?;
    let rate=piano.sample_rate();let mut pressure=Vec::with_capacity(rate as usize*6);
    let start=std::time::Instant::now();let mut peak:f64=0.0;
    for sample in 0..rate*6 {
        if sample==0{piano.note_on(69,0.6).map_err(|e|e.to_string())?;}
        if sample==rate{piano.note_on(69,2.0).map_err(|e|e.to_string())?;}
        if sample==rate*2{piano.note_on(69,4.5).map_err(|e|e.to_string())?;}
        if [rate/2,rate+rate/2,rate*2+rate/2].contains(&sample){piano.note_off(69).map_err(|e|e.to_string())?;}
        if sample==rate*3{for key in [48,60,64,67]{piano.note_on(key,2.5).map_err(|e|e.to_string())?;}}
        if sample==rate*4{for key in [48,60,64,67]{piano.note_off(key).map_err(|e|e.to_string())?;}}
        if sample==rate*4+rate/2{piano.set_sustain(0.5).map_err(|e|e.to_string())?;}
        if sample==rate*5{piano.set_sustain(0.0).map_err(|e|e.to_string())?;}
        // Explicit diagnostic observer gain [Pa / (m^3/s)]. This is NOT a
        // measured radiation transfer or a claim of calibrated acoustic SPL.
        let p=10_000.0*piano.step().map_err(|e|format!("sample {sample}: {e}"))?;
        peak=peak.max(p.abs());pressure.push(p);
    }
    let elapsed=start.elapsed().as_secs_f64();
    let (wav,clips)=fs_couple::pcm_wav::encode_pcm16_wav(&pressure,rate,2.0).map_err(|e|e.to_string())?;
    std::fs::write(path,wav).map_err(|e|e.to_string())?;
    println!("ESTIMATED scale/board; diagnostic volume-velocity observer, gain 10000 Pa/(m^3/s); no peak normalization.");
    println!("{} string modes; {} board modes; {} above-band duplex segments omitted from dynamic retention (static attachment retained).",
        piano.bank.modes.len(),piano.bank.board_count,piano.bank.omitted_duplex_modes);
    println!("6 s audio rendered in {elapsed:.6} s; wall/audio ratio {:.4}; peak {peak:.6} Pa-equivalent; {clips} PCM clips.",elapsed/6.0);
    println!("Input {:.9} J; stored {:.9} J; component losses {:.9} J; closure {:.3e} J; worst substep defect {:.3e} J.",
        piano.accounting.input_work_j,piano.energy_j(),piano.accounting.dissipated_j(),
        piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j(),piano.accounting.max_balance_error_j);
    Ok(())
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let scale = match args.as_slice() {
        [flag,path] if flag=="--render"=>return render(path),
        [flag,path] if flag=="--dump-board"=>{
            return std::fs::write(path,format!("# AUTHORED illustrative board, not measured Steinway geometry\n{}",board::write(&board::demonstration()))).map_err(|e|e.to_string());
        }
        [flag, path] if flag == "--scale" => geometry::read_scale(&std::fs::read_to_string(path).map_err(|e| e.to_string())?)?,
        [flag, path] if flag == "--dump-scale" => {
            let scale = geometry::demonstration_scale()?;
            let text = format!("# ESTIMATED demonstration scale; not measured Steinway geometry\n{}", geometry::write_scale(&scale));
            std::fs::write(path, text).map_err(|e| e.to_string())?;
            scale
        }
        [] => geometry::demonstration_scale()?,
        _ => return Err("usage: grand_piano [--render piano.wav | --scale measurements.csv | --dump-scale estimates.csv | --dump-board board.csv]".into()),
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
