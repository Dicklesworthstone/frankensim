//! Geometry-first grand-piano composition of existing FrankenSim owners.
//! `cargo run --release -p fs-couple --example grand_piano -- --preset steinway-d --render piano.wav`
//! Source-derived geometry is approximate, not a measured specimen digital twin.
mod geometry;
mod linear;
mod board;
mod board_geometry;
mod steinway_d;
mod performance;
mod felt;
mod engine;

const USAGE: &str = "grand_piano [--render piano.wav] [--scale strings.csv]
    [--preset steinway-d | --board board.csv | --board-geometry panel.fsb]
    [--mesh-divisions 4..24] [--dump-geometry panel.fsb] [--dump-obj soundboard.obj]
    [--board-band-hz Hz] [--performance events.csv] [--observer-gain Pa/(m^3/s)]
    [--note 21..108] [--velocity m/s] [--duration seconds]
    [--sample-rate Hz] [--substeps 1..16] [--modes 1..256]
    [--dump-scale strings.csv] [--dump-board board.csv]
--preset steinway-d reconstructs the published 17-rib Model D drawing, with
spruce panel, sugar-pine ribs, maple bridges, cut-off bar and 88 bridge stations.
--dump-geometry/--dump-obj export that same physical model; export alone skips
eigenanalysis. Thickness taper, material constants and key assignment include
explicit estimates. The string scale and hammer voicing remain replaceable estimates.
--note performs a single-key study; otherwise the demo also plays a chord of
available keys. --velocity overrides the three demo hammer launch speeds.
Velocity is POST-ESCAPEMENT hammer velocity, not MIDI velocity or key motion.
Geometric boards assemble a flat orthotropic plate before rendering. The explicit
frequency band admits at most 32 modes; mesh refinement is not a convergence claim.
--performance uses sample,event,key,value CSV instead of the demo and cannot be
combined with --note or --velocity. Note-on values are hammer velocity in m/s.
Observer gain is diagnostic, not a measured acoustic radiation transfer.";

#[derive(Debug)]
struct Options {
    render: Option<String>, scale: Option<String>, board: Option<String>,
    board_geometry: Option<String>, performance: Option<String>, preset: Option<String>,
    mesh_divisions: usize, dump_geometry: Option<String>, dump_obj: Option<String>,
    board_band_hz: f64, observer_gain: f64,
    dump_scale: Option<String>, dump_board: Option<String>,
    note: Option<u8>, velocity: Option<f64>, duration: f64,
    sample_rate: u32, substeps: usize, modes: usize, help: bool,
}
impl Default for Options {
    fn default() -> Self {
        Self { render: None, scale: None, board: None, board_geometry: None,
            performance: None, preset: None, mesh_divisions: 8, dump_geometry: None, dump_obj: None,
            board_band_hz: 400.0, observer_gain: 10_000.0, dump_scale: None,
            dump_board: None, note: None, velocity: None, duration: 6.0,
            sample_rate: 48_000, substeps: 4, modes: 24, help: false }
    }
}
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut options = Self::default();
        let mut seen = std::collections::BTreeSet::new();
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            if flag == "--help" || flag == "-h" { options.help = true; continue; }
            if !seen.insert(flag.as_str()) { return Err(format!("duplicate option {flag}")); }
            let value = args.next().ok_or_else(|| format!("missing value for {flag}"))?;
            let invalid = || format!("invalid value for {flag}: {value}");
            match flag.as_str() {
                "--render" => options.render = Some(value.clone()),
                "--scale" => options.scale = Some(value.clone()),
                "--board" => options.board = Some(value.clone()),
                "--board-geometry" => options.board_geometry = Some(value.clone()),
                "--preset" => options.preset = Some(value.clone()),
                "--mesh-divisions" => options.mesh_divisions = value.parse().map_err(|_| invalid())?,
                "--dump-geometry" => options.dump_geometry = Some(value.clone()),
                "--dump-obj" => options.dump_obj = Some(value.clone()),
                "--performance" => options.performance = Some(value.clone()),
                "--board-band-hz" => options.board_band_hz = value.parse().map_err(|_| invalid())?,
                "--observer-gain" => options.observer_gain = value.parse().map_err(|_| invalid())?,
                "--dump-scale" => options.dump_scale = Some(value.clone()),
                "--dump-board" => options.dump_board = Some(value.clone()),
                "--note" => options.note = Some(value.parse().map_err(|_| invalid())?),
                "--velocity" => options.velocity = Some(value.parse().map_err(|_| invalid())?),
                "--duration" => options.duration = value.parse().map_err(|_| invalid())?,
                "--sample-rate" => options.sample_rate = value.parse().map_err(|_| invalid())?,
                "--substeps" => options.substeps = value.parse().map_err(|_| invalid())?,
                "--modes" => options.modes = value.parse().map_err(|_| invalid())?,
                _ => return Err(format!("unknown option {flag}\n{USAGE}")),
            }
        }
        if options.note.is_some_and(|n| !(21..=108).contains(&n))
            || options.velocity.is_some_and(|v| !v.is_finite() || v <= 0.0 || v > 8.0)
            || !options.duration.is_finite() || !(0.001..=120.0).contains(&options.duration)
            || !(8_000..=192_000).contains(&options.sample_rate)
            || !(1..=16).contains(&options.substeps) || !(1..=256).contains(&options.modes) {
            return Err("render control outside its finite admitted range".into());
        }
        if usize::from(options.board.is_some()) + usize::from(options.board_geometry.is_some())
            + usize::from(options.preset.is_some()) > 1 {
            return Err("choose a preset, modal board or geometric board, not multiple board sources".into());
        }
        if options.preset.as_deref().is_some_and(|p| p != "steinway-d") {
            return Err("unknown piano preset; available: steinway-d".into());
        }
        if !(4..=24).contains(&options.mesh_divisions)
            || (options.preset.is_none() && (seen.contains("--mesh-divisions")
                || options.dump_geometry.is_some() || options.dump_obj.is_some())) {
            return Err("mesh/export controls require --preset steinway-d; divisions must be 4..24".into());
        }
        if seen.contains("--board-band-hz") && options.board_geometry.is_none() && options.preset.is_none() {
            return Err("--board-band-hz requires --board-geometry or --preset".into());
        }
        if !options.board_band_hz.is_finite() || options.board_band_hz <= 0.0
            || options.board_band_hz >= 0.45 * f64::from(options.sample_rate)
            || !options.observer_gain.is_finite() || options.observer_gain <= 0.0 {
            return Err("invalid board frequency band or diagnostic observer gain".into());
        }
        if options.performance.is_some()
            && (options.render.is_none() || options.note.is_some() || options.velocity.is_some()) {
            return Err("--performance requires --render and replaces --note/--velocity demo controls".into());
        }
        // Do not overwrite the very measurements that a render was asked to use.
        let inputs = [options.scale.as_ref(), options.board.as_ref(),
            options.board_geometry.as_ref(), options.performance.as_ref()];
        let outputs = [options.render.as_ref(), options.dump_scale.as_ref(), options.dump_board.as_ref(),
            options.dump_geometry.as_ref(), options.dump_obj.as_ref()];
        for (i, output) in outputs.iter().enumerate() {
            if let Some(path) = output {
                if path.is_empty() || inputs.iter().flatten().any(|input| input == path)
                    || outputs[..i].iter().flatten().any(|previous| previous == path) {
                    return Err("output paths must be distinct from inputs and each other".into());
                }
            }
        }
        Ok(options)
    }
}

fn load_scale(text: Option<&str>) -> Result<Vec<geometry::Course>, String> {
    match text { Some(text) => geometry::read_scale(text), None => geometry::demonstration_scale() }
}
fn load_board(text: Option<&str>, scale: &[geometry::Course]) -> Result<Vec<linear::BoardMode>, String> {
    match text {
        Some(text) => board::read(text, &scale.iter().map(|c| c.midi).collect::<Vec<_>>()),
        None => Ok(board::demonstration()),
    }
}
/// Export only the admitted keys. An absent measurement must not become an
/// apparently measured zero bridge coefficient when the table is re-imported.
fn write_board_for_scale(modes: &[linear::BoardMode], scale: &[geometry::Course]) -> String {
    let all = board::write(modes);
    let mut out = String::new();
    for row in all.lines() {
        let admitted = row == board::HEADER || row.split(',').nth(4)
            .and_then(|key| key.parse::<u8>().ok())
            .is_some_and(|key| scale.iter().any(|c| c.midi == key));
        if admitted { out.push_str(row); out.push('\n'); }
    }
    out
}
fn study_key(scale: &[geometry::Course], requested: Option<u8>) -> Result<u8, String> {
    if let Some(key) = requested {
        return scale.iter().any(|c| c.midi == key).then_some(key)
            .ok_or_else(|| format!("requested key {key} is absent from the input scale"));
    }
    scale.iter().min_by_key(|c| c.midi.abs_diff(69)).map(|c| c.midi)
        .ok_or_else(|| "empty input scale".into())
}

fn render(path: &str, scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options) -> Result<(), String> {
    study_key(&scale, options.note)?;
    let keys: Vec<u8> = scale.iter().map(|c| c.midi).collect();
    let rate = options.sample_rate;
    let count = (options.duration * f64::from(rate)).round() as u32;
    let mut score = match &options.performance {
        Some(path) => performance::Performance::read(
            &std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?,
            &keys, u64::from(count))?,
        None => performance::Performance::demonstration(&keys, rate, u64::from(count),
            options.note, options.velocity)?,
    };
    let mut piano = engine::Instrument::new(scale, modes, rate,
        options.substeps, options.modes, true)?;
    debug_assert_eq!(piano.sample_rate(), rate);
    let mut pressure = Vec::with_capacity(count as usize);
    let start = std::time::Instant::now();
    let mut peak: f64 = 0.0;
    for sample in 0..count {
        score.dispatch(u64::from(sample), &mut piano)?;
        // Explicit diagnostic observer gain [Pa / (m^3/s)]. This is NOT a
        // measured radiation transfer or a claim of calibrated acoustic SPL.
        let p = options.observer_gain * piano.step().map_err(|e| format!("sample {sample}: {e}"))?;
        if !p.is_finite() { return Err(format!("sample {sample}: diagnostic observer overflow")); }
        peak = peak.max(p.abs());
        pressure.push(p);
    }
    let elapsed = start.elapsed().as_secs_f64();
    let seconds = f64::from(count) / f64::from(rate);
    let (wav, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&pressure, rate, 2.0).map_err(|e| e.to_string())?;
    std::fs::write(path, wav).map_err(|e| format!("{path}: {e}"))?;
    println!("Diagnostic volume-velocity observer, gain {} Pa/(m^3/s); no peak normalization or calibrated SPL claim.", options.observer_gain);
    println!("{} string modes; {} board modes; {} above-band duplex segments omitted from dynamic retention (static attachment retained).",
        piano.bank.modes.len(), piano.bank.board_count, piano.bank.omitted_duplex_modes);
    println!("{seconds:.6} s audio rendered in {elapsed:.6} s; wall/audio ratio {:.4}; peak {peak:.6} Pa-equivalent; {clips} PCM clips.", elapsed / seconds);
    println!("Input {:.9} J; stored {:.9} J; component losses {:.9} J; closure {:.3e} J; worst substep defect {:.3e} J.",
        piano.accounting.input_work_j, piano.energy_j(), piano.accounting.dissipated_j(),
        piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j(), piano.accounting.max_balance_error_j);
    Ok(())
}

fn run() -> Result<(), String> {
    let options = Options::parse(&std::env::args().skip(1).collect::<Vec<_>>())?;
    if options.help { println!("{USAGE}"); return Ok(()); }
    let read = |path: &String| std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"));
    let scale_text = options.scale.as_ref().map(read).transpose()?;
    let board_text = options.board.as_ref().map(read).transpose()?;
    let scale = load_scale(scale_text.as_deref())?;
    let preset = options.preset.as_ref().map(|_| steinway_d::build(options.mesh_divisions)).transpose()?;
    if let Some(preset) = &preset {
        if let Some(path) = &options.dump_geometry { std::fs::write(path, &preset.geometry).map_err(|e| format!("{path}: {e}"))?; }
        if let Some(path) = &options.dump_obj { std::fs::write(path, &preset.obj).map_err(|e| format!("{path}: {e}"))?; }
        if options.render.is_none() && options.dump_board.is_none() && options.dump_scale.is_none()
            && (options.dump_geometry.is_some() || options.dump_obj.is_some()) {
            println!("Exported source-derived Model D: tapered panel, 17 ribs, maple bridges, cut-off bar and 88 bridge stations. No eigenanalysis was needed.");
            return Ok(());
        }
    }
    let geometry_text = match (&preset, &options.board_geometry) {
        (Some(p), _) => Some(p.geometry.clone()),
        (_, Some(path)) => Some(read(path)?),
        _ => None,
    };
    let (modes, board_source) = if let Some(text) = &geometry_text {
        let start = std::time::Instant::now();
        let geometry = board_geometry::BoardGeometry::read(text)?;
        let prepared = geometry.prepare(&scale.iter().map(|c| c.midi).collect::<Vec<_>>(),
            options.board_band_hz)?;
        println!("Flat plate {:.6} m^2, {:.6} kg (panel+ribs/bridges), {} free DOFs, {} modes in (0,{}] Hz; preparation {:.6} s.",
            prepared.area_m2, prepared.mass_kg, prepared.free_dofs, prepared.modes.len(),
            options.board_band_hz, start.elapsed().as_secs_f64());
        for (i, interval) in prepared.frequency_intervals_hz.iter().enumerate() {
            println!("board mode {i}: [{:.9}, {:.9}] Hz", interval.0, interval.1);
        }
        (prepared.modes, format!("GEOMETRY-DERIVED FLAT PLATE; {}; crown/rim compliance and acoustic radiation not modeled", prepared.provenance))
    } else {
        (load_board(board_text.as_deref(), &scale)?, options.board.as_deref()
            .unwrap_or("AUTHORED illustrative modes; not measured Steinway geometry").to_owned())
    };
    if modes.iter().any(|m| m.frequency_hz >= 0.45 * f64::from(options.sample_rate)) {
        return Err("soundboard mode at/above output retention ceiling; use an explicitly reduced board".into());
    }
    study_key(&scale, options.note)?;
    println!("String scale: {}.", options.scale.as_deref().unwrap_or("ESTIMATED demonstration"));
    println!("Soundboard: {board_source}.");
    println!("Source authority belongs to the inputs, not the model name; imported files are not independently certified measurements.");
    if let Some(path) = &options.dump_scale {
        let source = options.scale.as_deref().unwrap_or("ESTIMATED demonstration; not measured Steinway geometry");
        std::fs::write(path, format!("# Source: {source}\n{}", geometry::write_scale(&scale))).map_err(|e| e.to_string())?;
    }
    if let Some(path) = &options.dump_board {
        std::fs::write(path, format!("# Source: {board_source}\n{}", write_board_for_scale(&modes, &scale)))
            .map_err(|e| e.to_string())?;
    }
    if let Some(path) = &options.render { return render(path, scale, &modes, &options); }
    println!("Model D published envelope: {} x {} m; board {} -> {} m (center -> edge).",
        geometry::D_LENGTH_M, geometry::D_WIDTH_M, geometry::D_BOARD_CENTER_M, geometry::D_BOARD_EDGE_M);
    println!("{} courses; {} speaking strings; {} board modes.", scale.len(), scale.iter().map(|c| c.unison).sum::<usize>(), modes.len());
    for c in scale { println!("key {}: L={:.4} m, mu={:.6} kg/m, T={:.2} N, f1={:.3} Hz", c.midi, c.length_m, c.linear_density_kg_m, c.tension_n, c.partial_hz(1, c.tension_n)); }
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("grand_piano: {error}"); std::process::exit(1); }
}

#[cfg(test)]
mod render_tests {
    use super::*;
    fn options(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
    }
    #[test]
    fn rendering_and_both_physical_imports_are_composable() {
        let o = options(&["--scale", "strings.csv", "--render", "piano.wav", "--board", "board.csv",
            "--sample-rate", "44100", "--substeps", "2", "--modes", "40", "--note", "21"]).unwrap();
        assert_eq!(o.render.as_deref(), Some("piano.wav"));
        assert_eq!(o.scale.as_deref(), Some("strings.csv"));
        assert_eq!(o.board.as_deref(), Some("board.csv"));
        assert_eq!((o.sample_rate, o.substeps, o.modes, o.note), (44_100, 2, 40, Some(21)));
    }
    #[test]
    fn imported_subset_reaches_the_prepared_instrument_without_substitution() {
        let mut course = geometry::demonstration_scale().unwrap()[0];
        course.length_m *= 0.97;
        let scale = load_scale(Some(&geometry::write_scale(&[course]))).unwrap();
        assert_eq!(scale, vec![course]);
        assert_eq!(study_key(&scale, None).unwrap(), 21);
        assert!(study_key(&scale, Some(69)).is_err());
        let mut input = board::demonstration();
        input.truncate(1);
        input[0].frequency_hz = 123.0;
        let modes = load_board(Some(&board::write(&input)), &scale).unwrap();
        assert_eq!(modes[0].frequency_hz, 123.0);
        let mut piano = engine::Instrument::new(scale, &modes, 8_000, 4, 4, true).unwrap();
        assert_eq!(piano.bank.board_count, 1);
        piano.note_on(21, 1.0).unwrap();
        for _ in 0..16 { assert!(piano.step().unwrap().is_finite()); }
    }
    #[test]
    fn invalid_inputs_never_fall_back_to_estimates() {
        assert!(load_scale(Some("bad CSV")).is_err());
        let scale = geometry::demonstration_scale().unwrap();
        assert!(load_board(Some(board::HEADER), &scale).is_err());
        for args in [vec!["--velocity", "NaN"], vec!["--duration", "inf"], vec!["--modes", "0"],
            vec!["--scale", "input.csv", "--render", "input.csv"], vec!["--render"],
            vec!["--note", "69", "--note", "70"], vec!["--substeps", "17"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn geometric_board_and_performance_options_compose_without_replacing_existing_controls() {
        let o = options(&["--scale", "strings.csv", "--board-geometry", "panel.fsb",
            "--board-band-hz", "300", "--performance", "score.csv", "--render", "piano.wav",
            "--sample-rate", "44100", "--duration", "2", "--modes", "40"]).unwrap();
        assert_eq!(o.board_geometry.as_deref(), Some("panel.fsb"));
        assert_eq!(o.performance.as_deref(), Some("score.csv"));
        assert_eq!(o.board_band_hz, 300.0);
        assert_eq!((o.sample_rate, o.duration, o.modes), (44_100, 2.0, 40));
        for args in [vec!["--board", "a.csv", "--board-geometry", "b.fsb"],
            vec!["--board-band-hz", "300"], vec!["--performance", "score.csv"],
            vec!["--render", "a.wav", "--performance", "score.csv", "--note", "69"],
            vec!["--render", "a.wav", "--performance", "score.csv", "--velocity", "1"],
            vec!["--render", "score.csv", "--performance", "score.csv"],
            vec!["--observer-gain", "NaN"],
            vec!["--board-geometry", "a.fsb", "--board-band-hz", "21600"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn board_export_preserves_missing_measurements_and_reimports_the_admitted_subset() {
        let scale = vec![geometry::demonstration_scale().unwrap()[48]];
        let modes = board::demonstration();
        let text = write_board_for_scale(&modes, &scale);
        assert_eq!(text.lines().count(), modes.len() + 1);
        let roundtrip = board::read(&text, &[69]).unwrap();
        for (a, b) in roundtrip.iter().zip(&modes) {
            assert_eq!(a.frequency_hz, b.frequency_hz);
            assert_eq!(a.bridge[48], b.bridge[48]);
            assert_eq!(a.volume, b.volume);
        }
        assert!(board::read(&text, &[60, 69]).is_err());
    }
    #[test]
    fn source_preset_renders_and_exports_without_conflicting_board_inputs() {
        let o = options(&["--preset", "steinway-d", "--render", "d.wav", "--board-band-hz", "350",
            "--dump-geometry", "d.fsb", "--dump-obj", "d.obj", "--mesh-divisions", "12"]).unwrap();
        assert_eq!(o.preset.as_deref(), Some("steinway-d"));
        assert_eq!(o.mesh_divisions, 12);
        for args in [vec!["--preset", "unknown"],vec!["--preset", "steinway-d", "--board", "b.csv"],
            vec!["--preset", "steinway-d", "--board-geometry", "b.fsb"],vec!["--dump-obj", "d.obj"],
            vec!["--preset", "steinway-d", "--dump-obj", "d.obj", "--render", "d.obj"],
            vec!["--preset", "steinway-d", "--mesh-divisions", "3"]] {assert!(options(&args).is_err());}
    }
}
