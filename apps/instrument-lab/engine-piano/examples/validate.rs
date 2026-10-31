//! Kernel equivalence: the live kernel against the offline reference ENGINE
//! (`grand_piano/engine.rs::Instrument`, compiled from the same source file)
//! on an isolated single key, same board, same partial cap, one mechanical
//! step per sample. Nothing else hangs on the bridge in either model, so the
//! two must agree to roundoff: this checks the equations and discretization,
//! not the live reductions.
//!
//! usage: validate [board_band_hz] [partials] [seconds] [reduction_spec]
//!        [--rt0425-string-damping] [--rt0425-hammer-stiffness] [--rt0425-hammer-dissipation]
//!        [--felt-pack-hammers] [--string-waves]
//! With a reduction spec the board is the reference's bridge-informed
//! reduction of the (0, band] slice and both models use its full damping matrix.
//! The flags select the published laws in both models. `--felt-pack-hammers` plays the felt
//! pack's hammers (mass, felt card and rate) in both; the browser plays it with the string damping.
//! `--string-waves` stretches plain steel keys with longitudinal waves in both
//! (`stretch::plain_steel_waves_spec`, 10 kHz axial band): the reference iterates its
//! implicit discrete gradient, the live kernel kicks once per step, so the two then
//! agree to the time discretization, not to roundoff.
use fs_live_piano::refmods::{engine, geometry, hammer_felt, steinway_scale};
use fs_live_piano::{Config, Description, LivePiano};

fn main() {
    let all: Vec<String> = std::env::args().collect();
    const LAWS: [&str; 5] = ["--rt0425-string-damping", "--rt0425-hammer-stiffness", "--rt0425-hammer-dissipation", "--string-waves", "--felt-pack-hammers"];
    // A mistyped flag would otherwise validate the default law and report it as the published one.
    if let Some(bad) = all.iter().skip(1).find(|a| a.starts_with("--") && !LAWS.contains(&a.as_str())) {
        eprintln!("unknown option {bad}");
        std::process::exit(2);
    }
    let [string_damping, hammer_stiffness, hammer_dissipation, string_waves, felt_pack] = LAWS.map(|f| all.iter().any(|a| a == f));
    if felt_pack && (hammer_stiffness || hammer_dissipation) {
        eprintln!("--felt-pack-hammers excludes the RT-0425 hammer flags");
        std::process::exit(2);
    }
    let felt_card = felt_pack.then(|| hammer_felt::stulov_abel().expect("felt pack"));
    let args: Vec<String> = all.into_iter().filter(|a| !a.starts_with("--")).collect();
    let band: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(800.0);
    let partials: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(64);
    let seconds: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let reduction: Option<String> = args.get(4).cloned();
    let rate = 48_000u32;
    let scale = geometry::read_scale(fs_live_piano::kernel::SCALE_CSV).expect("scale");
    // The embedded dump must be the reference preset scale retuned to A4=440.
    let mut preset = steinway_scale::courses().expect("preset scale");
    for c in &mut preset {
        let target = 440.0 * fs_math::det::pow(2.0, (f64::from(c.midi) - 69.0) / 12.0);
        let cents = 1200.0 * fs_math::det::ln(target / c.partial_hz(1, c.tension_n)) / std::f64::consts::LN_2;
        c.tension_n = c.tension_at_cents(cents).expect("retune");
    }
    println!("embedded scale == steinway_scale::courses() retuned to A4=440: {}", preset == scale);
    println!("key,velocity_m_s,strings,partials,samples,max_bridge_v_ref,max_abs_dv,rel_dv,input_ref_j,input_live_j,felt_ref_j,felt_live_j,modal_ref_j,modal_live_j,stored_ref_j,stored_live_j,closure_ref_j,residual_live_j,contact_samples_live,peak_force_n");
    let steps = (seconds * f64::from(rate)) as usize;
    for &midi in &[24u8, 33, 45, 60, 72, 84, 96, 108] {
        for &velocity in &[0.5, 2.5, 4.5] {
            let mut course = scale[usize::from(midi - 21)];
            let pack_hammer = felt_card.as_ref().map(|card| card.hammer(midi).expect("pack hammer"));
            if let Some(h) = &pack_hammer {
                course.hammer_mass_kg = h.mass_kg;
            }
            let desc = match &reduction {
                Some(spec) => Description::with_reduction(vec![course], band, spec),
                None => Description::with_courses(vec![course], band),
            }.expect("board");
            let material = if let Some(h) = &pack_hammer { h.material(&course) }
                else if hammer_dissipation { steinway_scale::hammer_material_rt0425_damped(&course) }
                else if hammer_stiffness { steinway_scale::hammer_material_rt0425(&course) }
                else { steinway_scale::hammer_material(&course) }.expect("felt card");
            let mut reference = engine::Instrument::new_with_string_damping(vec![course], &desc.board, rate, 1,
                partials, true, vec![material], Some(engine::ShankGeometry::published()), None, None, string_damping)
                .expect("reference engine");
            if let Some(h) = &pack_hammer {
                reference.configure_source_hammer_dissipation(&[h.rate_per_string(&course)]).expect("reference R");
            } else if hammer_dissipation {
                let r_h = steinway_scale::hammer_relaxation_rt0425(midi).expect("R_H");
                reference.configure_source_hammer_dissipation(&[r_h]).expect("reference R_H");
            }
            if let Some(c) = &desc.damping {
                reference.configure_bare_board_damping(c).expect("reference board damping");
            }
            // Plain steel strings stretch with longitudinal waves in both models (wound ones stay linear).
            let waves = string_waves.then(|| fs_live_piano::stretch::plain_steel_waves_spec(&[course], 10_000.0, 0.01, 0.05));
            if let Some(text) = &waves {
                let mut spec = fs_live_piano::refmods::linear::string_stretching::Specification::read(text, &[course]).expect("stretching file");
                // The reference takes a per-mode axial row; to compare the same equations it gets the live
                // engine's own labelled scalar times its vertical bridge row (`stretch::LIVE_BRIDGE_AXIAL_RATIO`).
                let rows = spec.waves.keys().map(|&key| {
                    (key, desc.board.iter().map(|m| fs_live_piano::stretch::LIVE_BRIDGE_AXIAL_RATIO * m.bridge[usize::from(key - 21)]).collect())
                }).collect();
                spec.set_bridge_axial_rows(rows).expect("axial rows");
                reference.configure_string_stretching(&spec).expect("reference stretching");
            }
            let mut cfg = Config::live(rate);
            cfg.string_stretching = waves;
            cfg.board_band_hz = band;
            cfg.partial_band_hz = 0.45 * f64::from(rate);
            cfg.max_partials = partials;
            cfg.course_mode_budget = usize::MAX / 4;
            cfg.static_attachment = false;
            cfg.damper_points = 1;
            cfg.rt0425_string_damping = string_damping;
            cfg.rt0425_hammer_stiffness = hammer_stiffness;
            cfg.rt0425_hammer_dissipation = hammer_dissipation;
            cfg.felt_pack_hammers = felt_pack;
            let mut live = LivePiano::new(&desc, cfg).expect("live kernel");
            let rest_energy = reference.energy_j();
            reference.note_on(midi, velocity).expect("reference strike");
            live.note_on(midi, velocity).expect("live strike");
            live.apply_pending_notes();
            let row: Vec<f64> = reference.bank.strings[0].bridge.clone();
            let n = reference.bank.modes.len();
            let (mut peak, mut worst) = (0.0_f64, 0.0_f64);
            let (mut contact, mut peak_force) = (0usize, 0.0_f64);
            for i in 0..steps {
                reference.step().expect("reference step");
                live.step();
                let vr: f64 = row.iter().zip(&reference.bank.v[n..]).map(|(g, v)| g * v).sum();
                let vl = live.bridge_velocity(midi).expect("bridge");
                peak = peak.max(vr.abs());
                worst = worst.max((vr - vl).abs());
                let _ = i;
                let force = live.hammer_force(midi).unwrap_or(0.0);
                if force > 0.0 { contact += 1; }
                peak_force = peak_force.max(force);
            }
            let a = reference.accounting;
            let l = live.close_ledger();
            println!("{midi},{velocity},{},{},{steps},{peak:.6e},{worst:.3e},{:.3e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.3e},{:.3e},{contact},{peak_force:.3}",
                live.dynamic_strings(), live.dynamic_modes(), worst / peak,
                a.input_work_j, l.input_j, a.felt_loss_j, l.felt_loss_j, a.modal_loss_j, l.string_loss_j + l.board_loss_j,
                reference.energy_j() - rest_energy, l.stored_j,
                a.input_work_j - (reference.energy_j() - rest_energy) - a.dissipated_j(), l.residual_j);
            if let Some(e) = live.take_error() { println!("  live error: {e}"); }
        }
    }
}
