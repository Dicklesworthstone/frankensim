//! Kernel equivalence: the live kernel against the offline reference ENGINE
//! (`grand_piano/engine.rs::Instrument`, compiled from the same source file)
//! on an isolated single key, same board, same partial cap, one mechanical
//! step per sample. Nothing else hangs on the bridge in either model, so the
//! two must agree to roundoff: this checks the equations and discretization,
//! not the live reductions.
//!
//! usage: validate [board_band_hz] [partials] [seconds] [reduction_spec]
//! With a reduction spec the board is the reference's bridge-informed
//! reduction of the (0, band] slice and both models use its full damping matrix.
use fs_live_piano::refmods::{engine, geometry, steinway_scale};
use fs_live_piano::{Config, Description, LivePiano};

fn main() {
    let args: Vec<String> = std::env::args().collect();
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
            let course = scale[usize::from(midi - 21)];
            let desc = match &reduction {
                Some(spec) => Description::with_reduction(vec![course], band, spec),
                None => Description::with_courses(vec![course], band),
            }.expect("board");
            let material = steinway_scale::hammer_material(&course).expect("felt card");
            let mut reference = engine::Instrument::new_with_string_damping(vec![course], &desc.board, rate, 1,
                partials, true, vec![material], Some(engine::ShankGeometry::published()), None, None, false)
                .expect("reference engine");
            if let Some(c) = &desc.damping {
                reference.configure_bare_board_damping(c).expect("reference board damping");
            }
            let mut cfg = Config::live(rate);
            cfg.board_band_hz = band;
            cfg.partial_band_hz = 0.45 * f64::from(rate);
            cfg.max_partials = partials;
            cfg.course_mode_budget = usize::MAX / 4;
            cfg.static_attachment = false;
            cfg.damper_points = 1;
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
