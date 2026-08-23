//! Emergent-gate battery for the bowed-string fixture
//! (bead frankensim-music-v8-root-3ez8g.7.5).
//!
//! Every gate observes per-sample histories produced ONLY by the admitted
//! gesture acting through the friction island; the Helmholtz corner is never
//! injected. The viscous-only FALSIFIER proves the stick-slip gate detects
//! the MECHANISM, not mere oscillation.
//!
//! Determinism: ONE-HOST bitwise replay (see module docs of
//! `fs_couple::bowed_string`). No wall-clock, no RNG.

use fs_couple::bowed_string::{
    BowGesture, BowGestureError, BowedRunConfig, BowedStringCard, FrictionIsland, Termination,
    cents_deviation, classify, gate_metrics, run_bowed, run_log_hash,
};
use fs_couple::stribeck_friction::StribeckFriction;
use fs_couple::thin_plate::CompactBody;
use fs_scenario::RadiatingPlate;

/// Shared string card: same physical family as `bakeoff_string.rs`, with 16
/// retained modes so the Helmholtz sawtooth has spectral room to sharpen.
fn card() -> BowedStringCard {
    BowedStringCard {
        length_m: 0.65,
        tension_n: 60.0,
        linear_density_kg_m: 6.0e-4,
        mode_count: 16,
        zetas: (0..16).map(|k| 1.0e-3 + 1.5e-4 * k as f64).collect(),
        sample_rate_hz: 48_000,
    }
}

/// Authored rosin-class coefficients (no measured corpus exists; registry
/// rows carry Estimate accordingly). Meaning under the event-stiction
/// coupling: `mu_static` caps the hold force, and `stiction_m_s` is the
/// KINETIC decay width of the Stribeck drop. It must span cm/s — wide
/// enough that the sliding operating point samples the negative-slope
/// region and pumps the oscillation (a 1 mm/s drop is flat everywhere
/// reachable and the string dead-slides).
fn rosin() -> StribeckFriction {
    StribeckFriction {
        mu_static: 0.8,
        mu_dynamic: 0.4,
        stiction_m_s: 0.04,
    }
}

/// Violin-like beta = 0.11 from the bridge; moderate bow force and speed.
fn playable_gesture() -> BowGesture {
    BowGesture::admit(0.20, 1.5, 0.11).expect("canonical gesture admits")
}

fn config(gesture: BowGesture, island: FrictionIsland, steps: usize) -> BowedRunConfig {
    BowedRunConfig {
        card: card(),
        island,
        gesture,
        steps,
        subsamples: 4,
        termination: Termination::Rigid,
        listener_m: 1.0,
    }
}

#[test]
fn gesture_admission_refuses_nonphysical_inputs() {
    assert_eq!(
        BowGesture::admit(0.2, -1.0, 0.11),
        Err(BowGestureError::NonPositiveNormalForce)
    );
    assert_eq!(
        BowGesture::admit(f64::NAN, 1.0, 0.11),
        Err(BowGestureError::NonFinite {
            what: "v_bow_m_s"
        })
    );
    assert_eq!(
        BowGesture::admit(0.2, f64::INFINITY, 0.11),
        Err(BowGestureError::NonFinite {
            what: "normal_force_n"
        })
    );
    assert_eq!(
        BowGesture::admit(0.2, 1.0, 0.0),
        Err(BowGestureError::StationOutOfRange)
    );
    assert_eq!(
        BowGesture::admit(0.2, 1.0, 1.0),
        Err(BowGestureError::StationOutOfRange)
    );
}

#[test]
fn bowed_run_replays_bitwise_on_one_host() {
    let cfg = || config(playable_gesture(), FrictionIsland::Stribeck(rosin()), 6_000);
    let a = run_bowed(&cfg()).expect("run A");
    let b = run_bowed(&cfg()).expect("run B");
    assert_eq!(
        run_log_hash(&a),
        run_log_hash(&b),
        "identical configs must replay bitwise"
    );
    let longer = run_bowed(&config(
        playable_gesture(),
        FrictionIsland::Stribeck(rosin()),
        6_001,
    ))
    .expect("run C");
    assert_ne!(
        run_log_hash(&a),
        run_log_hash(&longer),
        "a different step count must not collide"
    );
}

/// GATE (a)+(b): the emergent limit cycle sits at the STRING's transverse
/// fundamental with one flyback interval per period at the bow point.
#[test]
fn stick_slip_and_helmholtz_corner_emerge_at_the_string_fundamental() {
    let c = card();
    // 500 ms total; the metrics window keeps the last 200 ms of steady motion.
    let log = run_bowed(
        &config(playable_gesture(), FrictionIsland::Stribeck(rosin()), 24_000),
    )
    .expect("playable run stays inside state budgets");
    let m = gate_metrics(&log, &c);
    println!(
        "gate a/b metrics: f1={:.2} Hz peak={:.2} Hz ratio={:.2} slip={:.3} intervals/period={:.3} E_peak={:.3e} J",
        m.fundamental_hz,
        m.peak_hz,
        m.peak_to_semitone_ratio,
        m.slip_frac,
        m.intervals_per_period,
        log.peak_total_energy_j
    );
    assert_eq!(
        classify(&m),
        "playable",
        "stick-slip/Helmholtz gate failed on metrics {m:?}"
    );
}

/// FALSIFIER: flatten the curve to purely viscous opposition. There is no
/// stiction window, so the SAME classifier must refuse the run. This proves
/// the gate discriminates the mechanism, not oscillation in general.
#[test]
fn falsifier_viscous_only_friction_fails_the_stick_slip_gate() {
    let c = card();
    let log = run_bowed(
        &config(
            playable_gesture(),
            FrictionIsland::ViscousOnly {
                viscous_n_s_per_m: 8.0,
            },
            24_000,
        ),
    )
    .expect("viscous falsifier run completes");
    let m = gate_metrics(&log, &c);
    println!(
        "falsifier metrics: f1={:.2} Hz peak={:.2} Hz ratio={:.2} slip={:.3} intervals/period={:.3}",
        m.fundamental_hz,
        m.peak_hz,
        m.peak_to_semitone_ratio,
        m.slip_frac,
        m.intervals_per_period
    );
    assert_ne!(
        classify(&m),
        "playable",
        "viscous-only MUST fail the stick-slip gate; classifier is broken if it passes"
    );
}

/// GATE (d): within the playable band, higher bow force flattens the pitch
/// (the classic effect): measured cents drift DOWN vs F_n, logged as a
/// table. Thresholds are this rig's declared tolerances, not literature
/// constants transcribed as truth.
#[test]
fn pitch_flattens_as_bow_force_rises_inside_the_playable_band() {
    let c = card();
    let sr = f64::from(c.sample_rate_hz);
    let forces = [0.75, 1.0, 1.5, 2.25, 3.0];
    let mut cents_rows = Vec::new();
    for force in forces {
        let gesture = BowGesture::admit(0.20, force, 0.11).expect("sweep gesture admits");
        let log = run_bowed(
            &config(gesture, FrictionIsland::Stribeck(rosin()), 14_400),
        )
        .expect("flattening sweep run stays bounded");
        let win_start = log.bow_point_velocity_m_s.len() - 4_800;
        let cents = cents_deviation(
            &log.bow_point_velocity_m_s[win_start..],
            c.fundamental_hz() * 0.92,
            c.fundamental_hz() * 1.10,
            c.fundamental_hz(),
            sr,
        );
        println!("F_n={force:>5.2} N -> pitch {cents:+7.1} cents vs f1");
        cents_rows.push((force, cents));
    }
    let first = cents_rows[0].1;
    let last = cents_rows[cents_rows.len() - 1].1;
    println!("flattening span: {first:+.1} -> {last:+.1} cents");
    assert!(
        last < first,
        "higher bow force must not SHARPEN the pitch in this rig (got {first:+.1} -> {last:+.1})"
    );
    assert!(
        first - last < 200.0,
        "flattening beyond two semitones means the sweep left the playable band"
    );
}

/// Second configuration: the rigid bridge transmits its force into a
/// one-port plate body; body motion and radiation stay finite and bounded.
#[test]
fn plate_one_port_configuration_runs_bounded_and_logs_body_motion() {
    let body_spec = RadiatingPlate {
        area_m2: 3.0e-3,
        mass_kg: 0.15,
        frequency_hz: 280.0,
        damping_ratio: 0.02,
    };
    let body = CompactBody::from_radiator(body_spec).expect("body spec admits");
    let cfg = BowedRunConfig {
        card: card(),
        island: FrictionIsland::Stribeck(rosin()),
        gesture: playable_gesture(),
        steps: 12_000,
        subsamples: 4,
        termination: Termination::PlateOnePort(Box::new(body)),
        listener_m: 1.0,
    };
    let log = run_bowed(&cfg).expect("plate configuration runs");
    assert_eq!(log.body_velocity_m_s.len(), cfg.steps);
    assert_eq!(log.radiated_pressure_pa.len(), cfg.steps);
    assert!(
        log.body_velocity_m_s.iter().all(|v| v.is_finite()),
        "body velocity must stay finite"
    );
    assert!(
        log.radiated_pressure_pa.iter().all(|p| p.is_finite()),
        "radiation must stay finite"
    );
    let body_speed_max = log
        .body_velocity_m_s
        .iter()
        .fold(0.0_f64, |a, v| a.max(v.abs()));
    println!("plate body max |volume velocity| = {body_speed_max:.3e} m^3/s");
    assert!(body_speed_max > 0.0, "the body must actually move");
}

#[test]
#[ignore = "temporary diagnostic"]
fn debug_probe_stick_dynamics() {
    let c = card();
    for &(f_n, v_bow) in &[(3.9_f64, 0.45_f64), (2.0, 0.10), (8.0, 0.20)] {
        let g = BowGesture::admit(v_bow, f_n, 0.11).unwrap();
        let log = run_bowed(&config(g, FrictionIsland::Stribeck(rosin()), 24_000)).unwrap();
        let rel = &log.relative_velocity_m_s[12_000..];
        let bp = &log.bow_point_velocity_m_s[12_000..];
        let stuck_samples = rel.iter().filter(|v| **v == 0.0).count();
        let crossings = rel.windows(2).filter(|w| w[0].signum() != w[1].signum()).count();
        println!(
            "F={f_n} v={v_bow}: stuck%={} zc={} vrel[min={},max={},mean={:.4}] vstr[max]={:.4}",
            100.0 * stuck_samples as f64 / rel.len() as f64,
            crossings,
            rel.iter().cloned().fold(f64::INFINITY, f64::min),
            rel.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            rel.iter().sum::<f64>() / rel.len() as f64,
            bp.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
    }
    let _ = c;
}

#[test]
#[ignore = "temporary diagnostic"]
fn debug_probe_force_pipeline() {
    let c = card();
    println!("f1={}", c.fundamental_hz());
    let law = rosin();
    println!("traction(0.45, 3.9) = {}", law.traction(0.45, 3.9));
    let g = BowGesture::admit(0.45, 3.9, 0.11).unwrap();
    let mut cfg = config(g, FrictionIsland::Stribeck(law), 8);
    cfg.subsamples = 1;
    let log = run_bowed(&cfg).unwrap();
    println!("bow_point_v[0..8]: {:?}", &log.bow_point_velocity_m_s[..]);
    println!("rel[0..8]: {:?}", &log.relative_velocity_m_s[..]);
    println!("E_final={:.3e}", log.final_total_energy_j);

    // Direct model probe bypassing run_bowed entirely.
    use fs_couple::bowed_string::BowedRunError;
    let _ : Option<BowedRunError> = None;
    let modes: Vec<fs_couple::modal_acoustic_time::ModalAcousticMode> = (0..c.mode_count)
        .map(|k| fs_couple::modal_acoustic_time::ModalAcousticMode {
            angular_frequency_rad_s: (k + 1) as f64 * std::f64::consts::PI
                * c.wave_speed_m_s() / c.length_m,
            damping_ratio: c.zetas[k],
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        })
        .collect();
    let mut model = fs_couple::modal_acoustic_time::ModalAcousticTimeModel::try_new(
        c.sample_rate_hz,
        modes,
        fs_couple::modal_acoustic_time::ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let phi1 = (std::f64::consts::PI * 0.11).sin() / (c.linear_density_kg_m * c.length_m * 0.5).sqrt();
    println!("phi1={phi1}");
    let q1 = 77.0;
    let mut f = vec![0.0; c.mode_count];
    f[0] = q1;
    for i in 0..6 {
        model.step_duration(&f, 1.0 / 48_000.0).unwrap();
        println!("direct step {i}: q1v={:+.6e} q1x={:+.6e}",
            model.states()[0].velocity_m_sqrt_kg_per_s,
            model.states()[0].displacement_m_sqrt_kg);
    }
}
