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
    BowGesture, BowGestureError, BowedRunConfig, BowedRunError, BowedStringCard, FrictionIsland,
    Termination, cents_deviation, classify, gate_metrics, run_bowed, run_log_hash,
};
use fs_couple::stribeck_friction::StribeckFriction;
use fs_couple::thin_plate::CompactBody;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::RadiatingPlate;

/// Shared string card: same physical family as `bakeoff_string.rs`, with 16
/// retained modes so the Helmholtz sawtooth has spectral room to sharpen.
fn card() -> BowedStringCard {
    BowedStringCard {
        length_m: 0.65,
        tension_n: 60.0,
        linear_density_kg_m: 6.0e-4,
        bending_stiffness_n_m2: 0.0,
        viscous_bending_n_m2_s: 0.0,
        mode_count: 16,
        zetas: (0..16)
            .map(|k| 1.0e-3 * (1.0 + 0.55 * f64::from(k)))
            .collect(),
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

fn interface_island(law: StribeckFriction, speed_min: f64, speed_max: f64) -> FrictionIsland {
    use fs_tribo::{
        ApplicabilityRange as Range, DryFrictionApplicability, DryInterfaceSystemCard, FrictionLaw,
        InputAuthority, InterfaceMedium, InterfaceSystemRef,
    };
    FrictionIsland::InterfaceStribeck {
        card: Box::new(
            DryInterfaceSystemCard::new(
                InterfaceSystemRef::new(
                    "bow/string",
                    "unworn",
                    "synthetic-test",
                    InputAuthority::SyntheticFixture,
                    InterfaceMedium::Dry,
                )
                .unwrap(),
                FrictionLaw::Stribeck {
                    static_mu: law.mu_static,
                    kinetic_mu: law.mu_dynamic,
                    characteristic_speed: law.stiction_m_s,
                    viscous_per_speed: 0.0,
                },
                DryFrictionApplicability::new(
                    Range::new(280.0, 320.0).unwrap(),
                    Range::new(5.0e5, 2.0e6).unwrap(),
                    Range::new(speed_min, speed_max).unwrap(),
                )
                .unwrap(),
            )
            .unwrap(),
        ),
        temperature_kelvin: 300.0,
        nominal_contact_area_m2: 1.0e-6,
    }
}

/// G1/G3: card coefficients drive the actual oscillator, not just admission.
#[test]
fn ordered_interface_drives_bowed_motion_from_its_declared_law() {
    let mut velocities = Vec::new();
    for kinetic in [0.2, 0.4] {
        let law = StribeckFriction::try_new(0.8, kinetic, 0.04).unwrap();
        let mut cfg = config(
            BowGesture::admit(0.2, 1.0, 0.2).unwrap(),
            interface_island(law, 0.0, 10.0),
            1,
        );
        cfg.card.mode_count = 1;
        cfg.card.zetas = vec![0.0];
        cfg.subsamples = 1;
        let omega = cfg.card.mode_omega_rad_s(1);
        let phi = (core::f64::consts::PI * 0.2).sin()
            / (cfg.card.linear_density_kg_m * cfg.card.length_m / 2.0).sqrt();
        let force = kinetic + (0.8 - kinetic) * (-25.0_f64).exp();
        let expected = phi * phi * force * (omega / 48_000.0).sin() / omega;
        let run = run_bowed(&cfg).unwrap();
        assert!((run.bow_point_velocity_m_s[0] - expected).abs() < 1.0e-12);
        velocities.push(run.bow_point_velocity_m_s[0]);
        let mut authored = cfg.clone();
        authored.island = FrictionIsland::Stribeck(law);
        assert_eq!(
            run_log_hash(&run),
            run_log_hash(&run_bowed(&authored).unwrap())
        );
    }
    assert!(velocities[1] > 1.9 * velocities[0]);
}

/// G0/G3: state applicability is checked at admission and after actual motion,
/// including the sticking branch and the final substep of a one-sample run.
#[test]
fn ordered_interface_checks_contact_state_throughout_bowed_run() {
    let mut cfg = config(
        BowGesture::admit(0.2, 1.0, 0.2).unwrap(),
        interface_island(rosin(), 0.0, 10.0),
        1,
    );
    cfg.card.mode_count = 1;
    cfg.card.zetas = vec![0.0];
    cfg.subsamples = 1;
    assert!(run_bowed(&cfg).is_ok());
    // Initial 0.2 m/s is admitted; acceleration leaves this narrow domain.
    cfg.island = interface_island(rosin(), 0.199999, 0.2);
    assert!(matches!(
        run_bowed(&cfg),
        Err(BowedRunError::Interface(
            fs_tribo::TriboError::OutsideApplicability {
                field: "slip_speed_mps",
                ..
            }
        ))
    ));
    cfg.island = interface_island(rosin(), 0.0, 10.0);
    cfg.gesture = BowGesture::admit(0.0, 1.0, 0.2).unwrap();
    assert!(run_bowed(&cfg).is_ok());
    if let FrictionIsland::InterfaceStribeck {
        temperature_kelvin, ..
    } = &mut cfg.island
    {
        *temperature_kelvin = 330.0;
    }
    assert!(matches!(
        run_bowed(&cfg),
        Err(BowedRunError::Interface(
            fs_tribo::TriboError::OutsideApplicability {
                field: "temperature_kelvin",
                ..
            }
        ))
    ));
    cfg.island = interface_island(rosin(), 0.0, 10.0);
    cfg.gesture.normal_force_n = 3.0;
    assert!(matches!(
        run_bowed(&cfg),
        Err(BowedRunError::Interface(
            fs_tribo::TriboError::OutsideApplicability {
                field: "contact_pressure_pa",
                ..
            }
        ))
    ));
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        cfg.island = interface_island(rosin(), 0.0, 10.0);
        if let FrictionIsland::InterfaceStribeck {
            nominal_contact_area_m2,
            ..
        } = &mut cfg.island
        {
            *nominal_contact_area_m2 = area;
        }
        assert!(matches!(run_bowed(&cfg), Err(BowedRunError::Friction(_))));
    }
    // Another law family cannot silently acquire a Stribeck capture model.
    cfg.island = interface_island(rosin(), 0.0, 10.0);
    if let FrictionIsland::InterfaceStribeck { card, .. } = &mut cfg.island {
        *card = Box::new(
            fs_tribo::DryInterfaceSystemCard::new(
                card.interface().clone(),
                fs_tribo::FrictionLaw::Coulomb {
                    static_mu: 0.8,
                    kinetic_mu: 0.4,
                },
                card.applicability(),
            )
            .unwrap(),
        );
    }
    assert!(matches!(run_bowed(&cfg), Err(BowedRunError::Friction(_))));
}

/// G1/G3: the complete card coefficient mu(v) drives signed modal motion.
/// The viscous slope is [s/m], so its force scales with both speed and load.
#[test]
fn ordered_interface_velocity_strengthening_reaches_bowed_motion() {
    for speed in [-0.4_f64, -0.2, 0.2, 0.4] {
        for normal in [0.5, 1.0, 2.0] {
            let mut responses = Vec::new();
            for slope in [0.0, 1.5] {
                let mut island = interface_island(rosin(), 0.0, 10.0);
                if let FrictionIsland::InterfaceStribeck { card, .. } = &mut island {
                    *card = Box::new(
                        fs_tribo::DryInterfaceSystemCard::new(
                            card.interface().clone(),
                            fs_tribo::FrictionLaw::Stribeck {
                                static_mu: 0.8,
                                kinetic_mu: 0.4,
                                characteristic_speed: 0.04,
                                viscous_per_speed: slope,
                            },
                            card.applicability(),
                        )
                        .unwrap(),
                    );
                }
                let mut cfg = config(BowGesture::admit(speed, normal, 0.2).unwrap(), island, 1);
                cfg.card.mode_count = 1;
                cfg.card.zetas = vec![0.0];
                cfg.subsamples = 1;
                let omega = cfg.card.mode_omega_rad_s(1);
                let phi = (core::f64::consts::PI * 0.2).sin()
                    / (cfg.card.linear_density_kg_m * cfg.card.length_m / 2.0).sqrt();
                let dt = 1.0 / f64::from(cfg.card.sample_rate_hz);
                let force = normal
                    * (0.4 + 0.4 * (-(speed / 0.04).powi(2)).exp() + slope * speed.abs())
                    * speed.signum();
                let expected = phi * phi * force * (omega * dt).sin() / omega;
                let run = run_bowed(&cfg).unwrap();
                assert!((run.bow_point_velocity_m_s[0] - expected).abs() < 1.0e-12);
                responses.push(run.bow_point_velocity_m_s[0]);
            }
            assert!(responses[1].abs() > 1.5 * responses[0].abs());
        }
    }
}

/// G1: independent pinned-beam dispersion, held-load response, energy and
/// support traction. A frequency-only change misses the bending shear here.
#[test]
fn bending_stiffness_reaches_bowed_motion_energy_and_bridge_force() {
    for ei in [0.0, 0.002, 0.02] {
        let mut cfg = config(
            BowGesture::admit(0.2, 1.0, 0.2).unwrap(),
            FrictionIsland::ViscousOnly {
                viscous_n_s_per_m: 0.01,
            },
            1,
        );
        cfg.card.bending_stiffness_n_m2 = ei;
        cfg.card.mode_count = 3;
        cfg.card.zetas = vec![0.0; 3];
        cfg.subsamples = 1;
        let dt = 1.0 / f64::from(cfg.card.sample_rate_hz);
        let norm = (cfg.card.linear_density_kg_m * cfg.card.length_m / 2.0).sqrt();
        let mut velocity = 0.0;
        let mut energy = 0.0;
        let mut bridge = 0.0;
        let mut tension_only_bridge = 0.0;
        for mode in 1..=3 {
            let kappa = mode as f64 * core::f64::consts::PI / cfg.card.length_m;
            let omega = ((cfg.card.tension_n * kappa.powi(2) + ei * kappa.powi(4))
                / cfg.card.linear_density_kg_m)
                .sqrt();
            assert!((cfg.card.mode_omega_rad_s(mode) / omega - 1.0).abs() < 1.0e-13);
            if mode == 1 {
                assert!(
                    (cfg.card.fundamental_hz() * core::f64::consts::TAU / omega - 1.0).abs()
                        < 1.0e-13
                );
            }
            let phi = (mode as f64 * core::f64::consts::PI * 0.2).sin() / norm;
            let force = phi * 0.002; // Initial rest, 0.01 N s/m times 0.2 m/s.
            let q = force * (1.0 - (omega * dt).cos()) / omega.powi(2);
            let v = force * (omega * dt).sin() / omega;
            velocity += phi * v;
            energy += 0.5 * (v * v + omega * omega * q * q);
            tension_only_bridge += cfg.card.tension_n * kappa / norm * q;
            bridge += (cfg.card.tension_n * kappa + ei * kappa.powi(3)) / norm * q;
        }
        let log = run_bowed(&cfg).unwrap();
        assert!((log.bow_point_velocity_m_s[0] - velocity).abs() < 1.0e-12);
        assert!((log.bridge_force_n[0] - bridge).abs() < 1.0e-12);
        assert!((log.final_total_energy_j - energy).abs() < 1.0e-18);
        assert!((log.peak_total_energy_j - energy).abs() < 1.0e-18);
        if ei > 0.0 {
            // Omitting bending shear must exceed the comparison tolerance
            // by a wide margin, including the lower-stiffness case.
            assert!((bridge - tension_only_bridge).abs() > 100.0 * 1.0e-12);
            assert!(cfg.card.mode_omega_rad_s(3) > 3.0 * cfg.card.mode_omega_rad_s(1));
        }
    }
}

#[test]
fn bowed_bending_stiffness_refuses_negative_or_nonfinite_values() {
    let mut cfg = config(playable_gesture(), FrictionIsland::Stribeck(rosin()), 0);
    for ei in [-0.01, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        cfg.card.bending_stiffness_n_m2 = ei;
        assert!(matches!(
            run_bowed(&cfg),
            Err(BowedRunError::InvalidCard { .. })
        ));
    }
}

/// G1: a viscous beam transmits velocity-dependent shear at the support.
#[test]
fn viscous_bending_reaches_bridge_shear_and_modal_decay() {
    let mut cfg = config(
        BowGesture::admit(0.2, 1.0, 0.2).unwrap(),
        FrictionIsland::ViscousOnly {
            viscous_n_s_per_m: 0.01,
        },
        1,
    );
    cfg.subsamples = 1;
    cfg.card.mode_count = 1;
    cfg.card.bending_stiffness_n_m2 = 0.02;
    cfg.card.viscous_bending_n_m2_s = 0.003;
    let card = &mut cfg.card;
    let kappa = core::f64::consts::PI / card.length_m;
    let omega = ((card.tension_n * kappa.powi(2) + card.bending_stiffness_n_m2 * kappa.powi(4))
        / card.linear_density_kg_m)
        .sqrt();
    let decay = 0.5 * card.viscous_bending_n_m2_s * kappa.powi(4) / card.linear_density_kg_m;
    card.zetas = vec![decay / card.mode_omega_rad_s(1)];
    let dt = 1.0 / f64::from(card.sample_rate_hz);
    let norm = (card.linear_density_kg_m * card.length_m / 2.0).sqrt();
    let phi = (core::f64::consts::PI * 0.2).sin() / norm;
    let force = 0.002 * phi;
    let wd = (omega * omega - decay * decay).sqrt();
    let q = force / omega.powi(2)
        * (1.0 - (-decay * dt).exp() * ((wd * dt).cos() + decay / wd * (wd * dt).sin()));
    let v = force * (-decay * dt).exp() * (wd * dt).sin() / wd;
    let elastic = (card.tension_n * kappa + card.bending_stiffness_n_m2 * kappa.powi(3)) / norm * q;
    let shear = card.viscous_bending_n_m2_s * kappa.powi(3) / norm * v;
    let log = run_bowed(&cfg).unwrap();
    assert!((log.bow_point_velocity_m_s[0] - phi * v).abs() < 1e-12);
    assert!((log.bridge_force_n[0] - elastic - shear).abs() < 1e-12);
    assert!(shear.abs() > 1e-8);
    cfg.card.zetas[0] = 0.0;
    assert!(matches!(
        run_bowed(&cfg),
        Err(BowedRunError::InvalidCard { .. })
    ));
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        cfg.card.viscous_bending_n_m2_s = invalid;
        assert!(matches!(
            run_bowed(&cfg),
            Err(BowedRunError::InvalidCard { .. })
        ));
    }
}

/// Verified playable point of THIS rig (see the flattening-sweep table):
/// violin-like beta = 0.11 from the bridge, forte bow speed, firm force.
fn playable_gesture() -> BowGesture {
    BowGesture::admit(0.45, 3.9, 0.11).expect("canonical gesture admits")
}

fn config(gesture: BowGesture, island: FrictionIsland, steps: usize) -> BowedRunConfig {
    BowedRunConfig {
        card: card(),
        island,
        gesture,
        steps,
        subsamples: 16,
        termination: Termination::Rigid,
        listener_m: 1.0,
    }
}

#[test]
fn empty_bowed_run_reports_rest_energy_without_samples() {
    let log = run_bowed(&config(
        playable_gesture(),
        FrictionIsland::Stribeck(rosin()),
        0,
    ))
    .expect("an admitted empty run stays at rest");
    assert!(log.bow_point_velocity_m_s.is_empty());
    assert!(log.relative_velocity_m_s.is_empty());
    assert!(log.bridge_force_n.is_empty());
    assert!(log.body_volume_velocity_m3_s.is_empty());
    assert!(log.radiated_pressure_pa.is_empty());
    assert_eq!(log.final_total_energy_j.to_bits(), 0.0_f64.to_bits());
    assert_eq!(log.peak_total_energy_j.to_bits(), 0.0_f64.to_bits());
}

#[test]
fn gesture_admission_refuses_nonphysical_inputs() {
    assert_eq!(
        BowGesture::admit(0.2, -1.0, 0.11),
        Err(BowGestureError::NonPositiveNormalForce)
    );
    assert_eq!(
        BowGesture::admit(f64::NAN, 1.0, 0.11),
        Err(BowGestureError::NonFinite { what: "v_bow_m_s" })
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
    let log = run_bowed(&config(
        playable_gesture(),
        FrictionIsland::Stribeck(rosin()),
        24_000,
    ))
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
    let log = run_bowed(&config(
        playable_gesture(),
        FrictionIsland::ViscousOnly {
            viscous_n_s_per_m: 8.0,
        },
        24_000,
    ))
    .expect("viscous falsifier run completes");
    let m = gate_metrics(&log, &c);
    println!(
        "falsifier metrics: f1={:.2} Hz peak={:.2} Hz ratio={:.2} slip={:.3} intervals/period={:.3}",
        m.fundamental_hz, m.peak_hz, m.peak_to_semitone_ratio, m.slip_frac, m.intervals_per_period
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
    let forces = [2.0, 2.75, 3.25, 3.9];
    let mut cents_rows = Vec::new();
    for force in forces {
        let gesture = BowGesture::admit(0.45, force, 0.11).expect("sweep gesture admits");
        let log = run_bowed(&config(gesture, FrictionIsland::Stribeck(rosin()), 14_400))
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
        subsamples: 16,
        termination: Termination::PlateOnePort {
            body: Box::new(body),
            ambient: GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0)
                .expect("declared ambient admits"),
        },
        listener_m: 1.0,
    };
    let log = run_bowed(&cfg).expect("plate configuration runs");
    assert_eq!(log.body_volume_velocity_m3_s.len(), cfg.steps);
    assert_eq!(log.radiated_pressure_pa.len(), cfg.steps);
    assert!(
        log.body_volume_velocity_m3_s.iter().all(|v| v.is_finite()),
        "body volume velocity must stay finite"
    );
    assert!(
        log.radiated_pressure_pa.iter().all(|p| p.is_finite()),
        "radiation must stay finite"
    );
    let body_speed_max = log
        .body_volume_velocity_m3_s
        .iter()
        .fold(0.0_f64, |a, v| a.max(v.abs()));
    println!("plate body max |volume velocity| = {body_speed_max:.3e} m^3/s");
    assert!(body_speed_max > 0.0, "the body must actually move");
}

/// A public `GasState` can be forged, so the plate boundary must refuse an
/// invalid density before it can publish a radiation trace.
#[test]
fn plate_one_port_refuses_nonphysical_declared_ambient_density() {
    let body_spec = RadiatingPlate {
        area_m2: 3.0e-3,
        mass_kg: 0.15,
        frequency_hz: 280.0,
        damping_ratio: 0.02,
    };
    for density_kg_m3 in [f64::NAN, 0.0, -1.0] {
        let mut ambient =
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
        ambient.density = density_kg_m3;
        let cfg = BowedRunConfig {
            card: card(),
            island: FrictionIsland::Stribeck(rosin()),
            gesture: playable_gesture(),
            steps: 1,
            subsamples: 16,
            termination: Termination::PlateOnePort {
                body: Box::new(CompactBody::from_radiator(body_spec).expect("body spec admits")),
                ambient,
            },
            listener_m: 1.0,
        };
        match run_bowed(&cfg) {
            Err(BowedRunError::InvalidAmbientDensity {
                density_kg_m3: refused,
            }) if refused.to_bits() == density_kg_m3.to_bits() => {}
            result => {
                panic!("invalid ambient density {density_kg_m3:?} must refuse, got {result:?}")
            }
        }
    }
}

/// G3: observer distance changes inverse-distance pressure, not the trajectory.
#[test]
fn plate_one_port_listener_distance_scales_pressure_without_changing_motion() {
    let mut near = config(playable_gesture(), FrictionIsland::Stribeck(rosin()), 1024);
    near.termination = Termination::PlateOnePort {
        body: Box::new(
            CompactBody::from_radiator(RadiatingPlate {
                area_m2: 3.0e-3,
                mass_kg: 0.15,
                frequency_hz: 280.0,
                damping_ratio: 0.02,
            })
            .unwrap(),
        ),
        ambient: GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap(),
    };
    near.listener_m = 1.0;
    let mut far = near.clone();
    far.listener_m = 2.0;
    let near_log = run_bowed(&near).unwrap();
    let far_log = run_bowed(&far).unwrap();
    assert_eq!(
        near_log.bow_point_velocity_m_s,
        far_log.bow_point_velocity_m_s
    );
    assert_eq!(near_log.bridge_force_n, far_log.bridge_force_n);
    assert_eq!(
        near_log.body_volume_velocity_m3_s,
        far_log.body_volume_velocity_m3_s
    );
    assert!(
        near_log
            .radiated_pressure_pa
            .iter()
            .any(|p| p.abs() > 1.0e-12)
    );
    for (near_pa, far_pa) in near_log
        .radiated_pressure_pa
        .iter()
        .zip(&far_log.radiated_pressure_pa)
    {
        assert_eq!(*near_pa, 2.0 * far_pa);
    }

    // Positive finite distance alone cannot guarantee representable pressure.
    far.listener_m = f64::from_bits(1);
    assert!(matches!(run_bowed(&far), Err(BowedRunError::Radiation(_))));
}

/// G0: observer admission runs even for an empty simulation, before stepping.
#[test]
fn plate_one_port_refuses_invalid_listener_distance() {
    let mut cfg = config(playable_gesture(), FrictionIsland::Stribeck(rosin()), 0);
    cfg.termination = Termination::PlateOnePort {
        body: Box::new(
            CompactBody::from_radiator(RadiatingPlate {
                area_m2: 3.0e-3,
                mass_kg: 0.15,
                frequency_hz: 280.0,
                damping_ratio: 0.02,
            })
            .unwrap(),
        ),
        ambient: GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap(),
    };
    for distance in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        cfg.listener_m = distance;
        match run_bowed(&cfg) {
            Err(BowedRunError::InvalidListenerDistance { distance_m }) => {
                assert_eq!(distance_m.to_bits(), distance.to_bits());
            }
            result => panic!("invalid listener distance {distance:?} must refuse: {result:?}"),
        }
    }
    // A rigid termination has no radiation observer to admit.
    cfg.termination = Termination::Rigid;
    assert!(run_bowed(&cfg).is_ok());
}

/// Declared ambient pressure scales only the one-way radiation conversion;
/// the string/body trajectory is independent of it.
#[test]
fn plate_one_port_radiation_scales_with_declared_ambient_density() {
    let air = GasSpec::dry_air_ussa1976();
    let low_pressure = GasState::try_new(&air, 293.15, 101_325.0).expect("low pressure admits");
    let high_pressure = GasState::try_new(&air, 293.15, 202_650.0).expect("high pressure admits");
    let body_spec = RadiatingPlate {
        area_m2: 3.0e-3,
        mass_kg: 0.15,
        frequency_hz: 280.0,
        damping_ratio: 0.02,
    };
    let config_with = |ambient| BowedRunConfig {
        card: card(),
        island: FrictionIsland::Stribeck(rosin()),
        gesture: playable_gesture(),
        steps: 12_000,
        subsamples: 16,
        termination: Termination::PlateOnePort {
            body: Box::new(CompactBody::from_radiator(body_spec).expect("body spec admits")),
            ambient,
        },
        listener_m: 1.0,
    };
    let low = run_bowed(&config_with(low_pressure)).expect("low-pressure run completes");
    let high = run_bowed(&config_with(high_pressure)).expect("high-pressure run completes");

    assert_eq!(
        low.body_volume_velocity_m3_s,
        high.body_volume_velocity_m3_s
    );
    assert_eq!(low.radiated_pressure_pa.len(), 12_000);
    assert_eq!(high.radiated_pressure_pa.len(), 12_000);
    assert!(
        low.radiated_pressure_pa.iter().all(|p| p.is_finite())
            && high.radiated_pressure_pa.iter().all(|p| p.is_finite()),
        "both radiation traces must stay finite"
    );
    let expected_ratio = high_pressure.density / low_pressure.density;
    let mut nonzero_samples = 0;
    for (sample, (low_pressure_pa, high_pressure_pa)) in low
        .radiated_pressure_pa
        .iter()
        .zip(&high.radiated_pressure_pa)
        .enumerate()
    {
        if *low_pressure_pa == 0.0 && *high_pressure_pa == 0.0 {
            continue;
        }
        assert!(
            *low_pressure_pa != 0.0 && *high_pressure_pa != 0.0,
            "density scaling must not introduce a one-sided zero at sample {sample}"
        );
        let observed_ratio = high_pressure_pa / low_pressure_pa;
        assert!(
            (observed_ratio - expected_ratio).abs() <= 1e-12 * expected_ratio,
            "radiation ratio {observed_ratio} must match density ratio {expected_ratio} at sample {sample}"
        );
        nonzero_samples += 1;
    }
    assert!(
        nonzero_samples > 0,
        "bridge-driven body must radiate at least one nonzero pressure sample"
    );
}
