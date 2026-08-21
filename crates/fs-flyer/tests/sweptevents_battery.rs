//! E3.4b / V-20 battery (bead wf-root-guzez.4.9): swept
//! point/segment/capsule event localization + phase-resolved
//! BladeCollisionProxyV1. Every V-20 hostile case EXECUTED:
//! disk-hits-but-blades-miss, hub-near-ground, high-Omega aliasing
//! (the naive endpoint check is run and shown blind), terrain-ridge,
//! event-bracketing. Cover certificate emitted with
//! BladeCollisionProxyArtifactId; generated radii (no hand-radius
//! API exists); hub exclusion; refinement 16 -> 24 -> typed refusal;
//! caps at cap AND beyond; determinism goldens (measure-then-pin).
//! Repro: cargo test -p fs-flyer --test sweptevents_battery

use fs_airscrew::{BladeStation, Rotor};
use fs_flyer::aircraft::wright_rotor_v1;
use fs_flyer::prelaunch::TerrainGrid;
use fs_flyer::simloop::TerminalEvent;
use fs_flyer::sweptevents::{
    BladeCollisionProxyV1, DISK_WARN_CLEARANCE_M, MARGIN_M, MAX_OMEGA_RAD_S, MAX_SEGMENT_FRAC,
    SweptOutcome, SweptPropMotion, build_blade_proxy, certs_at_resolutions, swept_feature_event,
    swept_prop_step,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-sweptevents\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 240.0;

/// 17x17 flat tile at 1 m spacing with an optional center spike.
fn tile(spike_m: f64) -> TerrainGrid {
    let mut rows = vec![vec![0.0f64; 17]; 17];
    rows[8][8] = spike_m;
    TerrainGrid::new(1.0, rows).unwrap()
}

fn constant_chord_rotor() -> Rotor {
    Rotor {
        radius_m: 1.2954,
        n_blades: 2,
        camber_ratio: 0.04,
        stations: vec![
            BladeStation {
                r_over_r: 0.30,
                chord_m: 0.15,
                beta_rad: 0.30,
            },
            BladeStation {
                r_over_r: 0.60,
                chord_m: 0.15,
                beta_rad: 0.30,
            },
            BladeStation {
                r_over_r: 0.95,
                chord_m: 0.15,
                beta_rad: 0.30,
            },
        ],
    }
}

/// Linear taper at slope 0.28 chord per r/R: measured excess 0.0106
/// at 16/blade (> the 0.0095 bound with 3 mm uncertainty) and 0.0086
/// at 24/blade — the refinement ladder's demonstrator.
fn taper_rotor() -> Rotor {
    Rotor {
        radius_m: 1.2954,
        n_blades: 2,
        camber_ratio: 0.04,
        stations: vec![
            BladeStation {
                r_over_r: 0.30,
                chord_m: 0.280,
                beta_rad: 0.30,
            },
            BladeStation {
                r_over_r: 0.60,
                chord_m: 0.196,
                beta_rad: 0.30,
            },
            BladeStation {
                r_over_r: 0.95,
                chord_m: 0.098,
                beta_rad: 0.30,
            },
        ],
    }
}

fn steep_chord_rotor() -> Rotor {
    Rotor {
        radius_m: 1.2954,
        n_blades: 2,
        camber_ratio: 0.04,
        stations: vec![
            BladeStation {
                r_over_r: 0.30,
                chord_m: 0.10,
                beta_rad: 0.30,
            },
            BladeStation {
                r_over_r: 0.60,
                chord_m: 0.55,
                beta_rad: 0.30,
            },
            BladeStation {
                r_over_r: 0.95,
                chord_m: 1.05,
                beta_rad: 0.30,
            },
        ],
    }
}

/// The registered rotor with its HONEST geometry uncertainty: E1.6 is
/// a reconstruction from the 1911 calibration tables — 1 cm class.
/// (Measured: tip-taper conservatism 0.0154 at 16/blade, 0.0118 at
/// 24/blade vs the max(5mm, unc)+margin bound of 0.0145 — the
/// refinement ladder is exercised by the real geometry.)
fn wright_proxy() -> BladeCollisionProxyV1 {
    build_blade_proxy(&wright_rotor_v1(), 0.010).unwrap()
}

#[test]
fn proxy_cover_certificate_and_artifact_id() {
    for (rotor, name) in [
        (constant_chord_rotor(), "const"),
        (wright_rotor_v1(), "wright"),
        (steep_chord_rotor(), "steep"),
    ] {
        for c in certs_at_resolutions(&rotor) {
            jlog(
                "cert-probe",
                &format!("\"rotor\":\"{name}\",\"cert\":\"{c:?}\""),
            );
        }
    }
    // Constant-chord rotor certifies at the 16/blade BASELINE.
    let slender = build_blade_proxy(&constant_chord_rotor(), 0.003).unwrap();
    assert_eq!(slender.certificate.capsules_per_blade, 16, "baseline");
    // The real 1903 rotor under its honest 1 cm reconstruction
    // uncertainty certifies at the baseline (measured excess 0.01175
    // vs the 0.0145 bound).
    let wright = wright_proxy();
    assert_eq!(wright.certificate.capsules_per_blade, 16, "wright");
    // The taper demonstrator EXERCISES the deterministic refinement:
    // 16/blade fails its 3 mm-uncertainty bound, 24/blade certifies.
    let tapered = build_blade_proxy(&taper_rotor(), 0.003).unwrap();
    assert_eq!(tapered.certificate.capsules_per_blade, 24, "refined");
    for proxy in [&slender, &wright, &tapered] {
        let cert = &proxy.certificate;
        assert!(cert.full_cover);
        assert!(cert.seg_len_max_frac <= MAX_SEGMENT_FRAC + 1e-12);
        // Hub void excluded: nothing reaches inboard of the first
        // registered station, and no capsule crosses the axis.
        assert!(cert.hub_min_r_frac >= 0.30 - 1e-12);
        for cap in &proxy.capsules {
            assert!(cap.r0_over_r >= 0.30 - 1e-12, "axis crossing");
            assert!(cap.r1_over_r > cap.r0_over_r);
            // Radii are GENERATED: strictly above the declared
            // margins, sane against the local chord.
            assert!(cap.radius_m > MARGIN_M);
            assert!(cap.radius_m < 0.30, "radius runaway: {}", cap.radius_m);
        }
    }
    // Artifact id: stable across rebuilds; geometry moves it.
    assert_eq!(wright.artifact_id, wright_proxy().artifact_id);
    let mut perturbed = wright_rotor_v1();
    perturbed.stations[2].chord_m += 0.001;
    assert_ne!(
        wright.artifact_id,
        build_blade_proxy(&perturbed, 0.010).unwrap().artifact_id,
        "geometry must move the artifact id"
    );
    // Adversarial steep chord: > 24 would be needed -> the blade
    // CLAIM refuses (typed) and only the disk warning remains.
    let err = build_blade_proxy(&steep_chord_rotor(), 0.003).unwrap_err();
    assert_eq!(err.code, "blade-cover-uncertifiable");
    // Uncertainty admission refusals.
    for bad in [f64::NAN, -1.0] {
        assert_eq!(
            build_blade_proxy(&wright_rotor_v1(), bad).unwrap_err().code,
            "blade-proxy-invalid"
        );
    }
    jlog(
        "cover",
        &format!(
            "\"baseline\":16,\"wright\":16,\"tapered\":24,\"wright_excess\":{},\"artifact_id\":\"{}\"",
            wright.certificate.excess_worst_m, wright.artifact_id
        ),
    );
}

#[test]
fn disk_hits_but_blades_miss() {
    // Spike under the hub tall enough to pierce the disk envelope;
    // both blades held horizontal (theta 0) and the phase barely
    // moves, so neither blade is anywhere near it.
    let proxy = wright_proxy();
    let terrain = tile(0.5);
    let motion = SweptPropMotion {
        hub0_m: [8.0, 8.0, 1.4],
        hub1_m: [8.0, 8.0, 1.4],
        theta0_rad: 0.0,
        omega_rad_s: 0.001,
        dt_s: DT,
    };
    let report = swept_prop_step(&wright_rotor_v1(), &proxy, &motion, &terrain).unwrap();
    assert_eq!(report.outcome, SweptOutcome::Clear, "no blade claim");
    assert!(report.outcome.terminal().is_none());
    let warn = report.disk_warning.expect("the disk envelope DOES hit");
    assert!(
        warn.min_clearance_m < 0.0,
        "envelope pierces: {}",
        warn.min_clearance_m
    );
    jlog(
        "disk-hits-blades-miss",
        &format!(
            "\"disk_clearance\":{},\"outcome\":\"clear\"",
            warn.min_clearance_m
        ),
    );
}

#[test]
fn hub_near_ground_is_the_hub_proxys_event() {
    // Hub descends onto a narrow pedestal; blades stay horizontal
    // over flat ground. The SEPARATE hub proxy owns this event — the
    // blade claim stays silent.
    let proxy = wright_proxy();
    let terrain = tile(0.45);
    let motion = SweptPropMotion {
        hub0_m: [8.0, 8.0, 0.62],
        hub1_m: [8.0, 8.0, 0.55],
        theta0_rad: 0.0,
        omega_rad_s: 0.001,
        dt_s: DT,
    };
    let report = swept_prop_step(&wright_rotor_v1(), &proxy, &motion, &terrain).unwrap();
    match report.outcome {
        SweptOutcome::HubStrike { t_event_s } => {
            assert!(t_event_s > 0.0 && t_event_s < DT, "localized inside");
            assert_eq!(
                report.outcome.terminal(),
                Some(TerminalEvent::DamageModelUnavailable)
            );
            jlog("hub-near-ground", &format!("\"t_event\":{t_event_s}"));
        }
        other => panic!("hub proxy must own this event, got {other:?}"),
    }
}

#[test]
fn high_omega_aliasing_caught_by_phase_resolved_sweep() {
    // Blades horizontal at BOTH step endpoints; omega spins the pair
    // through vertical mid-step, striking flat ground. The naive
    // endpoint check is EXECUTED first and is blind.
    let proxy = wright_proxy();
    let terrain = tile(0.0);
    let hub_z = 1.30;
    let omega = core::f64::consts::PI / DT; // theta sweeps exactly pi.
    // Naive endpoint check EXECUTED: a zero-omega sweep frozen at
    // each endpoint phase (blades horizontal both times) — blind.
    for theta0 in [0.0, core::f64::consts::PI] {
        let frozen = SweptPropMotion {
            hub0_m: [8.0, 8.0, hub_z],
            hub1_m: [8.0, 8.0, hub_z],
            theta0_rad: theta0,
            omega_rad_s: 0.0,
            dt_s: DT,
        };
        let r = swept_prop_step(&wright_rotor_v1(), &proxy, &frozen, &terrain).unwrap();
        assert_eq!(r.outcome, SweptOutcome::Clear, "endpoints ARE clear");
    }
    // Phase-resolved sweep catches the mid-step strike.
    let motion = SweptPropMotion {
        hub0_m: [8.0, 8.0, hub_z],
        hub1_m: [8.0, 8.0, hub_z],
        theta0_rad: 0.0,
        omega_rad_s: omega,
        dt_s: DT,
    };
    let report = swept_prop_step(&wright_rotor_v1(), &proxy, &motion, &terrain).unwrap();
    match report.outcome {
        SweptOutcome::BladeStrike {
            t_event_s,
            blade,
            capsule,
        } => {
            assert!(
                t_event_s > 0.3 * DT && t_event_s < 0.7 * DT,
                "strike near mid-step: {t_event_s}"
            );
            // Event bracketing: rerun bitwise-identically (fixed-count
            // bisection is deterministic).
            let again = swept_prop_step(&wright_rotor_v1(), &proxy, &motion, &terrain).unwrap();
            match again.outcome {
                SweptOutcome::BladeStrike { t_event_s: t2, .. } => {
                    assert_eq!(t_event_s.to_bits(), t2.to_bits(), "bit-identical twice");
                }
                other => panic!("determinism broke: {other:?}"),
            }
            assert_eq!(
                report.outcome.terminal(),
                Some(TerminalEvent::DamageModelUnavailable)
            );
            jlog(
                "high-omega",
                &format!(
                    "\"t_event\":{t_event_s},\"blade\":{blade},\"capsule\":{capsule},\"t_bits\":\"{:016x}\"",
                    t_event_s.to_bits()
                ),
            );
            // Golden over the localized time (measure-then-pin).
            assert_eq!(
                format!("{:016x}", t_event_s.to_bits()),
                "3f5ae4e93f983a76",
                "swept event time moved — determinism regression or an \
                 intentional kernel change requiring the golden-bump protocol"
            );
        }
        other => panic!("the swept pass must catch the mid-step strike, got {other:?}"),
    }
}

#[test]
fn terrain_ridge_and_feature_sweep() {
    // Temporal ridge: a skid segment glides over a spike — endpoints
    // of the STEP are clear, the middle of the step is not.
    let terrain = tile(0.8);
    let t = swept_feature_event(
        ([6.0, 7.0, 0.5], [6.0, 9.0, 0.5]),
        ([10.0, 7.0, 0.5], [10.0, 9.0, 0.5]),
        0.1,
        DT,
        &terrain,
    )
    .unwrap()
    .expect("ridge crossing must localize");
    assert!(t > 0.0 && t < DT);
    // Spatial ridge: the spike sits under the segment MIDDLE while
    // both endpoints fly over flat ground — interior axis samples
    // catch what endpoint-only sampling cannot.
    let t2 = swept_feature_event(
        ([8.0, 7.0, 1.0], [8.0, 9.0, 1.0]),
        ([8.0, 7.0, 0.6], [8.0, 9.0, 0.6]),
        0.05,
        DT,
        &terrain,
    )
    .unwrap()
    .expect("mid-segment ridge must localize");
    assert!(t2 > 0.0 && t2 < DT);
    // Bracketing witness: nudge the descent to stop just above the
    // spike — no event.
    let clear = swept_feature_event(
        ([8.0, 7.0, 1.0], [8.0, 9.0, 1.0]),
        ([8.0, 7.0, 0.90], [8.0, 9.0, 0.90]),
        0.05,
        DT,
        &terrain,
    )
    .unwrap();
    assert!(clear.is_none(), "stopping short must stay clear");
    // Refusals.
    assert_eq!(
        swept_feature_event(
            ([f64::NAN, 0.0, 1.0], [0.0, 1.0, 1.0]),
            ([0.0, 0.0, 1.0], [0.0, 1.0, 1.0]),
            0.1,
            DT,
            &terrain
        )
        .unwrap_err()
        .code,
        "swept-feature-invalid"
    );
    assert_eq!(
        swept_feature_event(
            ([0.0, 0.0, 1.0], [0.0, 1.0, 1.0]),
            ([0.0, 0.0, 1.0], [0.0, 1.0, 1.0]),
            -0.1,
            DT,
            &terrain
        )
        .unwrap_err()
        .code,
        "swept-feature-invalid"
    );
    // Outside-domain terrain refusals PROPAGATE (never invented).
    assert_eq!(
        swept_feature_event(
            ([-5.0, 0.0, 1.0], [-5.0, 1.0, 1.0]),
            ([-5.0, 0.0, 0.0], [-5.0, 1.0, 0.0]),
            0.1,
            DT,
            &terrain
        )
        .unwrap_err()
        .code,
        "terrain-query-outside-domain"
    );
    jlog(
        "terrain-ridge",
        &format!("\"t_temporal\":{t},\"t_spatial\":{t2}"),
    );
}

#[test]
fn motion_caps_at_cap_and_beyond() {
    let proxy = wright_proxy();
    let terrain = tile(0.0);
    let mk = |omega: f64, dt: f64| SweptPropMotion {
        hub0_m: [8.0, 8.0, 3.0],
        hub1_m: [8.0, 8.0, 3.0],
        theta0_rad: 0.0,
        omega_rad_s: omega,
        dt_s: dt,
    };
    // AT the omega cap admits; one part in 1e9 beyond refuses.
    assert!(
        swept_prop_step(
            &wright_rotor_v1(),
            &proxy,
            &mk(MAX_OMEGA_RAD_S, DT),
            &terrain
        )
        .is_ok()
    );
    assert_eq!(
        swept_prop_step(
            &wright_rotor_v1(),
            &proxy,
            &mk(MAX_OMEGA_RAD_S * (1.0 + 1e-9), DT),
            &terrain
        )
        .unwrap_err()
        .code,
        "swept-motion-invalid"
    );
    // Substep exhaustion: admitted motion whose phase rate demands
    // more than the cap refuses TYPED (never silently undersamples).
    assert_eq!(
        swept_prop_step(
            &wright_rotor_v1(),
            &proxy,
            &mk(MAX_OMEGA_RAD_S, 0.25),
            &terrain
        )
        .unwrap_err()
        .code,
        "swept-substeps-exhausted"
    );
    for bad in [mk(f64::NAN, DT), mk(100.0, 0.0), mk(100.0, -1.0)] {
        assert_eq!(
            swept_prop_step(&wright_rotor_v1(), &proxy, &bad, &terrain)
                .unwrap_err()
                .code,
            "swept-motion-invalid"
        );
    }
    // Disk warning threshold is a WARNING channel: hovering just
    // inside the threshold reports it without any outcome.
    let hover = mk(35.0, DT); // ~334 rpm, the historical band.
    let low = SweptPropMotion {
        hub0_m: [8.0, 8.0, 1.2954 + 0.02],
        hub1_m: [8.0, 8.0, 1.2954 + 0.02],
        ..hover
    };
    let report = swept_prop_step(&wright_rotor_v1(), &proxy, &low, &terrain).unwrap();
    let warn = report.disk_warning.expect("inside the warning band");
    assert!(warn.min_clearance_m > 0.0 && warn.min_clearance_m < DISK_WARN_CLEARANCE_M);
    jlog(
        "caps",
        &format!("\"warn_clearance\":{}", warn.min_clearance_m),
    );
}
