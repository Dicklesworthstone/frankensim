//! The 3-D slot-jet rig battery (bead frankensim-music-v8-root-3ez8g.10.1):
//! config refusals, determinism, fringe no-op law, momentum-exchange
//! sanity through a full smoke run, the shared-classifier wiring
//! (tone vs noise vs roundoff), and the MEASURED 3-D fringe
use fs_aeroac::regime::TONAL_FLATNESS_CEILING;
use fs_aeroac::slot_jet_3d::{
    FORCE_RMS_AMPLITUDE_FLOOR, Fringe3, SlotJet3dConfig, SlotJet3dDiagnostics, SlotJet3dRun,
    classify_rung, run_slot_jet_3d, run_slot_jet_3d_chunked,
};
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, CollisionModel3, FaceBoundary3};

/// Small admissible config for short smoke runs.
fn smoke_config() -> SlotJet3dConfig {
    SlotJet3dConfig {
        nx: 64,
        ny: 32,
        nz: 8,
        slot_half: 2.5,
        u_jet: 0.04,
        collision: CollisionModel3::CentralMoment {
            second_order_rate: 1.8,
            higher_order_rate: 1.9,
        },
        nozzle_thickness: 1,
        edge_distance: 8,
        plate_length: 2,
        fringe_width: 8,
        fringe_sigma: 0.4,
        seed_amplitude: 0.05,
        steps_settle: 20,
        steps_record: 64,
    }
}

#[test]
fn config_refusals() {
    let mut cfg = smoke_config();
    cfg.collision = CollisionModel3::Bgk { tau: 0.9 };
    assert!(run_slot_jet_3d(&cfg).is_err(), "BGK must be refused");

    let mut cfg = smoke_config();
    cfg.nz = 6;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "non-tile-aligned span must be refused"
    );

    let mut cfg = smoke_config();
    cfg.nozzle_thickness = 0;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "no-nozzle rig must be refused"
    );

    let mut cfg = smoke_config();
    cfg.fringe_width = 64;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "overflowing fringe must be refused"
    );

    let mut cfg = smoke_config();
    cfg.steps_record = 100;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "non-power-of-two record must be refused"
    );

    let mut cfg = smoke_config();
    cfg.seed_amplitude = 0.5;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "oversized seed must be refused"
    );
}

#[test]
fn smoke_run_classifies_and_discloses_bin() {
    let cfg = smoke_config();
    let run = run_slot_jet_3d(&cfg).expect("smoke run must execute");
    assert_eq!(run.force_series.len(), 64);
    assert!(!run.scope.is_empty());
    #[allow(clippy::cast_precision_loss)]
    let bin_width = (1.0f64 / 64.0) * 2.0 * cfg.slot_half / cfg.u_jet;
    assert!((run.diagnostics.strouhal_bin_width - bin_width).abs() < 1e-12);
    let rung = classify_rung(&run, &cfg).expect("classification must succeed");
    // A 64-step smoke cannot carry physics claims either way; the
    // receipt must still be fully populated.
    assert!(rung.flatness.is_finite());
    assert_eq!(
        rung.amplitude_qualified,
        rung.force_rms > FORCE_RMS_AMPLITUDE_FLOOR
    );
    let line = rung.to_jsonl();
    assert!(line.contains("\"schema\":\"fs-aeroac.slot-jet-3d.rung/v1\""));
    assert!(line.contains("\"strouhal_bin_width\":"));
}

#[test]
fn determinism_bitwise() {
    let cfg = smoke_config();
    let a = run_slot_jet_3d(&cfg).expect("run a");
    let b = run_slot_jet_3d(&cfg).expect("run b");
    assert_eq!(a.force_series.len(), b.force_series.len());
    for (fa, fb) in a.force_series.iter().zip(&b.force_series) {
        assert_eq!(fa[0].to_bits(), fb[0].to_bits());
        assert_eq!(fa[1].to_bits(), fb[1].to_bits());
    }
}

#[test]
fn fringe_is_noop_on_matched_state() {
    let mut grid = BoundaryGrid3::with_collision_model(
        32,
        16,
        8,
        CollisionModel3::CentralMoment {
            second_order_rate: 1.8,
            higher_order_rate: 1.9,
        },
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    let profile: Vec<(f64, [f64; 3])> = (0..16).map(|_| (1.0, [0.0, 0.0, 0.0])).collect();
    let fringe = Fringe3::with_profile(4, 0.5, &profile);
    let probe_cells: Vec<(usize, usize)> =
        (0..8).flat_map(|z| (0..16).map(move |y| (y, z))).collect();
    let snapshot: Vec<[f64; 19]> = probe_cells
        .iter()
        .map(|&(y, z)| grid.populations(3, y, z))
        .collect();
    fringe.apply(&mut grid);
    for (n, &(y, z)) in probe_cells.iter().enumerate() {
        let after = grid.populations(3, y, z);
        for q in 0..19 {
            assert!(
                (after[q] - snapshot[n][q]).abs() < 1.0e-14,
                "matched-state fringe must be a no-op"
            );
        }
    }
}

#[test]
fn classifier_tone_vs_noise_vs_roundoff() {
    // Deterministic synthetic records through the SAME pipeline the
    // rig feeds: a pure tone, seeded broadband noise, and a
    // roundoff-scale structured signal (the vacuous-oscillation trap).
    let n = 512usize;
    let cfg = smoke_config();
    let mk_run = |series: Vec<[f64; 2]>| SlotJet3dRun {
        force_series: series,
        diagnostics: SlotJet3dDiagnostics {
            mach_max_lattice: 0.05,
            flux_plate_plane: 0.01,
            flux_fringe_plane: 0.0099,
            reynolds: 100.0,
            record_len: n,
            #[allow(clippy::cast_precision_loss)]
            strouhal_bin_width: (1.0 / n as f64) * 2.0 * cfg.slot_half / cfg.u_jet,
        },
        scope: "test",
    };
    let mut tone: Vec<[f64; 2]> = Vec::with_capacity(n);
    for i in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let t = i as f64;
        tone.push([
            0.0,
            1e-3 * (2.0 * core::f64::consts::PI * 10.0 * t / n as f64).sin(),
        ]);
    }
    let mut lcg: u64 = 0x243f_6a88_85a3_08d3;
    let mut noise: Vec<[f64; 2]> = Vec::with_capacity(n);
    for _ in 0..n {
        lcg = lcg
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let unit = ((lcg >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0;
        noise.push([0.0, 1e-3 * unit]);
    }
    let mut roundoff: Vec<[f64; 2]> = Vec::with_capacity(n);
    for i in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let t = i as f64;
        roundoff.push([
            0.0,
            1e-16 * (2.0 * core::f64::consts::PI * 12.0 * t / n as f64).sin(),
        ]);
    }
    let tone_rung = classify_rung(&mk_run(tone), &cfg).expect("tone classifies");
    assert!(tone_rung.tonal, "pure tone must classify tonal");
    assert!(tone_rung.amplitude_qualified, "1e-3 tone clears the floor");
    let noise_rung = classify_rung(&mk_run(noise), &cfg).expect("noise classifies");
    assert!(!noise_rung.tonal, "white noise must NOT classify tonal");
    assert!(noise_rung.flatness > TONAL_FLATNESS_CEILING);
    let roundoff_rung = classify_rung(&mk_run(roundoff), &cfg).expect("roundoff classifies");
    assert!(
        !roundoff_rung.amplitude_qualified,
        "roundoff-scale force must fail the amplitude floor"
    );
}

/// MEASURED 3-D fringe reflection coefficient against a bounce-back
/// wall control (the port of the fs-lbm sponge battery recipe onto
/// BoundaryGrid3 + Fringe3). A Gaussian density pulse splits into two
/// acoustic packets at c_s = 1/sqrt(3); a probe between source and
/// layer records the incident and (later) reflected maxima at the
/// SAME location.
#[test]
fn fringe_reflection_measured_with_wall_control() {
    const NX: usize = 256;
    const NY: usize = 4;
    const NZ: usize = 4;
    const PULSE_X: usize = 64;
    const PROBE_X: usize = 112;
    const FRINGE_W: usize = 48;
    const TAU: f64 = 0.55;
    let pulse_grid = |solid_slab: Option<usize>| {
        let mut grid = BoundaryGrid3::new(
            NX,
            NY,
            NZ,
            TAU,
            [0.0; 3],
            BoundarySpec3::new([FaceBoundary3::Periodic; 6]),
        );
        if let Some(sx0) = solid_slab {
            grid.voxelize_sdf(move |[sx, _, _]| {
                if sx >= sx0 as f64 && sx < (sx0 + 8) as f64 {
                    -1.0
                } else {
                    1.0
                }
            });
        }
        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    if grid.is_solid(x, y, z) {
                        continue;
                    }
                    #[allow(clippy::cast_precision_loss)]
                    let dx = (x as f64 - PULSE_X as f64) / 4.0;
                    let rho = 1.0 + 1.0e-4 * (-dx * dx).exp();
                    let f = fs_lbm::d3q19::equilibrium3(rho, [0.0; 3]);
                    grid.set_populations(x, y, z, &f);
                }
            }
        }
        grid
    };
    let profile: Vec<(f64, [f64; 3])> = (0..NY).map(|_| (1.0, [0.0, 0.0, 0.0])).collect();
    let fringe = Fringe3::with_profile(FRINGE_W, 0.8, &profile);
    let measure = |mut grid: BoundaryGrid3, apply: bool| {
        let mut incident = 0.0f64;
        let mut reflected = 0.0f64;
        for t in 0..1000 {
            grid.step();
            if apply {
                fringe.apply(&mut grid);
            }
            let a = (grid.density(PROBE_X, 1, 1) - 1.0).abs();
            // Outbound packet passes the probe around t ~ 48*sqrt(3) ~ 83;
            // the layer inner edge is x = 208, so a reflection returns
            // around t ~ (48 + 2*(208-112))*sqrt(3) ~ 407.
            if (60..220).contains(&t) {
                incident = incident.max(a);
            }
            if (330..900).contains(&t) {
                reflected = reflected.max(a);
            }
        }
        reflected / incident
    };
    let r_fringe = measure(pulse_grid(None), true);
    let r_wall = measure(pulse_grid(Some(208)), false);
    assert!(
        r_fringe < 5.0e-2,
        "3-D fringe reflection must stay under the authored ceiling: {r_fringe:.4}"
    );
    assert!(
        r_wall > 10.0 * r_fringe.max(1e-6),
        "wall control must reflect far more than the fringe: wall={r_wall:.4} fringe={r_fringe:.4}"
    );
}

#[test]
fn chunked_resume_bitwise_matches_whole_run() {
    let mut cfg = smoke_config();
    cfg.steps_settle = 200;
    cfg.steps_record = 256;
    let whole = run_slot_jet_3d(&cfg).expect("whole run");

    // Unique per-process checkpoint dir; deliberately never deleted
    // (repository law: no file deletion without explicit approval).
    let dir = std::env::temp_dir().join(format!(
        "sj-ckpt-equiv-{}-{}",
        std::process::id(),
        cfg.steps_settle
    ));
    let first = run_slot_jet_3d_chunked(&cfg, &dir, 137).expect("chunk 1");
    assert!(matches!(first, fs_aeroac::slot_jet_3d::SweepProgress::Partial { .. }));
    let second = run_slot_jet_3d_chunked(&cfg, &dir, 10_000).expect("chunk 2");
    match second {
        fs_aeroac::slot_jet_3d::SweepProgress::Complete(run) => {
            assert_eq!(run.force_series.len(), whole.force_series.len());
            for (fa, fb) in run.force_series.iter().zip(&whole.force_series) {
                assert_eq!(fa[0].to_bits(), fb[0].to_bits());
                assert_eq!(fa[1].to_bits(), fb[1].to_bits());
            }
            assert_eq!(
                run.diagnostics.reynolds.to_bits(),
                whole.diagnostics.reynolds.to_bits()
            );
        }
        other => panic!("second chunk must complete, got {other:?}"),
    }
}

#[test]
fn checkpoint_refuses_foreign_configuration() {
    let cfg = smoke_config();
    let dir = std::env::temp_dir().join(format!("sj-ckpt-fp-{}", std::process::id()));
    let _ = run_slot_jet_3d_chunked(&cfg, &dir, 50).expect("seed chunk");
    let mut foreign = smoke_config();
    foreign.slot_half *= 2.0;
    let err = run_slot_jet_3d_chunked(&foreign, &dir, 50)
        .expect_err("fingerprint mismatch must refuse");
    assert!(
        format!("{err}").contains("different configuration"),
        "refusal must name the fingerprint guard: {err}"
    );
}
