//! The 3-D slot-jet rig battery (bead frankensim-music-v8-root-3ez8g.10.1):
//! config refusals, determinism, fringe no-op law, momentum-exchange
//! sanity through a full smoke run, the shared-classifier wiring
//! (short records refuse the admitted jet band; an in-band tone reads
//! tonal, a steep broadband record broadband, roundoff fails the
//! amplitude floor), and the MEASURED 3-D fringe
use fs_aeroac::jetlab::transverse_force_spectrum;
use fs_aeroac::regime::{TONAL_FLATNESS_CEILING, measure_spectral_flatness};
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

fn checkpoint_dir(label: &str) -> std::path::PathBuf {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("test clock is after the Unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("sj-ckpt-{label}-{}-{nonce}", std::process::id()))
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
    cfg.steps_record = 63;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "record below the 64-sample minimum must be refused"
    );

    let mut cfg = smoke_config();
    cfg.seed_amplitude = 0.5;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "oversized seed must be refused"
    );

    let mut cfg = smoke_config();
    cfg.fringe_width = 0;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "zero-width fringe must be refused without panicking"
    );

    let mut cfg = smoke_config();
    cfg.fringe_sigma = 0.0;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "out-of-range fringe strength must be refused without panicking"
    );

    let mut cfg = smoke_config();
    cfg.nozzle_thickness = usize::MAX;
    assert!(
        run_slot_jet_3d(&cfg).is_err(),
        "overflowing axial layout must be refused without panicking"
    );
}

/// The refusal every record too short for the admitted jet band receives.
const SHORT_BAND_REFUSAL: &str = "force record too short";

/// A synthetic force record wrapped as a run of the smoke rig's geometry.
fn synthetic_run(cfg: &SlotJet3dConfig, series: Vec<[f64; 2]>) -> SlotJet3dRun {
    let n = series.len();
    SlotJet3dRun {
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
    }
}

/// `n` deterministic uniform values in [-1, 1) from a 64-bit LCG.
fn lcg_units(n: usize, seed: u64) -> Vec<f64> {
    let mut lcg = seed;
    (0..n)
        .map(|_| {
            lcg = lcg
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((lcg >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        })
        .collect()
}

#[test]
fn smoke_run_classifies_and_discloses_bin() {
    let mut cfg = smoke_config();
    cfg.steps_record = 67;
    let run = run_slot_jet_3d(&cfg).expect("smoke run must execute");
    assert_eq!(run.force_series.len(), 67);
    assert!(!run.scope.is_empty());
    #[allow(clippy::cast_precision_loss)]
    let bin_width = (1.0 / cfg.steps_record as f64) * 2.0 * cfg.slot_half / cfg.u_jet;
    assert!((run.diagnostics.strouhal_bin_width - bin_width).abs() < 1e-12);
    // A 67-step smoke holds no bin of the admitted jet band (St_h 1 is
    // bin 0.5): the verdict refuses by name instead of reading noise.
    let err = classify_rung(&run, &cfg).expect_err("a 67-step record must refuse");
    assert!(
        err.to_string().contains(SHORT_BAND_REFUSAL),
        "refusal must name the short band: {err}"
    );
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
    // 512 steps put 4 bins under St_h 1, so the admitted band holds one:
    // each record refuses by name whatever it carries. The verdict on
    // band-holding records is exercised by the two regressions below.
    for (label, series) in [("tone", tone), ("noise", noise), ("roundoff", roundoff)] {
        let err = classify_rung(&mk_run(series), &cfg).expect_err("a 512-step record must refuse");
        assert!(
            err.to_string().contains(SHORT_BAND_REFUSAL),
            "{label}: refusal must name the short band: {err}"
        );
    }
}

/// Regression: a tone whose period (400 steps) sits well inside the jet
/// band reads tonal at its own Strouhal. The pre-fix verdict skipped the
/// first n/8 bins, so on this record it could find no period above 8
/// steps (its peak landed at St 56).
#[test]
fn rung_verdict_reads_period_400_tone_tonal() {
    // 16384 steps at u_jet 0.04, slot_half 2.5: St_h 1 is bin 131, so
    // the admitted band holds 128 bins.
    let n = 16_384usize;
    let cfg = smoke_config();
    let floor = lcg_units(n, 0x243f_6a88_85a3_08d3);
    let tone_at = |amplitude: f64, floor_amplitude: f64| -> Vec<[f64; 2]> {
        (0..n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let t = i as f64;
                [
                    0.0,
                    amplitude * (2.0 * core::f64::consts::PI * t / 400.0).sin()
                        + floor_amplitude * floor[i],
                ]
            })
            .collect()
    };
    let rung = classify_rung(&synthetic_run(&cfg, tone_at(1e-3, 1e-9)), &cfg)
        .expect("a band-holding record classifies");
    assert!(
        rung.tonal && rung.flatness <= TONAL_FLATNESS_CEILING,
        "a period-400 tone must read tonal: flatness {}",
        rung.flatness
    );
    let expected = 2.0 * cfg.slot_half / cfg.u_jet / 400.0;
    assert!(
        (rung.strouhal - expected).abs() <= rung.strouhal_bin_width,
        "peak Strouhal {} must sit within a bin ({}) of the tone's {expected}",
        rung.strouhal,
        rung.strouhal_bin_width
    );
    assert!(rung.amplitude_qualified, "1e-3 tone clears the floor");
    let line = rung.to_jsonl();
    assert!(line.contains("\"schema\":\"fs-aeroac.slot-jet-3d.rung/v1\""));
    assert!(line.contains("\"strouhal_bin_width\":"));
    // The same tone at roundoff scale still fails the amplitude floor.
    let roundoff = classify_rung(&synthetic_run(&cfg, tone_at(1e-16, 0.0)), &cfg)
        .expect("a roundoff-scale record classifies");
    assert!(
        !roundoff.amplitude_qualified && roundoff.force_rms <= FORCE_RMS_AMPLITUDE_FLOOR,
        "roundoff-scale force must fail the amplitude floor"
    );
}

/// Regression: a resolved broadband record whose power falls steeply
/// (white noise through four one-pole low-passes at 0.004 cycles/step,
/// about 190 dB across the half-spectrum) reads broadband. Its flatness
/// to Nyquist, the pre-fix verdict, reads it tonal.
#[test]
fn rung_verdict_reads_steep_broadband_record_broadband() {
    let n = 16_384usize;
    let cfg = smoke_config();
    let a = (-2.0 * core::f64::consts::PI * 0.004).exp();
    let mut y: Vec<f64> = lcg_units(n, 0x1319_8a2e_0370_7344)
        .into_iter()
        .map(|u| 1e-3 * u)
        .collect();
    for _ in 0..4 {
        let mut prev = 0.0;
        for v in &mut y {
            prev = a * prev + (1.0 - a) * *v;
            *v = prev;
        }
    }
    let series: Vec<[f64; 2]> = y.iter().map(|&v| [0.0, v]).collect();
    let (power, _) = transverse_force_spectrum(&series, cfg.slot_half, cfg.u_jet, n / 8)
        .expect("spectrum of a valid record");
    let to_nyquist = measure_spectral_flatness(&power).expect("full-band flatness");
    assert!(
        to_nyquist <= TONAL_FLATNESS_CEILING,
        "the record must exercise the defect: flatness to Nyquist {to_nyquist}"
    );
    let rung = classify_rung(&synthetic_run(&cfg, series), &cfg)
        .expect("a band-holding record classifies");
    assert!(
        !rung.tonal && rung.flatness > TONAL_FLATNESS_CEILING,
        "a steeply falling broadband record must read broadband: flatness {}",
        rung.flatness
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
    let dir = checkpoint_dir("equiv");
    let first = run_slot_jet_3d_chunked(&cfg, &dir, 137).expect("chunk 1");
    assert!(matches!(
        first,
        fs_aeroac::slot_jet_3d::SweepProgress::Partial { .. }
    ));
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
        other @ fs_aeroac::slot_jet_3d::SweepProgress::Partial { .. } => {
            panic!("second chunk must complete, got {other:?}")
        }
    }
    assert!(
        dir.join("state.bin").is_file(),
        "completion must retain its replayable terminal checkpoint"
    );
    let replay = run_slot_jet_3d_chunked(&cfg, &dir, 1).expect("replay terminal checkpoint");
    let fs_aeroac::slot_jet_3d::SweepProgress::Complete(replayed) = replay else {
        panic!("terminal checkpoint must replay as complete");
    };
    assert_eq!(replayed.force_series.len(), whole.force_series.len());
    for (actual, expected) in replayed.force_series.iter().zip(&whole.force_series) {
        assert_eq!(actual[0].to_bits(), expected[0].to_bits());
        assert_eq!(actual[1].to_bits(), expected[1].to_bits());
    }
}

#[test]
fn checkpoint_refuses_foreign_configuration() {
    let cfg = smoke_config();
    let dir = checkpoint_dir("fp");
    let _ = run_slot_jet_3d_chunked(&cfg, &dir, 50).expect("seed chunk");
    let mut foreign = smoke_config();
    foreign.slot_half *= 2.0;
    let err =
        run_slot_jet_3d_chunked(&foreign, &dir, 50).expect_err("fingerprint mismatch must refuse");
    assert!(
        format!("{err}").contains("different configuration"),
        "refusal must name the fingerprint guard: {err}"
    );

    let mut foreign_record_length = smoke_config();
    foreign_record_length.steps_record *= 2;
    let err = run_slot_jet_3d_chunked(&foreign_record_length, &dir, 50)
        .expect_err("record-length changes must invalidate checkpoint identity");
    assert!(
        format!("{err}").contains("different configuration"),
        "fingerprint must cover the recording horizon: {err}"
    );
}

#[test]
fn checkpoint_reader_refuses_malformed_state_before_resume() {
    let cfg = smoke_config();
    let dir = checkpoint_dir("malformed");
    let _ = run_slot_jet_3d_chunked(&cfg, &dir, 50).expect("seed checkpoint");
    let path = dir.join("state.bin");
    let original = std::fs::read(&path).expect("read seeded checkpoint");

    let refuse = |mutated: &[u8], expected: &str| {
        std::fs::write(&path, mutated).expect("write malformed checkpoint fixture");
        let error = run_slot_jet_3d_chunked(&cfg, &dir, 1)
            .expect_err("malformed checkpoint must refuse before stepping");
        assert!(
            error.to_string().contains(expected),
            "expected {expected:?} refusal, got {error}"
        );
    };

    let mut unknown_version = original.clone();
    unknown_version[8..12].copy_from_slice(&u32::MAX.to_le_bytes());
    refuse(&unknown_version, "magic/version mismatch");

    let mut impossible_progress = original.clone();
    let settle_past_end = u64::try_from(cfg.steps_settle)
        .expect("smoke settle count fits u64")
        .checked_add(1)
        .expect("smoke settle count can advance once");
    impossible_progress[20..28].copy_from_slice(&settle_past_end.to_le_bytes());
    refuse(&impossible_progress, "progress is inconsistent");

    let mut impossible_force_count = original.clone();
    impossible_force_count[36..44].copy_from_slice(&u64::MAX.to_le_bytes());
    refuse(&impossible_force_count, "force count disagrees");

    let mut truncated = original.clone();
    truncated.pop();
    refuse(&truncated, "length does not match");

    let mut trailing = original.clone();
    trailing.push(0);
    refuse(&trailing, "length does not match");

    let mut non_finite_force = original.clone();
    non_finite_force[44..52].copy_from_slice(&f64::NAN.to_le_bytes());
    refuse(&non_finite_force, "non-finite force history");

    let record_len = usize::try_from(u64::from_le_bytes(
        original[36..44].try_into().expect("fixed record length"),
    ))
    .expect("smoke record length fits usize");
    let first_population = 44usize
        .checked_add(
            16usize
                .checked_mul(record_len)
                .expect("smoke force byte count fits usize"),
        )
        .expect("smoke population offset fits usize");
    let mut zero_density = original.clone();
    zero_density[first_population..first_population + 19 * 8].fill(0);
    refuse(&zero_density, "non-positive fluid density");
}
