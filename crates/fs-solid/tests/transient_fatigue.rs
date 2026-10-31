//! G0/G1/G3/G4 native inertial stress histories and atomic fatigue observation.
//! The S-N curves are synthetic arithmetic references, not material life data.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::fatigue::{
    FatigueError, FatigueScope, MeanStressCorrection, RainflowLimits, RainflowMiner, SnCurve,
    SnPoint,
};
use fs_solid::fatigue::{TetFatigueError, TetFatigueGauge, TetFatigueHistoryConfig};
use fs_solid::linear3::*;
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};
use std::sync::Arc;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const FRAME: ContentHash = ContentHash([97; 32]);

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 97,
                kernel_id: 97,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(6.0, 60.0, 0.0, [97; 32]).unwrap()
}
fn thermal() -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        293.15,
        293.15,
        [0.0; 6],
        [97; 32],
        ContentHash([98; 32]),
        0.01,
    )
    .unwrap()
}
fn problem<'a>(
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::default(),
    }
}
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap()
}
fn initial() -> [[f64; 3]; 4] {
    let mut u = ZERO;
    u[1][0] = 0.001;
    u
}
fn curve(maximum: f64, temperature_k: f64) -> Arc<SnCurve> {
    Arc::new(
        SnCurve::new(
            FatigueScope {
                source_reference: "synthetic N(amplitude)=10/amplitude, not measured".into(),
                material: "synthetic oscillator".into(),
                process: "fixture".into(),
                surface_and_geometry: "smooth synthetic coupon".into(),
                environment: "fixture".into(),
                load_basis: "fully reversed normal stress amplitude".into(),
                temperature_k,
            },
            vec![
                SnPoint {
                    amplitude_pa: 1e-6,
                    cycles_to_failure: 1e7,
                },
                SnPoint {
                    amplitude_pa: maximum,
                    cycles_to_failure: 10.0 / maximum,
                },
            ],
            Some(1e-8),
        )
        .unwrap(),
    )
}
fn gauge(normal: [f64; 3], maximum: f64) -> TetFatigueGauge {
    TetFatigueGauge {
        element: 0,
        normal_in_mesh_frame: normal,
        material_state_identity: [97; 32],
        curve: curve(maximum, 293.15),
        mean_stress: MeanStressCorrection::Uncorrected,
    }
}
fn config() -> TetFatigueHistoryConfig {
    TetFatigueHistoryConfig {
        maximum_absolute_mandel_strain: 0.01,
        history: RainflowLimits {
            maximum_samples: 1000,
            maximum_residue: 32,
        },
        maximum_gauges: 4,
        maximum_total_residue: 128,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_unforced_vibration_counts_real_inertial_stresses_and_converges_in_time() {
    with_cx(false, |cx| {
        // Independent one-DOF model: V=1/6, K=E/6=10, M=rho/60=.1, omega=10.
        // A quasistatic zero-load replacement would return zero stress/damage.
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut period_errors = Vec::new();
        for per_cycle in [16, 32, 64] {
            // Newmark rotates by 2 atan(omega*h/2): put a full cycle on this
            // sample grid while refining its physical period towards 2*pi/omega.
            let h = (std::f64::consts::PI / f64::from(per_cycle)).tan() / 5.0;
            let mut native = load
                .begin_transient(
                    &problem,
                    Some(TetThermalStrainField::Uniform(&thermal)),
                    0.0,
                    &initial(),
                    &ZERO,
                    TransientConfig3::default(),
                    cx,
                )
                .unwrap();
            let gauges = vec![gauge([1.0, 0.0, 0.0], 1.0), gauge([1.0, 1.0, 0.0], 1.0)];
            let mut independent: Vec<_> = gauges
                .iter()
                .map(|g| {
                    RainflowMiner::new(g.curve.clone(), g.mean_stress, config().history).unwrap()
                })
                .collect();
            for (counter, factor) in independent.iter_mut().zip([1.0, 0.5]) {
                counter.push(0.06 * factor).unwrap();
            }
            let mut run = native
                .clone()
                .with_fatigue(FRAME, gauges, config(), cx)
                .unwrap();
            assert!(run.summaries().iter().all(|s| s.samples == 1));
            let initial_energy = native.state().total_energy_j;
            for step in 1..=2 * per_cycle {
                let original_step = native.advance(&load, h, cx).unwrap();
                let observed = run.advance(&load, h, cx).unwrap();
                assert_eq!(observed.mechanics, original_step);
                assert_eq!(run.mechanics().state(), native.state());
                assert_eq!(
                    observed.response,
                    native.recover_response(FRAME, 0.01, cx).unwrap()
                );
                assert!(
                    observed
                        .response
                        .external_force_n
                        .iter()
                        .flatten()
                        .all(|f| *f == 0.0)
                );
                close(observed.mechanics.quadratic_energy_j, initial_energy, 1e-16);
                let analytic = 0.06
                    * (2.0 * std::f64::consts::PI * f64::from(step) / f64::from(per_cycle)).cos();
                let signed = 60.0 * native.state().displacement[0];
                close(signed, analytic, 2e-14);
                for (i, factor) in [1.0, 0.5].into_iter().enumerate() {
                    close(observed.signed_stress_pa[i], signed * factor, 1e-14);
                    independent[i].push(signed * factor).unwrap();
                    assert_eq!(
                        observed.fatigue[i].samples,
                        u64::try_from(step + 1).unwrap()
                    );
                }
            }
            let fatigue = run.finish(cx).unwrap();
            for (i, expected_damage) in [0.012, 0.006].into_iter().enumerate() {
                close(fatigue[i].miner_damage, expected_damage, 1e-14);
                close(
                    fatigue[i].miner_damage,
                    independent[i].finish().unwrap().miner_damage,
                    1e-14,
                );
            }
            period_errors
                .push((run.mechanics().state().time_s - 4.0 * std::f64::consts::PI / 10.0).abs());
        }
        for pair in period_errors.windows(2) {
            assert!(
                (3.8..4.2).contains(&(pair[0] / pair[1])),
                "{period_errors:?}"
            );
        }
    });
}

#[test]
fn g3_processing_chunks_and_checkpoint_forks_preserve_open_cycles_bitwise() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let native = load
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial(),
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut whole = native
            .with_fatigue(FRAME, vec![gauge([1.0, 0.0, 0.0], 1.0)], config(), cx)
            .unwrap();
        for _ in 0..19 {
            whole.advance(&load, 0.01, cx).unwrap();
        }
        let mut chunked = whole.clone();
        for _ in 0..127 {
            whole.advance(&load, 0.01, cx).unwrap();
        }
        for chunk in [3, 17, 59, 48] {
            for _ in 0..chunk {
                chunked.advance(&load, 0.01, cx).unwrap();
            }
            let before = chunked.summaries();
            chunked.recover_response(cx).unwrap();
            chunked.recover_response(cx).unwrap();
            assert_eq!(chunked.summaries(), before); // observing is not sampling again
        }
        assert_eq!(chunked.mechanics().state(), whole.mechanics().state());
        assert_eq!(chunked.summaries(), whole.summaries());
        let final_whole = whole.finish(cx).unwrap();
        assert_eq!(chunked.finish(cx).unwrap(), final_whole);
        assert_eq!(chunked.finish(cx).unwrap(), final_whole);
        let state = chunked.mechanics().state().clone();
        assert!(matches!(
            chunked.advance(&load, 0.01, cx),
            Err(TetFatigueError::Finished)
        ));
        assert_eq!(chunked.mechanics().state(), &state);
    });
}

#[test]
fn g4_a_later_gauge_refusal_rolls_back_the_entire_native_mechanical_step() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let native = load
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial(),
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut run = native
            .with_fatigue(
                FRAME,
                vec![gauge([1.0, 0.0, 0.0], 1.0), gauge([1.0, 0.0, 0.0], 0.03)],
                config(),
                cx,
            )
            .unwrap();
        let mut rejected = false;
        for _ in 0..10 {
            let state = run.mechanics().state().clone();
            let response = run.recover_response(cx).unwrap();
            let summaries = run.summaries();
            match run.advance(&load, 0.2, cx) {
                Ok(_) => continue,
                Err(TetFatigueError::Gauge {
                    gauge: 1,
                    error: FatigueError::AmplitudeOutsideCurve { .. },
                }) => {
                    assert_eq!(run.mechanics().state(), &state);
                    assert_eq!(run.recover_response(cx).unwrap(), response);
                    assert_eq!(run.summaries(), summaries);
                    assert!(matches!(
                        run.finish(cx),
                        Err(TetFatigueError::Gauge { gauge: 1, .. })
                    ));
                    assert_eq!(run.summaries(), summaries);
                    assert!(run.summaries().iter().all(|s| !s.finished));
                    rejected = true;
                    break;
                }
                Err(error) => panic!("unexpected refusal: {error}"),
            }
        }
        assert!(rejected);
    });
}

#[test]
fn g4_cancellation_foreign_loads_and_budgets_preserve_the_joint_checkpoint() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let native = load
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial(),
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut cfg = config();
        cfg.history.maximum_samples = 2;
        let mut run = native
            .with_fatigue(
                FRAME,
                vec![gauge([1.0, 0.0, 0.0], 1.0), gauge([1.0, 1.0, 0.0], 1.0)],
                cfg,
                cx,
            )
            .unwrap();
        let state = run.mechanics().state().clone();
        let history = run.summaries();
        with_cx(true, |stopped| {
            assert!(matches!(
                run.advance(&load, 0.01, stopped),
                Err(TetFatigueError::Cancelled)
            ));
            assert!(matches!(
                run.finish(stopped),
                Err(TetFatigueError::Cancelled)
            ));
        });
        let other_nodes = NODES.map(|p| p.map(|x| x * 2.0));
        let other_geometry =
            TetLoadGeometry3::new(&other_nodes, &TETS, problem.budget, 4, cx).unwrap();
        let foreign = other_geometry
            .assemble(&[], VolumeLoad3::None, &[], cx)
            .unwrap();
        assert!(matches!(
            run.advance(&foreign, 0.01, cx),
            Err(TetFatigueError::Transient(TransientError3::Invalid(_)))
        ));
        assert!(run.advance(&load, 0.0, cx).is_err());
        assert_eq!(run.mechanics().state(), &state);
        assert_eq!(run.summaries(), history);
        run.advance(&load, 0.01, cx).unwrap();
        let state = run.mechanics().state().clone();
        let history = run.summaries();
        assert!(matches!(
            run.advance(&load, 0.01, cx),
            Err(TetFatigueError::Budget)
        ));
        assert_eq!(run.mechanics().state(), &state);
        assert_eq!(run.summaries(), history);
    });
}

#[test]
fn g0_temperature_source_frame_and_initial_history_are_bound_before_sampling() {
    with_cx(false, |cx| {
        let material = material();
        let mut source = [thermal()];
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut native = load
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::PerElement(&source)),
                0.0,
                &initial(),
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut wrong_source = gauge([1.0, 0.0, 0.0], 1.0);
        wrong_source.material_state_identity = [99; 32];
        let mut wrong_temperature = gauge([1.0, 0.0, 0.0], 1.0);
        wrong_temperature.curve = curve(1.0, 350.0);
        for wrong in [wrong_source, wrong_temperature] {
            assert!(matches!(
                native
                    .clone()
                    .with_fatigue(FRAME, vec![wrong], config(), cx),
                Err(TetFatigueError::SourceMismatch { gauge: 0, .. })
            ));
        }
        assert!(
            native
                .clone()
                .with_fatigue(
                    ContentHash([0; 32]),
                    vec![gauge([1.0, 0.0, 0.0], 1.0)],
                    config(),
                    cx
                )
                .is_err()
        );
        let no_temperature = load
            .begin_transient(
                &problem,
                None,
                0.0,
                &initial(),
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            no_temperature.with_fatigue(FRAME, vec![gauge([1.0, 0.0, 0.0], 1.0)], config(), cx),
            Err(TetFatigueError::Invalid(_))
        ));
        // The actual native field was frozen by begin_transient; a caller's
        // later source mutation cannot change the calibration binding.
        source[0] = TetThermalStrainState::try_new(
            350.0,
            293.15,
            [0.001; 6],
            [97; 32],
            ContentHash([98; 32]),
            0.01,
        )
        .unwrap();
        let run = native
            .clone()
            .with_fatigue(FRAME, vec![gauge([1.0, 0.0, 0.0], 1.0)], config(), cx)
            .unwrap();
        assert_eq!(run.curve(0).unwrap().scope().temperature_k, 293.15);
        close(
            run.recover_response(cx).unwrap().stress.elements()[0].stress_mandel_pa()[0],
            0.06,
            1e-14,
        );
        native.advance(&load, 0.01, cx).unwrap();
        assert!(matches!(
            native.with_fatigue(FRAME, vec![gauge([1.0, 0.0, 0.0], 1.0)], config(), cx),
            Err(TetFatigueError::Invalid(_))
        ));
    });
}

#[test]
fn g1_prescribed_support_cycles_keep_native_inertia_work_and_stress() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut native = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &ZERO,
                &ZERO,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let gauge = gauge([1.0, 0.0, 0.0], 1.0);
        let mut independent =
            RainflowMiner::new(gauge.curve.clone(), gauge.mean_stress, config().history).unwrap();
        independent.push(0.0).unwrap();
        let mut run = native
            .clone()
            .with_fatigue(FRAME, vec![gauge], config(), cx)
            .unwrap();
        assert!(run.advance(&load, 0.05, cx).is_err());
        assert_eq!(run.summaries()[0].samples, 1);
        for displacement in [0.001, -0.001, 0.001, -0.001, 0.0] {
            let original = native
                .advance_with_supports(&load, &[(3, displacement)], 0.05, cx)
                .unwrap();
            let step = run
                .advance_with_supports(&load, &[(3, displacement)], 0.05, cx)
                .unwrap();
            assert_eq!(step.mechanics, original);
            assert_eq!(
                step.response,
                native.recover_response(FRAME, 0.01, cx).unwrap()
            );
            close(step.signed_stress_pa[0], 60.0 * displacement, 1e-14);
            independent.push(60.0 * displacement).unwrap();
        }
        close(
            run.finish(cx).unwrap()[0].miner_damage,
            independent.finish().unwrap().miner_damage,
            1e-14,
        );
    });
}

#[test]
fn g4_stress_recovery_domain_refusal_does_not_commit_native_motion() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let initial_load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let next_load = geometry
            .assemble(&[], VolumeLoad3::None, &[(3, 0.04)], cx)
            .unwrap();
        let native = initial_load
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial(),
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut cfg = config();
        cfg.maximum_absolute_mandel_strain = 0.0015;
        let mut run = native
            .with_fatigue(FRAME, vec![gauge([1.0, 0.0, 0.0], 1.0)], cfg, cx)
            .unwrap();
        let before = run.mechanics().state().clone();
        let history = run.summaries();
        let mut mechanical_only = run.mechanics().clone();
        mechanical_only.advance(&next_load, 0.2, cx).unwrap();
        assert!(mechanical_only.state().displacement[0] > cfg.maximum_absolute_mandel_strain);
        assert!(matches!(
            run.advance(&next_load, 0.2, cx),
            Err(TetFatigueError::Transient(_))
        ));
        assert_eq!(run.mechanics().state(), &before);
        assert_eq!(run.summaries(), history);
    });
}
