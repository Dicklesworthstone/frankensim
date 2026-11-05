//! G0/G1/G4: actual native stress cycles, analytical growth and atomic refusal.
//! Every coefficient is a synthetic implementation fixture, not measured life data.

use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::FatigueScope;
use fs_material::fatigue::crack_growth::source::{
    ParisProperty, ParisSourceRequest, STRESS_INTENSITY_SQUARED_DIMS, resolve_paris,
};
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGrowthError, CrackGrowthLimits, ParisDomain, ParisLaw,
    ParisParameters,
};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use fs_solid::fatigue::TetFatigueError;
use fs_solid::fatigue::crack_growth::{
    CrackStressTolerance, TetCrackGrowthConfig, TetCrackGrowthError, TetCrackGrowthGauge,
    TetCrackGrowthPath,
};
use fs_solid::linear3::source::SolidStiffnessModel;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState,
};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::{TetLoadCase3, TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const FRAME: ContentHash = ContentHash([104; 32]);
const INITIAL_A: f64 = 1e-3;

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
                seed: 104,
                kernel_id: 104,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(7800.0, 200e9, 0.0, [104; 32]).unwrap()
}
fn thermal() -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        300.0,
        300.0,
        [0.0; 6],
        [104; 32],
        ContentHash([105; 32]),
        0.01,
    )
    .unwrap()
}
fn fixed() -> Vec<usize> {
    (0..12).filter(|d| *d != 3).collect()
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
fn load<'g, 'm>(
    geometry: &'g TetLoadGeometry3<'m>,
    stress: f64,
    cx: &Cx<'_>,
) -> TetLoadCase3<'g, 'm> {
    // Volume = 1/6, strain = displacement at node 1, stiffness = E/6.
    geometry
        .assemble(&[], VolumeLoad3::None, &[(3, stress / 6.0)], cx)
        .unwrap()
}
fn law(minimum_delta_k: f64, maximum_k: f64) -> Arc<ParisLaw> {
    Arc::new(
        ParisLaw::new(
            FatigueScope {
                source_reference: "synthetic Paris m=2 fixture; not measured".into(),
                material: "synthetic elastic solid".into(),
                process: "fixture".into(),
                surface_and_geometry:
                    "declared infinite-plate nominal stress; a is half crack length".into(),
                environment: "synthetic fixed conditions".into(),
                load_basis: "tensile Mode-I full stress-intensity range".into(),
                temperature_k: 300.0,
            },
            ParisParameters {
                reference_delta_k_pa_sqrt_m: 1e6,
                reference_rate_m_per_cycle: 1e-4,
                exponent: 2.0,
            },
            ParisDomain {
                minimum_delta_k_pa_sqrt_m: minimum_delta_k,
                maximum_delta_k_pa_sqrt_m: 1e8,
                minimum_load_ratio: 0.0,
                maximum_load_ratio: 0.95,
                maximum_k_pa_sqrt_m: maximum_k,
            },
        )
        .unwrap(),
    )
}
fn gauge() -> TetCrackGrowthGauge {
    TetCrackGrowthGauge {
        element: 0,
        normal_in_mesh_frame: [1.0, 0.0, 0.0],
        material_state_identity: [104; 32],
        material_source: None,
        law: law(1.0, 1e9),
        geometry: ConstantGeometryCrack {
            geometry_factor: 1.0,
            minimum_length_m: 1e-5,
            maximum_length_m: 0.01,
        },
        initial_length_m: INITIAL_A,
        limits: CrackGrowthLimits {
            maximum_blocks: 100,
            maximum_cycles: 100.0,
        },
        closure: CrackStressTolerance {
            absolute_stress_pa: 1e-5,
            relative: 1e-9,
        },
        mode_i: CrackStressTolerance {
            absolute_stress_pa: 1e-5,
            relative: 1e-9,
        },
    }
}
fn config() -> TetCrackGrowthConfig {
    TetCrackGrowthConfig {
        maximum_absolute_mandel_strain: 0.01,
        maximum_gauges: 4,
        maximum_samples: 1000,
    }
}
fn path<'p, 'm>(
    problem: &'p TetLinearElasticProblem<'m>,
    thermal: &'p TetThermalStrainState,
    gauges: Vec<TetCrackGrowthGauge>,
    cx: &Cx<'_>,
) -> TetCrackGrowthPath<'p, 'm> {
    TetCrackGrowthPath::new(
        problem,
        TetThermalStrainField::Uniform(thermal),
        FRAME,
        gauges,
        Default::default(),
        config(),
        cx,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual:e} != {expected:e}"
    );
}
fn analytical_length(initial: f64, stress_range: f64) -> f64 {
    // Independently integrated da/dN = (1e-4/1e12) pi (delta_sigma)^2 a.
    initial * (1e-16 * std::f64::consts::PI * stress_range * stress_range).exp()
}

#[test]
fn g1_native_nominal_cycles_match_analytical_growth_without_changing_mechanics() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut path = path(&problem, &thermal, vec![gauge()], cx);
        let mut time = 0.0;
        let mut expected = INITIAL_A;
        let first = load(&geometry, 2e6, cx);
        path.advance(time, &first, cx).unwrap();
        for maximum in [8e6, 5e6] {
            time += 1.0;
            let peak = load(&geometry, maximum, cx);
            let sample = path.advance(time, &peak, cx).unwrap();
            close(sample.normal_stress_pa[0], maximum, 1e-6);
            assert_eq!(sample.cracks[0].crack_length_m, expected);
            assert!(sample.completed_cycles.is_none());
            let original = peak
                .solve(
                    &problem,
                    TetThermalStrainField::Uniform(&thermal),
                    Default::default(),
                    FRAME,
                    0.01,
                    cx,
                )
                .unwrap();
            assert_eq!(sample.mechanics.stress, original.stress);
            assert_eq!(sample.mechanics.displacement, original.displacement);
            time += 1.0;
            let returned = path.advance_and_close_cycle(time, &first, cx).unwrap();
            expected = analytical_length(expected, maximum - 2e6);
            close(returned.cracks[0].crack_length_m, expected, 2e-16);
            assert_eq!(
                returned.completed_cycles.unwrap()[0].after,
                path.states()[0]
            );
            // Keep exact accepted arithmetic for the next native state check.
            expected = path.states()[0].crack_length_m;
        }
        assert_eq!(path.states()[0].cycles, 2.0);
        assert_eq!(path.states()[0].blocks, 2);
        assert_eq!(path.observed_samples(), 5);
        assert_eq!(path.finish(cx).unwrap(), path.states());
        assert_eq!(path.finish(cx).unwrap(), path.states());
        assert!(matches!(
            path.advance(5.0, &first, cx),
            Err(TetCrackGrowthError::Finished)
        ));
    });
}

#[test]
fn g0_nested_reversals_and_unclosed_windows_cannot_be_pooled_into_one_cycle() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut path = path(&problem, &thermal, vec![gauge()], cx);
        // Plateaus and intermediate monotone samples are retained without a trace.
        for (i, stress) in [2e6, 2e6, 5e6, 8e6, 4e6].into_iter().enumerate() {
            path.advance(i as f64, &load(&geometry, stress, cx), cx)
                .unwrap();
        }
        let before = path.states();
        assert!(matches!(
            path.close_cycle(cx),
            Err(TetCrackGrowthError::Cycle { .. })
        ));
        assert!(matches!(
            path.finish(cx),
            Err(TetCrackGrowthError::Cycle { .. })
        ));
        assert!(matches!(
            path.advance(5.0, &load(&geometry, 6e6, cx), cx),
            Err(TetCrackGrowthError::Cycle { .. })
        ));
        assert_eq!(path.last_time_s(), Some(4.0));
        assert_eq!(path.observed_samples(), 5);
        assert_eq!(path.states(), before);
        path.advance_and_close_cycle(5.0, &load(&geometry, 2e6, cx), cx)
            .unwrap();
        close(
            path.states()[0].crack_length_m,
            analytical_length(INITIAL_A, 6e6),
            2e-16,
        );
        assert!(matches!(
            path.close_cycle(cx),
            Err(TetCrackGrowthError::Cycle { .. })
        ));
    });
}

#[test]
fn g0_late_gauge_refusal_restores_all_growth_extrema_and_native_time() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut narrow = gauge();
        narrow.law = law(1e6, 1e9);
        let mut path = path(&problem, &thermal, vec![gauge(), narrow], cx);
        path.advance(0.0, &load(&geometry, 2e6, cx), cx).unwrap();
        path.advance(1.0, &load(&geometry, 8e6, cx), cx).unwrap();
        let before = path.states();
        assert!(matches!(
            path.advance_and_close_cycle(2.0, &load(&geometry, 2e6, cx), cx),
            Err(TetCrackGrowthError::Gauge {
                gauge: 1,
                error: CrackGrowthError::DeltaKOutsideDomain { .. }
            })
        ));
        assert_eq!(path.states(), before);
        assert_eq!(path.last_time_s(), Some(1.0));
        assert_eq!(path.observed_samples(), 2);
        // The refused return did not create a reversal: a higher peak can still
        // be proposed at the same time and the original single cycle completed.
        path.advance(2.0, &load(&geometry, 22e6, cx), cx).unwrap();
        path.advance_and_close_cycle(3.0, &load(&geometry, 2e6, cx), cx)
            .unwrap();
        for state in path.states() {
            close(
                state.crack_length_m,
                analytical_length(INITIAL_A, 20e6),
                2e-16,
            );
            assert_eq!(state.cycles, 1.0);
        }
    });
}

#[test]
fn g0_material_temperature_mode_mixity_and_sample_intensity_are_enforced() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut foreign = gauge();
        foreign.material_state_identity = [200; 32];
        let mut hot = gauge();
        let mut scope = hot.law.scope().clone();
        scope.temperature_k = 320.0;
        hot.law = Arc::new(ParisLaw::new(scope, hot.law.parameters(), hot.law.domain()).unwrap());
        for bad in [foreign, hot] {
            assert!(matches!(
                TetCrackGrowthPath::new(
                    &problem,
                    TetThermalStrainField::Uniform(&thermal),
                    FRAME,
                    vec![bad],
                    Default::default(),
                    config(),
                    cx
                ),
                Err(TetCrackGrowthError::Binding(_))
            ));
        }
        let mut inclined = gauge();
        inclined.normal_in_mesh_frame = [1.0, 1.0, 0.0];
        let mut mixed = path(&problem, &thermal, vec![inclined], cx);
        assert!(matches!(
            mixed.advance(0.0, &load(&geometry, 2e6, cx), cx),
            Err(TetCrackGrowthError::MixedMode { .. })
        ));
        assert_eq!(mixed.last_time_s(), None);
        let mut limited = gauge();
        limited.law = law(1.0, 2e6);
        let mut limited = path(&problem, &thermal, vec![limited], cx);
        assert!(matches!(
            limited.advance(0.0, &load(&geometry, 40e6, cx), cx),
            Err(TetCrackGrowthError::Gauge {
                error: CrackGrowthError::MaximumKExceeded { .. },
                ..
            })
        ));
        assert_eq!(limited.observed_samples(), 0);
        let mut tensile = path(&problem, &thermal, vec![gauge()], cx);
        assert!(matches!(
            tensile.advance(0.0, &load(&geometry, -2e6, cx), cx),
            Err(TetCrackGrowthError::Gauge {
                error: CrackGrowthError::StressOutsideDomain { .. },
                ..
            })
        ));
        assert_eq!(tensile.observed_samples(), 0);
    });
}

#[test]
fn g0_a_static_hold_is_not_growth_and_below_fit_excursions_are_not_arrest() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut narrow = gauge();
        narrow.law = law(1e6, 1e9);
        let mut path = path(&problem, &thermal, vec![narrow], cx);
        let low = load(&geometry, 2e6, cx);
        path.advance(0.0, &low, cx).unwrap();
        path.advance_and_close_cycle(1.0, &low, cx).unwrap();
        assert_eq!(path.states()[0].crack_length_m, INITIAL_A);
        path.advance(2.0, &load(&geometry, 3e6, cx), cx).unwrap();
        let before = path.states();
        assert!(matches!(
            path.advance_and_close_cycle(3.0, &low, cx),
            Err(TetCrackGrowthError::Gauge {
                error: CrackGrowthError::DeltaKOutsideDomain { .. },
                ..
            })
        ));
        assert_eq!(path.states(), before);
        assert_eq!(path.last_time_s(), Some(2.0));
    });
}

#[test]
fn g4_cancellation_and_global_sample_budget_do_not_publish_partial_state() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let low = load(&geometry, 2e6, cx);
        let high = load(&geometry, 8e6, cx);
        let mut config = config();
        config.maximum_samples = 2;
        let mut path = TetCrackGrowthPath::new(
            &problem,
            TetThermalStrainField::Uniform(&thermal),
            FRAME,
            vec![gauge()],
            Default::default(),
            config,
            cx,
        )
        .unwrap();
        path.advance(0.0, &low, cx).unwrap();
        let before = path.states();
        with_cx(true, |cancelled| {
            assert!(matches!(
                path.advance(1.0, &high, cancelled),
                Err(TetCrackGrowthError::Cancelled)
            ));
        });
        assert_eq!(path.states(), before);
        assert_eq!(path.last_time_s(), Some(0.0));
        path.advance(1.0, &high, cx).unwrap();
        assert!(matches!(
            path.advance_and_close_cycle(2.0, &low, cx),
            Err(TetCrackGrowthError::Budget)
        ));
        assert_eq!(path.states(), before);
        assert_eq!(path.observed_samples(), 2);
    });
}

#[test]
fn g1_inertial_vibration_drives_the_same_growth_owner_and_exact_checkpoint_replay() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let forcing = load(&geometry, 5e6, cx);
        let mut initial = ZERO;
        initial[1][0] = 2e6 / 200e9;
        let native = forcing
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut original = native.clone();
        let mut observed = native
            .with_crack_growth(FRAME, vec![gauge()], config(), cx)
            .unwrap();
        let mut checkpoint = observed.clone();
        // Independent oscillator: K=E/6, M=rho/60. This Newmark step rotates
        // by 2*pi/32 exactly up to floating-point solve error.
        let omega = (10.0_f64 * 200e9 / 7800.0).sqrt();
        let count = 32;
        let h = 2.0 * (std::f64::consts::PI / f64::from(count)).tan() / omega;
        let mut minimum: f64 = 2e6;
        let mut maximum: f64 = 2e6;
        for step in 1..=count {
            let expected_mechanics = original.advance(&forcing, h, cx).unwrap();
            let actual = if step == count {
                observed.advance_and_close_cycle(&forcing, h, cx).unwrap()
            } else {
                observed.advance(&forcing, h, cx).unwrap()
            };
            assert_eq!(actual.mechanics, expected_mechanics);
            assert_eq!(
                actual.response,
                original.recover_response(FRAME, 0.01, cx).unwrap()
            );
            assert_eq!(observed.mechanics().state(), original.state());
            minimum = minimum.min(actual.normal_stress_pa[0]);
            maximum = maximum.max(actual.normal_stress_pa[0]);
            if step != count {
                assert_eq!(actual.cracks[0].crack_length_m, INITIAL_A);
            }
        }
        close(minimum, 2e6, 0.01);
        close(maximum, 8e6, 0.01);
        close(
            observed.states()[0].crack_length_m,
            analytical_length(INITIAL_A, maximum - minimum),
            2e-16,
        );
        assert_eq!(observed.observed_samples(), 33);
        for step in 1..=count {
            if step == count {
                checkpoint.advance_and_close_cycle(&forcing, h, cx).unwrap();
            } else {
                checkpoint.advance(&forcing, h, cx).unwrap();
            }
        }
        assert_eq!(checkpoint.states(), observed.states());
        assert_eq!(checkpoint.mechanics().state(), observed.mechanics().state());
        assert_eq!(observed.finish(cx).unwrap(), observed.states());
    });
}

#[test]
fn g4_late_inertial_cycle_refusal_keeps_native_motion_and_all_cracks_at_checkpoint() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let forcing = load(&geometry, 5e6, cx);
        let mut initial = ZERO;
        initial[1][0] = 2e6 / 200e9;
        let native = forcing
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut narrow = gauge();
        narrow.law = law(1e6, 1e9);
        let mut observed = native
            .with_crack_growth(FRAME, vec![gauge(), narrow], config(), cx)
            .unwrap();
        let h = 2.0 * (std::f64::consts::PI / 32.0).tan() / (10.0_f64 * 200e9 / 7800.0).sqrt();
        for _ in 1..32 {
            observed.advance(&forcing, h, cx).unwrap();
        }
        let before = observed.clone();
        assert!(matches!(
            observed.advance_and_close_cycle(&forcing, h, cx),
            Err(TetCrackGrowthError::Gauge {
                gauge: 1,
                error: CrackGrowthError::DeltaKOutsideDomain { .. }
            })
        ));
        assert_eq!(observed.states(), before.states());
        assert_eq!(observed.mechanics().state(), before.mechanics().state());
        assert_eq!(observed.observed_samples(), before.observed_samples());
        assert_eq!(
            observed.recover_response(cx).unwrap(),
            before.recover_response(cx).unwrap()
        );
    });
}

#[test]
fn g1_peak_started_cycles_retain_actual_closure_drift_without_snapping() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut path = path(&problem, &thermal, vec![gauge()], cx);
        path.advance(0.0, &load(&geometry, 8e6, cx), cx).unwrap();
        path.advance(1.0, &load(&geometry, 2e6, cx), cx).unwrap();
        let before = path.states();
        assert!(matches!(
            path.advance_and_close_cycle(2.0, &load(&geometry, 8e6 + 1.0, cx), cx),
            Err(TetCrackGrowthError::Cycle { .. })
        ));
        assert_eq!(path.states(), before);
        assert_eq!(path.last_time_s(), Some(1.0));
        // This finite mismatch is within the declared 0.00801 Pa allowance.
        // The actual endpoint becomes the next anchor; it is never snapped.
        let returned_peak = 8e6 - 0.001;
        let first = path
            .advance_and_close_cycle(2.0, &load(&geometry, returned_peak, cx), cx)
            .unwrap();
        close(first.normal_stress_pa[0], returned_peak, 1e-6);
        let a1 = analytical_length(INITIAL_A, 6e6);
        close(first.cracks[0].crack_length_m, a1, 2e-16);
        path.advance(3.0, &load(&geometry, 2e6, cx), cx).unwrap();
        let second = path
            .advance_and_close_cycle(4.0, &load(&geometry, returned_peak, cx), cx)
            .unwrap();
        close(
            second.cracks[0].crack_length_m,
            analytical_length(a1, returned_peak - 2e6),
            2e-16,
        );
        assert_eq!(path.states()[0].cycles, 2.0);
        assert_eq!(path.observed_samples(), 5);
    });
}

#[test]
fn g1_moving_support_cycles_preserve_original_inertia_work_and_stress() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut initial = ZERO;
        initial[1][0] = 2e6 / 200e9;
        let mut original = load
            .begin_transient_with_supports(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                0.0,
                &initial,
                &ZERO,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut observed = original
            .clone()
            .with_crack_growth(FRAME, vec![gauge()], config(), cx)
            .unwrap();
        // A run with moving supports cannot silently switch to stationary support stepping.
        assert!(observed.advance(&load, 0.05, cx).is_err());
        assert_eq!(observed.observed_samples(), 1);
        for (index, stress) in [8e6, 2e6].into_iter().enumerate() {
            let supports = [(3, stress / 200e9)];
            let expected = original
                .advance_with_supports(&load, &supports, 0.05, cx)
                .unwrap();
            let actual = if index == 1 {
                observed
                    .advance_with_supports_and_close_cycle(&load, &supports, 0.05, cx)
                    .unwrap()
            } else {
                observed
                    .advance_with_supports(&load, &supports, 0.05, cx)
                    .unwrap()
            };
            assert_eq!(actual.mechanics, expected);
            assert_eq!(
                actual.response,
                original.recover_response(FRAME, 0.01, cx).unwrap()
            );
            assert_eq!(observed.mechanics().state(), original.state());
            close(actual.normal_stress_pa[0], stress, 1e-6);
        }
        close(
            observed.states()[0].crack_length_m,
            analytical_length(INITIAL_A, 6e6),
            2e-16,
        );
        assert_eq!(observed.states()[0].cycles, 1.0);
        assert_eq!(observed.observed_samples(), 3);
    });
}

#[test]
fn g1_declared_repeated_native_cycle_matches_analytical_growth_with_three_samples() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut path = path(&problem, &thermal, vec![gauge()], cx);
        assert!(matches!(
            path.cycles_to_length(0, 2.0 * INITIAL_A, cx),
            Err(TetCrackGrowthError::Cycle { gauge: 0, .. })
        ));
        for (index, stress) in [2e6, 8e6, 2e6].into_iter().enumerate() {
            path.advance(index as f64, &load(&geometry, stress, cx), cx)
                .unwrap();
        }
        let before = path.states();
        assert_eq!(before[0].crack_length_m, INITIAL_A);
        assert_eq!(before[0].cycles, 0.0);
        let forecast = path.cycles_to_length(0, 2.0 * INITIAL_A, cx).unwrap();
        let expected_exposure = std::f64::consts::LN_2 / (1e-16 * std::f64::consts::PI * 6e6 * 6e6);
        close(forecast.additional_cycles, expected_exposure, 2e-10);
        assert_eq!(forecast.before, before[0]);
        assert_eq!(forecast.target_length_m, 2.0 * INITIAL_A);
        assert_eq!(path.states(), before);
        assert_eq!(path.observed_samples(), 3);
        assert_eq!(path.last_time_s(), Some(2.0));
        let steps = path.close_repeated_cycle(100, cx).unwrap();
        // Independently integrated m=2 law over the declared 100 repetitions:
        // a(N) = a(0) exp((1e-4/1e12) pi (delta_sigma)^2 N).
        let expected = INITIAL_A * (1e-16 * std::f64::consts::PI * 6e6 * 6e6 * 100.0).exp();
        let after = path.states();
        close(after[0].crack_length_m, expected, 2e-16);
        assert_eq!(after[0].cycles, 100.0);
        assert_eq!(after[0].blocks, 1);
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].before, before[0]);
        assert_eq!(steps[0].after, after[0]);
        assert_eq!(path.observed_samples(), 3);
        assert_eq!(path.last_time_s(), Some(2.0));
        assert!(matches!(
            path.cycles_to_length(0, 4.0 * INITIAL_A, cx),
            Err(TetCrackGrowthError::Cycle { gauge: 0, .. })
        ));
        assert_eq!(path.finish(cx).unwrap(), after);
    });
}

#[test]
fn g4_repeated_cycle_late_gauge_budget_refusal_keeps_every_window_reusable() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut limited = gauge();
        limited.limits.maximum_cycles = 10.0;
        let mut path = path(&problem, &thermal, vec![gauge(), limited], cx);
        for (index, stress) in [2e6, 8e6, 2e6].into_iter().enumerate() {
            path.advance(index as f64, &load(&geometry, stress, cx), cx)
                .unwrap();
        }
        let before = path.states();
        assert!(matches!(
            path.close_repeated_cycle(100, cx),
            Err(TetCrackGrowthError::Gauge {
                gauge: 1,
                error: CrackGrowthError::Budget,
            })
        ));
        assert_eq!(path.states(), before);
        assert_eq!(path.observed_samples(), 3);
        assert_eq!(path.last_time_s(), Some(2.0));
        // Retry directly: neither the first gauge's growth nor either window
        // may have been consumed by the second gauge's late refusal.
        let steps = path.close_repeated_cycle(10, cx).unwrap();
        let expected = INITIAL_A * (1e-16 * std::f64::consts::PI * 6e6 * 6e6 * 10.0).exp();
        for (index, state) in path.states().into_iter().enumerate() {
            close(state.crack_length_m, expected, 2e-16);
            assert_eq!(state.cycles, 10.0);
            assert_eq!(state.blocks, 1);
            assert_eq!(steps[index].before, before[index]);
            assert_eq!(steps[index].after, state);
        }
        assert_eq!(path.observed_samples(), 3);
        assert_eq!(path.last_time_s(), Some(2.0));
    });
}

#[test]
fn g0_invalid_repetition_counts_refuse_atomically_without_consuming_the_cycle() {
    with_cx(false, |cx| {
        let material = material();
        let thermal = thermal();
        let fixed = fixed();
        let problem = problem(&material, &fixed);
        let geometry = geometry(cx);
        let mut wide_budget = gauge();
        wide_budget.limits.maximum_cycles = (1_u64 << 53) as f64;
        let mut path = path(&problem, &thermal, vec![wide_budget], cx);
        for (index, stress) in [2e6, 8e6, 2e6].into_iter().enumerate() {
            path.advance(index as f64, &load(&geometry, stress, cx), cx)
                .unwrap();
        }
        let before = path.states();
        // The first unrepresentable integer must not be rounded to 2^53.
        // Invalid input must be distinguished from later physical-domain or
        // exposure-budget refusals of an otherwise admitted count.
        for count in [0, (1_u64 << 53) + 1] {
            assert!(matches!(
                path.close_repeated_cycle(count, cx),
                Err(TetCrackGrowthError::Invalid(_))
            ));
            assert_eq!(path.states(), before);
            assert_eq!(path.observed_samples(), 3);
            assert_eq!(path.last_time_s(), Some(2.0));
        }
        let steps = path.close_repeated_cycle(1, cx).unwrap();
        let after = path.states();
        close(
            after[0].crack_length_m,
            analytical_length(INITIAL_A, 6e6),
            2e-16,
        );
        assert_eq!(after[0].cycles, 1.0);
        assert_eq!(after[0].blocks, 1);
        assert_eq!(steps[0].before, before[0]);
        assert_eq!(steps[0].after, after[0]);
        assert_eq!(path.observed_samples(), 3);
        assert_eq!(path.last_time_s(), Some(2.0));
    });
}

struct ParisSourceFixture {
    card: MaterialCard,
    fit_point: QueryPoint,
    bundle: Arc<ResolvedMaterialStatePoint>,
    material: TetElasticMaterial,
    thermal: TetThermalStrainState,
    law: Arc<ParisLaw>,
}

fn paris_source_request() -> ParisSourceRequest {
    ParisSourceRequest {
        temperature_axis: "absolute_temperature".into(),
        stress_ratio_axis: "calibration_load_ratio".into(),
        delta_k_squared_axis: "calibration_delta_k_squared".into(),
        surface_and_geometry: "synthetic constant-Y nominal-stress coupon".into(),
        environment: "synthetic fixed humidity context".into(),
        load_basis: "tensile full-range Mode-I Paris fit; metres per full cycle".into(),
    }
}

fn paris_source_fixture(revision: u32) -> ParisSourceFixture {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let request = paris_source_request();
    let external_domain = ValidityDomain::unconstrained()
        .with_quantity(&request.temperature_axis, temperature, 290.0, 330.0)
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
    let native_point = QueryPoint::new()
        .with_quantity(&request.temperature_axis, temperature, 300.0)
        .unwrap()
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.2)
        .unwrap();
    let fit_domain = external_domain
        .clone()
        .with_quantity(
            &request.stress_ratio_axis,
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            0.95,
        )
        .with_quantity(
            &request.delta_k_squared_axis,
            QuantitySpec::dimensional(STRESS_INTENSITY_SQUARED_DIMS),
            1.0,
            1e16,
        );
    let fit_point = native_point
        .clone()
        .with_quantity(
            &request.stress_ratio_axis,
            QuantitySpec::dimensional(Dims::NONE),
            0.25,
        )
        .unwrap()
        .with_quantity(
            &request.delta_k_squared_axis,
            QuantitySpec::dimensional(STRESS_INTENSITY_SQUARED_DIMS),
            1e12,
        )
        .unwrap();
    let mut claims = ClaimSet::new();
    for (key, value, validity) in [
        (
            PropertyKey::new("density", Density::DIMS),
            7800.0,
            external_domain.clone(),
        ),
        (
            PropertyKey::new("young_modulus", Pressure::DIMS),
            200e9,
            external_domain.clone(),
        ),
        (
            PropertyKey::new("poisson_ratio", Dims::NONE),
            0.0,
            external_domain,
        ),
        (
            ParisProperty::ReferenceDeltaKSquared.property_key(),
            1e12,
            fit_domain.clone(),
        ),
        (
            ParisProperty::ReferenceGrowthRate.property_key(),
            1e-4,
            fit_domain.clone(),
        ),
        (
            ParisProperty::Exponent.property_key(),
            2.0,
            fit_domain.clone(),
        ),
        (
            ParisProperty::MaximumKSquared.property_key(),
            1e18,
            fit_domain,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                value: PropertyValue::Scalar {
                    value,
                    dims: key.dims(),
                },
                key,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic native Paris source fixture; not measured".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic native Paris material".into(),
            phase: "solid".into(),
            process: "synthetic coupon process".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    }
    let bundle = Arc::new(
        SolidStiffnessModel::Isotropic
            .resolve_reference(
                &card,
                &native_point,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
    );
    let material = SolidStiffnessModel::Isotropic
        .material(
            &bundle,
            bundle.property("density").unwrap().value_si(),
            bundle.identity(),
            0.01,
        )
        .unwrap();
    let thermal = TetThermalStrainState::try_new(
        300.0,
        300.0,
        [0.0; 6],
        bundle.identity().0,
        ContentHash([106; 32]),
        0.01,
    )
    .unwrap();
    let law = Arc::new(
        resolve_paris(
            &card,
            &fit_point,
            &request,
            [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 4],
        )
        .unwrap(),
    );
    ParisSourceFixture {
        card,
        fit_point,
        bundle,
        material,
        thermal,
        law,
    }
}

#[test]
fn g0_sourced_native_paris_retains_original_card_and_refuses_foreign_operating_context() {
    with_cx(false, |cx| {
        let source = paris_source_fixture(0);
        let foreign = paris_source_fixture(1);
        let fixed = fixed();
        let problem = problem(&source.material, &fixed);
        let geometry = geometry(cx);
        let mut bound = gauge();
        bound.material_state_identity = source.bundle.identity().0;
        bound.material_source = Some(source.bundle.clone());
        bound.law = source.law.clone();
        assert_eq!(source.bundle.state_coordinate("T"), None);
        // Native material queries contain external conditions, not calibration
        // R/Delta-K coordinates. The original custom temperature axis is kept.
        assert_eq!(
            source.bundle.state_coordinate("absolute_temperature"),
            Some(300.0)
        );
        assert_eq!(
            source.bundle.state_coordinate("calibration_load_ratio"),
            None
        );
        let mut static_path = path(&problem, &source.thermal, vec![bound.clone()], cx);
        assert!(std::ptr::eq(
            static_path.material_source(0).unwrap(),
            source.bundle.as_ref()
        ));
        let retained = static_path.law(0).unwrap().source().unwrap();
        assert_eq!(retained.card_identity(), source.card.content_hash());
        assert_eq!(retained.query(), &source.fit_point);
        assert_eq!(
            retained
                .property(ParisProperty::ReferenceDeltaKSquared)
                .evidence
                .value
                .value,
            1e12
        );
        for (index, stress) in [2e6, 8e6, 2e6].into_iter().enumerate() {
            static_path
                .advance(index as f64, &load(&geometry, stress, cx), cx)
                .unwrap();
        }
        static_path.close_cycle(cx).unwrap();
        close(
            static_path.states()[0].crack_length_m,
            analytical_length(INITIAL_A, 6e6),
            2e-16,
        );

        let forcing = load(&geometry, 5e6, cx);
        let mut initial = ZERO;
        initial[1][0] = 2e6 / 200e9;
        let native = forcing
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&source.thermal)),
                0.0,
                &initial,
                &ZERO,
                TransientConfig3::default(),
                cx,
            )
            .unwrap();
        let mut observed = native
            .clone()
            .with_crack_growth(FRAME, vec![bound.clone()], config(), cx)
            .unwrap();
        assert!(std::ptr::eq(
            observed.material_source(0).unwrap(),
            source.bundle.as_ref()
        ));
        assert_eq!(
            observed.law(0).unwrap().source().unwrap().card_identity(),
            source.card.content_hash()
        );
        let h = 2.0 * (std::f64::consts::PI / 32.0).tan() / (10.0_f64 * 200e9 / 7800.0).sqrt();
        for step in 1..=32 {
            if step == 32 {
                observed.advance_and_close_cycle(&forcing, h, cx).unwrap();
            } else {
                observed.advance(&forcing, h, cx).unwrap();
            }
        }
        close(
            observed.states()[0].crack_length_m,
            static_path.states()[0].crack_length_m,
            2e-16,
        );
        assert_eq!(
            observed.law(0).unwrap().source().unwrap().query(),
            &source.fit_point
        );

        let mut missing_bundle = bound.clone();
        missing_bundle.material_source = None;
        let mut foreign_bundle = bound.clone();
        foreign_bundle.material_source = Some(foreign.bundle);
        let mut foreign_fit = bound.clone();
        foreign_fit.law = foreign.law;
        let other_context = source
            .fit_point
            .clone()
            .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.8)
            .unwrap();
        let mut foreign_context = bound;
        foreign_context.law = Arc::new(
            resolve_paris(
                &source.card,
                &other_context,
                &paris_source_request(),
                [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 4],
            )
            .unwrap(),
        );
        for bad in [missing_bundle, foreign_bundle, foreign_fit, foreign_context] {
            assert!(matches!(
                TetCrackGrowthPath::new(
                    &problem,
                    TetThermalStrainField::Uniform(&source.thermal),
                    FRAME,
                    vec![bad.clone()],
                    Default::default(),
                    config(),
                    cx,
                ),
                Err(TetCrackGrowthError::Binding(
                    TetFatigueError::SourceMismatch { .. }
                ))
            ));
            assert!(matches!(
                native
                    .clone()
                    .with_crack_growth(FRAME, vec![bad], config(), cx),
                Err(TetCrackGrowthError::Binding(
                    TetFatigueError::SourceMismatch { .. }
                ))
            ));
        }
    });
}
