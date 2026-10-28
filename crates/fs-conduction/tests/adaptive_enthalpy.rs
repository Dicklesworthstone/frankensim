//! G1/G3/G4: adaptive spatial heat transport against independent modal and
//! piecewise-linear latent-heat solutions, including rejected work and resume.

mod support;

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::adaptive::{
    AdaptiveEnthalpyConfig, AdaptiveEnthalpyStatus, EnthalpyAcceptedTrial,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ScalarField, ThermalBoundary,
    ThermalBoundaryBuilder,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;

fn mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    .unwrap()
}

fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    // h=5(T-280), rho0=2: volumetric heat capacity = 10 J/(m³ K).
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x64; 32]),
        SolidLiquidPhase::Solid,
        [(0.0, 280.0), (1100.0, 500.0)]
            .into_iter()
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2.0,
            })
            .collect(),
    )
    .unwrap()
}

fn latent() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x65; 32]),
        [
            (0.0, 300.0, 0.0),
            (100.0, 350.0, 0.0),
            (300.0, 350.0, 1.0),
            (600.0, 500.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 2.0,
        })
        .collect(),
    )
    .unwrap()
}

fn config(temperature_tolerance: f64) -> AdaptiveEnthalpyConfig {
    AdaptiveEnthalpyConfig {
        enthalpy_absolute_tolerance_j_kg: 5.0 * temperature_tolerance,
        temperature_absolute_tolerance_k: temperature_tolerance,
        liquid_fraction_absolute_tolerance: 1e-4,
        relative_tolerance: 0.0,
        minimum_step_s: 1e-8,
        maximum_step_s: 1.0,
        endpoint: EnthalpyStepConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-12,
                relative_tolerance: 1e-12,
                linear_restart: 4,
                max_linear_cycles: 8,
                forcing_minimum: 1e-12,
                forcing_maximum: 1e-3,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 32,
            energy_tolerance_j: 1e-9,
        },
    }
}

fn prepare<'m, 'c>(
    cx: &fs_exec::Cx<'_>,
    mesh: &'m ConductionMesh,
    curve: &'c EquilibriumEnthalpyPhaseCurve,
) -> EnthalpyBackwardEuler<'m, 'c> {
    EnthalpyBackwardEuler::uniform(
        cx,
        mesh,
        curve,
        2.0,
        EnthalpyBudget {
            max_vertices: 4,
            max_elements: 1,
        },
    )
    .unwrap()
}

fn problem<'a>(
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    material: &'a ConductivityModel,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem {
        mesh,
        boundary,
        material,
        source,
        element_materials: None,
    }
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn observation_bits(trial: EnthalpyAcceptedTrial<'_>) -> Vec<u64> {
    let mut result = vec![
        trial.accepted_trial as u64,
        trial.start_time_s.to_bits(),
        trial.midpoint_time_s.to_bits(),
        trial.end_time_s.to_bits(),
        trial.error.normalized.to_bits(),
    ];
    for half in [trial.first_half, trial.second_half] {
        result.extend(bits(&half.specific_enthalpy_j_kg));
        result.extend(bits(&half.temperature));
        result.extend(bits(&half.liquid_mass_fraction));
        result.push(half.stored_energy_change_j.to_bits());
    }
    result
}

#[test]
fn tighter_time_tolerance_resolves_an_independent_spatial_diffusion_mode() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let stepper = prepare(cx, &mesh, &curve);
        let p = problem(&mesh, &boundary, &material, &source);
        let initial = [200.0, 250.0, 150.0, 200.0];
        let mut errors = Vec::new();
        let mut steps = Vec::new();
        for tolerance in [1e-2, 1e-4] {
            let mut state = stepper
                .adaptive_state(cx, &initial, 0.0, 1.0, config(tolerance))
                .unwrap();
            let report = state.advance_to(cx, p, None, 1.6, 10000).unwrap();
            assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
            assert!(state.rejected_steps() > 0);
            assert_eq!(
                report.endpoint_solves,
                3 * (report.accepted + report.rejected)
            );
            // On this exact tet, [0,1,-1,0] is an eigenvector of the
            // independently assembled K/M with eigenvalue 4k/(rho Cp)=1.2/s.
            let amplitude = 10.0 * (-1.2_f64 * 1.6).exp();
            let reference = [320.0, 320.0 + amplitude, 320.0 - amplitude, 320.0];
            let error = state
                .temperature()
                .iter()
                .zip(reference)
                .map(|(got, expected)| (got - expected).abs())
                .fold(0.0_f64, f64::max);
            errors.push(error);
            steps.push(state.accepted_steps());
            assert!(report.integrated_net_input_j.abs() < 1e-14);
            assert!(report.stored_energy_change_j.abs() < 1e-8);
            assert!(report.energy_residual_j.abs() < 1e-8);
        }
        assert!(errors[1] < 0.02, "independent modal errors: {errors:?}");
        assert!(
            errors[1] < 0.2 * errors[0],
            "independent modal errors: {errors:?}"
        );
        assert!(steps[1] > steps[0]);
    });
}

#[test]
fn liquid_fraction_controls_time_error_while_latent_nodes_stay_isothermal() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = latent();
        let stepper = prepare(cx, &mesh, &curve);
        let mut policy = config(100.0);
        policy.enthalpy_absolute_tolerance_j_kg = 1000.0;
        policy.liquid_fraction_absolute_tolerance = 1e-6;
        let mut state = stepper
            .adaptive_state(cx, &[80.0, 200.0, 200.0, 320.0], 0.0, 0.15, policy)
            .unwrap();
        let report = state
            .advance_to(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                0.15,
                10000,
            )
            .unwrap();
        assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert!(state.rejected_steps() > 0 && state.accepted_steps() > 1);
        assert_eq!(&state.temperature()[1..3], &[350.0, 350.0]);
        // In these unchanged chart segments, y=(T0-350,T3-350) obeys
        // y'=[[-3,1],[1,-1]]y with y(0)=(-10,10). The exact matrix
        // exponential gives h1=h2=200-y0-y3 from conserved enthalpy.
        let t = 0.15_f64;
        let a = 2.0_f64.sqrt() * t;
        let y0 = (-2.0 * t).exp() * (-10.0 * a.cosh() + 10.0 * 2.0_f64.sqrt() * a.sinh());
        let y3 = (-2.0 * t).exp() * (10.0 * a.cosh());
        let reference = [
            100.0 + 2.0 * y0,
            200.0 - y0 - y3,
            200.0 - y0 - y3,
            300.0 + 2.0 * y3,
        ];
        for (got, expected) in state.specific_enthalpy_j_kg().iter().zip(reference) {
            assert!(
                (got - expected).abs() < 0.02,
                "latent modal h: {got} versus {expected}"
            );
        }
        assert!(report.energy_residual_j.abs() < 1e-8);
        assert!(report.last_error.unwrap().normalized <= 1.0);
    });
}

#[test]
fn events_and_accepted_half_step_heat_cross_a_true_latent_plateau() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let heating = ScalarField::Uniform(100.0);
        let cooling = ScalarField::Uniform(-50.0);
        let curve = latent();
        let stepper = prepare(cx, &mesh, &curve);
        let mut state = stepper
            .adaptive_state(cx, &[50.0; 4], 0.0, 0.4, config(1e-5))
            .unwrap();
        let hot = state
            .advance_to(
                cx,
                problem(&mesh, &boundary, &material, &heating),
                None,
                5.0,
                1000,
            )
            .unwrap();
        assert_eq!(hot.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(state.time_s().to_bits(), 5.0_f64.to_bits());
        for &h in state.specific_enthalpy_j_kg() {
            assert!((h - 300.0).abs() < 1e-8);
        }
        assert!((hot.integrated_net_input_j - 500.0 / 6.0).abs() < 1e-10);
        // A discontinuous workload is changed ONLY after reaching its event.
        let cold = state
            .advance_to(
                cx,
                problem(&mesh, &boundary, &material, &cooling),
                None,
                7.0,
                1000,
            )
            .unwrap();
        assert_eq!(cold.status, AdaptiveEnthalpyStatus::ReachedEnd);
        for &h in state.specific_enthalpy_j_kg() {
            assert!((h - 250.0).abs() < 1e-8);
        }
        assert_eq!(state.temperature(), &[350.0; 4]);
        for &fraction in state.liquid_mass_fraction() {
            assert!((fraction - 0.75).abs() < 1e-10);
        }
        assert!((cold.integrated_net_input_j - 400.0 / 6.0).abs() < 1e-10);
        assert!((cold.stored_energy_change_j - 400.0 / 6.0).abs() < 1e-8);
        assert!(cold.energy_residual_j.abs() < 1e-8);
    });
}

#[test]
fn adaptive_pointwise_radiation_resolves_independent_quartic_cooling() {
    use fs_conduction::{
        AmbientRadiationPatch, EMISSIVITY_DIMS, PointwiseAmbientRadiation,
        STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ThermalBc,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
    };
    support::with_cx(|cx| {
        // Equal surface areas and nodal masses preserve the uniform mode on a
        // regular tet. Its exact continuous equation is C T'=-epsilon sigma A
        // (T^4-Ta^4), independent of the production quadrature and time solver.
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
            vec![
                [0.0; 3],
                [1.0, 0.0, 0.0],
                [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
                [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
            ],
        )
        .unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("exterior", ThermalBc::adiabatic())
            .unwrap()
            .finish()
            .unwrap();
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
                value: PropertyValue::Scalar {
                    value: 0.8,
                    dims: EMISSIVITY_DIMS,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 500.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic adaptive radiation reference".into(),
                    license: "internal-test-use".into(),
                    artifact: None,
                },
            })
            .unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic gray surface".into(),
                phase: "solid".into(),
                process: "numerical reference".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap();
        let emissivity = SurfaceEmissivity::from_card(
            "exterior",
            &card,
            400.0,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        let patches = [AmbientRadiationPatch::new("exterior", emissivity, 300.0).unwrap()];
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let stepper = prepare(cx, &mesh, &curve);
        let mut policy = config(1e-4);
        policy.endpoint.newton.absolute_tolerance = 1e-10;
        policy.endpoint.energy_tolerance_j = 1e-8;
        let mut state = stepper
            .adaptive_state(cx, &[600.0; 4], 0.0, 0.02, policy)
            .unwrap();
        let mut observed_pairs = 0;
        let mut observed_input_j = 0.0;
        let report = state
            .advance_to_with_pointwise_ambient_radiation_observed(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                &radiation,
                0.02,
                10000,
                |trial| {
                    observed_pairs += 1;
                    assert_eq!(trial.accepted_trial, observed_pairs);
                    for (half, dt) in [
                        (trial.first_half, trial.midpoint_time_s - trial.start_time_s),
                        (trial.second_half, trial.end_time_s - trial.midpoint_time_s),
                    ] {
                        let radiation = half.pointwise_radiation.as_ref().unwrap();
                        observed_input_j -= dt * radiation.outward_heat_w;
                    }
                    ControlFlow::Continue(())
                },
            )
            .unwrap();
        assert_eq!(observed_pairs, report.accepted);
        assert!((observed_input_j - report.integrated_net_input_j).abs() < 1e-10);
        assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert!(state.rejected_steps() > 0);
        let capacity = 10.0 * 2.0_f64.sqrt() / 12.0;
        let alpha = 0.8 * STEFAN_BOLTZMANN_W_M2_K4 * 3.0_f64.sqrt() / capacity;
        // This analytic primitive differentiates to 1/(T^4-Ta^4).
        let primitive = |t: f64| {
            ((t - 300.0) / (t + 300.0)).ln() / (4.0 * 300.0_f64.powi(3))
                - (t / 300.0).atan() / (2.0 * 300.0_f64.powi(3))
        };
        let target = primitive(400.0) - alpha * 0.02;
        let (mut low, mut high) = (300.0, 400.0);
        for _ in 0..80 {
            let middle = 0.5 * (low + high);
            if primitive(middle) < target {
                low = middle;
            } else {
                high = middle;
            }
        }
        let exact = 0.5 * (low + high);
        for &t in state.temperature() {
            assert!(
                (t - exact).abs() < 0.04,
                "quartic cooling {t} versus {exact}"
            );
        }
        assert!(report.integrated_net_input_j < 0.0);
        assert!(report.energy_residual_j.abs() < 1e-7);
        let stored = capacity * (state.temperature().iter().sum::<f64>() / 4.0 - 400.0);
        assert!((report.stored_energy_change_j - stored).abs() < 1e-10);
        assert!(state.last_endpoint().unwrap().pointwise_radiation.is_some());
    });
}

#[test]
fn rejected_trials_and_budget_resume_preserve_the_same_accepted_trajectory() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let stepper = prepare(cx, &mesh, &curve);
        let p = problem(&mesh, &boundary, &material, &source);
        let initial = [200.0, 250.0, 150.0, 200.0];
        let mut straight = stepper
            .adaptive_state(cx, &initial, 0.0, 1.0, config(1e-3))
            .unwrap();
        let mut resumed = straight.clone();
        let mut full_trace = Vec::new();
        let full = straight
            .advance_to_observed(cx, p, None, 1.0, 10000, |trial| {
                full_trace.push(observation_bits(trial));
                ControlFlow::Continue(())
            })
            .unwrap();
        assert_eq!(full_trace.len(), straight.accepted_steps());
        let zero = resumed
            .advance_to_observed(cx, p, None, 1.0, 0, |_| {
                panic!("zero budget must not publish a fine pair")
            })
            .unwrap();
        assert_eq!(zero.status, AdaptiveEnthalpyStatus::TrialLimit);
        assert_eq!(zero.endpoint_solves, 0);
        let rejected = resumed
            .advance_to_observed(cx, p, None, 1.0, 1, |_| {
                panic!("rejected comparison must not publish either half step")
            })
            .unwrap();
        assert_eq!(rejected.status, AdaptiveEnthalpyStatus::TrialLimit);
        assert_eq!((rejected.accepted, rejected.rejected), (0, 1));
        assert_eq!(bits(resumed.specific_enthalpy_j_kg()), bits(&initial));
        assert_eq!(rejected.stored_energy_change_j, 0.0);
        assert_eq!(rejected.integrated_net_input_j, 0.0);
        let mut resumed_trace = Vec::new();
        let last = loop {
            let report = resumed
                .advance_to_observed(cx, p, None, 1.0, 1, |trial| {
                    resumed_trace.push(observation_bits(trial));
                    ControlFlow::Continue(())
                })
                .unwrap();
            if report.status == AdaptiveEnthalpyStatus::ReachedEnd {
                break report;
            }
            assert_eq!(report.status, AdaptiveEnthalpyStatus::TrialLimit);
        };
        assert_eq!(full_trace, resumed_trace);
        assert_eq!(
            bits(straight.specific_enthalpy_j_kg()),
            bits(resumed.specific_enthalpy_j_kg())
        );
        assert_eq!(
            straight.next_step_s().to_bits(),
            resumed.next_step_s().to_bits()
        );
        assert_eq!(straight.accepted_steps(), resumed.accepted_steps());
        assert_eq!(straight.rejected_steps(), resumed.rejected_steps());
        assert_eq!(straight.endpoint_solves(), resumed.endpoint_solves());
        assert_eq!(
            full.integrated_net_input_j.to_bits(),
            last.integrated_net_input_j.to_bits()
        );
        assert_eq!(
            full.energy_residual_j.to_bits(),
            last.energy_residual_j.to_bits()
        );
    });
}

#[test]
fn minimum_step_and_cancellation_preserve_accepted_physics() {
    support::with_gate(|gate, cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let stepper = prepare(cx, &mesh, &curve);
        let p = problem(&mesh, &boundary, &material, &source);
        let initial = [200.0, 250.0, 150.0, 200.0];
        let mut policy = config(1e-10);
        policy.minimum_step_s = 0.5;
        policy.maximum_step_s = 0.5;
        let mut state = stepper
            .adaptive_state(cx, &initial, 0.0, 0.5, policy)
            .unwrap();
        let refused = state.advance_to(cx, p, None, 1.0, 100).unwrap();
        assert_eq!(refused.status, AdaptiveEnthalpyStatus::MinimumStep);
        assert_eq!((refused.accepted, refused.rejected), (0, 1));
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(bits(state.specific_enthalpy_j_kg()), bits(&initial));
        assert_eq!(refused.integrated_net_input_j, 0.0);
        let mut progressed = stepper
            .adaptive_state(cx, &initial, 0.0, 0.1, config(1e-3))
            .unwrap();
        let accepted = progressed.advance_to(cx, p, None, 0.2, 1000).unwrap();
        assert_eq!(accepted.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert!(progressed.accepted_steps() > 0);
        let before = bits(progressed.specific_enthalpy_j_kg());
        let before_time = progressed.time_s().to_bits();
        let before_step = progressed.next_step_s().to_bits();
        let before_solves = progressed.endpoint_solves();
        gate.request();
        assert!(matches!(
            progressed.advance_to(cx, p, None, 1.0, 1),
            Err(
                fs_conduction::transient::enthalpy::EnthalpyError::Conduction(
                    fs_conduction::ConductionError::Cancelled { .. }
                )
            )
        ));
        assert_eq!(bits(progressed.specific_enthalpy_j_kg()), before);
        assert_eq!(progressed.time_s().to_bits(), before_time);
        assert_eq!(progressed.next_step_s().to_bits(), before_step);
        assert_eq!(progressed.endpoint_solves(), before_solves);
    });
}

#[test]
fn accepted_pair_exposes_a_local_half_step_peak_without_cloning_the_endpoint() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let stepper = prepare(cx, &mesh, &curve);
        let mut policy = config(10.0);
        policy.maximum_step_s = 4.0;
        // T(0) = (380,310,295,295). Independent star-Laplacian modes give
        // T1(t) = 320 + 10 exp(-1.2t) - 20 exp(-4.8t).
        let mut state = stepper
            .adaptive_state(cx, &[500.0, 150.0, 75.0, 75.0], 0.0, 4.0, policy)
            .unwrap();
        let initial_local_temperature = state.temperature()[1];
        let mut pairs = 0;
        let mut observed_peak = initial_local_temperature;
        let mut final_buffer = std::ptr::null();
        let report = state
            .advance_to_observed(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                4.0,
                1,
                |trial| {
                    pairs += 1;
                    assert_eq!(trial.accepted_trial, 1);
                    assert_eq!(
                        (trial.start_time_s, trial.midpoint_time_s, trial.end_time_s),
                        (0.0, 2.0, 4.0)
                    );
                    // Each two-second backward-Euler half step multiplies the
                    // independent modes by 1/3.4 and 1/10.6 respectively.
                    let first = 320.0 + 10.0 / 3.4 - 20.0 / 10.6;
                    let second = 320.0 + 10.0 / 3.4_f64.powi(2)
                        - 20.0 / 10.6_f64.powi(2);
                    assert!((trial.first_half.temperature[1] - first).abs() < 1e-8);
                    assert!((trial.second_half.temperature[1] - second).abs() < 1e-8);
                    assert!(first > second && first > initial_local_temperature);
                    observed_peak = observed_peak
                        .max(trial.first_half.temperature[1])
                        .max(trial.second_half.temperature[1]);
                    final_buffer = trial.second_half.specific_enthalpy_j_kg.as_ptr();
                    assert!(trial.error.normalized <= 1.0);
                    // Reaching the requested event takes precedence over a stop.
                    ControlFlow::Break(())
                },
            )
            .unwrap();
        assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!((report.accepted, report.rejected, pairs), (1, 0, 1));
        assert!(observed_peak > state.temperature()[1]);
        assert_eq!(final_buffer, state.specific_enthalpy_j_kg().as_ptr());
        assert!(report.energy_residual_j.abs() < 1e-8);
    });
}

#[test]
fn observer_stop_commits_a_pair_and_resume_does_not_repeat_it() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(10.0);
        let curve = sensible();
        let stepper = prepare(cx, &mesh, &curve);
        let p = problem(&mesh, &boundary, &material, &source);
        let mut straight = stepper
            .adaptive_state(cx, &[200.0; 4], 0.0, 0.25, config(1e-4))
            .unwrap();
        let mut resumed = straight.clone();
        let mut full_trace = Vec::new();
        let full = straight
            .advance_to_observed(cx, p, None, 1.0, 100, |trial| {
                full_trace.push(observation_bits(trial));
                ControlFlow::Continue(())
            })
            .unwrap();
        let mut resumed_trace = Vec::new();
        let stopped = resumed
            .advance_to_observed(cx, p, None, 1.0, 100, |trial| {
                resumed_trace.push(observation_bits(trial));
                ControlFlow::Break(())
            })
            .unwrap();
        assert_eq!(stopped.status, AdaptiveEnthalpyStatus::ObserverStop);
        assert_eq!((stopped.accepted, stopped.rejected, stopped.endpoint_solves), (1, 0, 3));
        assert_eq!(resumed.time_s(), 0.25);
        assert_eq!(resumed.accepted_steps(), 1);
        assert!(stopped.integrated_net_input_j > 0.0);
        let last = resumed
            .advance_to_observed(cx, p, None, 1.0, 100, |trial| {
                resumed_trace.push(observation_bits(trial));
                ControlFlow::Continue(())
            })
            .unwrap();
        assert_eq!(last.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(full_trace, resumed_trace);
        assert_eq!(
            bits(straight.specific_enthalpy_j_kg()),
            bits(resumed.specific_enthalpy_j_kg())
        );
        assert_eq!(full.integrated_net_input_j.to_bits(), last.integrated_net_input_j.to_bits());
        assert_eq!(full.energy_residual_j.to_bits(), last.energy_residual_j.to_bits());
        // A completed event emits nothing when resumed, even with a stop observer.
        let done = resumed
            .advance_to_observed(cx, p, None, 1.0, 100, |_| {
                panic!("an already reached event must not emit another pair")
            })
            .unwrap();
        assert_eq!(done.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(done.accepted, 0);
    });
}

#[test]
fn cancellation_requested_by_observer_preserves_the_published_pair() {
    support::with_gate(|gate, cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(10.0);
        let curve = sensible();
        let stepper = prepare(cx, &mesh, &curve);
        let mut state = stepper
            .adaptive_state(cx, &[200.0; 4], 0.0, 0.25, config(1e-4))
            .unwrap();
        let mut observed = None;
        let result = state.advance_to_observed(
            cx,
            problem(&mesh, &boundary, &material, &source),
            None,
            1.0,
            100,
            |trial| {
                observed = Some((
                    trial.end_time_s.to_bits(),
                    bits(&trial.second_half.specific_enthalpy_j_kg),
                ));
                gate.request();
                ControlFlow::Continue(())
            },
        );
        assert!(matches!(
            result,
            Err(
                fs_conduction::transient::enthalpy::EnthalpyError::Conduction(
                    fs_conduction::ConductionError::Cancelled { .. }
                )
            )
        ));
        let (time, field) = observed.expect("one complete accepted pair precedes cancellation");
        assert_eq!(state.time_s().to_bits(), time);
        assert_eq!(bits(state.specific_enthalpy_j_kg()), field);
        assert_eq!(state.accepted_steps(), 1);
        assert_eq!(state.endpoint_solves(), 3);
    });
}
