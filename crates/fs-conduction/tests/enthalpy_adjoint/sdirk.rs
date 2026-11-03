//! Two-stage derivatives checked against independently perturbed SDIRK solves.
use super::*;
use fs_conduction::ThermalBc;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA as G};

const Q2: [f64; 4] = [-5.0, 6.0, 2.0, 11.0];
const FIRST_H: [f64; 4] = [0.03, -0.01, 0.02, 0.04];
const LAST_H: [f64; 4] = [-0.02, 0.04, -0.03, 0.01];
const FIRST_T: [f64; 4] = [-0.2, 0.5, 0.3, -0.1];

fn input(problem: ConductionProblem<'_>) -> EnthalpySdirk2Input<'_> {
    EnthalpySdirk2Input {
        problem,
        interfaces: None,
        radiation: EnthalpyRadiation::None,
        contacts: None,
    }
}
fn nodal(values: &[f64]) -> ScalarField {
    ScalarField::nodal("stage source", 4, values.to_vec()).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance={tolerance:e}"
    );
}
fn both_objective(step: &fs_conduction::transient::enthalpy::sdirk::EnthalpySdirk2Step) -> f64 {
    dot(&step.first_stage.specific_enthalpy_j_kg, &FIRST_H)
        + dot(&step.second_stage.specific_enthalpy_j_kg, &LAST_H)
        + step
            .first_stage
            .temperature
            .iter()
            .zip(FIRST_T)
            .map(|(t, w)| (t - 350.0) * w)
            .sum::<f64>()
        + step
            .second_stage
            .temperature
            .iter()
            .zip(WEIGHTS)
            .map(|(t, w)| (t - 350.0) * w)
            .sum::<f64>()
}

#[test]
fn nonlinear_mixed_phase_stage_sources_and_history_match_finite_differences() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        let engine = fixture
            .stepper(cx)
            .with_phase_conductivity(
                cx,
                vec![LiquidMassFractionConductivity::declared(0.4, 2.5).unwrap()],
            )
            .unwrap();
        let sources = [nodal(&SOURCE), nodal(&Q2)];
        let inputs = sources.each_ref().map(|q| input(fixture.problem(q)));
        let response = engine
            .linearize_sdirk2(cx, inputs, &OLD, DT, config(), 0)
            .unwrap();
        let mut seeds = response
            .temperature_pullback(cx, [&FIRST_T, &WEIGHTS])
            .unwrap();
        for i in 0..4 {
            seeds[0][i] += FIRST_H[i];
            seeds[1][i] += LAST_H[i];
        }
        let gradient = response
            .pullback_stages(cx, [&seeds[0], &seeds[1]], linear())
            .unwrap();
        assert!(
            gradient
                .relative_residual
                .iter()
                .all(|r| *r < linear().tolerance)
        );
        assert!(
            gradient
                .iterations
                .iter()
                .all(|n| *n > 0 && *n <= linear().max_iterations)
        );
        let evaluate = |old: &[f64], q1: &[f64], q2: &[f64]| {
            let q = [nodal(q1), nodal(q2)];
            both_objective(
                &engine
                    .advance_sdirk2(
                        cx,
                        q.each_ref().map(|q| input(fixture.problem(q))),
                        old,
                        DT,
                        config(),
                    )
                    .unwrap(),
            )
        };
        for vertex in 0..4 {
            let old_fd = five_point(|delta| {
                let mut old = OLD;
                old[vertex] += delta;
                evaluate(&old, &SOURCE, &Q2)
            });
            near(gradient.previous_specific_enthalpy[vertex], old_fd, 2e-8);
            for stage in 0..2 {
                let source_fd = five_point(|delta| {
                    let mut q = [SOURCE, Q2];
                    q[stage][vertex] += delta;
                    evaluate(&OLD, &q[0], &q[1])
                });
                near(gradient.source_density[stage][vertex], source_fd, 2e-8);
            }
        }
        let v = [[0.7, -0.8, 0.2, 0.4], [-0.3, 0.6, 1.1, -0.5]];
        let w = [[-0.1, 0.4, 0.7, -0.3], [0.8, -0.5, 0.2, 0.6]];
        let action = response.apply_jacobian(cx, [&v[0], &v[1]]).unwrap();
        let transpose = response
            .apply_jacobian_transpose(cx, [&w[0], &w[1]])
            .unwrap();
        near(
            dot(&action[0], &w[0]) + dot(&action[1], &w[1]),
            dot(&v[0], &transpose[0]) + dot(&v[1], &transpose[1]),
            1e-13,
        );
        // Isolate the interstage block: it is -r*M, including latent nodes.
        let coupling = response.apply_jacobian(cx, [&v[0], &[0.0; 4]]).unwrap();
        for i in 0..4 {
            near(
                coupling[1][i],
                -(1.0 - G) / G * (2.0 / 24.0) * v[0][i],
                1e-14,
            );
        }
    });
}

fn bath(mesh: &ConductionMesh, axis: usize, temperature: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "bath",
            |face| face.centroid[axis] == 0.0,
            ThermalBc::dirichlet(temperature).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

#[test]
fn changing_dirichlet_masks_keep_original_history_and_separate_bath_gradients() {
    let fixture = Fixture::new();
    let old = [60.0, 70.0, 65.0, 55.0];
    let source = nodal(&SOURCE);
    with_cx(|cx| {
        let engine = fixture.stepper(cx);
        let boundaries = [bath(&fixture.mesh, 0, 330.0), bath(&fixture.mesh, 1, 335.0)];
        let inputs = boundaries.each_ref().map(|boundary| {
            input(ConductionProblem {
                boundary,
                ..fixture.problem(&source)
            })
        });
        let response = engine
            .linearize_sdirk2(cx, inputs, &old, DT, config(), 0)
            .unwrap();
        let gradient = response
            .pullback_stages(cx, [&FIRST_H, &WEIGHTS], linear())
            .unwrap();
        let objective = |previous: &[f64], temperatures: [f64; 2]| {
            let boundaries = [
                bath(&fixture.mesh, 0, temperatures[0]),
                bath(&fixture.mesh, 1, temperatures[1]),
            ];
            let inputs = boundaries.each_ref().map(|boundary| {
                input(ConductionProblem {
                    boundary,
                    ..fixture.problem(&source)
                })
            });
            let step = engine
                .advance_sdirk2(cx, inputs, previous, DT, config())
                .unwrap();
            dot(&step.first_stage.specific_enthalpy_j_kg, &FIRST_H)
                + dot(&step.second_stage.specific_enthalpy_j_kg, &WEIGHTS)
        };
        for vertex in 0..4 {
            let fd = five_point(|delta| {
                let mut previous = old;
                previous[vertex] += delta;
                objective(&previous, [330.0, 335.0])
            });
            near(gradient.previous_specific_enthalpy[vertex], fd, 2e-8);
        }
        // Node 2 is prescribed at stage one but free at stage two. The first
        // reservoir reaction carries its ORIGINAL history into the final state.
        assert!(gradient.previous_specific_enthalpy[2] < -0.1);
        for stage in 0..2 {
            let fd = five_point(|delta| {
                let mut temperatures = [330.0, 335.0];
                temperatures[stage] += delta;
                objective(&old, temperatures)
            });
            near(
                gradient.prescribed_temperature[stage].iter().sum(),
                fd,
                3e-8,
            );
        }
        assert_eq!(gradient.prescribed_temperature[0][1], 0.0);
        assert_eq!(gradient.prescribed_temperature[1][2], 0.0);
        assert_eq!(gradient.previous_specific_enthalpy[0], 0.0);
        assert_eq!(gradient.previous_specific_enthalpy[3], 0.0);
    });
}

#[test]
fn authentic_stages_admit_when_extrapolated_history_is_outside_the_chart() {
    let mut fixture = Fixture::new();
    fixture.curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x22; 32]),
        fs_material::phase::SolidLiquidPhase::Solid,
        [(0.0, 300.0), (1000.0, 500.0)]
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2.0,
            })
            .to_vec(),
    )
    .unwrap();
    fixture.material = ConductivityModel::isotropic_declared(3.0).unwrap();
    fixture.boundary = bath(&fixture.mesh, 0, 305.0);
    let source = ScalarField::Uniform(0.0);
    let old = [25.0, 100.0, 25.0, 25.0];
    with_cx(|cx| {
        let response = fixture
            .stepper(cx)
            .linearize_sdirk2(
                cx,
                [input(fixture.problem(&source)); 2],
                &old,
                5.0 / 1.2,
                config(),
                0,
            )
            .unwrap();
        let step = response.primal();
        assert!(100.0 + (1.0 - G) / G * (step.first_stage.specific_enthalpy_j_kg[1] - 100.0) < 0.0);
        assert!(step.second_stage.specific_enthalpy_j_kg[1] > 0.0);
        let gradient = response
            .pullback(cx, &[0.0, 1.0, 0.0, 0.0], linear())
            .unwrap();
        let z = -5.0;
        let stability = (1.0 + (1.0 - 2.0 * G) * z) / (1.0 - G * z).powi(2);
        near(gradient.previous_specific_enthalpy[1], stability, 1e-13);
        near(
            step.second_stage.specific_enthalpy_j_kg[1],
            25.0 + 75.0 * stability,
            1e-10,
        );
    });
}

#[test]
fn latent_storage_retains_tableau_weights_and_refused_adjoints_publish_no_gradient() {
    let fixture = Fixture::new();
    let source = ScalarField::Uniform(4.0);
    with_cx(|cx| {
        let engine = fixture.stepper(cx);
        let inputs = [input(fixture.problem(&source)); 2];
        let response = engine
            .linearize_sdirk2(cx, inputs, &[200.0; 4], DT, config(), 0)
            .unwrap();
        let seed = response
            .temperature_pullback(cx, [&WEIGHTS, &WEIGHTS])
            .unwrap();
        assert_eq!(seed, [vec![0.0; 4], vec![0.0; 4]]);
        let zero = response
            .pullback_stages(cx, [&seed[0], &seed[1]], linear())
            .unwrap();
        assert_eq!(zero.previous_specific_enthalpy, [0.0; 4]);
        assert_eq!(zero.iterations, [0; 2]);
        let gradient = response.pullback(cx, &[1.0; 4], linear()).unwrap();
        for i in 0..4 {
            near(gradient.previous_specific_enthalpy[i], 1.0, 2e-14);
            near(gradient.source_density[0][i], (1.0 - G) * DT / 2.0, 2e-14);
            near(gradient.source_density[1][i], G * DT / 2.0, 2e-14);
        }
        let saved = response
            .primal()
            .second_stage
            .specific_enthalpy_j_kg
            .clone();
        assert!(response.pullback(cx, &[1.0; 3], linear()).is_err());
        assert!(
            response
                .pullback_stages(cx, [&[f64::NAN; 4], &[0.0; 4]], linear())
                .is_err()
        );
        assert!(
            response
                .pullback(
                    cx,
                    &[1.0; 4],
                    LinearConfig {
                        max_iterations: 0,
                        ..linear()
                    }
                )
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(response.pullback(cancelled, &[1.0; 4], linear()).is_err());
            assert!(
                engine
                    .linearize_sdirk2(cancelled, inputs, &[200.0; 4], DT, config(), 0)
                    .is_err()
            );
        });
        assert_eq!(response.primal().second_stage.specific_enthalpy_j_kg, saved);
        let zero = ScalarField::Uniform(0.0);
        assert!(matches!(
            engine.linearize_sdirk2(
                cx,
                [input(fixture.problem(&zero)); 2],
                &[100.0; 4],
                DT,
                config(),
                0
            ),
            Err(EnthalpyAdjointError::ChartKink { .. })
        ));
    });
}

#[test]
fn reverse_fixed_grid_trajectory_matches_independent_history_and_forcing_perturbations() {
    let fixture = Fixture::new();
    let d0 = [0.3, -0.4, 0.2, 0.7];
    let dq = [[0.5, -0.1, 0.7, 0.2], [-0.3, 0.6, 0.2, -0.4]];
    with_cx(|cx| {
        let engine = fixture.stepper(cx);
        let sources = [nodal(&SOURCE), nodal(&Q2)];
        let inputs = sources.each_ref().map(|q| input(fixture.problem(q)));
        let mut h = OLD.to_vec();
        let mut tape = Vec::new();
        for dt in [0.01, 0.02, 0.03] {
            let response = engine
                .linearize_sdirk2(cx, inputs, &h, dt, config(), 0)
                .unwrap();
            h.clone_from(&response.primal().second_stage.specific_enthalpy_j_kg);
            tape.push(response);
        }
        let mut carry = WEIGHTS.to_vec();
        let mut forcing = 0.0;
        for response in tape.iter().rev() {
            let gradient = response.pullback(cx, &carry, linear()).unwrap();
            forcing +=
                dot(&gradient.source_density[0], &dq[0]) + dot(&gradient.source_density[1], &dq[1]);
            carry = gradient.previous_specific_enthalpy;
        }
        let evaluate = |initial_shift: f64, forcing_shift: f64| {
            let sources = [SOURCE, Q2]; // Distinct samples remain distinct throughout the history.
            let q = [
                nodal(&std::array::from_fn::<_, 4, _>(|i| {
                    sources[0][i] + forcing_shift * dq[0][i]
                })),
                nodal(&std::array::from_fn::<_, 4, _>(|i| {
                    sources[1][i] + forcing_shift * dq[1][i]
                })),
            ];
            let inputs = q.each_ref().map(|q| input(fixture.problem(q)));
            let mut h: Vec<_> = OLD
                .iter()
                .zip(d0)
                .map(|(h, d)| h + initial_shift * d)
                .collect();
            for dt in [0.01, 0.02, 0.03] {
                h = engine
                    .advance_sdirk2(cx, inputs, &h, dt, config())
                    .unwrap()
                    .second_stage
                    .specific_enthalpy_j_kg;
            }
            dot(&h, &WEIGHTS)
        };
        near(
            dot(&carry, &d0),
            five_point(|delta| evaluate(delta, 0.0)),
            2e-8,
        );
        near(forcing, five_point(|delta| evaluate(0.0, delta)), 2e-8);
    });
}
