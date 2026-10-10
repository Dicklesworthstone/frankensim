//! Discrete enthalpy sensitivities checked against actual perturbed solves.

mod support;

use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig,
};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_exec::Cx;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const OLD: [f64; 4] = [60.0, 180.0, 240.0, 420.0];
const SOURCE: [f64; 4] = [8.0, -3.0, 12.0, 5.0];
const WEIGHTS: [f64; 4] = [0.7, -0.4, 0.2, 1.1];
const DT: f64 = 0.03;

fn chart(knots: &[(f64, f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x21; 32]),
        knots
            .iter()
            .map(|&(h, temperature, liquid)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: temperature,
                liquid_mass_fraction: liquid,
                bulk_density_kg_m3: 2.0,
            })
            .collect(),
    )
    .unwrap()
}

struct Fixture {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    curve: EquilibriumEnthalpyPhaseCurve,
    material: ConductivityModel,
}

impl Fixture {
    fn new() -> Self {
        // Same reference-tetrahedron fixture as enthalpy_transport.rs.
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        )
        .unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        Self {
            mesh,
            boundary,
            // Solid, latent and liquid slopes are 0.5, 0 and 0.25 K/(J/kg).
            curve: chart(&[
                (0.0, 300.0, 0.0),
                (100.0, 350.0, 0.0),
                (300.0, 350.0, 1.0),
                (700.0, 450.0, 1.0),
            ]),
            material: ConductivityModel::isotropic(
                ConductivityTable::declared_curve(vec![(300.0, 2.0), (450.0, 5.0)]).unwrap(),
            ),
        }
    }

    fn stepper(&self, cx: &Cx<'_>) -> EnthalpyBackwardEuler<'_, '_> {
        EnthalpyBackwardEuler::uniform(
            cx,
            &self.mesh,
            &self.curve,
            2.0,
            EnthalpyBudget {
                max_vertices: 4,
                max_elements: 1,
            },
        )
        .unwrap()
    }

    fn problem<'a>(&'a self, source: &'a ScalarField) -> ConductionProblem<'a> {
        ConductionProblem {
            mesh: &self.mesh,
            boundary: &self.boundary,
            material: &self.material,
            element_materials: None,
            source,
        }
    }

    fn temperature_objective(&self, cx: &Cx<'_>, old: &[f64], q: &[f64]) -> f64 {
        let source = ScalarField::nodal("nodal deposition", 4, q.to_vec()).unwrap();
        let primal = self
            .stepper(cx)
            .advance(cx, self.problem(&source), None, old, DT, config())
            .unwrap();
        // A fixed offset avoids large cancellation in the finite differences.
        primal
            .temperature
            .iter()
            .zip(WEIGHTS)
            .map(|(t, w)| w * (t - 350.0))
            .sum()
    }
}

fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-10,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn five_point(mut evaluate: impl FnMut(f64) -> f64) -> f64 {
    let epsilon = 0.02;
    (-evaluate(2.0 * epsilon) + 8.0 * evaluate(epsilon) - 8.0 * evaluate(-epsilon)
        + evaluate(-2.0 * epsilon))
        / (12.0 * epsilon)
}

fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 64,
        restart: 4,
    }
}

#[test]
fn convective_response_keeps_latent_carry_separate_from_temperature_seeds() {
    let mut fixture = Fixture::new();
    fixture.boundary = ThermalBoundaryBuilder::new(&fixture.mesh)
        .region(
            "air",
            |_| true,
            fs_conduction::ThermalBc::robin(4.0, 330.0).unwrap(),
        )
        .unwrap()
        .finish()
        .unwrap();
    with_cx(|cx| {
        let source = ScalarField::Uniform(0.0);
        let step = fixture
            .stepper(cx)
            .linearize_step(
                cx,
                fixture.problem(&source),
                None,
                &[200.0; 4],
                DT,
                config(),
            )
            .unwrap();
        let response = step.robin_response(cx, &["air"], linear()).unwrap();
        assert_eq!(response.temperature(), &[350.0; 4]);
        let gradient = response
            .pullback(cx, &[0.1; 4], &WEIGHTS, &[0.7], &[0.2])
            .unwrap();
        for (&history, &source) in gradient
            .transport
            .previous_specific_enthalpy
            .iter()
            .zip(&gradient.transport.source_density)
        {
            assert!((history - 0.1).abs() < 1e-14);
            assert!((source - DT * 0.1 / 2.0).abs() < 1e-14);
        }
        let area = fixture
            .mesh
            .boundary()
            .iter()
            .map(|face| face.area)
            .sum::<f64>();
        let lambda = 0.1 / (2.0 / 24.0);
        let reference = (DT * lambda - 0.2) * 4.0 * area;
        assert!((gradient.references[0] - reference).abs() < 1e-12);
        assert!((gradient.log_htc[0] + 20.0 * reference).abs() < 1e-12);
        assert!(
            gradient
                .nodal_load
                .iter()
                .all(|v| (*v - DT * lambda).abs() < 1e-14)
        );
        assert!(step.robin_response(cx, &["air", "air"], linear()).is_err());
        assert!(step.robin_response(cx, &["missing"], linear()).is_err());
        assert!(
            response
                .pullback(cx, &[0.1; 3], &WEIGHTS, &[0.7], &[0.2])
                .is_err()
        );
        assert_eq!(
            gradient,
            response
                .pullback(cx, &[0.1; 4], &WEIGHTS, &[0.7], &[0.2])
                .unwrap()
        );
    });
}

#[test]
fn mixed_phase_jacobian_and_transpose_satisfy_the_dot_identity() {
    let fixture = Fixture::new();
    let source = ScalarField::nodal("nodal deposition", 4, SOURCE.to_vec()).unwrap();
    with_cx(|cx| {
        let stepper = fixture.stepper(cx);
        let response = stepper
            .linearize_step(cx, fixture.problem(&source), None, &OLD, DT, config())
            .unwrap();
        let h = &response.primal().specific_enthalpy_j_kg;
        let slopes: Vec<_> = h
            .iter()
            .map(|&value| {
                fixture
                    .curve
                    .temperature_derivative_at_specific_enthalpy(value)
                    .unwrap()
            })
            .collect();
        assert_eq!(slopes, [0.5, 0.0, 0.0, 0.25]);
        let u = [0.7, -1.2, 0.6, 0.3];
        let v = [-0.5, 0.8, 1.4, -0.2];
        let jv = response.apply_jacobian(cx, &v).unwrap();
        let jtu = response.apply_jacobian_transpose(cx, &u).unwrap();
        assert!((dot(&u, &jv) - dot(&jtu, &v)).abs() < 1e-13);
        // This fixture must expose the nonsymmetry from nonuniform dT/dh
        // and k'(T); passing the primal action off as its transpose cannot pass.
        let ju = response.apply_jacobian(cx, &u).unwrap();
        assert!((dot(&u, &jv) - dot(&ju, &v)).abs() > 1e-4);
    });
}

#[test]
fn nonlinear_history_and_nodal_source_gradients_match_five_point_differences() {
    let fixture = Fixture::new();
    let source = ScalarField::nodal("nodal deposition", 4, SOURCE.to_vec()).unwrap();
    with_cx(|cx| {
        let stepper = fixture.stepper(cx);
        // Linearize the actual accepted forward result, as a reverse time
        // traversal does, rather than differentiating Newton iterations.
        let accepted = stepper
            .advance(cx, fixture.problem(&source), None, &OLD, DT, config())
            .unwrap();
        let response = stepper
            .linearize_accepted(
                cx,
                fixture.problem(&source),
                None,
                &OLD,
                DT,
                config(),
                accepted,
            )
            .unwrap();
        let seed = response.temperature_pullback(cx, &WEIGHTS).unwrap();
        let gradient = response.pullback(cx, &seed, linear()).unwrap();
        assert!(gradient.iterations > 0 && gradient.relative_residual <= 1e-12);
        for vertex in 0..4 {
            let history_fd = five_point(|delta| {
                let mut old = OLD;
                old[vertex] += delta;
                fixture.temperature_objective(cx, &old, &SOURCE)
            });
            let source_fd = five_point(|delta| {
                let mut q = SOURCE;
                q[vertex] += delta;
                fixture.temperature_objective(cx, &OLD, &q)
            });
            assert!(
                (gradient.previous_specific_enthalpy[vertex] - history_fd).abs() < 1e-8,
                "history vertex {vertex}: adjoint={}, FD={history_fd}",
                gradient.previous_specific_enthalpy[vertex]
            );
            assert!(
                (gradient.source_density[vertex] - source_fd).abs() < 1e-8,
                "source vertex {vertex}: adjoint={}, FD={source_fd}",
                gradient.source_density[vertex]
            );
        }
        // Plateau history cannot affect temperature locally, but a source
        // basis there loads its sensible neighbors through consistent P1 mass.
        for vertex in [1, 2] {
            assert!(gradient.previous_specific_enthalpy[vertex].abs() < 1e-12);
            assert!(gradient.source_density[vertex].abs() > 1e-4);
        }
    });
}

#[test]
fn latent_interior_has_exactly_zero_temperature_sensitivity_but_stores_heat() {
    let fixture = Fixture::new();
    let source = ScalarField::uniform("latent heating", 4.0).unwrap();
    let old = [150.0, 180.0, 210.0, 240.0];
    with_cx(|cx| {
        let stepper = fixture.stepper(cx);
        let response = stepper
            .linearize_step(cx, fixture.problem(&source), None, &old, DT, config())
            .unwrap();
        let seed = response.temperature_pullback(cx, &WEIGHTS).unwrap();
        assert_eq!(seed, vec![0.0; 4]);
        let temperature = response.pullback(cx, &seed, linear()).unwrap();
        assert_eq!(temperature.previous_specific_enthalpy, vec![0.0; 4]);
        assert_eq!(temperature.source_density, vec![0.0; 4]);
        assert_eq!(temperature.iterations, 0);
        assert_eq!(
            five_point(|delta| fixture.temperature_objective(cx, &old, &[4.0 + delta; 4])),
            0.0
        );

        // The flat temperature is physical, not a missing derivative: sum(h)
        // still changes with history and source. With rho0=2, d(sum h)/dq_i=dt/2.
        let enthalpy = response.pullback(cx, &[1.0; 4], linear()).unwrap();
        for (previous, source) in enthalpy
            .previous_specific_enthalpy
            .iter()
            .zip(&enthalpy.source_density)
        {
            assert!((previous - 1.0).abs() < 1e-12);
            assert!((source - DT / 2.0).abs() < 1e-12);
        }
    });
}

#[test]
fn genuine_chart_and_conductivity_kinks_refuse_but_collinear_knots_admit() {
    let mut fixture = Fixture::new();
    let source = ScalarField::uniform("no heating", 0.0).unwrap();
    with_cx(|cx| {
        assert!(matches!(
            fixture.stepper(cx).linearize_step(
                cx,
                fixture.problem(&source),
                None,
                &[100.0; 4],
                DT,
                config()
            ),
            Err(EnthalpyAdjointError::ChartKink { .. })
        ));

        fixture.curve = chart(&[
            (0.0, 300.0, 0.0),
            (50.0, 325.0, 0.0),
            (100.0, 350.0, 0.0),
            (300.0, 350.0, 1.0),
            (700.0, 450.0, 1.0),
        ]);
        assert!(
            fixture
                .stepper(cx)
                .linearize_step(cx, fixture.problem(&source), None, &[50.0; 4], DT, config())
                .is_ok()
        );

        fixture.material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(300.0, 1.0), (325.0, 2.0), (400.0, 3.0)])
                .unwrap(),
        );
        assert!(matches!(
            fixture.stepper(cx).linearize_step(
                cx,
                fixture.problem(&source),
                None,
                &[50.0; 4],
                DT,
                config()
            ),
            Err(EnthalpyAdjointError::MaterialKink)
        ));
        fixture.material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(300.0, 1.0), (325.0, 2.0), (400.0, 5.0)])
                .unwrap(),
        );
        assert!(
            fixture
                .stepper(cx)
                .linearize_step(cx, fixture.problem(&source), None, &[50.0; 4], DT, config())
                .is_ok()
        );
    });
}

#[test]
fn failed_or_cancelled_pullbacks_leave_the_accepted_endpoint_reusable() {
    let fixture = Fixture::new();
    let source = ScalarField::nodal("nodal deposition", 4, SOURCE.to_vec()).unwrap();
    with_cx(|cx| {
        let stepper = fixture.stepper(cx);
        let response = stepper
            .linearize_step(cx, fixture.problem(&source), None, &OLD, DT, config())
            .unwrap();
        let saved: Vec<_> = response
            .primal()
            .specific_enthalpy_j_kg
            .iter()
            .map(|value| value.to_bits())
            .collect();
        let seed = response.temperature_pullback(cx, &WEIGHTS).unwrap();
        let expected = response.pullback(cx, &seed, linear()).unwrap();
        assert!(response.temperature_pullback(cx, &[f64::NAN; 4]).is_err());
        assert!(response.pullback(cx, &[1.0; 3], linear()).is_err());
        assert!(
            response
                .pullback(
                    cx,
                    &seed,
                    LinearConfig {
                        max_iterations: 0,
                        ..linear()
                    }
                )
                .is_err()
        );
        assert!(matches!(
            response.pullback(
                cx,
                &seed,
                LinearConfig {
                    max_iterations: 1,
                    restart: 1,
                    ..linear()
                }
            ),
            Err(EnthalpyAdjointError::NotConverged { .. })
        ));
        let mut tampered = response.primal().clone();
        tampered.specific_enthalpy_j_kg[0] += 1.0;
        assert!(
            stepper
                .linearize_accepted(
                    cx,
                    fixture.problem(&source),
                    None,
                    &OLD,
                    DT,
                    config(),
                    tampered
                )
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                response.pullback(cancelled, &seed, linear()),
                Err(EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
                    ConductionError::Cancelled { .. }
                )))
            ));
            assert!(
                response
                    .apply_jacobian_transpose(cancelled, &WEIGHTS)
                    .is_err()
            );
        });
        let replay = response.pullback(cx, &seed, linear()).unwrap();
        assert_eq!(
            replay.previous_specific_enthalpy,
            expected.previous_specific_enthalpy
        );
        assert_eq!(replay.source_density, expected.source_density);
        assert_eq!(
            response
                .primal()
                .specific_enthalpy_j_kg
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            saved
        );
    });
}
