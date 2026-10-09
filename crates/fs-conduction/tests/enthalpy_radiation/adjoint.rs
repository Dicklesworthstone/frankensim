//! Actual perturbed nonlinear solves, not differentiation of a frozen secant.

use super::*;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::{ConductivityTable, LinearConfig, ThermalBoundary};

const OLD: [f64; 4] = [600.0, 1800.0, 3400.0, 5000.0];
const SOURCE: [f64; 4] = [12.0, -3.0, 8.0, 5.0];
const GOAL: [f64; 4] = [0.7, -0.4, 0.2, 1.1];
const DT: f64 = 0.05;

fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-11,
        max_iterations: 64,
        restart: 4,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

struct Fixture {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    curve: EquilibriumEnthalpyPhaseCurve,
    material: ConductivityModel,
}

impl Fixture {
    fn new() -> Self {
        let mesh = regular_tet();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("exterior", |_| true, ThermalBc::robin(HTC, 300.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        // The endpoint has solid, latent and two different sensible slopes.
        // A scalar multiple of a temperature tangent cannot pass this case.
        let knots = [
            (0.0, 250.0, 0.0),
            (1000.0, 350.0, 0.0),
            (3000.0, 350.0, 1.0),
            (4000.0, 450.0, 1.0),
            (6000.0, 550.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: RHO,
        })
        .collect();
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x53; 32]), knots).unwrap();
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(240.0, 2.0), (560.0, 8.0)]).unwrap(),
        );
        Self {
            mesh,
            boundary,
            curve,
            material,
        }
    }

    fn problem<'a>(&'a self, source: &'a ScalarField) -> ConductionProblem<'a> {
        ConductionProblem {
            mesh: &self.mesh,
            boundary: &self.boundary,
            material: &self.material,
            source,
            element_materials: None,
        }
    }

    fn stepper(&self, cx: &fs_exec::Cx<'_>) -> EnthalpyBackwardEuler<'_, '_> {
        EnthalpyBackwardEuler::uniform(cx, &self.mesh, &self.curve, RHO, budget()).unwrap()
    }

    fn objective(
        &self,
        cx: &fs_exec::Cx<'_>,
        old: &[f64],
        q: &[f64],
        dt: f64,
        ambient: f64,
        epsilon: f64,
    ) -> f64 {
        let source = ScalarField::nodal("heat density", 4, q.to_vec()).unwrap();
        let solved = self
            .stepper(cx)
            .advance_with_ambient_radiation(
                cx,
                self.problem(&source),
                None,
                old,
                dt,
                config(),
                &[patch_with_emissivity(ambient, 520.0, epsilon)],
                radiation_config(),
            )
            .unwrap();
        solved
            .conduction
            .temperature
            .iter()
            .zip(GOAL)
            .map(|(t, w)| w * (t - 350.0))
            .sum()
    }
}

fn five_point(step: f64, f: impl Fn(f64) -> f64) -> f64 {
    (-f(2.0 * step) + 8.0 * f(step) - 8.0 * f(-step) + f(-2.0 * step)) / (12.0 * step)
}

fn agree(actual: f64, expected: f64, name: &str) {
    assert!(
        (actual - expected).abs() <= 2e-7 * (1.0 + expected.abs()),
        "{name}: implicit {actual:.12e}, independent perturbation {expected:.12e}"
    );
}

#[test]
fn total_radiative_enthalpy_gradients_match_perturbed_physics_and_transpose() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        for ambient in [300.0, 500.0] {
            let source = ScalarField::nodal("heat density", 4, SOURCE.to_vec()).unwrap();
            let stepper = fixture.stepper(cx);
            let patches = [patch(ambient, 520.0)];
            let linearization = stepper
                .linearize_step_with_ambient_radiation(
                    cx,
                    fixture.problem(&source),
                    None,
                    &OLD,
                    DT,
                    config(),
                    &patches,
                    radiation_config(),
                    8,
                )
                .unwrap();
            let u = [0.2, -0.5, 0.7, 0.3];
            let v = [-1.0, 0.1, 0.6, 0.8];
            let jv = linearization.apply_jacobian(cx, &v).unwrap();
            let jtu = linearization.apply_jacobian_transpose(cx, &u).unwrap();
            assert!((dot(&u, &jv) - dot(&jtu, &v)).abs() < 1e-12);
            let seed = linearization.temperature_pullback(cx, &GOAL).unwrap();
            let gradient = linearization.pullback(cx, &seed, linear()).unwrap();
            assert!(gradient.transport.relative_residual < linear().tolerance);
            assert_eq!(gradient.radiation_regions, ["exterior"]);
            let actual_dual = linearization
                .apply_jacobian_transpose(cx, &gradient.transport.adjoint)
                .unwrap();
            assert!(
                actual_dual
                    .iter()
                    .zip(&seed)
                    .all(|(a, b)| (a - b).abs() < 2e-12)
            );
            for i in 0..4 {
                let expected = five_point(0.03, |delta| {
                    let mut old = OLD;
                    old[i] += delta;
                    fixture.objective(cx, &old, &SOURCE, DT, ambient, EPSILON)
                });
                agree(
                    gradient.transport.previous_specific_enthalpy[i],
                    expected,
                    "history",
                );
                let expected = five_point(0.1, |delta| {
                    let mut q = SOURCE;
                    q[i] += delta;
                    fixture.objective(cx, &OLD, &q, DT, ambient, EPSILON)
                });
                agree(
                    gradient.transport.source_density[i],
                    expected,
                    "source density",
                );
            }
            agree(
                gradient.time_step_s,
                five_point(2e-4, |delta| {
                    fixture.objective(cx, &OLD, &SOURCE, DT + delta, ambient, EPSILON)
                }),
                "time step",
            );
            agree(
                gradient.reservoir_temperatures[0],
                five_point(0.05, |delta| {
                    fixture.objective(cx, &OLD, &SOURCE, DT, ambient + delta, EPSILON)
                }),
                "reservoir temperature",
            );
            agree(
                gradient.emissivities[0],
                five_point(1e-4, |delta| {
                    fixture.objective(cx, &OLD, &SOURCE, DT, ambient, EPSILON + delta)
                }),
                "emissivity",
            );
            let repeat = linearization.pullback(cx, &seed, linear()).unwrap();
            assert_eq!(gradient, repeat);
            let assigned = HeterogeneousEnthalpyBackwardEuler::new(
                cx,
                &fixture.mesh,
                &[ReferenceEnthalpyMaterial {
                    curve: &fixture.curve,
                    reference_density_kg_m3: RHO,
                }],
                &[0],
                budget(),
            )
            .unwrap();
            let facade = assigned
                .linearize_step_with_ambient_radiation(
                    cx,
                    fixture.problem(&source),
                    None,
                    &OLD,
                    DT,
                    config(),
                    &patches,
                    radiation_config(),
                    8,
                )
                .unwrap();
            assert_eq!(facade.pullback(cx, &seed, linear()).unwrap(), gradient);
            // Conversion for existing reverse-history tapes keeps every
            // radiative state derivative, rather than discarding the update.
            let transport = linearization.into_transport();
            assert_eq!(
                transport.pullback(cx, &seed, linear()).unwrap(),
                gradient.transport
            );
        }
    });
}

#[test]
fn latent_plateau_keeps_enthalpy_controls_and_zero_temperature_sensitivity() {
    let mesh = regular_tet();
    let curve = phase_curve(true);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |_| true, ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let area = 3.0_f64.sqrt();
    let mass = RHO * 2.0_f64.sqrt() / 12.0;
    with_cx(|cx| {
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        for ambient in [300.0, 400.0] {
            let dt = 0.2;
            let lin = stepper
                .linearize_step_with_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &[2000.0; 4],
                    dt,
                    config(),
                    &[patch(ambient, 520.0)],
                    radiation_config(),
                    8,
                )
                .unwrap();
            let gradient = lin.pullback(cx, &[0.25; 4], linear()).unwrap();
            let radiation_per_epsilon =
                area * STEFAN_BOLTZMANN_W_M2_K4 * (350.0_f64.powi(4) - ambient.powi(4));
            let power = area * HTC * 50.0 + EPSILON * radiation_per_epsilon;
            agree(gradient.time_step_s, -power / mass, "latent dt");
            agree(
                gradient.emissivities[0],
                -dt * radiation_per_epsilon / mass,
                "latent emissivity",
            );
            agree(
                gradient.reservoir_temperatures[0],
                dt * area * 4.0 * EPSILON * STEFAN_BOLTZMANN_W_M2_K4 * ambient.powi(3) / mass,
                "latent ambient",
            );
            for &value in &gradient.transport.previous_specific_enthalpy {
                agree(value, 0.25, "latent history");
            }
            agree(
                gradient.transport.source_density.iter().sum(),
                dt / RHO,
                "latent source",
            );
            let temperature_seed = lin.temperature_pullback(cx, &[0.25; 4]).unwrap();
            assert_eq!(temperature_seed, [0.0; 4]);
            let zero = lin.pullback(cx, &temperature_seed, linear()).unwrap();
            assert_eq!(zero.transport.adjoint, [0.0; 4]);
            assert_eq!(zero.time_step_s, 0.0);
            assert_eq!(zero.emissivities, [0.0]);
            assert_eq!(zero.reservoir_temperatures, [0.0]);
        }
    });
}

#[test]
fn radiative_adjoint_rechecks_primal_and_refuses_kinks_budgets_and_cancellation() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        let source = ScalarField::Uniform(0.0);
        let problem = fixture.problem(&source);
        let stepper = fixture.stepper(cx);
        let patches = [patch(300.0, 520.0)];
        assert!(
            stepper
                .linearize_step_with_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &OLD,
                    DT,
                    config(),
                    &patches,
                    radiation_config(),
                    7,
                )
                .is_err()
        );
        let mut accepted = stepper
            .advance_with_ambient_radiation(
                cx,
                problem,
                None,
                &OLD,
                DT,
                config(),
                &patches,
                radiation_config(),
            )
            .unwrap()
            .conduction;
        accepted.specific_enthalpy_j_kg[0] += 10.0;
        assert!(matches!(
            stepper.linearize_accepted_with_ambient_radiation(
                cx,
                problem,
                None,
                &OLD,
                DT,
                config(),
                &patches,
                accepted,
                8,
            ),
            Err(EnthalpyAdjointError::PrimalResidual { .. })
        ));
        let lin = stepper
            .linearize_step_with_ambient_radiation(
                cx,
                problem,
                None,
                &OLD,
                DT,
                config(),
                &patches,
                radiation_config(),
                8,
            )
            .unwrap();
        assert!(
            lin.pullback(
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
            assert!(matches!(
                lin.pullback(cancelled, &[1.0; 4], linear()),
                Err(EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
                    ConductionError::Cancelled { .. }
                )))
            ));
        });
        assert!(lin.pullback(cx, &[1.0; 4], linear()).is_ok());
        let equilibrium = ThermalBoundaryBuilder::new(&fixture.mesh)
            .region("exterior", |_| true, ThermalBc::robin(HTC, 350.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        let equilibrium_problem = ConductionProblem {
            boundary: &equilibrium,
            ..problem
        };
        assert!(matches!(
            stepper.linearize_step_with_ambient_radiation(
                cx,
                equilibrium_problem,
                None,
                &[1000.0; 4],
                DT,
                config(),
                &[patch(350.0, 520.0)],
                radiation_config(),
                8,
            ),
            Err(EnthalpyAdjointError::ChartKink { .. })
        ));
    });
}
