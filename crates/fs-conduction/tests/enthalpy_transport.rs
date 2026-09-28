//! G1/G3 spatial enthalpy checks against independent energy and P1 identities.

mod support;

use fs_blake3::ContentHash;
use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::lumped::{
    BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig, solve_lumped_enthalpy,
};
use fs_conduction::material::{ConductivityModel, ConductivityTable};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig,
};
use fs_conduction::{ConductionError, ConductionProblem, LinearConfig};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const REFERENCE_DENSITY: f64 = 2.0;
const VOLUME: f64 = 1.0 / 6.0;

fn tetrahedron() -> ConductionMesh {
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

fn adiabatic(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

fn curve(knots: &[(f64, f64, f64, f64)]) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x19; 32]),
        knots
            .iter()
            .map(|&(h, t, liquid, density)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: liquid,
                bulk_density_kg_m3: density,
            })
            .collect(),
    )
    .unwrap()
}

fn sensible_curve() -> EquilibriumEnthalpyPhaseCurve {
    // h = 5 (T - 280) J/kg, so rho0 Cp = 10 J/(m3 K).
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x20; 32]),
        SolidLiquidPhase::Solid,
        [(0.0, 280.0), (1100.0, 500.0)]
            .into_iter()
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: REFERENCE_DENSITY,
            })
            .collect(),
    )
    .unwrap()
}

fn latent_curve() -> EquilibriumEnthalpyPhaseCurve {
    // Density varies on purpose: it must not alter stationary reference mass.
    curve(&[
        (0.0, 300.0, 0.0, 4.0),
        (100.0, 350.0, 0.0, 3.0),
        (300.0, 350.0, 1.0, 2.0),
        (600.0, 500.0, 1.0, 1.0),
    ])
}

fn budget() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 4,
        max_elements: 1,
    }
}

fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-9,
    }
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
    values.iter().map(|value| value.to_bits()).collect()
}

#[test]
fn single_phase_transport_matches_existing_backward_euler() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::uniform("uniform source", 4.0).unwrap();
    let curve = sensible_curve();
    let old_t = [310.0, 330.0, 320.0, 300.0];
    let old_h: Vec<_> = old_t.iter().map(|t| 5.0 * (t - 280.0)).collect();
    let saved_h = bits(&old_h);
    with_cx(|cx| {
        let enthalpy =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
        let temperature =
            BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(10.0).unwrap())
                .unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        for dt in [0.01, 0.2, 2.0] {
            let got = enthalpy.advance(cx, p, None, &old_h, dt, config()).unwrap();
            let reference = temperature
                .advance(
                    cx,
                    p,
                    None,
                    &old_t,
                    dt,
                    StepConfig {
                        linear: LinearConfig {
                            tolerance: 1e-12,
                            max_iterations: 64,
                            restart: 4,
                        },
                        energy_tolerance_j: 1e-9,
                    },
                )
                .unwrap();
            for (actual, expected) in got.temperature.iter().zip(&reference.temperature) {
                assert!(
                    (actual - expected).abs() < 1e-9,
                    "dt={dt}: {actual} != {expected}"
                );
            }
            assert!((got.stored_energy_change_j - reference.stored_energy_change_j).abs() < 1e-9);
            assert!(got.newton.converged && got.energy_residual_j.abs() < 1e-9);
        }
        assert_eq!(bits(&old_h), saved_h);
    });
}

#[test]
fn uniform_latent_heating_matches_energy_and_the_lumped_limit() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::uniform("uniform deposition", 100.0).unwrap();
    let curve = latent_curve();
    let mass = REFERENCE_DENSITY * VOLUME;
    // Uniform initial conditions, uniform internal deposition and an insulated
    // boundary make the lumped limit exact: grad(T)=0 and Bi=0 at every step.
    let body = LumpedEnthalpyBody::try_new(
        "uniform synthetic tet",
        mass,
        1.5 + 0.5 * 3.0_f64.sqrt(),
        0.0,
        0.0,
        0.1,
        3.0,
        &curve,
    )
    .unwrap();
    with_cx(|cx| {
        let lumped = solve_lumped_enthalpy(
            cx,
            &body,
            BiotGate::corpus_default(),
            LumpedEnthalpyMarchConfig {
                initial_specific_enthalpy_j_kg: 50.0,
                ambient_temperature_k: 300.0,
                radiation_temperature_k: 300.0,
                internal_power_w: 100.0 * VOLUME,
                duration_s: 6.0,
                maximum_step_s: 1.0,
                maximum_steps: 6,
                enthalpy_tolerance_j_kg: 1e-9,
            },
        )
        .unwrap();
        let stepper =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
        for &nodal_mass in stepper.reference_nodal_masses_kg() {
            assert!((nodal_mass - mass / 4.0).abs() < 1e-15);
        }
        let mut h = vec![50.0; 4];
        let mut energy = 0.0;
        for step in 1..=6 {
            let got = stepper
                .advance(
                    cx,
                    problem(&mesh, &boundary, &material, &source),
                    None,
                    &h,
                    1.0,
                    config(),
                )
                .unwrap();
            let expected_h = 50.0 + 50.0 * step as f64;
            let expected_t = if expected_h <= 100.0 {
                300.0 + 0.5 * expected_h
            } else if expected_h <= 300.0 {
                350.0
            } else {
                350.0 + 0.5 * (expected_h - 300.0)
            };
            let expected_liquid = ((expected_h - 100.0) / 200.0).clamp(0.0, 1.0);
            for i in 0..4 {
                assert!((got.specific_enthalpy_j_kg[i] - expected_h).abs() < 1e-8);
                assert!((got.temperature[i] - expected_t).abs() < 1e-9);
                assert!((got.liquid_mass_fraction[i] - expected_liquid).abs() < 1e-10);
                assert!(
                    (got.specific_enthalpy_j_kg[i]
                        - lumped.samples()[step].phase_state.specific_enthalpy_j_kg())
                    .abs()
                        < 1e-8
                );
            }
            assert!((got.stored_energy_change_j - 100.0 * VOLUME).abs() < 1e-9);
            assert!(got.energy_residual_j.abs() < 1e-9);
            energy += got.stored_energy_change_j;
            h = got.specific_enthalpy_j_kg;
        }
        assert!((energy - 100.0).abs() < 1e-8);
        assert!((mass * (h.iter().sum::<f64>() / 4.0 - 50.0) - energy).abs() < 1e-9);
    });
}

#[test]
fn temperature_dependent_conductivity_matches_an_independent_p1_endpoint() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(300.0, 1.0), (400.0, 11.0)]).unwrap(),
    );
    let source = ScalarField::uniform("uniform deposition", 100.0).unwrap();
    let curve = sensible_curve();
    let old_h = [200.0, 100.0, 100.0, 100.0]; // T = [320, 300, 300, 300].
    with_cx(|cx| {
        let stepper =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
        let got = stepper
            .advance(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                &old_h,
                1.0,
                config(),
            )
            .unwrap();
        // Reference-tet gradients are (-1,-1,-1), e1, e2, e3. Equal cold
        // vertices therefore remain equal. Energy fixes mean(T_new)=315 K;
        // endpoint k=2.5, and the hot/cold contrast obeys
        // D_new = 20 / (1 + 16*k*dt/(rho Cp)) = 4 K.
        for (actual, expected) in got.temperature.iter().zip([318.0, 314.0, 314.0, 314.0]) {
            assert!((actual - expected).abs() < 1e-8, "{actual} != {expected}");
        }
        assert!((got.source_w - 100.0 * VOLUME).abs() < 1e-12);
        assert!((got.stored_energy_change_j - 100.0 * VOLUME).abs() < 1e-9);
        assert!(got.newton.converged && got.energy_residual_j.abs() < 1e-9);
    });
}

#[test]
fn zero_external_heat_preserves_latent_equilibrium_and_mixed_phase_energy() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::uniform("no heat", 0.0).unwrap();
    let curve = latent_curve();
    let old_h = vec![120.0, 180.0, 220.0, 260.0];
    let checkpoint = old_h.clone();
    with_cx(|cx| {
        let stepper =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        let first = stepper.advance(cx, p, None, &old_h, 2.0, config()).unwrap();
        let replay = stepper
            .advance(cx, p, None, &checkpoint, 2.0, config())
            .unwrap();
        assert_eq!(bits(&first.specific_enthalpy_j_kg), bits(&checkpoint));
        assert_eq!(
            bits(&first.specific_enthalpy_j_kg),
            bits(&replay.specific_enthalpy_j_kg)
        );
        assert_eq!(first.newton, replay.newton);
        assert_eq!(first.newton.iterations, 0);
        for (t, fraction) in first.temperature.iter().zip(&first.liquid_mass_fraction) {
            assert!((*t - 350.0).abs() < 1e-12 && *fraction > 0.0 && *fraction < 1.0);
        }
        assert!(first.stored_energy_change_j.abs() < 1e-12);
        assert!(first.energy_residual_j.abs() < 1e-10);

        // The hot vertex is above the plateau, while the three others remain
        // on it. Exact reference-tet balance gives h_hot=310 and h_cold=180:
        // (h_hot-400)/12 + (3*k*V)*0.5*(h_hot-300) = 0. Its sensible
        // dT/dh=0.5 couples INTO the latent rows, whose own dT/dh is zero.
        // One Newton correction suffices on these fixed constitutive branches;
        // scaling tangent rows instead of columns does not satisfy this check.
        let mixed = [400.0, 150.0, 150.0, 150.0];
        let mut one_attempt = config();
        one_attempt.max_newton_iterations = 1;
        let exchange = stepper
            .advance(cx, p, None, &mixed, 1.0, one_attempt)
            .unwrap();
        for (actual, expected) in exchange
            .specific_enthalpy_j_kg
            .iter()
            .zip([310.0, 180.0, 180.0, 180.0])
        {
            assert!((actual - expected).abs() < 1e-9);
        }
        assert_eq!(exchange.newton.iterations, 1);
        assert!(exchange.stored_energy_change_j.abs() < 1e-10);
        assert!(exchange.energy_residual_j.abs() < 1e-10);
    });
    assert_eq!(bits(&old_h), bits(&checkpoint));
}

#[test]
fn invalid_or_cancelled_steps_refuse_without_changing_history() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::uniform("uniform deposition", 100.0).unwrap();
    let curve = latent_curve();
    let old_h = vec![50.0; 4];
    let saved = bits(&old_h);
    with_cx(|cx| {
        assert!(matches!(
            EnthalpyBackwardEuler::uniform(
                cx,
                &mesh,
                &curve,
                REFERENCE_DENSITY,
                EnthalpyBudget {
                    max_vertices: 3,
                    max_elements: 1
                }
            ),
            Err(EnthalpyError::Budget { .. })
        ));
        assert!(EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 0.0, budget()).is_err());
        let stepper =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, REFERENCE_DENSITY, budget()).unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        for invalid in [
            vec![50.0; 3],
            vec![f64::NAN; 4],
            vec![-1.0; 4],
            vec![601.0; 4],
        ] {
            assert!(
                stepper
                    .advance(cx, p, None, &invalid, 1.0, config())
                    .is_err()
            );
        }
        assert!(stepper.advance(cx, p, None, &old_h, 0.0, config()).is_err());
        let mut exhausted = config();
        exhausted.max_newton_iterations = 0;
        assert!(
            stepper
                .advance(cx, p, None, &old_h, 1.0, exhausted)
                .is_err()
        );
        let prescribed = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "fixed temperature",
                |_| true,
                ThermalBc::dirichlet(300.0).unwrap(),
            )
            .unwrap()
            .finish()
            .unwrap();
        assert!(matches!(
            stepper.advance(
                cx,
                problem(&mesh, &prescribed, &material, &source),
                None,
                &old_h,
                1.0,
                config()
            ),
            Err(EnthalpyError::UnsupportedDirichlet)
        ));
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                stepper.advance(cancelled, p, None, &old_h, 1.0, config()),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ));
        });
        assert_eq!(bits(&old_h), saved);
    });
}
