//! Independent P1/scalar balances for implicit spatial radiation and latent heat.

mod support;

use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig,
};
use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, ConductionError, ConductionMesh,
    ConductionProblem, ConductivityModel, EMISSIVITY_DIMS, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, ScalarField, SurfaceEmissivity, ThermalBc, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const RHO: f64 = 10.0;
const CP: f64 = 10.0;
const K: f64 = 3.0;
const HTC: f64 = 3.0;
const EPSILON: f64 = 0.8;

fn regular_tet() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
        ],
    )
    .unwrap()
}

fn phase_curve(latent: bool) -> EquilibriumEnthalpyPhaseCurve {
    let knots: Vec<_> = if latent {
        vec![
            (0.0, 250.0, 0.0),
            (1000.0, 350.0, 0.0),
            (3000.0, 350.0, 1.0),
            (4500.0, 500.0, 1.0),
        ]
    } else {
        vec![(0.0, 200.0, 0.0), (4000.0, 600.0, 0.0)]
    }
    .into_iter()
    .map(|(h, t, liquid)| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: liquid,
        bulk_density_kg_m3: RHO,
    })
    .collect();
    if latent {
        EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x52; 32]), knots).unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x51; 32]),
            SolidLiquidPhase::Solid,
            knots,
        )
        .unwrap()
    }
}

fn patch(ambient: f64, maximum_k: f64) -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: EPSILON,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 240.0, maximum_k),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic implicit-radiation numerical reference".into(),
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
    let emissivity =
        SurfaceEmissivity::from_card("exterior", &card, 350.0, SelectionPolicy::SingleClaimOnly)
            .unwrap();
    AmbientRadiationPatch::new("exterior", emissivity, ambient).unwrap()
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
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-8,
    }
}

fn radiation_config() -> AmbientRadiationConfig {
    AmbientRadiationConfig {
        max_iterations: 128,
        ..AmbientRadiationConfig::default()
    }
}

// Equal vertex masses and equal face areas imply an independent scalar
// equation for the mean on a regular tetrahedron, including nonuniform data.
fn scalar_mean(old_mean: f64, ambient: f64, dt: f64) -> f64 {
    let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
    let area = 3.0_f64.sqrt();
    let mut lo = 240.0_f64;
    let mut hi = 520.0_f64;
    for _ in 0..80 {
        let t = (lo + hi) / 2.0;
        let residual = capacity * (t - old_mean)
            + dt * area
                * (HTC * (t - 300.0)
                    + EPSILON * STEFAN_BOLTZMANN_W_M2_K4 * (t.powi(4) - ambient.powi(4)));
        if residual > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    (lo + hi) / 2.0
}

#[test]
fn implicit_radiation_matches_independent_mean_and_spatial_p1_modes() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |_| true, ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    with_cx(|cx| {
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        for ambient in [300.0, 500.0] {
            for old_t in [[400.0; 4], [420.0, 380.0, 350.0, 450.0]] {
                let old_h: Vec<_> = old_t.iter().map(|t| CP * (t - 200.0)).collect();
                let history = old_h.clone();
                let patches = [patch(ambient, 520.0)];
                let solve = || {
                    stepper
                        .advance_with_ambient_radiation(
                            cx,
                            problem,
                            None,
                            &old_h,
                            0.5,
                            config(),
                            &patches,
                            radiation_config(),
                        )
                        .unwrap()
                };
                let got = solve();
                let mean = scalar_mean(400.0, ambient, 0.5);
                // Analytic stiffness eigenvalue 2*k*V and boundary mass
                // eigenvalue A/12 on the three zero-mean nodal modes.
                let volume = 2.0_f64.sqrt() / 12.0;
                let area = 3.0_f64.sqrt();
                let h_rad = EPSILON * STEFAN_BOLTZMANN_W_M2_K4 * (mean.powi(4) - ambient.powi(4))
                    / (mean - ambient);
                let nodal_capacity = RHO * CP * volume / 4.0;
                let attenuation = nodal_capacity
                    / (nodal_capacity + 0.5 * (2.0 * K * volume + (HTC + h_rad) * area / 12.0));
                for (i, actual) in got.conduction.temperature.iter().enumerate() {
                    let expected = mean + attenuation * (old_t[i] - 400.0);
                    assert!(
                        (actual - expected).abs() < 2e-8,
                        "ambient={ambient}, vertex={i}: {actual} versus {expected}"
                    );
                }
                let nonlinear_heat =
                    area * EPSILON * STEFAN_BOLTZMANN_W_M2_K4 * (mean.powi(4) - ambient.powi(4));
                let convective_heat = area * HTC * (mean - 300.0);
                assert!((got.radiation.nonlinear_radiation_out_w - nonlinear_heat).abs() < 2e-7);
                assert!((got.convective_out_w - convective_heat).abs() < 2e-7);
                assert!(
                    (got.conduction.stored_energy_change_j
                        + 0.5 * (nonlinear_heat + convective_heat))
                        .abs()
                        < 2e-7
                );
                assert!(got.physical_residual_norm_j <= got.physical_residual_tolerance_j);
                assert!(got.physical_energy_residual_j.abs() <= config().energy_tolerance_j);
                assert!(got.radiation.iterations > 1);
                assert_eq!(
                    got.conduction.specific_enthalpy_j_kg,
                    solve().conduction.specific_enthalpy_j_kg
                );
                assert_eq!(old_h, history, "radiative trials must not advance history");
            }
        }
    });
}

#[test]
fn latent_radiative_exchange_changes_phase_without_smoothing_temperature() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |_| true, ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(true);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    with_cx(|cx| {
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        for ambient in [300.0, 400.0] {
            let got = stepper
                .advance_with_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &[2000.0; 4],
                    0.5,
                    config(),
                    &[patch(ambient, 520.0)],
                    radiation_config(),
                )
                .unwrap();
            let outward = 3.0_f64.sqrt()
                * (HTC * 50.0
                    + EPSILON * STEFAN_BOLTZMANN_W_M2_K4 * (350.0_f64.powi(4) - ambient.powi(4)));
            let expected_h = 2000.0 - 0.5 * outward / (RHO * 2.0_f64.sqrt() / 12.0);
            assert!(expected_h > 1000.0 && expected_h < 3000.0);
            for i in 0..4 {
                assert!((got.conduction.temperature[i] - 350.0).abs() < 1e-12);
                assert!((got.conduction.specific_enthalpy_j_kg[i] - expected_h).abs() < 1e-7);
                assert!(
                    (got.conduction.liquid_mass_fraction[i] - (expected_h - 1000.0) / 2000.0).abs()
                        < 1e-10
                );
            }
            assert!(got.physical_energy_residual_j.abs() < 1e-8);
        }
    });
}

#[test]
fn residual_budget_domain_and_cancellation_refusals_leave_history_unchanged() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |_| true, ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let old_h = [2000.0; 4];
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    with_cx(|cx| {
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let loose_outer = AmbientRadiationConfig {
            max_iterations: 1,
            temperature_tolerance_k: 1e6,
            balance_tolerance_w: 1e20,
            ..radiation_config()
        };
        assert!(matches!(
            stepper.advance_with_ambient_radiation(
                cx,
                problem,
                None,
                &old_h,
                0.5,
                config(),
                &[patch(300.0, 520.0)],
                loose_outer,
            ),
            Err(EnthalpyError::Conduction(
                ConductionError::AmbientRadiationNotConverged { .. }
            ))
        ));
        // Even enormous watt/temperature allowances cannot bypass the real
        // nodal physical residual or the unchanged absolute joule budget.
        assert!(
            stepper
                .advance_with_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &old_h,
                    0.5,
                    config(),
                    &[patch(300.0, 375.0)],
                    radiation_config(),
                )
                .is_err()
        );
        let duplicate = [patch(300.0, 520.0), patch(300.0, 520.0)];
        for patches in [&duplicate[..], &[][..]] {
            assert!(
                stepper
                    .advance_with_ambient_radiation(
                        cx,
                        problem,
                        None,
                        &old_h,
                        0.5,
                        config(),
                        patches,
                        radiation_config(),
                    )
                    .is_err()
            );
        }
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                stepper.advance_with_ambient_radiation(
                    cancelled,
                    problem,
                    None,
                    &old_h,
                    0.5,
                    config(),
                    &[patch(300.0, 520.0)],
                    radiation_config(),
                ),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ));
        });
        assert!(
            stepper
                .advance_with_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &old_h,
                    0.5,
                    config(),
                    &[patch(300.0, 520.0)],
                    radiation_config(),
                )
                .is_ok(),
            "refused trials must leave a usable original history"
        );
        assert_eq!(old_h, [2000.0; 4]);
    });
}
