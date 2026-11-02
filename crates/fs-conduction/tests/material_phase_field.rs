//! Original source k(T,f) reaches the native P1 enthalpy field and its tangent.
//! Constitutive numbers below are manufactured; they are not measured material data.
mod support;

use fs_conduction::material::LiquidMassFractionConductivity;
use fs_conduction::transient::enthalpy::adaptive::{
    AdaptiveEnthalpyConfig, AdaptiveEnthalpyStatus,
};
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::material::phase::MaterialPhaseTransport;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyError, EnthalpyStepConfig};
use fs_conduction::{
    ConductionMesh, ConductivityModel, ElementMaterials, LinearConfig, MaterialId, MaterialTable,
    ScalarField, ThermalBc, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const OLD: [f64; 4] = [-20.0, 80.0, 20.0, 140.0];
const Q: [f64; 4] = [8.0, -3.0, 12.0, 5.0];
const DT: f64 = 0.01;
const MASS: f64 = 2.0 / 24.0;

fn card(values: [f64; 3]) -> MaterialCard {
    conductivity_card(PropertyValue::Curve {
        abscissa: "f".into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(0.0, values[0]), (0.5, values[1]), (1.0, values[2])],
        dims: KD,
    })
}

fn conductivity_card(value: PropertyValue) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("k", KD),
            value,
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    800.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
                .with_quantity(
                    "p",
                    QuantitySpec::dimensional(PRESSURE),
                    100_000.0,
                    200_000.0,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured source-bound phase-field fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic phase transport".into(),
            phase: "declared equilibrium chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn caloric(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (100.0, 400.0, 1.0),
            (500.0, 800.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            // Current bulk density is deliberately not reference rho0=2.
            bulk_density_kg_m3: 1000.0,
        })
        .collect(),
    )
    .unwrap()
}

fn prepare<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    range: [f64; 2],
) -> PreparedPhaseConductivity<'a> {
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(PRESSURE), 150_000.0)
        .unwrap();
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        range,
        || Ok(()),
    )
    .unwrap()
}

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

fn budget() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 4,
        max_elements: 1,
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

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn temperature(h: f64) -> f64 {
    if h < 0.0 {
        400.0 + h
    } else if h <= 100.0 {
        400.0
    } else {
        300.0 + h
    }
}

// Independent truncated-power integral for four DISTINCT nodal enthalpies:
// E[(h_h-c)+] = sum_i (h_i-c)+^4 / (4 product_{j!=i}(h_i-h_j)).
// This closed polynomial formula differs from the production B-spline/Gauss leaf.
fn average_conductivity(h: &[f64]) -> f64 {
    let hinge = |cut: f64| {
        (0..4)
            .map(|i| {
                let denominator = (0..4)
                    .filter(|&j| i != j)
                    .map(|j| h[i] - h[j])
                    .product::<f64>();
                assert_ne!(denominator, 0.0);
                (h[i] - cut).max(0.0).powi(4) / (4.0 * denominator)
            })
            .sum::<f64>()
    };
    // Original f knots: k=2 at f=0, k=8 at f=1/2, k=3 at f=1.
    2.0 + 0.12 * hinge(0.0) - 0.22 * hinge(50.0) + 0.1 * hinge(100.0)
}

fn reference_residual(h: &[f64], old: &[f64]) -> [f64; 4] {
    let t: Vec<_> = h.iter().map(|&h| temperature(h)).collect();
    let laplace = [
        3.0 * t[0] - t[1] - t[2] - t[3],
        t[1] - t[0],
        t[2] - t[0],
        t[3] - t[0],
    ];
    let k = average_conductivity(h);
    // Right tetrahedron has V=1/6; integral N_i N_j = (1+delta_ij)/120.
    std::array::from_fn(|i| {
        MASS * (h[i] - old[i])
            + DT * (k * laplace[i] / 6.0 - (Q.iter().sum::<f64>() + Q[i]) / 120.0)
    })
}

#[test]
fn source_knots_are_integrated_in_enthalpy_space_including_negative_enthalpy() {
    let mesh = mesh();
    let card = card([2.0, 8.0, 3.0]);
    let curve = caloric(&card);
    let source = prepare(&card, &curve, [-100.0, 500.0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        assert!(std::ptr::eq(law.conductivity_source(), &source));
        assert!(std::ptr::eq(
            law.conductivity_source().caloric_curve(),
            &curve
        ));
        assert_eq!(
            source.knots().iter().map(|k| k[0]).collect::<Vec<_>>(),
            vec![-100.0, 0.0, 50.0, 100.0, 500.0]
        );
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let h = [-20.0, 80.0, -20.0, -20.0];
        let endpoint = solver
            .evaluate_endpoint(
                cx,
                law.problem(&mesh, &boundary, &q),
                None,
                &h,
                &h,
                DT,
                config(),
                None,
            )
            .unwrap();
        // Here h_h=-20+100*x. The normalized volume density is 3*(1-x)^2.
        // Integrate the two source hinges exactly, across h=0 and h=50.
        let mean_k = 2.0 + 0.12 * 80.0_f64.powi(4) / (4.0 * 100.0_f64.powi(3))
            - 0.22 * 30.0_f64.powi(4) / (4.0 * 100.0_f64.powi(3));
        near(mean_k, 3.18425, 1e-14);
        assert!((mean_k - 2.75).abs() > 0.4); // mean of nodal k is wrong
        assert!((mean_k - 2.6).abs() > 0.5); // k at mean h is also wrong
        near(endpoint.residual_j[0], -DT * mean_k * 20.0 / 6.0, 1e-12);
        near(endpoint.residual_j[1], DT * mean_k * 20.0 / 6.0, 1e-12);
        near(endpoint.residual_j[2], 0.0, 1e-12);
        near(endpoint.residual_j[3], 0.0, 1e-12);
        near(endpoint.energy_residual_j, 0.0, 1e-12);
        let flux = law.element_heat_flux(cx, &mesh, &h).unwrap();
        near(flux[0][0], -mean_k * 20.0, 1e-11);
        near(flux[0][1], 0.0, 1e-11);
        near(flux[0][2], 0.0, 1e-11);
        assert!(law.element_heat_flux(cx, &mesh, &h[..3]).is_err());
        assert_eq!(solver.reference_nodal_masses_kg(), &[MASS; 4]);
    });
}

#[test]
fn latent_source_transport_balance_fresh_tangent_and_history_adjoint_agree() {
    let mesh = mesh();
    let card = card([2.0, 8.0, 3.0]);
    let curve = caloric(&card);
    let source = prepare(&card, &curve, [-100.0, 500.0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::nodal("manufactured deposition", 4, Q.to_vec()).unwrap();
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &q);
        let response = solver
            .linearize_step(cx, problem, None, &OLD, DT, config())
            .unwrap();
        let h = &response.primal().specific_enthalpy_j_kg;
        assert!(h[0] < 0.0 && h[3] > 100.0);
        assert!(h[1] > 50.0 && h[1] < 100.0 && h[2] > 0.0 && h[2] < 50.0);
        for residual in reference_residual(h, &OLD) {
            near(residual, 0.0, 2e-10);
        }
        near(
            response.primal().stored_energy_change_j,
            DT * Q.iter().sum::<f64>() / 24.0,
            1e-10,
        );
        for column in 0..4 {
            let mut direction = [0.0; 4];
            direction[column] = 1.0;
            let action = response.apply_jacobian(cx, &direction).unwrap();
            let mut plus = h.clone();
            let mut minus = h.clone();
            plus[column] += 0.01;
            minus[column] -= 0.01;
            // Reassemble physical residuals at both fresh enthalpy fields.
            let rp = solver
                .evaluate_endpoint(cx, problem, None, &OLD, &plus, DT, config(), None)
                .unwrap();
            let rm = solver
                .evaluate_endpoint(cx, problem, None, &OLD, &minus, DT, config(), None)
                .unwrap();
            let independent_plus = reference_residual(&plus, &OLD);
            let independent_minus = reference_residual(&minus, &OLD);
            for row in 0..4 {
                near(
                    action[row],
                    (rp.residual_j[row] - rm.residual_j[row]) / 0.02,
                    2e-9,
                );
                near(
                    action[row],
                    (independent_plus[row] - independent_minus[row]) / 0.02,
                    2e-9,
                );
            }
            if column == 1 || column == 2 {
                // On the latent plateau dT/dh=0, but dk/dh still carries heat sensitivity.
                assert!(action[0].abs() > 1e-6);
            }
        }
        let weights = [0.7, -0.4, 0.2, 1.1];
        let seed = response.temperature_pullback(cx, &weights).unwrap();
        assert_eq!([seed[1], seed[2]], [0.0, 0.0]);
        let gradient = response
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-12,
                    max_iterations: 64,
                    restart: 4,
                },
            )
            .unwrap();
        for column in 0..4 {
            let objective = |delta| {
                let mut old = OLD;
                old[column] += delta;
                let next = solver
                    .advance(cx, problem, None, &old, DT, config())
                    .unwrap();
                next.temperature
                    .iter()
                    .zip(weights)
                    .map(|(t, w)| (t - 400.0) * w)
                    .sum::<f64>()
            };
            let e = 0.02;
            let fd = (-objective(2.0 * e) + 8.0 * objective(e) - 8.0 * objective(-e)
                + objective(-2.0 * e))
                / (12.0 * e);
            near(gradient.previous_specific_enthalpy[column], fd, 1e-8);
        }
    });
}

#[test]
fn native_adaptive_storage_crosses_latent_heat_once_with_fixed_reference_mass() {
    let mesh = mesh();
    let card = card([2.0, 8.0, 3.0]);
    let curve = caloric(&card);
    let source = prepare(&card, &curve, [-100.0, 500.0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::Uniform(100.0);
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &q);
        let masses = solver.reference_nodal_masses_kg().to_vec();
        let mut state = solver
            .adaptive_state(
                cx,
                &[-20.0; 4],
                0.0,
                0.25,
                AdaptiveEnthalpyConfig {
                    enthalpy_absolute_tolerance_j_kg: 1e-6,
                    temperature_absolute_tolerance_k: 1e-6,
                    liquid_fraction_absolute_tolerance: 1e-8,
                    relative_tolerance: 0.0,
                    minimum_step_s: 1e-6,
                    maximum_step_s: 0.4,
                    endpoint: config(),
                },
            )
            .unwrap();
        // Uniform adiabatic heating gives h=-20+(q/rho0)*t exactly.
        for (time, expected_h, expected_t, expected_f) in
            [(1.4, 50.0, 400.0, 0.5), (3.0, 130.0, 430.0, 1.0)]
        {
            let report = state.advance_to(cx, problem, None, time, 100).unwrap();
            assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
            for vertex in 0..4 {
                near(state.specific_enthalpy_j_kg()[vertex], expected_h, 1e-8);
                near(state.temperature()[vertex], expected_t, 1e-8);
                near(state.liquid_mass_fraction()[vertex], expected_f, 1e-10);
            }
            near(state.last_endpoint().unwrap().energy_residual_j, 0.0, 1e-10);
        }
        assert_eq!(solver.reference_nodal_masses_kg(), masses.as_slice());
        near(
            masses.iter().sum::<f64>() * (130.0 + 20.0),
            100.0 * 3.0 / 6.0,
            1e-12,
        );
    });
}

#[test]
fn a_uniform_source_corner_is_forward_valid_but_has_no_claimed_classical_adjoint() {
    let mesh = mesh();
    let card = card([2.0, 8.0, 3.0]);
    let curve = caloric(&card);
    let source = prepare(&card, &curve, [-100.0, 500.0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &q);
        let old = [50.0; 4];
        let forward = solver
            .advance(cx, problem, None, &old, DT, config())
            .unwrap();
        assert_eq!(forward.specific_enthalpy_j_kg, old);
        assert!(matches!(
            solver.linearize_accepted(cx, problem, None, &old, DT, config(), forward),
            Err(EnthalpyAdjointError::MaterialKink)
        ));
    });
}

#[test]
fn source_domain_basis_double_scaling_and_execution_guards_are_preserved() {
    let mesh = mesh();
    let card = card([2.0, 8.0, 3.0]);
    let curve = caloric(&card);
    let source = prepare(&card, &curve, [-10.0, 200.0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        assert!(MaterialPhaseTransport::try_new_isotropic(cx, &source, 0.0).is_err());
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let other_mass = MaterialPhaseTransport::try_new_isotropic(cx, &source, 3.0).unwrap();
        assert_ne!(law.identity(), other_mass.identity());
        assert!(matches!(
            law.prepare_uniform(
                cx,
                &mesh,
                EnthalpyBudget {
                    max_vertices: 3,
                    max_elements: 1,
                }
            ),
            Err(EnthalpyError::Budget { .. })
        ));
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &q);
        let outside = [-20.0; 4]; // Caloric chart supports it; original k source does not.
        assert!(matches!(
            solver.advance(cx, problem, None, &outside, DT, config()),
            Err(EnthalpyError::SourceEnthalpyOutsideDomain { vertex: 0, .. })
        ));
        assert_eq!(outside, [-20.0; 4]);
        let valid = [25.0; 4];
        let unrelated = ConductivityModel::isotropic_declared(1.0).unwrap();
        let mut replaced = problem;
        replaced.material = &unrelated;
        assert!(matches!(
            solver.advance(cx, replaced, None, &valid, DT, config()),
            Err(EnthalpyError::InvalidInput(_))
        ));
        let overrides = ElementMaterials::new(
            MaterialTable::new([(MaterialId(0), unrelated.clone())]).unwrap(),
            vec![MaterialId(0)],
        )
        .unwrap();
        let mut replaced = problem;
        replaced.element_materials = Some(&overrides);
        assert!(
            solver
                .advance(cx, replaced, None, &valid, DT, config())
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(MaterialPhaseTransport::try_new_isotropic(cancelled, &source, 2.0).is_err());
            assert!(
                solver
                    .advance(cancelled, problem, None, &valid, DT, config())
                    .is_err()
            );
        });
        assert!(matches!(
            solver.with_phase_conductivity(
                cx,
                vec![LiquidMassFractionConductivity::declared(1.0, 2.0).unwrap()],
            ),
            Err(EnthalpyError::InvalidInput(_))
        ));
    });
    let zero_card = self::card([0.0; 3]);
    let zero_curve = caloric(&zero_card);
    let zero = prepare(&zero_card, &zero_curve, [-100.0, 500.0]);
    with_cx(|cx| {
        assert!(MaterialPhaseTransport::try_new_isotropic(cx, &zero, 2.0).is_err());
    });
}

#[test]
fn a_strongly_descending_source_preserves_its_positive_upper_endpoint() {
    let mesh = mesh();
    let card = card([1e20, 1e20, 1e-20]);
    let curve = caloric(&card);
    // Clip at the end of the descending source segment, so the uniform
    // endpoint uses its left segment rather than a constant liquid branch.
    let source = prepare(&card, &curve, [50.0, 100.0]);
    assert_eq!(
        source
            .state_at_specific_enthalpy(100.0)
            .unwrap()
            .conductivity_w_m_k()
            .to_bits(),
        1e-20_f64.to_bits(),
    );
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &q);
        let h = [100.0; 4];
        // Both residual-only and tangent-producing native paths must retain
        // the admitted positive k. Lower-end interpolation used to round it
        // to zero and spuriously refuse this original supported source state.
        let endpoint = solver
            .evaluate_endpoint(cx, problem, None, &h, &h, DT, config(), None)
            .unwrap();
        near(endpoint.residual_norm_j, 0.0, 1e-30);
        let next = solver.advance(cx, problem, None, &h, DT, config()).unwrap();
        assert_eq!(next.specific_enthalpy_j_kg, h);
        assert_eq!(
            law.element_heat_flux(cx, &mesh, &h).unwrap(),
            vec![[0.0; 3]]
        );
    });
}

#[test]
fn source_domain_overshoots_backtrack_without_extrapolating_the_physical_source() {
    let mesh = mesh();
    let card = conductivity_card(PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(400.0, 0.1), (401.0, 100.1)],
        dims: KD,
    });
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [(-1.0, 399.0), (5.0, 405.0)]
            .into_iter()
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2.0,
            })
            .collect(),
    )
    .unwrap();
    let source = prepare(&card, &curve, [0.0, 1.0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "fixed bath",
            |face| face.centroid[0] == 0.0,
            ThermalBc::dirichlet(400.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::Uniform(82.96);
    let old = [0.0, 0.1, 0.0, 0.0];
    // For the single free vertex, T=400+h and <k>=0.1+25h.
    // The exact scalar residual has supported root h=0.9, but its first
    // full Newton proposal leaves k's support while remaining caloric-valid.
    let residual = |h: f64| (h - 0.1) / 12.0 + h * (0.1 + 25.0 * h) / 6.0 - 82.96 / 24.0;
    let first_trial = 0.1 - residual(0.1) / (1.0 / 12.0 + (0.1 + 50.0 * 0.1) / 6.0);
    assert!(first_trial > 1.0 && first_trial < 5.0);
    curve.state_at_specific_enthalpy(first_trial).unwrap();
    near(residual(0.9), 0.0, 1e-14);
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &q);
        let unsupported = [0.0, first_trial, 0.0, 0.0];
        assert!(matches!(
            solver.evaluate_endpoint(cx, problem, None, &old, &unsupported, 1.0, config(), None),
            Err(EnthalpyError::SourceEnthalpyOutsideDomain { vertex: 1, .. })
        ));
        assert!(matches!(
            solver.advance(cx, problem, None, &unsupported, 1.0, config()),
            Err(EnthalpyError::SourceEnthalpyOutsideDomain { vertex: 1, .. })
        ));
        let next = solver
            .advance(cx, problem, None, &old, 1.0, config())
            .unwrap();
        near(next.specific_enthalpy_j_kg[1], 0.9, 1e-9);
        for vertex in [0, 2, 3] {
            near(next.specific_enthalpy_j_kg[vertex], 0.0, 1e-12);
        }
        assert!(
            next.newton
                .history
                .iter()
                .any(|iteration| iteration.step_length > 0.0 && iteration.step_length < 1.0)
        );
        near(next.energy_residual_j, 0.0, 1e-10);
        let endpoint = solver
            .evaluate_endpoint(
                cx,
                problem,
                None,
                &old,
                &next.specific_enthalpy_j_kg,
                1.0,
                config(),
                None,
            )
            .unwrap();
        near(endpoint.residual_norm_j, 0.0, 1e-10);
        assert_eq!(old, [0.0, 0.1, 0.0, 0.0]);
    });
}
