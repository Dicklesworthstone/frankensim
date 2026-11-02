//! Source emissivity changes latent storage, with independent P1 moment oracles.
mod support;
#[path = "material_ambient_radiation/adaptive.rs"]
mod adaptive;

use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::heterogeneous::{
    EnthalpyRadiation, HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{
    AmbientRadiationPatch, ConductionMesh, ConductionProblem, ConductivityModel,
    EMISSIVITY_DIMS, MaterialAmbientRadiationPatch, MaterialPointwiseAmbientRadiation,
    PointwiseAmbientRadiation, STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY,
    ScalarField, SurfaceEmissivity, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, Temperature};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const RHO: f64 = 24.0; // Unit right tetra: each reference nodal mass is 1 kg.
const AREA: f64 = 0.5;
const LATENT: f64 = 1000.0;

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0, 1, 2, 3]]), vec![
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    ]).unwrap()
}

fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("surface", |face| face.centroid[2] == 0.0, ThermalBc::adiabatic()).unwrap()
        .remainder("insulated", ThermalBc::adiabatic()).unwrap().finish().unwrap()
}

fn card(axis: Option<&str>) -> MaterialCard {
    let (value, interpolation) = match axis {
        Some(axis) => (PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: if axis == "T" { Temperature::DIMS } else { Dims::NONE },
            knots: if axis == "T" { vec![(300.0, 0.2), (500.0, 0.8)] }
                else { vec![(0.0, 0.2), (1.0, 0.8)] },
            dims: EMISSIVITY_DIMS,
        }, InterpolationPolicy::LinearInside),
        None => (PropertyValue::Scalar { value: 0.5, dims: EMISSIVITY_DIMS },
            InterpolationPolicy::ConstantWithinValidity),
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS), value,
        validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0)
            .with("melt", 0.0, 1.0).with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation, observations: Vec::new(),
        provenance: Provenance { source: "synthetic P1 material radiation oracle".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "chart".into(),
        process: "radiation reference".into(), revision: 0 }, claims, Vec::new()).unwrap()
}

fn curve(card: &MaterialCard, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), [
        (-100.0, 300.0, 0.0), (0.0, 400.0, 0.0),
        (LATENT, 400.0, 1.0), (LATENT + 200.0, 500.0, 1.0),
    ].map(|(h, t, f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: datum + h,
        temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: RHO }).to_vec()).unwrap()
}

fn context() -> QueryPoint {
    QueryPoint::new().with("T", 400.0).unwrap().with("melt", 0.5).unwrap()
        .with("lot", 7.0).unwrap().with("unrelated", 0.5).unwrap()
}

fn binding<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> ThermalMaterialPoint<'a> {
    ThermalMaterialPoint { material: card, context, temperature_axis: "T", liquid_fraction_axis: Some("melt") }
}

fn patch<'a>(card: &'a MaterialCard, curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint, bath: f64) -> MaterialAmbientRadiationPatch<'a>
{
    MaterialAmbientRadiationPatch { region: "surface", ambient_temperature_k: bath, curve,
        material: binding(card, context), selection: MaterialPropertySelection::SingleClaimOnly }
}

fn budget() -> EnthalpyBudget { EnthalpyBudget { max_vertices: 4, max_elements: 1 } }
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig { newton: NewtonKrylovConfig {
        absolute_tolerance: 1e-10, relative_tolerance: 0.0, linear_restart: 4,
        max_linear_cycles: 8, forcing_minimum: 1e-12, forcing_maximum: 1e-4,
        ..NewtonKrylovConfig::default()
    }, max_newton_iterations: 24, energy_tolerance_j: 1e-8 }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}");
}

// Independent exact simplex moments, not the implementation's Duffy rule.
// Expanding Ni * (sum Nj eps_j) * (sum Nk T_k)^4 has total degree six.
fn exact_face_load(temperature: [f64; 3], epsilon: [f64; 3], bath: f64) -> [f64; 3] {
    let factorial = |n: usize| (1..=n).map(|v| v as f64).product::<f64>();
    let mut result = [0.0; 3];
    for a in 0..=4 {
        for b in 0..=4 - a {
            let c = 4 - a - b;
            let alpha = [a, b, c];
            let multinomial = 24.0 / (factorial(a) * factorial(b) * factorial(c));
            let monomial = temperature[0].powi(a as i32) * temperature[1].powi(b as i32)
                * temperature[2].powi(c as i32);
            for (i, out) in result.iter_mut().enumerate() {
                for (j, emissivity) in epsilon.iter().enumerate() {
                    let mut power = alpha; power[i] += 1; power[j] += 1;
                    let moment = 2.0 * AREA * power.map(factorial).into_iter().product::<f64>() / factorial(8);
                    *out += emissivity * multinomial * monomial * moment;
                }
            }
        }
    }
    let sum: f64 = epsilon.iter().sum();
    std::array::from_fn(|i| STEFAN_BOLTZMANN_W_M2_K4
        * (result[i] - bath.powi(4) * AREA / 12.0 * (sum + epsilon[i])))
}

#[test]
fn latent_material_emissivity_changes_storage_at_fixed_kelvin_in_both_heat_directions() {
    let mesh = mesh(); let boundary = boundary(&mesh); let card = card(Some("melt")); let context = context();
    let material = ConductivityModel::isotropic_declared(3.0).unwrap(); let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    for datum in [0.0, 4096.0] {
        let curve = curve(&card, datum);
        for bath in [300.0_f64, 500.0] {
            let dt = 0.2;
            // Three face nodes store heat. The insulated apex stays at its old h.
            let coefficient = AREA * STEFAN_BOLTZMANN_W_M2_K4 * (400.0_f64.powi(4) - bath.powi(4));
            let expected = (3.0 * 500.0 - dt * coefficient * 0.2)
                / (3.0 + dt * coefficient * 0.6 / LATENT);
            let old = [datum + 500.0; 4];
            with_cx(|cx| {
                let radiation = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
                    &[patch(&card, &curve, &context, bath)]).unwrap();
                let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
                    &[ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: RHO }],
                    &[0], budget()).unwrap();
                let result = engine.advance_with_radiation(cx, problem, None, &old, dt, config(),
                    EnthalpyRadiation::MaterialAmbient(&radiation)).unwrap();
                let next = result.conduction();
                for vertex in 0..3 {
                    close(next.specific_enthalpy_j_kg[vertex] - datum, expected, 2e-9);
                    close(next.liquid_mass_fraction[vertex], expected / LATENT, 3e-12);
                }
                close(next.specific_enthalpy_j_kg[3], old[3], 2e-9);
                assert_eq!(next.temperature, vec![400.0; 4]);
                let heat = coefficient * (0.2 + 0.6 * expected / LATENT);
                let report = next.pointwise_radiation.as_ref().unwrap();
                close(report.outward_heat_w, heat, 2e-9);
                close(next.stored_energy_change_j, 3.0 * (expected - 500.0), 6e-9);
                close(next.stored_energy_change_j + dt * heat, 0.0, 6e-9);
                close(next.energy_residual_j, 0.0, config().energy_tolerance_j);
                assert_eq!(next.robin_out_w, 0.0);
                assert_eq!(report.patches[0].material_states.len(), 3);
                for state in &report.patches[0].material_states {
                    close(state.emissivity.value(), 0.2 + 0.6 * expected / LATENT, 3e-12);
                    assert_eq!(state.emissivity.card_identity(), card.content_hash());
                }
                let audit = engine.evaluate_endpoint_with_radiation(cx, problem, None, &old,
                    &next.specific_enthalpy_j_kg, dt, config(), EnthalpyRadiation::MaterialAmbient(&radiation)).unwrap();
                assert!(audit.residual_norm_j <= 1e-8);

                // A phase-dependent accepted report cannot lose its h-feedback by
                // rebinding an otherwise valid fixed-emissivity endpoint adjoint.
                let snapshot = AmbientRadiationPatch::new("surface",
                    SurfaceEmissivity::from_condition("surface", binding(&card, &context),
                        MaterialPropertySelection::SingleClaimOnly).unwrap().0, bath).unwrap();
                let fixed_patches = [snapshot];
                let fixed = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &fixed_patches).unwrap();
                let uniform = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
                assert!(matches!(uniform.linearize_accepted_with_pointwise_ambient_radiation(
                    cx, problem, None, &old, dt, config(), &fixed, next.clone(), 1e-8, 10000,
                ), Err(EnthalpyAdjointError::InvalidInput(message))
                    if message.contains("state-dependent material emissivity")));
            });
            assert_eq!(old, [datum + 500.0; 4]);
        }
    }
}

#[test]
fn nonuniform_face_load_and_complete_enthalpy_tangent_match_independent_moments() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    for axis in ["melt", "T"] {
        let card = card(Some(axis)); let curve = curve(&card, 0.0);
        let h = if axis == "melt" { [LATENT / 6.0, LATENT / 2.0, 5.0 * LATENT / 6.0, 500.0] }
            else { [-80.0, -30.0, 1100.0, 500.0] };
        let t = h.map(|h| curve.state_at_specific_enthalpy(h).unwrap().temperature_k());
        let epsilon: [f64; 3] = std::array::from_fn(|i| if axis == "melt" {
            0.2 + 0.6 * h[i] / LATENT
        } else { 0.2 + 0.003 * (t[i] - 300.0) });
        with_cx(|cx| {
            let radiation = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
                &[patch(&card, &curve, &context, 300.0)]).unwrap();
            let assembly = radiation.assemble(cx, &h).unwrap();
            let expected = exact_face_load([t[0], t[1], t[2]], epsilon, 300.0);
            for (actual, expected) in assembly.outward_load_w[..3].iter().zip(expected) { close(*actual, expected, 2e-10); }
            assert_eq!(assembly.outward_load_w[3], 0.0);
            assert!(assembly.temperature_diagonal_w_k[..3].iter().all(|value| *value > 0.0));
            let direction = [0.7, -0.9, 1.3, 2.0]; let weights = [1.1, 0.3, -0.8, 0.4];
            let mut action = [0.0; 4]; let mut transpose = [0.0; 4];
            assembly.enthalpy_tangent_apply(cx, &direction, &mut action).unwrap();
            assembly.enthalpy_tangent_apply_transpose(cx, &weights, &mut transpose).unwrap();
            close(action.iter().zip(weights).map(|(a, b)| a * b).sum(),
                transpose.iter().zip(direction).map(|(a, b)| a * b).sum(), 2e-12);
            let delta = 1e-3;
            let plus = radiation.assemble(cx, &std::array::from_fn::<_, 4, _>(|i| h[i] + delta * direction[i])).unwrap();
            let minus = radiation.assemble(cx, &std::array::from_fn::<_, 4, _>(|i| h[i] - delta * direction[i])).unwrap();
            for i in 0..4 { close(action[i], (plus.outward_load_w[i] - minus.outward_load_w[i]) / (2.0 * delta), 2e-7); }
            if axis == "melt" {
                let coefficient = AREA * STEFAN_BOLTZMANN_W_M2_K4
                    * (400.0_f64.powi(4) - 300.0_f64.powi(4)) * 0.6 / LATENT / 12.0;
                for i in 0..3 { close(action[i], coefficient * (direction[..3].iter().sum::<f64>() + direction[i]), 2e-13); }
                assert!(action[..3].iter().any(|value| value.abs() > 1e-3));
            }
            assert!(assembly.tangent_apply(cx, &direction, &mut [0.0; 4]).is_err());
            assert!(assembly.parameter_pullback(cx, &weights).is_err());
        });
    }
}

#[test]
fn source_chart_binding_and_original_condition_support_are_not_silently_replaced() {
    let unsupported = card(Some("unrelated")); let unsupported_curve = curve(&unsupported, 0.0);
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let card = card(None); let curve = curve(&card, 0.0); let other_curve = curve.clone();
    let conductivity = ConductivityModel::isotropic_declared(1.0).unwrap(); let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &conductivity,
        source: &source, element_materials: None };
    with_cx(|cx| {
        let radiation = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
            &[patch(&card, &other_curve, &context, 300.0)]).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        assert!(engine.advance_with_material_ambient_radiation(cx, problem, None, &[500.0; 4],
            0.1, config(), &radiation).is_err());
        let bound = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
            &[patch(&card, &curve, &context, 300.0)]).unwrap();
        let got = bound.assemble(cx, &[500.0; 4]).unwrap();
        close(got.report.outward_heat_w, AREA * 0.5 * STEFAN_BOLTZMANN_W_M2_K4
            * (400.0_f64.powi(4) - 300.0_f64.powi(4)), 2e-10);
        assert!(bound.assemble(cx, &[1201.0, 500.0, 500.0, 500.0]).is_err());
        assert!(MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
            &[patch(&unsupported, &unsupported_curve, &context, 300.0)]).is_err());
    });
    with_cancelled_cx(|cx| assert!(MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
        &[patch(&card, &curve, &context, 300.0)]).is_err()));
}

#[test]
fn a_foreign_larger_mesh_refuses_before_indexing_heterogeneous_caloric_bindings() {
    let small = mesh();
    let large = ConductionMesh::new(TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]), vec![
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0],
    ]).unwrap();
    let boundary = boundary(&large); let card = card(Some("melt"));
    let curve = curve(&card, 0.0); let context = context();
    let material = ConductivityModel::isotropic_declared(1.0).unwrap(); let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &large, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    with_cx(|cx| {
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &small,
            &[ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: RHO }], &[0], budget()).unwrap();
        let radiation = MaterialPointwiseAmbientRadiation::new(cx, &large, &boundary,
            &[patch(&card, &curve, &context, 300.0)]).unwrap();
        assert!(engine.advance_with_material_ambient_radiation(cx, problem, None, &[500.0; 8],
            0.1, config(), &radiation).is_err());
        assert!(engine.evaluate_endpoint_with_material_ambient_radiation(cx, problem, None, &[500.0; 8],
            &[500.0; 8], 0.1, config(), &radiation).is_err());
    });
}
