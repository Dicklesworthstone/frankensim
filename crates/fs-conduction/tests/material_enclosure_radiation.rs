//! Independent gray-patch equations with spatial material emissivity and latent h.
mod support;
#[path = "material_enclosure_radiation/adaptive.rs"]
mod adaptive;
#[path = "material_enclosure_radiation/adjoint.rs"]
mod adjoint;

use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::heterogeneous::{
    EnthalpyRadiation, HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, EMISSIVITY_DIMS, GrayDiffuseEnclosure,
    MaterialEnclosureRadiationPatch, MaterialPointwiseEnclosureRadiation, OpenDiffuseEnclosure,
    PointwiseEnclosureRadiation, PointwiseEnclosureTolerance, RadiationSurface, ScalarField,
    STEFAN_BOLTZMANN_W_M2_K4 as SIGMA, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    ThermalBc, ThermalBoundary, ThermalBoundaryBuilder, ViewFactorEvidence, ViewFactorMatrix,
    ViewFactorTolerance,
};
use fs_evidence::ValidityDomain;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const AREA: f64 = 0.5;
const RHO: f64 = 24.0; // Each right tetra has four 1-kg reference nodal masses.
const LATENT: f64 = 1000.0;

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]), vec![
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0],
    ]).unwrap()
}
fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("left", |face| face.centroid[2] == 0.0 && face.vertices[0] < 4, ThermalBc::adiabatic()).unwrap()
        .region("right", |face| face.centroid[2] == 0.0 && face.vertices[0] >= 4, ThermalBc::adiabatic()).unwrap()
        .remainder("insulated", ThermalBc::adiabatic()).unwrap().finish().unwrap()
}
fn material_card(name: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "melt".into(), abscissa_dims: Dims::NONE,
            knots: vec![(0.0, 0.2), (1.0, 0.8)], dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0)
            .with("melt", 0.0, 1.0).with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(), provenance: Provenance { source: "synthetic material enclosure oracle".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "chart".into(),
        process: "radiation fixture".into(), revision: 0 }, claims, Vec::new()).unwrap()
}
fn phase_curve(card: &MaterialCard, melting_k: f64, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), [
        (-100.0, 300.0, 0.0), (0.0, melting_k, 0.0),
        (LATENT, melting_k, 1.0), (LATENT + 200.0, 500.0, 1.0),
    ].map(|(h, t, f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: datum + h,
        temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: RHO }).to_vec()).unwrap()
}
fn context() -> QueryPoint {
    QueryPoint::new().with("T", 400.0).unwrap().with("melt", 0.5).unwrap().with("lot", 7.0).unwrap()
}
fn binding<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> ThermalMaterialPoint<'a> {
    ThermalMaterialPoint { material: card, context, temperature_axis: "T", liquid_fraction_axis: Some("melt") }
}
fn bindings<'a>(cards: &'a [MaterialCard; 2], curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    context: &'a QueryPoint) -> [MaterialEnclosureRadiationPatch<'a>; 2]
{
    [0, 1].map(|i| MaterialEnclosureRadiationPatch { region: ["left", "right"][i], curve: &curves[i],
        material: binding(&cards[i], context), selection: MaterialPropertySelection::SingleClaimOnly })
}
fn surfaces(mesh: &ConductionMesh, cards: &[MaterialCard; 2], context: &QueryPoint) -> Vec<RadiationSurface> {
    [0, 1].map(|body| {
        let name = ["left", "right"][body];
        let snapshot = SurfaceEmissivity::from_condition(name, binding(&cards[body], context),
            MaterialPropertySelection::SingleClaimOnly).unwrap().0;
        RadiationSurface::new(mesh, name, |face| face.centroid[2] == 0.0
            && (face.vertices[0] as usize / 4) == body, snapshot).unwrap()
    }).to_vec()
}
fn closed(mesh: &ConductionMesh, cards: &[MaterialCard; 2], context: &QueryPoint) -> GrayDiffuseEnclosure {
    // Algebraic parallel-patch view model; tetrahedra provide traces and storage.
    GrayDiffuseEnclosure::new(surfaces(mesh, cards, context), ViewFactorMatrix::infinite_parallel_plates(AREA).unwrap()).unwrap()
}
fn budget() -> EnthalpyBudget { EnthalpyBudget { max_vertices: 8, max_elements: 2 } }
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig { newton: NewtonKrylovConfig { absolute_tolerance: 1e-9, relative_tolerance: 0.0,
        linear_restart: 8, max_linear_cycles: 12, forcing_minimum: 1e-12, forcing_maximum: 1e-4,
        ..NewtonKrylovConfig::default() }, max_newton_iterations: 32, energy_tolerance_j: 1e-8 }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.16e} != {expected:.16e}, tolerance={tolerance:e}");
}

// Independent simplex moments of Ni * epsilon_h * T_h^4, not quadrature.
fn emitted_load(t: [f64; 3], epsilon: [f64; 3]) -> [f64; 3] {
    let factorial = |n: usize| (1..=n).map(|v| v as f64).product::<f64>();
    let mut result = [0.0; 3];
    for a in 0..=4 { for b in 0..=4 - a {
        let c = 4 - a - b; let alpha = [a, b, c];
        let term = 24.0 / (factorial(a) * factorial(b) * factorial(c))
            * t[0].powi(a as i32) * t[1].powi(b as i32) * t[2].powi(c as i32);
        for (i, out) in result.iter_mut().enumerate() { for (j, eps) in epsilon.iter().enumerate() {
            let mut power = alpha; power[i] += 1; power[j] += 1;
            *out += SIGMA * eps * term * 2.0 * AREA
                * power.map(factorial).into_iter().product::<f64>() / factorial(8);
        } }
    } }
    result
}

fn oracle(curves: &[EquilibriumEnthalpyPhaseCurve; 2], h: &[f64; 8], mutual: f64, bath: f64)
    -> ([f64; 8], [f64; 2], f64)
{
    let mut emitted = [[0.0; 3]; 2]; let mut weights = [[0.0; 3]; 2]; let mut means = [0.0; 2];
    for body in 0..2 {
        let phases = [0, 1, 2].map(|i| curves[body].state_at_specific_enthalpy(h[body * 4 + i]).unwrap());
        let t = phases.map(|p| p.temperature_k()); let epsilon = phases.map(|p| 0.2 + 0.6 * p.liquid_mass_fraction());
        emitted[body] = emitted_load(t, epsilon); means[body] = epsilon.iter().sum::<f64>() / 3.0;
        weights[body] = std::array::from_fn(|i| AREA / 12.0 * (epsilon.iter().sum::<f64>() + epsilon[i]));
    }
    let bath_power = SIGMA * bath.powi(4);
    let rhs: [f64; 2] = std::array::from_fn(|i| emitted[i].iter().sum::<f64>() / AREA
        + (1.0 - means[i]) * (1.0 - mutual) * bath_power);
    let a = (1.0 - means[0]) * mutual; let b = (1.0 - means[1]) * mutual;
    let radiosity = [(rhs[0] + a * rhs[1]) / (1.0 - a * b), (rhs[1] + b * rhs[0]) / (1.0 - a * b)];
    let irradiation = [mutual * radiosity[1] + (1.0 - mutual) * bath_power,
        mutual * radiosity[0] + (1.0 - mutual) * bath_power];
    let mut load = [0.0; 8]; let mut heat = [0.0; 2];
    for body in 0..2 { for i in 0..3 {
        load[body * 4 + i] = emitted[body][i] - weights[body][i] * irradiation[body];
        heat[body] += load[body * 4 + i];
    } }
    let bath_heat = 2.0 * AREA * (1.0 - mutual) * (bath_power - (radiosity[0] + radiosity[1]) / 2.0);
    (load, heat, bath_heat)
}

#[test]
fn material_reflection_matches_spatial_moments_and_the_full_nonsymmetric_h_action() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let cards = [material_card("left"), material_card("right")];
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 350.0, 0.0)];
    let declarations = bindings(&cards, &curves, &context); let enclosure = closed(&mesh, &cards, &context);
    let h = [-30.0, 500.0, 1100.0, 500.0, -60.0, 200.0, 1050.0, 500.0];
    with_cx(|cx| {
        for mutual in [1.0, 0.4] {
            let open = OpenDiffuseEnclosure::new(cx, surfaces(&mesh, &cards, &context),
                vec![vec![0.0, 0.4], vec![0.4, 0.0]],
                ViewFactorEvidence::Analytic { geometry: "synthetic reciprocal patch/bath algebra".into() },
                ViewFactorTolerance::default()).unwrap();
            let base = if mutual == 1.0 { PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &enclosure,
                PointwiseEnclosureTolerance::default()).unwrap() }
                else { PointwiseEnclosureRadiation::new_open(cx, &mesh, &boundary, &open, 300.0,
                    PointwiseEnclosureTolerance::default()).unwrap() };
            let owner = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
            assert_eq!(owner.is_closed(), mutual == 1.0);
            let actual = owner.assemble(cx, &h).unwrap(); let (load, heat, bath) = oracle(&curves, &h, mutual, 300.0);
            for i in 0..8 { close(actual.outward_load_w[i], load[i], 5e-10); }
            for (i, patch) in actual.report.patches.iter().enumerate() {
                close(patch.outward_heat_w, heat[i], 5e-10); assert_eq!(patch.material_states.len(), 3);
                close(patch.material_mean_emissivity.unwrap(), patch.material_states.iter()
                    .map(|state| state.emissivity.value()).sum::<f64>() / 3.0, 2e-15);
            }
            close(actual.report.outward_heat_w, heat.iter().sum(), 1e-9);
            close(actual.report.ambient.as_ref().map_or(0.0, |bath| bath.outward_heat_w), bath, 1e-9);
            close(actual.report.system_energy_closure_w, 0.0, 1e-9);
            let direction = [0.5, -0.7, 0.9, -0.2, -0.4, 1.1, 0.3, 0.6];
            let weights = [-0.3, 1.2, 0.7, 0.2, 0.4, -0.8, 0.5, -0.1];
            let mut action = [0.0; 8]; let mut transpose = [0.0; 8];
            actual.enthalpy_tangent_apply(cx, &direction, &mut action).unwrap();
            actual.enthalpy_tangent_apply_transpose(cx, &weights, &mut transpose).unwrap();
            close(action.iter().zip(weights).map(|(a, b)| a * b).sum(),
                transpose.iter().zip(direction).map(|(a, b)| a * b).sum(), 5e-12);
            let delta = 1e-3;
            let plus = owner.assemble(cx, &std::array::from_fn::<_, 8, _>(|i| h[i] + delta * direction[i])).unwrap();
            let minus = owner.assemble(cx, &std::array::from_fn::<_, 8, _>(|i| h[i] - delta * direction[i])).unwrap();
            for i in 0..8 { close(action[i], (plus.outward_load_w[i] - minus.outward_load_w[i]) / (2.0 * delta), 2e-7); }
        }
    });
}

#[test]
fn two_latent_materials_exchange_heat_through_changing_reflection_without_external_energy() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let cards = [material_card("left"), material_card("right")]; let enclosure = closed(&mesh, &cards, &context);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap(); let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material, source: &source, element_materials: None };
    let dt = 5.0; let driving = AREA * SIGMA * (400.0_f64.powi(4) - 350.0_f64.powi(4));
    let (mut lo, mut hi) = (0.0, driving / 3.0);
    for _ in 0..90 {
        let q = (lo + hi) / 2.0; let e_left = 0.5 - 0.6 * dt * q / (3.0 * LATENT);
        let e_right = 0.5 + 0.6 * dt * q / (3.0 * LATENT);
        if q - driving / (1.0 / e_left + 1.0 / e_right - 1.0) < 0.0 { lo = q; } else { hi = q; }
    }
    let heat = (lo + hi) / 2.0;
    for datum in [0.0, 4096.0] {
        let curves = [phase_curve(&cards[0], 400.0, datum), phase_curve(&cards[1], 350.0, datum)];
        let declarations = bindings(&cards, &curves, &context); let old = [datum + 500.0; 8];
        with_cx(|cx| {
            let base = PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &enclosure, PointwiseEnclosureTolerance::default()).unwrap();
            let owner = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
            let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
                &[ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: RHO },
                    ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: RHO }], &[0, 1], budget()).unwrap();
            let got = engine.advance_with_radiation(cx, problem, None, &old, dt, config(), EnthalpyRadiation::MaterialEnclosure(&owner)).unwrap();
            let next = got.conduction();
            for i in 0..3 {
                close(next.specific_enthalpy_j_kg[i] - datum, 500.0 - dt * heat / 3.0, 2e-8);
                close(next.specific_enthalpy_j_kg[i + 4] - datum, 500.0 + dt * heat / 3.0, 2e-8);
            }
            close(next.specific_enthalpy_j_kg[3], old[3], 2e-9); close(next.specific_enthalpy_j_kg[7], old[7], 2e-9);
            assert_eq!(&next.temperature[..4], &[400.0; 4]); assert_eq!(&next.temperature[4..], &[350.0; 4]);
            let report = next.pointwise_enclosure_radiation.as_ref().unwrap();
            close(report.patches[0].outward_heat_w, heat, 2e-8); close(report.patches[1].outward_heat_w, -heat, 2e-8);
            close(next.stored_energy_change_j, 0.0, 1e-8); close(next.energy_residual_j, 0.0, 1e-8);
            assert!(matches!(engine.linearize_accepted_with_pointwise_enclosure_radiation(cx, problem, None,
                &old, dt, config(), &base, next.clone(), 1e-8, 20000), Err(EnthalpyAdjointError::InvalidInput(_))));
        });
        assert_eq!(old, [datum + 500.0; 8]);
    }
}

#[test]
fn equal_temperature_latent_enclosure_has_exactly_zero_heat_despite_nonuniform_phase() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let cards = [material_card("left"), material_card("right")]; let enclosure = closed(&mesh, &cards, &context);
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 400.0, 0.0)];
    let declarations = bindings(&cards, &curves, &context);
    with_cx(|cx| {
        let base = PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &enclosure, PointwiseEnclosureTolerance::default()).unwrap();
        let owner = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
        let actual = owner.assemble(cx, &[100.0, 400.0, 900.0, 500.0, 200.0, 700.0, 800.0, 500.0]).unwrap();
        assert_eq!(actual.outward_load_w, vec![0.0; 8]); assert_eq!(actual.report.outward_heat_w, 0.0);
        let mut action = [0.0; 8]; actual.enthalpy_tangent_apply(cx, &[1.0; 8], &mut action).unwrap();
        assert_eq!(action, [0.0; 8]);
        assert!(owner.assemble(cx, &[-101.0; 8]).is_err());
    });
    with_cancelled_cx(|cx| assert!(PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &enclosure,
        PointwiseEnclosureTolerance::default()).is_err()));
}
