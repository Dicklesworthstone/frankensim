#![allow(dead_code)]
use fs_em::{
    planar::PlanarMagnetostatics,
    winding::{PlanarWindingBasis, WindingConfig, circuit::WindingCircuit,
        ohmic::{OhmicWindingControls, WindingConductor}},
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY},
    state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
        ScalarAdmissibility, ScalarPropertyRequirement, SPECIFIC_HEAT_CAPACITY_DIMS,
        resolve_material_state_point},
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, QuantitySpec, Temperature, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;

fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("T",
        QuantitySpec::dimensional(Temperature::DIMS), 200.0, 2000.0)
}
fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}
fn claim(key: PropertyKey, value: PropertyValue) -> PropertyClaim {
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else { InterpolationPolicy::LinearInside };
    let validity = if matches!(&value, PropertyValue::Curve { abscissa, .. } if abscissa == "pressure") {
        validity().with_quantity("pressure", pressure_quantity(), 1e5, 2e5)
    } else { validity() };
    PropertyClaim { key, value, validity, interpolation,
        uncertainty: UncertaintyModel::Unstated, observations: vec![],
        provenance: Provenance { source: "synthetic winding material regression".into(),
            license: "test-only".into(), artifact: None } }
}
fn assemble(rho: f64, alpha: f64, mu: PropertyValue, basis: ElectricalPropertyBasis) -> MaterialCard {
    let electrical = |temperature| {
        let value = rho * (1.0 + alpha * (temperature - 300.0));
        match basis { ElectricalPropertyBasis::Resistivity => value,
            ElectricalPropertyBasis::Conductivity => 1.0 / value }
    };
    let electrical = if alpha == 0.0 {
        PropertyValue::Scalar { value: electrical(300.0), dims: basis.dimensions() }
    } else {
        PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: [200.0, 300.0, 2000.0].map(|t| (t, electrical(t))).to_vec(), dims: basis.dimensions() }
    };
    assemble_values(electrical, mu, basis)
}
fn assemble_values(electrical: PropertyValue, mu: PropertyValue, basis: ElectricalPropertyBasis) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for item in [
        claim(PropertyKey::new(basis.property_name(), basis.dimensions()), electrical),
        claim(PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS), mu),
        claim(PropertyKey::new("density", Density::DIMS),
            PropertyValue::Scalar { value: 1.0, dims: Density::DIMS }),
        claim(PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            PropertyValue::Scalar { value: 4.0, dims: SPECIFIC_HEAT_CAPACITY_DIMS }),
    ] { claims.insert_claim(item).unwrap(); }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic winding".into(),
        phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap()
}
pub fn card(rho: f64, alpha: f64, mu: f64, basis: ElectricalPropertyBasis) -> MaterialCard {
    assemble(rho, alpha, PropertyValue::Scalar { value: mu, dims: MAGNETIC_PERMEABILITY_DIMS }, basis)
}
pub fn magnetic_curve_card(rho: f64, mu: f64, basis: ElectricalPropertyBasis) -> MaterialCard {
    // Same endpoint mu still cannot hide work along a temperature-dependent law.
    assemble(rho, 0.001, PropertyValue::Curve {
        abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(200.0, mu), (300.0, mu), (600.0, 2.0 * mu), (2000.0, mu)],
        dims: MAGNETIC_PERMEABILITY_DIMS,
    }, basis)
}
pub fn affine_magnetic_card(rho:f64, alpha:f64, mu:f64, slope:f64,
    basis:ElectricalPropertyBasis)->MaterialCard {
    assemble(rho,alpha,PropertyValue::Curve {
        abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
        knots:[200.0,300.0,2000.0].map(|t|(t,mu+slope*(t-300.0))).to_vec(),
        dims:MAGNETIC_PERMEABILITY_DIMS,
    },basis)
}
/// Explicit absolute-temperature source domain for thermodynamic entropy tests.
/// Dimensional-kelvin transport fixtures remain separately authored above.
pub fn absolute_affine_magnetic_card(rho:f64, alpha:f64, mu:f64, slope:f64,
    basis:ElectricalPropertyBasis)->MaterialCard {
    let original=affine_magnetic_card(rho,alpha,mu,slope,basis);
    let mut claims=ClaimSet::new();
    for (_,source) in original.claims().claims_ordered() {
        let mut claim=source.clone();
        claim.validity=claim.validity.with_quantity("T",QuantitySpec::semantic(
            SemanticType::new(QuantityKind::AbsoluteTemperature,ValueForm::Static)),200.0,2000.0);
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(),claims,vec![]).unwrap()
}
pub fn absolute_point(t:f64)->QueryPoint {
    QueryPoint::new().with_quantity("T",QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,ValueForm::Static)),t).unwrap()
}
pub fn pressure_magnetic_card(rho: f64, mu: f64, basis: ElectricalPropertyBasis) -> MaterialCard {
    assemble(rho, 0.001, PropertyValue::Curve {
        abscissa: "pressure".into(), abscissa_dims: fs_qty::Dims([-1, 1, -2, 0, 0, 0]),
        knots: vec![(1e5, mu), (2e5, 2.0 * mu)], dims: MAGNETIC_PERMEABILITY_DIMS,
    }, basis)
}
pub fn pressure_electrical_magnetic_card(rho: f64, mu: f64, basis: ElectricalPropertyBasis) -> MaterialCard {
    let value = |rho| match basis {
        ElectricalPropertyBasis::Resistivity => rho,
        ElectricalPropertyBasis::Conductivity => 1.0 / rho,
    };
    let pressure_dims = fs_qty::Dims([-1, 1, -2, 0, 0, 0]);
    assemble_values(
        PropertyValue::Curve {
            abscissa: "pressure".into(), abscissa_dims: pressure_dims,
            knots: vec![(1e5, value(rho)), (2e5, value(2.0 * rho))],
            dims: basis.dimensions(),
        },
        PropertyValue::Curve {
            abscissa: "pressure".into(), abscissa_dims: pressure_dims,
            knots: vec![(1e5, mu), (2e5, 2.0 * mu)], dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        basis,
    )
}
pub fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
}
pub fn pressure_point(t: f64, pressure: f64) -> QueryPoint {
    point(t).with_quantity("pressure", pressure_quantity(), pressure).unwrap()
}
pub fn bundle(card: &MaterialCard, t: f64, basis: ElectricalPropertyBasis) -> ResolvedMaterialStatePoint {
    bundle_at(card, &point(t), basis)
}
pub fn bundle_at(card: &MaterialCard, point: &QueryPoint, basis: ElectricalPropertyBasis) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(card, point, &[
        basis.requirement().unwrap(),
        ScalarPropertyRequirement::try_new(MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_new("density", Density::DIMS,
            ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS,
            ScalarAdmissibility::StrictlyPositive).unwrap(),
    ], MaterialPropertySelection::SingleClaimOnly).unwrap()
}
pub fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits { max_cells: 1000, max_dofs: 1000, max_triplets: 10000 }
}
pub fn linear() -> CheckedCgConfig {
    CheckedCgConfig { tolerance: 1e-11, max_iterations: 1000, max_corrections: 2 }
}
pub fn winding_config() -> WindingConfig {
    WindingConfig { max_windings: 4, max_basis_values: 10000,
        pairing_tolerance_h: 1e-10, current_balance_tolerance_turns: 1e-9,
        inverse_residual_tolerance: 1e-10 }
}
pub fn problem(cards: &[MaterialCard; 2], basis: ElectricalPropertyBasis, depth_m: f64)
    -> (PlanarMagnetostatics, Vec<Vec<f64>>, Vec<usize>)
{
    problem_at(cards, basis, depth_m, &point(300.0))
}
pub fn problem_at(cards: &[MaterialCard; 2], basis: ElectricalPropertyBasis, depth_m: f64, query: &QueryPoint)
    -> (PlanarMagnetostatics, Vec<Vec<f64>>, Vec<usize>)
{
    let n = 8_usize;
    let mut vertices = Vec::new();
    let mut boundary = Vec::new();
    for y in 0..=n { for x in 0..=n {
        if x == 0 || y == 0 || x == n || y == n { boundary.push(vertices.len()); }
        vertices.push([x as f64 / n as f64, y as f64 / n as f64]);
    } }
    let mut faces = Vec::new();
    for y in 0..n { for x in 0..n {
        let a = (y * (n + 1) + x) as u32;
        let b = a + 1;
        let d = a + (n + 1) as u32;
        let c = d + 1;
        faces.extend([[a, b, c], [a, c, d]]);
    } }
    let mesh = TriComplex2::from_indexed_triangles("material-stranded-winding", vertices,
        faces, Metric2::planar(depth_m).unwrap()).unwrap();
    let mut profiles = vec![Vec::new(), Vec::new()];
    let mut materials = Vec::new();
    let resolved = cards.each_ref().map(|c| bundle_at(c, query, basis));
    for face in mesh.faces() {
        let centre: [f64; 2] = [0, 1].map(|axis| face.iter()
            .map(|&v| mesh.vertices()[v as usize][axis] / 3.0).sum());
        let owner = usize::from(centre[1] >= 0.5);
        // Each disjoint half has 100 outgoing + 100 return turns, and volume
        // 0.25 m3 at depth 0.5 m. Magnetic mutual coupling crosses the split.
        let j = if centre[0] < 0.5 { 400.0 } else { -400.0 };
        for (winding, profile) in profiles.iter_mut().enumerate() {
            profile.push(if winding == owner { j } else { 0.0 });
        }
        materials.push(resolved[owner].clone());
    }
    (PlanarMagnetostatics::new(mesh, materials, limits(), &mut || false).unwrap(), profiles, boundary)
}
pub fn graph() -> CircuitGraph {
    CircuitGraph { node_count: 4, branches: vec![
        (3, 0, Branch::Capacitor { farads: 0.4 }),
        (1, 0, Branch::Resistor { ohms: 1.0 }),
        (2, 0, Branch::Inductor { henries: 7.0 }),
        (3, 0, Branch::Resistor { ohms: 2.0 }),
        (1, 0, Branch::Inductor { henries: 11.0 }),
        (2, 0, Branch::Resistor { ohms: 3.0 }),
        (1, 0, Branch::CurrentSource { port: 0 }),
    ], transformers: vec![] }
}
pub fn winding<'a>(problem: &'a PlanarMagnetostatics, profiles: Vec<Vec<f64>>, boundary: &[usize]) -> WindingCircuit<'a> {
    PlanarWindingBasis::new(problem, profiles, boundary, vec![0.002, 0.003],
        winding_config(), linear(), &mut || false).unwrap()
        .into_circuit(graph(), vec![4, 2], &mut || false).unwrap()
}
pub fn sources(magnetic: &WindingCircuit<'_>, basis: ElectricalPropertyBasis) -> Vec<WindingConductor> {
    magnetic.basis().problem().materials().iter().enumerate().filter(|(triangle, _)| {
        magnetic.basis().profiles_per_m2().iter().any(|profile| profile[*triangle] != 0.0)
    }).map(|(triangle, material)| WindingConductor { triangle,
        conductor: ResolvedConductorState::try_from_material_with_basis(material.clone(), basis).unwrap() }).collect()
}
pub fn controls() -> OhmicWindingControls {
    OhmicWindingControls { power_tolerance_w: 1e-8, energy_tolerance_j: 1e-8 }
}
pub fn two_by_two(a: [f64; 4], b: [f64; 2]) -> [f64; 2] {
    let det = a[0] * a[3] - a[1] * a[2];
    [(b[0] * a[3] - a[1] * b[1]) / det, (a[0] * b[1] - b[0] * a[2]) / det]
}
pub fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-9 * b.abs().max(1.0), "{a} != {b}");
}
pub fn refreshed(circuit: &fs_em::winding::ohmic::OhmicWindingCircuit<'_, '_>,
    cards: &[MaterialCard; 2], temperatures: [f64; 2]) -> Vec<WindingConductor>
{
    circuit.conductors().iter().zip(circuit.conductor_windings()).map(|(source, &winding)| {
        let snapshot = source.conductor.material().at_query_point(&cards[winding], &point(temperatures[winding])).unwrap();
        WindingConductor { triangle: source.triangle,
            conductor: ResolvedConductorState::try_from_material_with_basis(snapshot, source.conductor.source_basis()).unwrap() }
    }).collect()
}

