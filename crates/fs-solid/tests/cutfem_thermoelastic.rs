//! G1/G3: original thermal sources through native CutFEM and solid dynamics.
//! Authored curves verify implementation, not experimental calibration.

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::eigenstrain::PreparedEigenstrain3;
use fs_cutfem::elastic3::mass::MaterialMass3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3, MaterialRefreshLimits3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ivl::Interval;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::thermal::{
    AffineIsotropicThermoelastic, IsothermalThermoelasticInput, IsothermalThermoelasticState,
};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_material::tensor::{from_mandel, to_mandel};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solid::cutfront::dynamics3::{CutElasticityRef3, CutSecondOrderProblem3};
use fs_solver::LinearOp;
use fs_time::galpha::SecondOrderProblem;

const POINTS: usize = 100_000;
const TOTAL_STRAIN: f64 = 0.1;

fn poll() -> ControlFlow<()> { ControlFlow::Continue(()) }

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn control<T>(action: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 { depth: 0, ..Default::default() }, &mut callback,
    ).unwrap();
    action(&mut control)
}
fn domain() -> HexCell { HexCell::try_new([0.0; 3], [1.0; 3]).unwrap() }
fn origin(p: [f64; 3]) -> bool { p == [0.0; 3] }
fn options() -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation: ElasticityFormulation3::BBar,
        ghost_gamma: 0.3,
        ..Default::default()
    }
}
fn tree() -> Octree3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    tree.refined(&[*tree.leaves().iter().next().unwrap()], poll).unwrap()
}
fn cartesian(field: ElasticMaterialField3<'_>) -> CutElasticity3 {
    control(|control| CutElasticity3::build_with_material_field(
        domain(), [2; 3], &Top, field, &origin, options(), control,
    ).unwrap())
}
fn adaptive(field: ElasticMaterialField3<'_>, tree: &Octree3) -> AdaptiveElasticity3 {
    control(|control| AdaptiveElasticity3::build_with_material_field(
        domain(), tree, &Top, field, &origin, options(), control,
    ).unwrap())
}

fn law(strain_limit: f64) -> AffineIsotropicThermoelastic {
    let [poisson, young_slope, nu_slope, alpha, alpha_slope] =
        [0.2, -2.0, 0.0008, 0.0001, 0.000002];
    let temperature_quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0, young_slope),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), poisson, nu_slope),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            alpha,
            alpha_slope,
        ),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (280.0, value - 20.0 * slope),
                    (300.0, value),
                    (360.0, value + 60.0 * slope),
                ],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity, 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored CutFEM native thermal-state fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermoelastic source".into(),
            phase: "solid".into(),
            process: "numerical derivative verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("T", temperature_quantity, 300.0)
        .unwrap();
    let elastic = resolve_material_state_point(
        &card,
        &query,
        &[
            ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                POISSON_RATIO_PROPERTY, Dims([0; 6]), ScalarAdmissibility::Finite,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let expansion = resolve_material_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ScalarAdmissibility::Finite,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    AffineIsotropicThermoelastic::new(
        &card, &elastic, &expansion, [280.0, 360.0], strain_limit,
    )
    .unwrap()
}


fn frozen(law: &AffineIsotropicThermoelastic, temperature: f64) -> IsothermalThermoelasticState {
    law.at_temperature(temperature).unwrap()
        .isothermal_elasticity(IsothermalThermoelasticInput {
            // Explicit authored reference inventory, independent of T.
            reference_density_kg_m3: 3.0,
            mechanical_material_identity: ContentHash([82; 32]),
        }).unwrap()
}
fn temperature(bounds: HexCell) -> f64 {
    let [x, y, z] = bounds.lo();
    306.0 + 15.0 * x + 10.0 * y + 5.0 * z
}
fn affine(nodes: &[[f64; 3]], strain: f64) -> Vec<f64> {
    nodes.iter().flat_map(|&p| p.map(|value| strain * value)).collect()
}
fn field(nodes: &[[f64; 3]]) -> Vec<f64> {
    nodes.iter().flat_map(|&[x, y, z]| [
        0.006 * x * (1.0 + y),
        -0.003 * y * (1.0 + z),
        0.004 * z * (1.0 + x),
    ]).collect()
}
fn action(op: &dyn LinearOp, u: &[f64]) -> Vec<f64> {
    let mut force = vec![0.0; op.n()];
    op.apply(u, &mut force);
    force
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!((a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}");
}
fn close_vectors(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) { close(a, b, 3e-12, 3e-11); }
}
fn check_native_dynamics(
    source: CutElasticityRef3<'_>,
    mass: &MaterialMass3,
    prepared: &PreparedEigenstrain3,
    displacement: &[f64],
    cells: usize,
) {
    let occupancy = vec![1.0; cells];
    let gate = CancelGate::new_clock_free();
    let dynamics = CutSecondOrderProblem3::new_with_eigenstrain(
        source, mass, &occupancy, prepared, TOTAL_STRAIN, POINTS, 10_000, &gate,
    ).unwrap();
    let q: Vec<_> = dynamics.free_dofs().iter().map(|&dof| displacement[dof]).collect();
    let mut resistance = vec![0.0; q.len()];
    dynamics.internal_force(&q, &mut resistance);
    assert!(resistance.iter().all(|value| value.is_finite() && value.abs() < 3e-10));
}

#[test]
fn current_source_free_expansion_is_stress_free_on_both_bbar_spaces_and_native_dynamics() {
    let law = law(0.05);
    let current = frozen(&law, 330.0);
    let free = current.source().free_strain().unwrap();
    close(free, 0.0001 * 30.0 + 0.000001 * 30.0 * 30.0, 0.0, 2e-17);
    let op = cartesian(ElasticMaterialField3::ThermoelasticUniform(&current));
    let tree = tree();
    let adapted = adaptive(ElasticMaterialField3::ThermoelasticUniform(&current), &tree);
    assert!(adapted.physical_nodes().len() > adapted.nodes().len());
    let u = affine(op.nodes(), free);
    let ua = affine(adapted.nodes(), free);
    let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let prepared_a = adapted.prepare_thermoelastic(POINTS, poll).unwrap();
    close_vectors(&action(&op, &u), &op.eigenstrain_load(&prepared, poll).unwrap());
    close_vectors(&action(&adapted, &ua), &adapted.eigenstrain_load(&prepared_a, poll).unwrap());
    op.admit_thermoelastic(&u, TOTAL_STRAIN, POINTS, poll).unwrap();
    adapted.admit_thermoelastic(&ua, TOTAL_STRAIN, POINTS, poll).unwrap();
    assert_eq!(op.cell_thermoelastic_state(0).unwrap(), &current);
    assert_eq!(adapted.cell_thermoelastic_state(0).unwrap(), &current);
    for responses in [
        op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap(),
        adapted.thermoelastic_response(&ua, TOTAL_STRAIN, POINTS, poll).unwrap(),
    ] {
        assert!(!responses.is_empty());
        for response in responses {
            assert_eq!(response.temperature_k, 330.0);
            assert!(response.mechanical.stress.iter().all(|value| value.abs() < 2e-10));
            close(response.free_energy_j, 0.0, 0.0, 1e-18);
            close(response.entropy_j_k, 0.0, 0.0, 1e-16);
            close(response.internal_energy_j, 0.0, 0.0, 1e-13);
            // At fixed total strain, the elastic capacity contribution is
            // negative even at a stress-free state. It is not caloric Cp.
            assert!(response.heat_capacity_j_k < 0.0);
        }
    }
    let mass = op.prepare_material_mass(POINTS, poll).unwrap();
    let mass_a = adapted.prepare_material_mass(POINTS, poll).unwrap();
    check_native_dynamics(CutElasticityRef3::Cartesian(&op), &mass, &prepared, &u, op.cells());
    check_native_dynamics(CutElasticityRef3::Adaptive(&adapted), &mass_a, &prepared_a, &ua, adapted.cells());
}

#[test]
fn spatial_temperature_refresh_matches_fresh_operators_and_preserves_original_mass() {
    let law = law(0.05);
    let reference = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));

    let mut op = cartesian(ElasticMaterialField3::ThermoelasticUniform(&reference));
    let nodes = op.nodes().to_vec();
    let inventory = op.reference_cell_masses_kg().unwrap().to_vec();
    let scales: Vec<_> = (0..op.cells()).map(|cell| 0.55 + 0.02 * cell as f64).collect();
    op.set_scales(&scales).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let old_mass = op.prepare_material_mass(POINTS, poll).unwrap();
    let diagonal = old_mass.diagonal(&occupancy, poll).unwrap();
    let old_source = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let u = field(op.nodes());
    let cold_action = action(&op, &u);
    op.refresh_material_field(
        ElasticMaterialField3::ThermoelasticByCell(&hot),
        MaterialRefreshLimits3::default(), poll,
    ).unwrap();
    let mut fresh = cartesian(ElasticMaterialField3::ThermoelasticByCell(&hot));
    fresh.set_scales(&scales).unwrap();
    assert_eq!(op.nodes(), nodes.as_slice());
    assert_eq!(op.scales(), scales.as_slice());
    assert_eq!(op.reference_cell_masses_kg().unwrap(), inventory.as_slice());
    assert!(op.validate_material_mass(&old_mass).is_err());
    assert!(op.validate_eigenstrain(&old_source).is_err());
    let mass = op.prepare_material_mass(POINTS, poll).unwrap();
    assert_eq!(mass.diagonal(&occupancy, poll).unwrap(), diagonal);
    let image = action(&op, &u);
    close_vectors(&image, &action(&fresh, &u));
    assert!(image.iter().zip(&cold_action).any(|(a, b)| (a-b).abs() > 1e-4));
    let source = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let fresh_source = fresh.prepare_thermoelastic(POINTS, poll).unwrap();
    close_vectors(&op.eigenstrain_load(&source, poll).unwrap(),
        &fresh.eigenstrain_load(&fresh_source, poll).unwrap());
    let responses = op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap();
    assert!(responses.iter().any(|response| response.temperature_k != responses[0].temperature_k));
    assert!(responses.iter().any(|response| response.mechanical.stress.iter().any(|s| s.abs() > 1e-4)));
    for cell in 0..op.cells() {
        let current = op.cell_thermoelastic_state(cell).unwrap();
        assert_eq!(current, fresh.cell_thermoelastic_state(cell).unwrap());
        assert_eq!(current.input().reference_density_kg_m3, 3.0);
    }

    let tree = tree();
    let mut op = adaptive(ElasticMaterialField3::ThermoelasticUniform(&reference), &tree);
    let nodes = op.nodes().to_vec();
    let physical_nodes = op.physical_nodes().to_vec();
    let inventory = op.reference_cell_masses_kg().unwrap().to_vec();
    let scales: Vec<_> = (0..op.cells()).map(|cell| 0.55 + 0.01 * cell as f64).collect();
    op.set_scales(&scales).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let old_mass = op.prepare_material_mass(POINTS, poll).unwrap();
    let diagonal = old_mass.diagonal(&occupancy, poll).unwrap();
    let old_source = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let u = field(op.nodes());
    op.refresh_material_field(
        ElasticMaterialField3::ThermoelasticByCell(&hot),
        MaterialRefreshLimits3::default(), poll,
    ).unwrap();
    let mut fresh = adaptive(ElasticMaterialField3::ThermoelasticByCell(&hot), &tree);
    fresh.set_scales(&scales).unwrap();
    assert_eq!(op.nodes(), nodes.as_slice());
    assert_eq!(op.physical_nodes(), physical_nodes.as_slice());
    assert_eq!(op.scales(), scales.as_slice());
    assert_eq!(op.reference_cell_masses_kg().unwrap(), inventory.as_slice());
    assert!(op.validate_material_mass(&old_mass).is_err());
    assert!(op.validate_eigenstrain(&old_source).is_err());
    let mass = op.prepare_material_mass(POINTS, poll).unwrap();
    assert_eq!(mass.diagonal(&occupancy, poll).unwrap(), diagonal);
    close_vectors(&action(&op, &u), &action(&fresh, &u));
    let source = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let fresh_source = fresh.prepare_thermoelastic(POINTS, poll).unwrap();
    close_vectors(&op.eigenstrain_load(&source, poll).unwrap(),
        &fresh.eigenstrain_load(&fresh_source, poll).unwrap());
    for cell in 0..op.cells() {
        assert_eq!(op.cell_thermoelastic_state(cell).unwrap(), fresh.cell_thermoelastic_state(cell).unwrap());
    }
}

#[test]
fn native_bulk_potentials_and_mixed_stress_derivatives_follow_fixed_displacement_temperature_changes() {
    let law = law(0.05);
    let t = 321.0;
    let h = 0.02;
    let current = frozen(&law, t);
    let mut op = cartesian(ElasticMaterialField3::ThermoelasticUniform(&current));
    let scales: Vec<_> = (0..op.cells()).map(|cell| 0.6 + 0.04 * cell as f64).collect();
    op.set_scales(&scales).unwrap();
    let u = field(op.nodes());
    let mut samples = Vec::new();
    for offset in [-2.0, -1.0, 0.0, 1.0, 2.0] {
        let state = frozen(&law, t + offset * h);
        op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticUniform(&state),
            MaterialRefreshLimits3::default(), poll,
        ).unwrap();
        samples.push(op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll).unwrap());
    }
    let first = |values: [f64; 5]| {
        (values[0] - 8.0 * values[1] + 8.0 * values[3] - values[4]) / (12.0 * h)
    };
    let second = |values: [f64; 5]| {
        (-values[0] + 16.0 * values[1] - 30.0 * values[2]
            + 16.0 * values[3] - values[4]) / (12.0 * h * h)
    };
    for point in 0..samples[2].len() {
        let actual = &samples[2][point];
        let scale = scales[actual.mechanical.cell];
        let projected_total = std::array::from_fn(|i| {
            actual.mechanical.constitutive_elastic_strain[i] + actual.mechanical.free_strain[i]
        });
        let native = current.source().state(
            to_mandel(&projected_total), actual.mechanical.weight,
        ).unwrap();
        for (got, expected) in [
            (actual.free_energy_j, native.free_energy_j),
            (actual.entropy_j_k, native.entropy_j_k),
            (actual.internal_energy_j, native.internal_energy_j),
            (actual.heat_capacity_j_k, native.heat_capacity_j_k),
        ] {
            close(got, scale * expected, 3e-11, 2e-13);
        }
        let beta = from_mandel(&native.stress_temperature_derivative_mandel_pa_k);
        for i in 0..6 {
            close(actual.stress_temperature_derivative[i], scale * beta[i], 3e-11, 2e-11);
        }
        for sample in &samples {
            assert_eq!(sample.len(), samples[2].len());
            assert_eq!(sample[point].mechanical.cell, actual.mechanical.cell);
            assert_eq!(sample[point].mechanical.position, actual.mechanical.position);
            assert_eq!(sample[point].mechanical.weight, actual.mechanical.weight);
        }
        let psi = std::array::from_fn(|i| samples[i][point].free_energy_j);
        let internal = std::array::from_fn(|i| samples[i][point].internal_energy_j);
        close(actual.entropy_j_k, -first(psi), 2e-7, 2e-10);
        close(actual.heat_capacity_j_k, first(internal), 2e-7, 2e-8);
        close(actual.heat_capacity_j_k, -t * second(psi), 2e-5, 2e-7);
        close(actual.internal_energy_j,
            actual.free_energy_j + t * actual.entropy_j_k, 2e-12, 2e-13);
        for component in 0..6 {
            let stress = std::array::from_fn(|i| samples[i][point].mechanical.stress[component]);
            close(actual.stress_temperature_derivative[component], first(stress), 2e-7, 2e-8);
        }
    }
}

#[test]
fn original_total_component_domain_source_ownership_and_cancellation_remain_mandatory() {
    let law = law(0.01);
    let current = frozen(&law, 350.0);
    let free = current.source().free_strain().unwrap();
    let elastic = [0.012-free, 0.012-free, 0.012-free, 0.0, 0.0, 0.0];
    assert!(current.acoustic().admits_mandel_strain(&elastic));
    let op = cartesian(ElasticMaterialField3::ThermoelasticUniform(&current));
    let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let total_outside = affine(op.nodes(), 0.012);
    // The separate Frobenius budget is deliberately generous. Total native
    // components exceed 0.01 even though every elastic component is inside.
    assert!(op.admit_thermoelastic(&total_outside, TOTAL_STRAIN, POINTS, poll).is_err());
    assert!(op.thermoelastic_response(&total_outside, TOTAL_STRAIN, POINTS, poll).is_err());
    assert!(op.eigenstrain_response(&prepared, &total_outside, TOTAL_STRAIN, POINTS, poll).is_err());
    let tree = tree();
    let adapted = adaptive(ElasticMaterialField3::ThermoelasticUniform(&current), &tree);
    let total_outside = affine(adapted.nodes(), 0.012);
    assert!(adapted.admit_thermoelastic(&total_outside, TOTAL_STRAIN, POINTS, poll).is_err());
    assert!(adapted.thermoelastic_response(&total_outside, TOTAL_STRAIN, POINTS, poll).is_err());

    let valid = affine(op.nodes(), free);
    op.admit_thermoelastic(&valid, TOTAL_STRAIN, POINTS, poll).unwrap();
    let force_before = op.eigenstrain_load(&prepared, poll).unwrap();
    let mut calls = 0usize;
    op.thermoelastic_response(&valid, TOTAL_STRAIN, POINTS,
        || { calls += 1; poll() }).unwrap();
    assert!(calls > 4);
    for stop in [1, calls / 2, calls] {
        let mut count = 0usize;
        assert!(matches!(op.thermoelastic_response(&valid, TOTAL_STRAIN, POINTS,
            || {
                count += 1;
                if count == stop { ControlFlow::Break(()) } else { poll() }
            }), Err(ElasticityError3::Cancelled)));
    }
    assert!(op.prepare_thermoelastic(0, poll).is_err());
    assert!(op.admit_thermoelastic(&valid, TOTAL_STRAIN, 0, poll).is_err());
    op.validate_eigenstrain(&prepared).unwrap();
    assert_eq!(op.eigenstrain_load(&prepared, poll).unwrap(), force_before);

    // Equal numerical C and density cannot replace retained thermal sources.
    let plain = cartesian(ElasticMaterialField3::Uniform(current.acoustic()));
    assert!(plain.cell_thermoelastic_state(0).is_none());
    assert!(plain.prepare_thermoelastic(POINTS, poll).is_err());
    assert!(plain.thermoelastic_response(&valid, TOTAL_STRAIN, POINTS, poll).is_err());
}
