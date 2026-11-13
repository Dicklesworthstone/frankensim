//! G0/G1/G3 source data -> temperature refresh -> actual rod equilibrium.
//! Synthetic constitutive curves test wiring; they do not qualify a real alloy.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::LayerIsotropicParameters;
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::thermal::source::{RodLayerSource, RodMaterialSource, RodSegmentSource};

// Pressure is a fixed external condition only as a typed static pressure. Every rod layer
// property must be constant at fixed environment, and that rule refuses an untyped extra axis
// as a possibly evolving coordinate (only the legacy untyped "T" is exempt).
fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}
fn claim(name: &str, dims: Dims, knots: Vec<(f64, f64)>) -> PropertyClaim {
    let lo = knots[0].0;
    let hi = knots.last().unwrap().0;
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS, knots, dims },
        validity: ValidityDomain::unconstrained().with("T", lo, hi)
            .with_quantity("p", pressure_quantity(), 100_000.0, 100_000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside, observations: vec![],
        provenance: Provenance { source: "synthetic heated rod".into(), license: "test-only".into(), artifact: None },
    }
}
fn card(alpha_scale: f64, split_modulus: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, a, b) in [
        ("density", Density::DIMS, 1000.0, 800.0),
        ("young_modulus", Pressure::DIMS, 1e6, 5e5),
        ("shear_modulus", Pressure::DIMS, 4e5, 2e5),
        (ALPHA, INVERSE_TEMPERATURE_DIMS, alpha_scale*1e-5, alpha_scale*3e-5),
    ] {
        if split_modulus && name == "young_modulus" {
            claims.insert_claim(claim(name, dims, vec![(300.0, a), (350.0, 0.9*a)])).unwrap();
            claims.insert_claim(claim(name, dims, vec![(400.0, b), (500.0, b)])).unwrap();
        } else {
            claims.insert_claim(claim(name, dims, vec![(300.0, a), (500.0, b)])).unwrap();
        }
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
        .with_quantity("p", pressure_quantity(), 100_000.0).unwrap()
}
fn segment(card: &MaterialCard, radius: f64) -> RodSegmentSource<'_> {
    RodSegmentSource { inner_radius_m: 0.0, shear_correction: 0.8,
        layers: vec![RodLayerSource { outer_radius_m: radius, card, reference_point: point(),
            mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
            expansion_selection: MaterialPropertySelection::SingleClaimOnly }] }
}
fn source(card: &MaterialCard) -> RodMaterialSource<'_> {
    RodMaterialSource::try_new(0.2, vec![segment(card, 0.01)], "shear_modulus", 0.03,
        &CancelGate::new()).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1e-12), "{a} != {b}");
}

#[test]
fn nonuniform_thermal_refresh_preserves_inventory_and_actual_series_expansion() {
    use fs_solid::rod::material_field::motion::RodRootBoundary;
    let material = card(1.0, false);
    let lengths = vec![0.125, 0.375];
    let radii = [0.005_f64, 0.01];
    let gate = CancelGate::new();
    let source = RodMaterialSource::try_new_nonuniform(lengths.clone(),
        radii.iter().map(|&radius| segment(&material, radius)).collect(),
        "shear_modulus", 0.03, &gate).unwrap();
    let original = source.reference_rod();
    let model = original.clone().into_motion_model(RodRootBoundary::Clamped, &gate).unwrap();
    let temperatures = [350.0, 450.0];
    let mut hot = source.resolve_at(&temperatures, &gate).unwrap();
    assert_eq!(hot.reference_segment_lengths_m(), lengths);
    assert_eq!(hot.reference_mass(), original.reference_mass());
    let mut expected_extension = 0.0;
    let mut expected_mass = 0.0;
    for i in 0..2 {
        let delta = temperatures[i] - 300.0;
        let area = core::f64::consts::PI * radii[i].powi(2);
        let young = 1e6 - 2500.0 * delta;
        let eigenstrain = 1e-5 * delta + 0.5e-7 * delta * delta;
        expected_extension += lengths[i] * (eigenstrain + 0.1 / (young * area));
        expected_mass += lengths[i] * 1000.0 * area;
        near(hot.thermal_sections().unwrap()[i].zero_resultant_axial_strain(), eigenstrain, 1e-12);
    }
    hot.solve_static(&TipLoad { force: [0.1, 0.0, 0.0], moment: [0.0; 3] },
        2, 1e-10, &gate).unwrap();
    near(hot.positions_m()[2][0] - 0.5, expected_extension, 2e-7);
    near(source.reference_mass_kg(), expected_mass, 2e-15);
    let refreshed = model.rebind_material_temperatures(&source, &temperatures, &gate).unwrap();
    assert_eq!(refreshed.inertia(), model.inertia());
    assert_eq!(refreshed.source().reference_mass(), model.source().reference_mass());
    assert_eq!(refreshed.source().reference_segment_lengths_m(), lengths);
    let cold = source.resolve_at(&[300.0; 2], &gate).unwrap();
    assert_eq!(cold.energy_j().unwrap(), 0.0);
}

#[test]
fn nonuniform_isotropic_sources_retain_interval_identity_across_reference_geometry() {
    use fs_ga::So3;
    use fs_solid::rod::material_field::reference::RodReferenceGeometry;
    let material = bulk_shear_card();
    let lengths = vec![0.125, 0.375];
    let gate = CancelGate::new();
    let bindings = vec![vec![bulk_shear_binding()]; 2];
    let source = RodMaterialSource::try_new_nonuniform_with_isotropic_parameters(lengths.clone(),
        vec![segment(&material, 0.01); 2], &bindings, 0.03, &gate).unwrap();
    let geometry = RodReferenceGeometry::try_nonuniform(lengths.clone(),
        vec![[0.0; 3], [0.125, 0.0, 0.0], [0.5, 0.0, 0.0]], vec![So3::identity(); 3]).unwrap();
    let bound = source.with_reference_geometry(geometry, &gate).unwrap();
    let hot = bound.resolve_at(&[350.0, 450.0], &gate).unwrap();
    assert_eq!(hot.reference_geometry().unwrap().segment_lengths_m(), lengths);
    assert_eq!(hot.reference_mass(), source.reference_rod().reference_mass());
    let wrong_partition = RodReferenceGeometry::try_new(0.5,
        vec![[0.0; 3], [0.25, 0.0, 0.0], [0.5, 0.0, 0.0]], vec![So3::identity(); 3]).unwrap();
    assert!(source.with_reference_geometry(wrong_partition, &gate).is_err());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(source.resolve_at(&[350.0, 450.0], &cancelled).is_err());
}

#[test]
fn temperature_refresh_changes_stiffness_expansion_and_solved_extension_together() {
    let card = card(1.0, false);
    let source = source(&card);
    let gate = CancelGate::new();
    let original_mass = source.reference_mass_kg();
    let mut hot = source.resolve_at(&[400.0], &gate).unwrap();
    let area = std::f64::consts::PI*0.01_f64.powi(2);
    near(hot.sections()[0].elastic().axial_stiffness_n(), 750_000.0*area, 1e-14);
    near(hot.sections()[0].shear_stiffness_n(), 0.8*300_000.0*area, 1e-14);
    let thermal = &hot.thermal_sections().unwrap()[0];
    near(thermal.zero_resultant_axial_strain(), 0.0015, 1e-13);
    let load = TipLoad { force: [0.1, 0.0, 0.0], moment: [0.0; 3] };
    hot.solve_static(&load, 1, 1e-8, &gate).unwrap();
    let expected = 0.2*(0.0015 + 0.1/(750_000.0*area));
    near(hot.positions_m()[1][0]-0.2, expected, 2e-5);
    near(source.reference_mass_kg(), 1000.0*area*0.2, 1e-14);
    assert_eq!(source.reference_mass_kg(), original_mass);
    let current = &hot.sections()[0].elastic().layers()[0].material;
    assert_eq!(current.query_point(), &[("T".into(), 400.0), ("p".into(), 100_000.0)]);
    near(current.property("density").unwrap().value_si(), 900.0, 1e-14);
}

#[test]
fn paths_return_to_original_reference_without_accumulated_thermal_strain() {
    let card = card(1.0, false);
    let source = source(&card);
    let gate = CancelGate::new();
    for t in [400.0, 350.0, 450.0, 300.0, 400.0] {
        let rod = source.resolve_at(&[t], &gate).unwrap();
        let dt = t-300.0;
        let expected = 1e-5*dt + 0.5e-7*dt*dt;
        near(rod.thermal_sections().unwrap()[0].zero_resultant_axial_strain(), expected, 1e-12);
    }
    assert_eq!(source.reference_rod().positions_m(), &[[0.0; 3], [0.2, 0.0, 0.0]]);
}

#[test]
fn reference_pins_prevent_switching_to_a_different_fit_during_heating() {
    let card = card(1.0, true);
    let source = source(&card);
    assert!(source.resolve_at(&[325.0], &CancelGate::new()).is_ok());
    // An in-domain alternative exists at 425 K, but is not the reference claim.
    assert!(source.resolve_at(&[425.0], &CancelGate::new()).is_err());
}

#[test]
fn per_layer_and_per_segment_temperatures_keep_sources_and_spatial_order() {
    let a = card(1.0, false);
    let b = card(2.0, false);
    let mut first = segment(&a, 0.005);
    first.layers.push(RodLayerSource { outer_radius_m: 0.01, card: &b, reference_point: point(),
        mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
        expansion_selection: MaterialPropertySelection::SingleClaimOnly });
    let source = RodMaterialSource::try_new(0.4, vec![first, segment(&a, 0.01)],
        "shear_modulus", 0.03, &CancelGate::new()).unwrap();
    let rod = source.resolve_at(&[400.0, 350.0, 450.0], &CancelGate::new()).unwrap();
    let sections = rod.thermal_sections().unwrap();
    near(sections[0].expansion_paths()[0].free_linear_strain(), 0.0015, 1e-13);
    near(sections[0].expansion_paths()[1].free_linear_strain(), 0.00125, 1e-13);
    near(sections[1].expansion_paths()[0].free_linear_strain(), 0.002625, 1e-13);
    assert_eq!(sections[0].section().layers()[1].material.card_identity(), b.content_hash());
    assert_eq!(source.reference_temperatures_k(), &[300.0; 3]);
}

#[test]
fn incomplete_states_invalid_geometry_and_cancellation_publish_no_candidate() {
    let card = card(1.0, false);
    let source = source(&card);
    for values in [vec![], vec![300.0, 300.0], vec![f64::NAN], vec![0.0], vec![600.0]] {
        assert!(source.resolve_at(&values, &CancelGate::new()).is_err());
    }
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(source.resolve_at(&[350.0], &cancelled).is_err());
    let mut invalid = segment(&card, 0.01);
    invalid.layers[0].reference_point = QueryPoint::new().with("T", 300.0).unwrap();
    assert!(RodMaterialSource::try_new(0.2, vec![invalid], "shear_modulus", 0.03, &CancelGate::new()).is_err());
    assert!(RodMaterialSource::try_new(0.2, vec![segment(&card, 0.0)], "shear_modulus", 0.03, &CancelGate::new()).is_err());
    assert!(RodMaterialSource::try_new(0.2, vec![segment(&card, 0.01)], "young_modulus", 0.03, &CancelGate::new()).is_err());
    assert_eq!(source.reference_temperatures_k(), &[300.0]);
}

#[test]
fn curved_reference_survives_material_temperature_refresh_without_mass_or_stress_reset() {
    use fs_ga::{So3, So3Tangent, Vec3};
    use fs_math::det;
    use fs_solid::rod::material_field::reference::RodReferenceGeometry;
    let card = card(1.0, false);
    let original = source(&card);
    let gate = CancelGate::new();
    let angle = 0.2;
    let nodes = vec![[0.0; 3], [det::sin(angle), 1.0-det::cos(angle), 0.0]];
    let frames = vec![So3::identity(), So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, angle))).unwrap()];
    let source = original.with_reference_geometry(
        RodReferenceGeometry::try_new(0.2, nodes.clone(), frames.clone()).unwrap(), &gate).unwrap();
    assert_eq!(source.reference_rod().energy_j().unwrap(), 0.0);
    for temperature in [300.0, 400.0, 300.0] {
        let rod = source.resolve_at(&[temperature], &gate).unwrap()
            .with_gravity([0.0, -9.81, 0.0], &gate).unwrap();
        assert_eq!(rod.positions_m(), nodes);
        assert_eq!(rod.frames(), frames);
        assert_eq!(rod.segment_strains(0).unwrap(), ([0.0; 3], [0.0; 3]));
        assert_eq!(rod.reference_mass().total_mass_kg(), original.reference_mass_kg());
        assert_eq!(rod.gravity_potential_j().unwrap(), 0.0);
        let expansion = rod.thermal_sections().unwrap()[0].zero_resultant_axial_strain();
        near(rod.segment_resultants(0).unwrap().0[0],
            -rod.sections()[0].elastic().axial_stiffness_n() * expansion, 1e-12);
        assert_eq!(rod.segment_resultants(0).unwrap().1, [0.0; 3]);
    }
    assert_eq!(original.reference_rod().positions_m()[1], [0.2, 0.0, 0.0]);
}

fn scalar_claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    let mut property = claim(name, dims, vec![(300.0, value), (500.0, value)]);
    property.value = PropertyValue::Scalar { value, dims };
    property.interpolation = InterpolationPolicy::ConstantWithinValidity;
    property
}

fn material_from_claims(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic independent elastic sources".into(), phase: "solid".into(),
        process: "test".into(), revision: 0,
    }, claims, vec![]).unwrap()
}

fn positive_requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap()
}

fn bulk_shear_binding() -> LayerIsotropicParameters {
    LayerIsotropicParameters {
        density: positive_requirement("specimen_density", Density::DIMS),
        parameters: IsotropicElasticParameters::BulkShear {
            bulk_modulus: positive_requirement("bulk_response", Pressure::DIMS),
            shear_modulus: positive_requirement("shear_response", Pressure::DIMS),
        },
    }
}

fn bulk_shear_card() -> MaterialCard {
    material_from_claims(vec![
        claim("specimen_density", Density::DIMS, vec![(300.0, 1000.0), (500.0, 800.0)]),
        claim("bulk_response", Pressure::DIMS, vec![(300.0, 8e5), (500.0, 6e5)]),
        claim("shear_response", Pressure::DIMS, vec![(300.0, 4e5), (500.0, 2e5)]),
        claim(ALPHA, INVERSE_TEMPERATURE_DIMS, vec![(300.0, 1e-5), (500.0, 3e-5)]),
    ])
}

fn young_poisson_binding() -> LayerIsotropicParameters {
    LayerIsotropicParameters {
        density: positive_requirement("mass_density", Density::DIMS),
        parameters: IsotropicElasticParameters::YoungPoisson {
            young_modulus: positive_requirement("tensile_response", Pressure::DIMS),
            poisson_ratio: ScalarPropertyRequirement::try_new("transverse_ratio", Dims::NONE,
                ScalarAdmissibility::OpenInterval { lower: -1.0, upper: 0.5 }).unwrap(),
        },
    }
}

fn young_poisson_card() -> MaterialCard {
    material_from_claims(vec![
        scalar_claim("mass_density", Density::DIMS, 1200.0),
        claim("tensile_response", Pressure::DIMS, vec![(300.0, 1.2e6), (500.0, 6e5)]),
        scalar_claim("transverse_ratio", Dims::NONE, 0.25),
        scalar_claim(ALPHA, INVERSE_TEMPERATURE_DIMS, 2e-5),
    ])
}

#[test]
fn bulk_shear_material_drives_hot_rod_extension_and_torsion_without_young_claim() {
    let card = bulk_shear_card();
    let gate = CancelGate::new();
    let bindings = vec![vec![bulk_shear_binding()]];
    let source = RodMaterialSource::try_new_with_isotropic_parameters(
        0.2, vec![segment(&card, 0.01)], &bindings, 0.03, &gate).unwrap();
    let mut hot = source.resolve_at(&[400.0], &gate).unwrap();
    let area = std::f64::consts::PI * 0.01_f64.powi(2);
    let inertia = std::f64::consts::PI * 0.01_f64.powi(4) / 4.0;
    let young = 9.0 * 700_000.0 * 300_000.0 / (3.0 * 700_000.0 + 300_000.0);
    let section = &hot.sections()[0];
    near(section.elastic().axial_stiffness_n(), young * area, 1e-14);
    near(section.elastic().bending_stiffness_n_m2(), young * inertia, 1e-14);
    near(section.shear_stiffness_n(), 0.8 * 300_000.0 * area, 1e-14);
    near(section.torsional_stiffness_n_m2(), 300_000.0 * 2.0 * inertia, 1e-14);
    assert_eq!(section.shear_property(), None);
    assert_eq!(section.elastic().isotropic_parameters(), Some(bindings[0].as_slice()));
    let current = &section.elastic().layers()[0].material;
    for absent in ["density", "young_modulus", "shear_modulus"] {
        assert!(current.property(absent).is_none());
    }
    near(current.property("specimen_density").unwrap().value_si(), 900.0, 1e-14);
    near(hot.reference_mass().total_mass_kg(), 1000.0 * area * 0.2, 1e-14);
    assert_eq!(hot.reference_mass(), source.reference_rod().reference_mass());
    for property in current.properties() {
        let original = source.reference_rod().sections()[0].elastic().layers()[0].material
            .property(property.requirement().name()).unwrap();
        assert_eq!(property.answer().receipt.selected, original.answer().receipt.selected);
    }
    let load = TipLoad { force: [0.1, 0.0, 0.0], moment: [1e-4, 0.0, 0.0] };
    hot.solve_static(&load, 1, 1e-10, &gate).unwrap();
    near(hot.positions_m()[1][0] - 0.2,
        0.2 * (0.0015 + load.force[0] / (young * area)), 2e-5);
    let (_, curvature) = hot.segment_strains(0).unwrap();
    near(curvature[0], load.moment[0] / (300_000.0 * 2.0 * inertia), 2e-5);
    let (force, moment) = hot.segment_resultants(0).unwrap();
    near(force[0], load.force[0], 2e-5);
    near(moment[0], load.moment[0], 2e-5);
}

#[test]
fn mixed_elastic_pairs_remain_per_layer_and_per_segment_during_heating() {
    let bulk = bulk_shear_card();
    let poisson = young_poisson_card();
    let direct = card(1.0, false);
    let mut first = segment(&bulk, 0.005);
    first.layers.push(RodLayerSource {
        outer_radius_m: 0.01, card: &poisson, reference_point: point(),
        mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let bindings = vec![
        vec![bulk_shear_binding(), young_poisson_binding()],
        vec![LayerIsotropicParameters {
            density: positive_requirement("density", Density::DIMS),
            parameters: IsotropicElasticParameters::YoungShear {
                young_modulus: positive_requirement("young_modulus", Pressure::DIMS),
                shear_modulus: positive_requirement("shear_modulus", Pressure::DIMS),
            },
        }],
    ];
    let gate = CancelGate::new();
    let source = RodMaterialSource::try_new_with_isotropic_parameters(
        0.4, vec![first, segment(&direct, 0.01)], &bindings, 0.03, &gate).unwrap();
    let hot = source.resolve_at(&[400.0, 350.0, 450.0], &gate).unwrap();
    let core_area = std::f64::consts::PI * 0.005_f64.powi(2);
    let shell_area = std::f64::consts::PI * (0.01_f64.powi(2) - 0.005_f64.powi(2));
    near(hot.sections()[0].elastic().axial_stiffness_n(),
        787_500.0 * core_area + 1_050_000.0 * shell_area, 1e-14);
    near(hot.sections()[0].shear_stiffness_n(),
        0.8 * (300_000.0 * core_area + 420_000.0 * shell_area), 1e-14);
    near(hot.sections()[1].elastic().axial_stiffness_n(),
        625_000.0 * (core_area + shell_area), 1e-14);
    near(hot.sections()[1].shear_stiffness_n(),
        0.8 * 250_000.0 * (core_area + shell_area), 1e-14);
    let expansion = &hot.thermal_sections().unwrap()[0];
    near(expansion.zero_resultant_axial_strain(),
        (787_500.0 * core_area * 0.0015 + 1_050_000.0 * shell_area * 0.001)
            / hot.sections()[0].elastic().axial_stiffness_n(), 1e-13);
    for (section, binding) in hot.sections().iter().zip(&bindings) {
        assert_eq!(section.elastic().isotropic_parameters(), Some(binding.as_slice()));
    }
    let shell = &hot.sections()[0].elastic().layers()[1].material;
    assert_eq!(shell.card_identity(), poisson.content_hash());
    assert!(shell.property("shear_modulus").is_none());
    near(shell.property("tensile_response").unwrap().value_si(), 1_050_000.0, 1e-14);
    assert_eq!(hot.reference_mass(), source.reference_rod().reference_mass());
}

#[test]
fn source_thermoelastic_law_uses_derived_shear_in_the_same_energy_and_entropy() {
    let card = young_poisson_card();
    let gate = CancelGate::new();
    let source = RodMaterialSource::try_new_with_isotropic_parameters(
        0.2, vec![segment(&card, 0.01)], &[vec![young_poisson_binding()]], 0.03, &gate).unwrap();
    let laws = source.thermoelastic_laws(&gate).unwrap();
    let temperature = 400.0;
    let gamma = [0.005, 0.002, -0.001];
    let curvature = [0.1, 0.15, -0.2];
    let state = laws[0].at_strain(&[temperature], gamma, curvature, 0.03).unwrap();
    let area = std::f64::consts::PI * 0.01_f64.powi(2);
    let inertia = std::f64::consts::PI * 0.01_f64.powi(4) / 4.0;
    let young = 900_000.0;
    let shear = young / 2.5;
    let elastic_strain = gamma[0] - 2e-5 * (temperature - 300.0);
    let bending_geometry = inertia * (curvature[1].powi(2) + curvature[2].powi(2));
    let shear_geometry = 0.8 * area * (gamma[1].powi(2) + gamma[2].powi(2))
        + 2.0 * inertia * curvature[0].powi(2);
    let energy = 0.5 * (young * (area * elastic_strain.powi(2) + bending_geometry)
        + shear * shear_geometry);
    let entropy = young * area * 2e-5 * elastic_strain
        - 0.5 * (-3000.0 * (area * elastic_strain.powi(2) + bending_geometry)
            - 1200.0 * shear_geometry);
    near(state.layers[0].free_energy_j_m, energy, 1e-13);
    near(state.layers[0].entropy_j_m_k, entropy, 1e-13);
    near(state.layers[0].internal_energy_j_m, energy + temperature * entropy, 1e-13);
    near(state.layers[0].shear_modulus_pa, shear, 1e-14);
    near(state.layers[0].shear_modulus_derivative_pa_k, -1200.0, 1e-14);
    let h = 0.01;
    let plus = laws[0].at_strain(&[temperature + h], gamma, curvature, 0.03).unwrap();
    let minus = laws[0].at_strain(&[temperature - h], gamma, curvature, 0.03).unwrap();
    near(-(plus.layers[0].free_energy_j_m - minus.layers[0].free_energy_j_m) / (2.0 * h),
        entropy, 1e-7);
    let hot = source.resolve_at(&[temperature], &gate).unwrap();
    near(state.section.elastic().axial_stiffness_n(),
        hot.sections()[0].elastic().axial_stiffness_n(), 1e-14);
    near(state.section.torsional_stiffness_n_m2(),
        hot.sections()[0].torsional_stiffness_n_m2(), 1e-14);
    assert_eq!(state.section.shear_property(), None);
    assert!(state.section.elastic().layers()[0].material.property("shear_modulus").is_none());
}

#[test]
fn thermoelastic_source_admits_scalar_bulk_shear_but_refuses_unsupported_thermal_curvature() {
    let gate = CancelGate::new();
    for varying_moduli in [false, true] {
        let moduli = if varying_moduli {
            [claim("bulk_response", Pressure::DIMS, vec![(300.0, 8e5), (500.0, 6e5)]),
             claim("shear_response", Pressure::DIMS, vec![(300.0, 4e5), (500.0, 2e5)])]
        } else {
            [scalar_claim("bulk_response", Pressure::DIMS, 8e5),
             scalar_claim("shear_response", Pressure::DIMS, 4e5)]
        };
        let card = material_from_claims(vec![
            scalar_claim("specimen_density", Density::DIMS, 1000.0),
            scalar_claim(ALPHA, INVERSE_TEMPERATURE_DIMS, 1e-5),
        ].into_iter().chain(moduli).collect());
        let source = RodMaterialSource::try_new_with_isotropic_parameters(
            0.2, vec![segment(&card, 0.01)], &[vec![bulk_shear_binding()]], 0.03, &gate).unwrap();
        let hot = source.resolve_at(&[400.0], &gate).unwrap();
        let result = source.thermoelastic_laws(&gate);
        if varying_moduli {
            assert!(format!("{}", result.unwrap_err()).contains("derived nonlinear E"));
        } else {
            let laws = result.unwrap();
            let state = laws[0].at_strain(&[400.0], [0.002, 0.001, 0.0], [0.1, 0.2, 0.0], 0.03).unwrap();
            assert_eq!(state.layers[0].shear_modulus_derivative_pa_k, 0.0);
            assert_eq!(state.axial.layers[0].young_modulus_derivative_pa_k, 0.0);
            near(state.section.elastic().axial_stiffness_n(),
                hot.sections()[0].elastic().axial_stiffness_n(), 1e-14);
            assert!(state.section.elastic().layers()[0].material.property("young_modulus").is_none());
        }
        assert_eq!(hot.reference_mass(), source.reference_rod().reference_mass());
    }
}

#[test]
fn parameterized_sources_refuse_incomplete_roles_and_cannot_switch_bulk_claims() {
    let card = material_from_claims(vec![
        scalar_claim("specimen_density", Density::DIMS, 1000.0),
        claim("bulk_response", Pressure::DIMS, vec![(300.0, 8e5), (350.0, 7e5)]),
        claim("bulk_response", Pressure::DIMS, vec![(400.0, 6e5), (500.0, 5e5)]),
        scalar_claim("shear_response", Pressure::DIMS, 4e5),
        scalar_claim(ALPHA, INVERSE_TEMPERATURE_DIMS, 1e-5),
    ]);
    let gate = CancelGate::new();
    for bindings in [vec![], vec![vec![]], vec![vec![bulk_shear_binding(), bulk_shear_binding()]]] {
        assert!(RodMaterialSource::try_new_with_isotropic_parameters(
            0.2, vec![segment(&card, 0.01)], &bindings, 0.03, &gate).is_err());
    }
    let bindings = [vec![bulk_shear_binding()]];
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(RodMaterialSource::try_new_with_isotropic_parameters(
        0.2, vec![segment(&card, 0.01)], &bindings, 0.03, &cancelled).is_err());
    let source = RodMaterialSource::try_new_with_isotropic_parameters(
        0.2, vec![segment(&card, 0.01)], &bindings, 0.03, &gate).unwrap();
    assert!(source.resolve_at(&[325.0], &gate).is_ok());
    // The high-temperature claim exists, but was not selected at the reference.
    assert!(source.resolve_at(&[425.0], &gate).is_err());
    assert_eq!(source.reference_temperatures_k(), &[300.0]);
    assert!(source.reference_rod().sections()[0].elastic().layers()[0].material
        .property("young_modulus").is_none());
}
