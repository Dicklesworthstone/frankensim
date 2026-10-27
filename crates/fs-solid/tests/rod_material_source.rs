//! G0/G1/G3 source data -> temperature refresh -> actual rod equilibrium.
//! Synthetic constitutive curves test wiring; they do not qualify a real alloy.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection};
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::thermal::source::{RodLayerSource, RodMaterialSource, RodSegmentSource};

fn claim(name: &str, dims: Dims, knots: Vec<(f64, f64)>) -> PropertyClaim {
    let lo = knots[0].0;
    let hi = knots.last().unwrap().0;
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS, knots, dims },
        validity: ValidityDomain::unconstrained().with("T", lo, hi).with("p", 100_000.0, 100_000.0),
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
    QueryPoint::new().with("T", 300.0).unwrap().with("p", 100_000.0).unwrap()
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
