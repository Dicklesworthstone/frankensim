//! Tensile capacity and helical-wrap inertia of string sections. Synthetic
//! coefficients test the section physics; no measured alloy is implied.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, EvaluationDecision, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::wound::{HelicalWrap, HelicalWrapping, MAX_HELICAL_WRAPS, close_wound_pitch_m};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ResolvedScalarProperty, ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, Length, Pressure, QuantitySpec};
use std::f64::consts::PI;

const LIMIT: &str = "tensile_limit";

fn material(values: &[(&str, Dims, f64)]) -> ResolvedMaterialStatePoint { material_at(values, 293.15) }
fn material_at(values: &[(&str, Dims, f64)], temperature: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in values {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "synthetic string section".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "test".into(), phase: "solid".into(),
        process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", temperature).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn solid(rho: f64, young: f64, limit: f64) -> ResolvedMaterialStatePoint {
    material(&[("density", Density::DIMS, rho), ("young_modulus", Pressure::DIMS, young),
        (LIMIT, Pressure::DIMS, limit)])
}
fn rod(layers: &[(f64, ResolvedMaterialStatePoint)]) -> ConcentricElasticSection {
    ConcentricElasticSection::try_new(0.0, layers.iter().map(|(radius, material)|
        ConcentricMaterialLayer { outer_radius_m: *radius, material: material.clone() }).collect()).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs(), "{a:e} != {b:e} (relative {tolerance:e})");
}
/// A synthetic limit that falls linearly from 400 MPa at 0.8 mm to 300 MPa at
/// 1.2 mm wire diameter, resolved at `diameter_m`: an answer no layer's fixed
/// state bundle can hold.
fn limit_by_diameter(diameter_m: f64) -> ResolvedScalarProperty {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(LIMIT, Pressure::DIMS),
        value: PropertyValue::Curve { abscissa: "wire_diameter".into(), abscissa_dims: Length::DIMS,
            knots: vec![(0.8e-3, 4.0e8), (1.2e-3, 3.0e8)], dims: Pressure::DIMS },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0)
            .with_quantity("wire_diameter", QuantitySpec::dimensional(Length::DIMS), 0.8e-3, 1.2e-3),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic strength by diameter".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "test".into(), phase: "solid".into(),
        process: "drawn".into(), revision: 0 }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 293.15).unwrap()
        .with_quantity("wire_diameter", QuantitySpec::dimensional(Length::DIMS), diameter_m).unwrap();
    let requirement = ScalarPropertyRequirement::try_new(LIMIT, Pressure::DIMS,
        ScalarAdmissibility::StrictlyPositive).unwrap();
    resolve_material_state_point(&card, &point, &[requirement], MaterialPropertySelection::SingleClaimOnly)
        .unwrap().property(LIMIT).unwrap().clone()
}
fn wrap(d: f64, pitch: f64, rho: f64) -> HelicalWrap {
    HelicalWrap { wire_diameter_m: d, pitch_m: pitch, material: solid(rho, 1.2e11, 2e8) }
}

// G0: a limit resolved at another state point (a strength claim over wire
// diameter, evaluated at 1.0 mm) enters the capacity only as its complete
// fs-matdb answer, and the capacity keeps that answer with its claim and its
// interpolation bracket. The claim path keeps each layer's own answer the same
// way. A count other than the layer count, an answer for another property and
// a nonpositive limit refuse.
#[test]
fn g0_resolved_limits_keep_their_answers_and_match_the_claim_path() {
    let section = rod(&[(0.0004, solid(7850.0, 2.0e11, 2.0e9)), (0.0005, solid(19300.0, 8.0e10, 2.0e8))]);
    let claimed = section.axial_tensile_capacity(LIMIT).unwrap();
    let own: Vec<ResolvedScalarProperty> = section.layers().iter()
        .map(|layer| layer.material.property(LIMIT).unwrap().clone()).collect();
    assert_eq!(claimed.layer_limit_sources(), own.as_slice());
    assert_eq!(section.axial_tensile_capacity_from_resolved(LIMIT, &own).unwrap(), claimed);

    let by_diameter = limit_by_diameter(1.0e-3);
    let decision = &by_diameter.answer().receipt.decision;
    assert!(matches!(decision, EvaluationDecision::LinearInside { x_lo, x_hi }
        if *x_lo == 0.8e-3 && *x_hi == 1.2e-3), "{decision:?}");
    let mixed = section.axial_tensile_capacity_from_resolved(LIMIT, &[own[0].clone(), by_diameter.clone()])
        .unwrap();
    near(mixed.layer_limits_pa()[1], 3.5e8, 1e-12);
    assert_eq!(mixed.layer_limit_sources()[1], by_diameter);
    // limit/E: 2e9/2e11 = 1e-2 for the core, 3.5e8/8e10 = 4.4e-3 for the coat, which governs.
    assert_eq!(mixed.governing_layer(), 1);
    near(mixed.capacity_n(), mixed.layer_limits_pa()[1] * section.axial_stiffness_n() / 8.0e10, 2e-15);
    // A layer whose own bundle has no limit claim takes a resolved limit.
    let unrated = material(&[("density", Density::DIMS, 7850.0), ("young_modulus", Pressure::DIMS, 2e11)]);
    let bare = rod(&[(0.0005, unrated)]);
    let capacity = bare.axial_tensile_capacity_from_resolved(LIMIT, std::slice::from_ref(&by_diameter)).unwrap();
    near(capacity.capacity_n(), capacity.layer_limits_pa()[0] * PI * 0.0005 * 0.0005, 2e-15);

    assert!(section.axial_tensile_capacity_from_resolved(LIMIT, &own[..1]).is_err());
    assert!(section.axial_tensile_capacity_from_resolved(LIMIT, &[own[0].clone(), own[1].clone(), own[1].clone()])
        .is_err());
    let modulus = section.layers()[0].material.property("young_modulus").unwrap().clone();
    assert!(section.axial_tensile_capacity_from_resolved(LIMIT, &[own[0].clone(), modulus]).is_err());
    assert!(section.axial_tensile_capacity_from_resolved("another_limit", &own).is_err());
    let zero = material(&[(LIMIT, Pressure::DIMS, 0.0)]).property(LIMIT).unwrap().clone();
    assert!(section.axial_tensile_capacity_from_resolved(LIMIT, &[own[0].clone(), zero]).is_err());
}

// G0: one material: capacity = limit * area, and at a fixed pitch the stress
// 4 L^2 f^2 rho does not depend on the wire diameter.
#[test]
fn g0_single_layer_capacity_is_limit_times_area_and_pitch_fixes_stress() {
    let (rho, young, limit) = (7850.0, 2.0e11, 2.0e9);
    let section = rod(&[(0.0005, solid(rho, young, limit))]);
    let capacity = section.axial_tensile_capacity(LIMIT).unwrap();
    assert_eq!(capacity.limit_property(), LIMIT);
    assert_eq!(capacity.governing_layer(), 0);
    assert_eq!(capacity.layer_limits_pa(), &[limit]);
    near(capacity.capacity_n(), limit * PI * 0.0005 * 0.0005, 2e-15);
    near(capacity.layer_stress_pa(0, capacity.capacity_n()).unwrap(), limit, 2e-15);
    near(capacity.utilization(capacity.capacity_n()).unwrap(), 1.0, 2e-15);
    let (length, frequency): (f64, f64) = (0.404, 440.0);
    for radius in [0.0003, 0.0005, 0.0008] {
        let section = rod(&[(radius, solid(rho, young, limit))]);
        // Flexible-string tension that puts the fundamental at `frequency`.
        let tension = section.line_density_kg_m() * (2.0 * length * frequency).powi(2);
        let stress = section.axial_tensile_capacity(LIMIT).unwrap().layer_stress_pa(0, tension).unwrap();
        near(stress, 4.0 * length * length * frequency * frequency * rho, 4e-15);
    }
}

// G0: bonded layers share one strain; the layer with the smallest limit/E
// governs, and its stress equals its limit exactly at the capacity.
#[test]
fn g0_bonded_layers_share_one_strain_and_the_smallest_limit_ratio_governs() {
    let core = solid(7850.0, 2.0e11, 2.0e9);
    let coat = solid(19300.0, 8.0e10, 2.0e8);
    let section = rod(&[(0.0004, core.clone()), (0.0005, coat.clone())]);
    let capacity = section.axial_tensile_capacity(LIMIT).unwrap();
    let force = 100.0;
    let stresses = [0, 1].map(|i| capacity.layer_stress_pa(i, force).unwrap());
    let areas = section.layer_integrals().iter().map(|g| g.area_m2()).collect::<Vec<_>>();
    near(stresses[0] * areas[0] + stresses[1] * areas[1], force, 4e-15);
    near(stresses[0] / 2.0e11, stresses[1] / 8.0e10, 2e-15);
    assert_eq!(capacity.governing_layer(), 1);
    near(capacity.capacity_n(), 2.0e8 * section.axial_stiffness_n() / 8.0e10, 2e-15);
    near(capacity.layer_stress_pa(1, capacity.capacity_n()).unwrap(), 2.0e8, 4e-15);
    assert!(capacity.layer_stress_pa(0, capacity.capacity_n()).unwrap() < 2.0e9);
    // A stronger coat moves the governing layer to the core.
    let strong = rod(&[(0.0004, core), (0.0005, solid(19300.0, 8.0e10, 1.0e9))]);
    let capacity = strong.axial_tensile_capacity(LIMIT).unwrap();
    assert_eq!(capacity.governing_layer(), 0);
    near(capacity.capacity_n(), 2.0e9 * strong.axial_stiffness_n() / 2.0e11, 2e-15);
}

// G0: splitting one material into annuli is not physics; capacity is unchanged.
#[test]
fn g0_layer_subdivision_does_not_change_capacity() {
    let m = solid(7850.0, 2.0e11, 2.0e9);
    let whole = rod(&[(0.0005, m.clone())]).axial_tensile_capacity(LIMIT).unwrap();
    let split = rod(&[(0.0002, m.clone()), (0.00035, m.clone()), (0.0005, m)])
        .axial_tensile_capacity(LIMIT).unwrap();
    near(split.capacity_n(), whole.capacity_n(), 4e-15);
    assert_eq!(split.layer_limits_pa(), &[2.0e9; 3]);
}

// Planted negatives: absent strength data is missing data, never "unbreakable".
#[test]
fn g0_missing_or_invalid_limit_data_refuses() {
    let unrated = material(&[("density", Density::DIMS, 7850.0), ("young_modulus", Pressure::DIMS, 2e11)]);
    assert!(rod(&[(0.0005, unrated.clone())]).axial_tensile_capacity(LIMIT).is_err());
    // A rated core does not excuse an unrated coat.
    let mixed = rod(&[(0.0004, solid(7850.0, 2e11, 2e9)), (0.0005, unrated)]);
    assert!(mixed.axial_tensile_capacity(LIMIT).is_err());
    let wrong_quantity = material(&[("density", Density::DIMS, 7850.0),
        ("young_modulus", Pressure::DIMS, 2e11), (LIMIT, Density::DIMS, 2e9)]);
    assert!(rod(&[(0.0005, wrong_quantity)]).axial_tensile_capacity(LIMIT).is_err());
    for bad in [0.0, -2e9] {
        assert!(rod(&[(0.0005, solid(7850.0, 2e11, bad))]).axial_tensile_capacity(LIMIT).is_err());
    }
    let capacity = rod(&[(0.0005, solid(7850.0, 2e11, 2e9))]).axial_tensile_capacity(LIMIT).unwrap();
    assert!(capacity.layer_stress_pa(1, 10.0).is_err());
    for force in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(capacity.layer_stress_pa(0, force).is_err());
        assert!(capacity.utilization(force).is_err());
    }
}

// Independent closest approach of two different turns of a centre-line helix
// (radius R, pitch p), by dense sampling of the angular separation over
// [pi, 3 pi]; separations below half a turn belong to the same wire.
fn sampled_clearance(radius: f64, pitch: f64) -> f64 {
    let steps = 400_000;
    (0..=steps).map(|i| {
        let turn = PI + 2.0 * PI * i as f64 / steps as f64;
        (2.0 * radius * (0.5 * turn).sin()).hypot(pitch * turn / (2.0 * PI))
    }).fold(f64::INFINITY, f64::min)
}

// G0: the tightest admissible pitch is where different turns just touch; a
// pitch of exactly one wire diameter already overlaps them.
#[test]
fn g0_close_wound_pitch_is_where_neighbouring_turns_touch() {
    for (under, d) in [(0.0012, 0.0012), (0.0032, 0.0019), (0.0002, 0.0001), (0.0005, 0.0020)] {
        let pitch = close_wound_pitch_m(under, d).unwrap();
        assert!(pitch > d && pitch < 2.0 * d);
        let radius = 0.5 * (under + d);
        near(sampled_clearance(radius, pitch), d, 1e-9);
        assert!(sampled_clearance(radius, d) < d);
        let core = rod(&[(0.5 * under, solid(7850.0, 2e11, 2e9))]);
        assert!(HelicalWrapping::try_new(&core, vec![wrap(d, pitch, 8960.0)]).is_ok());
        assert!(HelicalWrapping::try_new(&core, vec![wrap(d, pitch * (1.0 - 1e-9), 8960.0)]).is_err());
        assert!(HelicalWrapping::try_new(&core, vec![wrap(d, d, 8960.0)]).is_err());
    }
    for (under, d) in [(0.0, 0.001), (0.001, 0.0), (f64::NAN, 0.001), (0.001, f64::INFINITY)] {
        assert!(close_wound_pitch_m(under, d).is_err());
    }
}

// G0: a wrap's mass is its wire volume along the centre-line helix.
#[test]
fn g0_wrap_mass_is_the_helix_integral() {
    let (r, d, rho) = (0.0006, 0.0012, 8960.0);
    let pitch = close_wound_pitch_m(2.0 * r, d).unwrap();
    let core = rod(&[(r, solid(7850.0, 2e11, 2e9))]);
    let wrapping = HelicalWrapping::try_new(&core, vec![wrap(d, pitch, rho)]).unwrap();
    let turns = 1.0 / pitch;
    let turn_length = ((PI * (2.0 * r + d)).powi(2) + pitch * pitch).sqrt();
    near(wrapping.line_density_kg_m(), rho * PI * d * d / 4.0 * turns * turn_length, 4e-15);
    assert_eq!(wrapping.wrap_line_density_kg_m(), &[wrapping.line_density_kg_m()]);
    assert_eq!(wrapping.core_outer_radius_m(), r);
    near(wrapping.outer_radius_m(), r + d, 2e-15);
    // Widening the pitch thins the helix toward straight axial wire.
    let open = HelicalWrapping::try_new(&core, vec![wrap(d, 1.0e3, rho)]).unwrap();
    assert!(open.line_density_kg_m() < wrapping.line_density_kg_m());
    near(open.line_density_kg_m(), rho * PI * d * d / 4.0, 1e-10);
    assert!(wrapping.line_density_kg_m() > core.line_density_kg_m());
}

// G0: a second wrap lies on the first, so its helix sits two first-wire
// diameters further out and the outer radius adds both wires.
#[test]
fn g0_second_wrap_lies_on_the_first() {
    let (r, d1, d2) = (0.0008, 0.0010, 0.0016);
    let (p1, p2) = (close_wound_pitch_m(2.0 * r, d1).unwrap(),
        close_wound_pitch_m(2.0 * r + 2.0 * d1, d2).unwrap());
    let core = rod(&[(r, solid(7850.0, 2e11, 2e9))]);
    let both = HelicalWrapping::try_new(&core, vec![wrap(d1, p1, 8960.0), wrap(d2, p2, 8960.0)]).unwrap();
    let helix = |diameter: f64, d: f64, p: f64| 8960.0 * PI * d * d / 4.0
        * (1.0 + (PI * (diameter + d) / p).powi(2)).sqrt();
    let masses = both.wrap_line_density_kg_m();
    near(masses[0], helix(2.0 * r, d1, p1), 4e-15);
    near(masses[1], helix(2.0 * r + 2.0 * d1, d2, p2), 4e-15);
    near(both.line_density_kg_m(), masses[0] + masses[1], 2e-15);
    near(both.outer_radius_m(), r + d1 + d2, 2e-15);
    // The second wire wound directly on the core is a different string: on
    // that smaller helix the pitch solved for the outer position overlaps.
    assert!(HelicalWrapping::try_new(&core, vec![wrap(d2, p2, 8960.0)]).is_err());
    let alone = HelicalWrapping::try_new(&core,
        vec![wrap(d2, close_wound_pitch_m(2.0 * r, d2).unwrap(), 8960.0)]).unwrap();
    assert!(alone.line_density_kg_m() < masses[1]);
}

#[test]
fn g0_empty_oversized_mismatched_or_unsourced_wraps_refuse() {
    let core = rod(&[(0.0006, solid(7850.0, 2e11, 2e9))]);
    assert!(HelicalWrapping::try_new(&core, vec![]).is_err());
    assert!(HelicalWrapping::try_new(&core, vec![wrap(0.001, 0.002, 8960.0); MAX_HELICAL_WRAPS + 1]).is_err());
    assert!(HelicalWrapping::try_new(&core, vec![wrap(0.001, 0.002, 8960.0); MAX_HELICAL_WRAPS]).is_ok());
    for d in [0.0, -0.001, f64::NAN, f64::INFINITY] {
        assert!(HelicalWrapping::try_new(&core, vec![wrap(d, 0.002, 8960.0)]).is_err());
    }
    for pitch in [0.0, -0.002, f64::NAN] {
        assert!(HelicalWrapping::try_new(&core, vec![wrap(0.001, pitch, 8960.0)]).is_err());
    }
    let no_density = HelicalWrap { wire_diameter_m: 0.001, pitch_m: 0.002,
        material: material(&[("young_modulus", Pressure::DIMS, 1.2e11)]) };
    assert!(HelicalWrapping::try_new(&core, vec![no_density]).is_err());
    // A wrap resolved at another temperature is a different physical state.
    let warm = HelicalWrap { wire_diameter_m: 0.001, pitch_m: 0.002,
        material: material_at(&[("density", Density::DIMS, 8960.0)], 300.0) };
    assert!(HelicalWrapping::try_new(&core, vec![warm]).is_err());
}

// ---- Loss-angle (hysteretic) bending loss ----

const TAN: &str = "bending_loss_tangent";
fn lossy(rho: f64, young: f64, tangent: f64, dims: Dims, band: Option<(f64, f64)>) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [("density", Density::DIMS, rho), ("young_modulus", Pressure::DIMS, young), (TAN, dims, tangent)] {
        let validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        let validity = match band { Some((lo, hi)) if name == TAN => validity.with("omega", lo, hi), _ => validity };
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims }, validity,
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "synthetic loss tangent".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "test".into(), phase: "solid".into(),
        process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 5000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}

// G0: H = sum(E_i I_i tan(delta_i)) from independent annular integrals; the
// bands intersect; subdivision does not change H; explicit zero is lossless.
#[test]
fn g0_loss_angle_bending_integrates_each_layer_without_averaging() {
    let (r1, r2) = (0.0004, 0.0006);
    let a = lossy(7850.0, 2.1e11, 1.1e-4, Dims::NONE, Some((50.0, 2.0e4)));
    let b = lossy(19300.0, 7.9e10, 5.0e-3, Dims::NONE, Some((10.0, 1.2e4)));
    let section = rod(&[(r1, a.clone()), (r2, b)]).with_loss_angle_bending(TAN).unwrap();
    let moments = [PI * r1.powi(4) / 4.0, PI * (r2.powi(4) - r1.powi(4)) / 4.0];
    let layers = [2.1e11 * moments[0] * 1.1e-4, 7.9e10 * moments[1] * 5.0e-3];
    near(section.layer_hysteretic_stiffness_n_m2()[0], layers[0], 1e-14);
    near(section.layer_hysteretic_stiffness_n_m2()[1], layers[1], 1e-14);
    near(section.hysteretic_stiffness_n_m2(), layers[0] + layers[1], 1e-14);
    assert_eq!(section.omega_band_rad_s(), (50.0, 1.2e4));
    assert_eq!(section.loss_property(), TAN);
    let whole = rod(&[(r2, a.clone())]).with_loss_angle_bending(TAN).unwrap();
    let split = rod(&[(r1, a.clone()), (r2, a)]).with_loss_angle_bending(TAN).unwrap();
    near(split.hysteretic_stiffness_n_m2(), whole.hysteretic_stiffness_n_m2(), 1e-14);
    let zero = rod(&[(r2, lossy(7850.0, 2.1e11, 0.0, Dims::NONE, Some((1.0, 1e5))))]);
    assert_eq!(zero.with_loss_angle_bending(TAN).unwrap().hysteretic_stiffness_n_m2(), 0.0);
}

// The modal reduction dissipates only the bending part of a tensioned mode.
// Added fluid mass changes its frequency/rate, not the source H or dry section.
#[test]
fn g1_modal_loss_angle_uses_actual_frequency_mass_and_tension_dilution() {
    use fs_material::section::modal_loss::ConcentricModalLossSection;
    let original = rod(&[(0.0004, lossy(7850.0, 2.1e11, 1.1e-4, Dims::NONE,
        Some((50.0, 2.0e4)))),
        (0.0006, lossy(19300.0, 7.9e10, 5e-3, Dims::NONE, Some((10.0, 1.2e4))))]);
    let source = original.clone().with_loss_angle_bending(TAN).unwrap();
    let h = source.hysteretic_stiffness_n_m2();
    let loss = ConcentricModalLossSection::from_loss_angle(source.clone());
    assert_eq!(loss.elastic(), &original);
    assert_eq!(loss.loss_angle(), Some(&source));
    assert!(loss.thermoelastic().is_none());
    assert_eq!(loss.omega_band_rad_s(), (50.0, 1.2e4));
    let wave: f64 = PI / 0.7;
    let stiffness = 100.0 * wave.powi(2) + original.bending_stiffness_n_m2() * wave.powi(4);
    let dry_mass = original.line_density_kg_m();
    let mut dry_rate = 0.0;
    for added_mass in [0.0, dry_mass] {
        let mass = dry_mass + added_mass;
        let omega = (stiffness / mass).sqrt();
        let rate = loss.damping_rate_per_s(omega, wave, mass).unwrap();
        near(rate / (2.0 * omega), h * wave.powi(4) / (2.0 * stiffness), 1e-14);
        if added_mass == 0.0 { dry_rate = rate; }
        else { near(rate / dry_rate, 1.0 / 2.0_f64.sqrt(), 1e-14); }
        assert!(rate / (2.0 * omega) < 0.5 * h / original.bending_stiffness_n_m2());
    }
    for omega in [0.0, 49.0, 12001.0, f64::NAN] {
        assert!(loss.damping_rate_per_s(omega, wave, dry_mass).is_err());
    }
    assert!(loss.damping_rate_per_s(1000.0, 0.0, dry_mass).is_err());
    assert!(loss.damping_rate_per_s(1000.0, wave, -dry_mass).is_err());
}

// Planted negatives: a viscosity in Pa s, a negative tangent, a band-less
// claim and an absent claim never become a loss tangent.
#[test]
fn g0_loss_angle_refuses_wrong_quantity_negative_bandless_or_missing_data() {
    let r = 0.0005;
    for bad in [lossy(7850.0, 2e11, 1e-4, Pressure::DIMS, Some((1.0, 1e5))),
        lossy(7850.0, 2e11, -1e-4, Dims::NONE, Some((1.0, 1e5))),
        lossy(7850.0, 2e11, 1e-4, Dims::NONE, None)] {
        assert!(rod(&[(r, bad)]).with_loss_angle_bending(TAN).is_err());
    }
    let rated = lossy(7850.0, 2e11, 1e-4, Dims::NONE, Some((1.0, 1e5)));
    assert!(rod(&[(r, rated)]).with_loss_angle_bending("absent_loss_tangent").is_err());
}
