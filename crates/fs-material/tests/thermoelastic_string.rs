//! Zener thermoelastic bending loss of solid circular members. Synthetic
//! coefficients test the law; no measured alloy is implied.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, ResolvedMaterialStatePoint,
    SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    resolve_material_state_point};
use fs_material::visco::circular_ted::{CIRCULAR_TED_MODES, circular_relaxation_shape, circular_thermal_modes};
use fs_math::special::bessel_j1;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use std::f64::consts::PI;

// A gold-like synthetic metal: E, alpha, rho, c_p, k.
const METAL: [f64; 5] = [71.0e9, 14.2e-6, 19_300.0, 129.7, 296.0];
const T0: f64 = 293.15;

fn kelvin() -> QuantitySpec { QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static)) }
fn state(values: &[(&str, Dims, f64)], typed_t: bool) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in values {
        let validity = if typed_t { ValidityDomain::unconstrained().with_quantity("T", kelvin(), 280.0, 320.0) }
            else { ValidityDomain::unconstrained().with("T", 280.0, 320.0) };
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims }, validity,
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "synthetic thermoelastic metal".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "test".into(), phase: "solid".into(),
        process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap();
    let point = if typed_t { QueryPoint::new().with_quantity("T", kelvin(), T0).unwrap() }
        else { QueryPoint::new().with("T", T0).unwrap() };
    resolve_material_state_point(&card, &point, &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn metal(values: [f64; 5]) -> ResolvedMaterialStatePoint {
    state(&[("young_modulus", Pressure::DIMS, values[0]),
        ("linear_thermal_expansion_coefficient", INVERSE_TEMPERATURE_DIMS, values[1]),
        ("density", Density::DIMS, values[2]),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, values[3]),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, values[4])], true)
}
fn rod(radius: f64, material: ResolvedMaterialStatePoint) -> ConcentricElasticSection {
    ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer { outer_radius_m: radius, material }]).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs(), "{a:e} != {b:e} (relative {tolerance:e})");
}

// G0: xi_n are the zeros of J1' (checked against fs-math's certified J1, not
// the quadrature that found them), and w_n equals the projection of the
// bending strain r sin(theta) onto J1(xi_n r) sin(theta), integrated
// independently; the weights sum to one.
#[test]
fn g0_thermal_modes_are_zeros_of_j1_prime_with_projection_weights() {
    let modes = circular_thermal_modes();
    assert_eq!(modes.len(), CIRCULAR_TED_MODES);
    // Abramowitz and Stegun Table 9.5: j'_{1,1} = 1.8411837813, j'_{1,2} = 5.3314427735.
    near(modes[0].0, 1.841_183_781_340_659_5, 1e-12);
    near(modes[1].0, 5.331_442_773_525_032_6, 1e-12);
    // Heptonstall et al. eq. 6: f_c = 2.16 k/(rho C d^2), i.e. 4 xi_1^2 / (2 pi).
    assert!((4.0 * modes[0].0.powi(2) / (2.0 * PI) - 2.16).abs() < 0.005);
    let mut previous = 0.0;
    for (index, &(xi, weight)) in modes.iter().enumerate() {
        let n = (index + 1) as f64;
        assert!(xi > previous && xi > (n - 0.5) * PI && xi < (n + 0.5) * PI);
        previous = xi;
        let h = 1e-5;
        let slope = (bessel_j1(xi + h) - bessel_j1(xi - h)) / (2.0 * h);
        assert!(slope.abs() < 1e-8, "J1'({xi}) = {slope:e}");
        if index < 4 {
            // <r, J1(xi r)> and <J1, J1> on the unit disk (weight r dr) by composite Simpson.
            let steps = 20_000;
            let simpson = |f: &dyn Fn(f64) -> f64| {
                let dx = 1.0 / steps as f64;
                (0..=steps).map(|i| {
                    let c = if i == 0 || i == steps { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
                    c * f(i as f64 * dx)
                }).sum::<f64>() * dx / 3.0
            };
            let overlap = simpson(&|r| r * bessel_j1(xi * r) * r);
            let norm = simpson(&|r| bessel_j1(xi * r).powi(2) * r);
            near(overlap * overlap / (0.25 * norm), weight, 1e-9);
        }
    }
    let total: f64 = modes.iter().map(|m| m.1).sum();
    assert!(total < 1.0 && total > 1.0 - 2e-6, "{total}");
    near(modes[0].1, 0.987_4, 2e-4);
}

// G1: Delta and the base time by hand, the low- and high-frequency limits of
// the series, and the single-mode peak bound.
#[test]
fn g1_relaxation_strength_time_and_limits_by_hand() {
    let radius = 0.0005;
    let section = rod(radius, metal(METAL)).with_thermoelastic_bending().unwrap();
    let [e, alpha, rho, cp, k] = METAL;
    let delta = e * alpha * alpha * T0 / (rho * cp);
    near(section.law().relaxation_strength(), delta, 1e-14);
    near(section.base_time_s(), rho * cp * radius * radius / k, 1e-14);
    near(section.hysteretic_scale_n_m2(), e * PI * radius.powi(4) / 4.0 * delta, 1e-13);
    assert_eq!(section.law().t0, T0);
    let modes = circular_thermal_modes();
    let base = section.base_time_s();
    // omega tau << 1: Delta omega sum(w tau); omega tau >> 1: Delta sum(w / tau) / omega.
    let low = 1e-6 / base;
    near(section.loss_tangent(low), delta * low * modes.iter().map(|m| m.1 * base / m.0.powi(2)).sum::<f64>(), 1e-9);
    let high = 1e7 / base;
    near(section.loss_tangent(high), delta / high * modes.iter().map(|m| m.1 * m.0.powi(2) / base).sum::<f64>(), 1e-6);
    // At the first mode's corner the first term is w_1/2; the rest is bounded by their weight.
    let corner = modes[0].0.powi(2) / base;
    let rest: f64 = modes[1..].iter().map(|m| m.1).sum();
    assert!((section.loss_tangent(corner) - delta * modes[0].1 / 2.0).abs() <= delta * rest);
    // Gold-like 1 mm wire: a peak of about 8e-4 near 0.25 kHz, as estimated in the U2 ruling.
    let peak = (0..2000).map(|i| 10f64.powf(1.0 + 3.0 * i as f64 / 2000.0))
        .map(|f| (f, section.loss_tangent(2.0 * PI * f))).fold((0.0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
    assert!(peak.0 > 200.0 && peak.0 < 320.0 && peak.1 > 7e-4 && peak.1 < 9e-4, "{peak:?}");
}

// G1: the law depends on omega only through omega R^2: doubling the radius
// moves the whole curve down by four in frequency, so the peak moves as 1/d^2.
#[test]
fn g1_doubling_the_radius_moves_the_curve_by_four() {
    let small = rod(0.0004, metal(METAL)).with_thermoelastic_bending().unwrap();
    let large = rod(0.0008, metal(METAL)).with_thermoelastic_bending().unwrap();
    near(large.base_time_s(), 4.0 * small.base_time_s(), 1e-14);
    for f in [10.0, 63.0, 250.0, 1000.0, 8000.0] {
        let omega = 2.0 * PI * f;
        near(large.loss_tangent(omega / 4.0), small.loss_tangent(omega), 1e-12);
    }
    assert_eq!(circular_relaxation_shape(0.0, 1.0), 0.0);
    assert_eq!(circular_relaxation_shape(1.0, f64::NAN), 0.0);
}

// Planted negatives: a missing conductivity or heat capacity, a wrong
// quantity, a layered or hollow section, and an untyped temperature refuse.
#[test]
fn g0_missing_wrong_layered_hollow_or_untyped_inputs_refuse() {
    let full = [("young_modulus", Pressure::DIMS, METAL[0]),
        ("linear_thermal_expansion_coefficient", INVERSE_TEMPERATURE_DIMS, METAL[1]),
        ("density", Density::DIMS, METAL[2]),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, METAL[3]),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, METAL[4])];
    for skip in [3, 4] {
        let partial: Vec<_> = full.iter().enumerate().filter(|(i, _)| *i != skip).map(|(_, v)| *v).collect();
        assert!(rod(0.0005, state(&partial, true)).with_thermoelastic_bending().is_err());
    }
    let mut wrong = full;
    wrong[3].1 = Pressure::DIMS;
    assert!(rod(0.0005, state(&wrong, true)).with_thermoelastic_bending().is_err());
    let layered = ConcentricElasticSection::try_new(0.0, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0003, material: metal(METAL) },
        ConcentricMaterialLayer { outer_radius_m: 0.0005, material: metal([7.0e10, 14.2e-6, 19_300.0, 129.7, 296.0]) },
    ]).unwrap();
    assert!(layered.with_thermoelastic_bending().is_err());
    // Splitting one material into annuli is not layering.
    let split = ConcentricElasticSection::try_new(0.0, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0003, material: metal(METAL) },
        ConcentricMaterialLayer { outer_radius_m: 0.0005, material: metal(METAL) },
    ]).unwrap().with_thermoelastic_bending().unwrap();
    near(split.base_time_s(), rod(0.0005, metal(METAL)).with_thermoelastic_bending().unwrap().base_time_s(), 1e-14);
    let hollow = ConcentricElasticSection::try_new(0.0001, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0005, material: metal(METAL) }]).unwrap();
    assert!(hollow.with_thermoelastic_bending().is_err());
    assert!(rod(0.0005, state(&full, false)).with_thermoelastic_bending().is_err());
}


// Typed held-moisture curves exercise the same environmental policy as every
// other section law. These manufactured values are not measured material data.
fn modal_environment(moisture_value: f64, typed: bool, extra_strain: bool)
    -> ResolvedMaterialStatePoint
{
    use fs_qty::semantic::MoistureBasis;
    let moisture = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static));
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value, varying, band) in [
        ("young_modulus", Pressure::DIMS, METAL[0], false, (50.0, 12000.0)),
        ("density", Density::DIMS, METAL[2], false, (50.0, 12000.0)),
        ("linear_thermal_expansion_coefficient", INVERSE_TEMPERATURE_DIMS, METAL[1], true, (50.0, 12000.0)),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, METAL[3], false, (100.0, 9000.0)),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, METAL[4], false, (50.0, 12000.0)),
        ("excess_loss_tangent", Dims::NONE, 1e-4, true, (200.0, 8000.0)),
    ] {
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity("T", kelvin(), 280.0, 320.0).with("omega", band.0, band.1);
        validity = if typed { validity.with_quantity("water", moisture, 0.0, 0.2) }
            else { validity.with("water", 0.0, 0.2) };
        if extra_strain { validity = validity.with("strain", 0.0, 0.1); }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if varying { PropertyValue::Curve {
                abscissa: "water".into(), abscissa_dims: Dims::NONE,
                knots: vec![(0.0, value), (0.2, 2.0 * value)], dims,
            }} else { PropertyValue::Scalar { value, dims } },
            validity, uncertainty: UncertaintyModel::Unstated,
            interpolation: if varying { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance {
                source: "manufactured modal-loss held environment".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims,
            ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "test".into(),
        phase: "solid".into(), process: "held environment".into(), revision: 0 },
        claims, vec![]).unwrap();
    let point = QueryPoint::new().with_quantity("T", kelvin(), T0).unwrap()
        .with("omega", 1000.0).unwrap();
    let mut point = if typed { point.with_quantity("water", moisture, moisture_value).unwrap() }
        else { point.with("water", moisture_value).unwrap() };
    if extra_strain { point = point.with("strain", 0.01).unwrap(); }
    resolve_material_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}

#[test]
fn g1_modal_thermoelastic_loss_retains_excess_support_and_loaded_inertia() {
    use fs_material::section::modal_loss::ConcentricModalLossSection;
    let material = modal_environment(0.0, true, false);
    let elastic = rod(0.0005, material.clone());
    let thermal = elastic.clone().with_thermoelastic_bending().unwrap();
    let excess = elastic.clone().with_loss_angle_bending("excess_loss_tangent").unwrap();
    assert_eq!(thermal.omega_band_rad_s(), (100.0, 9000.0));
    let only = ConcentricModalLossSection::from_thermoelastic(thermal.clone(), None).unwrap();
    let loss = ConcentricModalLossSection::from_thermoelastic(
        thermal.clone(), Some(excess.clone())).unwrap();
    assert_eq!(loss.omega_band_rad_s(), (200.0, 8000.0));
    assert_eq!(loss.thermoelastic(), Some(&thermal));
    assert_eq!(loss.loss_angle(), Some(&excess));
    assert_eq!(loss.elastic().layers()[0].material, material);
    let wave: f64 = PI / 0.7;
    let stiffness = 100.0 * wave.powi(2) + elastic.bending_stiffness_n_m2() * wave.powi(4);
    for mass in [elastic.line_density_kg_m(), 2.0 * elastic.line_density_kg_m()] {
        let omega = (stiffness / mass).sqrt();
        let rate = loss.damping_rate_per_s(omega, wave, mass).unwrap();
        let h = excess.hysteretic_stiffness_n_m2()
            + elastic.bending_stiffness_n_m2() * thermal.loss_tangent(omega);
        near(rate / (2.0 * omega), h * wave.powi(4) / (2.0 * stiffness), 1e-14);
        let intrinsic = only.damping_rate_per_s(omega, wave, mass).unwrap();
        near(rate - intrinsic, excess.hysteretic_stiffness_n_m2()
            * wave.powi(4) / mass / omega, 1e-12);
    }
    // A selected excess never disappears just because thermal data go further.
    assert!(only.hysteretic_stiffness_n_m2(150.0).is_ok());
    assert!(loss.hysteretic_stiffness_n_m2(150.0).is_err());
    assert!(only.hysteretic_stiffness_n_m2(9500.0).is_err());
    let wrong_geometry = rod(0.0006, material)
        .with_loss_angle_bending("excess_loss_tangent").unwrap();
    assert!(ConcentricModalLossSection::from_thermoelastic(
        thermal, Some(wrong_geometry)).is_err());
}

#[test]
fn g1_modal_loss_owners_share_held_environment_without_freezing_strain() {
    let dry = rod(0.0005, modal_environment(0.0, true, false));
    let wet = rod(0.0005, modal_environment(0.2, true, false));
    let a = dry.clone().with_thermoelastic_bending().unwrap();
    let b = wet.clone().with_thermoelastic_bending().unwrap();
    near(b.law().relaxation_strength(), 4.0 * a.law().relaxation_strength(), 1e-14);
    assert_eq!(a.base_time_s(), b.base_time_s());
    let ta = dry.with_loss_angle_bending("excess_loss_tangent").unwrap();
    let tb = wet.with_loss_angle_bending("excess_loss_tangent").unwrap();
    near(tb.hysteretic_stiffness_n_m2(), 2.0 * ta.hysteretic_stiffness_n_m2(), 1e-14);
    for state in [modal_environment(0.1, false, false), modal_environment(0.1, true, true)] {
        let section = rod(0.0005, state);
        assert!(section.clone().with_thermoelastic_bending().is_err());
        assert!(section.with_loss_angle_bending("excess_loss_tangent").is_err());
    }
}
