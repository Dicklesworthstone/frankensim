use super::*;
use crate::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use crate::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Density, Dims, Pressure, Temperature};

fn scalar(value: f64, dims: Dims) -> PropertyValue {
    PropertyValue::Scalar { value, dims }
}

fn curve(values: [f64; 3], dims: Dims) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
        knots: [280.0, 350.0, 420.0].into_iter().zip(values).collect(), dims,
    }
}

fn card(scale: f64, alpha: f64, defect: Option<&str>)
    -> (MaterialCard, Vec<ScalarPropertyRequirement>)
{
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, mut value, admissibility) in [
        (DENSITY_PROPERTY, Density::DIMS, scalar(4000.0 * scale, Density::DIMS),
            ScalarAdmissibility::StrictlyPositive),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS,
            curve([1.02e9 * scale, 0.95e9 * scale, 0.88e9 * scale], Pressure::DIMS),
            ScalarAdmissibility::StrictlyPositive),
        (ALPHA, INVERSE_TEMPERATURE_DIMS, scalar(alpha, INVERSE_TEMPERATURE_DIMS),
            ScalarAdmissibility::Finite),
        ("auxiliary", Dims::NONE, curve([1.8, 2.5, 3.2], Dims::NONE),
            ScalarAdmissibility::StrictlyPositive),
    ] {
        if defect == Some(name) {
            value = match name {
                YOUNG_MODULUS_PROPERTY => curve([1.02e9, 0.97e9, 0.88e9], dims),
                DENSITY_PROPERTY => curve([4000.0, 3990.0, 3980.0], dims),
                ALPHA => curve([alpha, alpha, alpha], dims),
                _ => unreachable!(),
            };
        }
        let interpolation = match &value {
            PropertyValue::Scalar { .. } => InterpolationPolicy::ConstantWithinValidity,
            _ => InterpolationPolicy::LinearInside,
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value,
            validity: ValidityDomain::unconstrained().with("T", 280.0, 420.0)
                .with("pressure", 100000.0, 100000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation, observations: vec![],
            provenance: Provenance {
                source: "synthetic affine layered thermoelasticity; no measurement claim".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, admissibility).unwrap());
    }
    (MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(),
        process: "layered thermoelastic test".into(), revision: 0,
    }, claims, vec![]).unwrap(), requirements)
}

fn source<'a>(cards: &[&'a MaterialCard], requirements: &[Vec<ScalarPropertyRequirement>])
    -> ConcentricSectionTemperatureLaw<'a>
{
    let point = QueryPoint::new().with("T", 300.0).unwrap()
        .with("pressure", 100000.0).unwrap();
    let layers = cards.iter().zip(requirements).enumerate().map(|(i, (card, requirements))| {
        ConcentricMaterialLayer {
            outer_radius_m: (i + 1) as f64 * 0.001,
            material: resolve_material_state_point(card, &point, requirements,
                MaterialPropertySelection::SingleClaimOnly).unwrap(),
        }
    }).collect();
    ConcentricSectionTemperatureLaw::try_new(
        ConcentricElasticSection::try_new(0.0, layers).unwrap(), cards,
        &vec!["T"; cards.len()],
    ).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    let scale = actual.abs().max(expected.abs()).max(1e-18);
    assert!((actual - expected).abs() <= relative * scale, "{actual} != {expected}");
}

#[test]
fn g1_unequal_layer_temperatures_share_force_free_energy_and_entropy() {
    let (a, ra) = card(1.0, 1e-5, None);
    let (b, rb) = card(2.0, 2e-5, None);
    let law = ConcentricThermoelasticLaw::try_new(source(&[&a, &b], &[ra, rb])).unwrap();
    let strain = 0.0006;
    let temperatures = [340.0, 370.0];
    let result = law.at_temperatures_k(&temperatures, strain, 0.01).unwrap();
    let response = result.section.response_at_axial_strain(strain, 0.01).unwrap();
    let mut force = 0.0;
    let mut energy = 0.0;
    for (i, layer) in result.layers.iter().enumerate() {
        let scale = (i + 1) as f64;
        let alpha = scale * 1e-5;
        let area = core::f64::consts::PI * ((i + 1).pow(2) - i.pow(2)) as f64 * 1e-6;
        let moment = core::f64::consts::PI * ((i + 1).pow(4) - i.pow(4)) as f64 * 1e-12 / 4.0;
        let e = strain - alpha * (temperatures[i] - 300.0);
        let modulus = scale * (1e9 - 1e6 * (temperatures[i] - 300.0));
        let slope = -scale * 1e6;
        let psi = 0.5 * modulus * area * e * e;
        let entropy = area * (modulus * alpha * e - 0.5 * slope * e * e);
        let curvature = area * (modulus * alpha * alpha - 2.0 * slope * alpha * e);
        close(layer.area_m2, area, 4e-13);
        close(layer.second_moment_m4, moment, 4e-13);
        close(layer.young_modulus_pa, modulus, 4e-13);
        close(layer.young_modulus_derivative_pa_k, slope, 4e-13);
        close(layer.alpha_per_k, alpha, 4e-13);
        close(layer.elastic_strain, e, 4e-13);
        close(layer.axial_free_energy_per_reference_length_j_m, psi, 4e-13);
        close(layer.entropy_per_reference_length_j_m_k, entropy, 4e-13);
        close(layer.free_energy_temperature_second_derivative_j_m_k2, curvature, 4e-13);
        assert_eq!(layer.temperature_k, temperatures[i]);
        assert_eq!(layer.density_kg_m3, scale * 4000.0);
        force += modulus * area * e;
        energy += psi;
        let original = &law.reference().layers()[i].material;
        let current = &result.section.section().layers()[i].material;
        for (old, new) in original.properties().iter().zip(current.properties()) {
            assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
        }
        close(current.property("auxiliary").unwrap().value_si(),
            2.0 + 0.01 * (temperatures[i] - 300.0), 4e-13);
        assert_eq!(current.query().axes()["pressure"], 100000.0);
    }
    assert!(response.layers()[0].force_n() > 0.0);
    assert!(response.layers()[1].force_n() < 0.0);
    close(response.force_n(), force, 4e-13);
    close(response.strain_energy_per_length_j_m(), energy, 4e-13);
}

#[test]
fn g1_temperature_derivatives_and_internal_energy_obey_helmholtz_relations() {
    // Negative expansion is data, and must not be clamped to zero.
    let (card, requirements) = card(1.0, -1.5e-5, None);
    let law = ConcentricThermoelasticLaw::try_new(source(&[&card], &[requirements])).unwrap();
    let temperature = 345.0;
    let delta = 0.01;
    let strain = 0.002;
    let low = law.at_temperatures_k(&[temperature - delta], strain, 0.01).unwrap();
    let mid = law.at_temperatures_k(&[temperature], strain, 0.01).unwrap();
    let high = law.at_temperatures_k(&[temperature + delta], strain, 0.01).unwrap();
    let (low, mid, high) = (&low.layers[0], &mid.layers[0], &high.layers[0]);
    let derivative = (high.axial_free_energy_per_reference_length_j_m
        - low.axial_free_energy_per_reference_length_j_m) / (2.0 * delta);
    close(-derivative, mid.entropy_per_reference_length_j_m_k, 2e-8);
    let entropy_derivative = (high.entropy_per_reference_length_j_m_k
        - low.entropy_per_reference_length_j_m_k) / (2.0 * delta);
    close(-entropy_derivative, mid.free_energy_temperature_second_derivative_j_m_k2, 2e-8);
    let internal = |layer: &LayerThermoelasticState| {
        layer.axial_free_energy_per_reference_length_j_m
            + layer.temperature_k * layer.entropy_per_reference_length_j_m_k
    };
    close((internal(high) - internal(low)) / (2.0 * delta),
        -temperature * mid.free_energy_temperature_second_derivative_j_m_k2, 2e-8);
    let reference = law.at_temperatures_k(&[300.0], 0.0, 0.01).unwrap();
    assert_eq!(reference.layers[0].axial_free_energy_per_reference_length_j_m, 0.0);
    assert_eq!(reference.layers[0].entropy_per_reference_length_j_m_k, 0.0);
}

#[test]
fn g0_nonaffine_or_variable_reference_laws_and_invalid_states_refuse() {
    for defect in [YOUNG_MODULUS_PROPERTY, ALPHA, DENSITY_PROPERTY] {
        let (card, requirements) = card(1.0, 1e-5, Some(defect));
        assert!(ConcentricThermoelasticLaw::try_new(source(&[&card], &[requirements])).is_err(),
            "unsupported law admitted: {defect}");
    }
    let (card, requirements) = card(1.0, 1e-5, None);
    let law = ConcentricThermoelasticLaw::try_new(source(&[&card], &[requirements])).unwrap();
    let reference = law.reference().clone();
    for temperatures in [vec![], vec![300.0, 300.0], vec![279.0], vec![421.0], vec![f64::NAN]] {
        assert!(law.at_temperatures_k(&temperatures, 0.0, 0.01).is_err());
    }
    assert!(law.at_temperatures_k(&[340.0], 0.02, 0.01).is_err());
    assert!(law.at_temperatures_k(&[340.0], 0.0, 1e-4).is_err());
    assert!(law.at_temperatures_k(&[340.0], f64::NAN, 0.01).is_err());
    assert_eq!(law.reference(), &reference);
    let recovered = law.at_temperatures_k(&[300.0], 0.0, 0.01).unwrap();
    assert_eq!(recovered.section.section(), &reference);
}
