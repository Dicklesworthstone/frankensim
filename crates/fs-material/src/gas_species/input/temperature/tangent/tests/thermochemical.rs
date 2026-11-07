use super::*;
use crate::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use crate::thermochemical_gas::tests::{point, species_in_domain};
use crate::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, DynViscosity, QuantitySpec};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

mod surface;

fn curve(knots: &[(f64, f64)], dims: Dims) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: knots.to_vec(),
        dims,
    }
}

fn source(
    name: &str,
    mu: PropertyValue,
    k: PropertyValue,
    interpolation: InterpolationPolicy,
    bounds: Option<[f64; 2]>,
) -> SpeciesTransportMaterial {
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        "pressure",
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
        50000.0,
        150000.0,
    );
    if let Some([lo, hi]) = bounds {
        validity = validity.with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            lo,
            hi,
        );
    }
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        ("dynamic_viscosity", mu, DynViscosity::DIMS),
        ("thermal_conductivity", k, THERMAL_CONDUCTIVITY_DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value,
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured transport".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "manufactured".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

fn gas(transport: Vec<SpeciesTransportMaterial>, fractions: Vec<f64>) -> ThermochemicalGasModel {
    let components = transport
        .iter()
        .map(|transport| {
            let (mass, cp, slope) = if transport.species.as_str() == "A" {
                (0.028, 3.5, 0.001)
            } else {
                (0.004, 2.5, 0.002)
            };
            species_in_domain(
                transport.species.as_str(),
                mass,
                cp,
                slope,
                2e7,
                [200.0, 700.0],
            )
        })
        .collect();
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            components,
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap()
}

#[test]
fn variable_pure_slopes_use_actual_wilke_wms_interactions_and_nasa_calorics() {
    let a = source(
        "A",
        curve(&[(250.0, 1.8e-5), (650.0, 3.4e-5)], DynViscosity::DIMS),
        curve(
            &[(250.0, 0.0275), (650.0, 0.0475)],
            THERMAL_CONDUCTIVITY_DIMS,
        ),
        InterpolationPolicy::LinearInside,
        Some([250.0, 650.0]),
    );
    let b = source(
        "B",
        curve(&[(250.0, 9.5e-6), (650.0, 1.35e-5)], DynViscosity::DIMS),
        curve(&[(250.0, 0.11), (650.0, 0.19)], THERMAL_CONDUCTIVITY_DIMS),
        InterpolationPolicy::LinearInside,
        Some([250.0, 650.0]),
    );
    let model = gas(vec![a.clone(), b.clone()], vec![0.75, 0.25]);
    let law = IsobaricGasLaw::try_from_model(model.clone(), point(300.0, 100000.0)).unwrap();
    let permuted =
        IsobaricGasLaw::try_from_model(gas(vec![b, a], vec![0.25, 0.75]), point(300.0, 100000.0))
            .unwrap();
    for t in [275.0, 330.0, 500.0, 625.0] {
        let response = law
            .transport_temperature_tangent(Temperature::new(t))
            .unwrap();
        let h = 0.01;
        let samples = [-2.0, -1.0, 1.0, 2.0]
            .map(|scale| values(&law.at_temperature(Temperature::new(t + scale * h)).unwrap()));
        for (i, exact) in derivatives(&response).iter().enumerate() {
            let finite_difference = (samples[0][i] - 8.0 * samples[1][i] + 8.0 * samples[2][i]
                - samples[3][i])
                / (12.0 * h);
            near(*exact, finite_difference, 3e-8);
        }
        near(
            response.specific_heat_derivative_j_kg_k2,
            R * (0.75 * 0.001 + 0.25 * 0.002) / (0.75 * 0.028 + 0.25 * 0.004),
            3e-14,
        );
        assert_eq!(
            response.gas,
            law.at_temperature(Temperature::new(t)).unwrap()
        );
        assert_eq!(
            response,
            permuted
                .transport_temperature_tangent(Temperature::new(t))
                .unwrap()
        );
    }
    let high_pressure = IsobaricGasLaw::try_from_model(model, point(300.0, 150000.0)).unwrap();
    let a = law
        .transport_temperature_tangent(Temperature::new(400.0))
        .unwrap();
    let b = high_pressure
        .transport_temperature_tangent(Temperature::new(400.0))
        .unwrap();
    near(
        1.5 * a.density_derivative_kg_m3_k,
        b.density_derivative_kg_m3_k,
        3e-14,
    );
    assert_eq!(
        a.specific_heat_derivative_j_kg_k2,
        b.specific_heat_derivative_j_kg_k2
    );
    assert_eq!(
        a.conductivity_derivative_w_m_k2,
        b.conductivity_derivative_w_m_k2
    );
    assert_eq!(
        a.dynamic_viscosity_derivative_pa_s_k,
        b.dynamic_viscosity_derivative_pa_s_k
    );
}

#[test]
fn source_corner_refuses_but_validity_endpoints_use_their_admitted_side() {
    let bind = |bounds| {
        IsobaricGasLaw::try_from_model(
            gas(
                vec![source(
                    "A",
                    curve(
                        &[(250.0, 1e-5), (300.0, 2e-5), (400.0, 5e-5)],
                        DynViscosity::DIMS,
                    ),
                    curve(&[(250.0, 0.02), (400.0, 0.05)], THERMAL_CONDUCTIVITY_DIMS),
                    InterpolationPolicy::LinearInside,
                    Some(bounds),
                )],
                vec![1.0],
            ),
            point(300.0, 100000.0),
        )
        .unwrap()
    };
    let complete = bind([250.0, 400.0]);
    complete.at_temperature(Temperature::new(300.0)).unwrap();
    assert!(
        matches!(complete.transport_temperature_tangent(Temperature::new(300.0)),
        Err(MaterialError::State { what }) if what.contains("no unique temperature derivative"))
    );
    let left = bind([250.0, 300.0])
        .transport_temperature_tangent(Temperature::new(300.0))
        .unwrap();
    let right = bind([300.0, 400.0])
        .transport_temperature_tangent(Temperature::new(300.0))
        .unwrap();
    near(left.dynamic_viscosity_derivative_pa_s_k, 2e-7, 3e-14);
    near(right.dynamic_viscosity_derivative_pa_s_k, 3e-7, 3e-14);
    assert_eq!(left.gas.state(), right.gas.state());
}

#[test]
fn explicit_constant_source_support_supplies_zero_tangents_without_promoting_discrete_data() {
    let bind = |interpolation, bounds| {
        let (mu, k) = if interpolation == InterpolationPolicy::ConstantWithinValidity {
            (
                PropertyValue::Scalar {
                    value: 2e-5,
                    dims: DynViscosity::DIMS,
                },
                PropertyValue::Scalar {
                    value: 0.03,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                },
            )
        } else {
            (
                curve(&[(300.0, 2e-5), (400.0, 3e-5)], DynViscosity::DIMS),
                curve(&[(300.0, 0.03), (400.0, 0.04)], THERMAL_CONDUCTIVITY_DIMS),
            )
        };
        IsobaricGasLaw::try_from_model(
            gas(vec![source("A", mu, k, interpolation, bounds)], vec![1.0]),
            point(300.0, 100000.0),
        )
        .unwrap()
    };
    let admitted = bind(
        InterpolationPolicy::ConstantWithinValidity,
        Some([250.0, 650.0]),
    );
    let response = admitted
        .transport_temperature_tangent(Temperature::new(315.0))
        .unwrap();
    assert_eq!(response.dynamic_viscosity_derivative_pa_s_k, 0.0);
    assert_eq!(response.conductivity_derivative_w_m_k2, 0.0);
    assert!(response.specific_heat_derivative_j_kg_k2 > 0.0);
    for law in [
        bind(InterpolationPolicy::ConstantWithinValidity, None),
        bind(InterpolationPolicy::TabulatedOnly, Some([250.0, 650.0])),
    ] {
        law.at_temperature(Temperature::new(300.0)).unwrap();
        assert!(
            law.transport_temperature_tangent(Temperature::new(300.0))
                .is_err()
        );
    }
}
