//! Independent kinetic-theory limits through the original NASA/material owners.
use super::*;
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    tests::{
        model, point, species, species_in_domain, species_with_energy_reference, transport,
        transport_in_domain,
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const AREA: f64 = 1e-10;
fn support() -> EffusiveSupport {
    EffusiveSupport {
        temperature_k: (200.0, 2000.0),
        pressure_pa: (1000.0, 1e7),
    }
}
fn regime() -> EffusiveRegime {
    EffusiveRegime::try_new(
        1e-9,
        1e-12,
        1e-7,
        "synthetic qualification for every mixture of the retained test gas catalog",
    )
    .unwrap()
}
fn interface(gas: ThermochemicalGasModel, probability: f64) -> EffusiveInterface {
    let rows = gas
        .thermodynamics()
        .components()
        .iter()
        .map(|s| {
            EffusiveSpeciesTransmission::try_new(
                s.species().clone(),
                probability,
                "synthetic reciprocal energy/internal-state-independent transmission",
            )
            .unwrap()
        })
        .collect();
    EffusiveInterface::try_new(gas, rows, support(), regime()).unwrap()
}
fn pure_model(
    name: &str,
    mass: f64,
    cp_r: f64,
    slope: f64,
    reference: f64,
) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species_with_energy_reference(
                name, mass, cp_r, slope, reference,
            )],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport(name, 1.7e-5, 0.025)],
        None,
    )
    .unwrap()
}
fn pure(mass: f64, cp_r: f64, slope: f64, reference: f64, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        pure_model("A", mass, cp_r, slope, reference),
        point(t, p),
        0.001,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-12 * actual.abs().max(expected.abs()).max(1e-100),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn one_way(area: f64, transmission: f64, partial_pressure: f64, mass: f64, t: f64) -> f64 {
    area * transmission * partial_pressure / (2.0 * core::f64::consts::PI * mass * R * t).sqrt()
}

#[test]
fn maxwellian_beams_recover_monatomic_and_polyatomic_molecular_energies() {
    for (mass, cp_r, slope) in [(0.004, 2.5, 0.0), (0.028, 3.5, 0.001)] {
        let (left, right) = (
            pure(mass, cp_r, slope, 0.0, 400.0, 110_000.0),
            pure(mass, cp_r, slope, 0.0, 900.0, 80_000.0),
        );
        let law = interface(left.model().clone(), 0.65);
        let flux = left.effusive_flux_to(&right, &law, AREA, || false).unwrap();
        assert_eq!(flux.species.len(), 1);
        let row = &flux.species[0];
        for (beam, t, p) in [
            (row.left_to_right, 400.0, 110_000.0),
            (row.right_to_left, 900.0, 80_000.0),
        ] {
            let n = one_way(AREA, 0.65, p, mass, t);
            let h = R * (cp_r * t + 0.5 * slope * t * t) / mass;
            let crossing = R * ((cp_r - 0.5) * t + 0.5 * slope * t * t) / mass;
            near(beam.molar_rate_mol_s, n);
            near(beam.mass_rate_kg_s, n * mass);
            near(beam.source_enthalpy_j_kg.unwrap(), h);
            near(beam.crossing_energy_j_kg.unwrap(), crossing);
            near(
                beam.additional_specific_energy_j_kg.unwrap(),
                -R * t / (2.0 * mass),
            );
            near(beam.energy_rate_w, n * mass * crossing);
            near(beam.additional_energy_rate_w, -n * R * t / 2.0);
        }
        near(
            flux.energy_rate_left_to_right_w,
            row.left_to_right.energy_rate_w - row.right_to_left.energy_rate_w,
        );
        near(
            flux.additional_energy_rate_left_to_right_w,
            row.left_to_right.additional_energy_rate_w - row.right_to_left.additional_energy_rate_w,
        );
    }
}

#[test]
fn different_molecular_masses_reciprocate_with_endpoint_and_transmission_order() {
    let left =
        ClosedGasParcel::new(model(vec![0.75, 0.25], 0.0), point(400.0, 100_000.0), 0.001).unwrap();
    let right =
        ClosedGasParcel::new(model(vec![0.25, 0.75], 0.0), point(900.0, 120_000.0), 0.002).unwrap();
    let rows = [("B", 0.3), ("A", 0.6)]
        .into_iter()
        .map(|(id, p)| {
            EffusiveSpeciesTransmission::try_new(
                SpeciesId::new(id).unwrap(),
                p,
                "synthetic reciprocal species selectivity",
            )
            .unwrap()
        })
        .collect();
    let law =
        EffusiveInterface::try_new(model(vec![0.5, 0.5], 0.0), rows, support(), regime()).unwrap();
    let flux = left.effusive_flux_to(&right, &law, AREA, || false).unwrap();
    let reverse = right.effusive_flux_to(&left, &law, AREA, || false).unwrap();
    assert_eq!(flux.species.len(), 2);
    for (i, (mass, tau, xl, xr)) in [(0.028, 0.6, 0.75, 0.25), (0.004, 0.3, 0.25, 0.75)]
        .into_iter()
        .enumerate()
    {
        near(
            flux.species[i].left_to_right.molar_rate_mol_s,
            one_way(AREA, tau, xl * 100_000.0, mass, 400.0),
        );
        near(
            flux.species[i].right_to_left.molar_rate_mol_s,
            one_way(AREA, tau, xr * 120_000.0, mass, 900.0),
        );
        assert_eq!(
            flux.species[i].left_to_right,
            reverse.species[i].right_to_left
        );
        assert_eq!(
            flux.species[i].right_to_left,
            reverse.species[i].left_to_right
        );
    }
    near(
        reverse.energy_rate_left_to_right_w,
        -flux.energy_rate_left_to_right_w,
    );
    near(
        reverse.additional_energy_rate_left_to_right_w,
        -flux.additional_energy_rate_left_to_right_w,
    );
    near(
        reverse.net_mass_rate_left_to_right_kg_s,
        -flux.net_mass_rate_left_to_right_kg_s,
    );
    let mut ordered = law.transmissions().to_vec();
    ordered.reverse();
    let permuted =
        EffusiveInterface::try_new(law.gas_catalog().clone(), ordered, support(), regime())
            .unwrap();
    assert_eq!(
        flux,
        left.effusive_flux_to(&right, &permuted, AREA, || false)
            .unwrap()
    );
}

#[test]
fn zero_net_mass_thermal_transpiration_keeps_heat_under_huge_energy_reference_shifts() {
    let mut reference_flux = None;
    for reference in [0.0, 1e24, -1e24] {
        // Square-root temperature ratio2 and pressure ratio2 give equal gross
        // rates exactly, avoiding a numerical near-equilibrium ambiguity.
        let left = pure(0.004, 2.5, 0.0, reference, 400.0, 100_000.0);
        let right = pure(0.004, 2.5, 0.0, reference, 1600.0, 200_000.0);
        let law = interface(left.model().clone(), 1.0);
        let flux = left.effusive_flux_to(&right, &law, AREA, || false).unwrap();
        let l = flux.species[0].left_to_right;
        let r = flux.species[0].right_to_left;
        assert_eq!(l.molar_rate_mol_s, r.molar_rate_mol_s);
        assert_eq!(l.mass_rate_kg_s, r.mass_rate_kg_s);
        assert_eq!(flux.net_molar_rate_left_to_right_mol_s, 0.0);
        assert_eq!(flux.net_mass_rate_left_to_right_kg_s, 0.0);
        let expected = l.molar_rate_mol_s * 2.0 * R * (400.0 - 1600.0);
        near(flux.energy_rate_left_to_right_w, expected);
        near(
            flux.additional_energy_rate_left_to_right_w,
            -l.molar_rate_mol_s * R * (400.0 - 1600.0) / 2.0,
        );
        assert!(flux.energy_rate_left_to_right_w < 0.0);
        assert!(flux.energy_rate_left_to_right_w * (1.0 / 1600.0 - 1.0 / 400.0) > 0.0);
        if let Some(previous) = reference_flux {
            near(flux.energy_rate_left_to_right_w, previous);
        }
        reference_flux = Some(flux.energy_rate_left_to_right_w);
        if reference < 0.0 {
            assert!(l.crossing_energy_j_kg.unwrap() < 0.0);
            assert!(l.energy_rate_w < 0.0);
        }
    }
}

#[test]
fn only_zero_area_bypasses_inactive_support_and_cancellation_is_observed() {
    let reference = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let foreign = pure(0.004, 2.5, 0.0, 0.0, 1000.0, 100_000.0);
    let law = EffusiveInterface::try_new(
        reference.model().clone(),
        vec![
            EffusiveSpeciesTransmission::try_new(
                SpeciesId::new("A").unwrap(),
                0.0,
                "synthetic sealed species only inside qualification",
            )
            .unwrap(),
        ],
        EffusiveSupport {
            temperature_k: (300.0, 500.0),
            pressure_pa: (90_000.0, 110_000.0),
        },
        regime(),
    )
    .unwrap();
    let closed = foreign
        .effusive_flux_to(&foreign, &law, 0.0, || false)
        .unwrap();
    assert!(closed.species.is_empty());
    assert_eq!(closed.energy_rate_left_to_right_w, 0.0);
    assert_eq!(closed.additional_energy_rate_left_to_right_w, 0.0);
    assert!(
        foreign
            .effusive_flux_to(&foreign, &law, AREA, || false)
            .is_err()
    );
    let admitted = reference
        .effusive_flux_to(&reference, &law, AREA, || false)
        .unwrap();
    assert_eq!(admitted.species[0].left_to_right, EffusiveBeam::default());
    assert_eq!(admitted.species[0].right_to_left, EffusiveBeam::default());
    assert!(matches!(
        reference.effusive_flux_to(&reference, &law, 0.0, || true),
        Err(GasEnergyError::Cancelled)
    ));
    let mut calls = 0;
    assert!(matches!(
        reference.effusive_flux_to(&reference, &law, AREA, || {
            calls += 1;
            calls >= 3
        }),
        Err(GasEnergyError::Cancelled)
    ));
}

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn transmission_source(
    value: PropertyValue,
    dims: Dims,
    policy: InterpolationPolicy,
    extra_axis: bool,
) -> ResolvedMaterialStatePoint {
    let mut validity = ValidityDomain::unconstrained()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            300.0,
            800.0,
        )
        .with_quantity(
            "pressure",
            typed(QuantityKind::Pressure),
            50_000.0,
            200_000.0,
        );
    let mut query = point(400.0, 100_000.0);
    if extra_axis {
        validity = validity.with("strain", 0.0, 1.0);
        query = query.with("strain", 0.2).unwrap();
    }
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(EFFUSIVE_TRANSMISSION_PROPERTY, dims),
            value,
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: policy,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic calibrated reciprocal transmission".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic A".into(),
            phase: "solid aperture".into(),
            process: "synthetic energy-independent transmission fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirement = ScalarPropertyRequirement::try_new(
        EFFUSIVE_TRANSMISSION_PROPERTY,
        dims,
        ScalarAdmissibility::NonNegative,
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &query,
        &[requirement],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn constant_source(probability: f64) -> ResolvedMaterialStatePoint {
    transmission_source(
        PropertyValue::Scalar {
            value: probability,
            dims: EFFUSIVE_TRANSMISSION_DIMS,
        },
        EFFUSIVE_TRANSMISSION_DIMS,
        InterpolationPolicy::ConstantWithinValidity,
        false,
    )
}

#[test]
fn material_receipts_and_actual_endpoint_support_are_retained_without_freezing_axes() {
    let source = constant_source(0.4);
    let row = EffusiveSpeciesTransmission::from_material(
        SpeciesId::new("A").unwrap(),
        source.clone(),
        EFFUSIVE_TRANSMISSION_PROPERTY,
        "synthetic calibration: reciprocal, independent of energy and internal state",
    )
    .unwrap();
    assert_eq!(row.source(), Some(&source));
    assert_eq!(row.source_property(), Some(EFFUSIVE_TRANSMISSION_PROPERTY));
    assert_eq!(row.probability(), 0.4);
    let left = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let right = pure(0.028, 3.5, 0.0, 0.0, 700.0, 100_000.0);
    let law =
        EffusiveInterface::try_new(left.model().clone(), vec![row], support(), regime()).unwrap();
    let flux = left.effusive_flux_to(&right, &law, AREA, || false).unwrap();
    assert_eq!(
        flux.interface.transmissions()[0]
            .source()
            .unwrap()
            .identity(),
        source.identity()
    );
    for (t, p) in [(900.0, 100_000.0), (400.0, 300_000.0)] {
        let outside = pure(0.028, 3.5, 0.0, 0.0, t, p);
        assert!(
            left.effusive_flux_to(&outside, &law, AREA, || false)
                .is_err()
        );
    }
    let closed = EffusiveSpeciesTransmission::from_material(
        SpeciesId::new("A").unwrap(),
        constant_source(0.0),
        EFFUSIVE_TRANSMISSION_PROPERTY,
        "calibrated zero in its source box",
    )
    .unwrap();
    let closed_law =
        EffusiveInterface::try_new(left.model().clone(), vec![closed], support(), regime())
            .unwrap();
    let outside = pure(0.028, 3.5, 0.0, 0.0, 900.0, 100_000.0);
    assert!(
        left.effusive_flux_to(&outside, &closed_law, AREA, || false)
            .is_err()
    );

    let wrong_dims = Dims([2, 0, 0, 0, 0, 0]);
    let wrong = transmission_source(
        PropertyValue::Scalar {
            value: 0.5,
            dims: wrong_dims,
        },
        wrong_dims,
        InterpolationPolicy::ConstantWithinValidity,
        false,
    );
    assert!(
        EffusiveSpeciesTransmission::from_material(
            SpeciesId::new("A").unwrap(),
            wrong,
            EFFUSIVE_TRANSMISSION_PROPERTY,
            "wrong intrinsic-permeability dimensions"
        )
        .is_err()
    );
    let varying = transmission_source(
        PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 0.4), (800.0, 0.6)],
            dims: EFFUSIVE_TRANSMISSION_DIMS,
        },
        EFFUSIVE_TRANSMISSION_DIMS,
        InterpolationPolicy::LinearInside,
        false,
    );
    assert!(
        EffusiveSpeciesTransmission::from_material(
            SpeciesId::new("A").unwrap(),
            varying,
            EFFUSIVE_TRANSMISSION_PROPERTY,
            "temperature-dependent transmission is not frozen"
        )
        .is_err()
    );
    let strained = transmission_source(
        PropertyValue::Scalar {
            value: 0.5,
            dims: EFFUSIVE_TRANSMISSION_DIMS,
        },
        EFFUSIVE_TRANSMISSION_DIMS,
        InterpolationPolicy::ConstantWithinValidity,
        true,
    );
    assert!(
        EffusiveSpeciesTransmission::from_material(
            SpeciesId::new("A").unwrap(),
            strained,
            EFFUSIVE_TRANSMISSION_PROPERTY,
            "strain requires a separate evolving owner"
        )
        .is_err()
    );
    assert!(
        EffusiveSpeciesTransmission::from_material(
            SpeciesId::new("A").unwrap(),
            constant_source(1.01),
            EFFUSIVE_TRANSMISSION_PROPERTY,
            "probability must not exceed one"
        )
        .is_err()
    );
}

#[test]
fn original_gas_catalog_and_physical_qualification_refuse_unsupported_substitution() {
    let left = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let law = interface(left.model().clone(), 1.0);
    let other_reference = pure(0.028, 3.5, 0.0, 1e6, 400.0, 100_000.0);
    assert!(
        other_reference
            .effusive_flux_to(&other_reference, &law, AREA, || false)
            .is_err()
    );
    let other_mass = pure(0.004, 2.5, 0.0, 0.0, 400.0, 100_000.0);
    assert!(
        other_mass
            .effusive_flux_to(&other_mass, &law, AREA, || false)
            .is_err()
    );
    assert!(EffusiveRegime::try_new(1e-6, 1e-10, 1e-7, "collisional").is_err());
    assert!(EffusiveRegime::try_new(1e-9, 1e-9, 1e-6, "long pore").is_err());
    assert!(EffusiveRegime::try_new(1e-9, 0.0, 1e-6, "").is_err());
    assert!(
        EffusiveSpeciesTransmission::try_new(SpeciesId::new("A").unwrap(), -0.1, "invalid")
            .is_err()
    );
    assert!(
        EffusiveSpeciesTransmission::try_new(SpeciesId::new("A").unwrap(), f64::NAN, "invalid")
            .is_err()
    );
    let rows = vec![
        law.transmissions()[0].clone(),
        law.transmissions()[0].clone(),
    ];
    assert!(EffusiveInterface::try_new(left.model().clone(), rows, support(), regime()).is_err());
    let invalid = EffusiveSupport {
        temperature_k: (500.0, 300.0),
        ..support()
    };
    assert!(
        EffusiveInterface::try_new(
            left.model().clone(),
            law.transmissions().to_vec(),
            invalid,
            regime()
        )
        .is_err()
    );
    // A larger explicitly qualified source catalog admits its actual subsets.
    let larger = interface(model(vec![0.5, 0.5], 0.0), 0.5);
    assert!(
        left.effusive_flux_to(&left, &larger, AREA, || false)
            .is_ok()
    );
}

#[test]
fn absent_side_species_never_query_the_other_gas_caloric_temperature() {
    let a = species_in_domain("A", 0.028, 3.5, 0.0, 0.0, [200.0, 500.0]);
    let b = species_in_domain("B", 0.004, 2.5, 0.0, 0.0, [800.0, 1800.0]);
    let ta = transport_in_domain("A", 1.7e-5, 0.025, [200.0, 500.0], [1000.0, 1e7]);
    let tb = transport_in_domain("B", 2e-5, 0.15, [800.0, 1800.0], [1000.0, 1e7]);
    let make = |s, t| {
        ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![s],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            )
            .unwrap(),
            vec![t],
            None,
        )
        .unwrap()
    };
    let left =
        ClosedGasParcel::new(make(a.clone(), ta.clone()), point(400.0, 100_000.0), 0.001).unwrap();
    let right =
        ClosedGasParcel::new(make(b.clone(), tb.clone()), point(1000.0, 100_000.0), 0.001).unwrap();
    let bound = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![ta, tb],
        None,
    )
    .unwrap();
    let flux = left
        .effusive_flux_to(&right, &interface(bound, 1.0), AREA, || false)
        .unwrap();
    assert!(flux.species[0].left_to_right.mass_rate_kg_s > 0.0);
    assert_eq!(flux.species[0].right_to_left, EffusiveBeam::default());
    assert_eq!(flux.species[1].left_to_right, EffusiveBeam::default());
    assert!(flux.species[1].right_to_left.mass_rate_kg_s > 0.0);
}

#[test]
fn positive_tiny_beams_are_not_closed_by_intermediate_underflow() {
    let gas = pure(0.028, 3.5, 0.0, 0.0, 400.0, 100_000.0);
    let law = interface(gas.model().clone(), 0.5);
    // area*tau alone rounds to zero; multiplication by the physical incident
    // number density makes the final beam representable and must occur first
    // in a scaled evaluation rather than silently closing the aperture.
    let tiny_area = f64::from_bits(1);
    assert_eq!(tiny_area * 0.5, 0.0);
    let tiny = gas
        .effusive_flux_to(&gas, &law, tiny_area, || false)
        .unwrap();
    assert!(tiny.species[0].left_to_right.molar_rate_mol_s > 0.0);
    assert!(tiny.species[0].left_to_right.mass_rate_kg_s > 0.0);
    assert!(tiny.species[0].left_to_right.additional_energy_rate_w < 0.0);
    assert_eq!(tiny.net_mass_rate_left_to_right_kg_s, 0.0);
    assert_eq!(tiny.energy_rate_left_to_right_w, 0.0);
    let unrepresentable = interface(gas.model().clone(), 1e-300);
    assert!(
        gas.effusive_flux_to(&gas, &unrepresentable, 1e-300, || false)
            .is_err()
    );
    assert!(
        gas.effusive_flux_to(&gas, &law, f64::MAX, || false)
            .is_err()
    );
    assert!(gas.effusive_flux_to(&gas, &law, -1.0, || false).is_err());
}

#[test]
fn ternary_calibration_covers_shifted_nonbinary_fractions_without_renormalizing_sources() {
    let make = |fractions| {
        ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![
                    species("A", 0.028, 3.5, 0.0),
                    species("B", 0.004, 2.5, 0.0),
                    species("C", 0.044, 4.0, 0.0),
                ],
                Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
            )
            .unwrap(),
            vec![
                transport("A", 1.7e-5, 0.025),
                transport("B", 2e-5, 0.15),
                transport("C", 1.5e-5, 0.02),
            ],
            None,
        )
        .unwrap()
    };
    let catalog = make(vec![0.1, 0.2, 0.7]);
    let law = interface(catalog.clone(), 0.8);
    let left = ClosedGasParcel::new(catalog, point(400.0, 100_000.0), 0.001).unwrap();
    let right =
        ClosedGasParcel::new(make(vec![0.7, 0.1, 0.2]), point(900.0, 140_000.0), 0.003).unwrap();
    let flux = left.effusive_flux_to(&right, &law, AREA, || false).unwrap();
    assert_eq!(flux.species.len(), 3);
    for (i, (molar_mass, xl, xr)) in [(0.028, 0.1, 0.7), (0.004, 0.2, 0.1), (0.044, 0.7, 0.2)]
        .into_iter()
        .enumerate()
    {
        near(
            flux.species[i].left_to_right.molar_rate_mol_s,
            one_way(AREA, 0.8, xl * 100_000.0, molar_mass, 400.0),
        );
        near(
            flux.species[i].right_to_left.molar_rate_mol_s,
            one_way(AREA, 0.8, xr * 140_000.0, molar_mass, 900.0),
        );
    }
    let pure_c = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("C", 0.044, 4.0, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport("C", 1.5e-5, 0.02)],
        None,
    )
    .unwrap();
    let subset = ClosedGasParcel::new(pure_c, point(500.0, 90_000.0), 0.001).unwrap();
    let mixed = subset
        .effusive_flux_to(&right, &law, AREA, || false)
        .unwrap();
    assert_eq!(mixed.species[0].left_to_right, EffusiveBeam::default());
    assert_eq!(mixed.species[1].left_to_right, EffusiveBeam::default());
    assert!(mixed.species[2].left_to_right.mass_rate_kg_s > 0.0);
}
