//! Manufactured constitutive/source checks, not measured metal properties.
use super::*;
use crate::conductor::ElectricalPropertyBasis;
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, Temperature};

fn card(rho: f64, liquid: bool, basis: LiquidViscosityBasis) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let viscosity = match basis {
        LiquidViscosityBasis::Dynamic => (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 0.02),
        LiquidViscosityBasis::Kinematic => {
            (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS, 2e-5)
        }
    };
    let electrical = ElectricalPropertyBasis::Resistivity;
    for (name, dims, value, liquid_only) in [
        (DENSITY_PROPERTY, Density::DIMS, rho, false),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            Dims([1, 1, -3, -1, 0, 0]),
            10.0,
            false,
        ),
        (
            electrical.property_name(),
            electrical.dimensions(),
            1e-6,
            false,
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            3000.0,
            true,
        ),
        (viscosity.0, viscosity.1, viscosity.2, true),
    ] {
        if liquid_only && !liquid {
            continue;
        }
        let conductivity = name == THERMAL_CONDUCTIVITY_PROPERTY;
        let quantity = if conductivity {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        let varying = conductivity || name == electrical.property_name();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: if varying {
                    PropertyValue::Curve {
                        abscissa: "liquid".into(),
                        abscissa_dims: Dims([0; 6]),
                        knots: vec![(0.0, value), (1.0, if conductivity { 2.0 } else { 3e-6 })],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                // Validity axes are typed exactly as point() types them (claims and
                // queries must agree on each axis quantity since 01ac273e9).
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        if liquid_only { 500.0 } else { 400.0 },
                        600.0,
                    )
                    // The smallest positive bound excludes exactly solid, without
                    // inventing a finite cutoff for the material source.
                    .with_quantity(
                        "liquid",
                        QuantitySpec::dimensional(Dims([0; 6])),
                        if liquid_only { f64::from_bits(1) } else { 0.0 },
                        1.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if varying {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "MANUFACTURED phase transport fixture; not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured equal-density phase material".into(),
            phase: "solid-liquid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn curve(card: &MaterialCard, rho: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 400.0, 0.0),
            (0.0, 500.0, 0.0),
            (100.0, 500.0, 1.0),
            (200.0, 600.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, phi)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: phi,
            bulk_density_kg_m3: rho,
        })
        .collect(),
    )
    .unwrap()
}

fn point() -> QueryPoint {
    // Deliberately stale thermal coordinates: only the ORIGINAL h chart may
    // choose the two values that reach the source query.
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            550.0,
        )
        .unwrap()
        .with_quantity("liquid", QuantitySpec::dimensional(Dims([0; 6])), 0.9)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100000.0,
        )
        .unwrap()
}

fn resolve(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
    basis: LiquidViscosityBasis,
) -> Result<ResolvedPhaseTransportState, MaterialError> {
    let point = point();
    let extras = [ElectricalPropertyBasis::Resistivity.requirement().unwrap()];
    resolve_phase_transport_state(
        law,
        h,
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material: card,
                context: &point,
                temperature_axis: "T",
                liquid_fraction_axis: Some("liquid"),
            },
            viscosity_basis: basis,
            extra_requirements: &extras,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    )
}

#[test]
fn exact_solid_keeps_thermal_and_electrical_sources_without_any_liquid_claim() {
    let card = card(1000.0, false, LiquidViscosityBasis::Dynamic);
    let curve = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1000.0, 1000.0, 1e-8).unwrap();
    for h in [-100.0, -50.0, 0.0] {
        let solid = resolve(&card, &law, h, LiquidViscosityBasis::Dynamic).unwrap();
        assert_eq!(solid.phase_state().phase(), SolidLiquidPhase::Solid);
        assert!(matches!(
            solid.mobility(),
            PhaseMobility::ImmobileSolid { .. }
        ));
        assert_eq!(solid.mobile_transport(), None);
        assert!(solid.transport_state().is_err());
        assert_eq!(solid.thermal_conductivity_w_m_k().unwrap(), 10.0);
        assert_eq!(solid.material().properties().len(), 3);
        assert!(
            solid
                .material()
                .property(DYNAMIC_VISCOSITY_PROPERTY)
                .is_none()
        );
        assert!(
            solid
                .material()
                .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                .is_none()
        );
        assert_eq!(
            solid
                .material()
                .property(ElectricalPropertyBasis::Resistivity.property_name())
                .unwrap()
                .value_si(),
            1e-6
        );
        assert_eq!(solid.material().card_identity(), card.content_hash());
        assert_eq!(solid.phase_state().phase_curve_identity(), curve.identity());
        assert_eq!(solid.mobility().law_identity(), law.identity());
    }
    // Melting requires the actual missing liquid sources; no fallback is used.
    assert!(resolve(&card, &law, 25.0, LiquidViscosityBasis::Dynamic).is_err());
}

#[test]
fn latent_plateau_uses_one_h_derived_fraction_for_heat_current_and_momentum() {
    let card = card(1000.0, true, LiquidViscosityBasis::Dynamic);
    let curve = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1000.0, 1000.0, 1e-8).unwrap();
    let a = resolve(&card, &law, 25.0, LiquidViscosityBasis::Dynamic).unwrap();
    let b = resolve(&card, &law, 75.0, LiquidViscosityBasis::Dynamic).unwrap();
    for (state, h, phi, k, resistivity) in
        [(&a, 25.0, 0.25, 8.0, 1.5e-6), (&b, 75.0, 0.75, 4.0, 2.5e-6)]
    {
        assert_eq!(state.phase_state().temperature_k(), 500.0);
        assert_eq!(state.phase_state().specific_enthalpy_j_kg(), h);
        assert_eq!(state.phase_state().liquid_mass_fraction(), phi);
        assert_eq!(state.mobility().volumes().liquid, phi);
        assert_eq!(state.thermal_conductivity_w_m_k().unwrap(), k);
        let transport = state.transport_state().unwrap();
        assert_eq!(transport.thermal_conductivity_w_m_k(), k);
        assert_eq!(transport.density_kg_m3(), 1000.0);
        assert_eq!(transport.dynamic_viscosity_pa_s(), 0.02);
        let electrical = state
            .material()
            .property(ElectricalPropertyBasis::Resistivity.property_name())
            .unwrap()
            .value_si();
        assert!((electrical - resistivity).abs() < 1e-20);
        for property in state.material().properties() {
            let receipt = &property.answer().receipt;
            assert!(receipt.query_point.contains(&("T".into(), 500.0)));
            assert!(receipt.query_point.contains(&("liquid".into(), phi)));
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
        }
    }
    assert!(
        a.mobility().mobile().unwrap().inverse_permeability_per_m2
            > b.mobility().mobile().unwrap().inverse_permeability_per_m2
    );
    for property in a.material().properties() {
        assert_eq!(
            property.answer().receipt.selected,
            b.material()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
}

#[test]
fn melting_and_freezing_keep_original_phase_and_source_identity() {
    let card = card(1000.0, true, LiquidViscosityBasis::Dynamic);
    let curve = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1000.0, 1000.0, 1e-8).unwrap();
    let solid_before = resolve(&card, &law, -50.0, LiquidViscosityBasis::Dynamic).unwrap();
    for h in [-50.0, 0.0, 25.0, 100.0, 150.0, 75.0, 0.0, -50.0] {
        let state = resolve(&card, &law, h, LiquidViscosityBasis::Dynamic).unwrap();
        assert_eq!(
            *state.phase_state(),
            curve.state_at_specific_enthalpy(h).unwrap()
        );
        assert_eq!(state.material().card_identity(), card.content_hash());
        assert_eq!(state.phase_state().bulk_density_kg_m3(), 1000.0);
        assert_eq!(state.mobile_transport().is_some(), h > 0.0);
        assert_eq!(
            state,
            resolve(&card, &law, h, LiquidViscosityBasis::Dynamic).unwrap()
        );
    }
    assert_eq!(
        solid_before,
        resolve(&card, &law, -50.0, LiquidViscosityBasis::Dynamic).unwrap()
    );
}

#[test]
fn equal_density_interpolation_preserves_thermal_mass_bits() {
    // A reciprocal round trip changes this density by 1 ULP.
    let rho = 15.0_f64 / 37.0;
    assert_ne!(rho.to_bits(), rho.recip().recip().to_bits());
    let card = card(rho, true, LiquidViscosityBasis::Dynamic);
    let curve = curve(&card, rho);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, rho, rho, 1e-8).unwrap();
    for h in [-50.0, 25.0, 75.0, 150.0] {
        let state = resolve(&card, &law, h, LiquidViscosityBasis::Dynamic).unwrap();
        assert_eq!(
            state.phase_state().bulk_density_kg_m3().to_bits(),
            rho.to_bits()
        );
        assert_eq!(
            state.mobility().volumes().bulk_density_kg_m3.to_bits(),
            rho.to_bits()
        );
    }
}

#[test]
fn kinematic_source_uses_same_density_without_inventing_a_dynamic_claim() {
    let card = card(1000.0, true, LiquidViscosityBasis::Kinematic);
    let curve = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1000.0, 1000.0, 1e-8).unwrap();
    let state = resolve(&card, &law, 50.0, LiquidViscosityBasis::Kinematic).unwrap();
    assert_eq!(
        state.transport_state().unwrap().dynamic_viscosity_pa_s(),
        0.02
    );
    assert!(
        state
            .material()
            .property(KINEMATIC_VISCOSITY_PROPERTY)
            .is_some()
    );
    assert!(
        state
            .material()
            .property(DYNAMIC_VISCOSITY_PROPERTY)
            .is_none()
    );
    assert!(resolve(&card, &law, 50.0, LiquidViscosityBasis::Dynamic).is_err());
    assert!(
        resolve(&card, &law, 0.0, LiquidViscosityBasis::Kinematic)
            .unwrap()
            .material()
            .property(KINEMATIC_VISCOSITY_PROPERTY)
            .is_none()
    );
}

#[test]
fn source_identity_density_coordinate_and_domain_errors_refuse() {
    let card = card(1000.0, true, LiquidViscosityBasis::Dynamic);
    // The chart is still bound to this card, but its declared mass density is
    // incompatible with the property source.
    let wrong_density = curve(&card, 999.0);
    let law =
        MushyPhaseLaw::try_new_with_immobile_solid(&wrong_density, 999.0, 999.0, 1e-8).unwrap();
    assert!(resolve(&card, &law, 50.0, LiquidViscosityBasis::Dynamic).is_err());
    let correct = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&correct, 1000.0, 1000.0, 1e-8).unwrap();
    let other = self::card(1001.0, true, LiquidViscosityBasis::Dynamic);
    assert!(resolve(&other, &law, 50.0, LiquidViscosityBasis::Dynamic).is_err());
    for h in [-101.0, 201.0, f64::NAN, f64::INFINITY] {
        assert!(resolve(&card, &law, h, LiquidViscosityBasis::Dynamic).is_err());
    }
    let point = point();
    // A source dependent on phase cannot silently retain an unrelated fraction
    // from context when the caller did not bind that coordinate to h.
    assert!(
        resolve_phase_transport_state(
            &law,
            50.0,
            PhaseTransportSourceQuery {
                source: ThermalMaterialPoint {
                    material: &card,
                    context: &point,
                    temperature_axis: "T",
                    liquid_fraction_axis: None
                },
                viscosity_basis: LiquidViscosityBasis::Dynamic,
                extra_requirements: &[],
                selection: MaterialPropertySelection::SingleClaimOnly,
            }
        )
        .is_err()
    );
    let bad_temperature = point
        .clone()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            500.0,
        )
        .unwrap();
    assert!(
        resolve_phase_transport_state(
            &law,
            50.0,
            PhaseTransportSourceQuery {
                source: ThermalMaterialPoint {
                    material: &card,
                    context: &bad_temperature,
                    temperature_axis: "T",
                    liquid_fraction_axis: Some("liquid")
                },
                viscosity_basis: LiquidViscosityBasis::Dynamic,
                extra_requirements: &[],
                selection: MaterialPropertySelection::SingleClaimOnly,
            }
        )
        .is_err()
    );
    let extras = [ScalarPropertyRequirement::try_new(
        DYNAMIC_VISCOSITY_PROPERTY,
        DYNAMIC_VISCOSITY_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )
    .unwrap()];
    assert!(
        resolve_phase_transport_state(
            &law,
            0.0,
            PhaseTransportSourceQuery {
                source: ThermalMaterialPoint {
                    material: &card,
                    context: &point,
                    temperature_axis: "T",
                    liquid_fraction_axis: Some("liquid")
                },
                viscosity_basis: LiquidViscosityBasis::Dynamic,
                extra_requirements: &extras,
                selection: MaterialPropertySelection::SingleClaimOnly,
            }
        )
        .is_err()
    );
    assert_eq!(Temperature::DIMS, Dims([0, 0, 0, 1, 0, 0]));
}
