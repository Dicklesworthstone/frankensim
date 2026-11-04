//! Original source density and phase transport for a physical-volume consumer.

use super::*;
use crate::phase::HeatCapacitySourceQuery;
use crate::phase::melting::{MeltingSourceQuery, ResolvedMeltingCurve};
use fs_matdb::SelectionPolicy;

fn query<'a>(
    card: &'a MaterialCard,
    context: &'a QueryPoint,
    basis: LiquidViscosityBasis,
) -> PhaseTransportSourceQuery<'a> {
    PhaseTransportSourceQuery {
        source: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("liquid"),
        },
        viscosity_basis: basis,
        extra_requirements: &[],
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

fn unequal_card(liquid: bool, basis: LiquidViscosityBasis, middle_density: f64) -> MaterialCard {
    let original = card(4.0, liquid, basis);
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut source = source.clone();
        if source.key.name() == DENSITY_PROPERTY {
            source.value = PropertyValue::Curve {
                abscissa: "liquid".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 4.0), (0.5, middle_density), (1.0, 2.0)],
                dims: Density::DIMS,
            };
            source.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(source).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn unequal_curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 400.0, 0.0, 4.0),
            (0.0, 500.0, 0.0, 4.0),
            (100.0, 500.0, 1.0, 2.0),
            (200.0, 600.0, 1.0, 2.0),
        ]
        .into_iter()
        .map(|(h, t, fraction, density)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: density,
        })
        .collect(),
    )
    .unwrap()
}

#[test]
fn unequal_phase_sources_retain_actual_density_viscosity_and_volume_fraction() {
    let card = unequal_card(true, LiquidViscosityBasis::Kinematic, 8.0 / 3.0);
    let curve = unequal_curve(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1e-8).unwrap();
    assert!(law.fixed_density_kg_m3().is_none());
    let context = point();
    let extras = [ElectricalPropertyBasis::Resistivity.requirement().unwrap()];
    let mut source = query(&card, &context, LiquidViscosityBasis::Kinematic);
    source.extra_requirements = &extras;
    for (h, density) in [
        (-50.0, 4.0_f64),
        (0.0, 4.0),
        (50.0, 8.0 / 3.0),
        (150.0, 2.0),
    ] {
        let state =
            resolve_phase_transport_state_for_volume_change(&law, h, source.clone()).unwrap();
        let phase = curve.state_at_specific_enthalpy(h).unwrap();
        assert_eq!(state.phase_state(), &phase);
        assert_eq!(state.source_density_kg_m3().to_bits(), density.to_bits());
        assert_eq!(
            state.transport_density_kg_m3().to_bits(),
            phase.bulk_density_kg_m3().to_bits()
        );
        assert_eq!(state.material().card_identity(), card.content_hash());
        if h > 0.0 {
            let transport = state.transport_state().unwrap();
            assert_eq!(transport.density_kg_m3(), density);
            assert_eq!(transport.dynamic_viscosity_pa_s(), density * 2e-5);
            assert!((transport.kinematic_viscosity_m2_s() - 2e-5).abs() < 1e-20);
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
        } else {
            assert!(state.mobile_transport().is_none());
            assert!(state.transport_state().is_err());
            assert!(
                state
                    .material()
                    .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                    .is_none()
            );
        }
        for property in state.material().properties() {
            let receipt = &property.answer().receipt;
            assert!(
                receipt
                    .query_point
                    .contains(&("T".into(), phase.temperature_k()))
            );
            assert!(
                receipt
                    .query_point
                    .contains(&("liquid".into(), phase.liquid_mass_fraction()))
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
        }
        assert!(resolve_phase_transport_state(&law, h, source.clone()).is_err());
        assert!(
            resolve_phase_transport_state_with_density_departure(&law, h, source.clone(), 0.5)
                .is_err()
        );
    }
    let mixed = resolve_phase_transport_state_for_volume_change(&law, 50.0, source).unwrap();
    assert_eq!(mixed.phase_state().liquid_mass_fraction(), 0.5);
    assert_eq!(mixed.mobility().volumes().liquid, 2.0 / 3.0);
    assert_eq!(mixed.thermal_conductivity_w_m_k().unwrap(), 6.0);
    // One kg occupies 1/4 m3 before melting, 3/8 m3 at half liquid MASS,
    // and 1/2 m3 when liquid. Source transport cannot freeze any of these.
    assert_eq!(1.0 / mixed.transport_density_kg_m3(), 3.0 / 8.0);
}

#[test]
fn volume_change_refuses_density_disagreement_without_repair_or_solid_fallback() {
    let card = unequal_card(false, LiquidViscosityBasis::Dynamic, 8.0 / 3.0);
    let curve = unequal_curve(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1e-8).unwrap();
    let context = point();
    let source = query(&card, &context, LiquidViscosityBasis::Dynamic);
    let solid =
        resolve_phase_transport_state_for_volume_change(&law, -50.0, source.clone()).unwrap();
    assert_eq!(solid.phase_state().phase(), SolidLiquidPhase::Solid);
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
    assert!(resolve_phase_transport_state_for_volume_change(&law, 50.0, source).is_err());

    for middle_density in [8.0 / 3.0, f64::from_bits((8.0_f64 / 3.0).to_bits() + 1)] {
        let card = unequal_card(true, LiquidViscosityBasis::Dynamic, middle_density);
        let curve = unequal_curve(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1e-8).unwrap();
        let source = query(&card, &context, LiquidViscosityBasis::Dynamic);
        // Source linear-rho interpolation is not the chart's additive volume
        // between the tabulated matching states. Even a one-ULP mismatch at
        // the middle knot refuses: no invented agreement tolerance is used.
        assert!(
            resolve_phase_transport_state_for_volume_change(&law, 25.0, source.clone()).is_err()
        );
        assert_eq!(
            resolve_phase_transport_state_for_volume_change(&law, 50.0, source.clone()).is_ok(),
            middle_density.to_bits() == (8.0_f64 / 3.0).to_bits()
        );
        for h in [-101.0, 201.0, f64::NAN] {
            assert!(
                resolve_phase_transport_state_for_volume_change(&law, h, source.clone()).is_err()
            );
        }
        let other = unequal_card(true, LiquidViscosityBasis::Kinematic, middle_density);
        assert!(
            resolve_phase_transport_state_for_volume_change(
                &law,
                50.0,
                query(&other, &context, LiquidViscosityBasis::Kinematic)
            )
            .is_err()
        );
    }
}

#[test]
fn equal_density_state_is_unchanged_and_boussinesq_departure_stays_explicit() {
    let context = point();
    let card = card(1000.0, true, LiquidViscosityBasis::Dynamic);
    let curve = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1000.0, 1000.0, 1e-8).unwrap();
    let source = query(&card, &context, LiquidViscosityBasis::Dynamic);
    for h in [-50.0, 50.0, 150.0] {
        assert_eq!(
            resolve_phase_transport_state(&law, h, source.clone()).unwrap(),
            resolve_phase_transport_state_for_volume_change(&law, h, source.clone()).unwrap()
        );
    }
    let departure_card = super::card(990.0, true, LiquidViscosityBasis::Dynamic);
    let departure_curve = super::curve(&departure_card, 1000.0);
    let departure_law =
        MushyPhaseLaw::try_new_with_immobile_solid(&departure_curve, 1000.0, 1000.0, 1e-8).unwrap();
    let source = query(&departure_card, &context, LiquidViscosityBasis::Dynamic);
    assert!(
        resolve_phase_transport_state_with_density_departure(
            &departure_law,
            50.0,
            source.clone(),
            0.02
        )
        .is_ok()
    );
    assert!(resolve_phase_transport_state_for_volume_change(&departure_law, 50.0, source).is_err());
}

#[test]
fn sensible_expansion_uses_original_melting_density_at_actual_temperature() {
    let original = card(4.0, true, LiquidViscosityBasis::Dynamic);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DENSITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(400.0, 4.0), (600.0, 3.0)],
                dims: Density::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        } else if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            claim.value = PropertyValue::Scalar {
                value: 2.0,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            };
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    400.0,
                    600.0,
                )
                .with_quantity("liquid", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
        }
        claims.insert_claim(claim).unwrap();
    }
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let latent_dims = Dims([2, 0, -2, 0, 0, 0]);
    for (name, quantity, value) in [
        ("melting-point", absolute, 500.0),
        (
            "latent-heat-fusion",
            QuantitySpec::dimensional(latent_dims),
            100.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured expanding phase source".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let context = point();
    let branch_point = |t, fraction| {
        context
            .clone()
            .with("T", t)
            .unwrap()
            .with("liquid", fraction)
            .unwrap()
    };
    let points = [
        branch_point(400.0, 0.0),
        branch_point(500.0, 0.0),
        branch_point(500.0, 1.0),
        branch_point(600.0, 1.0),
    ];
    let transition = QueryPoint::new()
        .with_quantity("T", absolute, 500.0)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100000.0,
        )
        .unwrap();
    let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    let tm = PropertyKey::with_quantity("melting-point", absolute);
    let latent = PropertyKey::new("latent-heat-fusion", latent_dims);
    let branch = |i| HeatCapacitySourceQuery {
        heat_capacity: &cp,
        density: &rho,
        temperature_axis: "T",
        lower: &points[i],
        upper: &points[i + 1],
        selection: SelectionPolicy::SingleClaimOnly,
    };
    let melting = ResolvedMeltingCurve::try_new(
        &card,
        MeltingSourceQuery {
            solid: branch(0),
            liquid: branch(2),
            melting_point: &tm,
            latent_heat_fusion: &latent,
            transition: &transition,
            transition_selection: SelectionPolicy::SingleClaimOnly,
            phase_selector: Some("liquid"),
        },
        0.0,
        0.1,
    )
    .unwrap();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    assert!(law.fixed_density_kg_m3().is_none());
    for (h, t, density) in [
        (100.0, 450.0, 3.75),
        (250.0, 500.0, 3.5),
        (400.0, 550.0, 3.25),
    ] {
        let source = query(&card, &context, LiquidViscosityBasis::Dynamic);
        let state =
            resolve_phase_transport_state_for_volume_change(&law, h, source.clone()).unwrap();
        assert_eq!(state.phase_state().temperature_k(), t);
        assert_eq!(state.phase_state().bulk_density_kg_m3(), density);
        assert_eq!(state.transport_density_kg_m3(), density);
        if h > 200.0 {
            assert_eq!(
                state.transport_state().unwrap().dynamic_viscosity_pa_s(),
                0.02
            );
        }
        assert!(resolve_phase_transport_state(&law, h, source).is_err());
    }
}
