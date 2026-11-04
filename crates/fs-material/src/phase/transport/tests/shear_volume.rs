//! Physical-rate source response on the same changing phase-density chart.

use super::*;

fn changing_card(basis: LiquidViscosityBasis) -> MaterialCard {
    let original = shear_card(None);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DENSITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "liquid".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 4.0), (0.5, 8.0 / 3.0), (1.0, 2.0)],
                dims: Density::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        } else if claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY
            && basis == LiquidViscosityBasis::Kinematic
        {
            claim.key = PropertyKey::new(KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS);
            claim.value = PropertyValue::Curve {
                abscissa: "rate".into(),
                abscissa_dims: Frequency::DIMS,
                knots: vec![(0.0, 2e-5), (4.0, 4e-5)],
                dims: KINEMATIC_VISCOSITY_DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn changing_curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 400.0, 0.0, 4.0),
            (0.0, 500.0, 0.0, 4.0),
            (100.0, 500.0, 1.0, 2.0),
            (200.0, 600.0, 1.0, 2.0),
        ]
        .into_iter()
        .map(|(h, t, f, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: rho,
        })
        .collect(),
    )
    .unwrap()
}

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

#[test]
fn shear_response_and_current_phase_density_share_original_pinned_sources() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let card = changing_card(basis);
        let curve = changing_curve(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1e-8).unwrap();
        let context = point()
            .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 3.0)
            .unwrap();
        let source = query(&card, &context, basis);
        assert!(
            resolve_phase_transport_state_for_volume_change(&law, 50.0, source.clone()).is_err()
        );
        assert!(
            resolve_phase_transport_state_with_shear_rate(&law, 50.0, source.clone(), "rate", 2.0)
                .is_err()
        );
        let first = resolve_phase_transport_state_with_shear_for_volume_change(
            &law,
            50.0,
            source.clone(),
            "rate",
            2.0,
        )
        .unwrap();
        let density = 8.0 / 3.0;
        assert_eq!(first.transport_density_kg_m3(), density);
        assert_eq!(first.phase_state().bulk_density_kg_m3(), density);
        assert_eq!(first.phase_state().liquid_mass_fraction(), 0.5);
        assert_eq!(first.mobility().volumes().liquid, 2.0 / 3.0);
        assert_eq!(first.thermal_conductivity_w_m_k().unwrap(), 7.0);
        let expected_mu = match basis {
            LiquidViscosityBasis::Dynamic => 0.03,
            LiquidViscosityBasis::Kinematic => density * 3e-5,
        };
        assert!(
            (first.transport_state().unwrap().dynamic_viscosity_pa_s() - expected_mu).abs() < 1e-16
        );
        let pins = first
            .material()
            .properties()
            .iter()
            .map(|property| {
                (
                    property.requirement().name().to_owned(),
                    property.answer().receipt.selected,
                )
            })
            .collect();
        let next = resolve_phase_transport_state_with_shear_for_volume_change(
            &law,
            150.0,
            PhaseTransportSourceQuery {
                selection: MaterialPropertySelection::PinnedByProperty(pins),
                ..source
            },
            "rate",
            4.0,
        )
        .unwrap();
        assert_eq!(next.transport_density_kg_m3(), 2.0);
        assert_eq!(next.phase_state().bulk_density_kg_m3(), 2.0);
        assert_eq!(next.thermal_conductivity_w_m_k().unwrap(), 4.0);
        assert_eq!(
            next.transport_state().unwrap().dynamic_viscosity_pa_s(),
            match basis {
                LiquidViscosityBasis::Dynamic => 0.04,
                LiquidViscosityBasis::Kinematic => 8e-5,
            }
        );
        for (state, t, f, rate) in [(&first, 500.0, 0.5, 2.0), (&next, 550.0, 1.0, 4.0)] {
            assert_eq!(state.material().card_identity(), card.content_hash());
            assert_eq!(state.phase_state().phase_curve_identity(), curve.identity());
            assert_eq!(state.mobility().law_identity(), law.identity());
            for property in state.material().properties() {
                let receipt = &property.answer().receipt;
                for (axis, value) in [
                    ("T", t),
                    ("liquid", f),
                    ("rate", rate),
                    ("pressure", 100000.0),
                ] {
                    assert!(receipt.query_point.contains(&(axis.into(), value)));
                }
                assert_eq!(
                    receipt.selected,
                    first
                        .material()
                        .property(property.requirement().name())
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
            }
        }
        // Authored rate replaces the stale query coordinate without mutating it.
        assert_eq!(context.axes().get("rate"), Some(&3.0));
        assert_eq!(context.axes().get("liquid"), Some(&0.9));
    }
}

#[test]
fn exact_solid_keeps_no_fluid_tuple_and_changing_density_does_not_relax_rate_admission() {
    let full = changing_card(LiquidViscosityBasis::Dynamic);
    let mut claims = ClaimSet::new();
    for (_, claim) in full.claims().claims_ordered() {
        if ![DYNAMIC_VISCOSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY]
            .contains(&claim.key.name())
        {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    let card = MaterialCard::assemble(full.id().clone(), claims, vec![]).unwrap();
    let curve = changing_curve(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1e-8).unwrap();
    let context = point()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 3.0)
        .unwrap();
    let source = query(&card, &context, LiquidViscosityBasis::Dynamic);
    let solid = resolve_phase_transport_state_with_shear_for_volume_change(
        &law,
        -50.0,
        source.clone(),
        "rate",
        0.0,
    )
    .unwrap();
    assert_eq!(solid.transport_density_kg_m3(), 4.0);
    assert_eq!(solid.phase_state().phase(), SolidLiquidPhase::Solid);
    assert!(solid.mobile_transport().is_none());
    assert!(solid.transport_state().is_err());
    assert_eq!(solid.material().properties().len(), 2);
    for (h, rate) in [(-50.0, 0.01), (50.0, 2.0)] {
        assert!(
            resolve_phase_transport_state_with_shear_for_volume_change(
                &law,
                h,
                source.clone(),
                "rate",
                rate
            )
            .is_err()
        );
    }
    let curve = changing_curve(&full);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1e-8).unwrap();
    let source = query(&full, &context, LiquidViscosityBasis::Dynamic);
    for rate in [-1.0, 4.1, f64::NAN, f64::INFINITY] {
        assert!(
            resolve_phase_transport_state_with_shear_for_volume_change(
                &law,
                50.0,
                source.clone(),
                "rate",
                rate
            )
            .is_err()
        );
    }
    for axis in ["absent", "T", "liquid", "pressure"] {
        assert!(
            resolve_phase_transport_state_with_shear_for_volume_change(
                &law,
                50.0,
                source.clone(),
                axis,
                0.0
            )
            .is_err()
        );
    }
    for context in [
        point(),
        point().with("rate", 2.0).unwrap(),
        point()
            .with_quantity("rate", QuantitySpec::dimensional(Temperature::DIMS), 2.0)
            .unwrap(),
    ] {
        assert!(
            resolve_phase_transport_state_with_shear_for_volume_change(
                &law,
                50.0,
                query(&full, &context, LiquidViscosityBasis::Dynamic),
                "rate",
                2.0,
            )
            .is_err()
        );
    }
    // A correctly typed rate cannot repair unsupported or inconsistent h/rho.
    for h in [-101.0, 25.0, 201.0, f64::NAN] {
        assert!(
            resolve_phase_transport_state_with_shear_for_volume_change(
                &law,
                h,
                source.clone(),
                "rate",
                2.0
            )
            .is_err()
        );
    }
    assert!(
        resolve_phase_transport_state_with_shear_for_volume_change(
            &law,
            50.0,
            query(&card, &context, LiquidViscosityBasis::Dynamic),
            "rate",
            2.0,
        )
        .is_err()
    );
}

#[test]
fn original_fixed_shear_states_match_and_shear_dependent_storage_still_refuses() {
    let context = point()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 3.0)
        .unwrap();
    let card = shear_card(None);
    let chart = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&chart, 1000.0, 1000.0, 1e-8).unwrap();
    let source = query(&card, &context, LiquidViscosityBasis::Dynamic);
    for (h, rate) in [(-50.0, 0.0), (25.0, 2.0), (150.0, 4.0)] {
        assert_eq!(
            resolve_phase_transport_state_with_shear_rate(&law, h, source.clone(), "rate", rate)
                .unwrap(),
            resolve_phase_transport_state_with_shear_for_volume_change(
                &law,
                h,
                source.clone(),
                "rate",
                rate
            )
            .unwrap(),
        );
    }
    for property in [DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY] {
        let card = shear_card(Some(property));
        let chart = curve(&card, 1000.0);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&chart, 1000.0, 1000.0, 1e-8).unwrap();
        let error = resolve_phase_transport_state_with_shear_for_volume_change(
            &law,
            25.0,
            query(&card, &context, LiquidViscosityBasis::Dynamic),
            "rate",
            0.0,
        )
        .unwrap_err();
        assert!(matches!(error, MaterialError::State { what }
            if what.contains("shear-dependent density or heat capacity")));
    }
}
