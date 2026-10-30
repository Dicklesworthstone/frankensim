//! Joint constitutive admission on the original h chart; manufactured data.
use super::*;

#[test]
fn phase_density_and_native_rate_share_original_source_without_changing_mass_storage() {
    let original = shear_card(None);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DENSITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(400.0, 1010.0), (500.0, 1000.0), (600.0, 990.0)],
                dims: Density::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = curve(&card, 1000.0);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1000.0, 1000.0, 1e-8).unwrap();
    let context = point()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 3.0)
        .unwrap();
    let query = PhaseTransportSourceQuery {
        source: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("liquid"),
        },
        viscosity_basis: LiquidViscosityBasis::Dynamic,
        extra_requirements: &[],
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let resolve = |h, rate, limit| {
        resolve_phase_transport_state_with_shear_rate_and_density_departure(
            &law,
            h,
            query.clone(),
            "rate",
            rate,
            limit,
        )
    };
    for (h, rate, rho, k, mu) in [
        (-50.0, 0.0, 1005.0, 10.0, None),
        (25.0, 2.0, 1000.0, 9.0, Some(0.03)),
        (150.0, 4.0, 995.0, 4.0, Some(0.04)),
    ] {
        let state = resolve(h, rate, 0.02).unwrap();
        assert!((state.source_density_kg_m3() - rho).abs() < 1e-12);
        assert_eq!(state.phase_state().bulk_density_kg_m3(), 1000.0);
        assert_eq!(
            state.phase_state(),
            &curve.state_at_specific_enthalpy(h).unwrap()
        );
        assert_eq!(state.material().card_identity(), card.content_hash());
        assert!((state.thermal_conductivity_w_m_k().unwrap() - k).abs() < 1e-12);
        if let Some(mu) = mu {
            let transport = state.transport_state().unwrap();
            assert_eq!(transport.density_kg_m3(), state.source_density_kg_m3());
            assert!((transport.dynamic_viscosity_pa_s() - mu).abs() < 1e-14);
        } else {
            assert!(state.mobile_transport().is_none());
            assert!(
                state
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .is_none()
            );
            assert!(
                state
                    .material()
                    .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
                    .is_none()
            );
        }
        for property in state.material().properties() {
            let receipt = &property.answer().receipt;
            assert!(receipt.query_point.contains(&("rate".into(), rate)));
            assert!(
                receipt
                    .query_point
                    .contains(&("T".into(), state.phase_state().temperature_k(),))
            );
            assert!(
                receipt
                    .query_point
                    .contains(&("liquid".into(), state.phase_state().liquid_mass_fraction(),))
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
        }
    }
    // Neither original one-coordinate entry can silently supply the missing owner.
    assert!(
        resolve_phase_transport_state_with_shear_rate(&law, 150.0, query.clone(), "rate", 2.0,)
            .is_err()
    );
    assert!(
        resolve_phase_transport_state_with_density_departure(&law, 150.0, query.clone(), 0.02,)
            .is_err()
    );
    assert!(resolve(-50.0, 0.01, 0.02).is_err());
    assert!(resolve(-50.0, 0.0, 0.001).is_err());
    for rate in [-1.0, 4.1, f64::NAN] {
        assert!(resolve(25.0, rate, 0.02).is_err());
    }
    for limit in [-0.1, 1.0, f64::NAN] {
        assert!(resolve(25.0, 0.0, limit).is_err());
    }
}

#[test]
fn density_departure_does_not_admit_shear_dependent_thermodynamic_storage() {
    for name in [DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY] {
        let card = shear_card(Some(name));
        let curve = curve(&card, 1000.0);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1000.0, 1000.0, 1e-8).unwrap();
        let context = point()
            .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
            .unwrap();
        assert!(
            resolve_phase_transport_state_with_shear_rate_and_density_departure(
                &law,
                25.0,
                PhaseTransportSourceQuery {
                    source: ThermalMaterialPoint {
                        material: &card,
                        context: &context,
                        temperature_axis: "T",
                        liquid_fraction_axis: Some("liquid"),
                    },
                    viscosity_basis: LiquidViscosityBasis::Dynamic,
                    extra_requirements: &[],
                    selection: MaterialPropertySelection::SingleClaimOnly,
                },
                "rate",
                0.0,
                0.02,
            )
            .is_err()
        );
    }
}
