//! Instantaneous original-source responses used by coupled process owners.
use super::*;

#[test]
fn loading_targets_keep_original_owner_and_actual_semantic_coordinates() {
    let model = equilibrium(1.0);
    for (concentration, loaded, vacant) in [
        (0.0, 0.0, 2.0),
        (1.0, 1.0, 1.0),
        (2.0, 4.0 / 3.0, 2.0 / 3.0),
    ] {
        let result = model
            .loading_at(
                AmountConcentration::new(concentration),
                Temperature::new(298.15),
                &conditions(),
            )
            .unwrap();
        assert_eq!(result.model(), &model);
        assert_eq!(result.model().source(), model.source());
        close(result.loading.value(), loaded, 3e-16);
        close(result.vacant_loading.value(), vacant, 3e-16);
        close(
            result.loading.value() + result.vacant_loading.value(),
            2.0,
            3e-16,
        );
        assert_eq!(result.concentration.value(), concentration);
        assert_eq!(result.temperature.value(), 298.15);
        assert_eq!(result.conditions.axes()[CONCENTRATION_AXIS], concentration);
        assert_eq!(result.conditions.axes()[TEMPERATURE_AXIS], 298.15);
        assert_eq!(
            result.conditions.axis_quantities()[TEMPERATURE_AXIS],
            temperature_spec()
        );
        assert_eq!(
            result.conditions.axis_quantities()[CONCENTRATION_AXIS],
            concentration_spec()
        );
    }
}

#[test]
fn local_rates_cover_uptake_desorption_stationary_and_disabled_transfer() {
    for (rate, loading, expected) in [
        (2.0, 0.25, 1.5),
        (2.0, 1.5, -1.0),
        (2.0, 1.0, 0.0),
        (0.0, 0.25, 0.0),
        (0.0, 1.5, 0.0),
    ] {
        let model = LocalLangmuirTransfer::from_card(equilibrium(1.0), rate_source(rate)).unwrap();
        let result = model
            .rate_at(
                AmountConcentration::new(1.0),
                AdsorbentLoading::new(loading),
                Volume::new(2.0),
                Temperature::new(298.15),
                &conditions(),
            )
            .unwrap();
        assert_eq!(result.model(), &model);
        assert_eq!(result.equilibrium.model(), model.equilibrium());
        assert_eq!(result.loading.value(), loading);
        assert_eq!(result.loading_rate.value(), expected);
        assert_eq!(result.conditions, result.equilibrium.conditions);
        assert_eq!(result.conditions.axes()[LOADING_AXIS], loading);
        assert_eq!(result.conditions.axes()[VOLUME_AXIS], 2.0);
        assert_eq!(
            result.conditions.axis_quantities()[LOADING_AXIS],
            QuantitySpec::dimensional(AdsorbentLoading::DIMS)
        );
        assert_eq!(
            result.conditions.axis_quantities()[VOLUME_AXIS],
            QuantitySpec::dimensional(Volume::DIMS)
        );
    }
    assert_eq!(AdsorbentLoadingRate::DIMS, Dims([0, -1, -1, 0, 0, 1]));
}

#[test]
fn trace_vacancy_and_signed_desorption_survive_saturation_but_range_loss_refuses() {
    let wide = |mut card: ConstitutiveModelCard| {
        card.validity =
            card.validity
                .with_quantity(CONCENTRATION_AXIS, concentration_spec(), 0.0, 1e100);
        card
    };
    let isotherm =
        AqueousLangmuir::from_card(ElementId::new("Au").unwrap(), wide(equilibrium_source(1.0)))
            .unwrap();
    let model = LocalLangmuirTransfer::from_card(isotherm, wide(rate_source(1.0))).unwrap();
    let result = model
        .rate_at(
            AmountConcentration::new(1e100),
            AdsorbentLoading::new(2.0),
            Volume::new(2.0),
            Temperature::new(298.15),
            &conditions(),
        )
        .unwrap();
    assert_eq!(result.equilibrium.loading.value(), 2.0);
    close(
        result.equilibrium.vacant_loading.value() / 2e-100,
        1.0,
        3e-16,
    );
    close(result.loading_rate.value() / -2e-100, 1.0, 3e-16);
    assert!(result.loading_rate.value() < 0.0);

    for (rate, concentration) in [(f64::MAX, 100.0), (1e-300, 1e-100)] {
        let model = LocalLangmuirTransfer::from_card(equilibrium(1.0), rate_source(rate)).unwrap();
        assert!(matches!(
            model.rate_at(
                AmountConcentration::new(concentration),
                AdsorbentLoading::new(0.0),
                Volume::new(2.0),
                Temperature::new(298.15),
                &conditions(),
            ),
            Err(AdsorptionFlowError::Source(
                AqueousAdsorptionError::Arithmetic(_)
            ))
        ));
    }
}

#[test]
fn instantaneous_sources_use_actual_axes_and_refuse_unsupported_or_forged_conditions() {
    let model = LocalLangmuirTransfer::from_card(equilibrium(1.0), rate_source(0.0)).unwrap();
    let evaluate = |c, q, v, t, extra: &QueryPoint| {
        model.rate_at(
            AmountConcentration::new(c),
            AdsorbentLoading::new(q),
            Volume::new(v),
            Temperature::new(t),
            extra,
        )
    };
    for (axis, quantity, _) in owned_axes() {
        let forged = conditions().with_quantity(axis, quantity, 1.0).unwrap();
        assert!(evaluate(1.0, 0.25, 2.0, 298.15, &forged).is_err());
    }
    for (axis, quantity) in [
        (CONCENTRATION_AXIS, concentration_spec()),
        (TEMPERATURE_AXIS, temperature_spec()),
    ] {
        let forged = conditions().with_quantity(axis, quantity, 1.0).unwrap();
        assert!(
            model
                .equilibrium()
                .loading_at(
                    AmountConcentration::new(1.0),
                    Temperature::new(298.15),
                    &forged
                )
                .is_err()
        );
    }
    let wrong_quantity = QueryPoint::new()
        .with_quantity("pH", concentration_spec(), 11.0)
        .unwrap();
    for extra in [QueryPoint::new(), wrong_quantity] {
        assert!(evaluate(1.0, 0.25, 2.0, 298.15, &extra).is_err());
    }
    for (c, q, v, t) in [
        (-1.0, 0.25, 2.0, 298.15),
        (f64::NAN, 0.25, 2.0, 298.15),
        (1.0, -1.0, 2.0, 298.15),
        (1.0, 2.1, 2.0, 298.15),
        (1.0, 0.25, 0.0, 298.15),
        (1.0, 0.25, 2.0, 0.0),
        (101.0, 0.25, 2.0, 298.15),
    ] {
        assert!(evaluate(c, q, v, t, &conditions()).is_err());
    }

    // Even a zero-rate card must admit both independent source domains.
    let mut narrow_rate = rate_source(0.0);
    narrow_rate.validity =
        narrow_rate
            .validity
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.5, 11.5);
    let narrow = LocalLangmuirTransfer::from_card(equilibrium(1.0), narrow_rate).unwrap();
    let outside_rate = QueryPoint::new()
        .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.25)
        .unwrap();
    assert!(
        narrow
            .rate_at(
                AmountConcentration::new(1.0),
                AdsorbentLoading::new(0.25),
                Volume::new(2.0),
                Temperature::new(298.15),
                &outside_rate
            )
            .is_err()
    );

    // Isotherm conditions for volume/loading receive the actual typed values
    // before admission, rather than caller substitutes or a borrowed source box.
    let mut bounded = equilibrium_source(1.0);
    bounded.validity = bounded
        .validity
        .with_quantity(
            VOLUME_AXIS,
            QuantitySpec::dimensional(Volume::DIMS),
            1.0,
            3.0,
        )
        .with_quantity(
            LOADING_AXIS,
            QuantitySpec::dimensional(AdsorbentLoading::DIMS),
            0.0,
            0.5,
        );
    let bounded = AqueousLangmuir::from_card(ElementId::new("Au").unwrap(), bounded).unwrap();
    let bounded = LocalLangmuirTransfer::from_card(bounded, rate_source(0.0)).unwrap();
    for (q, v, admitted) in [(0.25, 2.0, true), (0.25, 4.0, false), (0.75, 2.0, false)] {
        assert_eq!(
            bounded
                .rate_at(
                    AmountConcentration::new(1.0),
                    AdsorbentLoading::new(q),
                    Volume::new(v),
                    Temperature::new(298.15),
                    &conditions()
                )
                .is_ok(),
            admitted
        );
    }
}

#[test]
fn local_derivative_agrees_with_original_retained_sorbent_dynamics() {
    let retained = population("carbon", 3.0, 0.75, 2.0, 1.0);
    let rate = retained
        .model()
        .rate_at(
            AmountConcentration::new(1.0),
            AdsorbentLoading::new(0.25),
            Volume::new(2.0),
            Temperature::new(298.15),
            retained.conditions(),
        )
        .unwrap();
    let expected = retained.mass().value() * rate.loading_rate.value();
    assert_eq!(expected, 4.5);
    let mut tank = tank(vec![retained], 2.0);
    let h = 1.0 / 1_048_576.0;
    let result = advance(&mut tank, 0.0, 0.0, h, fixed_step(h));
    close(
        (result.populations[0].adsorbed.value() - 0.75) / h,
        expected,
        2e-5,
    );
    close((result.dissolved.value() - 2.0) / h, -expected, 2e-5);
    balances(&result, 2e-14);
}
