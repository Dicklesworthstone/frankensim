use super::*;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_material::thermal::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn resistance_key() -> PropertyKey {
    PropertyKey::new(
        AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    )
}

fn resistance_claim(r: f64, pressure_dependent: bool) -> PropertyClaim {
    let mut value = claim(resistance_key(), [r; 2]);
    if pressure_dependent {
        // Manufactured linear R(p) = r*(2-p/1 MPa), valid only on this interval.
        value.value = PropertyValue::Curve {
            abscissa: "loaded_pressure".into(),
            abscissa_dims: Dims([-1, 1, -2, 0, 0, 0]),
            knots: vec![(0.1e6, 1.9 * r), (1.1e6, 0.9 * r)],
            dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        };
        value.validity = value.validity.with("loaded_pressure", 0.1e6, 1.1e6);
    }
    value
}

fn resolved(data: &Fixture) -> ResolvedInterfaceStatePoint {
    let card = data
        .source
        .wear_source()
        .contact()
        .normal_source()
        .elastic()
        .interface()
        .card();
    let query = point(T0)
        .with_quantity(
            "loaded_pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e6,
        )
        .unwrap();
    resolve_interface_state_point(
        card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn finite_resistance_heat_partition_matches_independent_asymmetric_quadratic() {
    with_cx(|cx, _| {
        let r = 1e-4;
        let fraction = 0.25;
        let data = fixture_with_extra_claim(cx, Some(resistance_claim(r, false)));
        let thermal = network_with_resistance(&data.caloric, T0, r);
        assert!(FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).is_err());
        let resistances = [fraction * r, (1.0 - fraction) * r];
        let route = FrictionHeatRoute::from_link_with_resistance(
            &thermal,
            LINK,
            [CONDUCTIVITY; 2],
            resistances,
        )
        .unwrap();
        let source = resolved(&data);
        let contact = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx)
            .with_sourced_contact_resistance(source.clone(), "loaded_pressure".into(), fraction, cx)
            .unwrap();
        let dt = 0.1;
        let distance = 0.01;
        let step = contact
            .step_candidate(
                &thermal,
                dt,
                &motion(distance),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        let capacity = RHO * AREA * LENGTH * CP;
        let g = resistances.map(|r| AREA / (0.5 * LENGTH / CONDUCTIVITY + r));
        let b = g.map(|g| dt * g / (capacity + dt * g));
        let measure = AREA * 1e6 * distance;
        let qa = measure * 0.0001 * b[0];
        let qb = capacity * (b[0] + b[1]) - measure * (0.01 - 0.004 * b[0]);
        let theta = 2.0 * measure * 0.4 / (qb + (qb * qb + 4.0 * qa * measure * 0.4).sqrt());
        let report = step.report();
        close(report.heat.interface_temperature_k, T0 + theta, 0.0, 3e-6);
        for side in 0..2 {
            close(
                step.thermal().temperature_k(CELLS[side]).unwrap(),
                T0 + b[side] * theta,
                0.0,
                3e-6,
            );
            close(report.heat.conductances_w_k[side], g[side], 1e-14, 0.0);
            close(
                report.heat.body_heat_j[side],
                capacity * b[side] * theta,
                2e-6,
                1e-12,
            );
        }
        assert!(report.heat.body_heat_j[0] > report.heat.body_heat_j[1]);
        let expected_work = measure * (1.0 - 0.01 * b[0] * theta) * (0.4 + 0.01 * theta);
        close(report.frictional_work_j, expected_work, 2e-6, 1e-12);
        close(
            report.heat.thermal.stored_energy_change_j,
            expected_work,
            2e-6,
            3e-11,
        );
        let sample = report.contact_resistance.as_ref().unwrap();
        assert_eq!(sample.source.card_identity(), source.card_identity());
        assert_eq!(
            sample.source.properties()[0].answer().receipt.selected,
            source.properties()[0].answer().receipt.selected
        );
        assert_eq!(
            sample.source.state_coordinate("T"),
            Some(report.heat.interface_temperature_k)
        );
        assert_eq!(step.contact().route().area_resistance_m2_k_w(), resistances);
        assert_eq!(contact.route().area_resistance_m2_k_w(), resistances);
        assert_eq!(thermal.temperature_k(CELLS[0]).unwrap(), T0);
    });
}

#[test]
fn pressure_resistance_tracks_actual_worn_endpoint_and_survives_continuation() {
    with_cx(|cx, _| {
        let r = 1e-4;
        let data = fixture_with_extra_claim(cx, Some(resistance_claim(r, true)));
        let thermal = network_with_resistance(&data.caloric, T0, r);
        let route = FrictionHeatRoute::from_link_with_resistance(
            &thermal,
            LINK,
            [CONDUCTIVITY; 2],
            [0.5 * r; 2],
        )
        .unwrap();
        let contact = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx)
            .with_sourced_contact_resistance(resolved(&data), "loaded_pressure".into(), 0.5, cx)
            .unwrap();
        let first = contact
            .step_candidate(
                &thermal,
                0.1,
                &motion(0.01),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        let pressure = first
            .contact()
            .state()
            .wear_state()
            .contact_state()
            .contact_state()
            .pressure_pa()[12];
        let sample = first.report().contact_resistance.as_ref().unwrap();
        assert_eq!(sample.mean_loaded_pressure_pa, pressure);
        assert_eq!(
            sample.source.state_coordinate("loaded_pressure"),
            Some(pressure)
        );
        let expected = r * (2.0 - pressure / 1e6);
        close(
            sample.area_resistance_m2_k_w.iter().sum(),
            expected,
            1e-14,
            1e-20,
        );
        assert!(expected > r);
        assert_eq!(
            first.contact().route().area_resistance_m2_k_w(),
            sample.area_resistance_m2_k_w
        );
        let patch = first.thermal().contact_patch(LINK).unwrap();
        close(
            patch.link.conductance_w_k,
            AREA / (LENGTH / CONDUCTIVITY + expected),
            1e-14,
            1e-20,
        );
        assert_eq!(patch.covered_area_m2, [AREA; 2]);
        // Export removes actual inventory but cannot silently restore perfect
        // contact or erase the accepted resistance/source checkpoint.
        let (geometry, mut patches) =
            wear_export::export_geometry(first.contact(), first.thermal());
        let distances = geometry.map(|g| 0.5 * (g.end_m[0] - g.start_m[0]));
        patches[0]
            .replacement
            .as_mut()
            .unwrap()
            .link
            .conductance_w_k =
            series_contact_conductance_w_k(AREA, distances, [CONDUCTIVITY; 2], expected).unwrap();
        let exported = first
            .contact()
            .export_wear_candidate(
                first.thermal(),
                geometry,
                &patches,
                [CONDUCTIVITY; 2],
                fs_thermal::enthalpy::removal::MaterialRemovalTolerances {
                    mass_kg: 1e-20,
                    volume_m3: 1e-23,
                },
                cx,
            )
            .unwrap();
        assert_eq!(exported.contact().state(), first.contact().state());
        assert_eq!(
            exported.contact().route().area_resistance_m2_k_w(),
            sample.area_resistance_m2_k_w
        );
        assert_eq!(exported.contact().unexported_wear_volume_m3(), [0.0; 2]);
        let next = exported
            .contact()
            .step_candidate(
                exported.thermal(),
                0.2,
                &motion(0.001),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                first.report().heat.interface_temperature_k,
                transport,
                cx,
            )
            .unwrap();
        assert!(next.contact().contact_resistance().is_some());
        assert_eq!(
            next.report()
                .contact_resistance
                .as_ref()
                .unwrap()
                .source
                .card_identity(),
            sample.source.card_identity()
        );
        close(
            next.report().heat.thermal.stored_energy_change_j,
            next.report().frictional_work_j,
            0.0,
            3e-11,
        );
    });
}

#[test]
fn zero_resistance_replays_perfect_contact_without_double_heat() {
    with_cx(|cx, _| {
        let data = fixture_with_extra_claim(cx, Some(resistance_claim(0.0, false)));
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let plain = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx);
        let sourced = plain
            .clone()
            .with_sourced_contact_resistance(resolved(&data), "loaded_pressure".into(), 0.37, cx)
            .unwrap();
        let old = plain
            .step_candidate(
                &thermal,
                0.1,
                &motion(0.01),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        let new = sourced
            .step_candidate(
                &thermal,
                0.1,
                &motion(0.01),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                T0,
                transport,
                cx,
            )
            .unwrap();
        assert_eq!(old.contact().state(), new.contact().state());
        assert_eq!(old.thermal().cells(), new.thermal().cells());
        assert_eq!(old.report().heat.body_heat_j, new.report().heat.body_heat_j);
        assert_eq!(
            old.report().heat.interface_temperature_k,
            new.report().heat.interface_temperature_k
        );
        assert!(old.report().contact_resistance.is_none());
        assert_eq!(
            new.report()
                .contact_resistance
                .as_ref()
                .unwrap()
                .area_resistance_m2_k_w,
            [0.0; 2]
        );
    });
}

#[test]
fn resistance_binding_and_source_failure_leave_native_histories_untouched() {
    with_cx(|cx, _| {
        let r = 1e-4;
        let mut narrow = resistance_claim(r, false);
        narrow.validity = narrow.validity.with("T", 280.0, T0 + 1e-4);
        let data = fixture_with_extra_claim(cx, Some(narrow));
        let thermal = network_with_resistance(&data.caloric, T0, r);
        let route = FrictionHeatRoute::from_link_with_resistance(
            &thermal,
            LINK,
            [CONDUCTIVITY; 2],
            [0.5 * r; 2],
        )
        .unwrap();
        let plain = bind(&data, &thermal, &route, loaded(&data, cx), limits(), cx);
        let source = resolved(&data);
        assert!(
            plain
                .clone()
                .with_sourced_contact_resistance(source.clone(), "loaded_pressure".into(), 0.1, cx)
                .is_err()
        );
        assert!(
            plain
                .clone()
                .with_sourced_contact_resistance(source.clone(), "T".into(), 0.5, cx)
                .is_err()
        );
        let contact = plain
            .with_sourced_contact_resistance(source, "loaded_pressure".into(), 0.5, cx)
            .unwrap();
        let before = contact.state().clone();
        let before_cells = thermal.cells().clone();
        let failed = contact.step_candidate(
            &thermal,
            0.1,
            &motion(0.01),
            PlasticSlidingHeatBoundary::HeldApproach,
            controls(),
            &[0.0; 2],
            T0,
            transport,
            cx,
        );
        assert!(matches!(
            failed,
            Err(PlasticSlidingHeatError::ResistanceSource(_))
        ));
        assert_eq!(contact.state(), &before);
        assert_eq!(thermal.cells(), &before_cells);
        assert_eq!(contact.route().area_resistance_m2_k_w(), [0.5 * r; 2]);
    });
}
