//! Original R(T_interface, final loaded pressure) through depth wear and heat.

use super::*;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_material::thermal::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

pub(super) const PRESSURE: &str = "loaded_pressure";
pub(super) const R0: f64 = 1e-4;

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

pub(super) fn resistance_claim(r: f64, variable: bool) -> PropertyClaim {
    let mut property = claim(
        PropertyKey::new(
            AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        ),
        [r; 2],
    );
    property.validity =
        property
            .validity
            .with_quantity(PRESSURE, pressure_quantity(), 0.1e6, 1.1e6);
    if variable {
        // Affine manufactured source; both independent coordinates are active.
        property.value = PropertyValue::Surface2D {
            axes: [PRESSURE.into(), "T".into()],
            axis_dims: [Dims([-1, 1, -2, 0, 0, 0]), Temperature::DIMS],
            samples: [
                [0.1e6, 280.0],
                [1.1e6, 280.0],
                [1.1e6, 320.0],
                [0.1e6, 320.0],
            ]
            .into_iter()
            .map(|[p, t]| ([p, t], r * (2.0 - p / 1e6 + 0.03 * (t - T0))))
            .collect(),
            triangles: vec![[0, 1, 2], [0, 2, 3]],
            dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
        };
    }
    property
}

pub(super) fn resolved(data: &DepthFixture) -> ResolvedInterfaceStatePoint {
    let card = data
        .source
        .wear_state()
        .normal_source()
        .elastic()
        .interface()
        .card();
    let query = point(T0)
        .with_quantity(PRESSURE, pressure_quantity(), 1e6)
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

pub(super) fn bound<'a>(
    data: &'a DepthFixture,
    thermal: &DepositionSimulation<'_>,
    r: f64,
    fraction: f64,
    cx: &Cx<'_>,
) -> DepthSlidingHeatContact<'a> {
    let route = FrictionHeatRoute::from_link_with_resistance(
        thermal,
        LINK,
        [CONDUCTIVITY; 2],
        [fraction * r, (1.0 - fraction) * r],
    )
    .unwrap();
    bind_depth(data, thermal, &route, depth_loaded(data, cx), limits(), cx)
        .with_sourced_contact_resistance(resolved(data), PRESSURE.into(), fraction, cx)
        .unwrap()
}

#[test]
fn g1_depth_resistance_heat_matches_independent_asymmetric_quadratic() {
    with_cx(|cx, _| {
        let data = depth_fixture_with_extra_claim(320.0, cx, Some(resistance_claim(R0, false)));
        let thermal = network_with_resistance(&data.caloric, T0, R0);
        let fraction = 0.25;
        let contact = bound(&data, &thermal, R0, fraction, cx);
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
        // W = A H(T_A, z_old=0) s mu(T_interface). Native wear changes the
        // endpoint cap; default plastic heating is zero in this independent rung.
        let capacity = RHO * AREA * LENGTH * CP;
        let resistances = [fraction * R0, (1.0 - fraction) * R0];
        let g = resistances.map(|r| AREA / (0.5 * LENGTH / CONDUCTIVITY + r));
        let b = g.map(|value| dt * value / (capacity + dt * value));
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
        let work = measure * (1.0 - 0.01 * b[0] * theta) * (0.4 + 0.01 * theta);
        close(report.frictional_work_j, work, 2e-6, 1e-12);
        close(
            report.heat.thermal.stored_energy_change_j,
            work,
            2e-6,
            3e-11,
        );
        assert_eq!(report.normal_plastic_heat_j, [0.0; 2]);
        assert_eq!(step.contact().route().area_resistance_m2_k_w(), resistances);
        assert_eq!(thermal.temperature_k(CELLS[0]).unwrap(), T0);
    });
}

#[test]
fn g3_depth_resistance_samples_final_pressure_and_temperature_through_export() {
    with_cx(|cx, _| {
        let data = depth_fixture_with_extra_claim(320.0, cx, Some(resistance_claim(R0, true)));
        let thermal = network_with_resistance(&data.caloric, T0, R0);
        let contact = bound(&data, &thermal, R0, 0.5, cx);
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
        let report = first.report();
        let pressure = first
            .contact()
            .source()
            .wear_state()
            .contact_state()
            .pressure_pa()[12];
        let beginning_pressure = report
            .mechanics
            .temperature
            .candidate
            .wear_state()
            .contact_state()
            .pressure_pa()[12];
        assert!(pressure < beginning_pressure);
        let sample = report.contact_resistance.as_ref().unwrap();
        let temperature = report.heat.interface_temperature_k;
        assert_eq!(sample.mean_loaded_pressure_pa, pressure);
        assert_eq!(sample.source.state_coordinate(PRESSURE), Some(pressure));
        assert_eq!(sample.source.state_coordinate("T"), Some(temperature));
        let original = contact.contact_resistance().unwrap().original_source();
        assert_eq!(sample.source.card_identity(), original.card_identity());
        assert_eq!(
            sample.source.properties()[0].answer().receipt.selected,
            original.properties()[0].answer().receipt.selected
        );
        let expected = R0 * (2.0 - pressure / 1e6 + 0.03 * (temperature - T0));
        assert!(expected > R0);
        close(
            sample.area_resistance_m2_k_w.iter().sum(),
            expected,
            1e-14,
            1e-20,
        );
        assert_eq!(
            first.contact().route().area_resistance_m2_k_w(),
            sample.area_resistance_m2_k_w
        );
        close(
            first.thermal().links()[&LINK].conductance_w_k,
            AREA / (LENGTH / CONDUCTIVITY + expected),
            1e-14,
            1e-20,
        );
        for side in 0..2 {
            close(
                report.heat.conductances_w_k[side],
                AREA / (0.5 * LENGTH / CONDUCTIVITY + sample.area_resistance_m2_k_w[side]),
                1e-14,
                1e-20,
            );
        }
        close(
            report.heat.thermal.stored_energy_change_j,
            report.frictional_work_j,
            0.0,
            3e-11,
        );

        let (geometry, mut patches) =
            removal::export_geometry(first.contact().unexported_wear_volume_m3(), first.thermal());
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
        assert_eq!(exported.contact().source(), first.contact().source());
        assert_eq!(exported.contact().unexported_wear_volume_m3(), [0.0; 2]);
        assert_eq!(
            exported.contact().route().area_resistance_m2_k_w(),
            sample.area_resistance_m2_k_w
        );
        assert_eq!(
            exported
                .contact()
                .contact_resistance()
                .unwrap()
                .original_source(),
            original
        );
        let next = exported
            .contact()
            .step_candidate(
                exported.thermal(),
                0.2,
                &motion(-0.001),
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                &[0.0; 2],
                temperature,
                transport,
                cx,
            )
            .unwrap();
        let next_sample = next.report().contact_resistance.as_ref().unwrap();
        assert_eq!(next_sample.source.card_identity(), original.card_identity());
        assert_eq!(
            next_sample.mean_loaded_pressure_pa,
            next.contact()
                .source()
                .wear_state()
                .contact_state()
                .pressure_pa()[12]
        );
        assert_eq!(
            next.contact()
                .contact_resistance()
                .unwrap()
                .original_source(),
            original
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
fn g3_zero_resistance_replays_native_depth_history_and_heat() {
    with_cx(|cx, _| {
        let data = depth_fixture_with_extra_claim(320.0, cx, Some(resistance_claim(0.0, false)));
        let thermal = network(&data.caloric, T0);
        let route = FrictionHeatRoute::from_link(&thermal, LINK, [CONDUCTIVITY; 2]).unwrap();
        let plain = bind_depth(
            &data,
            &thermal,
            &route,
            depth_loaded(&data, cx),
            limits(),
            cx,
        );
        let sourced = plain
            .clone()
            .with_sourced_contact_resistance(resolved(&data), PRESSURE.into(), 0.37, cx)
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
        assert_eq!(old.contact().source(), new.contact().source());
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
fn g0_resistance_refuses_foreign_sources_and_unsupported_worn_endpoints_atomically() {
    with_cx(|cx, _| {
        for pressure_limited in [false, true] {
            let mut narrow = resistance_claim(R0, true);
            narrow.validity = if pressure_limited {
                narrow
                    .validity
                    .with_quantity(PRESSURE, pressure_quantity(), 0.99e6, 1.1e6)
            } else {
                narrow.validity.with("T", 280.0, T0 + 1e-4)
            };
            let data = depth_fixture_with_extra_claim(320.0, cx, Some(narrow));
            let thermal = network_with_resistance(&data.caloric, T0, R0);
            let route = FrictionHeatRoute::from_link_with_resistance(
                &thermal,
                LINK,
                [CONDUCTIVITY; 2],
                [0.5 * R0; 2],
            )
            .unwrap();
            let plain = bind_depth(
                &data,
                &thermal,
                &route,
                depth_loaded(&data, cx),
                limits(),
                cx,
            );
            let source = resolved(&data);
            assert!(
                plain
                    .clone()
                    .with_sourced_contact_resistance(source.clone(), PRESSURE.into(), 0.1, cx,)
                    .is_err()
            );
            assert!(
                plain
                    .clone()
                    .with_sourced_contact_resistance(source.clone(), "T".into(), 0.5, cx,)
                    .is_err()
            );
            let foreign =
                depth_fixture_with_extra_claim(320.0, cx, Some(resistance_claim(2.0 * R0, true)));
            assert!(matches!(
                plain.clone().with_sourced_contact_resistance(
                    resolved(&foreign),
                    PRESSURE.into(),
                    0.5,
                    cx,
                ),
                Err(PlasticSlidingHeatError::Invalid(_))
            ));
            let contact = plain
                .with_sourced_contact_resistance(source, PRESSURE.into(), 0.5, cx)
                .unwrap();
            let before = contact.source().clone();
            let before_cells = thermal.cells().clone();
            // The pressure-limited source accepts the starting 1 MPa but must
            // refuse pressure AFTER recession and exposed-depth cap refresh.
            assert!(matches!(
                contact.step_candidate(
                    &thermal,
                    0.1,
                    &motion(0.01),
                    PlasticSlidingHeatBoundary::HeldApproach,
                    controls(),
                    &[0.0; 2],
                    T0,
                    transport,
                    cx,
                ),
                Err(PlasticSlidingHeatError::ResistanceSource(_))
            ));
            assert_eq!(contact.source(), &before);
            assert_eq!(thermal.cells(), &before_cells);
            assert_eq!(contact.route().area_resistance_m2_k_w(), [0.5 * R0; 2]);
            assert_eq!(contact.unexported_wear_volume_m3(), [0.0; 2]);
        }
    });
}
