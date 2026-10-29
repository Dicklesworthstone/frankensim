use super::*;
use crate::multicomponent::CaloricDomain;
use crate::multicomponent::tests::{bank, source};
use fs_qty::{Area, Length};
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn mechanism(heat: f64, rate: f64, variable_cp: bool, reference: f64) -> IdealGasKineticsV1 {
    let gas = bank(vec![
        source("A", R, 3.5, if variable_cp { 0.05 } else { 0.0 }, reference),
        source(
            "B",
            R,
            if variable_cp { 4.0 } else { 3.5 },
            if variable_cp { 0.15 } else { 0.0 },
            reference - heat,
        ),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(rate),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    IdealGasKineticsV1::new(
        gas,
        "zone-kinetics-manufactured-isomers",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap()
}
fn controls() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-10,
        relative_tolerance: 1e-6,
        absolute_amount_tolerance_mol: 1e-14,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-4,
        maximum_step_s: 1.0,
        maximum_steps: 30000,
        ..IsochoricKineticsSolve::default()
    }
}
fn zone_controls() -> ZoneSolve {
    ZoneSolve {
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
        ..ZoneSolve::default()
    }
}
fn room(m: &IdealGasKineticsV1, fraction: f64, t: [f64; 2], y: [[f64; 2]; 2]) -> TwoZone {
    let gas = MixtureGas::new(
        m.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    TwoZone::new(
        gas,
        Area::new(0.5),
        Length::new(2.0),
        Length::new(2.0 * fraction),
        Pressure::new(4.0),
        t.map(Temperature::new),
        y.map(Vec::from),
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance {tolerance:e}"
    );
}
fn energy_by_original_species(room: &TwoZone) -> f64 {
    [Layer::Lower, Layer::Upper]
        .into_iter()
        .map(|l| {
            let state = room.layer(l);
            state
                .species_masses()
                .iter()
                .zip(room.gas().source().models())
                .map(|(mass, model)| {
                    mass / model.molar_mass().value()
                        * model
                            .evaluate_properties(state.temperature())
                            .unwrap()
                            .u()
                            .value()
                })
                .sum::<f64>()
        })
        .sum()
}
fn check_original_eos(room: &TwoZone) {
    let mut volume = 0.0;
    for layer in [Layer::Lower, Layer::Upper] {
        let state = room.layer(layer);
        volume += state.volume().value();
        assert!(state.species_masses().iter().all(|m| *m >= 0.0));
        if state.mass_kg() == 0.0 {
            assert_eq!(state.volume().value(), 0.0);
            assert_eq!(state.sensible_internal_energy().value(), 0.0);
        } else {
            let amounts: f64 = state
                .species_masses()
                .iter()
                .zip(room.gas().source().models())
                .map(|(m, source)| m / source.molar_mass().value())
                .sum();
            close(
                amounts * R * state.temperature().value() / state.volume().value(),
                room.pressure().value(),
                2e-9,
            );
        }
    }
    close(volume, room.volume().value(), 1e-14);
}

#[test]
fn finite_time_relaxation_matches_analytic_solution_without_equilibrium_jump() {
    let m = mechanism(0.0, 3.0, false, 0.0);
    let mut zone = room(&m, 1.0, [2.0; 2], [[1.0, 0.0]; 2]);
    let initial = zone.clone();
    let mut operator = ZoneKinetics::new(&zone, m, controls(), || false).unwrap();
    let report = operator
        .apply_split(
            &mut zone,
            Time::new(0.1),
            [true; 2],
            zone_controls(),
            || false,
        )
        .unwrap();
    let expected_b_fraction = 0.5 * (1.0 - (-0.6_f64).exp());
    let state = zone.layer(Layer::Lower);
    close(
        state.species_masses()[1] / state.mass_kg(),
        expected_b_fraction,
        1.2e-4,
    );
    assert!(state.species_masses()[1] / state.mass_kg() < 0.3);
    close(state.temperature().value(), 2.0, 1e-9);
    assert_eq!(zone.time(), initial.time());
    assert_eq!(report.chemistry.layers[1], None);
    assert_eq!(report.chemistry.chemical_rates[1], LayerRates::default());
    assert_eq!(zone.layer(Layer::Upper), initial.layer(Layer::Upper));
    let detail = report.chemistry.layers[0].unwrap();
    assert!(detail.reactor.entropy_change_j_k > 0.0);
    assert!(detail.reactor.accepted_steps > 0);
    assert!(detail.reactor.maximum_estimated_local_error <= 1.0);
    close(report.physical_energy_residual_j, 0.0, 1e-10);
    check_original_eos(&zone);
}

#[test]
fn both_layers_react_with_variable_calorics_and_exchange_real_interface_work() {
    let m = mechanism(8.0, 2.0, true, 0.0);
    let mut zone = room(&m, 0.4, [1.5, 3.0], [[0.95, 0.05], [0.2, 0.8]]);
    let initial = zone.clone();
    let original_energy = energy_by_original_species(&zone);
    let original_mass = zone.layers[0].mass_kg() + zone.layers[1].mass_kg();
    let mut operator = ZoneKinetics::new(&zone, m, controls(), || false).unwrap();
    let report = operator
        .apply_split(
            &mut zone,
            Time::new(0.15),
            [true; 2],
            zone_controls(),
            || false,
        )
        .unwrap();
    for k in 0..2 {
        let detail = report.chemistry.layers[k].unwrap();
        assert!(detail.reactor.entropy_change_j_k >= -controls().entropy_tolerance_j_k);
        assert!(
            detail.reactor.maximum_relative_element_residual
                <= controls().relative_conservation_tolerance
        );
        close(
            detail.physical_energy_source_j,
            0.0,
            controls().energy_tolerance_j,
        );
        close(detail.mass_source_kg, 0.0, 1e-13);
        for s in 0..2 {
            close(
                zone.layers[k].masses[s],
                initial.layers[k].masses[s]
                    + 0.15 * report.chemistry.chemical_rates[k].species_kg_s[s],
                1e-15,
            );
        }
    }
    assert!(
        (report.chemistry.layers[0]
            .unwrap()
            .isochoric_pressure
            .value()
            - report.chemistry.layers[1]
                .unwrap()
                .isochoric_pressure
                .value())
        .abs()
            > 0.01
    );
    assert!(report.relaxation.interface_work_j[0].abs() > 1e-4);
    assert_eq!(
        report.relaxation.interface_work_j[0],
        -report.relaxation.interface_work_j[1]
    );
    assert_ne!(zone.interface_height(), initial.interface_height());
    close(
        zone.layers[0].mass_kg() + zone.layers[1].mass_kg(),
        original_mass,
        1e-13,
    );
    close(energy_by_original_species(&zone), original_energy, 2e-9);
    close(
        zone.physical_internal_energy().value(),
        original_energy,
        2e-9,
    );
    assert_eq!(zone.time(), initial.time());
    check_original_eos(&zone);
}

#[test]
fn disabled_layer_preserves_species_and_receives_only_original_mechanical_work() {
    let m = mechanism(10.0, 3.0, false, 0.0);
    let mut zone = room(&m, 0.5, [2.0, 3.0], [[1.0, 0.0]; 2]);
    let before = zone.clone();
    let mut operator = ZoneKinetics::new(&zone, m, controls(), || false).unwrap();
    let report = operator
        .apply_split(
            &mut zone,
            Time::new(0.1),
            [false, true],
            zone_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(report.chemistry.layers[0], None);
    assert_eq!(report.chemistry.chemical_rates[0], LayerRates::default());
    assert_eq!(
        zone.layer(Layer::Lower).species_masses(),
        before.layer(Layer::Lower).species_masses()
    );
    assert!(
        zone.layer(Layer::Lower).temperature().value()
            > before.layer(Layer::Lower).temperature().value()
    );
    close(
        zone.layers[0].energy - before.layers[0].energy,
        -report.relaxation.interface_work_j[0],
        1e-9,
    );
    check_original_eos(&zone);
}

#[test]
fn repeated_pressure_trials_reuse_original_state_and_preserve_nonzero_clock() {
    let m = mechanism(4.0, 2.0, false, 0.0);
    let mut accepted = room(&m, 0.5, [2.0, 3.0], [[0.9, 0.1], [0.4, 0.6]]);
    accepted
        .advance(
            Time::new(0.123),
            [LayerRates::default(); 2],
            zone_controls(),
            || false,
        )
        .unwrap();
    let before = accepted.clone();
    let mut operator = ZoneKinetics::new(&accepted, m, controls(), || false).unwrap();
    let first = operator
        .apply_split(
            &mut accepted,
            Time::new(0.07),
            [true; 2],
            zone_controls(),
            || false,
        )
        .unwrap();
    let endpoint = accepted.clone();
    accepted.copy_state_from(&before).unwrap();
    let second = operator
        .apply_split(
            &mut accepted,
            Time::new(0.07),
            [true; 2],
            zone_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(first, second);
    assert_eq!(accepted, endpoint);
    assert_eq!(
        accepted.time().value().to_bits(),
        before.time().value().to_bits()
    );
}

#[test]
fn original_energy_gauge_shift_does_not_change_accepted_physics() {
    let run = |reference| {
        let m = mechanism(5.0, 2.0, true, reference);
        let mut zone = room(&m, 0.3, [2.0, 3.0], [[0.8, 0.2], [0.4, 0.6]]);
        let mut operator = ZoneKinetics::new(&zone, m, controls(), || false).unwrap();
        operator
            .apply_split(
                &mut zone,
                Time::new(0.1),
                [true; 2],
                zone_controls(),
                || false,
            )
            .unwrap();
        zone
    };
    let a = run(0.0);
    let b = run(100.0);
    close(a.pressure().value(), b.pressure().value(), 1e-9);
    close(
        a.interface_height().value(),
        b.interface_height().value(),
        1e-9,
    );
    for layer in [Layer::Lower, Layer::Upper] {
        close(
            a.layer(layer).temperature().value(),
            b.layer(layer).temperature().value(),
            1e-9,
        );
        for (am, bm) in a
            .layer(layer)
            .species_masses()
            .iter()
            .zip(b.layer(layer).species_masses())
        {
            close(*am, *bm, 1e-10);
        }
    }
}

#[test]
fn source_domain_geometry_and_work_refusals_leave_room_unchanged() {
    let m = mechanism(5.0, 2.0, false, 0.0);
    let mut zone = room(&m, 0.5, [2.0, 3.0], [[0.9, 0.1]; 2]);
    let initial = zone.clone();
    assert!(matches!(
        ZoneKinetics::new(&zone, mechanism(5.0, 2.0, false, 1.0), controls(), || false),
        Err(ZoneKineticsError::SourceMismatch)
    ));
    let mut bad = controls();
    bad.maximum_temperature_k = 51.0;
    assert!(matches!(
        ZoneKinetics::new(&zone, m.clone(), bad, || false),
        Err(ZoneKineticsError::Invalid(_))
    ));
    let mut operator = ZoneKinetics::new(&zone, m.clone(), controls(), || false).unwrap();
    assert!(
        operator
            .apply_split(
                &mut zone,
                Time::new(0.0),
                [true; 2],
                zone_controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(zone, initial);
    let mut other = room(
        &mechanism(5.0, 2.0, false, 1.0),
        0.5,
        [2.0, 3.0],
        [[0.9, 0.1]; 2],
    );
    let other_before = other.clone();
    assert_eq!(
        operator.apply_split(
            &mut other,
            Time::new(0.1),
            [true; 2],
            zone_controls(),
            || false
        ),
        Err(ZoneKineticsError::SourceMismatch)
    );
    assert_eq!(other, other_before);
    let mut geometry = zone.clone();
    geometry.floor_area *= 2.0;
    let geometry_before = geometry.clone();
    assert!(matches!(
        operator.apply_split(
            &mut geometry,
            Time::new(0.1),
            [true; 2],
            zone_controls(),
            || false
        ),
        Err(ZoneKineticsError::Zone(ZoneError::Invalid(_)))
    ));
    assert_eq!(geometry, geometry_before);
    let mut bounded = controls();
    bounded.maximum_steps = 1;
    let mut operator = ZoneKinetics::new(&zone, m.clone(), bounded, || false).unwrap();
    assert!(matches!(
        operator.apply_split(
            &mut zone,
            Time::new(0.5),
            [true; 2],
            zone_controls(),
            || false
        ),
        Err(ZoneKineticsError::Chemistry { .. })
    ));
    assert_eq!(zone, initial);
    let mut narrowed = controls();
    narrowed.maximum_temperature_k = 2.5;
    let mut operator = ZoneKinetics::new(&zone, m, narrowed, || false).unwrap();
    assert!(matches!(
        operator.apply_split(
            &mut zone,
            Time::new(0.1),
            [true; 2],
            zone_controls(),
            || false
        ),
        Err(ZoneKineticsError::TemperatureDomain {
            layer: Layer::Upper,
            ..
        })
    ));
    assert_eq!(zone, initial);
}

#[test]
fn pressure_relaxation_leaving_kinetic_interval_refuses_before_publication() {
    let m = mechanism(8.0, 3.0, false, 0.0);
    let upper_t = 3.99;
    // Exact equilibrium of these equal-Cp, equal-molecular-count isomers.
    let upper_a = 1.0 / (1.0 + (8.0_f64 / upper_t).exp());
    let mut zone = room(
        &m,
        0.5,
        [1.0, upper_t],
        [[1.0, 0.0], [upper_a, 1.0 - upper_a]],
    );
    let before = zone.clone();
    let mut narrow = controls();
    narrow.maximum_temperature_k = 4.0;
    let mut operator = ZoneKinetics::new(&zone, m, narrow, || false).unwrap();
    let homogeneous = operator
        .evaluate(&zone, Time::new(0.1), [true; 2], || false)
        .unwrap();
    assert!(
        homogeneous
            .layers
            .iter()
            .flatten()
            .all(|l| l.isochoric_temperature.value() <= 4.0)
    );
    let error = operator.apply_split(
        &mut zone,
        Time::new(0.1),
        [true; 2],
        zone_controls(),
        || false,
    );
    assert!(
        matches!(error, Err(ZoneKineticsError::TemperatureDomain { layer: Layer::Upper, temperature_k }) if temperature_k > 4.0),
        "{error:?}"
    );
    assert_eq!(zone, before);
    operator
        .apply_split(
            &mut zone,
            Time::new(1e-5),
            [true; 2],
            zone_controls(),
            || false,
        )
        .unwrap();
    assert!(zone.layer(Layer::Upper).temperature().value() < 4.0);
    assert_eq!(zone.time(), before.time());
    check_original_eos(&zone);
}

#[test]
fn late_cancellation_keeps_both_layers_and_prepared_operator_reusable() {
    let m = mechanism(3.0, 2.0, false, 0.0);
    let mut zone = room(&m, 0.5, [2.0, 3.0], [[0.8, 0.2]; 2]);
    let before = zone.clone();
    let mut operator = ZoneKinetics::new(&zone, m, controls(), || false).unwrap();
    let mut calls = 0usize;
    let expected = operator
        .apply_split(
            &mut zone,
            Time::new(0.04),
            [true; 2],
            zone_controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let endpoint = zone.clone();
    zone.copy_state_from(&before).unwrap();
    let mut seen = 0usize;
    assert_eq!(
        operator.apply_split(
            &mut zone,
            Time::new(0.04),
            [true; 2],
            zone_controls(),
            || {
                seen += 1;
                seen == calls
            }
        ),
        Err(ZoneKineticsError::Cancelled)
    );
    assert_eq!(zone, before);
    let repeated = operator
        .apply_split(
            &mut zone,
            Time::new(0.04),
            [true; 2],
            zone_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(repeated, expected);
    assert_eq!(zone, endpoint);
}
