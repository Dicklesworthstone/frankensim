use super::*;
use crate::multicomponent::CaloricDomain;
use crate::multicomponent::tests::{bank, source};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn gas(linear: f64, unequal: bool, formation: f64) -> MixtureGas {
    let b_mass = if unequal { 2.0 * R } else { R };
    let b_cp = if unequal { 4.0 } else { 3.5 };
    let b_linear = if unequal { 0.5 * linear } else { linear };
    MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, linear, formation),
            source("B", b_mass, b_cp, b_linear, -2.0 * formation),
        ]),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn room(
    gas: MixtureGas,
    volume: f64,
    fraction: f64,
    temperatures: [f64; 2],
    y: [[f64; 2]; 2],
) -> TwoZone {
    TwoZone::new(
        gas,
        Area::new(0.5 * volume),
        Length::new(2.0),
        Length::new(2.0 * fraction),
        Pressure::new(4.0),
        temperatures.map(Temperature::new),
        y.map(Vec::from),
        || false,
    )
    .unwrap()
}
fn solve() -> ZoneSolve {
    ZoneSolve {
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
        ..ZoneSolve::default()
    }
}
fn admit(room: &TwoZone) {
    let n = room.gas().species_count();
    let mut volumes = 0.0;
    for layer in [Layer::Lower, Layer::Upper] {
        let state = room.layer(layer);
        volumes += state.volume().value();
        assert!(
            state
                .species_masses()
                .iter()
                .all(|m| m.is_finite() && *m >= 0.0)
        );
        if state.mass_kg() == 0.0 {
            assert_eq!(state.volume().value(), 0.0);
            assert_eq!(state.sensible_internal_energy().value(), 0.0);
            continue;
        }
        assert!(state.volume().value() > 0.0);
        let fraction_sum = sum(state.species_masses().iter().map(|m| m / state.mass_kg()));
        assert!((fraction_sum - 1.0).abs() < 1e-14);
        let mut q = [0.0; MAX_COMPONENTS];
        for i in 0..n {
            q[i] = state.species_masses()[i] / state.volume().value();
        }
        q[n + 1] = state.sensible_internal_energy().value() / state.volume().value();
        let original = room.gas().decode(&q[..n + 2], None, &mut || false).unwrap();
        assert!((original.pressure - room.pressure().value()).abs() < 2e-8);
        assert!((original.temperature - state.temperature().value()).abs() < 2e-8);
    }
    assert!((volumes - room.volume().value()).abs() < 1e-14 * room.volume().value());
    assert!((0.0..=room.height().value()).contains(&room.interface_height().value()));
}
fn cfast_error(steps: usize) -> f64 {
    let mut zone = room(gas(0.0, false, 0.0), 3.0, 0.5, [2.0, 3.0], [[0.7, 0.3]; 2]);
    let initial_energy = zone.physical_internal_energy().value();
    let mut upper = LayerRates::default();
    upper.add_heat(Power::new(6.0)).unwrap();
    for _ in 0..steps {
        let report = zone
            .advance(
                Time::new(1.0 / steps as f64),
                [LayerRates::default(), upper],
                solve(),
                || false,
            )
            .unwrap();
        assert_eq!(report.interface_work_j[0], -report.interface_work_j[1]);
        assert!(report.interface_work_j[0] < 0.0);
        assert!(report.physical_energy_residual_j.abs() < 1e-10);
    }
    // CFAST constant-Cp equations give p=p0+(gamma-1)Qt/V. The unheated
    // lower layer is adiabatically compressed: Vl/Vl0=(p0/p)^(1/gamma).
    let expected_p: f64 = 4.0 + 0.4 * 6.0 / 3.0;
    let expected_volume = 1.5 * (4.0 / expected_p).powf(1.0 / 1.4);
    assert!((zone.pressure().value() - expected_p).abs() < 2e-8);
    assert!((zone.physical_internal_energy().value() - initial_energy - 6.0).abs() < 1e-9);
    assert_eq!(zone.time().value(), 1.0);
    admit(&zone);
    (zone.layer(Layer::Lower).volume().value() - expected_volume).abs()
}

#[test]
fn g1_constant_cp_recovers_cfast_pressure_and_refines_adiabatic_interface_motion() {
    let coarse = cfast_error(8);
    let fine = cfast_error(16);
    assert!(coarse < 1e-4);
    assert!(fine < 0.3 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn g1_variable_cp_heating_uses_the_original_caloric_integral() {
    let mut zone = room(gas(0.2, false, 0.0), 3.0, 0.4, [2.0; 2], [[0.7, 0.3]; 2]);
    let mut rates = [LayerRates::default(); 2];
    rates[0].add_heat(Power::new(2.4)).unwrap();
    rates[1].add_heat(Power::new(3.6)).unwrap();
    let report = zone
        .advance(Time::new(1.0), rates, solve(), || false)
        .unwrap();
    // Cv(T)=2.5+0.2T J/(kg K); equal heating per kg adds 1 J/kg.
    let expected = 2.0 + 2.0 / (2.9 + (2.9_f64 * 2.9 + 0.4).sqrt());
    for layer in [Layer::Lower, Layer::Upper] {
        assert!((zone.layer(layer).temperature().value() - expected).abs() < 1e-8);
    }
    assert!((expected - (2.0 + 1.0 / 2.9)).abs() > 1e-3);
    assert!((zone.pressure().value() - 2.0 * expected).abs() < 1e-8);
    assert!((zone.interface_height().value() - 0.8).abs() < 1e-8);
    assert!(report.physical_energy_residual_j.abs() < 1e-9);
    admit(&zone);
}

#[test]
fn g0_internal_enthalpy_transfer_conserves_each_species_and_physical_energy() {
    let mut zone = room(
        gas(0.05, true, 10.0),
        3.0,
        0.55,
        [2.0, 5.0],
        [[0.9, 0.1], [0.2, 0.8]],
    );
    let original = zone.clone();
    let carried = zone
        .carried_rate(Layer::Lower, MassFlowRate::new(0.1))
        .unwrap();
    let donor = zone.layer(Layer::Lower);
    let flow_work = carried.referenced_energy_w / 0.1
        - donor.sensible_internal_energy().value() / donor.mass_kg();
    assert!((flow_work - zone.pressure().value() / donor.density_kg_m3().unwrap()).abs() < 1e-12);
    let report = zone
        .advance(
            Time::new(0.2),
            [carried.opposite(), carried],
            solve(),
            || false,
        )
        .unwrap();
    for i in 0..zone.gas().species_count() {
        let before = original.layers[0].masses[i] + original.layers[1].masses[i];
        let after = zone.layers[0].masses[i] + zone.layers[1].masses[i];
        assert!((after - before).abs() < 1e-14);
    }
    assert!(
        (zone.physical_internal_energy().value() - original.physical_internal_energy().value())
            .abs()
            < 1e-10
    );
    assert!(report.mass_residual_kg.abs() < 1e-14);
    assert!(report.physical_energy_residual_j.abs() < 1e-10);
    admit(&zone);
}

#[test]
fn g0_entrainment_and_heat_create_an_absent_upper_layer_without_species_seeding() {
    let mut zone = room(
        gas(0.05, true, 0.0),
        3.0,
        1.0,
        [2.0; 2],
        [[1.0, 0.0], [0.0, 1.0]],
    );
    assert_eq!(zone.layer(Layer::Upper).species_masses(), &[0.0, 0.0]);
    assert_eq!(zone.layer(Layer::Upper).volume().value(), 0.0);
    let carried = zone
        .carried_rate(Layer::Lower, MassFlowRate::new(0.15))
        .unwrap();
    let mut upper = carried;
    upper.add_heat(Power::new(3.0)).unwrap();
    let report = zone
        .advance(Time::new(0.2), [carried.opposite(), upper], solve(), || {
            false
        })
        .unwrap();
    assert!((zone.layer(Layer::Upper).mass_kg() - 0.03).abs() < 1e-14);
    assert_eq!(zone.layer(Layer::Upper).species_masses()[1], 0.0);
    assert_eq!(zone.layer(Layer::Lower).species_masses()[1], 0.0);
    assert!(zone.layer(Layer::Upper).temperature().value() > 2.0);
    assert!(zone.interface_height().value() < zone.height().value());
    assert!(report.physical_energy_residual_j.abs() < 1e-9);
    admit(&zone);
}

#[test]
fn g4_late_cancellation_iteration_limit_and_depletion_preserve_every_inventory() {
    let mut zone = room(gas(0.2, true, 0.0), 3.0, 0.4, [2.0, 4.0], [[0.7, 0.3]; 2]);
    let original = zone.clone();
    let mut upper = LayerRates::default();
    upper.add_heat(Power::new(6.0)).unwrap();
    let rates = [LayerRates::default(), upper];
    let mut reference = zone.clone();
    let mut polls = 0;
    reference
        .advance(Time::new(0.1), rates, solve(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    let error = zone
        .advance(Time::new(0.1), rates, solve(), || {
            calls += 1;
            calls >= polls
        })
        .unwrap_err();
    assert_eq!(error, ZoneError::Cancelled);
    assert_eq!(zone, original);
    let error = zone
        .advance(
            Time::new(0.1),
            rates,
            ZoneSolve {
                maximum_iterations: 1,
                ..solve()
            },
            || false,
        )
        .unwrap_err();
    assert_eq!(error, ZoneError::IterationLimit);
    assert_eq!(zone, original);
    let mut drain = LayerRates::default();
    drain.species_kg_s[0] = -1e6;
    assert!(matches!(
        zone.advance(Time::new(1.0), [drain, upper], solve(), || false),
        Err(ZoneError::SpeciesDepletion { .. })
    ));
    assert_eq!(zone, original);
    zone.advance(Time::new(0.1), rates, solve(), || false)
        .unwrap();
    assert_eq!(zone, reference);
    admit(&zone);
}

#[test]
fn g0_interior_states_are_not_rejected_by_caloric_bracket_roundoff() {
    for volume in [0.3, 0.7, 1.3, 3.5, 16.6992] {
        let mut zone = room(
            gas(0.1, true, 0.0),
            volume,
            0.47,
            [2.3, 4.1],
            [[0.3, 0.7], [0.8, 0.2]],
        );
        let initial = zone.clone();
        let result = zone.advance(Time::new(0.01), [LayerRates::default(); 2], solve(), || {
            false
        });
        assert!(result.is_ok(), "volume={volume}: {result:?}");
        assert!((zone.pressure().value() - initial.pressure().value()).abs() < 1e-8);
        assert!(
            (zone.interface_height().value() - initial.interface_height().value()).abs() < 1e-8
        );
        admit(&zone);
    }
}

#[test]
fn g0_initialization_refuses_unrepresentable_extensive_inventories() {
    let result = TwoZone::new(
        gas(0.0, false, 0.0),
        Area::new(1e300),
        Length::new(1.0),
        Length::new(0.5),
        Pressure::new(1e20),
        [Temperature::new(2.0); 2],
        [vec![0.7, 0.3], vec![0.7, 0.3]],
        || false,
    );
    assert!(result.is_err());
}

#[test]
fn g1_split_stream_and_heat_preserve_clock_and_original_physical_energy() {
    let mut zone = room(gas(0.2, true, 30.0), 3.0, 0.4, [2.0, 3.0], [[0.7, 0.3]; 2]);
    zone.advance(Time::new(0.3), [LayerRates::default(); 2], solve(), || {
        false
    })
    .unwrap();
    let before = zone.clone();
    let dt = Time::new(0.1);
    let mut rates = [LayerRates::default(); 2];
    rates[0] = zone
        .inflow(MassFlowRate::new(0.2), &[0.1, 0.9], Temperature::new(4.0))
        .unwrap();
    rates[1] = zone
        .carried_rate(Layer::Upper, MassFlowRate::new(0.1))
        .unwrap()
        .opposite();
    rates[1].add_heat(Power::new(2.0)).unwrap();
    let mut advanced = before.clone();
    let reference = advanced.advance(dt, rates, solve(), || false).unwrap();
    let report = zone.apply_split(dt, rates, solve(), || false).unwrap();
    assert_eq!(report, reference);
    assert_eq!(zone.layers, advanced.layers);
    assert_eq!(zone.pressure(), advanced.pressure());
    assert_eq!(
        zone.time().value().to_bits(),
        before.time().value().to_bits()
    );
    let physical_rate: f64 = rates
        .iter()
        .map(|rate| {
            rate.referenced_energy_w
                + rate
                    .species_kg_s
                    .iter()
                    .zip(zone.gas().energy_references())
                    .map(|(m, e)| m * e)
                    .sum::<f64>()
        })
        .sum();
    assert!(
        (zone.physical_internal_energy().value()
            - before.physical_internal_energy().value()
            - dt.value() * physical_rate)
            .abs()
            < 1e-9
    );
    assert_ne!(zone.interface_height(), before.interface_height());
    let accepted = zone.clone();
    let mut drain = rates;
    drain[0].species_kg_s[0] = -1e6;
    assert!(matches!(
        zone.apply_split(dt, drain, solve(), || false),
        Err(ZoneError::SpeciesDepletion { .. })
    ));
    assert_eq!(zone, accepted);
    assert_eq!(
        zone.apply_split(dt, rates, solve(), || true),
        Err(ZoneError::Cancelled)
    );
    assert_eq!(zone, accepted);
}
