use super::*;
use fs_material::phase::water::WaterPhase;
use fs_qty::Pressure;
use fs_thermochem::water::if97::{If97Region, If97SaturationPair};

fn controls(minimum_temperature_k: f64, maximum_temperature_k: f64) -> WaterCellControls {
    WaterCellControls {
        equilibrium: IsochoricWaterControls {
            minimum_temperature_k,
            maximum_temperature_k,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 100e6,
            volume_tolerance_m3_kg: 1e-14,
            internal_energy_tolerance_j_kg: 1e-6,
            maximum_pressure_iterations: 96,
            maximum_temperature_iterations: 96,
        },
        maximum_relative_mass_exchange: 0.25,
        maximum_inflow_parcels: 8,
        mass_tolerance_kg: 1e-13,
        energy_tolerance_j: 1e-7,
        volume_tolerance_m3: 1e-13,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (+/- {tolerance})"
    );
}
fn compressed_cell() -> (EquilibriumWaterCell, If97State) {
    let source = If97State::at_pressure_temperature(
        Pressure::new(100_000.0),
        Temperature::new(300.0),
        If97Region::Region1,
    )
    .unwrap();
    let cell = EquilibriumWaterCell::new(
        source.specific_volume().value(),
        1.0,
        source.specific_internal_energy().value(),
        0.0,
        controls(299.9, 300.1),
        || false,
    )
    .unwrap();
    (cell, source)
}
fn receiver_cell() -> EquilibriumWaterCell {
    EquilibriumWaterCell::at_temperature(
        0.005,
        0.05,
        Temperature::new(300.0),
        0.0,
        controls(299.0, 301.0),
        || false,
    )
    .unwrap()
}

#[test]
fn conserved_mass_volume_energy_recover_original_finite_phase_inventory() {
    for (temperature, quality) in [(300.0, 0.001), (500.0, 0.4)] {
        let pair = If97SaturationPair::at_temperature(Temperature::new(temperature)).unwrap();
        let (liquid, vapor) = (pair.liquid(), pair.vapor());
        let mass = 0.2;
        let specific_volume = (1.0 - quality) * liquid.specific_volume().value()
            + quality * vapor.specific_volume().value();
        let specific_energy = (1.0 - quality) * liquid.specific_internal_energy().value()
            + quality * vapor.specific_internal_energy().value();
        let volume = mass * specific_volume;
        let energy = mass * specific_energy;
        let cell = EquilibriumWaterCell::new(
            volume,
            mass,
            energy,
            0.0,
            controls(temperature - 1.0, temperature + 1.0),
            || false,
        )
        .unwrap();
        let state = cell.state();
        assert_eq!(state.mass_kg(), mass);
        assert_eq!(state.internal_energy_j(), energy);
        near(state.source_state().temperature_k(), temperature, 1e-8);
        near(
            state.source_state().pressure_pa(),
            pair.pressure().value(),
            2e-3,
        );
        near(state.phase_masses_kg()[1], mass * quality, 1e-12);
        near(state.phase_masses_kg().iter().sum(), mass, 1e-15);
        near(
            state.phase_volumes_m3()[0],
            mass * (1.0 - quality) * liquid.specific_volume().value(),
            1e-12,
        );
        near(
            state.phase_volumes_m3()[1],
            mass * quality * vapor.specific_volume().value(),
            1e-12,
        );
        near(
            state.vapor_volume_fraction(),
            quality * vapor.specific_volume().value() / specific_volume,
            1e-10,
        );
        near(
            state.source_internal_energy_j() - energy,
            state.source_energy_residual_j(),
            0.0,
        );
        assert!(state.source_energy_residual_j().abs() <= cell.controls().energy_tolerance_j);
        assert!(state.source_volume_residual_m3().abs() <= cell.controls().volume_tolerance_m3);
        assert!(state.source_state().saturation_pair().is_some());
        assert!(state.source_state().single_phase_source().is_none());
    }
}

#[test]
fn withdrawal_cavitates_and_source_owned_refill_condenses_with_enthalpy_flow_work() {
    let (mut cell, reservoir) = compressed_cell();
    let initial = cell.state();
    assert_eq!(initial.source_state().phase(), WaterPhase::CompressedLiquid);
    assert_eq!(initial.phase_masses_kg()[1], 0.0);
    let out_mass = 0.01;
    let donor_h = initial.source_state().specific_enthalpy_j_kg();
    let drained = cell.advance(&[], out_mass, 0.0, 1.0, || false).unwrap();
    let parcel = drained.outgoing.unwrap();
    assert_eq!(parcel.mass_kg(), out_mass);
    assert_eq!(parcel.enthalpy_j(), out_mass * donor_h);
    assert_eq!(drained.state.mass_kg(), initial.mass_kg() - out_mass);
    assert_eq!(
        drained.state.internal_energy_j(),
        initial.internal_energy_j() - out_mass * donor_h
    );
    assert_eq!(drained.state.volume_m3(), initial.volume_m3());
    assert_eq!(drained.boundary_work_j, 0.0);
    assert_eq!(
        drained.state.source_state().phase(),
        WaterPhase::LiquidVapor
    );
    assert!(drained.state.phase_masses_kg()[1] > 0.0);
    assert!(
        drained.state.source_state().pressure_pa() < 0.1 * initial.source_state().pressure_pa()
    );
    // Pressure-flow work leaves with h, even though the physical boundary does
    // not move. Subtracting only parcel internal energy would miss this term.
    let flow_work = out_mass
        * initial.source_state().pressure_pa()
        * initial.source_state().specific_volume_m3_kg();
    near(
        parcel.enthalpy_j() - out_mass * initial.source_state().specific_internal_energy_j_kg(),
        flow_work,
        1e-9,
    );
    assert!(flow_work > 0.9);
    let incoming = WaterParcel::new(out_mass, WaterParcelSource::Single(reservoir)).unwrap();
    let refilled = cell.advance(&[incoming], 0.0, 0.0, 2.0, || false).unwrap();
    assert_eq!(refilled.state.mass_kg(), initial.mass_kg());
    near(
        refilled.state.internal_energy_j(),
        initial.internal_energy_j(),
        1e-7,
    );
    assert_eq!(
        refilled.state.source_state().phase(),
        WaterPhase::CompressedLiquid
    );
    assert_eq!(refilled.state.phase_masses_kg()[1], 0.0);
    near(
        refilled.state.source_state().pressure_pa(),
        initial.source_state().pressure_pa(),
        0.1,
    );
    for step in [drained, refilled] {
        assert!(step.mass_residual_kg.abs() <= cell.controls().mass_tolerance_kg);
        assert!(step.energy_residual_j.abs() <= cell.controls().energy_tolerance_j);
        near(
            step.source_first_law_residual_j,
            step.energy_residual_j + step.state.source_energy_residual_j()
                - step.initial.source_energy_residual_j(),
            2e-10,
        );
    }
}

#[test]
fn paired_transfer_keeps_both_finite_water_inventories_and_total_energy() {
    let (mut donor, _) = compressed_cell();
    let mut receiver = receiver_cell();
    let before = [donor.state(), receiver.state()];
    let change =
        transfer_water(&mut donor, &mut receiver, 0.001, [0.2, -0.1], 0.5, || false).unwrap();
    assert_eq!(change.parcel, change.donor.outgoing.unwrap());
    assert_eq!(
        change.receiver.incoming_enthalpy_j,
        change.parcel.enthalpy_j()
    );
    near(
        donor.state().mass_kg() + receiver.state().mass_kg(),
        before[0].mass_kg() + before[1].mass_kg(),
        1e-13,
    );
    near(
        donor.state().internal_energy_j() + receiver.state().internal_energy_j()
            - before[0].internal_energy_j()
            - before[1].internal_energy_j(),
        0.1,
        1e-7,
    );
    assert_eq!(donor.state().time_s(), 0.5);
    assert_eq!(receiver.state().time_s(), 0.5);
    assert_eq!(donor.state().volume_m3(), before[0].volume_m3());
    assert_eq!(receiver.state().volume_m3(), before[1].volume_m3());
    near(change.total_mass_residual_kg, 0.0, 1e-13);
    near(change.total_energy_residual_j, 0.0, 2e-7);
}

#[test]
fn late_cancellation_and_receiver_failure_publish_neither_candidate() {
    let (donor, _) = compressed_cell();
    let receiver = receiver_cell();
    let (mut accepted_donor, mut accepted_receiver) = (donor, receiver);
    let mut polls = 0;
    transfer_water(
        &mut accepted_donor,
        &mut accepted_receiver,
        0.001,
        [0.0; 2],
        1.0,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 4);
    let (mut cancelled_donor, mut cancelled_receiver) = (donor, receiver);
    let mut visits = 0;
    let cancelled = transfer_water(
        &mut cancelled_donor,
        &mut cancelled_receiver,
        0.001,
        [0.0; 2],
        1.0,
        || {
            visits += 1;
            visits == polls
        },
    );
    assert_eq!(cancelled, Err(WaterCellError::Cancelled));
    assert_eq!(visits, polls);
    assert_eq!(cancelled_donor, donor);
    assert_eq!(cancelled_receiver, receiver);
    transfer_water(
        &mut cancelled_donor,
        &mut cancelled_receiver,
        0.001,
        [0.0; 2],
        1.0,
        || false,
    )
    .unwrap();
    assert_eq!(cancelled_donor, accepted_donor);
    assert_eq!(cancelled_receiver, accepted_receiver);
    let mut rejected_donor = donor;
    let mut rejected_receiver = receiver;
    rejected_receiver.controls.maximum_relative_mass_exchange = 0.001;
    let before_receiver = rejected_receiver;
    assert!(matches!(
        transfer_water(
            &mut rejected_donor,
            &mut rejected_receiver,
            0.001,
            [0.0; 2],
            1.0,
            || false
        ),
        Err(WaterCellError::ExcessiveMassExchange { .. })
    ));
    assert_eq!(rejected_donor, donor);
    assert_eq!(rejected_receiver, before_receiver);
}

#[test]
fn conservative_energy_never_accumulates_source_inverse_error_and_gross_flow_is_bounded() {
    let mut cfg = controls(399.0, 401.0);
    cfg.energy_tolerance_j = 0.01;
    cfg.equilibrium.internal_energy_tolerance_j_kg = 0.2;
    let mut cell = EquilibriumWaterCell::at_temperature(
        0.005,
        0.05,
        Temperature::new(400.0),
        0.0,
        cfg,
        || false,
    )
    .unwrap();
    let original = cell;
    let mut ledger = cell.state().internal_energy_j();
    for step in 1..=20 {
        let heat = if step % 2 == 0 {
            -0.123456789
        } else {
            0.123456789
        };
        ledger += heat;
        cell.advance(&[], 0.0, heat, f64::from(step), || false)
            .unwrap();
        assert_eq!(cell.state().internal_energy_j().to_bits(), ledger.to_bits());
        assert!(cell.state().source_energy_residual_j().abs() <= cfg.energy_tolerance_j);
    }
    assert_eq!(
        cell.state().internal_energy_j(),
        original.state().internal_energy_j()
    );
    let before = cell;
    let parcel = WaterParcel::new(
        0.15 * cell.state().mass_kg(),
        WaterParcelSource::Equilibrium(cell.state().source_state()),
    )
    .unwrap();
    let withdrawn = parcel.mass_kg();
    assert!(matches!(
        cell.advance(&[parcel], withdrawn, 0.0, 21.0, || false),
        Err(WaterCellError::ExcessiveMassExchange { .. })
    ));
    assert_eq!(cell, before);
}
