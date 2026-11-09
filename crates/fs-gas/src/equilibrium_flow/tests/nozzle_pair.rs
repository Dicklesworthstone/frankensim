//! G1/G3 two-volume conservation, reacting receiving and reversible orientation.
use super::*;
use super::super::nozzle::{
    EquilibriumNozzleControl,
    pair::{EquilibriumGasPair, EquilibriumGasPairControl, EquilibriumGasPairDirection,
        EquilibriumGasPairError, EquilibriumTankInitialState},
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::equilibrium::blowdown::EquilibriumBlowdownSolve;
use fs_qty::Area;
use fs_thermochem::equilibrium::composition::volume::{
    EquilibriumInternalEnergySolveV1, EquilibriumVolumeSolveV1,
};
use std::sync::Arc;

const PAIR_AREA: f64 = 1e-5;

fn pair_controls() -> EquilibriumGasPairControl {
    let sp = controls().isentropic;
    let mut nozzle = EquilibriumNozzleControl::new(sp, Pressure::new(100.0));
    nozzle.relative_sonic_tolerance = 1e-10;
    nozzle.relative_pressure_tolerance = 1e-10;
    EquilibriumGasPairControl {
        reservoir: EquilibriumBlowdownSolve {
            minimum_pressure_pa: 100.0,
            relative_density_tolerance: 1e-12,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-5,
            maximum_density_iterations: 80,
        },
        nozzle,
        receiver: EquilibriumInternalEnergySolveV1 {
            volume: EquilibriumVolumeSolveV1 {
                equilibrium: sp.equilibrium,
                relative_volume_tolerance: 1e-12,
                maximum_pressure_iterations: 80,
            },
            minimum_temperature_k: sp.minimum_temperature_k,
            maximum_temperature_k: sp.maximum_temperature_k,
            energy_tolerance_j: 1e-7,
            maximum_temperature_iterations: 80,
        },
        entropy_tolerance_j_k: 1e-8,
        relative_element_tolerance: 1e-10,
        maximum_outgoing_fraction: 0.25,
        maximum_mass_iterations: 64,
    }
}
fn initial(temperature: f64, pressure: f64, volume: f64) -> EquilibriumTankInitialState {
    EquilibriumTankInitialState {
        volume_m3: volume, temperature: Temperature::new(temperature),
        pressure: Pressure::new(pressure),
    }
}
fn pair(
    chemistry: &IdealGasEquilibriumV1,
    left: EquilibriumTankInitialState, right: EquilibriumTankInitialState,
    controls: EquilibriumGasPairControl,
) -> EquilibriumGasPair {
    EquilibriumGasPair::new(
        Arc::new(chemistry.clone()),
        vec![(id(if chemistry.models().len() == 1 { "A" } else { "A2" }), Amount::new(1.0))],
        left, right, controls, || false,
    ).unwrap()
}
fn total_energy(pair: &EquilibriumGasPair) -> f64 {
    pair.tanks().iter().map(|tank| tank.internal_energy_j()).sum()
}
fn total_mass(pair: &EquilibriumGasPair) -> f64 {
    pair.tanks().iter().map(|tank| tank.mass_kg()).sum()
}

#[test]
fn g1_frozen_pair_receives_integrated_donor_energy_and_has_symmetric_direction() {
    let chemistry = single();
    let hot = initial(T0, 2.0 * P0, 0.01);
    let low = initial(700.0, 0.3 * P0, 0.02);
    let first = pair(&chemistry, hot, low, pair_controls());
    let swapped = pair(&chemistry, low, hot, pair_controls());
    let a = first.advance(Area::new(PAIR_AREA), 0.01, || false).unwrap();
    let b = swapped.advance(Area::new(PAIR_AREA), 0.01, || false).unwrap();
    assert_eq!(a.step().direction, Some(EquilibriumGasPairDirection::LeftToRight));
    assert_eq!(b.step().direction, Some(EquilibriumGasPairDirection::RightToLeft));
    let before = first.tanks();
    let next = a.endpoint().tanks();
    let gamma = 4.0 / 3.0;
    let cv = 3.0 * R / 0.01;
    let donor_t = T0 * (next[0].mass_kg() / before[0].mass_kg()).powf(gamma - 1.0);
    let heat = cv * (before[0].mass_kg() * T0 - next[0].mass_kg() * donor_t);
    let receiver_t = (before[1].mass_kg() * 700.0 + heat / cv) / next[1].mass_kg();
    close(next[0].isentropic_sample().temperature.value(), donor_t, 1e-6);
    close(next[1].isentropic_sample().temperature.value(), receiver_t, 1e-6);
    close(a.step().donor.as_ref().unwrap().exhausted_enthalpy_j, heat, 1e-5);
    close(total_mass(a.endpoint()), total_mass(&first), 1e-13);
    close(total_energy(a.endpoint()), total_energy(&first), 1e-5);
    assert!(a.step().entropy_production_j_k > 0.0);
    for index in 0..2 {
        close(next[index].mass_kg(), b.endpoint().tanks()[1-index].mass_kg(), 1e-13);
        close(next[index].isentropic_sample().temperature.value(),
            b.endpoint().tanks()[1-index].isentropic_sample().temperature.value(), 1e-6);
    }
    let flow = a.step().midpoint_nozzle.as_ref().unwrap();
    assert!(flow.back_pressure_pa > low.pressure.value());
    assert!(flow.back_pressure_pa < next[1].isentropic_sample().pressure.value());
    assert!(a.step().energy_change_since_reference_j.abs() <= pair_controls().reservoir.energy_tolerance_j);
}

#[test]
fn g3_choked_pair_keeps_second_order_time_and_closed_pair_energy() {
    let chemistry = single();
    let original = pair(&chemistry, initial(T0, 2.0*P0, 0.01),
        initial(T0, 0.2*P0, 0.02), pair_controls());
    let gamma = 4.0_f64 / 3.0;
    let r = R / 0.01;
    let rate = PAIR_AREA * 2.0*P0 / (r*T0).sqrt() * gamma.sqrt()
        * (2.0/(gamma+1.0)).powf((gamma+1.0)/(2.0*(gamma-1.0)));
    let duration = 0.2;
    let exponent = 0.5*(gamma-1.0);
    let m0 = original.tanks()[0].mass_kg();
    let exact = m0 * (1.0 + exponent*rate*duration/m0).powf(-1.0/exponent);
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let mut current = original.clone();
        for _ in 0..steps {
            let candidate = current.advance(
                Area::new(PAIR_AREA), duration/steps as f64, || false,
            ).unwrap();
            assert!(candidate.step().midpoint_nozzle.as_ref().unwrap().choked);
            assert!(candidate.step().mass_quadrature_residual_kg.abs()
                <= pair_controls().reservoir.mass_tolerance_kg);
            assert!(candidate.step().energy_change_since_reference_j.abs()
                <= pair_controls().reservoir.energy_tolerance_j);
            current = candidate.into_parts().0;
        }
        errors.push((current.tanks()[0].mass_kg()-exact).abs());
        close(total_energy(&current), total_energy(&original), 1e-5);
    }
    for e in errors.windows(2) {
        assert!(e[0]/e[1] > 3.5 && e[0]/e[1] < 4.5, "{errors:?}");
    }
}

#[test]
fn g1_reacting_receiver_changes_composition_while_both_element_and_energy_accounts_close() {
    let chemistry = dissociation();
    let original = pair(&chemistry, initial(1600.0, 2.0*P0, 0.01),
        initial(600.0, 0.3*P0, 0.02), pair_controls());
    let mut current = original.clone();
    let mut last = None;
    for _ in 0..8 {
        let candidate = current.advance(Area::new(PAIR_AREA), 0.002, || false).unwrap();
        assert!(candidate.step().entropy_production_j_k >= -pair_controls().entropy_tolerance_j_k);
        assert!(candidate.step().maximum_relative_element_residual
            <= pair_controls().relative_element_tolerance);
        for tank in candidate.endpoint().tanks() {
            check_products(tank.thermodynamics(), &chemistry);
        }
        let (next, step) = candidate.into_parts();
        current = next;
        last = Some(step);
    }
    let a = |pair: &EquilibriumGasPair| pair.tanks()[1].thermodynamics()
        .equilibrium().species().iter().find(|species| species.species == id("A"))
        .unwrap().mole_fraction;
    assert!(a(&current) > a(&original));
    close(total_energy(&current), total_energy(&original), 1e-5);
    close(total_mass(&current), total_mass(&original), 1e-13);
    assert!(last.unwrap().entropy_change_since_reference_j_k > 0.0);
    let elements = |pair: &EquilibriumGasPair| pair.tanks().iter().map(|tank| {
        tank.thermodynamics().equilibrium().species().iter().map(|species| {
            species.amount.value() * if species.species == id("A2") { 2.0 } else { 1.0 }
        }).sum::<f64>() * tank.physical_amount_scale()
    }).sum::<f64>();
    close(elements(&current), elements(&original), 1e-10);
}

#[test]
fn g3_pair_formation_reference_does_not_erase_small_received_energy() {
    let plain = single();
    let shifted = owner(vec![model("A", 0.01, 4.0, 1e12, 0.0)], vec![1]);
    let mut settings = pair_controls();
    settings.reservoir.energy_tolerance_j = 1e-4;
    let left = initial(T0, 2.0*P0, 0.01);
    let right = initial(900.0, 0.4*P0, 0.02);
    let a = pair(&plain, left, right, settings)
        .advance(Area::new(PAIR_AREA), 0.001, || false).unwrap();
    let b = pair(&shifted, left, right, settings)
        .advance(Area::new(PAIR_AREA), 0.001, || false).unwrap();
    for i in 0..2 {
        close(a.endpoint().tanks()[i].mass_kg(), b.endpoint().tanks()[i].mass_kg(), 1e-11);
        close(a.endpoint().tanks()[i].isentropic_sample().temperature.value(),
            b.endpoint().tanks()[i].isentropic_sample().temperature.value(), 2e-4);
    }
    assert!(b.endpoint().tanks()[1].isentropic_sample().temperature.value() > 900.0);
    assert!(b.step().energy_change_since_reference_j.abs() <= settings.reservoir.energy_tolerance_j);
}

#[test]
fn g4_pair_rest_crossing_work_exhaustion_and_final_cancellation_are_atomic() {
    let chemistry = single();
    let same = pair(&chemistry, initial(T0, P0, 0.01), initial(700.0, P0, 0.02), pair_controls());
    let rest = same.advance(Area::new(PAIR_AREA), 1.0, || false).unwrap();
    assert_eq!(rest.endpoint(), &same);
    assert!(rest.step().midpoint_nozzle.is_none());
    let original = pair(&chemistry, initial(T0, 1.01*P0, 0.01),
        initial(T0, P0, 0.01), pair_controls());
    let closed = original.advance(Area::new(0.0), 1.0, || false).unwrap();
    assert_eq!(closed.endpoint(), &original);
    let crossing = original.advance(Area::new(PAIR_AREA), 1.0, || false);
    assert!(matches!(crossing, Err(EquilibriumGasPairError::PressureCrossing { .. })
        | Err(EquilibriumGasPairError::NoConvergence { .. })), "{crossing:?}");
    let mut limited = pair_controls();
    limited.maximum_mass_iterations = 1;
    let limited = pair(&chemistry, initial(T0, 2.0*P0, 0.01),
        initial(T0, 0.2*P0, 0.02), limited);
    assert!(matches!(limited.advance(Area::new(PAIR_AREA), 0.01, || false),
        Err(EquilibriumGasPairError::NoConvergence { .. })));
    let source = pair(&chemistry, initial(T0, 2.0*P0, 0.01),
        initial(T0, 0.2*P0, 0.02), pair_controls());
    let before = source.clone();
    let mut polls = 0;
    let successful = source.advance(Area::new(PAIR_AREA), 0.01, || { polls += 1; false }).unwrap();
    assert!(polls > 10);
    let mut again = 0;
    let cancelled = source.advance(Area::new(PAIR_AREA), 0.01, || { again += 1; again == polls });
    assert!(matches!(cancelled, Err(EquilibriumGasPairError::Cancelled)));
    assert_eq!(source, before);
    assert_eq!(source.advance(Area::new(PAIR_AREA), 0.01, || false).unwrap(), successful);
}

#[test]
fn g4_pair_does_not_reuse_a_receivers_old_outflow_entropy() {
    let chemistry = dissociation();
    let original = pair(&chemistry, initial(1600.0, 2.0*P0, 0.01),
        initial(600.0, 0.3*P0, 0.02), pair_controls());
    let result = original.advance(Area::new(PAIR_AREA), 0.002, || false).unwrap();
    let receiver = &result.endpoint().tanks()[1];
    assert_eq!(receiver.isentropic_sample().entropy_change_j_k, 0.0);
    assert_ne!(receiver.caloric_anchor(), original.tanks()[1].caloric_anchor());
    let original_copy = receiver.clone();
    let mut prepared = receiver.prepare_depletion(pair_controls().reservoir, || false).unwrap();
    let foreign = fs_thermochem::equilibrium::composition::isentrope::EquilibriumIsentropeV1::new(
        &chemistry, Temperature::new(900.0), receiver.isentropic_sample().pressure,
        &[(id("A2"), Amount::new(1.0))], controls().isentropic, || false,
    ).unwrap();
    *prepared.isentrope() = foreign;
    assert!(prepared.sample_remaining_mass(0.99*receiver.mass_kg(), || false).is_err());
    assert!(prepared.finish_remaining_mass(0.99*receiver.mass_kg(), || false).is_err());
    assert_eq!(receiver, &original_copy);
}
