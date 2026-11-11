//! Adaptive finite-wall regression fixtures.
//! Synthetic capacities and gases exercise integration/accounting, not validation.

use super::*;
use super::super::tests::{channel, endpoints, link, material};
use super::super::super::super::tests::{controls, gas, near};
use crate::thermal::transport::volumetric::timeline::FixedVolumeThermalLaw;
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportPaths;
use std::collections::BTreeMap;

fn fixture() -> (
    crate::thermal::transport::volumetric::PreparedVolumetricThermalTransport,
    crate::thermochemical_gas::energy::parcel::transfer::species::effusion::channel::DiffuseMolecularChannel,
    BTreeMap<u64, ClosedGasParcel>,
) {
    let source = material(1.0e6, 700.0, [200.0, 1800.0]);
    let path = channel(2.0e-8);
    let gases = BTreeMap::from([
        (0, gas(650.0, 160_000.0, 0.0)),
        (1, gas(300.0, 80_000.0, 0.0)),
    ]);
    (source, path, gases)
}

fn accuracy(dt: f64, relative: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: dt,
        minimum_step_s: dt * 1.0e-5,
        maximum_step_s: dt,
        maximum_attempts: 4096,
        maximum_accepted_intervals: 1024,
        relative_tolerance: relative,
        temperature_tolerance_k: 1.0e-6,
        pressure_tolerance_pa: 1.0e-5,
        species_mass_tolerance_kg: 1.0e-16,
        stream_energy_tolerance_j: 1.0e-9,
    }
}

#[test]
fn adaptive_finite_wall_converges_to_small_fixed_steps_and_closes_window_energy() {
    let (source, path, gases) = fixture();
    let law = FixedVolumeThermalLaw::try_new(&source, 2800.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 2.0e-7, 500.0).unwrap(),
    )]);
    let wall_temperature_k = walls[&9].temperature_k();
    let links = [link(&path)];
    let solve = endpoints(&gases);
    let physical = controls();
    let duration = 0.02;

    let adaptive = advance_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &links,
        duration,
        &solve,
        physical,
        accuracy(duration, 1.0e-5),
        || false,
    )
    .unwrap();
    assert!(adaptive.step().accepted_intervals > 0);
    assert_eq!(
        adaptive.step().accepted_steps.len(),
        2 * adaptive.step().accepted_intervals
    );
    assert!(adaptive.step().maximum_accepted_error <= 1.0);
    assert!(
        adaptive.step().balance.energy_residual_j.abs()
            <= physical.exchange.energy_tolerance_j
    );
    assert!(adaptive.step().balance.total_entropy_change_j_k > 0.0);

    let mut reference_gas = gases.clone();
    let mut reference_walls = walls.clone();
    let n = 128;
    for _ in 0..n {
        let next = super::super::step_candidate(
            &reference_gas,
            &reference_walls,
            GasTransportPaths::default(),
            &links,
            duration / f64::from(n),
            &solve,
            physical,
            || false,
        )
        .unwrap();
        let (g, w, _) = next.into_parts();
        reference_gas = g;
        reference_walls = w;
    }

    for id in [0, 1] {
        near(
            adaptive.parcels()[&id].state().state().temperature,
            reference_gas[&id].state().state().temperature,
            2.0e-3,
        );
        near(
            adaptive.parcels()[&id].state().state().pressure,
            reference_gas[&id].state().state().pressure,
            0.5,
        );
    }
    near(
        adaptive.walls()[&9].temperature_k(),
        reference_walls[&9].temperature_k(),
        2.0e-3,
    );
    assert_eq!(gases[&0].state().state().temperature, 650.0);
    // The original wall is unchanged. try_new_isochoric stores the law's own
    // inverse of the stored energy (T to e to T), so the requested 500 K comes
    // back within rounding, not necessarily bit for bit.
    assert_eq!(walls[&9].temperature_k(), wall_temperature_k);
    assert!((wall_temperature_k - 500.0).abs() <= 2.0 * (500.0_f64.next_up() - 500.0));
}

#[test]
fn wall_accuracy_is_not_hidden_by_small_gas_state_error() {
    let (source, path, gases) = fixture();
    let law = FixedVolumeThermalLaw::try_new(&source, 1200.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 2.0e-9, 500.0).unwrap(),
    )]);
    let links = [link(&path)];
    let duration = 1.0e-4;
    let loose = advance_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &links,
        duration,
        &endpoints(&gases),
        controls(),
        EffusionAdvanceControls {
            temperature_tolerance_k: 1.0e6,
            pressure_tolerance_pa: 1.0e12,
            species_mass_tolerance_kg: 1.0,
            stream_energy_tolerance_j: 1.0e-12,
            ..accuracy(duration, 1.0e-7)
        },
        || false,
    )
    .unwrap();
    assert!(
        loose.step().rejected_intervals > 0 || loose.step().accepted_intervals > 1,
        "wall energy/temperature must participate in the adaptive estimate"
    );
}

#[test]
fn minimum_step_refusal_keeps_original_gas_and_wall_inventories_immutable() {
    let (source, path, gases) = fixture();
    let original_gas = gases.clone();
    let law = FixedVolumeThermalLaw::try_new(&source, 2800.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1.0e-9, 500.0).unwrap(),
    )]);
    let original_temperature = walls[&9].temperature_k();
    let links = [link(&path)];
    let dt = 1.0e-3;
    let strict = EffusionAdvanceControls {
        initial_step_s: dt,
        minimum_step_s: dt,
        maximum_step_s: dt,
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1.0e-15,
        pressure_tolerance_pa: 1.0e-12,
        species_mass_tolerance_kg: 1.0e-20,
        stream_energy_tolerance_j: 1.0e-18,
        ..accuracy(dt, 1.0e-6)
    };
    assert!(matches!(
        advance_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &links,
            dt,
            &endpoints(&gases),
            controls(),
            strict,
            || false,
        ),
        Err(FiniteChannelAdvanceError::MinimumStep { .. }
            | FiniteChannelAdvanceError::Step { .. })
    ));
    assert_eq!(gases, original_gas);
    assert_eq!(walls[&9].temperature_k(), original_temperature);
}
