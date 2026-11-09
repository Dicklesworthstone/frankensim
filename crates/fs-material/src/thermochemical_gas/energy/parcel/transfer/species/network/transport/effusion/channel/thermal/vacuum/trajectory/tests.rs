//! G1/G3/G4: analytic blowdown, nonlinear wall feedback and atomic refusal.
use super::super::super::super::super::tests::{controls, gas, interface, near};
use super::super::super::tests::{channel, endpoints, material};
use super::*;
use crate::thermal::transport::volumetric::timeline::FixedVolumeThermalLaw;
use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::channel::{
    CylindricalMolecularGeometry, DiffuseMolecularChannel,
};

fn link(channel: &DiffuseMolecularChannel) -> FiniteVacuumChannelLink<'_> {
    FiniteVacuumChannelLink {
        id: 7,
        parcel: 0,
        channel,
        area_m2: 1e-6,
        wall: 9,
    }
}
fn accuracy(dt: f64, relative: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: dt,
        minimum_step_s: dt * 1e-6,
        maximum_step_s: dt,
        maximum_attempts: 4096,
        maximum_accepted_intervals: 1024,
        relative_tolerance: relative,
        temperature_tolerance_k: 1e-6,
        pressure_tolerance_pa: 1e-4,
        species_mass_tolerance_kg: 1e-14,
        stream_energy_tolerance_j: 1e-5,
    }
}

#[test]
fn ballistic_limit_converges_to_independent_monatomic_blowdown() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 2e-8, 900.0).unwrap(),
    )]);
    let path = channel(0.0);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let dt = 1.0;
    // dN/dt=-lambda*N*sqrt(T/T0), T/T0=(N/N0)^(1/3), hence
    // N/N0=(1+lambda*t/6)^(-6), T/T0=(1+lambda*t/6)^(-2).
    let lambda =
        1e-6 / gases[&0].volume_m3() * (R * 400.0 / (2.0 * core::f64::consts::PI * 0.028)).sqrt();
    let factor = 1.0 + lambda * dt / 6.0;
    let expected_t = 400.0 / factor.powi(2);
    let expected_m = gases[&0].mass_kg() / factor.powi(6);
    let mut previous_error = f64::INFINITY;
    for relative in [1e-2, 1e-3] {
        let candidate = advance_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            dt,
            &endpoints(&gases),
            controls(),
            accuracy(dt, relative),
            || false,
        )
        .unwrap();
        let error = (candidate.parcels()[&0].state().state().temperature - expected_t).abs();
        assert!(error < previous_error);
        previous_error = error;
        assert_eq!(
            candidate.walls()[&9].temperature_k(),
            walls[&9].temperature_k()
        );
        assert_eq!(candidate.step().balance.wall_energy_change_j, 0.0);
        assert_eq!(
            candidate.step().accepted_steps.len(),
            2 * candidate.step().accepted_intervals
        );
        assert!(candidate.step().maximum_accepted_error <= 1.0);
        assert!(
            candidate.step().balance.energy_residual_j.abs()
                <= controls().exchange.energy_tolerance_j
        );
        if relative == 1e-3 {
            near(
                candidate.parcels()[&0].state().state().temperature,
                expected_t,
                0.1,
            );
            near(candidate.parcels()[&0].mass_kg(), expected_m, 1e-7);
        }
    }
}

#[test]
fn adaptive_hot_wall_matches_small_native_steps_and_whole_window_storage() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 2e-8, 900.0).unwrap(),
    )]);
    let path = channel(2e-8);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let dt = 0.02;
    let links = [link(&path)];
    let solve = endpoints(&gases);
    let candidate = advance_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &links,
        dt,
        &solve,
        controls(),
        accuracy(dt, 1e-3),
        || false,
    )
    .unwrap();
    let (mut reference_gas, mut reference_wall) = (gases.clone(), walls.clone());
    for _ in 0..512 {
        let next = step_candidate(
            &reference_gas,
            &reference_wall,
            GasTransportPaths::default(),
            &links,
            dt / 512.0,
            &solve,
            controls(),
            || false,
        )
        .unwrap();
        (reference_gas, reference_wall, _) = next.into_parts();
    }
    near(
        candidate.parcels()[&0].state().state().temperature,
        reference_gas[&0].state().state().temperature,
        0.04,
    );
    near(
        candidate.walls()[&9].temperature_k(),
        reference_wall[&9].temperature_k(),
        0.04,
    );
    assert!(candidate.walls()[&9].temperature_k() < 900.0);
    assert!(candidate.step().accepted_intervals > 1);
    let actual_wall = walls[&9].volume_m3()
        * source
            .storage_change_j_m3(900.0, candidate.walls()[&9].temperature_k())
            .unwrap();
    near(
        candidate.step().balance.wall_energy_change_j,
        actual_wall,
        1e-10,
    );
    assert!(candidate.step().balance.wall_entropy_change_j_k < 0.0);
    let gas_change = 1.5 * R / 0.028
        * (candidate.parcels()[&0].mass_kg() * candidate.parcels()[&0].state().state().temperature
            - gases[&0].mass_kg() * 400.0);
    near(
        gas_change + actual_wall + candidate.step().balance.gas.escaped_energy_j,
        0.0,
        1e-8,
    );
    near(
        candidate.step().balance.gas.escaped_mass_kg,
        gases[&0].mass_kg() - candidate.parcels()[&0].mass_kg(),
        1e-15,
    );
    assert_eq!(gases[&0].state().state().temperature, 400.0);
    assert_eq!(walls[&9].temperature_k(), 900.0);
}

#[test]
fn formation_and_wall_storage_datums_do_not_loosen_adaptive_accuracy() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let mut reference = None;
    for (formation, wall_datum) in [(0.0, 200.0), (1e20, 200.0), (0.0, 800.0)] {
        let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, wall_datum)
            .unwrap()
            .with_sensible_entropy();
        let walls = BTreeMap::from([(
            9,
            FixedVolumeChannelWall::try_new_isochoric(&law, 2e-8, 900.0).unwrap(),
        )]);
        let path = DiffuseMolecularChannel::try_new(
            interface(formation),
            CylindricalMolecularGeometry::try_new(1e-8, 2e-8, 1e-4, "synthetic fixed channel")
                .unwrap(),
        )
        .unwrap();
        let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, formation))]);
        let dt = 0.02;
        let candidate = advance_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            dt,
            &endpoints(&gases),
            controls(),
            accuracy(dt, 1e-3),
            || false,
        )
        .unwrap();
        let values = [
            candidate.parcels()[&0].state().state().temperature,
            candidate.walls()[&9].temperature_k(),
            candidate.step().balance.wall_energy_change_j,
        ];
        if let Some((original, count)) = reference {
            for (value, target) in values.into_iter().zip(original) {
                near(value, target, 1e-7);
            }
            assert_eq!(candidate.step().accepted_intervals, count);
        } else {
            reference = Some((values, candidate.step().accepted_intervals));
        }
    }
}

#[test]
fn work_minimum_step_and_late_cancellation_publish_neither_inventory() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 2e-8, 900.0).unwrap(),
    )]);
    let path = channel(2e-8);
    let gases = BTreeMap::from([(0, gas(400.0, 100_000.0, 0.0))]);
    let original = gases.clone();
    let dt = 0.02;
    let strict = EffusionAdvanceControls {
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1e-12,
        pressure_tolerance_pa: 1e-10,
        species_mass_tolerance_kg: 1e-20,
        stream_energy_tolerance_j: 1e-15,
        minimum_step_s: dt,
        ..accuracy(dt, 1e-3)
    };
    assert!(matches!(
        advance_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            dt,
            &endpoints(&gases),
            controls(),
            strict,
            || false
        ),
        Err(Error::MinimumStep { .. })
    ));
    let bounded = EffusionAdvanceControls {
        maximum_attempts: 3,
        minimum_step_s: dt * 1e-6,
        ..strict
    };
    assert!(matches!(
        advance_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            dt,
            &endpoints(&gases),
            controls(),
            bounded,
            || false
        ),
        Err(Error::WorkBudget { .. })
    ));
    let mut polls = 0usize;
    let completed = advance_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        dt,
        &endpoints(&gases),
        controls(),
        accuracy(dt, 1e-3),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(completed.step().accepted_intervals > 0);
    let mut seen = 0usize;
    assert!(matches!(
        advance_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            dt,
            &endpoints(&gases),
            controls(),
            accuracy(dt, 1e-3),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(gases, original);
    assert_eq!(walls[&9].temperature_k(), 900.0);
}
