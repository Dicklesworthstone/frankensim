//! Native gas changes retain resolved no-slip momentum and separate its work.
//! Synthetic material changes exercise the mechanism, not measured gas fidelity.

use fs_couple::bernoulli_aperture::slit_channel::{
    FixedChannel, SlitChannel, SlitChannelSpec, SlitExit,
    circular::{CircularChannel, CircularChannelSpec},
};
use fs_material::gas::{GasSpec, GasState};

const LENGTH: f64 = 0.012;
const WIDTH: f64 = 0.003;
const HEIGHT: f64 = 0.0004;
const RADIUS: f64 = 0.0003;
const DT: f64 = 1.0 / 48_000.0;

fn gas(viscosity_scale: f64, pressure: f64) -> GasState {
    let mut spec = GasSpec::dry_air_ussa1976();
    spec.sutherland_beta *= viscosity_scale;
    GasState::try_new(&spec, 293.15, pressure).unwrap()
}

fn channel(circular: bool, cells: usize, gas: &GasState, exit: SlitExit) -> FixedChannel {
    if circular {
        CircularChannel::from_gas(
            gas,
            CircularChannelSpec {
                length_m: LENGTH,
                radius_m: RADIUS,
                radial_cells: cells,
                sample_period_s: DT,
                exit,
            },
        )
        .unwrap()
        .into()
    } else {
        SlitChannel::from_gas(
            gas,
            SlitChannelSpec {
                length_m: LENGTH,
                width_m: WIDTH,
                height_m: HEIGHT,
                transverse_cells: cells,
                sample_period_s: DT,
                exit,
            },
        )
        .unwrap()
        .into()
    }
}

fn cell_volume(circular: bool, cells: usize, index: usize) -> f64 {
    if circular {
        let inner = RADIUS * index as f64 / cells as f64;
        let outer = RADIUS * (index + 1) as f64 / cells as f64;
        core::f64::consts::PI * LENGTH * (outer * outer - inner * inner)
    } else {
        LENGTH * WIDTH * HEIGHT / cells as f64
    }
}

fn independent_energy(circular: bool, density: f64, velocities: &[f64]) -> f64 {
    velocities
        .iter()
        .enumerate()
        .map(|(i, v)| 0.5 * density * cell_volume(circular, velocities.len(), i) * v * v)
        .sum()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale.max(1e-30),
        "actual {actual:e}, expected {expected:e}, scale {scale:e}"
    );
}

#[test]
fn g1_source_change_preserves_each_cell_momentum_and_accounts_actual_kinetic_work() {
    let before = gas(1.0, 101_325.0);
    let after = gas(1.7, 157_000.0);
    for circular in [false, true] {
        for exit in [SlitExit::StaticPressure, SlitExit::UniformKineticHead] {
            let mut owner = channel(circular, 17, &before, exit);
            for pressure in [90.0, 90.0, -20.0, 40.0, -60.0] {
                owner.step(pressure).unwrap();
            }
            let velocities = owner.velocities_m_s().to_vec();
            let last = owner.last_step();
            let old_energy = independent_energy(circular, before.density, &velocities);
            assert!(old_energy > 0.0);
            let report = owner
                .rebind_gas_controlled(&after, 1.0, &mut || false)
                .unwrap();
            assert_eq!(owner.gas(), &after);
            assert_eq!(owner.last_step(), last, "accepted exchange is historical");
            assert_eq!(owner.sample_period_s(), DT);
            assert_eq!(owner.velocities_m_s().len(), velocities.len());
            for (i, (&old, &new)) in velocities.iter().zip(owner.velocities_m_s()).enumerate() {
                let volume = cell_volume(circular, velocities.len(), i);
                let momentum = before.density * volume * old;
                near(after.density * volume * new, momentum, momentum.abs());
            }
            let new_energy = independent_energy(circular, after.density, owner.velocities_m_s());
            near(report.kinetic_energy_before_j, old_energy, old_energy);
            near(report.kinetic_energy_after_j, new_energy, new_energy);
            near(
                new_energy,
                before.density / after.density * old_energy,
                old_energy,
            );
            near(
                report.parameter_work_j,
                new_energy - old_energy,
                old_energy + new_energy,
            );
            assert_eq!(owner.stored_energy_j(), report.kinetic_energy_after_j);
            assert_ne!(owner.stored_energy_j(), last.kinetic_energy_j);
            // The next actual step starts from the transferred energy. No part
            // of the zero-duration work is relabeled as wall or numerical loss.
            let step = owner.step(0.0).unwrap();
            assert_eq!(step.pressure_work_j, 0.0);
            assert_eq!(owner.stored_energy_j(), step.kinetic_energy_j);
            near(
                step.stored_delta_j,
                step.kinetic_energy_j - new_energy,
                new_energy,
            );
            near(
                new_energy - step.kinetic_energy_j,
                step.viscous_dissipation_j + step.numerical_dissipation_j + step.exit_work_j,
                new_energy,
            );
        }
    }
}

#[test]
fn g3_viscosity_only_refresh_changes_decay_without_instantaneous_work_or_reset() {
    let before = gas(1.0, 101_325.0);
    let viscous = gas(3.0, 101_325.0);
    assert_eq!(before.density, viscous.density);
    for circular in [false, true] {
        let mut changed = channel(circular, 32, &before, SlitExit::StaticPressure);
        for _ in 0..12 {
            changed.step(70.0).unwrap();
        }
        let mut control = changed.clone();
        let previous = changed.last_step();
        let report = changed
            .rebind_gas_controlled(&viscous, 0.0, &mut || false)
            .unwrap();
        assert_eq!(report.parameter_work_j, 0.0);
        assert_eq!(changed.velocities_m_s(), control.velocities_m_s());
        assert_eq!(changed.stored_energy_j(), control.stored_energy_j());
        assert_eq!(changed.last_step(), previous);
        assert_ne!(
            changed.pressure_load().unwrap(),
            control.pressure_load().unwrap()
        );
        for _ in 0..8 {
            changed.step(0.0).unwrap();
            control.step(0.0).unwrap();
        }
        assert!(changed.stored_energy_j() < 0.99 * control.stored_energy_j());
        // Reapplying an identical complete source preserves the exact operator,
        // velocity, endpoint energy and historical report, including on retry.
        let mut unchanged = changed.clone();
        assert_eq!(
            changed
                .rebind_gas_controlled(&viscous, 0.0, &mut || false)
                .unwrap()
                .parameter_work_j,
            0.0
        );
        for pressure in [10.0, -30.0, 0.0] {
            assert_eq!(
                changed.step(pressure).unwrap(),
                unchanged.step(pressure).unwrap()
            );
            assert_eq!(changed.velocities_m_s(), unchanged.velocities_m_s());
        }
    }
}

#[test]
fn g4_rejected_gas_source_work_and_final_cancellation_preserve_exact_retry() {
    let before = gas(1.0, 101_325.0);
    let after = gas(2.0, 80_000.0);
    for circular in [false, true] {
        let mut owner = channel(circular, 130, &before, SlitExit::StaticPressure);
        for _ in 0..10 {
            owner.step(50.0).unwrap();
        }
        let mut clean = owner.clone();
        let velocities = owner.velocities_m_s().to_vec();
        let last = owner.last_step();
        let energy = owner.stored_energy_j();
        let load = owner.pressure_load().unwrap();
        let same = |owner: &FixedChannel| {
            assert_eq!(owner.gas(), &before);
            assert_eq!(owner.velocities_m_s(), velocities);
            assert_eq!(owner.last_step(), last);
            assert_eq!(owner.stored_energy_j(), energy);
            assert_eq!(owner.pressure_load().unwrap(), load);
        };
        let mut bad = after;
        bad.density *= 1.1;
        assert!(
            owner
                .rebind_gas_controlled(&bad, 1.0, &mut || false)
                .is_err()
        );
        same(&owner);
        for budget in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                owner
                    .rebind_gas_controlled(&after, budget, &mut || false)
                    .is_err()
            );
            same(&owner);
        }
        let mut polls = 0;
        let accepted = clean
            .rebind_gas_controlled(&after, 1.0, &mut || {
                polls += 1;
                false
            })
            .unwrap();
        assert!(polls >= 5 && accepted.parameter_work_j > 0.0);
        let mut calls = 0;
        assert!(
            owner
                .rebind_gas_controlled(&after, 1.0, &mut || {
                    calls += 1;
                    calls == polls
                })
                .is_err()
        );
        assert_eq!(calls, polls);
        same(&owner);
        assert_eq!(
            owner
                .rebind_gas_controlled(&after, 1.0, &mut || false)
                .unwrap(),
            accepted
        );
        for pressure in [-50.0, 0.0, 70.0] {
            assert_eq!(owner.step(pressure).unwrap(), clean.step(pressure).unwrap());
            assert_eq!(owner.velocities_m_s(), clean.velocities_m_s());
        }
    }
}
