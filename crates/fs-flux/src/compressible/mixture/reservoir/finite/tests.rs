use super::*;
use crate::compressible::mixture::tests::{close, standard};

fn held(side: BoundarySide, pressure: f64, shear: f64) -> MixtureEuler2d {
    let gas = standard(0.4);
    let grid = CartesianGrid::new(2, 2, 0.2, 0.3, 0.4).unwrap();
    let axis = side.axis();
    let mut pairs = [AxisBoundaries::Periodic; 2];
    pairs[axis] = AxisBoundaries::Ends {
        lower: if side.sign() < 0.0 { BoundaryCondition::Transmissive } else { BoundaryCondition::Reflecting },
        upper: if side.sign() > 0.0 { BoundaryCondition::Transmissive } else { BoundaryCondition::Reflecting },
    };
    let mut velocity = [0.0; 2];
    velocity[1 - axis] = shear;
    let mut row = [0.0; 6];
    state_into(&gas, &[0.7, 0.3, 0.0], velocity, 1.0, &mut row).unwrap();
    let mut states = Vec::with_capacity(6 * grid.len());
    for _ in 0..grid.len() { states.extend_from_slice(&row); }
    let bytes = MixtureEuler2d::required_workspace_bytes_with_reservoirs(&gas, &grid, 1).unwrap();
    MixtureEuler2d::new_with_reservoirs(
        gas, grid, states,
        Boundaries2d { x: pairs[0], y: pairs[1] },
        &[ReservoirBoundary {
            side,
            partial_densities: &[0.15, 0.85, 0.0],
            velocity: [0.0; 2],
            pressure_pa: pressure,
        }],
        bytes, || false,
    ).unwrap()
}

fn finite(side: BoundarySide, pressure: f64, volume_m3: f64, shear: f64) -> MixtureEuler2d {
    let flow = held(side, pressure, shear);
    let bytes = flow.required_workspace_bytes_with_finite_reservoirs(1).unwrap();
    let flow = flow.with_finite_reservoirs(&[FiniteReservoir {side, volume_m3}], bytes, || false).unwrap();
    assert_eq!(flow.workspace_bytes(), bytes);
    flow
}

fn storage(flow: &MixtureEuler2d, side: BoundarySide) -> [f64; 4] {
    let total = flow.totals();
    let tank = flow.finite_reservoir(side).unwrap();
    let q = tank.state();
    [
        total.species_masses[0] + q.partial_densities[0] * tank.volume_m3(),
        total.species_masses[1] + q.partial_densities[1] * tank.volume_m3(),
        total.species_masses[2] + q.partial_densities[2] * tank.volume_m3(),
        total.energy + q.energy * tank.volume_m3(),
    ]
}

#[test]
fn g0_planar_finite_mixture_tanks_exchange_species_energy_and_support_all_sides() {
    for side in SIDES {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            for pressure in [0.5, 2.0] {
                let shear = if pressure < 1.0 { -0.3 } else { 0.0 };
                let mut flow = finite(side, pressure, 0.04, shear);
                let original = storage(&flow, side);
                let initial_tank = flow.finite_reservoir(side).unwrap().clone();
                for _ in 0..3 {
                    let report = flow.step(Time::new(0.0002), StepControl {
                        reconstruction, ..StepControl::default()
                    }, || false).unwrap();
                    assert_eq!(report.halvings, 0);
                    for (got, want) in storage(&flow, side).into_iter().zip(original) {
                        close(got, want, 4e-12);
                    }
                    let tank = flow.finite_reservoir_report(side).unwrap();
                    for d in 0..3 {
                        close(tank.species_residual_kg[d], 0.0, 2e-14);
                        close(report.boundary_change.species_masses[d],
                            -tank.transfer.species_masses[d], 2e-14);
                        close(tank.after.species_masses[d],
                            flow.finite_reservoir(side).unwrap().state().partial_densities[d] * 0.04, 2e-14);
                    }
                    close(tank.energy_residual_j, 0.0, 2e-14);
                    close(tank.absolute_energy_residual_j, 0.0, 2e-12);
                    close(report.boundary_change.energy, -tank.transfer.energy, 2e-14);
                    assert!(if pressure < 1.0 {tank.transfer.mass() > 0.0} else {tank.transfer.mass() < 0.0});
                    for d in 0..2 {
                        close(tank.transfer.momentum[d],
                            flow.boundary_impulses(side).map(|q| q.momentum[d]).sum(), 2e-14);
                    }
                    if pressure < 1.0 {
                        close(tank.transfer.momentum[1 - side.axis()], shear * tank.transfer.mass(), 3e-12);
                        assert!(tank.after.species_masses[0] / tank.after.mass() > 0.15);
                    } else {
                        close(tank.after.species_masses[0] / tank.after.mass(), 0.15, 2e-13);
                    }
                    let current = flow.reservoir(side).unwrap();
                    close(current.pressure_pa, tank.after.pressure_pa, 2e-14);
                    assert_eq!(current.partial_densities[2], 0.0);
                    assert_eq!(current.velocity, [0.0; 2]);
                    for k in 0..flow.cell_count() { assert_eq!(flow.cell(k).unwrap().partial_densities[2], 0.0); }
                }
                let before = flow.gas().primitive(initial_tank.state(), || false).unwrap();
                let after = flow.gas().primitive(flow.finite_reservoir(side).unwrap().state(), || false).unwrap();
                if pressure > 1.0 {
                    assert!(after.pressure < before.pressure);
                    assert!(after.temperature < before.temperature);
                } else {
                    assert!(after.pressure > before.pressure);
                }
            }
        }
    }
}

#[test]
fn g1_planar_finite_mixture_pressure_feedback_refines_at_rk2_order() {
    let run = |steps: usize| {
        let side = BoundarySide::XMin;
        let mut flow = finite(side, 2.0, 0.006, 0.0);
        for _ in 0..steps {
            let report = flow.step(Time::new(0.01 / steps as f64), StepControl {
                reconstruction: Reconstruction::FirstOrder, ..StepControl::default()
            }, || false).unwrap();
            assert_eq!(report.halvings, 0);
            close(report.elapsed.value(), 0.01 / steps as f64, 1e-15);
        }
        flow.reservoir(side).unwrap().pressure_pa
    };
    let reference = run(128);
    let coarse = (run(4) - reference).abs();
    let fine = (run(8) - reference).abs();
    assert!(fine > 1e-12 && coarse > 3.2 * fine, "{coarse:e} {fine:e}");
}

#[test]
fn g0_planar_finite_mixture_tiny_tank_limits_initial_timestep() {
    for side in SIDES {
        let mut flow = finite(side, 2.0, 1e-10, 0.0);
        let tank = flow.finite_reservoir(side).unwrap();
        let mut flux = [0.0; 6];
        let receipt = flow.gas().reservoir_normal_flux(
            flow.cell(0).unwrap().normal(side.axis()), tank.state(),
            side.sign(), &mut flux, || false,
        ).unwrap();
        let rate = 0.5 * flow.grid().face_count(side) as f64
            * flow.grid().face_area(side) * receipt.signal_speed / tank.volume_m3();
        let report = flow.step(Time::new(0.01), StepControl {
            cfl: 0.05, reconstruction: Reconstruction::FirstOrder, ..StepControl::default()
        }, || false).unwrap();
        assert!(report.elapsed.value() <= 0.05 / rate * (1.0 + 1e-12));
        assert!(report.elapsed.value() < 1e-8);
        assert!(report.maximum_cfl <= 0.05 * (1.0 + 1e-12));
        assert_eq!(report.halvings, 0);
    }
}

#[test]
fn g0_planar_finite_mixture_final_poll_preserves_all_accepted_state_and_retry() {
    let side = BoundarySide::YMax;
    let make = || {
        let mut flow = finite(side, 0.5, 0.02, 0.2);
        flow.step(Time::new(1e-4), StepControl::default(), || false).unwrap();
        flow
    };
    let mut probe = make();
    let mut polls = 0;
    probe.step(Time::new(1e-4), StepControl::default(), || { polls += 1; false }).unwrap();
    let mut flow = make();
    let states = flow.states().to_vec();
    let tank = flow.finite_reservoir(side).unwrap().clone();
    let report = *flow.finite_reservoir_report(side).unwrap();
    let step = flow.last_step().copied();
    let time = flow.time();
    let boundary = flow.boundary_impulses.clone();
    let bytes = flow.workspace_bytes();
    let mut seen = 0;
    assert_eq!(flow.step(Time::new(1e-4), StepControl::default(), || {
        seen += 1; seen == polls
    }), Err(Error::Cancelled));
    assert_eq!(seen, polls);
    assert_eq!(flow.states(), states);
    assert_eq!(flow.finite_reservoir(side), Some(&tank));
    assert_eq!(flow.finite_reservoir_report(side), Some(&report));
    assert_eq!(flow.last_step().copied(), step);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.boundary_impulses, boundary);
    assert_eq!(flow.workspace_bytes(), bytes);
    flow.step(Time::new(1e-4), StepControl::default(), || false).unwrap();
    assert_eq!(flow.states(), probe.states());
    assert_eq!(flow.finite_reservoir(side), probe.finite_reservoir(side));
    assert_eq!(flow.finite_reservoir_report(side), probe.finite_reservoir_report(side));
}

#[test]
fn g0_planar_finite_mixture_binding_budget_clone_and_compact_law_are_retained() {
    let side = BoundarySide::XMin;
    let input = FiniteReservoir { side, volume_m3: 0.02 };
    let flow = held(side, 2.0, 0.0);
    let bytes = flow.required_workspace_bytes_with_finite_reservoirs(1).unwrap();
    assert!(flow.with_finite_reservoirs(&[input], bytes - 1, || false).is_err());
    assert!(held(side, 2.0, 0.0).with_finite_reservoirs(
        &[input, input], MAX_MIXTURE_WORKSPACE_BYTES, || false).is_err());
    let mut flow = finite(side, 2.0, 0.02, 0.0)
        .with_stagnation_reservoirs(&[], || false).unwrap();
    assert!(flow.uses_stagnation_reservoir(side));
    flow.step(Time::new(1e-4), StepControl::default(), || false).unwrap();
    let mut clone = flow.reservoirs.try_clone(flow.gas(), &mut || false).unwrap();
    assert_eq!(clone.finite.inventory(0), flow.finite_reservoir(side));
    assert_eq!(clone.finite.report(0), flow.finite_reservoir_report(side));
    assert_eq!(clone.workspace_bytes(), flow.reservoirs.workspace_bytes());
    let initial = finite(side, 0.5, 0.04, 0.0);
    clone.copy_finite_from(initial.gas(), &initial.reservoirs, &mut || false).unwrap();
    assert_eq!(clone.finite.inventory(0), initial.finite_reservoir(side));
    assert_eq!(clone.finite.report(0), None);
    let bytes = clone.workspace_bytes();
    let saved = clone.finite.inventory(0).unwrap().clone();
    let mut polls = 0;
    clone.copy_finite_from(flow.gas(), &flow.reservoirs, &mut || {polls += 1; false}).unwrap();
    let mut check = initial.reservoirs.try_clone(initial.gas(), &mut || false).unwrap();
    let mut seen = 0;
    assert_eq!(check.copy_finite_from(flow.gas(), &flow.reservoirs, &mut || {
        seen += 1; seen == polls
    }), Err(Error::Cancelled));
    assert_eq!(check.finite.inventory(0), Some(&saved));
    assert_eq!(check.workspace_bytes(), bytes);
}
