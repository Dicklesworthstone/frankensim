use super::*;
use crate::compressible::spatial::mixture::molecular::{MolecularWall3d, WallThermalCondition};

fn navier_reactor(fraction: f64) -> (ReactingEuler3d, MolecularTransport3d<ConstantTransport>) {
    let (gas, mechanism) = mechanism(2.0, 3.0, 0.0);
    let flow = coupled(
        field(
            gas,
            CartesianGrid::new(3, 3, 4, 0.3, 0.3, 0.25).unwrap(),
            Boundaries3d {
                x: AxisBoundaries::Periodic,
                y: AxisBoundaries::Periodic,
                z: WALLS,
            },
            |_| ([0.8, 0.2], [0.3, -0.1, 0.0], 2.0),
        ),
        mechanism,
    );
    let mut molecular = operator(flow.flow(), 0.01);
    for side in [BoundarySide::ZMin, BoundarySide::ZMax] {
        molecular = molecular
            .with_wall(
                side,
                MolecularWall3d::NavierSlip {
                    velocity: [0.0; 3],
                    thermal: WallThermalCondition::Adiabatic,
                    slip_length_m: 0.2,
                    fluid_heat_fraction: fraction,
                },
            )
            .unwrap();
    }
    (flow, molecular)
}

#[test]
fn g3_reactive_navier_walls_keep_chemistry_species_and_interface_heat_in_one_budget() {
    for fraction in [0.0, 1.0] {
        let (mut flow, mut molecular) = navier_reactor(fraction);
        let before = flow.flow().totals();
        let absolute_before = before.absolute_energy(flow.flow().gas()).unwrap();
        let bytes = budget(&flow, &molecular);
        let report = flow
            .step_with_dissipation(
                &mut molecular,
                Time::new(0.004),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || false,
            )
            .unwrap();
        balance(&report);
        assert!(report.after.species_masses[0] < before.species_masses[0]);
        assert!(report.after.species_masses[1] > before.species_masses[1]);
        let mut outward_energy = 0.0;
        for side in [BoundarySide::ZMin, BoundarySide::ZMax] {
            for impulse in flow.boundary_impulses(side) {
                assert_eq!(impulse.species_masses, [0.0; 2]);
                outward_energy += impulse.energy;
            }
        }
        let absolute_after = report.after.absolute_energy(flow.flow().gas()).unwrap();
        close(
            absolute_after - absolute_before + outward_energy,
            0.0,
            3e-10,
        );
        if fraction == 1.0 {
            assert_eq!(outward_energy, 0.0);
        } else {
            assert!(outward_energy > 0.0);
            assert!(absolute_after < absolute_before);
        }
        let split = report.dissipation.unwrap();
        assert_eq!(split.first.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(split.second.elapsed, split.first.elapsed);
        assert_eq!(molecular.last_step(), Some(&split.second));
    }
}

#[test]
fn g4_reactive_navier_cancellation_preserves_composition_wall_heat_and_replay() {
    let (mut reference, mut reference_molecular) = navier_reactor(0.3);
    let bytes = budget(&reference, &reference_molecular);
    let mut polls = 0;
    let accepted = reference
        .step_with_dissipation(
            &mut reference_molecular,
            Time::new(0.004),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut flow, mut molecular) = navier_reactor(0.3);
    let initial = flow.flow().states().to_vec();
    let initial_loads = loads(&flow);
    let initial_molecular = molecular_loads(&molecular);
    let mut calls = 0;
    assert!(
        flow.step_with_dissipation(
            &mut molecular,
            Time::new(0.004),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || {
                calls += 1;
                calls == polls - 1
            }
        )
        .is_err()
    );
    assert!(calls > polls / 2);
    assert_eq!(flow.flow().states(), initial);
    assert_eq!(flow.flow().time().value(), 0.0);
    assert_eq!(loads(&flow), initial_loads);
    assert_eq!(molecular_loads(&molecular), initial_molecular);
    assert!(molecular.last_step().is_none());
    let replay = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.004),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    assert_eq!(replay, accepted);
    assert_eq!(flow.flow().states(), reference.flow().states());
    assert_eq!(loads(&flow), loads(&reference));
    assert_eq!(
        molecular_loads(&molecular),
        molecular_loads(&reference_molecular)
    );
}
