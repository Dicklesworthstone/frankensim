use super::*;
use crate::compressible::mixture::molecular::DissipationControl;
use crate::compressible::mixture::tests::molecular::{
    impulses as molecular_impulses, law, prepare,
};
use crate::compressible::mixture::tests::reservoir::bind;

#[test]
fn g3_reacting_reservoir_exchange_and_source_metadata_publish_at_one_interval() {
    let mechanism = make_mechanism(5.0, 8.0, 3.5);
    let grid = CartesianGrid::new(3, 2, 1.0 / 3.0, 0.5, 0.6).unwrap();
    let boundaries = Boundaries2d {
        x: OPEN,
        y: AxisBoundaries::Periodic,
    };
    let baths = [
        ReservoirBoundary {
            side: BoundarySide::XMin,
            partial_densities: &[0.8, 0.2, 0.0],
            velocity: [0.2, 0.1],
            pressure_pa: 1.2,
        },
        ReservoirBoundary {
            side: BoundarySide::XMax,
            partial_densities: &[0.2, 0.8, 0.0],
            velocity: [0.1, -0.1],
            pressure_pa: 0.9,
        },
    ];
    let create = || {
        let flow = bind(
            raw_field(&mechanism, grid, boundaries, |x| {
                (
                    [0.5 + 0.1 * x[1], 0.5 - 0.1 * x[1], 0.0],
                    [0.1, 0.05 * x[0]],
                    1.0,
                )
            }),
            &baths,
        );
        let flow = reactive(flow, &mechanism);
        let operator = prepare(flow.flow(), law(flow.flow().gas(), 0.02, 0.03, 0.01, 0.02));
        for bath in baths {
            assert_eq!(flow.flow().reservoir(bath.side), Some(bath));
        }
        (flow, operator)
    };
    let (mut reference, mut reference_operator) = create();
    let bytes = reference.workspace_bytes() + reference_operator.workspace_bytes();
    let molecular = DissipationControl {
        maximum_step_s: 0.0005,
        ..DissipationControl::default()
    };
    reference
        .step_with_dissipation(
            &mut reference_operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            || false,
        )
        .unwrap();
    let mut polls = 0;
    let expected = reference
        .step_with_dissipation(
            &mut reference_operator,
            Time::new(0.002),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    coupled_balance(&expected);
    // Both accepted publications swap the private full-step trial into place;
    // its copied original-source metadata must survive those swaps exactly.
    for bath in baths {
        assert_eq!(reference.flow().reservoir(bath.side), Some(bath));
    }
    let split = expected.dissipation.unwrap();
    assert!(expected.chemistry.referenced_energy_change_j > 0.0);
    assert!(split.first.boundary_change.species_masses[0].abs() > 1e-12);
    assert!(split.second.boundary_change.energy.abs() > 1e-12);
    for k in 0..3 {
        close(
            expected.boundary_change.species_masses[k],
            expected.transport.boundary_change.species_masses[k]
                + split.first.boundary_change.species_masses[k]
                + split.second.boundary_change.species_masses[k],
            3e-12,
        );
    }
    for axis in 0..2 {
        close(
            expected.boundary_change.momentum[axis],
            expected.transport.boundary_change.momentum[axis]
                + split.first.boundary_change.momentum[axis]
                + split.second.boundary_change.momentum[axis],
            3e-12,
        );
    }
    close(
        expected.boundary_change.energy,
        expected.transport.boundary_change.energy
            + split.first.boundary_change.energy
            + split.second.boundary_change.energy,
        3e-12,
    );
    let mut outward = [0.0; 6];
    for (species, momentum, energy) in impulses(&reference) {
        for k in 0..3 {
            outward[k] += species[k];
        }
        for axis in 0..2 {
            outward[3 + axis] += momentum[axis];
        }
        outward[5] += energy;
    }
    for k in 0..3 {
        close(
            expected.boundary_change.species_masses[k] + outward[k],
            0.0,
            3e-12,
        );
    }
    for axis in 0..2 {
        close(
            expected.boundary_change.momentum[axis] + outward[3 + axis],
            0.0,
            3e-12,
        );
    }
    close(expected.boundary_change.energy + outward[5], 0.0, 3e-12);

    let (mut flow, mut operator) = create();
    flow.step_with_dissipation(
        &mut operator,
        Time::new(0.001),
        ReactiveStepControl::default(),
        molecular,
        bytes,
        || false,
    )
    .unwrap();
    let cells = flow.flow().states().to_vec();
    let time = flow.flow().time();
    let previous = *flow.last_step().unwrap();
    let advective = *flow.flow().last_step().unwrap();
    let molecular_previous = *operator.last_step().unwrap();
    let loads = impulses(&flow);
    let molecular_loads = molecular_impulses(&operator);
    let mut attempted = 0;
    let error = flow
        .step_with_dissipation(
            &mut operator,
            Time::new(0.002),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            || {
                attempted += 1;
                attempted >= polls
            },
        )
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(flow.flow().states(), cells);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(flow.last_step(), Some(&previous));
    assert_eq!(flow.flow().last_step(), Some(&advective));
    assert_eq!(operator.last_step(), Some(&molecular_previous));
    assert_eq!(impulses(&flow), loads);
    assert_eq!(molecular_impulses(&operator), molecular_loads);
    assert_eq!(flow.workspace_bytes() + operator.workspace_bytes(), bytes);
    for bath in baths {
        assert_eq!(flow.flow().reservoir(bath.side), Some(bath));
    }
    let replay = flow
        .step_with_dissipation(
            &mut operator,
            Time::new(0.002),
            ReactiveStepControl::default(),
            molecular,
            bytes,
            || false,
        )
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(flow.flow().states(), reference.flow().states());
    assert_eq!(impulses(&flow), impulses(&reference));
    assert_eq!(
        molecular_impulses(&operator),
        molecular_impulses(&reference_operator)
    );
}
