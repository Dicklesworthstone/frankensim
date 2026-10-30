use super::*;
use crate::compressible::mixture::molecular::{DissipationControl, MolecularTransport2d};
use crate::compressible::mixture::tests::molecular::{
    impulses as molecular_impulses, law, prepare,
};
use fs_gas::multicomponent::dissipation::ConstantTransport;

fn budget(flow: &ReactingEuler2d, operator: &MolecularTransport2d<ConstantTransport>) -> usize {
    flow.workspace_bytes() + operator.workspace_bytes()
}

#[test]
fn g3_uniform_molecular_reactions_recover_original_chemistry_and_shared_cfl_retries() {
    let mechanism = make_mechanism(40.0, 500.0, 3.5);
    let mut reference = uniform(&mechanism);
    let expected = reference
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    let mut flow = uniform(&mechanism);
    let mut operator = prepare(flow.flow(), law(flow.flow().gas(), 0.02, 0.03, 0.01, 0.03));
    let bytes = budget(&flow, &operator);
    let report = flow
        .step_with_dissipation(
            &mut operator,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    coupled_balance(&report);
    assert!(report.coupling_retries > 0);
    assert_eq!(report.elapsed, expected.elapsed);
    assert_eq!(report.transport, expected.transport);
    assert_eq!(report.chemistry, expected.chemistry);
    assert_eq!(flow.flow().states(), reference.flow().states());
    assert_eq!(flow.flow().time(), reference.flow().time());
    let split = report.dissipation.unwrap();
    assert_eq!(split.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(split.second.elapsed, split.first.elapsed);
    assert_eq!(operator.last_step(), Some(&split.second));
    assert_eq!(budget(&flow, &operator), bytes);
}

#[test]
fn g4_reacting_molecular_wall_loads_and_both_histories_publish_as_one_interval() {
    let mechanism = make_mechanism(5.0, 8.0, 3.5);
    let create = || {
        let flow = wall_field(&mechanism);
        let operator = prepare(flow.flow(), law(flow.flow().gas(), 0.02, 0.03, 0.01, 0.02));
        (flow, operator)
    };
    let (mut reference, mut reference_operator) = create();
    let bytes = budget(&reference, &reference_operator);
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
    let split = expected.dissipation.unwrap();
    let mut outward = [0.0; 2];
    for (species, momentum, energy) in impulses(&reference) {
        assert_eq!(species, [0.0; 3]);
        assert_eq!(energy, 0.0);
        for axis in 0..2 {
            outward[axis] += momentum[axis];
        }
    }
    for axis in 0..2 {
        close(
            expected.after.momentum[axis] - expected.before.momentum[axis] + outward[axis],
            0.0,
            3e-12,
        );
        close(
            expected.boundary_change.momentum[axis],
            expected.transport.boundary_change.momentum[axis]
                + split.first.boundary_change.momentum[axis]
                + split.second.boundary_change.momentum[axis],
            3e-12,
        );
    }
    assert_ne!(
        expected.boundary_change.momentum, expected.transport.boundary_change.momentum,
        "molecular wall traction must reach the complete support load"
    );

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
    let clock = flow.flow().time();
    let report = *flow.last_step().unwrap();
    let advection_report = *flow.flow().last_step().unwrap();
    let molecular_report = *operator.last_step().unwrap();
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
    assert_eq!(flow.flow().time(), clock);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.flow().last_step(), Some(&advection_report));
    assert_eq!(operator.last_step(), Some(&molecular_report));
    assert_eq!(impulses(&flow), loads);
    assert_eq!(molecular_impulses(&operator), molecular_loads);
    assert_eq!(budget(&flow, &operator), bytes);

    let refused = flow.step_with_dissipation(
        &mut operator,
        Time::new(0.001),
        ReactiveStepControl {
            maximum_coupling_retries: 0,
            ..ReactiveStepControl::default()
        },
        DissipationControl {
            maximum_step_s: 1e-6,
            maximum_steps: 1,
            ..DissipationControl::default()
        },
        bytes,
        || false,
    );
    assert!(matches!(refused, Err(ReactiveError::RetryLimit { .. })));
    assert_eq!(flow.flow().states(), cells);
    assert_eq!(operator.last_step(), Some(&molecular_report));
    assert_eq!(impulses(&flow), loads);
    assert_eq!(molecular_impulses(&operator), molecular_loads);
    assert!(matches!(
        flow.step_with_dissipation(
            &mut operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            molecular,
            bytes - 1,
            || false
        ),
        Err(ReactiveError::Invalid(_))
    ));

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
