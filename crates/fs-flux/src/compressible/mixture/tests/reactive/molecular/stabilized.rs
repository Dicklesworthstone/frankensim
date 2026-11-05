//! Original chemistry and molecular sources share one stabilized interval and work cap.
use super::*;
use crate::compressible::mixture::molecular::{
    DissipationError, MolecularReport, PlanarTransport, StabilizedDissipationControl,
};
use crate::compressible::mixture::{BodyForce2d, ForceControl2d};

fn prepared(flow: &MixtureEuler2d, diffusion: f64) -> MolecularTransport2d<ConstantTransport> {
    let source = PlanarTransport::new(
        flow.gas().clone(),
        law(flow.gas(), diffusion, 0.03, 0.01, 0.03),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport2d::required_stabilized_workspace_bytes(flow, &source).unwrap();
    MolecularTransport2d::new(flow, source, bytes, || false).unwrap()
}

#[derive(Debug, PartialEq)]
struct Published {
    states: Vec<f64>,
    time: Time,
    coupled: Option<ReactiveStepReport>,
    advection: Option<crate::compressible::mixture::StepReport>,
    molecular: Option<MolecularReport>,
    boundary: Vec<([f64; 3], [f64; 2], f64)>,
    molecular_boundary: Vec<([f64; 3], [f64; 2], f64)>,
}
fn published(
    flow: &ReactingEuler2d,
    molecular: &MolecularTransport2d<ConstantTransport>,
) -> Published {
    Published {
        states: flow.flow().states().to_vec(),
        time: flow.flow().time(),
        coupled: flow.last_step().copied(),
        advection: flow.flow().last_step().copied(),
        molecular: molecular.last_step().copied(),
        boundary: impulses(flow),
        molecular_boundary: molecular_impulses(molecular),
    }
}

#[test]
fn g1_stabilized_reacting_split_preserves_native_chemistry_diffusion_and_signed_body_work() {
    let mechanism = make_mechanism(6.0, 8.0, 4.5);
    let velocity = [0.7, -0.3];
    let make = || {
        raw_field(
            &mechanism,
            CartesianGrid::new(4, 3, 0.25, 1.0 / 3.0, 0.7).unwrap(),
            PERIODIC,
            |x| {
                let a = 0.7 + 0.1 * (2.0 * std::f64::consts::PI * (x[0] + x[1])).sin();
                ([a, 1.0 - a, 0.0], velocity, 1.0)
            },
        )
    };
    let control = ReactiveStepControl::default();
    let molecular_control = StabilizedDissipationControl::default();
    let cap = ForceControl2d { maximum_mach: 8.0 };
    let acceleration = [-3.0, 1.0];
    let force = BodyForce2d::UniformAcceleration(acceleration);
    let dt = 0.002;

    // Compose the public native D, C-A-C and D owners independently. The
    // middle owner retains the original source-bound SDIRK chemistry and force map.
    let mut manual = make();
    let mut native_operator = prepared(&manual, 40.0);
    let first = native_operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), molecular_control, || {
            false
        })
        .unwrap();
    let mut central = reactive(manual, &mechanism);
    let mut native_force = central.flow().forcing_workspace(usize::MAX).unwrap();
    let middle = central
        .step_with_forcing(
            Time::new(dt),
            control,
            force,
            cap,
            &mut native_force,
            || false,
        )
        .unwrap();
    assert_eq!(middle.elapsed.value(), dt);
    let mut manual = central.into_flow();
    let second = native_operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), molecular_control, || {
            false
        })
        .unwrap();

    let mut flow = reactive(make(), &mechanism);
    let mut operator = prepared(flow.flow(), 40.0);
    let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
    let bytes =
        flow.required_stabilized_workspace_bytes(&operator).unwrap() + work.workspace_bytes();
    let report = flow
        .step_with_stabilized_dissipation_and_forcing(
            &mut operator,
            Time::new(dt),
            control,
            molecular_control,
            bytes,
            force,
            cap,
            &mut work,
            || false,
        )
        .unwrap();
    let physical = report.coupled;
    let split = physical.dissipation.unwrap();
    coupled_balance(&physical);
    assert_eq!(physical.coupling_retries, 0);
    assert_eq!(split.first, first.transport);
    assert_eq!(split.second, second.transport);
    assert_eq!(physical.transport, middle.transport);
    assert_eq!(physical.chemistry, middle.chemistry);
    assert_eq!(flow.flow().states(), manual.states());
    assert_eq!(flow.flow().time(), manual.time());
    assert_eq!(flow.last_step(), Some(&physical));
    assert_eq!(operator.last_step(), Some(&second.transport));
    assert_eq!(
        molecular_impulses(&operator),
        molecular_impulses(&native_operator)
    );
    assert_eq!(
        report.operator_evaluations,
        first.operator_evaluations + second.operator_evaluations
    );
    assert_eq!(
        report.maximum_stages_used,
        first.maximum_stages_used.max(second.maximum_stages_used)
    );
    assert!(report.maximum_stages_used > 2);
    assert!(split.first.maximum_parabolic_number > molecular_control.transport.parabolic_cfl);
    assert!(
        report.maximum_scaled_parabolic_number
            <= molecular_control.transport.parabolic_cfl * (1.0 + 1e-13)
    );
    assert_eq!(budget(&flow, &operator) + work.workspace_bytes(), bytes);

    assert_eq!(physical.chemistry.cell_solves, 2 * flow.flow().cell_count());
    assert!(physical.chemistry.species_change_kg[0] < 0.0);
    assert!(physical.chemistry.species_change_kg[1] > 0.0);
    assert!(physical.chemistry.referenced_energy_change_j > 0.0);
    for index in 0..flow.flow().cell_count() {
        assert_eq!(flow.flow().cell(index).unwrap().partial_densities[2], 0.0);
    }
    let mass = physical.before.mass();
    let mut expected_work = 0.0;
    for axis in 0..2 {
        let dv = acceleration[axis] * dt;
        close(physical.source_change.impulse[axis], mass * dv, 3e-12);
        expected_work += mass * dv * (velocity[axis] + 0.5 * dv);
    }
    assert!(expected_work < 0.0);
    close(physical.source_change.work, expected_work, 3e-12);
    close(
        physical.after.absolute_energy(flow.flow().gas()).unwrap()
            - physical.before.absolute_energy(flow.flow().gas()).unwrap()
            - physical.source_change.work,
        0.0,
        3e-10,
    );
}

fn seeded(
    mechanism: &IdealGasKineticsV1,
) -> (
    ReactingEuler2d,
    MolecularTransport2d<ConstantTransport>,
    usize,
) {
    let mut flow = uniform(mechanism);
    let mut operator = prepared(flow.flow(), 0.02);
    let bytes = flow.required_stabilized_workspace_bytes(&operator).unwrap();
    flow.step_with_stabilized_dissipation(
        &mut operator,
        Time::new(1e-5),
        ReactiveStepControl::default(),
        StabilizedDissipationControl::default(),
        bytes,
        || false,
    )
    .unwrap();
    assert_eq!(budget(&flow, &operator), bytes);
    (flow, operator, bytes)
}

#[test]
fn g4_heat_driven_restarts_share_molecular_work_and_publish_neither_owner_on_refusal() {
    let mechanism = make_mechanism(40.0, 500.0, 3.5);
    let controls = ReactiveStepControl::default();
    let molecular = StabilizedDissipationControl::default();
    let (mut reference, mut reference_operator, bytes) = seeded(&mechanism);
    let mut polls = 0;
    let expected = reference
        .step_with_stabilized_dissipation(
            &mut reference_operator,
            Time::new(0.2),
            controls,
            molecular,
            bytes,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    coupled_balance(&expected.coupled);
    assert!(expected.coupled.coupling_retries > 0);
    let (mut replay, mut replay_operator, bytes) = seeded(&mechanism);
    let accepted = replay
        .step_with_stabilized_dissipation(
            &mut replay_operator,
            expected.coupled.elapsed,
            controls,
            molecular,
            bytes,
            || false,
        )
        .unwrap();
    assert_eq!(accepted.coupled.coupling_retries, 0);
    assert_eq!(accepted.coupled.dissipation, expected.coupled.dissipation);
    assert_eq!(accepted.coupled.transport, expected.coupled.transport);
    assert_eq!(accepted.coupled.chemistry, expected.coupled.chemistry);
    assert_eq!(replay.flow().states(), reference.flow().states());
    assert!(expected.operator_evaluations > accepted.operator_evaluations);

    let (mut flow, mut operator, bytes) = seeded(&mechanism);
    let before = published(&flow, &operator);
    for (duration, limit) in [
        (Time::new(0.2), accepted.operator_evaluations),
        (expected.coupled.elapsed, accepted.operator_evaluations - 1),
    ] {
        let mut bounded = molecular;
        bounded.maximum_operator_evaluations = limit;
        match flow.step_with_stabilized_dissipation(
            &mut operator,
            duration,
            controls,
            bounded,
            bytes,
            || false,
        ) {
            Err(ReactiveError::Dissipation(DissipationError::StepFailure {
                elapsed_time_s,
                last,
                ..
            })) => {
                assert_eq!(elapsed_time_s, 0.0);
                assert!(matches!(
                    *last,
                    DissipationError::OperatorWorkLimit { operator_evaluations, .. }
                        if operator_evaluations == limit
                ));
            }
            _ => panic!("the original work refusal must terminate the complete reactive interval"),
        }
        assert_eq!(published(&flow, &operator), before);
        assert_eq!(budget(&flow, &operator), bytes);
    }

    let (mut interrupted, mut interrupted_operator, bytes) = seeded(&mechanism);
    let before = published(&interrupted, &interrupted_operator);
    let mut count = 0;
    assert_eq!(
        interrupted.step_with_stabilized_dissipation(
            &mut interrupted_operator,
            Time::new(0.2),
            controls,
            molecular,
            bytes,
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(published(&interrupted, &interrupted_operator), before);
    let actual = interrupted
        .step_with_stabilized_dissipation(
            &mut interrupted_operator,
            Time::new(0.2),
            controls,
            molecular,
            bytes,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(
        published(&interrupted, &interrupted_operator),
        published(&reference, &reference_operator)
    );
}
