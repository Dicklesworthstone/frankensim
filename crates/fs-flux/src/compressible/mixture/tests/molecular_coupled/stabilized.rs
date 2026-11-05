//! The original stabilized molecular owners compose without resetting work or publication.
use super::*;
use crate::compressible::mixture::molecular::{
    MolecularReport, MolecularTransport2d, PlanarTransport, StabilizedMolecularStepControl,
};

fn prepared(
    flow: &MixtureEuler2d,
    diffusion: f64,
    shear: f64,
    bulk: f64,
    conductivity: f64,
) -> MolecularTransport2d<ConstantTransport> {
    let source = PlanarTransport::new(
        flow.gas().clone(),
        law(flow.gas(), diffusion, shear, bulk, conductivity),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport2d::required_stabilized_workspace_bytes(flow, &source).unwrap();
    MolecularTransport2d::new(flow, source, bytes, || false).unwrap()
}

fn owner(
    flow: MixtureEuler2d,
    molecular: MolecularTransport2d<ConstantTransport>,
) -> (Coupled, ForcingWorkspace2d, usize) {
    let scratch = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let bytes = Coupled::required_stabilized_workspace_bytes(&flow, &molecular).unwrap()
        + scratch.workspace_bytes();
    (
        Coupled::new(flow, molecular, bytes, || false).unwrap(),
        scratch,
        bytes,
    )
}

#[derive(Debug, PartialEq)]
struct Published {
    states: Vec<f64>,
    time: Time,
    split: Option<MolecularStepReport>,
    advection: Option<crate::compressible::mixture::StepReport>,
    molecular: Option<MolecularReport>,
    boundary: Vec<[f64; 6]>,
    molecular_boundary: Vec<[f64; 6]>,
}
fn published(flow: &Coupled) -> Published {
    Published {
        states: flow.flow().states().to_vec(),
        time: flow.flow().time(),
        split: flow.last_step().copied(),
        advection: flow.flow().last_step().copied(),
        molecular: flow.molecular().last_step().copied(),
        boundary: loads(flow),
        molecular_boundary: molecular_loads(flow),
    }
}

#[test]
fn g1_stabilized_mixture_split_preserves_the_original_stiff_diffusion_and_advection() {
    let make = || {
        let grid = CartesianGrid::new(8, 6, 0.1, 0.12, 0.6).unwrap();
        field(equal_calorics(), grid, PERIODIC, |x| {
            let a = 0.5 + 0.1 * (2.0 * std::f64::consts::PI * (x[0] / 0.8 + x[1] / 0.72)).sin();
            ([a, 1.0 - a, 0.0], [0.35, -0.2], 1.0)
        })
    };
    let controls = StabilizedMolecularStepControl::default();
    let dt = 0.004;
    let mut manual = make();
    let mut operator = prepared(&manual, 10.0, 0.0, 0.0, 0.0);
    let first = operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), controls.molecular, || {
            false
        })
        .unwrap();
    let middle = manual
        .step(Time::new(dt), controls.advection, || false)
        .unwrap();
    assert_eq!(middle.elapsed.value(), dt);
    let second = operator
        .advance_stabilized(&mut manual, Time::new(0.5 * dt), controls.molecular, || {
            false
        })
        .unwrap();

    let original = make();
    let molecular = prepared(&original, 10.0, 0.0, 0.0, 0.0);
    let (mut flow, scratch, bytes) = owner(original, molecular);
    let before = flow.flow().states().to_vec();
    let report = flow
        .step_stabilized(Time::new(dt), controls, || false)
        .unwrap();
    check(&report.transport);
    assert_eq!(report.transport.coupling_retries, 0);
    assert_eq!(report.transport.first, first.transport);
    assert_eq!(report.transport.advection, middle);
    assert_eq!(report.transport.second, second.transport);
    assert_eq!(flow.flow().states(), manual.states());
    assert_ne!(flow.flow().states(), before);
    assert_eq!(flow.flow().time(), manual.time());
    assert_eq!(flow.last_step(), Some(&report.transport));
    assert_eq!(flow.molecular().last_step(), Some(&second.transport));
    assert_eq!(
        report.operator_evaluations,
        first.operator_evaluations + second.operator_evaluations
    );
    assert_eq!(
        report.maximum_stages_used,
        first.maximum_stages_used.max(second.maximum_stages_used)
    );
    assert!(report.maximum_stages_used > 2);
    assert!(
        report.transport.first.maximum_parabolic_number
            > controls.molecular.transport.parabolic_cfl
    );
    assert!(
        report.maximum_scaled_parabolic_number
            <= controls.molecular.transport.parabolic_cfl * (1.0 + 1e-13)
    );
    assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
    for index in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.density, 1.0, 2e-11);
        close(p.temperature, 1.0, 2e-10);
        close(p.pressure, 1.0, 2e-10);
        close(p.velocity[0], 0.35, 2e-11);
        close(p.velocity[1], -0.2, 2e-11);
    }

    // An allowance for the ordinary explicit path cannot silently acquire
    // the additional stabilized faces, even with a zero body source.
    let original = make();
    let molecular = prepared(&original, 10.0, 0.0, 0.0, 0.0);
    let mut scratch = original.forcing_workspace(usize::MAX).unwrap();
    let tight = Coupled::required_workspace_bytes(&original, &molecular).unwrap()
        + scratch.workspace_bytes();
    let mut limited = Coupled::new(original, molecular, tight, || false).unwrap();
    let before = published(&limited);
    assert!(matches!(
        limited.step_stabilized_with_forcing(
            Time::new(dt),
            controls,
            BodyForce2d::UniformAcceleration([0.0; 2]),
            CAP,
            &mut scratch,
            || false,
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(published(&limited), before);
}

fn seeded() -> (Coupled, ForcingWorkspace2d) {
    let flow = reservoir_field(31);
    let molecular = prepared(&flow, 0.02, 0.03, 0.01, 0.04);
    let (mut owner, scratch, bytes) = owner(flow, molecular);
    owner
        .step_stabilized(
            Time::new(0.001),
            StabilizedMolecularStepControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(owner.workspace_bytes() + scratch.workspace_bytes(), bytes);
    (owner, scratch)
}

#[test]
fn g4_stabilized_budget_counts_both_halves_and_discarded_forced_intervals_atomically() {
    let controls = StabilizedMolecularStepControl::default();
    let force = BodyForce2d::UniformAcceleration([400.0, -200.0]);
    let cap = ForceControl2d { maximum_mach: 0.8 };
    let (mut reference, mut reference_scratch) = seeded();
    let mut polls = 0;
    let expected = reference
        .step_stabilized_with_forcing(
            Time::new(0.2),
            controls,
            force,
            cap,
            &mut reference_scratch,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    check(&expected.transport);
    final_admission(&reference, &expected.transport, cap);
    assert!(expected.transport.coupling_retries > 0);
    let (mut replay, mut replay_scratch) = seeded();
    let accepted = replay
        .step_stabilized_with_forcing(
            expected.transport.elapsed,
            controls,
            force,
            cap,
            &mut replay_scratch,
            || false,
        )
        .unwrap();
    assert_eq!(accepted.transport.coupling_retries, 0);
    assert_eq!(accepted.transport.first, expected.transport.first);
    assert_eq!(accepted.transport.advection, expected.transport.advection);
    assert_eq!(accepted.transport.second, expected.transport.second);
    assert_eq!(replay.flow().states(), reference.flow().states());
    assert_eq!(loads(&replay), loads(&reference));
    assert_eq!(molecular_loads(&replay), molecular_loads(&reference));
    assert!(expected.operator_evaluations > accepted.operator_evaluations);

    let (mut limited, mut scratch) = seeded();
    let before = published(&limited);
    let bytes = limited.workspace_bytes();
    for (duration, budget) in [
        // Enough work for the final pair, but not the discarded complete trials.
        (Time::new(0.2), accepted.operator_evaluations),
        // No complete retry: the second half still shares the first half's cap.
        (
            expected.transport.elapsed,
            accepted.operator_evaluations - 1,
        ),
    ] {
        let mut bounded = controls;
        bounded.molecular.maximum_operator_evaluations = budget;
        match limited.step_stabilized_with_forcing(
            duration,
            bounded,
            force,
            cap,
            &mut scratch,
            || false,
        ) {
            Err(DissipationError::StepFailure {
                elapsed_time_s,
                last,
                ..
            }) => {
                assert_eq!(elapsed_time_s, 0.0);
                assert!(matches!(
                    *last,
                    DissipationError::OperatorWorkLimit { operator_evaluations, .. }
                        if operator_evaluations == budget
                ));
            }
            _ => panic!("the complete interval must retain its exhausted molecular work cap"),
        }
        assert_eq!(published(&limited), before);
        assert_eq!(limited.workspace_bytes(), bytes);
    }

    // Refuse the final checkpoint after a complete successful private replay.
    let (mut interrupted, mut scratch) = seeded();
    let before = published(&interrupted);
    let mut count = 0;
    assert_eq!(
        interrupted.step_stabilized_with_forcing(
            Time::new(0.2),
            controls,
            force,
            cap,
            &mut scratch,
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(published(&interrupted), before);
    let report = interrupted
        .step_stabilized_with_forcing(Time::new(0.2), controls, force, cap, &mut scratch, || false)
        .unwrap();
    assert_eq!(report, expected);
    assert_eq!(published(&interrupted), published(&reference));
}
