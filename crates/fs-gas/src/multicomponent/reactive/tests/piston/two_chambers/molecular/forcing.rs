//! G1/G3/G4: two independently forced molecular chambers sharing one piston.
use super::*;
use crate::body_force::{BodyForce, ForceControl};
use crate::dynamics::piston::PistonSuspension;

const FORCE_CONTROL: ForceControl = ForceControl { maximum_mach: 80.0 };

fn small_pair() -> (ReactingEuler1d, ReactingEuler1d) {
    (
        chamber(&mechanism(8.0, 30.0, 100.0), 0.0, 0.8, 2, chemistry()),
        chamber(&mechanism(2.0, 8.0, -30.0), 0.8, 1.2, 2, chemistry()),
    )
}

#[allow(clippy::too_many_arguments)]
fn forced_step(
    left: &mut ReactingEuler1d,
    right: &mut ReactingEuler1d,
    dl: &mut DissipationOperator<ConstantTransport>,
    dr: &mut DissipationOperator<ConstantTransport>,
    piston: &mut Piston,
    forces: [BodyForce<'_>; 2],
    cancelled: impl FnMut() -> bool,
) -> Result<PistonStep<ReactiveSharedPistonGasReport>, ReactiveError> {
    left.step_with_shared_piston_dissipation_and_body_force(
        right,
        dl,
        dr,
        piston,
        Time::new(0.001),
        ReactiveStepControl::default(),
        [DissipationControl::default(); 2],
        0.01,
        piston_controls(),
        forces,
        [FORCE_CONTROL; 2],
        cancelled,
    )
}

#[test]
fn g1_two_forced_chambers_close_chemical_spring_and_kinetic_energy_with_distinct_heat() {
    let (mut left, mut right) = small_pair();
    let mut dl = molecular(left.flow());
    let mut dr = molecular(right.flow());
    let mut piston = Piston::new(0.2, 0.2).unwrap();
    let law = PistonSuspension::new(1.5, 0.3, 0.75).unwrap();
    let before = [left.flow().totals(), right.flow().totals()];
    let step = left
        .step_with_shared_suspended_piston_dissipation_and_body_force(
            &mut right,
            &mut dl,
            &mut dr,
            &mut piston,
            Time::new(0.001),
            ReactiveStepControl::default(),
            [DissipationControl::default(); 2],
            0.01,
            law,
            piston_controls(),
            [
                BodyForce::UniformAcceleration(3.0),
                BodyForce::UniformForceDensity(-2.0),
            ],
            [FORCE_CONTROL; 2],
            || false,
        )
        .unwrap();
    let report = &step.piston;
    let h = report.gas.left.elapsed.value();
    let x1 = left.flow().mesh().positions()[2];
    let displacement = x1 - 0.8;
    let spring_change = 0.5 * 1.5 * ((x1 - 0.75).powi(2) - (0.8_f64 - 0.75).powi(2));
    let midpoint = report.midpoint_velocity_m_s;
    let damping_heat = 0.3 * h * midpoint.powi(2);
    let kinetic = 0.5 * piston.mass_kg() * (piston.velocity_m_s().powi(2) - 0.2_f64.powi(2));
    let sources = [
        report.gas.left.transport.source_change,
        report.gas.right.transport.source_change,
    ];
    close(sources[0].impulse, before[0].mass() * 3.0 * h, 3e-12);
    let right_volume: f64 = right.flow().mesh().volumes().iter().sum();
    close(
        sources[1].impulse,
        0.5 * h * -2.0 * (AREA * 1.2 + right_volume),
        3e-12,
    );
    close(
        physical_change(left.flow(), before[0])
            + physical_change(right.flow(), before[1])
            + kinetic
            + spring_change
            + damping_heat
            - 0.01 * displacement
            - sources.iter().map(|s| s.work).sum::<f64>(),
        0.0,
        4e-8,
    );
    let suspension_impulse = -h * (1.5 * (0.8 + 0.5 * displacement - 0.75) + 0.3 * midpoint);
    close(
        left.flow().totals().momentum - before[0].momentum + right.flow().totals().momentum
            - before[1].momentum
            + piston.mass_kg() * (piston.velocity_m_s() - 0.2)
            + report.fixed_wall_impulse_n_s
            - 0.01 * h
            - suspension_impulse
            - sources.iter().map(|s| s.impulse).sum::<f64>(),
        0.0,
        4e-11,
    );
    for (side, flow) in [&report.gas.left, &report.gas.right]
        .into_iter()
        .enumerate()
    {
        close(
            report.gas.gas_energy_change_j[side],
            physical_change(
                if side == 0 { left.flow() } else { right.flow() },
                before[side],
            ),
            3e-11,
        );
        assert!(flow.chemistry.species_change_kg[1] > 0.0);
        let halves = flow.dissipation.unwrap();
        for half in [halves.first, halves.second] {
            let ends = [half.left_boundary_change, half.right_boundary_change];
            close(
                ends[1 - side].energy,
                midpoint * ends[1 - side].momentum,
                3e-13,
            );
        }
    }
    assert_eq!(left.flow().time(), right.flow().time());
    assert_eq!(x1, right.flow().mesh().positions()[0]);
    assert_eq!(dl.mesh(), left.flow().mesh());
    assert_eq!(dr.mesh(), right.flow().mesh());
}

#[test]
fn g3_zero_sources_replay_original_shared_molecular_root() {
    let (mut left, mut right) = small_pair();
    let (mut a, mut b) = small_pair();
    let mut dl = molecular(left.flow());
    let mut dr = molecular(right.flow());
    let mut da = molecular(a.flow());
    let mut db = molecular(b.flow());
    let mut p = Piston::new(0.2, 0.2).unwrap();
    let mut q = p;
    let actual = forced_step(
        &mut left,
        &mut right,
        &mut dl,
        &mut dr,
        &mut p,
        [
            BodyForce::Acceleration(&[0.0; 2]),
            BodyForce::UniformForceDensity(0.0),
        ],
        || false,
    )
    .unwrap();
    let expected = molecular_step(
        &mut a,
        &mut b,
        &mut da,
        &mut db,
        &mut q,
        Time::new(0.001),
        || false,
    )
    .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(left.accepted.states, a.accepted.states);
    assert_eq!(right.accepted.states, b.accepted.states);
    assert_eq!(p, q);
}

#[test]
fn g4_two_forced_fields_and_both_molecular_bindings_publish_atomically() {
    let forces = [
        BodyForce::UniformAcceleration(3.0),
        BodyForce::UniformAcceleration(-2.0),
    ];
    let (mut reference_left, mut reference_right) = small_pair();
    let mut dl = molecular(reference_left.flow());
    let mut dr = molecular(reference_right.flow());
    let mut p = Piston::new(0.2, 0.2).unwrap();
    let mut polls = 0;
    let expected = forced_step(
        &mut reference_left,
        &mut reference_right,
        &mut dl,
        &mut dr,
        &mut p,
        forces,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for cutoff in [polls / 2, polls] {
        let (mut left, mut right) = small_pair();
        let mut dl = molecular(left.flow());
        let mut dr = molecular(right.flow());
        let mut piston = Piston::new(0.2, 0.2).unwrap();
        let states = [left.accepted.states.clone(), right.accepted.states.clone()];
        let meshes = [left.flow().mesh().clone(), right.flow().mesh().clone()];
        let original = piston;
        let mut seen = 0;
        let error = forced_step(
            &mut left,
            &mut right,
            &mut dl,
            &mut dr,
            &mut piston,
            forces,
            || {
                seen += 1;
                seen == cutoff
            },
        )
        .unwrap_err();
        assert_eq!(error, ReactiveError::Cancelled);
        assert_eq!(seen, cutoff);
        assert_eq!(left.accepted.states, states[0]);
        assert_eq!(right.accepted.states, states[1]);
        assert_eq!(left.flow().mesh(), &meshes[0]);
        assert_eq!(right.flow().mesh(), &meshes[1]);
        assert_eq!(dl.mesh(), &meshes[0]);
        assert_eq!(dr.mesh(), &meshes[1]);
        assert_eq!(left.flow().time().value(), 0.0);
        assert_eq!(right.flow().time().value(), 0.0);
        assert_eq!(piston, original);
        let resumed = forced_step(
            &mut left,
            &mut right,
            &mut dl,
            &mut dr,
            &mut piston,
            forces,
            || false,
        )
        .unwrap();
        assert_eq!(resumed, expected);
        assert_eq!(left.accepted.states, reference_left.accepted.states);
        assert_eq!(right.accepted.states, reference_right.accepted.states);
        assert_eq!(piston, p);
    }
}
