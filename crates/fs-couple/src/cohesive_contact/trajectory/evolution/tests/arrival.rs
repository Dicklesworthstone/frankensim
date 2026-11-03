//! Actual native arrival states before impact and supported continuation.
use super::*;
use crate::cohesive_contact::trajectory::TrajectoryError;

const MAXIMUM_STEP_S: f64 = 0.375;
const DURATION_S: f64 = 0.5;

// The broad phase admits candidates a millimetre away. The trajectory must
// still resolve the physical endpoint to the separate, much tighter gap gate.
const ACTIVATION_GAP_M: f64 = 1e-3;

fn approach<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    gaps_m: &[f64],
    surface_normal_velocity_m_s: f64,
    controls: EvolutionControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EvolutionCandidate<'m, 'a>, TrajectoryError> {
    assert_eq!(motions.len(), gaps_m.len());
    let rows: Vec<[[TraceTerm; 1]; 3]> = (0..motions.len())
        .map(|body| {
            core::array::from_fn(|dof| {
                [TraceTerm {
                    body,
                    dof,
                    coefficient: 1.0,
                }]
            })
        })
        .collect();
    let contacts: Vec<_> = rows
        .iter()
        .zip(gaps_m)
        .enumerate()
        .map(|(body, (rows, &gap_m))| FrictionalContact {
            id: 101 + body as u64,
            gap_m,
            activation_tolerance_m: ACTIVATION_GAP_M,
            friction_coefficient: 0.25,
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
            // Contact convention is Bv+c, so c is minus world wall velocity.
            prescribed_velocity_m_s: [-surface_normal_velocity_m_s, 0.0, 0.0],
            trace: [&rows[0], &rows[1], &rows[2]],
        })
        .collect();
    let forces = vec![[-8.0, 0.0, 0.0, HELD_FORCE_N]; motions.len()];
    let offsets = vec![[[0.0; 3]]; motions.len()];
    let loads: Vec<_> = forces
        .iter()
        .zip(&offsets)
        .map(|(force, offsets)| BodyLoad {
            mean_force_n: force,
            endpoint_prescribed_jumps_m: offsets,
        })
        .collect();
    let preconditioner = IdentityPrecond;
    let preconditioners: Vec<&dyn fs_sparse::precond::Precond> =
        (0..motions.len()).map(|_| &preconditioner as _).collect();
    let scales = vec![
        NaturalMapScales {
            normal_n_per_m: 256.0,
            tangent_n_per_m: 256.0,
        };
        contacts.len()
    ];
    let initial_forces = vec![[0.0; 3]; contacts.len()];
    advance_candidate_controlled(
        motions,
        &contacts,
        &loads,
        &preconditioners,
        &scales,
        &initial_forces,
        DURATION_S,
        controls,
        cancelled,
    )
}

#[test]
fn positive_gap_accelerated_arrival_impacts_at_the_physical_surface_then_stays_supported() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [-0.5, 0.0, 0.0], MAXIMUM_STEP_S);
    let original = before.clone();
    let candidate = approach(
        &[&before],
        &[3.0 / 32.0],
        0.0,
        controls(MAXIMUM_STEP_S),
        &mut || false,
    )
    .unwrap();
    let after = &candidate.motions[0];
    let report = &candidate.report;
    // m=2, a=-4: g(t)=3/32-t/2-2t^2 first vanishes at t=1/8.
    // The incoming velocity is -1, hence J=2 and plastic impact loss=1.
    assert_eq!(report.trajectory.events.len(), 1);
    assert_eq!(report.arrival_events.len(), 1);
    let arrival = &report.arrival_events[0];
    assert_eq!(arrival.contact_id, 101);
    assert!(arrival.arrival.gap_m >= 0.0);
    assert!(arrival.arrival.gap_m <= controls(MAXIMUM_STEP_S).arrival.gap_tolerance_m);
    assert_eq!(
        arrival.time_s.to_bits(),
        report.trajectory.events[0].time_s.to_bits()
    );
    assert_eq!(arrival.segment_start_time_s.to_bits(), START_S.to_bits());
    let analytical_time_s = START_S + 0.125;
    assert!(arrival.bracket_time_s[0] <= analytical_time_s);
    assert!(arrival.bracket_time_s[1] >= analytical_time_s);
    assert!(
        arrival.bracket_time_s[1] - arrival.bracket_time_s[0]
            <= controls(MAXIMUM_STEP_S).arrival.time_tolerance_s
    );
    assert!(arrival.arrival.response_evaluations >= 1);
    close(arrival.arrival.normal_velocity_m_s, -1.0, 1e-8);
    close(report.trajectory.events[0].time_s, START_S + 0.125, 1e-9);
    close(
        report.trajectory.events[0].contact.impulses[0][0],
        2.0,
        1e-8,
    );
    close(after.state().displacements_m()[0], -3.0 / 32.0, 1e-10);
    close(after.velocity_m_s()[0], 0.0, 1e-8);
    close(report.trajectory.final_gaps_m[0], 0.0, 1e-10);
    assert!(report.trajectory.minimum_gap_m >= 0.0);
    assert_eq!(after.time_s().to_bits(), (START_S + DURATION_S).to_bits());
    assert!(report.sustained_steps > 0);
    close(report.warm_start_forces_n[0][0], 8.0, 1e-8);
    close(report.trajectory.external_work_j, 0.75, 1e-8);
    close(report.trajectory.kinetic_change_j, -0.25, 1e-8);
    close(report.trajectory.impact_loss_j, 1.0, 1e-8);
    close(report.trajectory.contact_prescribed_work_j, 0.0, 1e-10);
    close(report.contact_friction_loss_j, 0.0, 1e-10);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    retained(after, &before);
    unchanged(&before, &original);
}

#[test]
fn moving_surface_arrival_uses_relative_position_and_retains_world_momentum_and_work() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [-0.5, 0.0, 0.0], MAXIMUM_STEP_S);
    let original = before.clone();
    let candidate = approach(
        &[&before],
        &[5.0 / 32.0],
        0.5,
        controls(MAXIMUM_STEP_S),
        &mut || false,
    )
    .unwrap();
    let after = &candidate.motions[0];
    let report = &candidate.report;
    // g(t)=5/32-t-2t^2 vanishes at 1/8. The actual body then has
    // v=-1 and acquires J=3 to match the moving surface's +1/2.
    assert_eq!(report.trajectory.events.len(), 1);
    close(report.trajectory.events[0].time_s, START_S + 0.125, 1e-9);
    close(
        report.trajectory.events[0].contact.impulses[0][0],
        3.0,
        1e-8,
    );
    close(after.velocity_m_s()[0], 0.5, 1e-8);
    close(after.state().displacements_m()[0], 3.0 / 32.0, 1e-10);
    close(report.trajectory.final_gaps_m[0], 0.0, 1e-10);
    assert!(report.trajectory.minimum_gap_m >= 0.0);
    assert_eq!(after.time_s().to_bits(), (START_S + DURATION_S).to_bits());
    // The surface supplies 3/2 J at impact and 3/2 J during support.
    close(report.trajectory.contact_prescribed_work_j, 3.0, 1e-8);
    close(report.trajectory.external_work_j, -0.75, 1e-8);
    close(report.trajectory.impact_loss_j, 2.25, 1e-8);
    close(report.trajectory.kinetic_change_j, 0.0, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    close(report.contact_friction_loss_j, 0.0, 1e-10);
    close(report.warm_start_forces_n[0][0], 8.0, 1e-8);
    retained(after, &before);
    unchanged(&before, &original);
}

#[test]
fn second_body_arrives_while_the_first_contact_is_already_supported() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let first = motion(&system, &mass, [-0.5, 0.0, 0.0], MAXIMUM_STEP_S);
    let second = motion(&system, &mass, [-0.5, 0.0, 0.0], MAXIMUM_STEP_S);
    let first_original = first.clone();
    let second_original = second.clone();
    let candidate = approach(
        &[&first, &second],
        &[3.0 / 32.0, 0.25],
        0.0,
        controls(MAXIMUM_STEP_S),
        &mut || false,
    )
    .unwrap();
    let report = &candidate.report;
    assert_eq!(report.trajectory.events.len(), 2);
    close(report.trajectory.events[0].time_s, START_S + 0.125, 1e-9);
    close(report.trajectory.events[1].time_s, START_S + 0.25, 1e-9);
    close(
        report.trajectory.events[0].contact.impulses[0][0],
        2.0,
        1e-8,
    );
    close(
        report.trajectory.events[0].contact.impulses[1][0],
        0.0,
        1e-8,
    );
    close(
        report.trajectory.events[1].contact.impulses[0][0],
        0.0,
        1e-8,
    );
    close(
        report.trajectory.events[1].contact.impulses[1][0],
        3.0,
        1e-8,
    );
    for (body, gap) in [3.0 / 32.0, 0.25].into_iter().enumerate() {
        let after = &candidate.motions[body];
        close(after.velocity_m_s()[0], 0.0, 1e-8);
        close(after.state().displacements_m()[0], -gap, 1e-10);
        close(report.trajectory.final_gaps_m[body], 0.0, 1e-10);
        close(report.warm_start_forces_n[body][0], 8.0, 1e-8);
        assert_eq!(after.time_s().to_bits(), (START_S + DURATION_S).to_bits());
        retained(after, [&first, &second][body]);
    }
    assert!(report.trajectory.minimum_gap_m >= 0.0);
    close(report.trajectory.external_work_j, 2.75, 1e-8);
    close(report.trajectory.kinetic_change_j, -0.5, 1e-8);
    close(report.trajectory.impact_loss_j, 3.25, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    unchanged(&first, &first_original);
    unchanged(&second, &second_original);
}

#[test]
fn rejected_arrival_trials_and_late_cancellation_publish_no_partial_native_state() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [-0.5, 0.0, 0.0], MAXIMUM_STEP_S);
    let original = before.clone();
    let mut polls = 0;
    let complete = approach(
        &[&before],
        &[3.0 / 32.0],
        0.0,
        controls(MAXIMUM_STEP_S),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(complete.report.trajectory.body_trials > 1);
    assert!(polls > 10);

    let mut unfunded = controls(MAXIMUM_STEP_S);
    unfunded.trajectory.maximum_body_trials = 1;
    assert!(matches!(
        approach(&[&before], &[3.0 / 32.0], 0.0, unfunded, &mut || false),
        Err(TrajectoryError::Budget(_))
    ));
    unchanged(&before, &original);

    let mut seen = 0;
    let refused = approach(
        &[&before],
        &[3.0 / 32.0],
        0.0,
        controls(MAXIMUM_STEP_S),
        &mut || {
            seen += 1;
            seen >= polls - 2
        },
    );
    assert!(refused.is_err());
    assert!(seen >= polls - 2);
    unchanged(&before, &original);
}

#[test]
fn positive_gap_inside_broad_activation_still_flies_to_the_physical_endpoint() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [-0.5, 0.0, 0.0], MAXIMUM_STEP_S);
    let original = before.clone();
    let gap_m = 257.0 / 524_288.0;
    assert!(gap_m < ACTIVATION_GAP_M);
    let candidate = approach(
        &[&before],
        &[gap_m],
        0.0,
        controls(MAXIMUM_STEP_S),
        &mut || false,
    )
    .unwrap();
    let report = &candidate.report;
    // The millimetre broad-phase tolerance is not permission to apply an
    // impulse now: the true constant-acceleration arrival is at t=1/1024.
    assert_eq!(report.trajectory.events.len(), 1);
    close(
        report.trajectory.events[0].time_s,
        START_S + 1.0 / 1024.0,
        1e-9,
    );
    close(
        report.trajectory.events[0].contact.impulses[0][0],
        129.0 / 128.0,
        1e-8,
    );
    close(
        candidate.motions[0].state().displacements_m()[0],
        -gap_m,
        1e-10,
    );
    close(candidate.motions[0].velocity_m_s()[0], 0.0, 1e-8);
    close(report.trajectory.final_gaps_m[0], 0.0, 1e-10);
    close(report.trajectory.external_work_j, 257.0 / 65_536.0, 1e-8);
    close(report.trajectory.impact_loss_j, 16_641.0 / 65_536.0, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    assert_eq!(
        candidate.motions[0].time_s().to_bits(),
        (START_S + DURATION_S).to_bits()
    );
    retained(&candidate.motions[0], &before);
    unchanged(&before, &original);
}

#[test]
fn initially_stationary_positive_gap_accelerates_to_impact_before_supported_rest() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0; 3], MAXIMUM_STEP_S);
    let original = before.clone();
    let candidate = approach(
        &[&before],
        &[1.0 / 32.0],
        0.0,
        controls(MAXIMUM_STEP_S),
        &mut || false,
    )
    .unwrap();
    let report = &candidate.report;
    // g(t)=1/32-2t^2: a stationary body still reaches contact at 1/8.
    // Arrival cannot be inferred from the initial gap/closing-speed quotient.
    assert_eq!(report.trajectory.events.len(), 1);
    close(report.trajectory.events[0].time_s, START_S + 0.125, 1e-9);
    close(
        report.trajectory.events[0].contact.impulses[0][0],
        1.0,
        1e-8,
    );
    close(
        candidate.motions[0].state().displacements_m()[0],
        -1.0 / 32.0,
        1e-10,
    );
    close(candidate.motions[0].velocity_m_s()[0], 0.0, 1e-8);
    close(report.trajectory.final_gaps_m[0], 0.0, 1e-10);
    close(report.trajectory.external_work_j, 0.25, 1e-8);
    close(report.trajectory.impact_loss_j, 0.25, 1e-8);
    close(report.trajectory.kinetic_change_j, 0.0, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    assert_eq!(
        candidate.motions[0].time_s().to_bits(),
        (START_S + DURATION_S).to_bits()
    );
    retained(&candidate.motions[0], &before);
    unchanged(&before, &original);
}
