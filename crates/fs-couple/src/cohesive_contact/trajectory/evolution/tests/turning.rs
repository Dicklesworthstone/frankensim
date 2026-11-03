//! Actual two-dimensional Coulomb turning through the original native bodies.
use super::*;
use crate::cohesive_contact::trajectory::sustained::{
    SegmentCandidate, SustainedContact, segment_candidate_controlled,
};

fn trace_rows() -> [[TraceTerm; 1]; 3] {
    core::array::from_fn(|dof| {
        [TraceTerm {
            body: 0,
            dof,
            coefficient: 1.0,
        }]
    })
}
fn scales() -> NaturalMapScales {
    NaturalMapScales {
        normal_n_per_m: 256.0,
        tangent_n_per_m: 256.0,
    }
}

/// Independent closed-form midpoint mechanics for m=2 and Coulomb cap=2.
/// The production owner solves the original nonlinear force equations.
fn midpoint(
    velocity: [f64; 2],
    authored_force: [f64; 2],
    duration_s: f64,
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let w = [
        velocity[0] + 0.25 * duration_s * authored_force[0],
        velocity[1] + 0.25 * duration_s * authored_force[1],
    ];
    let speed = w[0].hypot(w[1]);
    let shrink = 1.0 - duration_s / (2.0 * speed);
    assert!(shrink > 0.0, "accepted segment must retain actual sliding");
    let mean = [shrink * w[0], shrink * w[1]];
    (
        [2.0 * mean[0] - velocity[0], 2.0 * mean[1] - velocity[1]],
        [duration_s * mean[0], duration_s * mean[1]],
        [-2.0 * w[0] / speed, -2.0 * w[1] / speed],
    )
}
fn direct<'m, 'a>(
    before: &CohesiveMotion<'m, 'a>,
    authored_force: [f64; 2],
    scalar_stop_trials: usize,
) -> SegmentCandidate<'m, 'a> {
    let rows = trace_rows();
    let mut segment = controls(2.0).sustained;
    segment.stopping.maximum_trials = scalar_stop_trials;
    segment_candidate_controlled(
        &[before],
        &[SustainedContact {
            id: 701,
            gap_m: 0.0,
            activation_tolerance_m: 1e-8,
            friction_coefficient: 0.1,
            prescribed_velocity_m_s: [0.0; 3],
            trace: [&rows[0], &rows[1], &rows[2]],
        }],
        &[BodyLoad {
            mean_force_n: &[-20.0, authored_force[0], authored_force[1], HELD_FORCE_N],
            endpoint_prescribed_jumps_m: &[[0.0; 3]],
        }],
        &[scales()],
        &[[20.0, -2.0, 0.0]],
        2.0,
        segment,
        &mut || false,
    )
    .unwrap()
}
fn check_direct(
    before: &CohesiveMotion<'_, '_>,
    candidate: &SegmentCandidate<'_, '_>,
    authored_force: [f64; 2],
) {
    let duration = candidate.report.accepted_duration_s;
    close(duration, 0.25, 1e-12);
    assert!(candidate.report.stopping.is_none());
    assert!(candidate.report.stopped_contact_ids.is_empty());
    let turning = candidate.report.turning.as_ref().unwrap();
    assert_eq!(turning.trigger_contact_id, 701);
    assert!(turning.refinements >= 1);
    close(turning.elapsed_s, duration, 0.0);
    let before_velocity = [before.velocity_m_s()[1], before.velocity_m_s()[2]];
    let (velocity, slip, friction) = midpoint(before_velocity, authored_force, duration);
    let after = &candidate.step.bodies[0].motion;
    for axis in 0..2 {
        close(after.velocity_m_s()[axis + 1], velocity[axis], 1e-8);
        close(after.state().displacements_m()[axis + 1], slip[axis], 1e-9);
        close(
            candidate.step.report.contact.forces_n[0][axis + 1],
            friction[axis],
            1e-8,
        );
        close(
            2.0 * (after.velocity_m_s()[axis + 1] - before_velocity[axis]),
            duration * (authored_force[axis] + friction[axis]),
            1e-8,
        );
    }
    let cross = before_velocity[0] * velocity[1] - before_velocity[1] * velocity[0];
    let dot = before_velocity[0] * velocity[0] + before_velocity[1] * velocity[1];
    let angle = cross.abs().atan2(dot);
    assert!(angle > 0.0 && angle <= core::f64::consts::FRAC_PI_4);
    close(turning.maximum_angle_rad, angle, 1e-8);
    assert!(after.velocity_m_s()[1].hypot(after.velocity_m_s()[2]) > 1e-3);
    close(candidate.step.report.contact.forces_n[0][0], 20.0, 1e-8);
    let work = authored_force[0] * slip[0] + authored_force[1] * slip[1];
    let loss = 2.0 * slip[0].hypot(slip[1]);
    let kinetic_change = velocity[0] * velocity[0] + velocity[1] * velocity[1]
        - before_velocity[0] * before_velocity[0]
        - before_velocity[1] * before_velocity[1];
    close(work, kinetic_change + loss, 1e-10);
    close(candidate.step.report.audit.friction_loss_j, loss, 1e-8);
    close(candidate.step.report.audit.energy_residual_j, 0.0, 1e-8);
    assert_eq!(after.time_s().to_bits(), (START_S + duration).to_bits());
    retained(after, before);
}

#[test]
fn transverse_authored_load_resolves_a_real_turn_without_inventing_a_stop() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0, 1.2, 0.0], 2.0);
    let original = before.clone();
    let candidate = direct(&before, [-1.2, 4.8], 64);
    check_direct(&before, &candidate, [-1.2, 4.8]);
    unchanged(&before, &original);
}

#[test]
fn a_transverse_intermediate_response_intercepts_an_apparent_half_turn_stop() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = motion(&system, &mass, [0.0, 1.0, 0.0], 2.0);
    let original = before.clone();
    // The full h=2 response has v+=(-1,0), yet its h=1 stopping
    // probe has v+=(-1/sqrt(2), 1-1/sqrt(2)), a resolved turn.
    // One scalar probe is enough only if the actual transverse response is
    // handed back to angular refinement instead of a fabricated zero speed.
    let candidate = direct(&before, [-2.0, 2.0], 1);
    check_direct(&before, &candidate, [-2.0, 2.0]);
    unchanged(&before, &original);
}

#[test]
fn whole_curved_sliding_windows_preserve_native_momentum_and_work_through_every_turn() {
    for (initial, authored_force) in [(1.2, [-1.2, 4.8]), (1.0, [-2.0, 2.0])] {
        let elastic = HeldElastic;
        let mass = Mass;
        let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
        let system = system(&elastic, &law);
        let before = motion(&system, &mass, [0.0, initial, 0.0], 2.0);
        let original = before.clone();
        let rows = trace_rows();
        let contact = FrictionalContact {
            id: 701,
            gap_m: 0.0,
            activation_tolerance_m: 1e-8,
            friction_coefficient: 0.1,
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
            prescribed_velocity_m_s: [0.0; 3],
            trace: [&rows[0], &rows[1], &rows[2]],
        };
        let preconditioner = IdentityPrecond;
        let candidate = advance_candidate_controlled(
            &[&before],
            &[contact],
            &[BodyLoad {
                mean_force_n: &[-20.0, authored_force[0], authored_force[1], HELD_FORCE_N],
                endpoint_prescribed_jumps_m: &[[0.0; 3]],
            }],
            &[&preconditioner],
            &[scales()],
            &[[20.0, -2.0, 0.0]],
            2.0,
            controls(2.0),
            &mut || false,
        )
        .unwrap();
        let report = &candidate.report;
        let after = &candidate.motions[0];
        assert!(!report.turning_segments.is_empty());
        assert!(report.stopping_events.is_empty());
        assert!(report.arrival_events.is_empty());
        assert!(report.trajectory.events.is_empty());
        // With maximum_step=the entire window, every nonfinal accepted
        // interval has an angular receipt. Reconstruct independent midpoint
        // mechanics at those actual times, then the final unshortened interval.
        let mut ends: Vec<f64> = report.turning_segments.iter().map(|s| s.time_s).collect();
        if ends.last().copied() != Some(START_S + 2.0) {
            ends.push(START_S + 2.0);
        }
        let mut time_s = START_S;
        let mut velocity = [initial, 0.0];
        let mut displacement = [0.0; 2];
        let mut work = 0.0;
        let mut loss = 0.0;
        for endpoint_s in ends {
            assert!(endpoint_s > time_s);
            let (next, slip, _) = midpoint(velocity, authored_force, endpoint_s - time_s);
            for axis in 0..2 {
                displacement[axis] += slip[axis];
            }
            work += authored_force[0] * slip[0] + authored_force[1] * slip[1];
            loss += 2.0 * slip[0].hypot(slip[1]);
            velocity = next;
            time_s = endpoint_s;
        }
        for turning in &report.turning_segments {
            assert_eq!(turning.turning.trigger_contact_id, 701);
            assert!(turning.turning.maximum_angle_rad <= core::f64::consts::FRAC_PI_4);
            assert!(turning.turning.elapsed_s > 0.0);
        }
        for axis in 0..2 {
            close(after.velocity_m_s()[axis + 1], velocity[axis], 1e-8);
            close(
                after.state().displacements_m()[axis + 1],
                displacement[axis],
                1e-8,
            );
        }
        close(report.trajectory.external_work_j, work, 1e-8);
        close(report.contact_friction_loss_j, loss, 1e-8);
        close(
            report.trajectory.kinetic_change_j,
            velocity[0] * velocity[0] + velocity[1] * velocity[1] - initial * initial,
            1e-8,
        );
        close(report.trajectory.energy_residual_j, 0.0, 1e-8);
        close(after.velocity_m_s()[0], 0.0, 1e-8);
        close(report.trajectory.final_gaps_m[0], 0.0, 1e-10);
        assert_eq!(after.time_s().to_bits(), (START_S + 2.0).to_bits());
        retained(after, &before);
        unchanged(&before, &original);
    }
}
