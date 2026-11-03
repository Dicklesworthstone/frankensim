//! Original native owners across coupled impacts and ordered contact events.
use super::*;

struct CoupledMass;
impl LinearOp for CoupledMass {
    fn n(&self) -> usize {
        4
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y[0] = 2.0 * x[0] + x[1];
        y[1] = x[0] + 2.0 * x[1];
        y[2] = 2.0 * x[2];
        y[3] = 2.0 * x[3];
    }
}
fn coupled_motion<'m, 'a>(
    system: &'m CohesiveSystem<'a>,
    mass: &'m CoupledMass,
    maximum_step_s: f64,
) -> CohesiveMotion<'m, 'a> {
    let state = system
        .bonded_state()
        .solve_controlled(
            &[0.0, 0.0, 0.0, HELD_FORCE_N],
            &[[0.0; 3]],
            dynamics(maximum_step_s).solve,
            &mut || false,
        )
        .unwrap()
        .state;
    assert!(state.points()[0].damage() > 0.0 && state.points()[0].damage() < 1.0);
    close(state.displacements_m()[3], HELD_OPENING_M, 1e-11);
    state
        .into_dynamics(
            mass,
            vec![0.0, -1.0, 0.0, 0.0],
            START_S,
            dynamics(maximum_step_s),
        )
        .unwrap()
}
fn rows(body: usize, dofs: [usize; 3]) -> [[TraceTerm; 1]; 3] {
    dofs.map(|dof| {
        [TraceTerm {
            body,
            dof,
            coefficient: 1.0,
        }]
    })
}
fn contact<'r>(
    id: u64,
    gap_m: f64,
    friction_coefficient: f64,
    rows: &'r [[TraceTerm; 1]; 3],
) -> FrictionalContact<'r> {
    FrictionalContact {
        id,
        gap_m,
        activation_tolerance_m: 1e-3,
        friction_coefficient,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
        prescribed_velocity_m_s: [0.0; 3],
        trace: [&rows[0], &rows[1], &rows[2]],
    }
}
fn run<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[FrictionalContact<'_>],
    forces: &[[f64; 4]],
    initial_forces_n: &[[f64; 3]],
    duration_s: f64,
    maximum_step_s: f64,
) -> EvolutionCandidate<'m, 'a> {
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
    advance_candidate_controlled(
        motions,
        contacts,
        &loads,
        &preconditioners,
        &scales,
        initial_forces_n,
        duration_s,
        controls(maximum_step_s),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn arriving_contact_resolves_simultaneous_momentum_in_an_existing_mass_coupled_support() {
    let elastic = HeldElastic;
    let mass = CoupledMass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let before = coupled_motion(&system, &mass, 0.25);
    let original = before.clone();
    let x = rows(0, [0, 1, 2]);
    let y = rows(0, [1, 0, 2]);
    let contacts = [contact(401, 0.0, 0.0, &x), contact(402, 0.125, 0.0, &y)];
    let candidate = run(
        &[&before],
        &contacts,
        &[[-20.0, 0.0, 0.0, HELD_FORCE_N]],
        &[[20.0, 0.0, 0.0], [0.0; 3]],
        0.5,
        0.25,
    );
    let after = &candidate.motions[0];
    let report = &candidate.report;
    assert_eq!(report.trajectory.events.len(), 1);
    let impact = &report.trajectory.events[0];
    close(impact.time_s, START_S + 0.125, 1e-9);
    // M_xy=[[2,1],[1,2]], v^-=[0,-1], v^+=[0,0].
    // M(v^+-v^-)=[1,2]: the already supported x contact must acquire
    // an impulse even though its pre-impact normal velocity is exactly zero.
    close(impact.contact.impulses[0][0], 1.0, 1e-8);
    close(impact.contact.impulses[1][0], 2.0, 1e-8);
    for impulse in &impact.contact.impulses {
        close(impulse[1], 0.0, 1e-10);
        close(impulse[2], 0.0, 1e-10);
    }
    close(after.velocity_m_s()[0], 0.0, 1e-8);
    close(after.velocity_m_s()[1], 0.0, 1e-8);
    close(after.state().displacements_m()[0], 0.0, 1e-10);
    close(after.state().displacements_m()[1], -0.125, 1e-10);
    for &gap in &report.trajectory.final_gaps_m {
        close(gap, 0.0, 1e-10);
    }
    close(report.warm_start_forces_n[0][0], 20.0, 1e-8);
    close(report.warm_start_forces_n[1][0], 0.0, 1e-8);
    close(report.trajectory.external_work_j, 0.0, 1e-8);
    close(report.trajectory.impact_loss_j, 1.0, 1e-8);
    close(report.trajectory.kinetic_change_j, -1.0, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    close(report.contact_friction_loss_j, 0.0, 1e-10);
    assert_eq!(after.time_s().to_bits(), (START_S + 0.5).to_bits());
    assert_eq!(report.arrival_events.len(), 1);
    assert_eq!(report.arrival_events[0].contact_id, 402);
    retained(after, &before);
    unchanged(&before, &original);
}

#[test]
fn actual_sliding_stop_precedes_an_independent_normal_arrival_without_skipping_either() {
    let elastic = HeldElastic;
    let mass = Mass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let slider = motion(&system, &mass, [0.0, 5.0 / 32.0, 0.0], 0.25);
    let incoming = motion(&system, &mass, [-1.0, 0.0, 0.0], 0.25);
    let slider_original = slider.clone();
    let incoming_original = incoming.clone();
    let first = rows(0, [0, 1, 2]);
    let second = rows(1, [0, 1, 2]);
    let contacts = [
        contact(501, 0.0, 0.25, &first),
        contact(502, 0.125, 0.0, &second),
    ];
    let candidate = run(
        &[&slider, &incoming],
        &contacts,
        &[
            [-20.0, 0.0, 0.0, HELD_FORCE_N],
            [0.0, 0.0, 0.0, HELD_FORCE_N],
        ],
        &[[20.0, -5.0, 0.0], [0.0; 3]],
        0.25,
        0.25,
    );
    let report = &candidate.report;
    assert_eq!(report.stopping_events.len(), 1);
    assert_eq!(report.stopping_events[0].contact_ids, vec![501]);
    close(report.stopping_events[0].time_s, START_S + 1.0 / 16.0, 1e-9);
    assert_eq!(report.arrival_events.len(), 1);
    assert_eq!(report.arrival_events[0].contact_id, 502);
    close(report.arrival_events[0].time_s, START_S + 0.125, 1e-9);
    assert_eq!(report.trajectory.events.len(), 1);
    close(report.trajectory.events[0].time_s, START_S + 0.125, 1e-9);
    assert!(report.stopping_events[0].time_s < report.arrival_events[0].time_s);
    close(
        report.trajectory.events[0].contact.impulses[0][0],
        0.0,
        1e-8,
    );
    close(
        report.trajectory.events[0].contact.impulses[1][0],
        2.0,
        1e-8,
    );
    close(
        candidate.motions[0].state().displacements_m()[1],
        5.0 / 1024.0,
        1e-10,
    );
    close(
        candidate.motions[1].state().displacements_m()[0],
        -0.125,
        1e-10,
    );
    for (body, after) in candidate.motions.iter().enumerate() {
        for &velocity in &after.velocity_m_s()[..3] {
            close(velocity, 0.0, 1e-8);
        }
        assert_eq!(after.time_s().to_bits(), (START_S + 0.25).to_bits());
        retained(after, [&slider, &incoming][body]);
    }
    close(report.contact_friction_loss_j, 25.0 / 1024.0, 1e-8);
    close(report.trajectory.impact_loss_j, 1.0, 1e-8);
    close(report.trajectory.kinetic_change_j, -1049.0 / 1024.0, 1e-8);
    close(report.trajectory.external_work_j, 0.0, 1e-8);
    close(report.trajectory.energy_residual_j, 0.0, 1e-8);
    unchanged(&slider, &slider_original);
    unchanged(&incoming, &incoming_original);
}

#[test]
fn cold_force_solve_supports_the_unmodified_original_mass_coupled_post_impact_state() {
    use crate::cohesive_contact::trajectory::sustained::{
        SustainedContact, step_candidate_controlled,
    };
    let elastic = HeldElastic;
    let mass = CoupledMass;
    let law = BilinearCohesive::new(1000.0, 1.0, 0.002, 1000.0, 0.02).unwrap();
    let system = system(&elastic, &law);
    let initial = coupled_motion(&system, &mass, 0.25);
    let x = rows(0, [0, 1, 2]);
    let y = rows(0, [1, 0, 2]);
    let contacts = [contact(401, 0.0, 0.0, &x), contact(402, 0.125, 0.0, &y)];
    // Keep the actual impulse solver's returned state, including every
    // floating-point component. The force solve below receives no warm start.
    let arrived = run(
        &[&initial],
        &contacts,
        &[[-20.0, 0.0, 0.0, HELD_FORCE_N]],
        &[[20.0, 0.0, 0.0], [0.0; 3]],
        0.125,
        0.25,
    );
    assert_eq!(arrived.report.trajectory.events.len(), 1);
    let before = &arrived.motions[0];
    let original = before.clone();
    let held: Vec<_> = contacts
        .iter()
        .zip(&arrived.report.trajectory.final_gaps_m)
        .map(|(contact, &gap_m)| SustainedContact {
            id: contact.id,
            gap_m,
            activation_tolerance_m: contact.activation_tolerance_m,
            friction_coefficient: contact.friction_coefficient,
            prescribed_velocity_m_s: contact.prescribed_velocity_m_s,
            trace: contact.trace,
        })
        .collect();
    let step = step_candidate_controlled(
        &[before],
        &held,
        &[BodyLoad {
            mean_force_n: &[-20.0, 0.0, 0.0, HELD_FORCE_N],
            endpoint_prescribed_jumps_m: &[[0.0; 3]],
        }],
        &[NaturalMapScales {
            normal_n_per_m: 256.0,
            tangent_n_per_m: 256.0,
        }; 2],
        &[[0.0; 3]; 2],
        0.25,
        controls(0.25).sustained.step,
        &mut || false,
    )
    .unwrap();
    let after = &step.bodies[0].motion;
    for (contact, dof) in [0, 1].into_iter().enumerate() {
        let actual_gap_m = held[contact].gap_m + after.state().displacements_m()[dof]
            - before.state().displacements_m()[dof];
        assert!(actual_gap_m >= 0.0, "actual gap {actual_gap_m:e}");
        assert!(actual_gap_m <= 1e-10);
        assert!(step.report.contact.paths_m[contact][0] >= 0.0);
        assert!(step.report.contact.forces_n[contact][0] >= 0.0);
        close(after.velocity_m_s()[dof], 0.0, 1e-8);
    }
    close(step.report.contact.forces_n[0][0], 20.0, 1e-8);
    close(step.report.contact.forces_n[1][0], 0.0, 1e-8);
    assert!(step.report.contact.maximum_residual_n <= 1e-8);
    close(step.report.audit.energy_residual_j, 0.0, 1e-8);
    assert_eq!(after.time_s().to_bits(), (START_S + 0.375).to_bits());
    retained(after, before);
    unchanged(before, &original);
}
