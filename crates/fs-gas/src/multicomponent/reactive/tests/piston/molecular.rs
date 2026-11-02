use super::*;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationOperator, TransportCoefficients,
    TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

pub(super) fn molecular(flow: &MixtureEuler1d) -> DissipationOperator<ConstantTransport> {
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "moving-molecular-piston-fixture",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.02),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.04),
        },
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    DissipationOperator::new(flow, law, || false).unwrap()
}

fn nonuniform(reference: f64) -> ReactingEuler1d {
    let mechanism = mechanism(8.0, 30.0, reference);
    let mut flow = field(&mechanism, 12, 0.0, |x| {
        let a = 0.7 + 0.1 * (2.0 * core::f64::consts::PI * x).sin();
        [a, 1.0 - a]
    });
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap()
}

fn step(
    flow: &mut ReactingEuler1d,
    operator: &mut DissipationOperator<ConstantTransport>,
    piston: &mut Piston,
    cancelled: impl FnMut() -> bool,
) -> Result<crate::dynamics::piston::PistonStep<ReactiveStepReport>, ReactiveError> {
    flow.step_with_piston_and_dissipation(
        operator,
        piston,
        Time::new(0.005),
        ReactiveStepControl::default(),
        DissipationControl::default(),
        load(),
        piston_controls(),
        cancelled,
    )
}

#[test]
fn g1_moving_molecular_wall_matches_independent_rk2_stress_and_work() {
    let mechanism = mechanism(0.0, 0.0, 0.0);
    let mut flow = field(&mechanism, 1, 0.2, |_| [0.6, 0.4]);
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    let mut operator = molecular(&flow);
    operator.admit_motion(&flow, &mut || false).unwrap();
    let initial = flow.totals();
    let dt = 0.01;
    let wall_speed = 0.7;
    // L=rho=1. Each wall is L/2 from the cell center. Newtonian longitudinal
    // viscosity eta=4*mu/3 gives du/dt=2*eta*(w-2*u). Independently integrate
    // this scalar equation with the two SSP stages, including each wall force.
    let eta = (4.0 / 3.0) * 0.02;
    let u0 = 0.2;
    let u_stage = u0 + dt * 2.0 * eta * (wall_speed - 2.0 * u0);
    let u1 = 0.5 * u0 + 0.5 * (u_stage + dt * 2.0 * eta * (wall_speed - 2.0 * u_stage));
    let right_impulse = dt * AREA * 2.0 * eta * (wall_speed - 0.5 * (u0 + u_stage));
    let left_impulse = -dt * AREA * 2.0 * eta * (0.5 * (u0 + u_stage));
    let report = operator
        .advance_moving_trial(
            &mut flow,
            Time::new(dt),
            DissipationControl::default(),
            [0.0, wall_speed],
            &mut || false,
        )
        .unwrap();
    assert_eq!(report.accepted_steps, 1);
    assert_eq!(report.rejected_steps, 0);
    close(
        flow.gas()
            .primitive(flow.cell(0).unwrap(), || false)
            .unwrap()
            .velocity,
        u1,
        2e-14,
    );
    close(report.right_boundary_change.momentum, right_impulse, 2e-15);
    close(report.left_boundary_change.momentum, left_impulse, 2e-15);
    close(
        report.right_boundary_change.energy,
        wall_speed * right_impulse,
        2e-15,
    );
    assert_eq!(report.left_boundary_change.energy, 0.0);
    close(
        flow.totals().energy - initial.energy,
        wall_speed * right_impulse,
        2e-14,
    );
    close(
        flow.totals().momentum - initial.momentum,
        right_impulse + left_impulse,
        2e-14,
    );
    assert!(report.right_boundary_change.energy > 0.0);
    assert_eq!(
        flow.time().value(),
        0.0,
        "the split molecular stage owns no clock motion"
    );
}

#[test]
fn g1_reacting_viscous_piston_closes_physical_balances_and_energy_gauge() {
    let mut a = nonuniform(0.0);
    let mut b = nonuniform(1000.0);
    let mut da = molecular(a.flow());
    let mut db = molecular(b.flow());
    let mut pa = Piston::new(0.2, 0.0).unwrap();
    let mut pb = pa;
    let before = a.flow().totals();
    let mut exterior_work = 0.0;
    let mut exterior_impulse = 0.0;
    let mut support = 0.0;
    let mut viscous_impulse = 0.0;
    for _ in 0..6 {
        let report = step(&mut a, &mut da, &mut pa, || false).unwrap();
        let shifted = step(&mut b, &mut db, &mut pb, || false).unwrap();
        let halves = report.gas.dissipation.unwrap();
        assert_eq!(
            halves.first.elapsed.value(),
            0.5 * report.gas.elapsed.value()
        );
        assert_eq!(halves.second.elapsed, halves.first.elapsed);
        assert_eq!(da.mesh(), a.flow().mesh());
        assert_eq!(db.mesh(), b.flow().mesh());
        close(
            report.total_energy_residual_j,
            0.0,
            piston_controls().energy_tolerance_j,
        );
        close(
            report.combined_momentum_residual_n_s,
            0.0,
            piston_controls().momentum_tolerance_n_s,
        );
        for half in [halves.first, halves.second] {
            assert_eq!(half.boundary_change.species_masses[..2], [0.0; 2]);
            assert_eq!(half.left_boundary_change.energy, 0.0);
            close(
                half.right_boundary_change.energy,
                report.midpoint_velocity_m_s * half.right_boundary_change.momentum,
                3e-14,
            );
            close(half.absolute_energy_residual_j, 0.0, 2e-12);
            viscous_impulse += half.right_boundary_change.momentum.abs();
        }
        close(
            report.gas_impulse_n_s,
            -report.gas.transport.right_boundary_change.momentum
                - halves.first.right_boundary_change.momentum
                - halves.second.right_boundary_change.momentum,
            3e-14,
        );
        close(report.gas.mass_residual_kg, 0.0, 3e-12);
        exterior_work += report.exterior_pressure_work_j;
        exterior_impulse -= AREA * report.gas.elapsed.value();
        support += report.fixed_wall_impulse_n_s;
        close(
            report.gas.elapsed.value(),
            shifted.gas.elapsed.value(),
            1e-15,
        );
        close(pa.velocity_m_s(), pb.velocity_m_s(), 3e-9);
        for i in 0..a.flow().cell_count() {
            let qa = a.flow().cell(i).unwrap();
            let qb = b.flow().cell(i).unwrap();
            for k in 0..2 {
                close(qa.partial_densities[k], qb.partial_densities[k], 3e-9);
            }
            let ua = a.flow().gas().primitive(qa, || false).unwrap();
            let ub = b.flow().gas().primitive(qb, || false).unwrap();
            close(ua.pressure, ub.pressure, 3e-9);
            close(ua.temperature, ub.temperature, 3e-9);
        }
    }
    let after = a.flow().totals();
    let physical_energy = after.energy - before.energy
        + (0..2)
            .map(|k| {
                (after.species_masses[k] - before.species_masses[k])
                    * a.flow().gas().energy_references()[k]
            })
            .sum::<f64>();
    close(
        physical_energy + 0.5 * pa.mass_kg() * pa.velocity_m_s().powi(2) + exterior_work,
        0.0,
        1e-7,
    );
    close(
        after.momentum - before.momentum + pa.mass_kg() * pa.velocity_m_s() + support
            - exterior_impulse,
        0.0,
        3e-10,
    );
    assert!(
        viscous_impulse > 1e-9,
        "the coupled root must include a nonzero molecular load"
    );
    assert!(
        pa.velocity_m_s() > 0.0,
        "reaction pressure must expand the loaded chamber"
    );
    assert!(a.flow().mesh().positions()[a.flow().cell_count()] > 1.0);
}

#[test]
fn g1_shortened_piston_replays_all_sources_at_one_speed_and_duration() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut flow = closed(&mechanism, 4, chemistry());
    let mut replay = closed(&mechanism, 4, chemistry());
    let mut operator = molecular(flow.flow());
    let mut replay_operator = molecular(replay.flow());
    let mut piston = Piston::new(1.0, 0.0).unwrap();
    let report = flow
        .step_with_piston_and_dissipation(
            &mut operator,
            &mut piston,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            load(),
            piston_controls(),
            || false,
        )
        .unwrap();
    assert!(report.gas.coupling_retries > 0);
    let dt = report.gas.elapsed.value();
    let speed = report.midpoint_velocity_m_s;
    replay_operator
        .admit_motion(&replay.accepted, &mut || false)
        .unwrap();
    let first = replay_operator
        .advance_moving_trial(
            &mut replay.accepted,
            Time::new(0.5 * dt),
            DissipationControl::default(),
            [0.0, speed],
            &mut || false,
        )
        .unwrap();
    let mut chemical = ChemistryReport::default();
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut chemical,
        &mut || false,
    )
    .unwrap();
    let transport = replay
        .accepted
        .step_with_moving_walls(
            Time::new(dt),
            ReactiveStepControl::default().transport,
            [0.0, speed],
            || false,
        )
        .unwrap();
    assert_eq!(transport.elapsed.value(), dt);
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut chemical,
        &mut || false,
    )
    .unwrap();
    let second = replay_operator
        .advance_moving_trial(
            &mut replay.accepted,
            Time::new(0.5 * dt),
            DissipationControl::default(),
            [0.0, speed],
            &mut || false,
        )
        .unwrap();
    assert_eq!(flow.flow().mesh(), replay.flow().mesh());
    assert_eq!(flow.flow().time(), replay.flow().time());
    for (actual, expected) in flow.accepted.states.iter().zip(&replay.accepted.states) {
        close(*actual, *expected, 3e-12);
    }
    let halves = report.gas.dissipation.unwrap();
    close(
        halves.first.right_boundary_change.energy,
        first.right_boundary_change.energy,
        3e-13,
    );
    close(
        halves.second.right_boundary_change.energy,
        second.right_boundary_change.energy,
        3e-13,
    );
    assert_ne!(
        first.before, second.before,
        "the final molecular stage uses the evolved state"
    );
}

#[test]
fn g4_late_cancellation_preserves_moved_bindings_and_allows_stationary_continuation() {
    let mut flow = nonuniform(0.0);
    let mut reference = nonuniform(0.0);
    let mut operator = molecular(flow.flow());
    let mut reference_operator = molecular(reference.flow());
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let mut reference_piston = piston;
    step(&mut flow, &mut operator, &mut piston, || false).unwrap();
    step(
        &mut reference,
        &mut reference_operator,
        &mut reference_piston,
        || false,
    )
    .unwrap();
    let original_states = flow.accepted.states.clone();
    let original_mesh = flow.flow().mesh().clone();
    let original_time = flow.flow().time();
    let original_piston = piston;
    let mut successful_polls = 0;
    let expected = step(
        &mut reference,
        &mut reference_operator,
        &mut reference_piston,
        || {
            successful_polls += 1;
            false
        },
    )
    .unwrap();
    assert!(successful_polls > 50);
    let mut polls = 0;
    let error = step(&mut flow, &mut operator, &mut piston, || {
        polls += 1;
        polls == successful_polls
    })
    .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(polls, successful_polls);
    assert_eq!(flow.accepted.states, original_states);
    assert_eq!(flow.flow().mesh(), &original_mesh);
    assert_eq!(operator.mesh(), &original_mesh);
    assert_eq!(flow.flow().time(), original_time);
    assert_eq!(piston, original_piston);
    let resumed = step(&mut flow, &mut operator, &mut piston, || false).unwrap();
    assert_eq!(resumed, expected);
    assert_eq!(flow.accepted.states, reference.accepted.states);
    assert_eq!(operator.mesh(), reference_operator.mesh());
    assert_eq!(piston, reference_piston);
    // Stopping prescribed motion does not invalidate the accepted source or
    // geometry. Scratch from the refused moving trial must not leak widths.
    let stationary = flow
        .step_with_dissipation(
            &mut operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let stationary_reference = reference
        .step_with_dissipation(
            &mut reference_operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(stationary, stationary_reference);
    assert_eq!(flow.accepted.states, reference.accepted.states);
    assert_eq!(operator.mesh(), flow.flow().mesh());
}
