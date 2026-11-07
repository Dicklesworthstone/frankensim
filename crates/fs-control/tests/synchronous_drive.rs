//! Sampled current/torque regulation against the original saturated native
//! finite-bus plant. These are authored numerical scenarios, not a universal
//! closed-loop stability, transient-current or experimentally validated claim.
use fs_control::synchronous::{
    CurrentControllerParameters, CurrentControllerState, SynchronousCurrentController,
    SynchronousMeasurement,
};
use fs_phs::circuit::inductance::InductorFluxCurve;
use fs_phs::circuit::synchronous::{
    SynchronousMotor, SynchronousMotorParameters,
    inverter::{AveragedSynchronousInverter, InverterParameters, PreparedSynchronousInverter},
};
use fs_phs::{PreparedStepError, SolveBudget};

const TS: f64 = 1e-4;
const CURRENT_LIMIT: f64 = 8.0;
const D_CURRENT: f64 = -4.0;
const D: [(f64, f64); 7] = [
    (-0.22, -50.0),
    (-0.10, -10.0),
    (-0.04, -2.0),
    (0.0, 0.0),
    (0.04, 2.0),
    (0.10, 10.0),
    (0.22, 50.0),
];
const Q: [(f64, f64); 7] = [
    (-0.30, -60.0),
    (-0.12, -12.0),
    (-0.04, -2.0),
    (0.0, 0.0),
    (0.04, 2.0),
    (0.12, 12.0),
    (0.30, 60.0),
];

fn machine() -> SynchronousMotorParameters {
    SynchronousMotorParameters {
        d_inductance_h: 0.02,
        q_inductance_h: 0.02,
        winding_resistance_ohm: 0.35,
        pm_flux_linkage_wb: 0.08,
        pole_pairs: 3,
        rotor_inertia_kg_m2: 0.025,
        stator_inertia_kg_m2: 0.07,
        bearing_drag_nm_s_per_rad: 0.002,
    }
}

fn bridge() -> InverterParameters {
    InverterParameters {
        dc_capacitance_f: 0.3,
        high_on_resistance_ohm: [0.03, 0.04, 0.02],
        low_on_resistance_ohm: [0.02, 0.03, 0.04],
    }
}

fn fixture() -> (
    PreparedSynchronousInverter,
    SynchronousCurrentController,
    [f64; 6],
) {
    let p = machine();
    let motor = SynchronousMotor::with_flux_curves(
        p,
        InductorFluxCurve::try_new(D.to_vec()).unwrap(),
        InductorFluxCurve::try_new(Q.to_vec()).unwrap(),
    )
    .unwrap();
    let state = motor.state([D_CURRENT, 0.0], 0.1, [20.0, -3.0]).unwrap();
    let inverter = AveragedSynchronousInverter::new(motor, bridge()).unwrap();
    let initial = inverter.state(24.0, &state).unwrap();
    let plant = inverter
        .prepare_controlled(
            SolveBudget {
                relative: 2e-11,
                absolute: 2e-12,
                max_iterations: 30,
            },
            8,
            || false,
        )
        .unwrap();
    let controller = SynchronousCurrentController::new(CurrentControllerParameters {
        sample_period_s: TS,
        pole_pairs: p.pole_pairs,
        proportional_gain_ohm: [2.0; 2],
        integral_gain_ohm_per_s: [200.0; 2],
        anti_windup_time_s: [0.002; 2],
        winding_resistance_ohm: p.winding_resistance_ohm,
        current_limit_a: CURRENT_LIMIT,
    })
    .unwrap();
    (plant, controller, initial)
}

fn measurement(plant: &PreparedSynchronousInverter, x: &[f64; 6]) -> SynchronousMeasurement {
    let point = plant.machine_operating_point(x).unwrap();
    SynchronousMeasurement {
        current_alpha_beta_a: point.current_alpha_beta_a,
        flux_alpha_beta_wb: [x[1], x[2]],
        relative_mechanical_angle_rad: x[3],
        shaft_speed_rad_s: point.shaft_speed_rad_s,
        dc_bus_voltage_v: x[0] / bridge().dc_capacitance_f,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

// Independently integrate/interpolate the authored physical source segments.
// No native storage, torque, DG or energy callback defines these oracles.
fn source_current(knots: &[(f64, f64)], flux: f64) -> f64 {
    let pair = knots
        .windows(2)
        .find(|k| k[0].0 <= flux && flux <= k[1].0)
        .unwrap();
    pair[0].1 + (flux - pair[0].0) * (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0)
}

fn source_energy(knots: &[(f64, f64)], flux: f64) -> f64 {
    let (lo, hi, sign) = if flux < 0.0 {
        (flux, 0.0, -1.0)
    } else {
        (0.0, flux, 1.0)
    };
    assert!(lo >= knots[0].0 && hi <= knots[knots.len() - 1].0);
    let mut energy = 0.0;
    for pair in knots.windows(2) {
        let a = lo.max(pair[0].0);
        let b = hi.min(pair[1].0);
        if b > a {
            energy += 0.5 * (source_current(knots, a) + source_current(knots, b)) * (b - a);
        }
    }
    sign * energy
}

fn physical_flux(x: &[f64; 6]) -> [f64; 2] {
    let (s, c) = (f64::from(machine().pole_pairs) * x[3]).sin_cos();
    [c * x[1] + s * x[2], -s * x[1] + c * x[2]]
}

fn energy(x: &[f64; 6]) -> f64 {
    let p = machine();
    let [d, q] = physical_flux(x);
    0.5 * x[0] * x[0] / bridge().dc_capacitance_f
        + source_energy(&D, d - p.pm_flux_linkage_wb)
        + source_energy(&Q, q)
        + 0.5 * (x[4] * x[4] / p.rotor_inertia_kg_m2 + x[5] * x[5] / p.stator_inertia_kg_m2)
}

fn torque(x: &[f64; 6]) -> f64 {
    let p = machine();
    let [d, q] = physical_flux(x);
    let id = source_current(&D, d - p.pm_flux_linkage_wb);
    let iq = source_current(&Q, q);
    f64::from(p.pole_pairs) * (d * iq - q * id)
}

#[test]
fn saturated_closed_loop_tracks_torque_and_regenerates_with_both_shaft_reactions() {
    let (mut plant, controller, mut state) = fixture();
    let initial = state;
    let mut accepted_pi = CurrentControllerState::default();
    let mut total_heat = 0.0;
    let mut total_em_impulse = 0.0;
    let mut total_bearing_impulse = 0.0;
    let mut limited_samples = 0;
    let mut duty_changes = 0;
    let mut prior_duty = plant.high_fraction();
    let mut tail_error = [0.0_f64; 3];
    let mut stage_bridge_work = [0.0; 3];
    let mut stage_bus_rise = [0.0; 3];
    for (stage, requested_torque) in [1.2, -1.2, 0.0].into_iter().enumerate() {
        let before_bus = state[0] / bridge().dc_capacitance_f;
        for sample_index in 0..500 {
            let reference = plant
                .current_reference_for_torque_controlled(
                    requested_torque,
                    D_CURRENT,
                    CURRENT_LIMIT,
                    || false,
                )
                .unwrap();
            assert!(!reference.torque_limited);
            close(reference.achievable_torque_nm, requested_torque, 1e-12);
            let proposal = controller
                .propose_controlled(
                    &accepted_pi,
                    &measurement(&plant, &state),
                    reference.current_dq_a,
                    || false,
                )
                .unwrap();
            assert!(!proposal.current_limited);
            assert!(
                proposal
                    .duty_high_fraction
                    .iter()
                    .all(|d| (0.0..=1.0).contains(d))
            );
            limited_samples += usize::from(proposal.voltage_limited);
            duty_changes += usize::from(proposal.duty_high_fraction != prior_duty);
            let (mut next, mut outputs) = ([0.0; 6], [0.0; 3]);
            let step = plant
                .step_controlled(
                    &state,
                    proposal.duty_high_fraction,
                    &[0.0; 3],
                    TS,
                    &mut next,
                    &mut outputs,
                    || false,
                )
                .unwrap();
            let channel_heat: f64 = step.channel_heat_j.iter().flatten().sum();
            assert!(step.copper_heat_j >= 0.0 && step.bearing_heat_j >= 0.0);
            assert!(step.channel_heat_j.iter().flatten().all(|h| *h >= 0.0));
            let heat = step.copper_heat_j + step.bearing_heat_j + channel_heat;
            close(heat, step.native.dissipated, 3e-14);
            close(
                step.bridge_work_j,
                step.machine_electrical_work_j + channel_heat,
                3e-14,
            );
            close(next[0] - state[0], -TS * step.dc_bridge_current_a, 5e-11);
            close(
                next[4] - state[4],
                step.electromagnetic_angular_impulse_nm_s + step.bearing_angular_impulse_nm_s,
                5e-12,
            );
            close(next[4] + next[5], initial[4] + initial[5], 5e-11);
            total_heat += heat;
            total_em_impulse += step.electromagnetic_angular_impulse_nm_s;
            total_bearing_impulse += step.bearing_angular_impulse_nm_s;
            stage_bridge_work[stage] += step.bridge_work_j;
            // The proposed controller state is published only AFTER the full
            // original finite-bus/magnetic/mechanical interval was accepted.
            accepted_pi = proposal.next_state;
            state = next;
            prior_duty = proposal.duty_high_fraction;
            let point = plant.machine_operating_point(&state).unwrap();
            close(point.electromagnetic_torque_nm, torque(&state), 1e-12);
            assert!(point.differential_inductance_dq_h[0] < 0.02);
            if sample_index >= 450 {
                tail_error[stage] =
                    tail_error[stage].max((torque(&state) - requested_torque).abs());
                if stage < 2 {
                    assert!(point.differential_inductance_dq_h[1] < 0.02);
                }
            }
        }
        stage_bus_rise[stage] = state[0] / bridge().dc_capacitance_f - before_bus;
    }
    assert!(
        tail_error.iter().all(|e| *e < 0.05),
        "tail torque error [N m]: {tail_error:?}"
    );
    assert!(
        limited_samples > 0,
        "reversal must exercise physical voltage saturation"
    );
    assert!(
        duty_changes > 1400,
        "rotating sample/hold commands must actually refresh modulation"
    );
    assert!(
        stage_bridge_work[0] > 0.5,
        "motoring bus work: {stage_bridge_work:?}"
    );
    assert!(
        stage_bridge_work[1] < -0.1,
        "negative torque must regenerate net bus energy: {stage_bridge_work:?}"
    );
    assert!(
        stage_bus_rise[1] > 0.005,
        "finite bus must gain voltage during regeneration: {stage_bus_rise:?}"
    );
    assert!(total_heat > 1.0);
    close(energy(&state) - energy(&initial), -total_heat, 2e-7);
    close(
        state[4] - initial[4],
        total_em_impulse + total_bearing_impulse,
        2e-10,
    );
    close(
        state[5] - initial[5],
        -total_em_impulse - total_bearing_impulse,
        2e-10,
    );
}

#[test]
fn rejected_or_cancelled_plant_interval_does_not_accept_pi_state_or_duties() {
    let (mut plant, controller, mut accepted_physical) = fixture();
    let physical_before = accepted_physical;
    let mut accepted_pi = CurrentControllerState::default();
    let pi_before = accepted_pi;
    let held_before = plant.high_fraction();
    let reference = plant
        .current_reference_for_torque_controlled(1.2, D_CURRENT, CURRENT_LIMIT, || false)
        .unwrap();
    let proposal = controller
        .propose_controlled(
            &accepted_pi,
            &measurement(&plant, &accepted_physical),
            reference.current_dq_a,
            || false,
        )
        .unwrap();
    assert_ne!(proposal.next_state, accepted_pi);
    assert_ne!(proposal.duty_high_fraction, held_before);
    let (mut candidate, mut outputs) = ([123.0; 6], [456.0; 3]);
    // Removing twice the stored capacitor charge forces a real physical bus
    // refusal. The controller proposal must not become accepted history.
    let impossible_source = -2.0 * accepted_physical[0] / TS;
    assert!(
        plant
            .step_controlled(
                &accepted_physical,
                proposal.duty_high_fraction,
                &[impossible_source, 0.0, 0.0],
                TS,
                &mut candidate,
                &mut outputs,
                || false,
            )
            .is_err()
    );
    assert_eq!(candidate, [123.0; 6]);
    assert_eq!(outputs, [456.0; 3]);
    assert_eq!(accepted_physical, physical_before);
    assert_eq!(accepted_pi, pi_before);
    assert_eq!(plant.high_fraction(), held_before);
    // Polls 1/2 precede the changed duty stamp. The third reaches the native
    // solve after J/R/witness replacement and exercises actual stamp rollback.
    let mut polls = 0;
    assert!(matches!(
        plant.step_controlled(
            &accepted_physical,
            proposal.duty_high_fraction,
            &[0.0; 3],
            TS,
            &mut candidate,
            &mut outputs,
            || {
                polls += 1;
                polls == 3
            },
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(polls, 3);
    assert_eq!(candidate, [123.0; 6]);
    assert_eq!(outputs, [456.0; 3]);
    assert_eq!(accepted_physical, physical_before);
    assert_eq!(accepted_pi, pi_before);
    assert_eq!(plant.high_fraction(), held_before);
    let retry_proposal = controller
        .propose_controlled(
            &accepted_pi,
            &measurement(&plant, &accepted_physical),
            reference.current_dq_a,
            || false,
        )
        .unwrap();
    assert_eq!(proposal, retry_proposal);
    let retried = plant
        .step_controlled(
            &accepted_physical,
            retry_proposal.duty_high_fraction,
            &[0.0; 3],
            TS,
            &mut candidate,
            &mut outputs,
            || false,
        )
        .unwrap();
    let (mut fresh, _, fresh_initial) = fixture();
    let (mut expected, mut expected_outputs) = ([0.0; 6], [0.0; 3]);
    let direct = fresh
        .step_controlled(
            &fresh_initial,
            proposal.duty_high_fraction,
            &[0.0; 3],
            TS,
            &mut expected,
            &mut expected_outputs,
            || false,
        )
        .unwrap();
    assert_eq!(candidate, expected);
    assert_eq!(outputs, expected_outputs);
    assert_eq!(retried.native, direct.native);
    accepted_physical = candidate;
    accepted_pi = retry_proposal.next_state;
    assert_ne!(accepted_physical, physical_before);
    assert_ne!(accepted_pi, pi_before);
    assert_eq!(plant.high_fraction(), proposal.duty_high_fraction);
}
