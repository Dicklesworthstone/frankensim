use super::*;
use fs_phs::linear_exact::LinearSplitBudget;
use fs_solid::membrane::viscosity::modal::{
    MembraneDampingTimeError, PreparedMembraneDamping,
};

fn flow_budget() -> LinearSplitBudget {
    LinearSplitBudget { max_modes: 16, max_omega_dt: 1.0 }
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-13 + 3e-10 * a.abs().max(b.abs()), "{a:e} != {b:e}");
}

// Independent RK4 reference in the original nodal constitutive owner. It never
// reads the projected C used by the prepared exponential. The last two state
// entries integrate physical heat and held-force work as separate ODEs.
fn reference_derivative(
    damping: &MembraneViscosity<'_>,
    basis: &MembraneModes<'_>,
    state: &[f64],
    force: &[f64],
) -> Vec<f64> {
    let n = basis.modes().len();
    let velocity = basis.reconstruct(&state[n..2 * n]).unwrap();
    let physical = damping.evaluate(
        basis.linearization_displacement_m(), &velocity, || false,
    ).unwrap();
    let resistance = basis.project_nodal_force(&physical.resisting_force_n).unwrap();
    let mut result = vec![0.0; 2 * n + 2];
    for i in 0..n {
        result[i] = state[n + i];
        result[n + i] = force[i] - basis.modes()[i].lambda * state[i] - resistance[i];
    }
    result[2 * n] = physical.power_w;
    result[2 * n + 1] = (0..n).map(|i| force[i] * state[n + i]).sum();
    result
}

fn rk4(mut state: Vec<f64>, duration: f64, steps: usize, rhs: impl Fn(&[f64]) -> Vec<f64>)
    -> Vec<f64>
{
    let h = duration / steps as f64;
    for _ in 0..steps {
        let a = rhs(&state);
        let stage = |slope: &[f64], scale: f64| state.iter().zip(slope)
            .map(|(x, slope)| x + scale * slope).collect::<Vec<_>>();
        let b = rhs(&stage(&a, 0.5 * h));
        let c = rhs(&stage(&b, 0.5 * h));
        let d = rhs(&stage(&c, h));
        for i in 0..state.len() {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn exact_coupled_motion_matches_nodal_material_loss_and_held_force_work() {
    let (p, laws) = patch([40_000.0, 30_000.0, 20_000.0, 10_000.0, 4_000.0, 2_000.0], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let basis = modes(&p, &deformed(), &[5, 8]);
    let dt = 0.07;
    let mut prepared = PreparedMembraneDamping::new(
        &damping, &basis, dt, flow_budget(), 100_000, || false,
    ).unwrap();
    assert_eq!(prepared.mode_count(), 2);
    assert_eq!(prepared.state_dimension(), 4);
    assert!(prepared.projected_viscosity().matrix_per_s()[1].abs() > 0.1);
    let initial = [0.03, -0.02, 0.4, -0.2];
    let force = [0.7, -0.4];
    let mut endpoint = [0.0; 4];
    let mut average_velocity = [0.0; 2];
    let record = prepared.step_into(
        &initial, &force, &mut endpoint, &mut average_velocity, || false,
    ).unwrap();
    let reference = rk4(
        vec![initial[0], initial[1], initial[2], initial[3], 0.0, 0.0],
        dt, 2048, |state| reference_derivative(&damping, &basis, state, &force),
    );
    for i in 0..4 { near(endpoint[i], reference[i]); }
    near(record.dissipated, reference[4]);
    near(record.supplied, reference[5]);
    assert!(record.dissipated > 0.0);
    near(record.delta_h + record.dissipated, record.supplied);
    for i in 0..2 {
        near(dt * average_velocity[i], endpoint[i] - initial[i]);
    }

    // Causal negative control: removing off-diagonal material loss changes the
    // actual trajectory, even when the same modal masses/frequencies are kept.
    let c = prepared.projected_viscosity().matrix_per_s();
    let diagonal = rk4(initial.to_vec(), dt, 2048, |state| vec![
        state[2], state[3],
        force[0] - basis.modes()[0].lambda * state[0] - c[0] * state[2],
        force[1] - basis.modes()[1].lambda * state[1] - c[3] * state[3],
    ]);
    assert!(endpoint.iter().zip(diagonal).map(|(a, b)| (a - b).abs()).sum::<f64>() > 1e-4);
}

#[test]
fn scalar_membrane_decay_matches_closed_form_heat_and_subdivision() {
    let (p, laws) = patch([20_000.0, 15_000.0, 10_000.0, 5_000.0, 2_000.0, 1_000.0], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let basis = modes(&p, &[[0.0; 3]; 3], &[3]);
    let dt = 0.04;
    let mut full = PreparedMembraneDamping::new(
        &damping, &basis, dt, flow_budget(), 100_000, || false,
    ).unwrap();
    let mut half = PreparedMembraneDamping::new(
        &damping, &basis, dt / 2.0, flow_budget(), 100_000, || false,
    ).unwrap();
    let lambda = basis.modes()[0].lambda;
    let c = full.projected_viscosity().matrix_per_s()[0];
    // Reference triangle A=.5,h=.002,m_node=1: K=600 and C=20.
    close(lambda, 600.0);
    close(c, 20.0);
    let sigma = c / 2.0;
    let wd = (lambda - sigma * sigma).sqrt();
    let initial = [0.03, -0.2];
    let force = [0.7];
    let shifted = initial[0] - force[0] / lambda;
    let velocity = |t: f64| (-sigma * t).exp() * (
        initial[1] * (wd * t).cos()
        - (sigma * initial[1] + lambda * shifted) * (wd * t).sin() / wd
    );
    let expected_q = force[0] / lambda + (-sigma * dt).exp() * (
        shifted * (wd * dt).cos()
        + (initial[1] + sigma * shifted) * (wd * dt).sin() / wd
    );
    let mut endpoint = [0.0; 2];
    let mut average = [0.0; 1];
    let record = full.step_into(&initial, &force, &mut endpoint, &mut average, || false).unwrap();
    near(endpoint[0], expected_q);
    near(endpoint[1], velocity(dt));
    // Composite Simpson integration of the independent analytic velocity.
    let intervals = 2048;
    let h = dt / intervals as f64;
    let mut integral = velocity(0.0).powi(2) + velocity(dt).powi(2);
    for i in 1..intervals {
        integral += (if i % 2 == 0 { 2.0 } else { 4.0 }) * velocity(i as f64 * h).powi(2);
    }
    near(record.dissipated, c * h * integral / 3.0);

    let mut middle = [0.0; 2];
    let mut split = [0.0; 2];
    let a = half.step_into(&initial, &force, &mut middle, &mut average, || false).unwrap();
    let b = half.step_into(&middle, &force, &mut split, &mut average, || false).unwrap();
    near(split[0], endpoint[0]);
    near(split[1], endpoint[1]);
    near(a.dissipated + b.dissipated, record.dissipated);
    near(a.supplied + b.supplied, record.supplied);
}

#[test]
fn flat_transverse_modes_keep_their_physical_zero_loss() {
    let (p, laws) = patch([4.0, 3.0, 2.0, 1.0, 0.4, 0.2], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let stretch = [[0.0; 3], [0.02, 0.0, 0.0], [0.0, 0.02, 0.0]];
    let basis = modes(&p, &stretch, &[5]);
    let dt = 0.17;
    let mut prepared = PreparedMembraneDamping::new(
        &damping, &basis, dt, flow_budget(), 100_000, || false,
    ).unwrap();
    let initial = [0.003, 0.04];
    let mut next = [0.0; 2];
    let mut output = [0.0; 1];
    let record = prepared.step_into(&initial, &[0.0], &mut next, &mut output, || false).unwrap();
    let omega = basis.modes()[0].lambda.sqrt();
    near(next[0], initial[0] * (omega * dt).cos() + initial[1] * (omega * dt).sin() / omega);
    near(next[1], initial[1] * (omega * dt).cos() - omega * initial[0] * (omega * dt).sin());
    assert_eq!(record.dissipated, 0.0);
    near(record.delta_h, 0.0);
}

#[test]
fn refused_or_late_cancelled_steps_publish_no_partial_membrane_state() {
    let (p, laws) = patch([4.0, 3.0, 2.0, 1.0, 0.4, 0.2], 0.002, 0.0);
    let damping = MembraneViscosity::new(&p, laws, (0.0, 1000.0), || false).unwrap();
    let basis = modes(&p, &[[0.0; 3]; 3], &[3]);
    let mut prepared = PreparedMembraneDamping::new(
        &damping, &basis, 0.04, flow_budget(), 100_000, || false,
    ).unwrap();
    let initial = [0.003, 0.04];
    let mut next = [91.0; 2];
    let mut output = [73.0; 1];
    let mut polls = 0;
    let refused = prepared.step_into(&initial, &[0.0], &mut next, &mut output, || {
        polls += 1;
        polls == 2
    });
    assert!(matches!(refused, Err(MembraneDampingTimeError::Membrane(MembraneError::Cancelled))));
    assert_eq!(next, [91.0; 2]);
    assert_eq!(output, [73.0; 1]);
    assert!(prepared.step_into(&[f64::NAN, 0.0], &[0.0], &mut next, &mut output, || false).is_err());
    assert_eq!(next, [91.0; 2]);
    assert_eq!(output, [73.0; 1]);
    assert!(PreparedMembraneDamping::new(
        &damping, &basis, 0.04,
        LinearSplitBudget { max_modes: 1, max_omega_dt: 1.0 }, 100_000, || false,
    ).is_err());
    prepared.step_into(&initial, &[0.0], &mut next, &mut output, || false).unwrap();
    assert!(next.iter().all(|x| x.is_finite()));
}
