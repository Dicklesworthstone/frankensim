//! Omitted-mode compliance among fixed obstacles (pexfp.299, MagentaBadger
//! 7592 gate 3): the joint solve sees each support in series with it, the
//! ledger stores f G f / 2, and the step-held quasi-static response books its
//! declared numerical loss dR G dR / 2.

use super::*;

// obstacle(): the synthetic support's linear (alpha 1) stiffness [N/m].
const SUPPORT_N_M: f64 = 2000.0;
// 2 K g = 1 for the opposed pair: its net stiffness 2K halves.
const PAIR_G: f64 = 2.5e-4;

fn moving(velocity: f64) -> CoupledModalSystem {
    CoupledModalSystem::new(vec![model(100.0, velocity)], vec![], config(1000.0), &CancelGate::new_clock_free()).unwrap()
}
fn pair(compliance: Option<f64>, velocity: f64) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let contacts = [1.0, -1.0].map(|sign| (obstacle(sign, -1e-4), normal_limits())).to_vec();
    let joint = multi(friction(moving(velocity))).with_obstacles(contacts, limits(), &gate).unwrap();
    match compliance {
        Some(g) => joint.with_obstacle_residual_compliance(vec![g, -g, -g, g], &gate).unwrap(),
        None => joint,
    }
}
fn single(compliance: Option<f64>) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let joint = multi(friction(moving(0.0)))
        .with_obstacles(vec![(obstacle(1.0, -1e-4), normal_limits())], limits(), &gate).unwrap();
    match compliance {
        Some(g) => joint.with_obstacle_residual_compliance(vec![g], &gate).unwrap(),
        None => joint,
    }
}

// Damped angular frequency [rad/s] from the modal velocity's zero crossings.
fn frequency(joint: &mut MultiFrictionModalSystem, steps: u32) -> f64 {
    let dt = joint.sample_period_s();
    let mut crossings = Vec::new();
    let mut last = joint.components()[0].states()[0].velocity_m_sqrt_kg_per_s;
    for step in 1..=steps {
        let frame = joint.step(&[0.0], &[released()]).unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let v = joint.components()[0].states()[0].velocity_m_sqrt_kg_per_s;
        if last != 0.0 && v.signum() != last.signum() { crossings.push(dt*(f64::from(step)-v/(v-last))); }
        last = v;
    }
    assert!(crossings.len() >= 4, "only {} crossings", crossings.len());
    core::f64::consts::PI*(crossings.len()-1) as f64/(crossings[crossings.len()-1]-crossings[0])
}

#[test]
fn supports_act_in_series_with_the_omitted_modes() {
    // Pair: net stiffness 2K/(1 + 2Kg) on a 100 rad/s mode; single: K/(1 + Kg).
    // The 0.02 modal damping moves these by under 2e-4.
    for (mut joint, expected, alone) in [
        (pair(Some(PAIR_G), 0.004), 1e4+2.0*SUPPORT_N_M/(1.0+2.0*SUPPORT_N_M*PAIR_G), false),
        (pair(None, 0.004), 1e4+2.0*SUPPORT_N_M, false),
        (single(Some(2.0*PAIR_G)), 1e4+SUPPORT_N_M/(1.0+2.0*SUPPORT_N_M*PAIR_G), true),
        (single(None), 1e4+SUPPORT_N_M, true),
    ] {
        let paired = joint.obstacle_pairs();
        let bound = joint.obstacle_residual_compliance().is_some();
        let measured = frequency(&mut joint, 12_000);
        assert!((measured/expected.sqrt()-1.0).abs() < 1e-3,
            "single {alone} residual {bound}: {measured} rad/s against {}", expected.sqrt());
        assert_eq!(paired, if bound && !alone { vec![(0, 1)] } else { vec![] });
        for n in joint.last_obstacle_frames().unwrap() { assert!(n.normal_force_n > 0.0, "a support released"); }
    }
}

#[test]
fn the_ledger_needs_the_omitted_modes_storage_and_books_their_numerical_loss() {
    let mut joint = pair(Some(PAIR_G), 0.004);
    let (mut stored, mut net) = (0.0, 0.0);
    let mut unbooked = 0;
    let steps = 4096;
    for _ in 0..steps {
        let f = joint.step(&[0.0], &[released()]).unwrap().clone();
        let n = joint.last_obstacle_frames().unwrap();
        // G = g [[1, -1], [-1, 1]] sees only the net force.
        let now = n[0].normal_force_n-n[1].normal_force_n;
        assert!((f.obstacle_residual_energy_j-0.5*PAIR_G*now*now).abs() <= 1e-18);
        let change = now-net;
        assert!((f.obstacle_residual_loss_j-0.5*PAIR_G*change*change).abs() <= 1e-18);
        assert!(f.obstacle_residual_loss_j >= 0.0);
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        // The same step's identity with f G f / 2 dropped from both ends.
        let dropped = f.energy_residual_j-(f.obstacle_residual_energy_j-stored);
        if dropped.abs() > f.energy_tolerance_j { unbooked += 1; }
        (stored, net) = (f.obstacle_residual_energy_j, now);
    }
    assert!(unbooked > 3*steps/4, "without f G f / 2 the identity still closed on {} of {steps} steps", steps-unbooked);
}

#[test]
fn residual_compliance_refuses_what_it_does_not_own() {
    let gate = CancelGate::new_clock_free();
    assert!(multi(friction(moving(0.0))).with_obstacle_residual_compliance(vec![1e-4], &gate).is_err());
    let pinch = || vec![PAIR_G, -PAIR_G, -PAIR_G, PAIR_G];
    for bad in [vec![1e-4], vec![1e-4, 0.0, 0.0], vec![-1e-4, 0.0, 0.0, 1e-4], vec![1e-4, 2e-4, 2e-4, 1e-4],
        vec![1e-4, 1e-5, 1.000_000_1e-5, 1e-4], vec![f64::NAN, 0.0, 0.0, 1e-4]] {
        assert!(pair(None, 0.0).with_obstacle_residual_compliance(bad, &gate).is_err());
    }
    assert!(pair(Some(PAIR_G), 0.0).with_obstacle_residual_compliance(pinch(), &gate).is_err());
    let mut played = pair(None, 0.0);
    played.step(&[0.0], &[released()]).unwrap();
    assert!(played.with_obstacle_residual_compliance(pinch(), &gate).is_err());
    let error = pair(Some(PAIR_G), 0.0).initialize_released_equilibrium(&[0.0], &gate).unwrap_err();
    assert!(format!("{error:?}").contains("residual compliance"), "{error:?}");
}
