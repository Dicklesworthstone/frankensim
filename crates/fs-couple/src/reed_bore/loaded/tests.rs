//! G0/G1: the complete load sees the same interval as mechanics and lay work.
use super::*;

fn reed() -> BeatingReed {
    BeatingReed {
        rest_opening_m: 4e-4,
        width_m: 0.013,
        closing_pressure_pa: 6000.0,
        blowing_pressure_pa: 300.0,
        attack_s: 0.0,
        mass_kg: 1e-5,
        stiffness_n_m: 500.0,
        damping_ratio: 0.04,
    }
}

fn assert_close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * scale.max(f64::MIN_POSITIVE),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn held_load_accounts_for_contact_entry_release_and_original_viscous_heat() {
    let reed = reed();
    let mass = reed.mass_kg;
    let stiffness = reed.stiffness_n_m;
    let damping = 2.0 * reed.damping_ratio * (stiffness * mass).sqrt();
    let face = stiffness * reed.rest_opening_m / reed.closing_pressure_pa;
    let contact_stiffness = 1e8;
    let potential = |y: f64| contact_stiffness * (-y).max(0.0).powi(3) / 3.0;
    let energy = |y: f64, v: f64| {
        0.5 * mass * v * v + 0.5 * stiffness * (y - reed.rest_opening_m).powi(2) + potential(y)
    };
    let mut transitions = [false; 4];
    for chi in [0.0, 5.0] {
        for dt in [1e-5, 1e-4] {
            let lay = crate::unilateral_contact::slit_lay(contact_stiffness, 2.0)
                .unwrap()
                .with_internal_loss(chi)
                .unwrap();
            let owner = PreparedReedLoad::new(reed, 1.2, dt, lay).unwrap();
            for (y, v) in [
                (-1e-4, -0.2),
                (-1e-4, 2.0),
                (1e-5, -2.0),
                (1e-5, 0.1),
                (4e-4, 0.0),
            ] {
                for mouth in [-500.0, 300.0, 10000.0] {
                    let free = 75.0;
                    let impedance = 1e6;
                    let step = owner.step(y, v, mouth, free, impedance).unwrap();
                    let y1 = step.opening_m;
                    let v1 = step.velocity_m_s;
                    let vm = (y1 - y) / dt;
                    let drop = mouth - step.pressure_pa;
                    let jet = reed.width_m
                        * (0.5 * (y + y1)).max(0.0)
                        * drop.signum()
                        * (2.0 * drop.abs() / 1.2).sqrt();
                    let flow = jet - face * vm;
                    let flow_scale = jet.abs() + (face * vm).abs() + flow.abs();
                    assert_close(step.mean_velocity_m_s, vm, vm.abs().max(v.abs()));
                    assert_close(step.jet_flow_m3_s, jet, flow_scale);
                    assert_close(step.volume_flow_m3_s, flow, flow_scale);
                    assert_close(
                        step.pressure_pa,
                        free + impedance * flow,
                        free.abs() + (impedance * flow).abs(),
                    );

                    // Independent secant of the original power-law potential,
                    // and its original nonadhesive Hunt--Crossley force.
                    let elastic = if y1 == y {
                        contact_stiffness * (-y).max(0.0).powi(2)
                    } else {
                        (potential(y) - potential(y1)) / (y1 - y)
                    };
                    let contact = (elastic * (1.0 - chi * vm)).max(0.0);
                    let force = contact - face * drop;
                    let sigma = damping / (2.0 * mass);
                    let omega2 = stiffness / mass;
                    let wd = (omega2 - sigma * sigma).sqrt();
                    let displacement = y - reed.rest_opening_m;
                    let equilibrium = force / stiffness;
                    let excursion = displacement - equilibrium;
                    let velocity = |t: f64| {
                        (-sigma * t).exp()
                            * (v * (wd * t).cos()
                                + (-sigma * v - omega2 * excursion) / wd * (wd * t).sin())
                    };
                    let expected_y = reed.rest_opening_m
                        + equilibrium
                        + (-sigma * dt).exp()
                            * (excursion * (wd * dt).cos()
                                + (v + sigma * excursion) / wd * (wd * dt).sin());
                    assert_close(y1, expected_y, y.abs() + y1.abs() + reed.rest_opening_m);
                    assert_close(
                        v1,
                        velocity(dt),
                        v.abs() + v1.abs() + (dt * force / mass).abs(),
                    );
                    // Simpson quadrature of the independent analytical velocity
                    // is a positive heat reference; no energy remainder is used.
                    let intervals = 256;
                    let h = dt / f64::from(intervals);
                    let mut integral = velocity(0.0).powi(2) + velocity(dt).powi(2);
                    for index in 1..intervals {
                        integral += if index % 2 == 0 { 2.0 } else { 4.0 }
                            * velocity(f64::from(index) * h).powi(2);
                    }
                    let viscous_heat = damping * h * integral / 3.0;
                    let contact_heat = (elastic - contact) * (y1 - y);
                    let jet_heat = dt * drop * jet;
                    let mouth_work = dt * mouth * flow;
                    let bore_work = dt * step.pressure_pa * flow;
                    let delta_h = energy(y1, v1) - energy(y, v);
                    let scale = energy(y, v).abs()
                        + energy(y1, v1).abs()
                        + mouth_work.abs()
                        + bore_work.abs()
                        + jet_heat.abs();
                    let account = step.account;
                    assert_close(account.delta_h_j, delta_h, scale);
                    assert_close(account.viscous_heat_j, viscous_heat, scale);
                    assert_close(account.contact_heat_j, contact_heat, scale);
                    assert_close(account.jet_heat_j, jet_heat, scale);
                    assert_close(account.mouth_work_j, mouth_work, scale);
                    assert_close(account.bore_work_j, bore_work, scale);
                    let residual =
                        delta_h + viscous_heat + contact_heat + jet_heat + bore_work - mouth_work;
                    assert_close(residual, 0.0, scale);
                    assert!(viscous_heat >= 0.0 && contact_heat >= 0.0 && jet_heat >= 0.0);
                    transitions[usize::from(y >= 0.0) * 2 + usize::from(y1 >= 0.0)] = true;
                }
            }
        }
    }
    assert!(
        transitions.iter().all(|seen| *seen),
        "entry/release coverage: {transitions:?}"
    );
}

#[test]
fn rejected_loaded_trials_leave_the_prepared_response_unchanged() {
    let reed = reed();
    let lay = crate::unilateral_contact::slit_lay(1e8, 2.0).unwrap();
    let owner = PreparedReedLoad::new(reed, 1.2, 1.0 / 48000.0, lay).unwrap();
    let accepted = owner
        .step(reed.rest_opening_m, 0.0, 300.0, 75.0, 1e6)
        .unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(owner.step(bad, 0.0, 300.0, 75.0, 1e6).is_err());
        assert!(
            owner
                .step(reed.rest_opening_m, bad, 300.0, 75.0, 1e6)
                .is_err()
        );
        assert!(
            owner
                .step(reed.rest_opening_m, 0.0, bad, 75.0, 1e6)
                .is_err()
        );
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            owner
                .step(reed.rest_opening_m, 0.0, 300.0, 75.0, bad)
                .is_err()
        );
    }
    let replay = owner
        .step(reed.rest_opening_m, 0.0, 300.0, 75.0, 1e6)
        .unwrap();
    assert_eq!(accepted.opening_m, replay.opening_m);
    assert_eq!(accepted.velocity_m_s, replay.velocity_m_s);
    assert_eq!(accepted.pressure_pa, replay.pressure_pa);
    assert_eq!(accepted.volume_flow_m3_s, replay.volume_flow_m3_s);
    assert_eq!(accepted.account.residual_j, replay.account.residual_j);
}

#[test]
fn quasistatic_load_pays_for_reverse_flow_and_respects_complete_closure() {
    let reed = BeatingReed {
        mass_kg: 0.0,
        ..reed()
    };
    let dt = 1.0 / 48_000.0;
    let impedance = 1e6;
    for (mouth, free) in [(0.0, 75.0), (300.0, -75.0), (12_000.0, 0.0)] {
        let step = solve_quasistatic_load(reed, 1.2, dt, mouth, free, impedance, 0.0).unwrap();
        let drop = mouth - step.pressure_pa;
        let opening = reed.rest_opening_m * (1.0 - drop / reed.closing_pressure_pa).clamp(0.0, 1.0);
        let expected_flow =
            reed.width_m * opening * drop.signum() * (2.0 * drop.abs() / 1.2).sqrt();
        assert_close(step.volume_flow_m3_s, expected_flow, expected_flow.abs());
        assert_close(
            step.pressure_pa,
            free + impedance * expected_flow,
            free.abs() + (impedance * expected_flow).abs(),
        );
        let mouth_work = dt * mouth * expected_flow;
        let bore_work = dt * step.pressure_pa * expected_flow;
        let jet_heat = dt * drop * expected_flow;
        let scale = mouth_work.abs() + bore_work.abs() + jet_heat.abs();
        assert_close(step.mouth_work_j, mouth_work, scale);
        assert_close(step.bore_work_j, bore_work, scale);
        assert_close(step.jet_heat_j, jet_heat, scale);
        assert_close(jet_heat + bore_work - mouth_work, 0.0, scale);
        assert!(step.jet_heat_j >= 0.0);
        if mouth == 0.0 {
            assert!(step.volume_flow_m3_s < 0.0 && step.bore_work_j < 0.0);
            assert_eq!(step.mouth_work_j, 0.0);
        }
        if mouth > reed.closing_pressure_pa {
            assert_eq!(step.volume_flow_m3_s, 0.0);
            assert_eq!(step.pressure_pa, 0.0);
        }
    }
}
