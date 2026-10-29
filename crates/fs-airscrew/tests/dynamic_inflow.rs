//! G0/G1 physical wake inertia, independent blade forces and analytic transients.
use fs_airscrew::dynamic_inflow::UniformAxialInflow;
use fs_airscrew::{BladeStation, Rotor};

const RADIUS: f64 = 1.1;
const RHO: f64 = 1.225;

fn rotor() -> Rotor {
    Rotor {
        radius_m: RADIUS,
        n_blades: 2,
        camber_ratio: 0.02,
        stations: [(0.3, 0.1, 0.6), (0.55, 0.15, 0.42), (0.8, 0.12, 0.32)]
            .map(|(r_over_r, chord_m, beta_rad)| BladeStation {
                r_over_r,
                chord_m,
                beta_rad,
            })
            .to_vec(),
    }
}

fn model() -> UniformAxialInflow {
    UniformAxialInflow::new(RADIUS, RHO).unwrap()
}

fn coefficient() -> f64 {
    2.0 * RHO * std::f64::consts::PI * RADIUS.powi(2)
}

fn independent_attached_loads(axial: f64, induced: f64, omega: f64) -> [f64; 2] {
    let rotor = rotor();
    let mut total = [0.0; 2];
    // Independent station-cell boundaries for this concrete radial fixture.
    for (station, width) in rotor.stations.iter().zip([0.185, 0.25, 0.325]) {
        let r = RADIUS * station.r_over_r;
        let phi = (axial + induced).atan2(omega * r);
        let alpha = station.beta_rad - phi;
        assert!(alpha.abs() < 0.25);
        let cl = std::f64::consts::TAU * (alpha + 0.04);
        let cd = 0.012 + 0.8 * alpha * alpha;
        let area = station.chord_m * width * RADIUS * 2.0;
        let pressure = 0.5 * RHO * ((axial + induced).powi(2) + (omega * r).powi(2));
        total[0] += pressure * area * (cl * phi.cos() - cd * phi.sin());
        total[1] += pressure * area * r * (cl * phi.sin() + cd * phi.cos());
    }
    total
}

#[test]
fn stored_induction_drives_actual_blade_forces_and_physical_air_inertia() {
    let model = model();
    let air_mass = (8.0 / 3.0) * RHO * RADIUS.powi(3);
    assert!((model.apparent_mass_kg() - air_mass).abs() < 1e-15);
    for induction in [0.0, 1.2, 2.0] {
        let result = model.evaluate(&rotor(), 6.0, induction, 45.0).unwrap();
        let expected = independent_attached_loads(6.0, induction, 45.0);
        assert!((result.thrust_n - expected[0]).abs() < 2e-12);
        assert!((result.torque_nm - expected[1]).abs() < 2e-12);
        let flux = coefficient() * induction * (6.0 + induction);
        assert!((result.inflow.momentum_flux_n - flux).abs() < 2e-13);
        assert!((air_mass * result.inflow.acceleration_mps2 + flux - expected[0]).abs() < 2e-12);
    }
    let fresh = model.evaluate(&rotor(), 6.0, 0.0, 45.0).unwrap();
    let developed = model.evaluate(&rotor(), 6.0, 2.0, 45.0).unwrap();
    assert!(fresh.thrust_n > developed.thrust_n + 20.0);
    assert!(fresh.inflow.acceleration_mps2 > 0.0);
}

#[test]
fn zero_and_signed_steady_limits_preserve_momentum_without_a_flow_floor() {
    let model = model();
    let zero = model.step_held(0.0, 0.0, 0.0, 0.01).unwrap();
    assert_eq!(zero.induced_mps, 0.0);
    assert_eq!(zero.momentum_residual_ns, 0.0);
    assert_eq!(zero.energy_residual_j, 0.0);
    for (axial, induction) in [(0.0, 3.0), (8.0, 2.0), (8.0, -2.0)] {
        let thrust = coefficient() * induction * (axial + induction);
        let rate = model.rate(axial, induction, thrust).unwrap();
        assert!(rate.acceleration_mps2.abs() < 1e-14);
        let next = model.step_held(induction, axial, thrust, 0.03).unwrap();
        assert!((next.induced_mps - induction).abs() < 2e-15);
        assert!(next.momentum_residual_ns.abs() < 1e-13);
        assert!(next.energy_residual_j.abs() < 1e-12);
        assert!(next.transport_work_j >= 0.0);
    }
}

#[test]
fn hover_thrust_step_matches_analytic_tanh_with_second_order_and_closed_work() {
    let model = model();
    let equilibrium = 4.0;
    let thrust = coefficient() * equilibrium * equilibrium;
    let end_time = 0.2;
    let exact =
        equilibrium * (coefficient() / model.apparent_mass_kg() * equilibrium * end_time).tanh();
    let mut errors = Vec::new();
    for count in [10, 20, 40] {
        let dt = end_time / f64::from(count);
        let mut induction = 0.0;
        let mut impulse = 0.0;
        let mut work = 0.0;
        for _ in 0..count {
            let next = model.step_held(induction, 0.0, thrust, dt).unwrap();
            assert!(next.induced_mps > induction && next.induced_mps < equilibrium);
            assert!(next.momentum_residual_ns.abs() < 2e-14);
            assert!(next.energy_residual_j.abs() < 1e-13);
            impulse += next.thrust_impulse_ns - next.transport_impulse_ns;
            work += next.induced_work_j - next.transport_work_j;
            induction = next.induced_mps;
        }
        assert!((model.apparent_mass_kg() * induction - impulse).abs() < 2e-13);
        assert!((0.5 * model.apparent_mass_kg() * induction * induction - work).abs() < 2e-12);
        errors.push((induction - exact).abs());
    }
    println!("hover analytic inflow errors: {errors:?}");
    assert!(errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8);
    assert!(errors[0] / errors[1] < 4.2 && errors[1] / errors[2] < 4.2);
}

#[test]
fn unpowered_wake_decays_by_momentum_transport_without_instant_reset() {
    let model = model();
    let initial = 4.0;
    let end_time = 0.2;
    let exact = initial / (1.0 + coefficient() / model.apparent_mass_kg() * initial * end_time);
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut induction = initial;
        for _ in 0..count {
            let next = model
                .step_held(induction, 0.0, 0.0, end_time / f64::from(count))
                .unwrap();
            assert!(next.induced_mps > 0.0 && next.induced_mps < induction);
            assert_eq!(next.induced_work_j, 0.0);
            assert!(next.transport_work_j > 0.0);
            assert!(next.kinetic_after_j < next.kinetic_before_j);
            induction = next.induced_mps;
        }
        errors.push((induction - exact).abs());
    }
    println!("free-wake analytic inflow errors: {errors:?}");
    assert!(errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9);
}

#[test]
fn actual_blade_equilibrium_matches_independent_uniform_disk_balance() {
    let model = model();
    let mut low = 0.0;
    let mut high = 3.0;
    for _ in 0..60 {
        let middle = 0.5 * (low + high);
        let thrust = independent_attached_loads(6.0, middle, 45.0)[0];
        if thrust > coefficient() * middle * (6.0 + middle) {
            low = middle;
        } else {
            high = middle;
        }
    }
    let induction = 0.5 * (low + high);
    assert!(induction > 0.5 && induction < 2.5);
    let actual = model.evaluate(&rotor(), 6.0, induction, 45.0).unwrap();
    assert!(actual.inflow.acceleration_mps2.abs() < 5e-13);
    let after = model
        .step_held(induction, 6.0, actual.thrust_n, 0.02)
        .unwrap();
    assert!((after.induced_mps - induction).abs() < 1e-14);
}

#[test]
fn invalid_or_excessive_braking_steps_refuse_instead_of_clipping() {
    let model = model();
    assert_eq!(
        model.rate(-1.0, 2.0, 1.0).unwrap_err().code,
        "dynamic-inflow-invalid"
    );
    assert_eq!(
        model.rate(10.0, -4.1, -50.0).unwrap_err().code,
        "dynamic-inflow-wake-domain"
    );
    assert!(model.step_held(0.0, 0.0, -1.0, 0.001).is_err());
    assert!(model.evaluate(&rotor(), 6.0, 1.0, 0.0).is_err());
    assert!(model.step_held(1.0, 6.0, 50.0, f64::NAN).is_err());
}
