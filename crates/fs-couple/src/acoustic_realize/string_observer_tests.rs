//! G1: original spatial source moments and accepted joined-clock histories.
use super::*;
use crate::thin_plate::CompactBody;

fn string() -> PrestressedString {
    PrestressedString {
        length_m: 0.5,
        tension_n: 40_000.0,
        lin_density_kg_m: 1.0,
        axial_stiffness_n: 0.0,
        width_m: 0.001,
        n_modes: 1,
        damping_ratio: 0.0,
        // An explicit viscous law keeps the independent oscillator's damping
        // known, without reproducing the cylinder-air-loss implementation.
        rayleigh: Some(RayleighParams { alpha_per_s: 12.0, beta_s: 0.0 }),
        bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        polarization_detune: 0.0,
        moving_end: true,
        tension_loss_factor: None,
        bending_loss_factor: None,
    }
}

#[test]
fn signed_strip_moments_match_the_original_sine_and_cosine_shapes() {
    let mut string = string();
    string.n_modes = 6;
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let range = 2.3;
    let scale = 1.25;
    let pi = core::f64::consts::PI;
    let phi0 = (2.0 / (string.lin_density_kg_m * string.length_m)).sqrt();
    let factor = scale * gas.density * string.width_m / (4.0 * pi * range);
    for moving in [false, true] {
        string.moving_end = moving;
        let observer = StringObserver::new(&string, &gas, range, scale);
        for (mode, &(area, moment)) in observer.weights.iter().enumerate() {
            // Independent Simpson integration of the actual normalized shape,
            // including its sign, not the observer's integrated formulas.
            let cells = 4096;
            let dx = string.length_m / cells as f64;
            let mut integrals = [0.0; 2];
            for i in 0..=cells {
                let x = i as f64 * dx;
                let phase = (mode as f64 + if moving { 0.5 } else { 1.0 }) * pi * x
                    / string.length_m;
                let shape = phi0 * if moving { phase.cos() } else { phase.sin() };
                let weight = if i == 0 || i == cells { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
                integrals[0] += weight * shape;
                integrals[1] += weight * (x - string.length_m / 2.0) * shape;
            }
            let expected_area = factor * integrals[0] * dx / 3.0;
            let expected_moment = factor * integrals[1] * dx / (3.0 * gas.sound_speed);
            assert!((area - expected_area).abs() < 1e-11 * factor * phi0 * string.length_m);
            assert!((moment - expected_moment).abs()
                < 1e-11 * factor * phi0 * string.length_m.powi(2) / gas.sound_speed);
        }
        if moving {
            let (area, moment) = observer.weights[0];
            let expected_area = factor * phi0 * 2.0 * string.length_m / pi;
            let expected_moment = factor * phi0 * string.length_m.powi(2)
                * (1.0 / pi - 4.0 / pi.powi(2)) / gas.sound_speed;
            assert!((area / expected_area - 1.0).abs() < 1e-14);
            assert!((moment / expected_moment - 1.0).abs() < 1e-14);
            assert!(moment < 0.0);
        }
    }
}

#[test]
fn moving_renderer_observes_the_independent_constrained_oscillator_history() {
    let string = string();
    let gas = gas_state(AmbientGas::sea_level()).unwrap();
    let pluck = Pluck { station_frac: 0.5, height_m: 1e-3 };
    let rate = 8_000;
    let dt = 1.0 / f64::from(rate);
    let pi = core::f64::consts::PI;
    let phi0 = (2.0 / (string.lin_density_kg_m * string.length_m)).sqrt();
    let omega_string_sq = string.tension_n / string.lin_density_kg_m
        * (pi / (2.0 * string.length_m)).powi(2);
    // Use the original release projection; the oracle below independently
    // solves the subsequent damped, constrained mechanical motion.
    let q0 = free_fixed_pluck_modal(pluck, &string, 0);
    let spec = RadiatingPlate {
        area_m2: 2e-4, mass_kg: 0.2, frequency_hz: 160.0, damping_ratio: 0.015,
    };
    for body in [false, true] {
        let mut bank = PlateBank::default();
        if body { bank.linear.push(CompactBody::from_radiator(spec).unwrap()); }
        let actual = realize_dirac_join(
            &string, Some(pluck), None, &[], None, &mut bank,
            None, None, None, &gas, 1.0, rate, 64,
        ).unwrap();
        // Common endpoint velocity gives q_body = a(q_string-q0), a=phi0 sqrt(M).
        // Eliminate that constraint in Newton's equations, including the
        // body's original equilibrium offset and both physical dampers.
        let a_sq = if body { phi0.powi(2) * spec.mass_kg } else { 0.0 };
        let omega_body = 2.0 * pi * spec.frequency_hz;
        let stiffness = omega_string_sq + a_sq * omega_body.powi(2);
        let omega_sq = stiffness / (1.0 + a_sq);
        let gamma = (12.0 + a_sq * 2.0 * spec.damping_ratio * omega_body)
            / (2.0 * (1.0 + a_sq));
        let omega_d = (omega_sq - gamma.powi(2)).sqrt();
        let amplitude = q0 * omega_string_sq / stiffness;
        let velocity = |time: f64| -amplitude * omega_sq / omega_d
            * (-gamma * time).exp() * (omega_d * time).sin();
        let area = phi0 * 2.0 * string.length_m / pi;
        let moment = phi0 * string.length_m.powi(2) * (1.0 / pi - 4.0 / pi.powi(2));
        let factor = gas.density * string.width_m / (4.0 * pi);
        let body_factor = if body { gas.density * spec.area_m2 * phi0 / (2.0 * pi) } else { 0.0 };
        let mut previous = None;
        let mut peak = 0.0_f64;
        let mut error = 0.0_f64;
        let mut dipole = 0.0_f64;
        for (i, pressure) in actual.iter().enumerate() {
            let t = i as f64 * dt;
            let acceleration = (velocity(t + dt) - velocity(t)) / dt;
            // The first accepted interval primes itself. Later samples use
            // consecutive interval means, without a release-history impulse.
            let jerk = previous.map_or(0.0, |a| (acceleration - a) / dt);
            let monopole = (factor * area + body_factor) * acceleration;
            let signed_moment = factor * moment / gas.sound_speed * jerk;
            let expected = monopole + signed_moment;
            peak = peak.max(expected.abs());
            error = error.max((pressure - expected).abs());
            dipole = dipole.max(signed_moment.abs());
            previous = Some(acceleration);
        }
        assert!(peak > 0.0 && dipole > 1e-3 * peak, "body={body}: missing-moment falsifier");
        assert!(error < 2e-9 * peak, "body={body}: error {error} vs peak {peak}");
    }
}
