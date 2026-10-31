use super::super::{BaffledPressure, CircularOutletReceiver};
use super::*;

fn medium() -> RayleighMedium {
    RayleighMedium {
        density: 1.2,
        sound_speed: 320.0,
    }
}
fn patch(x: f64, shape: Vec<f64>) -> SurfaceSample {
    SurfaceSample {
        position_m: [x, 0.0, 0.0],
        area_m2: 0.08,
        mode_shape: shape,
    }
}
fn phase(value: f64) -> C64 {
    C64::new(value.cos(), value.sin())
}
fn close(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn baffled_monopole_retains_reactive_velocity_and_cardioid_proximity() {
    let gas = medium();
    let radius = 0.4;
    let surface = [patch(0.0, vec![1.0 / 0.08])]; // unit total volume flow
    let source = C64::new(0.03, -0.02);
    let front = [0.0, 0.0, -1.0];
    let patterns = [
        FirstOrder::default(),
        FirstOrder::new(0.5, front).unwrap(),
        FirstOrder::new(0.0, front).unwrap(),
        FirstOrder::new(0.5, [0.0, 0.0, 1.0]).unwrap(),
    ];
    for i in 0..=64 {
        let kr = 0.05 * 1000_f64.powf(i as f64 / 64.0);
        let omega = kr * gas.sound_speed / radius;
        let prepared = Prepared::new(
            &surface,
            1,
            &[[0.0, 0.0, radius]; 4],
            omega / TAU,
            gas,
            Options::default(),
            || false,
        )
        .unwrap();
        let mut fields = [Field::default(); 4];
        prepared.evaluate_into(&[source], &mut fields).unwrap();
        let expected_pressure =
            source * phase(kr) * C64::new(0.0, -gas.density * omega / (TAU * radius));
        let expected_ratio = C64::new(1.0, 1.0 / kr);
        close(fields[0].pressure_pa, expected_pressure, 3e-14);
        close(
            fields[0].particle_velocity_m_s[2]
                / expected_pressure.scale(1.0 / (gas.density * gas.sound_speed)),
            expected_ratio,
            3e-14,
        );
        close(
            fields[0].pressure_gradient_pa_m[2],
            expected_pressure * C64::new(-1.0 / radius, kr / radius),
            3e-14,
        );
        assert_eq!(fields[0].particle_velocity_m_s[..2], [C64::ZERO; 2]);
        let mut observed = [C64::ZERO; 4];
        prepared
            .observe_into(&[source], &patterns, &mut observed)
            .unwrap();
        assert_eq!(observed[0], fields[0].pressure_pa);
        close(
            observed[1],
            expected_pressure * C64::new(1.0, 0.5 / kr),
            3e-14,
        );
        close(observed[2], expected_pressure * expected_ratio, 3e-14);
        // An ideal rear cardioid null is a plane-wave limit; a nearby source
        // retains this quadrature-phase contribution instead of disappearing.
        close(
            observed[3],
            expected_pressure * C64::new(0.0, -0.5 / kr),
            2e-12,
        );
    }
    let kr = 100_000.0;
    let prepared = Prepared::new(
        &surface,
        1,
        &[[0.0, 0.0, radius]],
        kr * gas.sound_speed / (TAU * radius),
        gas,
        Options::default(),
        || false,
    )
    .unwrap();
    let field = prepared.transfer(0, 0).unwrap();
    let ratio = field.particle_velocity_m_s[2]
        / field
            .pressure_pa
            .scale(1.0 / (gas.density * gas.sound_speed));
    assert!((ratio - C64::ONE).abs() < 1.00001e-5);
}

#[test]
fn spatial_gradient_differentiates_complex_pressure_in_every_direction() {
    let surface = [patch(-0.18, vec![1.0, -0.4]), patch(0.23, vec![-0.2, 0.8])];
    let center = [0.12, 0.17, 0.35];
    let step = 1e-6;
    let mut positions = vec![center];
    for coordinate in 0..3 {
        let mut plus = center;
        plus[coordinate] += step;
        let mut minus = center;
        minus[coordinate] -= step;
        positions.extend([plus, minus]);
    }
    for frequency in [70.0, 700.0, 7000.0] {
        let prepared = Prepared::new(
            &surface,
            2,
            &positions,
            frequency,
            medium(),
            Options::default(),
            || false,
        )
        .unwrap();
        let mut fields = vec![Field::default(); positions.len()];
        prepared
            .evaluate_into(&[C64::new(0.7, 0.2), C64::new(-0.3, 0.9)], &mut fields)
            .unwrap();
        for c in 0..3 {
            let independent =
                (fields[1 + 2 * c].pressure_pa - fields[2 + 2 * c].pressure_pa).scale(0.5 / step);
            close(fields[0].pressure_gradient_pa_m[c], independent, 2e-8);
        }
    }
}

#[test]
fn coherent_pressure_null_keeps_the_physical_transverse_velocity() {
    let surface = [patch(-0.2, vec![1.0, 0.0]), patch(0.2, vec![0.0, 1.0])];
    let prepared = Prepared::new(
        &surface,
        2,
        &[[0.0, 0.0, 0.3]; 2],
        2300.0,
        medium(),
        Options::default(),
        || false,
    )
    .unwrap();
    let source = C64::new(0.2, -0.7);
    let mut fields = [Field::default(); 2];
    prepared
        .evaluate_into(&[source, -source], &mut fields)
        .unwrap();
    assert_eq!(fields[0].pressure_pa, C64::ZERO);
    assert_eq!(fields[0].particle_velocity_m_s[1..], [C64::ZERO; 2]);
    assert!(fields[0].particle_velocity_m_s[0].abs() > 0.1);
    let dipole = FirstOrder::new(0.0, [1.0, 0.0, 0.0]).unwrap();
    let mut output = [C64::ZERO; 2];
    prepared
        .observe_into(
            &[source, -source],
            &[FirstOrder::default(), dipole],
            &mut output,
        )
        .unwrap();
    assert_eq!(output[0], C64::ZERO);
    close(
        output[1],
        fields[0].particle_velocity_m_s[0].scale(-medium().density * medium().sound_speed),
        1e-15,
    );
    let rotated = phase(1.234);
    let mut next = [Field::default(); 2];
    prepared
        .evaluate_into(&[source * rotated, -source * rotated], &mut next)
        .unwrap();
    close(
        next[0].particle_velocity_m_s[0],
        fields[0].particle_velocity_m_s[0] * rotated,
        1e-15,
    );
}

#[test]
fn existing_circular_quadrature_converges_to_analytic_piston_pressure_and_velocity() {
    let radius = 0.08_f64;
    let height = 0.2_f64;
    let frequency = 800.0;
    let gas = medium();
    let k = TAU * frequency / gas.sound_speed;
    let edge = (height * height + radius * radius).sqrt();
    let area = std::f64::consts::PI * radius * radius;
    let expected_p =
        (phase(k * height) - phase(k * edge)).scale(gas.density * gas.sound_speed / area);
    let expected_u = (phase(k * height) - phase(k * edge).scale(height / edge)).scale(1.0 / area);
    let mut previous_error = f64::INFINITY;
    for rings in [8, 32, 128] {
        let geometry = CircularOutletReceiver {
            position_m: [0.0, 0.0, height],
            radial_rings: rings,
            angular_points: 64,
            maximum_frequency_hz: frequency,
        };
        let surface = BaffledPressure::circular_surface(radius, geometry, gas).unwrap();
        let prepared = Prepared::new(
            &surface,
            1,
            &[geometry.position_m],
            frequency,
            gas,
            Options::default(),
            || false,
        )
        .unwrap();
        let field = prepared.transfer(0, 0).unwrap();
        let error = ((field.pressure_pa - expected_p).abs() / expected_p.abs())
            .max((field.particle_velocity_m_s[2] - expected_u).abs() / expected_u.abs());
        assert!(
            error < previous_error * 0.1,
            "rings={rings}, error={error}, previous={previous_error}"
        );
        previous_error = error;
        assert!(field.particle_velocity_m_s[0].abs() < 1e-13);
        assert!(field.particle_velocity_m_s[1].abs() < 1e-13);
    }
    assert!(
        previous_error < 2e-6,
        "finest relative error {previous_error}"
    );
}

#[test]
fn harmonic_phasors_match_the_original_causal_pressure_owner_through_upper_partials() {
    let rate = 48_000;
    let gas = medium();
    let position = [0.0, 0.0, 75.5 * gas.sound_speed / f64::from(rate)];
    let surface = [patch(0.0, vec![1.0, -0.45])];
    let sources = [C64::new(0.003, 0.002), C64::new(-0.002, 0.001)];
    for frequency in [200.0, 4800.0, 15000.0, 21600.0] {
        let harmonic = Prepared::new(
            &surface,
            2,
            &[position],
            frequency,
            gas,
            Options::default(),
            || false,
        )
        .unwrap();
        let mut fields = [Field::default()];
        harmonic.evaluate_into(&sources, &mut fields).unwrap();
        let mut causal =
            BaffledPressure::new_differentiating(&surface, 2, rate, position, gas).unwrap();
        let latency = causal.filter_delay_samples();
        let omega = TAU * frequency / f64::from(rate);
        let (mut error, mut signal) = (0.0, 0.0);
        for n in 0..2400 {
            let angle = omega * n as f64;
            let input = sources.map(|v| v.re * angle.cos() + v.im * angle.sin());
            let heard = causal.step(&input).unwrap();
            if n > 500 {
                let angle = omega * (n as f64 - latency);
                let expected =
                    fields[0].pressure_pa.re * angle.cos() + fields[0].pressure_pa.im * angle.sin();
                error += (heard - expected).powi(2);
                signal += expected.powi(2);
            }
        }
        assert!(
            (error / signal).sqrt() < 5e-4,
            "frequency={frequency}, error={}",
            (error / signal).sqrt()
        );
    }
}

#[test]
fn complete_receiver_arrays_and_budgets_keep_coherent_order_and_transactional_outputs() {
    let mut shape = vec![0.0; MAX_RAYLEIGH_INPUTS];
    shape[MAX_RAYLEIGH_INPUTS - 1] = -0.3;
    let surface = [patch(0.0, shape)];
    let positions: Vec<_> = (0..MAX_RECEIVERS)
        .map(|i| [i as f64 * 0.001, 0.0, 0.4])
        .collect();
    let prepared = Prepared::new(
        &surface,
        MAX_RAYLEIGH_INPUTS,
        &positions,
        4000.0,
        medium(),
        Options::default(),
        || false,
    )
    .unwrap();
    assert_eq!(
        (prepared.input_count(), prepared.receiver_count()),
        (MAX_RAYLEIGH_INPUTS, MAX_RECEIVERS)
    );
    assert_eq!(prepared.frequency_hz(), 4000.0);
    assert_eq!(prepared.medium(), medium());
    assert!(prepared.transfer(MAX_RECEIVERS, 0).is_none());
    let mut sources = vec![C64::ZERO; MAX_RAYLEIGH_INPUTS];
    sources[MAX_RAYLEIGH_INPUTS - 1] = C64::ONE;
    let mut output = vec![Field::default(); MAX_RECEIVERS];
    prepared.evaluate_into(&sources, &mut output).unwrap();
    for (i, field) in output.iter().enumerate() {
        assert_eq!(
            *field,
            prepared.transfer(i, MAX_RAYLEIGH_INPUTS - 1).unwrap()
        );
    }
    let reversed: Vec<_> = positions.iter().rev().copied().collect();
    let reverse = Prepared::new(
        &surface,
        MAX_RAYLEIGH_INPUTS,
        &reversed,
        4000.0,
        medium(),
        Options::default(),
        || false,
    )
    .unwrap();
    let mut reverse_output = vec![Field::default(); MAX_RECEIVERS];
    reverse
        .evaluate_into(&sources, &mut reverse_output)
        .unwrap();
    assert!(output.iter().eq(reverse_output.iter().rev()));
    let before = output.clone();
    sources[MAX_RAYLEIGH_INPUTS - 1] = C64::new(f64::MAX, 0.0);
    assert!(prepared.evaluate_into(&sources, &mut output).is_err());
    assert_eq!(output, before);
    sources[MAX_RAYLEIGH_INPUTS - 1] = C64::new(f64::NAN, 0.0);
    assert!(prepared.evaluate_into(&sources, &mut output).is_err());
    assert_eq!(output, before);
    let budget = Options {
        maximum_kernel_terms: MAX_RAYLEIGH_INPUTS * MAX_RECEIVERS - 1,
        ..Options::default()
    };
    assert!(
        Prepared::new(
            &surface,
            MAX_RAYLEIGH_INPUTS,
            &positions,
            4000.0,
            medium(),
            budget,
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(
        Prepared::new(
            &surface,
            MAX_RAYLEIGH_INPUTS,
            &positions,
            4000.0,
            medium(),
            Options::default(),
            || {
                polls += 1;
                polls == 5
            }
        )
        .is_err()
    );
    assert_eq!(polls, 5);
    assert!(
        Prepared::new(
            &surface,
            MAX_RAYLEIGH_INPUTS,
            &[[0.0, 0.0, 0.049]],
            4000.0,
            medium(),
            Options::default(),
            || false
        )
        .is_err()
    );
    assert!(
        Prepared::new(
            &surface,
            MAX_RAYLEIGH_INPUTS,
            &positions,
            0.0,
            medium(),
            Options::default(),
            || false
        )
        .is_err()
    );
}
