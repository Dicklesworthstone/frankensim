use super::super::aperture::{Aperture, Sample};
use super::*;
use fs_material::fluid_acoustics::{FluidAcousticState, thermoviscous::ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;

fn material(gamma: f64, scale: f64) -> ClassicalBulkAbsorption {
    let state = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(2.0, 10.0).unwrap(),
        FluidTransportState::try_new(2.0, 8.0, 0.6 * scale, 0.003 * scale).unwrap(),
        8.0 / gamma,
    )
    .unwrap();
    ClassicalBulkAbsorption::from_material(&state, 0.002 * scale).unwrap()
}

fn patch(x: f64, shape: Vec<f64>) -> SurfaceSample {
    SurfaceSample {
        position_m: [x, 0.0, 0.0],
        area_m2: 0.08,
        mode_shape: shape,
    }
}

fn exponential(k: C64, distance: f64) -> C64 {
    let phase = k.re * distance;
    C64::new(phase.cos(), phase.sin()).scale((-k.im * distance).exp())
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-200),
        "actual={actual:?}, expected={expected:?}"
    );
}

// Cp=Cv removes the energy equation's coupling exactly. Solve the ORIGINAL
// continuity + viscous momentum equations: s²=c²-i*omega*eta/rho,
// k=omega/s, Z_p=rho*c²/s. The square root below uses real arithmetic and
// retains its weak-loss imaginary part without subtracting adjacent values.
fn independent_isocaloric_wave(omega: f64, scale: f64) -> (C64, C64) {
    let imaginary = omega * 0.003 * scale;
    let real = (0.5 * (100.0_f64.hypot(imaginary) + 100.0)).sqrt();
    let speed = C64::new(real, -imaginary / (2.0 * real));
    (C64::from_re(omega) / speed, C64::from_re(200.0) / speed)
}

#[test]
fn thermoviscous_point_pressure_and_velocity_resolve_weak_to_strong_viscosity() {
    let source = C64::new(0.03, -0.02);
    let law = material(1.0, 1.0);
    let surface = [patch(0.0, vec![1.0 / 0.08])];
    for x in [1e-8, 1e-3, 1.0, 1e4, 1e12] {
        let omega = x * 100.0 / 0.003;
        let (k, zp) = independent_isocaloric_wave(omega, 1.0);
        let radius = 0.7 / k.re;
        let options = Options {
            minimum_clearance_m: radius * 0.5,
            ..Options::default()
        };
        let prepared = Prepared::new_thermoviscous(
            &surface,
            1,
            &[[0.0, 0.0, radius]],
            omega / TAU,
            law,
            options,
            || false,
        )
        .unwrap();
        assert_eq!(prepared.thermoviscous(), Some(&law));
        assert!(prepared.relaxation().is_none());
        close(prepared.wavenumber_per_m(), k, 2e-14);
        let green = source * exponential(k, radius).scale(1.0 / (TAU * radius));
        let k_zp = k * zp;
        let pressure = green * C64::new(k_zp.im, -k_zp.re);
        let velocity = green * C64::new(1.0 / radius + k.im, -k.re);
        let mut field = [Field::default()];
        prepared.evaluate_into(&[source], &mut field).unwrap();
        close(field[0].pressure_pa, pressure, 3e-14);
        close(field[0].particle_velocity_m_s[2], velocity, 3e-14);
        close(
            field[0].pressure_gradient_pa_m[2],
            pressure * C64::new(-1.0 / radius - k.im, k.re),
            5e-14,
        );
        assert_eq!(field[0].particle_velocity_m_s[..2], [C64::ZERO; 2]);
        let pattern = FirstOrder::new(0.5, [0.0, 0.0, -1.0]).unwrap();
        let mut observed = [C64::ZERO];
        prepared
            .observe_into(&[source], &[pattern], &mut observed)
            .unwrap();
        close(
            observed[0],
            pressure.scale(0.5) + velocity.scale(10.0),
            3e-14,
        );
        if x >= 1.0 {
            // rho*omega is the TRACTION coefficient; using it for pressure
            // overstates thermodynamic pressure by sqrt(1+x²).
            let traction_image = green * C64::new(0.0, -2.0 * omega);
            assert!((traction_image - pressure).abs() > 0.5 * pressure.abs());
        }
    }
}

#[test]
fn thermoviscous_gradient_closes_viscous_momentum_and_spatial_refinement() {
    let law = material(2.0, 1.0);
    let surface = [patch(-0.18, vec![1.0, -0.4]), patch(0.23, vec![-0.2, 0.8])];
    let source = [C64::new(0.7, 0.2), C64::new(-0.3, 0.9)];
    let center = [0.12, 0.17, 0.35];
    for omega in [0.3, 3.0, 30.0, 300.0] {
        let mut previous_error = f64::INFINITY;
        for step in [2e-3, 1e-3, 5e-4] {
            let mut positions = vec![center];
            for c in 0..3 {
                let mut plus = center;
                let mut minus = center;
                plus[c] += step;
                minus[c] -= step;
                positions.extend([plus, minus]);
            }
            let prepared = Prepared::new_thermoviscous(
                &surface,
                2,
                &positions,
                omega / TAU,
                law,
                Options::default(),
                || false,
            )
            .unwrap();
            let mut fields = vec![Field::default(); positions.len()];
            prepared.evaluate_into(&source, &mut fields).unwrap();
            let k = prepared.wavenumber_per_m();
            // Direct longitudinal momentum: grad(p)=(i*rho*omega-eta*k²)u.
            // eta=4*mu/3+zeta=0.006. Thermal effects enter k, not this identity.
            let momentum = C64::new(0.0, 2.0 * omega) - (k * k).scale(0.006);
            let mut error = 0.0_f64;
            for c in 0..3 {
                let gradient = fields[0].pressure_gradient_pa_m[c];
                close(
                    gradient,
                    momentum * fields[0].particle_velocity_m_s[c],
                    2e-13,
                );
                let difference = (fields[1 + 2 * c].pressure_pa - fields[2 + 2 * c].pressure_pa)
                    .scale(0.5 / step);
                error = error.max((gradient - difference).abs() / gradient.abs());
            }
            assert!(
                error < previous_error * 0.26,
                "omega={omega}, h={step}, error={error}"
            );
            previous_error = error;
        }
        assert!(previous_error < 1e-4);
    }
}

#[test]
fn thermoviscous_disk_area_average_converges_to_independent_complex_integrals() {
    let law = material(1.0, 100.0);
    let omega = 300.0;
    let (k, zp) = independent_isocaloric_wave(omega, 100.0);
    let (height, radius) = (0.08_f64, 0.06_f64);
    let outer = height.hypot(radius);
    let axial = exponential(k, height);
    let edge = exponential(k, outer);
    let area = std::f64::consts::PI * radius * radius;
    let pressure = zp * (axial - edge).scale(1.0 / area);
    let velocity = (axial - edge.scale(height / outer)).scale(1.0 / area);
    let surface = [patch(0.0, vec![1.0 / 0.08])];
    let mut previous_error = f64::INFINITY;
    for rings in [4, 8, 16, 32] {
        let aperture =
            Aperture::circular([0.0, 0.0, height], [0.0, 0.0, 1.0], radius, rings, 8).unwrap();
        let prepared = Prepared::new_thermoviscous_aperture(
            &surface,
            1,
            &aperture,
            omega / TAU,
            law,
            Options::default(),
            || false,
        )
        .unwrap();
        assert_eq!(prepared.receiver_count(), 1);
        let field = prepared.transfer(0, 0).unwrap();
        let error = ((field.pressure_pa - pressure).abs() / pressure.abs())
            .max((field.particle_velocity_m_s[2] - velocity).abs() / velocity.abs());
        assert!(
            error < previous_error * 0.26,
            "rings={rings}, error={error}"
        );
        previous_error = error;
        let k_zp = k * zp;
        close(
            field.pressure_gradient_pa_m[2],
            field.particle_velocity_m_s[2] * C64::new(-k_zp.im, k_zp.re),
            3e-14,
        );
        assert!(
            field.particle_velocity_m_s[..2]
                .iter()
                .all(|u| u.abs() < 1e-14 * velocity.abs())
        );
    }
    assert!(previous_error < 8e-5);
}

#[test]
fn thermoviscous_fields_converge_to_original_lossless_receiver() {
    let surface = [patch(-0.1, vec![0.8]), patch(0.2, vec![-0.3])];
    let position = [[0.1, 0.2, 0.4]];
    let medium = RayleighMedium {
        density: 2.0,
        sound_speed: 10.0,
    };
    let lossless = Prepared::new(
        &surface,
        1,
        &position,
        80.0 / TAU,
        medium,
        Options::default(),
        || false,
    )
    .unwrap()
    .transfer(0, 0)
    .unwrap();
    let mut previous_error = f64::INFINITY;
    for scale in [1e-2, 1e-4, 1e-6] {
        let prepared = Prepared::new_thermoviscous(
            &surface,
            1,
            &position,
            80.0 / TAU,
            material(2.0, scale),
            Options::default(),
            || false,
        )
        .unwrap();
        assert_eq!(prepared.medium(), medium);
        let field = prepared.transfer(0, 0).unwrap();
        let mut error =
            (field.pressure_pa - lossless.pressure_pa).abs() / lossless.pressure_pa.abs();
        for c in 0..3 {
            error = error.max(
                (field.particle_velocity_m_s[c] - lossless.particle_velocity_m_s[c]).abs()
                    / lossless.particle_velocity_m_s[c].abs(),
            );
        }
        assert!(error < previous_error * 0.011);
        previous_error = error;
    }
    assert!(previous_error < 1e-7);
}

#[test]
fn thermoviscous_point_aperture_identity_and_transactional_refusals() {
    let law = material(2.0, 1.0);
    let surface = [patch(0.0, vec![1.0])];
    let position = [0.1, 0.2, 0.4];
    let point = Prepared::new_thermoviscous(
        &surface,
        1,
        &[position],
        30.0,
        law,
        Options::default(),
        || false,
    )
    .unwrap();
    let aperture = Aperture::new(vec![Sample {
        position_m: position,
        area_m2: 0.001,
    }])
    .unwrap();
    let area = Prepared::new_thermoviscous_aperture(
        &surface,
        1,
        &aperture,
        30.0,
        law,
        Options::default(),
        || false,
    )
    .unwrap();
    assert_eq!(point.transfer(0, 0), area.transfer(0, 0));
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            Prepared::new_thermoviscous(
                &surface,
                1,
                &[position],
                frequency,
                law,
                Options::default(),
                || false,
            )
            .is_err()
        );
    }
    assert!(
        Prepared::new_thermoviscous(
            &surface,
            1,
            &[position],
            30.0,
            law,
            Options::default(),
            || true,
        )
        .is_err()
    );
    let mut output = [point.transfer(0, 0).unwrap()];
    let before = output;
    assert!(
        point
            .evaluate_into(&[C64::new(f64::NAN, 0.0)], &mut output)
            .is_err()
    );
    assert_eq!(output, before);
    let mut observed = [C64::new(13.0, -7.0)];
    let before = observed;
    assert!(point.observe_into(&[C64::ONE], &[], &mut observed).is_err());
    assert_eq!(observed, before);
    let aperture = Aperture::new(vec![
        Sample {
            position_m: position,
            area_m2: 0.001,
        },
        Sample {
            position_m: [0.2, 0.1, 0.4],
            area_m2: 0.002,
        },
    ])
    .unwrap();
    assert!(
        Prepared::new_thermoviscous_aperture(
            &surface,
            1,
            &aperture,
            30.0,
            law,
            Options {
                maximum_kernel_terms: 1,
                ..Options::default()
            },
            || false,
        )
        .is_err()
    );
}
