//! G1/G3: Snell paths, stationary phase, reciprocity and physical observation.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::{
    PlanarFluidInterface, PreparedRefractedMonopole, RefractedMonopoleLimits,
};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_aeroac::observation::{HarmonicObservationBudget, PreparedHarmonicField};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn fluid(rho: f64, speed: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, speed).unwrap()
}
fn limits() -> RefractedMonopoleLimits {
    RefractedMonopoleLimits {
        minimum_kh: 10.0,
        minimum_direction_cosine: 0.05,
        path_tolerance_m: 1e-12,
        maximum_iterations: 64,
    }
}
fn plane(media: [FluidAcousticState; 2]) -> PlanarFluidInterface {
    PlanarFluidInterface::new(media[0], media[1], [0.0; 3], [0.0, 0.0, 3.0]).unwrap()
}
fn prepare(
    media: [FluidAcousticState; 2],
    source: [f64; 3],
    receiver: [f64; 3],
    hz: f64,
) -> PreparedRefractedMonopole {
    PreparedRefractedMonopole::new(plane(media), hz, source, receiver, limits(), &mut || false)
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-10),
        "actual={actual:e}, expected={expected:e}, relative tolerance={tolerance:e}"
    );
}
fn complex_close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-10),
        "actual={actual:?}, expected={expected:?}, relative tolerance={tolerance:e}"
    );
}

#[test]
fn equal_speed_separates_impedance_from_bending_and_recovers_spherical_pressure() {
    let source = [-0.2, 0.3, -2.0];
    let receiver = [1.3, 1.0, 3.0];
    let q = C64::new(0.002, -0.001);
    let gas = fluid(1.2, 320.0);
    let frequency = 2000.0;
    let reference = PreparedCompactSourceObserver::new(gas, frequency, receiver, &[source], 1)
        .unwrap()
        .evaluate(&[CompactSourcePhasor::monopole(q)])
        .unwrap();
    let range = (1.5_f64.powi(2) + 0.7_f64.powi(2) + 5.0_f64.powi(2)).sqrt();
    for rho2 in [1.2, 4.8] {
        let observer = prepare([gas, fluid(rho2, 320.0)], source, receiver, frequency);
        let path = observer.path();
        let tp = 2.0 * rho2 / (1.2 + rho2);
        close(path.spreading_distance_m, range, 2e-13);
        close(path.delay_s, range / 320.0, 2e-13);
        close(path.pressure_transmission, tp, 2e-13);
        close(path.interface_position_m[0], 0.4, 5e-12);
        close(path.interface_position_m[1], 0.58, 5e-12);
        let actual = observer.evaluate(q).unwrap();
        complex_close(actual.pressure_pa, reference.pressure_pa.scale(tp), 3e-12);
        for axis in 0..3 {
            close(
                path.incident_direction[axis],
                path.receiving_direction[axis],
                3e-12,
            );
        }
        // Pressure is exact in this limit; the declared leading-order velocity
        // omits the known spherical 1/(kr) reactive term, never claims it away.
        if rho2 == 1.2 {
            let kr = TAU * frequency * range / 320.0;
            let velocity_error =
                (actual.particle_velocity_m_s[2] - reference.particle_velocity_m_s[2]).abs();
            close(
                velocity_error,
                actual.particle_velocity_m_s[2].abs() / kr,
                2e-9,
            );
        }
        assert_eq!(observer.interface().receiving_state(), fluid(rho2, 320.0));
        assert_eq!(observer.source_position_m(), source);
        assert_eq!(observer.observer_position_m(), receiver);
    }
}

#[test]
fn constructed_oblique_rays_obey_snell_and_normal_power_continuity() {
    for (c1, c2, rho1, rho2, s1) in [
        (300.0_f64, 450.0, 1.2, 2.5, 0.4),
        (450.0, 300.0, 2.5, 1.2, 0.7),
        (300.0, 450.0, 2.0, 4.0 / 3.0, 0.4),
    ] {
        let s2 = s1 * c2 / c1;
        let u1 = (1.0 - s1 * s1).sqrt();
        let u2 = (1.0 - s2 * s2).sqrt();
        let x = 2.0 * s1 / u1;
        let r = x + 3.0 * s2 / u2;
        let observer = prepare(
            [fluid(rho1, c1), fluid(rho2, c2)],
            [0.0, 0.0, -2.0],
            [r, 0.0, 3.0],
            2000.0,
        );
        let path = observer.path();
        close(path.interface_position_m[0], x, 3e-12);
        close(path.incident_direction[0], s1, 3e-12);
        close(path.receiving_direction[0], s2, 3e-12);
        close(path.delay_s, 2.0 / (u1 * c1) + 3.0 / (u2 * c2), 2e-13);
        let tp = path.pressure_transmission;
        let power = tp * tp * (rho1 * c1 / (rho2 * c2)) * (u2 / u1);
        close(path.normal_power_transmission, power, 3e-13);
        close((tp - 1.0).powi(2) + power, 1.0, 3e-13);
        assert!(path.crossing_half_width_m <= limits().path_tolerance_m);
        if rho1 == 2.0 {
            assert!(
                (tp - 1.0).abs() > 0.05,
                "equal scalar Z does not remove oblique reflection"
            );
        }
    }
}

#[test]
fn normal_incidence_uses_refracted_spreading_instead_of_total_distance() {
    let media = [fluid(1.2, 300.0), fluid(2.1, 600.0)];
    let observer = prepare(media, [0.0, 0.0, -2.0], [0.0, 0.0, 3.0], 1700.0);
    let path = observer.path();
    let d = 2.0 + (600.0 / 300.0) * 3.0;
    close(path.spreading_distance_m, d, 2e-15);
    let phase = TAU * 1700.0 * (2.0 / 300.0 + 3.0 / 600.0);
    let (sin, cos) = phase.sin_cos();
    let tp = 2.0 * (2.1 * 600.0) / (1.2 * 300.0 + 2.1 * 600.0);
    let expected = C64::new(sin, -cos).scale(TAU * 1700.0 * 1.2 * tp / (4.0 * PI * d));
    complex_close(
        observer.evaluate(C64::ONE).unwrap().pressure_pa,
        expected,
        2e-13,
    );
    assert_eq!(path.interface_position_m, [0.0; 3]);
    assert_eq!(path.crossing_half_width_m, 0.0);
}

#[test]
fn pressure_transfer_is_reciprocal_when_source_fluids_and_normal_are_reversed() {
    let media = [fluid(1.2, 300.0), fluid(4.3, 510.0)];
    let source = [-0.3, 0.4, -2.0];
    let receiver = [1.7, -0.5, 3.0];
    let forward = prepare(media, source, receiver, 1900.0);
    let reverse_plane =
        PlanarFluidInterface::new(media[1], media[0], [0.0; 3], [0.0, 0.0, -1.0]).unwrap();
    let reverse = PreparedRefractedMonopole::new(
        reverse_plane,
        1900.0,
        receiver,
        source,
        limits(),
        &mut || false,
    )
    .unwrap();
    close(forward.path().delay_s, reverse.path().delay_s, 2e-13);
    complex_close(
        forward.evaluate(C64::ONE).unwrap().pressure_pa,
        reverse.evaluate(C64::ONE).unwrap().pressure_pa,
        3e-12,
    );
}

#[test]
fn stationary_phase_hessian_independently_recovers_the_refracted_amplitude() {
    let (rho1, rho2, c1, c2, a, b, s1) = (1.2_f64, 2.7, 300.0, 500.0, 1.5, 2.0, 0.35);
    let omega = TAU * 1500.0;
    let (k1, k2) = (omega / c1, omega / c2);
    let q = k1 * s1;
    let (q1, q2) = ((k1 * k1 - q * q).sqrt(), (k2 * k2 - q * q).sqrt());
    let r = a * q / q1 + b * q / q2;
    let observer = prepare(
        [fluid(rho1, c1), fluid(rho2, c2)],
        [0.0, 0.0, -a],
        [r, 0.0, b],
        1500.0,
    );
    // Independently differentiate the Weyl phase in transverse wavevector,
    // not either principal-width expression used by the implementation.
    let phase = |qx: f64, qy: f64| {
        r * qx + a * (k1 * k1 - qx * qx - qy * qy).sqrt() + b * (k2 * k2 - qx * qx - qy * qy).sqrt()
    };
    let phi = phase(q, 0.0);
    let h = 1e-4 * k1;
    let hxx = (phase(q + h, 0.0) - 2.0 * phi + phase(q - h, 0.0)) / (h * h);
    let hyy = (phase(q, h) - 2.0 * phi + phase(q, -h)) / (h * h);
    let transmission = 2.0 * (q1 / rho1) / (q1 / rho1 + q2 / rho2);
    // 2D stationary phase: Gt = T exp(i phi)/(4 pi q1 sqrt(det(-H))).
    let green_amplitude = transmission / (4.0 * PI * q1 * (hxx * hyy).sqrt());
    let (sin, cos) = phi.sin_cos();
    let expected = C64::new(sin, -cos).scale(omega * rho1 * green_amplitude);
    complex_close(
        observer.evaluate(C64::ONE).unwrap().pressure_pa,
        expected,
        3e-6,
    );
}

#[test]
fn coordinate_rotation_and_distance_scaling_preserve_physical_fields() {
    let media = [fluid(1.2, 300.0), fluid(2.3, 450.0)];
    let source = [0.0, 0.0, -2.0];
    let receiver = [1.0, 0.5, 3.0];
    let original = prepare(media, source, receiver, 2000.0);
    let transform = |p: [f64; 3]| [10.0 + p[2], -2.0 + p[0], 3.0 + p[1]];
    let rotated_plane =
        PlanarFluidInterface::new(media[0], media[1], transform([0.0; 3]), [7.0, 0.0, 0.0])
            .unwrap();
    let rotated = PreparedRefractedMonopole::new(
        rotated_plane,
        2000.0,
        transform(source),
        transform(receiver),
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = original.evaluate(C64::ONE).unwrap();
    let rotated_field = rotated.evaluate(C64::ONE).unwrap();
    complex_close(rotated_field.pressure_pa, field.pressure_pa, 3e-12);
    for (new, old) in [(0, 2), (1, 0), (2, 1)] {
        complex_close(
            rotated_field.particle_velocity_m_s[new],
            field.particle_velocity_m_s[old],
            3e-12,
        );
    }
    let scaled = prepare(
        media,
        source.map(|v| 10.0 * v),
        receiver.map(|v| 10.0 * v),
        2000.0,
    );
    close(scaled.path().delay_s, 10.0 * original.path().delay_s, 2e-13);
    close(
        scaled.evaluate(C64::ONE).unwrap().pressure_pa.abs(),
        field.pressure_pa.abs() / 10.0,
        3e-13,
    );
}

#[test]
fn original_finite_window_observer_consumes_the_refracted_field_on_its_source_clock() {
    let observer = prepare(
        [fluid(1.2, 300.0), fluid(2.1, 440.0)],
        [0.0, 0.0, -2.0],
        [1.0, 0.5, 3.0],
        1800.0,
    );
    let q = C64::new(0.002, -0.001);
    let line = observer.evaluate_harmonic_line(q).unwrap();
    let field = observer.evaluate(q).unwrap();
    assert_eq!(line.pressure_pa, field.pressure_pa);
    let observed = PreparedHarmonicField::new(
        &[line],
        HarmonicObservationBudget {
            max_lines: 1,
            max_line_pairs: 1,
        },
    )
    .unwrap()
    .observe_window(HarmonicWindow::new(0.137, 1.0 / observer.frequency_hz()).unwrap())
    .unwrap();
    close(
        observed.mean_square_pressure_pa2,
        0.5 * field.pressure_pa.norm_sq(),
        3e-12,
    );
    for axis in 0..3 {
        close(
            observed.mean_intensity_w_m2[axis],
            field.active_intensity_w_m2[axis],
            3e-12,
        );
    }
}

#[test]
fn invalid_envelopes_wrong_sides_exhaustion_and_cancellation_refuse() {
    let media = [fluid(1.2, 300.0), fluid(2.1, 450.0)];
    let interface = plane(media);
    assert!(PlanarFluidInterface::new(media[0], media[1], [0.0; 3], [0.0; 3]).is_err());
    let attempt = |hz, source, receiver, limits| {
        PreparedRefractedMonopole::new(interface, hz, source, receiver, limits, &mut || false)
    };
    let (source, receiver) = ([0.0, 0.0, -2.0], [3.0, 0.0, 3.0]);
    for bad in [
        RefractedMonopoleLimits {
            minimum_kh: 0.0,
            ..limits()
        },
        RefractedMonopoleLimits {
            minimum_direction_cosine: f64::NAN,
            ..limits()
        },
        RefractedMonopoleLimits {
            path_tolerance_m: -1.0,
            ..limits()
        },
        RefractedMonopoleLimits {
            maximum_iterations: 129,
            ..limits()
        },
    ] {
        assert!(attempt(2000.0, source, receiver, bad).is_err());
    }
    assert!(attempt(1.0, source, receiver, limits()).is_err());
    assert!(attempt(2000.0, receiver, source, limits()).is_err());
    assert!(attempt(2000.0, source, [0.0; 3], limits()).is_err());
    assert!(
        attempt(
            2000.0,
            source,
            receiver,
            RefractedMonopoleLimits {
                minimum_direction_cosine: 0.999,
                ..limits()
            }
        )
        .is_err()
    );
    assert!(matches!(
        attempt(
            2000.0,
            source,
            receiver,
            RefractedMonopoleLimits {
                maximum_iterations: 1,
                ..limits()
            }
        ),
        Err(AeroacError::NotConverged { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedRefractedMonopole::new(interface, 2000.0, source, receiver, limits(), &mut || {
            polls += 1;
            polls == 3
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let observer = attempt(2000.0, source, receiver, limits()).unwrap();
    assert!(observer.evaluate(C64::new(f64::INFINITY, 0.0)).is_err());
    assert_eq!(observer.evaluate(C64::ZERO).unwrap().pressure_pa, C64::ZERO);
}
