//! G1/G3: physical force orientation, source derivatives and impedance-wall motion.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::CoherentCompactField;
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    PlanarImpedanceWall, PreparedImpedanceWallForce, PreparedImpedanceWallMonopole,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

const HZ: f64 = 700.0;
const RHO: f64 = 1.2;
const SPEED: f64 = 343.0;

fn limits() -> FullWaveForceLimits {
    FullWaveForceLimits {
        absolute_pressure_tolerance_pa_per_n: 1e-10,
        absolute_velocity_tolerance_m_s_per_n: 1e-13,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn wall(normalized: Option<C64>) -> PlanarImpedanceWall {
    PlanarImpedanceWall::new(
        FluidAcousticState::try_new(RHO, SPEED).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        normalized.map(|z| z.scale(RHO * SPEED)),
    )
    .unwrap()
}
fn unit(axis: usize) -> [C64; 3] {
    core::array::from_fn(|i| if i == axis { C64::ONE } else { C64::ZERO })
}
fn components(value: &CoherentCompactField) -> [C64; 4] {
    [
        value.pressure_pa,
        value.particle_velocity_m_s[0],
        value.particle_velocity_m_s[1],
        value.particle_velocity_m_s[2],
    ]
}
fn close(actual: C64, expected: C64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:?}, expected={expected:?}"
    );
}

// Independent spherical force Green function and its spatial Hessian.
fn force_green(source: [f64; 3], receiver: [f64; 3], force: [C64; 3]) -> [C64; 4] {
    let displacement: [f64; 3] = core::array::from_fn(|i| receiver[i] - source[i]);
    let radius = displacement.iter().map(|x| x * x).sum::<f64>().sqrt();
    let direction = displacement.map(|x| x / radius);
    let omega = TAU * HZ;
    let k = omega / SPEED;
    let (sin, cos) = (k * radius).sin_cos();
    let green = C64::new(cos, sin).scale(1.0 / (4.0 * PI * radius));
    let radial = green * C64::new(1.0 / radius, -k);
    let derivative = green * C64::new(k * k - 2.0 / (radius * radius), 2.0 * k / radius);
    let radial_force = (0..3).fold(C64::ZERO, |sum, i| sum + force[i].scale(direction[i]));
    let velocity: [C64; 3] = core::array::from_fn(|i| {
        (derivative * radial_force.scale(direction[i])
            + radial.scale(1.0 / radius) * (force[i] - radial_force.scale(direction[i])))
            * C64::new(0.0, -1.0 / (omega * RHO))
    });
    [radial * radial_force, velocity[0], velocity[1], velocity[2]]
}

#[test]
fn ideal_wall_force_images_retain_orientation_and_coherent_intensity() {
    let origin = [0.2, -0.1, 0.4];
    let normal = [0.6, 0.8, 0.0];
    let tangent = [-0.8, 0.6, 0.0];
    let source = core::array::from_fn(|i| origin[i] - 0.08 * normal[i] + 0.03 * tangent[i]);
    let receiver = [
        origin[0] - 0.13 * normal[0] + 0.2 * tangent[0],
        origin[1] - 0.13 * normal[1] + 0.2 * tangent[1],
        origin[2] + 0.04,
    ];
    let image = core::array::from_fn(|i| source[i] + 0.16 * normal[i]);
    let force = [C64::new(1.0, 0.3), C64::new(-0.2, 0.8), C64::new(0.7, -0.1)];
    let normal_force = (0..3).fold(C64::ZERO, |sum, i| sum + force[i].scale(normal[i]));
    let direct = force_green(source, receiver, force);
    for (impedance, sign) in [(None, 1.0), (Some(C64::ZERO), -1.0)] {
        let plane = PlanarImpedanceWall::new(
            FluidAcousticState::try_new(RHO, SPEED).unwrap(),
            origin,
            normal,
            impedance,
        )
        .unwrap();
        let owner = PreparedImpedanceWallForce::new(
            plane,
            HZ,
            source,
            receiver,
            FullWaveForceLimits {
                maximum_evaluations: 1,
                ..limits()
            },
            &mut || false,
        )
        .unwrap();
        let image_force =
            core::array::from_fn(|i| (force[i] - normal_force.scale(2.0 * normal[i])).scale(sign));
        let reflected = force_green(image, receiver, image_force);
        let expected_total = core::array::from_fn::<_, 4, _>(|i| direct[i] + reflected[i]);
        let total = owner.evaluate_total(force).unwrap();
        for (actual, expected) in [
            (
                components(&owner.evaluate_reflected(force).unwrap()),
                reflected,
            ),
            (components(&total), expected_total),
        ] {
            for i in 0..4 {
                close(actual[i], expected[i], 3e-12, 1e-12);
            }
        }
        for i in 0..3 {
            close(
                C64::new(
                    total.active_intensity_w_m2[i],
                    total.reactive_intensity_w_m2[i],
                ),
                expected_total[0].scale(0.5) * expected_total[i + 1].conj(),
                6e-12,
                1e-12,
            );
        }
        assert!(owner.diagnostics().exact_ideal_boundary);
        assert_eq!(owner.diagnostics().evaluations, 0);
    }
}

#[test]
fn finite_wall_force_is_source_derivative_of_all_monopole_field_components() {
    let plane = wall(Some(C64::new(0.3, 1.1)));
    let source = [0.01, -0.02, -0.045];
    let receiver = [0.12, 0.035, -0.035];
    let owner =
        PreparedImpedanceWallForce::new(plane, HZ, source, receiver, limits(), &mut || false)
            .unwrap();
    let mono_limits = FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-9,
        absolute_velocity_tolerance_per_m2: 1e-12,
        relative_tolerance: 1e-12,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    };
    let delta = 2e-4;
    for axis in 0..3 {
        let samples = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
            let mut shifted = source;
            shifted[axis] += offset * delta;
            let monopole = PreparedImpedanceWallMonopole::new(
                plane,
                HZ,
                shifted,
                receiver,
                mono_limits,
                &mut || false,
            )
            .unwrap();
            [
                components(&monopole.evaluate_reflected(C64::ONE).unwrap()),
                components(&monopole.evaluate_total(C64::ONE).unwrap()),
            ]
        });
        let actual = [
            components(&owner.evaluate_reflected(unit(axis)).unwrap()),
            components(&owner.evaluate_total(unit(axis)).unwrap()),
        ];
        for mode in 0..2 {
            for row in 0..4 {
                // Fourth-order original-source derivative, including velocity.
                let derivative = ((samples[2][mode][row] - samples[1][mode][row]).scale(8.0)
                    - (samples[3][mode][row] - samples[0][mode][row]))
                    .scale(1.0 / (12.0 * delta))
                    * C64::new(0.0, 1.0 / (TAU * HZ * RHO));
                close(actual[mode][row], derivative, 5e-7, 1e-9);
            }
        }
    }
}

#[test]
fn force_boundary_pressure_velocity_and_absorption_obey_the_rotated_robin_law() {
    let origin = [0.23, -0.17, 0.31];
    let normal = [0.6, 0.8, 0.0];
    let tangent = [-0.8, 0.6, 0.0];
    let source = [
        origin[0] - 0.07 * normal[0] + 0.12 * tangent[0],
        origin[1] - 0.07 * normal[1] + 0.12 * tangent[1],
        origin[2] + 0.03,
    ];
    for z in [
        C64::from_re(0.3),
        C64::new(3.0, -2.0),
        C64::new(0.02, 1.4),
        C64::new(0.0, -1.5),
    ] {
        let impedance = z.scale(RHO * SPEED);
        let plane = PlanarImpedanceWall::new(
            FluidAcousticState::try_new(RHO, SPEED).unwrap(),
            origin,
            normal,
            Some(impedance),
        )
        .unwrap();
        let owner =
            PreparedImpedanceWallForce::new(plane, HZ, source, origin, limits(), &mut || false)
                .unwrap();
        for force in [
            unit(0),
            unit(1),
            unit(2),
            [C64::new(0.3, 0.2), C64::new(-0.7, 0.1), C64::new(0.4, -0.6)],
        ] {
            let value = owner.evaluate_total(force).unwrap();
            let un = (0..3).fold(C64::ZERO, |sum, i| {
                sum + value.particle_velocity_m_s[i].scale(normal[i])
            });
            close(value.pressure_pa, impedance * un, 2e-8, 1e-8);
            let outward_flux = (0..3)
                .map(|i| value.active_intensity_w_m2[i] * normal[i])
                .sum::<f64>();
            let absorbed = 0.5 * impedance.re * un.abs().powi(2);
            assert!(outward_flux >= -1e-9);
            assert!((outward_flux - absorbed).abs() <= 2e-8 * absorbed.abs().max(1.0));
        }
    }
}

#[test]
fn all_force_channels_resolve_a_narrow_passive_surface_mode_with_bounded_work() {
    // z=.001+i has a passive evanescent amplitude peak near 2000. On axis,
    // independently differentiate the E1 expression F from the monopole test:
    // p_n=-F'/(4*pi), v_nn=-i*F''/(4*pi*rho*omega), and
    // v_tt=v_ss=-i*(k*k*F+F'')/(8*pi*rho*omega), by the Helmholtz equation.
    // Off-axis values use independent vertical-wavenumber adaptive quadrature.
    for (radius, expected) in [
        (
            0.0,
            [
                C64::ZERO,
                C64::new(35.5638177360633, 41.87343034999833),
                C64::new(0.08303818274805802, -0.14421949287333596),
                C64::new(0.08303818274805802, -0.14421949287333596),
                C64::ZERO,
                C64::ZERO,
                C64::new(0.10078699231219132, -0.298945881957846),
            ],
        ),
        (
            0.17,
            [
                C64::new(-21.122278195304883, 11.488377282009878),
                C64::new(-12.820503691755373, -8.769174824978263),
                C64::new(-0.07063000765580908, 0.03022518203580587),
                C64::new(0.012804145886109982, 0.02354142145759751),
                C64::new(-0.04314107866497867, -0.04016522703358564),
                C64::new(0.04314107866497867, 0.04016522703358564),
                C64::new(-0.021856858568925594, 0.029833367804835698),
            ],
        ),
    ] {
        let budget = FullWaveForceLimits {
            maximum_evaluations: 16_384,
            ..limits()
        };
        let owner = PreparedImpedanceWallForce::new(
            wall(Some(C64::new(0.001, 1.0))),
            HZ,
            [0.0, 0.0, -0.04],
            [radius, 0.0, -0.02],
            budget,
            &mut || false,
        )
        .unwrap();
        let x = owner.evaluate_reflected(unit(0)).unwrap();
        let y = owner.evaluate_reflected(unit(1)).unwrap();
        let z = owner.evaluate_reflected(unit(2)).unwrap();
        let actual = [
            x.pressure_pa,
            z.pressure_pa,
            x.particle_velocity_m_s[0],
            y.particle_velocity_m_s[1],
            z.particle_velocity_m_s[0],
            x.particle_velocity_m_s[2],
            z.particle_velocity_m_s[2],
        ];
        for (actual, expected) in actual.into_iter().zip(expected) {
            close(actual, expected, 2e-9, 1e-10);
        }
        close(
            z.particle_velocity_m_s[0],
            -x.particle_velocity_m_s[2],
            2e-13,
            1e-14,
        );
        assert!(!owner.diagnostics().exact_ideal_boundary);
        assert!(owner.diagnostics().evaluations <= budget.maximum_evaluations);
        if radius == 0.0 {
            assert_eq!(x.pressure_pa, C64::ZERO);
            assert!(x.particle_velocity_m_s[0].abs() > 0.1);
        }
    }
}

#[test]
fn force_cache_keeps_original_clock_self_field_and_bounded_refusals() {
    let plane = wall(Some(C64::new(0.3, 0.7)));
    let source = [0.0, 0.0, -0.04];
    let receiver = [0.12, 0.0, -0.03];
    let owner =
        PreparedImpedanceWallForce::new(plane, HZ, source, receiver, limits(), &mut || false)
            .unwrap();
    let force = [C64::new(0.2, -0.4), C64::new(0.1, 0.3), C64::new(-0.5, 0.2)];
    for (field, line) in [
        (
            owner.evaluate_reflected(force).unwrap(),
            owner.evaluate_reflected_harmonic_line(force).unwrap(),
        ),
        (
            owner.evaluate_total(force).unwrap(),
            owner.evaluate_total_harmonic_line(force).unwrap(),
        ),
    ] {
        assert_eq!(line.frequency_hz, HZ);
        assert_eq!(line.pressure_pa, field.pressure_pa);
        assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    }
    let self_field =
        PreparedImpedanceWallForce::new(plane, HZ, source, source, limits(), &mut || false)
            .unwrap();
    assert!(!self_field.total_field_available());
    assert!(self_field.evaluate_reflected(force).is_ok());
    assert!(self_field.evaluate_total(force).is_err());
    assert!(self_field.evaluate_total_harmonic_line(force).is_err());
    assert_eq!(
        components(&owner.evaluate_total([C64::ZERO; 3]).unwrap()),
        [C64::ZERO; 4]
    );
    assert!(
        owner
            .evaluate_reflected([C64::new(f64::NAN, 0.0); 3])
            .is_err()
    );
    assert_eq!(owner.source_position_m(), source);
    assert_eq!(owner.observer_position_m(), receiver);
    assert!(matches!(
        PreparedImpedanceWallForce::new(plane, HZ, source, receiver, limits(), &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let budget = FullWaveForceLimits {
        maximum_evaluations: 30,
        ..limits()
    };
    assert!(
        PreparedImpedanceWallForce::new(plane, HZ, source, receiver, budget, &mut || false)
            .is_err()
    );
}
