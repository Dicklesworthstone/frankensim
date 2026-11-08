//! G1/G3: impedance-boundary physics, ideal Green limits and reciprocity.

use fs_aeroac::curle3d::refraction::full_wave::FullWaveMonopoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    PlanarImpedanceWall, PreparedImpedanceWallMonopole,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn fluid() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 343.0).unwrap()
}
fn limits() -> FullWaveMonopoleLimits {
    FullWaveMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1e-7,
        absolute_velocity_tolerance_per_m2: 1e-10,
        relative_tolerance: 1e-10,
        maximum_evaluations: 262_144,
        maximum_depth: 40,
    }
}
fn wall(normalized: Option<C64>) -> PlanarImpedanceWall {
    PlanarImpedanceWall::new(
        fluid(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        normalized.map(|z| z.scale(1.2 * 343.0)),
    )
    .unwrap()
}
fn prepare(
    normalized: Option<C64>,
    source: [f64; 3],
    receiver: [f64; 3],
) -> PreparedImpedanceWallMonopole {
    PreparedImpedanceWallMonopole::new(
        wall(normalized),
        700.0,
        source,
        receiver,
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn close(a: C64, b: C64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * b.abs(),
        "actual={a:?}, expected={b:?}"
    );
}

// Independent outgoing spherical source, using ordinary complex closed forms.
fn green(source: [f64; 3], receiver: [f64; 3]) -> [C64; 4] {
    let d: [f64; 3] = core::array::from_fn(|i| receiver[i] - source[i]);
    let radius = d.iter().map(|x| x * x).sum::<f64>().sqrt();
    let omega = TAU * 700.0;
    let k = omega / 343.0;
    let (sin, cos) = (k * radius).sin_cos();
    let g = C64::new(cos, sin).scale(1.0 / (4.0 * PI * radius));
    let p = C64::new(0.0, -omega * 1.2) * g;
    let u = C64::new(1.0 / radius, -k) * g;
    [
        p,
        u.scale(d[0] / radius),
        u.scale(d[1] / radius),
        u.scale(d[2] / radius),
    ]
}

#[test]
fn ideal_limits_recover_independent_spherical_images_and_wall_conditions() {
    let source = [0.02, -0.01, -0.12];
    let receiver = [0.14, 0.05, -0.08];
    let image = [source[0], source[1], -source[2]];
    let reflected = green(image, receiver);
    let direct = green(source, receiver);
    for (impedance, sign) in [(None, 1.0), (Some(C64::ZERO), -1.0)] {
        let owner = prepare(impedance, source, receiver);
        assert!(owner.diagnostics().exact_ideal_boundary);
        assert_eq!(owner.diagnostics().evaluations, 0);
        let value = owner.evaluate_reflected_harmonic_line(C64::ONE).unwrap();
        let value = [
            value.pressure_pa,
            value.particle_velocity_m_s[0],
            value.particle_velocity_m_s[1],
            value.particle_velocity_m_s[2],
        ];
        let total = owner.evaluate_total_harmonic_line(C64::ONE).unwrap();
        let total = [
            total.pressure_pa,
            total.particle_velocity_m_s[0],
            total.particle_velocity_m_s[1],
            total.particle_velocity_m_s[2],
        ];
        for i in 0..4 {
            close(value[i], reflected[i].scale(sign), 3e-13, 1e-11);
            close(total[i], direct[i] + reflected[i].scale(sign), 3e-13, 1e-11);
        }
        let on_wall = prepare(impedance, source, [0.14, 0.05, 0.0])
            .evaluate_total(C64::ONE)
            .unwrap();
        if sign > 0.0 {
            assert!(on_wall.particle_velocity_m_s[2].abs() < 1e-13);
        } else {
            assert!(on_wall.pressure_pa.abs() < 1e-10);
        }
    }
}

#[test]
fn full_spectrum_satisfies_robin_condition_including_lossy_surface_mode() {
    let source = [0.0, 0.0, -0.07];
    // The compliant wall has evanescent amplification >1 despite passive resistance.
    for z in [
        C64::from_re(0.3),
        C64::ONE,
        C64::new(4.0, -3.0),
        C64::new(0.15, 2.0),
        C64::new(0.0, -1.5),
        C64::new(0.0, 0.2),
        C64::new(0.0, 1.0),
        C64::new(0.0, 2.0),
    ] {
        for radius in [0.0, 0.18] {
            let owner = prepare(Some(z), source, [radius, 0.0, 0.0]);
            assert!(!owner.diagnostics().exact_ideal_boundary);
            assert!(owner.diagnostics().evaluations > 0);
            let value = owner.evaluate_total(C64::new(0.002, -0.001)).unwrap();
            close(
                value.pressure_pa,
                z.scale(1.2 * 343.0) * value.particle_velocity_m_s[2],
                2e-8,
                2e-9,
            );
            // Outward normal intensity is wall absorption, with no source of work.
            assert!(value.active_intensity_w_m2[2] >= -1e-9);
            if z.re == 0.0 {
                assert!(value.active_intensity_w_m2[2].abs() < 1e-9);
            } else {
                assert!(value.active_intensity_w_m2[2] > 0.0);
            }
        }
    }
}

#[test]
fn narrow_passive_surface_mode_matches_independent_fields_with_bounded_work() {
    let source = [0.0, 0.0, -0.04];
    let z = C64::new(0.001, 1.0);
    let mut budget = limits();
    budget.maximum_evaluations = 16_384;
    // This wall's evanescent amplitude peaks near 2000. On axis the Weyl
    // integral has the independent closed form p=-i*omega*rho*F/(4*pi),
    // u_n=F'/(4*pi), where s=i*k/z and
    // F(H)=exp(i*k*H)/H + 2*s*exp(-s*H)*E1(-H*(s+i*k)).
    // Principal-branch E1 was evaluated independently of the production
    // spectrum. The off-axis reference uses separate adaptive quadrature in
    // the physical vertical wavenumbers, split around the complex pole.
    for (radius, expected) in [
        (
            0.0,
            [
                C64::new(11060.970430630407, 1780.0211985817225),
                C64::ZERO,
                C64::new(-35.5638177360633, -41.87343034999833),
            ],
        ),
        (
            0.17,
            [
                C64::new(-6093.672850897635, 4054.6385935549797),
                C64::new(-21.122278195304883, 11.488377282009878),
                C64::new(12.820503691755377, 8.769174824978263),
            ],
        ),
    ] {
        let owner = PreparedImpedanceWallMonopole::new(
            wall(Some(z)),
            700.0,
            source,
            [radius, 0.0, -0.02],
            budget,
            &mut || false,
        )
        .unwrap();
        let reflected = owner.evaluate_reflected(C64::ONE).unwrap();
        close(reflected.pressure_pa, expected[0], 2e-9, 1e-7);
        close(reflected.particle_velocity_m_s[0], expected[1], 2e-9, 1e-10);
        close(reflected.particle_velocity_m_s[2], expected[2], 2e-9, 1e-10);
        assert_eq!(reflected.particle_velocity_m_s[1], C64::ZERO);
        assert!(owner.diagnostics().evaluations <= budget.maximum_evaluations);
    }
}

#[test]
fn lossless_surface_wave_matches_independent_principal_value_and_residue() {
    // Independent 90-digit arithmetic: for a=k/X, H=.06, the evanescent
    // integrals are I0=1/H+2*a*I and I1=1/H²+2*a/H+2*a²*I, where
    // I=-exp(-a*H)*Ei(a*H)+i*pi*exp(-a*H). Their pressure and outward
    // velocity are -i*omega*rho*I0/(4*pi) and -I1/(4*pi). The remaining
    // propagating integral was evaluated with an independent entire power
    // series and rational moments, without the production quadrature.
    for (reactance, pressure, normal_velocity) in [
        (
            0.2,
            C64::new(-1650.1417773107182, 12307.677951028231),
            C64::new(60.7157079516139, -42.272360402093314),
        ),
        (
            1.0,
            C64::new(11070.152457546734, 1778.4778983626435),
            C64::new(-35.60627954605893, -41.887885006570245),
        ),
        (
            2.0,
            C64::new(9009.845676865571, -2592.240233489048),
            C64::new(-36.9677355464601, -20.02170778055919),
        ),
    ] {
        let budget = FullWaveMonopoleLimits {
            maximum_evaluations: 16_384,
            ..limits()
        };
        let owner = PreparedImpedanceWallMonopole::new(
            wall(Some(C64::new(0.0, reactance))),
            700.0,
            [0.0, 0.0, -0.04],
            [0.0, 0.0, -0.02],
            budget,
            &mut || false,
        )
        .unwrap();
        let field = owner.evaluate_reflected(C64::ONE).unwrap();
        close(field.pressure_pa, pressure, 2e-11, 1e-7);
        close(
            field.particle_velocity_m_s[2],
            normal_velocity,
            2e-11,
            1e-10,
        );
        assert_eq!(field.particle_velocity_m_s[..2], [C64::ZERO; 2]);
        assert!(owner.diagnostics().evaluations <= budget.maximum_evaluations);
        // Three original 15-node regions plus the four cached numerator and
        // analytic-derivative samples; every subdivision adds 30 samples.
        assert_eq!(owner.diagnostics().evaluations % 30, 19);
    }
}

#[test]
fn positive_resistance_tends_to_the_same_outgoing_lossless_field() {
    let source = [0.0, 0.0, -0.04];
    for radius in [0.0, 0.17] {
        let receiver = [radius, 0.0, -0.02];
        let limit = prepare(Some(C64::new(0.0, 1.0)), source, receiver)
            .evaluate_reflected(C64::ONE)
            .unwrap();
        let limit = [
            limit.pressure_pa,
            limit.particle_velocity_m_s[0],
            limit.particle_velocity_m_s[1],
            limit.particle_velocity_m_s[2],
        ];
        let mut previous = [f64::INFINITY; 4];
        for resistance in [0.01, 0.001, 0.0001] {
            let value = prepare(Some(C64::new(resistance, 1.0)), source, receiver)
                .evaluate_reflected(C64::ONE)
                .unwrap();
            let value = [
                value.pressure_pa,
                value.particle_velocity_m_s[0],
                value.particle_velocity_m_s[1],
                value.particle_velocity_m_s[2],
            ];
            for i in 0..4 {
                let error = (value[i] - limit[i]).abs();
                assert!(error <= 0.15 * previous[i] + 1e-8);
                previous[i] = error;
                if resistance == 0.0001 {
                    close(value[i], limit[i], 5e-4, 1e-8);
                }
            }
        }
    }
}

#[test]
fn matched_normal_incidence_retains_spherical_near_field_and_reciprocity() {
    let source = [0.0, 0.0, -0.018];
    let receiver = [0.0, 0.0, -0.026];
    let owner = prepare(Some(C64::ONE), source, receiver);
    // A normal-incidence image multiplier would incorrectly return zero here.
    let reflected = owner.evaluate_reflected(C64::ONE).unwrap();
    assert!(reflected.pressure_pa.abs() > 100.0);
    let a = [0.03, -0.01, -0.06];
    let b = [0.19, 0.08, -0.09];
    for z in [C64::new(0.2, 1.3), C64::new(0.0, 1.3)] {
        let forward = prepare(Some(z), a, b).evaluate_total(C64::ONE).unwrap();
        let reverse = prepare(Some(z), b, a).evaluate_total(C64::ONE).unwrap();
        close(forward.pressure_pa, reverse.pressure_pa, 2e-10, 1e-7);
        let coincident = prepare(Some(z), a, a);
        assert!(!coincident.total_field_available());
        assert!(coincident.evaluate_reflected(C64::ONE).is_ok());
        assert!(coincident.evaluate_total(C64::ONE).is_err());
    }
}

#[test]
fn refined_quadrature_converges_without_changing_the_wall_physics() {
    let source = [0.0, 0.0, -0.04];
    let receiver = [0.17, 0.0, -0.02];
    for resistance in [0.1, 0.0] {
        let mut values = Vec::new();
        let mut evidence = Vec::new();
        for relative in [1e-5, 1e-7, 1e-9] {
            let mut budget = limits();
            budget.relative_tolerance = relative;
            budget.absolute_pressure_tolerance_pa_per_m3_s = relative;
            budget.absolute_velocity_tolerance_per_m2 = relative * 1e-3;
            let owner = PreparedImpedanceWallMonopole::new(
                wall(Some(C64::new(resistance, 1.0))),
                700.0,
                source,
                receiver,
                budget,
                &mut || false,
            )
            .unwrap();
            values.push(owner.evaluate_reflected(C64::ONE).unwrap().pressure_pa);
            evidence.push(
                owner.diagnostics().pressure_quadrature_error_pa_per_m3_s
                    + owner.diagnostics().pressure_tail_bound_pa_per_m3_s,
            );
        }
        assert!((values[1] - values[2]).abs() <= evidence[1] + evidence[2] + 1e-7);
        assert!((values[0] - values[2]).abs() <= evidence[0] + evidence[2] + 1e-7);
        assert!(evidence[2] < evidence[0]);
    }
}

#[test]
fn cancellation_work_and_unsupported_wall_refusals_publish_no_response() {
    for z in [C64::new(-0.1, 0.0), C64::new(f64::NAN, 0.0)] {
        assert!(PlanarImpedanceWall::new(fluid(), [0.0; 3], [0.0, 0.0, 1.0], Some(z)).is_err());
    }
    let source = [0.0, 0.0, -0.1];
    let receiver = [0.2, 0.0, -0.1];
    // Pure compliance is admissible, but its four cached pole evaluations
    // and subsequent quadrature still obey the caller's original work limit.
    for maximum_evaluations in [4, 49] {
        let budget = FullWaveMonopoleLimits {
            maximum_evaluations,
            ..limits()
        };
        assert!(
            PreparedImpedanceWallMonopole::new(
                wall(Some(C64::new(0.0, 1.0))),
                700.0,
                source,
                receiver,
                budget,
                &mut || false,
            )
            .is_err()
        );
    }
    let mut polls = 0;
    assert!(
        PreparedImpedanceWallMonopole::new(
            wall(Some(C64::new(0.0, 1.0))),
            700.0,
            source,
            receiver,
            limits(),
            &mut || {
                polls += 1;
                polls == 8
            },
        )
        .is_err()
    );
    assert!(
        PreparedImpedanceWallMonopole::new(
            wall(Some(C64::ONE)),
            700.0,
            source,
            receiver,
            limits(),
            &mut || true
        )
        .is_err()
    );
    let mut budget = limits();
    budget.maximum_evaluations = 15;
    assert!(
        PreparedImpedanceWallMonopole::new(
            wall(Some(C64::ONE)),
            700.0,
            source,
            receiver,
            budget,
            &mut || false
        )
        .is_err()
    );
    assert!(
        PreparedImpedanceWallMonopole::new(
            wall(Some(C64::ONE)),
            700.0,
            source,
            [0.0, 0.0, 0.1],
            limits(),
            &mut || false
        )
        .is_err()
    );
}
