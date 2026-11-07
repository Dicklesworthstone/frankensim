//! G1/G3: exact Hertz sliding/spin integrals against independent spatial work.

use core::f64::consts::PI;

use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

#[path = "sliding_patch/clipped.rs"]
mod clipped;

fn patch(normal: f64, radius: f64, mu: f64) -> HertzSlidingPatch {
    HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: normal,
        contact_radius_m: radius,
        friction_coefficient: mu,
    })
    .unwrap()
}

fn rate_grid(radial_panels: usize, angular_panels: usize) -> HertzFrictionQuadrature {
    HertzFrictionQuadrature {
        radial_panels,
        angular_panels,
    }
}

fn rate_laws() -> [FrictionLaw; 4] {
    [
        FrictionLaw::VelocityDependent {
            static_mu: 0.85,
            mu_zero: 0.65,
            slope_per_speed: -0.9,
        },
        FrictionLaw::Stribeck {
            static_mu: 0.85,
            kinetic_mu: 0.16,
            characteristic_speed: 0.3,
            viscous_per_speed: 0.08,
        },
        FrictionLaw::Hyperbolic {
            static_mu: 0.85,
            kinetic_mu: 0.16,
            characteristic_speed: 0.2,
        },
        FrictionLaw::ShiftedPower {
            static_mu: 0.85,
            kinetic_mu: 0.16,
            characteristic_speed: 0.2,
            exponent: 0.4,
        },
    ]
}

/// Independent full-disk midpoint integration of original local traction and
/// heat. It uses neither the production Coulomb subtraction nor its Gauss rule.
fn spatial_rate_law(
    law: &FrictionLaw,
    geometry: [f64; 2],
    slip: [f64; 2],
    spin: f64,
    grid: [usize; 2],
) -> [f64; 4] {
    let [normal, radius] = geometry;
    let [radial, angular] = grid;
    let mut output = [0.0; 4];
    let dtheta = (PI / 2.0) / radial as f64;
    for i in 0..radial {
        let theta = (i as f64 + 0.5) * dtheta;
        let (s, c) = theta.sin_cos();
        let r = radius * s;
        let load = 3.0 * normal * s * c * c * dtheta / angular as f64;
        for j in 0..angular {
            let phi = 2.0 * PI * (j as f64 + 0.5) / angular as f64;
            let (s, c) = phi.sin_cos();
            let position = [r * c, r * s];
            let velocity = [slip[0] - spin * position[1], slip[1] + spin * position[0]];
            let speed = velocity[0].hypot(velocity[1]);
            if speed > 0.0 {
                let magnitude = load * law.kinetic_coefficient(speed).unwrap();
                let local = velocity.map(|v| -magnitude * (v / speed));
                output[0] += local[0];
                output[1] += local[1];
                output[2] += position[0] * local[1] - position[1] * local[0];
                output[3] += magnitude * speed;
            }
        }
    }
    output
}

#[test]
fn shared_curves_keep_exact_coulomb_and_pure_translation_limits() {
    let original = patch(7.0, 0.12, 0.3);
    let coulomb = HertzFrictionPatch::new(
        7.0,
        0.12,
        FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.3,
        },
    )
    .unwrap();
    for (slip, spin) in [
        ([0.3, -0.4], 4.0),
        ([0.0; 2], -4.0),
        ([0.3, -0.4], 0.0),
        ([0.0; 2], 0.0),
    ] {
        assert_eq!(
            coulomb.evaluate(slip, spin, rate_grid(1, 1)).unwrap(),
            original.evaluate(slip, spin).unwrap()
        );
    }
    for law in rate_laws() {
        let model = HertzFrictionPatch::new(7.0, 0.12, law.clone()).unwrap();
        assert_eq!(model.friction_law(), &law);
        let expected = patch(7.0, 0.12, law.kinetic_coefficient(0.5).unwrap())
            .evaluate([0.3, -0.4], 0.0)
            .unwrap();
        assert_eq!(
            model.evaluate([0.3, -0.4], 0.0, rate_grid(1, 1)).unwrap(),
            expected
        );
    }
}

#[test]
fn affine_local_speed_law_has_exact_hertz_viscous_force_torque_and_heat() {
    let normal = 7.0;
    let radius = 0.12;
    let slip = [0.3, -0.4];
    let spin = -4.0;
    for (mu_zero, slope) in [(0.3, 0.4), (0.3, -0.1), (0.0, 0.4)] {
        let law = FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero,
            slope_per_speed: slope,
        };
        let model = HertzFrictionPatch::new(normal, radius, law.clone()).unwrap();
        let base = patch(normal, radius, mu_zero).evaluate(slip, spin).unwrap();
        let result = model.evaluate(slip, spin, rate_grid(1, 1)).unwrap();
        for i in 0..2 {
            close(
                result.force_n[i],
                base.force_n[i] - normal * slope * slip[i],
                2e-15,
            );
        }
        let polar_moment = (2.0 / 5.0) * normal * radius * radius;
        close(
            result.spin_moment_nm,
            base.spin_moment_nm - slope * polar_moment * spin,
            2e-16,
        );
        close(
            result.dissipated_power_w,
            base.dissipated_power_w + slope * (normal * 0.25 + polar_moment * spin * spin),
            3e-15,
        );
        close(result.power_residual_w, 0.0, 2e-15);
        assert_eq!(
            result,
            model.evaluate(slip, spin, rate_grid(256, 256)).unwrap()
        );
        let reference = spatial_rate_law(&law, [normal, radius], slip, spin, [256, 512]);
        close(result.force_n[0], reference[0], 3e-5);
        close(result.force_n[1], reference[1], 3e-5);
        close(result.spin_moment_nm, reference[2], 3e-6);
        close(result.dissipated_power_w, reference[3], 3e-5);
    }
}

#[test]
fn all_shared_curves_resolve_local_sliding_spin_and_heat_over_one_patch() {
    for law in rate_laws() {
        let model = HertzFrictionPatch::new(3.2, 0.4, law.clone()).unwrap();
        for slip in [[0.36, -0.48], [0.0; 2]] {
            let spin = -2.5;
            let result = model.evaluate(slip, spin, rate_grid(32, 32)).unwrap();
            let reference = spatial_rate_law(&law, [3.2, 0.4], slip, spin, [512, 1024]);
            close(result.force_n[0], reference[0], 3e-5);
            close(result.force_n[1], reference[1], 3e-5);
            close(result.spin_moment_nm, reference[2], 1e-5);
            close(result.dissipated_power_w, reference[3], 3e-5);
            close(result.power_residual_w, 0.0, 3e-14);
            assert!(result.dissipated_power_w > 0.0);
            if slip == [0.0; 2] {
                assert_eq!(result.force_n, [0.0; 2]);
                // Zero centre speed does not make a spinning patch stick.
                assert!(result.spin_moment_nm > 0.0);
                let frozen = patch(3.2, 0.4, law.zero_speed_kinetic_coefficient().unwrap())
                    .evaluate(slip, spin)
                    .unwrap();
                assert!((result.spin_moment_nm - frozen.spin_moment_nm).abs() > 0.03);
            }
        }
    }
}

#[test]
fn nonaffine_patch_quadrature_has_live_three_rung_refinement() {
    let law = FrictionLaw::Stribeck {
        static_mu: 0.85,
        kinetic_mu: 0.16,
        characteristic_speed: 0.3,
        viscous_per_speed: 0.08,
    };
    let model = HertzFrictionPatch::new(1.0, 1.0, law).unwrap();
    let slip = [0.6, 0.0];
    let reference = model.evaluate(slip, 1.0, rate_grid(128, 128)).unwrap();
    let mut errors = Vec::new();
    for panels in [2, 4, 8] {
        let result = model
            .evaluate(slip, 1.0, rate_grid(panels, panels))
            .unwrap();
        errors.push(
            (result.force_n[0] - reference.force_n[0])
                .abs()
                .max((result.spin_moment_nm - reference.spin_moment_nm).abs())
                .max((result.dissipated_power_w - reference.dissipated_power_w).abs()),
        );
        close(result.power_residual_w, 0.0, 2e-15);
    }
    eprintln!("Local-speed Hertz patch refinement errors: {errors:?}");
    assert!(errors[1] < errors[0] / 4.0, "{errors:?}");
    assert!(errors[2] < errors[1] / 4.0, "{errors:?}");
    assert!(errors[2] < 2e-5, "{errors:?}");

    // Steep velocity weakening can feed translation energy into spin while
    // local friction and total generalized work remain dissipative. For weak
    // spin, M = -N a^2 omega [mu(U)/U + mu'(U)]/5 + O(omega^3).
    let curve = FrictionLaw::Stribeck {
        static_mu: 0.8,
        kinetic_mu: 0.0,
        characteristic_speed: 0.3,
        viscous_per_speed: 0.0,
    };
    let (mu, derivative) = curve.kinetic_coefficient_and_slope(0.3).unwrap();
    let weakening = HertzFrictionPatch::new(2.0, 0.1, curve)
        .unwrap()
        .evaluate([0.3, 0.0], 0.1, rate_grid(8, 8))
        .unwrap();
    let moment = -2.0 * 0.1 * 0.1 * 0.1 * (mu / 0.3 + derivative) / 5.0;
    assert!(weakening.spin_moment_nm > 0.0);
    close(weakening.spin_moment_nm, moment, 2e-6);
    assert!(weakening.dissipated_power_w > 0.0);
    close(weakening.power_residual_w, 0.0, 2e-15);
}

#[test]
fn rate_dependent_patch_is_objective_reversible_and_deterministic() {
    let slip = [0.3, -0.4];
    let spin = 4.0;
    let angle: f64 = 0.713;
    let (s, c) = angle.sin_cos();
    let rotate = |v: [f64; 2]| [c * v[0] - s * v[1], s * v[0] + c * v[1]];
    for law in rate_laws() {
        let model = HertzFrictionPatch::new(7.0, 0.12, law).unwrap();
        let grid = rate_grid(16, 16);
        let base = model.evaluate(slip, spin, grid).unwrap();
        assert_eq!(base, model.evaluate(slip, spin, grid).unwrap());
        let rotated = model.evaluate(rotate(slip), spin, grid).unwrap();
        for (actual, expected) in rotated.force_n.into_iter().zip(rotate(base.force_n)) {
            close(actual, expected, 3e-15);
        }
        close(rotated.spin_moment_nm, base.spin_moment_nm, 3e-16);
        close(rotated.dissipated_power_w, base.dissipated_power_w, 3e-15);
        let reversed = model.evaluate(slip.map(|v| -v), -spin, grid).unwrap();
        assert_eq!(reversed.force_n, base.force_n.map(|f| -f));
        assert_eq!(reversed.spin_moment_nm, -base.spin_moment_nm);
        assert_eq!(reversed.dissipated_power_w, base.dissipated_power_w);
        let reflected = model.evaluate([slip[0], -slip[1]], -spin, grid).unwrap();
        assert_eq!(reflected.force_n, [base.force_n[0], -base.force_n[1]]);
        assert_eq!(reflected.spin_moment_nm, -base.spin_moment_nm);
    }
}

#[test]
fn rate_dependent_patch_preserves_rest_lift_off_clamping_and_refusals() {
    let law = rate_laws()[0].clone();
    let model = HertzFrictionPatch::new(3.2, 0.4, law.clone()).unwrap();
    let rest = model.evaluate([0.0; 2], 0.0, rate_grid(1, 1)).unwrap();
    assert!(rest.at_rest);
    assert_eq!(rest.dissipated_power_w, 0.0);
    let no_load = HertzFrictionPatch::new(0.0, 0.4, law)
        .unwrap()
        .evaluate([0.1, 0.2], 2.5, rate_grid(1, 1))
        .unwrap();
    let clamped = model.evaluate([3.0, 4.0], 2.5, rate_grid(1, 1)).unwrap();
    for response in [no_load, clamped] {
        assert!(!response.at_rest);
        assert_eq!(response.force_n, [0.0; 2]);
        assert_eq!(response.spin_moment_nm, 0.0);
        assert_eq!(response.dissipated_power_w, 0.0);
    }
    // A fast-spinning clipped law still dissipates in a very small central
    // core. Even the first node of a 256-panel full-disk grid misses this core.
    // Integrating that core directly retains its positive, representable heat.
    let spin = 10_000.0;
    let cutoff: f64 = 0.65 / 0.9;
    let edge_speed = 0.4 * spin;
    let ratio = cutoff / edge_speed;
    // Independent small-core pressure expansion; the next relative term is
    // 3*ratio^4/112, below 3e-17 for these physical inputs.
    let expected_heat =
        3.2 * 0.65 * cutoff.powi(3) / (4.0 * edge_speed * edge_speed) * (1.0 - ratio * ratio / 5.0);
    for panels in [1, 8, 256] {
        let core = model
            .evaluate([0.0; 2], spin, rate_grid(panels, 1))
            .unwrap();
        close(core.dissipated_power_w, expected_heat, 2e-22);
        close(core.spin_moment_nm, -expected_heat / spin, 2e-26);
        close(core.power_residual_w, 0.0, 2e-22);
    }
    let moving_core = model.evaluate([0.1, 0.0], spin, rate_grid(1, 4)).unwrap();
    assert!(moving_core.force_n[0] < 0.0);
    assert!(moving_core.dissipated_power_w > 0.0);
    close(moving_core.dissipated_power_w, expected_heat, 1e-16);
    close(moving_core.power_residual_w, 0.0, 2e-22);
    for grid in [
        rate_grid(0, 1),
        rate_grid(1, 0),
        rate_grid(257, 1),
        rate_grid(1, usize::MAX),
    ] {
        assert!(model.evaluate([0.0; 2], 0.0, grid).is_err());
    }
    assert!(
        model
            .evaluate([f64::NAN, 0.0], 0.0, rate_grid(1, 1))
            .is_err()
    );
    assert!(
        model
            .evaluate([0.0; 2], f64::INFINITY, rate_grid(1, 1))
            .is_err()
    );
    assert!(HertzFrictionPatch::new(-1.0, 0.4, rate_laws()[0].clone()).is_err());
    assert!(HertzFrictionPatch::new(1.0, 0.0, rate_laws()[0].clone()).is_err());
    assert!(
        HertzFrictionPatch::new(
            1.0,
            0.4,
            FrictionLaw::Coulomb {
                static_mu: 0.2,
                kinetic_mu: 0.8
            }
        )
        .is_err()
    );
}

#[test]
fn shrinking_off_centre_friction_disks_keep_force_spin_transfer_and_heat() {
    // G1: independent Taylor integration of p about the stationary point.
    // For k=v_cut/(a*omega) << 1 and H=sqrt(1-d^2):
    // F=-N*mu*d*k^3/(8H), D=N*mu*v_cut*k^2*H/4,
    // M=-N*mu*a*k^3*(2-3d^2)/(8H), with relative O(k^2) error.
    // At d>sqrt(2/3), translation feeds spin while total work stays passive.
    let normal = 3.2;
    let radius = 0.4;
    let mu = 0.65;
    let cutoff = mu / 0.9;
    let k = 1e-6;
    let edge = cutoff / k;
    let spin = edge / radius;
    let model = HertzFrictionPatch::new(normal, radius, rate_laws()[0].clone()).unwrap();
    for d in [0.1_f64, 0.5, 0.7, 0.9] {
        let h = (1.0 - d * d).sqrt();
        let force = -normal * mu * d * k.powi(3) / (8.0 * h);
        let heat = normal * mu * cutoff * k * k * h / 4.0;
        let moment = -normal * mu * radius * k.powi(3) * (2.0 - 3.0 * d * d) / (8.0 * h);
        let response = model
            .evaluate([d * edge, 0.0], spin, rate_grid(2, 8))
            .unwrap();
        close(response.force_n[0], force, force.abs() * 1e-8);
        assert_eq!(response.force_n[1], 0.0);
        close(response.spin_moment_nm, moment, moment.abs() * 1e-8);
        close(response.dissipated_power_w, heat, heat * 1e-8);
        close(response.power_residual_w, 0.0, heat * 2e-14);
        assert_eq!(response.spin_moment_nm > 0.0, d > (2.0_f64 / 3.0).sqrt());
        assert_eq!(
            response,
            model
                .evaluate([d * edge, 0.0], spin, rate_grid(2, 8))
                .unwrap()
        );
    }
}

#[test]
fn contained_active_disk_agrees_with_independent_full_patch_and_frame_changes() {
    // G1/G3: ordinary-size disks, including tangency to the Hertz boundary,
    // checked in the original patch coordinates with an unrelated rule.
    for (d, cutoff) in [(0.35, 0.3), (0.6, 0.4)] {
        let law = FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.6,
            slope_per_speed: -0.6 / cutoff,
        };
        let model = HertzFrictionPatch::new(2.0, 0.4, law.clone()).unwrap();
        let slip = [0.6 * d, -0.8 * d];
        let response = model.evaluate(slip, 2.5, rate_grid(24, 24)).unwrap();
        let reference = spatial_rate_law(&law, [2.0, 0.4], slip, 2.5, [512, 1024]);
        close(response.force_n[0], reference[0], 2e-6);
        close(response.force_n[1], reference[1], 2e-6);
        close(response.spin_moment_nm, reference[2], 1e-6);
        close(response.dissipated_power_w, reference[3], 2e-6);
        close(response.power_residual_w, 0.0, 2e-16);
        let rotated = model
            .evaluate([-slip[1], slip[0]], 2.5, rate_grid(24, 24))
            .unwrap();
        assert_eq!(rotated.force_n, [-response.force_n[1], response.force_n[0]]);
        assert_eq!(rotated.spin_moment_nm, response.spin_moment_nm);
        assert_eq!(rotated.dissipated_power_w, response.dissipated_power_w);
        let reversed = model
            .evaluate(slip.map(|v| -v), -2.5, rate_grid(24, 24))
            .unwrap();
        assert_eq!(reversed.force_n, response.force_n.map(|f| -f));
        assert_eq!(reversed.spin_moment_nm, -response.spin_moment_nm);
        assert_eq!(reversed.dissipated_power_w, response.dissipated_power_w);
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}

/// Independent midpoint surface quadrature; the production owner uses no grid.
/// r=a*sin(theta) turns p(r)*r*dr into 3P sin(theta) cos(theta)^2 dtheta.
/// Force, moment, and local-slip heat are accumulated independently.
fn spatial(
    parameters: HertzSlidingPatchParameters,
    slip: [f64; 2],
    spin: f64,
    radial: usize,
    angular: usize,
) -> [f64; 4] {
    let mut output = [0.0; 4];
    let dr = (PI / 2.0) / radial as f64;
    for i in 0..radial {
        let theta = (i as f64 + 0.5) * dr;
        let (sine, cosine) = theta.sin_cos();
        let radius = parameters.contact_radius_m * sine;
        let force = 3.0
            * parameters.normal_force_n
            * parameters.friction_coefficient
            * sine
            * cosine
            * cosine
            * dr
            / angular as f64;
        for j in 0..angular {
            let angle = 2.0 * PI * (j as f64 + 0.5) / angular as f64;
            let (s, c) = angle.sin_cos();
            let x = radius * c;
            let y = radius * s;
            let velocity = [slip[0] - spin * y, slip[1] + spin * x];
            let speed = velocity[0].hypot(velocity[1]);
            if speed > 0.0 {
                let local = velocity.map(|v| -force * v / speed);
                output[0] += local[0];
                output[1] += local[1];
                output[2] += x * local[1] - y * local[0];
                output[3] += force * speed;
            }
        }
    }
    output
}

#[test]
fn pure_translation_and_spin_use_pressure_first_moment_not_rms_radius() {
    let law = patch(12.0, 0.02, 0.4);
    let translate = law.evaluate([3.0, -4.0], 0.0).unwrap();
    close(translate.force_n[0], -2.88, 2e-15);
    close(translate.force_n[1], 3.84, 2e-15);
    assert_eq!(translate.spin_moment_nm, 0.0);
    close(translate.dissipated_power_w, 24.0, 8e-15);
    let expected = (3.0 * PI / 16.0) * 4.8 * 0.02;
    for spin in [-17.0, 17.0] {
        let result = law.evaluate([0.0; 2], spin).unwrap();
        assert_eq!(result.force_n, [0.0; 2]);
        close(result.spin_moment_nm, -expected.copysign(spin), 3e-17);
        close(result.dissipated_power_w, expected * spin.abs(), 3e-16);
        close(result.power_residual_w, 0.0, 3e-16);
    }
    assert!((expected - 0.4_f64.sqrt() * 4.8 * 0.02).abs() > 1e-3);
}

#[test]
fn sliding_and_spinning_share_one_friction_capacity() {
    let law = patch(1.0, 1.0, 1.0);
    let result = law.evaluate([1.0, 0.0], 1.0).unwrap();
    close(result.force_n[0], -9.0 * PI / 32.0, 4e-16);
    close(result.spin_moment_nm, -9.0 * PI / 128.0, 2e-16);
    assert!(result.force_n[0].abs() < law.sliding_capacity_n());
    assert!(result.spin_moment_nm.abs() < law.spinning_capacity_nm());
    let faster_translation = law.evaluate([10.0, 0.0], 1.0).unwrap();
    assert!(faster_translation.force_n[0].abs() > result.force_n[0].abs());
    assert!(faster_translation.spin_moment_nm.abs() < result.spin_moment_nm.abs());
}

#[test]
fn exact_wrench_matches_independent_hertz_surface_integrals() {
    let law = patch(3.2, 0.4, 0.7);
    for epsilon in [0.2, 0.75, 1.0, 1.5, 2.5, 20.0] {
        let slip = [0.6 * epsilon, -0.8 * epsilon];
        let result = law.evaluate(slip, 2.5).unwrap();
        let reference = spatial(law.parameters(), slip, 2.5, 384, 1024);
        close(result.force_n[0], reference[0], 2.0e-5);
        close(result.force_n[1], reference[1], 2.0e-5);
        close(result.spin_moment_nm, reference[2], 8.0e-6);
        close(result.dissipated_power_w, reference[3], 1.0e-4);
        close(result.power_residual_w, 0.0, 2e-14);
        assert!(result.dissipated_power_w > 0.0);
    }
}

#[test]
fn independent_surface_reference_has_a_live_three_rung_refinement() {
    let law = patch(1.0, 1.0, 1.0);
    let slip = [1.5, -2.0];
    let result = law.evaluate(slip, 1.0).unwrap();
    let exact = [
        result.force_n[0],
        result.force_n[1],
        result.spin_moment_nm,
        result.dissipated_power_w,
    ];
    let mut errors = Vec::new();
    for radial in [16, 32, 64] {
        let reference = spatial(law.parameters(), slip, 1.0, radial, 128);
        let error = reference
            .into_iter()
            .zip(exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(error);
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order > 1.95 && order < 2.05,
            "radial errors {errors:?}, order {order}"
        );
    }
    // The angular axis is independently live near the interior stagnation point.
    let exact = law.evaluate([0.6, 0.0], 1.0).unwrap();
    let mut angular_errors = Vec::new();
    for angular in [16, 32, 64] {
        let reference = spatial(law.parameters(), [0.6, 0.0], 1.0, 512, angular);
        angular_errors.push((reference[0] - exact.force_n[0]).abs());
    }
    assert!(angular_errors[0] > angular_errors[1]);
    assert!(angular_errors[1] > angular_errors[2]);
    assert!(angular_errors[2] < angular_errors[0] / 8.0);
    eprintln!(
        "Hertz surface reference radial errors {errors:?}, angular errors {angular_errors:?}"
    );
}

#[test]
fn arbitrary_slip_directions_rotation_reflection_and_rate_scaling_are_covariant() {
    let law = patch(7.0, 0.12, 0.3);
    let slip = [0.3, -0.4];
    let spin = 4.0;
    let base = law.evaluate(slip, spin).unwrap();
    let angle: f64 = 0.713;
    let (s, c) = angle.sin_cos();
    let rotate = |v: [f64; 2]| [c * v[0] - s * v[1], s * v[0] + c * v[1]];
    let rotated = law.evaluate(rotate(slip), spin).unwrap();
    for (actual, expected) in rotated.force_n.into_iter().zip(rotate(base.force_n)) {
        close(actual, expected, 2e-15);
    }
    close(rotated.spin_moment_nm, base.spin_moment_nm, 2e-16);
    close(rotated.dissipated_power_w, base.dissipated_power_w, 2e-15);
    let reflected = law.evaluate([slip[0], -slip[1]], -spin).unwrap();
    assert_eq!(reflected.force_n, [base.force_n[0], -base.force_n[1]]);
    assert_eq!(reflected.spin_moment_nm, -base.spin_moment_nm);
    for scale in [1e-200, 1e-20, 1.0, 1e20, 1e200] {
        let scaled = law.evaluate(slip.map(|x| scale * x), scale * spin).unwrap();
        for (actual, expected) in scaled.force_n.into_iter().zip(base.force_n) {
            close(actual, expected, 2e-15);
        }
        close(scaled.spin_moment_nm, base.spin_moment_nm, 2e-16);
        close(
            scaled.dissipated_power_w / scale,
            base.dissipated_power_w,
            2e-15,
        );
    }
    assert_eq!(base, law.evaluate(slip, spin).unwrap());
}

#[test]
fn vanishing_spin_remains_resolved_without_large_ratio_cancellation() {
    let law = patch(1.0, 1.0, 1.0);
    for z in [1e-4, 1e-8, 1e-100, 1e-300] {
        let result = law.evaluate([1.0, 0.0], z).unwrap();
        // Independent leading moments of the Hertz pressure: <y^2>=a^2/5.
        close(-result.spin_moment_nm / z, 0.2, 2e-9);
        close(-result.force_n[0], 1.0 - z * z / 10.0, 4e-16);
    }
    for epsilon in [1.0, 2.0] {
        let before = law.evaluate([epsilon * (1.0 - 1e-9), 0.0], 1.0).unwrap();
        let after = law.evaluate([epsilon * (1.0 + 1e-9), 0.0], 1.0).unwrap();
        close(before.force_n[0], after.force_n[0], 2e-9);
        close(before.spin_moment_nm, after.spin_moment_nm, 2e-9);
    }
}

#[test]
fn rest_frictionless_limits_and_invalid_inputs_are_explicit() {
    let law = patch(1.0, 1.0, 0.5);
    let rest = law.evaluate([0.0; 2], 0.0).unwrap();
    assert!(rest.at_rest);
    assert_eq!(rest.force_n, [0.0; 2]);
    assert_eq!(rest.spin_moment_nm, 0.0);
    assert_eq!(rest.dissipated_power_w, 0.0);
    for no_friction in [patch(0.0, 1.0, 0.5), patch(1.0, 1.0, 0.0)] {
        let response = no_friction.evaluate([3.0, 4.0], 7.0).unwrap();
        assert!(!response.at_rest);
        assert_eq!(response.force_n, [0.0; 2]);
        assert_eq!(response.spin_moment_nm, 0.0);
        assert_eq!(response.dissipated_power_w, 0.0);
    }
    for (normal, radius, mu) in [
        (-1.0, 1.0, 0.5),
        (1.0, 0.0, 0.5),
        (1.0, 1.0, -0.5),
        (f64::INFINITY, 1.0, 0.5),
    ] {
        assert!(
            HertzSlidingPatch::new(HertzSlidingPatchParameters {
                normal_force_n: normal,
                contact_radius_m: radius,
                friction_coefficient: mu,
            })
            .is_err()
        );
    }
    assert!(law.evaluate([f64::NAN, 0.0], 0.0).is_err());
    assert!(law.evaluate([0.0; 2], f64::INFINITY).is_err());
    let enormous = patch(1e200, 1.0, 1.0);
    assert!(enormous.evaluate([1e200, 0.0], 0.0).is_err());
}
