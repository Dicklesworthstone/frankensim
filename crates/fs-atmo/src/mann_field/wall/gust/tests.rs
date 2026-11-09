use super::*;
use crate::mann::MANN_TARGET_V1;
use crate::mann_field::fixed_step::FixedStepMannField;

const GAUSS: [(f64, f64); 8] = [
    (-0.9602898564975363, 0.1012285362903763),
    (-0.7966664774136267, 0.2223810344533745),
    (-0.525532409916329, 0.3137066458778873),
    (-0.1834346424956498, 0.362683783378362),
    (0.1834346424956498, 0.362683783378362),
    (0.525532409916329, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];

fn source() -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        45,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 12,
            min_kl: 0.1,
            max_kl: 4.0,
        },
        [3.4, -0.8, 0.0],
    )
    .unwrap();
    HomogeneousMannField::at_anchor(model, 87, vec![0.7; 12]).unwrap()
}

fn packet(source: &HomogeneousMannField) -> MannWallGust {
    MannWallGust::new(
        MannPlaneWall::new(0.0).unwrap(),
        [0.0, 0.0, 0.8],
        [2.0, 1.5, 2.0],
        source.tick() as f64 / TICK_HZ,
        1.3,
    )
    .unwrap()
}

fn norm(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().map(f64::abs).fold(0.0, f64::max)
}

fn compare(a: FlowSample, b: FlowSample, tolerance: f64) {
    for (left, right) in
        a.u.iter()
            .chain(a.grad.iter().flatten())
            .zip(b.u.iter().chain(b.grad.iter().flatten()))
    {
        assert!(
            (left - right).abs() <= tolerance * (1.0 + right.abs()),
            "{left} != {right}"
        );
    }
}

// Integrate the PUBLIC wall velocity numerically, then localize its potential.
// This oracle uses neither the production potential nor its envelope jet.
fn independent_potential(
    source: &HomogeneousMannField,
    gust: &MannWallGust,
    p: [f64; 3],
) -> [f64; 3] {
    let height = p[2] - gust.wall().height_m();
    let mut integrated = [0.0; 2];
    for (node, weight) in GAUSS {
        let q = [
            p[0],
            p[1],
            gust.wall().height_m() + 0.5 * height * (1.0 + node),
        ];
        let velocity = source.sample_wall_blocked(q, gust.wall()).unwrap().u;
        for i in 0..2 {
            integrated[i] += 0.5 * height * weight * velocity[i];
        }
    }
    let time = source.tick() as f64 / TICK_HZ;
    let mut gain = gust.strength();
    for i in 0..3 {
        let center =
            gust.center_m()[i] + source.model.advection_mps[i] * (time - gust.reference_time_s());
        let s = (p[i] - center) / gust.half_width_m()[i];
        gain *= if s.abs() < 1.0 {
            (1.0 - s * s).powi(3)
        } else {
            0.0
        };
    }
    [gain * integrated[1], -gain * integrated[0], 0.0]
}

#[test]
fn g0_packet_is_the_curl_of_an_independently_integrated_wall_potential() {
    let field = source();
    let gust = packet(&field);
    for p in [[0.2, -0.3, 0.7], [-0.7, 0.4, 1.1], [0.6, 0.3, 0.01]] {
        let h = 1e-5;
        let mut derivative = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut plus = p;
            let mut minus = p;
            plus[j] += h;
            minus[j] -= h;
            let a = independent_potential(&field, &gust, plus);
            let b = independent_potential(&field, &gust, minus);
            for i in 0..3 {
                derivative[i][j] = (a[i] - b[i]) / (2.0 * h);
            }
        }
        let expected = [
            derivative[2][1] - derivative[1][2],
            derivative[0][2] - derivative[2][0],
            derivative[1][0] - derivative[0][1],
        ];
        let actual = field.sample_wall_gust(p, &gust).unwrap();
        for i in 0..3 {
            assert!(
                (actual.u[i] - expected[i]).abs() < 2e-9,
                "component {i}: {:?} vs {expected:?}",
                actual.u
            );
        }
    }
}

#[test]
fn g0_all_gradients_and_continuity_hold_where_naive_velocity_modulation_fails() {
    let field = source();
    let gust = packet(&field);
    let mut naive_divergence: f64 = 0.0;
    let mut correction: f64 = 0.0;
    for p in [[0.2, -0.3, 0.7], [-0.7, 0.4, 1.1], [1.7, 0.6, 0.03]] {
        let flow = field.sample_wall_gust(p, &gust).unwrap();
        let scale = norm(flow.grad.into_iter().flatten()).max(1.0);
        assert!((flow.grad[0][0] + flow.grad[1][1] + flow.grad[2][2]).abs() < 3e-13 * scale);
        for j in 0..3 {
            let mut plus = p;
            let mut minus = p;
            plus[j] += 1e-5;
            minus[j] -= 1e-5;
            let a = field.sample_wall_gust(plus, &gust).unwrap();
            let b = field.sample_wall_gust(minus, &gust).unwrap();
            for i in 0..3 {
                let difference = (a.u[i] - b.u[i]) / 2e-5;
                assert!(
                    (difference - flow.grad[i][j]).abs() < 2e-9 * scale,
                    "gradient [{i}][{j}]: {difference} vs {}",
                    flow.grad[i][j]
                );
            }
        }
        let base = field.sample_wall_blocked(p, gust.wall()).unwrap();
        let envelope = gust
            .envelope(p, field.tick() as f64 / TICK_HZ, field.model.advection_mps)
            .unwrap()
            .unwrap();
        naive_divergence = naive_divergence.max(dot(envelope.grad, base.u).abs());
        correction = correction.max(norm((0..3).map(|i| flow.u[i] - envelope.value * base.u[i])));
    }
    assert!(
        naive_divergence > 0.01,
        "fixture must expose mass created by e*v"
    );
    assert!(
        correction > 0.01,
        "curl correction must affect the sampled velocity"
    );
}

#[test]
fn g0_closed_surface_flux_vanishes_without_using_the_analytic_gradient() {
    let field = source();
    let gust = packet(&field);
    let low = [-0.6, -0.5, 0.0];
    let high = [0.8, 0.4, 0.9];
    let mut flux = 0.0;
    let mut absolute = 0.0;
    for axis in 0..3 {
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        for side in 0..2 {
            for (x, wx) in GAUSS {
                for (y, wy) in GAUSS {
                    let mut p = [0.0; 3];
                    p[axis] = if side == 0 { low[axis] } else { high[axis] };
                    p[a] = 0.5 * (low[a] + high[a] + x * (high[a] - low[a]));
                    p[b] = 0.5 * (low[b] + high[b] + y * (high[b] - low[b]));
                    let v = field.sample_wall_gust(p, &gust).unwrap().u[axis];
                    let contribution = 0.25 * (high[a] - low[a]) * (high[b] - low[b]) * wx * wy * v;
                    flux += if side == 0 {
                        -contribution
                    } else {
                        contribution
                    };
                    absolute += contribution.abs();
                }
            }
        }
    }
    assert!(absolute > 0.1);
    assert!(
        flux.abs() < 2e-12 * absolute,
        "net flux {flux}, absolute flux {absolute}"
    );
}

#[test]
fn g0_wall_and_compact_support_have_exact_zeros_and_continuous_gradients() {
    let field = source();
    let gust = packet(&field);
    let zero = FlowSample {
        u: [0.0; 3],
        grad: [[0.0; 3]; 3],
    };
    for p in [[0.2, -0.3, 0.0], [-1.7, 0.4, 0.0], [0.0, 0.0, 0.0]] {
        let flow = field.sample_wall_gust(p, &gust).unwrap();
        assert_eq!(flow.u[2], 0.0);
        assert_eq!(flow.grad[2][0], 0.0);
        assert_eq!(flow.grad[2][1], 0.0);
    }
    for p in [
        [2.0, 0.0, 0.8],
        [3.0, 0.1, 0.8],
        [0.1, -1.5, 0.8],
        [0.1, 0.2, 3.0],
    ] {
        assert_eq!(field.sample_wall_gust(p, &gust).unwrap(), zero);
    }
    let coarse = field
        .sample_wall_gust([2.0 - 1e-4, 0.1, 0.8], &gust)
        .unwrap();
    let fine = field
        .sample_wall_gust([2.0 - 1e-5, 0.1, 0.8], &gust)
        .unwrap();
    assert!(norm(fine.u) < 0.02 * norm(coarse.u));
    assert!(norm(fine.grad.into_iter().flatten()) < 0.12 * norm(coarse.grad.into_iter().flatten()));
    let shifted = MannWallGust::new(
        MannPlaneWall::new(-2.5).unwrap(),
        [0.0, 0.0, -1.7],
        [2.0; 3],
        gust.reference_time_s(),
        1.0,
    )
    .unwrap();
    let at_wall = field.sample_wall_gust([0.2, 0.1, -2.5], &shifted).unwrap();
    assert_eq!(at_wall.u[2], 0.0);
    assert_eq!(at_wall.grad[2][0], 0.0);
    assert_eq!(at_wall.grad[2][1], 0.0);
}

#[test]
fn g0_vertical_integrals_handle_zero_wavenumbers_and_thin_wall_distances() {
    let wall = MannPlaneWall::new(0.0).unwrap();
    let modes = [
        SpectralMode {
            k: [1.0, 0.0, 0.0],
            weight: 1.0,
            polarization: [[0.0, 1.0, 0.5], [0.0, -0.2, 0.8]],
        },
        SpectralMode {
            k: [0.0, 0.0, 2.0],
            weight: 1.0,
            polarization: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        },
        SpectralMode {
            k: [1e-300, 0.0, 1.0],
            weight: 1.0,
            polarization: [[0.0, 1.0, 0.0], [-1.0, 0.0, 1e-300]],
        },
        SpectralMode {
            k: [0.5, -1.0, 1e-14],
            weight: 1.0,
            polarization: [[1.0, 0.5, 0.0], [-2e-14, 0.0, 1.0]],
        },
    ];
    for mode in modes {
        for height in [0.0, 1e-30, 1e-12, 0.7] {
            let p = [0.4, -0.2, height];
            let a = mode.polarization[0];
            let b = mode.polarization[1];
            let geometry = ModePoint::new(&mode, [0.0; 3], 0.0, p, wall, height).unwrap();
            let phase = dot(mode.k, [p[0], p[1], 0.0]);
            let actual = potential(&geometry, height, phase, a, b).unwrap();
            let mut integral = [0.0; 2];
            for (node, weight) in GAUSS {
                let h = 0.5 * height * (1.0 + node);
                let g = ModePoint::new(&mode, [0.0; 3], 0.0, [p[0], p[1], h], wall, h).unwrap();
                let v = g.values(a, b);
                for i in 0..2 {
                    integral[i] += 0.5 * height * weight * v[i];
                }
            }
            assert!((actual.u[0] - integral[1]).abs() <= 4e-15 * height.max(1e-300));
            assert!((actual.u[1] + integral[0]).abs() <= 4e-15 * height.max(1e-300));
            let curl = [
                -actual.grad[1][2],
                actual.grad[0][2],
                actual.grad[1][0] - actual.grad[0][1],
            ];
            let expected = geometry.values(a, b);
            for i in 0..3 {
                assert!((curl[i] - expected[i]).abs() < 2e-15);
            }
        }
    }
}

#[test]
fn g0_packet_inherits_source_advection_and_both_checkpointed_ou_clocks() {
    let mut source = source();
    let gust = packet(&source);
    let time = source.tick() as f64 / TICK_HZ;
    let p = [0.2, -0.3, 0.7];
    let amplitudes = core::array::from_fn(|i| source.paths[i].amplitudes());
    let before = sample_modes(source.model(), time, amplitudes, p, &gust).unwrap();
    let later = time + 0.25;
    let moved = core::array::from_fn(|i| p[i] + source.model.advection_mps[i] * 0.25);
    compare(
        before,
        sample_modes(source.model(), later, amplitudes, moved, &gust).unwrap(),
        2e-13,
    );
    let checkpoint = source.checkpoint();
    let mut restored = HomogeneousMannField::resume(checkpoint.clone()).unwrap();
    for point in [p, moved, [-0.4, 0.5, 0.8]] {
        assert_eq!(
            source.sample_wall_gust(point, &gust),
            restored.sample_wall_gust(point, &gust)
        );
    }
    assert_eq!(source.checkpoint(), checkpoint);
    source.advance_to(source.tick() + 7).unwrap();
    restored.advance_to(restored.tick() + 7).unwrap();
    assert_eq!(source.checkpoint(), restored.checkpoint());
    assert_ne!(before, source.sample_wall_gust(p, &gust).unwrap());
    let mut fixed = FixedStepMannField::from_snapshot(source.clone(), 991, 0.013, 88 * 12).unwrap();
    assert_eq!(
        fixed.sample_wall_gust(p, &gust),
        source.sample_wall_gust(p, &gust)
    );
    fixed.advance_to(3, 4 * 12 * 3, || false).unwrap();
    let checkpoint = fixed.checkpoint();
    let workspace = fixed.workspace_bytes();
    let mut restored = FixedStepMannField::resume(checkpoint.clone(), workspace).unwrap();
    for point in [moved, p, moved] {
        assert_eq!(
            fixed.sample_wall_gust(point, &gust),
            restored.sample_wall_gust(point, &gust)
        );
    }
    assert_eq!(fixed.checkpoint(), checkpoint);
    assert_eq!(fixed.workspace_bytes(), workspace);
    fixed.advance_to(5, 4 * 12 * 2, || false).unwrap();
    restored.advance_to(5, 4 * 12 * 2, || false).unwrap();
    assert_eq!(fixed.checkpoint(), restored.checkpoint());
    assert_eq!(
        fixed.sample_wall_gust(p, &gust),
        restored.sample_wall_gust(p, &gust)
    );
}

#[test]
fn g0_strength_is_linear_at_all_points_and_center_velocity_retains_its_source() {
    let field = source();
    let one = MannWallGust::new(
        MannPlaneWall::new(0.0).unwrap(),
        [0.0, 0.0, 0.8],
        [2.0; 3],
        field.tick() as f64 / TICK_HZ,
        1.0,
    )
    .unwrap();
    let two = MannWallGust::new(
        one.wall(),
        one.center_m(),
        one.half_width_m(),
        one.reference_time_s(),
        2.0,
    )
    .unwrap();
    let zero = MannWallGust::new(
        one.wall(),
        one.center_m(),
        one.half_width_m(),
        one.reference_time_s(),
        0.0,
    )
    .unwrap();
    for p in [one.center_m(), [0.3, 0.4, 0.1], [-0.6, 0.9, 1.3]] {
        let a = field.sample_wall_gust(p, &one).unwrap();
        let b = field.sample_wall_gust(p, &two).unwrap();
        for (a, b) in
            a.u.iter()
                .chain(a.grad.iter().flatten())
                .zip(b.u.iter().chain(b.grad.iter().flatten()))
        {
            assert_eq!(2.0 * a, *b);
        }
        assert_eq!(
            field.sample_wall_gust(p, &zero).unwrap(),
            FlowSample {
                u: [0.0; 3],
                grad: [[0.0; 3]; 3]
            }
        );
    }
    assert_eq!(
        field.sample_wall_gust(one.center_m(), &one).unwrap().u,
        field
            .sample_wall_blocked(one.center_m(), one.wall())
            .unwrap()
            .u
    );
}

#[test]
fn g0_refusals_do_not_change_either_source_and_outside_support_skips_spectral_work() {
    let source = source();
    let gust = packet(&source);
    for width in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-300, 1e300] {
        assert!(MannWallGust::new(gust.wall(), [0.0; 3], [width, 1.0, 1.0], 0.0, 1.0).is_err());
    }
    for strength in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(MannWallGust::new(gust.wall(), [0.0; 3], [1.0; 3], 0.0, strength).is_err());
    }
    assert!(MannWallGust::new(gust.wall(), [f64::NAN, 0.0, 1.0], [1.0; 3], 0.0, 1.0).is_err());
    assert!(MannWallGust::new(gust.wall(), [0.0; 3], [1.0; 3], f64::NAN, 1.0).is_err());
    let fixed = FixedStepMannField::from_snapshot(source.clone(), 991, 0.013, 88 * 12).unwrap();
    let original = source.checkpoint();
    let fixed_original = fixed.checkpoint();
    for point in [
        [0.0, 0.0, -0.01],
        [f64::NAN, 0.0, 1.0],
        [0.0, f64::INFINITY, 1.0],
    ] {
        assert!(source.sample_wall_gust(point, &gust).is_err());
        assert!(fixed.sample_wall_gust(point, &gust).is_err());
    }
    let huge = MannWallGust::new(
        gust.wall(),
        gust.center_m(),
        [0.1; 3],
        gust.reference_time_s(),
        1e308,
    )
    .unwrap();
    assert!(source.sample_wall_gust(huge.center_m(), &huge).is_err());
    assert!(fixed.sample_wall_gust(huge.center_m(), &huge).is_err());
    let ancient =
        MannWallGust::new(gust.wall(), gust.center_m(), [1.0; 3], -f64::MAX, 1.0).unwrap();
    assert!(source.sample_wall_gust([0.0, 0.0, 1.0], &ancient).is_err());
    assert_eq!(source.checkpoint(), original);
    assert_eq!(fixed.checkpoint(), fixed_original);
    let invalid_amplitudes = [f64::NAN; 12];
    assert_eq!(
        sample_modes(
            source.model(),
            gust.reference_time_s(),
            [&invalid_amplitudes; 4],
            [100.0, 0.0, 1.0],
            &gust
        )
        .unwrap(),
        FlowSample {
            u: [0.0; 3],
            grad: [[0.0; 3]; 3]
        }
    );
    assert!(source.sample_wall_gust([0.3, 0.2, 0.8], &gust).is_ok());
    assert!(fixed.sample_wall_gust([0.3, 0.2, 0.8], &gust).is_ok());
}
