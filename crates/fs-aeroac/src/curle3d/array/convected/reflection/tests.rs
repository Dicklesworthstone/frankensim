//! G1: independent convected-potential images, wall work and source coherence.

use super::*;
use crate::curle3d::{CompactSourcePhasor, IdealPlaneBoundary, PreparedCompactSourceObserver};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.7, 310.0).unwrap()
}
fn prepared(
    flow: [f64; 3],
    point: [f64; 3],
    sources: &[[f64; 3]],
) -> PreparedConvectedMonopoleObserver {
    PreparedConvectedMonopoleObserver::new(
        medium(),
        flow,
        117.0,
        point,
        sources,
        sources.len(),
        ConvectedMonopoleLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 1e-6,
        },
    )
    .unwrap()
}
fn plane(axis: usize, boundary: IdealPlaneBoundary) -> IdealAcousticPlane {
    let mut normal = [0.0; 3];
    normal[axis] = 1.0;
    IdealAcousticPlane::new([0.0; 3], normal, boundary).unwrap()
}
fn close(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() < 3e-12 * b.abs().max(scale),
        "actual={a:?}, expected={b:?}"
    );
}
// Independent differentiated convected Green potential; no prepared owner,
// image geometry helper or moving-source pressure coefficients are used.
fn analytic(flow: [f64; 3], point: [f64; 3], source: [f64; 3], q: C64) -> [C64; 4] {
    let rho = 1.7;
    let c = 310.0;
    let omega = core::f64::consts::TAU * 117.0;
    let m = flow.map(|v| v / c);
    let r = core::array::from_fn::<_, 3, _>(|i| point[i] - source[i]);
    let mdotr: f64 = (0..3).map(|i| m[i] * r[i]).sum();
    let beta2 = 1.0 - m.iter().map(|v| v * v).sum::<f64>();
    let range = (beta2 * r.iter().map(|v| v * v).sum::<f64>() + mdotr * mdotr).sqrt();
    let d = (range - mdotr) / beta2;
    let gr = core::array::from_fn::<_, 3, _>(|i| (beta2 * r[i] + m[i] * mdotr) / range);
    let gd = core::array::from_fn::<_, 3, _>(|i| (gr[i] - m[i]) / beta2);
    let phase = omega * d / c;
    let gq =
        C64::new(phase.cos(), phase.sin()).scale(1.0 / (4.0 * core::f64::consts::PI * range)) * q;
    let conv_r: f64 = (0..3).map(|i| flow[i] * gr[i]).sum();
    let conv_d: f64 = (0..3).map(|i| flow[i] * gd[i]).sum();
    let mut result = [C64::ZERO; 4];
    result[0] = gq * C64::new(-rho * conv_r / range, rho * (-omega + omega * conv_d / c));
    for i in 0..3 {
        result[i + 1] = gq * C64::new(gr[i] / range, -omega * gd[i] / c);
    }
    result
}

#[test]
fn complete_corner_images_match_independent_convected_pressure_velocity_and_energy() {
    let flow = [0.0, 0.0, 124.0];
    let point = [0.8, 0.9, -0.4];
    let sources = [[0.3, 0.4, 0.2], [0.6, 0.2, -0.1]];
    let flows = [C64::new(0.002, -0.003), C64::new(-0.001, 0.004)];
    for b0 in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        for b1 in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let planes = [plane(0, b0), plane(1, b1)];
            let original = prepared(flow, point, &sources);
            let direct_delays = [original.source_delay_s(0), original.source_delay_s(1)];
            let reflected = original.with_orthogonal_planes(&planes, 8).unwrap();
            let result = reflected.evaluate(&flows).unwrap();
            let mut expected = [C64::ZERO; 4];
            for (&source, &q) in sources.iter().zip(&flows) {
                for mask in 0..4 {
                    let mut image = source;
                    let mut sign = 1.0;
                    for i in 0..2 {
                        if mask & (1 << i) != 0 {
                            image[i] = -image[i];
                            if planes[i].boundary() == IdealPlaneBoundary::PressureRelease {
                                sign = -sign;
                            }
                        }
                    }
                    let field = analytic(flow, point, image, q.scale(sign));
                    for i in 0..4 {
                        expected[i] = expected[i] + field[i];
                    }
                }
            }
            close(result.field.pressure_pa, expected[0], 1e-3);
            for i in 0..3 {
                close(result.field.particle_velocity_m_s[i], expected[i + 1], 1e-6);
            }
            let energy = expected[0].norm_sq() / (4.0 * 1.7 * 310.0 * 310.0)
                + 1.7 * expected[1..].iter().map(|v| v.norm_sq()).sum::<f64>() / 4.0;
            assert!((result.mean_perturbation_energy_density_j_m3 - energy).abs() < 1e-12 * energy);
            for i in 0..3 {
                let work = (expected[0] * expected[i + 1].conj()).re * 0.5;
                let expected_flux = work + flow[i] * energy;
                assert!(
                    (result.mean_perturbation_energy_flux_w_m2[i] - expected_flux).abs()
                        < 1e-12 * expected_flux.abs().max(1e-8)
                );
            }
            assert_eq!(reflected.acoustic(), medium());
            assert_eq!(reflected.flow_velocity_m_s(), flow);
            assert_eq!(reflected.reflection_planes(), &planes);
            assert_eq!(reflected.source_count(), 2);
            assert_eq!(
                [reflected.source_delay_s(0), reflected.source_delay_s(1)],
                direct_delays
            );
            assert_eq!(reflected.source_position_m(1), Some(sources[1]));
        }
    }
}

#[test]
fn rigid_and_release_walls_have_zero_normal_energy_flux_in_tangential_flow() {
    let flow = [85.0, 170.0, 0.0];
    let point = [0.4, -0.7, 0.2];
    let source = [1.0, -1.0, 0.8]; // point + 0.3*[2,-1,2]
    let normal = [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let wall = IdealAcousticPlane::new(point, [2.0, -1.0, 2.0], boundary).unwrap();
        let owner = prepared(flow, point, &[source])
            .with_orthogonal_planes(&[wall], 2)
            .unwrap();
        let result = owner.evaluate(&[C64::new(0.004, -0.003)]).unwrap();
        let un = (0..3).fold(C64::ZERO, |sum, i| {
            sum + result.field.particle_velocity_m_s[i].scale(normal[i])
        });
        match boundary {
            IdealPlaneBoundary::Rigid => assert!(un.abs() < 1e-15),
            IdealPlaneBoundary::PressureRelease => assert!(result.field.pressure_pa.abs() < 1e-12),
        }
        let flux: f64 = (0..3)
            .map(|i| result.mean_perturbation_energy_flux_w_m2[i] * normal[i])
            .sum();
        assert!(flux.abs() < 1e-14);
    }
    let walls = [
        plane(0, IdealPlaneBoundary::Rigid),
        plane(1, IdealPlaneBoundary::PressureRelease),
    ];
    for point in [[0.0, 0.8, 0.4], [0.7, 0.0, 0.4], [0.0, 0.0, 0.4]] {
        let field = prepared([0.0, 0.0, 124.0], point, &[[0.3, 0.5, 0.2]])
            .with_orthogonal_planes(&walls, 4)
            .unwrap()
            .evaluate(&[C64::ONE])
            .unwrap();
        if point[0] == 0.0 {
            assert!(field.field.particle_velocity_m_s[0].abs() < 1e-13);
        }
        if point[1] == 0.0 {
            assert!(field.field.pressure_pa.abs() < 1e-10);
        }
        for axis in 0..2 {
            if point[axis] == 0.0 {
                assert!(field.mean_perturbation_energy_flux_w_m2[axis].abs() < 1e-10);
            }
        }
    }
}

#[test]
fn quiescent_three_wall_limit_matches_original_compact_owner() {
    let walls = [
        plane(0, IdealPlaneBoundary::Rigid),
        plane(1, IdealPlaneBoundary::PressureRelease),
        plane(2, IdealPlaneBoundary::Rigid),
    ];
    let point = [0.9, 0.7, 0.6];
    let sources = [[0.2, 0.3, 0.4], [0.6, 0.1, 0.2]];
    let flows = [C64::new(0.003, 0.002), C64::new(-0.001, 0.001)];
    let convected = prepared([0.0; 3], point, &sources)
        .with_orthogonal_planes(&walls, 16)
        .unwrap();
    let original = PreparedCompactSourceObserver::new(medium(), 117.0, point, &sources, 2)
        .unwrap()
        .with_orthogonal_planes(&walls, 16)
        .unwrap();
    let expected = original
        .evaluate(&flows.map(CompactSourcePhasor::monopole))
        .unwrap();
    let actual = convected.evaluate(&flows).unwrap();
    close(actual.field.pressure_pa, expected.pressure_pa, 1e-4);
    for i in 0..3 {
        close(
            actual.field.particle_velocity_m_s[i],
            expected.particle_velocity_m_s[i],
            1e-7,
        );
        assert_eq!(
            actual.mean_perturbation_energy_flux_w_m2[i],
            actual.field.active_intensity_w_m2[i]
        );
    }
}

#[test]
fn image_folding_preserves_each_physical_sources_coherence_and_clock() {
    let walls = [
        plane(0, IdealPlaneBoundary::Rigid),
        plane(1, IdealPlaneBoundary::PressureRelease),
    ];
    let owner = prepared(
        [0.0, 0.0, 124.0],
        [0.8, 0.7, -0.3],
        &[[0.2, 0.3, 0.4], [0.5, 0.1, 0.2]],
    )
    .with_orthogonal_planes(&walls, 8)
    .unwrap();
    let a = [C64::new(0.003, 0.002), C64::ZERO];
    let b = [C64::ZERO, C64::new(-0.002, 0.004)];
    let budget = CoherenceBudget {
        max_modes: 2,
        max_source_mode_pairs: 4,
    };
    let statistics = owner
        .evaluate_coherence_modes(&[&a, &b], budget, || false)
        .unwrap();
    let fa = owner.evaluate(&a).unwrap();
    let fb = owner.evaluate(&b).unwrap();
    let expected = 0.5 * (fa.field.pressure_pa.norm_sq() + fb.field.pressure_pa.norm_sq());
    assert!((statistics.field.mean_square_pressure_pa2 - expected).abs() < 1e-13 * expected);
    assert_eq!(statistics.field.scope, owner.scope());
    let coherent = [a[0], b[1]];
    let line = owner.evaluate_harmonic_line(&coherent).unwrap();
    let total = owner.evaluate(&coherent).unwrap();
    assert_eq!(line.frequency_hz, 117.0);
    assert_eq!(line.pressure_pa, total.field.pressure_pa);
    assert_eq!(
        line.particle_velocity_m_s,
        total.field.particle_velocity_m_s
    );
    assert!((0.5 * total.field.pressure_pa.norm_sq() - expected).abs() > 0.01 * expected);
}

#[test]
fn boundary_admission_refuses_normal_flow_invalid_groups_and_partial_work() {
    let source = [[0.3, 0.4, 0.5]];
    let point = [0.7, 0.8, 0.9];
    let x = plane(0, IdealPlaneBoundary::Rigid);
    let y = plane(1, IdealPlaneBoundary::Rigid);
    assert!(
        prepared([1e-100, 0.0, 0.0], point, &source)
            .with_orthogonal_planes(&[x], 2)
            .is_err()
    );
    assert!(
        prepared([2.0, 0.0, 124.0], point, &source)
            .with_orthogonal_planes(&[x], 2)
            .is_err()
    );
    let free = prepared([0.0, 0.0, 124.0], point, &source);
    assert!(free.clone().with_orthogonal_planes(&[], 8).is_err());
    assert!(free.clone().with_orthogonal_planes(&[x, x], 4).is_err());
    assert!(free.clone().with_orthogonal_planes(&[x, y], 3).is_err());
    assert!(
        prepared([0.0, 0.0, 124.0], [-0.1, 0.8, 0.9], &source)
            .with_orthogonal_planes(&[x], 2)
            .is_err()
    );
    assert!(
        prepared([0.0, 0.0, 124.0], point, &[[0.0, 0.4, 0.5]])
            .with_orthogonal_planes(&[x], 2)
            .is_err()
    );
    let bound = free.clone().with_orthogonal_planes(&[x], 2).unwrap();
    assert!(bound.with_orthogonal_planes(&[y], 4).is_err());
    let mut polls = 0;
    let result = free.with_orthogonal_planes_cancellable(&[x, y], 4, || {
        polls += 1;
        polls >= 15
    });
    assert!(matches!(result, Err(AeroacError::Cancelled { .. })));
}
