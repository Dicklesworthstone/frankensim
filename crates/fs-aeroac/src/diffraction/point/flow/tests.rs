//! G1/G3: independent convected source and linearized-Euler checks.
use super::*;
use crate::diffraction::{
    HalfPlaneDiffraction, HalfPlaneFace, HalfPlaneFrame, HalfPlaneMonopoleLimits,
    HalfPlaneReceiver, PreparedHalfPlaneMonopole,
};
use core::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 340.0).unwrap()
}
fn frame() -> HalfPlaneFrame {
    HalfPlaneFrame::new([0.0; 3], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap()
}
fn limits(k: f64) -> HalfPlaneMonopoleLimits {
    HalfPlaneMonopoleLimits {
        absolute_pressure_tolerance_pa_per_m3_s: 1.3 * 340.0 * k * 1e-12,
        absolute_velocity_tolerance_per_m2: 1e-11,
        relative_tolerance: 1e-10,
        maximum_evaluations: 100_000,
        maximum_depth: 42,
    }
}
fn prepared(
    source: [f64; 3],
    point: HalfPlaneReceiver,
    k: f64,
    flow: f64,
    boundary: IdealPlaneBoundary,
) -> PreparedHalfPlaneMonopole {
    PreparedHalfPlaneMonopole::with_edge_parallel_flow(
        medium(),
        k * 340.0 / TAU,
        frame(),
        boundary,
        source,
        point,
        flow,
        limits(k),
        || false,
    )
    .unwrap_or_else(|e| panic!("source={source:?}, point={point:?}, flow={flow}: {e:?}"))
}
fn values(owner: &PreparedHalfPlaneMonopole, q: C64) -> [C64; 4] {
    let line = owner.evaluate_harmonic_line(q).unwrap();
    [
        line.pressure_pa,
        line.particle_velocity_m_s[0],
        line.particle_velocity_m_s[1],
        line.particle_velocity_m_s[2],
    ]
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-10),
        "{a:?} != {b:?}"
    );
}

// Independent physical free-space convected Green function, Eq.(2.2).
// Differentiating this closed expression gives the source's actual p and u.
fn free(source: [f64; 3], point: [f64; 3], k: f64, flow: f64) -> [C64; 4] {
    let d: [f64; 3] = core::array::from_fn(|i| point[i] - source[i]);
    let m = flow / 340.0;
    let beta2 = 1.0 - m * m;
    let radius = (beta2 * (d[0] * d[0] + d[1] * d[1]) + d[2] * d[2]).sqrt();
    let angle = k * (radius - m * d[2]) / beta2;
    let green = C64::new(angle.cos(), angle.sin()).scale(1.0 / (4.0 * PI * radius));
    let gradient: [C64; 3] = core::array::from_fn(|i| {
        let dr = if i == 2 {
            d[i] / radius
        } else {
            beta2 * d[i] / radius
        };
        green
            * C64::new(
                -dr / radius,
                k * (dr - if i == 2 { m } else { 0.0 }) / beta2,
            )
    });
    let pressure = C64::new(0.0, -k * 340.0 * 1.3) * green + gradient[2].scale(1.3 * flow);
    [pressure, -gradient[0], -gradient[1], -gradient[2]]
}

#[test]
fn zero_flow_preserves_original_transfer_and_energy_uses_original_material() {
    let source = [-0.4, 0.3, 0.2];
    let point = HalfPlaneReceiver::WorldPosition([0.3, -0.2, 0.1]);
    let old = PreparedHalfPlaneMonopole::new(
        medium(),
        7.0 * 340.0 / TAU,
        frame(),
        IdealPlaneBoundary::Rigid,
        source,
        point,
        limits(7.0),
        || false,
    )
    .unwrap();
    let zero = prepared(source, point, 7.0, 0.0, IdealPlaneBoundary::Rigid);
    let q = C64::new(0.3, -0.7);
    assert_eq!(values(&old, q), values(&zero, q));
    assert_eq!(
        old.diagnostics().evaluations,
        zero.diagnostics().evaluations
    );
    for flow in [0.0, 0.6 * 340.0] {
        let owner = prepared(source, point, 7.0, flow, IdealPlaneBoundary::Rigid);
        assert_eq!(owner.acoustic(), medium());
        assert_eq!(owner.edge_flow_velocity_m_s(), flow);
        assert_eq!(owner.mean_flow_velocity_m_s(), [0.0, 0.0, flow]);
        let v = values(&owner, q);
        let energy = v[0].norm_sq() / (4.0 * medium().adiabatic_bulk_modulus_pa())
            + 0.25 * 1.3 * v[1..].iter().map(|u| u.norm_sq()).sum::<f64>();
        let result = owner.evaluate_energy_transport(q).unwrap();
        assert!((result.energy_density_j_m3 - energy).abs() < 3e-14 * energy);
        for i in 0..3 {
            let work = (v[0] * v[i + 1].conj()).re * 0.5;
            let flux = work + if i == 2 { flow * energy } else { 0.0 };
            assert!((result.pressure_work_flux_w_m2[i] - work).abs() < 3e-13 * work.abs().max(1.0));
            assert!((result.energy_flux_w_m2[i] - flux).abs() < 3e-13 * flux.abs().max(1.0));
        }
    }
}

#[test]
fn open_continuation_source_recovers_independent_convected_spherical_field() {
    let source = [-0.4, 0.0, 0.2];
    for mach in [-0.7, 0.35, 0.85] {
        for point in [[0.3, -0.2, 0.1], [-0.2, 0.4, 0.8], [-0.7, -0.1, -0.5]] {
            let owner = prepared(
                source,
                HalfPlaneReceiver::WorldPosition(point),
                7.0,
                mach * 340.0,
                IdealPlaneBoundary::Rigid,
            );
            for (a, b) in
                values(&owner, C64::ONE)
                    .into_iter()
                    .zip(free(source, point, 7.0, mach * 340.0))
            {
                close(a, b, 2e-9);
            }
        }
    }
}

#[test]
fn diffracted_field_closes_convected_momentum_continuity_and_energy_transport() {
    let source = [-0.4, 0.3, 0.2];
    let point = [0.3, -0.2, -0.1];
    let k = 7.0;
    let omega = k * 340.0;
    let q = C64::new(0.7, -0.2);
    for mach in [-0.6, 0.5] {
        let flow = mach * 340.0;
        for boundary in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let center = prepared(
                source,
                HalfPlaneReceiver::WorldPosition(point),
                k,
                flow,
                boundary,
            );
            let value = values(&center, q);
            let h = 2e-4;
            let mut derivative = [[C64::ZERO; 4]; 3];
            let mut second = [C64::ZERO; 3];
            let mut divergence_flux = 0.0;
            for axis in 0..3 {
                let samples = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
                    let mut at = point;
                    at[axis] += offset * h;
                    let owner = prepared(
                        source,
                        HalfPlaneReceiver::WorldPosition(at),
                        k,
                        flow,
                        boundary,
                    );
                    (
                        values(&owner, q),
                        owner.evaluate_energy_transport(q).unwrap().energy_flux_w_m2,
                    )
                });
                for component in 0..4 {
                    derivative[axis][component] = (samples[0].0[component]
                        - samples[1].0[component].scale(8.0)
                        + samples[2].0[component].scale(8.0)
                        - samples[3].0[component])
                        .scale(1.0 / (12.0 * h));
                }
                second[axis] = (-samples[0].0[0] + samples[1].0[0].scale(16.0)
                    - value[0].scale(30.0)
                    + samples[2].0[0].scale(16.0)
                    - samples[3].0[0])
                    .scale(1.0 / (12.0 * h * h));
                divergence_flux += (samples[0].1[axis] - 8.0 * samples[1].1[axis]
                    + 8.0 * samples[2].1[axis]
                    - samples[3].1[axis])
                    / (12.0 * h);
            }
            for i in 0..3 {
                let material =
                    C64::new(0.0, -omega) * value[i + 1] + derivative[2][i + 1].scale(flow);
                let residual = material.scale(1.3) + derivative[i][0];
                let scale = (1.3 * material.abs() + derivative[i][0].abs()).max(1.0);
                assert!(
                    residual.abs() < 2e-7 * scale,
                    "momentum {mach} {boundary:?}: {residual:?}"
                );
            }
            let pressure_rate = C64::new(0.0, -omega) * value[0] + derivative[2][0].scale(flow);
            let divergence = derivative[0][1] + derivative[1][2] + derivative[2][3];
            let residual =
                pressure_rate.scale(1.0 / medium().adiabatic_bulk_modulus_pa()) + divergence;
            let scale = (pressure_rate.abs() / medium().adiabatic_bulk_modulus_pa()
                + divergence.abs())
            .max(1e-6);
            assert!(
                residual.abs() < 2e-7 * scale,
                "continuity {mach} {boundary:?}: {residual:?}"
            );
            let wave = second[0]
                + second[1]
                + second[2].scale(1.0 - mach * mach)
                + C64::new(0.0, 2.0 * k * mach) * derivative[2][0]
                + value[0].scale(k * k);
            assert!(
                wave.abs() < 3e-5 * (k * k * value[0].abs()).max(1.0),
                "wave {mach} {boundary:?}: {wave:?}"
            );
            let flux = center
                .evaluate_energy_transport(q)
                .unwrap()
                .energy_flux_w_m2;
            let scale = k * flux.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!(
                divergence_flux.abs() < 2e-7 * scale.max(1e-6),
                "energy {divergence_flux:e}"
            );
        }
    }
}

#[test]
fn both_sheet_faces_and_flow_reversed_reciprocity_hold() {
    let source = [0.4, 0.3, 0.2];
    let point = [-0.15, -0.25, -0.4];
    for flow in [-210.0, 180.0] {
        for boundary in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            for face in [HalfPlaneFace::PositiveNormal, HalfPlaneFace::NegativeNormal] {
                let owner = prepared(
                    source,
                    HalfPlaneReceiver::SheetFace {
                        distance_from_edge_m: 0.3,
                        along_edge_m: 0.1,
                        face,
                    },
                    7.0,
                    flow,
                    boundary,
                );
                let v = values(&owner, C64::ONE);
                match boundary {
                    IdealPlaneBoundary::Rigid => assert_eq!(v[2], C64::ZERO),
                    IdealPlaneBoundary::PressureRelease => assert_eq!(v[0], C64::ZERO),
                }
            }
            let forward = prepared(
                source,
                HalfPlaneReceiver::WorldPosition(point),
                7.0,
                flow,
                boundary,
            );
            let reverse = prepared(
                point,
                HalfPlaneReceiver::WorldPosition(source),
                7.0,
                -flow,
                boundary,
            );
            close(
                values(&forward, C64::ONE)[0],
                values(&reverse, C64::ONE)[0],
                2e-9,
            );
            let wrong = prepared(
                point,
                HalfPlaneReceiver::WorldPosition(source),
                7.0,
                flow,
                boundary,
            );
            assert!((values(&forward, C64::ONE)[0] - values(&wrong, C64::ONE)[0]).abs() > 0.01);
        }
    }
}

#[test]
fn finite_sources_converge_to_original_flow_plane_wave_in_pressure_and_velocity() {
    let direction: [f64; 3] = [-0.48, -0.64, 0.6];
    let point = [0.3, -0.2, 0.1];
    let k = 7.0;
    for mach in [-0.45, 0.6] {
        let flow = mach * 340.0;
        let group = [direction[0], direction[1], direction[2] + mach];
        let norm = group[0].hypot(group[1]).hypot(group[2]);
        for boundary in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let reference = HalfPlaneDiffraction::with_edge_parallel_flow(
                medium(),
                k * 340.0 / TAU,
                direction,
                boundary,
                flow,
            )
            .unwrap()
            .evaluate_harmonic_line(point, C64::ONE, || false)
            .unwrap();
            let reference = [
                reference.pressure_pa,
                reference.particle_velocity_m_s[0],
                reference.particle_velocity_m_s[1],
                reference.particle_velocity_m_s[2],
            ];
            let mut errors = Vec::new();
            for distance in [100.0, 1000.0, 10000.0] {
                let source = group.map(|d| -distance * d / norm);
                let amplitude = free(source, [0.0; 3], k, flow)[0].recip();
                let owner = prepared(
                    source,
                    HalfPlaneReceiver::WorldPosition(point),
                    k,
                    flow,
                    boundary,
                );
                let v = values(&owner, amplitude);
                let error = (0..4)
                    .map(|i| (v[i] - reference[i]).abs() * if i == 0 { 1.0 } else { 1.3 * 340.0 })
                    .fold(0.0_f64, f64::max);
                errors.push(error);
            }
            assert!(
                errors[1] < 0.12 * errors[0] && errors[2] < 0.12 * errors[1],
                "{mach} {boundary:?} {errors:?}"
            );
            assert!(errors[2] < 5e-5, "{mach} {boundary:?} {errors:?}");
        }
    }
}

#[test]
fn flow_world_frame_source_clock_and_bounded_refusals_are_retained() {
    let source = [-0.4, 0.3, 0.2];
    let point = [0.3, -0.2, 0.1];
    let offset = [1.25, -2.0, 0.5];
    let world = |v: [f64; 3]| [offset[0] + v[2], offset[1] + v[0], offset[2] + v[1]];
    let placed = HalfPlaneFrame::new(offset, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let expected = prepared(
        source,
        HalfPlaneReceiver::WorldPosition(point),
        7.0,
        180.0,
        IdealPlaneBoundary::Rigid,
    );
    let actual = PreparedHalfPlaneMonopole::with_edge_parallel_flow(
        medium(),
        7.0 * 340.0 / TAU,
        placed,
        IdealPlaneBoundary::Rigid,
        world(source),
        HalfPlaneReceiver::WorldPosition(world(point)),
        180.0,
        limits(7.0),
        || false,
    )
    .unwrap();
    assert_eq!(actual.mean_flow_velocity_m_s(), [180.0, 0.0, 0.0]);
    assert_eq!(actual.frequency_hz(), 7.0 * 340.0 / TAU);
    let q = C64::new(0.3, -0.7);
    let a = values(&actual, q);
    let b = values(&expected, q);
    close(a[0], b[0], 3e-10);
    for (i, j) in [(1, 3), (2, 1), (3, 2)] {
        close(a[i], b[j], 3e-10);
    }
    let phase = C64::new(0.41_f64.cos(), 0.41_f64.sin());
    for (a, b) in values(&actual, q * phase).into_iter().zip(a) {
        close(a, b * phase, 3e-14);
    }
    let energy = actual.evaluate_energy_transport(q).unwrap();
    let original = expected.evaluate_energy_transport(q).unwrap();
    for (i, j) in [(0, 2), (1, 0), (2, 1)] {
        assert!(
            (energy.energy_flux_w_m2[i] - original.energy_flux_w_m2[j]).abs()
                < 1e-10 * original.energy_flux_w_m2[j].abs().max(1.0)
        );
    }
    let geometric = (0.5 + 0.3_f64.hypot(0.2)).hypot(0.1);
    assert!((actual.diagnostics().minimum_edge_path_m - geometric).abs() < 1e-14);
    for flow in [-340.0, 340.0, 341.0, f64::NAN, f64::INFINITY] {
        assert!(
            PreparedHalfPlaneMonopole::with_edge_parallel_flow(
                medium(),
                7.0 * 340.0 / TAU,
                frame(),
                IdealPlaneBoundary::Rigid,
                source,
                HalfPlaneReceiver::WorldPosition(point),
                flow,
                limits(7.0),
                || false
            )
            .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        PreparedHalfPlaneMonopole::with_edge_parallel_flow(
            medium(),
            7.0 * 340.0 / TAU,
            frame(),
            IdealPlaneBoundary::Rigid,
            source,
            HalfPlaneReceiver::WorldPosition(point),
            180.0,
            limits(7.0),
            || {
                polls += 1;
                polls == 20
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut bounded = limits(7.0);
    bounded.maximum_evaluations = 1;
    assert!(
        PreparedHalfPlaneMonopole::with_edge_parallel_flow(
            medium(),
            7.0 * 340.0 / TAU,
            frame(),
            IdealPlaneBoundary::Rigid,
            source,
            HalfPlaneReceiver::WorldPosition(point),
            180.0,
            bounded,
            || false
        )
        .is_err()
    );
}
