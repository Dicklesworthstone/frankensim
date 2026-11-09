use super::*;
use crate::mann::MANN_TARGET_V1;
use crate::mann_field::{
    HomogeneousMannField, MannQuadrature, MannQuadratureRule, MannSpectralModel,
    fixed_step::FixedStepMannField,
};

const EDGES_0: [f64; 5] = [-1.0, -0.73, -0.12, 0.41, 1.0];
const EDGES_1: [f64; 4] = [-1.0, -0.47, 0.36, 1.0];

fn source(seed: u64) -> HomogeneousMannField {
    let model = MannSpectralModel::build_with_rule(
        seed,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 16,
            min_kl: 0.04,
            max_kl: 12.0,
        },
        [7.8, -0.9, 0.2],
        MannQuadratureRule::BalancedEnergyV1,
    )
    .unwrap();
    let mut source = HomogeneousMannField::at_anchor(model, 173, vec![0.7; 16]).unwrap();
    source.advance_to(source.tick() + 7).unwrap();
    source
}

fn patch() -> MannInflowPatch {
    MannInflowPatch::new([0.3, -0.7, 4.0], [[0.2, 0.7, 0.1], [-0.1, 0.15, 0.65]]).unwrap()
}

fn frame(source: &HomogeneousMannField) -> MannInflowFrame {
    source.inflow_frame(patch(), 168 * 16, || false).unwrap()
}

fn plane(frame: &MannInflowFrame, density: MannInflowDensity) -> MannInflowPlane {
    MannInflowPlane::prepare(
        frame,
        [&EDGES_0, &EDGES_1],
        density,
        MannInflowPlane::required_workspace_bytes(16, [4, 3]).unwrap(),
        || false,
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} vs {b:e}"
    );
}

// Independent numerical surface integration of the physical point field.
// Neither modal gains nor the production exponential-average kernel enters.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    let mut result = Vec::new();
    for i in 0..n {
        let mut x = (core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 0.0;
        for _ in 0..30 {
            let (mut a, mut b) = (1.0, x);
            for j in 2..=n {
                let c = ((2 * j - 1) as f64 * x * b - (j - 1) as f64 * a) / j as f64;
                a = b;
                b = c;
            }
            derivative = n as f64 * (x * b - a) / (x * x - 1.0);
            let delta = b / derivative;
            x -= delta;
            if delta.abs() < 2e-16 {
                break;
            }
        }
        result.push((x, 1.0 / ((1.0 - x * x) * derivative * derivative)));
    }
    result
}

fn integrate<const N: usize>(
    face: MannInflowFace,
    mut value: impl FnMut([f64; 3]) -> [f64; N],
) -> [f64; N] {
    let rule = gauss(40);
    let mut result = [0.0; N];
    for &(a, wa) in &rule {
        for &(b, wb) in &rule {
            let point = core::array::from_fn(|d| {
                face.center_m[d] + a * face.half_edges_m[0][d] + b * face.half_edges_m[1][d]
            });
            let sample = value(point);
            for d in 0..N {
                result[d] += wa * wb * sample[d];
            }
        }
    }
    result
}

#[test]
fn nonuniform_skew_faces_have_exact_geometry_order_and_area_averaged_velocity() {
    let source = source(43);
    let frame = frame(&source);
    let plane = plane(&frame, MannInflowDensity::ConstantKgM3(1.19));
    assert_eq!(plane.dimensions(), [4, 3]);
    assert_eq!(plane.samples().len(), 12);
    assert_eq!(plane.flat_index([4, 0]), None);
    assert_eq!(plane.flat_index([0, 3]), None);
    assert_eq!(plane.face(12), None);
    let patch = patch();
    let mut area = 0.0;
    let mut center_bias: f64 = 0.0;
    for j in 0..3 {
        for i in 0..4 {
            let index = plane.flat_index([i, j]).unwrap();
            assert_eq!(index, i + 4 * j);
            let geometry = plane.face(index).unwrap();
            let expected_area =
                patch.area_m2() * (EDGES_0[i + 1] - EDGES_0[i]) * (EDGES_1[j + 1] - EDGES_1[j])
                    / 4.0;
            close(geometry.area_m2, expected_area, 2e-15);
            for d in 0..3 {
                let expected_center = patch.center_m()[d]
                    + 0.5 * (EDGES_0[i + 1] + EDGES_0[i]) * patch.half_edges_m()[0][d]
                    + 0.5 * (EDGES_1[j + 1] + EDGES_1[j]) * patch.half_edges_m()[1][d];
                assert_eq!(geometry.center_m[d], expected_center);
            }
            let averaged = integrate(geometry, |x| frame.sample(x).unwrap().u);
            let sample = plane.samples()[index];
            for d in 0..3 {
                close(sample.favre_velocity_m_s[d], averaged[d], 3e-13);
                close(sample.mean_mass_flux_kg_m2_s[d], 1.19 * averaged[d], 3e-13);
                center_bias = center_bias
                    .max((averaged[d] - frame.sample(geometry.center_m).unwrap().u[d]).abs());
            }
            close(
                sample.normal_mass_flux_kg_m2_s,
                1.19 * dot(averaged, patch.normal()),
                3e-13,
            );
            area += geometry.area_m2;
        }
    }
    close(area, patch.area_m2(), 2e-15);
    assert!(
        center_bias > 1e-4,
        "fixture must distinguish face averages from center samples"
    );
}

#[test]
fn anelastic_faces_use_integrated_density_and_favre_velocity_over_large_contrasts() {
    let source = source(43);
    let frame = frame(&source);
    let constant = plane(&frame, MannInflowDensity::ConstantKgM3(1.19));
    let zero = ExponentialDensityProfile::new(1.19, 4.0, 0.0, [0.0, 8.0]).unwrap();
    let zero_plane = plane(&frame, MannInflowDensity::Exponential(zero));
    assert_eq!(constant, zero_plane);
    for kappa in [1e-10, 0.7, 20.0] {
        let density = ExponentialDensityProfile::new(1.19, 4.0, kappa, [0.0, 8.0]).unwrap();
        let plane = plane(&frame, MannInflowDensity::Exponential(density));
        let mut point_density_error: f64 = 0.0;
        for (index, sample) in plane.samples().iter().enumerate() {
            let geometry = plane.face(index).unwrap();
            let integral = integrate(geometry, |x| {
                let physical = frame.sample_anelastic(x, &density).unwrap();
                [
                    physical.density_kg_m3,
                    physical.mass_flux_kg_m2_s[0],
                    physical.mass_flux_kg_m2_s[1],
                    physical.mass_flux_kg_m2_s[2],
                ]
            });
            close(sample.mean_density_kg_m3, integral[0], 3e-12);
            for d in 0..3 {
                close(sample.mean_mass_flux_kg_m2_s[d], integral[d + 1], 3e-12);
                close(
                    sample.favre_velocity_m_s[d],
                    integral[d + 1] / integral[0],
                    3e-12,
                );
                close(
                    sample.mean_density_kg_m3 * sample.favre_velocity_m_s[d],
                    sample.mean_mass_flux_kg_m2_s[d],
                    3e-15,
                );
            }
            point_density_error = point_density_error.max(
                (sample.mean_density_kg_m3 / density.density_kg_m3(geometry.center_m[2]).unwrap()
                    - 1.0)
                    .abs(),
            );
        }
        if kappa == 20.0 {
            assert!(point_density_error > 1.0);
        }
        close(plane.turbulent_mass_flow_kg_s(), 0.0, 2e-13);
    }
}

#[test]
fn whole_patch_flux_survives_mesh_refinement_without_a_face_residual_correction() {
    for seed in [43, 87, 193] {
        let source = source(seed);
        let frame = frame(&source);
        let whole = MannInflowPlane::prepare(
            &frame,
            [&[-1.0, 1.0], &[-1.0, 1.0]],
            MannInflowDensity::ConstantKgM3(1.23),
            usize::MAX,
            || false,
        )
        .unwrap();
        for subdivisions in [1usize, 4, 13] {
            let a: Vec<_> = (0..=subdivisions)
                .map(|i| -1.0 + 2.0 * (i as f64 / subdivisions as f64).powi(2))
                .collect();
            let b: Vec<_> = (0..=subdivisions + 1)
                .map(|i| -1.0 + 2.0 * (i as f64 / (subdivisions + 1) as f64).powi(3))
                .collect();
            let refined = MannInflowPlane::prepare(
                &frame,
                [&a, &b],
                MannInflowDensity::ConstantKgM3(1.23),
                usize::MAX,
                || false,
            )
            .unwrap();
            let mut integral = [0.0; 3];
            let mut local_flux = 0.0;
            for (index, sample) in refined.samples().iter().enumerate() {
                let geometry = refined.face(index).unwrap();
                let filter = MannBoxFilter::prepare(
                    source.model(),
                    [geometry.half_edges_m[0], geometry.half_edges_m[1], [0.0; 3]],
                    16 * 8,
                    || false,
                )
                .unwrap();
                let original = frame.sample_filtered(geometry.center_m, &filter).unwrap();
                for d in 0..3 {
                    assert_eq!(
                        sample.favre_velocity_m_s[d].to_bits(),
                        original.u[d].to_bits()
                    );
                    integral[d] += geometry.area_m2 * sample.mean_mass_flux_kg_m2_s[d];
                }
                local_flux += geometry.area_m2 * sample.normal_mass_flux_kg_m2_s.abs();
            }
            for d in 0..3 {
                close(
                    integral[d],
                    patch().area_m2() * whole.samples()[0].mean_mass_flux_kg_m2_s[d],
                    3e-13,
                );
            }
            close(refined.turbulent_mass_flow_kg_s(), 0.0, 3e-13);
            assert!(local_flux > 1e-3, "local turbulent exchange is retained");
        }
    }
}

#[test]
fn both_clocks_refresh_atomically_reuse_buffers_and_replay_without_source_mutation() {
    let mut source = source(43);
    let mut frame = frame(&source);
    let mut plane = plane(&frame, MannInflowDensity::ConstantKgM3(1.19));
    let before = plane.clone();
    let addresses = [plane.samples.as_ptr(), plane.trial.as_ptr()];
    let faces_address = plane.faces.as_ptr();
    let bytes = plane.workspace_bytes();
    source.advance_to(source.tick() + 17).unwrap();
    source.refresh_inflow_frame(&mut frame, || false).unwrap();
    let source_checkpoint = source.checkpoint();
    let frame_before = frame.clone();
    for stop in [1, 7, 14] {
        let mut count = 0;
        let refusal = plane
            .refresh(&frame, || {
                count += 1;
                count == stop
            })
            .unwrap_err();
        assert_eq!(refusal.code, "mann-inlet-plane-cancelled");
        assert_eq!(plane, before);
        assert_eq!(frame, frame_before);
    }
    plane.refresh(&frame, || false).unwrap();
    assert_ne!(plane.samples(), before.samples());
    assert_eq!(plane.time_s(), frame.time_s());
    assert_eq!(source.checkpoint(), source_checkpoint);
    assert!(addresses.contains(&plane.samples.as_ptr()));
    assert_eq!(plane.faces.as_ptr(), faces_address);
    assert_eq!(plane.workspace_bytes(), bytes);

    let mut fixed = FixedStepMannField::from_snapshot(source, 239, 0.017, 88 * 16).unwrap();
    fixed.advance_to(5, 4 * 16 * 5, || false).unwrap();
    let checkpoint = fixed.checkpoint();
    fixed.refresh_inflow_frame(&mut frame, || false).unwrap();
    plane.refresh(&frame, || false).unwrap();
    assert_eq!(plane.time_s(), fixed.time_s());
    assert_eq!(fixed.checkpoint(), checkpoint);
    assert!(addresses.contains(&plane.samples.as_ptr()));
    assert_eq!(plane.faces.as_ptr(), faces_address);
    assert_eq!(plane.workspace_bytes(), bytes);
    let fresh = MannInflowPlane::prepare(
        &frame,
        [&EDGES_0, &EDGES_1],
        MannInflowDensity::ConstantKgM3(1.19),
        bytes,
        || false,
    )
    .unwrap();
    assert_eq!(plane, fresh);

    // A late numerical refusal must also preserve all previously accepted
    // results. Inject an unrepresentable last observation after geometry
    // admission, so earlier trial faces have already been evaluated.
    plane.faces.last_mut().unwrap().geometry.center_m[0] = f64::INFINITY;
    let accepted = plane.clone();
    let mut polls = 0;
    assert!(
        plane
            .refresh(&frame, || {
                polls += 1;
                false
            })
            .is_err()
    );
    assert!(polls >= 12);
    assert_eq!(plane, accepted);
    assert_eq!(frame.time_s(), fixed.time_s());
}

#[test]
fn partitions_source_density_and_actual_capacity_are_admitted_before_publication() {
    let source = source(43);
    let frame = frame(&source);
    let bytes = MannInflowPlane::required_workspace_bytes(16, [4, 3]).unwrap();
    let admitted = plane(&frame, MannInflowDensity::ConstantKgM3(1.19));
    assert_eq!(admitted.workspace_bytes(), bytes);
    assert!(
        MannInflowPlane::prepare(
            &frame,
            [&EDGES_0, &EDGES_1],
            MannInflowDensity::ConstantKgM3(1.19),
            bytes - 1,
            || false
        )
        .is_err()
    );
    assert!(MannInflowPlane::required_workspace_bytes(16, [usize::MAX, 2]).is_err());
    assert!(MannInflowPlane::required_workspace_bytes(16, [1, 0]).is_err());
    assert!(MannInflowPlane::required_workspace_bytes(257, [1, 1]).is_err());
    for edges in [
        &[][..],
        &[-1.0][..],
        &[-1.0, 0.0, 0.0, 1.0][..],
        &[-1.0, 0.2, -0.2, 1.0][..],
        &[-0.9, 1.0][..],
        &[-1.0, 0.9][..],
        &[-1.0, f64::NAN, 1.0][..],
    ] {
        assert!(
            MannInflowPlane::prepare(
                &frame,
                [edges, &EDGES_1],
                MannInflowDensity::ConstantKgM3(1.19),
                usize::MAX,
                || false
            )
            .is_err()
        );
    }
    for density in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            MannInflowPlane::prepare(
                &frame,
                [&EDGES_0, &EDGES_1],
                MannInflowDensity::ConstantKgM3(density),
                bytes,
                || false
            )
            .is_err()
        );
    }
    let excluded_corner = ExponentialDensityProfile::new(1.19, 4.0, 0.1, [3.5, 4.5]).unwrap();
    assert!(
        MannInflowPlane::prepare(
            &frame,
            [&EDGES_0, &EDGES_1],
            MannInflowDensity::Exponential(excluded_corner),
            bytes,
            || false
        )
        .is_err()
    );
    assert!(
        MannInflowPlane::prepare(
            &frame,
            [&EDGES_0, &EDGES_1],
            MannInflowDensity::ConstantKgM3(1.19),
            bytes,
            || true
        )
        .is_err()
    );
    let other_patch = MannInflowPatch::new([0.3, -0.7, 4.1], patch().half_edges_m()).unwrap();
    let wrong_patch = source
        .inflow_frame(other_patch, 168 * 16, || false)
        .unwrap();
    let wrong_source = self::frame(&self::source(44));
    for wrong in [&wrong_patch, &wrong_source] {
        let mut candidate = admitted.clone();
        assert!(candidate.refresh(wrong, || false).is_err());
        assert_eq!(candidate, admitted);
    }
}

#[test]
fn cartesian_reservoir_targets_add_mean_flow_without_losing_turbulent_mass_balance() {
    let source = source(43);
    // Normal +x: the lowest remaining axis y is the fast index, followed by z.
    let patch = MannInflowPatch::new([0.0, 0.0, 4.0], [[0.0, 0.8, 0.0], [0.0, 0.0, 0.6]]).unwrap();
    let frame = source.inflow_frame(patch, 168 * 16, || false).unwrap();
    let density = ExponentialDensityProfile::new(1.19, 4.0, 0.8, [0.0, 8.0]).unwrap();
    let plane = MannInflowPlane::prepare(
        &frame,
        [&EDGES_0, &EDGES_1],
        MannInflowDensity::Exponential(density),
        usize::MAX,
        || false,
    )
    .unwrap();
    let mean_velocity = [6.0, 0.2, 0.0];
    let mut target_mass_flow = 0.0;
    for z in 0..3 {
        for y in 0..4 {
            let index = plane.flat_index([y, z]).unwrap();
            let sample = plane.samples()[index];
            // These values map directly to NsfReservoirTileState / NsfOpenState;
            // temperature and the installed boundary's fixed cell axis are
            // caller-owned, so fs-atmo needs no dependency on that solver.
            let cell = [0, y, z];
            let velocity: [f64; 3] =
                core::array::from_fn(|d| mean_velocity[d] + sample.favre_velocity_m_s[d]);
            let temperature_k = 288.15;
            assert_eq!(index, cell[1] + 4 * cell[2]);
            assert!(temperature_k > 0.0);
            target_mass_flow +=
                plane.face(index).unwrap().area_m2 * sample.mean_density_kg_m3 * velocity[0];
        }
    }
    let whole = MannInflowFace {
        center_m: patch.center_m(),
        half_edges_m: patch.half_edges_m(),
        area_m2: patch.area_m2(),
    };
    let rho_bar = integrate(whole, |x| [density.density_kg_m3(x[2]).unwrap()])[0];
    close(
        target_mass_flow,
        patch.area_m2() * rho_bar * mean_velocity[0],
        4e-13,
    );
    // This is the prescribed reservoir target, not a promise about a
    // compressible solver's actual streamed mass at any individual step.
    close(plane.turbulent_mass_flow_kg_s(), 0.0, 2e-13);
}
