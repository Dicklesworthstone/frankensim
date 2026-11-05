//! G1/G3 quadratic motion, independent triangle moments and coherent work.

use super::*;
use crate::acoustics::helmholtz::heterogeneous::{BoundaryAdmittance, HarmonicLimits};
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use crate::fixtures::masked_cube_grid;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};
const LENGTHS: [f64; 3] = [1.0, 0.4, 0.2];
const FACE_AREA: f64 = 0.5 * LENGTHS[1] * LENGTHS[2];
const MIDSIDE: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];

fn geometry() -> (fs_rep_mesh::TetComplex, Vec<[f64; 3]>) {
    let (mesh, mut positions) = masked_cube_grid(1, 1, 1, &|_, _, _| true);
    for point in &mut positions {
        for (coordinate, length) in point.iter_mut().zip(LENGTHS) {
            *coordinate *= length;
        }
    }
    (mesh, positions)
}

fn model<'f, 'c>(fe: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        fe,
        fe.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn source_face(fe: &HelmholtzFe<'_>) -> [u32; 3] {
    fe.boundary_faces()
        .into_iter()
        .find(|face| face.iter().all(|&v| fe.positions()[v as usize][0] == 0.0))
        .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn close_complex(actual: C64, expected: C64, tolerance: f64) {
    close(actual.re, expected.re, tolerance);
    close(actual.im, expected.im, tolerance);
}

// Integral of a barycentric monomial over one physical triangle:
// 2 A a! b! c! / (a+b+c+2)!. Independent of Duffy points or FE integration.
fn triangle_moment(powers: [usize; 3]) -> f64 {
    const FACTORIAL: [f64; 8] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0];
    2.0 * FACE_AREA * powers.iter().map(|&p| FACTORIAL[p]).product::<f64>()
        / FACTORIAL[powers.iter().sum::<usize>() + 2]
}

fn midside_load(fe: &HelmholtzFe<'_>, face: [u32; 3]) -> Vec<f64> {
    let mut expected = vec![0.0; fe.ndof()];
    for (local, vertex) in face.into_iter().enumerate() {
        let mut powers = [1, 1, 0];
        powers[local] += 1;
        expected[vertex as usize] = 4.0 * triangle_moment(powers);
    }
    for (a, b) in [(0, 1), (1, 2), (2, 0)] {
        let mut edge = [face[a], face[b]];
        edge.sort_unstable();
        let global = fe.space.complex.edges.binary_search(&edge).unwrap();
        let mut powers = [1, 1, 0];
        powers[a] += 1;
        powers[b] += 1;
        expected[fe.space.edge_off + global] = 4.0 * triangle_moment(powers);
    }
    expected
}

#[test]
fn quadratic_boundary_midside_motion_retains_original_basis_and_exact_work() {
    let (mesh, positions) = geometry();
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let model = model(&fe);
    let face = source_face(&fe);
    let source = model
        .quadratic_velocity_source(&[(face, MIDSIDE)], &mut || false)
        .unwrap();
    let expected = midside_load(&fe, face);
    close(
        source.integrated_shape_m2(),
        FACE_AREA / 3.0,
        2e-14 * FACE_AREA,
    );
    for (&actual, &truth) in source.load.iter().zip(&expected) {
        close(actual, truth, 2e-14 * FACE_AREA);
    }
    // A corner-only interface would erase this motion completely.
    assert!(
        fe.boundary_linear_load(&[(face, [0.0; 3])])
            .unwrap()
            .iter()
            .all(|&value| value == 0.0)
    );
    assert!(source.load.iter().any(|&value| value > 0.0));
    let pressure: Vec<_> = (0..model.ndof())
        .map(|i| C64::new(0.2 * i as f64 - 1.3, 0.4 - 0.1 * i as f64))
        .collect();
    let exact = pressure
        .iter()
        .zip(&expected)
        .fold(C64::ZERO, |sum, (&p, &weight)| sum + p.scale(weight));
    close_complex(source.harmonic_pressure(&pressure), exact, 1e-14);
    let initial: Vec<_> = pressure.iter().map(|p| p.re).collect();
    let next: Vec<_> = pressure.iter().map(|p| p.im).collect();
    let midpoint: f64 = initial
        .iter()
        .zip(&next)
        .zip(expected)
        .map(|((&a, &b), weight)| 0.5 * (a + b) * weight)
        .sum();
    close(source.midpoint_pressure(&initial, &next), midpoint, 1e-14);
}

#[test]
fn quadratic_boundary_preserves_higher_order_face_tests() {
    let (mesh, positions) = geometry();
    for order in [3, 4] {
        let fe = HelmholtzFe::new(&mesh, &positions, order, MEDIUM).unwrap();
        let model = model(&fe);
        let face = source_face(&fe);
        let source = model
            .quadratic_velocity_source(&[(face, MIDSIDE)], &mut || false)
            .unwrap();
        let index = mesh.faces.binary_search(&face).unwrap();
        // The first face basis is lambda0 lambda1 lambda2 at either order.
        let dof = fe.space.face_off + index * fe.space.per_face;
        close(
            source.load[dof],
            4.0 * triangle_moment([2, 2, 1]),
            2e-14 * FACE_AREA,
        );
        assert!(source.load[dof] > 0.0);
    }
}

#[test]
fn quadratic_boundary_permutation_linear_recovery_and_zero_integral_shape() {
    let (mesh, positions) = geometry();
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let model = model(&fe);
    let face = source_face(&fe);
    let values = [0.4, -0.8, 1.2, 0.7, -0.6, 0.3];
    let source = model
        .quadratic_velocity_source(&[(face, values)], &mut || false)
        .unwrap();
    for permutation in [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [0, 2, 1],
        [2, 1, 0],
        [1, 0, 2],
    ] {
        let mut reordered = [0.0; 6];
        for i in 0..3 {
            reordered[i] = values[permutation[i]];
        }
        for (edge, (a, b)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
            let mut pair = [permutation[a], permutation[b]];
            pair.sort_unstable();
            reordered[3 + edge] = values[match pair {
                [0, 1] => 3,
                [1, 2] => 4,
                [0, 2] => 5,
                _ => unreachable!("a triangle permutation preserves its edges"),
            }];
        }
        let permuted = model
            .quadratic_velocity_source(&[(permutation.map(|i| face[i]), reordered)], &mut || false)
            .unwrap();
        assert_eq!(source.faces, permuted.faces);
        close(
            source.integrated_shape_m2(),
            permuted.integrated_shape_m2(),
            2e-14 * FACE_AREA,
        );
        for (&a, &b) in source.load.iter().zip(&permuted.load) {
            close(a, b, 2e-14 * FACE_AREA);
        }
    }
    let corners = [-0.5, 1.2, 0.8];
    let linear = model
        .quadratic_velocity_source(
            &[(
                face,
                [
                    corners[0],
                    corners[1],
                    corners[2],
                    f64::midpoint(corners[0], corners[1]),
                    f64::midpoint(corners[1], corners[2]),
                    f64::midpoint(corners[2], corners[0]),
                ],
            )],
            &mut || false,
        )
        .unwrap();
    for (&a, b) in linear
        .load
        .iter()
        .zip(fe.boundary_linear_load(&[(face, corners)]).unwrap())
    {
        close(a, b, 2e-14 * FACE_AREA);
    }
    // N0=2 lambda0^2-lambda0 has zero mean, but nonzero pressure coupling.
    let balanced = model
        .quadratic_velocity_source(&[(face, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0])], &mut || false)
        .unwrap();
    close(balanced.integrated_shape_m2(), 0.0, 2e-14 * FACE_AREA);
    assert!(balanced.load.iter().any(|v| v.abs() > FACE_AREA * 1e-3));
}

#[test]
fn quadratic_boundary_harmonic_midside_motion_and_signed_coherent_power() {
    let (mesh, positions) = geometry();
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let model = model(&fe);
    let face = source_face(&fe);
    let source = model
        .quadratic_velocity_source(&[(face, MIDSIDE)], &mut || false)
        .unwrap();
    let owner = model
        .prepare(
            0.31,
            &[BoundaryAdmittance {
                faces: fe.boundary_faces(),
                admittance_m_pa_s: C64::from_re(
                    1.0 / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s),
                ),
            }],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let amplitude = C64::new(0.03, -0.02);
    let field = owner
        .solve_with_boundary_sources(&[], &[(&source, amplitude)], &[], &mut || false)
        .unwrap();
    let original_work = field
        .pressure_pa
        .iter()
        .zip(midside_load(&fe, face))
        .fold(C64::ZERO, |sum, (&p, weight)| sum + p.scale(weight));
    let expected = 0.5 * (original_work.conj() * amplitude).re;
    assert!(expected > 0.0);
    let tolerance = expected * 2e-10;
    close(field.boundary_source_power_w, expected, tolerance);
    close(field.input_power_w, expected, tolerance);
    close(field.boundary_dissipation_w, expected, tolerance);
    assert_eq!(field.volume_source_power_w, 0.0);
    assert!(field.pressure_pa.iter().any(|p| p.abs() > 0.0));
    let cancelled = owner
        .solve_with_boundary_sources(
            &[],
            &[(&source, amplitude), (&source, -amplitude)],
            &[],
            &mut || false,
        )
        .unwrap();
    assert!(cancelled.pressure_pa.iter().all(|&p| p == C64::ZERO));
    assert_eq!(cancelled.boundary_source_power_w, 0.0);

    let uniform = model
        .quadratic_velocity_source(&[(face, [1.0; 6])], &mut || false)
        .unwrap();
    let absorbing = owner
        .solve_with_boundary_sources(
            &[(face, [amplitude; 3])],
            &[(&uniform, amplitude.scale(-0.25))],
            &[],
            &mut || false,
        )
        .unwrap();
    assert!(absorbing.input_power_w > 0.0);
    assert!(absorbing.boundary_source_power_w < 0.0);
    close(
        absorbing.boundary_source_power_w,
        -absorbing.input_power_w / 3.0,
        absorbing.input_power_w * 2e-10,
    );
    close(
        absorbing.power_balance_residual_w(),
        0.0,
        absorbing.input_power_w * 2e-10,
    );

    let volume = model
        .volume_velocity_source(
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / (LENGTHS[0] * LENGTHS[1] * LENGTHS[2])),
            &mut || false,
        )
        .unwrap();
    let volume_amplitude = C64::new(0.0002, 0.0001);
    let mixed = owner
        .solve_with_boundary_sources(
            &[],
            &[(&source, amplitude)],
            &[(&volume, volume_amplitude)],
            &mut || false,
        )
        .unwrap();
    let boundary_work = 0.5 * (source.harmonic_pressure(&mixed.pressure_pa).conj() * amplitude).re;
    let volume_work =
        0.5 * (volume.harmonic_pressure(&mixed.pressure_pa).conj() * volume_amplitude).re;
    let tolerance = (boundary_work.abs() + volume_work.abs()) * 2e-10;
    close(mixed.boundary_source_power_w, boundary_work, tolerance);
    close(mixed.volume_source_power_w, volume_work, tolerance);
    close(mixed.input_power_w, boundary_work + volume_work, tolerance);
    close(mixed.power_balance_residual_w(), 0.0, tolerance);

    let legacy = owner
        .solve(&[(face, [amplitude; 3])], &mut || false)
        .unwrap();
    let empty = owner
        .solve_with_boundary_sources(&[(face, [amplitude; 3])], &[], &[], &mut || false)
        .unwrap();
    assert_eq!(legacy.pressure_pa, empty.pressure_pa);
    assert_eq!(legacy.input_power_w, empty.input_power_w);
    assert_eq!(legacy.boundary_source_power_w, 0.0);
    assert_eq!(empty.boundary_source_power_w, 0.0);
}

#[test]
fn quadratic_boundary_support_identity_release_limits_and_cancellation_refuse() {
    let (mesh, positions) = geometry();
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let model = model(&fe);
    let face = source_face(&fe);
    assert!(model.quadratic_velocity_source(&[], &mut || false).is_err());
    assert!(
        model
            .quadratic_velocity_source(
                &[(face, MIDSIDE), ([face[1], face[0], face[2]], MIDSIDE)],
                &mut || false,
            )
            .is_err()
    );
    let exterior = fe.boundary_faces();
    let internal = *mesh
        .faces
        .iter()
        .find(|face| !exterior.contains(face))
        .unwrap();
    assert!(
        model
            .quadratic_velocity_source(&[(internal, MIDSIDE)], &mut || false)
            .is_err()
    );
    for invalid_values in [
        [f64::NAN; 6],
        [f64::INFINITY; 6],
        [0.0, 0.0, 0.0, f64::MAX, f64::MAX, f64::MAX],
    ] {
        assert!(
            model
                .quadratic_velocity_source(&[(face, invalid_values)], &mut || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        model.quadratic_velocity_source(&[(face, MIDSIDE)], &mut || {
            polls += 1;
            polls >= 4
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    let source = model
        .quadratic_velocity_source(&[(face, MIDSIDE)], &mut || false)
        .unwrap();
    let owner = model
        .prepare(0.31, &[], HarmonicLimits::default(), &mut || false)
        .unwrap();
    let other_model = HeterogeneousHelmholtzFe::new(
        &fe,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let wrong = other_model
        .quadratic_velocity_source(&[(face, MIDSIDE)], &mut || false)
        .unwrap();
    assert!(
        owner
            .solve_with_boundary_sources(&[], &[(&wrong, C64::ZERO)], &[], &mut || false)
            .is_err()
    );
    assert!(
        owner
            .solve_with_boundary_sources(
                &[],
                &[(&source, C64::new(f64::NAN, 0.0))],
                &[],
                &mut || false,
            )
            .is_err()
    );
    let too_many = vec![(&source, C64::ZERO); MAX_BOUNDARY_SOURCES + 1];
    assert!(matches!(
        owner.solve_with_boundary_sources(&[], &too_many, &[], &mut || false),
        Err(HeterogeneousError::Limit {
            maximum: MAX_BOUNDARY_SOURCES,
            ..
        })
    ));
    assert!(matches!(
        owner.solve_with_boundary_sources(&[], &[(&source, C64::ONE)], &[], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let released = model
        .prepare_with_pressure_release(0.31, &[], &[face], HarmonicLimits::default(), &mut || false)
        .unwrap();
    let zero_shape = model
        .quadratic_velocity_source(&[(face, [0.0; 6])], &mut || false)
        .unwrap();
    for blocked in [&source, &zero_shape] {
        assert!(
            released
                .solve_with_boundary_sources(&[], &[(blocked, C64::ZERO)], &[], &mut || false)
                .is_err()
        );
    }
}
