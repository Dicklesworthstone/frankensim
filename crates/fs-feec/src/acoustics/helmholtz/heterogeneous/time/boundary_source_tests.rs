use super::tests::{MEDIUM, brick, plane, uniform};
use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::MAX_BOUNDARY_SOURCES;

pub(super) fn quadratic_profile(
    geometry: &HelmholtzFe<'_>,
    faces: &[[u32; 3]],
    shape: impl Fn([f64; 3]) -> f64,
) -> Vec<([u32; 3], [f64; 6])> {
    faces
        .iter()
        .map(|&face| {
            let points = face.map(|vertex| geometry.positions()[vertex as usize]);
            let mut values = [0.0; 6];
            for i in 0..3 {
                values[i] = shape(points[i]);
            }
            for (i, [a, b]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
                values[3 + i] = shape(core::array::from_fn(|d| {
                    f64::midpoint(points[a][d], points[b][d])
                }));
            }
            (face, values)
        })
        .collect()
}

// Independent exact nodal P2 triangle mass, area/180 times this integer matrix.
// It integrates the product of original FE pressure and the supplied six-node
// velocity polynomial without reusing source quadrature or its prepared load.
pub(super) fn quadratic_work(
    geometry: &HelmholtzFe<'_>,
    profiles: &[([u32; 3], [f64; 6])],
    initial: &[f64],
    next: &[f64],
    amplitude: f64,
    duration: f64,
) -> f64 {
    assert_eq!(geometry.order(), 2);
    const MASS: [[f64; 6]; 6] = [
        [6., -1., -1., 0., -4., 0.],
        [-1., 6., -1., 0., 0., -4.],
        [-1., -1., 6., -4., 0., 0.],
        [0., 0., -4., 32., 16., 16.],
        [-4., 0., 0., 16., 32., 16.],
        [0., -4., 0., 16., 16., 32.],
    ];
    let mut result = 0.0;
    for &(face, values) in profiles {
        let element = geometry
            .tetrahedra()
            .iter()
            .position(|tet| face.iter().all(|v| tet.contains(v)))
            .unwrap();
        let tet = geometry.tetrahedra()[element];
        let local = face.map(|v| tet.iter().position(|&i| i == v).unwrap());
        let mut points = [[0.0; 4]; 6];
        for i in 0..3 {
            points[i][local[i]] = 1.0;
        }
        for (i, [a, b]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
            points[3 + i][local[a]] = 0.5;
            points[3 + i][local[b]] = 0.5;
        }
        let pressure = points.map(|barycentric| {
            geometry
                .space
                .element_dofs(element)
                .iter()
                .map(|&(dof, basis)| {
                    f64::midpoint(initial[dof], next[dof]) * basis.eval(barycentric, 2)
                })
                .sum::<f64>()
        });
        let vertices = face.map(|v| geometry.positions()[v as usize]);
        let a: [f64; 3] = core::array::from_fn(|d| vertices[1][d] - vertices[0][d]);
        let b: [f64; 3] = core::array::from_fn(|d| vertices[2][d] - vertices[0][d]);
        let cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let area = 0.5 * det::sqrt(cross.iter().map(|v| v * v).sum::<f64>());
        for i in 0..6 {
            for j in 0..6 {
                result +=
                    duration * amplitude * area / 180.0 * pressure[i] * MASS[i][j] * values[j];
            }
        }
    }
    result
}

#[test]
fn midside_boundary_motion_drives_resolved_transient_field_with_original_pressure_work() {
    let lengths = [1.0, 0.2, 0.3];
    let area = lengths[1] * lengths[2];
    let (complex, positions) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let inlet = plane(&geometry, 0.0);
    let profiles = quadratic_profile(&geometry, &inlet, |x| {
        let y = x[1] / lengths[1];
        4.0 * y * (1.0 - y)
    });
    assert!(profiles.iter().all(|(_, values)| values[..3] == [0.0; 3]));
    assert!(
        profiles
            .iter()
            .any(|(_, values)| values[3..].iter().any(|&v| v != 0.0))
    );
    let source = model
        .quadratic_velocity_source(&profiles, &mut || false)
        .unwrap();
    assert!((source.integrated_shape_m2() / (2.0 * area / 3.0) - 1.0).abs() < 2e-12);
    let owner = model
        .prepare_time_with_pressure_release(
            0.005,
            &[],
            &plane(&geometry, 1.0),
            HeterogeneousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let zero = owner.zero_state();
    let corner_only = inlet
        .iter()
        .map(|&face| (face, [0.0; 3]))
        .collect::<Vec<_>>();
    let silent = owner.step(&zero, &corner_only, &mut || false).unwrap();
    assert!(silent.state.pressure_pa().iter().all(|&p| p == 0.0));
    let mut state = zero;
    let mut input = 0.0;
    for _ in 0..6 {
        let next = owner
            .step_with_boundary_sources(&state, &[], &[(&source, 0.02)], &[], &mut || false)
            .unwrap();
        let expected = quadratic_work(
            &geometry,
            &profiles,
            state.pressure_pa(),
            next.state.pressure_pa(),
            0.02,
            0.005,
        );
        assert!(expected > 0.0);
        assert!((next.work.boundary_source_work_j / expected - 1.0).abs() < 2e-10);
        assert!((next.work.input_work_j / expected - 1.0).abs() < 2e-10);
        assert_eq!(next.work.volume_source_work_j, 0.0);
        input += expected;
        assert!((owner.energy(&next.state).unwrap().total_j() / input - 1.0).abs() < 3e-10);
        for (dof, &fixed) in owner.pressure_release.constrained.iter().enumerate() {
            if fixed {
                assert_eq!(next.state.pressure_pa()[dof], 0.0);
                assert_eq!(next.state.pressure_potential_pa_s()[dof], 0.0);
            }
        }
        state = next.state;
    }
    let reversed = owner
        .step_with_boundary_sources(&state, &[], &[(&source, -0.02)], &[], &mut || false)
        .unwrap();
    assert!(reversed.work.boundary_source_work_j < 0.0);
}

#[test]
fn prepared_boundary_sources_recover_linear_forcing_and_refuse_wrong_model_or_release_face() {
    let (complex, positions) = brick(1, [1.0, 0.2, 0.3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let foreign = uniform(&geometry);
    let inlet = plane(&geometry, 0.0);
    let profiles = quadratic_profile(&geometry, &inlet, |x| 1.0 + x[1]);
    let source = model
        .quadratic_velocity_source(&profiles, &mut || false)
        .unwrap();
    let foreign_source = foreign
        .quadratic_velocity_source(&profiles, &mut || false)
        .unwrap();
    let release_source = model
        .quadratic_velocity_source(
            &quadratic_profile(&geometry, &plane(&geometry, 1.0), |_| 0.0),
            &mut || false,
        )
        .unwrap();
    let volume_source = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let owner = model
        .prepare_time_with_pressure_release(
            0.005,
            &[],
            &plane(&geometry, 1.0),
            HeterogeneousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let zero = owner.zero_state();
    let amplitude = 0.02;
    let linear = profiles
        .iter()
        .map(|&(face, values)| {
            (
                face,
                [
                    amplitude * values[0],
                    amplitude * values[1],
                    amplitude * values[2],
                ],
            )
        })
        .collect::<Vec<_>>();
    let original = owner
        .step_with_volume_sources(&zero, &linear, &[(&volume_source, -0.01)], &mut || false)
        .unwrap();
    let prepared = owner
        .step_with_boundary_sources(
            &zero,
            &[],
            &[(&source, amplitude)],
            &[(&volume_source, -0.01)],
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in original
        .state
        .pressure_pa()
        .iter()
        .zip(prepared.state.pressure_pa())
    {
        assert!((a - b).abs() < 2e-10);
    }
    let scale = prepared.work.input_work_j.abs();
    assert!((original.work.input_work_j - prepared.work.input_work_j).abs() < 2e-10 * scale);
    assert!(
        (prepared.work.boundary_source_work_j + prepared.work.volume_source_work_j
            - prepared.work.input_work_j)
            .abs()
            < 2e-10 * scale
    );
    for invalid in [
        &[(&foreign_source, 0.0)][..],
        &[(&source, f64::NAN)][..],
        &[(&release_source, 0.0)][..],
    ] {
        assert!(
            owner
                .step_with_boundary_sources(&zero, &[], invalid, &[], &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        owner.step_with_boundary_sources(
            &zero,
            &[],
            &vec![(&source, 0.0); MAX_BOUNDARY_SOURCES + 1],
            &[],
            &mut || false
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    assert!(matches!(
        owner.step_with_boundary_sources(&zero, &[], &[(&source, 0.02)], &[], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let cancelled = owner
        .step_with_boundary_sources(
            &zero,
            &[],
            &[(&source, 0.02), (&source, -0.02)],
            &[],
            &mut || false,
        )
        .unwrap();
    assert_eq!(cancelled.work.input_work_j, 0.0);
    assert_eq!(cancelled.work.boundary_source_work_j, 0.0);
    assert!(cancelled.state.pressure_pa().iter().all(|&p| p == 0.0));
    let unchanged = owner
        .step_with_boundary_sources(&zero, &[], &[], &[], &mut || false)
        .unwrap();
    assert_eq!(unchanged.state.pressure_pa(), zero.pressure_pa());
    assert_eq!(zero.pressure_pa(), vec![0.0; model.ndof()]);
}
