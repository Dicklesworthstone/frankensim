use super::super::super::tests::{plane, state};
use super::*;

fn sphere(centre: f64) -> SpherePanels {
    let original = SpherePanels::icosphere(1.0, 0).unwrap();
    SpherePanels::from_triangles(
        original
            .triangles()
            .unwrap()
            .iter()
            .map(|triangle| triangle.map(|p| [p[0], p[1], p[2] + centre]))
            .collect(),
    )
    .unwrap()
}
fn close(actual: &VectorObservation, expected: &VectorObservation) {
    assert_eq!(actual.scalar.pressure.len(), expected.scalar.pressure.len());
    for (a, b) in actual
        .scalar
        .pressure
        .iter()
        .zip(&expected.scalar.pressure)
        .chain(
            actual
                .particle_velocity_m_s
                .iter()
                .flatten()
                .zip(expected.particle_velocity_m_s.iter().flatten()),
        )
    {
        assert!((*a - *b).abs() < 3e-11 * (1.0 + b.abs()), "{a:?} != {b:?}");
    }
    for (a, b) in actual
        .scalar
        .quadrature_error_estimate_pa
        .iter()
        .zip(&expected.scalar.quadrature_error_estimate_pa)
        .chain(
            actual
                .quadrature_error_estimate_m_s
                .iter()
                .flatten()
                .zip(expected.quadrature_error_estimate_m_s.iter().flatten()),
        )
    {
        assert!((*a - *b).abs() < 3e-11 * (1.0 + b.abs()));
    }
}

#[test]
fn one_body_adaptive_pair_recovers_single_owner_inside_and_outside() {
    let gate = CancelGate::new();
    let surface = sphere(0.0);
    let outside = state(293.15, 101325.0, Some(0.5));
    let inside = state(383.15, 101325.0, Some(5.0));
    let multiple = MolecularFluidBodies::new(
        &[Body {
            surface: &surface,
            material: &inside,
        }],
        &outside,
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let single = MolecularFluidInterface::new(
        &surface,
        &outside,
        &inside,
        200.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let actual = multiple.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let expected = single.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let options = MolecularNearFieldOptions::default();
    for (point, interior) in [([0.2, 0.1, 0.3], true), ([0.3, 0.2, 1.1], false)] {
        let a = if interior {
            actual.interior_pressure_velocity_adaptive_at(0, &[point], 0.01, 1000, options, &gate)
        } else {
            actual.total_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options, &gate)
        }
        .unwrap();
        let b = if interior {
            expected.interior_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options, &gate)
        } else {
            expected.total_pressure_velocity_adaptive_at(&[point], 0.01, 1000, options, &gate)
        }
        .unwrap();
        close(&a, &b);
    }
}

#[test]
fn adaptive_transparent_interior_receives_neighbors_with_one_global_work_budget() {
    let gate = CancelGate::new();
    let surfaces = [sphere(-1.6), sphere(1.6)];
    let outside = state(293.15, 101325.0, Some(0.5));
    let inside = state(383.15, 101325.0, Some(5.0));
    let owner = MolecularFluidBodies::new(
        &[
            Body {
                surface: &surfaces[0],
                material: &outside,
            },
            Body {
                surface: &surfaces[1],
                material: &inside,
            },
        ],
        &outside,
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let single = MolecularFluidInterface::new(
        &surfaces[1],
        &outside,
        &inside,
        200.0,
        Options::default(),
        &gate,
    )
    .unwrap();
    let actual = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let expected = single.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    let points = [[0.12, -0.08, -1.3]];
    let options = MolecularNearFieldOptions::default();
    let a = actual
        .interior_pressure_velocity_adaptive_at(0, &points, 0.01, 1000, options, &gate)
        .unwrap();
    let b = expected
        .total_pressure_velocity_adaptive_at(&points, 0.01, 1000, options, &gate)
        .unwrap();
    close(&a, &b);
    let incident = plane()
        .sample_wave(
            owner.wave,
            outside.gas().density * 200.0,
            None,
            points[0],
            &mut || false,
        )
        .unwrap();
    assert!((a.scalar.pressure[0] - incident.pressure_pa).abs() > 1e-5);
    assert!(
        actual
            .interior_pressure_velocity_adaptive_at(0, &points, 0.01, 39, options, &gate)
            .is_err()
    );
    for refused in [
        MolecularNearFieldOptions {
            maximum_row_bytes: crate::near_field::complex::required_row_bytes(20, 1).unwrap(),
            ..options
        },
        MolecularNearFieldOptions {
            quadrature: crate::near_field::Options {
                maximum_kernel_evaluations: 39 * 80,
                ..options.quadrature
            },
            ..options
        },
    ] {
        assert!(
            actual
                .interior_pressure_velocity_adaptive_at(0, &points, 0.01, 1000, refused, &gate)
                .is_err()
        );
    }
    assert!(
        actual
            .interior_pressure_velocity_adaptive_at(
                0,
                &[[0.0, 0.0, 1.5]],
                0.01,
                1000,
                options,
                &gate
            )
            .is_err()
    );
    let replay = actual
        .interior_pressure_velocity_adaptive_at(0, &points, 0.01, 1000, options, &gate)
        .unwrap();
    assert_eq!(a.scalar.pressure, replay.scalar.pressure);
    assert_eq!(a.particle_velocity_m_s, replay.particle_velocity_m_s);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        actual.interior_pressure_velocity_adaptive_at(0, &points, 0.01, 1000, options, &cancelled),
        Err(HelmholtzError::Cancelled)
    ));

    let transparent = MolecularFluidBodies::new(
        &[
            Body {
                surface: &surfaces[0],
                material: &outside,
            },
            Body {
                surface: &surfaces[1],
                material: &outside,
            },
        ],
        &outside,
        200.0,
        0.01,
        Options::default(),
        &gate,
    )
    .unwrap();
    let field = transparent
        .scatter_field(plane(), 0.01, 1000, &gate)
        .unwrap();
    let exact = field
        .interior_pressure_velocity_adaptive_at(0, &points, 0.01, 1000, options, &gate)
        .unwrap();
    assert_eq!(exact.scalar.pressure, [incident.pressure_pa]);
    assert_eq!(exact.particle_velocity_m_s, [incident.velocity_m_s]);
    assert_eq!(exact.scalar.quadrature_error_estimate_pa, [0.0]);
    assert_eq!(exact.quadrature_error_estimate_m_s, [[0.0; 3]]);
}
