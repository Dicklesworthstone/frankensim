//! Independent source-energy and endpoint checks for the normal path operator.
use super::*;

#[test]
fn force_and_endpoint_tangent_have_correct_contact_and_equal_depth_limits() {
    let law = source(1e7, None);
    let d = 2e-5;
    let dt = 1e-3;
    let endpoint = law.evaluate("endpoint", d, 0.0, dt).unwrap();
    let zero = law.elastic_path_average(0.0, 0.0, dt).unwrap();
    assert_eq!(zero.mean_force_n, 0.0);
    assert_eq!(zero.endpoint_derivative_n_per_m, 0.0);
    assert_eq!(zero.work_j, 0.0);
    let fixed = law.elastic_path_average(d, d, dt).unwrap();
    close(fixed.mean_force_n, endpoint.normal_force_n);
    close(
        fixed.endpoint_derivative_n_per_m,
        0.5 * endpoint.tangent_n_per_m,
    );
    assert_eq!(fixed.work_j, 0.0);
    let enter = law.elastic_path_average(0.0, d, dt).unwrap();
    close(enter.mean_force_n, 0.4 * endpoint.normal_force_n);
    close(
        enter.endpoint_derivative_n_per_m,
        0.4 * endpoint.tangent_n_per_m,
    );
    close(enter.work_j, endpoint.reversible_energy_j);
    let leave = law.elastic_path_average(d, 0.0, dt).unwrap();
    close(leave.mean_force_n, enter.mean_force_n);
    close(
        leave.endpoint_derivative_n_per_m,
        (4.0 / 15.0) * endpoint.tangent_n_per_m,
    );
    close(leave.work_j, -endpoint.reversible_energy_j);
    // One-ulp travel must not disappear in a difference of endpoint potentials.
    let adjacent = f64::from_bits(d.to_bits() + 1);
    let tiny = law.elastic_path_average(d, adjacent, dt).unwrap();
    close(tiny.mean_force_n, endpoint.normal_force_n);
    close(tiny.work_j, endpoint.normal_force_n * (adjacent - d));
    assert!(tiny.work_j > 0.0);
}

#[test]
fn unequal_depth_force_tangent_and_subdivision_match_original_source_energy() {
    let law = source(1e7, None);
    let dt = 1e-3;
    let d0 = 1e-5;
    let d1 = 7e-5;
    let dm = 3e-5;
    let energy = |d| {
        law.evaluate("energy-reference", d, 0.0, dt)
            .unwrap()
            .reversible_energy_j
    };
    let path = law.elastic_path_average(d0, d1, dt).unwrap();
    close(path.work_j, energy(d1) - energy(d0));
    close(path.mean_force_n, (energy(d1) - energy(d0)) / (d1 - d0));
    // Differentiate the independently evaluated source-energy secant.
    let h = d1 * 1e-5;
    let plus = (energy(d1 + h) - energy(d0)) / (d1 + h - d0);
    let minus = (energy(d1 - h) - energy(d0)) / (d1 - h - d0);
    let fd = (plus - minus) / (2.0 * h);
    assert!((fd - path.endpoint_derivative_n_per_m).abs() < 2e-8 * fd.abs());
    let first = law.elastic_path_average(d0, dm, dt / 3.0).unwrap();
    let second = law.elastic_path_average(dm, d1, 2.0 * dt / 3.0).unwrap();
    close(first.work_j + second.work_j, path.work_j);
    close(
        (first.mean_force_n + 2.0 * second.mean_force_n) / 3.0,
        path.mean_force_n,
    );
    let reverse = law.elastic_path_average(d1, d0, dt).unwrap();
    close(reverse.mean_force_n, path.mean_force_n);
    close(reverse.work_j, -path.work_j);
}

#[test]
fn elliptic_source_geometry_and_original_admission_limits_are_preserved() {
    let a = solid("A", 1000.0, 1e7);
    let b = solid("B", 2000.0, 2e7);
    let law = SourcePointContact::new(
        bind(&a, &b, None, 0.0),
        PointContactGeometry::EllipticParaboloid {
            maximum_curvature_m_inverse: 20.0,
            minimum_curvature_m_inverse: 5.0,
        },
        envelope(),
    )
    .unwrap();
    let d0 = 1e-5;
    let d1 = 6e-5;
    let r0 = law.evaluate("ellipse-start", d0, 0.0, 1e-3).unwrap();
    let r1 = law.evaluate("ellipse-end", d1, 0.0, 1e-3).unwrap();
    let average = law.elastic_path_average(d0, d1, 1e-3).unwrap();
    assert!(r1.elliptic_patch_axes.unwrap().aspect_ratio > 1.0);
    close(
        average.work_j,
        r1.reversible_energy_j - r0.reversible_energy_j,
    );
    for (start, end, dt) in [
        (d0, 1.01e-4, 1e-3),
        (1.01e-4, d0, 1e-3),
        (0.0, d1, 1e-5),
        (d1, 0.0, 1e-5),
        (d0, d1, 2e-3),
        (d0, d1, 0.0),
        (-d0, d1, 1e-3),
        (d0, -d1, 1e-3),
        (d0, f64::NAN, 1e-3),
    ] {
        assert!(law.elastic_path_average(start, end, dt).is_err());
    }
    for chi in [0.0, 0.4] {
        assert!(
            source(1e7, Some(chi))
                .elastic_path_average(d0, d1, 1e-3)
                .is_err()
        );
    }
}
