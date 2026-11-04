use super::*;
use fs_la::nonlinear_eig::PoleOptions;

fn options(n: usize) -> BeynOptions {
    BeynOptions {
        nodes: 128,
        probe_columns: n,
        rank_relative_tolerance: 1e-9,
        count_tolerance: 1e-7,
        quadrature_tolerance: 1e-9,
        contour_min_reciprocal_condition: 1e-12,
        linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 100,
        eigen_max_sweeps: 200,
        max_dimension: 8,
        max_dense_bytes: 1 << 20,
        max_work: 100_000_000,
        refinement: PoleOptions {
            spectral_scale: 4.0,
            step_tolerance: 1e-12,
            residual_tolerance: 1e-12,
            max_iterations: 24,
            max_dimension: 8,
            max_dense_bytes: 1 << 20,
            max_work: 1_000_000,
        },
    }
}

fn check_rows(p: &Pencil, result: &BeynReport) {
    assert_eq!(result.argument_count, result.rank);
    assert_eq!(result.rank, result.poles.len());
    for pole in &result.poles {
        assert_eq!(pole.vector.len(), p.dimension);
        assert!(
            p.relative_residual(pole.value, &pole.vector, &mut || false)
                .unwrap()
                < 5e-12
        );
    }
}

#[test]
fn two_loaded_quadratics_are_discovered_without_seeds_or_row_scale_bias() {
    // G1/G3: independent physical oscillators, 2*q''+4*q'+10*q=0 and
    // 2*q''+8*q'+26*q=0. The region selects one member of each conjugate pair.
    // Original equation units differ by 160 decades; no coordinate is removed.
    let mut p = Pencil::new(2, 4.0);
    p.constant = vec![10e80, 0.0, 0.0, 26e-80];
    p.linear = vec![4e80, 0.0, 0.0, 8e-80];
    p.quadratic_diagonal = vec![2e80, 2e-80];
    p.condition(&mut || false).unwrap();
    let region = Ellipse {
        center: C64::new(-1.5, -2.5),
        real_radius: 1.6,
        imaginary_radius: 1.6,
    };
    let result = p.discover(region, &options(2), || false).unwrap();
    check_rows(&p, &result);
    assert_eq!(result.argument_count, 2);
    for (pole, target) in result
        .poles
        .iter()
        .zip([C64::new(-2.0, -3.0), C64::new(-1.0, -2.0)])
    {
        assert!((pole.value - target).abs() < 2e-11);
    }
    let replay = p.discover(region, &options(2), || false).unwrap();
    assert_eq!(result.argument_integral, replay.argument_integral);
    for (first, second) in result.poles.iter().zip(replay.poles) {
        assert_eq!(first.value, second.value);
        assert_eq!(first.vector, second.vector);
    }
}

#[test]
fn contour_retains_coupled_maxwell_and_dark_relaxations() {
    // G1: q''+6*q'+6*q+(q-z)=0, z'=2*(q-z), w'=-5*w.
    // The physical characteristic polynomial is (s+1)(s+3)(s+4)(s+5).
    // Select -3,-4,-5; the last root has identically zero displacement.
    let mut p = Pencil::new(3, 4.0);
    p.constant = vec![7.0, -1.0, 0.0, -2.0, 2.0, 0.0, 0.0, 0.0, 5.0];
    p.linear = vec![6.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    p.quadratic_diagonal[0] = 1.0;
    p.condition(&mut || false).unwrap();
    let region = Ellipse {
        center: C64::from_re(-4.0),
        real_radius: 1.6,
        imaginary_radius: 0.9,
    };
    let result = p.discover(region, &options(3), || false).unwrap();
    check_rows(&p, &result);
    assert_eq!(result.argument_count, 3);
    for (pole, target) in result.poles.iter().zip([-5.0, -4.0, -3.0]) {
        assert!((pole.value - C64::from_re(target)).abs() < 2e-11);
    }
    let dark = &result.poles[0];
    assert!(dark.vector[0].abs() < 1e-11);
    assert!(dark.vector[1].abs() < 1e-11);
    assert!((dark.vector[2].abs() - 1.0).abs() < 1e-12);
    let mut narrow = options(3);
    narrow.probe_columns = 2;
    assert!(matches!(
        p.discover(region, &narrow, || false),
        Err(BeynError::ProbeRankSaturated {
            probes: 2,
            count: 3
        })
    ));
}

#[test]
fn contour_keeps_finite_gap_inertia_and_algebraic_pressure() {
    // G1: q''+q+dp=0, v'+2*v=dp, v=q'. Retain the gap's kinetic
    // coordinate and pressure multiplier; no inverse of the singular mass.
    // The finite roots are those of 2*s^2+2*s+1, with the actual added mass.
    let mut p = Pencil::new(3, 1.0);
    p.constant = vec![1.0, 0.0, 1.0, 0.0, 2.0, -1.0, 0.0, 1.0, 0.0];
    p.linear = vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0];
    p.quadratic_diagonal[0] = 1.0;
    p.condition(&mut || false).unwrap();
    let region = Ellipse {
        center: C64::from_re(-0.5),
        real_radius: 0.9,
        imaginary_radius: 1.0,
    };
    let result = p.discover(region, &options(3), || false).unwrap();
    check_rows(&p, &result);
    assert_eq!(result.argument_count, 2);
    for pole in &result.poles {
        let s = pole.value;
        assert!((s.re + 0.5).abs() < 2e-11);
        assert!((s.im.abs() - 0.5).abs() < 2e-11);
        assert!((pole.vector[1] - s * pole.vector[0]).abs() < 2e-11);
        assert!((pole.vector[2] - (s + C64::from_re(2.0)) * pole.vector[1]).abs() < 2e-11);
    }
}

#[test]
fn contour_refuses_unresolved_quadratic_count_and_bounds_without_partial_modes() {
    // G0/G4: one displacement has two finite poles. Two resolvent moments
    // cannot resolve both in one state direction; success must not lose one.
    let mut p = Pencil::new(1, 2.0);
    p.constant[0] = 5.0;
    p.linear[0] = 2.0;
    p.quadratic_diagonal[0] = 1.0;
    p.condition(&mut || false).unwrap();
    let region = Ellipse {
        center: C64::from_re(-1.0),
        real_radius: 3.0,
        imaginary_radius: 3.0,
    };
    assert!(matches!(
        p.discover(region, &options(1), || false),
        Err(BeynError::CountUnresolved { .. })
    ));
    let empty = Ellipse {
        center: C64::from_re(8.0),
        real_radius: 1.0,
        imaginary_radius: 1.0,
    };
    let result = p.discover(empty, &options(1), || false).unwrap();
    assert!(result.poles.is_empty());
    assert_eq!(result.argument_count, 0);
    assert_eq!(result.rank, 0);
    let mut bytes = options(1);
    bytes.max_dense_bytes = 1;
    assert!(matches!(
        p.discover(empty, &bytes, || false),
        Err(BeynError::Invalid { .. })
    ));
    let mut work = options(1);
    work.max_work = 1;
    assert!(matches!(
        p.discover(empty, &work, || false),
        Err(BeynError::Invalid { .. })
    ));
    assert!(matches!(
        p.discover(empty, &options(1), || true),
        Err(BeynError::Operator(NonlinearEigError::Cancelled))
    ));
    let mut polls = 0;
    assert!(matches!(
        p.discover(empty, &options(1), || {
            polls += 1;
            polls > 20
        }),
        Err(BeynError::Operator(NonlinearEigError::Cancelled))
    ));
}
