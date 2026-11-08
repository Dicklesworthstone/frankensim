use fs_solid::{
    SolidError,
    stability::{buckling_loads, buckling_loads_cancellable},
};
use fs_sparse::{Coo, Csr};
fn matrix(n: usize, f: impl Fn(usize, usize) -> f64) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            let v = f(i, j);
            if v != 0.0 {
                coo.push(i, j, v);
            }
        }
    }
    coo.assemble()
}
fn diagonal(v: &[f64]) -> Csr {
    matrix(v.len(), |i, j| if i == j { v[i] } else { 0.0 })
}
fn map(n: usize) -> Vec<Option<usize>> {
    (0..n + (n % 2))
        .map(|i| if i < n { Some(i) } else { None })
        .collect()
}
fn close(a: f64, b: f64) {
    assert!(
        a.is_finite() && (a - b).abs() < 2e-7 * b.abs().max(1.0),
        "{a:e} != {b:e}"
    );
}
#[test]
fn scaling_preserves_loads_and_geometric_work() {
    for unit in [1e-240, 1.0, 1e240] {
        for reference in [1e-24, 1.0, 1e24] {
            let k = diagonal(&[2., 6., 12., 20., 30., 42.].map(|v| v * unit));
            let g = diagonal(&[-1., -2., -3., -4., -5., -6.].map(|v| v * unit * reference));
            let r = buckling_loads(&k, &g, &map(6), 2, 150).unwrap();
            assert_eq!(r.loads.len(), 2);
            for i in 0..2 {
                close(r.loads[i] * reference, 2.0 + i as f64);
                let mut gp = vec![0.0; 6];
                g.spmv(&r.modes[i], &mut gp);
                close(
                    -r.modes[i].iter().zip(gp).map(|(a, b)| a * b).sum::<f64>(),
                    1.0,
                );
                assert!(r.relative_residuals[i] <= 1e-8);
            }
        }
    }
}
#[test]
fn tiny_full_blocks_and_noncompressive_modes_are_admitted() {
    let r = buckling_loads(&diagonal(&[2.]), &diagonal(&[-1.]), &map(1), 1, 150).unwrap();
    assert_eq!(r.loads, vec![2.]);
    let k = diagonal(&[2., 5., 7., 11.]);
    let r = buckling_loads(&k, &diagonal(&[-1., -1., 0., 1.]), &map(4), 4, 150).unwrap();
    assert_eq!(r.loads.len(), 2);
    close(r.loads[0], 2.);
    close(r.loads[1], 5.);
    let r = buckling_loads(&k, &diagonal(&[0.; 4]), &map(4), 4, 150).unwrap();
    assert!(r.loads.is_empty() && r.modes.is_empty());
}
#[test]
fn underresolved_iteration_refuses_then_converges() {
    let k = Csr::identity(20);
    let g = diagonal(&(1..=20).map(|i| -(i as f64)).collect::<Vec<_>>());
    assert!(matches!(
        buckling_loads(&k, &g, &map(20), 1, 1),
        Err(SolidError::SolveFailed { .. })
    ));
    let r = buckling_loads(&k, &g, &map(20), 1, 200).unwrap();
    close(r.loads[0], 1.0 / 20.0);
}
#[test]
fn noncommuting_congruence_recovers_reference_pencil() {
    let l = [
        [2., 0., 0., 0.],
        [1., 3., 0., 0.],
        [-0.5, 0.25, 1.5, 0.],
        [0.2, -0.3, 0.7, 1.],
    ];
    let d = [4., 2., 0.5, -1.];
    let k = matrix(4, |i, j| (0..4).map(|p| l[i][p] * l[j][p]).sum());
    let g = matrix(4, |i, j| {
        -(0..4).map(|p| l[i][p] * d[p] * l[j][p]).sum::<f64>()
    });
    let r = buckling_loads(&k, &g, &map(4), 4, 150).unwrap();
    assert_eq!(r.loads.len(), 3);
    for (a, b) in r.loads.iter().zip([0.25, 0.5, 2.]) {
        close(*a, b);
    }
}
#[test]
fn malformed_pencils_are_structured_refusals() {
    let k = Csr::identity(4);
    let g = diagonal(&[-1.; 4]);
    let bad_maps = [
        vec![Some(0), Some(0), Some(2), Some(3)],
        vec![Some(0), Some(1)],
        vec![Some(0), Some(1), Some(2), Some(4)],
    ];
    for m in bad_maps {
        assert!(matches!(
            buckling_loads(&k, &g, &m, 1, 10),
            Err(SolidError::InvalidInput { .. })
        ));
    }
    for bad in [
        diagonal(&[1., 2., f64::NAN, 4.]),
        matrix(4, |i, j| {
            if i == j {
                1.
            } else if i == 0 && j == 1 {
                0.1
            } else {
                0.
            }
        }),
    ] {
        assert!(matches!(
            buckling_loads(&bad, &g, &map(4), 1, 10),
            Err(SolidError::InvalidInput { .. })
        ));
    }
}

#[test]
fn component_independent_maps_preserve_odd_three_dimensional_layouts() {
    let map = [None, Some(0), None];
    let result = buckling_loads(&diagonal(&[2.0]), &diagonal(&[-1.0]), &map, 1, 10)
        .expect("one xyz node with one free component is a valid pencil");
    assert_eq!(result.loads, [2.0]);
    assert_eq!(result.dof_map, map);
}

#[test]
fn cancellation_refuses_before_and_during_numerical_work() {
    use std::cell::Cell;
    let k = Csr::identity(20);
    let g = diagonal(&(1..=20).map(|i| -(i as f64)).collect::<Vec<_>>());
    for stop_at in [1, 9, 35] {
        let polls = Cell::new(0);
        let result = buckling_loads_cancellable(&k, &g, &map(20), 1, 200, &|| {
            polls.set(polls.get() + 1);
            polls.get() >= stop_at
        });
        assert!(
            matches!(result, Err(SolidError::Cancelled)),
            "poll {stop_at}"
        );
        assert_eq!(polls.get(), stop_at);
    }
    let result = buckling_loads_cancellable(&k, &g, &map(20), 1, 200, &|| false).unwrap();
    close(result.loads[0], 0.05);
    assert!(
        result.iters < 200,
        "a resolved block stops spending its allowance"
    );
}
