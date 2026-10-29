//! G0/G1: a load sees its physical flow map, independently of external drives.
use fs_phs::linear_exact::{LinearSplitBudget, split_quadratic, step_render};
use fs_phs::{PortHamiltonian, QuadraticStorage, modal_bank_ports};

fn host() -> PortHamiltonian {
    // Two original unit inertias; external force acts only on the first.
    PortHamiltonian::new(
        4,
        1,
        vec![
            0., 1., 0., 0., -1., 0., 0., 0., 0., 0., 0., 1., 0., 0., -1., 0.,
        ],
        vec![0.; 16],
        vec![0., 1., 0., 0.],
        Box::new(
            QuadraticStorage::new(
                vec![
                    0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
                ],
                4,
            )
            .unwrap(),
        ),
    )
    .unwrap()
}
fn load() -> PortHamiltonian {
    modal_bank_ports(&[3.0], &[0.2], &[&[2.0]]).unwrap()
}

#[test]
fn signed_flow_projection_has_reciprocal_force_and_preserves_the_external_drive() {
    let s = host()
        .with_parallel_load_columns(load(), &[0., 0.3, 0., -0.8])
        .unwrap();
    assert_eq!(s.state_dim(), 6);
    assert_eq!(s.port_dim(), 1);
    let (j, _, g) = s.structure();
    assert_eq!(g, &[0., 1., 0., 0., 0., 0.]);
    assert_eq!(j[1 * 6 + 5], -0.6);
    assert_eq!(j[3 * 6 + 5], 1.6);
    assert_eq!(j[5 * 6 + 1], 0.6);
    assert_eq!(j[5 * 6 + 3], -1.6);
    let x = [0., 0.7, 0., -0.5, 0., 0.3];
    let effort = s.effort(&x);
    let power: f64 = (0..6)
        .map(|i| effort[i] * (0..6).map(|k| j[i * 6 + k] * effort[k]).sum::<f64>())
        .sum();
    assert!(power.abs() < 1e-14);
    assert_eq!(s.output(&x), vec![0.7]);
}

#[test]
fn initially_unforced_coordinate_moves_and_complete_work_closes() {
    let s = host()
        .with_parallel_load_columns(load(), &[0., 0.3, 0., -0.8])
        .unwrap();
    let split = split_quadratic(
        &s,
        LinearSplitBudget {
            max_modes: 6,
            max_omega_dt: 3.0,
        },
    )
    .unwrap();
    let mut x = vec![0.; 6];
    let (mut heat, mut work) = (0., 0.);
    for i in 0..100 {
        let f = if i < 50 { 1.0 } else { 0.0 };
        let record = step_render(&s, &split, &x, &[f], 0.01).unwrap();
        assert!((record.supplied - f * (record.x[0] - x[0])).abs() < 1e-12);
        heat += record.dissipated;
        work += record.supplied;
        x = record.x;
    }
    assert!(
        x[3].abs() > 1e-3,
        "the undriven velocity must feel mutual reaction"
    );
    assert!(heat > 0.0);
    assert!((s.hamiltonian(&x) + heat - work).abs() < 1e-11);
    // The old index API is exactly the corresponding column projection.
    let old = host().with_parallel_load(load(), &[0]).unwrap();
    let new = host()
        .with_parallel_load_columns(load(), &[0., 1., 0., 0.])
        .unwrap();
    assert_eq!(old.structure(), new.structure());
}

#[test]
fn malformed_columns_and_nonfinite_projections_refuse() {
    for columns in [
        vec![],
        vec![0.; 3],
        vec![0.; 5],
        vec![0., f64::NAN, 0., 1.],
        vec![0., f64::INFINITY, 0., 1.],
    ] {
        assert!(host().with_parallel_load_columns(load(), &columns).is_err());
    }
    let uncoupled = host().with_parallel_load_columns(load(), &[0.; 4]).unwrap();
    assert_eq!(uncoupled.port_dim(), 1);
}
