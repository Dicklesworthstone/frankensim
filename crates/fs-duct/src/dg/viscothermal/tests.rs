use super::*;
use crate::{LossModel, Termination, input_impedance};
use fs_la::factor::lu;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;

fn air(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn cylinder() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.05,
        }],
    }
}
fn duct(elements: usize, inlet: DgBoundary, outlet: DgBoundary) -> DgDuct {
    DgDuct::new(
        &cylinder(),
        &air(101325.0),
        DgMesh { elements, order: 2 },
        inlet,
        outlet,
    )
    .unwrap()
}

// Harmonic response of the ACTUAL assembled operator. Condense only its
// independent memory coordinates, then use fs-la's existing real LU on the
// doubled complex system. This tests the physical coupling and mass maps,
// independently of the pole producer's analytic response helper.
fn harmonic_port(duct: &DgDuct, omega: f64) -> C64 {
    let n = duct.total_dofs();
    let base = 2 * duct.field_dofs;
    let size = 2 * base;
    let mut matrix = vec![0.0; size * size];
    for i in 0..base {
        for j in 0..base {
            let mut value = C64::new(-duct.l_tilde[i * n + j], if i == j { omega } else { 0.0 });
            for k in base..n {
                value = value
                    - C64::from_re(duct.l_tilde[i * n + k] * duct.l_tilde[k * n + j])
                        / C64::new(-duct.l_tilde[k * n + k], omega);
            }
            matrix[i * size + j] = value.re;
            matrix[i * size + base + j] = -value.im;
            matrix[(base + i) * size + j] = value.im;
            matrix[(base + i) * size + base + j] = value.re;
        }
    }
    let port = &duct.ports[0].1;
    let mut solution = vec![0.0; size];
    solution[..base].copy_from_slice(&port[..base]);
    lu(&matrix, size).unwrap().solve(&mut solution);
    (0..base).fold(C64::ZERO, |sum, i| {
        sum + C64::new(solution[i], solution[base + i]).scale(port[i])
    })
}

#[test]
fn narrow_dg_viscothermal_impedance_converges_to_exact_bessel() {
    let gas = air(101325.0);
    let mut previous = f64::INFINITY;
    for elements in [1, 2, 4] {
        let dg = duct(elements, DgBoundary::Port(PortId(0)), DgBoundary::IdealOpen)
            .with_viscothermal(16, [50.0, 8000.0], 1e-4)
            .unwrap();
        let mut worst = 0.0_f64;
        for frequency in [100.0, 500.0, 1000.0] {
            let omega = core::f64::consts::TAU * frequency;
            let exact = input_impedance(
                &cylinder(),
                &gas,
                omega,
                LossModel::Bessel,
                Termination::IdealOpen,
            )
            .unwrap()
            .impedance
            .conj();
            let actual = harmonic_port(&dg, omega);
            assert!(actual.re > 0.0);
            worst = worst.max((actual - exact).abs() / exact.abs());
        }
        println!("DG elements={elements} exact-Bessel complex error={worst:e}");
        assert!(worst < previous);
        previous = worst;
    }
    assert!(previous < 2e-3);
}

#[test]
fn dg_viscothermal_memory_uses_one_passive_exponential_and_bound_state() {
    let original = duct(1, DgBoundary::Closed, DgBoundary::Closed);
    let dg = original
        .clone()
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap();
    assert_ne!(dg.operator_identity(), original.operator_identity());
    assert_eq!(dg.total_dofs(), 2 * dg.field_dofs() * 9);
    let initial = dg.project_initial_state(|x| det::sin(core::f64::consts::PI * x / 0.05), |_| 0.0);
    assert_eq!(initial.tilde.len(), dg.total_dofs());
    assert!(
        initial.tilde[2 * dg.field_dofs()..]
            .iter()
            .all(|v| *v == 0.0)
    );
    let coarse = dg.propagator(2e-5, 1e-11).unwrap();
    let fine = dg.propagator(1e-5, 1e-11).unwrap();
    let mut state = initial.clone();
    let mut refined = initial.clone();
    let mut previous = state.stored_energy_j();
    for _ in 0..32 {
        let row = coarse.step(&mut state, &[]).unwrap();
        fine.step(&mut refined, &[]).unwrap();
        fine.step(&mut refined, &[]).unwrap();
        assert!(!row.ledger.balance_checked); // independent loss integral is not claimed
        assert!(row.ledger.cumulative_dissipated_j >= 0.0);
        assert!(state.stored_energy_j() <= previous * (1.0 + 1e-12));
        previous = state.stored_energy_j();
    }
    assert!(previous < 0.99 * initial.stored_energy_j());
    let difference = state
        .tilde
        .iter()
        .zip(&refined.tilde)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>();
    assert!(difference < 1e-16 * initial.stored_energy_j());
    coarse.sync_nodal(&mut state).unwrap();
    assert!(state.pressure().iter().all(|p| p.is_finite()));
    let mut wrong = original.zero_state();
    let before = wrong.clone();
    assert!(coarse.step(&mut wrong, &[]).is_err());
    assert_eq!(wrong, before);
}

#[test]
fn material_pressure_changes_bessel_drag_and_the_complete_duct_response() {
    let a = duct(2, DgBoundary::Port(PortId(0)), DgBoundary::IdealOpen)
        .with_viscothermal(16, [50.0, 8000.0], 1e-4)
        .unwrap();
    let b = DgDuct::new(
        &cylinder(),
        &air(60000.0),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_viscothermal(16, [50.0, 8000.0], 1e-4)
    .unwrap();
    assert!(
        b.viscothermal_model().unwrap().series_drag_per_s()
            > a.viscothermal_model().unwrap().series_drag_per_s()
    );
    let za = harmonic_port(&a, core::f64::consts::TAU * 500.0);
    let zb = harmonic_port(&b, core::f64::consts::TAU * 500.0);
    assert!((za - zb).abs() / za.abs() > 0.1);
}

#[test]
fn dg_wall_memory_refuses_unsupported_sources_and_oversized_states() {
    let valid = duct(1, DgBoundary::Closed, DgBoundary::IdealOpen);
    assert!(
        valid
            .clone()
            .with_viscothermal(2, [50.0, 8000.0], 1e-4)
            .is_err()
    );
    let attached = valid.with_viscothermal(16, [50.0, 8000.0], 1e-4).unwrap();
    assert!(
        attached
            .with_viscothermal(16, [50.0, 8000.0], 1e-4)
            .is_err()
    );
    let cone = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.001,
            outlet_radius: 0.002,
            length: 0.05,
        }],
    };
    let cone = DgDuct::new(
        &cone,
        &air(101325.0),
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    )
    .unwrap();
    assert!(cone.with_viscothermal(16, [50.0, 8000.0], 1e-4).is_err());
    let large = duct(16, DgBoundary::Closed, DgBoundary::IdealOpen);
    assert!(large.with_viscothermal(32, [50.0, 8000.0], 1e-4).is_err());
    let field = crate::profile::gas_field::tests::field([290.0, 310.0], [0.8; 2], 101325.0);
    let geometry = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: field.length_m(),
        }],
    };
    let nonuniform = DgDuct::new_with_gas_field(
        &geometry,
        &field,
        DgMesh {
            elements: 4,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    )
    .unwrap();
    assert!(
        nonuniform
            .with_viscothermal(16, [50.0, 8000.0], 1e-4)
            .is_err()
    );
}
