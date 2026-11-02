//! Analytic link traction, exact old-path parity and transactional refusal.

use super::*;
use crate::d3q19::W3;
use std::cell::{Cell, RefCell};

fn isolated(coordinates: [usize; 3], force: [f64; 3], sigma: f64) -> FreeSurface3 {
    FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        force,
        sigma,
        ContactModel3::Neutral,
        |x, y, z| [x, y, z] == coordinates,
    )
    .unwrap()
    .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal)
}

#[test]
fn distinct_neighbors_reconstruct_each_population_with_the_correct_traction_sign() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut flow = isolated([3; 3], [0.0; 3], 0.0);
            let mass = flow.ledger_mass();
            let seen = RefCell::new(Vec::new());
            flow.step_with_gas_link_pressure(&gate, 0.02, |link| {
                assert_eq!(link.liquid_cell, [3; 3]);
                let (ex, ey, ez) = E3[link.incoming_direction];
                assert_eq!(
                    link.gas_cell,
                    [3 - i64::from(ex), 3 - i64::from(ey), 3 - i64::from(ez)]
                );
                seen.borrow_mut().push(link.incoming_direction);
                // Identify pressure by the gas neighbor, independently of the
                // population index: higher pressure on the negative-axis side.
                Ok(sign * 0.001 * (3 - link.gas_cell[axis]) as f64)
            })
            .unwrap();
            assert_eq!(*seen.borrow(), (1..Q3).collect::<Vec<_>>());
            let (slot, lane) = flow.grid.resolve_source(3, 3, 3).unwrap();
            let populations = flow.grid.populations(slot, lane);
            for q in 0..Q3 {
                let direction = [E3[q].0, E3[q].1, E3[q].2];
                // At rest rho_b=1+3p and f_q=2*w_q*rho_b-w_q.
                // No production equilibrium/collision helper supplies this oracle.
                let expected = W3[q] * (1.0 + 6.0 * sign * 0.001 * f64::from(direction[axis]));
                assert!(
                    (populations[q] - expected).abs() < 2e-15,
                    "axis={axis} q={q}"
                );
            }
            let state = flow.sample_at(3, 3, 3).unwrap();
            assert!((state.density - 1.0).abs() < 2e-14);
            // Sum_i w_i c_i c_i = cs² I gives delta momentum = 2*p e_axis.
            for d in 0..3 {
                let expected = if d == axis { sign * 0.002 } else { 0.0 };
                assert!((state.velocity[d] - expected).abs() < 2e-14);
            }
            assert_eq!(flow.ledger_mass().to_bits(), mass.to_bits());
        }
    }
}

#[test]
fn uniform_link_load_keeps_cell_load_bits_and_original_callback_cadence() {
    let gate = CancelGate::new_clock_free();
    for model in [
        ConversionRedistribution3::LegacyGlobal,
        ConversionRedistribution3::Local,
        ConversionRedistribution3::ConservativeLocal,
    ] {
        let mut cell =
            isolated([3; 3], [1e-6, -2e-6, 3e-6], 0.001).with_conversion_redistribution(model);
        let mut link = cell.clone();
        for _ in 0..3 {
            let cell_calls = Cell::new(0);
            let link_calls = Cell::new(0);
            cell.step_with_gas_pressure(&gate, 0.05, |_| {
                cell_calls.set(cell_calls.get() + 1);
                Ok(0.0005)
            })
            .unwrap();
            link.step_with_gas_link_pressure(&gate, 0.05, |_| {
                link_calls.set(link_calls.get() + 1);
                Ok(0.0005)
            })
            .unwrap();
            assert_eq!(cell_calls.get(), 1);
            assert_eq!(link_calls.get(), Q3 - 1);
            assert_eq!(cell.grid.state_bits(), link.grid.state_bits());
            assert_eq!(cell.mass_snapshot(), link.mass_snapshot());
            assert_eq!(cell.conversions(), link.conversions());
        }
    }
}

#[test]
fn wall_sources_are_not_gas_links_and_coordinates_do_not_underflow() {
    let gate = CancelGate::new_clock_free();
    let mut flow = isolated([0; 3], [0.0; 3], 0.0);
    let calls = Cell::new(0);
    flow.step_with_gas_link_pressure(&gate, 0.05, |link| {
        assert_eq!(link.liquid_cell, [0; 3]);
        assert!(
            link.gas_cell
                .iter()
                .all(|&coordinate| (0..=1).contains(&coordinate))
        );
        assert!(link.gas_cell.iter().any(|&coordinate| coordinate == 1));
        calls.set(calls.get() + 1);
        Ok(0.0)
    })
    .unwrap();
    // Three positive-axis sources plus three positive face diagonals.
    assert_eq!(calls.get(), 6);
}

#[test]
fn last_link_pressure_refusal_and_cancellation_preserve_the_complete_candidate() {
    let mut flow = isolated([3; 3], [1e-6, 0.0, 0.0], 0.001);
    let original = flow.clone();
    for cancel in [false, true] {
        let gate = CancelGate::new_clock_free();
        let calls = Cell::new(0);
        let result = flow.step_with_gas_link_pressure(&gate, 0.05, |_| {
            calls.set(calls.get() + 1);
            if calls.get() == Q3 - 1 {
                if cancel {
                    gate.request();
                } else {
                    return Ok(0.02);
                }
            }
            Ok(0.0005)
        });
        assert_eq!(calls.get(), Q3 - 1);
        if cancel {
            assert_eq!(result, Err(FreeSurfaceError3::Cancelled));
        } else {
            assert!(matches!(
                result,
                Err(FreeSurfaceError3::InvalidGasPressure {
                    coordinates: Some([3, 3, 3]),
                    ..
                })
            ));
        }
        assert_eq!(flow.grid.state_bits(), original.grid.state_bits());
        assert_eq!(flow.mass_snapshot(), original.mass_snapshot());
        assert_eq!(flow.carry.to_bits(), original.carry.to_bits());
        assert_eq!(flow.conversions(), original.conversions());
        assert_eq!(
            flow.worst_step_violation().to_bits(),
            original.worst_step_violation().to_bits()
        );
        assert_eq!(flow.geometry_epoch, original.geometry_epoch);
    }
}
