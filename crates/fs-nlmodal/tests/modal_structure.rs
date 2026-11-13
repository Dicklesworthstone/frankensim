//! fs-nlmodal's `ModalStructure` refuses a modal model whose modes carry an observer transfer: the
//! Component books its modes' work and heat only, so a model that also radiates to an observer
//! would lose that energy from the ledger. The same model without its transfers is admitted.
//! Below that: the structure's port-map tangent against its own proposals, its transpose, and the
//! end of its flows, told on signal outputs when asked.
//! Run with `cargo test -p fs-nlmodal --test modal_structure`.

// Exact float comparisons here assert bitwise agreement or exact zeros.
#![allow(clippy::float_cmp)]

use fs_math::c64::C64;
use fs_nlmodal::component::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalStructure,
    ModalStructureRefusal,
};
use fs_ports::component::Component;

#[test]
fn a_model_with_an_observer_transfer_is_refused() {
    let budget = ModalAcousticTimeBudget::audible_reference();
    let mode = |hz: f64, transfer: C64| ModalAcousticMode {
        angular_frequency_rad_s: core::f64::consts::TAU * hz,
        damping_ratio: 1.0e-4,
        pressure_per_modal_velocity: transfer,
    };
    let model = ModalAcousticTimeModel::try_new(
        44_100,
        vec![mode(300.0, C64::ZERO), mode(900.0, C64::new(0.0, 1.0e-3))],
        budget,
    )
    .expect("model");
    assert_eq!(
        ModalStructure::new(model.clone(), Vec::new()).unwrap_err(),
        ModalStructureRefusal::ObserverTransfer { mode: 1 }
    );
    let silent =
        ModalStructure::new(model.without_pressure_transfers(), Vec::new()).expect("no transfer");
    assert_eq!((silent.modes(), silent.ports().len()), (2, 1));
}

// The end of a step and the port map's tangent (bead rc-1109 .12.5, commit 2).

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_nlmodal::component::{ModalAcousticState, ModalAttachment};
use fs_ports::component::{EndOfStep, PortLayout, PortTangent, Proposal};
use fs_ports::{PortKind, StableId};

const RATE_HZ: u32 = 48_000;
const KELVIN: f64 = 300.0;
/// Two rows of unequal weights on three modes: a bridge and a finger.
const BRIDGE: [f64; 3] = [0.9, -0.4, 0.3];
const FINGER: [f64; 3] = [0.5, 1.2, -0.8];
/// The efforts the tangent is taken at [N].
const EFFORTS: [f64; 2] = [0.3, -0.7];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn modal(omega: f64) -> ModalAcousticMode {
    ModalAcousticMode {
        angular_frequency_rad_s: omega,
        damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO,
    }
}

/// A damped mode (2 % of critical at 200 Hz), an undamped one (700 Hz) and a free coordinate
/// damped at 3 /s, on the reference budget.
fn three_modes() -> ModalAcousticTimeModel {
    let tau = core::f64::consts::TAU;
    ModalAcousticTimeModel::try_new_with_damping_rates(
        RATE_HZ,
        vec![modal(tau * 200.0), modal(tau * 700.0), modal(0.0)],
        &[2.0 * 0.02 * tau * 200.0, 0.0, 3.0],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .expect("three modes")
}

fn row(id: &str, weights: [f64; 3]) -> ModalAttachment {
    ModalAttachment {
        id: StableId::new(id).expect("id"),
        kind: PortKind::MechanicalForceVelocity,
        width: NonZeroUsize::MIN,
        weights: weights.to_vec(),
    }
}

fn two_rows(model: ModalAcousticTimeModel) -> ModalStructure {
    ModalStructure::new(model, vec![row("bridge", BRIDGE), row("finger", FINGER)])
        .expect("structure")
        .publishing_ends()
        .expect("ends")
}

fn moving() -> Vec<ModalAcousticState> {
    [(1.0e-6, 2.0e-3), (-3.0e-7, 1.0e-3), (2.0e-6, -4.0e-4)]
        .iter()
        .map(|&(q, v)| ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        })
        .collect()
}

fn durations() -> [f64; 3] {
    let full = 1.0 / f64::from(RATE_HZ);
    [full, 0.5 * full, 0.37 * full]
}

/// One proposal's outputs and candidate. A proposal the structure refuses or shortens is never
/// compared as a column: it fails the test here.
fn proposed(
    structure: &ModalStructure,
    old: &[ModalAcousticState],
    dt_s: f64,
    inputs: &[f64],
) -> (Vec<f64>, Vec<ModalAcousticState>) {
    let layout = PortLayout::new(structure.ports()).expect("layout");
    let mut outputs = vec![f64::NAN; layout.output_len()];
    let mut new = Vec::new();
    let answer =
        with_cx(|cx| structure.propose(&old.to_vec(), dt_s, inputs, &mut new, &mut outputs, cx));
    match answer {
        Ok(Proposal::Step(_)) => (outputs, new),
        other => panic!("a proposal at {inputs:?} was not taken: {other:?}"),
    }
}

fn jvp(
    structure: &ModalStructure,
    old: &[ModalAcousticState],
    dt_s: f64,
    inputs: &[f64],
    direction: &[f64],
) -> Vec<f64> {
    let layout = PortLayout::new(structure.ports()).expect("layout");
    let mut d_outputs = vec![f64::NAN; layout.output_len()];
    with_cx(|cx| structure.jvp(&old.to_vec(), dt_s, inputs, direction, &mut d_outputs, cx))
        .expect("jvp");
    d_outputs
}

fn unit(index: usize) -> Vec<f64> {
    let mut direction = vec![0.0; 3];
    direction[index] = 1.0;
    direction
}

#[test]
fn each_column_of_the_tangent_is_the_change_a_newton_makes_to_a_proposal() {
    let structure = two_rows(three_modes());
    let budget = ModalAcousticTimeBudget::audible_reference();
    let old = moving();
    let base = [EFFORTS[0], EFFORTS[1], KELVIN];
    let rows = [BRIDGE, FINGER];
    let mut heat_worst = 0.0_f64;
    for dt_s in durations() {
        let (at, reached) = proposed(&structure, &old, dt_s, &base);
        for column in 0..2 {
            let (mut up, mut down) = (base, base);
            up[column] += 1.0;
            down[column] -= 1.0;
            let ((above, high), (below, low)) = (
                proposed(&structure, &old, dt_s, &up),
                proposed(&structure, &old, dt_s, &down),
            );
            // The budgets these proposals keep clear of: every candidate's displacement under
            // 1e-5 of the 1.0 m sqrt(kg) allowed, its velocity under 1e-7 of 1e5 m sqrt(kg)/s,
            // and its energy under 1e-14 of 1e9 J.
            for candidate in [&reached, &high, &low] {
                for s in candidate {
                    assert!(
                        s.displacement_m_sqrt_kg.abs()
                            <= 1.0e-5 * budget.maximum_abs_displacement_m_sqrt_kg
                    );
                    assert!(
                        s.velocity_m_sqrt_kg_per_s.abs()
                            <= 1.0e-7 * budget.maximum_abs_velocity_m_sqrt_kg_per_s
                    );
                }
                let energy = structure.stored_energy_j(candidate).expect("energy");
                assert!(energy <= 1.0e-14 * budget.maximum_total_energy_j);
            }
            let tangent = jvp(&structure, &old, dt_s, &base, &unit(column));
            for r in 0..2 {
                // A mean flow is a sum of endpoint differences over dt, and rounds at the size of
                // the endpoints; an end flow rounds at the size of its terms.
                let mut mean_scale = 0.0_f64;
                let mut end_scale = 0.0_f64;
                for k in 0..3 {
                    for candidate in [&reached, &high] {
                        mean_scale += rows[r][k].abs()
                            * (candidate[k].displacement_m_sqrt_kg.abs()
                                + old[k].displacement_m_sqrt_kg.abs())
                            / dt_s;
                        end_scale += rows[r][k].abs() * candidate[k].velocity_m_sqrt_kg_per_s.abs();
                    }
                }
                let mean = above[r] - at[r];
                let end = above[3 + r] - at[3 + r];
                assert!(
                    (tangent[r] - mean).abs() <= 1.0e-12 * mean_scale,
                    "dt {dt_s:e}, mean row {r}, column {column}: {:e} against {mean:e}",
                    tangent[r]
                );
                assert!(
                    (tangent[3 + r] - end).abs() <= 1.0e-12 * end_scale,
                    "dt {dt_s:e}, end row {r}, column {column}: {:e} against {end:e}",
                    tangent[3 + r]
                );
            }
            // The heat row is quadratic in the efforts, so a central difference is its slope.
            let central = 0.5 * (above[2] - below[2]);
            let off = (tangent[2] - central).abs() / tangent[2].abs().max(central.abs());
            eprintln!(
                "heat column: dt {dt_s:e}, column {column}: tangent {:e}, central difference {central:e}, relative {off:e}",
                tangent[2]
            );
            heat_worst = heat_worst.max(off);
        }
    }
    eprintln!("heat column: worst relative difference {heat_worst:e}");
    // A regression guard of a measured value, not an accuracy claim. The step is exact for a held
    // force and a central difference is exact for a quadratic row, so the difference is roundoff.
    // The owner's frame bounds the roundoff of its heat, not of this row, which the structure
    // forms from terms the size of the work, so no bound is derived from the integrator.
    // Measured worst: 4.007e-11 (ts2, ca204c7b16, 2026-11-13; the 0.37-sample step, column 1).
    // The bound is twice that, rounded up to one digit.
    assert!(heat_worst <= 9.0e-11, "heat column off by {heat_worst:e}");
}

#[test]
fn the_transposed_tangent_is_the_transpose() {
    let structure = two_rows(three_modes());
    let old = moving();
    let base = [EFFORTS[0], EFFORTS[1], KELVIN];
    // A fixed sequence of directions and weights (a small linear congruential generator).
    let mut seed = 0x2545_f491_4f6c_dd1d_u64;
    let mut draw = || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
    };
    for dt_s in durations() {
        for _ in 0..8 {
            let u: Vec<f64> = (0..3).map(|_| draw()).collect();
            let w: Vec<f64> = (0..5).map(|_| draw()).collect();
            let ju = jvp(&structure, &old, dt_s, &base, &u);
            let mut jtw = vec![f64::NAN; 3];
            with_cx(|cx| structure.vjp(&old, dt_s, &base, &w, &mut jtw, cx)).expect("vjp");
            let left: f64 = w.iter().zip(&ju).map(|(a, b)| a * b).sum();
            let right: f64 = jtw.iter().zip(&u).map(|(a, b)| a * b).sum();
            let size: f64 = w.iter().zip(&ju).map(|(a, b)| (a * b).abs()).sum::<f64>()
                + jtw.iter().zip(&u).map(|(a, b)| (a * b).abs()).sum::<f64>();
            assert!(
                (left - right).abs() <= 4.0 * 5.0 * f64::EPSILON * size,
                "dt {dt_s:e}: <w, J u> = {left:e}, <J^T w, u> = {right:e}"
            );
            assert_eq!(jtw[2], 0.0, "the temperature has no adjoint");
        }
    }
}

#[test]
fn the_temperature_moves_nothing_and_an_undamped_structure_has_no_heat_slope() {
    let structure = two_rows(three_modes());
    let old = moving();
    let base = [EFFORTS[0], EFFORTS[1], KELVIN];
    for dt_s in durations() {
        assert!(
            jvp(&structure, &old, dt_s, &base, &unit(2))
                .iter()
                .all(|d| *d == 0.0)
        );
    }
    let tau = core::f64::consts::TAU;
    let undamped = ModalAcousticTimeModel::try_new(
        RATE_HZ,
        vec![modal(tau * 200.0), modal(tau * 700.0), modal(tau * 1300.0)],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .expect("undamped");
    let structure = two_rows(undamped);
    for dt_s in durations() {
        for column in 0..2 {
            assert_eq!(jvp(&structure, &old, dt_s, &base, &unit(column))[2], 0.0);
        }
    }
}

#[test]
fn a_structure_tells_the_end_of_its_flows_only_when_asked_and_then_to_the_bit() {
    let plain = ModalStructure::new(
        three_modes(),
        vec![row("bridge", BRIDGE), row("finger", FINGER)],
    )
    .expect("structure");
    assert_eq!(plain.ports().len(), 3);
    for port in 0..3 {
        assert_eq!(plain.end_of_step(port), EndOfStep::Untold);
    }
    let telling = plain.clone().publishing_ends().expect("ends");
    assert_eq!(telling.ports().len(), 5);
    assert_eq!(telling.ports()[3].id().as_str(), "bridge.end");
    assert_eq!(telling.ports()[4].id().as_str(), "finger.end");
    assert_eq!(telling.end_of_step(0), EndOfStep::Told(3));
    assert_eq!(telling.end_of_step(1), EndOfStep::Told(4));
    assert_eq!(
        telling.end_of_step(2),
        EndOfStep::Untold,
        "heat is not a line's flow"
    );
    assert_eq!(
        telling
            .clone()
            .publishing_ends()
            .expect("again")
            .ports()
            .len(),
        5,
        "asking twice adds them once"
    );

    let old = moving();
    let base = [EFFORTS[0], EFFORTS[1], KELVIN];
    for dt_s in durations() {
        let (outputs, new) = proposed(&telling, &old, dt_s, &base);
        let (plain_outputs, plain_new) = proposed(&plain, &old, dt_s, &base);
        assert_eq!(new, plain_new);
        assert_eq!(
            outputs[..3],
            plain_outputs[..],
            "telling the end changes nothing else"
        );
        for (r, weights) in [BRIDGE, FINGER].iter().enumerate() {
            let mut end = 0.0_f64;
            for k in 0..3 {
                end += weights[k] * new[k].velocity_m_sqrt_kg_per_s;
            }
            assert_eq!(outputs[3 + r].to_bits(), end.to_bits());
        }
    }

    let clash = ModalStructure::new(three_modes(), vec![row("a", BRIDGE), row("a.end", FINGER)])
        .expect("two rows");
    assert!(matches!(
        clash.publishing_ends(),
        Err(ModalStructureRefusal::Port(_))
    ));
}
