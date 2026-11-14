//! The Hertz patch between two bodies on the composition contract
//! (`fs_tribo::sliding_patch::component`, feature `component`).
//!
//! Run with `cargo test -p fs-tribo --features component --test
//! hertz_patch_interface`. The owner is the crate's `HertzWrenchContact`, and
//! the line is a synthetic body response, not a material.

// Exact float comparisons here assert what the component fixes exactly: the
// owner's own wrench and heat handed on, a pair held at rest, a state that is
// another's to the bit.
#![allow(clippy::float_cmp)]

use std::array::from_fn;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_load_line};
use fs_ports::component::{Component, LoadLine, Proposal, Receipt};
use fs_tribo::TriboError;
use fs_tribo::sliding_patch::component::{HertzPatchInterface, HertzPatchLoad, HertzPatchState};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
    HertzWrenchParameters, HertzWrenchPhase, HertzWrenchStep,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters, HertzSlidingWrench};

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

const NORMAL_N: f64 = 20.0;
const RADIUS_M: f64 = 0.02;
const DT_S: f64 = 0.01;

/// The owner's budgets in the module's own Hertz test
/// (`fs-couple`, tests/articulated_friction_heat/hertz.rs).
fn solve() -> HertzWrenchConfig {
    HertzWrenchConfig {
        maximum_static_iterations: 96,
        maximum_iterations: 256,
        maximum_evaluations: 512,
        force_tolerance_n: 2e-9,
        moment_tolerance_nm: 2e-10,
        displacement_tolerance_m: 2e-13,
        energy_tolerance_j: 2e-11,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 100.0,
        maximum_spin_moment_nm: 10.0,
    }
}

fn load() -> HertzPatchLoad {
    HertzPatchLoad {
        normal_n: NORMAL_N,
        radius_m: RADIUS_M,
        solve: solve(),
    }
}

fn hertz(static_mu: f64, kinetic_mu: f64, radius_m: f64) -> Result<HertzWrenchContact, TriboError> {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: NORMAL_N,
        contact_radius_m: radius_m,
        static_mu,
        kinetic_mu,
    })
    .map_err(|error| match error {
        HertzWrenchError::Law(error) => error,
        _ => TriboError::InvalidInput { field: "test_law" },
    })
}

/// A coefficient that falls with temperature, `0.4 - 0.01 (T - 300 K)`, as in
/// the module's test.
fn mu(kelvin: f64) -> f64 {
    0.4 - 0.01 * (kelvin - 300.0)
}

fn weakening(kelvin: f64) -> Result<HertzWrenchContact, TriboError> {
    hertz(mu(kelvin), mu(kelvin), RADIUS_M)
}

/// A body pair's line over a step: half the step times a symmetric positive
/// definite mobility [1/kg], with slip and spin coupled.
fn body_line() -> [f64; 9] {
    [0.9, 0.1, 0.05, 0.1, 0.8, -0.02, 0.05, -0.02, 0.6].map(|m| 0.5 * DT_S * m)
}

/// A unit mass on every coefficient, uncoupled.
fn unit_line() -> [f64; 9] {
    let c = 0.5 * DT_S;
    [c, 0.0, 0.0, 0.0, c, 0.0, 0.0, 0.0, c]
}

/// The native Hertz patch at a slip rate and a spin rate.
fn native(friction_coefficient: f64, slip_m_s: [f64; 2], spin_rad_s: f64) -> HertzSlidingWrench {
    HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: NORMAL_N,
        contact_radius_m: RADIUS_M,
        friction_coefficient,
    })
    .unwrap()
    .evaluate(slip_m_s, spin_rad_s)
    .unwrap()
}

/// One step on a line outside a kernel: the candidate, the receipt and the
/// outputs (the three efforts, then the heat rate).
fn on_line<F>(
    patch: &HertzPatchInterface<F>,
    state: &HertzPatchState,
    inputs: [f64; 4],
    line: &[f64; 9],
) -> Result<(HertzPatchState, Receipt, [f64; 4]), TriboError>
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError>,
{
    let mut next = *state;
    let mut outputs = [f64::NAN; 4];
    let proposal = with_cx(|cx| {
        patch.propose_on_line(state, DT_S, &inputs, line, &mut next, &mut outputs, cx)
    })?;
    let Proposal::Step(receipt) = proposal else {
        panic!("the patch takes the step or refuses");
    };
    Ok((next, receipt, outputs))
}

fn field<T: core::fmt::Debug>(result: Result<T, TriboError>) -> &'static str {
    match result {
        Err(TriboError::InvalidInput { field }) => field,
        other => panic!("expected a refusal by field, got {other:?}"),
    }
}

/// `|a - b| <= tolerance`, with both printed to the digit.
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} vs {b:.17e}");
}

/// Under a prescribed motion (no line) the patch is the native Hertz patch at
/// that motion: its force, its moment and its power, for slip, spin and both.
#[test]
fn prescribed_slip_and_spin_are_the_native_patch() {
    let patch = HertzPatchInterface::new(weakening, load()).unwrap();
    assert_eq!(patch.line_ports(), &[0]);
    assert_eq!(patch.ports()[0].width(), 3);
    let stale = patch.start([9.0; 3]).unwrap();
    for released in [
        [0.3, -0.4, 0.0],
        [0.0, 0.0, RADIUS_M * 50.0],
        [0.3, -0.4, RADIUS_M * 50.0],
    ] {
        let state = patch.start(released).unwrap();
        let inputs = [released[0], released[1], released[2], 300.0];
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: DT_S,
            inputs: &inputs,
        };
        let observed = with_cx(|cx| check_component(&patch, &case, cx));
        let Ok(Observed::Stepped {
            receipt, outputs, ..
        }) = observed
        else {
            panic!("{observed:?}");
        };
        let wrench = native(
            mu(300.0),
            [released[0], released[1]],
            released[2] / RADIUS_M,
        );
        let force_scale = mu(300.0) * NORMAL_N;
        close(outputs[0], -wrench.force_n[0], 1e-12 * force_scale);
        close(outputs[1], -wrench.force_n[1], 1e-12 * force_scale);
        close(
            outputs[2],
            -wrench.spin_moment_nm / RADIUS_M,
            1e-12 * force_scale,
        );
        close(
            receipt.dissipated_j,
            DT_S * wrench.dissipated_power_w,
            1e-12 * DT_S * wrench.dissipated_power_w,
        );
        assert!(receipt.dissipated_j > 0.0);
        assert_eq!(outputs[3], -receipt.dissipated_j / DT_S);
        assert_eq!(
            receipt.entropy_generated_j_per_k,
            receipt.dissipated_j / 300.0
        );
        assert_eq!(receipt.stored_change_j, 0.0);
        // A given motion is held over the step.
        let (next, ..) = on_line(&patch, &state, inputs, &[0.0; 9]).unwrap();
        assert_eq!(next.motion_m_s(), released);
    }
}

/// On a body pair's line the effort is the owner's own step at the released
/// motion and the line, to the bit, and so is the heat. The port's flows are
/// the owner's motion to its compatibility tolerance, and the step ends at
/// twice the owner's mean less the start.
#[test]
fn on_a_line_the_effort_is_the_owners_step() {
    let patch = HertzPatchInterface::new(weakening, load()).unwrap();
    let line = body_line();
    let start = [0.3, -0.4, RADIUS_M * 20.0];
    let state = patch.start(start).unwrap();
    let stale = patch.start([0.0; 3]).unwrap();
    let inputs = [start[0], start[1], start[2], 310.0];
    let case = Case {
        state: &state,
        stale: &stale,
        dt_s: DT_S,
        inputs: &inputs,
    };
    let checked = with_cx(|cx| check_load_line(&patch, &case, &line, cx));
    let Ok(on) = checked else {
        panic!("{checked:?}");
    };
    let Observed::Stepped {
        receipt, outputs, ..
    } = on.observed
    else {
        panic!("{:?}", on.observed);
    };
    let owner = weakening(310.0)
        .unwrap()
        .advance_phase(
            &HertzWrenchInput {
                free_displacement_m: [DT_S * start[0], DT_S * start[1]],
                free_spin_angle_rad: DT_S * start[2] / RADIUS_M,
                scaled_compliance_m_per_n: from_fn(|i| from_fn(|j| DT_S * line[3 * i + j])),
            },
            HertzWrenchPhase::Kinetic,
            solve(),
            || false,
        )
        .unwrap();
    assert_eq!(owner.mode, HertzWrenchMode::Sliding);
    assert_eq!(
        outputs[..3],
        [
            -owner.force_n[0],
            -owner.force_n[1],
            -owner.spin_moment_nm / RADIUS_M
        ]
    );
    assert_eq!(receipt.dissipated_j, owner.dissipated_energy_j);
    assert_eq!(outputs[3], -owner.dissipated_energy_j / DT_S);
    assert!(receipt.dissipated_j > 0.0);
    assert!(receipt.numerical_loss_j.abs() <= solve().energy_tolerance_j);
    let moved = [
        owner.slip_displacement_m[0] / DT_S,
        owner.slip_displacement_m[1] / DT_S,
        RADIUS_M * owner.spin_angle_rad / DT_S,
    ];
    for (flow, moved) in on.flows.iter().zip(moved) {
        close(DT_S * flow, DT_S * moved, solve().displacement_tolerance_m);
    }
    let (next, ..) = on_line(&patch, &state, inputs, &line).unwrap();
    assert_eq!(next.motion_m_s(), from_fn(|i| 2.0 * moved[i] - start[i]));
}

/// The phase is the module's: a pair that starts at rest is held by the
/// static set and ends exactly at rest; the same push on a pair that starts
/// moving is beyond the kinetic set, and the pair slides at it.
#[test]
fn a_pair_at_rest_is_held_by_the_static_set_and_a_moving_one_is_not() {
    let patch = HertzPatchInterface::new(|_| hertz(0.5, 0.3, RADIUS_M), load()).unwrap();
    let line = unit_line();
    // Holding the pair takes 0.4 of the load: inside the static set (0.5),
    // beyond the kinetic one (0.3).
    let push_n = 0.4 * NORMAL_N;
    let inputs = [line[0] * push_n, 0.0, 0.0, 300.0];
    let (held, receipt, outputs) =
        on_line(&patch, &patch.start([0.0; 3]).unwrap(), inputs, &line).unwrap();
    assert_eq!(held.motion_m_s(), [0.0; 3]);
    close(outputs[0], push_n, 1e-12 * push_n);
    close(outputs[1], 0.0, 1e-12 * push_n);
    close(outputs[2], 0.0, 1e-12 * push_n);
    assert_eq!(receipt.dissipated_j, 0.0);
    let moving = patch.start([1e-3, 0.0, 0.0]).unwrap();
    let (slid, receipt, outputs) = on_line(&patch, &moving, inputs, &line).unwrap();
    close(outputs[0], 0.3 * NORMAL_N, 1e-12 * NORMAL_N);
    close(outputs[1], 0.0, 1e-12 * NORMAL_N);
    close(outputs[2], 0.0, 1e-12 * NORMAL_N);
    assert!(receipt.dissipated_j > 0.0);
    assert!(slid.motion_m_s()[0] > 0.0);
}

/// A patch with no friction takes no power and makes no heat; the motion
/// passes through it.
#[test]
fn a_frictionless_patch_takes_no_power() {
    let patch = HertzPatchInterface::new(|_| hertz(0.0, 0.0, RADIUS_M), load()).unwrap();
    let start = [0.3, 0.0, RADIUS_M * 10.0];
    let inputs = [start[0], start[1], start[2], 300.0];
    let (next, receipt, outputs) =
        on_line(&patch, &patch.start(start).unwrap(), inputs, &body_line()).unwrap();
    assert_eq!(outputs, [0.0; 4]);
    assert_eq!(receipt.dissipated_j, 0.0);
    assert_eq!(receipt.numerical_loss_j, 0.0);
    for (end, start) in next.motion_m_s().into_iter().zip(start) {
        close(end, start, 1e-13 * start.abs());
    }
}

/// The law is asked at the temperature the heat port receives.
#[test]
fn the_law_is_asked_at_the_interface_temperature() {
    let patch = HertzPatchInterface::new(weakening, load()).unwrap();
    let state = patch.start([0.5, 0.0, 0.0]).unwrap();
    for kelvin in [300.0, 310.0, 330.0] {
        let (_, receipt, outputs) =
            on_line(&patch, &state, [0.5, 0.0, 0.0, kelvin], &[0.0; 9]).unwrap();
        close(-outputs[3], mu(kelvin) * NORMAL_N * 0.5, 1e-12 * NORMAL_N);
        assert_eq!(
            receipt.entropy_generated_j_per_k,
            receipt.dissipated_j / kelvin
        );
    }
}

/// A refused step, by field, with nothing written to the candidate or the
/// outputs.
fn refusal<F>(
    patch: &HertzPatchInterface<F>,
    start: [f64; 3],
    inputs: [f64; 4],
    line: &[f64; 9],
) -> &'static str
where
    F: Fn(f64) -> Result<HertzWrenchContact, TriboError>,
{
    let state = patch.start(start).unwrap();
    let sentinel = patch.start([7.0; 3]).unwrap();
    let mut next = sentinel;
    let mut outputs = [f64::NAN; 4];
    let result = with_cx(|cx| {
        patch.propose_on_line(&state, DT_S, &inputs, line, &mut next, &mut outputs, cx)
    });
    assert_eq!(next, sentinel, "a refusal writes no candidate");
    assert!(
        outputs.iter().all(|v| v.is_nan()),
        "a refusal writes no output"
    );
    field(result)
}

#[test]
fn refusals_name_their_field_and_write_nothing() {
    let patch = HertzPatchInterface::new(weakening, load()).unwrap();
    let slide = [0.3, -0.4, 0.0, 300.0];
    let unit = unit_line();
    // Another radius would rescale the port's third coefficient.
    let wider = HertzPatchInterface::new(|_| hertz(0.4, 0.4, 0.03), load()).unwrap();
    assert_eq!(
        refusal(&wider, [0.3, -0.4, 0.0], slide, &unit),
        "hertz_patch_law_geometry"
    );
    // A moving pair that the kinetic set holds at no mean slip would end the
    // step sliding the other way: the module's coarse reversal.
    let stops = HertzPatchInterface::new(|_| hertz(0.5, 0.3, RADIUS_M), load()).unwrap();
    assert_eq!(
        refusal(&stops, [0.01, 0.0, 0.0], [0.01, 0.0, 0.0, 300.0], &unit),
        "hertz_patch_reversal"
    );
    // A slide on a line with an antisymmetric part, a thousandth of the
    // line's coefficient or of rounding size, is taken on the full line; a
    // part the owner's repeat cannot settle is refused (the precessing line,
    // below).
    let c = unit[0];
    let state = patch.start([0.3, -0.4, 0.0]).unwrap();
    for size in [1e-3, 1e-12] {
        let mut skewed = unit;
        skewed[1] = size * c;
        skewed[3] = -size * c;
        let (_, receipt, _) = on_line(&patch, &state, slide, &skewed).unwrap();
        assert!(receipt.dissipated_j > 0.0);
        assert!(receipt.numerical_loss_j.abs() <= solve().energy_tolerance_j);
    }
    // Inputs.
    assert_eq!(
        refusal(&patch, [0.3, -0.4, 0.0], [0.3, -0.4, 0.0, 0.0], &unit),
        "interface_temperature_k"
    );
    assert_eq!(
        refusal(
            &patch,
            [0.3, -0.4, 0.0],
            [f64::NAN, -0.4, 0.0, 300.0],
            &unit
        ),
        "hertz_patch_inputs"
    );
    // The owner's effort budget.
    let mut small = load();
    small.solve.maximum_force_n = 1.0;
    let budgeted = HertzPatchInterface::new(weakening, small).unwrap();
    assert_eq!(
        refusal(&budgeted, [0.3, -0.4, 0.0], slide, &[0.0; 9]),
        "hertz_patch_effort_budget"
    );
}

/// The owner's step at the free motion `dt (released - K effort)` on the
/// symmetric part of `line`, `K` its antisymmetric part: what the component
/// hands the owner, written here from the line itself.
fn owner_on(line: &[f64; 9], released: [f64; 3], effort: [f64; 3]) -> HertzWrenchStep {
    let skew: [f64; 3] = from_fn(|i| {
        (0..3)
            .map(|j| 0.5 * (line[3 * i + j] - line[3 * j + i]) * effort[j])
            .sum()
    });
    weakening(300.0)
        .unwrap()
        .advance_phase(
            &HertzWrenchInput {
                free_displacement_m: [
                    DT_S * (released[0] - skew[0]),
                    DT_S * (released[1] - skew[1]),
                ],
                free_spin_angle_rad: DT_S * (released[2] - skew[2]) / RADIUS_M,
                scaled_compliance_m_per_n: from_fn(|i| {
                    from_fn(|j| DT_S * 0.5 * (line[3 * i + j] + line[3 * j + i]))
                }),
            },
            HertzWrenchPhase::Kinetic,
            solve(),
            || false,
        )
        .unwrap()
}

/// A line with a large antisymmetric slip block, as the gyroscopic coupling
/// of a body spinning about the patch normal gives. A pure spin has no force,
/// so the block does nothing to it: one owner step, the symmetric line's to
/// the bit. A slide across it is taken on the full line: its motion is
/// `dt (released - C effort)` with the whole `C`, its effort is the owner's
/// fixed point there, and it costs more owner steps. A block so large that
/// the owner's repeat does not contract is refused, and writes nothing.
#[test]
fn a_precessing_line_takes_a_spin_and_a_slide_on_the_full_line() {
    let patch = HertzPatchInterface::new(weakening, load()).unwrap();
    let symmetric = unit_line();
    let c = symmetric[0];
    let mut precessing = symmetric;
    precessing[1] = 0.05 * c;
    precessing[3] = -0.05 * c;
    let spin = [0.0, 0.0, RADIUS_M * 20.0];
    let state = patch.start(spin).unwrap();
    let inputs = [spin[0], spin[1], spin[2], 300.0];
    let (next, receipt, outputs) = on_line(&patch, &state, inputs, &precessing).unwrap();
    let (plain_next, plain_receipt, plain_outputs) =
        on_line(&patch, &state, inputs, &symmetric).unwrap();
    assert_eq!(outputs[..2], [0.0; 2]);
    assert!(outputs[2].abs() > 0.0 && receipt.dissipated_j > 0.0);
    assert_eq!(outputs.map(f64::to_bits), plain_outputs.map(f64::to_bits));
    assert_eq!(next, plain_next);
    assert_eq!(receipt, plain_receipt);
    assert_eq!(next.owner_steps(), 1);

    // The slide is taken, after more than one owner step, and the block has
    // turned its effort from the symmetric line's.
    let released = [0.3, -0.4, 0.0];
    let state = patch.start(released).unwrap();
    let inputs = [released[0], released[1], released[2], 300.0];
    let (next, receipt, outputs) = on_line(&patch, &state, inputs, &precessing).unwrap();
    let (plain_next, _, plain_outputs) = on_line(&patch, &state, inputs, &symmetric).unwrap();
    let effort = [outputs[0], outputs[1], outputs[2]];
    println!(
        "precessing slide: effort {effort:?} against the symmetric line's {:?}, owner steps {} against {}",
        &plain_outputs[..3],
        next.owner_steps(),
        plain_next.owner_steps()
    );
    assert!(receipt.dissipated_j > 0.0 && receipt.residual <= receipt.tolerance);
    assert_eq!(plain_next.owner_steps(), 1);
    assert!(next.owner_steps() > 1);
    assert_ne!(effort[..2], plain_outputs[..2]);
    // Its motion over the step, the owner's (the step ends at twice it less
    // the start), is the full line's at that effort, written from the whole
    // `C` here.
    for (i, (end, start)) in next.motion_m_s().into_iter().zip(released).enumerate() {
        let owner_m = DT_S * 0.5 * (end + start);
        let full_m = DT_S
            * (start
                - (0..3)
                    .map(|j| precessing[3 * i + j] * effort[j])
                    .sum::<f64>());
        close(owner_m, full_m, solve().displacement_tolerance_m);
    }
    // One more owner step at that effort moves it by no more than the
    // owner's own tolerances.
    let again = owner_on(&precessing, released, effort);
    close(-again.force_n[0], effort[0], solve().force_tolerance_n);
    close(-again.force_n[1], effort[1], solve().force_tolerance_n);
    close(
        -again.spin_moment_nm,
        RADIUS_M * effort[2],
        solve().moment_tolerance_nm,
    );

    // A block twenty times the line's own coefficient: the owner's force
    // turns across the free motion by more than the free motion, and the
    // repeat's gain, the owner's force change per free motion across the
    // slide times dt |K|, is above one. Refused; nothing written.
    let mut spinning = symmetric;
    spinning[1] = 20.0 * c;
    spinning[3] = -20.0 * c;
    let at = owner_on(&symmetric, released, [0.0; 3]);
    let across = [0.8, 0.6, 0.0];
    let nudge_m = 1e-9;
    let nudged = owner_on(
        &symmetric,
        from_fn(|i| released[i] + across[i] * nudge_m / DT_S),
        [0.0; 3],
    );
    let stiffness_n_m = (0..2)
        .map(|i| (nudged.force_n[i] - at.force_n[i]) * across[i])
        .sum::<f64>()
        .abs()
        / nudge_m;
    let gain = stiffness_n_m * DT_S * 20.0 * c;
    println!("non-contracting block: gain {gain} ({stiffness_n_m} N/m across the slide)");
    assert!(gain > 1.0, "{gain}");
    assert_eq!(
        refusal(&patch, released, inputs, &spinning),
        "hertz_patch_line_skew"
    );
}

#[test]
fn a_step_the_caller_cancels_is_refused_and_writes_nothing() {
    let patch = HertzPatchInterface::new(weakening, load()).unwrap();
    let state = patch.start([0.3, -0.4, 0.0]).unwrap();
    let sentinel = patch.start([7.0; 3]).unwrap();
    let (mut next, mut outputs) = (sentinel, [f64::NAN; 4]);
    let gate = CancelGate::new_clock_free();
    gate.request();
    let pool = ArenaPool::new(ArenaConfig::default());
    let result = pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        patch.propose_on_line(
            &state,
            DT_S,
            &[0.3, -0.4, 0.0, 300.0],
            &unit_line(),
            &mut next,
            &mut outputs,
            &cx,
        )
    });
    assert_eq!(field(result), "hertz_patch_cancelled");
    assert_eq!(next, sentinel, "a refusal writes no candidate");
    assert!(
        outputs.iter().all(|v| v.is_nan()),
        "a refusal writes no output"
    );
}

#[test]
fn a_patch_is_admitted_with_a_load_and_budgets_the_owner_takes() {
    let with = |change: fn(&mut HertzPatchLoad)| {
        let mut load = load();
        change(&mut load);
        field(HertzPatchInterface::new(weakening, load))
    };
    assert_eq!(with(|load| load.normal_n = -1.0), "normal_force");
    assert_eq!(with(|load| load.radius_m = 0.0), "contact_radius_m");
    assert_eq!(
        with(|load| load.solve.maximum_iterations = 0),
        "iteration_budget"
    );
    assert_eq!(
        field(
            HertzPatchInterface::new(weakening, load())
                .unwrap()
                .start([f64::NAN, 0.0, 0.0])
        ),
        "motion_m_s"
    );
}
