//! Owner slice: a joint loaded modal basis on the composition kernel, with its loss in a resistor
//! N-port.
//!
//! `fs_nlmodal::loaded::LoadedModalBasis` solves the modes of coupled structures under the added
//! mass of their shared air as one pencil and runs them as one undamped `ModalStructure`. The loss
//! of the set, `Phi^T C Phi`, is `elements::Resistance` on the structure's modal port, joined by a
//! common-flow junction. Here the kernel runs:
//! - one mode under a known added mass and loss, checked against the held scheme's step worked by
//!   hand and against the decay law the basis states for it;
//! - two coupled modes through a free decay, checked against the same scheme summed by hand, with
//!   the ledger closing on the resistance's heat.
//!
//! Run with `cargo test -p fs-system --features modal-slices --test loaded_modal_slice`.

// An undamped mode's heat is compared with zero exactly; the held scheme is written index by index,
// as its equations are.
#![allow(clippy::float_cmp, clippy::needless_range_loop)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_nlmodal::component::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_nlmodal::loaded::LoadedModalBasis;
use fs_ports::{PortKind, StableId};
use fs_system::elements::Resistance;
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};

const TAU: f64 = core::f64::consts::TAU;
const RATE_HZ: u32 = 44_100;
const KIND: PortKind = PortKind::MechanicalForceVelocity;
const AIR_K: f64 = 293.15;
/// The kernel's interface solve stops at a residual of 1e-12 of the terms it sums
/// (`SolveConfig::default`); with the roundoff of one exact modal step, a few eps, one step of the
/// kernel stands within this fraction of the state's scale from the scheme worked by hand.
const PER_STEP: f64 = 1.0e-12 + 16.0 * f64::EPSILON;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 11,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// The block: the loaded structure, the resistance `Phi^T C Phi` on its modal port, and each heat
/// port routed to a sink at the air's temperature.
fn block(
    basis: &LoadedModalBasis,
    loss: &[f64],
    start: Vec<ModalAcousticState>,
) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let structure = basis
        .structure(
            RATE_HZ,
            ModalAcousticTimeBudget::audible_reference(),
            StableId::new("modal").expect("id"),
            KIND,
            Vec::new(),
        )
        .expect("structure");
    let modes = builder.add("modes", structure, start).expect("modes");
    let resistance = Resistance::new(KIND, basis.congruent(loss).expect("loss")).expect("definite");
    let loss_node = builder
        .add_with_tangent("loss", resistance, ())
        .expect("loss");
    let modal = builder.port(modes, "modal").expect("modal");
    let u = builder.port(loss_node, "u").expect("u");
    builder.common_flow(modal, &[u]).expect("junction");
    for (node, sink) in [(modes, "air-modes"), (loss_node, "air-loss")] {
        let heat = builder.port(node, "heat").expect("heat");
        builder
            .heat_route(
                heat,
                &[],
                Some(heat_sink(sink, AIR_K).expect("sink")),
                RouteTemperature::Fixed(AIR_K),
            )
            .expect("route");
    }
    (builder.seal().expect("sealed"), modes, loss_node)
}

fn states(q: &[f64], v: &[f64]) -> Vec<ModalAcousticState> {
    q.iter()
        .zip(v)
        .map(|(&q, &v)| ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        })
        .collect()
}

/// The held scheme summed by hand, apart from the kernel: each loaded mode advances exactly under
/// a force held over the step, `q1 = q0 cos + v0 sin / w + f (1 - cos) / w^2` and
/// `v1 = -q0 w sin + v0 cos + f sin / w`, and the force is the resistance's answer to the step's
/// mean velocity, `f = -C' (q1 - q0) / dt`, solved for `f` by Gaussian elimination.
fn held_step(w: &[f64], loss: &[f64], dt: f64, q: &mut [f64], v: &mut [f64]) {
    let n = w.len();
    let (cos, sin): (Vec<f64>, Vec<f64>) =
        w.iter().map(|&w| ((w * dt).cos(), (w * dt).sin())).unzip();
    // (1 - cos) / w^2 without the cancellation.
    let a: Vec<f64> = w
        .iter()
        .map(|&w| 2.0 * (0.5 * w * dt).sin().powi(2) / (w * w))
        .collect();
    let free: Vec<f64> = (0..n)
        .map(|k| q[k] * (cos[k] - 1.0) + v[k] * sin[k] / w[k])
        .collect();
    // (dt I + C' diag(a)) f = -C' free.
    let mut m = vec![0.0; n * n];
    let mut f = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = loss[i * n + j] * a[j] + if i == j { dt } else { 0.0 };
            f[i] -= loss[i * n + j] * free[j];
        }
    }
    for p in 0..n {
        for r in p + 1..n {
            let factor = m[r * n + p] / m[p * n + p];
            for c in p..n {
                m[r * n + c] -= factor * m[p * n + c];
            }
            f[r] -= factor * f[p];
        }
    }
    for p in (0..n).rev() {
        let mut sum = f[p];
        for c in p + 1..n {
            sum -= m[p * n + c] * f[c];
        }
        f[p] = sum / m[p * n + p];
    }
    for k in 0..n {
        let (q0, v0) = (q[k], v[k]);
        q[k] = q0 * cos[k] + v0 * sin[k] / w[k] + f[k] * a[k];
        v[k] = -q0 * w[k] * sin[k] + v0 * cos[k] + f[k] * sin[k] / w[k];
    }
}

fn kernel_step(system: &mut System, dt: f64, steps: usize) {
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).expect("kernel step");
        }
    });
}

/// One mode under added mass `mu` and loss `c`: the kernel's step is the held scheme's (checked on
/// both columns of its step matrix), and the determinant of that step, the factor by which the
/// energy envelope decays per step, is `(1 - e) / (1 + e)` with
/// `e = (c' dt / 2) sinc^2(w' dt / 2)`, from the factor the basis states.
#[test]
fn one_loaded_mode_decays_by_the_held_schemes_law() {
    let (w0, mu, c) = (TAU * 329.63, 0.012, 6.9);
    let basis = LoadedModalBasis::new(&[w0], &[mu]).expect("definite");
    let w = basis.angular_frequencies_rad_s()[0];
    let loss = basis.congruent(&[c]).expect("loss");
    let dt = 1.0 / f64::from(RATE_HZ);
    // Columns of the step matrix in the scaled state (q w, v): from (1, 0) and from (0, 1).
    let scale = 1.0e-4;
    let mut kernel = [[0.0; 2]; 2];
    let mut hand = [[0.0; 2]; 2];
    for (column, (q0, v0)) in [(scale / w, 0.0), (0.0, scale)].into_iter().enumerate() {
        let (mut system, modes, _) = block(&basis, &[c], states(&[q0], &[v0]));
        kernel_step(&mut system, dt, 1);
        let reached = &system
            .state::<Vec<ModalAcousticState>>(modes)
            .expect("state")[0];
        kernel[0][column] = reached.displacement_m_sqrt_kg * w / scale;
        kernel[1][column] = reached.velocity_m_sqrt_kg_per_s / scale;
        let (mut q, mut v) = ([q0], [v0]);
        held_step(&[w], &loss, dt, &mut q, &mut v);
        hand[0][column] = q[0] * w / scale;
        hand[1][column] = v[0] / scale;
    }
    for row in 0..2 {
        for column in 0..2 {
            assert!(
                (kernel[row][column] - hand[row][column]).abs() <= PER_STEP,
                "({row}, {column}): {} against {}",
                kernel[row][column],
                hand[row][column]
            );
        }
    }
    let det = kernel[0][0] * kernel[1][1] - kernel[0][1] * kernel[1][0];
    let factor = basis.held_loss_factors(dt)[0];
    let e = 0.5 * loss[0] * dt * factor;
    let law = (1.0 - e) / (1.0 + e);
    // The determinant is two products of entries of size at most one, each entry within PER_STEP.
    assert!(
        (det - law).abs() <= 4.0 * PER_STEP,
        "determinant {det} against {law}"
    );
    // Within one step the energy ratio |A x|^2 / |x|^2 depends on the phase of the start x: it
    // ranges over the eigenvalues of A^T A, which for the held scheme are 1 and law^2. Each entry
    // of A^T A is two products of entries within PER_STEP, so within 4 PER_STEP, and its
    // eigenvalues move by at most its 2-norm, below twice its largest entry (Weyl), with a few
    // roundings of the closed form beside.
    let gram = [
        kernel[0][0] * kernel[0][0] + kernel[1][0] * kernel[1][0],
        kernel[0][0] * kernel[0][1] + kernel[1][0] * kernel[1][1],
        kernel[0][1] * kernel[0][1] + kernel[1][1] * kernel[1][1],
    ];
    let (mid, half) = (
        0.5 * (gram[0] + gram[2]),
        (0.25 * (gram[0] - gram[2]).powi(2) + gram[1] * gram[1]).sqrt(),
    );
    let (largest, smallest) = (mid + half, mid - half);
    let spread = 8.0 * PER_STEP + 8.0 * f64::EPSILON;
    assert!(
        (largest - 1.0).abs() <= spread && (smallest - law * law).abs() <= spread,
        "energy ratio over a step: {smallest} to {largest}, law {} to 1",
        law * law
    );
    // The same determinant read as the decay rate of the energy envelope, against the continuous
    // rate c': the rate is c' times the factor times atanh(e) / e, so above c' times the factor by
    // at most e^2 / (3 (1 - e^2)). The determinant's error, 4 PER_STEP, moves the ratio by at most
    // 4 PER_STEP / (2 e det).
    let rate = -(det - 1.0).ln_1p() / dt;
    let over = rate / (loss[0] * factor) - 1.0;
    let bound = e * e / (3.0 * (1.0 - e * e)) + 2.0 * PER_STEP / (e * det);
    assert!(
        over.abs() <= bound,
        "rate {rate}: {over} over, bound {bound}"
    );
    // The factor itself: about 1.8e-4 short of one at 330 Hz and 44.1 kHz.
    assert!(factor < 1.0 && factor > 1.0 - 2.0e-4, "factor {factor}");
}

/// Two coupled modes, near each other as a bar's first mode and its tube's: a free decay through
/// the kernel against the held scheme summed by hand, and the ledger closing on the heat the
/// resistance books.
#[test]
fn two_coupled_modes_follow_the_held_scheme_and_close_the_ledger() {
    let w0 = [TAU * 329.63, TAU * 331.9];
    let added = [0.004, 0.0015, 0.0015, 0.02];
    let c = [5.0, 3.0, 3.0, 60.0];
    let basis = LoadedModalBasis::new(&w0, &added).expect("definite");
    let w = basis.angular_frequencies_rad_s().to_vec();
    let loss = basis.congruent(&c).expect("loss");
    let dt = 1.0 / f64::from(RATE_HZ);
    let (q0, v0) = ([0.0, 0.0], [1.0e-3, 0.0]);
    let (mut system, modes, loss_node) = block(&basis, &c, states(&q0, &v0));
    let (mut q, mut v) = (q0, v0);
    let steps = 4_410_usize;
    let start_j = 0.5 * v0[0] * v0[0];
    let scale = v0[0];
    for chunk in 0..10 {
        kernel_step(&mut system, dt, steps / 10);
        for _ in 0..steps / 10 {
            held_step(&w, &loss, dt, &mut q, &mut v);
        }
        let done = (chunk + 1) * steps / 10;
        let bound = done as f64 * PER_STEP * scale;
        let reached = system
            .state::<Vec<ModalAcousticState>>(modes)
            .expect("state");
        for k in 0..2 {
            let (rq, rv) = (
                reached[k].displacement_m_sqrt_kg,
                reached[k].velocity_m_sqrt_kg_per_s,
            );
            assert!(
                (rq * w[k] - q[k] * w[k]).abs() <= bound,
                "step {done} mode {k}: q {rq} against {}",
                q[k]
            );
            assert!(
                (rv - v[k]).abs() <= bound,
                "step {done} mode {k}: v {rv} against {}",
                v[k]
            );
        }
    }
    let reached = system
        .state::<Vec<ModalAcousticState>>(modes)
        .expect("state");
    let end_j: f64 = reached
        .iter()
        .zip(&w)
        .map(|(s, &w)| {
            0.5 * (s.velocity_m_sqrt_kg_per_s * s.velocity_m_sqrt_kg_per_s
                + w * w * s.displacement_m_sqrt_kg * s.displacement_m_sqrt_kg)
        })
        .sum();
    let ledger = system.ledger();
    let account = |node: NodeId| &ledger.nodes[node.index()];
    let heat_j = account(loss_node).dissipated_j.value();
    assert!(
        heat_j > 0.1 * start_j,
        "a decay worth measuring: {heat_j} of {start_j} J"
    );
    // Energy crossing the junction cancels to the solve tolerance each step, a fraction of the
    // energy that crossed; the modes' own account and the resistance's heat close on it.
    assert!(
        (start_j - end_j - heat_j).abs() <= 1.0e-10 * start_j,
        "{start_j} - {end_j} against {heat_j}"
    );
    assert!(
        account(modes).dissipated_j.value() == 0.0,
        "the modes run undamped"
    );
    assert!((ledger.sink_heat_j() - heat_j).abs() <= 1.0e-12 * heat_j);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10 * start_j);
    assert!(ledger.entropy_generated_j_per_k() > 0.0);
}
