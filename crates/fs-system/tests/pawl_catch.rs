//! The catch (`frankensim-rc-1109-root-3gc12.13.1`, O2b and O3b): a platform that falls on a rail
//! is caught by a pawl on a ladder of tooth plates, from generic elements, against hand
//! calculations; and the same platform hung on a rope, cut loose by a latch that frees the rope.
//!
//! The platform (`M = 2 kg`) slides down a rail whose shoes press it with `N = 5 N` (a `Slider`,
//! `mu_s = 0.3`, `mu_k = 0.2`) under its weight, so it falls at `a = g - mu_k N / M`. A pawl turns
//! about a pivot on the platform, its centre of mass at the pivot, pushed out by a leaf spring
//! (rest at turn 0) from its retracted turn 0.6 rad to a stop at 0.3 rad, where its tip is
//! straight below the pivot: a strut. The stop is a heel against a face parallel to the rail, so
//! its gap does not depend on the fall. Four thin tooth plates, each a bounded `PoseGap` under a
//! lossy Hertz contact, stand beside the tip's path. Every heat port joins one heat capacity.
//!
//! Hand calculations: the swing `theta(t) = 0.6 cos(omega t)` gives the steps at which the tip
//! enters the plates' span and the heel meets the stop; the fall gives the plate the tip is above
//! at entry (it passes under the first, behind its band) and the step of contact; the impact's
//! first integral `M v dv/d delta = M a - K delta^1.5 (1 + chi v)`, integrated in `delta`
//! (fourth-order Runge-Kutta, checked against the closed form at `a = 0`), gives the deepest
//! penetration, the peak load on the pawl and the plate's heat at the arrest. Until the catch the
//! platform falls at `a`; up to the arrest the rail's heat is `mu_k N` times the fall; at the
//! arrest the stores of every node and the heat sum to nothing.
//!
//! The rope (O3b) is tied to the pawl 5 mm beside the pivot, on the retracted side, and runs up to
//! a fitting held by a latch (an open mechanical `IdealSwitch`). It holds the platform's weight,
//! and its torque holds the pawl on a retract stop (a second heel) against the leaf. It is a
//! lossless one-sided linear spring on a `PoseGap`, `K = E A / L`: `E` = 8.1e6 psi, AISI's rated
//! modulus of 8x19 fibre-core elevator rope as US 4,887,422 quotes it (the area it is on is not
//! stated there), on the metallic area `0.357 d^2` of a 2 mm rope (PFEIFER, 6x19M-FC); the
//! fitting is the rope's steel. The length, 0.985 m, is chosen so that every event of the cut is
//! at least 0.1 of a step from a step's end at both step lengths, as asserted. The platform hangs
//! at rest for 32 steps; then a scheduled event closes the latch, which frees the fitting. The
//! release to the rope's slack is the same laws by RK4 at a hundredth of a step, the platform held
//! until the first step whose mean rope force puts the rail's load past `mu_s N`, as the `Slider`
//! decides; from the slack, the swing and the fall above start from the release's state. The
//! patent's toggles are assembly edits: without the leaf, or without the pawl (the rope tied to
//! the platform), no plate closes and the platform falls at `a` from the slack past the plates.
//!
//! No claim: the pawl's centre of mass is at its pivot (no coupled mass matrix), the tip has no
//! friction, the plates have no flanks or backs as faces, one pawl. The rope is vertical, its mass
//! is lumped in the fitting (no waves), the cut takes no time and the falling fitting meets
//! nothing (it stays above the platform through the arrest, asserted); the drum is O4.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::Obstacle;
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_system::elements::{
    Force, HeatCapacity, IdealSwitch, Mass, Motion, Pose, PoseGap, Slider, Slip, SourceAccount,
    Spring, Waveform,
};
use fs_system::{Ledger, NodeId, SolveConfig, System, SystemBuilder};

const G_M_S2: f64 = 9.806_65;
const MASS_KG: f64 = 2.0;
const NORMAL_N: f64 = 5.0;
const MU_S: f64 = 0.3;
const MU_K: f64 = 0.2;
/// Pivot to tip and pivot to heel [m].
const ARM_M: f64 = 0.02;
const HEEL_M: f64 = 0.005;
const PAWL_KG_M2: f64 = 2.0e-6;
const LEAF_N_M_RAD: f64 = 0.05;
/// The turn at the cut and at the stop, where the tip is straight below the pivot [rad].
const RETRACTED: f64 = 0.6;
const STRUT: f64 = 0.3;
/// The plates span from `INNER_M` beside the pivot's line, on the side the tip comes from, to
/// `OUTER_M` past it, with a band `DEPTH_M` deep.
const INNER_M: f64 = 2.0e-3;
const OUTER_M: f64 = 5.0e-3;
const DEPTH_M: f64 = 3.0e-4;
/// The first plate's top below the tip at the cut, and the pitch of the plates [m].
const FIRST_M: f64 = 2.0e-4;
const PITCH_M: f64 = 2.5e-3;
const TEETH: usize = 4;
/// Hertz stiffness [N/m^1.5] and Hunt-Crossley loss [s/m] of a plate and of the stops.
const TOOTH: [f64; 2] = [1.0e10, 2.0];
const STOP: [f64; 2] = [1.0e11, 50.0];
/// The rope's modulus [Pa] (8.1e6 psi) on its metallic area `METALLIC d^2`, its diameter and
/// length [m], its steel [kg/m^3], and its tie on the pawl from the pivot [m].
const ROPE_PA: f64 = 8.1e6 * 6_894.757_293_168_361;
const METALLIC: f64 = 0.357;
const ROPE_M: [f64; 2] = [2.0e-3, 0.985];
const STEEL_KG_M3: f64 = 7_850.0;
const LEVER_M: f64 = 5.0e-3;
/// A binary step (2^-20 s), so that the cut, after 32 of them, ends a step at `dt / 2` too.
const DT_S: f64 = 1.0 / 1_048_576.0;
const CUT_S: f64 = 32.0 * DT_S;
/// The release's RK4 steps in a step.
const SUBSTEPS: usize = 100;
/// A run's events: the rail's break-away, the tip's entry into the plates' span, then the first
/// change of each gap's contact from its start, open or closed: the stop (`2`), the plates, the
/// retract stop (the pawl leaves it) and the rope (slack).
const EVENTS: usize = 5 + TEETH;
const RETRACT: usize = 3 + TEETH;
const ROPE: usize = 4 + TEETH;
/// Steps of the reference's quadrature; half as many give its error.
const QUADRATURE: usize = 16_000;
const SOLVE: ContactSolve = ContactSolve {
    max_iterations: 100,
    maximum_force_n: 1.0e6,
    maximum_penetration_m: 1.0e-3,
    force_absolute_tolerance_n: 1.0e-9,
    force_relative_tolerance: 1.0e-12,
};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 8,
            tile: 9,
            iteration: 10,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// The tip from the pivot at the turn `theta` [m]: beside it (towards the retracted side) and up.
fn tip(theta: f64) -> [f64; 2] {
    let (sin, cos) = (theta - STRUT).sin_cos();
    [ARM_M * sin, -ARM_M * cos]
}

/// The height of plate `k`'s top, from the pivot at the cut [m].
fn top_m(k: usize) -> f64 {
    tip(RETRACTED)[1] - FIRST_M - PITCH_M * k as f64
}

fn contact(name: &str, alpha: f64, [stiffness, loss]: [f64; 2]) -> ContactInterface {
    let law = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        stiffness,
        alpha,
        name.into(),
    );
    let law = law.unwrap().with_internal_loss(loss).unwrap();
    ContactInterface::new(law, SOLVE).unwrap()
}

/// The peak of a smooth series sampled at equal steps: the parabola through its largest sample
/// and the two beside it.
fn vertex(y: &[f64]) -> f64 {
    let i = (1..y.len() - 1)
        .max_by(|&i, &j| y[i].total_cmp(&y[j]))
        .unwrap();
    let (a, b, c) = (y[i - 1], y[i], y[i + 1]);
    b + (c - a).powi(2) / (8.0 * (2.0 * b - a - c))
}

/// `w = v^2 / 2` against `delta`, `dw/d delta = (F - K delta^1.5 (1 + chi sqrt(2 w))) / M`, by
/// RK4 in `steps` steps over `span`, to `w = 0`: the deepest penetration and the peak load.
fn impact(force_n: f64, v0: f64, span: f64, steps: usize) -> [f64; 2] {
    let [k, chi] = TOOTH;
    let load = |d: f64, w: f64| k * d.powf(1.5) * (1.0 + chi * (2.0 * w.max(0.0)).sqrt());
    let slope = |d: f64, w: f64| (force_n - load(d, w)) / MASS_KG;
    let (h, mut d, mut w, mut loads) = (span / steps as f64, 0.0, 0.5 * v0 * v0, vec![0.0]);
    loop {
        let k1 = slope(d, w);
        let k2 = slope(d + 0.5 * h, w + 0.5 * h * k1);
        let k3 = slope(d + 0.5 * h, w + 0.5 * h * k2);
        let next = w + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + slope(d + h, w + h * k3));
        if next <= 0.0 {
            return [d + h * w / (w - next), vertex(&loads)];
        }
        (d, w) = (d + h, next);
        loads.push(load(d, w));
    }
}

/// The hand calculations.
struct Plan {
    accel_m_s2: f64,
    entry_s: f64,
    stop_s: f64,
    /// The plate the tip is above at entry.
    plate: usize,
    contact_s: f64,
    /// Deepest penetration [m], peak load [N] and the plate's heat at the arrest [J], with the
    /// quadrature's error of each.
    impact: [f64; 3],
    quadrature: [f64; 3],
}

/// The swing and the fall from `[t, theta, theta', s, s']`, the pawl free of the rope: `theta`
/// is `A cos(omega (t' - t) + phase)`, the fall is at `a`.
fn plan([t0, theta, rate, s0, u0]: [f64; 5]) -> Plan {
    let accel_m_s2 = G_M_S2 - MU_K * NORMAL_N / MASS_KG;
    let omega = (LEAF_N_M_RAD / PAWL_KG_M2).sqrt();
    let (amplitude, phase) = (theta.hypot(rate / omega), (-rate / omega).atan2(theta));
    let swing = |turn: f64| t0 + ((turn / amplitude).acos() - phase) / omega;
    let entry = STRUT + (INNER_M / ARM_M).asin();
    let entry_s = swing(entry);
    let fall = |t: f64| s0 + u0 * (t - t0) + 0.5 * accel_m_s2 * (t - t0).powi(2);
    let height = -fall(entry_s) + tip(entry)[1];
    let plate = (0..TEETH).find(|&k| height >= top_m(k)).unwrap();
    // The plates above it have the tip behind them, under their band.
    assert!((0..plate).all(|k| height < top_m(k) - DEPTH_M));
    let drop_m = -ARM_M - top_m(plate) - s0;
    let tau = (u0.mul_add(u0, 2.0 * accel_m_s2 * drop_m).sqrt() - u0) / accel_m_s2;
    let v0 = u0 + accel_m_s2 * tau;
    let [k, chi] = TOOTH;
    let free = (2.5 * MASS_KG * (chi * v0 - (chi * v0).ln_1p()) / (k * chi * chi)).powf(0.4);
    let reference = |force_n: f64, steps: usize| {
        let [deepest_m, peak_n] = impact(force_n, v0, 2.0 * free, steps);
        let heat_j = 0.5 * MASS_KG * v0 * v0 + force_n * deepest_m - k * deepest_m.powf(2.5) / 2.5;
        [deepest_m, peak_n, heat_j]
    };
    let [fine, coarse] = [QUADRATURE, QUADRATURE / 2].map(|n| reference(0.0, n)[0]);
    assert!(
        (fine - free).abs() <= (coarse - fine).abs(),
        "{fine} {coarse} {free}"
    );
    let weight = MASS_KG * accel_m_s2;
    let [fine, coarse] = [QUADRATURE, QUADRATURE / 2].map(|n| reference(weight, n));
    Plan {
        accel_m_s2,
        entry_s,
        stop_s: swing(STRUT),
        plate,
        contact_s: t0 + tau,
        impact: fine,
        quadrature: [0, 1, 2].map(|i| (coarse[i] - fine[i]).abs()),
    }
}

/// The platform hung at rest: the rope's stiffness [N/m], the fitting's mass [kg], the pawl's
/// turn on the retract stop, where the rope's torque meets the stop's and the leaf's (by
/// bisection) [rad], and the height above the pivot of the rope's face, at which the rope holds
/// the platform's weight [m].
fn hung(leaf: bool) -> [f64; 4] {
    let area_m2 = METALLIC * ROPE_M[0] * ROPE_M[0];
    let stiffness = ROPE_PA * area_m2 / ROPE_M[1];
    let (tension, spring) = (MASS_KG * G_M_S2, if leaf { LEAF_N_M_RAD } else { 0.0 });
    let torque = |theta: f64| {
        let (sin, cos) = (theta - RETRACTED).sin_cos();
        let stop = STOP[0] * (HEEL_M * sin).max(0.0).powf(1.5);
        (tension * LEVER_M - stop * HEEL_M) * cos - spring * theta
    };
    let (mut lo, mut hi) = (RETRACTED, RETRACTED + 1.0e-3);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        if torque(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let face_m = tension / stiffness + LEVER_M * (lo - RETRACTED).sin();
    [stiffness, STEEL_KG_M3 * area_m2 * ROPE_M[1], lo, face_m]
}

/// The release, from the cut to the rope's slack.
struct Release {
    /// The step of the break-away, and the fraction of it at which the force the rail would
    /// hold, as step means, passes `mu_s N`.
    step: usize,
    fraction: f64,
    /// When the pawl leaves the retract stop and the rope goes slack, from the cut [s].
    leave_s: f64,
    slack_s: f64,
    /// At the first substep with the rope slack: its time, the turn and the fall with rates.
    end: [f64; 5],
}

/// The release by RK4 at `dt / SUBSTEPS`: the fitting falls under the rope and its weight, the
/// pawl turns under the rope, the stop and the leaf, and the platform is held until the first
/// step whose mean rope force (the mean of its step-end tensions: the rope is linear) puts the
/// rail's load past `mu_s N`, then slides at `mu_k N` from that step's start, as the `Slider`
/// decides.
fn release(dt: f64) -> Release {
    let [stiffness, fitting_kg, theta, face_m] = hung(true);
    let h = dt / SUBSTEPS as f64;
    // The fitting's fall, the platform's, the turn, each with its rate.
    let closure = |y: &[f64; 6]| face_m + y[2] - LEVER_M * (y[4] - RETRACTED).sin() - y[0];
    let slope = |y: &[f64; 6], slides: bool| -> [f64; 6] {
        let tension = stiffness * closure(y).max(0.0);
        let (sin, cos) = (y[4] - RETRACTED).sin_cos();
        let ([k, chi], press) = (STOP, (HEEL_M * sin).max(0.0));
        let stop = k * press.powf(1.5) * (1.0 + chi * HEEL_M * cos * y[5]).max(0.0);
        let fall = G_M_S2 - (tension + MU_K * NORMAL_N) / MASS_KG;
        let fall = if slides { fall } else { 0.0 };
        let torque = (tension * LEVER_M - stop * HEEL_M) * cos - LEAF_N_M_RAD * y[4];
        let pulled = tension / fitting_kg + G_M_S2;
        [y[1], pulled, y[3], fall, y[5], torque / PAWL_KG_M2]
    };
    let advance = |y: &[f64; 6], slides: bool| -> [f64; 6] {
        let at = |k: &[f64; 6], by: f64| std::array::from_fn(|i| y[i] + by * k[i]);
        let k1 = slope(y, slides);
        let k2 = slope(&at(&k1, 0.5 * h), slides);
        let k3 = slope(&at(&k2, 0.5 * h), slides);
        let k4 = slope(&at(&k3, h), slides);
        std::array::from_fn(|i| y[i] + h / 6.0 * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]))
    };
    let load = |c0: f64, c1: f64| (MASS_KG * G_M_S2 - stiffness * 0.5 * (c0 + c1)).abs();
    let mut y = [0.0, 0.0, 0.0, 0.0, theta, 0.0];
    let (mut n, mut broke, mut leave_s) = (0, None, f64::NAN);
    let c0 = closure(&y);
    let (mut saved, mut before) = (y, [c0, load(c0, c0) - MU_S * NORMAL_N]);
    loop {
        let next = advance(&y, broke.is_some());
        n += 1;
        if broke.is_none() && n % SUBSTEPS == 0 {
            let c = closure(&next);
            let over = load(before[0], c) - MU_S * NORMAL_N;
            if over > 0.0 {
                broke = Some((n / SUBSTEPS, before[1] / (before[1] - over)));
                (y, n) = (saved, n - SUBSTEPS);
                continue;
            }
            (saved, before) = (next, [c, over]);
        }
        let t = n as f64 * h;
        let cross = |old: f64, new: f64| t - h * new / (new - old);
        let (old, new) = ((y[4] - RETRACTED).sin(), (next[4] - RETRACTED).sin());
        if leave_s.is_nan() && new <= 0.0 && old > 0.0 {
            leave_s = cross(old, new);
        }
        let (old, new) = (closure(&y), closure(&next));
        if new <= 0.0 {
            let (step, fraction) = broke.unwrap();
            let end = [t, next[4], next[5], next[2], next[3]];
            let slack_s = cross(old, new);
            return Release {
                step,
                fraction,
                leave_s,
                slack_s,
                end,
            };
        }
        y = next;
    }
}

/// What a run at one step length shows: its events; the deepest penetration, the peak load and
/// the plate's heat at the arrest; the largest penetration sampled and the tip's longest move in
/// a step [m]; before the cut, the platform's and the fitting's fastest speeds [m/s], the pawl's
/// largest turn off its rest [rad] and the rail's largest heat [J]; and the fitting's fall past
/// the platform's at the arrest [m].
#[derive(Default)]
struct Run {
    events: [Option<usize>; EVENTS],
    impact: [f64; 3],
    deepest_m: f64,
    travel_m: f64,
    held: [f64; 4],
    gain_m: f64,
}

/// The assembly's edits: the rope and its cut, the leaf, the pawl with its stops.
#[derive(Clone, Copy)]
struct Edits {
    rope: bool,
    leaf: bool,
    pawl: bool,
}

/// The sealed catch and its nodes: the gaps (with the pawl the heel's first, then the plates',
/// the retract stop's and the rope's) with their event slots and contacts, the shapes of the
/// gaps and the pose they start at.
struct Catch {
    system: System,
    bodies: [NodeId; 2],
    fitting: Option<NodeId>,
    nodes: Vec<(usize, NodeId, NodeId)>,
    shapes: Vec<PoseGap>,
    start: Pose,
}

/// The gaps: with the pawl the heel's and the plates', with the rope too the retract stop's;
/// with the rope, last, the rope's from its tie to its face at `face_m`.
fn shapes(edits: Edits, face_m: f64) -> Vec<PoseGap> {
    let (sin, cos) = STRUT.sin_cos();
    let carried = |point, face, normal| PoseGap::new([0.0; 2], [0.0, -1.0], point, face, normal);
    let mut shapes = Vec::new();
    if edits.pawl {
        shapes.push(carried([HEEL_M * sin, HEEL_M * cos], [0.0; 2], [-1.0, 0.0]).unwrap());
        for k in 0..TEETH {
            let plate = carried([-ARM_M * sin, -ARM_M * cos], [0.0, top_m(k)], [0.0, 1.0]).unwrap();
            shapes.push(plate.bounded([-INNER_M, OUTER_M], DEPTH_M).unwrap());
        }
    }
    let (sin, cos) = RETRACTED.sin_cos();
    if edits.pawl && edits.rope {
        shapes.push(carried([HEEL_M * sin, HEEL_M * cos], [0.0; 2], [1.0, 0.0]).unwrap());
    }
    if edits.rope {
        let tie = [LEVER_M * cos, -LEVER_M * sin];
        shapes.push(carried(tie, [0.0, face_m], [0.0, 1.0]).unwrap());
    }
    shapes
}

/// The rope's far end: the fitting under its weight, held by the latch until the cut frees it.
fn hang(b: &mut SystemBuilder, fitting_kg: f64, rope: NodeId, rest: Motion) -> NodeId {
    let p = |b: &SystemBuilder, node, name| b.port(node, name).unwrap();
    let body = b.add_with_tangent("fitting", Mass::new(fitting_kg).unwrap(), rest);
    let weight = Force::new(Waveform::constant(-fitting_kg * G_M_S2).unwrap());
    let weight = b.add_with_tangent("fitting weight", weight, SourceAccount::default());
    let latch = IdealSwitch::in_domain(PortKind::MechanicalForceVelocity).unwrap();
    let latch = b.add_on_line("latch", latch, ()).unwrap();
    let (body, weight) = (body.unwrap(), weight.unwrap());
    let held = [p(b, weight, "p"), p(b, latch, "p"), p(b, rope, "b")];
    b.common_flow(p(b, body, "p"), &held).unwrap();
    let cut = |time: f64, _: f64, out: &mut [f64]| out[0] = f64::from(time >= CUT_S);
    b.scheduled_boundary(p(b, latch, "closed"), Box::new(cut), vec![CUT_S])
        .unwrap();
    body
}

fn assemble(edits: Edits) -> Catch {
    let rotational = PortKind::RotationalTorqueAngularVelocity;
    let rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let [stiffness, fitting_kg, theta, face_m] = hung(edits.leaf);
    // The platform's velocity is up, and what is joined to it receives minus it: the pose's `s`
    // is the fall, along an axis that points down. The turn is the leaf's extension.
    let start = Pose {
        s_m: 0.0,
        theta_rad: if edits.rope { theta } else { RETRACTED },
    };
    let shapes = shapes(edits, face_m);
    let mut b = SystemBuilder::new();
    let body = Mass::new(MASS_KG).unwrap();
    let platform = b.add_with_tangent("platform", body, rest).unwrap();
    let weight = Force::new(Waveform::constant(-MASS_KG * G_M_S2).unwrap());
    let gravity = b
        .add_with_tangent("gravity", weight, SourceAccount::default())
        .unwrap();
    let shoes = Slider::new(NORMAL_N, MU_S, MU_K).unwrap();
    let rail = b.add_on_line("rail", shoes, Slip::stuck()).unwrap();
    let inertia = Mass::in_domain(rotational, PAWL_KG_M2).unwrap();
    let (pawl, leaf) = (edits.pawl, edits.leaf && edits.pawl);
    let pawl = pawl.then(|| b.add_with_tangent("pawl", inertia, rest).unwrap());
    let spring = Spring::in_domain(rotational, LEAF_N_M_RAD, 0.0).unwrap();
    let leaf = leaf.then(|| b.add_with_tangent("leaf", spring, start.theta_rad).unwrap());
    let capacity = HeatCapacity::new(1.0e3).unwrap();
    let sink = b.add_with_tangent("sink", capacity, 300.0).unwrap();
    let mut nodes = Vec::new();
    let first = if edits.pawl { 2 } else { ROPE };
    for (i, shape) in shapes.iter().enumerate() {
        let law = match first + i {
            ROPE => contact("rope", 1.0, [stiffness, 0.0]),
            2 | RETRACT => contact("face", 1.5, STOP),
            _ => contact("face", 1.5, TOOTH),
        };
        let closure = ContactState::at(-shape.gap_m(&start));
        let gap = b.add_with_tangent(&format!("gap {i}"), shape.clone(), start);
        let touch = b.add_on_line(&format!("contact {i}"), law, closure);
        nodes.push((first + i, gap.unwrap(), touch.unwrap()));
    }
    let p = |b: &SystemBuilder, node, name| b.port(node, name).unwrap();
    let still = || Box::new(|_: f64, _: f64, v: &mut [f64]| v[0] = 0.0);
    let mut fall = vec![p(&b, gravity, "p"), p(&b, rail, "a")];
    let mut turn: Vec<_> = leaf.iter().map(|&leaf| p(&b, leaf, "a")).collect();
    let mut heat = vec![p(&b, rail, "heat")];
    let mut fitting = None;
    for &(slot, gap, touch) in &nodes {
        fall.push(p(&b, gap, "s"));
        turn.push(p(&b, gap, "theta"));
        b.bond(p(&b, touch, "a"), p(&b, gap, "gap")).unwrap();
        if slot == ROPE {
            fitting = Some(hang(&mut b, fitting_kg, touch, rest));
        } else {
            heat.push(p(&b, touch, "heat"));
            b.boundary(p(&b, touch, "b"), still()).unwrap();
        }
    }
    b.common_flow(p(&b, platform, "p"), &fall).unwrap();
    if let Some(pawl) = pawl {
        b.common_flow(p(&b, pawl, "p"), &turn).unwrap();
    } else {
        // Without the pawl the rope is tied to the platform where its tie was.
        for port in turn {
            b.boundary(port, still()).unwrap();
        }
    }
    b.boundary(p(&b, rail, "b"), still()).unwrap();
    if let Some(leaf) = leaf {
        b.boundary(p(&b, leaf, "b"), still()).unwrap();
    }
    b.common_effort(p(&b, sink, "heat"), &heat).unwrap();
    let mut system = b.seal().unwrap();
    for &port in &heat {
        system.guess(port, &[300.0]);
    }
    // At one step of the impact at `dt / 2` the default solve leaves the sink's heat rate 6.7e-13
    // W open at 25 iterations and at 100, so each unknown is allowed 1e-12 in its own unit.
    if edits.rope {
        system.set_config(SolveConfig {
            absolute_tolerance: 1.0e-12,
            ..SolveConfig::default()
        });
    }
    Catch {
        system,
        bodies: [platform, rail],
        fitting,
        nodes,
        shapes,
        start,
    }
}

fn catch(plan: &Plan, dt: f64, edits: Edits) -> Run {
    let Catch {
        mut system,
        bodies: [platform, rail],
        fitting,
        nodes,
        shapes,
        start,
    } = assemble(edits);
    let caught = 3 + plan.plate;
    let held = if edits.rope { (CUT_S / dt) as u32 } else { 0 };
    // The rounding each step adds to a closure and to the gap evaluated at the pose, of the size
    // of the gap's terms; and to the kernel's energy, of the size of the energy moved.
    let scale_m = top_m(TEETH - 1).abs() + 2.0 * ARM_M + OUTER_M;
    let energy_j = MASS_KG * G_M_S2 * scale_m + 0.5 * LEAF_N_M_RAD * RETRACTED * RETRACTED;
    // A run with no catch ends when the tip is past the last plate's band.
    let past_m = FIRST_M + PITCH_M * (TEETH - 1) as f64 + DEPTH_M;
    let mut run = Run::default();
    let (mut closures, mut loads, mut arrest) = (Vec::new(), Vec::new(), None);
    // From the step that ends with the rope slack (from the start, with none) the platform's
    // velocity drops by `a dt` a step.
    let mut free = (!edits.rope).then_some((0.0, 0.0));
    let mut last = tip(start.theta_rad);
    with_cx(|cx| {
        for k in 1..=held + (4.0 * plan.contact_s / dt) as u32 {
            system.step(dt, cx).unwrap();
            let t = f64::from(k) * dt;
            let ledger = system.ledger();
            let account = |node: NodeId| &ledger.nodes[node.index()];
            let rounding = |size: f64, count: u32| 16.0 * f64::EPSILON * size * f64::from(count);
            let pose = *system.state::<Pose>(nodes[0].1).unwrap();
            // The fitting's height and velocity, up.
            let lift = fitting.map_or([0.0; 2], |body| {
                let motion = system.state::<Motion>(body).unwrap();
                [motion.position_m, motion.velocity_m_s]
            });
            if edits.pawl {
                let [x, y] = tip(pose.theta_rad);
                let here = [x, y - pose.s_m];
                let moved_m = (here[0] - last[0]).hypot(here[1] - last[1]);
                (run.travel_m, last) = (run.travel_m.max(moved_m), here);
                if x <= INNER_M {
                    run.events[1].get_or_insert(k as usize);
                }
            }
            for (&(slot, gap, touch), shape) in nodes.iter().zip(&shapes) {
                let own = system.state::<Pose>(gap).unwrap();
                let closure_m = system.state::<ContactState>(touch).unwrap().closure_m;
                // The rope's closure is its gap's less the fitting's fall.
                let far_m = if slot == ROPE { lift[0] } else { 0.0 };
                let off = closure_m + shape.gap_m(own) - far_m;
                assert!(off.abs() <= rounding(scale_m, k + 1), "{k} {slot}: {off}");
                if (closure_m > 0.0) == (shape.gap_m(&start) >= 0.0) {
                    run.events[slot].get_or_insert(k as usize);
                }
            }
            if !system.state::<Slip>(rail).unwrap().is_stuck() {
                run.events[0].get_or_insert(k as usize);
            }
            let v = system.state::<Motion>(platform).unwrap().velocity_m_s;
            let rail_j = account(rail).dissipated_j.value();
            if k <= held {
                let turned = (pose.theta_rad - start.theta_rad).abs();
                let moved = [v.abs(), lift[1].abs(), turned, rail_j.abs()];
                run.held = std::array::from_fn(|i| run.held[i].max(moved[i]));
            }
            if free.is_none() && run.events[ROPE].is_some() {
                free = Some((t, v));
            }
            if run.events[caught].is_none() {
                if let Some((from_s, from_m_s)) = free.filter(|&(from_s, _)| t > from_s) {
                    let off = v - from_m_s + plan.accel_m_s2 * (t - from_s);
                    let allowed = rounding(plan.accel_m_s2 * t, k);
                    assert!(off.abs() <= allowed, "step {k}: falls {off}");
                }
            } else {
                let plate = system
                    .state::<ContactState>(nodes[1 + plan.plate].2)
                    .unwrap();
                closures.push(plate.closure_m);
                loads.push(plate.reaction_n);
            }
            if arrest.is_none() && run.events[caught].is_some() && v >= 0.0 {
                arrest = Some(k);
                arrested(ledger, dt, k, energy_j);
                run.impact[2] = account(nodes[1 + plan.plate].2).dissipated_j.value();
                run.gain_m = -lift[0] - pose.s_m;
            } else if arrest.is_none() && run.events[0].is_some() {
                let off = rail_j - MU_K * NORMAL_N * pose.s_m;
                let allowed = rounding(MU_K * NORMAL_N * pose.s_m, k);
                assert!(off.abs() <= allowed, "step {k}: rail heat {off} J off");
            }
            if arrest.is_some_and(|step| k > step + 1) || pose.s_m > past_m {
                break;
            }
        }
    });
    if arrest.is_none() {
        return run;
    }
    run.deepest_m = closures.iter().copied().fold(0.0, f64::max);
    run.impact[0] = vertex(&closures);
    run.impact[1] = vertex(&loads);
    run
}

/// At the arrest the stores of every node and the heat sum to nothing, to the rounding of the
/// energy moved; each node's account is printed.
fn arrested(ledger: &Ledger, dt: f64, k: u32, energy_j: f64) {
    let total: f64 = ledger.nodes.iter().map(|n| n.stored_change_j.value()).sum();
    let allowed = 64.0 * f64::EPSILON * energy_j * f64::from(k);
    assert!(
        total.abs() <= allowed,
        "at the arrest: {total} J of {allowed}"
    );
    println!("dt {dt:e}, arrest at step {k}: stored change, heat [J]; sum {total:.3e}");
    for node in &ledger.nodes {
        let stored = node.stored_change_j.value();
        let heat = node.dissipated_j.value();
        println!("  {}: {stored:.6e}, {heat:.6e}", node.name);
    }
}

/// The kernel's step is second order: its error at `dt / 2` is a third of the change from `dt`.
/// The reference may change between the two by its own error too.
fn second_order(runs: [(&Plan, &Run); 2]) {
    let names = ["deepest penetration [m]", "peak load [N]", "plate heat [J]"];
    for (i, name) in names.iter().enumerate() {
        let [(plan_x, x), (plan_y, y)] = runs.map(|(plan, run)| (plan.impact[i], run.impact[i]));
        let allowed = (x - y).abs() + runs[1].0.quadrature[i] + (plan_x - plan_y).abs();
        assert!(
            (y - plan_y).abs() <= allowed,
            "{name}: {x} and {y} against {plan_y}"
        );
        println!("{name}: {x:.9e} at dt, {y:.9e} at dt/2, {plan_y:.9e} by hand");
    }
}

#[test]
fn a_falling_platform_is_caught_by_the_plate_the_hand_calculation_names() {
    let plan = plan([0.0, RETRACTED, 0.0, 0.0, 0.0]);
    let falls = Edits {
        rope: false,
        leaf: true,
        pawl: true,
    };
    let runs = [DT_S, DT_S / 2.0].map(|dt| (dt, catch(&plan, dt, falls)));
    for (dt, run) in &runs {
        let at = |t: f64| Some((t / dt).ceil() as usize);
        let mut events = [None; EVENTS];
        (events[0], events[1], events[2]) = (Some(1), at(plan.entry_s), at(plan.stop_s));
        events[3 + plan.plate] = at(plan.contact_s);
        assert_eq!(
            run.events, events,
            "dt {dt}: break-away, entry, stop and each plate's contact"
        );
        assert!(
            run.deepest_m <= DEPTH_M / 4.0,
            "dt {dt}: {} m",
            run.deepest_m
        );
        assert!(
            run.travel_m <= DEPTH_M / 100.0,
            "dt {dt}: {} m",
            run.travel_m
        );
    }
    let (at_dt, at_half) = (&runs[0].1, &runs[1].1);
    second_order([(&plan, at_dt), (&plan, at_half)]);
    println!(
        "plate {} closes at step {:?} (t_c / dt = {:.3}); deepest {:.3e} m and longest tip move \
         {:.3e} m in a step, against a band {DEPTH_M:e} m deep",
        plan.plate,
        at_dt.events[3 + plan.plate],
        plan.contact_s / DT_S,
        at_dt.deepest_m,
        at_dt.travel_m,
    );
}

#[test]
fn a_platform_cut_loose_from_its_rope_is_caught_where_the_release_and_the_fall_put_it() {
    let hung = Edits {
        rope: true,
        leaf: true,
        pawl: true,
    };
    let releases = [DT_S, DT_S / 2.0].map(release);
    // Both runs slide from the same instant (the break at `dt / 2` is in the first half of the
    // break at `dt`), so the two references are one up to their RK4 error.
    assert_eq!(2 * releases[0].step - 1, releases[1].step);
    let plans = releases.each_ref().map(|release| plan(release.end));
    let mut runs = Vec::new();
    for ((dt, release), plan) in [DT_S, DT_S / 2.0].into_iter().zip(&releases).zip(&plans) {
        // Each event is at least 0.1 of a step from a step's end.
        let (step, edge) = (release.step, release.fraction.min(1.0 - release.fraction));
        println!("dt {dt:e}: the break-away in step {step}, {edge:.3} of a step from either end");
        assert!(edge >= 0.1, "dt {dt}: the break at {edge}");
        let held = (CUT_S / dt) as usize;
        let at = |t: f64| {
            let (steps, into) = (t / dt, (t / dt).fract());
            println!("dt {dt:e}: an event {steps:.3} steps after the cut");
            assert!(into.min(1.0 - into) >= 0.1, "dt {dt}: {steps} steps");
            Some(held + steps.ceil() as usize)
        };
        let mut events = [None; EVENTS];
        events[0] = Some(held + release.step);
        (events[1], events[2]) = (at(plan.entry_s), at(plan.stop_s));
        events[3 + plan.plate] = at(plan.contact_s);
        (events[RETRACT], events[ROPE]) = (at(release.leave_s), at(release.slack_s));
        let run = catch(plan, dt, hung);
        // Before the cut nothing moves: the platform and the fitting within the solve's 1e-12
        // m/s, and the pawl not at all; the rail makes no heat.
        let moved = run.held[0].max(run.held[1]);
        assert!(moved <= 1.0e-12, "dt {dt}: {:?}", run.held);
        assert_eq!(run.held[2..], [0.0; 2], "dt {dt}");
        assert_eq!(run.events, events, "dt {dt}");
        assert!(run.deepest_m <= DEPTH_M / 4.0 && run.travel_m <= DEPTH_M / 100.0);
        // The fitting falls faster than the platform, and is above it at the arrest.
        let gain_m = run.gain_m;
        println!("dt {dt:e}: the fitting fell {gain_m:.6e} m further than the platform");
        assert!(gain_m > 0.0 && gain_m < ROPE_M[1]);
        runs.push(run);
    }
    second_order([(&plans[0], &runs[0]), (&plans[1], &runs[1])]);
    let no_leaf = Edits {
        leaf: false,
        ..hung
    };
    let no_pawl = Edits {
        pawl: false,
        ..no_leaf
    };
    for edits in [no_leaf, no_pawl] {
        let run = catch(&plans[0], DT_S, edits);
        let none = run.events[1..RETRACT].iter().all(Option::is_none);
        assert!(none && run.events[ROPE].is_some(), "{:?}", run.events);
    }
}
