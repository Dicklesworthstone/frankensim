//! The catch (`frankensim-rc-1109-root-3gc12.13.1`, O2b): a platform that falls on a rail is
//! caught by a pawl on a ladder of tooth plates, from generic elements, against hand
//! calculations.
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
//! No claim: the pawl's centre of mass is at its pivot (no coupled mass matrix), the tip has no
//! friction, the plates have no flanks or backs as faces, one pawl, and the rope that holds the
//! pawl retracted before the cut is not modelled (O3).

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::Obstacle;
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_system::elements::{
    Force, HeatCapacity, Mass, Motion, Pose, PoseGap, Slider, Slip, SourceAccount, Spring, Waveform,
};
use fs_system::{NodeId, System, SystemBuilder};

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
/// Hertz stiffness [N/m^1.5] and Hunt-Crossley loss [s/m] of a plate and of the stop.
const TOOTH: [f64; 2] = [1.0e10, 2.0];
const STOP: [f64; 2] = [1.0e11, 50.0];
const DT_S: f64 = 1.0e-6;
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

fn contact(name: &str, [stiffness, loss]: [f64; 2]) -> ContactInterface {
    let law = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        stiffness,
        1.5,
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

fn plan() -> Plan {
    let accel_m_s2 = G_M_S2 - MU_K * NORMAL_N / MASS_KG;
    let omega = (LEAF_N_M_RAD / PAWL_KG_M2).sqrt();
    let entry = STRUT + (INNER_M / ARM_M).asin();
    let entry_s = (entry / RETRACTED).acos() / omega;
    let height = -0.5 * accel_m_s2 * entry_s * entry_s + tip(entry)[1];
    let plate = (0..TEETH).find(|&k| height >= top_m(k)).unwrap();
    // The plates above it have the tip behind them, under their band.
    assert!((0..plate).all(|k| height < top_m(k) - DEPTH_M));
    let contact_s = (2.0 * (-ARM_M - top_m(plate)) / accel_m_s2).sqrt();
    let v0 = accel_m_s2 * contact_s;
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
        stop_s: (STRUT / RETRACTED).acos() / omega,
        plate,
        contact_s,
        impact: fine,
        quadrature: [0, 1, 2].map(|i| (coarse[i] - fine[i]).abs()),
    }
}

/// What a run at one step length shows: the steps of the entry, the stop and each plate's
/// contact; the deepest penetration, the peak load and the plate's heat at the arrest; the
/// largest penetration sampled and the tip's longest move in a step [m].
#[derive(Default)]
struct Run {
    events: [Option<usize>; 2 + TEETH],
    impact: [f64; 3],
    deepest_m: f64,
    travel_m: f64,
}

/// The sealed catch and its nodes: the gaps (the heel's first, then the plates') with their
/// contacts, and the shapes of the gaps.
struct Catch {
    system: System,
    bodies: [NodeId; 6],
    nodes: Vec<(NodeId, NodeId)>,
    shapes: Vec<PoseGap>,
}

fn assemble() -> Catch {
    let rotational = PortKind::RotationalTorqueAngularVelocity;
    let rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    // The platform's velocity is up, and what is joined to it receives minus it: the pose's `s`
    // is the fall, along an axis that points down. The turn is the leaf's extension.
    let start = Pose {
        s_m: 0.0,
        theta_rad: RETRACTED,
    };
    let (sin, cos) = STRUT.sin_cos();
    let carried = |point, face, normal| PoseGap::new([0.0; 2], [0.0, -1.0], point, face, normal);
    let mut shapes = vec![carried([HEEL_M * sin, HEEL_M * cos], [0.0; 2], [-1.0, 0.0]).unwrap()];
    for k in 0..TEETH {
        let plate = carried([-ARM_M * sin, -ARM_M * cos], [0.0, top_m(k)], [0.0, 1.0]).unwrap();
        shapes.push(plate.bounded([-INNER_M, OUTER_M], DEPTH_M).unwrap());
    }
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
    let pawl = b.add_with_tangent("pawl", inertia, rest).unwrap();
    let spring = Spring::in_domain(rotational, LEAF_N_M_RAD, 0.0).unwrap();
    let leaf = b.add_with_tangent("leaf", spring, RETRACTED).unwrap();
    let capacity = HeatCapacity::new(1.0e3).unwrap();
    let sink = b.add_with_tangent("sink", capacity, 300.0).unwrap();
    let mut nodes: Vec<(NodeId, NodeId)> = Vec::new();
    for (i, shape) in shapes.iter().enumerate() {
        let law = if i == 0 { STOP } else { TOOTH };
        let closure = ContactState::at(-shape.gap_m(&start));
        let gap = b.add_with_tangent(&format!("gap {i}"), shape.clone(), start);
        let touch = b.add_on_line(&format!("contact {i}"), contact("face", law), closure);
        nodes.push((gap.unwrap(), touch.unwrap()));
    }
    let p = |b: &SystemBuilder, node, name| b.port(node, name).unwrap();
    let still = || Box::new(|_: f64, _: f64, v: &mut [f64]| v[0] = 0.0);
    let mut fall = vec![p(&b, gravity, "p"), p(&b, rail, "a")];
    let (mut turn, mut heat) = (vec![p(&b, leaf, "a")], vec![p(&b, rail, "heat")]);
    for &(gap, touch) in &nodes {
        fall.push(p(&b, gap, "s"));
        turn.push(p(&b, gap, "theta"));
        heat.push(p(&b, touch, "heat"));
        b.bond(p(&b, touch, "a"), p(&b, gap, "gap")).unwrap();
        b.boundary(p(&b, touch, "b"), still()).unwrap();
    }
    b.common_flow(p(&b, platform, "p"), &fall).unwrap();
    b.common_flow(p(&b, pawl, "p"), &turn).unwrap();
    b.boundary(p(&b, rail, "b"), still()).unwrap();
    b.boundary(p(&b, leaf, "b"), still()).unwrap();
    b.common_effort(p(&b, sink, "heat"), &heat).unwrap();
    let mut system = b.seal().unwrap();
    for &port in &heat {
        system.guess(port, &[300.0]);
    }
    Catch {
        system,
        bodies: [platform, gravity, rail, pawl, leaf, sink],
        nodes,
        shapes,
    }
}

fn catch(plan: &Plan, dt: f64) -> Run {
    let Catch {
        mut system,
        bodies: [platform, gravity, rail, pawl, leaf, sink],
        nodes,
        shapes,
    } = assemble();
    let caught = 2 + plan.plate;
    // The rounding each step adds to a closure and to the gap evaluated at the pose, of the size
    // of the gap's terms; and to the kernel's energy, of the size of the energy moved.
    let scale_m = top_m(TEETH - 1).abs() + 2.0 * ARM_M + OUTER_M;
    let energy_j = MASS_KG * G_M_S2 * scale_m + 0.5 * LEAF_N_M_RAD * RETRACTED * RETRACTED;
    let mut run = Run::default();
    let (mut closures, mut loads, mut arrest) = (Vec::new(), Vec::new(), None);
    let mut last = tip(RETRACTED);
    with_cx(|cx| {
        for k in 1_u32.. {
            system.step(dt, cx).unwrap();
            let t = f64::from(k) * dt;
            let ledger = system.ledger();
            let account = |node: NodeId| &ledger.nodes[node.index()];
            let rounding = |size: f64, count: u32| 16.0 * f64::EPSILON * size * f64::from(count);
            let pose = *system.state::<Pose>(nodes[0].0).unwrap();
            let [x, y] = tip(pose.theta_rad);
            let here = [x, y - pose.s_m];
            run.travel_m = run
                .travel_m
                .max((here[0] - last[0]).hypot(here[1] - last[1]));
            last = here;
            if x <= INNER_M {
                run.events[0].get_or_insert(k as usize);
            }
            for (i, (&(gap, touch), shape)) in nodes.iter().zip(&shapes).enumerate() {
                let own = system.state::<Pose>(gap).unwrap();
                let closure_m = system.state::<ContactState>(touch).unwrap().closure_m;
                let off = closure_m + shape.gap_m(own);
                assert!(
                    off.abs() <= rounding(scale_m, k + 1),
                    "step {k}, gap {i}: {off}"
                );
                if closure_m > 0.0 {
                    run.events[1 + i].get_or_insert(k as usize);
                }
            }
            let v = system.state::<Motion>(platform).unwrap().velocity_m_s;
            if run.events[caught].is_none() {
                let off = v + plan.accel_m_s2 * t;
                assert!(
                    off.abs() <= rounding(plan.accel_m_s2 * t, k),
                    "step {k}: falls {off}"
                );
            } else {
                let plate = system
                    .state::<ContactState>(nodes[1 + plan.plate].1)
                    .unwrap();
                closures.push(plate.closure_m);
                loads.push(plate.reaction_n);
            }
            if arrest.is_none() && run.events[caught].is_some() && v >= 0.0 {
                arrest = Some(k);
                let total: f64 = ledger.nodes.iter().map(|n| n.stored_change_j.value()).sum();
                let allowed = 4.0 * rounding(energy_j, k);
                assert!(
                    total.abs() <= allowed,
                    "at the arrest: {total} J of {allowed}"
                );
                let plate = account(nodes[1 + plan.plate].1).dissipated_j.value();
                run.impact[2] = plate;
                let stored = |node: NodeId| account(node).stored_change_j.value();
                println!(
                    "dt {dt:e}, arrest at step {k}: stored changes [J] platform {:.6e}, gravity \
                     {:.6e}, pawl {:.3e}, leaf {:.6e}, stop {:.3e}, plate {:.6e}, sink {:.6e} \
                     (rail {:.6e}, stop {:.6e}, plate {plate:.6e} of heat); sum {total:.3e}",
                    stored(platform),
                    stored(gravity),
                    stored(pawl),
                    stored(leaf),
                    stored(nodes[0].1),
                    stored(nodes[1 + plan.plate].1),
                    stored(sink),
                    account(rail).dissipated_j.value(),
                    account(nodes[0].1).dissipated_j.value(),
                );
            } else if arrest.is_none() {
                let off = account(rail).dissipated_j.value() - MU_K * NORMAL_N * pose.s_m;
                let allowed = rounding(MU_K * NORMAL_N * pose.s_m, k);
                assert!(
                    off.abs() <= allowed,
                    "step {k}: the rail's heat is {off} J off"
                );
            }
            if arrest.is_some_and(|step| k > step + 1) {
                break;
            }
        }
    });
    run.deepest_m = closures.iter().copied().fold(0.0, f64::max);
    run.impact[0] = vertex(&closures);
    run.impact[1] = vertex(&loads);
    run
}

#[test]
fn a_falling_platform_is_caught_by_the_plate_the_hand_calculation_names() {
    let plan = plan();
    let runs = [DT_S, DT_S / 2.0].map(|dt| (dt, catch(&plan, dt)));
    for (dt, run) in &runs {
        let at = |t: f64| Some((t / dt).ceil() as usize);
        let mut events = [None; 2 + TEETH];
        (events[0], events[1]) = (at(plan.entry_s), at(plan.stop_s));
        events[2 + plan.plate] = at(plan.contact_s);
        assert_eq!(
            run.events, events,
            "dt {dt}: entry, stop and each plate's contact"
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
    // The kernel's step is second order: its error at dt/2 is a third of the change from dt.
    let (at_dt, at_half) = (&runs[0].1, &runs[1].1);
    let names = ["deepest penetration [m]", "peak load [N]", "plate heat [J]"];
    for (i, name) in names.iter().enumerate() {
        let (x, y, reference) = (at_dt.impact[i], at_half.impact[i], plan.impact[i]);
        let allowed = (x - y).abs() + plan.quadrature[i];
        assert!(
            (y - reference).abs() <= allowed,
            "{name}: {x} and {y} against {reference}"
        );
        println!("{name}: {x:.9e} at dt, {y:.9e} at dt/2, {reference:.9e} by hand");
    }
    println!(
        "plate {} closes at step {:?} (t_c / dt = {:.3}); deepest {:.3e} m and longest tip move \
         {:.3e} m in a step, against a band {DEPTH_M:e} m deep",
        plan.plate,
        at_dt.events[2 + plan.plate],
        plan.contact_s / DT_S,
        at_dt.deepest_m,
        at_dt.travel_m,
    );
}
