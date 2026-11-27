//! The Painleve rod: the load line beside a body of several rows. A uniform rod stands at
//! tan(theta) = 2 to a Hertz floor, its lower end sliding towards -x with Coulomb friction `mu`.
//! The rod is a `MultiportInertia` whose rows are its centre's (vx, vy, omega), so
//! `M = diag(m, m, m l^2 / 3)` for the half-length `l`. Its ports are the end's normal flow
//! `vy - l c omega` (the contact), the end's tangential flow `vx + l s omega` (the slider) and
//! the centre's vertical flow `vy` (gravity), with `c` and `s` the cosine and sine of theta. `G`
//! is held at the starting pose, so each run is a short window: the rod turns by under a degree.
//!
//! Friction pushes the end along +x with `mu N`. With `W = G^T M^-1 G`, the end's normal
//! acceleration is `(W_nn + mu W_nt) N + W_ng (-m g)`, and the coefficient of `N` is
//! `(1 + 3 c^2 - 3 mu c s) / m`. At mu = 1 it is 0.4 / m and the sliding end is held by
//! `N = m g / 0.4`; at mu = 1.5 it is -0.2 / m and no sliding contact is consistent (Painleve).
//!
//! A stuck end has no size. Under the default solve the kernel refuses a step soon after the end
//! sticks, `NotConverged` on `slider.a` at a residual of about 1e-16 m/s, the limit stated by
//! `load_line.rs`'s test
//! `an_end_held_at_exactly_nothing_has_no_size_and_is_solved_to_what_the_solve_is_told`. The
//! runs here tell the solve an absolute tolerance of 1e-15 m/s, as that test does, and the
//! mu = 1.5 test pins the default refusal too.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::Obstacle;
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::Component;
use fs_system::elements::{
    Force, HeatCapacity, MultiportInertia, Slider, Slip, SourceAccount, Waveform,
};
use fs_system::{SolveConfig, SystemBuilder, SystemError};

const MASS_KG: f64 = 1.0;
const GRAVITY_M_S2: f64 = 9.81;
const HALF_LENGTH_M: f64 = 0.5;
const DT_S: f64 = 1.0e-5;
/// The end's tangential flow at the start [m/s].
const V0: f64 = -1.0;
/// The absolute tolerance the runs tell the solve, in each unknown's own unit: m/s at the stuck
/// end, whose flow has no size (see the module docs).
const TOLD: f64 = 1.0e-15;
/// Hertz stiffness [N/m^1.5]: `4/3 E* sqrt(R)` with `E* = 0.8 GPa` and `R = 3 mm`.
const HERTZ: f64 = 4.0 / 3.0 * 0.8e9 * 0.054_772_255_750_516_61;
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

/// `G`, row-major: rows vx, vy, omega; columns the end's normal and tangential flows and the
/// centre's vertical flow. The omega row carries the lever arms [m].
fn map() -> [f64; 9] {
    let (c, s, l) = (1.0 / 5f64.sqrt(), 2.0 / 5f64.sqrt(), HALF_LENGTH_M);
    [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, -l * c, l * s, 0.0]
}

/// `(M^-1 G)_ij` for the diagonal `M`.
fn reach(i: usize, j: usize) -> f64 {
    let inverse = [1.0, 1.0, 3.0 / (HALF_LENGTH_M * HALF_LENGTH_M)][i] / MASS_KG;
    inverse * map()[3 * i + j]
}

/// The sliding closed form at `mu`: the coefficient `W_nn + mu W_nt`, the normal force that holds
/// the end on the floor, and the rows' accelerations under it, its friction and the weight.
fn sliding(mu: f64) -> (f64, f64, [f64; 3]) {
    let g = map();
    let w = |a: usize, b: usize| (0..3).map(|i| g[3 * i + a] * reach(i, b)).sum::<f64>();
    let coefficient = w(0, 0) + mu * w(0, 1);
    let weight = MASS_KG * GRAVITY_M_S2;
    let n = w(0, 2) * weight / coefficient;
    let a = std::array::from_fn(|i| (reach(i, 0) + mu * reach(i, 1)) * n - reach(i, 2) * weight);
    (coefficient, n, a)
}

/// A step's record. `mechanical_j` is the energy change since the start of the rod, gravity's
/// account and the contact's store.
#[derive(Clone, Copy)]
struct Step {
    v: [f64; 3],
    tangential: f64,
    closure_m: f64,
    reaction_n: f64,
    stuck: bool,
    mechanical_j: f64,
}

struct Run {
    steps: Vec<Step>,
    refused: Option<SystemError>,
    heat_j: f64,
}

/// `count` steps from the end sliding at `V0` on the floor closed to the Hertz closure of `n0`.
fn simulate(mu: f64, n0: f64, count: usize, told: bool) -> Run {
    let (m, inertia) = (MASS_KG, MASS_KG * HALF_LENGTH_M * HALF_LENGTH_M / 3.0);
    let mass = vec![m, 0.0, 0.0, 0.0, m, 0.0, 0.0, 0.0, inertia];
    let kind = PortKind::MechanicalForceVelocity;
    let body = MultiportInertia::new(kind, mass, map().to_vec(), 3).unwrap();
    // The omega row is not a velocity: the element declares no state dimension.
    assert!(body.state_dimensions().is_none());
    let weight = Force::new(Waveform::constant(-m * GRAVITY_M_S2).unwrap());
    let law = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        HERTZ,
        1.5,
        "floor".into(),
    );
    let floor = ContactInterface::new(law.unwrap(), SOLVE).unwrap();
    let mut b = SystemBuilder::new();
    let rod = b.add_with_tangent("rod", body, vec![V0, 0.0, 0.0]).unwrap();
    let gravity = b
        .add_with_tangent("gravity", weight, SourceAccount::default())
        .unwrap();
    let closure = ContactState::at((n0 / HERTZ).powf(2.0 / 3.0));
    let wall = b
        .add_on_line("wall", floor.publishing_reaction().unwrap(), closure)
        .unwrap();
    let friction = Slider::loaded_by_signal(mu, mu).unwrap();
    let slider = b
        .add_on_line("slider", friction, Slip::sliding(-V0))
        .unwrap();
    let sink = b
        .add_with_tangent("sink", HeatCapacity::new(1.0e3).unwrap(), 300.0)
        .unwrap();
    let p = |b: &SystemBuilder, node, name| b.port(node, name).unwrap();
    b.bond(p(&b, wall, "a"), p(&b, rod, "p0")).unwrap();
    b.bond(p(&b, slider, "a"), p(&b, rod, "p1")).unwrap();
    b.bond(p(&b, gravity, "p"), p(&b, rod, "p2")).unwrap();
    for node in [wall, slider] {
        b.boundary(p(&b, node, "b"), Box::new(|_, _, v: &mut [f64]| v[0] = 0.0))
            .unwrap();
    }
    b.signal(p(&b, wall, "reaction_n"), p(&b, slider, "normal_n"))
        .unwrap();
    let heat = p(&b, slider, "heat");
    b.common_effort(p(&b, sink, "heat"), &[heat]).unwrap();
    let mut system = b.seal().unwrap();
    system.guess(heat, &[300.0]);
    if told {
        system.set_config(SolveConfig {
            absolute_tolerance: TOLD,
            ..SolveConfig::default()
        });
    }
    let g = map();
    let mut steps = Vec::new();
    let refused = with_cx(|cx| {
        for _ in 0..count {
            if let Err(error) = system.step(DT_S, cx) {
                return Some(error);
            }
            let v = system.state::<Vec<f64>>(rod).unwrap();
            let v = [v[0], v[1], v[2]];
            let contact = system.state::<ContactState>(wall).unwrap();
            let ledger = system.ledger();
            let stored =
                |node: fs_system::NodeId| ledger.nodes[node.index()].stored_change_j.value();
            steps.push(Step {
                v,
                tangential: (0..3).map(|i| g[3 * i + 1] * v[i]).sum(),
                closure_m: contact.closure_m,
                reaction_n: contact.reaction_n,
                stuck: system.state::<Slip>(slider).unwrap().is_stuck(),
                mechanical_j: stored(rod) + stored(gravity) + stored(wall),
            });
        }
        None
    });
    let heat_j = system.ledger().nodes[sink.index()].stored_change_j.value();
    Run {
        steps,
        refused,
        heat_j,
    }
}

/// From the step the end sticks to the window's end: it stays stuck and on the floor, its
/// tangential flow is nothing to the solve's tolerance, and over the whole run the mechanical
/// energy never rises above the rounding of its accounts and the rod turns by under a degree.
/// Returns the index of the step that sticks.
fn holds_after_the_stick(run: &Run) -> usize {
    let stuck = run
        .steps
        .iter()
        .position(|step| step.stuck)
        .expect("the end sticks");
    let flow = SolveConfig::default().relative_tolerance * V0.abs() + TOLD;
    for (k, step) in (1..).zip(&run.steps).skip(stuck) {
        assert!(step.stuck && step.reaction_n > 0.0, "step {k} lets go");
        assert!(
            step.tangential.abs() <= flow,
            "step {k}: {:e} m/s",
            step.tangential
        );
    }
    let rounding = 64.0 * f64::EPSILON * 0.5 * MASS_KG * V0 * V0;
    let mut before = 0.0;
    for (k, step) in (1..).zip(&run.steps) {
        let rise = step.mechanical_j - before;
        assert!(
            rise <= rounding,
            "step {k}: the mechanical energy rises by {rise:e} J"
        );
        before = step.mechanical_j;
    }
    let turned: f64 = run.steps.iter().map(|step| step.v[2].abs() * DT_S).sum();
    assert!(turned < 1f64.to_radians(), "{turned} rad");
    stuck
}

#[test]
fn at_mu_one_the_end_slides_and_sticks_as_the_closed_form_says() {
    let mu = 1.0;
    let (coefficient, n, a) = sliding(mu);
    assert!(coefficient > 0.0, "{coefficient}");
    // The end's tangential acceleration; it stops in the step in which it would pass zero.
    let along: f64 = (0..3).map(|i| map()[3 * i + 1] * a[i]).sum();
    let stop = (V0.abs() / along / DT_S).ceil() as usize;
    let run = simulate(mu, n, stop + 100, true);
    assert!(run.refused.is_none(), "{:?}", run.refused);
    let stuck = holds_after_the_stick(&run);
    assert_eq!(stuck + 1, stop);
    // While it slides, each step changes the rows by the closed form to the step's own error: the
    // contact's force tolerance and the kernel's relative tolerance on the block's largest force,
    // carried by M^-1 G, and the rounding of a difference of two velocities.
    let largest = (n * mu.max(1.0)).max(MASS_KG * GRAVITY_M_S2);
    let kernel = SolveConfig::default().relative_tolerance * largest;
    let contact = SOLVE.force_absolute_tolerance_n + SOLVE.force_relative_tolerance * n;
    let (normal, tangential) = (contact + kernel, mu * contact + kernel);
    let mut start = [V0, 0.0, 0.0];
    let mut worst = 0.0_f64;
    for step in &run.steps[..stuck] {
        for i in 0..3 {
            let size = step.v[i].abs().max(start[i].abs());
            let allowed = reach(i, 0).abs() * normal
                + reach(i, 1).abs() * tangential
                + 4.0 * f64::EPSILON * size / DT_S;
            worst = worst.max(((step.v[i] - start[i]) / DT_S - a[i]).abs() / allowed);
        }
        start = step.v;
    }
    println!(
        "mu = 1: N = {n} N, a = {a:?}; sticks at step {stop} as the closed form says; the slide's \
         worst step is {worst:.3} of its own error; heat {:.9e} J",
        run.heat_j
    );
    assert!(worst <= 1.0, "{worst}");
}

#[test]
fn at_mu_one_and_a_half_where_no_slide_is_consistent_the_end_jams_and_holds() {
    let mu = 1.5;
    let (coefficient, ..) = sliding(mu);
    assert!(coefficient < 0.0, "{coefficient}");
    // From the weight on the floor the normal force grows on its own while the end slides, until
    // the end stops: the jam, on a compliant floor.
    let weight = MASS_KG * GRAVITY_M_S2;
    let run = simulate(mu, weight, 600, true);
    assert!(run.refused.is_none(), "{:?}", run.refused);
    let stuck = holds_after_the_stick(&run);
    assert!(stuck < 500, "sticks at step {}", stuck + 1);
    let slide = &run.steps[..=stuck];
    assert!(
        slide
            .windows(2)
            .all(|pair| pair[1].reaction_n > pair[0].reaction_n)
    );
    assert!(slide[stuck].reaction_n > 10.0 * weight);
    let deepest = run
        .steps
        .iter()
        .map(|step| step.closure_m)
        .fold(0.0, f64::max);
    assert!(deepest <= SOLVE.maximum_penetration_m, "{deepest} m");
    println!(
        "mu = 1.5: the reaction grows from {:.6} N to {:.6} N and the end sticks at step {}; \
         the deepest closure is {deepest:.6e} m; heat {:.9e} J",
        slide[0].reaction_n,
        slide[stuck].reaction_n,
        stuck + 1,
        run.heat_j
    );
    // Under the default solve the stuck end is refused by name, not answered wrongly.
    let default = simulate(mu, weight, 600, false);
    assert!(default.steps.len() > stuck, "refused before the stick");
    match default.refused {
        Some(SystemError::NotConverged {
            worst, residual, ..
        }) => {
            assert_eq!(worst, "slider.a");
            assert!(residual.abs() < 1.0e-15, "{residual}");
            println!("default solve: refused at step {}", default.steps.len() + 1);
        }
        other => panic!("{other:?}"),
    }
}
