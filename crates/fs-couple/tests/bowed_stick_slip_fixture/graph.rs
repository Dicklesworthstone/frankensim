//! The bowed string as a composition on the kernel (bead rc-1109 .12.5).
//!
//! **Graph.** The fixture's string is fs-nlmodal's `ModalStructure`, publishing
//! the end of its flow at the bow. The contact is fs-tribo's
//! `FrictionInterface` with the fixture's Stribeck law as `HystereticFriction`,
//! reading that end on the kernel's load line. The bow is a boundary at the bow
//! speed, and each node's heat goes to a sink at the air's temperature. Nothing
//! here is hand-coupled: the kernel steps the three parts, and the interface
//! cuts a step where the slip stops or the pair lets go.
//!
//! **Against the module.** The graph is judged by the module's own verdicts and
//! figures (Helmholtz or not, the sliding fraction, the period), not to
//! rounding. The two differ in rule, not only in rounding. The graph asks the
//! law at a sliding step's mean slip and locates its events inside a step; the
//! module asks at the end slip and takes whole steps. Held at no end slip, a
//! step's mean slip is not zero beside a string, and the work of that defect is
//! booked as numerical loss. It shrinks with the step on both sides (bead
//! comment, 2026-11-14).

use core::num::NonZeroUsize;
use std::time::Instant;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_nlmodal::component::{ModalAcousticTimeBudget, ModalAttachment, ModalStructure};
use fs_ports::{PortKind, StableId};
use fs_system::{RouteTemperature, SystemBuilder, heat_sink};
use fs_tribo::stribeck_friction::{HystereticFriction, HystereticFrictionState};
use fs_tribo::tangential::component::{FrictionInterface, FrictionState, InterfaceLoad};
use fs_tribo::{FrictionLaw, TangentialSolveConfig, TriboError};

use super::*;

/// The air's temperature, of both heat sinks and of the interface [K].
const AIR_K: f64 = 293.15;

type LawAt = fn(f64) -> Result<HystereticFriction, TriboError>;
type Rubbed = FrictionState<HystereticFrictionState>;

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

/// The fixture's friction, as the module lowers its Stribeck declaration.
fn dry_law() -> FrictionLaw {
    FrictionLaw::Stribeck {
        static_mu: MU_S,
        kinetic_mu: MU_D,
        characteristic_speed: DECAY_M_S,
        viscous_per_speed: 0.0,
    }
}

fn law_at(_: f64) -> Result<HystereticFriction, TriboError> {
    HystereticFriction::new(dry_law())
}

/// Unit-mass mode shapes at the bow station, as the module forms them.
fn shapes_at(modes: usize, beta: f64) -> Vec<f64> {
    let norm = (DENSITY_KG_M * length_m() * 0.5).sqrt();
    (1..=modes)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let kk = k as f64;
            (kk * core::f64::consts::PI * beta).sin() / norm
        })
        .collect()
}

/// The fixture's string as the modal structure, one attachment at the bow.
fn string(rate: u32, modes: usize, beta: f64) -> ModalStructure {
    let card = card(rate, modes);
    let pairs: Vec<(f64, f64)> = (1..=modes)
        .map(|k| (card.mode_omega_rad_s(k), card.zetas[k - 1]))
        .collect();
    ModalStructure::from_modes(
        rate,
        &pairs,
        ModalAcousticTimeBudget::audible_reference(),
        vec![ModalAttachment {
            id: StableId::new("bow").unwrap(),
            kind: PortKind::MechanicalForceVelocity,
            width: NonZeroUsize::MIN,
            weights: shapes_at(modes, beta),
        }],
    )
    .unwrap()
}

/// One step the graph took.
struct Taken {
    /// Its length and the time it ended at [s].
    dt_s: f64,
    end_s: f64,
    started_at_rest: bool,
    ended_at_rest: bool,
    /// The work the two ends gave up over it, force times mean slip times dt [J].
    work_j: f64,
}

/// The graph's run: its steps, the refusal that stopped it if any, and its wall.
struct Played {
    steps: Vec<Taken>,
    refusal: Option<String>,
    wall_s: f64,
}

/// The bowed graph from the fixture's Helmholtz start, stepped at `1/rate` until
/// `until_s`.
fn play(rate: u32, modes: usize, beta: f64, force_n: f64, until_s: f64) -> Played {
    let h = 1.0 / f64::from(rate);
    let structure = string(rate, modes, beta).publishing_ends().unwrap();
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let string_node = builder
        .add_with_tangent("string", structure, helmholtz_start(modes, beta))
        .unwrap();
    let load = InterfaceLoad {
        normal_n: force_n,
        nominal_pressure_pa: None,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 1.0e6,
            force_absolute_tolerance_n: 1.0e-12,
            force_relative_tolerance: 1.0e-12,
        },
        impulse_tolerance_n_s: 1.0e-10,
        work_tolerance_j: 1.0e-12,
    };
    // The bow's end is a boundary, so only the string's end is read.
    let interface = FrictionInterface::new(law_at as LawAt, load)
        .unwrap()
        .reading_ends([true, false])
        .unwrap();
    let stuck = interface.start(AIR_K, 0.0).unwrap();
    let rubbing = builder.add_on_line("rubbing", interface, stuck).unwrap();
    let (a, b) = (port(&builder, rubbing, "a"), port(&builder, rubbing, "b"));
    builder
        .common_flow(port(&builder, string_node, "bow"), &[a])
        .unwrap();
    builder
        .boundary(b, Box::new(|_, _, values| values[0] = -BOW_SPEED_M_S))
        .unwrap();
    for (node, sink) in [(string_node, "air"), (rubbing, "contact")] {
        builder
            .heat_route(
                port(&builder, node, "heat"),
                &[],
                Some(heat_sink(sink, AIR_K).unwrap()),
                RouteTemperature::Fixed(AIR_K),
            )
            .unwrap();
    }
    let mut system = builder.seal().unwrap();
    let mut steps = Vec::new();
    let mut refusal = None;
    let started = Instant::now();
    with_cx(|cx| {
        while system.time_s() < until_s {
            let started_at_rest = system.state::<Rubbed>(rubbing).unwrap().at_rest();
            let report = match system.step(h, cx) {
                Ok(report) => report,
                Err(error) => {
                    refusal = Some(format!(
                        "refused at step {} from {} s: {error}",
                        steps.len() + 1,
                        system.time_s()
                    ));
                    break;
                }
            };
            let mean_m_s = system.last_input(a)[0] - system.last_input(b)[0];
            let force_n = system.last_output(a)[0];
            steps.push(Taken {
                dt_s: report.dt_s,
                end_s: report.time_s,
                started_at_rest,
                ended_at_rest: system.state::<Rubbed>(rubbing).unwrap().at_rest(),
                work_j: report.dt_s * force_n * mean_m_s,
            });
        }
    });
    Played {
        steps,
        refusal,
        wall_s: started.elapsed().as_secs_f64(),
    }
}

/// What a run shows by the module's rules, over the last `TAIL_PERIODS` periods:
/// classify's Helmholtz verdict, the sliding fraction, the slip onsets and the
/// intervals between them [s].
struct Figures {
    helmholtz: bool,
    fraction: f64,
    onsets: usize,
    intervals_s: Vec<f64>,
}

impl Figures {
    fn period_s(&self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let n = self.intervals_s.len() as f64;
        self.intervals_s.iter().sum::<f64>() / n
    }

    /// The shortest and the longest interval [s].
    fn span_s(&self) -> (f64, f64) {
        let least = self
            .intervals_s
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let most = self
            .intervals_s
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        (least, most)
    }
}

/// The module's figures, from its receipts. Its onsets lie on the step grid.
#[allow(clippy::cast_precision_loss)]
fn module_figures(trace: &[BowContactReceipt], h: f64) -> Figures {
    let (helmholtz, fraction) = classify(trace, h);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let window = (TAIL_PERIODS as f64 / F0_HZ / h).round() as usize;
    let found = onsets(&trace[trace.len() - window..]);
    Figures {
        helmholtz,
        fraction,
        onsets: found.len(),
        intervals_s: found
            .windows(2)
            .map(|pair| (pair[1] - pair[0]) as f64 * h)
            .collect(),
    }
}

/// The graph's figures by the same rules, from the tail starting at
/// `tail_from_s`.
/// - An onset is a step that starts at rest and ends sliding, at that step's
///   start.
/// - The sliding time is that of the steps that END sliding, as a receipt's
///   branch is the one after its step.
fn graph_figures(steps: &[Taken], tail_from_s: f64) -> Figures {
    let releases: Vec<f64> = steps
        .iter()
        .filter(|step| step.started_at_rest && !step.ended_at_rest)
        .map(|step| step.end_s - step.dt_s)
        .filter(|&at_s| at_s >= tail_from_s)
        .collect();
    let intervals_s: Vec<f64> = releases.windows(2).map(|pair| pair[1] - pair[0]).collect();
    let helmholtz = releases.len() + 1 >= TAIL_PERIODS
        && releases.len() <= TAIL_PERIODS
        && intervals_s
            .iter()
            .all(|interval| (interval * F0_HZ - 1.0).abs() <= 0.05);
    let end_s = steps.last().map_or(tail_from_s, |step| step.end_s);
    let sliding_s: f64 = steps
        .iter()
        .filter(|step| !step.ended_at_rest)
        .map(|step| (step.end_s - (step.end_s - step.dt_s).max(tail_from_s)).max(0.0))
        .sum();
    Figures {
        helmholtz,
        fraction: sliding_s / (end_s - tail_from_s),
        onsets: releases.len(),
        intervals_s,
    }
}

/// The held steps' work against the force as a share of the heat of the rest,
/// in % (the graph's booking: held is started and ended at rest).
fn held_against_percent(steps: &[Taken]) -> f64 {
    let (mut against_j, mut slid_j) = (0.0_f64, 0.0_f64);
    for step in steps {
        if step.started_at_rest && step.ended_at_rest {
            against_j += step.work_j.min(0.0);
        } else {
            slid_j += step.work_j;
        }
    }
    100.0 * against_j / slid_j
}

#[test]
#[allow(clippy::cast_precision_loss)]
fn g2_the_bowed_graph_keeps_the_module_s_helmholtz_motion() {
    // The coarsest rung, at the geometric mean of Schelleng's bounds for beta
    // 1/10, from the Helmholtz start, as g1 runs the module.
    let (rate, modes) = RUNGS[0];
    let beta = 1.0 / 10.0;
    let force_n = det::sqrt(schelleng_fmin(beta) * schelleng_fmax(beta));
    let h = 1.0 / f64::from(rate);
    let until_s = RUN_PERIODS / F0_HZ;
    let tail_from_s = until_s - TAIL_PERIODS as f64 / F0_HZ;

    let graph = play(rate, modes, beta, force_n, until_s);
    let started = Instant::now();
    let mut state = voice(rate, modes, 1, beta, force_n);
    let trace = receipts(&mut state, rate, RUN_PERIODS);
    let module_wall_s = started.elapsed().as_secs_f64();

    let end_s = graph.steps.last().map_or(0.0, |step| step.end_s);
    let ours = graph_figures(&graph.steps, tail_from_s);
    let theirs = module_figures(&trace, h);
    let against_percent = held_against_percent(&graph.steps);
    let fraction_tolerance = 0.05 * beta + 2.0 * F0_HZ / f64::from(rate);
    println!(
        "{{\"test\":\"bowed_graph_against_module\",\"rung_hz\":{rate},\"modes\":{modes},\
         \"beta\":{beta},\"force_n\":{force_n},\"refusal\":\"{}\",\"end_s\":{end_s},\
         \"steps\":[{}, {}],\"us_a_step\":[{}, {}],\"helmholtz\":[{}, {}],\
         \"tail_onsets\":[{}, {}],\"sliding_fraction\":[{}, {}],\"fraction_tolerance\":{fraction_tolerance},\
         \"tail_period_s\":[{}, {}],\"tail_intervals_s\":[{:?}, {:?}],\
         \"period_off_in_steps\":{},\"held_against_percent_of_slid\":{against_percent}}}",
        graph.refusal.as_deref().unwrap_or("none"),
        graph.steps.len(),
        trace.len(),
        1.0e6 * graph.wall_s / graph.steps.len() as f64,
        1.0e6 * module_wall_s / trace.len() as f64,
        ours.helmholtz,
        theirs.helmholtz,
        ours.onsets,
        theirs.onsets,
        ours.fraction,
        theirs.fraction,
        ours.period_s(),
        theirs.period_s(),
        ours.span_s(),
        theirs.span_s(),
        (ours.period_s() - theirs.period_s()) / h
    );

    // 1. The graph plays the whole run.
    assert!(
        graph.refusal.is_none() && end_s >= until_s,
        "the graph stopped at {end_s} s of {until_s} s: {:?}",
        graph.refusal
    );
    // 2. Helmholtz on both sides, by classify's rule.
    assert!(theirs.helmholtz, "the module is not in Helmholtz motion");
    assert!(
        ours.helmholtz,
        "the graph is not in Helmholtz motion: {} onsets in the last {TAIL_PERIODS} periods, \
         intervals {:?} s against 1/f0 = {} s",
        ours.onsets,
        ours.intervals_s,
        1.0 / F0_HZ
    );
    // 3. The sliding fraction within g1's own tolerance of the module's.
    assert!(
        (ours.fraction - theirs.fraction).abs() <= fraction_tolerance,
        "sliding fraction {} against the module's {} (tolerance {fraction_tolerance})",
        ours.fraction,
        theirs.fraction
    );
    // 4. The period within one step of the module's, which resolves no finer.
    assert!(
        (ours.period_s() - theirs.period_s()).abs() <= h,
        "period {} s against the module's {} s (a step is {h} s)",
        ours.period_s(),
        theirs.period_s()
    );
    // 5. A regression guard, set above a measured figure and not derived: the
    // held steps' work against the force was 0.6217 % of the slid heat here
    // (the module stepped by hand on its own end line: 0.4123 %; both shrink
    // with the step, bead comment 2026-11-14).
    assert!(
        against_percent >= -1.0,
        "the held steps' work against the force is {against_percent} % of the slid heat"
    );
}
