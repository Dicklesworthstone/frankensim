//! The memory pair declared as a draft and lowered by the owner constructors
//! (`fs_system::lowering::owners`), against `memory_graph` above with no
//! actuator. It is 1b's first pair of 2 kg boxes closing at 2 m/s, its bushing
//! on a `relaxation-bushing` card: Kelvin viscosity 0.4 N s/m and one branch
//! of stiffness 1.69 N/m (factor 1.3) and time 0.2 s, all on x. The interface
//! declares the branch's elastic displacement, 0.6 / 1.3 m, and two loss
//! routes send 0.7 and 0.2 of its heat to the cells. Fixture values; not a
//! seed-dataset authority. With the module's Arrhenius shift on the card and
//! cell 0 as its `TemperatureSource`, against `memory_graph` with that shift
//! on cell 0's published temperature. With bushing0 holding 0.7 N or -0.7 N
//! on x (a constant source on its `drive_x`), against `memory_graph` with
//! that held force. With bushing0's damping on `module_curve(0.4)` instead,
//! against `memory_graph` with the test's `CurveLaw` on that claim. With an
//! equilibrium spring beside the branch, against `memory_graph` with that spring.

use super::*;
use crate::graph::lowering::{
    CurveLaw, admit, curve, curved, fixture_cards, fixture_laws, interface_card, lower,
    lowered_mass, module_curve, observe, pairing, sourced, temperature_bits, under,
};
use crate::graph::with_cx;
use fs_exec::CancelGate;
use fs_matdb::{InterpolationPolicy, PropertyValue};
use fs_material::visco::temperature::TemperatureShiftModel;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::ParameterKey;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftAssemblyInput, DraftInterfaceKind, SourceDecl, SourceTarget, TimeHistory,
    TimeScheme,
};
use fs_system::lowering::Observed;
use fs_system::lowering::owners::{
    RELAXATION_BUSHING_FAMILY, VISCOUS_BUSHING_FAMILY, register_owners,
};
use fs_system::{CardStore, Ledger, LoweringRegistry, ParameterSource};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const DAMPING: Dims = Dims([0, 1, -1, 0, 0, 0]);
const TORSIONAL_DAMPING: Dims = Dims([2, 1, -1, 0, 0, 0]);
const STIFFNESS: Dims = Dims([0, 1, -2, 0, 0, 0]);
const TORSIONAL_STIFFNESS: Dims = Dims([2, 1, -2, 0, 0, 0]);
const KELVIN: Dims = Dims([0, 0, 0, 1, 0, 0]);
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
/// The branch's elastic displacement [m]: 1.3 times it is the fixture's
/// 0.6 sqrt(J), to the bit.
const DISPLACEMENT_M: f64 = 0.461_538_461_538_461_5;

type State = Vec<(ParameterKey, QtyAny)>;

/// The branch's declared displacement: `x` on x, 0 on the other axes.
fn displacement(x: f64) -> State {
    let axes = [("x", LENGTH), ("y", LENGTH), ("z", LENGTH)];
    let axes = axes
        .into_iter()
        .chain(["rx", "ry", "rz"].map(|axis| (axis, Dims::NONE)));
    let declared = axes.enumerate().map(|(i, (axis, dims))| {
        let key = ParameterKey::new(format!("displacement_0_{axis}"), dims);
        (key, QtyAny::new(if i == 0 { x } else { 0.0 }, dims))
    });
    declared.collect()
}

/// `shift`'s keys on a `relaxation-bushing` card.
fn shift_keys(shift: &TemperatureShift) -> Vec<(&'static str, f64, Dims)> {
    let TemperatureShiftModel::Arrhenius {
        activation_energy_j_mol,
    } = shift.model()
    else {
        panic!("the module's shifts are Arrhenius");
    };
    let [lowest, highest] = shift.validity_k();
    vec![
        ("activation_energy", activation_energy_j_mol, MOLAR_ENERGY),
        (
            "shift_reference_temperature",
            shift.reference_temperature_k(),
            KELVIN,
        ),
        ("shift_lowest_temperature", lowest, KELVIN),
        ("shift_highest_temperature", highest, KELVIN),
    ]
}

/// 1b's pairing with bushing0 on the relaxation card, declaring `state`. The
/// bodies' work budget covers the dense solve's multiply-adds (about 3.4e8 a
/// step for this pair). `shift` adds its keys to the card.
fn with_relaxation(state: State, shift: &[(&str, f64, Dims)]) -> (DraftAssemblyInput, CardStore) {
    let (draft, mut cards) = pairing(fixture_cards(&laws()), 0.0);
    let mut input = draft.input().clone();
    let values = [relaxation_values().as_slice(), shift].concat();
    on_card(
        (&mut input, &mut cards),
        0,
        (RELAXATION_BUSHING_FAMILY, &values),
        state,
    );
    for body in &mut input.parts[..4] {
        body.fidelity.error_budget.max_work_units = Some(10_000_000_000);
    }
    (input, cards)
}

/// The relaxation card's values: the damping and one branch.
fn relaxation_values() -> Vec<(&'static str, f64, Dims)> {
    let mut values = vec![
        ("damping_x", 0.4, DAMPING),
        ("stiffness_0_x", 1.69, STIFFNESS),
        ("relaxation_time_0", 0.2, TIME),
    ];
    values.extend(["damping_y", "damping_z"].map(|name| (name, 0.0, DAMPING)));
    values.extend(["stiffness_0_y", "stiffness_0_z"].map(|name| (name, 0.0, STIFFNESS)));
    let turning = ["damping_rx", "damping_ry", "damping_rz"];
    values.extend(turning.map(|name| (name, 0.0, TORSIONAL_DAMPING)));
    let twisting = ["stiffness_0_rx", "stiffness_0_ry", "stiffness_0_rz"];
    values.extend(twisting.map(|name| (name, 0.0, TORSIONAL_STIFFNESS)));
    values
}

/// Puts interface `index` on a new `family` card of `values`, declaring `state`.
fn on_card(
    (input, cards): (&mut DraftAssemblyInput, &mut CardStore),
    index: usize,
    (family, values): (&str, &[(&str, f64, Dims)]),
    state: State,
) {
    let card = interface_card(family, values);
    let law_card = cards.insert_interface(family, card).unwrap();
    input.interfaces[index].kind = DraftInterfaceKind::MechanicalContact {
        law_card,
        third_body: None,
        state,
    };
}

/// The memory pair: `with_relaxation` without bodies 2 and 3, bushing1 and
/// its routes, and with bushing0's routes at the fixture's 0.7 and 0.2, and
/// `sources` temperature sources from bushing0's part to cell 0 on route00's faces.
fn memory_pair(
    state: State,
    shift: &[(&str, f64, Dims)],
    sources: usize,
) -> (DraftAssembly, CardStore) {
    let (mut input, cards) = with_relaxation(state, shift);
    input.parts.drain(2..4);
    // bushing0, the two beds, route00 and route01.
    input.interfaces = [0, 2, 3, 4, 5]
        .map(|i| input.interfaces[i].clone())
        .to_vec();
    for (route, f) in input.interfaces[3..].iter_mut().zip([0.7, 0.2]) {
        let state = vec![(
            ParameterKey::new("fraction", Dims::NONE),
            QtyAny::dimensionless(f),
        )];
        route.kind = DraftInterfaceKind::HeatRoute { state };
    }
    sourced((admit(input), cards), &vec![3; sources])
}

/// `pair` with bushing0 holding each `(port, history)` as a source.
fn holding(
    (draft, cards): (DraftAssembly, CardStore),
    held: &[(&str, TimeHistory)],
) -> (DraftAssembly, CardStore) {
    let mut input = draft.input().clone();
    let entity = input.interfaces[0].entity;
    for (n, (port, history)) in held.iter().enumerate() {
        input.sources.push(SourceDecl {
            name: format!("hold{n}"),
            target: SourceTarget::Port {
                entity,
                port: (*port).to_owned(),
            },
            history: history.clone(),
        });
    }
    (admit(input), cards)
}

#[test]
fn the_lowered_memory_draft_steps_with_the_memory_graph_bit_for_bit() {
    // M1 (MB 39274): the factor sqrt(1.69) is 1.3, and 1.3 times the declared
    // displacement is the fixture's 0.6 sqrt(J), each to the bit.
    assert_eq!(1.69_f64.sqrt().to_bits(), 1.3_f64.to_bits());
    assert_eq!((1.3 * DISPLACEMENT_M).to_bits(), 0.6_f64.to_bits());
    steps_with_the_memory_graph(None, None, false);
}

#[test]
fn a_shifted_memory_draft_follows_its_cell_as_the_shifted_memory_graph_bit_for_bit() {
    steps_with_the_memory_graph(Some(&Follows::Shift(arrhenius(300.0, 360.0))), None, false);
}

#[test]
fn a_memory_draft_on_a_damping_curve_steps_with_the_memory_graph_on_its_law_bit_for_bit() {
    let law = Follows::Damping(CurveLaw::new([module_curve(0.4)]));
    for held in [None, Some(0.7)] {
        steps_with_the_memory_graph(Some(&law), held, false);
    }
}

#[test]
fn a_held_memory_draft_steps_with_the_driven_memory_graph_bit_for_bit() {
    // The oracle's held forces (`HELD`), then 0.7 N with the shift.
    for force in [0.7, -0.7] {
        steps_with_the_memory_graph(None, Some(force), false);
    }
    let shift = Follows::Shift(arrhenius(300.0, 360.0));
    steps_with_the_memory_graph(Some(&shift), Some(0.7), false);
}

#[test]
fn an_equilibrium_spring_beside_the_branch_steps_with_the_memory_graph_bit_for_bit() {
    // The factor 0.5 times the 0.4 m stretch is the oracle's 0.2 sqrt(J).
    assert_eq!((0.25_f64.sqrt() * 0.4).to_bits(), 0.2_f64.to_bits());
    steps_with_the_memory_graph(None, None, true);
    // Held and shifted, as a measured mount is.
    let shift = Follows::Shift(arrhenius(300.0, 360.0));
    steps_with_the_memory_graph(Some(&shift), Some(0.7), true);
}

/// The memory pair following cell 0 as its temperature source by `follows`:
/// a shift's keys on its card, or bushing0's damping on `module_curve(0.4)`,
/// the claim the law must read; with `eq`, `memory(true)`'s spring.
fn following(follows: Option<&Follows>, eq: bool) -> (DraftAssembly, CardStore) {
    let (state, spring) = memory(eq);
    let pair = |keys: &[_]| memory_pair(state.clone(), &[&spring[..], keys].concat(), 1);
    match follows {
        None => memory_pair(state.clone(), &spring, 0),
        Some(Follows::Shift(shift)) => pair(&shift_keys(shift)),
        Some(Follows::Damping(_)) => curved(pair(&[]), 0, module_curve(0.4)),
    }
}

impl Follows {
    /// The kinetics' or the damping's factor at `kelvin` (for the curve,
    /// `ModuleLaw`'s, which it is to 2 ulp).
    fn factor(&self, kelvin: f64) -> f64 {
        match self {
            Self::Shift(shift) => shift.evaluate(kelvin).unwrap().factor(),
            Self::Damping(_) => 1.0 + 2.0 * (kelvin - 300.0),
        }
    }
}

/// `following(follows, eq)` with bushing0 holding `held` [N] on x, against
/// `memory_graph` following by `follows`, cell 0 published, that force and
/// with `eq` its spring.
fn steps_with_the_memory_graph(follows: Option<&Follows>, held: Option<f64>, eq: bool) {
    let pair = following(follows, eq);
    let force = held.map(|f| ("drive_x", TimeHistory::Constant(QtyAny::new(f, FORCE))));
    let (draft, cards) = holding(pair, force.as_slice());
    // The oracle takes the bodies' mass properties from the receipts, as 1b's.
    let mut fixture = oracle_pair(held, eq);
    let mass: [MassProperties; 2] = lowered_mass(&lower(&draft, &cards));
    fixture.bodies.properties[..2].copy_from_slice(&mass);
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(laws()));
    let controls = memory_controls();
    // Under an implicit scheme the draft drops its beds (`oracle_timeline`).
    for (scheme, target_s) in [
        (TimeScheme::Explicit, 0.25),
        (TimeScheme::Explicit, 0.3),
        (TimeScheme::BackwardEuler, 0.3),
        (TimeScheme::Sdirk2, 0.3),
    ] {
        let case = format!("{scheme:?} to {target_s} s, held {held:?}");
        let draft = under(&draft, scheme, scheme != TimeScheme::Explicit);
        let (mut timeline, heat) = oracle_timeline(laws, scheme);
        if follows.is_some() {
            timeline = timeline.publishing(vec![CellId(0)]).unwrap();
        }
        let follow = follows.cloned();
        let graph = memory_graph(&fixture, controls.mechanical, timeline, &heat, held, follow);
        let (mut oracle, network, cells) = graph;
        let mut lowered = lower(&draft, &cards);
        let nodes = ["relaxation-network", "thermal-timeline"].map(|name| lowered.node_map[name]);
        let system = &mut lowered.system;
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            while oracle.time_s() < target_s {
                let theirs = oracle.step(target_s - oracle.time_s(), cx).unwrap().dt_s;
                let step = system.step(target_s - system.time_s(), cx).unwrap().dt_s;
                assert_eq!(step.to_bits(), theirs.to_bits(), "{case}: the step");
                let ours = system.state::<RelaxationState>(nodes[0]);
                assert_eq!(ours, oracle.state(network), "{case}: bodies and memory");
                let bits = [
                    temperature_bits(system, nodes[1]),
                    temperature_bits(&oracle, cells),
                ];
                assert_eq!(bits[0], bits[1], "{case}: the cells");
            }
        });
        // The step histories are equal to the bit, so each account is too.
        let (got, want) = (system.ledger(), oracle.ledger());
        let dissipated = node(want, "network").dissipated_j.value();
        assert!(
            dissipated > 0.1,
            "{case}: enough heat to compare, {dissipated} J"
        );
        let account = |ledger: &Ledger, name: &str| {
            let account = node(ledger, name);
            [
                account.dissipated_j.value(),
                account.stored_change_j.value(),
            ]
            .map(f64::to_bits)
        };
        for (ours, theirs) in [
            ("relaxation-network", "network"),
            ("thermal-timeline", "cells"),
        ] {
            assert_eq!(account(got, ours), account(want, theirs), "{case}: {ours}");
        }
        for name in ["unallocated", "environment"] {
            assert_eq!(
                sink(got, name).to_bits(),
                sink(want, name).to_bits(),
                "{case}: {name}"
            );
        }
        // The actuator's work entered through its boundary, as the oracle's.
        let drive = |ledger: &Ledger, node: &str| {
            let port = format!("{node}.drive0");
            let account = ledger.boundaries.iter().find(|a| a.port == port);
            account.map(|account| account.energy_in_j.value().to_bits())
        };
        let work = drive(got, "relaxation-network");
        assert_eq!(work, drive(want, "network"), "{case}: the actuator's work");
        assert_eq!(work.is_some(), held.is_some(), "{case}: an actuator");
        // The residual is a sum of accounts that cancel: within 1e-12 of the
        // dissipated energy, as in 1b.
        let residual = got.energy_residual_j();
        assert!(
            residual.abs() <= 1e-12 * dissipated,
            "{case}: {residual:e} J"
        );
        // The kinetics or the damping moved: the factor at cell 0's last temperature.
        if let Some(follows) = follows {
            let kelvin = f64::from_bits(temperature_bits(system, nodes[1])[0]);
            let factor = follows.factor(kelvin);
            eprintln!("{case}: the factor is {factor} at cell 0's {kelvin} K");
            assert!(
                (factor - 1.0).abs() >= 0.05,
                "{case}: factor {factor} at {kelvin} K"
            );
        }
        // A body's reader reads its slot in the relaxation state.
        let state = system.state::<RelaxationState>(nodes[0]).unwrap();
        let v = state.bodies[1].linear_momentum_world().scale(0.5);
        let seen = observe(&lowered, draft.input().parts[1].entity);
        let Some(Observed::Body { velocity_m_s, .. }) = seen else {
            panic!("{case}: {seen:?}");
        };
        let want = [v.x, v.y, v.z].map(f64::to_bits);
        assert_eq!(velocity_m_s.map(f64::to_bits), want, "{case}");
    }
}

/// The branch's memory and, with `eq`, an equilibrium spring beside it: its
/// keys, 0.25 N/m on x (factor 0.5) and 0 on the other axes, and its memory,
/// stretched 0.4 m on x.
fn memory(eq: bool) -> (State, Vec<(&'static str, f64, Dims)>) {
    let mut state = displacement(DISPLACEMENT_M);
    if !eq {
        return (state, Vec::new());
    }
    let to_eq = |(key, q): (ParameterKey, _)| {
        let name = key.name.replace("_0_", "_eq_");
        (ParameterKey::new(name, key.dims), q)
    };
    state.extend(displacement(0.4).into_iter().map(to_eq));
    let mut spring = vec![("stiffness_eq_x", 0.25, STIFFNESS)];
    spring.extend(["stiffness_eq_y", "stiffness_eq_z"].map(|name| (name, 0.0, STIFFNESS)));
    let twisting = ["stiffness_eq_rx", "stiffness_eq_ry", "stiffness_eq_rz"];
    spring.extend(twisting.map(|name| (name, 0.0, TORSIONAL_STIFFNESS)));
    (state, spring)
}

/// `memory_graph`'s pair holding `held` [N] on x and, with `eq`, `memory(true)`'s
/// spring beside its branch: factor 0.5 on x, memory 0.2 sqrt(J).
fn oracle_pair(held: Option<f64>, eq: bool) -> MemoryPair {
    let mut pair = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, held.unwrap_or(0.0));
    if eq {
        let mut rows = [[0.0; 6]; 6];
        rows[0][0] = 0.5;
        let spring = SpatialRelaxationBranch::from_equilibrium_gram_rows(rows).unwrap();
        pair.branches[0].push(spring);
        pair.memory[0].push([0.2, 0.0, 0.0, 0.0, 0.0, 0.0]);
    }
    pair
}

/// The oracle's timeline and cells under `scheme`. The implicit one is a twin
/// of the timeline in `…reproduces_the_implicit_memory_candidate`, from the
/// same constructors, but its cells have no beds, which carry nothing in this
/// fixture: no reservoirs and no links. Its
/// controls are fs-thermal's defaults with the pairing's step, absolute
/// tolerance and work (39298 Q3).
fn oracle_timeline(
    laws: &'static [Caloric; 2],
    scheme: TimeScheme,
) -> (ThermalTimeline<'static>, DepositionSimulation<'static>) {
    let cells = vec![CellId(0), CellId(1)];
    let scheme = match scheme {
        TimeScheme::Explicit => {
            let timeline = ThermalTimeline::new(cells, true, memory_controls().thermal);
            return (timeline.unwrap(), thermal(laws, false, true));
        }
        TimeScheme::BackwardEuler => ImplicitScheme::BackwardEuler,
        TimeScheme::Sdirk2 => ImplicitScheme::Sdirk2,
    };
    let mut schedule = vec![birth(0, 0.0, 0.0), birth(1, 0.0, 0.0)];
    for event in &mut schedule {
        if let Event::Birth(cell) = &mut event.event {
            cell.boundary.reservoir = None;
        }
    }
    let heat = DepositionSimulation::new(
        schedule,
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        budget(),
    )
    .unwrap()
    .advance_to(0.0, config().thermal, || Ok(()))
    .unwrap()
    .0;
    let observer = ConductionObserver::new(
        &heat,
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        &[],
        &[],
        || Ok(()),
    )
    .unwrap();
    let graph = ImplicitConductionGraph::new(
        &heat,
        &observer,
        &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let controls = ImplicitConductionControls {
        maximum_dt_s: 0.125,
        maximum_steps: 10_000,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-10,
        ..ImplicitConductionControls::default()
    };
    let timeline = ThermalTimeline::implicit(cells, graph, controls, scheme);
    (timeline.unwrap(), heat)
}

/// Why lowering `draft` on `cards` refuses.
fn refusal((draft, cards): (DraftAssembly, CardStore)) -> String {
    let mut registry = LoweringRegistry::with_kernel_defaults();
    register_owners(&mut registry, fixture_laws());
    registry.lower(&draft, &cards).err().unwrap().to_string()
}

#[test]
fn a_shift_without_all_its_keys_one_model_or_one_source_refuses() {
    let keys = shift_keys(&arrhenius(300.0, 360.0));
    let wlf = [("wlf_c1", 8.86, Dims::NONE), ("wlf_c2", 101.6, KELVIN)];
    let both = [keys.as_slice(), &wlf].concat();
    let one = "a bushing that follows temperature reads exactly one `TemperatureSource`";
    let no_law = "its bushing's card states no temperature shift or damping curve";
    for (shift, sources, named) in [
        (&keys[..3], 1, "shift_highest_temperature"),
        (&both[..], 1, "a temperature shift has one model"),
        (&keys[..], 0, one),
        (&keys[..], 2, one),
        (&[][..], 1, no_law),
    ] {
        let why = refusal(memory_pair(displacement(DISPLACEMENT_M), shift, sources));
        assert!(why.contains(named), "{sources} sources: {why}");
    }
}

#[test]
fn shift_keys_on_a_viscous_card_or_in_a_partly_shifted_network_refuse() {
    let keys = shift_keys(&arrhenius(300.0, 360.0));
    // bushing0 on a viscous card with the keys: refused, not ignored.
    let (mut input, mut cards) = with_relaxation(Vec::new(), &[]);
    let damping = relaxation_values().into_iter();
    let damping = damping.filter(|(name, ..)| name.starts_with("damping"));
    let viscous = [damping.collect(), keys.clone()].concat();
    on_card(
        (&mut input, &mut cards),
        0,
        (VISCOUS_BUSHING_FAMILY, &viscous),
        Vec::new(),
    );
    let why = refusal((admit(input), cards));
    assert!(
        why.contains("a temperature shift needs a `relaxation-bushing` card"),
        "{why}"
    );
    // bushing1 on the relaxation card without the keys, beside bushing0 with them.
    let (mut input, mut cards) = with_relaxation(displacement(DISPLACEMENT_M), &keys);
    let unshifted = (RELAXATION_BUSHING_FAMILY, relaxation_values());
    let unshifted = (unshifted.0, unshifted.1.as_slice());
    on_card(
        (&mut input, &mut cards),
        1,
        unshifted,
        displacement(DISPLACEMENT_M),
    );
    let why = refusal((admit(input), cards));
    assert!(why.contains("all shifted or none"), "{why}");
}

#[test]
fn a_held_force_on_a_viscous_bushing_not_constant_in_other_units_or_twice_refuses() {
    let force = |f| TimeHistory::Constant(QtyAny::new(f, FORCE));
    let harmonic = TimeHistory::Harmonic {
        amplitude: QtyAny::new(0.7, FORCE),
        frequency: QtyAny::new(2.0, Dims([0, 0, -1, 0, 0, 0])),
        phase: QtyAny::dimensionless(0.0),
        offset: QtyAny::new(0.0, FORCE),
    };
    let couple = TimeHistory::Constant(QtyAny::new(0.7, TORSIONAL_STIFFNESS));
    let memory = || memory_pair(displacement(DISPLACEMENT_M), &[], 0);
    // 1b's pairing: both bushings on viscous cards.
    let viscous = pairing(fixture_cards(&laws()), 0.0);
    let twice = vec![force(0.7), force(0.1)];
    for (pair, held, named) in [
        (viscous, vec![force(0.7)], "needs a `relaxation-bushing`"),
        (memory(), vec![harmonic], "a held force is constant"),
        (memory(), vec![couple], "`drive_x` holds a force [N]"),
        (memory(), twice, "a bushing holds one source on `drive_x`"),
    ] {
        let held: Vec<_> = held.into_iter().map(|h| ("drive_x", h)).collect();
        // The last source is the one refused.
        let item = format!("refusal for source:hold{}: ", held.len() - 1);
        let why = refusal(holding(pair, &held));
        assert!(why.contains(&item) && why.contains(named), "{why}");
    }
}

#[test]
fn a_branch_without_its_memory_or_beside_a_viscous_bushing_refuses() {
    // The memory has no default: a missing key refuses by name.
    let mut state = displacement(DISPLACEMENT_M);
    state.pop();
    let why = refusal(memory_pair(state, &[], 0));
    assert!(
        why.contains("bushing0") && why.contains("declare `displacement_0_rz`"),
        "{why}"
    );
    // bushing1 stays on its viscous card: both sets and both families are named.
    let (input, cards) = with_relaxation(displacement(DISPLACEMENT_M), &[]);
    let why = refusal((admit(input), cards));
    for named in [
        "`bushing0` on `relaxation-bushing`",
        "`bushing1` on `viscous-bushing`",
    ] {
        assert!(why.contains(named), "{why}");
    }
}

#[test]
fn a_damping_curve_off_t_on_a_branch_beside_a_shift_or_a_fixed_bushing_refuses() {
    let keys = shift_keys(&arrhenius(300.0, 360.0));
    let shape = "`damping_x` follows temperature as a curve over `T`, `LinearInside`";
    let mut tabulated = module_curve(0.4);
    tabulated.interpolation = InterpolationPolicy::TabulatedOnly;
    let mut over_p = module_curve(0.4);
    over_p.value = PropertyValue::Curve {
        abscissa: "P".to_owned(),
        abscissa_dims: Dims([-1, 1, -2, 0, 0, 0]),
        knots: vec![(0.0, 0.4), (1e6, 0.4)],
        dims: DAMPING,
    };
    let flat = |name, dims, c| curve(name, dims, vec![(300.0, c), (303.0, c)]);
    let branch = |name| format!("`{name}` is a curve: a branch follows temperature by a shift");
    let none: &[_] = &[];
    let one = "a bushing that follows temperature reads exactly one `TemperatureSource`";
    let way = "a network follows temperature one way: a shift or damping curves";
    for (claim, shift, sources, named) in [
        (tabulated, none, 1, shape.to_owned()),
        (over_p, none, 1, shape.to_owned()),
        (
            flat("stiffness_0_x", STIFFNESS, 1.69),
            none,
            0,
            branch("stiffness_0_x"),
        ),
        (
            flat("relaxation_time_0", TIME, 0.2),
            none,
            0,
            branch("relaxation_time_0"),
        ),
        (module_curve(0.4), &keys, 1, way.to_owned()),
        (module_curve(0.4), none, 0, one.to_owned()),
        (module_curve(0.4), none, 2, one.to_owned()),
    ] {
        let pair = memory_pair(displacement(DISPLACEMENT_M), shift, sources);
        let why = refusal(curved(pair, 0, claim));
        assert!(why.contains(&named), "{sources} sources: {why}");
    }
    // 1b's pairing, bushing0's damping on a curve and bushing1's not.
    let pair = curved(pairing(fixture_cards(&laws()), 0.0), 0, module_curve(0.7));
    let why = refusal(sourced(pair, &[4]));
    let mixed = "a network's bushings all follow temperature or none";
    assert!(why.contains(mixed), "{why}");
}

#[test]
fn unstated_axes_step_as_explicit_zeros_and_a_partly_stated_card_refuses() {
    // bushing0 on a card of `values` alone, declaring `state`.
    let on = |values: &[_], state: State| {
        let (draft, mut cards) = memory_pair(displacement(DISPLACEMENT_M), &[], 0);
        let mut input = draft.input().clone();
        let family = (RELAXATION_BUSHING_FAMILY, values);
        on_card((&mut input, &mut cards), 0, family, state);
        (admit(input), cards)
    };
    // Only x stated: the damping, the branch and its memory.
    let x = relaxation_values()[..3].to_vec();
    let only_x = on(&x, displacement(DISPLACEMENT_M)[..1].to_vec());
    let zeros = memory_pair(displacement(DISPLACEMENT_M), &[], 0);
    let mut pairs = [only_x, zeros].map(|(draft, cards)| lower(&draft, &cards));
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        while pairs[1].system.time_s() < 0.3 {
            let seen = pairs.each_mut().map(|lowered| {
                let system = &mut lowered.system;
                let dt = system.step(0.3 - system.time_s(), cx).unwrap().dt_s;
                let node = lowered.node_map["relaxation-network"];
                (dt.to_bits(), system.state::<RelaxationState>(node).cloned())
            });
            assert_eq!(seen[0], seen[1], "at {} s", pairs[1].system.time_s());
        }
    });
    let sources = &pairs[0].receipt.interfaces[0].parameter_sources;
    let unstated = "not stated by the card: 0 for this bushing";
    for name in ["damping_y", "damping_rz", "stiffness_0_z", "stiffness_0_ry"] {
        let receipt = sources.iter().find(|p| p.parameter == name).unwrap();
        let ParameterSource::DeclaredAssemblyField { field, .. } = &receipt.source else {
            panic!("{name}: {receipt:?}");
        };
        assert_eq!(field, unstated, "{name}");
        assert_eq!(receipt.value.to_bits(), 0.0_f64.to_bits(), "{name}");
    }
    let (state, spring) = memory(true);
    let mut no_x = state.clone();
    no_x.retain(|(key, _)| key.name != "displacement_eq_x");
    let timed = [&spring[..], &[("relaxation_time_eq", 0.1, TIME)]].concat();
    let pairs = [
        memory_pair(state.clone(), &timed, 0),
        memory_pair(no_x, &spring, 0),
        memory_pair(state, &spring[..5], 0),
        on(&x[..1], Vec::new()),
        on(&x, displacement(DISPLACEMENT_M)),
    ];
    let named = [
        "`relaxation_time_eq`",
        "declare `displacement_eq_x`",
        "declare `stiffness_eq_rz`",
        "a branch's stiffness on some axis",
        "`displacement_0_y` is on an axis",
    ];
    for (pair, named) in pairs.into_iter().zip(named) {
        let why = refusal(pair);
        assert!(why.contains("bushing0") && why.contains(named), "{why}");
    }
}

#[test]
fn a_cell_past_its_damping_curves_last_knot_refuses_the_step_by_name() {
    // The knots end at 300.01 K, which cell 0 passes in the run.
    let knots = vec![(300.0, 0.4), (300.01, 0.4)];
    let pair = memory_pair(displacement(DISPLACEMENT_M), &[], 1);
    let (draft, cards) = curved(pair, 0, curve("damping_x", DAMPING, knots));
    let system = &mut lower(&draft, &cards).system;
    let gate = CancelGate::new_clock_free();
    let why = with_cx(&gate, |cx| {
        loop {
            assert!(system.time_s() < 0.5, "the run passed 300.01 K unrefused");
            if let Err(err) = system.step(0.5 - system.time_s(), cx) {
                break err.to_string();
            }
        }
    });
    eprintln!("the step refuses: {why}");
    assert!(why.contains("`damping_x`: OutsideKnotSpan"), "{why}");
}

#[test]
fn a_following_network_resumes_from_a_checkpoint_as_the_plain_run() {
    // As it is, then with an equilibrium spring beside the branch.
    for eq in [false, true] {
        let (state, spring) = memory(eq);
        resumes_as_the_plain_run(memory_pair(state, &spring, 1));
    }
}

/// `pair` on a damping curve, from a checkpoint two steps in, as the plain run.
fn resumes_as_the_plain_run(pair: (DraftAssembly, CardStore)) {
    let (draft, cards) = curved(pair, 0, module_curve(0.4));
    let mut lowered = lower(&draft, &cards);
    let nodes = ["relaxation-network", "thermal-timeline"].map(|name| lowered.node_map[name]);
    let system = &mut lowered.system;
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        // Two steps in: cell 0 is above 300 K, so the damping has moved.
        for _ in 0..2 {
            system.step(0.5 - system.time_s(), cx).unwrap();
        }
        let saved = system.checkpoint().unwrap();
        let run = |system: &mut System| {
            let mut seen = Vec::new();
            while system.time_s() < 0.5 {
                let dt = system.step(0.5 - system.time_s(), cx).unwrap().dt_s;
                let state = system.state::<RelaxationState>(nodes[0]).cloned();
                seen.push((dt.to_bits(), state, temperature_bits(system, nodes[1])));
            }
            (seen, format!("{:?}", system.ledger()))
        };
        let plain = run(system);
        assert!(
            plain.0.len() >= 2,
            "{} steps after the checkpoint",
            plain.0.len()
        );
        system.restore(&saved).unwrap();
        assert_eq!(run(system), plain);
    });
}
