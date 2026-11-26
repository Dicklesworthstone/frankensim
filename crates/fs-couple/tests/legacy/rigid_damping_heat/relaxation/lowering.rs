//! The memory pair declared as a draft and lowered by the owner constructors
//! (`fs_system::lowering::owners`), against `memory_graph` above with no
//! actuator. It is 1b's first pair of 2 kg boxes closing at 2 m/s, its bushing
//! on a `relaxation-bushing` card: Kelvin viscosity 0.4 N s/m and one branch
//! of stiffness 1.69 N/m (factor 1.3) and time 0.2 s, all on x. The interface
//! declares the branch's elastic displacement, 0.6 / 1.3 m, and two loss
//! routes send 0.7 and 0.2 of its heat to the cells. Fixture values; not a
//! seed-dataset authority. With the module's Arrhenius shift on the card and
//! cell 0 as its `TemperatureSource`, against `memory_graph` with that shift
//! on cell 0's published temperature.

use super::*;
use crate::graph::lowering::{
    fixture_cards, fixture_laws, interface_card, lower, lowered_mass, observe, pairing,
    temperature_bits, under,
};
use crate::graph::with_cx;
use fs_exec::CancelGate;
use fs_material::visco::temperature::TemperatureShiftModel;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::ParameterKey;
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftAssemblyInput, DraftInterfaceKind, TimeScheme,
};
use fs_scenario::{EntityDeclaration, InterfacePair, Scenario};
use fs_system::lowering::Observed;
use fs_system::lowering::owners::{
    RELAXATION_BUSHING_FAMILY, VISCOUS_BUSHING_FAMILY, register_owners,
};
use fs_system::{CardStore, Ledger, LoweringRegistry};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const DAMPING: Dims = Dims([0, 1, -1, 0, 0, 0]);
const TORSIONAL_DAMPING: Dims = Dims([2, 1, -1, 0, 0, 0]);
const STIFFNESS: Dims = Dims([0, 1, -2, 0, 0, 0]);
const TORSIONAL_STIFFNESS: Dims = Dims([2, 1, -2, 0, 0, 0]);
const KELVIN: Dims = Dims([0, 0, 0, 1, 0, 0]);
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
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

fn admit(input: DraftAssemblyInput) -> DraftAssembly {
    let scenario = Scenario::new("pairing", 112, input.environment.base.clone());
    DraftAssembly::admit(input, &scenario).unwrap()
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
    for n in 0..sources {
        let mut source = input.interfaces[3].clone();
        let pair = InterfacePair::ordered(source.feature_a, source.feature_b);
        source.name = format!("source{n}");
        let declared = EntityDeclaration::interface(input.root, &source.name, pair);
        source.entity = input.catalog.declare(declared).unwrap();
        source.kind = DraftInterfaceKind::TemperatureSource;
        input.interfaces.push(source);
    }
    (admit(input), cards)
}

#[test]
fn the_lowered_memory_draft_steps_with_the_memory_graph_bit_for_bit() {
    // M1 (MB 39274): the factor sqrt(1.69) is 1.3, and 1.3 times the declared
    // displacement is the fixture's 0.6 sqrt(J), each to the bit.
    assert_eq!(1.69_f64.sqrt().to_bits(), 1.3_f64.to_bits());
    assert_eq!((1.3 * DISPLACEMENT_M).to_bits(), 0.6_f64.to_bits());
    steps_with_the_memory_graph(None);
}

#[test]
fn a_shifted_memory_draft_follows_its_cell_as_the_shifted_memory_graph_bit_for_bit() {
    steps_with_the_memory_graph(Some(arrhenius(300.0, 360.0)));
}

/// The memory pair, with `shift` on its card and cell 0 as its temperature
/// source, against `memory_graph` with that shift and cell 0 published.
fn steps_with_the_memory_graph(shift: Option<TemperatureShift>) {
    let keys = shift.as_ref().map(shift_keys).unwrap_or_default();
    let sources = usize::from(shift.is_some());
    let (draft, cards) = memory_pair(displacement(DISPLACEMENT_M), &keys, sources);
    // The oracle takes the bodies' mass properties from the receipts, as 1b's.
    let mut fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, 0.0);
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
        let case = format!("{scheme:?} to {target_s} s");
        let draft = under(&draft, scheme, scheme != TimeScheme::Explicit);
        let (mut timeline, heat) = oracle_timeline(laws, scheme);
        if shift.is_some() {
            timeline = timeline.publishing(vec![CellId(0)]).unwrap();
        }
        let graph = memory_graph(&fixture, controls.mechanical, timeline, &heat, None, shift);
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
        // The residual is a sum of accounts that cancel: within 1e-12 of the
        // dissipated energy, as in 1b.
        let residual = got.energy_residual_j();
        assert!(
            residual.abs() <= 1e-12 * dissipated,
            "{case}: {residual:e} J"
        );
        // The kinetics moved: the factor at cell 0's last temperature.
        if let Some(shift) = shift {
            let kelvin = f64::from_bits(temperature_bits(system, nodes[1])[0]);
            let factor = shift.evaluate(kelvin).unwrap().factor();
            eprintln!("{case}: the shift factor is {factor} at cell 0's {kelvin} K");
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
    let one = "a shifted bushing reads exactly one `TemperatureSource`";
    for (shift, sources, named) in [
        (&keys[..3], 1, "shift_highest_temperature"),
        (&both[..], 1, "a temperature shift has one model"),
        (&keys[..], 0, one),
        (&keys[..], 2, one),
        (&[][..], 1, "its bushing's card states no temperature shift"),
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
