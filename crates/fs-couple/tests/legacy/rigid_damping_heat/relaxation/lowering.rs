//! The memory pair declared as a draft and lowered by the owner constructors
//! (`fs_system::lowering::owners`), against `memory_graph` above with no
//! actuator. It is 1b's first pair of 2 kg boxes closing at 2 m/s, its bushing
//! on a `relaxation-bushing` card: Kelvin viscosity 0.4 N s/m and one branch
//! of stiffness 1.69 N/m (factor 1.3) and time 0.2 s, all on x. The interface
//! declares the branch's elastic displacement, 0.6 / 1.3 m, and two loss
//! routes send 0.7 and 0.2 of its heat to the cells. Fixture values; not a
//! seed-dataset authority.

use super::*;
use crate::graph::lowering::{
    fixture_cards, fixture_laws, interface_card, lower, lowered_mass, observe, pairing,
    temperature_bits,
};
use crate::graph::with_cx;
use fs_exec::CancelGate;
use fs_qty::{Dims, QtyAny};
use fs_scenario::Scenario;
use fs_scenario::assembly::ParameterKey;
use fs_scenario::assembly::draft::{DraftAssembly, DraftAssemblyInput, DraftInterfaceKind};
use fs_system::lowering::Observed;
use fs_system::lowering::owners::{RELAXATION_BUSHING_FAMILY, register_owners};
use fs_system::{CardStore, Ledger, LoweringRegistry};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const DAMPING: Dims = Dims([0, 1, -1, 0, 0, 0]);
const TORSIONAL_DAMPING: Dims = Dims([2, 1, -1, 0, 0, 0]);
const STIFFNESS: Dims = Dims([0, 1, -2, 0, 0, 0]);
const TORSIONAL_STIFFNESS: Dims = Dims([2, 1, -2, 0, 0, 0]);
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

/// 1b's pairing with bushing0 on the relaxation card, declaring `state`. The
/// bodies' work budget covers the dense solve's multiply-adds (about 3.4e8 a
/// step for this pair).
fn with_relaxation(state: State) -> (DraftAssemblyInput, CardStore) {
    let (draft, mut cards) = pairing(fixture_cards(&laws()), 0.0);
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
    let card = interface_card(RELAXATION_BUSHING_FAMILY, &values);
    let law_card = cards
        .insert_interface(RELAXATION_BUSHING_FAMILY, card)
        .unwrap();
    let mut input = draft.input().clone();
    input.interfaces[0].kind = DraftInterfaceKind::MechanicalContact {
        law_card,
        third_body: None,
        state,
    };
    for body in &mut input.parts[..4] {
        body.fidelity.error_budget.max_work_units = Some(10_000_000_000);
    }
    (input, cards)
}

fn admit(input: DraftAssemblyInput) -> DraftAssembly {
    let scenario = Scenario::new("pairing", 112, input.environment.base.clone());
    DraftAssembly::admit(input, &scenario).unwrap()
}

/// The memory pair: `with_relaxation` without bodies 2 and 3, bushing1 and
/// its routes, and with bushing0's routes at the fixture's 0.7 and 0.2.
fn memory_pair(state: State) -> (DraftAssembly, CardStore) {
    let (mut input, cards) = with_relaxation(state);
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
    (admit(input), cards)
}

#[test]
fn the_lowered_memory_draft_steps_with_the_memory_graph_bit_for_bit() {
    // M1 (MB 39274): the factor sqrt(1.69) is 1.3, and 1.3 times the declared
    // displacement is the fixture's 0.6 sqrt(J), each to the bit.
    assert_eq!(1.69_f64.sqrt().to_bits(), 1.3_f64.to_bits());
    assert_eq!((1.3 * DISPLACEMENT_M).to_bits(), 0.6_f64.to_bits());
    let (draft, cards) = memory_pair(displacement(DISPLACEMENT_M));
    // The oracle takes the bodies' mass properties from the receipts, as 1b's.
    let mut fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, 0.0);
    let mass: [MassProperties; 2] = lowered_mass(&lower(&draft, &cards));
    fixture.bodies.properties[..2].copy_from_slice(&mass);
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(laws()));
    let (heat, controls) = (thermal(laws, false, true), memory_controls());
    for target_s in [0.25, 0.3] {
        let case = format!("to {target_s} s");
        let timeline = ThermalTimeline::new(vec![CellId(0), CellId(1)], true, controls.thermal);
        let timeline = timeline.unwrap();
        let graph = memory_graph(&fixture, controls.mechanical, timeline, &heat, None, None);
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

#[test]
fn a_branch_without_its_memory_or_beside_a_viscous_bushing_refuses() {
    let refusal = |(draft, cards): (DraftAssembly, CardStore)| {
        let mut registry = LoweringRegistry::with_kernel_defaults();
        register_owners(&mut registry, fixture_laws());
        registry.lower(&draft, &cards).err().unwrap().to_string()
    };
    // The memory has no default: a missing key refuses by name.
    let mut state = displacement(DISPLACEMENT_M);
    state.pop();
    let why = refusal(memory_pair(state));
    assert!(
        why.contains("bushing0") && why.contains("declare `displacement_0_rz`"),
        "{why}"
    );
    // bushing1 stays on its viscous card: both sets and both families are named.
    let (input, cards) = with_relaxation(displacement(DISPLACEMENT_M));
    let why = refusal((admit(input), cards));
    for named in [
        "`bushing0` on `relaxation-bushing`",
        "`bushing1` on `viscous-bushing`",
    ] {
        assert!(why.contains(named), "{why}");
    }
}
