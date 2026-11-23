//! The planar nonlinear winding owner as one Component: its ports, a step
//! that is the owner's own step bit for bit, heat that leaves through the
//! ports as the owner books it, the kernel's audit, conformance, refusals
//! and identity.
#![cfg(feature = "component")]

#[path = "nonlinear_winding_circuit/common.rs"]
mod common;

use common::{conductors, controls, profiles};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_em::planar::eddy::{EddyConductor, nonlinear::PlanarNonlinearEddy};
use fs_em::planar::nonlinear::PlanarNonlinearMagnetostatics;
use fs_em::winding::component::{
    HeatReceiver, PlanarWinding, PlanarWindingError, PlanarWindingState,
};
use fs_em::winding::nonlinear::{
    NonlinearWindingCircuit, NonlinearWindingInitial, NonlinearWindings,
};
use fs_em::winding::ohmic::WindingConductor;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{
    Carrier, Class, Component, EndOfStep, PortLayout, PortRole, PortSpec, Proposal, audit_step,
};

const DT: f64 = 0.2;
const VOLTS: f64 = 9.0;
const GRAPH_OHM: f64 = 0.3;
const TOLERANCE_J: f64 = 1e-9;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// Core conductors in triangles 2, 3, 6 and 7, at `kelvin`.
fn core(card: &MaterialCard, kelvin: f64) -> Vec<EddyConductor> {
    let material = common::material::conductor(card, kelvin, ElectricalPropertyBasis::Resistivity);
    [2, 3, 6, 7]
        .into_iter()
        .map(|triangle| EddyConductor {
            triangle,
            region: if triangle < 4 { 7 } else { 8 },
            conductor: material.clone(),
        })
        .collect()
}

/// Stranded conductors in triangles 0, 1, 4 and 5.
fn stranded(card: &MaterialCard) -> Vec<WindingConductor> {
    conductors(card, 300.0)
        .into_iter()
        .filter(|c| [0, 1, 4, 5].contains(&c.triangle))
        .collect()
}

/// A voltage source driving the winding into a graph resistor with a native
/// capacitor across it. The winding has an end inductance, and the section
/// carries core conductors, so every storage term and every kind of heat is
/// present.
fn series<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
) -> NonlinearWindingCircuit<'m> {
    let zeros = problem
        .solve(
            &[0.0; 8],
            &common::boundary(),
            &[0.0; 8],
            common::material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let field = PlanarNonlinearEddy::new_with_impressed_current(
        &zeros,
        core(card, 300.0),
        &[0, 1, 4, 5],
        &common::boundary(),
        &[0.0; 8],
        0.0,
        common::material::controls(),
        &mut || false,
    )
    .unwrap();
    NonlinearWindingCircuit::new(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: GRAPH_OHM }),
                (2, 0, Branch::Capacitor { farads: 0.05 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profiles(2).remove(0)],
            end_inductance_h: vec![0.05],
        },
        stranded(card),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![0.0],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap()
}

/// The left half, the right half and the graph resistor.
fn receivers() -> Vec<HeatReceiver> {
    vec![
        HeatReceiver {
            triangles: vec![0, 1, 2, 3],
            resistors: vec![],
        },
        HeatReceiver {
            triangles: vec![4, 5, 6, 7],
            resistors: vec![],
        },
        HeatReceiver {
            triangles: vec![],
            resistors: vec![2],
        },
    ]
}

fn component<'m>(circuit: &NonlinearWindingCircuit<'m>, card: &MaterialCard) -> PlanarWinding<'m> {
    PlanarWinding::new(
        circuit,
        core(card, 300.0),
        stranded(card),
        receivers(),
        controls(),
        TOLERANCE_J,
    )
    .unwrap()
}

#[test]
fn ports_follow_the_native_graph_and_no_port_tells_its_end() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let winding = component(&circuit, &card);
    let ports = winding.ports();
    let names: Vec<_> = ports.iter().map(|p| p.id().as_str().to_owned()).collect();
    assert_eq!(names, ["source0", "heat0", "heat1", "heat2"]);
    assert_eq!(
        ports[0].carrier(),
        Carrier::Power(PortKind::ElectricalVoltageCurrent)
    );
    assert_eq!(ports[0].role(), PortRole::EffortIn);
    for port in &ports[1..] {
        assert_eq!(port.carrier(), Carrier::Heat);
        assert_eq!(port.role(), PortRole::EffortIn);
    }
    assert_eq!(winding.class(), Class::Smooth);
    for port in 0..ports.len() {
        assert_eq!(winding.end_of_step(port), EndOfStep::Untold);
    }

    // A current source driving a resistor of its own is a `FlowIn` port.
    let profiles = profiles(1);
    let field = common::initial_field(&problem, &profiles, &[0.0]);
    let two = NonlinearWindingCircuit::new(
        field,
        CircuitGraph {
            node_count: 4,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: GRAPH_OHM }),
                (3, 0, Branch::CurrentSource { port: 1 }),
                (3, 0, Branch::Resistor { ohms: 2.0 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        conductors(&card, 300.0),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0, 0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap();
    let both = PlanarWinding::new(
        &two,
        vec![],
        conductors(&card, 300.0),
        vec![
            HeatReceiver {
                triangles: (0..8).collect(),
                resistors: vec![2],
            },
            HeatReceiver {
                triangles: vec![],
                resistors: vec![4],
            },
        ],
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let roles: Vec<_> = both.ports().iter().map(PortSpec::role).collect();
    assert_eq!(
        roles,
        [
            PortRole::EffortIn,
            PortRole::FlowIn,
            PortRole::EffortIn,
            PortRole::EffortIn
        ]
    );
    // The current source's port returns the voltage across it: 1.5 A into
    // 2 ohm. Its resistor's heat leaves through `heat1`.
    let mut next = PlanarWindingState::new(two.clone());
    let mut outputs = [0.0; 4];
    let state = PlanarWindingState::new(two);
    with_cx(|cx| {
        both.propose(
            &state,
            DT,
            &[VOLTS, 1.5, 300.0, 300.0],
            &mut next,
            &mut outputs,
            cx,
        )
    })
    .unwrap();
    let close = |a: f64, b: f64| assert!((a - b).abs() <= 1e-9 * (1.0 + b.abs()), "{a} != {b}");
    close(outputs[1].abs(), 3.0);
    close(outputs[3], -2.0 * 1.5 * 1.5);
}

#[test]
fn a_step_is_the_owners_own_step_and_its_heat_leaves_through_the_ports() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let winding = component(&circuit, &card);
    let layout = PortLayout::new(winding.ports()).unwrap();
    let kelvin = [300.0, 310.0, 320.0];
    let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
    let mut state = PlanarWindingState::new(circuit.clone());
    let mut owner = circuit;
    for n in 1..=3 {
        let mut next = state.clone();
        let mut outputs = [0.0; 4];
        let proposal =
            with_cx(|cx| winding.propose(&state, DT, &inputs, &mut next, &mut outputs, cx))
                .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("the winding shortened the step");
        };
        let expected = owner
            .step_to_candidate(
                owner.time_s() + DT,
                core(&card, 300.0),
                stranded(&card),
                &[VOLTS],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = expected.step();
        assert_eq!(next, PlanarWindingState::new(expected.circuit().clone()));
        assert_eq!(outputs[0].to_bits(), step.native.y[0].to_bits());
        // What the heat ports return is the owner's heat, to the rounding
        // floor the receipt claims.
        let returned: f64 = outputs[1..].iter().map(|rate| -rate * DT).sum();
        assert!(
            (returned - step.heat_j).abs() <= receipt.energy_tolerance_j,
            "ports return {returned} J, the owner booked {} J, floor {} J",
            step.heat_j,
            receipt.energy_tolerance_j
        );
        // The owner's records, booked to the receivers by hand.
        let mut heat = [0.0; 3];
        for cell in &step.cells {
            heat[usize::from(cell.source.triangle >= 4)] += cell.heat_j;
        }
        for cell in &step.field.cells {
            heat[usize::from(cell.source.triangle >= 4)] += cell.heat_j;
        }
        for resistor in &step.graph_resistors {
            assert_eq!(resistor.branch, 2);
            heat[2] += resistor.heat_j;
        }
        for (k, q) in heat.iter().enumerate() {
            assert!(*q > 0.0, "receiver {k} received no heat at step {n}");
            assert_eq!(outputs[1 + k].to_bits(), (-q / DT).to_bits());
        }
        assert!(!step.field.cells.is_empty() && !step.cells.is_empty());
        assert_eq!(
            receipt.stored_change_j.to_bits(),
            step.stored_change_j.to_bits()
        );
        let entropy: f64 = heat.iter().zip(kelvin).map(|(q, t)| q / t).sum();
        assert_eq!(
            receipt.entropy_generated_j_per_k.to_bits(),
            entropy.to_bits()
        );
        assert!(
            (receipt.numerical_loss_j - step.numerical_dissipation_j).abs()
                <= receipt.energy_tolerance_j + TOLERANCE_J
        );
        assert!(step.numerical_dissipation_j > 0.0);
        assert!(receipt.residual <= receipt.tolerance);
        // The kernel's first law, with the endpoints measured from the states.
        let old_j = winding.stored_energy_j(&state).unwrap();
        let new_j = winding.stored_energy_j(&next).unwrap();
        assert!((old_j - step.stored_before_j).abs() <= receipt.energy_tolerance_j);
        assert!((new_j - step.stored_after_j).abs() <= receipt.energy_tolerance_j);
        let audit = audit_step(
            winding.ports(),
            &layout,
            DT,
            &inputs,
            &outputs,
            &receipt,
            old_j,
            new_j,
        )
        .unwrap();
        println!(
            "component against owner: step {n}, current {:.9} A, heat {:?} J, stored {:.6e} J, \
             numerical {:.6e} J, balance defect {:.1e} J, endpoint defect {:.1e} J",
            next.circuit().currents_a()[0],
            heat,
            new_j,
            receipt.numerical_loss_j,
            audit.balance_defect_j,
            audit.endpoint_defect_j
        );
        owner = expected.into_parts().0;
        state = next;
    }
    assert!(owner.time_s() > 3.0 * DT - 1e-12);
}

#[test]
fn the_component_passes_the_conformance_kit() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let winding = component(&circuit, &card);
    let state = PlanarWindingState::new(circuit);
    let inputs = [VOLTS, 300.0, 310.0, 320.0];
    let mut stale = state.clone();
    let mut outputs = [0.0; 4];
    with_cx(|cx| winding.propose(&state, DT, &inputs, &mut stale, &mut outputs, cx)).unwrap();
    let case = Case {
        state: &state,
        stale: &stale,
        dt_s: DT,
        inputs: &inputs,
    };
    let observed =
        with_cx(|cx| check_component(&winding, &case, cx)).expect("component conformance");
    let Observed::Stepped { outputs, .. } = observed else {
        panic!("the winding shortened the step");
    };
    assert!(
        outputs[0] > 0.0 && outputs[1..].iter().all(|rate| *rate < 0.0),
        "{outputs:?}"
    );
}

#[test]
fn construction_and_steps_refuse_what_the_component_cannot_book() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let build = |receivers: Vec<HeatReceiver>, tolerance: f64| {
        PlanarWinding::new(
            &circuit,
            core(&card, 300.0),
            stranded(&card),
            receivers,
            controls(),
            tolerance,
        )
    };
    let refused = |result: Result<PlanarWinding<'_>, PlanarWindingError>, what: &str| match result {
        Err(PlanarWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    let mut missing = receivers();
    missing[2].resistors.clear();
    missing.truncate(2);
    refused(
        build(missing, TOLERANCE_J),
        "a dissipating element is mapped to no receiver",
    );
    let mut unmapped = receivers();
    unmapped[0].triangles.retain(|&triangle| triangle != 2);
    refused(
        build(unmapped, TOLERANCE_J),
        "a dissipating element is mapped to no receiver",
    );
    let mut twice = receivers();
    twice[1].triangles.push(0);
    refused(
        build(twice, TOLERANCE_J),
        "a dissipating element is mapped to two receivers",
    );
    let mut empty = receivers();
    empty.push(HeatReceiver::default());
    refused(build(empty, TOLERANCE_J), "a heat receiver maps no element");
    let mut inductor = receivers();
    inductor[2].resistors.push(1);
    refused(
        build(inductor, TOLERANCE_J),
        "a receiver names a branch that is not a resistor",
    );
    let mut foreign = receivers();
    foreign[0].triangles.push(8);
    refused(
        build(foreign, TOLERANCE_J),
        "a receiver names a triangle that carries no conductor",
    );
    refused(
        build(receivers(), 0.0),
        "energy tolerance must be positive and finite",
    );

    let winding = component(&circuit, &card);
    let state = PlanarWindingState::new(circuit.clone());
    let step = |inputs: &[f64], from: &PlanarWindingState<'_>| {
        let mut next = from.clone();
        let mut outputs = vec![0.0; inputs.len()];
        with_cx(|cx| winding.propose(from, DT, inputs, &mut next, &mut outputs, cx))
    };
    let invalid = |result: Result<Proposal, PlanarWindingError>, what: &str| match result {
        Err(PlanarWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    let ports = "positive duration, complete finite ports and positive temperatures required";
    invalid(step(&[VOLTS, 300.0, 0.0, 320.0], &state), ports);
    invalid(step(&[VOLTS, 300.0, 310.0], &state), ports);
    invalid(step(&[f64::NAN, 300.0, 310.0, 320.0], &state), ports);
    let other = PlanarWindingState::new(common::single(&problem, &card, GRAPH_OHM, 0.0, 0.0));
    invalid(
        step(&[VOLTS, 300.0, 310.0, 320.0], &other),
        "the state's circuit is not this Component's topology",
    );
    // A refused step leaves the accepted state as it was.
    assert_eq!(state, PlanarWindingState::new(circuit));
}

#[test]
fn identity_words_tell_conductors_receivers_and_tolerances_apart() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let words = |winding: &PlanarWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let base = words(&component(&circuit, &card));
    assert_eq!(base, words(&component(&circuit, &card)));
    let hot = PlanarWinding::new(
        &circuit,
        core(&card, 400.0),
        stranded(&card),
        receivers(),
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let mut swapped = receivers();
    swapped.swap(0, 1);
    let moved = PlanarWinding::new(
        &circuit,
        core(&card, 300.0),
        stranded(&card),
        swapped,
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let mut looser = controls();
    looser.power_tolerance_w *= 2.0;
    let loose = PlanarWinding::new(
        &circuit,
        core(&card, 300.0),
        stranded(&card),
        receivers(),
        looser,
        TOLERANCE_J,
    )
    .unwrap();
    let wider = PlanarWinding::new(
        &circuit,
        core(&card, 300.0),
        stranded(&card),
        receivers(),
        controls(),
        2.0 * TOLERANCE_J,
    )
    .unwrap();
    for other in [&hot, &moved, &loose, &wider] {
        assert_ne!(base, words(other));
    }
}
