//! The play-memory owner (`SpatialHystereticWindingCircuit`) as the 3D
//! Component's circuit, on the one-mode cube of
//! `tests/spatial_induction_winding/hysteretic_common.rs`: a step is the
//! owner's own, bit for bit, with every magnetic cell's pinning and viscous
//! heat at its own receiver, conducting or not; steps whose ports leave the
//! fits' temperature support refuse; refused and cancelled steps leave the
//! play memory as it was; a state whose field has another mesh is another
//! topology; the owner has a word of its own; the magnetic lane refuses it.

// Lints of the fixture's own lines, which this module does not edit.
#[allow(
    dead_code,
    clippy::elidable_lifetime_names,
    clippy::float_cmp,
    clippy::needless_range_loop
)]
#[path = "../spatial_induction_winding/hysteretic_common.rs"]
mod play;

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_em::spatial::vector_potential::induction::SpatialInductionLimits;
use fs_em::spatial::vector_potential::induction::hysteretic::winding::{
    SpatialHystereticWindingCircuit, SpatialHystereticWindingStep,
};
use fs_em::spatial::vector_potential::induction::hysteretic::{
    SpatialHystereticInduction, SpatialHystereticInductionControls, SpatialHystereticProblem,
};
use fs_em::spatial::vector_potential::induction::nonlinear::winding::SpatialNonlinearWindings;
use fs_em::winding::component::spatial::{
    SpatialConductorLaw, SpatialHeatReceiver, SpatialLawElement, SpatialWinding,
    SpatialWindingError, SpatialWindingState,
};
use fs_em::winding::end::EndWindingMode;
use fs_em::winding::nonlinear::NonlinearWindingPort;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::MaterialCard;
use fs_material::conductor::ResolvedConductorState;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, PortLayout, Proposal, Receipt, audit_step};
use play::{GRAPH_R, Mode, at_field, config, controls, initial, material, near, profile, strands};

use super::{DT, TOLERANCE_J, with_cx};

/// The Component's energy allowance: the owner's own, the
/// `energy_tolerance_j` of the fixture's controls. This owner closes its
/// account only to that, well beyond the saturating tests' 1e-9 J.
const ALLOWANCE_J: f64 = 4e-7;
/// The fits' intrinsic branch viscosity, as the 3D rate tests take it.
const ETA: f64 = 0.6;
/// Every fit holds from 200 K to this temperature; all are resolved at
/// 300 K.
const FIT_UPPER: f64 = 400.0;
const KELVIN: [f64; 4] = [300.0, 310.0, 320.0, 330.0];
/// The field \[A/m\] the fixture's one-mode reference drives each step to:
/// two minor loops.
const FIELDS: [f64; 7] = [12.0, 0.0, 6.0, 1.0, 6.0, -12.0, 0.0];
/// The core conducts in two of the six cells around the free diagonal. The
/// strands cross two others (24 and 25), and 28 and 29 are magnetic only.
const EDDY: [usize; 2] = [26, 27];

type PlayWinding<'m> = SpatialWinding<
    'm,
    ResolvedConductorState,
    SpatialHystereticWindingCircuit<'m>,
    SpatialHystereticInductionControls,
>;
type PlayState<'m> =
    SpatialWindingState<'m, ResolvedConductorState, SpatialHystereticWindingCircuit<'m>>;

fn card() -> Arc<MaterialCard> {
    material::card_with_viscosity(FIT_UPPER, ETA)
}

/// The fixture's cube, 2 x 2 x `layers` sub-cubes of six tetrahedra, every
/// cell's fit switching by pinning and viscosity. For `layers` from 2,
/// vertex 13 and the body diagonal [9, 22] are the fixture's. Vertices 12,
/// 18 and 19 of the diagonal's star leave the grid, so that its six cells,
/// congruent on the grid, each heat by their own amount.
fn problem(card: &Arc<MaterialCard>, layers: u32) -> SpatialHystereticProblem {
    let mut positions = Vec::new();
    for k in 0..=layers {
        for j in 0..=2 {
            for i in 0..=2 {
                positions.push([f64::from(i) / 2.0, f64::from(j) / 2.0, f64::from(k) / 2.0]);
            }
        }
    }
    positions[12][2] += 0.04;
    positions[18][2] += 0.07;
    positions[19][1] -= 0.03;
    let index = |i: u32, j: u32, k: u32| (k * 3 + j) * 3 + i;
    let mut tets = Vec::new();
    for k in 0..layers {
        for j in 0..2 {
            for i in 0..2 {
                let [a, b, c, d] = [
                    index(i, j, k),
                    index(i + 1, j, k),
                    index(i, j + 1, k),
                    index(i + 1, j + 1, k),
                ];
                let [e, f, g, h] = [
                    index(i, j, k + 1),
                    index(i + 1, j, k + 1),
                    index(i, j + 1, k + 1),
                    index(i + 1, j + 1, k + 1),
                ];
                tets.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    material::problem_with_viscosity(card, positions, tets)
}

/// The fixture's winding (a voltage source, the graph resistor as branch 1,
/// the winding as branch 2, end modes) on its one-mode field, with every
/// cochain but the body diagonal prescribed and the core conducting in
/// [`EDDY`].
fn circuit<'m>(
    problem: &'m SpatialHystereticProblem,
    card: &MaterialCard,
) -> SpatialHystereticWindingCircuit<'m> {
    let source = material::conductor(card, 300.0);
    let diagonal = problem.mesh().edges.binary_search(&[9, 22]).unwrap();
    let prescribed: Vec<usize> = (0..problem.mesh().edges.len())
        .filter(|&edge| edge != diagonal)
        .collect();
    let core = (0..problem.mesh().tets.len())
        .map(|cell| EDDY.contains(&cell).then(|| source.clone()))
        .collect();
    let field = SpatialHystereticInduction::new(
        problem,
        core,
        &prescribed,
        SpatialInductionLimits {
            maximum_unknowns: 2000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 100_000,
        },
        1e-10,
        &mut || false,
    )
    .unwrap();
    let profiles = vec![profile(problem)];
    SpatialHystereticWindingCircuit::new_with_end_modes(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: GRAPH_R }),
                (2, 0, Branch::Inductor { henries: f64::NAN }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![2],
            profiles_per_m2: profiles.clone(),
            end_inductance_h: vec![0.03],
        },
        vec![NonlinearWindingPort::CurrentInput],
        vec![EndWindingMode {
            inductance_h: 0.02,
            turns: vec![2.0],
        }],
        strands(&profiles, &source, 1.0),
        initial(1),
        config(),
        &mut || false,
    )
    .unwrap()
}

/// Source voltages \[V\]: the fixture's one-mode reference, whose core
/// conducts in four cells and whose fits have no viscosity, driven through
/// [`FIELDS`]. Here they only need to switch the play.
fn voltages() -> Vec<f64> {
    let mut old = Mode::default();
    FIELDS
        .iter()
        .map(|&h| {
            let next = at_field(old, 300.0, h);
            let volts = play::voltage(old, next, 300.0, GRAPH_R);
            old = next;
            volts
        })
        .collect()
}

fn stranded(circuit: &SpatialHystereticWindingCircuit<'_>) -> Vec<usize> {
    circuit
        .conductors()
        .iter()
        .map(|strand| strand.cell)
        .collect()
}

/// The strands' cells, the eddy cells, every other cell (magnetic only),
/// and the graph resistor.
fn receivers(circuit: &SpatialHystereticWindingCircuit<'_>) -> Vec<SpatialHeatReceiver> {
    let strands = stranded(circuit);
    let only = (0..circuit.field().problem().mesh().tets.len())
        .filter(|cell| !strands.contains(cell) && !EDDY.contains(cell))
        .collect();
    vec![
        SpatialHeatReceiver {
            cells: strands,
            resistors: vec![],
        },
        SpatialHeatReceiver {
            cells: EDDY.to_vec(),
            resistors: vec![],
        },
        SpatialHeatReceiver {
            cells: only,
            resistors: vec![],
        },
        SpatialHeatReceiver {
            cells: vec![],
            resistors: vec![1],
        },
    ]
}

fn component<'m>(
    circuit: &SpatialHystereticWindingCircuit<'m>,
    tolerance_j: f64,
) -> PlayWinding<'m> {
    SpatialWinding::new(
        circuit,
        circuit.field().conductors().to_vec(),
        circuit.conductors().to_vec(),
        receivers(circuit),
        controls(),
        tolerance_j,
    )
    .unwrap()
}

/// The owner's heat records booked to the four receivers by hand: stranded
/// cells, eddy cells, then every magnetic cell's pinning plus viscous heat,
/// by cell; the graph resistor to the fourth.
fn booked(step: &SpatialHystereticWindingStep, strands: &[usize]) -> [f64; 4] {
    let receiver = |cell: usize| {
        if strands.contains(&cell) {
            0
        } else if EDDY.contains(&cell) {
            1
        } else {
            2
        }
    };
    let mut heat = [0.0; 4];
    for cell in &step.cells {
        heat[receiver(cell.cell)] += cell.heat_j;
    }
    for cell in &step.field.cells {
        heat[receiver(cell.cell)] += cell.heat_j;
    }
    for cell in &step.field.magnetic_cells {
        heat[receiver(cell.cell)] += cell.irreversible_heat_j;
    }
    for resistor in &step.graph_resistors {
        assert_eq!(resistor.branch, 1);
        heat[3] += resistor.heat_j;
    }
    heat
}

/// A step's account and every heat record, as bits: the step record has no
/// `Debug`.
fn record(step: &SpatialHystereticWindingStep) -> Vec<u64> {
    [
        step.stored_before_j,
        step.stored_after_j,
        step.stored_change_j,
        step.heat_j,
        step.supplied_work_j,
        step.numerical_dissipation_j,
        step.field.stored_after_j,
        step.field.pinning_heat_j,
        step.field.viscous_heat_j,
    ]
    .into_iter()
    .chain(step.native.y.iter().copied())
    .chain(step.cells.iter().map(|cell| cell.heat_j))
    .chain(step.field.cells.iter().map(|cell| cell.heat_j))
    .chain(
        step.field
            .magnetic_cells
            .iter()
            .map(|cell| cell.irreversible_heat_j),
    )
    .map(f64::to_bits)
    .collect()
}

/// A context whose cancellation was requested before the step.
fn cancelled<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// `winding`'s accepted state after `volts`, each a step with every port at
/// 300 K.
fn driven<'m>(
    winding: &PlayWinding<'m>,
    circuit: SpatialHystereticWindingCircuit<'m>,
    volts: &[f64],
) -> PlayState<'m> {
    let mut state = SpatialWindingState::new(circuit).unwrap();
    for &v in volts {
        let mut next = state.clone();
        let mut outputs = [0.0; 5];
        with_cx(|cx| {
            winding.propose(
                &state,
                DT,
                &[v, 300.0, 300.0, 300.0, 300.0],
                &mut next,
                &mut outputs,
                cx,
            )
        })
        .unwrap();
        state = next;
    }
    state
}

/// Step `n`'s heat ports each return minus their receiver's heat rate as
/// booked by hand, bit for bit; each magnetic cell, in mesh order, heats by
/// its pinning plus viscous heat; the generated entropy is the ports' sum,
/// bit for bit; and the ports return the owner's heat within the allowance.
fn heat_leaves_as_booked(
    n: usize,
    step: &SpatialHystereticWindingStep,
    outputs: &[f64; 5],
    receipt: &Receipt,
    strands: &[usize],
) -> [f64; 4] {
    let heat = booked(step, strands);
    for (k, q) in heat.iter().enumerate() {
        assert_eq!(
            outputs[1 + k].to_bits(),
            (-q / DT).to_bits(),
            "receiver {k} at step {n}"
        );
    }
    for (index, cell) in step.field.magnetic_cells.iter().enumerate() {
        assert_eq!(cell.cell, index);
        near(
            cell.irreversible_heat_j,
            cell.pinning_heat_j + cell.viscous_heat_j,
            1e-15,
        );
    }
    let entropy: f64 = heat.iter().zip(KELVIN).map(|(q, t)| q / t).sum();
    assert_eq!(
        receipt.entropy_generated_j_per_k.to_bits(),
        entropy.to_bits()
    );
    let returned: f64 = heat.iter().sum();
    assert!(
        (returned - step.heat_j).abs() <= receipt.energy_tolerance_j,
        "ports return {returned} J, the owner booked {} J",
        step.heat_j
    );
    heat
}

#[test]
fn a_play_memory_step_is_the_owners_own_and_each_cell_heats_its_own_receiver() {
    let card = card();
    let p = problem(&card, 2);
    let circuit = circuit(&p, &card);
    let strands = stranded(&circuit);
    let winding = component(&circuit, ALLOWANCE_J);
    let layout = PortLayout::new(winding.ports()).unwrap();
    assert_eq!(
        controls().energy_tolerance_j.to_bits(),
        ALLOWANCE_J.to_bits()
    );
    let mut state = SpatialWindingState::new(circuit.clone()).unwrap();
    assert!(state.step().is_none());
    let mut native = circuit;
    let (mut insulating, mut viscous) = (0.0, 0.0);
    for (n, volts) in voltages().into_iter().enumerate() {
        let inputs = [volts, KELVIN[0], KELVIN[1], KELVIN[2], KELVIN[3]];
        let mut next = state.clone();
        let mut outputs = [0.0; 5];
        let proposal =
            with_cx(|cx| winding.propose(&state, DT, &inputs, &mut next, &mut outputs, cx))
                .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("the winding shortened the step");
        };
        let expected = native
            .step_with_conductors_candidate(
                native.time_s() + DT,
                native.field().conductors().to_vec(),
                native.conductors().to_vec(),
                &[volts],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = expected.step();
        let rebuilt = SpatialWindingState::new(expected.circuit().clone()).unwrap();
        assert_eq!(next, rebuilt);
        assert_eq!(
            rebuilt.stored_energy_j().to_bits(),
            step.stored_after_j.to_bits()
        );
        assert_eq!(record(next.step().unwrap()), record(step));
        assert_eq!(
            state.stored_energy_j().to_bits(),
            step.stored_before_j.to_bits()
        );
        assert_eq!(outputs[0].to_bits(), step.native.y[0].to_bits());
        let heat = heat_leaves_as_booked(n, step, &outputs, &receipt, &strands);
        // Off the grid, no two cells of the star heat alike, so a heat
        // booked at another cell's receiver shows in the ports.
        let star: Vec<u64> = step.field.magnetic_cells[24..30]
            .iter()
            .map(|cell| cell.irreversible_heat_j.to_bits())
            .collect();
        assert!(
            (1..star.len()).all(|i| !star[..i].contains(&star[i])),
            "step {n}: two cells of the star heat alike, {star:?}"
        );
        let only = &step.field.magnetic_cells[28..30];
        insulating += only
            .iter()
            .map(|cell| cell.irreversible_heat_j)
            .sum::<f64>();
        viscous += step.field.viscous_heat_j;
        let old_j = winding.stored_energy_j(&state).unwrap();
        let new_j = winding.stored_energy_j(&next).unwrap();
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
            "3D play memory against owner: step {n}, {volts:.6} V, current {:.9} A, heat \
             {heat:?} J, cells 28 and 29 {:.6e} J, viscous {:.6e} J, owner defect {:.3e} J \
             (the owner's residual {:.3e} J, floor {:.1e} J), balance defect {:.1e} J, \
             endpoint defect {:.1e} J",
            next.circuit().currents_a()[0],
            only.iter()
                .map(|cell| cell.irreversible_heat_j)
                .sum::<f64>(),
            step.field.viscous_heat_j,
            receipt.numerical_loss_j - step.numerical_dissipation_j,
            step.energy_residual_j,
            receipt.energy_tolerance_j,
            audit.balance_defect_j,
            audit.endpoint_defect_j
        );
        native = expected.into_parts().0;
        state = next;
    }
    assert!(
        insulating > 0.0 && viscous > 0.0,
        "cells 28 and 29 {insulating} J, viscous {viscous} J"
    );
}

#[test]
fn the_play_memory_component_passes_the_conformance_kit() {
    let card = card();
    let p = problem(&card, 2);
    let circuit = circuit(&p, &card);
    let winding = component(&circuit, ALLOWANCE_J);
    let volts = voltages();
    let state = driven(&winding, circuit, &volts[..2]);
    let inputs = [volts[2], KELVIN[0], KELVIN[1], KELVIN[2], KELVIN[3]];
    let mut stale = state.clone();
    let mut outputs = [0.0; 5];
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
    assert!(outputs[1..].iter().all(|rate| *rate < 0.0), "{outputs:?}");
}

#[test]
fn ports_beyond_the_fits_support_refuse_and_refused_steps_leave_the_memory_as_it_was() {
    let card = card();
    let p = problem(&card, 2);
    let start = circuit(&p, &card);
    let only = receivers(&start)[2].cells[0];
    let winding = component(&start, ALLOWANCE_J);
    let volts = voltages();
    let state = driven(&winding, start.clone(), &volts[..2]);
    let memory = state.circuit().field().cell_states();
    assert!(
        memory[28].polarizations_t()[0]
            .iter()
            .any(|j| j.abs() > 0.0)
    );
    let attempt = |winding: &PlayWinding<'_>, inputs: [f64; 5], cancel: bool| {
        let mut next = state.clone();
        let mut outputs = [0.0; 5];
        let step = |cx: &Cx<'_>| winding.propose(&state, DT, &inputs, &mut next, &mut outputs, cx);
        let result = if cancel {
            cancelled(step)
        } else {
            with_cx(step)
        };
        if result.is_err() {
            assert_eq!(
                next, state,
                "a refused step changed the state it was handed"
            );
        }
        result
    };
    let law = |result: Result<Proposal, SpatialWindingError>, port: usize, cell: usize| match result
    {
        Err(SpatialWindingError::Law {
            receiver,
            element: SpatialLawElement::Field(at),
            kelvin,
            message,
        }) if receiver == port && at == cell => {
            println!("field cell {at} at {kelvin} K refused: {message}");
        }
        other => panic!("expected cell {cell}'s fit to refuse, got {other:?}"),
    };
    let v = volts[2];
    // The support is inclusive on both sides. Each refusal names the first
    // cell, in mesh order, of the receiver whose port leaves it.
    attempt(&winding, [v, 400.0, 200.0, 300.0, 300.0], false).unwrap();
    law(
        attempt(&winding, [v, 420.0, 300.0, 300.0, 300.0], false),
        0,
        0,
    );
    law(
        attempt(&winding, [v, 300.0, 150.0, 300.0, 300.0], false),
        1,
        26,
    );
    law(
        attempt(&winding, [v, 300.0, 300.0, 400.5, 300.0], false),
        2,
        only,
    );
    // The graph resistor's port reaches no fit.
    attempt(&winding, [v, 300.0, 300.0, 300.0, 600.0], false).unwrap();
    assert!(matches!(
        attempt(&winding, [v, 300.0, 300.0, 300.0, 300.0], true),
        Err(SpatialWindingError::Cancelled)
    ));
    // An allowance at the rounding floor alone refuses the owner's own
    // closure defect after the step was solved: from rest, the first step's
    // exceeds its floor.
    let rest = SpatialWindingState::new(start).unwrap();
    let strict = component(rest.circuit(), f64::MIN_POSITIVE);
    let mut next = rest.clone();
    let mut outputs = [0.0; 5];
    let inputs = [volts[0], 300.0, 300.0, 300.0, 300.0];
    match with_cx(|cx| strict.propose(&rest, DT, &inputs, &mut next, &mut outputs, cx)) {
        Err(SpatialWindingError::Balance {
            defect_j,
            tolerance_j,
        }) => println!("floor alone: defect {defect_j:.3e} J beyond {tolerance_j:.3e} J"),
        other => panic!("expected the balance refusal, got {other:?}"),
    }
    assert_eq!(next, rest, "a refused step changed the state it was handed");

    // A state whose field has another mesh: three layers, 72 cells.
    let tall = problem(&card, 3);
    let foreign = SpatialWindingState::new(circuit(&tall, &card)).unwrap();
    let mut next = foreign.clone();
    let mut outputs = [0.0; 5];
    match with_cx(|cx| {
        winding.propose(
            &foreign,
            DT,
            &[v, 300.0, 300.0, 300.0, 300.0],
            &mut next,
            &mut outputs,
            cx,
        )
    }) {
        Err(SpatialWindingError::Invalid { what }) => {
            assert_eq!(what, "the state's circuit is not this Component's topology");
        }
        other => panic!("expected the topology refusal, got {other:?}"),
    }
}

#[test]
fn every_magnetic_cell_maps_to_one_receiver_and_the_owner_has_its_own_word() {
    let card = card();
    let p = problem(&card, 2);
    let circuit = circuit(&p, &card);
    let build = |receivers: Vec<SpatialHeatReceiver>| {
        SpatialWinding::new(
            &circuit,
            circuit.field().conductors().to_vec(),
            circuit.conductors().to_vec(),
            receivers,
            controls(),
            ALLOWANCE_J,
        )
    };
    let refused = |result: Result<PlayWinding<'_>, SpatialWindingError>, what: &str| match result {
        Err(SpatialWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    // Cells 28 and 29 carry no conductor; their fits still dissipate.
    let mut without = receivers(&circuit);
    without[2].cells.retain(|&cell| cell != 28);
    refused(
        build(without),
        "a dissipating element is mapped to no receiver",
    );
    let mut twice = receivers(&circuit);
    twice[1].cells.push(28);
    refused(
        build(twice),
        "a dissipating element is mapped to two receivers",
    );
    let mut beyond = receivers(&circuit);
    beyond[2].cells.push(48);
    refused(
        build(beyond),
        "a receiver names a cell that carries no conductor",
    );

    let words = |winding: &PlayWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let base = words(&component(&circuit, ALLOWANCE_J));
    assert_eq!(base, words(&component(&circuit, ALLOWANCE_J)));
    // The tolerance, then the owner's word.
    let owner = u64::from_le_bytes(*b"vecplay1");
    assert_eq!(base[base.len() - 2..], [ALLOWANCE_J.to_bits(), owner]);
    let saturating = super::card(18.0);
    let q = super::cube(&saturating);
    let mut other = Vec::new();
    super::component(&super::series(&q, &saturating), &q, &saturating).identity_words(&mut other);
    assert_eq!(other.last(), Some(&TOLERANCE_J.to_bits()));
    // An eddy cell and a magnetic-only cell in each other's receiver.
    let mut swapped = receivers(&circuit);
    let at = swapped[2]
        .cells
        .iter()
        .position(|&cell| cell == 28)
        .unwrap();
    swapped[1].cells[0] = 28;
    swapped[2].cells[at] = EDDY[0];
    assert_ne!(base, words(&build(swapped).unwrap()));
}

/// The circuit's own conductor in each core cell and strand, at every
/// temperature.
#[derive(Debug)]
struct Own {
    core: BTreeMap<usize, ResolvedConductorState>,
    strands: BTreeMap<(usize, usize), ResolvedConductorState>,
}

impl SpatialConductorLaw<ResolvedConductorState> for Own {
    fn core(&self, cell: usize, _: f64) -> Result<ResolvedConductorState, String> {
        self.core
            .get(&cell)
            .cloned()
            .ok_or_else(|| format!("no core conductor in cell {cell}"))
    }

    fn strand(
        &self,
        winding: usize,
        cell: usize,
        _: f64,
    ) -> Result<ResolvedConductorState, String> {
        self.strands
            .get(&(winding, cell))
            .cloned()
            .ok_or_else(|| format!("no strand of winding {winding} in cell {cell}"))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([self.core.len() as u64, self.strands.len() as u64]);
    }
}

#[test]
fn laws_that_answer_the_owners_conductors_step_as_fixed_conductors_do() {
    let card = card();
    let p = problem(&card, 2);
    let circuit = circuit(&p, &card);
    let fixed = component(&circuit, ALLOWANCE_J);
    let own: Arc<dyn SpatialConductorLaw<ResolvedConductorState>> = Arc::new(Own {
        core: circuit
            .field()
            .conductors()
            .iter()
            .enumerate()
            .filter_map(|(cell, conductor)| Some((cell, conductor.clone()?)))
            .collect(),
        strands: circuit
            .conductors()
            .iter()
            .map(|strand| ((strand.winding, strand.cell), strand.conductor.clone()))
            .collect(),
    });
    let following: PlayWinding<'_> = SpatialWinding::at_port_temperatures(
        &circuit,
        vec![own; 4],
        BTreeMap::new(),
        receivers(&circuit),
        controls(),
        ALLOWANCE_J,
    )
    .unwrap();
    let (mut a, mut b) = (
        SpatialWindingState::new(circuit.clone()).unwrap(),
        SpatialWindingState::new(circuit).unwrap(),
    );
    for volts in voltages().into_iter().take(3) {
        let inputs = [volts, KELVIN[0], KELVIN[1], KELVIN[2], KELVIN[3]];
        let (mut next_a, mut next_b) = (a.clone(), b.clone());
        let (mut out_a, mut out_b) = ([0.0; 5], [0.0; 5]);
        with_cx(|cx| fixed.propose(&a, DT, &inputs, &mut next_a, &mut out_a, cx)).unwrap();
        with_cx(|cx| following.propose(&b, DT, &inputs, &mut next_b, &mut out_b, cx)).unwrap();
        assert_eq!(next_a, next_b);
        assert_eq!(out_a.map(f64::to_bits), out_b.map(f64::to_bits));
        (a, b) = (next_a, next_b);
    }
}

#[test]
fn the_play_memory_owner_has_no_magnetic_lane() {
    let card = card();
    let p = problem(&card, 2);
    let circuit = circuit(&p, &card);
    let refused = PlayWinding::magnetic_at_port_temperatures(
        &circuit,
        vec![],
        BTreeMap::new(),
        vec![],
        receivers(&circuit),
        controls(),
        ALLOWANCE_J,
    );
    match refused {
        Err(SpatialWindingError::Invalid { what }) => {
            assert_eq!(
                what,
                "this owner refreshes no field materials with temperature"
            );
        }
        other => panic!("expected the lane to refuse the play-memory owner, got {other:?}"),
    }
}
