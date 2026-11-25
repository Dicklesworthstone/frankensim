//! The play-memory owner (`HystereticWindingCircuit`) as the Component's
//! circuit: a step is the owner's own, bit for bit, with every magnetic
//! triangle's pinning and viscous heat at its own receiver, conducting or
//! not; steps whose ports leave the fits' temperature support refuse;
//! refused and cancelled steps leave the vector memory as it was; the lane
//! has a word of its own; the magnetic lane refuses the play-memory owner.

#[path = "../hysteretic_eddy/common.rs"]
mod fixture;

use std::collections::BTreeMap;
use std::sync::Arc;

use fixture::{
    MU, boundary, controls, near, paired_squares, problem_from_materials, sources_at, square,
    terminal_controls,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_em::planar::eddy::EddyConductor;
use fs_em::planar::eddy::hysteretic::{HystereticEddyControls, PlanarHystereticEddy};
use fs_em::planar::nonlinear::hysteretic::PlanarHystereticMagnetostatics;
use fs_em::winding::component::{
    ConductorLaw, HeatReceiver, LawElement, PlanarWinding, PlanarWindingError, PlanarWindingState,
};
use fs_em::winding::end::EndWindingMode;
use fs_em::winding::hysteretic::{HystereticWindingCircuit, HystereticWindingStep};
use fs_em::winding::nonlinear::{
    NonlinearWindingConfig, NonlinearWindingInitial, NonlinearWindingPort, NonlinearWindings,
};
use fs_em::winding::ohmic::WindingConductor;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::conductor::ResolvedConductorState;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, PortLayout, Proposal, audit_step};

use super::{DT, GRAPH_OHM, TOLERANCE_J, with_cx};

const DEPTH: f64 = 0.3;
const SIGMA: f64 = 18.0;
const ETA: f64 = 0.6;
const VOLTAGES: [f64; 6] = [20.0, 0.0, 8.0, 0.0, -20.0, 0.0];
const KELVIN: [f64; 3] = [300.0, 310.0, 320.0];
/// Fits resolved at 300 K, every claim valid from 200 K to 500 K.
const SUPPORT: [f64; 3] = [300.0, 200.0, 500.0];

type PlayWinding<'m> = PlanarWinding<'m, HystereticWindingCircuit<'m>, HystereticEddyControls>;

/// The paired squares. Every triangle's fit switches with pinning and
/// viscosity but triangle 6's, which is linear and never dissipates, so the
/// two squares heat differently. With `support`, every claim is bounded in
/// temperature (`[at, low, high]`); without it, the fits have no temperature.
fn problem(support: Option<[f64; 3]>) -> (PlanarHystereticMagnetostatics, Vec<EddyConductor>) {
    let (play, conductor) = sources_at(MU, true, SIGMA, Some(ETA), support);
    let (linear, _) = sources_at(MU, false, SIGMA, None, support);
    let mut materials = vec![play; 8];
    materials[6] = linear;
    problem_from_materials(paired_squares(DEPTH), materials, conductor)
}

/// A voltage source drives the winding through the graph resistor (branch
/// 2). The core conducts in triangles 1 and 7, the copper is stranded in 0
/// and 3 (8 turns per m²) and 4 and 5 (-8), and triangles 2 and 6 are
/// magnetic only.
fn circuit<'m>(
    owner: &'m PlanarHystereticMagnetostatics,
    all: &[EddyConductor],
) -> HystereticWindingCircuit<'m> {
    circuit_on(
        owner,
        all,
        &[1, 7],
        &[(0, 8.0), (3, 8.0), (4, -8.0), (5, -8.0)],
        6,
    )
}

/// The same drive on any mesh: the core triangles, the stranded copper's
/// turns per m² by triangle, and the number of boundary vertices, held at
/// zero potential. Every other triangle is magnetic only.
fn circuit_on<'m>(
    owner: &'m PlanarHystereticMagnetostatics,
    all: &[EddyConductor],
    core: &[usize],
    turns: &[(usize, f64)],
    boundary_vertices: usize,
) -> HystereticWindingCircuit<'m> {
    let cores: Vec<_> = all
        .iter()
        .filter(|s| core.contains(&s.triangle))
        .cloned()
        .collect();
    let mut profile = vec![0.0; all.len()];
    let mut copper = Vec::new();
    for &(triangle, density) in turns {
        profile[triangle] = density;
        copper.push(WindingConductor {
            triangle,
            conductor: all[triangle].conductor.clone(),
        });
    }
    let insulators: Vec<_> = (0..all.len()).filter(|t| !core.contains(t)).collect();
    let initial = owner.demagnetized(0.0, &mut || false).unwrap();
    let field = PlanarHystereticEddy::new(
        &initial,
        cores,
        &insulators,
        &boundary(boundary_vertices, 0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    HystereticWindingCircuit::new_with_end_modes(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 999.0 }),
                (2, 0, Branch::Resistor { ohms: GRAPH_OHM }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profile],
            end_inductance_h: vec![0.02],
        },
        vec![NonlinearWindingPort::CurrentInput],
        vec![EndWindingMode {
            inductance_h: 0.003,
            turns: vec![2.0],
        }],
        copper,
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        NonlinearWindingConfig {
            max_windings: 4,
            max_profile_values: 32,
            current_balance_tolerance_turns: 1e-10,
        },
        &mut || false,
    )
    .unwrap()
}

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

fn component<'m>(circuit: &HystereticWindingCircuit<'m>, tolerance_j: f64) -> PlayWinding<'m> {
    PlanarWinding::new(
        circuit,
        circuit.field().conductors().to_vec(),
        circuit.conductors().to_vec(),
        receivers(),
        terminal_controls(),
        tolerance_j,
    )
    .unwrap()
}

/// The owner's heat records, booked to the three receivers by hand:
/// stranded cells, eddy cells, then every magnetic cell's pinning plus
/// viscous heat, by triangle; the graph resistor to the third.
fn booked(step: &HystereticWindingStep) -> [f64; 3] {
    let side = |triangle: usize| usize::from(triangle >= 4);
    let mut heat = [0.0; 3];
    for cell in &step.cells {
        heat[side(cell.source.triangle)] += cell.heat_j;
    }
    for cell in &step.field.cells {
        heat[side(cell.source.triangle)] += cell.heat_j;
    }
    for cell in &step.field.magnetic.cells {
        heat[side(cell.triangle)] += cell.irreversible_heat_j;
    }
    for resistor in &step.graph_resistors {
        assert_eq!(resistor.branch, 2);
        heat[2] += resistor.heat_j;
    }
    heat
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
                kernel_id: 23,
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
    circuit: HystereticWindingCircuit<'m>,
    volts: &[f64],
) -> PlanarWindingState<'m, HystereticWindingCircuit<'m>> {
    let mut state = PlanarWindingState::new(circuit);
    for &v in volts {
        let mut next = state.clone();
        let mut outputs = [0.0; 4];
        with_cx(|cx| {
            winding.propose(
                &state,
                DT,
                &[v, 300.0, 300.0, 300.0],
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

#[test]
fn a_play_memory_step_is_the_owners_own_and_each_triangle_heats_its_own_receiver() {
    let (owner, all) = problem(Some(SUPPORT));
    let circuit = circuit(&owner, &all);
    let winding = component(&circuit, TOLERANCE_J);
    let layout = PortLayout::new(winding.ports()).unwrap();
    let mut state = PlanarWindingState::new(circuit.clone());
    assert!(state.step().is_none());
    let mut native = circuit;
    let (mut insulating, mut viscous) = (0.0, 0.0);
    for (n, volts) in VOLTAGES.into_iter().enumerate() {
        let inputs = [volts, KELVIN[0], KELVIN[1], KELVIN[2]];
        let mut next = state.clone();
        let mut outputs = [0.0; 4];
        let proposal =
            with_cx(|cx| winding.propose(&state, DT, &inputs, &mut next, &mut outputs, cx))
                .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("the winding shortened the step");
        };
        let expected = native
            .step_to_candidate(
                native.time_s() + DT,
                native.field().conductors().to_vec(),
                native.conductors().to_vec(),
                &[volts],
                terminal_controls(),
                &mut || false,
            )
            .unwrap();
        let step = expected.step();
        assert_eq!(next, PlanarWindingState::new(expected.circuit().clone()));
        assert_eq!(format!("{:?}", next.step()), format!("{:?}", Some(step)));
        assert_eq!(outputs[0].to_bits(), step.native.y[0].to_bits());
        let heat = booked(step);
        for (k, q) in heat.iter().enumerate() {
            assert_eq!(
                outputs[1 + k].to_bits(),
                (-q / DT).to_bits(),
                "receiver {k} at step {n}"
            );
        }
        for (triangle, cell) in step.field.magnetic.cells.iter().enumerate() {
            assert_eq!(cell.triangle, triangle);
            near(
                cell.irreversible_heat_j,
                cell.pinning_heat_j + cell.viscous_heat_j,
            );
            if triangle == 6 {
                assert_eq!(cell.irreversible_heat_j.to_bits(), 0);
            }
        }
        insulating += step.field.magnetic.cells[2].irreversible_heat_j;
        viscous += step.field.magnetic.viscous_heat_j;
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
            "play memory against owner: step {n}, {volts} V, current {:.9} A, heat {heat:?} J, \
             triangle 2 {:.6e} J, viscous {:.6e} J, owner defect {:.1e} J (floor {:.1e} J), \
             balance defect {:.1e} J, endpoint defect {:.1e} J",
            next.circuit().currents_a()[0],
            step.field.magnetic.cells[2].irreversible_heat_j,
            step.field.magnetic.viscous_heat_j,
            receipt.numerical_loss_j - step.numerical_dissipation_j,
            receipt.energy_tolerance_j,
            audit.balance_defect_j,
            audit.endpoint_defect_j
        );
        native = expected.into_parts().0;
        state = next;
    }
    assert!(
        insulating > 0.0 && viscous > 0.0,
        "triangle 2 {insulating} J, viscous {viscous} J"
    );
}

#[test]
fn the_play_memory_component_passes_the_conformance_kit() {
    let (owner, all) = problem(Some(SUPPORT));
    let circuit = circuit(&owner, &all);
    let winding = component(&circuit, TOLERANCE_J);
    let state = driven(&winding, circuit, &[20.0, 0.0]);
    let inputs = [8.0, KELVIN[0], KELVIN[1], KELVIN[2]];
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
    assert!(outputs[1..].iter().all(|rate| *rate < 0.0), "{outputs:?}");
}

#[test]
fn ports_beyond_the_fits_support_refuse_and_refused_steps_leave_the_memory_as_it_was() {
    let (owner, all) = problem(Some(SUPPORT));
    let start = circuit(&owner, &all);
    let winding = component(&start, TOLERANCE_J);
    let state = driven(&winding, start, &[20.0, 0.0]);
    let memory = state.circuit().field().magnetic_state().cell_states();
    assert!(memory[0].polarizations_t()[0][0] != 0.0);
    let attempt = |winding: &PlayWinding<'_>, inputs: [f64; 4], cancel: bool| {
        let mut next = state.clone();
        let mut outputs = [0.0; 4];
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
    let law =
        |result: Result<Proposal, PlanarWindingError>, port: usize, triangle: usize| match result {
            Err(PlanarWindingError::Law {
                receiver,
                element: LawElement::Field(at),
                kelvin,
                message,
            }) if receiver == port && at == triangle => {
                println!("field triangle {at} at {kelvin} K refused: {message}");
                message
            }
            other => panic!("expected triangle {triangle}'s fit to refuse, got {other:?}"),
        };
    // The support is inclusive on both sides.
    attempt(&winding, [8.0, 500.0, 200.0, 300.0], false).unwrap();
    law(attempt(&winding, [8.0, 520.0, 300.0, 300.0], false), 0, 0);
    law(attempt(&winding, [8.0, 300.0, 150.0, 300.0], false), 1, 4);
    assert!(matches!(
        attempt(&winding, [8.0, 300.0, 300.0, 300.0], true),
        Err(PlanarWindingError::Cancelled)
    ));
    // An allowance at the rounding floor alone refuses the owner's own
    // closure defect after the step was solved.
    let strict = component(state.circuit(), f64::MIN_POSITIVE);
    assert!(matches!(
        attempt(&strict, [8.0, 300.0, 300.0, 300.0], false),
        Err(PlanarWindingError::Balance { .. })
    ));

    // A state whose field has another mesh: one square, four triangles.
    let (play, conductor) = sources_at(MU, true, SIGMA, Some(ETA), Some(SUPPORT));
    let (small, small_all) = problem_from_materials(square(DEPTH), vec![play; 4], conductor);
    let foreign = PlanarWindingState::new(circuit_on(
        &small,
        &small_all,
        &[],
        &[(0, 8.0), (2, -8.0)],
        4,
    ));
    let mut next = foreign.clone();
    let mut outputs = [0.0; 4];
    match with_cx(|cx| {
        winding.propose(
            &foreign,
            DT,
            &[8.0, 300.0, 300.0, 300.0],
            &mut next,
            &mut outputs,
            cx,
        )
    }) {
        Err(PlanarWindingError::Invalid { what }) => {
            assert_eq!(what, "the state's circuit is not this Component's topology");
        }
        other => panic!("expected the topology refusal, got {other:?}"),
    }

    // Fits with no temperature of their own hold at none.
    let (bare, bare_all) = problem(None);
    let unbounded = circuit(&bare, &bare_all);
    let winding = component(&unbounded, TOLERANCE_J);
    let state = PlanarWindingState::new(unbounded);
    let mut next = state.clone();
    let mut outputs = [0.0; 4];
    let message = law(
        with_cx(|cx| {
            winding.propose(
                &state,
                DT,
                &[8.0, 300.0, 300.0, 300.0],
                &mut next,
                &mut outputs,
                cx,
            )
        }),
        0,
        0,
    );
    assert_eq!(
        message,
        "the magnetic fit was resolved at no single absolute temperature"
    );
    assert_eq!(next, state);
}

#[test]
fn every_magnetic_triangle_maps_to_one_receiver_and_the_lane_has_its_own_word() {
    let (owner, all) = problem(Some(SUPPORT));
    let circuit = circuit(&owner, &all);
    let build = |receivers: Vec<HeatReceiver>| {
        PlanarWinding::new(
            &circuit,
            circuit.field().conductors().to_vec(),
            circuit.conductors().to_vec(),
            receivers,
            terminal_controls(),
            TOLERANCE_J,
        )
    };
    let refused = |result: Result<PlayWinding<'_>, PlanarWindingError>, what: &str| match result {
        Err(PlanarWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    // Triangles 2 and 6 carry no conductor; their fits still dissipate.
    let mut without = receivers();
    without[0].triangles.retain(|&triangle| triangle != 2);
    refused(
        build(without),
        "a dissipating element is mapped to no receiver",
    );
    let mut twice = receivers();
    twice[1].triangles.push(2);
    refused(
        build(twice),
        "a dissipating element is mapped to two receivers",
    );
    let mut beyond = receivers();
    beyond[1].triangles.push(8);
    refused(
        build(beyond),
        "a receiver names a triangle that carries no conductor",
    );

    let words = |winding: &PlayWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let base = words(&component(&circuit, TOLERANCE_J));
    assert_eq!(base, words(&component(&circuit, TOLERANCE_J)));
    // The tolerance, then the lane's word.
    let lane = u64::from_le_bytes(*b"vecplay1");
    assert_eq!(base[base.len() - 2..], [TOLERANCE_J.to_bits(), lane]);
    let card = super::common::source(true, 18.0);
    let saturating = super::common::problem(&card, 0.3);
    let mut other = Vec::new();
    super::component(&super::series(&saturating, &card), &card).identity_words(&mut other);
    assert_eq!(other.last(), Some(&TOLERANCE_J.to_bits()));
    let mut swapped = receivers();
    swapped[0].triangles[2] = 6;
    swapped[1].triangles[2] = 2;
    assert_ne!(base, words(&build(swapped).unwrap()));
}

/// The circuit's own conductor in each triangle, at every temperature.
#[derive(Debug)]
struct Own(BTreeMap<usize, ResolvedConductorState>);

impl ConductorLaw for Own {
    fn conductor(&self, triangle: usize, _: f64) -> Result<ResolvedConductorState, String> {
        self.0
            .get(&triangle)
            .cloned()
            .ok_or_else(|| format!("no conductor in triangle {triangle}"))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.0.len() as u64);
    }
}

#[test]
fn laws_that_answer_the_owners_conductors_step_as_fixed_conductors_do() {
    let (owner, all) = problem(Some(SUPPORT));
    let circuit = circuit(&owner, &all);
    let fixed = component(&circuit, TOLERANCE_J);
    let own: Arc<dyn ConductorLaw> = Arc::new(Own(circuit
        .field()
        .conductors()
        .iter()
        .map(|c| (c.triangle, c.conductor.clone()))
        .chain(
            circuit
                .conductors()
                .iter()
                .map(|c| (c.triangle, c.conductor.clone())),
        )
        .collect()));
    let following: PlayWinding<'_> = PlanarWinding::at_port_temperatures(
        &circuit,
        vec![own.clone(), own.clone(), own],
        BTreeMap::new(),
        receivers(),
        terminal_controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let (mut a, mut b) = (
        PlanarWindingState::new(circuit.clone()),
        PlanarWindingState::new(circuit),
    );
    for volts in [20.0, 0.0, -8.0] {
        let inputs = [volts, KELVIN[0], KELVIN[1], KELVIN[2]];
        let (mut next_a, mut next_b) = (a.clone(), b.clone());
        let (mut out_a, mut out_b) = ([0.0; 4], [0.0; 4]);
        with_cx(|cx| fixed.propose(&a, DT, &inputs, &mut next_a, &mut out_a, cx)).unwrap();
        with_cx(|cx| following.propose(&b, DT, &inputs, &mut next_b, &mut out_b, cx)).unwrap();
        assert_eq!(next_a, next_b);
        assert_eq!(out_a.map(f64::to_bits), out_b.map(f64::to_bits));
        (a, b) = (next_a, next_b);
    }
}

#[test]
fn the_play_memory_owner_has_no_magnetic_lane() {
    let (owner, all) = problem(Some(SUPPORT));
    let circuit = circuit(&owner, &all);
    let laws = super::magnetic::laws(&super::magnetic::card(0.002));
    let refused = PlayWinding::magnetic_at_port_temperatures(
        &circuit,
        vec![],
        BTreeMap::new(),
        laws,
        receivers(),
        terminal_controls(),
        TOLERANCE_J,
    );
    let Err(PlanarWindingError::Invalid { what }) = refused else {
        panic!("{refused:?}");
    };
    assert_eq!(
        what,
        "this owner refreshes no field materials with temperature"
    );
}
