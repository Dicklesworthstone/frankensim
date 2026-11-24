//! The planar nonlinear winding owner as one Component: its ports, a step
//! that is the owner's own step bit for bit, heat that leaves through the
//! ports as the owner books it, the kernel's audit, conformance, refusals
//! and identity; then the same with conductor and resistor laws at the
//! temperatures the heat ports receive. With `magnetic-hysteresis`, the
//! play-memory owner in `winding_component/hysteretic.rs`.
#![cfg(feature = "component")]

#[path = "nonlinear_winding_circuit/common.rs"]
mod common;
#[cfg(feature = "magnetic-hysteresis")]
#[path = "winding_component/hysteretic.rs"]
mod hysteretic;

use std::collections::BTreeMap;
use std::sync::Arc;

use common::{conductors, controls, profiles};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_em::planar::eddy::{EddyConductor, nonlinear::PlanarNonlinearEddy};
use fs_em::planar::nonlinear::PlanarNonlinearMagnetostatics;
use fs_em::winding::component::{
    ConductorLaw, HeatReceiver, LawElement, PlanarWinding, PlanarWindingError, PlanarWindingState,
    ResistorLaw,
};
use fs_em::winding::nonlinear::{
    NonlinearWindingCandidate, NonlinearWindingCircuit, NonlinearWindingInitial,
    NonlinearWindingStep, NonlinearWindings, resistors::GraphResistorMaterial,
};
use fs_em::winding::ohmic::WindingConductor;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::MaterialCard;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductor, ResolvedConductorState};
use fs_material::field_storage::{
    nonlinear::REVERSIBLE_FIELD_STRENGTH_PROPERTY, permanent::FIELD_STRENGTH_DIMS,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{
    Carrier, Class, Component, EndOfStep, PortLayout, PortRole, PortSpec, Proposal, audit_step,
};
use fs_qty::{Density, Dims};

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
    series_on(
        problem,
        card,
        core(card, 300.0),
        Branch::Capacitor { farads: 0.05 },
    )
}

/// [`series`] with these core conductors and this element across the graph
/// resistor.
fn series_on<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
    core: Vec<EddyConductor>,
    across: Branch,
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
        core,
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
                (2, 0, across),
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
            charges_c: if matches!(across, Branch::Capacitor { .. }) {
                vec![0.0]
            } else {
                vec![]
            },
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

/// The owner's heat records, booked to the three receivers by hand.
fn booked(step: &NonlinearWindingStep) -> [f64; 3] {
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
    heat
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
        let heat = booked(step);
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

/// The fixture's conductor at `kelvin`, as `common::material::conductor`
/// resolves it, but refused rather than unwrapped outside the card's
/// validity.
fn resolve(card: &MaterialCard, kelvin: f64) -> Result<ResolvedConductorState, String> {
    let basis = ElectricalPropertyBasis::Resistivity;
    let scalar = |name: &str, dims: Dims, admissibility| {
        ScalarPropertyRequirement::try_new(name, dims, admissibility).map_err(|e| format!("{e:?}"))
    };
    let requirements = [
        basis.requirement().map_err(|e| format!("{e:?}"))?,
        scalar(
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            ScalarAdmissibility::NonNegative,
        )?,
        scalar(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        scalar(
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        scalar(
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    let point = common::material::point(kelvin);
    let bundle = resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|e| format!("{e:?}"))?;
    ResolvedConductorState::try_from_material_with_basis(bundle, basis)
        .map_err(|e| format!("{e:?}"))
}

/// The card's identity and resistivity at 300 K, as law words.
fn card_words(card: &MaterialCard) -> Vec<u64> {
    let state = resolve(card, 300.0).unwrap();
    let mut words = vec![state.resistivity_ohm_m().to_bits()];
    for chunk in state.material().card_identity().0.chunks(8) {
        words.push(u64::from_le_bytes(chunk.try_into().unwrap()));
    }
    words
}

/// The length of a unit-area path of the card with `GRAPH_OHM` at 300 K.
fn path_length(card: &MaterialCard) -> f64 {
    GRAPH_OHM / resolve(card, 300.0).unwrap().resistivity_ohm_m()
}

/// Every triangle's conductor from one card, at the port's temperature or,
/// when frozen, at that fixed temperature.
#[derive(Debug)]
struct CardLaw {
    card: MaterialCard,
    frozen: Option<f64>,
}

impl ConductorLaw for CardLaw {
    fn conductor(&self, _triangle: usize, kelvin: f64) -> Result<ResolvedConductorState, String> {
        resolve(&self.card, self.frozen.unwrap_or(kelvin))
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend(card_words(&self.card));
        words.push(self.frozen.map_or(0, f64::to_bits));
    }
}

/// A uniform unit-area path of one card.
#[derive(Debug)]
struct PathLaw {
    card: MaterialCard,
    length_m: f64,
    frozen: Option<f64>,
}

impl ResistorLaw for PathLaw {
    fn conductor(&self, kelvin: f64) -> Result<ResolvedConductor, String> {
        resolve(&self.card, self.frozen.unwrap_or(kelvin))?
            .into_uniform_conductor(self.length_m, 1.0)
            .map_err(|e| format!("{e:?}"))
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend(card_words(&self.card));
        words.extend([self.length_m.to_bits(), self.frozen.map_or(0, f64::to_bits)]);
    }
}

/// Conductor laws by receiver, and resistor laws by branch.
type Laws = (
    Vec<Arc<dyn ConductorLaw>>,
    BTreeMap<usize, Arc<dyn ResistorLaw>>,
);

/// One card law per receiver and a path law for the graph resistor.
fn laws(card: &MaterialCard, frozen: Option<f64>) -> Laws {
    let cells: Arc<dyn ConductorLaw> = Arc::new(CardLaw {
        card: card.clone(),
        frozen,
    });
    let path: Arc<dyn ResistorLaw> = Arc::new(PathLaw {
        card: card.clone(),
        length_m: path_length(card),
        frozen,
    });
    (
        vec![cells.clone(), cells.clone(), cells],
        BTreeMap::from([(2, path)]),
    )
}

fn following<'m>(circuit: &NonlinearWindingCircuit<'m>, card: &MaterialCard) -> PlanarWinding<'m> {
    let (cells, paths) = laws(card, None);
    PlanarWinding::at_port_temperatures(circuit, cells, paths, receivers(), controls(), TOLERANCE_J)
        .unwrap()
}

/// The core and stranded conductors, each resolved by hand at its
/// receiver's temperature.
fn conductors_at(
    card: &MaterialCard,
    kelvin: [f64; 3],
) -> (Vec<EddyConductor>, Vec<WindingConductor>) {
    let at = |triangle: usize| resolve(card, kelvin[usize::from(triangle >= 4)]).unwrap();
    let core = [2, 3, 6, 7]
        .into_iter()
        .map(|triangle| EddyConductor {
            triangle,
            region: if triangle < 4 { 7 } else { 8 },
            conductor: at(triangle),
        })
        .collect();
    let windings = [0, 1, 4, 5]
        .into_iter()
        .map(|triangle| WindingConductor {
            triangle,
            conductor: at(triangle),
        })
        .collect();
    (core, windings)
}

/// The owner's own step from `owner`, with every element resolved by hand at
/// its receiver's temperature, and the graph resistor's resistance.
fn owner_at<'m>(
    owner: &NonlinearWindingCircuit<'m>,
    card: &MaterialCard,
    kelvin: [f64; 3],
) -> (NonlinearWindingCandidate<'m>, f64) {
    let (core, windings) = conductors_at(card, kelvin);
    let path = resolve(card, kelvin[2])
        .unwrap()
        .into_uniform_conductor(path_length(card), 1.0)
        .unwrap();
    let ohms = path.resistance_ohm();
    let driven = owner
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: 2,
                conductor: path,
            }],
            &mut || false,
        )
        .unwrap();
    let candidate = driven
        .step_to_candidate(
            owner.time_s() + DT,
            core,
            windings,
            &[VOLTS],
            controls(),
            &mut || false,
        )
        .unwrap();
    (candidate, ohms)
}

#[test]
fn laws_give_each_element_at_its_ports_temperature_as_the_owner_would() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    // The fallible resolution is the fixture's own inside the card's validity,
    // and the conductivity moves by more than 5 % across 300-400 K.
    assert_eq!(
        resolve(&card, 300.0).unwrap(),
        common::material::conductor(&card, 300.0, ElectricalPropertyBasis::Resistivity)
    );
    let ratio = resolve(&card, 300.0).unwrap().resistivity_ohm_m()
        / resolve(&card, 400.0).unwrap().resistivity_ohm_m();
    println!("laws: conductivity at 400 K is {ratio} of its value at 300 K");
    assert!(ratio < 0.95);
    let winding = following(&circuit, &card);
    let layout = PortLayout::new(winding.ports()).unwrap();
    let mut state = PlanarWindingState::new(circuit.clone());
    let mut owner = circuit;
    let steps = [
        [300.0, 350.0, 400.0],
        [330.0, 390.0, 310.0],
        [400.0, 300.0, 360.0],
    ];
    for (n, kelvin) in steps.into_iter().enumerate() {
        let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
        let mut next = state.clone();
        let mut outputs = [0.0; 4];
        let proposal =
            with_cx(|cx| winding.propose(&state, DT, &inputs, &mut next, &mut outputs, cx))
                .unwrap();
        let Proposal::Step(receipt) = proposal else {
            panic!("the winding shortened the step");
        };
        let (expected, ohms) = owner_at(&owner, &card, kelvin);
        let step = expected.step();
        assert_eq!(next, PlanarWindingState::new(expected.circuit().clone()));
        assert_eq!(
            next.circuit().graph().branches[2].2,
            Branch::Resistor { ohms }
        );
        assert_eq!(outputs[0].to_bits(), step.native.y[0].to_bits());
        let heat = booked(step);
        for (k, q) in heat.iter().enumerate() {
            assert_eq!(outputs[1 + k].to_bits(), (-q / DT).to_bits());
        }
        assert_eq!(
            receipt.stored_change_j.to_bits(),
            step.stored_change_j.to_bits()
        );
        let entropy: f64 = heat.iter().zip(kelvin).map(|(q, t)| q / t).sum();
        assert_eq!(
            receipt.entropy_generated_j_per_k.to_bits(),
            entropy.to_bits()
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
            "laws against owner: step {}, ports at {kelvin:?} K, current {:.9} A, graph \
             resistor {ohms:.9} ohm, heat {heat:?} J, balance defect {:.1e} J, endpoint \
             defect {:.1e} J",
            n + 1,
            next.circuit().currents_a()[0],
            audit.balance_defect_j,
            audit.endpoint_defect_j
        );
        owner = expected.into_parts().0;
        state = next;
    }
}

#[test]
fn a_law_that_ignores_temperature_steps_as_fixed_conductors_do() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let (cells, paths) = laws(&card, Some(300.0));
    let frozen = PlanarWinding::at_port_temperatures(
        &circuit,
        cells,
        paths,
        receivers(),
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    // The fixed Component steps the circuit whose graph resistor already
    // carries the law's resistance.
    let path = resolve(&card, 300.0)
        .unwrap()
        .into_uniform_conductor(path_length(&card), 1.0)
        .unwrap();
    let bound = circuit
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: 2,
                conductor: path,
            }],
            &mut || false,
        )
        .unwrap();
    let fixed = component(&bound, &card);
    let mut a = PlanarWindingState::new(circuit);
    let mut b = PlanarWindingState::new(bound);
    for kelvin in [[300.0, 350.0, 400.0], [400.0, 330.0, 310.0]] {
        let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
        let (mut na, mut nb) = (a.clone(), b.clone());
        let (mut oa, mut ob) = ([0.0; 4], [0.0; 4]);
        with_cx(|cx| frozen.propose(&a, DT, &inputs, &mut na, &mut oa, cx)).unwrap();
        with_cx(|cx| fixed.propose(&b, DT, &inputs, &mut nb, &mut ob, cx)).unwrap();
        assert_eq!(na, nb);
        assert_eq!(oa.map(f64::to_bits), ob.map(f64::to_bits));
        (a, b) = (na, nb);
    }
}

#[test]
fn with_no_resistor_law_the_graph_resistor_keeps_its_authored_value() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let (cells, paths) = laws(&card, None);
    // The owner rebinds every graph resistor or none.
    let parallel = series_on(
        &problem,
        &card,
        core(&card, 300.0),
        Branch::Resistor { ohms: 1.0 },
    );
    let mut both = receivers();
    both[2].resistors.push(3);
    match PlanarWinding::at_port_temperatures(
        &parallel,
        cells.clone(),
        paths.clone(),
        both,
        controls(),
        TOLERANCE_J,
    ) {
        Err(PlanarWindingError::Invalid { what }) => {
            assert_eq!(what, "every graph resistor or none has a law");
        }
        other => panic!("expected the partial-law refusal, got {other:?}"),
    }
    let winding = PlanarWinding::at_port_temperatures(
        &circuit,
        cells,
        BTreeMap::new(),
        receivers(),
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let kelvin = [300.0, 350.0, 400.0];
    let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
    let state = PlanarWindingState::new(circuit.clone());
    let mut next = state.clone();
    let mut outputs = [0.0; 4];
    with_cx(|cx| winding.propose(&state, DT, &inputs, &mut next, &mut outputs, cx)).unwrap();
    let (core, windings) = conductors_at(&card, kelvin);
    let expected = circuit
        .step_to_candidate(DT, core, windings, &[VOLTS], controls(), &mut || false)
        .unwrap();
    assert_eq!(next, PlanarWindingState::new(expected.circuit().clone()));
    assert_eq!(next.circuit().graph(), circuit.graph());
    // The authored value is then the Component's own.
    let hot = circuit
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: 2,
                conductor: paths[&2].conductor(400.0).unwrap(),
            }],
            &mut || false,
        )
        .unwrap();
    let hot = PlanarWindingState::new(hot);
    let mut ignored = hot.clone();
    match with_cx(|cx| winding.propose(&hot, DT, &inputs, &mut ignored, &mut outputs, cx)) {
        Err(PlanarWindingError::Invalid { what }) => {
            assert_eq!(what, "the state's circuit is not this Component's topology");
        }
        other => panic!("expected the topology refusal, got {other:?}"),
    }
}

#[test]
fn the_component_with_laws_passes_the_conformance_kit() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let winding = following(&circuit, &card);
    let state = PlanarWindingState::new(circuit);
    let inputs = [VOLTS, 300.0, 350.0, 400.0];
    let mut stale = state.clone();
    let mut outputs = [0.0; 4];
    with_cx(|cx| winding.propose(&state, DT, &inputs, &mut stale, &mut outputs, cx)).unwrap();
    let case = Case {
        state: &state,
        stale: &stale,
        dt_s: DT,
        inputs: &inputs,
    };
    with_cx(|cx| check_component(&winding, &case, cx)).expect("component conformance");
}

#[test]
fn laws_refuse_what_they_cannot_answer_and_constructions_they_cannot_drive() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let (cells, paths) = laws(&card, None);
    let build = |cells: Vec<Arc<dyn ConductorLaw>>,
                 paths: BTreeMap<usize, Arc<dyn ResistorLaw>>| {
        PlanarWinding::at_port_temperatures(
            &circuit,
            cells,
            paths,
            receivers(),
            controls(),
            TOLERANCE_J,
        )
    };
    let refused = |result: Result<PlanarWinding<'_>, PlanarWindingError>, what: &str| match result {
        Err(PlanarWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        other => panic!("expected the refusal {what:?}, got {other:?}"),
    };
    refused(
        build(cells[..2].to_vec(), paths.clone()),
        "one conductor law per heat receiver required",
    );
    let mut inductor = paths.clone();
    inductor.insert(1, paths[&2].clone());
    refused(
        build(cells.clone(), inductor),
        "a resistor law names a branch that is not a resistor",
    );

    // The card's validity ends at 500 K. A law asked above it refuses with
    // the element, its receiver and the temperature; nothing is clamped.
    let winding = following(&circuit, &card);
    let state = PlanarWindingState::new(circuit.clone());
    let step = |winding: &PlanarWinding<'_>, inputs: &[f64], from: &PlanarWindingState<'_>| {
        let mut next = from.clone();
        let mut outputs = vec![0.0; inputs.len()];
        with_cx(|cx| winding.propose(from, DT, inputs, &mut next, &mut outputs, cx))
    };
    match step(&winding, &[VOLTS, 300.0, 520.0, 310.0], &state) {
        Err(PlanarWindingError::Law {
            receiver: 1,
            element: LawElement::Triangle(4),
            kelvin,
            ..
        }) if kelvin.to_bits() == 520.0_f64.to_bits() => {}
        other => panic!("expected the stranded law's refusal, got {other:?}"),
    }
    match step(&winding, &[VOLTS, 300.0, 310.0, 520.0], &state) {
        Err(PlanarWindingError::Law {
            receiver: 2,
            element: LawElement::Resistor(2),
            kelvin,
            ..
        }) if kelvin.to_bits() == 520.0_f64.to_bits() => {}
        other => panic!("expected the resistor law's refusal, got {other:?}"),
    }
    // A law of another card meets the owner's own continuity refusal.
    let (other, _) = laws(&common::source(true, 9.0), None);
    let foreign = build(other, paths.clone()).unwrap();
    match step(&foreign, &[VOLTS, 300.0, 310.0, 320.0], &state) {
        Err(PlanarWindingError::Owner(error)) => println!("laws: another card: {error}"),
        other => panic!("expected the owner's refusal, got {other:?}"),
    }
    // Resistor values may differ, the graph's shape may not.
    let parallel = series_on(
        &problem,
        &card,
        core(&card, 300.0),
        Branch::Resistor { ohms: 1.0 },
    );
    match step(
        &winding,
        &[VOLTS, 300.0, 310.0, 320.0],
        &PlanarWindingState::new(parallel),
    ) {
        Err(PlanarWindingError::Invalid { what }) => {
            assert_eq!(what, "the state's circuit is not this Component's topology");
        }
        other => panic!("expected the topology refusal, got {other:?}"),
    }
    // Fixed conductors fix the resistor values too.
    let hot = circuit
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: 2,
                conductor: paths[&2].conductor(400.0).unwrap(),
            }],
            &mut || false,
        )
        .unwrap();
    match step(
        &component(&circuit, &card),
        &[VOLTS, 300.0, 310.0, 320.0],
        &PlanarWindingState::new(hot),
    ) {
        Err(PlanarWindingError::Invalid { what }) => {
            assert_eq!(what, "the state's circuit is not this Component's topology");
        }
        other => panic!("expected the topology refusal, got {other:?}"),
    }
    assert_eq!(state, PlanarWindingState::new(circuit));
}

#[test]
fn identity_words_tell_laws_and_cells_apart() {
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = series(&problem, &card);
    let words = |winding: &PlanarWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let at = |cells, paths| {
        PlanarWinding::at_port_temperatures(
            &circuit,
            cells,
            paths,
            receivers(),
            controls(),
            TOLERANCE_J,
        )
        .unwrap()
    };
    let base = words(&following(&circuit, &card));
    assert_eq!(base, words(&following(&circuit, &card)));
    assert_ne!(base, words(&component(&circuit, &card)));
    let (cells, paths) = laws(&card, Some(300.0));
    assert_ne!(base, words(&at(cells, paths)));
    let other = common::source(true, 9.0);
    let (cells, paths) = laws(&card, None);
    let (foreign_cells, foreign_paths) = laws(&other, None);
    assert_ne!(base, words(&at(foreign_cells, paths.clone())));
    assert_ne!(base, words(&at(cells, foreign_paths)));
    // The same graph with its core regions labelled 9 and 10.
    let relabelled = core(&card, 300.0)
        .into_iter()
        .map(|c| EddyConductor {
            region: c.region + 2,
            ..c
        })
        .collect();
    let moved = series_on(
        &problem,
        &card,
        relabelled,
        Branch::Capacitor { farads: 0.05 },
    );
    assert_ne!(base, words(&following(&moved, &card)));
}
