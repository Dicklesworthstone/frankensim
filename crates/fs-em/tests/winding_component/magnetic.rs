//! The magnetic lane: each triangle's permeability follows the temperature
//! its receiver's port gives (fs-material's affine `ThermalMagneticMaterial`),
//! the state carries the temperatures its materials were evaluated at, and
//! the Component declares its stored entropy. Synthetic constitutive data
//! test the discrete accounts; they are not a measured core.

#[path = "../winding_ohmic/common.rs"]
mod source;

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_em::planar::eddy::nonlinear::NonlinearEddyMagneticCell;
use fs_em::planar::nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics};
use fs_em::winding::component::{
    ConductorLaw, LawElement, MagneticLaw, MagneticLaws, PlanarWinding, PlanarWindingError,
    PlanarWindingState,
};
use fs_em::winding::nonlinear::{
    NonlinearWindingCandidate, NonlinearWindingCircuit, NonlinearWindingInitial, NonlinearWindings,
};
use fs_em::winding::ohmic::WindingConductor;
use fs_matdb::MaterialCard;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_material::field_storage::thermal::magnetic::ThermalMagneticMaterial;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, PortLayout, Proposal, Receipt, audit_step};

use super::{DT, GRAPH_OHM, HeatReceiver, TOLERANCE_J, VOLTS, booked, common, receivers, with_cx};

const DEPTH: f64 = 0.3;
/// The permeability's slope \[H/(m K)\]: 1 H/m at 300 K, 1.2 H/m at 400 K.
const SLOPE: f64 = 0.002;
const START_K: f64 = 300.0;
const RANGE_K: [f64; 2] = [250.0, 400.0];
/// The source's volts, and every port at 300 K.
const HELD: [f64; 4] = [VOLTS, START_K, START_K, START_K];
const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;

/// A conductor with a resistivity of 1/180 Ohm m at 300 K, rising 0.4 % per
/// K, and a permeability affine in absolute temperature with `slope`.
pub(super) fn card(slope: f64) -> MaterialCard {
    source::absolute_affine_magnetic_card(1.0 / 180.0, 0.004, 1.0, slope, BASIS)
}

fn material(card: &MaterialCard, kelvin: f64) -> ResolvedMaterialStatePoint {
    source::bundle_at(card, &source::absolute_point(kelvin), BASIS)
}

fn conductor(card: &MaterialCard, kelvin: f64) -> Result<ResolvedConductorState, String> {
    ResolvedConductorState::try_from_material_with_basis(material(card, kelvin), BASIS)
        .map_err(|e| format!("{e:?}"))
}

/// The card's law over `RANGE_K`, pinned at 300 K.
fn law(card: &MaterialCard) -> ThermalMagneticMaterial {
    ThermalMagneticMaterial::new(card, &material(card, START_K), "T", RANGE_K).unwrap()
}

/// The card's law in every triangle.
pub(super) fn laws(card: &MaterialCard) -> MagneticLaws {
    let law = law(card);
    MagneticLaws {
        triangles: (0..8)
            .map(|triangle| MagneticLaw {
                triangle,
                law: law.clone(),
            })
            .collect(),
        assembly: source::limits(),
    }
}

/// Every triangle linear, as the card gives it at 300 K.
fn problem(card: &MaterialCard) -> PlanarNonlinearMagnetostatics {
    PlanarNonlinearMagnetostatics::new(
        common::mesh(DEPTH),
        vec![IsotropicMagneticMaterial::Linear(material(card, START_K)); 8],
        source::limits(),
        &mut || false,
    )
    .unwrap()
}

/// The stranded copper of [`circuit`] at its receivers' temperatures.
fn windings(card: &MaterialCard, kelvin: [f64; 3]) -> Vec<WindingConductor> {
    [0, 1, 4, 5]
        .into_iter()
        .map(|triangle| WindingConductor {
            triangle,
            conductor: conductor(card, kelvin[usize::from(triangle >= 4)]).unwrap(),
        })
        .collect()
}

/// A voltage source drives the winding, carrying `amperes`, through the
/// graph resistor (branch 2). The copper is stranded in triangles 0 and 1
/// (12 turns per m²) and 4 and 5 (-12); triangles 2, 3, 6 and 7 are
/// magnetic only.
fn circuit<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    card: &MaterialCard,
    amperes: f64,
) -> NonlinearWindingCircuit<'m> {
    let profiles = vec![common::profiles(2).remove(0)];
    NonlinearWindingCircuit::new(
        common::initial_field(problem, &profiles, &[amperes]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: GRAPH_OHM }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        windings(card, [START_K; 3]),
        NonlinearWindingInitial {
            currents_a: vec![amperes],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![VOLTS],
        },
        common::config(),
        &mut || false,
    )
    .unwrap()
}

/// Each receiver's copper from the card at the port's temperature.
#[derive(Debug)]
struct Copper(MaterialCard);

impl ConductorLaw for Copper {
    fn conductor(&self, _triangle: usize, kelvin: f64) -> Result<ResolvedConductorState, String> {
        conductor(&self.0, kelvin)
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        for chunk in material(&self.0, START_K).card_identity().0.chunks(8) {
            words.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
    }
}

fn copper(card: &MaterialCard) -> Vec<Arc<dyn ConductorLaw>> {
    let law: Arc<dyn ConductorLaw> = Arc::new(Copper(card.clone()));
    vec![law.clone(), law.clone(), law]
}

fn lane_of<'m>(
    circuit: &NonlinearWindingCircuit<'m>,
    card: &MaterialCard,
    laws: MagneticLaws,
    receivers: Vec<HeatReceiver>,
) -> Result<PlanarWinding<'m>, PlanarWindingError> {
    PlanarWinding::magnetic_at_port_temperatures(
        circuit,
        copper(card),
        BTreeMap::new(),
        laws,
        receivers,
        common::controls(),
        TOLERANCE_J,
    )
}

fn lane<'m>(circuit: &NonlinearWindingCircuit<'m>, card: &MaterialCard) -> PlanarWinding<'m> {
    lane_of(circuit, card, laws(card), receivers()).unwrap()
}

/// The owner's own step from `owner` with every element and every
/// triangle's material resolved by hand at its receiver's temperature.
fn owner_at<'m>(
    owner: &NonlinearWindingCircuit<'m>,
    card: &MaterialCard,
    kelvin: [f64; 3],
) -> NonlinearWindingCandidate<'m> {
    let law = law(card);
    let materials = (0..8)
        .map(|triangle| {
            let at = law
                .at_temperature(kelvin[usize::from(triangle >= 4)])
                .unwrap();
            IsotropicMagneticMaterial::Linear(at.material().clone())
        })
        .collect();
    owner
        .step_with_materials_candidate(
            owner.time_s() + DT,
            materials,
            vec![],
            windings(card, kelvin),
            &[VOLTS],
            source::limits(),
            common::controls(),
            &mut || false,
        )
        .unwrap()
}

/// Each receiver's `s = -dpsi/dT` at one end of a step, from the owner's
/// cell records and the card's law, by hand.
fn entropy(
    card: &MaterialCard,
    cells: &[NonlinearEddyMagneticCell],
    end: usize,
    kelvin: [f64; 3],
) -> [f64; 3] {
    let law = law(card);
    let mut s = [0.0; 3];
    for cell in cells {
        let receiver = usize::from(cell.triangle >= 4);
        let [x, y] = cell.flux_density_t[end];
        s[receiver] += law
            .at_temperature(kelvin[receiver])
            .unwrap()
            .state(x.hypot(y), cell.volume_m3)
            .unwrap()
            .entropy_j_k;
    }
    s
}

type Proposed<'m> = (PlanarWindingState<'m>, [f64; 4], Receipt);

fn proposed<'m>(
    winding: &PlanarWinding<'m>,
    state: &PlanarWindingState<'m>,
    inputs: &[f64],
) -> Result<Proposed<'m>, PlanarWindingError> {
    let mut new = state.clone();
    let mut outputs = [0.0; 4];
    match with_cx(|cx| winding.propose(state, DT, inputs, &mut new, &mut outputs, cx))? {
        Proposal::Step(receipt) => Ok((new, outputs, receipt)),
        Proposal::Shorten { dt_s } => panic!("the winding shortened the step to {dt_s} s"),
    }
}

/// A proposal from `state`, passed by the kernel's `audit_step` with the
/// Component's stored energies.
fn audited<'m>(
    winding: &PlanarWinding<'m>,
    state: &PlanarWindingState<'m>,
    inputs: &[f64],
) -> Proposed<'m> {
    let (new, outputs, receipt) = proposed(winding, state, inputs).unwrap();
    let layout = PortLayout::new(winding.ports()).unwrap();
    let stored = |state: &PlanarWindingState<'m>| winding.stored_energy_j(state).unwrap();
    audit_step(
        winding.ports(),
        &layout,
        DT,
        inputs,
        &outputs,
        &receipt,
        stored(state),
        stored(&new),
    )
    .expect("audited");
    (new, outputs, receipt)
}

/// The conformance kit's case from `state`, with a stale candidate slot;
/// its energy and entropy audits run inside.
fn kit(winding: &PlanarWinding<'_>, state: &PlanarWindingState<'_>, inputs: &[f64]) -> Receipt {
    let (stale, _, _) = proposed(winding, state, inputs).unwrap();
    let case = Case {
        state,
        stale: &stale,
        dt_s: DT,
        inputs,
    };
    match with_cx(|cx| check_component(winding, &case, cx)).expect("component conformance") {
        Observed::Stepped { receipt, .. } => receipt,
        Observed::Shortened { dt_s } => panic!("the winding shortened the step to {dt_s} s"),
    }
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn invalid<T>(result: Result<T, PlanarWindingError>, what: &str) {
    match result {
        Err(PlanarWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
        Err(other) => panic!("expected the refusal {what:?}, got {other:?}"),
        Ok(_) => panic!("expected the refusal {what:?}"),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} is not {expected:e} within {tolerance:e}"
    );
}

#[test]
fn at_held_temperatures_a_step_is_the_owners_own_and_its_heat_is_joule_less_t_ds() {
    let card = card(SLOPE);
    let problem = problem(&card);
    let circuit = circuit(&problem, &card, 1.0);
    let winding = lane(&circuit, &card);
    let kelvin = [START_K; 3];
    let state = PlanarWindingState::at_temperatures(circuit.clone(), kelvin.to_vec());
    let inputs = HELD;
    let (new, outputs, receipt) = audited(&winding, &state, &inputs);
    let owner = owner_at(&circuit, &card, kelvin);
    let step = owner.step();
    assert_eq!(
        bits(new.circuit().currents_a()),
        bits(owner.circuit().currents_a())
    );
    assert_eq!(
        bits(new.circuit().field().potential_wb_m()),
        bits(owner.circuit().field().potential_wb_m())
    );
    assert_eq!(new.temperatures(), Some(&kelvin[..]));
    // No change of temperature: no parameter work, and nothing is booked
    // beyond the owner's own numerical dissipation.
    assert_eq!(step.parameter_work_j.to_bits(), 0.0_f64.to_bits());
    let cells = &step.field.magnetic_cells;
    assert_eq!(cells.len(), 8);
    assert!(cells.iter().all(|cell| cell.parameter_work_j == 0.0));
    close(
        receipt.numerical_loss_j,
        step.numerical_dissipation_j,
        1e-14,
    );
    // Each port's heat is its Joule heat less T dS of its triangles; the
    // field rises, so dS > 0 and the reversible part flows in.
    let joule = booked(step);
    let (s0, s1) = (
        entropy(&card, cells, 0, kelvin),
        entropy(&card, cells, 1, kelvin),
    );
    for receiver in 0..3 {
        let net = joule[receiver] - kelvin[receiver] * (s1[receiver] - s0[receiver]);
        close(outputs[1 + receiver], -net / DT, 1e-12);
    }
    assert!(s1[0] > s0[0] && s1[1] > s0[1], "{s0:?} -> {s1:?}");
    assert_eq!(s0[2].to_bits(), 0.0_f64.to_bits());
    close(receipt.dissipated_j, joule.iter().sum(), 1e-14);
    close(
        receipt.entropy_generated_j_per_k,
        joule.iter().sum::<f64>() / START_K,
        1e-15,
    );
}

#[test]
fn ramped_temperatures_pass_the_energy_audit_and_book_a_second_order_parameter_term() {
    let card = card(SLOPE);
    let problem = problem(&card);
    let circuit = circuit(&problem, &card, 1.0);
    let winding = lane(&circuit, &card);
    let start = PlanarWindingState::at_temperatures(circuit.clone(), vec![START_K; 3]);
    // Ten steps, the halves rising 2 K per step: `audit_step` closes each.
    let mut state = start.clone();
    for k in 1..=10 {
        let t = 2.0_f64.mul_add(f64::from(k), START_K);
        let inputs = [VOLTS, t, t, START_K];
        state = audited(&winding, &state, &inputs).0;
        assert_eq!(state.temperatures(), Some(&inputs[1..]));
    }
    // One step from the start, the rise halved each time. The term booked
    // as numerical loss beyond the owner's dissipation, -(P + dT s0), is
    // measured from the receipt and the owner's own step.
    let mut terms = Vec::new();
    for rise in [16.0, 8.0, 4.0, 2.0] {
        let kelvin = [START_K + rise, START_K + rise, START_K];
        let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
        let (_, _, receipt) = proposed(&winding, &start, &inputs).unwrap();
        let owner = owner_at(&circuit, &card, kelvin);
        let step = owner.step();
        let measured = receipt.numerical_loss_j - step.numerical_dissipation_j;
        let s0 = entropy(&card, &step.field.magnetic_cells, 0, [START_K; 3]);
        let parameter: f64 = step
            .field
            .magnetic_cells
            .iter()
            .map(|cell| cell.parameter_work_j)
            .sum();
        let by_hand = -(parameter + rise * (s0[0] + s0[1]));
        println!("rise {rise} K: booked {measured:e} J, by hand {by_hand:e} J");
        close(measured, by_hand, 1e-12);
        // psi is convex in T, so P > -s0 dT and the term is negative.
        assert!(measured < 0.0, "{measured:e}");
        terms.push(measured);
    }
    for pair in terms.windows(2) {
        let ratio = pair[0] / pair[1];
        println!("halving the rise divides the term by {ratio}");
        assert!((3.6..4.4).contains(&ratio), "{ratio}");
    }
}

#[test]
fn a_rising_field_cools_its_receivers_while_joule_heat_and_generation_stay_positive() {
    let card = card(SLOPE);
    let problem = problem(&card);
    let circuit = circuit(&problem, &card, 0.0);
    let winding = lane(&circuit, &card);
    let state = PlanarWindingState::at_temperatures(circuit, vec![START_K; 3]);
    let inputs = HELD;
    let (_, outputs, receipt) = audited(&winding, &state, &inputs);
    println!("from rest, heat rates in [W]: {:?}", &outputs[1..]);
    // Heat flows into the winding from both halves; the graph resistor's
    // receiver still takes its Joule heat.
    assert!(outputs[1] > 0.0 && outputs[2] > 0.0, "{outputs:?}");
    assert!(outputs[3] < 0.0, "{outputs:?}");
    assert!(receipt.dissipated_j > 0.0);
    assert!(receipt.entropy_generated_j_per_k > 0.0);
}

#[test]
fn the_lane_refuses_what_it_cannot_book_and_never_clamps_a_temperature() {
    let card = card(SLOPE);
    let problem = problem(&card);
    let circuit = circuit(&problem, &card, 1.0);
    let winding = lane(&circuit, &card);
    let state = PlanarWindingState::at_temperatures(circuit.clone(), vec![START_K; 3]);
    // A port beyond the law's interval refuses; the temperature is never
    // clamped.
    let beyond = [VOLTS, 450.0, START_K, START_K];
    match proposed(&winding, &state, &beyond) {
        Err(PlanarWindingError::Law {
            receiver: 0,
            element: LawElement::Permeability(0),
            kelvin,
            ..
        }) if kelvin.to_bits() == 450.0_f64.to_bits() => {}
        other => panic!("{other:?}"),
    }
    // A state whose temperatures are not the field's: the law's free energy
    // at them is not the field's.
    let hotter = PlanarWindingState::at_temperatures(circuit.clone(), vec![340.0, 340.0, START_K]);
    match proposed(&winding, &hotter, &[VOLTS, 340.0, 340.0, START_K]) {
        Err(PlanarWindingError::Law {
            receiver: 0,
            element: LawElement::Permeability(0),
            kelvin,
            message,
        }) => {
            assert_eq!(kelvin.to_bits(), 340.0_f64.to_bits());
            assert!(message.contains("is not the field's"), "{message}");
        }
        other => panic!("{other:?}"),
    }
    // A state without temperatures, or with the wrong number of them.
    let lane_only = "a state carries temperatures in the magnetic lane only, one positive \
                     temperature per receiver";
    let bare = PlanarWindingState::new(circuit.clone());
    invalid(proposed(&winding, &bare, &HELD), lane_only);
    let short = PlanarWindingState::at_temperatures(circuit.clone(), vec![START_K; 2]);
    invalid(proposed(&winding, &short, &HELD), lane_only);
    // Outside the lane, a state with temperatures refuses.
    let mut conducting = receivers();
    for receiver in &mut conducting {
        receiver.triangles.retain(|t| [0, 1, 4, 5].contains(t));
    }
    let plain = PlanarWinding::at_port_temperatures(
        &circuit,
        copper(&card),
        BTreeMap::new(),
        conducting.clone(),
        common::controls(),
        TOLERANCE_J,
    )
    .unwrap();
    invalid(proposed(&plain, &state, &HELD), lane_only);
    // Only law triangles may be receivers' triangles without a conductor.
    invalid(
        PlanarWinding::at_port_temperatures(
            &circuit,
            copper(&card),
            BTreeMap::new(),
            receivers(),
            common::controls(),
            TOLERANCE_J,
        ),
        "a receiver names a triangle that carries no conductor",
    );
    // Every law triangle maps to a receiver.
    invalid(
        lane_of(&circuit, &card, laws(&card), conducting),
        "a dissipating element is mapped to no receiver",
    );
    // Laws of another card, two laws in one triangle, a triangle beyond the
    // mesh, and no law.
    invalid(
        lane_of(&circuit, &card, laws(&self::card(0.001)), receivers()),
        "a magnetic law's source is not its triangle's field material",
    );
    let with = |edit: fn(&mut Vec<MagneticLaw>)| {
        let mut laws = laws(&card);
        edit(&mut laws.triangles);
        lane_of(&circuit, &card, laws, receivers())
    };
    invalid(
        with(|laws| laws[1].triangle = 0),
        "two magnetic laws name one triangle",
    );
    invalid(
        with(|laws| laws[7].triangle = 8),
        "a magnetic law names a triangle whose field material is not linear",
    );
    invalid(with(Vec::clear), "the magnetic lane has no law");
}

#[test]
fn proposals_build_their_own_problems_and_states_compare_them_by_content() {
    let card = card(SLOPE);
    let problem = problem(&card);
    let circuit = circuit(&problem, &card, 1.0);
    let winding = lane(&circuit, &card);
    let state = PlanarWindingState::at_temperatures(circuit.clone(), vec![START_K; 3]);
    let inputs = [VOLTS, 310.0, 320.0, START_K];
    let (a, _, _) = proposed(&winding, &state, &inputs).unwrap();
    let (b, _, _) = proposed(&winding, &state, &inputs).unwrap();
    assert!(!core::ptr::eq(
        a.circuit().field().problem(),
        b.circuit().field().problem()
    ));
    assert_eq!(a, b);
    let cooler =
        PlanarWindingState::at_temperatures(a.circuit().clone(), vec![310.0, 319.0, START_K]);
    assert_ne!(a, cooler);
    assert_ne!(a, PlanarWindingState::new(a.circuit().clone()));
    assert_ne!(state, PlanarWindingState::new(circuit));
}

#[test]
fn the_lane_has_words_of_its_own() {
    let card = card(SLOPE);
    let problem = problem(&card);
    let circuit = circuit(&problem, &card, 1.0);
    let words = |winding: &PlanarWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let lane = words(&lane(&circuit, &card));
    assert_eq!(lane, words(&self::lane(&circuit, &card)));
    let mut narrower = laws(&card);
    narrower.triangles[3].law =
        ThermalMagneticMaterial::new(&card, &material(&card, START_K), "T", [250.0, 390.0])
            .unwrap();
    let mut bounded = laws(&card);
    bounded.assembly.max_triplets -= 1;
    // Triangle 0 conducts, so it may go without a law.
    let mut fewer = laws(&card);
    fewer.triangles.remove(0);
    let mut conducting = receivers();
    conducting[0].triangles.retain(|&t| t != 3);
    conducting[1].triangles.push(3);
    for (laws, receivers) in [
        (narrower, receivers()),
        (bounded, receivers()),
        (fewer, receivers()),
        (self::laws(&card), conducting),
    ] {
        assert_ne!(
            lane,
            words(&lane_of(&circuit, &card, laws, receivers).unwrap())
        );
    }
}

#[test]
fn the_lane_declares_the_stored_entropy_of_its_state_and_other_lanes_none() {
    let card = card(SLOPE);
    let problem = problem(&card);
    let circuit = circuit(&problem, &card, 1.0);
    let winding = lane(&circuit, &card);
    let kelvin = [310.0, 320.0, START_K];
    let state = PlanarWindingState::at_temperatures(circuit.clone(), vec![START_K; 3]);
    let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
    let (new, _, _) = proposed(&winding, &state, &inputs).unwrap();
    // Measured from each state alone, it is the sum of `s` by hand from the
    // owner's own step and the card's law, at the state's temperatures.
    let owner = owner_at(&circuit, &card, kelvin);
    let cells = &owner.step().field.magnetic_cells;
    let measured = |state: &PlanarWindingState<'_>| {
        winding
            .stored_entropy_j_per_k(state)
            .unwrap()
            .expect("declared")
    };
    let s0 = entropy(&card, cells, 0, [START_K; 3]);
    let s1 = entropy(&card, cells, 1, kelvin);
    close(measured(&state), s0.iter().sum(), 1e-15);
    close(measured(&new), s1.iter().sum(), 1e-15);
    invalid(
        winding.stored_entropy_j_per_k(&PlanarWindingState::new(circuit)),
        "a state of the magnetic lane carries its temperatures",
    );
    // The saturating lanes declare no stored entropy.
    let card = common::source(true, 18.0);
    let problem = common::problem(&card, 0.3);
    let circuit = super::series(&problem, &card);
    let state = PlanarWindingState::new(circuit.clone());
    for winding in [
        super::component(&circuit, &card),
        super::following(&circuit, &card),
    ] {
        assert_eq!(winding.stored_entropy_j_per_k(&state).unwrap(), None);
    }
}

#[test]
fn the_lane_passes_the_conformance_kit_and_its_entropy_audit() {
    let card = card(SLOPE);
    let problem = problem(&card);
    // Held, rising 2 K per step, and from rest: ten steps each, the kit's
    // energy and entropy audits closing at every one.
    for (amperes, rise) in [(1.0, 0.0_f64), (1.0, 2.0), (0.0, 0.0)] {
        let circuit = circuit(&problem, &card, amperes);
        let winding = lane(&circuit, &card);
        let mut state = PlanarWindingState::at_temperatures(circuit, vec![START_K; 3]);
        for k in 1..=10 {
            let t = rise.mul_add(f64::from(k), START_K);
            let inputs = [VOLTS, t, t, START_K];
            let receipt = kit(&winding, &state, &inputs);
            assert!(receipt.entropy_generated_j_per_k > 0.0);
            state = proposed(&winding, &state, &inputs).unwrap().0;
        }
    }
}
