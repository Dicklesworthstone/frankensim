//! The magnetic lane on the 3D owner: each law cell's permeability follows
//! the temperature its receiver's port gives (fs-material's affine
//! `ThermalMagneticMaterial`), the state carries the temperatures its
//! materials were evaluated at, and the Component declares its stored
//! entropy. Every cell of the cube is linear and has a law; the strands cross
//! 24 cells, core conductors fill the other 12 of the lower half, and the
//! other 12 of the upper half are magnetic only. Receiver 0 takes the
//! conducting cells, receiver 1 the magnetic-only cells and receiver 2 the
//! graph resistor. Synthetic constitutive data test the discrete accounts;
//! they are not a measured core.

#[path = "../winding_ohmic/common.rs"]
mod source;

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_em::spatial::nonlinear::SpatialNonlinearMaterial;
use fs_em::spatial::vector_potential::induction::nonlinear::SpatialNonlinearInductionMagneticCell;
use fs_em::spatial::vector_potential::induction::nonlinear::winding::{
    SpatialNonlinearWindingCandidate, SpatialNonlinearWindingCircuit, SpatialNonlinearWindingStep,
    SpatialWindingConductor,
};
use fs_em::spatial::vector_potential::nonlinear::SpatialNonlinearVectorMagnetostatics;
use fs_em::winding::component::spatial::{
    SpatialConductorLaw, SpatialHeatReceiver, SpatialLawElement, SpatialMagneticLaw,
    SpatialWinding, SpatialWindingError, SpatialWindingState,
};
use fs_matdb::MaterialCard;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_material::field_storage::thermal::magnetic::ThermalMagneticMaterial;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_phs::circuit::Branch;
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, PortLayout, Proposal, Receipt, audit_step};

use super::{DT, FILL, LIMITS, TOLERANCE_J, VOLTS, cells, controls, cube_of, series_of, with_cx};

/// The permeability's slope \[H/(m K)\]: 1 H/m at 300 K, 1.2 H/m at 400 K.
const SLOPE: f64 = 0.002;
const START_K: f64 = 300.0;
const RANGE_K: [f64; 2] = [250.0, 400.0];
/// The source's volts, and every port at 300 K.
const HELD: [f64; 4] = [VOLTS, START_K, START_K, START_K];
const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;
const CELLS: usize = 48;

/// A conductor with a resistivity of 1/18 Ohm m at 300 K, rising 0.4 % per
/// K, and a permeability affine in absolute temperature with `slope`.
fn card(slope: f64) -> MaterialCard {
    source::absolute_affine_magnetic_card(1.0 / 18.0, 0.004, 1.0, slope, BASIS)
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

/// The card's law in every cell.
fn laws(card: &MaterialCard) -> Vec<SpatialMagneticLaw> {
    let law = law(card);
    (0..CELLS)
        .map(|cell| SpatialMagneticLaw {
            cell,
            law: law.clone(),
        })
        .collect()
}

/// The cube with every cell linear, as the card gives it at 300 K.
fn problem(card: &MaterialCard) -> SpatialNonlinearVectorMagnetostatics {
    cube_of(SpatialNonlinearMaterial::Linear(material(card, START_K)))
}

/// Whether a cell carries no conductor: the strands do not cross it and it
/// lies in the upper half (cells 24 to 47).
fn magnetic_only(crossed: &[usize], cell: usize) -> bool {
    cell >= CELLS / 2 && !crossed.contains(&cell)
}

/// A cell's receiver: 1 for a magnetic-only cell, otherwise 0.
fn receiver(crossed: &[usize], cell: usize) -> usize {
    usize::from(magnetic_only(crossed, cell))
}

/// Core conductors at `kelvin` in the cells of the lower half (0 to 23) that
/// the strands do not cross.
fn core(
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    kelvin: f64,
) -> Vec<Option<ResolvedConductorState>> {
    let conductor = conductor(card, kelvin).unwrap();
    let mut core = vec![None; CELLS];
    for cell in cells(p, false) {
        if cell < CELLS / 2 {
            core[cell] = Some(conductor.clone());
        }
    }
    core
}

/// The strands at `kelvin` in every cell they cross.
fn strands(
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    kelvin: f64,
) -> Vec<SpatialWindingConductor> {
    let conductor = conductor(card, kelvin).unwrap();
    cells(p, true)
        .into_iter()
        .map(|cell| SpatialWindingConductor {
            winding: 0,
            cell,
            fill_fraction: FILL,
            conductor: conductor.clone(),
        })
        .collect()
}

/// The parent's series circuit at rest, its conductors at 300 K.
fn circuit<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
) -> SpatialNonlinearWindingCircuit<'m> {
    series_of(
        p,
        core(p, card, START_K),
        strands(p, card, START_K),
        Branch::Capacitor { farads: 0.05 },
    )
}

/// Each receiver's conductors from the card at the port's temperature.
#[derive(Debug)]
struct Copper(MaterialCard);

impl SpatialConductorLaw<ResolvedConductorState> for Copper {
    fn core(&self, _cell: usize, kelvin: f64) -> Result<ResolvedConductorState, String> {
        conductor(&self.0, kelvin)
    }
    fn strand(
        &self,
        _winding: usize,
        _cell: usize,
        kelvin: f64,
    ) -> Result<ResolvedConductorState, String> {
        conductor(&self.0, kelvin)
    }
    fn identity_words(&self, words: &mut Vec<u64>) {
        for chunk in material(&self.0, START_K).card_identity().0.chunks(8) {
            words.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
    }
}

fn copper(card: &MaterialCard) -> Vec<Arc<dyn SpatialConductorLaw<ResolvedConductorState>>> {
    let law: Arc<dyn SpatialConductorLaw<ResolvedConductorState>> = Arc::new(Copper(card.clone()));
    vec![law.clone(), law.clone(), law]
}

/// The conducting cells; the magnetic-only cells; the graph resistor.
fn receivers(p: &SpatialNonlinearVectorMagnetostatics) -> Vec<SpatialHeatReceiver> {
    let crossed = cells(p, true);
    let (only, conducting) = (0..CELLS).partition(|&cell| magnetic_only(&crossed, cell));
    vec![
        SpatialHeatReceiver {
            cells: conducting,
            resistors: vec![],
        },
        SpatialHeatReceiver {
            cells: only,
            resistors: vec![],
        },
        SpatialHeatReceiver {
            cells: vec![],
            resistors: vec![2],
        },
    ]
}

/// Receivers for a Component outside the lane: every cell they name
/// conducts. The strands' cells; the core cells; the graph resistor.
fn conducting(p: &SpatialNonlinearVectorMagnetostatics) -> Vec<SpatialHeatReceiver> {
    let mut receivers = super::receivers(p);
    receivers[1].cells.retain(|&cell| cell < CELLS / 2);
    receivers
}

/// The owner's heat records, booked to the lane's receivers by hand: the
/// strands' and the core's to receiver 0, the graph resistor's to 2.
fn booked(step: &SpatialNonlinearWindingStep) -> [f64; 3] {
    let mut heat = [0.0; 3];
    for cell in &step.cells {
        heat[0] += cell.heat_j;
    }
    for cell in &step.field.cells {
        heat[0] += cell.heat_j;
    }
    for resistor in &step.graph_resistors {
        heat[2] += resistor.heat_j;
    }
    heat
}

fn lane_of<'m>(
    circuit: &SpatialNonlinearWindingCircuit<'m>,
    card: &MaterialCard,
    laws: Vec<SpatialMagneticLaw>,
    receivers: Vec<SpatialHeatReceiver>,
) -> Result<SpatialWinding<'m>, SpatialWindingError> {
    SpatialWinding::magnetic_at_port_temperatures(
        circuit,
        copper(card),
        BTreeMap::new(),
        laws,
        receivers,
        controls(),
        TOLERANCE_J,
    )
}

fn lane<'m>(
    circuit: &SpatialNonlinearWindingCircuit<'m>,
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
) -> SpatialWinding<'m> {
    lane_of(circuit, card, laws(card), receivers(p)).unwrap()
}

/// The owner's own step from `owner` with every conductor and every cell's
/// material resolved by hand at its receiver's temperature.
fn owner_at<'m>(
    owner: &SpatialNonlinearWindingCircuit<'m>,
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    kelvin: [f64; 3],
) -> SpatialNonlinearWindingCandidate<'m> {
    let law = law(card);
    let crossed = cells(p, true);
    let materials = (0..CELLS)
        .map(|cell| {
            let at = law
                .at_temperature(kelvin[receiver(&crossed, cell)])
                .unwrap();
            SpatialNonlinearMaterial::Linear(at.material().clone())
        })
        .collect();
    owner
        .step_with_materials_candidate(
            owner.time_s() + DT,
            &[[0.0; 3]; 27],
            materials,
            core(p, card, kelvin[0]),
            strands(p, card, kelvin[0]),
            &[VOLTS],
            controls(),
            &mut || false,
        )
        .unwrap()
}

/// Each receiver's `s = -dpsi/dT` at one end of a step, from the owner's
/// cell records and the card's law, by hand.
fn entropy(
    p: &SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    records: &[SpatialNonlinearInductionMagneticCell],
    end: usize,
    kelvin: [f64; 3],
) -> [f64; 3] {
    let law = law(card);
    let crossed = cells(p, true);
    let mut s = [0.0; 3];
    for record in records {
        let receiver = receiver(&crossed, record.cell);
        let [x, y, z] = record.flux_density_t[end];
        s[receiver] += law
            .at_temperature(kelvin[receiver])
            .unwrap()
            .state(x.hypot(y).hypot(z), record.volume_m3)
            .unwrap()
            .entropy_j_k;
    }
    s
}

type Proposed<'m> = (SpatialWindingState<'m>, [f64; 4], Receipt);

fn proposed<'m>(
    winding: &SpatialWinding<'m>,
    state: &SpatialWindingState<'m>,
    inputs: &[f64],
) -> Result<Proposed<'m>, SpatialWindingError> {
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
    winding: &SpatialWinding<'m>,
    state: &SpatialWindingState<'m>,
    inputs: &[f64],
) -> Proposed<'m> {
    let (new, outputs, receipt) = proposed(winding, state, inputs).unwrap();
    let layout = PortLayout::new(winding.ports()).unwrap();
    let stored = |state: &SpatialWindingState<'m>| winding.stored_energy_j(state).unwrap();
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

/// The lane's state after two held steps from rest at 300 K: the field has
/// risen, so every law cell stores entropy.
fn charged<'m>(
    winding: &SpatialWinding<'m>,
    circuit: &SpatialNonlinearWindingCircuit<'m>,
) -> SpatialWindingState<'m> {
    let mut state =
        SpatialWindingState::at_temperatures(circuit.clone(), vec![START_K; 3]).unwrap();
    for _ in 0..2 {
        state = audited(winding, &state, &HELD).0;
    }
    state
}

/// The conformance kit's case from `state`, with a stale candidate slot;
/// its energy and entropy audits run inside.
fn kit(winding: &SpatialWinding<'_>, state: &SpatialWindingState<'_>, inputs: &[f64]) -> Receipt {
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

fn invalid<T>(result: Result<T, SpatialWindingError>, what: &str) {
    match result {
        Err(SpatialWindingError::Invalid { what: actual }) => assert_eq!(actual, what),
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
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let winding = lane(&circuit, &p, &card);
    let kelvin = [START_K; 3];
    let state = charged(&winding, &circuit);
    assert_eq!(state.temperatures(), Some(&kelvin[..]));
    let (new, outputs, receipt) = audited(&winding, &state, &HELD);
    let owner = owner_at(state.circuit(), &p, &card, kelvin);
    let step = owner.step();
    let rebuilt =
        SpatialWindingState::at_temperatures(owner.circuit().clone(), kelvin.to_vec()).unwrap();
    assert_eq!(new, rebuilt);
    let kept = new.step().unwrap();
    assert_eq!(kept.stored_after_j.to_bits(), step.stored_after_j.to_bits());
    assert_eq!(outputs[0].to_bits(), step.native.y[0].to_bits());
    // No change of temperature: no parameter work, and nothing is booked
    // beyond the owner's own numerical dissipation and the residual of its
    // balance, here within the rounding floor.
    assert_eq!(step.parameter_work_j.to_bits(), 0.0_f64.to_bits());
    let records = &step.field.magnetic_cells;
    assert_eq!(records.len(), CELLS);
    assert!(records.iter().all(|record| record.parameter_work_j == 0.0));
    let gap = receipt.numerical_loss_j - step.numerical_dissipation_j;
    println!(
        "held: numerical loss {:e} J, the owner's {:e} J, gap {gap:e} J, floor {:e} J",
        receipt.numerical_loss_j, step.numerical_dissipation_j, receipt.energy_tolerance_j
    );
    close(
        receipt.numerical_loss_j,
        step.numerical_dissipation_j,
        receipt.energy_tolerance_j,
    );
    // Each port's heat is its Joule heat less T dS of its cells; the field
    // rises, so dS > 0 and the reversible part flows in.
    let joule = booked(step);
    let (s0, s1) = (
        entropy(&p, &card, records, 0, kelvin),
        entropy(&p, &card, records, 1, kelvin),
    );
    println!("held: Joule {joule:?} J, s {s0:?} -> {s1:?} J/K, rates {outputs:?}");
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
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let winding = lane(&circuit, &p, &card);
    let start = charged(&winding, &circuit);
    // Five steps, both cell receivers rising 2 K per step: `audit_step`
    // closes each.
    let mut state = start.clone();
    for k in 1..=5 {
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
        let owner = owner_at(start.circuit(), &p, &card, kelvin);
        let step = owner.step();
        let measured = receipt.numerical_loss_j - step.numerical_dissipation_j;
        let records = &step.field.magnetic_cells;
        let s0 = entropy(&p, &card, records, 0, [START_K; 3]);
        let parameter: f64 = records.iter().map(|record| record.parameter_work_j).sum();
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
fn a_rising_field_cools_the_cells_without_conductors_while_generation_stays_positive() {
    let card = card(SLOPE);
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let winding = lane(&circuit, &p, &card);
    let state = SpatialWindingState::at_temperatures(circuit, vec![START_K; 3]).unwrap();
    let (_, outputs, receipt) = audited(&winding, &state, &HELD);
    println!("from rest, heat rates in [W]: {:?}", &outputs[1..]);
    // Heat flows into the magnetic-only cells, which have no Joule heat;
    // the conducting cells and the graph resistor still give theirs out.
    assert!(outputs[2] > 0.0, "{outputs:?}");
    assert!(outputs[1] < 0.0 && outputs[3] < 0.0, "{outputs:?}");
    assert!(receipt.dissipated_j > 0.0);
    assert!(receipt.entropy_generated_j_per_k > 0.0);
}

#[test]
fn the_lane_refuses_what_it_cannot_book_and_never_clamps_a_temperature() {
    let card = card(SLOPE);
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let winding = lane(&circuit, &p, &card);
    let state = charged(&winding, &circuit);
    // A port beyond the law's interval refuses; the temperature is never
    // clamped. Cell 0 is a core cell.
    let beyond = [VOLTS, 450.0, START_K, START_K];
    match proposed(&winding, &state, &beyond) {
        Err(SpatialWindingError::Law {
            receiver: 0,
            element: SpatialLawElement::Permeability(0),
            kelvin,
            ..
        }) if kelvin.to_bits() == 450.0_f64.to_bits() => {}
        other => panic!("{other:?}"),
    }
    // A state whose temperatures are not the field's: the law's free energy
    // at them is not the field's.
    let hotter =
        SpatialWindingState::at_temperatures(state.circuit().clone(), vec![340.0, 340.0, START_K])
            .unwrap();
    match proposed(&winding, &hotter, &[VOLTS, 340.0, 340.0, START_K]) {
        Err(SpatialWindingError::Law {
            receiver: 0,
            element: SpatialLawElement::Permeability(0),
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
    let bare = SpatialWindingState::new(circuit.clone()).unwrap();
    invalid(proposed(&winding, &bare, &HELD), lane_only);
    let short = SpatialWindingState::at_temperatures(circuit.clone(), vec![START_K; 2]).unwrap();
    invalid(proposed(&winding, &short, &HELD), lane_only);
    // Outside the lane, a state with temperatures refuses.
    let plain = SpatialWinding::at_port_temperatures(
        &circuit,
        copper(&card),
        BTreeMap::new(),
        conducting(&p),
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    invalid(proposed(&plain, &state, &HELD), lane_only);
    // Only law cells may be receivers' cells without a conductor.
    invalid(
        SpatialWinding::at_port_temperatures(
            &circuit,
            copper(&card),
            BTreeMap::new(),
            receivers(&p),
            controls(),
            TOLERANCE_J,
        ),
        "a receiver names a cell that carries no conductor",
    );
    // Every law cell maps to a receiver.
    invalid(
        lane_of(&circuit, &card, laws(&card), conducting(&p)),
        "a dissipating element is mapped to no receiver",
    );
    // Laws of another card, two laws in one cell, a cell beyond the mesh,
    // and no law.
    invalid(
        lane_of(&circuit, &card, laws(&self::card(0.001)), receivers(&p)),
        "a magnetic law's source is not its cell's field material",
    );
    let with = |edit: fn(&mut Vec<SpatialMagneticLaw>)| {
        let mut laws = laws(&card);
        edit(&mut laws);
        lane_of(&circuit, &card, laws, receivers(&p))
    };
    invalid(
        with(|laws| laws[1].cell = 0),
        "two magnetic laws name one cell",
    );
    let not_linear = "a magnetic law names a cell whose field material is not linear";
    invalid(with(|laws| laws[47].cell = CELLS), not_linear);
    invalid(with(Vec::clear), "the magnetic lane has no law");
    // A law in a cell of the saturating cube, whose law is reversible.
    let saturating = super::card(18.0);
    let q = super::cube(&saturating);
    let (conductors, paths) = super::laws(&saturating, None);
    invalid(
        SpatialWinding::magnetic_at_port_temperatures(
            &super::series(&q, &saturating),
            conductors,
            paths,
            laws(&card).split_off(47),
            super::receivers(&q),
            controls(),
            TOLERANCE_J,
        ),
        not_linear,
    );
}

#[test]
fn proposals_build_their_own_problems_and_states_compare_them_by_content() {
    let card = card(SLOPE);
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let winding = lane(&circuit, &p, &card);
    let state = charged(&winding, &circuit);
    let inputs = [VOLTS, 310.0, 320.0, START_K];
    let (a, _, _) = proposed(&winding, &state, &inputs).unwrap();
    let (b, _, _) = proposed(&winding, &state, &inputs).unwrap();
    assert!(!core::ptr::eq(
        a.circuit().field().problem(),
        b.circuit().field().problem()
    ));
    assert_eq!(a, b);
    let cooler =
        SpatialWindingState::at_temperatures(a.circuit().clone(), vec![310.0, 319.0, START_K])
            .unwrap();
    assert_ne!(a, cooler);
    assert_ne!(a, SpatialWindingState::new(a.circuit().clone()).unwrap());
    assert_ne!(
        state,
        SpatialWindingState::new(state.circuit().clone()).unwrap()
    );
}

/// `p` rebuilt with `edit` applied to its positions, tetrahedra and
/// materials.
fn rebuilt(
    p: &SpatialNonlinearVectorMagnetostatics,
    edit: impl FnOnce(&mut [[f64; 3]], &mut [[u32; 4]], &mut [SpatialNonlinearMaterial]),
) -> SpatialNonlinearVectorMagnetostatics {
    let mut positions = p.positions_m().to_vec();
    let mut tets = p.mesh().tets.clone();
    let mut materials = p.materials().to_vec();
    edit(&mut positions, &mut tets, &mut materials);
    SpatialNonlinearVectorMagnetostatics::new(positions, tets, materials, LIMITS, &mut || false)
        .unwrap()
}

/// The lane's state of [`circuit`] on `p`, at rest at 300 K.
fn at_rest<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
) -> SpatialWindingState<'m> {
    SpatialWindingState::at_temperatures(circuit(p, card), vec![START_K; 3]).unwrap()
}

#[test]
fn states_of_the_lane_compare_their_problems_by_mesh_and_materials() {
    let card = card(SLOPE);
    let p = problem(&card);
    let rest = at_rest(&p, &card);
    // Another allocation of the same problem.
    let same = rebuilt(&p, |_, _, _| {});
    assert_eq!(rest, at_rest(&same, &card));
    // At rest the potential, currents and linkages are zero, and none of
    // these edits touches a conducting cell (moving a strand's vertex would
    // change the winding's resistance), so only the problem tells them
    // apart: a magnetic-only cell's material at 310 K, two magnetic-only
    // cells in each other's place, and the corner vertex (0, 0, 1), which
    // only two magnetic-only cells touch, moved in its plane.
    let only = receivers(&p)[1].cells.clone();
    let (x, y) = (only[0], only[1]);
    let warmer = rebuilt(&p, |_, _, materials| {
        materials[x] = SpatialNonlinearMaterial::Linear(material(&card, 310.0));
    });
    let swapped = rebuilt(&p, |_, tets, _| tets.swap(x, y));
    let corner = p
        .positions_m()
        .iter()
        .position(|x| x.map(f64::to_bits) == [0.0, 0.0, 1.0].map(f64::to_bits))
        .unwrap();
    let touching: Vec<usize> = (0..CELLS)
        .filter(|&cell| p.mesh().tets[cell].contains(&(corner as u32)))
        .collect();
    assert!(
        touching.iter().all(|cell| only.contains(cell)),
        "{touching:?}"
    );
    let moved = rebuilt(&p, |positions, _, _| {
        positions[corner] = [-0.01, -0.01, 1.0];
    });
    for other in [&warmer, &swapped, &moved] {
        assert_ne!(rest, at_rest(other, &card));
    }
}

#[test]
fn the_lane_has_words_of_its_own() {
    let card = card(SLOPE);
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let words = |winding: &SpatialWinding<'_>| {
        let mut words = Vec::new();
        winding.identity_words(&mut words);
        words
    };
    let lane = words(&lane(&circuit, &p, &card));
    assert_eq!(lane, words(&self::lane(&circuit, &p, &card)));
    let mut narrower = laws(&card);
    narrower[30].law =
        ThermalMagneticMaterial::new(&card, &material(&card, START_K), "T", [250.0, 390.0])
            .unwrap();
    // Cell 0 has a core conductor, so it may go without a law.
    let mut fewer = laws(&card);
    fewer.remove(0);
    // Move a magnetic-only cell to the conducting cells' receiver.
    let mut regrouped = receivers(&p);
    let moved = regrouped[1].cells.remove(0);
    regrouped[0].cells.push(moved);
    for (laws, receivers) in [
        (narrower, receivers(&p)),
        (fewer, receivers(&p)),
        (self::laws(&card), regrouped),
    ] {
        assert_ne!(
            lane,
            words(&lane_of(&circuit, &card, laws, receivers).unwrap())
        );
    }
    // Outside the lane the words are the Component's with conductor laws.
    let plain = |receivers| {
        SpatialWinding::at_port_temperatures(
            &circuit,
            copper(&card),
            BTreeMap::new(),
            receivers,
            controls(),
            TOLERANCE_J,
        )
        .unwrap()
    };
    assert_ne!(lane, words(&plain(conducting(&p))));
}

#[test]
fn the_lane_declares_the_stored_entropy_of_its_state_and_other_lanes_none() {
    let card = card(SLOPE);
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let winding = lane(&circuit, &p, &card);
    let kelvin = [310.0, 320.0, START_K];
    let state = charged(&winding, &circuit);
    let inputs = [VOLTS, kelvin[0], kelvin[1], kelvin[2]];
    let (new, _, _) = proposed(&winding, &state, &inputs).unwrap();
    // Measured from each state alone, it is the sum of `s` by hand from the
    // owner's own step and the card's law, at the state's temperatures; the
    // stored energy adds `T s` to the energy the state carries.
    let owner = owner_at(state.circuit(), &p, &card, kelvin);
    let records = &owner.step().field.magnetic_cells;
    let measured = |state: &SpatialWindingState<'_>| {
        winding
            .stored_entropy_j_per_k(state)
            .unwrap()
            .expect("declared")
    };
    let s0 = entropy(&p, &card, records, 0, [START_K; 3]);
    let s1 = entropy(&p, &card, records, 1, kelvin);
    close(measured(&state), s0.iter().sum(), 1e-15);
    close(measured(&new), s1.iter().sum(), 1e-15);
    let ts: f64 = s1.iter().zip(kelvin).map(|(s, t)| t * s).sum();
    close(
        winding.stored_energy_j(&new).unwrap() - new.stored_energy_j(),
        ts,
        1e-15,
    );
    let carries = "a state of the magnetic lane carries its temperatures";
    let short = SpatialWindingState::at_temperatures(circuit.clone(), vec![START_K; 2]).unwrap();
    invalid(winding.stored_entropy_j_per_k(&short), carries);
    let bare = SpatialWindingState::new(circuit).unwrap();
    invalid(winding.stored_entropy_j_per_k(&bare), carries);
    invalid(winding.stored_energy_j(&bare), carries);
    // The fixed and conductor-law Components declare no stored entropy.
    let card = super::card(18.0);
    let p = super::cube(&card);
    let circuit = super::series(&p, &card);
    let state = super::state(circuit.clone());
    for winding in [
        super::component(&circuit, &p, &card),
        super::following(&circuit, &p, &card),
    ] {
        assert_eq!(winding.stored_entropy_j_per_k(&state).unwrap(), None);
    }
}

#[test]
fn the_lane_passes_the_conformance_kit_and_its_entropy_audit() {
    let card = card(SLOPE);
    let p = problem(&card);
    let circuit = circuit(&p, &card);
    let winding = lane(&circuit, &p, &card);
    // Charged and held, charged and rising 2 K per step, and from rest:
    // three steps each, the kit's energy and entropy audits closing at every
    // one.
    let rest = SpatialWindingState::at_temperatures(circuit.clone(), vec![START_K; 3]).unwrap();
    let charged = charged(&winding, &circuit);
    for (start, rise) in [(&charged, 0.0_f64), (&charged, 2.0), (&rest, 0.0)] {
        let mut state = start.clone();
        for k in 1..=3 {
            let t = rise.mul_add(f64::from(k), START_K);
            let inputs = [VOLTS, t, t, START_K];
            let receipt = kit(&winding, &state, &inputs);
            assert!(receipt.entropy_generated_j_per_k > 0.0);
            state = proposed(&winding, &state, &inputs).unwrap().0;
        }
    }
}
