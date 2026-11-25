//! Graphs of tests 1 and 3: the fixture's winding as fs-em's `SpatialWinding`
//! (feature `component`), built with `magnetic_at_port_temperatures`, beside
//! the shared kernel's backward-Euler cell (`winding_graph/kernel.rs`):
//! - one heat port receives the heat of all 48 tetrahedra, strands and core
//!   alike, and of the authored graph resistor;
//! - every conductor follows the cell's card at the temperature the port
//!   receives ([`CardConductors`]), rebound as the module rebinds it;
//! - every tetrahedron's permeability follows the card's `mu(T)` law there
//!   ([`mu_laws`]): the card's permeability has a temperature axis, so the
//!   module and the Component both refuse a temperature refresh of a field
//!   that holds it fixed;
//! - the source port holds 10 V.

use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_em::spatial::vector_potential::induction::nonlinear::winding::SpatialNonlinearWindingStep;
use fs_em::winding::component::spatial::{
    SpatialConductorLaw, SpatialHeatReceiver, SpatialLawElement, SpatialMagneticLaw,
    SpatialWinding, SpatialWindingError, SpatialWindingState,
};
use fs_matdb::QueryPoint;
use fs_system::{NodeId, Phase, System, SystemError};
use fs_thermal::component::TimelineState;

use super::kernel::{card_words, conduction, graph, ledger_closes, timeline, with_cx, with_gate};
use super::*;

const STEPS: usize = 2;
const DRIVE: f64 = 10.0;
const TOLERANCE_J: f64 = 1e-9;

/// Each tetrahedron's conductor in the construction circuit, strand and core
/// alike, rebound to its material at the cell's query point with `kelvin` on
/// the `T` axis, as the module rebinds them. It counts its answers; with
/// `trip`, it requests the gate's cancellation at the answer named.
#[derive(Debug)]
struct CardConductors {
    card: MaterialCard,
    context: QueryPoint,
    core: BTreeMap<usize, ResolvedConductorState>,
    strands: BTreeMap<(usize, usize), ResolvedConductorState>,
    answers: AtomicUsize,
    trip: Option<(usize, &'static CancelGate)>,
}

impl CardConductors {
    fn new(
        card: &MaterialCard,
        circuit: &SpatialNonlinearWindingCircuit<'_>,
        trip: Option<(usize, &'static CancelGate)>,
    ) -> Self {
        Self {
            card: card.clone(),
            context: source::absolute_point(300.0),
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
            answers: AtomicUsize::new(0),
            trip,
        }
    }

    fn at(
        &self,
        original: Option<&ResolvedConductorState>,
        kelvin: f64,
    ) -> Result<ResolvedConductorState, String> {
        let answer = self.answers.fetch_add(1, Ordering::Relaxed) + 1;
        if let Some((at, gate)) = self.trip
            && at == answer
        {
            gate.request();
        }
        let original = original.ok_or("no conductor there")?;
        let point = self
            .context
            .clone()
            .with("T", kelvin)
            .map_err(|e| format!("{e:?}"))?;
        let material = original
            .material()
            .at_query_point(&self.card, &point)
            .map_err(|e| format!("{e:?}"))?;
        original.with_material(material).map_err(|e| e.to_string())
    }
}

impl SpatialConductorLaw<ResolvedConductorState> for CardConductors {
    fn core(&self, cell: usize, kelvin: f64) -> Result<ResolvedConductorState, String> {
        self.at(self.core.get(&cell), kelvin)
    }

    fn strand(
        &self,
        winding: usize,
        cell: usize,
        kelvin: f64,
    ) -> Result<ResolvedConductorState, String> {
        self.at(self.strands.get(&(winding, cell)), kelvin)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        card_words(&self.card, words);
    }
}

/// One accepted interval: the cell's end temperature \[K\]; the circuit's
/// currents \[A\], state and field potential \[Wb\]; the resistivity
/// \[ohm m\] of its strands, then of its core; and the step's strand, eddy,
/// graph and total heat, its field's and its own numerical dissipation, its
/// stored change and its supplied work \[J\].
#[derive(Debug)]
struct Row {
    kelvin: f64,
    currents_a: Vec<f64>,
    state: Vec<f64>,
    potential_wb: Vec<f64>,
    resistivity: Vec<f64>,
    energies_j: [f64; 8],
}

impl Row {
    fn new(
        kelvin: f64,
        circuit: &SpatialNonlinearWindingCircuit<'_>,
        step: &SpatialNonlinearWindingStep,
    ) -> Self {
        let strands = circuit.conductors().iter().map(|s| &s.conductor);
        Self {
            kelvin,
            currents_a: circuit.currents_a().to_vec(),
            state: circuit.state().to_vec(),
            potential_wb: circuit.field().potential_wb().to_vec(),
            resistivity: strands
                .chain(circuit.field().conductors().iter().flatten())
                .map(ResolvedConductorState::resistivity_ohm_m)
                .collect(),
            energies_j: [
                step.cells.iter().map(|cell| cell.heat_j).sum(),
                step.field.heat_j,
                step.graph_resistors.iter().map(|r| r.heat_j).sum(),
                step.heat_j,
                step.field.numerical_dissipation_j,
                step.numerical_dissipation_j,
                step.stored_change_j,
                step.supplied_work_j,
            ],
        }
    }

    /// Every value's bits, in field order.
    fn bits(&self) -> Vec<u64> {
        [self.kelvin]
            .iter()
            .chain(&self.currents_a)
            .chain(&self.state)
            .chain(&self.potential_wb)
            .chain(&self.resistivity)
            .chain(&self.energies_j)
            .map(|value| value.to_bits())
            .collect()
    }
}

/// The card's linear `mu(T)` law over 280-350 K, from each tetrahedron's
/// field material, as test 2 builds it.
fn mu_laws(
    card: &MaterialCard,
    p: &SpatialNonlinearVectorMagnetostatics,
) -> Vec<ThermalMagneticMaterial> {
    p.materials()
        .iter()
        .map(|m| ThermalMagneticMaterial::new(card, m.material(), "T", [280.0, 350.0]).unwrap())
        .collect()
}

/// The fixture's winding as the Component, its one receiver taking every
/// tetrahedron and the graph resistor, beside the cell's timeline.
fn kernel(
    curve: &'static ResolvedHeatCapacityCurve,
    circuit: &SpatialNonlinearWindingCircuit<'static>,
    law: Arc<CardConductors>,
    mu: &[ThermalMagneticMaterial],
) -> (System, NodeId, NodeId) {
    let winding = SpatialWinding::magnetic_at_port_temperatures(
        circuit,
        vec![law as Arc<dyn SpatialConductorLaw<ResolvedConductorState>>],
        BTreeMap::new(),
        mu.iter()
            .enumerate()
            .map(|(cell, law)| SpatialMagneticLaw {
                cell,
                law: law.clone(),
            })
            .collect(),
        vec![SpatialHeatReceiver {
            cells: (0..48).collect(),
            resistors: vec![ROUTES[0].branch],
        }],
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let start = thermal(curve, 0.0);
    graph(
        winding,
        SpatialWindingState::at_temperatures(circuit.clone(), vec![300.0]).unwrap(),
        timeline(&start, conduction(&start, &[(MAT, curve)], &[])),
        &[CELL],
        Some(|_| DRIVE),
    )
}

/// The accepted interval the system holds.
fn row(system: &System, winding: NodeId, cells: NodeId) -> Row {
    let state = system
        .state::<SpatialWindingState<'static>>(winding)
        .unwrap();
    let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
    Row::new(
        timeline.simulation().temperature_k(CELL).unwrap(),
        state.circuit(),
        state.step().unwrap(),
    )
}

fn kernel_run(system: &mut System, winding: NodeId, cells: NodeId) -> Vec<Row> {
    with_cx(|cx| {
        (0..STEPS)
            .map(|_| {
                let report = system.step(DT, cx).unwrap();
                assert_eq!(
                    report.dt_s.to_bits(),
                    DT.to_bits(),
                    "the kernel took the whole interval"
                );
                row(system, winding, cells)
            })
            .collect()
    })
}

/// The module's run: `step_candidate` from `circuit` with the `mu` laws, its
/// conductors held to `law` at its own end temperature, bit for bit.
fn module_run(
    card: &MaterialCard,
    curve: &ResolvedHeatCapacityCurve,
    circuit: &SpatialNonlinearWindingCircuit<'_>,
    law: &CardConductors,
    mu: &[ThermalMagneticMaterial],
) -> Vec<Row> {
    let mu: Vec<_> = mu
        .iter()
        .enumerate()
        .map(|(tetrahedron, law)| SpatialMagneticThermalLaw {
            tetrahedron,
            law: MagneticConstitutiveThermalLaw::Linear(law.clone()),
            maximum_flux_density_t: 10.0,
        })
        .collect();
    let query = source::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CELL,
        material_id: MAT,
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context: &query,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions = regions();
    let network = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &ROUTES,
        fill_tolerance: 1e-12,
    };
    let (mut heat, mut electrical) = (thermal(curve, 0.0), circuit.clone());
    let mut rows = Vec::new();
    for interval in 1..=STEPS {
        let next = network
            .step_candidate(
                &heat,
                DT * interval as f64,
                &electrical,
                &[DRIVE],
                Some(&mu),
                thermal_controls(),
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let t = next.thermal().temperature_k(CELL).unwrap();
        for strand in next.circuit().conductors() {
            let (winding, cell) = (strand.winding, strand.cell);
            assert_eq!(
                law.strand(winding, cell, t).unwrap(),
                strand.conductor,
                "interval {interval}, winding {winding}, cell {cell}"
            );
        }
        for (cell, conductor) in next.circuit().field().conductors().iter().enumerate() {
            if let Some(conductor) = conductor {
                assert_eq!(
                    &law.core(cell, t).unwrap(),
                    conductor,
                    "interval {interval}, cell {cell}"
                );
            }
        }
        rows.push(Row::new(
            t,
            next.circuit(),
            &next.report().coupled.electrical,
        ));
        (heat, electrical, _) = next.into_parts();
    }
    rows
}

/// Test 1's equations on each interval: the strand heat is `DT i^2 rho R`,
/// with `R` the profile's geometric resistance per unit resistivity; the
/// graph heat `DT i^2 GRAPH_R`; the heat their sum with the eddy heat; the
/// cell's `4 dT` that heat; and every conductor's resistivity `rho(T)` at
/// the end temperature.
fn closed(rows: &[Row], alpha: f64, geometric_r: f64) {
    let mut before = 300.0;
    for (n, row) in rows.iter().enumerate() {
        let (t, i) = (row.kelvin, row.currents_a[0]);
        let [strand, eddy, resistor, heat, field_numerical, ..] = row.energies_j;
        let rho = 0.2 * (1.0 + alpha * (t - 300.0));
        println!(
            "alpha {alpha}, interval {n}: {t:.9} K, {i:.9} A, strand {strand:.6e} J, eddy {eddy:.6e} J"
        );
        near(strand, DT * i * i * rho * geometric_r, 3e-7);
        near(resistor, DT * i * i * GRAPH_R, 3e-7);
        near(heat, strand + eddy + resistor, 3e-7);
        near(4.0 * (t - before), heat, 5e-7);
        assert!(field_numerical > 0.0, "interval {n}");
        for &value in &row.resistivity {
            near(value, rho, 1e-12);
        }
        before = t;
    }
}

/// The largest differences over the intervals: temperature \[K\], current
/// \[A\], circuit state, field potential \[Wb\] and energy \[J\].
fn largest(kernel: &[Row], module: &[Row]) -> [f64; 5] {
    assert_eq!((kernel.len(), module.len()), (STEPS, STEPS));
    let mut most = [0.0_f64; 5];
    for (a, b) in kernel.iter().zip(module) {
        let mut take = |slot: usize, x: &[f64], y: &[f64]| {
            assert_eq!(x.len(), y.len());
            for (p, q) in x.iter().zip(y) {
                most[slot] = most[slot].max((p - q).abs());
            }
        };
        take(0, &[a.kelvin], &[b.kelvin]);
        take(1, &a.currents_a, &b.currents_a);
        take(2, &a.state, &b.state);
        take(3, &a.potential_wb, &b.potential_wb);
        take(4, &a.energies_j, &b.energies_j);
    }
    most
}

/// Bounds of [`largest`], in its order. On the export of 128e8869d7 the
/// kernel and the module agree within 1.308e-10 K, 1.553e-12 A, 3.553e-15,
/// 6.539e-14 Wb and 3.105e-12 J over both cards (the fixed card's currents,
/// state, potential and energies to the bit); the bounds sit about fifty
/// times above.
const MODULE: [f64; 5] = [7e-9, 8e-11, 2e-13, 3.5e-12, 1.6e-10];

/// Graph 1: test 1's two intervals at a fixed and a responsive card.
#[test]
fn the_declared_graph_reproduces_the_original_temperature_feedback_over_two_steps() {
    let mut runs = Vec::new();
    for alpha in [0.0, 0.005] {
        let card: &'static MaterialCard = Box::leak(Box::new(
            source::absolute_affine_magnetic_card(0.2, alpha, 2.0, 0.0, BASIS),
        ));
        let caloric: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(curve(card)));
        let p: &'static SpatialNonlinearVectorMagnetostatics = Box::leak(Box::new(problem(card)));
        let start = circuit(p, 1.0);
        let law = Arc::new(CardConductors::new(card, &start, None));
        let mu = mu_laws(card, p);
        let (mut system, winding, cells) = kernel(caloric, &start, law.clone(), &mu);
        let rows = kernel_run(&mut system, winding, cells);
        ledger_closes(&system, winding, cells, 0.0);
        let geometric_r = profile(p).iter().map(|&j| dot(j, j) / 48.0).sum::<f64>();
        closed(&rows, alpha, geometric_r);
        let module = module_run(card, caloric, &start, &law, &mu);
        let most = largest(&rows, &module);
        println!(
            "alpha {alpha}: kernel against module: {:.3e} K, {:.3e} A, {:.3e}, {:.3e} Wb, {:.3e} J",
            most[0], most[1], most[2], most[3], most[4]
        );
        for (i, (value, bound)) in most.iter().zip(MODULE).enumerate() {
            assert!(*value <= bound, "alpha {alpha}, quantity {i}: {most:?}");
        }
        runs.push(rows);
    }
    // The feedback is real: the responsive cell warms and its currents stay
    // below the fixed card's.
    let [fixed, responsive] = [&runs[0], &runs[1]];
    assert!(responsive[1].kelvin > responsive[0].kelvin);
    for n in 0..STEPS {
        assert!(responsive[n].currents_a[0] < fixed[n].currents_a[0], "{n}");
    }
}

/// Graph 3: test 3's refusals on the responsive card. A law over another
/// card has no answer; a cancellation requested at the law's last
/// answer in a step is the kernel's `Cancelled`. Neither publishes anything,
/// and a retry equals the good step, bit for bit. The module's partial-fill
/// refusal and circuit budget have no counterpart here: its caloric owner
/// maps mesh volumes to cells and budgets its own thermal Newton, while the
/// graph's cell is the timeline's.
#[test]
fn the_declared_graph_refuses_a_foreign_card_and_a_late_cancellation_and_keeps_every_history() {
    let card: &'static MaterialCard = Box::leak(Box::new(source::absolute_affine_magnetic_card(
        0.2, 0.005, 2.0, 0.0, BASIS,
    )));
    let caloric: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(curve(card)));
    let p: &'static SpatialNonlinearVectorMagnetostatics = Box::leak(Box::new(problem(card)));
    let start = circuit(p, 1.0);
    let mu = mu_laws(card, p);
    let unchanged = |system: &System, winding: NodeId, cells: NodeId| {
        assert_eq!(system.time_s().to_bits(), 0.0_f64.to_bits());
        let state = system
            .state::<SpatialWindingState<'static>>(winding)
            .unwrap();
        let constructed = SpatialWindingState::at_temperatures(start.clone(), vec![300.0]).unwrap();
        assert!(*state == constructed && state.step().is_none());
        let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
        let kelvin = timeline.simulation().temperature_k(CELL).unwrap();
        assert_eq!(kelvin.to_bits(), 300.0_f64.to_bits());
    };

    // A law over the 0.3 ohm m card has no answer at the first core cell:
    // fs-material rebinds a material only to its own card, as in the
    // module's refresh.
    let foreign: &'static MaterialCard = Box::leak(Box::new(
        source::absolute_affine_magnetic_card(0.3, 0.005, 2.0, 0.0, BASIS),
    ));
    let law = Arc::new(CardConductors::new(foreign, &start, None));
    let (mut system, winding, cells) = kernel(caloric, &start, law, &mu);
    let error = with_cx(|cx| system.step(DT, cx)).unwrap_err();
    let SystemError::Component {
        node,
        phase: Phase::Propose,
        error: owner,
    } = &error
    else {
        panic!("{error:?}");
    };
    assert!(
        node == "winding"
            && matches!(
                owner.downcast_ref::<SpatialWindingError>(),
                Some(SpatialWindingError::Law {
                    receiver: 0,
                    element: SpatialLawElement::Core(0),
                    message,
                    ..
                }) if message.contains("original immutable card")
            ),
        "{error:?}"
    );
    unchanged(&system, winding, cells);

    // The good step, its law's answers counted.
    let law = Arc::new(CardConductors::new(card, &start, None));
    let (mut system, winding, cells) = kernel(caloric, &start, law.clone(), &mu);
    with_cx(|cx| system.step(DT, cx)).unwrap();
    let good = row(&system, winding, cells);
    let answers = law.answers.load(Ordering::Relaxed);
    // Each proposal asks once per tetrahedron, and the last comes after the
    // cell has answered the winding's heat at least once.
    assert!(answers.is_multiple_of(48) && answers >= 2 * 48, "{answers}");

    // The same step, cancelled at the law's last answer: the owner's field
    // solve notices; nothing is published, and the retry is the good step.
    let gate: &'static CancelGate = Box::leak(Box::new(CancelGate::new_clock_free()));
    let law = Arc::new(CardConductors::new(card, &start, Some((answers, gate))));
    let (mut system, winding, cells) = kernel(caloric, &start, law.clone(), &mu);
    let error = with_gate(gate, |cx| system.step(DT, cx)).unwrap_err();
    assert!(matches!(error, SystemError::Cancelled), "{error:?}");
    assert_eq!(law.answers.load(Ordering::Relaxed), answers);
    unchanged(&system, winding, cells);
    with_cx(|cx| system.step(DT, cx)).unwrap();
    assert_eq!(row(&system, winding, cells).bits(), good.bits());
    println!(
        "graph 3: {} proposals in the step; the cancelled one was the last",
        answers / 48
    );
}
