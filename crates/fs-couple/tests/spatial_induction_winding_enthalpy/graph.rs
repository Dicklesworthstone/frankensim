//! Graphs of tests 1 to 4: the fixture's winding as fs-em's `SpatialWinding`
//! (feature `component`), built with `magnetic_at_port_temperatures`, beside
//! the shared kernel's backward-Euler cells (`winding_graph/kernel.rs`):
//! - one heat port receives the heat of all 48 tetrahedra, strands and core
//!   alike, and of the authored graph resistor; in graph 4 a second port
//!   receives that resistor's heat in its own body ([`Bodies`]);
//! - every conductor follows the cell's card at the temperature the port
//!   receives ([`CardConductors`]), rebound as the module rebinds it, and
//!   graph 4's resistor its body's card ([`BodyResistor`]);
//! - every tetrahedron's permeability follows the card's `mu(T)` law there
//!   ([`mu_laws`]): the card's permeability has a temperature axis, so the
//!   module and the Component both refuse a temperature refresh of a field
//!   that holds it fixed;
//! - the source port holds 10 V.
//!
//! Graphs 5 and 6, of the play-memory owner, are in [`hysteretic`].

use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_em::spatial::vector_potential::induction::nonlinear::winding::SpatialNonlinearWindingStep;
use fs_em::winding::component::ResistorLaw;
use fs_em::winding::component::spatial::{
    SpatialConductorLaw, SpatialHeatReceiver, SpatialLawElement, SpatialMagneticLaw,
    SpatialWinding, SpatialWindingError, SpatialWindingState,
};
use fs_matdb::QueryPoint;
use fs_material::conductor::ResolvedConductor;
use fs_material::state_point::{ScalarPropertyRequirement, resolve_material_state_point};
use fs_system::{NodeId, Phase, System, SystemError};
use fs_thermal::component::TimelineState;

use super::kernel::{card_words, conduction, graph, ledger_closes, timeline, with_cx, with_gate};
use super::*;

#[cfg(feature = "magnetic-hysteresis")]
#[path = "hysteretic_graph.rs"]
mod hysteretic;

const STEPS: usize = 2;
const DRIVE: f64 = 10.0;
const TOLERANCE_J: f64 = 1e-9;
/// Graph 2's boundary conductance to 290 K \[W/K\], as test 2 authors it.
const COOLING_W_K: f64 = 3.0;
/// Graph 4's resistor body, as test 4 authors it.
const RES_CELL: CellId = CellId(8);
const RES_MAT: MaterialId = MaterialId(10);

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
        let context = source::absolute_point(300.0);
        let (core, strands) = (circuit.field().conductors(), circuit.conductors());
        Self::on(card, context, core, strands, trip)
    }

    /// The law over any owner's `core` and `strands`, each cell's query
    /// point `context`.
    fn on(
        card: &MaterialCard,
        context: QueryPoint,
        core: &[Option<ResolvedConductorState>],
        strands: &[SpatialWindingConductor],
        trip: Option<(usize, &'static CancelGate)>,
    ) -> Self {
        Self {
            card: card.clone(),
            context,
            core: core
                .iter()
                .enumerate()
                .filter_map(|(cell, conductor)| Some((cell, conductor.clone()?)))
                .collect(),
            strands: strands
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

/// A graph resistor's law: its body's card at the temperature its port
/// receives, at that body's query `point`, as a `length_m` long conductor of
/// `area_m2`. As the module's `CellResistor` does, it resolves the
/// resistivity alone, from the claim its single-claim selection picks at
/// 300 K.
#[derive(Debug)]
struct BodyResistor {
    card: MaterialCard,
    point: fn(f64) -> QueryPoint,
    length_m: f64,
    area_m2: f64,
    requirement: ScalarPropertyRequirement,
    pins: MaterialPropertySelection,
}

impl BodyResistor {
    /// Graph 4's: 1 m long, 1 m^2, as test 4 authors it.
    fn new(card: &MaterialCard) -> Self {
        Self::sized(card, source::absolute_point, 1.0, 1.0)
    }

    fn sized(
        card: &MaterialCard,
        point: fn(f64) -> QueryPoint,
        length_m: f64,
        area_m2: f64,
    ) -> Self {
        let requirement = BASIS.requirement().unwrap();
        let selected = resolve_material_state_point(
            card,
            &point(300.0),
            std::slice::from_ref(&requirement),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let claim = selected.properties()[0].answer().receipt.selected;
        Self {
            card: card.clone(),
            point,
            length_m,
            area_m2,
            pins: MaterialPropertySelection::PinnedByProperty(
                [(requirement.name().to_owned(), claim)]
                    .into_iter()
                    .collect(),
            ),
            requirement,
        }
    }
}

impl ResistorLaw for BodyResistor {
    fn conductor(&self, kelvin: f64) -> Result<ResolvedConductor, String> {
        let material = resolve_material_state_point(
            &self.card,
            &(self.point)(kelvin),
            std::slice::from_ref(&self.requirement),
            self.pins.clone(),
        )
        .map_err(|e| format!("{e:?}"))?;
        ResolvedConductorState::try_from_material_with_basis(material, BASIS)
            .and_then(|state| state.into_uniform_conductor(self.length_m, self.area_m2))
            .map_err(|e| e.to_string())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        card_words(&self.card, words);
    }
}

/// One accepted interval: each cell's end temperature \[K\], in cell order;
/// the circuit's currents \[A\], state and field potential \[Wb\]; the
/// resistivity \[ohm m\] of its strands, then of its core; the resistance
/// \[ohm\] of each graph resistor the step bound to a material; the field's
/// end flux density \[T\] in each tetrahedron; and the step's strand, eddy,
/// graph and total heat, its field's and its own numerical dissipation, its
/// stored change and its supplied work \[J\].
#[derive(Debug)]
struct Row {
    kelvin: Vec<f64>,
    currents_a: Vec<f64>,
    state: Vec<f64>,
    potential_wb: Vec<f64>,
    resistivity: Vec<f64>,
    resistance_ohm: Vec<f64>,
    flux_density_t: Vec<[f64; 3]>,
    energies_j: [f64; 8],
}

impl Row {
    fn new(
        kelvin: Vec<f64>,
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
            resistance_ohm: step
                .graph_materials
                .iter()
                .map(|m| m.conductor.resistance_ohm())
                .collect(),
            flux_density_t: step.field.endpoint_flux_density_t.clone(),
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
        self.kelvin
            .iter()
            .chain(&self.currents_a)
            .chain(&self.state)
            .chain(&self.potential_wb)
            .chain(&self.resistivity)
            .chain(&self.resistance_ohm)
            .chain(self.flux_density_t.iter().flatten())
            .chain(&self.energies_j)
            .map(|value| value.to_bits())
            .collect()
    }
}

/// Each cell's temperature \[K\], in cell order.
fn kelvins(cells: &DepositionSimulation<'_>) -> Vec<f64> {
    cells
        .cells()
        .keys()
        .map(|&cell| cells.temperature_k(cell).unwrap())
        .collect()
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

/// The `mu` laws as the module takes them.
fn thermal_laws(
    mu: &[ThermalMagneticMaterial],
) -> Vec<SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>> {
    mu.iter()
        .enumerate()
        .map(|(tetrahedron, law)| SpatialMagneticThermalLaw {
            tetrahedron,
            law: MagneticConstitutiveThermalLaw::Linear(law.clone()),
            maximum_flux_density_t: 10.0,
        })
        .collect()
}

/// A module cell on `card`'s curve, its temperature on the `T` axis of
/// `context`.
fn material_cell<'a>(
    (cell, material_id): (CellId, MaterialId),
    curve: &'a ResolvedHeatCapacityCurve,
    card: &'a MaterialCard,
    context: &'a QueryPoint,
) -> MaterialCell<'a> {
    MaterialCell {
        cell,
        material_id,
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }
}

/// The cells beside the winding: the coil's cell, which receives every
/// tetrahedron's heat and, without `resistor`, the graph resistor's; with
/// `resistor`, the graph resistor's own body, which receives that heat and
/// sets the resistor through the law.
struct Bodies {
    start: DepositionSimulation<'static>,
    materials: Vec<(MaterialId, &'static ResolvedHeatCapacityCurve)>,
    resistor: Option<Arc<dyn ResistorLaw>>,
}

impl Bodies {
    /// The coil's cell alone, its boundary taking `cooling_w_k` \[W/K\] to
    /// 290 K.
    fn coil(curve: &'static ResolvedHeatCapacityCurve, cooling_w_k: f64) -> Self {
        Self {
            start: thermal(curve, cooling_w_k),
            materials: vec![(MAT, curve)],
            resistor: None,
        }
    }

    /// The coil's cell and the resistor's body, 1 kg each at 300 K, both
    /// insulated, as test 4 authors them.
    fn pair(
        coil: &'static ResolvedHeatCapacityCurve,
        body: &'static ResolvedHeatCapacityCurve,
        resistor: Arc<dyn ResistorLaw>,
    ) -> Self {
        let birth = |id: CellId, material_id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: 0,
                body_id: id.0,
                material_id,
                mass_kg: 1.0,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        };
        let start = DepositionSimulation::new(
            vec![birth(CELL, MAT), birth(RES_CELL, RES_MAT)],
            &[
                (MAT, coil.curve() as &dyn ThermalLaw),
                (RES_MAT, body.curve() as &dyn ThermalLaw),
            ],
            DepositionBudget {
                max_cells: 2,
                max_links: 2,
                max_events: 2,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        Self {
            start,
            materials: vec![(MAT, coil), (RES_MAT, body)],
            resistor: Some(resistor),
        }
    }
}

/// The Component's heat receivers, their cells' ports and its resistor laws.
type Layout = (
    Vec<SpatialHeatReceiver>,
    Vec<CellId>,
    BTreeMap<usize, Arc<dyn ResistorLaw>>,
);

/// The coil's receiver takes all 48 tetrahedra and, without `resistor`, the
/// graph resistor's heat; with `resistor`, a second receiver at the
/// resistor's body takes that heat.
fn layout(resistor: Option<Arc<dyn ResistorLaw>>) -> Layout {
    let branch = ROUTES[0].branch;
    let coil = |resistors| SpatialHeatReceiver {
        cells: (0..48).collect(),
        resistors,
    };
    match resistor {
        None => (vec![coil(vec![branch])], vec![CELL], BTreeMap::new()),
        Some(resistor) => (
            vec![
                coil(Vec::new()),
                SpatialHeatReceiver {
                    cells: Vec::new(),
                    resistors: vec![branch],
                },
            ],
            vec![CELL, RES_CELL],
            BTreeMap::from([(branch, resistor)]),
        ),
    }
}

/// The fixture's winding as the Component beside `bodies`' timeline: one
/// receiver per body, each at its cell's port.
fn kernel(
    bodies: Bodies,
    circuit: &SpatialNonlinearWindingCircuit<'static>,
    law: Arc<CardConductors>,
    mu: &[ThermalMagneticMaterial],
) -> (System, NodeId, NodeId) {
    let (receivers, ports, resistors) = layout(bodies.resistor);
    let temperatures = vec![300.0; ports.len()];
    // One law per receiver; the resistor body's has no cells to ask it.
    let law = law as Arc<dyn SpatialConductorLaw<ResolvedConductorState>>;
    let winding = SpatialWinding::magnetic_at_port_temperatures(
        circuit,
        vec![law; receivers.len()],
        resistors,
        mu.iter()
            .enumerate()
            .map(|(cell, law)| SpatialMagneticLaw {
                cell,
                law: law.clone(),
            })
            .collect(),
        receivers,
        controls(),
        TOLERANCE_J,
    )
    .unwrap();
    let start = &bodies.start;
    graph(
        winding,
        SpatialWindingState::at_temperatures(circuit.clone(), temperatures).unwrap(),
        timeline(start, conduction(start, &bodies.materials, &[])),
        &ports,
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
        kelvins(timeline.simulation()),
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

/// The module's conductors at `t` \[K\] are `law`'s there, bit for bit.
fn follows(
    law: &CardConductors,
    (strands, core): (
        &[SpatialWindingConductor],
        &[Option<ResolvedConductorState>],
    ),
    t: f64,
    interval: usize,
) {
    for strand in strands {
        let (winding, cell) = (strand.winding, strand.cell);
        assert_eq!(
            law.strand(winding, cell, t).unwrap(),
            strand.conductor,
            "interval {interval}, winding {winding}, cell {cell}"
        );
    }
    for (cell, conductor) in core.iter().enumerate() {
        if let Some(conductor) = conductor {
            assert_eq!(
                &law.core(cell, t).unwrap(),
                conductor,
                "interval {interval}, cell {cell}"
            );
        }
    }
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
    let mu = thermal_laws(mu);
    let query = source::absolute_point(300.0);
    let cells = [material_cell((CELL, MAT), curve, card, &query)];
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
        let (t, c) = (kelvins(next.thermal()), next.circuit());
        let conductors = (c.conductors(), c.field().conductors());
        follows(law, conductors, t[0], interval);
        let (thermal, circuit, report) = next.into_parts();
        rows.push(Row::new(t, &circuit, &report.coupled.electrical));
        (heat, electrical) = (thermal, circuit);
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
        let (t, i) = (row.kelvin[0], row.currents_a[0]);
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
        take(0, &a.kelvin, &b.kelvin);
        take(1, &a.currents_a, &b.currents_a);
        take(2, &a.state, &b.state);
        take(3, &a.potential_wb, &b.potential_wb);
        take(4, &a.energies_j, &b.energies_j);
    }
    most
}

/// Bounds of [`largest`], in its order, for graphs 1 and 2. On the export of
/// 128e8869d7 graph 1's kernel and module agree within 1.308e-10 K,
/// 1.553e-12 A, 3.553e-15, 6.539e-14 Wb and 3.105e-12 J over both cards (the
/// fixed card's currents, state, potential and energies to the bit); on
/// e53a1fccd9's, graph 2's within 2.022e-10 K, 2.457e-12 A, 7.105e-15,
/// 1.001e-13 Wb and 4.913e-12 J. The bounds sit 51 to 56 times above graph
/// 1's values and 28 to 35 times above graph 2's.
const MODULE: [f64; 5] = [7e-9, 8e-11, 2e-13, 3.5e-12, 1.6e-10];

/// Graph 2's bound on its reversible heat against the module's \[J\]: they
/// agree within 1.346e-9 J (each about 7e-3 to 9e-3 J), about the cell's
/// 4 J/K times their temperature difference; the bound sits 52 times above.
const REVERSIBLE_J: f64 = 7e-8;

/// Graph 4's bounds, in [`largest`]'s order: on e53a1fccd9's export the
/// kernel and the module agree within 2.083e-10 K, 1.234e-11 A, 3.553e-15,
/// 5.760e-13 Wb and 4.989e-11 J over both resistor cards (the fixed card's
/// currents, state, potential and energies to the bit). The bounds sit 48
/// to 56 times above.
const SEPARATE: [f64; 5] = [1e-8, 6e-10, 2e-13, 3e-11, 2.5e-9];

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
        let (mut system, winding, cells) =
            kernel(Bodies::coil(caloric, 0.0), &start, law.clone(), &mu);
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
    assert!(responsive[1].kelvin[0] > responsive[0].kelvin[0]);
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
    let (mut system, winding, cells) = kernel(Bodies::coil(caloric, 0.0), &start, law, &mu);
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
    let (mut system, winding, cells) = kernel(Bodies::coil(caloric, 0.0), &start, law.clone(), &mu);
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
    let (mut system, winding, cells) = kernel(Bodies::coil(caloric, 0.0), &start, law.clone(), &mu);
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

/// The module's run of graph 2: test 2's heat assembly from `circuit` and
/// the cooled cell with the `mu` laws, its conductors held to `law` at its
/// own end temperature, bit for bit, and its external heat the boundary's;
/// each interval's row and reversible heat \[J\].
fn assembly_run(
    card: &MaterialCard,
    curve: &ResolvedHeatCapacityCurve,
    circuit: &SpatialNonlinearWindingCircuit<'_>,
    law: &CardConductors,
    mu: &[ThermalMagneticMaterial],
) -> (Vec<Row>, Vec<f64>) {
    let laws = thermal_laws(mu);
    let query = source::absolute_point(300.0);
    let cells = [material_cell((CELL, MAT), curve, card, &query)];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let regions = regions();
    let (mut heat, mut now) = (thermal(curve, COOLING_W_K), circuit.clone());
    let (mut rows, mut reversible) = (Vec::new(), Vec::new());
    with_cx(|cx| {
        for interval in 1..=STEPS {
            let prepared = assembly.prepare(cx, &heat).unwrap();
            let next = heat_assembly_candidate(
                &prepared,
                cx,
                &heat,
                DT * interval as f64,
                &now,
                &[DRIVE],
                &regions,
                &ROUTES,
                1e-12,
                1e-12,
                Some(&laws),
                thermal_controls(),
                controls(),
                &[0.0],
            )
            .unwrap();
            let (t, c) = (kelvins(next.thermal()), next.circuit());
            let conductors = (c.conductors(), c.field().conductors());
            follows(law, conductors, t[0], interval);
            let (thermal, circuit, report) = next.into_parts();
            near(
                report.external_heat_j,
                -DT * COOLING_W_K * (t[0] - 290.0),
                5e-7,
            );
            reversible.push(report.coupled.source.reversible_heat_j);
            rows.push(Row::new(t, &circuit, &report.coupled.electrical));
            (heat, now) = (thermal, circuit);
        }
    });
    (rows, reversible)
}

/// Graph 2: test 2's cooled cell, its permeability rising 0.0004 per kelvin,
/// over two intervals. On each, the boundary gives the cell
/// `-DT 3 (t - 290)` J, and the cell's `4 dT` is that, the owner's heat and
/// the law's reversible heat, which is not zero. The winding's stored change
/// is the field's internal energy, `b^2 / (2 mu) (1 + T mu' / mu)` over the
/// tetrahedra, and the end inductance's `L i^2 / 2`. The module's heat
/// assembly, test 2's own path, agrees within [`MODULE`] and
/// [`REVERSIBLE_J`].
#[test]
fn the_declared_graph_closes_cooling_and_the_permeability_s_reversible_heat() {
    let slope = 0.0004;
    let card: &'static MaterialCard = Box::leak(Box::new(source::absolute_affine_magnetic_card(
        0.2, 0.005, 2.0, slope, BASIS,
    )));
    let caloric: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(curve(card)));
    let p: &'static SpatialNonlinearVectorMagnetostatics = Box::leak(Box::new(problem(card)));
    let start = circuit(p, 1.0);
    let law = Arc::new(CardConductors::new(card, &start, None));
    let mu = mu_laws(card, p);
    let (mut system, winding, cells) =
        kernel(Bodies::coil(caloric, COOLING_W_K), &start, law.clone(), &mu);
    let rows = kernel_run(&mut system, winding, cells);
    let (mut before, mut cooled, mut reversible) = (300.0, 0.0, Vec::new());
    for row in &rows {
        let t = row.kelvin[0];
        let cooling = -DT * COOLING_W_K * (t - 290.0);
        reversible.push(4.0 * (t - before) - row.energies_j[3] - cooling);
        (before, cooled) = (t, cooled + cooling);
    }
    ledger_closes(&system, winding, cells, cooled);
    let last = &rows[STEPS - 1];
    let (t, i) = (last.kelvin[0], last.currents_a[0]);
    let mu_t = 2.0 + slope * (t - 300.0);
    let field_u: f64 = last
        .flux_density_t
        .iter()
        .map(|&b| dot(b, b) / (96.0 * mu_t) * (1.0 + t * slope / mu_t))
        .sum();
    let stored = system.ledger().nodes[winding.index()]
        .stored_change_j
        .value();
    near(stored, field_u + 0.5 * END_L * i * i, 5e-7);

    let (module, module_reversible) = assembly_run(card, caloric, &start, &law, &mu);
    let most = largest(&rows, &module);
    let apart = reversible
        .iter()
        .zip(&module_reversible)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!(
        "graph 2: reversible heat {reversible:?} J; kernel against module: {:.3e} K, {:.3e} A, {:.3e}, {:.3e} Wb, {:.3e} J, reversible {apart:.3e} J",
        most[0], most[1], most[2], most[3], most[4]
    );
    for (n, r) in reversible.iter().enumerate() {
        assert!(r.abs() > 1e-10, "interval {n}: {r}");
    }
    for (i, (value, bound)) in most.iter().zip(MODULE).enumerate() {
        assert!(*value <= bound, "quantity {i}: {most:?}");
    }
    assert!(apart <= REVERSIBLE_J, "{apart}");
}

/// Test 4's balance on each interval: the body's resistance is
/// `GRAPH_R (1 + alpha (t_r - 300))` at its end temperature; its heat
/// `DT r i^2` is the body's `4 dT`; the coil cell's `4 dT` is the strand and
/// eddy heat.
fn separate(rows: &[Row], alpha: f64) {
    let mut before = [300.0; 2];
    for (n, row) in rows.iter().enumerate() {
        let ([coil, body], i) = ([row.kelvin[0], row.kelvin[1]], row.currents_a[0]);
        let [strand, eddy, resistor, ..] = row.energies_j;
        let r = GRAPH_R * (1.0 + alpha * (body - 300.0));
        assert_eq!(row.resistance_ohm.len(), 1, "interval {n}");
        near(row.resistance_ohm[0], r, 1e-12);
        near(resistor, DT * r * i * i, 3e-7);
        near(4.0 * (body - before[1]), resistor, 5e-7);
        near(4.0 * (coil - before[0]), strand + eddy, 5e-7);
        before = [coil, body];
    }
}

/// The module's run of graph 4: the material step from the coil's and the
/// body's cells (`start`, on these cards and curves) and `circuit` with the
/// `mu` laws, its conductors held to the coil's law and its resistor to the
/// body's at its own end temperatures, bit for bit.
fn material_run(
    [(coil, coil_curve), (body, body_curve)]: [(&MaterialCard, &ResolvedHeatCapacityCurve); 2],
    start: DepositionSimulation<'_>,
    circuit: &SpatialNonlinearWindingCircuit<'_>,
    (law, resistor): (&CardConductors, &BodyResistor),
    mu: &[ThermalMagneticMaterial],
) -> Vec<Row> {
    let query = source::absolute_point(300.0);
    let cells = [
        material_cell((CELL, MAT), coil_curve, coil, &query),
        material_cell((RES_CELL, RES_MAT), body_curve, body, &query),
    ];
    let routes = [ResistorCell {
        branch: 1,
        cell: RES_CELL,
    }];
    let resistors = [CellResistor {
        destination: routes[0],
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }];
    let regions = regions();
    let layout = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &routes,
        fill_tolerance: 1e-12,
    };
    let laws = thermal_laws(mu);
    let (mut heat, mut now) = (start, circuit.clone());
    let mut rows = Vec::new();
    for interval in 1..=STEPS {
        let next = layout
            .step_material_candidate(
                &heat,
                DT * interval as f64,
                &now,
                &[DRIVE],
                &resistors,
                &[(1, BASIS)],
                Some(&laws),
                thermal_controls(),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let (t, c) = (kelvins(next.thermal()), next.circuit());
        let conductors = (c.conductors(), c.field().conductors());
        follows(law, conductors, t[0], interval);
        let (thermal, circuit, report) = next.into_parts();
        let step = &report.coupled.electrical;
        assert_eq!(
            resistor.conductor(t[1]).unwrap(),
            step.graph_materials[0].conductor,
            "interval {interval}"
        );
        rows.push(Row::new(t, &circuit, step));
        (heat, now) = (thermal, circuit);
    }
    rows
}

/// Graph 4: test 4's resistor in its own body, on a fixed and a responsive
/// resistor card, over two intervals. Both runs hold [`separate`]; the
/// module's material step agrees within [`SEPARATE`], its resistor its
/// body's card at the body's end temperature, bit for bit; a construction
/// circuit holding the resistor at 301 K steps as the original does; the
/// responsive body's currents stay below the fixed one's.
#[test]
fn the_declared_graph_heats_a_separate_resistor_body_that_sets_its_own_resistance() {
    let mut currents = Vec::new();
    for alpha in [0.0, 0.005] {
        let coil: &'static MaterialCard = Box::leak(Box::new(
            source::absolute_affine_magnetic_card(0.2, 0.0, 2.0, 0.0, BASIS),
        ));
        let body: &'static MaterialCard = Box::leak(Box::new(
            source::absolute_affine_magnetic_card(GRAPH_R, alpha, 2.0, 0.0, BASIS),
        ));
        let coil_curve: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(curve(coil)));
        let body_curve: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(curve(body)));
        let p: &'static SpatialNonlinearVectorMagnetostatics = Box::leak(Box::new(problem(coil)));
        let start = circuit(p, 1.0);
        let law = Arc::new(CardConductors::new(coil, &start, None));
        let resistor = Arc::new(BodyResistor::new(body));
        let mu = mu_laws(coil, p);
        let bodies = || Bodies::pair(coil_curve, body_curve, resistor.clone());
        let (mut system, winding, cells) = kernel(bodies(), &start, law.clone(), &mu);
        let rows = kernel_run(&mut system, winding, cells);
        ledger_closes(&system, winding, cells, 0.0);
        separate(&rows, alpha);

        // The module's `admit_initial_graph_sources` refuses a circuit whose
        // resistor disagrees with its body's thermal state. The kernel has
        // no input that holds a resistor against its law: the law rebinds it
        // before the solve, so a construction circuit carrying the body's
        // card at 301 K steps to the original's bits.
        let wrong = resistor.conductor(301.0).unwrap();
        let mismatched = start
            .with_graph_resistors(
                vec![GraphResistorMaterial {
                    branch: 1,
                    conductor: wrong,
                }],
                &mut || false,
            )
            .unwrap();
        let (mut other, other_winding, other_cells) =
            kernel(bodies(), &mismatched, law.clone(), &mu);
        with_cx(|cx| other.step(DT, cx)).unwrap();
        assert_eq!(
            row(&other, other_winding, other_cells).bits(),
            rows[0].bits()
        );

        let module = material_run(
            [(coil, coil_curve), (body, body_curve)],
            bodies().start,
            &start,
            (&law, &resistor),
            &mu,
        );
        separate(&module, alpha);
        let most = largest(&rows, &module);
        println!(
            "graph 4, alpha {alpha}: kernel against module: {:.3e} K, {:.3e} A, {:.3e}, {:.3e} Wb, {:.3e} J",
            most[0], most[1], most[2], most[3], most[4]
        );
        for (i, (value, bound)) in most.iter().zip(SEPARATE).enumerate() {
            assert!(*value <= bound, "alpha {alpha}, quantity {i}: {most:?}");
        }
        currents.push(rows.iter().map(|row| row.currents_a[0]).collect::<Vec<_>>());
    }
    let [fixed, responsive] = [&currents[0], &currents[1]];
    for (n, (r, f)) in responsive.iter().zip(fixed).enumerate() {
        assert!(r < f, "interval {n}");
    }
}
