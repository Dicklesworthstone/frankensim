//! The module's nonlinear winding trials as declared graphs on the
//! composition kernel.
//!
//! `WindingThermalNetwork::step_nonlinear_candidate` and
//! `step_nonlinear_material_candidate` couple fs-em's nonlinear winding
//! circuit to the enthalpy cells by hand: the circuit's backward-Euler step
//! with every conductor at its cell's END temperature, Newton on the cell
//! energies. The graphs say the same with no code for the pair:
//!
//! - `fs_em::winding::component::PlanarWinding::at_port_temperatures`: the
//!   same circuit, one heat port per destination cell, the winding's
//!   conductors ([`CellConductors`]) and a material graph resistor
//!   ([`BodyResistor`]) following their cell's card at the temperature the
//!   port receives; an authored graph resistor keeps its ohms;
//! - `fs_thermal::component::ThermalTimeline::implicit` with backward
//!   Euler, returning end-of-step temperatures;
//! - one heat route per port into its cell, and the cells' `environment`
//!   port to a sink at 300 K.
//!
//! The laws come from public calls only. At each of the module's end
//! temperatures they give its own conductors, bit for bit, before the graph
//! runs. The two sides solve the same equations with different solvers (the
//! module's Newton on cell energies; the kernel's interface iteration and
//! the implicit conduction graph), so they agree to within the solves'
//! tolerances.

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::circuit_enthalpy::material::CellResistor;
use fs_em::planar::nonlinear::PlanarNonlinearMagnetostatics;
use fs_em::winding::component::{
    ConductorLaw, HeatReceiver, PlanarWinding, PlanarWindingState, ResistorLaw,
};
use fs_em::winding::nonlinear::{NonlinearWindingCircuit, NonlinearWindingStep};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductor, ResolvedConductorState};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_system::{NodeId, RouteTemperature, SolveConfig, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};

use super::graph_material::{
    DESTINATIONS, LINKS, RESISTOR_CELL, RESISTOR_MATERIAL, bindings, graph_card, native, resistors,
    solve_controls, two_bodies,
};
use super::*;

pub(super) const STEPS: usize = 4;
pub(super) const TOLERANCE_J: f64 = 1e-9;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
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

fn card_words(card: &MaterialCard, words: &mut Vec<u64>) {
    for chunk in card.content_hash().0.chunks(8) {
        words.push(u64::from_le_bytes(chunk.try_into().unwrap()));
    }
}

/// A cell's winding conductors at the cell's temperature, as the module
/// re-resolves them: each triangle's conductor in the construction circuit,
/// rebound to its material at the cell's query point with `kelvin` on the
/// `T` axis. That is the point `ThermalMaterialPoint::at` writes for a cell
/// with no liquid fraction axis.
#[derive(Debug)]
pub(super) struct CellConductors {
    card: MaterialCard,
    pub(super) context: QueryPoint,
    originals: BTreeMap<usize, ResolvedConductorState>,
}

impl CellConductors {
    pub(super) fn new(card: &MaterialCard, circuit: &NonlinearWindingCircuit<'_>) -> Self {
        Self {
            card: card.clone(),
            context: field::material::point(300.0),
            originals: circuit
                .conductors()
                .iter()
                .map(|source| (source.triangle, source.conductor.clone()))
                .collect(),
        }
    }
}

impl ConductorLaw for CellConductors {
    fn conductor(&self, triangle: usize, kelvin: f64) -> Result<ResolvedConductorState, String> {
        let original = self
            .originals
            .get(&triangle)
            .ok_or_else(|| format!("triangle {triangle} has no conductor"))?;
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

    fn identity_words(&self, words: &mut Vec<u64>) {
        card_words(&self.card, words);
    }
}

/// A material graph resistor at its cell's temperature, as the module
/// resolves it: the cell's card at its query point with `kelvin` on the `T`
/// axis, the cell's conductivity and the resistor's resistivity pinned to
/// the claims they select at 300 K, over the `CellResistor`'s path.
#[derive(Debug)]
struct BodyResistor {
    card: MaterialCard,
    context: QueryPoint,
    requirements: Vec<ScalarPropertyRequirement>,
    pins: MaterialPropertySelection,
    length_m: f64,
    area_m2: f64,
}

impl BodyResistor {
    fn new(card: &MaterialCard, conductivity: &PropertyKey, resistor: &CellResistor) -> Self {
        let context = field::material::point(300.0);
        let requirements = vec![
            ScalarPropertyRequirement::try_with_key(conductivity, ScalarAdmissibility::NonNegative)
                .unwrap(),
            ElectricalPropertyBasis::Resistivity.requirement().unwrap(),
        ];
        let pins = requirements
            .iter()
            .map(|requirement| {
                let selected = resolve_material_state_point(
                    card,
                    &context,
                    std::slice::from_ref(requirement),
                    resistor.selection.clone(),
                )
                .unwrap();
                (
                    requirement.name().to_owned(),
                    selected.properties()[0].answer().receipt.selected,
                )
            })
            .collect();
        Self {
            card: card.clone(),
            context,
            requirements,
            pins: MaterialPropertySelection::PinnedByProperty(pins),
            length_m: resistor.length_m,
            area_m2: resistor.area_m2,
        }
    }
}

impl ResistorLaw for BodyResistor {
    fn conductor(&self, kelvin: f64) -> Result<ResolvedConductor, String> {
        let point = self
            .context
            .clone()
            .with("T", kelvin)
            .map_err(|e| format!("{e:?}"))?;
        let material =
            resolve_material_state_point(&self.card, &point, &self.requirements, self.pins.clone())
                .map_err(|e| format!("{e:?}"))?;
        ResolvedConductorState::try_from_material_with_basis(
            material,
            ElectricalPropertyBasis::Resistivity,
        )
        .and_then(|state| state.into_uniform_conductor(self.length_m, self.area_m2))
        .map_err(|e| e.to_string())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        card_words(&self.card, words);
        words.extend([self.length_m.to_bits(), self.area_m2.to_bits()]);
    }
}

/// A caloric curve of constant specific heat with the cards' conductivity
/// of 1 W/(m K), for the implicit conduction graph. The temperature is the
/// curve's own.
struct Conduction {
    curve: &'static ResolvedHeatCapacityCurve,
}

impl ConductivityLaw for Conduction {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(1.0)
    }
}

impl ImplicitConductionLaw for Conduction {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.curve.curve()
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        let knots = self.curve.curve().knots();
        [
            knots[0].specific_enthalpy_j_kg,
            knots[knots.len() - 1].specific_enthalpy_j_kg,
        ]
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let knots = self.curve.curve().knots();
        let (first, last) = (&knots[0], &knots[knots.len() - 1]);
        Ok(ImplicitConductionState {
            temperature_k: self.curve.curve().temperature_k(h)?,
            temperature_slope_kg_k_j: (last.temperature_k - first.temperature_k)
                / (last.specific_enthalpy_j_kg - first.specific_enthalpy_j_kg),
            conductivity_w_m_k: 1.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// The cells of `start` as one backward-Euler timeline whose cell ports
/// return end-of-step temperatures, joined by `links`.
fn timeline(
    start: &DepositionSimulation<'static>,
    materials: &[(MaterialId, &'static ResolvedHeatCapacityCurve)],
    links: &[ConductiveLink],
) -> (ThermalTimeline<'static>, TimelineState<'static>) {
    let laws: &'static [Conduction] = Box::leak(
        materials
            .iter()
            .map(|&(_, curve)| Conduction { curve })
            .collect(),
    );
    let ids = materials.iter().map(|&(id, _)| id);
    let conductivity: Vec<_> = ids
        .clone()
        .zip(laws)
        .map(|(id, law)| (id, law as &dyn ConductivityLaw))
        .collect();
    let observer = ConductionObserver::new(start, &conductivity, links, &[], || Ok(())).unwrap();
    let implicit: Vec<_> = ids
        .zip(laws)
        .map(|(id, law)| (id, law as &dyn ImplicitConductionLaw))
        .collect();
    let conduction = ImplicitConductionGraph::new(
        start,
        &observer,
        &implicit,
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let controls = ImplicitConductionControls {
        maximum_dt_s: DT,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-9,
        ..Default::default()
    };
    let timeline = ThermalTimeline::implicit(
        start.cells().keys().copied().collect(),
        conduction,
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .unwrap()
    .returning_end_temperatures()
    .unwrap();
    let state = timeline.start(start).unwrap();
    (timeline, state)
}

/// `winding` from `state` beside the cells of `start`, heat port `k`
/// routed into `ports[k]`, its source port (if any) held at `drive` volts.
pub(super) fn graph(
    winding: PlanarWinding<'static>,
    state: PlanarWindingState<'static>,
    (start, materials, links): (
        &DepositionSimulation<'static>,
        &[(MaterialId, &'static ResolvedHeatCapacityCurve)],
        &[ConductiveLink],
    ),
    ports: &[CellId],
    drive: Option<f64>,
) -> (System, NodeId, NodeId) {
    let (timeline, cells_start) = timeline(start, materials, links);
    let mut builder = SystemBuilder::new();
    let electrical = builder.add("winding", winding, state).unwrap();
    let cells = builder.add("cells", timeline, cells_start).unwrap();
    if let Some(volts) = drive {
        builder
            .boundary(
                builder.port(electrical, "source0").unwrap(),
                Box::new(move |_, _, out: &mut [f64]| out.fill(volts)),
            )
            .unwrap();
    }
    for (k, cell) in ports.iter().enumerate() {
        builder
            .heat_route(
                builder.port(electrical, &format!("heat{k}")).unwrap(),
                &[(
                    builder.port(cells, &format!("cell{}", cell.0)).unwrap(),
                    1.0,
                )],
                None,
                RouteTemperature::SinkEquivalent,
            )
            .unwrap();
    }
    builder
        .heat_route(
            builder.port(cells, "environment").unwrap(),
            &[],
            Some(heat_sink("environment", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    // The winding's heat comes out of its owner's nested solves, whose field
    // Newton stops at a relative residual of 1e-10 (`field::controls()`). The
    // interface closes to that, not to the kernel's default 1e-12, which a
    // heat answered to 1e-10 cannot meet.
    system.set_config(SolveConfig {
        relative_tolerance: 1e-10,
        ..SolveConfig::default()
    });
    for k in 0..ports.len() {
        let port = system.port(electrical, &format!("heat{k}")).unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    (system, electrical, cells)
}

/// One accepted interval: the cells' temperatures \[K\], the winding current
/// \[A\] and circuit state, and over the interval the circuit's stored
/// change, Joule heat, numerical dissipation and supplied work \[J\], the
/// entropy \[J/K\] its losses generate at their cells and the entropy
/// conduction among the cells generates.
#[derive(Debug, Clone)]
pub(super) struct Row {
    pub(super) kelvin: Vec<f64>,
    pub(super) current_a: f64,
    pub(super) state: Vec<f64>,
    stored_j: f64,
    joule_j: f64,
    pub(super) numerical_j: f64,
    supplied_j: f64,
    entropy_j_k: [f64; 2],
}

pub(super) fn module_row(
    kelvin: Vec<f64>,
    circuit: &NonlinearWindingCircuit<'_>,
    step: &NonlinearWindingStep,
    entropy_j_k: [f64; 2],
) -> Row {
    Row {
        kelvin,
        current_a: circuit.currents_a()[0],
        state: circuit.state().to_vec(),
        stored_j: step.stored_change_j,
        joule_j: step.heat_j,
        numerical_j: step.numerical_dissipation_j,
        supplied_j: step.supplied_work_j,
        entropy_j_k,
    }
}

pub(super) fn graph_run(
    system: &mut System,
    winding: NodeId,
    cells: NodeId,
    order: &[CellId],
) -> Vec<Row> {
    let mut rows = Vec::new();
    let mut before = [0.0; 6];
    with_cx(|cx| {
        for _ in 0..STEPS {
            let report = system.step(DT, cx).unwrap();
            assert_eq!(
                report.dt_s.to_bits(),
                DT.to_bits(),
                "the kernel took the whole interval"
            );
            let ledger = system.ledger();
            let node = &ledger.nodes[winding.index()];
            let now = [
                node.stored_change_j.value(),
                node.dissipated_j.value(),
                node.numerical_loss_j.value(),
                ledger
                    .boundaries
                    .iter()
                    .map(|b| b.energy_in_j.value())
                    .sum(),
                node.entropy_generated_j_per_k.value(),
                ledger.nodes[cells.index()]
                    .entropy_generated_j_per_k
                    .value(),
            ];
            let circuit = system
                .state::<PlanarWindingState<'static>>(winding)
                .unwrap()
                .circuit();
            let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
            rows.push(Row {
                kelvin: order
                    .iter()
                    .map(|&id| timeline.simulation().temperature_k(id).unwrap())
                    .collect(),
                current_a: circuit.currents_a()[0],
                state: circuit.state().to_vec(),
                stored_j: now[0] - before[0],
                joule_j: now[1] - before[1],
                numerical_j: now[2] - before[2],
                supplied_j: now[3] - before[3],
                entropy_j_k: [now[4] - before[4], now[5] - before[5]],
            });
            before = now;
        }
    });
    rows
}

/// The largest differences over the intervals: temperature \[K\], current
/// \[A\], state coordinate, the four energies \[J\] and the two entropy
/// terms \[J/K\].
pub(super) fn largest(kernel: &[Row], module: &[Row]) -> [f64; 5] {
    assert_eq!((kernel.len(), module.len()), (STEPS, STEPS));
    let mut most = [0.0_f64; 5];
    for (ours, theirs) in kernel.iter().zip(module) {
        assert_eq!(ours.kelvin.len(), theirs.kelvin.len());
        assert_eq!(ours.state.len(), theirs.state.len());
        let pairs = [
            (
                0,
                ours.kelvin.iter().zip(&theirs.kelvin).collect::<Vec<_>>(),
            ),
            (2, ours.state.iter().zip(&theirs.state).collect()),
            (
                3,
                vec![
                    (&ours.stored_j, &theirs.stored_j),
                    (&ours.joule_j, &theirs.joule_j),
                    (&ours.numerical_j, &theirs.numerical_j),
                    (&ours.supplied_j, &theirs.supplied_j),
                ],
            ),
            (
                4,
                ours.entropy_j_k.iter().zip(&theirs.entropy_j_k).collect(),
            ),
        ];
        most[1] = most[1].max((ours.current_a - theirs.current_a).abs());
        for (slot, values) in pairs {
            for (x, y) in values {
                most[slot] = most[slot].max((x - y).abs());
            }
        }
    }
    most
}

/// Bounds of [`largest`], in its order. On the export of d4e7a28481 the
/// kernel and the module agree within 6.6e-11 K, 2.3e-11 A, 3.0e-12 C,
/// 6.8e-11 J and 2.3e-13 J/K over the four intervals of every graph here;
/// the bounds sit about fifty times above that, far below the module's own
/// 1e-8 J cell stop.
const BOUNDS: [f64; 5] = [3e-9, 1.2e-9, 1.5e-10, 3.5e-9, 1.2e-11];

pub(super) fn hold(most: [f64; 5], what: &str) {
    for (i, (value, bound)) in most.iter().zip(BOUNDS).enumerate() {
        assert!(
            *value <= bound,
            "{what}, quantity {i}: {value:e} above {bound:e}"
        );
    }
}

/// The kernel's ledger closes: the energy the source supplied, less the
/// circuit's stored change and numerical dissipation, is the cells' stored
/// change; the environment sink is idle; the winding and the cells generate
/// the ledger's entropy, the routes carrying only the interface iteration's
/// residue.
pub(super) fn ledger_closes(system: &System, winding: NodeId, cells: NodeId) {
    let ledger = system.ledger();
    assert!(
        ledger.energy_residual_j().abs() <= 1e-8,
        "{}",
        ledger.energy_residual_j()
    );
    assert_eq!(ledger.sinks.len(), 1);
    assert!(ledger.sinks[0].heat_j.value().abs() <= 1e-12);
    let node = &ledger.nodes[winding.index()];
    let supplied: f64 = ledger
        .boundaries
        .iter()
        .map(|b| b.energy_in_j.value())
        .sum();
    let left = supplied - node.stored_change_j.value() - node.numerical_loss_j.value();
    let stored = ledger.nodes[cells.index()].stored_change_j.value();
    assert!(
        left > 1.0 && (stored - left).abs() <= 1e-8,
        "{stored} vs {left}"
    );
    let nodes = [winding, cells].map(|n| ledger.nodes[n.index()].entropy_generated_j_per_k.value());
    let total = ledger.entropy_generated_j_per_k();
    // A route's residue is its heat times the difference of the reciprocal
    // temperatures its two ends last saw, a relative 1e-10 of the heat's
    // entropy at most: the interface tolerance.
    assert!(
        (total - nodes[0] - nodes[1]).abs() <= 1e-10 * total,
        "{total} vs {nodes:?}"
    );
}

/// The module's run of graphs 1 and 2: `step_nonlinear_candidate` from
/// `circuit`, its conductors held to `law` at its own end temperature, bit
/// for bit.
fn fixed_module(
    card: &MaterialCard,
    curve: &ResolvedHeatCapacityCurve,
    circuit: &NonlinearWindingCircuit<'_>,
    law: &CellConductors,
) -> Vec<Row> {
    let mut module = Vec::new();
    let mut thermal = thermal(curve);
    let mut circuit = circuit.clone();
    let context = field::material::point(300.0);
    let cells = cells(card, curve, &context);
    let regions = routes();
    let layout = network(&cells, &regions);
    for step in 1..=STEPS {
        let candidate = layout
            .step_nonlinear_candidate(
                &thermal,
                step as f64 * DT,
                &circuit,
                &[DRIVE],
                controls(),
                field::controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let t = candidate.thermal().temperature_k(CELL).unwrap();
        let electrical = &candidate.report().electrical;
        for loss in &electrical.cells {
            let triangle = loss.source.triangle;
            assert_eq!(
                law.conductor(triangle, t).unwrap(),
                loss.source.conductor,
                "step {step}, triangle {triangle}"
            );
        }
        module.push(module_row(
            vec![t],
            candidate.circuit(),
            electrical,
            [electrical.heat_j / t, 0.0],
        ));
        let (next_thermal, next_circuit, _) = candidate.into_parts();
        thermal = next_thermal;
        circuit = next_circuit;
    }
    module
}

/// The module's run of graph 3: `step_nonlinear_material_candidate` from
/// `circuit`, its conductors and graph resistor held to `law` and `resistor`
/// at its own end temperatures, bit for bit.
fn material_module(
    (cards, curves): (&[MaterialCard; 2], &[ResolvedHeatCapacityCurve; 2]),
    circuit: &NonlinearWindingCircuit<'_>,
    law: &CellConductors,
    resistor: &BodyResistor,
    key: &PropertyKey,
) -> Vec<Row> {
    let mut module = Vec::new();
    let mut thermal = two_bodies(curves);
    let mut circuit = circuit.clone();
    let context = field::material::point(300.0);
    let cells = bindings(cards, curves, &context, key);
    let regions = routes();
    let resistors = resistors();
    let layout = WindingThermalNetwork {
        cells: &cells,
        links: &LINKS,
        magnetic_regions: &regions,
        graph_resistors: &DESTINATIONS,
        volume_tolerance_m3: 1e-12,
    };
    let bases = [(2, ElectricalPropertyBasis::Resistivity)];
    for step in 1..=STEPS {
        let candidate = layout
            .step_nonlinear_material_candidate(
                &thermal,
                step as f64 * DT,
                &circuit,
                &[],
                &resistors,
                &bases,
                solve_controls(),
                field::controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let t = [CELL, RESISTOR_CELL].map(|id| candidate.thermal().temperature_k(id).unwrap());
        let report = candidate.report();
        let electrical = &report.electrical;
        for loss in &electrical.cells {
            let triangle = loss.source.triangle;
            assert_eq!(
                law.conductor(triangle, t[0]).unwrap(),
                loss.source.conductor,
                "step {step}, triangle {triangle}"
            );
        }
        assert_eq!(
            resistor.conductor(t[1]).unwrap(),
            electrical.graph_materials[0].conductor,
            "step {step}"
        );
        let winding_j: f64 = electrical.cells.iter().map(|loss| loss.heat_j).sum();
        let link_j = DT * report.source.network.links[0].conductance_w_k * (t[0] - t[1]);
        module.push(module_row(
            t.to_vec(),
            candidate.circuit(),
            electrical,
            [
                winding_j / t[0] + electrical.graph_resistors[0].heat_j / t[1],
                link_j * (1.0 / t[1] - 1.0 / t[0]),
            ],
        ));
        let (next_thermal, next_circuit, _) = candidate.into_parts();
        thermal = next_thermal;
        circuit = next_circuit;
    }
    module
}

/// Graphs 1 and 2: the fixture of
/// `data_only_saturation_changes_actual_voltage_current_and_finite_thermal_feedback`.
/// One cell holds every triangle and the authored 1-ohm graph resistor; no
/// links, no sources.
#[test]
fn the_declared_graph_reproduces_the_fixed_resistor_winding_and_the_bisection_fixed_point() {
    for saturating in [false, true] {
        let card: &'static MaterialCard = Box::leak(Box::new(field::source(saturating, 18.0)));
        let curve: &'static ResolvedHeatCapacityCurve = Box::leak(Box::new(caloric(card)));
        let problem: &'static PlanarNonlinearMagnetostatics =
            Box::leak(Box::new(field::problem(card, DEPTH)));
        let circuit = field::single(problem, card, SERIES_OHM, 0.0, DRIVE);
        let law = Arc::new(CellConductors::new(card, &circuit));

        // The module, its conductors checked against the law at its own end
        // temperature, bit for bit.
        let module = fixed_module(card, curve, &circuit, &law);

        let winding = PlanarWinding::at_port_temperatures(
            &circuit,
            vec![law],
            BTreeMap::new(),
            vec![HeatReceiver {
                triangles: (0..8).collect(),
                resistors: vec![2],
            }],
            field::controls(),
            TOLERANCE_J,
        )
        .unwrap();
        let start = thermal(curve);
        let (mut system, electrical, cells) = graph(
            winding,
            PlanarWindingState::new(circuit),
            (&start, &[(MATERIAL, curve)], &[]),
            &[CELL],
            Some(DRIVE),
        );
        let kernel = graph_run(&mut system, electrical, cells, &[CELL]);
        ledger_closes(&system, electrical, cells);

        // Graph 1: each interval of the kernel against the bisection fixed
        // point from the kernel's own start of interval.
        let (mut from, mut oracle) = ((0.0, 300.0), [0.0_f64; 3]);
        for row in &kernel {
            let expected = reference(from.0, from.1, saturating);
            oracle[0] = oracle[0].max((row.kelvin[0] - expected.temperature).abs());
            oracle[1] = oracle[1].max((row.current_a - expected.current).abs());
            oracle[2] = oracle[2].max((row.joule_j - expected.heat).abs());
            from = (row.current_a, row.kelvin[0]);
        }
        println!(
            "graph 1 (saturating {saturating}): kernel against bisection: {:.1e} K, {:.1e} A, {:.1e} J",
            oracle[0], oracle[1], oracle[2]
        );
        // Within 3.1e-11 K, 1.6e-11 A and 6.0e-11 J on the export of
        // d4e7a28481; the bisection itself closes to rounding.
        for (value, bound) in oracle.iter().zip([2e-9, 1e-9, 3e-9]) {
            assert!(*value <= bound, "graph 1: {oracle:?}");
        }
        let most = largest(&kernel, &module);
        println!(
            "graph 2 (saturating {saturating}): kernel against module: {:.1e} K, {:.1e} A, {:.1e} J, {:.1e} J/K; end {:.6} K",
            most[0],
            most[1],
            most[3],
            most[4],
            kernel[STEPS - 1].kelvin[0]
        );
        hold(most, "graph 2");
        // The feedback is real: the cell warms by more than 6 K, so the
        // winding's resistivity rises by more than 3 % of itself.
        assert!(kernel[STEPS - 1].kelvin[0] > 306.0, "{kernel:?}");
    }
}

/// Graph 3: the fixture of
/// `finite_capacitor_and_nonlinear_field_use_ordinary_resistor_temperature_and_its_own_heat_cell`.
/// The winding cell holds every triangle, a second cell the material graph
/// resistor; the two share link 17 of 1 W/K.
#[test]
fn the_declared_graph_reproduces_the_material_resistor_winding_and_its_heat_cell() {
    let mut currents = Vec::new();
    for alpha in [0.0, 0.08] {
        let cards: &'static [MaterialCard; 2] = Box::leak(Box::new([
            field::source(true, 18.0),
            graph_card(alpha, 500.0),
        ]));
        let curves: &'static [ResolvedHeatCapacityCurve; 2] =
            Box::leak(Box::new(cards.each_ref().map(caloric)));
        let problem: &'static PlanarNonlinearMagnetostatics =
            Box::leak(Box::new(field::problem(&cards[0], DEPTH)));
        let circuit = native(problem, &cards[0]);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let law = Arc::new(CellConductors::new(&cards[0], &circuit));
        let resistor = Arc::new(BodyResistor::new(&cards[1], &key, &resistors()[0]));

        let module = material_module((cards, curves), &circuit, &law, &resistor, &key);

        let winding = PlanarWinding::at_port_temperatures(
            &circuit,
            vec![law.clone(), law],
            BTreeMap::from([(2, resistor as Arc<dyn ResistorLaw>)]),
            vec![
                HeatReceiver {
                    triangles: (0..8).collect(),
                    resistors: vec![],
                },
                HeatReceiver {
                    triangles: vec![],
                    resistors: vec![2],
                },
            ],
            field::controls(),
            TOLERANCE_J,
        )
        .unwrap();
        let start = two_bodies(curves);
        let (mut system, electrical, cells) = graph(
            winding,
            PlanarWindingState::new(circuit),
            (
                &start,
                &[(MATERIAL, &curves[0]), (RESISTOR_MATERIAL, &curves[1])],
                &[ConductiveLink {
                    link_id: 17,
                    half_distances_m: [0.5; 2],
                    contact_resistance_m2_k_w: 0.0,
                }],
            ),
            &[CELL, RESISTOR_CELL],
            None,
        );
        let kernel = graph_run(&mut system, electrical, cells, &[CELL, RESISTOR_CELL]);
        ledger_closes(&system, electrical, cells);
        let most = largest(&kernel, &module);
        println!(
            "graph 3 (alpha {alpha}): kernel against module: {:.1e} K, {:.1e} A, state {:.1e}, {:.1e} J, {:.1e} J/K; end {:.6} K, {:.6} K",
            most[0],
            most[1],
            most[2],
            most[3],
            most[4],
            kernel[STEPS - 1].kelvin[0],
            kernel[STEPS - 1].kelvin[1]
        );
        hold(most, "graph 3");
        currents.push(kernel[STEPS - 1].current_a);
    }
    // The resistor's R(T) reacts back on the circuit, as in the module's test.
    assert!(currents[1] < currents[0] - 1e-3, "{currents:?}");
}
