//! The composition-kernel plumbing of the winding graphs: the law that
//! rebinds each triangle's conductor to its cell's card at its port's
//! temperature, the cells as one backward-Euler timeline whose cell ports
//! return end-of-step temperatures, a `PlanarWinding` beside them with one
//! heat route per heat port, and the closure of the kernel's ledger. Each
//! includer's root supplies `DT`.

use std::collections::BTreeMap;
use std::fmt;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_em::planar::eddy::EddyConductor;
use fs_em::winding::component::{ConductorLaw, PlanarWinding, PlanarWindingState, WindingOwner};
use fs_em::winding::ohmic::WindingConductor;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::conductor::ResolvedConductorState;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_system::{NodeId, RouteTemperature, SolveConfig, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, MaterialId, ThermalLaw};

use super::DT;

pub(super) fn card_words(card: &MaterialCard, words: &mut Vec<u64>) {
    for chunk in card.content_hash().0.chunks(8) {
        words.push(u64::from_le_bytes(chunk.try_into().unwrap()));
    }
}

/// Each triangle's conductor in the construction circuit, core and winding
/// alike, rebound to its material at the cell's query point with `kelvin` on
/// the `T` axis, as the modules rebind them. That is the point
/// `ThermalMaterialPoint::at` writes for a cell with no liquid fraction axis.
#[derive(Debug)]
pub(super) struct CellConductors {
    card: MaterialCard,
    context: QueryPoint,
    originals: BTreeMap<usize, ResolvedConductorState>,
}

impl CellConductors {
    pub(super) fn new(
        card: &MaterialCard,
        context: QueryPoint,
        core: &[EddyConductor],
        windings: &[WindingConductor],
    ) -> Self {
        let core = core
            .iter()
            .map(|source| (source.triangle, &source.conductor));
        let windings = windings
            .iter()
            .map(|source| (source.triangle, &source.conductor));
        Self {
            card: card.clone(),
            context,
            originals: core
                .chain(windings)
                .map(|(triangle, conductor)| (triangle, conductor.clone()))
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

pub(super) fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
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

/// The cells of `start` as one implicit conduction graph, joined by
/// `links`.
pub(super) fn conduction(
    start: &DepositionSimulation<'static>,
    materials: &[(MaterialId, &'static ResolvedHeatCapacityCurve)],
    links: &[ConductiveLink],
) -> ImplicitConductionGraph<'static> {
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
    ImplicitConductionGraph::new(
        start,
        &observer,
        &implicit,
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap()
}

/// `conduction` over the cells of `start` as one backward-Euler timeline
/// whose cell ports return end-of-step temperatures.
pub(super) fn timeline(
    start: &DepositionSimulation<'static>,
    conduction: ImplicitConductionGraph<'static>,
) -> (ThermalTimeline<'static>, TimelineState<'static>) {
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

/// `winding` from `state` beside the cells' `timeline` from `cells_start`,
/// heat port `k` routed into `ports[k]`, its source port (if any) held at `drive(t0)` volts
/// over the step that starts at `t0`.
pub(super) fn graph<C, F>(
    winding: PlanarWinding<'static, C, F>,
    state: PlanarWindingState<'static, C>,
    (timeline, cells_start): (ThermalTimeline<'static>, TimelineState<'static>),
    ports: &[CellId],
    drive: Option<fn(f64) -> f64>,
) -> (System, NodeId, NodeId)
where
    C: WindingOwner<'static, FieldControls = F> + 'static,
    F: Copy + fmt::Debug + 'static,
{
    let mut builder = SystemBuilder::new();
    let electrical = builder.add("winding", winding, state).unwrap();
    let cells = builder.add("cells", timeline, cells_start).unwrap();
    if let Some(volts) = drive {
        builder
            .boundary(
                builder.port(electrical, "source0").unwrap(),
                Box::new(move |t0, _, out: &mut [f64]| out.fill(volts(t0))),
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

/// The kernel's ledger closes: the energy the source supplied, less the
/// circuit's stored change and numerical dissipation, plus `environment_j`,
/// the heat the cells' authored boundaries take in, is the cells' stored
/// change; the environment sink gives exactly `environment_j`; the winding
/// and the cells generate the ledger's entropy, the routes carrying only the
/// interface iteration's residue.
pub(super) fn ledger_closes(system: &System, winding: NodeId, cells: NodeId, environment_j: f64) {
    let ledger = system.ledger();
    assert!(
        ledger.energy_residual_j().abs() <= 1e-8,
        "{}",
        ledger.energy_residual_j()
    );
    assert_eq!(ledger.sinks.len(), 1);
    let sink = ledger.sinks[0].heat_j.value();
    assert!(
        (sink + environment_j).abs() <= 1e-12 * (1.0 + environment_j.abs()),
        "{sink} vs {environment_j}"
    );
    let node = &ledger.nodes[winding.index()];
    let supplied: f64 = ledger
        .boundaries
        .iter()
        .map(|b| b.energy_in_j.value())
        .sum();
    let left =
        supplied - node.stored_change_j.value() - node.numerical_loss_j.value() + environment_j;
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
