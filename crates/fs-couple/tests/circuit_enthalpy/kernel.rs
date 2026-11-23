//! The kernel side shared by the declared-graph comparisons: the fixture's
//! two cells (`common`) as one backward-Euler timeline whose cell ports
//! return end-of-step temperatures, and the execution context the kernel
//! steps in.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, MaterialId, ThermalLaw};

use crate::common::sourced;

/// The fixture's cells, cell 30 then cell 10.
pub const CELLS: [CellId; 2] = [CellId(30), CellId(10)];

pub fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 11,
            kernel_id: 29,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// The specific enthalpy [J/kg] at `kelvin` on a caloric curve whose knots
/// rise in temperature, by the curve's own linear interpolation.
pub fn enthalpy_at(curve: &ResolvedHeatCapacityCurve, kelvin: f64) -> Result<f64, String> {
    let knots = curve.curve().knots();
    for pair in knots.windows(2) {
        let (low, high) = (&pair[0], &pair[1]);
        if (low.temperature_k..=high.temperature_k).contains(&kelvin)
            && high.temperature_k > low.temperature_k
        {
            let fraction = (kelvin - low.temperature_k) / (high.temperature_k - low.temperature_k);
            return Ok(low.specific_enthalpy_j_kg
                + fraction * (high.specific_enthalpy_j_kg - low.specific_enthalpy_j_kg));
        }
    }
    Err(format!("{kelvin} K is outside the caloric curve"))
}

/// The fixture's caloric curve with the card's constant conductivity `k`,
/// for the implicit conduction graph. The temperature is the curve's own.
struct CardConduction {
    curve: &'static ResolvedHeatCapacityCurve,
    k_w_m_k: f64,
}

impl ConductivityLaw for CardConduction {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.k_w_m_k)
    }
}

impl ImplicitConductionLaw for CardConduction {
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
        // The fixture's curve is one linear segment: dT/dh is its inverse slope.
        let slope = (last.temperature_k - first.temperature_k)
            / (last.specific_enthalpy_j_kg - first.specific_enthalpy_j_kg);
        Ok(ImplicitConductionState {
            temperature_k: self.curve.curve().temperature_k(h)?,
            temperature_slope_kg_k_j: slope,
            conductivity_w_m_k: self.k_w_m_k,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// The fixture's two cells as one backward-Euler timeline whose cell ports
/// return end-of-step temperatures, and its start, each cell's boundary
/// carrying its source of `sources_w` (the timeline books that heat on its
/// `environment` port). The link's conductance
/// comes from the cards' conductivity `k` across two half-cells of 0.5 m and
/// the link's 1 m^2: `1 / (0.5 / k + 0.5 / k)`, the fixture's authored
/// conductance when `k` is that conductance.
pub fn timeline(
    curves: &'static [ResolvedHeatCapacityCurve; 2],
    conductance_w_k: f64,
    sources_w: [f64; 2],
    dt_s: f64,
) -> (ThermalTimeline<'static>, TimelineState<'static>) {
    let start: DepositionSimulation<'static> = sourced(curves, conductance_w_k, None, sources_w);
    let laws: &'static [CardConduction; 2] =
        Box::leak(Box::new(core::array::from_fn(|i| CardConduction {
            curve: &curves[i],
            k_w_m_k: conductance_w_k,
        })));
    let observer = ConductionObserver::new(
        &start,
        &[
            (MaterialId(0), &laws[0] as &dyn ConductivityLaw),
            (MaterialId(1), &laws[1]),
        ],
        &[ConductiveLink {
            link_id: 7,
            half_distances_m: [0.5; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let conduction = ImplicitConductionGraph::new(
        &start,
        &observer,
        &[
            (MaterialId(0), &laws[0] as &dyn ImplicitConductionLaw),
            (MaterialId(1), &laws[1]),
        ],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let controls = ImplicitConductionControls {
        maximum_dt_s: dt_s,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-9,
        ..Default::default()
    };
    let timeline = ThermalTimeline::implicit(
        CELLS.to_vec(),
        conduction,
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .unwrap()
    .returning_end_temperatures()
    .unwrap();
    let cells_start = timeline.start(&start).unwrap();
    (timeline, cells_start)
}
