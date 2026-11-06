//! Equal-energy photon recycling on the original thermal observer.
//!
//! Resolved reabsorption creates a new excitation in its receiving cell. All
//! cells evolve simultaneously through the existing positive source propagator;
//! there is no lagged next-step injection, finite generation truncation or
//! steady-state inverse. Photon flight is quasistatic; waiting in an excited
//! state retains its physical lifetime. Broad fluorescence spectra require an
//! authored excitation-dependent redistribution law and are not inferred here.

use super::{
    AbsorbingLuminescenceEnvironment, Cell, CheckedSum, LuminescenceSource,
    LuminescentOpticalHeating, OpticalHeatingError, RecyclingLimits, finite, invalid,
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_thermal::luminescence::ExcitationDecay;
use fs_thermal::luminescence::recycling::{
    RecyclingEdge, RecyclingError, RecyclingStep, ResonantRecyclingPlan,
};
use std::sync::Arc;

#[derive(Debug)]
struct Topology {
    edges: Vec<RecyclingEdge>,
    escape: Vec<f64>,
    unresolved: Vec<f64>,
    fraction_roundoff: f64,
    limits: RecyclingLimits,
}

#[derive(Debug)]
struct Prepared {
    decays: Vec<ExcitationDecay>,
    duration_s: f64,
    plan: ResonantRecyclingPlan,
}

/// Topology and dense interval maps are shared across transactional clones.
/// Only actual accepted stores/accounts are copied with the thermal observer.
#[derive(Clone, Debug)]
pub(super) struct ResonantFeedback {
    topology: Arc<Topology>,
    prepared: Option<Arc<Prepared>>,
}

impl<const N: usize> LuminescentOpticalHeating<N, 1> {
    /// Explicitly declare that every emission-band absorption creates an
    /// excitation of the same sourced emitting state. The single emission line
    /// must equal the stored-state wavelength and lie inside the source's
    /// original excitation domain. No unsupported spectrum is clipped or
    /// broadened. Intrinsic quantum yield and lifetime still govern every
    /// subsequent excitation, including returns to the same optical cell.
    ///
    /// All original material cells become possible emitters, including those
    /// never reached by the external pump, and require directed transport rows.
    /// Call before advancing; switching a running passive history would change
    /// its physical meaning. Quenching is frozen at each native start state.
    /// Numerical work obeys the existing bounded reaction owner's limits.
    pub fn with_resonant_recycling(
        mut self,
        thermal: &DepositionSimulation<'_>,
        limits: RecyclingLimits,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        checkpoint()?;
        let wavelength = self.source.fluorescence().wavelengths_nm()[0];
        let (lo, hi) = self.source.fluorescence().excitation_bounds_nm();
        if self.pending.is_some()
            || self.energy.duration_s != 0.0
            || self.recycling.is_some()
            || wavelength != self.source.excitation_wavelength_nm()
            || !(lo..=hi).contains(&wavelength)
        {
            return Err(invalid(
                "resonant recycling requires an unadvanced supported equal-energy source",
            ));
        }
        let environment = self.emission_transport.as_ref().ok_or_else(|| {
            invalid("resonant recycling requires directed emitted-light transport")
        })?;
        let count = environment.operator().material_cells.len();
        if count == 0 || count > limits.maximum_nodes {
            return Err(invalid(
                "resonant recycling state count exceeds the supplied budget",
            ));
        }
        self.cells.sort_unstable_by_key(|cell| cell.grid_cell);
        let mut additions = Vec::new();
        additions
            .try_reserve_exact(count.saturating_sub(self.cells.len()))
            .map_err(|_| invalid("recycling state allocation refused"))?;
        for (grid_cell, thermal_cell) in environment.bindings() {
            checkpoint()?;
            environment.check_emitter(grid_cell, thermal_cell)?;
            if self
                .cells
                .binary_search_by_key(&grid_cell, |cell| cell.grid_cell)
                .is_err()
            {
                additions.push(Cell {
                    grid_cell,
                    thermal_cell,
                    reservoir_pump_w: [0.0; N],
                    prompt_heat_w: [0.0; N],
                    absorbed_w: [0.0; N],
                    stored_j: 0.0,
                });
            }
        }
        self.cells
            .try_reserve_exact(additions.len())
            .map_err(|_| invalid("recycling state allocation refused"))?;
        self.cells.extend(additions);
        self.cells.sort_unstable_by_key(|cell| cell.grid_cell);
        self.recycling = Some(ResonantFeedback::try_new(
            environment,
            &self.cells,
            limits,
            &mut checkpoint,
        )?);
        self.cell_emission_w = self
            .cells
            .iter()
            .map(|cell| (cell.grid_cell, [0.0]))
            .collect();
        self.mean_cell_emission_w = self.cell_emission_w.clone();
        self.check_state(thermal, true, &mut checkpoint)?;
        Ok(self)
    }
}

impl ResonantFeedback {
    fn try_new<const N: usize>(
        environment: &AbsorbingLuminescenceEnvironment<1>,
        cells: &[Cell<N>],
        limits: RecyclingLimits,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        let mut topology = Topology {
            edges: Vec::new(),
            escape: Vec::with_capacity(cells.len()),
            unresolved: Vec::with_capacity(cells.len()),
            fraction_roundoff: environment.tolerance().fraction_roundoff,
            limits,
        };
        for (source_cell, cell) in cells.iter().enumerate() {
            checkpoint()?;
            let rows = &environment.operator().rows;
            let row_index = rows
                .binary_search_by_key(&cell.grid_cell, |row| row.emitter_cell)
                .map_err(|_| invalid("recycling requires every material's emission row"))?;
            let row = &rows[row_index];
            topology.escape.push(row.escape.mean[0]);
            topology.unresolved.push(row.unresolved.mean[0]);
            for link in &row.absorption {
                checkpoint()?;
                if link.fraction.mean[0] == 0.0 {
                    continue;
                }
                if topology.edges.len() == limits.maximum_edges {
                    return Err(invalid("recycling links exceed the supplied budget"));
                }
                let receiver_cell = cells
                    .binary_search_by_key(&link.cell_index, |cell| cell.grid_cell)
                    .map_err(|_| invalid("reabsorbed photons lack an original excitation store"))?;
                topology
                    .edges
                    .try_reserve(1)
                    .map_err(|_| invalid("recycling link allocation refused"))?;
                topology.edges.push(RecyclingEdge {
                    source_cell,
                    receiver_cell,
                    absorbed_fraction: link.fraction.mean[0],
                });
            }
        }
        Ok(Self {
            topology: Arc::new(topology),
            prepared: None,
        })
    }

    pub(super) fn advance<const N: usize, const M: usize>(
        &mut self,
        source: &LuminescenceSource<M>,
        cells: &[Cell<N>],
        scale: &[f64; N],
        duration_s: f64,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<RecyclingStep, OpticalHeatingError> {
        let mut decays = Vec::with_capacity(cells.len());
        let mut stored = Vec::with_capacity(cells.len());
        let mut pumps = Vec::with_capacity(cells.len());
        for cell in cells {
            checkpoint()?;
            decays.push(source.decay_at_temperature(thermal.temperature_k(cell.thermal_cell)?)?);
            stored.push(cell.stored_j);
            let mut pump = CheckedSum::default();
            for (lane, factor) in scale.iter().enumerate() {
                pump.add(finite(cell.reservoir_pump_w[lane] * factor)?)?;
            }
            pumps.push(pump.value()?);
        }
        if self
            .prepared
            .as_ref()
            .is_none_or(|cached| cached.duration_s != duration_s || cached.decays != decays)
        {
            let topology = &self.topology;
            let plan = preserve_checkpoint(checkpoint, |cancelled| {
                ResonantRecyclingPlan::prepare(
                    &decays,
                    &topology.edges,
                    &topology.escape,
                    &topology.unresolved,
                    duration_s,
                    topology.fraction_roundoff,
                    topology.limits,
                    cancelled,
                )
            })?;
            self.prepared = Some(Arc::new(Prepared {
                decays,
                duration_s,
                plan,
            }));
        }
        let prepared = self.prepared.as_ref().expect("prepared recycling interval");
        preserve_checkpoint(checkpoint, |cancelled| {
            prepared
                .plan
                .advance_with_cancel(&stored, &pumps, cancelled)
        })
    }
}

fn preserve_checkpoint<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, RecyclingError>,
) -> Result<T, OpticalHeatingError> {
    let mut refused = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            refused = Some(error);
            true
        }
    });
    if let Some(error) = refused {
        return Err(error.into());
    }
    result.map_err(|_| invalid("coupled resonant photon recycling refused"))
}
