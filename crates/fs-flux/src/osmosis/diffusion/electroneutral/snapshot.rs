//! Accepted-inventory electrical observations through the original native solve.

use super::{ElectroneutralDiffusionReport, OsmoticNetwork, OsmoticNetworkError, filled};
use crate::osmosis::poll;
use fs_material::osmosis::diffusion::electrical::PrescribedPotentialDiffusionReport;

/// One active molecular connection in its original left-to-right orientation.
/// Fluxes are in the network's declared original solute order; zero-area links
/// are omitted and can be identified through their original network indices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonicDiffusionConnection {
    pub link_index: usize,
    pub endpoints: [usize; 2],
    pub area_m2: f64,
    pub distance_m: f64,
    /// Native force/friction, current, voltage and entropy diagnostics.
    pub transport: PrescribedPotentialDiffusionReport,
}

/// A newly solved physical electrical state at the accepted network time.
/// This is not the last coarse, fine or Newton probe from a prior advance.
#[derive(Clone, Debug)]
pub struct IonicDiffusionSnapshot {
    pub time_s: f64,
    pub temperature_k: f64,
    pub species_count: usize,
    /// Work belongs only to this readout; it consumes no later advance's budget.
    pub work: ElectroneutralDiffusionReport,
    potentials: Vec<f64>,
    gauges: Vec<usize>,
    charge_rates: Vec<f64>,
    connections: Vec<IonicDiffusionConnection>,
    fluxes: Vec<f64>,
}

impl IonicDiffusionSnapshot {
    /// `F*phi/(RT)` per compartment, with one zero reference per conducting
    /// component. Absolute differences BETWEEN disconnected components have
    /// no physical meaning. Native connection reports also carry volts.
    pub fn dimensionless_potentials(&self) -> &[f64] {
        &self.potentials
    }

    /// Reference-compartment index for each node's conducting component.
    /// An isolated or exactly carrier-free node supplies its own reference.
    pub fn potential_gauges(&self) -> &[usize] {
        &self.gauges
    }

    /// Actual signed net OUTGOING elementary-charge equivalents [mol/s].
    /// Every entry passed the configured absolute node-current budget.
    pub fn node_charge_rates_mol_s(&self) -> &[f64] {
        &self.charge_rates
    }

    pub fn connections(&self) -> &[IonicDiffusionConnection] {
        &self.connections
    }

    /// Native stationary-solvent species fluxes [mol/(m² s)] for a connection
    /// in `connections()`, positive from its left endpoint to its right.
    pub fn solute_fluxes_mol_m2_s(&self, connection: usize) -> Option<&[f64]> {
        if connection >= self.connections.len() {
            return None;
        }
        let start = connection * self.species_count;
        Some(&self.fluxes[start..start + self.species_count])
    }
}

impl OsmoticNetwork {
    /// Solve and observe ionic transport at CURRENT accepted inventories.
    /// Uses the same original material/source responses, gauge-fixed native CG
    /// and combined-force replays as evolution, with fresh source qualification.
    /// No amount, solvent volume, accepted time or advance work counter changes.
    /// Works for both explicit and implicit global electrical integration.
    ///
    /// `maximum_snapshot_bytes` bounds additional retained output-array payload,
    /// beyond the already admitted network numerical workspace. Source-query
    /// receipt strings and allocator overhead retain their existing exclusions.
    /// The global electrical controls also bound this independent solve's work.
    /// Cancellation or any source/current refusal returns no partial snapshot.
    pub fn ionic_diffusion_snapshot(
        &mut self,
        maximum_snapshot_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonicDiffusionSnapshot, OsmoticNetworkError> {
        poll(&mut cancelled)?;
        let physical = self.diffusion_controls.ok_or(OsmoticNetworkError::Invalid(
            "ionic snapshot requires prepared molecular transport",
        ))?;
        let electrical =
            self.electroneutral_diffusion
                .as_mut()
                .ok_or(OsmoticNetworkError::Invalid(
                    "ionic snapshot requires the global potential/current closure",
                ))?;
        let nodes = self.cells.len();
        let species = self.species.len();
        let active = self
            .diffusion
            .iter()
            .filter(|link| link.area_m2 != 0.0)
            .count();
        let entries = active
            .checked_mul(species)
            .ok_or(OsmoticNetworkError::Invalid("ionic snapshot size overflow"))?;
        let required = nodes
            .checked_mul(2 * core::mem::size_of::<f64>() + core::mem::size_of::<usize>())
            .and_then(|bytes| {
                bytes.checked_add(
                    active.checked_mul(core::mem::size_of::<IonicDiffusionConnection>())?,
                )
            })
            .and_then(|bytes| bytes.checked_add(entries.checked_mul(core::mem::size_of::<f64>())?))
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<IonicDiffusionSnapshot>()))
            .ok_or(OsmoticNetworkError::Invalid("ionic snapshot size overflow"))?;
        if required > maximum_snapshot_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_snapshot_bytes,
            });
        }
        let mut snapshot = IonicDiffusionSnapshot {
            time_s: self.time_s,
            temperature_k: self.cells[0].solution.temperature_k(),
            species_count: species,
            work: ElectroneutralDiffusionReport::default(),
            potentials: filled(nodes, 0.0)?,
            gauges: filled(nodes, 0)?,
            charge_rates: filled(nodes, 0.0)?,
            connections: Vec::new(),
            fluxes: Vec::new(),
        };
        snapshot
            .connections
            .try_reserve_exact(active)
            .map_err(|_| {
                OsmoticNetworkError::Invalid("ionic connection readout allocation failed")
            })?;
        snapshot
            .fluxes
            .try_reserve_exact(entries)
            .map_err(|_| OsmoticNetworkError::Invalid("ionic flux readout allocation failed"))?;
        // Restored on BOTH ordinary success and failure; numerical scratch may
        // change, but the following advance always prepares its own fresh state.
        let prior_report = electrical.report;
        electrical.begin_window();
        let result = (|| {
            let mut evaluations = 0;
            electrical.evaluate_transfers_observed(
                &self.cells,
                &mut self.diffusion,
                1.0,
                physical,
                &mut evaluations,
                &mut cancelled,
                &mut |link_index, endpoints, area_m2, distance_m, flux, transport| {
                    snapshot.connections.push(IonicDiffusionConnection {
                        link_index,
                        endpoints,
                        area_m2,
                        distance_m,
                        transport,
                    });
                    snapshot.fluxes.extend_from_slice(flux);
                    Ok(())
                },
            )?;
            snapshot.potentials.copy_from_slice(&electrical.potentials);
            snapshot.gauges.copy_from_slice(&electrical.parents);
            snapshot
                .charge_rates
                .copy_from_slice(&electrical.charge_rate);
            snapshot.work = electrical.report;
            poll(&mut cancelled)?;
            Ok(snapshot)
        })();
        electrical.report = prior_report;
        result
    }
}
