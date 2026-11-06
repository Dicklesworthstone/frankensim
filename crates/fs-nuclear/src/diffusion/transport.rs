//! Explicit finite-volume transport of delayed-neutron precursor counts.
//!
//! A directed volumetric flow Q carries each family from its upstream cell at
//! rate Q/V_upstream. Transfers conserve tracked precursor count; an outflow
//! removes it to a separately reported exterior. A precursor continues to
//! decay while moving through the admitted network, emitting its neutron in
//! the cell where it decays, not necessarily where it was produced.
//!
//! These are prescribed tracer flows through fixed control volumes. They do
//! not evolve bulk fluid mass, dissolved fuel, temperature or enthalpy, infer
//! unlisted return paths, or solve the carrier's momentum/continuity equations.

use super::{
    DiffusionData, DiffusionError, DiffusionLimits, FissionError, nonnegative, quotient, reserved,
};

/// Expected precursor transfer from one well-mixed cell to another.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrecursorTransfer {
    pub from_cell: usize,
    pub to_cell: usize,
    pub volumetric_flow_m3_s: f64,
}

/// Removal to an untracked exterior. Future exterior decays are not credited
/// to any retained cell or treated as a neutron leakage event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrecursorOutflow {
    pub cell: usize,
    pub volumetric_flow_m3_s: f64,
}

/// Ordered authored flows; empty transport keeps every precursor stationary.
/// Repeated endpoints represent distinct additive flow paths.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PrecursorTransport {
    pub transfers: Vec<PrecursorTransfer>,
    pub outflows: Vec<PrecursorOutflow>,
}

pub(super) struct TransportRates {
    pub transfers: Vec<f64>,
    pub outflows: Vec<f64>,
}

impl PrecursorTransport {
    pub(super) fn check_limits(
        &self,
        data: &DiffusionData,
        limits: DiffusionLimits,
    ) -> Result<(), DiffusionError> {
        let transfers =
            self.transfers
                .len()
                .checked_add(data.links.len())
                .ok_or(DiffusionError::Budget(
                    "combined diffusion and precursor interfaces",
                ))?;
        let outflows =
            self.outflows
                .len()
                .checked_add(data.leakage.len())
                .ok_or(DiffusionError::Budget(
                    "combined diffusion and precursor boundaries",
                ))?;
        if transfers > limits.max_links {
            return Err(DiffusionError::Budget(
                "combined diffusion and precursor interfaces",
            ));
        }
        if outflows > limits.max_leakage_boundaries {
            return Err(DiffusionError::Budget(
                "combined diffusion and precursor boundaries",
            ));
        }
        Ok(())
    }

    pub(super) fn rates(
        &self,
        data: &DiffusionData,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<TransportRates, DiffusionError> {
        let volumes = &data.cell_volumes_m3;
        let mut transfers = reserved(self.transfers.len())?;
        for transfer in &self.transfers {
            checkpoint()?;
            if transfer.from_cell >= volumes.len()
                || transfer.to_cell >= volumes.len()
                || transfer.from_cell == transfer.to_cell
                || !nonnegative(transfer.volumetric_flow_m3_s)
            {
                return Err(DiffusionError::Invalid(
                    "invalid directed precursor transfer",
                ));
            }
            transfers.push(quotient(
                transfer.volumetric_flow_m3_s,
                volumes[transfer.from_cell],
            )?);
        }
        let mut outflows = reserved(self.outflows.len())?;
        for outflow in &self.outflows {
            checkpoint()?;
            if outflow.cell >= volumes.len() || !nonnegative(outflow.volumetric_flow_m3_s) {
                return Err(DiffusionError::Invalid("invalid precursor outflow"));
            }
            outflows.push(quotient(
                outflow.volumetric_flow_m3_s,
                volumes[outflow.cell],
            )?);
        }
        Ok(TransportRates {
            transfers,
            outflows,
        })
    }

    pub(super) fn same_topology(&self, other: &Self) -> bool {
        self.transfers.len() == other.transfers.len()
            && self
                .transfers
                .iter()
                .zip(&other.transfers)
                .all(|(a, b)| (a.from_cell, a.to_cell) == (b.from_cell, b.to_cell))
            && self.outflows.len() == other.outflows.len()
            && self
                .outflows
                .iter()
                .zip(&other.outflows)
                .all(|(a, b)| a.cell == b.cell)
    }
}

#[cfg(test)]
mod tests;
