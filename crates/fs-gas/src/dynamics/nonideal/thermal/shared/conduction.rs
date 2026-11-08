//! Conservative heat exchange between the original finite wall regions.

use super::*;
use crate::multicomponent::wall_exchange::shared::links;
pub(super) use links::WallLinks;
pub use links::{SharedWallLink, SharedWallLinkReport};

impl RealGasWallOperator<'_> {
    /// Connect distinct finite wall regions without duplicating their storage.
    /// Bind after `with_shared_walls`, before any accepted wall interval.
    ///
    /// Each gas/wall group first solves its original common endpoint. These
    /// links then advance in the supplied order, each solving
    /// `Q = dt G (Ta(-Q) - Tb(Q))` on the same nonlinear or latent wall laws.
    /// Parallel links and cycles are permitted. Every stage conserves its two
    /// original inventories; the ordered composition is first order and needs
    /// time refinement. The final wall temperature generally differs from its
    /// earlier gas-contact endpoint, retained in the link-stage reports.
    ///
    /// Conductances are complete authored W/K values; no geometry, material
    /// conductivity or spatial PDE is inferred. All source calls spend the
    /// shared attempt's work budget. Failure of any link prevents publication
    /// of the entire gas/wall interval, including the flow clock and reports.
    pub fn with_shared_wall_links(
        mut self,
        links: Vec<SharedWallLink>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RealGasFlowError> {
        poll(&mut cancelled)?;
        let shared = self.shared.as_mut().ok_or(RealGasFlowError::Invalid(
            "real-gas wall links require shared inventories",
        ))?;
        if self.last_report.is_some()
            || shared.links.is_some()
            || links.is_empty()
            || links.len() > shared.controls.maximum_source_evaluations / 2
        {
            return Err(RealGasFlowError::Invalid(
                "real-gas wall link binding state or work budget",
            ));
        }
        let mut reports = Vec::new();
        reports
            .try_reserve_exact(links.len())
            .map_err(|_| crate::dynamics::Error::Allocation)?;
        for link in &links {
            poll(&mut cancelled)?;
            if link.walls[0] == link.walls[1]
                || link.walls.iter().any(|&i| i >= shared.groups.len())
                || !link.conductance_w_per_k.is_finite()
                || link.conductance_w_per_k < 0.0
            {
                return Err(RealGasFlowError::Invalid(
                    "real-gas wall link indices or conductance",
                ));
            }
            let state = link.walls.map(|i| shared.groups[i].candidate);
            reports.push(SharedWallLinkReport {
                elapsed: Time::new(0.0),
                initial_state: state,
                endpoint_state: state,
                initial_temperature_k: [0.0; 2],
                endpoint_temperature_k: [0.0; 2],
                heat_a_to_b_j: 0.0,
                wall_energy_change_j: [0.0; 2],
                energy_residual_j: 0.0,
                flux_residual_j: 0.0,
                iterations: 0,
            });
        }
        let mut candidate = Vec::new();
        candidate
            .try_reserve_exact(reports.len())
            .map_err(|_| crate::dynamics::Error::Allocation)?;
        candidate.extend_from_slice(&reports);
        shared.links = Some(WallLinks {
            ports: links,
            candidate,
            accepted: reports,
            conduction: None,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Original authored link order and conductances.
    pub fn shared_wall_links(&self) -> Option<&[SharedWallLink]> {
        Some(&self.shared.as_ref()?.links.as_ref()?.ports)
    }

    /// Last accepted stage states and heat balances in the authored link order.
    /// A refusal or cancellation retains the preceding accepted reports.
    /// A thermal piston interval has two source stages: these receipts describe
    /// the second half, while `last_shared_wall_report` covers both halves.
    pub fn shared_wall_link_reports(&self) -> Option<&[SharedWallLinkReport]> {
        let shared = self.shared.as_ref()?;
        shared.last_report?;
        Some(&shared.links.as_ref()?.accepted)
    }
}

pub(super) fn prepare(
    groups: &mut [Group<'_>],
    ports: &mut WallLinks,
    duration: Time,
    controls: WallControl,
    work: &Work,
    total: &mut HeatSum,
    maximum_wall: &mut f64,
    iterations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), RealGasFlowError> {
    for index in 0..ports.ports.len() {
        poll(cancelled)?;
        let link = ports.ports[index];
        let receipt = links::exchange(
            link.walls.map(|i| &groups[i]),
            ports,
            index,
            duration,
            controls,
            work,
            cancelled,
        )?;
        ports.candidate[index] = receipt;
        for (i, &wall) in link.walls.iter().enumerate() {
            groups[wall].candidate = receipt.endpoint_state[i];
        }
        *maximum_wall = maximum_wall.max(receipt.energy_residual_j.abs());
        *iterations =
            iterations
                .checked_add(receipt.iterations)
                .ok_or(RealGasFlowError::Invalid(
                    "real-gas wall link iteration count overflow",
                ))?;
        total.add(WallCellReport {
            elapsed: duration,
            wall_impulse_ns: 0.0,
            mechanical_work_to_gas_j: 0.0,
            drag_heating_j: 0.0,
            heat_to_gas_j: 0.0,
            // Internal transfer cancels; retain the represented balance error
            // instead of adding either directed heat to the gas or a reservoir.
            wall_energy_change_j: receipt.energy_residual_j,
            external_heat_to_gas_j: 0.0,
            energy_residual_j: receipt.energy_residual_j,
            thermal_residual_j: receipt.flux_residual_j,
            iterations: receipt.iterations,
        })?;
    }
    poll(cancelled)?;
    Ok(())
}

