//! Conservative finite thermal ports between original shared-wall inventories.
//!
//! Gas/shared-wall solves precede these ordered backward-Euler link stages.
//! Each stage uses the original scalar owner and source storage. This is a
//! first-order composition, not a global graph endpoint or geometric PDE solve.

use super::*;
use crate::multicomponent::wall_exchange::shared::links;
pub(super) use links::WallLinks;
pub use links::{SharedWallLink, SharedWallLinkReport};

impl NetworkWallOperator<'_> {
    /// Connect original shared-wall inventories through finite thermal ports.
    /// Bind after `with_shared_walls`, before any accepted network interval.
    ///
    /// Each gas/shared-wall group first solves its original common endpoint.
    /// The supplied wall links then execute in order, each solving
    /// `Q = dt G (Ta(-Q) - Tb(Q))` on the original h/u or capacity laws.
    /// Each link is implicit and conservative; the full composition is first
    /// order and requires time refinement, including link-order comparison.
    /// The final wall state does not generally obey the earlier gas-stage flux
    /// equation. `shared_wall_link_reports` retains the actual stage endpoints.
    ///
    /// No second wall inventory, global nonlinear solver, material conductivity,
    /// geometry, or spatial-discretization claim is introduced. Resolved spatial
    /// diffusion belongs to the original fs-thermal graph owner. All source
    /// queries spend the existing whole-attempt `SharedWallControl` budget.
    /// Any failed link cancels publication of the entire gas/wall interval.
    pub fn with_shared_wall_links(
        mut self,
        links: Vec<SharedWallLink>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let shared = self
            .shared
            .as_mut()
            .ok_or(Error::Invalid("wall links require shared inventories"))?;
        if self.last_report.is_some()
            || shared.links.is_some()
            || links.is_empty()
            || links.len() > shared.controls.maximum_source_evaluations / 2
        {
            return Err(Error::Invalid("wall link binding state or work budget"));
        }
        let mut reports = Vec::new();
        reports
            .try_reserve_exact(links.len())
            .map_err(|_| Error::Allocation)?;
        for link in &links {
            poll(&mut cancelled)?;
            if link.walls[0] == link.walls[1]
                || link.walls.iter().any(|&i| i >= shared.groups.len())
                || !link.conductance_w_per_k.is_finite()
                || link.conductance_w_per_k < 0.0
            {
                return Err(Error::Invalid("wall link indices or conductance"));
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
            .map_err(|_| Error::Allocation)?;
        candidate.extend_from_slice(&reports);
        shared.links = Some(WallLinks {
            ports: links,
            candidate,
            accepted: reports,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Original authored port order and conductances, unchanged by evolution.
    pub fn shared_wall_links(&self) -> Option<&[SharedWallLink]> {
        Some(&self.shared.as_ref()?.links.as_ref()?.ports)
    }

    /// Last fully accepted link-stage receipts, in the authored link order.
    /// A failed step retains the previous states and these same receipts.
    pub fn shared_wall_link_reports(&self) -> Option<&[SharedWallLinkReport]> {
        let shared = self.shared.as_ref()?;
        shared.last_report?;
        Some(&shared.links.as_ref()?.accepted)
    }
}

pub(super) fn exchange(
    groups: &[Group<'_>],
    link: SharedWallLink,
    duration: Time,
    control: WallControl,
    work: &Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SharedWallLinkReport, WallExchangeError> {
    links::exchange(
        link.walls.map(|i| &groups[i]),
        link,
        duration,
        control,
        work,
        cancelled,
    )
}

#[cfg(test)]
mod tests;
