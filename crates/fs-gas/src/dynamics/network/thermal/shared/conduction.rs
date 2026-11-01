//! Conservative finite thermal ports between original shared-wall inventories.
//!
//! Gas/shared-wall solves precede these ordered backward-Euler link stages.
//! Each stage uses the original scalar owner and source storage. This is a
//! first-order composition, not a global graph endpoint or geometric PDE solve.

use super::*;

/// One stationary, passive thermal port between two original wall inventories.
/// Order is significant at finite timestep: links execute in the supplied
/// order. Parallel interfaces and cycles are admitted without merging them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedWallLink {
    /// Distinct indices in the bound shared-wall array, oriented a then b.
    pub walls: [usize; 2],
    /// Complete authored conductance [W/K], finite and nonnegative. Geometry
    /// or temperature-dependent conductivity is not inferred from this value.
    pub conductance_w_per_k: f64,
}

/// Original states and physical accounts of one accepted link stage.
/// Later links may change these walls again; these temperatures are deliberately
/// retained rather than relabelled as a simultaneous whole-network endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedWallLinkReport {
    /// Actual accepted gas-network interval, used once by each link stage.
    pub elapsed: Time,
    /// Original states entering this particular composition stage.
    pub initial_state: [WallThermal; 2],
    /// Accepted original states leaving this particular composition stage.
    pub endpoint_state: [WallThermal; 2],
    /// Initial temperatures in the link's a,b order [K].
    pub initial_temperature_k: [f64; 2],
    /// Endpoint temperatures in the same order [K].
    pub endpoint_temperature_k: [f64; 2],
    /// Actual represented energy received by wall b [J]; negative reverses it.
    pub heat_a_to_b_j: f64,
    /// Actual original storage changes of wall a and b [J].
    pub wall_energy_change_j: [f64; 2],
    /// Sum of the two actual storage changes [J].
    pub energy_residual_j: f64,
    /// Q - dt G (Ta_new - Tb_new) at this stage [J].
    pub flux_residual_j: f64,
    /// Safeguarded scalar iterations, excluding bracket queries.
    pub iterations: usize,
}

pub(super) struct WallLinks {
    pub(super) ports: Vec<SharedWallLink>,
    pub(super) candidate: Vec<SharedWallLinkReport>,
    pub(super) accepted: Vec<SharedWallLinkReport>,
}

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

// Heat leaving this original wall before either constitutive endpoint. An
// authored constant-capacity wall uses the same passive initial temperature
// hull as the existing gas/wall exchange. Source h/u retains its latent range.
fn heat_limits(
    group: &Group<'_>,
    state: WallThermal,
    other_temperature: f64,
    work: &Work,
) -> Result<[f64; 2], WallExchangeError> {
    work.charge()?;
    if let Some(source) = &group.accepted.storage.heat_storage {
        return source.heat_limits(
            WallContact {
                drag: WallDrag::None,
                thermal: state,
            },
            other_temperature,
        );
    }
    let WallThermal::Lumped {
        temperature_k,
        heat_capacity_j_per_k,
        ..
    } = state
    else {
        return Err(WallExchangeError::Invalid("wall link finite storage"));
    };
    let heat = heat_capacity_j_per_k * (temperature_k - other_temperature);
    if !heat.is_finite() {
        return Err(WallExchangeError::Invalid("wall link storage bound"));
    }
    Ok([heat.min(0.0), heat.max(0.0)])
}

pub(super) fn exchange(
    groups: &[Group<'_>],
    link: SharedWallLink,
    duration: Time,
    control: WallControl,
    work: &Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SharedWallLinkReport, WallExchangeError> {
    poll(cancelled)?;
    let [a, b] = link.walls.map(|i| &groups[i]);
    let states = [a.candidate, b.candidate];
    let initial = [
        a.after_state(states[0], 0.0, work)?,
        b.after_state(states[1], 0.0, work)?,
    ];
    let temperatures = initial.each_ref().map(|p| p.temperature);
    let integrated = duration.value() * link.conductance_w_per_k;
    let initial_flux = integrated * (temperatures[0] - temperatures[1]);
    if !integrated.is_finite() || !initial_flux.is_finite() {
        return Err(WallExchangeError::Invalid(
            "wall link integrated conductance",
        ));
    }
    let mut receipt = SharedWallLinkReport {
        elapsed: duration,
        initial_state: states,
        endpoint_state: states,
        initial_temperature_k: temperatures,
        endpoint_temperature_k: temperatures,
        heat_a_to_b_j: 0.0,
        wall_energy_change_j: [0.0; 2],
        energy_residual_j: 0.0,
        flux_residual_j: -initial_flux,
        iterations: 0,
    };
    if initial_flux.abs() <= control.tolerance(initial_flux.abs()) {
        poll(cancelled)?;
        return Ok(receipt);
    }
    let la = heat_limits(a, states[0], temperatures[1], work)?;
    let lb = heat_limits(b, states[1], temperatures[0], work)?;
    let (mut low, mut high) = if initial_flux > 0.0 {
        (0.0, initial_flux.min(la[1]).min(-lb[0]))
    } else {
        (initial_flux.max(la[0]).max(-lb[1]), 0.0)
    };
    if low < 0.0 {
        low = low.next_up();
    }
    if high > 0.0 {
        high = high.next_down();
    }
    if !low.is_finite() || !high.is_finite() || low >= high {
        return Err(WallExchangeError::Invalid(
            "wall link supported heat interval",
        ));
    }
    let (heat, iterations) = solve_scalar(
        low,
        high,
        0.0,
        control,
        |heat, checkpoint| {
            if checkpoint() {
                return Err(WallExchangeError::Cancelled);
            }
            let wa = a.after_state(states[0], -heat, work)?;
            let wb = b.after_state(states[1], heat, work)?;
            let flux = integrated * (wa.temperature - wb.temperature);
            let residual = heat - flux;
            let slope = wa
                .slope
                .zip(wb.slope)
                .map(|(a, b)| 1.0 + integrated * (a + b));
            if !residual.is_finite() || slope.is_some_and(|s| !positive(s)) {
                return Err(WallExchangeError::Invalid("wall link heat equation"));
            }
            Ok((residual, slope, heat.abs().max(flux.abs())))
        },
        cancelled,
    )?;
    let endpoint = [
        a.after_state(states[0], -heat, work)?,
        b.after_state(states[1], heat, work)?,
    ];
    let actual_heat = endpoint[1].change;
    let energy_residual = endpoint[0].change + actual_heat;
    let flux = integrated * (endpoint[0].temperature - endpoint[1].temperature);
    let flux_residual = actual_heat - flux;
    let scale = endpoint[0]
        .change
        .abs()
        .max(actual_heat.abs())
        .max(flux.abs());
    if !energy_residual.is_finite()
        || !flux_residual.is_finite()
        || energy_residual.abs() > control.tolerance(scale)
        || flux_residual.abs() > control.tolerance(scale)
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j: energy_residual,
            thermal_residual_j: flux_residual,
        });
    }
    receipt.endpoint_state = endpoint.each_ref().map(|p| p.thermal);
    receipt.endpoint_temperature_k = endpoint.each_ref().map(|p| p.temperature);
    receipt.heat_a_to_b_j = actual_heat;
    receipt.wall_energy_change_j = endpoint.each_ref().map(|p| p.change);
    receipt.energy_residual_j = energy_residual;
    receipt.flux_residual_j = flux_residual;
    receipt.iterations = iterations;
    poll(cancelled)?;
    Ok(receipt)
}

#[cfg(test)]
mod tests;
