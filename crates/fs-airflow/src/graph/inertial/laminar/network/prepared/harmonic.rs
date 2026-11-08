//! Continuous pressure/head response of the already projected laminar network.
//!
//! The original port-Hamiltonian coordinates retain every independent mean-flow
//! mode and every pipe's viscous-memory arm. Harmonic queries reuse that system
//! directly under `exp(-i omega t)`; no steady resistance network or isolated
//! impedance assembly replaces it. Pressure nodes are gauge-normalized exactly
//! as in time-domain stepping. Outputs are pressure-node flows first, followed by
//! one original branch flow per explicit head port.

mod audit;
mod batch;
mod settled;
#[cfg(test)]
mod tests;

use super::*;
use fs_math::c64::C64;
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel, HarmonicResponse};

impl<M: FluidTransport> LaminarNetwork<M> {
    /// Prepare the complete native continuous-frequency model once. External
    /// composition may reuse it, but that low-level path owns every original
    /// frequency-band, peak-Reynolds and continuity admission itself. The
    /// checked response methods below perform all of those admissions.
    pub fn harmonic_model(
        &self,
        cx: &Cx<'_>,
        budget: HarmonicBudget,
    ) -> Result<HarmonicModel, HarmonicError> {
        poll(cx)?;
        let model = HarmonicModel::from_ode(&self.system, budget, || cx.checkpoint().is_err())?;
        poll(cx)?;
        Ok(model)
    }

    /// Gauge-normalize reservoir pressure phasors and append one signed pressure
    /// rise per branch. Positive heads drive each authored branch direction.
    /// The returned port order exactly matches [`Self::system`].
    pub fn harmonic_input_phasors(
        &self,
        reservoir_pressures_pa: &[C64],
        branch_heads_pa: Option<&[C64]>,
    ) -> Result<Vec<C64>, HarmonicError> {
        self.admit_harmonic_input(reservoir_pressures_pa, branch_heads_pa)?;
        let boundary = self.topology.pressure_nodes.len();
        let edges = self.branches.len();
        let mut input = vec![C64::ZERO; boundary + edges];
        for port in 0..boundary {
            input[port] = reservoir_pressures_pa[port]
                - reservoir_pressures_pa[self.basis.pressure_reference[port]];
        }
        for edge in 0..edges {
            input[boundary + edge] = branch_heads_pa.map_or(C64::ZERO, |heads| heads[edge]);
        }
        if input.iter().any(|value| !finite(*value)) {
            return Err(invalid(
                "gauge-normalized pressure phasors are not representable",
            ));
        }
        Ok(input)
    }

    /// Solve one reservoir-pressure excitation with zero explicit branch heads.
    ///
    /// The solved peak flow in every original pipe is admitted against that
    /// pipe's Reynolds ceiling. Aggregate storage and heat are independently
    /// reconstructed from the original per-pipe coordinates and source squares.
    pub fn harmonic_response(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressures_pa: &[C64],
        budget: HarmonicBudget,
    ) -> Result<HarmonicResponse, HarmonicError> {
        self.harmonic_response_impl(cx, omega_rad_s, reservoir_pressures_pa, None, budget)
    }

    /// Solve one reservoir-pressure and explicit branch-head excitation.
    pub fn harmonic_response_with_heads(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressures_pa: &[C64],
        branch_heads_pa: &[C64],
        budget: HarmonicBudget,
    ) -> Result<HarmonicResponse, HarmonicError> {
        self.harmonic_response_impl(
            cx,
            omega_rad_s,
            reservoir_pressures_pa,
            Some(branch_heads_pa),
            budget,
        )
    }

    fn harmonic_response_impl(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressures_pa: &[C64],
        branch_heads_pa: Option<&[C64]>,
        budget: HarmonicBudget,
    ) -> Result<HarmonicResponse, HarmonicError> {
        poll(cx)?;
        self.admit_harmonic_input(reservoir_pressures_pa, branch_heads_pa)?;
        self.admit_harmonic_frequency(cx, omega_rad_s)?;
        let input = self.harmonic_input_phasors(reservoir_pressures_pa, branch_heads_pa)?;
        let tolerance = harmonic_tolerance(self, input.len(), budget)?;
        let model = self.harmonic_model(cx, budget)?;
        let factor = model.factor(omega_rad_s, || cx.checkpoint().is_err())?;
        self.response_with_factor(cx, &factor, &input, tolerance)
    }

    fn admit_harmonic_input(
        &self,
        reservoir_pressures_pa: &[C64],
        branch_heads_pa: Option<&[C64]>,
    ) -> Result<(), HarmonicError> {
        let boundary = self.topology.pressure_nodes.len();
        let edges = self.branches.len();
        if reservoir_pressures_pa.len() != boundary
            || branch_heads_pa.is_some_and(|heads| heads.len() != edges)
            || reservoir_pressures_pa.iter().any(|value| !finite(*value))
            || branch_heads_pa.is_some_and(|heads| heads.iter().any(|value| !finite(*value)))
        {
            return Err(invalid(
                "complete finite reservoir-pressure and branch-head phasors required",
            ));
        }
        Ok(())
    }

    fn admit_harmonic_frequency(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
    ) -> Result<(), HarmonicError> {
        poll(cx)?;
        for pipe in &self.branches {
            pipe.impedance_at(omega_rad_s).map_err(|_| {
                invalid("frequency lies outside an original pipe wall-memory admission")
            })?;
            poll(cx)?;
        }
        Ok(())
    }
}

fn harmonic_tolerance<M>(
    network: &LaminarNetwork<M>,
    ports: usize,
    budget: HarmonicBudget,
) -> Result<f64, HarmonicError> {
    let tolerance = budget.relative_residual_tolerance
        + ROUNDOFF * network.system.state_dim().max(ports) as f64;
    if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0 {
        return Err(invalid("harmonic network residual tolerance is invalid"));
    }
    Ok(tolerance)
}

fn poll(cx: &Cx<'_>) -> Result<(), HarmonicError> {
    cx.checkpoint().map_err(|_| HarmonicError::Cancelled)
}
fn invalid(what: &'static str) -> HarmonicError {
    HarmonicError::Invalid { what }
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
