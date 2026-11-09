//! Local physical viscous loss using the SAME positive RK quadrature as flow.
//!
//! Cell i owns pi*A_i*dz and 3*mu_i*(du/dz)_i². This is produced heat,
//! not an evolved caloric inventory. Momentum upwinding and implicit numerical
//! damping never enter it. A consuming thermal owner must also transport the
//! changing cell mass and its energy; a fixed-mass deposition is not equivalent.

use super::*;

/// Borrowed local account from the last successful whole advance call.
#[derive(Clone, Copy, Debug)]
pub struct SlenderJetViscousHeat<'a> {
    /// Nonnegative physical heat generated in each original axial cell [J].
    pub cell_heat_j: &'a [f64],
    /// The identical aggregate returned by SlenderJetAdvance [J].
    pub total_heat_j: f64,
    /// Compensated sum(cell_heat_j) minus total_heat_j [J].
    /// This retains quadrature accumulation roundoff; no cell is corrected.
    pub partition_residual_j: f64,
}

#[derive(Clone, Debug)]
pub(super) struct LocalHeat {
    power: Vec<f64>,
    step: Vec<f64>,
    candidate: Vec<f64>,
    accepted: Vec<f64>,
    total: f64,
    residual: f64,
}
impl LocalHeat {
    fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            power: zeros(n, gate)?,
            step: zeros(n, gate)?,
            candidate: zeros(n, gate)?,
            accepted: zeros(n, gate)?,
            total: 0.0,
            residual: 0.0,
        })
    }
    pub(super) fn set_power(&mut self, i: usize, power_w: f64) {
        self.power[i] = power_w;
    }
}

impl PeriodicSlenderJet {
    /// Prepare local physical-heat accounting once, retaining the native flow.
    ///
    /// Four reusable cell arrays hold current power, a trial RK increment,
    /// the whole-call candidate and the last accepted heat account. Successful
    /// advances allocate nothing here. Errors and cancellation preserve the
    /// preceding account. A domain stop publishes only completed substeps;
    /// a zero-duration successful call publishes zero heat.
    ///
    /// This does not evolve temperature or deposit energy. In a jet the cell
    /// masses change with conservative axial flow; a thermal consumer must
    /// account for advective energy transport as well as this generated heat.
    pub fn enable_cell_viscous_heat(&mut self, gate: &CancelGate) -> Result<(), NonlinearJetError> {
        poll(gate)?;
        if self.spatial.heat.is_none() {
            let heat = LocalHeat::new(self.accepted.a.len(), gate)?;
            poll(gate)?;
            self.spatial.heat = Some(heat);
        }
        Ok(())
    }

    /// Local heat and aggregate closure for the last successful advance.
    /// None means local accounting has not been enabled; zero arrays precede
    /// the first successful advance. Temperature changes do not rewrite history.
    #[must_use]
    pub fn cell_viscous_heat(&self) -> Option<SlenderJetViscousHeat<'_>> {
        self.spatial
            .heat
            .as_ref()
            .map(|heat| SlenderJetViscousHeat {
                cell_heat_j: &heat.accepted,
                total_heat_j: heat.total,
                partition_residual_j: heat.residual,
            })
    }
}

impl Spatial {
    #[cfg(feature = "enthalpy-coupling")]
    pub(super) fn cell_viscous_power_w(&self) -> Result<&[f64], NonlinearJetError> {
        self.heat
            .as_ref()
            .map(|heat| heat.power.as_slice())
            .ok_or(NonlinearJetError::InvalidState(
                "thermal jet heat workspace absent",
            ))
    }

    /// Current trial quadrature only; never part of published heat until accepted.
    #[cfg(feature = "implicit-enthalpy-coupling")]
    pub(super) fn trial_viscous_heat_j(&self) -> Result<&[f64], NonlinearJetError> {
        self.heat
            .as_ref()
            .map(|heat| heat.step.as_slice())
            .ok_or(NonlinearJetError::InvalidState(
                "thermal jet heat workspace absent",
            ))
    }

    pub(super) fn reset_viscous_heat(
        &mut self,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        if let Some(heat) = &mut self.heat {
            for i in 0..heat.candidate.len() {
                poll_at(gate, i)?;
                heat.candidate[i] = 0.0;
            }
        }
        Ok(())
    }

    pub(super) fn start_viscous_heat(
        &mut self,
        weight_s: f64,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        if let Some(heat) = &mut self.heat {
            for i in 0..heat.step.len() {
                poll_at(gate, i)?;
                heat.step[i] = finite(weight_s * heat.power[i])?;
            }
        }
        Ok(())
    }

    pub(super) fn add_viscous_heat(
        &mut self,
        weight_s: f64,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        if let Some(heat) = &mut self.heat {
            for i in 0..heat.step.len() {
                poll_at(gate, i)?;
                heat.step[i] = finite(heat.step[i] + weight_s * heat.power[i])?;
            }
        }
        Ok(())
    }

    pub(super) fn accept_viscous_heat_step(
        &mut self,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        if let Some(heat) = &mut self.heat {
            for i in 0..heat.step.len() {
                poll_at(gate, i)?;
                heat.candidate[i] = finite(heat.candidate[i] + heat.step[i])?;
            }
        }
        Ok(())
    }

    pub(super) fn viscous_heat_partition(
        &self,
        total: f64,
        gate: &CancelGate,
    ) -> Result<f64, NonlinearJetError> {
        let Some(heat) = &self.heat else {
            return Ok(0.0);
        };
        let mut sum = Sum::default();
        for (i, &cell) in heat.candidate.iter().enumerate() {
            poll_at(gate, i)?;
            sum.add(cell);
        }
        finite(sum.total() - total)
    }

    pub(super) fn publish_viscous_heat(&mut self, total: f64, residual: f64) {
        if let Some(heat) = &mut self.heat {
            core::mem::swap(&mut heat.accepted, &mut heat.candidate);
            heat.total = total;
            heat.residual = residual;
        }
    }
}

#[cfg(test)]
mod tests;
