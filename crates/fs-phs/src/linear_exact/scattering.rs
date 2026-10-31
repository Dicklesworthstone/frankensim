//! One held incident wave at an explicitly selected original physical port.
//! The reference impedance is a boundary relation, never material viscosity.
use super::{LinearSplitBudget, PreparedMultiportScatteringZoh, ScatteringPortRecord};
use crate::{PhsError, PortHamiltonian};

/// Independent physical accounts over the complete prepared duration.
#[derive(Clone, Copy, Debug)]
pub struct ScatteringStepRecord {
    /// Original constitutive storage change [J].
    pub delta_h: f64,
    /// Integral of the ORIGINAL physical resistance, excluding reference loading [J].
    pub physical_heat: f64,
    /// Work through all unselected original ports [J].
    pub external_work: f64,
    /// Net incoming boundary work, including outgoing-wave exchange [J].
    pub boundary_work: f64,
    /// Integral of a²/Z for the held incident wave [J].
    pub incident_energy: f64,
    /// Independently integrated (y-a)²/Z, with no energy subtraction [J].
    pub outgoing_energy: f64,
    /// Interval mean outgoing wave b = y-a, in original port-output units.
    pub outgoing_wave: f64,
}

/// Complete endpoint and ORIGINAL mean port outputs, in their original order.
pub struct ScatteringStep {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub account: ScatteringStepRecord,
}

/// Compatible one-port view of the shared exact multiport scattering owner.
/// Original input/output semantics and independent energy accounts are unchanged.
pub struct PreparedScatteringZoh {
    inner: PreparedMultiportScatteringZoh,
}

impl PreparedScatteringZoh {
    pub fn new_owned(
        system: PortHamiltonian,
        port: usize,
        impedance: f64,
        duration: f64,
        budget: LinearSplitBudget,
    ) -> Result<Self, PhsError> {
        Ok(Self {
            inner: PreparedMultiportScatteringZoh::new_owned(
                system, &[(port, impedance)], duration, budget,
            )?,
        })
    }

    pub fn system(&self) -> &PortHamiltonian { self.inner.system() }
    pub fn duration(&self) -> f64 { self.inner.duration() }
    pub fn selected_port(&self) -> usize { self.inner.ports()[0].0 }

    /// Current-window relation y_bar=free+beta*a in original output units.
    /// The selected input value is ignored; all other authored inputs enter free.
    pub fn average_relation(&self, state: &[f64], inputs: &[f64]) -> Result<(f64, f64), PhsError> {
        let (mut free, mut direct) = ([0.0], [0.0]);
        self.inner.average_relation_into(state, inputs, &mut free, &mut direct)?;
        Ok((free[0], direct[0]))
    }

    /// Allocating complete candidate; accepted source state is never modified.
    pub fn step(&mut self, state: &[f64], inputs: &[f64]) -> Result<ScatteringStep, PhsError> {
        let mut x = vec![0.0; self.system().state_dim()];
        let mut y = vec![0.0; self.system().port_dim()];
        let account = self.step_into(state, inputs, &mut x, &mut y)?;
        Ok(ScatteringStep { x, y, account })
    }

    /// Reuse all prepared scratch; refusal changes neither caller output.
    pub fn step_into(
        &mut self,
        state: &[f64],
        inputs: &[f64],
        next: &mut [f64],
        average_outputs: &mut [f64],
    ) -> Result<ScatteringStepRecord, PhsError> {
        let mut ports = [ScatteringPortRecord::default()];
        let account = self.inner.step_into(state, inputs, next, average_outputs, &mut ports)?;
        Ok(ScatteringStepRecord {
            delta_h: account.delta_h,
            physical_heat: account.physical_heat,
            external_work: account.external_work,
            boundary_work: account.boundary_work,
            incident_energy: account.incident_energy,
            outgoing_energy: account.outgoing_energy,
            outgoing_wave: ports[0].outgoing_wave,
        })
    }
}
