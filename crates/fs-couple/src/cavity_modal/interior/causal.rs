//! Causal exterior reaction on the complete original acoustic pressure field.
//!
//! The FE owner advances every free pressure/potential coordinate exactly for
//! held inward velocities. The existing rational or retarded MIMO owner advances
//! its original states exactly for held OUTWARD velocities. Matching the two
//! interval-mean conjugate forces closes the same-interval feedback, including
//! every mutual term. No pressure-release substitute, lagged boundary force,
//! modal truncation or second exterior integrator is introduced.
//!
//! The impedance uses native velocity/force coordinates: its input is the
//! outward amplitude [m/s] of each declared dimensionless opening shape and
//! its output is integral(shape * pressure) dA [N]. Thus Z has units N s/m,
//! just like the native BEM aperture matrix. Signed and zero-net-volume modes
//! need no division by area or signed volume gain. Supplied rational models
//! already use exp(+i omega t); outgoing exp(-i omega t) BEM data must have
//! been conjugated once when that SAME source realization was constructed.
//!
//! This module interconnects admitted source owners; it does not fit a model or
//! establish that a supplied model describes the FE openings, exterior gas or
//! physical geometry. Those are caller obligations, including source-band,
//! spatial, aperture-mode, realization and exchange-clock convergence. The
//! native continuous/sampled passivity reports retain their original classes;
//! GridOnly is never promoted to a storage certificate. Exterior exchange work
//! is signed boundary work, not inferred radiation, irreversible heat or a
//! separately known exterior energy. A returning wave may make it negative.

use fs_feec::acoustics::helmholtz::heterogeneous::time::exact::{
    ExactTimeLimits, ExactVelocityPort, PreparedExactTime,
};
use fs_feec::acoustics::helmholtz::heterogeneous::time::{
    HeterogeneousTimeState, ResistiveBoundary,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{
    HeterogeneousError, HeterogeneousHelmholtzFe, VolumeVelocitySource,
};
use fs_la::factor::{FactorError, Lu, lu};
use fs_math::dd::Dd;
use fs_vfit::mimo::MimoPassivityReport;
use fs_vfit::mimo::delayed::PreparedDelayedMimoLoad;
use fs_vfit::mimo::load::{MimoLoadError, PreparedMimoLoad};
use std::collections::BTreeSet;

/// The original admitted impedance and accepted history, without a new model
/// conversion. Both variants use exact held-input interval-mean outputs.
#[derive(Clone, Debug)]
pub enum CausalExteriorLoad {
    Rational(PreparedMimoLoad),
    /// Fractional arrival times retain the native two-interval analytic rule;
    /// integer delays are not rounded and pole/history states are not dropped.
    Delayed(PreparedDelayedMimoLoad),
}

impl CausalExteriorLoad {
    #[must_use]
    pub fn ports(&self) -> usize {
        match self {
            Self::Rational(load) => load.ports(),
            Self::Delayed(load) => load.ports(),
        }
    }

    #[must_use]
    pub fn duration_s(&self) -> f64 {
        match self {
            Self::Rational(load) => load.time_step_s(),
            Self::Delayed(load) => load.time_step_s(),
        }
    }

    /// Unmodified continuous and sampled admission reports, in that order.
    #[must_use]
    pub fn reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        match self {
            Self::Rational(load) => load.reports(),
            Self::Delayed(load) => load.reports(),
        }
    }

    fn mean_direct_term(&self) -> &[f64] {
        match self {
            Self::Rational(load) => load.mean_direct_term(),
            Self::Delayed(load) => load.mean_direct_term(),
        }
    }

    fn mean_history_into(&self, effort: &mut [f64]) -> Result<(), MimoLoadError> {
        match self {
            Self::Rational(load) => load.mean_history_output_into(effort),
            Self::Delayed(load) => load.mean_history_output_into(effort),
        }
    }

    fn step_mean_into(&mut self, flow: &[f64], effort: &mut [f64]) -> Result<(), MimoLoadError> {
        match self {
            Self::Rational(load) => load.step_mean_into(flow, effort),
            Self::Delayed(load) => load.step_mean_into(flow, effort),
        }
    }
}

/// The FE exponential and source runtime keep their own original budgets.
#[derive(Clone, Copy, Debug)]
pub struct CausalTimeLimits {
    pub interior: ExactTimeLimits,
    /// Explicit complete opening-mode count, in 1..=128. No mode is dropped.
    pub maximum_exterior_ports: usize,
    /// Relative original mean-force constraint residual, strictly in (0,1).
    /// This is a numerical closure tolerance, not a coupled-model error bound.
    pub interface_relative_tolerance: f64,
}

/// Refusal leaves the accepted interior field and every exterior history intact.
#[derive(Debug)]
pub enum CausalTimeError {
    Invalid(&'static str),
    Interior(HeterogeneousError),
    Exterior(MimoLoadError),
    Closure(FactorError),
    InterfaceResidual { relative: f64, tolerance: f64 },
    Cancelled,
}
impl core::fmt::Display for CausalTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "causal interior/exterior acoustics: {self:?}")
    }
}
impl std::error::Error for CausalTimeError {}
impl From<HeterogeneousError> for CausalTimeError {
    fn from(error: HeterogeneousError) -> Self {
        Self::Interior(error)
    }
}
impl From<MimoLoadError> for CausalTimeError {
    fn from(error: MimoLoadError) -> Self {
        Self::Exterior(error)
    }
}
fn invalid(what: &'static str) -> CausalTimeError {
    CausalTimeError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CausalTimeError> {
    if cancelled() {
        Err(CausalTimeError::Cancelled)
    } else {
        Ok(())
    }
}
fn face_key(mut face: [u32; 3]) -> [u32; 3] {
    face.sort_unstable();
    face
}

/// Same-trajectory physical fluid accounts, with exterior transfer counted once.
#[derive(Clone, Copy, Debug)]
pub struct CausalTimeWork {
    pub initial_fluid_energy_j: f64,
    pub final_fluid_energy_j: f64,
    pub delta_fluid_energy_j: f64,
    /// Original FE integral of p^T C p; never an inferred energy deficit.
    pub boundary_dissipation_j: f64,
    /// Signed prescribed feed work plus original volume-source work [J].
    pub source_work_j: f64,
    /// Signed prescribed velocity-feed subtotal, already in source_work_j.
    pub boundary_drive_work_j: f64,
    /// Signed volume-source subtotal, already in source_work_j.
    pub volume_source_work_j: f64,
    /// Work delivered to the exterior by its held velocities and exact mean
    /// force. This includes reversible storage and returning-wave contributions.
    pub exterior_exchange_work_j: f64,
    /// Exterior work minus work leaving the fluid at the same opening modes.
    /// This is the measured coupling defect, never a loss term.
    pub interface_work_defect_j: f64,
}
impl CausalTimeWork {
    #[must_use]
    pub fn balance_residual_j(self) -> f64 {
        self.delta_fluid_energy_j + self.boundary_dissipation_j + self.exterior_exchange_work_j
            - self.source_work_j
    }
}

/// Accepted same-interval opening motion and the two independent force readings.
#[derive(Clone, Debug)]
pub struct CausalTimeStep {
    /// Native modal velocity coordinates, in exterior_ports() order [m/s].
    /// A signed shape's coordinate is not its net volume flow.
    pub outward_velocity_m_s: Vec<f64>,
    /// Original FE pressure integral divided by this interval's duration [N].
    pub interior_mean_force_n: Vec<f64>,
    /// Original rational/delayed owner's exact mean output [N].
    pub exterior_mean_force_n: Vec<f64>,
    pub interface_relative_residual: f64,
    pub work: CausalTimeWork,
}

/// Full original pressure state and original exterior histories on one clock.
/// Preparation caches one small port factor; warm steps have no factorization,
/// BEM/material query, exponential, delay interpolation or history cloning.
pub struct CausalInteriorTime<'m, 'f, 'c> {
    prepared: PreparedExactTime<'m, 'f, 'c>,
    state: HeterogeneousTimeState<'m, 'f, 'c>,
    exterior: CausalExteriorLoad,
    exterior_ports: Vec<usize>,
    driven_ports: Vec<usize>,
    /// H maps held inward velocities to original interval-mean FE forces.
    interior_direct: Vec<f64>,
    /// (H + D_mean) maps outward velocity to unforced FE force minus history.
    closure: Lu,
    limits: CausalTimeLimits,
    accepted_steps: u64,
}

impl<'m, 'f, 'c> CausalInteriorTime<'m, 'f, 'c> {
    /// Bind complete original spatial ports to an already admitted exterior
    /// realization. `ports` contains dimensionless P1 normal-velocity shapes;
    /// all FE pressure coordinates, of any admitted order, remain present.
    /// `exterior_ports` maps the source matrix's rows/columns to those shapes in
    /// source order. All remaining ports are prescribed feeds in original order.
    /// Signed and overlapping exterior modes on the same opening are allowed.
    /// A feed or local resistance cannot share an exterior opening face, and
    /// release faces retain the original exact owner's stricter admission.
    ///
    /// The existing source band must already cover this clock's Nyquist limit.
    /// Original load history is retained, even when nonzero. The fluid starts
    /// at zero; initialize_field permits an explicit coherent initial pair.
    /// No original source state is silently reset. Spatial, realization and
    /// same-interval held-flow exchange errors require independent ladders.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
        duration_s: f64,
        boundaries: &[ResistiveBoundary],
        ports: &[ExactVelocityPort],
        exterior_ports: &[usize],
        volume_sources: &[&VolumeVelocitySource<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        exterior: CausalExteriorLoad,
        limits: CausalTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CausalTimeError> {
        poll(cancelled)?;
        let n = exterior_ports.len();
        if !(1..=128).contains(&limits.maximum_exterior_ports)
            || n == 0
            || n > limits.maximum_exterior_ports
            || n != exterior.ports()
            || ports.len() > 128
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || duration_s.to_bits() != exterior.duration_s().to_bits()
            || !limits.interface_relative_tolerance.is_finite()
            || !(0.0 < limits.interface_relative_tolerance
                && limits.interface_relative_tolerance < 1.0)
        {
            return Err(invalid(
                "complete bounded source ports, matching physical clocks and a finite closure tolerance are required",
            ));
        }
        let mut selected = vec![false; ports.len()];
        for &index in exterior_ports {
            if index >= ports.len() || selected[index] {
                return Err(invalid(
                    "exterior indices must uniquely identify original velocity ports",
                ));
            }
            if ports[index]
                .profile
                .iter()
                .all(|(_, values)| values.iter().all(|&v| v == 0.0))
            {
                return Err(invalid(
                    "each exterior mode needs a nonzero original physical velocity shape",
                ));
            }
            selected[index] = true;
        }
        let driven_ports: Vec<_> = (0..ports.len()).filter(|&index| !selected[index]).collect();
        let exterior_faces: BTreeSet<_> = exterior_ports
            .iter()
            .flat_map(|&index| ports[index].profile.iter().map(|(face, _)| face_key(*face)))
            .collect();
        for &index in &driven_ports {
            poll(cancelled)?;
            if ports[index]
                .profile
                .iter()
                .any(|(face, _)| exterior_faces.contains(&face_key(*face)))
            {
                return Err(invalid(
                    "all velocity modes on an exterior face must belong to that exterior",
                ));
            }
        }
        for boundary in boundaries {
            poll(cancelled)?;
            if boundary
                .faces
                .iter()
                .any(|&face| exterior_faces.contains(&face_key(face)))
            {
                return Err(invalid(
                    "a local resistance and exterior reaction cannot own the same face",
                ));
            }
        }
        let prepared = model.prepare_exact_time_with_pressure_release(
            duration_s,
            boundaries,
            ports,
            volume_sources,
            pressure_release_faces,
            limits.interior,
            cancelled,
        )?;
        let state = prepared.zero_state();
        let zeros = vec![0.0; prepared.volume_source_count()];
        let mut amplitudes = vec![0.0; prepared.port_count()];
        let mut interior_direct = vec![0.0; n * n];
        // Unit columns are cold contractions of the SAME cached exact operator.
        // No second propagator, field reduction or separate pressure quadrature.
        for (column, &port) in exterior_ports.iter().enumerate() {
            poll(cancelled)?;
            amplitudes[port] = 1.0;
            let response =
                prepared.step_with_volume_sources(&state, &amplitudes, &zeros, cancelled)?;
            amplitudes[port] = 0.0;
            for (row, &observation) in exterior_ports.iter().enumerate() {
                interior_direct[row * n + column] =
                    response.integrated_port_pressure_pa_m2_s[observation] / duration_s;
            }
        }
        let combined: Vec<_> = interior_direct
            .iter()
            .zip(exterior.mean_direct_term())
            .map(|(&interior, &outside)| interior + outside)
            .collect();
        if combined
            .iter()
            .chain(&interior_direct)
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "original same-interval port response is not representable",
            ));
        }
        // This <=128-port solve is bounded independently of full FE/exterior
        // state counts. The native factor is reused for every accepted interval.
        poll(cancelled)?;
        let closure = lu(&combined, n).map_err(CausalTimeError::Closure)?;
        poll(cancelled)?;
        Ok(Self {
            prepared,
            state,
            exterior,
            exterior_ports: exterior_ports.to_vec(),
            driven_ports,
            interior_direct,
            closure,
            limits,
            accepted_steps: 0,
        })
    }

    #[must_use]
    pub fn state(&self) -> &HeterogeneousTimeState<'m, 'f, 'c> {
        &self.state
    }
    #[must_use]
    pub fn interior(&self) -> &PreparedExactTime<'m, 'f, 'c> {
        &self.prepared
    }
    /// Borrow the unchanged source model, reports and accepted histories.
    #[must_use]
    pub fn exterior(&self) -> &CausalExteriorLoad {
        &self.exterior
    }
    #[must_use]
    pub fn exterior_ports(&self) -> &[usize] {
        &self.exterior_ports
    }
    #[must_use]
    pub fn driven_ports(&self) -> &[usize] {
        &self.driven_ports
    }
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Set the complete original fluid field before the first accepted interval.
    /// Source histories are retained. The caller supplies the initial physical
    /// pair deliberately; no modal projection or exterior-history reset occurs.
    pub fn initialize_field(
        &mut self,
        potential_pa_s: &[f64],
        pressure_pa: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CausalTimeError> {
        poll(cancelled)?;
        if self.accepted_steps != 0 {
            return Err(invalid(
                "accepted coupled histories cannot be reinitialized independently",
            ));
        }
        let candidate = self.prepared.state(potential_pa_s, pressure_pa)?;
        self.prepared.energy(&candidate)?;
        poll(cancelled)?;
        self.state = candidate;
        Ok(())
    }

    /// Advance both original histories with held feed velocities and held volume
    /// injections. Every prescribed source requires an explicit finite amplitude,
    /// including zero. The original force matching uses means over this SAME
    /// interval, not previous/endpoint pressures. The two owners are exact for
    /// their held inputs; their finite-clock exchange is a discretization of the
    /// continuously coupled problem, and is not claimed exact for that problem.
    ///
    /// All fallible fluid, source, residual, work and cancellation checks precede
    /// acceptance. The native exterior step performs its own transactional finite
    /// state check last; only its success commits the fluid field and step count.
    /// No exterior history clone or per-step factor/exponential is needed.
    pub fn step(
        &mut self,
        inward_feed_velocity_m_s: &[f64],
        volume_flow_m3_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CausalTimeStep, CausalTimeError> {
        poll(cancelled)?;
        if inward_feed_velocity_m_s.len() != self.driven_ports.len()
            || volume_flow_m3_s.len() != self.prepared.volume_source_count()
            || inward_feed_velocity_m_s
                .iter()
                .chain(volume_flow_m3_s)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite held amplitude for every original prescribed source is required",
            ));
        }
        let next_step = self
            .accepted_steps
            .checked_add(1)
            .ok_or_else(|| invalid("coupled acoustic step count overflow"))?;
        let n = self.exterior_ports.len();
        let dt = self.prepared.duration_s();
        let mut amplitudes = vec![0.0; self.prepared.port_count()];
        for (&port, &velocity) in self.driven_ports.iter().zip(inward_feed_velocity_m_s) {
            amplitudes[port] = velocity;
        }
        let baseline = self.prepared.step_with_volume_sources(
            &self.state,
            &amplitudes,
            volume_flow_m3_s,
            cancelled,
        )?;
        let base_forces: Vec<_> = self
            .exterior_ports
            .iter()
            .map(|&port| baseline.integrated_port_pressure_pa_m2_s[port] / dt)
            .collect();
        let mut history = vec![0.0; n];
        self.exterior.mean_history_into(&mut history)?;
        let mut velocities: Vec<_> = base_forces
            .iter()
            .zip(&history)
            .map(|(&p, &h)| p - h)
            .collect();
        poll(cancelled)?;
        self.closure.solve(&mut velocities);
        if velocities.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "same-interval exterior velocity is not representable",
            ));
        }
        for (&port, &velocity) in self.exterior_ports.iter().zip(&velocities) {
            amplitudes[port] = -velocity;
        }
        let candidate = self.prepared.step_with_volume_sources(
            &self.state,
            &amplitudes,
            volume_flow_m3_s,
            cancelled,
        )?;
        let interior_forces: Vec<_> = self
            .exterior_ports
            .iter()
            .map(|&port| candidate.integrated_port_pressure_pa_m2_s[port] / dt)
            .collect();
        // Match the original load's arithmetic exactly. The final native step
        // repeats this bounded observation and checks its candidate pole/history
        // states, before changing any accepted state. Only that owner advances it.
        let mut exterior_forces = history.clone();
        let direct = self.exterior.mean_direct_term();
        for (row, effort) in exterior_forces.iter_mut().enumerate() {
            *effort += direct[row * n..(row + 1) * n]
                .iter()
                .zip(&velocities)
                .map(|(d, u)| d * u)
                .sum::<f64>();
        }
        let mut scale = 0.0f64;
        let mut error = 0.0f64;
        for row in 0..n {
            poll(cancelled)?;
            scale = scale
                .max(interior_forces[row].abs())
                .max(exterior_forces[row].abs())
                .max(base_forces[row].abs())
                .max(history[row].abs());
            error = error.max((interior_forces[row] - exterior_forces[row]).abs());
            for column in 0..n {
                scale = scale
                    .max((self.interior_direct[row * n + column] * velocities[column]).abs())
                    .max((direct[row * n + column] * velocities[column]).abs());
            }
        }
        if !scale.is_finite()
            || !error.is_finite()
            || interior_forces
                .iter()
                .chain(&exterior_forces)
                .any(|v| !v.is_finite())
        {
            return Err(invalid("original interval-mean interface force overflow"));
        }
        let relative = if scale == 0.0 { 0.0 } else { error / scale };
        if relative > self.limits.interface_relative_tolerance {
            return Err(CausalTimeError::InterfaceResidual {
                relative,
                tolerance: self.limits.interface_relative_tolerance,
            });
        }
        let mut feed_work = Dd::from_f64(0.0);
        for (&port, &velocity) in self.driven_ports.iter().zip(inward_feed_velocity_m_s) {
            feed_work = feed_work
                + Dd::from_f64(velocity)
                    * Dd::from_f64(candidate.integrated_port_pressure_pa_m2_s[port]);
        }
        let mut exterior_work = Dd::from_f64(0.0);
        let mut fluid_outgoing_work = Dd::from_f64(0.0);
        for (row, &port) in self.exterior_ports.iter().enumerate() {
            exterior_work = exterior_work
                + Dd::from_f64(dt)
                    * Dd::from_f64(velocities[row])
                    * Dd::from_f64(exterior_forces[row]);
            fluid_outgoing_work = fluid_outgoing_work
                + Dd::from_f64(velocities[row])
                    * Dd::from_f64(candidate.integrated_port_pressure_pa_m2_s[port]);
        }
        let work = CausalTimeWork {
            initial_fluid_energy_j: candidate.work.initial_energy_j,
            final_fluid_energy_j: candidate.work.final_energy_j,
            delta_fluid_energy_j: candidate.work.delta_energy_j,
            boundary_dissipation_j: candidate.work.boundary_dissipation_j,
            source_work_j: (feed_work + Dd::from_f64(candidate.work.volume_source_work_j)).to_f64(),
            boundary_drive_work_j: feed_work.to_f64(),
            volume_source_work_j: candidate.work.volume_source_work_j,
            exterior_exchange_work_j: exterior_work.to_f64(),
            interface_work_defect_j: (exterior_work - fluid_outgoing_work).to_f64(),
        };
        if [
            work.source_work_j,
            work.boundary_drive_work_j,
            work.volume_source_work_j,
            work.exterior_exchange_work_j,
            work.interface_work_defect_j,
            work.balance_residual_j(),
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "original coupled acoustic work is not representable",
            ));
        }
        poll(cancelled)?;
        // Final atomic acceptance: the native load changes nothing on failure.
        // No fallible operation or cancellation checkpoint follows its success.
        self.exterior
            .step_mean_into(&velocities, &mut exterior_forces)?;
        self.state = candidate.state;
        self.accepted_steps = next_step;
        Ok(CausalTimeStep {
            outward_velocity_m_s: velocities,
            interior_mean_force_n: interior_forces,
            exterior_mean_force_n: exterior_forces,
            interface_relative_residual: relative,
            work,
        })
    }
}

#[cfg(test)]
mod tests;
