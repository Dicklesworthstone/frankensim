//! Reciprocal coupling of complete resolved thermoviscous fields to a DG bore.
//!
//! The original FE holds pressure traction over one physical interval; the
//! original DG holds the equal and opposite mean volume flow. Their cached
//! exact mean relations determine one common pressure before either candidate
//! is propagated. Every fluid, molecular, finite-solid and structural state is
//! retained. Coupling-clock refinement remains necessary independently of the
//! original spatial discretizations. This is not a joint continuous exponential
//! or a nonlinear mean-flow/jet model.

use super::ThermoviscousInterior;
use crate::acoustic_realize::AcousticRealizeError;
use crate::dg_port::{DgPortLoad, sum_terminal_work};
use fs_duct::dg::DgStepOutput;
use fs_exec::CancelGate;
use fs_feec::acoustics::thermoviscous::time::{
    PreparedThermoviscousTime, ThermoviscousTimeError, ThermoviscousTimeState,
    ThermoviscousTimeStep,
};

/// Physical-pressure playback through the shared observation owners.
pub mod performance;

/// Refusal by an original participant or the common pressure/work admission.
#[derive(Debug)]
pub enum ThermoviscousDgError {
    Time(ThermoviscousTimeError),
    Duct(AcousticRealizeError),
    Invalid(&'static str),
    Cancelled,
}
impl From<ThermoviscousTimeError> for ThermoviscousDgError {
    fn from(error: ThermoviscousTimeError) -> Self { Self::Time(error) }
}
impl From<AcousticRealizeError> for ThermoviscousDgError {
    fn from(error: AcousticRealizeError) -> Self { Self::Duct(error) }
}
impl core::fmt::Display for ThermoviscousDgError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Time(error) => error.fmt(f),
            Self::Duct(error) => error.fmt(f),
            Self::Invalid(what) => f.write_str(what),
            Self::Cancelled => f.write_str("resolved field/DG coupling cancelled"),
        }
    }
}
impl std::error::Error for ThermoviscousDgError {}
fn invalid(what: &'static str) -> ThermoviscousDgError {
    ThermoviscousDgError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), ThermoviscousDgError> {
    if gate.is_requested() { Err(ThermoviscousDgError::Cancelled) } else { Ok(()) }
}

/// Explicit physical and numerical bounds, additional to both owners' budgets.
#[derive(Clone, Copy, Debug)]
pub struct ThermoviscousDgLimits {
    pub maximum_steps: u64,
    /// Mean interface, authored traction and FE nodal pressure ceiling [Pa].
    pub maximum_abs_pressure_pa: f64,
    /// Complete FE and original DG/branch/gap storage ceiling [J].
    /// An attached rational terminal supplies work, not certified storage.
    pub maximum_total_energy_j: f64,
    pub energy_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub interface_relative_tolerance: f64,
}
impl Default for ThermoviscousDgLimits {
    fn default() -> Self {
        Self {
            maximum_steps: 10_000_000,
            maximum_abs_pressure_pa: 1_000_000.0,
            maximum_total_energy_j: 1.0,
            energy_relative_tolerance: 1e-8,
            energy_absolute_tolerance_j: 1e-24,
            interface_relative_tolerance: 1e-10,
        }
    }
}
impl ThermoviscousDgLimits {
    fn validate(self) -> Result<(), ThermoviscousDgError> {
        if self.maximum_steps == 0
            || !self.maximum_abs_pressure_pa.is_finite() || self.maximum_abs_pressure_pa <= 0.0
            || !self.maximum_total_energy_j.is_finite() || self.maximum_total_energy_j <= 0.0
            || !self.energy_relative_tolerance.is_finite() || self.energy_relative_tolerance <= 0.0
            || self.energy_relative_tolerance > 1e-6
            || !self.energy_absolute_tolerance_j.is_finite() || self.energy_absolute_tolerance_j < 0.0
            || !self.interface_relative_tolerance.is_finite() || self.interface_relative_tolerance <= 0.0
            || self.interface_relative_tolerance > 1e-6
        {
            return Err(invalid("finite positive coupled field/DG budgets and bounded tolerances required"));
        }
        Ok(())
    }
}

/// One fully accepted common interval. Interface work cancels internally.
#[derive(Clone, Debug)]
pub struct ThermoviscousDgFrame {
    pub step: u64,
    /// Original complete fields and separate integrated heat/work channels.
    /// Its pressure work includes the internal DG reaction.
    pub field: ThermoviscousTimeStep,
    pub duct: DgStepOutput,
    /// Shared mean inlet pressure [Pa]; this is not microphone pressure.
    pub mean_pressure_pa: f64,
    /// Mean flow INTO the DG, opposite the selected FE inward volume flow.
    pub mean_flow_m3_s: f64,
    pub stored_energy_j: f64,
    pub storage_change_j: f64,
    /// Authored external FE tractions and structural forces only.
    pub authored_external_work_j: f64,
    /// Independently integrated FE heat plus the original DG loss receipt,
    /// including any matched-boundary absorption. DG loss may be inferred by
    /// its owner; this does not turn an unchecked DG ledger into an independent
    /// dissipation certificate or label every DG energy sink as material heat.
    pub material_dissipation_j: f64,
    pub terminal_work_j: f64,
    pub interface_work_residual_j: f64,
    pub balance_residual_j: f64,
    pub balance_tolerance_j: f64,
}

/// One immutable source/preparation and two privately owned accepted histories.
pub struct ThermoviscousDg<'a, 'c> {
    source: &'a ThermoviscousInterior<'c>,
    prepared: PreparedThermoviscousTime<'a, 'c>,
    initial_state: ThermoviscousTimeState,
    load: DgPortLoad,
    port: usize,
    mobility_m3_pa_s: f64,
    limits: ThermoviscousDgLimits,
    accepted_steps: u64,
    frame: Option<ThermoviscousDgFrame>,
}
impl<'a, 'c> ThermoviscousDg<'a, 'c> {
    /// Connect an original unit-pressure FE port to an additional DG inlet.
    /// The prepared fields must borrow this exact source FE. Gas, clock and
    /// planar interface area must agree; no material, end correction or direct
    /// mouth radiation is substituted at the internal connection. Downstream
    /// branches, material walls, finite gaps and terminal radiation remain
    /// owned by the supplied original duct.
    pub fn new(
        source: &'a ThermoviscousInterior<'c>,
        prepared: PreparedThermoviscousTime<'a, 'c>,
        load: DgPortLoad,
        port: usize,
        limits: ThermoviscousDgLimits,
        gate: &CancelGate,
    ) -> Result<Self, ThermoviscousDgError> {
        poll(gate)?;
        limits.validate()?;
        let area = source.inner().plane_wave_port_area_m2(port)
            .ok_or_else(|| invalid("DG connection requires an exclusive planar unit-pressure FE port"))?;
        // Pointer identity alone does not prove material identity: a caller can
        // prepare the same FE with equilibrium Cp or invented reservoirs. Keep
        // the actual frozen calorics and every original source process here.
        let calorics_match = match source.relaxation() {
            Some(original) => {
                prepared.frozen_cp_j_kg_k().to_bits()
                    == (original.frozen_cv_j_kg_k() + source.gas().specific_gas_constant).to_bits()
                    && prepared.caloric_reservoirs().len() == original.processes().len()
                    && prepared.caloric_reservoirs().iter().zip(original.processes()).all(|(a, b)| {
                        a.heat_capacity_j_kg_k.to_bits() == b.heat_capacity_j_kg_k.to_bits()
                            && a.relaxation_time_s.to_bits() == b.relaxation_time_s.to_bits()
                    })
            }
            None => prepared.caloric_reservoirs().is_empty()
                && prepared.frozen_cp_j_kg_k().to_bits() == source.inner().medium().cp_j_kg_k.to_bits(),
        };
        if !core::ptr::eq(source.inner(), prepared.inner())
            || !calorics_match
            || prepared.duration_s() * source.band_hz()[1] > 0.5
            || prepared.conductive_wall_layouts().len() != source.wall_sources().len()
            || prepared.maximum_wall_frequency_hz().is_some_and(|frequency|
                frequency.to_bits() != source.band_hz()[1].to_bits())
            || source.gas() != load.duct().inlet_gas()
            || prepared.duration_s().to_bits() != load.propagator().dt().to_bits()
            || load.state().step_count() != 0 || load.has_driven_mouth()
            || (area - load.duct().inlet_area_m2()).abs() > 64.0 * f64::EPSILON * area
        {
            return Err(invalid("resolved field/DG connection requires its original source, equal gas/area/clock and an unadvanced ordinary inlet"));
        }
        let initial_state = prepared.zero_state();
        let pressure = vec![0.0; source.inner().port_mode_count()];
        let forces = vec![0.0; initial_state.wall_velocity_m_sqrt_kg_per_s.len()];
        let (free, mobility_m3_pa_s) = prepared.average_relation(
            &initial_state, &pressure, &forces, port, &mut || gate.is_requested(),
        )?;
        let relation = load.port_response()?;
        if free != 0.0 || !mobility_m3_pa_s.is_finite() || mobility_m3_pa_s <= 0.0
            || !relation.impedance_pa_s_m3.is_finite() || relation.impedance_pa_s_m3 <= 0.0
        {
            return Err(invalid("positive original mean-port mobility and impedance required"));
        }
        let result = Self { source, prepared, initial_state, load, port, mobility_m3_pa_s,
            limits, accepted_steps: 0, frame: None };
        result.admit_state(result.state(), result.load.stored_energy_j())?;
        poll(gate)?;
        Ok(result)
    }

    /// Retain an authored complete initial field, including original molecular
    /// and finite-wall temperatures. This never resets an accepted trajectory.
    pub fn with_initial_state(mut self, state: ThermoviscousTimeState)
        -> Result<Self, ThermoviscousDgError>
    {
        if self.accepted_steps != 0 {
            return Err(invalid("initial resolved fields can only be supplied before the first coupled step"));
        }
        self.admit_state(&state, self.load.stored_energy_j())?;
        self.initial_state = state;
        Ok(self)
    }

    pub const fn source(&self) -> &'a ThermoviscousInterior<'c> { self.source }
    pub const fn prepared(&self) -> &PreparedThermoviscousTime<'a, 'c> { &self.prepared }
    pub const fn load(&self) -> &DgPortLoad { &self.load }
    pub fn state(&self) -> &ThermoviscousTimeState {
        self.frame.as_ref().map_or(&self.initial_state, |frame| &frame.field.state)
    }
    pub const fn port(&self) -> usize { self.port }
    pub const fn limits(&self) -> ThermoviscousDgLimits { self.limits }
    pub const fn accepted_steps(&self) -> u64 { self.accepted_steps }
    pub const fn last_frame(&self) -> Option<&ThermoviscousDgFrame> { self.frame.as_ref() }
    pub fn stored_energy_j(&self) -> Result<f64, ThermoviscousDgError> {
        Ok(self.prepared.stored_energy_j(self.state())? + self.load.stored_energy_j())
    }

    fn admit_state(&self, state: &ThermoviscousTimeState, duct_energy: f64)
        -> Result<f64, ThermoviscousDgError>
    {
        let energy = self.prepared.stored_energy_j(state)? + duct_energy;
        if !energy.is_finite() || energy > self.limits.maximum_total_energy_j
            || state.pressure_pa.iter().any(|p| p.abs() > self.limits.maximum_abs_pressure_pa)
        {
            return Err(invalid("coupled resolved pressure or complete stored energy exceeds its physical budget"));
        }
        Ok(energy)
    }

    /// Advance one common interval from complete held external drives. The
    /// connected pressure slot MUST be zero: its reaction is solved here.
    /// Neither participant nor the accepted frame changes after any refusal,
    /// including cancellation after either original candidate is propagated.
    pub fn step(
        &mut self,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<&ThermoviscousDgFrame, ThermoviscousDgError> {
        poll(gate)?;
        if self.accepted_steps >= self.limits.maximum_steps
            || traction_pa.len() != self.source.inner().port_mode_count()
            || traction_pa.get(self.port).copied() != Some(0.0)
            || traction_pa.iter().any(|p| !p.is_finite() || p.abs() > self.limits.maximum_abs_pressure_pa)
        {
            return Err(invalid("bounded coupled clock and complete external tractions with zero internal pressure slot required"));
        }
        let (free, mobility) = self.prepared.average_relation(
            self.state(), traction_pa, wall_forces_n_per_sqrt_kg,
            self.port, &mut || gate.is_requested(),
        )?;
        if mobility.to_bits() != self.mobility_m3_pa_s.to_bits() {
            return Err(invalid("original resolved port mobility changed"));
        }
        let relation = self.load.port_response()?;
        let denominator = 1.0 + mobility * relation.impedance_pa_s_m3;
        // Q_FE = free + mobility*p and Q_DG = -Q_FE; p = h + Z*Q_DG.
        let pressure = (relation.history_pressure_pa - relation.impedance_pa_s_m3 * free) / denominator;
        let flow = -(free + mobility * relation.history_pressure_pa) / denominator;
        if !denominator.is_finite() || denominator <= 0.0 || !pressure.is_finite()
            || !flow.is_finite() || pressure.abs() > self.limits.maximum_abs_pressure_pa
        {
            return Err(invalid("common resolved field/DG mean relation left its finite pressure domain"));
        }
        let mut applied = traction_pa.to_vec();
        applied[self.port] = pressure;
        let field = self.prepared.step(self.state(), &applied, wall_forces_n_per_sqrt_kg,
            &mut || gate.is_requested())?;
        let actual_flow = -field.mean_inward_flow_m3_s[self.port];
        let flow_scale = free.abs() + (mobility * pressure).abs() + actual_flow.abs();
        if (actual_flow - flow).abs() > self.limits.interface_relative_tolerance * flow_scale.max(1e-300) {
            return Err(invalid("resolved field and original port map disagree on swept volume"));
        }
        let dt = self.prepared.duration_s();
        let external_pressure_work: f64 = traction_pa.iter().zip(&field.mean_inward_flow_m3_s)
            .map(|(p, q)| p * (q * dt)).sum();
        let authored_external_work_j = external_pressure_work + field.work.structural_work_j;
        let inlet = self.load.inlet();
        let trial = self.load.preview_step(actual_flow)?;
        let duct = trial.output();
        let port = duct.port_outputs.iter().find(|p| p.id == inlet)
            .ok_or_else(|| invalid("original DG candidate omitted its connected inlet"))?;
        let pressure_scale = pressure.abs() + relation.history_pressure_pa.abs()
            + (relation.impedance_pa_s_m3 * actual_flow).abs();
        if !port.conjugate_rate.is_finite() || (port.conjugate_rate - pressure).abs()
            > self.limits.interface_relative_tolerance * pressure_scale.max(1e-300)
        {
            return Err(invalid("resolved field and DG candidates disagree on mean pressure"));
        }
        let stored_energy_j = field.work.final_energy_j + duct.ledger.stored_energy_j;
        let before = field.work.initial_energy_j + trial.previous_stored_energy_j();
        let storage_change_j = field.work.delta_energy_j
            + (duct.ledger.stored_energy_j - trial.previous_stored_energy_j());
        let material_dissipation_j = field.work.dissipated_j()
            + (duct.ledger.cumulative_dissipated_j - trial.previous_cumulative_dissipated_j());
        let terminal_work_j = sum_terminal_work(trial.terminal_frames());
        let interface_work_residual_j = trial.inlet_work_j()
            + pressure * (field.mean_inward_flow_m3_s[self.port] * dt);
        let balance_residual_j = storage_change_j + material_dissipation_j
            + terminal_work_j - authored_external_work_j;
        let scale = before.abs() + stored_energy_j.abs() + authored_external_work_j.abs()
            + material_dissipation_j.abs() + terminal_work_j.abs();
        let balance_tolerance_j = self.limits.energy_absolute_tolerance_j
            + self.limits.energy_relative_tolerance * scale;
        if [stored_energy_j, before, storage_change_j, material_dissipation_j,
            terminal_work_j, interface_work_residual_j, balance_residual_j, scale, balance_tolerance_j]
            .iter().any(|value| !value.is_finite())
            || stored_energy_j > self.limits.maximum_total_energy_j
            || field.state.pressure_pa.iter().any(|p| p.abs() > self.limits.maximum_abs_pressure_pa)
            || material_dissipation_j < -balance_tolerance_j
            || interface_work_residual_j.abs() > balance_tolerance_j
            || balance_residual_j.abs() > balance_tolerance_j
        {
            return Err(invalid("original resolved field/DG pressure, storage or same-window work admission refused"));
        }
        poll(gate)?;
        // Both publications are infallible. No cancellation/error separates them.
        let duct = trial.accept();
        self.accepted_steps += 1;
        self.frame = Some(ThermoviscousDgFrame { step: self.accepted_steps, field, duct,
            mean_pressure_pa: pressure, mean_flow_m3_s: actual_flow, stored_energy_j,
            storage_change_j, authored_external_work_j, material_dissipation_j, terminal_work_j,
            interface_work_residual_j, balance_residual_j, balance_tolerance_j });
        Ok(self.frame.as_ref().expect("joint field/DG step accepted"))
    }
}
