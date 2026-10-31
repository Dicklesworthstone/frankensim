//! Fixed Hertz-preload material microslip in the original simultaneous solve.
//! All constitutive history, traction quadrature and loss belong to the original
//! source-bound Cattaneo--Mindlin owner. This adapter supplies its actual modal
//! compliance and publishes its candidate only with the complete mechanical step.

use super::*;
use fs_contact::interface_binding::mindlin::{SourceMindlinContact, SourceMindlinState};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput};

/// Source-derived elastic contact storage and physical microslip/gross-slip heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MindlinFrame {
    /// Native constitutive branch. This owner's Sliding includes dissipative
    /// partial slip; it must not be interpreted as proof of gross sliding.
    pub regime: fs_tribo::FrictionRegime,
    /// Accepted endpoint force, distinct from the work-average port traction [N].
    pub endpoint_force_n: f64,
    /// Accepted signed macro tangential displacement [m].
    pub displacement_m: f64,
    /// Exact retained history intervals, bounded by the original source law.
    pub history_intervals: usize,
    /// Reversible tangential contact energy [J]; fixed normal preload is separate.
    pub stored_energy_j: f64,
    /// Change of reversible tangential energy in this accepted interval [J].
    pub stored_change_j: f64,
    /// Physical microslip and gross-slip heat, never total traction work [J].
    pub dissipated_j: f64,
    /// Original constitutive storage + heat - work residual [J].
    pub closure_residual_j: f64,
}

#[derive(Clone, Debug)]
pub(super) struct MaterialMindlinContact {
    source: SourceMindlinContact,
    state: SourceMindlinState,
    staged: Option<TangentialStep<SourceMindlinState>>,
}
impl MaterialMindlinContact {
    fn new(source: SourceMindlinContact) -> Self {
        let state = source.initial_state();
        Self { source, state, staged: None }
    }
    pub(super) fn normal_force_n(&self) -> f64 { self.source.normal_response().normal_force_n }
    pub(super) fn energy_j(&self) -> f64 { self.source.stored_energy_j(&self.state) }
    pub(super) fn staged_energy_j(&self) -> f64 {
        self.staged.as_ref().map_or(self.energy_j(), |step| self.source.stored_energy_j(&step.state))
    }
    fn trial(&self, free: f64, compliance: f64, normal: f64, dt: f64,
        config: ModalFrictionConfig) -> Result<TangentialStep<SourceMindlinState>, ModalFrictionError>
    {
        Ok(self.source.advance(&self.state, &TangentialStepInput {
            normal_n: normal, slip_increment_m: [free, 0.0], dt_s: dt,
            contact_temperature_k: Some(self.source.temperature_kelvin()),
            nominal_pressure_pa: None, compliance_m_per_n: compliance,
            solve: TangentialSolveConfig { max_iterations: config.max_iterations,
                maximum_force_n: config.maximum_force_n,
                force_absolute_tolerance_n: config.force_absolute_tolerance_n,
                force_relative_tolerance: config.force_relative_tolerance },
        })?)
    }
    pub(super) fn solve(&self, free: f64, compliance: f64, normal: f64, dt: f64,
        config: ModalFrictionConfig)
        -> Result<(f64, ModalFrictionRegime, usize, Option<FrictionRoots>), ModalFrictionError>
    {
        let step = self.trial(free, compliance, normal, dt, config)?;
        Ok((step.traction_n[0], ModalFrictionRegime::Microslip,
            step.iterations, None))
    }
    fn actual(&self, slip: f64, traction: f64, normal: f64, dt: f64,
        config: ModalFrictionConfig, iterations: usize)
        -> Result<(TangentialStep<SourceMindlinState>, f64), ModalFrictionError>
    {
        limit("material Mindlin slip speed", finite(slip / dt)?.abs(), config.maximum_slip_speed_m_s)?;
        limit("material Mindlin traction", traction.abs(), config.maximum_force_n)?;
        // Re-evaluate the original history on the ACTUAL jointly solved slip.
        // Zero compliance means body motion is already known, not a rigid body.
        let step = self.trial(slip, 0.0, normal, dt, config)?;
        let residual = finite(traction - step.traction_n[0])?;
        let tolerance = scalar::force_tolerance(config, traction, step.traction_n[0])?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::ContactSolve { residual_n: residual,
                tolerance_n: tolerance, iterations }.into());
        }
        Ok((step, residual))
    }
    pub(super) fn check(&self, slip: f64, traction: f64, normal: f64, dt: f64,
        config: ModalFrictionConfig, iterations: usize) -> Result<(), ModalFrictionError>
    {
        self.actual(slip, traction, normal, dt, config, iterations).map(|_| ())
    }
    pub(super) fn stage(&mut self, slip: f64, traction: f64, normal: f64, dt: f64,
        config: ModalFrictionConfig, iterations: usize) -> Result<(MindlinFrame, f64), ModalFrictionError>
    {
        let (step, residual) = self.actual(slip, traction, normal, dt, config, iterations)?;
        let frame = MindlinFrame { regime: step.regime, endpoint_force_n: step.state.endpoint_force_n(),
            displacement_m: step.state.slip_m(), history_intervals: step.state.history_intervals(),
            stored_energy_j: self.source.stored_energy_j(&step.state),
            stored_change_j: step.stored_change_j, dissipated_j: step.dissipated_j,
            closure_residual_j: step.closure_residual_j };
        self.staged = Some(step);
        Ok((frame, residual))
    }
    pub(super) fn publish(&mut self) {
        if let Some(step) = self.staged.take() { self.state = step.state; }
    }
}

impl MultiFrictionModalSystem {
    /// Bind original material-derived fixed-load Hertz--Mindlin sources once,
    /// before playback. Other ports retain their existing mechanisms. The full
    /// modal response, series-drive compliance and nonlinear strain iteration
    /// include the retained contact storage and physical microslip heat.
    ///
    /// Each bound drive MUST prescribe zero additional normal force. Its normal
    /// force comes from the source's original Hertz indentation. Changing load,
    /// opening, migration or temperature would require an explicit work/history
    /// transfer and is not approximated by reinitializing the contact. The fixed
    /// normal preload was supplied before this window and is never charged as
    /// tangential work or heat. This is an isothermal, proportional-tangent rung.
    pub fn with_material_mindlin_contacts(mut self,
        bindings: Vec<(usize, SourceMindlinContact)>, gate: &CancelGate)
        -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if bindings.is_empty() { return Ok(self); }
        if self.samples_rendered() != 0 || self.mindlin_contacts.iter().any(Option::is_some)
            || self.friction_residual.is_some() {
            return Err(invalid("Mindlin contacts require one complete binding before playback, without residual-tail compliance").into());
        }
        for (port, source) in bindings {
            poll(Some(gate))?;
            if port >= self.port_count() || self.mindlin_contacts[port].is_some()
                || self.shear_layers[port].is_some() || self.interface(port).is_some()
                || self.contact_loads[port].is_some() {
                return Err(invalid("one fixed-normal Mindlin source per actual port, without a competing interface or normal-load owner").into());
            }
            let contact = MaterialMindlinContact::new(source);
            // Admit the original source's clock/envelope before retaining it.
            contact.trial(0.0, 0.0, contact.normal_force_n(), self.sample_period_s(), self.points[port].config)?;
            self.mindlin_contacts[port] = Some(contact);
        }
        poll(Some(gate))?;
        Ok(self)
    }
    /// Original bulk, normal-model and ordered-interface receipts at this port.
    #[must_use]
    pub fn material_mindlin_contact(&self, port: usize) -> Option<&SourceMindlinContact> {
        self.mindlin_contacts.get(port).and_then(Option::as_ref).map(|contact| &contact.source)
    }
    /// Accepted source-owned history; no rejected candidate is exposed.
    #[must_use]
    pub fn material_mindlin_state(&self, port: usize) -> Option<&SourceMindlinState> {
        self.mindlin_contacts.get(port).and_then(Option::as_ref).map(|contact| &contact.state)
    }
}
