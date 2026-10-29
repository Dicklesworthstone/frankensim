//! Real adhesive junction shear in the original finite-body modal transaction.
use super::*;
use fs_contact::interface_binding::adhesive_shear::{SourceAdhesiveContactUpdate, SourceAdhesiveTangential};
use fs_tribo::adhesive::{AdhesiveContactState, AdhesiveTangentialContact};
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

/// Current real area/load plus the complete junction law. Optional source binding
/// retains both original material owners rather than copying fitted strengths.
/// Geometry is held over each sample; the external normal owner supplies updates.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveModalBinding {
    law: AdhesiveTangentialContact,
    source: Option<SourceAdhesiveTangential>,
}

/// Geometry supplied by the existing normal owner for one mechanical sample.
/// A sourced binding requires a response-derived update with the SAME original
/// normal identity. Authored geometry cannot overwrite a source-backed patch.
/// These inputs do not transfer normal work into the tangential energy account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AdhesiveModalContactUpdate {
    /// Actual JKR, DMT or Maugis solid contact geometry with its source identity.
    Source(SourceAdhesiveContactUpdate),
    /// Explicit geometry for an originally unsourced, caller-authored binding.
    Declared(AdhesiveContactState),
}
impl AdhesiveModalBinding {
    /// Bind explicitly declared physical contact geometry and coefficients.
    /// Coulomb and affine laws use the held-effort solve. Stribeck/hyperbolic
    /// pressure laws use the native network's retained branch continuation,
    /// including the independent static/kinetic junction-strength drop.
    pub fn new(law: AdhesiveTangentialContact) -> Result<Self, ModalFrictionError> {
        Ok(Self { law, source: None })
    }
    /// Retain the original normal/shear sources and their exact fixed query.
    pub fn from_source(source: SourceAdhesiveTangential) -> Result<Self, ModalFrictionError> {
        let mut binding = Self::new(source.tangential().clone())?;
        binding.source = Some(source);
        Ok(binding)
    }
    /// Actual complete junction/pressure law and last accepted normal patch.
    #[must_use] pub const fn tangential(&self) -> &AdhesiveTangentialContact { &self.law }
    /// Original source owners when explicitly supplied; never fabricated.
    #[must_use] pub const fn source(&self) -> Option<&SourceAdhesiveTangential> { self.source.as_ref() }

    fn prepare_contact(
        &self,
        update: AdhesiveModalContactUpdate,
    ) -> Result<AdhesiveTangentialContact, ModalFrictionError> {
        match (self.source.as_ref(), update) {
            (Some(source), AdhesiveModalContactUpdate::Source(update)) => {
                source.tangential_for_update(update).map_err(Into::into)
            }
            (None, AdhesiveModalContactUpdate::Declared(contact)) => {
                AdhesiveTangentialContact::new(self.law.law().clone(), contact).map_err(Into::into)
            }
            _ => Err(invalid(
                "adhesive geometry update must retain the original source or declared authority",
            )
            .into()),
        }
    }

    pub(super) fn solve(&self, config: ModalFrictionConfig, compliance: f64, dt: f64,
        free_slip: f64, normal: f64, gate: Option<&CancelGate>)
        -> Result<(f64, ModalFrictionRegime, usize), ModalFrictionError>
    {
        poll(gate)?;
        let x = TangentialStepInput { normal_n: normal, slip_increment_m: [free_slip,0.0],
            dt_s: dt, contact_temperature_k: self.source.as_ref().and_then(|s| s.temperature_kelvin()),
            nominal_pressure_pa: None, compliance_m_per_n: compliance,
            solve: TangentialSolveConfig { max_iterations: config.max_iterations,
                maximum_force_n: config.maximum_force_n,
                force_absolute_tolerance_n: config.force_absolute_tolerance_n,
                force_relative_tolerance: config.force_relative_tolerance } };
        // The staged leaf keeps the original source's complete law. Its fixed
        // source temperature above is unchanged; only actual area/load evolves.
        // Calling the source's accepted leaf here would use stale geometry.
        let step = self.law.advance(&self.law.initial_state(), &x)?;
        let regime = if self.law.contact().area_m2() == 0.0 { ModalFrictionRegime::Released }
            else if step.regime == FrictionRegime::Sticking { ModalFrictionRegime::Sticking }
            else { ModalFrictionRegime::Sliding };
        poll(gate)?;
        Ok((step.traction_n[0], regime, step.iterations))
    }
    pub(super) fn check(&self, config: ModalFrictionConfig, dt: f64,
        traction: f64, slip: f64, regime: ModalFrictionRegime, iterations: usize)
        -> Result<(f64,f64), ModalFrictionError>
    {
        match regime {
            ModalFrictionRegime::Released => {
                if self.law.contact().area_m2() != 0.0 || traction != 0.0 {
                    return Err(invalid("released adhesive patch must have zero real area and traction").into());
                }
                Ok((0.0,0.0))
            }
            ModalFrictionRegime::Sticking => {
                if slip.abs() > config.sticking_tolerance_m {
                    return Err(ModalFrictionError::Sticking { residual_m: slip, tolerance_m: config.sticking_tolerance_m });
                }
                let cap = self.law.law().static_capacity_n(self.law.contact())?;
                if traction.abs() > cap + scalar::force_tolerance(config, traction, cap)? {
                    return Err(invalid("actual static adhesive reaction exceeds its area/load capacity").into());
                }
                Ok((0.0,0.0))
            }
            ModalFrictionRegime::Sliding => {
                let speed = finite(slip/dt)?;
                limit("friction slip speed", speed.abs(), config.maximum_slip_speed_m_s)?;
                if speed == 0.0 || traction != 0.0 && traction.signum() != speed.signum() {
                    return Err(invalid("actual adhesive traction must oppose nonzero relative slip").into());
                }
                let expected = finite(speed.signum() * self.law.law().sliding_capacity_n(self.law.contact(), speed.abs())?)?;
                let residual = finite(traction - expected)?;
                let tolerance = scalar::force_tolerance(config, traction, expected)?;
                if residual.abs() > tolerance {
                    return Err(ModalCouplingError::ContactSolve { residual_n: residual, tolerance_n: tolerance, iterations }.into());
                }
                Ok((finite(traction*slip)?,residual))
            }
            ModalFrictionRegime::Shearing {..} => Err(invalid("adhesive junction friction has no elastic shear store").into()),
        }
    }
}
impl FrictionModalSystem {
    /// Construct the original network owner with a complete adhesive junction.
    pub fn new_with_adhesion(network: CoupledModalSystem, left: ModalAttachment,
        right: Option<ModalAttachment>, binding: AdhesiveModalBinding,
        config: ModalFrictionConfig, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        let port = ModalFrictionPort { left, right, law: binding.law.law().pressure_law().clone() };
        Self::new(network,port,config,gate)?.with_adhesion(binding)
    }
    /// Retain an explicit patch around the unchanged pressure law. Existing
    /// mechanical states, body compliance, observers and energy transaction stay
    /// intact. Replacing an existing source or patch refuses.
    pub fn with_adhesion(mut self, binding: AdhesiveModalBinding) -> Result<Self, ModalFrictionError> {
        if self.adhesive.is_some() || self.interface.is_some()
            || self.port.law != *binding.law.law().pressure_law() {
            return Err(invalid("adhesive binding must match the unchanged pressure law and cannot replace a retained source").into());
        }
        self.adhesive = Some(binding);
        // A former dry load line lacks the new junction force. Mechanical
        // motion is retained, but its old friction branch is not transferable.
        self.load_line = None;
        self.branch = ContactBranch::Fresh;
        self.staged_branch = ContactBranch::Fresh;
        Ok(self)
    }
    /// Complete retained adhesive geometry, law and original sources.
    #[must_use] pub const fn adhesive_binding(&self) -> Option<&AdhesiveModalBinding> { self.adhesive.as_ref() }

    /// Advance with normal-owner geometry held over this tangential sample.
    /// The SAME attached junction retains its occupied sliding branch as area
    /// or repulsive load changes. Actual zero-area detachment and reattachment
    /// start a fresh junction. Source query and temperature remain unchanged.
    /// Geometry, cache, branch and native mechanics publish together; refusal
    /// restores the accepted binding and history without copying source owners.
    /// Normal motion/work remain external, and no shear peeling is inferred.
    pub fn step_with_adhesive_contact(
        &mut self,
        external: &[f64],
        drive: FrictionDrive,
        update: AdhesiveModalContactUpdate,
    ) -> Result<&FrictionModalFrame, ModalFrictionError> {
        self.step_with_adhesive_contact_inner(external, drive, update, None)
    }

    /// The same geometry/mechanics transaction with bounded cancellation polls.
    pub fn step_with_adhesive_contact_under_gate(
        &mut self,
        external: &[f64],
        drive: FrictionDrive,
        update: AdhesiveModalContactUpdate,
        gate: &CancelGate,
    ) -> Result<&FrictionModalFrame, ModalFrictionError> {
        self.step_with_adhesive_contact_inner(external, drive, update, Some(gate))
    }

    fn step_with_adhesive_contact_inner(
        &mut self,
        external: &[f64],
        drive: FrictionDrive,
        update: AdhesiveModalContactUpdate,
        gate: Option<&CancelGate>,
    ) -> Result<&FrictionModalFrame, ModalFrictionError> {
        poll(gate)?;
        let binding = self.adhesive.as_mut().ok_or_else(|| {
            invalid("changing adhesive geometry requires an existing adhesive binding")
        })?;
        let candidate = binding.prepare_contact(update)?;
        let before_contact = binding.law.contact();
        let after_contact = candidate.contact();
        // Move only the small constitutive leaf. Original source allocations
        // remain in place throughout the trial, including late refusals.
        let previous_law = core::mem::replace(&mut binding.law, candidate);
        let previous_line = self.load_line.clone();
        let previous_branch = self.branch;
        let previous_staged_branch = self.staged_branch;
        if before_contact != after_contact {
            // Area changes the junction offset even at unchanged normal load.
            self.load_line = None;
        }
        if before_contact.area_m2() == 0.0 || after_contact.area_m2() == 0.0 {
            self.branch = ContactBranch::Fresh;
        }
        let trial = (|| {
            let candidate = self.stage_candidate(external, drive, gate)?;
            self.network.validate_staged_material_points(gate)?;
            poll(gate)?;
            // Last fallible operation: source admission is repeated against the
            // unchanged owner before any accepted mechanics can publish.
            if let AdhesiveModalContactUpdate::Source(update) = update {
                self.adhesive
                    .as_mut()
                    .and_then(|binding| binding.source.as_mut())
                    .ok_or_else(|| invalid("missing retained adhesive source"))?
                    .update_contact(update)?;
            }
            Ok::<_, ModalFrictionError>(candidate)
        })();
        self.branch = previous_branch;
        match trial {
            Ok((accepted, nonlinear_frame)) => {
                Ok(self.publish_candidate(accepted, nonlinear_frame))
            }
            Err(error) => {
                if let Some(binding) = &mut self.adhesive {
                    binding.law = previous_law;
                }
                self.load_line = previous_line;
                self.staged_branch = previous_staged_branch;
                Err(error)
            }
        }
    }
}
