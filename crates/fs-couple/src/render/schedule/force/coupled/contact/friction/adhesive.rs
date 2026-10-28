//! Real adhesive junction shear in the original finite-body modal transaction.
use super::*;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveTangential;
use fs_tribo::adhesive::AdhesiveTangentialContact;
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

/// Fixed real area/load plus the complete junction law. Optional source binding
/// retains both original material owners rather than copying fitted strengths.
/// The normal patch is explicit frozen data, not a new normal-contact solver.
#[derive(Clone, Debug, PartialEq)]
pub struct AdhesiveModalBinding {
    law: AdhesiveTangentialContact,
    source: Option<SourceAdhesiveTangential>,
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
    /// Actual complete junction/pressure law and frozen normal patch.
    #[must_use] pub const fn tangential(&self) -> &AdhesiveTangentialContact { &self.law }
    /// Original source owners when explicitly supplied; never fabricated.
    #[must_use] pub const fn source(&self) -> Option<&SourceAdhesiveTangential> { self.source.as_ref() }

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
        let step = if let Some(source) = &self.source { source.advance(&source.initial_state(), &x)? }
            else { self.law.advance(&self.law.initial_state(), &x)? };
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
    /// Retain an explicit fixed patch around the unchanged pressure law. Existing
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
}
