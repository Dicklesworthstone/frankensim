//! Source-bound dry-interface applicability on the actual reciprocal motion.
//!
//! fs-tribo owns the law, ordered surface/history identity and its state domain.
//! This binding adds a fixed contact temperature and explicit nominal area.
//! It does not turn bulk elasticity into friction, infer a contact patch, or
//! upgrade caller-declared/synthetic data into measured material evidence.

use super::*;
use fs_tribo::{DryFrictionState, DryInterfaceSystemCard};

/// Immutable interface data and the geometry/state needed to query its domain.
/// The nominal area is NOT a Hertz patch and does not determine normal load.
/// Temperature and area are frozen during a realization; no heat/wear model is
/// implied. The producer remains responsible for mapping the ordered surfaces
/// to the actual participants of the mechanical port.
#[derive(Clone, Debug, PartialEq)]
pub struct DryFrictionBinding {
    card: DryInterfaceSystemCard,
    temperature_kelvin: f64,
    nominal_area_m2: f64,
}
impl DryFrictionBinding {
    /// Bind the complete original interface, without copying selected mu values.
    /// Source-temperature applicability is checked even before contact exists.
    /// A positive normal load is checked at runtime using pressure N/area.
    pub fn new(card: DryInterfaceSystemCard, temperature_kelvin: f64, nominal_area_m2: f64)
        -> Result<Self, ModalFrictionError>
    {
        if !nominal_area_m2.is_finite() || nominal_area_m2 <= 0.0 {
            return Err(invalid("dry-interface binding requires a finite positive nominal contact area").into());
        }
        // Exercise the ORIGINAL owner's allocation-free state admission at
        // in-domain pressure/speed, without inventing a physical sample.
        let domain = card.applicability();
        let speed = domain.slip_speed_mps().minimum();
        card.validate_state(DryFrictionState::new(temperature_kelvin,
            domain.contact_pressure_pa().minimum(), speed)?)?;
        Ok(Self { card, temperature_kelvin, nominal_area_m2 })
    }
    /// Complete ordered surface/history, law, authority and applicability.
    #[must_use] pub const fn card(&self) -> &DryInterfaceSystemCard { &self.card }
    /// Explicit fixed interface temperature, independent of ambient gas [K].
    #[must_use] pub const fn temperature_kelvin(&self) -> f64 { self.temperature_kelvin }
    /// Original nominal contact area used to map prescribed load to pressure.
    #[must_use] pub const fn nominal_area_m2(&self) -> f64 { self.nominal_area_m2 }

    fn check_speed(&self, normal: f64, relative_speed: f64) -> Result<(), ModalFrictionError> {
        // No interface traction exists while released. A card may legitimately
        // require strictly positive pressure; release must not query it at zero.
        if normal == 0.0 { return Ok(()); }
        let pressure = normal / self.nominal_area_m2;
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(invalid("positive-load nominal contact pressure must be representable").into());
        }
        let state = DryFrictionState::new(self.temperature_kelvin, pressure, relative_speed.abs())?;
        // No per-sample source/receipt clone and no second traction evaluation.
        // The SAME admission is also used by the original receipt-producing query.
        self.card.validate_state(state)?;
        Ok(())
    }
}
/// Validate the true contact speed in the unchanged component/basis order.
/// This includes torsion or a moving opposite body when present in `column`.
pub(super) fn check_endpoint(binding: Option<&DryFrictionBinding>, column: &[f64],
    models: &[ModalAcousticTimeModel], drive: FrictionDrive) -> Result<(), ModalFrictionError>
{
    let Some(binding) = binding else { return Ok(()); };
    if drive.normal_force_n == 0.0 { return Ok(()); }
    let mut velocity = 0.0;
    let mut index = 0;
    for model in models {
        for state in model.states() {
            velocity = finite(velocity + column[index] * state.velocity_m_sqrt_kg_per_s)?;
            index += 1;
        }
    }
    binding.check_speed(drive.normal_force_n, finite(drive.speed_m_s - velocity)?)
}

/// The actual work-conjugate rate is checked separately from endpoint speeds.
/// Trial scalar root probes are algebraic values, not accepted physical states.
/// This does not certify all continuous-time extrema between sampled endpoints.
pub(super) fn check_average(binding: Option<&DryFrictionBinding>, drive: FrictionDrive,
    slip: f64, dt: f64) -> Result<(), ModalFrictionError>
{
    if let Some(binding) = binding {
        if drive.normal_force_n != 0.0 {
            binding.check_speed(drive.normal_force_n, finite(slip / dt)?)?;
        }
    }
    Ok(())
}

impl FrictionModalSystem {
    /// Construct the existing solver with its force law taken from this source.
    /// Geometry, original mechanical state and pressure observations are retained.
    /// No fallback coefficient, regularized traction or new integrator is used.
    pub fn new_with_interface(network: CoupledModalSystem, left: ModalAttachment,
        right: Option<ModalAttachment>, binding: DryFrictionBinding,
        config: ModalFrictionConfig, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        let port = ModalFrictionPort { left, right, law: binding.card.friction_law().clone() };
        Self::new(network, port, config, gate)?.with_interface(binding)
    }

    /// Attach source applicability to an IDENTICAL retained constitutive law.
    /// This may bind an already vibrating system: it changes no mechanics or
    /// history, and it makes no retroactive claim about previously accepted steps.
    /// A different law or an already bound source refuses rather than silently
    /// replacing an active interface. The next positive-load step checks its
    /// initial, actual step-average and endpoint rates before state publication.
    pub fn with_interface(mut self, binding: DryFrictionBinding) -> Result<Self, ModalFrictionError> {
        if self.interface.is_some() || self.adhesive.is_some() || self.port.law != *binding.card.friction_law() {
            return Err(invalid("dry-interface binding must match the unchanged law and cannot replace a retained source").into());
        }
        self.interface = Some(binding);
        Ok(self)
    }
    /// Original source and applicability attached to subsequent physical steps.
    #[must_use] pub const fn interface(&self) -> Option<&DryFrictionBinding> { self.interface.as_ref() }
}
