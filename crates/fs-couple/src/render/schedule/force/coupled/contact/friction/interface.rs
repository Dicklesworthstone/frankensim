//! Source-bound dry-interface applicability on the actual reciprocal motion.
//!
//! fs-tribo owns the law, ordered surface/history identity and its state domain.
//! This binding adds a contact temperature and explicit nominal area.
//! It does not turn bulk elasticity into friction, infer a contact patch, or
//! upgrade caller-declared/synthetic data into measured material evidence.

use super::*;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxisKind, SourceThermalCoulomb, SourceThermalCoulombError,
};
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_qty::Temperature;
use fs_tribo::{DryFrictionApplicability, DryFrictionState, DryInterfaceSystemCard};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
struct TemperatureSource {
    law: Arc<SourceThermalCoulomb>,
    resolved: ResolvedInterfaceStatePoint,
}

/// Immutable interface data and the geometry/state needed to query its domain.
/// The nominal area is NOT a Hertz patch and does not determine normal load.
/// Ordinary stepping holds temperature and area. An explicit temperature source
/// instead requires the shared finite-cell heat transaction, which refreshes
/// its original coefficient claims at the jointly solved interface temperature.
/// The producer remains responsible for mapping the ordered surfaces
/// to the actual participants of the mechanical port.
#[derive(Clone, Debug, PartialEq)]
pub struct DryFrictionBinding {
    card: DryInterfaceSystemCard,
    temperature_kelvin: f64,
    nominal_area_m2: f64,
    temperature_source: Option<TemperatureSource>,
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
        Ok(Self { card, temperature_kelvin, nominal_area_m2, temperature_source: None })
    }

    /// Bind existing source-pinned Coulomb coefficients to finite contact heat.
    /// Only temperature may change coefficient VALUES here. Explicit pressure
    /// and sliding-speed applicability are still checked on actual mechanics;
    /// pressure-, speed-, path- or surface-temperature-dependent responses need
    /// their own joint constitutive unknowns and refuse this owner.
    pub fn new_temperature_source(
        source: SourceThermalCoulomb,
        applicability: DryFrictionApplicability,
        nominal_area_m2: f64,
    ) -> Result<Self, ModalFrictionError> {
        if source.query_axes().iter().any(|axis| axis.kind != FrictionQueryAxisKind::SlidingSpeed)
            || source.initial_source().properties().iter().any(|property|
                property.source_axes().iter().any(|(axis, dims)|
                    axis != &source.axes().temperature || *dims != Temperature::DIMS))
        {
            return Err(invalid("modal thermal Coulomb coefficients may respond only to their owned contact temperature").into());
        }
        let point = source.initial_source().query();
        let temperature = *point.axes().get(&source.axes().temperature)
            .ok_or_else(|| invalid("thermal friction lost its original temperature coordinate"))?;
        let (resolved, law) = source.coefficients_at_query(&point)?;
        let card = DryInterfaceSystemCard::new(source.interface().interface().clone(), law, applicability)?;
        let mut binding = Self::new(card, temperature, nominal_area_m2)?;
        binding.temperature_source = Some(TemperatureSource { law: Arc::new(source), resolved });
        Ok(binding)
    }

    /// Original coefficient owner, never a fitted replacement for its claims.
    pub fn temperature_source(&self) -> Option<&SourceThermalCoulomb> {
        self.temperature_source.as_ref().map(|source| source.law.as_ref())
    }

    /// Latest accepted coefficient receipts at actual average slip and pressure.
    pub fn resolved_temperature_source(&self) -> Option<&ResolvedInterfaceStatePoint> {
        self.temperature_source.as_ref().map(|source| &source.resolved)
    }
    /// Complete ordered surface/history, law, authority and applicability.
    #[must_use] pub const fn card(&self) -> &DryInterfaceSystemCard { &self.card }
    /// Explicit fixed interface temperature, independent of ambient gas [K].
    #[must_use] pub const fn temperature_kelvin(&self) -> f64 { self.temperature_kelvin }
    /// Original nominal contact area used to map prescribed load to pressure.
    #[must_use] pub const fn nominal_area_m2(&self) -> f64 { self.nominal_area_m2 }

    /// Re-admit the same interface and area at a thermal trial. An explicit
    /// temperature owner refreshes its pinned coefficients; a fixed law stays
    /// fixed. Actual pressure/rates are admitted separately by the staged solve.
    pub(crate) fn at_temperature(&self, temperature_kelvin: f64) -> Result<Self, ModalFrictionError> {
        let Some(source) = &self.temperature_source else {
            return Self::new(self.card.clone(), temperature_kelvin, self.nominal_area_m2);
        };
        let point = source.resolved.query().with(&source.law.axes().temperature, temperature_kelvin)
            .map_err(SourceThermalCoulombError::from)?;
        let (resolved, law) = source.law.coefficients_at_query(&point)?;
        let card = DryInterfaceSystemCard::new(
            self.card.interface().clone(), law, self.card.applicability(),
        )?;
        let mut binding = Self::new(card, temperature_kelvin, self.nominal_area_m2)?;
        binding.temperature_source = Some(TemperatureSource { law: source.law.clone(), resolved });
        Ok(binding)
    }

    fn source_at_motion(&self, pressure: f64, speed: f64)
        -> Result<Option<ResolvedInterfaceStatePoint>, ModalFrictionError>
    {
        let Some(source) = &self.temperature_source else { return Ok(None); };
        let mut point = source.resolved.query();
        if let Some(axis) = &source.law.axes().nominal_pressure {
            point = point.with(axis, pressure).map_err(SourceThermalCoulombError::from)?;
        }
        for axis in source.law.query_axes() {
            point = point.with(&axis.name, speed.abs()).map_err(SourceThermalCoulombError::from)?;
        }
        let (resolved, law) = source.law.coefficients_at_query(&point)?;
        if law != *self.card.friction_law() {
            return Err(invalid("actual pressure/speed changed a temperature-only friction coefficient").into());
        }
        Ok(Some(resolved))
    }

    pub(crate) fn at_average_motion(&self, normal: f64, speed: f64)
        -> Result<Self, ModalFrictionError>
    {
        let mut binding = self.clone();
        if normal != 0.0 {
            self.check_speed(normal, speed)?;
            if let Some(resolved) = self.source_at_motion(normal / self.nominal_area_m2, speed)? {
                binding.temperature_source.as_mut().expect("retained source").resolved = resolved;
            }
        }
        Ok(binding)
    }

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
        self.source_at_motion(pressure, relative_speed)?;
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
