//! Sealed adiabatic chemistry at one prescribed pressure and total enthalpy.
//!
//! The original extent integrator, prepared Newton solve, step doubling and
//! physical admission gates remain the only evolution owner. Actual volume is
//! V=nRT/p at every trial; concentration dilution and volume-scaled rates enter
//! its analytic Jacobian. Formation energy occurs once in H. Boundary work is
//! p(Vfinal-Vinitial), positive outward, so U changes even though H is conserved.
//! This homogeneous quasistatic model does not resolve piston inertia, acoustics,
//! spatial transport, heat exchange, inflow or outflow.

use super::*;

pub use super::{
    IsochoricKineticsBalance as IsobaricKineticsBalance,
    IsochoricKineticsError as IsobaricKineticsError,
    IsochoricKineticsMethod as IsobaricKineticsMethod,
    IsochoricKineticsSolve as IsobaricKineticsSolve,
};

/// Prepared constant-H,p adapter over the original constrained reactor.
/// Successful advances allocate no storage and never mutate caller histories.
pub struct IsobaricReactor {
    inner: IsochoricReactor,
}

impl IsobaricReactor {
    pub fn new(
        mechanism: IdealGasKineticsV1,
        controls: IsobaricKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IsobaricKineticsError> {
        let mut inner = IsochoricReactor::new(mechanism, controls, &mut cancelled)?;
        inner.reactor.prepare_pressure()?;
        poll(&mut cancelled)?;
        Ok(Self { inner })
    }

    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.inner.mechanism()
    }
    pub const fn controls(&self) -> IsobaricKineticsSolve {
        self.inner.controls()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.inner.workspace_bytes()
    }

    /// Use the original feed enthalpy at the supplied physical temperature.
    pub fn advance_adiabatic(
        &mut self,
        amounts: &[Amount],
        pressure: Pressure,
        initial_temperature: Temperature,
        duration: Time,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricKineticsView<'_>, IsobaricKineticsError> {
        self.advance(
            amounts,
            None,
            pressure,
            initial_temperature,
            duration,
            cancelled,
        )
    }

    /// Resolve the unchanged feed to one original-source H before kinetic time.
    /// The supplied temperature is only a root hint; initial_volume() is the
    /// actual volume AFTER this initial enthalpy resolution. Absolute targets
    /// whose f64 spacing exceeds the Joule tolerance retain the shared refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adiabatic_enthalpy(
        &mut self,
        amounts: &[Amount],
        target_enthalpy: Energy,
        pressure: Pressure,
        initial_temperature: Temperature,
        duration: Time,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricKineticsView<'_>, IsobaricKineticsError> {
        self.advance(
            amounts,
            Some(target_enthalpy),
            pressure,
            initial_temperature,
            duration,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance(
        &mut self,
        amounts: &[Amount],
        target: Option<Energy>,
        pressure: Pressure,
        initial_temperature: Temperature,
        duration: Time,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricKineticsView<'_>, IsobaricKineticsError> {
        poll(&mut cancelled)?;
        if !positive(pressure.value()) || !positive(initial_temperature.value()) {
            return Err(IsobaricKineticsError::Invalid(
                "invalid isobaric pressure or temperature",
            ));
        }
        let mut count = Sum::default();
        for amount in amounts {
            if !amount.value().is_finite() || amount.value() < 0.0 {
                return Err(IsobaricKineticsError::Invalid("invalid dense feed amount"));
            }
            count.add(amount.value());
        }
        let reference_volume = Volume::new(
            count.value()
                * (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * initial_temperature.value()
                    / pressure.value()),
        );
        let inner = advance_prepared_constraint(
            &mut self.inner,
            amounts,
            target,
            reference_volume,
            Some(pressure),
            initial_temperature,
            duration,
            &mut cancelled,
        )?;
        let view = IsobaricKineticsView { inner, pressure };
        if !view.boundary_work().value().is_finite() {
            return Err(IsobaricKineticsError::Invalid(
                "unrepresentable isobaric boundary work",
            ));
        }
        poll(&mut cancelled)?;
        Ok(view)
    }
}

/// Fully accepted endpoint borrowing reusable storage. The balance energy
/// residuals refer to H; mass, elemental, charge and entropy gates are unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricKineticsView<'a> {
    inner: IsochoricKineticsView<'a>,
    pressure: Pressure,
}

impl IsobaricKineticsView<'_> {
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.inner.mechanism()
    }
    pub const fn amounts(&self) -> &[Amount] {
        self.inner.amounts()
    }
    pub const fn temperature(&self) -> Temperature {
        self.inner.temperature()
    }
    /// Prescribed physical pressure; pressure_residual_pa() retains EOS roundoff.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    pub fn pressure_residual_pa(&self) -> f64 {
        self.inner.pressure().value() - self.pressure.value()
    }
    pub const fn enthalpy(&self) -> Energy {
        self.inner.constrained_energy
    }
    pub const fn target_enthalpy(&self) -> Energy {
        self.inner.target_energy
    }
    pub const fn internal_energy(&self) -> Energy {
        self.inner.internal_energy()
    }
    pub const fn mass(&self) -> Mass {
        self.inner.mass()
    }
    pub const fn volume(&self) -> Volume {
        self.inner.volume()
    }
    pub const fn initial_volume(&self) -> Volume {
        self.inner.initial_volume
    }
    /// Outward mechanical work against the prescribed pressure, J.
    pub fn boundary_work(&self) -> Energy {
        Energy::new(self.pressure.value() * (self.volume().value() - self.initial_volume().value()))
    }
    pub const fn elapsed_time(&self) -> Time {
        self.inner.elapsed_time()
    }
    pub const fn balance(&self) -> IsobaricKineticsBalance {
        self.inner.balance()
    }
    pub fn frozen_state(
        &self,
        relative_charge_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FrozenIonizedGasStateV1, FrozenIonizedGasStateErrorV1> {
        self.inner
            .frozen_state(relative_charge_tolerance, cancelled)
    }
    /// Explicit owning conversion, including original endpoint source receipts.
    pub fn to_owned(
        self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricKineticsResult, IsobaricKineticsError> {
        Ok(IsobaricKineticsResult {
            inner: self.inner.to_owned(cancelled)?,
            pressure: self.pressure,
        })
    }
}

/// Independent accepted endpoint retaining the complete original mechanism.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricKineticsResult {
    inner: IsochoricKineticsResult,
    pressure: Pressure,
}

impl IsobaricKineticsResult {
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.inner.mechanism()
    }
    pub fn amounts(&self) -> &[Amount] {
        self.inner.amounts()
    }
    pub const fn temperature(&self) -> Temperature {
        self.inner.temperature()
    }
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    pub fn pressure_residual_pa(&self) -> f64 {
        self.inner.pressure().value() - self.pressure.value()
    }
    pub const fn enthalpy(&self) -> Energy {
        self.inner.constrained_energy
    }
    pub const fn target_enthalpy(&self) -> Energy {
        self.inner.target_energy
    }
    pub const fn internal_energy(&self) -> Energy {
        self.inner.internal_energy()
    }
    pub const fn mass(&self) -> Mass {
        self.inner.mass()
    }
    pub const fn volume(&self) -> Volume {
        self.inner.volume()
    }
    pub const fn initial_volume(&self) -> Volume {
        self.inner.initial_volume
    }
    pub fn boundary_work(&self) -> Energy {
        Energy::new(self.pressure.value() * (self.volume().value() - self.initial_volume().value()))
    }
    pub const fn elapsed_time(&self) -> Time {
        self.inner.elapsed_time()
    }
    pub const fn balance(&self) -> IsobaricKineticsBalance {
        self.inner.balance()
    }
    pub fn sources(&self) -> &[Nasa9EvaluationReceiptV1] {
        self.inner.sources()
    }
    pub fn frozen_state(
        &self,
        relative_charge_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FrozenIonizedGasStateV1, FrozenIonizedGasStateErrorV1> {
        self.inner
            .frozen_state(relative_charge_tolerance, cancelled)
    }
}

#[allow(clippy::too_many_arguments)]
pub fn advance_adiabatic(
    mechanism: &IdealGasKineticsV1,
    feed: &[(SpeciesId, Amount)],
    pressure: Pressure,
    initial_temperature: Temperature,
    duration: Time,
    controls: IsobaricKineticsSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<IsobaricKineticsResult, IsobaricKineticsError> {
    poll(&mut cancelled)?;
    let amounts = feed_amounts(mechanism, feed)?;
    let mut reactor = IsobaricReactor::new(mechanism.clone(), controls, &mut cancelled)?;
    reactor
        .advance_adiabatic(
            &amounts,
            pressure,
            initial_temperature,
            duration,
            &mut cancelled,
        )?
        .to_owned(cancelled)
}

#[allow(clippy::too_many_arguments)]
pub fn advance_adiabatic_enthalpy(
    mechanism: &IdealGasKineticsV1,
    feed: &[(SpeciesId, Amount)],
    target: Energy,
    pressure: Pressure,
    initial_temperature: Temperature,
    duration: Time,
    controls: IsobaricKineticsSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<IsobaricKineticsResult, IsobaricKineticsError> {
    poll(&mut cancelled)?;
    let amounts = feed_amounts(mechanism, feed)?;
    let mut reactor = IsobaricReactor::new(mechanism.clone(), controls, &mut cancelled)?;
    reactor
        .advance_adiabatic_enthalpy(
            &amounts,
            target,
            pressure,
            initial_temperature,
            duration,
            &mut cancelled,
        )?
        .to_owned(cancelled)
}
