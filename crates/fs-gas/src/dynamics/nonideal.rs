//! Conservative homogeneous real-gas fields on the original material owner.
//!
//! NASA calorics and Peng--Robinson mechanics supply one pressure, energy,
//! heat capacity and sound speed at each conserved density. One immutable
//! density/temperature reference defines the linear energy gauge for ALL cells:
//! `E = rho * (u(rho,T) - u(rho_ref,T_ref) + velocity^2/2)`.
//! The finite original-law increment cancels formation datums before rounding;
//! a density-dependent per-cell reference would violate energy conservation.
//!
//! The existing HLLC/HLL kernel, physical mesh, conservative update and balance
//! types are shared with the ideal-gas field, including smooth area variation
//! and the stationary side-wall pressure source. Only constitutive admission and
//! the bounded homogeneous field owner differ. First-order or source-admitted
//! MUSCL reconstruction and SSP-RK2 time stepping are an Estimated baseline. Local
//! source/stability checks do not establish a globally single-phase rectangle,
//! entropy stability, positivity for arbitrary EOS, or a shock-thickness model.

mod evolve;
pub use evolve::{RealGasBoundaries, RealGasBoundary, RealGasEuler1d};
mod thermal;
pub use thermal::{
    RealGasWallContact, RealGasWallOperator, RealGasWallReport, SharedWall, SharedWallControl,
    SharedWallGas, SharedWallLink, SharedWallLinkReport, SharedWallPort, SharedWallReport,
};
mod wall;
pub use wall::RealGasWallControl;

use super::{Conserved, Error, FaceFlux, Primitive, poll, positive};
use fs_matdb::QueryPoint;
use fs_material::MaterialError;
use fs_material::thermochemical_gas::energy::entropy::GasEntropyError;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::nonideal::{
    ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel,
};
use fs_qty::{Density, Temperature, Velocity};

/// Explicit density range and native caloric inversion work/temperature range.
///
/// The caller declares one continuous mechanically/calorically stable
/// homogeneous branch over this rectangle. Every actual state is checked by
/// the native material owner; this declaration does not perform a phase flash.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasDomain {
    /// Inclusive minimum/maximum mass densities [kg/m³].
    pub density_bounds_kg_m3: [f64; 2],
    /// Original source temperature bracket, absolute energy tolerance and work.
    pub energy: GasEnergySolve,
}

impl RealGasDomain {
    pub fn validate(self) -> Result<(), RealGasFlowError> {
        self.energy.validate()?;
        let [lo, hi] = self.density_bounds_kg_m3;
        if !positive(lo) || !positive(hi) || lo >= hi {
            return Err(RealGasFlowError::Invalid(
                "positive strictly increasing real-gas density bounds required",
            ));
        }
        Ok(())
    }

    fn admits_density(self, rho: f64) -> bool {
        rho.is_finite()
            && (self.density_bounds_kg_m3[0]..=self.density_bounds_kg_m3[1]).contains(&rho)
    }

    fn admits_temperature(self, t: f64) -> bool {
        t.is_finite()
            && (self.energy.minimum_temperature_k..=self.energy.maximum_temperature_k).contains(&t)
    }
}

/// Atomic material, numerical, work or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum RealGasFlowError {
    Invalid(&'static str),
    Material(MaterialError),
    Energy(GasEnergyError),
    /// Original density/entropy path for a physical wall rarefaction.
    Entropy(GasEntropyError),
    /// Original Rankine--Hugoniot pressure/energy/entropy wall admission.
    Shock(Box<crate::material::nonideal_shock::RealGasShockError>),
    /// Existing geometry, shared numerical flux or CFL admission refused.
    Numerical(Error),
    /// Original implicit gas/wall source refused before publication.
    Thermal(Box<crate::multicomponent::wall_exchange::WallExchangeError>),
    StepRetries {
        attempts: usize,
        last: Box<RealGasFlowError>,
    },
    Cancelled,
}

impl core::fmt::Display for RealGasFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "homogeneous real-gas field: {self:?}")
    }
}
impl std::error::Error for RealGasFlowError {}
impl From<MaterialError> for RealGasFlowError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<GasEnergyError> for RealGasFlowError {
    fn from(error: GasEnergyError) -> Self {
        match error {
            GasEnergyError::Cancelled => Self::Cancelled,
            other => Self::Energy(other),
        }
    }
}
impl From<GasEntropyError> for RealGasFlowError {
    fn from(error: GasEntropyError) -> Self {
        match error {
            GasEntropyError::Cancelled => Self::Cancelled,
            other => Self::Entropy(other),
        }
    }
}
impl From<crate::material::nonideal_shock::RealGasShockError> for RealGasFlowError {
    fn from(error: crate::material::nonideal_shock::RealGasShockError) -> Self {
        match error {
            crate::material::nonideal_shock::RealGasShockError::Cancelled => Self::Cancelled,
            other => Self::Shock(Box::new(other)),
        }
    }
}
impl From<Error> for RealGasFlowError {
    fn from(error: Error) -> Self {
        match error {
            Error::Cancelled => Self::Cancelled,
            other => Self::Numerical(other),
        }
    }
}
impl From<crate::multicomponent::wall_exchange::WallExchangeError> for RealGasFlowError {
    fn from(error: crate::multicomponent::wall_exchange::WallExchangeError) -> Self {
        match error {
            crate::multicomponent::wall_exchange::WallExchangeError::Cancelled => Self::Cancelled,
            other => Self::Thermal(Box::new(other)),
        }
    }
}

/// Original material model, fixed mixture/context and one common energy gauge.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGas {
    model: ThermochemicalPengRobinsonModel,
    context: QueryPoint,
    reference: ResolvedThermochemicalRealGasState,
    domain: RealGasDomain,
}

impl RealGas {
    /// Bind one immutable material mixture and a caller-declared homogeneous
    /// domain. The context temperature and explicit density choose the common
    /// reference; context pressure is replaced by the original EOS pressure.
    pub fn new(
        model: ThermochemicalPengRobinsonModel,
        context: QueryPoint,
        reference_density: Density,
        domain: RealGasDomain,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RealGasFlowError> {
        domain.validate()?;
        poll(&mut cancelled)?;
        if !domain.admits_density(reference_density.value()) {
            return Err(RealGasFlowError::Invalid(
                "reference density outside declared real-gas domain",
            ));
        }
        let reference = model.evaluate_density_temperature(&context, reference_density.value())?;
        if !domain.admits_temperature(reference.temperature_k()) {
            return Err(RealGasFlowError::Invalid(
                "reference temperature outside declared caloric interval",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            model,
            context,
            reference,
            domain,
        })
    }

    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        &self.model
    }
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }
    pub const fn reference_state(&self) -> &ResolvedThermochemicalRealGasState {
        &self.reference
    }
    pub const fn domain(&self) -> RealGasDomain {
        self.domain
    }
    /// Original absolute internal energy [J/kg] at the one retained reference.
    pub fn energy_reference(&self) -> f64 {
        self.reference.specific_internal_energy_j_kg()
    }

    /// Construct conservative data from explicit density, velocity and
    /// temperature. Pressure follows the original EOS, never an ideal-gas fit.
    pub fn conserved(
        &self,
        density: Density,
        velocity: Velocity,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Conserved, RealGasFlowError> {
        self.state_at(
            density.value(),
            velocity.value(),
            temperature.value(),
            &mut cancelled,
        )
        .map(|(q, _)| q)
    }

    pub(super) fn state_at(
        &self,
        density: f64,
        velocity: f64,
        temperature: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Conserved, Primitive), RealGasFlowError> {
        poll(cancelled)?;
        if !self.domain.admits_density(density)
            || !self.domain.admits_temperature(temperature)
            || !velocity.is_finite()
        {
            return Err(RealGasFlowError::Invalid(
                "primitive density/temperature outside real-gas domain or nonfinite velocity",
            ));
        }
        // Preserve all original axis descriptors, including the temperature kind.
        let point = self
            .context
            .clone()
            .with("temperature", temperature)
            .map_err(|error| MaterialError::State {
                what: format!("real-gas field temperature coordinate: {error}"),
            })?;
        let resolved = self.model.evaluate_density_temperature(&point, density)?;
        let relative_energy = self
            .model
            .specific_internal_energy_change_between_densities(
                self.reference.density_kg_m3(),
                self.reference.temperature_k(),
                density,
                temperature,
            )?;
        let q = Conserved {
            density,
            momentum: density * velocity,
            energy: density * (relative_energy + (0.5 * velocity) * velocity),
        };
        let recovered_velocity = q.momentum / density;
        let recovered_energy = q.energy / density - (0.5 * recovered_velocity) * recovered_velocity;
        if !q.finite()
            || !recovered_energy.is_finite()
            || (recovered_energy - relative_energy).abs() > self.domain.energy.energy_tolerance_j_kg
            || (relative_energy != 0.0 && recovered_energy == 0.0)
            || resolved.density_kg_m3().to_bits() != density.to_bits()
        {
            return Err(RealGasFlowError::Invalid(
                "conserved state cannot retain original density or declared caloric precision",
            ));
        }
        poll(cancelled)?;
        Ok((q, primitive(&resolved, recovered_velocity)))
    }

    /// Resolve a conservative cell through the original density/energy owner.
    pub fn thermodynamics(
        &self,
        state: Conserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalRealGasState, RealGasFlowError> {
        self.resolve(state, &mut cancelled)
    }

    pub fn primitive(
        &self,
        state: Conserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, RealGasFlowError> {
        self.decode(state, &mut cancelled)
    }

    /// Recover the physical energy density, including the original source datum.
    pub fn absolute_energy_density(&self, state: Conserved) -> Result<f64, RealGasFlowError> {
        let energy = state.energy + state.density * self.energy_reference();
        if !energy.is_finite() {
            return Err(RealGasFlowError::Invalid(
                "unrepresentable absolute energy density",
            ));
        }
        Ok(energy)
    }

    fn resolve(
        &self,
        state: Conserved,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalRealGasState, RealGasFlowError> {
        poll(cancelled)?;
        if !state.finite() || !self.domain.admits_density(state.density) {
            return Err(RealGasFlowError::Invalid(
                "conserved density outside real-gas domain or nonfinite state",
            ));
        }
        let velocity = state.momentum / state.density;
        let specific_total = state.energy / state.density;
        let kinetic = (0.5 * velocity) * velocity;
        let relative_energy = specific_total - kinetic;
        if !velocity.is_finite() || !relative_energy.is_finite() {
            return Err(RealGasFlowError::Invalid(
                "unrepresentable real-gas specific energy",
            ));
        }
        // At the exact reference density, removing kinetic energy can create
        // a tiny signed remainder even for data constructed at T_ref. The
        // native heat-increment API deliberately refuses a nonzero heat whose
        // temperature is unrepresentable. Here the input is a total-energy
        // cell, not an authored heat input: admit the retained reference only
        // inside BOTH the subtraction-roundoff scale and the caller's absolute
        // energy tolerance. The conserved energy is never altered or clipped.
        let roundoff = (8.0 * f64::EPSILON * specific_total.abs()
            + 8.0 * f64::EPSILON * kinetic.abs())
        .min(self.domain.energy.energy_tolerance_j_kg);
        if state.density.to_bits() == self.reference.density_kg_m3().to_bits()
            && relative_energy.abs() <= roundoff
        {
            poll(cancelled)?;
            return Ok(self.reference.clone());
        }
        Ok(self
            .model
            .reconstruct_density_internal_energy_increment_between_densities(
                &self.context,
                self.reference.density_kg_m3(),
                state.density,
                relative_energy,
                self.domain.energy,
                &mut *cancelled,
            )?
            .into_state())
    }

    fn decode(
        &self,
        state: Conserved,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, RealGasFlowError> {
        let resolved = self.resolve(state, cancelled)?;
        Ok(primitive(&resolved, state.momentum / state.density))
    }

    /// Original shared HLLC/HLL face algebra with actual real-gas acoustic
    /// speeds. Both HLLC star states and every HLL middle state must resolve
    /// through the same density/energy/material owner. Source-invalid stars
    /// trigger the bounded HLL fallback; cancellation always propagates.
    pub fn hllc_flux(
        &self,
        left: Conserved,
        right: Conserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        let pl = self.decode(left, &mut cancelled)?;
        let pr = self.decode(right, &mut cancelled)?;
        self.interface(left, pl, right, pr, &mut cancelled)
    }

    fn interface(
        &self,
        left: Conserved,
        pl: Primitive,
        right: Conserved,
        pr: Primitive,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        super::flux::checked_euler_flux(
            left,
            pl,
            right,
            pr,
            |q, cancelled| match self.decode(q, cancelled) {
                Ok(_) => Ok(true),
                Err(RealGasFlowError::Cancelled) => Err(RealGasFlowError::Cancelled),
                Err(_) => Ok(false),
            },
            cancelled,
        )
    }
}

fn primitive(state: &ResolvedThermochemicalRealGasState, velocity: f64) -> Primitive {
    Primitive {
        density: state.density_kg_m3(),
        velocity,
        pressure: state.pressure_pa(),
        temperature: state.temperature_k(),
        sound_speed: state.sound_speed_m_s(),
    }
}
