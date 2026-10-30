//! Source-owned fluid calorics in resolved three-dimensional thermoviscous flow.
//!
//! The original mixed FE owner resolves velocity, pressure and temperature.
//! This adapter binds its equilibrium and transport coefficients to one
//! `GasState` or a shared material fluid, evaluates original gas molecular reservoirs once per frequency,
//! and exports the same work-conjugate ports to the existing duct algebra.
//! There is no additional attenuation coefficient or fitted wall resistance.
//!
//! The FE fields use exp(+i omega t), while material calorics and ducts use
//! exp(-i omega t). Both changes of convention are explicit conjugations.
//! A traction shape psi means imposed stress -P psi n and inward flow
//! Q=-integral(psi v.n). Thus P is a total normal-traction coefficient, not a
//! point sample of thermodynamic pressure. For psi=1, Q is net volume flow.
//! Connecting nonuniform shapes requires the caller's matching dual traces;
//! these coordinates are not automatically evanescent traveling-wave modes.

use fs_duct::{DuctError, compose::TwoPort};
use fs_feec::acoustics::thermoviscous::time::{
    PreparedThermoviscousTime, ThermoviscousCaloricReservoir, ThermoviscousTimeError,
    ThermoviscousTimeLimits,
};
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousError, ThermoviscousFe, ThermoviscousField, ThermoviscousMedium,
    ThermoviscousThermodynamics,
};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_material::thermal::transport::ThermalTransportSample;
use fs_material::{MaterialError, gas::GasState, relaxation::RelaxationState};
use fs_math::c64::C64;

/// Finite solid-wall patches evaluated by the original material owner.
pub mod wall;

/// Retains the original material, finite-element and duct-composition refusals.
#[derive(Debug)]
pub enum ThermoviscousSourceError {
    Material(MaterialError),
    Interior(ThermoviscousError),
    Duct(DuctError),
    Time(ThermoviscousTimeError),
}

impl core::fmt::Display for ThermoviscousSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => e.fmt(f),
            Self::Interior(e) => e.fmt(f),
            Self::Duct(e) => e.fmt(f),
            Self::Time(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for ThermoviscousSourceError {}
impl From<MaterialError> for ThermoviscousSourceError {
    fn from(e: MaterialError) -> Self {
        Self::Material(e)
    }
}
impl From<ThermoviscousError> for ThermoviscousSourceError {
    fn from(e: ThermoviscousError) -> Self {
        Self::Interior(e)
    }
}
impl From<DuctError> for ThermoviscousSourceError {
    fn from(e: DuctError) -> Self {
        Self::Duct(e)
    }
}
impl From<ThermoviscousTimeError> for ThermoviscousSourceError {
    fn from(e: ThermoviscousTimeError) -> Self {
        Self::Time(e)
    }
}

fn invalid(what: &'static str) -> ThermoviscousSourceError {
    ThermoviscousError::InvalidInput(what).into()
}

/// Numerical view of one coherent source gas. Bulk viscosity is an explicit
/// independent dissipative channel [Pa s]; it must not include the molecular
/// reservoirs also supplied through `from_relaxation`. No value is inferred.
pub fn source_medium(
    gas: &GasState,
    bulk_viscosity_pa_s: f64,
) -> Result<ThermoviscousMedium, ThermoviscousSourceError> {
    gas.validate_coherence()?;
    if !bulk_viscosity_pa_s.is_finite() || bulk_viscosity_pa_s < 0.0 {
        return Err(invalid(
            "resolved gas requires finite nonnegative bulk viscosity",
        ));
    }
    Ok(ThermoviscousMedium {
        density_kg_m3: gas.density,
        pressure_pa: gas.pressure,
        temperature_k: gas.temperature,
        cp_j_kg_k: gas.specific_heat_cp,
        dynamic_viscosity_pa_s: gas.dynamic_viscosity,
        thermal_conductivity_w_m_k: gas.thermal_conductivity,
        bulk_viscosity_pa_s,
    })
}

/// The actual EOS derivatives and transport of one generic material fluid.
/// `pressure_pa` is the declared equilibrium pressure (metadata for this local
/// linearization), never an effective ideal-gas pressure. Bulk viscosity is
/// an explicit independent input; no relaxation law is invented.
///
/// Cp, Cv and adiabatic c fix beta_T=(Cp/Cv)/(rho*c^2). Signed alpha_p must
/// come from the SAME source's thermal view. Its thermodynamic identity
/// beta_S=beta_T-alpha_p^2*T/(rho*Cp) must agree to floating-point accuracy.
/// Unknown expansion or independently inconsistent property tuples refuse;
/// the sign cannot be recovered from Cp/Cv and is never guessed.
pub fn source_fluid_medium<M: ThermoviscousFluid + ThermalFluidTransport + ?Sized>(
    material: &M,
    pressure_pa: f64,
    bulk_viscosity_pa_s: f64,
) -> Result<(ThermoviscousMedium, ThermoviscousThermodynamics), ThermoviscousSourceError> {
    let fluid = material.thermoviscous_state()?;
    let thermal = material.thermal_transport_state()?;
    if fluid.transport() != thermal.transport() {
        return Err(invalid(
            "resolved fluid acoustic and thermal transport must describe one state",
        ));
    }
    let alpha = thermal.volumetric_expansion_per_k().ok_or_else(|| {
        invalid("resolved fluid acoustics requires known signed isobaric expansion")
    })?;
    let acoustic = fluid.acoustic();
    let transport = fluid.transport();
    let medium = ThermoviscousMedium {
        density_kg_m3: acoustic.density_kg_m3(),
        pressure_pa,
        temperature_k: thermal.temperature_k(),
        cp_j_kg_k: transport.specific_heat_j_kg_k(),
        dynamic_viscosity_pa_s: transport.dynamic_viscosity_pa_s(),
        thermal_conductivity_w_m_k: transport.thermal_conductivity_w_m_k(),
        bulk_viscosity_pa_s,
    };
    let beta_s = 1.0 / acoustic.adiabatic_bulk_modulus_pa();
    let thermodynamics = ThermoviscousThermodynamics {
        isothermal_compressibility_pa_inv: fluid.heat_capacity_ratio() * beta_s,
        isobaric_expansion_per_k: alpha,
    };
    let actual = thermodynamics.adiabatic_compressibility_pa_inv(medium)?;
    if (actual - beta_s).abs() > 256.0 * f64::EPSILON * beta_s {
        return Err(invalid(
            "resolved fluid Cp/Cv, sound speed and signed expansion violate one thermodynamic closure",
        ));
    }
    Ok((medium, thermodynamics))
}

/// Resolved mixed FE with an original source fluid and declared frequency band.
/// The band is a computation domain, not a source-accuracy certificate or a
/// substitute for resolving the thermal and viscous penetration depths.
pub struct ThermoviscousInterior<'c, S = GasState> {
    fe: ThermoviscousFe<'c>,
    source: S,
    relaxation: Option<RelaxationState>,
    band_hz: [f64; 2],
    wall_sources: Vec<ThermalTransportSample>,
}

impl<'c> ThermoviscousInterior<'c> {
    /// Bind an original FE whose numerical coefficients came from this gas.
    /// All equilibrium and transport coefficients must agree bitwise.
    pub fn from_gas(
        fe: ThermoviscousFe<'c>,
        gas: GasState,
        band_hz: [f64; 2],
    ) -> Result<Self, ThermoviscousSourceError> {
        if !fe.conductive_walls().is_empty() {
            return Err(invalid(
                "bind the source gas before attaching source-evaluated conductive walls",
            ));
        }
        let actual = fe.medium();
        let expected = source_medium(&gas, actual.bulk_viscosity_pa_s)?;
        if fe.thermodynamics() != ThermoviscousThermodynamics::ideal_gas(expected)? {
            return Err(invalid(
                "resolved gas cannot bind a different fluid equation of state",
            ));
        }
        for (a, b) in [
            (actual.density_kg_m3, expected.density_kg_m3),
            (actual.pressure_pa, expected.pressure_pa),
            (actual.temperature_k, expected.temperature_k),
            (actual.cp_j_kg_k, expected.cp_j_kg_k),
            (
                actual.dynamic_viscosity_pa_s,
                expected.dynamic_viscosity_pa_s,
            ),
            (
                actual.thermal_conductivity_w_m_k,
                expected.thermal_conductivity_w_m_k,
            ),
        ] {
            if a.to_bits() != b.to_bits() {
                return Err(invalid(
                    "resolved thermoviscous FE differs from its source gas",
                ));
            }
        }
        if !band_hz.iter().all(|f| f.is_finite() && *f > 0.0) || band_hz[1] <= band_hz[0] {
            return Err(invalid(
                "resolved gas needs a positive increasing frequency band",
            ));
        }
        Ok(Self {
            fe,
            source: gas,
            relaxation: None,
            band_hz,
            wall_sources: Vec::new(),
        })
    }

    /// Retain the original evaluated reservoirs. Dynamic heat capacity enters
    /// the temperature equation itself, so molecular and conduction losses use
    /// the same solved temperature field rather than independent loss overlays.
    pub fn from_relaxation(
        fe: ThermoviscousFe<'c>,
        source: RelaxationState,
        band_hz: [f64; 2],
    ) -> Result<Self, ThermoviscousSourceError> {
        let mut model = Self::from_gas(fe, *source.gas(), band_hz)?;
        model.relaxation = Some(source);
        Ok(model)
    }

    pub fn gas(&self) -> &GasState {
        &self.source
    }
}

impl<'c, S: ThermoviscousFluid + ThermalFluidTransport> ThermoviscousInterior<'c, S> {
    /// Retain the original material source beside its full mixed FE. All
    /// coefficients and EOS derivatives must match exactly, and neither wall
    /// sources nor caloric reservoirs can be silently inherited from a gas.
    pub fn from_material(
        fe: ThermoviscousFe<'c>,
        source: S,
        band_hz: [f64; 2],
    ) -> Result<Self, ThermoviscousSourceError> {
        if !fe.conductive_walls().is_empty() {
            return Err(invalid(
                "bind the original fluid before its source-evaluated conductive walls",
            ));
        }
        let actual = fe.medium();
        let (expected, thermodynamics) =
            source_fluid_medium(&source, actual.pressure_pa, actual.bulk_viscosity_pa_s)?;
        if actual != expected || fe.thermodynamics() != thermodynamics {
            return Err(invalid(
                "resolved thermoviscous FE differs from its original fluid source",
            ));
        }
        if !band_hz.iter().all(|f| f.is_finite() && *f > 0.0) || band_hz[1] <= band_hz[0] {
            return Err(invalid(
                "resolved fluid needs a positive increasing frequency band",
            ));
        }
        Ok(Self {
            fe,
            source,
            relaxation: None,
            band_hz,
            wall_sources: Vec::new(),
        })
    }
}

impl<'c, S> ThermoviscousInterior<'c, S> {
    pub fn inner(&self) -> &ThermoviscousFe<'c> {
        &self.fe
    }
    /// Original material/source owner; the immutable numerical snapshot lives
    /// in `inner()`. A source change requires an explicit new assembly.
    pub fn material_source(&self) -> &S {
        &self.source
    }
    pub fn relaxation(&self) -> Option<&RelaxationState> {
        self.relaxation.as_ref()
    }
    pub const fn band_hz(&self) -> [f64; 2] {
        self.band_hz
    }
    /// Original source samples in the same order as the finite wall regions.
    /// Their conductivity direction and caloric constraint are selected when
    /// preparing the material, not inferred by the acoustic adapter.
    pub fn wall_sources(&self) -> &[ThermalTransportSample] {
        &self.wall_sources
    }

    /// Causal full-field motion and work with this SAME immutable fluid source.
    /// The borrow retains the source, FE and structural model for every step.
    /// Original molecular capacities and rates own full P2 internal-temperature
    /// fields in the same solve. Finite-wall impedances still require their
    /// actual time states and are refused.
    /// The clock must sample the declared source band; no frequency truncation
    /// or proof of spatial/input-waveform resolution follows from that check.
    pub fn prepare_time<'a>(
        &'a self,
        duration_s: f64,
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousTime<'a, 'c>, ThermoviscousSourceError> {
        if cancelled() {
            return Err(ThermoviscousError::Cancelled.into());
        }
        if duration_s * self.band_hz[1] > 0.5 {
            return Err(invalid(
                "resolved transient clock does not sample the declared source band",
            ));
        }
        match &self.relaxation {
            Some(source) => {
                let reservoirs: Vec<_> = source
                    .processes()
                    .iter()
                    .map(|process| ThermoviscousCaloricReservoir {
                        heat_capacity_j_kg_k: process.heat_capacity_j_kg_k,
                        relaxation_time_s: process.relaxation_time_s,
                    })
                    .collect();
                Ok(self.fe.prepare_time_with_caloric_reservoirs(
                    duration_s,
                    source.frozen_cv_j_kg_k() + source.gas().specific_gas_constant,
                    &reservoirs,
                    limits,
                    cancelled,
                )?)
            }
            None => Ok(self.fe.prepare_time(duration_s, limits, cancelled)?),
        }
    }

    fn heat_capacity(
        &self,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<C64, ThermoviscousSourceError> {
        if cancelled() {
            return Err(invalid("resolved source-fluid evaluation cancelled"));
        }
        if !frequency_hz.is_finite()
            || frequency_hz < self.band_hz[0]
            || frequency_hz > self.band_hz[1]
        {
            return Err(invalid(
                "resolved source-fluid frequency outside declared band",
            ));
        }
        match &self.relaxation {
            Some(source) => Ok(source
                .caloric_response(core::f64::consts::TAU * frequency_hz)?
                .cp_j_kg_k
                .conj()),
            None => Ok(C64::from_re(self.fe.medium().cp_j_kg_k)),
        }
    }

    /// Simultaneous traction excitation with all fields and work terms from
    /// the original mixed solve. Every pressure coefficient has units Pa.
    pub fn response(
        &self,
        frequency_hz: f64,
        traction_pa: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousField, ThermoviscousSourceError> {
        let cp = self.heat_capacity(frequency_hz, cancelled)?;
        Ok(self
            .fe
            .prepare_with_heat_capacity(frequency_hz, cp, cancelled)?
            .response(traction_pa, cancelled)?)
    }

    /// Simultaneous source-bound fluid traction and original structural force.
    /// The FE must have been assembled with its original modal moving wall.
    /// Forces use N/sqrt(kg) and both drives use exp(+i omega t); gas molecular
    /// calorics, material fluid EOS and finite wall heat remain in this same solve.
    pub fn response_with_wall_forces(
        &self,
        frequency_hz: f64,
        traction_pa: &[C64],
        wall_forces_n_per_sqrt_kg: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousField, ThermoviscousSourceError> {
        let cp = self.heat_capacity(frequency_hz, cancelled)?;
        Ok(self
            .fe
            .prepare_with_heat_capacity(frequency_hz, cp, cancelled)?
            .response_with_wall_forces(traction_pa, wall_forces_n_per_sqrt_kg, cancelled)?)
    }

    /// Row-major Q=Y P [m3/(Pa s)] in the FE's exp(+i omega t) convention.
    /// The original owner factors once for all excitation columns.
    pub fn admittance(
        &self,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, ThermoviscousSourceError> {
        let cp = self.heat_capacity(frequency_hz, cancelled)?;
        Ok(self
            .fe
            .admittance_with_heat_capacity(frequency_hz, cp, cancelled)?)
    }

    /// Direct power scattering in exp(-i omega t), ready for the original
    /// Redheffer composition. References are positive real pressure/flow
    /// resistances. Connected coordinates, shapes and references must match.
    pub fn duct_scattering(
        &self,
        frequency_hz: f64,
        reference_resistance: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, ThermoviscousSourceError> {
        if self.fe.port_mode_count() == 0 {
            return Err(invalid(
                "a closed moving-wall cavity has no duct scattering coordinates",
            ));
        }
        let cp = self.heat_capacity(frequency_hz, cancelled)?;
        let mut s = self.fe.scattering_with_heat_capacity(
            frequency_hz,
            cp,
            reference_resistance,
            cancelled,
        )?;
        for row in s.chunks_mut(self.fe.port_mode_count()) {
            if cancelled() {
                return Err(invalid("resolved thermoviscous duct conversion cancelled"));
            }
            for value in row {
                *value = value.conj();
            }
        }
        Ok(s)
    }

    /// Exactly two retained work coordinates in the original duct `TwoPort`.
    /// No modes are dropped and no transfer-matrix inversion is introduced.
    pub fn duct_two_port(
        &self,
        frequency_hz: f64,
        reference_impedance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TwoPort, ThermoviscousSourceError> {
        if self.fe.port_mode_count() != 2 {
            return Err(invalid(
                "resolved duct TwoPort needs exactly two port coordinates",
            ));
        }
        let s = self.duct_scattering(frequency_hz, &[reference_impedance; 2], cancelled)?;
        Ok(TwoPort::new([s[0], s[1], s[2], s[3]], reference_impedance)?)
    }
}

#[cfg(test)]
mod tests;
