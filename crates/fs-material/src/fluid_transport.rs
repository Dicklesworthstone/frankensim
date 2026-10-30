//! One checked transport snapshot shared by gas and liquid consumers.
//!
//! This adapts existing material evaluations; it is not another material
//! database, equation of state, or phase selector. Four independent SI values
//! determine momentum diffusion, thermal diffusion and transport groups at ONE
//! state. Source receipts remain on the original input; this numerical view
//! alone does not convey source or constitutive-model authority.

use crate::MaterialError;
use crate::gas::{GasState, ResolvedGasState};
use crate::gas_mixture::ResolvedGasMixtureState;
use crate::liquid::ResolvedLiquidState;
use crate::moist_air::ResolvedMoistAirState;

mod thermal;
pub use thermal::{ThermalFluidTransport, ThermalFluidTransportState};

/// Immutable transport coefficients for one evaluated material state point.
/// Derived quantities cannot be mutated independently of their inputs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluidTransportState {
    density: f64,
    specific_heat: f64,
    conductivity: f64,
    viscosity: f64,
    kinematic_viscosity: f64,
    thermal_diffusivity: f64,
    prandtl: f64,
}

impl FluidTransportState {
    /// Admit coherent-SI `(rho [kg/m³], cp [J/(kg K)], k [W/(m K)], mu [Pa s])`.
    /// Raw values carry no material-card authority. This is Newtonian,
    /// single-state transport, not a constant-property claim over an evolution.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite inputs and unrepresentable derived values.
    pub fn try_new(
        density: f64,
        specific_heat: f64,
        conductivity: f64,
        viscosity: f64,
    ) -> Result<Self, MaterialError> {
        if ![density, specific_heat, conductivity, viscosity]
            .into_iter()
            .all(positive)
        {
            return Err(invalid("fluid density, heat capacity, conductivity and viscosity must be finite and positive"));
        }
        let kinematic_viscosity = viscosity / density;
        let thermal_diffusivity = (conductivity / density) / specific_heat;
        let prandtl = kinematic_viscosity / thermal_diffusivity;
        if ![kinematic_viscosity, thermal_diffusivity, prandtl]
            .into_iter()
            .all(positive)
        {
            return Err(invalid("derived fluid transport is not representable"));
        }
        Ok(Self {
            density,
            specific_heat,
            conductivity,
            viscosity,
            kinematic_viscosity,
            thermal_diffusivity,
            prandtl,
        })
    }

    /// Density [kg/m³].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 { self.density }
    /// Isobaric specific heat [J/(kg K)].
    #[must_use]
    pub const fn specific_heat_j_kg_k(&self) -> f64 { self.specific_heat }
    /// Thermal conductivity [W/(m K)].
    #[must_use]
    pub const fn thermal_conductivity_w_m_k(&self) -> f64 { self.conductivity }
    /// Dynamic viscosity [Pa s].
    #[must_use]
    pub const fn dynamic_viscosity_pa_s(&self) -> f64 { self.viscosity }
    /// Kinematic viscosity `nu = mu/rho` [m²/s].
    #[must_use]
    pub const fn kinematic_viscosity_m2_s(&self) -> f64 { self.kinematic_viscosity }
    /// Thermal diffusivity `alpha = k/(rho cp)` [m²/s].
    #[must_use]
    pub const fn thermal_diffusivity_m2_s(&self) -> f64 { self.thermal_diffusivity }
    /// `Pr = nu/alpha = mu cp/k`.
    #[must_use]
    pub const fn prandtl(&self) -> f64 { self.prandtl }

    /// Reynolds number for nonnegative speed [m/s] and positive length [m].
    ///
    /// # Errors
    /// Refuses invalid inputs and an unrepresentable dimensionless group.
    pub fn reynolds(&self, speed_m_s: f64, length_m: f64) -> Result<f64, MaterialError> {
        dimensionless_transport(speed_m_s, length_m, self.kinematic_viscosity)
    }

    /// Peclet number for nonnegative speed [m/s] and positive length [m].
    ///
    /// # Errors
    /// Refuses invalid inputs and an unrepresentable dimensionless group.
    pub fn peclet(&self, speed_m_s: f64, length_m: f64) -> Result<f64, MaterialError> {
        dimensionless_transport(speed_m_s, length_m, self.thermal_diffusivity)
    }

    /// Signed Rayleigh number `g beta delta_T L³/(nu alpha)`.
    /// Gravity is a nonnegative magnitude [m/s²]. Volumetric expansion [1/K]
    /// and wall-minus-reference temperature difference [K] remain explicit:
    /// liquid expansion is NOT inferred as `1/T`. A correlation may further
    /// refuse a nonpositive Rayleigh number.
    ///
    /// # Errors
    /// Refuses invalid inputs and an unrepresentable nonzero driving term.
    pub fn rayleigh(
        &self,
        gravity_m_s2: f64,
        expansion_per_k: f64,
        delta_temperature_k: f64,
        length_m: f64,
    ) -> Result<f64, MaterialError> {
        if !gravity_m_s2.is_finite() || gravity_m_s2 < 0.0
            || !expansion_per_k.is_finite() || !delta_temperature_k.is_finite()
            || !positive(length_m)
        {
            return Err(invalid("Rayleigh inputs must be finite, with nonnegative gravity and positive length"));
        }
        if gravity_m_s2 == 0.0 || expansion_per_k == 0.0 || delta_temperature_k == 0.0 {
            return Ok(0.0);
        }
        let ra = gravity_m_s2 * expansion_per_k * delta_temperature_k
            * (length_m / self.kinematic_viscosity)
            * (length_m / self.thermal_diffusivity) * length_m;
        if !ra.is_finite() || ra == 0.0 {
            return Err(invalid("Rayleigh number is not representable"));
        }
        Ok(ra)
    }
}

/// Shared view of one material evaluation. The original input remains the
/// owner of source receipts, query conditions and constitutive-model choices.
pub trait FluidTransport {
    /// Checked numerical snapshot at the input's evaluated state point.
    ///
    /// # Errors
    /// Refuses missing or inadmissible transport coefficients. Gas inputs also
    /// require the SAME coherent EOS, caloric and acoustic state; stale public
    /// fields are not silently reinterpreted as an unrelated raw fluid tuple.
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError>;
}

/// Borrow the original owner rather than copying its coefficients into another
/// material. Supports concrete sources and heterogeneous `&dyn FluidTransport`
/// inputs; the selected source's complete admission and errors are unchanged.
impl<T: FluidTransport + ?Sized> FluidTransport for &T {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        T::transport_state(*self)
    }
}
/// The borrowed thermal view delegates to the SAME source, including unknown
/// expansion and refusal of incoherent temperature or transport state.
impl<T: ThermalFluidTransport + ?Sized> ThermalFluidTransport for &T {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        T::thermal_transport_state(*self)
    }
}

impl FluidTransport for FluidTransportState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> { Ok(*self) }
}

impl FluidTransport for GasState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        // One admission and one coefficient derivation for flow, heat and
        // acoustics. Raw transport-only tuples still have their explicit
        // FluidTransportState constructor; a GasState must remain a gas state.
        self.coherent_transport()
    }
}

impl FluidTransport for ResolvedGasState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.state().transport_state()
    }
}
impl FluidTransport for ResolvedMoistAirState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.state().transport_state()
    }
}
impl FluidTransport for ResolvedGasMixtureState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.state().transport_state()
    }
}
impl FluidTransport for ResolvedLiquidState {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        FluidTransportState::try_new(
            self.density_kg_m3(), self.specific_heat_j_kg_k(),
            self.thermal_conductivity_w_m_k(), self.dynamic_viscosity_pa_s(),
        )
    }
}

fn dimensionless_transport(speed: f64, length: f64, diffusivity: f64) -> Result<f64, MaterialError> {
    if !speed.is_finite() || speed < 0.0 || !positive(length) {
        return Err(invalid("transport groups need finite nonnegative speed and positive length"));
    }
    if speed == 0.0 { return Ok(0.0); }
    let group = (speed / diffusivity) * length;
    if !positive(group) {
        return Err(invalid("transport dimensionless group is not representable"));
    }
    Ok(group)
}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn invalid(what: &str) -> MaterialError { MaterialError::Parameters { what: what.into() } }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gas::GasSpec;
    use fs_qty::semantic::{Composition, CompositionBasis};

    #[test]
    fn transport_identity_and_derived_groups_close() {
        let state = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        assert_eq!(state.kinematic_viscosity_m2_s(), 1e-6);
        assert_eq!(state.thermal_diffusivity_m2_s(), 1.25e-7);
        assert_eq!(state.prandtl(), 8.0);
        assert!((state.peclet(2.0, 0.1).unwrap() / state.reynolds(2.0, 0.1).unwrap() - 8.0).abs() < 1e-12);
        assert_eq!(state.reynolds(0.0, 1.0).unwrap(), 0.0);
    }

    #[test]
    fn gas_and_liquid_use_the_same_view_without_air_defaults() {
        let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
        let gas = air.transport_state().unwrap();
        assert_eq!(gas.dynamic_viscosity_pa_s(), air.dynamic_viscosity);
        let liquid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        assert!(liquid.reynolds(1.0, 0.1).unwrap() > gas.reynolds(1.0, 0.1).unwrap());
    }

    #[test]
    fn pressure_changes_both_momentum_and_thermal_diffusion() {
        let spec = GasSpec::dry_air_ussa1976();
        let a = GasState::try_new(&spec, 300.0, 100000.0).unwrap().transport_state().unwrap();
        let b = GasState::try_new(&spec, 300.0, 200000.0).unwrap().transport_state().unwrap();
        assert_eq!(b.kinematic_viscosity_m2_s(), a.kinematic_viscosity_m2_s() / 2.0);
        assert_eq!(b.thermal_diffusivity_m2_s(), a.thermal_diffusivity_m2_s() / 2.0);
        assert_eq!(b.prandtl(), a.prandtl());
    }

    #[test]
    fn arbitrary_composition_uses_the_existing_gas_mixture() {
        let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.98, 0.02]).unwrap();
        let gas = GasState::try_new_mixture(
            &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()],
            &composition, 300.0, 100000.0, Some(1),
        ).unwrap();
        let state = gas.transport_state().unwrap();
        assert_eq!(state.density_kg_m3(), gas.density);
        assert_eq!(state.thermal_conductivity_w_m_k(), gas.thermal_conductivity);
        assert!((state.prandtl() / gas.prandtl - 1.0).abs() < 1e-12);
    }

    #[test]
    fn incoherent_gas_is_not_silently_reinterpreted_as_raw_transport() {
        let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
        let expected = gas.transport_state().unwrap();
        gas.prandtl = 100.0;
        assert!(gas.transport_state().is_err());
        // A caller declaring only a transport tuple may use that explicit type,
        // but it has neither gas EOS semantics nor an acoustic-state claim.
        assert_eq!(FluidTransportState::try_new(gas.density, gas.specific_heat_cp,
            gas.thermal_conductivity, gas.dynamic_viscosity).unwrap(), expected);
    }

    #[test]
    fn invalid_materials_and_overflows_refuse() {
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(FluidTransportState::try_new(value, 1000.0, 1.0, 0.001).is_err());
            assert!(FluidTransportState::try_new(1.0, value, 1.0, 0.001).is_err());
            assert!(FluidTransportState::try_new(1.0, 1000.0, value, 0.001).is_err());
            assert!(FluidTransportState::try_new(1.0, 1000.0, 1.0, value).is_err());
        }
        assert!(FluidTransportState::try_new(f64::MIN_POSITIVE, f64::MIN_POSITIVE, f64::MAX, 1.0).is_err());
        let state = FluidTransportState::try_new(1.0, 1000.0, 1.0, 0.001).unwrap();
        assert!(state.reynolds(-1.0, 1.0).is_err());
        assert!(state.reynolds(1.0, 0.0).is_err());
        assert!(state.peclet(f64::MAX, f64::MAX).is_err());
    }

    #[test]
    fn thermal_expansion_is_explicit_and_sign_is_not_discarded() {
        let state = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let hot = state.rayleigh(9.81, 2e-4, 10.0, 0.1).unwrap();
        let cold = state.rayleigh(9.81, 2e-4, -10.0, 0.1).unwrap();
        assert_eq!(hot, -cold);
        assert_eq!(state.rayleigh(9.81, 0.0, 10.0, 0.1).unwrap(), 0.0);
        assert!(state.rayleigh(-9.81, 2e-4, 10.0, 0.1).is_err());
    }
}
