//! Bounded material-density feedback without changing conserved VOF mass.

use super::*;
use crate::material::transport::boussinesq::{BoussinesqConfig, BoussinesqMaterial};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

#[derive(Clone, Copy)]
pub(super) struct BoussinesqLiquid {
    config: BoussinesqConfig<3>,
    reference: FluidTransportState,
}

impl BoussinesqLiquid {
    pub(super) fn try_new(
        flow: &MaterialFreeSurface3<'_>,
        law: &LiquidTemperatureLaw<'_>,
        config: BoussinesqConfig<3>,
    ) -> Result<Self, Error> {
        config
            .validate()
            .map_err(|_| Error::Invalid("invalid liquid Boussinesq bounds"))?;
        let reference = law
            .reference_state()
            .transport_state()
            .map_err(MaterialLbmError::from)?;
        let existing = flow.lattice.reference_force();
        for axis in 0..3 {
            if flow.scaling.gravity_lattice(config.gravity_m_s2[axis])? != existing[axis] {
                return Err(Error::Invalid(
                    "Boussinesq gravity must equal the original free-surface reference gravity",
                ));
            }
        }
        Ok(Self { config, reference })
    }

    pub(super) fn force(
        self,
        flow: &MaterialFreeSurface3<'_>,
        state: &ResolvedLiquidState,
    ) -> Result<[f64; 3], Error> {
        let departure = (state.density_kg_m3() - self.reference.density_kg_m3())
            / self.reference.density_kg_m3();
        let mut total = flow.lattice.reference_force();
        for axis in 0..3 {
            let acceleration = departure * self.config.gravity_m_s2[axis];
            let perturbation = flow.scaling.gravity_lattice(acceleration)?;
            if acceleration != 0.0 && perturbation == 0.0 {
                return Err(Error::Invalid(
                    "liquid buoyancy force underflows lattice units",
                ));
            }
            total[axis] += perturbation;
            if !total[axis].is_finite() {
                return Err(Error::Invalid(
                    "liquid buoyancy force overflows lattice units",
                ));
            }
        }
        Ok(total)
    }
}

/// Resolve once from the original h and source owner. Conduction and native
/// momentum share this admission, including isolated heated cells and ghosts.
pub(super) fn resolve_liquid(
    law: &LiquidTemperatureLaw<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    density: f64,
    buoyancy: Option<BoussinesqLiquid>,
) -> Result<(ResolvedLiquidState, FluidTransportState), Error> {
    let temperature = liquid_temperature(caloric, h, density)?;
    let state = law.at_temperature(Temperature::new(temperature))?;
    let transport = liquid_transport(&state, density, buoyancy)?;
    Ok((state, transport))
}

pub(super) fn liquid_temperature(
    caloric: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    density: f64,
) -> Result<f64, Error> {
    let phase = caloric.state_at_specific_enthalpy(h)?;
    if phase.phase() != SolidLiquidPhase::Liquid
        || phase.bulk_density_kg_m3().to_bits() != density.to_bits()
    {
        return Err(Error::Invalid(
            "free-surface heat needs pure liquid at the original fixed caloric density",
        ));
    }
    Ok(phase.temperature_k())
}

pub(super) fn liquid_transport(
    state: &ResolvedLiquidState,
    density: f64,
    buoyancy: Option<BoussinesqLiquid>,
) -> Result<FluidTransportState, Error> {
    let actual = state.transport_state().map_err(MaterialLbmError::from)?;
    let transport = if let Some(binding) = buoyancy {
        // The original source and its receipts remain in `state`; this is the
        // shared numerical projection of that SAME already-resolved tuple.
        *BoussinesqMaterial::try_new(actual, &binding.reference, binding.config)
            .map_err(MaterialLbmError::from)?
            .reference_storage_transport()
    } else {
        if state.density_kg_m3().to_bits() != density.to_bits() {
            return Err(Error::Invalid(
                "temperature-dependent density needs an admitted reference-mass projection or volume-change owner",
            ));
        }
        actual
    };
    Ok(transport)
}

impl<'s> MaterialFreeSurface3<'s> {
    /// Bind the original liquid h chart with source-density-driven buoyancy.
    ///
    /// VOF mass, cell volume and the caloric chart retain the original reference
    /// density. The shared Boussinesq owner bounds actual source rho/Cp departures;
    /// local dynamic viscosity divided by rho0 supplies collision tau. The
    /// original h chart still owns thermal storage, including its actual slope.
    /// `material_at` retains the unprojected source density and all receipts.
    ///
    /// Gravity must match the original domain. Its rho0*g reference load stays
    /// present at the free surface; only (rho(T)-rho0)*g is added. No ideal-gas
    /// expansion coefficient or material-name rule replaces the density source.
    /// Attachment changes only the Guo population frame, preserving initial
    /// physical velocity and force-corrected stress even for an evolved flow.
    /// Each step holds the post-thermal force during native momentum, then
    /// rebinds the endpoint source frame without an extra half-force impulse.
    ///
    /// This is a first-order held-force thermal/flow split, not a moving-volume
    /// approximation. It adds no phase conversion, pressure work, gravitational
    /// energy balance or physical viscous heating. Original stationary walls,
    /// constant sourced capillarity, native pressure-boundary accuracy and low-Mach limits remain.
    pub fn with_boussinesq_enthalpy<'m>(
        self,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        config: BoussinesqConfig<3>,
        initial_h: impl Fn([usize; 3]) -> f64,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        self.with_enthalpy_model(law, caloric, controls, initial_h, Some(config), None, None, gate)
    }
}
