//! Existing transverse no-slip momentum coupled to exact DG pressure ports.
//!
//! A fixed passage is distinct from the moving aperture and from either duct.
//! Its original finite-volume law is dp = Z Q - H. The shared series owner
//! eliminates this law together with both optional tract and bore relations,
//! then advances every participant at the same accepted through-flow. No
//! fitted inertance, hydraulic-diameter substitution or delayed pressure enters.

use super::{ApertureDg, ApertureDgLoad, invalid, series_inlet_state};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::paired::BladePropagator;
use crate::bernoulli_aperture::slit_channel::FixedChannel;
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use fs_material::gas::GasState;

/// Admit the complete inlet material and SUM of the exact discrete pressure
/// loads before constructing a single, paired or material-modal valve.
///
/// Both original acoustic owners must be cold. A pre-existing channel velocity
/// profile is permitted and retains its explicit initial kinetic energy. The
/// complete channel gas and clock must equal the inlet's; StaticPressure is
/// required because a kinetic-head exit is not an affine pressure port.
/// Geometry must not include this same fixed passage a second time in a duct
/// or in the moving slit. The caller owns that physical assembly declaration.
pub fn inlet_state_with_supply(
    upstream: Option<&ApertureDgLoad>,
    downstream: &ApertureDgLoad,
    channel: &FixedChannel,
) -> Result<(GasState, f64), AcousticRealizeError> {
    let (gas, acoustic) = match upstream {
        Some(upstream) => series_inlet_state(upstream, downstream)?,
        None => {
            if downstream.has_driven_mouth() || downstream.state().step_count() != 0 {
                return Err(invalid(
                    "DG supply requires a cold acoustic load without a driven exterior mouth",
                ));
            }
            (
                *downstream.duct().inlet_gas(),
                downstream.impedance_pa_s_m3(),
            )
        }
    };
    if channel.gas() != &gas
        || channel.sample_period_s().to_bits() != downstream.propagator().dt().to_bits()
    {
        return Err(invalid(
            "DG supply passage and valve inlets must share one complete material gas and clock",
        ));
    }
    let impedance = acoustic + channel.pressure_load()?.impedance_pa_s_m3;
    if !impedance.is_finite() || impedance <= acoustic {
        return Err(invalid(
            "DG supply impedance cannot be represented beside its acoustic loads",
        ));
    }
    Ok((gas, impedance))
}

impl ApertureDg<LoadedAperture> {
    /// Bind the original valve, optional upstream tract, bore and fixed no-slip
    /// passage on one physical clock. Construct the valve with the gas and
    /// impedance returned by `inlet_state_with_supply`; no source is retuned.
    ///
    /// The passage can be planar or circular, with a refinable transverse mesh.
    /// Its momentum survives attacks, releases and flow reversals. Original
    /// viscous heat and backward-Euler numerical dissipation stay separate in
    /// each frame's `supply` account and enter the total balance exactly once.
    /// All nonlinear, acoustic and channel candidates pass their pressure/work
    /// checks before any participant publishes a step.
    ///
    /// This is fixed-geometry, fully developed incompressible channel flow at
    /// a prescribed complete material state. Axial development, moving-channel
    /// squeeze flow, turbulence and channel heating are not supplied. Existing
    /// finite-background/live-material operations refuse this attachment until
    /// they include its separate inventory and parameter-work account.
    pub fn new_loaded_with_supply(
        aperture: impl Into<LoadedAperture>,
        upstream: Option<ApertureDgLoad>,
        load: ApertureDgLoad,
        channel: impl Into<FixedChannel>,
    ) -> Result<Self, AcousticRealizeError> {
        let channel = channel.into();
        let (gas, impedance) = inlet_state_with_supply(upstream.as_ref(), &load, &channel)?;
        let aperture = aperture.into();
        let (density, declared_impedance) = aperture.load();
        if aperture.accepted_steps() != 0
            || aperture.blade_propagator() != BladePropagator::ExactZoh
            || aperture.time_step_s().to_bits() != load.propagator().dt().to_bits()
            || density.to_bits() != gas.density.to_bits()
            || declared_impedance.to_bits() != impedance.to_bits()
            || matches!(&aperture, LoadedAperture::Paired(a) if a.gas().is_some_and(|g| g != &gas))
            || matches!(&aperture, LoadedAperture::Modal(a) if a.gas().is_some_and(|g| g != &gas))
        {
            return Err(invalid(
                "DG supply mechanics require cold ExactZoh state and the actual combined inlet impedance, complete gas and clock",
            ));
        }
        Ok(Self {
            aperture,
            load,
            upstream,
            supply_channel: Some(channel),
        })
    }
}

#[cfg(test)]
mod tests;
