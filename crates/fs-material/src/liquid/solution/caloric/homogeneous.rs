//! Source-resolved homogeneous liquid feeds, without assuming phase equilibrium.
//!
//! An initially homogeneous mixture can be metastable or unstable. Its original
//! excess enthalpy must remain in the heat account until separation occurs:
//! replacing it by an isothermal equilibrium before a thermal solve discards
//! the heat of demixing. The existing L1 caloric owner supplies this initial
//! state and its subsequent energy-balanced equilibrium. No additional heater,
//! Cp integration or phase solver is introduced here.
//!
//! A distinct feed type retains its actual material properties and component
//! inventory without claiming equilibrium or entering equilibrium-only tank
//! consumers. Relaxation returns the established equilibrium/balance types.
//! This is an equilibrium endpoint calculation at maintained pressure, not a
//! separation rate, nucleation model, spinodal dynamics or interfacial-energy
//! account. Transport Cp retains the same independent source qualification as
//! the parent caloric adapter.

use super::super::{BinaryLiquidError, positive_value};
use super::{
    BinaryCaloricErrorV1, BinaryCaloricSolveV1, BinaryCaloricStateV1, BinaryLiquidCaloricBalance,
    BinaryLiquidCaloricError, BinaryLiquidCaloricLaw, SourceOwner, poll,
};
use crate::MaterialError;
use crate::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use crate::liquid::ResolvedLiquidState;
use fs_qty::Temperature;
use fs_thermochem::solution::BinarySolutionStateV1;
use std::sync::Arc;

/// One explicitly homogeneous finite feed at its original temperature and
/// composition. Local thermodynamic instability does not silently equilibrate
/// it. The source owner and original enthalpy remain available for relaxation.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedHomogeneousBinaryLiquidFeed {
    owner: Arc<SourceOwner>,
    caloric: BinaryCaloricStateV1,
    homogeneous: BinarySolutionStateV1,
    material: ResolvedLiquidState,
    mass_kg: f64,
    volume_m3: f64,
}

impl ResolvedHomogeneousBinaryLiquidFeed {
    /// Original homogeneous enthalpy. `equilibrium()` on this caloric state is None.
    pub const fn caloric_state(&self) -> &BinaryCaloricStateV1 {
        &self.caloric
    }
    /// Original free-energy response, including the actual thermodynamic factor.
    /// A negative factor identifies a spinodal instability, not negative diffusion.
    pub const fn homogeneous_state(&self) -> &BinarySolutionStateV1 {
        &self.homogeneous
    }
    /// Original transport claims at the actual initial temperature/composition.
    pub const fn material(&self) -> &ResolvedLiquidState {
        &self.material
    }
    pub const fn temperature(&self) -> Temperature {
        self.caloric.temperature()
    }
    /// Conserved amounts [mol], in the original pair's species order.
    pub const fn component_moles(&self) -> [f64; 2] {
        self.caloric.component_moles()
    }
    pub const fn mole_fractions(&self) -> [f64; 2] {
        self.homogeneous.mole_fractions()
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Source-derived volume, not an additive pure-component volume estimate.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
}

impl FluidTransport for ResolvedHomogeneousBinaryLiquidFeed {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.material.transport_state()
    }
}
impl ThermalFluidTransport for ResolvedHomogeneousBinaryLiquidFeed {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        ThermalFluidTransportState::try_new(
            self.transport_state()?,
            self.temperature().value(),
            self.material.volumetric_expansion_per_k(),
        )
    }
}

impl BinaryLiquidCaloricLaw<'_, '_> {
    /// Resolve an explicitly homogeneous feed, including metastable or unstable
    /// compositions. Both caloric and transport sources must admit the original
    /// state. No phase equilibrium, heat exchange or inventory change occurs.
    /// Refusal or cancellation publishes no partially resolved feed.
    pub fn homogeneous(
        &self,
        temperature: Temperature,
        component_moles: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedHomogeneousBinaryLiquidFeed, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        let caloric = self.owner.model.homogeneous(temperature, component_moles)?;
        let total_moles = positive_value(component_moles[0] + component_moles[1], "total amount")?;
        // Use the same original pair owner and its two-fraction normalization
        // as L1 homogeneous initialization. This is a local constitutive sample,
        // not an equilibrium solve or a competing caloric representation.
        let homogeneous = self
            .owner
            .model
            .thermodynamics()
            .evaluate_composition(temperature, component_moles.map(|n| n / total_moles))
            .map_err(BinaryCaloricErrorV1::from)?;
        let fractions = homogeneous.mole_fractions();
        if (component_moles[0] > 0.0 && fractions[0] == 0.0)
            || (component_moles[1] > 0.0 && fractions[0] == 1.0)
        {
            return Err(BinaryLiquidError::Arithmetic(
                "nonzero component lost in the source's first-species coordinate",
            )
            .into());
        }
        poll(&mut cancelled)?;
        let point = self
            .liquid
            .point
            .clone()
            .with(&self.liquid.temperature_axis, temperature.value())
            .and_then(|point| point.with(&self.liquid.composition_axis, fractions[0]))
            .map_err(|_| {
                BinaryLiquidError::Invalid(
                    "homogeneous feed violates source coordinate conventions",
                )
            })?;
        let material = self
            .liquid
            .reference
            .at_query_point(self.liquid.card, &point)
            .map_err(BinaryLiquidError::from)?;
        let component_mass = self.liquid.component_masses(component_moles)?;
        let mass_kg = positive_value(
            component_mass[0] + component_mass[1],
            "homogeneous feed mass",
        )?;
        let volume_m3 = positive_value(
            mass_kg / material.density_kg_m3(),
            "homogeneous feed volume",
        )?;
        poll(&mut cancelled)?;
        Ok(ResolvedHomogeneousBinaryLiquidFeed {
            owner: Arc::clone(&self.owner),
            caloric,
            homogeneous,
            material,
            mass_kg,
            volume_m3,
        })
    }

    /// Relax the original homogeneous feed with signed external heat [J].
    /// Zero heat resolves adiabatic demixing: released excess enthalpy changes
    /// the temperature and therefore the accepted phase compositions. Positive
    /// heat can instead carry the feed through its miscibility boundary.
    ///
    /// The existing energy solver retains the initial homogeneous enthalpy and
    /// both component inventories. Only its accepted equilibrium refreshes the
    /// original mixture sources, and the returned balance remains compatible
    /// with ordinary `add_heat`, `mix` and equilibrium material consumers.
    /// Inputs stay unchanged on failed convergence, source loss or cancellation.
    pub fn relax_homogeneous(
        &self,
        initial: &ResolvedHomogeneousBinaryLiquidFeed,
        heat_input_j: f64,
        controls: BinaryCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidCaloricBalance, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.owner, &initial.owner) && self.owner != initial.owner {
            return Err(BinaryLiquidCaloricError::Invalid(
                "homogeneous feed requires its original caloric and mixture source binding",
            ));
        }
        let balance =
            self.owner
                .model
                .add_heat(&initial.caloric, heat_input_j, controls, &mut cancelled)?;
        self.resolve_balance(balance, &mut cancelled)
    }
}
