//! Spherical evaporation with the original finite gas and liquid calorics.
//!
//! Every interval derives radius/area from its single mass/enthalpy inventory
//! and the original liquid density, then uses the stagnant-sphere `D_AB/r`
//! coefficient or an explicitly supplied diameter-based Sherwood number with
//! the existing Stefan/kinetic phase owner. Accepted size is
//! reconstructed from final mass AND final material density. No radius state,
//! surrogate saturation source, constant-cp gas or latent heat source is added.
//!
//! Geometry and supplied molecular diffusivity are frozen within each interval;
//! the optional original heat-first split has the same first-order accuracy.
//! Refine the complete timestep, including film-source changes. The gas remains
//! its original fixed external control volume: droplet volume is not subtracted
//! from it and no sealed-vessel displacement or pressure-work model is implied.
//! Positive seeds, original inventory limits and all source/conservation gates
//! remain authoritative; complete dry-out, nucleation, deformation, droplet
//! motion, breakup, a forced-convection correlation and rarefied transport are
//! not supplied. Correlation selection/admission belongs to the transport owner.

/// Spherical inventories sharing one original finite-gas donor and caloric solve.
pub mod chamber;

use super::{
    GasEvaporationError, GasEvaporationSolve, GasEvaporationStep, GasFilmHeatSolve,
    GasFilmHeatStep, GasLiquidBinding, advance_surface_with_gas, advance_surface_with_gas_and_heat,
    poll,
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
pub use fs_thermal::evaporation::droplet::{SphericalDropletGeometry, SphericalDropletState};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, LiquidSurfaceState};

/// The actual geometric change and the original heat/phase energy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletReport {
    pub geometry: [SphericalDropletGeometry; 2],
    pub binary_diffusivity_m2_s: f64,
    /// Explicit diameter-based transfer ratio; 2 is the stagnant sphere.
    pub sherwood_number: f64,
    /// Initial Sh*D/diameter over this interval; refreshed on the next call.
    pub mass_transfer_coefficient_m_s: f64,
    /// Phase substep, after optional sensible heat.
    pub phase: GasEvaporationStep,
    /// Complete original-to-final exchange. Vapor energy already includes its
    /// latent subaccount; sensible heat is reported separately below.
    pub exchange: GasEvaporationStep,
    pub heat: Option<GasFilmHeatStep>,
}

/// One accepted sphere and finite gas. Neither original input is mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletCandidate<'catalog> {
    droplet: SphericalDropletState,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    report: GasDropletReport,
}

impl<'catalog> GasDropletCandidate<'catalog> {
    pub const fn droplet(&self) -> SphericalDropletState {
        self.droplet
    }
    pub const fn geometry(&self) -> &[SphericalDropletGeometry; 2] {
        &self.report.geometry
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &GasDropletReport {
        &self.report
    }
    /// Complete original-to-final storage and temperature changes.
    pub const fn step(&self) -> &GasEvaporationStep {
        &self.report.exchange
    }
    pub const fn phase_step(&self) -> &GasEvaporationStep {
        &self.report.phase
    }
    pub const fn heat(&self) -> Option<&GasFilmHeatStep> {
        self.report.heat.as_ref()
    }
    pub fn into_parts(
        self,
    ) -> (
        SphericalDropletState,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasDropletReport,
    ) {
        (self.droplet, self.gas, self.catalog, self.report)
    }
}

fn initial_surface<L: EvaporationLaw + ?Sized>(
    law: &L,
    droplet: SphericalDropletState,
    temperature: &impl Fn(f64) -> Result<f64, EvaporationError>,
    binary_diffusivity_m2_s: f64,
    sherwood_number: f64,
) -> Result<(LiquidSurfaceState, SphericalDropletGeometry, f64), GasEvaporationError> {
    let temperature_k = temperature(droplet.specific_enthalpy_j_kg)?;
    let geometry = droplet.geometry(law, temperature_k)?;
    let coefficient = geometry
        .mass_transfer_coefficient_at_sherwood_m_s(binary_diffusivity_m2_s, sherwood_number)?;
    Ok((
        LiquidSurfaceState {
            mass_kg: droplet.mass_kg,
            specific_enthalpy_j_kg: droplet.specific_enthalpy_j_kg,
            area_m2: geometry.area_m2(),
        },
        geometry,
        coefficient,
    ))
}

/// Advance one isolated spherical liquid inventory in the actual finite gas.
///
/// `binary_diffusivity_m2_s` is the caller's molecular transport source for the
/// interval, not thermal diffusivity or an inferred heat/mass analogy. Initial
/// gas temperature, pressure, humidity and donor inventory come from `gas`.
/// The shared phase law supplies the nonlinear Stefan film and selected finite
/// interface kinetics. The original caloric owner supplies carried vapor heat.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplet_with_gas<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    droplet: SphericalDropletState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
    binary_diffusivity_m2_s: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasDropletCandidate<'catalog>, GasEvaporationError> {
    advance_droplet_with_gas_at_sherwood(
        law,
        binding,
        droplet,
        gas,
        catalog,
        temperature,
        binary_diffusivity_m2_s,
        2.0,
        duration_s,
        controls,
        cancelled,
    )
}

/// Original finite-gas phase owner with explicit diameter-based transfer.
///
/// `k_m=Sh*D_AB/diameter`; molecular diffusivity and Sherwood retain separate
/// receipts. The caller admits the correlation, representative properties and
/// any dilute/blowing limits. No correlation is inferred here. `Sh=2` is the
/// exact original stagnant path. Area and coefficient freeze for this interval.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplet_with_gas_at_sherwood<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    droplet: SphericalDropletState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
    binary_diffusivity_m2_s: f64,
    sherwood_number: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasDropletCandidate<'catalog>, GasEvaporationError> {
    poll(&mut cancelled)?;
    let (surface, geometry, coefficient) =
        initial_surface(law, droplet, &temperature, binary_diffusivity_m2_s, sherwood_number)?;
    let candidate = advance_surface_with_gas(
        law,
        binding,
        surface,
        gas,
        catalog,
        temperature,
        coefficient,
        duration_s,
        controls,
        &mut cancelled,
    )?;
    let (liquid, gas, catalog, exchange) = candidate.into_parts();
    finish(
        law,
        liquid,
        gas,
        catalog,
        GasDropletReport {
            geometry: [geometry; 2],
            binary_diffusivity_m2_s,
            sherwood_number,
            mass_transfer_coefficient_m_s: coefficient,
            phase: exchange.clone(),
            exchange,
            heat: None,
        },
        &mut cancelled,
    )
}

/// Reciprocal sensible heat, phase transfer and accepted spherical geometry.
///
/// This uses the original heat-then-phase owner over one common interval. The
/// sensible coefficient excludes vapor-carried and latent energy. Initial area
/// and D/r stay fixed during this split; even with D=0, thermal density changes
/// update the accepted radius. No source or state is published on a refusal.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplet_with_gas_and_heat<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    droplet: SphericalDropletState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
    binary_diffusivity_m2_s: f64,
    heat_transfer_coefficient_w_m2_k: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    heat_controls: GasFilmHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasDropletCandidate<'catalog>, GasEvaporationError> {
    advance_droplet_with_gas_and_heat_at_sherwood(
        law,
        binding,
        droplet,
        gas,
        catalog,
        temperature,
        binary_diffusivity_m2_s,
        2.0,
        heat_transfer_coefficient_w_m2_k,
        duration_s,
        controls,
        heat_controls,
        cancelled,
    )
}

/// Original heat-first/phase split with explicit molecular D and Sherwood.
///
/// The sensible heat coefficient is independent, excludes vapor/latent energy,
/// and is never inferred from Sherwood. The same initially admitted area and
/// `Sh*D_AB/diameter` are frozen across both original substeps. Transport-model
/// validity remains caller authority; `Sh=2` preserves the original path.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplet_with_gas_and_heat_at_sherwood<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    droplet: SphericalDropletState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
    binary_diffusivity_m2_s: f64,
    sherwood_number: f64,
    heat_transfer_coefficient_w_m2_k: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    heat_controls: GasFilmHeatSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasDropletCandidate<'catalog>, GasEvaporationError> {
    poll(&mut cancelled)?;
    let (surface, geometry, coefficient) =
        initial_surface(law, droplet, &temperature, binary_diffusivity_m2_s, sherwood_number)?;
    let candidate = advance_surface_with_gas_and_heat(
        law,
        binding,
        surface,
        gas,
        catalog,
        temperature,
        coefficient,
        heat_transfer_coefficient_w_m2_k,
        duration_s,
        controls,
        heat_controls,
        &mut cancelled,
    )?;
    let phase = candidate.phase().step().clone();
    let (liquid, gas, catalog, exchange, heat) = candidate.into_parts();
    finish(
        law,
        liquid,
        gas,
        catalog,
        GasDropletReport {
            geometry: [geometry; 2],
            binary_diffusivity_m2_s,
            sherwood_number,
            mass_transfer_coefficient_m_s: coefficient,
            phase,
            exchange,
            heat: Some(heat),
        },
        &mut cancelled,
    )
}

fn finish<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    liquid: LiquidSurfaceState,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    mut report: GasDropletReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasDropletCandidate<'catalog>, GasEvaporationError> {
    poll(cancelled)?;
    let droplet = SphericalDropletState {
        mass_kg: liquid.mass_kg,
        specific_enthalpy_j_kg: liquid.specific_enthalpy_j_kg,
    };
    report.geometry[1] = droplet.geometry(law, report.exchange.surface_temperature_k[1])?;
    poll(cancelled)?;
    Ok(GasDropletCandidate {
        droplet,
        gas,
        catalog,
        report,
    })
}

#[cfg(test)]
mod tests;
