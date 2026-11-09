//! Joint spherical liquid inventories exchanging with one original finite gas.
//!
//! Every droplet derives area and diameter from its own positive mass, original
//! enthalpy chart and original liquid density. Its molecular diffusivity and
//! diameter-based Sherwood number supply k_m=Sh*D/diameter. The existing chamber
//! owns the single shared vapor donor, species update and gas caloric solve;
//! invoking independent scalar gas updates is never part of this adapter.
//!
//! Geometry, D and Sh freeze for one first-order interval, including optional
//! reciprocal sensible heat before phase transfer. Final geometry uses each
//! actual accepted mass AND final source density. There is no separate radius
//! inventory, statistical multiplicity, inferred sensible heat coefficient,
//! droplet interaction/shielding model, motion, breakup, coalescence, nucleation
//! or complete dry-out. The prescribed gas volume stays fixed and external;
//! droplet volume displacement and pressure work are not inferred.

/// Conservative finite-duration evolution with refreshed spherical geometry.
pub mod adaptive;

use super::super::chamber::{
    GasLiquidChamberCandidate, GasLiquidChamberSolve, GasLiquidChamberStep,
    GasLiquidChamberSurface, GasLiquidChamberSurfaceStep, Total, advance_chamber_with_gas,
    heat::{
        GasLiquidChamberHeatSolve, GasLiquidChamberHeatStep, GasLiquidChamberHeatSurface,
        GasLiquidChamberSurfaceHeatStep, advance_chamber_with_gas_and_heat,
    },
};
use super::{
    ClosedGasParcel, EvaporationError, EvaporationLaw, GasEvaporationError, GasFilmHeatSolve,
    GasLiquidBinding, SphericalDropletGeometry, SphericalDropletState, ThermochemicalGasModel,
    initial_surface, poll,
};

/// One physical droplet inventory and its original phase/caloric sources.
/// Stable IDs distinguish physical droplets; duplicate IDs refuse.
/// Each row can retain its own source law, volatile binding and enthalpy chart.
#[derive(Clone, Copy)]
pub struct GasDropletSurface<'a> {
    pub id: u64,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub droplet: SphericalDropletState,
    pub temperature_from_enthalpy: &'a dyn Fn(f64) -> Result<f64, EvaporationError>,
    pub binary_diffusivity_m2_s: f64,
    /// Explicit diameter-based transfer ratio; 2 preserves stagnant diffusion.
    /// Correlation selection, applicability and film-property sampling remain
    /// the caller's responsibility.
    pub sherwood_number: f64,
}

/// One physical sphere with a separately prescribed sensible heat law.
/// This coefficient excludes vapor-carried sensible energy and latent heat.
#[derive(Clone, Copy)]
pub struct GasDropletHeatSurface<'a> {
    pub surface: GasDropletSurface<'a>,
    pub heat_transfer_coefficient_w_m2_k: f64,
    pub heat_controls: GasFilmHeatSolve,
}

/// Actual endpoint geometry and the original chamber receipts for one droplet.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletSurfaceStep {
    pub id: u64,
    pub droplet: SphericalDropletState,
    pub geometry: [SphericalDropletGeometry; 2],
    pub binary_diffusivity_m2_s: f64,
    pub sherwood_number: f64,
    /// Initial Sh*D/diameter used throughout this interval.
    pub mass_transfer_coefficient_m_s: f64,
    /// Original phase-only account, after optional sensible heat.
    pub phase: GasLiquidChamberSurfaceStep,
    /// Original-to-final account; its liquid area is the FROZEN initial area.
    /// Use droplet plus geometry[1] for the actual accepted size.
    pub exchange: GasLiquidChamberSurfaceStep,
    pub heat: Option<GasLiquidChamberSurfaceHeatStep>,
}

/// One shared original donor/energy budget, never multiplied by droplet count.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletChamberReport {
    pub phase: GasLiquidChamberStep,
    pub exchange: GasLiquidChamberStep,
    pub heat: Option<GasLiquidChamberHeatStep>,
    /// Sum of absolute original liquid mass defects and the sole gas boundary
    /// defect, at least the original global and maximum-species defect magnitudes.
    /// Sensible heat has no mass transfer and does not duplicate this account.
    pub sum_absolute_mass_residual_kg: f64,
}

/// All positive droplet endpoints and their common finite gas publish together.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletChamberCandidate<'catalog> {
    droplets: Vec<GasDropletSurfaceStep>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    report: GasDropletChamberReport,
}
impl<'catalog> GasDropletChamberCandidate<'catalog> {
    /// Accepted droplets and original receipts in canonical stable-ID order.
    pub fn droplets(&self) -> &[GasDropletSurfaceStep] {
        &self.droplets
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &GasDropletChamberReport {
        &self.report
    }
    pub const fn step(&self) -> &GasLiquidChamberStep {
        &self.report.exchange
    }
    pub const fn phase_step(&self) -> &GasLiquidChamberStep {
        &self.report.phase
    }
    pub const fn heat(&self) -> Option<&GasLiquidChamberHeatStep> {
        self.report.heat.as_ref()
    }
    pub fn into_parts(
        self,
    ) -> (
        Vec<GasDropletSurfaceStep>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasDropletChamberReport,
    ) {
        (self.droplets, self.gas, self.catalog, self.report)
    }
}

/// Prepare all spherical films against the SAME original atmosphere and call
/// the original joint chamber owner once. Opposing streams retain their full
/// carried energy; gross condensation shares the old vapor inventory and
/// cannot borrow same-interval evaporation from another droplet.
///
/// Original kinetics, liquid/gas caloric references, endpoint equilibrium
/// checks, cancellation and shared work/conservation budgets remain active.
/// Empty rosters preserve the original chamber's empty-interval behavior.
/// Inputs stay unchanged if any droplet or the shared gas refuses.
pub fn advance_droplets_with_gas<'catalog>(
    surfaces: &[GasDropletSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasDropletChamberCandidate<'catalog>, GasEvaporationError> {
    advance(
        surfaces.iter().copied().map(|surface| (surface, None)),
        gas,
        catalog,
        duration_s,
        controls,
        None,
        cancelled,
    )
}

/// One simultaneous reciprocal sensible-heat solve followed by one joint phase
/// interval. Every droplet sees the same evolving gas heat capacity, then the
/// same heated gas composition/temperature in the original phase owner.
/// Initial sphere area and Sh*D/diameter are retained through both stages;
/// endpoint size includes thermal density changes even when D=0.
pub fn advance_droplets_with_gas_and_heat<'catalog>(
    surfaces: &[GasDropletHeatSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasDropletChamberCandidate<'catalog>, GasEvaporationError> {
    advance(
        surfaces.iter().map(|surface| {
            (
                surface.surface,
                Some((
                    surface.heat_transfer_coefficient_w_m2_k,
                    surface.heat_controls,
                )),
            )
        }),
        gas,
        catalog,
        duration_s,
        controls,
        Some(heat_controls),
        cancelled,
    )
}

struct Prepared<'a> {
    input: GasDropletSurface<'a>,
    geometry: SphericalDropletGeometry,
    phase: GasLiquidChamberSurface<'a>,
    heat: Option<(f64, GasFilmHeatSolve)>,
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<'surface, 'catalog>(
    surfaces: impl ExactSizeIterator<
        Item = (GasDropletSurface<'surface>, Option<(f64, GasFilmHeatSolve)>),
    >,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    heat_controls: Option<GasLiquidChamberHeatSolve>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasDropletChamberCandidate<'catalog>, GasEvaporationError> {
    poll(&mut cancelled)?;
    if controls.maximum_surfaces == 0 || surfaces.len() > controls.maximum_surfaces {
        return Err(GasEvaporationError::Invalid(
            "spherical droplet count exceeds its positive chamber work ceiling",
        ));
    }
    let mut ordered: Vec<_> = surfaces.collect();
    ordered.sort_unstable_by_key(|(surface, _)| surface.id);
    if ordered.windows(2).any(|pair| pair[0].0.id == pair[1].0.id) {
        return Err(GasEvaporationError::Invalid(
            "spherical droplets require unique stable IDs",
        ));
    }
    let mut prepared = Vec::with_capacity(ordered.len());
    for (input, heat) in ordered {
        poll(&mut cancelled)?;
        let (liquid, geometry, coefficient) = initial_surface(
            input.law,
            input.droplet,
            &input.temperature_from_enthalpy,
            input.binary_diffusivity_m2_s,
            input.sherwood_number,
        )?;
        prepared.push(Prepared {
            input,
            geometry,
            phase: GasLiquidChamberSurface {
                id: input.id,
                law: input.law,
                binding: input.binding,
                liquid,
                temperature_from_enthalpy: input.temperature_from_enthalpy,
                mass_transfer_coefficient_m_s: coefficient,
            },
            heat,
        });
    }
    let (phase_surfaces, exchanges, gas, catalog, report, heat_surfaces) =
        if let Some(heat_controls) = heat_controls {
            let original: Vec<_> = prepared
                .iter()
                .map(|row| {
                    let (coefficient, controls) =
                        row.heat.expect("heat entry point supplies heat rows");
                    GasLiquidChamberHeatSurface {
                        surface: row.phase,
                        heat_transfer_coefficient_w_m2_k: coefficient,
                        heat_controls: controls,
                    }
                })
                .collect();
            let next = advance_chamber_with_gas_and_heat(
                &original,
                gas,
                catalog,
                duration_s,
                controls,
                heat_controls,
                &mut cancelled,
            )?;
            let sum_absolute_mass_residual_kg = mass_defect(next.phase(), &mut cancelled)?;
            let phase = next.phase().step().clone();
            let phase_surfaces = next.phase().surfaces().to_vec();
            let (exchanges, gas, catalog, exchange, heat, heat_surfaces) = next.into_parts();
            (
                phase_surfaces,
                exchanges,
                gas,
                catalog,
                GasDropletChamberReport {
                    phase,
                    exchange,
                    heat: Some(heat),
                    sum_absolute_mass_residual_kg,
                },
                heat_surfaces,
            )
        } else {
            let original: Vec<_> = prepared.iter().map(|row| row.phase).collect();
            let next = advance_chamber_with_gas(
                &original,
                gas,
                catalog,
                duration_s,
                controls,
                &mut cancelled,
            )?;
            let sum_absolute_mass_residual_kg = mass_defect(&next, &mut cancelled)?;
            let phase = next.step().clone();
            let (exchanges, gas, catalog, exchange) = next.into_parts();
            (
                exchanges.clone(),
                exchanges,
                gas,
                catalog,
                GasDropletChamberReport {
                    phase,
                    exchange,
                    heat: None,
                    sum_absolute_mass_residual_kg,
                },
                Vec::new(),
            )
        };
    if phase_surfaces.len() != prepared.len()
        || exchanges.len() != prepared.len()
        || (report.heat.is_some() && heat_surfaces.len() != prepared.len())
    {
        return Err(GasEvaporationError::Invalid(
            "original spherical chamber omitted a droplet receipt",
        ));
    }
    let mut heat_surfaces = heat_surfaces.into_iter();
    let mut droplets = Vec::with_capacity(prepared.len());
    for ((row, phase), exchange) in prepared.into_iter().zip(phase_surfaces).zip(exchanges) {
        poll(&mut cancelled)?;
        let heat = heat_surfaces.next();
        if phase.id != row.input.id
            || exchange.id != row.input.id
            || heat.as_ref().is_some_and(|heat| heat.id != row.input.id)
        {
            return Err(GasEvaporationError::Invalid(
                "original spherical chamber changed a droplet identity",
            ));
        }
        let droplet = SphericalDropletState {
            mass_kg: exchange.liquid.mass_kg,
            specific_enthalpy_j_kg: exchange.liquid.specific_enthalpy_j_kg,
        };
        let final_geometry = droplet.geometry(row.input.law, exchange.surface_temperature_k[1])?;
        droplets.push(GasDropletSurfaceStep {
            id: row.input.id,
            droplet,
            geometry: [row.geometry, final_geometry],
            binary_diffusivity_m2_s: row.input.binary_diffusivity_m2_s,
            sherwood_number: row.input.sherwood_number,
            mass_transfer_coefficient_m_s: row.phase.mass_transfer_coefficient_m_s,
            phase,
            exchange,
            heat,
        });
    }
    poll(&mut cancelled)?;
    Ok(GasDropletChamberCandidate {
        droplets,
        gas,
        catalog,
        report,
    })
}

fn mass_defect(
    phase: &GasLiquidChamberCandidate<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasEvaporationError> {
    let mut absolute = Total::default();
    absolute.add(phase.gas_exchange().step().mass_residual_kg.abs())?;
    for surface in phase.surfaces() {
        poll(cancelled)?;
        absolute.add(surface.mass_residual_kg.abs())?;
    }
    Ok(absolute
        .value()?
        .max(phase.step().mass_residual_kg.abs())
        .max(phase.step().maximum_species_mass_residual_kg.abs()))
}

#[cfg(test)]
mod tests;
