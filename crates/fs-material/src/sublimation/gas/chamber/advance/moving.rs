//! Moving phase area in the ORIGINAL adaptive nonisothermal chamber.
//!
//! Geometry is a function of conserved mass, not another integrated field.
//! Every coarse and fine half step refreshes it; both sticking and reflected
//! heat use that step's actual starting area. The original gas/solid caloric
//! owners, shared inventory allowance and error controller remain unchanged.
//!
//! This is explicitly a dilute-body approximation: gas volume is held fixed,
//! and total solid volume must remain below a declared fraction of it, at most
//! 1%. This bounds geometric displacement relative to gas volume, NOT omitted
//! pressure-work energy or total physical error. No gas-volume work, surface
//! free energy, Kelvin pressure, mesh/contact update or new phase is inferred.
//! Original laws must describe constant-density, uniform, shape-preserving
//! bodies. Complete finite-energy extinction is still refused by the molecular
//! kernel; the separate thermostatted front has an explicit extinction event.

use super::{
    GasChamberAdvanceCandidate, GasChamberAdvanceError, GasChamberAdvanceSolve,
    GasChamberCandidate, GasChamberSurface, GasSurfaceError, ClosedGasParcel,
    ThermochemicalGasModel, ThermalAccommodation, advance_history, error_ratio,
};
use fs_thermal::phase_geometry::{FrontGeometry, RetreatGeometry};
use std::collections::BTreeMap;

/// Explicit geometry and approximation/error admission, borrowed across trials.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingChamberGeometry {
    pub surfaces: BTreeMap<u64, RetreatGeometry>,
    /// Sum of every selected solid volume / fixed gas volume, in (0,0.01].
    /// This is an approximation-domain choice, not measured validation.
    pub maximum_solid_volume_fraction: f64,
    /// Absolute radius/thickness error scale [m], added to the same relative
    /// tolerance used for mass and temperature in the original controller.
    pub absolute_extent_tolerance_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingSurfaceHistory {
    pub initial: FrontGeometry,
    pub final_geometry: FrontGeometry,
    /// Exact mass-based shell displacement, not dm divided by a frozen area.
    pub normal_recession_m: f64,
}

/// Existing complete conservative history plus mass-derived endpoint geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingChamberCandidate<'catalog> {
    history: GasChamberAdvanceCandidate<'catalog>,
    geometry: BTreeMap<u64, MovingSurfaceHistory>,
}
impl<'catalog> MovingChamberCandidate<'catalog> {
    #[must_use]
    pub const fn history(&self) -> &GasChamberAdvanceCandidate<'catalog> { &self.history }
    #[must_use]
    pub fn geometry(&self) -> &BTreeMap<u64, MovingSurfaceHistory> { &self.geometry }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel { self.history.gas() }
    #[must_use]
    pub fn into_parts(self) -> (GasChamberAdvanceCandidate<'catalog>, BTreeMap<u64, MovingSurfaceHistory>) {
        (self.history, self.geometry)
    }
}

impl MovingChamberGeometry {
    pub(super) fn validate(&self, surfaces: &[GasChamberSurface<'_>]) -> Result<(), GasSurfaceError> {
        if !self.maximum_solid_volume_fraction.is_finite()
            || self.maximum_solid_volume_fraction <= 0.0 || self.maximum_solid_volume_fraction > 0.01
            || !self.absolute_extent_tolerance_m.is_finite() || self.absolute_extent_tolerance_m <= 0.0
            || self.surfaces.len() != surfaces.len()
            || surfaces.iter().any(|surface| !self.surfaces.contains_key(&surface.id)) {
            return Err(GasSurfaceError::Invalid("complete moving geometry, positive extent tolerance and dilute volume cap required"));
        }
        Ok(())
    }
    fn volume_admission(&self, volume: f64, gas: &ClosedGasParcel) -> Result<(), GasSurfaceError> {
        let fraction = volume / gas.volume_m3();
        if !fraction.is_finite() || fraction < 0.0 || fraction > self.maximum_solid_volume_fraction {
            Err(GasSurfaceError::Invalid("phase volume exceeded declared fixed-gas-volume approximation"))
        } else { Ok(()) }
    }
    pub(super) fn admit_interval(&self, surfaces: &[GasChamberSurface<'_>], gas: &ClosedGasParcel,
        cancelled: &mut impl FnMut() -> bool) -> Result<(), GasSurfaceError> {
        let mut volumes = super::super::Total::default();
        for surface in surfaces {
            super::super::super::poll(cancelled)?;
            let geometry = self.surfaces.get(&surface.id)
                .ok_or(GasSurfaceError::Invalid("missing moving phase geometry"))?;
            let physical = geometry.at_mass(surface.solid.mass_kg)?;
            let area = surface.solid.area_m2;
            if !area.is_finite() || area <= 0.0
                || (area / physical.area_m2 - 1.0).abs() > 64.0 * f64::EPSILON {
                return Err(GasSurfaceError::Invalid("surface area disagrees with its mass-owned moving geometry"));
            }
            let t = (surface.temperature_from_enthalpy)(surface.solid.specific_enthalpy_j_kg)?;
            geometry.admit_density(surface.law.properties(t)?.solid_density_kg_m3)?;
            volumes.add(physical.volume_m3)?;
        }
        self.volume_admission(volumes.value()?, gas)
    }
    pub(super) fn refresh_endpoint(&self, originals: &[GasChamberSurface<'_>],
        endpoint: &mut GasChamberCandidate<'_>, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasSurfaceError> {
        if originals.len() != endpoint.surfaces.len() {
            return Err(GasSurfaceError::Invalid("moving phase endpoint roster changed"));
        }
        let mut volumes = super::super::Total::default();
        for (original, row) in originals.iter().zip(&mut endpoint.surfaces) {
            super::super::super::poll(cancelled)?;
            if original.id != row.id { return Err(GasSurfaceError::Invalid("moving phase identity changed")); }
            let geometry = self.surfaces.get(&row.id)
                .ok_or(GasSurfaceError::Invalid("missing endpoint phase geometry"))?;
            geometry.admit_density(original.law.properties(row.surface_temperature_k[1])?.solid_density_kg_m3)?;
            let physical = geometry.at_mass(row.solid.mass_kg)?;
            if physical.mass_kg == 0.0 {
                return Err(GasSurfaceError::Invalid("finite-energy phase extinction requires its own event owner"));
            }
            row.solid.area_m2 = physical.area_m2;
            row.recession_m = geometry.recession_m(original.solid.mass_kg, row.solid.mass_kg)?;
            volumes.add(physical.volume_m3)?;
        }
        self.volume_admission(volumes.value()?, endpoint.gas())
    }
    pub(super) fn compare(&self, originals: &[GasChamberSurface<'_>], coarse: &GasChamberCandidate<'_>,
        fine: &GasChamberCandidate<'_>, relative: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, GasChamberAdvanceError> {
        let mut maximum: f64 = 0.0;
        for ((initial, a), b) in originals.iter().zip(coarse.surfaces()).zip(fine.surfaces()) {
            super::poll(cancelled)?;
            let law = &self.surfaces[&initial.id];
            let extent = |mass| law.at_mass(mass).map(|state| state.extent_m).map_err(GasSurfaceError::from);
            maximum = maximum.max(error_ratio(extent(a.solid.mass_kg)?, extent(b.solid.mass_kg)?,
                extent(initial.solid.mass_kg)?, self.absolute_extent_tolerance_m, relative)?);
        }
        Ok(maximum)
    }
}

/// Adaptive finite-energy sublimation/frost growth with actual changing area.
///
/// A caller supplies mass-consistent initial areas and exactly one geometry for
/// each surface. The original controller compares complete moving-geometry
/// trials, additionally controlling physical radius/thickness error. Geometry
/// and optional reflection tables are borrowed, never cloned on each half step.
/// Returned solid areas can be passed directly into the next continuation call.
/// No accepted gas/solid/geometry is published on a late failure or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_adaptive_moving<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    geometry: &MovingChamberGeometry,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasChamberAdvanceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MovingChamberCandidate<'catalog>, GasChamberAdvanceError> {
    let history = advance_history(surfaces, accommodation, Some(geometry), gas, catalog,
        duration_s, controls, &mut cancelled)?;
    let originals: BTreeMap<_, _> = surfaces.iter().map(|surface| (surface.id, surface.solid.mass_kg)).collect();
    let mut reports = BTreeMap::new();
    for endpoint in history.endpoint().surfaces() {
        super::poll(&mut cancelled)?;
        let shape = geometry.surfaces[&endpoint.id];
        let mass = originals[&endpoint.id];
        reports.insert(endpoint.id, MovingSurfaceHistory {
            initial: shape.at_mass(mass).map_err(GasSurfaceError::from)?,
            final_geometry: shape.at_mass(endpoint.solid.mass_kg).map_err(GasSurfaceError::from)?,
            normal_recession_m: shape.recession_m(mass, endpoint.solid.mass_kg).map_err(GasSurfaceError::from)?,
        });
    }
    Ok(MovingChamberCandidate { history, geometry: reports })
}

#[cfg(test)]
mod tests;
