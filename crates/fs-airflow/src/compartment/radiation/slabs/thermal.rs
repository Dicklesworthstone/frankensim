//! Non-isothermal absorbing slabs heat the original finite thermal cells.
//!
//! Each optical layer maps to one existing material-backed cell. Its true
//! absorption-minus-emission power is added once to that cell; exterior slab
//! ports never own a separate surface heat source. The original conductive
//! links, latent storage, authored sources and event timeline remain active.
//!
//! The sheet adapter's single bounded backward-Euler/Newton driver owns all
//! numerical stepping. This adapter supplies only the endpoint optical source
//! resolver. Gray or genuine finite spectral bands retain the slab's matched-
//! index, plane-parallel, nonscattering model and diffuse enclosure-port closure.
//! Gas temperatures and optical geometry remain prescribed over the interval.

use std::collections::{BTreeMap, BTreeSet};

use fs_conduction::ConductionError;
use fs_conduction::radiation::SurfaceEmissivity;
use fs_math::dd::Dd;
use fs_qty::{Power, Temperature, Time};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};

use super::SlabPair;
use super::super::{
    LayerAbsorption, ParticipatingEnclosure, RadiationSolve, RadiationWorkspace,
    sheets::thermal::{self as shared},
    spectral::{SpectralOptics, SpectralSolve, SpectralWorkspace},
    walls::WallRadiation,
};

/// The original thermal driver's controls, reports and error semantics apply.
pub use shared::{
    MAX_SHEET_THERMAL_CELLS as MAX_SLAB_THERMAL_CELLS,
    SheetThermalCandidate as SlabThermalCandidate,
    SheetThermalError as SlabThermalError,
    SheetThermalReport as SlabThermalReport,
    SheetThermalSolve as SlabThermalSolve,
    SheetThermalSurface as SlabThermalSurface,
};

/// Optical source held over one original thermal candidate interval.
#[derive(Clone, Copy)]
pub enum SlabThermalOptics<'a, 's> {
    Gray {
        absorption: &'a [LayerAbsorption; 2],
        slabs: &'a [SlabPair],
        solve: RadiationSolve,
    },
    Spectral {
        optics: &'a SpectralOptics<'s>,
        solve: SpectralSolve,
    },
}
impl SlabThermalOptics<'_, '_> {
    fn slabs(&self) -> Result<&[SlabPair], SlabThermalError> {
        match self {
            Self::Gray { slabs, .. } => Ok(slabs),
            Self::Spectral { optics, .. } => optics.slab_bands()
                .and_then(|bands| bands.first()).map(Vec::as_slice)
                .ok_or(SlabThermalError::Invalid("explicit slab spectral bands are required")),
        }
    }
}

/// Reusable scratch, retaining no accepted material state.
pub struct SlabThermalWorkspace {
    thermal: shared::ThermalWorkspace,
    optical: OpticalWorkspace,
}
struct OpticalWorkspace {
    gray: RadiationWorkspace,
    spectral: SpectralWorkspace,
}

impl ParticipatingEnclosure {
    /// Prepare the shared bounded thermal Jacobian and original optical buffers.
    pub fn slab_thermal_workspace(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<SlabThermalWorkspace, SlabThermalError> {
        Ok(SlabThermalWorkspace {
            thermal: shared::ThermalWorkspace::new(thermal)?,
            optical: OpticalWorkspace {
                gray: self.workspace()?,
                spectral: self.spectral_workspace()?,
            },
        })
    }

    /// Advance ALL original thermal cells with layer-resolved radiative heating.
    ///
    /// Layer cells are supplied in slab order, then front-to-back layer order.
    /// Each must be active and distinct within its slab. The two optical ports
    /// must bind to the first/last layer's Stored cell, which may have different
    /// temperatures. Other faces retain ordinary Stored or Prescribed ownership.
    /// A cell shared by different optical stacks receives their summed power,
    /// with its original mass counted only once.
    ///
    /// The gray callback queries only non-port opaque faces at their actual
    /// endpoint temperatures. Spectral bands skip gray callbacks. Each spectral
    /// band's pairing and layer axis must agree, as admitted by SpectralOptics.
    /// No port G-J balance is deposited alongside the actual layer heat.
    ///
    /// Every probe starts from the same original history, including all native
    /// conductive links and sources. Process start-time events first; an interior
    /// event refuses, and accepted target-time events remain pending until all
    /// coupled owners accept. Gas/bath inputs must describe this same interval.
    /// Original local/total joule gates and the combined energy gate remain
    /// independent of timestep accuracy; refine physical time separately.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_slab_thermal<'m, 's>(
        &self,
        thermal: &DepositionSimulation<'m>,
        duration: Time,
        surfaces: &[SlabThermalSurface],
        layer_cells: &[Vec<CellId>],
        gas_temperatures: [Temperature; 2],
        optics: SlabThermalOptics<'_, 's>,
        mut emissivity: impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
        controls: SlabThermalSolve,
        work: &mut SlabThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SlabThermalCandidate<'m, 's>, SlabThermalError> {
        shared::checkpoint(&mut cancelled)?;
        if surfaces.len() != self.geometry().len()
            || gas_temperatures.iter().any(|t| !shared::positive(t.value()))
        {
            return Err(SlabThermalError::Invalid("slab thermal optical axes or gas temperature"));
        }
        let ports = admit_layers(thermal, surfaces, layer_cells, optics.slabs()?, &mut cancelled)?;
        let optical = &mut work.optical;
        shared::candidate_thermal(
            thermal, duration, controls, &mut work.thermal,
            |state, cancel| resolve(
                self, state, surfaces, layer_cells, &ports, gas_temperatures, optics,
                &mut emissivity, optical, cancel,
            ),
            cancelled,
        )
    }
}

fn admit_layers(
    thermal: &DepositionSimulation<'_>,
    surfaces: &[SlabThermalSurface],
    cells: &[Vec<CellId>],
    slabs: &[SlabPair],
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Vec<bool>, SlabThermalError> {
    if slabs.is_empty() || cells.len() != slabs.len() {
        return Err(SlabThermalError::Invalid("one original cell list per optical slab is required"));
    }
    for surface in surfaces {
        shared::checkpoint(cancelled)?;
        match *surface {
            SlabThermalSurface::Stored(id) if thermal.cells().contains_key(&id) => {}
            SlabThermalSurface::Prescribed(t) if shared::positive(t.value()) => {}
            _ => return Err(SlabThermalError::Invalid("inactive cell or invalid bath temperature")),
        }
    }
    let mut ports = vec![false; surfaces.len()];
    for (slab, ids) in slabs.iter().zip(cells) {
        shared::checkpoint(cancelled)?;
        let [front, back] = slab.surfaces;
        if front == back || front >= surfaces.len() || back >= surfaces.len()
            || ports[front] || ports[back]
            || ids.is_empty() || ids.len() != slab.optics.layer_count()
        {
            return Err(SlabThermalError::Invalid("invalid slab ports or thermal layer axis"));
        }
        let mut seen = BTreeSet::new();
        for &id in ids {
            shared::checkpoint(cancelled)?;
            if !thermal.cells().contains_key(&id) || !seen.insert(id) {
                return Err(SlabThermalError::Invalid(
                    "each slab layer needs a distinct active original thermal cell",
                ));
            }
        }
        if surfaces[front] != SlabThermalSurface::Stored(ids[0])
            || surfaces[back] != SlabThermalSurface::Stored(ids[ids.len() - 1])
        {
            return Err(SlabThermalError::Invalid(
                "slab ports must bind to the first and last original layer cells",
            ));
        }
        ports[front] = true;
        ports[back] = true;
    }
    Ok(ports)
}

#[allow(clippy::too_many_arguments)]
fn resolve<'s>(
    enclosure: &ParticipatingEnclosure,
    state: &DepositionSimulation<'_>,
    surfaces: &[SlabThermalSurface],
    layer_cells: &[Vec<CellId>],
    ports: &[bool],
    gas_temperatures: [Temperature; 2],
    optics: SlabThermalOptics<'_, 's>,
    emissivity: &mut impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
    work: &mut OpticalWorkspace,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<shared::Applied<'s>, SlabThermalError> {
    let mut temperatures = Vec::with_capacity(surfaces.len());
    for surface in surfaces {
        shared::checkpoint(cancelled)?;
        temperatures.push(match *surface {
            SlabThermalSurface::Stored(id) => Temperature::new(state.temperature_k(id)?),
            SlabThermalSurface::Prescribed(t) => t,
        });
    }
    let mut layer_temperatures = Vec::with_capacity(layer_cells.len());
    for ids in layer_cells {
        let mut values = Vec::with_capacity(ids.len());
        for &id in ids {
            shared::checkpoint(cancelled)?;
            values.push(Temperature::new(state.temperature_k(id)?));
        }
        layer_temperatures.push(values);
    }
    let mut cards = Vec::new();
    if matches!(optics, SlabThermalOptics::Gray { .. }) {
        for (surface, &temperature) in temperatures.iter().enumerate() {
            if ports[surface] { continue; }
            shared::checkpoint(cancelled)?;
            let card = emissivity(surface, temperature)
                .map_err(|error| SlabThermalError::Emissivity { surface, error })?;
            if card.temperature_k() != temperature.value() {
                return Err(SlabThermalError::Invalid(
                    "emissivity must query the actual endpoint temperature",
                ));
            }
            cards.push(card);
        }
    }
    let radiation = match optics {
        SlabThermalOptics::Gray { absorption, slabs, solve } => {
            WallRadiation::Gray(enclosure.evaluate_with_slabs(
                &temperatures, &cards, slabs, &layer_temperatures, gas_temperatures,
                absorption, solve, &mut work.gray, &mut *cancelled,
            )?)
        }
        SlabThermalOptics::Spectral { optics, solve } => {
            WallRadiation::Spectral(enclosure.evaluate_spectral_with_slabs(
                &temperatures, &cards, &layer_temperatures, gas_temperatures,
                optics, solve, &mut work.spectral, &mut *cancelled,
            )?)
        }
    };
    let layer_heat = match &radiation {
        WallRadiation::Gray(report) => &report.layer_heat_power,
        WallRadiation::Spectral(report) => &report.layer_heat_power,
    };
    if layer_heat.len() != layer_cells.len() {
        return Err(SlabThermalError::Invalid("returned slab heat has a different layer axis"));
    }
    let mut sums = BTreeMap::new();
    for (powers, ids) in layer_heat.iter().zip(layer_cells) {
        if powers.len() != ids.len() {
            return Err(SlabThermalError::Invalid("returned slab heat has a different layer axis"));
        }
        for (&power, &id) in powers.iter().zip(ids) {
            shared::checkpoint(cancelled)?;
            let value = sums.entry(id).or_insert(Dd::from_f64(0.0));
            *value = *value + Dd::from_f64(power.value());
        }
    }
    let mut bath_power = vec![Power::new(0.0); surfaces.len()];
    for (i, (&surface, &power)) in surfaces.iter()
        .zip(radiation.surface_heat_power()).enumerate()
    {
        shared::checkpoint(cancelled)?;
        if ports[i] {
            if power.value() != 0.0 {
                return Err(SlabThermalError::Invalid("slab ports cannot own surface heat"));
            }
            continue;
        }
        match surface {
            SlabThermalSurface::Stored(id) => {
                let value = sums.entry(id).or_insert(Dd::from_f64(0.0));
                *value = *value + Dd::from_f64(power.value());
            }
            SlabThermalSurface::Prescribed(_) => bath_power[i] = power,
        }
    }
    shared::pack_applied(radiation, temperatures, sums, bath_power)
}

#[cfg(test)]
mod tests;
