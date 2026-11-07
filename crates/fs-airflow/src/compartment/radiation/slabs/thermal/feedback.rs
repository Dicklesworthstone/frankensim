//! State-dependent optics in the original coupled thermal energy residual.

use fs_conduction::radiation::slab::{AbsorbingLayer, AbsorbingSlab};
use fs_qty::Length;
use fs_thermal::enthalpy::DepositedCell;

use super::*;

/// Actual native endpoint state for one original source-bearing optical layer.
///
/// The reference layer always comes from the caller's original optics, never
/// from a previous Newton probe. The native cell retains specific enthalpy,
/// mass and material identity, so a phase-dependent law can change absorption
/// on a latent plateau without inventing a temperature-dependent heat capacity.
#[derive(Clone, Copy, Debug)]
pub struct SlabLayerState<'a> {
    /// `None` for genuine total-gray optics; otherwise the original band index.
    pub band: Option<usize>,
    /// Original finite band endpoints; never an inferred representative wavelength.
    pub wavelengths: Option<[Length; 2]>,
    pub slab: usize,
    pub layer: usize,
    pub cell: &'a DepositedCell,
    pub temperature: Temperature,
    pub reference: &'a AbsorbingLayer,
}

impl ParticipatingEnclosure {
    /// Couple temperature/enthalpy-dependent layer optics to native storage.
    ///
    /// Every backward-Euler residual, Jacobian and globalization trial resolves
    /// `layer_optics` at its actual material state. Return the sourced gray or
    /// band-constant absorption and real refractive index for THAT state. A law
    /// must be deterministic and must not advance its own accepted history.
    /// The input thermal state and original optics remain unchanged on refusal.
    ///
    /// Thickness, layer order, cell identity and port geometry stay fixed.
    /// Changing those requires a geometry/remap owner. Changed optical layers
    /// rebuild the ORIGINAL angular/Fresnel transport; identical responses
    /// reuse its immutable coefficients. Each preparation retains the original
    /// slab's quadrature controls, under the bounded shared thermal trial count.
    /// Fresnel reflection, total internal reflection, n² internal emission and
    /// all repeated returns therefore respond to the same endpoint state.
    ///
    /// Gray opaque-face emissivity remains its original endpoint callback;
    /// spectral opaque-face/gas/soot inputs and band coverage stay unchanged.
    /// This is an implicit layer-optics update, not a gas, deformation, optical
    /// coherence or nonequilibrium phase-kinetics solver. The returned candidate
    /// already includes its reported layer heat; accept gas/bath owners once.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_slab_thermal_with_feedback<'m, 's>(
        &self,
        thermal: &DepositionSimulation<'m>,
        duration: Time,
        surfaces: &[SlabThermalSurface],
        layer_cells: &[Vec<CellId>],
        gas_temperatures: [Temperature; 2],
        optics: SlabThermalOptics<'_, 's>,
        mut layer_optics: impl FnMut(SlabLayerState<'_>) -> Result<AbsorbingLayer, ConductionError>,
        mut emissivity: impl FnMut(usize, Temperature) -> Result<SurfaceEmissivity, ConductionError>,
        controls: SlabThermalSolve,
        work: &mut SlabThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SlabThermalCandidate<'m, 's>, SlabThermalError> {
        shared::checkpoint(&mut cancelled)?;
        if surfaces.len() != self.geometry().len()
            || gas_temperatures
                .iter()
                .any(|t| !shared::positive(t.value()))
        {
            return Err(SlabThermalError::Invalid(
                "slab thermal optical axes or gas temperature",
            ));
        }
        let ports = admit_layers(
            thermal,
            surfaces,
            layer_cells,
            optics.slabs()?,
            &mut cancelled,
        )?;
        let mut current = FeedbackOptics::new(optics);
        let optical = &mut work.optical;
        shared::candidate_thermal(
            thermal,
            duration,
            controls,
            &mut work.thermal,
            |state, cancel| {
                current.refresh(optics, state, layer_cells, &mut layer_optics, cancel)?;
                resolve(
                    self,
                    state,
                    surfaces,
                    layer_cells,
                    &ports,
                    gas_temperatures,
                    current.as_optics(optics),
                    &mut emissivity,
                    optical,
                    cancel,
                )
            },
            cancelled,
        )
    }
}

enum FeedbackOptics<'s> {
    Gray(Vec<SlabPair>),
    Spectral(SpectralOptics<'s>),
}
impl<'s> FeedbackOptics<'s> {
    fn new(optics: SlabThermalOptics<'_, 's>) -> Self {
        match optics {
            SlabThermalOptics::Gray { slabs, .. } => Self::Gray(slabs.to_vec()),
            SlabThermalOptics::Spectral { optics, .. } => Self::Spectral(optics.clone()),
        }
    }

    fn as_optics<'a>(&'a self, original: SlabThermalOptics<'a, 's>) -> SlabThermalOptics<'a, 's> {
        match (self, original) {
            (
                Self::Gray(slabs),
                SlabThermalOptics::Gray {
                    absorption, solve, ..
                },
            ) => SlabThermalOptics::Gray {
                absorption,
                slabs,
                solve,
            },
            (Self::Spectral(optics), SlabThermalOptics::Spectral { solve, .. }) => {
                SlabThermalOptics::Spectral { optics, solve }
            }
            _ => unreachable!("feedback optics retain the admitted original kind"),
        }
    }

    fn refresh(
        &mut self,
        original: SlabThermalOptics<'_, 's>,
        state: &DepositionSimulation<'_>,
        cells: &[Vec<CellId>],
        law: &mut impl FnMut(SlabLayerState<'_>) -> Result<AbsorbingLayer, ConductionError>,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), SlabThermalError> {
        match (self, original) {
            (Self::Gray(current), SlabThermalOptics::Gray { slabs, .. }) => {
                refresh_pairs(current, slabs, None, state, cells, law, cancelled)?;
            }
            (Self::Spectral(current), SlabThermalOptics::Spectral { optics, .. }) => {
                let mut bands = current
                    .slab_bands()
                    .ok_or(SlabThermalError::Invalid(
                        "explicit slab spectral bands are required",
                    ))?
                    .to_vec();
                let mut changed = false;
                for (band, ((pairs, reference), source)) in bands
                    .iter_mut()
                    .zip(optics.slab_bands().ok_or(SlabThermalError::Invalid(
                        "explicit slab spectral bands are required",
                    ))?)
                    .zip(optics.bands())
                    .enumerate()
                {
                    changed |= refresh_pairs(
                        pairs,
                        reference,
                        Some((band, source.wavelengths)),
                        state,
                        cells,
                        law,
                        cancelled,
                    )?;
                }
                if changed {
                    *current = optics.with_slab_feedback(bands)?;
                }
            }
            _ => unreachable!("feedback optics retain the admitted original kind"),
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn refresh_pairs(
    current: &mut [SlabPair],
    reference: &[SlabPair],
    band: Option<(usize, [Length; 2])>,
    state: &DepositionSimulation<'_>,
    cells: &[Vec<CellId>],
    law: &mut impl FnMut(SlabLayerState<'_>) -> Result<AbsorbingLayer, ConductionError>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<bool, SlabThermalError> {
    let mut changed = false;
    for (slab, ((current, original), ids)) in
        current.iter_mut().zip(reference).zip(cells).enumerate()
    {
        let mut layers = Vec::with_capacity(ids.len());
        for (layer, (&id, reference)) in ids.iter().zip(original.optics.layers()).enumerate() {
            shared::checkpoint(cancelled)?;
            let cell = state
                .cells()
                .get(&id)
                .ok_or(SlabThermalError::Invalid("inactive slab cell"))?;
            let temperature = Temperature::new(state.temperature_k(id)?);
            let updated = law(SlabLayerState {
                band: band.map(|b| b.0),
                wavelengths: band.map(|b| b.1),
                slab,
                layer,
                cell,
                temperature,
                reference,
            })
            .map_err(|error| match error {
                ConductionError::Cancelled { .. } => SlabThermalError::Cancelled,
                error => SlabThermalError::LayerOptics {
                    band: band.map(|b| b.0),
                    slab,
                    layer: Some(layer),
                    error,
                },
            })?;
            let support = updated.temperature_support_k();
            if updated.thickness_m() != reference.thickness_m()
                || !(support[0]..=support[1]).contains(&temperature.value())
            {
                return Err(SlabThermalError::Invalid(
                    "layer feedback must retain thickness and cover its actual endpoint temperature",
                ));
            }
            layers.push(updated);
        }
        if layers != current.optics.layers() {
            current.optics = AbsorbingSlab::new(
                layers,
                original.optics.quadrature().controls,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                ConductionError::Cancelled { .. } => SlabThermalError::Cancelled,
                error => SlabThermalError::LayerOptics {
                    band: band.map(|b| b.0),
                    slab,
                    layer: None,
                    error,
                },
            })?;
            changed = true;
        }
    }
    Ok(changed)
}
