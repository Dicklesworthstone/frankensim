//! Implicit density feedback in a fixed-footprint planar liquid column.
//!
//! Fixed material inventories occupy ordered, unmixed liquid layers. Their
//! original rho(h) determines current thickness at the explicitly fixed area;
//! the first face stays at coordinate zero and the far face is free to move.
//! The source must apply under the caller's maintained pressure/composition.
//! This is a quasistatic one-dimensional volume closure, not a constrained
//! solid expansion, momentum, mixing, sidewall, or pressure-displacement solve.
//!
//! Native material-density layers supply the same current half-thickness and
//! exact d(1/rho)/dh to Fourier/contact/bath fluxes and their Newton derivatives.
//! Storage remains the original m*h; no separate pressure-work heat is added.

use super::*;
use fs_thermal::enthalpy::conduction::{
    PlanarDensityLayer,
    implicit::{ImplicitConductionObservables, ImplicitConductionState},
};

/// Current physical liquid interval in the column's normal coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarLiquidCellGeometry {
    pub cell_id: CellId,
    pub mass_kg: f64,
    pub density_kg_m3: f64,
    pub volume_m3: f64,
    pub interval_m: [f64; 2],
}

/// Geometry observed from one accepted original thermal state.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanarLiquidGeometry {
    pub cells: Vec<PlanarLiquidCellGeometry>,
    pub footprint_area_m2: f64,
    pub total_mass_kg: f64,
    pub total_volume_m3: f64,
    pub total_thickness_m: f64,
}

/// Source-bound liquid volume and native implicit heat transport.
///
/// Prepare once at the coherent reference state, then retain this binding while
/// advancing its original timeline. No physical history lives in this adapter.
/// Geometry is observed from accepted enthalpies; reference cells stay immutable.
pub struct PlanarLiquidDensityColumn<'s, 'b, 'a> {
    bound: &'s BoundPhaseLayeredThermalBody<'b, 'a>,
    footprint_area_m2: f64,
    conduction: ConductionObserver<'a>,
    laws: Vec<LiquidTangent<'a>>,
}

#[derive(Clone, Copy)]
struct LiquidTangent<'a> {
    source: &'a PreparedPhaseConductivity<'a>,
}

impl<'b, 'a> BoundPhaseLayeredThermalBody<'b, 'a> {
    /// Bind a complete planar liquid column to the native implicit time owners.
    ///
    /// Every initial mass/rho(h) must agree with its original reference volume.
    /// Cross-section stays fixed while normal thickness changes; all cells must
    /// remain fully liquid. Bonded solids and coexistence need another closure.
    ///
    /// Only normal faces exchange heat. Original contacts occupy their full
    /// intercell faces, sidewalls are insulated, and selected end baths follow
    /// the moving ends at fixed area. Authored ambient htc applies only to the
    /// remaining uncovered end area. A dual-bath single cell retains the native
    /// independent reservoir/ambient channels without counting a bath twice.
    ///
    /// The original source pressure convention remains a caller obligation.
    /// This does not solve fluid motion or classify enthalpy as internal energy.
    pub fn planar_liquid_density_column<'s>(
        &'s self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PlanarLiquidDensityColumn<'s, 'b, 'a>, DepositionError> {
        checkpoint()?;
        let LayerStackGeometry::Planar { area_m2 } = self.body.geometry else {
            return Err(bad(
                "liquid column requires an explicitly fixed planar footprint",
            ));
        };
        observe_geometry(self, area_m2, thermal, true, &mut checkpoint)?;
        let observer = self.conduction_observer(thermal, &mut checkpoint)?;
        // An isolated single cell has no geometric transport channel. Its
        // native source/ambient-only solve still uses the same phase gate and
        // original caloric object; no dummy link or bath is introduced.
        let layers = if self.body.cell_count() > 1 || self.reservoirs.iter().any(Option::is_some) {
            self.cell_ids()
                .into_iter()
                .map(|cell_id| PlanarDensityLayer {
                    cell_id,
                    footprint_area_m2: area_m2,
                })
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let conduction =
            observer.with_material_density_layers(thermal, &layers, &mut checkpoint)?;
        let laws = self
            .body
            .layers
            .iter()
            .map(|layer| LiquidTangent {
                source: layer.conductivity,
            })
            .collect();
        checkpoint()?;
        Ok(PlanarLiquidDensityColumn {
            bound: self,
            footprint_area_m2: area_m2,
            conduction,
            laws,
        })
    }
}

impl PlanarLiquidDensityColumn<'_, '_, '_> {
    /// Build the original native BE/SDIRK/adaptive chain with exact source tangents.
    ///
    /// Both local and whole-call joule gates remain authoritative. Every trial,
    /// stage and accepted endpoint must remain inside the original fully liquid
    /// domain. A refusal or cancellation publishes no changed thermal state.
    /// Temporal and spatial refinement remain separate accuracy requirements.
    pub fn implicit_chain<'c>(
        &'c self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitConductionChain<'c>, ImplicitConductionError> {
        self.current_geometry(thermal, &mut checkpoint)?;
        let laws: Vec<_> = self
            .laws
            .iter()
            .enumerate()
            .map(|(index, law)| {
                (
                    MaterialId(self.bound.ids.first_material + index as u64),
                    law as &'c dyn ImplicitConductionLaw,
                )
            })
            .collect();
        ImplicitConductionChain::new(thermal, &self.conduction, &laws, checkpoint)
    }

    /// Resolve current positions and volumes without changing reference geometry.
    pub fn current_geometry(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PlanarLiquidGeometry, DepositionError> {
        observe_geometry(
            self.bound,
            self.footprint_area_m2,
            thermal,
            false,
            &mut checkpoint,
        )
    }
}

fn observe_geometry(
    bound: &BoundPhaseLayeredThermalBody<'_, '_>,
    area_m2: f64,
    thermal: &DepositionSimulation<'_>,
    initial: bool,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<PlanarLiquidGeometry, DepositionError> {
    checkpoint()?;
    if thermal.cells().len() != bound.body.cell_count()
        || thermal.links().len() != bound.body.links.len()
    {
        return Err(bad(
            "liquid column requires its complete original thermal body",
        ));
    }
    bound.admit_timeline(thermal, checkpoint)?;
    let mut cells = Vec::with_capacity(bound.body.cell_count());
    let (mut mass, mut volume, mut position) = (0.0, 0.0, 0.0);
    for (index, reference) in bound.body.cells.iter().enumerate() {
        checkpoint()?;
        let id = bound.cell_id(index);
        let cell = &thermal.cells()[&id];
        let state = bound.body.layers[reference.layer]
            .conductivity
            .caloric_curve()
            .values_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)?;
        if state.liquid_mass_fraction != 1.0 {
            return Err(DepositionError::MaterialDomain);
        }
        let cell_volume = finite(cell.mass_kg / state.bulk_density_kg_m3)?;
        if !positive(cell_volume)
            || initial
                && (cell_volume - reference.geometry.volume_m3).abs()
                    > 128.0 * f64::EPSILON * cell_volume.max(reference.geometry.volume_m3)
        {
            return Err(bad(
                "liquid mass and source density must match reference volume",
            ));
        }
        let thickness = finite(cell_volume / area_m2)?;
        let end = finite(position + thickness)?;
        if !positive(thickness) || end <= position {
            return Err(bad(
                "liquid column thickness is not geometrically representable",
            ));
        }
        mass = finite(mass + cell.mass_kg)?;
        volume = finite(volume + cell_volume)?;
        cells.push(PlanarLiquidCellGeometry {
            cell_id: id,
            mass_kg: cell.mass_kg,
            density_kg_m3: state.bulk_density_kg_m3,
            volume_m3: cell_volume,
            interval_m: [position, end],
        });
        position = end;
    }
    checkpoint()?;
    Ok(PlanarLiquidGeometry {
        cells,
        footprint_area_m2: area_m2,
        total_mass_kg: mass,
        total_volume_m3: volume,
        total_thickness_m: position,
    })
}

impl LiquidTangent<'_> {
    fn sample(
        &self,
        h: f64,
    ) -> Result<crate::thermal::conduction::phase::PhaseConductivitySample, DepositionError> {
        let sample = self
            .source
            .state_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        if sample.liquid_mass_fraction() != 1.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(sample)
    }
}

impl ImplicitConductionLaw for LiquidTangent<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.source.caloric_curve()
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self.source
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.source.enthalpy_bounds_j_kg()
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let sample = self.sample(h)?;
        Ok(ImplicitConductionState {
            temperature_k: sample.temperature_k(),
            temperature_slope_kg_k_j: sample.temperature_enthalpy_derivative_kg_k_j(),
            conductivity_w_m_k: sample.conductivity_w_m_k(),
            conductivity_slope_kg_w_m_k_j: sample.conductivity_enthalpy_derivative_kg_w_m_k_j(),
        })
    }
    fn specific_volume_derivative_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.sample(h)?;
        ImplicitConductionLaw::specific_volume_derivative_at_specific_enthalpy(self.source, h)
    }
    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        let sample = self.sample(h)?;
        Ok(ImplicitConductionObservables {
            temperature_k: sample.temperature_k(),
            phase_mass_fraction: Some(sample.liquid_mass_fraction()),
        })
    }
}
