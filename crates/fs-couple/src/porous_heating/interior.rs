//! Full-volume harmonic absorption on the original finite thermal timeline.
//!
//! Every original acoustic tetrahedron maps explicitly to one existing thermal
//! control volume. The prepared equivalent-fluid operator supplies its physical
//! density and compressional loss forms; no boundary loss or numerical residual
//! becomes heat. These watts already include volume and peak-to-mean scaling.
//!
//! The map is caller-authored: equal volumes do not establish spatial overlap.
//! Each destination must represent that complete physical element and its own
//! thermally lumped caloric law. The acoustic split does not assign separate
//! fluid/skeleton temperatures. Coefficients and coherent pressure remain frozen
//! during an advance; refresh the acoustic/material owner when that approximation
//! ceases to hold. No acoustic startup or temperature feedback is inferred.

use fs_feec::acoustics::helmholtz::heterogeneous::{
    HeterogeneousError, PreparedEquivalentFluid,
};
use fs_math::{c64::C64, dd::Dd};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
};
use std::collections::BTreeMap;

/// Hard storage/work ceiling; a caller may request a smaller cell budget.
pub const MAX_INTERIOR_HEATING_CELLS: usize = 4096;

/// Existing thermal destination, in the original acoustic tetrahedron order.
#[derive(Clone, Copy)]
pub struct InteriorHeatCell<'t> {
    /// Stable native thermal control-volume identity.
    pub cell: CellId,
    /// The exact borrowed caloric object used by this destination's timeline.
    pub caloric: &'t dyn ThermalLaw,
}

/// Physical absorbed powers for one original element and destination.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InteriorCellHeating {
    pub cell: CellId,
    /// Original tetrahedron index; no sorting or spatial projection is applied.
    pub element: usize,
    pub volume_m3: f64,
    /// Density-loss work, already integrated and cycle averaged [W].
    pub density_dissipation_w: f64,
    /// Compressional-loss work, already integrated and cycle averaged [W].
    pub compressional_dissipation_w: f64,
}

impl InteriorCellHeating {
    #[must_use]
    pub fn absorbed_power_w(&self) -> f64 {
        self.density_dissipation_w + self.compressional_dissipation_w
    }
}

/// Original physics refusal or an inadmissible acoustic/thermal correspondence.
#[derive(Debug)]
pub enum InteriorHeatingError {
    Acoustic(HeterogeneousError),
    Thermal(DepositionError),
    Invalid(&'static str),
}

impl core::fmt::Display for InteriorHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "interior acoustic heating refused: {self:?}")
    }
}
impl std::error::Error for InteriorHeatingError {}
impl From<HeterogeneousError> for InteriorHeatingError {
    fn from(value: HeterogeneousError) -> Self {
        Self::Acoustic(value)
    }
}
impl From<DepositionError> for InteriorHeatingError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}

struct BoundCell<'t> {
    heat: InteriorCellHeating,
    caloric: &'t dyn ThermalLaw,
    geometry: CellGeometry,
    material: MaterialId,
    mass_kg: f64,
    body_id: u64,
    segment_id: u64,
    storage_basis: ThermalStorageBasis,
    caloric_identity: Option<[u8; 32]>,
}

/// Cached physical sources retaining the original operator and pressure field.
///
/// Preparation evaluates the original local loss forms once. Thermal substeps
/// copy cached watts and accumulate their accepted native duration; they never
/// resample material, solve acoustics, or rebuild the source map by insertion.
pub struct InteriorAcousticHeating<'p, 'f, 'c, 't> {
    prepared: &'p PreparedEquivalentFluid<'f, 'c>,
    pressure_pa: &'p [C64],
    cells: Vec<BoundCell<'t>>,
    sources_w: BTreeMap<CellId, f64>,
}

/// Acoustic-only work over the accepted native thermal substeps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InteriorHeatingEnergy {
    pub density_dissipation_j: f64,
    pub compressional_dissipation_j: f64,
    /// Already included once in the native thermal external-heat report.
    pub deposited_j: f64,
}

/// One native thermal candidate and its independent acoustic heat receipt.
pub struct InteriorHeatingStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    /// All original authored thermal exchanges plus acoustic heat exactly once.
    pub heat_report: EnergyAudit,
    pub acoustic: InteriorHeatingEnergy,
    /// Acoustic heat only, in original tetrahedron order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

fn finite(value: f64) -> Result<f64, InteriorHeatingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}

impl<'p, 'f, 'c, 't> InteriorAcousticHeating<'p, 'f, 'c, 't> {
    /// Bind the complete pressure space and every original element exactly once.
    ///
    /// The pressure may come from the prepared solve or an externally composed
    /// coherent field on this same space and frequency. This adapter evaluates
    /// its physical loss; it does not establish that an arbitrary supplied field
    /// satisfies the acoustic equations. Heat receipts cannot replace pressure.
    ///
    /// Destination volumes must match the original quadrature-integrated element
    /// volumes to floating-point roundoff, including curved geometry. The caller
    /// owns the spatial correspondence and uniform-temperature justification.
    /// Native mass, caloric identity, storage basis and body/segment bindings are
    /// retained. No porosity, density, volume or factor of one half is added.
    pub fn try_new(
        prepared: &'p PreparedEquivalentFluid<'f, 'c>,
        pressure_pa: &'p [C64],
        thermal: &DepositionSimulation<'_>,
        destinations: &[InteriorHeatCell<'t>],
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, InteriorHeatingError> {
        checkpoint()?;
        let count = prepared.geometry().tetrahedra().len();
        if maximum_cells == 0
            || maximum_cells > MAX_INTERIOR_HEATING_CELLS
            || count > maximum_cells
        {
            return Err(DepositionError::Budget("interior acoustic heating cells").into());
        }
        if count == 0 || destinations.len() != count {
            return Err(InteriorHeatingError::Invalid(
                "interior heat needs exactly one destination per original tetrahedron",
            ));
        }
        let mut refusal = None;
        let powers = prepared.element_dissipation(pressure_pa, &mut || match checkpoint() {
            Ok(()) => false,
            Err(error) => {
                refusal = Some(error);
                true
            }
        });
        if let Some(error) = refusal {
            return Err(error.into());
        }
        let powers = powers?;
        let mut cells = Vec::with_capacity(count);
        let mut sources_w = BTreeMap::new();
        for (element, (route, power)) in destinations.iter().zip(powers).enumerate() {
            checkpoint()?;
            if sources_w.contains_key(&route.cell) {
                return Err(InteriorHeatingError::Invalid(
                    "interior acoustic destinations must be distinct",
                ));
            }
            let cell = thermal.cells().get(&route.cell).ok_or(
                InteriorHeatingError::Invalid("interior acoustic destination is not active"),
            )?;
            let native_volume = cell.geometry.volume_m3;
            let tolerance = 512.0 * f64::EPSILON * native_volume.abs().max(power.volume_m3.abs());
            if !native_volume.is_finite()
                || native_volume <= 0.0
                || !power.volume_m3.is_finite()
                || power.volume_m3 <= 0.0
                || (native_volume - power.volume_m3).abs() > tolerance
                || !thermal.uses_material_law(cell.material_id, route.caloric)
            {
                return Err(InteriorHeatingError::Invalid(
                    "interior heat needs the complete original volume and caloric object",
                ));
            }
            let heat = InteriorCellHeating {
                cell: route.cell,
                element,
                volume_m3: power.volume_m3,
                density_dissipation_w: power.density_dissipation_w,
                compressional_dissipation_w: power.compressional_dissipation_w,
            };
            let watts = finite(heat.absorbed_power_w())?;
            if heat.density_dissipation_w < 0.0 || heat.compressional_dissipation_w < 0.0 {
                return Err(InteriorHeatingError::Invalid(
                    "interior acoustic element heat must be passive",
                ));
            }
            sources_w.insert(route.cell, watts);
            cells.push(BoundCell {
                heat,
                caloric: route.caloric,
                geometry: cell.geometry,
                material: cell.material_id,
                mass_kg: cell.mass_kg,
                body_id: cell.body_id,
                segment_id: cell.segment_id,
                storage_basis: thermal.material_storage_basis(cell.material_id)?,
                caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
            });
        }
        checkpoint()?;
        Ok(Self {
            prepared,
            pressure_pa,
            cells,
            sources_w,
        })
    }

    #[must_use]
    pub const fn prepared(&self) -> &'p PreparedEquivalentFluid<'f, 'c> {
        self.prepared
    }

    #[must_use]
    pub const fn pressure_pa(&self) -> &'p [C64] {
        self.pressure_pa
    }

    /// Cached original element powers in caller-authored destination order.
    pub fn cells(&self) -> impl ExactSizeIterator<Item = &InteriorCellHeating> {
        self.cells.iter().map(|cell| &cell.heat)
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), InteriorHeatingError> {
        for bound in &self.cells {
            checkpoint()?;
            let cell = thermal.cells().get(&bound.heat.cell).ok_or(
                InteriorHeatingError::Invalid("interior acoustic destination is no longer active"),
            )?;
            if cell.geometry != bound.geometry
                || cell.material_id != bound.material
                || cell.mass_kg != bound.mass_kg
                || cell.body_id != bound.body_id
                || cell.segment_id != bound.segment_id
                || !thermal.uses_material_law(cell.material_id, bound.caloric)
                || thermal.material_storage_basis(cell.material_id)? != bound.storage_basis
                || thermal.material_caloric_identity(cell.material_id)? != bound.caloric_identity
            {
                return Err(InteriorHeatingError::Invalid(
                    "interior heat binding changed geometry, mass or original caloric law",
                ));
            }
        }
        Ok(())
    }

    /// Start a composable native observer with a zero acoustic work account.
    /// Accepted observers retain accumulated energy across subsequent advances.
    pub fn observer(&self) -> InteriorHeatingObserver<'_, 'p, 'f, 'c, 't> {
        InteriorHeatingObserver {
            owner: self,
            cell_energy: vec![Dd::ZERO; self.cells.len()],
            mechanisms: [Dd::ZERO; 2],
        }
    }

    /// Advance the original thermal timeline with this held cycle-mean source.
    ///
    /// Original contacts, reservoirs, authored sources, events, caloric domains,
    /// latent storage and timestep admission remain in the native thermal owner.
    /// Do not also author these watts as boundary sources. A refusal publishes
    /// neither the thermal candidate nor its observer's accumulated energy.
    pub fn advance_to<'a>(
        &self,
        initial: &DepositionSimulation<'a>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<InteriorHeatingStep<'a>, InteriorHeatingError> {
        self.validate(initial, &mut checkpoint)?;
        let (thermal, observer, heat_report) =
            initial.advance_observed(&self.observer(), target_s, config, &mut checkpoint)?;
        self.validate(&thermal, &mut checkpoint)?;
        let acoustic = observer.energy()?;
        let cell_energy_j = observer.cell_energy_j()?;
        checkpoint()?;
        Ok(InteriorHeatingStep {
            thermal,
            heat_report,
            acoustic,
            cell_energy_j,
        })
    }
}

/// Original native observer, suitable for composition with other source owners.
///
/// Cloning copies only bounded energy accounts and retains the immutable source
/// owner. The native timeline clones it transactionally; failed work never
/// modifies the input observer or its thermal state.
#[derive(Clone)]
pub struct InteriorHeatingObserver<'a, 'p, 'f, 'c, 't> {
    owner: &'a InteriorAcousticHeating<'p, 'f, 'c, 't>,
    cell_energy: Vec<Dd>,
    mechanisms: [Dd; 2],
}

impl InteriorHeatingObserver<'_, '_, '_, '_, '_> {
    pub fn energy(&self) -> Result<InteriorHeatingEnergy, InteriorHeatingError> {
        let deposited = self.cell_energy.iter().copied().fold(Dd::ZERO, |a, b| a + b);
        Ok(InteriorHeatingEnergy {
            density_dissipation_j: finite(self.mechanisms[0].to_f64())?,
            compressional_dissipation_j: finite(self.mechanisms[1].to_f64())?,
            deposited_j: finite(deposited.to_f64())?,
        })
    }

    pub fn cell_energy_j(&self) -> Result<Vec<(CellId, f64)>, InteriorHeatingError> {
        self.owner
            .cells
            .iter()
            .zip(&self.cell_energy)
            .map(|(cell, energy)| Ok((cell.heat.cell, finite(energy.to_f64())?)))
            .collect()
    }
}

impl DepositionObserver for InteriorHeatingObserver<'_, '_, '_, '_, '_> {
    type Error = InteriorHeatingError;

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        Ok(self.owner.sources_w.clone())
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.owner.validate(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        for (cell, energy) in self.owner.cells.iter().zip(&mut self.cell_energy) {
            checkpoint()?;
            *energy = *energy + Dd::from_f64(finite(cell.heat.absorbed_power_w() * dt)?);
            for (sum, power) in self.mechanisms.iter_mut().zip([
                cell.heat.density_dissipation_w,
                cell.heat.compressional_dissipation_w,
            ]) {
                *sum = *sum + Dd::from_f64(finite(power * dt)?);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
