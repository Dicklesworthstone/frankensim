//! Volume and mechanical-wall absorption on the original thermal timeline.
//!
//! Every original acoustic tetrahedron maps explicitly to one existing thermal
//! control volume. The prepared equivalent-fluid operator supplies its physical
//! density and compressional loss forms. Optional mechanical-wall receivers use
//! the original resistance and Maxwell loss forms; pressure-transfer work,
//! unspecified boundary loss and numerical residuals never become heat. These
//! watts already include their physical measure and peak-to-mean scaling.
//!
//! The map is caller-authored: equal volumes do not establish spatial overlap.
//! Each destination must represent that complete physical element and its own
//! thermally lumped caloric law. The acoustic split does not assign separate
//! fluid/skeleton temperatures. [`InteriorAcousticHeating`] holds coefficients
//! and coherent pressure fixed during an advance. The [`feedback`] owner refreshes
//! the original material and complete acoustic solve at accepted temperatures.
//! Both use an explicit quasisteady closure, without inferred acoustic startup.

use fs_feec::acoustics::helmholtz::heterogeneous::{
    HeterogeneousError, MechanicalWallPhasors, PreparedEquivalentFluid,
};
use fs_math::c64::C64;
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Original molecular cylinder/cone viscous losses on finite gas heat cells.
pub mod duct;
/// Quasisteady full-volume acoustics refreshed from accepted thermal states.
pub mod feedback;
mod heat;
/// Resolved mixed-field viscous heating with separate thermal availability.
pub mod thermoviscous;
use heat::{InteriorHeatAccount, InteriorHeatState};

/// Hard storage/work ceiling; a caller may request a smaller cell budget.
pub const MAX_INTERIOR_HEATING_CELLS: usize = 4096;

/// Existing thermal destination, in the acoustic owner's declared cell order.
#[derive(Clone, Copy)]
pub struct InteriorHeatCell<'t> {
    /// Stable native thermal control-volume identity.
    pub cell: CellId,
    /// The exact borrowed caloric object used by this destination's timeline.
    pub caloric: &'t dyn ThermalLaw,
}

/// Existing thermally lumped wall destination, in original mechanical-wall order.
#[derive(Clone, Copy)]
pub struct InteriorWallHeatCell<'t> {
    pub cell: CellId,
    /// The exact caloric object used by this destination's native timeline.
    pub caloric: &'t dyn ThermalLaw,
    /// Complete acoustic contact area represented by this cell [m²].
    ///
    /// Must equal the original FE trace area. This is distinct from the native
    /// cell's exposed convection area, which may include both faces and edges.
    pub acoustic_area_m2: f64,
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

/// Independently evaluated heat from one original mechanical wall patch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InteriorWallHeating {
    pub cell: CellId,
    pub wall: usize,
    pub area_m2: f64,
    pub mass_kg: f64,
    pub resistance_dissipation_w: f64,
    pub relaxation_dissipation_w: f64,
}

impl InteriorWallHeating {
    #[must_use]
    pub fn absorbed_power_w(&self) -> f64 {
        self.resistance_dissipation_w + self.relaxation_dissipation_w
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

struct CellBinding<'t> {
    cell: CellId,
    caloric: &'t dyn ThermalLaw,
    geometry: CellGeometry,
    material: MaterialId,
    mass_kg: f64,
    body_id: u64,
    segment_id: u64,
    storage_basis: ThermalStorageBasis,
    caloric_identity: Option<[u8; 32]>,
}

impl<'t> CellBinding<'t> {
    fn try_new(
        thermal: &DepositionSimulation<'_>,
        id: CellId,
        caloric: &'t dyn ThermalLaw,
    ) -> Result<Self, InteriorHeatingError> {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(InteriorHeatingError::Invalid(
                "interior acoustic destination is not active",
            ))?;
        if !thermal.uses_material_law(cell.material_id, caloric) {
            return Err(InteriorHeatingError::Invalid(
                "interior heat needs the original native caloric object",
            ));
        }
        Ok(Self {
            cell: id,
            caloric,
            geometry: cell.geometry,
            material: cell.material_id,
            mass_kg: cell.mass_kg,
            body_id: cell.body_id,
            segment_id: cell.segment_id,
            storage_basis: thermal.material_storage_basis(cell.material_id)?,
            caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
        })
    }

    fn validate(&self, thermal: &DepositionSimulation<'_>) -> Result<(), InteriorHeatingError> {
        let cell = thermal
            .cells()
            .get(&self.cell)
            .ok_or(InteriorHeatingError::Invalid(
                "interior acoustic destination is no longer active",
            ))?;
        if cell.geometry != self.geometry
            || cell.material_id != self.material
            || cell.mass_kg != self.mass_kg
            || cell.body_id != self.body_id
            || cell.segment_id != self.segment_id
            || !thermal.uses_material_law(cell.material_id, self.caloric)
            || thermal.material_storage_basis(cell.material_id)? != self.storage_basis
            || thermal.material_caloric_identity(cell.material_id)? != self.caloric_identity
        {
            return Err(InteriorHeatingError::Invalid(
                "interior heat binding changed geometry, mass or original caloric law",
            ));
        }
        Ok(())
    }
}

struct BoundCell<'t> {
    heat: InteriorCellHeating,
    binding: CellBinding<'t>,
}

struct BoundWall<'t> {
    heat: InteriorWallHeating,
    binding: CellBinding<'t>,
}

/// Cached physical sources retaining the original operator and pressure field.
///
/// Preparation evaluates the original local loss forms once. Thermal substeps
/// copy cached watts and accumulate their accepted native duration; they never
/// resample material, solve acoustics, or rebuild the source map by insertion.
pub struct InteriorAcousticHeating<'p, 'f, 'c, 't> {
    prepared: &'p PreparedEquivalentFluid<'f, 'c>,
    pressure_pa: &'p [C64],
    wall_phasors: Option<&'p [MechanicalWallPhasors]>,
    heat: Arc<InteriorHeatState<'t>>,
}

/// Acoustic-only work over the accepted native thermal substeps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InteriorHeatingEnergy {
    pub density_dissipation_j: f64,
    pub compressional_dissipation_j: f64,
    pub wall_resistance_dissipation_j: f64,
    pub wall_relaxation_dissipation_j: f64,
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
    /// Mechanical heat only, in original wall-patch order.
    pub wall_energy_j: Vec<(CellId, f64)>,
}

fn finite(value: f64) -> Result<f64, InteriorHeatingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}

fn same_measure(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 512.0 * f64::EPSILON * a.abs().max(b.abs())
}

impl<'p, 'f, 'c, 't> InteriorAcousticHeating<'p, 'f, 'c, 't> {
    /// Bind bulk absorption from every original element exactly once.
    ///
    /// Mechanical-wall absorption is included only by
    /// [`Self::try_new_with_mechanical_walls`], with explicit wall receivers.
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
        if maximum_cells == 0 || maximum_cells > MAX_INTERIOR_HEATING_CELLS || count > maximum_cells
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
            let binding = CellBinding::try_new(thermal, route.cell, route.caloric)?;
            if !same_measure(binding.geometry.volume_m3, power.volume_m3) {
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
            cells.push(BoundCell { heat, binding });
        }
        checkpoint()?;
        Ok(Self {
            prepared,
            pressure_pa,
            wall_phasors: None,
            heat: Arc::new(InteriorHeatState {
                cells,
                walls: Vec::new(),
                sources_w,
            }),
        })
    }

    /// Bind bulk and mechanical-wall absorption to distinct native receivers.
    ///
    /// Supply the complete coherent wall state from the same original operator,
    /// frequency and drive as `pressure_pa`, including its retained Maxwell
    /// extensions. The supplied field's equations are not recertified here.
    /// Every original wall maps once, in original order, to a separate cell whose
    /// mass equals the FE surface mass times its actual integrated trace area.
    /// The caller owns spatial/material correspondence and the justification for
    /// one temperature through each complete wall patch; area/mass agreement
    /// alone does not establish that correspondence. Native exposed area remains
    /// unchanged and need not equal the one-sided acoustic contact area.
    ///
    /// Wall heat is recomputed from the original positive resistance and Maxwell
    /// forms, never from pressure work or a balance residual. No second area,
    /// volume or peak-to-mean factor is applied. The budget counts both volume
    /// and wall destinations before either loss evaluation begins.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_mechanical_walls(
        prepared: &'p PreparedEquivalentFluid<'f, 'c>,
        pressure_pa: &'p [C64],
        wall_phasors: &'p [MechanicalWallPhasors],
        thermal: &DepositionSimulation<'_>,
        destinations: &[InteriorHeatCell<'t>],
        wall_destinations: &[InteriorWallHeatCell<'t>],
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, InteriorHeatingError> {
        checkpoint()?;
        let wall_count = prepared.wall_count();
        let count = prepared
            .geometry()
            .tetrahedra()
            .len()
            .checked_add(wall_count);
        if maximum_cells == 0
            || maximum_cells > MAX_INTERIOR_HEATING_CELLS
            || count.is_none_or(|count| count > maximum_cells)
        {
            return Err(DepositionError::Budget("interior acoustic heating cells").into());
        }
        if wall_destinations.len() != wall_count || wall_phasors.len() != wall_count {
            return Err(InteriorHeatingError::Invalid(
                "interior wall heat needs every original wall state and destination",
            ));
        }
        let mut source = Self::try_new(
            prepared,
            pressure_pa,
            thermal,
            destinations,
            maximum_cells,
            &mut checkpoint,
        )?;
        let mut refusal = None;
        let powers = prepared.wall_dissipation(wall_phasors, &mut || match checkpoint() {
            Ok(()) => false,
            Err(error) => {
                refusal = Some(error);
                true
            }
        });
        if let Some(error) = refusal {
            return Err(error.into());
        }
        let state =
            Arc::get_mut(&mut source.heat).expect("newly created interior heat is uniquely owned");
        state.walls.reserve(wall_count);
        for (wall, (route, power)) in wall_destinations.iter().zip(powers?).enumerate() {
            checkpoint()?;
            if state.sources_w.contains_key(&route.cell) {
                return Err(InteriorHeatingError::Invalid(
                    "interior volume and wall destinations must be distinct",
                ));
            }
            let binding = CellBinding::try_new(thermal, route.cell, route.caloric)?;
            if !same_measure(route.acoustic_area_m2, power.area_m2)
                || !same_measure(binding.mass_kg, power.mass_kg)
            {
                return Err(InteriorHeatingError::Invalid(
                    "wall heat needs the complete original acoustic area and wall mass",
                ));
            }
            let heat = InteriorWallHeating {
                cell: route.cell,
                wall,
                area_m2: power.area_m2,
                mass_kg: power.mass_kg,
                resistance_dissipation_w: power.resistance_dissipation_w,
                relaxation_dissipation_w: power.relaxation_dissipation_w,
            };
            let watts = finite(heat.absorbed_power_w())?;
            if heat.resistance_dissipation_w < 0.0 || heat.relaxation_dissipation_w < 0.0 {
                return Err(InteriorHeatingError::Invalid("wall heat must be passive"));
            }
            state.sources_w.insert(route.cell, watts);
            state.walls.push(BoundWall { heat, binding });
        }
        source.wall_phasors = Some(wall_phasors);
        checkpoint()?;
        Ok(source)
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
        self.heat.cells()
    }

    /// Cached original mechanical-wall powers, separate from volume heat.
    pub fn walls(&self) -> impl ExactSizeIterator<Item = &InteriorWallHeating> {
        self.heat.walls()
    }

    #[must_use]
    pub const fn wall_phasors(&self) -> Option<&'p [MechanicalWallPhasors]> {
        self.wall_phasors
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), InteriorHeatingError> {
        self.heat.validate(thermal, checkpoint)
    }

    /// Start a composable native observer with a zero acoustic work account.
    /// Accepted observers retain accumulated energy across subsequent advances.
    pub fn observer(&self) -> InteriorHeatingObserver<'_, 'p, 'f, 'c, 't> {
        InteriorHeatingObserver {
            owner: self,
            account: InteriorHeatAccount::new(&self.heat),
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
        observer
            .account
            .finish(&self.heat, thermal, heat_report, &mut checkpoint)
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
    account: InteriorHeatAccount,
}

impl InteriorHeatingObserver<'_, '_, '_, '_, '_> {
    pub fn energy(&self) -> Result<InteriorHeatingEnergy, InteriorHeatingError> {
        self.account.energy()
    }

    pub fn cell_energy_j(&self) -> Result<Vec<(CellId, f64)>, InteriorHeatingError> {
        self.account.cell_energy_j(&self.owner.heat)
    }

    pub fn wall_energy_j(&self) -> Result<Vec<(CellId, f64)>, InteriorHeatingError> {
        self.account.wall_energy_j(&self.owner.heat)
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
        Ok(self.owner.heat.sources())
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
        self.account.record(&self.owner.heat, dt, checkpoint)
    }
}

#[cfg(test)]
mod tests;
