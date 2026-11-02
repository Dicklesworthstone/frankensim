//! Coherent porous acoustic loss on the original finite thermal timeline.
//!
//! An explicitly ordered, complete depth partition receives the original
//! `PorousLayerExcitation`'s analytic viscous and thermal work. Thermal cells
//! are adjoining uniform prisms on that layer's normal axis, with the SAME full
//! face area and specimen volume. Their masses, caloric laws, conduction,
//! reservoirs, authored sources and events remain owned by `fs-thermal`.
//!
//! Each cell is an explicitly authored thermally lumped specimen volume. The
//! caller must justify one temperature/storage law for its fluid and skeleton:
//! the acoustic mechanism split does NOT partition heat between constituents.
//! The fixed-source adapters hold acoustic coefficients and the coherent
//! peak-pressure drive over an advance. They provide one-way cycle-mean heat,
//! not acoustic startup,
//! temperature feedback, frame motion or a new material-validity claim. Refresh
//! the original acoustic/material owner explicitly when that approximation ends.
//! [`feedback`] instead refreshes the original gas/porous response from accepted
//! native cell temperatures under an explicit isobaric background.

use fs_duct::DuctError;
use fs_duct::porous::{PorousLayerDissipation, PorousLayerExcitation};
use fs_math::{dd::Dd, det};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalStorageBasis,
};
use std::collections::BTreeMap;

mod stack;
pub use stack::PorousStackHeating;
/// Accepted-state gas-temperature feedback through the original porous stack.
pub mod feedback;

/// Hard preparation/storage ceiling; callers may declare a smaller cell budget.
pub const MAX_POROUS_HEATING_CELLS: usize = 4096;

/// One explicit thermal destination in incident-face-first physical order.
/// Its interval starts at the preceding endpoint, or zero for the first cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousHeatCell {
    /// Existing stable thermal cell identity.
    pub cell: CellId,
    /// Far face depth [m]. The final endpoint must equal the layer thickness.
    pub end_depth_m: f64,
}

/// Analytic physical heat cached once for one original thermal control volume.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousCellHeating {
    /// Original thermal destination, in authored depth order.
    pub cell: CellId,
    /// Near face depth [m] from the acoustic incident side.
    pub start_depth_m: f64,
    /// Far face depth [m] from the same origin.
    pub end_depth_m: f64,
    /// Already cycle averaged and integrated over the full specimen volume.
    pub dissipation: PorousLayerDissipation,
    geometry: CellGeometry,
    material: MaterialId,
    mass_kg: f64,
    body_id: u64,
    segment_id: u64,
    storage_basis: ThermalStorageBasis,
    caloric_identity: Option<[u8; 32]>,
}

/// Original acoustic or thermal refusal, or an invalid physical correspondence.
#[derive(Clone, Debug, PartialEq)]
pub enum PorousHeatingError {
    /// Unchanged original acoustic refusal.
    Acoustic(DuctError),
    /// Original gas or porous constitutive/source-domain refusal.
    Material(fs_material::MaterialError),
    /// Unchanged native thermal, domain, budget or cancellation refusal.
    Thermal(DepositionError),
    /// Physical mapping or retained binding is inadmissible.
    Invalid(&'static str),
    /// Numerical power closure failed; no cell heat is rescaled to conceal it.
    PowerBalance { residual_w: f64, tolerance_w: f64 },
}
impl core::fmt::Display for PorousHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "porous acoustic heating refused: {self:?}")
    }
}
impl std::error::Error for PorousHeatingError {}
impl From<DuctError> for PorousHeatingError {
    fn from(value: DuctError) -> Self {
        Self::Acoustic(value)
    }
}
impl From<DepositionError> for PorousHeatingError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<fs_material::MaterialError> for PorousHeatingError {
    fn from(value: fs_material::MaterialError) -> Self {
        Self::Material(value)
    }
}
impl From<fs_duct::fluid::FluidDuctError> for PorousHeatingError {
    fn from(value: fs_duct::fluid::FluidDuctError) -> Self {
        match value {
            fs_duct::fluid::FluidDuctError::Material(error) => Self::Material(error),
            fs_duct::fluid::FluidDuctError::Duct(error) => Self::Acoustic(error),
        }
    }
}

/// Reusable source bound to the original acoustic solution and existing cells.
/// Preparation uses one constant-work analytic integral per cell and bounded
/// cell-ID map construction. Thermal substeps copy cached watts; they never
/// solve or sample acoustics, and never rebuild the source map by insertion.
pub struct PorousHeating<'e, 'l, 'm> {
    excitation: &'e PorousLayerExcitation<'l, 'm>,
    heat: PreparedHeat,
}

/// The same cached native sources and observer serve a layer or a whole stack.
#[derive(Clone)]
struct PreparedHeat {
    cells: Vec<PorousCellHeating>,
    sources_w: BTreeMap<CellId, f64>,
    dissipation: PorousLayerDissipation,
    power_residual_w: f64,
    boundary_power_w: [f64; 2],
}

/// Acoustic energy over the exact accepted native thermal substeps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousHeatingEnergy {
    /// Integrated prescribed incoming acoustic power [J].
    pub incident_j: f64,
    /// Integrated coherent outgoing acoustic power [J]; never deposited as heat.
    pub outgoing_j: f64,
    /// Original pore-viscous work [J].
    pub viscous_j: f64,
    /// Original thermal-relaxation work [J].
    pub thermal_j: f64,
    /// Positive heat delivered once to the bound cells, already in heat_report.
    pub deposited_j: f64,
    /// Incident minus outgoing minus deposited; never used to manufacture heat.
    pub residual_j: f64,
}

/// Native candidate and both distinct energy accounts. No accepted input state
/// changes on geometry, material, energy, work-budget or cancellation refusal.
pub struct PorousHeatingStep<'t> {
    /// Accepted native thermal candidate, retaining its authored boundaries.
    pub thermal: DepositionSimulation<'t>,
    /// Includes acoustic heat AND all original authored thermal exchanges.
    pub heat_report: EnergyAudit,
    /// Additional held acoustic drive's independent energy account.
    pub acoustic: PorousHeatingEnergy,
    /// Acoustic heat [J] only, in the original authored depth order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

fn finite(value: f64) -> Result<f64, PorousHeatingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn norm(value: [f64; 3]) -> f64 {
    det::hypot(det::hypot(value[0], value[1]), value[2])
}
fn close_geometry(a: f64, b: f64) -> bool {
    (a - b).abs() <= 256.0 * f64::EPSILON * a.abs().max(b.abs())
}
fn close_power(a: f64, b: f64, scale: f64, count: usize) -> Result<f64, PorousHeatingError> {
    let residual_w = finite(a - b)?;
    let tolerance_w = finite(512.0 * f64::EPSILON * (count + 1) as f64 * scale)?;
    if residual_w.abs() > tolerance_w {
        Err(PorousHeatingError::PowerBalance {
            residual_w,
            tolerance_w,
        })
    } else {
        Ok(residual_w)
    }
}

impl<'e, 'l, 'm> PorousHeating<'e, 'l, 'm> {
    /// Bind every depth exactly once to distinct already-active thermal cells.
    ///
    /// Cell centreline endpoints are the centres of opposite prismatic faces:
    /// `start_m` faces the acoustic incident side and `end_m` the far side.
    /// Adjacent faces must coincide, all axes must be parallel, and each length
    /// and volume must match its depth interval and the layer's full face area
    /// to floating-point roundoff. No porosity factor or inferred density enters.
    /// Cells share one body; their original material IDs/masses remain explicit.
    /// The stored caloric identity is also checked whenever the law supplies it.
    ///
    /// The caller declares a positive preparation budget no larger than the
    /// hard ceiling. Geometry and acoustic power admission precede publication.
    pub fn try_new(
        excitation: &'e PorousLayerExcitation<'l, 'm>,
        thermal: &DepositionSimulation<'_>,
        partition: &[PorousHeatCell],
        maximum_cells: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PorousHeatingError> {
        checkpoint()?;
        if maximum_cells == 0
            || maximum_cells > MAX_POROUS_HEATING_CELLS
            || partition.is_empty()
            || partition.len() > maximum_cells
        {
            return Err(DepositionError::Budget("porous heating cells").into());
        }
        let length = excitation.layer().thickness_m();
        let area = excitation.layer().area_m2();
        if partition.last().expect("nonempty partition").end_depth_m != length {
            return Err(PorousHeatingError::Invalid(
                "porous heat must cover the complete layer",
            ));
        }
        let mut cells: Vec<PorousCellHeating> = Vec::with_capacity(partition.len());
        let mut sources_w = BTreeMap::new();
        let mut start = 0.0;
        let mut axis = [0.0; 3];
        let (mut viscous, mut relaxation, mut absorbed) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for route in partition {
            checkpoint()?;
            if !route.end_depth_m.is_finite()
                || route.end_depth_m <= start
                || route.end_depth_m > length
                || sources_w.contains_key(&route.cell)
            {
                return Err(PorousHeatingError::Invalid(
                    "porous heat needs ordered nonoverlapping depths and unique cells",
                ));
            }
            let cell = thermal
                .cells()
                .get(&route.cell)
                .ok_or(PorousHeatingError::Invalid(
                    "porous heat destination is not active",
                ))?;
            let width = route.end_depth_m - start;
            let delta = std::array::from_fn(|i| cell.geometry.end_m[i] - cell.geometry.start_m[i]);
            let physical_length = finite(norm(delta))?;
            let volume = finite(area * width)?;
            if physical_length <= 0.0
                || volume <= 0.0
                || !close_geometry(physical_length, width)
                || !close_geometry(cell.geometry.volume_m3, volume)
                || cell.geometry.surface_area_m2 < 2.0 * area * (1.0 - 256.0 * f64::EPSILON)
            {
                return Err(PorousHeatingError::Invalid(
                    "porous heat requires the original full-area prismatic cell geometry",
                ));
            }
            if let Some(previous) = cells.last() {
                if previous.geometry.end_m != cell.geometry.start_m
                    || previous.body_id != cell.body_id
                {
                    return Err(PorousHeatingError::Invalid(
                        "porous heat cells must adjoin in one physical body",
                    ));
                }
                let misalignment = norm(std::array::from_fn(|i| delta[i] - axis[i] * width));
                if misalignment > 256.0 * f64::EPSILON * width {
                    return Err(PorousHeatingError::Invalid(
                        "porous heat cells must share the layer normal",
                    ));
                }
            } else {
                axis = delta.map(|v| v / physical_length);
            }
            let dissipation = excitation.dissipation_between(start, route.end_depth_m)?;
            sources_w.insert(route.cell, dissipation.absorbed_power_w);
            viscous = viscous + Dd::from_f64(dissipation.viscous_power_w);
            relaxation = relaxation + Dd::from_f64(dissipation.thermal_power_w);
            absorbed = absorbed + Dd::from_f64(dissipation.absorbed_power_w);
            cells.push(PorousCellHeating {
                cell: route.cell,
                start_depth_m: start,
                end_depth_m: route.end_depth_m,
                dissipation,
                geometry: cell.geometry,
                material: cell.material_id,
                mass_kg: cell.mass_kg,
                body_id: cell.body_id,
                segment_id: cell.segment_id,
                storage_basis: thermal.material_storage_basis(cell.material_id)?,
                caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
            });
            start = route.end_depth_m;
        }
        let dissipation = PorousLayerDissipation {
            viscous_power_w: finite(viscous.to_f64())?,
            thermal_power_w: finite(relaxation.to_f64())?,
            absorbed_power_w: finite(absorbed.to_f64())?,
        };
        let whole = excitation.dissipation_between(0.0, length)?;
        close_power(
            whole.absorbed_power_w,
            dissipation.absorbed_power_w,
            whole.absorbed_power_w.max(dissipation.absorbed_power_w),
            cells.len(),
        )?;
        let boundary = excitation.boundary_power();
        let power_residual_w = close_power(
            boundary.absorbed_power_w,
            dissipation.absorbed_power_w,
            boundary
                .incident_power_w
                .max(boundary.outgoing_power_w)
                .max(dissipation.absorbed_power_w),
            cells.len(),
        )?;
        checkpoint()?;
        Ok(Self {
            excitation,
            heat: PreparedHeat {
                cells,
                sources_w,
                dissipation,
                power_residual_w,
                boundary_power_w: [boundary.incident_power_w, boundary.outgoing_power_w],
            },
        })
    }

    /// Original coherent acoustic field and material response, unchanged.
    pub const fn excitation(&self) -> &'e PorousLayerExcitation<'l, 'm> {
        self.excitation
    }
    /// Cached analytic sources in the complete authored depth order.
    pub fn cells(&self) -> &[PorousCellHeating] {
        &self.heat.cells
    }
    /// Compensated sum of physical cell powers, with the original mechanism split.
    pub const fn dissipation(&self) -> PorousLayerDissipation {
        self.heat.dissipation
    }
    /// Incoming minus outgoing minus summed physical cell heating [W].
    pub const fn power_residual_w(&self) -> f64 {
        self.heat.power_residual_w
    }

    /// Apply this frozen cycle-mean drive over the native accepted substeps.
    /// Generated watts are additional sources, never persistent boundary data.
    /// Their harmonic factor 1/2 and specimen volume are already included; the
    /// native thermal account includes them exactly once. Existing conduction,
    /// source/boundary events, caloric domains, CFL and latent heat are unchanged.
    /// Do not also author these same acoustic watts as a thermal boundary source.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousHeatingStep<'t>, PorousHeatingError> {
        self.heat.advance_to(initial, target_s, config, checkpoint)
    }
}

impl PreparedHeat {
    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        for bound in &self.cells {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&bound.cell)
                .ok_or(PorousHeatingError::Invalid(
                    "porous heat destination is no longer active",
                ))?;
            if cell.geometry != bound.geometry
                || cell.material_id != bound.material
                || cell.mass_kg != bound.mass_kg
                || cell.body_id != bound.body_id
                || cell.segment_id != bound.segment_id
                || thermal.material_storage_basis(cell.material_id)? != bound.storage_basis
                || thermal.material_caloric_identity(cell.material_id)? != bound.caloric_identity
            {
                return Err(PorousHeatingError::Invalid(
                    "porous heat binding changed geometry, mass or original thermal material",
                ));
            }
        }
        Ok(())
    }

    /// Apply this frozen cycle-mean drive over the native accepted substeps.
    /// Generated watts are additional sources, never persistent boundary data.
    /// Their harmonic factor 1/2 and specimen volume are already included; the
    /// native thermal account includes them exactly once. Existing conduction,
    /// source/boundary events, caloric domains, CFL and latent heat are unchanged.
    /// Do not also author these same acoustic watts as a thermal boundary source.
    fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousHeatingStep<'t>, PorousHeatingError> {
        checkpoint()?;
        self.validate(initial, &mut checkpoint)?;
        let observer = HeldHeat {
            owner: self,
            account: HeatAccount::new(self.cells.len()),
        };
        let (thermal, observer, heat_report) =
            initial.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        observer
            .account
            .finish(self, thermal, heat_report, &mut checkpoint)
    }
}

/// One physical energy accumulator serves frozen and temperature-refreshed heat.
#[derive(Clone)]
struct HeatAccount {
    energy: Vec<Dd>,
    mechanisms: [Dd; 2],
    boundary: [Dd; 2],
}
impl HeatAccount {
    fn new(cells: usize) -> Self {
        Self {
            energy: vec![Dd::ZERO; cells],
            mechanisms: [Dd::ZERO; 2],
            boundary: [Dd::ZERO; 2],
        }
    }
    fn record(
        &mut self,
        heat: &PreparedHeat,
        dt: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        for (cell, energy) in heat.cells.iter().zip(&mut self.energy) {
            checkpoint()?;
            *energy = *energy + Dd::from_f64(finite(cell.dissipation.absorbed_power_w * dt)?);
            for (sum, power) in self.mechanisms.iter_mut().zip([
                cell.dissipation.viscous_power_w,
                cell.dissipation.thermal_power_w,
            ]) {
                *sum = *sum + Dd::from_f64(finite(power * dt)?);
            }
        }
        for (sum, power) in self.boundary.iter_mut().zip(heat.boundary_power_w) {
            *sum = *sum + Dd::from_f64(finite(power * dt)?);
        }
        Ok(())
    }
    fn acoustic_energy(&self) -> Result<PorousHeatingEnergy, PorousHeatingError> {
        let deposited = self.energy.iter().copied().fold(Dd::ZERO, |a, b| a + b);
        Ok(PorousHeatingEnergy {
            incident_j: finite(self.boundary[0].to_f64())?,
            outgoing_j: finite(self.boundary[1].to_f64())?,
            viscous_j: finite(self.mechanisms[0].to_f64())?,
            thermal_j: finite(self.mechanisms[1].to_f64())?,
            deposited_j: finite(deposited.to_f64())?,
            residual_j: finite((self.boundary[0] - self.boundary[1] - deposited).to_f64())?,
        })
    }
    fn finish<'t>(
        self,
        heat: &PreparedHeat,
        thermal: DepositionSimulation<'t>,
        heat_report: EnergyAudit,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<PorousHeatingStep<'t>, PorousHeatingError> {
        let acoustic = self.acoustic_energy()?;
        let cell_energy_j = heat
            .cells
            .iter()
            .zip(self.energy)
            .map(|(cell, energy)| (cell.cell, energy.to_f64()))
            .collect();
        checkpoint()?;
        Ok(PorousHeatingStep {
            thermal,
            heat_report,
            acoustic,
            cell_energy_j,
        })
    }
}

#[derive(Clone)]
struct HeldHeat<'p> {
    owner: &'p PreparedHeat,
    account: HeatAccount,
}
impl DepositionObserver for HeldHeat<'_> {
    type Error = PorousHeatingError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.owner.validate(thermal, checkpoint)?;
        Ok(self.owner.sources_w.clone())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.account.record(self.owner, dt, checkpoint)
    }
}
