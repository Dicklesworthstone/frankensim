//! Coherent angular interference deposited in the original thermal cells.
//!
//! The original prepared layer field supplies cycle-mean viscous and thermal
//! heat AFTER coherent pressure/velocity summation. Positive tensor quadrature
//! integrates these densities over an explicitly positioned rectangular column.
//! Preparation caches watts; native thermal advances never resample acoustics.
//! The caller must resolve both transverse interference and depth attenuation
//! with quadrature refinement. Positive weights guarantee nonnegative heat,
//! not spatial accuracy. No field or quadrature weight is renormalized.
//!
//! This is frozen, one-way heating of a finite volume in the original infinite
//! planar field. It adds neither finite-aperture scattering nor temperature
//! feedback. A cropped column can exchange acoustic power through its sides;
//! deposited heat therefore does not claim a two-face/global power balance.
//! Each native cell retains its own mass and single-temperature caloric law.
//! The mechanism split does not partition heat between fluid and skeleton.

use super::{MAX_POROUS_HEATING_CELLS, PorousCellHeating, PorousHeatCell, close_geometry, norm};
use fs_duct::porous::PorousLayerDissipation;
use fs_material::porous::layers::{AcousticSpectrumError, Layer, PreparedLayerAcousticSpectrum};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig,
};
use std::collections::BTreeMap;

/// Hard bound on component evaluations during preparation, including all cells.
pub const MAX_COHERENT_HEAT_COMPONENT_EVALUATIONS: usize = 16_777_216;

/// One positive unit-interval quadrature node. Weights integrate a unit length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticHeatNode {
    /// Coordinate in [0,1]; mapped to the actual axis interval in metres.
    pub coordinate: f64,
    /// Strictly positive, at most one; each axis's weights must sum to one.
    pub weight: f64,
}

/// Caller-authored integration rule, shared by every native depth cell.
#[derive(Clone, Copy, Debug)]
pub struct CoherentHeatQuadrature<'q> {
    /// Unit-interval nodes for acoustic x, y and cell-local depth, respectively.
    /// All three slices must be nonempty. Ordering is retained deterministically.
    pub axes: [&'q [AcousticHeatNode]; 3],
}

/// Explicit allowance for cached cells and aggregate preparation work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoherentHeatingBudget {
    pub max_cells: usize,
    /// `cells * x_nodes * y_nodes * depth_nodes * acoustic_components`.
    pub max_component_evaluations: usize,
}

/// Rigid correspondence between the acoustic interface frame and native cells.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentThermalColumn {
    /// Native position [m] of acoustic `[0,0,0]` at the SELECTED layer's front.
    /// Acoustic x/y retain the original stack's transverse phase origin.
    pub origin_m: [f64; 3],
    /// Orthonormal native directions of acoustic x and y. Their cross product
    /// is the layer's positive depth direction; no scale or reflection enters.
    pub transverse_axes: [[f64; 3]; 2],
    /// Lower acoustic x/y coordinates [m] of the rectangular footprint.
    pub lower_m: [f64; 2],
    /// Positive acoustic x/y extents [m]. Full geometric volume is used once.
    pub extent_m: [f64; 2],
}

/// Original acoustic/thermal refusal or inadmissible physical correspondence.
#[derive(Clone, Debug, PartialEq)]
pub enum CoherentHeatingError {
    Acoustic(AcousticSpectrumError),
    Thermal(DepositionError),
    Invalid(&'static str),
}
impl core::fmt::Display for CoherentHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coherent layer heating refused: {self:?}")
    }
}
impl std::error::Error for CoherentHeatingError {}
impl From<DepositionError> for CoherentHeatingError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<AcousticSpectrumError> for CoherentHeatingError {
    fn from(error: AcousticSpectrumError) -> Self {
        Self::Acoustic(error)
    }
}

/// Cached physical work bound to the unchanged original field and native cells.
pub struct CoherentLayerHeating<'f, 's, 'p, 'a> {
    field: &'f PreparedLayerAcousticSpectrum<'s, 'p, 'a>,
    cells: Vec<PorousCellHeating>,
    sources_w: BTreeMap<CellId, f64>,
}

/// Acoustic work actually credited over the accepted native thermal substeps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentHeatingEnergy {
    pub viscous_j: f64,
    pub thermal_j: f64,
    /// Already included exactly once in `heat_report.external_heat_j`.
    pub deposited_j: f64,
}

/// Native thermal candidate and separate diagnostic accounting of acoustic heat.
pub struct CoherentHeatingStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    /// Acoustic heat plus all native authored sources, contacts and reservoirs.
    pub heat_report: EnergyAudit,
    pub acoustic: CoherentHeatingEnergy,
    /// Deposited acoustic Joules in the caller's original cell/depth order.
    pub cell_energy_j: Vec<(CellId, f64)>,
}

fn finite(value: f64) -> Result<f64, CoherentHeatingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}

impl CoherentThermalColumn {
    fn admit(self) -> Result<([f64; 3], f64), CoherentHeatingError> {
        let [x, y] = self.transverse_axes;
        if self
            .origin_m
            .iter()
            .chain(x.iter())
            .chain(y.iter())
            .chain(self.lower_m.iter())
            .chain(self.extent_m.iter())
            .any(|v| !v.is_finite())
            || self.extent_m.iter().any(|&v| v <= 0.0)
            || !close_geometry(norm(x), 1.0)
            || !close_geometry(norm(y), 1.0)
            || (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]).abs() > 256.0 * f64::EPSILON
        {
            return Err(CoherentHeatingError::Invalid(
                "finite orthonormal column frame and positive footprint required",
            ));
        }
        for axis in 0..2 {
            if finite(self.lower_m[axis] + self.extent_m[axis])? <= self.lower_m[axis] {
                return Err(CoherentHeatingError::Invalid(
                    "unrepresentable footprint extent",
                ));
            }
        }
        let area = finite(self.extent_m[0] * self.extent_m[1])?;
        if area <= 0.0 {
            return Err(CoherentHeatingError::Invalid(
                "unrepresentable footprint area",
            ));
        }
        Ok((
            [
                x[1] * y[2] - x[2] * y[1],
                x[2] * y[0] - x[0] * y[2],
                x[0] * y[1] - x[1] * y[0],
            ],
            area,
        ))
    }

    fn center(self, normal: [f64; 3], depth: f64) -> [f64; 3] {
        let xy = [
            self.lower_m[0] + 0.5 * self.extent_m[0],
            self.lower_m[1] + 0.5 * self.extent_m[1],
        ];
        core::array::from_fn(|i| {
            self.origin_m[i]
                + xy[0] * self.transverse_axes[0][i]
                + xy[1] * self.transverse_axes[1][i]
                + depth * normal[i]
        })
    }
}

impl<'f, 's, 'p, 'a> CoherentLayerHeating<'f, 's, 'p, 'a> {
    /// Integrate the original coherent heat over a complete layer-depth column.
    ///
    /// Existing cell centrelines must be the footprint centres at each authored
    /// interval's endpoints. Their actual volumes must equal footprint area
    /// times depth; no porosity or inferred density multiplies heat. Every cell
    /// belongs to the same body and appears once. Geometry, work and normalized
    /// positive axis weights are admitted before any acoustic field observation.
    /// The caller owns spatial resolution; refine against the desired heat QoI.
    pub fn try_new(
        field: &'f PreparedLayerAcousticSpectrum<'s, 'p, 'a>,
        thermal: &DepositionSimulation<'_>,
        column: CoherentThermalColumn,
        partition: &[PorousHeatCell],
        quadrature: CoherentHeatQuadrature<'_>,
        budget: CoherentHeatingBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, CoherentHeatingError> {
        checkpoint()?;
        if budget.max_cells == 0
            || budget.max_cells > MAX_POROUS_HEATING_CELLS
            || partition.is_empty()
            || partition.len() > budget.max_cells
            || budget.max_component_evaluations == 0
            || budget.max_component_evaluations > MAX_COHERENT_HEAT_COMPONENT_EVALUATIONS
            || quadrature.axes.iter().any(|axis| axis.is_empty())
        {
            return Err(CoherentHeatingError::Invalid(
                "coherent heat cell or work budget exceeded",
            ));
        }
        let work = quadrature
            .axes
            .iter()
            .try_fold(partition.len(), |n, axis| n.checked_mul(axis.len()))
            .and_then(|n| n.checked_mul(field.spectrum().component_count()));
        if work.is_none_or(|work| work > budget.max_component_evaluations) {
            return Err(CoherentHeatingError::Invalid(
                "coherent heat component evaluation budget exceeded",
            ));
        }
        for axis in quadrature.axes {
            let mut sum = Dd::ZERO;
            for node in axis {
                checkpoint()?;
                if !node.coordinate.is_finite()
                    || !(0.0..=1.0).contains(&node.coordinate)
                    || !node.weight.is_finite()
                    || node.weight <= 0.0
                    || node.weight > 1.0
                {
                    return Err(CoherentHeatingError::Invalid(
                        "finite unit-interval nodes and positive weights required",
                    ));
                }
                sum = sum + Dd::from_f64(node.weight);
            }
            if (sum.to_f64() - 1.0).abs() > 64.0 * f64::EPSILON {
                return Err(CoherentHeatingError::Invalid(
                    "quadrature axis must integrate unit length without renormalization",
                ));
            }
        }
        let (normal, area) = column.admit()?;
        let length = match field
            .spectrum()
            .stack()
            .prepared_layers()
            .layer(field.layer_index())
        {
            Some(
                Layer::Porous { thickness_m, .. }
                | Layer::RelaxingPorous { thickness_m, .. }
                | Layer::FluidGap { thickness_m }
                | Layer::FluidLayer { thickness_m, .. }
                | Layer::RelaxingFluid { thickness_m, .. },
            ) => thickness_m,
            None => {
                return Err(CoherentHeatingError::Invalid(
                    "missing original acoustic layer",
                ));
            }
        };
        if length <= 0.0 || partition.last().unwrap().end_depth_m != length {
            return Err(CoherentHeatingError::Invalid(
                "complete positive original layer depth required",
            ));
        }
        let mut cells = Vec::with_capacity(partition.len());
        let mut sources_w = BTreeMap::new();
        let mut start = 0.0;
        let mut body = None;
        for route in partition {
            checkpoint()?;
            let end = route.end_depth_m;
            if !end.is_finite()
                || end <= start
                || end > length
                || sources_w.insert(route.cell, 0.0).is_some()
            {
                return Err(CoherentHeatingError::Invalid(
                    "distinct cells and increasing complete depth partition required",
                ));
            }
            let cell = thermal
                .cells()
                .get(&route.cell)
                .ok_or(CoherentHeatingError::Invalid(
                    "inactive coherent heat destination",
                ))?;
            let width = end - start;
            let volume = finite(area * width)?;
            let axis = core::array::from_fn(|i| cell.geometry.end_m[i] - cell.geometry.start_m[i]);
            if volume <= 0.0
                || !close_geometry(cell.geometry.volume_m3, volume)
                || !close_geometry(norm(axis), width)
                || cell.geometry.surface_area_m2 < 2.0 * area
                || body.is_some_and(|body| body != cell.body_id)
            {
                return Err(CoherentHeatingError::Invalid(
                    "column needs original prism volumes, faces and one body",
                ));
            }
            body = Some(cell.body_id);
            for (actual, depth) in [(cell.geometry.start_m, start), (cell.geometry.end_m, end)] {
                let expected = column.center(normal, depth);
                let delta = norm(core::array::from_fn(|i| actual[i] - expected[i]));
                let scale = norm(actual)
                    .max(norm(expected))
                    .max(width)
                    .max(column.extent_m[0])
                    .max(column.extent_m[1]);
                if !delta.is_finite() || delta > 256.0 * f64::EPSILON * scale {
                    return Err(CoherentHeatingError::Invalid(
                        "native cell endpoints do not match the acoustic column frame",
                    ));
                }
            }
            cells.push(PorousCellHeating {
                cell: route.cell,
                start_depth_m: start,
                end_depth_m: end,
                dissipation: PorousLayerDissipation {
                    viscous_power_w: 0.0,
                    thermal_power_w: 0.0,
                    absorbed_power_w: 0.0,
                },
                geometry: cell.geometry,
                material: cell.material_id,
                mass_kg: cell.mass_kg,
                body_id: cell.body_id,
                segment_id: cell.segment_id,
                storage_basis: thermal.material_storage_basis(cell.material_id)?,
                caloric_identity: thermal.material_caloric_identity(cell.material_id)?,
            });
            start = end;
        }
        for cell in &mut cells {
            let mut powers = [Dd::ZERO; 2];
            let mut represented_volume = Dd::ZERO;
            for x in quadrature.axes[0] {
                for y in quadrature.axes[1] {
                    for z in quadrature.axes[2] {
                        checkpoint()?;
                        let position = [
                            column.lower_m[0] + column.extent_m[0] * x.coordinate,
                            column.lower_m[1] + column.extent_m[1] * y.coordinate,
                            if z.coordinate <= 0.5 {
                                cell.start_depth_m
                                    + (cell.end_depth_m - cell.start_depth_m) * z.coordinate
                            } else {
                                cell.end_depth_m
                                    - (cell.end_depth_m - cell.start_depth_m) * (1.0 - z.coordinate)
                            },
                        ];
                        let dv = cell.geometry.volume_m3 * x.weight * y.weight * z.weight;
                        if !dv.is_finite() || dv <= 0.0 {
                            return Err(CoherentHeatingError::Invalid(
                                "unrepresentable positive quadrature volume",
                            ));
                        }
                        represented_volume = represented_volume + Dd::from_f64(dv);
                        let mut cancelled = None;
                        let sample = field.field_with_heat(position, &mut || {
                            if cancelled.is_none() {
                                cancelled = checkpoint().err();
                            }
                            cancelled.is_some()
                        });
                        if let Some(error) = cancelled {
                            return Err(error.into());
                        }
                        let sample = sample?;
                        for (sum, density) in powers
                            .iter_mut()
                            .zip([sample.viscous_heating_w_m3, sample.thermal_heating_w_m3])
                        {
                            let watts = finite(density * dv)?;
                            if watts < 0.0 || (density > 0.0 && watts == 0.0) {
                                return Err(CoherentHeatingError::Invalid(
                                    "unrepresentable positive quadrature heat",
                                ));
                            }
                            *sum = *sum + Dd::from_f64(watts);
                        }
                    }
                }
            }
            if !close_geometry(represented_volume.to_f64(), cell.geometry.volume_m3) {
                return Err(CoherentHeatingError::Invalid(
                    "quadrature does not preserve the original cell volume",
                ));
            }
            cell.dissipation = PorousLayerDissipation {
                viscous_power_w: finite(powers[0].to_f64())?,
                thermal_power_w: finite(powers[1].to_f64())?,
                absorbed_power_w: finite((powers[0] + powers[1]).to_f64())?,
            };
            *sources_w.get_mut(&cell.cell).unwrap() = cell.dissipation.absorbed_power_w;
        }
        checkpoint()?;
        Ok(Self {
            field,
            cells,
            sources_w,
        })
    }

    /// Exact original prepared field; no reconstructed material/source identity.
    pub const fn field(&self) -> &'f PreparedLayerAcousticSpectrum<'s, 'p, 'a> {
        self.field
    }

    /// Cached mechanism watts in the authored native cell/depth order.
    pub fn cells(&self) -> &[PorousCellHeating] {
        &self.cells
    }

    /// Hold these frozen cycle-mean watts through the original native timeline.
    /// Thermal contacts, reservoirs, authored sources, events, caloric domains
    /// and rollback remain native. Acoustic heat already has the phasor factor
    /// 1/2 and actual cell volume; do not also author it as a boundary source.
    /// Refresh the original acoustic/material owner when temperature changes
    /// invalidate its frozen state. This method performs no acoustic queries.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<CoherentHeatingStep<'t>, CoherentHeatingError> {
        checkpoint()?;
        validate(&self.cells, initial, &mut checkpoint)?;
        let observer = HeldHeat {
            cells: &self.cells,
            sources_w: &self.sources_w,
            energy: vec![Dd::ZERO; self.cells.len()],
            mechanisms: [Dd::ZERO; 2],
        };
        let (thermal, observer, heat_report) =
            initial.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        validate(&self.cells, &thermal, &mut checkpoint)?;
        let deposited = observer.energy.iter().copied().fold(Dd::ZERO, |a, b| a + b);
        let acoustic = CoherentHeatingEnergy {
            viscous_j: finite(observer.mechanisms[0].to_f64())?,
            thermal_j: finite(observer.mechanisms[1].to_f64())?,
            deposited_j: finite(deposited.to_f64())?,
        };
        let cell_energy_j = self
            .cells
            .iter()
            .zip(observer.energy)
            .map(|(cell, energy)| (cell.cell, energy.to_f64()))
            .collect();
        checkpoint()?;
        Ok(CoherentHeatingStep {
            thermal,
            heat_report,
            acoustic,
            cell_energy_j,
        })
    }
}

fn validate(
    cells: &[PorousCellHeating],
    thermal: &DepositionSimulation<'_>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), CoherentHeatingError> {
    for bound in cells {
        checkpoint()?;
        let cell = thermal
            .cells()
            .get(&bound.cell)
            .ok_or(CoherentHeatingError::Invalid(
                "coherent heat destination is no longer active",
            ))?;
        if cell.geometry != bound.geometry
            || cell.material_id != bound.material
            || cell.mass_kg != bound.mass_kg
            || cell.body_id != bound.body_id
            || cell.segment_id != bound.segment_id
            || thermal.material_storage_basis(cell.material_id)? != bound.storage_basis
            || thermal.material_caloric_identity(cell.material_id)? != bound.caloric_identity
        {
            return Err(CoherentHeatingError::Invalid(
                "coherent heat binding changed geometry, mass or original thermal material",
            ));
        }
    }
    Ok(())
}

#[derive(Clone)]
struct HeldHeat<'h> {
    cells: &'h [PorousCellHeating],
    sources_w: &'h BTreeMap<CellId, f64>,
    energy: Vec<Dd>,
    mechanisms: [Dd; 2],
}
impl DepositionObserver for HeldHeat<'_> {
    type Error = CoherentHeatingError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        validate(self.cells, thermal, checkpoint)?;
        Ok(self.sources_w.clone())
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
        for (cell, energy) in self.cells.iter().zip(&mut self.energy) {
            checkpoint()?;
            *energy = *energy + Dd::from_f64(finite(cell.dissipation.absorbed_power_w * dt)?);
            for (sum, watts) in self.mechanisms.iter_mut().zip([
                cell.dissipation.viscous_power_w,
                cell.dissipation.thermal_power_w,
            ]) {
                *sum = *sum + Dd::from_f64(finite(watts * dt)?);
            }
        }
        Ok(())
    }
}
