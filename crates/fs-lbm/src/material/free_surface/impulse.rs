//! Zero-duration physical impulses on the original free-surface liquid owners.
//! Inputs are impulse density [N s/m³] per whole geometric cell, including
//! partially occupied cells. Native retained liquid mass supplies the inertia;
//! neither the timestep nor an assumed full-cell mass is used for the kick.

use super::{
    MaterialFreeSurface2, MaterialFreeSurface3, MaterialFreeSurfaceScaling,
    MaterialLbmError, finite, positive,
};
use crate::impulse::{ImpulseControls, ImpulseError, ImpulseReport, poll};
use fs_exec::CancelGate;

/// Explicit physical admission limits for a three-dimensional liquid impulse.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceImpulseControls3 {
    /// Original domain cell count, including gas, inactive cells and walls.
    pub max_cells: usize,
    /// Pre/post speed bound [m/s], no larger than the original scenario bound.
    pub maximum_speed_m_s: f64,
    /// Sum of local and global absolute momentum defects [kg m/s].
    pub maximum_momentum_defect_kg_m_s: f64,
    /// Sum of local and global absolute signed-work defects [J].
    pub maximum_work_defect_j: f64,
}

/// Explicit physical limits in two dimensions, per unit out-of-plane depth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceImpulseControls2 {
    /// Original domain cell count, including gas and wall rows.
    pub max_cells: usize,
    /// Pre/post speed bound [m/s], no larger than the original scenario bound.
    pub maximum_speed_m_s: f64,
    /// Sum of local and global absolute momentum defects [kg/s = N s/m].
    pub maximum_momentum_defect_kg_s: f64,
    /// Sum of local and global absolute signed-work defects [J/m].
    pub maximum_work_defect_j_m: f64,
}

/// Actual liquid impulse and signed work in three-dimensional SI units.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct MaterialFreeSurfaceImpulseReport3 {
    /// Prescribed impulse integrated over the affected geometric cells [N s].
    pub applied_impulse_kg_m_s: [f64; 3],
    /// Independently observed retained-liquid momentum change [kg m/s].
    pub momentum_change_kg_m_s: [f64; 3],
    /// Actual mass change [kg]; the impulse adds no mass.
    pub mass_change_kg: f64,
    /// Signed external work from impulse and mean endpoint velocity [J].
    pub work_j: f64,
    /// Independently observed liquid kinetic-energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Sum of local and global absolute momentum defects [kg m/s].
    pub absolute_momentum_defect_kg_m_s: f64,
    /// Sum of local and global absolute signed-work defects [J].
    pub absolute_work_defect_j: f64,
    /// Number of originally wet cells receiving a nonzero impulse.
    pub affected_cells: usize,
}

/// Actual liquid impulse and signed work per unit out-of-plane depth.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct MaterialFreeSurfaceImpulseReport2 {
    /// Prescribed impulse per depth [kg/s = N s/m].
    pub applied_impulse_kg_s: [f64; 2],
    /// Independently observed liquid momentum change per depth [kg/s].
    pub momentum_change_kg_s: [f64; 2],
    /// Actual mass change per depth [kg/m]; the impulse adds no mass.
    pub mass_change_kg_m: f64,
    /// Signed external work per depth [J/m].
    pub work_j_m: f64,
    /// Independently observed liquid kinetic-energy change per depth [J/m].
    pub kinetic_energy_change_j_m: f64,
    /// Sum of local and global absolute momentum defects [kg/s].
    pub absolute_momentum_defect_kg_s: f64,
    /// Sum of local and global absolute signed-work defects [J/m].
    pub absolute_work_defect_j_m: f64,
    /// Number of originally wet cells receiving a nonzero impulse.
    pub affected_cells: usize,
}

/// Physical-unit or native-state refusal; no populations have been published.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialFreeSurfaceImpulseError {
    /// Invalid or unrepresentable physical input, output or admission limit.
    Units(MaterialLbmError),
    /// Native shape, population, balance, allocation or cancellation refusal.
    Native(ImpulseError),
}
impl core::fmt::Display for MaterialFreeSurfaceImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Units(error) => write!(f, "physical liquid impulse: {error}"),
            Self::Native(error) => write!(f, "physical liquid impulse: {error}"),
        }
    }
}
impl std::error::Error for MaterialFreeSurfaceImpulseError {}
impl From<MaterialLbmError> for MaterialFreeSurfaceImpulseError {
    fn from(error: MaterialLbmError) -> Self { Self::Units(error) }
}
impl From<ImpulseError> for MaterialFreeSurfaceImpulseError {
    fn from(error: ImpulseError) -> Self { Self::Native(error) }
}

fn represented(input: f64, output: f64) -> Result<f64, MaterialLbmError> {
    finite(output, "physical impulse conversion is unrepresentable")?;
    if input != 0.0 && output == 0.0 {
        Err(MaterialLbmError::Invalid("physical impulse conversion underflows"))
    } else {
        Ok(output)
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare<const D: usize>(
    scaling: MaterialFreeSurfaceScaling<'_>,
    mass_unit: f64,
    impulses: &[[f64; D]],
    count: usize,
    max_cells: usize,
    limits: [f64; 3],
    gate: &CancelGate,
) -> Result<(Vec<[f64; D]>, ImpulseControls, [f64; 3]), MaterialFreeSurfaceImpulseError> {
    poll(gate)?;
    if impulses.len() != count {
        return Err(ImpulseError::Invalid("physical impulse field has the wrong shape").into());
    }
    let [speed, momentum_budget, work_budget] = limits;
    positive(speed, "physical impulse speed must be positive and finite")?;
    if !momentum_budget.is_finite() || momentum_budget < 0.0
        || !work_budget.is_finite() || work_budget < 0.0
    {
        return Err(MaterialLbmError::Invalid("physical impulse defects require finite nonnegative budgets").into());
    }
    let maximum_speed = represented(speed, speed / scaling.velocity_scale)?;
    if maximum_speed > scaling.flow.u_lattice {
        return Err(MaterialLbmError::Invalid("impulse speed exceeds the original scenario bound").into());
    }
    positive(mass_unit, "liquid mass unit is unrepresentable")?;
    let momentum_unit = positive(mass_unit * scaling.velocity_scale, "liquid momentum unit is unrepresentable")?;
    let work_unit = positive(momentum_unit * scaling.velocity_scale, "liquid work unit is unrepresentable")?;
    let native = ImpulseControls {
        max_cells,
        maximum_speed,
        maximum_momentum_defect: represented(momentum_budget, momentum_budget / momentum_unit)?,
        maximum_work_defect: represented(work_budget, work_budget / work_unit)?,
    };
    native.admit(count)?;
    let density_unit = positive(
        scaling.material.density_kg_m3() * scaling.velocity_scale,
        "liquid impulse-density unit is unrepresentable",
    )?;
    let mut field = Vec::new();
    field.try_reserve_exact(count)
        .map_err(|_| ImpulseError::Invalid("physical impulse allocation refused"))?;
    for (index, impulse) in impulses.iter().enumerate() {
        if index % 64 == 0 { poll(gate)?; }
        let mut converted = [0.0; D];
        for axis in 0..D {
            converted[axis] = represented(impulse[axis], impulse[axis] / density_unit)?;
        }
        field.push(converted);
    }
    Ok((field, native, [mass_unit, momentum_unit, work_unit]))
}

/// Private conversion container only: this native-shaped value never escapes
/// the adapter in physical units. Conversion and SI balance admission run inside
/// the native transaction, before its final cancellation poll and publication.
fn receipt_si<const D: usize>(
    mut report: ImpulseReport<D>,
    units: [f64; 3],
    momentum_budget: f64,
    work_budget: f64,
) -> Result<ImpulseReport<D>, MaterialLbmError> {
    let [mass, momentum, work] = units;
    let native_global_momentum = (0..D)
        .map(|axis| (report.momentum_change[axis] - report.applied_impulse[axis]).abs())
        .sum::<f64>();
    let local_momentum = report.absolute_momentum_defect - native_global_momentum;
    let local_work = report.absolute_work_defect - (report.kinetic_energy_change - report.work).abs();
    for axis in 0..D {
        report.applied_impulse[axis] = represented(report.applied_impulse[axis], report.applied_impulse[axis] * momentum)?;
        report.momentum_change[axis] = represented(report.momentum_change[axis], report.momentum_change[axis] * momentum)?;
    }
    report.mass_change = represented(report.mass_change, report.mass_change * mass)?;
    report.work = represented(report.work, report.work * work)?;
    report.kinetic_energy_change = represented(report.kinetic_energy_change, report.kinetic_energy_change * work)?;
    report.absolute_momentum_defect = finite(
        represented(local_momentum, local_momentum * momentum)?
            + (0..D).map(|axis| (report.momentum_change[axis] - report.applied_impulse[axis]).abs()).sum::<f64>(),
        "physical momentum readback is unrepresentable",
    )?;
    report.absolute_work_defect = finite(
        represented(local_work, local_work * work)?
            + (report.kinetic_energy_change - report.work).abs(),
        "physical work readback is unrepresentable",
    )?;
    if report.absolute_momentum_defect > momentum_budget
        || report.absolute_work_defect > work_budget
    {
        return Err(MaterialLbmError::Invalid("physical impulse readback exceeds its SI budget"));
    }
    Ok(report)
}

impl MaterialFreeSurface3<'_> {
    /// Apply one instantaneous impulse density [N s/m³] per original x-fastest
    /// geometric cell. Gas, inactive and wall entries must be zero. Partial
    /// liquid occupancy uses the actual retained mass, so the same impulse
    /// density causes a larger velocity jump at smaller fill.
    ///
    /// Preserves density, fill, held force, corrected viscous stress, source
    /// binding and elapsed time. No heat is added. The actual signed-work and
    /// momentum receipt is admitted in SI before native publication; any error
    /// or cancellation leaves the original state intact. Allocates one converted
    /// field and the native sparse proposal, without cloning the liquid owner.
    pub fn apply_impulse_density(
        &mut self,
        impulse_n_s_m3: &[[f64; 3]],
        controls: MaterialFreeSurfaceImpulseControls3,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceImpulseReport3, MaterialFreeSurfaceImpulseError> {
        let count = self.lattice.domain_cells().into_iter()
            .try_fold(1usize, |n, edge| n.checked_mul(edge))
            .ok_or(ImpulseError::Invalid("physical impulse cell count overflows"))?;
        let (field, native, units) = prepare(
            self.scaling, self.scaling.mass_unit_kg, impulse_n_s_m3, count,
            controls.max_cells,
            [controls.maximum_speed_m_s, controls.maximum_momentum_defect_kg_m_s, controls.maximum_work_defect_j],
            gate,
        )?;
        let mut refused_units = None;
        let result = self.lattice.apply_impulse_density_with_receipt(&field, native, gate, |report| {
            receipt_si(report, units, controls.maximum_momentum_defect_kg_m_s, controls.maximum_work_defect_j)
                .map_err(|error| {
                    refused_units = Some(error);
                    ImpulseError::Invalid("physical impulse receipt refused")
                })
        });
        if let Some(error) = refused_units { return Err(error.into()); }
        let report = result?;
        Ok(MaterialFreeSurfaceImpulseReport3 {
            applied_impulse_kg_m_s: report.applied_impulse,
            momentum_change_kg_m_s: report.momentum_change,
            mass_change_kg: report.mass_change,
            work_j: report.work,
            kinetic_energy_change_j: report.kinetic_energy_change,
            absolute_momentum_defect_kg_m_s: report.absolute_momentum_defect,
            absolute_work_defect_j: report.absolute_work_defect,
            affected_cells: report.affected_cells,
        })
    }
}

impl MaterialFreeSurface2<'_> {
    /// Apply an instantaneous impulse density [N s/m³] on the original
    /// x-fastest domain, including its wall rows. Non-liquid entries must be
    /// zero. The native held mass supplies partial-fill inertia.
    ///
    /// The receipt uses the original two-dimensional mass per depth; impulse
    /// is [N s/m] and work is [J/m]. No timestep is applied. Density, liquid
    /// mass, fill, held forces, corrected stress and source binding stay fixed.
    /// Physical conversion, balances and cancellation are admitted before the
    /// native transaction publishes populations, using one temporary input field.
    pub fn apply_impulse_density(
        &mut self,
        impulse_n_s_m3: &[[f64; 2]],
        controls: MaterialFreeSurfaceImpulseControls2,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceImpulseReport2, MaterialFreeSurfaceImpulseError> {
        let count = self.lattice.grid.nx.checked_mul(self.lattice.grid.ny)
            .ok_or(ImpulseError::Invalid("physical impulse cell count overflows"))?;
        let (field, native, units) = prepare(
            *self.scaling(), self.mass_per_depth_unit_kg_m, impulse_n_s_m3, count,
            controls.max_cells,
            [controls.maximum_speed_m_s, controls.maximum_momentum_defect_kg_s, controls.maximum_work_defect_j_m],
            gate,
        )?;
        let mut refused_units = None;
        let result = self.lattice.apply_impulse_density_with_receipt(&field, native, gate, |report| {
            receipt_si(report, units, controls.maximum_momentum_defect_kg_s, controls.maximum_work_defect_j_m)
                .map_err(|error| {
                    refused_units = Some(error);
                    ImpulseError::Invalid("physical impulse receipt refused")
                })
        });
        if let Some(error) = refused_units { return Err(error.into()); }
        let report = result?;
        Ok(MaterialFreeSurfaceImpulseReport2 {
            applied_impulse_kg_s: report.applied_impulse,
            momentum_change_kg_s: report.momentum_change,
            mass_change_kg_m: report.mass_change,
            work_j_m: report.work,
            kinetic_energy_change_j_m: report.kinetic_energy_change,
            absolute_momentum_defect_kg_s: report.absolute_momentum_defect,
            absolute_work_defect_j_m: report.absolute_work_defect,
            affected_cells: report.affected_cells,
        })
    }
}
