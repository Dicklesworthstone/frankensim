//! Physical units for uniform dyadic refinement of the original frozen liquid.

use super::{MaterialFreeSurface3, MaterialLbmError, finite, positive};
use crate::d3q19::freesurface3::{
    DyadicRefinementControls3, DyadicRefinementInventory3, DyadicRefinementReport3,
};
use fs_exec::CancelGate;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialDyadicRefinementControls3 {
    pub max_cells: usize,
    pub max_active_tiles: usize,
    pub maximum_speed_m_s: f64,
    pub maximum_mass_defect_kg: f64,
    pub maximum_first_moment_defect_kg_m: f64,
    pub maximum_momentum_defect_kg_m_s: f64,
    pub maximum_angular_momentum_defect_kg_m2_s: f64,
    pub maximum_kinetic_energy_defect_j: f64,
    pub maximum_viscous_stress_defect_pa: f64,
}

/// Located liquid moments about the fixed physical domain's lower corner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialDyadicRefinementInventory3 {
    pub mass_kg: f64,
    pub first_mass_moment_kg_m: [f64; 3],
    pub momentum_kg_m_s: [f64; 3],
    pub angular_momentum_kg_m2_s: [f64; 3],
    pub kinetic_energy_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct MaterialDyadicRefinementReport3 {
    pub source_cells: [usize; 3],
    pub target_cells: [usize; 3],
    pub source_active_tiles: usize,
    pub target_active_tiles: usize,
    pub before: MaterialDyadicRefinementInventory3,
    pub after: MaterialDyadicRefinementInventory3,
    pub absolute_mass_defect_kg: f64,
    pub absolute_first_moment_defect_kg_m: f64,
    pub absolute_momentum_defect_kg_m_s: f64,
    pub absolute_angular_momentum_defect_kg_m2_s: f64,
    pub kinetic_energy_change_j: f64,
    pub absolute_kinetic_energy_defect_j: f64,
    pub maximum_viscous_stress_defect_pa: f64,
    pub source_geometry_epoch: u64,
    pub target_geometry_epoch: u64,
}

fn represented(input: f64, output: f64) -> Result<f64, MaterialLbmError> {
    finite(output, "physical dyadic refinement is unrepresentable")?;
    if input != 0.0 && output == 0.0 {
        return Err(MaterialLbmError::Invalid(
            "physical dyadic refinement underflows",
        ));
    }
    Ok(output)
}

fn times(value: f64, factor: f64) -> Result<f64, MaterialLbmError> {
    represented(value, value * factor)
}

fn scaled_vector(value: [f64; 3], factor: f64) -> Result<[f64; 3], MaterialLbmError> {
    Ok([
        times(value[0], factor)?,
        times(value[1], factor)?,
        times(value[2], factor)?,
    ])
}

impl<'a> MaterialFreeSurface3<'a> {
    /// Refine an evolved frozen-material VOF state on the SAME physical domain.
    ///
    /// Cell size and timestep both halve, using the original material source,
    /// viscosity, pressure datum and velocity/pressure units. Native BGK
    /// relaxation and force/capillary scaling are inherited from the original
    /// accepted lattice. Mass per lattice cell becomes one eighth. The source
    /// is not re-queried or replaced and no physical time advances.
    ///
    /// The underlying native transfer preserves the actual voxel-average
    /// liquid mass, first moment, momenta and kinetic storage within the stated
    /// SI error budgets, and preserves physical viscous stress. It resolves no
    /// additional subcell interface shape. The following finer timesteps have
    /// their own spatial/temporal error; trajectory identity is not implied.
    ///
    /// This entry carries frozen isothermal source state, held/acceleration
    /// force ownership and supported fixed geometry. Evolving thermal/rheology,
    /// source-local bulk/capillarity and dynamic contact-line histories require
    /// their own transfer owners. All native and SI checks, including a final
    /// cancellation poll, precede returning a new owner; the original remains
    /// unchanged on success, refusal and cancellation.
    pub fn refine_dyadic(
        &self,
        controls: MaterialDyadicRefinementControls3,
        gate: &CancelGate,
    ) -> Result<(Self, MaterialDyadicRefinementReport3), MaterialLbmError> {
        crate::impulse::poll(gate)
            .map_err(|_| crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled)?;
        self.scaling.require_frozen_source()?;
        if self.scaling.bulk_viscosity.is_some() {
            return Err(MaterialLbmError::Invalid(
                "dyadic refinement needs the original BGK frozen source",
            ));
        }
        let s = self.scaling;
        let mass = positive(s.mass_unit_kg, "invalid refinement mass unit")?;
        let first = positive(mass * s.cell_size_m, "invalid refinement first-moment unit")?;
        let momentum = positive(mass * s.velocity_scale, "invalid refinement momentum unit")?;
        let angular = positive(
            momentum * s.cell_size_m,
            "invalid refinement angular-momentum unit",
        )?;
        let energy = positive(
            momentum * s.velocity_scale,
            "invalid refinement kinetic-energy unit",
        )?;
        let stress = positive(s.pressure_unit_pa, "invalid refinement stress unit")?;
        let units = [mass, first, momentum, angular, energy, stress];
        let physical = [
            controls.maximum_mass_defect_kg,
            controls.maximum_first_moment_defect_kg_m,
            controls.maximum_momentum_defect_kg_m_s,
            controls.maximum_angular_momentum_defect_kg_m2_s,
            controls.maximum_kinetic_energy_defect_j,
            controls.maximum_viscous_stress_defect_pa,
        ];
        let mut budgets = [0.0; 6];
        for i in 0..6 {
            if !physical[i].is_finite() || physical[i] < 0.0 {
                return Err(MaterialLbmError::Invalid(
                    "physical refinement needs finite nonnegative budgets",
                ));
            }
            budgets[i] = represented(physical[i], physical[i] / units[i])?;
        }
        positive(
            controls.maximum_speed_m_s,
            "physical refinement speed must be positive",
        )?;
        let maximum_speed = represented(
            controls.maximum_speed_m_s,
            controls.maximum_speed_m_s / s.velocity_scale,
        )?;
        if maximum_speed > s.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "refinement speed exceeds the original scenario bound",
            ));
        }
        let mut scaling = s;
        scaling.cell_size_m = positive(
            times(s.cell_size_m, 0.5)?,
            "refined cell width must be positive",
        )?;
        scaling.timestep_s = positive(
            times(s.timestep_s, 0.5)?,
            "refined timestep must be positive",
        )?;
        scaling.mass_unit_kg = positive(
            times(s.mass_unit_kg, 0.125)?,
            "refined mass unit must be positive",
        )?;
        scaling.surface_tension_lattice = times(s.surface_tension_lattice, 2.0)?;
        scaling.flow.viscosity = positive(
            times(s.flow.viscosity, 2.0)?,
            "refined lattice viscosity must be positive",
        )?;
        scaling.flow.tau = finite(
            2.0_f64.mul_add(s.flow.tau, -0.5),
            "refined relaxation is unrepresentable",
        )?;
        scaling.flow.tau_margin = positive(
            scaling.flow.tau - 0.5,
            "refined relaxation has no positive margin",
        )?;
        let (lattice, native) = self.lattice.refine_dyadic(
            DyadicRefinementControls3 {
                max_cells: controls.max_cells,
                max_active_tiles: controls.max_active_tiles,
                maximum_speed,
                maximum_mass_defect: budgets[0],
                maximum_first_moment_defect: budgets[1],
                maximum_momentum_defect: budgets[2],
                maximum_angular_momentum_defect: budgets[3],
                maximum_kinetic_energy_defect: budgets[4],
                maximum_viscous_stress_defect: budgets[5],
            },
            gate,
        )?;
        let report = receipt(native, units, physical)?;
        crate::impulse::poll(gate)
            .map_err(|_| crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled)?;
        Ok((Self { lattice, scaling }, report))
    }
}

fn inventory(
    value: DyadicRefinementInventory3,
    units: [f64; 6],
) -> Result<MaterialDyadicRefinementInventory3, MaterialLbmError> {
    Ok(MaterialDyadicRefinementInventory3 {
        mass_kg: times(value.mass, units[0])?,
        first_mass_moment_kg_m: scaled_vector(value.first_mass_moment, units[1])?,
        momentum_kg_m_s: scaled_vector(value.momentum, units[2])?,
        angular_momentum_kg_m2_s: scaled_vector(value.angular_momentum, units[3])?,
        kinetic_energy_j: times(value.kinetic_energy, units[4])?,
    })
}

fn difference(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).abs()).sum()
}

fn receipt(
    native: DyadicRefinementReport3,
    units: [f64; 6],
    budgets: [f64; 6],
) -> Result<MaterialDyadicRefinementReport3, MaterialLbmError> {
    let before = inventory(native.before, units)?;
    let after = inventory(native.after, units)?;
    // Retain the stable local native differences AND admit independently
    // rounded SI endpoint totals; conversion must not hide a physical defect.
    let extra = [
        (after.mass_kg - before.mass_kg).abs(),
        difference(after.first_mass_moment_kg_m, before.first_mass_moment_kg_m),
        difference(after.momentum_kg_m_s, before.momentum_kg_m_s),
        difference(
            after.angular_momentum_kg_m2_s,
            before.angular_momentum_kg_m2_s,
        ),
        (after.kinetic_energy_j - before.kinetic_energy_j).abs(),
        0.0,
    ];
    let values = [
        native.absolute_mass_defect,
        native.absolute_first_moment_defect,
        native.absolute_momentum_defect,
        native.absolute_angular_momentum_defect,
        native.absolute_kinetic_energy_defect,
        native.maximum_viscous_stress_defect,
    ];
    let mut defects = [0.0; 6];
    for i in 0..6 {
        defects[i] = finite(
            times(values[i], units[i])? + extra[i],
            "physical refinement readback is unrepresentable",
        )?;
        if defects[i] > budgets[i] {
            return Err(MaterialLbmError::Invalid(
                "physical refinement exceeds its SI error budget",
            ));
        }
    }
    Ok(MaterialDyadicRefinementReport3 {
        source_cells: native.source_cells,
        target_cells: native.target_cells,
        source_active_tiles: native.source_active_tiles,
        target_active_tiles: native.target_active_tiles,
        before,
        after,
        absolute_mass_defect_kg: defects[0],
        absolute_first_moment_defect_kg_m: defects[1],
        absolute_momentum_defect_kg_m_s: defects[2],
        absolute_angular_momentum_defect_kg_m2_s: defects[3],
        kinetic_energy_change_j: times(native.kinetic_energy_change, units[4])?,
        absolute_kinetic_energy_defect_j: defects[4],
        maximum_viscous_stress_defect_pa: defects[5],
        source_geometry_epoch: native.source_geometry_epoch,
        target_geometry_epoch: native.target_geometry_epoch,
    })
}

#[cfg(test)]
mod tests;
