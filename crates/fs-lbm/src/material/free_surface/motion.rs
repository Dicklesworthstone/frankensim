//! SI wall motion and reaction through the original frozen liquid trajectory.

use super::*;
use crate::d3q19::freesurface3::{FreeSurfaceError3, PlanarWallMotionReport3};

/// One accepted isothermal wall-driven free-surface timestep in SI units.
///
/// The paired raw impulses retain the numerical LBM pressure datum for momentum
/// accounting. Physical wall loads are reported separately about the declared
/// reference; an absolute load is available only with with_pressure_reference.
/// Wall work is tangential and therefore independent of the pressure datum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPlanarWallMotionReport3 {
    /// Domain plane whose material moves tangentially.
    pub wall: PlanarWall3,
    /// Prescribed velocity in the existing liquid inertial frame [m/s].
    pub velocity_m_s: [f64; 3],
    /// Actual native population impulse into liquid [kg m/s], including datum.
    pub fluid_impulse_kg_m_s: [f64; 3],
    /// Exactly reciprocal raw impulse on the wall [kg m/s], including datum.
    pub wall_impulse_kg_m_s: [f64; 3],
    /// Numerical density-one pressure impulse into liquid [kg m/s].
    pub reference_pressure_fluid_impulse_kg_m_s: [f64; 3],
    /// Wall force about the declared reference pressure [N].
    pub pressure_perturbation_wall_force_n: [f64; 3],
    /// Wall force including the explicitly declared absolute reference [N].
    /// None means no absolute pressure was supplied; no pressure is invented.
    pub absolute_wall_force_n: Option<[f64; 3]>,
    /// Mechanical energy delivered by the prescribed wall to liquid [J].
    pub fluid_work_j: f64,
    /// Equal and opposite energy received by the external wall actuator [J].
    pub wall_work_j: f64,
    /// Complete flooded plane area [m²].
    pub wetted_area_m2: f64,
    /// Native total liquid mass before and after the accepted step [kg].
    pub mass_before_kg: f64,
    /// Native total liquid mass after the accepted step [kg].
    pub mass_after_kg: f64,
    /// Absolute relative native liquid mass defect.
    pub relative_mass_defect: f64,
    /// Maximum physical liquid speed at the two admitted endpoints [m/s].
    pub maximum_speed_m_s: f64,
}

fn vector_si(
    value: [f64; 3],
    unit: f64,
    reason: &'static str,
) -> Result<[f64; 3], MaterialLbmError> {
    let mut result = [0.0; 3];
    for (out, value) in result.iter_mut().zip(value) {
        *out = finite(value * unit, reason)?;
    }
    Ok(result)
}

impl MaterialFreeSurface3<'_> {
    /// Advance the original liquid with one fully flooded moving planar wall.
    ///
    /// Velocity is prescribed tangentially in the same inertial frame as
    /// sample_at_si, constant over the original physical timestep. A moving
    /// bottom below a free surface is supported. Normal motion, cell
    /// cover/vacate, gas or partial fill on the selected plane and a moving
    /// contact line refuse. The caller owns the substrate/actuator dynamics.
    ///
    /// Original source density, viscosity, capillarity, force and native VOF
    /// transport remain authoritative. Actual link reflection supplies
    /// reciprocal impulse and mechanical wall work; no second fluid integrator
    /// or thermal dissipation model is introduced. This frozen isothermal
    /// method does not deposit wall work as caloric heat.
    ///
    /// The same private native candidate is committed only after source-scale,
    /// endpoint speed, SI representability, mass and cancellation checks.
    ///
    /// # Errors
    /// Refuses wall or liquid speed above the original scenario bound,
    /// unsupported native wall state, excessive mass defect, native failure,
    /// unrepresentable output, or observed cancellation without publication.
    pub fn step_with_planar_wall_motion_si(
        &mut self,
        wall: PlanarWall3,
        velocity_m_s: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarWallMotionReport3, MaterialLbmError> {
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        let scale = &self.scaling;
        let velocity = velocity_m_s.map(|v| v / scale.velocity_scale);
        let wall_speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
        if !wall_speed.is_finite() || wall_speed > scale.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "moving-wall velocity exceeds the original material scenario speed",
            ));
        }
        let impulse_unit = positive(
            scale.mass_unit_kg * scale.velocity_scale,
            "moving-wall SI impulse unit is unrepresentable",
        )?;
        let force_unit = positive(
            impulse_unit / scale.timestep_s,
            "moving-wall SI force unit is unrepresentable",
        )?;
        let work_unit = positive(
            impulse_unit * scale.velocity_scale,
            "moving-wall SI work unit is unrepresentable",
        )?;
        let (candidate, raw): (_, PlanarWallMotionReport3) = self
            .lattice
            .planar_wall_motion_candidate(wall, velocity, maximum_relative_mass_defect, gate)?;
        if raw.maximum_speed > scale.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "wall-driven liquid exceeds the original material scenario speed",
            ));
        }
        let fluid_impulse_kg_m_s = vector_si(
            raw.fluid_impulse,
            impulse_unit,
            "moving-wall SI fluid impulse is unrepresentable",
        )?;
        let wall_impulse_kg_m_s = fluid_impulse_kg_m_s.map(|v| -v);
        let reference_pressure_fluid_impulse_kg_m_s = vector_si(
            raw.reference_pressure_fluid_impulse,
            impulse_unit,
            "moving-wall SI reference impulse is unrepresentable",
        )?;
        let pressure_perturbation_wall_force_n = vector_si(
            core::array::from_fn(|d| raw.wall_impulse[d] + raw.reference_pressure_fluid_impulse[d]),
            force_unit,
            "moving-wall SI pressure-perturbation force is unrepresentable",
        )?;
        let wetted_area_m2 = positive(
            (raw.wetted_area * scale.cell_size_m) * scale.cell_size_m,
            "moving-wall SI wetted area is unrepresentable",
        )?;
        let absolute_wall_force_n = scale
            .reference_pressure_pa
            .map(|pressure| {
                let reference_force = positive(
                    pressure * wetted_area_m2,
                    "moving-wall absolute reference force is unrepresentable",
                )?;
                let normal = raw.reference_pressure_fluid_impulse.map(|v| v.signum());
                let mut absolute = pressure_perturbation_wall_force_n;
                for d in 0..3 {
                    // signum(0.0) is +1, so zero reference components must remain zero.
                    if raw.reference_pressure_fluid_impulse[d] != 0.0 {
                        absolute[d] = finite(
                            absolute[d] - normal[d] * reference_force,
                            "moving-wall absolute force is unrepresentable",
                        )?;
                    }
                }
                Ok::<_, MaterialLbmError>(absolute)
            })
            .transpose()?;
        let fluid_work_j = finite(
            raw.fluid_work * work_unit,
            "moving-wall SI actuator work is unrepresentable",
        )?;
        let report = MaterialPlanarWallMotionReport3 {
            wall,
            velocity_m_s,
            fluid_impulse_kg_m_s,
            wall_impulse_kg_m_s,
            reference_pressure_fluid_impulse_kg_m_s,
            pressure_perturbation_wall_force_n,
            absolute_wall_force_n,
            fluid_work_j,
            wall_work_j: -fluid_work_j,
            wetted_area_m2,
            mass_before_kg: positive(
                raw.mass_before * scale.mass_unit_kg,
                "moving-wall SI initial mass is unrepresentable",
            )?,
            mass_after_kg: positive(
                raw.mass_after * scale.mass_unit_kg,
                "moving-wall SI final mass is unrepresentable",
            )?,
            relative_mass_defect: raw.relative_mass_defect,
            maximum_speed_m_s: finite(
                raw.maximum_speed * scale.velocity_scale,
                "moving-wall SI liquid speed is unrepresentable",
            )?,
        };
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        self.lattice = candidate;
        Ok(report)
    }
}
