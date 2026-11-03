//! Physical rigid-wall wrench through the original material and VOF owners.

use super::*;
use crate::d3q19::freesurface3::{PlanarRigidWallMotion3, PlanarRigidWallReport3};

/// Rigid material motion along one fixed, fully flooded domain plane.
/// The liquid domain's physical origin is zero. Angular speed is right handed
/// about the positive coordinate axis normal to the selected plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPlanarRigidWallMotion3 {
    /// Unchanged plane; normal motion and tilting are unsupported.
    pub wall: PlanarWall3,
    /// Fixed reference point for local velocity and all moments [m].
    pub origin_m: [f64; 3],
    /// Tangential velocity at that reference point [m/s].
    pub velocity_at_origin_m_s: [f64; 3],
    /// Angular speed about the positive normal coordinate axis [rad/s].
    pub angular_speed_rad_s: f64,
}

impl MaterialPlanarRigidWallMotion3 {
    /// Bound the local rigid speed over a complete finite rectangular plane.
    /// Checks all four corners. This does not bound an infinite spinning wall.
    pub fn maximum_speed_m_s(
        self,
        cells: [usize; 3],
        cell_size_m: f64,
    ) -> Result<f64, MaterialLbmError> {
        positive(
            cell_size_m,
            "rigid wall cell size must be finite and positive",
        )?;
        let geometric = PlanarRigidWallMotion3 {
            wall: self.wall,
            origin: self.origin_m.map(|v| v / cell_size_m),
            velocity_at_origin: self.velocity_at_origin_m_s,
            angular_speed: self.angular_speed_rad_s * cell_size_m,
        };
        Ok(geometric.maximum_speed(cells)?)
    }
}

/// One accepted original free-surface step and its reciprocal physical wrench.
///
/// Raw momentum includes the native numerical pressure datum. Physical support
/// loads are separately relative to, or inclusive of, the explicitly declared
/// reference pressure. All moments use `motion.origin_m`. Mechanical wall work
/// is not caloric heat or a total energy certificate for the liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPlanarRigidWallReport3 {
    /// Prescribed local motion and common moment origin.
    pub motion: MaterialPlanarRigidWallMotion3,
    /// Actual raw liquid and reciprocal wall impulses [kg m/s].
    pub fluid_impulse_kg_m_s: [f64; 3],
    pub wall_impulse_kg_m_s: [f64; 3],
    /// Actual link angular impulse into liquid [kg m²/s].
    pub fluid_angular_impulse_kg_m2_s: [f64; 3],
    /// Exactly reciprocal wall angular impulse [kg m²/s].
    pub wall_angular_impulse_kg_m2_s: [f64; 3],
    /// Numerical density-one pressure impulse into liquid [kg m/s].
    pub reference_pressure_fluid_impulse_kg_m_s: [f64; 3],
    /// Moment of that numerical reference impulse [kg m²/s].
    pub reference_pressure_fluid_angular_impulse_kg_m2_s: [f64; 3],
    /// Wall force about the declared reference [N].
    pub pressure_perturbation_wall_force_n: [f64; 3],
    /// Wall moment about the declared reference [N m].
    pub pressure_perturbation_wall_torque_nm: [f64; 3],
    /// Wall force including an explicitly declared absolute reference [N].
    pub absolute_wall_force_n: Option<[f64; 3]>,
    /// Wall torque including that same absolute reference [N m].
    pub absolute_wall_torque_nm: Option<[f64; 3]>,
    /// Derivative of `(J_a,J_b,L_axis)` into liquid with respect to
    /// `(V_a,V_b,omega)` at this fixed origin and original input fluid.
    /// `a=(axis+1)%3`, `b=(axis+2)%3`. The 2×2 translation block has units kg,
    /// the cross blocks kg m and the spin diagonal kg m². This is measured
    /// in the original reflection, not a substitute for final native replay.
    pub fluid_response: [[f64; 3]; 3],
    /// Mechanical work into liquid and its reciprocal wall work [J].
    pub fluid_work_j: f64,
    pub wall_work_j: f64,
    /// Complete flooded area [m²] and its geometric centroid [m].
    pub wetted_area_m2: f64,
    pub wetted_centroid_m: [f64; 3],
    /// Original total liquid mass around the accepted step [kg].
    pub mass_before_kg: f64,
    pub mass_after_kg: f64,
    /// Absolute relative native liquid mass defect.
    pub relative_mass_defect: f64,
    /// Maximum physical wet-cell speed at both endpoints [m/s].
    pub maximum_speed_m_s: f64,
    /// Maximum local wall speed over the finite plane [m/s].
    pub maximum_wall_speed_m_s: f64,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

impl MaterialFreeSurface3<'_> {
    /// Propose two tangential translations and axial spin on the original liquid.
    /// Local corner, link-relative and endpoint speeds retain the original
    /// native and material scenario limits. All outputs retain the original
    /// source/scaling; failure or cancellation changes no accepted state.
    /// This fixed-plane operator has no cover/vacate, moving contact line,
    /// tilted geometry, finite plate edges or heat deposition.
    pub fn planar_rigid_wall_motion_candidate_si(
        &self,
        motion: MaterialPlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(Self, MaterialPlanarRigidWallReport3), MaterialLbmError> {
        let native_motion = self.prepare_planar_rigid_wall_motion_si(motion, gate)?;
        let (candidate, raw) = self.lattice.planar_rigid_wall_motion_candidate(
            native_motion,
            maximum_relative_mass_defect,
            gate,
        )?;
        let report = self.planar_rigid_wall_report_si(motion, raw, gate)?;
        Ok((
            Self {
                lattice: candidate,
                scaling: self.scaling,
            },
            report,
        ))
    }

    /// Shared source/scenario admission before either enclosing owner's clone.
    pub(crate) fn prepare_planar_rigid_wall_motion_si(
        &self,
        motion: MaterialPlanarRigidWallMotion3,
        gate: &CancelGate,
    ) -> Result<PlanarRigidWallMotion3, MaterialLbmError> {
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        let scale = &self.scaling;
        let native_motion = PlanarRigidWallMotion3 {
            wall: motion.wall,
            origin: motion.origin_m.map(|v| v / scale.cell_size_m),
            velocity_at_origin: motion
                .velocity_at_origin_m_s
                .map(|v| v / scale.velocity_scale),
            angular_speed: motion.angular_speed_rad_s * scale.timestep_s,
        };
        if native_motion.maximum_speed(self.lattice.domain_cells())? > scale.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "local rigid wall speed exceeds the original material scenario",
            ));
        }
        Ok(native_motion)
    }

    /// Read the actual final native receipt through the original SI/pressure owner.
    pub(crate) fn planar_rigid_wall_report_si(
        &self,
        motion: MaterialPlanarRigidWallMotion3,
        raw: PlanarRigidWallReport3,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarRigidWallReport3, MaterialLbmError> {
        let scale = &self.scaling;
        let units = self.planar_wall_units()?;
        let angular_unit = positive(
            units[0] * scale.cell_size_m,
            "rigid wall angular-impulse unit is unrepresentable",
        )?;
        let torque_unit = positive(
            angular_unit / scale.timestep_s,
            "rigid wall torque unit is unrepresentable",
        )?;
        let response_units = [
            scale.mass_unit_kg,
            positive(
                scale.mass_unit_kg * scale.cell_size_m,
                "rigid wall cross-response unit is unrepresentable",
            )?,
            positive(
                (scale.mass_unit_kg * scale.cell_size_m) * scale.cell_size_m,
                "rigid wall spin-response unit is unrepresentable",
            )?,
        ];
        if raw.maximum_speed > scale.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "rigid wall-driven liquid exceeds the original material scenario speed",
            ));
        }
        let linear = self.planar_wall_report_si(
            PlanarWallMotionReport3 {
                wall: motion.wall,
                velocity: raw.motion.velocity_at_origin,
                fluid_impulse: raw.fluid_impulse,
                wall_impulse: raw.wall_impulse,
                reference_pressure_fluid_impulse: raw.reference_pressure_fluid_impulse,
                fluid_work: raw.fluid_work,
                wall_work: raw.wall_work,
                wetted_area: raw.wetted_area,
                mass_before: raw.mass_before,
                mass_after: raw.mass_after,
                relative_mass_defect: raw.relative_mass_defect,
                maximum_speed: raw.maximum_speed,
            },
            motion.velocity_at_origin_m_s,
            units,
        )?;
        let fluid_angular_impulse_kg_m2_s = vector_si(
            raw.fluid_angular_impulse,
            angular_unit,
            "rigid wall SI angular impulse is unrepresentable",
        )?;
        let reference_pressure_fluid_angular_impulse_kg_m2_s = vector_si(
            raw.reference_pressure_fluid_angular_impulse,
            angular_unit,
            "rigid wall SI reference moment is unrepresentable",
        )?;
        let pressure_perturbation_wall_torque_nm = vector_si(
            core::array::from_fn(|d| {
                raw.wall_angular_impulse[d] + raw.reference_pressure_fluid_angular_impulse[d]
            }),
            torque_unit,
            "rigid wall SI perturbation torque is unrepresentable",
        )?;
        let wetted_centroid_m = vector_si(
            raw.wetted_centroid,
            scale.cell_size_m,
            "rigid wall centroid is unrepresentable",
        )?;
        let absolute_wall_torque_nm = scale
            .reference_pressure_pa
            .map(|pressure| {
                let reference_force = positive(
                    pressure * linear.wetted_area_m2,
                    "rigid wall absolute pressure force is unrepresentable",
                )?;
                let force = raw.reference_pressure_fluid_impulse.map(|value| {
                    if value == 0.0 {
                        0.0
                    } else {
                        -value.signum() * reference_force
                    }
                });
                let moment = cross(
                    core::array::from_fn(|d| wetted_centroid_m[d] - motion.origin_m[d]),
                    force,
                );
                vector_si(
                    core::array::from_fn(|d| pressure_perturbation_wall_torque_nm[d] + moment[d]),
                    1.0,
                    "rigid wall absolute torque is unrepresentable",
                )
            })
            .transpose()?;
        let mut fluid_response = [[0.0; 3]; 3];
        for (i, row) in fluid_response.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = finite(
                    raw.fluid_response[i][j]
                        * response_units[usize::from(i == 2) + usize::from(j == 2)],
                    "rigid wall SI response is unrepresentable",
                )?;
            }
        }
        let report = MaterialPlanarRigidWallReport3 {
            motion,
            fluid_impulse_kg_m_s: linear.fluid_impulse_kg_m_s,
            wall_impulse_kg_m_s: linear.wall_impulse_kg_m_s,
            fluid_angular_impulse_kg_m2_s,
            wall_angular_impulse_kg_m2_s: fluid_angular_impulse_kg_m2_s.map(|v| -v),
            reference_pressure_fluid_impulse_kg_m_s: linear.reference_pressure_fluid_impulse_kg_m_s,
            reference_pressure_fluid_angular_impulse_kg_m2_s,
            pressure_perturbation_wall_force_n: linear.pressure_perturbation_wall_force_n,
            pressure_perturbation_wall_torque_nm,
            absolute_wall_force_n: linear.absolute_wall_force_n,
            absolute_wall_torque_nm,
            fluid_response,
            fluid_work_j: linear.fluid_work_j,
            wall_work_j: linear.wall_work_j,
            wetted_area_m2: linear.wetted_area_m2,
            wetted_centroid_m,
            mass_before_kg: linear.mass_before_kg,
            mass_after_kg: linear.mass_after_kg,
            relative_mass_defect: linear.relative_mass_defect,
            maximum_speed_m_s: linear.maximum_speed_m_s,
            maximum_wall_speed_m_s: finite(
                raw.maximum_wall_speed * scale.velocity_scale,
                "rigid wall SI speed is unrepresentable",
            )?,
        };
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        Ok(report)
    }

    /// Atomically publish the same fully admitted original rigid-wall candidate.
    pub fn step_with_planar_rigid_wall_motion_si(
        &mut self,
        motion: MaterialPlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarRigidWallReport3, MaterialLbmError> {
        let (candidate, report) =
            self.planar_rigid_wall_motion_candidate_si(motion, maximum_relative_mass_defect, gate)?;
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        *self = candidate;
        Ok(report)
    }
}
