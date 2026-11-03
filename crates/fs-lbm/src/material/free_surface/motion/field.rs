//! Sourced physical units for the original spatial wall-face operator.

use super::*;
use crate::d3q19::freesurface3::{
    PlanarWallFieldControls3, PlanarWallFieldReport3, PlanarWallVelocityField3,
};

/// Prescribed tangential face velocities along one fixed flooded physical plane.
///
/// The cyclic tangent axes and row-major order are those of
/// [`PlanarWallVelocityField3`]. These are Cartesian face-quadrature values;
/// every reflected link at a face uses the same value. Geometry does not move
/// normally, and a continuous rigid velocity is not resampled at link arms.
#[derive(Clone, Copy, Debug)]
pub struct MaterialPlanarWallVelocityField3<'a> {
    pub wall: PlanarWall3,
    /// Common origin for all angular impulse readbacks [m].
    pub moment_origin_m: [f64; 3],
    /// Tangent a/b velocities [m/s], one pair per Cartesian wall face.
    pub tangential_velocity_m_s: &'a [[f64; 2]],
}

/// Actual reciprocal physical wrench/work, with explicit pressure references.
///
/// Raw impulses close native population momentum. The perturbation and optional
/// absolute loads are physical load interpretations and do not replace that
/// raw account. Moments use the containing report's common origin. Boundary
/// work is mechanical exchange, with no inferred thermal deposit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialWallWrenchExchange3 {
    pub fluid_impulse_kg_m_s: [f64; 3],
    pub wall_impulse_kg_m_s: [f64; 3],
    pub fluid_angular_impulse_kg_m2_s: [f64; 3],
    pub wall_angular_impulse_kg_m2_s: [f64; 3],
    pub reference_pressure_fluid_impulse_kg_m_s: [f64; 3],
    pub reference_pressure_fluid_angular_impulse_kg_m2_s: [f64; 3],
    pub pressure_perturbation_wall_force_n: [f64; 3],
    pub pressure_perturbation_wall_torque_nm: [f64; 3],
    pub absolute_wall_force_n: Option<[f64; 3]>,
    pub absolute_wall_torque_nm: Option<[f64; 3]>,
    pub fluid_work_j: f64,
    pub wall_work_j: f64,
}

/// Resolved force, couple and work for one deforming-wall quadrature face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPlanarWallFaceExchange3 {
    pub liquid_cell: [usize; 3],
    pub center_m: [f64; 3],
    pub velocity_m_s: [f64; 3],
    pub exchange: MaterialWallWrenchExchange3,
    /// Retain this couple as well as force at `center_m` when transferring
    /// the native angular impulse to a solid. Face-center force alone is
    /// insufficient. The piecewise translational face has no angular speed.
    pub fluid_couple_at_face_kg_m2_s: [f64; 3],
    pub wall_couple_at_face_kg_m2_s: [f64; 3],
    /// Diagonal dJ_tangent/du_tangent [kg], from this original input liquid.
    /// Different face/tangent cross terms are zero for this native reflection;
    /// this is not a steady-state or later-timestep fluid response.
    pub fluid_tangential_response_kg: [f64; 2],
    pub reflected_mass_defect_kg: f64,
    pub relative_mass_defect: f64,
}

/// One original material-fluid tick and its distributed mechanical exchange.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialPlanarWallFieldReport3 {
    pub wall: PlanarWall3,
    pub moment_origin_m: [f64; 3],
    /// Same ordering and length as the supplied tangential velocity field.
    pub faces: Vec<MaterialPlanarWallFaceExchange3>,
    pub exchange: MaterialWallWrenchExchange3,
    pub wetted_area_m2: f64,
    pub mass_before_kg: f64,
    pub mass_after_kg: f64,
    pub relative_mass_defect: f64,
    pub maximum_relative_face_mass_defect: f64,
    pub maximum_speed_m_s: f64,
    pub maximum_wall_speed_m_s: f64,
}

// One admitted conversion shared by the global sum and every local wrench.
struct WrenchUnits {
    impulse: f64,
    force: f64,
    work: f64,
    angular_impulse: f64,
    torque: f64,
    // Convert the native cs²-weighted reference, preserving its exact geometry.
    absolute_reference_force: Option<f64>,
    absolute_reference_torque: Option<f64>,
}

impl WrenchUnits {
    fn convert(
        &self,
        fluid_impulse: [f64; 3],
        fluid_angular: [f64; 3],
        reference_impulse: [f64; 3],
        reference_angular: [f64; 3],
        fluid_work: f64,
    ) -> Result<MaterialWallWrenchExchange3, MaterialLbmError> {
        let fluid_impulse_kg_m_s = vector_si(
            fluid_impulse,
            self.impulse,
            "wall field SI impulse is unrepresentable",
        )?;
        let fluid_angular_impulse_kg_m2_s = vector_si(
            fluid_angular,
            self.angular_impulse,
            "wall field SI angular impulse is unrepresentable",
        )?;
        let pressure_perturbation_wall_force_n = vector_si(
            core::array::from_fn(|d| -fluid_impulse[d] + reference_impulse[d]),
            self.force,
            "wall field SI pressure-perturbation force is unrepresentable",
        )?;
        let pressure_perturbation_wall_torque_nm = vector_si(
            core::array::from_fn(|d| -fluid_angular[d] + reference_angular[d]),
            self.torque,
            "wall field SI pressure-perturbation torque is unrepresentable",
        )?;
        let absolute_wall_force_n = self
            .absolute_reference_force
            .map(|unit| {
                vector_si(
                    core::array::from_fn(|d| {
                        pressure_perturbation_wall_force_n[d] - reference_impulse[d] * unit
                    }),
                    1.0,
                    "wall field absolute force is unrepresentable",
                )
            })
            .transpose()?;
        let absolute_wall_torque_nm = self
            .absolute_reference_torque
            .map(|unit| {
                vector_si(
                    core::array::from_fn(|d| {
                        pressure_perturbation_wall_torque_nm[d] - reference_angular[d] * unit
                    }),
                    1.0,
                    "wall field absolute torque is unrepresentable",
                )
            })
            .transpose()?;
        let fluid_work_j = finite(
            fluid_work * self.work,
            "wall field SI work is unrepresentable",
        )?;
        Ok(MaterialWallWrenchExchange3 {
            fluid_impulse_kg_m_s,
            wall_impulse_kg_m_s: fluid_impulse_kg_m_s.map(|v| -v),
            fluid_angular_impulse_kg_m2_s,
            wall_angular_impulse_kg_m2_s: fluid_angular_impulse_kg_m2_s.map(|v| -v),
            reference_pressure_fluid_impulse_kg_m_s: vector_si(
                reference_impulse,
                self.impulse,
                "wall field SI reference impulse is unrepresentable",
            )?,
            reference_pressure_fluid_angular_impulse_kg_m2_s: vector_si(
                reference_angular,
                self.angular_impulse,
                "wall field SI reference moment is unrepresentable",
            )?,
            pressure_perturbation_wall_force_n,
            pressure_perturbation_wall_torque_nm,
            absolute_wall_force_n,
            absolute_wall_torque_nm,
            fluid_work_j,
            wall_work_j: -fluid_work_j,
        })
    }
}

impl MaterialFreeSurface3<'_> {
    /// Propose nonuniform in-plane wall motion on the original material liquid.
    ///
    /// Every face and both native wet-cell endpoints obey the original source
    /// scenario speed bound. The same native reflection supplies locally
    /// mass-neutral momentum and actual link couples; physical unit/pressure
    /// conversion completes before a candidate can escape. Refusal or
    /// cancellation preserves the original population history and source.
    pub fn planar_wall_velocity_field_candidate_si(
        &self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<(Self, MaterialPlanarWallFieldReport3), MaterialLbmError> {
        let (velocities, origin) =
            self.prepare_planar_wall_velocity_field_si(motion, controls, gate)?;
        let (candidate, raw) = self.lattice.planar_wall_velocity_field_candidate(
            PlanarWallVelocityField3 {
                wall: motion.wall,
                moment_origin: origin,
                tangential_velocity: &velocities,
            },
            controls,
            gate,
        )?;
        let report = self.planar_wall_field_report_si(raw, motion.moment_origin_m, gate)?;
        Ok((
            Self {
                lattice: candidate,
                scaling: self.scaling,
            },
            report,
        ))
    }

    // Thermal and isothermal owners admit the same field before a native clone.
    pub(crate) fn prepare_planar_wall_velocity_field_si(
        &self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<(Vec<[f64; 2]>, [f64; 3]), MaterialLbmError> {
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        let count = controls.admit_shape(
            motion.wall,
            self.lattice.domain_cells(),
            motion.tangential_velocity_m_s.len(),
        )?;
        let scale = &self.scaling;
        let mut velocities = Vec::new();
        velocities
            .try_reserve_exact(count)
            .map_err(|_| MaterialLbmError::Invalid("physical wall field allocation refused"))?;
        for (index, pair) in motion.tangential_velocity_m_s.iter().enumerate() {
            if index % 64 == 0 && gate.is_requested() {
                return Err(FreeSurfaceError3::Cancelled.into());
            }
            let native = pair.map(|v| v / scale.velocity_scale);
            let speed = native[0].hypot(native[1]);
            if !speed.is_finite() || speed > scale.flow.u_lattice {
                return Err(MaterialLbmError::Invalid(
                    "local wall field velocity exceeds the original material scenario",
                ));
            }
            velocities.push(native);
        }
        let origin = vector_si(
            motion.moment_origin_m,
            scale.cell_size_m.recip(),
            "wall field lattice moment origin is unrepresentable",
        )?;
        Ok((velocities, origin))
    }

    // Convert the actual native receipt before either enclosing owner accepts.
    pub(crate) fn planar_wall_field_report_si(
        &self,
        raw: PlanarWallFieldReport3,
        moment_origin_m: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MaterialPlanarWallFieldReport3, MaterialLbmError> {
        let scale = &self.scaling;
        let [impulse, force, work] = self.planar_wall_units()?;
        let angular_impulse = positive(
            impulse * scale.cell_size_m,
            "wall field angular impulse unit is unrepresentable",
        )?;
        let area_unit = positive(
            scale.cell_size_m * scale.cell_size_m,
            "wall field physical face area is unrepresentable",
        )?;
        let absolute_reference_force = scale
            .reference_pressure_pa
            .map(|p| {
                positive(
                    (p * scale.cell_size_m) * scale.cell_size_m / CS2,
                    "wall field physical pressure unit is unrepresentable",
                )
            })
            .transpose()?;
        let units = WrenchUnits {
            impulse,
            force,
            work,
            angular_impulse,
            torque: positive(
                angular_impulse / scale.timestep_s,
                "wall field torque unit is unrepresentable",
            )?,
            absolute_reference_force,
            absolute_reference_torque: absolute_reference_force
                .map(|f| {
                    positive(
                        f * scale.cell_size_m,
                        "wall field physical pressure moment unit is unrepresentable",
                    )
                })
                .transpose()?,
        };
        if raw.maximum_speed > scale.flow.u_lattice {
            return Err(MaterialLbmError::Invalid(
                "wall field-driven liquid exceeds the original material scenario speed",
            ));
        }
        let exchange = units.convert(
            raw.fluid_impulse,
            raw.fluid_angular_impulse,
            raw.reference_pressure_fluid_impulse,
            raw.reference_pressure_fluid_angular_impulse,
            raw.fluid_work,
        )?;
        let mut faces = Vec::new();
        faces.try_reserve_exact(raw.faces.len()).map_err(|_| {
            MaterialLbmError::Invalid("physical wall field report allocation refused")
        })?;
        for (index, face) in raw.faces.into_iter().enumerate() {
            if index % 64 == 0 && gate.is_requested() {
                return Err(FreeSurfaceError3::Cancelled.into());
            }
            let fluid_couple_at_face_kg_m2_s = vector_si(
                face.fluid_couple_at_face,
                angular_impulse,
                "wall field SI local couple is unrepresentable",
            )?;
            let fluid_tangential_response_kg = [
                finite(
                    face.fluid_tangential_response[0] * scale.mass_unit_kg,
                    "wall field SI tangent response is unrepresentable",
                )?,
                finite(
                    face.fluid_tangential_response[1] * scale.mass_unit_kg,
                    "wall field SI tangent response is unrepresentable",
                )?,
            ];
            faces.push(MaterialPlanarWallFaceExchange3 {
                liquid_cell: face.liquid_cell,
                center_m: vector_si(
                    face.center,
                    scale.cell_size_m,
                    "wall field SI face center is unrepresentable",
                )?,
                velocity_m_s: vector_si(
                    face.velocity,
                    scale.velocity_scale,
                    "wall field SI velocity is unrepresentable",
                )?,
                exchange: units.convert(
                    face.fluid_impulse,
                    face.fluid_angular_impulse,
                    face.reference_pressure_fluid_impulse,
                    face.reference_pressure_fluid_angular_impulse,
                    face.fluid_work,
                )?,
                fluid_couple_at_face_kg_m2_s,
                wall_couple_at_face_kg_m2_s: fluid_couple_at_face_kg_m2_s.map(|v| -v),
                fluid_tangential_response_kg,
                reflected_mass_defect_kg: finite(
                    face.reflected_mass_defect * scale.mass_unit_kg,
                    "wall field SI local mass defect is unrepresentable",
                )?,
                relative_mass_defect: face.relative_mass_defect,
            });
        }
        let report = MaterialPlanarWallFieldReport3 {
            wall: raw.wall,
            moment_origin_m,
            faces,
            exchange,
            wetted_area_m2: positive(
                raw.wetted_area * area_unit,
                "wall field SI total area is unrepresentable",
            )?,
            mass_before_kg: positive(
                raw.mass_before * scale.mass_unit_kg,
                "wall field SI initial mass is unrepresentable",
            )?,
            mass_after_kg: positive(
                raw.mass_after * scale.mass_unit_kg,
                "wall field SI final mass is unrepresentable",
            )?,
            relative_mass_defect: raw.relative_mass_defect,
            maximum_relative_face_mass_defect: raw.maximum_relative_face_mass_defect,
            maximum_speed_m_s: finite(
                raw.maximum_speed * scale.velocity_scale,
                "wall field SI liquid speed is unrepresentable",
            )?,
            maximum_wall_speed_m_s: finite(
                raw.maximum_wall_speed * scale.velocity_scale,
                "wall field SI maximum speed is unrepresentable",
            )?,
        };
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        Ok(report)
    }

    /// Publish the same admitted native populations and physical field receipt.
    pub fn step_with_planar_wall_velocity_field_si(
        &mut self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarWallFieldReport3, MaterialLbmError> {
        let (candidate, report) =
            self.planar_wall_velocity_field_candidate_si(motion, controls, gate)?;
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        *self = candidate;
        Ok(report)
    }
}
