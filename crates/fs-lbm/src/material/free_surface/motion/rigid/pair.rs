//! Two original rigid reflections sharing one source-bound liquid timestep.

use super::*;
use crate::d3q19::freesurface3::PlanarRigidWallPairReport3;
use fs_math::dd::Dd;

/// Both opposing flooded wall exchanges from one accepted native VOF tick.
///
/// Each wall retains its own velocity and moment origin in caller order. The
/// mass and speed fields describe the common liquid once; the corresponding
/// per-wall fields repeat that same whole-domain observation and must not be
/// added. Raw impulse includes the numerical pressure datum. Physical support
/// loads remain in each wall's original pressure-reference account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPlanarRigidWallPairReport3 {
    /// Original physical wrench, affine response and pressure account per wall.
    pub walls: [MaterialPlanarRigidWallReport3; 2],
    /// Total raw impulse delivered by these two walls to liquid [kg m/s].
    pub fluid_impulse_kg_m_s: [f64; 3],
    /// Exactly reciprocal total impulse on both walls [kg m/s].
    pub wall_impulse_kg_m_s: [f64; 3],
    /// Sum of actual link velocity dotted with liquid impulse [J].
    pub fluid_work_j: f64,
    /// Exactly reciprocal mechanical work on both walls [J].
    pub wall_work_j: f64,
    /// Common whole-domain liquid mass before this single tick [kg].
    pub mass_before_kg: f64,
    /// Common whole-domain liquid mass after this single tick [kg].
    pub mass_after_kg: f64,
    /// Native absolute relative mass defect of the common candidate.
    pub relative_mass_defect: f64,
    /// Maximum wet-cell endpoint speed over the common liquid [m/s].
    pub maximum_speed_m_s: f64,
    /// Maximum local rigid speed over either finite wall plane [m/s].
    pub maximum_wall_speed_m_s: f64,
}

impl MaterialPlanarRigidWallPairReport3 {
    /// Combine both actual raw angular impulses about an explicitly chosen
    /// physical reference point [m], in the liquid domain's inertial frame.
    /// Each original moment is shifted by `(wall_origin - origin) cross J`.
    /// This includes orbital moment and does not add moments about unlike
    /// origins. Double-double accumulation retains opposing pressure terms.
    pub fn fluid_angular_impulse_at_si(
        &self,
        origin_m: [f64; 3],
    ) -> Result<[f64; 3], MaterialLbmError> {
        if !origin_m.into_iter().all(f64::is_finite) {
            return Err(MaterialLbmError::Invalid(
                "wall-pair angular-impulse origin must be finite",
            ));
        }
        let mut result = [Dd::ZERO; 3];
        for wall in &self.walls {
            let arm: [Dd; 3] = core::array::from_fn(|d| {
                Dd::from_f64(wall.motion.origin_m[d]) - Dd::from_f64(origin_m[d])
            });
            let impulse = wall.fluid_impulse_kg_m_s.map(Dd::from_f64);
            for d in 0..3 {
                let a = (d + 1) % 3;
                let b = (d + 2) % 3;
                result[d] = result[d]
                    + Dd::from_f64(wall.fluid_angular_impulse_kg_m2_s[d])
                    + arm[a] * impulse[b]
                    - arm[b] * impulse[a];
            }
        }
        vector_si(
            result.map(Dd::to_f64),
            1.0,
            "wall-pair SI angular impulse is unrepresentable",
        )
    }

    /// Reciprocal raw angular impulse on both walls about the SAME point,
    /// in kg m²/s. Guides and other domain boundaries are separate accounts.
    pub fn wall_angular_impulse_at_si(
        &self,
        origin_m: [f64; 3],
    ) -> Result<[f64; 3], MaterialLbmError> {
        Ok(self.fluid_angular_impulse_at_si(origin_m)?.map(|v| -v))
    }
}

impl MaterialFreeSurface3<'_> {
    /// Propose both opposing rigid wall motions in ONE original liquid tick.
    ///
    /// The walls must be distinct opposite faces of one coordinate axis; caller
    /// order is retained. Both complete planes remain flooded and fixed, with
    /// independent tangent translations and normal-axis spin. Each native link
    /// has one wall owner, including edge links. Adjacent planes, normal motion,
    /// tilt, finite edges and changing wall coverage remain unsupported.
    ///
    /// One candidate runs the original collision, both reflections and native
    /// VOF conversion. Each wall's 3-by-3 response describes that same incoming
    /// fluid. The one-tick cross-wall response blocks are zero because original
    /// collision consumes accepted populations; later liquid evolution carries
    /// interaction between the walls. This is not two sequential liquid steps.
    ///
    /// Original material/scaling, pressure conventions, body forcing and mass
    /// admission remain authoritative. Both local wall speeds and the common
    /// final liquid speed obey the original material scenario. Signed work is
    /// mechanical wall exchange and is not deposited as caloric heat. Failure
    /// or cancellation returns no partially admitted liquid or wall report.
    pub fn planar_rigid_wall_pair_candidate_si(
        &self,
        motions: [MaterialPlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(Self, MaterialPlanarRigidWallPairReport3), MaterialLbmError> {
        let native = [
            self.prepare_planar_rigid_wall_motion_si(motions[0], gate)?,
            self.prepare_planar_rigid_wall_motion_si(motions[1], gate)?,
        ];
        let (candidate, raw) = self.lattice.planar_rigid_wall_pair_candidate(
            native,
            maximum_relative_mass_defect,
            gate,
        )?;
        let report = self.planar_rigid_wall_pair_report_si(motions, raw, gate)?;
        Ok((
            Self {
                lattice: candidate,
                scaling: self.scaling,
            },
            report,
        ))
    }

    /// Shared original units, pressure reference and endpoint-speed admission
    /// for the cold pair and the complete thermal owner's final native pair.
    /// The measured boundary work stays distinct from any caloric source.
    pub(crate) fn planar_rigid_wall_pair_report_si(
        &self,
        motions: [MaterialPlanarRigidWallMotion3; 2],
        raw: PlanarRigidWallPairReport3,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarRigidWallPairReport3, MaterialLbmError> {
        let units = self.planar_wall_units()?;
        let walls = [
            self.planar_rigid_wall_report_si(motions[0], raw.walls[0], gate)?,
            self.planar_rigid_wall_report_si(motions[1], raw.walls[1], gate)?,
        ];
        let fluid_impulse_kg_m_s = vector_si(
            raw.fluid_impulse,
            units[0],
            "wall-pair SI total impulse is unrepresentable",
        )?;
        let fluid_work_j = finite(
            raw.fluid_work * units[2],
            "wall-pair SI total work is unrepresentable",
        )?;
        let report = MaterialPlanarRigidWallPairReport3 {
            walls,
            fluid_impulse_kg_m_s,
            wall_impulse_kg_m_s: fluid_impulse_kg_m_s.map(|v| -v),
            fluid_work_j,
            wall_work_j: -fluid_work_j,
            mass_before_kg: finite(
                raw.mass_before * self.scaling.mass_unit_kg,
                "wall-pair SI initial mass is unrepresentable",
            )?,
            mass_after_kg: finite(
                raw.mass_after * self.scaling.mass_unit_kg,
                "wall-pair SI final mass is unrepresentable",
            )?,
            relative_mass_defect: raw.relative_mass_defect,
            maximum_speed_m_s: finite(
                raw.maximum_speed * self.scaling.velocity_scale,
                "wall-pair SI liquid speed is unrepresentable",
            )?,
            maximum_wall_speed_m_s: finite(
                raw.maximum_wall_speed * self.scaling.velocity_scale,
                "wall-pair SI local wall speed is unrepresentable",
            )?,
        };
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        Ok(report)
    }

    /// Publish the single fully admitted liquid candidate and both physical
    /// wall exchanges together, after the final cancellation checkpoint.
    pub fn step_with_planar_rigid_wall_pair_si(
        &mut self,
        motions: [MaterialPlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarRigidWallPairReport3, MaterialLbmError> {
        let (candidate, report) =
            self.planar_rigid_wall_pair_candidate_si(motions, maximum_relative_mass_defect, gate)?;
        if gate.is_requested() {
            return Err(FreeSurfaceError3::Cancelled.into());
        }
        *self = candidate;
        Ok(report)
    }
}
