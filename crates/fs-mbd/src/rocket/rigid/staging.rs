//! Instantaneous separation of two native rocket bodies.
//!
//! Both pieces share the parent principal axes at separation. Their supplied
//! masses, centres and full parallel-axis assembly must reproduce the parent.
//! Equal and opposite impulses act at one common physical point. Each piece
//! inherits v_com + omega cross r before recoil; neither spin, orbital angular
//! momentum nor a discarded stage disappears from the returned result.
//! This event models a prescribed separation impulse, not pyrotechnic pressure,
//! contact/reimpact, fracture or plume interaction.

use super::{RigidRocketError, RigidRocketInertia, RigidRocketState, poll};
use crate::assembly::{
    AssemblyError, AssemblyPart, ImpulseTotals, SeparationImpulse, release_accounting,
    release_components_into,
};
use crate::{MassProperties, Pose, UnitQuaternion, Vec3, rocket::RocketState};

/// One piece, with axes aligned to the parent body and its own COM origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketStage {
    /// Dry mass of this piece [kg].
    pub dry_mass_kg: f64,
    /// Remaining propellant retained by this piece [kg].
    pub propellant_mass_kg: f64,
    /// Piece COM relative to the parent COM, in common body axes [m].
    pub center_of_mass_body_m: Vec3,
    /// Its own COM-referenced dry and propellant inertia law.
    pub inertia: RigidRocketInertia,
}

/// Both surviving owners and actual event accounting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StageSeparation {
    /// Continuing vehicle after recoil and transfer to its new COM.
    pub retained: RigidRocketState,
    /// Released vehicle after the opposite recoil; it can continue independently.
    pub released: RigidRocketState,
    /// Actual piece mass sum minus parent mass [kg].
    pub mass_residual_kg: f64,
    /// Actual sum of world linear momenta minus the parent [kg m/s].
    pub linear_momentum_residual_world_kg_m_s: Vec3,
    /// Spin plus orbital angular momentum residual about the parent COM [kg m²/s].
    pub angular_momentum_residual_world_kg_m2_s: Vec3,
    /// Kinetic energy supplied by the prescribed equal/opposite impulse [J].
    pub separation_work_j: f64,
    /// Relative kinetic-energy change minus separation work [J].
    /// Evaluated in the inertial frame translating with the original parent COM.
    pub energy_balance_residual_j: f64,
}

impl RocketStage {
    fn properties(self) -> Result<MassProperties, RigidRocketError> {
        let translation = RocketState {
            position_m: Vec3::ZERO,
            velocity_m_s: Vec3::ZERO,
            dry_mass_kg: self.dry_mass_kg,
            propellant_mass_kg: self.propellant_mass_kg,
        };
        translation.advance(super::zero_forcing(), 0.0, || false)?;
        if !self.center_of_mass_body_m.is_finite() {
            return Err(RigidRocketError::Invalid("finite stage COM required"));
        }
        self.inertia.mass_properties(translation)
    }
}
fn magnitude(vector: Vec3) -> f64 {
    vector.x.abs().max(vector.y.abs()).max(vector.z.abs())
}

impl RigidRocketState {
    /// Separate two supplied mass distributions at a shared impulse point.
    ///
    /// impulse_on_retained_body_n_s is applied to the retained piece and its
    /// negative to the released piece, at application_point_body_m measured
    /// from the original parent COM. These are body-frame input coordinates.
    /// Checks dry/propellant mass sums, COM and all six entries of the
    /// parallel-axis inertia assembly before publishing either piece.
    /// relative_assembly_tolerance must be positive and at most 1e-6; accepted
    /// discrepancies are reported, never repaired by changing a supplied mass.
    pub fn separate(
        &self,
        retained: RocketStage,
        released: RocketStage,
        impulse_on_retained_body_n_s: Vec3,
        application_point_body_m: Vec3,
        relative_assembly_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StageSeparation, RigidRocketError> {
        poll(&mut cancelled)?;
        let tolerance = relative_assembly_tolerance;
        if !tolerance.is_finite()
            || tolerance <= 0.0
            || tolerance > 1e-6
            || !impulse_on_retained_body_n_s.is_finite()
            || !application_point_body_m.is_finite()
        {
            return Err(RigidRocketError::Invalid(
                "finite staging inputs and assembly tolerance in (0, 1e-6] required",
            ));
        }
        let properties = self.mass_properties()?;
        let pieces = [retained, released];
        let parts = [retained.properties()?, released.properties()?];
        let total_mass = parts[0].mass() + parts[1].mass();
        let dry_residual =
            retained.dry_mass_kg + released.dry_mass_kg - self.translation.dry_mass_kg;
        let propellant_residual = retained.propellant_mass_kg + released.propellant_mass_kg
            - self.translation.propellant_mass_kg;
        if !total_mass.is_finite()
            || dry_residual.abs() > tolerance * self.translation.dry_mass_kg
            || propellant_residual.abs() > tolerance * self.translation.propellant_mass_kg
        {
            return Err(RigidRocketError::Invalid(
                "stages must partition both dry and propellant masses",
            ));
        }
        let first_moment = retained
            .center_of_mass_body_m
            .scale(parts[0].mass())
            .add(released.center_of_mass_body_m.scale(parts[1].mass()));
        let radius = magnitude(retained.center_of_mass_body_m)
            .max(magnitude(released.center_of_mass_body_m));
        if !first_moment.is_finite() || magnitude(first_moment) > (tolerance * total_mass) * radius
        {
            return Err(RigidRocketError::Invalid(
                "stage centres must reproduce the parent COM",
            ));
        }
        let mut assembled = [[0.0; 3]; 3];
        for (piece, part) in pieces.iter().zip(parts) {
            let r = piece.center_of_mass_body_m;
            let d = part.principal_inertia_body();
            let m = part.mass();
            let values = [r.x, r.y, r.z];
            let diagonal = [d.x, d.y, d.z];
            for row in 0..3 {
                for col in 0..3 {
                    assembled[row][col] += if row == col {
                        let a = values[(row + 1) % 3];
                        let b = values[(row + 2) % 3];
                        diagonal[row] + m * (a * a + b * b)
                    } else {
                        -m * values[row] * values[col]
                    };
                }
            }
        }
        let parent_inertia = properties.principal_inertia_body();
        let diagonal = [parent_inertia.x, parent_inertia.y, parent_inertia.z];
        let scale = magnitude(parent_inertia);
        for row in 0..3 {
            for col in 0..3 {
                let expected = if row == col { diagonal[row] } else { 0.0 };
                if !assembled[row][col].is_finite()
                    || (assembled[row][col] - expected).abs() > tolerance * scale
                {
                    return Err(RigidRocketError::Invalid(
                        "stage inertia must reproduce all parent parallel-axis entries",
                    ));
                }
            }
        }
        poll(&mut cancelled)?;
        let components = [
            AssemblyPart {
                properties: parts[0],
                pose_in_reference: Pose::new(
                    retained.center_of_mass_body_m,
                    UnitQuaternion::IDENTITY,
                )?,
            },
            AssemblyPart {
                properties: parts[1],
                pose_in_reference: Pose::new(
                    released.center_of_mass_body_m,
                    UnitQuaternion::IDENTITY,
                )?,
            },
        ];
        let rigid = self.rigid_body()?;
        let mut bodies = [rigid; 2];
        let mut deltas = [ImpulseTotals::default(); 2];
        let release = release_components_into(
            &components,
            properties,
            rigid,
            &[SeparationImpulse {
                part_a: 0,
                part_b: 1,
                point_assembly_m: application_point_body_m,
                impulse_on_a_assembly_ns: impulse_on_retained_body_n_s,
                couple_on_a_assembly_nms: Vec3::ZERO,
            }],
            &mut deltas,
            &mut bodies,
            &mut cancelled,
        )
        .map_err(separation_error)?;
        let mut states = [*self; 2];
        for (i, body) in bodies.iter().enumerate() {
            states[i] = Self::new(
                RocketState {
                    position_m: body.pose().position_world(),
                    velocity_m_s: body.center_of_mass_velocity_world(parts[i])?,
                    dry_mass_kg: pieces[i].dry_mass_kg,
                    propellant_mass_kg: pieces[i].propellant_mass_kg,
                },
                body.pose().orientation(),
                body.angular_momentum_body(),
                pieces[i].inertia,
            )?;
        }
        let actual = [states[0].rigid_body()?, states[1].rigid_body()?];
        let release = release_accounting(
            &components,
            properties,
            rigid,
            &actual,
            release.supplied_work_j,
            &mut cancelled,
        )
        .map_err(separation_error)?;
        poll(&mut cancelled)?;
        Ok(StageSeparation {
            retained: states[0],
            released: states[1],
            mass_residual_kg: total_mass - properties.mass(),
            linear_momentum_residual_world_kg_m_s: release
                .balance
                .linear_momentum_residual_world_ns,
            angular_momentum_residual_world_kg_m2_s: release
                .balance
                .angular_momentum_residual_world_nms,
            separation_work_j: release.supplied_work_j,
            energy_balance_residual_j: release.balance.energy_residual_j,
        })
    }
}

fn separation_error(error: AssemblyError) -> RigidRocketError {
    match error {
        AssemblyError::Dynamics(error) => RigidRocketError::Dynamics(error),
        AssemblyError::Cancelled => RigidRocketError::Cancelled,
        AssemblyError::Invalid(reason) => RigidRocketError::Invalid(reason),
        _ => RigidRocketError::Invalid("native rigid separation refused"),
    }
}
