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
use crate::{MassProperties, Vec3, rocket::RocketState};

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
fn multiply_diagonal(inertia: Vec3, vector: Vec3) -> Vec3 {
    Vec3::new(
        inertia.x * vector.x,
        inertia.y * vector.y,
        inertia.z * vector.z,
    )
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
        let omega_body =
            properties.angular_velocity_body_checked(self.angular_momentum_body_kg_m2_s)?;
        let mut states = [*self; 2];
        let mut relative_velocities = [Vec3::ZERO; 2];
        let mut separation_work = 0.0;
        for i in 0..2 {
            let sign = if i == 0 { 1.0 } else { -1.0 };
            let impulse = impulse_on_retained_body_n_s.scale(sign);
            let piece = pieces[i];
            let part = parts[i];
            let arm = application_point_body_m.sub(piece.center_of_mass_body_m);
            let torque_impulse = arm.cross(impulse);
            let spin =
                multiply_diagonal(part.principal_inertia_body(), omega_body).add(torque_impulse);
            let velocity_body = omega_body
                .cross(piece.center_of_mass_body_m)
                .add(impulse.scale(part.mass().recip()));
            relative_velocities[i] = self.orientation.rotate_body_to_world(velocity_body);
            let translation = RocketState {
                position_m: self.translation.position_m.add(
                    self.orientation
                        .rotate_body_to_world(piece.center_of_mass_body_m),
                ),
                velocity_m_s: self.translation.velocity_m_s.add(relative_velocities[i]),
                dry_mass_kg: piece.dry_mass_kg,
                propellant_mass_kg: piece.propellant_mass_kg,
            };
            states[i] = Self::new(translation, self.orientation, spin, piece.inertia)?;
            separation_work += 0.5 * impulse.dot(impulse.scale(part.mass().recip()))
                + 0.5 * torque_impulse.dot(part.angular_velocity_body_checked(torque_impulse)?);
        }
        let mut linear = Vec3::ZERO;
        let mut angular = Vec3::ZERO;
        let mut kinetic = 0.0;
        for i in 0..2 {
            let state = states[i];
            let part = parts[i];
            let relative_momentum = relative_velocities[i].scale(part.mass());
            linear = linear.add(state.translation.velocity_m_s.scale(part.mass()));
            angular = angular
                .add(
                    self.orientation
                        .rotate_body_to_world(state.angular_momentum_body_kg_m2_s),
                )
                .add(
                    self.orientation
                        .rotate_body_to_world(pieces[i].center_of_mass_body_m)
                        .cross(relative_momentum),
                );
            kinetic += 0.5 * relative_momentum.dot(relative_velocities[i])
                + 0.5
                    * state.angular_momentum_body_kg_m2_s.dot(
                        part.angular_velocity_body_checked(state.angular_momentum_body_kg_m2_s)?,
                    );
        }
        let linear_residual = linear.sub(self.translation.velocity_m_s.scale(properties.mass()));
        let angular_residual = angular.sub(
            self.orientation
                .rotate_body_to_world(self.angular_momentum_body_kg_m2_s),
        );
        let energy_residual =
            kinetic - 0.5 * self.angular_momentum_body_kg_m2_s.dot(omega_body) - separation_work;
        if !linear_residual.is_finite()
            || !angular_residual.is_finite()
            || !separation_work.is_finite()
            || !energy_residual.is_finite()
        {
            return Err(RigidRocketError::Invalid(
                "stage momenta or energy are not representable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(StageSeparation {
            retained: states[0],
            released: states[1],
            mass_residual_kg: total_mass - properties.mass(),
            linear_momentum_residual_world_kg_m_s: linear_residual,
            angular_momentum_residual_world_kg_m2_s: angular_residual,
            separation_work_j: separation_work,
            energy_balance_residual_j: energy_residual,
        })
    }
}
