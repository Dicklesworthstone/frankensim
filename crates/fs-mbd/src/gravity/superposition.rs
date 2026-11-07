//! Bounded superposition of the original stationary conservative field owners.

use super::{
    GravityError, GravityLoad, PointGravityField, PointGravitySample, RigidGravityField, finite,
    same_frame,
};
use crate::{MassProperties, Pose, Vec3, Wrench};
use fs_ga::Mat3;
use fs_math::dd::Dd;

/// Multiple stationary sources in one explicitly identified inertial frame.
///
/// Newtonian potentials add, together with their forces, torques and spatial
/// derivatives. This borrowed compositor preserves each source's own model,
/// datum and domain admission. It implements [`PointGravityField`] and/or
/// [`RigidGravityField`] according to the borrowed field type. A slice of
/// `&dyn PointGravityField` or `&dyn RigidGravityField` admits heterogeneous
/// sources without allocating or introducing another mass/integration owner.
///
/// Every query performs O(source count) work in the supplied deterministic order,
/// with fixed-size double-double accumulators and no heap allocation. The
/// caller explicitly caps that count at construction; cancellation remains with
/// the existing consumer around each bounded field query. Nested compositors
/// must budget their own children. Empty superpositions are zero fields, but
/// still validate the requested frame, coordinates and native mass properties.
///
/// Potential increments sum each source's direct increment before rounding the
/// total. Only a source returning `None` is sampled again at both endpoints for
/// its fallback difference. Thus an older custom source does not discard another
/// source's accurate work. Callbacks must be repeatable and side-effect-free;
/// fallback accuracy remains limited by that source's endpoint potentials.
///
/// Sources remain prescribed and stationary: this is not a recoiling N-body
/// system, an interior-gravity model, or an admission of the intervening path.
#[derive(Debug)]
pub struct SuperposedGravity<'a, F: ?Sized> {
    frame: u64,
    fields: &'a [&'a F],
}

impl<'a, F: ?Sized> SuperposedGravity<'a, F> {
    /// Bind the source list and its explicit maximum number of field evaluations.
    ///
    /// The cap counts children, not nested descendants. An increment query may
    /// additionally evaluate two endpoint loads per child without a direct
    /// increment. Frame and physical domains are checked at the actual query.
    pub fn new(
        frame: u64,
        fields: &'a [&'a F],
        maximum_fields: usize,
    ) -> Result<Self, GravityError> {
        if fields.len() > maximum_fields {
            return Err(GravityError::Invalid("gravity source budget exceeded"));
        }
        Ok(Self { frame, fields })
    }

    /// Number of directly borrowed sources, including nested compositors.
    pub const fn source_count(&self) -> usize {
        self.fields.len()
    }

    fn admit_position(&self, frame: u64, position: Vec3) -> Result<(), GravityError> {
        same_frame(frame, self.frame)?;
        if !position.is_finite() {
            return Err(GravityError::Invalid("nonfinite gravity position"));
        }
        Ok(())
    }
}

impl<F: PointGravityField + ?Sized> PointGravityField for SuperposedGravity<'_, F> {
    fn sample_point(
        &self,
        frame: u64,
        position_world_m: Vec3,
    ) -> Result<PointGravitySample, GravityError> {
        self.admit_position(frame, position_world_m)?;
        let mut acceleration = [Dd::ZERO; 3];
        let mut potential = Dd::ZERO;
        let mut gradient = [Dd::ZERO; 9];
        for field in self.fields {
            let sample = field.sample_point(frame, position_world_m)?.validate()?;
            let g = sample.acceleration_world_m_s2;
            accumulate(&mut acceleration, [g.x, g.y, g.z]);
            potential = potential + Dd::from_f64(sample.specific_potential_j_kg);
            accumulate(&mut gradient, sample.acceleration_gradient_s2.m);
        }
        PointGravitySample {
            acceleration_world_m_s2: vector(acceleration),
            specific_potential_j_kg: potential.to_f64(),
            acceleration_gradient_s2: Mat3 {
                m: gradient.map(Dd::to_f64),
            },
        }
        .validate()
    }

    fn specific_potential_increment(
        &self,
        frame: u64,
        before_world_m: Vec3,
        after_world_m: Vec3,
    ) -> Result<Option<f64>, GravityError> {
        self.admit_position(frame, before_world_m)?;
        self.admit_position(frame, after_world_m)?;
        let mut change = Dd::ZERO;
        for field in self.fields {
            let increment =
                match field.specific_potential_increment(frame, before_world_m, after_world_m)? {
                    Some(value) => Dd::from_f64(finite(value)?),
                    None => {
                        let before = field.sample_point(frame, before_world_m)?.validate()?;
                        let after = field.sample_point(frame, after_world_m)?.validate()?;
                        Dd::from_f64(after.specific_potential_j_kg)
                            - Dd::from_f64(before.specific_potential_j_kg)
                    }
                };
            change = change + increment;
        }
        Ok(Some(finite(change.to_f64())?))
    }
}

impl<F: RigidGravityField + ?Sized> RigidGravityField for SuperposedGravity<'_, F> {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        self.admit_position(frame, pose.position_world())?;
        properties.validate()?;
        let mut force = [Dd::ZERO; 3];
        let mut torque = [Dd::ZERO; 3];
        let mut potential = Dd::ZERO;
        for field in self.fields {
            let load = field.load(frame, pose, properties)?.validate()?;
            let f = load.wrench.force_world;
            let t = load.wrench.torque_body;
            accumulate(&mut force, [f.x, f.y, f.z]);
            accumulate(&mut torque, [t.x, t.y, t.z]);
            potential = potential + Dd::from_f64(load.potential_energy_j);
        }
        GravityLoad {
            wrench: Wrench {
                force_world: vector(force),
                torque_body: vector(torque),
            },
            potential_energy_j: potential.to_f64(),
        }
        .validate()
    }

    fn potential_energy_increment(
        &self,
        frame: u64,
        before: Pose,
        after: Pose,
        properties: MassProperties,
    ) -> Result<Option<f64>, GravityError> {
        self.admit_position(frame, before.position_world())?;
        self.admit_position(frame, after.position_world())?;
        properties.validate()?;
        let mut change = Dd::ZERO;
        for field in self.fields {
            let increment = match field
                .potential_energy_increment(frame, before, after, properties)?
            {
                Some(value) => Dd::from_f64(finite(value)?),
                None => {
                    let before = field.load(frame, before, properties)?.validate()?;
                    let after = field.load(frame, after, properties)?.validate()?;
                    Dd::from_f64(after.potential_energy_j) - Dd::from_f64(before.potential_energy_j)
                }
            };
            change = change + increment;
        }
        Ok(Some(finite(change.to_f64())?))
    }
}

fn accumulate<const N: usize>(sum: &mut [Dd; N], values: [f64; N]) {
    for (total, value) in sum.iter_mut().zip(values) {
        *total = *total + Dd::from_f64(value);
    }
}

fn vector(sum: [Dd; 3]) -> Vec3 {
    let [x, y, z] = sum.map(Dd::to_f64);
    Vec3::new(x, y, z)
}
