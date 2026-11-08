//! Exact Newtonian gravity of a prescribed sphere of uniform density.

use super::{PointGravityField, PointGravitySample};
use crate::gravity::{
    CentralGravity, GravityError, GravityLoad, RigidGravityField, finite, same_frame,
};
use crate::{MassProperties, Pose, Vec3, Wrench};
use fs_ga::Mat3;
use fs_math::dd::Dd;

/// A finite, uniform-density spherical source with potential zero at infinity.
///
/// Outside the sphere, the original [`CentralGravity`] owner supplies the exact
/// shell-theorem field. Inside, acceleration is linear in displacement and the
/// potential is quadratic, including the nonsingular center. This is a physical
/// density model, not point-source softening. Source position, G M and radius are
/// explicit; querying the field does not move or recoil the source.
///
/// The density jumps at the surface: potential and acceleration are continuous,
/// but the acceleration derivative has two one-sided limits. Exactly on the
/// surface, [`PointGravityField`] selects the exterior tangent. The rigid-body
/// trait uses the target's original COM mass as a monopole and supplies no tidal
/// torque or integrated finite-target energy. Continuum consumers can instead
/// sample their original material quadrature through the point-field trait.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformSphereGravity {
    central: CentralGravity,
    radius_m: f64,
}

impl UniformSphereGravity {
    /// Declare a positive finite G M [m³/s²] and source radius [m].
    pub fn new(
        frame: u64,
        center_world_m: Vec3,
        gravitational_parameter_m3_s2: f64,
        radius_m: f64,
    ) -> Result<Self, GravityError> {
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(GravityError::Invalid(
                "uniform sphere needs a positive finite radius",
            ));
        }
        Ok(Self {
            central: CentralGravity::new(frame, center_world_m, gravitational_parameter_m3_s2)?,
            radius_m,
        })
    }

    /// Acceleration [m/s²] and specific potential [J/kg], including the center.
    /// An unrepresentable tangent does not prevent this acceleration-only query.
    pub fn sample(&self, frame: u64, position_world_m: Vec3) -> Result<(Vec3, f64), GravityError> {
        let (displacement, radius) = self.geometry(frame, position_world_m)?;
        if radius >= self.radius_m {
            return self.central.sample(frame, position_world_m);
        }
        self.interior_sample(displacement, radius)
    }

    /// Specific potential change between snapshots of the same finite source.
    ///
    /// Source and sampled point may both translate. G M, radius and inertial
    /// frame must remain unchanged. The direct interior quadratic contraction
    /// retains small work hidden below either absolute potential's rounding.
    /// This does not advance the source or provide its reaction wrench.
    pub fn specific_potential_increment_to(
        &self,
        source_after: &Self,
        frame: u64,
        before_world_m: Vec3,
        after_world_m: Vec3,
    ) -> Result<f64, GravityError> {
        finite(
            self.specific_potential_change_to(source_after, frame, before_world_m, after_world_m)?
                .to_f64(),
        )
    }

    fn geometry(&self, frame: u64, position: Vec3) -> Result<(Vec3, f64), GravityError> {
        same_frame(frame, self.central.frame)?;
        let displacement = position.sub(self.central.center_world_m);
        Ok((
            displacement,
            displacement.stable_norm("uniform sphere radius")?,
        ))
    }

    fn interior_sample(
        &self,
        displacement: Vec3,
        radius: f64,
    ) -> Result<(Vec3, f64), GravityError> {
        let mu = self.central.gravitational_parameter_m3_s2;
        let q = radius / self.radius_m;
        let potential = finite(-(mu / self.radius_m) * (1.5 - 0.5 * q * q))?;
        let acceleration = |component| {
            finite(
                -scaled_product(mu, Dd::from_f64(component), Dd::ONE, self.radius_m, 3, 0)?
                    .to_f64(),
            )
        };
        Ok((
            Vec3::new(
                acceleration(displacement.x)?,
                acceleration(displacement.y)?,
                acceleration(displacement.z)?,
            ),
            potential,
        ))
    }

    fn specific_potential_change_to(
        &self,
        source_after: &Self,
        frame: u64,
        before: Vec3,
        after: Vec3,
    ) -> Result<Dd, GravityError> {
        if self.central.gravitational_parameter_m3_s2
            != source_after.central.gravitational_parameter_m3_s2
            || self.radius_m != source_after.radius_m
        {
            return Err(GravityError::Invalid(
                "moving uniform sphere must retain original G M and radius",
            ));
        }
        let (_, before_radius) = self.geometry(frame, before)?;
        let (_, after_radius) = source_after.geometry(frame, after)?;
        let before_inside = before_radius < self.radius_m;
        let after_inside = after_radius < self.radius_m;
        if !before_inside && !after_inside {
            return self.central.specific_potential_change_to(
                &source_after.central,
                frame,
                before,
                after,
            );
        }
        // Admit the same absolute potentials, independently of derivative range.
        self.sample(frame, before)?;
        source_after.sample(frame, after)?;
        let mu = self.central.gravitational_parameter_m3_s2;
        if before_inside && after_inside {
            let mut change = Dd::ZERO;
            for (a, ca, b, cb) in coordinates(
                after,
                source_after.central.center_world_m,
                before,
                self.central.center_world_m,
            ) {
                let [a, ca, b, cb] = [a, ca, b, cb].map(Dd::from_f64);
                let relative_after = a - ca;
                let relative_before = b - cb;
                let mut displacement = (a - b) - (ca - cb);
                let mut exponent = -1;
                if !displacement.hi.is_finite() {
                    // Both endpoint separations can be finite while a common
                    // source/point translation, or their difference, overflows.
                    displacement = (relative_after * Dd::from_f64(0.5))
                        - (relative_before * Dd::from_f64(0.5));
                    exponent = 0;
                }
                // Separate the endpoint products: their sum can overflow even
                // when the physical quadratic difference is representable.
                change = change
                    + scaled_product(mu, displacement, relative_after, self.radius_m, 3, exponent)?
                    + scaled_product(
                        mu,
                        displacement,
                        relative_before,
                        self.radius_m,
                        3,
                        exponent,
                    )?;
            }
            finite(change.to_f64())?;
            return Ok(change);
        }

        // Split a surface crossing at the common surface potential. Both terms
        // have the same sign; no large absolute potential is subtracted.
        let (inner, inner_point, outer, outer_point, outer_radius, sign) = if before_inside {
            (self, before, source_after, after, after_radius, 1.0)
        } else {
            (source_after, after, self, before, before_radius, -1.0)
        };
        let inner_deficit = Dd::ONE - inner.normalized_radius_squared(inner_point, self.radius_m);
        let inner_change = scaled_product(mu, inner_deficit, Dd::ONE, self.radius_m, 1, -1)?;
        let surface_ratio = Dd::from_f64(self.radius_m) / Dd::from_f64(outer_radius);
        let outer_excess = outer.normalized_radius_squared(outer_point, outer_radius)
            - surface_ratio * surface_ratio;
        let outer_factor = outer_excess / (Dd::ONE + surface_ratio);
        let outer_change = scaled_product(mu, outer_factor, Dd::ONE, self.radius_m, 1, 0)?;
        let change = (inner_change + outer_change) * Dd::from_f64(sign);
        finite(change.to_f64())?;
        Ok(change)
    }

    fn normalized_radius_squared(&self, point: Vec3, scale: f64) -> Dd {
        let divisor = Dd::from_f64(scale);
        let center = self.central.center_world_m;
        [
            (point.x, center.x),
            (point.y, center.y),
            (point.z, center.z),
        ]
        .into_iter()
        .map(|(point, center)| {
            let q = (Dd::from_f64(point) - Dd::from_f64(center)) / divisor;
            q * q
        })
        .fold(Dd::ZERO, |sum, term| sum + term)
    }
}

impl PointGravityField for UniformSphereGravity {
    fn sample_point(
        &self,
        frame: u64,
        position_world_m: Vec3,
    ) -> Result<PointGravitySample, GravityError> {
        let (displacement, radius) = self.geometry(frame, position_world_m)?;
        if radius >= self.radius_m {
            return self.central.sample_point(frame, position_world_m);
        }
        let (acceleration, potential) = self.interior_sample(displacement, radius)?;
        let slope = finite(
            -scaled_product(
                self.central.gravitational_parameter_m3_s2,
                Dd::ONE,
                Dd::ONE,
                self.radius_m,
                3,
                0,
            )?
            .to_f64(),
        )?;
        PointGravitySample {
            acceleration_world_m_s2: acceleration,
            specific_potential_j_kg: potential,
            acceleration_gradient_s2: Mat3 {
                m: [slope, 0.0, 0.0, 0.0, slope, 0.0, 0.0, 0.0, slope],
            },
        }
        .validate()
    }

    fn specific_potential_increment(
        &self,
        frame: u64,
        before: Vec3,
        after: Vec3,
    ) -> Result<Option<f64>, GravityError> {
        Ok(Some(self.specific_potential_increment_to(
            self, frame, before, after,
        )?))
    }
}

impl RigidGravityField for UniformSphereGravity {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        properties.validate()?;
        let (acceleration, potential) = self.sample(frame, pose.position_world())?;
        GravityLoad {
            wrench: Wrench {
                force_world: acceleration.scale(properties.mass()),
                torque_body: Vec3::ZERO,
            },
            potential_energy_j: properties.mass() * potential,
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
        properties.validate()?;
        let change = self.specific_potential_change_to(
            self,
            frame,
            before.position_world(),
            after.position_world(),
        )?;
        Ok(Some(finite(
            (Dd::from_f64(properties.mass()) * change).to_f64(),
        )?))
    }
}

fn coordinates(a: Vec3, ca: Vec3, b: Vec3, cb: Vec3) -> [(f64, f64, f64, f64); 3] {
    [
        (a.x, ca.x, b.x, cb.x),
        (a.y, ca.y, b.y, cb.y),
        (a.z, ca.z, b.z, cb.z),
    ]
}

// Evaluate mu*a*b/R^power * 2^adjustment without forming radius powers,
// overflowing products, or flushing a tiny factor before the final scaling.
pub(super) fn scaled_product(
    mu: f64,
    a: Dd,
    b: Dd,
    radius: f64,
    power: i32,
    adjustment: i32,
) -> Result<Dd, GravityError> {
    if a == Dd::ZERO || b == Dd::ZERO {
        return Ok(Dd::ZERO);
    }
    let (mu, mu_exp) = parts(mu);
    let (radius, radius_exp) = parts(radius);
    let (a, a_exp) = dd_parts(a);
    let (b, b_exp) = dd_parts(b);
    let mut value = Dd::from_f64(mu) * a * b;
    for _ in 0..power {
        value = value / Dd::from_f64(radius);
    }
    let exponent = mu_exp + a_exp + b_exp - power * radius_exp + adjustment;
    let hi = finite(scale_binary(value.hi, exponent))?;
    let lo = finite(scale_binary(value.lo, exponent))?;
    let result = Dd::from_pair(hi, lo);
    finite(result.to_f64())?;
    Ok(result)
}

fn dd_parts(value: Dd) -> (Dd, i32) {
    let (hi, exponent) = parts(value.hi);
    (
        Dd::from_pair(hi, scale_binary(value.lo, -exponent)),
        exponent,
    )
}

fn parts(value: f64) -> (f64, i32) {
    if value == 0.0 {
        return (value, 0);
    }
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent == 0 {
        let (mantissa, exponent) = parts(value * 18_014_398_509_481_984.0);
        return (mantissa, exponent - 54);
    }
    (
        f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52)),
        exponent - 1023,
    )
}

fn scale_binary(value: f64, exponent: i32) -> f64 {
    if value == 0.0 {
        return value;
    }
    let (mantissa, original_exponent) = parts(value);
    let exponent = exponent + original_exponent;
    if exponent > 1023 {
        return f64::INFINITY.copysign(value);
    }
    if exponent < -1075 {
        return 0.0_f64.copysign(value);
    }
    if exponent < -1022 {
        // The first product remains normal; only the final operation rounds
        // to the subnormal lattice, avoiding premature underflow.
        return (mantissa * f64::from_bits(((exponent + 1074 + 1023) as u64) << 52))
            * f64::from_bits(1);
    }
    mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
}
