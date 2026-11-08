//! Exact spherical-shell gravity for a prescribed piecewise uniform density.

use super::sphere::scaled_product;
use super::{PointGravityField, PointGravitySample};
use crate::gravity::{
    CentralGravity, GravityError, GravityLoad, RigidGravityField, finite, same_frame,
};
use crate::{MassProperties, Pose, Vec3, Wrench};
use fs_ga::Mat3;
use fs_math::dd::Dd;

/// Mass uniformly distributed between the preceding and this outer radius.
///
/// The first layer starts at radius zero. A zero mass declares a vacuum layer;
/// a positive mass declares density `3 m / (4 pi (b^3-a^3))`. Radii are meters
/// and masses kilograms. Construction never needs to form either radius cube.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialGravityLayer {
    pub outer_radius_m: f64,
    pub mass_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
struct Layer {
    inner: f64,
    outer: f64,
    mu: f64,
    inner_mu: f64,
    // (b^2 + ab + a^2)/b^2 and (b^3-a^3)/b^3, factored before subtraction.
    shape: f64,
    volume: f64,
    inner_potential: Dd,
    outer_potential: Dd,
}

/// Newtonian gravity of concentric layers with their actual uniform densities.
///
/// Enclosed masses and exterior-shell potentials are prepared in compensated
/// order. Queries locate a layer by binary search and allocate nothing. Inside
/// a cavity the surrounding shells supply constant potential and exactly zero
/// acceleration; outside the support the original point-source owner applies.
/// No artificial softening or negative-mass sphere subtraction is used.
///
/// Potential and acceleration are continuous at each density interface. The
/// point-field derivative selects the exterior-side density exactly at an
/// interface; its two one-sided limits generally differ. The rigid-body trait
/// applies the field to the target's original COM mass, with no inferred tidal
/// torque or finite-target averaging. Continuum consumers retain their own
/// material quadrature. The source is prescribed and does not recoil.
#[derive(Clone, Debug, PartialEq)]
pub struct RadialGravity {
    frame: u64,
    center_world_m: Vec3,
    layers: Vec<Layer>,
    exterior: Option<CentralGravity>,
}

impl RadialGravity {
    /// Prepare a bounded, nonempty sequence of increasing positive radii.
    ///
    /// G is positive and finite [m^3 kg^-1 s^-2]; masses are finite and
    /// nonnegative. All-zero masses describe a zero field on the declared
    /// support. The count cap is checked before allocation. Refuses an
    /// unrepresentable layer/total G M or central potential. Preparation costs
    /// O(layers) work and storage; point queries cost O(log layers).
    pub fn new(
        frame: u64,
        center_world_m: Vec3,
        gravitational_constant: f64,
        layers: &[RadialGravityLayer],
        maximum_layers: usize,
    ) -> Result<Self, GravityError> {
        if !center_world_m.is_finite()
            || !gravitational_constant.is_finite()
            || gravitational_constant <= 0.0
            || layers.is_empty()
            || layers.len() > maximum_layers
        {
            return Err(GravityError::Invalid(
                "radial gravity needs finite center, positive G and a bounded nonempty profile",
            ));
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(layers.len())
            .map_err(|_| GravityError::Invalid("radial gravity allocation failed"))?;
        let mut inner = 0.0;
        let mut enclosed = Dd::ZERO;
        for layer in layers {
            let outer = layer.outer_radius_m;
            if !outer.is_finite()
                || outer <= inner
                || !layer.mass_kg.is_finite()
                || layer.mass_kg < 0.0
            {
                return Err(GravityError::Invalid(
                    "radial gravity needs increasing positive radii and nonnegative masses",
                ));
            }
            let mu = scaled_product(
                gravitational_constant,
                Dd::from_f64(layer.mass_kg),
                Dd::ONE,
                1.0,
                0,
                0,
            )?
            .to_f64();
            if layer.mass_kg > 0.0 && mu == 0.0 {
                return Err(GravityError::Invalid("unrepresentable layer G M"));
            }
            let q = inner / outer;
            let shape = 1.0 + q + q * q;
            let volume = ((outer - inner) / outer) * shape;
            let inner_potential = scaled_product(
                mu,
                Dd::from_f64((1.0 + q) / shape),
                Dd::from_f64(3.0),
                outer,
                1,
                -1,
            )?;
            prepared.push(Layer {
                inner,
                outer,
                mu,
                inner_mu: finite(enclosed.to_f64())?,
                shape,
                volume,
                inner_potential,
                outer_potential: Dd::ZERO,
            });
            enclosed = enclosed + Dd::from_f64(mu);
            finite(enclosed.to_f64())?;
            inner = outer;
        }
        let mut outside = Dd::ZERO;
        for layer in prepared.iter_mut().rev() {
            layer.outer_potential = outside;
            outside = outside + layer.inner_potential;
            finite(outside.to_f64())?;
        }
        let total_mu = enclosed.to_f64();
        Ok(Self {
            frame,
            center_world_m,
            layers: prepared,
            exterior: if total_mu == 0.0 {
                None
            } else {
                Some(CentralGravity::new(frame, center_world_m, total_mu)?)
            },
        })
    }

    /// Acceleration [m/s^2] and specific potential [J/kg], zero at infinity.
    /// This admits finite acceleration even when its tangent is unrepresentable.
    pub fn sample(&self, frame: u64, point: Vec3) -> Result<(Vec3, f64), GravityError> {
        let (offset, radius) = self.geometry(frame, point)?;
        let index = self.layer_index(radius);
        if index == self.layers.len() {
            return match self.exterior {
                Some(source) => source.sample(frame, point),
                None => Ok((Vec3::ZERO, 0.0)),
            };
        }
        self.interior_sample(index, offset, radius)
    }

    fn interior_sample(
        &self,
        index: usize,
        offset: Vec3,
        radius: f64,
    ) -> Result<(Vec3, f64), GravityError> {
        let layer = &self.layers[index];
        if radius == 0.0 {
            return Ok((
                Vec3::ZERO,
                finite(-(layer.inner_potential + layer.outer_potential).to_f64())?,
            ));
        }
        let q = layer.inner / radius;
        let thickness = (radius - layer.inner) / (layer.outer - layer.inner);
        let fraction = (1.0 + q + q * q) / layer.shape;
        let acceleration = scaled_product(layer.inner_mu, Dd::ONE, Dd::ONE, radius, 2, 0)?
            + scaled_product(
                layer.mu,
                Dd::from_f64(thickness),
                Dd::from_f64(fraction),
                layer.outer,
                2,
                0,
            )?;
        let enclosed_potential = scaled_product(layer.inner_mu, Dd::ONE, Dd::ONE, radius, 1, 0)?
            + scaled_product(
                layer.mu,
                Dd::from_f64(thickness),
                Dd::from_f64((radius / layer.outer) * fraction),
                layer.outer,
                1,
                0,
            )?;
        let remaining = (layer.outer - radius) / (layer.outer - layer.inner);
        let remaining_potential = scaled_product(
            layer.mu,
            Dd::from_f64(remaining),
            Dd::from_f64(3.0 * (1.0 + radius / layer.outer) / layer.shape),
            layer.outer,
            1,
            -1,
        )?;
        let magnitude = finite(acceleration.to_f64())?;
        let acceleration = Vec3::new(
            -magnitude * (offset.x / radius),
            -magnitude * (offset.y / radius),
            -magnitude * (offset.z / radius),
        );
        let potential =
            finite(-(enclosed_potential + remaining_potential + layer.outer_potential).to_f64())?;
        if !acceleration.is_finite() {
            return Err(GravityError::Invalid("unrepresentable radial acceleration"));
        }
        Ok((acceleration, potential))
    }

    fn geometry(&self, frame: u64, point: Vec3) -> Result<(Vec3, f64), GravityError> {
        same_frame(frame, self.frame)?;
        let offset = point.sub(self.center_world_m);
        Ok((offset, offset.stable_norm("radial gravity radius")?))
    }

    fn layer_index(&self, radius: f64) -> usize {
        self.layers.partition_point(|layer| radius >= layer.outer)
    }

    fn precise_layer_index(&self, radius: Dd) -> usize {
        self.layers
            .partition_point(|layer| !radius.lt(Dd::from_f64(layer.outer)))
    }

    // Positive G M_enclosed/r, never a signed point deficit plus a sphere.
    fn enclosed_potential(&self, index: usize, radius: Dd) -> Result<Dd, GravityError> {
        if radius == Dd::ZERO {
            return Ok(Dd::ZERO);
        }
        if index == self.layers.len() {
            return match self.exterior {
                Some(source) => scaled_product(
                    source.gravitational_parameter_m3_s2,
                    Dd::from_f64(radius.hi) / radius,
                    Dd::ONE,
                    radius.hi,
                    1,
                    0,
                ),
                None => Ok(Dd::ZERO),
            };
        }
        let layer = &self.layers[index];
        let inner = Dd::from_f64(layer.inner);
        let outer = Dd::from_f64(layer.outer);
        let q = inner / radius;
        let thickness = (radius - inner) / (outer - inner);
        let fraction = (Dd::ONE + q + q * q) / Dd::from_f64(layer.shape);
        Ok(scaled_product(
            layer.inner_mu,
            Dd::from_f64(radius.hi) / radius,
            Dd::ONE,
            radius.hi,
            1,
            0,
        )? + scaled_product(
            layer.mu,
            thickness,
            (radius / outer) * fraction,
            layer.outer,
            1,
            0,
        )?)
    }

    fn precise_radius(&self, point: Vec3, radius: f64) -> Result<Dd, GravityError> {
        if radius == 0.0 {
            return Ok(Dd::ZERO);
        }
        let mut square = Dd::ZERO;
        for (point, center) in components(point, self.center_world_m) {
            let q = (Dd::from_f64(point) - Dd::from_f64(center)) / Dd::from_f64(radius);
            square = square + q * q;
        }
        scaled_product(radius, square.sqrt(), Dd::ONE, 1.0, 0, 0)
    }

    fn potential_change(&self, frame: u64, before: Vec3, after: Vec3) -> Result<Dd, GravityError> {
        let (_, rb) = self.geometry(frame, before)?;
        let (_, ra) = self.geometry(frame, after)?;
        self.sample(frame, before)?;
        self.sample(frame, after)?;
        if self.exterior.is_none() {
            return Ok(Dd::ZERO);
        }
        if rb >= self.layers.last().unwrap().outer && ra >= self.layers.last().unwrap().outer {
            return self
                .exterior
                .unwrap()
                .specific_potential_change(frame, before, after);
        }
        let before_radius = self.precise_radius(before, rb)?;
        let after_radius = self.precise_radius(after, ra)?;
        let scale = rb.max(ra);
        let difference = if scale == 0.0 {
            Dd::ZERO
        } else if rb < 0.5 * ra || ra < 0.5 * rb {
            after_radius - before_radius
        } else {
            let mut square_change = Dd::ZERO;
            let divisor = Dd::from_f64(scale);
            for ((a, center), (b, _)) in components(after, self.center_world_m)
                .into_iter()
                .zip(components(before, self.center_world_m))
            {
                let a = Dd::from_f64(a);
                let b = Dd::from_f64(b);
                let center = Dd::from_f64(center);
                let qa = (a - center) / divisor;
                let qb = (b - center) / divisor;
                let travel = a - b;
                let travel = if travel.hi.is_finite() {
                    travel / divisor
                } else {
                    qa - qb
                };
                square_change = square_change + travel * (qa + qb);
            }
            scaled_product(
                scale,
                square_change / (before_radius / divisor + after_radius / divisor),
                Dd::ONE,
                1.0,
                0,
                0,
            )?
        };
        if difference == Dd::ZERO {
            return Ok(Dd::ZERO);
        }
        let (first, last, difference, sign) = if difference.lt(Dd::ZERO) {
            (after_radius, before_radius, -difference, -1.0)
        } else {
            (before_radius, after_radius, difference, 1.0)
        };
        let first_index = self.precise_layer_index(first);
        let last_index = self.precise_layer_index(last);
        let mut change = Dd::ZERO;
        let mut start = first;
        for index in first_index..=last_index {
            let end = if index == last_index {
                last
            } else {
                Dd::from_f64(self.layers[index].outer)
            };
            let distance = if first_index == last_index {
                difference
            } else {
                end - start
            };
            change = change + self.layer_change(index, start, end, distance)?;
            start = end;
        }
        finite(change.to_f64())?;
        Ok(change * Dd::from_f64(sign))
    }

    // For an outward interval d: Delta phi = mu(r)*d/[r(r+d)]
    // + k*d^2*(3r+d)/[2(r+d)]. Both terms are nonnegative, even at the
    // inner boundary of a hollow shell where the linear contribution is zero.
    fn layer_change(
        &self,
        index: usize,
        start: Dd,
        end: Dd,
        distance: Dd,
    ) -> Result<Dd, GravityError> {
        if distance == Dd::ZERO {
            return Ok(Dd::ZERO);
        }
        if distance.lt(Dd::ZERO) || end.hi <= 0.0 {
            return Err(GravityError::Invalid("inconsistent radial work interval"));
        }
        let enclosed = self.enclosed_potential(index, start)?.to_f64();
        let first = scaled_product(enclosed, distance, Dd::from_f64(end.hi) / end, end.hi, 1, 0)?;
        if index == self.layers.len() {
            return Ok(first);
        }
        let layer = &self.layers[index];
        let width = Dd::from_f64(layer.outer) - Dd::from_f64(layer.inner);
        let factor = (distance / width) * (Dd::ONE + Dd::from_f64(2.0) * (start / end))
            / Dd::from_f64(layer.shape);
        Ok(first + scaled_product(layer.mu, distance, factor, layer.outer, 2, -1)?)
    }
}

impl PointGravityField for RadialGravity {
    fn sample_point(&self, frame: u64, point: Vec3) -> Result<PointGravitySample, GravityError> {
        let (offset, radius) = self.geometry(frame, point)?;
        let index = self.layer_index(radius);
        if index == self.layers.len() {
            return match self.exterior {
                Some(source) => source.sample_point(frame, point),
                None => Ok(PointGravitySample {
                    acceleration_world_m_s2: Vec3::ZERO,
                    specific_potential_j_kg: 0.0,
                    acceleration_gradient_s2: Mat3 { m: [0.0; 9] },
                }),
            };
        }
        let (acceleration, potential) = self.interior_sample(index, offset, radius)?;
        let layer = &self.layers[index];
        let inner_tide = if radius == 0.0 {
            0.0
        } else {
            scaled_product(layer.inner_mu, Dd::ONE, Dd::ONE, radius, 3, 0)?.to_f64()
        };
        let density = scaled_product(
            layer.mu,
            Dd::ONE / Dd::from_f64(layer.volume),
            Dd::ONE,
            layer.outer,
            3,
            0,
        )?
        .to_f64();
        let n = if radius == 0.0 {
            [0.0; 3]
        } else {
            [offset.x / radius, offset.y / radius, offset.z / radius]
        };
        let q = if radius == 0.0 {
            0.0
        } else {
            layer.inner / radius
        };
        let q3 = q * q * q;
        let shell_tangent = if radius == 0.0 {
            1.0
        } else {
            ((radius - layer.inner) / radius) * (1.0 + q + q * q)
        };
        let scale = inner_tide.max(density);
        let mut gradient = [0.0; 9];
        if scale > 0.0 {
            for i in 0..3 {
                for j in 0..3 {
                    let identity = if i == j { 1.0 } else { 0.0 };
                    let point_factor = 3.0 * n[i] * n[j] - identity;
                    let shell_factor = shell_tangent * identity + 3.0 * q3 * n[i] * n[j];
                    let value = (Dd::from_f64(inner_tide) / Dd::from_f64(scale))
                        * Dd::from_f64(point_factor)
                        - (Dd::from_f64(density) / Dd::from_f64(scale))
                            * Dd::from_f64(shell_factor);
                    gradient[3 * i + j] =
                        scaled_product(scale, value, Dd::ONE, 1.0, 0, 0)?.to_f64();
                }
            }
        }
        PointGravitySample {
            acceleration_world_m_s2: acceleration,
            specific_potential_j_kg: potential,
            acceleration_gradient_s2: Mat3 { m: gradient },
        }
        .validate()
    }

    fn specific_potential_increment(
        &self,
        frame: u64,
        before: Vec3,
        after: Vec3,
    ) -> Result<Option<f64>, GravityError> {
        Ok(Some(finite(
            self.potential_change(frame, before, after)?.to_f64(),
        )?))
    }
}

impl RigidGravityField for RadialGravity {
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
        let change =
            self.potential_change(frame, before.position_world(), after.position_world())?;
        Ok(Some(finite(
            (Dd::from_f64(properties.mass()) * change).to_f64(),
        )?))
    }
}

fn components(point: Vec3, center: Vec3) -> [(f64, f64); 3] {
    [
        (point.x, center.x),
        (point.y, center.y),
        (point.z, center.z),
    ]
}
