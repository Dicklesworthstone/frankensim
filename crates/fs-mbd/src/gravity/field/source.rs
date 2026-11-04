//! A native extended source sampled by point and continuum consumers.

use super::{PointGravityField, PointGravitySample, central_sample};
use crate::gravity::{
    CentralGravity, GravityError, finite, quadrupole_correction, quadrupole_inertia_trace,
};
use crate::{MassProperties, Pose, Vec3};
use fs_ga::Mat3;

/// Caller-declared COM support radius and far-field approximation domain.
/// The inertia/radius check is necessary, not a reconstruction of mass support
/// or a certified bound on neglected higher multipoles.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadrupoleSourceExtent {
    /// Bounding support radius measured from the source COM [m].
    pub radius_m: f64,
    /// Positive and strictly less than one; each sample requires radius/r <= it.
    pub maximum_radius_ratio: f64,
}

impl QuadrupoleSourceExtent {
    /// Admit the declaration and return its minimum permitted COM distance.
    pub fn minimum_sample_radius_m(self) -> Result<f64, GravityError> {
        if !self.radius_m.is_finite()
            || self.radius_m <= 0.0
            || !self.maximum_radius_ratio.is_finite()
            || self.maximum_radius_ratio <= 0.0
            || self.maximum_radius_ratio >= 1.0
        {
            return Err(GravityError::Invalid(
                "quadrupole source needs positive radius and ratio in (0,1)",
            ));
        }
        finite(self.radius_m / self.maximum_radius_ratio)
    }
}

/// Monopole plus the original rigid source's second mass moment.
///
/// With source COM separation r*n, full world inertia I and source mass M,
/// phi = -G M/r + G (3 n.I.n - tr I)/(2 r³). The original native quadrupole
/// owner evaluates the correction to potential and acceleration. Its analytic
/// spatial derivative supplies the same potential's Hessian to solid tangents.
/// Isotropic inertia gives exactly the original central-field arithmetic.
///
/// This is a stationary pose snapshot. A recoiling/rotating interaction must
/// rebuild it at the actual source pose and apply the opposite force AND the
/// source COM torque. Sampling alone does not advance or torque the source.
/// No allocation, eigensolve, spatial differencing, higher multipoles, resolved
/// near field, interior gravity, self-gravity or quadrupole-quadrupole product
/// is introduced. The declared radius/ratio is not a truncation-error bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadrupoleSourceGravity {
    central: CentralGravity,
    pose: Pose,
    properties: MassProperties,
    gravitational_constant_m3_kg_s2: f64,
    extent: QuadrupoleSourceExtent,
    trace: f64,
    inertia_world: Mat3,
    isotropic: bool,
}

impl QuadrupoleSourceGravity {
    /// Bind the original source mass, full inertia, COM pose and explicit G.
    pub fn new(
        frame: u64,
        pose: Pose,
        properties: MassProperties,
        gravitational_constant_m3_kg_s2: f64,
        extent: QuadrupoleSourceExtent,
    ) -> Result<Self, GravityError> {
        properties.validate()?;
        extent.minimum_sample_radius_m()?;
        if !gravitational_constant_m3_kg_s2.is_finite() || gravitational_constant_m3_kg_s2 <= 0.0 {
            return Err(GravityError::Invalid(
                "positive finite Newtonian constant required",
            ));
        }
        let trace = quadrupole_inertia_trace(properties, extent.radius_m)?;
        let central = CentralGravity::new(
            frame,
            pose.position_world(),
            gravitational_constant_m3_kg_s2 * properties.mass(),
        )?;
        let principal = properties.principal_inertia_body();
        let moments = [principal.x, principal.y, principal.z];
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
        .map(|axis| {
            pose.orientation().rotate_body_to_world(
                properties
                    .principal_axes_to_body()
                    .rotate_body_to_world(axis),
            )
        })
        .map(|axis| [axis.x, axis.y, axis.z]);
        let mut matrix = [0.0; 9];
        for i in 0..3 {
            for j in i..3 {
                let entry = finite((0..3).map(|k| moments[k] * axes[k][i] * axes[k][j]).sum())?;
                matrix[3 * i + j] = entry;
                matrix[3 * j + i] = entry;
            }
        }
        Ok(Self {
            central,
            pose,
            properties,
            gravitational_constant_m3_kg_s2,
            extent,
            trace,
            inertia_world: Mat3 { m: matrix },
            isotropic: principal.x == principal.y && principal.y == principal.z,
        })
    }
}

impl PointGravityField for QuadrupoleSourceGravity {
    fn sample_point(
        &self,
        frame: u64,
        position_world_m: Vec3,
    ) -> Result<PointGravitySample, GravityError> {
        let (direction, radius, potential) = self.central.radial(frame, position_world_m)?;
        if self.extent.radius_m / radius > self.extent.maximum_radius_ratio {
            return Err(GravityError::OutsideExpansionDomain);
        }
        let mut sample = central_sample(direction, radius, potential)?;
        if self.isotropic {
            return Ok(sample);
        }
        // Unit point mass: swapping source and target reverses the force, while
        // reversing n reverses the native correction a second time.
        let g_over_r3 = finite((self.gravitational_constant_m3_kg_s2 / radius) / radius / radius)?;
        let correction = quadrupole_correction(
            self.pose,
            self.properties,
            direction,
            radius,
            g_over_r3,
            self.trace,
        )?;
        sample.acceleration_world_m_s2 = sample
            .acceleration_world_m_s2
            .add(correction.wrench.force_world);
        sample.specific_potential_j_kg += correction.potential_energy_j;
        let n = [direction.x, direction.y, direction.z];
        let inertia = self.inertia_world.m;
        let b: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| inertia[3 * i + j] * n[j]).sum());
        let a: f64 = (0..3).map(|i| n[i] * b[i]).sum();
        let scale = finite(1.5 * (g_over_r3 / radius) / radius)?;
        // Dg_Q = 3G/(2r^5) [(5a-trI)1 -2I +10(nb'+bn')
        //                         +5(trI-7a)nn']; it is symmetric and trace-free.
        for i in 0..3 {
            for j in i..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                let derivative = scale
                    * ((5.0 * a - self.trace) * delta - 2.0 * inertia[3 * i + j]
                        + 10.0 * (n[i] * b[j] + b[i] * n[j])
                        + 5.0 * (self.trace - 7.0 * a) * n[i] * n[j]);
                let entry = sample.acceleration_gradient_s2.m[3 * i + j] + derivative;
                sample.acceleration_gradient_s2.m[3 * i + j] = entry;
                sample.acceleration_gradient_s2.m[3 * j + i] = entry;
            }
        }
        sample.validate()
    }
}
