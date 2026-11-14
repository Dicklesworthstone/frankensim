//! One fixed-origin, constant-angular-velocity boundary in the original GA owner.
use fs_ga::{GaError, Se3, So3, So3Tangent, Vec3};

/// Admitted rigid rotation, anchored at its original pose and clock origin.
/// The intended origin is unchanged; samples retain the original SE(3)
/// representation, roundoff and geometric admission. No support orbit is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformRotationFrame {
    initial_pose: Se3,
    angular_velocity_body: Vec3,
    angular_velocity_world: Vec3,
}

#[derive(Clone, Debug, PartialEq)]
pub enum UniformRotationError {
    Invalid(&'static str),
    Geometry(GaError),
}
impl From<GaError> for UniformRotationError {
    fn from(error: GaError) -> Self {
        Self::Geometry(error)
    }
}
impl core::fmt::Display for UniformRotationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "uniform supported rotation: {self:?}")
    }
}
impl std::error::Error for UniformRotationError {}

impl UniformRotationFrame {
    /// A constant body angular velocity with a fixed support origin.
    /// Zero spin is admitted and retains the initial pose bitwise.
    pub fn new(
        initial_pose: Se3,
        angular_velocity_body: Vec3,
    ) -> Result<Self, UniformRotationError> {
        if !finite(angular_velocity_body) {
            return Err(UniformRotationError::Invalid(
                "finite angular velocity required",
            ));
        }
        let angular_velocity_world = initial_pose.rotation().rotate(angular_velocity_body)?;
        if !finite(angular_velocity_world) {
            return Err(UniformRotationError::Invalid(
                "world angular velocity is not representable",
            ));
        }
        Ok(Self {
            initial_pose,
            angular_velocity_body,
            angular_velocity_world,
        })
    }

    pub fn angular_velocity_body(&self) -> Vec3 {
        self.angular_velocity_body
    }
    pub fn angular_velocity_world(&self) -> Vec3 {
        self.angular_velocity_world
    }

    /// Original pose after nonnegative physical elapsed seconds. Rotation uses
    /// R(t)=R(0) Exp(omega_body*t), at the original represented support origin.
    /// Native SE(3) encoding can round the returned translation and can refuse
    /// poorly represented distant poses. No orbital state is integrated here.
    pub fn pose_at_elapsed(&self, elapsed_s: f64) -> Result<Se3, UniformRotationError> {
        if !elapsed_s.is_finite() || elapsed_s < 0.0 {
            return Err(UniformRotationError::Invalid(
                "finite nonnegative elapsed time required",
            ));
        }
        if elapsed_s == 0.0 || self.angular_velocity_body == Vec3::default() {
            return Ok(self.initial_pose);
        }
        let increment = self.angular_velocity_body.scale(elapsed_s);
        if !finite(increment) || magnitude(increment) == 0.0 {
            return Err(UniformRotationError::Invalid(
                "rotation increment is not representable",
            ));
        }
        let rotation = self
            .initial_pose
            .rotation()
            .compose(So3::exp(So3Tangent::new(increment))?)?;
        Ok(Se3::from_parts(rotation, self.initial_pose.translation())?)
    }

    /// Whether a held world vector is invariant under this rotation. Inputs are
    /// never projected or modified. Scaled cross products avoid overflow and
    /// allow only 64 machine epsilons for encoded-axis roundoff, independently
    /// of physical force units. A zero vector or zero spin is invariant.
    pub fn admits_world_vector(&self, value: Vec3) -> bool {
        if !finite(value) {
            return false;
        }
        let a = magnitude(self.angular_velocity_world);
        let b = magnitude(value);
        if a == 0.0 || b == 0.0 {
            return true;
        }
        // Component division also works for subnormal vectors, whose inverse
        // scale would overflow. The original compensated cross product is used.
        let axis = divide(self.angular_velocity_world, a);
        let vector = divide(value, b);
        magnitude(axis.cross(vector)) <= 64.0 * f64::EPSILON
    }
}

fn finite(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}
fn magnitude(v: Vec3) -> f64 {
    v.x.abs().max(v.y.abs()).max(v.z.abs())
}
fn divide(v: Vec3, scale: f64) -> Vec3 {
    Vec3::new(v.x / scale, v.y / scale, v.z / scale)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: Vec3, b: Vec3, tolerance: f64) {
        assert!(magnitude(a - b) <= tolerance, "{a:?} != {b:?}");
    }

    #[test]
    fn arbitrary_initial_orientation_retains_the_world_axis_and_fixed_origin() {
        let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.7, 0.2))).unwrap();
        let origin = Vec3::new(3.0, -2.0, 5.0);
        let frame = UniformRotationFrame::new(
            Se3::from_parts(rotation, origin).unwrap(),
            Vec3::new(0.4, -0.2, 0.7),
        )
        .unwrap();
        let point = Vec3::new(0.7, -0.3, 0.1);
        let t = 0.8;
        let eps = 1e-5;
        let pose = frame.pose_at_elapsed(t).unwrap();
        close(pose.translation(), origin, 2e-13);
        close(
            pose.rotation()
                .rotate(frame.angular_velocity_body())
                .unwrap(),
            frame.angular_velocity_world(),
            5e-15,
        );
        // Independent circular-trajectory velocity, formed relative to the
        // actual support origin rather than about the arbitrary world origin.
        let plus = frame
            .pose_at_elapsed(t + eps)
            .unwrap()
            .rotation()
            .rotate(point)
            .unwrap();
        let minus = frame
            .pose_at_elapsed(t - eps)
            .unwrap()
            .rotation()
            .rotate(point)
            .unwrap();
        let arm = pose.rotation().rotate(point).unwrap();
        close(
            (plus - minus).scale(0.5 / eps),
            frame.angular_velocity_world().cross(arm),
            2e-10,
        );
    }

    #[test]
    fn distant_rotation_keeps_its_origin_and_turn() {
        let origin = Vec3::new(0.0, 0.0, 1_125_899_906_842_624.0);
        let frame = UniformRotationFrame::new(
            Se3::from_parts(So3::identity(), origin).unwrap(),
            Vec3::new(0.0, 0.0, 3.0),
        )
        .unwrap();
        assert_eq!(frame.pose_at_elapsed(0.0).unwrap().translation(), origin);
        // 2^50 m out, the rotated sample's ideal motor coordinates round in
        // proportion to the translation, whose own ulp is 0.25 m. fs-ga's SE(3)
        // unit check takes them relative to the translation, so this rigid
        // motion is admitted. No normalization or origin movement is
        // introduced: the sample keeps its origin (to 4 ulps) and its 0.3 rad
        // turn about the axis.
        let pose = frame.pose_at_elapsed(0.1).unwrap();
        close(pose.translation(), origin, 1.0);
        close(
            pose.rotation().rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            Vec3::new(0.3_f64.cos(), 0.3_f64.sin(), 0.0),
            1e-15,
        );
    }

    #[test]
    fn invariant_load_admission_is_scale_independent_and_preserves_small_nonaxial_fields() {
        let frame = UniformRotationFrame::new(Se3::identity(), Vec3::new(0.0, 0.0, 3.0)).unwrap();
        for scale in [1e-300, 1.0, 1e300] {
            assert!(frame.admits_world_vector(Vec3::new(0.0, 0.0, -scale)));
            assert!(!frame.admits_world_vector(Vec3::new(scale, 0.0, scale)));
            assert!(!frame.admits_world_vector(Vec3::new(1e-9 * scale, 0.0, scale)));
        }
        assert!(frame.admits_world_vector(Vec3::default()));
        assert!(!frame.admits_world_vector(Vec3::new(f64::NAN, 0.0, 0.0)));
        let tiny =
            UniformRotationFrame::new(Se3::identity(), Vec3::new(0.0, 0.0, f64::from_bits(1)))
                .unwrap();
        assert!(tiny.admits_world_vector(Vec3::new(0.0, 0.0, f64::from_bits(1))));
        assert!(!tiny.admits_world_vector(Vec3::new(f64::from_bits(1), 0.0, 0.0)));
    }

    #[test]
    fn zero_spin_and_zero_elapsed_preserve_pose_and_bad_motion_refuses() {
        let pose = Se3::from_parts(So3::identity(), Vec3::new(1e15, 2.0, 3.0)).unwrap();
        let stationary = UniformRotationFrame::new(pose, Vec3::default()).unwrap();
        assert_eq!(stationary.pose_at_elapsed(f64::MAX).unwrap(), pose);
        assert!(stationary.admits_world_vector(Vec3::new(1.0, 2.0, 3.0)));
        let frame = UniformRotationFrame::new(pose, Vec3::new(0.0, 0.0, 3.0)).unwrap();
        assert_eq!(frame.pose_at_elapsed(0.0).unwrap(), pose);
        assert!(frame.pose_at_elapsed(-1.0).is_err());
        assert!(frame.pose_at_elapsed(f64::INFINITY).is_err());
        assert!(frame.pose_at_elapsed(f64::MAX).is_err());
        assert!(UniformRotationFrame::new(pose, Vec3::new(f64::NAN, 0.0, 0.0)).is_err());
    }
}
