//! Shared material inertia in the centre-of-mass principal frame required by
//! the existing free-body integrator and point-contact/impulse APIs.
//!
//! Eigenvectors come from the shared L1 Jacobi kernel, not another eigensolver.
//! The actual quaternion-frame tensor reconstruction is checked before use.
//! Degenerate eigenvalues have non-unique axes: no smooth frame continuation
//! or derivative through eigenvector selection is claimed.
use crate::ResolvedBody;
use fs_ga::{GaError, Mat3, Quat, Se3, So3, Twist, Vec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::{DynamicsError, MassProperties, Pose, RigidBodyState, UnitQuaternion};

/// Momentum-preserving changes of an existing closed body's configuration.
pub mod transfer;

/// Failure never replaces the source body or an existing dynamics state.
#[derive(Clone, Debug, PartialEq)]
pub enum PrincipalFrameError {
    /// Nonfinite input, unrepresentable arithmetic or failed tensor reconstruction.
    Invalid(&'static str),
    /// Existing pose/rotation owner refused a frame.
    Frame(GaError),
    /// Existing free-body owner refused properties or state.
    Dynamics(DynamicsError),
    /// Cancellation before or after the bounded 3x3 eigendecomposition.
    Cancelled,
}
impl core::fmt::Display for PrincipalFrameError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "material principal frame: {s}"),
            Self::Frame(e) => write!(f, "material principal frame: {e}"),
            Self::Dynamics(e) => write!(f, "material principal dynamics: {e}"),
            Self::Cancelled => f.write_str("material principal frame: cancelled"),
        }
    }
}
impl std::error::Error for PrincipalFrameError {}
impl From<GaError> for PrincipalFrameError {
    fn from(e: GaError) -> Self {
        Self::Frame(e)
    }
}
impl From<DynamicsError> for PrincipalFrameError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}

/// One source body's exact mass and numerically admitted principal-frame view.
/// Original geometry and material receipts remain reachable through `source`.
#[derive(Clone, Debug)]
pub struct PrincipalBody<'a> {
    source: &'a ResolvedBody,
    frame: PrincipalFrame,
}

/// Principal-frame view of one complete, already admitted spatial inertia.
///
/// The same lowering serves reference bodies and their thermally or mechanically
/// deformed configurations. It carries the actual centroid and full tensor;
/// material/geometry lineage remains with the configuration that supplied it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrincipalFrame {
    inertia: SpatialInertia,
    properties: MassProperties,
    body_from_principal: Se3,
    relative_tensor_residual: f64,
}
impl ResolvedBody {
    /// Lower arbitrary full inertia to the existing free-body/contact owner.
    /// The caller supplies an entrywise tensor residual budget relative to the
    /// largest original tensor entry, in (0,1e-8]. This is a numerical residual,
    /// not an experimental error bound or a propagated source uncertainty.
    /// Cancellation surrounds the fixed-size, at-most-60-sweep shared kernel.
    pub fn principal_body(
        &self,
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PrincipalBody<'_>, PrincipalFrameError> {
        let frame =
            PrincipalFrame::from_inertia(self.inertia(), relative_tolerance, &mut cancelled)?;
        Ok(PrincipalBody {
            source: self,
            frame,
        })
    }
}
impl PrincipalFrame {
    /// Lower the same complete inertia used by articulated dynamics. This is
    /// not material resolution or a live-state update. The tensor residual
    /// tolerance and bounded eigendecomposition match `principal_body`.
    pub fn from_inertia(
        inertia: SpatialInertia,
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PrincipalFrameError> {
        poll(&mut cancelled)?;
        if !relative_tolerance.is_finite() || relative_tolerance <= 0.0 || relative_tolerance > 1e-8
        {
            return Err(PrincipalFrameError::Invalid(
                "principal tensor tolerance must be in (0,1e-8]",
            ));
        }
        let tensor = inertia.inertia_com();
        let scale = tensor.max_abs();
        let normalized = tensor.m.map(|x| x / scale);
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&normalized, 3);
        let mut rotation = Mat3 { m: [0.0; 9] };
        rotation.m.copy_from_slice(&vectors);
        let determinant = rotation.m[0]
            * (rotation.m[4] * rotation.m[8] - rotation.m[5] * rotation.m[7])
            - rotation.m[1] * (rotation.m[3] * rotation.m[8] - rotation.m[5] * rotation.m[6])
            + rotation.m[2] * (rotation.m[3] * rotation.m[7] - rotation.m[4] * rotation.m[6]);
        if !determinant.is_finite() || (determinant.abs() - 1.0).abs() > 1e-10 {
            return Err(PrincipalFrameError::Invalid(
                "principal axes are not an orthonormal frame",
            ));
        }
        // One eigenvector sign may change without changing the tensor. A pose
        // must be right-handed, unlike an arbitrary eigensolver basis.
        if determinant < 0.0 {
            for row in 0..3 {
                rotation.m[3 * row] = -rotation.m[3 * row];
            }
        }
        let rotation = rotation_from_matrix(rotation)?;
        let r = rotation.matrix();
        let diagonal = Mat3 {
            m: [
                values[0], 0.0, 0.0, 0.0, values[1], 0.0, 0.0, 0.0, values[2],
            ],
        };
        let reconstructed = r.compose(diagonal).compose(r.transpose());
        let mut residual = 0.0_f64;
        for (actual, expected) in reconstructed.m.into_iter().zip(normalized) {
            let error = (actual - expected).abs();
            if !error.is_finite() {
                return Err(PrincipalFrameError::Invalid(
                    "nonfinite principal tensor reconstruction",
                ));
            }
            residual = residual.max(error);
        }
        if residual > relative_tolerance {
            return Err(PrincipalFrameError::Invalid(
                "principal tensor reconstruction exceeds the declared budget",
            ));
        }
        let properties = MassProperties::new(
            inertia.mass(),
            fs_mbd::Vec3::ZERO,
            fs_mbd::Vec3::new(scale * values[0], scale * values[1], scale * values[2]),
        )?;
        let body_from_principal = Se3::from_parts(rotation, inertia.center_of_mass())?;
        poll(&mut cancelled)?;
        Ok(Self {
            inertia,
            properties,
            body_from_principal,
            relative_tensor_residual: residual,
        })
    }
}
impl<'a> PrincipalBody<'a> {
    /// Original region geometry and sourced material states, not a copied preset.
    pub const fn source(&self) -> &'a ResolvedBody {
        self.source
    }
    /// The same frame lowering used by continuing deformed-body motion.
    pub const fn frame(&self) -> &PrincipalFrame {
        &self.frame
    }
    /// Same MassProperties accepted by free-body stepping and point impulses.
    pub const fn properties(&self) -> MassProperties {
        self.frame.properties()
    }
    /// Maps principal-centre coordinates into the original geometry/body frame.
    pub const fn body_from_principal(&self) -> Se3 {
        self.frame.body_from_principal()
    }
    /// Maximum actual frame-reconstructed tensor error / largest input entry.
    pub const fn relative_tensor_residual(&self) -> f64 {
        self.frame.relative_tensor_residual()
    }
    /// Initialize this source body's pose and twist through the shared lowering.
    pub fn state_from_kinematics(
        &self,
        world_from_body: Se3,
        twist_body: Twist,
    ) -> Result<RigidBodyState, PrincipalFrameError> {
        self.frame
            .state_from_kinematics(world_from_body, twist_body)
    }
    /// Recover the original geometry's world pose after free-body stepping.
    pub fn world_from_body(&self, state: RigidBodyState) -> Result<Se3, PrincipalFrameError> {
        self.frame.world_from_body(state)
    }
    /// Convert an authored contact point to the principal centre-relative arm.
    pub fn arm_for_body_point(
        &self,
        point_body: Vec3,
    ) -> Result<fs_mbd::Vec3, PrincipalFrameError> {
        self.frame.arm_for_body_point(point_body)
    }
}
impl PrincipalFrame {
    /// Actual complete inertia in the configuration's authored body frame.
    pub const fn inertia(&self) -> SpatialInertia {
        self.inertia
    }
    /// Same MassProperties accepted by free-body stepping and point impulses.
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    /// Maps principal-centre coordinates into the original geometry/body frame.
    pub const fn body_from_principal(&self) -> Se3 {
        self.body_from_principal
    }
    /// Maximum actual frame-reconstructed tensor error / largest input entry.
    pub const fn relative_tensor_residual(&self) -> f64 {
        self.relative_tensor_residual
    }

    /// Initialize a principal-frame dynamics state from an original body pose
    /// and body-frame twist at its authored origin. The centre velocity includes
    /// omega cross centre_of_mass; angular momentum uses the same principal
    /// tensor. This is initialization, NOT momentum-conserving live replacement.
    pub fn state_from_kinematics(
        &self,
        world_from_body: Se3,
        twist_body: Twist,
    ) -> Result<RigidBodyState, PrincipalFrameError> {
        let w = twist_body.angular;
        let v = twist_body.linear;
        if ![w.x, w.y, w.z, v.x, v.y, v.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(PrincipalFrameError::Invalid("body twist must be finite"));
        }
        let world_from_principal = world_from_body.compose(self.body_from_principal)?;
        let q = *world_from_principal.rotation().as_quat();
        let pose = Pose::new(
            mbd(world_from_principal.translation()),
            UnitQuaternion::new(q.w, q.x, q.y, q.z)?,
        )?;
        let center_velocity = world_from_body
            .rotation()
            .rotate(v + w.cross(self.inertia.center_of_mass()))?;
        let omega = self.body_from_principal.rotation().inverse().rotate(w)?;
        let moments = self.properties.principal_inertia_body();
        Ok(RigidBodyState::new(
            pose,
            mbd(center_velocity.scale(self.properties.mass())),
            fs_mbd::Vec3::new(
                moments.x * omega.x,
                moments.y * omega.y,
                moments.z * omega.z,
            ),
        )?)
    }

    /// Recover the original geometry's world pose after the existing integrator
    /// advances the principal-frame state. The geometry does not jump to its
    /// mass centroid or rotate with an arbitrary eigenvector relabeling.
    pub fn world_from_body(&self, state: RigidBodyState) -> Result<Se3, PrincipalFrameError> {
        let pose = state.pose();
        let [w, x, y, z] = pose.orientation().components();
        let rotation = So3::try_from_quat(Quat { w, x, y, z })?;
        let world_from_principal = Se3::from_parts(rotation, ga(pose.position_world()))?;
        Ok(world_from_principal.compose(self.body_from_principal.inverse()?)?)
    }

    /// Convert an authored body-frame contact point to the centre-relative
    /// principal arm expected by point_kinematics, effective mass and impulse
    /// methods. No impulse, contact law or friction model is inferred here.
    pub fn arm_for_body_point(
        &self,
        point_body: Vec3,
    ) -> Result<fs_mbd::Vec3, PrincipalFrameError> {
        Ok(mbd(self
            .body_from_principal
            .rotation()
            .inverse()
            .rotate(point_body - self.inertia.center_of_mass())?))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PrincipalFrameError> {
    if cancelled() {
        Err(PrincipalFrameError::Cancelled)
    } else {
        Ok(())
    }
}
fn mbd(v: Vec3) -> fs_mbd::Vec3 {
    fs_mbd::Vec3::new(v.x, v.y, v.z)
}
fn ga(v: fs_mbd::Vec3) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}

// Matrix-to-quaternion boundary conversion. Select the largest quaternion
// coordinate to avoid a small divisor at a half-turn. Group admission remains
// with fs-ga, and the actual returned rotation is tensor-checked above.
fn rotation_from_matrix(matrix: Mat3) -> Result<So3, PrincipalFrameError> {
    let m = matrix.m;
    let candidates = [
        1.0 + m[0] + m[4] + m[8],
        1.0 + m[0] - m[4] - m[8],
        1.0 - m[0] + m[4] - m[8],
        1.0 - m[0] - m[4] + m[8],
    ];
    let mut index = 0;
    for i in 1..4 {
        if candidates[i] > candidates[index] {
            index = i;
        }
    }
    let s = 2.0 * candidates[index].sqrt();
    let [w, x, y, z] = match index {
        0 => [
            0.25 * s,
            (m[7] - m[5]) / s,
            (m[2] - m[6]) / s,
            (m[3] - m[1]) / s,
        ],
        1 => [
            (m[7] - m[5]) / s,
            0.25 * s,
            (m[1] + m[3]) / s,
            (m[2] + m[6]) / s,
        ],
        2 => [
            (m[2] - m[6]) / s,
            (m[1] + m[3]) / s,
            0.25 * s,
            (m[5] + m[7]) / s,
        ],
        _ => [
            (m[3] - m[1]) / s,
            (m[2] + m[6]) / s,
            (m[5] + m[7]) / s,
            0.25 * s,
        ],
    };
    Ok(So3::try_from_quat(Quat { w, x, y, z })?)
}
