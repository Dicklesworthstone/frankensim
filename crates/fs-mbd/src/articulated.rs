//! Spatial and articulated rigid-body dynamics over the `fs-ga` Lie owner.
//!
//! Poses, twists, wrenches, adjoints, and coadjoints remain owned by `fs-ga`.
//! This module adds their physical meaning: inertias, joint subspaces, a
//! topologically ordered body tree, forward kinematics, recursive Newton-Euler
//! inverse dynamics, and Featherstone's articulated-body forward dynamics for
//! prescribed-base and unconstrained free-flight boundaries.
//! Spatial coordinates use `[angular, linear]`; dual force coordinates use
//! `[torque, force]`.

use core::fmt;
use fs_ga::{GaError, Mat3, Mat6, Se3, Twist, Vec3, Wrench};

const SYMMETRY_TOLERANCE: f64 = 64.0 * f64::EPSILON;
const PHYSICAL_TOLERANCE: f64 = 256.0 * f64::EPSILON;
const ARTICULATED_PIVOT_TOLERANCE: f64 = 1.0e-14;
const FLOATING_BASE_SYMMETRY_TOLERANCE: f64 = 1024.0 * f64::EPSILON;
const FLOATING_BASE_PIVOT_TOLERANCE: f64 = 1024.0 * f64::EPSILON;
const FLOATING_BASE_CONDITION_LIMIT: f64 = 1.0e12;

/// Refusal channel for articulated model construction and evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum ArticulatedError {
    /// A named scalar or coordinate is NaN or infinite.
    NonFinite {
        /// Input or derived quantity that refused.
        field: &'static str,
        /// Scalar index in the documented coordinate ordering.
        index: usize,
    },
    /// A mass must be strictly positive.
    InvalidMass,
    /// The rotational inertia is not symmetric within the declared tolerance.
    NonSymmetricInertia {
        /// Largest mirrored-entry disagreement.
        defect: f64,
        /// Accepted absolute tolerance after scale normalization.
        tolerance: f64,
    },
    /// The inertia does not describe a physically consistent three-dimensional
    /// mass distribution.
    NonPhysicalInertia {
        /// Deterministic failing principal-minor or pseudo-inertia indicator.
        measure: f64,
    },
    /// A joint axis must be finite and nonzero.
    InvalidJointAxis,
    /// Joint lower/upper, speed, or effort bounds are invalid.
    InvalidJointLimits,
    /// A supplied joint coordinate violates its declared position limits.
    JointLimitViolation {
        /// Link carrying the joint.
        link: usize,
        /// Supplied generalized coordinate.
        value: f64,
        /// Inclusive lower limit.
        lower: f64,
        /// Inclusive upper limit.
        upper: f64,
    },
    /// A supplied joint speed exceeds its declared symmetric limit.
    JointVelocityLimitViolation {
        /// Link carrying the joint.
        link: usize,
        /// Supplied generalized speed.
        value: f64,
        /// Maximum admitted absolute speed.
        maximum: f64,
    },
    /// A supplied actuator effort exceeds its declared symmetric limit.
    JointEffortLimitViolation {
        /// Link carrying the joint.
        link: usize,
        /// Supplied generalized effort.
        value: f64,
        /// Maximum admitted absolute effort.
        maximum: f64,
    },
    /// An articulated model must contain at least one link.
    EmptyModel,
    /// A link name must be nonempty and unique.
    InvalidLinkName {
        /// Link with the invalid or repeated name.
        link: usize,
    },
    /// Link zero is the only root; every later parent must precede its child.
    InvalidParent {
        /// Link whose parent relation refused.
        link: usize,
        /// Supplied parent, or `None` for an unexpected additional root.
        parent: Option<usize>,
    },
    /// A generalized or per-link payload has the wrong length.
    LengthMismatch {
        /// Payload being checked.
        field: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        got: usize,
    },
    /// A one-degree-of-freedom articulated inertia became singular or
    /// non-positive.
    SingularArticulatedInertia {
        /// Link whose scalar pivot refused.
        link: usize,
        /// Observed motion-subspace pivot.
        pivot: f64,
    },
    /// A floating-base solve requires the model root to be rigidly attached to
    /// the six-DoF base frame; a scalar root joint would create a redundant
    /// generalized-coordinate gauge.
    FloatingBaseRootJointNotFixed,
    /// The fixed-size floating-base articulated inertia has no safely positive
    /// pivot after scale normalization.
    SingularFloatingBaseInertia {
        /// Zero-based pivot in `[angular, linear]` ordering.
        index: usize,
        /// Scale-normalized pivot value.
        pivot: f64,
    },
    /// The accumulated floating-base articulated inertia lost the symmetry
    /// required by its physical energy form beyond the rounding allowance.
    NonSymmetricFloatingBaseInertia {
        /// Largest mirrored-entry disagreement.
        defect: f64,
        /// Accepted absolute tolerance at the observed matrix scale.
        tolerance: f64,
    },
    /// The fixed-size floating-base articulated inertia is too poorly
    /// conditioned for a trustworthy acceleration solve.
    IllConditionedFloatingBaseInertia {
        /// Deterministic infinity-norm condition estimate.
        condition_estimate: f64,
        /// Maximum admitted estimate.
        limit: f64,
    },
    /// An `fs-ga` group operation refused.
    Geometry(GaError),
}

impl fmt::Display for ArticulatedError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite { field, index } => {
                write!(formatter, "{field} coordinate {index} must be finite")
            }
            Self::InvalidMass => {
                formatter.write_str("spatial inertia mass must be finite and positive")
            }
            Self::NonSymmetricInertia { defect, tolerance } => write!(
                formatter,
                "rotational inertia symmetry defect {defect} exceeds tolerance {tolerance}"
            ),
            Self::NonPhysicalInertia { measure } => write!(
                formatter,
                "rotational inertia is not physically consistent (indicator {measure})"
            ),
            Self::InvalidJointAxis => formatter.write_str("joint axis must be finite and nonzero"),
            Self::InvalidJointLimits => formatter.write_str(
                "joint limits require finite lower <= upper and positive finite speed/effort",
            ),
            Self::JointLimitViolation {
                link,
                value,
                lower,
                upper,
            } => write!(
                formatter,
                "joint coordinate {value} for link {link} is outside [{lower}, {upper}]"
            ),
            Self::JointVelocityLimitViolation {
                link,
                value,
                maximum,
            } => write!(
                formatter,
                "joint speed {value} for link {link} exceeds absolute limit {maximum}"
            ),
            Self::JointEffortLimitViolation {
                link,
                value,
                maximum,
            } => write!(
                formatter,
                "joint effort {value} for link {link} exceeds absolute limit {maximum}"
            ),
            Self::EmptyModel => formatter.write_str("an articulated model needs at least one link"),
            Self::InvalidLinkName { link } => {
                write!(formatter, "link {link} has an empty or duplicate name")
            }
            Self::InvalidParent { link, parent } => {
                write!(formatter, "link {link} has invalid parent {parent:?}")
            }
            Self::LengthMismatch {
                field,
                expected,
                got,
            } => write!(formatter, "{field} has length {got}, expected {expected}"),
            Self::SingularArticulatedInertia { link, pivot } => write!(
                formatter,
                "articulated inertia pivot for link {link} is singular or non-positive ({pivot})"
            ),
            Self::FloatingBaseRootJointNotFixed => formatter.write_str(
                "a floating-base model root must use a fixed joint to avoid a redundant gauge",
            ),
            Self::SingularFloatingBaseInertia { index, pivot } => write!(
                formatter,
                "floating-base articulated inertia pivot {index} is singular or non-positive ({pivot})"
            ),
            Self::NonSymmetricFloatingBaseInertia { defect, tolerance } => write!(
                formatter,
                "floating-base articulated inertia symmetry defect {defect} exceeds {tolerance}"
            ),
            Self::IllConditionedFloatingBaseInertia {
                condition_estimate,
                limit,
            } => write!(
                formatter,
                "floating-base articulated inertia condition estimate {condition_estimate} exceeds {limit}"
            ),
            Self::Geometry(source) => write!(formatter, "Lie-group operation refused: {source}"),
        }
    }
}

impl std::error::Error for ArticulatedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Geometry(source) => Some(source),
            _ => None,
        }
    }
}

impl From<GaError> for ArticulatedError {
    fn from(value: GaError) -> Self {
        Self::Geometry(value)
    }
}

/// General rigid-body inertia expressed at an arbitrary body-frame origin.
///
/// `inertia_com` is the symmetric 3×3 rotational inertia about the centre of
/// mass, expressed in the same body frame. Physical consistency is checked via
/// positive rotational inertia and the positive-semidefinite covariance
/// `0.5 tr(I) I - I`, which is equivalent to the principal-moment triangle
/// inequalities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialInertia {
    mass: f64,
    center_of_mass: Vec3,
    inertia_com: Mat3,
}

impl SpatialInertia {
    /// Validate a mass, centre of mass, and full rotational inertia.
    pub fn new(
        mass: f64,
        center_of_mass: Vec3,
        inertia_com: Mat3,
    ) -> Result<Self, ArticulatedError> {
        if !mass.is_finite() || mass <= 0.0 {
            return Err(ArticulatedError::InvalidMass);
        }
        validate_vec3(center_of_mass, "center_of_mass")?;
        validate_slice(&inertia_com.m, "inertia_com")?;

        let scale = inertia_com.max_abs();
        if scale == 0.0 {
            return Err(ArticulatedError::NonPhysicalInertia { measure: 0.0 });
        }
        let tolerance = SYMMETRY_TOLERANCE * scale;
        let symmetry_defect = (inertia_com.m[1] - inertia_com.m[3])
            .abs()
            .max((inertia_com.m[2] - inertia_com.m[6]).abs())
            .max((inertia_com.m[5] - inertia_com.m[7]).abs());
        if symmetry_defect > tolerance {
            return Err(ArticulatedError::NonSymmetricInertia {
                defect: symmetry_defect,
                tolerance,
            });
        }

        let symmetric = symmetrize(inertia_com);
        let normalized = divide_mat3(symmetric, scale);
        if !is_positive_definite(normalized, PHYSICAL_TOLERANCE) {
            return Err(ArticulatedError::NonPhysicalInertia {
                measure: minimum_leading_principal_minor(normalized),
            });
        }
        let trace = normalized.m[0] + normalized.m[4] + normalized.m[8];
        let pseudo = Mat3 {
            m: [
                0.5 * trace - normalized.m[0],
                -normalized.m[1],
                -normalized.m[2],
                -normalized.m[3],
                0.5 * trace - normalized.m[4],
                -normalized.m[5],
                -normalized.m[6],
                -normalized.m[7],
                0.5 * trace - normalized.m[8],
            ],
        };
        if !is_positive_semidefinite(pseudo, PHYSICAL_TOLERANCE) {
            return Err(ArticulatedError::NonPhysicalInertia {
                measure: minimum_principal_minor(pseudo),
            });
        }

        Ok(Self {
            mass,
            center_of_mass,
            inertia_com: symmetric,
        })
    }

    /// Mass in kilograms.
    #[must_use]
    pub const fn mass(self) -> f64 {
        self.mass
    }

    /// Centre-of-mass offset from the body-frame origin.
    #[must_use]
    pub const fn center_of_mass(self) -> Vec3 {
        self.center_of_mass
    }

    /// Rotational inertia about the centre of mass.
    #[must_use]
    pub const fn inertia_com(self) -> Mat3 {
        self.inertia_com
    }

    /// Full 6×6 spatial inertia in `[angular, linear]` ordering.
    ///
    /// # Errors
    /// Refuses if finite admitted properties overflow while assembling the
    /// offset-reference matrix.
    pub fn matrix(self) -> Result<Mat6, ArticulatedError> {
        let cross = hat(self.center_of_mass);
        let cross_squared = cross.compose(cross);
        let mut matrix = Mat6::zero();
        let mut row = 0;
        while row < 3 {
            let mut column = 0;
            while column < 3 {
                let index = row * 3 + column;
                matrix.m[row * 6 + column] =
                    self.inertia_com.m[index] - self.mass * cross_squared.m[index];
                matrix.m[row * 6 + column + 3] = self.mass * cross.m[index];
                matrix.m[(row + 3) * 6 + column] = -self.mass * cross.m[index];
                matrix.m[(row + 3) * 6 + column + 3] = if row == column { self.mass } else { 0.0 };
                column += 1;
            }
            row += 1;
        }
        validate_mat6(&matrix, "spatial_inertia.matrix")?;
        Ok(matrix)
    }

    /// Spatial momentum `I * velocity` as a dual wrench.
    ///
    /// # Errors
    /// Refuses non-finite velocity or an unrepresentable derived momentum.
    pub fn momentum(self, velocity: Twist) -> Result<Wrench, ArticulatedError> {
        validate_twist(velocity, "spatial_inertia.velocity")?;
        let center_velocity = velocity.linear + velocity.angular.cross(self.center_of_mass);
        let linear_momentum = center_velocity.scale(self.mass);
        let angular_momentum =
            self.inertia_com.apply(velocity.angular) + self.center_of_mass.cross(linear_momentum);
        let momentum = Wrench::new(angular_momentum, linear_momentum);
        validate_wrench(momentum, "spatial_inertia.momentum")?;
        Ok(momentum)
    }

    /// Kinetic energy `0.5 * velocityᵀ I velocity`.
    ///
    /// # Errors
    /// Refuses non-finite velocity or an unrepresentable derived energy.
    pub fn kinetic_energy(self, velocity: Twist) -> Result<f64, ArticulatedError> {
        let energy = 0.5 * self.momentum(velocity)?.pairing(velocity);
        ensure_finite_scalar(energy, "spatial_inertia.kinetic_energy", 0)?;
        Ok(energy)
    }
}

/// Inclusive position, symmetric speed, and symmetric actuator-effort limits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JointLimits {
    /// Minimum generalized position.
    lower: f64,
    /// Maximum generalized position.
    upper: f64,
    /// Maximum absolute generalized speed.
    velocity: f64,
    /// Maximum absolute generalized effort.
    effort: f64,
}

impl JointLimits {
    /// Validate scalar joint limits.
    pub fn new(
        lower: f64,
        upper: f64,
        velocity: f64,
        effort: f64,
    ) -> Result<Self, ArticulatedError> {
        if !lower.is_finite()
            || !upper.is_finite()
            || lower > upper
            || !velocity.is_finite()
            || velocity <= 0.0
            || !effort.is_finite()
            || effort <= 0.0
        {
            return Err(ArticulatedError::InvalidJointLimits);
        }
        Ok(Self {
            lower,
            upper,
            velocity,
            effort,
        })
    }

    /// Inclusive minimum generalized position.
    #[must_use]
    pub const fn lower(self) -> f64 {
        self.lower
    }

    /// Inclusive maximum generalized position.
    #[must_use]
    pub const fn upper(self) -> f64 {
        self.upper
    }

    /// Maximum absolute generalized speed.
    #[must_use]
    pub const fn velocity(self) -> f64 {
        self.velocity
    }

    /// Maximum absolute generalized effort.
    #[must_use]
    pub const fn effort(self) -> f64 {
        self.effort
    }
}

/// Joint category exposed without making validated axes forgeable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JointKind {
    /// No relative motion.
    Fixed,
    /// One rotational coordinate.
    Revolute,
    /// One translational coordinate.
    Prismatic,
    /// Coupled rotation and axial translation.
    Helical,
}

/// Validated joint motion subspace.
///
/// The articulated solver deliberately specializes its hot path to zero- and
/// one-degree-of-freedom joints. The frame above the root has a prescribed
/// `SE(3)` pose and body twist; [`free_floating_forward_dynamics`] solves its
/// physical six-DoF acceleration as a distinct boundary rather than creating
/// six fake scalar joints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JointModel {
    kind: JointKind,
    motion_subspace: Twist,
    limits: Option<JointLimits>,
}

impl JointModel {
    /// Fixed joint with no generalized coordinate.
    pub const FIXED: Self = Self {
        kind: JointKind::Fixed,
        motion_subspace: Twist::zero(),
        limits: None,
    };

    /// Construct a validated revolute joint.
    pub fn revolute(axis: Vec3, limits: Option<JointLimits>) -> Result<Self, ArticulatedError> {
        let axis = normalized_axis(axis)?;
        Ok(Self {
            kind: JointKind::Revolute,
            motion_subspace: Twist::new(axis, Vec3::new(0.0, 0.0, 0.0)),
            limits,
        })
    }

    /// Construct a validated prismatic joint.
    pub fn prismatic(axis: Vec3, limits: Option<JointLimits>) -> Result<Self, ArticulatedError> {
        let axis = normalized_axis(axis)?;
        Ok(Self {
            kind: JointKind::Prismatic,
            motion_subspace: Twist::new(Vec3::new(0.0, 0.0, 0.0), axis),
            limits,
        })
    }

    /// Construct a validated helical joint.
    pub fn helical(
        axis: Vec3,
        pitch: f64,
        limits: Option<JointLimits>,
    ) -> Result<Self, ArticulatedError> {
        if !pitch.is_finite() {
            return Err(ArticulatedError::NonFinite {
                field: "joint.pitch",
                index: 0,
            });
        }
        let axis = normalized_axis(axis)?;
        Ok(Self {
            kind: JointKind::Helical,
            motion_subspace: Twist::new(axis, axis.scale(pitch)),
            limits,
        })
    }

    /// Joint category.
    #[must_use]
    pub const fn kind(self) -> JointKind {
        self.kind
    }

    /// Number of scalar generalized coordinates owned by this joint.
    #[must_use]
    pub const fn dof_count(self) -> usize {
        match self.kind {
            JointKind::Fixed => 0,
            JointKind::Revolute | JointKind::Prismatic | JointKind::Helical => 1,
        }
    }

    /// Constant body-coordinate motion subspace.
    #[must_use]
    pub const fn motion_subspace(self) -> Twist {
        self.motion_subspace
    }

    fn limits(self) -> Option<JointLimits> {
        self.limits
    }

    fn motion(self, coordinate: f64) -> Result<Se3, ArticulatedError> {
        if !coordinate.is_finite() {
            return Err(ArticulatedError::NonFinite {
                field: "joint.coordinate",
                index: 0,
            });
        }
        Se3::exp(self.motion_subspace.scale(coordinate)).map_err(Into::into)
    }
}

/// One rigid link and the joint connecting it to its parent.
#[derive(Debug, Clone, PartialEq)]
pub struct Link {
    name: String,
    parent: Option<usize>,
    parent_from_child_zero: Se3,
    joint: JointModel,
    inertia: SpatialInertia,
}

impl Link {
    /// Construct a link. Tree ordering and name uniqueness are checked by
    /// [`ArticulatedModel::new`].
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        parent: Option<usize>,
        parent_from_child_zero: Se3,
        joint: JointModel,
        inertia: SpatialInertia,
    ) -> Self {
        Self {
            name: name.into(),
            parent,
            parent_from_child_zero,
            joint,
            inertia,
        }
    }

    /// Stable human-readable link name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Parent link, or `None` for link zero.
    #[must_use]
    pub const fn parent(&self) -> Option<usize> {
        self.parent
    }

    /// Zero-coordinate transform mapping child coordinates to parent.
    #[must_use]
    pub const fn parent_from_child_zero(&self) -> Se3 {
        self.parent_from_child_zero
    }

    /// Joint motion model.
    #[must_use]
    pub const fn joint(&self) -> JointModel {
        self.joint
    }

    /// Link spatial inertia in child coordinates.
    #[must_use]
    pub const fn inertia(&self) -> SpatialInertia {
        self.inertia
    }
}

/// Validated topologically ordered articulated body tree.
#[derive(Debug, Clone, PartialEq)]
pub struct ArticulatedModel {
    links: Vec<Link>,
    dof_index: Vec<Option<usize>>,
    dof_count: usize,
}

impl ArticulatedModel {
    /// Validate a root-first body tree and assign compact generalized indices.
    pub fn new(links: Vec<Link>) -> Result<Self, ArticulatedError> {
        if links.is_empty() {
            return Err(ArticulatedError::EmptyModel);
        }
        let mut dof_index = Vec::with_capacity(links.len());
        let mut dof_count = 0usize;
        for (index, link) in links.iter().enumerate() {
            if link.name.trim().is_empty()
                || links[..index]
                    .iter()
                    .any(|preceding| preceding.name == link.name)
            {
                return Err(ArticulatedError::InvalidLinkName { link: index });
            }
            let parent_is_valid = if index == 0 {
                link.parent.is_none()
            } else {
                link.parent.is_some_and(|parent| parent < index)
            };
            if !parent_is_valid {
                return Err(ArticulatedError::InvalidParent {
                    link: index,
                    parent: link.parent,
                });
            }
            if link.joint.dof_count() == 1 {
                dof_index.push(Some(dof_count));
                dof_count += 1;
            } else {
                dof_index.push(None);
            }
        }
        Ok(Self {
            links,
            dof_index,
            dof_count,
        })
    }

    /// Links in root-first topological order.
    #[must_use]
    pub fn links(&self) -> &[Link] {
        &self.links
    }

    /// Number of links, including fixed links.
    #[must_use]
    pub fn link_count(&self) -> usize {
        self.links.len()
    }

    /// Number of scalar generalized coordinates.
    #[must_use]
    pub const fn dof_count(&self) -> usize {
        self.dof_count
    }

    /// Generalized-coordinate index for one link, or `None` for a fixed joint.
    #[must_use]
    pub fn dof_index(&self, link: usize) -> Option<usize> {
        self.dof_index.get(link).copied().flatten()
    }

    /// Linear working-set claim for the articulated-body implementation.
    #[must_use]
    pub fn complexity(&self) -> DynamicsComplexity {
        DynamicsComplexity {
            links: self.links.len(),
            degrees_of_freedom: self.dof_count,
            dense_generalized_matrix_entries: 0,
            spatial_matrix_entries_per_link: 36,
        }
    }

    /// Linear working-set claim for the unconstrained free-floating solver.
    pub fn free_floating_complexity(
        &self,
    ) -> Result<FreeFloatingDynamicsComplexity, ArticulatedError> {
        if self.links[0].joint.kind() != JointKind::Fixed {
            return Err(ArticulatedError::FloatingBaseRootJointNotFixed);
        }
        Ok(FreeFloatingDynamicsComplexity {
            tree: self.complexity(),
            base_degrees_of_freedom: 6,
            fixed_root_solve_matrix_entries: 36,
        })
    }

    fn validate_generalized(
        &self,
        field: &'static str,
        values: &[f64],
    ) -> Result<(), ArticulatedError> {
        require_len(field, values.len(), self.dof_count)?;
        validate_slice(values, field)
    }

    fn validate_configuration(&self, position: &[f64]) -> Result<(), ArticulatedError> {
        self.validate_generalized("joint_position", position)?;
        for (link_index, link) in self.links.iter().enumerate() {
            if let Some(dof) = self.dof_index[link_index]
                && let Some(limits) = link.joint.limits()
            {
                let value = position[dof];
                if value < limits.lower() || value > limits.upper() {
                    return Err(ArticulatedError::JointLimitViolation {
                        link: link_index,
                        value,
                        lower: limits.lower(),
                        upper: limits.upper(),
                    });
                }
            }
        }
        Ok(())
    }

    fn validate_velocity(&self, velocity: &[f64]) -> Result<(), ArticulatedError> {
        self.validate_generalized("joint_velocity", velocity)?;
        for (link_index, link) in self.links.iter().enumerate() {
            if let Some(dof) = self.dof_index[link_index]
                && let Some(limits) = link.joint.limits()
            {
                let value = velocity[dof];
                if value.abs() > limits.velocity() {
                    return Err(ArticulatedError::JointVelocityLimitViolation {
                        link: link_index,
                        value,
                        maximum: limits.velocity(),
                    });
                }
            }
        }
        Ok(())
    }

    fn validate_applied_effort(&self, effort: &[f64]) -> Result<(), ArticulatedError> {
        self.validate_generalized("generalized_force", effort)?;
        for (link_index, link) in self.links.iter().enumerate() {
            if let Some(dof) = self.dof_index[link_index]
                && let Some(limits) = link.joint.limits()
            {
                let value = effort[dof];
                if value.abs() > limits.effort() {
                    return Err(ArticulatedError::JointEffortLimitViolation {
                        link: link_index,
                        value,
                        maximum: limits.effort(),
                    });
                }
            }
        }
        Ok(())
    }
}

/// Complexity metadata surfaced to browser and optimizer policy layers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DynamicsComplexity {
    /// Number of links traversed by each recursion.
    pub links: usize,
    /// Number of scalar joints.
    pub degrees_of_freedom: usize,
    /// Dense generalized mass-matrix entries allocated by this lane (zero).
    pub dense_generalized_matrix_entries: usize,
    /// Fixed-size spatial inertia entries retained per link.
    pub spatial_matrix_entries_per_link: usize,
}

/// Complexity metadata for unconstrained free-floating articulated dynamics.
///
/// The only dense solve is one fixed-size 6×6 root system. No matrix whose
/// extent grows with the generalized-coordinate count is formed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FreeFloatingDynamicsComplexity {
    /// Shared linear tree-pass metadata.
    pub tree: DynamicsComplexity,
    /// Unactuated base coordinates solved at the root.
    pub base_degrees_of_freedom: usize,
    /// Entries in the single fixed-size root articulated-inertia solve.
    pub fixed_root_solve_matrix_entries: usize,
}

/// Prescribed motion of the frame above the model root.
///
/// This is a fixed or externally driven base boundary. The articulated solver
/// using this input does not infer the base acceleration or solve free-floating
/// equilibrium; use [`free_floating_forward_dynamics`] for that distinct
/// boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaseState {
    /// Pose mapping base coordinates into world coordinates.
    pub world_from_base: Se3,
    /// Base twist expressed in base coordinates.
    pub twist_body: Twist,
    /// Prescribed non-gravity base acceleration in base coordinates.
    ///
    /// Dynamics subtract world gravity, expressed in the base frame, from
    /// this value to form Featherstone's root acceleration boundary.
    pub acceleration_body: Twist,
}

impl BaseState {
    /// A stationary base at a declared world pose.
    #[must_use]
    pub fn stationary(world_from_base: Se3) -> Self {
        Self {
            world_from_base,
            twist_body: Twist::zero(),
            acceleration_body: Twist::zero(),
        }
    }

    /// Construct an externally prescribed base motion boundary.
    #[must_use]
    pub const fn prescribed(
        world_from_base: Se3,
        twist_body: Twist,
        acceleration_body: Twist,
    ) -> Self {
        Self {
            world_from_base,
            twist_body,
            acceleration_body,
        }
    }
}

/// Pose and body twist of an unconstrained six-DoF frame above the model root.
///
/// This record reuses the canonical `fs-ga` pose and twist owners. It contains
/// no acceleration because [`free_floating_forward_dynamics`] solves that
/// quantity from force balance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeFloatingBaseState {
    /// Pose mapping floating-base coordinates into world coordinates.
    pub world_from_base: Se3,
    /// Floating-base twist expressed in base coordinates.
    pub twist_body: Twist,
}

impl FreeFloatingBaseState {
    /// Construct a free-floating base state from canonical `fs-ga` values.
    #[must_use]
    pub const fn new(world_from_base: Se3, twist_body: Twist) -> Self {
        Self {
            world_from_base,
            twist_body,
        }
    }

    /// A stationary free-floating base at a declared world pose.
    #[must_use]
    pub fn stationary(world_from_base: Se3) -> Self {
        Self::new(world_from_base, Twist::zero())
    }
}

/// Complete forward-kinematics receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct Kinematics {
    /// World-from-link poses in model order.
    pub world_from_link: Vec<Se3>,
    /// Body-coordinate link twists in model order.
    pub body_twist: Vec<Twist>,
}

/// Recursive inverse-dynamics receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct InverseDynamics {
    /// Generalized forces in compact joint order.
    pub generalized_force: Vec<f64>,
    /// Body-coordinate link accelerations, including the gravity convention.
    pub body_acceleration: Vec<Twist>,
    /// Net body-coordinate wrenches propagated through each link.
    pub body_wrench: Vec<Wrench>,
}

/// Articulated-body forward-dynamics receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct ForwardDynamics {
    /// Generalized accelerations in compact joint order.
    pub generalized_acceleration: Vec<f64>,
    /// Body-coordinate link accelerations, including the gravity convention.
    pub body_acceleration: Vec<Twist>,
}

/// Unconstrained free-flight articulated-body forward-dynamics receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct FreeFloatingForwardDynamics {
    /// Featherstone spatial acceleration of the floating base, in base
    /// coordinates. Uniform gravity is included rather than hidden in an
    /// apparent-acceleration convention.
    ///
    /// Its linear component is not the ordinary Cartesian acceleration of the
    /// base origin when the base twist is nonzero. Use
    /// [`origin_linear_acceleration_body`] for that conversion.
    pub base_spatial_acceleration_body: Twist,
    /// Actuated generalized accelerations in compact scalar-joint order.
    pub generalized_acceleration: Vec<f64>,
    /// Featherstone spatial acceleration of every retained link, expressed in
    /// that link's body coordinates. Convert each linear component with
    /// [`origin_linear_acceleration_body`] and the corresponding body twist
    /// before integrating a Cartesian origin velocity.
    pub body_spatial_acceleration: Vec<Twist>,
}

/// Convert a body-coordinate Featherstone spatial acceleration into the
/// ordinary Cartesian acceleration of that frame's origin, expressed in the
/// same body coordinates.
///
/// For body twist `[omega, v]` and spatial acceleration `[alpha, a]`, the
/// origin acceleration is `a + omega x v`. This explicit boundary prevents a
/// time integrator from treating the spatial acceleration's linear coordinate
/// as an ordinary vector derivative.
pub fn origin_linear_acceleration_body(
    spatial_acceleration_body: Twist,
    body_twist: Twist,
) -> Result<Vec3, ArticulatedError> {
    validate_twist(
        spatial_acceleration_body,
        "origin_linear_acceleration.spatial_acceleration_body",
    )?;
    validate_twist(body_twist, "origin_linear_acceleration.body_twist")?;
    let acceleration =
        spatial_acceleration_body.linear + body_twist.angular.cross(body_twist.linear);
    validate_vec3(acceleration, "origin_linear_acceleration.result")?;
    Ok(acceleration)
}

/// Compute root-first poses and body-coordinate twists.
pub fn forward_kinematics(
    model: &ArticulatedModel,
    base: BaseState,
    joint_position: &[f64],
    joint_velocity: &[f64],
) -> Result<Kinematics, ArticulatedError> {
    model.validate_configuration(joint_position)?;
    model.validate_velocity(joint_velocity)?;
    validate_twist(base.twist_body, "base.twist_body")?;
    validate_twist(base.acceleration_body, "base.acceleration_body")?;

    let mut poses = Vec::with_capacity(model.link_count());
    let mut velocities = Vec::with_capacity(model.link_count());
    for (link_index, link) in model.links.iter().enumerate() {
        let coordinate = model.dof_index[link_index].map_or(0.0, |dof| joint_position[dof]);
        let rate = model.dof_index[link_index].map_or(0.0, |dof| joint_velocity[dof]);
        let parent_from_child = link
            .parent_from_child_zero
            .compose(link.joint.motion(coordinate)?)?;
        let motion_parent_to_child = parent_from_child.inverse()?.adjoint();
        let joint_velocity_body = link.joint.motion_subspace().scale(rate);

        let (world_from_parent, parent_velocity) = match link.parent {
            Some(parent) => (poses[parent], velocities[parent]),
            None => (base.world_from_base, base.twist_body),
        };
        poses.push(world_from_parent.compose(parent_from_child)?);
        let velocity = motion_parent_to_child
            .apply_twist(parent_velocity)
            .plus(joint_velocity_body);
        validate_twist(velocity, "forward_kinematics.body_twist")?;
        velocities.push(velocity);
    }
    Ok(Kinematics {
        world_from_link: poses,
        body_twist: velocities,
    })
}

/// Recursive Newton-Euler inverse dynamics for zero- and one-DoF trees.
///
/// `external_body_wrench` is expressed in each link's body coordinates and is
/// subtracted from the wrench the actuators must supply. Gravity is supplied in
/// world coordinates. The base pose and twist are known; floating-base
/// equilibrium and non-gravity base acceleration are intentionally separate
/// solve inputs.
pub fn inverse_dynamics(
    model: &ArticulatedModel,
    base: BaseState,
    joint_position: &[f64],
    joint_velocity: &[f64],
    joint_acceleration: &[f64],
    gravity_world: Vec3,
    external_body_wrench: &[Wrench],
) -> Result<InverseDynamics, ArticulatedError> {
    model.validate_configuration(joint_position)?;
    model.validate_velocity(joint_velocity)?;
    model.validate_generalized("joint_acceleration", joint_acceleration)?;
    require_len(
        "external_body_wrench",
        external_body_wrench.len(),
        model.link_count(),
    )?;
    validate_vec3(gravity_world, "gravity_world")?;
    validate_twist(base.twist_body, "base.twist_body")?;
    validate_twist(base.acceleration_body, "base.acceleration_body")?;
    for wrench in external_body_wrench {
        validate_wrench(*wrench, "external_body_wrench")?;
    }

    let link_count = model.link_count();
    let mut transforms = vec![Mat6::zero(); link_count];
    let mut velocity = vec![Twist::zero(); link_count];
    let mut acceleration = vec![Twist::zero(); link_count];
    let mut wrench =
        vec![Wrench::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)); link_count];

    let gravity_base = base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(gravity_world)?;
    let base_acceleration = base.acceleration_body.plus(Twist::new(
        Vec3::new(0.0, 0.0, 0.0),
        gravity_base.scale(-1.0),
    ));
    validate_twist(base_acceleration, "inverse_dynamics.base_acceleration")?;

    for (link_index, link) in model.links.iter().enumerate() {
        let dof = model.dof_index[link_index];
        let coordinate = dof.map_or(0.0, |index| joint_position[index]);
        let rate = dof.map_or(0.0, |index| joint_velocity[index]);
        let rate_dot = dof.map_or(0.0, |index| joint_acceleration[index]);
        let parent_from_child = link
            .parent_from_child_zero
            .compose(link.joint.motion(coordinate)?)?;
        let x_up = parent_from_child.inverse()?.adjoint();
        transforms[link_index] = x_up;
        let subspace = link.joint.motion_subspace();
        let joint_velocity_body = subspace.scale(rate);
        let (parent_velocity, parent_acceleration) = match link.parent {
            Some(parent) => (velocity[parent], acceleration[parent]),
            None => (base.twist_body, base_acceleration),
        };
        velocity[link_index] = x_up.apply_twist(parent_velocity).plus(joint_velocity_body);
        acceleration[link_index] = x_up
            .apply_twist(parent_acceleration)
            .plus(subspace.scale(rate_dot))
            .plus(velocity[link_index].bracket(joint_velocity_body));
        validate_twist(velocity[link_index], "inverse_dynamics.body_velocity")?;
        validate_twist(
            acceleration[link_index],
            "inverse_dynamics.body_acceleration",
        )?;
        let momentum = link.inertia.momentum(velocity[link_index])?;
        wrench[link_index] = wrench_minus(
            wrench_plus(
                link.inertia.momentum(acceleration[link_index])?,
                cross_force(velocity[link_index], momentum),
            ),
            external_body_wrench[link_index],
        );
        validate_wrench(wrench[link_index], "inverse_dynamics.body_wrench")?;
    }

    let mut generalized_force = vec![0.0; model.dof_count];
    for link_index in (0..link_count).rev() {
        if let Some(dof) = model.dof_index[link_index] {
            generalized_force[dof] =
                wrench[link_index].pairing(model.links[link_index].joint.motion_subspace());
        }
        if let Some(parent) = model.links[link_index].parent {
            let parent_wrench = transforms[link_index]
                .transpose()
                .apply_wrench(wrench[link_index]);
            wrench[parent] = wrench_plus(wrench[parent], parent_wrench);
            validate_wrench(wrench[parent], "inverse_dynamics.propagated_wrench")?;
        }
    }
    validate_slice(&generalized_force, "inverse_dynamics.generalized_force")?;
    Ok(InverseDynamics {
        generalized_force,
        body_acceleration: acceleration,
        body_wrench: wrench,
    })
}

/// Featherstone articulated-body forward dynamics in linear time and storage.
///
/// This implementation never forms or factors a dense generalized mass
/// matrix. Every link contributes one fixed 6×6 articulated inertia and a
/// bounded set of six-vectors.
pub fn forward_dynamics(
    model: &ArticulatedModel,
    base: BaseState,
    joint_position: &[f64],
    joint_velocity: &[f64],
    generalized_force: &[f64],
    gravity_world: Vec3,
    external_body_wrench: &[Wrench],
) -> Result<ForwardDynamics, ArticulatedError> {
    validate_forward_dynamics_header(
        model,
        base.twist_body,
        joint_position,
        joint_velocity,
        generalized_force,
        gravity_world,
        external_body_wrench,
    )?;
    validate_twist(base.acceleration_body, "base.acceleration_body")?;
    validate_external_wrenches(external_body_wrench)?;
    let pass = articulated_body_pass(
        model,
        base.twist_body,
        joint_position,
        joint_velocity,
        generalized_force,
        external_body_wrench,
        None,
        false,
        PRESCRIBED_ABA_FIELDS,
    )?;

    let gravity_base = base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(gravity_world)?;
    let base_acceleration = base.acceleration_body.plus(Twist::new(
        Vec3::new(0.0, 0.0, 0.0),
        gravity_base.scale(-1.0),
    ));
    validate_twist(base_acceleration, "forward_dynamics.base_acceleration")?;
    let (generalized_acceleration, body_acceleration) = propagate_aba_accelerations(
        model,
        &pass,
        base_acceleration,
        "forward_dynamics.generalized_acceleration",
        "forward_dynamics.body_acceleration",
    )?;
    Ok(ForwardDynamics {
        generalized_acceleration,
        body_acceleration,
    })
}

/// Featherstone articulated-body forward dynamics with a physical free base.
///
/// The six base coordinates are unactuated and solved together with the scalar
/// joint accelerations. Gravity is a physical uniform body force in this API,
/// so the returned base and link accelerations are inertial-frame physical
/// accelerations expressed in their respective body coordinates.
///
/// This is unconstrained free flight. It does not impose support, contact,
/// impact, friction, or ground reaction forces, and it does not synthesize
/// controller effort or advance the state in time.
///
/// # Errors
/// Refuses invalid dimensions, limits, non-finite input or derived state, a
/// non-fixed root joint, scalar articulated-inertia singularity, or a singular
/// or ill-conditioned fixed-size root system.
pub fn free_floating_forward_dynamics(
    model: &ArticulatedModel,
    base: FreeFloatingBaseState,
    joint_position: &[f64],
    joint_velocity: &[f64],
    generalized_force: &[f64],
    gravity_world: Vec3,
    external_body_wrench: &[Wrench],
) -> Result<FreeFloatingForwardDynamics, ArticulatedError> {
    validate_forward_dynamics_header(
        model,
        base.twist_body,
        joint_position,
        joint_velocity,
        generalized_force,
        gravity_world,
        external_body_wrench,
    )?;
    validate_external_wrenches(external_body_wrench)?;
    if model.links[0].joint.kind() != JointKind::Fixed {
        return Err(ArticulatedError::FloatingBaseRootJointNotFixed);
    }

    let gravity_base = base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(gravity_world)?;
    let gravity_acceleration_base = Twist::new(Vec3::new(0.0, 0.0, 0.0), gravity_base);
    let pass = articulated_body_pass(
        model,
        base.twist_body,
        joint_position,
        joint_velocity,
        generalized_force,
        external_body_wrench,
        Some(gravity_acceleration_base),
        true,
        FREE_FLOATING_ABA_FIELDS,
    )?;
    let base_spatial_acceleration_body = solve_floating_base_system(
        &pass.base_articulated_inertia,
        wrench_from_array(scale6(pass.base_articulated_bias.to_array(), -1.0)),
    )?;
    validate_twist(
        base_spatial_acceleration_body,
        "free_floating_forward_dynamics.base_acceleration",
    )?;
    let (generalized_acceleration, body_spatial_acceleration) = propagate_aba_accelerations(
        model,
        &pass,
        base_spatial_acceleration_body,
        "free_floating_forward_dynamics.generalized_acceleration",
        "free_floating_forward_dynamics.body_acceleration",
    )?;
    Ok(FreeFloatingForwardDynamics {
        base_spatial_acceleration_body,
        generalized_acceleration,
        body_spatial_acceleration,
    })
}

struct ArticulatedBodyPass {
    x_up: Vec<Mat6>,
    bias_acceleration: Vec<Twist>,
    u: Vec<[f64; 6]>,
    d: Vec<f64>,
    scalar_u: Vec<f64>,
    base_articulated_inertia: Mat6,
    base_articulated_bias: Wrench,
}

#[derive(Clone, Copy)]
struct AbaFieldNames {
    body_velocity: &'static str,
    bias_acceleration: &'static str,
    articulated_bias: &'static str,
    reduced_inertia: &'static str,
    reduced_bias: &'static str,
    propagated_inertia: &'static str,
    propagated_bias: &'static str,
}

const PRESCRIBED_ABA_FIELDS: AbaFieldNames = AbaFieldNames {
    body_velocity: "forward_dynamics.body_velocity",
    bias_acceleration: "forward_dynamics.bias_acceleration",
    articulated_bias: "forward_dynamics.articulated_bias",
    reduced_inertia: "forward_dynamics.reduced_inertia",
    reduced_bias: "forward_dynamics.reduced_bias",
    propagated_inertia: "forward_dynamics.propagated_inertia",
    propagated_bias: "forward_dynamics.propagated_bias",
};

const FREE_FLOATING_ABA_FIELDS: AbaFieldNames = AbaFieldNames {
    body_velocity: "free_floating_forward_dynamics.body_velocity",
    bias_acceleration: "free_floating_forward_dynamics.bias_acceleration",
    articulated_bias: "free_floating_forward_dynamics.articulated_bias",
    reduced_inertia: "free_floating_forward_dynamics.reduced_inertia",
    reduced_bias: "free_floating_forward_dynamics.reduced_bias",
    propagated_inertia: "free_floating_forward_dynamics.propagated_inertia",
    propagated_bias: "free_floating_forward_dynamics.propagated_bias",
};

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn articulated_body_pass(
    model: &ArticulatedModel,
    base_twist_body: Twist,
    joint_position: &[f64],
    joint_velocity: &[f64],
    generalized_force: &[f64],
    external_body_wrench: &[Wrench],
    gravity_acceleration_base: Option<Twist>,
    collect_base_system: bool,
    fields: AbaFieldNames,
) -> Result<ArticulatedBodyPass, ArticulatedError> {
    let link_count = model.link_count();
    let zero_wrench = Wrench::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
    let mut x_up = vec![Mat6::zero(); link_count];
    let mut velocity = vec![Twist::zero(); link_count];
    let mut bias_acceleration = vec![Twist::zero(); link_count];
    let mut articulated_inertia = vec![Mat6::zero(); link_count];
    let mut articulated_bias = vec![zero_wrench; link_count];
    let mut gravity_acceleration =
        gravity_acceleration_base.map(|_| vec![Twist::zero(); link_count]);
    let mut u = vec![[0.0; 6]; link_count];
    let mut d = vec![0.0; link_count];
    let mut scalar_u = vec![0.0; link_count];

    for (link_index, link) in model.links.iter().enumerate() {
        let dof = model.dof_index[link_index];
        let coordinate = dof.map_or(0.0, |index| joint_position[index]);
        let rate = dof.map_or(0.0, |index| joint_velocity[index]);
        let parent_from_child = link
            .parent_from_child_zero
            .compose(link.joint.motion(coordinate)?)?;
        x_up[link_index] = parent_from_child.inverse()?.adjoint();
        let joint_velocity_body = link.joint.motion_subspace().scale(rate);
        let parent_velocity = link
            .parent
            .map_or(base_twist_body, |parent| velocity[parent]);
        velocity[link_index] = x_up[link_index]
            .apply_twist(parent_velocity)
            .plus(joint_velocity_body);
        bias_acceleration[link_index] = velocity[link_index].bracket(joint_velocity_body);
        validate_twist(velocity[link_index], fields.body_velocity)?;
        validate_twist(bias_acceleration[link_index], fields.bias_acceleration)?;
        articulated_inertia[link_index] = link.inertia.matrix()?;
        let velocity_bias = cross_force(
            velocity[link_index],
            link.inertia.momentum(velocity[link_index])?,
        );
        let bias_without_gravity = wrench_minus(velocity_bias, external_body_wrench[link_index]);
        articulated_bias[link_index] =
            match (gravity_acceleration_base, gravity_acceleration.as_mut()) {
                (Some(base_gravity), Some(link_gravity)) => {
                    let parent_gravity = link
                        .parent
                        .map_or(base_gravity, |parent| link_gravity[parent]);
                    link_gravity[link_index] = x_up[link_index].apply_twist(parent_gravity);
                    validate_twist(
                        link_gravity[link_index],
                        "free_floating_forward_dynamics.gravity_acceleration",
                    )?;
                    wrench_minus(
                        bias_without_gravity,
                        link.inertia.momentum(link_gravity[link_index])?,
                    )
                }
                _ => bias_without_gravity,
            };
        validate_wrench(articulated_bias[link_index], fields.articulated_bias)?;
    }

    let mut base_articulated_inertia = Mat6::zero();
    let mut base_articulated_bias = zero_wrench;
    for link_index in (0..link_count).rev() {
        let mut reduced_inertia = articulated_inertia[link_index];
        let mut reduced_bias = wrench_plus(
            articulated_bias[link_index],
            mat6_apply_twist(&reduced_inertia, bias_acceleration[link_index]),
        );
        if let Some(dof_index) = model.dof_index[link_index] {
            let subspace = model.links[link_index].joint.motion_subspace();
            u[link_index] = mat6_apply_array(&reduced_inertia, subspace.to_array());
            d[link_index] = dot6(subspace.to_array(), u[link_index]);
            let pivot_scale = dot_abs6(subspace.to_array(), u[link_index]);
            if !d[link_index].is_finite()
                || !pivot_scale.is_finite()
                || d[link_index] <= ARTICULATED_PIVOT_TOLERANCE * pivot_scale
            {
                return Err(ArticulatedError::SingularArticulatedInertia {
                    link: link_index,
                    pivot: d[link_index],
                });
            }
            scalar_u[link_index] =
                generalized_force[dof_index] - articulated_bias[link_index].pairing(subspace);
            reduced_inertia = mat6_sub_outer_scaled(
                &reduced_inertia,
                u[link_index],
                u[link_index],
                d[link_index].recip(),
            );
            validate_mat6(&reduced_inertia, fields.reduced_inertia)?;
            reduced_bias = wrench_plus(
                wrench_plus(
                    articulated_bias[link_index],
                    mat6_apply_twist(&reduced_inertia, bias_acceleration[link_index]),
                ),
                wrench_from_array(scale6(u[link_index], scalar_u[link_index] / d[link_index])),
            );
            validate_wrench(reduced_bias, fields.reduced_bias)?;
        }

        if let Some(parent) = model.links[link_index].parent {
            let child_to_parent = x_up[link_index].transpose();
            let inertia_in_child = mat6_compose(&reduced_inertia, &x_up[link_index]);
            let transformed_inertia = mat6_compose(&child_to_parent, &inertia_in_child);
            articulated_inertia[parent] =
                mat6_add(&articulated_inertia[parent], &transformed_inertia);
            validate_mat6(&articulated_inertia[parent], fields.propagated_inertia)?;
            let transformed_bias = child_to_parent.apply_wrench(reduced_bias);
            articulated_bias[parent] = wrench_plus(articulated_bias[parent], transformed_bias);
            validate_wrench(articulated_bias[parent], fields.propagated_bias)?;
        } else if collect_base_system {
            let root_to_base = x_up[link_index].transpose();
            base_articulated_inertia = mat6_compose(
                &root_to_base,
                &mat6_compose(&reduced_inertia, &x_up[link_index]),
            );
            base_articulated_bias = root_to_base.apply_wrench(reduced_bias);
            validate_mat6(
                &base_articulated_inertia,
                "free_floating_forward_dynamics.base_articulated_inertia",
            )?;
            validate_wrench(
                base_articulated_bias,
                "free_floating_forward_dynamics.base_articulated_bias",
            )?;
        }
    }

    Ok(ArticulatedBodyPass {
        x_up,
        bias_acceleration,
        u,
        d,
        scalar_u,
        base_articulated_inertia,
        base_articulated_bias,
    })
}

fn propagate_aba_accelerations(
    model: &ArticulatedModel,
    pass: &ArticulatedBodyPass,
    base_acceleration: Twist,
    generalized_field: &'static str,
    body_field: &'static str,
) -> Result<(Vec<f64>, Vec<Twist>), ArticulatedError> {
    let mut acceleration = vec![Twist::zero(); model.link_count()];
    let mut generalized_acceleration = vec![0.0; model.dof_count];
    for link_index in 0..model.link_count() {
        let parent_acceleration = model.links[link_index]
            .parent
            .map_or(base_acceleration, |parent| acceleration[parent]);
        let mut link_acceleration = pass.x_up[link_index]
            .apply_twist(parent_acceleration)
            .plus(pass.bias_acceleration[link_index]);
        if let Some(dof_index) = model.dof_index[link_index] {
            let qdd = (pass.scalar_u[link_index]
                - dot6(pass.u[link_index], link_acceleration.to_array()))
                / pass.d[link_index];
            if !qdd.is_finite() {
                return Err(ArticulatedError::NonFinite {
                    field: generalized_field,
                    index: dof_index,
                });
            }
            generalized_acceleration[dof_index] = qdd;
            link_acceleration =
                link_acceleration.plus(model.links[link_index].joint.motion_subspace().scale(qdd));
        }
        validate_twist(link_acceleration, body_field)?;
        acceleration[link_index] = link_acceleration;
    }
    Ok((generalized_acceleration, acceleration))
}

fn validate_forward_dynamics_header(
    model: &ArticulatedModel,
    base_twist_body: Twist,
    joint_position: &[f64],
    joint_velocity: &[f64],
    generalized_force: &[f64],
    gravity_world: Vec3,
    external_body_wrench: &[Wrench],
) -> Result<(), ArticulatedError> {
    model.validate_configuration(joint_position)?;
    model.validate_velocity(joint_velocity)?;
    model.validate_applied_effort(generalized_force)?;
    require_len(
        "external_body_wrench",
        external_body_wrench.len(),
        model.link_count(),
    )?;
    validate_vec3(gravity_world, "gravity_world")?;
    validate_twist(base_twist_body, "base.twist_body")?;
    Ok(())
}

fn validate_external_wrenches(external_body_wrench: &[Wrench]) -> Result<(), ArticulatedError> {
    for wrench in external_body_wrench {
        validate_wrench(*wrench, "external_body_wrench")?;
    }
    Ok(())
}

fn solve_floating_base_system(
    matrix: &Mat6,
    right_hand_side: Wrench,
) -> Result<Twist, ArticulatedError> {
    validate_mat6(
        matrix,
        "free_floating_forward_dynamics.base_articulated_inertia",
    )?;
    validate_wrench(
        right_hand_side,
        "free_floating_forward_dynamics.base_right_hand_side",
    )?;
    let scale = matrix.max_abs();
    if scale == 0.0 {
        return Err(ArticulatedError::SingularFloatingBaseInertia {
            index: 0,
            pivot: 0.0,
        });
    }

    let normalized = normalize_symmetric6(matrix, scale)?;
    let lower = cholesky6(&normalized)?;
    let condition_estimate = cholesky_condition_estimate6(&normalized, &lower);
    if !condition_estimate.is_finite() || condition_estimate > FLOATING_BASE_CONDITION_LIMIT {
        return Err(ArticulatedError::IllConditionedFloatingBaseInertia {
            condition_estimate,
            limit: FLOATING_BASE_CONDITION_LIMIT,
        });
    }

    let mut normalized_right_hand_side = right_hand_side.to_array();
    for value in &mut normalized_right_hand_side {
        *value /= scale;
    }
    validate_slice(
        &normalized_right_hand_side,
        "free_floating_forward_dynamics.normalized_base_right_hand_side",
    )?;
    let solution = solve_cholesky6(&lower, normalized_right_hand_side);
    validate_slice(
        &solution,
        "free_floating_forward_dynamics.base_acceleration",
    )?;
    Ok(twist_from_array(solution))
}

fn normalize_symmetric6(matrix: &Mat6, scale: f64) -> Result<[f64; 36], ArticulatedError> {
    let symmetry_tolerance = FLOATING_BASE_SYMMETRY_TOLERANCE * scale;
    let mut symmetry_defect = 0.0_f64;
    let mut row = 0;
    while row < 6 {
        let mut column = row + 1;
        while column < 6 {
            symmetry_defect = symmetry_defect
                .max((matrix.m[row * 6 + column] - matrix.m[column * 6 + row]).abs());
            column += 1;
        }
        row += 1;
    }
    if symmetry_defect > symmetry_tolerance {
        return Err(ArticulatedError::NonSymmetricFloatingBaseInertia {
            defect: symmetry_defect,
            tolerance: symmetry_tolerance,
        });
    }

    let mut normalized = [0.0; 36];
    row = 0;
    while row < 6 {
        let mut column = 0;
        while column < 6 {
            normalized[row * 6 + column] =
                (0.5 * matrix.m[row * 6 + column] + 0.5 * matrix.m[column * 6 + row]) / scale;
            column += 1;
        }
        row += 1;
    }
    validate_slice(
        &normalized,
        "free_floating_forward_dynamics.normalized_base_inertia",
    )?;
    Ok(normalized)
}

fn cholesky6(normalized: &[f64; 36]) -> Result<[f64; 36], ArticulatedError> {
    let mut lower = [0.0_f64; 36];
    let mut row = 0;
    while row < 6 {
        let mut column = 0;
        while column <= row {
            let mut pivot = normalized[row * 6 + column];
            let mut inner = 0;
            while inner < column {
                pivot = (-lower[row * 6 + inner]).mul_add(lower[column * 6 + inner], pivot);
                inner += 1;
            }
            if row == column {
                if !pivot.is_finite() || pivot <= FLOATING_BASE_PIVOT_TOLERANCE {
                    return Err(ArticulatedError::SingularFloatingBaseInertia {
                        index: row,
                        pivot,
                    });
                }
                lower[row * 6 + column] = pivot.sqrt();
            } else {
                lower[row * 6 + column] = pivot / lower[column * 6 + column];
            }
            column += 1;
        }
        row += 1;
    }
    validate_slice(
        &lower,
        "free_floating_forward_dynamics.base_cholesky_factor",
    )?;
    Ok(lower)
}

fn cholesky_condition_estimate6(normalized: &[f64; 36], lower: &[f64; 36]) -> f64 {
    let mut inverse_row_sums = [0.0; 6];
    let mut unit_column = 0;
    while unit_column < 6 {
        let mut unit = [0.0; 6];
        unit[unit_column] = 1.0;
        let inverse_column = solve_cholesky6(lower, unit);
        let mut inverse_row = 0;
        while inverse_row < 6 {
            inverse_row_sums[inverse_row] += inverse_column[inverse_row].abs();
            inverse_row += 1;
        }
        unit_column += 1;
    }
    let mut matrix_norm = 0.0_f64;
    let mut row = 0;
    while row < 6 {
        let mut row_sum = 0.0;
        let mut column = 0;
        while column < 6 {
            row_sum += normalized[row * 6 + column].abs();
            column += 1;
        }
        matrix_norm = matrix_norm.max(row_sum);
        row += 1;
    }
    let inverse_norm = inverse_row_sums.into_iter().fold(0.0_f64, f64::max);
    matrix_norm * inverse_norm
}

fn solve_cholesky6(lower: &[f64; 36], mut right_hand_side: [f64; 6]) -> [f64; 6] {
    let mut row = 0;
    while row < 6 {
        let mut value = right_hand_side[row];
        let mut column = 0;
        while column < row {
            value = (-lower[row * 6 + column]).mul_add(right_hand_side[column], value);
            column += 1;
        }
        right_hand_side[row] = value / lower[row * 6 + row];
        row += 1;
    }
    row = 6;
    while row > 0 {
        row -= 1;
        let mut value = right_hand_side[row];
        let mut column = row + 1;
        while column < 6 {
            value = (-lower[column * 6 + row]).mul_add(right_hand_side[column], value);
            column += 1;
        }
        right_hand_side[row] = value / lower[row * 6 + row];
    }
    right_hand_side
}

fn validate_slice(values: &[f64], field: &'static str) -> Result<(), ArticulatedError> {
    for (index, value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(ArticulatedError::NonFinite { field, index });
        }
    }
    Ok(())
}

fn validate_vec3(value: Vec3, field: &'static str) -> Result<(), ArticulatedError> {
    validate_slice(&[value.x, value.y, value.z], field)
}

fn validate_twist(value: Twist, field: &'static str) -> Result<(), ArticulatedError> {
    validate_slice(&value.to_array(), field)
}

fn validate_wrench(value: Wrench, field: &'static str) -> Result<(), ArticulatedError> {
    validate_slice(&value.to_array(), field)
}

fn validate_mat6(value: &Mat6, field: &'static str) -> Result<(), ArticulatedError> {
    validate_slice(&value.m, field)
}

fn ensure_finite_scalar(
    value: f64,
    field: &'static str,
    index: usize,
) -> Result<(), ArticulatedError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ArticulatedError::NonFinite { field, index })
    }
}

fn require_len(field: &'static str, got: usize, expected: usize) -> Result<(), ArticulatedError> {
    if got == expected {
        Ok(())
    } else {
        Err(ArticulatedError::LengthMismatch {
            field,
            expected,
            got,
        })
    }
}

fn normalized_axis(axis: Vec3) -> Result<Vec3, ArticulatedError> {
    validate_vec3(axis, "joint.axis")?;
    let scale = axis.x.abs().max(axis.y.abs()).max(axis.z.abs());
    if scale == 0.0 {
        return Err(ArticulatedError::InvalidJointAxis);
    }
    let scaled = Vec3::new(axis.x / scale, axis.y / scale, axis.z / scale);
    let norm = scaled.dot(scaled).sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return Err(ArticulatedError::InvalidJointAxis);
    }
    Ok(scaled.scale(norm.recip()))
}

fn symmetrize(matrix: Mat3) -> Mat3 {
    let xy = 0.5 * (matrix.m[1] + matrix.m[3]);
    let xz = 0.5 * (matrix.m[2] + matrix.m[6]);
    let yz = 0.5 * (matrix.m[5] + matrix.m[7]);
    Mat3 {
        m: [
            matrix.m[0],
            xy,
            xz,
            xy,
            matrix.m[4],
            yz,
            xz,
            yz,
            matrix.m[8],
        ],
    }
}

fn divide_mat3(mut matrix: Mat3, divisor: f64) -> Mat3 {
    for entry in &mut matrix.m {
        *entry /= divisor;
    }
    matrix
}

fn determinant3(matrix: Mat3) -> f64 {
    matrix.m[0] * (matrix.m[4] * matrix.m[8] - matrix.m[5] * matrix.m[7])
        - matrix.m[1] * (matrix.m[3] * matrix.m[8] - matrix.m[5] * matrix.m[6])
        + matrix.m[2] * (matrix.m[3] * matrix.m[7] - matrix.m[4] * matrix.m[6])
}

fn principal_minors(matrix: Mat3) -> [f64; 7] {
    [
        matrix.m[0],
        matrix.m[4],
        matrix.m[8],
        matrix.m[0] * matrix.m[4] - matrix.m[1] * matrix.m[3],
        matrix.m[0] * matrix.m[8] - matrix.m[2] * matrix.m[6],
        matrix.m[4] * matrix.m[8] - matrix.m[5] * matrix.m[7],
        determinant3(matrix),
    ]
}

fn minimum_principal_minor(matrix: Mat3) -> f64 {
    principal_minors(matrix)
        .into_iter()
        .fold(f64::INFINITY, f64::min)
}

fn minimum_leading_principal_minor(matrix: Mat3) -> f64 {
    matrix.m[0]
        .min(matrix.m[0] * matrix.m[4] - matrix.m[1] * matrix.m[3])
        .min(determinant3(matrix))
}

fn is_positive_semidefinite(matrix: Mat3, tolerance: f64) -> bool {
    principal_minors(matrix)
        .into_iter()
        .all(|minor| minor >= -tolerance)
}

fn is_positive_definite(matrix: Mat3, tolerance: f64) -> bool {
    matrix.m[0] > tolerance
        && matrix.m[0] * matrix.m[4] - matrix.m[1] * matrix.m[3] > tolerance
        && determinant3(matrix) > tolerance
}

fn hat(vector: Vec3) -> Mat3 {
    Mat3 {
        m: [
            0.0, -vector.z, vector.y, vector.z, 0.0, -vector.x, -vector.y, vector.x, 0.0,
        ],
    }
}

fn cross_force(motion: Twist, force: Wrench) -> Wrench {
    Wrench::new(
        motion.angular.cross(force.torque) + motion.linear.cross(force.force),
        motion.angular.cross(force.force),
    )
}

fn wrench_plus(lhs: Wrench, rhs: Wrench) -> Wrench {
    Wrench::new(lhs.torque + rhs.torque, lhs.force + rhs.force)
}

fn wrench_minus(lhs: Wrench, rhs: Wrench) -> Wrench {
    Wrench::new(lhs.torque - rhs.torque, lhs.force - rhs.force)
}

fn wrench_from_array(value: [f64; 6]) -> Wrench {
    Wrench::new(
        Vec3::new(value[0], value[1], value[2]),
        Vec3::new(value[3], value[4], value[5]),
    )
}

fn twist_from_array(value: [f64; 6]) -> Twist {
    Twist::new(
        Vec3::new(value[0], value[1], value[2]),
        Vec3::new(value[3], value[4], value[5]),
    )
}

fn mat6_apply_twist(matrix: &Mat6, value: Twist) -> Wrench {
    wrench_from_array(mat6_apply_array(matrix, value.to_array()))
}

fn mat6_apply_array(matrix: &Mat6, value: [f64; 6]) -> [f64; 6] {
    let mut output = [0.0; 6];
    let mut row = 0;
    while row < 6 {
        let mut column = 0;
        while column < 6 {
            output[row] += matrix.m[row * 6 + column] * value[column];
            column += 1;
        }
        row += 1;
    }
    output
}

fn mat6_compose(lhs: &Mat6, rhs: &Mat6) -> Mat6 {
    (*lhs).compose(*rhs)
}

fn mat6_add(lhs: &Mat6, rhs: &Mat6) -> Mat6 {
    let mut output = *lhs;
    let mut index = 0;
    while index < 36 {
        output.m[index] += rhs.m[index];
        index += 1;
    }
    output
}

fn mat6_sub_outer_scaled(matrix: &Mat6, lhs: [f64; 6], rhs: [f64; 6], scale: f64) -> Mat6 {
    let mut output = *matrix;
    let mut row = 0;
    while row < 6 {
        let mut column = 0;
        while column < 6 {
            output.m[row * 6 + column] -= scale * lhs[row] * rhs[column];
            column += 1;
        }
        row += 1;
    }
    output
}

fn dot6(lhs: [f64; 6], rhs: [f64; 6]) -> f64 {
    let mut result = 0.0;
    let mut index = 0;
    while index < 6 {
        result += lhs[index] * rhs[index];
        index += 1;
    }
    result
}

fn dot_abs6(lhs: [f64; 6], rhs: [f64; 6]) -> f64 {
    let mut result = 0.0;
    let mut index = 0;
    while index < 6 {
        result += (lhs[index] * rhs[index]).abs();
        index += 1;
    }
    result
}

fn scale6(mut value: [f64; 6], scale: f64) -> [f64; 6] {
    for coordinate in &mut value {
        *coordinate *= scale;
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 2.0e-10;

    fn diagonal(xx: f64, yy: f64, zz: f64) -> Mat3 {
        Mat3 {
            m: [xx, 0.0, 0.0, 0.0, yy, 0.0, 0.0, 0.0, zz],
        }
    }

    fn test_inertia(center: Vec3) -> SpatialInertia {
        SpatialInertia::new(1.0, center, diagonal(0.02, 0.26, 0.26)).unwrap()
    }

    fn centered_inertia(mass: f64, xx: f64, yy: f64, zz: f64) -> SpatialInertia {
        SpatialInertia::new(mass, Vec3::new(0.0, 0.0, 0.0), diagonal(xx, yy, zz)).unwrap()
    }

    fn zero_wrenches(count: usize) -> Vec<Wrench> {
        vec![Wrench::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)); count]
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:.16e}, got {actual:.16e}"
        );
    }

    fn assert_twist_close(actual: Twist, expected: Twist, tolerance: f64) {
        for (actual, expected) in actual.to_array().into_iter().zip(expected.to_array()) {
            assert_close(actual, expected, tolerance);
        }
    }

    fn floating_pendulum_model() -> ArticulatedModel {
        let root = Link::new(
            "floating_root",
            None,
            Se3::identity(),
            JointModel::FIXED,
            centered_inertia(2.0, 0.7, 0.8, 0.9),
        );
        let offset = Se3::from_parts(fs_ga::So3::identity(), Vec3::new(0.8, 0.0, 0.0)).unwrap();
        let pendulum = Link::new(
            "pendulum",
            Some(0),
            offset,
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            SpatialInertia::new(1.5, Vec3::new(0.4, 0.0, 0.0), diagonal(0.2, 0.3, 0.4)).unwrap(),
        );
        ArticulatedModel::new(vec![root, pendulum]).unwrap()
    }

    fn two_link_model() -> ArticulatedModel {
        let revolute = JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap();
        let root = Link::new(
            "shoulder",
            None,
            Se3::identity(),
            revolute,
            test_inertia(Vec3::new(0.5, 0.0, 0.0)),
        );
        let offset = Se3::from_parts(fs_ga::So3::identity(), Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let child = Link::new(
            "elbow",
            Some(0),
            offset,
            revolute,
            test_inertia(Vec3::new(0.4, 0.0, 0.0)),
        );
        ArticulatedModel::new(vec![root, child]).unwrap()
    }

    #[test]
    fn spatial_inertia_matrix_matches_direct_momentum_and_energy() {
        let inertia = SpatialInertia::new(
            3.0,
            Vec3::new(0.2, -0.1, 0.3),
            Mat3 {
                m: [0.9, 0.04, -0.02, 0.04, 1.1, 0.03, -0.02, 0.03, 1.2],
            },
        )
        .unwrap();
        let velocity = Twist::new(Vec3::new(0.7, -0.4, 0.2), Vec3::new(1.0, -0.5, 0.3));
        let direct = inertia.momentum(velocity).unwrap();
        let matrix = mat6_apply_twist(&inertia.matrix().unwrap(), velocity);
        for (actual, expected) in direct.to_array().into_iter().zip(matrix.to_array()) {
            assert_close(actual, expected, EPSILON);
        }
        assert_close(
            inertia.kinetic_energy(velocity).unwrap(),
            0.5 * matrix.pairing(velocity),
            EPSILON,
        );
    }

    #[test]
    fn physical_inertia_validation_enforces_symmetry_and_triangle_inequalities() {
        let nonsymmetric = Mat3 {
            m: [1.0, 0.2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        };
        assert!(matches!(
            SpatialInertia::new(1.0, Vec3::new(0.0, 0.0, 0.0), nonsymmetric),
            Err(ArticulatedError::NonSymmetricInertia { .. })
        ));
        assert!(matches!(
            SpatialInertia::new(1.0, Vec3::new(0.0, 0.0, 0.0), diagonal(3.0, 1.0, 1.0)),
            Err(ArticulatedError::NonPhysicalInertia { .. })
        ));
        assert_eq!(
            SpatialInertia::new(1.0, Vec3::new(0.0, 0.0, 0.0), diagonal(-1.0, -1.0, 1.0),),
            Err(ArticulatedError::NonPhysicalInertia { measure: -1.0 })
        );
        let subnormal = f64::from_bits(1);
        assert!(JointModel::revolute(Vec3::new(subnormal, 0.0, 0.0), None).is_ok());
    }

    #[test]
    fn spatial_boundaries_refuse_nonfinite_or_overflowing_derived_values() {
        let inertia = SpatialInertia::new(
            f64::MAX,
            Vec3::new(f64::MAX, 0.0, 0.0),
            diagonal(1.0, 1.0, 1.0),
        )
        .unwrap();
        assert!(matches!(
            inertia.matrix(),
            Err(ArticulatedError::NonFinite {
                field: "spatial_inertia.matrix",
                ..
            })
        ));
        assert!(matches!(
            inertia.momentum(Twist::new(
                Vec3::new(f64::NAN, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 0.0),
            )),
            Err(ArticulatedError::NonFinite {
                field: "spatial_inertia.velocity",
                ..
            })
        ));

        let prismatic = ArticulatedModel::new(vec![Link::new(
            "overflowing",
            None,
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            test_inertia(Vec3::new(0.0, 0.0, 0.0)),
        )])
        .unwrap();
        let base = BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(f64::MAX, 0.0, 0.0)),
            Twist::zero(),
        );
        assert!(matches!(
            forward_kinematics(&prismatic, base, &[0.0], &[f64::MAX]),
            Err(ArticulatedError::NonFinite {
                field: "forward_kinematics.body_twist",
                ..
            })
        ));
    }

    #[test]
    fn forward_kinematics_composes_two_revolute_links_in_the_declared_order() {
        let model = two_link_model();
        let kinematics = forward_kinematics(
            &model,
            BaseState::stationary(Se3::identity()),
            &[core::f64::consts::FRAC_PI_2, 0.0],
            &[0.0, 0.0],
        )
        .unwrap();
        let elbow_origin = kinematics.world_from_link[1]
            .transform_point(Vec3::new(0.0, 0.0, 0.0))
            .unwrap();
        assert_close(elbow_origin.x, 0.0, EPSILON);
        assert_close(elbow_origin.y, 1.0, EPSILON);
        assert_close(elbow_origin.z, 0.0, EPSILON);
    }

    #[test]
    fn joint_motion_rotates_about_its_declared_nonzero_parent_origin() {
        let offset = Se3::from_parts(fs_ga::So3::identity(), Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let model = ArticulatedModel::new(vec![Link::new(
            "offset_joint",
            None,
            offset,
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            test_inertia(Vec3::new(0.0, 0.0, 0.0)),
        )])
        .unwrap();
        let kinematics = forward_kinematics(
            &model,
            BaseState::stationary(Se3::identity()),
            &[core::f64::consts::FRAC_PI_2],
            &[0.0],
        )
        .unwrap();
        let origin = kinematics.world_from_link[0]
            .transform_point(Vec3::new(0.0, 0.0, 0.0))
            .unwrap();
        let local_x = kinematics.world_from_link[0]
            .transform_point(Vec3::new(1.0, 0.0, 0.0))
            .unwrap();
        assert_close(origin.x, 1.0, EPSILON);
        assert_close(origin.y, 0.0, EPSILON);
        assert_close(local_x.x, 1.0, EPSILON);
        assert_close(local_x.y, 1.0, EPSILON);
    }

    #[test]
    fn single_pendulum_inverse_dynamics_matches_closed_form_gravity_torque() {
        let model = ArticulatedModel::new(vec![Link::new(
            "pendulum",
            None,
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            test_inertia(Vec3::new(0.5, 0.0, 0.0)),
        )])
        .unwrap();
        let result = inverse_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &[0.0],
            &[0.0],
            &[0.0],
            Vec3::new(0.0, -9.81, 0.0),
            &zero_wrenches(1),
        )
        .unwrap();
        assert_close(result.generalized_force[0], 4.905, EPSILON);
    }

    #[test]
    fn inverse_dynamics_includes_prescribed_non_gravity_base_acceleration() {
        let model = ArticulatedModel::new(vec![Link::new(
            "pendulum",
            None,
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            test_inertia(Vec3::new(0.5, 0.0, 0.0)),
        )])
        .unwrap();
        let base = BaseState::prescribed(
            Se3::identity(),
            Twist::zero(),
            Twist::new(Vec3::new(0.0, 0.0, 2.0), Vec3::new(0.0, 0.0, 0.0)),
        );
        let result = inverse_dynamics(
            &model,
            base,
            &[0.0],
            &[0.0],
            &[0.0],
            Vec3::new(0.0, 0.0, 0.0),
            &zero_wrenches(1),
        )
        .unwrap();
        assert_close(result.generalized_force[0], 1.02, EPSILON);
    }

    #[test]
    fn inverse_dynamics_mixed_branch_matches_force_and_external_torque_oracles() {
        let root = Link::new(
            "base_link",
            None,
            Se3::identity(),
            JointModel::FIXED,
            test_inertia(Vec3::new(0.0, 0.0, 0.0)),
        );
        let slider = Link::new(
            "vertical_slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            SpatialInertia::new(2.0, Vec3::new(0.0, 0.0, 0.0), diagonal(0.1, 0.1, 0.1)).unwrap(),
        );
        let rotor = Link::new(
            "torque_rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            test_inertia(Vec3::new(0.0, 0.0, 0.0)),
        );
        let model = ArticulatedModel::new(vec![root, slider, rotor]).unwrap();
        let mut external = zero_wrenches(3);
        external[2] = Wrench::new(Vec3::new(0.0, 0.0, 3.0), Vec3::new(0.0, 0.0, 0.0));
        let result = inverse_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[0.0, 0.0],
            Vec3::new(0.0, -9.81, 0.0),
            &external,
        )
        .unwrap();
        assert_close(result.generalized_force[0], 19.62, EPSILON);
        assert_close(result.generalized_force[1], -3.0, EPSILON);
    }

    #[test]
    fn rotated_base_changes_world_gravity_in_the_declared_frame() {
        let model = ArticulatedModel::new(vec![Link::new(
            "pendulum",
            None,
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            test_inertia(Vec3::new(0.5, 0.0, 0.0)),
        )])
        .unwrap();
        let world_from_base = Se3::exp(Twist::new(
            Vec3::new(0.0, 0.0, core::f64::consts::FRAC_PI_2),
            Vec3::new(0.0, 0.0, 0.0),
        ))
        .unwrap();
        let result = inverse_dynamics(
            &model,
            BaseState::stationary(world_from_base),
            &[0.0],
            &[0.0],
            &[0.0],
            Vec3::new(0.0, -9.81, 0.0),
            &zero_wrenches(1),
        )
        .unwrap();
        assert_close(result.generalized_force[0], 0.0, EPSILON);
    }

    #[test]
    fn articulated_forward_dynamics_matches_single_pendulum_oracle() {
        let model = ArticulatedModel::new(vec![Link::new(
            "pendulum",
            None,
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            test_inertia(Vec3::new(0.5, 0.0, 0.0)),
        )])
        .unwrap();
        let result = forward_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &[0.0],
            &[0.0],
            &[0.0],
            Vec3::new(0.0, -9.81, 0.0),
            &zero_wrenches(1),
        )
        .unwrap();
        assert_close(result.generalized_acceleration[0], -4.905 / 0.51, 2.0e-9);
    }

    #[test]
    fn tiny_but_valid_spatial_inertia_does_not_trip_a_unit_dependent_pivot_floor() {
        let inertia = SpatialInertia::new(
            1.0e-15,
            Vec3::new(0.0, 0.0, 0.0),
            diagonal(1.0e-15, 1.0e-15, 1.0e-15),
        )
        .unwrap();
        let model = ArticulatedModel::new(vec![Link::new(
            "tiny_slider",
            None,
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia,
        )])
        .unwrap();
        let result = forward_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &[0.0],
            &[0.0],
            &[0.0],
            Vec3::new(0.0, 0.0, 0.0),
            &zero_wrenches(1),
        )
        .unwrap();
        assert_eq!(result.generalized_acceleration, vec![0.0]);
    }

    #[test]
    fn inverse_and_forward_dynamics_round_trip_on_a_coupled_tree() {
        let model = two_link_model();
        let q = [0.3, -0.7];
        let qd = [0.4, -0.2];
        let qdd = [1.1, -0.6];
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let external = zero_wrenches(model.link_count());
        let inverse = inverse_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &q,
            &qd,
            &qdd,
            gravity,
            &external,
        )
        .unwrap();
        let forward = forward_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &q,
            &qd,
            &inverse.generalized_force,
            gravity,
            &external,
        )
        .unwrap();
        for (actual, expected) in forward.generalized_acceleration.iter().zip(qdd.iter()) {
            assert_close(*actual, *expected, 2.0e-9);
        }
    }

    #[test]
    fn free_rigid_body_accelerates_with_uniform_gravity() {
        let inertia = centered_inertia(2.0, 0.4, 0.5, 0.6);
        let model = ArticulatedModel::new(vec![Link::new(
            "free_body",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia,
        )])
        .unwrap();
        let gravity = Vec3::new(0.4, -9.7, 1.2);
        let result = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::stationary(Se3::identity()),
            &[],
            &[],
            &[],
            gravity,
            &zero_wrenches(1),
        )
        .unwrap();
        let expected = Twist::new(Vec3::new(0.0, 0.0, 0.0), gravity);
        assert_twist_close(result.base_spatial_acceleration_body, expected, EPSILON);
        assert_twist_close(result.body_spatial_acceleration[0], expected, EPSILON);
        assert!(result.generalized_acceleration.is_empty());
    }

    #[test]
    fn free_rigid_body_no_force_obeys_instantaneous_conservation_identity() {
        let inertia = centered_inertia(2.0, 0.4, 0.5, 0.6);
        let model = ArticulatedModel::new(vec![Link::new(
            "free_body",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia,
        )])
        .unwrap();
        let velocity = Twist::new(Vec3::new(0.7, -0.4, 0.2), Vec3::new(1.0, -0.5, 0.3));
        let result = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::new(Se3::identity(), velocity),
            &[],
            &[],
            &[],
            Vec3::new(0.0, 0.0, 0.0),
            &zero_wrenches(1),
        )
        .unwrap();
        let momentum = inertia.momentum(velocity).unwrap();
        let balance = wrench_plus(
            inertia
                .momentum(result.base_spatial_acceleration_body)
                .unwrap(),
            cross_force(velocity, momentum),
        );
        for residual in balance.to_array() {
            assert_close(residual, 0.0, EPSILON);
        }
        assert_close(
            momentum.pairing(result.base_spatial_acceleration_body),
            0.0,
            EPSILON,
        );
    }

    #[test]
    fn free_rigid_body_external_wrench_matches_closed_form_acceleration() {
        let model = ArticulatedModel::new(vec![Link::new(
            "free_body",
            None,
            Se3::identity(),
            JointModel::FIXED,
            centered_inertia(2.0, 0.4, 0.5, 0.6),
        )])
        .unwrap();
        let external = [Wrench::new(
            Vec3::new(0.8, -1.0, 1.8),
            Vec3::new(4.0, -3.0, 2.0),
        )];
        let result = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::stationary(Se3::identity()),
            &[],
            &[],
            &[],
            Vec3::new(0.0, 0.0, 0.0),
            &external,
        )
        .unwrap();
        assert_twist_close(
            result.base_spatial_acceleration_body,
            Twist::new(Vec3::new(2.0, -2.0, 3.0), Vec3::new(2.0, -1.5, 1.0)),
            EPSILON,
        );
    }

    #[test]
    fn free_base_and_one_link_pendulum_have_the_expected_reaction_coupling() {
        let base_mass = 2.0;
        let base_inertia_z = 0.8;
        let pendulum_mass = 1.0;
        let pendulum_inertia_z = 0.4;
        let center_distance = 0.5;
        let root = Link::new(
            "base_rotor",
            None,
            Se3::identity(),
            JointModel::FIXED,
            centered_inertia(base_mass, 0.6, 0.7, base_inertia_z),
        );
        let child = Link::new(
            "pendulum",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            SpatialInertia::new(
                pendulum_mass,
                Vec3::new(center_distance, 0.0, 0.0),
                diagonal(0.2, 0.3, pendulum_inertia_z),
            )
            .unwrap(),
        );
        let model = ArticulatedModel::new(vec![root, child]).unwrap();
        let torque = 1.2;
        let result = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::stationary(Se3::identity()),
            &[0.0],
            &[0.0],
            &[torque],
            Vec3::new(0.0, 0.0, 0.0),
            &zero_wrenches(2),
        )
        .unwrap();
        let absolute_pendulum_inertia = pendulum_inertia_z
            + pendulum_mass * center_distance * center_distance * base_mass
                / (base_mass + pendulum_mass);
        let absolute_pendulum_acceleration = torque / absolute_pendulum_inertia;
        assert_close(
            result.base_spatial_acceleration_body.angular.z,
            -torque / base_inertia_z,
            EPSILON,
        );
        assert_close(
            result.base_spatial_acceleration_body.linear.y,
            -pendulum_mass * center_distance * absolute_pendulum_acceleration
                / (base_mass + pendulum_mass),
            EPSILON,
        );
        assert_close(
            result.generalized_acceleration[0],
            absolute_pendulum_acceleration + torque / base_inertia_z,
            EPSILON,
        );
    }

    #[test]
    fn free_flight_is_equivariant_under_a_common_world_rotation() {
        let inertia =
            SpatialInertia::new(2.0, Vec3::new(0.2, -0.1, 0.3), diagonal(0.4, 0.5, 0.6)).unwrap();
        let model = ArticulatedModel::new(vec![Link::new(
            "offset_free_body",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia,
        )])
        .unwrap();
        let velocity = Twist::new(Vec3::new(0.3, -0.2, 0.4), Vec3::new(0.5, 0.1, -0.3));
        let gravity = Vec3::new(0.2, -9.8, 0.7);
        let external = [Wrench::new(
            Vec3::new(0.4, -0.3, 0.2),
            Vec3::new(1.0, -0.5, 0.7),
        )];
        let reference = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::new(Se3::identity(), velocity),
            &[],
            &[],
            &[],
            gravity,
            &external,
        )
        .unwrap();
        let rotated_pose = Se3::exp(Twist::new(
            Vec3::new(0.3, -0.4, 0.7),
            Vec3::new(1.2, -0.6, 0.8),
        ))
        .unwrap();
        let rotated_gravity = rotated_pose.transform_vector(gravity).unwrap();
        let rotated = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::new(rotated_pose, velocity),
            &[],
            &[],
            &[],
            rotated_gravity,
            &external,
        )
        .unwrap();
        assert_twist_close(
            rotated.base_spatial_acceleration_body,
            reference.base_spatial_acceleration_body,
            2.0e-9,
        );
        assert_twist_close(
            rotated.body_spatial_acceleration[0],
            reference.body_spatial_acceleration[0],
            2.0e-9,
        );
    }

    #[test]
    fn free_floating_aba_round_trips_through_prescribed_base_rnea() {
        let model = floating_pendulum_model();
        let base = FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.4, 0.2, -0.1)),
        );
        let q = [0.3];
        let qd = [0.2];
        let applied = [0.7];
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let mut external = zero_wrenches(model.link_count());
        external[1] = Wrench::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.4, 0.2, -0.1));
        let forward =
            free_floating_forward_dynamics(&model, base, &q, &qd, &applied, gravity, &external)
                .unwrap();
        let inverse = inverse_dynamics(
            &model,
            BaseState::prescribed(
                base.world_from_base,
                base.twist_body,
                forward.base_spatial_acceleration_body,
            ),
            &q,
            &qd,
            &forward.generalized_acceleration,
            gravity,
            &external,
        )
        .unwrap();
        assert_close(inverse.generalized_force[0], applied[0], 3.0e-9);
        for residual in inverse.body_wrench[0].to_array() {
            assert_close(residual, 0.0, 3.0e-9);
        }
    }

    #[test]
    fn free_floating_boundary_refuses_nonfinite_and_redundant_inputs() {
        let model = ArticulatedModel::new(vec![Link::new(
            "finite_root",
            None,
            Se3::identity(),
            JointModel::FIXED,
            centered_inertia(1.0, 0.4, 0.5, 0.6),
        )])
        .unwrap();
        let invalid_base = FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(f64::NAN, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)),
        );
        assert!(matches!(
            free_floating_forward_dynamics(
                &model,
                invalid_base,
                &[],
                &[],
                &[],
                Vec3::new(0.0, 0.0, 0.0),
                &zero_wrenches(1),
            ),
            Err(ArticulatedError::NonFinite {
                field: "base.twist_body",
                index: 0,
            })
        ));

        let redundant = ArticulatedModel::new(vec![Link::new(
            "redundant_root",
            None,
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            centered_inertia(1.0, 0.4, 0.5, 0.6),
        )])
        .unwrap();
        assert_eq!(
            redundant.free_floating_complexity(),
            Err(ArticulatedError::FloatingBaseRootJointNotFixed)
        );
        assert_eq!(
            free_floating_forward_dynamics(
                &redundant,
                FreeFloatingBaseState::stationary(Se3::identity()),
                &[0.0],
                &[0.0],
                &[0.0],
                Vec3::new(0.0, 0.0, 0.0),
                &zero_wrenches(1),
            ),
            Err(ArticulatedError::FloatingBaseRootJointNotFixed)
        );
    }

    #[test]
    fn floating_base_root_solver_refuses_singular_and_ill_conditioned_systems() {
        let mut nonfinite = Mat6::identity();
        nonfinite.m[10] = f64::INFINITY;
        assert!(matches!(
            solve_floating_base_system(
                &nonfinite,
                Wrench::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)),
            ),
            Err(ArticulatedError::NonFinite {
                field: "free_floating_forward_dynamics.base_articulated_inertia",
                index: 10,
            })
        ));

        let mut nonsymmetric = Mat6::identity();
        nonsymmetric.m[1] = 1.0e-6;
        assert!(matches!(
            solve_floating_base_system(
                &nonsymmetric,
                Wrench::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)),
            ),
            Err(ArticulatedError::NonSymmetricFloatingBaseInertia { .. })
        ));

        let ill_scaled = ArticulatedModel::new(vec![Link::new(
            "ill_scaled_root",
            None,
            Se3::identity(),
            JointModel::FIXED,
            SpatialInertia::new(1.0, Vec3::new(1.0e8, 0.0, 0.0), diagonal(1.0, 1.0, 1.0)).unwrap(),
        )])
        .unwrap();
        assert!(matches!(
            free_floating_forward_dynamics(
                &ill_scaled,
                FreeFloatingBaseState::stationary(Se3::identity()),
                &[],
                &[],
                &[],
                Vec3::new(0.0, 0.0, 0.0),
                &zero_wrenches(1),
            ),
            Err(ArticulatedError::SingularFloatingBaseInertia { .. })
        ));

        let mut ill_conditioned = Mat6::identity();
        ill_conditioned.m[35] = 5.0e-13;
        assert!(matches!(
            solve_floating_base_system(
                &ill_conditioned,
                Wrench::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)),
            ),
            Err(ArticulatedError::IllConditionedFloatingBaseInertia { .. })
        ));
    }

    #[test]
    fn floating_base_root_solver_handles_coupling_and_origin_acceleration_is_explicit() {
        let mut coupled = Mat6::identity();
        coupled.m[1] = 0.2;
        coupled.m[6] = 0.2;
        coupled.m[22] = -0.3;
        coupled.m[27] = -0.3;
        let expected = Twist::new(Vec3::new(1.0, -2.0, 0.5), Vec3::new(0.25, -0.75, 1.5));
        let right_hand_side = mat6_apply_twist(&coupled, expected);
        let actual = solve_floating_base_system(&coupled, right_hand_side).unwrap();
        assert_twist_close(actual, expected, EPSILON);

        let body_twist = Twist::new(Vec3::new(0.0, 0.0, 2.0), Vec3::new(3.0, 4.0, 0.0));
        let spatial_acceleration = Twist::new(Vec3::new(0.1, 0.2, 0.3), Vec3::new(1.0, 2.0, 3.0));
        let origin = origin_linear_acceleration_body(spatial_acceleration, body_twist).unwrap();
        assert_close(origin.x, -7.0, EPSILON);
        assert_close(origin.y, 8.0, EPSILON);
        assert_close(origin.z, 3.0, EPSILON);
    }

    #[test]
    fn g1_catalog_runs_a_linear_storage_free_fall_solve_without_internal_actuation() {
        let catalog = crate::robot_models::unitree_g1_lower_body_waist_15dof().unwrap();
        let model = catalog.model();
        let complexity = model.free_floating_complexity().unwrap();
        assert_eq!(complexity.tree.links, 16);
        assert_eq!(complexity.tree.degrees_of_freedom, 15);
        assert_eq!(complexity.tree.dense_generalized_matrix_entries, 0);
        assert_eq!(complexity.tree.spatial_matrix_entries_per_link, 36);
        assert_eq!(complexity.base_degrees_of_freedom, 6);
        assert_eq!(complexity.fixed_root_solve_matrix_entries, 36);

        let zero = vec![0.0; model.dof_count()];
        let gravity = Vec3::new(0.0, 0.0, -9.81);
        let result = free_floating_forward_dynamics(
            model,
            FreeFloatingBaseState::stationary(Se3::identity()),
            &zero,
            &zero,
            &zero,
            gravity,
            &zero_wrenches(model.link_count()),
        )
        .unwrap();
        assert_twist_close(
            result.base_spatial_acceleration_body,
            Twist::new(Vec3::new(0.0, 0.0, 0.0), gravity),
            3.0e-9,
        );
        for acceleration in &result.generalized_acceleration {
            assert_close(*acceleration, 0.0, 3.0e-9);
        }

        // Uniform gravity cannot excite internal coordinates in unconstrained
        // free fall. Each link sees the same world acceleration expressed in
        // its own body frame; this catches a gravity sign or transform error
        // that a mere finite-value smoke assertion would miss.
        let kinematics =
            forward_kinematics(model, BaseState::stationary(Se3::identity()), &zero, &zero)
                .unwrap();
        for (link_index, acceleration) in result.body_spatial_acceleration.iter().enumerate() {
            let gravity_body = kinematics.world_from_link[link_index]
                .rotation()
                .inverse()
                .rotate(gravity)
                .unwrap();
            assert_twist_close(
                *acceleration,
                Twist::new(Vec3::new(0.0, 0.0, 0.0), gravity_body),
                3.0e-9,
            );
        }
    }

    #[test]
    fn model_reports_linear_articulated_working_set_and_checks_limits() {
        let limits = JointLimits::new(-1.0, 1.0, 2.0, 20.0).unwrap();
        let model = ArticulatedModel::new(vec![Link::new(
            "limited",
            None,
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 2.0), Some(limits)).unwrap(),
            test_inertia(Vec3::new(0.5, 0.0, 0.0)),
        )])
        .unwrap();
        let complexity = model.complexity();
        assert_eq!(complexity.links, 1);
        assert_eq!(complexity.degrees_of_freedom, 1);
        assert_eq!(complexity.dense_generalized_matrix_entries, 0);
        assert_eq!(complexity.spatial_matrix_entries_per_link, 36);
        assert!(matches!(
            forward_kinematics(
                &model,
                BaseState::stationary(Se3::identity()),
                &[1.1],
                &[0.0]
            ),
            Err(ArticulatedError::JointLimitViolation { link: 0, .. })
        ));
        assert!(matches!(
            forward_kinematics(
                &model,
                BaseState::stationary(Se3::identity()),
                &[0.0],
                &[2.1]
            ),
            Err(ArticulatedError::JointVelocityLimitViolation { link: 0, .. })
        ));
        assert!(matches!(
            forward_dynamics(
                &model,
                BaseState::stationary(Se3::identity()),
                &[0.0],
                &[0.0],
                &[20.1],
                Vec3::new(0.0, -9.81, 0.0),
                &zero_wrenches(1),
            ),
            Err(ArticulatedError::JointEffortLimitViolation { link: 0, .. })
        ));
    }
}
