//! Native generalized inertial terms for coupling an articulated tree.
//!
//! Coordinates are `[omega_base, v_base, joint_rate]`, with the base twist in
//! base axes. The equation is `M(q) udot + bias(q,u) = applied_force`; its base
//! linear acceleration is the derivative of BODY velocity. The world origin
//! acceleration is `R (udot_linear + omega cross v)`. Gravity is an explicit
//! projected load, never included in `bias`.
//!
//! Composite rigid-body projection supplies the full mass matrix. The original
//! Newton-Euler recursion supplies bias and each link's zero-generalized-
//! acceleration motion; no second dynamics or time integration engine is used.

use super::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, Kinematics,
    forward_kinematics, inverse_dynamics, mat6_add, mat6_apply_twist, require_len,
    twist_from_array, validate_mat6, validate_slice, validate_twist, validate_vec3,
    validate_wrench,
};
use core::fmt;
use fs_ga::{GaError, Mat6, Se3, Twist, Vec3, Wrench};

/// Link cap bounding the original indivisible FK and Newton-Euler passes.
pub const MAX_GENERALIZED_LINKS: usize = 4096;
/// Scalar-joint cap for the explicit generalized matrix and link Jacobians.
pub const MAX_GENERALIZED_JOINTS: usize = 64;

/// Atomic refusal of a generalized projection.
#[derive(Clone, Debug, PartialEq)]
pub enum GeneralizedError {
    /// Explicit projection exceeds its documented working-set bound.
    SizeLimit {
        /// Supplied links, including fixed links.
        links: usize,
        /// Supplied scalar joints.
        joints: usize,
    },
    /// Cancellation was observed before publishing the result.
    Cancelled,
    /// Original model, geometry, state or finite-output admission refused.
    Dynamics(ArticulatedError),
}

impl fmt::Display for GeneralizedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SizeLimit { links, joints } => write!(
                f,
                "generalized projection has {links} links/{joints} joints; maxima are {MAX_GENERALIZED_LINKS}/{MAX_GENERALIZED_JOINTS}"
            ),
            Self::Cancelled => f.write_str("generalized projection cancelled"),
            Self::Dynamics(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for GeneralizedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(error) => Some(error),
            _ => None,
        }
    }
}
impl From<ArticulatedError> for GeneralizedError {
    fn from(value: ArticulatedError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<GaError> for GeneralizedError {
    fn from(value: GaError) -> Self {
        Self::Dynamics(ArticulatedError::Geometry(value))
    }
}

/// Admitted link motion at one configuration and generalized velocity.
///
/// For link `i`, `V_i = J_i u` and `A_i = J_i udot + c_i`, all in that link's
/// body axes. Columns include six base coordinates followed by compact joints.
/// This same map projects structural, fluid and applied-load operators.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneralizedKinematics {
    kinematics: Kinematics,
    body_jacobian: Vec<Vec<Twist>>,
    body_bias_acceleration: Vec<Twist>,
    mass_and_center: Vec<(f64, Vec3)>,
    joints: usize,
}

impl GeneralizedKinematics {
    /// Borrow the original FK poses/twists for shared native contact geometry.
    #[must_use]
    pub fn native_kinematics(&self) -> &Kinematics {
        &self.kinematics
    }

    /// Number of scalar joints, excluding the six base coordinates.
    #[must_use]
    pub const fn dof_count(&self) -> usize {
        self.joints
    }
    /// Number of generalized columns, including the free base.
    #[must_use]
    pub const fn coordinate_count(&self) -> usize {
        6 + self.joints
    }
    /// Number of links in the original model order.
    #[must_use]
    pub fn link_count(&self) -> usize {
        self.kinematics.body_twist.len()
    }
    /// Actual twist in link axes; `link` must be below [`Self::link_count`].
    #[must_use]
    pub fn body_twist(&self, link: usize) -> Twist {
        self.kinematics.body_twist[link]
    }
    /// Spatial acceleration with every generalized derivative zero.
    /// `link` must be below [`Self::link_count`].
    #[must_use]
    pub fn body_bias_acceleration(&self, link: usize) -> Twist {
        self.body_bias_acceleration[link]
    }
    /// Link-body twist columns in `[omega_base, v_base, joint_rate]` order.
    /// `link` must be below [`Self::link_count`].
    #[must_use]
    pub fn body_jacobian(&self, link: usize) -> &[Twist] {
        &self.body_jacobian[link]
    }
    /// World pose of an original link; `link` must be below [`Self::link_count`].
    #[must_use]
    pub fn world_from_link(&self, link: usize) -> Se3 {
        self.kinematics.world_from_link[link]
    }

    /// Project per-link wrenches by virtual work, `Q = sum J_i^T W_i`.
    /// Both force and torque about each link origin are in that link's axes.
    /// These are external loads; internal joint efforts are added only to
    /// their compact joint entries by the caller.
    pub fn project_body_wrenches(&self, wrenches: &[Wrench]) -> Result<Vec<f64>, GeneralizedError> {
        require_len(
            "generalized.body_wrenches",
            wrenches.len(),
            self.link_count(),
        )?;
        let mut force = vec![0.0; self.coordinate_count()];
        for (jacobian, wrench) in self.body_jacobian.iter().zip(wrenches) {
            validate_wrench(*wrench, "generalized.body_wrench")?;
            for (entry, column) in force.iter_mut().zip(jacobian) {
                *entry += wrench.pairing(*column);
            }
        }
        validate_slice(&force, "generalized.projected_force")?;
        Ok(force)
    }

    /// Project world-axis forces and torques about EACH LINK ORIGIN.
    /// A torque about the world origin must first be shifted to the link
    /// origin; this API rotates axes and does not silently shift its moment.
    pub fn project_world_wrenches(
        &self,
        wrenches: &[Wrench],
    ) -> Result<Vec<f64>, GeneralizedError> {
        require_len(
            "generalized.world_wrenches",
            wrenches.len(),
            self.link_count(),
        )?;
        let local = wrenches
            .iter()
            .enumerate()
            .map(|(i, wrench)| {
                validate_wrench(*wrench, "generalized.world_wrench")?;
                let inverse = self.world_from_link(i).rotation().inverse();
                Ok(Wrench::new(
                    inverse.rotate(wrench.torque)?,
                    inverse.rotate(wrench.force)?,
                ))
            })
            .collect::<Result<Vec<_>, GeneralizedError>>()?;
        self.project_body_wrenches(&local)
    }

    /// Apply a declared uniform world acceleration to every admitted link COM.
    /// The result includes its actual base moment and scalar-joint loads.
    /// It must be added exactly once to the applied-force side of the equation.
    pub fn project_uniform_gravity(
        &self,
        gravity_world: Vec3,
    ) -> Result<Vec<f64>, GeneralizedError> {
        validate_vec3(gravity_world, "generalized.gravity_world")?;
        let local = self
            .mass_and_center
            .iter()
            .enumerate()
            .map(|(i, (mass, center))| {
                let force = self
                    .world_from_link(i)
                    .rotation()
                    .inverse()
                    .rotate(gravity_world)?
                    .scale(*mass);
                Ok(Wrench::new(center.cross(force), force))
            })
            .collect::<Result<Vec<_>, GeneralizedError>>()?;
        self.project_body_wrenches(&local)
    }
}

/// Complete dry articulated terms at one admitted state.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeFloatingTerms {
    /// Full symmetric matrix in `[omega_base, v_base, joint_rate]` ordering.
    /// All locked-body, base/joint and joint/joint couplings are retained.
    pub mass: Vec<Vec<f64>>,
    /// Left-hand inertial bias with gravity and applied forces excluded.
    pub bias: Vec<f64>,
    /// Shared native motion map for loads and additional kinetic operators.
    pub kinematics: GeneralizedKinematics,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GeneralizedError> {
    if cancelled() {
        Err(GeneralizedError::Cancelled)
    } else {
        Ok(())
    }
}

/// Extract the full free-base mass law and analytic velocity bias.
///
/// A composite rigid-body pass retains the locked mass; it never substitutes
/// the articulated Schur complement for the base block. Jacobian propagation
/// uses admitted joint subspaces, so no fictitious unit rate is checked against
/// a physical velocity limit. The existing RNEA runs once at zero generalized
/// acceleration to supply both bias and link `c` terms. No numerical derivative
/// or extra dynamics solver is used.
///
/// Work/storage are bounded by O(links*(6+joints) + (6+joints)^2), with at most
/// [`MAX_GENERALIZED_LINKS`] links and [`MAX_GENERALIZED_JOINTS`] scalar joints.
/// Cancellation is checked around the original indivisible FK/RNEA passes,
/// every 64 links in projection, and before publishing an immutable result.
/// The original configuration/rate limits and fixed-root requirement apply.
pub fn free_floating_terms(
    model: &ArticulatedModel,
    base: FreeFloatingBaseState,
    joint_position: &[f64],
    joint_velocity: &[f64],
    mut cancelled: impl FnMut() -> bool,
) -> Result<FreeFloatingTerms, GeneralizedError> {
    poll(&mut cancelled)?;
    let links = model.link_count();
    let joints = model.dof_count();
    if links > MAX_GENERALIZED_LINKS || joints > MAX_GENERALIZED_JOINTS {
        return Err(GeneralizedError::SizeLimit { links, joints });
    }
    model.free_floating_complexity()?;
    let prescribed = BaseState::prescribed(base.world_from_base, base.twist_body, Twist::zero());
    let kinematics = forward_kinematics(model, prescribed, joint_position, joint_velocity)?;
    poll(&mut cancelled)?;
    let inverse = inverse_dynamics(
        model,
        prescribed,
        joint_position,
        joint_velocity,
        &vec![0.0; joints],
        Vec3::default(),
        &vec![Wrench::default(); links],
    )?;
    poll(&mut cancelled)?;
    let columns = 6 + joints;
    let mut x_up = vec![Mat6::zero(); links];
    let mut composite = Vec::with_capacity(links);
    let mut jacobians: Vec<Vec<Twist>> = Vec::with_capacity(links);
    for (i, link) in model.links().iter().enumerate() {
        if i % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let coordinate = model.dof_index(i).map_or(0.0, |dof| joint_position[dof]);
        let parent_from_child = link
            .parent_from_child_zero()
            .compose(link.joint().motion(coordinate)?)?;
        x_up[i] = parent_from_child.inverse()?.adjoint();
        let mut jacobian = vec![Twist::zero(); columns];
        if let Some(parent) = link.parent() {
            for (column, parent_column) in jacobian.iter_mut().zip(&jacobians[parent]) {
                *column = x_up[i].apply_twist(*parent_column);
            }
        } else {
            for (column, entry) in jacobian[..6].iter_mut().enumerate() {
                let mut unit = [0.0; 6];
                unit[column] = 1.0;
                *entry = x_up[i].apply_twist(twist_from_array(unit));
            }
        }
        if let Some(dof) = model.dof_index(i) {
            jacobian[6 + dof] = link.joint().motion_subspace();
        }
        for column in &jacobian {
            validate_twist(*column, "generalized.body_jacobian")?;
        }
        jacobians.push(jacobian);
        composite.push(link.inertia().matrix()?);
    }
    let mut mass = vec![vec![0.0; columns]; columns];
    for i in (0..links).rev() {
        if i % 64 == 0 {
            poll(&mut cancelled)?;
        }
        if let Some(dof) = model.dof_index(i) {
            let column = 6 + dof;
            let mut force =
                mat6_apply_twist(&composite[i], model.links()[i].joint().motion_subspace());
            mass[column][column] = force.pairing(model.links()[i].joint().motion_subspace());
            let mut current = i;
            loop {
                force = x_up[current].transpose().apply_wrench(force);
                match model.links()[current].parent() {
                    Some(parent) => {
                        if let Some(parent_dof) = model.dof_index(parent) {
                            let value =
                                force.pairing(model.links()[parent].joint().motion_subspace());
                            mass[6 + parent_dof][column] = value;
                            mass[column][6 + parent_dof] = value;
                        }
                        current = parent;
                    }
                    None => break,
                }
            }
            for (row, value) in force.to_array().into_iter().enumerate() {
                mass[row][column] = value;
                mass[column][row] = value;
            }
        }
        let transformed = x_up[i].transpose().compose(composite[i]).compose(x_up[i]);
        validate_mat6(&transformed, "generalized.composite_inertia")?;
        if let Some(parent) = model.links()[i].parent() {
            composite[parent] = mat6_add(&composite[parent], &transformed);
            validate_mat6(&composite[parent], "generalized.composite_sum")?;
        } else {
            for row in 0..6 {
                for column in 0..=row {
                    // Both evaluations are the same symmetric bilinear form;
                    // remove only its matrix-product rounding asymmetry.
                    let value = 0.5 * transformed.m[6 * row + column]
                        + 0.5 * transformed.m[6 * column + row];
                    mass[row][column] = value;
                    mass[column][row] = value;
                }
            }
        }
    }
    for row in &mass {
        validate_slice(row, "generalized.mass")?;
    }
    let mut bias = x_up[0]
        .transpose()
        .apply_wrench(inverse.body_wrench[0])
        .to_array()
        .to_vec();
    bias.extend(inverse.generalized_force);
    validate_slice(&bias, "generalized.bias")?;
    poll(&mut cancelled)?;
    Ok(FreeFloatingTerms {
        mass,
        bias,
        kinematics: GeneralizedKinematics {
            kinematics,
            body_jacobian: jacobians,
            body_bias_acceleration: inverse.body_acceleration,
            mass_and_center: model
                .links()
                .iter()
                .map(|link| (link.inertia().mass(), link.inertia().center_of_mass()))
                .collect(),
            joints,
        },
    })
}
