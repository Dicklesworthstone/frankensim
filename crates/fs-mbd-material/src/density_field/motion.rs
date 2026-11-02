//! Exact P1 material momentum handed to the existing rigid-body state owner.

use fs_ga::{Se3, Twist, Vec3};
use fs_mbd::RigidBodyState;
use fs_mbd::articulated::SpatialInertia;

use super::{DeformedDensityBody, TetrahedralDensityBody};
use crate::principal::{PrincipalFrame, PrincipalFrameError};

/// Numerical admission and the explicit cost of discarding nonrigid motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidProjectionControls {
    /// Relative tensor, spin-frame and centred kinetic-energy residual limit,
    /// in (0, 1e-8]. This is a numerical tolerance, not material uncertainty.
    pub relative_tolerance: f64,
    /// Maximum unresolved internal kinetic energy [J]. This nonnegative energy
    /// is returned separately; projection does not make it caloric heat.
    pub maximum_unresolved_kinetic_energy_j: f64,
}

/// A proposed native rigid state retaining its original material inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMotionProjection<'a> {
    reference: &'a TetrahedralDensityBody,
    /// Current full inertia and its original geometry-to-principal transform.
    pub frame: PrincipalFrame,
    /// Actual native state accepted by free-body stepping and point impulses.
    pub state: RigidBodyState,
    /// Kinetic energy of the supplied nodal velocity field [J].
    pub nodal_kinetic_energy_j: f64,
    /// Kinetic energy of the returned native rigid state [J].
    pub rigid_kinetic_energy_j: f64,
    /// Integral of rho |u - V_COM - omega cross r|^2 / 2 [J].
    /// A composing physical handoff must retain or explicitly dispose of this
    /// account. It is neither numerical damping nor an inferred heat transfer.
    pub unresolved_kinetic_energy_j: f64,
    /// Relative spin error through the actual returned native pose.
    pub relative_momentum_residual: f64,
    /// Centred kinetic-energy split error divided by its largest term.
    pub relative_energy_residual: f64,
}

impl<'a> RigidMotionProjection<'a> {
    /// Original source states, density receipts, reference mesh and inventory.
    pub const fn reference(&self) -> &'a TetrahedralDensityBody {
        self.reference
    }
}

impl TetrahedralDensityBody {
    /// Project an actual P1 nodal velocity field into native free-body motion.
    ///
    /// Velocities are absolute physical velocities [m/s], expressed in the
    /// geometry/body axes. `world_from_body` places that geometry and rotates
    /// those components; no additional frame velocity or twist is inferred.
    /// The density-weighted linear and all three angular momenta are retained.
    /// The best rigid velocity is V = P/M and omega = I_COM^-1 L_COM, using the
    /// existing full-inertia principal-frame owner.
    ///
    /// Exact positive cubic tetrahedral quadrature integrates the P1 density
    /// times P1 position/velocity products. Internal kinetic energy is measured
    /// directly from the residual velocity and must fit the caller's bound.
    /// This is an explicit model projection, not a solidification/fracture law,
    /// a finite-strain solve, or an elastic/thermal energy transfer. Geometry
    /// remains at the supplied configuration; no existing state is mutated.
    /// Work is O(vertices + cells), with constant scratch space and cancellation
    /// before work, per vertex/cell and before returning the complete candidate.
    pub fn project_rigid_motion(
        &self,
        world_from_body: Se3,
        nodal_velocity_in_body_axes: &[Vec3],
        controls: RigidProjectionControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RigidMotionProjection<'_>, PrincipalFrameError> {
        project(
            self,
            self.inertia,
            |vertex| self.vertices[vertex].position_m,
            world_from_body,
            nodal_velocity_in_body_axes,
            controls,
            cancelled,
        )
    }
}

impl<'a> DeformedDensityBody<'a> {
    /// The same momentum projection on the actual deformed configuration.
    ///
    /// Original per-cell masses and density fractions integrate the current
    /// P1 positions and velocities. Thus rho_current = rho_reference/J enters
    /// without re-querying a source or averaging unequal cell Jacobians. The
    /// velocity convention, explicit unresolved-energy bound and no-claim
    /// limits are those of [`TetrahedralDensityBody::project_rigid_motion`].
    pub fn project_rigid_motion(
        &self,
        world_from_body: Se3,
        nodal_velocity_in_body_axes: &[Vec3],
        controls: RigidProjectionControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RigidMotionProjection<'a>, PrincipalFrameError> {
        project(
            self.reference,
            self.inertia,
            |vertex| self.vertices_m[vertex],
            world_from_body,
            nodal_velocity_in_body_axes,
            controls,
            cancelled,
        )
    }
}

fn project<'a>(
    reference: &'a TetrahedralDensityBody,
    inertia: SpatialInertia,
    position: impl Fn(usize) -> Vec3,
    world_from_body: Se3,
    velocities: &[Vec3],
    controls: RigidProjectionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RigidMotionProjection<'a>, PrincipalFrameError> {
    poll(&mut cancelled)?;
    if !controls.relative_tolerance.is_finite()
        || controls.relative_tolerance <= 0.0
        || controls.relative_tolerance > 1e-8
        || !controls.maximum_unresolved_kinetic_energy_j.is_finite()
        || controls.maximum_unresolved_kinetic_energy_j < 0.0
        || velocities.len() != reference.vertices.len()
    {
        return Err(invalid(
            "one velocity per vertex, tolerance in (0,1e-8] and finite nonnegative energy bound required",
        ));
    }
    for velocity in velocities {
        poll(&mut cancelled)?;
        if !finite(*velocity) {
            return Err(invalid("nodal velocities must be finite"));
        }
    }
    let frame = PrincipalFrame::from_inertia(inertia, controls.relative_tolerance, &mut cancelled)?;
    // Remove a common boost before taking moments. In particular, a large
    // translational energy cannot erase a small unresolved internal energy.
    let anchor = velocities[reference.tetrahedra[0][0]];
    let mut mean_delta = Vec3::new(0.0, 0.0, 0.0);
    samples(
        reference,
        &position,
        velocities,
        anchor,
        inertia,
        &mut cancelled,
        |mass, _, delta| mean_delta = mean_delta + delta.scale(mass / inertia.mass()),
    )?;
    let velocity_com = anchor + mean_delta;
    if !finite(velocity_com) {
        return Err(invalid("centre velocity is unrepresentable"));
    }
    let mut spin_body = Vec3::new(0.0, 0.0, 0.0);
    let mut centred_energy = 0.0;
    samples(
        reference,
        &position,
        velocities,
        anchor,
        inertia,
        &mut cancelled,
        |mass, arm, delta| {
            let weighted_velocity = (delta - mean_delta).scale(mass.sqrt());
            spin_body = spin_body + arm.scale(mass.sqrt()).cross(weighted_velocity);
            centred_energy += 0.5 * weighted_velocity.dot(weighted_velocity);
        },
    )?;
    let pose = frame
        .state_from_kinematics(world_from_body, Twist::zero())?
        .pose();
    let spin_principal = frame
        .body_from_principal()
        .rotation()
        .inverse()
        .rotate(spin_body)?;
    let momentum_world = world_from_body
        .rotation()
        .rotate(velocity_com)?
        .scale(inertia.mass());
    let state = RigidBodyState::new(pose, mbd(momentum_world), mbd(spin_principal))?;
    let omega_principal = frame
        .properties()
        .angular_velocity_body_checked(state.angular_momentum_body())?;
    let omega_body = frame
        .body_from_principal()
        .rotation()
        .rotate(ga(omega_principal))?;
    let rotational_energy = 0.5 * state.angular_momentum_body().dot(omega_principal);
    let rigid_kinetic_energy_j = state.kinetic_energy(frame.properties())?;
    let mut unresolved_kinetic_energy_j = 0.0;
    samples(
        reference,
        &position,
        velocities,
        anchor,
        inertia,
        &mut cancelled,
        |mass, arm, delta| {
            let residual = (delta - mean_delta - omega_body.cross(arm)).scale(mass.sqrt());
            unresolved_kinetic_energy_j += 0.5 * residual.dot(residual);
        },
    )?;
    let expected_spin = world_from_body.rotation().rotate(spin_body)?;
    let actual_spin = ga(pose
        .orientation()
        .rotate_body_to_world(state.angular_momentum_body()));
    let spin_norm = norm(expected_spin);
    let relative_momentum_residual = if spin_norm == 0.0 {
        norm(actual_spin)
    } else {
        norm(actual_spin - expected_spin) / spin_norm
    };
    let energy_scale = centred_energy
        .max(rotational_energy)
        .max(unresolved_kinetic_energy_j);
    let relative_energy_residual = if energy_scale == 0.0 {
        0.0
    } else {
        (centred_energy / energy_scale
            - rotational_energy / energy_scale
            - unresolved_kinetic_energy_j / energy_scale)
            .abs()
    };
    let translation = 0.5
        * state
            .linear_momentum_world()
            .dot(state.center_of_mass_velocity_world(frame.properties())?);
    let nodal_kinetic_energy_j = translation + centred_energy;
    if !centred_energy.is_finite()
        || !nodal_kinetic_energy_j.is_finite()
        || !relative_momentum_residual.is_finite()
        || relative_momentum_residual > controls.relative_tolerance
        || !relative_energy_residual.is_finite()
        || relative_energy_residual > controls.relative_tolerance
    {
        return Err(invalid(
            "rigid projection exceeded its numerical residual bound",
        ));
    }
    if !unresolved_kinetic_energy_j.is_finite()
        || unresolved_kinetic_energy_j > controls.maximum_unresolved_kinetic_energy_j
    {
        return Err(invalid(
            "rigid projection exceeded its unresolved kinetic-energy bound",
        ));
    }
    poll(&mut cancelled)?;
    Ok(RigidMotionProjection {
        reference,
        frame,
        state,
        nodal_kinetic_energy_j,
        rigid_kinetic_energy_j,
        unresolved_kinetic_energy_j,
        relative_momentum_residual,
        relative_energy_residual,
    })
}

// Cubic-exact positive rule: vertices have volume weights 1/40 and opposite
// face centroids 9/40. With f_i = rho_i / sum(rho), their MASS weights are
// m*f_i/10 and 3*m*sum_{j != i}(f_j)/10. Integrands after absorbing density
// are degree <= 2; no negative kinetic-energy terms or centroid-density alias.
fn samples(
    reference: &TetrahedralDensityBody,
    position: &impl Fn(usize) -> Vec3,
    velocities: &[Vec3],
    anchor: Vec3,
    inertia: SpatialInertia,
    cancelled: &mut impl FnMut() -> bool,
    mut consume: impl FnMut(f64, Vec3, Vec3),
) -> Result<(), PrincipalFrameError> {
    for (tet, cell) in reference.tetrahedra.iter().zip(&reference.cells) {
        poll(cancelled)?;
        let scale = cell
            .vertex_density_kg_m3
            .into_iter()
            .fold(0.0_f64, f64::max);
        let scaled = cell.vertex_density_kg_m3.map(|rho| rho / scale);
        let sum = scaled.into_iter().sum::<f64>();
        let fractions = scaled.map(|rho| rho / sum);
        let arms = tet.map(|vertex| position(vertex) - inertia.center_of_mass());
        let deltas = tet.map(|vertex| velocities[vertex] - anchor);
        let mass = cell.inertia.mass();
        for i in 0..4 {
            consume((mass * 0.1) * fractions[i], arms[i], deltas[i]);
            let [a, b, c] = [(i + 1) % 4, (i + 2) % 4, (i + 3) % 4];
            let arm =
                arms[a].scale(1.0 / 3.0) + arms[b].scale(1.0 / 3.0) + arms[c].scale(1.0 / 3.0);
            let delta = deltas[a].scale(1.0 / 3.0)
                + deltas[b].scale(1.0 / 3.0)
                + deltas[c].scale(1.0 / 3.0);
            consume(
                (mass * 0.3) * (fractions[a] + fractions[b] + fractions[c]),
                arm,
                delta,
            );
        }
    }
    Ok(())
}

fn mbd(v: Vec3) -> fs_mbd::Vec3 {
    fs_mbd::Vec3::new(v.x, v.y, v.z)
}
fn ga(v: fs_mbd::Vec3) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn finite(v: Vec3) -> bool {
    [v.x, v.y, v.z].into_iter().all(f64::is_finite)
}
fn invalid(message: &'static str) -> PrincipalFrameError {
    PrincipalFrameError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PrincipalFrameError> {
    if cancelled() {
        Err(PrincipalFrameError::Cancelled)
    } else {
        Ok(())
    }
}
