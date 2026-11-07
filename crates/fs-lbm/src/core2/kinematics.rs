//! Resolved planar rotation from the native force-corrected velocity field.

use super::{Cell, Grid, MAX_REGULARIZED_BOUNDARY_SPEED_SQ};
use crate::kinematics::{FlowKinematics3, KinematicsError, KinematicsWorkspace3, checkpoint};
use fs_exec::CancelGate;

impl Grid {
    /// Full spatial velocity gradient of the accepted D2Q9 state. The planar
    /// tensor is embedded in three dimensions with homogeneous z: curl has
    /// only a z component and vortex stretching is exactly zero. Unlike the
    /// non-equilibrium stress estimate, this retains antisymmetric rotation.
    ///
    /// Allocate the reusable workspace for `[nx, ny, 1]`. Output is x-fastest,
    /// with `None` at walls. Units are lattice units; `to_physical(dt_s)` uses
    /// this state's physical timestep. Each fluid velocity is sampled once,
    /// including both gravity and spatial-force half-step corrections.
    ///
    /// Interior/periodic differences are centered. Stationary halfway walls,
    /// including nonperiodic outer faces, use a quadratic reconstruction
    /// through the wall, fluid center and opposite fluid neighbor. A one-cell
    /// periodic axis is a homogeneous extrusion; a two-cell periodic axis
    /// refuses because its two neighbors alias and cannot resolve a gradient.
    ///
    /// This adapter covers the native Fluid/Wall topology, not free surfaces
    /// or curved-link bounce-back. Gas/interface cells refuse rather than
    /// assigning fictitious zero velocity to gas or assuming a wall there.
    /// Refusals and cancellation never modify this flow owner.
    pub fn kinematics<'w>(
        &self,
        gate: &CancelGate,
        workspace: &'w mut KinematicsWorkspace3,
    ) -> Result<&'w [Option<FlowKinematics3>], KinematicsError> {
        self.kinematics_impl(gate, None, workspace)
    }

    /// Sample with the same full-grid, wall-only velocity array supplied to
    /// [`Self::stream_from_with_moving_wall_momentum`]. The prescribed value
    /// belongs to the halfway wall, not to the solid-node center. Nonperiodic
    /// exterior faces remain stationary, as in the native streaming operator.
    /// The caller must supply the wall velocities belonging to this state.
    pub fn kinematics_with_moving_walls<'w>(
        &self,
        gate: &CancelGate,
        wall_velocities: &[[f64; 2]],
        workspace: &'w mut KinematicsWorkspace3,
    ) -> Result<&'w [Option<FlowKinematics3>], KinematicsError> {
        self.kinematics_impl(gate, Some(wall_velocities), workspace)
    }

    fn kinematics_impl<'w>(
        &self,
        gate: &CancelGate,
        wall_velocities: Option<&[[f64; 2]]>,
        workspace: &'w mut KinematicsWorkspace3,
    ) -> Result<&'w [Option<FlowKinematics3>], KinematicsError> {
        workspace.admit([self.nx, self.ny, 1], gate)?;
        let cells = workspace.velocity.len();
        if self.flags.len() != cells || self.f.len() != cells || self.fext.len() != cells {
            return Err(KinematicsError::Input("planar flow field shape differs"));
        }
        if (self.periodic_x && self.nx == 2) || (self.periodic_y && self.ny == 2) {
            return Err(KinematicsError::Input(
                "two-cell periodic velocity derivative is unresolved",
            ));
        }
        if wall_velocities.is_some_and(|v| v.len() != cells) {
            return Err(KinematicsError::Input(
                "moving-wall velocity field shape differs",
            ));
        }
        if !self.g.into_iter().all(f64::is_finite) {
            return Err(KinematicsError::Input("gravity must be finite"));
        }
        for i in 0..cells {
            if i % 64 == 0 {
                checkpoint(gate)?;
            }
            if !matches!(self.flags[i], Cell::Fluid | Cell::Wall) {
                return Err(KinematicsError::Input(
                    "planar kinematics requires a Fluid/Wall topology",
                ));
            }
            if let Some(walls) = wall_velocities {
                let u = walls[i];
                if !u.into_iter().all(f64::is_finite)
                    || u[0].mul_add(u[0], u[1] * u[1]) >= MAX_REGULARIZED_BOUNDARY_SPEED_SQ
                    || (self.flags[i] != Cell::Wall && u != [0.0; 2])
                {
                    return Err(KinematicsError::Input("invalid moving-wall velocity"));
                }
            }
            if self.flags[i] == Cell::Wall {
                workspace.velocity[i] = [0.0; 3];
                continue;
            }
            // The public Grid buffers can be edited by callers. Check the
            // preconditions of moments() before entering its asserting path.
            let rho: f64 = self.f[i].iter().sum();
            if !self.f[i].into_iter().all(f64::is_finite)
                || !(rho.is_finite() && rho > 0.0)
                || !self.fext[i].into_iter().all(f64::is_finite)
            {
                return Err(KinematicsError::InvalidState {
                    cell: [i % self.nx, i / self.nx, 0],
                });
            }
            let u = self.moments(i).u;
            if !u.into_iter().all(f64::is_finite) {
                return Err(KinematicsError::InvalidState {
                    cell: [i % self.nx, i / self.nx, 0],
                });
            }
            workspace.velocity[i] = [u[0], u[1], 0.0];
        }
        for i in 0..cells {
            if i % 64 == 0 {
                checkpoint(gate)?;
            }
            if self.flags[i] == Cell::Wall {
                workspace.field[i] = None;
                continue;
            }
            let mut gradient = [[0.0; 3]; 3];
            for (axis, (minus, plus)) in [(1, 3), (2, 4)].into_iter().enumerate() {
                // Reuse the streaming topology, including periodic wrapping.
                let low = self.source(i % self.nx, i / self.nx, minus);
                let high = self.source(i % self.nx, i / self.nx, plus);
                let fluid = |j: Option<usize>| j.filter(|&j| self.flags[j] == Cell::Fluid);
                let wall = |j: Option<usize>, component: usize| {
                    j.and_then(|j| wall_velocities.map(|walls| walls[j][component]))
                        .unwrap_or(0.0)
                };
                for (component, row) in gradient.iter_mut().enumerate().take(2) {
                    let center = workspace.velocity[i][component];
                    row[axis] = match (fluid(low), fluid(high)) {
                        (Some(lo), Some(hi)) => {
                            0.5 * (workspace.velocity[hi][component]
                                - workspace.velocity[lo][component])
                        }
                        (None, Some(hi)) => {
                            let u_wall = wall(low, component);
                            (center - u_wall) + (workspace.velocity[hi][component] - u_wall) / 3.0
                        }
                        (Some(lo), None) => {
                            let u_wall = wall(high, component);
                            -(center - u_wall) - (workspace.velocity[lo][component] - u_wall) / 3.0
                        }
                        (None, None) => wall(high, component) - wall(low, component),
                    };
                }
            }
            workspace.field[i] = Some(FlowKinematics3::from_gradient(gradient)?);
        }
        checkpoint(gate)?;
        Ok(&workspace.field)
    }
}
