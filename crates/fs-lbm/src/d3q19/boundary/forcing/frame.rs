//! Prescribed rigid-frame acceleration on the original D3Q19 Guo source path.
//!
//! In frame coordinates, `a = g - a_origin - alpha cross r
//! - omega cross (omega cross r) - 2 omega cross u`. All quantities use
//! lattice units, `dx = dt = 1`. The velocity-dependent Coriolis term uses
//! an analytic three-dimensional midpoint solve; its isolated source step
//! rotates velocity without adding kinetic energy. Spatial transport, viscous
//! stress and wall momentum remain in the existing collision/stream operator.

use super::super::{BoundaryGrid3, CancelGate, Face3, Q3, TILE_CELLS, cell_force};
use super::{BoundaryForcingError3, admit_populations, poll};

/// Prescribed translation and arbitrary-axis rotation of a co-moving grid.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct RigidFrame3 {
    /// Rotation origin in lattice coordinates. Native D3Q19 cell centres are
    /// `[x+0.5, y+0.5, z+0.5]`, matching voxel geometry and wall torque origins.
    pub origin: [f64; 3],
    /// Right-handed angular velocity in radians per lattice time.
    pub angular_velocity: [f64; 3],
    /// Right-handed angular acceleration in radians per lattice time squared.
    pub angular_acceleration: [f64; 3],
    /// Inertial acceleration of the frame origin, resolved in frame axes.
    pub origin_acceleration: [f64; 3],
}

impl RigidFrame3 {
    /// Replace the owned force field for the next original lattice step.
    ///
    /// `acceleration` is physical body acceleration, such as inertial gravity,
    /// already resolved in this moving frame. Supply midpoint values of the
    /// prescribed frame motion and hold them for one unit lattice step. Grid
    /// velocities and wall velocities are RELATIVE to the frame; a co-moving
    /// stationary vessel uses the existing stationary-wall boundary conditions.
    /// This does not transform an inertial initial velocity into frame axes.
    ///
    /// Replaces BOTH uniform and per-cell forces; pressure, viscosity and wall
    /// reactions must not be added again here. The Guo half-source transform
    /// preserves the current physical mass, momentum and corrected stress.
    /// Preparation applies no impulse and advances no time. Refresh this owner
    /// before each subsequent native step to update its velocity-dependent load.
    ///
    /// Closed finite domains are admitted. A periodic direction is admitted
    /// only when both angular velocity and acceleration are parallel to it,
    /// making the position-dependent inertial force invariant under that
    /// periodic translation. Open pressure/velocity boundaries are unsupported.
    /// A rotating torus with incompatible periodic axes is refused explicitly.
    ///
    /// Existing collision scratch holds the population proposal. Two linear
    /// passes reuse an existing per-cell force allocation, allocating it only
    /// on first use. Cancellation is polled per proposal tile and immediately
    /// before the final publication pass; that deterministic, admitted pass
    /// completes without further cancellation or fallible state changes.
    ///
    /// This is prescribed-frame flow, not finite-inertia vessel dynamics or
    /// exact full-flow energy conservation. Original low-Mach, spatial/time
    /// resolution and collision-regime requirements still apply.
    pub fn prepare_step(
        self,
        grid: &mut BoundaryGrid3,
        acceleration: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        poll(gate)?;
        let denominator = self.admit(grid, acceleration)?;
        let mut initial_force_field = if grid.force_field.is_none() {
            Some(vec![[0.0; 3]; grid.f[0].len() * TILE_CELLS])
        } else {
            None
        };
        for tile in 0..grid.f[0].len() {
            poll(gate)?;
            for lane in 0..TILE_CELLS {
                let mut populations: [f64; Q3] = core::array::from_fn(|q| grid.f[q][tile].0[lane]);
                if grid.solid[tile] & (1u64 << lane) == 0 {
                    let cell = grid.coords(tile, lane);
                    admit_populations(&populations, cell)?;
                    let (previous, next, velocity) =
                        self.cell_force(grid, tile, lane, acceleration, denominator)?;
                    if previous != next {
                        crate::thermal::forcing::reframe3(
                            &mut populations,
                            velocity,
                            previous,
                            next,
                        );
                        admit_populations(&populations, cell)?;
                    }
                }
                for (proposed, value) in grid.post.iter_mut().zip(populations) {
                    proposed[tile].0[lane] = value;
                }
            }
        }
        poll(gate)?;
        // Each proposal depends only on the still-accepted populations and
        // previous force at THAT cell. Recomputing after complete admission
        // lets this commit reuse the force vector without another dense buffer.
        // The uniform force and all original populations remain fixed until
        // the entire force pass completes; earlier cells cannot affect later ones.
        for tile in 0..grid.f[0].len() {
            for lane in 0..TILE_CELLS {
                let next = if grid.solid[tile] & (1u64 << lane) == 0 {
                    self.cell_force(grid, tile, lane, acceleration, denominator)
                        .expect("rigid-frame cell was admitted against unchanged local state")
                        .1
                } else {
                    [0.0; 3]
                };
                let field = grid
                    .force_field
                    .as_mut()
                    .or(initial_force_field.as_mut())
                    .expect("original or first-use force allocation exists");
                field[tile * TILE_CELLS + lane] = next;
            }
        }
        core::mem::swap(&mut grid.f, &mut grid.post);
        if let Some(field) = initial_force_field {
            grid.force_field = Some(field);
        }
        grid.force = [0.0; 3];
        grid.force_field_shape = None;
        Ok(())
    }

    fn admit(
        self,
        grid: &BoundaryGrid3,
        acceleration: [f64; 3],
    ) -> Result<f64, BoundaryForcingError3> {
        use BoundaryForcingError3::Invalid;
        let denominator = self.midpoint_denominator(acceleration).map_err(Invalid)?;
        if !grid.collision_model.supports_body_force() {
            return Err(Invalid(
                "rigid-frame flow requires a forcing-compatible collision model",
            ));
        }
        for face in Face3::ALL {
            let boundary = grid.boundaries.face(face);
            if boundary.is_open() {
                return Err(Invalid(
                    "rigid-frame flow requires closed or compatible periodic boundaries",
                ));
            }
            if boundary.is_periodic() {
                let axis = face.axis();
                if (0..3).any(|other| {
                    other != axis
                        && (self.angular_velocity[other] != 0.0
                            || self.angular_acceleration[other] != 0.0)
                }) {
                    return Err(Invalid(
                        "rigid-frame inertial forces are not periodic along this axis",
                    ));
                }
            }
        }
        Ok(denominator)
    }

    fn cell_force(
        self,
        grid: &BoundaryGrid3,
        tile: usize,
        lane: usize,
        acceleration: [f64; 3],
        denominator: f64,
    ) -> Result<([f64; 3], [f64; 3], [f64; 3]), BoundaryForcingError3> {
        let cell = grid.coords(tile, lane);
        // Admission precedes raw_moments' positive-density assertion. During
        // publication these original populations have not changed.
        let (rho, momentum) = grid.raw_moments(cell[0], cell[1], cell[2]);
        let previous = cell_force(grid.force, grid.force_field.as_deref(), tile, lane);
        let velocity = core::array::from_fn(|axis| (momentum[axis] + 0.5 * previous[axis]) / rho);
        let position = cell.map(|coordinate| coordinate as f64 + 0.5);
        let effective = self.effective_acceleration(position, velocity, acceleration, denominator);
        let next = effective.map(|value| rho * value);
        if !previous
            .into_iter()
            .chain(velocity)
            .chain(next)
            .all(f64::is_finite)
        {
            return Err(BoundaryForcingError3::Unphysical { cell });
        }
        Ok((previous, next, velocity))
    }

    // Shared by dense and sparse free-surface owners. Admission is once per
    // source update; each caller separately validates its local state/result.
    pub(crate) fn midpoint_denominator(self, acceleration: [f64; 3]) -> Result<f64, &'static str> {
        let denominator = 1.0 + dot(self.angular_velocity, self.angular_velocity);
        if !self
            .origin
            .into_iter()
            .chain(self.angular_velocity)
            .chain(self.angular_acceleration)
            .chain(self.origin_acceleration)
            .chain(acceleration)
            .all(f64::is_finite)
            || !denominator.is_finite()
        {
            return Err("rigid-frame motion must be finite and representable");
        }
        Ok(denominator)
    }

    pub(crate) fn effective_acceleration(
        self,
        position: [f64; 3],
        velocity: [f64; 3],
        acceleration: [f64; 3],
        denominator: f64,
    ) -> [f64; 3] {
        let r = core::array::from_fn(|axis| position[axis] - self.origin[axis]);
        let euler = cross(self.angular_acceleration, r);
        let centrifugal = cross(self.angular_velocity, cross(self.angular_velocity, r));
        let base = core::array::from_fn(|axis| {
            acceleration[axis] - self.origin_acceleration[axis] - euler[axis] - centrifugal[axis]
        });
        midpoint_acceleration(self.angular_velocity, velocity, base, denominator)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}

fn midpoint_acceleration(
    omega: [f64; 3],
    velocity: [f64; 3],
    base: [f64; 3],
    denominator: f64,
) -> [f64; 3] {
    // (I + [omega]x) u_mid = u_start + a_non_coriolis/2.
    // Its inverse is (I - [omega]x + omega*omega^T)/(1 + |omega|^2).
    // Unlike an explicit Coriolis kick, this skew solve does no midpoint work.
    let b = core::array::from_fn(|axis| velocity[axis] + 0.5 * base[axis]);
    let skew = cross(omega, b);
    let projection = dot(omega, b);
    let midpoint = core::array::from_fn(|axis| {
        omega[axis].mul_add(projection, b[axis] - skew[axis]) / denominator
    });
    let coriolis = cross(omega, midpoint);
    core::array::from_fn(|axis| base[axis] - 2.0 * coriolis[axis])
}
