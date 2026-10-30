//! Moving immersed geometry on the native six-velocity fluid/body operator.
//!
//! Free translation and the native momentum-preserving SO(3) map advance half
//! a lattice tick, geometry is sampled in that pose, and the existing implicit
//! midpoint exchange advances velocities. A second free half-step completes
//! the pose. All endpoint admission and diagnostics precede fluid publication.

use fs_ga::{So3, Vec3};
use fs_time::lie::{LieStepError, gyrostat::InertiaTensor3};

use super::partial_saturation::{checkpoint, cross, dot, invalid, low_mach, point_velocity};
use super::{
    BoundaryGrid3, CancelGate, PartialSaturationBody3, PartialSaturationBodyExchange3,
    PartialSaturationError3, TILE_CELLS,
};

type Result<T> = core::result::Result<T, PartialSaturationError3>;

mod multiple;
pub use multiple::{MovingLubricatedPlanesExchange3, MovingPartialSaturationBodiesExchange3};

/// Fixed geometry expressed relative to the body's center in BODY axes.
///
/// Both queries must be deterministic and unchanged during a step. All points
/// for which `contains` is true must lie within the declared bounding sphere.
/// The bound is a caller obligation for custom shapes; it is used to refuse
/// unsupported wall contact and periodic body seam crossing conservatively.
pub trait BodyFrameShape3 {
    /// Finite positive radius enclosing the shape about the body origin.
    fn bounding_radius(&self) -> f64;
    /// Membership at a finite point in body coordinates, including the surface.
    fn contains(&self, point: Vec3) -> bool;
}

/// Centered ellipsoid with fixed positive semiaxes in the body frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyFrameEllipsoid3 {
    semiaxes: Vec3,
}

impl BodyFrameEllipsoid3 {
    /// Admit finite positive semiaxes in lattice-length units.
    ///
    /// # Errors
    /// Refuses non-finite or non-positive semiaxes.
    pub fn new(semiaxes: Vec3) -> Result<Self> {
        if !array(semiaxes)
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
        {
            return Err(invalid("ellipsoid semiaxes must be positive and finite"));
        }
        Ok(Self { semiaxes })
    }

    /// Admit a sphere; the same geometry sampler is used for both primitives.
    ///
    /// # Errors
    /// Refuses a non-finite or non-positive radius.
    pub fn sphere(radius: f64) -> Result<Self> {
        Self::new(Vec3::new(radius, radius, radius))
    }

    /// Fixed semiaxes in body coordinates.
    #[must_use]
    pub const fn semiaxes(self) -> Vec3 {
        self.semiaxes
    }
}

impl BodyFrameShape3 for BodyFrameEllipsoid3 {
    fn bounding_radius(&self) -> f64 {
        self.semiaxes.x.max(self.semiaxes.y).max(self.semiaxes.z)
    }

    fn contains(&self, point: Vec3) -> bool {
        let coordinates = array(point);
        let axes = array(self.semiaxes);
        if (0..3).any(|i| !coordinates[i].is_finite() || coordinates[i].abs() > axes[i]) {
            return false;
        }
        let scaled = core::array::from_fn(|i| coordinates[i] / axes[i]);
        dot(scaled, scaled) <= 1.0
    }
}

/// Reusable deterministic subcell quadrature; one fraction buffer per body.
///
/// Each cell uses `samples_per_axis^3` equally spaced midpoint samples. Values
/// are raw occupancy counts, without volume renormalization or density repair.
/// This is a finite-resolution volume estimate, not an exact intersection or
/// a geometric-conservation certificate. Only this scratch buffer changes on
/// rejected motion steps; the next call clears and rebuilds it.
#[derive(Debug)]
pub struct BodyFrameSampler3 {
    dims: [usize; 3],
    samples_per_axis: usize,
    fractions: Vec<f64>,
}

impl BodyFrameSampler3 {
    /// Allocate one reusable field for this grid; admit 1 through 8 samples
    /// along each cell axis (at most 512 membership queries per cell).
    ///
    /// # Errors
    /// Refuses an unsupported sample count, size overflow or allocation failure.
    pub fn new(grid: &BoundaryGrid3, samples_per_axis: usize) -> Result<Self> {
        if !(1..=8).contains(&samples_per_axis) {
            return Err(invalid(
                "body quadrature requires 1 through 8 samples per axis",
            ));
        }
        let dims = grid.dimensions();
        let count = dims[0]
            .checked_mul(dims[1])
            .and_then(|xy| xy.checked_mul(dims[2]))
            .ok_or_else(|| invalid("body quadrature field size overflowed"))?;
        let mut fractions = Vec::new();
        fractions
            .try_reserve_exact(count)
            .map_err(|_| invalid("body quadrature allocation refused"))?;
        fractions.resize(count, 0.0);
        Ok(Self {
            dims,
            samples_per_axis,
            fractions,
        })
    }

    /// Last proposed row-major fractions, including scratch from a refusal.
    #[must_use]
    pub fn fractions(&self) -> &[f64] {
        &self.fractions
    }

    /// Sample a pose independently of motion. Geometry outside the domain is
    /// clipped; the motion step separately refuses support near the domain halo.
    ///
    /// # Errors
    /// Refuses non-finite geometry, invalid bounds, arithmetic or cancellation.
    pub fn sample(
        &mut self,
        center: [f64; 3],
        orientation: So3,
        shape: &impl BodyFrameShape3,
        gate: &CancelGate,
    ) -> Result<f64> {
        sample_geometry(
            self.dims,
            self.samples_per_axis,
            center,
            orientation,
            shape,
            None,
            Some(&mut self.fractions),
            gate,
        )
    }
}

/// One moving rigid body whose inertia is ADDITIONAL to retained interior fluid.
/// All lengths and times are lattice units, `dx = dt = 1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MovingPartialSaturationBody3 {
    /// World center of mass; cell centers are `[x+0.5,y+0.5,z+0.5]`.
    pub center: [f64; 3],
    /// Proper rotation from the fixed geometry/body frame to WORLD axes.
    pub orientation: So3,
    /// Finite positive additional solid mass, held fixed throughout the step.
    pub mass: f64,
    /// Fixed body-frame inertia, admitted once by the native rotation owner.
    pub inertia: InertiaTensor3,
    /// World translational velocity.
    pub velocity: [f64; 3],
    /// Spin angular momentum about the center in the rotating BODY frame.
    /// No hidden internal rotor momentum is included in this body model.
    pub angular_momentum_body: Vec3,
}

impl MovingPartialSaturationBody3 {
    /// World angular velocity from the current body momentum and orientation.
    ///
    /// # Errors
    /// Refuses non-finite or unrepresentable native inertia/group operations.
    pub fn angular_velocity_world(self) -> Result<[f64; 3]> {
        Ok(array(native(self.orientation.rotate(native(
            self.inertia.velocity(self.angular_momentum_body),
        )?))?))
    }
}

/// Actual accepted body trajectory and native fluid/body exchange for one tick.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct MovingPartialSaturationExchange3 {
    /// Existing midpoint exchange, including separately measured fixed walls.
    pub coupling: PartialSaturationBodyExchange3,
    /// Configuration at which geometry, world inertia and force arms were used.
    pub midpoint_center: [f64; 3],
    /// Body-to-world rotation at the coupling configuration.
    pub midpoint_orientation: So3,
    /// Raw sampled solid volumes at the initial, midpoint and final poses.
    /// Variation is quadrature error; no mass is inserted or removed to hide it.
    pub sampled_solid_volume: [f64; 3],
    /// External angular impulse about the requested origin, including the
    /// supplied center force's arm at the midpoint configuration.
    pub external_angular_impulse: [f64; 3],
    /// Actual endpoint total body linear/angular momentum change minus external
    /// and measured hydrodynamic impulses, about the requested fixed origin.
    /// Periodic fluid angular transport is separate from this body account.
    pub momentum_residual: [f64; 6],
    /// Actual endpoint kinetic-energy change, including native free drift.
    pub kinetic_energy_change: f64,
    /// Rotational energy change in both native free half-steps. The free map
    /// preserves spatial momentum; its finite-step kinetic energy is not exact.
    pub free_rotation_energy_change: f64,
    /// Endpoint kinetic change minus hydrodynamic/external work and separately
    /// reported free-rotation energy change. Measures frame/roundoff defects.
    pub work_residual: f64,
}

impl BoundaryGrid3 {
    /// Advance one body's pose, momenta and actual native populations together.
    ///
    /// The split is free half-drift, implicit six-velocity midpoint exchange,
    /// free half-drift. Only midpoint fractions drive collision. External input
    /// is WORLD linear impulse at the midpoint center, then WORLD center couple.
    /// The body's stored angular momentum is transported natively rather than
    /// repeatedly reconstructed from its angular velocity.
    ///
    /// # Errors
    /// Inherits the BGK/held-Guo-force and low-Mach admission of the native coupling.
    /// Both swept enclosing balls must stay in the interior domain halo and
    /// avoid fixed solid voxels; contact and periodic body wrapping are refused.
    /// Native drift, endpoint geometry, receipt and cancellation failures leave
    /// fluid, body and topology admission unchanged. Scratch may change.
    ///
    /// This is one rigid immersed body with retained interior fluid, not a
    /// displaced-fluid density closure, moving-mask remap, GCL or calibrated drag
    /// model. No full coupled temporal-order claim is made.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)] // One atomic split step.
    pub fn step_with_moving_partial_saturation_body(
        &mut self,
        sampler: &mut BodyFrameSampler3,
        shape: &impl BodyFrameShape3,
        body: &mut MovingPartialSaturationBody3,
        external_impulse: [f64; 6],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MovingPartialSaturationExchange3> {
        checkpoint(gate)?;
        if !(sampler.dims == self.dimensions()
            && body.mass.is_finite()
            && body.mass > 0.0
            && body
                .center
                .into_iter()
                .chain(external_impulse)
                .chain(moment_origin)
                .all(f64::is_finite)
            && low_mach(body.velocity))
        {
            return Err(invalid(
                "moving body, sampler, impulse or origin is inadmissible",
            ));
        }
        let initial = *body;
        let radius = shape_radius(shape)?;
        let midpoint_center =
            core::array::from_fn(|i| initial.center[i] + 0.5 * initial.velocity[i]);
        admit_sweep(self, initial.center, midpoint_center, radius, gate)?;
        let (midpoint_orientation, midpoint_momentum) = native(initial.inertia.momentum_step(
            initial.orientation,
            initial.angular_momentum_body,
            Vec3::new(0.0, 0.0, 0.0),
            0.5,
        ))?;
        let omega_initial = initial.angular_velocity_world()?;
        let omega_before_exchange = array(native(
            midpoint_orientation.rotate(native(initial.inertia.velocity(midpoint_momentum))?),
        )?);
        let volume_before = sample_geometry(
            sampler.dims,
            sampler.samples_per_axis,
            initial.center,
            initial.orientation,
            shape,
            Some((initial.velocity, omega_initial)),
            None,
            gate,
        )?;
        let volume_midpoint = sample_geometry(
            sampler.dims,
            sampler.samples_per_axis,
            midpoint_center,
            midpoint_orientation,
            shape,
            Some((initial.velocity, omega_before_exchange)),
            Some(&mut sampler.fractions),
            gate,
        )?;
        let rotation = midpoint_orientation.matrix();
        let rotated = rotation
            .compose(initial.inertia.matrix())
            .compose(rotation.transpose());
        // Native coupling admits exact symmetry; compute one triangle and
        // mirror its rounded values, rather than hiding an arbitrary asymmetry.
        let inertia_world =
            core::array::from_fn(|i| core::array::from_fn(|j| rotated.m[3 * i.max(j) + i.min(j)]));
        let mut frozen = PartialSaturationBody3 {
            center: midpoint_center,
            mass: initial.mass,
            inertia_world,
            velocity: initial.velocity,
            angular_velocity: omega_before_exchange,
        };
        let dims = sampler.dims;
        let samples = sampler.samples_per_axis;
        let (_, (accepted, receipt)) = self.step_with_partial_saturation_body_checked(
            &sampler.fractions,
            &mut frozen,
            external_impulse,
            moment_origin,
            gate,
            |grid, exchanged, coupling| {
                let dw = core::array::from_fn(|i| {
                    exchanged.angular_velocity[i] - omega_before_exchange[i]
                });
                let spin_impulse = inertia_world.map(|row| dot(row, dw));
                let after_exchange_momentum = midpoint_momentum
                    + native(midpoint_orientation.inverse().rotate(vector(spin_impulse)))?;
                let (orientation, angular_momentum_body) = native(initial.inertia.momentum_step(
                    midpoint_orientation,
                    after_exchange_momentum,
                    Vec3::new(0.0, 0.0, 0.0),
                    0.5,
                ))?;
                let accepted = MovingPartialSaturationBody3 {
                    center: core::array::from_fn(|i| {
                        midpoint_center[i] + 0.5 * exchanged.velocity[i]
                    }),
                    orientation,
                    velocity: exchanged.velocity,
                    angular_momentum_body,
                    ..initial
                };
                admit_sweep(grid, midpoint_center, accepted.center, radius, gate)?;
                let volume_after = sample_geometry(
                    dims,
                    samples,
                    accepted.center,
                    accepted.orientation,
                    shape,
                    Some((accepted.velocity, accepted.angular_velocity_world()?)),
                    None,
                    gate,
                )?;
                let external_linear = external_impulse[..3].try_into().unwrap();
                let arm = core::array::from_fn(|i| midpoint_center[i] - moment_origin[i]);
                let external_orbital = cross(arm, external_linear);
                let external_angular_impulse =
                    core::array::from_fn(|i| external_impulse[i + 3] + external_orbital[i]);
                let before_angular = total_angular_momentum(initial, moment_origin)?;
                let after_angular = total_angular_momentum(accepted, moment_origin)?;
                let momentum_residual = core::array::from_fn(|i| {
                    if i < 3 {
                        initial.mass * (accepted.velocity[i] - initial.velocity[i])
                            - external_impulse[i]
                            - coupling.hydrodynamic.body_impulse[i]
                    } else {
                        after_angular[i - 3]
                            - before_angular[i - 3]
                            - external_angular_impulse[i - 3]
                            - coupling.hydrodynamic.body_angular_impulse[i - 3]
                    }
                });
                let spin_before =
                    rotational_energy(initial.inertia, initial.angular_momentum_body)?;
                let spin_midpoint = rotational_energy(initial.inertia, midpoint_momentum)?;
                let spin_exchanged = rotational_energy(initial.inertia, after_exchange_momentum)?;
                let spin_after = rotational_energy(initial.inertia, angular_momentum_body)?;
                let kinetic_energy_change = 0.5
                    * initial.mass
                    * (dot(accepted.velocity, accepted.velocity)
                        - dot(initial.velocity, initial.velocity))
                    + spin_after
                    - spin_before;
                let free_rotation_energy_change =
                    (spin_midpoint - spin_before) + (spin_after - spin_exchanged);
                let work_residual = kinetic_energy_change
                    - coupling.external_work
                    - coupling.hydrodynamic.body_work
                    - free_rotation_energy_change;
                if ![
                    volume_before,
                    volume_midpoint,
                    volume_after,
                    kinetic_energy_change,
                    free_rotation_energy_change,
                    work_residual,
                ]
                .into_iter()
                .chain(external_angular_impulse)
                .chain(momentum_residual)
                .all(f64::is_finite)
                {
                    return Err(invalid("moving body momentum or energy account overflowed"));
                }
                Ok((
                    accepted,
                    MovingPartialSaturationExchange3 {
                        coupling: *coupling,
                        midpoint_center,
                        midpoint_orientation,
                        sampled_solid_volume: [volume_before, volume_midpoint, volume_after],
                        external_angular_impulse,
                        momentum_residual,
                        kinetic_energy_change,
                        free_rotation_energy_change,
                        work_residual,
                    },
                ))
            },
        )?;
        *body = accepted;
        Ok(receipt)
    }
}

fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn native<T>(value: core::result::Result<T, impl Into<LieStepError>>) -> Result<T> {
    value.map_err(|source| PartialSaturationError3::Motion {
        source: source.into(),
    })
}

fn shape_radius(shape: &(impl BodyFrameShape3 + ?Sized)) -> Result<f64> {
    let radius = shape.bounding_radius();
    if !(radius.is_finite() && radius > 0.0) {
        return Err(invalid(
            "body shape requires a positive finite bounding radius",
        ));
    }
    Ok(radius)
}

fn rotational_energy(inertia: InertiaTensor3, momentum: Vec3) -> Result<f64> {
    let energy = 0.5 * dot(array(momentum), array(native(inertia.velocity(momentum))?));
    if !energy.is_finite() {
        return Err(invalid("body rotational energy overflowed"));
    }
    Ok(energy)
}

fn total_angular_momentum(
    body: MovingPartialSaturationBody3,
    origin: [f64; 3],
) -> Result<[f64; 3]> {
    let spin = array(native(body.orientation.rotate(body.angular_momentum_body))?);
    let arm = core::array::from_fn(|i| body.center[i] - origin[i]);
    let orbital = cross(arm, body.velocity.map(|v| body.mass * v));
    Ok(core::array::from_fn(|i| spin[i] + orbital[i]))
}

fn admit_sweep(
    grid: &BoundaryGrid3,
    before: [f64; 3],
    after: [f64; 3],
    radius: f64,
    gate: &CancelGate,
) -> Result<()> {
    let dims = grid.dimensions();
    let lower: [f64; 3] = core::array::from_fn(|i| before[i].min(after[i]) - radius);
    let upper: [f64; 3] = core::array::from_fn(|i| before[i].max(after[i]) + radius);
    if (0..3).any(|i| {
        !before[i].is_finite()
            || !after[i].is_finite()
            || !lower[i].is_finite()
            || !upper[i].is_finite()
            || lower[i] < 1.0
            || upper[i] > dims[i] as f64 - 1.0
    }) {
        return Err(invalid("moving body swept support reaches the domain halo"));
    }
    let low = lower.map(|x| x.floor() as usize);
    let high = upper.map(|x| x.ceil() as usize);
    let mut visited = 0usize;
    for z in low[2]..high[2] {
        for y in low[1]..high[1] {
            for x in low[0]..high[0] {
                if visited.is_multiple_of(TILE_CELLS) {
                    checkpoint(gate)?;
                }
                visited += 1;
                if grid.is_solid(x, y, z) {
                    return Err(PartialSaturationError3::Cell {
                        cell: [x, y, z],
                        reason: "moving body swept support meets a fixed solid cell",
                    });
                }
            }
        }
    }
    checkpoint(gate)
}

#[allow(clippy::too_many_arguments)] // Shared quadrature for three poses and public sampling.
fn sample_geometry(
    dims: [usize; 3],
    samples: usize,
    center: [f64; 3],
    orientation: So3,
    shape: &(impl BodyFrameShape3 + ?Sized),
    velocity: Option<([f64; 3], [f64; 3])>,
    mut fractions: Option<&mut [f64]>,
    gate: &CancelGate,
) -> Result<f64> {
    checkpoint(gate)?;
    let radius = shape_radius(shape)?;
    if center
        .into_iter()
        .any(|c| !c.is_finite() || !(c - radius).is_finite() || !(c + radius).is_finite())
    {
        return Err(invalid("body sampling center or support overflowed"));
    }
    if let Some(field) = fractions.as_mut() {
        for tile in field.chunks_mut(TILE_CELLS) {
            checkpoint(gate)?;
            tile.fill(0.0);
        }
    }
    let low: [usize; 3] =
        core::array::from_fn(|i| (center[i] - radius).floor().clamp(0.0, dims[i] as f64) as usize);
    let high: [usize; 3] =
        core::array::from_fn(|i| (center[i] + radius).ceil().clamp(0.0, dims[i] as f64) as usize);
    let inverse = orientation.matrix().transpose();
    let denominator = (samples * samples * samples) as f64;
    let mut volume = 0.0;
    let mut visited = 0usize;
    for z in low[2]..high[2] {
        for y in low[1]..high[1] {
            for x in low[0]..high[0] {
                if visited.is_multiple_of(TILE_CELLS) {
                    checkpoint(gate)?;
                }
                visited += 1;
                let cell = [x, y, z];
                let mut covered = 0usize;
                for sz in 0..samples {
                    for sy in 0..samples {
                        for sx in 0..samples {
                            let offset = [sx, sy, sz];
                            let relative = core::array::from_fn(|i| {
                                cell[i] as f64 + (offset[i] as f64 + 0.5) / samples as f64
                                    - center[i]
                            });
                            let local = inverse.apply(vector(relative));
                            if !array(local).into_iter().all(f64::is_finite) {
                                return Err(PartialSaturationError3::Cell {
                                    cell,
                                    reason: "body-frame sample overflowed",
                                });
                            }
                            covered += usize::from(shape.contains(local));
                        }
                    }
                }
                let epsilon = covered as f64 / denominator;
                if let Some((linear, angular)) = velocity {
                    let r = core::array::from_fn(|i| cell[i] as f64 + 0.5 - center[i]);
                    if covered > 0 && !low_mach(point_velocity(linear, angular, r)) {
                        return Err(PartialSaturationError3::Cell {
                            cell,
                            reason: "moving body sampled point exceeds low-Mach envelope",
                        });
                    }
                }
                if let Some(field) = fractions.as_mut() {
                    field[(z * dims[1] + y) * dims[0] + x] = epsilon;
                }
                volume += epsilon;
            }
        }
    }
    checkpoint(gate)?;
    Ok(volume)
}
