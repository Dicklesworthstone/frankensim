//! Finite-inertia planar rotor coupled to the native D3Q19 wall links.
//!
//! The fluid stays in an inertial frame. The circular wall patch has a fixed
//! center and axis; only its spin responds to the fluid. Its bearing carries
//! translation and transverse torque. No rotating-frame force or fluid-inertia
//! correction is needed for this fixed geometry.

use fs_exec::CancelGate;

use super::{
    BoundaryGrid3, BoundaryLink3, BoundaryStepError3, E3, Face3, FaceBoundary3, OPP3, Q3,
    StreamView, TILE_CELLS, WallMomentumExchange3, cell_force, collision_kernel3,
    planar_face_cells, wall_population,
};

/// Refusal of a finite-inertia rotor proposal; accepted fluid and rotor state
/// remain unchanged.
#[derive(Debug)]
pub enum WallRotorError3 {
    /// Geometry, inertia, applied impulse, or solved motion is inadmissible.
    Invalid(&'static str),
    /// The original fluid collision, wall measurement, or cancellation refused.
    Boundary(BoundaryStepError3),
}

impl core::fmt::Display for WallRotorError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "D3Q19 wall rotor refused: {reason}"),
            Self::Boundary(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for WallRotorError3 {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Boundary(error) => Some(error),
            Self::Invalid(_) => None,
        }
    }
}

impl From<BoundaryStepError3> for WallRotorError3 {
    fn from(error: BoundaryStepError3) -> Self {
        Self::Boundary(error)
    }
}

#[derive(Debug, Clone, Copy)]
struct RotorCell3 {
    tile: usize,
    lane: usize,
    cell: [usize; 3],
    /// Velocity at unit spin: positive Cartesian axis crossed with radius.
    unit_velocity: [f64; 3],
}

/// A circular, finite-inertia rotating wall patch on a closed [BoundaryGrid3].
///
/// Coordinates, inertia, and impulses use lattice units (dx = dt = 1). Cell
/// centers are at (x + 0.5, y + 0.5, z + 0.5); wall planes are at 0 and the
/// corresponding grid dimension. Positive spin is about the positive Cartesian
/// axis normal to face, on both the minimum and maximum face.
///
/// The disk is discretized by selecting face-adjacent cells whose tangential
/// centers lie within radius. Every selected cell must have exactly the five
/// D3Q19 links into that face, with no side-wall seam or adjacent solid link.
/// Whole cells own their complete link set: an arbitrary cropped link selection
/// would inject mass and invalidate the scalar inertia response.
///
/// The remaining exterior and voxel walls are stationary. Open and periodic
/// domains are refused. This owner models a flush disk and its bearing, not a
/// vessel rotating around fluid or a translating/tilting immersed body.
#[derive(Debug, Clone)]
pub struct WallRotor3 {
    dimensions: [usize; 3],
    face: Face3,
    center: [f64; 3],
    radius: f64,
    inertia: f64,
    angular_velocity: f64,
    angle: f64,
    face_mask: u32,
    cells: Vec<RotorCell3>,
}

impl WallRotor3 {
    /// Prepare a reusable rotor patch without changing the grid.
    ///
    /// maximum_face_cells bounds the complete planar face traversal and the
    /// prepared allocation. The patch must fit the face and contain a resolved
    /// off-axis cell. Spin must keep every selected wall speed within the native
    /// low-Mach envelope (squared wall speed < 0.03). Patch topology and all face
    /// rules are revalidated on every coupled step.
    ///
    /// # Errors
    /// Invalid geometry/inertia/spin, unsupported face rules, incomplete patch
    /// cells, exceeded traversal budget, or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        grid: &BoundaryGrid3,
        face: Face3,
        center: [f64; 3],
        radius: f64,
        inertia: f64,
        angular_velocity: f64,
        maximum_face_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, WallRotorError3> {
        checkpoint(gate)?;
        stationary_closed(grid)?;
        let dimensions = [grid.nx, grid.ny, grid.nz];
        let axis = face.axis();
        let plane = if face.is_min() {
            0.0
        } else {
            dimensions[axis] as f64
        };
        let radius_sq = radius * radius;
        if !center.into_iter().all(f64::is_finite)
            || center[axis] != plane
            || !(radius.is_finite() && radius > 0.0 && radius_sq.is_finite())
            || !(inertia.is_finite() && inertia > 0.0)
            || !angular_velocity.is_finite()
        {
            return Err(WallRotorError3::Invalid(
                "invalid center, radius, inertia, or spin",
            ));
        }
        for tangent in (0..3).filter(|&index| index != axis) {
            if center[tangent] - radius < 0.0
                || center[tangent] + radius > dimensions[tangent] as f64
            {
                return Err(WallRotorError3::Invalid(
                    "rotor disk must fit inside its wall face",
                ));
            }
        }
        let face_cells = dimensions[(axis + 1) % 3]
            .checked_mul(dimensions[(axis + 2) % 3])
            .ok_or(WallRotorError3::Invalid("face cell count overflow"))?;
        if face_cells > maximum_face_cells {
            return Err(WallRotorError3::Invalid(
                "face exceeds the preparation budget",
            ));
        }
        let incoming_sign = if face.is_min() { 1 } else { -1 };
        let face_mask = E3.iter().enumerate().fold(0u32, |mask, (q, e)| {
            if [e.0, e.1, e.2][axis] == incoming_sign {
                mask | (1u32 << q)
            } else {
                mask
            }
        });
        let mut cells = Vec::new();
        for (index, cell) in planar_face_cells(dimensions, face).enumerate() {
            if index.is_multiple_of(TILE_CELLS) {
                checkpoint(gate)?;
            }
            let mut offset = core::array::from_fn(|i| cell[i] as f64 + 0.5 - center[i]);
            offset[axis] = 0.0;
            let r_sq = dot(offset, offset);
            if r_sq > radius_sq {
                continue;
            }
            let (tile, lane) = grid.addr(cell[0], cell[1], cell[2]);
            let mut axis_vector = [0.0; 3];
            axis_vector[axis] = 1.0;
            let candidate = RotorCell3 {
                tile,
                lane,
                cell,
                unit_velocity: cross(axis_vector, offset),
            };
            validate_cell(grid, candidate, face_mask)?;
            validate_speed(candidate.unit_velocity, angular_velocity)?;
            cells.push(candidate);
        }
        if !cells
            .iter()
            .any(|cell| dot(cell.unit_velocity, cell.unit_velocity) > 0.0)
        {
            return Err(WallRotorError3::Invalid(
                "rotor patch needs a resolved off-axis cell",
            ));
        }
        // The stream visits tile/lane order. A single cursor then recognizes
        // the sparse patch without a dense mask or a per-cell search.
        cells.sort_unstable_by_key(|cell| (cell.tile, cell.lane));
        checkpoint(gate)?;
        Ok(Self {
            dimensions,
            face,
            center,
            radius,
            inertia,
            angular_velocity,
            angle: 0.0,
            face_mask,
            cells,
        })
    }

    /// Face containing the disk.
    #[must_use]
    pub const fn face(&self) -> Face3 {
        self.face
    }

    /// Fixed center on the physical wall plane.
    #[must_use]
    pub const fn center(&self) -> [f64; 3] {
        self.center
    }

    /// Radius used to select complete face cells.
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// Axial moment of inertia in lattice mass-length-squared units.
    #[must_use]
    pub const fn inertia(&self) -> f64 {
        self.inertia
    }

    /// Accepted angular velocity about the positive Cartesian face axis.
    #[must_use]
    pub const fn angular_velocity(&self) -> f64 {
        self.angular_velocity
    }

    /// Accumulated rotation in radians, advanced with midpoint angular velocity.
    /// The circular patch itself is independent of this phase.
    #[must_use]
    pub const fn angle(&self) -> f64 {
        self.angle
    }

    /// Number of complete face cells in the prepared disk.
    #[must_use]
    pub fn patch_cell_count(&self) -> usize {
        self.cells.len()
    }

    fn validate(&self, grid: &BoundaryGrid3, gate: &CancelGate) -> Result<(), WallRotorError3> {
        checkpoint(gate)?;
        if self.dimensions != [grid.nx, grid.ny, grid.nz] {
            return Err(WallRotorError3::Invalid("rotor and grid dimensions differ"));
        }
        stationary_closed(grid)?;
        for (index, &cell) in self.cells.iter().enumerate() {
            if index.is_multiple_of(TILE_CELLS) {
                checkpoint(gate)?;
            }
            validate_cell(grid, cell, self.face_mask)?;
        }
        Ok(())
    }
}

/// Accepted fluid, rotor, and housing exchange for one native lattice step.
///
/// Wall exchanges and bearing angular impulse use the requested moment origin.
/// Bearing impulses act **on the rotor**; the housing receives their negatives.
/// The rotor's only momentum increment is its axial inertia times spin change.
/// Consequently fluid + rotor + housing momentum closes against the applied
/// rotor couple and held fluid body force, including about an arbitrary origin.
///
/// Mechanical work is the rotor's midpoint impulse work. Viscous and numerical
/// fluid dissipation remain with the original LBM owner; this receipt does not
/// assert conservation of fluid kinetic energy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallRotorExchange3 {
    /// Accepted spin before this step.
    pub angular_velocity_before: f64,
    /// Held wall spin during the native bounce.
    pub angular_velocity_midpoint: f64,
    /// Accepted spin after the fluid reaction.
    pub angular_velocity_after: f64,
    /// Native moving-wall exchange on the complete rotating patch.
    pub hydrodynamic: WallMomentumExchange3,
    /// Native exchange on every other exterior or voxel wall.
    pub stationary_walls: WallMomentumExchange3,
    /// Fluid angular impulse on the rotor about its center, projected on its axis.
    pub hydrodynamic_angular_impulse: f64,
    /// Bearing impulse on the rotor, cancelling its hydrodynamic force.
    pub bearing_impulse: [f64; 3],
    /// Bearing angular impulse on the rotor about the requested moment origin.
    /// It cancels transverse motion and the moment of the constrained force.
    pub bearing_angular_impulse: [f64; 3],
    /// Held fluid body-force impulse summed over the actual collided fluid cells.
    pub body_force_impulse: [f64; 3],
    /// Moment of that body force at native cell centers about the requested origin.
    pub body_force_angular_impulse: [f64; 3],
    /// Applied axial rotor couple integrated over this lattice step.
    pub external_angular_impulse: f64,
    /// Rotor kinetic energy change, evaluated in its stable polarized form.
    pub kinetic_energy_change: f64,
    /// Applied couple work: external angular impulse times midpoint spin.
    pub external_work: f64,
    /// Kinetic energy change minus measured fluid-on-rotor work and applied work.
    pub mechanical_work_residual: f64,
}

/// One rotor's accepted response in a simultaneous multi-rotor step.
///
/// All impulses act on this rotor. Angular impulses use the requested moment
/// origin, except the scalar hydrodynamic angular impulse, which uses the
/// rotor center and its positive Cartesian axis. Global housing/body-force
/// exchanges appear only once in [WallRotorsExchange3].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallRotorResponse3 {
    /// Accepted spin before this step.
    pub angular_velocity_before: f64,
    /// Held wall spin during the native bounce.
    pub angular_velocity_midpoint: f64,
    /// Accepted spin after the fluid reaction.
    pub angular_velocity_after: f64,
    /// Native moving-wall exchange on this rotor's complete patch.
    pub hydrodynamic: WallMomentumExchange3,
    /// Fluid angular impulse about this rotor's center, projected on its axis.
    pub hydrodynamic_angular_impulse: f64,
    /// Bearing impulse on this rotor, cancelling its hydrodynamic force.
    pub bearing_impulse: [f64; 3],
    /// Bearing angular impulse on this rotor about the requested moment origin.
    pub bearing_angular_impulse: [f64; 3],
    /// Applied axial couple integrated over this lattice step.
    pub external_angular_impulse: f64,
    /// Rotor kinetic energy change in stable polarized form.
    pub kinetic_energy_change: f64,
    /// Applied axial couple times midpoint spin.
    pub external_work: f64,
    /// Kinetic change minus measured fluid-on-rotor work and applied work.
    pub mechanical_work_residual: f64,
}

/// Simultaneous native fluid coupling to a fixed array of independent rotors.
///
/// Each rotor response appears in input-array order. The stationary housing
/// and held fluid body-force exchanges are global and are counted exactly once.
/// Summing fluid and rotor momentum changes plus housing impulses minus all
/// bearing impulses closes against the fluid body force and applied couples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallRotorsExchange3<const N: usize> {
    /// Individual rotor responses in the same order as the input array.
    pub rotors: [WallRotorResponse3; N],
    /// Exchange on all walls outside the union of rotating patches.
    pub stationary_walls: WallMomentumExchange3,
    /// Held body-force impulse over all collided fluid cells, counted once.
    pub body_force_impulse: [f64; 3],
    /// Moment of the held fluid body force about the requested origin.
    pub body_force_angular_impulse: [f64; 3],
}

#[derive(Clone, Copy, Default)]
struct RotorMotion3 {
    before: f64,
    midpoint: f64,
    after: f64,
    angle: f64,
}

impl BoundaryGrid3 {
    /// Advance native fluid collision/streaming and a finite-inertia wall rotor.
    ///
    /// The post-collision wall response is affine in the held spin:
    /// J_h = A - D * omega_mid, with
    /// D = sum_patch(rho * |axis cross radius|^2 / 3) >= 0.
    /// Solve the scalar midpoint inertia equation analytically, then use that
    /// same midpoint spin in the original moving-wall bounce and link exchange.
    /// The response has no explicit added-inertia stability restriction.
    ///
    /// One common velocity owns every link of a selected face cell. Tangential
    /// opposite links then cancel moving-wall mass injection and the Wen
    /// moving-mass correction. Its projection into each link is also the exact
    /// rigid rotation at the half-link intersection, since
    /// e dot (omega cross e) = 0.
    ///
    /// The external impulse is a pure axial couple. Bearings do no work because
    /// the center is fixed and their torque has no allowed axial component.
    /// All other walls are stationary; no periodic/open face is admitted.
    /// Use [Self::step_with_finite_inertia_rotors] for several simultaneous rotors.
    /// Call this method for each coupled step: ordinary [Self::step] has no
    /// rotor owner and treats the patch as the underlying stationary wall.
    ///
    /// Cancellation is checked at least once per 64-cell tile/pass and before
    /// the final population swap. Refusal leaves accepted populations, forces,
    /// topology lock, and rotor motion intact. Prepared cells and the grid's
    /// reusable next-state buffer avoid per-step dense allocation after warmup.

    pub fn step_with_finite_inertia_rotor(
        &mut self,
        rotor: &mut WallRotor3,
        external_angular_impulse: f64,
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<WallRotorExchange3, WallRotorError3> {
        let batch = self.step_with_finite_inertia_rotors(
            core::array::from_mut(rotor),
            [external_angular_impulse],
            moment_origin,
            gate,
        )?;
        let [response] = batch.rotors;
        Ok(WallRotorExchange3 {
            angular_velocity_before: response.angular_velocity_before,
            angular_velocity_midpoint: response.angular_velocity_midpoint,
            angular_velocity_after: response.angular_velocity_after,
            hydrodynamic: response.hydrodynamic,
            stationary_walls: batch.stationary_walls,
            hydrodynamic_angular_impulse: response.hydrodynamic_angular_impulse,
            bearing_impulse: response.bearing_impulse,
            bearing_angular_impulse: response.bearing_angular_impulse,
            body_force_impulse: batch.body_force_impulse,
            body_force_angular_impulse: batch.body_force_angular_impulse,
            external_angular_impulse: response.external_angular_impulse,
            kinetic_energy_change: response.kinetic_energy_change,
            external_work: response.external_work,
            mechanical_work_residual: response.mechanical_work_residual,
        })
    }

    /// Advance several independent finite-inertia wall rotors in one fluid step.
    ///
    /// All rotors see the same original post-collision fluid state. Their
    /// complete patches must be disjoint; repeated or intersecting patch cells
    /// refuse before collision. Each rotor's affine wall response is solved at
    /// its own midpoint spin, and one native stream publishes all rotor motion
    /// and fluid populations together. Calling the single-rotor step repeatedly
    /// would instead advance the fluid multiple times.
    ///
    /// Rotors can occupy the same face, opposite faces, or different Cartesian
    /// axes. The fixed-size arrays keep per-rotor workspaces bounded. A sparse
    /// ordered merge recognizes patch cells without a dense per-cell owner map:
    /// ordinary interior cells compare one cached address, and rotor cursors are
    /// searched only when a patch cell is consumed.
    ///
    /// Global housing and body-force exchanges appear exactly once in the
    /// returned receipt; each rotor retains its own bearing and work exchange.
    /// Geometry/low-Mach rules match [Self::step_with_finite_inertia_rotor].
    /// An empty array refuses. Cancellation or any late rotor/fluid refusal
    /// publishes none of the proposed fluid or rotor states.
    pub fn step_with_finite_inertia_rotors<const N: usize>(
        &mut self,
        rotors: &mut [WallRotor3; N],
        external_angular_impulses: [f64; N],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<WallRotorsExchange3<N>, WallRotorError3> {
        checkpoint(gate)?;
        if N == 0 {
            return Err(WallRotorError3::Invalid("at least one rotor is required"));
        }
        if !external_angular_impulses.into_iter().all(f64::is_finite)
            || !moment_origin.into_iter().all(f64::is_finite)
        {
            return Err(WallRotorError3::Invalid(
                "applied impulses and moment origin must be finite",
            ));
        }
        for rotor in rotors.iter() {
            rotor.validate(self, gate)?;
        }
        validate_disjoint(rotors, gate)?;

        let mut body_force_impulse = [0.0; 3];
        let mut body_force_angular_impulse = [0.0; 3];
        let collide = collision_kernel3();
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                let post = collide(populations, self.collision_model, force).map_err(|source| {
                    BoundaryStepError3::Collision {
                        cell: self.coords(tile, lane),
                        source,
                    }
                })?;
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
                let cell = self.coords(tile, lane);
                let offset = core::array::from_fn(|i| cell[i] as f64 + 0.5 - moment_origin[i]);
                let moment = cross(offset, force);
                for i in 0..3 {
                    body_force_impulse[i] += force[i];
                    body_force_angular_impulse[i] += moment[i];
                }
            }
        }
        let mut motions = [RotorMotion3::default(); N];
        for i in 0..N {
            motions[i] = propose_motion(self, &rotors[i], external_angular_impulses[i], gate)?;
        }

        let mut next = self.next.take().unwrap_or_else(|| self.f.clone());
        let outcome = (|| {
            let view = StreamView {
                post: &self.post,
                pull: &self.pull,
                link_masks: &self.link_masks,
                open_link_masks: &self.open_link_masks,
                stationary_link_masks: &self.stationary_link_masks,
                boundaries: self.boundaries,
                dims: [self.nx, self.ny, self.nz],
            };
            let mut hydrodynamic = [WallMomentumExchange3::default(); N];
            let mut stationary_walls = WallMomentumExchange3::default();
            let mut cursors = [0; N];
            let mut next_patch = next_rotor_cell(rotors, &cursors);
            for tile in 0..self.f[0].len() {
                checkpoint(gate)?;
                for lane in 0..TILE_CELLS {
                    if self.solid[tile] & (1u64 << lane) != 0 {
                        continue;
                    }
                    let patch =
                        next_patch.filter(|(_, cell)| cell.tile == tile && cell.lane == lane);
                    let velocity = patch.map_or([0.0; 3], |(owner, cell)| {
                        cell.unit_velocity.map(|u| u * motions[owner].midpoint)
                    });
                    let mut values = view.pull_cell(tile, lane);
                    let mask = self.link_masks[tile].0[lane];
                    if mask != 0 {
                        let cell = self.coords(tile, lane);
                        let rho = (0..Q3).map(|q| self.post[q][tile].0[lane]).sum::<f64>();
                        for (q, value) in values.iter_mut().enumerate() {
                            if mask & (1u32 << q) == 0 {
                                continue;
                            }
                            let outgoing = self.post[OPP3[q]][tile].0[lane];
                            let link = BoundaryLink3 {
                                tile,
                                lane,
                                cell,
                                direction: q,
                            };
                            if let Some((owner, _)) = patch {
                                *value = wall_population(outgoing, rho, q, velocity);
                                hydrodynamic[owner].add_link(
                                    link,
                                    outgoing,
                                    rho,
                                    velocity,
                                    moment_origin,
                                );
                            } else {
                                stationary_walls.add_link(
                                    link,
                                    outgoing,
                                    rho,
                                    [0.0; 3],
                                    moment_origin,
                                );
                            }
                        }
                    }
                    let rho_next = values.iter().sum::<f64>();
                    if !(values.iter().all(|value| value.is_finite())
                        && rho_next.is_finite()
                        && rho_next > 0.0)
                    {
                        return Err(BoundaryStepError3::Unphysical {
                            cell: self.coords(tile, lane),
                        }
                        .into());
                    }
                    for (field, value) in next.iter_mut().zip(values) {
                        field[tile].0[lane] = value;
                    }
                    if let Some((owner, _)) = patch {
                        cursors[owner] += 1;
                        next_patch = next_rotor_cell(rotors, &cursors);
                    }
                }
            }
            debug_assert!(next_patch.is_none());
            for exchange in &mut hydrodynamic {
                *exchange = exchange.checked()?;
            }
            let stationary_walls = stationary_walls.checked()?;
            let responses = core::array::from_fn(|i| {
                rotor_response(
                    &rotors[i],
                    motions[i],
                    hydrodynamic[i],
                    external_angular_impulses[i],
                    moment_origin,
                )
            });
            if ![body_force_impulse, body_force_angular_impulse]
                .into_iter()
                .flatten()
                .all(f64::is_finite)
                || !responses.iter().all(response_is_finite)
            {
                return Err(WallRotorError3::Invalid(
                    "rotor momentum or work is non-finite",
                ));
            }
            checkpoint(gate)?;
            Ok(WallRotorsExchange3 {
                rotors: responses,
                stationary_walls,
                body_force_impulse,
                body_force_angular_impulse,
            })
        })();
        match outcome {
            Ok(exchange) => {
                core::mem::swap(&mut self.f, &mut next);
                self.next = Some(next);
                self.topology_locked = true;
                for (rotor, motion) in rotors.iter_mut().zip(motions) {
                    rotor.angular_velocity = motion.after;
                    rotor.angle = motion.angle;
                }
                Ok(exchange)
            }
            Err(error) => {
                self.next = Some(next);
                Err(error)
            }
        }
    }
}

fn propose_motion(
    grid: &BoundaryGrid3,
    rotor: &WallRotor3,
    external_angular_impulse: f64,
    gate: &CancelGate,
) -> Result<RotorMotion3, WallRotorError3> {
    let mut stationary_angular_impulse = 0.0;
    let mut damping = 0.0;
    for (index, cell) in rotor.cells.iter().enumerate() {
        if index.is_multiple_of(TILE_CELLS) {
            checkpoint(gate)?;
        }
        let rho = (0..Q3)
            .map(|q| grid.post[q][cell.tile].0[cell.lane])
            .sum::<f64>();
        if !(rho.is_finite() && rho > 0.0) {
            return Err(BoundaryStepError3::Unphysical { cell: cell.cell }.into());
        }
        damping += rho * dot(cell.unit_velocity, cell.unit_velocity) / 3.0;
        for (q, e) in E3.iter().enumerate() {
            if rotor.face_mask & (1u32 << q) == 0 {
                continue;
            }
            let projection = dot(
                [f64::from(e.0), f64::from(e.1), f64::from(e.2)],
                cell.unit_velocity,
            );
            stationary_angular_impulse -=
                2.0 * projection * grid.post[OPP3[q]][cell.tile].0[cell.lane];
        }
    }
    let before = rotor.angular_velocity;
    let denominator = rotor.inertia + 0.5 * damping;
    let delta =
        (external_angular_impulse + stationary_angular_impulse - damping * before) / denominator;
    let midpoint = before + 0.5 * delta;
    let after = before + delta;
    let angle = rotor.angle + midpoint;
    if ![
        stationary_angular_impulse,
        damping,
        denominator,
        delta,
        midpoint,
        after,
        angle,
    ]
    .into_iter()
    .all(f64::is_finite)
        || denominator <= 0.0
    {
        return Err(WallRotorError3::Invalid(
            "rotor midpoint response is non-finite",
        ));
    }
    for (index, cell) in rotor.cells.iter().enumerate() {
        if index.is_multiple_of(TILE_CELLS) {
            checkpoint(gate)?;
        }
        validate_speed(cell.unit_velocity, midpoint)?;
        validate_speed(cell.unit_velocity, after)?;
    }
    Ok(RotorMotion3 {
        before,
        midpoint,
        after,
        angle,
    })
}

fn rotor_response(
    rotor: &WallRotor3,
    motion: RotorMotion3,
    hydrodynamic: WallMomentumExchange3,
    external_angular_impulse: f64,
    origin: [f64; 3],
) -> WallRotorResponse3 {
    let center_offset = core::array::from_fn(|i| rotor.center[i] - origin[i]);
    let force_moment = cross(center_offset, hydrodynamic.wall_impulse);
    let axis = rotor.face.axis();
    let hydrodynamic_angular_impulse = hydrodynamic.wall_angular_impulse[axis] - force_moment[axis];
    let bearing_impulse = hydrodynamic.wall_impulse.map(|impulse| -impulse);
    let bearing_angular_impulse = core::array::from_fn(|i| {
        if i == axis {
            hydrodynamic_angular_impulse - hydrodynamic.wall_angular_impulse[i]
        } else {
            -hydrodynamic.wall_angular_impulse[i]
        }
    });
    let kinetic_energy_change = rotor.inertia * (motion.after - motion.before) * motion.midpoint;
    let external_work = external_angular_impulse * motion.midpoint;
    let mechanical_work_residual = kinetic_energy_change - external_work - hydrodynamic.wall_work;
    WallRotorResponse3 {
        angular_velocity_before: motion.before,
        angular_velocity_midpoint: motion.midpoint,
        angular_velocity_after: motion.after,
        hydrodynamic,
        hydrodynamic_angular_impulse,
        bearing_impulse,
        bearing_angular_impulse,
        external_angular_impulse,
        kinetic_energy_change,
        external_work,
        mechanical_work_residual,
    }
}

fn response_is_finite(response: &WallRotorResponse3) -> bool {
    [response.bearing_impulse, response.bearing_angular_impulse]
        .into_iter()
        .flatten()
        .chain([
            response.hydrodynamic_angular_impulse,
            response.kinetic_energy_change,
            response.external_work,
            response.mechanical_work_residual,
        ])
        .all(f64::is_finite)
}

// Lists are prepared in native tile/lane order. Search only when a sparse
// patch entry is consumed; this helper is never called for an ordinary cell.
fn next_rotor_cell<const N: usize>(
    rotors: &[WallRotor3; N],
    cursors: &[usize; N],
) -> Option<(usize, RotorCell3)> {
    rotors
        .iter()
        .enumerate()
        .filter_map(|(owner, rotor)| {
            rotor
                .cells
                .get(cursors[owner])
                .copied()
                .map(|cell| (owner, cell))
        })
        .min_by_key(|(_, cell)| (cell.tile, cell.lane))
}

fn validate_disjoint<const N: usize>(
    rotors: &[WallRotor3; N],
    gate: &CancelGate,
) -> Result<(), WallRotorError3> {
    let mut cursors = [0; N];
    let mut previous = None;
    let mut visited = 0usize;
    while let Some((owner, cell)) = next_rotor_cell(rotors, &cursors) {
        if visited.is_multiple_of(TILE_CELLS) {
            checkpoint(gate)?;
        }
        let address = (cell.tile, cell.lane);
        if previous == Some(address) {
            return Err(WallRotorError3::Invalid("rotor patches overlap"));
        }
        previous = Some(address);
        cursors[owner] += 1;
        visited += 1;
    }
    checkpoint(gate)
}

fn stationary_closed(grid: &BoundaryGrid3) -> Result<(), WallRotorError3> {
    if Face3::ALL
        .into_iter()
        .any(|face| grid.boundaries.face(face) != FaceBoundary3::stationary_wall())
    {
        return Err(WallRotorError3::Invalid(
            "rotor requires a closed grid with stationary exterior walls",
        ));
    }
    Ok(())
}

fn validate_cell(grid: &BoundaryGrid3, cell: RotorCell3, mask: u32) -> Result<(), WallRotorError3> {
    if grid.solid[cell.tile] & (1u64 << cell.lane) != 0
        || grid.link_masks[cell.tile].0[cell.lane] != mask
        || grid.open_link_masks[cell.tile].0[cell.lane] != 0
        || grid.stationary_link_masks[cell.tile].0[cell.lane] != 0
    {
        return Err(WallRotorError3::Invalid(
            "rotor cell must own a complete unblocked planar face stencil",
        ));
    }
    Ok(())
}

fn validate_speed(unit_velocity: [f64; 3], spin: f64) -> Result<(), WallRotorError3> {
    let velocity = unit_velocity.map(|u| u * spin);
    let speed_sq = dot(velocity, velocity);
    if !speed_sq.is_finite() || speed_sq >= 0.03 {
        return Err(WallRotorError3::Invalid(
            "rotor wall speed exceeds the low-Mach admission envelope",
        ));
    }
    Ok(())
}

fn checkpoint(gate: &CancelGate) -> Result<(), WallRotorError3> {
    if gate.is_requested() {
        Err(BoundaryStepError3::Cancelled.into())
    } else {
        Ok(())
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
