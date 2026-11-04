//! Atomic flow preparation and held-force changes in the original Guo frame.

use super::*;

mod frame;
pub use frame::RigidFrame3;
mod impulse;
mod resistance;
pub use resistance::{LinearResistance3, ResistanceMomentumExchange3};
pub(super) use resistance::{PreparedResistance3, ResistanceStream3};

// Own a sampled endpoint until the shared publication point. Existing scalar
// histories retain their selected shape and original arithmetic.
pub(super) enum EndpointForcing3 {
    Scale(f64),
    Field {
        uniform: [f64; 3],
        local: Vec<[f64; 3]>,
    },
}

impl EndpointForcing3 {
    pub(super) fn sampled_field(&self) -> Option<([f64; 3], &[[f64; 3]])> {
        match self {
            Self::Field { uniform, local } => Some((*uniform, local)),
            Self::Scale(_) => None,
        }
    }
}

/// Refusal before a flow or force-frame proposal is published.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BoundaryForcingError3 {
    /// Cancellation at a lattice-tile boundary or before publication.
    Cancelled,
    /// Invalid physical data or unsupported boundary/collision combination.
    Invalid(&'static str),
    /// Nonpositive density or nonfinite force, velocity or proposed populations.
    Unphysical {
        /// Original cell coordinates.
        cell: [usize; 3],
    },
}
impl core::fmt::Display for BoundaryForcingError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "D3Q19 flow preparation: {self:?}")
    }
}
impl std::error::Error for BoundaryForcingError3 {}

impl BoundaryGrid3 {
    /// Advance one lattice tick with a prescribed endpoint force history.
    ///
    /// The held force belongs to `t_n`; `next_scale` multiplies the field's
    /// selected spatial shape at `t_(n+1)`. Collision and streaming use the
    /// OLD force, then accepted populations and the NEW force are published
    /// together. Since physical momentum is `sum(f_q c_q) + F/2`, the net
    /// body-force impulse is `(F_n + F_(n+1))/2`, including the unchanged
    /// uniform force. In a periodic domain, internal transport cancels in
    /// this balance; wall reactions and open-boundary flux remain separate.
    ///
    /// Initialize `F_0` with [`Self::rebind_force_field`] to retain the chosen
    /// initial physical velocity. This is trapezoidal integration of a smooth,
    /// prescribed fixed-shape force history, not an implicit solve for a force
    /// depending on the unknown endpoint velocity. A discontinuous held load
    /// instead uses `rebind_force_field` followed by an ordinary step.
    /// See Guo et al., PRE 65, 046308 (2002), doi:10.1103/PhysRevE.65.046308.
    ///
    /// Uses the original collision/stream functions and reusable publication
    /// scratch. Invalid endpoint data, collision refusal or cancellation leave
    /// both the accepted populations and force unchanged. Checkpoints occur
    /// at lattice-tile boundaries and immediately before publication.
    pub fn step_with_endpoint_force_scale(
        &mut self,
        next_scale: f64,
        gate: &CancelGate,
    ) -> Result<(), BoundaryStepError3> {
        self.validate_endpoint_force_scale(next_scale, gate)?;
        self.step_with_endpoint_forcing_inner(EndpointForcing3::Scale(next_scale), None, gate)
            .map(|_| ())
    }

    /// Advance with a prescribed uniform and spatial force at the next time.
    ///
    /// Both inputs are lattice force densities, not accelerations. `local`
    /// samples fluid cells in x-fastest order and is added to `uniform`. The
    /// force may change direction and spatial shape independently at every
    /// step. Static solid cells are skipped. Collision uses the accepted old
    /// force; open faces reconstruct physical velocity, density and stress in
    /// the NEW endpoint frame before populations and forces publish together.
    /// Thus body forcing supplies `(F_n + F_(n+1))/2` per lattice tick, with
    /// wall reaction and open-boundary transport accounted separately.
    ///
    /// The sampled field becomes the shape for later scalar histories. Use
    /// `rebind_uniform_force` and `rebind_force_field_with_gate` for an
    /// instantaneous change of a held load that supplies no elapsed impulse.
    /// This prescribed-source step does not solve velocity-dependent forcing.
    ///
    /// One sampled force vector is allocated; native population scratch is
    /// reused. Cancellation, invalid data or a late reconstruction refusal
    /// publishes neither populations nor either force component.
    pub fn step_with_endpoint_forcing(
        &mut self,
        uniform: [f64; 3],
        local: impl FnMut([usize; 3]) -> [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryStepError3> {
        let endpoint = self.sample_endpoint_forcing(uniform, local, gate)?;
        self.step_with_endpoint_forcing_inner(endpoint, None, gate)
            .map(|_| ())
    }

    /// Endpoint-force step with the native impulse, torque and work on selected
    /// walls. Selection and moving-mass conventions match
    /// [`Self::step_with_wall_momentum`]; the loads come from the same collision
    /// used by the accepted stream, without a second collision or grid clone.
    pub fn step_with_endpoint_forcing_and_wall_momentum(
        &mut self,
        uniform: [f64; 3],
        local: impl FnMut([usize; 3]) -> [f64; 3],
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        self.validate_wall_measurement(links, moment_origin)?;
        let endpoint = self.sample_endpoint_forcing(uniform, local, gate)?;
        self.step_with_endpoint_forcing_inner(endpoint, Some((links, moment_origin)), gate)
    }

    fn step_with_endpoint_forcing_inner(
        &mut self,
        endpoint: EndpointForcing3,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
        gate: &CancelGate,
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        let collide = collision_kernel3();
        for tile in 0..self.f[0].len() {
            endpoint_checkpoint(gate)?;
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
            }
        }
        let exchange = match measurement {
            Some((links, origin)) => self.wall_momentum(links, origin, Some(gate))?,
            None => WallMomentumExchange3::default(),
        };
        let mut next = self.next.take().unwrap_or_else(|| self.f.clone());
        let view = StreamView {
            post: &self.post,
            pull: &self.pull,
            link_masks: &self.link_masks,
            open_link_masks: &self.open_link_masks,
            stationary_link_masks: &self.stationary_link_masks,
            boundaries: self.boundaries,
            dims: [self.nx, self.ny, self.nz],
        };
        let endpoint_field = endpoint.sampled_field();
        for tile in 0..self.f[0].len() {
            if let Err(error) = endpoint_checkpoint(gate) {
                self.next = Some(next);
                return Err(error);
            }
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let values = view.pull_cell(tile, lane);
                if !view.endpoint_is_physical(&values, tile, lane, endpoint_field) {
                    self.next = Some(next);
                    return Err(BoundaryStepError3::Unphysical {
                        cell: self.coords(tile, lane),
                    });
                }
                for (field, value) in next.iter_mut().zip(values) {
                    field[tile].0[lane] = value;
                }
            }
        }
        self.publish_next_with_forcing(next, gate, true, Some(endpoint))?;
        self.topology_locked = true;
        Ok(exchange)
    }

    /// Pooled counterpart of [`Self::step_with_endpoint_force_scale`]. Both
    /// methods call the same per-cell kernels, so worker count changes neither
    /// accepted population bits nor the endpoint force. The existing atomic
    /// pooled publication path owns collision, streaming and cancellation.
    pub fn step_pooled_with_endpoint_force_scale<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        next_scale: f64,
    ) -> Result<(), BoundaryStepError3> {
        self.validate_endpoint_force_scale(next_scale, gate)?;
        // An empty selection uses the existing final cancellation checkpoint
        // without scanning wall links or introducing another stepping kernel.
        self.step_pooled_coupling_inner(
            pool,
            gate,
            Some((&[], [0.0; 3])),
            None,
            None,
            Some(EndpointForcing3::Scale(next_scale)),
        )
        .map(|_| ())
    }

    /// Pooled counterpart of [`Self::step_with_endpoint_forcing`], with the
    /// same sampled field, native cell kernels and atomic publication.
    pub fn step_pooled_with_endpoint_forcing<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        uniform: [f64; 3],
        local: impl FnMut([usize; 3]) -> [f64; 3],
    ) -> Result<(), BoundaryStepError3> {
        let endpoint = self.sample_endpoint_forcing(uniform, local, gate)?;
        self.step_pooled_coupling_inner(
            pool,
            gate,
            Some((&[], [0.0; 3])),
            None,
            None,
            Some(endpoint),
        )
        .map(|_| ())
    }

    /// Pooled endpoint-force step returning the same selected-link loads as
    /// [`Self::step_with_endpoint_forcing_and_wall_momentum`].
    pub fn step_pooled_with_endpoint_forcing_and_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        uniform: [f64; 3],
        local: impl FnMut([usize; 3]) -> [f64; 3],
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        self.validate_wall_measurement(links, moment_origin)?;
        let endpoint = self.sample_endpoint_forcing(uniform, local, gate)?;
        self.step_pooled_coupling_inner(
            pool,
            gate,
            Some((links, moment_origin)),
            None,
            None,
            Some(endpoint),
        )
        .map(|(exchange, _)| exchange)
    }

    fn sample_endpoint_forcing(
        &self,
        uniform: [f64; 3],
        local: impl FnMut([usize; 3]) -> [f64; 3],
        gate: &CancelGate,
    ) -> Result<EndpointForcing3, BoundaryStepError3> {
        self.sample_force_field_for_uniform(uniform, local, gate)
            .map(|local| EndpointForcing3::Field { uniform, local })
            .map_err(|error| match error {
                BoundaryForcingError3::Cancelled => BoundaryStepError3::Cancelled,
                BoundaryForcingError3::Invalid(reason) => {
                    BoundaryStepError3::EndpointForce { reason }
                }
                BoundaryForcingError3::Unphysical { cell } => {
                    BoundaryStepError3::Unphysical { cell }
                }
            })
    }

    pub(super) fn validate_endpoint_force_scale(
        &self,
        scale: f64,
        gate: &CancelGate,
    ) -> Result<(), BoundaryStepError3> {
        endpoint_checkpoint(gate)?;
        let refuse = |reason| BoundaryStepError3::EndpointForce { reason };
        if !scale.is_finite() {
            return Err(refuse("endpoint force scale must be finite"));
        }
        let shape = self
            .force_field_shape
            .as_deref()
            .or(self.force_field.as_deref())
            .ok_or_else(|| refuse("select a force field before advancing its history"))?;
        for tile in shape.chunks(TILE_CELLS) {
            endpoint_checkpoint(gate)?;
            if tile.iter().any(|local| {
                local
                    .iter()
                    .zip(self.force)
                    .any(|(value, uniform)| !(scale * value + uniform).is_finite())
            }) {
                return Err(refuse("scaled endpoint total force must be finite"));
            }
        }
        Ok(())
    }

    // The complete shape was admitted before stepping. Do not poll or refuse
    // between population publication and its matching endpoint force frame.
    pub(super) fn publish_endpoint_force_scale(&mut self, scale: f64) {
        let field = self.force_field.as_mut().expect("admitted force field");
        let shape = self.force_field_shape.get_or_insert_with(|| field.clone());
        for (held, local) in field.iter_mut().zip(shape.iter()) {
            *held = local.map(|value| scale * value);
        }
    }

    /// Held uniform volume-force density, excluding any per-cell contribution.
    pub const fn uniform_force_density(&self) -> [f64; 3] {
        self.force
    }

    /// Replace a spatial held volume-force density without an instantaneous impulse.
    ///
    /// Samples fluid cells in x-fastest order, skipping static solids. The field
    /// is ADDED to the unchanged uniform force. Partially occupied immersed-body
    /// cells remain fluid cells; their existing coupling owns any fictitious
    /// interior-force subtraction. A zero field disables the spatial contribution.
    ///
    /// The original Guo source transform preserves each cell's physical density,
    /// momentum and full force-corrected stress. No time advances and no work is
    /// supplied. The following collision applies the newly held load. This also
    /// selects the new spatial shape used by `set_force_field_scale`.
    ///
    /// One force-field allocation holds the sampled proposal; existing collision
    /// scratch holds populations. Cancellation or any late refusal publishes
    /// neither. An unchanged cell preserves its population bits. Callbacks may
    /// still panic, but even then accepted fields and populations remain intact.
    pub fn rebind_force_field_with_gate(
        &mut self,
        force: impl FnMut([usize; 3]) -> [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        let field = self.sample_force_field_with_gate(force, gate)?;
        for tile in 0..self.f[0].len() {
            poll(gate)?;
            for lane in 0..TILE_CELLS {
                let mut populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                if self.solid[tile] & (1u64 << lane) == 0 {
                    let previous = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                    let next = cell_force(self.force, Some(&field), tile, lane);
                    if previous != next {
                        let cell = self.coords(tile, lane);
                        admit_populations(&populations, cell)?;
                        let (rho, momentum) = self.raw_moments(cell[0], cell[1], cell[2]);
                        let velocity = core::array::from_fn(|axis| {
                            (momentum[axis] + 0.5 * previous[axis]) / rho
                        });
                        if !velocity.iter().all(|value| value.is_finite()) {
                            return Err(BoundaryForcingError3::Unphysical { cell });
                        }
                        crate::thermal::forcing::reframe3(
                            &mut populations,
                            velocity,
                            previous,
                            next,
                        );
                        admit_populations(&populations, cell)?;
                    }
                }
                for (proposed, value) in self.post.iter_mut().zip(populations) {
                    proposed[tile].0[lane] = value;
                }
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.f, &mut self.post);
        self.force_field = Some(field);
        self.force_field_shape = None;
        Ok(())
    }

    pub(super) fn sample_force_field_with_gate(
        &self,
        force: impl FnMut([usize; 3]) -> [f64; 3],
        gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, BoundaryForcingError3> {
        self.sample_force_field_for_uniform(self.force, force, gate)
    }

    fn sample_force_field_for_uniform(
        &self,
        uniform: [f64; 3],
        mut force: impl FnMut([usize; 3]) -> [f64; 3],
        gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, BoundaryForcingError3> {
        poll(gate)?;
        if !uniform.iter().all(|value| value.is_finite()) {
            return Err(BoundaryForcingError3::Invalid(
                "uniform force must be finite",
            ));
        }
        if !self.collision_model.supports_body_force() {
            return Err(BoundaryForcingError3::Invalid(
                "the selected collision model does not admit body forcing",
            ));
        }
        let mut field = vec![[0.0; 3]; self.f[0].len() * TILE_CELLS];
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    if ((z * self.ny + y) * self.nx + x).is_multiple_of(TILE_CELLS) {
                        poll(gate)?;
                    }
                    if self.is_solid(x, y, z) {
                        continue;
                    }
                    let cell = [x, y, z];
                    let value = force(cell);
                    if !value
                        .iter()
                        .zip(uniform)
                        .all(|(local, uniform)| local.is_finite() && (local + uniform).is_finite())
                    {
                        return Err(BoundaryForcingError3::Unphysical { cell });
                    }
                    let (tile, lane) = self.addr(x, y, z);
                    field[tile * TILE_CELLS + lane] = value;
                }
            }
        }
        poll(gate)?;
        Ok(field)
    }

    /// Replace the uniform part of a held volume-force density, in lattice units.
    ///
    /// The original Guo source transform preserves physical mass, momentum and
    /// force-corrected stress at each fluid cell. This changes the force held for
    /// subsequent steps; it advances no time and supplies no impulse or work.
    /// Existing per-cell forces and their selected spatial shape remain intact.
    /// A discontinuous pressure-gradient change therefore creates no artificial
    /// half-force impulse. `set_force_field_scale` preserves the same physical
    /// state when selecting a held force from a prescribed spatial shape;
    /// smooth endpoint histories instead use `step_with_endpoint_force_scale`,
    /// which owns the old-collision/new-endpoint trapezoidal timing.
    ///
    /// Uses existing collision scratch with no new population or force-field
    /// allocation. Cancellation or a late-cell refusal leaves accepted fields
    /// and forcing unchanged. An identical force preserves population bits.
    pub fn rebind_uniform_force(
        &mut self,
        force: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        poll(gate)?;
        if !force.iter().all(|value| value.is_finite()) {
            return Err(BoundaryForcingError3::Invalid(
                "uniform force must be finite",
            ));
        }
        if force != [0.0; 3] && !self.collision_model.supports_body_force() {
            return Err(BoundaryForcingError3::Invalid(
                "nonzero uniform force requires a forcing-compatible collision model",
            ));
        }
        if self.force == force {
            return Ok(());
        }
        for tile in 0..self.f[0].len() {
            poll(gate)?;
            for lane in 0..TILE_CELLS {
                let mut populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                if self.solid[tile] & (1u64 << lane) == 0 {
                    let previous = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                    let next = cell_force(force, self.force_field.as_deref(), tile, lane);
                    let cell = self.coords(tile, lane);
                    if !previous.iter().chain(&next).all(|value| value.is_finite()) {
                        return Err(BoundaryForcingError3::Unphysical { cell });
                    }
                    if previous != next {
                        // The original moment observer asserts positive density;
                        // turn a late invalid cell into this atomic typed refusal.
                        admit_populations(&populations, cell)?;
                        let (rho, momentum) = self.raw_moments(cell[0], cell[1], cell[2]);
                        let velocity = core::array::from_fn(|axis| {
                            (momentum[axis] + 0.5 * previous[axis]) / rho
                        });
                        if !velocity.iter().all(|value| value.is_finite()) {
                            return Err(BoundaryForcingError3::Unphysical { cell });
                        }
                        crate::thermal::forcing::reframe3(
                            &mut populations,
                            velocity,
                            previous,
                            next,
                        );
                        admit_populations(&populations, cell)?;
                    }
                }
                for (proposed, value) in self.post.iter_mut().zip(populations) {
                    proposed[tile].0[lane] = value;
                }
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.f, &mut self.post);
        self.force = force;
        Ok(())
    }

    /// Prepare uniform PHYSICAL density/velocity using the original equilibrium.
    ///
    /// Existing uniform and per-cell forces select each cell's Guo half-source
    /// frame. Thus the requested velocity is the force-corrected physical
    /// velocity even when a force is already held. Solid cells are unchanged.
    ///
    /// This deliberately replaces the fluid initial state and its stress; it
    /// is not a momentum-conserving operation on an evolved flow. No body state,
    /// boundary topology, force, viscosity or clock changes. Material adapters
    /// restrict this initializer to their initial accepted time. Wall mismatch
    /// describes a startup transient, not an already solved steady flow.
    /// The caller supplies physical Mach/regime admission for its scenario.
    pub fn initialize_uniform_equilibrium(
        &mut self,
        density: f64,
        velocity: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryForcingError3> {
        poll(gate)?;
        if !density.is_finite() || density <= 0.0 || !velocity.iter().all(|value| value.is_finite())
        {
            return Err(BoundaryForcingError3::Invalid(
                "initial density must be positive and velocity finite",
            ));
        }
        let equilibrium = equilibrium3(density, velocity);
        if !equilibrium.iter().all(|value| value.is_finite()) {
            return Err(BoundaryForcingError3::Invalid(
                "unrepresentable initial equilibrium",
            ));
        }
        for tile in 0..self.f[0].len() {
            poll(gate)?;
            for lane in 0..TILE_CELLS {
                let mut populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                if self.solid[tile] & (1u64 << lane) == 0 {
                    populations = equilibrium;
                    let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                    if force != [0.0; 3] {
                        crate::thermal::forcing::reframe3(
                            &mut populations,
                            velocity,
                            [0.0; 3],
                            force,
                        );
                    }
                    admit_populations(&populations, self.coords(tile, lane))?;
                }
                for (proposed, value) in self.post.iter_mut().zip(populations) {
                    proposed[tile].0[lane] = value;
                }
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.f, &mut self.post);
        Ok(())
    }
}

fn endpoint_checkpoint(gate: &CancelGate) -> Result<(), BoundaryStepError3> {
    if gate.is_requested() {
        Err(BoundaryStepError3::Cancelled)
    } else {
        Ok(())
    }
}

fn admit_populations(
    populations: &[f64; Q3],
    cell: [usize; 3],
) -> Result<(), BoundaryForcingError3> {
    let density = populations.iter().sum::<f64>();
    if populations.iter().all(|value| value.is_finite()) && density.is_finite() && density > 0.0 {
        Ok(())
    } else {
        Err(BoundaryForcingError3::Unphysical { cell })
    }
}
fn poll(gate: &CancelGate) -> Result<(), BoundaryForcingError3> {
    if gate.is_requested() {
        Err(BoundaryForcingError3::Cancelled)
    } else {
        Ok(())
    }
}
