//! Actual-density acceleration through the original open reservoir owner.

use super::*;
use fs_exec::KernelRunner;

impl<L: EquilibriumLattice3<Q>, const Q: usize> ForcedOpenBoundaryGrid3<L, Q> {
    /// Attach acceleration without changing accepted density, velocity or stress.
    /// Values use lattice acceleration (`a_SI*dt²/dx`) in x-fastest cell order;
    /// solid entries must be zero. The native source transform handles the
    /// selected D3Q19/D3Q27 collision family and preserves existing time/profiles.
    pub fn try_new_with_body_acceleration_preserving_velocity(
        inner: OpenBoundaryGrid3<L, Q>,
        cx: &Cx<'_>,
        mut acceleration: Vec<[f64; 3]>,
    ) -> Result<Self, GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        admit_forces(&inner.grid, &acceleration, None, &mut poll)?;
        for (cell, value) in acceleration.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll()?;
            }
            if !inner.grid.solid[cell] {
                let rho = inner.grid.f[cell].iter().sum::<f64>();
                *value = value.map(|a| rho * a);
            }
        }
        Self::try_new_preserving_velocity_with_poll(inner, acceleration, &mut poll)
    }

    /// Advance with a prescribed acceleration field at the NEXT accepted time.
    /// Collision uses retained `F_n`. Donors and the final reconstructed cells
    /// use their own `F_(n+1)=rho_(n+1)*a_(n+1)`, including wall/open seam density.
    /// The native source is linear in density, so seam mass closes directly;
    /// no fixed-point iteration, frozen reference density or extra kick occurs.
    ///
    /// The field is borrowed and never altered. A reusable force scratch vector
    /// permits native stream workers to resolve endpoints while neighboring
    /// donors still read the authored acceleration. Body impulse is trapezoidal;
    /// the returned reservoir transport and `last_wall_exchange` remain separate.
    pub fn advance_serial_with_body_acceleration(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        acceleration_at_end: &[[f64; 3]],
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let mut poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        poll()?;
        self.prepare_force_workspace()?;
        let (receipt, step) = self.inner.stage_serial_with_source(
            self.inner.plan,
            collisions,
            Some(&self.forces),
            None,
            Some((acceleration_at_end, &mut self.pending_forces)),
            origin,
            &mut poll,
        )?;
        self.publish_acceleration(step);
        Ok(receipt)
    }

    /// Same source timing and closure on the original two executor passes.
    /// Physical endpoint admission and force resolution stay inside each stream
    /// tile. No extra full-grid density sweep or population buffer is introduced.
    pub fn advance_pooled_with_body_acceleration<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
        acceleration_at_end: &[[f64; 3]],
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        if gate.is_requested() {
            return Err(GridError3::Cancelled);
        }
        self.prepare_force_workspace()?;
        let (receipt, step) = self.inner.stage_pooled_with_source(
            pool,
            gate,
            self.inner.plan,
            collisions,
            Some(&self.forces),
            None,
            Some((acceleration_at_end, &mut self.pending_forces)),
            origin,
        )?;
        self.publish_acceleration(step);
        Ok(receipt)
    }

    fn prepare_force_workspace(&mut self) -> Result<(), GridError3<L::Error>> {
        let cells = self.inner.grid.f.len();
        if self.pending_forces.len() != cells {
            self.pending_forces
                .try_reserve_exact(cells - self.pending_forces.len())
                .map_err(|_| GridError3::Allocation)?;
            self.pending_forces.resize(cells, [0.0; 3]);
        }
        Ok(())
    }

    fn publish_acceleration(&mut self, step: usize) {
        core::mem::swap(&mut self.inner.grid.f, &mut self.inner.grid.next);
        core::mem::swap(&mut self.forces, &mut self.pending_forces);
        self.inner.grid.completed_steps = step;
        self.inner.publish_wall_exchange();
    }

    /// Executor counterpart of `advance_serial` for authored force densities.
    /// Uses the same native reservoir/wall transaction as acceleration endpoints.
    pub fn advance_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
        mut forces_at_end: Vec<[f64; 3]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let (receipt, step) = self.inner.stage_pooled_with_source(
            pool,
            gate,
            self.inner.plan,
            collisions,
            Some(&self.forces),
            Some(&forces_at_end),
            None,
            origin,
        )?;
        self.publish_with_force_swap(step, &mut forces_at_end);
        Ok(receipt)
    }

    /// Hold the accepted spatial force density over one pooled native interval.
    /// A held acceleration instead requires an explicit next acceleration field.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<OpenBoundaryMomentumTransport3, GridError3<L::Error>> {
        let (receipt, step) = self.inner.stage_pooled_with_source(
            pool,
            gate,
            self.inner.plan,
            collisions,
            Some(&self.forces),
            Some(&self.forces),
            None,
            origin,
        )?;
        core::mem::swap(&mut self.inner.grid.f, &mut self.inner.grid.next);
        self.inner.grid.completed_steps = step;
        self.inner.publish_wall_exchange();
        Ok(receipt)
    }
}
