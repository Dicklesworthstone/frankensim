//! Exact mechanical kicks on the original density-dependent caloric owner.

use super::*;
use crate::body_force::{Failure, ForceChange, kick};

impl Failure for RealGasFlowError {
    fn invalid(reason: &'static str) -> Self {
        Self::Invalid(reason)
    }
    fn cancelled() -> Self {
        Self::Cancelled
    }
}

pub(super) struct Workspace {
    pub(super) states: Vec<Conserved>,
    pub(super) primitive: Vec<Primitive>,
}

impl Workspace {
    pub(super) fn new(cells: usize) -> Result<Self, RealGasFlowError> {
        Ok(Self {
            states: allocated(cells, Conserved::default())?,
            primitive: allocated(cells, Primitive::default())?,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        &mut self,
        gas: &RealGas,
        states: &[Conserved],
        primitive: &[Primitive],
        mesh: &Mesh1d,
        force: BodyForce<'_>,
        controls: ForceControl,
        dt: f64,
        cfl: f64,
        exchange: &mut Accumulator,
        rate: &mut f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RealGasFlowError> {
        for (i, &q) in states.iter().enumerate() {
            (self.states[i], self.primitive[i]) = update(
                gas,
                q,
                primitive[i],
                force,
                controls,
                i,
                0.5 * dt,
                mesh.volumes()[i],
                exchange,
                cancelled,
            )?;
            admit_cfl(self.primitive[i], mesh, i, dt, cfl, rate)?;
        }
        Ok(())
    }
}

impl RealGasEuler1d {
    /// Apply a held acceleration or force density without advancing the clock.
    ///
    /// The original shared kick changes momentum and total energy together:
    /// `delta E = delta(rho*u) * (u_before + u_after)/2`. Density and caloric
    /// reference stay fixed, so this source supplies kinetic work, not heat.
    /// Every candidate resolves pressure, temperature and sound speed through
    /// the retained NASA/Peng--Robinson model and the caller's Mach limit.
    /// All accepted cells and cached primitives publish together or not at all.
    pub fn apply_body_force(
        &mut self,
        duration: Time,
        force: BodyForce<'_>,
        controls: ForceControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForceChange, RealGasFlowError> {
        poll(&mut cancelled)?;
        controls.admit::<RealGasFlowError>()?;
        if !positive(duration.value()) {
            return Err(RealGasFlowError::Invalid(
                "body-force duration must be positive and finite",
            ));
        }
        if !force.admit::<RealGasFlowError>(self.states.len(), &mut cancelled)? {
            return Ok(ForceChange::default());
        }
        let mut exchange = Accumulator::default();
        for (i, &q) in self.states.iter().enumerate() {
            (self.work.candidate[i], self.work.candidate_primitive[i]) = update(
                &self.gas,
                q,
                self.primitive[i],
                force,
                controls,
                i,
                duration.value(),
                self.mesh.volumes()[i],
                &mut exchange,
                &mut cancelled,
            )?;
        }
        let total = totals(&self.work.candidate, self.mesh.volumes(), &mut cancelled)?;
        if total.array().iter().any(|x| !x.is_finite()) {
            return Err(RealGasFlowError::Invalid(
                "nonfinite forced real-gas storage",
            ));
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.states, &mut self.work.candidate);
        core::mem::swap(&mut self.primitive, &mut self.work.candidate_primitive);
        Ok(exchange.finish())
    }

    /// Compose exact held-source half steps around the original SSP-RK2 waves.
    ///
    /// Acceleration uses density at EACH source stage; force per volume does
    /// not acquire an extra density factor. Original EOS, endpoint Mach and
    /// full-interval CFL admission cover both kicks. A retry restarts from the
    /// accepted field, and commits gas, primitive cache, source report and clock
    /// together. Two cell buffers are allocated on first nonzero use and reused.
    /// Zero force follows the ordinary wave path without these allocations.
    ///
    /// The source is held over the complete attempted interval. This does not
    /// provide a state-dependent force law, hydrostatic well balancing, moving
    /// walls or a new spatial reconstruction. Time refinement remains required.
    pub fn step_with_body_force(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        self.step_inner(
            maximum_duration,
            controls,
            None,
            Some((force, force_controls)),
            cancelled,
        )
    }

    /// Add original finite-wall heat/drag after the forced wave interval.
    ///
    /// The two force halves surround transport; wall exchange then uses the
    /// actual accepted duration. This wall/source composition is first order.
    /// Body impulse/work, wall impulse/heat/drive work and boundary transport
    /// stay separate and enter the gas balance once. The gas plus finite-wall
    /// storage balance subtracts external body work as well as wall drive work.
    /// Refusal or cancellation commits neither participant nor either report.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_walls_and_body_force(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        walls: &mut RealGasWallOperator<'_>,
        wall_controls: WallControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, RealGasFlowError> {
        self.step_inner(
            maximum_duration,
            controls,
            Some((walls, wall_controls)),
            Some((force, force_controls)),
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn update(
    gas: &RealGas,
    q: Conserved,
    before: Primitive,
    force: BodyForce<'_>,
    controls: ForceControl,
    cell: usize,
    duration: f64,
    volume: f64,
    exchange: &mut Accumulator,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Conserved, Primitive), RealGasFlowError> {
    poll(cancelled)?;
    controls.admit_velocity::<RealGasFlowError>(before.velocity, before.sound_speed)?;
    let source = kick::<RealGasFlowError>(force, cell, q.density, q.momentum, q.energy, duration)?;
    let next = Conserved {
        density: q.density,
        momentum: source.momentum,
        energy: source.energy,
    };
    let primitive = if next == q {
        before
    } else {
        gas.decode(next, cancelled)?
    };
    controls.admit_velocity::<RealGasFlowError>(primitive.velocity, primitive.sound_speed)?;
    exchange.add::<RealGasFlowError>(&source, volume)?;
    Ok((next, primitive))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn finish(
    gas: &RealGas,
    states: &mut [Conserved],
    primitive: &mut [Primitive],
    mesh: &Mesh1d,
    force: BodyForce<'_>,
    controls: ForceControl,
    dt: f64,
    cfl: f64,
    exchange: &mut Accumulator,
    rate: &mut f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), RealGasFlowError> {
    for (i, q) in states.iter_mut().enumerate() {
        (*q, primitive[i]) = update(
            gas,
            *q,
            primitive[i],
            force,
            controls,
            i,
            0.5 * dt,
            mesh.volumes()[i],
            exchange,
            cancelled,
        )?;
        admit_cfl(primitive[i], mesh, i, dt, cfl, rate)?;
    }
    Ok(())
}

fn admit_cfl(
    primitive: Primitive,
    mesh: &Mesh1d,
    cell: usize,
    dt: f64,
    maximum: f64,
    rate: &mut f64,
) -> Result<(), RealGasFlowError> {
    let local = (primitive.velocity.abs() + primitive.sound_speed)
        * (mesh.face_areas()[cell] / mesh.volumes()[cell]);
    cfl_gate(local * dt, maximum)?;
    *rate = rate.max(local);
    Ok(())
}
