//! Mechanical forcing retains all native species and formation-energy axes.

use super::*;
use crate::body_force::{Accumulator, kick};

impl MixtureEuler1d {
    /// Apply a held mechanical source alone, without advancing the flow clock.
    ///
    /// All partial densities remain bit-identical. Source-stage mixture mass
    /// converts acceleration into force; a force density instead supplies its
    /// authored impulse directly. Actual momentum increments and midpoint
    /// velocities supply total-energy work without adding chemical energy or
    /// heat. Original caloric and Mach admission precede atomic publication.
    /// Cancellation or any refusal preserves every accepted species and cell.
    /// Use [`Self::step_with_body_force`] for a complete coupled time interval.
    pub fn apply_body_force(
        &mut self,
        duration: Time,
        force: BodyForce<'_>,
        controls: ForceControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForceChange, Error> {
        poll(&mut cancelled)?;
        controls.admit::<Error>()?;
        if !positive(duration.value()) {
            return Err(Error::Invalid(
                "body-force duration must be positive and finite",
            ));
        }
        if !force.admit::<Error>(self.cell_count(), &mut cancelled)? {
            return Ok(ForceChange::default());
        }
        let m = self.gas.species_count() + 2;
        let mut exchange = Accumulator::default();
        for (i, (q, next)) in self
            .states
            .chunks_exact(m)
            .zip(self.work.candidate.chunks_exact_mut(m))
            .enumerate()
        {
            next.copy_from_slice(q);
            update(
                &self.gas,
                next,
                &mut self.work.primitive[i],
                i,
                force,
                duration.value(),
                controls,
                self.mesh.volumes()[i],
                &mut exchange,
                &mut cancelled,
            )?;
        }
        totals(
            &self.work.candidate,
            self.mesh.volumes(),
            m - 2,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.states, &mut self.work.candidate);
        Ok(exchange.finish())
    }

    /// Advance held body forcing and native species transport atomically.
    ///
    /// Exact force half-steps surround the original MUSCL/HLLC SSP-RK2 map.
    /// Both source maps use their own actual mixture density; no frozen
    /// composition or separately evolved total-density field is introduced.
    /// Source endpoints retain the original EOS, explicit Mach bound and full
    /// interval CFL. All retries restart from the accepted species field, and
    /// only a complete admitted field, clock and source report are published.
    /// One flat source buffer is allocated on first nonzero use and reused.
    /// Zero forcing executes the ordinary step exactly without that allocation.
    ///
    /// `source_change` enters the returned first-law/momentum residual once.
    /// Pressure flow work remains in the original enthalpy boundary flux. This
    /// stationary-mesh method does not include moving walls, chemistry, diffusion
    /// or finite network feedback. State- or time-dependent forcing is held for
    /// this complete accepted interval and requires temporal refinement.
    pub fn step_with_body_force(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        poll(&mut cancelled)?;
        force_controls.admit::<Error>()?;
        if !force.admit::<Error>(self.cell_count(), &mut cancelled)? {
            return self.step(maximum_duration, controls, cancelled);
        }
        if self.body_force_stage.is_none() {
            self.body_force_stage = Some(allocated(self.states.len(), 0.0)?);
        }
        self.advance(
            maximum_duration,
            controls,
            [None; 2],
            Some((force, force_controls)),
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
pub(in crate::multicomponent) fn update(
    gas: &MixtureGas,
    state: &mut [f64],
    primitive: &mut Primitive,
    cell: usize,
    force: BodyForce<'_>,
    duration: f64,
    controls: ForceControl,
    volume: f64,
    exchange: &mut Accumulator,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    poll(cancelled)?;
    let n = gas.species_count();
    let before = gas.decode(state, Some(primitive.temperature), cancelled)?;
    controls.admit_velocity::<Error>(before.velocity, before.sound_speed)?;
    let source = kick::<Error>(
        force,
        cell,
        before.density,
        state[n],
        state[n + 1],
        duration,
    )?;
    let unchanged = source.momentum == state[n] && source.energy == state[n + 1];
    state[n] = source.momentum;
    state[n + 1] = source.energy;
    *primitive = if unchanged {
        before
    } else {
        gas.decode(state, Some(before.temperature), cancelled)?
    };
    controls.admit_velocity::<Error>(primitive.velocity, primitive.sound_speed)?;
    exchange.add::<Error>(&source, volume)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt(
    gas: &MixtureGas,
    states: &[f64],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    reservoirs: ReservoirEnds<'_>,
    work: &mut Workspace,
    forced: &mut [f64],
    force: BodyForce<'_>,
    controls: ForceControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage, ForceChange, f64), Error> {
    let m = gas.species_count() + 2;
    let mut exchange = Accumulator::default();
    let mut signal_rate = 0.0_f64;
    for (i, (q, next)) in states
        .chunks_exact(m)
        .zip(forced.chunks_exact_mut(m))
        .enumerate()
    {
        next.copy_from_slice(q);
        update(
            gas,
            next,
            &mut work.primitive[i],
            i,
            force,
            0.5 * dt,
            controls,
            mesh.volumes()[i],
            &mut exchange,
            cancelled,
        )?;
        admit_cfl(mesh, i, work.primitive[i], dt, cfl, &mut signal_rate)?;
    }
    let (a, b) = super::attempt(
        gas,
        forced,
        mesh,
        boundaries,
        dt,
        cfl,
        reconstruction,
        reservoirs,
        work,
        cancelled,
    )?;
    for (i, q) in work.candidate.chunks_exact_mut(m).enumerate() {
        update(
            gas,
            q,
            &mut work.primitive[i],
            i,
            force,
            0.5 * dt,
            controls,
            mesh.volumes()[i],
            &mut exchange,
            cancelled,
        )?;
        admit_cfl(mesh, i, work.primitive[i], dt, cfl, &mut signal_rate)?;
    }
    Ok((a, b, exchange.finish(), signal_rate))
}

fn admit_cfl(
    mesh: &Mesh1d,
    cell: usize,
    state: Primitive,
    dt: f64,
    maximum: f64,
    signal_rate: &mut f64,
) -> Result<(), Error> {
    let area = 0.5 * mesh.face_areas()[cell] + 0.5 * mesh.face_areas()[cell + 1];
    let rate = (state.velocity.abs() + state.sound_speed) * (area / mesh.volumes()[cell]);
    let observed = rate * dt;
    if !observed.is_finite() || observed > maximum * (1.0 + 64.0 * f64::EPSILON) {
        return Err(Error::CflExceeded { observed, maximum });
    }
    *signal_rate = signal_rate.max(rate);
    Ok(())
}
