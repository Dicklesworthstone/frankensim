//! Source splitting shares the original transport operator and retry owner.

use super::*;
use crate::body_force::{Accumulator, kick};

impl Euler1d {
    /// Apply the exact flow of a held distributed mechanical source alone.
    ///
    /// Mass and the original caloric reference are unchanged. The represented
    /// momentum increment supplies midpoint kinetic work to total energy;
    /// internal energy and temperature are preserved to floating-point accuracy.
    /// The returned energy residual exposes representation loss, never heat.
    /// Every candidate uses the original EOS and caller's Mach bound. Failure
    /// or cancellation leaves every accepted cell unchanged.
    ///
    /// Time does not advance: this is a source substep for a composing owner.
    /// Use [`Self::step_with_body_force`] for atomic source/transport evolution.
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
        if !force.admit::<Error>(self.states.len(), &mut cancelled)? {
            return Ok(ForceChange::default());
        }
        let mut exchange = Accumulator::default();
        for (i, &q) in self.states.iter().enumerate() {
            self.work.candidate[i] = update(
                &self.gas,
                q,
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
        // Also admit represented extensive storage before an infallible totals query.
        totals(&self.work.candidate, self.mesh.volumes(), &mut cancelled)?;
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.states, &mut self.work.candidate);
        Ok(exchange.finish())
    }

    /// Advance held distributed forcing and native Euler transport together.
    ///
    /// Symmetric force/advection/force splitting applies two exact half-duration
    /// source maps around the original MUSCL/HLLC SSP-RK2 update. The same retry
    /// restarts both source halves and transport from the accepted field. Source
    /// endpoints obey the original caloric domain, authored Mach bound and full
    /// interval CFL. Gas, source accounting and clock publish only together.
    /// The first nonzero call allocates one bounded cell buffer, reused thereafter;
    /// zero forcing follows the ordinary step exactly without that allocation.
    ///
    /// `source_change` reports actual volume-integrated impulse and conjugate
    /// work; they enter `conservation_residual` exactly once. A bound atmosphere's
    /// gravity remains a distinct, additive source: do not supply it again here.
    /// Boundary pressure work is already included in the native boundary flux.
    /// This stationary-mesh method does not compose moving-wall, wall-heat or
    /// network updates. Time-varying or state-dependent forces are held at the
    /// authored value for this complete accepted interval and require refinement.
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
        if !force.admit::<Error>(self.states.len(), &mut cancelled)? {
            return self.step(maximum_duration, controls, cancelled);
        }
        if self.body_force_stage.is_none() {
            self.body_force_stage = Some(allocated(self.states.len(), Conserved::default())?);
        }
        self.advance(
            maximum_duration,
            controls,
            None::<&mut NoPorts>,
            Some((force, force_controls)),
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn update(
    gas: &FrozenGas,
    q: Conserved,
    primitive: &mut Primitive,
    cell: usize,
    force: BodyForce<'_>,
    duration: f64,
    controls: ForceControl,
    volume: f64,
    exchange: &mut Accumulator,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Conserved, Error> {
    poll(cancelled)?;
    let before = gas.decode(q, Some(primitive.temperature), cancelled)?;
    controls.admit_velocity::<Error>(before.velocity, before.sound_speed)?;
    let source = kick::<Error>(force, cell, q.density, q.momentum, q.energy, duration)?;
    let next = Conserved {
        density: q.density,
        momentum: source.momentum,
        energy: source.energy,
    };
    *primitive = if next == q {
        before
    } else {
        gas.decode(next, Some(before.temperature), cancelled)?
    };
    controls.admit_velocity::<Error>(primitive.velocity, primitive.sound_speed)?;
    exchange.add::<Error>(&source, volume)?;
    Ok(next)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attempt<P>(
    gas: &FrozenGas,
    states: &[Conserved],
    mesh: &Mesh1d,
    boundaries: Boundaries,
    time: [f64; 2],
    cfl: f64,
    reconstruction: Reconstruction,
    hydrostatic: Option<&HydrostaticBalance>,
    work: &mut Workspace,
    forced: &mut [Conserved],
    force: BodyForce<'_>,
    controls: ForceControl,
    ports: Option<&mut P>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Stage, Stage, ForceChange, f64), Error>
where
    P: FnMut(&PortStage) -> Result<[Option<FaceFlux>; 2], Error>,
{
    let dt = time[1];
    let mut exchange = Accumulator::default();
    let mut signal_rate = 0.0_f64;
    for (i, &q) in states.iter().enumerate() {
        forced[i] = update(
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
    let (a, b) = super::attempt(
        gas,
        forced,
        mesh,
        boundaries,
        time,
        cfl,
        reconstruction,
        hydrostatic,
        work,
        ports,
        cancelled,
    )?;
    for (i, q) in work.candidate.iter_mut().enumerate() {
        *q = update(
            gas,
            *q,
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
