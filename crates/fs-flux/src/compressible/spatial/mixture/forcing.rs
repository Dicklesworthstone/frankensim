//! Mechanical source kicks on the original conservative mixture transaction.
//!
//! Held acceleration uses the actual density of each source stage; a force
//! density is already per volume. Neither changes species or formation-energy
//! references. Two exact kinetic-work kicks surround the original SSP RK2
//! advection, and every failed trial restarts from accepted cells. This is
//! prescribed forcing, without a hydrostatic well-balanced discretization or
//! an implicit velocity-dependent force solve.

use super::super::forcing::kick;
use super::*;

pub(super) type Context<'a, 'b> = (BodyForce3d<'a>, ForceControl3d, &'b mut ForcingWorkspace3d);

/// Reusable original state after the first half kick; no second flow owner.
/// The numerical budget excludes fixed metadata and allocator overhead.
pub struct ForcingWorkspace3d {
    shape: [usize; 3],
    species: usize,
    pub(super) kicked: Vec<f64>,
}

impl ForcingWorkspace3d {
    pub fn workspace_bytes(&self) -> usize {
        self.kicked.capacity() * size_of::<f64>()
    }

    pub(super) fn admit(&self, flow: &MixtureEuler3d) -> Result<(), Error> {
        if self.shape != flow.grid.shape
            || self.species != flow.gas.species_count()
            || self.kicked.len() != flow.states.len()
            || flow
                .workspace_bytes()
                .checked_add(self.workspace_bytes())
                .is_none_or(|bytes| bytes > MAX_MIXTURE_WORKSPACE_BYTES)
        {
            return Err(Error::Invalid(
                "mixture forcing workspace geometry or total budget",
            ));
        }
        Ok(())
    }
}

impl MixtureEuler3d {
    /// Additional numerical bytes: one species/momentum/energy row per cell.
    pub fn required_forcing_workspace_bytes(&self) -> Result<usize, Error> {
        self.states
            .len()
            .checked_mul(size_of::<f64>())
            .ok_or(Error::Invalid("mixture forcing workspace size overflow"))
    }

    /// Allocate once within the explicit additional-byte budget. The original
    /// field plus this scratch must also obey the mixture's total memory cap.
    pub fn forcing_workspace(&self, maximum_bytes: usize) -> Result<ForcingWorkspace3d, Error> {
        let required = self.required_forcing_workspace_bytes()?;
        if required > maximum_bytes
            || self
                .workspace_bytes()
                .checked_add(required)
                .is_none_or(|bytes| bytes > MAX_MIXTURE_WORKSPACE_BYTES)
        {
            return Err(Error::Invalid("mixture forcing workspace byte budget"));
        }
        let workspace = ForcingWorkspace3d {
            shape: self.grid.shape,
            species: self.gas.species_count(),
            kicked: allocated(self.states.len(), 0.0)?,
        };
        if workspace.workspace_bytes() > maximum_bytes {
            return Err(Error::Invalid(
                "mixture forcing capacity exceeds byte budget",
            ));
        }
        workspace.admit(self)?;
        Ok(workspace)
    }

    /// Advance held physical forcing with the original advection and clock.
    /// A zero field follows the original numerical path exactly. Nonzero
    /// forcing admits the authored Mach bound at every source and RK stage;
    /// its endpoint acoustic cell rate supplements the original face CFL.
    /// Source work and impulse are reported separately from boundary exchange.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        force_controls.admit::<Error>()?;
        workspace.admit(self)?;
        if !force.admit::<Error>(self.cell_count(), &mut cancelled)? {
            return self.step(maximum_duration, controls, cancelled);
        }
        self.step_impl(
            maximum_duration,
            controls,
            Some((force, force_controls, workspace)),
            cancelled,
        )
    }
}

/// Fixed cell order, with one compensated sum for each mechanical source term.
#[derive(Default)]
pub(super) struct Accumulator {
    values: [f64; 4],
    corrections: [f64; 4],
}

impl Accumulator {
    fn add(&mut self, impulse: [f64; 3], work: f64, volume: f64) -> Result<(), Error> {
        for (k, density) in impulse.into_iter().chain([work]).enumerate() {
            let value = density * volume;
            if !value.is_finite() {
                return Err(Error::Invalid(
                    "unrepresentable extensive mixture force source",
                ));
            }
            let y = value - self.corrections[k];
            let next = self.values[k] + y;
            self.corrections[k] = (next - self.values[k]) - y;
            self.values[k] = next;
        }
        Ok(())
    }

    pub(super) fn finish(self) -> Result<ForceChange3d, Error> {
        if !self
            .values
            .into_iter()
            .chain(self.corrections)
            .all(f64::is_finite)
        {
            return Err(Error::Invalid("mixture force source reduction overflow"));
        }
        Ok(ForceChange3d {
            impulse: [self.values[0], self.values[1], self.values[2]],
            work: self.values[3],
        })
    }
}

pub(super) fn acoustic_rate(p: Primitive, grid: &CartesianGrid) -> Result<f64, Error> {
    let rate: f64 = (0..3)
        .map(|axis| (p.velocity[axis].abs() + p.sound_speed) / grid.spacing[axis])
        .sum();
    if !rate.is_finite() || rate <= 0.0 {
        return Err(Error::Invalid(
            "unrepresentable forced mixture acoustic rate",
        ));
    }
    Ok(rate)
}

/// Admit an enclosing composition's final physical state after chemistry or
/// molecular heating, rather than only its central forced Euler endpoint.
pub(super) fn admit_field(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    states: &[f64],
    control: ForceControl3d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut rate: f64 = 0.0;
    for (i, row) in states.chunks_exact(gas.species_count() + 4).enumerate() {
        if i % 64 == 0 {
            stepping::poll::<Error>(cancelled)?;
        }
        let p = cell(row, gas.species_count()).primitive(gas, &mut *cancelled)?;
        control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
        rate = rate.max(acoustic_rate(p, grid)?);
    }
    stepping::poll::<Error>(cancelled)?;
    Ok(rate)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn kick_field(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    states: &mut [f64],
    force: BodyForce3d<'_>,
    control: ForceControl3d,
    duration: f64,
    source: &mut Accumulator,
    admit_endpoint: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let n = gas.species_count();
    let mut rate: f64 = 0.0;
    for (i, row) in states.chunks_exact_mut(n + 4).enumerate() {
        if i % 64 == 0 {
            stepping::poll::<Error>(cancelled)?;
        }
        // Inputs are original or already source-admitted RK states. The mass
        // remains the compensated species sum, never a separately evolved rho.
        let rho = sum(row[..n].iter().copied());
        let next = kick::<Error, 3>(
            force,
            i,
            rho,
            [row[n], row[n + 1], row[n + 2]],
            row[n + 3],
            duration,
        )?;
        row[n..n + 3].copy_from_slice(&next.momentum);
        row[n + 3] = next.energy;
        if admit_endpoint {
            let p = cell(row, n).primitive(gas, &mut *cancelled)?;
            control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            rate = rate.max(acoustic_rate(p, grid)?);
        }
        source.add(next.impulse_density, next.work_density, grid.cell_volume())?;
    }
    stepping::poll::<Error>(cancelled)?;
    Ok(rate)
}

#[cfg(test)]
#[path = "forcing_tests.rs"]
mod tests;
