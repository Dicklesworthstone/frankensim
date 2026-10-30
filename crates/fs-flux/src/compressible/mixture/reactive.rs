//! Atomic finite-rate chemistry in the original evolving-species Cartesian gas.
//!
//! The original `fs-gas` owner advances each cell at fixed volume, both momenta
//! and full chemical energy. This adapter only composes its two half intervals
//! with the shared conservative Cartesian transport. No reaction, thermodynamic
//! conversion, heat-release law or second spatial discretization is introduced.

use super::{
    BoundaryImpulse, BoundarySide, MAX_MIXTURE_WORKSPACE_BYTES, MixtureEuler2d, StepControl,
    StepReport as TransportStepReport, Totals, WIDTH_LIMIT, allocated, cell, checked_totals, sum,
};
pub use fs_gas::multicomponent::reactive::{CellChemistry, ChemistryReport, ReactiveError};
use fs_qty::{Time, Volume};

/// Original advection controls and bounded retries of the complete split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveStepControl {
    /// Original Cartesian CFL, reconstruction and advection retry policy.
    pub transport: StepControl,
    /// Complete-step restarts after heat release changes the acoustic interval.
    /// At most 32; chemistry itself retains its separately admitted work limits.
    pub maximum_coupling_retries: usize,
}
impl Default for ReactiveStepControl {
    fn default() -> Self {
        Self {
            transport: StepControl::default(),
            maximum_coupling_retries: 16,
        }
    }
}

/// One complete chemistry/advection/chemistry interval and its physical balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveStepReport {
    /// Actual common duration of all three phases.
    pub elapsed: Time,
    /// Rejected complete-step attempts before the accepted one.
    pub coupling_retries: usize,
    /// Original advection phase; its before/after states lie between chemistry.
    pub transport: TransportStepReport,
    /// Accepted chemistry changes, excluding every rejected complete attempt.
    pub chemistry: ChemistryReport,
    /// Complete original field.
    pub before: Totals,
    /// Complete final field.
    pub after: Totals,
    /// Inward physical boundary transfer; chemistry adds no boundary flux.
    pub boundary_change: Totals,
    /// `after-before-boundary-chemistry` for each conserved component. Chemistry
    /// changes species and the referenced energy, but neither momentum component.
    pub conservation_residual: Totals,
    /// Total mass change minus boundary mass transfer [kg].
    pub mass_residual_kg: f64,
    /// Full chemical/thermal/kinetic energy change minus boundary transfer [J].
    /// Uses differences plus original species references, avoiding subtraction
    /// of two large absolute formation-energy totals.
    pub absolute_energy_residual_j: f64,
}

/// Source-owned reactions coupled to the shared two-dimensional gas field.
///
/// All buffers are prepared once. Complete-step retries restore the original
/// accepted field before repeating chemistry, transport and chemistry at one
/// duration. Failure and cancellation leave cells, time, loads and reports
/// unchanged. The chemistry owner may modify its private evaluation scratch.
///
/// The original reactor's selected backward Euler or SDIRK2 method is retained.
/// Symmetric composition alone is not a coupled temporal-order certificate;
/// limiting and first-order advection retries may also reduce observed order.
/// No global accuracy, entropy, flame or experimental-validation certificate is
/// implied. This stationary inviscid lane does not include species diffusion.
pub struct ReactingEuler2d {
    accepted: MixtureEuler2d,
    trial: MixtureEuler2d,
    chemistry: CellChemistry,
    workspace_bytes: usize,
    last_step: Option<ReactiveStepReport>,
}
impl ReactingEuler2d {
    /// Total numerical heap requirement: both full Cartesian workspaces,
    /// including existing accepted-vector excess capacity, and the prepared
    /// chemistry workspace. Immutable sources/prepared caloric coefficients and
    /// bounded stack scratch are excluded, as in the original owners.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler2d,
        chemistry: &CellChemistry,
    ) -> Result<usize, ReactiveError> {
        flow.workspace_bytes()
            .checked_add(MixtureEuler2d::required_workspace_bytes(
                flow.gas(),
                flow.grid(),
            )?)
            .and_then(|bytes| bytes.checked_add(chemistry.workspace_bytes()))
            .ok_or(ReactiveError::Invalid(
                "reacting Cartesian workspace size overflow",
            ))
    }

    /// Bind the original field and exact source/domain-bound prepared chemistry.
    /// The caller's complete numerical byte budget and the shared 1 GiB ceiling
    /// are checked before allocating the duplicate Cartesian workspace.
    pub fn new(
        flow: MixtureEuler2d,
        chemistry: CellChemistry,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        poll(&mut cancelled)?;
        if flow.gas() != chemistry.gas() {
            return Err(ReactiveError::Invalid(
                "Cartesian flow and chemistry have different source or caloric paths",
            ));
        }
        let bytes = Self::required_workspace_bytes(&flow, &chemistry)?;
        if bytes > maximum_workspace_bytes || bytes > MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(ReactiveError::Invalid(
                "reacting Cartesian total workspace budget",
            ));
        }
        let controls = chemistry.controls();
        for index in 0..flow.cell_count() {
            let primitive = flow
                .cell(index)
                .expect("admitted index")
                .primitive(flow.gas(), &mut cancelled)?;
            if !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
                .contains(&primitive.temperature)
            {
                return Err(ReactiveError::Invalid(
                    "initial Cartesian temperature is outside chemistry interval",
                ));
            }
        }
        let mut states = allocated(flow.states.len(), 0.0)?;
        copy_values(&flow.states, &mut states, &mut cancelled)?;
        let mut trial = MixtureEuler2d::new(
            flow.gas.clone(),
            flow.grid,
            states,
            flow.boundaries,
            MixtureEuler2d::required_workspace_bytes(flow.gas(), flow.grid())?,
            &mut cancelled,
        )?;
        trial.time = flow.time;
        trial.last_step = flow.last_step;
        copy_values(
            &flow.boundary_impulses,
            &mut trial.boundary_impulses,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            accepted: flow,
            trial,
            chemistry,
            workspace_bytes: bytes,
            last_step: None,
        })
    }

    /// Complete accepted field, original sources, grid and physical clock.
    /// Its Euler metadata describes only the middle advective phase; this
    /// owner's [`Self::last_step`] includes both chemistry half steps.
    pub const fn flow(&self) -> &MixtureEuler2d {
        &self.accepted
    }
    /// Original source-bound chemistry, accessible without mutable scratch.
    pub const fn chemistry(&self) -> &CellChemistry {
        &self.chemistry
    }
    /// Actual combined prepared numerical heap requirement in bytes.
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    /// Most recent completely accepted chemistry/transport interval.
    pub const fn last_step(&self) -> Option<&ReactiveStepReport> {
        self.last_step.as_ref()
    }
    /// Complete outward exchange through each physical face. Chemistry has no
    /// boundary flux, so these are the original integrated advective transfers.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        self.accepted.boundary_impulses(side)
    }
    /// Release the accepted conservative field after all completed reactions.
    pub fn into_flow(self) -> MixtureEuler2d {
        self.accepted
    }

    /// Advance at most the requested duration, atomically publishing a complete
    /// chemistry/advection/chemistry split. If heat release shortens the Euler
    /// interval, both chemistry half steps restart at that same shorter time.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl)
            || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        {
            return Err(ReactiveError::Invalid(
                "reacting Cartesian duration or controls",
            ));
        }
        let n = self.accepted.gas.species_count();
        let before = checked_totals(
            &self.accepted.states,
            &self.accepted.grid,
            n,
            &mut cancelled,
        )?;
        let mut rate: f64 = 0.0;
        for row in self.accepted.states.chunks_exact(n + 3) {
            let p = cell(row, n).primitive(&self.accepted.gas, &mut cancelled)?;
            rate = rate.max(
                (p.velocity[0].abs() + p.sound_speed) / self.accepted.grid.spacing[0]
                    + (p.velocity[1].abs() + p.sound_speed) / self.accepted.grid.spacing[1],
            );
        }
        // A modest pre-chemistry margin avoids paying for predictably rejected
        // acoustic intervals. The actual post-reaction faces still own CFL.
        let mut dt = maximum_duration
            .value()
            .min(0.8 * controls.transport.cfl / rate);
        let mut last = ReactiveError::Invalid("no representable reacting Cartesian increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(0.5 * dt)
                || !self.accepted.time.is_finite()
                || self.accepted.time + dt <= self.accepted.time
                || !(self.accepted.time + dt).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "unrepresentable reacting Cartesian increment",
                ));
            }
            copy_physical(&self.accepted, &mut self.trial, &mut cancelled)?;
            let mut chemistry = ChemistryReport::default();
            let mut shortened = None;
            let attempt = (|| {
                react_cells(
                    &mut self.trial,
                    &mut self.chemistry,
                    0.5 * dt,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let transport =
                    self.trial
                        .step(Time::new(dt), controls.transport, &mut cancelled)?;
                if transport.elapsed.value() != dt {
                    shortened = Some(transport.elapsed.value());
                    return Err(ReactiveError::Invalid(
                        "Euler shortened the reacting interval",
                    ));
                }
                react_cells(
                    &mut self.trial,
                    &mut self.chemistry,
                    0.5 * dt,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let after =
                    checked_totals(&self.trial.states, &self.trial.grid, n, &mut cancelled)?;
                let boundary_change = transport.boundary_change;
                let mut residual = [0.0; WIDTH_LIMIT];
                for k in 0..n {
                    residual[k] = after.species_masses[k]
                        - before.species_masses[k]
                        - boundary_change.species_masses[k]
                        - chemistry.species_change_kg[k];
                }
                for axis in 0..2 {
                    residual[n + axis] = after.momentum[axis]
                        - before.momentum[axis]
                        - boundary_change.momentum[axis];
                }
                residual[n + 2] = after.energy
                    - before.energy
                    - boundary_change.energy
                    - chemistry.referenced_energy_change_j;
                let mass_residual_kg = sum((0..n).map(|k| {
                    after.species_masses[k]
                        - before.species_masses[k]
                        - boundary_change.species_masses[k]
                }));
                let absolute_energy_residual_j = sum(std::iter::once(
                    after.energy - before.energy - boundary_change.energy,
                )
                .chain((0..n).map(|k| {
                    (after.species_masses[k]
                        - before.species_masses[k]
                        - boundary_change.species_masses[k])
                        * self.trial.gas.energy_references()[k]
                })));
                if !residual[..n + 3].iter().all(|value| value.is_finite())
                    || !mass_residual_kg.is_finite()
                    || !absolute_energy_residual_j.is_finite()
                    || !chemistry.entropy_change_j_k.is_finite()
                {
                    return Err(ReactiveError::Invalid(
                        "reacting Cartesian balance overflow",
                    ));
                }
                poll(&mut cancelled)?;
                Ok(ReactiveStepReport {
                    elapsed: Time::new(dt),
                    coupling_retries: retry,
                    transport,
                    chemistry,
                    before,
                    after,
                    boundary_change,
                    conservation_residual: Totals::from_row(&residual[..n + 3], n),
                    mass_residual_kg,
                    absolute_energy_residual_j,
                })
            })();
            match attempt {
                Ok(report) => {
                    std::mem::swap(&mut self.accepted, &mut self.trial);
                    self.last_step = Some(report);
                    return Ok(report);
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(error) => {
                    last = error;
                    dt = shortened.unwrap_or(0.5 * dt);
                }
            }
        }
        Err(ReactiveError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}

fn copy_physical(
    source: &MixtureEuler2d,
    target: &mut MixtureEuler2d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    copy_values(&source.states, &mut target.states, cancelled)?;
    copy_values(
        &source.boundary_impulses,
        &mut target.boundary_impulses,
        cancelled,
    )?;
    target.time = source.time;
    target.last_step = source.last_step;
    Ok(())
}

fn copy_values(
    source: &[f64],
    target: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    for (source, target) in source.chunks(4096).zip(target.chunks_mut(4096)) {
        poll(cancelled)?;
        target.copy_from_slice(source);
    }
    Ok(())
}

fn react_cells(
    flow: &mut MixtureEuler2d,
    chemistry: &mut CellChemistry,
    duration: f64,
    report: &mut ChemistryReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    let width = flow.gas.species_count() + 3;
    let volume = Volume::new(flow.grid.cell_volume());
    for (index, row) in flow.states.chunks_exact_mut(width).enumerate() {
        chemistry.advance_planar(
            row,
            volume,
            Time::new(duration),
            index,
            report,
            &mut *cancelled,
        )?;
    }
    // Extensive admission is independent of individual cell caloric admission.
    checked_totals(&flow.states, &flow.grid, width - 3, cancelled)?;
    poll(cancelled)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReactiveError> {
    if cancelled() {
        Err(ReactiveError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
