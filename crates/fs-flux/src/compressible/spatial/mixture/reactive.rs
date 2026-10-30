//! Finite-rate chemistry coupled to the original volumetric species field.
//!
//! The source-owned prepared reactor advances every cell at fixed volume and
//! all three momenta. Chemistry/waves/chemistry use one actual accepted time
//! interval; heat-release-induced CFL reductions restart the complete trial.
//! This stationary inviscid engineering lane does not model diffusion, moving
//! boundaries, physical detonation thickness or a fully discrete entropy law.

use super::super::{BoundarySide, StepControl};
use super::{BoundaryImpulse, MixtureEuler3d, StepReport, Totals, allocated, cell, checked_totals};
pub use fs_gas::chemistry::IsochoricKineticsSolve;
use fs_gas::multicomponent::MAX_SPECIES;
pub use fs_gas::multicomponent::reactive::{
    CellChemistry, ChemistryReport, IdealGasKineticsV1, ReactiveError,
};
use fs_qty::{Time, Volume};

/// Bounded whole-interval restarts and original wave controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveStepControl {
    /// Original unsplit three-dimensional transport controls.
    pub transport: StepControl,
    /// Restarts permitted after the initial attempt; at most 32.
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

/// Complete accepted chemistry/wave/chemistry interval and physical balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveStepReport3d {
    /// Actual duration shared by both chemistry halves and the wave step.
    pub elapsed: Time,
    /// Rejected complete trials before this accepted interval.
    pub coupling_retries: usize,
    /// Middle transport report; its before/after totals exclude outer chemistry.
    pub transport: StepReport,
    /// Accepted source changes and work over both chemistry halves.
    pub chemistry: ChemistryReport,
    /// Original complete field totals.
    pub before: Totals,
    /// Final complete field totals.
    pub after: Totals,
    /// Inward boundary transfer for the complete physical interval.
    /// Chemistry is local and contributes no face flux.
    pub boundary_change: Totals,
    /// Species, momentum and referenced-energy residual after source and boundary terms.
    pub conservation_residual: Totals,
    /// Total mass change minus boundary mass transfer [kg].
    pub mass_residual_kg: f64,
    /// Chemical-plus-kinetic energy change minus boundary energy transfer [J].
    /// Differences are formed before applying the original species references.
    pub absolute_energy_residual_j: f64,
}

/// Atomic conservative three-dimensional reacting compressible gas.
///
/// Two preallocated fields and one source-owned prepared cell reactor are
/// reused for every step. All accepted cells, clock, last report and boundary
/// transfers remain unchanged until the complete split passes admission.
/// The selected chemistry method controls source accuracy; symmetric ordering
/// alone does not certify second-order shock convergence or entropy stability.
pub struct ReactingEuler3d {
    accepted: MixtureEuler3d,
    trial: MixtureEuler3d,
    reactor: CellChemistry,
    field_workspace_bytes: usize,
    workspace_bytes: usize,
    last_step: Option<ReactiveStepReport3d>,
}

impl ReactingEuler3d {
    /// Buffer bytes for both accepted/trial fields, including the accepted
    /// input's excess capacity. The original immutable gas/mechanism models,
    /// bounded reactor workspace and fixed-size reports are excluded.
    pub fn required_field_workspace_bytes(flow: &MixtureEuler3d) -> Result<usize, ReactiveError> {
        flow.workspace_bytes()
            .checked_add(MixtureEuler3d::required_workspace_bytes(
                flow.gas(),
                flow.grid(),
            )?)
            .ok_or(ReactiveError::Invalid(
                "reacting spatial workspace size overflow",
            ))
    }

    /// Complete numerical heap capacity for accepted/trial fields and the
    /// original prepared chemistry workspace. Immutable material/source data,
    /// prepared caloric coefficients and bounded stack scratch are excluded.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler3d,
        chemistry: &CellChemistry,
    ) -> Result<usize, ReactiveError> {
        Self::required_field_workspace_bytes(flow)?
            .checked_add(chemistry.workspace_bytes())
            .ok_or(ReactiveError::Invalid(
                "reacting spatial workspace size overflow",
            ))
    }

    /// Bind the exact original mechanism while retaining the accepted field.
    ///
    /// `maximum_workspace_bytes` bounds both field workspaces and the original
    /// prepared chemistry buffers together, under the original hard ceiling.
    /// Chemistry also retains its original species/reaction and work bounds.
    pub fn new(
        flow: MixtureEuler3d,
        mechanism: IdealGasKineticsV1,
        chemistry: IsochoricKineticsSolve,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        poll(&mut cancelled)?;
        if Self::required_field_workspace_bytes(&flow)? > maximum_workspace_bytes {
            return Err(ReactiveError::Invalid("reacting spatial workspace budget"));
        }
        let reactor = CellChemistry::new(flow.gas().clone(), mechanism, chemistry, &mut cancelled)?;
        Self::with_chemistry(flow, reactor, maximum_workspace_bytes, cancelled)
    }

    /// Compose an already prepared source without rebuilding its mechanism or
    /// reactor. Both owners must retain the same original gas and caloric gauge.
    pub fn with_chemistry(
        flow: MixtureEuler3d,
        reactor: CellChemistry,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        poll(&mut cancelled)?;
        let field_workspace_bytes = Self::required_field_workspace_bytes(&flow)?;
        let workspace_bytes = Self::required_workspace_bytes(&flow, &reactor)?;
        if workspace_bytes > maximum_workspace_bytes
            || workspace_bytes > super::MAX_MIXTURE_WORKSPACE_BYTES
        {
            return Err(ReactiveError::Invalid("reacting spatial workspace budget"));
        }
        if flow.gas() != reactor.gas() {
            return Err(ReactiveError::Invalid(
                "reacting field and cell source differ",
            ));
        }
        let chemistry = reactor.controls();
        let n = flow.gas().species_count();
        let width = n + 4;
        for q in flow.states.chunks_exact(width) {
            let p = cell(q, n).primitive(flow.gas(), &mut cancelled)?;
            if !(chemistry.minimum_temperature_k..=chemistry.maximum_temperature_k)
                .contains(&p.temperature)
            {
                return Err(ReactiveError::Invalid(
                    "initial spatial temperature is outside chemistry interval",
                ));
            }
        }
        let mut cells = allocated(flow.states.len(), 0.0)?;
        copy_field(&flow.states, &mut cells, &mut cancelled)?;
        let mut trial = MixtureEuler3d::new(
            flow.gas().clone(),
            *flow.grid(),
            cells,
            flow.boundaries(),
            maximum_workspace_bytes - flow.workspace_bytes() - reactor.workspace_bytes(),
            &mut cancelled,
        )?;
        copy_physical(&flow, &mut trial, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            accepted: flow,
            trial,
            reactor,
            field_workspace_bytes,
            workspace_bytes,
            last_step: None,
        })
    }

    /// Accepted physical field and clock. Its transport report describes only
    /// the middle wave stage; this owner's report covers the complete interval.
    pub const fn flow(&self) -> &MixtureEuler3d {
        &self.accepted
    }

    /// Source-owned finite-rate mechanism retained without reconstruction.
    pub fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.reactor.mechanism()
    }

    /// Unchanged source tolerances, time method, interval and work limits.
    pub fn chemistry_controls(&self) -> IsochoricKineticsSolve {
        self.reactor.controls()
    }

    /// Admitted accepted-plus-trial field-buffer bytes.
    pub const fn field_workspace_bytes(&self) -> usize {
        self.field_workspace_bytes
    }

    /// Complete admitted numerical heap capacity, including chemistry buffers.
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Last completely published reacting interval.
    pub const fn last_step(&self) -> Option<&ReactiveStepReport3d> {
        self.last_step.as_ref()
    }

    /// Complete outward face transfer, in the original physical face order.
    /// The two local chemistry stages have exactly zero boundary flux.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        self.accepted.boundary_impulses(side)
    }

    /// Release the complete accepted conservative field.
    pub fn into_flow(self) -> MixtureEuler3d {
        self.accepted
    }

    /// Advance at most the requested duration with an atomic C/A/C composition.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport3d, ReactiveError> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || controls.maximum_coupling_retries > 32
            || !positive(controls.transport.cfl)
            || controls.transport.cfl > 0.45
            || controls.transport.maximum_halvings > 32
        {
            return Err(ReactiveError::Invalid(
                "reacting duration, CFL or retry controls",
            ));
        }
        let n = self.accepted.gas.species_count();
        let width = n + 4;
        let before = checked_totals(
            &self.accepted.states,
            &self.accepted.grid,
            n,
            &mut cancelled,
        )?;
        let spacing = self.accepted.grid.spacing();
        let mut rate: f64 = 0.0;
        for q in self.accepted.states.chunks_exact(width) {
            let p = cell(q, n).primitive(&self.accepted.gas, &mut cancelled)?;
            let local = (0..3)
                .map(|axis| (p.velocity[axis].abs() + p.sound_speed) / spacing[axis])
                .sum::<f64>();
            rate = rate.max(local);
        }
        // Leave acoustic headroom before paying for the first reactive half.
        // Heat release can still increase the true face CFL; that is retried.
        let mut dt = maximum_duration
            .value()
            .min(0.8 * controls.transport.cfl / rate);
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(0.5 * dt)
                || !(self.accepted.time + dt).is_finite()
                || !positive(self.accepted.time + dt - self.accepted.time)
            {
                return Err(ReactiveError::Invalid(
                    "unrepresentable reacting time increment",
                ));
            }
            copy_physical(&self.accepted, &mut self.trial, &mut cancelled)?;
            let mut chemistry = ChemistryAccumulator::default();
            let mut shortened = None;
            let attempt = (|| {
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
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
                        "wave CFL shortened the reacting interval",
                    ));
                }
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    0.5 * dt,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let after =
                    checked_totals(&self.trial.states, &self.trial.grid, n, &mut cancelled)?;
                let chemistry = chemistry.report;
                let boundary_change = transport.boundary_change;
                let mut delta = [0.0; MAX_SPECIES + 4];
                for d in 0..width {
                    delta[d] =
                        after.component(d) - before.component(d) - boundary_change.component(d);
                }
                let mass_residual_kg = sum(delta[..n].iter().copied());
                let absolute_energy_residual_j = sum(std::iter::once(delta[n + 3])
                    .chain((0..n).map(|k| delta[k] * self.trial.gas.energy_references()[k])));
                for k in 0..n {
                    delta[k] -= chemistry.species_change_kg[k];
                }
                delta[n + 3] -= chemistry.referenced_energy_change_j;
                if !delta[..width].iter().all(|v| v.is_finite())
                    || !mass_residual_kg.is_finite()
                    || !absolute_energy_residual_j.is_finite()
                    || !chemistry.entropy_change_j_k.is_finite()
                {
                    return Err(ReactiveError::Invalid(
                        "non-finite spatial reacting balance",
                    ));
                }
                let report = ReactiveStepReport3d {
                    elapsed: Time::new(dt),
                    coupling_retries: retry,
                    transport,
                    chemistry,
                    before,
                    after,
                    boundary_change,
                    conservation_residual: Totals::from_row(&delta[..width], n),
                    mass_residual_kg,
                    absolute_energy_residual_j,
                };
                poll(&mut cancelled)?;
                Ok(report)
            })();
            match attempt {
                Ok(report) => {
                    std::mem::swap(&mut self.accepted, &mut self.trial);
                    self.last_step = Some(report);
                    return Ok(report);
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(last) => {
                    if retry == controls.maximum_coupling_retries {
                        return Err(ReactiveError::RetryLimit {
                            attempts: retry + 1,
                            last: Box::new(last),
                        });
                    }
                    dt = shortened.unwrap_or(0.5 * dt);
                }
            }
        }
        unreachable!("bounded retries return on their last attempt")
    }
}

fn copy_field(
    source: &[f64],
    target: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    for (a, b) in source.chunks(16384).zip(target.chunks_mut(16384)) {
        poll(cancelled)?;
        b.copy_from_slice(a);
    }
    Ok(())
}

fn copy_physical(
    source: &MixtureEuler3d,
    target: &mut MixtureEuler3d,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    copy_field(&source.states, &mut target.states, cancelled)?;
    copy_field(
        &source.boundary_impulses,
        &mut target.boundary_impulses,
        cancelled,
    )?;
    target.time = source.time;
    target.last_step = source.last_step;
    Ok(())
}

struct ChemistryAccumulator {
    report: ChemistryReport,
    correction: [f64; MAX_SPECIES + 2],
}

impl Default for ChemistryAccumulator {
    fn default() -> Self {
        Self {
            report: ChemistryReport::default(),
            correction: [0.0; MAX_SPECIES + 2],
        }
    }
}

impl ChemistryAccumulator {
    fn add(&mut self, next: ChemistryReport, n: usize) -> Result<(), ReactiveError> {
        let report = &mut self.report;
        for (total, increment) in [
            (&mut report.cell_solves, next.cell_solves),
            (&mut report.accepted_steps, next.accepted_steps),
            (&mut report.rejected_steps, next.rejected_steps),
            (&mut report.nonlinear_iterations, next.nonlinear_iterations),
        ] {
            *total = total.checked_add(increment).ok_or(ReactiveError::Invalid(
                "spatial chemistry work counter overflow",
            ))?;
        }
        report.maximum_energy_residual_j = report
            .maximum_energy_residual_j
            .max(next.maximum_energy_residual_j);
        report.maximum_reconstruction_residual_j = report
            .maximum_reconstruction_residual_j
            .max(next.maximum_reconstruction_residual_j);
        for k in 0..n {
            add(
                &mut report.species_change_kg[k],
                &mut self.correction[k],
                next.species_change_kg[k],
            );
        }
        add(
            &mut report.entropy_change_j_k,
            &mut self.correction[n],
            next.entropy_change_j_k,
        );
        add(
            &mut report.referenced_energy_change_j,
            &mut self.correction[n + 1],
            next.referenced_energy_change_j,
        );
        Ok(())
    }
}

fn react_cells(
    flow: &mut MixtureEuler3d,
    reactor: &mut CellChemistry,
    duration: f64,
    report: &mut ChemistryAccumulator,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    let n = flow.gas.species_count();
    let volume = Volume::new(flow.grid.cell_volume());
    for (i, q) in flow.states.chunks_exact_mut(n + 4).enumerate() {
        poll(cancelled)?;
        let mut next = ChemistryReport::default();
        reactor.advance_spatial(
            q,
            volume,
            Time::new(duration),
            i,
            &mut next,
            &mut *cancelled,
        )?;
        report.add(next, n)?;
    }
    poll(cancelled)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReactiveError> {
    if cancelled() {
        Err(ReactiveError::Cancelled)
    } else {
        Ok(())
    }
}

fn add(total: &mut f64, correction: &mut f64, value: f64) {
    let y = value - *correction;
    let next = *total + y;
    *correction = (next - *total) - y;
    *total = next;
}

fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        add(&mut total, &mut correction, value);
    }
    total
}

#[cfg(test)]
mod tests;
