//! Allocation-free temporal error control for finite bulk and film inventories.
//!
//! A trial compares one backward-Euler step with two half steps from the SAME
//! initial inventories. Their difference estimates the fine solution's local
//! error (the first-order Richardson divisor is one for equal halves). Rounded
//! unequal halves use the corresponding first-order error factor. Only the positive,
//! conservative fine state is retained; extrapolation would lose those bounds.
//! This is a local temporal estimate, not a spatial or global error certificate.

use super::*;

/// Temporal tolerances, separate from the algebraic controls of each step.
#[derive(Clone, Copy, Debug)]
pub struct TransientAdvanceControls {
    /// Largest proposed interval; endpoint clock rounding may extend it by one
    /// ULP of the target time to avoid an unsplittable terminal remainder.
    pub maximum_step_s: f64,
    /// A final remainder may be shorter; it must still admit two distinct halves.
    pub minimum_step_s: f64,
    pub relative_tolerance: f64,
    /// Applied independently to EVERY bulk/film solute amount, including zeros.
    pub absolute_amount_tolerance_mol: f64,
    /// Applied independently to both finite bulk solvent volumes.
    pub absolute_water_tolerance_m3: f64,
    /// Also controls both reconstructed membrane surface concentrations.
    pub absolute_face_concentration_tolerance_mol_m3: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    /// Explicit caps on inner root/linear tolerances and root work. `advance`
    /// further tightens them to 1% of temporal budgets, distributing absolute
    /// allowances over the window and converting face concentration to amounts
    /// with the smallest admitted cell/bulk volume. No tolerance is relaxed.
    pub step: TransientStepControls,
}
impl Default for TransientAdvanceControls {
    fn default() -> Self {
        Self {
            maximum_step_s: 1.0,
            minimum_step_s: 1e-12,
            relative_tolerance: 1e-4,
            absolute_amount_tolerance_mol: 1e-14,
            absolute_water_tolerance_m3: 1e-15,
            absolute_face_concentration_tolerance_mol_m3: 1e-6,
            maximum_steps: 10_000,
            maximum_rejections: 128,
            step: TransientStepControls::default(),
        }
    }
}

/// Accepted window and balances computed directly from original/final state.
#[derive(Clone, Debug)]
pub struct TransientAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    /// All implicit calls, including coarse and refused trials.
    pub step_evaluations: usize,
    /// All constitutive root evaluations, including refused trials.
    pub root_evaluations: usize,
    pub maximum_scaled_error: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub solvent_transfer_m3: f64,
    pub pressure_work_into_solution_j: f64,
    pub water_balance_residual_m3: f64,
    pub solute_balance_residual_mol: [f64; MAX_OSMOTIC_SPECIES],
    pub species_count: usize,
    /// Last accepted half step: its actual endpoint flux, membrane faces and
    /// algebraic residuals. Its interval covers only that final half step.
    pub endpoint: TransientStepReport,
}

struct State {
    bulks: [OsmoticCompartment; 2],
    film_amounts: Box<[f64]>,
    kinetic: Option<kinetic::Pools>,
    time_s: f64,
}
impl State {
    fn new(pair: &TransientOsmoticPair) -> Result<Self> {
        Ok(Self {
            bulks: pair.bulks.clone(),
            film_amounts: zeros(pair.film_amounts.len())?,
            kinetic: pair
                .kinetic
                .as_ref()
                .map(|k| kinetic::Pools::new(k.state.trapped.len()))
                .transpose()?,
            time_s: pair.time_s,
        })
    }

    fn capture(
        &mut self,
        pair: &TransientOsmoticPair,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        copy_bulks(&pair.bulks, &mut self.bulks)?;
        copy_amounts(&pair.film_amounts, &mut self.film_amounts, cancelled)?;
        if let (Some(source), Some(target)) = (&pair.kinetic, &mut self.kinetic) {
            copy_amounts(
                &source.state.equilibrium,
                &mut target.equilibrium,
                cancelled,
            )?;
            copy_amounts(&source.state.trapped, &mut target.trapped, cancelled)?;
        }
        self.time_s = pair.time_s;
        Ok(())
    }

    fn restore(
        &self,
        pair: &mut TransientOsmoticPair,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        copy_bulks(&self.bulks, &mut pair.bulks)?;
        copy_amounts(&self.film_amounts, &mut pair.film_amounts, cancelled)?;
        if let (Some(source), Some(target)) = (&self.kinetic, &mut pair.kinetic) {
            copy_amounts(
                &source.equilibrium,
                &mut target.state.equilibrium,
                cancelled,
            )?;
            copy_amounts(&source.trapped, &mut target.state.trapped, cancelled)?;
        }
        pair.time_s = self.time_s;
        Ok(())
    }

    fn amount(&self, row: usize, species: usize, rows: usize, count: usize) -> f64 {
        if row == 0 {
            self.bulks[0].solution().solute_moles()[species]
        } else if row == rows - 1 {
            self.bulks[1].solution().solute_moles()[species]
        } else {
            self.film_amounts[(row - 1) * count + species]
        }
    }
    fn pools(&self, row: usize, species: usize, rows: usize, count: usize) -> [f64; 2] {
        match &self.kinetic {
            Some(k) if row > 0 && row + 1 < rows => {
                let index = (row - 1) * count + species;
                [k.equilibrium[index], k.trapped[index]]
            }
            _ => [self.amount(row, species, rows, count), 0.0],
        }
    }

    /// Infallible, allocation-free rollback even if cancellation interrupted a
    /// trial restore or occurred after the final accepted internal substep.
    fn swap_with(&mut self, pair: &mut TransientOsmoticPair) {
        core::mem::swap(&mut self.bulks, &mut pair.bulks);
        core::mem::swap(&mut self.film_amounts, &mut pair.film_amounts);
        if let (Some(snapshot), Some(current)) = (&mut self.kinetic, &mut pair.kinetic) {
            core::mem::swap(snapshot, &mut current.state);
        }
        core::mem::swap(&mut self.time_s, &mut pair.time_s);
    }
}

fn copy_bulks(from: &[OsmoticCompartment; 2], to: &mut [OsmoticCompartment; 2]) -> Result<()> {
    for (source, target) in from.iter().zip(to) {
        target.solution.replace_inventories(
            source.solution().solvent_volume_m3(),
            source.solution().solute_moles(),
        )?;
        target.pressure = source.pressure;
        target.minimum_volume_m3 = source.minimum_volume_m3;
    }
    Ok(())
}
fn copy_amounts(from: &[f64], to: &mut [f64], cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    for (source, target) in from.chunks(4096).zip(to.chunks_mut(4096)) {
        poll(cancelled)?;
        target.copy_from_slice(source);
    }
    Ok(())
}

/// Three prepared snapshots; the pair's existing root/tridiagonal storage is
/// reused. Preparation is opt-in so single-step callers pay no extra storage.
pub struct TransientAdvanceWorkspace {
    original: State,
    start: State,
    coarse: State,
    rows: usize,
    species_count: usize,
    bytes: usize,
}
impl TransientAdvanceWorkspace {
    /// Additional owned storage, including the workspace itself.
    pub const fn workspace_bytes(&self) -> usize {
        self.bytes
    }
    fn compatible(&self, pair: &TransientOsmoticPair) -> bool {
        self.rows == pair.rows
            && self.species_count == pair.species_count
            && self.original.kinetic.is_some() == pair.kinetic.is_some()
            && self.original.bulks.iter().zip(&pair.bulks).all(|(a, b)| {
                let a = a.solution();
                let b = b.solution();
                a.temperature_k().to_bits() == b.temperature_k().to_bits()
                    && a.solvent_molar_volume_m3_mol().to_bits()
                        == b.solvent_molar_volume_m3_mol().to_bits()
            })
    }
}

impl TransientAdvanceControls {
    fn validate(&self) -> Result<()> {
        if !positive(self.maximum_step_s)
            || !positive(self.minimum_step_s)
            || self.minimum_step_s > self.maximum_step_s
            || !positive(self.relative_tolerance)
            || self.relative_tolerance >= 1.0
            || !positive(self.absolute_amount_tolerance_mol)
            || !positive(self.absolute_water_tolerance_m3)
            || !positive(self.absolute_face_concentration_tolerance_mol_m3)
            || self.maximum_steps == 0
            || !positive(self.step.relative_root_tolerance)
            || self.step.relative_root_tolerance >= 1.0
            || !positive(self.step.absolute_water_tolerance_m3)
            || !positive(self.step.relative_amount_tolerance)
            || self.step.relative_amount_tolerance >= 1.0
            || !positive(self.step.absolute_amount_tolerance_mol)
            || !(2..=4096).contains(&self.step.maximum_root_evaluations)
            || !(1..=4096).contains(&self.step.maximum_sorption_iterations)
            || !(1..=128).contains(&self.step.maximum_sorption_backtracks)
        {
            return Err(invalid("invalid transient adaptive controls"));
        }
        Ok(())
    }

    fn inner(&self, pair: &TransientOsmoticPair, fraction: f64) -> Result<TransientStepControls> {
        let minimum_volume = pair
            .films
            .iter()
            .flatten()
            .flat_map(|film| film.cell_volumes_m3.iter().copied())
            .chain(pair.bulks.iter().map(|bulk| bulk.minimum_volume_m3))
            .fold(f64::INFINITY, f64::min);
        let mut inner = self.step;
        inner.relative_root_tolerance = inner
            .relative_root_tolerance
            .min(0.01 * self.relative_tolerance);
        inner.relative_amount_tolerance = inner
            .relative_amount_tolerance
            .min(0.01 * self.relative_tolerance * fraction);
        inner.absolute_water_tolerance_m3 = inner
            .absolute_water_tolerance_m3
            .min(0.01 * self.absolute_water_tolerance_m3 * fraction);
        inner.absolute_amount_tolerance_mol = inner
            .absolute_amount_tolerance_mol
            .min(0.01 * self.absolute_amount_tolerance_mol * fraction)
            .min(
                0.01 * self.absolute_face_concentration_tolerance_mol_m3
                    * minimum_volume
                    * fraction,
            );
        if !positive(inner.relative_amount_tolerance)
            || !positive(inner.absolute_amount_tolerance_mol)
            || !positive(inner.absolute_water_tolerance_m3)
        {
            return Err(invalid(
                "adaptive algebraic tolerances are not representable",
            ));
        }
        Ok(inner)
    }
}

#[derive(Default)]
pub(super) struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    pub(super) fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    pub(super) fn total(&self) -> Result<f64> {
        finite(self.value + self.correction)
    }
}

fn scaled_error(start: f64, coarse: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64> {
    let scale = start.abs().max(coarse.abs()).max(fine.abs());
    let allowance = finite(absolute + relative * scale)?;
    finite((coarse - fine).abs() / allowance)
}

fn retryable(error: &TransientOsmosisError) -> bool {
    matches!(
        error,
        TransientOsmosisError::Network(
            OsmoticNetworkError::Unbracketed
                | OsmoticNetworkError::NonlinearConvergence
                | OsmoticNetworkError::WorkBudget
        ) | TransientOsmosisError::LinearResidual { .. }
            | TransientOsmosisError::InventoryBalance { .. }
            | TransientOsmosisError::SorptionConvergence { .. }
    )
}

impl TransientOsmoticPair {
    /// Prepare optional temporal-control storage under a separate byte budget.
    /// It may be reused with a pair of the same state shape and ideal solvent
    /// model; geometry and membrane coefficients remain owned by that pair.
    pub fn adaptive_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<TransientAdvanceWorkspace> {
        let bytes = core::mem::size_of::<TransientAdvanceWorkspace>()
            + 3 * self.rows * self.species_count * core::mem::size_of::<f64>()
            + self
                .kinetic
                .as_ref()
                .map_or(0, |k| 6 * k.state.trapped.len() * size_of::<f64>());
        if bytes > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: bytes,
                budget_bytes: maximum_workspace_bytes,
            }
            .into());
        }
        Ok(TransientAdvanceWorkspace {
            original: State::new(self)?,
            start: State::new(self)?,
            coarse: State::new(self)?,
            rows: self.rows,
            species_count: self.species_count,
            bytes,
        })
    }

    /// Advance exactly to `time_s() + duration_s` using step doubling. Every
    /// amount (including mobile/linear and kinetic pools separately), both bulk volumes and both membrane faces enter the maximum
    /// scaled local-error norm. Root work, accepted steps and rejections are
    /// bounded. Any refusal/cancellation restores the complete original window.
    /// No allocation or extrapolation occurs during prepared advancement.
    pub fn advance(
        &mut self,
        duration_s: f64,
        controls: TransientAdvanceControls,
        workspace: &mut TransientAdvanceWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransientAdvanceReport> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !positive(duration_s) || !workspace.compatible(self) {
            return Err(invalid(
                "invalid adaptive duration or incompatible prepared workspace",
            ));
        }
        let target = finite(self.time_s + duration_s)?;
        if target <= self.time_s {
            return Err(OsmoticNetworkError::StepUnderflow.into());
        }
        workspace.original.capture(self, &mut cancelled)?;
        let result = self.advance_inner(target, controls, workspace, &mut cancelled);
        if result.is_err() {
            workspace.original.swap_with(self);
        }
        result
    }

    fn advance_inner(
        &mut self,
        target: f64,
        controls: TransientAdvanceControls,
        workspace: &mut TransientAdvanceWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransientAdvanceReport> {
        let start_time = self.time_s;
        let duration = target - start_time;
        let mut step = controls.maximum_step_s.min(duration);
        let mut accepted = 0;
        let mut rejected = 0;
        let mut evaluations = 0usize;
        let mut roots = 0;
        let mut maximum_error: f64 = 0.0;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step: f64 = 0.0;
        let mut pressure_work = Sum::default();
        let mut endpoint = None;
        while self.time_s < target {
            poll(cancelled)?;
            if accepted >= controls.maximum_steps {
                return Err(OsmoticNetworkError::WorkBudget.into());
            }
            workspace.start.capture(self, cancelled)?;
            let remaining = target - self.time_s;
            let planned = if step >= remaining {
                target
            } else {
                self.time_s + step
            };
            // Fold a one-ULP terminal remainder before either solve. Merely
            // changing the published clock would omit physical transport;
            // attempting a separate remainder cannot form two distinct halves.
            let end = if planned < target && planned.next_up() >= target {
                target
            } else {
                planned
            };
            let dt = end - self.time_s;
            let midpoint = self.time_s + 0.5 * dt;
            if midpoint <= self.time_s || midpoint >= end || end > target {
                return Err(OsmoticNetworkError::StepUnderflow.into());
            }
            let first_dt = midpoint - self.time_s;
            let second_dt = end - midpoint;
            // Use actual representable durations: near the clock's resolution
            // limit a midpoint can split the interval into unequal halves.
            let coarse_controls = controls.inner(self, dt / duration)?;
            let first_controls = controls.inner(self, first_dt / duration)?;
            let second_controls = controls.inner(self, second_dt / duration)?;
            let u = first_dt / dt;
            let v = second_dt / dt;
            let fine_error_factor = (u * u + v * v) / (2.0 * u * v);
            let trial = (|| -> Result<_> {
                evaluations = evaluations
                    .checked_add(1)
                    .ok_or(OsmoticNetworkError::WorkBudget)?;
                let coarse = self.step_counted(dt, coarse_controls, &mut roots, &mut *cancelled)?;
                workspace.coarse.capture(self, cancelled)?;
                workspace.start.restore(self, cancelled)?;
                evaluations = evaluations
                    .checked_add(1)
                    .ok_or(OsmoticNetworkError::WorkBudget)?;
                let first =
                    self.step_counted(first_dt, first_controls, &mut roots, &mut *cancelled)?;
                if self.time_s != midpoint {
                    return Err(OsmoticNetworkError::StepUnderflow.into());
                }
                evaluations = evaluations
                    .checked_add(1)
                    .ok_or(OsmoticNetworkError::WorkBudget)?;
                let second =
                    self.step_counted(second_dt, second_controls, &mut roots, &mut *cancelled)?;
                if self.time_s != end {
                    return Err(OsmoticNetworkError::StepUnderflow.into());
                }
                let error = self.temporal_error(
                    &workspace.start,
                    &workspace.coarse,
                    &coarse,
                    &second,
                    controls,
                    cancelled,
                )? * fine_error_factor;
                Ok((error, first.pressure_work_into_solution_j, second))
            })();
            let (error, first_work, second) = match trial {
                Ok(value) => value,
                Err(error) if retryable(&error) => {
                    workspace.start.restore(self, cancelled)?;
                    if rejected >= controls.maximum_rejections || 0.5 * dt < controls.minimum_step_s
                    {
                        return Err(error);
                    }
                    rejected += 1;
                    step = 0.5 * dt;
                    continue;
                }
                Err(error) => return Err(error),
            };
            if error > 1.0 {
                workspace.start.restore(self, cancelled)?;
                if 0.5 * dt < controls.minimum_step_s {
                    return Err(TransientOsmosisError::TemporalAccuracy {
                        step_s: dt,
                        scaled_error: error,
                    });
                }
                if rejected >= controls.maximum_rejections {
                    return Err(OsmoticNetworkError::WorkBudget.into());
                }
                rejected += 1;
                step = 0.5 * dt;
                continue;
            }
            accepted += 1;
            maximum_error = maximum_error.max(error);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            pressure_work.add(first_work);
            pressure_work.add(second.pressure_work_into_solution_j);
            endpoint = Some(second);
            step = if error < 0.125 {
                (2.0 * dt).min(controls.maximum_step_s)
            } else {
                dt
            };
        }
        let original = &workspace.original;
        let mut water_balance = Sum::default();
        for side in 0..2 {
            water_balance.add(
                self.bulks[side].solution().solvent_volume_m3()
                    - original.bulks[side].solution().solvent_volume_m3(),
            );
        }
        let water_residual = water_balance.total()?;
        if water_residual.abs()
            > controls.step.absolute_water_tolerance_m3
                + 32.0 * f64::EPSILON * self.total_solvent_volume_m3()
        {
            return Err(invalid("adaptive bulk-plus-film water balance failed"));
        }
        let mut balance = [0.0; MAX_OSMOTIC_SPECIES];
        for (species, residual) in balance.iter_mut().enumerate().take(self.species_count) {
            poll(cancelled)?;
            let mut difference = Sum::default();
            let mut before = Sum::default();
            for row in 0..self.rows {
                let pools = original.pools(row, species, self.rows, self.species_count);
                difference.add(self.equilibrium_amount(row, species) - pools[0]);
                difference.add(self.trapped_amount(row, species) - pools[1]);
                before.add(pools[0]);
                before.add(pools[1]);
            }
            *residual = difference.total()?;
            if residual.abs()
                > controls.step.absolute_amount_tolerance_mol
                    + controls.step.relative_amount_tolerance * before.total()?
            {
                return Err(TransientOsmosisError::InventoryBalance {
                    species,
                    residual_mol: *residual,
                });
            }
        }
        let report = TransientAdvanceReport {
            interval_s: [start_time, target],
            accepted_steps: accepted,
            rejected_steps: rejected,
            step_evaluations: evaluations,
            root_evaluations: roots,
            maximum_scaled_error: maximum_error,
            minimum_accepted_step_s: minimum_step,
            maximum_accepted_step_s: maximum_step,
            solvent_transfer_m3: original.bulks[0].solution().solvent_volume_m3()
                - self.bulks[0].solution().solvent_volume_m3(),
            pressure_work_into_solution_j: pressure_work.total()?,
            water_balance_residual_m3: water_residual,
            solute_balance_residual_mol: balance,
            species_count: self.species_count,
            endpoint: endpoint.ok_or(OsmoticNetworkError::StepUnderflow)?,
        };
        poll(cancelled)?;
        Ok(report)
    }

    fn temporal_error(
        &self,
        start: &State,
        coarse: &State,
        coarse_report: &TransientStepReport,
        fine_report: &TransientStepReport,
        controls: TransientAdvanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        let mut error: f64 = 0.0;
        for row in 0..self.rows {
            if row % 64 == 0 {
                poll(cancelled)?;
            }
            for species in 0..self.species_count {
                error = error.max(scaled_error(
                    start.amount(row, species, self.rows, self.species_count),
                    coarse.amount(row, species, self.rows, self.species_count),
                    self.amount(row, species),
                    controls.absolute_amount_tolerance_mol,
                    controls.relative_tolerance,
                )?);
                if self.kinetic.is_some() {
                    let start = start.pools(row, species, self.rows, self.species_count);
                    let coarse = coarse.pools(row, species, self.rows, self.species_count);
                    for (index, fine) in [
                        self.equilibrium_amount(row, species),
                        self.trapped_amount(row, species),
                    ]
                    .into_iter()
                    .enumerate()
                    {
                        error = error.max(scaled_error(
                            start[index],
                            coarse[index],
                            fine,
                            controls.absolute_amount_tolerance_mol,
                            controls.relative_tolerance,
                        )?);
                    }
                }
            }
        }
        for side in 0..2 {
            error = error.max(scaled_error(
                start.bulks[side].solution().solvent_volume_m3(),
                coarse.bulks[side].solution().solvent_volume_m3(),
                self.bulks[side].solution().solvent_volume_m3(),
                controls.absolute_water_tolerance_m3,
                controls.relative_tolerance,
            )?);
        }
        for (coarse_faces, fine_faces) in [
            (
                &coarse_report.membrane_left_concentration_mol_m3,
                &fine_report.membrane_left_concentration_mol_m3,
            ),
            (
                &coarse_report.membrane_right_concentration_mol_m3,
                &fine_report.membrane_right_concentration_mol_m3,
            ),
        ] {
            for species in 0..self.species_count {
                error = error.max(scaled_error(
                    0.0,
                    coarse_faces[species],
                    fine_faces[species],
                    controls.absolute_face_concentration_tolerance_mol_m3,
                    controls.relative_tolerance,
                )?);
            }
        }
        Ok(error)
    }
}
