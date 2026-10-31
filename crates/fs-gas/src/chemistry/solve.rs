use super::*;

#[derive(Default)]
pub(super) struct Statistics {
    pub nonlinear_iterations: usize,
    pub temperature_evaluations: usize,
    pub direct_fallbacks: usize,
}

struct StateWorkspace {
    concentrations: Vec<f64>,
    kinetics: KineticWorkspaceV1,
}

/// One allocation per bounded vector/matrix at reactor construction. Stages,
/// rejected trials, temperature roots and Newton iterations reuse that storage.
pub(super) struct Reactor {
    pub mechanism: IdealGasKineticsV1,
    pub controls: IsochoricKineticsSolve,
    volume: f64,
    target_density: f64,
    /// Untouched whole-call feed for elemental and charge inventory audits.
    pub initial: Vec<f64>,
    /// Accepted composition about which the next local reaction extents act.
    origin: Vec<f64>,
    nu: Vec<f64>,
    state: StateWorkspace,
    trial: Vec<f64>,
    residual: Vec<f64>,
    delta: Vec<f64>,
    newton: NewtonWorkspace,
    temperature_extent_derivatives: Vec<f64>,
    pub stats: Statistics,
}
impl Reactor {
    pub(super) fn workspace_bytes(&self) -> usize {
        let scalars = [
            &self.initial,
            &self.origin,
            &self.nu,
            &self.state.concentrations,
            &self.trial,
            &self.residual,
            &self.delta,
            &self.temperature_extent_derivatives,
            &self.newton.matrix,
            &self.newton.rhs,
            &self.newton.row_scale,
            &self.newton.refinement_residual,
            &self.newton.refinement_delta,
        ]
        .iter()
        .map(|values| values.capacity())
        .sum::<usize>();
        let sparse = self.newton.stoichiometric_rows.as_ref().map_or(0, |rows| {
            (rows.offsets.capacity() + rows.column_offsets.capacity())
                * core::mem::size_of::<usize>()
                + (rows.entries.capacity() + rows.column_entries.capacity())
                    * core::mem::size_of::<(usize, f64)>()
        });
        let direct = self.newton.direct_fallback.as_ref().map_or(0, |work| {
            (work.matrix.capacity() + work.rhs.capacity() + work.row_scale.capacity())
                * core::mem::size_of::<f64>()
        });
        scalars * core::mem::size_of::<f64>()
            + self.state.kinetics.workspace_bytes()
            + sparse
            + direct
    }

    pub fn new(
        mechanism: IdealGasKineticsV1,
        controls: IsochoricKineticsSolve,
    ) -> Result<Self, IsochoricKineticsError> {
        let ns = mechanism.thermodynamics().models().len();
        let nr = mechanism.reactions().len();
        let kinetics = KineticWorkspaceV1::new(&mechanism)?;
        let mut nu = zeros(ns * nr)?;
        for i in 0..ns {
            for j in 0..nr {
                nu[i * nr + j] = mechanism
                    .stoichiometry()
                    .get(i, j)
                    .expect("admitted stoichiometry") as f64;
            }
        }
        let newton = NewtonWorkspace::new(ns, nr, &nu)?;
        Ok(Self {
            mechanism,
            controls,
            volume: 1.0,
            target_density: 0.0,
            initial: zeros(ns)?,
            origin: zeros(ns)?,
            nu,
            state: StateWorkspace {
                concentrations: zeros(ns)?,
                kinetics,
            },
            trial: zeros(nr)?,
            residual: zeros(nr)?,
            delta: zeros(nr)?,
            newton,
            temperature_extent_derivatives: zeros(nr)?,
            stats: Statistics::default(),
        })
    }
    /// Replace every cell-dependent input without reallocating admitted storage.
    pub fn initialize(
        &mut self,
        amounts: &[Amount],
        volume: Volume,
        temperature: Temperature,
        target: Option<Energy>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Energy, IsochoricKineticsError> {
        if amounts.len() != self.initial.len() {
            return Err(IsochoricKineticsError::Invalid(
                "dense feed does not match canonical mechanism species",
            ));
        }
        for (concentration, amount) in self.initial.iter_mut().zip(amounts) {
            if !amount.value().is_finite() || amount.value() < 0.0 {
                return Err(IsochoricKineticsError::Invalid("invalid dense feed amount"));
            }
            *concentration = amount.value() / volume.value();
            if !concentration.is_finite() || (amount.value() > 0.0 && *concentration == 0.0) {
                return Err(IsochoricKineticsError::Invalid(
                    "feed concentration is not representable",
                ));
            }
        }
        let initial = self.mechanism.evaluate_thermodynamics(
            temperature,
            &self.initial,
            &mut self.state.kinetics,
            cancelled,
        )?;
        if !positive(initial.mass_density_kg_m3 * volume.value()) {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable original parcel mass",
            ));
        }
        let target = target
            .unwrap_or_else(|| Energy::new(initial.internal_energy_density_j_m3 * volume.value()));
        let target_density = target.value() / volume.value();
        if !target.value().is_finite()
            || !target_density.is_finite()
            || !positive(self.controls.energy_tolerance_j / volume.value())
            || !positive(self.controls.absolute_amount_tolerance_mol / volume.value())
        {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable reactor energy or concentration tolerance",
            ));
        }
        self.volume = volume.value();
        self.target_density = target_density;
        self.origin.copy_from_slice(&self.initial);
        self.stats = Statistics::default();
        Ok(target)
    }
    pub fn concentrations(&self) -> &[f64] {
        &self.state.concentrations
    }

    /// Change coordinates after an accepted step without changing any physical
    /// amount or the conserved energy. Accumulating every reaction cycle over
    /// the complete trajectory makes tiny intermediates the difference of
    /// macroscopic extents, destroying their representable positive domain.
    /// Local extents retain the same stoichiometric update and final original-
    /// feed conservation gates without clipping a species or seeding radicals.
    pub fn rebase(&mut self, concentrations: &[f64]) {
        self.origin.copy_from_slice(concentrations);
    }
    pub fn resolve(
        &mut self,
        extents: &[f64],
        hint: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<KineticThermodynamicsV1, IsochoricKineticsError> {
        resolve(
            &self.mechanism,
            self.controls,
            self.volume,
            self.target_density,
            &self.origin,
            &self.nu,
            extents,
            hint,
            &mut self.state,
            &mut self.stats,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn implicit_step(
        &mut self,
        base: &[f64],
        step: f64,
        out: &mut [f64],
        stage_forcing: &mut [f64],
        base_state: KineticThermodynamicsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<KineticThermodynamicsV1, IsochoricKineticsError> {
        out.copy_from_slice(base);
        if self.controls.method == IsochoricKineticsMethod::BackwardEuler {
            return self.implicit_stage(
                base,
                step,
                None,
                out,
                base_state.temperature.value(),
                cancelled,
            );
        }
        // Alexander SDIRK2: A = [[gamma,0],[1-gamma,gamma]], b = A[1].
        // The accepted second stage is stiffly accurate. Keep the forcing
        // separate from `base`: base+h*(1-gamma)*r(stage1) need not describe
        // any admissible composition or caloric state.
        const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let stage = self.implicit_stage(
            base,
            GAMMA * step,
            None,
            out,
            base_state.temperature.value(),
            cancelled,
        )?;
        self.stage_entropy(stage, base_state, step)?;
        // Every successful implicit_stage exit has evaluated reaction rates
        // at exactly the returned state. Reserve this sole extra vector once.
        for (forcing, rate) in stage_forcing
            .iter_mut()
            .zip(self.state.kinetics.reaction_rates())
        {
            *forcing = ((1.0 - GAMMA) * step) * rate;
            if !forcing.is_finite() {
                return Err(failed(step, "unrepresentable SDIRK stage forcing"));
            }
        }
        let endpoint = self.implicit_stage(
            base,
            GAMMA * step,
            Some(stage_forcing),
            out,
            stage.temperature.value(),
            cancelled,
        )?;
        self.stage_entropy(endpoint, base_state, step)?;
        Ok(endpoint)
    }

    fn stage_entropy(
        &self,
        stage: KineticThermodynamicsV1,
        base: KineticThermodynamicsV1,
        step: f64,
    ) -> Result<(), IsochoricKineticsError> {
        let change = (stage.entropy_density_j_m3_k - base.entropy_density_j_m3_k) * self.volume;
        if !change.is_finite() || change < -self.controls.entropy_tolerance_j_k {
            return Err(failed(
                step,
                "SDIRK stage entropy decreases from its physical baseline",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn implicit_stage(
        &mut self,
        base: &[f64],
        step: f64,
        forcing: Option<&[f64]>,
        out: &mut [f64],
        hint: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<KineticThermodynamicsV1, IsochoricKineticsError> {
        if !positive(step) {
            return Err(failed(step, "unrepresentable implicit stage duration"));
        }
        let nr = base.len();
        let ns = self.initial.len();
        let mut state = self.resolve(out, hint, cancelled)?;
        for iteration in 0..self.controls.maximum_newton_iterations {
            poll(cancelled)?;
            self.stats.nonlinear_iterations += 1;
            self.mechanism.evaluate(
                state.temperature,
                &self.state.concentrations,
                &mut self.state.kinetics,
                &mut *cancelled,
            )?;
            for j in 0..nr {
                self.residual[j] = (out[j] - base[j])
                    - forcing.map_or(0.0, |f| f[j])
                    - step * self.state.kinetics.reaction_rates()[j];
            }
            let norm = residual_norm(
                self.controls,
                self.volume,
                &self.nu,
                &self.state.concentrations,
                &self.residual,
            );
            // A nonzero source must receive at least one Newton correction.
            // Otherwise many small admitted steps could silently freeze slow
            // chemistry merely because each individual increment is small.
            if norm <= 0.02 && (iteration > 0 || self.residual.iter().all(|r| *r == 0.0)) {
                return Ok(state);
            }
            if !norm.is_finite() {
                return Err(failed(step, "non-finite implicit residual"));
            }
            let u = self.state.kinetics.species_internal_energies();
            for k in 0..nr {
                let mut du = Sum::default();
                for i in 0..ns {
                    du.add(u[i] * self.nu[i * nr + k]);
                }
                self.temperature_extent_derivatives[k] =
                    -du.value() / state.heat_capacity_v_density_j_m3_k;
            }
            let direct_fallback = self.newton.solve(
                &self.nu,
                self.state.kinetics.rate_jacobian(),
                self.state.kinetics.temperature_derivatives(),
                &self.temperature_extent_derivatives,
                state.temperature.value(),
                step,
                &self.residual,
                &mut self.delta,
                cancelled,
            )?;
            self.stats.direct_fallbacks += usize::from(direct_fallback);
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..self.controls.maximum_backtracks {
                poll(cancelled)?;
                for j in 0..nr {
                    self.trial[j] = out[j] + fraction * self.delta[j];
                }
                let candidate = resolve(
                    &self.mechanism,
                    self.controls,
                    self.volume,
                    self.target_density,
                    &self.origin,
                    &self.nu,
                    &self.trial,
                    state.temperature.value(),
                    &mut self.state,
                    &mut self.stats,
                    cancelled,
                );
                match candidate {
                    Err(IsochoricKineticsError::Cancelled) => {
                        return Err(IsochoricKineticsError::Cancelled);
                    }
                    Ok(candidate) => {
                        let rates = self.mechanism.evaluate(
                            candidate.temperature,
                            &self.state.concentrations,
                            &mut self.state.kinetics,
                            &mut *cancelled,
                        );
                        match rates {
                            Err(KineticsErrorV1::Cancelled) => {
                                return Err(IsochoricKineticsError::Cancelled);
                            }
                            Ok(_) => {
                                for j in 0..nr {
                                    self.residual[j] = (self.trial[j] - base[j])
                                        - forcing.map_or(0.0, |f| f[j])
                                        - step * self.state.kinetics.reaction_rates()[j];
                                }
                                let trial_norm = residual_norm(
                                    self.controls,
                                    self.volume,
                                    &self.nu,
                                    &self.state.concentrations,
                                    &self.residual,
                                );
                                if trial_norm.is_finite()
                                    && (trial_norm <= 0.02
                                        || trial_norm < norm * (1.0 - 1e-4 * fraction))
                                {
                                    accepted = Some((candidate, trial_norm));
                                    break;
                                }
                            }
                            Err(_) => {}
                        }
                    }
                    Err(_) => {}
                }
                fraction *= 0.5;
            }
            let Some((candidate, trial_norm)) = accepted else {
                return Err(failed(
                    step,
                    "positivity/energy Newton line search exhausted",
                ));
            };
            out.copy_from_slice(&self.trial);
            state = candidate;
            if trial_norm <= 0.02 {
                return Ok(state);
            }
        }
        Err(failed(step, "implicit Newton iteration budget exhausted"))
    }
}

/// Bounded factorization of the smaller equivalent Newton system. State is
/// always reconstructed from every original reaction extent, including cycles
/// and duplicate reaction directions; no stoichiometric rank cutoff is used.
pub(super) struct NewtonWorkspace {
    dimension: usize,
    stoichiometric_rows: Option<StoichiometricRows>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    row_scale: Vec<f64>,
    refinement_residual: Vec<f64>,
    refinement_delta: Vec<f64>,
    direct_fallback: Option<DirectNewtonWorkspace>,
}

/// Cold-prepared storage for the original extent equation. It is used only
/// when recovery from the smaller system fails its unchanged backward gate.
struct DirectNewtonWorkspace {
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    row_scale: Vec<f64>,
}

impl NewtonWorkspace {
    pub(super) fn new(ns: usize, nr: usize, nu: &[f64]) -> Result<Self, IsochoricKineticsError> {
        let mut workspace = Self::with_dimension(nr.min(ns + 1), nr)?;
        if workspace.dimension < nr {
            workspace.stoichiometric_rows = Some(StoichiometricRows::new(ns, nr, nu)?);
        }
        Ok(workspace)
    }

    fn with_dimension(dimension: usize, nr: usize) -> Result<Self, IsochoricKineticsError> {
        let refinement = if dimension < nr { nr } else { 0 };
        Ok(Self {
            dimension,
            stoichiometric_rows: None,
            matrix: zeros(dimension * dimension)?,
            rhs: zeros(dimension)?,
            row_scale: zeros(dimension)?,
            refinement_residual: zeros(refinement)?,
            refinement_delta: zeros(refinement)?,
            direct_fallback: if dimension < nr {
                Some(DirectNewtonWorkspace {
                    matrix: zeros(nr * nr)?,
                    rhs: zeros(nr)?,
                    row_scale: zeros(nr)?,
                })
            } else {
                None
            },
        })
    }

    #[cfg(test)]
    pub(super) fn full_extent(nr: usize) -> Result<Self, IsochoricKineticsError> {
        Self::with_dimension(nr, nr)
    }

    /// Dense reference to the same reduced equation, retained only for tests.
    #[cfg(test)]
    pub(super) fn dense_reduced(ns: usize, nr: usize) -> Result<Self, IsochoricKineticsError> {
        Self::with_dimension(nr.min(ns + 1), nr)
    }

    #[cfg(test)]
    pub(super) fn nonzero_stoichiometry(&self) -> usize {
        self.stoichiometric_rows
            .as_ref()
            .map_or(0, |rows| rows.entries.len())
    }

    #[cfg(test)]
    pub(super) fn dimension(&self) -> usize {
        self.dimension
    }

    #[allow(clippy::too_many_arguments)]
    /// Return whether a direct fallback was needed after reduced recovery.
    pub(super) fn solve(
        &mut self,
        nu: &[f64],
        concentration_derivatives: &[f64],
        temperature_derivatives: &[f64],
        temperature_extent_derivatives: &[f64],
        temperature_scale: f64,
        step: f64,
        residual: &[f64],
        delta: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, IsochoricKineticsError> {
        let equation = NewtonEquation {
            nu,
            stoichiometric_rows: self.stoichiometric_rows.as_ref(),
            concentration_derivatives,
            temperature_derivatives,
            temperature_extent_derivatives,
            temperature_scale,
            step,
            residual,
        };
        if self.dimension == residual.len() {
            equation.direct(
                &mut self.matrix,
                &mut self.rhs,
                &mut self.row_scale,
                delta,
                cancelled,
            )?;
            return Ok(false);
        }
        equation.reduced(
            &mut self.matrix,
            &mut self.rhs,
            &mut self.row_scale,
            delta,
            cancelled,
        )?;
        // Recovering extent corrections as b + h U y can cancel large stiff
        // terms. Check the original extent equation, with at most four residual
        // corrections. This uses the same prepared storage and never accepts a
        // reduced solve solely because its smaller matrix was nonsingular.
        for refinement in 0..=4 {
            if equation.check_correction(
                delta,
                &mut self.rhs,
                &mut self.row_scale,
                &mut self.refinement_residual,
                cancelled,
            )? {
                return Ok(false);
            }
            if refinement == 4 {
                break;
            }
            NewtonEquation {
                residual: &self.refinement_residual,
                ..equation
            }
            .reduced(
                &mut self.matrix,
                &mut self.rhs,
                &mut self.row_scale,
                &mut self.refinement_delta,
                cancelled,
            )?;
            for (value, correction) in delta.iter_mut().zip(&self.refinement_delta) {
                *value += correction;
            }
        }
        // Tiny source channels may require a componentwise correction that
        // the Schur recovery cannot represent accurately, even with a modest
        // matrix condition number. Reuse the original scaled-partial-pivot
        // solver on I-hUV, then apply EXACTLY the same unscaled residual gate.
        // No extent is clipped and no physical or numerical tolerance changes.
        let fallback = self
            .direct_fallback
            .as_mut()
            .expect("prepared reduced fallback");
        equation.direct(
            &mut fallback.matrix,
            &mut fallback.rhs,
            &mut fallback.row_scale,
            delta,
            cancelled,
        )?;
        for refinement in 0..=2 {
            if equation.check_correction(
                delta,
                &mut self.rhs,
                &mut self.row_scale,
                &mut self.refinement_residual,
                cancelled,
            )? {
                return Ok(true);
            }
            if refinement == 2 {
                break;
            }
            NewtonEquation {
                residual: &self.refinement_residual,
                ..equation
            }
            .direct(
                &mut fallback.matrix,
                &mut fallback.rhs,
                &mut fallback.row_scale,
                &mut self.refinement_delta,
                cancelled,
            )?;
            for (value, correction) in delta.iter_mut().zip(&self.refinement_delta) {
                *value += correction;
            }
        }
        Err(failed(
            step,
            "direct Newton fallback exceeds the original linear residual gate",
        ))
    }
}

/// Prepared exact nonzero structure of N. Entries retain increasing source
/// reaction order, so every compensated reduction has the same term order as
/// the dense equation. There is no numerical threshold or rank approximation.
struct StoichiometricRows {
    offsets: Vec<usize>,
    entries: Vec<(usize, f64)>,
    column_offsets: Vec<usize>,
    column_entries: Vec<(usize, f64)>,
}

impl StoichiometricRows {
    fn new(ns: usize, nr: usize, nu: &[f64]) -> Result<Self, IsochoricKineticsError> {
        if nr == 0 || ns.checked_mul(nr) != Some(nu.len()) || nu.iter().any(|v| !v.is_finite()) {
            return Err(IsochoricKineticsError::Invalid(
                "invalid prepared stoichiometric matrix",
            ));
        }
        let mut offsets = Vec::new();
        offsets
            .try_reserve_exact(ns + 1)
            .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(nu.iter().filter(|&&v| v != 0.0).count())
            .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
        offsets.push(0);
        for row in nu.chunks_exact(nr) {
            for (j, &coefficient) in row.iter().enumerate() {
                if coefficient != 0.0 {
                    entries.push((j, coefficient));
                }
            }
            offsets.push(entries.len());
        }
        // The uncommon direct fallback uses the same exact sparsity in column
        // order. Its sums retain increasing species order from the dense owner.
        let mut column_offsets = Vec::new();
        column_offsets
            .try_reserve_exact(nr + 1)
            .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
        let mut column_entries = Vec::new();
        column_entries
            .try_reserve_exact(entries.len())
            .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
        column_offsets.push(0);
        for j in 0..nr {
            for i in 0..ns {
                let coefficient = nu[i * nr + j];
                if coefficient != 0.0 {
                    column_entries.push((i, coefficient));
                }
            }
            column_offsets.push(column_entries.len());
        }
        Ok(Self {
            offsets,
            entries,
            column_offsets,
            column_entries,
        })
    }

    fn row(&self, species: usize) -> &[(usize, f64)] {
        &self.entries[self.offsets[species]..self.offsets[species + 1]]
    }

    fn column(&self, reaction: usize) -> &[(usize, f64)] {
        &self.column_entries[self.column_offsets[reaction]..self.column_offsets[reaction + 1]]
    }
}

#[derive(Clone, Copy)]
struct NewtonEquation<'a> {
    nu: &'a [f64],
    stoichiometric_rows: Option<&'a StoichiometricRows>,
    concentration_derivatives: &'a [f64],
    temperature_derivatives: &'a [f64],
    temperature_extent_derivatives: &'a [f64],
    temperature_scale: f64,
    step: f64,
    residual: &'a [f64],
}

impl NewtonEquation<'_> {
    fn direct(
        self,
        matrix: &mut [f64],
        rhs: &mut [f64],
        row_scale: &mut [f64],
        delta: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IsochoricKineticsError> {
        let nr = self.residual.len();
        let ns = self.nu.len() / nr;
        for j in 0..nr {
            poll(cancelled)?;
            for k in 0..nr {
                let mut derivative = Sum::default();
                if let Some(sparse) = self.stoichiometric_rows {
                    for &(i, coefficient) in sparse.column(k) {
                        derivative.add(self.concentration_derivatives[j * ns + i] * coefficient);
                    }
                } else {
                    for i in 0..ns {
                        derivative
                            .add(self.concentration_derivatives[j * ns + i] * self.nu[i * nr + k]);
                    }
                }
                derivative
                    .add(self.temperature_derivatives[j] * self.temperature_extent_derivatives[k]);
                matrix[j * nr + k] =
                    if j == k { 1.0 } else { 0.0 } - self.step * derivative.value();
            }
            rhs[j] = -self.residual[j];
        }
        solve_linear(matrix, rhs, row_scale, nr, cancelled)?;
        delta.copy_from_slice(rhs);
        Ok(())
    }

    fn reduced(
        self,
        matrix: &mut [f64],
        rhs: &mut [f64],
        row_scale: &mut [f64],
        delta: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IsochoricKineticsError> {
        let nr = self.residual.len();
        let ns = self.nu.len() / nr;
        let n = ns + 1;
        // J = I - h U V, U = [dr/dc, T dr/dT],
        // V = [N; (dT/dxi)/T]. Solve (I - h V U)y = -V F,
        // then delta_xi = -F + h U y. The final coordinate is scaled
        // temperature, not an independently advanced state. Keeping dT/dxi
        // from the original compensated stoichiometric energy sum avoids
        // differentiating arbitrary elemental reference energies off the
        // conserved composition manifold. No inverse or rank test on N occurs.
        for i in 0..ns {
            poll(cancelled)?;
            let row = &mut matrix[i * n..(i + 1) * n];
            rhs[i] = if let Some(sparse) = self.stoichiometric_rows {
                self.assemble_row(i, sparse.row(i).iter().copied(), row)
            } else {
                self.assemble_row(
                    i,
                    self.nu[i * nr..(i + 1) * nr].iter().copied().enumerate(),
                    row,
                )
            };
        }
        poll(cancelled)?;
        // The thermal row remains dense: dropping its small values would alter
        // the original energy-constrained Newton equation and datum invariance.
        rhs[ns] = self.assemble_row(
            ns,
            self.temperature_extent_derivatives
                .iter()
                .map(|v| v / self.temperature_scale)
                .enumerate(),
            &mut matrix[ns * n..],
        );
        solve_linear(matrix, rhs, row_scale, n, cancelled)?;
        for j in 0..nr {
            poll(cancelled)?;
            let mut derivative = Sum::default();
            for i in 0..ns {
                derivative.add(self.concentration_derivatives[j * ns + i] * rhs[i]);
            }
            derivative.add(self.temperature_derivatives[j] * (self.temperature_scale * rhs[ns]));
            let mut correction = Sum::default();
            correction.add(-self.residual[j]);
            correction.add(self.step * derivative.value());
            delta[j] = correction.value();
        }
        Ok(())
    }

    fn assemble_row(
        self,
        row: usize,
        terms: impl Iterator<Item = (usize, f64)> + Clone,
        matrix_row: &mut [f64],
    ) -> f64 {
        let ns = matrix_row.len() - 1;
        let mut projected = Sum::default();
        for (j, coefficient) in terms.clone() {
            projected.add(-coefficient * self.residual[j]);
        }
        for (column, value) in matrix_row.iter_mut().enumerate() {
            let mut derivative = Sum::default();
            for (j, coefficient) in terms.clone() {
                let u = if column < ns {
                    self.concentration_derivatives[j * ns + column]
                } else {
                    self.temperature_derivatives[j] * self.temperature_scale
                };
                derivative.add(coefficient * u);
            }
            *value = if row == column { 1.0 } else { 0.0 } - self.step * derivative.value();
        }
        projected.value()
    }

    /// Componentwise backward residual in the original unscaled extent
    /// equation. This is a linear arithmetic gate, not a nonlinear error bound.
    fn check_correction(
        self,
        delta: &[f64],
        projected: &mut [f64],
        projected_magnitude: &mut [f64],
        defect: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, IsochoricKineticsError> {
        let nr = self.residual.len();
        let ns = self.nu.len() / nr;
        for i in 0..ns {
            poll(cancelled)?;
            (projected[i], projected_magnitude[i]) = if let Some(sparse) = self.stoichiometric_rows
            {
                Self::project_delta(sparse.row(i).iter().copied(), delta)
            } else {
                Self::project_delta(
                    self.nu[i * nr..(i + 1) * nr].iter().copied().enumerate(),
                    delta,
                )
            };
        }
        poll(cancelled)?;
        (projected[ns], projected_magnitude[ns]) = Self::project_delta(
            self.temperature_extent_derivatives
                .iter()
                .copied()
                .enumerate(),
            delta,
        );
        let tolerance = 64.0 * f64::EPSILON * (ns + nr + 1) as f64;
        let mut closed = true;
        for j in 0..nr {
            poll(cancelled)?;
            let mut derivative = Sum::default();
            let mut magnitude = Sum::default();
            for i in 0..ns {
                let term = self.concentration_derivatives[j * ns + i] * projected[i];
                derivative.add(term);
                // Retain |N| |delta|, not |N delta|: reaction cycles may
                // cancel in species space while their rounding scale remains
                // large in the original extent equation.
                magnitude
                    .add(self.concentration_derivatives[j * ns + i].abs() * projected_magnitude[i]);
            }
            let thermal = self.temperature_derivatives[j] * projected[ns];
            derivative.add(thermal);
            magnitude.add(self.temperature_derivatives[j].abs() * projected_magnitude[ns]);
            let mut residual = Sum::default();
            residual.add(self.residual[j]);
            residual.add(delta[j]);
            residual.add(-self.step * derivative.value());
            defect[j] = residual.value();
            let scale = self.residual[j].abs() + delta[j].abs() + self.step * magnitude.value();
            if !defect[j].is_finite() || !scale.is_finite() {
                return Err(failed(self.step, "non-finite reduced Newton correction"));
            }
            closed &= defect[j].abs() <= tolerance * scale;
        }
        Ok(closed)
    }

    fn project_delta(terms: impl Iterator<Item = (usize, f64)>, delta: &[f64]) -> (f64, f64) {
        let mut value = Sum::default();
        let mut magnitude = Sum::default();
        for (j, coefficient) in terms {
            let term = coefficient * delta[j];
            value.add(term);
            magnitude.add(term.abs());
        }
        (value.value(), magnitude.value())
    }
}

fn failed(step: f64, reason: &'static str) -> IsochoricKineticsError {
    IsochoricKineticsError::StepFailure {
        elapsed_time_s: 0.0,
        step_s: step,
        reason,
    }
}

fn residual_norm(
    controls: IsochoricKineticsSolve,
    volume: f64,
    nu: &[f64],
    c: &[f64],
    residual: &[f64],
) -> f64 {
    let nr = residual.len();
    let absolute = controls.absolute_amount_tolerance_mol / volume;
    let mut norm = 0.0_f64;
    for j in 0..nr {
        if !residual[j].is_finite() {
            return f64::INFINITY;
        }
        let mut scale = f64::INFINITY;
        for i in 0..c.len() {
            let coefficient = nu[i * nr + j].abs();
            if coefficient > 0.0 {
                scale = scale.min((absolute + controls.relative_tolerance * c[i]) / coefficient);
            }
        }
        let ratio = residual[j].abs() / (scale / nr as f64);
        if !positive(scale) || !ratio.is_finite() {
            return f64::INFINITY;
        }
        norm = norm.max(ratio);
    }
    norm
}

#[allow(clippy::too_many_arguments)]
fn resolve(
    mechanism: &IdealGasKineticsV1,
    controls: IsochoricKineticsSolve,
    volume: f64,
    target_density: f64,
    initial: &[f64],
    nu: &[f64],
    extents: &[f64],
    hint: f64,
    work: &mut StateWorkspace,
    statistics: &mut Statistics,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<KineticThermodynamicsV1, IsochoricKineticsError> {
    poll(cancelled)?;
    let nr = extents.len();
    for (i, &value) in initial.iter().enumerate() {
        let mut c = Sum::default();
        c.add(value);
        for j in 0..nr {
            c.add(nu[i * nr + j] * extents[j]);
        }
        let c = c.value();
        if !c.is_finite() || c < 0.0 {
            return Err(failed(
                0.0,
                "reaction extent leaves the nonnegative species domain",
            ));
        }
        work.concentrations[i] = c;
    }
    let tolerance = controls.energy_tolerance_j / volume;
    let evaluate = |temperature: f64,
                    work: &mut StateWorkspace,
                    statistics: &mut Statistics,
                    cancelled: &mut _| {
        statistics.temperature_evaluations += 1;
        mechanism
            .evaluate_thermodynamics(
                Temperature::new(temperature),
                &work.concentrations,
                &mut work.kinetics,
                cancelled,
            )
            .map_err(IsochoricKineticsError::from)
    };
    let mut lo = controls.minimum_temperature_k;
    let mut hi = controls.maximum_temperature_k;
    let mut t = hint.clamp(lo, hi);
    let first = evaluate(t, work, statistics, &mut *cancelled)?;
    let mut residual = first.internal_energy_density_j_m3 - target_density;
    if residual.abs() <= tolerance {
        return Ok(first);
    }
    let lower = evaluate(lo, work, statistics, &mut *cancelled)?;
    let upper = evaluate(hi, work, statistics, &mut *cancelled)?;
    if (lower.internal_energy_density_j_m3 - target_density).abs() <= tolerance {
        return evaluate(lo, work, statistics, &mut *cancelled);
    }
    if (upper.internal_energy_density_j_m3 - target_density).abs() <= tolerance {
        return Ok(upper);
    }
    if !(lower.internal_energy_density_j_m3 < target_density
        && target_density < upper.internal_energy_density_j_m3)
    {
        return Err(IsochoricKineticsError::EnergyNotBracketed {
            lower_energy_j: lower.internal_energy_density_j_m3 * volume,
            upper_energy_j: upper.internal_energy_density_j_m3 * volume,
        });
    }
    let mut state = first;
    for iteration in 0..controls.maximum_temperature_iterations {
        poll(cancelled)?;
        if residual < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        let next = t - residual / state.heat_capacity_v_density_j_m3_k;
        t = if next > lo && next < hi {
            next
        } else {
            lo + 0.5 * (hi - lo)
        };
        state = evaluate(t, work, statistics, &mut *cancelled)?;
        residual = state.internal_energy_density_j_m3 - target_density;
        if residual.abs() <= tolerance {
            return Ok(state);
        }
        if t == lo || t == hi {
            return Err(IsochoricKineticsError::TemperatureNoConvergence {
                iterations: iteration + 1,
                residual_j: residual * volume,
            });
        }
    }
    Err(IsochoricKineticsError::TemperatureNoConvergence {
        iterations: controls.maximum_temperature_iterations,
        residual_j: residual * volume,
    })
}

// Purpose-specific bounded dense Newton solve (at most512 reaction unknowns).
// Scaled partial pivoting is deterministic; it refuses unresolved pivots and
// nonfinite arithmetic. The nonlinear residual remains the acceptance owner.
fn solve_linear(
    a: &mut [f64],
    b: &mut [f64],
    scale: &mut [f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), IsochoricKineticsError> {
    for i in 0..n {
        scale[i] = a[i * n..(i + 1) * n]
            .iter()
            .fold(0.0_f64, |s, x| s.max(x.abs()));
        if !positive(scale[i])
            || !b[i].is_finite()
            || a[i * n..(i + 1) * n].iter().any(|x| !x.is_finite())
        {
            return Err(failed(0.0, "non-finite Newton linear system"));
        }
    }
    for k in 0..n {
        poll(cancelled)?;
        let mut pivot = k;
        let mut best = a[k * n + k].abs() / scale[k];
        for i in k + 1..n {
            let value = a[i * n + k].abs() / scale[i];
            if value > best {
                best = value;
                pivot = i;
            }
        }
        if !best.is_finite() || best <= f64::EPSILON * n as f64 {
            return Err(failed(0.0, "numerically singular Newton matrix"));
        }
        if pivot != k {
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            b.swap(k, pivot);
            scale.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.0;
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    for i in (0..n).rev() {
        let mut value = Sum::default();
        value.add(b[i]);
        for j in i + 1..n {
            value.add(-a[i * n + j] * b[j]);
        }
        b[i] = value.value() / a[i * n + i];
        if !b[i].is_finite() {
            return Err(failed(0.0, "non-finite Newton correction"));
        }
    }
    Ok(())
}
