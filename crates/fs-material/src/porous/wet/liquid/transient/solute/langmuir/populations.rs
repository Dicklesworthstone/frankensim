//! Sparse independent Langmuir populations sharing the original mobile solve.
//!
//! Every occupied/vacant stock and local source remains separate. The existing
//! site response eliminates each LDF loading; only its amount change and
//! nonnegative concentration tangent are summed into the aqueous equation.
//! This supplies heterogeneous sorption memory, not competitive shared sites.

use super::*;
use crate::aqueous_adsorption::kinetics::MAX_POPULATIONS;

/// One independent physical sorbent population in a cell/constituent slot.
/// Identity is unique within that slot; the same name may recur in other cells.
/// Distinct populations must describe disjoint dry-sorbent/site inventories.
/// All species slots remain independent, with no inference of competition.
#[derive(Clone, Copy, Debug)]
pub struct LangmuirSolutePopulation<'a> {
    pub cell: usize,
    pub species: usize,
    pub id: &'a str,
    pub sorption: LangmuirSoluteSorption<'a>,
    pub sorbed_amount_mol_m2: f64,
    pub vacant_amount_mol_m2: f64,
}

/// The original nonlinear controls plus an explicit total population work cap.
/// At most the original batch owner's 128 populations occupy any one slot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiLangmuirSoluteControls<const N: usize> {
    pub langmuir: LangmuirSoluteControls<N>,
    pub maximum_populations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MultiLangmuirSoluteError {
    Transport(DiluteSoluteError),
    /// Index into the original population slice, including its retained ID.
    Population {
        population: usize,
        error: LangmuirSoluteError,
    },
    PopulationWork {
        required_populations: usize,
        maximum_populations: usize,
    },
}
impl From<DiluteSoluteError> for MultiLangmuirSoluteError {
    fn from(value: DiluteSoluteError) -> Self {
        Self::Transport(value)
    }
}
impl core::fmt::Display for MultiLangmuirSoluteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "porous multiple-population Langmuir transport: {self:?}")
    }
}
impl std::error::Error for MultiLangmuirSoluteError {}

/// Accepted stocks and represented exchange for ONE original population.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LangmuirPopulationState {
    pub sorbed_amount_mol_m2: f64,
    pub vacant_amount_mol_m2: f64,
    pub exchange_into_sorbed_mol_m2: f64,
    pub equation_defect_mol_m2: f64,
    pub site_balance_defect_mol_m2: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MultiLangmuirSoluteStep<const N: usize> {
    /// The original dissolved and summed-bound accounts and physical face fluxes.
    pub sorption: SorbingDiluteSoluteStep<N>,
    /// Same order as the input population slice, preserving each source identity.
    /// Use these separate stocks, not a split of the aggregate, in the next step.
    pub populations: Vec<LangmuirPopulationState>,
    pub maximum_site_balance_defect_mol_m2: [f64; N],
    pub source_evaluations: usize,
}

pub(in super::super) struct MultiLangmuirStorage<'a, const N: usize> {
    bindings: &'a [LangmuirSolutePopulation<'a>],
    offsets: Vec<usize>,
    indices: Vec<usize>,
    states: Vec<LangmuirPopulationState>,
    original: Vec<[f64; N]>,
    bound: Vec<[f64; N]>,
    exchange: Vec<[f64; N]>,
    initial_concentration: Vec<[f64; N]>,
    final_water: Vec<f64>,
    limits: Vec<[[f64; 2]; N]>,
    active: [bool; N],
    area: f64,
    duration: f64,
    controls: MultiLangmuirSoluteControls<N>,
    source_evaluations: usize,
    linear_solves: usize,
    source_error: Option<MultiLangmuirSoluteError>,
}

impl<'a, const N: usize> MultiLangmuirStorage<'a, N> {
    pub(in super::super) fn new(
        n: usize,
        duration: f64,
        bindings: &'a [LangmuirSolutePopulation<'a>],
        area: f64,
        controls: MultiLangmuirSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<Self, MultiLangmuirSoluteError> {
        poll(cancelled)?;
        let ctl = controls.langmuir;
        let (slots, base) = transport::admit_extent(n, ctl.transport)?;
        if !positive(area)
            || !positive(duration)
            || !(1..=256).contains(&ctl.maximum_iterations)
            || !(1..=64).contains(&ctl.maximum_line_searches)
            || !(1..=8192).contains(&ctl.maximum_linear_solves)
            || !(1..=10_000_000).contains(&ctl.maximum_source_evaluations)
            || controls.maximum_populations == 0
            || controls.maximum_populations > slots * MAX_POPULATIONS
        {
            return Err(DiluteSoluteError::Invalid(
                "positive physical area/time and bounded population/nonlinear work required",
            )
            .into());
        }
        if bindings.len() > controls.maximum_populations {
            return Err(MultiLangmuirSoluteError::PopulationWork {
                required_populations: bindings.len(),
                maximum_populations: controls.maximum_populations,
            });
        }
        let bytes = base + Self::extra_bytes(n, bindings.len());
        if bytes > ctl.transport.maximum_storage_bytes {
            return Err(DiluteSoluteError::Workspace {
                required_bytes: bytes,
                maximum_bytes: ctl.transport.maximum_storage_bytes,
            }
            .into());
        }
        let mut offsets = vec![0; slots + 1];
        let mut limits = vec![[[0.0, f64::MAX]; N]; n];
        let mut elements = [None; N];
        let mut active = [false; N];
        for (p, population) in bindings.iter().enumerate() {
            if p % 128 == 0 {
                poll(cancelled)?;
            }
            let (i, k, b) = (population.cell, population.species, population.sorption);
            if i >= n
                || k >= N
                || population.id.is_empty()
                || population.id.len() > 256
                || population.id.chars().any(char::is_control)
            {
                return Err(Self::invalid_population(
                    p,
                    "valid cell/species and printable population identity required",
                ));
            }
            let old = population.sorbed_amount_mol_m2;
            let vacant = population.vacant_amount_mol_m2;
            if !nonnegative(old)
                || !nonnegative(vacant)
                || !positive(b.dry_mass_kg_m2)
                || !positive(b.temperature.value())
                || b.law.limits()[0] > b.law.limits()[1]
                || [
                    CONCENTRATION_AXIS,
                    TEMPERATURE_AXIS,
                    LOADING_AXIS,
                    VOLUME_AXIS,
                ]
                .iter()
                .any(|axis| b.conditions.axes().contains_key(*axis))
            {
                return Err(Self::invalid_population(
                    p,
                    "nonnegative stocks, positive dry mass/temperature and original source coordinates required",
                ));
            }
            let capacity =
                integrated_flux(b.dry_mass_kg_m2, b.law.equilibrium().capacity().value())?;
            if old > capacity
                || vacant > capacity
                || site_balance(old, vacant, capacity)?.abs()
                    > ctl.transport.absolute_equation_tolerance_mol_m2[k]
            {
                return Err(Self::invalid_population(
                    p,
                    "occupied plus vacant stocks must equal this population's original finite capacity",
                ));
            }
            let element = b.law.equilibrium().element();
            if elements[k].is_some_and(|e| e != element) {
                return Err(Self::invalid_population(
                    p,
                    "all populations and cells in one constituent must name the same analytical element",
                ));
            }
            elements[k] = Some(element);
            let support = b.law.limits();
            limits[i][k][0] = limits[i][k][0].max(support[0]);
            limits[i][k][1] = limits[i][k][1].min(support[1]);
            if limits[i][k][0] > limits[i][k][1] {
                return Err(Self::invalid_population(
                    p,
                    "coexisting populations need intersecting original concentration support",
                ));
            }
            active[k] |= b.law.rate() != Some(0.0);
            let slot = i * N + k;
            offsets[slot + 1] += 1;
            if offsets[slot + 1] > MAX_POPULATIONS {
                return Err(Self::invalid_population(
                    p,
                    "at most 128 independent populations per cell/constituent required",
                ));
            }
        }
        for k in 0..N {
            if elements[k].is_some() && elements[..k].contains(&elements[k]) {
                return Err(DiluteSoluteError::Invalid(
                    "duplicate analytical element in independent constituent slots",
                )
                .into());
            }
        }
        for slot in 0..slots {
            offsets[slot + 1] += offsets[slot];
        }
        let mut indices = vec![0; bindings.len()];
        // Reuse the starts as scatter cursors, then shift the resulting ends.
        // The original input order is retained within every sparse slot.
        for (p, population) in bindings.iter().enumerate() {
            if p % 128 == 0 {
                poll(cancelled)?;
            }
            let slot = population.cell * N + population.species;
            indices[offsets[slot]] = p;
            offsets[slot] += 1;
        }
        for slot in (1..=slots).rev() {
            offsets[slot] = offsets[slot - 1];
        }
        offsets[0] = 0;
        let mut original = vec![[0.0; N]; n];
        for i in 0..n {
            poll(cancelled)?;
            let cell_first = indices[offsets[i * N]..offsets[(i + 1) * N]]
                .first()
                .copied();
            for k in 0..N {
                let slot = i * N + k;
                let group = &indices[offsets[slot]..offsets[slot + 1]];
                for (j, &p) in group.iter().enumerate() {
                    if group[..j].iter().any(|&q| bindings[p].id == bindings[q].id) {
                        return Err(Self::invalid_population(
                            p,
                            "duplicate population identity in one cell/constituent",
                        ));
                    }
                    if cell_first.is_some_and(|q| {
                        bindings[p].sorption.temperature != bindings[q].sorption.temperature
                    }) {
                        return Err(Self::invalid_population(
                            p,
                            "all populations and constituents in a cell must use the same prescribed liquid temperature",
                        ));
                    }
                }
                original[i][k] =
                    compensated_sum(group.iter().map(|&p| bindings[p].sorbed_amount_mol_m2))?;
            }
        }
        Ok(Self {
            bindings,
            offsets,
            indices,
            states: vec![LangmuirPopulationState::default(); bindings.len()],
            original,
            bound: vec![[0.0; N]; n],
            exchange: vec![[0.0; N]; n],
            initial_concentration: vec![[0.0; N]; n],
            final_water: vec![0.0; n],
            limits,
            active,
            area,
            duration,
            controls,
            source_evaluations: 0,
            linear_solves: 0,
            source_error: None,
        })
    }

    fn extra_bytes(n: usize, populations: usize) -> usize {
        // Six cell/species fields (including two source limits), final water,
        // per-species guess/limits and the original eight Newton arrays, CSR,
        // each output state, stack tangents, and overlapping phase reports.
        (6 * n * N + 12 * n + 24 * N + MAX_POPULATIONS) * core::mem::size_of::<f64>()
            + (n * N + 1 + populations) * core::mem::size_of::<usize>()
            + populations * core::mem::size_of::<LangmuirPopulationState>()
            + N * core::mem::size_of::<bool>()
            + 1536
    }

    fn invalid_population(p: usize, message: &'static str) -> MultiLangmuirSoluteError {
        MultiLangmuirSoluteError::Population {
            population: p,
            error: DiluteSoluteError::Invalid(message).into(),
        }
    }

    fn refuse(&mut self, p: usize, error: LangmuirSoluteError) -> DiluteSoluteError {
        self.source_error = Some(MultiLangmuirSoluteError::Population {
            population: p,
            error,
        });
        DiluteSoluteError::Invalid(
            "original population source, equation or evaluation budget refused",
        )
    }

    fn tick(&mut self, p: usize) -> Result<()> {
        let maximum = self.controls.langmuir.maximum_source_evaluations;
        if self.source_evaluations >= maximum {
            return Err(self.refuse(
                p,
                LangmuirSoluteError::SourceBudget {
                    evaluations: self.source_evaluations,
                    maximum_evaluations: maximum,
                },
            ));
        }
        self.source_evaluations += 1;
        Ok(())
    }

    fn response(&mut self, i: usize, k: usize, c: f64) -> Result<(f64, f64)> {
        let slot = i * N + k;
        let (start, end) = (self.offsets[slot], self.offsets[slot + 1]);
        let mut tangents = [0.0; MAX_POPULATIONS];
        for j in start..end {
            let p = self.indices[j];
            let population = self.bindings[p];
            self.tick(p)?;
            let response = site_response(
                population.sorption,
                i,
                k,
                population.sorbed_amount_mol_m2,
                population.vacant_amount_mol_m2,
                c,
                self.final_water[i],
                self.area,
                self.duration,
                self.controls.langmuir.transport,
            )
            .map_err(|error| self.refuse(p, error))?;
            self.states[p] = LangmuirPopulationState {
                sorbed_amount_mol_m2: response.bound,
                vacant_amount_mol_m2: response.vacant,
                exchange_into_sorbed_mol_m2: response.exchange,
                equation_defect_mol_m2: response.equation_defect,
                site_balance_defect_mol_m2: response.site_defect,
            };
            tangents[j - start] = response.tangent;
        }
        let group = &self.indices[start..end];
        self.bound[i][k] =
            compensated_sum(group.iter().map(|&p| self.states[p].sorbed_amount_mol_m2))?;
        self.exchange[i][k] = compensated_sum(
            group
                .iter()
                .map(|&p| self.states[p].exchange_into_sorbed_mol_m2),
        )?;
        Ok((
            self.exchange[i][k],
            compensated_sum(tangents[..end - start].iter().copied())?,
        ))
    }

    pub(in super::super) fn recover<T>(
        &mut self,
        result: Result<T>,
    ) -> core::result::Result<T, MultiLangmuirSoluteError> {
        result.map_err(|error| self.source_error.take().unwrap_or(error.into()))
    }

    pub(in super::super) fn finish(
        self,
        initial: &[[f64; N]],
        step: DiluteSoluteStep<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<MultiLangmuirSoluteStep<N>, MultiLangmuirSoluteError> {
        let mut maximum_defect = [0.0_f64; N];
        let mut maximum_site_defect = [0.0_f64; N];
        for (p, state) in self.states.iter().enumerate() {
            if p % 128 == 0 {
                poll(cancelled)?;
            }
            let k = self.bindings[p].species;
            maximum_defect[k] = maximum_defect[k].max(state.equation_defect_mol_m2.abs());
            maximum_site_defect[k] =
                maximum_site_defect[k].max(state.site_balance_defect_mol_m2.abs());
        }
        let sorption = finish_sorption(
            &self.original,
            self.bound,
            self.exchange,
            maximum_defect,
            initial,
            step,
            cancelled,
        )?;
        Ok(MultiLangmuirSoluteStep {
            sorption,
            populations: self.states,
            maximum_site_balance_defect_mol_m2: maximum_site_defect,
            source_evaluations: self.source_evaluations,
        })
    }
}

impl<const N: usize> ImmobileStorage<N> for MultiLangmuirStorage<'_, N> {
    fn additional_bytes(&self) -> usize {
        Self::extra_bytes(self.original.len(), self.bindings.len())
    }

    fn admit(
        &mut self,
        mobile: &[[f64; N]],
        before: &[f64],
        after: &[f64],
        controls: DiluteSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        for i in 0..mobile.len() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            self.final_water[i] = after[i];
            for k in 0..N {
                let c = finite(mobile[i][k] / before[i])?;
                if mobile[i][k] > 0.0 && c == 0.0 {
                    return Err(DiluteSoluteError::Invalid(
                        "positive initial concentration is unrepresentable",
                    ));
                }
                self.initial_concentration[i][k] = c;
            }
        }
        for p in 0..self.bindings.len() {
            if p % 128 == 0 {
                poll(cancelled)?;
            }
            let population = self.bindings[p];
            let (i, k, b) = (population.cell, population.species, population.sorption);
            let c = self.initial_concentration[i][k];
            self.tick(p)?;
            admit_site_source(
                b,
                i,
                k,
                c,
                before[i],
                self.area,
                population.sorbed_amount_mol_m2,
            )
            .map_err(|error| self.refuse(p, error))?;
            if matches!(b.law, LangmuirSoluteLaw::Equilibrium(_)) {
                let (q, vacant_q, _) =
                    b.law.equilibrium().loading_components(c).map_err(|error| {
                        self.refuse(
                            p,
                            LangmuirSoluteError::EquilibriumSource {
                                cell: i,
                                species: k,
                                error,
                            },
                        )
                    })?;
                for residual in [
                    population.sorbed_amount_mol_m2 - integrated_flux(b.dry_mass_kg_m2, q)?,
                    population.vacant_amount_mol_m2 - integrated_flux(b.dry_mass_kg_m2, vacant_q)?,
                ] {
                    if residual.abs() > controls.absolute_equation_tolerance_mol_m2[k] {
                        return Err(self.refuse(
                            p,
                            DiluteSoluteError::Sorption {
                                cell: i,
                                species: k,
                                residual_mol_m2: residual,
                                allowance_mol_m2: controls.absolute_equation_tolerance_mol_m2[k],
                            }
                            .into(),
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    fn solve_species(
        &mut self,
        k: usize,
        lower: &[f64],
        diagonal: &[f64],
        upper: &[f64],
        rhs: &mut [f64],
        solver: &mut TridiagonalWorkspace,
        controls: DiluteSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize> {
        let ctl = self.controls.langmuir;
        if !self.active[k] {
            if self.linear_solves == ctl.maximum_linear_solves {
                return Err(DiluteSoluteError::NonlinearWork {
                    required_linear_solves: self.linear_solves + 1,
                    maximum_linear_solves: ctl.maximum_linear_solves,
                });
            }
            solver.factor(lower, diagonal, upper, &mut *cancelled)?;
            solver.solve_factored_in_place(rhs, &mut *cancelled)?;
            self.linear_solves += 1;
            return Ok(1);
        }
        let initial: Vec<_> = self.initial_concentration.iter().map(|c| c[k]).collect();
        let limits: Vec<_> = self.limits.iter().map(|row| row[k]).collect();
        let tolerance = controls.absolute_equation_tolerance_mol_m2[k]
            .min(controls.amount_balance_tolerance_mol_m2[k] / rhs.len() as f64)
            * 0.25;
        let count = transport::nonlinear::solve(
            k,
            lower,
            diagonal,
            upper,
            rhs,
            &initial,
            |i| limits[i],
            solver,
            transport::nonlinear::Controls {
                tolerance_mol_m2: tolerance,
                maximum_iterations: ctl.maximum_iterations,
                maximum_line_searches: ctl.maximum_line_searches,
                remaining_linear_solves: ctl.maximum_linear_solves - self.linear_solves,
            },
            |i, c| self.response(i, k, c),
            cancelled,
        )?;
        self.linear_solves += count;
        Ok(count)
    }

    fn finish_cell(
        &mut self,
        i: usize,
        k: usize,
        c: f64,
        _: DiluteSoluteControls<N>,
    ) -> Result<()> {
        self.response(i, k, c).map(|_| ())
    }

    fn amounts(&self, i: usize, k: usize) -> [f64; 3] {
        [self.original[i][k], self.bound[i][k], self.exchange[i][k]]
    }
}

#[cfg(test)]
mod tests;
