//! Finite immobile storage, coupled to the original mobile transport matrix.
use super::*;
use fs_material::adsorption::competitive::sites::scalar::{ScalarSiteStep, backward_euler_site};

/// One independent Langmuir population for one cell/species. Its fixed capacity
/// is in MOLES, supplied on an explicit solid skeleton basis by the caller.
/// `dn_bound/dt = k_a*c_mobile*(capacity-n_bound) - k_d*n_bound`.
/// No density, partition coefficient, cross-species competition or material
/// calibration is inferred. Zero rates admit irreversible/frozen limiting cases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmLangmuirSite {
    capacity_mol: f64,
    adsorption_m3_mol_s: f64,
    desorption_s_inverse: f64,
    initial_trapped_moles: f64,
}
impl FilmLangmuirSite {
    pub fn new(
        capacity_mol: f64,
        adsorption_m3_mol_s: f64,
        desorption_s_inverse: f64,
        initial_trapped_moles: f64,
    ) -> Result<Self> {
        if !positive(capacity_mol)
            || !nonnegative(adsorption_m3_mol_s)
            || !nonnegative(desorption_s_inverse)
            || !nonnegative(initial_trapped_moles)
            || initial_trapped_moles > capacity_mol
        {
            return Err(invalid(
                "invalid finite film site capacity, rates or loading",
            ));
        }
        Ok(Self {
            capacity_mol,
            adsorption_m3_mol_s,
            desorption_s_inverse,
            initial_trapped_moles,
        })
    }
    pub const fn capacity_mol(&self) -> f64 {
        self.capacity_mol
    }
    pub const fn adsorption_m3_mol_s(&self) -> f64 {
        self.adsorption_m3_mol_s
    }
    pub const fn desorption_s_inverse(&self) -> f64 {
        self.desorption_s_inverse
    }
}

impl TransientFilm {
    /// Add one optional finite kinetic population per cell/species, cell-major.
    /// Its initial trapped moles are ADDITIONAL to the initial mobile profile.
    /// If linear R is also supplied, the three independent pools are V*c,
    /// (R-1)*V*c and this finite loading. Neither bound pool contributes directly
    /// to osmotic pressure or face transport. All-None retains the original path.
    pub fn with_langmuir_sorption(mut self, sites: Vec<Option<FilmLangmuirSite>>) -> Result<Self> {
        if sites.len() != self.initial_concentrations_mol_m3.len() {
            return Err(invalid(
                "finite film sites require the cell-major species shape",
            ));
        }
        self.kinetic_sites = if sites.iter().all(Option::is_none) {
            None
        } else {
            Some(sites.into_boxed_slice())
        };
        Ok(self)
    }
}

/// Keep the mobile-plus-linear amount independently of finite bound inventory:
/// subtraction from a much larger total would erase a small mobile pool.
pub(super) struct Pools {
    pub equilibrium: Box<[f64]>,
    pub trapped: Box<[f64]>,
}
impl Pools {
    pub fn new(len: usize) -> Result<Self> {
        Ok(Self {
            equilibrium: zeros(len)?,
            trapped: zeros(len)?,
        })
    }
}
pub(super) struct KineticStorage {
    sites: Box<[Option<FilmLangmuirSite>]>,
    pub state: Pools,
    pub trial: Pools,
    pub capacities: Box<[f64]>,
    candidate: Box<[f64]>,
    diagonal: Box<[f64]>,
    residual: Box<[f64]>,
    correction: Box<[f64]>,
}
impl KineticStorage {
    pub fn bytes(entries: usize, rows: usize) -> usize {
        // Include caller definitions during preparation, both retained pools,
        // both staged pools, and five reusable nonlinear arrays.
        2 * entries * size_of::<Option<FilmLangmuirSite>>()
            + (4 * entries + 5 * rows) * size_of::<f64>()
            + size_of::<Self>()
    }
    pub fn new(
        definitions: &[Option<TransientFilm>; 2],
        amounts: &[f64],
        rows: usize,
    ) -> Result<Self> {
        let mut sites = Vec::new();
        sites
            .try_reserve_exact(amounts.len())
            .map_err(|_| invalid("film site allocation failed"))?;
        for film in definitions.iter().flatten() {
            match &film.kinetic_sites {
                Some(values) => sites.extend_from_slice(values),
                None => sites.resize(sites.len() + film.initial_concentrations_mol_m3.len(), None),
            }
        }
        let mut state = Pools::new(amounts.len())?;
        state.equilibrium.copy_from_slice(amounts);
        for (trapped, site) in state.trapped.iter_mut().zip(&sites) {
            *trapped = site.map_or(0.0, |s| s.initial_trapped_moles);
        }
        Ok(Self {
            sites: sites.into_boxed_slice(),
            state,
            trial: Pools::new(amounts.len())?,
            capacities: zeros(rows)?,
            candidate: zeros(rows)?,
            diagonal: zeros(rows)?,
            residual: zeros(rows)?,
            correction: zeros(rows)?,
        })
    }

    pub fn solve(
        &mut self,
        work: &mut Workspace,
        species: usize,
        count: usize,
        dt: f64,
        controls: TransientStepControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64)> {
        let equation = Equation {
            sites: &self.sites,
            old_trapped: &self.state.trapped,
            capacities: &self.capacities,
            lower: &work.lower,
            upper: &work.upper,
            base_diagonal: &work.diagonal,
            old: &work.rhs,
            species,
            count,
            dt,
            controls,
        };
        work.solution.copy_from_slice(&work.rhs);
        let mut norm = equation.evaluate(
            &work.solution,
            &mut self.diagonal,
            &mut self.residual,
            cancelled,
        )?;
        if norm.0 > 0.5 {
            // Zero is a subsolution: desorption can only add mobile inventory.
            // Langmuir endpoint loading is increasing and concave in c; the
            // positive diagonal derivative preserves the tridiagonal M-matrix.
            work.solution.fill(0.0);
            norm = equation.evaluate(
                &work.solution,
                &mut self.diagonal,
                &mut self.residual,
                cancelled,
            )?;
        }
        for _ in 0..controls.maximum_sorption_iterations {
            if norm.0 <= 0.5 {
                break;
            }
            for residual in &mut self.residual {
                *residual = -*residual;
            }
            work.linear.solve_into(
                &work.lower,
                &self.diagonal,
                &work.upper,
                &self.residual,
                &mut self.correction,
                &mut *cancelled,
            )?;
            let mut fraction: f64 = 1.0;
            let mut accepted = false;
            for _ in 0..controls.maximum_sorption_backtracks {
                poll(cancelled)?;
                for ((next, old), delta) in self
                    .candidate
                    .iter_mut()
                    .zip(&work.solution)
                    .zip(&self.correction)
                {
                    *next = fraction.mul_add(*delta, *old);
                }
                if self.candidate.iter().all(|&x| nonnegative(x)) {
                    let next = equation.evaluate(
                        &self.candidate,
                        &mut self.diagonal,
                        &mut self.residual,
                        cancelled,
                    )?;
                    if next.0 <= 0.5 || next.1 < norm.1 {
                        work.solution.copy_from_slice(&self.candidate);
                        norm = next;
                        accepted = true;
                        break;
                    }
                }
                fraction *= 0.5;
            }
            if !accepted {
                return Err(TransientOsmosisError::SorptionConvergence { species });
            }
        }
        if norm.0 > 0.5 {
            return Err(TransientOsmosisError::SorptionConvergence { species });
        }
        let mut balance = adaptive::Sum::default();
        let mut before = adaptive::Sum::default();
        let mut maximum = norm.1;
        for row in 0..work.solution.len() {
            if row % 128 == 0 {
                poll(cancelled)?;
            }
            let value = work.solution[row];
            balance.add(value - work.rhs[row]);
            before.add(work.rhs[row]);
            let mut trapped = 0.0;
            if row > 0 && row + 1 < work.solution.len() {
                let index = (row - 1) * count + species;
                if let Some(step) = equation.site(row, value)? {
                    // Stored amounts can balance even when a sub-ULP change
                    // suppresses physical uptake in a large preloaded pool.
                    // Audit the owner's INDEPENDENT kinetic equation at the
                    // converged candidate, using the mobile amount budget.
                    // Its separate occupied+vacant capacity closure retains
                    // the scalar owner's capacity-scaled roundoff gate.
                    let tolerance = controls.absolute_amount_tolerance_mol
                        + controls.relative_amount_tolerance * value.max(work.rhs[row]);
                    if step.kinetic_balance_residual_mol > tolerance {
                        return Err(TransientOsmosisError::SorptionConvergence { species });
                    }
                    maximum = maximum.max(step.kinetic_balance_residual_mol);
                    trapped = step.occupied_moles;
                }
                balance.add(trapped - self.state.trapped[index]);
                before.add(self.state.trapped[index]);
                self.trial.equilibrium[index] = value;
                self.trial.trapped[index] = trapped;
            }
            work.trial_amounts[row * count + species] = finite(value + trapped)?;
        }
        let balance = balance.total()?;
        if balance.abs()
            > controls.absolute_amount_tolerance_mol
                + controls.relative_amount_tolerance * before.total()?
        {
            return Err(TransientOsmosisError::InventoryBalance {
                species,
                residual_mol: balance,
            });
        }
        Ok((maximum, balance))
    }
}

struct Equation<'a> {
    sites: &'a [Option<FilmLangmuirSite>],
    old_trapped: &'a [f64],
    capacities: &'a [f64],
    lower: &'a [f64],
    upper: &'a [f64],
    base_diagonal: &'a [f64],
    old: &'a [f64],
    species: usize,
    count: usize,
    dt: f64,
    controls: TransientStepControls,
}
impl Equation<'_> {
    fn site(&self, row: usize, value: f64) -> Result<Option<ScalarSiteStep>> {
        if row == 0 || row + 1 == self.old.len() {
            return Ok(None);
        }
        let index = (row - 1) * self.count + self.species;
        self.sites[index]
            .map(|site| {
                let rate = finite(site.adsorption_m3_mol_s * (value / self.capacities[row]))?;
                if value > 0.0 && site.adsorption_m3_mol_s > 0.0 && rate == 0.0 {
                    return Err(invalid("film adsorption rate underflow"));
                }
                backward_euler_site(
                    self.old_trapped[index],
                    site.capacity_mol,
                    rate,
                    site.desorption_s_inverse,
                    self.dt,
                )
                .map_err(TransientOsmosisError::Sorption)
            })
            .transpose()
    }
    fn evaluate(
        &self,
        values: &[f64],
        diagonal: &mut [f64],
        residual: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64)> {
        let mut scaled: f64 = 0.0;
        let mut maximum: f64 = 0.0;
        for row in 0..values.len() {
            if row % 128 == 0 {
                poll(cancelled)?;
            }
            let mut difference = values[row] - self.old[row];
            // Use the independently retained mobile/linear inventory, never a
            // large bound capacity or gross throughflow, to scale algebraic
            // error. Otherwise those quantities could hide mobile-state error.
            let scale = values[row].max(self.old[row]);
            diagonal[row] = self.base_diagonal[row];
            if let Some(step) = self.site(row, values[row])? {
                let index = (row - 1) * self.count + self.species;
                let change = step.occupied_moles - self.old_trapped[index];
                difference += change;
                let site = self.sites[index].expect("evaluated finite site");
                diagonal[row] = finite(
                    diagonal[row]
                        + step.adsorption_rate_derivative_mol_s
                            * (site.adsorption_m3_mol_s / self.capacities[row]),
                )?;
            }
            if row > 0 {
                let incoming = (-self.lower[row - 1])
                    .mul_add(values[row - 1], self.upper[row - 1] * values[row]);
                difference -= incoming;
            }
            if row + 1 < values.len() {
                let outgoing =
                    (-self.lower[row]).mul_add(values[row], self.upper[row] * values[row + 1]);
                difference += outgoing;
            }
            residual[row] = finite(difference)?;
            maximum = maximum.max(difference.abs());
            scaled = scaled.max(
                difference.abs()
                    / (self.controls.absolute_amount_tolerance_mol
                        + self.controls.relative_amount_tolerance * scale),
            );
        }
        Ok((finite(scaled)?, maximum))
    }
}

impl TransientOsmoticPair {
    /// Finite kinetic loading only; the instantaneous linear pool is excluded.
    pub fn film_trapped_solute_moles(
        &self,
        side: usize,
        cell: usize,
        species: usize,
    ) -> Option<f64> {
        self.film_solute_moles(side, cell)?.get(species)?;
        Some(self.trapped_amount(
            1 + cell + if side == 0 { 0 } else { self.left_cells() },
            species,
        ))
    }
    /// Authored capacity and rates; loading evolves separately. None means the
    /// specified cell/species has no finite kinetic population (or is invalid).
    pub fn film_langmuir_site(
        &self,
        side: usize,
        cell: usize,
        species: usize,
    ) -> Option<&FilmLangmuirSite> {
        self.film_solute_moles(side, cell)?.get(species)?;
        let index =
            (cell + if side == 0 { 0 } else { self.left_cells() }) * self.species_count + species;
        self.kinetic.as_ref()?.sites[index].as_ref()
    }
}
