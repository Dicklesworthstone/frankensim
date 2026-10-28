//! Reusable log-amount Newton workspace. The ordinary solve uses an element
//! system of at most 17 unknowns. A bounded full KKT buffer is allocated lazily
//! once if traces make the reduced Gram system numerically singular.
//! Source receipts are allocated only at accepted publication.

use super::*;
use CompositionEquilibriumErrorV1 as E;
use fs_math::det::{exp, ln};

pub(super) mod expansion;

pub(super) struct Workspace {
    pressure: Pressure,
    log_pressure_ratio: f64,
    scale: f64,
    initial_mass: f64,
    inventory: Vec<f64>,
    active: Vec<usize>,
    a: Vec<f64>,
    rank: usize,
    log_n: Vec<f64>,
    n: Vec<f64>,
    log_total: f64,
    lambda: Vec<f64>,
    trial_log_n: Vec<f64>,
    trial_n: Vec<f64>,
    trial_lambda: Vec<f64>,
    g: Vec<f64>,
    h: Vec<f64>,
    cp: Vec<f64>,
    mu: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    step: Vec<f64>,
    full_matrix: Vec<f64>,
    full_rhs: Vec<f64>,
}
impl Workspace {
    pub(super) fn new(
        owner: &IdealGasEquilibriumV1,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
    ) -> Result<Self, E> {
        let p = pressure.value();
        if !p.is_finite() || p <= 0.0 {
            return Err(E::Invalid(
                "equilibrium pressure must be positive and finite",
            ));
        }
        let ns = owner.models.len();
        let ne = owner.elemental.row_count();
        let mut initial = owner.feed_amounts(feed)?;
        let mut sum = Sum::default();
        for &x in &initial {
            sum.add(x);
        }
        let scale = sum.value();
        let mut initial_mass = Sum::default();
        for (i, x) in initial.iter_mut().enumerate() {
            let value = *x / scale;
            if *x > 0.0 && value == 0.0 {
                return Err(E::Invalid("feed dynamic range cannot be normalized"));
            }
            *x = value;
            initial_mass.add(value * owner.models[i].molar_mass().value());
        }
        let initial_mass = initial_mass.value();
        if !initial_mass.is_finite() || initial_mass <= 0.0 {
            return Err(E::Arithmetic("initial mass"));
        }
        let mut inventory = allocated(ne, 0.0)?;
        for (k, b) in inventory.iter_mut().enumerate() {
            let mut value = Sum::default();
            for (i, &x) in initial.iter().enumerate() {
                value.add(owner.coefficients[k * ns + i] * x);
            }
            *b = value.value();
            if !b.is_finite() {
                return Err(E::Arithmetic("element inventory"));
            }
        }
        let mut active = Vec::new();
        active
            .try_reserve_exact(ns)
            .map_err(|_| E::AllocationRefused)?;
        for i in 0..ns {
            if (0..ne).all(|k| inventory[k] > 0.0 || owner.coefficients[k * ns + i] == 0.0) {
                active.push(i);
            }
        }
        if active.is_empty() {
            return Err(E::Invalid("no products represent the feed elements"));
        }
        let rows = independent_rows(owner, &active, &inventory)?;
        let rank = rows.len();
        if rank == 0 {
            return Err(E::SingularSystem);
        }
        // Every elemental equation uses its own inventory as scale.
        let mut a = allocated(rank * ns, 0.0)?;
        for (k, &row) in rows.iter().enumerate() {
            for &i in &active {
                let value = owner.coefficients[row * ns + i] / inventory[row];
                if !value.is_finite() {
                    return Err(E::Arithmetic("scaled element coefficients"));
                }
                a[k * ns + i] = value;
            }
        }
        let mut log_n = allocated(ns, f64::NEG_INFINITY)?;
        let mut n = allocated(ns, 0.0)?;
        let mut total = Sum::default();
        for &i in &active {
            let mut allowance = f64::INFINITY;
            for (k, &b) in inventory.iter().enumerate() {
                let count = owner.coefficients[k * ns + i];
                if count > 0.0 {
                    allowance = allowance.min(b / count);
                }
            }
            // Initialization only: all seeds together consume <=10% of each
            // element. No seed amount is retained as a physical concentration.
            let x = 0.9 * initial[i] + 0.1 * allowance / (active.len() as f64);
            if !x.is_finite() || x <= 0.0 {
                return Err(E::Arithmetic("positive interior initialization"));
            }
            n[i] = x;
            log_n[i] = ln(x);
            total.add(x);
        }
        let dim = rank + 1;
        Ok(Self {
            pressure,
            log_pressure_ratio: ln(p) - ln(owner.convention.reference_pressure().value()),
            scale,
            initial_mass,
            inventory,
            active,
            a,
            rank,
            log_n,
            n,
            log_total: ln(total.value()),
            lambda: allocated(rank, 0.0)?,
            trial_log_n: allocated(ns, f64::NEG_INFINITY)?,
            trial_n: allocated(ns, 0.0)?,
            trial_lambda: allocated(rank, 0.0)?,
            g: allocated(ns, 0.0)?,
            h: allocated(ns, 0.0)?,
            cp: allocated(ns, 0.0)?,
            mu: allocated(ns, 0.0)?,
            matrix: allocated(dim * dim, 0.0)?,
            rhs: allocated(dim, 0.0)?,
            step: allocated(ns, 0.0)?,
            full_matrix: Vec::new(),
            full_rhs: Vec::new(),
        })
    }

    pub(super) fn solve(
        &mut self,
        owner: &IdealGasEquilibriumV1,
        temperature: Temperature,
        controls: EquilibriumSolveV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize, E> {
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value();
        if !rt.is_finite() || rt <= 0.0 {
            return Err(E::Invalid("temperature must be positive and finite"));
        }
        for &i in &self.active {
            poll(cancelled)?;
            let p = owner.models[i]
                .evaluate_properties(temperature)
                .map_err(|source| owner.species_error(i, source))?;
            if p.cp().value() <= UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K {
                return Err(E::Invalid(
                    "eligible ideal-gas species require positive cv=cp-R",
                ));
            }
            self.g[i] = p.h().value() / rt - p.s().value() / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
            self.h[i] = p.h().value();
            self.cp[i] = p.cp().value();
            if !self.g[i].is_finite() {
                return Err(E::Arithmetic("dimensionless standard Gibbs energy"));
            }
        }
        for iteration in 0..=controls.maximum_iterations {
            poll(cancelled)?;
            let residuals = self.residuals(owner, false, self.log_total)?;
            if converged(residuals, controls) {
                return Ok(iteration);
            }
            if iteration == controls.maximum_iterations {
                return Err(E::NoConvergence {
                    iterations: iteration,
                    residuals,
                });
            }
            self.newton_system()?;
            let ns = self.n.len();
            if linear_solve(
                &mut self.matrix,
                &mut self.rhs,
                self.rank + 1,
                64.0 * f64::EPSILON,
            )
            .is_ok()
            {
                for &i in &self.active {
                    let mut value = Sum::default();
                    value.add(self.rhs[self.rank] - self.mu[i]);
                    for k in 0..self.rank {
                        value.add(self.a[k * ns + i] * self.rhs[k]);
                    }
                    self.step[i] = value.value();
                }
            } else {
                // Keep the species equations explicit when forming A diag(n)
                // A^T loses the distinctions carried by trace species.
                self.solve_full(None, cancelled)?;
            }
            let mut largest = self.rhs[self.rank].abs();
            for &i in &self.active {
                if !self.step[i].is_finite() {
                    return Err(E::Arithmetic("log-amount increment"));
                }
                // Damping growth prevents overshoot. Downward log changes may
                // be large; positivity and the merit search still guard them.
                largest = largest.max(self.step[i]);
            }
            let current_merit = merit(residuals, controls);
            let mut alpha = if largest > 2.0 { 2.0 / largest } else { 1.0 };
            let mut accepted = false;
            for _ in 0..controls.maximum_backtracks {
                poll(cancelled)?;
                let trial_total = self.log_total + alpha * self.rhs[self.rank];
                for k in 0..self.rank {
                    self.trial_lambda[k] = self.lambda[k] + alpha * self.rhs[k];
                }
                for &i in &self.active {
                    self.trial_log_n[i] = self.log_n[i] + alpha * self.step[i];
                    self.trial_n[i] = exp(self.trial_log_n[i]);
                }
                if let Ok(trial) = self.residuals(owner, true, trial_total) {
                    if converged(trial, controls)
                        || merit(trial, controls) <= (1.0 - 1e-4 * alpha) * current_merit
                    {
                        core::mem::swap(&mut self.log_n, &mut self.trial_log_n);
                        core::mem::swap(&mut self.n, &mut self.trial_n);
                        core::mem::swap(&mut self.lambda, &mut self.trial_lambda);
                        self.log_total = trial_total;
                        accepted = true;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if !accepted {
                return Err(E::NoConvergence {
                    iterations: iteration + 1,
                    residuals,
                });
            }
        }
        unreachable!("bounded Newton loop returns on its final iteration")
    }

    fn residuals(
        &self,
        owner: &IdealGasEquilibriumV1,
        trial: bool,
        log_total: f64,
    ) -> Result<EquilibriumResidualsV1, E> {
        let (amounts, logs, lambda) = if trial {
            (&self.trial_n, &self.trial_log_n, &self.trial_lambda)
        } else {
            (&self.n, &self.log_n, &self.lambda)
        };
        let ns = amounts.len();
        let mut sum = Sum::default();
        for &i in &self.active {
            if !amounts[i].is_finite() || !logs[i].is_finite() {
                return Err(E::Arithmetic("trial amount"));
            }
            sum.add(amounts[i]);
        }
        let sum = sum.value();
        let total = exp(log_total);
        if !sum.is_finite() || sum <= 0.0 || !total.is_finite() || total <= 0.0 {
            return Err(E::Arithmetic("total moles"));
        }
        let mut maximum_element = 0.0_f64;
        for (k, &b) in self.inventory.iter().enumerate() {
            let mut error = Sum::default();
            error.add(-b);
            for &i in &self.active {
                error.add(owner.coefficients[k * ns + i] * amounts[i]);
            }
            let error = error.value();
            if !error.is_finite() {
                return Err(E::Arithmetic("element residual"));
            }
            if b == 0.0 {
                if error != 0.0 {
                    return Err(E::Arithmetic("absent element appeared"));
                }
            } else {
                maximum_element = maximum_element.max((error / b).abs());
            }
        }
        let mut maximum_potential = 0.0_f64;
        let actual_log_total = ln(sum);
        for &i in &self.active {
            let mut value = Sum::default();
            value.add(self.g[i]);
            value.add(logs[i] - log_total);
            value.add(self.log_pressure_ratio);
            for k in 0..self.rank {
                value.add(-self.a[k * ns + i] * lambda[k]);
            }
            let value = value.value();
            if !value.is_finite() {
                return Err(E::Arithmetic("chemical potential"));
            }
            maximum_potential = maximum_potential
                .max(value.abs())
                .max((value + log_total - actual_log_total).abs());
        }
        Ok(EquilibriumResidualsV1 {
            maximum_relative_element_residual: maximum_element,
            relative_total_moles_residual: ((sum - total) / total).abs(),
            maximum_dimensionless_potential_residual: maximum_potential,
        })
    }

    fn fill_matrix(&mut self) -> Result<(f64, f64), E> {
        let ns = self.n.len();
        let dim = self.rank + 1;
        let mut sum = Sum::default();
        for &i in &self.active {
            sum.add(self.n[i]);
        }
        let sum = sum.value();
        let total = exp(self.log_total);
        for k in 0..self.rank {
            let mut an = Sum::default();
            for &i in &self.active {
                an.add(self.a[k * ns + i] * self.n[i]);
            }
            let an = an.value();
            self.matrix[k * dim + self.rank] = an;
            self.matrix[self.rank * dim + k] = an;
            for l in 0..=k {
                let mut value = Sum::default();
                for &i in &self.active {
                    value.add((self.a[k * ns + i] * self.n[i]) * self.a[l * ns + i]);
                }
                self.matrix[k * dim + l] = value.value();
                self.matrix[l * dim + k] = self.matrix[k * dim + l];
            }
        }
        self.matrix[self.rank * dim + self.rank] = sum - total;
        if self.matrix.iter().any(|x| !x.is_finite()) {
            return Err(E::Arithmetic("element-space matrix"));
        }
        Ok((sum, total))
    }

    fn newton_system(&mut self) -> Result<(), E> {
        let ns = self.n.len();
        let (sum, total) = self.fill_matrix()?;
        for &i in &self.active {
            let mut mu = Sum::default();
            mu.add(self.g[i]);
            mu.add(self.log_n[i] - self.log_total);
            mu.add(self.log_pressure_ratio);
            for k in 0..self.rank {
                mu.add(-self.a[k * ns + i] * self.lambda[k]);
            }
            self.mu[i] = mu.value();
        }
        for k in 0..self.rank {
            let mut right = Sum::default();
            right.add(1.0);
            for &i in &self.active {
                right.add((self.a[k * ns + i] * self.n[i]) * (self.mu[i] - 1.0));
            }
            self.rhs[k] = right.value();
        }
        let mut right = Sum::default();
        right.add(total - sum);
        for &i in &self.active {
            right.add(self.n[i] * self.mu[i]);
        }
        self.rhs[self.rank] = right.value();
        Ok(())
    }

    // Unknown order [element potentials, log-total, species log-amounts]
    // avoids first eliminating species into a numerically rank-lost Gram matrix.
    fn solve_full(
        &mut self,
        thermal_rt: Option<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), E> {
        let ns = self.n.len();
        let na = self.active.len();
        let dim = self.rank + 1 + na;
        if self.full_matrix.len() != dim * dim {
            self.full_matrix = allocated(dim * dim, 0.0)?;
            self.full_rhs = allocated(dim, 0.0)?;
        } else {
            self.full_matrix.fill(0.0);
            self.full_rhs.fill(0.0);
        }
        for (z, &i) in self.active.iter().enumerate() {
            poll(cancelled)?;
            for k in 0..self.rank {
                self.full_matrix[z * dim + k] = -self.a[k * ns + i];
            }
            self.full_matrix[z * dim + self.rank] = -1.0;
            self.full_matrix[z * dim + self.rank + 1 + z] = 1.0;
            self.full_rhs[z] = match thermal_rt {
                Some(rt) => self.h[i] / rt,
                None => -self.mu[i],
            };
        }
        for k in 0..self.rank {
            let mut residual = Sum::default();
            residual.add(1.0);
            for (z, &i) in self.active.iter().enumerate() {
                let value = self.a[k * ns + i] * self.n[i];
                self.full_matrix[(na + k) * dim + self.rank + 1 + z] = value;
                residual.add(-value);
            }
            self.full_rhs[na + k] = if thermal_rt.is_some() {
                0.0
            } else {
                residual.value()
            };
        }
        let total = exp(self.log_total);
        let mut residual = Sum::default();
        residual.add(total);
        self.full_matrix[(dim - 1) * dim + self.rank] = -total;
        for (z, &i) in self.active.iter().enumerate() {
            self.full_matrix[(dim - 1) * dim + self.rank + 1 + z] = self.n[i];
            residual.add(-self.n[i]);
        }
        self.full_rhs[dim - 1] = if thermal_rt.is_some() {
            0.0
        } else {
            residual.value()
        };
        // The larger fallback retains very small pivots; actual nonlinear
        // residual admission, finite arithmetic and backtracking remain binding.
        linear_solve(&mut self.full_matrix, &mut self.full_rhs, dim, 0.0)?;
        self.rhs.copy_from_slice(&self.full_rhs[..self.rank + 1]);
        for (z, &i) in self.active.iter().enumerate() {
            self.step[i] = self.full_rhs[self.rank + 1 + z];
        }
        poll(cancelled)
    }

    fn enthalpy(&self) -> Result<f64, E> {
        let mut value = Sum::default();
        for &i in &self.active {
            value.add(self.n[i] * self.h[i]);
        }
        let value = value.value() * self.scale;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(E::Arithmetic("total enthalpy"))
        }
    }

    // q_i=dln(n_i)/dln(T), holding p and every element fixed. The chemical
    // contribution is the positive Gibbs-Hessian form R sum n_i(q_i-mean(q))².
    // This avoids cancellation between large species formation enthalpies.
    fn heat_capacity(
        &mut self,
        temperature: Temperature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, E> {
        poll(cancelled)?;
        let ns = self.n.len();
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value();
        self.fill_matrix()?;
        for k in 0..self.rank {
            let mut right = Sum::default();
            for &i in &self.active {
                right.add(-(self.a[k * ns + i] * self.n[i]) * (self.h[i] / rt));
            }
            self.rhs[k] = right.value();
        }
        let mut right = Sum::default();
        let mut sum = Sum::default();
        for &i in &self.active {
            right.add(-self.n[i] * (self.h[i] / rt));
            sum.add(self.n[i]);
        }
        self.rhs[self.rank] = right.value();
        if linear_solve(
            &mut self.matrix,
            &mut self.rhs,
            self.rank + 1,
            64.0 * f64::EPSILON,
        )
        .is_ok()
        {
            for &i in &self.active {
                let mut value = Sum::default();
                value.add(self.h[i] / rt);
                value.add(self.rhs[self.rank]);
                for k in 0..self.rank {
                    value.add(self.a[k * ns + i] * self.rhs[k]);
                }
                self.step[i] = value.value();
            }
        } else {
            self.solve_full(Some(rt), cancelled)?;
        }
        let mut mean = Sum::default();
        for &i in &self.active {
            mean.add(self.n[i] * self.step[i]);
        }
        let mean = mean.value() / sum.value();
        let mut capacity = Sum::default();
        for &i in &self.active {
            let q = self.step[i] - mean;
            capacity.add(self.n[i] * (self.cp[i] + UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * q * q));
        }
        let capacity = capacity.value() * self.scale;
        if capacity.is_finite() && capacity > 0.0 {
            Ok(capacity)
        } else {
            Err(E::Arithmetic("equilibrium heat capacity"))
        }
    }

    pub(super) fn solve_hp(
        &mut self,
        owner: &IdealGasEquilibriumV1,
        target: Energy,
        controls: EquilibriumEnthalpySolveV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IdealGasEquilibriumSolutionV1, E> {
        let (mut low, mut high) = (
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        );
        let mut iterations = self.solve(
            owner,
            Temperature::new(low),
            controls.equilibrium,
            cancelled,
        )?;
        let lower_h = self.enthalpy()?;
        if (lower_h - target.value()).abs() <= controls.enthalpy_tolerance_j {
            return self.publish(
                owner,
                Temperature::new(low),
                controls.equilibrium,
                iterations,
                0,
                Some(target),
                cancelled,
            );
        }
        iterations += self.solve(
            owner,
            Temperature::new(high),
            controls.equilibrium,
            cancelled,
        )?;
        let upper_h = self.enthalpy()?;
        if (upper_h - target.value()).abs() <= controls.enthalpy_tolerance_j {
            return self.publish(
                owner,
                Temperature::new(high),
                controls.equilibrium,
                iterations,
                0,
                Some(target),
                cancelled,
            );
        }
        if !(lower_h < target.value() && target.value() < upper_h) {
            return Err(E::EnthalpyNotBracketed {
                lower_enthalpy_j: lower_h,
                upper_enthalpy_j: upper_h,
            });
        }
        let mut t = low + (high - low) * 0.5;
        let mut last_residual = upper_h - target.value();
        let mut temperature_iterations = 0;
        for step in 1..=controls.maximum_temperature_iterations {
            poll(cancelled)?;
            temperature_iterations = step;
            iterations +=
                self.solve(owner, Temperature::new(t), controls.equilibrium, cancelled)?;
            let residual = self.enthalpy()? - target.value();
            if residual.abs() <= controls.enthalpy_tolerance_j {
                return self.publish(
                    owner,
                    Temperature::new(t),
                    controls.equilibrium,
                    iterations,
                    step,
                    Some(target),
                    cancelled,
                );
            }
            if residual < 0.0 {
                low = t;
            } else {
                high = t;
            }
            let midpoint = low + (high - low) * 0.5;
            let capacity = match self.heat_capacity(Temperature::new(t), cancelled) {
                Ok(cp) => Some(cp),
                Err(E::Cancelled) => return Err(E::Cancelled),
                Err(_) => None,
            };
            let newton = capacity.map(|cp| t - residual / cp);
            let margin = (high - low) * 0.05;
            t = match newton {
                Some(next) if next.is_finite() && next > low + margin && next < high - margin => {
                    next
                }
                _ => midpoint,
            };
            last_residual = residual;
            if t == low || t == high {
                break;
            }
        }
        Err(E::EnthalpyNoConvergence {
            iterations: temperature_iterations,
            residual_j: last_residual,
        })
    }

    pub(super) fn publish(
        &mut self,
        owner: &IdealGasEquilibriumV1,
        temperature: Temperature,
        controls: EquilibriumSolveV1,
        composition_iterations: usize,
        temperature_iterations: usize,
        target_enthalpy: Option<Energy>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IdealGasEquilibriumSolutionV1, E> {
        poll(cancelled)?;
        let residuals = self.residuals(owner, false, self.log_total)?;
        if !converged(residuals, controls) {
            return Err(E::NoConvergence {
                iterations: composition_iterations,
                residuals,
            });
        }
        let mut total = Sum::default();
        let mut mass = Sum::default();
        let mut gibbs = Sum::default();
        for &i in &self.active {
            total.add(self.n[i]);
            mass.add(self.n[i] * owner.models[i].molar_mass().value());
        }
        let total = total.value();
        let mass = mass.value();
        let relative_mass_residual = (mass - self.initial_mass) / self.initial_mass;
        if !relative_mass_residual.is_finite()
            || relative_mass_residual.abs() > controls.relative_balance_tolerance
        {
            return Err(E::MassBalance {
                relative_residual: relative_mass_residual,
            });
        }
        let physical_total = total * self.scale;
        let physical_mass = mass * self.scale;
        if !physical_total.is_finite()
            || physical_total <= 0.0
            || !physical_mass.is_finite()
            || physical_mass <= 0.0
        {
            return Err(E::Arithmetic("published total amount or mass"));
        }
        let log_total = ln(total);
        let mut species = Vec::new();
        species
            .try_reserve_exact(self.n.len())
            .map_err(|_| E::AllocationRefused)?;
        for i in 0..self.n.len() {
            poll(cancelled)?;
            let (log_fraction, fraction, amount, source) = if self.log_n[i].is_finite() {
                let log_fraction = self.log_n[i] - log_total;
                let fraction = exp(log_fraction);
                let amount = if self.n[i] > 0.0 {
                    self.n[i] * self.scale
                } else {
                    exp(self.log_n[i] + ln(self.scale))
                };
                if !fraction.is_finite() || !amount.is_finite() {
                    return Err(E::Arithmetic("published species amount"));
                }
                let evaluation = owner.models[i]
                    .evaluate(temperature)
                    .map_err(|source| owner.species_error(i, source))?;
                gibbs.add(self.n[i] * (self.g[i] + log_fraction + self.log_pressure_ratio));
                (
                    Some(log_fraction),
                    fraction,
                    amount,
                    Some(evaluation.receipt().clone()),
                )
            } else {
                (None, 0.0, 0.0, None)
            };
            species.push(EquilibriumSpeciesAmountV1 {
                species: owner.models[i].species().clone(),
                amount: Amount::new(amount),
                mole_fraction: fraction,
                log_mole_fraction: log_fraction,
                source,
            });
        }
        let gibbs = (gibbs.value() * self.scale)
            * (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * temperature.value());
        if !gibbs.is_finite() {
            return Err(E::Arithmetic("total Gibbs energy"));
        }
        let enthalpy = self.enthalpy()?;
        let equilibrium_heat_capacity_p_j_k = match self.heat_capacity(temperature, cancelled) {
            Ok(cp) => Some(cp),
            Err(E::Cancelled) => return Err(E::Cancelled),
            Err(_) => None,
        };
        poll(cancelled)?;
        Ok(IdealGasEquilibriumSolutionV1 {
            temperature,
            pressure: self.pressure,
            species,
            total_amount: Amount::new(physical_total),
            mass: Mass::new(physical_mass),
            enthalpy: Energy::new(enthalpy),
            gibbs_energy: Energy::new(gibbs),
            equilibrium_heat_capacity_p_j_k,
            residuals,
            relative_mass_residual,
            composition_iterations,
            temperature_iterations,
            target_enthalpy,
        })
    }
}
fn converged(r: EquilibriumResidualsV1, c: EquilibriumSolveV1) -> bool {
    r.maximum_relative_element_residual <= c.relative_balance_tolerance
        && r.relative_total_moles_residual <= c.relative_balance_tolerance
        && r.maximum_dimensionless_potential_residual <= c.dimensionless_potential_tolerance
}
fn merit(r: EquilibriumResidualsV1, c: EquilibriumSolveV1) -> f64 {
    (r.maximum_relative_element_residual / c.relative_balance_tolerance)
        .max(r.relative_total_moles_residual / c.relative_balance_tolerance)
        .max(r.maximum_dimensionless_potential_residual / c.dimensionless_potential_tolerance)
}

// Numerical basis admission; every original elemental row is still tested.
// This is not an exact symbolic conservation or rank certificate.
fn independent_rows(
    owner: &IdealGasEquilibriumV1,
    active: &[usize],
    inventory: &[f64],
) -> Result<Vec<usize>, E> {
    let ns = owner.models.len();
    let ne = inventory.len();
    let mut basis = allocated(ne * ns, 0.0)?;
    let mut row = allocated(ns, 0.0)?;
    let mut pivots: Vec<usize> = Vec::new();
    let mut rows = Vec::new();
    pivots
        .try_reserve_exact(ne)
        .map_err(|_| E::AllocationRefused)?;
    rows.try_reserve_exact(ne)
        .map_err(|_| E::AllocationRefused)?;
    for (k, &b) in inventory.iter().enumerate() {
        if b == 0.0 {
            continue;
        }
        row.fill(0.0);
        let scale = active
            .iter()
            .map(|&i| owner.coefficients[k * ns + i])
            .fold(0.0_f64, f64::max);
        for &i in active {
            row[i] = owner.coefficients[k * ns + i] / scale;
        }
        for (r, &pivot) in pivots.iter().enumerate() {
            let factor = row[pivot];
            for &i in active {
                row[i] -= factor * basis[r * ns + i];
            }
        }
        let mut pivot = active[0];
        for &i in active {
            if row[i].abs() > row[pivot].abs() {
                pivot = i;
            }
        }
        if row[pivot].abs() <= 128.0 * f64::EPSILON {
            continue;
        }
        let r = rows.len();
        for &i in active {
            basis[r * ns + i] = row[i] / row[pivot];
        }
        rows.push(k);
        pivots.push(pivot);
    }
    Ok(rows)
}

// Row-equilibrated, partial-pivoted dense solve in deterministic tie order.
// Fast element-space path is <=17; fallback is <=145. No inner allocation.
fn linear_solve(a: &mut [f64], rhs: &mut [f64], n: usize, pivot_floor: f64) -> Result<(), E> {
    for i in 0..n {
        let scale = a[i * n..(i + 1) * n]
            .iter()
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
        if !scale.is_finite() || scale == 0.0 || !rhs[i].is_finite() {
            return Err(E::SingularSystem);
        }
        for j in 0..n {
            a[i * n + j] /= scale;
        }
        rhs[i] /= scale;
    }
    for k in 0..n {
        let mut pivot = k;
        for i in k + 1..n {
            if a[i * n + k].abs() > a[pivot * n + k].abs() {
                pivot = i;
            }
        }
        if !a[pivot * n + k].is_finite() || a[pivot * n + k].abs() <= pivot_floor {
            return Err(E::SingularSystem);
        }
        if pivot != k {
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            rhs.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.0;
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for i in (0..n).rev() {
        let mut value = Sum::default();
        value.add(rhs[i]);
        for j in i + 1..n {
            value.add(-a[i * n + j] * rhs[j]);
        }
        rhs[i] = value.value() / a[i * n + i];
        if !rhs[i].is_finite() {
            return Err(E::Arithmetic("linear system solution"));
        }
    }
    Ok(())
}
