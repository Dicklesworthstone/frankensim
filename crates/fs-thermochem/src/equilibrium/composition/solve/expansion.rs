//! Equilibrium entropy inversion and reacting thermodynamic response.
//!
//! The original element-constrained Gibbs solver owns composition. Its
//! differentiated KKT system, not a frozen cp/cv approximation, supplies
//! expansion, compressibility and equilibrium sound speed. Temperature trials
//! reuse its bounded workspace; only accepted endpoints allocate receipts.
//!
//! For fixed element inventory, let a=dln(V)/dln(T)|p and
//! b=-dln(V)/dln(p)|T. Then Cv=Cp-N*R*a*a/b and
//! c_eq^2=(p/rho)/(b-N*R*a*a/Cp). In general Cp-Cv != N*R.
//! These are the zero-frequency chemical-equilibrium derivatives, not the
//! frozen or finite-rate acoustic response. The neutral ideal-gas and source
//! restrictions of the original composition owner remain in force.

use super::*;

/// Temperature bracket and total-entropy tolerance for an equilibrium SP solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumEntropySolveV1 {
    /// Original composition tolerances and per-temperature work bounds.
    pub equilibrium: EquilibriumSolveV1,
    /// Positive source-supported lower temperature [K].
    pub minimum_temperature_k: f64,
    /// Source-supported upper temperature, greater than the lower bound [K].
    pub maximum_temperature_k: f64,
    /// Absolute total entropy residual allowed for the supplied feed [J/K].
    pub entropy_tolerance_j_k: f64,
    /// At most 256 temperature iterations, excluding bracket endpoints.
    pub maximum_temperature_iterations: usize,
}

/// Refusal of an equilibrium derivative or entropy-constrained endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumExpansionErrorV1 {
    /// Original chemistry, source, arithmetic, allocation or cancellation error.
    Equilibrium(CompositionEquilibriumErrorV1),
    /// Invalid controls or unsupported continuous temperature path.
    Invalid(&'static str),
    /// Both admitted endpoints lie on the same side of the requested entropy.
    EntropyNotBracketed {
        lower_entropy_j_k: f64,
        upper_entropy_j_k: f64,
    },
    /// Temperature convergence cannot substitute for entropy closure.
    EntropyNoConvergence {
        iterations: usize,
        residual_j_k: f64,
    },
    /// A differentiated conservation or stationarity equation did not close.
    DerivativeResidual(&'static str),
    /// A required finite positive response could not be established.
    UnstableResponse(&'static str),
}
impl From<CompositionEquilibriumErrorV1> for EquilibriumExpansionErrorV1 {
    fn from(value: CompositionEquilibriumErrorV1) -> Self {
        Self::Equilibrium(value)
    }
}
impl core::fmt::Display for EquilibriumExpansionErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium gas expansion: {self:?}")
    }
}
impl std::error::Error for EquilibriumExpansionErrorV1 {}

/// Local equilibrium derivatives at fixed elemental inventory.
/// Numerical estimates, not source uncertainty or interval certificates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumThermodynamicResponseV1 {
    /// Ideal-gas density using the actual equilibrium molecular amount [kg/m3].
    pub density_kg_m3: f64,
    /// dH/dT at fixed p and elements, including composition changes [J/K].
    pub heat_capacity_p_j_k: f64,
    /// dU/dT at fixed V and elements, including composition changes [J/K].
    pub heat_capacity_v_j_k: f64,
    /// (1/V) dV/dT at fixed p and elements [1/K]; may be signed.
    pub thermal_expansion_per_k: f64,
    /// -(1/V) dV/dp at fixed T and elements [1/Pa].
    pub isothermal_compressibility_per_pa: f64,
    /// dln(T)/dln(p) at fixed entropy and elements; may be signed.
    pub isentropic_temperature_exponent: f64,
    /// dln(p)/dln(rho) at fixed entropy and elements. Not Cp/Cv in general.
    pub isentropic_pressure_exponent: f64,
    /// sqrt(dp/drho at fixed entropy and elements) [m/s].
    pub sound_speed_m_s: f64,
}

/// An accepted equilibrium state and its consistent reacting response.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumThermodynamicStateV1 {
    equilibrium: IdealGasEquilibriumSolutionV1,
    entropy_j_k: f64,
    response: EquilibriumThermodynamicResponseV1,
    target_entropy_j_k: Option<f64>,
}
impl EquilibriumThermodynamicStateV1 {
    /// Original composition solution, including every source receipt.
    pub const fn equilibrium(&self) -> &IdealGasEquilibriumSolutionV1 {
        &self.equilibrium
    }
    /// Direct species-plus-mixing total entropy [J/K], not (H-G)/T.
    pub const fn entropy_j_k(&self) -> f64 {
        self.entropy_j_k
    }
    /// Coupled equilibrium derivatives at this same endpoint.
    pub const fn response(&self) -> EquilibriumThermodynamicResponseV1 {
        self.response
    }
    /// Actual total entropy minus the requested SP target [J/K].
    pub fn entropy_residual_j_k(&self) -> Option<f64> {
        self.target_entropy_j_k.map(|s| self.entropy_j_k - s)
    }
}

impl EquilibriumEntropySolveV1 {
    fn validate(self) -> Result<(), EquilibriumExpansionErrorV1> {
        self.equilibrium.validate()?;
        if !self.minimum_temperature_k.is_finite()
            || self.minimum_temperature_k <= 0.0
            || !self.maximum_temperature_k.is_finite()
            || self.maximum_temperature_k <= self.minimum_temperature_k
            || !self.entropy_tolerance_j_k.is_finite()
            || self.entropy_tolerance_j_k <= 0.0
            || !(1..=256).contains(&self.maximum_temperature_iterations)
        {
            return Err(EquilibriumExpansionErrorV1::Invalid(
                "invalid SP temperature bracket, entropy tolerance or work limit",
            ));
        }
        Ok(())
    }
}

impl IdealGasEquilibriumV1 {
    /// TP equilibrium with analytic reacting compressibility and sound speed.
    ///
    /// The two sensitivity right-hand sides use the same bounded element-space
    /// and full-KKT linear owners as composition and equilibrium Cp. No chemical
    /// finite differences, frozen-gamma substitution or new reaction model is used.
    pub fn equilibrate_tp_with_response(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
        controls: EquilibriumSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumThermodynamicStateV1, EquilibriumExpansionErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        let mut work = Workspace::new(self, pressure, feed)?;
        let iterations = work.solve(self, temperature, controls, &mut cancelled)?;
        work.publish_response(self, temperature, controls, iterations, 0, None, &mut cancelled)
    }

    /// Equilibrate at fixed total entropy [J/K], pressure and feed elements.
    ///
    /// Every trial re-equilibrates composition. The safeguarded Newton slope is
    /// dS/dT=Cp_equilibrium/T, not frozen Cp/T. One workspace is reused through
    /// all temperature trials, including the trace-species KKT fallback.
    /// Only a complete converged endpoint receives source receipts.
    ///
    /// Each eligible species must support the complete bracket. Crossing its
    /// NASA-9 regions requires its explicit continuous-caloric-reference policy;
    /// no source discontinuity or gap is silently treated as an isentrope.
    pub fn equilibrate_sp(
        &self,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
        target_entropy_j_k: f64,
        controls: EquilibriumEntropySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumThermodynamicStateV1, EquilibriumExpansionErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !target_entropy_j_k.is_finite() {
            return Err(EquilibriumExpansionErrorV1::Invalid("target entropy must be finite"));
        }
        let mut work = Workspace::new(self, pressure, feed)?;
        work.admit_temperature_path(self, controls)?;
        let (temperature, iterations, thermal_iterations) =
            work.solve_sp(self, target_entropy_j_k, controls, &mut cancelled)?;
        work.publish_response(
            self, temperature, controls.equilibrium, iterations, thermal_iterations,
            Some(target_entropy_j_k), &mut cancelled,
        )
    }
}

impl Workspace {
    fn admit_temperature_path(
        &self,
        owner: &IdealGasEquilibriumV1,
        controls: EquilibriumEntropySolveV1,
    ) -> Result<(), EquilibriumExpansionErrorV1> {
        for &i in &self.active {
            let model = &owner.models[i];
            let lower = model.region_for_temperature(controls.minimum_temperature_k)
                .ok_or(EquilibriumExpansionErrorV1::Invalid("SP lower temperature lacks source support"))?.0;
            let upper = model.region_for_temperature(controls.maximum_temperature_k)
                .ok_or(EquilibriumExpansionErrorV1::Invalid("SP upper temperature lacks source support"))?.0;
            if lower != upper && model.caloric_reference.is_none() {
                return Err(EquilibriumExpansionErrorV1::Invalid(
                    "cross-region equilibrium expansion requires explicit caloric continuity",
                ));
            }
        }
        Ok(())
    }

    fn entropy(&self, owner: &IdealGasEquilibriumV1, temperature: Temperature,
        cancelled: &mut impl FnMut() -> bool) -> Result<f64, E> {
        let mut total = Sum::default();
        for &i in &self.active { total.add(self.n[i]); }
        let log_total = ln(total.value());
        let mut entropy = Sum::default();
        for &i in &self.active {
            poll(cancelled)?;
            let properties = owner.models[i].evaluate_properties(temperature)
                .map_err(|source| owner.species_error(i, source))?;
            let log_activity = self.log_n[i] - log_total + self.log_pressure_ratio;
            entropy.add(self.n[i] * (properties.s().value()
                - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * log_activity));
        }
        let entropy = entropy.value() * self.scale;
        if entropy.is_finite() { Ok(entropy) } else { Err(E::Arithmetic("equilibrium total entropy")) }
    }

    // Reuse the original element matrix and dense solver for pressure response.
    // q_i=dln(n_i)/dln(p): q_i-q_total-A_i.lambda_p=-1, A diag(n) q=0.
    fn pressure_sensitivity(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), E> {
        poll(cancelled)?;
        let ns = self.n.len();
        let (total, _) = self.fill_matrix()?;
        let dim = self.rank + 1;
        self.matrix[self.rank * dim + self.rank] = 0.0;
        for k in 0..self.rank { self.rhs[k] = self.matrix[k * dim + self.rank]; }
        self.rhs[self.rank] = total;
        if linear_solve(&mut self.matrix, &mut self.rhs, dim, 64.0 * f64::EPSILON).is_ok() {
            for &i in &self.active {
                let mut q = Sum::default();
                q.add(-1.0);
                q.add(self.rhs[self.rank]);
                for k in 0..self.rank { q.add(self.a[k * ns + i] * self.rhs[k]); }
                self.step[i] = q.value();
            }
        } else {
            // Same full-KKT owner and buffers as solve_full, with the pressure
            // forcing (-1) in the species rows and zero conservation forcing.
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
                    self.full_matrix[(na + k) * dim + self.rank + 1 + z] = self.a[k * ns + i] * self.n[i];
                }
                self.full_matrix[z * dim + self.rank] = -1.0;
                self.full_matrix[z * dim + self.rank + 1 + z] = 1.0;
                self.full_rhs[z] = -1.0;
                self.full_matrix[(dim - 1) * dim + self.rank + 1 + z] = self.n[i];
            }
            self.full_matrix[(dim - 1) * dim + self.rank] = -total;
            linear_solve(&mut self.full_matrix, &mut self.full_rhs, dim, 0.0)?;
            self.rhs.copy_from_slice(&self.full_rhs[..self.rank + 1]);
            for (z, &i) in self.active.iter().enumerate() { self.step[i] = self.full_rhs[self.rank + 1 + z]; }
        }
        poll(cancelled)
    }

    fn admit_sensitivity(&self, owner: &IdealGasEquilibriumV1, thermal_rt: Option<f64>,
        controls: EquilibriumSolveV1) -> Result<f64, EquilibriumExpansionErrorV1> {
        let ns = self.n.len();
        let mut total = Sum::default();
        let mut weighted = Sum::default();
        let mut mass_change = Sum::default();
        let mut mass_scale = Sum::default();
        let mut largest = 1.0_f64;
        for &i in &self.active {
            let q = self.step[i];
            if !q.is_finite() { return Err(EquilibriumExpansionErrorV1::DerivativeResidual("nonfinite species response")); }
            total.add(self.n[i]);
            weighted.add(self.n[i] * q);
            let mass = self.n[i] * owner.models[i].molar_mass().value();
            mass_change.add(mass * q);
            mass_scale.add(mass * (1.0 + q.abs()));
            largest = largest.max(q.abs());
            let force = thermal_rt.map_or(-1.0, |rt| self.h[i] / rt);
            let mut residual = Sum::default();
            residual.add(q - self.rhs[self.rank]);
            residual.add(-force);
            let mut scale = 1.0 + q.abs() + self.rhs[self.rank].abs() + force.abs();
            for k in 0..self.rank {
                let v = self.a[k * ns + i] * self.rhs[k];
                residual.add(-v);
                scale += v.abs();
            }
            if !scale.is_finite() || !residual.value().is_finite()
                || residual.value().abs() > controls.dimensionless_potential_tolerance * scale {
                return Err(EquilibriumExpansionErrorV1::DerivativeResidual("chemical-potential response"));
            }
        }
        let tolerance = controls.relative_balance_tolerance + 512.0 * f64::EPSILON;
        for (k, &inventory) in self.inventory.iter().enumerate() {
            if inventory == 0.0 { continue; }
            let mut residual = Sum::default();
            let mut scale = Sum::default();
            scale.add(inventory);
            for &i in &self.active {
                let v = owner.coefficients[k * ns + i] * self.n[i] * self.step[i];
                residual.add(v);
                scale.add(v.abs());
            }
            if !scale.value().is_finite() || !residual.value().is_finite()
                || residual.value().abs() > tolerance * scale.value() {
                return Err(EquilibriumExpansionErrorV1::DerivativeResidual("element response"));
            }
        }
        let mean = weighted.value() / total.value();
        if !mean.is_finite() || (mean - self.rhs[self.rank]).abs() > tolerance * largest
            || !mass_change.value().is_finite() || !mass_scale.value().is_finite()
            || mass_change.value().abs() > tolerance * mass_scale.value() {
            return Err(EquilibriumExpansionErrorV1::DerivativeResidual("total amount or mass response"));
        }
        Ok(mean)
    }

    fn response(&mut self, owner: &IdealGasEquilibriumV1, temperature: Temperature,
        controls: EquilibriumSolveV1, cancelled: &mut impl FnMut() -> bool)
        -> Result<EquilibriumThermodynamicResponseV1, EquilibriumExpansionErrorV1> {
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let t = temperature.value();
        let p = self.pressure.value();
        let cp = self.heat_capacity(temperature, cancelled)?;
        let a = 1.0 + self.admit_sensitivity(owner, Some(r * t), controls)?;
        self.pressure_sensitivity(cancelled)?;
        let b = 1.0 - self.admit_sensitivity(owner, None, controls)?;
        let mut total = Sum::default();
        let mut mass = Sum::default();
        for &i in &self.active {
            total.add(self.n[i]);
            mass.add(self.n[i] * owner.models[i].molar_mass().value());
        }
        // Use normalized amounts for intensive properties, preserving feed-scale invariance.
        let nr = total.value() * r;
        let cp_normalized = cp / self.scale;
        if !nr.is_finite() || nr <= 0.0 || !cp_normalized.is_finite() || cp_normalized <= 0.0 {
            return Err(EquilibriumExpansionErrorV1::UnstableResponse("unrepresentable intensive heat capacity"));
        }
        let exponent = nr / cp_normalized * a;
        let adiabatic = b - a * exponent;
        if !b.is_finite() || b <= 0.0 || !adiabatic.is_finite() || adiabatic <= 0.0 {
            return Err(EquilibriumExpansionErrorV1::UnstableResponse("nonpositive equilibrium compressibility or Cv"));
        }
        let density = mass.value() * (p / (nr * t));
        let cv = cp * (adiabatic / b);
        let sound = fs_math::det::sqrt((p / density) / adiabatic);
        let response = EquilibriumThermodynamicResponseV1 {
            density_kg_m3: density,
            heat_capacity_p_j_k: cp,
            heat_capacity_v_j_k: cv,
            thermal_expansion_per_k: a / t,
            isothermal_compressibility_per_pa: b / p,
            isentropic_temperature_exponent: exponent,
            isentropic_pressure_exponent: 1.0 / adiabatic,
            sound_speed_m_s: sound,
        };
        for value in [density, cp, cv, b / p, 1.0 / adiabatic, sound] {
            if !value.is_finite() || value <= 0.0 {
                return Err(EquilibriumExpansionErrorV1::UnstableResponse("unrepresentable positive thermodynamic response"));
            }
        }
        if !(a / t).is_finite() || !exponent.is_finite() {
            return Err(EquilibriumExpansionErrorV1::UnstableResponse("unrepresentable thermal response"));
        }
        poll(cancelled)?;
        Ok(response)
    }

    fn solve_sp(&mut self, owner: &IdealGasEquilibriumV1, target: f64,
        controls: EquilibriumEntropySolveV1, cancelled: &mut impl FnMut() -> bool)
        -> Result<(Temperature, usize, usize), EquilibriumExpansionErrorV1> {
        let (mut low, mut high) = (controls.minimum_temperature_k, controls.maximum_temperature_k);
        let mut iterations = self.solve(owner, Temperature::new(low), controls.equilibrium, cancelled)?;
        let lower = self.entropy(owner, Temperature::new(low), cancelled)?;
        if (lower - target).abs() <= controls.entropy_tolerance_j_k {
            return Ok((Temperature::new(low), iterations, 0));
        }
        iterations += self.solve(owner, Temperature::new(high), controls.equilibrium, cancelled)?;
        let upper = self.entropy(owner, Temperature::new(high), cancelled)?;
        if (upper - target).abs() <= controls.entropy_tolerance_j_k {
            return Ok((Temperature::new(high), iterations, 0));
        }
        if !(lower < target && target < upper) {
            return Err(EquilibriumExpansionErrorV1::EntropyNotBracketed {
                lower_entropy_j_k: lower, upper_entropy_j_k: upper,
            });
        }
        let mut t = low + 0.5 * (high - low);
        let mut residual = upper - target;
        let mut steps = 0;
        for step in 1..=controls.maximum_temperature_iterations {
            poll(cancelled)?;
            steps = step;
            iterations += self.solve(owner, Temperature::new(t), controls.equilibrium, cancelled)?;
            residual = self.entropy(owner, Temperature::new(t), cancelled)? - target;
            if residual.abs() <= controls.entropy_tolerance_j_k {
                return Ok((Temperature::new(t), iterations, step));
            }
            if residual < 0.0 { low = t; } else { high = t; }
            let newton = match self.heat_capacity(Temperature::new(t), cancelled) {
                Ok(cp) => t - (residual / cp) * t,
                Err(E::Cancelled) => return Err(E::Cancelled.into()),
                Err(_) => f64::NAN,
            };
            let margin = 0.05 * (high - low);
            t = if newton.is_finite() && newton > low + margin && newton < high - margin {
                newton
            } else { low + 0.5 * (high - low) };
            if t == low || t == high { break; }
        }
        Err(EquilibriumExpansionErrorV1::EntropyNoConvergence { iterations: steps, residual_j_k: residual })
    }

    fn publish_response(&mut self, owner: &IdealGasEquilibriumV1, temperature: Temperature,
        controls: EquilibriumSolveV1, iterations: usize, thermal_iterations: usize,
        target_entropy_j_k: Option<f64>, cancelled: &mut impl FnMut() -> bool)
        -> Result<EquilibriumThermodynamicStateV1, EquilibriumExpansionErrorV1> {
        let entropy_j_k = self.entropy(owner, temperature, cancelled)?;
        let response = self.response(owner, temperature, controls, cancelled)?;
        let equilibrium = self.publish(owner, temperature, controls, iterations, thermal_iterations, None, cancelled)?;
        poll(cancelled)?;
        Ok(EquilibriumThermodynamicStateV1 { equilibrium, entropy_j_k, response, target_entropy_j_k })
    }
}

#[cfg(test)]
mod tests;
