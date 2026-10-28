//! Adiabatic steady-flow compressor/turbine closure on original gas enthalpy.
//!
//! Uniform low-speed ports with negligible kinetic/potential energy change,
//! one unchanged ideal-gas composition, no heat leak, no shaft inertia. The
//! caller explicitly selects overall isentropic or differential polytropic
//! efficiency; neither is inferred or calibrated.
//! Overall isentropic law: compressor Delta h = Delta h_s / eta,
//! turbine Delta h = eta Delta h_s.
//! Shaft work delivered BY gas is -Delta h; no p Delta V term is added again.
//! References: NASA Glenn compressor and power-turbine thermodynamics,
//! https://www.grc.nasa.gov/www/k-12/airplane/compth.html and
//! https://www.grc.nasa.gov/WWW/K-12/airplane/powtrbth.html.
//! Polytropic variable-Cp relation: Drela, TASOPT (2010), Appendix D.3,
//! https://web.mit.edu/drela/Public/N%2B3/Final_Report_App.pdf.
//! Constant local efficiency integrates Cp(T)/T, not a fitted constant exponent.

use super::{GasCaloricBranch, GasEnergyError, ResolvedThermochemicalGasState, coordinates,
    failure, poll, positive, updated_point};
use super::super::entropy::GasEntropyError;
use fs_matdb::QueryPoint;
use fs_qty::{Pressure, Temperature};

/// Explicit physical direction; a compressor cannot masquerade as a turbine.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GasMachineMode { Compressor, Turbine }

/// The declared efficiency law, not interchangeable numerical conventions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GasMachineEfficiency {
    /// Overall ratio of actual and isentropic enthalpy changes for a finite stage.
    Isentropic(f64),
    /// Constant differential efficiency: v dp / dh for a compressor, or
    /// dh / (v dp) for a turbine, integrated through the original Cp(T).
    Polytropic(f64),
}
impl GasMachineEfficiency {
    pub fn value(self) -> f64 {
        match self { Self::Isentropic(value) | Self::Polytropic(value) => value }
    }
}

/// No accepted state or shaft load is returned on any failed stage.
#[derive(Clone, Debug, PartialEq)]
pub enum GasMachineError {
    Invalid(&'static str),
    Energy(GasEnergyError),
    Entropy(GasEntropyError),
}
impl core::fmt::Display for GasMachineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f, "gas machine: {self:?}") }
}
impl std::error::Error for GasMachineError {}
impl From<GasEnergyError> for GasMachineError {
    fn from(error: GasEnergyError) -> Self { Self::Energy(error) }
}
impl From<GasEntropyError> for GasMachineError {
    fn from(error: GasEntropyError) -> Self { Self::Entropy(error) }
}

/// Analytic physical derivatives at fixed inlet state and composition.
/// Efficiency derivatives refer to the DECLARED efficiency law.
/// At eta=1 the efficiency derivative is the admissible one-sided limit.
/// Temperature derivatives use the selected Cp chart. At a Cp-discontinuous
/// source seam they are that chart's one-sided tangent, not a claim of two-sided
/// smoothness. Work/pressure derivatives do not divide by the ideal-state Cp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasMachineDerivatives {
    pub outlet_temperature_pressure_k_pa: f64,
    pub outlet_temperature_efficiency_k: f64,
    pub shaft_work_pressure_j_kg_pa: f64,
    pub shaft_work_efficiency_j_kg: f64,
}

/// Source-bearing steady-flow state pair with explicit work and entropy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct GasMachineStep {
    inlet: ResolvedThermochemicalGasState,
    outlet: ResolvedThermochemicalGasState,
    outlet_point: QueryPoint,
    pub mode: GasMachineMode,
    /// Authored law and its parameter; derivatives refer to this parameter.
    pub declared_efficiency: GasMachineEfficiency,
    /// Overall finite-stage isentropic efficiency. For a polytropic stage this
    /// is DERIVED from its enthalpy changes, not copied from the local efficiency.
    /// At equal pressures the declared efficiency is its zero-rise limit.
    pub isentropic_efficiency: f64,
    /// Caloric comparison only: no hypothetical transport state is fabricated.
    pub isentropic_temperature_k: f64,
    pub isentropic_enthalpy_change_j_kg: f64,
    pub isentropic_entropy_residual_j_kg_k: f64,
    pub enthalpy_change_j_kg: f64,
    /// Positive leaves gas through its shaft; negative is required shaft input.
    /// Isentropic law: from its prescribed enthalpy target. Polytropic law:
    /// inferred from the accepted original-law enthalpy rise/drop, not an
    /// independent imposed shaft-work datum or rotor-dynamics measurement.
    pub shaft_work_by_gas_j_kg: f64,
    /// Accepted Delta h plus shaft work [J/kg]. This is an independent endpoint
    /// residual for the isentropic law. For the polytropic law work is derived
    /// from Delta h, so this identity is not an independent convergence proof;
    /// `polytropic_entropy_residual_j_kg_k` provides its constitutive solve audit.
    pub energy_residual_j_kg: f64,
    /// Raw gas entropy increase [J/(kg K)]. It equals generated entropy under
    /// this adiabatic reduction, up to the declared absolute numerical tolerance.
    /// Small negative residuals within tolerance are reported, never clipped.
    pub entropy_change_j_kg_k: f64,
    pub entropy_tolerance_j_kg_k: f64,
    /// Actual Delta s minus the integrated local-efficiency entropy target.
    /// Present only for the polytropic law; never clipped to hide dissipation.
    pub polytropic_entropy_residual_j_kg_k: Option<f64>,
    pub derivatives: GasMachineDerivatives,
    pub iterations: usize,
}
impl GasMachineStep {
    pub fn inlet(&self) -> &ResolvedThermochemicalGasState { &self.inlet }
    pub fn outlet(&self) -> &ResolvedThermochemicalGasState { &self.outlet }
    pub fn outlet_point(&self) -> &QueryPoint { &self.outlet_point }
    pub fn into_outlet(self) -> (ResolvedThermochemicalGasState, QueryPoint) { (self.outlet, self.outlet_point) }

    /// Scale the specific work and entropy accounts at one nonnegative mass
    /// flow [kg/s]. This supplies a shared shaft/thermal-network load, not a
    /// second gas-state update or an extra heat source. Zero flow is exact zero.
    pub fn rates(&self, mass_flow_kg_s: f64) -> Result<GasMachineRates, GasMachineError> {
        if !mass_flow_kg_s.is_finite() || mass_flow_kg_s < 0.0 {
            return Err(GasMachineError::Invalid("finite nonnegative machine mass flow required"));
        }
        let product = |value: f64| {
            let result = mass_flow_kg_s * value;
            if !result.is_finite() || (mass_flow_kg_s != 0.0 && value != 0.0 && result == 0.0) {
                Err(GasMachineError::Invalid("machine rate is unrepresentable"))
            } else { Ok(result) }
        };
        Ok(GasMachineRates {
            shaft_power_by_gas_w: product(self.shaft_work_by_gas_j_kg)?,
            enthalpy_flow_change_w: product(self.enthalpy_change_j_kg)?,
            entropy_flow_increase_w_k: product(self.entropy_change_j_kg_k)?,
            energy_residual_w: product(self.energy_residual_j_kg)?,
        })
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasMachineRates {
    pub shaft_power_by_gas_w: f64,
    pub enthalpy_flow_change_w: f64,
    pub entropy_flow_increase_w_k: f64,
    pub energy_residual_w: f64,
}

impl GasCaloricBranch<'_> {
    /// Solve a complete adiabatic compressor/turbine stage at a prescribed
    /// outlet pressure and 0<eta<=1. Exact equal pressures give zero shaft work.
    ///
    /// The original relative entropy law finds the isentropic temperature;
    /// original relative enthalpy applies efficiency and finds the real outlet.
    /// Cp(T), mixture composition, reference policy and all material query axes
    /// survive. Only inlet and actual outlet require full transport. The ideal
    /// comparison remains a thermodynamic state, not an invented transport card.
    ///
    /// Source admission is reused, and trial roots allocate no receipts. The
    /// entropy root target is also bounded by the energy tolerance, preventing
    /// an inaccurate ideal-work reference from hiding behind an accurate final
    /// inversion. Independent signed-work, energy and second-law gates precede
    /// return. This is not a surge/choke map, rotor dynamics, real gas or a heat
    /// transfer model; no extra "efficiency loss heat" may be added to its gas.
    pub fn adiabatic_machine(
        &self, context: &QueryPoint, outlet_pressure: Pressure,
        mode: GasMachineMode, efficiency: f64, entropy_tolerance_j_kg_k: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMachineStep, GasMachineError> {
        self.solve_machine(context, outlet_pressure, mode,
            GasMachineEfficiency::Isentropic(efficiency), entropy_tolerance_j_kg_k, cancelled)
    }

    /// Adiabatic stage with constant DIFFERENTIAL polytropic efficiency.
    ///
    /// Compressor: integral Cp(T)/T dT = R ln(p_out/p_in)/eta.
    /// Turbine: integral Cp(T)/T dT = eta R ln(p_out/p_in).
    /// Original relative-entropy arithmetic cancels references before rounding;
    /// no constant gamma/exponent or pressure substepping is used. Chaining
    /// pressure subdivisions with the SAME local efficiency recovers the same
    /// endpoint/work within the declared numerical tolerances, unlike reusing
    /// one overall isentropic efficiency on arbitrarily subdivided stages.
    ///
    /// Shaft work is inferred from actual Delta h. The independently evaluated
    /// local-efficiency entropy residual, source support and work-direction
    /// checks govern acceptance. The report distinguishes authored polytropic
    /// efficiency from the resulting overall isentropic efficiency. Low-speed,
    /// frozen-gas, no-heat-leak assumptions match `adiabatic_machine`.
    pub fn polytropic_machine(
        &self, context: &QueryPoint, outlet_pressure: Pressure,
        mode: GasMachineMode, efficiency: f64, entropy_tolerance_j_kg_k: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMachineStep, GasMachineError> {
        self.solve_machine(context, outlet_pressure, mode,
            GasMachineEfficiency::Polytropic(efficiency), entropy_tolerance_j_kg_k, cancelled)
    }

    fn solve_machine(
        &self, context: &QueryPoint, outlet_pressure: Pressure,
        mode: GasMachineMode, efficiency_law: GasMachineEfficiency,
        entropy_tolerance_j_kg_k: f64, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMachineStep, GasMachineError> {
        let efficiency = efficiency_law.value();
        poll(&mut cancelled)?;
        let (t0, p0) = coordinates(context).map_err(GasEnergyError::from)?;
        let p1 = outlet_pressure.value();
        if !positive(p1) || !positive(efficiency) || efficiency > 1.0
            || !positive(entropy_tolerance_j_kg_k)
            || (mode == GasMachineMode::Compressor && p1 < p0.value())
            || (mode == GasMachineMode::Turbine && p1 > p0.value())
        {
            return Err(GasMachineError::Invalid("machine pressure direction, 0<efficiency<=1 and positive entropy tolerance required"));
        }
        let inlet = self.model.evaluate(context).map_err(GasEnergyError::from)?;
        let mut actual_branch = *self;
        actual_branch.controls.energy_tolerance_j_kg *= 0.5;
        if !positive(actual_branch.controls.energy_tolerance_j_kg) {
            return Err(GasMachineError::Invalid("machine energy accuracy target is unrepresentable"));
        }
        let entropy_target = entropy_tolerance_j_kg_k.min(
            (0.25 * self.controls.energy_tolerance_j_kg * efficiency) / self.controls.maximum_temperature_k);
        if !positive(entropy_target) {
            return Err(GasMachineError::Invalid("coupled entropy/enthalpy accuracy target is unrepresentable"));
        }
        let ideal = self.pressure_entropy_temperature(context, outlet_pressure, 0.0,
            entropy_target, &mut cancelled)?;
        let ideal_h = self.model.thermodynamics().specific_enthalpy_change(t0,
            Temperature::new(ideal.temperature_k))
            .map_err(|error| GasEnergyError::from(failure("machine ideal enthalpy", error)))?.value();
        let factor = match mode {
            GasMachineMode::Compressor => 1.0 / efficiency,
            GasMachineMode::Turbine => efficiency,
        };
        if !factor.is_finite() || !ideal_h.is_finite()
            || (p1 != p0.value() && (ideal_h == 0.0
                || ideal_h.is_sign_positive() != (mode == GasMachineMode::Compressor)))
        {
            return Err(GasMachineError::Invalid("ideal machine enthalpy change is unrepresentable"));
        }
        let pressure_entropy = if let GasMachineEfficiency::Polytropic(_) = efficiency_law {
            self.model.thermodynamics().specific_entropy_change(t0, p0, t0, outlet_pressure)
                .map_err(|error| GasEnergyError::from(failure("machine pressure entropy", error)))?.value()
        } else { 0.0 };
        let polytropic_target = if let GasMachineEfficiency::Polytropic(_) = efficiency_law {
            // Avoid factor-1 cancellation when efficiency is close to one.
            let loss = match mode {
                GasMachineMode::Compressor => -pressure_entropy * ((1.0 - efficiency) / efficiency),
                GasMachineMode::Turbine => pressure_entropy * (1.0 - efficiency),
            };
            if !loss.is_finite() || loss < 0.0
                || (p1 != p0.value() && efficiency != 1.0 && loss == 0.0)
            {
                return Err(GasMachineError::Invalid("polytropic entropy generation is unrepresentable"));
            }
            Some(loss)
        } else { None };
        let (outlet, actual_iterations, target_h) = if p1 == p0.value() {
            (inlet.clone(), 0, 0.0)
        } else if efficiency == 1.0 {
            poll(&mut cancelled)?;
            (self.model.evaluate(&updated_point(context, ideal.temperature_k, p1)?)
                .map_err(GasEnergyError::from)?, 0, ideal_h)
        } else if let Some(entropy_change) = polytropic_target {
            let actual = self.reconstruct_pressure_entropy_change(context, outlet_pressure,
                entropy_change, entropy_target, &mut cancelled)?;
            let actual_h = self.model.thermodynamics().specific_enthalpy_change(t0,
                Temperature::new(actual.state().state().temperature))
                .map_err(|error| GasEnergyError::from(failure("polytropic shaft enthalpy", error)))?.value();
            let iterations = actual.iterations();
            (actual.into_state(), iterations, actual_h)
        } else {
            let target_h = ideal_h * factor;
            if !target_h.is_finite() || target_h == 0.0 {
                return Err(GasMachineError::Invalid("actual machine enthalpy target is unrepresentable"));
            }
            let context_at_pressure = updated_point(context, t0.value(), p1)?;
            let actual = actual_branch.reconstruct_pressure_enthalpy_increment(&context_at_pressure,
                target_h, &mut cancelled)?;
            let iterations = actual.iterations();
            (actual.into_state(), iterations, target_h)
        };
        let t1 = outlet.state().temperature;
        let point = updated_point(context, t1, p1)?;
        let actual_h = self.model.thermodynamics().specific_enthalpy_change(t0, Temperature::new(t1))
            .map_err(|error| GasEnergyError::from(failure("machine accepted enthalpy", error)))?.value();
        let entropy = self.model.specific_entropy_change_j_kg_k(context, &point).map_err(GasEnergyError::from)?;
        let work = -target_h;
        let residual = actual_h + work;
        let polytropic_residual = polytropic_target.map(|target| entropy - target);
        if !residual.is_finite() || residual.abs() > self.controls.energy_tolerance_j_kg
            || !entropy.is_finite() || entropy < -entropy_tolerance_j_kg_k
            || polytropic_residual.is_some_and(|value| !value.is_finite() || value.abs() > entropy_target)
            || (p1 != p0.value() && (actual_h == 0.0
                || actual_h.is_sign_positive() != (mode == GasMachineMode::Compressor)))
        {
            return Err(GasMachineError::Invalid("machine endpoint violates energy, work direction or entropy tolerance"));
        }
        let overall_efficiency = if polytropic_target.is_none() || p1 == p0.value() || efficiency == 1.0 {
            efficiency
        } else {
            match mode {
                GasMachineMode::Compressor => ideal_h / actual_h,
                GasMachineMode::Turbine => actual_h / ideal_h,
            }
        };
        if !positive(overall_efficiency) || overall_efficiency > 1.0 {
            return Err(GasMachineError::Invalid("overall machine efficiency is outside represented physical bounds"));
        }
        // Isentropic finite-stage law uses Ts; local polytropic law uses the
        // actual T. Both follow dh = T ds + v dp without a frozen Cp/gamma.
        let work_temperature = if polytropic_target.is_some() { t1 } else { ideal.temperature_k };
        let dhdp = factor * inlet.state().specific_gas_constant * (work_temperature / p1);
        let dhdeta = if polytropic_target.is_some() {
            match mode {
                GasMachineMode::Compressor => (pressure_entropy / efficiency) * (t1 / efficiency),
                GasMachineMode::Turbine => -pressure_entropy * t1,
            }
        } else {
            match mode {
                GasMachineMode::Compressor => -target_h / efficiency,
                GasMachineMode::Turbine => ideal_h,
            }
        };
        let cp = outlet.thermodynamics().properties().mass_specific_cp().value();
        let derivatives = GasMachineDerivatives {
            outlet_temperature_pressure_k_pa: dhdp / cp,
            outlet_temperature_efficiency_k: dhdeta / cp,
            shaft_work_pressure_j_kg_pa: -dhdp,
            shaft_work_efficiency_j_kg: -dhdeta,
        };
        if !positive(cp) || !positive(dhdp) || !positive(derivatives.outlet_temperature_pressure_k_pa)
            || !dhdeta.is_finite() || !derivatives.outlet_temperature_efficiency_k.is_finite()
            || (dhdeta != 0.0 && derivatives.outlet_temperature_efficiency_k == 0.0)
        {
            return Err(GasMachineError::Invalid("machine tangent is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(GasMachineStep { inlet, outlet, outlet_point: point, mode,
            declared_efficiency: efficiency_law,
            isentropic_efficiency: overall_efficiency, isentropic_temperature_k: ideal.temperature_k,
            isentropic_enthalpy_change_j_kg: ideal_h,
            isentropic_entropy_residual_j_kg_k: ideal.sample.residual,
            enthalpy_change_j_kg: actual_h, shaft_work_by_gas_j_kg: work,
            energy_residual_j_kg: residual, entropy_change_j_kg_k: entropy,
            entropy_tolerance_j_kg_k, polytropic_entropy_residual_j_kg_k: polytropic_residual,
            derivatives, iterations: ideal.iterations + actual_iterations })
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod polytropic_tests;
