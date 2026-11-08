//! Adiabatic absorption/exsolution in a finite gas and finite liquid solvent.
//!
//! Both rigid, well-mixed phases share one temperature: thermal equilibration
//! is assumed fast compared with uptake. The caller supplies the solvent mass
//! and original single-phase liquid calorics. Negligible thermal expansion and
//! dissolved partial molar volume identify liquid caloric changes with internal
//! energy changes; chart density does not infer a solvent mass or a liquid EOS.
//!
//! Each stage reuses the exact finite-inventory isothermal film exchange, then
//! solves the CLOSED combined energy balance for its new composition. Original
//! NASA-9 gas Cv and original solvent/solute Cp integrals set the temperature.
//! Henry, diffusivity and optional gas-film sources are continued at that actual
//! temperature before further exchange. No thermostat supplies or removes heat.
//!
//! This is a first-order split integrator, not an exact nonisothermal kinetic
//! solution. Full/two-half-step differences in BOTH loading and temperature
//! drive bounded refinement. Their accumulated absolute estimates are reported;
//! they do not certify global error or source uncertainty. Every thermal solve
//! closes against the beginning of the requested interval, so many substeps do
//! not accumulate independent energy residuals or caloric-reference drift.
//!
//! All dilute, quasi-steady film assumptions of the isothermal owner remain.
//! No finite thermal-contact resistance, solvent evaporation, reaction, phase
//! change, thermal expansion or entropy-production theorem is supplied. The
//! author must align gas/solute caloric datums and supply thermodynamically
//! consistent Henry and caloric laws; this owner does not infer either law.

use super::*;

/// Explicit solvent thermal storage; dissolved-solute calorics stay in the law.
#[derive(Clone, Debug, PartialEq)]
pub struct DissolutionSolventCalorics {
    /// Conserved solvent mass [kg], excluding the exchanging solute.
    pub mass_kg: f64,
    /// Original liquid Cp(T), under the declared negligible-expansion convention.
    pub heat_capacity: ResolvedHeatCapacityCurve,
}

/// Absolute estimated errors over one requested interval and bounded work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdiabaticDissolutionControls {
    /// Existing source, film, dilute, energy and mass admission. Its maximum
    /// step bounds each trial; its gas temperature bracket bounds thermal solves.
    pub dissolution: DissolvedGasControls,
    pub temperature_tolerance_k: f64,
    pub dissolved_moles_tolerance: f64,
    /// Refinement refuses below this duration [s]; a shorter final remainder
    /// may be accepted if it already satisfies the local error estimate.
    pub minimum_step_s: f64,
    /// Each trial evaluates one full and two half stages, including rejected trials.
    pub maximum_trials: usize,
}

/// Accepted closed-system evolution. Positive compartment energy changes warm
/// or add inventory to that compartment; their sum is zero within the budget.
#[derive(Clone, Debug, PartialEq)]
pub struct AdiabaticDissolutionStep {
    pub initial: DissolvedGasInventory,
    pub state: DissolvedGasInventory,
    pub initial_temperature_k: f64,
    pub temperature_k: f64,
    pub sources: [DissolutionSample; 2],
    pub dissolved_moles_change: f64,
    pub gas_internal_energy_change_j: f64,
    pub dissolved_internal_energy_change_j: f64,
    pub solvent_internal_energy_change_j: f64,
    pub mass_residual_kg: f64,
    /// Source-relative closed energy residual; no external heat or work enters.
    pub energy_residual_j: f64,
    pub accepted_steps: usize,
    pub trials: usize,
    pub estimated_temperature_error_k: f64,
    pub estimated_dissolved_moles_error: f64,
}

#[derive(Clone)]
struct State {
    law: HenryDissolution,
    gas: ClosedGasParcel,
    liquid: DissolvedGasInventory,
    sample: DissolutionSample,
}

struct EnergyAnchor<'a> {
    law: &'a HenryDissolution,
    solvent: &'a DissolutionSolventCalorics,
    initial_moles: f64,
    gas_specific_internal_energy_j_kg: f64,
}

struct EnergyChange {
    gas: f64,
    dissolved: f64,
    solvent: f64,
    residual: f64,
    capacity: f64,
}

fn gas_source(error: impl core::fmt::Display) -> DissolvedGasError {
    GasEnergyError::Material(crate::MaterialError::Parameters {
        what: format!("adiabatic Henry gas calorics: {error}"),
    })
    .into()
}

impl EnergyAnchor<'_> {
    // Integrate relative calorics before adding the transferred inventory.
    // Formation references cancel in (h_d - u_g), not after multiplying two
    // large whole-compartment absolute stored energies.
    fn energy(&self, state: &State, temperature: f64) -> Result<EnergyChange> {
        let initial_t = self.law.temperature_k;
        let dm = finite(
            (state.liquid.dissolved_moles - self.initial_moles) * self.law.molar_mass_kg_mol,
        )?;
        let dissolved_mass = finite(state.liquid.dissolved_moles * self.law.molar_mass_kg_mol)?;
        let thermo = state.gas.model().thermodynamics();
        let gas_sensible = finite(
            state.gas.mass_kg()
                * thermo
                    .specific_internal_energy_change(
                        Temperature::new(initial_t),
                        Temperature::new(temperature),
                    )
                    .map_err(gas_source)?
                    .value(),
        )?;
        let dissolved_sensible = finite(
            dissolved_mass
                * self
                    .law
                    .calorics
                    .partial_enthalpy
                    .specific_enthalpy_change(initial_t, temperature)?,
        )?;
        let solvent = finite(
            self.solvent.mass_kg
                * self
                    .solvent
                    .heat_capacity
                    .specific_enthalpy_change(initial_t, temperature)?,
        )?;
        let capacity = finite(
            state.gas.mass_kg()
                * thermo
                    .evaluate_heat_capacities(Temperature::new(temperature))
                    .map_err(gas_source)?
                    .mass_specific_cv()
                    .value()
                + dissolved_mass
                    * self
                        .law
                        .calorics
                        .partial_enthalpy
                        .secant_heat_capacity(temperature, temperature)?
                + self.solvent.mass_kg
                    * self
                        .solvent
                        .heat_capacity
                        .secant_heat_capacity(temperature, temperature)?,
        )?;
        if !positive(capacity) {
            return Err(invalid(
                "adiabatic gas/solution needs positive thermal capacity",
            ));
        }
        Ok(EnergyChange {
            gas: finite(gas_sensible - dm * self.gas_specific_internal_energy_j_kg)?,
            dissolved: finite(dissolved_sensible + dm * self.law.dissolved_enthalpy_j_kg)?,
            solvent,
            residual: finite(
                dm * (self.law.dissolved_enthalpy_j_kg - self.gas_specific_internal_energy_j_kg)
                    + gas_sensible
                    + dissolved_sensible
                    + solvent,
            )?,
            capacity,
        })
    }

    fn temperature(
        &self,
        state: &State,
        controls: DissolvedGasControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        let solve = controls.gas.gas;
        let mut lo = solve.minimum_temperature_k;
        let mut hi = solve.maximum_temperature_k;
        let mut t = state.gas.state().state().temperature;
        let tolerance = controls.gas.energy_tolerance_j;
        let initial = self.energy(state, t)?;
        if initial.residual.abs() <= tolerance {
            return Ok(t);
        }
        // Positive original heat capacities make this root unique. Only the
        // endpoint in the required thermal direction needs evaluating.
        if initial.residual > 0.0 {
            hi = t;
            if self.energy(state, lo)?.residual > 0.0 {
                return Err(invalid(
                    "adiabatic dissolution cools below the thermal bracket",
                ));
            }
        } else {
            lo = t;
            if self.energy(state, hi)?.residual < 0.0 {
                return Err(invalid(
                    "adiabatic dissolution heats above the thermal bracket",
                ));
            }
        }
        for _ in 0..solve.maximum_iterations {
            poll(cancelled)?;
            let energy = self.energy(state, t)?;
            if energy.residual.abs() <= tolerance {
                return Ok(t);
            }
            if energy.residual > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
            let newton = t - energy.residual / energy.capacity;
            let next = if newton > lo && newton < hi {
                newton
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == t {
                return Err(invalid(
                    "adiabatic thermal solve stagnated above its energy budget",
                ));
            }
            t = next;
        }
        Err(invalid(
            "adiabatic thermal solve exhausted its iteration budget",
        ))
    }

    fn stage(
        &self,
        initial: &State,
        end_time_s: f64,
        film: DissolutionFilm,
        controls: DissolvedGasControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<State> {
        poll(cancelled)?;
        let mut next = initial.clone();
        next.law.advance(
            &mut next.gas,
            &mut next.liquid,
            film,
            end_time_s,
            controls,
            &mut *cancelled,
        )?;
        let temperature = self.temperature(&next, controls, cancelled)?;
        (next.gas, _) = next.gas.at_temperature(temperature)?;
        (next.law, next.sample) = next
            .law
            .at_conditions(&next.gas, next.liquid, film, controls)?;
        poll(cancelled)?;
        Ok(next)
    }
}

impl HenryDissolution {
    /// Advance a closed gas/solution at one common, evolving temperature.
    ///
    /// The solvent, Henry/diffusion and dissolved caloric sources stay original.
    /// Sources, gas and liquid publish together only after all requested time,
    /// energy, mass and error-estimate checks pass, including final cancellation.
    /// The supplied gas bracket must support every original species and both
    /// liquid caloric charts on a continuous positive-capacity branch.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adiabatic(
        &mut self,
        gas: &mut ClosedGasParcel,
        liquid: &mut DissolvedGasInventory,
        solvent: &DissolutionSolventCalorics,
        film: DissolutionFilm,
        end_time_s: f64,
        controls: AdiabaticDissolutionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdiabaticDissolutionStep> {
        poll(&mut cancelled)?;
        let exchange = controls.dissolution;
        let duration = end_time_s - liquid.time_s;
        let initial_temperature = gas.state().state().temperature;
        if !positive(duration)
            || !end_time_s.is_finite()
            || !positive(solvent.mass_kg)
            || !positive(controls.temperature_tolerance_k)
            || !positive(controls.dissolved_moles_tolerance)
            || !positive(controls.minimum_step_s)
            || !positive(exchange.maximum_step_s)
            || controls.minimum_step_s > exchange.maximum_step_s
            || controls.maximum_trials == 0
            || !positive(exchange.gas.energy_tolerance_j)
            || initial_temperature < exchange.gas.gas.minimum_temperature_k
            || initial_temperature > exchange.gas.gas.maximum_temperature_k
            || solvent
                .heat_capacity
                .curve()
                .knots()
                .iter()
                .any(|knot| knot.liquid_mass_fraction != 1.0)
        {
            return Err(invalid(
                "adiabatic dissolution needs finite liquid storage and bounded controls",
            ));
        }
        // Admit continuous positive Cv once for the complete original catalog,
        // including a species which may initially be absent from the gas.
        self.catalog
            .maximum_isochoric_temperature_slope(exchange.gas.gas)?;
        for chart in [&solvent.heat_capacity, &self.calorics.partial_enthalpy] {
            chart.specific_enthalpy_change(
                exchange.gas.gas.minimum_temperature_k,
                exchange.gas.gas.maximum_temperature_k,
            )?;
        }
        let (initial_law, initial_sample) = self.at_conditions(gas, *liquid, film, exchange)?;
        let anchor = EnergyAnchor {
            initial_moles: liquid.dissolved_moles,
            gas_specific_internal_energy_j_kg: finite(
                self.catalog
                    .species_specific_enthalpy_j_kg(&self.species, initial_temperature)?
                    - R * initial_temperature / self.molar_mass_kg_mol,
            )?,
            law: &initial_law,
            solvent,
        };
        let mut state = State {
            law: initial_law.clone(),
            gas: gas.clone(),
            liquid: *liquid,
            sample: initial_sample.clone(),
        };
        let mut step = duration.min(exchange.maximum_step_s);
        let mut trials = 0;
        let mut accepted = 0;
        let mut temperature_error = 0.0;
        let mut moles_error = 0.0;
        while state.liquid.time_s < end_time_s {
            poll(&mut cancelled)?;
            if trials == controls.maximum_trials {
                return Err(invalid(
                    "adiabatic dissolution exhausted its refinement budget",
                ));
            }
            trials += 1;
            let end = (state.liquid.time_s + step).min(end_time_s);
            let dt = end - state.liquid.time_s;
            let mid = state.liquid.time_s + 0.5 * dt;
            if !positive(dt) || mid == state.liquid.time_s || mid == end {
                return Err(invalid(
                    "adiabatic refinement has no representable half step",
                ));
            }
            let full = anchor.stage(&state, end, film, exchange, &mut cancelled)?;
            let half = anchor.stage(&state, mid, film, exchange, &mut cancelled)?;
            let fine = anchor.stage(&half, end, film, exchange, &mut cancelled)?;
            let dt_error = (full.law.temperature_k - fine.law.temperature_k).abs();
            let dn_error = (full.liquid.dissolved_moles - fine.liquid.dissolved_moles).abs();
            let fraction = dt / duration;
            let estimate = (dt_error / (controls.temperature_tolerance_k * fraction))
                .max(dn_error / (controls.dissolved_moles_tolerance * fraction));
            if !estimate.is_finite() {
                return Err(invalid(
                    "adiabatic local error allowance is unrepresentable",
                ));
            }
            if estimate <= 1.0 {
                state = fine;
                accepted += 1;
                temperature_error = finite(temperature_error + dt_error)?;
                moles_error = finite(moles_error + dn_error)?;
                if estimate < 0.125 {
                    step = (2.0 * dt).min(exchange.maximum_step_s);
                }
            } else {
                step = 0.5 * dt;
                if step < controls.minimum_step_s {
                    return Err(invalid(
                        "adiabatic error estimate requires a smaller permitted step",
                    ));
                }
            }
        }
        let energy = anchor.energy(&state, state.law.temperature_k)?;
        let dn = finite(state.liquid.dissolved_moles - liquid.dissolved_moles)?;
        let mass_residual = finite(
            state.gas.species_mass_kg(&self.species)? - gas.species_mass_kg(&self.species)?
                + dn * self.molar_mass_kg_mol,
        )?;
        budget(
            "adiabatic selected-species mass",
            mass_residual,
            exchange.gas.mass_tolerance_kg,
        )?;
        budget(
            "adiabatic closed system energy",
            energy.residual,
            exchange.gas.energy_tolerance_j,
        )?;
        let result = AdiabaticDissolutionStep {
            initial: *liquid,
            state: state.liquid,
            initial_temperature_k: initial_temperature,
            temperature_k: state.law.temperature_k,
            sources: [initial_sample, state.sample],
            dissolved_moles_change: dn,
            gas_internal_energy_change_j: energy.gas,
            dissolved_internal_energy_change_j: energy.dissolved,
            solvent_internal_energy_change_j: energy.solvent,
            mass_residual_kg: mass_residual,
            energy_residual_j: energy.residual,
            accepted_steps: accepted,
            trials,
            estimated_temperature_error_k: temperature_error,
            estimated_dissolved_moles_error: moles_error,
        };
        poll(&mut cancelled)?;
        *gas = state.gas;
        *liquid = state.liquid;
        *self = state.law;
        Ok(result)
    }
}
