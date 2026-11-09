//! Prepared reacting enthalpy changes at prescribed pressure.
//!
//! Reuse the original Gibbs composition owner and its safeguarded HP inverse.
//! Original source integrals resolve H-H_reference without subtracting rounded
//! absolute formation-energy totals. Chemical redistribution and its reacting
//! Cp remain active at every trial. Zero enthalpy change at a changed pressure
//! supplies an equilibrium isenthalpic endpoint; it is not an isentrope.
//!
//! This is a thermodynamic endpoint, not a finite-rate, throttle-flow or heat-
//! transfer law. A consumer must account for heat, shaft work and kinetic energy
//! before specifying the total enthalpy change. Entropy change is returned but
//! no reservoir entropy or process admissibility is inferred.

use super::expansion::path::Reference;
use super::*;

/// Allocation-free result of one original-source pressure/enthalpy query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumEnthalpyPressureSampleV1 {
    /// Actual absolute temperature [K].
    pub temperature: Temperature,
    /// Caller-prescribed pressure [Pa].
    pub pressure: Pressure,
    /// Reacting caloric, compressibility and acoustic response at this endpoint.
    pub response: EquilibriumThermodynamicResponseV1,
    /// Requested extensive H-H_reference [J], not heat unless the process allows it.
    pub target_enthalpy_change_j: f64,
    /// Actual extensive H-H_reference from the original source integrals [J].
    pub enthalpy_change_j: f64,
    /// Actual minus requested enthalpy change [J].
    pub enthalpy_residual_j: f64,
    /// Actual S-S_reference, including chemical redistribution and mixing [J/K].
    pub entropy_change_j_k: f64,
    /// Actual product-minus-feed mass divided by conserved feed mass.
    pub relative_mass_residual: f64,
    /// All Gibbs iterations, including the temperature bracket trials.
    pub composition_iterations: usize,
    /// Interior temperature trials, excluding the two bracket endpoints.
    pub temperature_iterations: usize,
}

/// One immutable equilibrium reference with a bounded, reusable HP workspace.
///
/// No source receipt is allocated during temperature trials. The original
/// trace-species KKT fallback may allocate once if needed. Failed or cancelled
/// queries invalidate publication; their numerical warm starts do not change
/// the original source, feed inventory or enthalpy datum.
pub struct EquilibriumEnthalpyPressurePathV1<'a> {
    owner: &'a IdealGasEquilibriumV1,
    work: Workspace,
    reference: Reference,
    reference_state: EquilibriumThermodynamicStateV1,
    controls: EquilibriumEnthalpySolveV1,
    current: Option<EquilibriumEnthalpyPressureSampleV1>,
}

impl<'a> EquilibriumEnthalpyPressurePathV1<'a> {
    /// Establish an equilibrium reference at explicit T,p and feed elements.
    ///
    /// The feed specifies conserved inventory, not an unreacted caloric inlet.
    /// The energy reference is the accepted INITIAL EQUILIBRIUM state. All
    /// eligible sources must cover the complete declared temperature bracket;
    /// crossing their NASA-9 charts requires explicit caloric continuity.
    pub fn new(
        owner: &'a IdealGasEquilibriumV1,
        temperature: Temperature,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
        controls: EquilibriumEnthalpySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, EquilibriumExpansionErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
            .contains(&temperature.value())
        {
            return Err(E::Invalid("enthalpy reference is outside the temperature bracket").into());
        }
        let mut work = Workspace::new(owner, pressure, feed)?;
        work.source_temperature_range(
            owner,
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        )?;
        let iterations = work.solve(owner, temperature, controls.equilibrium, &mut cancelled)?;
        let reference = Reference::new(owner, &work, temperature, &mut cancelled)?;
        let reference_state = work.publish_response(
            owner,
            temperature,
            controls.equilibrium,
            iterations,
            0,
            None,
            &mut cancelled,
        )?;
        let mass = work.initial_mass * work.scale;
        if !mass.is_finite() || mass <= 0.0 {
            return Err(E::Arithmetic("original enthalpy-pressure feed mass").into());
        }
        let current = Some(EquilibriumEnthalpyPressureSampleV1 {
            temperature,
            pressure,
            response: reference_state.response(),
            target_enthalpy_change_j: 0.0,
            enthalpy_change_j: 0.0,
            enthalpy_residual_j: 0.0,
            entropy_change_j_k: 0.0,
            relative_mass_residual: reference_state.equilibrium().relative_mass_residual(),
            composition_iterations: 0,
            temperature_iterations: 0,
        });
        poll(&mut cancelled)?;
        Ok(Self {
            owner,
            work,
            reference,
            reference_state,
            controls,
            current,
        })
    }

    /// Immutable reference with original species sources and reacting response.
    pub const fn reference_state(&self) -> &EquilibriumThermodynamicStateV1 {
        &self.reference_state
    }

    /// Conserved original feed mass [kg], not a rounded product-mass estimate.
    pub fn feed_mass_kg(&self) -> f64 {
        self.work.initial_mass * self.work.scale
    }

    /// Re-equilibrate at prescribed p and total H-H_reference [J].
    ///
    /// At fixed pressure with only pressure-volume work, this change is the
    /// signed supplied heat. At a different pressure it can represent an
    /// isenthalpic throttle or a separately closed heat/work balance; this
    /// owner never equates enthalpy change with heat without that assumption.
    /// Small changes retain their own absolute tolerance even beside a large
    /// formation reference. Reacting-composition residual sensitivity to that
    /// reference remains, just as on the original UV and isentropic paths.
    /// Bracketing errors report enthalpy CHANGE endpoints for this path.
    pub fn sample(
        &mut self,
        pressure: Pressure,
        enthalpy_change: Energy,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumEnthalpyPressureSampleV1, EquilibriumExpansionErrorV1> {
        let previous = self.current.take();
        poll(&mut cancelled)?;
        if !pressure.value().is_finite()
            || pressure.value() <= 0.0
            || !enthalpy_change.value().is_finite()
        {
            return Err(E::Invalid("positive pressure and finite enthalpy change required").into());
        }
        if let Some(mut same) = previous.filter(|s| {
            s.pressure == pressure && s.target_enthalpy_change_j == enthalpy_change.value()
        }) {
            same.composition_iterations = 0;
            same.temperature_iterations = 0;
            poll(&mut cancelled)?;
            self.current = Some(same);
            return Ok(same);
        }
        self.work.pressure = pressure;
        self.work.log_pressure_ratio = crate::entropy::log_ratio(
            pressure.value(),
            self.owner.convention.reference_pressure().value(),
        );
        let preferred = previous.map(|s| {
            // dH=Cp*dT + V*(1-T*alpha)*dp at fixed elemental inventory.
            // This predictor uses reacting expansion, including molecular-count
            // changes. Only the actual original-source residual accepts a root.
            let volume =
                self.feed_mass_kg() * (1.0 + s.relative_mass_residual) / s.response.density_kg_m3;
            let pressure_work = volume
                * (1.0 - s.temperature.value() * s.response.thermal_expansion_per_k)
                * (pressure.value() - s.pressure.value());
            s.temperature.value()
                + (enthalpy_change.value() - s.enthalpy_change_j - pressure_work)
                    / s.response.heat_capacity_p_j_k
        });
        let owner = self.owner;
        let reference = &self.reference;
        let (temperature, composition_iterations, temperature_iterations) =
            self.work.solve_enthalpy(
                owner,
                enthalpy_change.value(),
                self.controls,
                preferred,
                &mut cancelled,
                |work, t, cancel| reference.enthalpy_change(owner, work, t, cancel),
            )?;
        let response = self.work.response(
            owner,
            temperature,
            self.controls.equilibrium,
            &mut cancelled,
        )?;
        let dh = reference.enthalpy_change(owner, &self.work, temperature, &mut cancelled)?;
        let ds = reference.entropy_change(owner, &self.work, temperature, &mut cancelled)?;
        let mut mass = Sum::default();
        for &i in &self.work.active {
            mass.add(self.work.n[i] * owner.models[i].molar_mass().value());
        }
        let mass_residual = (mass.value() - self.work.initial_mass) / self.work.initial_mass;
        if !mass_residual.is_finite()
            || mass_residual.abs() > self.controls.equilibrium.relative_balance_tolerance
        {
            return Err(E::MassBalance {
                relative_residual: mass_residual,
            }
            .into());
        }
        let residual = dh - enthalpy_change.value();
        if !residual.is_finite() || residual.abs() > self.controls.enthalpy_tolerance_j {
            return Err(E::EnthalpyNoConvergence {
                iterations: temperature_iterations,
                residual_j: residual,
            }
            .into());
        }
        let sample = EquilibriumEnthalpyPressureSampleV1 {
            temperature,
            pressure,
            response,
            target_enthalpy_change_j: enthalpy_change.value(),
            enthalpy_change_j: dh,
            enthalpy_residual_j: residual,
            entropy_change_j_k: ds,
            relative_mass_residual: mass_residual,
            composition_iterations,
            temperature_iterations,
        };
        poll(&mut cancelled)?;
        self.current = Some(sample);
        Ok(sample)
    }

    /// Publish only the last complete query, with original per-species receipts.
    ///
    /// The target and residual remain in the sample: adding a small change to
    /// the absolute formation datum must not silently redefine this HP query.
    pub fn publish_current(
        &mut self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumThermodynamicStateV1, EquilibriumExpansionErrorV1> {
        poll(&mut cancelled)?;
        let sample = self
            .current
            .ok_or(E::Invalid("no accepted enthalpy-pressure query to publish"))?;
        let state = self.work.publish_response(
            self.owner,
            sample.temperature,
            self.controls.equilibrium,
            sample.composition_iterations,
            sample.temperature_iterations,
            None,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(state)
    }
}
