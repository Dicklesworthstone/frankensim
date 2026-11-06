//! Entropy admission of the existing axial molecular SSP-RK2 evolution.
//!
//! This measures the actual ideal-mixture state increment, including conversion
//! of kinetic energy through the original total-energy inversion. It does not
//! identify that finite-step increment with a constitutive production integral.
//! Spatial and temporal discretization remain part of the measured result.
//!
//! For an arbitrary fixed positive T* and species mass densities rho_i*, write
//! phi_i(T) = s_i^0(T)-s_i^0(T*)-R ln(T/T*). The entropy density, apart from a
//! constant linear in each species inventory, is
//!   sum_i [rho_i phi_i(T)-R(rho_i ln(rho_i/rho_i*)-rho_i)] / M_i.
//! The omitted constants multiply only numerical closed-inventory residuals;
//! they are reported separately, never hidden in the second-law admission.
//! Density differences and the relative-entropy remainder are evaluated before
//! rounding, preserving weak mixing and cancellation of source reference data.

use super::super::wall::log_ratio;
use super::super::{MAX_SPECIES, MixtureEuler1d, MixtureGas, Primitive, allocated};
use super::{
    Boundaries, Boundary, DissipationControl, DissipationError, DissipationOperator,
    DissipationReport, TransportLaw, poll, positive,
};
use fs_math::dd::Dd;
use fs_qty::{Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

mod reservoir;
use reservoir::ReservoirEntropyAccount;

/// Absolute entropy and relative species-continuity allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularEntropyControl {
    /// Maximum admitted negative state-minus-reservoir-input increment, J/K,
    /// at every accepted substep AND over the complete original-to-final interval.
    /// The reservoir input is zero for the closed entry point.
    pub absolute_tolerance_j_k: f64,
    /// Closed entry: maximum absolute mass drift divided by original species mass.
    /// With reservoirs, subtract actual species transfer and scale by original
    /// mass plus accepted absolute stage/end throughput. Exact-zero inventory
    /// and throughput require zero drift. The continuity defect is independent
    /// of its arbitrary entropy-reference effect.
    pub relative_species_tolerance: f64,
}

impl Default for MolecularEntropyControl {
    fn default() -> Self {
        Self {
            absolute_tolerance_j_k: 1e-10,
            relative_species_tolerance: 1e-12,
        }
    }
}

impl MolecularEntropyControl {
    fn validate(self) -> Result<(), DissipationError> {
        if !self.absolute_tolerance_j_k.is_finite()
            || self.absolute_tolerance_j_k < 0.0
            || !self.relative_species_tolerance.is_finite()
            || !(0.0..=1e-3).contains(&self.relative_species_tolerance)
        {
            return Err(DissipationError::Invalid(
                "closed molecular entropy controls",
            ));
        }
        Ok(())
    }
}

/// Actual closed-interval state entropy; no physical-production theorem is implied.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MolecularEntropyReport {
    /// Original-to-final thermal and mixing entropy change [J/K], with the
    /// explicitly reported species-reference continuity defect removed.
    /// This is the quantity admitted at every accepted substep and at completion.
    pub entropy_change_j_k: f64,
    /// Same original source-state entropy change INCLUDING its linear inventory
    /// defect [J/K]. Large arbitrary source datums can dominate this quantity.
    pub source_entropy_change_j_k: f64,
    /// Sum_i alpha_i * delta(m_i)/M_i [J/K], where
    /// alpha_i=s_i^0(T*)-R ln((rho_i*/M_i) R T*/p0)-R.
    /// T* is the original first-cell temperature and rho_i* the original
    /// volume-mean density (a positive original trace if that mean underflows).
    pub reference_entropy_residual_j_k: f64,
    /// Smallest admitted substep entropy change [J/K]; zero for zero duration.
    pub minimum_step_entropy_change_j_k: f64,
    /// Largest per-species relative inventory defect over accepted substeps
    /// and the complete interval. This is not an element-only mass check.
    pub maximum_relative_species_residual: f64,
    /// Attempts rejected specifically by entropy or species-inventory admission.
    /// These attempts are also counted by the native dissipation report.
    pub rejected_entropy_steps: usize,
}


/// Entropy of an open molecular interval with explicit stationary source baths.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReservoirMolecularEntropyReport {
    /// Original thermochemical source-state entropy change [J/K].
    pub source_entropy_change_j_k: f64,
    /// Actual entropy entering the gas from the left and right baths [J/K].
    /// Heat and partial-species entropy use their own actual bath states.
    pub boundary_entropy_input_j_k: [f64; 2],
    /// State change minus bath input, with species-reference contributions
    /// from numerical continuity residuals removed [J/K]. This is admitted.
    pub entropy_generation_j_k: f64,
    /// Same generation INCLUDING the explicitly bounded continuity defect [J/K].
    pub source_entropy_generation_j_k: f64,
    /// Source-reference entropy of (species change - boundary transfer) [J/K].
    pub reference_entropy_residual_j_k: f64,
    /// Smallest admitted substep state-minus-input value [J/K].
    pub minimum_step_entropy_generation_j_k: f64,
    /// Largest species residual divided by original mass plus accepted absolute
    /// stage/end throughput. Exact-zero inventory and throughput require zero drift.
    pub maximum_relative_species_residual: f64,
    pub rejected_entropy_steps: usize,
}

impl<L: TransportLaw> DissipationOperator<L> {


    /// Advance the closed molecular field with source-conjugate interior faces
    /// and admit every actual SSP-RK2 candidate against the original entropy.
    ///
    /// Interior diffusion uses logarithmic chemical efforts and divided
    /// enthalpies derived from each original caloric model. Fourier heat,
    /// longitudinal stress and authored reciprocal Soret/Dufour terms share
    /// one conservative energy flux. This corrects the spatial caloric/force
    /// pairing; the native finite-step positivity and entropy gates remain
    /// necessary. It does not establish a fully discrete stability theorem.
    ///
    /// Every active species must be positive at both sides of each interior
    /// face; absence on both sides is admitted. One-sided zeros refuse without
    /// flooring composition or selecting another flux. Periodic and stationary
    /// zero-heat/species ends are supported. Reservoir, moving-wall, RKL2 and
    /// advective/reactive integration are separate entry points.
    pub fn advance_entropy_consistent_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DissipationReport, MolecularEntropyReport), DissipationError> {
        controls.validate()?;
        entropy.validate()?;
        let mut account = EntropyTracker::new(flow, entropy, &mut cancelled)?;
        let report = super::evolve::advance_entropy_consistent(
            self, flow, duration, controls, &mut account, &mut cancelled,
        )?;
        Ok((report, account.report))
    }

    /// Advance source-conjugate molecular transport against stationary baths.
    ///
    /// Both interior and bath faces use the original divided caloric enthalpies,
    /// logarithmic chemical efforts, Fourier heat and authored Soret/Dufour
    /// pair. The cell-to-bath layer has the physical half-cell distance and
    /// arithmetic velocity in stress work. This gives the same paired
    /// semidiscrete entropy production as an interior face, with the bath
    /// entropy removed from the gas's accounting boundary.
    ///
    /// Bath entropy uses the ACTUAL conservative energy and species transfers:
    /// E_ref/T_b - sum_i(mu_i,ref/T_b * dm_i/M_i). Equivalently, bath heat is
    /// E_ref minus the bath partial-enthalpy transfer. It includes the divided-
    /// enthalpy correction and the bath share of viscous conversion; substituting
    /// the interior Fourier/Dufour heat here would omit those contributions.
    /// The stationary bath has no external traction-work input.
    ///
    /// Every actual SSP-RK2 substep and the complete original-to-final interval
    /// must pass the native positivity, species-continuity and state-minus-bath
    /// entropy checks before publication. One-sided zero species activities
    /// refuse at interior AND bath faces, without a composition floor or an
    /// engineering-flux fallback. Common absence is admitted. The flow clock
    /// is unchanged, and two absent baths retain the closed conjugate faces.
    ///
    /// This is an ideal-mixture, axial molecular closure on fixed geometry,
    /// not a fully discrete stability theorem or advective/reactive boundary
    /// condition. Bath cells carry the original species ordering and energy
    /// references, as required by `advance_with_reservoirs_for`.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_entropy_consistent_with_reservoirs_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(DissipationReport, ReservoirMolecularEntropyReport), DissipationError> {
        self.advance_reservoir_entropy::<true>(
            flow, duration, controls, entropy, reservoirs, cancelled,
        )
    }

    /// Admit the original SSP-RK2 interval with stationary molecular reservoirs.
    ///
    /// Entropy entering each end is q_reduced/T_b + sum_i(s_i,b * j_i/M_i),
    /// with the original partial enthalpies subtracted from the SAME native
    /// energy transfer. Source chemical entropy, reciprocal heat and Fourier
    /// heat are counted once. Actual stage/end species transfers enter the
    /// continuity check; open exchange is not mistaken for a mass defect.
    ///
    /// Both native substeps and the complete interval must pass one declared
    /// state-minus-input allowance. No fields or clock publish on final refusal.
    /// A zero-activity bath species may carry no nonzero molecular flux; its
    /// infinite ideal-mixture chemical potential is never replaced by a floor.
    /// Reservoirs must be stationary and use the original source species/gauges,
    /// as required by advance_with_reservoirs_for. No advective or reactive
    /// entropy is included. Two absent baths also admit periodic closed fields.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_reservoirs_and_entropy_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(DissipationReport, ReservoirMolecularEntropyReport), DissipationError> {
        self.advance_reservoir_entropy::<false>(
            flow, duration, controls, entropy, reservoirs, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_reservoir_entropy<const ENTROPY_FACES: bool>(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DissipationReport, ReservoirMolecularEntropyReport), DissipationError> {
        controls.validate()?;
        entropy.validate()?;
        let reservoirs = if reservoirs.iter().all(Option::is_none) {
            [None; 2]
        } else {
            super::reservoir::admit(&self.gas, self.boundaries, reservoirs, &mut cancelled)?
        };
        let mut account = EntropyTracker::new(flow, entropy, &mut cancelled)?;
        account.boundary = Some(ReservoirEntropyAccount::new(
            &self.gas, account.reference_temperature, &account.reference_density,
            &account.reference_entropy, reservoirs, &mut cancelled,
        )?);
        let report = super::evolve::advance_with_reservoirs_and_entropy::<ENTROPY_FACES, L>(
            self, flow, duration, controls, reservoirs, &mut account, &mut cancelled,
        )?;
        Ok((report, account.completed_reservoir_report
            .expect("native final admission retained reservoir entropy")))
    }

    /// Evolve the existing molecular operator with an additional closed-system
    /// entropy gate. Periodic, stationary adiabatic reflecting, and zero-
    /// molecular-flux transmissive ends are admitted. Prescribed reservoirs,
    /// moving walls, chemistry and advective transport need their own entropy
    /// exchange and are not included by this method.
    ///
    /// Each actual SSP-RK2 candidate must preserve every species inventory and
    /// pass the state-entropy allowance; rejected candidates halve the native
    /// step without clipping species, temperature or energy. The complete
    /// original-to-final account is checked again before caller-state publication.
    /// The flow clock is unchanged. Admitted trials retain the original
    /// molecular fluxes and energy trajectory bit for bit.
    ///
    /// A spatial discretization with negative entropy rate may refuse at every
    /// step size. Refusal does not manufacture an entropy-stable face flux or
    /// grant a fully discrete entropy/invariant-domain theorem to this solver.
    pub fn advance_for_with_entropy(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DissipationReport, MolecularEntropyReport), DissipationError> {
        controls.validate()?;
        entropy.validate()?;
        let mut account = EntropyTracker::new(flow, entropy, &mut cancelled)?;
        let report = super::evolve::advance_with_entropy(
            self,
            flow,
            duration,
            controls,
            &mut account,
            &mut cancelled,
        )?;
        Ok((report, account.report))
    }
}

/// The fixed basis is shared by every actual trial and the final direct check.
pub(super) struct EntropyTracker {
    controls: MolecularEntropyControl,
    reference_temperature: f64,
    reference_density: [f64; MAX_SPECIES],
    reference_entropy: [f64; MAX_SPECIES],
    original_mass: [f64; MAX_SPECIES],
    initial_temperature: Vec<f64>,
    current_temperature: Vec<f64>,
    trial_entropy: Dd,
    trial_mass: [Dd; MAX_SPECIES],
    trial_relative: f64,
    accepted: usize,
    report: MolecularEntropyReport,
    boundary: Option<ReservoirEntropyAccount>,
    completed_reservoir_report: Option<ReservoirMolecularEntropyReport>,
}

impl EntropyTracker {
    pub(super) fn new(
        flow: &MixtureEuler1d,
        controls: MolecularEntropyControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        controls.validate()?;
        poll(cancelled)?;
        match flow.boundaries {
            Boundaries::Periodic => {}
            Boundaries::Ends { left, right }
                if matches!(left, Boundary::Reflecting | Boundary::Transmissive)
                    && matches!(right, Boundary::Reflecting | Boundary::Transmissive) => {}
            _ => {
                return Err(DissipationError::Invalid(
                    "closed molecular entropy requires periodic or stationary zero-heat/species ends",
                ));
            }
        }
        let gas = flow.gas();
        let n = gas.species_count();
        let m = n + 2;
        let mut temperature = allocated(flow.cell_count(), 0.0)?;
        let mut masses = [Dd::ZERO; MAX_SPECIES];
        let mut volume = Dd::ZERO;
        let mut reference_density = [1.0; MAX_SPECIES];
        for (i, q) in flow.states.chunks_exact(m).enumerate() {
            poll(cancelled)?;
            temperature[i] = gas.decode(q, None, cancelled)?.temperature;
            let v = Dd::from_f64(flow.mesh.volumes()[i]);
            volume = volume + v;
            for k in 0..n {
                masses[k] = masses[k] + v * Dd::from_f64(q[k]);
                if q[k] > 0.0 {
                    reference_density[k] = q[k];
                }
            }
        }
        let volume = scalar(volume)?;
        if !positive(volume) {
            return Err(DissipationError::Invalid("closed molecular entropy volume"));
        }
        let reference_temperature = temperature[0];
        let t_ref = Temperature::new(reference_temperature);
        let mut reference_entropy = [0.0; MAX_SPECIES];
        let mut original_mass = [0.0; MAX_SPECIES];
        for (k, model) in gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            original_mass[k] = scalar(masses[k])?;
            let mean = original_mass[k] / volume;
            if positive(mean) {
                reference_density[k] = mean;
            }
            let log_partial_pressure = log_ratio(reference_density[k], model.molar_mass().value())
                + fs_math::det::ln(R)
                + fs_math::det::ln(reference_temperature)
                - fs_math::det::ln(model.convention().reference_pressure().value());
            reference_entropy[k] = model
                .evaluate_properties(t_ref)
                .map_err(super::FlowError::from)?
                .s()
                .value()
                - R * (log_partial_pressure + 1.0);
            if !reference_entropy[k].is_finite() {
                return Err(DissipationError::Invalid(
                    "unrepresentable source entropy reference",
                ));
            }
        }
        let mut current_temperature = allocated(temperature.len(), 0.0)?;
        current_temperature.copy_from_slice(&temperature);
        Ok(Self {
            controls,
            reference_temperature,
            reference_density,
            reference_entropy,
            original_mass,
            initial_temperature: temperature,
            current_temperature,
            trial_entropy: Dd::ZERO,
            trial_mass: [Dd::ZERO; MAX_SPECIES],
            trial_relative: 0.0,
            accepted: 0,
            report: MolecularEntropyReport::default(),
            boundary: None,
            completed_reservoir_report: None,
        })
    }

    pub(super) fn begin_attempt(&mut self, primitive: &[Primitive]) {
        for (i, (temperature, state)) in self
            .current_temperature
            .iter_mut()
            .zip(primitive)
            .enumerate()
        {
            *temperature = state.temperature;
            if self.accepted == 0 {
                self.initial_temperature[i] = state.temperature;
            }
        }
        self.trial_entropy = Dd::ZERO;
        self.trial_mass.fill(Dd::ZERO);
        self.trial_relative = 0.0;
    }

    pub(super) fn candidate_cell(
        &mut self,
        gas: &MixtureGas,
        cell: usize,
        old: &[f64],
        new: &[f64],
        temperature: f64,
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let increment = cell_change(
            gas,
            &self.reference_density,
            self.reference_temperature,
            old,
            new,
            self.current_temperature[cell],
            temperature,
            cancelled,
        )?;
        self.trial_entropy = self.trial_entropy + Dd::from_f64(volume) * increment;
        // Preserve the exact inversion already admitted by the native candidate
        // loop. No second EOS solve selects a different within-tolerance root.
        self.current_temperature[cell] = temperature;
        for k in 0..gas.species_count() {
            self.trial_mass[k] = self.trial_mass[k]
                + Dd::from_f64(volume) * (Dd::from_f64(new[k]) - Dd::from_f64(old[k]));
        }
        Ok(())
    }

    fn inventory(&self, n: usize, whole: bool) -> Result<f64, DissipationError> {
        let mut maximum = 0.0_f64;
        for k in 0..n {
            let transfer = self.boundary.as_ref().map_or(Dd::ZERO, |b| {
                if whole { b.total.species[k] } else { b.trial.species[k] }
            });
            let throughput = self.boundary.as_ref().map_or(Dd::ZERO, |b| {
                if whole { b.total.throughput[k] }
                else { b.total.throughput[k] + b.trial.throughput[k] }
            });
            let residual = scalar(self.trial_mass[k] - transfer)?;
            let scale = scalar(Dd::from_f64(self.original_mass[k]) + throughput)?;
            let tolerance = self.controls.relative_species_tolerance * scale;
            if !tolerance.is_finite() || residual.abs() > tolerance {
                return Err(DissipationError::SpeciesInventoryDrift {
                    species: k, residual_kg: residual, tolerance_kg: tolerance,
                });
            }
            if scale > 0.0 {
                maximum = maximum.max(residual.abs() / scale);
            }
        }
        Ok(maximum)
    }

    pub(super) fn check_attempt(
        &mut self,
        n: usize,
        a: &super::evolve::Stage,
        b: &super::evolve::Stage,
        dt: f64,
        areas: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let check = (|| {
            if let Some(boundary) = &mut self.boundary {
                boundary.stage_pair(a, b, dt, areas, n, cancelled)?;
            }
            self.trial_relative = self.inventory(n, false)?;
            let input = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.trial.input);
            let change = scalar(self.trial_entropy - input)?;
            if change < -self.controls.absolute_tolerance_j_k {
                return Err(DissipationError::EntropyDecrease {
                    change_j_k: change,
                    tolerance_j_k: self.controls.absolute_tolerance_j_k,
                });
            }
            Ok(())
        })();
        if check.is_err() {
            self.report.rejected_entropy_steps += 1;
        }
        check
    }

    pub(super) fn accept_attempt(&mut self, n: usize) -> Result<(), DissipationError> {
        let input = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.trial.input);
        let change = scalar(self.trial_entropy - input)?;
        if let Some(boundary) = &mut self.boundary { boundary.accept(n)?; }
        if self.accepted == 0 {
            self.report.minimum_step_entropy_change_j_k = change;
        } else {
            self.report.minimum_step_entropy_change_j_k =
                self.report.minimum_step_entropy_change_j_k.min(change);
        }
        self.report.maximum_relative_species_residual = self
            .report.maximum_relative_species_residual.max(self.trial_relative);
        self.accepted += 1;
        Ok(())
    }

    fn reservoir_report(&self) -> Result<ReservoirMolecularEntropyReport, DissipationError> {
        let boundary = self.boundary.as_ref().ok_or(DissipationError::Invalid(
            "missing molecular reservoir entropy account",
        ))?;
        let generation = self.trial_entropy - boundary.total.input;
        Ok(ReservoirMolecularEntropyReport {
            source_entropy_change_j_k: self.report.source_entropy_change_j_k,
            boundary_entropy_input_j_k: [
                scalar(boundary.total.source_input[0])?,
                scalar(boundary.total.source_input[1])?,
            ],
            entropy_generation_j_k: scalar(generation)?,
            source_entropy_generation_j_k: scalar(generation
                + Dd::from_f64(self.report.reference_entropy_residual_j_k))?,
            reference_entropy_residual_j_k: self.report.reference_entropy_residual_j_k,
            minimum_step_entropy_generation_j_k: self.report.minimum_step_entropy_change_j_k,
            maximum_relative_species_residual: self.report.maximum_relative_species_residual,
            rejected_entropy_steps: self.report.rejected_entropy_steps,
        })
    }

    /// Recompute the original-to-final state difference, not a sum of accepted
    /// allowances or rounded entropy-density differences. Caller state is intact.
    pub(super) fn finish(
        &mut self,
        gas: &MixtureGas,
        volumes: &[f64],
        original: &[f64],
        final_state: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let n = gas.species_count();
        let m = n + 2;
        self.trial_entropy = Dd::ZERO;
        self.trial_mass.fill(Dd::ZERO);
        for (i, (old, new)) in original
            .chunks_exact(m)
            .zip(final_state.chunks_exact(m))
            .enumerate()
        {
            poll(cancelled)?;
            let temperature = self.current_temperature[i];
            let change = cell_change(
                gas,
                &self.reference_density,
                self.reference_temperature,
                old,
                new,
                self.initial_temperature[i],
                temperature,
                cancelled,
            )?;
            let volume = Dd::from_f64(volumes[i]);
            self.trial_entropy = self.trial_entropy + volume * change;
            for k in 0..n {
                self.trial_mass[k] =
                    self.trial_mass[k] + volume * (Dd::from_f64(new[k]) - Dd::from_f64(old[k]));
            }
        }
        let relative = self.inventory(n, true)?;
        let input = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.total.input);
        let change = scalar(self.trial_entropy - input)?;
        if change < -self.controls.absolute_tolerance_j_k {
            return Err(DissipationError::EntropyDecrease {
                change_j_k: change,
                tolerance_j_k: self.controls.absolute_tolerance_j_k,
            });
        }
        let mut reference = Dd::ZERO;
        let mut reference_residual = Dd::ZERO;
        for (k, model) in gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            reference = reference
                + self.trial_mass[k] / Dd::from_f64(model.molar_mass().value())
                    * Dd::from_f64(self.reference_entropy[k]);
            let transfer = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.total.species[k]);
            reference_residual = reference_residual
                + (self.trial_mass[k] - transfer) / Dd::from_f64(model.molar_mass().value())
                    * Dd::from_f64(self.reference_entropy[k]);
        }
        self.report.entropy_change_j_k = scalar(self.trial_entropy)?;
        self.report.reference_entropy_residual_j_k = scalar(reference_residual)?;
        self.report.source_entropy_change_j_k = scalar(self.trial_entropy + reference)?;
        self.report.maximum_relative_species_residual =
            self.report.maximum_relative_species_residual.max(relative);
        // Every fallible report reduction precedes the native caller-state copy.
        if self.boundary.is_some() {
            self.completed_reservoir_report = Some(self.reservoir_report()?);
        }
        Ok(())
    }
}

fn scalar(value: Dd) -> Result<f64, DissipationError> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(DissipationError::Invalid(
            "non-finite closed molecular entropy",
        ))
    }
}

/// rho1 ln(rho1/rho0) - (rho1-rho0), including its finite zero-density limit.
/// For weak changes the leading quadratic is formed before any subtraction.
fn relative_density_remainder(rho0: f64, rho1: f64) -> Dd {
    if rho1 == 0.0 {
        return Dd::from_f64(rho0);
    }
    let delta = Dd::from_f64(rho1) - Dd::from_f64(rho0);
    let z = delta.to_f64() / rho0;
    if z.abs() <= 0.125 {
        let mut term = 0.5 * z * z;
        let mut result = Dd::from_f64(term);
        for k in 2..24 {
            term *= -z * (k - 1) as f64 / (k + 1) as f64;
            result = result + Dd::from_f64(term);
            if term.abs() <= 0.25 * f64::EPSILON * result.to_f64().abs() {
                break;
            }
        }
        Dd::from_f64(rho0) * result
    } else {
        Dd::from_f64(rho1) * Dd::from_f64(log_ratio(rho1, rho0)) - delta
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn cell_change(
    gas: &MixtureGas,
    reference_density: &[f64; MAX_SPECIES],
    reference_temperature: f64,
    old: &[f64],
    new: &[f64],
    old_temperature: f64,
    new_temperature: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, DissipationError> {
    let mut entropy = Dd::ZERO;
    let reference = Temperature::new(reference_temperature);
    let initial = Temperature::new(old_temperature);
    let final_t = Temperature::new(new_temperature);
    let log_t = log_ratio(new_temperature, old_temperature);
    for (k, model) in gas.source().models().iter().enumerate() {
        poll(cancelled)?;
        let rho0 = old[k];
        let rho1 = new[k];
        if (rho0 == rho1 && old_temperature == new_temperature) || (rho0 == 0.0 && rho1 == 0.0) {
            continue;
        }
        let delta = Dd::from_f64(rho1) - Dd::from_f64(rho0);
        let mut change = Dd::ZERO;
        if rho1 > 0.0 && old_temperature != new_temperature {
            let ds = model
                .molar_entropy_change(initial, final_t)
                .map_err(super::FlowError::from)?
                .value();
            change =
                Dd::from_f64(rho1) * (Dd::from_f64(ds) - Dd::from_f64(R) * Dd::from_f64(log_t));
        }
        if rho0 != rho1 {
            let ds0 = model
                .molar_entropy_change(reference, initial)
                .map_err(super::FlowError::from)?
                .value();
            let phi0 = Dd::from_f64(ds0)
                - Dd::from_f64(R) * Dd::from_f64(log_ratio(old_temperature, reference_temperature));
            let mixing = if rho0 == 0.0 {
                Dd::from_f64(rho1) * Dd::from_f64(log_ratio(rho1, reference_density[k]) - 1.0)
            } else {
                delta * Dd::from_f64(log_ratio(rho0, reference_density[k]))
                    + relative_density_remainder(rho0, rho1)
            };
            change = change + delta * phi0 - Dd::from_f64(R) * mixing;
        }
        entropy = entropy + change / Dd::from_f64(model.molar_mass().value());
    }
    scalar(entropy)?;
    Ok(entropy)
}
