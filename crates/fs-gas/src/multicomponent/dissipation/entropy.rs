//! Original-source state entropy shared by axial gas transport owners.
//!
//! Molecular SSP-RK2/RKL2 and the inviscid step use the same source intervals,
//! reference basis, species continuity checks and final report arithmetic.
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
use super::super::{
    CaloricDomain, Cell, Error as FlowError, MAX_SPECIES, MixtureEuler1d, MixtureGas, Primitive,
    allocated, poll,
};
use super::{
    Boundaries, Boundary, DissipationControl, DissipationError, DissipationOperator,
    DissipationReport, TransportLaw, positive,
};
use fs_math::dd::Dd;
use fs_qty::{Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;

mod advection;
mod reservoir;
pub use advection::ExteriorMixtureEntropyReport;
pub(in crate::multicomponent) use advection::ExteriorTransportInventory;
pub(in crate::multicomponent) mod reactive;
mod stabilized;
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
    fn validate(self) -> Result<(), FlowError> {
        if !self.absolute_tolerance_j_k.is_finite()
            || self.absolute_tolerance_j_k < 0.0
            || !self.relative_species_tolerance.is_finite()
            || !(0.0..=1e-3).contains(&self.relative_species_tolerance)
        {
            return Err(FlowError::EntropyDomain(
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
    /// Smallest admitted physical substep or RKL2 macrostep entropy change [J/K].
    /// Zero for zero duration; internal Legendre stages are not physical steps.
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

    /// Original source-conjugate SSP owner with its actual complete boundary
    /// inventory for source compositions. The caller receives no reconstructed
    /// entropy or approximate transfer in place of the accepted native axes.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent) fn advance_entropy_consistent_with_reservoir_inventory_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(
        DissipationReport,
        ReservoirMolecularEntropyReport,
        ExteriorTransportInventory,
    ), DissipationError> {
        self.advance_reservoir_entropy_with_inventory::<true>(
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<(DissipationReport, ReservoirMolecularEntropyReport), DissipationError> {
        let (report, measured, _) = self.advance_reservoir_entropy_with_inventory::<ENTROPY_FACES>(
            flow, duration, controls, entropy, reservoirs, cancelled,
        )?;
        Ok((report, measured))
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_reservoir_entropy_with_inventory<const ENTROPY_FACES: bool>(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(
        DissipationReport,
        ReservoirMolecularEntropyReport,
        ExteriorTransportInventory,
    ), DissipationError> {
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
        // All physical/source checks precede the original publication above.
        // Export only already-admitted native transfers and stage throughput.
        let boundary = account.boundary.as_ref()
            .expect("native reservoir interval retained its admitted bath account");
        let inventory = ExteriorTransportInventory::from_molecular_report(
            &report, boundary.total.throughput,
        );
        Ok((report, account.completed_reservoir_report
            .expect("native final admission retained reservoir entropy"), inventory))
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
pub(in crate::multicomponent) struct EntropyTracker {
    controls: MolecularEntropyControl,
    // The original source owns immutable shared models. Retaining it verifies
    // borrowed basis identity without another field clone or caloric query.
    source_owner: IdealGasEquilibriumV1,
    source_domain: CaloricDomain,
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
    // One pair per accounted native branch. Keep cached exterior efforts off
    // the closed-lane stack; open networks retain their common source basis.
    exterior: Vec<advection::ExteriorEntropyAccount>,
    completed_exterior_report: Option<ExteriorMixtureEntropyReport>,
}

impl EntropyTracker {
    pub(in crate::multicomponent) fn new(
        flow: &MixtureEuler1d,
        controls: MolecularEntropyControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        Self::new_with_basis(flow, controls, None, cancelled)
    }

    /// Keep fresh LOCAL initial temperatures and inventories while optionally
    /// sharing an immutable complete-window source/reference basis. The borrowed
    /// basis changes no state, density, temperature, source or numerical allowance.
    pub(in crate::multicomponent) fn new_with_basis(
        flow: &MixtureEuler1d,
        controls: MolecularEntropyControl,
        basis: Option<&Self>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        controls.validate()?;
        poll(cancelled)?;
        match flow.boundaries {
            Boundaries::Periodic => {}
            Boundaries::Ends { left, right }
                if matches!(left, Boundary::Reflecting | Boundary::Transmissive)
                    && matches!(right, Boundary::Reflecting | Boundary::Transmissive) => {}
            _ => {
                return Err(FlowError::EntropyDomain(
                    "closed molecular entropy requires periodic or stationary zero-heat/species ends",
                ));
            }
        }
        let gas = flow.gas();
        let m = gas.species_count() + 2;
        Self::from_cells(
            gas,
            flow.cell_count(),
            controls,
            basis,
            flow.states
                .chunks_exact(m)
                .zip(flow.mesh.volumes())
                .map(|(q, &volume)| (super::super::cell(q, gas.species_count()), volume)),
            cancelled,
        )
    }

    /// One source/reference basis for heterogeneous native storage owners.
    /// Cells retain their actual local density and physical volume; no flattened
    /// state copy, composition averaging or extra thermodynamic model is needed.
    pub(in crate::multicomponent) fn from_cells<'a>(
        gas: &MixtureGas,
        count: usize,
        controls: MolecularEntropyControl,
        basis: Option<&Self>,
        cells: impl IntoIterator<Item = (Cell<'a>, f64)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        controls.validate()?;
        poll(cancelled)?;
        if !(1..=crate::dynamics::MAX_EULER_CELLS).contains(&count) {
            return Err(FlowError::EntropyDomain("entropy storage count"));
        }
        if basis.is_some_and(|basis| !basis.matches_source(gas)) {
            return Err(FlowError::EntropyDomain(
                "entropy basis source or caloric domain mismatch",
            ));
        }
        let n = gas.species_count();
        let mut temperature = allocated(count, 0.0)?;
        let mut masses = [Dd::ZERO; MAX_SPECIES];
        let mut volume = Dd::ZERO;
        let mut reference_density = [1.0; MAX_SPECIES];
        let mut seen = 0;
        for (i, (q, cell_volume)) in cells.into_iter().enumerate() {
            poll(cancelled)?;
            if i >= count || !positive(cell_volume) {
                return Err(FlowError::EntropyDomain(
                    "entropy storage count or cell volume",
                ));
            }
            temperature[i] = gas.decode_cell(q, None, cancelled)?.temperature;
            let v = Dd::from_f64(cell_volume);
            volume = volume + v;
            for k in 0..n {
                masses[k] = masses[k] + v * Dd::from_f64(q.partial_densities[k]);
                if q.partial_densities[k] > 0.0 {
                    reference_density[k] = q.partial_densities[k];
                }
            }
            seen += 1;
        }
        if seen != count {
            return Err(FlowError::EntropyDomain("entropy storage count"));
        }
        let volume = scalar(volume)?;
        if !positive(volume) {
            return Err(FlowError::EntropyDomain("closed molecular entropy volume"));
        }
        let reference_temperature = basis.map_or(temperature[0], |basis| basis.reference_temperature);
        let t_ref = Temperature::new(reference_temperature);
        let mut reference_entropy = [0.0; MAX_SPECIES];
        let mut original_mass = [0.0; MAX_SPECIES];
        if let Some(basis) = basis {
            reference_density = basis.reference_density;
            reference_entropy = basis.reference_entropy;
            for k in 0..n {
                poll(cancelled)?;
                original_mass[k] = scalar(masses[k])?;
            }
        } else {
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
                    return Err(FlowError::EntropyDomain(
                        "unrepresentable source entropy reference",
                    ));
                }
            }
        }
        let mut current_temperature = allocated(temperature.len(), 0.0)?;
        current_temperature.copy_from_slice(&temperature);
        Ok(Self {
            controls,
            source_owner: gas.source().clone(),
            source_domain: gas.domain(),
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
            exterior: Vec::new(),
            completed_exterior_report: None,
        })
    }

    fn matches_source(&self, gas: &MixtureGas) -> bool {
        self.source_domain == gas.domain() && &self.source_owner == gas.source()
    }

    /// A normalized transfer can be summed directly only in the identical
    /// source chart and bit-exact fixed reference basis. No report rebasing occurs.
    fn same_basis(&self, other: &Self) -> bool {
        self.source_domain == other.source_domain
            && self.source_owner == other.source_owner
            && self.reference_temperature.to_bits() == other.reference_temperature.to_bits()
            && self.reference_density.iter().zip(&other.reference_density)
                .all(|(a, b)| a.to_bits() == b.to_bits())
            && self.reference_entropy.iter().zip(&other.reference_entropy)
                .all(|(a, b)| a.to_bits() == b.to_bits())
    }

    pub(in crate::multicomponent) fn begin_attempt(&mut self, primitive: &[Primitive]) {
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

    /// Retry a complete closed native step from its original accepted storage.
    /// Network primitive scratch may contain a rejected stage; initial admitted
    /// temperatures remain the entropy anchors until publication.
    pub(in crate::multicomponent) fn begin_network_attempt(&mut self) {
        self.current_temperature
            .copy_from_slice(&self.initial_temperature);
        self.trial_entropy = Dd::ZERO;
        self.trial_mass.fill(Dd::ZERO);
        self.trial_relative = 0.0;
    }

    pub(in crate::multicomponent) fn candidate_cell(
        &mut self,
        gas: &MixtureGas,
        cell: usize,
        old: &[f64],
        new: &[f64],
        temperature: f64,
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
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

    /// Account actual extensive storage during prescribed chamber-volume motion.
    /// Mechanical pressure work changes the caloric state but carries no entropy.
    /// The unchanged-volume path keeps the existing field arithmetic exactly.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent) fn candidate_moving_cell(
        &mut self,
        gas: &MixtureGas,
        cell: usize,
        old: &[f64],
        new: &[f64],
        temperature: f64,
        old_volume: f64,
        new_volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        if old_volume == new_volume {
            return self.candidate_cell(gas, cell, old, new, temperature, old_volume, cancelled);
        }
        let increment = moving_cell_change(
            gas,
            &self.reference_density,
            self.reference_temperature,
            old,
            new,
            self.current_temperature[cell],
            temperature,
            old_volume,
            new_volume,
            cancelled,
        )?;
        self.trial_entropy = self.trial_entropy + increment;
        self.current_temperature[cell] = temperature;
        for k in 0..gas.species_count() {
            self.trial_mass[k] = self.trial_mass[k]
                + (Dd::from_f64(new_volume) * Dd::from_f64(new[k])
                    - Dd::from_f64(old_volume) * Dd::from_f64(old[k]));
        }
        Ok(())
    }

    fn inventory(&self, n: usize, whole: bool) -> Result<f64, FlowError> {
        self.inventory_with_reaction(n, whole, None)
    }

    fn inventory_with_reaction(
        &self,
        n: usize,
        whole: bool,
        reaction: Option<(&[Dd; MAX_SPECIES], &[Dd; MAX_SPECIES])>,
    ) -> Result<f64, FlowError> {
        let mut maximum = 0.0_f64;
        for k in 0..n {
            let mut transfer = self.boundary.as_ref().map_or(Dd::ZERO, |b| {
                if whole { b.total.species[k] } else { b.trial.species[k] }
            });
            let mut throughput = self.boundary.as_ref().map_or(Dd::ZERO, |b| {
                if whole { b.total.throughput[k] }
                else { b.total.throughput[k] + b.trial.throughput[k] }
            });
            if let Some((change, absolute)) = reaction {
                transfer = transfer + change[k];
                throughput = throughput + absolute[k];
            }
            let residual = scalar(self.trial_mass[k] - transfer)?;
            let scale = scalar(Dd::from_f64(self.original_mass[k]) + throughput)?;
            let tolerance = self.controls.relative_species_tolerance * scale;
            if !tolerance.is_finite() || residual.abs() > tolerance {
                return Err(FlowError::SpeciesInventoryDrift {
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
    ) -> Result<(), FlowError> {
        let check = (|| {
            if let Some(boundary) = &mut self.boundary {
                boundary.stage_pair(a, b, dt, areas, n, cancelled)?;
            }
            self.trial_relative = self.inventory(n, false)?;
            let input = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.trial.input);
            let change = scalar(self.trial_entropy - input)?;
            if change < -self.controls.absolute_tolerance_j_k {
                return Err(FlowError::EntropyDecrease {
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

    pub(super) fn accept_attempt(&mut self, n: usize) -> Result<(), FlowError> {
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

    fn reservoir_report(
        &self,
        report: MolecularEntropyReport,
    ) -> Result<ReservoirMolecularEntropyReport, FlowError> {
        let boundary = self.boundary.as_ref().ok_or(FlowError::EntropyDomain(
            "missing molecular reservoir entropy account",
        ))?;
        let generation = self.trial_entropy - boundary.total.input;
        Ok(ReservoirMolecularEntropyReport {
            source_entropy_change_j_k: report.source_entropy_change_j_k,
            boundary_entropy_input_j_k: [
                scalar(boundary.total.source_input[0])?,
                scalar(boundary.total.source_input[1])?,
            ],
            entropy_generation_j_k: scalar(generation)?,
            source_entropy_generation_j_k: scalar(generation
                + Dd::from_f64(report.reference_entropy_residual_j_k))?,
            reference_entropy_residual_j_k: report.reference_entropy_residual_j_k,
            minimum_step_entropy_generation_j_k: report.minimum_step_entropy_change_j_k,
            maximum_relative_species_residual: report.maximum_relative_species_residual,
            rejected_entropy_steps: report.rejected_entropy_steps,
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
    ) -> Result<(), FlowError> {
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
        let report = self.complete_report(gas, cancelled)?;
        let reservoir_report = if self.boundary.is_some() {
            Some(self.reservoir_report(report)?)
        } else {
            None
        };
        // Every fallible report reduction precedes both account and field publication.
        self.report = report;
        self.completed_reservoir_report = reservoir_report;
        Ok(())
    }

    /// Admit one closed native step using the original-to-candidate increments
    /// already collected by candidate_cell. No second source interval or EOS
    /// evaluation is performed. Failed report arithmetic cannot accept the trial.
    pub(in crate::multicomponent) fn admit_single_step(
        &mut self,
        gas: &MixtureGas,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        if self.boundary.is_some() || self.accepted != 0 {
            return Err(FlowError::EntropyDomain(
                "single-step entropy admission requires a fresh closed account",
            ));
        }
        let completion = self.complete_report(gas, cancelled);
        match completion {
            Ok(mut report) => {
                report.minimum_step_entropy_change_j_k = report.entropy_change_j_k;
                self.report = report;
                self.accepted = 1;
                Ok(())
            }
            Err(error) => {
                if matches!(
                    &error,
                    FlowError::EntropyDecrease { .. } | FlowError::SpeciesInventoryDrift { .. }
                ) {
                    self.report.rejected_entropy_steps += 1;
                }
                Err(error)
            }
        }
    }

    /// The report is copied only after the native caller has admitted its step.
    pub(in crate::multicomponent) fn report(&self) -> MolecularEntropyReport {
        self.report
    }

    /// Shared original-to-final reduction for one native step and a complete
    /// molecular interval. The input increments are already in the fixed basis;
    /// all reference and continuity reductions complete before changing a report.
    fn complete_report(
        &self,
        gas: &MixtureGas,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularEntropyReport, FlowError> {
        let relative = self.inventory(gas.species_count(), true)?;
        let input = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.total.input);
        let change = scalar(self.trial_entropy - input)?;
        if change < -self.controls.absolute_tolerance_j_k {
            return Err(FlowError::EntropyDecrease {
                change_j_k: change,
                tolerance_j_k: self.controls.absolute_tolerance_j_k,
            });
        }
        let (reference, reference_residual, _) =
            self.reference_changes(gas, None, cancelled)?;
        let mut report = self.report;
        report.entropy_change_j_k = scalar(self.trial_entropy)?;
        report.reference_entropy_residual_j_k = scalar(reference_residual)?;
        report.source_entropy_change_j_k = scalar(self.trial_entropy + reference)?;
        report.maximum_relative_species_residual =
            report.maximum_relative_species_residual.max(relative);
        Ok(report)
    }

    /// Original source constants on the total change, continuity defect and
    /// actual chemical source, respectively. A chemical species change is
    /// physical and must not be removed with the numerical continuity defect.
    fn reference_changes(
        &self,
        gas: &MixtureGas,
        reaction: Option<&[Dd; MAX_SPECIES]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Dd, Dd, Dd), FlowError> {
        self.reference_changes_with_transfer(gas, reaction, None, cancelled)
    }

    /// The same original source basis with an additional actual boundary
    /// transfer. Subtract transfers at the species level before applying alpha;
    /// genuine chemical source entropy remains a separate retained quantity.
    fn reference_changes_with_transfer(
        &self,
        gas: &MixtureGas,
        reaction: Option<&[Dd; MAX_SPECIES]>,
        boundary_transfer: Option<&[Dd; MAX_SPECIES]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Dd, Dd, Dd), FlowError> {
        let mut reference = Dd::ZERO;
        let mut reference_residual = Dd::ZERO;
        let mut reaction_reference = Dd::ZERO;
        for (k, model) in gas.source().models().iter().enumerate() {
            poll(cancelled)?;
            reference = reference
                + self.trial_mass[k] / Dd::from_f64(model.molar_mass().value())
                    * Dd::from_f64(self.reference_entropy[k]);
            let transfer = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.total.species[k]);
            let mut residual = self.trial_mass[k] - transfer;
            if let Some(boundary) = boundary_transfer {
                residual = residual - boundary[k];
            }
            if let Some(change) = reaction {
                residual = residual - change[k];
                reaction_reference = reaction_reference
                    + change[k] / Dd::from_f64(model.molar_mass().value())
                        * Dd::from_f64(self.reference_entropy[k]);
            }
            reference_residual = reference_residual
                + residual / Dd::from_f64(model.molar_mass().value())
                    * Dd::from_f64(self.reference_entropy[k]);
        }
        Ok((reference, reference_residual, reaction_reference))
    }
}

fn scalar(value: Dd) -> Result<f64, FlowError> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(FlowError::EntropyDomain(
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
) -> Result<Dd, FlowError> {
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

/// Extensive version of the same original-mixture entropy primitive, with
/// actual endpoint volumes. Write m=V*rho before differencing; rescaling a final
/// density into the old volume can invent a species-roundoff contribution.
#[allow(clippy::too_many_arguments)]
fn moving_cell_change(
    gas: &MixtureGas,
    reference_density: &[f64; MAX_SPECIES],
    reference_temperature: f64,
    old: &[f64],
    new: &[f64],
    old_temperature: f64,
    new_temperature: f64,
    old_volume: f64,
    new_volume: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, FlowError> {
    let mut entropy = Dd::ZERO;
    let reference = Temperature::new(reference_temperature);
    let initial = Temperature::new(old_temperature);
    let final_t = Temperature::new(new_temperature);
    let log_t = log_ratio(new_temperature, old_temperature);
    let v0 = Dd::from_f64(old_volume);
    let v1 = Dd::from_f64(new_volume);
    for (k, model) in gas.source().models().iter().enumerate() {
        poll(cancelled)?;
        let (rho0, rho1) = (old[k], new[k]);
        if rho0 == 0.0 && rho1 == 0.0 {
            continue;
        }
        let m0 = v0 * Dd::from_f64(rho0);
        let m1 = v1 * Dd::from_f64(rho1);
        let delta = m1 - m0;
        let mut change = Dd::ZERO;
        if rho1 > 0.0 && old_temperature != new_temperature {
            let ds = model
                .molar_entropy_change(initial, final_t)
                .map_err(super::FlowError::from)?
                .value();
            change = m1 * (Dd::from_f64(ds) - Dd::from_f64(R) * Dd::from_f64(log_t));
        }
        let ds0 = model
            .molar_entropy_change(reference, initial)
            .map_err(super::FlowError::from)?
            .value();
        let phi0 = Dd::from_f64(ds0)
            - Dd::from_f64(R) * Dd::from_f64(log_ratio(old_temperature, reference_temperature));
        // Delta[m*(ln(rho/rho*)-1)], retaining its quadratic weak-change
        // remainder and the finite zero-population limit. The last term is
        // essential: expansion changes density without changing species mass.
        let mixing = if rho0 == 0.0 {
            m1 * Dd::from_f64(log_ratio(rho1, reference_density[k]) - 1.0)
        } else {
            delta * Dd::from_f64(log_ratio(rho0, reference_density[k]))
                + v1 * relative_density_remainder(rho0, rho1)
                - (v1 - v0) * Dd::from_f64(rho0)
        };
        change = change + delta * phi0 - Dd::from_f64(R) * mixing;
        entropy = entropy + change / Dd::from_f64(model.molar_mass().value());
    }
    scalar(entropy)?;
    Ok(entropy)
}
