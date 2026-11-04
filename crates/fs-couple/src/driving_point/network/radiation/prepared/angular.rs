//! Loaded generator coordinates reach the original positive BEM angular owner.

use super::*;
use fs_bem::observation_transfers::RadiatedPowerMoments;

impl<'a> PreparedRadiatingNetwork<'a> {
    /// Available only after `prepare_sources_with_radiation`. The returned view
    /// borrows this exact loaded network and its retained opening angular bank.
    pub fn far_field(&self) -> Option<PreparedRadiatingFarField<'_, 'a>> {
        self.radiation
            .as_ref()
            .map(|radiation| PreparedRadiatingFarField {
                owner: self,
                radiation,
            })
    }
}

/// Allocation-free outgoing observations of the original loaded source bank.
/// Inputs retain canonical open-circuit pressure [Pa] or graph-inflow [m³/s]
/// coordinates; they are mapped into actual outward opening flows before any
/// field products. Independent generators generally excite correlated openings.
///
/// All amplitudes use peak exp(-iωt). Far-field F [Pa m] means p ~ F exp(ikr)/r;
/// dP/dΩ = E[|F|²]/(2ρc) [W/sr]. Returned full-sphere powers/moments are line
/// quantities, not per-Hz densities. They observe the existing radiation load
/// and never add a second debit to network work. Angular convergence does not
/// establish mesh/port-basis convergence or replace the signed surface-work
/// diagnostic. Discard the caller's output slice after any error/cancellation;
/// no converged power result is returned on failure.
pub struct PreparedRadiatingFarField<'p, 'a> {
    owner: &'p PreparedRadiatingNetwork<'a>,
    radiation: &'p PreparedRadiationObservation,
}

impl PreparedRadiatingFarField<'_, '_> {
    /// Canonical physical source labels and their original, distinct SI units.
    pub fn source_coordinates(&self) -> &[NetworkSourceCoordinate] {
        self.owner.source_coordinates()
    }

    /// Requested unit directions, independent of the power quadrature ladder.
    pub fn directions(&self) -> &[[f64; 3]] {
        self.radiation.directions()
    }

    /// Actual angular kernel work performed once during cold preparation.
    pub fn kernel_evaluations(&self) -> usize {
        self.radiation.kernel_evaluations()
    }

    /// Full cached positive quadrature ladder, before drive-specific convergence.
    pub fn cached_angular_directions(&self) -> usize {
        self.radiation.cached_angular_directions()
    }

    /// Measured generator bin dP/dΩ [W/sr] at the original requested directions.
    /// Source-to-opening composition was admitted and performed once by the
    /// measured bin; this runtime allowance bounds the original angular owner.
    pub fn observe_welch_bin(
        &self,
        measured: &PreparedRadiatingWelchBin<'_, '_>,
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        if !std::ptr::eq(self.owner, measured.owner) {
            return Err(bad("measured angular bin belongs to another preparation"));
        }
        poll(gate)?;
        Ok(self.radiation.observe_coherence_modes(
            &measured.opening_factors,
            measured.segments,
            budget,
            output_w_per_sr,
            || gate.is_requested(),
        )?)
    }

    /// Measured directional powers and total outgoing power/moments [W].
    /// All original segment factors enter the SAME positive quadrature ladder;
    /// convergence is tested against this actual measured covariance. Summing
    /// separately converged segment powers could miss a different required
    /// angular resolution. This observation adds no radiation-work debit.
    pub fn evaluate_welch_bin(
        &self,
        measured: &PreparedRadiatingWelchBin<'_, '_>,
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        gate: &CancelGate,
    ) -> Result<RadiatedPowerMoments, ExteriorNetworkError> {
        if !std::ptr::eq(self.owner, measured.owner) {
            return Err(bad("measured angular bin belongs to another preparation"));
        }
        poll(gate)?;
        Ok(self.radiation.evaluate_coherence_modes(
            &measured.opening_factors,
            measured.segments,
            budget,
            output_w_per_sr,
            || gate.is_requested(),
        )?)
    }

    /// Requested-direction coherent amplitudes only; no power integration.
    pub fn observe_coherent(
        &self,
        excitation: &[C64],
        budget: HarmonicWorkBudget,
        output: &mut [C64],
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        let budget = self.remaining_budget(1, false, budget)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .opening_factors(excitation, 1, &mut scratch, gate)?;
        Ok(self
            .radiation
            .observe_coherent(openings, budget, output, || gate.is_requested())?)
    }

    /// Peak F at requested directions plus full positive radiated power/moments.
    /// Two successive quadrature refinements must agree for this actual drive.
    pub fn evaluate_coherent(
        &self,
        excitation: &[C64],
        budget: HarmonicWorkBudget,
        output: &mut [C64],
        gate: &CancelGate,
    ) -> Result<RadiatedPowerMoments, ExteriorNetworkError> {
        let budget = self.remaining_budget(1, false, budget)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .opening_factors(excitation, 1, &mut scratch, gate)?;
        Ok(self
            .radiation
            .evaluate_coherent(openings, budget, output, || gate.is_requested())?)
    }

    /// Requested dP/dΩ for factor-major C=B Bᴴ in canonical mixed source units.
    /// Sources within a factor interfere; distinct factors contribute powers.
    pub fn observe_coherence_modes(
        &self,
        factors: &[C64],
        rank: usize,
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        let budget = self.remaining_budget(rank, false, budget)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .opening_factors(factors, rank, &mut scratch, gate)?;
        Ok(self.radiation.observe_coherence_modes(
            openings,
            rank,
            budget,
            output_w_per_sr,
            || gate.is_requested(),
        )?)
    }

    /// Partial-coherence dP/dΩ and total power/moments; zero factors are silence.
    pub fn evaluate_coherence_modes(
        &self,
        factors: &[C64],
        rank: usize,
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        gate: &CancelGate,
    ) -> Result<RadiatedPowerMoments, ExteriorNetworkError> {
        let budget = self.remaining_budget(rank, false, budget)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .opening_factors(factors, rank, &mut scratch, gate)?;
        Ok(self.radiation.evaluate_coherence_modes(
            openings,
            rank,
            budget,
            output_w_per_sr,
            || gate.is_requested(),
        )?)
    }

    /// Independent canonical generator second moments E[|x_i|²], in squared SI
    /// units. Each generator's correlated opening field is summed before |F|².
    pub fn observe_independent_sources(
        &self,
        second_moments: &[f64],
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        gate: &CancelGate,
    ) -> Result<(), ExteriorNetworkError> {
        let rank = self.source_coordinates().len();
        let budget = self.remaining_budget(rank, true, budget)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .independent_openings(second_moments, &mut scratch, gate)?;
        Ok(self.radiation.observe_coherence_modes(
            openings,
            rank,
            budget,
            output_w_per_sr,
            || gate.is_requested(),
        )?)
    }

    /// Independent generators' dP/dΩ and total power/moments. Using diagonal
    /// opening powers here would discard their physically induced coherence.
    pub fn evaluate_independent_sources(
        &self,
        second_moments: &[f64],
        budget: HarmonicWorkBudget,
        output_w_per_sr: &mut [f64],
        gate: &CancelGate,
    ) -> Result<RadiatedPowerMoments, ExteriorNetworkError> {
        let rank = self.source_coordinates().len();
        let budget = self.remaining_budget(rank, true, budget)?;
        let mut scratch = [C64::ZERO; FACTOR_ENTRIES];
        let openings = self
            .owner
            .independent_openings(second_moments, &mut scratch, gate)?;
        Ok(self.radiation.evaluate_coherence_modes(
            openings,
            rank,
            budget,
            output_w_per_sr,
            || gate.is_requested(),
        )?)
    }

    // The caller's one runtime allowance covers BOTH loaded factor mapping and
    // the original BEM angular contraction. Storage is fixed at <=32² complex
    // values; mapping polls between factors and costs at most 32³ products.
    fn remaining_budget(
        &self,
        rank: usize,
        independent: bool,
        budget: HarmonicWorkBudget,
    ) -> Result<HarmonicWorkBudget, ExteriorNetworkError> {
        if rank > MAX_NETWORK_PORTS {
            return Err(bad("radiating source coherence rank exceeds 32"));
        }
        let sources = self.source_coordinates().len();
        let products =
            self.owner.opening_indices.len() * if independent { sources } else { rank * sources };
        let maximum_products = budget
            .maximum_products
            .checked_sub(products)
            .ok_or_else(|| bad("loaded source mapping exceeds angular work budget"))?;
        Ok(HarmonicWorkBudget { maximum_products })
    }
}
