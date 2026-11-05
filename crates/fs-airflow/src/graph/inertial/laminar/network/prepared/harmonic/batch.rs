//! Atomic multi-excitation reuse of one admitted complex factorization.

use super::*;

impl<M: FluidTransport> LaminarNetwork<M> {
    /// Solve several actual pressure excitations at one frequency while reusing
    /// ONE native complex factorization. Every case is preflighted before the
    /// factorization; every solved peak flow, local heat/storage account and
    /// continuity identity is checked before any response vector is returned.
    /// A refusal discards the complete batch. The caller-owned case slice is the
    /// explicit batch-size and result-allocation bound.
    pub fn harmonic_responses_batch(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressure_cases_pa: &[Vec<C64>],
        budget: HarmonicBudget,
    ) -> Result<Vec<HarmonicResponse>, HarmonicError> {
        self.harmonic_responses_batch_impl(
            cx,
            omega_rad_s,
            reservoir_pressure_cases_pa,
            None,
            budget,
        )
    }

    /// Batched checked response with one explicit branch-head vector per case.
    /// The two case arrays must have identical nonzero length.
    pub fn harmonic_responses_with_heads_batch(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressure_cases_pa: &[Vec<C64>],
        branch_head_cases_pa: &[Vec<C64>],
        budget: HarmonicBudget,
    ) -> Result<Vec<HarmonicResponse>, HarmonicError> {
        self.harmonic_responses_batch_impl(
            cx,
            omega_rad_s,
            reservoir_pressure_cases_pa,
            Some(branch_head_cases_pa),
            budget,
        )
    }

    fn harmonic_responses_batch_impl(
        &self,
        cx: &Cx<'_>,
        omega_rad_s: f64,
        reservoir_pressure_cases_pa: &[Vec<C64>],
        branch_head_cases_pa: Option<&[Vec<C64>]>,
        budget: HarmonicBudget,
    ) -> Result<Vec<HarmonicResponse>, HarmonicError> {
        poll(cx)?;
        if reservoir_pressure_cases_pa.is_empty()
            || branch_head_cases_pa
                .is_some_and(|cases| cases.len() != reservoir_pressure_cases_pa.len())
        {
            return Err(invalid(
                "harmonic batch needs matching nonempty pressure and head cases",
            ));
        }
        for case in 0..reservoir_pressure_cases_pa.len() {
            let heads = branch_head_cases_pa.map(|cases| cases[case].as_slice());
            self.admit_harmonic_input(&reservoir_pressure_cases_pa[case], heads)?;
            poll(cx)?;
        }
        self.admit_harmonic_frequency(cx, omega_rad_s)?;
        let ports = self.port_count();
        let tolerance = harmonic_tolerance(self, ports, budget)?;
        let model = self.harmonic_model(cx, budget)?;
        let factor = model.factor(omega_rad_s, || cx.checkpoint().is_err())?;
        let mut responses = Vec::with_capacity(reservoir_pressure_cases_pa.len());
        for case in 0..reservoir_pressure_cases_pa.len() {
            poll(cx)?;
            let heads = branch_head_cases_pa.map(|cases| cases[case].as_slice());
            let input =
                self.harmonic_input_phasors(&reservoir_pressure_cases_pa[case], heads)?;
            responses.push(self.response_with_factor(cx, &factor, &input, tolerance)?);
        }
        poll(cx)?;
        Ok(responses)
    }
}
