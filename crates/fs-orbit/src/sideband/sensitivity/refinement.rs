//! Joint resolution of a response and its exact carrier-frequency derivative.

use super::super::refinement::{
    EquationAssessment, SidebandRefinementCheck, assess_equations, state_change,
};
use super::super::{
    HbOrbit, PreparedPeriodicLinearization, Result, SidebandControls, SidebandError,
    SidebandScales, admit, bound, finite, invalid, poll, positive, resource_costs,
};
use super::{SidebandCarrierDerivativeControls, SidebandCarrierDerivativeResponse};
use crate::OrbitProblem;
use fs_math::c64::C64;

/// Independent physical allowances for the original and differentiated equations.
#[derive(Clone, Debug)]
pub struct SidebandCarrierRefinementTolerances {
    /// Complex-envelope L2 balance allowance per original equation row.
    pub equation_tolerances: Vec<f64>,
    /// Complex-envelope L2 change per original state between actual solves.
    pub state_tolerances: Vec<f64>,
    /// Original equation units PER angular frequency, including P_nu X - F_nu.
    pub derivative_equation_tolerances: Vec<f64>,
    /// Original state units PER angular frequency for actual inter-order change.
    pub derivative_state_tolerances: Vec<f64>,
}

impl SidebandCarrierRefinementTolerances {
    fn admit(&self, dimension: usize) -> Result<()> {
        let vectors = [
            &self.equation_tolerances,
            &self.state_tolerances,
            &self.derivative_equation_tolerances,
            &self.derivative_state_tolerances,
        ];
        if vectors
            .iter()
            .any(|v| v.len() != dimension || v.iter().any(|&x| !positive(x)))
        {
            return Err(invalid(
                "positive physical response and derivative tolerances for every original row/state required",
            ));
        }
        Ok(())
    }
}

/// Joint adaptive retained-order controls, independent of numerical chart scales.
#[derive(Clone, Debug)]
pub struct SidebandCarrierRefinementControls {
    /// Largest attempted order, strictly above the starting order and at most 64.
    pub maximum_harmonics: usize,
    pub tolerances: SidebandCarrierRefinementTolerances,
}

/// Full original-equation measurements of BOTH fields at one retained order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SidebandCarrierRefinementCheck {
    pub response: SidebandRefinementCheck,
    pub angular_derivative: SidebandRefinementCheck,
}

impl SidebandCarrierRefinementCheck {
    /// Both full sampled equations, grid checks and finite omitted bands pass.
    /// This alone does not provide actual inter-order convergence evidence.
    pub fn equations_resolved(&self) -> bool {
        [&self.response, &self.angular_derivative].iter().all(|c| {
            c.equation_residual_ratio <= 1.0
                && c.sampling_defect_ratio <= 1.0
                && c.omitted_residual_ratio <= 1.0
        })
    }
    /// All equation checks AND both actually measured state changes pass.
    pub fn is_resolved(&self) -> bool {
        self.equations_resolved()
            && self.response.state_change_ratio.is_some_and(|x| x <= 1.0)
            && self
                .angular_derivative
                .state_change_ratio
                .is_some_and(|x| x <= 1.0)
    }
}

/// An assessed candidate, which may still require more sidebands.
///
/// Inspect `check`; this type deliberately does not imply spectral acceptance.
/// The original retained systems are always solved and backward-error checked.
#[derive(Clone, Debug)]
pub struct AssessedSidebandCarrierDerivativeResponse {
    pub response: SidebandCarrierDerivativeResponse,
    pub harmonics: usize,
    pub carrier_angular_frequency: f64,
    pub orbital_angular_frequency: f64,
    pub base_residual: f64,
    pub tangent_sampling_defect: f64,
    pub check: SidebandCarrierRefinementCheck,
}

impl AssessedSidebandCarrierDerivativeResponse {
    /// Compare two ACTUALLY solved orders at the same carrier and base clock.
    /// The caller must use the same original source, base, numerical chart and
    /// forcing family. Missing lower-order bands are zero only in this norm;
    /// the current response itself is never obtained by padding or clipping.
    /// Both original-state and derivative units have independent tolerances.
    ///
    /// # Errors
    /// Invalid tolerances, incompatible shapes/frequencies or nonincreasing order.
    pub fn compare_with(
        &mut self,
        previous: &Self,
        tolerances: &SidebandCarrierRefinementTolerances,
    ) -> Result<()> {
        let d = tolerances.state_tolerances.len();
        tolerances.admit(d)?;
        if d == 0
            || self.harmonics > 64
            || previous.harmonics >= self.harmonics
            || self.carrier_angular_frequency.to_bits()
                != previous.carrier_angular_frequency.to_bits()
            || self.orbital_angular_frequency.to_bits()
                != previous.orbital_angular_frequency.to_bits()
            || self.response.response.coefficients.len() != d * (2 * self.harmonics + 1)
            || self.response.angular_derivative.coefficients.len() != d * (2 * self.harmonics + 1)
            || previous.response.response.coefficients.len() != d * (2 * previous.harmonics + 1)
            || previous.response.angular_derivative.coefficients.len()
                != d * (2 * previous.harmonics + 1)
        {
            return Err(invalid(
                "actual increasing orders at the same carrier, clock and original state dimension required",
            ));
        }
        self.check.response.state_change_ratio = Some(state_change(
            &previous.response.response.coefficients,
            previous.harmonics,
            &self.response.response.coefficients,
            self.harmonics,
            &tolerances.state_tolerances,
        ));
        self.check.angular_derivative.state_change_ratio = Some(state_change(
            &previous.response.angular_derivative.coefficients,
            previous.harmonics,
            &self.response.angular_derivative.coefficients,
            self.harmonics,
            &tolerances.derivative_state_tolerances,
        ));
        Ok(())
    }
}

/// A jointly resolved finite response/derivative and its actual refinement chain.
#[derive(Clone, Debug)]
pub struct AdaptiveSidebandCarrierDerivativeResponse {
    pub response: SidebandCarrierDerivativeResponse,
    pub harmonics: usize,
    pub base_residual: f64,
    pub tangent_sampling_defect: f64,
    pub refinement_trace: Vec<SidebandCarrierRefinementCheck>,
}

impl<P: OrbitProblem> PreparedPeriodicLinearization<'_, P> {
    /// Actually solve and independently assess a response/derivative at an order
    /// no larger than this preparation, reusing its SAME cached original tangent.
    ///
    /// One LU serves X and X_nu. Its exact original port derivative is sampled
    /// once per retained signed band and retained in original physical units for
    /// the independent equation check. After the factor is released, BOTH full
    /// original equations are sampled together on shifted 2*M+1 and 4*M+3 grids:
    /// `P X - A(t) X - F` and `P X_nu - A(t) X_nu + P_nu X - F_nu`.
    /// One original tangent callback per sample serves both fields. The original
    /// prepared physical base is frozen in a bounded snapshot, without diagnostics.
    ///
    /// Only the finite omitted-band ACCEPTANCE gates are deferred. Their
    /// measurements, full original equations and separate physical units remain
    /// in `check`. Singularities, source refusals and retained backward-error
    /// failures still refuse. A candidate is not an accepted adaptive response.
    /// `compare_with` supplies actual inter-order changes before `is_resolved`.
    ///
    /// Lower and upper orders use the preparation's common fine quadrature. This
    /// allows a resonance campaign to compare independently solved low/high
    /// systems at EVERY queried carrier while keeping one tangent and one live
    /// LU. It measures truncation convergence, not a second base/quadrature proof.
    ///
    /// Additional derivative limits now also admit this pair's full assessments,
    /// padded inputs, one retained comparator pair and bounded refinement trace.
    /// Core and additional payload/work allowances add. Callback internals remain
    /// outside the numerical work model, as in the fixed-order owner.
    ///
    /// # Errors
    /// Invalid order/data/limits, additional budget exhaustion before new source
    /// calls, any original retained numerical/source refusal, or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_carrier_derivative_assessed_at_order(
        &self,
        harmonics: usize,
        carrier: f64,
        forcing: &[C64],
        forcing_angular_derivative: &[C64],
        derivative_controls: SidebandCarrierDerivativeControls,
        tolerances: &SidebandCarrierRefinementTolerances,
        derivative: impl FnMut(f64, &mut [C64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AssessedSidebandCarrierDerivativeResponse> {
        poll(&mut cancelled)?;
        if harmonics > self.harmonics() {
            return Err(invalid("assessed order exceeds the prepared tangent"));
        }
        tolerances.admit(self.dimension)?;
        admit_data(
            self.dimension,
            harmonics,
            carrier,
            forcing,
            forcing_angular_derivative,
        )?;
        admit_derivative(derivative_controls)?;
        let controls = SidebandControls {
            harmonics,
            ..self.controls
        };
        let (payload, work) =
            additional_costs(self.dimension, self.orbit.coeffs[0].len(), controls);
        bound(
            "carrier refinement payload bytes",
            payload,
            derivative_controls.maximum_payload_bytes as u128,
        )?;
        bound(
            "carrier refinement work",
            work,
            derivative_controls.maximum_work as u128,
        )?;
        let callback = self
            .problem
            .island_state_jacobian()
            .ok_or(SidebandError::TangentUnavailable)?;
        let (response, original_derivative) = {
            let factor = self.factor_at_order(harmonics, carrier, &mut cancelled)?;
            let sensitivity = factor.prepare_carrier_derivative(
                derivative_controls,
                derivative,
                &mut cancelled,
            )?;
            let response = sensitivity.solve_candidate(
                forcing,
                forcing_angular_derivative,
                false,
                &mut cancelled,
            )?;
            (response, sensitivity.original_derivative)
        }; // No dense factor is retained during full original-equation sampling.
        let measured = assess_equations(
            self.problem,
            callback,
            &self.orbit,
            carrier,
            controls,
            &self.scales,
            &[
                EquationAssessment {
                    forcing,
                    response: &response.response,
                    tolerances: &tolerances.equation_tolerances,
                    additional_port: None,
                },
                EquationAssessment {
                    forcing: forcing_angular_derivative,
                    response: &response.angular_derivative,
                    tolerances: &tolerances.derivative_equation_tolerances,
                    additional_port: Some((&original_derivative, &response.response.coefficients)),
                },
            ],
            &mut cancelled,
        )?;
        let check = |index: usize, omitted, allowance| SidebandRefinementCheck {
            harmonics,
            assessment_samples: 4 * controls.samples_per_period + 3,
            equation_residual_ratio: measured[index].0,
            sampling_defect_ratio: measured[index].1,
            state_change_ratio: None,
            omitted_residual_ratio: omitted / allowance,
        };
        let check = SidebandCarrierRefinementCheck {
            response: check(
                0,
                response.response.omitted_residual,
                controls.maximum_omitted_residual,
            ),
            angular_derivative: check(
                1,
                response.angular_derivative.omitted_residual,
                derivative_controls.maximum_omitted_residual,
            ),
        };
        poll(&mut cancelled)?;
        Ok(AssessedSidebandCarrierDerivativeResponse {
            response,
            harmonics,
            carrier_angular_frequency: carrier,
            orbital_angular_frequency: self.omega,
            base_residual: self.base_residual,
            tangent_sampling_defect: self.tangent_sampling_defect,
            check,
        })
    }
}

/// Resolve response and exact carrier sensitivity together at a fixed original
/// periodic operating orbit. Orders grow as min(2*N+1, maximum_harmonics).
/// Both supplied forcing arrays declare their COMPLETE initial signed support;
/// they are zero-padded exactly, never folded, conjugated or clipped.
///
/// Success requires independently sampled full equations, finite omitted bands,
/// grid agreement, and actual inter-order coefficient changes for BOTH fields.
/// A smooth-looking response cannot hide an unresolved derivative. The base,
/// clock and phase are fixed: this is not differentiation along an orbit branch.
///
/// Core `maximum_work` bounds the cumulative core costs of all attempted rounds;
/// derivative `maximum_work` separately bounds the cumulative added derivative
/// solves and joint assessments. Unattempted orders are not charged. Their
/// payload caps add and include the live comparator and bounded history. Finite
/// sampled convergence is measured evidence, not an infinite-tail, forward-error
/// or Floquet stability enclosure; original base resolution remains independent.
///
/// # Errors
/// Original source/numerical refusals, cumulative budget exhaustion, cancellation
/// or `CarrierHarmonicsUnresolved` with the last measurements and no partial pair.
#[allow(clippy::too_many_arguments)]
pub fn solve_sideband_carrier_derivative_adaptive<P: OrbitProblem>(
    problem: &P,
    orbit: &HbOrbit,
    carrier: f64,
    forcing: &[C64],
    forcing_angular_derivative: &[C64],
    controls: SidebandControls,
    scales: &SidebandScales,
    derivative_controls: SidebandCarrierDerivativeControls,
    refinement: &SidebandCarrierRefinementControls,
    mut derivative: impl FnMut(f64, &mut [C64]) -> bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveSidebandCarrierDerivativeResponse> {
    poll(&mut cancelled)?;
    let d = problem.dim();
    admit(d, orbit, controls, scales)?;
    admit_data(
        d,
        controls.harmonics,
        carrier,
        forcing,
        forcing_angular_derivative,
    )?;
    admit_derivative(derivative_controls)?;
    refinement.tolerances.admit(d)?;
    if refinement.maximum_harmonics <= controls.harmonics || refinement.maximum_harmonics > 64 {
        return Err(invalid(
            "joint refinement requires distinct orders bounded by 64",
        ));
    }
    let mut previous: Option<AssessedSidebandCarrierDerivativeResponse> = None;
    let mut trace = Vec::new();
    let mut harmonics = controls.harmonics;
    let mut core_spent = 0_u128;
    let mut derivative_spent = 0_u128;
    loop {
        poll(&mut cancelled)?;
        let current = SidebandControls {
            harmonics,
            samples_per_period: controls
                .samples_per_period
                .max(2 * (3 * harmonics + 2).max(orbit.coeffs[0].len() - 1) + 1),
            ..controls
        };
        admit(d, orbit, current, scales)?;
        let (_, work) = resource_costs(d, orbit.coeffs[0].len(), current);
        let (extra_payload, extra_work) = additional_costs(d, orbit.coeffs[0].len(), current);
        core_spent += work;
        derivative_spent += extra_work;
        // Admit the ENTIRE next round before allocation, tangent or port calls.
        bound(
            "adaptive carrier core work",
            core_spent,
            controls.maximum_work as u128,
        )?;
        bound(
            "adaptive carrier derivative work",
            derivative_spent,
            derivative_controls.maximum_work as u128,
        )?;
        bound(
            "carrier refinement payload bytes",
            extra_payload,
            derivative_controls.maximum_payload_bytes as u128,
        )?;
        let n = d * (2 * harmonics + 1);
        let offset = (harmonics - controls.harmonics) * d;
        let mut padded = vec![C64::ZERO; n];
        let mut padded_derivative = vec![C64::ZERO; n];
        padded[offset..offset + forcing.len()].copy_from_slice(forcing);
        padded_derivative[offset..offset + forcing.len()]
            .copy_from_slice(forcing_angular_derivative);
        let mut candidate = {
            let prepared = PreparedPeriodicLinearization::new(
                problem,
                orbit,
                current,
                scales,
                &mut cancelled,
            )?;
            prepared.solve_carrier_derivative_assessed_at_order(
                harmonics,
                carrier,
                &padded,
                &padded_derivative,
                derivative_controls,
                &refinement.tolerances,
                &mut derivative,
                &mut cancelled,
            )?
        };
        if let Some(previous) = &previous {
            candidate.compare_with(previous, &refinement.tolerances)?;
        }
        trace.push(candidate.check);
        poll(&mut cancelled)?;
        if candidate.check.is_resolved() {
            return Ok(AdaptiveSidebandCarrierDerivativeResponse {
                response: candidate.response,
                harmonics,
                base_residual: candidate.base_residual,
                tangent_sampling_defect: candidate.tangent_sampling_defect,
                refinement_trace: trace,
            });
        }
        if harmonics == refinement.maximum_harmonics {
            return Err(SidebandError::CarrierHarmonicsUnresolved {
                last: candidate.check,
            });
        }
        previous = Some(candidate);
        harmonics = (2 * harmonics + 1).min(refinement.maximum_harmonics);
    }
}

fn admit_data(d: usize, h: usize, carrier: f64, forcing: &[C64], slope: &[C64]) -> Result<()> {
    let n = d * (2 * h + 1);
    if !carrier.is_finite()
        || forcing.len() != n
        || slope.len() != n
        || forcing.iter().chain(slope).any(|&z| !finite(z))
    {
        return Err(invalid(
            "finite carrier and complete original forcing/derivative arrays required",
        ));
    }
    Ok(())
}

fn admit_derivative(controls: SidebandCarrierDerivativeControls) -> Result<()> {
    if !positive(controls.maximum_relative_residual)
        || controls.maximum_relative_residual >= 1.0
        || !positive(controls.maximum_omitted_residual)
    {
        return Err(invalid(
            "positive derivative limits and backward error below one required",
        ));
    }
    Ok(())
}

fn additional_costs(
    dimension: usize,
    base_width: usize,
    controls: SidebandControls,
) -> (u128, u128) {
    let d = dimension as u128;
    let n = d * (2 * controls.harmonics + 1) as u128;
    let base = d * base_width as u128;
    let samples = (6 * controls.samples_per_period + 4) as u128;
    let (payload, work) = super::resource_costs(dimension, controls.harmonics);
    // Original port derivatives, two normalized full-equation envelopes, raw
    // tangent scratch, both padded sources, a prior pair and <=8 trace entries.
    (
        payload + 192 * n + 256 * d * d + 64 * base + 256 * d + 4096,
        work + samples * 512 * (d * d + 4 * n + base) + 512 * n * d + 256 * n + 2048,
    )
}
