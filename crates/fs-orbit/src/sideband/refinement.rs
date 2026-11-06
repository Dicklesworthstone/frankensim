//! Measured sideband-order refinement against the original perturbation equation.

use super::{
    ComplexSum, HbOrbit, PreparedPeriodicLinearization, Result, SidebandControls, SidebandError,
    SidebandResponse, SidebandScales, TAU, admit, bound, check_port, finite, invalid, poll,
    positive, resource_costs, scale_entry,
};
use crate::{AftSum, IslandStateJacobian, OrbitProblem};
use fs_math::{c64::C64, dd::Dd, det};

/// Physical resolution requirements in addition to the fixed-order controls.
#[derive(Clone, Debug)]
pub struct SidebandRefinementControls {
    /// Largest solved sideband order, strictly above the starting order and at
    /// most 64. Acceptance requires two distinct orders actually to be solved.
    pub maximum_harmonics: usize,
    /// Absolute complex-envelope L2 tolerance per ORIGINAL equation row.
    /// These physical units are independent of numerical conditioning scales.
    pub equation_tolerances: Vec<f64>,
    /// Absolute complex-envelope L2 change per ORIGINAL state coordinate,
    /// including all newly added bands, between successive solved orders.
    pub state_tolerances: Vec<f64>,
}

/// Independent measurements at one actually solved retained order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SidebandRefinementCheck {
    /// Retained signed indices are -harmonics through +harmonics.
    pub harmonics: usize,
    /// Samples on the finer of two independent shifted assessment grids.
    pub assessment_samples: usize,
    /// Worst per-row full sampled equation L2 / physical tolerance, taking
    /// the larger measurement on the two grids, including distant harmonics.
    pub equation_residual_ratio: f64,
    /// Worst per-row difference between the two L2 measurements / tolerance.
    pub sampling_defect_ratio: f64,
    /// Worst per-state response change / physical tolerance. Absent at the
    /// first solve; zero-padding the same response never supplies this evidence.
    pub state_change_ratio: Option<f64>,
    /// Existing finite omitted-band defect / its fixed-order allowance.
    pub omitted_residual_ratio: f64,
}

/// A resolved finite response with measured, not certified, convergence evidence.
#[derive(Clone, Debug)]
pub struct AdaptiveSidebandResponse {
    /// Complete original-state perturbation, at the accepted retained order.
    pub response: SidebandResponse,
    /// Sideband order of `response`; its indexing uses this order, not the seed.
    pub harmonics: usize,
    /// Independently recomputed original base-equation residual at acceptance.
    pub base_residual: f64,
    /// Independent tangent Fourier sampling difference at acceptance.
    pub tangent_sampling_defect: f64,
    /// Every actually solved and independently assessed retained order.
    pub refinement_trace: Vec<SidebandRefinementCheck>,
}

/// Resolve frequency conversion without mistaking low-order agreement for
/// convergence. Retained orders grow as `min(2*N+1, maximum_harmonics)`.
///
/// `forcing` declares the COMPLETE input at the initial order, in the same
/// signed peak-amplitude convention as [`super::PreparedSidebandResponse::solve`].
/// Its support is preserved exactly by zero-padding, never clipped or folded.
/// Each candidate uses the existing complete port/tangent operator and native
/// complex factorization. An intermediate omitted-band defect requests more
/// bands; singularities and retained backward-error failures still refuse.
///
/// Independently of Fourier guard bands, the original equation envelope
/// `sum_k [P(i*(nu+k*omega))*X_k-F_k]*exp(i*k*theta)
///  - A(theta)*sum_k X_k*exp(i*k*theta)` is sampled on two shifted, unequal
/// grids. This detects distant generated bands even when two truncated solves
/// agree. Acceptance requires both full equation measurements, finite omitted
/// bands, and a comparison of two actually solved state envelopes to pass.
///
/// L2 here means the COMPLEX envelope norm: signed coefficients have weight
/// one. It is not the real physical RMS when signed physical frequencies
/// coincide; consumers still combine those contributions coherently.
///
/// The initial fixed-order controls must be valid. Their `samples_per_period`
/// is then a minimum and grows to resolve each subsequent retained order.
/// `maximum_work` is a CUMULATIVE allowance for all attempted preparations,
/// factorizations, solves and assessments in this call. `maximum_payload_bytes`
/// covers all simultaneously live adaptive numerical buffers, including the
/// previous response. Both use the fixed owner's resource model and exclude
/// source/caller storage, callback internals, allocator metadata and stacks.
/// Unused higher orders are not charged. No samples-by-state history is stored.
///
/// Finite-grid agreement is not an infinite-tail, forward-error or stability
/// certificate. The original base and tangent must repeat at the stated period.
/// Arbitrary frequencies beyond BOTH assessment grids can remain unresolved.
/// Tighten the sampling and base-orbit resolution independently when needed.
///
/// # Errors
/// The original fixed-order refusals, cumulative resource exhaustion, or
/// [`SidebandError::HarmonicsUnresolved`] with the last measurements. No partial
/// accepted response is published on any refusal or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn solve_sideband_adaptive<P: OrbitProblem>(
    problem: &P,
    orbit: &HbOrbit,
    carrier_angular_frequency: f64,
    forcing: &[C64],
    controls: SidebandControls,
    scales: &SidebandScales,
    refinement: &SidebandRefinementControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveSidebandResponse> {
    poll(&mut cancelled)?;
    let d = problem.dim();
    admit(d, orbit, controls, scales)?;
    if !carrier_angular_frequency.is_finite()
        || refinement.maximum_harmonics <= controls.harmonics
        || refinement.maximum_harmonics > 64
        || refinement.equation_tolerances.len() != d
        || refinement.state_tolerances.len() != d
        || refinement
            .equation_tolerances
            .iter()
            .chain(&refinement.state_tolerances)
            .any(|&v| !positive(v))
        || forcing.len() != d * (2 * controls.harmonics + 1)
        || forcing.iter().any(|&z| !finite(z))
    {
        return Err(invalid(
            "adaptive sidebands need distinct bounded orders, complete finite forcing, a finite carrier and positive physical tolerances",
        ));
    }
    let callback = problem
        .island_state_jacobian()
        .ok_or(SidebandError::TangentUnavailable)?;
    let mut previous: Option<(usize, SidebandResponse)> = None;
    let mut trace = Vec::new();
    let mut harmonics = controls.harmonics;
    let mut spent = 0_u128;
    loop {
        poll(&mut cancelled)?;
        let mut current = controls;
        current.harmonics = harmonics;
        current.samples_per_period = controls
            .samples_per_period
            .max(2 * (3 * harmonics + 2).max(orbit.coeffs[0].len() - 1) + 1);
        admit(d, orbit, current, scales)?;
        let n = (d * (2 * harmonics + 1)) as u128;
        let dimension = d as u128;
        let base = dimension * orbit.coeffs[0].len() as u128;
        let (core_payload, core_work) = resource_costs(d, orbit.coeffs[0].len(), current);
        // Conservative simultaneous upper bound: previous/current responses,
        // padded forcing, tolerance-scaled assessment coefficients, streamed
        // physical waveform/Jacobian scratch and at most eight history entries.
        let payload = core_payload
            + 128 * n
            + 256 * dimension * dimension
            + 64 * base
            + 128 * dimension
            + 2048;
        bound(
            "adaptive numerical payload bytes",
            payload,
            controls.maximum_payload_bytes as u128,
        )?;
        let samples = (6 * current.samples_per_period + 4) as u128;
        let assessment_work = samples * 256 * (dimension * dimension + 2 * n + base)
            + 256 * n * dimension
            + 128 * n
            + 1024;
        spent += core_work + assessment_work;
        bound(
            "adaptive numerical work",
            spent,
            controls.maximum_work as u128,
        )?;

        let mut padded = vec![C64::ZERO; n as usize];
        let offset = (harmonics - controls.harmonics) * d;
        padded[offset..offset + forcing.len()].copy_from_slice(forcing);
        let (response, base_residual, tangent_sampling_defect) = {
            let prepared = PreparedPeriodicLinearization::new(
                problem,
                orbit,
                current,
                scales,
                &mut cancelled,
            )?;
            let factor = prepared.factor(carrier_angular_frequency, &mut cancelled)?;
            let response = factor.solve_candidate(&padded, false, &mut cancelled)?;
            (
                response,
                prepared.base_residual,
                prepared.tangent_sampling_defect,
            )
        }; // Release the dense factor before the independent assessment.
        let (equation_residual_ratio, sampling_defect_ratio) = assess(
            problem,
            callback,
            orbit,
            carrier_angular_frequency,
            &padded,
            &response,
            current,
            scales,
            &refinement.equation_tolerances,
            &mut cancelled,
        )?;
        let state_change_ratio = previous.as_ref().map(|(old_h, old)| {
            state_change(
                &old.coefficients,
                *old_h,
                &response.coefficients,
                harmonics,
                &refinement.state_tolerances,
            )
        });
        let check = SidebandRefinementCheck {
            harmonics,
            assessment_samples: 4 * current.samples_per_period + 3,
            equation_residual_ratio,
            sampling_defect_ratio,
            state_change_ratio,
            omitted_residual_ratio: response.omitted_residual / controls.maximum_omitted_residual,
        };
        trace.push(check);
        poll(&mut cancelled)?;
        if check.equation_residual_ratio <= 1.0
            && check.sampling_defect_ratio <= 1.0
            && check.omitted_residual_ratio <= 1.0
            && check.state_change_ratio.is_some_and(|v| v <= 1.0)
        {
            return Ok(AdaptiveSidebandResponse {
                response,
                harmonics,
                base_residual,
                tangent_sampling_defect,
                refinement_trace: trace,
            });
        }
        if harmonics == refinement.maximum_harmonics {
            return Err(SidebandError::HarmonicsUnresolved { last: check });
        }
        previous = Some((harmonics, response));
        harmonics = (2 * harmonics + 1).min(refinement.maximum_harmonics);
    }
}

pub(super) fn state_change(
    old: &[C64],
    old_h: usize,
    new: &[C64],
    new_h: usize,
    tolerances: &[f64],
) -> f64 {
    let d = tolerances.len();
    let offset = (new_h - old_h) * d;
    let mut worst = 0.0_f64;
    for (row, &tolerance) in tolerances.iter().enumerate() {
        let mut norm = 0.0;
        for band in 0..2 * new_h + 1 {
            let index = band * d + row;
            let prior = index
                .checked_sub(offset)
                .and_then(|i| old.get(i))
                .copied()
                .unwrap_or(C64::ZERO);
            let re = (Dd::from_f64(new[index].re) - Dd::from_f64(prior.re)).to_f64() / tolerance;
            let im = (Dd::from_f64(new[index].im) - Dd::from_f64(prior.im)).to_f64() / tolerance;
            norm = det::hypot(det::hypot(norm, re), im);
        }
        if !norm.is_finite() {
            return f64::INFINITY;
        }
        worst = worst.max(norm);
    }
    worst
}

/// One complete original equation envelope, with an optional linear chain rule.
pub(super) struct EquationAssessment<'a> {
    pub forcing: &'a [C64],
    pub response: &'a SidebandResponse,
    pub tolerances: &'a [f64],
    // Raw original dP/dnu blocks and the primal coefficients. The derivative
    // equation includes +P_nu X, independently of the rounded solve RHS.
    pub additional_port: Option<(&'a [C64], &'a [C64])>,
}

struct EquationEnvelope {
    normalized: Vec<C64>,
    linear: Vec<C64>,
    chart: SidebandScales,
}

#[allow(clippy::too_many_arguments)]
fn assess<P: OrbitProblem>(
    problem: &P,
    callback: IslandStateJacobian<P>,
    orbit: &HbOrbit,
    carrier: f64,
    forcing: &[C64],
    response: &SidebandResponse,
    controls: SidebandControls,
    scales: &SidebandScales,
    tolerances: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64)> {
    Ok(assess_equations(
        problem,
        callback,
        orbit,
        carrier,
        controls,
        scales,
        &[EquationAssessment {
            forcing,
            response,
            tolerances,
            additional_port: None,
        }],
        cancelled,
    )?[0])
}

/// Stream several envelopes through the SAME original tangent samples. In
/// particular the response and its derivative do not double physical Hessian
/// callbacks or retain a samples-by-state/matrix history.
#[allow(clippy::too_many_arguments)]
pub(super) fn assess_equations<P: OrbitProblem>(
    problem: &P,
    callback: IslandStateJacobian<P>,
    orbit: &HbOrbit,
    carrier: f64,
    controls: SidebandControls,
    scales: &SidebandScales,
    fields: &[EquationAssessment<'_>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<(f64, f64)>> {
    let d = problem.dim();
    let h = controls.harmonics;
    let n = d * (2 * h + 1);
    let mut envelopes = Vec::with_capacity(fields.len());
    let mut additional = Vec::with_capacity(fields.len());
    for field in fields {
        // Divide by physical tolerances BEFORE norms: squaring tiny original
        // equation units must not make a nonzero defect disappear.
        let chart = SidebandScales {
            states: scales.states.clone(),
            equations: field.tolerances.to_vec(),
        };
        let normalized = field
            .response
            .coefficients
            .iter()
            .enumerate()
            .map(|(i, &z)| divide(z, scales.states[i % d]))
            .collect::<Result<Vec<_>>>()?;
        additional.push(match field.additional_port {
            Some((_, primal)) => Some(
                primal
                    .iter()
                    .enumerate()
                    .map(|(i, &z)| divide(z, scales.states[i % d]))
                    .collect::<Result<Vec<_>>>()?,
            ),
            None => None,
        });
        envelopes.push(EquationEnvelope {
            normalized,
            linear: vec![C64::ZERO; n],
            chart,
        });
    }
    for band in 0..2 * h + 1 {
        poll(cancelled)?;
        let frequency = (band as f64 - h as f64).mul_add(orbit.omega, carrier);
        if !frequency.is_finite() {
            return Err(invalid("shifted assessment frequency is not representable"));
        }
        let port = problem.port(C64::new(0.0, frequency));
        poll(cancelled)?;
        check_port(&port, d)?;
        for (index, (field, envelope)) in fields.iter().zip(&mut envelopes).enumerate() {
            for row in 0..d {
                poll(cancelled)?;
                let mut sum = ComplexSum::new();
                for col in 0..d {
                    sum.product(
                        scale_entry(port[row * d + col], row, col, &envelope.chart)?,
                        envelope.normalized[band * d + col],
                    );
                    if let Some((blocks, _)) = field.additional_port {
                        sum.product(
                            scale_entry(
                                blocks[(band * d + row) * d + col],
                                row,
                                col,
                                &envelope.chart,
                            )?,
                            additional[index].as_ref().unwrap()[band * d + col],
                        );
                    }
                }
                let input = divide(field.forcing[band * d + row], field.tolerances[row])?;
                sum.re = sum.re - Dd::from_f64(input.re);
                sum.im = sum.im - Dd::from_f64(input.im);
                let value = sum.value();
                if !finite(value) {
                    return Err(invalid("linear perturbation envelope is not representable"));
                }
                envelope.linear[band * d + row] = value;
            }
        }
    }
    drop(additional);
    let coarse = sample_equations(
        problem,
        callback,
        orbit,
        h,
        &envelopes,
        2 * controls.samples_per_period + 1,
        0.25,
        cancelled,
    )?;
    let fine = sample_equations(
        problem,
        callback,
        orbit,
        h,
        &envelopes,
        4 * controls.samples_per_period + 3,
        0.5,
        cancelled,
    )?;
    let mut measurements = Vec::with_capacity(fields.len());
    for (coarse_rows, fine_rows) in coarse.chunks_exact(d).zip(fine.chunks_exact(d)) {
        let mut total = 0.0_f64;
        let mut sampling = 0.0_f64;
        for (&a, &b) in coarse_rows.iter().zip(fine_rows) {
            total = total.max(a).max(b);
            sampling = sampling.max((a - b).abs());
        }
        measurements.push((total, sampling));
    }
    Ok(measurements)
}

fn divide(value: C64, scale: f64) -> Result<C64> {
    let result = C64::new(value.re / scale, value.im / scale);
    if !finite(result)
        || (value.re != 0.0 && result.re == 0.0)
        || (value.im != 0.0 && result.im == 0.0)
    {
        return Err(invalid(
            "normalized perturbation coefficient is not representable",
        ));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn sample_equations<P: OrbitProblem>(
    problem: &P,
    callback: IslandStateJacobian<P>,
    orbit: &HbOrbit,
    h: usize,
    fields: &[EquationEnvelope],
    samples: usize,
    shift: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    let d = problem.dim();
    let base_order = orbit.coeffs[0].len() - 1;
    let mut phases = vec![C64::ZERO; h.max(base_order) + 1];
    let mut base = vec![0.0; d];
    let mut envelope = vec![C64::ZERO; d];
    let mut jacobian = vec![0.0; d * d];
    let mut norms = vec![0.0; d * fields.len()];
    let weight = 1.0 / det::sqrt(samples as f64);
    let period = TAU / orbit.omega;
    for sample in 0..samples {
        poll(cancelled)?;
        let fraction = (sample as f64 + shift) / samples as f64;
        let phase = TAU * fraction;
        for (order, value) in phases.iter_mut().enumerate() {
            let (sine, cosine) = det::sin_cos(order as f64 * phase);
            *value = C64::new(cosine, sine);
        }
        for row in 0..d {
            let mut sum = AftSum::new(orbit.coeffs[row][0].re);
            for (order, &value) in orbit.coeffs[row].iter().enumerate().skip(1) {
                sum.add(2.0 * (value.re * phases[order].re - value.im * phases[order].im));
            }
            base[row] = sum.total();
        }
        if base.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "perturbation assessment waveform is not representable",
            ));
        }
        jacobian.fill(f64::NAN);
        if !callback(problem, period * fraction, &base, &mut jacobian)
            || jacobian.iter().any(|v| !v.is_finite())
        {
            return Err(SidebandError::TangentRefused { sample, samples });
        }
        poll(cancelled)?;
        for (field_index, field) in fields.iter().enumerate() {
            for row in 0..d {
                let mut state = ComplexSum::new();
                for band in 0..2 * h + 1 {
                    let k = band as isize - h as isize;
                    let p = if k < 0 {
                        phases[k.unsigned_abs()].conj()
                    } else {
                        phases[k as usize]
                    };
                    state.product(field.normalized[band * d + row], p);
                }
                envelope[row] = state.value();
            }
            if envelope.iter().any(|&z| !finite(z)) {
                return Err(invalid(
                    "perturbation assessment waveform is not representable",
                ));
            }
            for row in 0..d {
                poll(cancelled)?;
                let mut balance = ComplexSum::new();
                for band in 0..2 * h + 1 {
                    let k = band as isize - h as isize;
                    let p = if k < 0 {
                        phases[k.unsigned_abs()].conj()
                    } else {
                        phases[k as usize]
                    };
                    balance.product(field.linear[band * d + row], p);
                }
                for col in 0..d {
                    let entry = scale_entry(
                        C64::from_re(jacobian[row * d + col]),
                        row,
                        col,
                        &field.chart,
                    )?;
                    balance.product(-entry, envelope[col]);
                }
                let defect = balance.value();
                if !finite(defect) {
                    return Err(invalid("full perturbation equation is not representable"));
                }
                let norm = &mut norms[field_index * d + row];
                *norm = det::hypot(det::hypot(*norm, weight * defect.re), weight * defect.im);
                if !norm.is_finite() {
                    return Err(invalid("full perturbation norm is not representable"));
                }
            }
        }
    }
    Ok(norms)
}
