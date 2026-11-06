//! Exact carrier sensitivities at a fixed complete periodic operating orbit.

use super::{
    AftSum, ComplexSum, PreparedSidebandResponse, Result, SidebandError, SidebandResponse, bound,
    finite, invalid, poll, positive, scale_entry,
};
use fs_math::{c64::C64, dd::Dd};

pub mod refinement;

/// Additional limits for a prepared carrier derivative and each sensitivity solve.
#[derive(Clone, Copy, Debug)]
pub struct SidebandCarrierDerivativeControls {
    /// Additional live numerical payload beyond the original factor/solve cap.
    /// Includes cached full port derivatives, the retained primal response,
    /// normalized derivative forcing and original-equation assessment scratch.
    /// Source/caller storage, callback internals and allocator metadata are excluded.
    pub maximum_payload_bytes: usize,
    /// Additional numerical work for derivative preparation plus one differentiated
    /// solve and its assessments. The primal operation keeps its original cap;
    /// these two allowances add. Each later query has the same allowance.
    pub maximum_work: u64,
    /// Componentwise backward-error limit for the differentiated retained system,
    /// including an independent check of D X' + P' X = F' in original units.
    pub maximum_relative_residual: f64,
    /// Omitted derivative-equation norm in equation-scale units PER angular
    /// frequency, for N < |k| <= 2*N+2. Independent of the primal allowance.
    pub maximum_omitted_residual: f64,
}

/// Response and its derivative with respect to the signed carrier frequency.
#[derive(Clone, Debug)]
pub struct SidebandCarrierDerivativeResponse {
    /// Original peak complex response coefficients and their own error checks.
    pub response: SidebandResponse,
    /// Derivatives of those coefficients, in state units per (radian/unit time).
    /// This is NOT the time derivative of the perturbation waveform. Base orbit,
    /// base frequency, phase gauge and sideband indices are held fixed.
    /// The retained residual includes the complete differentiated original equation.
    pub angular_derivative: SidebandResponse,
}

/// Cached exact port derivatives borrowing one unchanged original operator/LU.
///
/// The complete physical derivative `dP(i*w)/dw` is sampled at every signed
/// shifted frequency, normalized in the SAME numerical chart as the factor,
/// and reused across right-hand sides. This works with full mass, damping,
/// gyroscopic and dispersive/memory ports without inferring a time realization.
/// The original source and the derivative declaration must remain unchanged.
pub struct PreparedSidebandCarrierDerivative<'factor, 'prepared> {
    factor: &'factor PreparedSidebandResponse<'prepared>,
    derivative: Vec<C64>,
    // Original entries support a separate physical-tolerance assessment without
    // rescaling the rounded numerical chart or repeating derivative callbacks.
    original_derivative: Vec<C64>,
    controls: SidebandCarrierDerivativeControls,
}

impl<'prepared> PreparedSidebandResponse<'prepared> {
    /// Prepare an EXACT derivative of this factor's original frequency port.
    ///
    /// `derivative(w, out)` must write every row-major entry of `dP(i*w)/dw`
    /// in original equation/state units and return true. This is the derivative
    /// along the imaginary axis, including the factor i; it is NOT `dP(s)/ds`.
    /// Signed shifted frequencies are passed unchanged. Native/authored exact
    /// derivative actions may be used; no finite differences or model fitting
    /// are performed. NaN prefills detect incomplete declarations.
    ///
    /// The extra payload/work is admitted before any derivative callback or
    /// allocation. No new factorization is needed. A complete derivative matrix
    /// is cached once per sideband, independently of the number of excitations.
    ///
    /// # Errors
    /// Invalid limits, exhausted extra budgets, refused/incomplete/nonfinite
    /// declared derivatives, unrepresentable scaling, or cancellation.
    pub fn prepare_carrier_derivative<'factor>(
        &'factor self,
        controls: SidebandCarrierDerivativeControls,
        mut derivative: impl FnMut(f64, &mut [C64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedSidebandCarrierDerivative<'factor, 'prepared>> {
        poll(&mut cancelled)?;
        if !positive(controls.maximum_relative_residual)
            || controls.maximum_relative_residual >= 1.0
            || !positive(controls.maximum_omitted_residual)
        {
            return Err(invalid(
                "positive carrier derivative limits and backward error below one required",
            ));
        }
        let d = self.dimension;
        let h = self.controls.harmonics;
        let entries = d * d * (2 * h + 1);
        let (payload, work) = resource_costs(d, h);
        bound(
            "carrier derivative payload bytes",
            payload,
            controls.maximum_payload_bytes as u128,
        )?;
        // The native correction owner permits up to five residual evaluations
        // and five triangular solves. Include both original derivative products
        // and the finite omitted-band contraction with double-double sums.
        bound(
            "carrier derivative work",
            work,
            controls.maximum_work as u128,
        )?;
        let mut blocks = vec![C64::ZERO; entries];
        let mut original_derivative = vec![C64::ZERO; entries];
        let mut port = vec![C64::ZERO; d * d];
        for band in 0..2 * h + 1 {
            poll(&mut cancelled)?;
            let k = band as isize - h as isize;
            let frequency = (k as f64).mul_add(self.omega, self.carrier);
            port.fill(C64::new(f64::NAN, f64::NAN));
            if !derivative(frequency, &mut port) || port.iter().any(|&z| !finite(z)) {
                return Err(SidebandError::CarrierDerivativeRefused {
                    sideband: k,
                    angular_frequency: frequency,
                });
            }
            poll(&mut cancelled)?;
            original_derivative[band * d * d..(band + 1) * d * d].copy_from_slice(&port);
            for row in 0..d {
                poll(&mut cancelled)?;
                for col in 0..d {
                    blocks[band * d * d + row * d + col] =
                        scale_entry(port[row * d + col], row, col, self.scales)?;
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedSidebandCarrierDerivative {
            factor: self,
            derivative: blocks,
            original_derivative,
            controls,
        })
    }
}

impl PreparedSidebandCarrierDerivative<'_, '_> {
    /// Complete original state dimension.
    pub fn dimension(&self) -> usize {
        self.factor.dimension
    }
    /// Retained signed sideband order.
    pub fn harmonics(&self) -> usize {
        self.factor.controls.harmonics
    }
    /// Signed carrier with respect to which coefficients are differentiated.
    pub fn carrier_angular_frequency(&self) -> f64 {
        self.factor.carrier
    }

    /// Solve the primal and exact differentiated retained equations with ONE LU.
    ///
    /// Both supplied arrays use `(k+N)*dimension+row` ordering. The second is
    /// the complete physical `dF_k/dnu` for this same source family; supply zeros
    /// only for a source whose complex phasors are actually held fixed. The base
    /// orbit, clock and tangent A(t) stay fixed, hence
    /// `D(nu) X' = F' - blockdiag(P'(i*(nu+k*omega))) X`.
    /// A total derivative along a changing nonlinear operating orbit is different
    /// and is not supplied by this operation.
    ///
    /// The existing native solve/refinement kernel serves both systems. A second
    /// original-row contraction checks `D X' + P' X - F'` independently of the
    /// rounded formed RHS. The differentiated omitted bands receive their own
    /// unit-bearing allowance. No tangent, port or derivative callbacks occur
    /// after preparation; exact singularities retain their original refusal.
    ///
    /// These are derivatives of the FINITE retained response. Small backward
    /// errors do not establish a forward-sensitivity bound near a resonance.
    /// Resolve sideband order and tangent/base sampling independently when
    /// using the derivative for peak location, bandwidth or optimization.
    ///
    /// # Errors
    /// Invalid full forcing/derivative arrays, either solve's numerical or
    /// omitted-band refusal, unrepresentable arithmetic, or cancellation. Neither
    /// a partial pair nor an accepted primal alone is published on failure.
    pub fn solve(
        &self,
        forcing: &[C64],
        forcing_angular_derivative: &[C64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SidebandCarrierDerivativeResponse> {
        self.solve_candidate(forcing, forcing_angular_derivative, true, cancelled)
    }

    fn solve_candidate(
        &self,
        forcing: &[C64],
        forcing_angular_derivative: &[C64],
        require_omitted_resolution: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SidebandCarrierDerivativeResponse> {
        poll(&mut cancelled)?;
        let d = self.dimension();
        let n = d * (2 * self.harmonics() + 1);
        if forcing.len() != n
            || forcing_angular_derivative.len() != n
            || forcing
                .iter()
                .chain(forcing_angular_derivative)
                .any(|&z| !finite(z))
        {
            return Err(invalid(
                "complete finite original forcing and carrier derivative required",
            ));
        }
        let response =
            self.factor
                .solve_candidate(forcing, require_omitted_resolution, &mut cancelled)?;
        let mut primal = Vec::with_capacity(n);
        let mut input = Vec::with_capacity(n);
        for i in 0..n {
            primal.push(normalize(
                response.coefficients[i],
                self.factor.scales.states[i % d],
            )?);
            input.push(normalize(
                forcing_angular_derivative[i],
                self.factor.scales.equations[i % d],
            )?);
        }
        let mut rhs = vec![C64::ZERO; n];
        for row in 0..n {
            poll(&mut cancelled)?;
            let band = row / d;
            let mut sum = ComplexSum::new();
            sum.re = Dd::from_f64(input[row].re);
            sum.im = Dd::from_f64(input[row].im);
            for col in 0..d {
                sum.product(-self.derivative[row * d + col], primal[band * d + col]);
            }
            rhs[row] = sum.value();
        }
        let mut angular_derivative = self.factor.solve_normalized(
            rhs,
            self.controls.maximum_relative_residual,
            require_omitted_resolution.then_some(self.controls.maximum_omitted_residual),
            &mut cancelled,
        )?;
        let mut slope = Vec::with_capacity(n);
        for (i, &value) in angular_derivative.coefficients.iter().enumerate() {
            slope.push(normalize(value, self.factor.scales.states[i % d])?);
        }
        let measured = original_residual(self, &primal, &slope, &input, &mut cancelled)?;
        angular_derivative.relative_residual = angular_derivative.relative_residual.max(measured);
        if measured > self.controls.maximum_relative_residual {
            return Err(SidebandError::Residual {
                measured,
                limit: self.controls.maximum_relative_residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(SidebandCarrierDerivativeResponse {
            response,
            angular_derivative,
        })
    }
}

fn resource_costs(dimension: usize, harmonics: usize) -> (u128, u128) {
    let d = dimension as u128;
    let n = d * (2 * harmonics + 1) as u128;
    let entries = n * d;
    (
        48 * entries + 256 * n + 64 * d * d + 512,
        768 * entries + 4096 * n * n + 1024 * n,
    )
}

fn normalize(value: C64, scale: f64) -> Result<C64> {
    let normalized = C64::new(value.re / scale, value.im / scale);
    if !finite(normalized)
        || (value.re != 0.0 && normalized.re == 0.0)
        || (value.im != 0.0 && normalized.im == 0.0)
    {
        return Err(invalid(
            "normalized carrier derivative data is not representable",
        ));
    }
    Ok(normalized)
}

fn original_residual(
    owner: &PreparedSidebandCarrierDerivative<'_, '_>,
    primal: &[C64],
    slope: &[C64],
    input: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let n = primal.len();
    let d = owner.dimension();
    let mut maximum = 0.0_f64;
    for row in 0..n {
        poll(cancelled)?;
        let band = row / d;
        let mut sum = ComplexSum::new();
        let mut scale = AftSum::new(input[row].abs());
        for col in 0..n {
            let a = owner.factor.operator[row * n + col];
            sum.product(a, slope[col]);
            scale.add(a.abs() * slope[col].abs());
        }
        for col in 0..d {
            let a = owner.derivative[row * d + col];
            let x = primal[band * d + col];
            sum.product(a, x);
            scale.add(a.abs() * x.abs());
        }
        sum.re = sum.re - Dd::from_f64(input[row].re);
        sum.im = sum.im - Dd::from_f64(input[row].im);
        let defect = sum.value().abs();
        let scale = scale.total();
        if !defect.is_finite() || !scale.is_finite() {
            return Err(invalid(
                "original differentiated equation is not representable",
            ));
        }
        let ratio = if scale == 0.0 {
            if defect == 0.0 { 0.0 } else { f64::INFINITY }
        } else {
            defect / scale
        };
        maximum = maximum.max(ratio);
    }
    Ok(maximum)
}
