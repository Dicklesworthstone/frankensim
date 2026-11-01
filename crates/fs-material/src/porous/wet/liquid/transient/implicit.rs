//! Bounded backward Euler on the original column's source-qualified Darcy flux.

use super::{
    Hydraulics, LiquidCapillaryPressure, LiquidSaturationBoundaryStepInput,
    LiquidSaturationCandidate, LiquidSaturationColumn, LiquidSaturationStep,
    LiquidSaturationStepInput, WetPoreLiquidState, compensated_sum, hydraulics, interface_mobility,
    invalid, positive,
};
use crate::MaterialError;
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};
use fs_math::det;

/// Independent endpoint, whole-column conservation and nonlinear work budgets.
/// No tolerance expands a constitutive source's saturation/T/p/speed support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidSaturationImplicitControls {
    /// Absolute componentwise tolerance on the dimensionless saturation equation.
    pub absolute_saturation_residual_tolerance: f64,
    /// Relative to the larger actual storage increment or integrated flux term,
    /// never to the background saturation or a sum over other cells.
    pub relative_saturation_residual_tolerance: f64,
    /// One fixed original-to-final liquid-volume budget per transverse area [m].
    pub volume_balance_tolerance_m: f64,
    pub maximum_iterations: usize,
    /// Includes the full correction; successive trials halve its length.
    pub maximum_line_search_trials: usize,
}
impl LiquidSaturationImplicitControls {
    fn validate(self) -> Result<(), MaterialError> {
        if !positive(self.absolute_saturation_residual_tolerance)
            || self.absolute_saturation_residual_tolerance >= 1.0
            || !self.relative_saturation_residual_tolerance.is_finite()
            || self.relative_saturation_residual_tolerance < 0.0
            || self.relative_saturation_residual_tolerance >= 1.0
            || !positive(self.volume_balance_tolerance_m)
            || !(1..=256).contains(&self.maximum_iterations)
            || !(1..=64).contains(&self.maximum_line_search_trials)
        {
            return Err(invalid(
                "implicit liquid step needs positive absolute saturation/volume budgets, relative tolerance in [0,1), 1..256 iterations and 1..64 line-search trials",
            ));
        }
        Ok(())
    }
}

/// Actual bounded work and independently checked original endpoint equations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidSaturationImplicitReport {
    /// Accepted nonlinear corrections; an exact equilibrium needs none.
    pub iterations: usize,
    pub linear_solves: usize,
    /// Includes the original state and rejected or source-invalid trial states.
    pub state_evaluations: usize,
    pub line_search_backtracks: usize,
    /// Directions using only the capillary-pressure tangent when a full mobility
    /// tangent loses the diffusion signs. Accepted fluxes always use full laws.
    pub frozen_mobility_iterations: usize,
    pub maximum_absolute_residual_saturation: f64,
    /// Maximum componentwise residual divided by its own permitted tolerance.
    pub saturation_residual_ratio: f64,
    pub volume_balance_residual_m: f64,
}

/// Failure preserves the original column and its clock; no prefix is published.
#[derive(Debug)]
pub enum LiquidSaturationImplicitError {
    /// Invalid original input/source coordinates or unrepresentable arithmetic.
    Material(MaterialError),
    Tridiagonal(TridiagonalError),
    /// A supported iterate failed to converge within the explicit work budgets.
    NotConverged {
        report: LiquidSaturationImplicitReport,
        reason: &'static str,
    },
    /// Every candidate on a correction left the original source support or
    /// could not be represented. The last candidate refusal remains available.
    NoAdmissibleTrial {
        report: LiquidSaturationImplicitReport,
        error: MaterialError,
    },
    Cancelled,
}
impl core::fmt::Display for LiquidSaturationImplicitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit liquid saturation: {self:?}")
    }
}
impl std::error::Error for LiquidSaturationImplicitError {}
impl From<MaterialError> for LiquidSaturationImplicitError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<TridiagonalError> for LiquidSaturationImplicitError {
    fn from(error: TridiagonalError) -> Self {
        if error == TridiagonalError::Cancelled {
            Self::Cancelled
        } else {
            Self::Tridiagonal(error)
        }
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LiquidSaturationImplicitError> {
    if cancelled() {
        Err(LiquidSaturationImplicitError::Cancelled)
    } else {
        Ok(())
    }
}

struct Endpoint {
    saturations: Vec<f64>,
    states: Vec<WetPoreLiquidState>,
    flow: Hydraulics,
    residual: Vec<f64>,
    maximum_residual: f64,
    residual_ratio: f64,
    boundary_volume_m: f64,
    balance_residual_m: f64,
    merit: f64,
}
impl Endpoint {
    fn report_into(&self, report: &mut LiquidSaturationImplicitReport) {
        report.maximum_absolute_residual_saturation = self.maximum_residual;
        report.saturation_residual_ratio = self.residual_ratio;
        report.volume_balance_residual_m = self.balance_residual_m;
    }
}

fn endpoint(
    original: &LiquidSaturationColumn,
    saturations: Vec<f64>,
    input: LiquidSaturationBoundaryStepInput,
    controls: LiquidSaturationImplicitControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Endpoint, LiquidSaturationImplicitError> {
    poll(cancelled)?;
    let mut states = Vec::with_capacity(saturations.len());
    for (i, (cell, &s)) in original.cells.iter().zip(&saturations).enumerate() {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        states.push(cell.law.sample(s)?);
    }
    let flow = hydraulics(&original.cells, &states, input, || poll(cancelled))?;
    let mut residual = Vec::with_capacity(saturations.len());
    let mut maximum_residual: f64 = 0.0;
    let mut residual_ratio: f64 = 0.0;
    for i in 0..saturations.len() {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        let increment = saturations[i] - original.saturations[i];
        let transported =
            input.duration_s * (flow.faces[i] - flow.faces[i + 1]) / original.cells[i].capacity_m();
        let value = increment - transported;
        let tolerance = controls.absolute_saturation_residual_tolerance
            + controls.relative_saturation_residual_tolerance
                * increment.abs().max(transported.abs());
        let ratio = value.abs() / tolerance;
        if !value.is_finite() || !positive(tolerance) || !ratio.is_finite() {
            return Err(invalid("implicit liquid saturation residual is unrepresentable").into());
        }
        maximum_residual = maximum_residual.max(value.abs());
        residual_ratio = residual_ratio.max(ratio);
        residual.push(value);
    }
    // Sum actual changes directly: subtracting two large absolute inventories
    // would unnecessarily lose the small transfer being budgeted here.
    let change = compensated_sum(
        original
            .cells
            .iter()
            .enumerate()
            .map(|(i, cell)| cell.capacity_m() * (saturations[i] - original.saturations[i])),
    )?;
    let boundary_volume_m = input.duration_s * (flow.faces[0] - flow.faces[saturations.len()]);
    let balance_residual_m = change - boundary_volume_m;
    let merit = residual_ratio.max(balance_residual_m.abs() / controls.volume_balance_tolerance_m);
    if !boundary_volume_m.is_finite() || !balance_residual_m.is_finite() || !merit.is_finite() {
        return Err(invalid("implicit liquid volume balance is unrepresentable").into());
    }
    poll(cancelled)?;
    Ok(Endpoint {
        saturations,
        states,
        flow,
        residual,
        maximum_residual,
        residual_ratio,
        boundary_volume_m,
        balance_residual_m,
        merit,
    })
}

/// d(pl)/dS and d(ln(k*kr/mu))/dS on the original calibrated branch.
/// The latter may be unrepresentable near a wet endpoint; the pressure-only
/// iteration direction remains available without modifying the physical law.
fn tangents(state: &WetPoreLiquidState, mobility: f64) -> Result<(f64, f64), MaterialError> {
    let law = state.law();
    let s = state.liquid_saturation();
    let maximum = 1.0 - law.residual_gas_saturation();
    if mobility == 0.0 || s <= law.residual_liquid_saturation() || s >= maximum {
        return Ok((0.0, 0.0));
    }
    let span = maximum - law.residual_liquid_saturation();
    let remaining = (maximum - s) / span;
    let log_effective = if remaining < 0.5 {
        super::super::log_one_minus(remaining)
    } else {
        det::ln((s - law.residual_liquid_saturation()) / span)
    };
    let m = law.van_genuchten_m();
    let log_root = log_effective / m;
    let root = det::exp(log_root);
    let log_complement = if root < 0.5 {
        super::super::log_one_minus(root)
    } else {
        det::ln(-det::expm1(log_root))
    };
    let LiquidCapillaryPressure::Finite(pc) = state.capillary_pressure() else {
        return Err(invalid("mobile liquid pressure tangent is missing"));
    };
    let pressure = det::exp(
        det::ln(pc) + det::ln(1.0 - m)
            - det::ln(m)
            - det::ln(span)
            - log_effective
            - log_complement,
    );
    if !positive(pressure) {
        return Err(invalid(
            "liquid capillary-pressure tangent is unrepresentable",
        ));
    }
    let integral = -det::expm1(m * log_complement);
    let logarithmic_mobility = (law.connectivity_exponent()
        + 2.0 * root * det::exp((m - 1.0) * log_complement) / integral)
        / det::exp(log_effective)
        / span;
    Ok((pressure, logarithmic_mobility))
}

struct Linearization {
    lower: Vec<f64>,
    diagonal: Vec<f64>,
    upper: Vec<f64>,
    pressure_tangent: Vec<f64>,
    mobility_tangent: Vec<f64>,
}

fn storage_tangent(duration: f64, capacity: f64, flux_tangent: f64) -> f64 {
    // A sealed face remains exactly zero even when dt/capacity alone would
    // overflow. Preserve the ordinary scaling order for nonzero derivatives.
    if flux_tangent == 0.0 {
        0.0
    } else {
        duration / capacity * flux_tangent
    }
}

impl Linearization {
    fn new(n: usize) -> Self {
        Self {
            lower: vec![0.0; n - 1],
            diagonal: vec![1.0; n],
            upper: vec![0.0; n - 1],
            pressure_tangent: vec![0.0; n],
            mobility_tangent: vec![0.0; n],
        }
    }

    /// Full analytic Newton where the flux tangent retains diffusion signs;
    /// otherwise a monotone frozen-mobility pressure direction. Both directions
    /// solve the SAME complete nonlinear residual, with endpoint line search.
    fn assemble(
        &mut self,
        original: &LiquidSaturationColumn,
        point: &Endpoint,
        input: LiquidSaturationBoundaryStepInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, LiquidSaturationImplicitError> {
        let mut full = true;
        for i in 0..original.cells.len() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let (pressure, mobility) = tangents(&point.states[i], point.flow.mobility[i])?;
            self.pressure_tangent[i] = pressure;
            self.mobility_tangent[i] = mobility;
            full &= mobility.is_finite();
        }
        if full && self.fill(original, point, input, true, cancelled)? {
            return Ok(false);
        }
        self.fill(original, point, input, false, cancelled)?;
        Ok(true)
    }

    fn fill(
        &mut self,
        original: &LiquidSaturationColumn,
        point: &Endpoint,
        input: LiquidSaturationBoundaryStepInput,
        full: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, LiquidSaturationImplicitError> {
        self.lower.fill(0.0);
        self.diagonal.fill(1.0);
        self.upper.fill(0.0);
        let n = original.cells.len();
        let mut monotone = true;
        for face in 0..=n {
            if face % 128 == 0 {
                poll(cancelled)?;
            }
            if face == 0 || face == n {
                let i = if face == 0 { 0 } else { n - 1 };
                let sign = if face == 0 { -1.0 } else { 1.0 };
                let boundary = &point.flow.boundaries[usize::from(face == n)];
                let mut derivative = if boundary.pressure_conductance == 0.0 {
                    0.0
                } else {
                    sign * boundary.pressure_conductance * self.pressure_tangent[i]
                };
                if full && boundary.mobility_fraction != 0.0 {
                    derivative += point.flow.faces[face]
                        * boundary.mobility_fraction
                        * self.mobility_tangent[i];
                }
                monotone &= derivative.is_finite() && sign * derivative >= 0.0;
                self.diagonal[i] += sign
                    * storage_tangent(input.duration_s, original.cells[i].capacity_m(), derivative);
            } else {
                let left = face - 1;
                let right = face;
                let left_half = 0.5 * original.cells[left].thickness_m;
                let right_half = 0.5 * original.cells[right].thickness_m;
                let (mobility, shares) = interface_mobility(
                    point.flow.mobility[left],
                    point.flow.mobility[right],
                    left_half,
                    right_half,
                )?;
                // An exact sealed face has no pressure or mobility tangent.
                let pressure_factor = if mobility == 0.0 {
                    0.0
                } else {
                    mobility / (left_half + right_half)
                };
                let mut derivative_left = pressure_factor * self.pressure_tangent[left];
                let mut derivative_right = -pressure_factor * self.pressure_tangent[right];
                if full {
                    derivative_left +=
                        point.flow.faces[face] * shares[0] * self.mobility_tangent[left];
                    derivative_right +=
                        point.flow.faces[face] * shares[1] * self.mobility_tangent[right];
                }
                monotone &= derivative_left.is_finite()
                    && derivative_right.is_finite()
                    && derivative_left >= 0.0
                    && derivative_right <= 0.0;
                let capacity_left = original.cells[left].capacity_m();
                let capacity_right = original.cells[right].capacity_m();
                self.diagonal[left] +=
                    storage_tangent(input.duration_s, capacity_left, derivative_left);
                self.upper[left] +=
                    storage_tangent(input.duration_s, capacity_left, derivative_right);
                self.lower[left] -=
                    storage_tangent(input.duration_s, capacity_right, derivative_left);
                self.diagonal[right] -=
                    storage_tangent(input.duration_s, capacity_right, derivative_right);
            }
        }
        Ok(monotone)
    }
}

impl LiquidSaturationColumn {
    /// Solve one complete backward-Euler interval using the original nonlinear
    /// retention and mobility at the accepted endpoint. O(n) analytic assembly
    /// and adjacent-pivot tridiagonal LU avoid a dense or finite-difference solve.
    ///
    /// Source-supported backtracking bounds work explicitly. When the full
    /// mobility tangent is nonmonotone, a frozen-mobility pressure direction is
    /// used, but every trial resamples the full endpoint law and exact flux.
    /// Acceptance requires EVERY local saturation residual and one fixed whole-
    /// column volume balance. Dry/intrinsic seals remain exact. No clipping,
    /// epsilon mobility, hidden substeps, gas evolution or heat solve is applied.
    pub fn step_candidate_implicit(
        &self,
        input: LiquidSaturationStepInput,
        controls: LiquidSaturationImplicitControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        (LiquidSaturationCandidate, LiquidSaturationImplicitReport),
        LiquidSaturationImplicitError,
    > {
        self.step_candidate_implicit_with_boundaries(input.into(), controls, cancelled)
    }

    /// Solve the same fully implicit interval with sealed, prescribed-flow or
    /// pressure/resistive-pressure boundaries. Surface and porous half-cell
    /// resistances add in series; the exact interface mobility sensitivity
    /// enters the native tridiagonal Jacobian. Prescribed flux has zero
    /// boundary tangent but retains actual pressure/speed source admission.
    /// Local residual, whole-column volume and bounded-work controls are
    /// unchanged; no rejected trial or partial interval is published.
    pub fn step_candidate_implicit_with_boundaries(
        &self,
        input: LiquidSaturationBoundaryStepInput,
        controls: LiquidSaturationImplicitControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (LiquidSaturationCandidate, LiquidSaturationImplicitReport),
        LiquidSaturationImplicitError,
    > {
        poll(&mut cancelled)?;
        input.validate()?;
        controls.validate()?;
        let elapsed_s = self.elapsed_s + input.duration_s;
        if !elapsed_s.is_finite() || elapsed_s <= self.elapsed_s {
            return Err(invalid("transient liquid time must advance representably").into());
        }
        let mut report = LiquidSaturationImplicitReport {
            iterations: 0,
            linear_solves: 0,
            state_evaluations: 1,
            line_search_backtracks: 0,
            frozen_mobility_iterations: 0,
            maximum_absolute_residual_saturation: 0.0,
            saturation_residual_ratio: 0.0,
            volume_balance_residual_m: 0.0,
        };
        // An unsupported original source state is terminal, unlike an
        // unsupported trial direction that can be shortened without clipping.
        let mut point = endpoint(
            self,
            self.saturations.clone(),
            input,
            controls,
            &mut cancelled,
        )?;
        let before_pressures = point.flow.pressures.clone();
        let before = self.liquid_volume_per_area_m()?;
        point.report_into(&mut report);
        let n = self.cells.len();
        let mut matrix = Linearization::new(n);
        let mut workspace = TridiagonalWorkspace::new(n)?;
        let mut rhs = vec![0.0; n];
        let mut correction = vec![0.0; n];
        while point.merit > 1.0 {
            poll(&mut cancelled)?;
            if report.iterations == controls.maximum_iterations {
                return Err(LiquidSaturationImplicitError::NotConverged {
                    report,
                    reason: "nonlinear iteration budget exhausted",
                });
            }
            if matrix.assemble(self, &point, input, &mut cancelled)? {
                report.frozen_mobility_iterations += 1;
            }
            workspace.factor(
                &matrix.lower,
                &matrix.diagonal,
                &matrix.upper,
                &mut cancelled,
            )?;
            for (value, &residual) in rhs.iter_mut().zip(&point.residual) {
                *value = -residual;
            }
            workspace.solve_factored_into(&rhs, &mut correction, &mut cancelled)?;
            report.linear_solves += 1;
            let mut accepted = None;
            let mut last_refusal = None;
            let mut admissible = false;
            let mut alpha: f64 = 1.0;
            for trial in 0..controls.maximum_line_search_trials {
                poll(&mut cancelled)?;
                if trial != 0 {
                    report.line_search_backtracks += 1;
                }
                let trial_saturation = point
                    .saturations
                    .iter()
                    .zip(&correction)
                    .map(|(&s, &delta)| alpha.mul_add(delta, s))
                    .collect();
                report.state_evaluations += 1;
                match endpoint(self, trial_saturation, input, controls, &mut cancelled) {
                    Ok(candidate) => {
                        admissible = true;
                        if candidate.merit <= 1.0 || candidate.merit < point.merit {
                            accepted = Some(candidate);
                            break;
                        }
                    }
                    Err(LiquidSaturationImplicitError::Material(error)) => {
                        last_refusal = Some(error);
                    }
                    Err(error) => return Err(error),
                }
                alpha *= 0.5;
            }
            let Some(next) = accepted else {
                if !admissible {
                    if let Some(error) = last_refusal {
                        return Err(LiquidSaturationImplicitError::NoAdmissibleTrial {
                            report,
                            error,
                        });
                    }
                }
                return Err(LiquidSaturationImplicitError::NotConverged {
                    report,
                    reason: "bounded source-supported line search made no residual progress",
                });
            };
            point = next;
            report.iterations += 1;
            point.report_into(&mut report);
        }
        let next = Self {
            cells: self.cells.clone(),
            saturations: point.saturations,
            states: point.states,
            elapsed_s,
        };
        let after = next.liquid_volume_per_area_m()?;
        poll(&mut cancelled)?;
        Ok((
            LiquidSaturationCandidate {
                state: next,
                step: LiquidSaturationStep {
                    duration_s: input.duration_s,
                    face_flux_m_s: point.flow.faces,
                    liquid_pressure_pa_before: before_pressures,
                    liquid_pressure_pa_after: point.flow.pressures,
                    liquid_volume_per_area_before_m: before,
                    liquid_volume_per_area_after_m: after,
                    boundary_volume_per_area_m: point.boundary_volume_m,
                    balance_residual_m: point.balance_residual_m,
                },
            },
            report,
        ))
    }
}
