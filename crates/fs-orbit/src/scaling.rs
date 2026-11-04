//! Physical conditioning of the original harmonic-balance Newton solve.
//!
//! State and equation units are independent: an acoustic pressure, a mechanical
//! displacement and an electrical current must not share a numerical tolerance.

use super::{HbAnchor, HbBudget, HbOrbit, OrbitError, OrbitProblem, Packing, det};

/// Characteristic physical magnitudes for a harmonic-balance solve.
///
/// These values condition the equations; they do not change the physical model.
/// With tolerance `epsilon`, acceptance bounds the combined Parseval RMS norm
/// of equation defects divided by `equations` by `epsilon`. Large DC states
/// never relax this absolute normalized bound. Choose each equation magnitude
/// as its desired physical RMS tolerance divided by `epsilon`.
#[derive(Debug, Clone, PartialEq)]
pub struct HbScales {
    /// Positive finite magnitude in each original state coordinate's units.
    /// Conditions Newton columns and their finite-difference perturbations.
    pub states: Vec<f64>,
    /// Positive finite magnitude in each original balance equation's units.
    /// A general port's row units need not be state units divided by time.
    pub equations: Vec<f64>,
    /// Positive finite angular-frequency magnitude [rad/s], used for an
    /// unknown-frequency column. It never rescales the physical clock.
    pub angular_frequency: f64,
}

/// One physically scaled solve, retaining the original physical orbit.
#[derive(Debug, Clone)]
pub struct ScaledHbOrbit {
    /// Physical coefficients/frequency and unscaled residual diagnostics.
    pub orbit: HbOrbit,
    /// Independently sampled normalized Parseval residual. For a backbone,
    /// includes its relative first-harmonic amplitude error as an extra row.
    pub scaled_residual: f64,
    /// Normalized Parseval change of the balance on the finer AFT grid.
    pub scaled_sampling_defect: f64,
}

/// Solve from a physical packed seed using independently declared state and
/// equation magnitudes. The original AFT, Newton, line search, phase condition,
/// backbone constraint and finer-grid check remain the numerical owners.
///
/// DC has Parseval weight one; positive-frequency real/imaginary coefficients
/// have weight two. Both residual and sampling acceptance are absolute in this
/// normalized chart; neither static preload nor an unrelated coordinate can
/// hide a missing harmonic. Newton is row/column conditioned, and perturbations
/// use each coordinate's own magnitude rather than one implicit SI unit.
/// Returned Fourier coefficients, angular frequency and the supplied seed all
/// retain the original physical units. This entry adds no regularization to
/// the original solver's existing singular-Jacobian retry.
///
/// Unanchored autonomous equilibrium detection uses the normalized non-DC
/// state norm and the original 1e-9 threshold. Backbone amplitude admission is
/// relative to its prescribed positive physical amplitude. This is fixed-order
/// discrete balance, not an existence, infinite-tail or forward-error bound.
///
/// # Errors
/// Malformed, nonpositive, nonfinite or unrepresentable scales/seed/budgets
/// refuse before source callbacks; subsequent numerical failures retain the
/// original [`OrbitError`] and physical residual diagnostics.
pub fn solve_hb_scaled_seeded<P: OrbitProblem>(
    problem: &P,
    anchor: HbAnchor,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
) -> Result<ScaledHbOrbit, OrbitError> {
    let admitted = AdmittedScales::new(problem.dim(), anchor, &seed, budget, scales)?;
    super::solve_hb_seeded_impl(problem, anchor, seed, budget, Some(&admitted))
}

pub(super) struct AdmittedScales {
    pub(super) columns: Vec<f64>,
    pub(super) rows: Vec<f64>,
    coefficient_count: usize,
}

// Raw diagnostics also retain resolved small/large physical units: squaring
// before taking a norm could otherwise report zero or infinity for finite data.
pub(super) fn length(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0.0, det::hypot)
}

impl AdmittedScales {
    fn new(
        d: usize,
        anchor: HbAnchor,
        seed: &[f64],
        budget: &HbBudget,
        scales: &HbScales,
    ) -> Result<Self, OrbitError> {
        let invalid = || OrbitError::BadParameter {
            what: "scaled HB requires finite physical seed, positive state/equation/frequency scales and valid budgets",
        };
        let positive = |v: f64| v.is_finite() && v > 0.0;
        let nc = budget
            .harmonics
            .checked_mul(2)
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(d))
            .ok_or_else(invalid)?;
        let unknown = !matches!(anchor, HbAnchor::Forced { .. });
        let count = nc.checked_add(usize::from(unknown)).ok_or_else(invalid)?;
        if d == 0
            || budget.harmonics == 0
            || budget.max_newton == 0
            || !positive(budget.tolerance)
            || budget.tolerance >= 1.0
            || count.checked_mul(count).is_none()
            || scales.states.len() != d
            || scales.equations.len() != d
            || scales
                .states
                .iter()
                .chain(&scales.equations)
                .any(|&v| !positive(v))
            || !positive(scales.angular_frequency)
            || seed.len() != count
            || seed.iter().any(|v| !v.is_finite())
        {
            return Err(invalid());
        }
        super::admitted_backbone_amplitude(anchor, budget.tolerance)?;
        let omega = if unknown {
            seed[nc]
        } else if let HbAnchor::Forced { omega } = anchor {
            omega
        } else {
            unreachable!()
        };
        if !positive(omega) || !positive(super::TAU / omega) {
            return Err(invalid());
        }
        let pack = Packing {
            d,
            n_harm: budget.harmonics,
            omega_slot: unknown.then_some(nc),
        };
        let mut columns = vec![0.0; count];
        let mut rows = vec![0.0; nc];
        for k in 0..d {
            columns[k] = scales.states[k];
            rows[k] = scales.equations[k];
            for n in 1..=budget.harmonics {
                for i in [pack.re_index(n, k), pack.im_index(n, k)] {
                    columns[i] = scales.states[k] / core::f64::consts::SQRT_2;
                    rows[i] = scales.equations[k] / core::f64::consts::SQRT_2;
                }
            }
        }
        if unknown {
            columns[nc] = scales.angular_frequency;
        }
        if let HbAnchor::Backbone { amplitude, .. } = anchor {
            rows.push(0.5 * amplitude);
        }
        if columns.iter().chain(&rows).any(|&v| !positive(v)) {
            return Err(invalid());
        }
        Ok(Self {
            columns,
            rows,
            coefficient_count: nc,
        })
    }

    pub(super) fn row_norm(&self, values: &[f64]) -> f64 {
        values
            .iter()
            .zip(&self.rows)
            .fold(0.0, |norm, (&value, &unit)| det::hypot(norm, value / unit))
    }

    pub(super) fn balance_norm(&self, balance: &[f64], amplitude: f64) -> f64 {
        let norm = self.row_norm(balance);
        self.rows
            .get(self.coefficient_count)
            .map_or(norm, |unit| det::hypot(norm, amplitude / unit))
    }

    pub(super) fn difference_norm(&self, a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b)
            .zip(&self.rows)
            .fold(0.0, |norm, ((a, b), unit)| det::hypot(norm, (a - b) / unit))
    }

    pub(super) fn oscillatory_norm(&self, u: &[f64], pack: &Packing) -> f64 {
        u[pack.d..pack.coeff_len()]
            .iter()
            .zip(&self.columns[pack.d..pack.coeff_len()])
            .fold(0.0, |norm, (value, unit)| det::hypot(norm, value / unit))
    }
}
