//! A full-waveform phase chart for autonomous harmonic balance.

use super::{
    AftSum, HbAnchor, HbBudget, HbScales, OrbitError, OrbitProblem, Packing, ScaledHbOrbit,
    scaling::AdmittedScales,
};

/// Solve an autonomous cycle from all physical Fourier coefficients and an
/// unknown angular frequency, without selecting a privileged state or harmonic.
///
/// The packed seed contains all DC states, then real and imaginary state blocks
/// for each positive harmonic, followed by angular frequency in rad/s. Its
/// complete oscillatory waveform defines a fixed phase hyperplane: displacement
/// from the seed is orthogonal to the seed's phase derivative in the normalized
/// Parseval state metric. Every original balance equation and every coefficient
/// remain free. Static, nodal or absent first-state fundamentals are admissible.
///
/// This uses the same AFT/Newton, exact-island-tangent, deflation and independent
/// sampling checks as [`super::solve_hb_scaled_seeded`]. The appended phase row
/// has an exact derivative and is never deflated. Its normalized defect enters
/// `scaled_residual`; the returned orbit's unscaled diagnostics contain only
/// the original physical balance. Coefficients and frequency retain their units.
///
/// The seed must contain a nonzero oscillatory reference. This is a local phase
/// chart for an isolated autonomous cycle, not an amplitude constraint, an
/// existence proof, a primitive-period test or a stability determination.
///
/// # Errors
/// A nonautonomous problem, invalid physical scales/seed/budget, or an
/// unrepresentable reference refuses before source callbacks. A DC-only seed is
/// [`OrbitError::TrivialCollapse`]. Subsequent failures retain the shared solver's
/// typed refusals, including an explicitly supplied tangent's refusal.
pub fn solve_hb_autonomous_scaled_seeded<P: OrbitProblem>(
    problem: &P,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
) -> Result<ScaledHbOrbit, OrbitError> {
    let reference = PhaseCondition::new(problem, &seed, budget, scales)?;
    solve_with_reference(problem, seed, budget, scales, &reference)
}

pub(super) fn solve_with_reference<P: OrbitProblem>(
    problem: &P,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
    reference: &PhaseCondition,
) -> Result<ScaledHbOrbit, OrbitError> {
    let anchor = HbAnchor::Autonomous {
        omega_guess: seed.last().copied().unwrap_or(f64::NAN),
    };
    let mut admitted = AdmittedScales::new(problem.dim(), anchor, &seed, budget, scales)?;
    admitted.rows.push(1.0);
    super::solve_hb_seeded_impl(
        problem,
        anchor,
        seed,
        budget,
        Some(&admitted),
        Some(reference),
    )
}

// Keep the initial reference across every adaptive order. New harmonics extend
// the waveform with zero reference entries, and the frequency is never a phase
// coordinate. Storage/work is linear in the initial coefficient count.
pub(super) struct PhaseCondition {
    reference: Vec<f64>,
    direction: Vec<f64>,
    units: Vec<f64>,
}

impl PhaseCondition {
    pub(super) fn new<P: OrbitProblem>(
        problem: &P,
        seed: &[f64],
        budget: &HbBudget,
        scales: &HbScales,
    ) -> Result<Self, OrbitError> {
        if !problem.autonomous() {
            return Err(OrbitError::BadParameter {
                what: "full-waveform phase requires an autonomous problem",
            });
        }
        let anchor = HbAnchor::Autonomous {
            omega_guess: seed.last().copied().unwrap_or(f64::NAN),
        };
        let admitted = AdmittedScales::new(problem.dim(), anchor, seed, budget, scales)?;
        let pack = Packing {
            d: problem.dim(),
            n_harm: budget.harmonics,
            omega_slot: Some(seed.len() - 1),
        };
        let nc = pack.coeff_len();
        let units = admitted.columns[..nc].to_vec();
        let mut reference = vec![0.0; nc];
        let mut largest = 0.0_f64;
        for i in pack.d..nc {
            let value = seed[i] / units[i];
            if !value.is_finite() || (value == 0.0 && seed[i] != 0.0) {
                return Err(OrbitError::BadParameter {
                    what: "normalized autonomous phase reference must be representable",
                });
            }
            reference[i] = value;
            largest = largest.max(value.abs());
        }
        if largest == 0.0 {
            return Err(OrbitError::TrivialCollapse);
        }
        let mut direction = vec![0.0; nc];
        for n in 1..=pack.n_harm {
            for k in 0..pack.d {
                let re = pack.re_index(n, k);
                let im = pack.im_index(n, k);
                // Normalize before multiplying by harmonic order, so a finite
                // large physical reference cannot overflow its phase tangent.
                direction[re] = -(reference[im] / largest) * n as f64;
                direction[im] = (reference[re] / largest) * n as f64;
            }
        }
        let length = super::scaling::length(direction.iter().copied());
        for value in &mut direction {
            *value /= length;
        }
        Ok(Self {
            reference,
            direction,
            units,
        })
    }

    pub(super) fn value(&self, coefficients: &[f64]) -> f64 {
        let mut sum = AftSum::new(0.0);
        for (i, &direction) in self.direction.iter().enumerate() {
            if direction != 0.0 {
                sum.add(direction * (coefficients[i] / self.units[i] - self.reference[i]));
            }
        }
        sum.total()
    }

    pub(super) fn scaled_derivative(&self, column: usize, column_unit: f64) -> f64 {
        self.direction.get(column).map_or(0.0, |direction| {
            // Taking the unit ratio first avoids overflowing 1/unit for a small
            // but valid physical coordinate. Same-chart ratios are exactly one.
            direction * (column_unit / self.units[column])
        })
    }
}
