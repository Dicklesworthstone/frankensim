//! Kirchhoff/van't Hoff integration of the original phase Cp difference.
//!
//! For a neutral ideal solute, integrating `d ln(c)/dT = ΔH(T)/(R T²)` gives
//! `R ln(c/c_ref) = ΔH_ref (1/T_ref - 1/T)
//!                 + ∫[T_ref,T] ΔCp(u) (1/u - 1/T) du`.
//! This is the variable-enthalpy correction, not a second fitted solubility law.
//! The original piecewise-linear curves are merged once, preserving all knots.
//! Local analytic kernels integrate their difference without subtracting two
//! large phase entropies or cancelling first-order terms near the reference.
//! Kirchhoff and temperature correction convention:
//! <https://www.oecd-nea.org/upload/docs/application/pdf/2020-09/tdb4new_2020-09-02_11-54-56_843.pdf>.

use super::{MolarHeatCapacityCurveV1, SolubilityErrorV1, SpeciesId, Temperature};
use fs_math::det::ln;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct PhaseHeatCapacities {
    curves: [MolarHeatCapacityCurveV1; 2],
    difference: Box<[(f64, f64)]>,
    identical: bool,
}
impl PhaseHeatCapacities {
    pub(super) fn new(
        solute: &SpeciesId,
        domain: [Temperature; 2],
        solid: MolarHeatCapacityCurveV1,
        dissolved: MolarHeatCapacityCurveV1,
    ) -> Result<Self, SolubilityErrorV1> {
        for curve in [&solid, &dissolved] {
            let span = curve.temperature_domain();
            if curve.species() != solute
                || span[0].value() > domain[0].value()
                || span[1].value() < domain[1].value()
            {
                return Err(SolubilityErrorV1::Invalid(
                    "phase Cp species and coverage of the complete solubility domain required",
                ));
            }
        }
        let mut temperatures =
            Vec::with_capacity(solid.knots_si().len() + dissolved.knots_si().len() + 2);
        temperatures.extend(domain.map(Temperature::value));
        for curve in [&solid, &dissolved] {
            temperatures.extend(curve.knots_si().iter().filter_map(|&(t, _)| {
                (t > domain[0].value() && t < domain[1].value()).then_some(t)
            }));
        }
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        let mut difference = Vec::with_capacity(temperatures.len());
        for t in temperatures {
            let cp = |curve: &MolarHeatCapacityCurveV1| {
                curve
                    .heat_capacity_at(Temperature::new(t))
                    .map(|cp| cp.value())
                    .map_err(|_| SolubilityErrorV1::Unrepresentable("original phase heat capacity"))
            };
            difference.push((t, cp(&dissolved)? - cp(&solid)?));
        }
        let identical = difference.iter().all(|&(_, cp)| cp == 0.0);
        Ok(Self {
            curves: [solid, dissolved],
            difference: difference.into_boxed_slice(),
            identical,
        })
    }

    pub(super) const fn curves(&self) -> &[MolarHeatCapacityCurveV1; 2] {
        &self.curves
    }

    pub(super) fn heap_bytes(&self) -> usize {
        core::mem::size_of_val(self.difference.as_ref())
            + self
                .curves
                .iter()
                .map(MolarHeatCapacityCurveV1::heap_bytes)
                .sum::<usize>()
    }

    pub(super) fn enthalpy_increment(&self, from: f64, to: f64) -> Result<f64, SolubilityErrorV1> {
        if from == to || self.identical {
            return Ok(0.0);
        }
        let mut sum = Sum::default();
        self.segments(from.min(to), from.max(to), |a, b, cp| {
            sum.add(mean(cp) * (b - a))
        })?;
        Ok(if to >= from {
            sum.value()?
        } else {
            -sum.value()?
        })
    }

    /// The integrated Cp correction to R*ln(c/c_ref), in J/(mol K).
    pub(super) fn log_saturation_correction(
        &self,
        from: f64,
        to: f64,
    ) -> Result<f64, SolubilityErrorV1> {
        if from == to || self.identical {
            return Ok(0.0);
        }
        let warming = to > from;
        let mut sum = Sum::default();
        self.segments(from.min(to), from.max(to), |a, b, cp| {
            let ratio = (b - a) / a;
            let weights = local_weights(ratio, warming);
            sum.add(cp[0] * weights[0])?;
            sum.add(cp[1] * weights[1])?;
            // Complete this segment's kernel to the ACTUAL terminal T. For
            // cooling the reversed integration sign makes the kernel positive.
            // Group the ratios to avoid forming a large sensible enthalpy only
            // to divide it by temperature again.
            let span = if warming {
                ((b - a) / b) * ((to - b) / to)
            } else {
                ratio * ((a - to) / to)
            };
            sum.add(mean(cp) * span)
        })?;
        sum.value()
    }

    fn segments(
        &self,
        low: f64,
        high: f64,
        mut visit: impl FnMut(f64, f64, [f64; 2]) -> Result<(), SolubilityErrorV1>,
    ) -> Result<(), SolubilityErrorV1> {
        let start = self
            .difference
            .partition_point(|&(t, _)| t <= low)
            .saturating_sub(1);
        for pair in self.difference[start..].windows(2) {
            let a = low.max(pair[0].0);
            let b = high.min(pair[1].0);
            if b > a {
                let at = |t: f64| {
                    let weight = (t - pair[0].0) / (pair[1].0 - pair[0].0);
                    // Convex weighting also handles opposite-sign ΔCp values
                    // whose direct difference could exceed f64's range.
                    if pair[0].1 == pair[1].1 {
                        pair[0].1
                    } else {
                        (1.0 - weight) * pair[0].1 + weight * pair[1].1
                    }
                };
                visit(a, b, [at(a), at(b)])?;
            }
            if pair[1].0 >= high {
                break;
            }
        }
        Ok(())
    }
}

// The same-sign form retains subnormal differences without adding two large
// capacities. Opposite signs can be added safely before taking their mean.
fn mean(cp: [f64; 2]) -> f64 {
    if cp[0].is_sign_positive() == cp[1].is_sign_positive() {
        cp[0] + 0.5 * (cp[1] - cp[0])
    } else {
        0.5 * (cp[0] + cp[1])
    }
}

/// Integral weights of the two linear Cp endpoints over [a,b], with q=(b-a)/a.
/// The local terminal temperature is b for warming and a for cooling. Both
/// weights are O(q²): computing entropy minus enthalpy/T would erase them near
/// the reference, especially when ΔH_ref is zero.
fn local_weights(q: f64, warming: bool) -> [f64; 2] {
    if q <= 0.125 {
        // The omitted term is < q^25; at q=1/8 this is below 3e-23.
        // Horner evaluation factors out q², retaining near-reference curvature.
        let mut weights = [0.0; 2];
        for n in (2..=24).rev() {
            let k = f64::from(n);
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            let coefficients = if warming {
                [
                    (k - 1.0) * (k + 2.0) / (2.0 * k * (k + 1.0)),
                    (k - 1.0) / (2.0 * (k + 1.0)),
                ]
            } else {
                [1.0 / (k * (k + 1.0)), 1.0 / (k + 1.0)]
            };
            for i in 0..2 {
                weights[i] = weights[i] * q + sign * coefficients[i];
            }
        }
        return weights.map(|weight| (weight * q) * q);
    }
    let logarithm = ln(1.0 + q);
    let (constant, linear) = if warming {
        let fraction = q / (1.0 + q);
        (logarithm - fraction, 1.0 - logarithm / q - 0.5 * fraction)
    } else {
        (q - logarithm, 0.5 * q - 1.0 + logarithm / q)
    };
    [constant - linear, linear]
}

#[derive(Default)]
struct Sum {
    main: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), SolubilityErrorV1> {
        let next = self.main + value;
        if !value.is_finite() || !next.is_finite() {
            return Err(SolubilityErrorV1::Unrepresentable("phase Cp integral"));
        }
        self.correction += if self.main.abs() >= value.abs() {
            (self.main - next) + value
        } else {
            (value - next) + self.main
        };
        self.main = next;
        Ok(())
    }

    fn value(&self) -> Result<f64, SolubilityErrorV1> {
        let value = self.main + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(SolubilityErrorV1::Unrepresentable(
                "phase Cp sum correction",
            ))
        }
    }
}

#[cfg(test)]
mod tests;
