//! Intrinsic exposure coordinate for one fixed outer-radius population.
//!
//! F(u) = f X(u) + s u + d u²(3-2u), with X(u)=1-(1-u)³.
//! The coefficients contain the original radius, solid density and a/b ratio;
//! exposure has units mol s/m³. Concentration, inventories and physical time
//! remain the parent batch's responsibility.

use super::super::{ShrinkingCoreError, ShrinkingCoreReaction};
use fs_qty::Length;

type Result<T> = core::result::Result<T, ShrinkingCoreError>;

#[derive(Clone, Debug)]
pub(super) struct CoreExposure {
    film: f64,
    surface: f64,
    diffusion: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct CoreProgress {
    pub(super) penetration: f64,
    pub(super) conversion_increment: f64,
    pub(super) remaining_fraction: f64,
    pub(super) complete: bool,
}

impl CoreExposure {
    pub(super) fn new(model: &ShrinkingCoreReaction, radius: Length) -> Result<Self> {
        let scale = model.scale(radius)?;
        let coefficients = [
            scale * model.resistances[0] / 3.0,
            scale * model.resistances[1],
            scale * radius.value() * model.resistances[2] / 6.0,
        ];
        for (coefficient, resistance) in coefficients.iter().zip(model.resistances) {
            if !coefficient.is_finite()
                || *coefficient < 0.0
                || (resistance > 0.0 && *coefficient == 0.0)
            {
                return Err(ShrinkingCoreError::Arithmetic(
                    "intrinsic exposure coefficient range loss",
                ));
            }
        }
        let result = Self {
            film: coefficients[0],
            surface: coefficients[1],
            diffusion: coefficients[2],
        };
        result.remaining(0.0)?;
        Ok(result)
    }

    /// Evaluate the positive tail directly, never F(1)-F(u).
    pub(super) fn remaining(&self, u0: f64) -> Result<f64> {
        validate_penetration(u0)?;
        if u0 == 1.0 {
            return Ok(0.0);
        }
        let y = 1.0 - u0;
        let tail =
            self.film * y * y * y + self.surface * y + self.diffusion * y * y * (3.0 - 2.0 * y);
        if !tail.is_finite() || tail <= 0.0 {
            return Err(ShrinkingCoreError::Arithmetic(
                "positive remaining exposure range loss",
            ));
        }
        Ok(tail)
    }

    /// Invert only the incremental exposure. Every potential evaluation polls
    /// the parent's shared work/cancellation counter. At most max_iterations+3
    /// evaluations are needed; zero dose and completion need only one.
    pub(super) fn advance(
        &self,
        u0: f64,
        dose: f64,
        max_iterations: usize,
        tick: &mut impl FnMut() -> Result<()>,
    ) -> Result<CoreProgress> {
        validate_penetration(u0)?;
        if !dose.is_finite() || dose < 0.0 || !(1..=128).contains(&max_iterations) {
            return Err(ShrinkingCoreError::Invalid(
                "finite nonnegative exposure and 1..=128 inverse iterations required",
            ));
        }
        tick()?;
        let remaining = self.remaining(u0)?;
        let y0 = 1.0 - u0;
        let fraction0 = y0 * y0 * y0;
        if dose == 0.0 {
            return Ok(CoreProgress {
                penetration: u0,
                conversion_increment: 0.0,
                remaining_fraction: fraction0,
                complete: u0 == 1.0,
            });
        }
        if dose >= remaining {
            return Ok(CoreProgress {
                penetration: 1.0,
                conversion_increment: fraction0,
                remaining_fraction: 0.0,
                complete: true,
            });
        }

        let film = self.film * y0 * y0 * y0;
        let surface = self.surface * y0;
        let diffusion = self.diffusion * y0 * y0;
        // v is either fractional penetration of the remaining core, or its
        // remaining radius fraction. Switch to the latter near completion so
        // the inverse never subtracts nearly equal cumulative potentials.
        let tail = dose > remaining * 0.5;
        let target = if tail { remaining - dose } else { dose };
        let a = diffusion * u0 * 3.0;
        let b = diffusion * y0;
        let bounds = if tail {
            [
                (film, 3, 1.0),
                (surface, 1, 1.0),
                (diffusion * (3.0 - 2.0 * y0), 2, 3.0),
                (0.0, 1, 1.0),
            ]
        } else {
            [(film, 1, 3.0), (surface, 1, 1.0), (a, 1, 2.0), (b, 2, 3.0)]
        };
        let mut evaluate = |v: f64| -> Result<(f64, f64)> {
            tick()?;
            let (value, derivative) = if tail {
                let y = y0 * v;
                (
                    film * v * v * v + surface * v + diffusion * v * v * (3.0 - 2.0 * y),
                    film * v * v * 3.0 + surface + diffusion * v * (1.0 - y) * 6.0,
                )
            } else {
                let w = 1.0 - v;
                (
                    film * v * (3.0 - v * (3.0 - v))
                        + surface * v
                        + a * v * (2.0 - v)
                        + b * v * v * (3.0 - 2.0 * v),
                    film * w * w * 3.0 + surface + a * w * 2.0 + b * v * w * 6.0,
                )
            };
            if !value.is_finite() || value < 0.0 {
                return Err(ShrinkingCoreError::Arithmetic(
                    "intrinsic exposure evaluation range loss",
                ));
            }
            // An overflowing derivative only disables Newton. The positive
            // potential and its bracket still admit bounded bisection.
            Ok((value, derivative))
        };
        let v = invert(target, &bounds, max_iterations, &mut evaluate)?;
        let u1 = if tail {
            1.0 - y0 * v
        } else {
            y0.mul_add(v, u0)
        };
        if u1 <= u0 || u1 >= 1.0 {
            return Err(ShrinkingCoreError::Arithmetic(
                "positive partial core progress lost penetration resolution",
            ));
        }
        let y1 = 1.0 - u1;
        let increment = (u1 - u0) * (y0 * y0 + y0 * y1 + y1 * y1);
        let fraction1 = y1 * y1 * y1;
        if increment <= 0.0 || fraction1 <= 0.0 {
            return Err(ShrinkingCoreError::Arithmetic(
                "positive conversion or remaining core fraction range loss",
            ));
        }
        Ok(CoreProgress {
            penetration: u1,
            conversion_increment: increment,
            remaining_fraction: fraction1,
            complete: false,
        })
    }
}

fn validate_penetration(u: f64) -> Result<()> {
    if u.is_finite() && (0.0..=1.0).contains(&u) {
        Ok(())
    } else {
        Err(ShrinkingCoreError::Invalid("penetration must lie in [0,1]"))
    }
}

// Each positive basis is bounded by c*v^p and factor*c*v^p. These
// bounds start the inverse near its scale even for extremely small exposure;
// bisecting [0,1] alone would spend hundreds of iterations reaching it.
fn invert(
    target: f64,
    terms: &[(f64, u32, f64); 4],
    iterations: usize,
    evaluate: &mut impl FnMut(f64) -> Result<(f64, f64)>,
) -> Result<f64> {
    let count = terms.iter().filter(|(c, _, _)| *c > 0.0).count() as f64;
    let (mut lo, mut hi) = (1.0_f64, 1.0_f64);
    for &(coefficient, degree, factor) in terms {
        if coefficient == 0.0 {
            continue;
        }
        let root = if target >= coefficient {
            1.0
        } else {
            match degree {
                1 => target / coefficient,
                2 => target.sqrt() / coefficient.sqrt(),
                _ => target.cbrt() / coefficient.cbrt(),
            }
        };
        let divisor = match degree {
            1 => count * factor,
            2 => (count * factor).sqrt(),
            _ => (count * factor).cbrt(),
        };
        lo = lo.min(root / divisor);
        hi = hi.min(root);
    }
    lo *= 0.5;
    hi = (2.0 * hi).min(1.0);
    if hi == 0.0 || !hi.is_finite() {
        return Err(ShrinkingCoreError::Arithmetic(
            "positive exposure inverse range loss",
        ));
    }
    if evaluate(lo)?.0 > target || evaluate(hi)?.0 < target {
        return Err(ShrinkingCoreError::Arithmetic(
            "exposure inverse could not establish its monotone bracket",
        ));
    }
    let tolerance = 4.0 * (target.next_up() - target);
    let mut x = lo + (hi - lo) * 0.5;
    let (mut best, mut error) = (x, f64::INFINITY);
    for _ in 0..iterations {
        let (value, derivative) = evaluate(x)?;
        let residual = value - target;
        if residual.abs() < error {
            (best, error) = (x, residual.abs());
        }
        if residual.abs() <= tolerance {
            return Ok(x);
        }
        if residual > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
        let midpoint = lo + (hi - lo) * 0.5;
        if midpoint == lo || midpoint == hi {
            return Ok(best);
        }
        let newton = x - residual / derivative;
        x = if derivative.is_finite()
            && derivative > 0.0
            && newton > lo
            && newton < hi
            && newton != x
        {
            newton
        } else {
            midpoint
        };
    }
    Err(ShrinkingCoreError::Budget)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(value: f64, expected: f64) {
        assert!((value - expected).abs() <= 4e-14 * expected.abs().max(1e-300));
    }

    #[test]
    fn mixed_resistances_keep_original_core_history() {
        let curve = CoreExposure {
            film: 1.0 / 3.0,
            surface: 2.0,
            diffusion: 0.5,
        };
        close(curve.remaining(0.25).unwrap(), 33.0 / 16.0);
        let direct = curve
            .advance(0.25, 71.0 / 48.0, 96, &mut || Ok(()))
            .unwrap();
        close(direct.penetration, 0.75);
        close(direct.conversion_increment, 13.0 / 32.0);
        close(direct.remaining_fraction, 1.0 / 64.0);
        assert!(!direct.complete);
        let first = curve
            .advance(0.25, 37.0 / 48.0, 96, &mut || Ok(()))
            .unwrap();
        let resumed = curve
            .advance(first.penetration, 17.0 / 24.0, 96, &mut || Ok(()))
            .unwrap();
        close(resumed.penetration, direct.penetration);
        close(
            first.conversion_increment + resumed.conversion_increment,
            13.0 / 32.0,
        );
        let done = curve
            .advance(0.25, 33.0 / 16.0, 96, &mut || Ok(()))
            .unwrap();
        assert!(done.complete);
        assert_eq!(done.remaining_fraction, 0.0);
        close(done.conversion_increment, 27.0 / 64.0);
    }

    #[test]
    fn pure_diffusion_retains_tiny_initial_progress_and_film_keeps_positive_tail() {
        let diffusion = CoreExposure {
            film: 0.0,
            surface: 0.0,
            diffusion: 1.0,
        };
        let tiny = diffusion.advance(0.0, 1e-240, 96, &mut || Ok(())).unwrap();
        close(tiny.penetration, 1e-120 / 3.0_f64.sqrt());
        assert!(tiny.conversion_increment > 0.0);
        assert!(!tiny.complete);
        let film = CoreExposure {
            film: 1.0,
            surface: 0.0,
            diffusion: 0.0,
        };
        let almost = film
            .advance(0.0, 1.0_f64.next_down(), 96, &mut || Ok(()))
            .unwrap();
        assert!(!almost.complete);
        assert!(almost.remaining_fraction > 0.0);
        // The returned penetration has an absolute binary64 spacing: its
        // subtraction from one limits relative accuracy in this tiny core.
        let tail = 1.0 - 1.0_f64.next_down();
        assert!((almost.remaining_fraction / tail - 1.0).abs() < 1e-10);
        let identity = film.advance(0.5, 0.0, 96, &mut || Ok(())).unwrap();
        assert_eq!(identity.penetration, 0.5);
        assert_eq!(identity.conversion_increment, 0.0);
        assert_eq!(identity.remaining_fraction, 0.125);
        assert!(matches!(
            film.advance(0.5, 1e-240, 96, &mut || Ok(())),
            Err(ShrinkingCoreError::Arithmetic(_))
        ));
        assert_eq!(
            film.advance(0.0, 0.5, 96, &mut || Err(ShrinkingCoreError::Cancelled))
                .unwrap_err(),
            ShrinkingCoreError::Cancelled
        );
    }
}
