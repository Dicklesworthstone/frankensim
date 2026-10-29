//! One stopped-wall Riemann solver over a borrowed, fixed-composition caloric law.
//! Owners retain source admission, composition and typed source failures. No gas
//! owner, species vector or numerical workspace is allocated by this kernel.

use crate::dynamics::{CaloricDomain, Primitive};
use fs_math::det;

/// The owner's mass-specific caloric values, from the same source evaluation.
#[derive(Clone, Copy)]
pub(crate) struct HeatCapacities {
    pub cp: f64,
    pub cv: f64,
}

/// Minimal borrowed caloric interface needed by the physical wall problem.
/// Composition is fixed for one solve, while Cp and Cv evolve with temperature.
pub(crate) trait CaloricLaw {
    type Error;
    fn domain(&self) -> CaloricDomain;
    fn gas_constant(&self) -> f64;
    fn heat_capacities(
        &self,
        temperature: f64,
        pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeatCapacities, Self::Error>;
    fn internal_energy_change(
        &self,
        initial_temperature: f64,
        final_temperature: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Self::Error>;
    fn entropy_change(
        &self,
        initial_temperature: f64,
        initial_pressure: f64,
        final_temperature: f64,
        final_pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Self::Error>;
    fn sound_speed(
        &self,
        density: f64,
        pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Self::Error>;
    /// Construct a physical primitive with the owner's original conservative
    /// precision and source-path admission, retaining this local composition.
    fn primitive(
        &self,
        density: f64,
        velocity: f64,
        pressure: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Self::Error>;
    /// First original caloric seam strictly after `start`, or `end`.
    fn next_seam(
        &self,
        start: f64,
        end: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Self::Error>;
}

/// Source failures are returned unchanged to their original owner.
#[derive(Debug)]
pub(crate) enum WallError<E> {
    Source(E),
    Invalid(&'static str),
    CaloricPath(&'static str),
    Iterations,
    Cancelled,
}

#[derive(Clone, Copy)]
pub(crate) struct WallState {
    pub pressure: f64,
    pub temperature: f64,
    pub density: f64,
    pub signal_speed: f64,
}

pub(crate) struct WallSolver<'a, L: CaloricLaw> {
    law: &'a L,
}
impl<'a, L: CaloricLaw> WallSolver<'a, L> {
    pub const fn new(law: &'a L) -> Self {
        Self { law }
    }

    pub fn stopped_wall(
        &self,
        p: Primitive,
        toward_wall: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallState, WallError<L::Error>> {
        if toward_wall == 0.0 {
            return Ok(WallState {
                pressure: p.pressure,
                temperature: p.temperature,
                density: p.density,
                signal_speed: p.sound_speed,
            });
        }
        let tolerance = self.law.domain().absolute_energy_tolerance
            / (p.sound_speed + toward_wall.abs())
            + self.law.domain().relative_energy_tolerance * (p.sound_speed + toward_wall.abs())
            + 32.0 * f64::EPSILON * (p.sound_speed + toward_wall.abs());
        let (temperature, pressure, wave_speed) = if toward_wall > 0.0 {
            self.compressed_wall(p, toward_wall, tolerance, cancelled)?
        } else {
            self.expanded_wall(p, -toward_wall, tolerance, cancelled)?
        };
        let density = (pressure / self.law.gas_constant()) / temperature;
        let wall_sound_speed = self
            .law
            .sound_speed(density, pressure, cancelled)
            .map_err(WallError::Source)?;
        let signal_speed = (p.velocity.abs() + p.sound_speed)
            .max(wall_sound_speed)
            .max(wave_speed);
        if !positive(pressure) || !positive(density) || !positive(signal_speed) {
            return Err(WallError::Invalid("unrepresentable wall Riemann state"));
        }
        Ok(WallState {
            pressure,
            temperature,
            density,
            signal_speed,
        })
    }

    fn compressed_wall(
        &self,
        p: Primitive,
        q: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64, f64), WallError<L::Error>> {
        let r = self.law.gas_constant();
        let t1 = p.temperature;
        let maximum = self.law.domain().maximum_temperature;
        if t1 >= maximum {
            return Err(WallError::CaloricPath(
                "wall compression exceeds the source temperature path",
            ));
        }
        let cv1 = self
            .law
            .heat_capacities(t1, p.pressure, cancelled)
            .map_err(WallError::Source)?
            .cv;
        // w is the positive shock speed INTO the domain. Factoring out q
        // prevents the unchanged upstream state becoming a spurious root.
        let evaluate =
            |w: f64, cancelled: &mut _| -> Result<(f64, f64, f64), WallError<L::Error>> {
                let shape = w / r - t1 / (w + q);
                let t = (t1 + q * shape).max(t1);
                if t > maximum * (1.0 + 4.0 * f64::EPSILON) {
                    return Err(WallError::CaloricPath(
                        "wall shock temperature outside source path",
                    ));
                }
                let t = t.min(maximum); // endpoint rounding in the bounded parameterization
                let actual_dt = t - t1;
                let cv_secant = if actual_dt == 0.0 {
                    cv1
                } else {
                    self.law
                        .internal_energy_change(t1, t, cancelled)
                        .map_err(WallError::Source)?
                        / actual_dt
                };
                let cv = self
                    .law
                    .heat_capacities(t, p.pressure, cancelled)
                    .map_err(WallError::Source)?
                    .cv;
                let f = cv_secant * shape - 0.5 * q - r * t1 / (w + q);
                let derivative = cv / r + (cv + r) * t1 / ((w + q) * (w + q));
                if !positive(cv_secant) || !positive(derivative) || !f.is_finite() {
                    return Err(WallError::CaloricPath("wall shock caloric derivative"));
                }
                Ok((f, derivative, t))
            };
        let mut lo = 2.0 * r * t1 / (hypot(q, 2.0 * det::sqrt(r * t1)) + q);
        // The upper w bound maps exactly to the declared maximum temperature.
        let b = q - r * (maximum - t1) / q;
        let upper = if b == f64::NEG_INFINITY {
            f64::MAX
        } else {
            let h = hypot(b, 2.0 * det::sqrt(r * maximum));
            if b >= 0.0 {
                2.0 * r * maximum / (b + h)
            } else {
                0.5 * h - 0.5 * b
            }
        };
        if !positive(lo) || !positive(upper) || upper <= lo {
            return Err(WallError::CaloricPath(
                "wall shock has no supported compression bracket",
            ));
        }
        let mut hi = (p.sound_speed + q).min(upper);
        for _ in 0..self.law.domain().maximum_iterations {
            poll(cancelled)?;
            if evaluate(hi, cancelled)?.0 >= 0.0 {
                break;
            }
            if hi == upper {
                return Err(WallError::CaloricPath(
                    "stopped shock state exceeds source temperature maximum",
                ));
            }
            hi = (2.0 * hi).min(upper);
        }
        if evaluate(hi, cancelled)?.0 < 0.0 {
            return Err(WallError::Iterations);
        }
        let mut w = (p.sound_speed + 0.25 * q).clamp(lo, hi);
        for _ in 0..self.law.domain().maximum_iterations {
            poll(cancelled)?;
            let (f, derivative, t) = evaluate(w, cancelled)?;
            if f.abs() <= tolerance {
                let pressure = p.pressure + p.density * q * (q + w);
                let entropy = self
                    .law
                    .entropy_change(t1, p.pressure, t, pressure, cancelled)
                    .map_err(WallError::Source)?;
                let entropy_tolerance = 8.0 * q * tolerance / t1 + 128.0 * f64::EPSILON * r;
                if !entropy.is_finite() || entropy < -entropy_tolerance {
                    return Err(WallError::CaloricPath(
                        "wall compression has inadmissible entropy production",
                    ));
                }
                return Ok((t, pressure, w));
            }
            if f < 0.0 {
                lo = w;
            } else {
                hi = w;
            }
            let trial = w - f / derivative;
            w = if trial > lo && trial < hi {
                trial
            } else {
                lo + 0.5 * (hi - lo)
            };
        }
        Err(WallError::Iterations)
    }

    fn expanded_wall(
        &self,
        p: Primitive,
        speed: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64, f64), WallError<L::Error>> {
        let mut lo = det::sqrt(self.law.domain().minimum_temperature);
        let mut hi = det::sqrt(p.temperature);
        if self.rarefaction_integral(
            self.law.domain().minimum_temperature,
            p.temperature,
            0.25 * tolerance,
            cancelled,
        )? + tolerance
            < speed
        {
            return Err(WallError::CaloricPath(
                "wall rarefaction leaves the supported cold source domain",
            ));
        }
        let mut z = (hi - speed / self.rarefaction_integrand(p.temperature, cancelled)?).max(lo);
        for _ in 0..self.law.domain().maximum_iterations {
            poll(cancelled)?;
            let t = (z * z).clamp(self.law.domain().minimum_temperature, p.temperature);
            let residual =
                self.rarefaction_integral(t, p.temperature, 0.25 * tolerance, cancelled)? - speed;
            if residual.abs() <= tolerance {
                let ds = self
                    .law
                    .entropy_change(p.temperature, p.pressure, t, p.pressure, cancelled)
                    .map_err(WallError::Source)?;
                let pressure = p.pressure * det::exp(ds / self.law.gas_constant());
                return Ok((t, pressure, p.sound_speed + speed));
            }
            if residual > 0.0 {
                lo = z;
            } else {
                hi = z;
            }
            let trial = z + residual / self.rarefaction_integrand(t, cancelled)?;
            z = if trial > lo && trial < hi {
                trial
            } else {
                lo + 0.5 * (hi - lo)
            };
        }
        Err(WallError::Iterations)
    }

    // T=z² makes (Cp/c)dT = 2 sqrt(Cp Cv/R) dz. The constant-Cp
    // limit is integrated exactly and the cold 1/sqrt(T) singularity disappears.
    fn rarefaction_integrand(
        &self,
        t: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, WallError<L::Error>> {
        let properties = self
            .law
            .heat_capacities(t, 1.0, cancelled)
            .map_err(WallError::Source)?;
        let value = 2.0 * det::sqrt(properties.cp * (properties.cv / self.law.gas_constant()));
        if !positive(value) {
            return Err(WallError::CaloricPath(
                "invalid rarefaction caloric integrand",
            ));
        }
        Ok(value)
    }

    pub fn rarefaction_integral(
        &self,
        low: f64,
        high: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, WallError<L::Error>> {
        if low == high {
            return Ok(0.0);
        }
        let width = det::sqrt(high) - det::sqrt(low);
        if width == 0.0 {
            return Ok(0.0);
        }
        let mut start = low;
        let mut total = 0.0;
        while start < high {
            poll(cancelled)?;
            // Cp may jump at a reconciled h/s seam: never smooth across it.
            let end = self
                .law
                .next_seam(start, high, cancelled)
                .map_err(WallError::Source)?;
            let a = det::sqrt(start);
            let b = det::sqrt(end);
            let mid = a + 0.5 * (b - a);
            let fa = self.rarefaction_integrand(start, cancelled)?;
            let fm = self.rarefaction_integrand((mid * mid).clamp(start, end), cancelled)?;
            let fb = self
                .rarefaction_integrand(if end < high { end.next_down() } else { end }, cancelled)?;
            let estimate = (b - a) * (fa + 4.0 * fm + fb) / 6.0;
            total += self.rarefaction_panel(
                a,
                b,
                [fa, fm, fb],
                estimate,
                tolerance * ((b - a) / width),
                0,
                cancelled,
            )?;
            start = end;
        }
        if !total.is_finite() {
            return Err(WallError::Invalid("wall rarefaction integral overflow"));
        }
        Ok(total)
    }

    #[allow(clippy::too_many_arguments)]
    fn rarefaction_panel(
        &self,
        a: f64,
        b: f64,
        values: [f64; 3],
        estimate: f64,
        tolerance: f64,
        depth: u8,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, WallError<L::Error>> {
        poll(cancelled)?;
        let m = a + 0.5 * (b - a);
        let l = a + 0.5 * (m - a);
        let r = m + 0.5 * (b - m);
        let fl = self.rarefaction_integrand(l * l, cancelled)?;
        let fr = self.rarefaction_integrand(r * r, cancelled)?;
        let left = (m - a) * (values[0] + 4.0 * fl + values[1]) / 6.0;
        let right = (b - m) * (values[1] + 4.0 * fr + values[2]) / 6.0;
        let defect = left + right - estimate;
        if defect.abs() <= 15.0 * tolerance {
            return Ok(left + right + defect / 15.0);
        }
        if depth == 12 || m <= a || m >= b {
            return Err(WallError::CaloricPath(
                "wall rarefaction quadrature budget exhausted",
            ));
        }
        Ok(self.rarefaction_panel(
            a,
            m,
            [values[0], fl, values[1]],
            left,
            0.5 * tolerance,
            depth + 1,
            cancelled,
        )? + self.rarefaction_panel(
            m,
            b,
            [values[1], fr, values[2]],
            right,
            0.5 * tolerance,
            depth + 1,
            cancelled,
        )?)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), WallError<E>> {
    if cancelled() {
        Err(WallError::Cancelled)
    } else {
        Ok(())
    }
}
fn hypot(a: f64, b: f64) -> f64 {
    let scale = a.abs().max(b.abs());
    if scale == 0.0 {
        0.0
    } else {
        scale * det::sqrt((a / scale).powi(2) + (b / scale).powi(2))
    }
}
