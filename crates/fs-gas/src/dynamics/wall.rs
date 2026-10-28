//! One-sided thermally perfect Euler wall Riemann problem.
//!
//! A stopped state is connected to incoming gas by a Hugoniot shock and to
//! retreating gas by an isentropic rarefaction. Its pressure is the wall stress;
//! no HLL artificial tensile stress, caloric extrapolation, or state floor is
//! substituted. The fixed, adiabatic wall transfers neither mass nor energy.

use super::{Conserved, Error, FaceFlux, FrozenGas, Primitive, poll, positive};
use fs_math::det;
use fs_qty::{Pressure, Temperature};

#[derive(Clone, Copy)]
struct WallState {
    pressure: f64,
    temperature: f64,
    density: f64,
    signal_speed: f64,
}

impl FrozenGas {
    /// Reflect one interior trace against a stationary adiabatic wall.
    /// `outward_normal` is -1 at the left face and +1 at the right face.
    pub(super) fn wall_flux(
        &self,
        inside: Conserved,
        primitive: Primitive,
        outward_normal: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FaceFlux, Error> {
        poll(cancelled)?;
        if !self.admissible_state(inside)
            || !matches!(outward_normal, -1.0 | 1.0)
            || !positive(primitive.density)
            || !positive(primitive.pressure)
            || !positive(primitive.temperature)
            || !positive(primitive.sound_speed)
            || !primitive.velocity.is_finite()
        {
            return Err(Error::Invalid("wall normal or interior trace"));
        }
        let wall = self.stopped_wall(primitive, outward_normal * primitive.velocity, cancelled)?;
        if !positive(wall.density) || !positive(wall.temperature) {
            return Err(Error::Invalid(
                "unrepresentable wall density or temperature",
            ));
        }
        poll(cancelled)?;
        Ok(FaceFlux {
            mass: 0.0,
            momentum: wall.pressure,
            energy: 0.0,
            signal_speed: wall.signal_speed,
            used_hll_fallback: false,
        })
    }

    fn stopped_wall(
        &self,
        p: Primitive,
        toward_wall: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallState, Error> {
        if toward_wall == 0.0 {
            return Ok(WallState {
                pressure: p.pressure,
                temperature: p.temperature,
                density: p.density,
                signal_speed: p.sound_speed,
            });
        }
        let tolerance = self.domain().absolute_energy_tolerance
            / (p.sound_speed + toward_wall.abs())
            + self.domain().relative_energy_tolerance * (p.sound_speed + toward_wall.abs())
            + 32.0 * f64::EPSILON * (p.sound_speed + toward_wall.abs());
        let (temperature, pressure, wave_speed) = if toward_wall > 0.0 {
            self.compressed_wall(p, toward_wall, tolerance, cancelled)?
        } else {
            self.expanded_wall(p, -toward_wall, tolerance, cancelled)?
        };
        let density = (pressure / self.gas_constant()) / temperature;
        let wall = self.from_rho_u_p(density, 0.0, pressure)?;
        let signal_speed = (p.velocity.abs() + p.sound_speed)
            .max(wall.sound_speed)
            .max(wave_speed);
        if !positive(pressure) || !positive(density) || !positive(signal_speed) {
            return Err(Error::Invalid("unrepresentable wall Riemann state"));
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
    ) -> Result<(f64, f64, f64), Error> {
        let r = self.gas_constant();
        let t1 = p.temperature;
        let maximum = self.domain().maximum_temperature;
        if t1 >= maximum {
            return Err(Error::CaloricPath(
                "wall compression exceeds the source temperature path",
            ));
        }
        let cp1 = self
            .source()
            .evaluate_properties(Temperature::new(t1), Pressure::new(p.pressure))?;
        let cv1 = cp1.mass_specific_cv().value();
        // w is the positive shock speed INTO the domain. Factoring out q
        // prevents the unchanged upstream state becoming a spurious root.
        let evaluate = |w: f64| -> Result<(f64, f64, f64), Error> {
            let shape = w / r - t1 / (w + q);
            let t = (t1 + q * shape).max(t1);
            if t > maximum * (1.0 + 4.0 * f64::EPSILON) {
                return Err(Error::CaloricPath(
                    "wall shock temperature outside source path",
                ));
            }
            let t = t.min(maximum); // endpoint rounding in the bounded parameterization
            let actual_dt = t - t1;
            let cv_secant = if actual_dt == 0.0 {
                cv1
            } else {
                self.source()
                    .specific_internal_energy_change(Temperature::new(t1), Temperature::new(t))?
                    .value()
                    / actual_dt
            };
            let cv = self
                .source()
                .evaluate_properties(Temperature::new(t), Pressure::new(p.pressure))?
                .mass_specific_cv()
                .value();
            let f = cv_secant * shape - 0.5 * q - r * t1 / (w + q);
            let derivative = cv / r + (cv + r) * t1 / ((w + q) * (w + q));
            if !positive(cv_secant) || !positive(derivative) || !f.is_finite() {
                return Err(Error::CaloricPath("wall shock caloric derivative"));
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
            return Err(Error::CaloricPath(
                "wall shock has no supported compression bracket",
            ));
        }
        let mut hi = (p.sound_speed + q).min(upper);
        for _ in 0..self.domain().maximum_iterations {
            poll(cancelled)?;
            if evaluate(hi)?.0 >= 0.0 {
                break;
            }
            if hi == upper {
                return Err(Error::CaloricPath(
                    "stopped shock state exceeds source temperature maximum",
                ));
            }
            hi = (2.0 * hi).min(upper);
        }
        if evaluate(hi)?.0 < 0.0 {
            return Err(Error::TemperatureIterations);
        }
        let mut w = (p.sound_speed + 0.25 * q).clamp(lo, hi);
        for _ in 0..self.domain().maximum_iterations {
            poll(cancelled)?;
            let (f, derivative, t) = evaluate(w)?;
            if f.abs() <= tolerance {
                let pressure = p.pressure + p.density * q * (q + w);
                let entropy = self
                    .source()
                    .specific_entropy_change(
                        Temperature::new(t1),
                        Pressure::new(p.pressure),
                        Temperature::new(t),
                        Pressure::new(pressure),
                    )?
                    .value();
                let entropy_tolerance = 8.0 * q * tolerance / t1 + 128.0 * f64::EPSILON * r;
                if !entropy.is_finite() || entropy < -entropy_tolerance {
                    return Err(Error::CaloricPath(
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
        Err(Error::TemperatureIterations)
    }

    fn expanded_wall(
        &self,
        p: Primitive,
        speed: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64, f64), Error> {
        let mut lo = det::sqrt(self.domain().minimum_temperature);
        let mut hi = det::sqrt(p.temperature);
        if self.rarefaction_integral(
            self.domain().minimum_temperature,
            p.temperature,
            0.25 * tolerance,
            cancelled,
        )? + tolerance
            < speed
        {
            // A positive Tmin cannot authorize extrapolation to vacuum.
            return Err(Error::CaloricPath(
                "wall rarefaction leaves the supported cold source domain",
            ));
        }
        let mut z = (hi - speed / self.rarefaction_integrand(p.temperature)?).max(lo);
        for _ in 0..self.domain().maximum_iterations {
            poll(cancelled)?;
            let t = (z * z).clamp(self.domain().minimum_temperature, p.temperature);
            let residual =
                self.rarefaction_integral(t, p.temperature, 0.25 * tolerance, cancelled)? - speed;
            if residual.abs() <= tolerance {
                let ds = self
                    .source()
                    .specific_entropy_change(
                        Temperature::new(p.temperature),
                        Pressure::new(p.pressure),
                        Temperature::new(t),
                        Pressure::new(p.pressure),
                    )?
                    .value();
                let pressure = p.pressure * det::exp(ds / self.gas_constant());
                return Ok((t, pressure, p.sound_speed + speed));
            }
            if residual > 0.0 {
                lo = z;
            } else {
                hi = z;
            }
            let trial = z + residual / self.rarefaction_integrand(t)?;
            z = if trial > lo && trial < hi {
                trial
            } else {
                lo + 0.5 * (hi - lo)
            };
        }
        Err(Error::TemperatureIterations)
    }

    // T=z² makes (Cp/c)dT = 2 sqrt(Cp Cv/R) dz. The constant-Cp
    // limit is integrated exactly and the cold 1/sqrt(T) singularity disappears.
    fn rarefaction_integrand(&self, t: f64) -> Result<f64, Error> {
        let properties = self
            .source()
            .evaluate_properties(Temperature::new(t), Pressure::new(1.0))?;
        let value = 2.0
            * det::sqrt(
                properties.mass_specific_cp().value()
                    * (properties.mass_specific_cv().value() / self.gas_constant()),
            );
        if !positive(value) {
            return Err(Error::CaloricPath("invalid rarefaction caloric integrand"));
        }
        Ok(value)
    }

    fn rarefaction_integral(
        &self,
        low: f64,
        high: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
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
            let mut end = high;
            // Cp may jump at a reconciled h/s seam: never smooth across it.
            for component in self.source().components() {
                for region in component.regions() {
                    let seam = region.temperature_min().value();
                    if seam > start && seam < end {
                        end = seam;
                    }
                }
            }
            let a = det::sqrt(start);
            let b = det::sqrt(end);
            let mid = a + 0.5 * (b - a);
            let fa = self.rarefaction_integrand(start)?;
            let fm = self.rarefaction_integrand((mid * mid).clamp(start, end))?;
            let fb = self.rarefaction_integrand(if end < high { end.next_down() } else { end })?;
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
            return Err(Error::Invalid("wall rarefaction integral overflow"));
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
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        let m = a + 0.5 * (b - a);
        let l = a + 0.5 * (m - a);
        let r = m + 0.5 * (b - m);
        let fl = self.rarefaction_integrand(l * l)?;
        let fr = self.rarefaction_integrand(r * r)?;
        let left = (m - a) * (values[0] + 4.0 * fl + values[1]) / 6.0;
        let right = (b - m) * (values[1] + 4.0 * fr + values[2]) / 6.0;
        let defect = left + right - estimate;
        if defect.abs() <= 15.0 * tolerance {
            return Ok(left + right + defect / 15.0);
        }
        if depth == 12 || m <= a || m >= b {
            return Err(Error::CaloricPath(
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

fn hypot(a: f64, b: f64) -> f64 {
    let scale = a.abs().max(b.abs());
    if scale == 0.0 {
        0.0
    } else {
        scale * det::sqrt((a / scale).powi(2) + (b / scale).powi(2))
    }
}

#[cfg(test)]
#[path = "wall_tests.rs"]
mod tests;
