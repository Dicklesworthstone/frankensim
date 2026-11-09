//! Bounded normal-wave scattering and positive local work.
//!
//! For exp(+i omega t), k=b-i*a, a>=0. Inward pressure amplitudes are
//! p(z)=A exp(-i k z)+B exp(-i k (d-z)); no growing exponential is formed.
//! Normal strain is -p/M, M=K+4G/3. Its cycle-mean bulk and shear heat are
//! omega*K''*|strain|²/2 and 2*omega*G''*|strain|²/3.

use fs_math::{c64::C64, dd::Dd, det};

#[path = "kernel/layers.rs"]
pub(crate) mod layers;

pub(crate) type Result<T> = core::result::Result<T, &'static str>;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Slab {
    pub(crate) omega: f64,
    pub(crate) length: f64,
    pub(crate) modulus: C64,
    pub(crate) moduli: [C64; 2],
    pub(crate) k: C64,
    pub(crate) impedance: C64,
    exterior: [f64; 2],
    pub(crate) scattering: [C64; 4],
    inward: [[C64; 2]; 2],
    parity_factor: C64,
    normalized_impedance: C64,
    normalized_exterior: [f64; 2],
    missing_one_way: C64,
    log_gram: [f64; 2],
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Drive {
    pub(crate) inward: [C64; 2],
    pub(crate) parity: [C64; 2],
    pub(crate) outgoing: [C64; 2],
    pub(crate) incident_power: [f64; 2],
    pub(crate) outgoing_power: [f64; 2],
    pub(crate) heat: [f64; 2],
    pub(crate) residual: f64,
    pub(crate) strain_bound: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Field {
    pub(crate) pressure: C64,
    pub(crate) velocity: C64,
    pub(crate) strain: C64,
    pub(crate) heat_density: [f64; 2],
}

pub(super) fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn divided(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}

fn from_log(log: f64) -> Result<f64> {
    let value = det::exp(log);
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err("viscoelastic wave power is unrepresentable")
    }
}

fn product_ratio(a: f64, b: f64, c: f64, d: f64) -> Result<f64> {
    if a == 0.0 || b == 0.0 || c == 0.0 {
        return Ok(0.0);
    }
    let ab = a * b;
    let quotient = ab / d;
    let value = quotient * c;
    if ab.is_normal() && quotient.is_normal() && value.is_normal() {
        Ok(value)
    } else {
        from_log(det::ln(a) + det::ln(b) + det::ln(c) - det::ln(d))
    }
}

fn propagation(k: C64, distance: f64) -> Result<C64> {
    let magnitude = det::exp(k.im * distance);
    if magnitude == 0.0 {
        return Ok(C64::ZERO);
    }
    let phase = k.re * distance;
    if !phase.is_finite() {
        return Err("viscoelastic wave phase is unrepresentable");
    }
    let (sine, cosine) = det::sin_cos(phase);
    Ok(C64::new(cosine, -sine).scale(magnitude))
}

// 1-exp(-i*k*d), retaining the small real loss and phase independently.
fn missing_propagation(k: C64, distance: f64) -> Result<C64> {
    let attenuation = det::exp(k.im * distance);
    if attenuation == 0.0 {
        return Ok(C64::ONE);
    }
    let phase = k.re * distance;
    if !phase.is_finite() {
        return Err("viscoelastic wave phase is unrepresentable");
    }
    let half_sine = det::sin_cos(0.5 * phase).0;
    Ok(C64::new(
        -det::expm1(k.im * distance) + 2.0 * attenuation * half_sine * half_sine,
        attenuation * det::sin_cos(phase).0,
    ))
}

fn weighted(a: C64, x: f64, b: C64, y: f64) -> C64 {
    let (x, y) = (Dd::from_f64(x), Dd::from_f64(y));
    C64::new(
        (Dd::from_f64(a.re) * x + Dd::from_f64(b.re) * y).to_f64(),
        (Dd::from_f64(a.im) * x + Dd::from_f64(b.im) * y).to_f64(),
    )
}

// Eigenvalues of the exact pressure Gram integral. For thin layers, evaluate
// sinhc(x)-1 and 1-sinc(y) separately so an odd standing wave retains O(d³).
fn wave_gram_logs(k: C64, d: f64) -> Result<[f64; 2]> {
    if d == 0.0 {
        return Ok([f64::NEG_INFINITY; 2]);
    }
    let x = -k.im * d;
    let y = k.re * d;
    let attenuation = det::exp(-x);
    let values = if x < 0.1 && y.abs() < 0.1 {
        let x2 = x * x;
        let y2 = y * y;
        let sh = x2
            * (1.0 / 6.0
                + x2 * (1.0 / 120.0
                    + x2 * (1.0 / 5040.0 + x2 * (1.0 / 362880.0 + x2 / 39916800.0))));
        let si = y2
            * (1.0 / 6.0
                - y2 * (1.0 / 120.0
                    - y2 * (1.0 / 5040.0 - y2 * (1.0 / 362880.0 - y2 / 39916800.0))));
        let radius = det::hypot(x, y);
        let odd = if radius == 0.0 {
            f64::NEG_INFINITY
        } else {
            let xn = x / radius;
            let yn = y / radius;
            let coefficient = xn
                * xn
                * (1.0 / 6.0
                    + x2 * (1.0 / 120.0
                        + x2 * (1.0 / 5040.0 + x2 * (1.0 / 362880.0 + x2 / 39916800.0))))
                + yn * yn
                    * (1.0 / 6.0
                        - y2 * (1.0 / 120.0
                            - y2 * (1.0 / 5040.0 - y2 * (1.0 / 362880.0 - y2 / 39916800.0))));
            det::ln(d) - x - det::ln(2.0) + 2.0 * det::ln(radius) + det::ln(coefficient)
        };
        [det::ln(d) - x + det::ln(1.0 + 0.5 * sh - 0.5 * si), odd]
    } else {
        let j = if x < 1.0e-4 {
            d * if x == 0.0 {
                1.0
            } else {
                -det::expm1(-2.0 * x) / (2.0 * x)
            }
        } else {
            (0.5 * -det::expm1(-2.0 * x)) / -k.im
        };
        let l = if attenuation == 0.0 {
            0.0
        } else if y.abs() < 1.0e-4 {
            (d * attenuation) * if y == 0.0 { 1.0 } else { det::sin_cos(y).0 / y }
        } else if y.is_finite() {
            attenuation * (det::sin_cos(y).0 / k.re)
        } else {
            return Err("viscoelastic wave integral phase is unrepresentable");
        };
        [det::ln(0.5 * j + 0.5 * l), det::ln(0.5 * j - 0.5 * l)]
    };
    if values
        .iter()
        .all(|v| v.is_finite() || *v == f64::NEG_INFINITY)
    {
        Ok(values)
    } else {
        Err("viscoelastic positive wave integral is unrepresentable")
    }
}

fn integral_log(parity: [C64; 2], log_gram: [f64; 2]) -> f64 {
    let terms = [
        log_gram[0] + 2.0 * det::ln(parity[0].abs()),
        log_gram[1] + 2.0 * det::ln(parity[1].abs()),
    ];
    let larger = terms[0].max(terms[1]);
    if larger == f64::NEG_INFINITY {
        return larger;
    }
    larger + det::ln(1.0 + det::exp(terms[0].min(terms[1]) - larger))
}

impl Slab {
    pub(crate) fn new(
        density: f64,
        moduli: [C64; 2],
        omega: f64,
        length: f64,
        exterior: [f64; 2],
    ) -> Result<Self> {
        if !density.is_finite()
            || density <= 0.0
            || !omega.is_finite()
            || omega <= 0.0
            || !length.is_finite()
            || length < 0.0
            || exterior.iter().any(|z| !z.is_finite() || *z <= 0.0)
            || moduli
                .iter()
                .any(|m| !finite(*m) || m.re <= 0.0 || m.im < 0.0)
        {
            return Err(
                "viscoelastic slab needs positive density, moduli, frequency and impedances",
            );
        }
        let modulus = moduli[0] + moduli[1].scale(4.0 / 3.0);
        let magnitude = modulus.abs();
        if !finite(modulus) || !magnitude.is_finite() {
            return Err("longitudinal viscoelastic modulus is unrepresentable");
        }
        let root = modulus.sqrt();
        let density_root = det::sqrt(density);
        let impedance = root.scale(density_root);
        let k = C64::new(
            product_ratio(omega, density_root, root.re, magnitude)?,
            -product_ratio(omega, density_root, root.im, magnitude)?,
        );
        if !finite(k) || k.re <= 0.0 || !finite(impedance) || impedance.re <= 0.0 {
            return Err("viscoelastic wavenumber or impedance is unrepresentable");
        }
        let g = propagation(k, length)?;
        let scale = impedance.abs().max(exterior[0]).max(exterior[1]);
        let z = divided(impedance, scale);
        let [l, r] = exterior.map(|x| x / scale);
        if l <= 0.0 || r <= 0.0 || z.re <= 0.0 {
            return Err("viscoelastic exterior contrast is unrepresentable");
        }
        let phase = k.re * length;
        let x = -k.im * length;
        let missing_round_trip = if g == C64::ZERO {
            C64::ONE
        } else {
            let (sine, cosine) = det::sin_cos(phase);
            let decay = det::exp(-2.0 * x);
            C64::new(
                -det::expm1(-2.0 * x) + 2.0 * decay * sine * sine,
                2.0 * decay * sine * cosine,
            )
        };
        let (zl, zr) = (C64::from_re(l), C64::from_re(r));
        // Factored determinant retains a thin lossy slab and high-Q round trips.
        let determinant = z.scale(2.0 * (l + r)) + (z - zl) * (z - zr) * missing_round_trip;
        if !finite(determinant) || determinant == C64::ZERO {
            return Err("viscoelastic boundary solve is singular or unrepresentable");
        }
        // Cancel the common impedance BEFORE multiplying a tiny propagation
        // phase. Otherwise Im(Z) can erase a representable Im(g), even when
        // the normalized transmission is simply g to working precision.
        let reduced = C64::from_re(2.0 * (l + r)) + ((z - zl) * (z - zr) * missing_round_trip) / z;
        let use_reduced = finite(reduced) && reduced != C64::ZERO;
        let factor = if use_reduced {
            C64::from_re(2.0) / reduced
        } else {
            z.scale(2.0) / determinant
        };
        let inward = [
            [factor * (z + zr), -(factor * (z - zl) * g)],
            [-(factor * (z - zr) * g), factor * (z + zl)],
        ];
        let scattering = if use_reduced {
            [
                (C64::from_re(2.0 * (r - l)) + ((z + zl) * (z - zr) * missing_round_trip) / z)
                    / reduced,
                (C64::from_re(4.0 * l) / reduced) * g,
                (C64::from_re(4.0 * r) / reduced) * g,
                (C64::from_re(2.0 * (l - r)) + ((z + zr) * (z - zl) * missing_round_trip) / z)
                    / reduced,
            ]
        } else {
            [
                (z.scale(2.0 * (r - l)) + (z + zl) * (z - zr) * missing_round_trip) / determinant,
                (z.scale(4.0 * l) / determinant) * g,
                (z.scale(4.0 * r) / determinant) * g,
                (z.scale(2.0 * (l - r)) + (z + zr) * (z - zl) * missing_round_trip) / determinant,
            ]
        };
        if inward
            .iter()
            .flatten()
            .chain(scattering.iter())
            .any(|z| !finite(*z))
        {
            return Err("viscoelastic boundary coefficients are unrepresentable");
        }
        Ok(Self {
            omega,
            length,
            modulus,
            moduli,
            k,
            impedance,
            exterior,
            scattering,
            inward,
            parity_factor: factor,
            normalized_impedance: z,
            normalized_exterior: [l, r],
            missing_one_way: missing_propagation(k, length)?,
            log_gram: wave_gram_logs(k, length)?,
        })
    }

    fn parity(&self, incoming: [C64; 2]) -> [C64; 2] {
        // Evaluate A+B and A-B before their large individual wave terms round.
        // In a thin standing wave, subtracting already-rounded A and B can
        // erase the entire even-mode contribution to positive local heating.
        let [a, b] = incoming;
        let [l, r] = self.normalized_exterior;
        let z = self.normalized_impedance;
        let sum = weighted(a, r, b, l);
        let difference = weighted(a, r, b, -l);
        let z_difference = z * (a - b);
        [
            self.parity_factor * (sum.scale(2.0) + (z * (a + b) - sum) * self.missing_one_way),
            self.parity_factor
                * (z_difference.scale(2.0) - (z_difference - difference) * self.missing_one_way),
        ]
    }

    fn heat(&self, log_pressure_square: f64, volume_factor: f64) -> Result<[f64; 2]> {
        let common = det::ln(self.omega) + log_pressure_square - 2.0 * det::ln(self.modulus.abs())
            + det::ln(volume_factor);
        Ok([
            if self.moduli[0].im == 0.0 {
                0.0
            } else {
                from_log(common + det::ln(self.moduli[0].im) - det::ln(2.0))?
            },
            if self.moduli[1].im == 0.0 {
                0.0
            } else {
                from_log(common + det::ln(self.moduli[1].im) + det::ln(2.0 / 3.0))?
            },
        ])
    }

    pub(crate) fn drive(&self, incoming: [C64; 2], area: f64) -> Result<Drive> {
        if incoming.iter().any(|z| !finite(*z) || !z.abs().is_finite())
            || !area.is_finite()
            || area <= 0.0
        {
            return Err("viscoelastic drive needs finite peak pressures and positive area");
        }
        let inward = self
            .inward
            .map(|row| row[0] * incoming[0] + row[1] * incoming[1]);
        let parity = self.parity(incoming);
        let outgoing = [
            self.scattering[0] * incoming[0] + self.scattering[1] * incoming[1],
            self.scattering[2] * incoming[0] + self.scattering[3] * incoming[1],
        ];
        if inward
            .iter()
            .chain(outgoing.iter())
            .chain(parity.iter())
            .any(|z| !finite(*z) || !z.abs().is_finite())
        {
            return Err("viscoelastic driven field is unrepresentable");
        }
        let power = |p: C64, z: f64| {
            from_log(2.0 * det::ln(p.abs()) + det::ln(area) - det::ln(2.0) - det::ln(z))
        };
        let incident_power = [
            power(incoming[0], self.exterior[0])?,
            power(incoming[1], self.exterior[1])?,
        ];
        let outgoing_power = [
            power(outgoing[0], self.exterior[0])?,
            power(outgoing[1], self.exterior[1])?,
        ];
        let heat = self.heat(integral_log(parity, self.log_gram), area)?;
        let positive_sum = incident_power
            .into_iter()
            .chain(outgoing_power)
            .chain(heat)
            .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value))
            .to_f64();
        let residual = incident_power
            .into_iter()
            .fold(Dd::ZERO, |sum, x| sum + Dd::from_f64(x))
            - outgoing_power
                .into_iter()
                .chain(heat)
                .fold(Dd::ZERO, |sum, x| sum + Dd::from_f64(x));
        let residual = residual.to_f64();
        if !positive_sum.is_finite()
            || !residual.is_finite()
            || residual.abs() > 4096.0 * f64::EPSILON * positive_sum.max(f64::MIN_POSITIVE)
        {
            return Err("viscoelastic independent boundary and local heat accounts disagree");
        }
        let a = inward[0].abs();
        let b = inward[1].abs();
        let scale = a.max(b);
        let strain_bound = if scale == 0.0 || self.length == 0.0 {
            0.0
        } else {
            from_log(det::ln(scale) + det::ln(a / scale + b / scale) - det::ln(self.modulus.abs()))?
        };
        Ok(Drive {
            inward,
            parity,
            outgoing,
            incident_power,
            outgoing_power,
            heat,
            residual,
            strain_bound,
        })
    }

    /// Positive original strain-work integral on a physical depth interval.
    /// Restriction changes wave coordinates only; no new interface is inserted.
    pub(crate) fn interval_heat(
        &self,
        inward: [C64; 2],
        parity: [C64; 2],
        depths: [f64; 2],
        area: f64,
    ) -> Result<[f64; 2]> {
        let [start, end] = depths;
        if !start.is_finite()
            || !end.is_finite()
            || start < 0.0
            || end < start
            || end > self.length
            || !area.is_finite()
            || area <= 0.0
        {
            return Err(
                "viscoelastic heat interval needs ordered admitted depths and positive area",
            );
        }
        if start == end {
            return Ok([0.0; 2]);
        }
        if start == 0.0 && end == self.length {
            return self.heat(integral_log(parity, self.log_gram), area);
        }
        let back = self.length - end;
        let difference = (start - back).abs();
        let (local, log_scale) = if self.k.abs() * difference < 0.1 {
            // Keep the ORIGINAL accurately prepared A+B and A-B when the
            // restriction has almost equal propagation on its two sides.
            // A common phase cancels from |p|^2. Keep its attenuation in log
            // form so a small field times a large area can remain representable.
            let missing = missing_propagation(self.k, difference)?;
            let even = C64::ONE - missing.scale(0.5);
            let odd = missing.scale(if start <= back { 0.5 } else { -0.5 });
            (
                [
                    parity[0] * even + parity[1] * odd,
                    parity[1] * even + parity[0] * odd,
                ],
                self.k.im * start.min(back),
            )
        } else {
            // A far opaque restriction must NOT recover a tiny forward wave
            // by subtracting (A+B) and (A-B). Scale each bounded inward wave
            // by its own logarithmic amplitude first, then combine them.
            let magnitudes = inward.map(C64::abs);
            let logs = [
                det::ln(magnitudes[0]) + self.k.im * start,
                det::ln(magnitudes[1]) + self.k.im * back,
            ];
            let dominant = usize::from(logs[1] > logs[0]);
            let scale = logs[dominant];
            if scale == f64::NEG_INFINITY {
                return Ok([0.0; 2]);
            }
            let distances = [start, back];
            let mut waves = [C64::ZERO; 2];
            for i in 0..2 {
                if magnitudes[i] == 0.0 {
                    continue;
                }
                let relative = det::exp(logs[i] - scale);
                if relative == 0.0 {
                    continue;
                }
                let phase = self.k.re * (distances[i] - distances[dominant]);
                if !phase.is_finite() {
                    return Err("viscoelastic interval phase is unrepresentable");
                }
                let (sine, cosine) = det::sin_cos(phase);
                waves[i] =
                    divided(inward[i], magnitudes[i]) * C64::new(cosine, -sine).scale(relative);
            }
            ([waves[0] + waves[1], waves[0] - waves[1]], scale)
        };
        if local.iter().any(|v| !finite(*v) || !v.abs().is_finite()) {
            return Err("viscoelastic restricted wave is unrepresentable");
        }
        self.heat(
            integral_log(local, wave_gram_logs(self.k, end - start)?) + 2.0 * log_scale,
            area,
        )
    }

    pub(crate) fn field(&self, inward: [C64; 2], parity: [C64; 2], depth: f64) -> Result<Field> {
        if !depth.is_finite() || depth < 0.0 || depth > self.length || self.length == 0.0 {
            return Err(
                "viscoelastic interior observation needs a positive slab and admitted depth",
            );
        }
        let (pressure, velocity) = if self.k.abs() * self.length < 0.1 {
            let back_depth = self.length - depth;
            let forward_nearer = depth <= back_depth;
            let near_depth = depth.min(back_depth);
            let distance = depth.max(back_depth) - near_depth;
            let g = propagation(self.k, near_depth)?;
            let missing = missing_propagation(self.k, distance)?;
            let even = g * (C64::ONE - missing.scale(0.5));
            let odd = (g * missing).scale(if forward_nearer { 0.5 } else { -0.5 });
            (
                parity[0] * even + parity[1] * odd,
                (parity[1] * even + parity[0] * odd) / self.impedance,
            )
        } else {
            let a = inward[0] * propagation(self.k, depth)?;
            let b = inward[1] * propagation(self.k, self.length - depth)?;
            (a + b, (a - b) / self.impedance)
        };
        self.field_from_trace(pressure, velocity)
    }

    /// Reuse the original constitutive observation at an independently retained
    /// boundary trace, without recovering tiny velocity from opposite waves.
    pub(crate) fn field_from_trace(&self, pressure: C64, velocity: C64) -> Result<Field> {
        let strain = -(pressure / self.modulus);
        if [pressure, velocity, strain].iter().any(|z| !finite(*z)) {
            return Err("viscoelastic interior field is unrepresentable");
        }
        Ok(Field {
            pressure,
            velocity,
            strain,
            heat_density: self.heat(2.0 * det::ln(pressure.abs()), 1.0)?,
        })
    }
}

#[cfg(test)]
#[path = "kernel/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "kernel/interval_tests.rs"]
mod interval_tests;
