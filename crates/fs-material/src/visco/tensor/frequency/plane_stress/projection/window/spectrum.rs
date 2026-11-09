//! Coherent real-field products with the original binary material scales.

use super::*;

pub(crate) const MAX_SPECTRUM_LINES: usize = 64;
pub(crate) type StrainParts = [[Dd; 2]; 6];
type Spectrum = [(f64, i32); 2];
type Channel = [Sum; 2];

/// Shared temporal cache. Original material factors are supplied per point.
pub(crate) struct HarmonicEnergySpectrum {
    duration: f64,
    frequencies: Vec<f64>,
    start: Vec<[Dd; 2]>,
    center: Vec<[Dd; 2]>,
    end: Vec<[Dd; 2]>,
    gram: Vec<[Dd; 2]>,
    means: Vec<Dd>,
}

fn project_parts(row: &[f64; 6], strain: &StrainParts, spectrum: Spectrum) -> Channel {
    let [a, b] = spectrum;
    let (mut re, mut im) = (Sum::default(), Sum::default());
    for (&map, value) in row.iter().zip(strain) {
        for [r, i] in [[value[0].hi, value[1].hi], [value[0].lo, value[1].lo]] {
            re.add_binary_ratio(&[map, r, a.0], &[], a.1);
            re.add_binary_ratio(&[map, i, b.0], &[], b.1);
            im.add_binary_ratio(&[map, i, a.0], &[], a.1);
            im.add_binary_ratio(&[-map, r, b.0], &[], b.1);
        }
    }
    [re, im]
}

impl HarmonicEnergySpectrum {
    pub(crate) fn new(frequencies: &[f64], window: HarmonicWindow) -> Result<Self, WindowError> {
        Self::at_clock(frequencies, window, window.start_clock_dd())
    }

    pub(crate) fn new_at_clock(
        frequencies: &[f64], start_s: Dd, duration_s: f64,
    ) -> Result<Self, WindowError> {
        if !start_s.hi.is_finite() || !start_s.lo.is_finite() {
            return Err(WindowError::InvalidWindow);
        }
        let start_s = Dd::from_pair(start_s.hi, start_s.lo);
        Self::at_clock(frequencies, HarmonicWindow::from_clock_dd(start_s, duration_s)?, start_s)
    }

    fn at_clock(frequencies: &[f64], window: HarmonicWindow, start_s: Dd) -> Result<Self, WindowError> {
        if frequencies.len() > MAX_SPECTRUM_LINES {
            return Err(WindowError::Unrepresentable);
        }
        let mut start = Vec::with_capacity(frequencies.len());
        let mut center = Vec::with_capacity(frequencies.len());
        let mut end = Vec::with_capacity(frequencies.len());
        let mut means = Vec::with_capacity(frequencies.len());
        for &omega in frequencies {
            // Same separated clock/width rotations as the single-tone owner.
            let mut temporal = HarmonicEnergyWindow::new(omega, window)?;
            if start_s.lo != 0.0 {
                // The accepted state clock keeps a short previous interval in
                // its low word. Reuse the original separated width rotations.
                let w = Dd::from_f64(omega);
                let half_width = if omega.abs() >= window.duration_s() {
                    (w * Dd::from_f64(0.5)) * Dd::from_f64(window.duration_s())
                } else {
                    w * (Dd::from_f64(window.duration_s()) * Dd::from_f64(0.5))
                };
                let half = phase(half_width)?;
                temporal.start = phase(w * start_s)?;
                temporal.center = phase_product(temporal.start, half);
                temporal.end = phase_product(temporal.center, half);
            }
            start.push(temporal.start);
            center.push(temporal.center);
            end.push(temporal.end);
            // Integral of the centered carrier cosine. The sine mean is zero.
            // The original zero-frequency product retains arbitrarily short
            // windows without subtracting two rounded endpoint phases.
            means.push(window.product_kernel(omega, 0.0)?.centered_gram_dd()[0]);
        }
        let mut gram = Vec::with_capacity(frequencies.len() * frequencies.len());
        for &a in frequencies {
            for &b in frequencies {
                gram.push(window.product_kernel(a, b)?.centered_gram_dd());
            }
        }
        Ok(Self {
            duration: window.duration_s(),
            frequencies: frequencies.to_vec(),
            start,
            center,
            end,
            gram,
            means,
        })
    }

    pub(crate) fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>()
            + self.frequencies.len()
                * (core::mem::size_of::<f64>()
                    + 3 * core::mem::size_of::<[Dd; 2]>()
                    + core::mem::size_of::<Dd>())
            + self.gram.len() * core::mem::size_of::<[Dd; 2]>()
    }

    fn centered(
        &self,
        row: &[f64; 6],
        strains: &[StrainParts],
        spectrum: &[Spectrum],
    ) -> [Channel; MAX_SPECTRUM_LINES] {
        let mut values = [[Sum::default(); 2]; MAX_SPECTRUM_LINES];
        for i in 0..strains.len() {
            values[i] = rotate(project_parts(row, &strains[i], spectrum[i]), self.center[i]);
        }
        values
    }

    fn product_integral(&self, a: &[Channel], b: &[Channel], denominator: &[f64], total: &mut Sum) {
        let n = self.frequencies.len();
        for i in 0..n {
            for j in 0..n {
                for k in 0..2 {
                    weighted(
                        total,
                        a[i][k],
                        b[j][k],
                        self.gram[i * n + j][k],
                        &[self.duration],
                        denominator,
                    );
                }
            }
        }
    }

    /// Original spring energy, direct storage-rate integral, and physical heat.
    pub(crate) fn spring(
        &self,
        gram: &[[f64; 6]; 6],
        strains: &[StrainParts],
        spectrum: &[Spectrum],
        tau: Option<f64>,
    ) -> Result<[f64; 4], WindowError> {
        self.spring_with_bias(gram, strains, spectrum, tau, &[Dd::ZERO; 6])
    }

    // A relaxed constant strain belongs only to the equilibrium spring. The
    // caller supplies zero bias for every Maxwell branch; no DC memory is made.
    pub(crate) fn spring_with_bias(
        &self,
        gram: &[[f64; 6]; 6],
        strains: &[StrainParts],
        spectrum: &[Spectrum],
        tau: Option<f64>,
        bias: &[Dd; 6],
    ) -> Result<[f64; 4], WindowError> {
        let n = self.frequencies.len();
        let mut totals = [Sum::default(); 4];
        for row in gram {
            let mut constant = Sum::default();
            for (&map, value) in row.iter().zip(bias) {
                constant.add([map, value.hi, 1.0, 1.0]);
                constant.add([map, value.lo, 1.0, 1.0]);
            }
            let mut endpoints = [constant; 2];
            for i in 0..n {
                let z = project_parts(row, &strains[i], spectrum[i]);
                add(&mut endpoints[0], rotate(z, self.start[i])[0]);
                add(&mut endpoints[1], rotate(z, self.end[i])[0]);
            }
            for i in 0..2 {
                weighted(
                    &mut totals[i],
                    endpoints[i],
                    endpoints[i],
                    Dd::ONE,
                    &[0.5],
                    &[],
                );
            }
            let values = self.centered(row, strains, spectrum);
            for i in 0..n {
                if constant.fraction != Dd::ZERO {
                    weighted(
                        &mut totals[2],
                        constant,
                        values[i][1],
                        self.means[i],
                        &[self.frequencies[i], self.duration],
                        &[],
                    );
                }
                let g = self.gram[i * n + i];
                weighted(
                    &mut totals[2],
                    values[i][0],
                    values[i][1],
                    g[0] - g[1],
                    &[self.frequencies[i], self.duration],
                    &[],
                );
                for j in i + 1..n {
                    let [cc, ss] = self.gram[i * n + j];
                    // Differentiate each cross product as a sum- and a beat-
                    // frequency contribution. Subtract frequencies before a
                    // large carrier multiplication can erase the slow beat.
                    let mut sum_frequency = Sum::default();
                    sum_frequency.add([self.frequencies[i], 1.0, 1.0, 1.0]);
                    sum_frequency.add([self.frequencies[j], 1.0, 1.0, 1.0]);
                    let difference =
                        dd(Dd::from_f64(self.frequencies[j]) - Dd::from_f64(self.frequencies[i]));
                    let xy = product(values[i][0], values[j][1]);
                    let yx = product(values[i][1], values[j][0]);
                    let mut symmetric = Sum::default();
                    add(&mut symmetric, xy);
                    add(&mut symmetric, yx);
                    let mut antisymmetric = Sum::default();
                    add(&mut antisymmetric, xy);
                    add(&mut antisymmetric, scaled(yx, -Dd::ONE));
                    weighted(
                        &mut totals[2],
                        symmetric,
                        sum_frequency,
                        cc - ss,
                        &[0.5, self.duration],
                        &[],
                    );
                    weighted(
                        &mut totals[2],
                        antisymmetric,
                        difference,
                        cc + ss,
                        &[0.5, self.duration],
                        &[],
                    );
                }
            }
            if let Some(tau) = tau {
                self.product_integral(&values[..n], &values[..n], &[tau], &mut totals[3]);
            }
        }
        Ok([
            positive(totals[0].value())?,
            positive(totals[1].value())?,
            finite(totals[2].value())?,
            positive(totals[3].value())?,
        ])
    }

    pub(crate) fn kelvin_heat(
        &self,
        gram: &[[f64; 6]; 6],
        strains: &[StrainParts],
    ) -> Result<f64, WindowError> {
        let n = self.frequencies.len();
        let mut rates = [[(0.0, 0); 2]; MAX_SPECTRUM_LINES];
        for (rate, &omega) in rates.iter_mut().zip(&self.frequencies) {
            *rate = [(0.0, 0), (omega, 0)];
        }
        let mut total = Sum::default();
        for row in gram {
            let values = self.centered(row, strains, &rates[..n]);
            self.product_integral(&values[..n], &values[..n], &[], &mut total);
        }
        positive(total.value())
    }

    /// Original stress : strain_rate, independent of storage and loss factors.
    pub(crate) fn add_stress_work(
        &self,
        matrix: &[[f64; 6]; 6],
        strains: &[StrainParts],
        spectrum: &[Spectrum],
        total: &mut Sum,
    ) {
        self.add_stress_work_with_bias(matrix, strains, spectrum, &[Dd::ZERO; 6], total);
    }

    /// The equilibrium DC stress works on the actual AC strain increment.
    /// Its original constitutive tensor remains independent of storage factors.
    pub(crate) fn add_stress_work_with_bias(
        &self,
        matrix: &[[f64; 6]; 6],
        strains: &[StrainParts],
        spectrum: &[Spectrum],
        bias: &[Dd; 6],
        total: &mut Sum,
    ) {
        let n = self.frequencies.len();
        for (i, row) in matrix.iter().enumerate() {
            let mut constant = Sum::default();
            for (&map, value) in row.iter().zip(bias) {
                constant.add([map, value.hi, 1.0, 1.0]);
                constant.add([map, value.lo, 1.0, 1.0]);
            }
            let stress = self.centered(row, strains, spectrum);
            let mut rate = [[Sum::default(); 2]; MAX_SPECTRUM_LINES];
            let mut unit = [0.0; 6];
            unit[i] = 1.0;
            for j in 0..n {
                rate[j] = rotate(
                    project_parts(&unit, &strains[j], [(0.0, 0), (self.frequencies[j], 0)]),
                    self.center[j],
                );
                if constant.fraction != Dd::ZERO {
                    weighted(
                        total,
                        constant,
                        rate[j][0],
                        self.means[j],
                        &[self.duration],
                        &[],
                    );
                }
            }
            self.product_integral(&stress[..n], &rate[..n], &[], total);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() < 3e-11 * scale,
            "{actual:e} != {expected:e}, scale={scale:e}"
        );
    }

    #[test]
    fn independent_real_maxwell_odes_close_multitone_storage_heat_and_work() {
        let frequencies = [2.7, 3.1, 11.0];
        let tau = 0.3;
        let strains: [StrainParts; 3] = core::array::from_fn(|f| {
            core::array::from_fn(|i| {
                [
                    Dd::from_f64((f + 1) as f64 * (i as f64 - 2.0) * 0.001),
                    Dd::from_f64((2.0 - f as f64) * (3.0 - i as f64) * 0.0007),
                ]
            })
        });
        let spectrum: [Spectrum; 3] = frequencies.map(|w| {
            let x = w * tau;
            [(x * x / (1.0 + x * x), 0), (x / (1.0 + x * x), 0)]
        });
        let mut gram = [[0.0; 6]; 6];
        gram[0] = [1.0, -2.0, 0.3, 0.0, 0.5, -0.1];
        gram[1] = [0.2, 0.0, -0.6, 0.7, 0.1, -0.2];
        let matrix = core::array::from_fn(|i| {
            core::array::from_fn(|j| (0..6).map(|k| gram[k][i] * gram[k][j]).sum())
        });
        let window = HarmonicWindow::new(0.37, 1.9).unwrap();
        let temporal = HarmonicEnergySpectrum::new(&frequencies, window).unwrap();
        let response = temporal
            .spring(&gram, &strains, &spectrum, Some(tau))
            .unwrap();
        let mut work = Sum::default();
        temporal.add_stress_work(&matrix, &strains, &spectrum, &mut work);
        // Independently solve tau*q_dot+q=tau*epsilon_dot in real fields.
        let instant = |t: f64| {
            let (mut q, mut rate) = ([0.0; 6], [0.0; 6]);
            for f in 0..3 {
                let w = frequencies[f];
                let x = w * tau;
                let (s, c) = (w * t).sin_cos();
                for i in 0..6 {
                    let [re, im] = strains[f][i].map(Dd::to_f64);
                    let a = (x * im + x * x * re) / (1.0 + x * x);
                    let b = (x * x * im - x * re) / (1.0 + x * x);
                    q[i] += a * c + b * s;
                    rate[i] += w * (im * c - re * s);
                }
            }
            let stress: [f64; 6] =
                core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * q[j]).sum());
            let square: f64 = q.iter().zip(stress).map(|(a, b)| a * b).sum();
            [
                0.5 * square,
                square / tau,
                rate.iter().zip(stress).map(|(a, b)| a * b).sum(),
            ]
        };
        let (mut heat, mut direct_work) = (0.0, 0.0);
        let panels = 8192;
        for i in 0..=panels {
            let t = window.start_s() + window.duration_s() * i as f64 / panels as f64;
            let [_, h, w] = instant(t);
            let weight = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            heat += weight * h;
            direct_work += weight * w;
        }
        heat *= window.duration_s() / (3.0 * panels as f64);
        direct_work *= window.duration_s() / (3.0 * panels as f64);
        let start = instant(window.start_s())[0];
        let end = instant(window.start_s() + window.duration_s())[0];
        let scale = heat.max(start).max(end);
        near(response[0], start, scale);
        near(response[1], end, scale);
        near(response[2], end - start, scale);
        near(response[3], heat, scale);
        near(work.value(), direct_work, scale);
        near(work.value() - response[2] - response[3], 0.0, scale);
    }

    #[test]
    fn one_line_retains_tiny_window_increment_and_binary_constitutive_scale() {
        let omega = 3.0;
        let e = [C64::new(1e-150, -0.4e-150); 6];
        let strains: [StrainParts; 1] = [e.map(|z| [Dd::from_f64(z.re), Dd::from_f64(z.im)])];
        let gram =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1e150 } else { 0.0 }));
        let transfer = [(0.5, -1100), (0.7, -550)];
        let tau = 1e-300;
        let window = HarmonicWindow::new(0.4, 1e-12).unwrap();
        let single = HarmonicEnergyWindow::new(omega, window)
            .unwrap()
            .spring(&gram, &e, transfer, Some(tau))
            .unwrap();
        let multiple = HarmonicEnergySpectrum::new(&[omega], window)
            .unwrap()
            .spring(&gram, &strains, &[transfer], Some(tau))
            .unwrap();
        assert!(single[3] > 0.0);
        for (a, b) in multiple.into_iter().zip(single) {
            if b == 0.0 {
                assert_eq!(a, b);
            } else {
                near(a, b, b.abs());
            }
        }
    }

    #[test]
    fn coherent_close_tone_kelvin_heat_retains_physical_beat() {
        let frequencies = [10.0, 10.2];
        let mut e = [[Dd::ZERO; 2]; 6];
        e[0][0] = Dd::from_f64(0.01);
        let gram =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 2.0 } else { 0.0 }));
        let window = HarmonicWindow::new(0.0, 0.7).unwrap();
        let actual = HarmonicEnergySpectrum::new(&frequencies, window)
            .unwrap()
            .kelvin_heat(&gram, &[e; 2])
            .unwrap();
        let diagonal: f64 = frequencies
            .iter()
            .map(|&w| {
                HarmonicEnergySpectrum::new(&[w], window)
                    .unwrap()
                    .kelvin_heat(&gram, &[e])
                    .unwrap()
            })
            .sum();
        assert!(actual > 1.9 * diagonal);
        let mut integral = 0.0;
        let panels = 8192;
        for i in 0..=panels {
            let t = window.duration_s() * i as f64 / panels as f64;
            let rate: f64 = frequencies.iter().map(|&w| -w * 0.01 * (w * t).sin()).sum();
            let weight = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            integral += weight * 4.0 * rate * rate;
        }
        integral *= window.duration_s() / (3.0 * panels as f64);
        near(actual, integral, integral);
    }
}
