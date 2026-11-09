//! Coherent multitone material heat and work, retaining the original memories.

use super::super::plane_stress::projection::{
    Sum,
    window::{HarmonicEnergySpectrum, MAX_SPECTRUM_LINES, StrainParts},
};
use super::{TensorProny, TensorPronyHarmonicWindowPoint, ViscoError, invalid, positive, scalar};
use crate::visco::prony_spectrum::PronySpectrumWeights;
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_math::{c64::C64, dd::Dd};

/// Bounded authored positive-frequency lines, including equal-frequency entries.
pub const MAX_HARMONIC_MATERIAL_LINES: usize = MAX_SPECTRUM_LINES;

/// Material-admitted settled strains under one original exp(-i omega t) clock.
///
/// Each original frequency has its own full anisotropic Maxwell response.
/// Frequencies need not be commensurate. Exact equal-frequency strains combine
/// before physical admission or energy products, preserving their interference.
/// This is neither a fitted damping law nor a causal startup history.
pub struct TensorPronyHarmonicSpectrum<'a> {
    law: &'a TensorProny,
    authored_frequencies: Vec<f64>,
    frequencies: Vec<f64>,
    groups: Vec<usize>,
    weights: Vec<Vec<PronySpectrumWeights>>,
}

/// Cached exact temporal pair moments, reusable across original material points.
/// Preparation takes O(F²) work and storage for F distinct frequencies. A point
/// uses O(F²*(branches+1)*6) contractions with fixed stack scratch and no heap
/// allocation. Source geometry, frequency/strain domain and coefficients stay fixed.
pub struct TensorPronySpectrumWindow<'p, 'a> {
    source: &'p TensorPronyHarmonicSpectrum<'a>,
    window: HarmonicWindow,
    temporal: HarmonicEnergySpectrum,
}

impl TensorProny {
    /// Prepare at most 64 positive finite source-admitted angular frequencies.
    /// The input order is retained for all point strain arrays. Empty input is
    /// the zero AC field. DC/prestrain and startup histories are not inferred.
    ///
    /// Retained heap payload is at most size_of::<TensorPronyHarmonicSpectrum>()
    /// + 128*F*(branches+1) bytes, excluding source law and allocator metadata.
    pub fn prepare_harmonic_spectrum(
        &self,
        angular_frequencies_rad_s: &[f64],
    ) -> Result<TensorPronyHarmonicSpectrum<'_>, ViscoError> {
        if angular_frequencies_rad_s.len() > MAX_HARMONIC_MATERIAL_LINES {
            return Err(ViscoError::Parameters {
                what: "coherent material observation exceeds its 64-line limit",
            });
        }
        let mut frequencies = Vec::with_capacity(angular_frequencies_rad_s.len());
        let mut groups = Vec::with_capacity(angular_frequencies_rad_s.len());
        for &omega in angular_frequencies_rad_s {
            if !omega.is_finite() || omega <= 0.0 {
                return Err(ViscoError::Parameters {
                    what: "material spectrum requires positive finite angular frequencies",
                });
            }
            super::super::admit_frequency(self, omega)?;
            let group = frequencies
                .iter()
                .position(|&w| w == omega)
                .unwrap_or_else(|| {
                    frequencies.push(omega);
                    frequencies.len() - 1
                });
            groups.push(group);
        }
        let weights = frequencies
            .iter()
            .map(|&omega| {
                self.branches
                    .iter()
                    .map(|b| PronySpectrumWeights::new(omega, b.tau_s).ok_or_else(invalid))
                    .collect()
            })
            .collect::<Result<Vec<Vec<_>>, _>>()?;
        Ok(TensorPronyHarmonicSpectrum {
            law: self,
            authored_frequencies: angular_frequencies_rad_s.to_vec(),
            frequencies,
            groups,
            weights,
        })
    }
}

impl<'a> TensorPronyHarmonicSpectrum<'a> {
    pub const fn law(&self) -> &'a TensorProny {
        self.law
    }
    /// Original caller order, including repeated frequencies [rad/s].
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.authored_frequencies
    }
    /// Retained numeric payload, excluding the separately owned source law.
    pub fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>()
            + self.authored_frequencies.len() * core::mem::size_of::<f64>()
            + self.groups.len() * core::mem::size_of::<usize>()
            + self.frequencies.capacity() * core::mem::size_of::<f64>()
            + self.weights.len() * core::mem::size_of::<Vec<PronySpectrumWeights>>()
            + self
                .weights
                .iter()
                .map(|w| w.len() * core::mem::size_of::<PronySpectrumWeights>())
                .sum::<usize>()
    }

    /// Admit the complete coherent strain waveform, not just its sampled window.
    /// Equal-frequency lines first combine with retained compensated arithmetic.
    /// For distinct frequencies the sum of each ellipse's exact peak norm is a
    /// conservative all-time upper bound; no common period or phase sampling is
    /// assumed. This bound can refuse a valid waveform whose separate peaks can
    /// never coincide. It is not a certified supremum or spatial error bound.
    pub fn validate_strains(&self, strains: &[[C64; 6]]) -> Result<(), ViscoError> {
        let grouped = self.group_strains(strains)?;
        self.admit_grouped_strains(&grouped[..self.frequencies.len()])
    }

    /// Admit settled AC strains superposed on an explicitly relaxed DC strain.
    /// A nonzero bias also requires zero frequency in the original source band.
    /// The all-time bound is the DC norm plus the distinct AC ellipse bounds;
    /// initial Maxwell memory or geometric prestress effects are not inferred.
    pub fn validate_strains_with_bias(
        &self,
        strains: &[[C64; 6]],
        bias: &[f64; 6],
    ) -> Result<(), ViscoError> {
        let grouped = self.group_strains(strains)?;
        self.admit_grouped_strains_with_bias(
            &grouped[..self.frequencies.len()],
            &bias.map(Dd::from_f64),
        )
    }

    pub(in crate::visco::tensor::frequency) fn distinct_frequencies(&self) -> &[f64] {
        &self.frequencies
    }

    // Equal-frequency authored coordinates must combine before a physical
    // constraint is recovered. In particular, separately admitting thin-body
    // strains would reject a cancelling field or overflow its transverse map.
    pub(in crate::visco::tensor::frequency) fn group_strains<const N: usize>(
        &self,
        strains: &[[C64; N]],
    ) -> Result<[[[Dd; 2]; N]; MAX_SPECTRUM_LINES], ViscoError> {
        if strains.len() != self.authored_frequencies.len()
            || strains
                .iter()
                .flatten()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(ViscoError::Parameters {
                what: "one finite strain phasor per original spectral line required",
            });
        }
        let mut sums = [[[Sum::default(); 2]; N]; MAX_SPECTRUM_LINES];
        for (strain, &group) in strains.iter().zip(&self.groups) {
            for i in 0..N {
                sums[group][i][0].add([strain[i].re, 1.0, 1.0, 1.0]);
                sums[group][i][1].add([strain[i].im, 1.0, 1.0, 1.0]);
            }
        }
        Ok(sums.map(|strain| strain.map(|z| z.map(Sum::retained_value))))
    }

    pub(in crate::visco::tensor::frequency) fn admit_grouped_strains(
        &self,
        strains: &[StrainParts],
    ) -> Result<(), ViscoError> {
        self.admit_grouped_strains_with_bias(strains, &[Dd::ZERO; 6])
    }

    pub(in crate::visco::tensor::frequency) fn admit_grouped_strains_with_bias(
        &self,
        strains: &[StrainParts],
        bias: &[Dd; 6],
    ) -> Result<(), ViscoError> {
        if strains.len() != self.frequencies.len() {
            return Err(invalid());
        }
        if bias
            .iter()
            .any(|value| !value.hi.is_finite() || !value.lo.is_finite())
        {
            return Err(ViscoError::Parameters {
                what: "relaxed material bias must be a finite Mandel strain",
            });
        }
        if bias.iter().any(|&value| value != Dd::ZERO) {
            super::super::admit_frequency(self.law, 0.0)?;
        }
        let mut envelope = Dd::ZERO;
        for word in 0..2 {
            let phasor =
                bias.map(|value| C64::from_re(if word == 0 { value.hi } else { value.lo }));
            envelope = envelope
                + Dd::from_f64(super::super::strain::peak_norm(&phasor).ok_or_else(invalid)?);
        }
        for strain in strains {
            // Retain the low words in subsequent material contractions. Their
            // separate norm contributes conservatively to source admission.
            for word in 0..2 {
                let phasor: [C64; 6] = core::array::from_fn(|i| {
                    let z = strain[i];
                    if word == 0 {
                        C64::new(z[0].hi, z[1].hi)
                    } else {
                        C64::new(z[0].lo, z[1].lo)
                    }
                });
                let peak = super::super::strain::peak_norm(&phasor).ok_or_else(invalid)?;
                envelope = envelope + Dd::from_f64(peak);
            }
        }
        let excess = if self.law.strain_limit.is_finite() {
            envelope - Dd::from_f64(self.law.strain_limit)
        } else {
            Dd::ZERO
        };
        if !envelope.hi.is_finite()
            || !envelope.lo.is_finite()
            || (self.law.strain_limit.is_finite()
                && (excess.hi > 0.0 || (excess.hi == 0.0 && excess.lo > 0.0)))
        {
            return Err(ViscoError::Parameters {
                what: "coherent waveform strain bound exceeds the original tensor material domain",
            });
        }
        Ok(())
    }

    /// Cache unweighted physical interval products once, for every source point.
    /// Additional temporal payload is bounded by size_of::<TensorPronySpectrumWindow>()
    /// + 256*F + 64*F² bytes. Source spectrum payload is retained separately.
    pub fn prepare_window(
        &self,
        window: HarmonicWindow,
    ) -> Result<TensorPronySpectrumWindow<'_, 'a>, ViscoError> {
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(ViscoError::Parameters {
                what: "physical material energy requires a uniform observation window",
            });
        }
        Ok(TensorPronySpectrumWindow {
            source: self,
            window,
            temporal: HarmonicEnergySpectrum::new(&self.frequencies, window)
                .map_err(|_| invalid())?,
        })
    }
}

impl<'p, 'a> TensorPronySpectrumWindow<'p, 'a> {
    pub const fn source(&self) -> &'p TensorPronyHarmonicSpectrum<'a> {
        self.source
    }
    pub const fn window(&self) -> HarmonicWindow {
        self.window
    }
    /// Additional temporal payload; the separately borrowed source is excluded.
    pub fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>() + self.temporal.numerical_payload_bytes()
            - core::mem::size_of::<HarmonicEnergySpectrum>()
    }

    /// Exact real-field pair moments retain carrier and beat interference.
    /// Heat comes from original positive Gram channels; work comes independently
    /// from stress : strain_rate. Storage change integrates the original spring
    /// derivatives instead of subtracting two rounded endpoint energies.
    pub fn evaluate(
        &self,
        strains: &[[C64; 6]],
    ) -> Result<TensorPronyHarmonicWindowPoint, ViscoError> {
        self.evaluate_with_checkpoint(strains, || Ok(()))
    }

    /// Same observation with a caller-owned cancellation/error channel. Polls
    /// before strain admission and every Maxwell branch, Kelvin and work group.
    /// Each bounded group uses at most 12*64² real pair contractions. No accepted
    /// history or partially returned field is mutated on refusal.
    pub fn evaluate_with_checkpoint<E: From<ViscoError>>(
        &self,
        strains: &[[C64; 6]],
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyHarmonicWindowPoint, E> {
        checkpoint()?;
        let grouped = self.source.group_strains(strains)?;
        self.evaluate_grouped_with_checkpoint(&grouped[..self.source.frequencies.len()], checkpoint)
    }

    /// Observe a relaxed static strain plus the original coherent AC waveform.
    ///
    /// Only the equilibrium spring carries the DC bias. Its cross terms enter
    /// actual endpoint storage, directly integrated storage change and the
    /// independent stress-work integral. Maxwell and Kelvin heat remain their
    /// original AC values. Nonzero bias requires original source admission at
    /// zero frequency; exactly zero preserves the existing AC-only contract.
    /// This is a settled constitutive observation, not a startup/state advance,
    /// finite-strain or geometric-stiffness/prestress calculation.
    pub fn evaluate_with_bias(
        &self,
        strains: &[[C64; 6]],
        bias: &[f64; 6],
    ) -> Result<TensorPronyHarmonicWindowPoint, ViscoError> {
        self.evaluate_with_bias_and_checkpoint(strains, bias, || Ok(()))
    }

    /// Same biased observation with the original bounded cancellation points.
    pub fn evaluate_with_bias_and_checkpoint<E: From<ViscoError>>(
        &self,
        strains: &[[C64; 6]],
        bias: &[f64; 6],
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyHarmonicWindowPoint, E> {
        checkpoint()?;
        let grouped = self.source.group_strains(strains)?;
        self.evaluate_grouped_with_bias_and_checkpoint(
            &grouped[..self.source.frequencies.len()],
            &bias.map(Dd::from_f64),
            checkpoint,
        )
    }

    // Shared by source-native constrained materials. The retained low words
    // survive their transverse recovery and every original energy contraction.
    pub(in crate::visco::tensor::frequency) fn evaluate_grouped_with_checkpoint<
        E: From<ViscoError>,
    >(
        &self,
        strains: &[StrainParts],
        checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyHarmonicWindowPoint, E> {
        self.evaluate_grouped_with_bias_and_checkpoint(strains, &[Dd::ZERO; 6], checkpoint)
    }

    pub(in crate::visco::tensor::frequency) fn evaluate_grouped_with_bias_and_checkpoint<
        E: From<ViscoError>,
    >(
        &self,
        strains: &[StrainParts],
        bias: &[Dd; 6],
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyHarmonicWindowPoint, E> {
        checkpoint()?;
        self.source.admit_grouped_strains_with_bias(strains, bias)?;
        let n = strains.len();
        let law = self.source.law;
        let mut spectrum = [[(1.0, 0), (0.0, 0)]; MAX_SPECTRUM_LINES];
        let equilibrium = self
            .temporal
            .spring_with_bias(&law.equilibrium_gram, strains, &spectrum[..n], None, bias)
            .map_err(|_| invalid())?;
        let mut storage = [equilibrium[0], equilibrium[1], equilibrium[2]].map(Dd::from_f64);
        let mut heat = Dd::ZERO;
        let mut work = Sum::default();
        checkpoint()?;
        self.temporal.add_stress_work_with_bias(
            &law.c_eq,
            strains,
            &spectrum[..n],
            bias,
            &mut work,
        );
        for (index, branch) in law.branches.iter().enumerate() {
            checkpoint()?;
            for i in 0..n {
                spectrum[i] = self.source.weights[i][index].primal_binary();
            }
            let gram = branch.resolved_gram_rows(&law.factors[index]);
            let response = self
                .temporal
                .spring(&gram, strains, &spectrum[..n], Some(branch.tau_s))
                .map_err(|_| invalid())?;
            for j in 0..3 {
                storage[j] = storage[j] + Dd::from_f64(response[j]);
            }
            heat = heat + Dd::from_f64(response[3]);
            checkpoint()?;
            self.temporal
                .add_stress_work(&branch.c_pa, strains, &spectrum[..n], &mut work);
        }
        checkpoint()?;
        let kelvin = if let Some(factor) = &law.kv_factor {
            let gram = core::array::from_fn(|i| core::array::from_fn(|j| factor[j][i]));
            self.temporal
                .kelvin_heat(&gram, strains)
                .map_err(|_| invalid())?
        } else {
            0.0
        };
        if let Some(viscosity) = &law.kv_pa_s {
            checkpoint()?;
            for (spectrum, &omega) in spectrum.iter_mut().zip(&self.source.frequencies) {
                *spectrum = [(0.0, 0), (omega, 0)];
            }
            self.temporal
                .add_stress_work(viscosity, strains, &spectrum[..n], &mut work);
        }
        let work = scalar(work.value())?;
        let total_heat = heat + Dd::from_f64(kelvin);
        let result = TensorPronyHarmonicWindowPoint {
            start_stored_energy_j_m3: positive(storage[0].to_f64())?,
            end_stored_energy_j_m3: positive(storage[1].to_f64())?,
            stored_energy_change_j_m3: scalar(storage[2].to_f64())?,
            kelvin_heat_j_m3: kelvin,
            maxwell_heat_j_m3: positive(heat.to_f64())?,
            material_heat_j_m3: positive(total_heat.to_f64())?,
            stress_work_j_m3: work,
            energy_balance_residual_j_m3: scalar(
                (Dd::from_f64(work) - storage[2] - total_heat).to_f64(),
            )?,
        };
        checkpoint()?;
        Ok(result)
    }
}
