//! Explicit mixed real/proper forcing through one original mechanical pencil.
//!
//! Factor amplitudes carry their square-root statistical weights. One factor
//! spans all its coordinates and frequencies; different factors have zero
//! cross moments. No covariance, independence distribution, Gaussian law or
//! stationary spectrum is inferred from a list of frequencies.

use super::{
    HarmonicForceTone, HarmonicInterferenceLoad, HarmonicInterferenceMoments, HarmonicResponse,
    HarmonicResponseControls, HarmonicResponseError, MAX_INTERFERENCE_TONES, Tone, admit_size,
    apply, cap, checked, finish_moments, finite, invalid, load_payload, observe_sums, parts, poll,
    prepare_response,
};
use fs_exec::Cx;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::{c64::C64, dd::Dd};

/// Hard ceiling on explicitly authored random factors, excluding the mean.
pub const MAX_INTERFERENCE_FACTORS: usize = 64;
/// Hard ceiling on authored tones summed across the mean and every factor.
pub const MAX_ENSEMBLE_FORCE_TONES: usize = 4096;

/// Authored joint second-moment law for one factor's complete force signal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HarmonicFactorStatistics {
    /// A zero-mean real variable of unit variance multiplies the real signal.
    /// Different real factors have E[xi_r*xi_s] = delta_rs. Both sum-frequency
    /// and difference-frequency products survive within each factor.
    Real,
    /// A zero-mean proper complex variable multiplies the complete analytic
    /// signal before taking its real part. Jointly, E[z_r*conj(z_s)] = delta_rs
    /// and E[z_r*z_s] = 0. Only difference-frequency products survive.
    ProperComplex,
}

/// One covariance factor, retaining coherence across coordinates and tones.
///
/// The physical force is its deterministic mean plus each real factor times
/// the real tone sum, plus the real part of each proper factor times its
/// complex tone sum. All real/proper cross moments are zero. Coefficients
/// contain physical force units and square-root statistical weights. Factors
/// with shared source dependence must be encoded jointly, not declared as
/// separate uncorrelated factors. Empty factors are explicit zero columns.
///
/// Only these second moments are used. No Gaussian or fourth-moment law,
/// realization generator, energy variance, or source PSD normalization is added.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicForceFactor<'a> {
    /// Peak force phasors on the common source clock, in original coordinates.
    pub tones: &'a [HarmonicForceTone<'a>],
    /// Explicit real or jointly proper complex factor law.
    pub statistics: HarmonicFactorStatistics,
}

/// Aggregate admission for deterministic bias plus all random force factors.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicEnsembleControls {
    /// Maximum random factors, capped by MAX_INTERFERENCE_FACTORS.
    pub maximum_factors: usize,
    /// Maximum authored tones across all banks, including exact duplicates.
    /// Capped by MAX_ENSEMBLE_FORCE_TONES; each bank also has the original
    /// MAX_INTERFERENCE_TONES ceiling.
    pub maximum_tones: usize,
    /// Conservative peak scalar arrays and headers: retained responses/actions,
    /// ordering records, bank headers and one existing dense response owner.
    /// Caller inputs, allocator overhead and thread stacks are excluded.
    pub maximum_payload_bytes: usize,
    /// Maximum products: `(6*n + loaded_ports)*sum_b(distinct_tones_in_bank_b^2)`.
    /// Different factors have no cross contractions. Checked before any solve.
    pub maximum_product_terms: usize,
}

#[derive(Clone, Debug)]
struct Bank {
    statistics: HarmonicFactorStatistics,
    tones: Vec<Tone>,
}

#[derive(Clone, Copy)]
struct Input {
    omega: f64,
    bank: usize,
    index: usize,
}

/// Expected finite-window mechanical energy and work under authored coherence.
///
/// One original response owner is prepared per distinct frequency globally and
/// reused for every factor RHS at that frequency. Exact duplicates are summed
/// in DD separately inside each bank before solving; independent factors are
/// never combined. Full coupled M/C/K and the source clock remain unchanged.
///
/// Each observation reuses original-coordinate M*v, K*q and D*v actions, with
/// D=(C+C^T)/2. It allocates nothing, performs no solve or matrix application,
/// and accumulates the shared physical contractions in DD across banks,
/// including signed boundary work when a physical impedance is supplied.
///
/// Results are ensemble expectations of the actual sensitivity-weighted real
/// products. Real factors retain sum-frequency terms; proper factors retain
/// every within-factor beat. Singular frequencies, unrepresentable arithmetic
/// and response residuals retain the original owner's refusal. This is the
/// particular linear harmonic motion, not startup, nonlinear mixing, inferred
/// source statistics, or a near-null relative-error certificate.
#[derive(Clone, Debug)]
pub struct PreparedHarmonicInterferenceEnsemble {
    n: usize,
    banks: Vec<Bank>,
    frequency_count: usize,
    maximum_relative_response_residual: f64,
}

impl PreparedHarmonicInterferenceEnsemble {
    /// Prepare explicit deterministic mean and random covariance factors.
    ///
    /// Mean-only and factor-only sources are supported. At least one authored
    /// tone is required; no dummy frequency is invented for an empty source.
    /// Zero-amplitude authored tones still receive normal pencil admission.
    #[allow(clippy::too_many_arguments)] // Original pencil plus explicit statistical model and limits.
    pub fn try_new<'a>(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        mean_tones: &[HarmonicForceTone<'a>],
        factors: &[HarmonicForceFactor<'a>],
        response_controls: HarmonicResponseControls,
        controls: HarmonicEnsembleControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        Self::prepare(
            mass,
            damping,
            stiffness,
            n,
            mean_tones,
            factors,
            None,
            response_controls,
            controls,
            cx,
        )
    }

    /// Prepare every bank through the same physical frequency-dependent load.
    /// One loaded factorization per distinct frequency serves every RHS there;
    /// cross-frequency boundary work follows each factor's authored statistics.
    #[allow(clippy::too_many_arguments)]
    pub fn try_with_port_impedance<'a>(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        mean_tones: &[HarmonicForceTone<'a>],
        factors: &[HarmonicForceFactor<'a>],
        load: HarmonicInterferenceLoad<'_>,
        response_controls: HarmonicResponseControls,
        controls: HarmonicEnsembleControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        Self::prepare(
            mass,
            damping,
            stiffness,
            n,
            mean_tones,
            factors,
            Some(load),
            response_controls,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare<'a>(
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        n: usize,
        mean_tones: &[HarmonicForceTone<'a>],
        factors: &[HarmonicForceFactor<'a>],
        load: Option<HarmonicInterferenceLoad<'_>>,
        response_controls: HarmonicResponseControls,
        controls: HarmonicEnsembleControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        admit_size(n, response_controls)?;
        cap(
            "interference factors",
            factors.len(),
            controls.maximum_factors.min(MAX_INTERFERENCE_FACTORS),
        )?;
        let extent = n.checked_mul(n).ok_or(invalid("pencil extent overflow"))?;
        if [mass, damping, stiffness].iter().any(|a| a.len() != extent) {
            return Err(invalid("one complete M/C/K pencil is required"));
        }
        let bank_count = factors.len() + 1;
        let input = |bank: usize| {
            if bank == 0 {
                mean_tones
            } else {
                factors[bank - 1].tones
            }
        };
        let mut authored = 0usize;
        for bank in 0..bank_count {
            poll(cx)?;
            let source = input(bank);
            cap(
                "interference bank tones",
                source.len(),
                MAX_INTERFERENCE_TONES,
            )?;
            authored = authored
                .checked_add(source.len())
                .ok_or(invalid("ensemble tone count overflow"))?;
            cap(
                "interference ensemble tones",
                authored,
                controls.maximum_tones.min(MAX_ENSEMBLE_FORCE_TONES),
            )?;
            for tone in source {
                if !tone.angular_frequency_rad_s.is_finite()
                    || tone.angular_frequency_rad_s < 0.0
                    || tone.force.len() != n
                    || tone.force.iter().any(|&z| !finite(z))
                {
                    return Err(invalid(
                        "finite nonnegative tone frequencies and complete forces required",
                    ));
                }
            }
        }
        if authored == 0 {
            return Err(invalid("at least one ensemble force tone is required"));
        }
        if let Some(load) = load {
            load.validate(n, authored, cx)?;
        }
        // Admit sorting storage before its first allocation; retain the exact
        // response count only after duplicates have been identified.
        admit_payload(n, authored, bank_count, 0, load, controls)?;
        let mut order = Vec::with_capacity(authored);
        for bank in 0..bank_count {
            poll(cx)?;
            for (index, tone) in input(bank).iter().enumerate() {
                order.push(Input {
                    // Positive and negative zero are the same physical DC.
                    omega: if tone.angular_frequency_rad_s == 0.0 {
                        0.0
                    } else {
                        tone.angular_frequency_rad_s
                    },
                    bank,
                    index,
                });
            }
        }
        // Nonallocating sort, bounded by the hard authored-tone ceiling.
        order.sort_unstable_by(|a, b| {
            a.omega
                .total_cmp(&b.omega)
                .then_with(|| a.bank.cmp(&b.bank))
                .then_with(|| a.index.cmp(&b.index))
        });
        let mut counts = [0usize; MAX_INTERFERENCE_FACTORS + 1];
        let mut retained = 0usize;
        let mut frequency_count = 0usize;
        for (i, item) in order.iter().enumerate() {
            if i % 64 == 0 {
                poll(cx)?;
            }
            if i == 0 || item.omega != order[i - 1].omega {
                frequency_count += 1;
            }
            if i == 0 || item.omega != order[i - 1].omega || item.bank != order[i - 1].bank {
                counts[item.bank] += 1;
                retained += 1;
            }
        }
        if let Some(load) = load {
            if load.samples.len() != frequency_count {
                return Err(invalid(
                    "exactly one load sample per distinct source frequency required",
                ));
            }
            for item in &order {
                load.impedance(item.omega)?;
            }
        }
        admit_payload(n, authored, bank_count, retained, load, controls)?;
        let pair_count = counts[..bank_count].iter().try_fold(0usize, |sum, &count| {
            sum.checked_add(count.checked_mul(count)?)
        });
        let products = pair_count
            .and_then(|pairs| {
                pairs.checked_mul(n.checked_mul(6)?.checked_add(load.map_or(0, |l| l.ports))?)
            })
            .ok_or(invalid("ensemble product extent overflow"))?;
        cap(
            "interference ensemble product terms",
            products,
            controls.maximum_product_terms,
        )?;
        let mut banks = Vec::with_capacity(bank_count);
        for (bank, &count) in counts[..bank_count].iter().enumerate() {
            banks.push(Bank {
                statistics: if bank == 0 {
                    HarmonicFactorStatistics::Real
                } else {
                    factors[bank - 1].statistics
                },
                tones: Vec::with_capacity(count),
            });
        }

        let mut maximum_relative_response_residual = 0.0_f64;
        let mut first = 0;
        while first < authored {
            poll(cx)?;
            let omega = order[first].omega;
            let mut frequency_end = first + 1;
            while frequency_end < authored && order[frequency_end].omega == omega {
                frequency_end += 1;
            }
            let prepared = prepare_response(
                mass,
                damping,
                stiffness,
                n,
                omega,
                load,
                response_controls,
                cx,
            )?;
            let mut group = first;
            while group < frequency_end {
                poll(cx)?;
                let bank = order[group].bank;
                let source = input(bank);
                let mut end = group + 1;
                while end < frequency_end && order[end].bank == bank {
                    end += 1;
                }
                let mut force = Vec::with_capacity(n);
                for i in 0..n {
                    if i % 64 == 0 {
                        poll(cx)?;
                    }
                    let mut sum = [Dd::ZERO; 2];
                    for item in &order[group..end] {
                        let value = parts(source[item.index].force[i]);
                        sum[0] = checked(sum[0] + value[0])?;
                        sum[1] = checked(sum[1] + value[1])?;
                    }
                    force.push(C64::new(sum[0].to_f64(), sum[1].to_f64()));
                }
                let response = prepared.solve(&force, cx)?;
                maximum_relative_response_residual =
                    maximum_relative_response_residual.max(response.relative_residual);
                let mass_velocity = apply(mass, &response.velocity, n, false, cx)?;
                let stiffness_displacement =
                    apply(stiffness, &response.displacement, n, false, cx)?;
                let viscous_velocity = apply(damping, &response.velocity, n, true, cx)?;
                banks[bank].tones.push(Tone {
                    omega,
                    force,
                    response,
                    mass_velocity,
                    stiffness_displacement,
                    viscous_velocity,
                });
                group = end;
            }
            first = frequency_end;
        }
        poll(cx)?;
        Ok(Self {
            n,
            banks,
            frequency_count,
            maximum_relative_response_residual,
        })
    }

    /// Number of authored factors, including any explicitly empty zero columns.
    #[must_use]
    pub fn factor_count(&self) -> usize {
        self.banks.len() - 1
    }

    /// Distinct frequencies globally; one response preparation is reused at each.
    #[must_use]
    pub const fn frequency_count(&self) -> usize {
        self.frequency_count
    }

    /// Distinct frequencies in the deterministic mean, possibly zero.
    #[must_use]
    pub fn mean_tone_count(&self) -> usize {
        self.banks[0].tones.len()
    }

    /// Distinct frequencies in one random factor, preserving its original index.
    #[must_use]
    pub fn factor_tone_count(&self, factor: usize) -> Option<usize> {
        self.factor(factor).map(|bank| bank.tones.len())
    }

    /// Original authored factor law; the deterministic mean has no factor index.
    #[must_use]
    pub fn factor_statistics(&self, factor: usize) -> Option<HarmonicFactorStatistics> {
        self.factor(factor).map(|bank| bank.statistics)
    }

    /// Immutable deterministic-mean response, in increasing frequency order.
    #[must_use]
    pub fn mean_tone_response(&self, tone: usize) -> Option<(f64, &HarmonicResponse)> {
        self.banks[0]
            .tones
            .get(tone)
            .map(|t| (t.omega, &t.response))
    }

    /// Immutable response of one factor, in increasing frequency order.
    #[must_use]
    pub fn factor_tone_response(
        &self,
        factor: usize,
        tone: usize,
    ) -> Option<(f64, &HarmonicResponse)> {
        self.factor(factor)?
            .tones
            .get(tone)
            .map(|t| (t.omega, &t.response))
    }

    /// Maximum original-pencil backward residual over the mean and all factors.
    #[must_use]
    pub const fn maximum_relative_response_residual(&self) -> f64 {
        self.maximum_relative_response_residual
    }

    /// Expected window-averaged energy/work of the complete random real motion.
    ///
    /// Cross-frequency coherence stays inside its original factor. Different
    /// factors have no cross terms under the declared joint moments. Raw DD
    /// contributions accumulate before one final rounding; covariance is not
    /// recovered by subtracting a large deterministic mean. No allocation,
    /// matrix application, sampling, or solve occurs in this observation.
    pub fn observe(
        &self,
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<HarmonicInterferenceMoments, HarmonicResponseError> {
        poll(cx)?;
        let mut sums = [Dd::ZERO; 7];
        for bank in &self.banks {
            poll(cx)?;
            let values = observe_sums(self.n, &bank.tones, window, bank.statistics, cx)?;
            for (sum, value) in sums.iter_mut().zip(values) {
                *sum = checked(*sum + value)?;
            }
        }
        let result = finish_moments(sums)?;
        poll(cx)?;
        Ok(result)
    }

    fn factor(&self, factor: usize) -> Option<&Bank> {
        self.banks.get(factor.checked_add(1)?)
    }
}

fn admit_payload(
    n: usize,
    authored: usize,
    banks: usize,
    retained: usize,
    load: Option<HarmonicInterferenceLoad<'_>>,
    controls: HarmonicEnsembleControls,
) -> Result<(), HarmonicResponseError> {
    let additional_load = load_payload(n, retained, load)?;
    let bytes = n
        .checked_mul(n)
        .and_then(|x| x.checked_mul(144))
        .and_then(|x| x.checked_add(n.checked_mul(192)?))
        .and_then(|x| x.checked_add(authored.checked_mul(core::mem::size_of::<Input>())?))
        .and_then(|x| x.checked_add(banks.checked_mul(core::mem::size_of::<Bank>())?))
        .and_then(|x| {
            x.checked_add(
                retained.checked_mul(
                    n.checked_mul(144)?
                        .checked_add(core::mem::size_of::<Tone>())?,
                )?,
            )
        })
        .and_then(|x| x.checked_add(additional_load))
        .ok_or(invalid("ensemble payload extent overflow"))?;
    cap(
        "interference ensemble payload bytes",
        bytes,
        controls.maximum_payload_bytes,
    )
}
