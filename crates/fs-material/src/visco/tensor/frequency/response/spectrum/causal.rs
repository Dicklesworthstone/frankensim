//! Original Maxwell memories driven by a prescribed coherent harmonic strain.

use super::*;
use crate::visco::TensorPronyState;
use fs_fft::harmonic_window::{
    HarmonicDecayIntervalLimits, HarmonicPhaseInterval, PreparedHarmonicDecayInterval,
    ProductSpectrumError, phase_angular_interval_from_clock_dd,
};

/// Transactional endpoint of the original raw-Mandel Maxwell history.
/// The accepted input is never changed. `material` keeps interval increments
/// before accumulation into the original f64 cumulative state ledgers.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyCausalHarmonicPoint {
    pub state: TensorPronyState,
    pub material: TensorPronyHarmonicWindowPoint,
}

/// Fixed-source harmonic forcing of the original, arbitrarily retained branch
/// memories. This is a full-tensor strain-controlled interval; a settled
/// plane-stress transverse map does not solve a constrained startup history.
pub struct TensorPronyCausalSpectrumWindow<'p, 'a> {
    source: &'p TensorPronyHarmonicSpectrum<'a>,
    start_s: Dd,
    duration_s: f64,
    temporal: HarmonicEnergySpectrum,
    phases: Vec<HarmonicPhaseInterval>,
    branches: Vec<PreparedHarmonicDecayInterval>,
}

type Channels = [[Dd; 2]; MAX_SPECTRUM_LINES];

fn finite(value: Dd) -> Result<Dd, ViscoError> {
    if value.hi.is_finite() && value.lo.is_finite() { Ok(value) } else { Err(invalid()) }
}

// Preserve the caller's cancellation/error value through the shared moment
// owner's bool cancellation channel, including its bounded inner contractions.
fn numerical<E: From<ViscoError>, T>(
    checkpoint: &mut impl FnMut() -> Result<(), E>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, ProductSpectrumError>,
) -> Result<T, E> {
    let mut failure = None;
    let result = operation(&mut || {
        if failure.is_some() { return true; }
        if let Err(error) = checkpoint() { failure = Some(error); true } else { false }
    });
    if let Some(error) = failure { return Err(error); }
    result.map_err(|_| invalid().into())
}

fn projection(row: &[f64; 6], values: &[f64; 6]) -> Result<Dd, ViscoError> {
    let mut sum = Sum::default();
    for (&a, &b) in row.iter().zip(values) { sum.add([a, b, 1.0, 1.0]); }
    finite(sum.retained_value())
}

fn rate_channels(row: &[f64; 6], strains: &[StrainParts], frequencies: &[f64])
    -> Result<Channels, ViscoError>
{
    let mut result = [[Dd::ZERO; 2]; MAX_SPECTRUM_LINES];
    for (i, (&omega, strain)) in frequencies.iter().zip(strains).enumerate() {
        let (mut re, mut im) = (Sum::default(), Sum::default());
        for (&a, e) in row.iter().zip(strain) {
            for [r, s] in [[e[0].hi, e[1].hi], [e[0].lo, e[1].lo]] {
                re.add([a, omega, s, 1.0]);
                im.add([-a, omega, r, 1.0]);
            }
        }
        result[i] = [finite(re.retained_value())?, finite(im.retained_value())?];
    }
    Ok(result)
}

fn add_scaled(sum: &mut Sum, value: Dd, numerator: f64, denominator: f64) {
    sum.add_ratio(&[value.hi, numerator], &[denominator]);
    sum.add_ratio(&[value.lo, numerator], &[denominator]);
}

impl<'a> TensorPronyHarmonicSpectrum<'a> {
    /// Conservative additional temporal/scratch payload, before allocation.
    /// Includes a 128-KiB bound for fixed point/grouping stack scratch. Excludes
    /// the borrowed source and returned history (48 bytes per branch).
    pub fn causal_window_payload_bound_bytes(&self) -> Result<usize, ViscoError> {
        let n = self.frequencies.len();
        let branch = PreparedHarmonicDecayInterval::required_workspace_bytes(n)
            .map_err(|_| invalid())?;
        self.law.branches.len().checked_mul(branch)
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<TensorPronyCausalSpectrumWindow>()))
            .and_then(|bytes| bytes.checked_add(128 * 1024 + 512 * n + 64 * n * n))
            .ok_or_else(invalid)
    }

    /// Conservative shared-moment work for preparation plus one complete point.
    /// Source limits remain 64 authored lines and the original branch inventory.
    pub fn causal_window_work_bound(&self) -> Result<usize, ViscoError> {
        let n = self.frequencies.len();
        PreparedHarmonicDecayInterval::required_work(n).map_err(|_| invalid())?
            .checked_mul(19 * self.law.branches.len())
            .and_then(|work| work.checked_add(64 * n * n + 1024 * n + 1024))
            .ok_or_else(invalid)
    }

    pub fn prepare_causal_window(&self, window: HarmonicWindow)
        -> Result<TensorPronyCausalSpectrumWindow<'_, 'a>, ViscoError>
    {
        self.prepare_causal_window_with_checkpoint(window, || Ok(()))
    }

    pub fn prepare_causal_window_with_checkpoint<E: From<ViscoError>>(
        &self, window: HarmonicWindow, checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyCausalSpectrumWindow<'_, 'a>, E> {
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(ViscoError::Parameters {
                what: "causal material evolution requires a uniform physical interval",
            }.into());
        }
        self.prepare_causal_window_from_clock_dd_with_checkpoint(
            window.start_clock_dd(), window.duration_s(), checkpoint,
        )
    }

    /// Preserve the accepted source clock across intervals shorter than its
    /// f64 spacing. Frequencies, temperature and original material stay fixed.
    pub fn prepare_causal_window_from_clock_dd(&self, start_s: Dd, duration_s: f64)
        -> Result<TensorPronyCausalSpectrumWindow<'_, 'a>, ViscoError>
    {
        self.prepare_causal_window_from_clock_dd_with_checkpoint(start_s, duration_s, || Ok(()))
    }

    pub fn prepare_causal_window_from_clock_dd_with_checkpoint<E: From<ViscoError>>(
        &self, start_s: Dd, duration_s: f64, mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyCausalSpectrumWindow<'_, 'a>, E> {
        checkpoint()?;
        let start_s = finite(Dd::from_pair(start_s.hi, start_s.lo))?;
        self.causal_window_payload_bound_bytes()?;
        let n = self.frequencies.len();
        let limits = HarmonicDecayIntervalLimits {
            maximum_tones: MAX_SPECTRUM_LINES,
            maximum_workspace_bytes: PreparedHarmonicDecayInterval::required_workspace_bytes(n)
                .map_err(|_| invalid())?,
            maximum_work: PreparedHarmonicDecayInterval::required_work(n).map_err(|_| invalid())?,
        };
        let temporal = HarmonicEnergySpectrum::new_at_clock(&self.frequencies, start_s, duration_s)
            .map_err(|_| invalid())?;
        let mut phases = Vec::with_capacity(n);
        for &omega in &self.frequencies {
            checkpoint()?;
            phases.push(phase_angular_interval_from_clock_dd(omega, start_s, duration_s)
                .map_err(|_| invalid())?);
        }
        let mut branches = Vec::with_capacity(self.law.branches.len());
        for branch in &self.law.branches {
            checkpoint()?;
            branches.push(numerical(&mut checkpoint, |cancelled| {
                PreparedHarmonicDecayInterval::new_cancellable(
                    &self.frequencies, start_s, duration_s, branch.tau_s, limits,
                    &mut || cancelled(),
                )
            })?);
        }
        checkpoint()?;
        Ok(TensorPronyCausalSpectrumWindow { source: self, start_s, duration_s, temporal, phases, branches })
    }
}

impl<'p, 'a> TensorPronyCausalSpectrumWindow<'p, 'a> {
    pub const fn source(&self) -> &'p TensorPronyHarmonicSpectrum<'a> { self.source }
    pub const fn start_time_s_dd(&self) -> Dd { self.start_s }
    pub const fn duration_s(&self) -> f64 { self.duration_s }
    /// Retained temporal payload only; point queries allocate the endpoint q
    /// inventory and otherwise use fixed 64-line stack scratch.
    pub fn numerical_payload_bytes(&self) -> usize {
        core::mem::size_of::<Self>() + self.temporal.numerical_payload_bytes()
            - core::mem::size_of::<HarmonicEnergySpectrum>()
            + self.phases.capacity() * core::mem::size_of::<HarmonicPhaseInterval>()
            + self.branches.iter().map(|branch| branch.numerical_payload_bytes()).sum::<usize>()
            + (self.branches.capacity() - self.branches.len()) * core::mem::size_of::<PreparedHarmonicDecayInterval>()
    }

    pub fn evaluate(&self, strains: &[[C64; 6]], initial: &TensorPronyState)
        -> Result<TensorPronyCausalHarmonicPoint, ViscoError>
    { self.evaluate_with_bias(strains, &[0.0; 6], initial) }

    pub fn evaluate_with_bias(&self, strains: &[[C64; 6]], bias: &[f64; 6], initial: &TensorPronyState)
        -> Result<TensorPronyCausalHarmonicPoint, ViscoError>
    { self.evaluate_with_bias_and_checkpoint(strains, bias, initial, || Ok(())) }

    pub fn evaluate_with_checkpoint<E: From<ViscoError>>(
        &self, strains: &[[C64; 6]], initial: &TensorPronyState,
        checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyCausalHarmonicPoint, E> {
        self.evaluate_with_bias_and_checkpoint(strains, &[0.0; 6], initial, checkpoint)
    }

    /// Advance arbitrary original raw-Mandel spring strains without resetting
    /// them to the settled spectrum. Nonzero DC is a prescribed constant total
    /// strain, not an instruction to erase existing Maxwell memory.
    ///
    /// The public state has no model identifier: callers retain the same law
    /// and branch ordering. Shape, finite ledgers, physical cached storage and
    /// continuity of strain at the interval start are checked. Continuity uses
    /// 128*EPS times the component's original initial/bias/phasor magnitude,
    /// without an absolute-unit tolerance. A Kelvin stress change caused by a
    /// changed strain rate is allowed. Initial transients are not claimed to
    /// have a band-limited spectrum; authored AC/DC still requires source admission.
    pub fn evaluate_with_bias_and_checkpoint<E: From<ViscoError>>(
        &self, strains: &[[C64; 6]], bias: &[f64; 6], initial: &TensorPronyState,
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<TensorPronyCausalHarmonicPoint, E> {
        checkpoint()?;
        let grouped = self.source.group_strains(strains)?;
        let n = self.source.frequencies.len();
        let strains = &grouped[..n];
        let bias = bias.map(Dd::from_f64);
        self.source.admit_grouped_strains_with_bias(strains, &bias)?;
        let law = self.source.law;
        if initial.q.len() != law.branches.len()
            || initial.q.iter().flatten().chain(&initial.strain).chain(&initial.stress)
                .any(|value| !value.is_finite())
            || !initial.work.is_finite() || !initial.stored.is_finite() || initial.stored < 0.0
            || !initial.dissipated.is_finite() || initial.dissipated < 0.0
        {
            return Err(ViscoError::Parameters { what: "invalid original Mandel Maxwell history" }.into());
        }
        law.validate_mandel_strain(&initial.strain)?;
        let mut endpoints = [[Dd::ZERO; 6]; 2];
        let mut end_rate = [Dd::ZERO; 6];
        for component in 0..6 {
            let mut sums = [Sum::default(); 3];
            for sum in &mut sums[..2] { add_scaled(sum, bias[component], 1.0, 1.0); }
            let mut tolerance = 128.0 * f64::EPSILON * initial.strain[component].abs()
                + 128.0 * f64::EPSILON * bias[component].to_f64().abs();
            for i in 0..n {
                let z = strains[i][component];
                for (slot, phase) in [self.phases[i].start, self.phases[i].end].iter().enumerate() {
                    add_scaled(&mut sums[slot], z[0] * phase[0] - z[1] * phase[1], 1.0, 1.0);
                }
                let phase = self.phases[i].end;
                add_scaled(&mut sums[2], z[1] * phase[0] + z[0] * phase[1], self.source.frequencies[i], 1.0);
                for value in z { tolerance += 128.0 * f64::EPSILON * (value.hi.abs() + value.lo.abs()); }
            }
            endpoints[0][component] = finite(sums[0].retained_value())?;
            endpoints[1][component] = finite(sums[1].retained_value())?;
            end_rate[component] = finite(sums[2].retained_value())?;
            if (endpoints[0][component] - Dd::from_f64(initial.strain[component])).to_f64().abs() > tolerance {
                return Err(ViscoError::Parameters { what: "harmonic strain is discontinuous at the accepted history clock" }.into());
            }
        }
        let next_strain = endpoints[1].map(Dd::to_f64);
        law.validate_mandel_strain(&next_strain)?;
        let identity = [[(1.0, 0), (0.0, 0)]; MAX_SPECTRUM_LINES];
        let equilibrium = self.temporal.spring_with_bias(
            &law.equilibrium_gram, strains, &identity[..n], None, &bias,
        ).map_err(|_| invalid())?;
        let mut storage = [equilibrium[0], equilibrium[1], equilibrium[2]].map(Dd::from_f64);
        let mut work = Sum::default();
        self.temporal.add_stress_work_with_bias(&law.c_eq, strains, &identity[..n], &bias, &mut work);
        let mut stress = [Sum::default(); 6];
        for i in 0..6 { for j in 0..6 { add_scaled(&mut stress[i], endpoints[1][j], law.c_eq[i][j], 1.0); } }
        let mut heat = Sum::default();
        let zero = [[Dd::ZERO; 2]; MAX_SPECTRUM_LINES];
        let mut next_q = Vec::with_capacity(law.branches.len());
        for (b, (branch, kernel)) in law.branches.iter().zip(&self.branches).enumerate() {
            checkpoint()?;
            let mut q = [0.0; 6];
            for i in 0..6 {
                let mut axis = [0.0; 6]; axis[i] = 1.0;
                let rate = rate_channels(&axis, strains, &self.source.frequencies)?;
                // q_i * (C^T strain_rate)_i is the original stress work.
                // Use the source column explicitly and share the raw-q query;
                // squaring Cq solely to obtain work would impose an unrelated
                // stress-squared range restriction on a valid material.
                let column = core::array::from_fn(|j| branch.c_pa[j][i]);
                let work_rate = rate_channels(&column, strains, &self.source.frequencies)?;
                let moment = numerical(&mut checkpoint, |cancelled| kernel.evaluate_cancellable(
                    Dd::from_f64(initial.q[b][i]), &rate[..n], &work_rate[..n], &mut || cancelled(),
                ))?;
                q[i] = scalar(moment.end.to_f64())?;
                add_scaled(&mut work, moment.mean_mixed, self.duration_s, 1.0);
            }
            for row in branch.resolved_gram_rows(&law.factors[b]) {
                let initial_channel = projection(&row, &initial.q[b])?;
                let drive = rate_channels(&row, strains, &self.source.frequencies)?;
                let moment = numerical(&mut checkpoint, |cancelled| kernel.evaluate_cancellable(
                    initial_channel, &drive[..n], &zero[..n], &mut || cancelled(),
                ))?;
                storage[0] = storage[0] + Dd::from_f64(0.5) * initial_channel * initial_channel;
                storage[1] = storage[1] + Dd::from_f64(0.5) * moment.end * moment.end;
                storage[2] = storage[2] + initial_channel * moment.increment
                    + Dd::from_f64(0.5) * moment.increment * moment.increment;
                add_scaled(&mut heat, moment.mean_square, self.duration_s, branch.tau_s);
            }
            for i in 0..6 { for j in 0..6 { stress[i].add([branch.c_pa[i][j], q[j], 1.0, 1.0]); } }
            next_q.push(q);
        }
        checkpoint()?;
        let kelvin = if let Some(factor) = &law.kv_factor {
            let gram = core::array::from_fn(|i| core::array::from_fn(|j| factor[j][i]));
            self.temporal.kelvin_heat(&gram, strains).map_err(|_| invalid())?
        } else { 0.0 };
        if let Some(viscosity) = &law.kv_pa_s {
            let mut spectrum = [[(0.0, 0), (0.0, 0)]; MAX_SPECTRUM_LINES];
            for (weight, &omega) in spectrum.iter_mut().zip(&self.source.frequencies) { weight[1] = (omega, 0); }
            self.temporal.add_stress_work(viscosity, strains, &spectrum[..n], &mut work);
            for i in 0..6 { for j in 0..6 { add_scaled(&mut stress[i], end_rate[j], viscosity[i][j], 1.0); } }
        }
        let start_stored = positive(storage[0].to_f64())?;
        if (start_stored - initial.stored).abs()
            > 128.0 * f64::EPSILON * start_stored + 128.0 * f64::EPSILON * initial.stored
        {
            return Err(ViscoError::Parameters { what: "accepted Maxwell history has inconsistent stored energy" }.into());
        }
        let total_heat = heat.retained_value() + Dd::from_f64(kelvin);
        let work = finite(work.retained_value())?;
        let material = TensorPronyHarmonicWindowPoint {
            start_stored_energy_j_m3: start_stored,
            end_stored_energy_j_m3: positive(storage[1].to_f64())?,
            stored_energy_change_j_m3: scalar(storage[2].to_f64())?,
            kelvin_heat_j_m3: kelvin,
            maxwell_heat_j_m3: positive(heat.value())?,
            material_heat_j_m3: positive(total_heat.to_f64())?,
            stress_work_j_m3: scalar(work.to_f64())?,
            energy_balance_residual_j_m3: scalar((work - storage[2] - total_heat).to_f64())?,
        };
        let state = TensorPronyState {
            q: next_q, strain: next_strain,
            stress: stress.map(Sum::value), stored: material.end_stored_energy_j_m3,
            work: scalar((Dd::from_f64(initial.work) + work).to_f64())?,
            dissipated: positive((Dd::from_f64(initial.dissipated) + total_heat).to_f64())?,
        };
        if state.stress.iter().any(|value| !value.is_finite()) { return Err(invalid().into()); }
        checkpoint()?;
        Ok(TensorPronyCausalHarmonicPoint { state, material })
    }
}
