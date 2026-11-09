//! Harmonic driving of the original coupled plane-stress memory realization.

use super::*;
use crate::plane_stress_viscosity::MaterialPlane;
use crate::visco::tensor::{HarmonicEnergySpectrum, HarmonicWorkSum as Sum, StrainParts};
use crate::visco::{
    MAX_HARMONIC_MATERIAL_LINES, TensorPronyHarmonicSpectrum, TensorPronyHarmonicWindowPoint,
};
use fs_fft::harmonic_window::{
    HarmonicDecayIntervalLimits, HarmonicPhaseInterval, PreparedHarmonicDecayInterval,
    ProductSpectrumError, phase_angular_interval_from_clock_dd,
};
use fs_la::factor::{Qr, qr};
use fs_math::c64::C64;

const LINES: usize = MAX_HARMONIC_MATERIAL_LINES;
const ENGINEERING: [f64; 3] = [1.0, 1.0, core::f64::consts::FRAC_1_SQRT_2];
type PlaneParts = [[Dd; 2]; 3];

/// Additional storage and preparation-plus-one-point work, after the original
/// mixed controller's independently budgeted QR/SVD preparation.
#[derive(Clone, Copy, Debug)]
pub struct TensorPlaneStressHarmonicBudget {
    pub maximum_payload_bytes: usize,
    pub maximum_scalar_work: usize,
}

/// Separate candidate in the original source-bound, Gram-coordinate history.
/// No accepted history changes until its caller adopts this result.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPlaneStressHarmonicPoint {
    pub history: TensorPronyHistory,
    pub material: TensorPronyHarmonicWindowPoint,
    /// Direct full strain increment, in tensor-shear Voigt coordinates.
    pub strain_increment: Voigt,
    /// Conservative full-tensor bounds over the complete observed interval.
    pub strain_bounds: [Voigt; 2],
}

/// Reusable causal loading on one original plane-stress partition, with fixed
/// material and temperature. Input phasors use exp(-i omega t) and ordered
/// engineering coordinates [epsilon_aa, epsilon_bb, gamma_ab]. Transverse
/// strain follows the coupled ORIGINAL mixed-control modes, including initial
/// memory, rather than a settled frequency-specific transverse map.
pub struct TensorPlaneStressHarmonicWindow<'p> {
    control: &'p PreparedTensorStressControl,
    source: TensorPronyHarmonicSpectrum<'p>,
    plane: MaterialPlane,
    inside: [usize; 3],
    constraint: EquilibriumConstraint,
    map: [Voigt; 3],
    transverse: Vec<Voigt>,
    work_rows: Vec<[f64; 3]>,
    kernels: Vec<PreparedHarmonicDecayInterval>,
    temporal: HarmonicEnergySpectrum,
    phases: Vec<HarmonicPhaseInterval>,
    start_s: Dd,
    duration_s: f64,
    payload_bound: usize,
}

struct EquilibriumConstraint {
    axes: [usize; 3],
    factor: Qr,
    scale: f64,
}

fn moment_error(error: ProductSpectrumError) -> TensorPronyError {
    match error {
        ProductSpectrumError::Cancelled => TensorPronyError::Cancelled,
        _ => invalid("causal constrained harmonic moment is unresolved"),
    }
}
fn number(value: Dd) -> Result<Dd> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable plane-stress harmonic response"))
    }
}
fn add(sum: &mut Sum, value: Dd, factor: f64) {
    sum.add([value.hi, factor, 1.0, 1.0]);
    sum.add([value.lo, factor, 1.0, 1.0]);
}
fn l1(value: [Dd; 2]) -> f64 {
    value
        .iter()
        .map(|value| value.hi.abs() + value.lo.abs())
        .sum()
}
fn equilibrium_stress(model: &TensorProny, strain: &[Dd; 6]) -> [Dd; 6] {
    let channels = model
        .equilibrium_gram
        .map(|row| (0..6).fold(Dd::ZERO, |sum, i| sum + Dd::from_f64(row[i]) * strain[i]));
    core::array::from_fn(|i| {
        (0..6).fold(Dd::ZERO, |sum, row| {
            sum + Dd::from_f64(model.equilibrium_gram[row][i]) * channels[row]
        })
    })
}

impl EquilibriumConstraint {
    fn new(model: &TensorProny, axes: [usize; 3]) -> Result<Self> {
        let scale = model
            .equilibrium_gram
            .iter()
            .flat_map(|row| axes.map(|i| row[i].abs()))
            .fold(0.0_f64, f64::max);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid("unresolved constrained equilibrium factor"));
        }
        let values: Vec<_> = model
            .equilibrium_gram
            .iter()
            .flat_map(|row| axes.map(|i| row[i] / scale))
            .collect();
        let factor = qr(&values, 6, 3);
        if (0..3).any(|i| !factor.r(i, i).is_finite() || factor.r(i, i) == 0.0) {
            return Err(invalid("singular constrained equilibrium factor"));
        }
        Ok(Self {
            axes,
            factor,
            scale,
        })
    }
    fn solve(&self, rhs: &[Dd; 6]) -> Result<[Dd; 6]> {
        let mut x = [Dd::ZERO; 3];
        for i in 0..3 {
            let mut value = rhs[self.axes[i]] / Dd::from_f64(self.scale);
            for j in 0..i {
                value = value - Dd::from_f64(self.factor.r(j, i)) * x[j];
            }
            x[i] = number(value / Dd::from_f64(self.factor.r(i, i)))?;
        }
        for i in (0..3).rev() {
            let mut value = x[i] / Dd::from_f64(self.scale);
            for j in i + 1..3 {
                value = value - Dd::from_f64(self.factor.r(i, j)) * x[j];
            }
            x[i] = number(value / Dd::from_f64(self.factor.r(i, i)))?;
        }
        let mut full = [Dd::ZERO; 6];
        for i in 0..3 {
            full[self.axes[i]] = x[i];
        }
        Ok(full)
    }
}

impl<'p> TensorPlaneStressHarmonicWindow<'p> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        control: &'p PreparedTensorStressControl,
        plane: MaterialPlane,
        frequencies: &[f64],
        start_s: Dd,
        duration_s: f64,
        budget: TensorPlaneStressHarmonicBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let (inside, outside) = match plane {
            MaterialPlane::Axes12 => ([0, 1, 3], [2, 4, 5]),
            MaterialPlane::Axes21 => ([1, 0, 3], [2, 4, 5]),
            MaterialPlane::Axes13 => ([0, 2, 5], [1, 3, 4]),
            MaterialPlane::Axes23 => ([1, 2, 4], [0, 3, 5]),
        };
        if (0..6).any(|i| control.stress_components[i] != outside.contains(&i)) {
            return Err(invalid(
                "harmonic plane stress needs its matching original mixed-control partition",
            ));
        }
        let f = frequencies.len();
        let n = control.memory_coordinates();
        let per_kernel =
            PreparedHarmonicDecayInterval::required_workspace_bytes(f).map_err(moment_error)?;
        let kernel_work = PreparedHarmonicDecayInterval::required_work(f).map_err(moment_error)?;
        let payload_bound = n
            .checked_mul(per_kernel + 256)
            .and_then(|v| {
                v.checked_add(
                    core::mem::size_of::<Self>()
                        + 128 * 1024
                        + 128 * f * (control.model.branches().len() + 1)
                        + 512 * f
                        + 64 * f * f,
                )
            })
            .ok_or_else(|| invalid("harmonic plane-stress payload extent overflow"))?;
        let work = n
            .checked_mul(kernel_work)
            .and_then(|v| v.checked_add(64 * f * f + 1024 * f * (n + 1) + 32 * n * n + 4096))
            .ok_or_else(|| invalid("harmonic plane-stress work extent overflow"))?;
        if payload_bound > budget.maximum_payload_bytes || work > budget.maximum_scalar_work {
            return Err(invalid(
                "harmonic plane stress exceeds its explicit payload or work budget",
            ));
        }
        let source = control
            .model
            .prepare_harmonic_spectrum(frequencies)
            .map_err(|_| invalid("harmonic carriers leave the original tensor material domain"))?;
        let start_s = number(Dd::from_pair(start_s.hi, start_s.lo))?;
        let temporal = HarmonicEnergySpectrum::new_at_clock(
            source.distinct_frequencies(),
            start_s,
            duration_s,
        )
        .map_err(|_| invalid("invalid physical harmonic interval"))?;
        let constraint = EquilibriumConstraint::new(&control.model, outside)?;
        let mut map = [[0.0; 6]; 3];
        for j in 0..3 {
            let mut unit = [Dd::ZERO; 6];
            unit[inside[j]] = Dd::from_f64(ENGINEERING[j]);
            let rhs = equilibrium_stress(&control.model, &unit).map(|v| -v);
            let extension = constraint.solve(&rhs)?;
            map[j] = core::array::from_fn(|i| (unit[i] + extension[i]).to_f64());
        }
        let mut transverse = Vec::with_capacity(n);
        let mut work_rows = Vec::with_capacity(n);
        let mut memory = vec![[0.0; 6]; control.model.branches().len()];
        for mode in 0..n {
            poll(&mut cancelled)?;
            for row in 0..n {
                memory[row / 6][row % 6] =
                    control.modes[row * n + mode] / control.branch_rate_roots[row / 6];
            }
            let (stress, _) = control.source_stress(&[0.0; 6], &memory)?;
            let free = constraint
                .solve(&stress.map(|v| Dd::from_f64(-v)))?
                .map(Dd::to_f64);
            let (reaction, scale) = control.source_stress(&free, &memory)?;
            if outside
                .iter()
                .any(|&i| reaction[i].abs() > 8192.0 * f64::EPSILON * (n + 6) as f64 * scale[i])
            {
                return Err(invalid(
                    "constrained harmonic memory reconstruction is unresolved",
                ));
            }
            work_rows.push(core::array::from_fn(|j| {
                reaction[inside[j]] * ENGINEERING[j]
            }));
            transverse.push(free);
        }
        let limits = HarmonicDecayIntervalLimits {
            maximum_tones: LINES,
            maximum_workspace_bytes: per_kernel,
            maximum_work: kernel_work,
        };
        let mut kernels = Vec::with_capacity(n);
        for &rate in &control.rates {
            poll(&mut cancelled)?;
            kernels.push(
                PreparedHarmonicDecayInterval::new_cancellable(
                    source.distinct_frequencies(),
                    start_s,
                    duration_s,
                    1.0 / rate,
                    limits,
                    &mut cancelled,
                )
                .map_err(moment_error)?,
            );
        }
        let phases = source
            .distinct_frequencies()
            .iter()
            .map(|&omega| {
                phase_angular_interval_from_clock_dd(omega, start_s, duration_s)
                    .map_err(|_| invalid("unrepresentable original harmonic phase"))
            })
            .collect::<Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(Self {
            control,
            source,
            plane,
            inside,
            constraint,
            map,
            transverse,
            work_rows,
            kernels,
            temporal,
            phases,
            start_s,
            duration_s,
            payload_bound,
        })
    }

    pub const fn plane(&self) -> MaterialPlane {
        self.plane
    }
    pub const fn start_time_s_dd(&self) -> Dd {
        self.start_s
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.source.angular_frequencies_rad_s()
    }
    /// Conservative numeric storage including point scratch and returned history;
    /// the separately owned original mixed-controller factors are excluded.
    pub const fn numerical_payload_bound_bytes(&self) -> usize {
        self.payload_bound
    }

    /// Drive an already constrained history with engineering-strain harmonics.
    /// Equal carriers combine before recovery/admission. The full imposed
    /// waveform and a conservative full-tensor transient envelope are admitted.
    /// This is not a maximum-strain certificate or a temperature evolution.
    /// Initial transverse stress must already vanish; no initial jump is hidden.
    /// Refusal/cancellation leaves the original source-bound history unchanged.
    pub fn evaluate(
        &self,
        strains: &[[C64; 3]],
        bias: &[f64; 3],
        history: &TensorPronyHistory,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorPlaneStressHarmonicPoint> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.control.model, &history.model) || bias.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "harmonic plane stress requires its original history and finite strain bias",
            ));
        }
        let grouped = self
            .source
            .group_strains(strains)
            .map_err(|_| invalid("invalid coherent in-plane harmonic strains"))?;
        let frequencies = self.source.distinct_frequencies();
        let f = frequencies.len();
        let n = self.control.memory_coordinates();
        let tolerance = 8192.0 * f64::EPSILON * (n + 6) as f64;
        let (initial_stress, initial_scale) =
            self.control.source_stress(&history.strain, &history.q)?;
        if self
            .constraint
            .axes
            .iter()
            .any(|&i| initial_stress[i].abs() > tolerance * initial_scale[i])
        {
            return Err(invalid(
                "initial history is not in the prescribed plane-stress state",
            ));
        }
        let mut mapped: [StrainParts; LINES] = [[[Dd::ZERO; 2]; 6]; LINES];
        let mut mapped_bias = [Dd::ZERO; 6];
        for axis in 0..6 {
            let mut sum = Sum::default();
            for j in 0..3 {
                sum.add([self.map[j][axis], bias[j], 1.0, 1.0]);
            }
            mapped_bias[axis] = number(sum.retained_value())?;
            for k in 0..f {
                for part in 0..2 {
                    let mut sum = Sum::default();
                    for j in 0..3 {
                        add(&mut sum, grouped[k][j][part], self.map[j][axis]);
                    }
                    mapped[k][axis][part] = number(sum.retained_value())?;
                }
            }
        }
        self.source
            .admit_grouped_strains_with_bias(&mapped[..f], &mapped_bias)
            .map_err(|_| {
                invalid("full imposed harmonic waveform leaves the original material domain")
            })?;
        let mut delta = [Dd::ZERO; 6];
        let mut target = history.strain;
        for j in 0..3 {
            let mut start = Dd::from_f64(bias[j]);
            let mut end = start;
            let mut change = Dd::ZERO;
            let mut magnitude =
                bias[j].abs() + history.strain[self.inside[j]].abs() / ENGINEERING[j];
            for k in 0..f {
                let [a, b] = grouped[k][j];
                let phase = self.phases[k];
                start = start + a * phase.start[0] - b * phase.start[1];
                end = end + a * phase.end[0] - b * phase.end[1];
                change = change + a * phase.increment[0] - b * phase.increment[1];
                magnitude += l1(grouped[k][j]);
            }
            if (start * Dd::from_f64(ENGINEERING[j]) - Dd::from_f64(history.strain[self.inside[j]]))
                .to_f64()
                .abs()
                > 128.0 * f64::EPSILON * ENGINEERING[j] * magnitude
            {
                return Err(invalid(
                    "in-plane harmonic strain is discontinuous at the accepted clock",
                ));
            }
            delta[self.inside[j]] = number(change * Dd::from_f64(ENGINEERING[j]))?;
            target[self.inside[j]] = number(end * Dd::from_f64(ENGINEERING[j]))?.to_f64();
        }
        let identity = [[(1.0, 0), (0.0, 0)]; LINES];
        let model = &self.control.model;
        let mut work = Sum::default();
        self.temporal.add_stress_work_with_bias(
            &model.c_eq,
            &mapped[..f],
            &identity[..f],
            &mapped_bias,
            &mut work,
        );
        let mut heat = Sum::default();
        let mut mode_changes = Vec::with_capacity(n);
        let mut absolute = mapped_bias.map(|v| v.hi.abs() + v.lo.abs());
        let mut excursion = history.strain.map(f64::abs);
        for k in 0..f {
            for i in 0..6 {
                absolute[i] += l1(mapped[k][i]);
                excursion[i] += (frequencies[k] * self.duration_s).min(2.0) * l1(mapped[k][i]);
            }
        }
        for mode in 0..n {
            poll(&mut cancelled)?;
            let mut initial = Sum::default();
            for row in 0..n {
                initial.add([
                    self.control.modes[row * n + mode],
                    self.control.branch_rate_roots[row / 6],
                    history.q[row / 6][row % 6],
                    1.0,
                ]);
            }
            let initial = number(initial.retained_value())?;
            let drive_row = core::array::from_fn(|j| {
                self.control.forcing[mode][self.inside[j]] * ENGINEERING[j]
            });
            let drive = rate_fields(&drive_row, &grouped[..f], frequencies)?;
            let port = rate_fields(&self.work_rows[mode], &grouped[..f], frequencies)?;
            let response = self.kernels[mode]
                .evaluate_cancellable(initial, &drive[..f], &port[..f], &mut cancelled)
                .map_err(moment_error)?;
            add(&mut heat, response.mean_square, self.duration_s);
            add(&mut work, response.mean_mixed, self.duration_s);
            mode_changes.push(response.increment);
            let tau = 1.0 / self.control.rates[mode];
            let forcing: f64 = drive[..f].iter().copied().map(l1).sum();
            let driven = self.duration_s.min(tau) * forcing;
            let loss = if self.duration_s >= tau {
                1.0
            } else {
                self.duration_s / tau
            };
            for i in 0..6 {
                absolute[i] += self.transverse[mode][i].abs() * (initial.abs().to_f64() + driven);
                excursion[i] +=
                    self.transverse[mode][i].abs() * (loss * initial.abs().to_f64() + driven);
            }
        }
        let maximum = core::array::from_fn(|i| absolute[i].min(excursion[i]));
        model.validate_mandel_strain(&maximum).map_err(|_| {
            invalid("transient plane-stress strain envelope exceeds the original material domain")
        })?;
        let mut dq = vec![[Dd::ZERO; 6]; history.q.len()];
        let mut q = history.q.clone();
        let mut rhs = equilibrium_stress(model, &delta).map(|v| -v);
        for row in 0..n {
            poll(&mut cancelled)?;
            let mut sum = Sum::default();
            for (mode, &change) in mode_changes.iter().enumerate() {
                for word in [change.hi, change.lo] {
                    sum.add_ratio(
                        &[self.control.modes[row * n + mode], word],
                        &[self.control.branch_rate_roots[row / 6]],
                    );
                }
            }
            dq[row / 6][row % 6] = number(sum.retained_value())?;
            q[row / 6][row % 6] =
                (Dd::from_f64(history.q[row / 6][row % 6]) + dq[row / 6][row % 6]).to_f64();
        }
        for (j, branch) in model.branches().iter().enumerate() {
            let gram = branch.resolved_gram_rows(&model.factors[j]);
            for i in 0..6 {
                for row in 0..6 {
                    rhs[i] = rhs[i] - Dd::from_f64(gram[row][i]) * dq[j][row];
                }
            }
        }
        let transverse = self.constraint.solve(&rhs)?;
        for &i in &self.constraint.axes {
            delta[i] = transverse[i];
            target[i] = (Dd::from_f64(history.strain[i]) + delta[i]).to_f64();
        }
        model
            .validate_mandel_strain(&target)
            .map_err(|_| invalid("plane-stress endpoint leaves original strain domain"))?;
        let (stress, scale) = self.control.source_stress(&target, &q)?;
        if self
            .constraint
            .axes
            .iter()
            .any(|&i| stress[i].abs() > tolerance * scale[i])
        {
            return Err(invalid(
                "causal plane-stress endpoint constraint is unresolved",
            ));
        }
        let eq0 = apply(&model.equilibrium_gram, &history.strain);
        let eq1 = apply(&model.equilibrium_gram, &target);
        let eq_delta = model
            .equilibrium_gram
            .map(|row| (0..6).fold(Dd::ZERO, |s, i| s + Dd::from_f64(row[i]) * delta[i]));
        let mut change = Dd::ZERO;
        let mut exchange = 0.0;
        let mut stored = 0.5 * dot(&eq1, &eq1);
        for (old, difference) in
            core::iter::once((&eq0, &eq_delta)).chain(history.q.iter().zip(&dq))
        {
            for i in 0..6 {
                let term =
                    (Dd::from_f64(old[i]) + Dd::from_f64(0.5) * difference[i]) * difference[i];
                change = change + term;
                exchange += term.abs().to_f64();
            }
        }
        for memory in &q {
            stored += 0.5 * dot(memory, memory);
        }
        let work = number(work.retained_value())?;
        let heat = number(heat.retained_value())?;
        let residual = number(work - change - heat)?.to_f64();
        let total_work = (Dd::from_f64(history.work) + work).to_f64();
        let total_heat = (Dd::from_f64(history.heat) + heat).to_f64();
        if q.iter().any(|v| !finite(v))
            || ![stored, exchange, total_work, total_heat]
                .iter()
                .all(|v| v.is_finite())
            || stored < 0.0
            || heat.lt(Dd::ZERO)
            || residual.abs() > tolerance * (work.abs().to_f64() + heat.to_f64() + exchange)
        {
            return Err(invalid(
                "constrained harmonic work/storage/heat balance is unresolved",
            ));
        }
        let material = TensorPronyHarmonicWindowPoint {
            start_stored_energy_j_m3: history.stored,
            end_stored_energy_j_m3: stored,
            stored_energy_change_j_m3: number(change)?.to_f64(),
            kelvin_heat_j_m3: 0.0,
            maxwell_heat_j_m3: heat.to_f64(),
            material_heat_j_m3: heat.to_f64(),
            stress_work_j_m3: work.to_f64(),
            energy_balance_residual_j_m3: residual,
        };
        poll(&mut cancelled)?;
        let mut next = history.clone();
        next.strain = target;
        next.q = q;
        next.stress = core::array::from_fn(|i| {
            if self.control.stress_components[i] {
                0.0
            } else {
                stress[i]
            }
        });
        next.stored = stored;
        next.work = total_work;
        next.heat = total_heat;
        Ok(TensorPlaneStressHarmonicPoint {
            history: next,
            material,
            strain_increment: tensor(&delta.map(Dd::to_f64)),
            strain_bounds: [tensor(&maximum.map(|v| -v)), tensor(&maximum)],
        })
    }
}

fn rate_fields(
    row: &[f64; 3],
    strains: &[PlaneParts],
    frequencies: &[f64],
) -> Result<[[Dd; 2]; LINES]> {
    let mut values = [[Dd::ZERO; 2]; LINES];
    for (i, (&omega, strain)) in frequencies.iter().zip(strains).enumerate() {
        let mut sum = [Sum::default(); 2];
        for j in 0..3 {
            for [re, im] in [
                [strain[j][0].hi, strain[j][1].hi],
                [strain[j][0].lo, strain[j][1].lo],
            ] {
                sum[0].add([row[j], omega, im, 1.0]);
                sum[1].add([-row[j], omega, re, 1.0]);
            }
        }
        values[i] = [
            number(sum[0].retained_value())?,
            number(sum[1].retained_value())?,
        ];
    }
    Ok(values)
}
