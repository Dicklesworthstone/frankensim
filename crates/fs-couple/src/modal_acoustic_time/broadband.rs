//! One-way broadband receiver observation. The lower rows of the augmented
//! [mode, continuous filter, held force] exponential integrate the filter under
//! the actual within-sample modal trajectory. The original cached modal owner
//! advances the mechanical states, retaining its exact energy/work bits.
use super::{ModalAcousticFrame, ModalAcousticTimeModel, ModalAcousticWorkspace};
pub use crate::broadband_radiation::modal::{
    BroadbandModalControls, BroadbandModalError, BroadbandModalReport, ModalTransferSamples,
};
use crate::broadband_radiation::modal::{admit_modal_samples, fit_modal_velocity_transfers};
use fs_math::{c64::C64, det};

/// Maximum lag in a receiver's explicit sampled delay history.
pub const MAX_MODAL_OBSERVER_DELAY_SAMPLES: usize = 1_048_576;
/// Aggregate dense exponential-entry budget across every retained mode.
pub const MAX_MODAL_AUGMENTED_MATRIX_ENTRIES: usize = 4_194_304;

struct ExactFilter {
    n: usize,
    transition: Vec<f64>,
    state: Vec<f64>,
    candidate: Vec<f64>,
    c: Vec<f64>,
    d: f64,
    q_scale: f64,
    model: fs_vfit::RationalModel,
}
struct Delay {
    taps: Vec<(usize, f64)>,
    history: Vec<f64>,
    cursor: usize,
}
impl Delay {
    fn new(samples: f64, half_len: usize) -> Result<Self, BroadbandModalError> {
        if !samples.is_finite()
            || samples < 0.0
            || samples > MAX_MODAL_OBSERVER_DELAY_SAMPLES as f64
        {
            return Err(BroadbandModalError::Invalid(
                "observer propagation delay exceeds its sample envelope",
            ));
        }
        let nearest = samples.round();
        let integer = (samples - nearest).abs() <= 4.0 * f64::EPSILON * samples.max(1.0);
        let n = if integer {
            nearest as usize
        } else {
            samples.floor() as usize
        };
        let mut taps = Vec::new();
        if integer {
            taps.push((n, 1.0));
        } else {
            if n < half_len {
                return Err(BroadbandModalError::Invalid(
                    "physical delay is too short for the declared centred fractional kernel; raise the observation clock or refine with a shorter kernel",
                ));
            }
            let fraction = samples - n as f64;
            let l = half_len as isize;
            for j in -l..=l {
                let x = j as f64 - fraction;
                let u = x / (half_len + 1) as f64;
                let sinc = if x == 0.0 {
                    1.0
                } else {
                    det::sin(core::f64::consts::PI * x) / (core::f64::consts::PI * x)
                };
                let window = 0.42
                    + 0.5 * det::cos(core::f64::consts::PI * u)
                    + 0.08 * det::cos(core::f64::consts::TAU * u);
                taps.push(((n as isize + j) as usize, sinc * window));
            }
            let sum = taps.iter().map(|(_, v)| v).sum::<f64>();
            if !sum.is_finite() || sum == 0.0 {
                return Err(BroadbandModalError::NonFinite {
                    what: "fractional-delay normalization",
                });
            }
            for (_, v) in &mut taps {
                *v /= sum;
            }
        }
        let last = taps
            .iter()
            .map(|(lag, _)| *lag)
            .max()
            .ok_or(BroadbandModalError::Invalid("empty delay"))?;
        if last > MAX_MODAL_OBSERVER_DELAY_SAMPLES {
            return Err(BroadbandModalError::Cap {
                what: "observer delay history",
                found: last,
                cap: MAX_MODAL_OBSERVER_DELAY_SAMPLES,
            });
        }
        Ok(Self {
            taps,
            history: vec![0.0; last + 1],
            cursor: 0,
        })
    }
    fn response(&self, w: f64, dt: f64) -> C64 {
        self.taps.iter().fold(C64::ZERO, |sum, &(lag, c)| {
            let phase = w * dt * lag as f64;
            sum + C64::new(det::cos(phase), det::sin(phase)).scale(c)
        })
    }
    fn preview(&self, new: f64) -> Result<f64, BroadbandModalError> {
        let len = self.history.len();
        let mut y = 0.0;
        for &(lag, c) in &self.taps {
            y += c * if lag == 0 {
                new
            } else {
                self.history[(self.cursor + len - lag) % len]
            };
        }
        if !y.is_finite() {
            return Err(BroadbandModalError::NonFinite {
                what: "delayed receiver pressure",
            });
        }
        Ok(y)
    }
    fn commit(&mut self, new: f64) {
        self.history[self.cursor] = new;
        self.cursor = (self.cursor + 1) % self.history.len();
    }
}

/// One-way exact held-force receiver; the original mechanical owner is retained.
pub struct BroadbandModalObserver {
    model: ModalAcousticTimeModel,
    workspace: ModalAcousticWorkspace,
    filters: Vec<ExactFilter>,
    constant: Option<Vec<C64>>,
    delay: Delay,
    report: BroadbandModalReport,
}
impl BroadbandModalObserver {
    /// Fit/gate the complete declared source and prepare exact augmented ZOH.
    /// Initial acoustic state and delay history are zero. Every mode is retained.
    ///
    /// # Errors
    /// Refuses malformed/over-cap inputs, unstable/improper fits, held-out
    /// budget violations, noncausal centered delay taps and exponential refusal.
    pub fn build(
        model: &ModalAcousticTimeModel,
        s: &ModalTransferSamples,
        c: BroadbandModalControls,
        rate: u32,
    ) -> Result<(Self, BroadbandModalReport), BroadbandModalError> {
        admit_modal_samples(s.per_mode.len(), s.omega_rad_s.len(), s.held_out.len())?;
        if rate == 0
            || model.sample_period_s().to_bits() != (1.0 / f64::from(rate)).to_bits()
            || model.modes().len() != s.per_mode.len()
        {
            return Err(BroadbandModalError::Invalid(
                "observer model, source basis and declared sample clock must match",
            ));
        }
        let bank = fit_modal_velocity_transfers(s, c, rate)?;
        let dt = model.sample_period_s();
        let delay = Delay::new(bank.report.delay_samples, c.delay_kernel_half_len)?;
        let mut report = bank.report;
        if bank.constant.is_none() {
            let mut maximum = 0.0_f64;
            for (w, row) in &s.held_out {
                let d = delay.response(*w, dt);
                for (fit, &reference) in bank.models.iter().zip(row) {
                    let actual = fit.eval_iw(*w).conj() * d;
                    let error = if reference == C64::ZERO {
                        if actual == C64::ZERO {
                            0.0
                        } else {
                            f64::INFINITY
                        }
                    } else {
                        (actual - reference).abs() / reference.abs()
                    };
                    if !error.is_finite() {
                        return Err(BroadbandModalError::HeldOut {
                            max_rel: error,
                            budget: c.max_held_out_rel_error,
                        });
                    }
                    maximum = maximum.max(error);
                }
            }
            report.held_out_max_rel_error = maximum;
            if maximum > c.max_held_out_rel_error {
                return Err(BroadbandModalError::HeldOut {
                    max_rel: maximum,
                    budget: c.max_held_out_rel_error,
                });
            }
        }
        let entries = bank
            .models
            .iter()
            .try_fold(0usize, |sum, m| {
                let n = m.order() + 3;
                sum.checked_add(n * n)
            })
            .ok_or(BroadbandModalError::Invalid(
                "augmented observer shape overflow",
            ))?;
        if entries > MAX_MODAL_AUGMENTED_MATRIX_ENTRIES {
            return Err(BroadbandModalError::Cap {
                what: "augmented observer matrix entries",
                found: entries,
                cap: MAX_MODAL_AUGMENTED_MATRIX_ENTRIES,
            });
        }
        let mut filters = Vec::with_capacity(bank.models.len());
        for (mode, fit) in model.modes().iter().zip(bank.models) {
            let ss = fit.state_space();
            let n = ss.n + 3;
            let mut a = vec![0.0; n * n];
            let w = mode.angular_frequency_rad_s;
            let q_scale = if w == 0.0 { 1.0 } else { w };
            a[1] = q_scale * dt;
            a[n] = -w * w / q_scale * dt;
            a[n + 1] = -2.0 * mode.damping_ratio * w * dt;
            a[n + n - 1] = dt;
            for i in 0..ss.n {
                a[(i + 2) * n + 1] = ss.b[i] * dt;
                for j in 0..ss.n {
                    a[(i + 2) * n + j + 2] = ss.a[i * ss.n + j] * dt;
                }
            }
            let transition =
                fs_phs::linear_exact::exponential(&a, n).map_err(BroadbandModalError::Exact)?;
            filters.push(ExactFilter {
                n,
                transition,
                state: vec![0.0; ss.n],
                candidate: vec![0.0; ss.n],
                c: ss.c,
                d: ss.d,
                q_scale,
                model: fit,
            });
        }
        let mut mechanics = model.clone();
        // The new receiver owns the pressure ceiling. Retaining an unrelated
        // old narrow-band receiver here would admit/refuse on a second pressure.
        // This read-only transfer edit changes no mechanical coefficient/state.
        for mode in &mut mechanics.modes {
            mode.pressure_per_modal_velocity = C64::ZERO;
        }
        let observer = Self {
            model: mechanics,
            workspace: ModalAcousticWorkspace::new(model),
            filters,
            constant: bank.constant,
            delay,
            report: report.clone(),
        };
        Ok((observer, report))
    }
    pub const fn model(&self) -> &ModalAcousticTimeModel {
        &self.model
    }
    pub const fn report(&self) -> &BroadbandModalReport {
        &self.report
    }
    pub fn last_modal_frame(&self) -> &ModalAcousticFrame {
        &self.workspace.frame
    }
    /// Physical exp(-i*w*t) fitted operator INCLUDING the sampled delay kernel.
    /// This is a sampled held-out diagnostic, not a continuum error bound or
    /// the force-to-pressure transfer of the held mechanical oscillator.
    pub fn transfer_at(&self, mode: usize, omega: f64) -> Result<C64, BroadbandModalError> {
        if mode >= self.model.modes().len()
            || !omega.is_finite()
            || omega <= 0.0
            || omega * self.model.sample_period_s() >= core::f64::consts::PI
        {
            return Err(BroadbandModalError::Invalid(
                "operator query needs an existing mode and sub-Nyquist positive frequency",
            ));
        }
        let h = if let Some(constant) = &self.constant {
            constant[mode]
        } else {
            self.filters[mode].model.eval_iw(omega).conj()
        };
        Ok(h * self.delay.response(omega, self.model.sample_period_s()))
    }
    /// Exact-ZOH step under the supplied held modal forces. Filters and delay
    /// commit only after all pressure, force and mechanical budget checks pass.
    /// Initial acoustic filter/delay history is explicitly zero.
    pub fn step(&mut self, force: &[f64]) -> Result<f64, BroadbandModalError> {
        if force.len() != self.model.modes().len() || force.iter().any(|v| !v.is_finite()) {
            return Err(BroadbandModalError::Invalid(
                "finite held force required for every retained mode",
            ));
        }
        let mut raw = 0.0;
        if let Some(h) = &self.constant {
            for (i, (&transfer, &f)) in h.iter().zip(force).enumerate() {
                let state = self
                    .model
                    .advance_exact_zoh_step(i, self.model.states()[i], f);
                let w = self.model.modes()[i].angular_frequency_rad_s;
                if w == 0.0 {
                    if transfer != C64::ZERO {
                        return Err(BroadbandModalError::Invalid(
                            "free coordinates have no constant narrow-band quadrature",
                        ));
                    }
                } else {
                    let q = state.displacement_m_sqrt_kg - f / (w * w);
                    raw += transfer.re * state.velocity_m_sqrt_kg_per_s + transfer.im * w * q;
                }
            }
        } else {
            for (i, (filter, &f)) in self.filters.iter_mut().zip(force).enumerate() {
                let old = self.model.states()[i];
                let next = self.model.advance_exact_zoh_step(i, old, f);
                for row in 0..filter.state.len() {
                    let r = (row + 2) * filter.n;
                    let mut x = filter.transition[r] * filter.q_scale * old.displacement_m_sqrt_kg
                        + filter.transition[r + 1] * old.velocity_m_sqrt_kg_per_s;
                    for (j, &value) in filter.state.iter().enumerate() {
                        x += filter.transition[r + j + 2] * value;
                    }
                    x += filter.transition[r + filter.n - 1] * f;
                    if !x.is_finite() {
                        return Err(BroadbandModalError::NonFinite {
                            what: "augmented filter state",
                        });
                    }
                    filter.candidate[row] = x;
                }
                let mut y = filter.d * next.velocity_m_sqrt_kg_per_s;
                for (&c, &x) in filter.c.iter().zip(&filter.candidate) {
                    y += c * x;
                }
                raw += y;
            }
        }
        if !raw.is_finite() {
            return Err(BroadbandModalError::NonFinite {
                what: "modal pressure superposition",
            });
        }
        let pressure = self.delay.preview(raw)?;
        if pressure.abs() > self.model.budget.maximum_abs_pressure_pa {
            return Err(BroadbandModalError::Invalid(
                "broadband observer exceeds the declared pressure ceiling",
            ));
        }
        self.model
            .step_into(force, &mut self.workspace)
            .map_err(BroadbandModalError::Modal)?;
        for filter in &mut self.filters {
            filter.state.copy_from_slice(&filter.candidate);
        }
        self.delay.commit(raw);
        Ok(pressure)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn integer_delay_copies_bits_and_history_cap_refuses_plus_one() {
        let mut delay = Delay::new(7.0, 1).unwrap();
        let mut values = Vec::new();
        for n in 0..73 {
            let value = det::cos(n as f64 * 0.137);
            let expected = if n < 7 { 0.0_f64 } else { values[n - 7] };
            assert_eq!(delay.preview(value).unwrap().to_bits(), expected.to_bits());
            delay.commit(value);
            values.push(value);
        }
        assert!(Delay::new(MAX_MODAL_OBSERVER_DELAY_SAMPLES as f64, 1).is_ok());
        assert!(Delay::new((MAX_MODAL_OBSERVER_DELAY_SAMPLES + 1) as f64, 1).is_err());
        assert!(Delay::new(1.25, 8).is_err());
    }
}
