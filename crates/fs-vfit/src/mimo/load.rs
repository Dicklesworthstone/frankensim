//! Bounded admission and held-flow stepping for proper causal MIMO loads.
//!
//! This is the common admission used by waveguide and field propagators.
//! It retains the original MimoZoh coordinates and exact interval-mean
//! output map; it does not construct another rational integrator. Continuous
//! passivity and the actual sampled mean-output check remain separate reports.
//! The sampled check is GridOnly, not a coupled discrete-storage certificate.

use super::{MimoModel, MimoPassivityReport, MimoZoh, check_mimo_passivity};
use crate::model::PoleTerm;
use fs_math::c64::C64;

/// Source-bound continuous access to the original admitted pole realization.
pub mod continuous;
pub use continuous::MimoContinuous;

/// Cold resource and frequency limits for one proper rational load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MimoLoadLimits {
    /// Declared source/operator band [rad/s], covering the clock's Nyquist limit.
    pub band_rad_s: [f64; 2],
    /// Actual interval-mean ZOH passivity samples from DC through Nyquist.
    /// Must lie in 17..=65537; increasing this count does not prove passivity.
    pub sampled_passivity_points: usize,
    /// Maximum real state dimension, at most 256 for this dense owner.
    pub max_states: usize,
    /// Conservative retained matrix/state/report payload allowance [bytes].
    /// Excludes caller vectors, allocator overhead and cold eigensolver scratch.
    pub max_memory_bytes: usize,
}

/// A named admission or finite-step refusal; accepted pole states stay unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MimoLoadError(pub &'static str);

impl core::fmt::Display for MimoLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for MimoLoadError {}

/// An admitted source model, exact ZOH runtime and its two distinct reports.
///
/// Clone this owner when a larger interconnection needs a candidate transaction.
/// A refused step changes neither the accepted pole coordinates nor the caller's
/// output. Successful steps preserve the source model, clock and admission limits.
#[derive(Clone, Debug)]
pub struct PreparedMimoLoad {
    source: MimoModel,
    limits: MimoLoadLimits,
    time_step_s: f64,
    runtime: MimoZoh,
    continuous: MimoPassivityReport,
    sampled: MimoPassivityReport,
}

impl PreparedMimoLoad {
    /// Admit a finite stable proper impedance before allocating runtime states.
    ///
    /// Ports are bounded by 128, pole terms by 256, and real state dimension by
    /// the caller's at-most-256 budget. Exact matrix shapes, real-pole residues,
    /// conjugate-pair orientation and zero improper E are required. The complete
    /// source band must reach Nyquist at this clock; both continuous and actual
    /// sampled mean-output passivity checks must pass.
    pub fn new(
        model: &MimoModel, time_step_s: f64, limits: MimoLoadLimits,
    ) -> Result<Self, MimoLoadError> {
        if !time_step_s.is_finite() || time_step_s <= 0.0 {
            return Err(MimoLoadError("MIMO load requires a finite positive physical clock"));
        }
        let n = model.ports;
        let terms = model.poles.len();
        if !(1..=128).contains(&n)
            || terms > 256
            || limits.max_states > 256
            || !(17..=65_537).contains(&limits.sampled_passivity_points)
        {
            return Err(MimoLoadError("MIMO load exceeds port, pole or sampling limits"));
        }
        let states = model.order().checked_mul(n)
            .ok_or(MimoLoadError("MIMO state size overflow"))?;
        let entries = n * n;
        if states > limits.max_states
            || model.d.len() != entries
            || model.e.len() != entries
            || model.residues.len() != terms
            || model.residues.iter().any(|r| {
                r.len() != entries || r.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
            || model.d.iter().any(|v| !v.is_finite())
            || model.e.iter().any(|v| *v != 0.0)
            || !model.is_stable()
            || model.poles.iter().zip(&model.residues).any(|(p, r)| match p {
                PoleTerm::Real { pole, .. } => {
                    !pole.is_finite() || r.iter().any(|z| z.im != 0.0)
                }
                PoleTerm::Pair { pole, .. } => {
                    !pole.re.is_finite() || !pole.im.is_finite() || pole.im <= 0.0
                }
            })
        {
            return Err(MimoLoadError(
                "MIMO load requires a finite stable proper model with exact shapes",
            ));
        }
        let band = limits.band_rad_s;
        if !band.iter().all(|v| v.is_finite())
            || band[0] <= 0.0 || band[1] <= band[0]
            || band[1] < core::f64::consts::PI / time_step_s
        {
            return Err(MimoLoadError(
                "MIMO operator band must cover the physical clock's Nyquist limit",
            ));
        }
        // The original network owner's conservative retained-payload envelope.
        // Finite dense limits above keep this arithmetic representable on wasm32.
        let required = entries * (64 * terms + 64) + states * 64 + n * 256
            + (limits.sampled_passivity_points + 6400 + terms * 128) * 8;
        if required > limits.max_memory_bytes {
            return Err(MimoLoadError("MIMO load payload exceeds its separate memory allowance"));
        }
        let continuous = check_mimo_passivity(model, (band[0], band[1]))
            .map_err(|_| MimoLoadError("MIMO continuous passivity admission failed"))?;
        if !continuous.passive {
            return Err(MimoLoadError("MIMO continuous impedance has a passivity violation"));
        }
        let runtime = model.zoh(time_step_s)
            .map_err(|_| MimoLoadError("MIMO exact ZOH realization failed"))?;
        let sampled = runtime.mean_sampled_passivity(limits.sampled_passivity_points);
        if !sampled.passive {
            return Err(MimoLoadError("actual sampled MIMO impedance has a passivity violation"));
        }
        Ok(Self { source: model.clone(), limits, time_step_s, runtime, continuous, sampled })
    }

    /// Original source model, including port order and mutual matrix entries.
    #[must_use]
    pub fn source_model(&self) -> &MimoModel { &self.source }
    /// Limits actually admitted; no default band or sampling policy is inserted.
    #[must_use]
    pub fn limits(&self) -> &MimoLoadLimits { &self.limits }
    #[must_use]
    pub fn ports(&self) -> usize { self.source.ports }
    #[must_use]
    pub fn time_step_s(&self) -> f64 { self.time_step_s }
    /// Continuous and actual sampled mean-output reports, in that order.
    /// The second report retains its original GridOnly classification.
    #[must_use]
    pub fn reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        (&self.continuous, &self.sampled)
    }
    /// Row-major interval-mean direct matrix for a same-interval coupled solve.
    #[must_use]
    pub fn mean_direct_term(&self) -> &[f64] { self.runtime.mean_direct_term() }
    /// Accepted original ZOH pole coordinates, in term then input-port order.
    #[must_use]
    pub fn pole_states(&self) -> &[C64] { self.runtime.pole_states() }

    /// Write the accepted history's interval-mean effort, without advancing it.
    /// The output must have one entry per port; refusal leaves it unchanged.
    pub fn mean_history_output_into(&self, output: &mut [f64]) -> Result<(), MimoLoadError> {
        let n = self.ports();
        if output.len() != n {
            return Err(MimoLoadError("MIMO history output must have one entry per port"));
        }
        let mut candidate = [0.0; 128];
        self.runtime.mean_history_output_into(&mut candidate[..n]);
        if candidate[..n].iter().any(|v| !v.is_finite()) {
            return Err(MimoLoadError("MIMO history output left the finite set"));
        }
        output.copy_from_slice(&candidate[..n]);
        Ok(())
    }

    /// Advance the original exact ZOH states for one held input, writing the
    /// exact interval-mean effort. Both slices must contain one value per port.
    ///
    /// No allocation occurs. Invalid shapes, nonfinite input, output or candidate
    /// state refuse before accepted history or caller output changes. The caller
    /// remains responsible for checking its coupled work and accepting the whole
    /// interconnection transaction; clone this owner for that candidate state.
    pub fn step_mean_into(
        &mut self, input: &[f64], output: &mut [f64],
    ) -> Result<(), MimoLoadError> {
        let n = self.ports();
        if input.len() != n || output.len() != n {
            return Err(MimoLoadError("MIMO step requires one input and output per port"));
        }
        if input.iter().any(|v| !v.is_finite()) {
            return Err(MimoLoadError("MIMO held input must be finite"));
        }
        let mut candidate = [0.0; 128];
        self.runtime.mean_history_output_into(&mut candidate[..n]);
        for (i, effort) in candidate[..n].iter_mut().enumerate() {
            *effort += self.runtime.mean_direct_term()[i * n..(i + 1) * n]
                .iter().zip(input).map(|(d, u)| d * u).sum::<f64>();
        }
        if candidate[..n].iter().any(|v| !v.is_finite()) || !self.runtime.preview_state(input) {
            return Err(MimoLoadError("MIMO candidate output or state left the finite set"));
        }
        self.runtime.accept_state(input);
        output.copy_from_slice(&candidate[..n]);
        Ok(())
    }

    /// Transfer the same admitted runtime and reports to an existing owner that
    /// already retains the original source specification and manages preview/accept.
    #[must_use]
    pub(crate) fn into_parts(self) -> (MimoZoh, MimoPassivityReport, MimoPassivityReport) {
        (self.runtime, self.continuous, self.sampled)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::passivity::CertificateClass;

    const DT: f64 = 1.0 / 48_000.0;
    fn limits() -> MimoLoadLimits {
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 33, max_states: 4, max_memory_bytes: 1 << 20,
        }
    }
    fn model() -> MimoModel {
        MimoModel {
            ports: 2,
            poles: vec![PoleTerm::Real { pole: -400.0, residue: 0.0 }],
            residues: vec![vec![
                C64::from_re(200.0), C64::from_re(50.0),
                C64::from_re(50.0), C64::from_re(180.0),
            ]],
            d: vec![2.0, 0.5, 0.5, 3.0], e: vec![0.0; 4], reciprocal: true,
        }
    }

    #[test]
    fn prepared_load_preserves_original_mean_output_states_and_reports() {
        let source = model();
        let mut prepared = PreparedMimoLoad::new(&source, DT, limits()).unwrap();
        let mut original = source.zoh(DT).unwrap();
        assert_eq!(prepared.source_model(), &source);
        assert_eq!(prepared.limits(), &limits());
        assert_eq!(prepared.time_step_s().to_bits(), DT.to_bits());
        assert!(prepared.reports().0.passive);
        assert!(prepared.reports().1.passive);
        assert_eq!(prepared.reports().1.class, CertificateClass::GridOnly);
        for input in [[1.0, -0.5], [0.0, 0.75], [-0.3, 0.0], [0.0, 0.0]] {
            let mut history = [0.0; 2];
            prepared.mean_history_output_into(&mut history).unwrap();
            let expected = original.step_mean(&input);
            let mut actual = [0.0; 2];
            prepared.step_mean_into(&input, &mut actual).unwrap();
            for i in 0..2 {
                let solve_pressure = history[i]
                    + prepared.mean_direct_term()[i * 2..(i + 1) * 2]
                        .iter().zip(&input).map(|(d, u)| d * u).sum::<f64>();
                assert_eq!(actual[i], solve_pressure);
                assert!((actual[i] - expected[i]).abs()
                    <= 16.0 * f64::EPSILON * (actual[i].abs() + expected[i].abs()));
            }
            assert_eq!(prepared.pole_states(), original.pole_states());
        }
    }

    #[test]
    fn malformed_raw_loads_refuse_before_realization() {
        let mut cases = Vec::new();
        let mut bad = model();
        bad.d.pop();
        cases.push(bad);
        let mut bad = model();
        bad.e[0] = 1.0;
        cases.push(bad);
        let mut bad = model();
        bad.residues[0][0].im = 1.0;
        cases.push(bad);
        let mut bad = model();
        bad.residues[0][0].re = f64::NAN;
        cases.push(bad);
        let mut bad = model();
        bad.poles[0] = PoleTerm::Real { pole: 1.0, residue: 0.0 };
        cases.push(bad);
        let mut bad = model();
        bad.poles[0] = PoleTerm::Pair {
            pole: C64::new(-1.0, -1.0), residue: C64::ZERO,
        };
        cases.push(bad);
        let mut bad = model();
        bad.ports = usize::MAX;
        cases.push(bad);
        for bad in cases {
            assert!(PreparedMimoLoad::new(&bad, DT, limits()).is_err());
        }
        for dt in [0.0, -DT, f64::NAN, f64::INFINITY] {
            assert!(PreparedMimoLoad::new(&model(), dt, limits()).is_err());
        }
        let mut bad = limits();
        bad.band_rad_s[1] *= 1.0 - 1e-9;
        assert!(PreparedMimoLoad::new(&model(), DT, bad).is_err());
        let mut bad = limits();
        bad.sampled_passivity_points = 16;
        assert!(PreparedMimoLoad::new(&model(), DT, bad).is_err());
        let mut bad = limits();
        bad.max_states = 1;
        assert!(PreparedMimoLoad::new(&model(), DT, bad).is_err());
        let mut bad = limits();
        bad.max_memory_bytes = 0;
        assert!(PreparedMimoLoad::new(&model(), DT, bad).is_err());
        let mut active = model();
        active.d = vec![-2.0, 0.0, 0.0, -3.0];
        assert!(PreparedMimoLoad::new(&active, DT, limits()).is_err());
    }

    #[test]
    fn finite_effort_cannot_hide_overflowing_candidate_coordinates() {
        let source = MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real { pole: -1.0e-4, residue: 0.0 }],
            residues: vec![vec![C64::from_re(1.0e-300)]],
            d: vec![0.0], e: vec![0.0], reciprocal: true,
        };
        let dt = 1.0e4;
        let limits = MimoLoadLimits {
            band_rad_s: [1.0e-9, core::f64::consts::PI / dt],
            sampled_passivity_points: 17, max_states: 1, max_memory_bytes: 1 << 20,
        };
        let mut load = PreparedMimoLoad::new(&source, dt, limits).unwrap();
        let input = 1.0e307;
        assert!((load.mean_direct_term()[0] * input).is_finite());
        let before = load.pole_states().to_vec();
        let mut output = [42.0];
        assert!(load.step_mean_into(&[input], &mut output).is_err());
        assert_eq!(output, [42.0]);
        assert_eq!(load.pole_states(), before);
        load.step_mean_into(&[1.0], &mut output).unwrap();
        assert!(load.pole_states()[0].re.is_finite());
    }

    #[test]
    fn failed_steps_and_abandoned_clones_do_not_advance_accepted_history() {
        let mut load = PreparedMimoLoad::new(&model(), DT, limits()).unwrap();
        let mut output = [0.0; 2];
        load.step_mean_into(&[0.3, -0.2], &mut output).unwrap();
        let before = load.pole_states().to_vec();
        let mut reference = load.clone();
        let mut candidate = load.clone();
        candidate.step_mean_into(&[0.7, 0.2], &mut output).unwrap();
        assert_ne!(candidate.pole_states(), load.pole_states());
        assert_eq!(load.pole_states(), before);
        let mut unchanged = [42.0; 2];
        assert!(load.step_mean_into(&[f64::NAN, 0.0], &mut unchanged).is_err());
        assert!(load.step_mean_into(&[f64::MAX, f64::MAX], &mut unchanged).is_err());
        assert!(load.step_mean_into(&[1.0], &mut unchanged).is_err());
        assert!(load.step_mean_into(&[1.0, 0.0], &mut unchanged[..1]).is_err());
        assert!(load.mean_history_output_into(&mut unchanged[..1]).is_err());
        assert_eq!(unchanged, [42.0; 2]);
        assert_eq!(load.pole_states(), before);
        let mut expected = [0.0; 2];
        reference.step_mean_into(&[-0.1, 0.4], &mut expected).unwrap();
        load.step_mean_into(&[-0.1, 0.4], &mut output).unwrap();
        assert_eq!(output, expected);
        assert_eq!(load.pole_states(), reference.pole_states());
    }
}
