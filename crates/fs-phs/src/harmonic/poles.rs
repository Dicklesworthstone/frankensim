//! Isolated finite modes of the complete continuous pHS descriptor pencil.
//! A supplied seed selects one local pole; no spectral completeness or transfer
//! observability claim is inferred from a single mode.

use super::{C64, HarmonicModel, finite};
use fs_la::nonlinear_eig::{
    NonlinearEigError, NonlinearPole, PoleEvaluation, newton_pole_with_cancel,
};
use fs_math::det;

pub use fs_la::nonlinear_eig::PoleOptions;

mod spectrum;
pub use spectrum::{BeynOptions, HarmonicSpectrum, HarmonicSpectrumError, ResonanceRegion};

/// Physical port participation reuses the original pencil and mode observer.
pub mod response;

/// A local mode of exp(s t), including every original constraint and memory state.
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicPole {
    /// Laplace pole s [1/s]. Re(s)<0 means decay; no unstable pole is clamped.
    pub laplace_pole_per_s: C64,
    /// omega=i*s [rad/s] under exp(-i omega t). Decay has Im(omega)<0.
    pub angular_frequency_rad_s: C64,
    /// Original differential state followed by algebraic efforts, unit numerical
    /// 2-norm. This is a coordinate normalization, not unit physical energy.
    pub state: Vec<C64>,
    /// Original storage gradient and algebraic efforts.
    pub effort: Vec<C64>,
    /// Original G^T effort in port order. A zero output can be an unobservable
    /// internal mode; finite pencil poles are not automatically transfer poles.
    pub port_output: Vec<C64>,
    /// x* Q x / 2 for the returned amplitude; constraints store no energy.
    pub stored_energy_j: f64,
    /// Re(e* R e + e* T x) for the same amplitude. T is present only for
    /// explicit nonlinear equilibrium models. Then this is signed incremental
    /// work, not total biased-source heat; negative values are not clipped.
    pub dissipated_power_w: f64,
    /// 2 Re(s) H + D, independently evaluated from physical storage and loss.
    pub decay_power_residual_w: f64,
    /// Maximum componentwise residual of the ORIGINAL split pencil sE-A.
    pub relative_residual: f64,
    pub iterations: usize,
}

/// Pole-specific failures leave the model and any reusable harmonic LU intact.
#[derive(Debug, Clone, PartialEq)]
pub enum HarmonicPoleError {
    /// Local nonlinear iteration, data, capacity or cancellation refusal.
    Solve(NonlinearEigError),
    /// The original descriptor equations failed their existing model tolerance.
    Residual { relative: f64, tolerance: f64 },
    /// A physical observation or conditioning scale is not representable.
    Invalid { what: &'static str },
}
impl From<NonlinearEigError> for HarmonicPoleError {
    fn from(error: NonlinearEigError) -> Self {
        Self::Solve(error)
    }
}
impl core::fmt::Display for HarmonicPoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "harmonic pole refused: {self:?}")
    }
}
impl std::error::Error for HarmonicPoleError {}

fn poll(cancelled: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), NonlinearEigError> {
    if cancelled() {
        Err(NonlinearEigError::Cancelled)
    } else {
        Ok(())
    }
}
fn divide(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}

// Local and contour solves share precisely this original descriptor pencil.
// Preparation stays lazy so each numerical owner admits its resources first.
struct PolePencil<'a> {
    model: &'a HarmonicModel,
    spectral_scale: f64,
    conditioning: Option<(Vec<f64>, f64, f64)>,
}

impl<'a> PolePencil<'a> {
    fn new(model: &'a HarmonicModel, spectral_scale: f64) -> Self {
        Self {
            model,
            spectral_scale,
            conditioning: None,
        }
    }

    fn evaluate(
        &mut self,
        s: C64,
        cancel: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, NonlinearEigError> {
        let model = self.model;
        let n = model.n;
        if self.conditioning.is_none() {
            let mut rows = Vec::with_capacity(n);
            let mut generator_sq = 0.0;
            let mut descriptor_norm = 0.0;
            for i in 0..n {
                poll(cancel)?;
                let floor = if i < model.n_diff {
                    self.spectral_scale
                } else {
                    0.0
                };
                let scale = model.generator[i * n..(i + 1) * n]
                    .iter()
                    .map(|v| v.abs())
                    .fold(floor, f64::max);
                if !scale.is_finite() || scale <= 0.0 {
                    return Err(NonlinearEigError::BadInput {
                        what: "a zero descriptor row cannot define an isolated pole",
                    });
                }
                for &a in &model.generator[i * n..(i + 1) * n] {
                    generator_sq += (a / scale) * (a / scale);
                }
                if i < model.n_diff {
                    descriptor_norm = det::hypot(descriptor_norm, 1.0 / scale);
                }
                rows.push(scale);
            }
            self.conditioning = Some((rows, det::sqrt(generator_sq), descriptor_norm));
        }
        let (rows, a_norm, e_norm) = self.conditioning.as_ref().expect("prepared above");
        let residual_scale = *a_norm + s.abs().max(self.spectral_scale) * *e_norm;
        if !residual_scale.is_finite() || residual_scale <= 0.0 {
            return Err(NonlinearEigError::FunctionRefused { at: s });
        }
        let mut matrix = vec![C64::ZERO; n * n];
        let mut derivative = vec![C64::ZERO; n * n];
        for i in 0..n {
            poll(cancel)?;
            for j in 0..n {
                matrix[i * n + j] = C64::from_re(-model.generator[i * n + j] / rows[i]);
            }
            if i < model.n_diff {
                matrix[i * n + i] = matrix[i * n + i] + divide(s, rows[i]);
                derivative[i * n + i] = C64::from_re(1.0 / rows[i]);
            }
        }
        Ok(PoleEvaluation {
            matrix,
            derivative,
            residual_scale,
        })
    }
}

impl HarmonicModel {
    /// Refine one isolated finite free mode of the ORIGINAL continuous model.
    ///
    /// The pencil is s*diag(I_diff,0) - (J-R)*diag(Q_diff,I_alg), plus the
    /// retained nonlinear resisting tangent on explicit equilibrium images.
    /// All Maxwell states, mutual storage and algebraic constraints remain in
    /// the solve; the singular descriptor matrix is never inverted. A nearby
    /// harmonic response's complete state is a useful seed for observable modes.
    ///
    /// Constant row equilibration is prepared once using the declared spectral
    /// scale, so the derivative remains analytic and the eigenvector keeps its
    /// original coordinates. Acceptance also checks the unscaled split pencil
    /// with this model's existing componentwise residual tolerance.
    ///
    /// This is not a spectrum search, nonlinear limit cycle, complete mode count,
    /// or certificate of a transfer-function pole. Infinite descriptor roots,
    /// unresolved repeated/defective roots and a singular pencil can refuse.
    /// The separate options bound the pole workspace/work; the already prepared
    /// model retains its original allocation. No accepted history is changed.
    pub fn resonance_pole(
        &self,
        laplace_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicPole, HarmonicPoleError> {
        // The numerical owner admits sizes/resources BEFORE this lazy callback
        // allocates or assembles any additional matrices.
        let mut pencil = PolePencil::new(self, options.spectral_scale);
        let n = self.n;
        let pole = newton_pole_with_cancel(
            n,
            laplace_seed_per_s,
            state_seed,
            options,
            |s, cancel| pencil.evaluate(s, cancel),
            &mut cancelled,
        )?;
        self.observe_pole(pole, &mut cancelled)
    }

    fn observe_pole(
        &self,
        pole: NonlinearPole,
        mut cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicPole, HarmonicPoleError> {
        let n = self.n;
        let s = pole.value;
        let state = pole.vector;
        let mut relative = 0.0_f64;
        for i in 0..n {
            poll(&mut cancelled)?;
            let scale = self.generator[i * n..(i + 1) * n]
                .iter()
                .map(|v| v.abs())
                .fold(if i < self.n_diff { s.abs() } else { 0.0 }, f64::max);
            if scale == 0.0 {
                continue;
            }
            let diagonal = if i < self.n_diff {
                divide(s, scale)
            } else {
                C64::ZERO
            };
            let mut residual = diagonal * state[i];
            let mut denominator = diagonal.abs() * state[i].abs();
            for (j, &x) in state.iter().enumerate() {
                let a = self.generator[i * n + j] / scale;
                residual = residual - x.scale(a);
                denominator += a.abs() * x.abs();
            }
            let row = if denominator == 0.0 {
                if residual == C64::ZERO {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                residual.abs() / denominator
            };
            relative = relative.max(row);
        }
        if !relative.is_finite() || relative > self.budget.relative_residual_tolerance {
            return Err(HarmonicPoleError::Residual {
                relative,
                tolerance: self.budget.relative_residual_tolerance,
            });
        }
        let mut effort = vec![C64::ZERO; n];
        for i in 0..n {
            poll(&mut cancelled)?;
            if i < self.n_diff {
                for (j, &x) in state.iter().enumerate().take(self.n_diff) {
                    effort[i] = effort[i] + x.scale(self.q[i * n + j]);
                }
            } else {
                effort[i] = state[i];
            }
        }
        let stored_energy_j = 0.5
            * (0..self.n_diff)
                .map(|i| (state[i].conj() * effort[i]).re)
                .sum::<f64>();
        let mut dissipated_power_w = 0.0;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut row = C64::ZERO;
            for (j, &e) in effort.iter().enumerate() {
                row = row + e.scale(self.r[i * n + j]);
            }
            if self.resisting_tangent.is_some() {
                row = row + self.nonlinear_loss_row(i, &state);
            }
            dissipated_power_w += (effort[i].conj() * row).re;
        }
        let mut port_output = vec![C64::ZERO; self.m];
        for (p, output) in port_output.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            for (i, &e) in effort.iter().enumerate() {
                *output = *output + e.scale(self.g[i * self.m + p]);
            }
        }
        let decay_power_residual_w = 2.0 * s.re * stored_energy_j + dissipated_power_w;
        if effort.iter().chain(&port_output).any(|&z| !finite(z))
            || [stored_energy_j, dissipated_power_w, decay_power_residual_w]
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(HarmonicPoleError::Invalid {
                what: "physical pole storage, loss or port output is not representable",
            });
        }
        poll(&mut cancelled)?;
        Ok(HarmonicPole {
            laplace_pole_per_s: s,
            angular_frequency_rad_s: C64::new(-s.im, s.re),
            state,
            effort,
            port_output,
            stored_energy_j,
            dissipated_power_w,
            decay_power_residual_w,
            relative_residual: relative,
            iterations: pole.iterations,
        })
    }
}
