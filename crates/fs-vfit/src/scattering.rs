//! Same-sample FIR scattering ports for implicit physical junctions.
//!
//! A sampled reflection is `b[n] = r[0] a[n] + sum(k>=1, r[k] a[n-k])`.
//! Feeding a junction `DelayedFilter::incoming()` after the previous push
//! instead evaluates `b[n-1]`, adding a sample of phase and delaying r[0].
//! Here the existing FIR owner evaluates the one-sample-advanced STRICT tail;
//! the instantaneous coefficient is kept in the junction's algebraic equation.
//! No new convolution, fit, oscillator, gain correction or constitutive law.
//!
//! Supplying a kernel is not a passivity certificate. A convex crossfade of
//! frozen responses is not a proof of a physically moving or passive boundary.
//! All kernels use one characteristic-wave normalization supplied by the caller.

use crate::discretize::{DelayedFilter, DiscretizeError};

/// Refusal of a same-sample port or its bounded bank.
#[derive(Debug, Clone, PartialEq)]
pub enum ScatteringError {
    /// Invalid geometry-independent port input or work/storage bound.
    Invalid(&'static str),
    /// The shared discrete filter refused a coefficient, history or sample.
    Discrete(DiscretizeError),
    /// A bank participant refused after another participant advanced.
    /// Reconstruct the bank; do not continue a partially accepted sample.
    Poisoned,
}
impl core::fmt::Display for ScatteringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "FIR scattering port: {what}"),
            Self::Discrete(error) => write!(f, "FIR scattering port: {error}"),
            Self::Poisoned => write!(f, "FIR scattering bank has a partially accepted sample"),
        }
    }
}
impl std::error::Error for ScatteringError {}
impl From<DiscretizeError> for ScatteringError {
    fn from(error: DiscretizeError) -> Self { Self::Discrete(error) }
}

/// Immutable affine relation for the NEXT accepted sample, not the last output.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScatteringResponse {
    /// Dimensionless coefficient of the as-yet-unknown outgoing wave.
    pub direct: f64,
    /// Incoming contribution from already accepted outgoing samples.
    /// Has the caller's wave units, e.g. Pa for pressure waves.
    pub history: f64,
}
impl ScatteringResponse {
    /// Evaluate without changing any state. The source's implicit solve must
    /// determine the outgoing wave before either participant commits the sample.
    pub fn incoming(self, outgoing: f64) -> Result<f64, ScatteringError> {
        let value = self.direct * outgoing + self.history;
        if !outgoing.is_finite() || !value.is_finite() {
            return Err(ScatteringError::Invalid("incoming/outgoing wave is not finite"));
        }
        Ok(value)
    }
}

/// Exact FIR port with the direct reflection exposed to a physical junction.
/// The tail and its history are stepped only by `DelayedFilter`.
#[derive(Debug, Clone)]
pub struct FirScatteringPort {
    direct: f64,
    tail: DelayedFilter,
}
impl FirScatteringPort {
    /// Preserve every coefficient, including nonzero r[0]. At least four taps
    /// are required by the shared FIR owner. No passivity repair is performed.
    pub fn new(sample_period_s: f64, impulse: Vec<f64>) -> Result<Self, ScatteringError> {
        if impulse.len() < 4 || impulse.iter().any(|x| !x.is_finite()) {
            return Err(ScatteringError::Invalid("a finite reflection kernel needs at least four taps"));
        }
        let direct = impulse[0];
        // At the end of sample n, this filter stores the tail needed at n+1.
        // The last zero keeps exactly the original history capacity.
        let mut tail: Vec<_> = impulse.into_iter().skip(1).collect();
        tail.push(0.0);
        Ok(Self { direct, tail: DelayedFilter::from_impulse_response(sample_period_s, tail)? })
    }

    /// Prepare a new kernel from prior outgoing waves, oldest first. Replaying
    /// the history also prepares the pending tail; merely loading a ring would
    /// leave `incoming()` at zero until a sample had incorrectly been consumed.
    /// This is control/offline work and allocates; `push` and `response` do not.
    pub fn with_history(sample_period_s: f64, impulse: Vec<f64>, history: &[f64])
        -> Result<Self, ScatteringError>
    {
        if history.iter().any(|x| !x.is_finite()) {
            return Err(ScatteringError::Invalid("outgoing history must be finite"));
        }
        let keep = impulse.len();
        let mut port = Self::new(sample_period_s, impulse)?;
        for &outgoing in &history[history.len().saturating_sub(keep)..] {
            port.push(outgoing)?;
        }
        Ok(port)
    }

    /// Next-sample affine response. Reading never advances the port.
    #[must_use]
    pub fn response(&self) -> ScatteringResponse {
        ScatteringResponse { direct: self.direct, history: self.tail.incoming() }
    }

    /// Commit one outgoing wave and return its SAME-sample reflection. A refusal
    /// leaves the port untouched. No callback allocation or kernel copying.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, ScatteringError> {
        let reflected = self.response().incoming(outgoing)?;
        self.tail.push(outgoing)?;
        Ok(reflected)
    }

    /// Accepted outgoing history, oldest first; allocated only for a control lift.
    #[must_use]
    pub fn history(&self) -> Vec<f64> { self.tail.history() }
}

#[derive(Debug, Clone)]
struct Fade {
    old: FirScatteringPort,
    kernel: Vec<f64>,
    remaining: usize,
    total: usize,
}

/// A bounded menu of FIR responses sharing one wave normalization and history.
/// Only the active port (and one fading predecessor) performs sample work.
/// Preparation and switching allocate; steady stepping does not.
#[derive(Debug, Clone)]
pub struct FirScatteringBank {
    kernels: Vec<Vec<f64>>,
    active: FirScatteringPort,
    index: usize,
    fade: Option<Fade>,
    sample_period_s: f64,
    poisoned: bool,
}
impl FirScatteringBank {
    /// Prepare every response with the SAME history capacity so changing to a
    /// longer kernel cannot discard still-relevant in-flight waves. The caller's
    /// `maximum_total_taps` bounds the complete padded menu before allocation.
    pub fn new(sample_period_s: f64, mut kernels: Vec<Vec<f64>>, maximum_total_taps: usize)
        -> Result<Self, ScatteringError>
    {
        let longest = kernels.iter().map(Vec::len).max()
            .ok_or(ScatteringError::Invalid("reflection menu must not be empty"))?;
        let total = longest.checked_mul(kernels.len())
            .ok_or(ScatteringError::Invalid("reflection menu size overflow"))?;
        if total > maximum_total_taps || kernels.iter().any(|k| k.len() < 4 || k.iter().any(|x| !x.is_finite())) {
            return Err(ScatteringError::Invalid("reflection menu is invalid or exceeds its padded tap budget"));
        }
        for kernel in &mut kernels { kernel.resize(longest, 0.0); }
        let active = FirScatteringPort::new(sample_period_s, kernels[0].clone())?;
        Ok(Self { kernels, active, index: 0, fade: None, sample_period_s, poisoned: false })
    }

    /// Index of the selected target response.
    #[must_use]
    pub const fn active_index(&self) -> usize { self.index }

    /// Number of accepted samples remaining in the current transition.
    #[must_use]
    pub fn fade_remaining(&self) -> usize { self.fade.as_ref().map_or(0, |f| f.remaining) }

    fn check_live(&self) -> Result<(), ScatteringError> {
        if self.poisoned { Err(ScatteringError::Poisoned) } else { Ok(()) }
    }

    /// The affine response used by the very NEXT `push`, including that sample's
    /// crossfade weight. Using the previous weight produces a different boundary.
    pub fn response(&self) -> Result<ScatteringResponse, ScatteringError> {
        self.check_live()?;
        let new = self.active.response();
        let response = match &self.fade {
            Some(fade) => {
                let w = (fade.total - fade.remaining + 1) as f64 / fade.total as f64;
                let old = fade.old.response();
                ScatteringResponse { direct: (1.0-w)*old.direct+w*new.direct,
                    history: (1.0-w)*old.history+w*new.history }
            }
            None => new,
        };
        if !response.direct.is_finite() || !response.history.is_finite() {
            return Err(ScatteringError::Invalid("crossfaded response is not finite"));
        }
        Ok(response)
    }

    /// Accept one outgoing wave. Returns exactly the response supplied to the
    /// junction, rather than advancing one participant and recomputing its phase.
    /// An underlying filter refusal poisons the bank; no partial sample may resume.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, ScatteringError> {
        let incoming = self.response()?.incoming(outgoing)?;
        let advanced = (|| -> Result<(), ScatteringError> {
            self.active.push(outgoing)?;
            if let Some(fade) = &mut self.fade { fade.old.push(outgoing)?; }
            Ok(())
        })();
        if let Err(error) = advanced { self.poisoned = true; return Err(error); }
        if let Some(fade) = &mut self.fade {
            fade.remaining -= 1;
            if fade.remaining == 0 { self.fade = None; }
        }
        Ok(incoming)
    }

    /// Select a kernel without resetting acoustic history. During a second
    /// switch, the predecessor is the CURRENT blended kernel, not either endpoint
    /// of the interrupted transition. Zero fade explicitly changes the boundary
    /// immediately. This is coefficient interpolation, not moving-hole mechanics.
    /// A refused control changes neither port nor its transition/history.
    pub fn switch(&mut self, index: usize, fade_samples: usize) -> Result<(), ScatteringError> {
        self.check_live()?;
        let target = self.kernels.get(index)
            .ok_or(ScatteringError::Invalid("reflection index is outside the prepared menu"))?;
        let history = self.active.history();
        let active = FirScatteringPort::with_history(self.sample_period_s, target.clone(), &history)?;
        let fade = if fade_samples == 0 { None } else {
            let kernel = match &self.fade {
                Some(fade) => {
                    let w = (fade.total-fade.remaining) as f64 / fade.total as f64;
                    fade.kernel.iter().zip(&self.kernels[self.index])
                        .map(|(old, new)| (1.0-w)*old+w*new).collect()
                }
                None => self.kernels[self.index].clone(),
            };
            let old = FirScatteringPort::with_history(self.sample_period_s, kernel.clone(), &history)?;
            Some(Fade { old, kernel, remaining: fade_samples, total: fade_samples })
        };
        self.active = active;
        self.index = index;
        self.fade = fade;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(a: f64, b: f64) {
        assert!((a-b).abs() <= 3e-12*a.abs().max(b.abs()).max(1.0), "{a} != {b}");
    }
    fn signal(i: usize) -> f64 { ((i*17 % 43) as f64 - 21.0)/13.0 }
    fn convolution(kernel: &[f64], history: &[f64]) -> f64 {
        kernel.iter().zip(history.iter().rev()).map(|(a,b)|a*b).sum()
    }

    #[test]
    fn instantaneous_term_and_delayed_echoes_match_the_existing_full_fir() {
        for direct in [-1.0, -0.4, 0.0, 0.7, 1.0] {
            let mut kernel = vec![0.0; 19];
            kernel[0]=direct; kernel[3]=-0.2; kernel[18]=0.04;
            let mut original = DelayedFilter::from_impulse_response(1.0/48000.0, kernel.clone()).unwrap();
            let mut port = FirScatteringPort::new(1.0/48000.0, kernel.clone()).unwrap();
            let mut history = Vec::new();
            for i in 0..200 {
                let x=signal(i); let preview=port.response();
                history.push(x);
                let y=port.push(x).unwrap();
                assert_eq!(y, preview.incoming(x).unwrap());
                close(y, original.push(x).unwrap());
                close(y, convolution(&kernel, &history));
            }
        }
    }

    #[test]
    fn lifted_tail_is_ready_before_the_first_new_sample_and_refusal_is_atomic() {
        let kernel=vec![0.3,-0.2,0.4,0.1,0.05];
        let history:Vec<_>=(0..51).map(signal).collect();
        let mut port=FirScatteringPort::with_history(0.001,kernel.clone(),&history).unwrap();
        let mut expected=history.clone(); expected.push(0.8);
        close(port.response().incoming(0.8).unwrap(),convolution(&kernel,&expected));
        let before=port.response();let remembered=port.history();
        assert!(port.push(f64::NAN).is_err());
        assert_eq!(port.response(),before);assert_eq!(port.history(),remembered);
        close(port.push(0.8).unwrap(),convolution(&kernel,&expected));
    }

    #[test]
    fn switched_kernels_use_the_current_sample_weight_and_complete_history() {
        let kernels=vec![vec![0.2,-0.1,0.0,0.3],vec![-0.4,0.0,0.1,0.2,0.0,0.3,0.0,-0.1]];
        let mut bank=FirScatteringBank::new(0.001,kernels.clone(),16).unwrap();
        let mut history=Vec::new();
        for i in 0..40 {let x=signal(i);history.push(x);bank.push(x).unwrap();}
        bank.switch(1,5).unwrap();
        for i in 0..9 {
            let x=signal(i+40);history.push(x);
            let w=((i+1) as f64/5.0).min(1.0);
            let expected=(1.0-w)*convolution(&kernels[0],&history)+w*convolution(&kernels[1],&history);
            let response=bank.response().unwrap();
            close(response.incoming(x).unwrap(),expected);
            assert_eq!(bank.push(x).unwrap(),response.incoming(x).unwrap());
        }
        assert_eq!(bank.fade_remaining(),0);
        let before=bank.response().unwrap();
        assert!(bank.switch(usize::MAX,4).is_err());assert_eq!(bank.response().unwrap(),before);
    }

    #[test]
    fn interrupted_fade_starts_from_the_actual_mixture_not_the_old_target() {
        let kernels=vec![vec![0.3,0.1,0.0,0.2],vec![-0.2,0.0,0.4,0.0],vec![0.1,-0.2,0.3,0.1]];
        let mut bank=FirScatteringBank::new(0.001,kernels.clone(),12).unwrap();
        let mut history=Vec::new();
        for i in 0..13 {let x=signal(i);history.push(x);bank.push(x).unwrap();}
        bank.switch(1,4).unwrap();
        for i in 13..15 {let x=signal(i);history.push(x);bank.push(x).unwrap();}
        bank.switch(2,3).unwrap();
        for i in 0..3 {
            let x=signal(15+i);history.push(x);let w=(i+1) as f64/3.0;
            let previous=0.5*convolution(&kernels[0],&history)+0.5*convolution(&kernels[1],&history);
            close(bank.push(x).unwrap(),(1.0-w)*previous+w*convolution(&kernels[2],&history));
        }
    }

    #[test]
    fn constant_reflectance_closes_the_boundary_without_a_spurious_delay() {
        // Matched resistive source p_drive = p_in + R_source*q, with Zc=R_source=1:
        // a=p_drive/2, b=r*a. r=+/-1 must be rigid/open THIS sample, not next.
        for r in [-1.0,0.0,1.0] {
            let mut port=FirScatteringPort::new(0.001,vec![r,0.0,0.0,0.0]).unwrap();
            for i in 0..16 {
                let source=signal(i);let a=0.5*source;let b=port.push(a).unwrap();
                close(a+b,0.5*(1.0+r)*source);
                close(a-b,0.5*(1.0-r)*source);
            }
        }
    }

    #[test]
    fn padded_menu_budget_and_nonfinite_inputs_refuse() {
        assert!(FirScatteringBank::new(0.001,vec![],10).is_err());
        assert!(FirScatteringBank::new(0.001,vec![vec![0.0;4],vec![0.0;8]],15).is_err());
        assert!(FirScatteringPort::new(0.0,vec![0.0;4]).is_err());
        assert!(FirScatteringPort::new(0.001,vec![f64::NAN;4]).is_err());
        let mut bank=FirScatteringBank::new(0.001,vec![vec![0.2,0.0,0.0,0.0]],4).unwrap();
        let before=bank.response().unwrap();
        assert!(bank.push(f64::INFINITY).is_err());assert_eq!(bank.response().unwrap(),before);
    }
}
