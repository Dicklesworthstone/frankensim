//! A characteristic line whose taps move between frozen shapes during one render (pexfp.327.1).
//!
//! Each shape of a moving tract is realized at build time as its own exact FIR. A switch carries the outgoing-wave
//! history into the next shape's line with [`DelayedFilter::from_impulse_response_with_history`], the O(N) lift of
//! the brass valve bank (`MmLineBank::switch`): the ring buffer is preloaded, so waves already in flight persist and
//! the new taps govern every later reflection. The incoming waves then crossfade linearly from the old line over the
//! fade. Replaying the history through the new FIR instead (as `FirScatteringBank::switch` does) costs N^2 per
//! switch, which at the voice's 131072-tap glottal line is about 1.7e10 multiply-adds.
//!
//! Every shape is padded with zeros to the longest, so a switch to a longer shape cannot drop a wave the longer
//! taps still need. The padding is bit-neutral: [`DelayedFilter::push`] accumulates in tap order, the zeros come
//! last, and a sum that starts at `+0.0` is never `-0.0`, so adding `0.0 * x` leaves it bit for bit unchanged.
//!
//! WHAT THIS IS NOT: a crossfade between frozen linear shapes is coefficient interpolation, not a moving boundary.
//! A tract whose walls move also injects the volume they sweep; that source is absent here.

use fs_vfit::discretize::{DelayedFilter, DiscretizeError};

/// Refusal of a keyframe line, by name.
#[derive(Debug, Clone, PartialEq)]
pub enum KeyframeLineError {
    /// Bad shapes or an invalid switch.
    Invalid(&'static str),
    /// The underlying FIR refused a tap, a history sample or a wave.
    Discrete(DiscretizeError),
    /// An earlier push failed after one of the two lines of a fade had advanced; rebuild the line.
    Poisoned,
}

impl core::fmt::Display for KeyframeLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "keyframe line: {what}"),
            Self::Discrete(error) => write!(f, "keyframe line: {error:?}"),
            Self::Poisoned => write!(f, "keyframe line has a partially accepted sample"),
        }
    }
}

impl std::error::Error for KeyframeLineError {}

impl From<DiscretizeError> for KeyframeLineError {
    fn from(error: DiscretizeError) -> Self {
        Self::Discrete(error)
    }
}

#[derive(Debug, Clone)]
struct Fade {
    old: DelayedFilter,
    remaining: usize,
    total: usize,
}

/// A menu of frozen FIR shapes sharing one outgoing-wave history, one of them active.
#[derive(Debug, Clone)]
pub struct KeyframeLine {
    taps: Vec<Vec<f64>>,
    active: DelayedFilter,
    index: usize,
    fade: Option<Fade>,
    sample_period_s: f64,
    poisoned: bool,
}

impl KeyframeLine {
    /// Shape 0 starts active with an empty history.
    ///
    /// # Errors
    /// No shapes, a shape under four taps, a non-finite tap, or a bad sample period.
    pub fn new(sample_period_s: f64, mut taps: Vec<Vec<f64>>) -> Result<Self, KeyframeLineError> {
        let longest = taps
            .iter()
            .map(Vec::len)
            .max()
            .ok_or(KeyframeLineError::Invalid("a keyframe line needs at least one shape"))?;
        if taps.iter().any(|t| t.len() < 4) {
            return Err(KeyframeLineError::Invalid("every shape needs at least four taps"));
        }
        if taps.iter().flatten().any(|x| !x.is_finite()) {
            return Err(KeyframeLineError::Invalid("every tap must be finite"));
        }
        for t in &mut taps {
            t.resize(longest, 0.0);
        }
        let active = DelayedFilter::from_impulse_response(sample_period_s, taps[0].clone())?;
        Ok(Self {
            taps,
            active,
            index: 0,
            fade: None,
            sample_period_s,
            poisoned: false,
        })
    }

    /// The shape the line is moving to, or holds.
    #[must_use]
    pub const fn active_index(&self) -> usize {
        self.index
    }

    /// Number of shapes.
    #[must_use]
    pub fn shapes(&self) -> usize {
        self.taps.len()
    }

    /// Taps of every shape after padding.
    #[must_use]
    pub fn padded_taps(&self) -> usize {
        self.taps[0].len()
    }

    /// Samples left in the present crossfade.
    #[must_use]
    pub fn fade_remaining(&self) -> usize {
        self.fade.as_ref().map_or(0, |f| f.remaining)
    }

    /// Move to shape `index`. The outgoing-wave history is carried into the new shape's line; over the next
    /// `fade_samples` pushes the incoming wave goes linearly from the old line's to the new one's, and `0` swaps at
    /// once. A refused switch changes nothing.
    ///
    /// # Errors
    /// An unknown shape, a switch while a crossfade is still running, or a poisoned line.
    pub fn switch(&mut self, index: usize, fade_samples: usize) -> Result<(), KeyframeLineError> {
        if self.poisoned {
            return Err(KeyframeLineError::Poisoned);
        }
        let taps = self
            .taps
            .get(index)
            .ok_or(KeyframeLineError::Invalid("shape index is outside the prepared menu"))?;
        if self.fade.is_some() {
            return Err(KeyframeLineError::Invalid("a switch must wait for the running crossfade to finish"));
        }
        let lifted =
            DelayedFilter::from_impulse_response_with_history(self.sample_period_s, taps.clone(), &self.active.history())?;
        let old = core::mem::replace(&mut self.active, lifted);
        self.fade = (fade_samples > 0).then_some(Fade {
            old,
            remaining: fade_samples,
            total: fade_samples,
        });
        self.index = index;
        Ok(())
    }

    /// Push one outgoing wave and return the incoming one, the same-sample convolution of
    /// [`DelayedFilter::push`]. During a crossfade it is `(1 - w) old + w new` with `w = (total - remaining + 1) /
    /// total`, the weight of this sample, so the last sample of a fade is the new line's alone.
    ///
    /// # Errors
    /// A refused wave; one refused after the other line of a fade advanced poisons the line.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, KeyframeLineError> {
        if self.poisoned {
            return Err(KeyframeLineError::Poisoned);
        }
        let Some(fade) = self.fade.as_mut() else {
            return Ok(self.active.push(outgoing)?);
        };
        let old = fade.old.push(outgoing)?;
        let new = match self.active.push(outgoing) {
            Ok(new) => new,
            Err(error) => {
                self.poisoned = true;
                return Err(error.into());
            }
        };
        let w = (fade.total - fade.remaining + 1) as f64 / fade.total as f64;
        fade.remaining -= 1;
        if fade.remaining == 0 {
            self.fade = None;
        }
        Ok((1.0 - w) * old + w * new)
    }
}

/// A tract line as an owner holds it: one shape on the plain FIR path every line took before moving tracts, or a
/// menu of shapes it moves through.
#[derive(Debug, Clone)]
pub enum ShapeLine {
    /// One shape, pushed exactly as before.
    Fixed(DelayedFilter),
    /// Shapes switched with the history carried.
    Moving(KeyframeLine),
}

impl ShapeLine {
    /// Push one outgoing wave and return the incoming one.
    ///
    /// # Errors
    /// The line's own refusal.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, KeyframeLineError> {
        match self {
            Self::Fixed(line) => Ok(line.push(outgoing)?),
            Self::Moving(line) => line.push(outgoing),
        }
    }

    /// The moving line made of a fixed line's own taps (shape 0) and `further` shapes.
    ///
    /// # Errors
    /// A line that already moves, one that is not a plain un-stepped FIR, or refused shapes.
    pub fn to_moving(&self, further: Vec<Vec<f64>>) -> Result<Self, KeyframeLineError> {
        let Self::Fixed(line) = self else {
            return Err(KeyframeLineError::Invalid("the line already moves between shapes"));
        };
        let (dt, own) = line
            .quiescent_fir()
            .ok_or(KeyframeLineError::Invalid("a moving tract needs every shape's line to be a plain FIR"))?;
        let mut menu = Vec::with_capacity(further.len() + 1);
        menu.push(own.to_vec());
        menu.extend(further);
        Ok(Self::Moving(KeyframeLine::new(dt, menu)?))
    }

    /// Switch a moving line ([`KeyframeLine::switch`]).
    ///
    /// # Errors
    /// A fixed line, or the switch's own refusal.
    pub fn switch(&mut self, index: usize, fade_samples: usize) -> Result<(), KeyframeLineError> {
        match self {
            Self::Fixed(_) => Err(KeyframeLineError::Invalid("a fixed line has no shapes to switch between")),
            Self::Moving(line) => line.switch(index, fade_samples),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DT: f64 = 1.0 / 48_000.0;

    fn signal(i: usize) -> f64 {
        ((i * 17 % 43) as f64 - 21.0) / 13.0
    }

    fn convolution(taps: &[f64], history: &[f64]) -> f64 {
        taps.iter().zip(history.iter().rev()).map(|(a, b)| a * b).sum()
    }

    fn shapes() -> Vec<Vec<f64>> {
        vec![
            vec![0.2, -0.1, 0.0, 0.3, 0.05],
            vec![-0.4, 0.0, 0.1, 0.2, 0.0, 0.3, 0.0, -0.1, 0.02],
            vec![0.1, 0.1, -0.2, 0.0],
        ]
    }

    #[test]
    fn one_shape_is_bit_identical_to_its_plain_line() {
        let taps = shapes()[1].clone();
        let mut plain = DelayedFilter::from_impulse_response(DT, taps.clone()).expect("plain");
        let mut line = KeyframeLine::new(DT, vec![taps]).expect("one shape");
        for i in 0..2000 {
            let x = signal(i);
            assert_eq!(line.push(x).expect("line").to_bits(), plain.push(x).expect("plain").to_bits());
        }
    }

    #[test]
    fn zero_padding_changes_no_output_bit() {
        // Shape 0 is the 5-tap kernel padded to 9 by the longer shape 1; it must play as its own plain line.
        let mut plain = DelayedFilter::from_impulse_response(DT, shapes()[0].clone()).expect("plain");
        let mut line = KeyframeLine::new(DT, shapes()).expect("menu");
        assert_eq!(line.padded_taps(), 9);
        for i in 0..2000 {
            let x = signal(i) * if i % 7 == 0 { -1.0 } else { 1.0 };
            assert_eq!(line.push(x).expect("line").to_bits(), plain.push(x).expect("plain").to_bits());
        }
    }

    #[test]
    fn a_switch_carries_the_whole_history_and_the_crossfade_weights_are_exact() {
        let menu = shapes();
        let longest = menu.iter().map(Vec::len).max().expect("menu");
        let padded = |k: usize| {
            let mut t = menu[k].clone();
            t.resize(longest, 0.0);
            t
        };
        let mut line = KeyframeLine::new(DT, menu.clone()).expect("menu");
        let mut history = Vec::new();
        for i in 0..40 {
            let x = signal(i);
            history.push(x);
            line.push(x).expect("push");
        }
        // Hard swap: the very next output is the new taps over the whole history.
        line.switch(1, 0).expect("hard");
        assert_eq!(line.fade_remaining(), 0);
        let x = signal(40);
        history.push(x);
        let y = line.push(x).expect("after hard");
        assert!((y - convolution(&padded(1), &history)).abs() < 1e-12, "{y}");
        // Crossfade over 5 samples to shape 2: w = 1/5 .. 5/5, the last sample is the new line alone.
        line.switch(2, 5).expect("fade");
        assert!(matches!(
            line.switch(0, 3),
            Err(KeyframeLineError::Invalid("a switch must wait for the running crossfade to finish"))
        ));
        for k in 0..7 {
            let x = signal(41 + k);
            history.push(x);
            let y = line.push(x).expect("fading");
            let w = ((k + 1) as f64 / 5.0).min(1.0);
            let want = (1.0 - w) * convolution(&padded(1), &history) + w * convolution(&padded(2), &history);
            assert!((y - want).abs() < 1e-12, "sample {k}: {y} vs {want}");
        }
        assert_eq!(line.fade_remaining(), 0);
        assert_eq!(line.active_index(), 2);
    }

    #[test]
    fn refusals_fire_by_name_and_change_nothing() {
        assert!(matches!(
            KeyframeLine::new(DT, vec![]),
            Err(KeyframeLineError::Invalid("a keyframe line needs at least one shape"))
        ));
        assert!(matches!(
            KeyframeLine::new(DT, vec![vec![0.1, 0.2, 0.3]]),
            Err(KeyframeLineError::Invalid("every shape needs at least four taps"))
        ));
        assert!(matches!(
            KeyframeLine::new(DT, vec![vec![0.1, f64::NAN, 0.3, 0.0]]),
            Err(KeyframeLineError::Invalid("every tap must be finite"))
        ));
        let mut line = KeyframeLine::new(DT, shapes()).expect("menu");
        assert!(matches!(
            line.switch(3, 0),
            Err(KeyframeLineError::Invalid("shape index is outside the prepared menu"))
        ));
        assert_eq!(line.active_index(), 0);
        assert!(line.push(f64::NAN).is_err());
        assert!(line.push(0.5).is_ok());
    }
}
