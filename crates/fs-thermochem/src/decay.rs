//! Constant-hazard, first-order decay shared by physical and chemical owners.
//!
//! `lambda = ln(2)/t_half`, survival is `exp(-lambda*t)`, and activity is
//! `lambda*N` for a population of N nuclei. See IAEA TECDOC 1741,
//! <https://www-pub.iaea.org/MTCD/Publications/PDF/TE-1741_web.pdf>.
//! Half-life describes a population's exponential survival, not a deterministic
//! lifetime assigned to each constituent. A stable species is explicit and
//! never encoded as an infinite input. Rates and their applicable conditions
//! are caller-authored; this module supplies no isotope catalog or calibration.
//!
//! This immutable L1 closure evaluates one law. Inventory evolution, daughter
//! branching, thermal coupling and changing-rate integration belong to consumers.

use core::f64::consts::LN_2;
use fs_math::decay::{decayed_fraction, remaining_fraction};
use fs_qty::{Frequency, Time};

/// Invalid physical input or a required nonzero value outside binary64 range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DecayErrorV1 {
    /// Half-life must be finite and positive; rate/time finite and nonnegative.
    Invalid(&'static str),
    /// A positive rate, lifetime or decay fraction cannot be represented.
    Arithmetic(&'static str),
}

impl core::fmt::Display for DecayErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "first-order decay: {self:?}")
    }
}
impl std::error::Error for DecayErrorV1 {}

/// Separately evaluated survival and loss fractions over one interval.
///
/// `decayed` uses `-expm1(-lambda*t)`, preserving losses much smaller than one
/// ULP of `remaining`. Consumers must use it directly for produced daughters or
/// released energy, rather than subtracting `remaining` from one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayFractionsV1 {
    /// Undecayed fraction, between zero and one.
    pub remaining: f64,
    /// Decayed fraction, between zero and one.
    pub decayed: f64,
}

/// An admitted, immutable constant first-order removal law.
///
/// Both parameterizations are retained in coherent SI. A supplied half-life is
/// preserved exactly; its derived inverse-second rate is rounded once. The
/// constructor refuses an unrepresentable positive rate or finite half-life.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayLawV1 {
    rate: Frequency,
    half_life: Option<Time>,
}

impl DecayLawV1 {
    /// Explicitly stable population: zero rate and no finite half-life.
    #[must_use]
    pub const fn stable() -> Self {
        Self {
            rate: Frequency::new(0.0),
            half_life: None,
        }
    }

    /// Admit a positive finite half-life in seconds.
    pub fn from_half_life(half_life: Time) -> Result<Self, DecayErrorV1> {
        if !half_life.value().is_finite() || half_life.value() <= 0.0 {
            return Err(DecayErrorV1::Invalid(
                "half-life must be finite and positive",
            ));
        }
        let rate = LN_2 / half_life.value();
        if !rate.is_finite() || rate == 0.0 {
            return Err(DecayErrorV1::Arithmetic("decay constant from half-life"));
        }
        Ok(Self {
            rate: Frequency::new(rate),
            half_life: Some(half_life),
        })
    }

    /// Admit a finite nonnegative rate in inverse seconds; zero is stable.
    pub fn from_decay_constant(rate: Frequency) -> Result<Self, DecayErrorV1> {
        if !rate.value().is_finite() || rate.value() < 0.0 {
            return Err(DecayErrorV1::Invalid(
                "decay rate must be finite and nonnegative",
            ));
        }
        if rate.value() == 0.0 {
            return Ok(Self::stable());
        }
        let half_life = LN_2 / rate.value();
        if !half_life.is_finite() || half_life == 0.0 {
            return Err(DecayErrorV1::Arithmetic("half-life from decay constant"));
        }
        Ok(Self {
            rate,
            half_life: Some(Time::new(half_life)),
        })
    }

    /// Admitted inverse-second rate; this is an event hazard, not an angular rate.
    #[must_use]
    pub const fn decay_constant(self) -> Frequency {
        self.rate
    }

    /// Finite half-life in seconds; `None` means the explicitly stable law.
    #[must_use]
    pub const fn half_life(self) -> Option<Time> {
        self.half_life
    }

    /// Mean survival time `1/lambda`, distinct from the half-life.
    /// Stable populations return `None`; a finite law whose mean lifetime
    /// exceeds binary64 range returns an arithmetic refusal.
    pub fn mean_lifetime(self) -> Result<Option<Time>, DecayErrorV1> {
        if self.rate.value() == 0.0 {
            return Ok(None);
        }
        let lifetime = 1.0 / self.rate.value();
        if !lifetime.is_finite() || lifetime == 0.0 {
            return Err(DecayErrorV1::Arithmetic("mean lifetime"));
        }
        Ok(Some(Time::new(lifetime)))
    }

    /// Combine independent competing removal processes by adding their hazards.
    /// This gives `1/t_effective = 1/t_1 + 1/t_2`; it does not describe
    /// consecutive daughter decays or erase the attribution of either channel.
    pub fn competing_with(self, other: Self) -> Result<Self, DecayErrorV1> {
        if self.rate.value() == 0.0 {
            return Ok(other);
        }
        if other.rate.value() == 0.0 {
            return Ok(self);
        }
        let rate = self.rate.value() + other.rate.value();
        if !rate.is_finite() {
            return Err(DecayErrorV1::Arithmetic("sum of competing decay rates"));
        }
        Self::from_decay_constant(Frequency::new(rate))
    }

    /// Evaluate survival and loss with no allocation or timestep approximation.
    ///
    /// A positive exposure below binary64 range refuses rather than reporting
    /// exactly zero loss. Very large exposures saturate to full decay;
    /// an unrepresentably small survivor rounds to zero. Separate loss and
    /// survival preserve both tails whenever each is representable.
    pub fn fractions(self, elapsed: Time) -> Result<DecayFractionsV1, DecayErrorV1> {
        if !elapsed.value().is_finite() || elapsed.value() < 0.0 {
            return Err(DecayErrorV1::Invalid(
                "elapsed time must be finite and nonnegative",
            ));
        }
        if elapsed.value() == 0.0 || self.rate.value() == 0.0 {
            return Ok(DecayFractionsV1 {
                remaining: 1.0,
                decayed: 0.0,
            });
        }
        let exposure = self.rate.value() * elapsed.value();
        if exposure == 0.0 {
            return Err(DecayErrorV1::Arithmetic(
                "positive integrated decay exposure",
            ));
        }
        Ok(DecayFractionsV1 {
            remaining: remaining_fraction(exposure),
            decayed: decayed_fraction(exposure),
        })
    }
}

#[cfg(test)]
mod tests;
