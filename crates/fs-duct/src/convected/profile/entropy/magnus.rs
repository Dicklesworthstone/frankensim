//! Fixed-size upper-block Magnus action for one advected entropy channel.
//! All three variables have the same reference scale before exponentiation.

use crate::convected::{bad, finite, phase};
use crate::{DuctError, mul2};
use fs_math::c64::C64;

#[derive(Clone, Copy)]
pub(super) struct Block {
    pub acoustic: [C64; 4],
    pub forcing: [C64; 2],
    pub entropy: C64,
}

impl Block {
    const IDENTITY: Self = Self {
        acoustic: [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE],
        forcing: [C64::ZERO; 2],
        entropy: C64::ONE,
    };

    pub fn is_finite(self) -> bool {
        self.acoustic
            .into_iter()
            .chain(self.forcing)
            .chain([self.entropy])
            .all(finite)
    }

    fn add(self, other: Self) -> Self {
        Self {
            acoustic: core::array::from_fn(|i| self.acoustic[i] + other.acoustic[i]),
            forcing: core::array::from_fn(|i| self.forcing[i] + other.forcing[i]),
            entropy: self.entropy + other.entropy,
        }
    }

    fn scale(self, value: f64) -> Self {
        Self {
            acoustic: self.acoustic.map(|z| z.scale(value)),
            forcing: self.forcing.map(|z| z.scale(value)),
            entropy: self.entropy.scale(value),
        }
    }

    fn mul(self, other: Self) -> Self {
        Self {
            acoustic: mul2(self.acoustic, other.acoustic),
            forcing: [
                self.acoustic[0] * other.forcing[0]
                    + self.acoustic[1] * other.forcing[1]
                    + self.forcing[0] * other.entropy,
                self.acoustic[2] * other.forcing[0]
                    + self.acoustic[3] * other.forcing[1]
                    + self.forcing[1] * other.entropy,
            ],
            entropy: self.entropy * other.entropy,
        }
    }

    fn norm(self) -> f64 {
        (self.acoustic[0].abs() + self.acoustic[1].abs() + self.forcing[0].abs())
            .max(self.acoustic[2].abs() + self.acoustic[3].abs() + self.forcing[1].abs())
            .max(self.entropy.abs())
    }
}

pub(super) fn step(a: Block, b: Block, h: f64, root3: f64) -> Result<Block, DuctError> {
    let a = a.scale(h);
    let b = b.scale(h);
    let logarithm = a
        .add(b)
        .scale(0.5)
        .add(a.mul(b).add(b.mul(a).scale(-1.0)).scale(-root3 / 12.0));
    if !logarithm.is_finite() {
        return Err(bad(
            "entropy Magnus generator overflow; refine the original panels",
        ));
    }
    // Bounded scaling and squaring. For norm <= 1/2, 24 Taylor terms leave a
    // truncation tail < 4e-33 before floating-point rounding. This controls the
    // LOCAL exponential only, not the spatial Magnus/travel-phase error.
    let mut scaled = logarithm;
    let mut squares = 0;
    while scaled.norm() > 0.5 {
        if squares == 32 {
            return Err(bad(
                "entropy exponential work bound exceeded; refine the original panels",
            ));
        }
        scaled = scaled.scale(0.5);
        squares += 1;
    }
    let mut sum = Block::IDENTITY;
    let mut term = Block::IDENTITY;
    for order in 1..=24 {
        term = term.mul(scaled).scale(1.0 / f64::from(order));
        sum = sum.add(term);
    }
    for _ in 0..squares {
        sum = sum.mul(sum);
    }
    // The scalar entropy block has an exactly imaginary logarithm; evaluating
    // its phase directly avoids artificial transport gain from squaring.
    sum.entropy = phase(logarithm.entropy.im);
    if !sum.is_finite() {
        return Err(bad("entropy Magnus exponential is unrepresentable"));
    }
    Ok(sum)
}
