//! Equivalent lengths of a compact tone hole: the engine's laws and the printed thick-wall set.
//!
//! A hole of radius `b`, chimney height `t` and bore radius `a` (`d = b/a`) is a symmetric T:
//! a shunt branch (the chimney with its inner matching, and for an open hole its radiating
//! mouth) between two halves of a series impedance `j k Z0 t_a` that acts on the bore.
//!
//! `Engine` is the model as it was: inner length from the Dubos et al. 1999 polynomial
//! (`fs_phs::side_hole_inner_length`), the series length `-0.37 b^2/a` on an open hole only, a
//! closed hole a pure shunt. `Printed` takes the Lefebvre and Scavone, J. Acoust. Soc. Am.
//! 131(4), 3153-3163 (2012) set instead, each term with its printed equation:
//!
//! - `t_i` open (Eq. 31 times `G` of Eq. 32): `t_i/b = 0.822 - 0.095 d - 1.566 d^2 + 2.138 d^3 -
//!   1.640 d^4 + 0.502 d^5`, `G = 1 + H(d) I(ka)`, `H = 1 - 4.56 d + 6.55 d^2`,
//!   `I = 0.17 ka + 0.92 (ka)^2 + 0.16 (ka)^3 - 0.29 (ka)^4`. Printed for `ka <= 1` and `d <= 1`:
//!   above that `I` is held at its `ka = 1` value.
//! - `t_m` on open and closed holes (Eq. 8, Nederveen et al. 1998): `t_m = b d (1 + 0.207 d^3)/8`.
//! - `t_r` (Eq. 11, Dalmont and Nederveen 2001): `t_r = (0.822 - 0.47 [b/(a+t)]^0.8) b`. The flanged
//!   load keeps its radiation resistance and its frequency dependence; only its low-frequency mass
//!   `0.8216 b` moves, as a constant length `t_r - 0.8216 b`.
//! - `t_a` open (Eq. 33) and closed (Eq. 34): `t_a/(b d^2) = -0.35 + 0.06 tanh(2.7 t/b)` and
//!   `-0.12 - 0.17 tanh(2.4 t/b)`. A closed hole gains a series inertance the engine never had.
//!
//! The open shunt is then Eq. 7, `j Z0h [k t_i + tan k(t + t_m + t_r)]`, and the closed one Eq. 16,
//! `-j Z0h cot k(t + t_m)`: the closed-hole `k t_i` mass of Eq. 16 is not applied, the paper calls
//! it negligible at low frequency and the scratch model puts it under 0.01 cent at D3, F4 and G4.
//! The `-t_a/4` of Eq. 13 is not added to `t_i`: Eq. 31 is the fit used. A hole with `b >= a`
//! keeps the engine terms (the fits are printed for `d <= 1`). A partly open hole weights the
//! open and closed series lengths by `sigma`, as the engine weights its open series length. The
//! Nederveen near-field mutual length between open holes is kept as it is.

use fs_math::det;

/// The flanged-pipe end correction `0.8216 b` that the engine's radiation load carries at low frequency.
const FLANGED_END_CORRECTION: f64 = 0.8216;

/// Which equivalent lengths a compact tone hole carries.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum HoleLengths {
    /// No bore-junction corrections: the chimney and the radiation load alone.
    Uncorrected,
    /// The engine's laws, unchanged.
    #[default]
    Engine,
    /// The Lefebvre and Scavone 2012 thick-wall set described in the module documentation.
    Printed,
}

impl HoleLengths {
    pub(crate) const fn from_compact(compact: bool) -> Self {
        if compact { Self::Engine } else { Self::Uncorrected }
    }

    /// `Printed` needs the hole inside the bore (`b < a`); otherwise the engine's terms apply.
    fn effective(self, hole_radius: f64, bore_radius: f64) -> Self {
        if self == Self::Printed && !(hole_radius > 0.0 && bore_radius > hole_radius && bore_radius.is_finite()) {
            Self::Engine
        } else {
            self
        }
    }

    /// Length added to the physical height of the OPEN part of the chimney, before its flanged mouth.
    /// `wavenumber` is `omega / c` of the hole's gas.
    pub(crate) fn open_extra_m(self, hole_radius: f64, chimney_height: f64, bore_radius: f64, wavenumber: f64) -> f64 {
        match self.effective(hole_radius, bore_radius) {
            Self::Uncorrected => 0.0,
            Self::Engine => fs_phs::side_hole_inner_length(hole_radius, bore_radius),
            Self::Printed => {
                let (b, a) = (hole_radius, bore_radius);
                let d = b / a;
                let inner = inner_over_b(d) * frequency_factor(d, wavenumber * a) * b;
                inner + matching_m(b, d) + (radiation_m(b, a, chimney_height) - FLANGED_END_CORRECTION * b)
            }
        }
    }

    /// Length added to the height of the SEALED part of the chimney.
    pub(crate) fn sealed_extra_m(self, hole_radius: f64, bore_radius: f64) -> f64 {
        match self.effective(hole_radius, bore_radius) {
            Self::Printed => matching_m(hole_radius, hole_radius / bore_radius),
            Self::Uncorrected | Self::Engine => 0.0,
        }
    }

    /// Signed series length referred to the bore (its inertance is `rho * length / (pi a^2)`), the
    /// open and closed parts weighted by `sigma`, plus the mutual length of the open holes.
    pub(crate) fn series_m(
        self, hole_radius: f64, chimney_height: f64, bore_radius: f64, sigma: f64, mutual_m: f64,
    ) -> f64 {
        match self.effective(hole_radius, bore_radius) {
            Self::Uncorrected => 0.0,
            Self::Engine => {
                if sigma > 0.0 { fs_phs::side_hole_series_length(hole_radius, bore_radius) * sigma + mutual_m } else { 0.0 }
            }
            Self::Printed => {
                let (b, a, t) = (hole_radius, bore_radius, chimney_height);
                let d = b / a;
                let weighted = sigma * ((-0.35 + 0.06 * det::tanh(2.7 * t / b)) * b * d * d)
                    + (1.0 - sigma) * ((-0.12 - 0.17 * det::tanh(2.4 * t / b)) * b * d * d);
                if sigma > 0.0 { weighted + mutual_m } else { weighted }
            }
        }
    }
}

/// Eq. 31: the inner length over `b`.
fn inner_over_b(d: f64) -> f64 {
    0.822 + d * (-0.095 + d * (-1.566 + d * (2.138 + d * (-1.640 + d * 0.502))))
}

/// Eq. 32: `G = 1 + H(d) I(ka)`, `ka` held inside the printed range `0..=1`.
fn frequency_factor(d: f64, ka: f64) -> f64 {
    let x = ka.clamp(0.0, 1.0);
    let h = 1.0 + d * (-4.56 + d * 6.55);
    let i = x * (0.17 + x * (0.92 + x * (0.16 - 0.29 * x)));
    1.0 + h * i
}

/// Eq. 8: the matching-volume length.
fn matching_m(b: f64, d: f64) -> f64 {
    b * d * (1.0 + 0.207 * d * d * d) / 8.0
}

/// Eq. 11: the radiation length of a hole drilled through a thick pipe of outer radius `a + t`.
fn radiation_m(b: f64, a: f64, t: f64) -> f64 {
    (0.822 - 0.47 * det::pow(b / (a + t), 0.8)) * b
}

#[cfg(test)]
mod tests;
