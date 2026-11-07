//! Bounded homogeneous E/H evolution; exp(+i omega t), E'=-i A H, H'=-i B E.
use fs_math::{c64::C64, dd::Dd, det};

#[derive(Clone, Copy)]
enum Interior {
    Entire {
        a: C64,
        b: C64,
    },
    Decaying {
        admittance: C64,
        forward: C64,
        backward: C64,
    },
}

#[derive(Clone, Copy)]
pub(in crate::plane_wave::layers) struct LayerFieldKernel {
    thickness: f64,
    kz: C64,
    // Retained real q*k for phase; kz still owns admittance and attenuation.
    phase_wavenumber: Dd,
    left: [C64; 2],
    right: [C64; 2],
    interior: Interior,
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

impl LayerFieldKernel {
    pub(in crate::plane_wave::layers) fn new(
        thickness: f64,
        kz: C64,
        phase_wavenumber: Dd,
        a: C64,
        b: C64,
        left: [C64; 2],
        right: [C64; 2],
    ) -> Option<Self> {
        if !thickness.is_finite()
            || thickness < 0.0
            || kz.re < 0.0
            || kz.im > 0.0
            || phase_wavenumber.hi != kz.re
            || !phase_wavenumber.lo.is_finite()
            || ![kz, a, b]
                .iter()
                .chain(&left)
                .chain(&right)
                .all(|&v| finite(v))
        {
            return None;
        }
        let interior = if kz.abs() * thickness <= 0.5 {
            Interior::Entire { a, b }
        } else {
            // TE: Y=kz/A, TM: Y=B/kz. Pick the larger coefficient to
            // avoid dividing by a tiny A or B near its critical limit.
            let admittance = if a.abs() >= b.abs() { kz / a } else { b / kz };
            if admittance == C64::ZERO || !finite(admittance) {
                return None;
            }
            let forward = (left[0] + left[1] / admittance).scale(0.5);
            let backward = (right[0] - right[1] / admittance).scale(0.5);
            if !finite(forward) || !finite(backward) {
                return None;
            }
            Interior::Decaying {
                admittance,
                forward,
                backward,
            }
        };
        Some(Self {
            thickness,
            kz,
            phase_wavenumber,
            left,
            right,
            interior,
        })
    }

    /// Bind this original linear solution to one actual coherent amplitude.
    /// Material coefficients, physical phases and both face locations stay fixed.
    pub(super) fn scaled(self, amplitude: C64) -> Option<Self> {
        if !finite(amplitude) {
            return None;
        }
        let left = self.left.map(|value| value * amplitude);
        let right = self.right.map(|value| value * amplitude);
        let interior = match self.interior {
            Interior::Entire { a, b } => Interior::Entire { a, b },
            Interior::Decaying {
                admittance,
                forward,
                backward,
            } => {
                let forward = forward * amplitude;
                let backward = backward * amplitude;
                if !finite(forward) || !finite(backward) {
                    return None;
                }
                Interior::Decaying {
                    admittance,
                    forward,
                    backward,
                }
            }
        };
        left.iter()
            .chain(&right)
            .all(|&value| finite(value))
            .then_some(Self {
                left,
                right,
                interior,
                ..self
            })
    }

    /// Conservative whole-film magnitudes of the original tangential E/H.
    /// This analytic envelope is for material-domain admission, not an
    /// interval certificate or a claim that its maximum is attained.
    pub(super) fn amplitude_bounds(self) -> Option<[f64; 2]> {
        let endpoints =
            core::array::from_fn::<_, 2, _>(|i| self.left[i].abs().max(self.right[i].abs()));
        let values = if self.thickness == 0.0 || self.kz == C64::ZERO {
            // At criticality the entire field is affine, so its norm is
            // convex and is bounded by its two original endpoint norms.
            endpoints
        } else {
            match self.interior {
                Interior::Entire { a, b } => {
                    // |cos(kz*z)| and |sinc(kz*z)| <= exp(|Im(kz)|*d).
                    // Propagate from either original endpoint; both bounds
                    // are valid, so the smaller one retains tighter admission.
                    let growth = det::exp(-self.kz.im * self.thickness);
                    let coupling = [a.abs(), b.abs()];
                    core::array::from_fn(|i| {
                        let bound = |trace: [C64; 2]| {
                            growth
                                * (trace[i].abs()
                                    + (coupling[i] * self.thickness) * trace[1 - i].abs())
                        };
                        bound(self.left).min(bound(self.right)).max(endpoints[i])
                    })
                }
                Interior::Decaying {
                    admittance,
                    forward,
                    backward,
                } => {
                    // Each wave decays from its own face throughout the film.
                    let electric = forward.abs() + backward.abs();
                    [
                        electric.max(endpoints[0]),
                        (admittance.abs() * electric).max(endpoints[1]),
                    ]
                }
            }
        };
        values
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0)
            .then_some(values)
    }

    pub(in crate::plane_wave::layers) fn field(self, depth: f64) -> Option<[C64; 2]> {
        if !depth.is_finite() || depth < 0.0 || depth > self.thickness {
            return None;
        }
        if depth == 0.0 {
            return Some(self.left);
        }
        if depth == self.thickness {
            return Some(self.right);
        }
        let value = match self.interior {
            Interior::Entire { a, b } => {
                let (cosine, sinc) = entire(self.kz.scale(depth));
                let sine_over_k = sinc.scale(depth);
                [
                    cosine * self.left[0] - C64::new(0.0, 1.0) * a * sine_over_k * self.left[1],
                    cosine * self.left[1] - C64::new(0.0, 1.0) * b * sine_over_k * self.left[0],
                ]
            }
            Interior::Decaying {
                admittance,
                forward,
                backward,
            } => {
                let phase_depth = Dd::from_f64(depth);
                let f = forward * decay(self.kz, self.phase_wavenumber, depth, phase_depth)?;
                let b = backward
                    * decay(
                        self.kz,
                        self.phase_wavenumber,
                        self.thickness - depth,
                        Dd::from_f64(self.thickness) - phase_depth,
                    )?;
                [f + b, admittance * (f - b)]
            }
        };
        value.iter().all(|&v| finite(v)).then_some(value)
    }
}

// Entire cos(z), sinc(z) series is used only for |z| <= 1/2. Twelve terms
// make the truncation far smaller than binary64 rounding, including z=0.
fn entire(z: C64) -> (C64, C64) {
    let negative_squared = -(z * z);
    let mut cos_term = C64::ONE;
    let mut sinc_term = C64::ONE;
    let mut cosine = C64::ONE;
    let mut sinc = C64::ONE;
    for n in 1..=12 {
        let twice = (2 * n) as f64;
        cos_term = (cos_term * negative_squared).scale(1.0 / (twice * (twice - 1.0)));
        sinc_term = (sinc_term * negative_squared).scale(1.0 / (twice * (twice + 1.0)));
        cosine = cosine + cos_term;
        sinc = sinc + sinc_term;
    }
    (cosine, sinc)
}

pub(super) fn decay(
    kz: C64,
    phase_wavenumber: Dd,
    distance: f64,
    phase_distance: Dd,
) -> Option<C64> {
    let amplitude = det::exp(kz.im * distance);
    if amplitude == 0.0 {
        return Some(C64::ZERO);
    }
    let phase = kz.re * distance;
    if !phase.is_finite() {
        return None;
    }
    // Preserve an endpoint difference as well as q*k before reducing phase.
    // The original rounded distance continues to own attenuation/admission.
    let (sine, cosine) = det::sin_cos_dd(phase_wavenumber * phase_distance);
    Some(C64::new(cosine.to_f64(), -sine.to_f64()).scale(amplitude))
}

#[cfg(test)]
#[path = "kernel/tests.rs"]
mod tests;
