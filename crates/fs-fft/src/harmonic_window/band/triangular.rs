//! Orthogonal field projection for normalized triangular sensitivity.
//!
//! On x in [-1,1], dmu=(1-|x|)dx has moments mu_(2k)=2/((2k+1)(2k+2)).
//! Its monic polynomials satisfy x*P_m=P_(m+1)+beta_m*P_(m-1), with
//! beta_m=h_m/h_(m-1), h_m=integral(P_m^2 dmu). Symmetry removes the
//! diagonal recurrence coefficient. See DLMF 18.2.11_4 and 18.2.11_5:
//! https://dlmf.nist.gov/18.2.E11_4 and https://dlmf.nist.gov/18.2.E11_5.
//!
//! Expanding x^n in P_m gives positive coefficients C_(n,m). Thus the
//! exponential's mth coefficient, with its common (-i)^m omitted, is
//! a^m * sum_k (-a^2)^k C_(m+2k,m)/(m+2k)!.
//! Scale by sqrt(h_m/(2m+1)) so every existing positive modal contraction
//! keeps its weight 2m+1. Mode zero is exactly the weighted field mean.
//! This projects fields before products, never a matrix of bright pair powers.
//!
//! The exact triangular projection minimizes its weighted residual. Its tail
//! is therefore no larger than the uniform Legendre residual under dmu, at
//! most twice that uniform tail because dmu <= 2*(dx/2). The existing mode
//! admission gives omitted power below 1e-49*(sum |A_i|)^2. Twenty Taylor
//! terms per mode retain at least degree 39; at |a|<=0.5 the exponential
//! remainder bound is below 1.84e-60. Neither bound covers floating arithmetic,
//! source phase uncertainty, or arbitrary-null relative accuracy.

use super::{Dd, MAX_UNIFORM_HARMONIC_BAND_MODES};

const TERMS: usize = 20;
const HEIGHT: usize = MAX_UNIFORM_HARMONIC_BAND_MODES + TERMS - 1;

// Basis construction keeps its result and two bounded recurrence vectors.
// Count the mode metadata as well as the 6,864-byte numerical payload.
pub(super) const PREPARATION_SCRATCH_BYTES: usize =
    core::mem::size_of::<Basis>() + 2 * core::mem::size_of::<[Dd; HEIGHT]>();

/// beta[0] is unused. Constants are DD roundings of exact rationals obtained
/// from the displayed moments: P0=1, P1=x; P_(n+1)=x*P_n-beta_n*P_(n-1);
/// h_n=sum_(i,j) P_n[i]*P_n[j]*mu_(i+j), beta_n=h_n/h_(n-1).
/// No fitted weight, quadrature or numerical Gram factorization is involved.
const BETA: [Dd; HEIGHT] = [
    Dd::ZERO,
    Dd {
        hi: 0.16666666666666666,
        lo: 9.25185853854297e-18,
    },
    Dd {
        hi: 0.23333333333333334,
        lo: -3.7007434154171884e-18,
    },
    Dd {
        hi: 0.23265306122448978,
        lo: 1.3707855712208565e-17,
    },
    Dd {
        hi: 0.24453992123165055,
        lo: 1.6794458571576392e-18,
    },
    Dd {
        hi: 0.24253282614228905,
        lo: -1.0383687618373551e-17,
    },
    Dd {
        hi: 0.24734048652178317,
        lo: -4.396195522435889e-19,
    },
    Dd {
        hi: 0.24584455698969224,
        lo: 9.203705389845228e-18,
    },
    Dd {
        hi: 0.24843705244991143,
        lo: 9.510321725024748e-18,
    },
    Dd {
        hi: 0.24735248607766946,
        lo: -1.374037523751846e-17,
    },
    Dd {
        hi: 0.2489746603424117,
        lo: 1.311155788748039e-17,
    },
    Dd {
        hi: 0.2481652625939299,
        lo: -7.707576643889179e-18,
    },
    Dd {
        hi: 0.249276872138724,
        lo: 1.0060389398784347e-17,
    },
    Dd {
        hi: 0.24865332472773147,
        lo: -5.058914532921021e-18,
    },
    Dd {
        hi: 0.24946320940416367,
        lo: -1.3483985148946863e-17,
    },
    Dd {
        hi: 0.24896936222417373,
        lo: 1.1266049111302432e-17,
    },
    Dd {
        hi: 0.24958603295918674,
        lo: -5.4412576604934134e-18,
    },
    Dd {
        hi: 0.2491857471865592,
        lo: -4.23388513691319e-18,
    },
    Dd {
        hi: 0.2496711872300111,
        lo: 7.347283934868221e-18,
    },
    Dd {
        hi: 0.2493404078326574,
        lo: -1.2760751759088584e-17,
    },
    Dd {
        hi: 0.24973260927953286,
        lo: 1.147230907709552e-17,
    },
    Dd {
        hi: 0.2494547917176547,
        lo: 1.308621006777582e-17,
    },
    Dd {
        hi: 0.24977834763425336,
        lo: -2.8036559875648746e-18,
    },
    Dd {
        hi: 0.2495417731015751,
        lo: -1.3115935150027751e-17,
    },
    Dd {
        hi: 0.24981331126551504,
        lo: 9.826910093212699e-18,
    },
    Dd {
        hi: 0.24960946154104752,
        lo: 4.003392554719851e-19,
    },
    Dd {
        hi: 0.24984063097865555,
        lo: -7.271839768681094e-18,
    },
    Dd {
        hi: 0.2496631724749628,
        lo: 8.50881403251657e-18,
    },
    Dd {
        hi: 0.24986237918274679,
        lo: 1.049606272390046e-18,
    },
    Dd {
        hi: 0.2497065081991481,
        lo: 6.211917040048296e-18,
    },
    Dd {
        hi: 0.24987997155873026,
        lo: 9.291811446455081e-18,
    },
    Dd {
        hi: 0.24974198043531665,
        lo: -1.123024992931898e-17,
    },
    Dd {
        hi: 0.24989440162911028,
        lo: -5.066101920802153e-18,
    },
    Dd {
        hi: 0.24977138398661322,
        lo: -9.241822941624209e-18,
    },
    Dd {
        hi: 0.24990638315895625,
        lo: -2.6294771942445705e-18,
    },
    Dd {
        hi: 0.24979602908120746,
        lo: 2.9142901324343603e-18,
    },
];

pub(super) struct Basis {
    series: [[Dd; TERMS]; MAX_UNIFORM_HARMONIC_BAND_MODES],
    normalization: [Dd; MAX_UNIFORM_HARMONIC_BAND_MODES],
    modes: usize,
}

impl Basis {
    pub(super) fn new(modes: usize) -> Self {
        let mut result = Self {
            series: [[Dd::ZERO; TERMS]; MAX_UNIFORM_HARMONIC_BAND_MODES],
            normalization: [Dd::ZERO; MAX_UNIFORM_HARMONIC_BAND_MODES],
            modes,
        };
        let degree = modes - 1 + 2 * (TERMS - 1);
        let height = modes - 1 + TERMS;
        let mut power = [Dd::ZERO; HEIGHT];
        power[0] = Dd::ONE;
        let mut factorial = Dd::ONE;
        for n in 0..=degree {
            for m in (n % 2..modes.min(n + 1)).step_by(2) {
                let k = (n - m) / 2;
                if k < TERMS {
                    result.series[m][k] = power[m] / factorial;
                }
            }
            if n == degree {
                break;
            }
            let mut next = [Dd::ZERO; HEIGHT];
            for m in (n % 2..height.min(n + 1)).step_by(2) {
                if m + 1 < height {
                    next[m + 1] = next[m + 1] + power[m];
                }
                if m != 0 {
                    next[m - 1] = next[m - 1] + BETA[m] * power[m];
                }
            }
            power = next;
            factorial = factorial * Dd::from_f64((n + 1) as f64);
        }
        let mut norm = Dd::ONE;
        for m in 0..modes {
            if m != 0 {
                norm = norm * BETA[m];
            }
            result.normalization[m] = (norm / Dd::from_f64((2 * m + 1) as f64)).sqrt();
        }
        result
    }

    pub(super) fn fill(&self, a: Dd, source: usize, stride: usize, output: &mut [Dd]) {
        let negative_square = -(a * a);
        let mut leading = Dd::ONE;
        for m in 0..self.modes {
            let mut series = self.series[m][TERMS - 1];
            for k in (0..TERMS - 1).rev() {
                series = series * negative_square + self.series[m][k];
            }
            output[m * stride + source] = (leading * self.normalization[m]) * series;
            leading = leading * a;
        }
    }
}
