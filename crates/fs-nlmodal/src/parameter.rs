//! Original coefficient derivatives of the quartic sum-of-squares effort law.

use super::{NlModalError, SYM_TOL, SosModalStorage};

/// One original stress channel's derivative per caller-owned real parameter.
/// Signed directions need not be positive: they are local derivatives, not a
/// new admitted storage. None holds the complete coupling matrix fixed.
#[derive(Clone, Copy, Debug, Default)]
pub struct StressChannelParameterDirection<'a> {
    pub coefficient: f64,
    /// Original row-major symmetric mode-by-mode coupling derivative.
    pub coupling: Option<&'a [f64]>,
}

/// A complete local coefficient direction in the ORIGINAL modal state chart.
/// None holds that whole field fixed; supplied arrays cover every original
/// mode/channel. Frequencies are angular frequencies, not squared stiffnesses.
/// The mass-normalized coordinates, modal basis and momentum law stay fixed.
/// A changing basis or mass normalization needs its own coordinate chain rule.
#[derive(Clone, Copy, Debug, Default)]
pub struct SosModalParameterDirection<'a> {
    pub omegas: Option<&'a [f64]>,
    pub channels: Option<&'a [StressChannelParameterDirection<'a>]>,
}

fn invalid(what: &'static str) -> NlModalError {
    NlModalError::Parameter { what }
}

impl SosModalStorage {
    /// Differentiate the original effort at FIXED physical state.
    ///
    /// For each channel, `s=q^T E q`, `e=E q` and the effort is `c s e`.
    /// Its parameter derivative is `dc s e + c (q^T dE q) e + c s dE q`.
    /// All signed/off-diagonal channels remain; a zero baseline coefficient
    /// still contributes when its direction is nonzero. Linear rows add
    /// `2 omega domega q`, while momentum derivatives are zero.
    ///
    /// This is distinct from [`Self::hessian_action`], which differentiates
    /// the state. A total response derivative needs BOTH actions. The result
    /// has original effort units per declared parameter unit. No perturbed
    /// model, state, modal reduction or finite-difference step is formed.
    /// Two matrix-vector passes per channel use O(channels*modes^2) work and
    /// constant temporary memory. Source admission is rechecked because its
    /// fields are public. Roundoff-sized admitted asymmetry is differentiated
    /// as originally represented rather than symmetrized.
    ///
    /// # Errors
    /// Invalid source/direction or incomplete/nonfinite state/output shape
    /// refuses before writing output. Unrepresentable arithmetic refuses with
    /// unspecified scratch output. No finite parameter interval is admitted.
    pub fn gradient_parameter_derivative(
        &self,
        state: &[f64],
        direction: SosModalParameterDirection<'_>,
        out: &mut [f64],
    ) -> Result<(), NlModalError> {
        self.validate()?;
        let n = self.omegas.len();
        let dimension = n
            .checked_mul(2)
            .ok_or_else(|| invalid("modal parameter dimension overflow"))?;
        if state.len() != dimension
            || out.len() != dimension
            || state.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "modal parameter derivative needs complete finite state and matching output",
            ));
        }
        if let Some(omegas) = direction.omegas {
            if omegas.len() != n || omegas.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "complete finite modal frequency direction required",
                ));
            }
        }
        if let Some(channels) = direction.channels {
            if channels.len() != self.channels.len() {
                return Err(invalid(
                    "direction must cover every original stress channel",
                ));
            }
            for channel in channels {
                if !channel.coefficient.is_finite() {
                    return Err(invalid("finite stress coefficient direction required"));
                }
                if let Some(matrix) = channel.coupling {
                    if matrix.len() != n * n || matrix.iter().any(|v| !v.is_finite()) {
                        return Err(invalid(
                            "complete finite stress coupling direction required",
                        ));
                    }
                    let scale = matrix
                        .iter()
                        .fold(f64::MIN_POSITIVE, |a, &v| a.max(v.abs()));
                    let mut worst = 0.0_f64;
                    for i in 0..n {
                        for j in 0..i {
                            worst = worst.max((matrix[i * n + j] - matrix[j * n + i]).abs());
                        }
                    }
                    if worst > SYM_TOL * scale {
                        return Err(NlModalError::AsymmetricCoupling {
                            residual: worst / scale,
                        });
                    }
                }
            }
        }
        out.fill(0.0);
        if let Some(omegas) = direction.omegas {
            for i in 0..n {
                out[2 * i] = 2.0 * self.omegas[i] * omegas[i] * state[2 * i];
            }
        }
        for (index, channel) in self.channels.iter().enumerate() {
            let Some(slope) = direction.channels.map(|channels| channels[index]) else {
                continue;
            };
            if slope.coefficient == 0.0 && slope.coupling.is_none() {
                continue;
            }
            let mut form = 0.0;
            let mut form_derivative = 0.0;
            for i in 0..n {
                let mut eq = 0.0;
                let mut deq = 0.0;
                for j in 0..n {
                    eq += channel.coupling[i * n + j] * state[2 * j];
                    deq += slope.coupling.map_or(0.0, |matrix| matrix[i * n + j]) * state[2 * j];
                }
                form += state[2 * i] * eq;
                form_derivative += state[2 * i] * deq;
            }
            if !form.is_finite() || !form_derivative.is_finite() {
                return Err(invalid("modal parameter channel contraction overflow"));
            }
            for i in 0..n {
                let mut eq = 0.0;
                let mut deq = 0.0;
                for j in 0..n {
                    eq += channel.coupling[i * n + j] * state[2 * j];
                    deq += slope.coupling.map_or(0.0, |matrix| matrix[i * n + j]) * state[2 * j];
                }
                out[2 * i] += slope.coefficient * form * eq
                    + channel.coefficient * (form_derivative * eq + form * deq);
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid("modal parameter effort derivative overflow"));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StressChannel;
    use fs_phs::Storage;

    #[test]
    fn complete_signed_channel_direction_matches_converging_original_gradients() {
        let source = SosModalStorage {
            omegas: vec![2.0, 3.0],
            channels: vec![
                StressChannel {
                    coefficient: 0.8,
                    coupling: vec![1.0, -0.2, -0.2, 0.3],
                },
                StressChannel {
                    coefficient: 0.2,
                    coupling: vec![-0.4, 0.8, 0.8, 1.0],
                },
            ],
        };
        let omegas = [0.3, -0.7];
        let first = [0.2, 0.3, 0.3, -0.1];
        let second = [-0.4, -0.1, -0.1, 0.5];
        let channels = [
            StressChannelParameterDirection {
                coefficient: -0.4,
                coupling: Some(&first),
            },
            StressChannelParameterDirection {
                coefficient: 0.3,
                coupling: Some(&second),
            },
        ];
        let state = [0.4, -0.2, -0.3, 0.7];
        let direction = SosModalParameterDirection {
            omegas: Some(&omegas),
            channels: Some(&channels),
        };
        let mut actual = [0.0; 4];
        source
            .gradient_parameter_derivative(&state, direction, &mut actual)
            .unwrap();
        assert_eq!(actual[1], 0.0);
        assert_eq!(actual[3], 0.0);
        let mut errors = Vec::new();
        for step in [0.04, 0.02, 0.01] {
            let gradient = |delta: f64| {
                let perturbed = SosModalStorage {
                    omegas: source
                        .omegas
                        .iter()
                        .zip(omegas)
                        .map(|(w, d)| w + delta * d)
                        .collect(),
                    channels: source
                        .channels
                        .iter()
                        .zip(channels)
                        .map(|(channel, direction)| StressChannel {
                            coefficient: channel.coefficient + delta * direction.coefficient,
                            coupling: channel
                                .coupling
                                .iter()
                                .zip(direction.coupling.unwrap())
                                .map(|(e, de)| e + delta * de)
                                .collect(),
                        })
                        .collect(),
                };
                let mut out = [0.0; 4];
                perturbed.gradient(&state, &mut out);
                out
            };
            let plus = gradient(step);
            let minus = gradient(-step);
            let error = actual
                .iter()
                .zip(plus)
                .zip(minus)
                .map(|((a, p), m)| (a - (p - m) / (2.0 * step)).abs())
                .fold(0.0, f64::max);
            errors.push(error);
        }
        assert!(errors[0] > 3.9 * errors[1]);
        assert!(errors[1] > 3.9 * errors[2]);
        assert!(errors[2] < 1e-5);
    }

    #[test]
    fn zero_baseline_channel_and_linear_small_amplitude_limits_remain_physical() {
        let source = SosModalStorage {
            omegas: vec![3.0],
            channels: vec![StressChannel {
                coefficient: 0.0,
                coupling: vec![-2.0],
            }],
        };
        let channels = [StressChannelParameterDirection {
            coefficient: 0.5,
            coupling: None,
        }];
        let mut out = [0.0; 2];
        for amplitude in [0.4_f64, 0.2, 1e-3] {
            source
                .gradient_parameter_derivative(
                    &[amplitude, 0.7],
                    SosModalParameterDirection {
                        omegas: None,
                        channels: Some(&channels),
                    },
                    &mut out,
                )
                .unwrap();
            assert!((out[0] - 2.0 * amplitude.powi(3)).abs() < 1e-15);
            assert_eq!(out[1], 0.0);
        }
        source
            .gradient_parameter_derivative(
                &[0.4, 0.7],
                SosModalParameterDirection {
                    omegas: Some(&[0.2]),
                    channels: None,
                },
                &mut out,
            )
            .unwrap();
        assert!((out[0] - 0.48).abs() < 1e-15);
    }

    #[test]
    fn malformed_direction_refuses_before_writing_output() {
        let source = SosModalStorage {
            omegas: vec![2.0, 3.0],
            channels: vec![StressChannel {
                coefficient: 1.0,
                coupling: vec![1.0, 0.0, 0.0, 1.0],
            }],
        };
        let bad = [StressChannelParameterDirection {
            coefficient: 0.0,
            coupling: Some(&[1.0, 0.2, 0.0, 1.0]),
        }];
        let mut out = [7.0; 4];
        for direction in [
            SosModalParameterDirection {
                omegas: Some(&[f64::NAN, 0.0]),
                channels: None,
            },
            SosModalParameterDirection {
                omegas: None,
                channels: Some(&[]),
            },
            SosModalParameterDirection {
                omegas: None,
                channels: Some(&bad),
            },
        ] {
            assert!(
                source
                    .gradient_parameter_derivative(&[0.3; 4], direction, &mut out)
                    .is_err()
            );
            assert_eq!(out, [7.0; 4]);
        }
    }
}
