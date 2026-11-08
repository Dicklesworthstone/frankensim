//! Factored frequency actions shared by the first and second response jets.
//! This numerical kernel takes a cancellation callback, not an executor.

use super::super::{HarmonicResponseError, MatrixLayout, finite, invalid, scale, scale_complex};
use super::add;
use crate::quadratic_fraction;
use fs_math::c64::C64;

pub(super) struct FrequencyAction {
    /// Minus the scaled pencil derivative applied to the supplied state,
    /// divided by 2^exponent. Derivatives are with respect to scaled frequency.
    pub(super) value: Vec<C64>,
    pub(super) envelope: Vec<f64>,
    pub(super) exponent: i32,
}

#[derive(Clone, Copy)]
pub(super) enum FrequencyOrder<'a> {
    First(Option<&'a [C64]>),
    Second {
        slope: Option<&'a [C64]>,
        curvature: Option<&'a [C64]>,
    },
}

pub(super) struct FrequencyPort<'a> {
    pub(super) ports: usize,
    pub(super) force_map: &'a [f64],
    pub(super) impedance: &'a [C64],
    pub(super) exponents: &'a [i32],
    pub(super) minimum_time_exponent: i32,
}

pub(super) struct FrequencyPencil<'a> {
    pub(super) n: usize,
    pub(super) layout: MatrixLayout,
    pub(super) scaled_frequency: f64,
    pub(super) time_exponent: i32,
    pub(super) coordinate_exponents: &'a [i32],
    pub(super) scaled_mass: &'a [f64],
    pub(super) scaled_damping: &'a [f64],
    pub(super) port: Option<FrequencyPort<'a>>,
}

impl FrequencyPencil<'_> {
    /// Apply -A_mu or -A_mumu without forming an n-by-n port operator.
    /// The parent has admitted all source extents, domains and work budgets.
    pub(super) fn apply(
        &self,
        state: &[C64],
        order: FrequencyOrder<'_>,
        mut poll: impl FnMut() -> Result<(), HarmonicResponseError>,
    ) -> Result<FrequencyAction, HarmonicResponseError> {
        let mut exponent = 0;
        let mu = self.scaled_frequency;
        let (mu_fraction, mu_exponent) = if mu == 0.0 {
            (0.0, 0)
        } else {
            quadratic_fraction(mu)
        };
        if let Some(load) = &self.port {
            let headroom = 2 * load.ports.next_power_of_two().ilog2() as i32 + 2;
            match order {
                FrequencyOrder::First(Some(slope)) => {
                    if load.minimum_time_exponent != i32::MIN {
                        exponent = exponent.max(load.minimum_time_exponent - self.time_exponent);
                    }
                    if mu != 0.0 {
                        for a in 0..load.ports {
                            poll()?;
                            for b in 0..load.ports {
                                if load.exponents[a] == i32::MIN || load.exponents[b] == i32::MIN {
                                    continue;
                                }
                                for value in
                                    [slope[a * load.ports + b].re, slope[a * load.ports + b].im]
                                {
                                    if value != 0.0 {
                                        exponent = exponent.max(
                                            load.exponents[a]
                                                + load.exponents[b]
                                                + quadratic_fraction(value).1
                                                + mu_exponent
                                                + headroom,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
                FrequencyOrder::Second {
                    slope: Some(slope),
                    curvature: Some(curvature),
                } => {
                    for a in 0..load.ports {
                        poll()?;
                        for b in 0..load.ports {
                            if load.exponents[a] == i32::MIN || load.exponents[b] == i32::MIN {
                                continue;
                            }
                            let at = a * load.ports + b;
                            for (value, shift, active) in [
                                (slope[at].re, 1, true),
                                (slope[at].im, 1, true),
                                (
                                    curvature[at].re,
                                    mu_exponent + self.time_exponent,
                                    mu != 0.0,
                                ),
                                (
                                    curvature[at].im,
                                    mu_exponent + self.time_exponent,
                                    mu != 0.0,
                                ),
                            ] {
                                if value != 0.0 && active {
                                    exponent = exponent.max(
                                        load.exponents[a]
                                            + load.exponents[b]
                                            + quadratic_fraction(value).1
                                            + shift
                                            + headroom,
                                    );
                                }
                            }
                        }
                    }
                }
                _ => return Err(invalid("complete load frequency derivatives required")),
            }
        }
        let mut value = vec![C64::ZERO; self.n];
        let mut envelope = vec![0.0; self.n];
        for i in 0..self.n {
            poll()?;
            let mut correction = C64::ZERO;
            for j in self.layout.columns(i, self.n) {
                let q = state[j];
                if q == C64::ZERO {
                    continue;
                }
                let at = self.layout.index(i, j, self.n);
                let coefficient = match order {
                    FrequencyOrder::First(_) => {
                        let mass = if mu == 0.0 {
                            0.0
                        } else {
                            2.0 * mu * scale(self.scaled_mass[at], -exponent)?
                        };
                        C64::new(mass, scale(self.scaled_damping[at], -exponent)?)
                    }
                    FrequencyOrder::Second { .. } => {
                        C64::from_re(scale(self.scaled_mass[at], 1 - exponent)?)
                    }
                };
                add(&mut value[i], &mut correction, coefficient * q);
                envelope[i] += coefficient.abs() * q.abs();
            }
        }
        if let Some(load) = &self.port {
            let p = load.ports;
            let mut b = vec![0.0; self.n * p];
            for i in 0..self.n {
                poll()?;
                for a in 0..p {
                    if load.exponents[a] != i32::MIN {
                        b[i * p + a] = scale(
                            load.force_map[i * p + a],
                            -self.coordinate_exponents[i] - load.exponents[a],
                        )?;
                    }
                }
            }
            let mut flow = vec![C64::ZERO; p];
            let mut flow_envelope = vec![0.0; p];
            for a in 0..p {
                poll()?;
                let mut correction = C64::ZERO;
                for (i, &q) in state.iter().enumerate() {
                    add(&mut flow[a], &mut correction, q.scale(b[i * p + a]));
                    flow_envelope[a] += q.abs() * b[i * p + a].abs();
                }
            }
            let mut effort = vec![C64::ZERO; p];
            let mut effort_envelope = vec![0.0; p];
            for a in 0..p {
                poll()?;
                if load.exponents[a] == i32::MIN {
                    continue;
                }
                let mut correction = C64::ZERO;
                for c in 0..p {
                    if load.exponents[c] == i32::MIN {
                        continue;
                    }
                    let common = load.exponents[a] + load.exponents[c] - exponent;
                    let (first, second) = match order {
                        FrequencyOrder::First(Some(slope)) => (
                            scale_complex(load.impedance[a * p + c], common - self.time_exponent)?,
                            if mu == 0.0 {
                                C64::ZERO
                            } else {
                                scale_complex(slope[a * p + c], common + mu_exponent)?
                                    .scale(mu_fraction)
                            },
                        ),
                        FrequencyOrder::Second {
                            slope: Some(slope),
                            curvature: Some(curvature),
                        } => (
                            scale_complex(slope[a * p + c], common + 1)?,
                            if mu == 0.0 {
                                C64::ZERO
                            } else {
                                scale_complex(
                                    curvature[a * p + c],
                                    common + mu_exponent + self.time_exponent,
                                )?
                                .scale(mu_fraction)
                            },
                        ),
                        _ => return Err(invalid("complete load frequency derivatives required")),
                    };
                    add(&mut effort[a], &mut correction, first * flow[c]);
                    add(&mut effort[a], &mut correction, second * flow[c]);
                    effort_envelope[a] += (first.abs() + second.abs()) * flow_envelope[c];
                }
            }
            for i in 0..self.n {
                poll()?;
                let mut correction = C64::ZERO;
                for a in 0..p {
                    let term = effort[a].scale(b[i * p + a]);
                    add(&mut value[i], &mut correction, C64::new(-term.im, term.re));
                    envelope[i] += effort_envelope[a] * b[i * p + a].abs();
                }
            }
        }
        if value.iter().any(|&z| !finite(z)) || envelope.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete frequency action is not representable"));
        }
        Ok(FrequencyAction {
            value,
            envelope,
            exponent,
        })
    }
}

#[cfg(test)]
#[path = "action_tests.rs"]
mod tests;
