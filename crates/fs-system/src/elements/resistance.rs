//! A memoryless resistor N-port: through one power port of width `n` it receives the flows `u`
//! and returns the efforts `R u`, with `R` symmetric positive definite. All the power it absorbs,
//! `u^T R u`, leaves as heat through port `heat` at the temperature that port receives.
//!
//! It is the many-port [`super::Damper`]: a matrix of damping coefficients between `n`
//! coordinates of one power domain, such as the self and mutual radiation resistance of a set of
//! modes. A diagonal `R` gives each row the damper's effort `r u` exactly.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};
use fs_ports::{PortKind, StableId};

use super::{ElementError, heat, inputs_finite, roundoff, temperature};

/// A linear resistor N-port with port `u` (flow in, `n` rows) and port `heat`.
///
/// `R` must be positive definite, not only semidefinite: a passive device whose coordinates each
/// carry some loss of their own has a definite `R`, and the definite factor `R = L L^T` lets the
/// absorbed power be computed as `|L^T u|^2`, which is never negative, so the generated entropy
/// is never negative either.
#[derive(Debug, Clone)]
pub struct Resistance {
    n: usize,
    /// `R`, row-major.
    matrix: Vec<f64>,
    /// The lower Cholesky factor of `R`, row-major.
    factor: Vec<f64>,
    ports: [PortSpec; 2],
}

impl Resistance {
    /// A resistor N-port in `kind` with the row-major `n x n` matrix `matrix`, effort per unit flow
    /// in that domain.
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when the entry count is not the square of a positive `n`, or an
    /// entry is not finite; [`ElementError::Asymmetric`] when an entry and its mirror differ (a
    /// signed zero equals its mirror); [`ElementError::NotPositiveDefinite`] when the Cholesky
    /// factorization meets a pivot that is not positive and finite.
    #[expect(
        clippy::float_cmp,
        reason = "symmetry is tested exactly: an entry must equal its mirror, not lie near it"
    )]
    pub fn new(kind: PortKind, matrix: Vec<f64>) -> Result<Self, ElementError> {
        let n = square_side(matrix.len()).ok_or(ElementError::Parameter {
            name: "matrix entries (n x n, n >= 1)",
            value: matrix.len() as f64,
        })?;
        if let Some(&value) = matrix.iter().find(|value| !value.is_finite()) {
            return Err(ElementError::Parameter {
                name: "matrix entry",
                value,
            });
        }
        for row in 0..n {
            for column in 0..row {
                if matrix[row * n + column] != matrix[column * n + row] {
                    return Err(ElementError::Asymmetric { row, column });
                }
            }
        }
        let factor = cholesky(&matrix, n)?;
        let width = NonZeroUsize::new(n).expect("n is at least 1");
        let flow = PortSpec::power(
            StableId::new("u").expect("a valid stable id"),
            kind,
            PortRole::FlowIn,
            width,
        )
        .expect("a valid power port");
        Ok(Self {
            n,
            matrix,
            factor,
            ports: [flow, heat("heat", PortRole::EffortIn)],
        })
    }

    /// The number of rows `n`.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.n
    }

    /// `R`, row-major.
    #[must_use]
    pub fn matrix(&self) -> &[f64] {
        &self.matrix
    }

    /// `efforts = R flows`, row by row in index order.
    fn apply(&self, flows: &[f64], efforts: &mut [f64]) {
        for (row, effort) in self.matrix.chunks_exact(self.n).zip(efforts.iter_mut()) {
            *effort = row.iter().zip(flows).map(|(r, u)| r * u).sum();
        }
    }

    /// The absorbed power `|L^T u|^2` [W], never negative.
    fn power_w(&self, flows: &[f64]) -> f64 {
        let n = self.n;
        (0..n)
            .map(|j| {
                let y: f64 = (j..n).map(|i| self.factor[i * n + j] * flows[i]).sum();
                y * y
            })
            .sum()
    }

    /// The size of the energy balance's terms: `sum |R_ij u_i u_j|` plus `(sum sqrt(R_ii) |u_i|)^2`,
    /// which bounds `|L| |L|^T` against `|u|` (each row of `L` has the norm `sqrt(R_ii)`).
    fn balance_scale(&self, flows: &[f64]) -> f64 {
        let (mut products, mut diagonal) = (0.0, 0.0);
        for (i, row) in self.matrix.chunks_exact(self.n).enumerate() {
            for (r, v) in row.iter().zip(flows) {
                products += (r * flows[i] * v).abs();
            }
            diagonal += row[i].sqrt() * flows[i].abs();
        }
        products + diagonal * diagonal
    }
}

/// `n` when `len` is `n * n` with `n >= 1`.
fn square_side(len: usize) -> Option<usize> {
    let n = len.isqrt();
    (n >= 1 && n * n == len).then_some(n)
}

/// The lower Cholesky factor of the symmetric `m`, or the first pivot that is not positive and
/// finite (a sum of finite entries can still overflow).
fn cholesky(m: &[f64], n: usize) -> Result<Vec<f64>, ElementError> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = m[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if !(sum > 0.0 && sum.is_finite()) {
                    return Err(ElementError::NotPositiveDefinite {
                        pivot: i,
                        value: sum,
                    });
                }
                l[i * n + i] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    Ok(l)
}

impl Component for Resistance {
    type State = ();
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        let n = self.n;
        let flows = &inputs[..n];
        let kelvin = temperature(inputs[n])?;
        self.apply(flows, &mut outputs[..n]);
        let loss_w = self.power_w(flows);
        outputs[n] = -loss_w;
        let dissipated_j = dt_s * loss_w;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dissipated_j / kelvin,
            // The work in, `dt u . R u`, and the heat out, `dt |L^T u|^2`, are one quantity summed
            // two ways; each errs by a few `n eps` times the balance scale (Higham, Accuracy and
            // Stability of Numerical Algorithms, 2nd ed., Thm 10.3, for the factor).
            energy_tolerance_j: roundoff((n as f64 + 2.0) * dt_s * self.balance_scale(flows)),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for Resistance {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.n;
        self.apply(&d_inputs[..n], &mut d_outputs[..n]);
        // The heat rate is -u^T R u; its derivative is -2 (R u) . du.
        let mut pushed = vec![0.0; n];
        self.apply(&inputs[..n], &mut pushed);
        d_outputs[n] = -2.0
            * pushed
                .iter()
                .zip(&d_inputs[..n])
                .map(|(a, b)| a * b)
                .sum::<f64>();
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.n;
        // R is symmetric, so R^T bar_e is R bar_e.
        self.apply(&bar_outputs[..n], &mut bar_inputs[..n]);
        let mut pushed = vec![0.0; n];
        self.apply(&inputs[..n], &mut pushed);
        for (bar, along) in bar_inputs[..n].iter_mut().zip(&pushed) {
            *bar -= 2.0 * along * bar_outputs[n];
        }
        bar_inputs[n] = 0.0;
        Ok(())
    }
}
