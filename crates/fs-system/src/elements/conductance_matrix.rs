//! A linear conductance between two heat traces. Ports `a` and `b` are field heat ports: they
//! receive the temperatures at `n_a` and `n_b` points and return the heat load `Q = K T` [W] into
//! the element at each point, with `K = D^T W D`.
//!
//! `D` has `q` rows over the `n_a + n_b` points, `a`'s first: row `r` is a temperature jump across
//! the interface, a combination of the points whose weights sum to zero, so a uniform temperature
//! makes no jump and the loads sum to zero. `W` is the interface's conductance over those jumps
//! [W/K], symmetric positive definite. A contact between two matching P1 traces is
//! `D = [I, -P]` with `P` pairing each of `a`'s points with `b`'s point at the same place, and `W`
//! the face mass matrix of the patch times `1 / R''`.
//!
//! It is the many-point [`super::Conductance`]: one point on each side and `D = [1, -1]`,
//! `W = [G]` give `G (T_a - T_b)` entering at `a` and leaving at `b`.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use super::resistance::{cholesky, square_side};
use super::{ElementError, roundoff, temperature};

/// A conductance `K = D^T W D` between the field heat ports `a` and `b`.
///
/// It stores nothing and books no dissipation: what enters at the hot points leaves at the cold
/// ones. It generates `-dt (D u)^T W (D T)` of entropy per step, with `u = 1 / T` at each point,
/// the sum of `Q_i / T_i` the kernel takes over both ports.
///
/// No claim: that sum is not positive for every `D`, `W` and temperature, and the kernel
/// refuses a step in which it is negative. For a matching contact (`D = [I, -P]`) it is
/// `dt sum_rs j_r W_rs j_s / (T_a T_b)_r` over the jumps `j`: no term is negative when the jumps
/// have one sign and `W` has no negative entry. Joined points with opposed jumps at very
/// different temperatures can make it negative.
#[derive(Debug, Clone)]
pub struct ConductanceMatrix {
    widths: [usize; 2],
    rows: usize,
    /// `D`, row-major, `rows` rows of `n_a + n_b` entries.
    jump: Vec<f64>,
    /// `W`, row-major.
    weights: Vec<f64>,
    /// `sum_j |D_rj|` of each row.
    row_size: Vec<f64>,
    ports: [PortSpec; 2],
}

impl ConductanceMatrix {
    /// A conductance between a trace of `width_a` points at port `a` and one of `width_b` points
    /// at port `b`, through the row-major jumps `jump` (`q x (width_a + width_b)`) and the
    /// row-major conductance `weights` (`q x q`) [W/K].
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when a width is zero, the weights are not the square of a
    /// positive `q`, the jumps are not `q` rows of `width_a + width_b`, an entry is not finite, or
    /// a jump row's sum is not zero within `4 n eps` of the sum of its sizes;
    /// [`ElementError::Asymmetric`] when a weight and its mirror differ (a signed zero equals its
    /// mirror); [`ElementError::NotPositiveDefinite`] when the Cholesky factorization of the
    /// weights meets a pivot that is not positive and finite.
    #[expect(
        clippy::float_cmp,
        reason = "symmetry is tested exactly: an entry must equal its mirror, not lie near it"
    )]
    pub fn new(
        width_a: usize,
        width_b: usize,
        jump: Vec<f64>,
        weights: Vec<f64>,
    ) -> Result<Self, ElementError> {
        let ports = [("a", width_a), ("b", width_b)].map(|(name, width)| {
            NonZeroUsize::new(width).map(|width| {
                let id = StableId::new(name).expect("a valid stable id");
                PortSpec::heat_field(id, PortRole::EffortIn, width).expect("a valid heat port")
            })
        });
        let [Some(port_a), Some(port_b)] = ports else {
            return Err(ElementError::Parameter {
                name: "port width (points, at least 1)",
                value: width_a.min(width_b) as f64,
            });
        };
        let rows = square_side(weights.len()).ok_or(ElementError::Parameter {
            name: "weights entries (q x q, q >= 1)",
            value: weights.len() as f64,
        })?;
        let n = width_a + width_b;
        if jump.len() != rows * n {
            return Err(ElementError::Parameter {
                name: "jump entries (q x (n_a + n_b))",
                value: jump.len() as f64,
            });
        }
        if let Some(&value) = jump.iter().find(|value| !value.is_finite()) {
            return Err(ElementError::Parameter {
                name: "jump entry",
                value,
            });
        }
        if let Some(&value) = weights.iter().find(|value| !value.is_finite()) {
            return Err(ElementError::Parameter {
                name: "weights entry",
                value,
            });
        }
        let mut row_size = Vec::with_capacity(rows);
        for row in jump.chunks_exact(n) {
            let (sum, size) = row.iter().fold((0.0_f64, 0.0_f64), |(sum, size), d| {
                (sum + d, size + d.abs())
            });
            if !(sum.abs() <= 4.0 * n as f64 * f64::EPSILON * size) {
                return Err(ElementError::Parameter {
                    name: "jump row sum (zero within 4 n eps of its size)",
                    value: sum,
                });
            }
            row_size.push(size);
        }
        for row in 0..rows {
            for column in 0..row {
                if weights[row * rows + column] != weights[column * rows + row] {
                    return Err(ElementError::Asymmetric { row, column });
                }
            }
        }
        cholesky(&weights, rows)?;
        Ok(Self {
            widths: [width_a, width_b],
            rows,
            jump,
            weights,
            row_size,
            ports: [port_a, port_b],
        })
    }

    /// The points of ports `a` and `b`.
    #[must_use]
    pub fn widths(&self) -> [usize; 2] {
        self.widths
    }

    /// The number of jumps `q`.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// `K = D^T W D`, row-major over the `n_a + n_b` points, `a`'s first.
    #[must_use]
    pub fn matrix(&self) -> Vec<f64> {
        let n = self.widths[0] + self.widths[1];
        let mut column = vec![0.0; n];
        let mut matrix = vec![0.0; n * n];
        let mut unit = vec![0.0; n];
        for j in 0..n {
            unit[j] = 1.0;
            self.apply(&unit, &mut column);
            unit[j] = 0.0;
            for (i, value) in column.iter().enumerate() {
                matrix[i * n + j] = *value;
            }
        }
        matrix
    }

    /// `D x`, one jump per row.
    fn jumps(&self, values: &[f64]) -> Vec<f64> {
        let n = self.widths[0] + self.widths[1];
        self.jump
            .chunks_exact(n)
            .map(|row| row.iter().zip(values).map(|(d, x)| d * x).sum())
            .collect()
    }

    /// `W D x`, the conductance over the jumps of `values`.
    fn weighed(&self, values: &[f64]) -> Vec<f64> {
        let jumps = self.jumps(values);
        self.weights
            .chunks_exact(self.rows)
            .map(|row| row.iter().zip(&jumps).map(|(w, j)| w * j).sum())
            .collect()
    }

    /// `D^T y` into `out`.
    fn spread(&self, loads: &[f64], out: &mut [f64]) {
        let n = self.widths[0] + self.widths[1];
        out[..n].fill(0.0);
        for (row, &load) in self.jump.chunks_exact(n).zip(loads) {
            for (value, d) in out[..n].iter_mut().zip(row) {
                *value = d.mul_add(load, *value);
            }
        }
    }

    /// `K x` into `out`.
    fn apply(&self, values: &[f64], out: &mut [f64]) {
        let loads = self.weighed(values);
        self.spread(&loads, out);
    }
}

impl Component for ConductanceMatrix {
    type State = ();
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::ConductanceMatrix",
        state_family: "stateless",
        schema_grammar: "(): a memoryless conductance between two traces stores nothing",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        // The widths are part of the system's port description; the number of jumps is not.
        words.push(self.rows as u64);
        words.extend(self.jump.iter().map(|value| value.to_bits()));
        words.extend(self.weights.iter().map(|value| value.to_bits()));
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
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
        let n = self.widths[0] + self.widths[1];
        for &kelvin in &inputs[..n] {
            temperature(kelvin)?;
        }
        let loads = self.weighed(&inputs[..n]);
        self.spread(&loads, outputs);
        // The entropy into the element, dt sum Q_i / T_i, is minus what it generates. Summed
        // point by point it cancels between the two sides; `(D u)^T (W D T)` is the same sum
        // over the jumps, whose terms in a matching contact are each `1 / T_a - 1 / T_b` times
        // the conductance over `T_a - T_b`, and so of one sign.
        let inverse: Vec<f64> = inputs[..n].iter().map(|kelvin| 1.0 / kelvin).collect();
        let inverse_jumps = self.jumps(&inverse);
        let generated = -dt_s
            * inverse_jumps
                .iter()
                .zip(&loads)
                .map(|(u, y)| u * y)
                .sum::<f64>();
        // The loads sum to zero because each row of D does: the sum of `D^T y` errs by a few
        // `n eps` of `sum_r |y_r| sum_j |D_rj|`, and the row sums add at most `4 n eps` of the
        // same (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., sec. 3.1).
        let scale: f64 = loads
            .iter()
            .zip(&self.row_size)
            .map(|(y, size)| y.abs() * size)
            .sum();
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: generated,
            energy_tolerance_j: roundoff((n + self.rows + 2) as f64 * dt_s * scale),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for ConductanceMatrix {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let n = self.widths[0] + self.widths[1];
        self.apply(&d_inputs[..n], d_outputs);
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // K = D^T W D with W symmetric, so K^T bar is K bar.
        let n = self.widths[0] + self.widths[1];
        self.apply(&bar_outputs[..n], bar_inputs);
        Ok(())
    }
}
