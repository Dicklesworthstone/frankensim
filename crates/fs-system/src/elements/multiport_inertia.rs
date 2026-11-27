//! A multi-port inertia: the many-port [`super::Mass`]. Its state is `n` generalized velocities
//! `v`; its momentum is `p = M v` with a constant symmetric positive definite mass matrix `M`, and
//! its kinetic energy `v^T M v / 2 = p^T M^-1 p / 2`. Each of its `m` ports `p0`, `p1`, ...
//! receives an effort and returns a flow; a constant port map `G` (`n x m`) takes the efforts `e`
//! to the generalized force `G e` and the velocities to the port flows `G^T v`.
//!
//! The units are the caller's. A row of `v` is in the port flow's unit over its entries of `G`:
//! a row whose entries carry a lever arm (metres) is an angular velocity, and `M` carries the
//! matching mixed units. So the state declares no dimension (`state_dimensions` is `None`).
//!
//! The step is the mass's in matrix form, the implicit midpoint rule under held efforts:
//! `M (v1 - v0) = dt G e`, solved with the square-root-free Cholesky factor `M = L D L^T`, and each
//! port returns the step-mean flow `G^T (v0 + v1) / 2`. A port's flow is linear in `v`, so it ends
//! its step by the midpoint relation, and a load line (a contact, a slider) sits beside it as it
//! sits beside a `Mass`. The port work `dt e . G^T (v0 + v1) / 2` is the kinetic energy's change
//! `(v0 + v1)^T M (v1 - v0) / 2`, so the receipt is lossless. One row and one port with `G = [1]`
//! give the `Mass` step to the bit.
//!
//! No claim: `M` and `G` are constants. That is exact for a body whose configuration does not
//! change over the window it is used for, such as a rigid body held at one pose through a short
//! impact or slip. It is not a rigid body that rotates, whose mass matrix and port map move with
//! its pose. It keeps no positions and offers no staged residual.

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::PortKind;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use super::{ElementError, inputs_finite, port, roundoff};

/// An inertia of `n` generalized velocities seen through `m` effort-input ports `p0` to
/// `p<m-1>`, with a constant mass matrix and a constant port map (see the module's no-claim).
#[derive(Debug, Clone)]
pub struct MultiportInertia {
    n: usize,
    /// `M`, row-major `n x n`.
    mass: Vec<f64>,
    /// The unit lower factor `L` of `M = L D L^T`, row-major; only its strict lower part is read.
    lower: Vec<f64>,
    /// The pivots `D`.
    pivots: Vec<f64>,
    /// `G`, row-major `n x m`.
    map: Vec<f64>,
    ports: Vec<PortSpec>,
}

impl MultiportInertia {
    /// The inertia with the row-major `n x n` mass matrix `mass` and the row-major `n x ports`
    /// port map `map`, both in the rows' units (see the module docs). Each port receives an
    /// effort and returns a flow in `kind`.
    ///
    /// # Errors
    /// [`ElementError::Parameter`] when `mass` has not `n x n` entries for some `n >= 1`, `map`
    /// has not `n x ports` entries with `ports >= 1`, or an entry of either is not finite;
    /// [`ElementError::Asymmetric`] when an entry of `mass` differs from its mirror (a signed zero
    /// equals its mirror); [`ElementError::NotPositiveDefinite`] when a pivot of `L D L^T` (the
    /// Cholesky pivot before its square root) is not positive and finite.
    #[expect(
        clippy::float_cmp,
        reason = "symmetry is tested exactly: an entry must equal its mirror, not lie near it"
    )]
    pub fn new(
        kind: PortKind,
        mass: Vec<f64>,
        map: Vec<f64>,
        ports: usize,
    ) -> Result<Self, ElementError> {
        let n = mass.len().isqrt();
        if n == 0 || n * n != mass.len() {
            return Err(ElementError::Parameter {
                name: "mass matrix entries (n x n, n >= 1)",
                value: mass.len() as f64,
            });
        }
        if ports == 0 || n.checked_mul(ports) != Some(map.len()) {
            return Err(ElementError::Parameter {
                name: "port map entries (n x ports, ports >= 1)",
                value: map.len() as f64,
            });
        }
        for (name, values) in [("mass matrix entry", &mass), ("port map entry", &map)] {
            if let Some(&value) = values.iter().find(|value| !value.is_finite()) {
                return Err(ElementError::Parameter { name, value });
            }
        }
        for row in 0..n {
            for column in 0..row {
                if mass[row * n + column] != mass[column * n + row] {
                    return Err(ElementError::Asymmetric { row, column });
                }
            }
        }
        let (lower, pivots) = factor(&mass, n)?;
        let ports = (0..ports)
            .map(|j| port(&format!("p{j}"), Some(kind), PortRole::EffortIn))
            .collect();
        Ok(Self {
            n,
            mass,
            lower,
            pivots,
            map,
            ports,
        })
    }

    /// `x = M^-1 x` by `L`, `D` and `L^T` in turn; with one row, `x / M`.
    fn solve(&self, x: &mut [f64]) {
        let n = self.n;
        for i in 0..n {
            for k in 0..i {
                x[i] -= self.lower[i * n + k] * x[k];
            }
        }
        for (value, pivot) in x.iter_mut().zip(&self.pivots) {
            *value /= pivot;
        }
        for i in (0..n).rev() {
            for k in i + 1..n {
                x[i] -= self.lower[k * n + i] * x[k];
            }
        }
    }

    /// `scale G e` into `out`, each row's sum taken in port order before the scaling.
    fn force(&self, scale: f64, efforts: &[f64], out: &mut Vec<f64>) {
        let rows = self.map.chunks_exact(self.ports.len());
        out.clear();
        out.extend(rows.map(|row| scale * dot(row, efforts)));
    }

    /// Port `j`'s flow `(G^T w)_j`, summed in row order.
    fn flow(&self, j: usize, w: impl Fn(usize) -> f64) -> f64 {
        let m = self.ports.len();
        (0..self.n).fold(-0.0, |sum, i| sum + self.map[i * m + j] * w(i))
    }

    /// `sum_ij (M_ij / 2) v_i v_j`, or with every factor taken by its size, the size of the
    /// balance's terms. With one row both are `Mass`'s `(m / 2) v v`, to the bit.
    fn kinetic(&self, v: &[f64], size: bool) -> f64 {
        let mut total = -0.0;
        for (row, &vi) in self.mass.chunks_exact(self.n).zip(v) {
            for (&mij, &vj) in row.iter().zip(v) {
                total += if size {
                    0.5 * mij.abs() * vi.abs() * vj.abs()
                } else {
                    0.5 * mij * vi * vj
                };
            }
        }
        total
    }

    /// One finite velocity per row, with a finite kinetic energy.
    fn admits(&self, v: &[f64]) -> bool {
        v.len() == self.n && v.iter().all(|x| x.is_finite()) && self.kinetic(v, false).is_finite()
    }

    /// One velocity per row.
    fn fits(&self, v: &[f64]) -> Result<(), ElementError> {
        let (name, value) = ("state velocities (one per row of M)", v.len() as f64);
        (v.len() == self.n)
            .then_some(())
            .ok_or(ElementError::Parameter { name, value })
    }

    /// `out = (dt / 2) G^T M^-1 G d`, the same at every state and input.
    fn slope(&self, dt_s: f64, d: &[f64], out: &mut [f64]) {
        let mut x = Vec::with_capacity(self.n);
        self.force(0.5 * dt_s, d, &mut x);
        self.solve(&mut x);
        for (j, value) in out.iter_mut().enumerate() {
            *value = self.flow(j, |i| x[i]);
        }
    }
}

/// A sum of products from a negative zero, which is the identity of addition: a single product
/// keeps its own bits.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(-0.0, |sum, (x, y)| sum + x * y)
}

/// `M = L D L^T` with `L` unit lower triangular: `L` row-major and the pivots `D`, or the first
/// pivot that is not positive and finite. A symmetric matrix is positive definite exactly when
/// every pivot is positive; pivot `j` is the Cholesky pivot of row `j` before its square root.
fn factor(m: &[f64], n: usize) -> Result<(Vec<f64>, Vec<f64>), ElementError> {
    let (mut lower, mut pivots) = (vec![0.0; n * n], vec![0.0; n]);
    for j in 0..n {
        let mut pivot = m[j * n + j];
        for k in 0..j {
            pivot -= lower[j * n + k] * lower[j * n + k] * pivots[k];
        }
        if !(pivot > 0.0 && pivot.is_finite()) {
            return Err(ElementError::NotPositiveDefinite {
                pivot: j,
                value: pivot,
            });
        }
        pivots[j] = pivot;
        for i in j + 1..n {
            let mut sum = m[i * n + j];
            for k in 0..j {
                sum -= lower[i * n + k] * lower[j * n + k] * pivots[k];
            }
            lower[i * n + j] = sum / pivot;
        }
    }
    Ok((lower, pivots))
}

impl Component for MultiportInertia {
    type State = Vec<f64>;
    type Error = ElementError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::elements::MultiportInertia",
        state_family: "generalized-velocity",
        schema_grammar: "n f64 velocities, one per row of the mass matrix; finite, with finite \
                         kinetic energy",
        codec_grammar: "n u64 words: the bits of each velocity in row order",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Vec<f64>, words: &mut Vec<u64>) -> bool {
        if !self.admits(state) {
            return false;
        }
        words.extend(state.iter().map(|value| value.to_bits()));
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Vec<f64>> {
        let state: Vec<f64> = words.iter().map(|&word| f64::from_bits(word)).collect();
        self.admits(&state).then_some(state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        // The factor is derived from the admitted matrix; the domain is part of the ports.
        words.extend([self.n as u64, self.ports.len() as u64]);
        words.extend(self.mass.iter().map(|value| value.to_bits()));
        words.extend(self.map.iter().map(|value| value.to_bits()));
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, ElementError> {
        self.fits(state)?;
        Ok(self.kinetic(state, false))
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        inputs_finite(inputs)?;
        self.fits(old)?;
        let n = self.n;
        self.force(dt_s, inputs, new);
        self.solve(new);
        for (velocity, start) in new.iter_mut().zip(old) {
            *velocity += start;
        }
        let mean = |i: usize| 0.5 * (old[i] + new[i]);
        for (j, out) in outputs.iter_mut().enumerate() {
            *out = self.flow(j, mean);
        }
        let mut stored_change_j = -0.0;
        for (i, row) in self.mass.chunks_exact(n).enumerate() {
            for (j, mij) in row.iter().enumerate() {
                stored_change_j += mij * mean(i) * (new[j] - old[j]);
            }
        }
        let work = inputs
            .iter()
            .zip(outputs.iter())
            .fold(-0.0, |sum, (e, y)| sum + (dt_s * e * y).abs());
        let scale = self.kinetic(old, true) + self.kinetic(new, true) + work;
        // The solve errs by a few n units of |L| |D| |L^T|, at most n times the scale for a
        // positive definite M; the port sums by m units. One row and one port: the Mass's
        // allowance.
        let units = (n * n * self.ports.len()) as f64;
        let receipt = Receipt::lossless(stored_change_j, roundoff(units * scale));
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for MultiportInertia {
    fn jvp(
        &self,
        _old: &Vec<f64>,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        self.slope(dt_s, d_inputs, d_outputs);
        Ok(())
    }

    fn vjp(
        &self,
        _old: &Vec<f64>,
        dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // `(dt / 2) G^T M^-1 G` is symmetric: its own adjoint.
        self.slope(dt_s, bar_outputs, bar_inputs);
        Ok(())
    }
}
