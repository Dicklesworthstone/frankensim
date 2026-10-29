//! Two-way exterior air on the prepared owner.
//!
//! A passive radiation load (`radiation::Binding`) is a set of unit-mass acoustic oscillators. Each one is driven by a
//! signed row l of the source modes' velocities and reacts on them with F = -l p. Here the poles run as a bank beside
//! the dense mechanical transition, not inside its exponential. A travel-time-resolved skin load has thousands of
//! poles, and folding them into the dense state would multiply its dimension and the bytes every step streams.
//!
//! Over one step each side is integrated exactly while the other side's step average is held:
//! - the mechanics feel the held port force u = -L^T dq_air / dt;
//! - each pole feels the held drive w = L dq_ports / dt.
//!
//! Both are exact step averages, because q' = v on both sides, so the interconnection's work cancels:
//! u . dq_ports = -(dq_air . L dq_ports) / dt = -w . dq_air. This is the implicit midpoint rule applied to the
//! interconnection alone: second order in dt, A-stable and lossless. Each side's own step stays exact.
//!
//! The step is linear. Let B be the transition's own port block and D each pole's held-drive displacement. Then the
//! exchange is one solve of the constant ports x ports operator I + L^T D L B / dt^2, factored once. Contacts keep
//! their held reactions around the air-loaded affine response, which `predict_free` and `response` return.
//!
//! A pole's state is (q, p) with storage 1/2 ((omega q)^2 + p^2), as on the exact block. Its radiated loss over a step
//! is its own Van Loan quadratic of 2 zeta omega p^2.

use super::*;
use fs_la::factor::{Lu, lu};

/// Hard ceiling on the poles of one prepared air bank: the shared passive model's own limit.
pub(super) const MAX_PREPARED_AIR_POLES: usize = crate::render::plate::impact::radiation::MAX_POLES;

/// One pole's exact step under a held drive, in the scaled coordinates z = (omega q, p, w).
struct PoleStep {
    omega: f64,
    /// The end (omega q, p) rows of the exact step from z.
    rows: [[f64; 3]; 2],
    /// The step's radiated loss is z^T loss z.
    loss: [f64; 9],
}

/// The air bank of a prepared network and its constant exchange operator.
pub(super) struct PreparedAir {
    sources: Vec<usize>,
    /// Pole-major coupling L, poles x ports [1/s].
    coupling: Vec<f64>,
    poles: Vec<PoleStep>,
    /// The transition's own q rows of every source mode, ports x 3 count, before the mode weight: row . z / weight is the
    /// end displacement, the very sum the transition's product forms.
    port_rows: Vec<f64>,
    /// End displacement of every mode per unit held force on each port, count x ports.
    port_response: Vec<f64>,
    /// L^T D L, ports x ports.
    exchange: Vec<f64>,
    /// I + L^T D L B / dt^2.
    factor: Lu,
    dt: f64,
}

/// One accepted air step.
pub(super) struct AirStep {
    /// End (q, p) of every pole.
    pub(super) state: Vec<f64>,
    pub(super) loss: f64,
    pub(super) delta_h: f64,
    /// Work the mechanics delivered to the air, w . dq_air.
    pub(super) work: f64,
    /// The held port force the step implies, -L^T dq_air / dt.
    pub(super) force: Vec<f64>,
}

impl PreparedAir {
    /// The bytes this bank keeps beside the network: port rows and responses, the exchange and its factor, the coupling
    /// and every pole's step.
    pub(super) fn dense_entries(count: usize, memory: usize, ports: usize, poles: usize) -> Option<usize> {
        let rows = ports.checked_mul(count.checked_mul(3)?.checked_add(memory)?)?;
        let response = count.checked_mul(ports)?;
        let square = ports.checked_mul(ports)?.checked_mul(4)?;
        let coupling = poles.checked_mul(ports)?;
        let steps = poles.checked_mul(15)?;
        rows.checked_add(response)?.checked_add(square)?.checked_add(coupling)?.checked_add(steps)
    }

    /// The multiply-adds of the bank's preparation: every pole's 3 x 3 exponentials and loss integral, L^T D L, the
    /// exchange product and its factorization.
    pub(super) fn setup_terms(ports: usize, poles: usize) -> Option<u64> {
        let (m, p) = (ports as u64, poles as u64);
        let pole = 27_u64.checked_mul(80 * 8 + 2 * MAX_SCALES as u64 + 2)?;
        p.checked_mul(pole)?.checked_add(p.checked_mul(m.checked_mul(m)?)?)?.checked_add(m.checked_pow(3)?.checked_mul(2)?)
    }

    pub(super) fn new(
        binding: &radiation::Binding, count: usize, memory: usize, dt: f64, weights: &[f64], transition: &[f64], gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        let model = &binding.model;
        let (ports, d) = (model.ports, 3 * count + memory);
        if model.poles.is_empty() || model.poles.len() > MAX_PREPARED_AIR_POLES || binding.source_modes.len() != ports
            || binding.source_modes.iter().any(|&mode| mode >= count) || weights.len() != count
            || transition.len() != d * d {
            return Err(invalid("prepared air needs a complete pole bank on original modes"));
        }
        let budget = DenseExponentialBudget { max_dimension: 6, max_dense_bytes: 1 << 16, max_multiply_terms: 1 << 24 };
        let mut poles = Vec::with_capacity(model.poles.len());
        let mut coupling = Vec::with_capacity(model.poles.len() * ports);
        for pole in &model.poles {
            poll(Some(gate))?;
            coupling.extend_from_slice(&pole.coupling);
            let (omega, rate) = (pole.omega, finite(2.0 * pole.zeta * pole.omega)?);
            let generator = [0.0, omega, 0.0, -omega, -rate, 1.0, 0.0, 0.0, 0.0];
            let (levels, h) = scaled_levels(&generator, 3, dt, budget, gate)?;
            let resistance = [0.0, 0.0, 0.0, 0.0, rate, 0.0, 0.0, 0.0, 0.0];
            let loss = gramian(&generator, &resistance, &levels, h, 3, budget, gate)?;
            let step = levels.last().ok_or_else(|| invalid("prepared air step"))?;
            let pole = PoleStep {
                omega,
                rows: [[step[0], step[1], step[2]], [step[3], step[4], step[5]]],
                loss: loss.try_into().map_err(|_| invalid("prepared air loss"))?,
            };
            pole.validate_energy_identity()?;
            poles.push(pole);
        }
        let sources = binding.source_modes.clone();
        let mut port_rows = Vec::with_capacity(ports * d);
        for &mode in &sources {
            port_rows.extend_from_slice(&transition[2 * mode * d..(2 * mode + 1) * d]);
        }
        let mut port_response = vec![0.0; count * ports];
        for i in 0..count {
            for (k, &mode) in sources.iter().enumerate() {
                port_response[i * ports + k] = finite(transition[2 * i * d + 2 * count + memory + mode] / weights[i])?;
            }
        }
        // D: each pole's displacement per unit held drive over the step.
        let drive: Vec<f64> = poles.iter().map(|pole| finite(pole.rows[0][2] / pole.omega)).collect::<Result<_, _>>()?;
        let mut exchange = vec![0.0; ports * ports];
        for (row, &d_a) in coupling.chunks_exact(ports).zip(&drive) {
            poll(Some(gate))?;
            for i in 0..ports {
                for j in 0..ports {
                    exchange[i * ports + j] = finite(exchange[i * ports + j] + row[i] * d_a * row[j])?;
                }
            }
        }
        let dt2 = finite(dt * dt)?;
        let mut operator = vec![0.0; ports * ports];
        for i in 0..ports {
            for j in 0..ports {
                let gb = (0..ports).try_fold(0.0, |s, k| finite(s + exchange[i * ports + k] * port_response[sources[k] * ports + j]))?;
                operator[i * ports + j] = finite(if i == j { 1.0 } else { 0.0 } + gb / dt2)?;
            }
        }
        let factor = lu(&operator, ports).map_err(|_| invalid("prepared air exchange operator is singular"))?;
        poll(Some(gate))?;
        Ok(Self { sources, coupling, poles, port_rows, port_response, exchange, factor, dt })
    }

    pub(super) fn sources(&self) -> &[usize] { &self.sources }

    pub(super) fn state_len(&self) -> usize { 2 * self.poles.len() }

    /// Stored acoustic energy of a complete (q, p) state.
    pub(super) fn energy(&self, air: &[f64]) -> Result<f64, ModalCouplingError> {
        if air.len() != self.state_len() {
            return Err(invalid("prepared air state layout"));
        }
        self.poles.iter().zip(air.chunks_exact(2))
            .try_fold(0.0, |sum, (pole, s)| finite(sum + 0.5 * ((pole.omega * s[0]).powi(2) + s[1] * s[1])))
    }

    /// Each source mode's end-displacement increment for the mechanical state vector z. Every row's sum is formed in the
    /// transition product's own order, so these are the bits of the endpoint the same z propagates to.
    pub(super) fn port_increments(&self, z: &[f64], x: &[f64], weights: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        let d = z.len();
        self.sources.iter().zip(self.port_rows.chunks_exact(d)).map(|(&mode, row)| {
            let end = finite(row.iter().zip(z).fold(0.0, |s, (a, b)| s + a * b))?;
            finite(finite(end / weights[mode])? - x[2 * mode])
        }).collect()
    }

    /// Each pole's displacement increment without a drive.
    pub(super) fn free_increments(&self, air: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        if air.len() != self.state_len() || air.iter().any(|v| !v.is_finite()) {
            return Err(invalid("prepared air state layout"));
        }
        self.poles.iter().zip(air.chunks_exact(2)).map(|(pole, s)| {
            let end = finite(pole.rows[0][0] * (pole.omega * s[0]) + pole.rows[0][1] * s[1])?;
            finite(finite(end / pole.omega)? - s[0])
        }).collect()
    }

    /// The held port force of a step whose ports move dq_ports and whose poles move dq_air before the exchange:
    /// (I + G B / dt^2) u = -L^T dq_air / dt - G dq_ports / dt^2.
    pub(super) fn exchange_force(&self, dq_ports: &[f64], dq_air: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        let ports = self.sources.len();
        if dq_ports.len() != ports || dq_air.len() != self.poles.len() {
            return Err(invalid("prepared air exchange layout"));
        }
        let mut rhs = vec![0.0; ports];
        for (row, &dq) in self.coupling.chunks_exact(ports).zip(dq_air) {
            for (value, l) in rhs.iter_mut().zip(row) {
                *value -= l * dq;
            }
        }
        let dt2 = self.dt * self.dt;
        for (value, g) in rhs.iter_mut().zip(self.exchange.chunks_exact(ports)) {
            let mixed = g.iter().zip(dq_ports).fold(0.0, |s, (a, b)| s + a * b);
            *value = finite(*value / self.dt - mixed / dt2)?;
        }
        self.factor.solve(&mut rhs);
        rhs.into_iter().map(finite).collect()
    }

    /// Add the air-loaded part of a response: `response` holds every mode's end displacement for a held force; the ports'
    /// own displacements then drive the exchange with the air at rest.
    pub(super) fn load_response(&self, response: &mut [f64]) -> Result<(), ModalCouplingError> {
        let ports = self.sources.len();
        let dq_ports: Vec<_> = self.sources.iter().map(|&mode| response[mode]).collect();
        let force = self.exchange_force(&dq_ports, &vec![0.0; self.poles.len()])?;
        for (value, row) in response.iter_mut().zip(self.port_response.chunks_exact(ports)) {
            *value = finite(*value + row.iter().zip(&force).fold(0.0, |s, (a, b)| s + a * b))?;
        }
        Ok(())
    }

    /// Each port's first-order rounding bound on the difference of its two held forces, `exchange_force`'s and
    /// `advance`'s. Both are formed from increments that are differences of end and start coordinates: a port's row sum
    /// over z less its start displacement, a pole's end less its start. The rounding of a difference is relative to its
    /// operands, not to the difference, so a mode resting on a preload resolves its first increments only to the last
    /// bits of its displacement, far above the force they imply: a violin body mode resting on the strings' preload at
    /// x = -1.117e-7 moves 5.37e-20 in its first step at 48 kHz, and the two forces of that step differ by 1.8e-17 N,
    /// ulp(x) through G / dt^2, on a force of 1.5e-13 N. The operands reach a port's force through
    /// |l| |D| |l| / dt^2 (ports) and |l| / dt (poles). gamma_n = n u / (1 - n u) is taken over the longest chain of
    /// roundings an increment passes: the held input's sum and the displacement's weight (2), the row sum over z (d),
    /// its division and the start's subtraction (2), the pole drive (ports + 1), the pole's row, division and
    /// subtraction (5) and the force's sum over the poles (poles + 1). Each of the two evaluations carries it once.
    pub(super) fn force_rounding(&self, z: &[f64], x: &[f64], weights: &[f64], air: &[f64], dq_ports: &[f64])
        -> Result<Vec<f64>, ModalCouplingError>
    {
        let ports = self.sources.len();
        let d = z.len();
        if air.len() != self.state_len() || air.iter().any(|v| !v.is_finite()) || dq_ports.len() != ports {
            return Err(invalid("prepared air state layout"));
        }
        let displacement: Vec<f64> = self.sources.iter().zip(self.port_rows.chunks_exact(d)).map(|(&mode, row)| {
            let terms = row.iter().zip(z).fold(0.0, |s, (a, b)| s + (a * b).abs());
            finite(terms / weights[mode] + x[2 * mode].abs())
        }).collect::<Result<_, _>>()?;
        let drives = self.drives(dq_ports)?;
        let dt2 = self.dt * self.dt;
        let mut operands = vec![0.0; ports];
        for ((pole, s), (&w, row)) in self.poles.iter().zip(air.chunks_exact(2)).zip(drives.iter().zip(self.coupling.chunks_exact(ports))) {
            let through = row.iter().zip(&displacement).fold(0.0, |acc, (l, m)| acc + (l * m).abs());
            let terms = pole.rows[0].iter().zip([pole.omega * s[0], s[1], w]).fold(0.0, |acc, (a, b)| acc + (a * b).abs());
            let pole_side = (terms / pole.omega + s[0].abs()) / self.dt;
            let port_side = (pole.rows[0][2] / pole.omega).abs() * through / dt2;
            for (value, l) in operands.iter_mut().zip(row) {
                *value = finite(*value + l.abs() * (port_side + pole_side))?;
            }
        }
        let unit = 0.5 * f64::EPSILON;
        let n = (d + ports + self.poles.len() + 11) as f64;
        let gamma = n * unit / (1.0 - n * unit);
        operands.into_iter().map(|value| finite(2.0 * gamma * value)).collect()
    }

    /// The pole drives L dq_ports / dt of accepted port increments.
    fn drives(&self, dq_ports: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
        if dq_ports.len() != self.sources.len() {
            return Err(invalid("prepared air drive layout"));
        }
        self.coupling.chunks_exact(self.sources.len())
            .map(|row| finite(row.iter().zip(dq_ports).fold(0.0, |s, (a, b)| s + a * b) / self.dt)).collect()
    }

    /// The air's exact step under the accepted port increments.
    pub(super) fn advance(&self, air: &[f64], dq_ports: &[f64]) -> Result<AirStep, ModalCouplingError> {
        if air.len() != self.state_len() || air.iter().any(|v| !v.is_finite()) {
            return Err(invalid("prepared air state layout"));
        }
        let ports = self.sources.len();
        let drives = self.drives(dq_ports)?;
        let mut state = Vec::with_capacity(air.len());
        let (mut loss, mut delta_h, mut work) = (0.0, 0.0, 0.0);
        let mut force = vec![0.0; ports];
        for ((pole, s), (&w, row)) in self.poles.iter().zip(air.chunks_exact(2)).zip(drives.iter().zip(self.coupling.chunks_exact(ports))) {
            let z = [pole.omega * s[0], s[1], w];
            let end = [dot3(&pole.rows[0], &z)?, dot3(&pole.rows[1], &z)?];
            let q = finite(end[0] / pole.omega)?;
            let dq = finite(q - s[0])?;
            state.extend([q, end[1]]);
            loss = finite(loss + quadratic(&pole.loss, &z)?)?;
            delta_h = finite(delta_h + 0.5 * (end[0] * end[0] + end[1] * end[1]) - 0.5 * (z[0] * z[0] + z[1] * z[1]))?;
            work = finite(work + w * dq)?;
            for (f, l) in force.iter_mut().zip(row) {
                *f = finite(*f - l * dq / self.dt)?;
            }
        }
        let scale = delta_h.abs().max(loss).max(work.abs()).max(self.energy(air)?).max(f64::MIN_POSITIVE);
        if (delta_h + loss - work).abs() > MATRIX_TOL * scale {
            return Err(invalid("prepared air step energy residual"));
        }
        Ok(AirStep { state, loss, delta_h, work, force })
    }

    /// The radiated loss of the step from `air` under the accepted port increments, the same value `advance` returns.
    pub(super) fn loss(&self, air: &[f64], dq_ports: &[f64]) -> Result<f64, ModalCouplingError> {
        if air.len() != self.state_len() {
            return Err(invalid("prepared air state layout"));
        }
        let drives = self.drives(dq_ports)?;
        self.poles.iter().zip(air.chunks_exact(2)).zip(&drives)
            .try_fold(0.0, |sum, ((pole, s), &w)| finite(sum + quadratic(&pole.loss, &[pole.omega * s[0], s[1], w])?))
    }
}

impl PoleStep {
    /// 1/2 (T^T S T - S) + W - work = 0 as a quadratic form in z = (omega q, p, w), with S = diag(1, 1, 0) and the work
    /// w (q(dt) - q(0)): the step's storage change and loss are exactly the drive's work.
    fn validate_energy_identity(&self) -> Result<(), ModalCouplingError> {
        let t = [self.rows[0], self.rows[1], [0.0, 0.0, 1.0]];
        let mut defect = [0.0; 9];
        let mut scale = [0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                let propagated = t[0][i] * t[0][j] + t[1][i] * t[1][j];
                let storage = if i == j && i < 2 { 1.0 } else { 0.0 };
                defect[i * 3 + j] = finite(0.5 * (propagated - storage) + self.loss[i * 3 + j])?;
                scale[i * 3 + j] = finite(0.5 * (propagated.abs() + storage) + self.loss[i * 3 + j].abs())?;
            }
        }
        for j in 0..3 {
            let initial = if j == 0 { 1.0 } else { 0.0 };
            let work = finite(0.5 * (self.rows[0][j] - initial) / self.omega)?;
            defect[2 * 3 + j] -= work;
            defect[j * 3 + 2] -= work;
            scale[2 * 3 + j] += work.abs();
            scale[j * 3 + 2] += work.abs();
        }
        if matrix_norm(&defect, 3)? > MATRIX_TOL * matrix_norm(&scale, 3)?.max(f64::MIN_POSITIVE) {
            return Err(invalid("prepared air pole energy identity"));
        }
        Ok(())
    }
}

fn dot3(row: &[f64; 3], z: &[f64; 3]) -> Result<f64, ModalCouplingError> {
    finite(row[0] * z[0] + row[1] * z[1] + row[2] * z[2])
}
