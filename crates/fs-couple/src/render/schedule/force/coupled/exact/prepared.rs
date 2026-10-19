//! Cold exact mechanical propagation with full per-resistance loss integrals.
//!
//! z = [D q, v, u] contains every original mass-normalized coordinate and
//! every held generalized force. D only balances units (sqrt(K_ii), or one
//! on free coordinates); it is invertible and removes no state. The constant
//! augmented generator is the same q'=v, v'=-Kq-Cv+u pencil as ExactBlock.
//! Short-interval Van Loan integrals avoid exp(-A^T dt) overflow for stiff
//! viscous memory: W(2h)=W(h)+F(h)^T W(h) F(h). Each physical dashpot keeps
//! its OWN integral. This is a floating-point matrix exponential, not a
//! nonlinear/continuum accuracy certificate or a real-time timing promise.

use super::*;
use fs_phs::linear_exact::{DenseExponentialBudget, exponential_with_budget};

/// Explicit upper bound on complete mechanical modes in the cold dense path.
/// Each mode retains displacement, velocity and its independent held input.
pub const MAX_PREPARED_MODAL_MODES: usize = 512;

/// Additional cold dense resources; original modal/connection/work and
/// physical/energy admissions remain binding. There is deliberately no default.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PreparedModalBudget {
    /// Maximum TOTAL original mechanical modes, up to 512; none are dropped.
    pub max_modes: usize,
    /// Conservative peak bytes of additional dense preparation/cache storage.
    /// Includes up to 48 scale levels and all separately retained loss matrices.
    pub max_dense_bytes: usize,
    /// Conservative scalar multiply-add term ceiling, including every shared
    /// exponential's 32 Taylor products and at most 48 squarings.
    pub max_multiply_terms: u64,
}

pub(super) struct PreparedNetwork {
    count: usize,
    dt: f64,
    weights: Vec<f64>,
    stiffness: Vec<f64>,
    rest_force: Vec<f64>,
    transition: Vec<f64>,
    component_loss: Option<Vec<f64>>,
    connection_loss: Vec<Option<Vec<f64>>>,
}

const MATRIX_TOL: f64 = 1e-11;
const MAX_SCALES: usize = 48;

impl PreparedModalBudget {
    /// Validate fixed capacity fields without allocating or choosing a backend.
    pub fn validate(self) -> Result<(), ModalCouplingError> {
        if self.max_modes == 0 || self.max_modes > MAX_PREPARED_MODAL_MODES
            || self.max_dense_bytes == 0 || self.max_multiply_terms == 0 {
            return Err(invalid("prepared exact modal byte/work/mode budget"));
        }
        Ok(())
    }

    fn admit(self, count: usize, losses: usize, links: usize) -> Result<DenseExponentialBudget, ModalCouplingError> {
        self.validate()?;
        if count == 0 || count > self.max_modes {
            return Err(invalid("complete prepared modal count exceeds max_modes; no modes were dropped"));
        }
        // The hard mode limit makes these small, but every resource product
        // remains checked before the first prepared dense allocation.
        let d = count.checked_mul(3).ok_or_else(|| invalid("prepared state extent"))?;
        let extent = d.checked_mul(d).ok_or_else(|| invalid("prepared matrix extent"))?;
        let matrices = losses.checked_add(128).ok_or_else(|| invalid("prepared workspace count"))?;
        let bytes = extent.checked_mul(matrices)
            .and_then(|n| n.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|n| n.checked_add(links.checked_mul(core::mem::size_of::<Option<Vec<f64>>>())?))
            .ok_or_else(|| invalid("prepared dense byte overflow"))?;
        let cube = (d as u64).checked_pow(3).ok_or_else(|| invalid("prepared dense work overflow"))?;
        // A 2d Van Loan exponential per nonzero channel; two d transition
        // exponentials for the independent h versus h/2 check; Gramian
        // doubling, energy-identity validation and all physical assembly.
        let products = (losses as u64).checked_mul(80 * 8 + 2 * MAX_SCALES as u64 + 2)
            .and_then(|n| n.checked_add(2 * 80 + 2 * MAX_SCALES as u64 + 16))
            .ok_or_else(|| invalid("prepared product work overflow"))?;
        let terms = cube.checked_mul(products)
            .and_then(|n| n.checked_add((links as u64).checked_mul((count as u64).checked_pow(2)?)?.checked_mul(8)?))
            .ok_or_else(|| invalid("prepared complete work overflow"))?;
        if bytes > self.max_dense_bytes {
            return Err(invalid("prepared exact dense bytes exceed max_dense_bytes"));
        }
        if terms > self.max_multiply_terms {
            return Err(invalid("prepared exact multiply terms exceed max_multiply_terms"));
        }
        Ok(DenseExponentialBudget {
            max_dimension: 2 * d,
            max_dense_bytes: self.max_dense_bytes,
            max_multiply_terms: self.max_multiply_terms,
        })
    }
}

impl PreparedNetwork {
    pub(super) fn new(
        models: &[ModalAcousticTimeModel],
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        count: usize,
        dt: f64,
        budget: PreparedModalBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if !dt.is_finite() || dt <= 0.0 || columns.len() != links.len()
            || columns.iter().any(|c| c.len() != count)
            || models.iter().map(|m| m.modes().len()).sum::<usize>() != count {
            return Err(invalid("prepared complete mechanical layout/clock"));
        }
        let has_component_loss = models.iter().flat_map(|m| m.modes())
            .any(|m| m.damping_ratio != 0.0 && m.angular_frequency_rad_s != 0.0);
        let losses = usize::from(has_component_loss)
            + links.iter().filter(|c| c.damping_n_s_m != 0.0).count();
        let exp_budget = budget.admit(count, losses, links.len())?;
        let d = 3 * count;
        let mut stiffness = vec![0.0; count * count];
        let mut damping = vec![0.0; count * count];
        let mut component_rates = vec![0.0; count];
        let mut rest_force = vec![0.0; count];
        for (i, mode) in models.iter().flat_map(|m| m.modes()).enumerate() {
            stiffness[i * count + i] = finite(mode.angular_frequency_rad_s.powi(2))?;
            component_rates[i] = finite(2.0 * mode.damping_ratio * mode.angular_frequency_rad_s)?;
            damping[i * count + i] = component_rates[i];
        }
        for (column, link) in columns.iter().zip(links) {
            poll(Some(gate))?;
            for i in 0..count {
                rest_force[i] = finite(rest_force[i] + link.stiffness_n_m * link.rest_extension_m * column[i])?;
                for j in 0..count {
                    stiffness[i * count + j] = finite(stiffness[i * count + j] + link.stiffness_n_m * column[i] * column[j])?;
                    damping[i * count + j] = finite(damping[i * count + j] + link.damping_n_s_m * column[i] * column[j])?;
                }
            }
        }
        let weights: Vec<_> = (0..count).map(|i| {
            let k = stiffness[i * count + i];
            if k > 0.0 { det::sqrt(k) } else { 1.0 }
        }).collect();
        let mut generator = vec![0.0; d * d];
        for i in 0..count {
            generator[(2 * i) * d + 2 * i + 1] = weights[i];
            generator[(2 * i + 1) * d + 2 * count + i] = 1.0;
            for j in 0..count {
                generator[(2 * i + 1) * d + 2 * j] = finite(-stiffness[i * count + j] / weights[j])?;
                generator[(2 * i + 1) * d + 2 * j + 1] = -damping[i * count + j];
            }
        }
        let norm = matrix_norm(&generator, d)?;
        let mut h = dt;
        let mut scales = 0;
        while norm * h > 0.5 {
            if scales == MAX_SCALES {
                return Err(invalid("prepared exact scaling exceeds 48 doublings"));
            }
            h *= 0.5;
            scales += 1;
        }
        let scaled: Vec<_> = generator.iter().map(|v| v * h).collect();
        let mut levels = Vec::with_capacity(scales + 1);
        levels.push(exponential_with_budget(&scaled, d, exp_budget).map_err(ModalCouplingError::Propagation)?);
        for i in 0..scales {
            levels.push(product(&levels[i], &levels[i], d, gate)?);
        }
        // Compare the complete transition against one finer independent series
        // evaluation. Extra squaring roundoff is admitted, not silently trusted.
        let finer: Vec<_> = scaled.iter().map(|v| v * 0.5).collect();
        let mut check = exponential_with_budget(&finer, d, exp_budget).map_err(ModalCouplingError::Propagation)?;
        for _ in 0..=scales { check = product(&check, &check, d, gate)?; }
        matrix_agrees(levels.last().unwrap(), &check, d, "prepared transition refinement residual")?;

        let mut component_loss = None;
        if has_component_loss {
            let mut resistance = vec![0.0; d * d];
            for (i, c) in component_rates.iter().enumerate() {
                resistance[(2 * i + 1) * d + 2 * i + 1] = *c;
            }
            component_loss = Some(gramian(&generator, &resistance, &levels, h, d, exp_budget, gate)?);
        }
        let mut connection_loss = Vec::with_capacity(links.len());
        for (column, link) in columns.iter().zip(links) {
            poll(Some(gate))?;
            if link.damping_n_s_m == 0.0 {
                connection_loss.push(None);
                continue;
            }
            let mut resistance = vec![0.0; d * d];
            for i in 0..count {
                for j in 0..count {
                    resistance[(2 * i + 1) * d + 2 * j + 1] =
                        finite(link.damping_n_s_m * column[i] * column[j])?;
                }
            }
            connection_loss.push(Some(gramian(&generator, &resistance, &levels, h, d, exp_budget, gate)?));
        }
        let transition = levels.pop().unwrap();
        let result = Self { count, dt, weights, stiffness, rest_force, transition, component_loss, connection_loss };
        result.validate_energy_identity(gate)?;
        poll(Some(gate))?;
        Ok(result)
    }

    fn state_input(&self, x: &[f64], u: &[f64], dt: f64) -> Result<Vec<f64>, ModalCouplingError> {
        if dt.to_bits() != self.dt.to_bits() || x.len() != 2 * self.count || u.len() != self.count
            || x.iter().chain(u).any(|v| !v.is_finite()) {
            return Err(invalid("prepared exact clock or complete finite state/input changed"));
        }
        let mut z = vec![0.0; 3 * self.count];
        for i in 0..self.count {
            z[2 * i] = finite(self.weights[i] * x[2 * i])?;
            z[2 * i + 1] = x[2 * i + 1];
            z[2 * self.count + i] = u[i];
        }
        Ok(z)
    }

    pub(super) fn hamiltonian(&self, x: &[f64]) -> f64 {
        let n = self.count;
        let mut energy = 0.0;
        for i in 0..n {
            energy += 0.5 * x[2 * i + 1].powi(2);
            for j in 0..n {
                energy += 0.5 * x[2 * i] * self.stiffness[i * n + j] * x[2 * j];
            }
        }
        energy
    }

    pub(super) fn advance(&self, models: &[ModalAcousticTimeModel], external: &[f64], dt: f64)
        -> Result<(StepRecord, Vec<f64>, Vec<f64>), ModalCouplingError>
    {
        if external.len() != self.count {
            return Err(invalid("prepared external force layout"));
        }
        let x = flatten(models);
        let u: Vec<_> = external.iter().zip(&self.rest_force)
            .map(|(f, r)| finite(f + r)).collect::<Result<_, _>>()?;
        let z = self.state_input(&x, &u, dt)?;
        let end = apply(&self.transition, &z)?;
        let mut endpoint = vec![0.0; 2 * self.count];
        let mut y = vec![0.0; self.count];
        let mut supplied = 0.0;
        for i in 0..self.count {
            endpoint[2 * i] = finite(end[2 * i] / self.weights[i])?;
            endpoint[2 * i + 1] = end[2 * i + 1];
            let dq = finite(endpoint[2 * i] - x[2 * i])?;
            y[i] = finite(dq / dt)?;
            supplied = finite(supplied + u[i] * dq)?;
        }
        let (component, connection) = self.loss_sums(&z)?;
        let dissipated = finite(component + connection)?;
        let before = finite(self.hamiltonian(&x))?;
        let after = finite(self.hamiltonian(&endpoint))?;
        let delta_h = finite(after - before)?;
        let scale = before.abs().max(after.abs()).max(supplied.abs()).max(dissipated).max(f64::MIN_POSITIVE);
        let residual = finite(delta_h + dissipated - supplied)?;
        if residual.abs() > MATRIX_TOL * scale {
            return Err(invalid("prepared exact trajectory energy residual"));
        }
        Ok((StepRecord {
            x: endpoint, y, delta_h, dissipated, supplied,
            newton_iters: 0, solver_residual: residual.abs(),
            solver_certificate: None,
        }, x, u))
    }

    pub(super) fn response(&self, column: &[f64], dt: f64) -> Result<Vec<f64>, ModalCouplingError> {
        let z = self.state_input(&vec![0.0; 2 * self.count], column, dt)?;
        let end = apply(&self.transition, &z)?;
        (0..self.count).map(|i| finite(end[2 * i] / self.weights[i])).collect()
    }

    fn loss_sums(&self, z: &[f64]) -> Result<(f64, f64), ModalCouplingError> {
        let component = match &self.component_loss { Some(w) => quadratic(w, z)?, None => 0.0 };
        let mut connection = 0.0;
        for matrix in self.connection_loss.iter().flatten() {
            connection = finite(connection + quadratic(matrix, z)?)?;
        }
        Ok((component, connection))
    }

    pub(super) fn losses(&self, x: &[f64], u: &[f64], dt: f64, total: f64)
        -> Result<(f64, f64), ModalCouplingError>
    {
        let result = self.loss_sums(&self.state_input(x, u, dt)?)?;
        let sum = finite(result.0 + result.1)?;
        if !total.is_finite() || (sum - total).abs() > MATRIX_TOL * sum.abs().max(total.abs()).max(f64::MIN_POSITIVE) {
            return Err(invalid("prepared independent loss channels disagree with total"));
        }
        Ok(result)
    }

    pub(super) fn connection_work(
        &self, x: &[f64], u: &[f64], endpoint: &[f64], dt: f64,
        columns: &[Vec<f64>], links: &[ConnectionLaw], connection_loss: f64,
        gate: Option<&CancelGate>,
    ) -> Result<Vec<crate::port_adoption::PortTransfer>, ModalCouplingError> {
        if links.len() != self.connection_loss.len() || columns.len() != links.len()
            || endpoint.len() != x.len() {
            return Err(invalid("prepared connection observation layout"));
        }
        let z = self.state_input(x, u, dt)?;
        let mut transfers = Vec::with_capacity(links.len());
        let mut total = 0.0;
        for ((column, link), matrix) in columns.iter().zip(links).zip(&self.connection_loss) {
            poll(gate)?;
            let coordinate = |state: &[f64]| -> Result<f64, ModalCouplingError> {
                column.iter().enumerate().try_fold(-link.rest_extension_m,
                    |s, (i, b)| finite(s + b * state[2 * i]))
            };
            let before = coordinate(x)?;
            let after = coordinate(endpoint)?;
            let loss = match matrix { Some(w) => quadratic(w, &z)?, None => 0.0 };
            total = finite(total + loss)?;
            let work = finite(-0.5 * link.stiffness_n_m * (after - before) * (after + before) - loss)?;
            transfers.push(crate::port_adoption::PortTransfer::work(work, false));
        }
        if (total - connection_loss).abs() > MATRIX_TOL * total.abs().max(connection_loss.abs()).max(f64::MIN_POSITIVE) {
            return Err(invalid("prepared per-port losses do not sum to connection loss"));
        }
        Ok(transfers)
    }

    fn validate_energy_identity(&self, gate: &CancelGate) -> Result<(), ModalCouplingError> {
        let (n, d) = (self.count, 3 * self.count);
        let mut storage = vec![0.0; d * d];
        for i in 0..n {
            storage[(2 * i + 1) * d + 2 * i + 1] = 1.0;
            for j in 0..n {
                storage[(2 * i) * d + 2 * j] = finite(self.stiffness[i * n + j] / self.weights[i] / self.weights[j])?;
            }
        }
        let propagated = sandwich(&self.transition, &storage, d, gate)?;
        let mut defect = vec![0.0; d * d];
        let mut scale = vec![0.0; d * d];
        for i in 0..d * d {
            defect[i] = finite(0.5 * (propagated[i] - storage[i]))?;
            scale[i] = finite(0.5 * (propagated[i].abs() + storage[i].abs()))?;
        }
        for matrix in self.component_loss.iter().chain(self.connection_loss.iter().flatten()) {
            for i in 0..d * d {
                defect[i] = finite(defect[i] + matrix[i])?;
                scale[i] = finite(scale[i] + matrix[i].abs())?;
            }
        }
        // u^T (q(dt)-q(0)) is an EXACT held-force work quadratic. Symmetric
        // input/state cross terms count once in z^T work z.
        for i in 0..n {
            let qi = 2 * i;
            let ui = 2 * n + i;
            for j in 0..d {
                let initial = if qi == j { 1.0 } else { 0.0 };
                let work = finite(0.5 * (self.transition[qi * d + j] - initial) / self.weights[i])?;
                defect[ui * d + j] = finite(defect[ui * d + j] - work)?;
                defect[j * d + ui] = finite(defect[j * d + ui] - work)?;
                scale[ui * d + j] = finite(scale[ui * d + j] + work.abs())?;
                scale[j * d + ui] = finite(scale[j * d + ui] + work.abs())?;
            }
        }
        if matrix_norm(&defect, d)? > MATRIX_TOL * matrix_norm(&scale, d)?.max(f64::MIN_POSITIVE) {
            return Err(invalid("prepared complete transition/loss matrix energy identity"));
        }
        Ok(())
    }
}

fn matrix_norm(a: &[f64], n: usize) -> Result<f64, ModalCouplingError> {
    let mut norm = 0.0_f64;
    for row in a.chunks_exact(n) {
        norm = norm.max(row.iter().try_fold(0.0, |s, v| finite(s + v.abs()))?);
    }
    Ok(norm)
}

fn matrix_agrees(a: &[f64], b: &[f64], n: usize, what: &'static str) -> Result<(), ModalCouplingError> {
    let mut residual = 0.0_f64;
    for (left, right) in a.chunks_exact(n).zip(b.chunks_exact(n)) {
        residual = residual.max(left.iter().zip(right).try_fold(0.0,
            |s, (x, y)| finite(s + (x - y).abs()))?);
    }
    let scale = matrix_norm(a, n)?.max(matrix_norm(b, n)?).max(f64::MIN_POSITIVE);
    if residual > MATRIX_TOL * scale { return Err(invalid(what)); }
    Ok(())
}

fn product(a: &[f64], b: &[f64], n: usize, gate: &CancelGate) -> Result<Vec<f64>, ModalCouplingError> {
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        poll(Some(gate))?;
        for k in 0..n {
            for j in 0..n {
                result[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
    if result.iter().any(|v| !v.is_finite()) { return Err(invalid("prepared dense product overflow")); }
    Ok(result)
}

fn sandwich(f: &[f64], w: &[f64], n: usize, gate: &CancelGate) -> Result<Vec<f64>, ModalCouplingError> {
    let wf = product(w, f, n, gate)?;
    let transpose: Vec<_> = (0..n * n).map(|i| f[(i % n) * n + i / n]).collect();
    product(&transpose, &wf, n, gate)
}

fn gramian(
    generator: &[f64], resistance: &[f64], levels: &[Vec<f64>], h: f64,
    d: usize, budget: DenseExponentialBudget, gate: &CancelGate,
) -> Result<Vec<f64>, ModalCouplingError> {
    poll(Some(gate))?;
    let block = 2 * d;
    let mut van_loan = vec![0.0; block * block];
    for i in 0..d {
        for j in 0..d {
            van_loan[i * block + j] = finite(-generator[j * d + i] * h)?;
            van_loan[(d + i) * block + d + j] = finite(generator[i * d + j] * h)?;
            van_loan[i * block + d + j] = finite(resistance[i * d + j] * h)?;
        }
    }
    let exponential = exponential_with_budget(&van_loan, block, budget).map_err(ModalCouplingError::Propagation)?;
    let e: Vec<_> = (0..d * d).map(|i| exponential[(d + i / d) * block + d + i % d]).collect();
    matrix_agrees(&e, &levels[0], d, "prepared loss/transition short-interval disagreement")?;
    let cross: Vec<_> = (0..d * d).map(|i| exponential[(i / d) * block + d + i % d]).collect();
    let transpose: Vec<_> = (0..d * d).map(|i| e[(i % d) * d + i / d]).collect();
    let mut w = product(&transpose, &cross, d, gate)?;
    for level in levels.iter().take(levels.len() - 1) {
        let extra = sandwich(level, &w, d, gate)?;
        for (value, increment) in w.iter_mut().zip(extra) { *value = finite(*value + increment)?; }
    }
    let transpose: Vec<_> = (0..d * d).map(|i| w[(i % d) * d + i / d]).collect();
    matrix_agrees(&w, &transpose, d, "prepared loss Gramian symmetry residual")?;
    // Only antisymmetric floating roundoff is removed. No spectral projection,
    // eigenvalue clamp or omitted channel is used to make a loss look positive.
    for i in 0..d {
        for j in i + 1..d {
            let value = finite(0.5 * w[i * d + j] + 0.5 * w[j * d + i])?;
            w[i * d + j] = value;
            w[j * d + i] = value;
        }
    }
    Ok(w)
}

fn apply(matrix: &[f64], z: &[f64]) -> Result<Vec<f64>, ModalCouplingError> {
    matrix.chunks_exact(z.len()).map(|row| row.iter().zip(z)
        .try_fold(0.0, |s, (a, b)| finite(s + a * b))).collect()
}

fn quadratic(matrix: &[f64], z: &[f64]) -> Result<f64, ModalCouplingError> {
    let mut sum = 0.0;
    let mut magnitude = 0.0;
    for (i, row) in matrix.chunks_exact(z.len()).enumerate() {
        for (j, value) in row.iter().enumerate() {
            let term = finite(z[i] * value * z[j])?;
            sum = finite(sum + term)?;
            magnitude = finite(magnitude + term.abs())?;
        }
    }
    // A Gramian is PSD in exact arithmetic. Only a negative value below its
    // explicit summation-roundoff allowance is reported as zero physical loss.
    // A larger negative quadratic refuses the trial and preserves accepted state.
    let roundoff = 4096.0 * f64::EPSILON * magnitude.max(f64::MIN_POSITIVE);
    if sum < -roundoff { return Err(invalid("prepared resistance integral is negative beyond roundoff")); }
    Ok(sum.max(0.0))
}

#[cfg(test)]
mod tests;
