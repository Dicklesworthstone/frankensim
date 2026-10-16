//! Work-conjugate real velocity ports on the existing P1/P0 acoustic solution.
//!
//! Unit generalized VELOCITY drives one source column. Z_ij = integral b_i p_j
//! is the resisting impedance under exp(-i*w*t); body force is -Z*v. The same
//! returned solutions supply receivers, without substituting panel mean pressure
//! into a P0 observer or relabeling the auxiliary density as normal velocity.
use super::{HelmholtzError, Medium, Prepared, Solution, C64, bad, finite, poll, zeroed};

/// One complete velocity-driven source batch and its unchanged projected work.
/// Column j of the row-major impedance and solution j have the same source.
/// All densities stay bound to the source geometry, frequency and medium.
pub struct VelocityPortSolution<'p, 's> {
    solutions: Vec<Solution<'p, 's>>,
    impedance: Vec<C64>,
}
impl<'p, 's> VelocityPortSolution<'p, 's> {
    pub fn solutions(&self) -> &[Solution<'p, 's>] { &self.solutions }
    pub fn impedance(&self) -> &[C64] { &self.impedance }
    pub fn ports(&self) -> usize { self.solutions.len() }
}

impl<'s> Prepared<'s> {
    /// Bound the original source solve AND the full mechanical work projection.
    /// This is a cold-work estimate, not wall time. The inherited source-count
    /// ceiling still applies; dependent/zero/signed columns are never removed.
    pub fn velocity_port_work(&self, ports: usize) -> Result<u64, HelmholtzError> {
        let solve = self.frequency_work(ports)?;
        let projection = (ports as u128) * (ports as u128)
            * (self.facets.len() as u128) * 8;
        let extra = u64::try_from(projection).map_err(|_| bad("velocity-port projection work overflow"))?;
        solve.checked_add(extra).ok_or_else(|| bad("velocity-port work overflow"))
    }

    /// Solve and project one immutable real source basis, once per frequency.
    /// `rows[j][panel]` converts generalized velocity j to outward panel velocity.
    /// Rows are physical shape coefficients, not normalized acoustic amplitudes.
    /// Work uses the same rows and actual panel areas, with NO pressure repair,
    /// Hermitian symmetrization, positive diagonal loading or discarded modes.
    ///
    /// An acceleration response is obtained by dividing a solution's pressure
    /// by -i*w AFTER observation; it must not be mistaken for this velocity Z.
    /// The inherited flux/backward-residual gates remain authoritative. The
    /// discrete formulation's passivity does not certify continuum accuracy.
    /// All row/budget admission precedes solving; late cancellation or arithmetic
    /// refusal returns no partially projected batch and does not mutate self.
    pub fn solve_velocity_ports<'p>(&'p self, k: f64, medium: Medium,
        rows: &[&[f64]], mut cancelled: impl FnMut() -> bool)
        -> Result<VelocityPortSolution<'p, 's>, HelmholtzError>
    {
        poll(&mut cancelled)?;
        let n = rows.len();
        if self.velocity_port_work(n)? > self.budget.maximum_work {
            return Err(bad("velocity-port solve/projection exceeds declared work budget"));
        }
        if rows.iter().any(|row| row.len() != self.facets.len()
            || row.iter().any(|v| !v.is_finite())) {
            return Err(bad("velocity ports require complete finite real panel rows"));
        }
        let mut fields = Vec::new();
        fields.try_reserve_exact(n).map_err(|_| HelmholtzError::AllocationFailed {
            what: "velocity-port source rows" })?;
        for row in rows {
            poll(&mut cancelled)?;
            let mut field = zeroed(row.len(), C64::ZERO)?;
            for (v, &b) in field.iter_mut().zip(*row) { *v = C64::from_re(b); }
            fields.push(field);
        }
        let refs: Vec<_> = fields.iter().map(Vec::as_slice).collect();
        let solutions = self.solve_batch(k, medium, &refs, &mut cancelled)?;
        let mut impedance = zeroed(n*n, C64::ZERO)?;
        for (j, solution) in solutions.iter().enumerate() {
            for (i, row) in rows.iter().enumerate() {
                poll(&mut cancelled)?;
                let mut force = C64::ZERO;
                for ((&b, &p), &area) in row.iter().zip(&solution.pressure).zip(self.surface.areas()) {
                    // Preserve the transpose of the physical velocity map, not
                    // complex conjugation of the acoustic pressure solution.
                    force = force + p.scale(area*b);
                }
                if !finite(force) { return Err(bad("nonfinite velocity-port work projection")); }
                impedance[i*n+j] = force;
            }
        }
        poll(&mut cancelled)?;
        Ok(VelocityPortSolution { solutions, impedance })
    }
}

#[cfg(test)]
mod tests;
