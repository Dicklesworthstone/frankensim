//! Original material relaxation in the two-sided pressure-jump screen.
//!
//! The material supplies k(omega), while inviscid momentum still gives the
//! physical drive H(k)*[p] = i*rho*omega*G^T*rate. The same exact P1 pairing G
//! returns every self/mutual effort. No complex sound speed replaces omega,
//! and no extra face, phenomenological damping or absorption filter is added.
//!
//! This uses the original bounded cosine/antipodal continuation of H. Its
//! finite quadrature must converge independently; the signed boundary work is
//! neither repaired nor certified positive for every complex k. In a relaxing
//! medium it contains exterior molecular heat as well as outgoing energy, and
//! is not far-field flux or a causal time-domain realization.
use super::*;
use crate::helmholtz::molecular::{MolecularError, MolecularExterior, PreparedMolecularIncident};

mod receivers;
pub use receivers::PreparedMolecularScreenReceivers;
mod incident;
pub use incident::{
    MolecularScreenIllumination, MolecularScreenIlluminationLimits,
    PreparedMolecularScreenIlluminationReceivers,
};

/// Whole source-batch limits in addition to the original geometry admission.
#[derive(Debug, Clone, Copy)]
pub struct MolecularScreenPortLimits {
    /// Complete impedance/jump extents; maximum_work admits the ENTIRE
    /// values assembly, LU, source solves, residuals and work projection.
    pub ports: ScreenPortLimits,
    /// One H, its centered growth chart and values-only assembly scratch.
    pub operator: ScreenComplexValueLimits,
    /// Live solve payload after assembly: H+LU, exact weak pairing, all jumps,
    /// impedance, RHS/solution/pivot vectors and blocked-LU descriptors.
    /// Assembly scratch has its separate nonoverlapping `operator` bound.
    /// Excludes geometry/source ownership, later receiver output, allocator
    /// metadata and thread stacks; these are not process RSS ceilings.
    pub maximum_workspace_bytes: usize,
}

/// Complete original pressure-jump columns at one real material frequency.
/// Both the immutable source and the original prepared screen remain bound.
pub struct MolecularScreenVelocityPortSolution<'p, 's> {
    prepared: &'p PreparedScreen<'s>,
    source: MolecularExterior,
    omega: f64,
    k: C64,
    jumps: Vec<Vec<C64>>,
    impedance: Vec<C64>,
    maximum_backward_residual: f64,
    frequency_work: u64,
}

impl<'s> PreparedScreen<'s> {
    /// Conservative scalar-work units for this complete complex-k batch.
    /// Values-only assembly uses four times the original non-LU assembly
    /// allowance; all remaining terms include original pressure residuals and
    /// exact full P1 work projection. This is not a measured FLOP count.
    pub fn molecular_velocity_port_work(&self, ports: usize) -> Result<u64, HelmholtzError> {
        if ports == 0 || ports > crate::helmholtz::MAX_RADIATION_FIELDS_PER_BATCH {
            return Err(bad("molecular screen port count must be in 1..=256"));
        }
        let (p, n, f, v) = (
            ports as u128,
            self.n as u128,
            self.facets.len() as u128,
            self.surface.vertices.len() as u128,
        );
        let work = 4 * (self.frequency_work() - n * n * n)
            + 32 * n * n * n
            + 128 * p * n * n
            + 128 * p * f
            + 64 * p * v
            + 64 * p * p * v;
        u64::try_from(work).map_err(|_| bad("molecular screen batch work overflows"))
    }

    /// Solve every facet-local P1 velocity column against one original H and
    /// one LU. The source material is evaluated once at the physical omega;
    /// pressure has the original zero free-edge values, velocity does not.
    ///
    /// Both sides contact the SAME homogeneous relaxing gas. Shear viscosity,
    /// thermal layers, mean flow, a cavity backing and independent face coats
    /// require other owners. Failed admission, cancellation or a bad residual
    /// returns no partial source fields or impedance matrix.
    pub fn solve_molecular_velocity_ports<'p>(
        &'p self,
        source: &MolecularExterior,
        omega: f64,
        rows: &[&[[f64; 3]]],
        limits: MolecularScreenPortLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularScreenVelocityPortSolution<'p, 's>, MolecularError> {
        self.solve_molecular_ports_with(
            source,
            omega,
            rows,
            limits,
            None,
            |_, _, _, _, _| Ok(()),
            cancelled,
        )
        .map(|(motion, ())| motion)
    }

    /// Keep the original H/LU alive while an admitted companion batch solves
    /// its own weak drives. All callbacks finish before either result escapes.
    #[allow(clippy::too_many_arguments)]
    fn solve_molecular_ports_with<'p, T, C: FnMut() -> bool>(
        &'p self,
        source: &MolecularExterior,
        omega: f64,
        rows: &[&[[f64; 3]]],
        limits: MolecularScreenPortLimits,
        prepared_incident: Option<&PreparedMolecularIncident<'_>>,
        mut finish: impl FnMut(
            Vec<Vec<f64>>,
            &[C64],
            &fs_la::eigen_complex::LuComplex,
            f64,
            &mut C,
        ) -> Result<T, MolecularError>,
        mut cancelled: C,
    ) -> Result<(MolecularScreenVelocityPortSolution<'p, 's>, T), MolecularError> {
        poll(&mut cancelled)?;
        let ports = rows.len();
        let count = limits.ports;
        let band = source.band_rad_s();
        if !omega.is_finite()
            || omega < band[0]
            || omega > band[1]
            || count.maximum_ports == 0
            || count.maximum_ports > crate::helmholtz::MAX_RADIATION_FIELDS_PER_BATCH
            || ports > count.maximum_ports
        {
            return Err(bad(
                "molecular screen frequency or source count is outside its declaration",
            )
            .into());
        }
        let frequency_work = self.molecular_velocity_port_work(ports)?;
        let entries = ports
            .checked_mul(ports)
            .ok_or_else(|| bad("molecular screen impedance extent overflows"))?;
        let jumps = ports
            .checked_mul(self.surface.vertices.len())
            .ok_or_else(|| bad("molecular screen jump extent overflows"))?;
        let bytes = self.molecular_port_payload(ports)?;
        if entries > count.maximum_matrix_entries
            || jumps > count.maximum_jump_entries
            || frequency_work > count.maximum_work
            || frequency_work > self.options.maximum_work
            || bytes > limits.maximum_workspace_bytes as u128
        {
            return Err(bad("complete molecular screen batch exceeds declared resources").into());
        }
        for row in rows {
            poll(&mut cancelled)?;
            if row.len() != self.facets.len() || row.iter().flatten().any(|x| !x.is_finite()) {
                return Err(
                    bad("molecular screen needs complete finite facet-corner velocities").into(),
                );
            }
        }
        let rho_omega = source.source().gas().density * omega;
        let k = if let Some(incident) = prepared_incident {
            // Only the original source's private illumination call supplies
            // this owner; its evaluated gas and frequency remain immutable.
            if incident.angular_frequency_rad_s().to_bits() != omega.to_bits() {
                return Err(bad("molecular incident frequency differs from its screen").into());
            }
            incident.wavenumber_per_m()
        } else {
            source.source().relaxing_wavenumber(omega)?
        };
        if !finite(k) || k.re <= 0.0 || k.im < 0.0 || !rho_omega.is_finite() || rho_omega <= 0.0 {
            return Err(bad(
                "molecular screen requires passive outgoing k and representable rho*omega",
            )
            .into());
        }
        // Assembly scratch is released before weak source columns or LU exist.
        // The value-only specialization keeps the original real-axis H bits.
        let h = self.complex_matrix(k, limits.operator, &mut cancelled)?;
        let weak = self.velocity_port_pairing(
            rows,
            count.maximum_jump_entries,
            count.maximum_work,
            &mut cancelled,
        )?;
        let lu = fs_la::eigen_complex::with_lu_workers(self.options.workers, || {
            lu_complex_with_cancel(&h, self.n, &mut cancelled)
        })?;
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(ports)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "molecular screen pressure-jump columns",
            })?;
        let mut maximum_backward_residual = 0.0_f64;
        for row in &weak {
            poll(&mut cancelled)?;
            let mut rhs = zeroed(self.n, C64::ZERO)?;
            for (i, &value) in row.iter().enumerate() {
                if i % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                rhs[i] = C64::new(0.0, rho_omega * value);
                if !finite(rhs[i]) || (value != 0.0 && rhs[i].im == 0.0) {
                    return Err(bad("molecular screen physical drive is unrepresentable").into());
                }
            }
            let (jump, backward) = self.solve_molecular_rhs(&rhs, &h, &lu, &mut cancelled)?;
            maximum_backward_residual = maximum_backward_residual.max(backward);
            fields.push(jump);
        }
        let mut impedance = zeroed(entries, C64::ZERO)?;
        for (j, jump) in fields.iter().enumerate() {
            for (i, row) in weak.iter().enumerate() {
                poll(&mut cancelled)?;
                impedance[i * ports + j] = self.molecular_jump_effort(row, jump, &mut cancelled)?;
            }
        }
        let extra = finish(weak, &h, &lu, rho_omega, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok((
            MolecularScreenVelocityPortSolution {
                prepared: self,
                source: source.clone(),
                omega,
                k,
                jumps: fields,
                impedance,
                maximum_backward_residual,
                frequency_work,
            },
            extra,
        ))
    }

    fn molecular_port_payload(&self, ports: usize) -> Result<u128, HelmholtzError> {
        let (n, p, v) = (
            self.n as u128,
            ports as u128,
            self.surface.vertices.len() as u128,
        );
        let scratch =
            fs_la::eigen_complex::lu_complex_blocked_workspace_bytes(self.n, self.options.workers)
                .ok_or_else(|| bad("molecular screen LU scratch extent overflows"))?;
        Ok(32 * n * n
            + 8 * p * n
            + 16 * p * v
            + 16 * p * p
            + 40 * n
            + 128 * p
            + 4096
            + scratch as u128)
    }

    /// Original physical residual and nodal reconstruction shared by the
    /// real mechanical columns and complex blocked-illumination drives.
    fn solve_molecular_rhs(
        &self,
        rhs: &[C64],
        h: &[C64],
        lu: &fs_la::eigen_complex::LuComplex,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, f64), MolecularError> {
        let mut x = rhs.to_vec();
        lu.solve_with_cancel(&mut x, &mut *cancelled)?;
        let mut residual = 0.0_f64;
        let mut scale = 0.0_f64;
        for i in 0..self.n {
            poll(cancelled)?;
            let mut lhs = C64::ZERO;
            let mut envelope = rhs[i].abs();
            for j in 0..self.n {
                lhs = lhs + h[i * self.n + j] * x[j];
                envelope += h[i * self.n + j].abs() * x[j].abs();
            }
            if !finite(lhs) || !envelope.is_finite() {
                return Err(bad("molecular screen original residual is unrepresentable").into());
            }
            residual = residual.max((lhs - rhs[i]).abs());
            scale = scale.max(envelope);
        }
        let backward = residual / scale.max(f64::MIN_POSITIVE);
        if !backward.is_finite() || backward > 1e-10 || x.iter().any(|&z| !finite(z)) {
            return Err(bad("molecular screen original pressure residual failed").into());
        }
        let mut jump = zeroed(self.surface.vertices.len(), C64::ZERO)?;
        for (i, dof) in self.dof.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            if let Some(u) = dof {
                jump[i] = x[*u];
            }
        }
        Ok((jump, backward))
    }

    fn molecular_jump_effort(
        &self,
        row: &[f64],
        jump: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<C64, HelmholtzError> {
        let mut effort = Sum::default();
        for (vertex, dof) in self.dof.iter().enumerate() {
            if vertex % 256 == 0 {
                poll(cancelled)?;
            }
            if let Some(u) = dof {
                effort.add(jump[vertex].scale(row[*u]));
            }
        }
        let value = effort.value();
        if !finite(value) {
            return Err(bad("molecular screen full work projection overflows"));
        }
        Ok(value)
    }
}

impl MolecularScreenVelocityPortSolution<'_, '_> {
    pub fn source(&self) -> &MolecularExterior {
        &self.source
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn wavenumber_per_m(&self) -> C64 {
        self.k
    }
    pub fn ports(&self) -> usize {
        self.jumps.len()
    }
    pub fn impedance(&self) -> &[C64] {
        &self.impedance
    }
    pub fn pressure_jumps_pa(&self) -> &[Vec<C64>] {
        &self.jumps
    }
    pub fn frequency_work(&self) -> u64 {
        self.frequency_work
    }
    /// Original column solve residual only: no spatial error or coherent
    /// cancellation-relative field accuracy is inferred from this number.
    pub fn maximum_backward_residual(&self) -> f64 {
        self.maximum_backward_residual
    }

    fn admit_rates(&self, rates: &[C64]) -> Result<(), HelmholtzError> {
        if rates.len() != self.ports() || rates.iter().any(|&z| !finite(z)) {
            return Err(bad(
                "molecular screen requires every finite source coordinate",
            ));
        }
        Ok(())
    }

    /// Every original modal resisting effort Z*rate; fluid reaction is its
    /// negative. This bounded matrix has at most 256 squared entries.
    pub fn resisting_efforts(&self, rates: &[C64]) -> Result<Vec<C64>, HelmholtzError> {
        self.admit_rates(rates)?;
        let mut efforts = zeroed(self.ports(), C64::ZERO)?;
        for (effort, row) in efforts.iter_mut().zip(self.impedance.chunks(self.ports())) {
            let mut sum = Sum::default();
            for (&z, &rate) in row.iter().zip(rates) {
                sum.add(z * rate);
            }
            *effort = sum.value();
        }
        if efforts.iter().any(|&z| !finite(z)) {
            return Err(bad("molecular screen coherent efforts overflow"));
        }
        Ok(efforts)
    }

    /// Signed cycle-mean work entering the exterior, .5 Re(rate^H Z rate).
    /// Gas absorption and outgoing radiation are not split or counted twice.
    pub fn boundary_work_w(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, HelmholtzError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        let mut work = Sum::default();
        for (i, row) in self.impedance.chunks(self.ports()).enumerate() {
            poll(&mut cancelled)?;
            let mut effort = Sum::default();
            for (&z, &rate) in row.iter().zip(rates) {
                effort.add(z * rate);
            }
            work.add(rates[i].conj() * effort.value());
        }
        let complex_work = work.value();
        let work = 0.5 * complex_work.re;
        if !finite(complex_work) || !work.is_finite() {
            return Err(bad("molecular screen coherent work overflows"));
        }
        poll(&mut cancelled)?;
        Ok(work)
    }

    /// Combine the actual complex rates BEFORE observing the two-sided field.
    /// Free-edge pressure remains exactly zero with nonzero edge velocity.
    pub fn pressure_jump(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        let mut jump = zeroed(self.prepared.surface.vertices.len(), C64::ZERO)?;
        for (v, out) in jump.iter_mut().enumerate() {
            if v % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let mut sum = Sum::default();
            for (field, &rate) in self.jumps.iter().zip(rates) {
                sum.add(field[v] * rate);
            }
            *out = sum.value();
            if !finite(*out) {
                return Err(bad("molecular screen coherent pressure jump overflows"));
            }
        }
        poll(&mut cancelled)?;
        Ok(jump)
    }
}

#[derive(Default)]
struct Sum {
    value: C64,
    correction: C64,
}
impl Sum {
    fn add(&mut self, value: C64) {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
    }
    fn value(self) -> C64 {
        self.value
    }
}

#[cfg(test)]
pub(super) mod tests;
