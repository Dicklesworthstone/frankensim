//! Complete work-conjugate velocity ports on the original open P1 screen.
//!
//! A curved facet mesh needs facet-local normal traces: a shared nodal vector
//! has a different dot product with each adjacent facet normal. For real trace
//! b_r and pressure jump p_s, Z_rs = integral b_r p_s dS. The exact P1 triangle
//! mass matrix A/12 * [[2,1,1],[1,2,1],[1,1,2]] supplies BOTH the acoustic RHS
//! and its work transpose. Only pressure has zero free-edge nodal values.
use super::*;

/// Additional complete-batch limits. The prepared geometry's limits still apply.
/// These are logical payload/work bounds, not allocator overhead or process RSS.
#[derive(Debug, Clone, Copy)]
pub struct ScreenPortLimits {
    /// Complete source count, at most MAX_RADIATION_FIELDS_PER_BATCH (256).
    pub maximum_ports: usize,
    /// Every row-major impedance entry, including all self and mutual terms.
    pub maximum_matrix_entries: usize,
    /// Every retained nodal pressure jump across all source columns. The real
    /// weak-drive scratch has no more entries than this same bound.
    pub maximum_jump_entries: usize,
    /// One assembly/factor, all source solves/audits and full work projection.
    /// The same total must also fit OpenOptions::maximum_work.
    pub maximum_work: u64,
}

/// Original source columns and their unchanged work matrix under exp(-i*w*t).
/// Fluid force on the structure is -Z*rate. No modes, signs or cross terms are
/// removed, and no symmetrization or passive fitting changes this matrix.
pub struct ScreenVelocityPortSolution<'p, 's> {
    solutions: Vec<JumpSolution<'p, 's>>,
    impedance: Vec<C64>,
    frequency_work: u64,
}

impl<'s> PreparedScreen<'s> {
    /// Conservative scalar-work units for the entire rigid source batch,
    /// including its exact P1 work projection. This is not a timing estimate.
    pub fn velocity_port_work(&self, ports: usize) -> Result<u64, HelmholtzError> {
        if ports == 0 || ports > crate::helmholtz::MAX_RADIATION_FIELDS_PER_BATCH {
            return Err(bad("open velocity-port count must be in 1..=256"));
        }
        let (p, n, f, v) = (
            ports as u128,
            self.n as u128,
            self.facets.len() as u128,
            self.surface.vertices.len() as u128,
        );
        let work =
            self.frequency_work() + 128 * p * n * n + 128 * p * f + 64 * p * v + 8 * p * p * v;
        u64::try_from(work).map_err(|_| bad("open velocity-port work overflow"))
    }

    /// Solve every physical facet-local P1 velocity against ONE original
    /// hypersingular assembly and LU. rows[mode][facet][corner] is the signed
    /// normal-velocity coefficient on that oriented facet. It need not agree
    /// across a folded edge, and it need not vanish at the free boundary.
    ///
    /// The structural screen is impermeable, with both faces in this medium.
    /// A sheet admittance would be an additional mechanical law and is not
    /// silently attached here. The pressure jump already represents both faces;
    /// no factor two, coincident second surface or artificial thickness enters.
    /// Source residuals and receiver quadrature are those of the original owner.
    /// A batch failure returns no partial source fields or impedance matrix.
    pub fn solve_velocity_ports<'p>(
        &'p self,
        k: f64,
        medium: Medium,
        rows: &[&[[f64; 3]]],
        limits: ScreenPortLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScreenVelocityPortSolution<'p, 's>, HelmholtzError> {
        poll(&mut cancelled)?;
        let ports = rows.len();
        let frequency_work = self.velocity_port_work(ports)?;
        if !k.is_finite()
            || k <= 0.0
            || k > self.highest_k
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || limits.maximum_ports == 0
            || limits.maximum_ports > crate::helmholtz::MAX_RADIATION_FIELDS_PER_BATCH
            || ports > limits.maximum_ports
        {
            return Err(bad(
                "invalid open velocity-port frequency, medium or source limit",
            ));
        }
        let entries = ports
            .checked_mul(ports)
            .ok_or_else(|| bad("open velocity-port matrix extent overflow"))?;
        let jumps = ports
            .checked_mul(self.surface.vertices.len())
            .ok_or_else(|| bad("open velocity-port jump extent overflow"))?;
        if entries > limits.maximum_matrix_entries
            || jumps > limits.maximum_jump_entries
            || frequency_work > limits.maximum_work
            || frequency_work > self.options.maximum_work
        {
            return Err(bad(
                "complete open velocity-port batch exceeds declared resources",
            ));
        }
        for row in rows {
            poll(&mut cancelled)?;
            if row.len() != self.facets.len() || row.iter().flatten().any(|v| !v.is_finite()) {
                return Err(bad(
                    "open velocity ports require complete finite facet-corner rows",
                ));
            }
        }
        // Retain the actual test/drive pairing once, rather than repeating nine
        // facet products for every mode pair or inventing a lumped mass rule.
        let mut weak = Vec::new();
        weak.try_reserve_exact(ports)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "open velocity-port weak drives",
            })?;
        for row in rows {
            poll(&mut cancelled)?;
            let mut drive = zeroed(self.n, 0.0)?;
            for (index, (f, values)) in self.facets.iter().zip(*row).enumerate() {
                if index % 32 == 0 {
                    poll(&mut cancelled)?;
                }
                for a in 0..3 {
                    if let Some(u) = self.dof[f.ids[a]] {
                        for b in 0..3 {
                            drive[u] +=
                                values[b] * (f.area * if a == b { 1.0 / 6.0 } else { 1.0 / 12.0 });
                        }
                    }
                }
            }
            if drive.iter().any(|v| !v.is_finite()) {
                return Err(bad("open velocity-port weak projection overflow"));
            }
            weak.push(drive);
        }
        poll(&mut cancelled)?;
        let (h, lu) = self.factor_wave(k, true, Admittance::Rigid, &mut cancelled)?;
        let mut solutions = Vec::new();
        solutions
            .try_reserve_exact(ports)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "open velocity-port source fields",
            })?;
        let mut drive = zeroed(self.n, C64::ZERO)?;
        for row in &weak {
            poll(&mut cancelled)?;
            for (v, &value) in drive.iter_mut().zip(row) {
                *v = C64::from_re(value);
            }
            solutions.push(self.finish_weak_drive(
                k,
                true,
                medium,
                Admittance::Rigid,
                &drive,
                &h,
                &lu,
                &mut cancelled,
            )?);
        }
        let mut impedance = zeroed(entries, C64::ZERO)?;
        for (j, solution) in solutions.iter().enumerate() {
            for (i, row) in weak.iter().enumerate() {
                poll(&mut cancelled)?;
                let mut force = C64::ZERO;
                for (v, &dof) in self.dof.iter().enumerate() {
                    if let Some(u) = dof {
                        force = force + solution.jump[v].scale(row[u]);
                    }
                }
                if !finite(force) {
                    return Err(bad("open velocity-port work projection overflow"));
                }
                impedance[i * ports + j] = force;
            }
        }
        poll(&mut cancelled)?;
        Ok(ScreenVelocityPortSolution {
            solutions,
            impedance,
            frequency_work,
        })
    }
}

impl<'p, 's> ScreenVelocityPortSolution<'p, 's> {
    pub fn solutions(&self) -> &[JumpSolution<'p, 's>] {
        &self.solutions
    }
    pub fn impedance(&self) -> &[C64] {
        &self.impedance
    }
    pub fn ports(&self) -> usize {
        self.solutions.len()
    }
    pub fn frequency_work(&self) -> u64 {
        self.frequency_work
    }
    /// Largest original per-source backward residual; not an error bound for
    /// a nearly cancelling coherent combination or the continuum field.
    pub fn maximum_backward_residual(&self) -> f64 {
        self.solutions
            .iter()
            .map(|s| s.backward_residual)
            .fold(0.0, f64::max)
    }

    fn admit_rates(&self, rates: &[C64]) -> Result<(), HelmholtzError> {
        if rates.len() != self.ports() || rates.iter().any(|v| !finite(*v)) {
            return Err(bad(
                "open velocity-port rates must match every finite source coordinate",
            ));
        }
        Ok(())
    }

    /// Complete resisting efforts Z*rate in the original source coordinates.
    pub fn resisting_efforts(&self, rates: &[C64]) -> Result<Vec<C64>, HelmholtzError> {
        self.admit_rates(rates)?;
        let mut result = zeroed(self.ports(), C64::ZERO)?;
        for (out, row) in result.iter_mut().zip(self.impedance.chunks(self.ports())) {
            for (&z, &rate) in row.iter().zip(rates) {
                *out = *out + z * rate;
            }
        }
        if result.iter().any(|v| !finite(*v)) {
            return Err(bad("open velocity-port effort overflow"));
        }
        Ok(result)
    }

    /// Peak-amplitude cycle-average source work, .5 Re(rate^H Z rate) [W].
    /// This includes coherent cross terms and retains its original signed value.
    pub fn input_power_w(&self, rates: &[C64]) -> Result<f64, HelmholtzError> {
        let efforts = self.resisting_efforts(rates)?;
        let power = 0.5
            * rates
                .iter()
                .zip(efforts)
                .map(|(&rate, force)| (rate.conj() * force).re)
                .sum::<f64>();
        if !power.is_finite() {
            return Err(bad("open velocity-port work overflow"));
        }
        Ok(power)
    }

    /// Coherent nodal pressure jump, in the original geometry's vertex order.
    /// Free-edge values stay zero even when the physical velocity there is not.
    pub fn pressure_jump(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        let mut jump = zeroed(self.solutions[0].jump.len(), C64::ZERO)?;
        for (solution, &rate) in self.solutions.iter().zip(rates) {
            for (i, (out, &source)) in jump.iter_mut().zip(&solution.jump).enumerate() {
                if i % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                *out = *out + source * rate;
            }
        }
        if jump.iter().any(|v| !finite(*v)) {
            return Err(bad("coherent open pressure jump overflow"));
        }
        poll(&mut cancelled)?;
        Ok(jump)
    }

    /// Observe the combined field once, with the original adaptive P1 kernels.
    /// Complex source amplitudes combine before quadrature, never as powers.
    pub fn pressure_at(
        &self,
        rates: &[C64],
        points: &[Point],
        separation: f64,
        options: crate::near_field::Options,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        poll(&mut cancelled)?;
        self.admit_rates(rates)?;
        options.validate()?;
        let source = &self.solutions[0];
        if !separation.is_finite()
            || separation <= 0.0
            || points.iter().flatten().any(|v| !v.is_finite())
            || points
                .len()
                .checked_mul(source.prepared.facets.len())
                .and_then(|v| v.checked_mul(80))
                .is_none_or(|v| v > options.maximum_kernel_evaluations)
        {
            return Err(bad(
                "invalid or over-budget coherent open-screen receiver request",
            ));
        }
        let jump = self.pressure_jump(rates, &mut cancelled)?;
        source.pressure_at_with_jump(&jump, points, separation, options, cancelled)
    }

    /// Coherent r*exp(-ikr)*p(r*d), observed from either side of the screen.
    pub fn far_field(
        &self,
        rates: &[C64],
        direction: Point,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<C64, HelmholtzError> {
        let jump = self.pressure_jump(rates, &mut cancelled)?;
        let field = self.solutions[0].far_field_with_jump(&jump, direction)?;
        poll(&mut cancelled)?;
        if !finite(field) {
            return Err(bad("coherent open far-field overflow"));
        }
        Ok(field)
    }

    /// Original antipodal outgoing-power rule applied to the SAME coherent
    /// field, including all interference. Its angular chart remains unchanged.
    pub fn radiated_power_w(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, HelmholtzError> {
        let jump = self.pressure_jump(rates, &mut cancelled)?;
        let source = &self.solutions[0];
        let mut sum = 0.0;
        for (direction, weight) in directions(source.prepared.options.angular_refinement) {
            poll(&mut cancelled)?;
            sum += weight
                * (source.far_field_with_jump(&jump, direction)?.abs().powi(2)
                    + source
                        .far_field_with_jump(&jump, direction.map(|v| -v))?
                        .abs()
                        .powi(2));
        }
        let power = sum / (2.0 * source.medium.density * source.medium.sound_speed);
        poll(&mut cancelled)?;
        if !power.is_finite() {
            return Err(bad("coherent open radiated power overflow"));
        }
        Ok(power)
    }
}
