//! Sparse pressure-driven compressibility memories on original volume forms.
use fs_math::{dd::Dd, det};
use fs_sparse::direct::{
    DirectLimits, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt,
};
use fs_sparse::{Coo, Csr};
use std::cell::RefCell;

/// Original beta(s)=beta_inf+sum delta_beta/(1+s*tau) data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkMemoryArm {
    /// Positive relaxed compressibility increment [1/Pa].
    pub delta_compressibility_pa_inv: f64,
    /// Positive pressure-driven time [s], in tau*q_dot+q=p.
    pub relaxation_time_s: f64,
}
/// Bounds original fields, retained region forms and native sparse preparation.
#[derive(Clone, Copy, Debug)]
pub struct BulkTimeLimits {
    /// All original pressure coordinates, including constrained ones.
    pub maximum_dofs: usize,
    /// Disjoint original material regions carrying memories.
    pub maximum_regions: usize,
    /// Total pressure-memory arms across the regions.
    pub maximum_arms: usize,
    /// Sum of the complete original region coordinate count for every arm.
    pub maximum_memory_dofs: usize,
    /// Aggregate retained physical CSR entries; independently also bounds
    /// each region/effective-operator COO staging buffer. This is not a peak
    /// allocator-byte cap. Original geometry and caller state are additional.
    pub maximum_matrix_entries: usize,
    /// Original LDLT factor/front/scratch payload cap. Retained CSR, symbolic
    /// analysis, field/work vectors and allocator overhead are additional.
    pub maximum_factor_bytes: usize,
    /// Original sparse symbolic-work, numeric-work and front-dimension caps.
    pub sparse: DirectLimits,
    /// Componentwise original continuity, kinematic and memory ODE tolerance.
    pub maximum_relative_residual: f64,
    /// Relative physical work allowance, with evaluated roundoff added.
    pub energy_balance_tolerance: f64,
}
impl Default for BulkTimeLimits {
    fn default() -> Self {
        Self {
            maximum_dofs: 16384,
            maximum_regions: 256,
            maximum_arms: 8192,
            maximum_memory_dofs: 2_000_000,
            maximum_matrix_entries: 8_000_000,
            maximum_factor_bytes: 128 * 1024 * 1024,
            sparse: DirectLimits {
                max_symbolic_entries: 8_000_000,
                max_analysis_work: 100_000_000,
                max_numeric_work: 1_000_000_000,
                max_front_dim: 2048,
            },
            maximum_relative_residual: 5e-11,
            energy_balance_tolerance: 5e-10,
        }
    }
}
#[derive(Debug)]
pub enum BulkTimeError {
    Invalid(&'static str),
    Limit {
        what: &'static str,
        requested: usize,
        maximum: usize,
    },
    Factor(LdltError),
    Cancelled,
    Residual {
        observed: f64,
        maximum: f64,
    },
    EnergyBalance {
        account: &'static str,
        residual_j: f64,
        maximum_j: f64,
    },
}
impl core::fmt::Display for BulkTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bulk pressure memory: {self:?}")
    }
}
impl std::error::Error for BulkTimeError {}
impl From<LdltError> for BulkTimeError {
    fn from(e: LdltError) -> Self {
        if matches!(e, LdltError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Factor(e)
        }
    }
}
type Result<T> = core::result::Result<T, BulkTimeError>;
const ROUND: f64 = 4096.0 * f64::EPSILON;

pub(super) struct Region {
    pub(super) dofs: Vec<usize>,
    /// One consistent unit-weight volume mass, shared by all region arms.
    pub(super) mass: Csr,
    pub(super) arms: Vec<BulkMemoryArm>,
}
struct Arm {
    region: usize,
    source: BulkMemoryArm,
    retention: f64,
    response: f64,
}
pub(super) struct PreparedBulkDynamics {
    mass: Csr,
    stiffness: Csr,
    conductance: Csr,
    regions: Vec<Region>,
    arms: Vec<Arm>,
    free: Vec<usize>,
    factor: LdltFactor,
    scale: Vec<f64>,
    duration: f64,
    limits: BulkTimeLimits,
}
/// Complete pressure/potential coefficients and every original volume history.
#[derive(Clone)]
pub struct BulkPressureState<'p> {
    owner: &'p PreparedBulkDynamics,
    potential: Vec<f64>,
    pressure: Vec<f64>,
    memories: Vec<Vec<f64>>,
}
impl core::fmt::Debug for BulkPressureState<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BulkPressureState")
            .field("pressure_dofs", &self.pressure.len())
            .field("memory_arms", &self.memories.len())
            .finish_non_exhaustive()
    }
}
impl BulkPressureState<'_> {
    /// Original hierarchical pressure coefficients [Pa], including bubbles.
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure
    }
    /// Original potential coefficients [Pa s]; velocity is -grad(psi)/rho.
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        &self.potential
    }
    /// Complete region coefficients [Pa]. Arms are flattened in input region
    /// order, then authored arm order; coefficients follow `arm_dofs`.
    pub fn memory_pressure_pa(&self, arm: usize) -> Option<&[f64]> {
        self.memories.get(arm).map(Vec::as_slice)
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct BulkTimeEnergy {
    /// Original physical psi^T*K_rho*psi/2.
    pub kinetic_j: f64,
    /// Instantaneous pressure storage p^T*M_beta_inf*p/2.
    pub frozen_compression_j: f64,
    /// Sum of delta_beta*q^T*M_region*q/2 over all original volume memories.
    pub molecular_storage_j: f64,
}
impl BulkTimeEnergy {
    pub fn total_j(self) -> f64 {
        self.kinetic_j + self.frozen_compression_j + self.molecular_storage_j
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct BulkMemoryWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub energy_change_j: f64,
    /// Work from the physical pressure into this original memory.
    pub pressure_work_j: f64,
    /// Independently integrated h*delta_beta*(p_mid-q_mid)^2/tau.
    pub dissipation_j: f64,
    pub maximum_relative_residual: f64,
}
impl BulkMemoryWork {
    pub fn balance_residual_j(self) -> f64 {
        difference(
            self.energy_change_j,
            self.dissipation_j,
            self.pressure_work_j,
        )
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct BulkTimeWork {
    pub initial_energy: BulkTimeEnergy,
    pub final_energy: BulkTimeEnergy,
    pub acoustic_energy_change_j: f64,
    pub memory_energy_change_j: f64,
    pub input_work_j: f64,
    pub boundary_dissipation_j: f64,
    /// Subtotals already included in input_work_j.
    pub boundary_source_work_j: f64,
    pub volume_source_work_j: f64,
    pub molecular_dissipation_j: f64,
    pub pressure_to_memory_work_j: f64,
    pub maximum_relative_residual: f64,
}
impl BulkTimeWork {
    pub fn acoustic_balance_residual_j(self) -> f64 {
        difference(
            self.acoustic_energy_change_j,
            self.boundary_dissipation_j + self.pressure_to_memory_work_j,
            self.input_work_j,
        )
    }
    pub fn memory_balance_residual_j(self) -> f64 {
        difference(
            self.memory_energy_change_j,
            self.molecular_dissipation_j,
            self.pressure_to_memory_work_j,
        )
    }
    pub fn balance_residual_j(self) -> f64 {
        difference(
            self.acoustic_energy_change_j + self.memory_energy_change_j,
            self.boundary_dissipation_j + self.molecular_dissipation_j,
            self.input_work_j,
        )
    }
}
#[derive(Debug)]
pub struct BulkTimeStep<'p> {
    pub state: BulkPressureState<'p>,
    pub work: BulkTimeWork,
    pub memories: Vec<BulkMemoryWork>,
}

impl PreparedBulkDynamics {
    pub(super) fn new(
        mass: Csr,
        stiffness: Csr,
        conductance: Csr,
        regions: Vec<Region>,
        free: Vec<usize>,
        duration: f64,
        limits: BulkTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let n = mass.nrows();
        let t = 0.5 * duration;
        if n == 0
            || !positive(duration)
            || !positive(t * t)
            || limits.maximum_factor_bytes == 0
            || !tolerance(limits.maximum_relative_residual)
            || !tolerance(limits.energy_balance_tolerance)
        {
            return Err(BulkTimeError::Invalid(
                "positive representable step, dimensions and physical limits required",
            ));
        }
        bound("pressure coordinates", n, limits.maximum_dofs)?;
        bound("memory regions", regions.len(), limits.maximum_regions)?;
        let mut entries = 0usize;
        for m in [&mass, &stiffness, &conductance] {
            if m.nrows() != n || m.ncols() != n {
                return Err(BulkTimeError::Invalid(
                    "complete original matrix shapes required",
                ));
            }
            entries = entries
                .checked_add(m.nnz())
                .ok_or(BulkTimeError::Invalid("matrix extent overflow"))?;
        }
        let mut map = vec![usize::MAX; n];
        for (i, &dof) in free.iter().enumerate() {
            if dof >= n || map[dof] != usize::MAX {
                return Err(BulkTimeError::Invalid(
                    "unique original free coordinates required",
                ));
            }
            map[dof] = i;
        }
        if free.is_empty() {
            return Err(BulkTimeError::Invalid(
                "at least one free pressure coordinate required",
            ));
        }
        let mut arms = Vec::new();
        let mut memory_dofs = 0usize;
        for (region, r) in regions.iter().enumerate() {
            poll(cancelled)?;
            if r.dofs.is_empty()
                || r.mass.nrows() != r.dofs.len()
                || r.mass.ncols() != r.dofs.len()
                || r.dofs.iter().any(|&d| d >= n)
                || r.dofs.windows(2).any(|w| w[0] >= w[1])
            {
                return Err(BulkTimeError::Invalid(
                    "original sorted region coordinates and consistent volume form required",
                ));
            }
            entries = entries
                .checked_add(r.mass.nnz())
                .ok_or(BulkTimeError::Invalid("region matrix extent overflow"))?;
            for &source in &r.arms {
                if !positive(source.delta_compressibility_pa_inv)
                    || !positive(source.relaxation_time_s)
                {
                    return Err(BulkTimeError::Invalid(
                        "positive original compressibility strength and pressure-memory time required",
                    ));
                }
                bound("memory arms", arms.len() + 1, limits.maximum_arms)?;
                memory_dofs = memory_dofs
                    .checked_add(r.dofs.len())
                    .ok_or(BulkTimeError::Invalid("memory extent overflow"))?;
                bound(
                    "memory coordinates",
                    memory_dofs,
                    limits.maximum_memory_dofs,
                )?;
                let tau = source.relaxation_time_s;
                let (retention, response) = if tau >= t {
                    let a = t / tau;
                    (1.0 / (1.0 + a), a / (1.0 + a))
                } else {
                    let a = tau / t;
                    (a / (1.0 + a), 1.0 / (1.0 + a))
                };
                arms.push(Arm {
                    region,
                    source,
                    retention,
                    response,
                });
            }
        }
        bound(
            "retained sparse entries",
            entries,
            limits.maximum_matrix_entries,
        )?;
        let nf = free.len();
        let mut coo = Coo::new(nf, nf);
        for (m, w) in [(&mass, 1.0), (&conductance, t), (&stiffness, t * t)] {
            for &row in &free {
                poll(cancelled)?;
                let (cols, values) = m.row(row);
                for (&col, &v) in cols.iter().zip(values) {
                    if map[col] != usize::MAX {
                        push(&mut coo, map[row], map[col], w * v, &limits)?;
                    }
                }
            }
        }
        // Sum arm coefficients once per shared region mass, retaining sparse structure.
        for (index, r) in regions.iter().enumerate() {
            let weight = arms
                .iter()
                .filter(|a| a.region == index)
                .fold(Dd::ZERO, |s, a| {
                    s + Dd::from_f64(a.response)
                        * Dd::from_f64(a.source.delta_compressibility_pa_inv)
                })
                .to_f64();
            for (i, &row) in r.dofs.iter().enumerate() {
                poll(cancelled)?;
                if map[row] == usize::MAX {
                    continue;
                }
                let (cols, values) = r.mass.row(i);
                for (&j, &v) in cols.iter().zip(values) {
                    let col = r.dofs[j];
                    if map[col] != usize::MAX {
                        push(&mut coo, map[row], map[col], weight * v, &limits)?;
                    }
                }
            }
        }
        let mut effective = coo.assemble();
        let mut scale = Vec::with_capacity(nf);
        for i in 0..nf {
            let d = effective.get(i, i);
            if !positive(d) {
                return Err(BulkTimeError::Invalid(
                    "positive finite physical diagonal required",
                ));
            }
            scale.push(1.0 / det::sqrt(d));
        }
        for row in 0..nf {
            poll(cancelled)?;
            let cols = effective.row(row).0.to_vec();
            for (&col, v) in cols.iter().zip(effective.row_values_mut(row)) {
                let (a, b) = if row <= col { (row, col) } else { (col, row) };
                *v = (*v * scale[a]) * scale[b];
                if !v.is_finite() {
                    return Err(BulkTimeError::Invalid("pressure congruence overflow"));
                }
            }
        }
        let gate = RefCell::new(&mut *cancelled);
        let check = || (*gate.borrow_mut())();
        let symbolic =
            SymbolicLdlt::analyze_bounded(&effective, DirectOrdering::Amd, &limits.sparse, &check)?;
        let factor = symbolic.factor_bounded(
            &effective,
            &LdltOptions {
                max_factor_bytes: Some(limits.maximum_factor_bytes),
                ..LdltOptions::default()
            },
            &limits.sparse,
            &check,
        )?;
        if factor.inertia().positive != nf || factor.inertia().negative != 0 {
            return Err(BulkTimeError::Invalid(
                "condensed original pressure operator is not positive definite",
            ));
        }
        Ok(Self {
            mass,
            stiffness,
            conductance,
            regions,
            arms,
            free,
            factor,
            scale,
            duration,
            limits,
        })
    }
    pub(super) fn duration(&self) -> f64 {
        self.duration
    }
    pub(super) fn arm_count(&self) -> usize {
        self.arms.len()
    }
    pub(super) fn arm_dofs(&self, arm: usize) -> Option<&[usize]> {
        self.arms
            .get(arm)
            .map(|a| self.regions[a.region].dofs.as_slice())
    }
    pub(super) fn zero_state(&self) -> BulkPressureState<'_> {
        BulkPressureState {
            owner: self,
            potential: vec![0.0; self.mass.nrows()],
            pressure: vec![0.0; self.mass.nrows()],
            memories: self
                .arms
                .iter()
                .map(|a| vec![0.0; self.regions[a.region].dofs.len()])
                .collect(),
        }
    }
    pub(super) fn state(
        &self,
        potential: &[f64],
        pressure: &[f64],
        memories: &[Vec<f64>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BulkPressureState<'_>> {
        let n = self.mass.nrows();
        poll(cancelled)?;
        if potential.len() != n || pressure.len() != n || memories.len() != self.arms.len() {
            return Err(BulkTimeError::Invalid(
                "all original pressure and memory coordinates required",
            ));
        }
        finite(potential)?;
        finite(pressure)?;
        let mut fixed = vec![true; n];
        for &i in &self.free {
            fixed[i] = false;
        }
        if (0..n).any(|i| fixed[i] && (potential[i] != 0.0 || pressure[i] != 0.0)) {
            return Err(BulkTimeError::Invalid(
                "exact zero pressure and potential on release trace required",
            ));
        }
        for (a, q) in self.arms.iter().zip(memories) {
            if q.len() != self.regions[a.region].dofs.len() {
                return Err(BulkTimeError::Invalid(
                    "original complete region memory order required",
                ));
            }
            finite(q)?;
        }
        let result = BulkPressureState {
            owner: self,
            potential: potential.to_vec(),
            pressure: pressure.to_vec(),
            memories: memories.to_vec(),
        };
        self.energy(&result, cancelled)?;
        Ok(result)
    }
    pub(super) fn check(&self, state: &BulkPressureState<'_>) -> Result<()> {
        if core::ptr::eq(self, state.owner) {
            Ok(())
        } else {
            Err(BulkTimeError::Invalid(
                "state belongs to another material, geometry or timestep",
            ))
        }
    }
    pub(super) fn energy(
        &self,
        state: &BulkPressureState<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BulkTimeEnergy> {
        self.check(state)?;
        let mut result = BulkTimeEnergy {
            kinetic_j: 0.5 * form(&self.stiffness, &state.potential, cancelled)?,
            frozen_compression_j: 0.5 * form(&self.mass, &state.pressure, cancelled)?,
            molecular_storage_j: 0.0,
        };
        let mut memory = Dd::ZERO;
        for (a, q) in self.arms.iter().zip(&state.memories) {
            memory = memory
                + Dd::from_f64(0.5 * a.source.delta_compressibility_pa_inv)
                    * Dd::from_f64(form(&self.regions[a.region].mass, q, cancelled)?);
        }
        result.molecular_storage_j = memory.to_f64();
        finite(&[result.total_j()])?;
        Ok(result)
    }
    pub(super) fn step<'p>(
        &'p self,
        initial: &BulkPressureState<'p>,
        load: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BulkTimeStep<'p>> {
        poll(cancelled)?;
        self.check(initial)?;
        let n = self.mass.nrows();
        if load.len() != n {
            return Err(BulkTimeError::Invalid(
                "complete original pressure load required",
            ));
        }
        finite(load)?;
        let h = self.duration;
        let t = 0.5 * h;
        let sqrt_h = det::sqrt(h);
        let mut rhs = vec![Dd::ZERO; n];
        for i in 0..n {
            poll(cancelled)?;
            rhs[i] = row(&self.mass, i, &initial.pressure)
                + Dd::from_f64(t)
                    * (Dd::from_f64(load[i]) - row(&self.stiffness, i, &initial.potential));
        }
        for (a, q) in self.arms.iter().zip(&initial.memories) {
            let r = &self.regions[a.region];
            for (i, &dof) in r.dofs.iter().enumerate() {
                poll(cancelled)?;
                rhs[dof] = rhs[dof]
                    + Dd::from_f64(a.response)
                        * Dd::from_f64(a.source.delta_compressibility_pa_inv)
                        * row(&r.mass, i, q);
            }
        }
        let mut solved: Vec<_> = self
            .free
            .iter()
            .zip(&self.scale)
            .map(|(&i, &s)| rhs[i].to_f64() * s)
            .collect();
        finite(&solved)?;
        poll(cancelled)?;
        self.factor.solve_in_place(&mut solved);
        poll(cancelled)?;
        let mut midpoint = vec![0.0; n];
        for (i, &dof) in self.free.iter().enumerate() {
            midpoint[dof] = solved[i] * self.scale[i];
        }
        finite(&midpoint)?;
        let mut next = self.zero_state();
        for i in 0..n {
            next.pressure[i] = (Dd::from_f64(2.0) * Dd::from_f64(midpoint[i])
                - Dd::from_f64(initial.pressure[i]))
            .to_f64();
            next.potential[i] = h.mul_add(midpoint[i], initial.potential[i]);
        }
        finite(&next.pressure)?;
        finite(&next.potential)?;
        for (index, a) in self.arms.iter().enumerate() {
            let r = &self.regions[a.region];
            for (i, &dof) in r.dofs.iter().enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                let q0 = Dd::from_f64(initial.memories[index][i]);
                next.memories[index][i] = (q0
                    + Dd::from_f64(2.0 * a.response) * (Dd::from_f64(midpoint[dof]) - q0))
                    .to_f64();
            }
            finite(&next.memories[index])?;
        }
        // Reconstruct every audited midpoint from the actual published endpoints.
        let pm: Vec<_> = (0..n)
            .map(|i| mid(initial.pressure[i], next.pressure[i]))
            .collect();
        let psim: Vec<_> = (0..n)
            .map(|i| mid(initial.potential[i], next.potential[i]))
            .collect();
        let mut memory_change_load = vec![Dd::ZERO; n];
        let mut accounts = Vec::with_capacity(self.arms.len());
        let mut work = BulkTimeWork {
            initial_energy: self.energy(initial, cancelled)?,
            final_energy: self.energy(&next, cancelled)?,
            ..BulkTimeWork::default()
        };
        for (index, a) in self.arms.iter().enumerate() {
            poll(cancelled)?;
            let r = &self.regions[a.region];
            let q0 = &initial.memories[index];
            let q1 = &next.memories[index];
            let mut delta = vec![0.0; r.dofs.len()];
            let mut pressure = vec![0.0; r.dofs.len()];
            let mut heat_field = vec![0.0; r.dofs.len()];
            let mut residual = 0.0_f64;
            let inverse_sqrt_tau = 1.0 / det::sqrt(a.source.relaxation_time_s);
            for (i, &dof) in r.dofs.iter().enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                let qm = mid(q0[i], q1[i]);
                delta[i] = (Dd::from_f64(q1[i]) - Dd::from_f64(q0[i])).to_f64();
                pressure[i] = pm[dof];
                let lag = (Dd::from_f64(pm[dof]) - Dd::from_f64(qm)).to_f64();
                let defect = (Dd::from_f64(a.retention) * Dd::from_f64(delta[i])
                    - Dd::from_f64(2.0 * a.response) * Dd::from_f64(lag))
                .to_f64();
                let scale = a.retention * (q0[i].abs() + q1[i].abs())
                    + 2.0 * a.response * (pm[dof].abs() + qm.abs());
                residual = residual.max(relative(defect, scale));
                heat_field[i] = product(lag, sqrt_h, inverse_sqrt_tau).to_f64();
            }
            finite(&heat_field)?;
            let strength = a.source.delta_compressibility_pa_inv;
            for (i, &dof) in r.dofs.iter().enumerate() {
                poll(cancelled)?;
                memory_change_load[dof] =
                    memory_change_load[dof] + Dd::from_f64(strength) * row(&r.mass, i, &delta);
            }
            let account = BulkMemoryWork {
                initial_energy_j: 0.5 * strength * form(&r.mass, q0, cancelled)?,
                final_energy_j: 0.5 * strength * form(&r.mass, q1, cancelled)?,
                energy_change_j: strength * change(&r.mass, q0, q1, cancelled)?,
                pressure_work_j: strength * pair(&r.mass, &pressure, &delta, cancelled)?,
                dissipation_j: strength * form(&r.mass, &heat_field, cancelled)?,
                maximum_relative_residual: residual,
            };
            residual_gate(residual, &self.limits)?;
            balance(
                "original volume memory",
                account.balance_residual_j(),
                account.energy_change_j.abs()
                    + account.pressure_work_j.abs()
                    + account.dissipation_j,
                account.initial_energy_j + account.final_energy_j,
                &self.limits,
            )?;
            work.memory_energy_change_j += account.energy_change_j;
            work.molecular_dissipation_j += account.dissipation_j;
            work.pressure_to_memory_work_j += account.pressure_work_j;
            work.maximum_relative_residual = work.maximum_relative_residual.max(residual);
            accounts.push(account);
        }
        let dp: Vec<_> = (0..n)
            .map(|i| (Dd::from_f64(next.pressure[i]) - Dd::from_f64(initial.pressure[i])).to_f64())
            .collect();
        for &i in &self.free {
            poll(cancelled)?;
            let inertia = row(&self.mass, i, &dp);
            let drag = row(&self.conductance, i, &pm);
            let spring = row(&self.stiffness, i, &psim);
            let defect = (inertia
                + memory_change_load[i]
                + Dd::from_f64(h) * (drag + spring - Dd::from_f64(load[i])))
            .to_f64();
            let scale = row_abs(&self.mass, i, &initial.pressure)
                + row_abs(&self.mass, i, &next.pressure)
                + memory_change_load[i].to_f64().abs()
                + h * (row_abs(&self.conductance, i, &pm)
                    + row_abs(&self.stiffness, i, &psim)
                    + load[i].abs());
            let kinematic = (Dd::from_f64(next.potential[i])
                - Dd::from_f64(initial.potential[i])
                - Dd::from_f64(h) * Dd::from_f64(pm[i]))
            .to_f64();
            work.maximum_relative_residual = work
                .maximum_relative_residual
                .max(relative(defect, scale))
                .max(relative(
                    kinematic,
                    next.potential[i].abs() + initial.potential[i].abs() + h * pm[i].abs(),
                ));
        }
        residual_gate(work.maximum_relative_residual, &self.limits)?;
        work.acoustic_energy_change_j =
            change(&self.mass, &initial.pressure, &next.pressure, cancelled)?
                + change(
                    &self.stiffness,
                    &initial.potential,
                    &next.potential,
                    cancelled,
                )?;
        work.boundary_dissipation_j = h * form(&self.conductance, &pm, cancelled)?;
        let mut input = Dd::ZERO;
        for (&p, &b) in pm.iter().zip(load) {
            input = input + Dd::from_f64(h) * Dd::from_f64(p) * Dd::from_f64(b);
        }
        work.input_work_j = input.to_f64();
        let scale = work.initial_energy.total_j() + work.final_energy.total_j();
        for (account, defect, physical) in [
            (
                "frozen fluid",
                work.acoustic_balance_residual_j(),
                work.acoustic_energy_change_j.abs()
                    + work.input_work_j.abs()
                    + work.boundary_dissipation_j
                    + work.pressure_to_memory_work_j.abs(),
            ),
            (
                "molecular volume",
                work.memory_balance_residual_j(),
                work.memory_energy_change_j.abs()
                    + work.molecular_dissipation_j
                    + work.pressure_to_memory_work_j.abs(),
            ),
            (
                "complete acoustic volume",
                work.balance_residual_j(),
                work.acoustic_energy_change_j.abs()
                    + work.memory_energy_change_j.abs()
                    + work.input_work_j.abs()
                    + work.boundary_dissipation_j
                    + work.molecular_dissipation_j,
            ),
        ] {
            balance(account, defect, physical, scale, &self.limits)?;
        }
        poll(cancelled)?;
        Ok(BulkTimeStep {
            state: next,
            work,
            memories: accounts,
        })
    }
}
pub(super) fn bound(what: &'static str, requested: usize, maximum: usize) -> Result<()> {
    if requested > maximum {
        Err(BulkTimeError::Limit {
            what,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}
pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(BulkTimeError::Cancelled)
    } else {
        Ok(())
    }
}
fn push(c: &mut Coo, i: usize, j: usize, v: f64, l: &BulkTimeLimits) -> Result<()> {
    if !v.is_finite() {
        return Err(BulkTimeError::Invalid(
            "physical sparse coefficient overflow",
        ));
    }
    if v != 0.0 {
        bound(
            "sparse staging entries",
            c.staged() + 1,
            l.maximum_matrix_entries,
        )?;
        c.push(i, j, v);
    }
    Ok(())
}
fn finite(v: &[f64]) -> Result<()> {
    if v.iter().any(|v| !v.is_finite()) {
        Err(BulkTimeError::Invalid("nonfinite original field or work"))
    } else {
        Ok(())
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn tolerance(v: f64) -> bool {
    positive(v) && v <= 1e-6
}
fn relative(defect: f64, scale: f64) -> f64 {
    if !defect.is_finite() || !scale.is_finite() {
        f64::INFINITY
    } else if defect == 0.0 {
        0.0
    } else if scale == 0.0 {
        f64::INFINITY
    } else {
        defect.abs() / scale
    }
}
fn residual_gate(observed: f64, limits: &BulkTimeLimits) -> Result<()> {
    if observed > limits.maximum_relative_residual {
        Err(BulkTimeError::Residual {
            observed,
            maximum: limits.maximum_relative_residual,
        })
    } else {
        Ok(())
    }
}
fn balance(
    account: &'static str,
    residual: f64,
    physical: f64,
    scale: f64,
    limits: &BulkTimeLimits,
) -> Result<()> {
    let maximum = limits.energy_balance_tolerance * physical + ROUND * (physical + scale);
    finite(&[residual, maximum])?;
    if residual.abs() > maximum {
        Err(BulkTimeError::EnergyBalance {
            account,
            residual_j: residual,
            maximum_j: maximum,
        })
    } else {
        Ok(())
    }
}
fn mid(a: f64, b: f64) -> f64 {
    ((Dd::from_f64(a) + Dd::from_f64(b)) * Dd::from_f64(0.5)).to_f64()
}
fn difference(a: f64, b: f64, c: f64) -> f64 {
    (Dd::from_f64(a) + Dd::from_f64(b) - Dd::from_f64(c)).to_f64()
}
fn product(a: f64, b: f64, c: f64) -> Dd {
    let mut f = [a, b, c];
    f.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
    Dd::from_f64(f[0]) * Dd::from_f64(f[2]) * Dd::from_f64(f[1])
}
fn row(a: &Csr, i: usize, x: &[f64]) -> Dd {
    let (c, v) = a.row(i);
    c.iter().zip(v).fold(Dd::ZERO, |s, (&j, &v)| {
        s + Dd::from_f64(v) * Dd::from_f64(x[j])
    })
}
fn row_abs(a: &Csr, i: usize, x: &[f64]) -> f64 {
    let (c, v) = a.row(i);
    c.iter().zip(v).map(|(&j, &v)| (v * x[j]).abs()).sum()
}
fn pair(a: &Csr, x: &[f64], y: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
    let mut sum = Dd::ZERO;
    for (i, &x) in x.iter().enumerate() {
        poll(cancelled)?;
        sum = sum + Dd::from_f64(x) * row(a, i, y);
    }
    let r = sum.to_f64();
    finite(&[r])?;
    Ok(r)
}
fn form(a: &Csr, x: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
    let r = pair(a, x, x, cancelled)?;
    let mut s = 0.0;
    for (i, &v) in x.iter().enumerate() {
        poll(cancelled)?;
        s += v.abs() * row_abs(a, i, x);
    }
    finite(&[s])?;
    if r < -ROUND * s {
        return Err(BulkTimeError::Invalid(
            "negative original storage or heat form",
        ));
    }
    Ok(r.max(0.0))
}
fn change(a: &Csr, x: &[f64], y: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
    let sum: Vec<_> = x.iter().zip(y).map(|(&x, &y)| mid(x, y)).collect();
    let mut value = Dd::ZERO;
    for (i, (&x, &y)) in x.iter().zip(y).enumerate() {
        poll(cancelled)?;
        value = value + (Dd::from_f64(y) - Dd::from_f64(x)) * row(a, i, &sum);
    }
    let r = value.to_f64();
    finite(&[r])?;
    Ok(r)
}

#[cfg(test)]
mod tests;
