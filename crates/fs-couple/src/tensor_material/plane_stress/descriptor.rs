//! Rank-free plane stress through the existing port-Hamiltonian descriptor.
//!
//! Keep all normal strains e_n as stored coordinates, and introduce algebraic
//! efforts lambda = e_n'. Their lossless pair gives e_n'=lambda and
//! 0=-H_en-D_nu*v-D_nn*lambda. The ORIGINAL full Kelvin matrix is a PSD block
//! over [v,lambda]. It is never inverted, ranked, clipped or regularized.
//! Maxwell spring coordinates z retain every original factor row and satisfy
//! z'=A_u*v+A_n*lambda-z/tau. The shared descriptor owns the complete step.
use super::*;
use fs_phs::DescriptorPortHamiltonian;

const TRACTION_RELATIVE_TOLERANCE: f64 = 1e-9;
const ENERGY_RELATIVE_TOLERANCE: f64 = 1e-10;

/// Original 3D material under local z-plane stress, with all three normal
/// strains retained. This admits oblique singular Kelvin viscosity without
/// guessing its rank. Columns are [epsilon_xx,epsilon_yy,gamma_xy] per metre.
#[derive(Clone, Debug)]
pub struct DescriptorPlaneStressPoint {
    source: ProjectedTensorProny,
}

impl DescriptorPlaneStressPoint {
    /// Project the original law once. The coordinate budget includes three
    /// local normal strains; every six-row Maxwell factor remains complete.
    pub fn new(
        material: &TensorProny,
        columns: &[[f64; 3]],
        volume_m3: f64,
        budget: TensorProjectionBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        if columns.is_empty()
            || columns
                .len()
                .checked_add(3)
                .is_none_or(|n| n > budget.max_coordinates)
        {
            return Err(bad(
                "descriptor plane stress needs bounded nonempty strain columns",
            ));
        }
        let mut full: Vec<Voigt> = columns
            .iter()
            .map(|s| [s[0], s[1], 0.0, 0.5 * s[2], 0.0, 0.0])
            .collect();
        full.extend([
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        let source = material
            .project_strain(&full, volume_m3, budget, || gate.is_requested())
            .map_err(TensorMaterialError::Material)?;
        Ok(Self { source })
    }
    /// Original physical displacement-coordinate count.
    pub fn coordinates(&self) -> usize {
        self.source.coordinates() - 3
    }
    /// Six spring-memory coordinates per original Maxwell branch.
    pub fn memory_coordinates(&self) -> usize {
        6 * self.source.branches().len()
    }
    /// Original immutable constitutive law and applicability domain.
    pub fn material(&self) -> &TensorProny {
        self.source.material()
    }
    /// Original full strain/factor projection; no reduced material is fitted.
    pub fn projection(&self) -> &ProjectedTensorProny {
        &self.source
    }
}

/// Complete physical inertia and rank-free plane-stress constraints.
///
/// Unlike the condensed exact-ZOH body, this explicit API uses the existing
/// second-order descriptor/Gonzalez step. Energy and traction acceptance DO NOT
/// certify frequency accuracy. Callers must refine time; this owner is not
/// automatically admitted to an audio render clock. It adds no normal inertia.
pub struct DescriptorPlaneStressBody {
    system: DescriptorPortHamiltonian,
    inertia: Inertia,
    points: Vec<DescriptorPlaneStressPoint>,
    offsets: Vec<usize>,
}

fn extent(
    n: usize,
    points: &[DescriptorPlaneStressPoint],
    budget: TensorMaterialBudget,
) -> Result<(usize, usize)> {
    if n == 0
        || n > 128
        || points.is_empty()
        || points.len() > budget.max_material_points
        || budget.max_material_points == 0
        || budget.max_material_points > 4096
        || budget.max_state_coordinates == 0
        || budget.max_state_coordinates > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        || points.iter().any(|p| p.coordinates() != n)
    {
        return Err(bad(
            "bounded matching descriptor plane-stress points and inertia required",
        ));
    }
    let differential = points
        .iter()
        .try_fold(2 * n, |s, p| {
            s.checked_add(p.memory_coordinates())
                .and_then(|s| s.checked_add(3))
        })
        .ok_or_else(|| bad("descriptor material differential extent overflow"))?;
    let total = points
        .len()
        .checked_mul(3)
        .and_then(|k| differential.checked_add(k))
        .ok_or_else(|| bad("descriptor material algebraic extent overflow"))?;
    if total > budget.max_state_coordinates {
        return Err(bad(
            "complete material, normal-strain and rate coordinates exceed capacity",
        ));
    }
    Ok((differential, total))
}

impl DescriptorPlaneStressBody {
    pub fn new(
        masses_kg: Vec<f64>,
        points: Vec<DescriptorPlaneStressPoint>,
        budget: TensorMaterialBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        let dimensions = extent(masses_kg.len(), &points, budget)?;
        Self::assemble(Inertia::diagonal(masses_kg)?, points, dimensions, gate)
    }
    /// Keep the entire SPD mass matrix and p=M*v in the original coordinates.
    pub fn from_mass_matrix(
        coordinates: usize,
        mass_kg: Vec<f64>,
        points: Vec<DescriptorPlaneStressPoint>,
        budget: TensorMaterialBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        let dimensions = extent(coordinates, &points, budget)?;
        Self::assemble(
            Inertia::consistent(coordinates, mass_kg, gate)?,
            points,
            dimensions,
            gate,
        )
    }

    fn assemble(
        inertia: Inertia,
        points: Vec<DescriptorPlaneStressPoint>,
        (differential, total): (usize, usize),
        gate: &CancelGate,
    ) -> Result<Self> {
        let n = inertia.dimension();
        let mut q = vec![0.0; differential * differential];
        let mut j = vec![0.0; total * total];
        let mut r = vec![0.0; total * total];
        let mut g = vec![0.0; total * n];
        for i in 0..n {
            for k in 0..n {
                q[(n + i) * differential + n + k] = inertia.inverse_entry(i, k);
            }
            j[i * total + n + i] = 1.0;
            j[(n + i) * total + i] = -1.0;
            g[(n + i) * n + i] = 1.0;
        }
        let mut offsets = Vec::with_capacity(points.len());
        let mut offset = 2 * n;
        for (index, point) in points.iter().enumerate() {
            poll(gate)?;
            offsets.push(offset);
            let memory = point.memory_coordinates();
            let normal = offset + memory;
            let rate = differential + 3 * index;
            let d = n + 3;
            let strain_index = |i| if i < n { i } else { normal + i - n };
            let rate_index = |i| if i < n { n + i } else { rate + i - n };
            for a in 0..d {
                for b in 0..d {
                    q[strain_index(a) * differential + strain_index(b)] +=
                        point.source.equilibrium_stiffness()[a * d + b];
                    r[rate_index(a) * total + rate_index(b)] +=
                        point.source.viscosity_matrix()[a * d + b];
                }
            }
            for k in 0..3 {
                j[(normal + k) * total + rate + k] = 1.0;
                j[(rate + k) * total + normal + k] = -1.0;
            }
            for (branch_index, branch) in point.source.branches().iter().enumerate() {
                let relaxation = branch.relaxation_time_s().recip();
                if !relaxation.is_finite() || relaxation <= 0.0 {
                    return Err(bad("unrepresentable Maxwell rate"));
                }
                for (row_index, row) in branch.rows().iter().enumerate() {
                    let z = offset + 6 * branch_index + row_index;
                    q[z * differential + z] = 1.0;
                    r[z * total + z] = relaxation;
                    for k in 0..d {
                        j[z * total + rate_index(k)] = row[k];
                        j[rate_index(k) * total + z] = -row[k];
                    }
                }
            }
            offset = normal + 3;
        }
        poll(gate)?;
        let storage =
            QuadraticStorage::new(q, differential).map_err(TensorMaterialError::Propagation)?;
        let system =
            DescriptorPortHamiltonian::new(total, differential, n, j, r, g, Box::new(storage))
                .map_err(TensorMaterialError::Propagation)?;
        poll(gate)?;
        Ok(Self {
            system,
            inertia,
            points,
            offsets,
        })
    }
    pub fn coordinates(&self) -> usize {
        self.inertia.dimension()
    }
    pub fn state_dimension(&self) -> usize {
        self.system.state_dim()
    }
    pub fn differential_dimension(&self) -> usize {
        self.system.differential_dim()
    }
    pub fn mass_matrix_kg(&self) -> Vec<f64> {
        self.inertia.matrix()
    }
    pub fn material_points(&self) -> &[DescriptorPlaneStressPoint] {
        &self.points
    }
    /// Existing descriptor operator with the original force/velocity ports.
    pub fn port_system(&self) -> &DescriptorPortHamiltonian {
        &self.system
    }

    fn point_coordinates(&self, x: &[f64], point: usize) -> Vec<f64> {
        let n = self.coordinates();
        let start = self.offsets[point] + self.points[point].memory_coordinates();
        let mut full = x[..n].to_vec();
        full.extend_from_slice(&x[start..start + 3]);
        full
    }
    fn traction(&self, x: &[f64], point: usize) -> ([f64; 3], [f64; 3]) {
        let p = &self.points[point];
        let n = self.coordinates();
        let d = n + 3;
        let full = self.point_coordinates(x, point);
        let velocity = self.inertia.velocity(&x[n..2 * n]);
        let rate_start = self.differential_dimension() + 3 * point;
        let mut rates = velocity;
        rates.extend_from_slice(&x[rate_start..rate_start + 3]);
        let mut force = [0.0; 3];
        let mut magnitude = [0.0; 3];
        for i in 0..3 {
            for k in 0..d {
                for value in [
                    p.source.equilibrium_stiffness()[(n + i) * d + k] * full[k],
                    p.source.viscosity_matrix()[(n + i) * d + k] * rates[k],
                ] {
                    force[i] += value;
                    magnitude[i] += value.abs();
                }
            }
            for (k, row) in p
                .source
                .branches()
                .iter()
                .flat_map(|b| b.rows())
                .enumerate()
            {
                let value = row[n + i] * x[self.offsets[point] + k];
                force[i] += value;
                magnitude[i] += value.abs();
            }
        }
        (force, magnitude)
    }
    fn admit_state(&self, x: &[f64]) -> Result<()> {
        if x.len() != self.state_dimension() || !finite(x) {
            return Err(bad("complete finite descriptor material state required"));
        }
        for (index, point) in self.points.iter().enumerate() {
            point
                .source
                .strain(&self.point_coordinates(x, index))
                .map_err(TensorMaterialError::Material)?;
            let (traction, scale) = self.traction(x, index);
            if (0..3).any(|i| {
                !traction[i].is_finite()
                    || !scale[i].is_finite()
                    || traction[i].abs()
                        > TRACTION_RELATIVE_TOLERANCE * scale[i].max(f64::MIN_POSITIVE)
            }) {
                return Err(bad(
                    "original elastic, Maxwell and Kelvin normal traction is inconsistent",
                ));
            }
        }
        if !self.system.hamiltonian(x).is_finite() {
            return Err(bad("descriptor material energy is nonfinite"));
        }
        Ok(())
    }

    /// Explicit relaxed prehistory at rest. Equilibrium normal strains solve
    /// the ORIGINAL SPD elastic block, all Maxwell spring strains and all rate
    /// efforts start at zero. Moving preloads require a complete checked restart.
    pub fn history_relaxed_at_rest(
        &self,
        displacement_m: &[f64],
        gate: &CancelGate,
    ) -> Result<DescriptorPlaneStressHistory<'_>> {
        poll(gate)?;
        let n = self.coordinates();
        if displacement_m.len() != n || !finite(displacement_m) {
            return Err(bad("finite matching initial displacement required"));
        }
        let mut state = vec![0.0; self.state_dimension()];
        state[..n].copy_from_slice(displacement_m);
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            let d = n + 3;
            let stiffness = point.source.equilibrium_stiffness();
            let rhs: [f64; 3] = core::array::from_fn(|a| {
                -(0..n)
                    .map(|k| stiffness[(n + a) * d + k] * displacement_m[k])
                    .sum::<f64>()
            });
            let normal = normal_solve(stiffness, d, n, &rhs)?;
            let at = self.offsets[i] + point.memory_coordinates();
            state[at..at + 3].copy_from_slice(&normal);
        }
        self.history_from_state(state, gate)
    }

    /// Restart the SAME operator from a complete accepted state. State order:
    /// [u,p], then each point's [six Maxwell spring coordinates per branch,
    /// three normal tensor strains], then three algebraic rate representatives
    /// per point. No missing history is inferred. Original normal traction and
    /// strain-domain admission must hold; counters begin at this explicit state.
    /// Algebraic representatives have no energy and are not an interpolant of
    /// instantaneous rates in undamped constraint directions.
    pub fn history_from_state(
        &self,
        state: Vec<f64>,
        gate: &CancelGate,
    ) -> Result<DescriptorPlaneStressHistory<'_>> {
        poll(gate)?;
        self.admit_state(&state)?;
        let initial_energy = self.system.hamiltonian(&state);
        poll(gate)?;
        Ok(DescriptorPlaneStressHistory {
            body: self,
            state,
            initial_energy,
            work: 0.0,
            heat: 0.0,
            time: 0.0,
            steps: 0,
        })
    }
}

/// Accepted complete state and separate physical work/heat accounts.
pub struct DescriptorPlaneStressHistory<'a> {
    body: &'a DescriptorPlaneStressBody,
    state: Vec<f64>,
    initial_energy: f64,
    work: f64,
    heat: f64,
    time: f64,
    steps: u64,
}
impl DescriptorPlaneStressHistory<'_> {
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    pub fn displacement_m(&self) -> &[f64] {
        &self.state[..self.body.coordinates()]
    }
    pub fn velocity_m_s(&self) -> Vec<f64> {
        let n = self.body.coordinates();
        self.body.inertia.velocity(&self.state[n..2 * n])
    }
    pub fn time_s(&self) -> f64 {
        self.time
    }
    pub fn steps(&self) -> u64 {
        self.steps
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.body.system.hamiltonian(&self.state)
    }
    pub fn supplied_work_j(&self) -> f64 {
        self.work
    }
    pub fn dissipated_heat_j(&self) -> f64 {
        self.heat
    }
    pub fn balance_residual_j(&self) -> f64 {
        self.stored_energy_j() - self.initial_energy + self.heat - self.work
    }
    pub fn plane_stress_strain(&self, point: usize) -> Result<Voigt> {
        let p = self
            .body
            .points
            .get(point)
            .ok_or_else(|| bad("unknown descriptor material point"))?;
        p.source
            .strain(&self.body.point_coordinates(&self.state, point))
            .map_err(TensorMaterialError::Material)
    }
    pub fn plane_stress_normal_stress_pa(&self, point: usize) -> Result<[f64; 3]> {
        let p = self
            .body
            .points
            .get(point)
            .ok_or_else(|| bad("unknown descriptor material point"))?;
        let (f, _) = self.body.traction(&self.state, point);
        let volume = p.source.volume_m3();
        Ok([f[0] / volume, 0.5 * f[1] / volume, 0.5 * f[2] / volume])
    }
    pub fn advance(&mut self, force_n: &[f64], dt_s: f64, gate: &CancelGate) -> Result<StepRecord> {
        self.advance_controlled(force_n, dt_s, || gate.is_requested())
    }
    /// Existing bounded descriptor solve; cancellation is polled around that
    /// solve and immediately before atomic publication, not inside its Newton
    /// work. No frequency/phase certificate is inferred from energy balance.
    pub fn advance_controlled(
        &mut self,
        force_n: &[f64],
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepRecord> {
        if cancelled() {
            return Err(TensorMaterialError::Cancelled);
        }
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || force_n.len() != self.body.coordinates()
            || !finite(force_n)
        {
            return Err(bad(
                "positive finite clock and complete finite force required",
            ));
        }
        let time = self.time + dt_s;
        let steps = self
            .steps
            .checked_add(1)
            .ok_or_else(|| bad("descriptor sample count overflow"))?;
        if !time.is_finite() || time <= self.time {
            return Err(bad("descriptor clock is unrepresentable"));
        }
        let record = fs_phs::step_descriptor(&self.body.system, &self.state, force_n, dt_s)
            .map_err(TensorMaterialError::Propagation)?;
        if cancelled() {
            return Err(TensorMaterialError::Cancelled);
        }
        self.body.admit_state(&record.x)?;
        let before = self.stored_energy_j();
        let after = self.body.system.hamiltonian(&record.x);
        let scale = before
            .abs()
            .max(after.abs())
            .max(record.supplied.abs())
            .max(record.dissipated.abs())
            .max(f64::MIN_POSITIVE);
        let balance = after - before + record.dissipated - record.supplied;
        if !balance.is_finite()
            || record.dissipated < 0.0
            || balance.abs() > ENERGY_RELATIVE_TOLERANCE * scale
        {
            return Err(bad(
                "descriptor material physical heat or original energy balance refused",
            ));
        }
        let work = self.work + record.supplied;
        let heat = self.heat + record.dissipated;
        if !work.is_finite() || !heat.is_finite() {
            return Err(bad("descriptor material counters overflow"));
        }
        if cancelled() {
            return Err(TensorMaterialError::Cancelled);
        }
        self.state = record.x.clone();
        self.time = time;
        self.steps = steps;
        self.work = work;
        self.heat = heat;
        Ok(record)
    }
}

#[cfg(test)]
mod tests;
