//! Exact forced motion in declared windows of the accepted loaded body.
//!
//! The original mass and tangent are projected in full, including finite
//! eigenvector residuals and mass-orthogonality error. The existing fs-phs
//! coupled ZOH owner advances q,p with p=M_r*v and retains every glue damping
//! entry. Forces are increments about the already balanced static loads.

use super::{AssemblyVibration, checkpoint};
use fs_exec::Cx;
use fs_phs::linear_exact::{LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, PreparedExactZoh};
use fs_phs::{PortHamiltonian, QuadraticStorage, ResistanceWitness};
use fs_sparse::{Coo, DirectOrdering, LdltFactor, LdltOptions, SymbolicLdlt};

/// Explicit source-projection ceilings and the original exact owner's cap.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyResponseBudget {
    /// Conservative sum of mass/tangent projection, physical glue-source
    /// projection and rectangular resistance-witness QR multiply-add terms.
    pub maximum_projection_terms: u64,
    /// Combined projection, retained glue-source coefficient and resistance-
    /// witness buffer payload ceiling and, separately, each numeric factor's
    /// payload ceiling. Borrowed modes, container metadata and the exact owner's
    /// cache are excluded; the latter retains its original bounded preparation.
    pub maximum_projection_bytes: usize,
    /// Original exact owner budget: two scalar states per retained mode. No
    /// modes are discarded. With one independent force per mode its augmented
    /// cap admits at most 42 modes. The matching band is unused by exact ZOH.
    pub linear: LinearSplitBudget,
}

/// Physical output and independent energy accounts of one accepted interval.
#[derive(Clone, Debug)]
pub struct AssemblyResponseStep {
    /// Full component displacement/rotation, including the loaded mean shape.
    pub displacements: Vec<Vec<f64>>,
    /// Full component endpoint translation/angular velocity.
    pub velocities: Vec<Vec<f64>>,
    /// Full component displacement increment divided by the fixed duration.
    /// This is the conjugate flow for the held incremental forces/moments.
    pub mean_velocities: Vec<Vec<f64>>,
    /// Quadratic perturbation storage, excluding constant static preload energy.
    pub perturbation_energy_j: f64,
    /// Original quadratic owner's stable change of perturbation storage [J].
    pub energy_change_j: f64,
    /// Held incremental force/moment times physical displacement increment [J].
    pub supplied_work_j: f64,
    /// Independent exact-trajectory integral of v^T C_glue v [J].
    pub glue_dissipation_j: f64,
    /// Energy change plus glue heat minus supplied work [J].
    pub energy_balance_residual_j: f64,
}

/// Fixed-clock response borrowing the unchanged equilibrium and its mode basis.
///
/// This advances the complete supplied window basis, not an independently
/// certified full-body response. Omitted modes, material/radiation losses,
/// nonlinear motion and changes of tether/strut/contact branch are outside this
/// linear perturbation model. Their validity and convergence remain required
/// before using it as a physical render. No static force is applied twice.
///
/// The original bounded cold preparation and each warm pH step are indivisible
/// cancellation units. Projection/physical maps poll through the caller's Cx,
/// and no state is published before the final cancellation check.
pub struct AssemblyResponse<'v, 'body> {
    vibration: &'v AssemblyVibration<'body>,
    exact: PreparedExactZoh<'static>,
    mass: Vec<f64>,
    stiffness: Vec<f64>,
    inverse_mass: Vec<f64>,
    state: Vec<f64>,
}

impl<'v, 'body> AssemblyResponse<'v, 'body> {
    /// Prepare once; the initial state is exactly the accepted static body.
    /// Full projected M/K must be SPD. Singular/unstable pencils, non-PSD glue,
    /// inaccurate inversion, nonfinite data or exhausted budgets refuse.
    pub fn prepare(
        vibration: &'v AssemblyVibration<'body>,
        duration_s: f64,
        budget: AssemblyResponseBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, String> {
        checkpoint(cx)?;
        let n = vibration.mode_count();
        let scalar = n.checked_mul(2).ok_or("response state size overflow")?;
        let augmented = scalar
            .checked_add(n)
            .and_then(|v| v.checked_add(1))
            .and_then(|v| v.checked_mul(2))
            .ok_or("response operator size overflow")?;
        if n == 0
            || scalar > budget.linear.max_modes
            || augmented > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || !duration_s.is_finite()
            || duration_s <= 0.
            || !budget.linear.max_omega_dt.is_finite()
            || budget.linear.max_omega_dt <= 0.
        {
            return Err(
                "assembled response requires a positive clock and complete exact-owner capacity"
                    .into(),
            );
        }
        let square = n * n; // the original operator cap bounds n before allocation
        let assembly = &vibration.equilibrium.assembly;
        let model = &assembly.model;
        let free = model.free;
        let sparse = model
            .k
            .nnz()
            .checked_add(model.m.nnz())
            .ok_or("response projection size overflow")?;
        let dense_terms = (2 * square as u64)
            .checked_mul(free as u64)
            .ok_or("response projection work overflow")?;
        let terms = (n as u64)
            .checked_mul(sparse as u64)
            .and_then(|v| v.checked_add(dense_terms))
            .ok_or("response projection work overflow")?;
        let source_height = assembly.glue_damping_source_row_count();
        // The original witness keeps its source rows while normalized rows and
        // Householder storage coexist. Reserve all three rectangular buffers,
        // plus square roots/matrix copies in the common dense allowance. Its
        // QR and reconstruction work is admitted before allocating any rows.
        let witness_entries = source_height
            .checked_mul(scalar)
            .ok_or("response resistance witness size overflow")?;
        let witness_terms = (source_height as u64)
            .checked_mul(scalar as u64)
            .and_then(|v| v.checked_mul(scalar as u64))
            .and_then(|v| v.checked_mul(8))
            .and_then(|v| {
                v.checked_add(if source_height == 0 {
                    0
                } else {
                    2 * (scalar as u64).pow(3)
                })
            })
            .ok_or("response resistance witness work overflow")?;
        let reserved_terms = terms
            .checked_add(witness_terms)
            .ok_or("response projection work overflow")?;
        let bytes = free
            .checked_mul(2)
            .and_then(|v| v.checked_add(40 * square + 24 * n))
            .and_then(|v| {
                witness_entries
                    .checked_mul(3)
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| v.checked_mul(8))
            .ok_or("response projection bytes overflow")?;
        if reserved_terms > budget.maximum_projection_terms
            || bytes > budget.maximum_projection_bytes
        {
            return Err("assembled response source projection budget exceeded".into());
        }
        // Passing only the remainder makes the public junction owner's row
        // staging/source payload and work part of these same total ceilings.
        let source_bytes = budget.maximum_projection_bytes - bytes;
        let source_terms = budget.maximum_projection_terms - reserved_terms;
        if vibration.glue_damping().len() != square {
            return Err("assembled response basis and source dimensions disagree".into());
        }
        let (mass, stiffness) = vibration.projected_pencil(
            budget.maximum_projection_terms,
            budget.maximum_projection_bytes,
            cx,
        )?;
        let factor = positive_factor(&mass, n, budget.maximum_projection_bytes, cx)?;
        positive_factor(&stiffness, n, budget.maximum_projection_bytes, cx)?;
        let mut inverse_mass = vec![0.; square];
        let mut unit = vec![0.; n];
        for j in 0..n {
            checkpoint(cx)?;
            unit.fill(0.);
            unit[j] = 1.;
            let column = factor.solve(&unit);
            for i in 0..n {
                inverse_mass[i * n + j] = column[i];
            }
        }
        for i in 0..n {
            for j in 0..i {
                let value = 0.5 * inverse_mass[i * n + j] + 0.5 * inverse_mass[j * n + i];
                inverse_mass[i * n + j] = value;
                inverse_mass[j * n + i] = value;
            }
        }
        // Verify the actual symmetric inverse used in the physical generator.
        for i in 0..n {
            checkpoint(cx)?;
            for j in 0..n {
                let (mut value, mut scale) = (0., 1.);
                for k in 0..n {
                    let term = mass[i * n + k] * inverse_mass[k * n + j];
                    value += term;
                    scale += term.abs();
                }
                let expected = if i == j { 1. } else { 0. };
                if !value.is_finite()
                    || !scale.is_finite()
                    || (value - expected).abs() > 64. * n as f64 * f64::EPSILON * scale
                {
                    return Err("assembled response mass inverse residual refused".into());
                }
            }
        }
        drop((factor, unit));
        let (mut q, mut j, mut r) = (
            vec![0.; scalar * scalar],
            vec![0.; scalar * scalar],
            vec![0.; scalar * scalar],
        );
        let mut g = vec![0.; scalar * n];
        for a in 0..n {
            checkpoint(cx)?;
            j[a * scalar + n + a] = 1.;
            j[(n + a) * scalar + a] = -1.;
            g[(n + a) * n + a] = 1.;
            for b in 0..n {
                q[a * scalar + b] = stiffness[a * n + b];
                q[(n + a) * scalar + n + b] = inverse_mass[a * n + b];
                r[(n + a) * scalar + n + b] = vibration.glue_damping()[a * n + b];
            }
        }
        let witness = if source_height == 0 {
            None
        } else {
            let modes: Vec<_> = vibration.modes().map(|mode| mode.phi.as_slice()).collect();
            let rows = assembly
                .project_glue_damping_source_rows(&modes, source_bytes, source_terms, || {
                    checkpoint(cx).is_err()
                })
                .map_err(|e| e.to_string())?;
            let mut canonical_rows = vec![0.; witness_entries];
            for a in 0..source_height {
                checkpoint(cx)?;
                canonical_rows[a * scalar + n..(a + 1) * scalar]
                    .copy_from_slice(&rows[a * n..(a + 1) * n]);
            }
            drop((rows, modes));
            checkpoint(cx)?;
            Some(
                ResistanceWitness::try_from_source_rows(scalar, r.clone(), canonical_rows)
                    .map_err(|e| e.to_string())?,
            )
        };
        let storage = QuadraticStorage::new(q, scalar).map_err(|e| e.to_string())?;
        let mut system = PortHamiltonian::new(scalar, n, j, r, g, Box::new(storage))
            .map_err(|e| e.to_string())?;
        if let Some(witness) = witness {
            system = system
                .with_resistance_witness(witness)
                .map_err(|e| e.to_string())?;
        }
        checkpoint(cx)?;
        let exact = PreparedExactZoh::new_owned(system, duration_s, budget.linear)
            .map_err(|e| e.to_string())?;
        checkpoint(cx)?;
        Ok(Self {
            vibration,
            exact,
            mass,
            stiffness,
            inverse_mass,
            state: vec![0.; scalar],
        })
    }

    /// Fixed held-input interval [s]. Preparation is reused on every step.
    pub fn duration_s(&self) -> f64 {
        self.exact.duration()
    }

    /// Original projected mass Phi^T M Phi, row-major; never replaced by I.
    pub fn mass_matrix(&self) -> &[f64] {
        &self.mass
    }

    /// Original projected tangent Phi^T K_t Phi, including its off-diagonals.
    pub fn stiffness_matrix(&self) -> &[f64] {
        &self.stiffness
    }

    /// Complete projected glue damping, row-major, counted exactly once.
    pub fn damping_matrix(&self) -> &[f64] {
        self.vibration.glue_damping()
    }

    /// Canonical [q, p] state, p=M_r*v. Copy this for an exact checkpoint.
    /// Restoring requires this same prepared body, mode order and clock.
    pub fn canonical_state(&self) -> &[f64] {
        &self.state
    }

    /// Restore complete canonical state atomically, without changing the model.
    pub fn restore_state(&mut self, state: &[f64], cx: &Cx<'_>) -> Result<(), String> {
        checkpoint(cx)?;
        if state.len() != self.state.len()
            || state.iter().any(|v| !v.is_finite())
            || !self.exact.system().hamiltonian(state).is_finite()
        {
            return Err("assembled response checkpoint must be complete and finite".into());
        }
        let state = state.to_vec();
        checkpoint(cx)?;
        self.state = state;
        Ok(())
    }

    /// Set modal displacement and velocity increments about loaded equilibrium.
    pub fn initialize(&mut self, q: &[f64], v: &[f64], cx: &Cx<'_>) -> Result<(), String> {
        checkpoint(cx)?;
        let n = self.vibration.mode_count();
        if q.len() != n || v.len() != n || q.iter().chain(v).any(|x| !x.is_finite()) {
            return Err(
                "assembled response initialization must match the finite modal basis".into(),
            );
        }
        let momentum = apply_dense(&self.mass, v, cx)?;
        let mut state = q.to_vec();
        state.extend(momentum);
        self.restore_state(&state, cx)
    }

    /// Advance one held incremental full-component force/moment vector.
    /// The static preload is already balanced; zero input and zero state stay
    /// exactly at that equilibrium. Every failure leaves the physical state
    /// unchanged. Warm work uses the original cached transition/heat owner.
    pub fn step(&mut self, loads: &[&[f64]], cx: &Cx<'_>) -> Result<AssemblyResponseStep, String> {
        self.step_with_checkpoint(loads, cx, &mut || checkpoint(cx))
    }

    fn step_with_checkpoint(
        &mut self,
        loads: &[&[f64]],
        cx: &Cx<'_>,
        check: &mut impl FnMut() -> Result<(), String>,
    ) -> Result<AssemblyResponseStep, String> {
        check()?;
        let n = self.vibration.mode_count();
        let force = self.vibration.project_loads(loads, cx)?;
        let mut next = vec![0.; 2 * n];
        let mut mean = vec![0.; n];
        let record = self
            .exact
            .step_into(&self.state, &force, &mut next, &mut mean)
            .map_err(|e| e.to_string())?;
        check()?;
        let velocity = apply_dense(&self.inverse_mass, &next[n..], cx)?;
        let delta: Vec<_> = next[..n]
            .iter()
            .zip(&self.state[..n])
            .map(|(a, b)| a - b)
            .collect();
        let work: f64 = force.iter().zip(&delta).map(|(a, b)| a * b).sum();
        for (value, increment) in mean.iter_mut().zip(&delta) {
            *value = increment / self.duration_s();
        }
        let energy = self.exact.system().hamiltonian(&next);
        let scale = energy
            .abs()
            .max(self.exact.system().hamiltonian(&self.state).abs())
            .max(record.supplied.abs())
            .max(work.abs())
            .max(record.dissipated.abs());
        let residual = record.delta_h + record.dissipated - work;
        if !work.is_finite()
            || !energy.is_finite()
            || !residual.is_finite()
            || record.dissipated < 0.
            || (work - record.supplied).abs() > 1e-12 * scale
            || residual.abs() > 1e-12 * scale
        {
            return Err("assembled response held-force work/energy balance refused".into());
        }
        let displacements = self.vibration.displacements(&next[..n], cx)?;
        let velocities = self.vibration.perturbation_displacements(&velocity, cx)?;
        let mean_velocities = self.vibration.perturbation_displacements(&mean, cx)?;
        if displacements
            .iter()
            .chain(&velocities)
            .chain(&mean_velocities)
            .flatten()
            .any(|v| !v.is_finite())
        {
            return Err("assembled response physical output overflow".into());
        }
        check()?;
        self.state = next;
        Ok(AssemblyResponseStep {
            displacements,
            velocities,
            mean_velocities,
            perturbation_energy_j: energy,
            energy_change_j: record.delta_h,
            supplied_work_j: work,
            glue_dissipation_j: record.dissipated,
            energy_balance_residual_j: residual,
        })
    }
}

fn positive_factor(
    matrix: &[f64],
    n: usize,
    bytes: usize,
    cx: &Cx<'_>,
) -> Result<LdltFactor, String> {
    checkpoint(cx)?;
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if matrix[i * n + j] != 0. {
                coo.push(i, j, matrix[i * n + j]);
            }
        }
    }
    let csr = coo.assemble();
    let factor = SymbolicLdlt::analyze(&csr, DirectOrdering::Amd)
        .and_then(|s| {
            s.factor(
                &csr,
                &LdltOptions {
                    max_factor_bytes: Some(bytes),
                    ..LdltOptions::default()
                },
            )
        })
        .map_err(|e| e.to_string())?;
    if factor.inertia().positive != n || factor.inertia().negative != 0 {
        return Err(
            "assembled response requires positive definite projected mass and tangent".into(),
        );
    }
    checkpoint(cx)?;
    Ok(factor)
}

fn apply_dense(matrix: &[f64], vector: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, String> {
    let n = vector.len();
    let mut out = vec![0.; n];
    for i in 0..n {
        checkpoint(cx)?;
        out[i] = matrix[i * n..(i + 1) * n]
            .iter()
            .zip(vector)
            .map(|(a, b)| a * b)
            .sum();
        if !out[i].is_finite() {
            return Err("assembled response state transformation overflow".into());
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
