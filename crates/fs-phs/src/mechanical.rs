//! Original-coordinate linear mechanics through the existing pHS owner.
//!
//! `M q_ddot + C q_dot + K_relaxed q + sum(b_j s_j) = P u`,
//! `s_j_dot = k_j b_j^T q_dot - s_j/tau_j`. Consistent mass, full
//! nonproportional resistance, gyroscopic coupling and material histories
//! remain. No modal diagonalization, mass lumping or new integrator is used.

use crate::interconnection::accumulation::Sum;
use crate::linear_exact::positive_heat::factor_resistance;
use crate::{
    PhsError, PortHamiltonian, QuadraticStorage, RelaxationBranch, ResistanceWitness,
    require_symmetric,
};
use fs_la::factor::cholesky;
use fs_math::det;

mod impulse;
pub use impulse::{MechanicalImpulse, MechanicalImpulseControls};

/// Maximum complete state dimension of this bounded dense source adapter.
/// Exact-propagation owners impose their own additional preparation budgets.
pub const MAX_MECHANICAL_STATES: usize = crate::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM;

/// Original structural matrices in canonical physical coordinates.
///
/// State order is `[q_0..q_n, p_0..p_n, z_0..z_b]`, with `p=M v` and
/// `z_j=sqrt(k_j)*viscous_strain_j`. Each branch stores
/// `(sqrt(k_j)*b_j^T q-z_j)^2/2`. K is the RELAXED stiffness; branch stiffness
/// is added by the existing Maxwell owner exactly once. Input forces and
/// conjugate output velocities satisfy `y=P^T v` and power `u^T y`.
///
/// `system()` supports the existing discrete-gradient, prepared exact ZOH,
/// continuous-tone, harmonic and power-port interconnection APIs. Their
/// state, work, heat and admission rules apply unchanged. The mass inverse
/// and physical resistance factor are prepared once; warm stepping adds no
/// mass solve or factorization. This adapter performs no time integration.
pub struct LinearMechanicalSystem {
    system: PortHamiltonian,
    coordinates: usize,
    mass: Vec<f64>,
    inverse_mass: Vec<f64>,
    branches: Vec<RelaxationBranch>,
}

impl LinearMechanicalSystem {
    /// Admit row-major n-by-n M, C and relaxed K, and n-by-ports P.
    ///
    /// M must be symmetric positive definite, K symmetric PSD, and the
    /// symmetric part of C PSD. The skew part of C remains gyroscopic and
    /// produces no heat. Zero stiffness/damping and zero strain directions
    /// retain free or independently relaxing states. Each branch projection
    /// has n entries in ORIGINAL displacement order; equal-time arms remain
    /// distinct. The retained resistance witness includes every loss channel.
    ///
    /// State and port counts are bounded by `MAX_MECHANICAL_STATES`. Native
    /// floating-point structural/factor admission applies; this is not an
    /// interval conditioning certificate or a material-validity claim.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mass: &[f64],
        damping: &[f64],
        relaxed_stiffness: &[f64],
        force_map: &[f64],
        coordinates: usize,
        ports: usize,
        branches: Vec<RelaxationBranch>,
    ) -> Result<Self, PhsError> {
        let n = coordinates;
        let state_dim = n
            .checked_mul(2)
            .and_then(|v| v.checked_add(branches.len()))
            .ok_or(PhsError::LinearSplitBudget {
                what: "mechanical state extent overflow",
            })?;
        if n == 0 || state_dim > MAX_MECHANICAL_STATES || ports > MAX_MECHANICAL_STATES {
            return Err(PhsError::LinearSplitBudget {
                what: "mechanical state or port capacity",
            });
        }
        // Q and the appended memory use native dense spectral admission.
        // Preflight its work ceiling before its infallible numerical entry.
        fs_la::eigen::admit_jacobi_eigh(state_dim).map_err(|_| PhsError::LinearSplitBudget {
            what: "mechanical dense admission capacity",
        })?;
        let square = n * n;
        if mass.len() != square
            || damping.len() != square
            || relaxed_stiffness.len() != square
            || force_map.len() != n * ports
            || branches.iter().any(|b| b.projection.len() != n)
        {
            return Err(PhsError::Dimension {
                what: "original mechanical matrices and strain maps",
            });
        }
        if mass
            .iter()
            .chain(damping)
            .chain(relaxed_stiffness)
            .chain(force_map)
            .any(|v| !v.is_finite())
        {
            return Err(PhsError::Unadmitted {
                what: "finite mechanical source matrices",
            });
        }
        require_symmetric(mass, n, "mechanical mass")?;
        require_symmetric(relaxed_stiffness, n, "mechanical relaxed stiffness")?;
        // Resolve each original elastic direction before unrelated kinetic or
        // branch storage can enlarge the aggregate quadratic PSD allowance.
        factor_resistance(relaxed_stiffness, n).map_err(|_| PhsError::Unadmitted {
            what: "resolved positive semidefinite mechanical relaxed stiffness",
        })?;
        let factor = cholesky(mass, n).map_err(|_| PhsError::Unadmitted {
            what: "positive definite mechanical mass",
        })?;
        let mut inverse = vec![0.0; square];
        let mut column = vec![0.0; n];
        for j in 0..n {
            column.fill(0.0);
            column[j] = 1.0;
            factor.solve(&mut column);
            for i in 0..n {
                inverse[i * n + j] = column[i];
            }
        }
        // Triangular solves can round reciprocal inverse entries differently.
        for i in 0..n {
            for j in 0..i {
                let value = 0.5 * inverse[i * n + j] + 0.5 * inverse[j * n + i];
                inverse[i * n + j] = value;
                inverse[j * n + i] = value;
            }
        }
        if inverse.iter().any(|v| !v.is_finite()) {
            return Err(PhsError::Unadmitted {
                what: "representable mechanical inverse mass",
            });
        }
        let symmetric: Vec<_> = (0..square)
            .map(|index| {
                let (a, b) = (index / n, index % n);
                0.5 * damping[a * n + b] + 0.5 * damping[b * n + a]
            })
            .collect();
        // Reuse native positive-factor admission: no clipping or fitting C.
        let damping_root = factor_resistance(&symmetric, n)?;
        let d = 2 * n;
        let mut q = vec![0.0; d * d];
        let mut j = vec![0.0; d * d];
        let mut resistance = vec![0.0; d * d];
        let mut root = vec![0.0; d * d];
        let mut g = vec![0.0; d * ports];
        for a in 0..n {
            j[a * d + n + a] = 1.0;
            j[(n + a) * d + a] = -1.0;
            for b in 0..n {
                q[a * d + b] = relaxed_stiffness[a * n + b];
                q[(n + a) * d + n + b] = inverse[a * n + b];
                resistance[(n + a) * d + n + b] = symmetric[a * n + b];
                root[(n + a) * d + n + b] = damping_root[a * n + b];
                j[(n + a) * d + n + b] = -(0.5 * damping[a * n + b] - 0.5 * damping[b * n + a]);
            }
            for p in 0..ports {
                g[(n + a) * ports + p] = force_map[a * ports + p];
            }
        }
        let witness = ResistanceWitness {
            matrix: resistance.clone(),
            root,
        };
        let storage = Box::new(QuadraticStorage::new(q, d)?);
        let base = PortHamiltonian::new(d, ports, j, resistance, g, storage)?
            .with_resistance_witness(witness)?;
        let canonical_branches = branches
            .iter()
            .map(|branch| {
                let mut projection = branch.projection.clone();
                projection.resize(d, 0.0);
                RelaxationBranch {
                    projection,
                    stiffness: branch.stiffness,
                    relaxation_time_s: branch.relaxation_time_s,
                }
            })
            .collect();
        let system = base.with_relaxation_branches(canonical_branches)?;
        Ok(Self {
            system,
            coordinates: n,
            mass: mass.to_vec(),
            inverse_mass: inverse,
            branches,
        })
    }

    pub const fn coordinates(&self) -> usize {
        self.coordinates
    }
    pub fn branches(&self) -> &[RelaxationBranch] {
        &self.branches
    }
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }

    /// Initialize original displacement, velocity and per-arm viscous strain.
    /// Zero viscous strain retains a loaded arm when b^T q is nonzero; use
    /// `relaxed_state` for explicitly fully relaxed initial material.
    pub fn state(
        &self,
        displacement: &[f64],
        velocity: &[f64],
        viscous_strain: &[f64],
    ) -> Result<Vec<f64>, PhsError> {
        if viscous_strain.len() != self.branches.len() || !finite(viscous_strain) {
            return Err(PhsError::Dimension {
                what: "mechanical viscous strain state",
            });
        }
        let mut state = self.prefix(displacement, velocity)?;
        for (i, branch) in self.branches.iter().enumerate() {
            state[2 * self.coordinates + i] = det::sqrt(branch.stiffness) * viscous_strain[i];
        }
        self.admit_state(&state)?;
        Ok(state)
    }

    /// Explicit zero-stress branch initialization on the original displacement.
    pub fn relaxed_state(
        &self,
        displacement: &[f64],
        velocity: &[f64],
    ) -> Result<Vec<f64>, PhsError> {
        let mut state = self.prefix(displacement, velocity)?;
        for (i, branch) in self.branches.iter().enumerate() {
            let root = det::sqrt(branch.stiffness);
            state[2 * self.coordinates + i] = branch
                .projection
                .iter()
                .zip(displacement)
                .map(|(b, q)| (root * b) * q)
                .sum();
        }
        self.admit_state(&state)?;
        Ok(state)
    }

    /// Physical velocities on the original coordinates, with all memory retained.
    pub fn velocities(&self, state: &[f64]) -> Result<Vec<f64>, PhsError> {
        let mut velocity = vec![0.0; self.coordinates];
        self.velocities_into(state, &mut velocity)?;
        Ok(velocity)
    }

    /// Allocation-free full-inertia observation; no elastic or memory effort
    /// evaluation is needed. On refusal, output is scratch; state is unchanged.
    pub fn velocities_into(&self, state: &[f64], velocity: &mut [f64]) -> Result<(), PhsError> {
        self.admit_state(state)?;
        let n = self.coordinates;
        if velocity.len() != n {
            return Err(PhsError::Dimension {
                what: "mechanical velocity output extent",
            });
        }
        for (i, value) in velocity.iter_mut().enumerate() {
            *value = Sum::dot(&self.inverse_mass[i * n..(i + 1) * n], &state[n..2 * n])
                .and_then(Sum::value)
                .ok_or(PhsError::Unadmitted {
                    what: "representable original mechanical velocity",
                })?;
        }
        Ok(())
    }

    /// Physical stresses k_j*(b_j^T q-viscous_strain_j), in original arm order.
    pub fn branch_stresses(&self, state: &[f64]) -> Result<Vec<f64>, PhsError> {
        self.admit_state(state)?;
        let effort = self.system.effort(state);
        let stresses: Vec<_> = self
            .branches
            .iter()
            .enumerate()
            .map(|(i, branch)| -det::sqrt(branch.stiffness) * effort[2 * self.coordinates + i])
            .collect();
        if !finite(&stresses) {
            return Err(PhsError::Unadmitted {
                what: "finite mechanical stress",
            });
        }
        Ok(stresses)
    }

    fn prefix(&self, displacement: &[f64], velocity: &[f64]) -> Result<Vec<f64>, PhsError> {
        let n = self.coordinates;
        if displacement.len() != n
            || velocity.len() != n
            || !finite(displacement)
            || !finite(velocity)
        {
            return Err(PhsError::Dimension {
                what: "mechanical displacement and velocity",
            });
        }
        let mut state = vec![0.0; self.system.state_dim()];
        state[..n].copy_from_slice(displacement);
        for i in 0..n {
            state[n + i] = Sum::dot(&self.mass[i * n..(i + 1) * n], velocity)
                .and_then(Sum::value)
                .ok_or(PhsError::Unadmitted {
                    what: "representable original mechanical momentum",
                })?;
        }
        Ok(state)
    }
    fn admit_state(&self, state: &[f64]) -> Result<(), PhsError> {
        if state.len() != self.system.state_dim() || !finite(state) {
            return Err(PhsError::Dimension {
                what: "complete finite mechanical state",
            });
        }
        Ok(())
    }
}
fn finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}
