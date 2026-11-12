//! One-time bounded reference-mass factorization for the original Verlet map.

use super::super::{MembraneError, MembranePatch, inertia::MembraneMassImage, poll};
use fs_la::factor::{Cholesky, FactorRunError, cholesky_cancellable};
use fs_solver::op::LinearOp;
use fs_sparse::direct::{DirectOrdering, LdltFactor};
use std::{ops::Range, sync::Arc};

mod sparse;
pub use sparse::{SparseConsistentMassBudget, SparseMassStatistics};

/// Bound on each non-interruptible triangular-solve leaf. Factorization itself
/// polls the existing cancellable Cholesky owner. This dense entry does not
/// silently allocate dense 3N²; larger meshes can explicitly select sparse
/// consistent inertia with its separate storage/work allowances.
pub const MAX_CONSISTENT_AXIS_NODES: usize = 512;

/// Explicit dense work bound for consistent membrane inertia.
#[derive(Clone, Copy, Debug)]
pub struct ConsistentMassBudget {
    /// Maximum free scalar nodes in EACH Cartesian component, additionally
    /// bounded by `MAX_CONSISTENT_AXIS_NODES`. Factor storage is at most three
    /// n×n scalar matrices, and identical support patterns share one factor.
    pub max_nodes_per_axis: usize,
}

#[derive(Clone, Copy, Debug)]
pub(super) enum ConsistentPreparation {
    Dense(ConsistentMassBudget),
    Sparse(SparseConsistentMassBudget, DirectOrdering),
}

#[derive(Clone, Debug)]
enum MassFactor {
    Dense(Cholesky),
    Sparse(LdltFactor),
}

impl MassFactor {
    fn solve(&self, values: &mut [f64]) {
        match self {
            Self::Dense(factor) => factor.solve(values),
            Self::Sparse(factor) => factor.solve_in_place(values),
        }
    }
}

#[derive(Clone, Debug)]
struct ScaledFactor {
    // M_ff = D^(1/2) C D^(1/2), with D the physical positive row sums.
    // The diagonal changes solver coordinates; it never replaces the mass.
    root_mass: Vec<f64>,
    factor: MassFactor,
}

#[derive(Clone, Debug)]
pub(super) struct AxisInverse {
    axis: usize,
    range: Range<usize>,
    factor: Arc<ScaledFactor>,
}

#[derive(Clone, Debug)]
pub(super) enum InverseMass {
    Lumped(Vec<f64>),
    Consistent {
        n: usize,
        axes: Vec<AxisInverse>,
        // Sum over all prescribed test rows of M_sf, per free coordinate.
        // This retains net support inertia without storing a dense cross block.
        support_mass: Vec<f64>,
    },
}

impl InverseMass {
    pub(super) fn image(&self) -> MembraneMassImage {
        match self {
            Self::Lumped(_) => MembraneMassImage::Lumped,
            Self::Consistent { .. } => MembraneMassImage::Consistent,
        }
    }

    /// `free` is grouped by Cartesian axis and ordered by original node within
    /// each group. This permits direct allocation-free solves on each output
    /// slice, with no strided gather or interior-mutable work buffer.
    pub(super) fn consistent(
        patch: &MembranePatch,
        free: &[usize],
        fixed: &[bool],
        preparation: ConsistentPreparation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MembraneError> {
        poll(cancelled)?;
        let counts: [usize; 3] =
            std::array::from_fn(|axis| free.iter().filter(|&&dof| dof % 3 == axis).count());
        let cap = match preparation {
            ConsistentPreparation::Dense(budget) => {
                budget.max_nodes_per_axis.min(MAX_CONSISTENT_AXIS_NODES)
            }
            ConsistentPreparation::Sparse(budget, _) => {
                budget.admit()?;
                budget.max_nodes_per_axis
            }
        };
        if cap == 0 || counts.iter().any(|&n| n > cap) {
            return Err(MembraneError::Budget(
                "bounded free scalar nodes per consistent membrane inertia axis",
            ));
        }
        let mut axes: Vec<AxisInverse> = Vec::with_capacity(3);
        let mut first = 0;
        for (axis, &n) in counts.iter().enumerate() {
            poll(cancelled)?;
            let range = first..first + n;
            first += n;
            if n == 0 {
                continue;
            }
            let reused = axes.iter().find(|previous| {
                previous.range.len() == n
                    && free[previous.range.clone()]
                        .iter()
                        .zip(&free[range.clone()])
                        .all(|(a, b)| a / 3 == b / 3)
            });
            let factor = if let Some(previous) = reused {
                Arc::clone(&previous.factor)
            } else if let ConsistentPreparation::Sparse(budget, ordering) = preparation {
                Arc::new(sparse::factor(
                    patch,
                    &free[range.clone()],
                    budget,
                    ordering,
                    cancelled,
                )?)
            } else {
                let mut mass = patch.mass_on_dofs(
                    MembraneMassImage::Consistent,
                    &free[range.clone()],
                    cancelled,
                )?;
                let root_mass: Vec<_> = free[range.clone()]
                    .iter()
                    .map(|&dof| patch.nodal_mass_kg[dof / 3].sqrt())
                    .collect();
                for a in 0..n {
                    poll(cancelled)?;
                    for b in 0..=a {
                        // Divide by the smaller root first, retaining a weak
                        // coupling across widely separated physical masses.
                        let lo = root_mass[a].min(root_mass[b]);
                        let hi = root_mass[a].max(root_mass[b]);
                        let value = (mass[a * n + b] / lo) / hi;
                        mass[a * n + b] = value;
                        mass[b * n + a] = value;
                    }
                }
                let factor =
                    cholesky_cancellable(&mass, n, &mut *cancelled).map_err(
                        |error| match error {
                            FactorRunError::Cancelled => MembraneError::Cancelled,
                            FactorRunError::Factor(_) => MembraneError::Invalid(
                                "consistent membrane inertia is not numerically positive definite",
                            ),
                        },
                    )?;
                Arc::new(ScaledFactor {
                    root_mass,
                    factor: MassFactor::Dense(factor),
                })
            };
            axes.push(AxisInverse {
                axis,
                range,
                factor,
            });
        }
        let mut support_mass = vec![0.0; free.len()];
        if fixed.iter().any(|&held| held) {
            let mut index = vec![usize::MAX; fixed.len()];
            for (i, &dof) in free.iter().enumerate() {
                index[dof] = i;
            }
            for (element_index, element) in patch.elements.iter().enumerate() {
                poll(cancelled)?;
                let mass = patch.element_consistent_mass(element_index)?;
                for axis in 0..3 {
                    for a in 0..3 {
                        if !fixed[3 * element.nodes[a] + axis] {
                            continue;
                        }
                        for b in 0..3 {
                            let i = index[3 * element.nodes[b] + axis];
                            if i != usize::MAX {
                                support_mass[i] += mass[a][b];
                            }
                        }
                    }
                }
            }
        }
        if support_mass.iter().any(|m| !m.is_finite()) {
            return Err(MembraneError::Invalid(
                "unrepresentable consistent support inertia",
            ));
        }
        poll(cancelled)?;
        Ok(Self::Consistent {
            n: free.len(),
            axes,
            support_mass,
        })
    }

    pub(super) fn has_support_coupling(&self) -> bool {
        match self {
            Self::Lumped(_) => false,
            Self::Consistent { support_mass, .. } => support_mass.iter().any(|&m| m != 0.0),
        }
    }

    pub(super) fn sparse_statistics(&self) -> Option<SparseMassStatistics> {
        match self {
            Self::Lumped(_) => None,
            Self::Consistent { axes, .. } => sparse::statistics(axes),
        }
    }

    /// M_sf * delta_v, summed over supported rows and kept by spatial axis.
    /// The caller computes delta_v from M_ff^-1 * delta_p, preserving a small
    /// accepted momentum increment even on a large absolute velocity background.
    pub(super) fn support_impulse(&self, delta_velocity: &[f64]) -> [f64; 3] {
        let mut impulse = [0.0; 3];
        if let Self::Consistent {
            axes, support_mass, ..
        } = self
        {
            for block in axes {
                impulse[block.axis] = support_mass[block.range.clone()]
                    .iter()
                    .zip(&delta_velocity[block.range.clone()])
                    .map(|(mass, delta)| mass * delta)
                    .sum();
            }
        }
        impulse
    }
}

impl LinearOp for InverseMass {
    fn n(&self) -> usize {
        match self {
            Self::Lumped(diagonal) => diagonal.len(),
            Self::Consistent { n, .. } => *n,
        }
    }

    fn apply(&self, momentum: &[f64], velocity: &mut [f64]) {
        match self {
            Self::Lumped(diagonal) => {
                for ((out, p), inverse) in velocity.iter_mut().zip(momentum).zip(diagonal) {
                    *out = p * inverse;
                }
            }
            Self::Consistent { axes, .. } => {
                velocity.copy_from_slice(momentum);
                for block in axes {
                    let values = &mut velocity[block.range.clone()];
                    for (value, root) in values.iter_mut().zip(&block.factor.root_mass) {
                        *value /= root;
                    }
                    block.factor.factor.solve(values);
                    for (value, root) in values.iter_mut().zip(&block.factor.root_mass) {
                        *value /= root;
                    }
                }
            }
        }
    }

    // The inverse of a symmetric mass (lumped, or D^-1/2 M~^-1 D^-1/2 for the
    // consistent one) is symmetric.
    fn apply_transpose(&self, momentum: &[f64], velocity: &mut [f64]) {
        self.apply(momentum, velocity);
    }
}
