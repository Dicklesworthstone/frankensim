//! Original Maxwell histories in the sparse signed midpoint equations.
//!
//! A branch obeys sigma_dot = k b^T q_dot - sigma/tau and stores
//! sigma^2/(2k). Structural-only and pressure-potential-only maps preserve the
//! real symmetric signed system; a mixed map does not and is refused intact.
//! Each retained scalar history is an augmented unknown, never a dense b*b^T
//! update. With nu=tau/(tau+h/2), the unknown is sigma_mid/sqrt(k*nu).
//! Its diagonal is -1 for a structural branch, +1 for an acoustic branch.
//! This preserves physical storage, including dark branches with zero strain.

use super::{
    Coo, Dd, ROUND, Result, SparseStructuralTimeError, SparseStructuralTimeLimits, bound, det,
    finite, poll, push, relative,
};

/// Original dense authored strain maps are read once and retained sparsely.
/// The caller's arrays, allocator metadata and accepted states are additional
/// to these entry counts and the original native factor payload bound.
#[derive(Clone, Copy, Debug)]
pub struct SparseStructuralRelaxationLimits {
    pub maximum_branches: usize,
    /// Bounds scanning the original full [u,psi] arrays before sparse lowering.
    pub maximum_source_entries: usize,
}
impl Default for SparseStructuralRelaxationLimits {
    fn default() -> Self {
        Self {
            maximum_branches: 4096,
            maximum_source_entries: 4_000_000,
        }
    }
}

/// The public adapter borrows the unchanged native MaxwellDampingBranch data.
pub struct BranchSource<'a> {
    pub strain: &'a [f64],
    pub stiffness: f64,
    pub relaxation_time_s: f64,
}

/// One independently contracted physical memory account, in source order.
#[derive(Clone, Copy, Debug, Default)]
pub struct SparseStructuralMemoryWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub energy_change_j: f64,
    pub structural_work_j: f64,
    pub acoustic_work_j: f64,
    /// h*sigma_mid^2/(k*tau), independently evaluated and nonnegative.
    pub dissipation_j: f64,
    pub maximum_relative_residual: f64,
}
impl SparseStructuralMemoryWork {
    pub fn balance_residual_j(self) -> f64 {
        (Dd::from_f64(self.energy_change_j) + Dd::from_f64(self.dissipation_j)
            - Dd::from_f64(self.structural_work_j)
            - Dd::from_f64(self.acoustic_work_j))
        .to_f64()
    }
}

struct Branch {
    direction: Vec<(usize, f64)>,
    domain: usize,
    stiffness: f64,
    root_stiffness: f64,
    relaxation_time_s: f64,
    root_retention: f64,
    retention: f64,
    decay: f64,
    integrated_retention_s: f64,
    coupling: f64,
}

pub(super) struct Branches {
    entries: Vec<Branch>,
    pub(super) acoustic: usize,
}
pub(super) struct Audit {
    pub(super) reaction: Vec<f64>,
    pub(super) memory: Vec<SparseStructuralMemoryWork>,
    pub(super) work: [f64; 2],
    pub(super) energy_change: f64,
    pub(super) dissipation: f64,
    pub(super) relative_residual: f64,
    pub(super) work_scale: [f64; 2],
    pub(super) energy_scale: f64,
}

impl Branches {
    pub(super) fn prepare(
        source: &[BranchSource<'_>],
        n: usize,
        ns: usize,
        half_step: f64,
        source_matrix_entries: usize,
        limits: &SparseStructuralTimeLimits,
        memory_limits: SparseStructuralRelaxationLimits,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        bound(
            "material branches",
            source.len(),
            memory_limits.maximum_branches,
        )?;
        let source_entries =
            n.checked_mul(source.len())
                .ok_or(SparseStructuralTimeError::Invalid(
                    "material source extent overflow",
                ))?;
        bound(
            "material source coefficients",
            source_entries,
            memory_limits.maximum_source_entries,
        )?;
        let mut entries = Vec::with_capacity(source.len());
        let mut retained = source_matrix_entries;
        let mut acoustic = 0;
        for branch in source {
            poll(cancelled)?;
            if branch.strain.len() != n
                || !branch.stiffness.is_finite()
                || branch.stiffness <= 0.0
                || !branch.relaxation_time_s.is_finite()
                || branch.relaxation_time_s <= 0.0
            {
                return Err(SparseStructuralTimeError::Invalid(
                    "original full strain, positive stiffness and physical relaxation time required",
                ));
            }
            let mut direction = Vec::new();
            let mut present = [false; 2];
            for (i, &value) in branch.strain.iter().enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                if !value.is_finite() {
                    return Err(SparseStructuralTimeError::Invalid(
                        "finite original material strain coefficients required",
                    ));
                }
                if value != 0.0 {
                    retained =
                        retained
                            .checked_add(1)
                            .ok_or(SparseStructuralTimeError::Invalid(
                                "material entry extent overflow",
                            ))?;
                    bound(
                        "retained physical entries",
                        retained,
                        limits.maximum_matrix_entries,
                    )?;
                    present[usize::from(i >= ns)] = true;
                    direction.push((i, value));
                }
            }
            if present[0] && present[1] {
                return Err(SparseStructuralTimeError::Invalid(
                    "mixed structural/acoustic Maxwell maps require the original general exact owner",
                ));
            }
            let domain = usize::from(present[1]);
            acoustic += domain;
            // Do not obtain the weak decay by subtracting nearly equal numbers.
            let tau = branch.relaxation_time_s;
            let (retention, decay, integrated, root_retention) = if tau >= half_step {
                let ratio = half_step / tau;
                let den = 1.0 + ratio;
                (
                    1.0 / den,
                    ratio / den,
                    half_step / den,
                    1.0 / det::sqrt(den),
                )
            } else {
                let ratio = tau / half_step;
                let den = 1.0 + ratio;
                (
                    ratio / den,
                    1.0 / den,
                    tau / den,
                    (det::sqrt(tau) / det::sqrt(half_step)) / det::sqrt(den),
                )
            };
            let root_stiffness = det::sqrt(branch.stiffness);
            let coupling = product(half_step, root_stiffness, root_retention).to_f64();
            if !root_retention.is_finite()
                || root_retention <= 0.0
                || !coupling.is_finite()
                || coupling < 0.0
                || (coupling == 0.0 && !direction.is_empty())
            {
                return Err(SparseStructuralTimeError::Invalid(
                    "material midpoint scaling is unrepresentable",
                ));
            }
            entries.push(Branch {
                direction,
                domain,
                stiffness: branch.stiffness,
                root_stiffness,
                relaxation_time_s: tau,
                root_retention,
                retention,
                decay,
                integrated_retention_s: integrated,
                coupling,
            });
        }
        poll(cancelled)?;
        Ok(Self { entries, acoustic })
    }
    pub(super) fn len(&self) -> usize {
        self.entries.len()
    }
    pub(super) fn sign(&self, branch: usize) -> f64 {
        if self.entries[branch].domain == 0 {
            -1.0
        } else {
            1.0
        }
    }
    pub(super) fn assemble(
        &self,
        matrix: &mut Coo,
        spatial: usize,
        limits: &SparseStructuralTimeLimits,
        cancelled: &impl Fn() -> bool,
    ) -> Result<()> {
        for (j, branch) in self.entries.iter().enumerate() {
            poll(cancelled)?;
            let sign = self.sign(j);
            push(matrix, spatial + j, spatial + j, sign, limits)?;
            for (slot, &(i, b)) in branch.direction.iter().enumerate() {
                if slot % 64 == 0 {
                    poll(cancelled)?;
                }
                let value = -sign * branch.coupling * b;
                if value == 0.0 {
                    return Err(SparseStructuralTimeError::Invalid(
                        "nonzero original material coupling underflowed",
                    ));
                }
                push(matrix, i, spatial + j, value, limits)?;
                push(matrix, spatial + j, i, value, limits)?;
            }
        }
        Ok(())
    }
    pub(super) fn rhs(
        &self,
        stresses: &[f64],
        values: &mut [f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<()> {
        for (j, branch) in self.entries.iter().enumerate() {
            if j % 64 == 0 {
                poll(cancelled)?;
            }
            values[j] = self.sign(j)
                * product(
                    branch.root_retention,
                    stresses[j],
                    1.0 / branch.root_stiffness,
                )
                .to_f64();
        }
        poll(cancelled)?;
        Ok(())
    }
    pub(super) fn next_stresses(
        &self,
        initial: &[f64],
        solution: &[f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<Vec<f64>> {
        let mut next = Vec::with_capacity(self.len());
        for (j, branch) in self.entries.iter().enumerate() {
            if j % 64 == 0 {
                poll(cancelled)?;
            }
            let middle = product(branch.root_stiffness, branch.root_retention, solution[j]);
            next.push((Dd::from_f64(2.0) * middle - Dd::from_f64(initial[j])).to_f64());
        }
        finite(&next)?;
        Ok(next)
    }
    pub(super) fn energy(&self, stresses: &[f64], cancelled: &impl Fn() -> bool) -> Result<f64> {
        if stresses.len() != self.len() {
            return Err(SparseStructuralTimeError::Invalid(
                "every original Maxwell stress must be supplied in source order",
            ));
        }
        finite(stresses)?;
        let mut result = Dd::ZERO;
        for (j, branch) in self.entries.iter().enumerate() {
            poll(cancelled)?;
            let z = stresses[j] / branch.root_stiffness;
            result = result + Dd::from_f64(0.5) * Dd::from_f64(z) * Dd::from_f64(z);
        }
        let energy = result.to_f64();
        finite(&[energy])?;
        Ok(energy)
    }
    pub(super) fn audit(
        &self,
        initial: &[f64],
        next: &[f64],
        midpoint_velocity: &[f64],
        duration_s: f64,
        limits: &SparseStructuralTimeLimits,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Audit> {
        let mut result = Audit {
            reaction: Vec::new(),
            memory: Vec::with_capacity(self.len()),
            work: [0.0; 2],
            energy_change: 0.0,
            dissipation: 0.0,
            relative_residual: 0.0,
            work_scale: [0.0; 2],
            energy_scale: 0.0,
        };
        if self.len() == 0 {
            return Ok(result);
        }
        let mut reaction = vec![Dd::ZERO; midpoint_velocity.len()];
        let mut total_change = Dd::ZERO;
        let mut total_heat = Dd::ZERO;
        let mut work = [Dd::ZERO; 2];
        for (j, branch) in self.entries.iter().enumerate() {
            poll(cancelled)?;
            let sigma0 = Dd::from_f64(initial[j]);
            let sigma1 = Dd::from_f64(next[j]);
            let midpoint = (sigma0 + sigma1) * Dd::from_f64(0.5);
            let mut rate = Dd::ZERO;
            let mut rate_scale = 0.0;
            for (slot, &(i, b)) in branch.direction.iter().enumerate() {
                if slot % 64 == 0 {
                    poll(cancelled)?;
                }
                let term = Dd::from_f64(b) * Dd::from_f64(midpoint_velocity[i]);
                rate = rate + term;
                rate_scale += term.to_f64().abs();
                reaction[i] = reaction[i] + Dd::from_f64(b) * midpoint;
            }
            // The ORIGINAL branch ODE, scaled by tau/(tau+h/2). This avoids
            // forming an overflowed h/tau while retaining every physical term.
            let inertia = (sigma1 - sigma0) * Dd::from_f64(branch.retention);
            let decay = midpoint * Dd::from_f64(2.0 * branch.decay);
            let drive = product(
                2.0 * branch.integrated_retention_s,
                branch.stiffness,
                rate.to_f64(),
            );
            let defect = (inertia + decay - drive).to_f64();
            let scale = branch.retention * (initial[j].abs() + next[j].abs())
                + 2.0 * branch.decay * midpoint.to_f64().abs()
                + product(
                    2.0 * branch.integrated_retention_s,
                    branch.stiffness,
                    rate_scale,
                )
                .to_f64();
            let residual = relative(defect, scale)?;
            if residual > limits.maximum_relative_residual {
                return Err(SparseStructuralTimeError::Residual {
                    observed: residual,
                    maximum: limits.maximum_relative_residual,
                });
            }
            let root_k = Dd::from_f64(branch.root_stiffness);
            let z0 = sigma0 / root_k;
            let z1 = sigma1 / root_k;
            let initial_energy = (Dd::from_f64(0.5) * z0 * z0).to_f64();
            let final_energy = (Dd::from_f64(0.5) * z1 * z1).to_f64();
            let change =
                (Dd::from_f64(0.5) * ((sigma1 - sigma0) / root_k) * ((sigma1 + sigma0) / root_k))
                    .to_f64();
            // Independently contract actual sigma^2/(k*tau), not a work deficit.
            let root_heat = product(
                (midpoint / root_k).to_f64(),
                det::sqrt(duration_s),
                1.0 / det::sqrt(branch.relaxation_time_s),
            );
            let heat = (root_heat * root_heat).to_f64();
            let input = product(duration_s, midpoint.to_f64(), rate.to_f64()).to_f64();
            let input_scale = product(duration_s, midpoint.to_f64().abs(), rate_scale).to_f64();
            finite(&[
                initial_energy,
                final_energy,
                change,
                heat,
                input,
                input_scale,
            ])?;
            let mut account = SparseStructuralMemoryWork {
                initial_energy_j: initial_energy,
                final_energy_j: final_energy,
                energy_change_j: change,
                dissipation_j: heat,
                maximum_relative_residual: residual,
                ..SparseStructuralMemoryWork::default()
            };
            if branch.domain == 0 {
                account.structural_work_j = input;
            } else {
                account.acoustic_work_j = input;
            }
            let physical = change.abs() + heat + input.abs();
            let roundoff = initial_energy + final_energy + physical + input_scale;
            let maximum = limits.energy_balance_tolerance * physical + ROUND * roundoff;
            let defect = account.balance_residual_j();
            finite(&[maximum, defect])?;
            if defect.abs() > maximum {
                return Err(SparseStructuralTimeError::EnergyBalance {
                    account: "original Maxwell branch",
                    residual_j: defect,
                    maximum_j: maximum,
                });
            }
            total_change = total_change + Dd::from_f64(change);
            total_heat = total_heat + Dd::from_f64(heat);
            work[branch.domain] = work[branch.domain] + Dd::from_f64(input);
            result.work_scale[branch.domain] += input_scale;
            result.energy_scale += roundoff;
            result.relative_residual = result.relative_residual.max(residual);
            result.memory.push(account);
        }
        result.reaction = reaction.into_iter().map(Dd::to_f64).collect();
        result.work = work.map(Dd::to_f64);
        result.energy_change = total_change.to_f64();
        result.dissipation = total_heat.to_f64();
        finite(&result.reaction)?;
        finite(&result.work)?;
        finite(&result.work_scale)?;
        finite(&[
            result.energy_change,
            result.dissipation,
            result.energy_scale,
        ])?;
        poll(cancelled)?;
        Ok(result)
    }
}

// Three finite factors: multiply the largest and smallest first so an
// avoidable intermediate overflow/underflow does not discard a physical term.
fn product(a: f64, b: f64, c: f64) -> Dd {
    let mut factors = [a, b, c];
    factors.sort_by(|x, y| x.abs().total_cmp(&y.abs()));
    Dd::from_f64(factors[0]) * Dd::from_f64(factors[2]) * Dd::from_f64(factors[1])
}

#[cfg(test)]
mod tests;
