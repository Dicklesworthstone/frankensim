//! Integrated material heat on the original complete held-force trajectory.
//!
//! Each source row is (sqrt(k) b.q-z)/sqrt(tau), not a share of total loss.
//! The original modal/connection/Maxwell coefficients define the observation
//! generator; `PreparedSourceLoss` owns its positive integral. No physical state
//! or second mechanical time stepper is stored here. External acoustic owners
//! and nonlinear/contact trajectories need their own complete source integral.
use super::super::*;
use core::ops::Range;
use fs_phs::{PreparedStepError, RelaxationBranch};
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedSourceLoss};

/// Explicit complete observation resources, aggregated over ALL groups.
/// The 2*n source-integral dimension is at most 256. This bounds observation
/// work, not the original mechanical propagator's independent resources.
#[derive(Clone, Copy, Debug)]
pub struct MaterialHeatBudget {
    pub max_groups: usize,
    pub max_dimension: usize,
    pub max_dense_bytes: usize,
    pub max_multiply_terms: u64,
}

/// Immutable source-group heat on one fixed original operator and clock.
/// Reuse through continuing motion; rate, spring, damping or topology changes
/// invalidate the cache. Evaluation reads the accepted START state only.
pub struct MaterialHeatIntegral {
    signature: Signature,
    ranges: Vec<Range<usize>>,
    integrals: Vec<PreparedSourceLoss>,
    rest_force: Vec<f64>,
}

struct Signature {
    offsets: Vec<usize>,
    modes: Vec<(f64, f64)>,
    columns: Vec<Vec<f64>>,
    damping: Vec<Vec<f64>>,
    laws: Vec<(f64, f64, f64)>,
    branches: Vec<RelaxationBranch>,
    dt: f64,
}
fn damping_column(network: &CoupledModalSystem, i: usize) -> &[f64] {
    network.exact.as_ref().and_then(|e| e.damping_column(i))
        .or_else(|| network.harmonic_damping_columns.get(i).and_then(Option::as_deref))
        .unwrap_or(&network.columns[i])
}
fn prepared_error(e: PreparedStepError) -> ModalCouplingError {
    match e {
        PreparedStepError::Cancelled => ModalCouplingError::Cancelled,
        PreparedStepError::Solver(e) => ModalCouplingError::Propagation(e),
    }
}
fn same_branch(a: &RelaxationBranch, b: &RelaxationBranch) -> bool {
    a.stiffness.to_bits() == b.stiffness.to_bits()
        && a.relaxation_time_s.to_bits() == b.relaxation_time_s.to_bits()
        && a.projection == b.projection
}
impl Signature {
    fn matches(&self, network: &CoupledModalSystem, gate: &CancelGate) -> Result<bool, ModalCouplingError> {
        if network.radiation.is_some() || network.harmonic_only
            || network.offsets != self.offsets || network.dt.to_bits() != self.dt.to_bits()
            || network.mode_count() != self.modes.len()
            || network.connections.len() != self.laws.len()
            || network.relaxation_branches.len() != self.branches.len() {
            return Ok(false);
        }
        for (mode, &(omega, zeta)) in network.models.iter().flat_map(|m| m.modes()).zip(&self.modes) {
            poll(Some(gate))?;
            if mode.angular_frequency_rad_s.to_bits() != omega.to_bits()
                || mode.damping_ratio.to_bits() != zeta.to_bits() { return Ok(false); }
        }
        for (i, law) in network.connections.iter().enumerate() {
            poll(Some(gate))?;
            if (law.stiffness_n_m, law.damping_n_s_m, law.rest_extension_m) != self.laws[i]
                || network.columns[i] != self.columns[i]
                || damping_column(network, i) != self.damping[i] { return Ok(false); }
        }
        for (a, b) in network.relaxation_branches.iter().zip(&self.branches) {
            poll(Some(gate))?;
            if !same_branch(a,b) { return Ok(false); }
        }
        Ok(true)
    }
}

impl CoupledModalSystem {
    /// Prepare disjoint, ordered original material-memory ranges. Empty ranges
    /// are retained as exact zero heat (e.g. an explicitly elastic source layer).
    /// Unselected memories still affect the COMPLETE trajectory but not this
    /// heat account. Every fluid mode, connection and signed cross term survives.
    ///
    /// This observes the same frozen linear generator as ordinary/prepared
    /// playback, including constant-loss damping columns and spring preload.
    /// Radiation/DG/contact owners are NOT approximated by their endpoint forces.
    /// Their attachment requires a different complete-trajectory observation.
    pub fn prepare_material_heat(
        &self, ranges: &[Range<usize>], budget: MaterialHeatBudget, gate: &CancelGate,
    ) -> Result<MaterialHeatIntegral, ModalCouplingError> {
        poll(Some(gate))?;
        self.require_time_domain()?;
        if self.radiation.is_some() || ranges.is_empty() || ranges.len() > budget.max_groups {
            return Err(invalid("material heat needs bounded groups on an original modal/Maxwell trajectory without external acoustic state"));
        }
        let m = self.mode_count();
        let arms = self.relaxation_branches.len();
        let n = m.checked_mul(2).and_then(|n| n.checked_add(arms))
            .ok_or_else(|| invalid("material heat state extent overflow"))?;
        let d = n.checked_mul(2).ok_or_else(|| invalid("material heat observation extent overflow"))?;
        let mut previous_end = 0;
        for range in ranges {
            poll(Some(gate))?;
            if range.start < previous_end || range.start > range.end || range.end > arms {
                return Err(invalid("material heat groups must be ordered disjoint original branch ranges"));
            }
            previous_end = range.end;
        }
        // Include every cache, peak native positive-integral workspace, this
        // complete pencil and its source signature before constructing arrays.
        let signatures = m.checked_mul(2)
            .and_then(|v| self.connections.len().checked_mul(2*m+9).and_then(|k| v.checked_add(k)))
            .and_then(|v| arms.checked_mul(2*m+5).and_then(|k| v.checked_add(k)))
            .and_then(|v| v.checked_add(self.offsets.len()));
        let bytes = d.checked_mul(d).and_then(|v| v.checked_mul(40))
            .and_then(|v| v.checked_mul(ranges.len()))
            .and_then(|v| v.checked_add(n.checked_mul(n)?.checked_mul(8)?))
            .and_then(|v| v.checked_add(signatures?))
            .and_then(|v| v.checked_mul(8));
        let work = (d as u64).checked_pow(3)
            .and_then(|v| v.checked_mul(4096))
            .and_then(|v| v.checked_mul(ranges.len() as u64))
            .and_then(|v| v.checked_add((m as u64).checked_pow(2)?
                .checked_mul((self.connections.len()+arms+1) as u64)?.checked_mul(16)?));
        if n == 0 || d > 256 || d > budget.max_dimension
            || bytes.is_none_or(|b| b > budget.max_dense_bytes)
            || work.is_none_or(|w| w > budget.max_multiply_terms) {
            return Err(invalid("complete material heat exceeds aggregate dimension/byte/work budget"));
        }
        let modes: Vec<_> = self.models.iter().flat_map(|model| model.modes())
            .map(|mode| (mode.angular_frequency_rad_s,mode.damping_ratio)).collect();
        let mut flow = vec![0.0; n*n];
        let mut input = vec![0.0; n*m];
        let mut rest_force = vec![0.0; m];
        for (i, &(omega,zeta)) in modes.iter().enumerate() {
            flow[2*i*n+2*i+1] = 1.0;
            flow[(2*i+1)*n+2*i] = finite(-omega*omega)?;
            flow[(2*i+1)*n+2*i+1] = finite(-2.0*zeta*omega)?;
            input[(2*i+1)*m+i] = 1.0;
        }
        // These are the ORIGINAL complete signed source columns; no stiffness
        // eigendecomposition, fitted modal decay, or loss allocation is used.
        for (index, (column, law)) in self.columns.iter().zip(&self.connections).enumerate() {
            poll(Some(gate))?;
            let damping = damping_column(self,index);
            for i in 0..m {
                rest_force[i] = finite(rest_force[i]+law.stiffness_n_m*law.rest_extension_m*column[i])?;
                for j in 0..m {
                    let row = (2*i+1)*n;
                    flow[row+2*j] = finite(flow[row+2*j]-law.stiffness_n_m*column[i]*column[j])?;
                    flow[row+2*j+1] = finite(flow[row+2*j+1]-law.damping_n_s_m*damping[i]*damping[j])?;
                }
            }
        }
        for (arm, branch) in self.relaxation_branches.iter().enumerate() {
            poll(Some(gate))?;
            let z = 2*m+arm;
            let root = fs_math::det::sqrt(branch.stiffness);
            for i in 0..m {
                let c = finite(root*branch.projection[2*i])?;
                flow[(2*i+1)*n+z] = c;
                flow[z*n+2*i] = finite(c/branch.relaxation_time_s)?;
                for j in 0..m {
                    let slot = (2*i+1)*n+2*j;
                    flow[slot] = finite(flow[slot]-c*root*branch.projection[2*j])?;
                }
            }
            flow[z*n+z] = finite(-1.0/branch.relaxation_time_s)?;
        }
        let native_budget = DenseExponentialBudget { max_dimension: budget.max_dimension,
            max_dense_bytes: budget.max_dense_bytes/ranges.len(),
            max_multiply_terms: budget.max_multiply_terms/ranges.len() as u64 };
        let mut integrals = Vec::with_capacity(ranges.len());
        for range in ranges {
            poll(Some(gate))?;
            let mut rows = vec![0.0; n*range.len()];
            for (row, arm) in range.clone().enumerate() {
                let branch = &self.relaxation_branches[arm];
                let scale = fs_math::det::sqrt(branch.relaxation_time_s);
                let root = fs_math::det::sqrt(branch.stiffness);
                for i in 0..m { rows[row*n+2*i] = finite(root*branch.projection[2*i]/scale)?; }
                rows[row*n+2*m+arm] = finite(-1.0/scale)?;
            }
            integrals.push(PreparedSourceLoss::new_controlled(n,m,&flow,&input,&rows,self.dt,
                native_budget,||gate.is_requested()).map_err(prepared_error)?);
        }
        let signature = Signature { offsets: self.offsets.clone(), modes, columns: self.columns.clone(),
            damping: (0..self.columns.len()).map(|i| damping_column(self,i).to_vec()).collect(),
            laws: self.connections.iter().map(|l| (l.stiffness_n_m,l.damping_n_s_m,l.rest_extension_m)).collect(),
            branches: self.relaxation_branches.clone(), dt: self.dt };
        poll(Some(gate))?;
        Ok(MaterialHeatIntegral { signature, ranges: ranges.to_vec(), integrals, rest_force })
    }
}

impl MaterialHeatIntegral {
    pub fn groups(&self) -> &[Range<usize>] { &self.ranges }
    pub const fn duration_s(&self) -> f64 { self.signature.dt }

    /// True only while the original complete linear coefficients remain intact.
    /// Mechanical/material states and sample clocks may advance normally.
    pub fn is_current(&self, network: &CoupledModalSystem, gate: &CancelGate)
        -> Result<bool, ModalCouplingError> {
        poll(Some(gate))?;
        self.signature.matches(network,gate)
    }

    /// Physical joules in original group order from one accepted START state
    /// under held generalized forces. No clock/state is advanced. A thermal
    /// caller must admit this account with the corresponding native step, not
    /// deposit it and then attempt mechanics. Evaluation has bounded O(G*n²)
    /// scratch/contractions and performs no exponential or nonlinear solve.
    pub fn evaluate(&self, network: &CoupledModalSystem, external: &[f64], gate: &CancelGate)
        -> Result<Vec<f64>, ModalCouplingError> {
        poll(Some(gate))?;
        if !self.signature.matches(network,gate)? {
            return Err(invalid("material heat cache no longer matches complete native coefficients"));
        }
        if external.len() != self.rest_force.len() || external.iter().any(|f| !f.is_finite()) {
            return Err(invalid("material heat requires every finite held modal force"));
        }
        let mut state: Vec<_> = network.models.iter().flat_map(|model| model.states())
            .flat_map(|s| [s.displacement_m_sqrt_kg,s.velocity_m_sqrt_kg_per_s]).collect();
        state.extend_from_slice(&network.material_memory);
        let input = external.iter().zip(&self.rest_force).map(|(f,r)| finite(f+r))
            .collect::<Result<Vec<_>,_>>()?;
        let mut heat = Vec::with_capacity(self.integrals.len());
        for integral in &self.integrals {
            poll(Some(gate))?;
            heat.push(integral.evaluate_controlled(&state,&input,||gate.is_requested()).map_err(prepared_error)?);
        }
        poll(Some(gate))?;
        Ok(heat)
    }
}

#[cfg(test)]
mod tests;
