//! Exact spectral image of the existing fs-phs Maxwell storage/resistance law.
//! In elastic-memory coordinates e=c.x-z, e_dot=c.v-e/tau. This is a change
//! of physical coordinates, not a fit, time discretization or branch elimination.
use super::*;
use fs_phs::RelaxationBranch;

#[derive(Default)]
pub(super) struct RetainedMemory {
    pub(super) rows: Vec<Vec<f64>>,
    pub(super) times: Vec<f64>,
}

impl BemLoadedModes {
    /// Attach the SAME declared Maxwell arms used by fs-phs. Projections here
    /// have exactly `mode_count()` entries in mass-normalized displacement order,
    /// not the interleaved q/p layout of a time-domain model. The equilibrium K
    /// must not already contain the arms' excess stiffness. A Kelvin viscosity
    /// may coexist only when it represents a distinct physical loss mechanism.
    ///
    /// Each c=sqrt(stiffness)*projection retains its full signed coupling.
    /// No channel merging or rank truncation: even equal-time, cancelling-force
    /// histories keep their internal states. The caller owns matching the source
    /// basis and constitutive applicability, as for the existing viscosity input.
    ///
    /// One cold attachment, at most 256 branches, within both the supplied branch
    /// ceiling and the original FULL dimension/payload/work limits. Additional
    /// states are never traded for boundary panels. Empty input is a bitwise no-op.
    /// Failed admission/cancellation returns no partially attached model.
    pub fn with_relaxation_branches(
        mut self,
        branches: &[RelaxationBranch],
        maximum_branches: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BemModesError> {
        poll(&mut cancelled)?;
        if maximum_branches > 256 || branches.len() > maximum_branches {
            return Err(bad("material memory exceeds its explicit 256-branch ceiling"));
        }
        if branches.is_empty() { return Ok(self); }
        if self.memory_count() != 0 {
            return Err(bad("material memory requires one complete cold attachment"));
        }
        let m = self.mode_count();
        Self::admit_dimensions(m, self.panel_count(), branches.len(), &self.settings)?;
        let mut rows = Vec::with_capacity(branches.len());
        let mut times = Vec::with_capacity(branches.len());
        for branch in branches {
            poll(&mut cancelled)?;
            let tau = branch.relaxation_time_s;
            if branch.projection.len() != m
                || branch.projection.iter().any(|x| !x.is_finite())
                || !branch.stiffness.is_finite() || branch.stiffness <= 0.0
                || !tau.is_finite() || tau <= 0.0 || !(1.0/tau).is_finite()
            {
                return Err(bad("Maxwell arm needs a finite displacement map and positive representable stiffness/time"));
            }
            let root = det::sqrt(branch.stiffness);
            let row: Vec<_> = branch.projection.iter().map(|&b| b*root).collect();
            if row.iter().any(|x| !x.is_finite()) || row.iter().all(|x| *x == 0.0)
                || row.iter().zip(&branch.projection).any(|(&c,&b)| b != 0.0 && c == 0.0)
            {
                return Err(bad("Maxwell arm scale loses a nonzero source component"));
            }
            rows.push(row);
            times.push(tau);
        }
        self.memory = RetainedMemory { rows, times };
        self.rebalance(&mut cancelled)?;
        Ok(self)
    }

    /// c_j [1/s], in the SAME displacement basis; H_j=|c_j.x-z_j|^2/2.
    pub fn relaxation_rows_per_s(&self) -> &[Vec<f64>] { &self.memory.rows }
    pub fn relaxation_times_s(&self) -> &[f64] { &self.memory.times }
}

impl RetainedMemory {
    // O(branches*modes) per evaluation. In particular no dense C^T C and no
    // rational dynamic stiffness are formed. With no arms, not a bit is touched.
    pub(super) fn add_to_pencil(
        &self, omega: C64, m: usize, p: usize,
        matrix: &mut [C64], derivative: &mut [C64], envelope: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(), BemModesError> {
        let n = m+p+self.times.len();
        for (j,(row,&tau)) in self.rows.iter().zip(&self.times).enumerate() {
            poll(cancelled)?;
            let at = m+p+j;
            let iw_tau = C64::new(-omega.im*tau, omega.re*tau);
            for i in 0..m {
                matrix[i*n+at] = C64::from_re(row[i]);
                envelope[i*n+at] = row[i].abs();
                matrix[at*n+i] = iw_tau.scale(row[i]);
                derivative[at*n+i] = C64::new(0.0, tau*row[i]);
                envelope[at*n+i] = omega.abs()*tau*row[i].abs();
            }
            matrix[at*n+at] = C64::ONE-iw_tau;
            derivative[at*n+at] = C64::new(0.0, -tau);
            envelope[at*n+at] = 1.0+omega.abs()*tau;
        }
        Ok(())
    }

    pub(super) fn seed(
        &self, omega: C64, displacement: &[C64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<C64>, BemModesError> {
        let mut values = Vec::with_capacity(self.times.len());
        for (branch,(row,&tau)) in self.rows.iter().zip(&self.times).enumerate() {
            poll(cancelled)?;
            let s_tau = C64::new(omega.im*tau, -omega.re*tau);
            let denominator = C64::ONE+s_tau;
            if denominator.abs() == 0.0 {
                return Err(BemModesError::SeedMemorySingular { branch });
            }
            let strain = row.iter().zip(displacement)
                .fold(C64::ZERO, |sum,(&c,&x)| sum+x.scale(c));
            let value = (s_tau/denominator)*strain;
            if !finite(value) {
                return Err(bad("initial material memory is not representable"));
            }
            values.push(value);
        }
        Ok(values)
    }
}

#[cfg(test)]
#[path = "memory_tests.rs"]
mod tests;
