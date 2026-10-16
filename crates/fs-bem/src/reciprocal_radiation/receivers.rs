//! Close pressure/velocity observation of the bound P1/P0 solution.
//! No pressure-mean conversion, auxiliary-density/velocity substitution, new
//! boundary solve, or approximate source direction. The existing near-field
//! owner integrates each original P1 shape and all selected derivative kernels.
use super::*;
use crate::near_field::{Geometry,Options,Observation,VectorObservation,linear::Rows};

/// Reusable exterior points, tied to one immutable variational source geometry.
pub struct ReceiverGeometry<'p,'s> {
    owner:&'p Prepared<'s>,geometry:Geometry<'s>,
}
impl<'s> Prepared<'s> {
    /// Admit actual point-to-surface clearances and exterior winding before
    /// quadrature or solving any source. No enclosing-sphere exclusion is used.
    /// # Errors
    /// The existing near-field geometry, clearance and cardinality refusals.
    pub fn receivers(&self,points:&[Point],minimum_clearance_m:f64)
        ->Result<ReceiverGeometry<'_,'s>,HelmholtzError> {
        Ok(ReceiverGeometry {owner:self,geometry:Geometry::new(self.surface,points,minimum_clearance_m)?})
    }
}
impl<'p,'s> ReceiverGeometry<'p,'s> {
    /// Actual shortest receiver-to-triangle distances [m], in input order.
    pub fn clearances_m(&self)->&[f64] {self.geometry.clearances_m()}
    /// Prepare shared linear-density rows for every source at this frequency.
    /// Set `particle_velocity` for a directional microphone or vector probe.
    /// Each derivative component has its own discrepancy check; pressure-only
    /// convergence does not establish particle-velocity accuracy. Work is
    /// shared across the complete receiver set and polled per integration tile.
    /// # Errors
    /// Frequency/medium, quadrature depth/work, finite arithmetic or cancellation.
    pub fn prepare(&self,k:f64,medium:Medium,options:Options,particle_velocity:bool,
        mut cancelled:impl FnMut()->bool)->Result<PreparedReceivers<'p,'s>,HelmholtzError> {
        poll(&mut cancelled)?;
        if k>self.owner.highest_k {return Err(bad("receiver frequency exceeds the prepared source band"));}
        let rows=self.geometry.prepare_linear(k,medium,options,particle_velocity,cancelled)?;
        Ok(PreparedReceivers {owner:self.owner,rows,k,medium})
    }
}

/// A source-independent frequency row bank, reusable across all modal RHSs.
pub struct PreparedReceivers<'p,'s> {
    owner:&'p Prepared<'s>,rows:Rows,k:f64,medium:Medium,
}
impl PreparedReceivers<'_,'_> {
    /// Actual quadrature points visited during preparation, including rejected
    /// coarse tiles. No additional quadrature is performed by `evaluate`.
    pub fn kernel_evaluations(&self)->usize {self.rows.kernel_evaluations}
    fn evaluate_fields(&self,solution:&Solution<'_,'_>,velocity:bool)
        ->Result<VectorObservation,HelmholtzError> {
        if !core::ptr::eq(self.owner,solution.owner) || self.k.to_bits()!=solution.k.to_bits()
            || self.medium!=solution.medium || solution.coefficients.len()!=self.owner.unknowns() {
            return Err(bad("linear receiver and solution source/frequency/medium disagree"));
        }
        self.rows.evaluate(|i| {
            let f=&self.owner.facets[i];
            (f.ids.map(|a|solution.coefficients[a]),solution.coefficients[self.owner.vertices+i])
        },velocity)
    }
    /// Complex pressure [Pa] from D[p_P1]+S[lambda_P0], with propagated local
    /// integration discrepancies. These are estimates, not continuum/BEM bounds.
    /// # Errors
    /// Cross-wired source/frequency/medium, nonfinite coefficients or overflow.
    pub fn evaluate(&self,solution:&Solution<'_,'_>)->Result<Observation,HelmholtzError> {
        Ok(self.evaluate_fields(solution,false)?.scalar)
    }
    /// Pressure and fluid particle velocity [m/s] of the SAME Green field.
    /// Feed these to `near_field::FirstOrder` for an ideal cardioid or figure
    /// eight; near-field reactive velocity remains in the response.
    /// # Errors
    /// The scalar refusals, missing derivative preparation, or overflow.
    pub fn evaluate_velocity(&self,solution:&Solution<'_,'_>)->Result<VectorObservation,HelmholtzError> {
        self.evaluate_fields(solution,true)
    }
}

#[cfg(test)]
mod tests;
