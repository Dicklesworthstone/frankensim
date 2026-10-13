//! Cold selection of exact kinetic-coordinate elimination, not a reduced model.
use super::{ImpactError, ImpactSystem, PreparedImpactSystem, StepWorkspace, invalid};

pub(super) fn prepare_plan(inner: &ImpactSystem, workspace: &mut StepWorkspace) -> Result<(), ImpactError> {
    if let Some(air)=&inner.radiation {
        // Mechanical momentum exchanges power with the acoustic pairs. They
        // are not mutually independent leaves; keep the existing loaded image.
        workspace.set_condensed_pairs(&air.state_pairs()).map_err(ImpactError::PreparedSolve)?;
    } else if inner.modes>1 && !inner.has_nonlinear_dissipation() && diagonal_momentum_flow(inner) {
        select(inner,workspace)?;
    }
    Ok(())
}

fn diagonal_momentum_flow(inner: &ImpactSystem) -> bool {
    let sys=&inner.system;let n=sys.state_dim();let (j,r,_)=sys.structure();
    for row in 0..inner.modes {for col in 0..inner.modes {
        let index=(2*row+1)*n+2*col+1;
        if row!=col && j[index]-r[index]!=0.0 {return false;}
    }}
    true
}

fn select(inner: &ImpactSystem, workspace: &mut StepWorkspace) -> Result<(), ImpactError> {
    let momenta:Vec<_>=(0..inner.modes).map(|i|2*i+1).collect();
    workspace.set_condensed_leaves(&[],&momenta).map_err(ImpactError::PreparedSolve)
}

impl PreparedImpactSystem {
    /// Select exact canonical-momentum elimination for a model without radiation memory.
    ///
    /// The complete displacement/creep/relaxation block and the scalar Gonzalez
    /// correction remain coupled. All momenta are recovered before the original
    /// physical acceptance gate. No mass, contact law, shell nonlinearity, felt
    /// history, cavity overlap or acoustic observation is replaced or discarded.
    ///
    /// Preparation selects this by default when there are at least two modes,
    /// no nonlinear dissipative port and no off-diagonal momentum J-R entries.
    /// That cold screen avoids known unsuitable layouts; the ACTUAL analytic
    /// Jacobian is still checked on every update. Even the smallest nonzero
    /// inter-leaf term requires the original full dense solve. Explicit enabling
    /// may attempt other unloaded layouts, but cannot bypass that exact check.
    /// False selects full dense Newton. Finite differences ignore the selection.
    ///
    /// This method allocates and belongs outside the callback. Switching changes
    /// numerical scratch only, not accepted state, history, time or force ports.
    /// Rounding can differ from dense Newton. No runtime speedup is certified.
    ///
    /// # Errors
    /// Refuses a radiation-loaded model without changing its existing plan:
    /// the acoustic-pair plan remains controlled by `set_radiation_condensation`.
    /// Also returns the shared workspace's capacity/dimension refusal.
    pub fn set_momentum_condensation(&mut self, enabled: bool) -> Result<(), ImpactError> {
        if self.inner.radiation.is_some() {
            return Err(invalid("momentum condensation cannot replace a radiation-memory plan"));
        }
        if enabled {select(&self.inner,&mut self.workspace)}
        else {self.workspace.set_condensed_leaves(&[],&[]).map_err(ImpactError::PreparedSolve)}
    }
}

#[cfg(test)]
#[path="momentum_tests.rs"]
mod tests;
