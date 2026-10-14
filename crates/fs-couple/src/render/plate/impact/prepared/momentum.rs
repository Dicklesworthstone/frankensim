//! Cold selection of exact kinetic-coordinate elimination, not a reduced model.
use super::{ImpactError, ImpactSystem, PreparedImpactSystem, StepWorkspace, invalid};

pub(super) fn prepare_plan(inner: &ImpactSystem, workspace: &mut StepWorkspace) -> Result<(), ImpactError> {
    if let Some(air)=&inner.radiation {
        let pairs=air.state_pairs();
        if inner.modes>1 {select_kinematic(inner,workspace,&pairs)?;}
        else {workspace.set_condensed_pairs(&pairs).map_err(ImpactError::PreparedSolve)?;}
    } else if inner.modes>1 {
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
    if !inner.has_nonlinear_dissipation() && diagonal_momentum_flow(inner) {
        // Keep the established scalar traversal when its independent momentum
        // block is suitable. The actual Jacobian is still checked at each step.
        let momenta:Vec<_>=(0..inner.modes).map(|i|2*i+1).collect();
        workspace.set_condensed_leaves(&[],&momenta).map_err(ImpactError::PreparedSolve)
    } else {
        select_kinematic(inner,workspace,&[])
    }
}

fn select_kinematic(inner: &ImpactSystem, workspace: &mut StepWorkspace, pairs: &[[usize;2]]) -> Result<(), ImpactError> {
    // Displacement equations eliminate momenta; momentum equations stay in the
    // border. Thus off-diagonal resistance and nonlinear contact loss are NOT
    // discarded. Acoustic equations may depend on earlier momentum leaves;
    // the shared owner performs exact forward substitution, not independent
    // elimination or a staggered mechanical/acoustic step.
    let pivots:Vec<_>=(0..inner.modes).map(|i|[2*i,2*i+1]).collect();
    workspace.set_condensed_pivots(&pivots,pairs).map_err(ImpactError::PreparedSolve)
}

impl PreparedImpactSystem {
    /// Select exact canonical-momentum elimination for a model without radiation memory.
    ///
    /// The complete displacement/creep/relaxation block and the scalar Gonzalez
    /// correction remain coupled. All momenta are recovered before the original
    /// physical acceptance gate. No mass, contact law, shell nonlinearity, felt
    /// history, cavity overlap or acoustic observation is replaced or discarded.
    ///
    /// Preparation selects this by default with at least two modes. Models with
    /// diagonal momentum resistance and no nonlinear loss keep scalar momentum
    /// leaves. Other models use displacement-equation/momentum-unknown pivots:
    /// spatial dampers and nonlinear contact loss remain fully coupled.
    /// The ACTUAL base Jacobian must have the selected exact block structure on
    /// every update. Singular pivots, unexpected structure or failed original
    /// backward error use full LU of the SAME complete Jacobian.
    /// False selects full dense Newton. Finite differences ignore the selection.
    ///
    /// This method allocates and belongs outside the callback. Switching changes
    /// numerical scratch only, not accepted state, history, time or force ports.
    /// Rounding can differ from dense Newton. No runtime speedup is certified.
    ///
    /// # Errors
    /// Refuses a radiation-loaded model without changing its existing plan:
    /// its joint kinematic/acoustic plan is controlled by `set_radiation_condensation`.
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
