//! Continuous joint propagation of the unchanged cavity and downstream DG duct.
use super::*;
pub use crate::render::schedule::force::coupled::dg::ExactDgBudget;

impl SharedCavity {
    /// Prepare one continuous full modal/Maxwell/DG exponential. Every original
    /// coordinate and distal closed/open/matched boundary is retained, including
    /// optional viscothermal wall states. Rational terminals, branches, gaps,
    /// driven mouths, direct neck radiation and constant-loss connections refuse.
    ///
    /// Original neck geometry and inlet gas/area/clock admissions remain binding.
    /// This path never prepares the mean-port mobility. The explicit cold budget
    /// covers the complete states + held generalized forces + one constant row,
    /// all cached quadratic work/loss integrals and the existing dense kernel.
    /// Warm work is O(d^2) per retained work/loss channel, with no exponential or
    /// root solve. Original physical/energy ceilings remain binding; a passive
    /// energy bound also limits DG inlet pressure over the complete interval.
    ///
    /// Only authored forces are held. The internal pressure and flow evolve
    /// continuously together, and their true integral p(t)Q(t) is exchanged.
    /// Splitting a constant-drive interval preserves the same trajectory to
    /// numerical error; varying-drive reconstruction, spatial DG convergence and
    /// constitutive-band validity are still separate accuracy obligations.
    pub fn with_exact_downstream_duct(
        self,
        load: DgPortLoad,
        budget: ExactDgBudget,
        gate: &CancelGate,
    ) -> Result<SharedCavityDg, SharedCavityDgError> {
        let column = self.downstream_duct_column(&load, gate)?;
        let bridge = PreparedJointDg::new(&self.system, &load, column, budget, gate)?;
        Ok(SharedCavityDg {
            cavity: self,
            load,
            evolution: Evolution::Exact(bridge),
            frame: None,
        })
    }
}
impl SharedCavityDg {
    /// Switch a cold connection to the complete continuous exact backend. No
    /// accepted physical history may have advanced, and none is reset or fitted.
    pub fn with_exact_evolution(
        mut self,
        budget: ExactDgBudget,
        gate: &CancelGate,
    ) -> Result<Self, SharedCavityDgError> {
        let column = self.cavity.downstream_duct_column(&self.load, gate)?;
        self.evolution = Evolution::Exact(PreparedJointDg::new(
            &self.cavity.system,
            &self.load,
            column,
            budget,
            gate,
        )?);
        Ok(self)
    }
    pub fn uses_exact_evolution(&self) -> bool {
        matches!(&self.evolution, Evolution::Exact(_))
    }
}

#[cfg(test)]
mod tests;
