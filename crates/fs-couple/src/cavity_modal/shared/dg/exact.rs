//! Continuous joint propagation of the unchanged cavity and downstream DG duct.
use super::*;
pub use crate::render::schedule::force::coupled::dg::ExactDgBudget;

impl SharedCavity {
    /// Prepare one continuous full modal/Maxwell/DG exponential, including the
    /// original rational terminal group and optional viscothermal wall states.
    /// Every original coordinate and mutual terminal entry is retained. Ordinary
    /// closed/open/matched distal boundaries remain supported. Rational terminals
    /// require nonsingular, well-conditioned instantaneous D; a singular D needs
    /// a descriptor treatment and is refused. Branches, gaps, driven mouths,
    /// direct neck radiation and constant-loss connections remain unsupported.
    ///
    /// Original neck geometry and inlet gas/area/clock admissions remain binding.
    /// This path never prepares the mean-port mobility. The explicit cold budget
    /// covers the complete states + held generalized forces + one constant row,
    /// all cached quadratic work/loss integrals and the existing dense kernel.
    /// Warm work is O(d^2) per retained work/loss channel, with no exponential or
    /// root solve. Original physical and energy ceilings remain binding. Rational
    /// terminal work is signed boundary exchange; its original reports provide
    /// no generic storage law. The complete conditioned generator supplies a
    /// conservative numerical interval envelope for pressure/reaction guards,
    /// without assigning stored energy or an all-frequency passivity certificate
    /// to terminal coordinates.
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
