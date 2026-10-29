//! Continuous joint propagation of the unchanged cavity and downstream DG duct.
use super::*;
pub use crate::render::schedule::force::coupled::dg::ExactDgBudget;

impl SharedCavity {
    /// Prepare one continuous full modal/Maxwell/DG exponential, including the
    /// original finite side branches, rational terminal group and optional
    /// viscothermal wall states on every retained duct.
    /// Every original coordinate and mutual terminal entry is retained. Ordinary
    /// closed/open/matched distal boundaries remain supported. Rational terminals
    /// admit well-conditioned invertible instantaneous D or exactly zero D.
    /// Zero-D terminals (including pressure release and strictly proper dynamic
    /// loads) retain every source state and original pressure constraint; one
    /// differentiated constraint eliminates only the zero-storage flow. Initial
    /// and endpoint pressures and hidden derivative identities must agree.
    /// Mixed-rank singular D, singular/ill-conditioned differentiated constraints
    /// and inconsistent initial states refuse; no epsilon resistance or state
    /// projection is inserted. Gaps, driven mouths,
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
    /// Only authored forces are held. Every original branch state is retained;
    /// pressure continuity and opposite junction flow hold continuously. The
    /// internal pressure and flow evolve together, and their true integral
    /// p(t)Q(t) is evaluated independently on both sides of every junction.
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
