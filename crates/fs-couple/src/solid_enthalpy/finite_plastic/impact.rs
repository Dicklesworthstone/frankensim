//! Zero-duration contact on the original source-bound caloric/mechanical pair.
use super::{FiniteThermoplasticError, FiniteThermoplasticSolid, poll};
use crate::inelastic_impact::multiple::{
    InelasticMixedContact, InelasticMixedContactCandidate, InelasticMixedContactControls,
    InelasticMixedContactError, thermoplastic_tet_mixed_contact_candidate,
};
use fs_exec::Cx;
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solid::hyper3::plastic::thermal::{FiniteThermoplasticTetError, TetImpulsePoint};
use fs_thermal::enthalpy::DepositionSimulation;

/// Original source/caloric admission or complete native contact refusal.
#[derive(Debug)]
pub enum FiniteThermoplasticImpactError {
    Source(FiniteThermoplasticError),
    Contact(InelasticMixedContactError<FiniteThermoplasticTetError>),
}
impl std::fmt::Display for FiniteThermoplasticImpactError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "source-bound thermoplastic impact refused: {self:?}")
    }
}
impl std::error::Error for FiniteThermoplasticImpactError {}

impl<'a> FiniteThermoplasticSolid<'a> {
    /// Propose simultaneous rigid/solid and rigid/rigid contacts at this clock.
    ///
    /// The original thermal snapshot and exact caloric sources must still match.
    /// The returned source-bound owner retains them, every material history and
    /// the original reference mass, so `advance` continues on the same timeline.
    /// Thermal state is borrowed and unchanged: reported impact loss remains a
    /// mechanical receipt, without an inferred heat-deposition partition.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_contact_candidate(
        &self,
        thermal: &DepositionSimulation<'_>,
        held: &[usize],
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[InelasticMixedContact<TetImpulsePoint>],
        origin: [f64; 3],
        controls: InelasticMixedContactControls,
        cx: &Cx<'_>,
    ) -> Result<InelasticMixedContactCandidate<Self>, FiniteThermoplasticImpactError> {
        poll(cx).map_err(FiniteThermoplasticImpactError::Source)?;
        if thermal.time_s() != self.time_s {
            return Err(FiniteThermoplasticImpactError::Source(
                FiniteThermoplasticError::Invalid("matching impact and thermal clocks required"),
            ));
        }
        self.enthalpies(thermal, true, cx)
            .map_err(FiniteThermoplasticImpactError::Source)?;
        let event = thermoplastic_tet_mixed_contact_candidate(
            &self.solid,
            held,
            bodies,
            properties,
            contacts,
            origin,
            controls,
            cx,
        )
        .map_err(FiniteThermoplasticImpactError::Contact)?;
        let solid = Self {
            solid: event.solid,
            cells: self.cells.clone(),
            time_s: self.time_s,
        };
        poll(cx).map_err(FiniteThermoplasticImpactError::Source)?;
        Ok(InelasticMixedContactCandidate {
            solid,
            rigid_bodies: event.rigid_bodies,
            report: event.report,
        })
    }
}
