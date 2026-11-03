//! Zero-duration contact on the original source-bound caloric/mechanical pair.
use super::{
    FiniteThermoplasticError, FiniteThermoplasticMaterial, FiniteThermoplasticSolid, poll,
};
use crate::finite_solid::contact::mixed::FiniteSolidMixedContactKinematics;
use crate::inelastic_impact::directional::{
    SourceInelasticMixedContact, SourceInelasticMixedContactBinding,
    SourceInelasticMixedContactCandidate, SourceInelasticMixedContactError,
    thermoplastic_tet_source_elliptic_contact_candidate,
};
use crate::inelastic_impact::multiple::elliptic::{
    InelasticMixedEllipticCandidate, InelasticMixedEllipticContact,
    InelasticMixedEllipticContactError, thermoplastic_tet_mixed_elliptic_contact_candidate,
};
use crate::inelastic_impact::multiple::{
    InelasticMixedContact, InelasticMixedContactCandidate, InelasticMixedContactControls,
    InelasticMixedContactError, thermoplastic_tet_mixed_contact_candidate,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState,
};
use fs_exec::Cx;
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticTetError, FiniteThermoplasticTetSolid, TetImpulsePoint,
};
use fs_thermal::enthalpy::DepositionSimulation;

/// Original source/caloric admission or complete native contact refusal.
#[derive(Debug)]
pub enum FiniteThermoplasticImpactError<C = InelasticMixedContactError<FiniteThermoplasticTetError>>
{
    Source(FiniteThermoplasticError),
    Contact(C),
}
impl<C: std::fmt::Debug> std::fmt::Display for FiniteThermoplasticImpactError<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "source-bound thermoplastic impact refused: {self:?}")
    }
}
impl<C: std::fmt::Debug> std::error::Error for FiniteThermoplasticImpactError<C> {}

pub type FiniteThermoplasticEllipticImpactError =
    FiniteThermoplasticImpactError<InelasticMixedEllipticContactError<FiniteThermoplasticTetError>>;
pub type FiniteThermoplasticSourceImpactError<E> = FiniteThermoplasticImpactError<
    SourceInelasticMixedContactError<FiniteThermoplasticTetError, E>,
>;

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
        self.admit_impact(thermal, cx)
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
        let candidate = self.retain_impact(event);
        poll(cx).map_err(FiniteThermoplasticImpactError::Source)?;
        Ok(candidate)
    }

    /// Directional contact on the same accepted source/caloric trajectory.
    /// The supplied ellipse acts in each contact's admitted tangent frame;
    /// zero principal coefficients retain the original exact free-slip law.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_elliptic_contact_candidate(
        &self,
        thermal: &DepositionSimulation<'_>,
        held: &[usize],
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[InelasticMixedEllipticContact],
        origin: [f64; 3],
        controls: InelasticMixedContactControls,
        cx: &Cx<'_>,
    ) -> Result<InelasticMixedEllipticCandidate<Self>, FiniteThermoplasticEllipticImpactError> {
        self.admit_impact(thermal, cx)
            .map_err(FiniteThermoplasticImpactError::Source)?;
        let event = thermoplastic_tet_mixed_elliptic_contact_candidate(
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
        let candidate = self.retain_impact(event);
        poll(cx).map_err(FiniteThermoplasticImpactError::Source)?;
        Ok(candidate)
    }

    /// Resolve original interface sources beside the accepted caloric state.
    ///
    /// Source coefficients use actual incoming contact slip and an explicit
    /// coefficient role. The callback supplies the admitted local conditions;
    /// impulse does not define pressure or duration. Returned interface samples
    /// and mechanical state belong to one candidate. Thermal/material history,
    /// elapsed time and the source's accepted slip path remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_source_elliptic_contact_candidate<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        held: &[usize],
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[SourceInelasticMixedContact],
        histories: &[SourceThermalEllipticState],
        bindings: &[SourceInelasticMixedContactBinding<'_>],
        origin: [f64; 3],
        controls: InelasticMixedContactControls,
        conditions: impl FnMut(
            &SourceThermalEllipticState,
            &FiniteSolidMixedContactKinematics,
        ) -> Result<SourceFrictionConditions, E>,
        cx: &Cx<'_>,
    ) -> Result<SourceInelasticMixedContactCandidate<Self>, FiniteThermoplasticSourceImpactError<E>>
    {
        self.admit_impact(thermal, cx)
            .map_err(FiniteThermoplasticImpactError::Source)?;
        let event = thermoplastic_tet_source_elliptic_contact_candidate(
            &self.solid,
            held,
            bodies,
            properties,
            contacts,
            histories,
            bindings,
            origin,
            controls,
            conditions,
            cx,
        )
        .map_err(FiniteThermoplasticImpactError::Contact)?;
        let candidate = SourceInelasticMixedContactCandidate {
            mechanical: self.retain_impact(event.mechanical),
            sources: event.sources,
        };
        poll(cx).map_err(FiniteThermoplasticImpactError::Source)?;
        Ok(candidate)
    }

    fn admit_impact(
        &self,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> Result<(), FiniteThermoplasticError> {
        poll(cx)?;
        if thermal.time_s() != self.time_s {
            return Err(FiniteThermoplasticError::Invalid(
                "matching impact and thermal clocks required",
            ));
        }
        self.enthalpies(thermal, true, cx)?;
        Ok(())
    }

    fn retain_impact<R>(
        &self,
        event: InelasticMixedContactCandidate<
            FiniteThermoplasticTetSolid<'a, FiniteThermoplasticMaterial<'a>>,
            R,
        >,
    ) -> InelasticMixedContactCandidate<Self, R> {
        InelasticMixedContactCandidate {
            solid: Self {
                solid: event.solid,
                cells: self.cells.clone(),
                time_s: self.time_s,
            },
            rigid_bodies: event.rigid_bodies,
            report: event.report,
        }
    }
}
