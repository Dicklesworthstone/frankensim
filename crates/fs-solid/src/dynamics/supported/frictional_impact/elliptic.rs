//! Original directional material friction with prescribed support reactions.

use super::*;
use crate::dynamics::frictional_impact::elliptic::{SolidEllipticContact, source};
use fs_contact::frictional_multicontact::elliptic::{
    EllipticContactResiduals, EllipticFriction, EllipticFrictionalContactResult,
};
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticState;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseError;
use source::{
    SolidEllipticSourceBinding, SolidEllipticSourceError, SolidImpactGeometry,
    SolidImpactKinematics, SourceProposal,
};

/// Original physical ellipse and complete supported momentum/work response.
/// Only the law owner's determinacy witness uses its explicitly mapped basis.
pub type SupportedSolidEllipticImpactReport =
    SupportedSolidFrictionalImpactReport<EllipticFrictionalContactResult, EllipticContactResiduals>;

/// Source and native support failures retain their original typed cause.
/// Neither the body nor any accepted material history changes on refusal.
#[derive(Debug)]
pub enum SupportedSolidEllipticSourceError<E = SourceImpulseError> {
    Contact(SolidEllipticSourceError<E>),
    Solid(SupportedSolidError),
}
impl<E> From<SolidEllipticSourceError<E>> for SupportedSolidEllipticSourceError<E> {
    fn from(error: SolidEllipticSourceError<E>) -> Self {
        Self::Contact(error)
    }
}
impl<E> From<SupportedSolidFrictionalImpactError> for SupportedSolidEllipticSourceError<E> {
    fn from(error: SupportedSolidFrictionalImpactError) -> Self {
        match error {
            SupportedSolidFrictionalImpactError::Contact(error) => {
                Self::Contact(SolidEllipticSourceError::Impact(error))
            }
            SupportedSolidFrictionalImpactError::Solid(error) => Self::Solid(error),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SupportedSolidEllipticSourceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source-bound supported directional impact refused: {self:?}"
        )
    }
}
impl<E: core::error::Error + 'static> core::error::Error for SupportedSolidEllipticSourceError<E> {}

impl SupportedSpectralSolid {
    /// Resolve simultaneous material-directional impacts while supports retain
    /// their actual accepted instantaneous velocity. Geometry covectors include
    /// every represented coordinate, including the prescribed coordinates.
    ///
    /// The original full-mass response keeps normal/shear and intercontact
    /// coupling. The physical ellipse's impulse returns through the original
    /// transpose; direct support loads and M_cf Delta-v_f both contribute to
    /// the actuator reaction. A zero principal coefficient retains exact free
    /// slip, and equal coefficients recover the original circular owner.
    ///
    /// The actual native endpoint passes the original unscaled physical law,
    /// momentum and full work gates before either trajectory is accepted.
    /// Displacement, time, acceleration, preparation and material history stay
    /// unchanged. This is an event law, not collision discovery, a sustained
    /// contact force, a static/kinetic transition or a physical heat partition.
    pub fn apply_elliptic_impacts(
        &mut self,
        contacts: &[SolidEllipticContact<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> core::result::Result<SupportedSolidEllipticImpactReport, SupportedSolidFrictionalImpactError>
    {
        let proposal = frictional_impact::propose(
            &self.velocity,
            contacts,
            controls,
            |impulse, cx| self.momentum_response(impulse, cx),
            cx,
        )?;
        let (contact, trial) = self.trial_frictional_contact(proposal, contacts, controls, cx)?;
        let impulse = self.accept_momentum_impulse(trial);
        Ok(SupportedSolidEllipticImpactReport { contact, impulse })
    }

    /// Sample each active source at the ACTUAL incoming relative velocity,
    /// including prescribed-support motion, in canonical contact-ID order.
    ///
    /// The original source transaction retains geometry/history input order,
    /// never queries inactive contacts, and stages histories without Clone.
    /// Both histories and all native momentum/work state publish only after
    /// the complete supported candidate passes. Callback-owned allocation,
    /// evaluation and external effects retain their original ownership; no
    /// impulse duration, pressure, sliding travel or thermal state is inferred.
    pub fn apply_elliptic_impacts_with<S, E>(
        &mut self,
        contacts: &[SolidImpactGeometry<'_>],
        histories: &mut [S],
        controls: SolidFrictionalImpactControls,
        resolve: impl FnMut(
            &S,
            &SolidImpactKinematics,
        ) -> core::result::Result<(EllipticFriction, S), E>,
        cx: &Cx<'_>,
    ) -> core::result::Result<
        SupportedSolidEllipticImpactReport,
        SupportedSolidEllipticSourceError<E>,
    > {
        let SourceProposal {
            impulse,
            contacts,
            mut candidates,
        } = source::prepare(
            &self.velocity,
            contacts,
            histories,
            controls,
            |row, cx| self.momentum_response(row, cx),
            resolve,
            cx,
        )?;
        let (contact, trial) = self.trial_frictional_contact(impulse, &contacts, controls, cx)?;
        let impulse = self.accept_momentum_impulse(trial);
        source::publish(histories, &mut candidates);
        Ok(SupportedSolidEllipticImpactReport { contact, impulse })
    }

    /// Use the original ordered-interface law and its retained material axes.
    ///
    /// Both principal sources, explicit temperature/pressure/external query
    /// and impact coefficient role retain the original source sampler's
    /// applicability and authority. Actual incoming slip replaces only speed
    /// query coordinates. Accepted travel, slip, transported texture and the
    /// last finite-step temperature survive this instantaneous refresh.
    pub fn apply_sourced_elliptic_impacts(
        &mut self,
        contacts: &[SolidImpactGeometry<'_>],
        histories: &mut [SourceThermalEllipticState],
        sources: &[SolidEllipticSourceBinding<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> core::result::Result<SupportedSolidEllipticImpactReport, SupportedSolidEllipticSourceError>
    {
        source::admit_bindings(contacts, sources, cx)?;
        self.apply_elliptic_impacts_with(
            contacts,
            histories,
            controls,
            |history, query| source::sample(sources[query.original_index], history, query, cx),
            cx,
        )
    }
}
