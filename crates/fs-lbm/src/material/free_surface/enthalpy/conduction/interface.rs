//! Ordered, source-bound thermal resistance on the existing liquid/wall route.

use super::*;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{MaterialPropertySelection, ResolvedInterfaceStatePoint};
use fs_material::thermal::interface::ThermalInterfacePoint;

pub use fs_material::thermal::interface::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY;

/// One original contact source, pinned before any thermal/flow advancement.
///
/// The contact index names the matching entry in the supplied geometry slice.
/// Its interface's A/B order is liquid/solid. The two bulk cards and each
/// caloric chart must agree; matching material names alone do not suffice.
/// Pressure, gap, surface finish and medium remain explicit fixed source data.
/// This is a negligible-storage thermal contact, not a wetting or force law.
#[derive(Clone, Debug)]
pub struct FreeSurfaceBodyInterface<'a> {
    contact_index: usize,
    point: ThermalInterfacePoint<'a>,
    reference: ResolvedInterfaceStatePoint,
}

impl<'a> FreeSurfaceBodyInterface<'a> {
    /// Select once at supported initial liquid/solid states. Every subsequent
    /// substep and endpoint refreshes this exact claim; a support exit refuses
    /// the complete fluid/solid transaction instead of selecting another law.
    pub fn try_new(
        contact_index: usize,
        point: ThermalInterfacePoint<'a>,
        initial_states: [EquilibriumPhaseState; 2],
        selection: MaterialPropertySelection,
    ) -> Result<Self, Error> {
        let reference = point.resolve_resistance(initial_states, selection)?;
        Ok(Self {
            contact_index,
            point,
            reference,
        })
    }

    /// Index into the current step's ordered contact-geometry slice.
    pub const fn contact_index(&self) -> usize {
        self.contact_index
    }

    /// Original selected claim and full ordered surface/state authority.
    pub const fn reference(&self) -> &ResolvedInterfaceStatePoint {
        &self.reference
    }

    pub(super) fn at(
        &self,
        states: [EquilibriumPhaseState; 2],
    ) -> Result<ResolvedInterfaceStatePoint, Error> {
        let point = self.point.at(states)?;
        Ok(self
            .reference
            .at_query_point(self.point.interface, &point)?)
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Evolve sourced interface resistance, finite solid storage, conducting
    /// liquid and native VOF motion on the same accepted timeline.
    ///
    /// Each interface adds R''/A IN SERIES with the original liquid and solid
    /// half-cell resistances. Omitted interface bindings explicitly retain the
    /// zero-contact-resistance approximation. Undeclared geometric contacts
    /// remain insulated. Both participants use their own evolving temperature;
    /// no interface temperature or pressure is guessed. The resulting total
    /// conductance supplies both thermal flux and the original thermal CFL.
    /// Optional environmental gas patches use the existing reservoir source
    /// and combined thermal step bound; their heat is external, while solid
    /// contact exchange stays internal. Pass an empty slice to omit gas exchange.
    /// All source, thermal and native-flow refusals preserve both participants.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_body_interfaces(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        environment: &[FreeSurfaceGasContact<'_>],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        self.advance_with_interfaces_and_environment(
            faces,
            heat,
            controls,
            &[],
            bodies,
            contacts,
            interfaces,
            environment,
            gate,
        )
    }

    /// The same ordered interface coupling including actual molecular stress
    /// power from the original viscous-heating owner. No numerical loss is heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_body_interfaces_and_viscous_heating(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        environment: &[FreeSurfaceGasContact<'_>],
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        let applied = self.viscous_sources(gate)?;
        let conduction = self.advance_with_interfaces_and_environment(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            interfaces,
            environment,
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}
