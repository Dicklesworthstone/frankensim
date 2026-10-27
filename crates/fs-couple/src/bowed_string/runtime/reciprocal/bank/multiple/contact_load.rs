//! Geometry-checked composition of material-string normal and bow contacts.
use super::*;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::contact::friction::multiple::normal_load::ContactNormalLoad;

/// The same physical surface's tangent and normal in the retained instrument.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BowContactLoad {
    /// Existing bow/contact index in `source().bows()` order.
    pub bow: usize,
    /// Existing normal-stop index in `source().stops()` order.
    pub stop: usize,
}
fn error(what: &'static str) -> RenderError { RenderError::Control { what } }

impl MultiBowStringBankSource {
    /// Original geometry-checked tangent-to-normal binding, retained in playback.
    #[must_use]
    pub fn bow_contact_loads(&self) -> &[BowContactLoad] { &self.bow_contact_loads }
}
impl MultiBowStringBank {
    /// Drive bow friction from the corresponding solved normal contact.
    ///
    /// Both contacts must name the SAME material string and station. A two-plane
    /// string and an explicitly perpendicular `StringStopSide::Along` normal
    /// are required; a missing normal displacement component cannot be guessed.
    /// The normal points from string centre towards the prescribed surface.
    /// If torsion is retained, `normal x tangent = +1` in the source's transverse
    /// frame must match the existing positive-radius surface torque. Opposite
    /// torque orientation refuses, rather than silently changing a live port.
    ///
    /// The original normal law supplies its held reaction to the original
    /// friction law in the SAME joint solve. Contact separation releases traction
    /// automatically. Bow drives for mapped contacts must declare zero normal
    /// force (no extra prescribed load), while speed remains physical input.
    /// Prescribed normal translation uses the existing stop-motion path and its
    /// explicit actuator work. Unmapped bows retain their original semantics.
    ///
    /// This adds no material query, modal reduction, bow-hair structure or time
    /// loop. A normal receiver footprint may be attached to an existing body;
    /// tangential surface motion is externally prescribed unless that same
    /// receiving body is bound through `with_bow_body_surfaces`. Normal and
    /// tangential receiver maps may differ, but must name the same body when
    /// both are dynamically resolved. The original law producers retain
    /// responsibility for material correspondence, contact area and validity.
    ///
    /// Bind once at sample zero, possibly after contact-aware preload. All
    /// accepted vibration, gas/body/cavity sources, interfaces and observers
    /// survive. One normal may supply only one scalar tangential friction law.
    pub fn with_bow_contact_loads(mut self, bindings: Vec<BowContactLoad>, gate: &CancelGate)
        -> Result<Self, RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        if bindings.is_empty() { return Ok(self); }
        if self.bow_station_motion_limits().is_some() {
            return Err(error("contact-loaded bows need jointly transported normal patches before bow migration"));
        }
        if !self.source.bow_contact_loads.is_empty() || self.system.samples_rendered() != 0
            || bindings.len() > self.source.bows.len() {
            return Err(error("bow contact loads require one complete binding before playback"));
        }
        let mut resolved = Vec::with_capacity(bindings.len());
        for binding in &bindings {
            if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
            let bow = self.source.bows.get(binding.bow)
                .ok_or_else(|| error("normal-load binding names an absent bow"))?;
            let stop = self.source.stops.get(binding.stop)
                .ok_or_else(|| error("normal-load binding names an absent string stop"))?;
            if bow.location.string != stop.location.string
                || bow.location.station_fraction.to_bits() != stop.location.station_fraction.to_bits() {
                return Err(error("normal and tangential contacts must occupy the same string station"));
            }
            if let Some(surface) = self.system.port(binding.bow).and_then(|port| port.right.as_ref()) {
                if stop.body.as_ref().is_none_or(|body|
                    body.component.checked_add(self.cards.len()) != Some(surface.component)) {
                    return Err(error("resolved normal and tangential surfaces must belong to the same receiving body"));
                }
            }
            let source = &self.source.bank.strings()[bow.location.string].source;
            let MaterialStringSource::WithPolarizations(polarized) = source else {
                return Err(error("contact-driven bowing requires both normal and tangential transverse motion"));
            };
            let stops::StringStopSide::Along(normal) = stop.side else {
                return Err(error("a contact-driven bow requires an explicit normal perpendicular to its tangent"));
            };
            let n = normal.components();
            let t = polarized.directions().bow();
            let dot = n[0]*t[0] + n[1]*t[1];
            if dot.abs() > 1e-12 {
                return Err(error("normal and bow tangent must be perpendicular in the same transverse frame"));
            }
            if source.torsion_source().is_some() {
                let torque_sign = n[0]*t[1] - n[1]*t[0];
                if (torque_sign-1.0).abs() > 1e-12 {
                    return Err(error("normal surface orientation disagrees with the existing torsional bow port"));
                }
            }
            resolved.push(ContactNormalLoad { friction_port: binding.bow, obstacle: binding.stop });
        }
        self.system = self.system.with_contact_normal_loads(resolved, gate).map_err(RenderError::Friction)?;
        self.source.bow_contact_loads = bindings;
        Ok(self)
    }
}
