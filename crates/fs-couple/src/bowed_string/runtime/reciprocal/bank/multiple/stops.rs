//! Material-derived string clearance/contact, not a speaking-length edit.
//! These fixed normal contacts use the existing transverse centreline force
//! map, with zero first-order torque on a circular string. Bow slip retains its
//! independently selected surface/twist map. Unilateral engagement and release
//! are solved together with all bows, other strings, body and cavity reactions.

use super::*;
use super::performance::force::{StringForceCoupling, StringForcePort};
use crate::render::RenderError;
use crate::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_dcontact::Obstacle;

/// Side of the source's centreline/bow axis on which the stop is situated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StringStopSide {
    /// Positive relative transverse motion closes the contact.
    Positive,
    /// Negative relative transverse motion closes the contact.
    Negative,
}
impl StringStopSide {
    fn sign(self) -> f64 { match self { Self::Positive => 1.0, Self::Negative => -1.0 } }
}

/// A fixed normal contact; its compliance is a physical law, not a pitch knob.
#[derive(Clone, Debug)]
pub struct MaterialStringStop {
    /// Actual string and interior station along its full supported span.
    pub location: StringBowLocation,
    /// Closure direction. The same sign acts on string and optional receiver.
    pub side: StringStopSide,
    /// Optional receiver displacement map in its ORIGINAL supplied-body basis.
    /// Component zero means the first body AFTER all strings. Shapes must come
    /// from that receiver's actual mass-normalized basis. None is a fixed stop.
    /// This adapter does not invent a finger-body mode or a new receiving mass.
    pub body: Option<ModalAttachment>,
    /// Original one-point unit-opening [-1] obstacle, including any receipt.
    /// Its gap is signed unloaded SURFACE clearance [m], not string radius or
    /// desired speaking length. The law producer owns contact-material matching;
    /// bulk E/density alone do not certify tissue/fret compliance or loss.
    pub law: Obstacle,
    /// Existing normal-law solve, penetration and reaction limits.
    pub limits: ModalContactConfig,
}
fn error(what: &'static str) -> RenderError { RenderError::Control { what } }

impl MultiBowStringBankSource {
    /// Original physical stop descriptions, retained through scheduling/output.
    #[must_use]
    pub fn stops(&self) -> &[MaterialStringStop] { &self.stops }
}
impl MultiBowStringBank {
    /// Add fixed unilateral contact at actual material-string locations.
    ///
    /// Reuses `string_force_columns` in centreline mode. Density, geometry and
    /// layer changes therefore alter the physical contact map automatically;
    /// retained twist coordinates receive no invented radial-contact torque.
    /// Nonzero receiver shapes carry reactions back through the same network.
    /// No string length, tension, eigenfrequency or observer gain is reassigned.
    /// Every stop, bow and body shares one accepted mechanical step and clock.
    ///
    /// Stops bind once before any sample, retaining initial motion/potential.
    /// Drive indentation and load-release with the existing physical-force
    /// performer. Its linear static-preload shortcut refuses once stops exist;
    /// it cannot falsely claim equilibrium while ignoring nonlinear contacts.
    /// A vibrating string may leave and re-enter a stop without a state reset.
    ///
    /// This is fixed-geometry contact along the source's centreline/bow axis.
    /// A two-plane source projects that same physical direction into both planes;
    /// no independent stop-normal direction is inferred. It does not solve
    /// finger motion, changing contact area, contact friction, geometric
    /// tension, a moving fret, rigid constraints or contact-induced adaptivity.
    /// Material source domains and spatial/time refinement remain necessary.
    pub fn with_string_stops(mut self, stops: Vec<MaterialStringStop>,
        limits: FrictionObstacleConfig, gate: &CancelGate) -> Result<Self, RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        if stops.is_empty() { return Ok(self); }
        if !self.source.stops.is_empty() || self.system.samples_rendered() != 0
            || !(1..=32).contains(&limits.max_contacts) || stops.len() > limits.max_contacts
            || stops.len()+self.system.port_count() > 64 {
            return Err(error("string stops require one bounded binding before the first physical sample"));
        }
        let ports: Vec<_> = stops.iter().map(|stop| StringForcePort {
            location: stop.location, coupling: StringForceCoupling::Centreline,
        }).collect();
        let columns = self.string_force_columns(&ports, limits.max_setup_terms, gate)?;
        let string_count = self.cards.len();
        let mut contacts = Vec::with_capacity(stops.len());
        for (stop, column) in stops.iter().zip(columns) {
            if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
            let member = stop.location.string;
            let modes = self.system.components()[member].modes().len();
            let offset: usize = self.system.components()[..member].iter().map(|m| m.modes().len()).sum();
            let sign = stop.side.sign();
            let left = ModalAttachment { component: member,
                shapes: column[offset..offset+modes].iter().map(|x| sign*x).collect() };
            let right = match &stop.body {
                None => ModalAttachment { component: member, shapes: vec![0.0; modes] },
                Some(body) => {
                    let component = string_count.checked_add(body.component)
                        .ok_or_else(|| error("stop receiver index overflows"))?;
                    let receiver = self.system.components().get(component)
                        .ok_or_else(|| error("stop names an absent receiving body"))?;
                    if body.shapes.len() != receiver.modes().len() || body.shapes.iter().any(|x| !x.is_finite())
                        || !body.shapes.iter().any(|x| *x != 0.0) {
                        return Err(error("stop receiver needs a nonzero finite shape in its actual body basis"));
                    }
                    ModalAttachment { component, shapes: body.shapes.iter().map(|x| sign*x).collect() }
                },
            };
            contacts.push((ModalContact { left, right, law: stop.law.clone() }, stop.limits));
        }
        self.system = self.system.with_obstacles(contacts, limits, gate).map_err(RenderError::Friction)?;
        self.source.stops = stops;
        Ok(self)
    }
}
