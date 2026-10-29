//! Material-derived string clearance/contact, not a speaking-length edit.
//! These fixed normal contacts use the existing transverse centreline force
//! map, with zero first-order torque on a circular string. Bow slip retains its
//! independently selected surface/twist map. Unilateral engagement and release
//! are solved together with all bows, other strings, body and cavity reactions.

use super::*;
use super::performance::force::{StringForceCoupling, StringForcePort};
use crate::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirection;
use crate::bowed_string::runtime::reciprocal::bank::residual::{ResidualStop, directions, residual_bridges, stop_residual};
use crate::render::RenderError;
use crate::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::geometry::ObstacleGeometry;
use fs_dcontact::Obstacle;

/// Prescribed migration along a retained material string, with contact work.
pub mod stations;

/// Normal closure direction, either relative to the bow axis or independently declared.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StringStopSide {
    /// Positive relative transverse motion closes the contact.
    Positive,
    /// Negative relative transverse motion closes the contact.
    Negative,
    /// Independent unit closure normal in the string's transverse frame. The
    /// receiver shape, if any, must measure displacement along this SAME normal;
    /// its supplied signs are retained. No receiver orientation is inferred.
    Along(TransverseDirection),
}
impl StringStopSide {
    pub(super) fn sign(self) -> f64 {
        match self { Self::Positive | Self::Along(_) => 1.0, Self::Negative => -1.0 }
    }
}

/// A source normal contact; its compliance is a physical law, not a pitch knob.
#[derive(Clone, Debug)]
pub struct MaterialStringStop {
    /// Actual string and interior station along its full supported span.
    pub location: StringBowLocation,
    /// Closure direction. Positive/Negative apply the same sign to both sides;
    /// Along requires the receiver's original map along that declared normal.
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
fn source_error(error: ReciprocalStringError) -> RenderError {
    match error {
        ReciprocalStringError::Material(e) => RenderError::Voice(e),
        ReciprocalStringError::Network(e) => RenderError::Coupled(e),
        ReciprocalStringError::Friction(e) => RenderError::Friction(e),
    }
}

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
    /// performer. Static preload includes the original normal contact potentials
    /// in the full-network equilibrium, with every bow released before the window;
    /// bound geometric/material storage participates in the same static objective.
    /// A vibrating string may leave and re-enter a stop without a state reset.
    ///
    /// A declared normal can differ from the bow tangent or plucking direction.
    /// The same source basis projects it into both selected transverse planes;
    /// a missing plane refuses, rather than silently dropping force. Already
    /// bound geometric stretching stays in the same joint constitutive solve.
    /// Prescribed support translation follows the declared fixed normals.
    /// `with_stop_station_motion` also admits migration along the original span,
    /// retaining the complete source basis and accounting for attachment work.
    /// This does not solve finger-body motion, changing area, contact friction,
    /// rigid constraints or contact-induced adaptivity.
    /// Material source domains and spatial/time refinement remain necessary.
    pub fn with_string_stops(mut self, stops: Vec<MaterialStringStop>,
        limits: FrictionObstacleConfig, gate: &CancelGate) -> Result<Self, RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        if stops.is_empty() { return Ok(self); }
        if self.fiber_support_motion_limits().is_some() {
            return Err(error("normal stops bind before fiber support motion"));
        }
        if !self.source.stops.is_empty() || self.system.samples_rendered() != 0
            || !(1..=32).contains(&limits.max_contacts) || stops.len() > limits.max_contacts
            || stops.len()+self.system.port_count() > 64 {
            return Err(error("string stops require one bounded binding before the first physical sample"));
        }
        let mut geometry = self.stop_geometry(&stops, None, limits.max_setup_terms, gate)?;
        let residual = if self.source.bank().residual_flexibility() {
            Some(self.residual_stop_geometry(&stops, &mut geometry, gate)?)
        } else { None };
        let contacts = stops.iter().zip(geometry).map(|(stop, geometry)|
            (ModalContact { left: geometry.left, right: geometry.right, law: stop.law.clone() }, stop.limits)).collect();
        self.system = self.system.with_obstacles(contacts, limits, gate).map_err(RenderError::Friction)?;
        if let Some(compliance) = residual {
            self.system = self.system.with_obstacle_residual_compliance(compliance, gate).map_err(RenderError::Friction)?;
        }
        self.source.stops = stops;
        Ok(self)
    }

    /// Retain complete normal-contact material sources in original stop order.
    /// Each source must match its already bound unit-point obstacle exactly.
    /// The shared joint owner checks its real closure/rate envelope before every
    /// accepted step, including prescribed support translation and migration.
    pub fn with_string_stop_sources(mut self,
        sources: Vec<(usize, fs_contact::normal_patch::source::SourcePointContact)>,
        gate: &CancelGate) -> Result<Self, RenderError>
    {
        self.system = self.system.with_obstacle_sources(sources, gate)
            .map_err(RenderError::Friction)?;
        Ok(self)
    }

    // With residual flexibility (pexfp.299): each fixed stop's closure gains its
    // bridge's cross term, minus k_eff G_sb times the bridge link's string and
    // body columns, and the stops share the condensed omitted-mode compliance
    // (bank/residual.rs). Receiver stops are not wired. Returns that compliance.
    fn residual_stop_geometry(&self, stops: &[MaterialStringStop], geometry: &mut [ObstacleGeometry],
        gate: &CancelGate) -> Result<Vec<f64>, RenderError>
    {
        let strings = self.source.bank().strings();
        let string_count = self.cards.len();
        let mut ports = Vec::with_capacity(stops.len());
        for stop in stops {
            if stop.body.is_some() {
                return Err(error("a receiver stop is not wired with residual flexibility"));
            }
            let member = strings.get(stop.location.string).ok_or_else(|| error("stop names an unknown string"))?;
            let unsigned = match stop.side {
                StringStopSide::Along(direction) => direction.components(),
                _ => directions(&member.source).0,
            };
            let sign = stop.side.sign();
            // A directed bridge condenses on its own (bank/polarization.rs); a stop it is not
            // orthogonal to would need their omitted-mode cross term, which is not wired.
            if self.source.bank().translational_body_connections().iter()
                .filter(|t| t.string == stop.location.string).map(|t| t.direction.components())
                .any(|[tx, ty]| (tx*unsigned[0] + ty*unsigned[1]).abs() > 1e-12) {
                return Err(error("with residual flexibility a stop must be orthogonal to its string's directed bridges"));
            }
            ports.push(ResidualStop { string: stop.location.string, station_fraction: stop.location.station_fraction,
                direction: [sign*unsigned[0], sign*unsigned[1]] });
        }
        let bridges = residual_bridges(strings, &self.cards).map_err(source_error)?;
        let (cross, compliance) = stop_residual(&bridges, &ports).map_err(source_error)?;
        for ((stop, geometry), factor) in stops.iter().zip(geometry.iter_mut()).zip(cross) {
            if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
            let index = stop.location.string;
            let connection = &strings[index].connections[0];
            let bridge = strings[index].source.bridge_shape(&self.cards[index], connection.station_fraction)
                .map_err(source_error)?;
            if bridge.len() != geometry.left.shapes.len() {
                return Err(error("stop and bridge maps differ from the string's actual basis"));
            }
            for (left, b) in geometry.left.shapes.iter_mut().zip(&bridge) { *left -= factor*b; }
            // Closure is left minus right, so the body part +factor*B enters as right = -factor*B.
            geometry.right = ModalAttachment { component: string_count + connection.body.component,
                shapes: connection.body.shapes.iter().map(|b| -factor*b).collect() };
        }
        Ok(compliance)
    }

    // One source projection for initial and moving stops. Material/gas-loaded
    // mass, both transverse planes and the original receiver basis own these
    // maps; no sine basis, radius/torque or receiving mass is invented here.
    pub(super) fn stop_geometry(&self, stops: &[MaterialStringStop], stations: Option<&[f64]>,
        max_terms: usize, gate: &CancelGate) -> Result<Vec<ObstacleGeometry>, RenderError>
    {
        let ports: Vec<_> = stops.iter().enumerate().map(|(i, stop)| StringForcePort {
            location: StringBowLocation { string: stop.location.string,
                station_fraction: stations.map_or(stop.location.station_fraction, |s| s[i]) },
            coupling: match stop.side {
                StringStopSide::Along(direction) => StringForceCoupling::CentrelineAlong(direction),
                _ => StringForceCoupling::Centreline,
            },
        }).collect();
        let columns = self.string_force_columns(&ports, max_terms, gate)?;
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
            contacts.push(ObstacleGeometry { left, right });
        }
        Ok(contacts)
    }
}
