//! Two transverse planes of ONE circular material string.
//!
//! Both planes solve mu*y_tt = T*y_xx - EI*y_xxxx with the same material/gas
//! losses. For a unit transverse direction d, a physical port has row
//! [d_x*phi, d_y*phi]; its transpose applies force to those same coordinates.
//! Thus rotation changes neither mass nor work. A bridge spring stores
//! k*(d dot y - y_body)^2/2, so body compliance couples the planes and shifts
//! their assembled resonances without an authored detuning or output gain.
//!
//! See Mansour, Woodhouse & Scavone (2016), "Enhanced wave-based modelling of
//! musical strings. Part 2: Bowed strings", for two-polarization/body coupling:
//! https://api.repository.cam.ac.uk/server/api/core/bitstreams/bc2b1bfc-8b61-4499-999e-3de6dd494927/content
//!
//! The section is circular and frozen-state, with identical flexural laws in
//! both planes. Translational bridge ports follow their supplied direction;
//! optional torque/rotation connections use the original material torsion in
//! the same joint solve. The explicit `MaterialStringSource::into_nonlinear`
//! path retains both planes in one Kirchhoff--Carrier elongation/tension channel;
//! the linear bank continues to use its original quadratic storage. A six-DOF
//! termination remains outside this image. Independently directed translation
//! connections can load both planes through distinct body displacement maps;
//! the bridge need not be a rank-one termination.
//! Optional twist carries no fictitious direct acoustic output.

use super::*;
use fs_material::state_point::MaterialStatePointError;

/// A physical unit direction in the string's orthonormal transverse frame.
/// The first component addresses plane one and the second addresses plane two;
/// it is not a modal gain. Admission retains the supplied numbers verbatim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseDirection([f64; 2]);

// Construction excludes NaN, so equality is reflexive for every admitted value.
impl Eq for TransverseDirection {}

impl TransverseDirection {
    /// Direction along the first retained transverse plane.
    pub const X: Self = Self([1.0, 0.0]);
    /// Direction along the second retained transverse plane.
    pub const Y: Self = Self([0.0, 1.0]);

    /// Admit finite components with unit squared norm within 1e-12, without
    /// normalizing or discarding a component. A one-plane source subsequently
    /// refuses any nonzero second component instead of silently projecting it.
    pub fn try_new(direction: [f64; 2]) -> Result<Self, MaterialStatePointError> {
        let norm2 = direction[0] * direction[0] + direction[1] * direction[1];
        if direction.iter().any(|x| !x.is_finite())
            || !norm2.is_finite() || (norm2 - 1.0).abs() > 1e-12
        {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "transverse direction must be a finite unit vector",
            });
        }
        Ok(Self(direction))
    }

    /// Original dimensionless components in the declared transverse frame.
    #[must_use]
    pub const fn components(self) -> [f64; 2] { self.0 }
}

/// A centreline spring/dashpot along an independently declared transverse axis.
///
/// Its relative displacement is `direction dot y_string - x_body - rest`.
/// The same projection applies the equal and opposite force in the coupled
/// solve. Several such laws express a bridge with different compliances in its
/// two transverse directions, without duplicating a string or changing its
/// material. Torque is independent and receives exactly zero from this port.
#[derive(Clone, Debug, PartialEq)]
pub struct StringBodyTranslationConnection {
    /// Original material-string index.
    pub string: usize,
    /// Interior station on the full supported string, including afterlength.
    pub station_fraction: f64,
    /// Signed unit direction in this string's declared transverse frame.
    pub direction: TransverseDirection,
    /// Displacement along that SAME physical axis in the original body basis.
    /// Component indices address the supplied bodies, before strings are added.
    /// A material producer must derive this map from its actual geometry/basis;
    /// matching its length alone does not establish provenance or orientation.
    pub body: ModalAttachment,
    /// Nonnegative stiffness [N/m], independently sourced from the joint.
    pub stiffness_n_m: f64,
    /// Nonnegative viscous resistance [N s/m], independently sourced from the joint.
    pub damping_n_s_m: f64,
    /// Signed stress-free relative displacement [m].
    pub rest_extension_m: f64,
}

impl StringBodyTranslationConnection {
    pub(super) fn lower(
        &self, cards: &[BowedStringCard], members: &[MaterialStringMember],
        bodies: &[ModalAcousticTimeModel],
    ) -> Result<ModalConnection, ReciprocalStringError> {
        let member = members.get(self.string)
            .ok_or_else(|| invalid("directed bridge names an unknown material string"))?;
        let card = cards.get(self.string)
            .ok_or_else(|| invalid("directed bridge has no aligned string realization"))?;
        let body = bodies.get(self.body.component)
            .ok_or_else(|| invalid("directed bridge names an unknown original body"))?;
        if !interior(self.station_fraction)
            || self.body.shapes.len() != body.modes().len()
            || self.body.shapes.iter().any(|x| !x.is_finite())
            || self.body.shapes.iter().all(|x| *x == 0.0)
            || !self.stiffness_n_m.is_finite() || self.stiffness_n_m < 0.0
            || !self.damping_n_s_m.is_finite() || self.damping_n_s_m < 0.0
            || (self.stiffness_n_m == 0.0 && self.damping_n_s_m == 0.0)
            || !self.rest_extension_m.is_finite() {
            return Err(invalid("directed bridge requires an interior station, a nonzero finite body map and a passive law"));
        }
        let shapes = member.source.transverse_shape_along(card, self.station_fraction, self.direction)?;
        Ok(ModalConnection {
            left: ModalAttachment { component: self.string, shapes },
            right: ModalAttachment { component: cards.len() + self.body.component,
                shapes: self.body.shapes.clone() },
            stiffness_n_m: self.stiffness_n_m, damping_n_s_m: self.damping_n_s_m,
            rest_extension_m: self.rest_extension_m,
        })
    }
}

impl ReciprocalStringBank {
    /// Add reciprocal string/body translations with an independent axis per law.
    ///
    /// The existing bridge connections remain. These extra ports can carry the
    /// missing second-plane reaction, or a spatially separated support. Their
    /// string rows use the original material and actual gas-loaded transverse
    /// mass; the supplied body rows stay in that body's retained basis. Thus
    /// normal and lateral bridge motion, sympathetic strings and any torsion
    /// share one force solve and one work ledger. No modal frequency is retuned.
    /// Original source identities, memory, gas, radiation, observer and dry bow
    /// interface survive; subsequent gas/cavity/rotational bindings retain all
    /// directions and recompute their projections from the same sources.
    ///
    /// Bind one complete set before string motion or playback and before adding
    /// nonlinear storage. An empty set is a no-op. Omitted-mode residual bridge
    /// flexibility currently admits one axis only and refuses this composition.
    /// A one-plane string refuses any direction with a nonzero second component.
    /// This does not infer a bridge's anisotropic law, rigid frame, geometry or
    /// unilateral contact from the bulk string material.
    pub fn with_translational_body_connections(
        mut self, connections: Vec<StringBodyTranslationConnection>, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if connections.is_empty() { return Ok(self); }
        if self.system.has_nonlinear_potentials() || self.source.residual_flexibility {
            return Err(invalid("directed bridges bind before nonlinear storage and without scalar residual flexibility"));
        }
        if self.system.samples_rendered() != 0 || !self.source.translational_body_connections.is_empty()
            || self.system.components()[..self.cards.len()].iter().flat_map(|c| c.states())
                .any(|s| s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
            return Err(invalid("directed bridges require one complete binding before string motion or playback"));
        }
        let body_end = self.source.fluid_loading.as_ref().map_or(
            self.system.components().len(), |fluid| self.cards.len() + fluid.body_components,
        );
        let bodies = self.system.components()[self.cards.len()..body_end].to_vec();
        let link_count = self.source.body_connections.len()
            .checked_add(self.source.rotational_body_connections.len())
            .and_then(|count| count.checked_add(connections.len()))
            .ok_or_else(|| invalid("directed bridge connection count overflow"))?;
        preflight(&self.source.strings, &bodies, &self.source.options, link_count)?;
        let bridges: Vec<_> = self.source.strings.iter().map(|s| s.connections.as_slice()).collect();
        let mut system = assemble_sources(&self.cards, Some(&self.source.strings), bodies,
            &bridges, &self.source.body_connections, &connections,
            &self.source.rotational_body_connections, self.source.fluid_loading.as_ref(),
            self.source.prepared_linear, &self.source.constant_loss_bridges(),
            false, self.system.port().law.clone(), &self.source.options, gate)?;
        if let Some(weights) = self.system.acceleration_pressure_weights() {
            system = system.with_acceleration_pressure(weights.to_vec())
                .map_err(ReciprocalStringError::Friction)?;
        }
        if let Some(binding) = self.system.interface() {
            system = system.with_interface(binding.clone()).map_err(ReciprocalStringError::Friction)?;
        }
        system = self.source.restore_material_memory(system, &self.cards,
            Some(self.system.material_memory().to_vec()), gate)?;
        self.source.validate_body_memory_frequencies(&system, gate)?;
        self.source.validate_material_bridge_frequencies(&system, gate)?;
        poll(gate)?;
        self.system = system;
        self.source.translational_body_connections = connections;
        Ok(self)
    }
}

/// Physical directions expressed in one common orthonormal transverse frame.
/// No modal weights, gains, frequency offsets or implicit normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseDirections {
    bow: [f64; 2],
    bridge: [f64; 2],
}
impl TransverseDirections {
    /// Admit the tangent of the bow stroke and the string/body connector axis.
    /// Each vector must have unit squared norm within 1e-12. This tolerance is
    /// only floating-point direction admission; supplied vectors are retained.
    pub fn try_new(bow: [f64; 2], bridge: [f64; 2]) -> Result<Self, MaterialStatePointError> {
        for direction in [bow, bridge] {
            TransverseDirection::try_new(direction)?;
        }
        Ok(Self { bow, bridge })
    }
    /// Bow tangent; default centreline loads use it, explicitly directed loads do not.
    #[must_use]
    pub const fn bow(&self) -> [f64; 2] {
        self.bow
    }
    /// Direction of the physical scalar connection to the common receiver.
    #[must_use]
    pub const fn bridge(&self) -> [f64; 2] {
        self.bridge
    }
}

/// Original circular section, optional original torsion and physical directions.
/// Coordinate order is first transverse plane, second plane, then any twist.
#[derive(Clone, Debug)]
pub struct PolarizedStringSource {
    source: MaterialStringSource,
    directions: TransverseDirections,
    total_modes: usize,
}
impl PolarizedStringSource {
    /// Original source; no repeated material query or copied material identity.
    #[must_use]
    pub const fn source(&self) -> &MaterialStringSource {
        &self.source
    }
    /// Original directions, independent of the arbitrary transverse basis choice.
    #[must_use]
    pub const fn directions(&self) -> TransverseDirections {
        self.directions
    }
    pub(super) const fn total_modes(&self) -> usize {
        self.total_modes
    }

    pub(super) fn add_torsional_relaxation(
        self,
        properties: &[fs_material::section::relaxation::LayerRelaxationProperties],
    ) -> Result<MaterialStringSource, MaterialStatePointError> {
        self.source.with_torsional_relaxation(properties)?
            .with_transverse_polarizations(self.directions)
    }

    pub(super) fn add_torsion(
        self,
        shear: &str,
        viscosity: &str,
        modes: usize,
    ) -> Result<MaterialStringSource, MaterialStatePointError> {
        self.source
            .with_torsion(shear, viscosity, modes)?
            .with_transverse_polarizations(self.directions)
    }
}
impl MaterialStringSource {
    /// Retain both transverse displacement components from the SAME source.
    /// No second material, mass, damping curve or artificial pitch split exists.
    /// The body/friction network acts on the declared directions and determines
    /// its own polarization mixing. Adding torsion before or after this operation
    /// gives the same coordinate order and physical source.
    ///
    /// # Errors
    /// Already polarized/nested sources, zero or oversized transverse basis,
    /// or mode-count overflow. Full network and material-band limits still apply
    /// when the source is lowered into the bank.
    pub fn with_transverse_polarizations(
        self,
        directions: TransverseDirections,
    ) -> Result<Self, MaterialStatePointError> {
        let bending = match &self {
            Self::WithTorsion(t) => t.bending(),
            other => other,
        };
        let count = match bending {
            Self::Homogeneous(s) => s.string().n_modes,
            Self::Wound(s) => s.string().n_modes,
            Self::Concentric { template, .. } | Self::ConcentricProny { template, .. } => template.n_modes,
            Self::WithPolarizations(_) | Self::WithTorsion(_) => {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "two transverse planes require one original unpolarized source",
                });
            }
        };
        if count == 0 || count > 4096 {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "each transverse plane requires 1..=4096 modes",
            });
        }
        let total_modes =
            self.modes()
                .checked_add(count)
                .ok_or(MaterialStatePointError::InvalidDerived {
                    quantity: "combined transverse and torsional dimension overflow",
                })?;
        Ok(Self::WithPolarizations(Box::new(PolarizedStringSource {
            source: self,
            directions,
            total_modes,
        })))
    }

    pub(super) fn torsion_source(&self) -> Option<&torsion::TorsionalStringSource> {
        match self {
            Self::WithPolarizations(p) => p.source.torsion_source(),
            Self::WithTorsion(t) => Some(t),
            Self::Homogeneous(_) | Self::Wound(_) | Self::Concentric { .. } | Self::ConcentricProny { .. } => None,
        }
    }

    pub(super) fn append_extra_modes(
        &self,
        card: &BowedStringCard,
        modes: &mut Vec<ModalAcousticMode>,
        gate: &CancelGate,
    ) -> Result<(), ReciprocalStringError> {
        if matches!(self, Self::WithPolarizations(_)) {
            // An independent displacement component has the SAME physical mass,
            // stiffness and loss; this does not double the specimen's line mass.
            for k in 0..card.mode_count {
                poll(gate)?;
                modes.push(modes[k]);
            }
        }
        if let Some(twist) = self.torsion_source() {
            twist.append_modes(card, modes, gate)?;
        }
        Ok(())
    }

    pub(super) fn transverse_shape(
        &self,
        card: &BowedStringCard,
        station: f64,
    ) -> Result<Vec<f64>, ReciprocalStringError> {
        self.projected_shape(card, station, false)
    }
    pub(super) fn bridge_shape(
        &self,
        card: &BowedStringCard,
        station: f64,
    ) -> Result<Vec<f64>, ReciprocalStringError> {
        self.projected_shape(card, station, true)
    }
    fn projected_shape(
        &self,
        card: &BowedStringCard,
        station: f64,
        bridge: bool,
    ) -> Result<Vec<f64>, ReciprocalStringError> {
        let direction = match self {
            Self::WithPolarizations(p) => {
                if bridge { p.directions.bridge } else { p.directions.bow }
            },
            _ => [1.0, 0.0],
        };
        // Both sources of these components were admitted above. The shared
        // projector keeps existing default-port arithmetic and coordinate order.
        self.transverse_shape_along(card, station, TransverseDirection(direction))
    }

    /// Conjugate centreline displacement/force row for an independent physical
    /// direction. Twist entries are exactly zero: radial normal force has no
    /// first-order torque about a circular string's axis.
    pub(super) fn transverse_shape_along(
        &self,
        card: &BowedStringCard,
        station: f64,
        direction: TransverseDirection,
    ) -> Result<Vec<f64>, ReciprocalStringError> {
        if !interior(station) {
            return Err(invalid("directed string port requires an interior station"));
        }
        let direction = direction.components();
        let polarized = matches!(self, Self::WithPolarizations(_));
        if !polarized && direction[1] != 0.0 {
            return Err(invalid("a second-plane load requires a two-plane material string"));
        }
        let mut row = shape(card, station);
        if polarized {
            for k in 0..card.mode_count {
                row.push(direction[1] * row[k]);
            }
        }
        for value in &mut row[..card.mode_count] {
            *value *= direction[0];
        }
        row.resize(self.modes(), 0.0);
        if row.iter().any(|x| !x.is_finite()) || row.iter().all(|x| *x == 0.0) {
            return Err(invalid(
                "transverse string direction has no finite nonzero physical port",
            ));
        }
        Ok(row)
    }
}
