//! Twist and transverse motion at the SAME material-string surface contact.
//!
//! For rotary line inertia Ip and circular stiffness GJ, fixed-rotation ends
//! give omega_n = (n*pi/L)*sqrt(GJ/Ip). Shear loss gives 2*zeta*omega =
//! (eta_s J/Ip)*(n*pi/L)^2. The original modal owner advances these coordinates.
//! Bow slip sees v + r*theta_dot; its transpose applies BOTH F and r*F. The
//! An explicitly supplied rotational bridge additionally exchanges torque with
//! body rotation through that same mechanical solve, including after bow release.
//! No acoustic gain is assigned to twist, and no extra friction work or loss is
//! counted separately.
//!
//! This is small-angle twist of bonded circular layers, not helical winding,
//! interlayer slip, unilateral bridge-foot contact, rotary gas drag, or large
//! axial strain. All bows use the same declared positive surface side. The
//! material band screens the retained linear modes, not nonlinear harmonics.
//! See Bavu, Smith & Wolfe, 'Torsional waves in a bowed string' (2005),
//! arXiv:physics/0505059, for the importance of surface rather than centre speed.

use super::*;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::torsion::ConcentricTorsionalSection;
use fs_material::state_point::MaterialStatePointError;

/// Angular displacement/torque participation in the supplied body's SAME basis.
/// Unlike a translational `ModalAttachment`, entries have units
/// rad / (m sqrt(kg)); multiplying the row by torque gives generalized force.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalRotationAttachment {
    /// Component in the original body list, before material strings are added.
    pub component: usize,
    /// Small-angle rotation about one declared signed axis, one entry per mode.
    pub angular_shapes: Vec<f64>,
}

impl ModalRotationAttachment {
    /// Derive rotation from two physical displacement footprints and their
    /// positive separation [m], without independently authored modal constants.
    ///
    /// The returned row is `(positive - negative) / separation`. It is the
    /// small-angle chord rotation for that baseline, exact for a rigid rocking
    /// member. For a deforming plate it is a secant approximation, not a local
    /// slope certificate; the footprint/baseline needs spatial refinement.
    /// Both footprints must belong to the SAME modal reduction and displacement
    /// direction. Equal dimensions alone do not establish that provenance.
    /// A zero footprint can explicitly represent an immobile supported point.
    ///
    /// # Errors
    /// Different bodies/dimensions, invalid separation, nonfinite/overflowing
    /// maps, or a baseline with no retained angular participation.
    pub fn from_displacement_footprints(
        positive: &ModalAttachment,
        negative: &ModalAttachment,
        separation_m: f64,
    ) -> Result<Self, ReciprocalStringError> {
        if positive.component != negative.component || positive.shapes.is_empty()
            || positive.shapes.len() != negative.shapes.len()
            || !separation_m.is_finite() || separation_m <= 0.0
            || positive.shapes.iter().chain(&negative.shapes).any(|x| !x.is_finite()) {
            return Err(invalid("rotation needs two same-body footprints and a finite positive separation"));
        }
        let angular_shapes: Vec<_> = positive.shapes.iter().zip(&negative.shapes)
            .map(|(p, n)| (p - n) / separation_m).collect();
        if angular_shapes.iter().any(|x| !x.is_finite())
            || angular_shapes.iter().all(|x| *x == 0.0)
            || positive.shapes.iter().zip(&negative.shapes).zip(&angular_shapes)
                .any(|((p, n), angle)| p != n && *angle == 0.0) {
            return Err(invalid("body rotation is detached or not representable in the retained basis"));
        }
        Ok(Self { component: positive.component, angular_shapes })
    }
}

/// A retained translational bridge footprint in the original source bank.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StringBridgeLocation {
    /// Original material-string index.
    pub string: usize,
    /// Index within that member's `connections` list.
    pub connection: usize,
}

impl MaterialStringBankSource {
    /// Form a body rotation from two actually retained bridge footprints.
    /// Material plate/chart producers already projected both rows onto ONE
    /// body basis, so a material swap changes the rotation map automatically.
    /// The same small-angle secant assumptions as `from_displacement_footprints`
    /// apply. This operation neither adds bridge mass nor infers its geometry.
    pub fn rotation_between_bridges(
        &self,
        positive: StringBridgeLocation,
        negative: StringBridgeLocation,
        separation_m: f64,
    ) -> Result<ModalRotationAttachment, ReciprocalStringError> {
        let lookup = |at: StringBridgeLocation| {
            self.strings.get(at.string).and_then(|s| s.connections.get(at.connection))
                .map(|c| &c.body)
                .ok_or_else(|| invalid("body rotation names an unknown material string bridge"))
        };
        ModalRotationAttachment::from_displacement_footprints(
            lookup(positive)?, lookup(negative)?, separation_m,
        )
    }
}

/// A reciprocal small-angle spring/dashpot between string twist and body rotation.
/// The relative coordinate is `theta_string - theta_body - rest_angle_rad`.
/// Both angles must use the same physical signed axis. Law parameters retain
/// caller authority; bulk string properties do not determine a bridge joint.
#[derive(Clone, Debug, PartialEq)]
pub struct StringBodyRotationConnection {
    /// Original material-string index; its source must explicitly include torsion.
    pub string: usize,
    /// Interior station on the full supported string, including afterlength.
    pub station_fraction: f64,
    /// Actual angular participation of the receiving material/body reduction.
    pub body: ModalRotationAttachment,
    /// Nonnegative rotational stiffness [N m / rad].
    pub stiffness_n_m_per_rad: f64,
    /// Nonnegative rotational resistance [N m s / rad].
    pub damping_n_m_s_per_rad: f64,
    /// Signed stress-free relative rotation [rad].
    pub rest_angle_rad: f64,
}

impl StringBodyRotationConnection {
    pub(super) fn lower(
        &self,
        cards: &[BowedStringCard],
        members: &[MaterialStringMember],
        bodies: &[ModalAcousticTimeModel],
    ) -> Result<ModalConnection, ReciprocalStringError> {
        let member = members.get(self.string)
            .ok_or_else(|| invalid("rotational connection names an unknown material string"))?;
        let twist = member.source.torsion_source()
            .ok_or_else(|| invalid("a rotational connection requires material-derived string torsion"))?;
        let card = cards.get(self.string)
            .ok_or_else(|| invalid("rotational connection has no aligned string realization"))?;
        let body = bodies.get(self.body.component)
            .ok_or_else(|| invalid("rotational connection names an unknown original body"))?;
        if !interior(self.station_fraction)
            || self.body.angular_shapes.len() != body.modes().len()
            || self.body.angular_shapes.iter().any(|x| !x.is_finite())
            || self.body.angular_shapes.iter().all(|x| *x == 0.0)
            || !self.stiffness_n_m_per_rad.is_finite() || self.stiffness_n_m_per_rad < 0.0
            || !self.damping_n_m_s_per_rad.is_finite() || self.damping_n_m_s_per_rad < 0.0
            || (self.stiffness_n_m_per_rad == 0.0 && self.damping_n_m_s_per_rad == 0.0)
            || !self.rest_angle_rad.is_finite() {
            return Err(invalid("rotational bridge requires an interior point, a finite angular body map and a nonzero passive law"));
        }
        // The original owner uses metres/Newtons. This exact coordinate change
        // preserves energy and work: x=r*theta, k_x=k_theta/r², F_x=tau/r.
        // Use the actual material section radius, never an arbitrary modal gain.
        let radius = twist.section.shear().elastic().outer_radius_m();
        let stiffness = self.stiffness_n_m_per_rad / radius / radius;
        let damping = self.damping_n_m_s_per_rad / radius / radius;
        let rest = self.rest_angle_rad * radius;
        let right: Vec<_> = self.body.angular_shapes.iter().map(|x| radius * x).collect();
        if !stiffness.is_finite() || !damping.is_finite() || !rest.is_finite()
            || (self.stiffness_n_m_per_rad != 0.0 && stiffness == 0.0)
            || (self.damping_n_m_s_per_rad != 0.0 && damping == 0.0)
            || (self.rest_angle_rad != 0.0 && rest == 0.0)
            || right.iter().zip(&self.body.angular_shapes)
                .any(|(scaled, original)| !scaled.is_finite() || (*original != 0.0 && *scaled == 0.0)) {
            return Err(invalid("rotational bridge's work-conjugate metre/Newton lowering is not representable"));
        }
        let transverse_modes = member.source.modes() - twist.mode_count;
        let mut left = twist.angular_shape(card, self.station_fraction, transverse_modes)?;
        for value in &mut left {
            let scaled = *value * radius;
            if !scaled.is_finite() || (*value != 0.0 && scaled == 0.0) {
                return Err(invalid("rotational bridge's string surface map is not representable"));
            }
            *value = scaled;
        }
        if left.iter().any(|x| !x.is_finite()) || left[transverse_modes..].iter().all(|x| *x == 0.0) {
            return Err(invalid("rotational bridge's string surface map is not representable"));
        }
        Ok(ModalConnection {
            left: ModalAttachment { component: self.string, shapes: left },
            right: ModalAttachment { component: cards.len() + self.body.component, shapes: right },
            stiffness_n_m: stiffness, damping_n_s_m: damping, rest_extension_m: rest,
        })
    }
}

impl ReciprocalStringBank {
    /// Add a complete set of reciprocal torque/rotation laws before playback.
    ///
    /// Material bending/torsion, translational bridges, shared body/cavity and
    /// material damping remain in ONE existing exact linear/friction solve.
    /// The original pressure weights and dry-interface binding survive. Body
    /// modal initial states are retained; string states must still be at rest.
    /// All source identities and rotational inputs survive multibow scheduling.
    ///
    /// This prepares the complete network again, without redoing any material
    /// queries or body eigenproblem. It cannot modify an ongoing performance.
    /// An empty list returns the original bank unchanged. A second nonempty
    /// binding refuses; supply all rotational links together.
    ///
    /// The connection owner's force ceiling bounds torque divided by the
    /// sourced outer string radius (the corresponding tangential surface force).
    /// No contact-separation, large-angle, bridge-inertia, radiation-loading or
    /// measured-instrument claim follows from an admitted rotational spring.
    ///
    /// # Errors
    /// Missing torsional/body sources, invalid angular maps or law/geometry,
    /// already-started/bound bank, complete-network limits, or cancellation.
    pub fn with_rotational_body_connections(
        mut self,
        connections: Vec<StringBodyRotationConnection>,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if connections.is_empty() { return Ok(self); }
        if self.system.has_nonlinear_potentials() {
            return Err(invalid("bind rotational bridges before geometric stretching; a rebuild must not discard its storage"));
        }
        if self.system.samples_rendered() != 0 || !self.source.rotational_body_connections.is_empty()
            || self.system.components()[..self.cards.len()].iter().flat_map(|c| c.states())
                .any(|s| s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
            return Err(invalid("rotational bridges require one complete binding before string motion or playback"));
        }
        let body_end = self.source.fluid_loading.as_ref().map_or(
            self.system.components().len(), |fluid| self.cards.len() + fluid.body_components,
        );
        let bodies = self.system.components()[self.cards.len()..body_end].to_vec();
        let link_count = self.source.body_connections.len().checked_add(connections.len())
            .ok_or_else(|| invalid("rotational bridge connection count overflow"))?;
        preflight(&self.source.strings, &bodies, &self.source.options, link_count)?;
        let bridges: Vec<_> = self.source.strings.iter().map(|s| s.connections.as_slice()).collect();
        let mut system = assemble_sources(&self.cards, Some(&self.source.strings), bodies,
            &bridges, &self.source.body_connections, &connections, self.source.fluid_loading.as_ref(),
            self.source.prepared_linear,
            self.system.port().law.clone(),
            &self.source.options, gate)?;
        if let Some(weights) = self.system.acceleration_pressure_weights() {
            system = system.with_acceleration_pressure(weights.to_vec())
                .map_err(ReciprocalStringError::Friction)?;
        }
        if let Some(binding) = self.system.interface() {
            system = system.with_interface(binding.clone()).map_err(ReciprocalStringError::Friction)?;
        }
        self.source.validate_material_bridge_frequencies(&system, gate)?;
        poll(gate)?;
        self.system = system;
        self.source.rotational_body_connections = connections;
        Ok(self)
    }
}

/// Frozen original bending specimen plus torsion from its SAME material states.
/// Private fields prevent attaching an unrelated section or nesting reductions.
#[derive(Clone, Debug)]
pub struct TorsionalStringSource {
    bending: MaterialStringSource,
    section: ConcentricTorsionalSection,
    mode_count: usize,
    total_modes: usize,
}
impl MaterialStringSource {
    /// Add small-angle torsion using named G [Pa] and shear viscosity [Pa s].
    ///
    /// The original homogeneous or concentric material bundles supply all new
    /// coefficients; there is no additional material query or independent mass.
    /// Source values must already be present. The original bending model and
    /// prestress are retained, not retuned to match a torsional resonance.
    /// All torsional endpoints are fixed in rotation over the complete span.
    /// Gas continues to load transverse motion; rotational gas loading is absent.
    ///
    /// # Errors
    /// Zero/oversized/nested mode selection, unrepresentable total dimension, or
    /// original source-law refusal. Frequency/clock admission occurs at assembly.
    pub fn with_torsion(self, shear_property: &str, viscosity_property: &str, modes: usize)
        -> Result<Self, MaterialStatePointError>
    {
        if let Self::WithPolarizations(p) = self {
            return p.add_torsion(shear_property, viscosity_property, modes);
        }
        let invalid = |quantity| MaterialStatePointError::InvalidDerived { quantity };
        if modes == 0 || modes > 4096 || matches!(&self, Self::WithTorsion(_)) {
            return Err(invalid("torsional extension requires 1..=4096 modes and a non-nested source"));
        }
        let total_modes = self.modes().checked_add(modes)
            .ok_or_else(|| invalid("combined bending/torsion mode count overflow"))?;
        let elastic = match &self {
            Self::Homogeneous(s) => ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
                outer_radius_m: s.radius_m(), material: s.material().clone(),
            }])?,
            Self::Concentric { section, .. } => section.elastic().clone(),
            Self::WithTorsion(_) => return Err(invalid("nested torsion is not a new physical section")),
            Self::WithPolarizations(_) => return Err(invalid("polarized torsion must retain its original source")),
        };
        // The correction multiplies transverse GA only; it never changes GJ.
        let section = elastic.with_shear_modulus(shear_property, 1.0)?
            .with_torsional_viscosity(viscosity_property)?;
        Ok(Self::WithTorsion(Box::new(TorsionalStringSource {
            bending: self, section, mode_count: modes, total_modes,
        })))
    }

    pub(super) fn surface_shape(&self, card: &BowedStringCard, station: f64)
        -> Result<Vec<f64>, ReciprocalStringError>
    {
        let mut shapes = self.transverse_shape(card, station)?;
        if let Some(twist) = self.torsion_source() {
            let inertia = twist.section.rotary_inertia_kg_m();
            let radius = twist.section.shear().elastic().outer_radius_m();
            let offset = shapes.len() - twist.mode_count;
            for k in 0..twist.mode_count {
                let value = radius * unit_shape(k, station, card.length_m, inertia);
                if !value.is_finite() { return Err(invalid("torsional surface shape is not representable")); }
                shapes[offset + k] = value;
            }
        }
        Ok(shapes)
    }
}
impl TorsionalStringSource {
    /// Unmodified original specimen/layers, bending loss and support prescription.
    #[must_use]
    pub const fn bending(&self) -> &MaterialStringSource { &self.bending }
    /// Original material section with independent shear loss and rotary inertia.
    #[must_use]
    pub const fn section(&self) -> &ConcentricTorsionalSection { &self.section }
    /// Number of appended twist modes; bending coordinates always come first.
    #[must_use]
    pub const fn mode_count(&self) -> usize { self.mode_count }
    pub(super) const fn total_modes(&self) -> usize { self.total_modes }

    pub(super) fn angular_shape(&self, card: &BowedStringCard, station: f64, transverse_modes: usize)
        -> Result<Vec<f64>, ReciprocalStringError>
    {
        let inertia = self.section.rotary_inertia_kg_m();
        let mut row = vec![0.0; transverse_modes];
        for k in 0..self.mode_count {
            let value = unit_shape(k, station, card.length_m, inertia);
            if !value.is_finite() { return Err(invalid("torsional angular shape is not representable")); }
            row.push(value);
        }
        Ok(row)
    }

    pub(super) fn append_modes(&self, card: &BowedStringCard, modes: &mut Vec<ModalAcousticMode>, gate: &CancelGate)
        -> Result<(), ReciprocalStringError>
    {
        let inertia = self.section.rotary_inertia_kg_m();
        let speed = (self.section.stiffness_n_m2() / inertia).sqrt();
        let norm = (inertia * card.length_m * 0.5).sqrt();
        let radius = self.section.shear().elastic().outer_radius_m();
        if !speed.is_finite() || speed <= 0.0 || !norm.is_finite() || norm <= 0.0
            || !(radius / norm).is_finite() || radius / norm == 0.0 {
            return Err(invalid("torsion frequency and angular normalization must be representable"));
        }
        let (lo, hi) = self.section.omega_band_rad_s();
        for n in 1..=self.mode_count {
            poll(gate)?;
            let k = n as f64 * core::f64::consts::PI / card.length_m;
            let omega = k * speed;
            let zeta = (self.section.viscous_stiffness_n_m2_s() / inertia) * k * k / (2.0 * omega);
            if !omega.is_finite() || omega <= 0.0 || omega < lo || omega > hi
                || !zeta.is_finite() || zeta < 0.0
                || (self.section.viscous_stiffness_n_m2_s() > 0.0 && zeta == 0.0) {
                return Err(invalid("actual torsional modes exceed material support or coefficient representability"));
            }
            modes.push(ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: zeta,
                pressure_per_modal_velocity: fs_math::c64::C64::ZERO });
        }
        Ok(())
    }
}
