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
//! interlayer slip, unilateral bridge-foot contact, or large axial strain.
//! Optional cylinder gas loading supplies rotational drag in the same coupled
//! solve. All bows use the same declared positive surface side. The
//! material band screens the retained linear modes, not nonlinear harmonics.
//! See Bavu, Smith & Wolfe, 'Torsional waves in a bowed string' (2005),
//! arXiv:physics/0505059, for the importance of surface rather than centre speed.
//!
//! A measured whole string has no layers to derive torsion from. It can take
//! its own measured torsion instead ([`MeasuredStringTorsion`]): fundamental,
//! surface impedance and one Q, on the nut-to-bridge span. See that type.

use super::*;
use crate::shear_shell::constant_loss::ConstantLossSpectrum;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::torsion::{ConcentricTorsionalSection, ConcentricPronyTorsionalSection};
use fs_material::section::relaxation::LayerRelaxationProperties;
use fs_material::state_point::MaterialStatePointError;
use fs_qty::{Dims, QuantitySpec};

/// Optional layer property: the angular frequency [rad/s] at which a pack wrote its Kelvin-Voigt shear viscosity as
/// `eta_s = eta_t G / omega_ref` from a measured torsional loss factor `eta_t`. A layer that states it has the constant
/// loss stiffness `omega_ref eta_s`, and every twist mode gets that layer's `eta_t` at its own frequency, times the
/// shell's unit constant-loss spectrum ([`crate::shear_shell::constant_loss`]) fitted over
/// [`STRING_LOSS_SPECTRUM_BAND_RAD_S`]. A layer without it stays Kelvin-Voigt, bit for bit.
pub const SHEAR_VISCOSITY_REFERENCE_OMEGA: &str = "shear_viscosity_constant_loss_reference_omega";
/// Dimensions of [`SHEAR_VISCOSITY_REFERENCE_OMEGA`]: an angular frequency, 1/s.
pub const SHEAR_VISCOSITY_REFERENCE_OMEGA_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);
/// The span of that spectrum's branch rates and fit [rad/s], 50 Hz to 40 kHz: a bowed string keeps twist modes to
/// 30 kHz. A constant-loss twist mode outside it refuses.
pub const STRING_LOSS_SPECTRUM_BAND_RAD_S: (f64, f64) =
    (core::f64::consts::TAU * 50.0, core::f64::consts::TAU * 40_000.0);

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
        let Some(section) = twist.section() else {
            return Err(invalid("a rotational bridge requires material-derived string torsion; measured whole-string torsion has no polar inertia or radius"));
        };
        let radius = section.shear().elastic().outer_radius_m();
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
        let link_count = self.source.body_connections.len()
            .checked_add(self.source.translational_body_connections.len())
            .and_then(|count| count.checked_add(connections.len()))
            .ok_or_else(|| invalid("rotational bridge connection count overflow"))?;
        preflight(&self.source.strings, &bodies, &self.source.options, link_count)?;
        let bridges: Vec<_> = self.source.strings.iter().map(|s| s.connections.as_slice()).collect();
        let mut system = assemble_sources(&self.cards, Some(&self.source.strings), bodies,
            &bridges, &self.source.body_connections, &self.source.translational_body_connections,
            &connections, self.source.fluid_loading.as_ref(),
            self.source.prepared_linear, &self.source.constant_loss_bridges(),
            self.source.residual_flexibility, self.system.port().law.clone(),
            &self.source.options, gate)?;
        if let Some(weights) = self.system.acceleration_pressure_weights() {
            system = system.with_acceleration_pressure(weights.to_vec())
                .map_err(ReciprocalStringError::Friction)?;
        }
        if let Some(binding) = self.system.interface() {
            system = system.with_interface(binding.clone()).map_err(ReciprocalStringError::Friction)?;
        }
        system = self.source.restore_material_memory(system, &self.cards,
            Some(self.system.material_memory().to_vec()), gate)?;
        self.source.validate_material_bridge_frequencies(&system, gate)?;
        poll(gate)?;
        self.system = system;
        self.source.rotational_body_connections = connections;
        Ok(self)
    }
}

/// Torsion measured on one whole string: its torsional fundamental on the
/// speaking length, its characteristic torsional impedance at the surface,
/// and one quality factor. It belongs to a measured whole string
/// ([`MaterialStringSource::Measured`]); a material section derives its own.
///
/// The twist is a non-dispersive torsional line from the nut (station 0) to
/// the support station (the bridge). Rotation is fixed at both ends, and the
/// afterlength beyond the support carries no twist. In surface coordinates
/// `u = r theta`, the line inertia is `Z / c` [kg/m] and the stiffness `Z c`
/// [N], with `c = 2 L_s f_1`. The bow therefore sees exactly the measured
/// `Z = sqrt(K_R I_R) / r^2`, and mode `n` lies at `n f_1`. The source does not
/// state a radius, polar inertia or shear modulus, and none is inferred.
///
/// Every twist mode gets the measured Q, as `zeta = 1 / (2 Q)`. That follows a
/// measured Q that stays almost constant over the modes it was measured on; no
/// material law is used to vary it. No gas drag is added, because the measured
/// Q is the string's own loss in air. Consumers that need layers, polar inertia
/// or a radius (material heat, shear memory, rotational bridges, angular ports,
/// cylinder gas loading) refuse a measured torsion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeasuredStringTorsion {
    fundamental_hz: f64,
    impedance_kg_s: f64,
    quality_factor: f64,
    support_fraction: f64,
}

impl MeasuredStringTorsion {
    /// Admit the measured torsional fundamental [Hz] of the span from the nut to
    /// `support_fraction` of the full string, the surface impedance [kg/s], and
    /// the quality factor of every twist mode.
    ///
    /// # Errors
    /// A nonpositive or nonfinite fundamental or impedance, a quality factor that
    /// is not finite and above 1/2, or a support station outside (0, 1).
    pub fn try_new(fundamental_hz: f64, impedance_kg_s: f64, quality_factor: f64, support_fraction: f64)
        -> Result<Self, MaterialStatePointError>
    {
        if !fundamental_hz.is_finite() || fundamental_hz <= 0.0
            || !impedance_kg_s.is_finite() || impedance_kg_s <= 0.0
            || !quality_factor.is_finite() || quality_factor <= 0.5
            || !support_fraction.is_finite() || support_fraction <= 0.0 || support_fraction >= 1.0 {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "measured torsion needs a positive fundamental and impedance, a quality factor above 1/2 and an interior support station",
            });
        }
        Ok(Self { fundamental_hz, impedance_kg_s, quality_factor, support_fraction })
    }
    /// Measured torsional fundamental of the twisting span [Hz].
    #[must_use]
    pub const fn fundamental_hz(&self) -> f64 { self.fundamental_hz }
    /// Measured characteristic torsional impedance at the surface, `sqrt(K_R I_R) / r^2` [kg/s].
    #[must_use]
    pub const fn impedance_kg_s(&self) -> f64 { self.impedance_kg_s }
    /// Measured quality factor, the same for every twist mode.
    #[must_use]
    pub const fn quality_factor(&self) -> f64 { self.quality_factor }
    /// Station of the twist's far end (the bridge) as a fraction of the full span.
    #[must_use]
    pub const fn support_fraction(&self) -> f64 { self.support_fraction }
    /// Twisting span [m], nut to support, on a string whose full span is `length_m`.
    #[must_use]
    pub fn span_m(&self, length_m: f64) -> f64 { length_m * self.support_fraction }
    /// Torsional wave speed `c = 2 L_s f_1` [m/s].
    #[must_use]
    pub fn wave_speed_m_s(&self, length_m: f64) -> f64 { 2.0 * self.span_m(length_m) * self.fundamental_hz }
    /// Line inertia in surface coordinates, `I_R / r^2 = Z / c` [kg/m].
    #[must_use]
    pub fn surface_inertia_kg_m(&self, length_m: f64) -> f64 {
        self.impedance_kg_s / self.wave_speed_m_s(length_m)
    }
}

/// Where a twist source's coefficients come from.
#[derive(Clone, Debug)]
enum TwistLaw {
    /// The original material section's shear stiffness, polar inertia and shear loss.
    Section(ConcentricTorsionalSection),
    /// A measurement of the whole string's own torsion.
    Measured(MeasuredStringTorsion),
}

/// Frozen original bending specimen plus torsion from its SAME material states,
/// or from the same measured string's own torsion measurement.
/// Private fields prevent attaching an unrelated section or nesting reductions.
#[derive(Clone, Debug)]
pub struct TorsionalStringSource {
    bending: MaterialStringSource,
    law: TwistLaw,
    relaxation: Option<ConcentricPronyTorsionalSection>,
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
    /// Optional cylinder gas loading adds rotational drag using the same gas
    /// and section radius, while preserving the original material polar inertia.
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
        self.with_torsion_binding(Some(shear_property), viscosity_property, modes)
    }

    /// Add torsion using G from each layer's original explicit isotropic pair.
    /// E/nu, E/G and K/G share the material section's existing projection; no
    /// separate shear-modulus claim, material query or nominal value is added.
    /// The named shear viscosity remains an independent source requirement.
    ///
    /// # Errors
    /// Requires an explicitly parameterized concentric source and the same
    /// geometry, mode-count, loss and frequency admission as `with_torsion`.
    pub fn with_isotropic_torsion(
        self,
        viscosity_property: &str,
        modes: usize,
    ) -> Result<Self, MaterialStatePointError> {
        if let Self::WithPolarizations(p) = self {
            return p.add_isotropic_torsion(viscosity_property, modes);
        }
        self.with_torsion_binding(None, viscosity_property, modes)
    }

    fn with_torsion_binding(
        self,
        shear_property: Option<&str>,
        viscosity_property: &str,
        modes: usize,
    ) -> Result<Self, MaterialStatePointError> {
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
            Self::ConcentricProny { section, .. } => section.elastic().clone(),
            Self::ConcentricModalLoss { section, .. } => section.elastic().clone(),
            Self::Wound(_) => return Err(invalid("close-coiled winding torsion requires a resolved core/wrap model")),
            Self::Measured(_) => return Err(invalid("whole-string bending measurements do not determine torsional stiffness or polar inertia")),
            Self::WithTorsion(_) => return Err(invalid("nested torsion is not a new physical section")),
            Self::WithPolarizations(_) => return Err(invalid("polarized torsion must retain its original source")),
        };
        // The correction multiplies transverse GA only; it never changes GJ.
        let shear = match shear_property {
            Some(property) => elastic.with_shear_modulus(property, 1.0)?,
            None => elastic.with_isotropic_shear(1.0)?,
        };
        let section = shear.with_torsional_viscosity(viscosity_property)?;
        Ok(Self::WithTorsion(Box::new(TorsionalStringSource {
            bending: self,
            law: TwistLaw::Section(section),
            relaxation: None,
            mode_count: modes,
            total_modes,
        })))
    }

    /// Add a measured whole string's own measured torsion as `modes` twist modes
    /// on its nut-to-support span ([`MeasuredStringTorsion`]).
    ///
    /// Only [`Self::Measured`] takes it, with or without its two transverse
    /// planes, and either polarization order gives the same source. Bending
    /// coordinates come first, then the twist modes `1..=modes` at `n f_1`.
    ///
    /// # Errors
    /// A source that is not a measured whole string, already twisted, or a mode
    /// selection outside 1..=4096. Clock admission occurs at assembly.
    pub fn with_measured_torsion(self, torsion: MeasuredStringTorsion, modes: usize)
        -> Result<Self, MaterialStatePointError>
    {
        if let Self::WithPolarizations(p) = self {
            return p.add_measured_torsion(torsion, modes);
        }
        let invalid = |quantity| MaterialStatePointError::InvalidDerived { quantity };
        if modes == 0 || modes > 4096 || matches!(&self, Self::WithTorsion(_)) {
            return Err(invalid("torsional extension requires 1..=4096 modes and a non-nested source"));
        }
        if !matches!(&self, Self::Measured(_)) {
            return Err(invalid("measured torsion belongs to a measured whole string; a material section derives its own"));
        }
        let total_modes = self.modes().checked_add(modes)
            .ok_or_else(|| invalid("combined bending/torsion mode count overflow"))?;
        Ok(Self::WithTorsion(Box::new(TorsionalStringSource {
            bending: self, law: TwistLaw::Measured(torsion), relaxation: None, mode_count: modes, total_modes,
        })))
    }

    /// Add independently selected per-layer shear memory to original torsion.
    ///
    /// Call after `with_torsion`; either transverse-polarization order is valid.
    /// All source moduli, times and equilibrium roles must already be resolved.
    /// Kelvin shear viscosity remains a separately sourced parallel dashpot.
    /// The bank automatically projects these histories onto its actual twist
    /// modes, in the same friction/bridge/body solve and energy account.
    /// Source bending histories precede torsional histories; regional body and
    /// explicitly added layered bending histories follow both.
    ///
    /// # Errors
    /// Absent torsion, a second shear-memory selection, or shared material-law
    /// admission failure. Other consumers refuse if they cannot retain history.
    pub fn with_torsional_relaxation(
        self, properties: &[LayerRelaxationProperties],
    ) -> Result<Self, MaterialStatePointError> {
        match self {
            Self::WithPolarizations(p) => p.add_torsional_relaxation(properties),
            Self::WithTorsion(mut twist) => {
                if twist.relaxation.is_some() {
                    return Err(MaterialStatePointError::InvalidDerived {
                        quantity: "torsional relaxation requires one original source binding",
                    });
                }
                let TwistLaw::Section(section) = &twist.law else {
                    return Err(MaterialStatePointError::InvalidDerived {
                        quantity: "shear relaxation requires a material torsional section; measured torsion has no layers",
                    });
                };
                twist.relaxation = Some(section.clone().with_prony_torsion(properties)?);
                Ok(Self::WithTorsion(twist))
            }
            Self::Homogeneous(_) | Self::Measured(_) | Self::Wound(_) | Self::Concentric { .. }
                | Self::ConcentricProny { .. } | Self::ConcentricModalLoss { .. } =>
                Err(MaterialStatePointError::InvalidDerived {
                    quantity: "shear relaxation requires the original torsional source",
                }),
        }
    }

    pub(super) fn surface_shape(&self, card: &BowedStringCard, station: f64)
        -> Result<Vec<f64>, ReciprocalStringError>
    {
        let mut shapes = self.transverse_shape(card, station)?;
        if let Some(twist) = self.torsion_source() {
            let offset = shapes.len() - twist.mode_count;
            match &twist.law {
                TwistLaw::Section(section) => {
                    let inertia = section.rotary_inertia_kg_m();
                    let radius = section.shear().elastic().outer_radius_m();
                    for k in 0..twist.mode_count {
                        let value = radius * unit_shape(k, station, card.length_m, inertia);
                        if !value.is_finite() { return Err(invalid("torsional surface shape is not representable")); }
                        shapes[offset + k] = value;
                    }
                }
                TwistLaw::Measured(measured) => {
                    // Surface displacement per modal coordinate on the nut-to-support span;
                    // the afterlength does not twist.
                    let support = measured.support_fraction();
                    let span = measured.span_m(card.length_m);
                    let inertia = measured.surface_inertia_kg_m(card.length_m);
                    for k in 0..twist.mode_count {
                        let value = if station < support { unit_shape(k, station / support, span, inertia) } else { 0.0 };
                        if !value.is_finite() { return Err(invalid("torsional surface shape is not representable")); }
                        shapes[offset + k] = value;
                    }
                }
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
    /// `None` for a measured whole string's measured torsion, which has no section.
    #[must_use]
    pub const fn section(&self) -> Option<&ConcentricTorsionalSection> {
        match &self.law {
            TwistLaw::Section(section) => Some(section),
            TwistLaw::Measured(_) => None,
        }
    }
    /// The measured whole string's own torsion, when that is the source.
    #[must_use]
    pub const fn measured(&self) -> Option<&MeasuredStringTorsion> {
        match &self.law {
            TwistLaw::Section(_) => None,
            TwistLaw::Measured(measured) => Some(measured),
        }
    }
    /// Optional independent shear spectra, retained on the original source.
    #[must_use]
    pub const fn relaxation(&self) -> Option<&ConcentricPronyTorsionalSection> {
        self.relaxation.as_ref()
    }
    /// Number of appended twist modes; bending coordinates always come first.
    #[must_use]
    pub const fn mode_count(&self) -> usize { self.mode_count }
    pub(super) const fn total_modes(&self) -> usize { self.total_modes }

    pub(super) fn angular_shape(&self, card: &BowedStringCard, station: f64, transverse_modes: usize)
        -> Result<Vec<f64>, ReciprocalStringError>
    {
        let TwistLaw::Section(section) = &self.law else {
            return Err(invalid("an angular map requires material torsion; measured whole-string torsion has no polar inertia"));
        };
        let inertia = section.rotary_inertia_kg_m();
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
        let section = match &self.law {
            TwistLaw::Section(section) => section,
            TwistLaw::Measured(measured) => return self.append_measured_modes(measured, card, modes, gate),
        };
        let inertia = section.rotary_inertia_kg_m();
        let speed = (section.stiffness_n_m2() / inertia).sqrt();
        let norm = (inertia * card.length_m * 0.5).sqrt();
        let radius = section.shear().elastic().outer_radius_m();
        if !speed.is_finite() || speed <= 0.0 || !norm.is_finite() || norm <= 0.0
            || !(radius / norm).is_finite() || radius / norm == 0.0 {
            return Err(invalid("torsion frequency and angular normalization must be representable"));
        }
        let constant_loss = Self::constant_loss(section)?;
        let (lo, hi) = match &constant_loss {
            None => section.omega_band_rad_s(),
            Some(law) => law.band,
        };
        for n in 1..=self.mode_count {
            poll(gate)?;
            let k = n as f64 * core::f64::consts::PI / card.length_m;
            let omega = k * speed;
            let viscous = match &constant_loss {
                None => section.viscous_stiffness_n_m2_s(),
                Some(law) => law.viscous_stiffness_at(omega),
            };
            let zeta = (viscous / inertia) * k * k / (2.0 * omega);
            if !omega.is_finite() || omega <= 0.0 || omega < lo || omega > hi
                || !zeta.is_finite() || zeta < 0.0
                || (section.viscous_stiffness_n_m2_s() > 0.0 && zeta == 0.0) {
                return Err(invalid("actual torsional modes exceed material support or coefficient representability"));
            }
            modes.push(ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: zeta,
                pressure_per_modal_velocity: fs_math::c64::C64::ZERO });
        }
        Ok(())
    }

    /// Measured twist modes: `omega_n = 2 pi n f_1` and the measured `zeta = 1 / (2 Q)` on every mode.
    fn append_measured_modes(&self, measured: &MeasuredStringTorsion, card: &BowedStringCard,
        modes: &mut Vec<ModalAcousticMode>, gate: &CancelGate) -> Result<(), ReciprocalStringError>
    {
        let span = measured.span_m(card.length_m);
        let inertia = measured.surface_inertia_kg_m(card.length_m);
        let norm = (inertia * span * 0.5).sqrt();
        let zeta = 0.5 / measured.quality_factor();
        if !span.is_finite() || span <= 0.0 || !inertia.is_finite() || inertia <= 0.0
            || !norm.is_finite() || norm <= 0.0 || !(1.0 / norm).is_finite()
            || !zeta.is_finite() || zeta <= 0.0 {
            return Err(invalid("measured torsion's span, surface inertia and loss must be representable"));
        }
        for n in 1..=self.mode_count {
            poll(gate)?;
            let omega = core::f64::consts::TAU * measured.fundamental_hz() * n as f64;
            if !omega.is_finite() || omega <= 0.0 {
                return Err(invalid("measured torsional mode frequency is not representable"));
            }
            modes.push(ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: zeta,
                pressure_per_modal_velocity: fs_math::c64::C64::ZERO });
        }
        Ok(())
    }

    /// The per-layer shear loss when any layer states [`SHEAR_VISCOSITY_REFERENCE_OMEGA`], else `None`. Its band is
    /// the section's material band within the spectrum's fit band and each stated anchor's own validity.
    fn constant_loss(section: &ConcentricTorsionalSection) -> Result<Option<TorsionalConstantLoss>, ReciprocalStringError> {
        let elastic = section.shear().elastic();
        if elastic.layers().iter().all(|layer| layer.material.property(SHEAR_VISCOSITY_REFERENCE_OMEGA).is_none()) {
            return Ok(None);
        }
        let (lo, hi) = section.omega_band_rad_s();
        let mut band = (lo.max(STRING_LOSS_SPECTRUM_BAND_RAD_S.0), hi.min(STRING_LOSS_SPECTRUM_BAND_RAD_S.1));
        let mut layers = Vec::with_capacity(elastic.layers().len());
        for (layer, geometry) in elastic.layers().iter().zip(elastic.layer_integrals()) {
            // The same per-layer term as the section's sum: eta_s times the polar moment J = 2 I.
            let eta = layer.material.property(section.viscosity_property())
                .ok_or_else(|| invalid("constant-loss torsion lost a layer's shear viscosity"))?.value_si();
            let viscous = eta * (2.0 * geometry.second_moment_m4());
            let reference = match layer.material.property(SHEAR_VISCOSITY_REFERENCE_OMEGA) {
                None => None,
                Some(p) => {
                    let w = p.value_si();
                    if p.requirement().quantity() != QuantitySpec::dimensional(SHEAR_VISCOSITY_REFERENCE_OMEGA_DIMS)
                        || !w.is_finite() || w <= 0.0 || !p.is_constant_at_fixed_temperature() {
                        return Err(invalid("a shear-viscosity reference angular frequency is a positive constant in rad/s"));
                    }
                    if let Some((a, b)) = p.angular_frequency_band_rad_s()
                        .map_err(|_| invalid("invalid shear-viscosity reference frequency semantics"))? {
                        band = (band.0.max(a), band.1.min(b));
                    }
                    Some(w)
                }
            };
            layers.push((viscous, reference));
        }
        if band.0.is_nan() || band.1.is_nan() || band.0 > band.1 {
            return Err(invalid("constant-loss torsion has no frequency band inside the loss spectrum's fit band"));
        }
        Ok(Some(TorsionalConstantLoss { layers, spectrum: ConstantLossSpectrum::fit_over(STRING_LOSS_SPECTRUM_BAND_RAD_S),
            band }))
    }
}

/// A twist source's shear loss with at least one constant-loss layer.
struct TorsionalConstantLoss {
    /// Per layer: `eta_s J` [N m2 s] and the layer's reference angular frequency, `None` for a Kelvin-Voigt layer.
    layers: Vec<(f64, Option<f64>)>,
    spectrum: ConstantLossSpectrum,
    /// Where the law is claimed [rad/s]; a twist mode outside it refuses.
    band: (f64, f64),
}

impl TorsionalConstantLoss {
    /// The viscous stiffness that gives a twist mode at `omega` its loss: a constant-loss layer's `eta_s J` times
    /// `omega_ref G(omega) / omega` (loss factor `eta_t G(omega)` at every frequency), a Kelvin-Voigt layer's unchanged.
    fn viscous_stiffness_at(&self, omega: f64) -> f64 {
        let shape = self.spectrum.shape(omega);
        self.layers.iter().map(|&(viscous, reference)| reference.map_or(viscous, |w| viscous * w * shape / omega)).sum()
    }
}
