//! The SAME source-bound strings in the existing nonlinear contact mechanics.
//!
//! fs-nlmodal owns geometric stretching, fs-dcontact owns unilateral contact,
//! and fs-phs owns time integration. This module only compiles material inputs,
//! mass-normalized initial states and physical point maps. No assigned pitch
//! glide, contact force waveform, alternative integrator or sound gain is added.
//!
//! A frozen, moderate-slope Kirchhoff--Carrier model: axial strain equilibrates
//! instantaneously along the full fixed-end span. Both transverse planes of a
//! polarized source contribute to one shared elongation and tension. Optional torsion
//! remains a separate quadratic storage in that SAME joint mechanical solve.
//! This impact route supplies no axial wave inertia, evolving temperature,
//! plasticity, winding slip, changing topology or automatic transfer from a
//! running bow bank. The `longitudinal` module separately retains dynamic
//! material-derived axial coordinates in the playable friction bank.

use super::{MaterialStringSource, BowedStringCard, ModalAcousticMode,
    ModalAcousticTimeBudget, ModalAcousticTimeModel, GasState, Cx, CancelGate};
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::plate::impact::{BodyPotential, ImpactBody, ImpactError, MAX_IMPACT_MODES};
use crate::render::plate::impact::linear::wire::WireSpan;
use crate::render::plate::impact::string::StringStretching;

/// Geometric stretching in the original playable material bank and scheduler.
pub mod bank;

/// Dynamic material-derived axial strain in the original playable bank.
pub mod longitudinal;

fn invalid(what: &'static str) -> ImpactError { ImpactError::Invalid(what) }
fn poll(cx: &Cx<'_>, gate: &CancelGate) -> Result<(), ImpactError> {
    if gate.is_requested() { return Err(ImpactError::Cancelled); }
    cx.checkpoint().map_err(|_| ImpactError::Cancelled)
}

/// One immutable material realization, ready to join an existing ImpactSystem.
/// Bending is body zero; an explicitly requested torsional body follows. This
/// partition changes neither the coordinate ordering nor the total energy.
/// The original source remains attached; cloning bodies shares their storage.
#[derive(Clone, Debug)]
pub struct NonlinearMaterialString {
    source: MaterialStringSource,
    gas: GasState,
    card: BowedStringCard,
    wire: WireSpan,
    bodies: Vec<ImpactBody>,
    axial_rigidity_n: f64,
    sample_rate_hz: u32,
}
impl MaterialStringSource {
    fn axial_rigidity_n(&self) -> Result<f64, ImpactError> {
        let value = match self {
            Self::WithPolarizations(p) => return p.source().axial_rigidity_n(),
            Self::WithTorsion(t) => return t.bending().axial_rigidity_n(),
            Self::Homogeneous(s) => s.string().axial_stiffness_n,
            Self::Measured(_) => return Err(invalid("whole-string bending measurements do not determine axial rigidity for stretching")),
            Self::Wound(s) => s.string().axial_stiffness_n,
            Self::Concentric { section, .. } => section.elastic().axial_stiffness_n(),
            Self::ConcentricProny { section, .. } => section.elastic().axial_stiffness_n(),
            Self::ConcentricModalLoss { section, .. } => section.elastic().axial_stiffness_n(),
        };
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid("material string stretching requires its positive source-derived EA"));
        }
        Ok(value)
    }

    /// Derive linear bending, geometric stretching and optional twist together.
    ///
    /// Mass, T0, EA, EI, bending loss, shear loss and gas drag all come through
    /// the same source/section owners used by the bow bank. `maximum_slope`
    /// selects the original nonlinear owner's explicit moderate-slope envelope;
    /// it is not a force clamp or a replacement for the material's applicability.
    /// Nonlinear harmonics need temporal/modal refinement beyond a linear-mode
    /// Nyquist check. The material is a frozen linear-elastic constitutive image.
    ///
    /// `initial` contains every first-plane bending coordinate, then the complete
    /// second plane when selected, then optional twist. Both planes share EA,
    /// geometric elongation and one quartic storage; they are not separate strings.
    /// It specifies actual stored energy; no equilibrium or pluck is fabricated.
    /// Geometry is the full supported span, including afterlength. Interior
    /// fingers, frets, picks and bridges must be physical contacts/attachments.
    ///
    /// # Errors
    /// Missing source loss, inconsistent gas, invalid state/dimensions/slope,
    /// frequency-domain or clock refusal, unrepresentable data, cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn into_nonlinear(
        self, cx: &Cx<'_>, gas: GasState, sample_rate_hz: u32,
        maximum_slope: f64, initial: Vec<ModalAcousticState>, gate: &CancelGate,
    ) -> Result<NonlinearMaterialString, ImpactError> {
        poll(cx, gate)?;
        let count = self.modes();
        if count == 0 || count > MAX_IMPACT_MODES || initial.len() != count
            || initial.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite()
                || !s.velocity_m_sqrt_kg_per_s.is_finite()) {
            return Err(invalid("nonlinear material string needs bounded complete finite initial coordinates"));
        }
        let planes = if matches!(&self, Self::WithPolarizations(_)) { 2 } else { 1 };
        let axial_rigidity_n = self.axial_rigidity_n()?;
        let stretching = StringStretching { axial_rigidity_n, maximum_slope };
        stretching.validate()?;
        let card = self.lower(cx, &gas, sample_rate_hz).map_err(|e| ImpactError::Owner(e.to_string()))?;
        let wire = WireSpan {
            endpoints_m: [[0.0, 0.0], [card.length_m, 0.0]],
            linear_density_kg_m: card.linear_density_kg_m, tension_n: card.tension_n,
            bending_stiffness_n_m2: card.bending_stiffness_n_m2,
            damping_per_s: (0..card.mode_count)
                .map(|i| 2.0 * card.zetas[i] * card.mode_omega_rad_s(i + 1)).collect(),
        };
        // The original storage owns the quartic channel, its tangent, the exact
        // tension diagnostic and continuous-span slope admission at each step.
        let transverse_modes = card.mode_count.checked_mul(planes)
            .ok_or_else(|| invalid("nonlinear transverse coordinate count overflow"))?;
        if transverse_modes > count {
            return Err(invalid("nonlinear source has an inconsistent transverse prefix"));
        }
        let transverse = initial[..transverse_modes].to_vec();
        let body = if planes == 2 {
            wire.stretching_body_two_planes(transverse, stretching)?
        } else {
            wire.stretching_body(transverse, stretching)?
        };
        let mut bodies = vec![body];
        if let Some(twist) = self.torsion_source() {
            let mut modes: Vec<ModalAcousticMode> = Vec::with_capacity(twist.mode_count());
            twist.append_modes(&card, &mut modes, gate).map_err(|e| ImpactError::Owner(e.to_string()))?;
            // Reuse actual torsional source-band and modal-clock admission. No
            // fake twist frequency or borrowed transverse damping is supplied.
            ModalAcousticTimeModel::try_new(sample_rate_hz, modes.clone(), ModalAcousticTimeBudget::audible_reference())
                .map_err(|e| ImpactError::Owner(format!("torsional modal admission: {e:?}")))?;
            bodies.push(ImpactBody {
                potential: BodyPotential::Linear(modes.iter().map(|m| m.angular_frequency_rad_s).collect()),
                initial: initial[transverse_modes..].to_vec(),
                damping_per_s: modes.iter().map(|m| 2.0 * m.damping_ratio * m.angular_frequency_rad_s).collect(),
            });
        }
        poll(cx, gate)?;
        Ok(NonlinearMaterialString { source: self, gas, card, wire, bodies, axial_rigidity_n, sample_rate_hz })
    }
}
impl NonlinearMaterialString {
    /// Original material/section/rotary source, including all property receipts.
    #[must_use] pub const fn source(&self) -> &MaterialStringSource { &self.source }
    /// Same admitted medium that supplies the transverse drag.
    #[must_use] pub const fn gas(&self) -> &GasState { &self.gas }
    /// Complete original clock; the receiving time owner must use its reciprocal.
    #[must_use] pub const fn sample_rate_hz(&self) -> u32 { self.sample_rate_hz }
    /// Small-amplitude transverse dispersion and material/gas losses, not live tension.
    #[must_use] pub const fn card(&self) -> &BowedStringCard { &self.card }
    /// Actual span and damping produced by this same material realization.
    /// Consumers need not reconstruct EI, mass or material/gas loss from a preset.
    #[must_use] pub const fn wire_span(&self) -> &WireSpan { &self.wire }
    /// EA from the same original bulk/section data as mass and bending rigidity.
    #[must_use] pub const fn axial_rigidity_n(&self) -> f64 { self.axial_rigidity_n }
    /// Bending first, optional torsion second, for the SAME nonlinear time owner.
    #[must_use] pub fn bodies(&self) -> &[ImpactBody] { &self.bodies }
    /// Flattened coordinate count, including all selected torsional coordinates.
    #[must_use] pub fn mode_count(&self) -> usize { self.source.modes() }
    /// Complete transverse prefix; both planes belong to body zero.
    #[must_use] pub fn transverse_mode_count(&self) -> usize { self.bodies[0].initial.len() }

    /// Actual interior transverse displacement/force map [1/sqrt(kg)].
    /// Torsion has zero normal displacement to first order; it is not silently
    /// driven by a finger pressing radially on a circular string.
    pub fn transverse_shape(&self, station_fraction: f64) -> Result<Vec<f64>, ImpactError> {
        if !super::interior(station_fraction) {
            return Err(invalid("a fixed string endpoint is not a moving point contact"));
        }
        self.source.transverse_shape(&self.card, station_fraction)
            .map_err(|e| ImpactError::Owner(e.to_string()))
    }

    /// Centreline displacement/force along the original bridge direction.
    /// The same row and its transpose observe motion and apply reaction, across
    /// both transverse planes. Torsion has no centreline translation.
    pub fn bridge_shape(&self, station_fraction: f64) -> Result<Vec<f64>, ImpactError> {
        if !super::interior(station_fraction) {
            return Err(invalid("a fixed string endpoint is not a moving bridge port"));
        }
        self.source.bridge_shape(&self.card, station_fraction)
            .map_err(|e| ImpactError::Owner(e.to_string()))
    }

    /// Original material angular displacement/torque map, with a zero transverse
    /// prefix. Torque in N*m is work-conjugate to this angle in radians; no
    /// surface radius is inserted into an independently supplied torque.
    pub fn angular_shape(&self, station_fraction: f64) -> Result<Vec<f64>, ImpactError> {
        if !super::interior(station_fraction) {
            return Err(invalid("a fixed string endpoint is not a moving angular port"));
        }
        let twist = self.source.torsion_source()
            .ok_or_else(|| invalid("an angular port requires the original material torsion source"))?;
        twist.angular_shape(&self.card, station_fraction, self.transverse_mode_count())
            .map_err(|e| ImpactError::Owner(e.to_string()))
    }

    /// Actual tangential surface map, including r*theta for requested torsion.
    /// Its transpose applies F and rF with one work conjugacy. Choosing an
    /// external force here does not itself supply a bow/contact/friction law.
    pub fn surface_shape(&self, station_fraction: f64) -> Result<Vec<f64>, ImpactError> {
        if !super::interior(station_fraction) {
            return Err(invalid("a fixed string endpoint is not a moving surface port"));
        }
        self.source.surface_shape(&self.card, station_fraction)
            .map_err(|e| ImpactError::Owner(e.to_string()))
    }
}

/// Source-retaining preparation for the original shared mechanics and pressure path.
pub mod runtime;

/// Source-bound strings on the original inertially coupled translating support.
pub mod supported;
