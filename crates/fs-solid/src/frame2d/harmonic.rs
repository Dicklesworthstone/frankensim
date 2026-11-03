//! Forced resonance of the original linear frame, without a modal truncation.
//!
//! The native reference mass, elastic stiffness and complete material viscosity
//! enter `K - omega^2 M - i omega C` under `exp(-i omega t)`. Authored Rayleigh
//! damping is optional through explicit zero coefficients; it is not inferred
//! from a material name or substituted for the original section viscosity.
//! Stationary supports use the frame's original free-coordinate map.

use super::linear::{FrameLinearError, FrameLinearPencil};
use super::{Frame2d, Rayleigh};
use crate::SolidError;
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicPortControls, HarmonicResponse, HarmonicResponseControls, HarmonicResponseError,
    PreparedHarmonicResponse,
};

/// Original frame/source refusal or shared forced-response refusal.
#[derive(Debug)]
pub enum FrameHarmonicError {
    Frame(SolidError),
    Harmonic(HarmonicResponseError),
}
impl core::fmt::Display for FrameHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-FRAME-HARMONIC: {self:?}")
    }
}
impl std::error::Error for FrameHarmonicError {}
impl From<SolidError> for FrameHarmonicError {
    fn from(error: SolidError) -> Self {
        match error {
            SolidError::Cancelled => Self::Harmonic(HarmonicResponseError::Cancelled),
            error => Self::Frame(error),
        }
    }
}
impl From<HarmonicResponseError> for FrameHarmonicError {
    fn from(error: HarmonicResponseError) -> Self {
        Self::Harmonic(error)
    }
}

impl From<FrameLinearError> for FrameHarmonicError {
    fn from(error: FrameLinearError) -> Self {
        match error {
            FrameLinearError::Frame(error) => Self::Frame(error),
            FrameLinearError::Invalid(what) => HarmonicResponseError::Invalid(what).into(),
            FrameLinearError::Budget {
                quantity,
                requested,
                cap,
            } => HarmonicResponseError::Budget {
                quantity,
                requested,
                cap,
            }
            .into(),
            FrameLinearError::Cancelled => HarmonicResponseError::Cancelled.into(),
        }
    }
}

fn invalid(what: &'static str) -> FrameHarmonicError {
    HarmonicResponseError::Invalid(what).into()
}
fn poll(cx: &Cx<'_>) -> Result<(), FrameHarmonicError> {
    cx.checkpoint()
        .map_err(|_| HarmonicResponseError::Cancelled.into())
}

/// A physical effort/flow load evaluated at the requested frequency.
///
/// The real row-major force map has `frame.ndof()` rows and `ports` columns,
/// in the ORIGINAL node-major `[ux, uy, theta]` order, including fixed rows.
/// A column maps one load effort into nodal forces and moments; its transpose
/// maps nodal velocities into the conjugate port flow. Lever arms and surface
/// quadrature therefore belong in this same map. No bare-frame solve is needed
/// to specify a load at an otherwise singular resonance.
///
/// `impedance` is the complete row-major `ports * ports` matrix Z, mapping
/// flow to opposing effort under `exp(-i omega t)`. For example a resistance,
/// attached inertance and spring have Z = r - i*omega*m + i*k/omega at omega>0.
/// Mutual and reactive entries are retained. Fixed rows have zero velocity;
/// their load efforts can be returned, but support reactions are not computed.
///
/// The load's owner must supply physical port units, frequency-band admission
/// and broadband causality/passivity. Shared pointwise Hermitian-PSD admission
/// cannot establish those properties. Finite Z at DC cannot encode a spring's
/// 1/omega pole; put static boundary stiffness in the frame or query omega>0.
#[derive(Clone, Copy, Debug)]
pub struct FrameHarmonicPortLoad<'a> {
    /// Full native-coordinate force map, never a truncated modal image.
    pub force_map: &'a [f64],
    /// Number of physical ports, including any motionless supported ports.
    pub ports: usize,
    /// Original frequency-dependent complex impedance in conjugate port units.
    pub impedance: &'a [C64],
    /// Separate load budgets; payload includes reduction of the full force map
    /// plus the shared physical-load admission, projection and response arrays.
    pub controls: HarmonicPortControls,
}

impl FrameHarmonicPortLoad<'_> {
    fn reduced_controls(&self, n: usize) -> Result<HarmonicPortControls, FrameHarmonicError> {
        if self.ports == 0 {
            return Err(invalid("frame load requires at least one physical port"));
        }
        if self.ports > self.controls.maximum_ports {
            return Err(HarmonicResponseError::Budget {
                quantity: "frame load ports",
                requested: self.ports,
                cap: self.controls.maximum_ports,
            }
            .into());
        }
        let entries = n
            .checked_mul(self.ports)
            .ok_or_else(|| invalid("full frame load map extent overflow"))?;
        if self.force_map.len() != entries
            || self.ports.checked_mul(self.ports) != Some(self.impedance.len())
            || self.force_map.iter().any(|x| !x.is_finite())
            || self
                .impedance
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid(
                "finite full-frame force map and complete port impedance required",
            ));
        }
        // Reserve the full map size before assembly; the retained map is never
        // larger. Shared admission budgets all its own additional arrays.
        let bytes = entries
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or_else(|| invalid("frame load reduction payload overflow"))?;
        if bytes > self.controls.maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "frame load reduction bytes",
                requested: bytes,
                cap: self.controls.maximum_payload_bytes,
            }
            .into());
        }
        Ok(HarmonicPortControls {
            maximum_payload_bytes: self.controls.maximum_payload_bytes - bytes,
            ..self.controls
        })
    }
}

/// One immutable frequency factor, retaining the original frame and supports.
/// Repeated force/receiver queries reuse the same shared native factorization.
pub struct FrameHarmonicOperator<'a> {
    source: &'a Frame2d,
    free_dofs: Vec<usize>,
    response: PreparedHarmonicResponse,
}

impl FrameHarmonicOperator<'_> {
    /// Source geometry, material bindings, supports and accepted state, unchanged.
    pub fn source(&self) -> &Frame2d {
        self.source
    }

    /// Original node-major `[ux, uy, theta]` indices of every retained coordinate.
    /// Translation is in metres and rotation in radians; conjugate forces are
    /// newtons and newton-metres respectively. No massive coordinate is removed.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Solve for peak applied force phasors in `free_dofs()` order. Returned
    /// displacement and velocity use that same map and `exp(-i omega t)`.
    ///
    /// The power account includes the complete original material viscosity AND
    /// the explicitly supplied Rayleigh model once. It is dissipated mechanical
    /// power, not an automatic heat deposit. Fixed coordinates have zero motion;
    /// loads applied there affect reactions, which this response does not report.
    /// This particular periodic solution does not erase an existing transient.
    pub fn solve(
        &self,
        forces: &[C64],
        cx: &Cx<'_>,
    ) -> Result<HarmonicResponse, FrameHarmonicError> {
        Ok(self.response.solve(forces, cx)?)
    }
}

impl Frame2d {
    /// Prepare the full original frame's forced response at one frequency.
    ///
    /// Only fixed linear geometry and elastic sections are admitted. Every free
    /// coordinate requires positive reference inertia: use consistent member
    /// mass or explicit section spin rather than inventing rotational mass or
    /// statically condensing a viscous first-order coordinate. Source viscosity
    /// must admit the actual frequency, including DC when requested. Its full
    /// off-diagonal matrix survives; no fitted modal damping is introduced.
    /// Original rigid modes remain intact at positive forcing frequency. This
    /// query reuses source viscosity assembly without a transient Nyquist test
    /// or an unnecessary natural-mode eigensolve.
    /// Only roundoff-sized native stiffness skew is symmetrized in scratch;
    /// larger asymmetry refuses and the source frame is never rewritten.
    ///
    /// `maximum_coordinates` bounds the FULL frame, including fixed coordinates,
    /// before dense assembly. The payload ceiling covers a conservative envelope
    /// of numerical arrays for native assembly/admission and the shared solve;
    /// source material metadata, caller inputs and allocator overhead are excluded.
    /// Assembly polls per member; existing bounded mass/material preparations are
    /// indivisible calls between checkpoints. No accepted frame state is mutated.
    ///
    /// # Errors
    /// Unsupported geometry/section/mass, source-frequency refusal, resource
    /// limits, nonfinite coefficients, nonpassive operators, singular lossless
    /// poles, cancellation or the original shared solve's numerical refusal.
    pub fn prepare_harmonic_response(
        &self,
        omega_rad_s: f64,
        rayleigh: Rayleigh,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<FrameHarmonicOperator<'_>, FrameHarmonicError> {
        self.prepare_harmonic(omega_rad_s, rayleigh, controls, None, cx)
    }

    /// Prepare the complete native frame with a physical boundary load.
    ///
    /// Assemble `K - omega^2 M - i omega C - i omega B Z B^T` BEFORE
    /// factorization. A physical load can resolve a bare-frame pole and move
    /// the coupled resonance. All native material viscosity, authored Rayleigh
    /// damping and mutual load impedance remain in their original owners.
    /// Returned `port_load` contains the original flows, opposing efforts and
    /// net absorbed boundary work; that work is not automatically material heat
    /// or acoustic radiation. The same geometry, inertia, source-frequency and
    /// cancellation conditions as `prepare_harmonic_response` apply.
    pub fn prepare_harmonic_response_with_port_impedance(
        &self,
        omega_rad_s: f64,
        rayleigh: Rayleigh,
        load: FrameHarmonicPortLoad<'_>,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<FrameHarmonicOperator<'_>, FrameHarmonicError> {
        self.prepare_harmonic(omega_rad_s, rayleigh, controls, Some(load), cx)
    }

    fn prepare_harmonic(
        &self,
        omega_rad_s: f64,
        rayleigh: Rayleigh,
        controls: HarmonicResponseControls,
        load: Option<FrameHarmonicPortLoad<'_>>,
        cx: &Cx<'_>,
    ) -> Result<FrameHarmonicOperator<'_>, FrameHarmonicError> {
        poll(cx)?;
        if !omega_rad_s.is_finite()
            || omega_rad_s < 0.0
            || !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(invalid(
                "finite nonnegative frequency and residual budget in (0,1) required",
            ));
        }
        let n = self.ndof();
        let load = load
            .map(|load| load.reduced_controls(n).map(|controls| (load, controls)))
            .transpose()?;
        let FrameLinearPencil {
            free_dofs,
            mass: m,
            stiffness: k,
            damping: c,
            material_omega_band_rad_s,
        } = self.prepare_linear_pencil(
            rayleigh,
            controls.maximum_coordinates,
            controls.maximum_payload_bytes,
            cx,
        )?;
        if let Some((low, high)) = material_omega_band_rad_s {
            if omega_rad_s < low || omega_rad_s > high {
                return Err(invalid(
                    "harmonic frequency lies outside original frame viscosity support",
                ));
            }
        }
        let nf = free_dofs.len();
        let response = if let Some((load, port_controls)) = load {
            let mut force_map = Vec::with_capacity(nf * load.ports);
            for &row in &free_dofs {
                poll(cx)?;
                force_map
                    .extend_from_slice(&load.force_map[row * load.ports..(row + 1) * load.ports]);
            }
            PreparedHarmonicResponse::try_with_port_impedance(
                &m,
                &c,
                &k,
                nf,
                &force_map,
                load.ports,
                load.impedance,
                omega_rad_s,
                controls,
                port_controls,
                cx,
            )?
        } else {
            PreparedHarmonicResponse::try_new(&m, &c, &k, nf, omega_rad_s, controls, cx)?
        };
        poll(cx)?;
        Ok(FrameHarmonicOperator {
            source: self,
            free_dofs,
            response,
        })
    }
}
