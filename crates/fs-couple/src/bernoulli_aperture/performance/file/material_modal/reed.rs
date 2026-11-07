//! One typed, sourced lumped reed in place of the material plate of the DG-valve file.
//!
//! `reed HEIGHT WIDTH CLOSING_PA MASS STIFFNESS ZETA SOURCE_TOKEN` carries the six numbers of the `reed`
//! row of the reed-performance file (`render::schedule::reed`, lowered from a wind description by
//! `lower_reed_performance.py`) and the token of the description rows that sourced them, so a DG bore can
//! play the reed a modal-bore page plays and the bore is the only thing that changes.
//!
//! It is ONE structural mode of the existing `ModalAperture` (exact zero-order hold, uniform-slit contact,
//! Bernoulli slit flow, swept volume, the same mechanical limits as a plate): the gap is HEIGHT + q, the
//! pressure area is A = (k / P_c) H exactly as the live reed derives its face, and the damping is
//! c = 2 zeta sqrt(k m). Nothing is inferred that the row does not state: the reed-performance file's
//! zero-means-derived mass and stiffness are refused here, HEIGHT must equal the `valve` REST_GAP, the
//! contact is the file's own `contact` row (the lay of a lumped reed is the writer's to state), and no
//! viscous slit is selected. The row replaces `solid` or `loss`/`blades` and every plate record; it
//! carries no mesh, section, material memory, traction, closure or force port.

use super::{PlateValveInputError, Reader, bad};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use crate::bernoulli_aperture::paired::{BladePropagator, BladeState, PairedApertureLimits};
use fs_dcontact::Obstacle;

/// Longest source token, in bytes; the token names the sourcing rows, it is not prose.
const MAX_SOURCE_TOKEN_BYTES: usize = 256;

/// The authored reed exactly as the `reed` record states it.
#[derive(Clone, Debug, PartialEq)]
pub struct ReedRecord {
    rest_opening_m: f64,
    width_m: f64,
    closing_pressure_pa: f64,
    mass_kg: f64,
    stiffness_n_m: f64,
    damping_ratio: f64,
    source: String,
}

/// The host quantities the DG file supplies to every valve it builds.
pub(super) struct ReedHost {
    pub(super) density_kg_m3: f64,
    pub(super) impedance_pa_s_m3: f64,
    pub(super) time_step_s: f64,
    pub(super) max_steps: u64,
    pub(super) limits: PairedApertureLimits,
}

impl ReedRecord {
    /// Rest opening HEIGHT [m]; equal to the file's `valve` REST_GAP.
    pub const fn rest_opening_m(&self) -> f64 { self.rest_opening_m }
    /// Slit width [m].
    pub const fn width_m(&self) -> f64 { self.width_m }
    /// Static closing pressure [Pa] that, with the stiffness, fixes the pressure area.
    pub const fn closing_pressure_pa(&self) -> f64 { self.closing_pressure_pa }
    /// Moving mass [kg].
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Restoring stiffness [N/m].
    pub const fn stiffness_n_m(&self) -> f64 { self.stiffness_n_m }
    /// Viscous damping ratio.
    pub const fn damping_ratio(&self) -> f64 { self.damping_ratio }
    /// The token naming the description rows that sourced the six numbers.
    pub fn source(&self) -> &str { &self.source }
    /// Effective pressure area A = (k / P_c) H [m^2], associated as the live reed's own
    /// `reed_bore::reed_pressure_face` associates it, so both routes carry the same bits.
    pub fn pressure_area_m2(&self) -> f64 { (self.stiffness_n_m / self.closing_pressure_pa) * self.rest_opening_m }
    /// Viscous damping c = 2 zeta sqrt(k m) [N s/m].
    pub fn damping_n_s_m(&self) -> f64 {
        2.0 * self.damping_ratio * fs_math::det::sqrt(self.stiffness_n_m * self.mass_kg)
    }

    /// The reed as one `ModalAperture` mode on the DG host, at rest.
    pub(super) fn into_aperture(
        &self, host: ReedHost, contact: Obstacle, propagator: BladePropagator,
    ) -> Result<ModalAperture, AcousticRealizeError> {
        ModalAperture::new_with_propagator(
            ModalApertureSpec {
                rest_opening_m: self.rest_opening_m,
                width_m: self.width_m,
                modes: vec![ApertureMode {
                    mass_kg: self.mass_kg,
                    stiffness_n_m: self.stiffness_n_m,
                    damping_n_s_m: self.damping_n_s_m(),
                    pressure_area_m2: self.pressure_area_m2(),
                    gap_projection: 1.0,
                    force_projection: [0.0, 0.0],
                }],
                density_kg_m3: host.density_kg_m3,
                impedance_pa_s_m3: host.impedance_pa_s_m3,
                time_step_s: host.time_step_s,
                max_steps: host.max_steps,
                limits: host.limits,
            },
            vec![BladeState { displacement_m: 0.0, velocity_m_s: 0.0 }],
            contact,
            propagator,
        )
    }
}

/// Read and admit the `reed` record. `rest_opening` is the already read `valve` REST_GAP.
pub(super) fn read(rows: &mut Reader<'_>, rest_opening: f64) -> Result<ReedRecord, PlateValveInputError> {
    let mut row = rows.row("reed")?;
    let line = row.line;
    let rest_opening_m = row.scalar()?;
    let width_m = row.scalar()?;
    let closing_pressure_pa = row.scalar()?;
    let mass_kg = row.scalar()?;
    let stiffness_n_m = row.scalar()?;
    let damping_ratio = row.scalar()?;
    let source = row.word()?;
    row.finish()?;
    if !(rest_opening_m > 0.0 && width_m > 0.0 && closing_pressure_pa > 0.0
        && mass_kg > 0.0 && stiffness_n_m > 0.0)
    {
        return Err(bad(line,
            "a typed reed needs positive opening, width, closing pressure, mass and stiffness; none is derived"));
    }
    if damping_ratio < 0.0 {
        return Err(bad(line, "a typed reed needs a nonnegative damping ratio"));
    }
    if rest_opening_m != rest_opening {
        return Err(bad(line, "reed HEIGHT must equal the valve REST_GAP"));
    }
    if source.len() > MAX_SOURCE_TOKEN_BYTES || !source.bytes().all(|b| b.is_ascii_graphic()) {
        return Err(bad(line, "a typed reed needs a printable source token of at most 256 bytes"));
    }
    let record = ReedRecord {
        rest_opening_m, width_m, closing_pressure_pa, mass_kg, stiffness_n_m, damping_ratio,
        source: source.to_string(),
    };
    let area = record.pressure_area_m2();
    let damping = record.damping_n_s_m();
    if !(area.is_finite() && area > 0.0 && damping.is_finite()) {
        return Err(bad(line, "a typed reed's pressure area k H / P_c and damping must be finite"));
    }
    Ok(record)
}
