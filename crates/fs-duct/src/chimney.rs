//! Profiled tone-hole chimneys: a side branch whose chimney is a chain of
//! duct sections from the bore wall outward, for example an undercut cone
//! under a cylinder (IRP M5.5-W2 seam qualification).
//!
//! The model is one-dimensional. The chimney is the plane-wave chain of its
//! sections, terminated as [`tone_hole_shunt`](crate::tone_hole_shunt)
//! terminates a cylindrical one (flanged open, rigid closed, or the vent
//! admittance blend). The junction keeps the cylindrical hole's Dalmont inner
//! length and Nederveen series length, evaluated at the chimney's bore-side
//! radius. The evanescent near field that an undercut reshapes at the
//! junction is omitted; [`PROFILED_CHIMNEY_LIMITATION`] states this, and a
//! 3-D multiport tone hole is the reference for any real omission claim.
//! A single cylindrical section reproduces `tone_hole_shunt` exactly.

use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_phs::{WallPin, side_hole_inner_length, side_hole_series_length};

use crate::{
    Duct, DuctError, HoleState, LossModel, Segment, Termination, hole_sigma, input_impedance_wall,
    mul2,
};

/// Model-form caveat every comparison through a profiled chimney carries.
pub const PROFILED_CHIMNEY_LIMITATION: &str = "1-D profiled chimney: plane-wave sections with \
    cylindrical-hole junction lengths at the bore-side radius; the evanescent near field at the \
    junction is omitted";

fn chimney_radius(sections: &[Segment], bore_radius: f64) -> Result<f64, DuctError> {
    let bad = |what| DuctError::BadParameter { what };
    let first = sections
        .first()
        .ok_or(bad("a profiled chimney needs at least one section"))?;
    for s in sections {
        if matches!(s, Segment::ToneHole { .. }) {
            return Err(bad("chimney sections are cylinders, cones or profiles"));
        }
        s.validate()?;
    }
    let radius = first.inlet_radius();
    if !(bore_radius.is_finite() && radius < bore_radius) {
        return Err(bad(
            "the chimney's bore-side radius must be below the bore radius",
        ));
    }
    Ok(radius)
}

/// Shunt impedance of a tone hole whose chimney is `sections`, bore side
/// first, with the same [`WallPin`] the bore uses (`None` is rigid).
///
/// # Errors
/// [`DuctError`]: no sections, a tone-hole section, an invalid section, a
/// bore-side radius not below `bore_radius`, or the chimney chain's own
/// refusals.
pub fn profiled_tone_hole_shunt(
    state: &GasState,
    sections: &[Segment],
    hole_state: HoleState,
    omega: f64,
    loss: LossModel,
    bore_radius: f64,
    wall: Option<&WallPin>,
) -> Result<C64, DuctError> {
    let radius = chimney_radius(sections, bore_radius)?;
    let sigma = hole_sigma(hole_state);
    if sigma > 0.0 && sigma < 1.0 {
        let shunt =
            |s| profiled_tone_hole_shunt(state, sections, s, omega, loss, bore_radius, wall);
        let (z_o, z_c) = (shunt(HoleState::Open)?, shunt(HoleState::Closed)?);
        return Ok((z_o.recip().scale(sigma) + z_c.recip().scale(1.0 - sigma)).recip());
    }
    let chimney_loss = match loss {
        LossModel::WideTube => LossModel::Bessel,
        other => other,
    };
    let mut segments = sections.to_vec();
    let termination = if sigma <= 0.0 {
        Termination::Closed
    } else {
        // The inner length is extra tube at the bore-side radius: it extends
        // a cylindrical first section, otherwise it is prepended.
        let inner = side_hole_inner_length(radius, bore_radius);
        match segments[0] {
            Segment::Cylinder { radius, length } => {
                segments[0] = Segment::Cylinder {
                    radius,
                    length: length + inner,
                };
            }
            _ => segments.insert(
                0,
                Segment::Cylinder {
                    radius,
                    length: inner,
                },
            ),
        }
        Termination::FlangedOpen
    };
    let duct = Duct { segments };
    Ok(input_impedance_wall(&duct, state, omega, chimney_loss, termination, wall)?.impedance)
}

/// T-junction 2-port `[p_in, U_in] = J [p_out, U_out]` of a profiled tone
/// hole: series `Z_s/2`, shunt, series `Z_s/2`, with the Nederveen series
/// length of an open hole at the bore-side radius (a closed pad is a pure
/// shunt). There is no mutual term from other holes.
///
/// # Errors
/// As [`profiled_tone_hole_shunt`].
pub fn profiled_tone_hole_junction(
    state: &GasState,
    sections: &[Segment],
    hole_state: HoleState,
    omega: f64,
    loss: LossModel,
    bore_radius: f64,
    wall: Option<&WallPin>,
) -> Result<[C64; 4], DuctError> {
    let radius = chimney_radius(sections, bore_radius)?;
    let zh = profiled_tone_hole_shunt(state, sections, hole_state, omega, loss, bore_radius, wall)?;
    let shunt = [C64::ONE, C64::ZERO, zh.recip(), C64::ONE];
    let sigma = hole_sigma(hole_state);
    if sigma <= 0.0 {
        return Ok(shunt);
    }
    let ts = side_hole_series_length(radius, bore_radius) * sigma;
    if ts == 0.0 {
        return Ok(shunt);
    }
    let area = core::f64::consts::PI * bore_radius * bore_radius;
    let half = C64::new(0.0, -omega * state.density * ts / area).scale(0.5);
    let series = [C64::ONE, half, C64::ZERO, C64::ONE];
    Ok(mul2(mul2(series, shunt), series))
}
