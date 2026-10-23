//! Helical wraps on a tension-carrying core: added inertia and outer radius.
//!
//! A wrap of wire diameter d = 2a and pitch p lies on the diameter D beneath
//! it, so its centre-line helix has radius R = D/2 + a. One axial metre holds
//! `hypot(2 pi R, p) / p` metres of wire, and the wrap adds the wire's density
//! times its area times that length to the line density. The next wrap lies on
//! D + 2d. Neighbouring turns may not interpenetrate: their closest
//! centre-line approach, solved exactly for any helix angle, must be at least d.
//! A pitch of exactly d already overlaps slightly; [`close_wound_pitch_m`]
//! gives the tightest admissible pitch.
//!
//! This is the mass-only interim. The wraps carry none of the string tension
//! and add no axial or bending stiffness: the core section alone owns EA, EI,
//! its bending loss and the tensile check. Wrap stiffness, inter-turn contact,
//! slip friction, wrap loss and wrap tension from manufacture are absent. For a
//! single thin wrap, [`super::winding::CloseWoundSection`] adds the close-coiled
//! spring stiffness inside its admission domain; thick wire/coil ratios and
//! multiple wraps (piano bass strings) are outside it and use this owner.

use super::{ConcentricElasticSection, MaterialStatePointError, invalid, positive};
use crate::state_point::{DENSITY_PROPERTY, ResolvedMaterialStatePoint};
use fs_qty::{Density, QuantitySpec};
use std::f64::consts::{PI, TAU};

/// Hard cold-path bound on wrap layers over one core.
pub const MAX_HELICAL_WRAPS: usize = 8;

/// One helical round-wire wrap, innermost first.
#[derive(Clone, Debug, PartialEq)]
pub struct HelicalWrap {
    /// Round wire diameter [m].
    pub wire_diameter_m: f64,
    /// Axial advance per turn [m].
    pub pitch_m: f64,
    /// Exact source bundle containing the wire's dimensioned density.
    pub material: ResolvedMaterialStatePoint,
}

/// Wrap layers bound to one core, with their source bundles retained.
#[derive(Clone, Debug, PartialEq)]
pub struct HelicalWrapping {
    core_outer_radius_m: f64,
    wraps: Vec<HelicalWrap>,
    wrap_line_density_kg_m: Vec<f64>,
    line_density_kg_m: f64,
    outer_radius_m: f64,
}

/// Closest centre-line distance between different turns of a helix with
/// centre-line radius `radius` and axial `pitch`, over separations of at least
/// half a turn. Stationary separations solve `sin t + (c/R)^2 t = 0`, c = p/2pi.
fn turn_clearance(radius: f64, pitch: f64) -> f64 {
    let advance = pitch / TAU;
    let k = (advance / radius).powi(2);
    let distance = |turn: f64| (2.0 * radius * (0.5 * turn).sin()).hypot(advance * turn);
    let half_turn = distance(PI);
    if k >= 1.0 { return half_turn; }
    // f(t) = sin t + k t falls to its minimum at t* = pi + acos k and rises
    // monotonically on [t*, 2 pi]; the distance minimum is f's root there.
    let lowest = PI + k.acos();
    if lowest.sin() + k * lowest >= 0.0 { return half_turn; }
    let (mut lo, mut hi) = (lowest, TAU);
    for _ in 0..64 {
        let turn = 0.5 * (lo + hi);
        if turn.sin() + k * turn < 0.0 { lo = turn; } else { hi = turn; }
    }
    distance(0.5 * (lo + hi)).min(half_turn)
}

/// Tightest pitch [m] at which a wire of `wire_diameter_m` wound over a
/// cylinder of `under_diameter_m` does not interpenetrate its neighbouring turn.
///
/// # Errors
/// Nonpositive or nonfinite dimensions refuse.
pub fn close_wound_pitch_m(under_diameter_m: f64, wire_diameter_m: f64)
    -> Result<f64, MaterialStatePointError> {
    let under = positive(under_diameter_m, "positive finite wound-over diameter")?;
    let d = positive(wire_diameter_m, "positive finite wrap wire diameter")?;
    let radius = 0.5 * (under + d);
    // Clearance grows with pitch; it is below d at p = d and at least d at 2d.
    let (mut lo, mut hi) = (d, 2.0 * d);
    for _ in 0..64 {
        let pitch = 0.5 * (lo + hi);
        if turn_clearance(radius, pitch) < d { lo = pitch; } else { hi = pitch; }
    }
    Ok(hi)
}

impl HelicalWrapping {
    /// Wind `wraps` (innermost first) over the outer surface of `core`.
    /// Every wrap must resolve at the core's physical query coordinates.
    ///
    /// # Errors
    /// An empty or oversized wrap list, nonpositive or interpenetrating
    /// geometry, a different material state, missing/wrong-quantity/
    /// state-varying density, and unrepresentable arithmetic refuse atomically.
    pub fn try_new(
        core: &ConcentricElasticSection,
        wraps: Vec<HelicalWrap>,
    ) -> Result<Self, MaterialStatePointError> {
        if wraps.is_empty() || wraps.len() > MAX_HELICAL_WRAPS {
            return Err(invalid("helical wrapping requires 1..=8 wraps"));
        }
        let state = &core.layers()[0].material;
        let core_outer_radius_m = core.outer_radius_m();
        let mut diameter = 2.0 * core_outer_radius_m;
        let mut total = 0.0;
        let mut wrap_line_density_kg_m = Vec::with_capacity(wraps.len());
        for wrap in &wraps {
            let d = positive(wrap.wire_diameter_m, "positive finite wrap wire diameter")?;
            let pitch = positive(wrap.pitch_m, "positive finite wrap pitch")?;
            let radius = 0.5 * (diameter + d);
            if !(turn_clearance(radius, pitch) >= d) {
                return Err(invalid("helical wrap turns interpenetrate at the declared pitch"));
            }
            if wrap.material.query_point() != state.query_point()
                || wrap.material.axis_quantities() != state.axis_quantities() {
                return Err(invalid("core and wraps must resolve at identical physical coordinates"));
            }
            let density = wrap.material.property(DENSITY_PROPERTY)
                .ok_or_else(|| invalid("every wrap needs an unambiguous wire density"))?;
            if density.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS) {
                return Err(invalid("wrap density must be dimension-only kg/m3"));
            }
            if !density.is_constant_at_fixed_temperature() {
                return Err(invalid("wrap density must be a constant or T-only curve at fixed positive absolute T"));
            }
            let rho = positive(density.value_si(), "positive finite wrap density")?;
            let wire_per_axial_length = positive((TAU * radius).hypot(pitch) / pitch, "helix length ratio")?;
            let mass = positive(rho * (0.25 * PI * d * d) * wire_per_axial_length, "wrap line density")?;
            total = positive(total + mass, "summed wrap line density")?;
            wrap_line_density_kg_m.push(mass);
            diameter = positive(diameter + 2.0 * d, "wrapped diameter")?;
        }
        Ok(Self { core_outer_radius_m, wraps, wrap_line_density_kg_m,
            line_density_kg_m: total, outer_radius_m: 0.5 * diameter })
    }

    /// Outer radius [m] of the core these wraps were wound on.
    #[must_use]
    pub const fn core_outer_radius_m(&self) -> f64 { self.core_outer_radius_m }
    /// Original wrap declarations, innermost first.
    #[must_use]
    pub fn wraps(&self) -> &[HelicalWrap] { &self.wraps }
    /// Mass per axial length [kg/m] of each wrap, in declaration order.
    #[must_use]
    pub fn wrap_line_density_kg_m(&self) -> &[f64] { &self.wrap_line_density_kg_m }
    /// Mass per axial length [kg/m] of all wraps together, core excluded.
    #[must_use]
    pub const fn line_density_kg_m(&self) -> f64 { self.line_density_kg_m }
    /// Outer radius [m] of the outermost wrap.
    #[must_use]
    pub const fn outer_radius_m(&self) -> f64 { self.outer_radius_m }
}
