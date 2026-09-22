//! Moving, finite-mass compliant contacts, compiled into the existing impact owner.
//!
//! Each jaw moves along one declared axis and carries a conforming felt pad.
//! Positive jaw travel closes its gap. Separate footprint sites see their own
//! surface displacement and material/creep history, but the SAME jaw inertia.
//! Contact is compression-only; neither release nor jaw motion is an audio fade.
//! fs-material owns felt and fs-phs owns time integration, including reactions,
//! energy loss, analytic tangents and rollback. This module adds neither law.
use super::{ImpactBody, ImpactError, MAX_IMPACT_MODES, felt::{FeltPad, KelvinBranch}, invalid};
use fs_material::fiber::WoolFelt;

/// A bounded footprint quadrature, not a continuum convergence certificate.
pub const MAX_PAD_SITES: usize = 4;
/// One moving pad or two opposed jaws in a single shared mechanical solve.
pub const MAX_MOVING_JAWS: usize = 2;

/// Side of the surface, relative to the axis used by every supplied shape row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PadSide {
    /// Jaw is on the positive side; positive surface travel closes this gap.
    Positive,
    /// Jaw is on the negative side; negative surface travel closes this gap.
    Negative,
}
impl PadSide {
    /// Signed surface contribution to pad compression.
    #[must_use]
    pub const fn sign(self) -> f64 {
        match self { Self::Positive => 1.0, Self::Negative => -1.0 }
    }
}

/// A positive-area site in a physical pad footprint.
#[derive(Debug, Clone)]
pub struct PadSite {
    /// Actual, unnormalized displacement participation along the declared axis
    /// [1/sqrt(kg)] in the complete original body basis. Unattached entries are zero.
    pub weights: Vec<f64>,
    /// Physical contact quadrature area [m^2]. Sites partition the loaded area.
    pub area_m2: f64,
}

/// Physical data for one translating jaw. No hand/skin calibration is inferred.
#[derive(Debug, Clone)]
pub struct CompliantJaw {
    /// Surface side; jaw travel and applied force are positive inward on either side.
    pub side: PadSide,
    /// Translating effective mass [kg], not a stiffness or frequency surrogate.
    pub mass_kg: f64,
    /// Grounded jaw drag [N s/m]; zero leaves a genuinely free inertia.
    pub drag_n_s_m: f64,
    /// Uncompressed gap at zero jaw/surface travel [m], nonnegative.
    pub initial_gap_m: f64,
    /// Inward launch velocity [m/s]; zero starts at rest, never auto-closes the pad.
    pub initial_velocity_m_s: f64,
    /// Uncompressed pad thickness [m]. Faces conform to the reference surface.
    pub thickness_m: f64,
    /// The existing material law, with caller-supplied physical coefficients.
    pub law: WoolFelt,
    /// Prior conditioning strain. Each site retains independent subsequent history.
    pub prior_maximum_strain: f64,
    /// Series Kelvin elements for the WHOLE jaw footprint. Each site's stiffness
    /// and viscosity are scaled by its area fraction, preserving physical totals.
    pub creep: Vec<KelvinBranch>,
}

/// Work-conjugate force/observation port of an appended jaw.
#[derive(Debug, Clone, Copy)]
pub struct JawPort {
    /// Generalized coordinate in the complete body basis, before private creep state.
    pub coordinate: usize,
    /// Multiply physical inward force by this to obtain generalized force;
    /// multiply q or p by this to observe physical inward travel or velocity.
    pub inverse_sqrt_mass: f64,
    /// Reference surface side; global axial jaw travel has the opposite sign.
    pub side: PadSide,
}

/// Cold-compiled attachment. Append its bodies and pads before constructing the
/// existing ImpactSystem; extend all old contact/pad/volume/drag rows with zeros.
/// The original coordinates never move. No mechanical state is reset mid-playback.
#[derive(Debug)]
pub struct MovingPads {
    /// One independent inertia per jaw, shared by all that jaw's footprint sites.
    pub bodies: Vec<ImpactBody>,
    /// Compression-only reciprocal contacts in the COMPLETE extended body basis.
    pub pads: Vec<FeltPad>,
    /// Actual physical inputs/observations; not microphone channels.
    pub ports: Vec<JawPort>,
}
impl MovingPads {
    /// Compile supplied geometry/materials. Bounds and every derived pad are
    /// checked here; the final impact owner still checks total assembly budgets.
    pub fn new(structural_modes: usize, sites: &[PadSite], jaws: &[CompliantJaw])
        -> Result<Self, ImpactError>
    {
        if structural_modes == 0 || sites.is_empty() || sites.len() > MAX_PAD_SITES
            || jaws.is_empty() || jaws.len() > MAX_MOVING_JAWS
            || structural_modes.checked_add(jaws.len()).is_none_or(|n| n > MAX_IMPACT_MODES)
        { return Err(invalid("moving pads exceed the declared sites, jaws or mechanical mode budget")); }
        let total = structural_modes + jaws.len();
        let mut area = 0.0;
        for site in sites {
            if !site.area_m2.is_finite() || site.area_m2 <= 0.0
                || site.weights.len() != structural_modes
                || site.weights.iter().any(|b| !b.is_finite())
                || site.weights.iter().all(|b| *b == 0.0)
            { return Err(invalid("pad site needs positive area and an actual finite moving-surface row")); }
            area += site.area_m2;
        }
        if !area.is_finite() { return Err(invalid("moving-pad total area is not representable")); }
        let mut result = Self { bodies: Vec::with_capacity(jaws.len()),
            pads: Vec::with_capacity(jaws.len()*sites.len()), ports: Vec::with_capacity(jaws.len()) };
        for (index, jaw) in jaws.iter().enumerate() {
            if !jaw.initial_gap_m.is_finite() || jaw.initial_gap_m < 0.0
                || !jaw.drag_n_s_m.is_finite() || jaw.drag_n_s_m < 0.0
                || jaws[..index].iter().any(|j| j.side == jaw.side)
            { return Err(invalid("moving jaws need distinct sides, finite gaps and nonnegative physical drag")); }
            let (mut body, weight) = ImpactBody::free_mass(jaw.mass_kg, 0.0, jaw.initial_velocity_m_s)?;
            let drag = jaw.drag_n_s_m / jaw.mass_kg;
            if !weight.is_finite() || !drag.is_finite()
                || !body.initial[0].velocity_m_sqrt_kg_per_s.is_finite()
                || (jaw.drag_n_s_m > 0.0 && drag == 0.0)
            { return Err(invalid("moving jaw mass/drag cannot be represented in the mechanical basis")); }
            body.damping_per_s[0] = drag;
            let coordinate = structural_modes + index;
            for site in sites {
                let fraction = site.area_m2/area;
                let mut weights: Vec<_> = site.weights.iter().map(|b| jaw.side.sign()*b).collect();
                weights.resize(total, 0.0); weights[coordinate] = weight;
                let pad = FeltPad { area_m2: site.area_m2, thickness_m: jaw.thickness_m,
                    precompression_m: -jaw.initial_gap_m, weights, law: jaw.law.clone(),
                    prior_maximum_strain: jaw.prior_maximum_strain,
                    creep: jaw.creep.iter().map(|b| KelvinBranch {
                        stiffness_n_m: fraction*b.stiffness_n_m,
                        viscosity_n_s_m: fraction*b.viscosity_n_s_m,
                    }).collect() };
                pad.validate(total)?;
                result.pads.push(pad);
            }
            result.bodies.push(body);
            result.ports.push(JawPort { coordinate, inverse_sqrt_mass: weight, side: jaw.side });
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
