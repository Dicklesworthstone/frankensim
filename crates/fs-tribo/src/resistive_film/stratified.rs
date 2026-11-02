//! Exact pressure-flow reduction for a prescribed viscosity profile across a gap.
//!
//! Let z run from 0 to 1 across the instantaneous aperture, with stationary
//! tangential no-slip walls and continuous velocity/shear stress at interfaces.
//! Integrating d(mu du/dy)/dy = dp/dx gives the zero-shear plane z0 = B/A and
//! pressure mobility D = C - B*B/A, where A, B, C integrate 1/mu, z/mu, z*z/mu.
//! Thus Q = -width * h^3 * D * dp/dx and hydraulic-equivalent mu = 1/(12*D).
//! These are the variable-viscosity Reynolds integrals; see D. Dowson (1962),
//! "A generalized Reynolds equation for fluid-film lubrication",
//! <https://doi.org/10.1016/S0020-7403(62)80038-1>.
//!
//! This implementation evaluates D as a scaled weighted CENTRAL variance,
//! including each constant layer's exact thickness-squared/12 variance. It
//! never subtracts nearly equal raw second moments. The layer contribution
//! integral((z-z0)^2/mu dz) divided by D is its physical pressure-work fraction.
//!
//! The original pressure graph, analytic tangent and work accounting remain
//! the only flow/mechanical owner. Normalized interfaces and viscosities are
//! prescribed and held fixed within a solve. This is not a temperature solver,
//! layer-mass/advection model, mixing rule, slip/Couette law or moving-interface
//! conservation model. Axial profile variation needs a separate reduction.

use super::{FilmError, ResistiveFilm, fail};

/// Maximum constant-viscosity layers in one original channel.
/// A film has at most MAX_CHANNELS profiles; construction work is bounded.
pub const MAX_LAYERS: usize = 64;

/// One contiguous layer in original lower-wall to upper-wall order.
///
/// The lower bound is zero for the first layer and the preceding upper bound
/// thereafter. The final upper bound must be exactly one. Bounds are neither
/// sorted nor renormalized; gaps, overlaps and zero-thickness layers refuse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscosityLayer {
    /// Upper interface, normalized by the current hydraulic aperture, in (0, 1].
    pub upper_fraction: f64,
    /// Newtonian dynamic viscosity at zero gauge pressure [Pa s].
    pub dynamic_viscosity_pa_s: f64,
}

/// Retained layers and their exact piecewise-constant Reynolds reduction.
///
/// The hydraulic equivalent is a pressure-flow coefficient, not an arithmetic
/// material mean. Original layer values remain available without alteration.
#[derive(Debug, Clone, PartialEq)]
pub struct GapViscosityProfile {
    layers: Vec<ViscosityLayer>,
    hydraulic_equivalent_viscosity_pa_s: f64,
    neutral_plane_fraction: f64,
    layer_pressure_work_fractions: Vec<f64>,
}

impl GapViscosityProfile {
    /// Admit 1..MAX_LAYERS complete, ordered layers across the normalized gap.
    ///
    /// Positive finite viscosity and strictly increasing finite bounds are
    /// required. Computation is scaled by the least layer viscosity, then
    /// centered near the largest fluidity weight. Positive layer fluidity,
    /// variance contributions, work fractions and equivalent viscosity must be
    /// representable: underflow/overflow refuses rather than dropping a layer,
    /// flooring mobility or replacing a constitutive value.
    ///
    /// Equal viscosities recover that exact input viscosity and z0 = 0.5,
    /// independently of how the uniform profile was partitioned. Nonuniform
    /// profiles and work fractions are evaluated to floating-point roundoff.
    pub fn new(layers: Vec<ViscosityLayer>) -> Result<Self, FilmError> {
        if layers.is_empty() || layers.len() > MAX_LAYERS {
            return fail("stratified viscosity needs one to 64 complete gap layers");
        }
        let mut lower = [0.0; MAX_LAYERS];
        let mut thickness = [0.0; MAX_LAYERS];
        let mut minimum_viscosity = f64::INFINITY;
        let mut previous = 0.0;
        for (i, layer) in layers.iter().enumerate() {
            if !layer.upper_fraction.is_finite()
                || layer.upper_fraction <= previous
                || layer.upper_fraction > 1.0
                || !layer.dynamic_viscosity_pa_s.is_finite()
                || layer.dynamic_viscosity_pa_s <= 0.0
            {
                return fail("invalid ordered gap layer or positive viscosity");
            }
            lower[i] = previous;
            thickness[i] = layer.upper_fraction - previous;
            previous = layer.upper_fraction;
            minimum_viscosity = minimum_viscosity.min(layer.dynamic_viscosity_pa_s);
        }
        if previous != 1.0 {
            return fail("stratified viscosity must cover the complete normalized gap");
        }
        let uniform = layers
            .iter()
            .all(|layer| layer.dynamic_viscosity_pa_s == layers[0].dynamic_viscosity_pa_s);
        let mut weights = [0.0; MAX_LAYERS];
        let mut largest = 0;
        for (i, layer) in layers.iter().enumerate() {
            let relative_fluidity = minimum_viscosity / layer.dynamic_viscosity_pa_s;
            weights[i] = thickness[i] * relative_fluidity;
            if !weights[i].is_finite() || weights[i] <= 0.0 {
                return fail("scaled layer fluidity is not representable");
            }
            if weights[i] > weights[largest] {
                largest = i;
            }
        }
        let count = layers.len();
        let total_weight = compensated_sum(weights[..count].iter().copied());
        let anchor = lower[largest];
        let mut centers = [0.0; MAX_LAYERS];
        for i in 0..count {
            // Do not round a very thin layer's global midpoint before
            // centering: this preserves its local coordinate near either wall.
            centers[i] = (lower[i] - anchor) + 0.5 * thickness[i];
        }
        let mean_offset = if uniform {
            0.5 - anchor
        } else {
            compensated_sum((0..count).map(|i| weights[i] * centers[i])) / total_weight
        };
        let neutral_plane_fraction = if uniform { 0.5 } else { anchor + mean_offset };
        if !neutral_plane_fraction.is_finite() || !(0.0..=1.0).contains(&neutral_plane_fraction) {
            return fail("stratified neutral plane is not representable");
        }
        let mut contributions = [0.0; MAX_LAYERS];
        for i in 0..count {
            let distance = centers[i] - mean_offset;
            let within_layer = thickness[i] * thickness[i] / 12.0;
            contributions[i] = weights[i] * (distance * distance + within_layer);
            if !contributions[i].is_finite() || contributions[i] <= 0.0 {
                return fail("positive layer pressure mobility is not representable");
            }
        }
        let scaled_mobility = compensated_sum(contributions[..count].iter().copied());
        let hydraulic_equivalent_viscosity_pa_s = if uniform {
            layers[0].dynamic_viscosity_pa_s
        } else {
            minimum_viscosity / (12.0 * scaled_mobility)
        };
        if !hydraulic_equivalent_viscosity_pa_s.is_finite()
            || hydraulic_equivalent_viscosity_pa_s <= 0.0
        {
            return fail("stratified hydraulic-equivalent viscosity is not representable");
        }
        let mut layer_pressure_work_fractions = Vec::with_capacity(count);
        for contribution in &contributions[..count] {
            let fraction = contribution / scaled_mobility;
            if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
                return fail("stratified pressure-work fraction is not representable");
            }
            layer_pressure_work_fractions.push(fraction);
        }
        Ok(Self {
            layers,
            hydraulic_equivalent_viscosity_pa_s,
            neutral_plane_fraction,
            layer_pressure_work_fractions,
        })
    }

    /// Exact authored data, in original lower-wall to upper-wall order.
    #[must_use]
    pub fn layers(&self) -> &[ViscosityLayer] {
        &self.layers
    }

    /// The zero-gauge coefficient satisfying Q = -width*h^3/(12*mu)*dp/dx.
    /// This is a hydraulic equivalent, not a bulk or arithmetic material value.
    #[must_use]
    pub const fn hydraulic_equivalent_viscosity_pa_s(&self) -> f64 {
        self.hydraulic_equivalent_viscosity_pa_s
    }

    /// Normalized zero-shear location B/A, rounded to the global gap coordinate.
    /// The internal central variance retains local offsets before that rounding.
    #[must_use]
    pub const fn neutral_plane_fraction(&self) -> f64 {
        self.neutral_plane_fraction
    }

    /// Physical layer dissipation divided by total pressure-flow dissipation.
    ///
    /// Fractions are positive and sum to one to roundoff. Multiply by the
    /// original FilmChannelReport::dissipated_power_w, or its accepted integral,
    /// to locate pressure work. A closed/idle passage still generates zero work.
    /// This does not choose thermal destinations or a fluid/wall heat split.
    #[must_use]
    pub fn layer_pressure_work_fractions(&self) -> &[f64] {
        &self.layer_pressure_work_fractions
    }
}

/// One retained through-gap profile per original pressure-graph passage.
///
/// Only hydraulic equivalents are installed in the native ResistiveFilm.
/// Its ordinary flow, pressure, tangent and channel-work solver remain the
/// owner. Keep this binding to distinguish those equivalents from the actual
/// constituent viscosities; no layer averaging or second flow solver is hidden.
///
/// A common native Barus coefficient is retained. Within the original Reynolds
/// approximation p is uniform across the gap. Multiplying EVERY layer viscosity
/// by exp(alpha*p) multiplies all fluidities/mobility by exp(-alpha*p); the
/// neutral plane and work fractions stay fixed. The native transformed-pressure
/// solve and Q*delta_p channel work therefore remain valid. Different layer
/// pressure coefficients, changing profiles along a passage or non-Newtonian
/// layers require their own derivation and are not represented here.
#[derive(Debug, Clone)]
pub struct StratifiedFilm {
    profiles: Vec<GapViscosityProfile>,
    film: ResistiveFilm,
}

impl StratifiedFilm {
    /// Compile complete profiles in ORIGINAL channel order, including shut ones.
    /// Original geometry, force ports, pressure law and validity bounds survive.
    /// No admission failure mutates the supplied native operator.
    pub fn new(
        geometry: &ResistiveFilm,
        profiles: Vec<GapViscosityProfile>,
    ) -> Result<Self, FilmError> {
        if profiles.len() != geometry.channel_count() {
            return fail("stratified profiles must cover every original film channel");
        }
        let equivalent: Vec<_> = profiles
            .iter()
            .map(GapViscosityProfile::hydraulic_equivalent_viscosity_pa_s)
            .collect();
        let film = geometry.with_channel_viscosities(&equivalent)?;
        Ok(Self { profiles, film })
    }

    /// Exact retained profiles in the original channel insertion order.
    #[must_use]
    pub fn profiles(&self) -> &[GapViscosityProfile] {
        &self.profiles
    }

    /// Native flow/tangent/work owner with hydraulic-equivalent coefficients.
    /// Its viscosity getters describe these equivalents; layers() on each
    /// retained profile supplies the actual prescribed constituent values.
    #[must_use]
    pub const fn film(&self) -> &ResistiveFilm {
        &self.film
    }

    /// Release the compiled operator together with every original gap profile.
    pub fn into_parts(self) -> (ResistiveFilm, Vec<GapViscosityProfile>) {
        (self.film, self.profiles)
    }
}

fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum: f64 = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

#[cfg(test)]
mod tests;
