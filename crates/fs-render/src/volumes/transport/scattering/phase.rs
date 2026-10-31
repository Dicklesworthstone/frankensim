//! Reciprocal phase functions using the propagation-direction convention.

use super::{VolumeTransportError, unit_variate};
use fs_blake3::DomainHasher;
use fs_geom::Vec3;
use fs_math::det;

mod table;
pub use table::{MAX_PHASE_KNOTS, TabulatedPhase};
mod spectral;
pub use spectral::{
    MAX_PHASE_WAVELENGTHS, PhasePacketEvaluation, PhasePacketSample, SpectralPhaseTable,
};

const INV_FOUR_PI: f64 = 1.0 / (4.0 * core::f64::consts::PI);

/// Normalized phase probability per solid angle. Both directions describe
/// forward path travel: positive HG asymmetry scatters along the incoming ray,
/// unlike a surface convention with both vectors pointing away from a vertex.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseFunction {
    /// Uniform direction, probability `1 / (4 pi)`.
    Isotropic,
    /// HG anisotropy `-1 < g < 1`; its mean scattering cosine is exactly `g`.
    HenyeyGreenstein {
        /// Signed mean forward cosine, independent of wavelength.
        g: f64,
    },
    /// Unpolarized Rayleigh law `3 (1 + cos²(theta)) / (16 pi)`.
    Rayleigh,
    /// Molecular Rayleigh law with an explicit unpolarized 90-degree
    /// depolarization ratio, `I_parallel / I_perpendicular`.
    RayleighDepolarized {
        /// Finite ratio in `[0, 6/7)`; zero is the ordinary Rayleigh law.
        depolarization: f64,
    },
    /// Reciprocal, azimuthally symmetric tabulated law with a matched sampler.
    Tabulated(TabulatedPhase),
    /// Convex interpolation of normalized angular laws at vacuum wavelengths.
    /// Use the wavelength/packet methods; gray evaluation refuses this law.
    SpectralTable(SpectralPhaseTable),
}

/// An importance-sampled unit travel direction and its solid-angle density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseSample {
    /// Unit outgoing direction in world space.
    pub direction: Vec3,
    /// Matched phase probability per steradian.
    pub pdf: f64,
}

impl PhaseFunction {
    /// Construct the HG model with explicit finite, nonsingular asymmetry.
    pub fn henyey_greenstein(g: f64) -> Result<Self, VolumeTransportError> {
        let result = Self::HenyeyGreenstein { g };
        result.validate()?;
        Ok(result)
    }

    /// Normalize an angular table and construct its matched phase sampler.
    pub fn tabulated(knots: &[(f64, f64)]) -> Result<Self, VolumeTransportError> {
        TabulatedPhase::try_new(knots).map(Self::Tabulated)
    }

    /// Molecular depolarization, using the same ratio as the King correction
    /// `(6 + 3 delta)/(6 - 7 delta)` to the total scattering cross section.
    /// This scalar law describes initially unpolarized light, not Stokes transport.
    pub fn rayleigh_depolarized(depolarization: f64) -> Result<Self, VolumeTransportError> {
        let result = Self::RayleighDepolarized { depolarization };
        result.validate()?;
        Ok(if depolarization == 0.0 {
            Self::Rayleigh
        } else {
            result
        })
    }

    pub(super) fn validate(&self) -> Result<(), VolumeTransportError> {
        if let Self::HenyeyGreenstein { g } = self
            && (!g.is_finite() || g.abs() >= 1.0)
        {
            return Err(VolumeTransportError::Phase);
        }
        if let Self::RayleighDepolarized { depolarization } = self
            && (!depolarization.is_finite() || !(0.0..6.0 / 7.0).contains(depolarization))
        {
            return Err(VolumeTransportError::Phase);
        }
        Ok(())
    }

    pub(super) fn hash_into(&self, hash: &mut DomainHasher) {
        match self {
            Self::Isotropic => hash.update(&[0]),
            Self::HenyeyGreenstein { g } => {
                hash.update(&[1]);
                hash.update(&g.to_bits().to_le_bytes());
            }
            Self::Rayleigh => hash.update(&[2]),
            Self::Tabulated(table) => {
                hash.update(&[3]);
                table.hash_into(hash);
            }
            Self::RayleighDepolarized { depolarization } => {
                hash.update(&[4]);
                hash.update(&depolarization.to_bits().to_le_bytes());
            }
            Self::SpectralTable(table) => {
                hash.update(&[5]);
                table.hash_into(hash);
            }
        }
    }

    /// Evaluate using `incoming_travel.dot(outgoing_travel)`; the result is
    /// both the physical phase value and the sampling PDF. Tiny dot-product
    /// roundoff outside the unit interval clamps, while invalid angles refuse.
    pub fn eval_cos(&self, cos_theta: f64) -> Result<f64, VolumeTransportError> {
        self.validate()?;
        if !cos_theta.is_finite() || cos_theta.abs() > 1.0 + 2.0e-10 {
            return Err(VolumeTransportError::Phase);
        }
        let c = cos_theta.clamp(-1.0, 1.0);
        Ok(match self {
            Self::Isotropic => INV_FOUR_PI,
            Self::HenyeyGreenstein { g } => {
                let a = g.abs();
                // Keep the forward/backward peak stable as |g| approaches one.
                let denominator = if *g >= 0.0 {
                    (1.0 - a) * (1.0 - a) + 2.0 * a * (1.0 - c)
                } else {
                    (1.0 - a) * (1.0 - a) + 2.0 * a * (1.0 + c)
                };
                (1.0 - a) * (1.0 + a) * INV_FOUR_PI / (denominator * det::sqrt(denominator))
            }
            Self::Rayleigh => 0.75 * INV_FOUR_PI * (1.0 + c * c),
            Self::RayleighDepolarized { depolarization: d } => {
                1.5 * INV_FOUR_PI * ((1.0 + d) + (1.0 - d) * c * c) / (2.0 + d)
            }
            Self::Tabulated(table) => table.eval_clamped(c),
            Self::SpectralTable(_) => return Err(VolumeTransportError::Wavelength),
        })
    }

    /// Sample the exact angular CDF and a uniform azimuth from `[0,1)` variates.
    /// No small-anisotropy isotropic approximation is used.
    pub fn sample(
        &self,
        incoming_travel: Vec3,
        u_cos: f64,
        u_azimuth: f64,
    ) -> Result<PhaseSample, VolumeTransportError> {
        self.validate()?;
        let norm2 = incoming_travel.dot(incoming_travel);
        if !norm2.is_finite()
            || (norm2 - 1.0).abs() > 2.0e-10
            || !unit_variate(u_cos)
            || !unit_variate(u_azimuth)
        {
            return Err(VolumeTransportError::Phase);
        }
        let w = incoming_travel.scale(1.0 / det::sqrt(norm2));
        let cosine = match self {
            Self::Isotropic => 2.0f64.mul_add(u_cos, -1.0),
            Self::HenyeyGreenstein { g } => {
                // Algebraic HG inversion without subtracting nearly equal
                // terms then dividing by g; exact also at g=0 and small |g|.
                let a = g.abs();
                let u = if *g < 0.0 { 1.0 - u_cos } else { u_cos };
                let q = (1.0 - a) / ((1.0 - a) + 2.0 * a * u);
                let c = 1.0 - 2.0 * q * q * (1.0 - u) * (1.0 + a * u);
                if *g < 0.0 { -c } else { c }
            }
            Self::Rayleigh => rayleigh_cosine(u_cos),
            Self::RayleighDepolarized { depolarization: d } => {
                // Exact normalized Rayleigh/isotropic mixture. Conditional
                // remapping reuses the cosine variate without an extra stream.
                let isotropic = 3.0 * d / (2.0 + d);
                if u_cos < isotropic {
                    2.0 * (u_cos / isotropic) - 1.0
                } else {
                    rayleigh_cosine((u_cos - isotropic) / (1.0 - isotropic))
                }
            }
            Self::Tabulated(table) => table.sample_cosine(u_cos)?,
            Self::SpectralTable(_) => return Err(VolumeTransportError::Wavelength),
        }
        .clamp(-1.0, 1.0);
        let sine = det::sqrt(((1.0 - cosine) * (1.0 + cosine)).max(0.0));
        let (sin_phi, cos_phi) = det::sin_cos(core::f64::consts::TAU * u_azimuth);
        // Stable orthonormal frame, including both poles (Duff et al. 2017).
        let sign = 1.0f64.copysign(w.z);
        let a = -1.0 / (sign + w.z);
        let b = w.x * w.y * a;
        let tangent = Vec3::new(1.0 + sign * w.x * w.x * a, sign * b, -sign * w.x);
        let bitangent = Vec3::new(b, sign + w.y * w.y * a, -w.y);
        let x = sine * cos_phi;
        let y = sine * sin_phi;
        let raw = Vec3::new(
            tangent.x * x + bitangent.x * y + w.x * cosine,
            tangent.y * x + bitangent.y * y + w.y * cosine,
            tangent.z * x + bitangent.z * y + w.z * cosine,
        );
        let direction = raw.scale(1.0 / det::sqrt(raw.dot(raw)));
        let pdf = self.eval_cos(w.dot(direction))?;
        if let Self::Tabulated(table) = self
            && pdf == 0.0
        {
            // A zero-density CDF boundary has measure zero. A frame rotation
            // may also round an edge sample just outside very narrow support.
            // Use an actual interior direction in that same positive segment;
            // never replace the physical PDF by a nonzero floor.
            let cosine = table.interior_cosine(u_cos);
            let sine = det::sqrt(((1.0 - cosine) * (1.0 + cosine)).max(0.0));
            let x = sine * cos_phi;
            let y = sine * sin_phi;
            let raw = Vec3::new(
                tangent.x * x + bitangent.x * y + w.x * cosine,
                tangent.y * x + bitangent.y * y + w.y * cosine,
                tangent.z * x + bitangent.z * y + w.z * cosine,
            );
            let direction = raw.scale(1.0 / det::sqrt(raw.dot(raw)));
            let pdf = self.eval_cos(w.dot(direction))?;
            if !(pdf > 0.0 && pdf.is_finite()) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            return Ok(PhaseSample { direction, pdf });
        }
        Ok(PhaseSample { direction, pdf })
    }
}

fn rayleigh_cosine(u_cos: f64) -> f64 {
    // Invert c³+3c=8u-4 using 2*sinh(asinh(4u-2)/3).
    // Deterministic expm1 avoids Cardano cancellation near zero.
    let r = 4.0f64.mul_add(u_cos, -2.0);
    let t = det::ln(r.abs() + det::sqrt(r * r + 1.0)) / 3.0;
    let e = det::expm1(t);
    (e * (2.0 + e) / (1.0 + e)).copysign(r)
}
