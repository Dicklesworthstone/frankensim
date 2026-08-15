//! Anisotropic/layered surface detail and visualization-only spin cues
//! (bead frankensim-h7xu5.4.3, core slice).
//!
//! A perfectly axisymmetric, uniformly shaded disc reveals precession but
//! not spin about its own axis. Brushed-metal anisotropy and an engraved
//! radial mark make spin readable — as SHADING detail only:
//!
//! - The tangent frame, brushing rotation, anisotropic alpha pair, normal
//!   perturbation, and spin mark all live strictly on the shading side.
//!   Nothing here touches intersection, contact, mass, or motion vectors,
//!   and the geometric normal is carried unmodified next to the perturbed
//!   shading normal so a consumer cannot lose the distinction.
//! - Every parameter binds into a material-detail identity; changing the
//!   brushing pattern, anisotropy ratio, or mark parameters mints a new
//!   identity while the PHYSICAL chart/mass fingerprints (owned elsewhere)
//!   stay untouched by construction — this module cannot reach them.
//! - The zero-detail limit is exact: `SurfaceDetail::none()` reproduces
//!   the isotropic base surface bit-for-bit and perturbs nothing.
//!
//! No-claims: spin cues are presentation. The anisotropic distribution is
//! energy-consistent by construction (an alpha PAIR in the admitted GGX
//! domain reshapes the lobe; it neither gains energy nor rewrites optical
//! constants), but no radiometric validation is claimed here.

use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};

/// Admitted GGX alpha domain (mirrors the isotropic conductor surface).
pub const MIN_DETAIL_ALPHA: f64 = 1.0e-4;
/// Upper admitted GGX alpha.
pub const MAX_DETAIL_ALPHA: f64 = 1.0;
/// Hard cap on shading-normal perturbation (radians). Bounded so a shading
/// normal can never flip below the geometric horizon.
pub const MAX_NORMAL_PERTURBATION_RAD: f64 = 0.35;
/// Versioned identity domain for surface-detail parameter sets.
pub const SURFACE_DETAIL_IDENTITY_DOMAIN: &str = "org.frankensim.fs-render.surface-detail.v1";

/// Typed refusals of the surface-detail boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceDetailError {
    /// An alpha component is outside the admitted GGX domain.
    InvalidAlpha,
    /// The anisotropy ratio must be finite and >= 1 (the ratio stretches
    /// the tangential alpha; direction handles orientation).
    InvalidAnisotropyRatio,
    /// The normal perturbation exceeds the bounded cap or is non-finite.
    UnboundedNormalPerturbation,
    /// The mark angular width must lie in (0, pi/4].
    InvalidMarkWidth,
    /// The query point sits on the spin axis where no radial tangent
    /// exists; callers fall back to the geometric frame.
    PoleTangentUndefined,
    /// A non-finite local position or normal was supplied.
    NonFiniteInput,
}

impl core::fmt::Display for SurfaceDetailError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "surface-detail refusal: {self:?}")
    }
}

impl std::error::Error for SurfaceDetailError {}

/// Procedural brushing pattern over the axisymmetric chart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BrushingPattern {
    /// Grooves along the radial direction (spokes).
    Radial,
    /// Grooves along the circular direction (turntable finish).
    Circular,
}

impl BrushingPattern {
    const fn tag(self) -> u8 {
        match self {
            Self::Radial => 0,
            Self::Circular => 1,
        }
    }
}

/// Optional engraved/printed radial reference mark: a shading-layer-only
/// spin cue at one azimuth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialMark {
    /// Mark azimuth in the local chart frame [rad].
    pub azimuth_rad: f64,
    /// Angular half-width of the mark [rad], in (0, pi/4].
    pub half_width_rad: f64,
}

/// The complete, validated shading-detail parameter set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDetail {
    base_alpha: f64,
    anisotropy_ratio: f64,
    pattern: BrushingPattern,
    normal_perturbation_rad: f64,
    mark: Option<RadialMark>,
}

impl SurfaceDetail {
    /// The exact zero-detail limit: isotropic base alpha, no perturbation,
    /// no mark. Guaranteed to reproduce the base surface bit-for-bit.
    #[must_use]
    pub const fn none(base_alpha: f64) -> Self {
        Self {
            base_alpha,
            anisotropy_ratio: 1.0,
            pattern: BrushingPattern::Circular,
            normal_perturbation_rad: 0.0,
            mark: None,
        }
    }

    /// Validate a detail parameter set.
    ///
    /// # Errors
    /// Typed refusals for out-of-domain alpha, ratio below one,
    /// perturbations above the bounded cap, and malformed marks.
    pub fn try_new(
        base_alpha: f64,
        anisotropy_ratio: f64,
        pattern: BrushingPattern,
        normal_perturbation_rad: f64,
        mark: Option<RadialMark>,
    ) -> Result<Self, SurfaceDetailError> {
        if !base_alpha.is_finite() || !(MIN_DETAIL_ALPHA..=MAX_DETAIL_ALPHA).contains(&base_alpha) {
            return Err(SurfaceDetailError::InvalidAlpha);
        }
        if !anisotropy_ratio.is_finite() || anisotropy_ratio < 1.0 {
            return Err(SurfaceDetailError::InvalidAnisotropyRatio);
        }
        // The stretched alpha must stay inside the admitted domain too:
        // energy consistency comes from the domain, not from hope.
        if base_alpha * anisotropy_ratio > MAX_DETAIL_ALPHA {
            return Err(SurfaceDetailError::InvalidAlpha);
        }
        if !normal_perturbation_rad.is_finite()
            || normal_perturbation_rad < 0.0
            || normal_perturbation_rad > MAX_NORMAL_PERTURBATION_RAD
        {
            return Err(SurfaceDetailError::UnboundedNormalPerturbation);
        }
        if let Some(mark) = &mark {
            if !mark.azimuth_rad.is_finite()
                || !mark.half_width_rad.is_finite()
                || mark.half_width_rad <= 0.0
                || mark.half_width_rad > core::f64::consts::FRAC_PI_4
            {
                return Err(SurfaceDetailError::InvalidMarkWidth);
            }
        }
        Ok(Self {
            base_alpha,
            anisotropy_ratio,
            pattern,
            normal_perturbation_rad,
            mark,
        })
    }

    /// Anisotropic GGX alpha pair `(alpha_along_grooves, alpha_across)`.
    /// The lobe is TIGHT along the brushing grooves and wide across them;
    /// ratio 1 collapses to the isotropic base exactly.
    #[must_use]
    pub fn alpha_pair(&self) -> (f64, f64) {
        (self.base_alpha, self.base_alpha * self.anisotropy_ratio)
    }

    /// Whether this set is exactly the zero-detail limit.
    #[must_use]
    pub fn is_zero_detail(&self) -> bool {
        self.anisotropy_ratio == 1.0 && self.normal_perturbation_rad == 0.0 && self.mark.is_none()
    }

    /// Versioned identity of the complete parameter set. Every semantic
    /// field participates, so an anisotropy swap or mark edit mints a new
    /// material-detail identity.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut payload = Vec::new();
        payload.extend_from_slice(&self.base_alpha.to_bits().to_le_bytes());
        payload.extend_from_slice(&self.anisotropy_ratio.to_bits().to_le_bytes());
        payload.push(self.pattern.tag());
        payload.extend_from_slice(&self.normal_perturbation_rad.to_bits().to_le_bytes());
        match &self.mark {
            None => payload.push(0),
            Some(mark) => {
                payload.push(1);
                payload.extend_from_slice(&mark.azimuth_rad.to_bits().to_le_bytes());
                payload.extend_from_slice(&mark.half_width_rad.to_bits().to_le_bytes());
            }
        }
        fs_blake3::hash_domain(SURFACE_DETAIL_IDENTITY_DOMAIN, &payload)
    }

    /// Whether the local chart point falls inside the radial mark.
    #[must_use]
    pub fn mark_covers(&self, local: Point3) -> bool {
        let Some(mark) = &self.mark else {
            return false;
        };
        let azimuth = local.y.atan2(local.x);
        let mut delta = (azimuth - mark.azimuth_rad) % core::f64::consts::TAU;
        if delta > core::f64::consts::PI {
            delta -= core::f64::consts::TAU;
        }
        if delta < -core::f64::consts::PI {
            delta += core::f64::consts::TAU;
        }
        delta.abs() <= mark.half_width_rad
    }
}

/// The shading frame at one local chart point: geometric normal carried
/// UNMODIFIED next to the derived tangent basis and the (bounded) shading
/// normal, so downstream code cannot conflate them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShadingFrame {
    /// Geometric normal, exactly as supplied (unit).
    pub geometric_normal: Vec3,
    /// Bounded shading normal (unit; equals geometric under zero detail).
    pub shading_normal: Vec3,
    /// First tangent: along the brushing grooves (unit).
    pub tangent: Vec3,
    /// Second tangent: `geometric_normal x tangent` (unit; right-handed).
    pub bitangent: Vec3,
}

/// Derive the deterministic tangent frame for the axisymmetric chart at a
/// local point with the given geometric normal.
///
/// The radial direction projects the local position onto the plane
/// perpendicular to the normal; [`BrushingPattern::Radial`] grooves run
/// along it, [`BrushingPattern::Circular`] grooves run along the circular
/// direction (`normal x radial`). The frame is right-handed by
/// construction and rotates WITH the body: it is a pure function of local
/// coordinates, so a rigid rotation of (position, normal) rotates the
/// frame identically.
///
/// # Errors
/// Pole refusal where no radial tangent exists (callers fall back to the
/// geometric frame), non-finite inputs.
pub fn shading_frame(
    detail: &SurfaceDetail,
    local: Point3,
    geometric_normal: Vec3,
) -> Result<ShadingFrame, SurfaceDetailError> {
    if !local.x.is_finite()
        || !local.y.is_finite()
        || !local.z.is_finite()
        || !geometric_normal.x.is_finite()
        || !geometric_normal.y.is_finite()
        || !geometric_normal.z.is_finite()
    {
        return Err(SurfaceDetailError::NonFiniteInput);
    }
    let normal = normalize(geometric_normal).ok_or(SurfaceDetailError::NonFiniteInput)?;
    // Radial component of the local position perpendicular to the normal.
    let position = Vec3::new(local.x, local.y, local.z);
    let along_normal = dot(position, normal);
    let radial = Vec3::new(
        position.x - along_normal * normal.x,
        position.y - along_normal * normal.y,
        position.z - along_normal * normal.z,
    );
    let Some(radial) = normalize(radial) else {
        return Err(SurfaceDetailError::PoleTangentUndefined);
    };
    let circular = cross(normal, radial);
    let (tangent, bitangent) = match detail.pattern {
        BrushingPattern::Radial => (radial, circular),
        BrushingPattern::Circular => (circular, negate(radial)),
    };
    // Bounded shading-normal perturbation: tilt toward the bitangent by
    // the fixed admitted angle. Deterministic (no noise source), bounded
    // by admission, and exactly the geometric normal at zero detail.
    let shading_normal = if detail.normal_perturbation_rad == 0.0 {
        normal
    } else {
        let (sine, cosine) = detail.normal_perturbation_rad.sin_cos();
        normalize(Vec3::new(
            cosine * normal.x + sine * bitangent.x,
            cosine * normal.y + sine * bitangent.y,
            cosine * normal.z + sine * bitangent.z,
        ))
        .ok_or(SurfaceDetailError::NonFiniteInput)?
    };
    Ok(ShadingFrame {
        geometric_normal: normal,
        shading_normal,
        tangent,
        bitangent,
    })
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn negate(v: Vec3) -> Vec3 {
    Vec3::new(-v.x, -v.y, -v.z)
}

fn normalize(v: Vec3) -> Option<Vec3> {
    let length = dot(v, v).sqrt();
    if !length.is_finite() || length < 1.0e-12 {
        return None;
    }
    Some(Vec3::new(v.x / length, v.y / length, v.z / length))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn detail() -> SurfaceDetail {
        SurfaceDetail::try_new(
            0.08,
            4.0,
            BrushingPattern::Circular,
            0.05,
            Some(RadialMark {
                azimuth_rad: 0.0,
                half_width_rad: 0.1,
            }),
        )
        .expect("admits")
    }

    #[test]
    fn admission_fails_closed() {
        assert_eq!(
            SurfaceDetail::try_new(0.0, 1.0, BrushingPattern::Radial, 0.0, None).unwrap_err(),
            SurfaceDetailError::InvalidAlpha
        );
        assert_eq!(
            SurfaceDetail::try_new(0.5, 0.5, BrushingPattern::Radial, 0.0, None).unwrap_err(),
            SurfaceDetailError::InvalidAnisotropyRatio
        );
        // Stretched alpha escaping the domain refuses (energy discipline).
        assert_eq!(
            SurfaceDetail::try_new(0.5, 4.0, BrushingPattern::Radial, 0.0, None).unwrap_err(),
            SurfaceDetailError::InvalidAlpha
        );
        assert_eq!(
            SurfaceDetail::try_new(0.1, 2.0, BrushingPattern::Radial, 1.0, None).unwrap_err(),
            SurfaceDetailError::UnboundedNormalPerturbation
        );
        assert_eq!(
            SurfaceDetail::try_new(
                0.1,
                2.0,
                BrushingPattern::Radial,
                0.0,
                Some(RadialMark {
                    azimuth_rad: 0.0,
                    half_width_rad: 2.0,
                })
            )
            .unwrap_err(),
            SurfaceDetailError::InvalidMarkWidth
        );
    }

    #[test]
    fn zero_detail_limit_is_exact() {
        let base = SurfaceDetail::none(0.08);
        assert!(base.is_zero_detail());
        assert_eq!(base.alpha_pair(), (0.08, 0.08), "isotropic recovery");
        let frame = shading_frame(&base, Point3::new(0.5, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
            .expect("frame");
        assert_eq!(
            frame.shading_normal, frame.geometric_normal,
            "zero detail perturbs nothing"
        );
        assert!(!base.mark_covers(Point3::new(0.5, 0.0, 0.0)));
    }

    #[test]
    fn tangent_frame_is_right_handed_and_rotates_with_the_body() {
        let detail = detail();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let frame = shading_frame(&detail, Point3::new(0.5, 0.0, 0.0), normal).expect("frame");
        // Right-handedness: tangent x bitangent points along the normal.
        let handed = cross(frame.tangent, frame.bitangent);
        assert!(
            dot(handed, frame.geometric_normal) > 0.99,
            "right-handed frame"
        );
        // Circular brushing at +x: grooves run along +y.
        assert!((frame.tangent.y - 1.0).abs() < 1e-12);

        // Rigid rotation covariance: rotate the local point 90 degrees
        // about z; the tangent rotates identically (spin readability -
        // the frame follows the body).
        let rotated = shading_frame(&detail, Point3::new(0.0, 0.5, 0.0), normal).expect("frame");
        assert!(
            (rotated.tangent.x + 1.0).abs() < 1e-12,
            "tangent follows rotation"
        );
        // Bounded perturbation: shading normal within the declared cone.
        let cosine = dot(frame.shading_normal, frame.geometric_normal);
        assert!(cosine >= (MAX_NORMAL_PERTURBATION_RAD).cos() - 1e-12);
    }

    #[test]
    fn pole_and_nonfinite_inputs_refuse() {
        let detail = detail();
        assert_eq!(
            shading_frame(
                &detail,
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0)
            )
            .unwrap_err(),
            SurfaceDetailError::PoleTangentUndefined
        );
        assert_eq!(
            shading_frame(
                &detail,
                Point3::new(f64::NAN, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0)
            )
            .unwrap_err(),
            SurfaceDetailError::NonFiniteInput
        );
    }

    #[test]
    fn identity_binds_every_parameter() {
        let base = detail();
        let base_id = base.identity();
        let anisotropy_swap =
            SurfaceDetail::try_new(0.08, 5.0, BrushingPattern::Circular, 0.05, base_mark())
                .expect("admits");
        assert_ne!(
            anisotropy_swap.identity(),
            base_id,
            "anisotropy swap moves identity"
        );
        let pattern_swap =
            SurfaceDetail::try_new(0.08, 4.0, BrushingPattern::Radial, 0.05, base_mark())
                .expect("admits");
        assert_ne!(
            pattern_swap.identity(),
            base_id,
            "pattern swap moves identity"
        );
        let mark_removed = SurfaceDetail::try_new(0.08, 4.0, BrushingPattern::Circular, 0.05, None)
            .expect("admits");
        assert_ne!(mark_removed.identity(), base_id, "mark edit moves identity");
        // And re-deriving the same parameters reproduces the identity.
        assert_eq!(detail().identity(), base_id);
    }

    fn base_mark() -> Option<RadialMark> {
        Some(RadialMark {
            azimuth_rad: 0.0,
            half_width_rad: 0.1,
        })
    }

    #[test]
    fn the_radial_mark_reads_spin_deterministically() {
        let detail = detail();
        // The mark sits at azimuth 0 with half-width 0.1 rad.
        assert!(detail.mark_covers(Point3::new(1.0, 0.0, 0.0)));
        assert!(detail.mark_covers(Point3::new(1.0, 0.09, 0.0)));
        assert!(!detail.mark_covers(Point3::new(1.0, 0.2, 0.0)));
        assert!(!detail.mark_covers(Point3::new(-1.0, 0.0, 0.0)));
        // Wraparound: an azimuth just below +pi is far from the mark; the
        // delta normalization must not misclassify it.
        assert!(!detail.mark_covers(Point3::new(-1.0, -1.0e-6, 0.0)));
    }
}
