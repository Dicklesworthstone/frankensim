//! Solid mass properties for validated axisymmetric line/arc charts.
//!
//! The formulas here are boundary integrals, not a tessellation or sampled
//! volumetric quadrature.  For a counter-clockwise meridian boundary `C` in
//! `(rho, z)`, Green's theorem gives
//!
//! ```text
//! V       = pi * integral_C rho^2 dz
//! V zbar  = pi * integral_C rho^2 z dz
//! I_z     = rho_m pi/2 * integral_C rho^4 dz
//! I_perp  = rho_m pi * integral_C (rho^4/4 + rho^2 z^2) dz.
//! ```
//!
//! Lines are integrated as ordinary polynomials.  Arcs expand into bounded
//! powers of sine and cosine and use their closed antiderivatives.  Thus an
//! admitted profile stays an exact geometric profile throughout mechanics
//! admission; no chordal approximation is introduced here.

use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_geom::Point3;

use crate::{AxisymmetricChart, AxisymmetricError, MeridianSegment};

const PI: f64 = core::f64::consts::PI;
const TAU: f64 = core::f64::consts::TAU;
const POLL_STRIDE: usize = 16;
const CERTIFICATE_ULPS: usize = 4_096;

/// Axisymmetric principal moments about either the center of mass or the
/// world origin.  `transverse` is both `I_x` and `I_y`; `axial` is `I_z`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricPrincipalInertia {
    /// The repeated principal moment about an axis perpendicular to z.
    pub transverse: f64,
    /// The principal moment about the revolution (z) axis.
    pub axial: f64,
}

/// Explicit outward bands for the deterministic floating-point evaluation.
///
/// The analytic integration is exact over the represented real line/arc
/// geometry.  These bands conservatively widen the fixed operation sequence,
/// including the declared strict-trigonometry ULP budget used for arcs; they
/// do not turn binary64 input coordinates into exact real measurements.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricMassErrorBounds {
    /// Enclosure for solid volume.
    pub volume: NumericalCertificate,
    /// Enclosure for mass after density multiplication.
    pub mass: NumericalCertificate,
    /// Enclosure for the axial center-of-mass coordinate.
    pub center_axial: NumericalCertificate,
    /// Enclosure for transverse inertia about the world origin.
    pub origin_transverse: NumericalCertificate,
    /// Enclosure for axial inertia about the world origin.
    pub origin_axial: NumericalCertificate,
    /// Enclosure for transverse centroidal principal inertia.
    pub principal_transverse: NumericalCertificate,
    /// Enclosure for axial centroidal principal inertia.
    pub principal_axial: NumericalCertificate,
}

/// Production mass properties of a homogeneous solid of revolution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricMassProperties {
    /// Positive enclosed volume.
    pub volume: f64,
    /// `density * volume`.
    pub mass: f64,
    /// Center of mass.  Axisymmetry makes `x == y == 0` exactly by model.
    pub center_of_mass: Point3,
    /// Principal inertia about the center of mass.
    pub principal_inertia: AxisymmetricPrincipalInertia,
    /// Inertia about the world origin, retained for parallel-axis consumers.
    pub origin_inertia: AxisymmetricPrincipalInertia,
    /// Explicit floating-point output bands for each reported quantity.
    pub error: AxisymmetricMassErrorBounds,
}

/// Refusal from [`AxisymmetricChart::mass_properties`].
#[derive(Debug, Clone, PartialEq)]
pub enum AxisymmetricMassError {
    /// Density is not a finite, strictly positive volumetric density.
    InvalidDensity {
        /// Rejected value.
        density: f64,
    },
    /// The retained meridian no longer meets its construction obligations.
    InvalidChart(AxisymmetricError),
    /// Cancellation was observed before all features were integrated.
    Cancelled,
    /// A fixed exact-formula evaluation overflowed or otherwise became
    /// non-finite; no partial mass properties are published.
    NonFinite {
        /// Quantity that failed to remain finite.
        quantity: &'static str,
    },
    /// A validated CCW chart reconstructed a non-positive volume, so it is
    /// unsafe to divide by it for a center of mass.
    NonPositiveVolume {
        /// Reconstructed volume.
        volume: f64,
    },
    /// A central inertia violated non-negativity beyond rounding tolerance.
    InvalidInertia {
        /// Reconstructed transverse central inertia.
        transverse: f64,
        /// Reconstructed axial central inertia.
        axial: f64,
    },
}

impl core::fmt::Display for AxisymmetricMassError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDensity { density } => write!(
                f,
                "axisymmetric mass density must be finite and strictly positive, got {density}"
            ),
            Self::InvalidChart(error) => write!(
                f,
                "axisymmetric mass requires a still-valid closed CCW meridian: {error}"
            ),
            Self::Cancelled => write!(
                f,
                "axisymmetric mass integration was cancelled before publishing properties"
            ),
            Self::NonFinite { quantity } => write!(
                f,
                "axisymmetric mass exact-formula evaluation became non-finite at {quantity}"
            ),
            Self::NonPositiveVolume { volume } => write!(
                f,
                "axisymmetric mass reconstructed non-positive volume {volume} from its CCW meridian"
            ),
            Self::InvalidInertia { transverse, axial } => write!(
                f,
                "axisymmetric mass reconstructed invalid central inertia (transverse={transverse}, axial={axial})"
            ),
        }
    }
}

impl std::error::Error for AxisymmetricMassError {}

impl AxisymmetricChart {
    /// Compute homogeneous-solid mass properties from this chart's exact
    /// meridian features.
    ///
    /// The retained construction certificate is independently revalidated
    /// before integration.  Work is deterministic in segment order and polls
    /// cancellation at bounded feature strides.  The returned principal
    /// inertia is centroidal; [`AxisymmetricMassProperties::origin_inertia`]
    /// is supplied so callers can independently check/apply the parallel-axis
    /// theorem without hand-entering a disc inertia.
    pub fn mass_properties(
        &self,
        density: f64,
        cx: &Cx<'_>,
    ) -> Result<AxisymmetricMassProperties, AxisymmetricMassError> {
        if !density.is_finite() || density <= 0.0 {
            return Err(AxisymmetricMassError::InvalidDensity { density });
        }
        self.verify_construction()
            .map_err(AxisymmetricMassError::InvalidChart)?;

        let mut moments = BoundaryMoments::default();
        for (index, segment) in self.segments().iter().copied().enumerate() {
            if index % POLL_STRIDE == 0 && cx.checkpoint().is_err() {
                return Err(AxisymmetricMassError::Cancelled);
            }
            moments.add(segment)?;
        }
        if cx.checkpoint().is_err() {
            return Err(AxisymmetricMassError::Cancelled);
        }

        let volume = PI * moments.r2;
        let first_axial = PI * moments.r2_z;
        let origin_axial = density * (0.5 * PI * moments.r4);
        let origin_transverse = density * PI * (0.25 * moments.r4 + moments.r2_z2);
        ensure_finite(volume, "volume")?;
        ensure_finite(first_axial, "axial first moment")?;
        ensure_finite(origin_axial, "origin axial inertia")?;
        ensure_finite(origin_transverse, "origin transverse inertia")?;
        if volume <= 0.0 {
            return Err(AxisymmetricMassError::NonPositiveVolume { volume });
        }

        let mass = density * volume;
        let center_axial = first_axial / volume;
        let mut principal_transverse = origin_transverse - mass * center_axial * center_axial;
        let principal_axial = origin_axial;
        ensure_finite(mass, "mass")?;
        ensure_finite(center_axial, "center axial coordinate")?;
        ensure_finite(principal_transverse, "principal transverse inertia")?;
        ensure_finite(principal_axial, "principal axial inertia")?;

        // A centered subtraction can leave a few negative ulps for an
        // otherwise valid solid.  Only canonicalize that representational
        // residue; a material negative moment remains a hard refusal.
        let inertia_scale = origin_transverse
            .abs()
            .max(mass * center_axial * center_axial);
        let roundoff = inertia_scale * 64.0 * f64::EPSILON;
        if principal_transverse < 0.0 && principal_transverse >= -roundoff {
            principal_transverse = 0.0;
        }
        if principal_transverse < 0.0 || principal_axial < 0.0 {
            return Err(AxisymmetricMassError::InvalidInertia {
                transverse: principal_transverse,
                axial: principal_axial,
            });
        }

        Ok(AxisymmetricMassProperties {
            volume,
            mass,
            center_of_mass: Point3::new(0.0, 0.0, center_axial),
            principal_inertia: AxisymmetricPrincipalInertia {
                transverse: principal_transverse,
                axial: principal_axial,
            },
            origin_inertia: AxisymmetricPrincipalInertia {
                transverse: origin_transverse,
                axial: origin_axial,
            },
            error: AxisymmetricMassErrorBounds {
                volume: widened(volume),
                mass: widened(mass),
                center_axial: widened(center_axial),
                origin_transverse: widened(origin_transverse),
                origin_axial: widened(origin_axial),
                principal_transverse: widened(principal_transverse),
                principal_axial: widened(principal_axial),
            },
        })
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct BoundaryMoments {
    r2: f64,
    r2_z: f64,
    r4: f64,
    r2_z2: f64,
}

impl BoundaryMoments {
    fn add(&mut self, segment: MeridianSegment) -> Result<(), AxisymmetricMassError> {
        self.r2 += boundary_integral(segment, 2, 0);
        self.r2_z += boundary_integral(segment, 2, 1);
        self.r4 += boundary_integral(segment, 4, 0);
        self.r2_z2 += boundary_integral(segment, 2, 2);
        ensure_finite(self.r2, "volume boundary moment")?;
        ensure_finite(self.r2_z, "axial first boundary moment")?;
        ensure_finite(self.r4, "axial inertia boundary moment")?;
        ensure_finite(self.r2_z2, "transverse inertia boundary moment")?;
        Ok(())
    }
}

/// Evaluate `integral r^radius_power z^axial_power dz` over one feature.
fn boundary_integral(segment: MeridianSegment, radius_power: u32, axial_power: u32) -> f64 {
    match segment {
        MeridianSegment::Line { start, end } => line_integral(
            start.radius,
            start.axial,
            end.radius - start.radius,
            end.axial - start.axial,
            radius_power,
            axial_power,
        ),
        MeridianSegment::Arc {
            start,
            end,
            center,
            clockwise,
        } => arc_integral(
            start.radius,
            start.axial,
            end.radius,
            end.axial,
            center.radius,
            center.axial,
            clockwise,
            radius_power,
            axial_power,
        ),
    }
}

fn line_integral(
    radius: f64,
    axial: f64,
    delta_radius: f64,
    delta_axial: f64,
    radius_power: u32,
    axial_power: u32,
) -> f64 {
    let mut sum = 0.0;
    for radius_term in 0..=radius_power {
        let radius_coefficient = binomial(radius_power, radius_term)
            * integer_power(radius, radius_power - radius_term)
            * integer_power(delta_radius, radius_term);
        for axial_term in 0..=axial_power {
            let axial_coefficient = binomial(axial_power, axial_term)
                * integer_power(axial, axial_power - axial_term)
                * integer_power(delta_axial, axial_term);
            sum += delta_axial * radius_coefficient * axial_coefficient
                / f64::from(radius_term + axial_term + 1);
        }
    }
    sum
}

#[allow(clippy::too_many_arguments)]
fn arc_integral(
    start_radius: f64,
    start_axial: f64,
    end_radius: f64,
    end_axial: f64,
    center_radius: f64,
    center_axial: f64,
    clockwise: bool,
    radius_power: u32,
    axial_power: u32,
) -> f64 {
    let start_angle = (start_axial - center_axial).atan2(start_radius - center_radius);
    let end_angle = (end_axial - center_axial).atan2(end_radius - center_radius);
    let sweep = if clockwise {
        -(start_angle - end_angle).rem_euclid(TAU)
    } else {
        (end_angle - start_angle).rem_euclid(TAU)
    };
    let finish_angle = start_angle + sweep;
    let arc_radius = (start_radius - center_radius).hypot(start_axial - center_axial);

    let mut sum = 0.0;
    for radius_term in 0..=radius_power {
        let radius_coefficient = binomial(radius_power, radius_term)
            * integer_power(center_radius, radius_power - radius_term)
            * integer_power(arc_radius, radius_term);
        for axial_term in 0..=axial_power {
            let axial_coefficient = binomial(axial_power, axial_term)
                * integer_power(center_axial, axial_power - axial_term)
                * integer_power(arc_radius, axial_term);
            sum += radius_coefficient
                * axial_coefficient
                * arc_radius
                * trig_power_integral(radius_term + 1, axial_term, start_angle, finish_angle);
        }
    }
    sum
}

/// Exact antiderivative evaluation for `cos(theta)^m sin(theta)^n`.
///
/// The only exponents needed by the mechanics formulas are `m <= 5`,
/// `n <= 2`; expansion into complex exponentials keeps that bounded and
/// avoids a quadrature tolerance or an input-dependent refinement path.
fn trig_power_integral(cos_power: u32, sin_power: u32, start: f64, end: f64) -> f64 {
    let mut sum = 0.0;
    let denominator = 2.0_f64.powi((cos_power + sin_power) as i32);
    for cos_pick in 0..=cos_power {
        for sin_pick in 0..=sin_power {
            let sign = if sin_pick % 2 == 0 { 1.0 } else { -1.0 };
            let scalar =
                sign * binomial(cos_power, cos_pick) * binomial(sin_power, sin_pick) / denominator;
            let frequency = (cos_power + sin_power) as i32 - 2 * (cos_pick + sin_pick) as i32;
            let (integral_real, integral_imaginary) = if frequency == 0 {
                (end - start, 0.0)
            } else {
                let frequency = f64::from(frequency);
                (
                    (fs_math::det::sin(frequency * end) - fs_math::det::sin(frequency * start))
                        / frequency,
                    (fs_math::det::cos(frequency * start) - fs_math::det::cos(frequency * end))
                        / frequency,
                )
            };
            let (coefficient_real, coefficient_imaginary) = match sin_power % 4 {
                0 => (scalar, 0.0),
                1 => (0.0, -scalar),
                2 => (-scalar, 0.0),
                _ => (0.0, scalar),
            };
            sum += coefficient_real * integral_real - coefficient_imaginary * integral_imaginary;
        }
    }
    sum
}

fn integer_power(base: f64, exponent: u32) -> f64 {
    (0..exponent).fold(1.0, |product, _| product * base)
}

fn binomial(n: u32, k: u32) -> f64 {
    let mut numerator = 1_u64;
    let mut denominator = 1_u64;
    for factor in 1..=k.min(n - k) {
        numerator *= u64::from(n + 1 - factor);
        denominator *= u64::from(factor);
    }
    (numerator / denominator) as f64
}

fn ensure_finite(value: f64, quantity: &'static str) -> Result<(), AxisymmetricMassError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(AxisymmetricMassError::NonFinite { quantity })
    }
}

fn widened(value: f64) -> NumericalCertificate {
    let mut lo = value;
    let mut hi = value;
    for _ in 0..CERTIFICATE_ULPS {
        lo = lo.next_down();
        hi = hi.next_up();
    }
    NumericalCertificate::enclosure(lo, hi)
}
