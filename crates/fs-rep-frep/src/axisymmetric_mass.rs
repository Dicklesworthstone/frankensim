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
//! The exposed surface area is evaluated from the same retained meridian:
//! `A = 2 pi integral_C rho ds`. Lines use the frustum formula and circular
//! arcs use its closed trigonometric antiderivative, so thermal and acoustic
//! consumers do not need a separately typed or tessellated area.
//!
//! Lines are integrated as ordinary polynomials.  Arcs expand into bounded
//! powers of sine and cosine and use their closed antiderivatives.  Thus an
//! admitted profile stays an exact geometric profile throughout mechanics
//! admission; no chordal approximation is introduced here.

use fs_exec::Cx;
use fs_geom::Point3;

use crate::{AxisymmetricChart, AxisymmetricError, MeridianSegment};

const PI: f64 = core::f64::consts::PI;
const TAU: f64 = core::f64::consts::TAU;
const POLL_STRIDE: usize = 16;

/// Axisymmetric principal moments about either the center of mass or the
/// world origin.  `transverse` is both `I_x` and `I_y`; `axial` is `I_z`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricPrincipalInertia {
    /// The repeated principal moment about an axis perpendicular to z.
    pub transverse: f64,
    /// The principal moment about the revolution (z) axis.
    pub axial: f64,
}

/// Non-authoritative numerical telemetry for a mass-property evaluation.
///
/// These are absolute term-magnitude scales accumulated in deterministic
/// feature order. They help a caller identify cancellation-prone inputs, but
/// are **not error bounds, intervals, or certificates** and must not be used
/// to admit a mechanics result. The analytic formulas are exact over the
/// represented real line/arc geometry; floating-point certification needs a
/// separate directed-rounding implementation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricMassRoundoffDiagnostics {
    /// Sum of absolute volume-boundary terms, after multiplication by pi.
    pub volume_term_scale: f64,
    /// Sum of absolute axial-first-moment boundary terms, after multiplication
    /// by pi.
    pub axial_first_moment_term_scale: f64,
    /// Sum of absolute centroidal transverse-inertia boundary terms.
    pub centroidal_transverse_term_scale: f64,
    /// Sum of absolute axial-inertia boundary terms.
    pub axial_inertia_term_scale: f64,
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
    /// Non-authoritative term-magnitude telemetry; never an error bound.
    pub roundoff_diagnostics: AxisymmetricMassRoundoffDiagnostics,
}

/// Analytic exposed area of one validated solid of revolution.
///
/// This is the complete boundary area, including inner bores and planar caps.
/// Axis-closing meridian segments contribute exactly zero because their
/// revolution is not a surface. The value is deterministic binary64
/// evaluation of exact line/arc formulas, not a directed-rounding certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricSurfaceArea {
    /// Complete exposed boundary area [length^2].
    pub area: f64,
}

/// Refusal from [`AxisymmetricChart::surface_area`].
#[derive(Debug, Clone, PartialEq)]
pub enum AxisymmetricSurfaceAreaError {
    /// The retained meridian no longer meets its construction obligations.
    InvalidChart(AxisymmetricError),
    /// Cancellation was observed before all features were integrated.
    Cancelled,
    /// An exact-formula evaluation overflowed or otherwise became non-finite.
    NonFinite {
        /// Retained meridian feature that failed.
        source_feature: usize,
    },
    /// A validated closed solid reconstructed a non-positive exposed area.
    NonPositiveArea {
        /// Reconstructed total area.
        area: f64,
    },
}

impl core::fmt::Display for AxisymmetricSurfaceAreaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidChart(error) => write!(
                f,
                "axisymmetric area requires a still-valid closed meridian: {error}"
            ),
            Self::Cancelled => write!(
                f,
                "axisymmetric area integration was cancelled before publishing a value"
            ),
            Self::NonFinite { source_feature } => write!(
                f,
                "axisymmetric area became non-finite at meridian feature {source_feature}"
            ),
            Self::NonPositiveArea { area } => write!(
                f,
                "axisymmetric area reconstructed non-positive boundary area {area}"
            ),
        }
    }
}

impl std::error::Error for AxisymmetricSurfaceAreaError {}

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
    /// A positive finite density and volume underflowed to a non-positive
    /// representable mass, so no mechanics properties can be published.
    NonPositiveMass {
        /// Reconstructed representable mass.
        mass: f64,
    },
    /// A positive-volume 3D solid reconstructed a zero principal moment after
    /// finite binary64 evaluation, so it cannot be published as a rigid body.
    NonPositivePrincipalInertia {
        /// Reconstructed transverse centroidal moment.
        transverse: f64,
        /// Reconstructed axial centroidal moment.
        axial: f64,
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
            Self::NonPositiveMass { mass } => write!(
                f,
                "axisymmetric mass underflowed to non-positive representable mass {mass}"
            ),
            Self::NonPositivePrincipalInertia { transverse, axial } => write!(
                f,
                "axisymmetric mass reconstructed non-positive principal inertia for a positive-volume solid (transverse={transverse}, axial={axial})"
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
    /// Compute complete exposed area from this chart's exact meridian.
    ///
    /// The construction certificate is revalidated, work follows retained
    /// feature order, and cancellation is polled at the same bounded stride as
    /// mass integration. No render mesh or caller-provided area participates.
    pub fn surface_area(
        &self,
        cx: &Cx<'_>,
    ) -> Result<AxisymmetricSurfaceArea, AxisymmetricSurfaceAreaError> {
        match self.verify_construction_with_cx(cx) {
            Ok(_) => {}
            Err(AxisymmetricError::Cancelled) => {
                return Err(AxisymmetricSurfaceAreaError::Cancelled);
            }
            Err(error) => return Err(AxisymmetricSurfaceAreaError::InvalidChart(error)),
        }

        let mut area = 0.0;
        for (index, segment) in self.segments().iter().copied().enumerate() {
            if index % POLL_STRIDE == 0 && cx.checkpoint().is_err() {
                return Err(AxisymmetricSurfaceAreaError::Cancelled);
            }
            let contribution = revolved_segment_area(segment);
            if !contribution.is_finite() {
                return Err(AxisymmetricSurfaceAreaError::NonFinite {
                    source_feature: index,
                });
            }
            area += contribution;
            if !area.is_finite() {
                return Err(AxisymmetricSurfaceAreaError::NonFinite {
                    source_feature: index,
                });
            }
        }
        if cx.checkpoint().is_err() {
            return Err(AxisymmetricSurfaceAreaError::Cancelled);
        }
        if area <= 0.0 {
            return Err(AxisymmetricSurfaceAreaError::NonPositiveArea { area });
        }
        Ok(AxisymmetricSurfaceArea { area })
    }

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
        match self.verify_construction_with_cx(cx) {
            Ok(_) => {}
            Err(AxisymmetricError::Cancelled) => return Err(AxisymmetricMassError::Cancelled),
            Err(error) => return Err(AxisymmetricMassError::InvalidChart(error)),
        }

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
        ensure_finite(volume, "volume")?;
        ensure_finite(first_axial, "axial first moment")?;
        ensure_finite(origin_axial, "origin axial inertia")?;
        if volume <= 0.0 {
            return Err(AxisymmetricMassError::NonPositiveVolume { volume });
        }

        let mass = density * volume;
        ensure_finite(mass, "mass")?;
        if mass <= 0.0 {
            return Err(AxisymmetricMassError::NonPositiveMass { mass });
        }
        let center_axial = first_axial / volume;
        ensure_finite(center_axial, "center axial coordinate")?;
        let mut centered_r2_z2 = 0.0;
        let mut centered_r2_z2_abs = 0.0;
        for (index, segment) in self.segments().iter().copied().enumerate() {
            if index % POLL_STRIDE == 0 && cx.checkpoint().is_err() {
                return Err(AxisymmetricMassError::Cancelled);
            }
            let term = boundary_integral_centered_axial(segment, 2, 2, center_axial);
            centered_r2_z2 += term;
            centered_r2_z2_abs += term.abs();
            ensure_finite(
                centered_r2_z2,
                "centered transverse inertia boundary moment",
            )?;
            ensure_finite(
                centered_r2_z2_abs,
                "centered transverse inertia absolute boundary moment",
            )?;
        }
        if cx.checkpoint().is_err() {
            return Err(AxisymmetricMassError::Cancelled);
        }
        let mut principal_transverse = density * PI * (0.25 * moments.r4 + centered_r2_z2);
        let principal_axial = origin_axial;
        ensure_finite(principal_transverse, "principal transverse inertia")?;
        ensure_finite(principal_axial, "principal axial inertia")?;

        // The centroidal moment is evaluated in a centered second boundary
        // pass. Only the final boundary sum can leave a few negative ulps;
        // a material negative moment remains a hard refusal.
        let inertia_scale = density * PI * (0.25 * moments.r4_abs + centered_r2_z2_abs);
        ensure_finite(inertia_scale, "principal transverse inertia scale")?;
        let roundoff = inertia_scale * 64.0 * f64::EPSILON;
        ensure_finite(roundoff, "principal transverse roundoff diagnostic")?;
        if principal_transverse < 0.0 && principal_transverse >= -roundoff {
            principal_transverse = 0.0;
        }
        if principal_transverse < 0.0 || principal_axial < 0.0 {
            return Err(AxisymmetricMassError::InvalidInertia {
                transverse: principal_transverse,
                axial: principal_axial,
            });
        }
        if principal_transverse == 0.0 || principal_axial == 0.0 {
            return Err(AxisymmetricMassError::NonPositivePrincipalInertia {
                transverse: principal_transverse,
                axial: principal_axial,
            });
        }

        // Compute this only after the final centroidal value has been
        // accepted, so the published pair obeys the parallel-axis identity.
        let origin_transverse = principal_transverse + mass * center_axial * center_axial;
        ensure_finite(origin_transverse, "origin transverse inertia")?;

        // These term scales are diagnostics only, but publication remains
        // transactional: every returned scalar must be finite after its own
        // density/pi multipliers, not merely before them.
        let volume_term_scale = PI * moments.r2_abs;
        let axial_first_moment_term_scale = PI * moments.r2_z_abs;
        let centroidal_transverse_term_scale = inertia_scale;
        let axial_inertia_term_scale = density * 0.5 * PI * moments.r4_abs;
        ensure_finite(volume_term_scale, "volume roundoff diagnostic")?;
        ensure_finite(
            axial_first_moment_term_scale,
            "axial first-moment roundoff diagnostic",
        )?;
        ensure_finite(
            centroidal_transverse_term_scale,
            "centroidal transverse roundoff diagnostic",
        )?;
        ensure_finite(
            axial_inertia_term_scale,
            "axial inertia roundoff diagnostic",
        )?;

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
            roundoff_diagnostics: AxisymmetricMassRoundoffDiagnostics {
                volume_term_scale,
                axial_first_moment_term_scale,
                centroidal_transverse_term_scale,
                axial_inertia_term_scale,
            },
        })
    }
}

/// Evaluate `2 pi integral radius ds` over one retained meridian feature.
fn revolved_segment_area(segment: MeridianSegment) -> f64 {
    match segment {
        MeridianSegment::Line { start, end } => {
            let slant = (end.radius - start.radius).hypot(end.axial - start.axial);
            PI * (start.radius + end.radius) * slant
        }
        MeridianSegment::Arc {
            start,
            end,
            center,
            clockwise,
        } => {
            let start_angle = (start.axial - center.axial).atan2(start.radius - center.radius);
            let end_angle = (end.axial - center.axial).atan2(end.radius - center.radius);
            let sweep = if clockwise {
                -(start_angle - end_angle).rem_euclid(TAU)
            } else {
                (end_angle - start_angle).rem_euclid(TAU)
            };
            let sweep_magnitude = sweep.abs();
            let middle_angle = start_angle + 0.5 * sweep;
            let arc_radius = (start.radius - center.radius).hypot(start.axial - center.axial);

            // Parameterizing by unsigned arc length avoids subtracting the
            // endpoint sines. The bracket is the exact integral of radius
            // over the admitted sweep and remains well conditioned for short
            // fillets and either orientation.
            let radial_integral = center.radius * sweep_magnitude
                + 2.0
                    * arc_radius
                    * fs_math::det::cos(middle_angle)
                    * fs_math::det::sin(0.5 * sweep_magnitude);
            TAU * arc_radius * radial_integral
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct BoundaryMoments {
    r2: f64,
    r2_z: f64,
    r4: f64,
    r2_abs: f64,
    r2_z_abs: f64,
    r4_abs: f64,
}

impl BoundaryMoments {
    fn add(&mut self, segment: MeridianSegment) -> Result<(), AxisymmetricMassError> {
        let r2 = boundary_integral(segment, 2, 0);
        let r2_z = boundary_integral(segment, 2, 1);
        let r4 = boundary_integral(segment, 4, 0);
        self.r2 += r2;
        self.r2_z += r2_z;
        self.r4 += r4;
        self.r2_abs += r2.abs();
        self.r2_z_abs += r2_z.abs();
        self.r4_abs += r4.abs();
        ensure_finite(self.r2, "volume boundary moment")?;
        ensure_finite(self.r2_z, "axial first boundary moment")?;
        ensure_finite(self.r4, "axial inertia boundary moment")?;
        ensure_finite(self.r2_abs, "volume absolute boundary moment")?;
        ensure_finite(self.r2_z_abs, "axial first absolute boundary moment")?;
        ensure_finite(self.r4_abs, "axial inertia absolute boundary moment")?;
        Ok(())
    }
}

/// Evaluate a boundary integral after translating the axial coordinate to the
/// known center of mass. This is a second pass on the exact features rather
/// than a parallel-axis subtraction of two large origin moments.
fn boundary_integral_centered_axial(
    segment: MeridianSegment,
    radius_power: u32,
    axial_power: u32,
    axial_center: f64,
) -> f64 {
    match segment {
        MeridianSegment::Line { start, end } => line_integral(
            start.radius,
            start.axial - axial_center,
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
            start.axial - axial_center,
            end.radius,
            end.axial - axial_center,
            center.radius,
            center.axial - axial_center,
            clockwise,
            radius_power,
            axial_power,
        ),
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
