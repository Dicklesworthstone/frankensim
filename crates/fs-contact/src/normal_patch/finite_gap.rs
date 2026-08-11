//! Finite-geometry frictionless normal contact on elastic half spaces.
//!
//! Unlike the local Hertz rungs in [`super::NormalPatchLaw`], this solver does
//! not replace the undeformed separation by a quadratic tangent.  The caller
//! supplies the complete relative gap sampled on a uniform tangent-plane grid.
//! Constant-pressure cells are coupled by the exact rectangular-cell integral
//! of the Boussinesq surface-compliance kernel, and a deterministic active-set
//! solve enforces
//!
//! `pressure >= 0`, `separation >= 0`, and `pressure * separation = 0`.
//!
//! This remains a two-elastic-half-space model.  A contact patch reaching the
//! sampled boundary is refused rather than truncated, and callers remain
//! responsible for half-space depth, finite-layer, yield, adhesion, thermal,
//! and geometry-sampling applicability.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Chart, Point3, TraceStepClaim, Vec3};
use fs_la::factor::{FactorError, cholesky};

/// Stable model identity for the finite-gap half-space discretization.
pub const FINITE_GAP_HALF_SPACE_MODEL_ID: &str =
    "org.frankensim.fs-contact.finite-gap-half-space-bem.v1";

/// Stable identity for certified chart-to-gap sampling.
pub const FINITE_GAP_CHART_SAMPLING_MODEL_ID: &str =
    "org.frankensim.fs-contact.finite-gap-chart-sampling.v1";

/// Public work bound for one contact grid.
pub const MAX_FINITE_GAP_CELLS: usize = 16_384;

/// Uniform cell-centred tangent-plane grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteGapGrid {
    /// Number of cells along local tangent direction `x`.
    pub cells_x: usize,
    /// Number of cells along local tangent direction `y`.
    pub cells_y: usize,
    /// Cell width along `x` [m].
    pub cell_width_m: f64,
    /// Cell width along `y` [m].
    pub cell_depth_m: f64,
}

impl FiniteGapGrid {
    fn cell_count(self) -> Result<usize, FiniteGapContactError> {
        let count =
            self.cells_x
                .checked_mul(self.cells_y)
                .ok_or(FiniteGapContactError::InvalidInput {
                    field: "grid cell count",
                })?;
        if self.cells_x < 3
            || self.cells_y < 3
            || count > MAX_FINITE_GAP_CELLS
            || !(self.cell_width_m.is_finite() && self.cell_width_m > 0.0)
            || !(self.cell_depth_m.is_finite() && self.cell_depth_m > 0.0)
        {
            return Err(FiniteGapContactError::InvalidInput {
                field: "finite-gap grid",
            });
        }
        Ok(count)
    }

    fn center(self, index: usize) -> (f64, f64) {
        let ix = index % self.cells_x;
        let iy = index / self.cells_x;
        (
            (ix as f64 + 0.5 - 0.5 * self.cells_x as f64) * self.cell_width_m,
            (iy as f64 + 0.5 - 0.5 * self.cells_y as f64) * self.cell_depth_m,
        )
    }
}

/// Orthonormal contact frame anchored at a nominal surface support point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteGapContactFrame {
    /// Nominal point on the undeformed body surface [m].
    pub surface_point_m: Point3,
    /// Unit normal pointing out of the body at `surface_point_m`.
    pub outward_normal: Vec3,
    /// First unit tangent direction.
    pub tangent_x: Vec3,
    /// Second unit tangent direction.
    pub tangent_y: Vec3,
}

/// Minimum chart evidence admitted by a chart-to-gap request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiniteGapChartEvidenceRequirement {
    /// Require a rigorous field enclosure for every sign decision.
    RigorousEnclosure,
    /// Permit finite nominal values when rigorous trace evidence is absent;
    /// the resulting receipt is explicitly downgraded to an estimate.
    AllowEstimate,
}

/// Authority actually achieved by chart-to-gap sampling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiniteGapChartSamplingAuthority {
    /// Every retained interval is a rigorous outside-to-inside root bracket.
    Enclosure,
    /// At least one sign or surface location used a nominal chart value.
    Estimate,
}

/// Bounded request to sample a solid chart into a finite contact gap.
pub struct FiniteGapChartSamplingRequest<'a> {
    /// Solid chart whose negative set is the contacting body.
    pub chart: &'a dyn Chart,
    /// Caller-supplied immutable geometry identity for provenance.
    pub source_geometry_id: &'a str,
    /// Cell-centred tangent-plane grid.
    pub grid: FiniteGapGrid,
    /// Contact frame and nominal support point.
    pub frame: FiniteGapContactFrame,
    /// Minimum evidence required from the source chart.
    pub evidence_requirement: FiniteGapChartEvidenceRequirement,
    /// Distance probed outward from the tangent plane to establish an
    /// outside endpoint for every cell [m].
    pub outside_probe_m: f64,
    /// Maximum distance searched inward from the tangent plane [m].
    pub maximum_inward_search_m: f64,
    /// Maximum retained root bracket width [m].
    pub root_tolerance_m: f64,
    /// Maximum bisection iterations per cell.
    pub maximum_bisection_steps: usize,
}

/// Certified cell-wise surface-height brackets sampled from a solid chart.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapChartSamplingReceipt {
    /// Stable identity of the complete request and retained brackets.
    pub identity: ContentHash,
    /// Caller-supplied source geometry identity.
    pub source_geometry_id: String,
    /// Chart kind reported by [`Chart::name`].
    pub chart_name: &'static str,
    /// Sample grid.
    pub grid: FiniteGapGrid,
    /// Authority achieved across every cell.
    pub authority: FiniteGapChartSamplingAuthority,
    /// Midpoint of every retained surface-height interval [m].
    pub undeformed_gap_m: Vec<f64>,
    /// Lower endpoint of every retained interval [m]. It is a rigorous bound
    /// only when `authority == Enclosure`.
    pub surface_height_lower_m: Vec<f64>,
    /// Upper endpoint of every retained interval [m]. It is a rigorous bound
    /// only when `authority == Enclosure`.
    pub surface_height_upper_m: Vec<f64>,
    /// Largest retained bracket width [m].
    pub maximum_surface_bracket_m: f64,
    /// Largest number of bisection steps used by any cell.
    pub maximum_steps_used: usize,
}

/// Typed refusal from certified chart-to-gap sampling.
#[derive(Debug, Clone, PartialEq)]
pub enum FiniteGapChartSamplingError {
    /// The frame, grid, provenance, search range, or work bound is invalid.
    InvalidInput { field: &'static str },
    /// The chart does not provide the continuity and sign evidence needed for
    /// a certified line root.
    UnsupportedTraceClaim { chart_name: &'static str },
    /// A chart evaluation did not return a finite rigorous enclosure that
    /// contains its nominal field value.
    UncertifiedSample {
        cell_x: usize,
        cell_y: usize,
        inward_parameter_m: f64,
    },
    /// The outward probe was not certified outside the solid.
    OutsideEndpointNotOutside { cell_x: usize, cell_y: usize },
    /// The inward search endpoint was not certified inside the solid.
    InsideEndpointNotInside { cell_x: usize, cell_y: usize },
    /// A field enclosure straddled zero before the positional bracket met the
    /// requested tolerance.
    AmbiguousSurfaceSample {
        cell_x: usize,
        cell_y: usize,
        bracket_width_m: f64,
    },
    /// The requested bisection work bound was insufficient.
    BisectionDidNotConverge {
        cell_x: usize,
        cell_y: usize,
        bracket_width_m: f64,
    },
    /// Cancellation was observed between deterministic cells or bisections.
    Cancelled,
}

impl fmt::Display for FiniteGapChartSamplingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { field } => {
                write!(f, "finite-gap chart sampling received invalid {field}")
            }
            Self::UnsupportedTraceClaim { chart_name } => write!(
                f,
                "chart {chart_name} does not provide a certified continuous trace claim"
            ),
            Self::UncertifiedSample {
                cell_x,
                cell_y,
                inward_parameter_m,
            } => write!(
                f,
                "chart sample at cell ({cell_x}, {cell_y}) and inward parameter {inward_parameter_m:.9e} m lacks a valid rigorous enclosure"
            ),
            Self::OutsideEndpointNotOutside { cell_x, cell_y } => write!(
                f,
                "outward endpoint at cell ({cell_x}, {cell_y}) is not certified outside the chart"
            ),
            Self::InsideEndpointNotInside { cell_x, cell_y } => write!(
                f,
                "inward endpoint at cell ({cell_x}, {cell_y}) is not certified inside the chart"
            ),
            Self::AmbiguousSurfaceSample {
                cell_x,
                cell_y,
                bracket_width_m,
            } => write!(
                f,
                "chart field at cell ({cell_x}, {cell_y}) straddles zero while the surface bracket remains {bracket_width_m:.9e} m wide"
            ),
            Self::BisectionDidNotConverge {
                cell_x,
                cell_y,
                bracket_width_m,
            } => write!(
                f,
                "surface bisection at cell ({cell_x}, {cell_y}) stopped with a {bracket_width_m:.9e} m bracket"
            ),
            Self::Cancelled => write!(f, "finite-gap chart sampling cancelled"),
        }
    }
}

impl core::error::Error for FiniteGapChartSamplingError {}

/// Sample the first outside-to-inside crossing of a solid chart along the
/// contact-frame inward normal at every tangent-grid cell.
///
/// The output retains lower and upper surface-height brackets.  The midpoint
/// is a convenient nominal gap for [`solve_finite_gap_half_space`], while the
/// brackets remain available for refinement and response sensitivity studies.
pub fn sample_finite_gap_from_chart(
    request: &FiniteGapChartSamplingRequest<'_>,
    cx: &Cx<'_>,
) -> Result<FiniteGapChartSamplingReceipt, FiniteGapChartSamplingError> {
    let cell_count =
        request
            .grid
            .cell_count()
            .map_err(|_| FiniteGapChartSamplingError::InvalidInput {
                field: "finite-gap grid",
            })?;
    if request.source_geometry_id.trim().is_empty()
        || !finite_point(request.frame.surface_point_m)
        || !orthonormal_contact_frame(request.frame)
        || !(request.outside_probe_m.is_finite() && request.outside_probe_m > 0.0)
        || !(request.maximum_inward_search_m.is_finite() && request.maximum_inward_search_m > 0.0)
        || !(request.root_tolerance_m.is_finite() && request.root_tolerance_m > 0.0)
        || request.maximum_bisection_steps == 0
    {
        return Err(FiniteGapChartSamplingError::InvalidInput {
            field: "chart sampling request",
        });
    }
    if request.evidence_requirement == FiniteGapChartEvidenceRequirement::RigorousEnclosure
        && request.chart.trace_step_claim() == TraceStepClaim::NoClaim
    {
        return Err(FiniteGapChartSamplingError::UnsupportedTraceClaim {
            chart_name: request.chart.name(),
        });
    }
    cx.checkpoint()
        .map_err(|_| FiniteGapChartSamplingError::Cancelled)?;

    let inward = request.frame.outward_normal.scale(-1.0);
    let mut midpoints = Vec::with_capacity(cell_count);
    let mut lower = Vec::with_capacity(cell_count);
    let mut upper = Vec::with_capacity(cell_count);
    let mut maximum_surface_bracket_m = 0.0_f64;
    let mut maximum_steps_used = 0_usize;
    let mut authority = FiniteGapChartSamplingAuthority::Enclosure;
    for index in 0..cell_count {
        cx.checkpoint()
            .map_err(|_| FiniteGapChartSamplingError::Cancelled)?;
        let cell_x = index % request.grid.cells_x;
        let cell_y = index / request.grid.cells_x;
        let (x, y) = request.grid.center(index);
        let plane_point = request
            .frame
            .surface_point_m
            .offset(request.frame.tangent_x.scale(x))
            .offset(request.frame.tangent_y.scale(y));
        let mut lo = -request.outside_probe_m;
        let mut hi = request.maximum_inward_search_m;
        let lo_observation = chart_sign(request, plane_point, inward, lo, cell_x, cell_y, cx)?;
        authority = weakest_chart_authority(authority, lo_observation.authority);
        if lo_observation.sign != ChartSign::Outside {
            return Err(FiniteGapChartSamplingError::OutsideEndpointNotOutside { cell_x, cell_y });
        }
        let hi_observation = chart_sign(request, plane_point, inward, hi, cell_x, cell_y, cx)?;
        authority = weakest_chart_authority(authority, hi_observation.authority);
        if hi_observation.sign != ChartSign::Inside {
            return Err(FiniteGapChartSamplingError::InsideEndpointNotInside { cell_x, cell_y });
        }

        let mut steps_used = 0_usize;
        while hi - lo > request.root_tolerance_m && steps_used < request.maximum_bisection_steps {
            cx.checkpoint()
                .map_err(|_| FiniteGapChartSamplingError::Cancelled)?;
            let mid = f64::midpoint(lo, hi);
            let observation = chart_sign(request, plane_point, inward, mid, cell_x, cell_y, cx)?;
            authority = weakest_chart_authority(authority, observation.authority);
            match observation.sign {
                ChartSign::Outside => lo = mid,
                ChartSign::Inside => hi = mid,
                ChartSign::ExactBoundary | ChartSign::NominalBoundary => {
                    lo = mid;
                    hi = mid;
                }
                ChartSign::BoundaryAmbiguous => {
                    if hi - lo > request.root_tolerance_m {
                        return Err(FiniteGapChartSamplingError::AmbiguousSurfaceSample {
                            cell_x,
                            cell_y,
                            bracket_width_m: hi - lo,
                        });
                    }
                }
            }
            steps_used += 1;
        }
        let width = hi - lo;
        if width > request.root_tolerance_m {
            return Err(FiniteGapChartSamplingError::BisectionDidNotConverge {
                cell_x,
                cell_y,
                bracket_width_m: width,
            });
        }
        maximum_surface_bracket_m = maximum_surface_bracket_m.max(width);
        maximum_steps_used = maximum_steps_used.max(steps_used);
        lower.push(lo);
        upper.push(hi);
        midpoints.push(f64::midpoint(lo, hi));
    }

    let identity = chart_sampling_identity(request, &lower, &upper);
    Ok(FiniteGapChartSamplingReceipt {
        identity,
        source_geometry_id: request.source_geometry_id.to_owned(),
        chart_name: request.chart.name(),
        grid: request.grid,
        authority,
        undeformed_gap_m: midpoints,
        surface_height_lower_m: lower,
        surface_height_upper_m: upper,
        maximum_surface_bracket_m,
        maximum_steps_used,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChartSign {
    Outside,
    Inside,
    ExactBoundary,
    NominalBoundary,
    BoundaryAmbiguous,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ChartSignObservation {
    sign: ChartSign,
    authority: FiniteGapChartSamplingAuthority,
}

fn chart_sign(
    request: &FiniteGapChartSamplingRequest<'_>,
    plane_point: Point3,
    inward: Vec3,
    parameter_m: f64,
    cell_x: usize,
    cell_y: usize,
    cx: &Cx<'_>,
) -> Result<ChartSignObservation, FiniteGapChartSamplingError> {
    let point = plane_point.offset(inward.scale(parameter_m));
    let sample = request.chart.eval(point, cx);
    let enclosure = request.chart.trace_value_enclosure(point, &sample, cx);
    if valid_rigorous_enclosure(enclosure, sample.signed_distance) {
        let sign =
            if enclosure.kind == NumericalKind::Exact && enclosure.lo == 0.0 && enclosure.hi == 0.0
            {
                ChartSign::ExactBoundary
            } else if enclosure.lo > 0.0 {
                ChartSign::Outside
            } else if enclosure.hi < 0.0 {
                ChartSign::Inside
            } else {
                ChartSign::BoundaryAmbiguous
            };
        return Ok(ChartSignObservation {
            sign,
            authority: FiniteGapChartSamplingAuthority::Enclosure,
        });
    }
    if request.evidence_requirement == FiniteGapChartEvidenceRequirement::AllowEstimate
        && sample.signed_distance.is_finite()
    {
        let sign = if sample.signed_distance > 0.0 {
            ChartSign::Outside
        } else if sample.signed_distance < 0.0 {
            ChartSign::Inside
        } else {
            ChartSign::NominalBoundary
        };
        Ok(ChartSignObservation {
            sign,
            authority: FiniteGapChartSamplingAuthority::Estimate,
        })
    } else {
        Err(FiniteGapChartSamplingError::UncertifiedSample {
            cell_x,
            cell_y,
            inward_parameter_m: parameter_m,
        })
    }
}

fn weakest_chart_authority(
    left: FiniteGapChartSamplingAuthority,
    right: FiniteGapChartSamplingAuthority,
) -> FiniteGapChartSamplingAuthority {
    if left == FiniteGapChartSamplingAuthority::Estimate
        || right == FiniteGapChartSamplingAuthority::Estimate
    {
        FiniteGapChartSamplingAuthority::Estimate
    } else {
        FiniteGapChartSamplingAuthority::Enclosure
    }
}

fn valid_rigorous_enclosure(enclosure: NumericalCertificate, nominal: f64) -> bool {
    matches!(
        enclosure.kind,
        NumericalKind::Exact | NumericalKind::Enclosure
    ) && nominal.is_finite()
        && enclosure.lo.is_finite()
        && enclosure.hi.is_finite()
        && enclosure.lo <= nominal
        && nominal <= enclosure.hi
}

fn finite_point(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn orthonormal_contact_frame(frame: FiniteGapContactFrame) -> bool {
    let vectors = [frame.outward_normal, frame.tangent_x, frame.tangent_y];
    vectors
        .iter()
        .all(|vector| (vector.norm() - 1.0).abs() <= 2.0e-10)
        && frame.outward_normal.dot(frame.tangent_x).abs() <= 2.0e-10
        && frame.outward_normal.dot(frame.tangent_y).abs() <= 2.0e-10
        && frame.tangent_x.dot(frame.tangent_y).abs() <= 2.0e-10
}

fn chart_sampling_identity(
    request: &FiniteGapChartSamplingRequest<'_>,
    lower: &[f64],
    upper: &[f64],
) -> ContentHash {
    let mut bytes = Vec::with_capacity(256 + 16 * lower.len());
    bytes.extend_from_slice(&(request.source_geometry_id.len() as u64).to_le_bytes());
    bytes.extend_from_slice(request.source_geometry_id.as_bytes());
    bytes.extend_from_slice(&(request.chart.name().len() as u64).to_le_bytes());
    bytes.extend_from_slice(request.chart.name().as_bytes());
    bytes.push(match request.evidence_requirement {
        FiniteGapChartEvidenceRequirement::RigorousEnclosure => 0,
        FiniteGapChartEvidenceRequirement::AllowEstimate => 1,
    });
    bytes.extend_from_slice(&(request.grid.cells_x as u64).to_le_bytes());
    bytes.extend_from_slice(&(request.grid.cells_y as u64).to_le_bytes());
    bytes.extend_from_slice(&request.grid.cell_width_m.to_bits().to_le_bytes());
    bytes.extend_from_slice(&request.grid.cell_depth_m.to_bits().to_le_bytes());
    for value in [
        request.frame.surface_point_m.x,
        request.frame.surface_point_m.y,
        request.frame.surface_point_m.z,
        request.frame.outward_normal.x,
        request.frame.outward_normal.y,
        request.frame.outward_normal.z,
        request.frame.tangent_x.x,
        request.frame.tangent_x.y,
        request.frame.tangent_x.z,
        request.frame.tangent_y.x,
        request.frame.tangent_y.y,
        request.frame.tangent_y.z,
        request.outside_probe_m,
        request.maximum_inward_search_m,
        request.root_tolerance_m,
    ] {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(&(request.maximum_bisection_steps as u64).to_le_bytes());
    for (&lo, &hi) in lower.iter().zip(upper) {
        bytes.extend_from_slice(&lo.to_bits().to_le_bytes());
        bytes.extend_from_slice(&hi.to_bits().to_le_bytes());
    }
    hash_domain(FINITE_GAP_CHART_SAMPLING_MODEL_ID, &bytes)
}

/// One bounded prescribed-approach contact request.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapContactRequest {
    /// Uniform tangent-plane grid.
    pub grid: FiniteGapGrid,
    /// Undeformed relative separation at cell centres [m].
    ///
    /// Its finite minimum is treated as the zero-separation datum.  This makes
    /// the request invariant to a constant height offset without silently
    /// changing `approach_m`.
    pub undeformed_gap_m: Vec<f64>,
    /// Rigid relative approach from the minimum undeformed gap [m].
    pub approach_m: f64,
    /// Classical reduced modulus `E*` of the two contacting bodies [Pa].
    pub reduced_modulus_pa: f64,
    /// Maximum deterministic active-set changes.
    pub maximum_active_set_iterations: usize,
    /// Absolute separation/complementarity tolerance [m].
    pub complementarity_tolerance_m: f64,
    /// Required ring of inactive cells between the patch and every grid edge.
    pub boundary_clearance_cells: usize,
}

/// Converged pressure and integral response for one finite gap.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapContactReceipt {
    /// Stable identity of the complete request and retained pressure field.
    pub identity: ContentHash,
    /// Pressure at every cell centre under the constant-cell basis [Pa].
    pub pressure_pa: Vec<f64>,
    /// Relative elastic surface displacement at every cell centre [m].
    pub elastic_opening_m: Vec<f64>,
    /// Number of positive-pressure cells.
    pub active_cells: usize,
    /// Active-set changes required for convergence.
    pub active_set_iterations: usize,
    /// Integrated compressive resultant [N].
    pub normal_force_n: f64,
    /// Maximum cell pressure [Pa].
    pub peak_pressure_pa: f64,
    /// Pressure centroid in the supplied tangent frame [m].
    pub pressure_centroid_m: [f64; 2],
    /// Equivalent Hertz semiaxes reconstructed from pressure second moments
    /// (`<x^2>=a^2/5`, `<y^2>=b^2/5`) [m].
    pub equivalent_pressure_semiaxes_m: [f64; 2],
    /// Reversible half-space strain energy `1/2 integral(p u dA)` [J].
    pub reversible_energy_j: f64,
    /// Maximum active closure or inactive penetration residual [m].
    pub complementarity_residual_m: f64,
    /// Constant height datum subtracted from the supplied gap [m].
    pub undeformed_gap_datum_m: f64,
}

/// Typed refusal from the finite-gap half-space solve.
#[derive(Debug, Clone, PartialEq)]
pub enum FiniteGapContactError {
    /// A scalar, shape, tolerance, or work bound is invalid.
    InvalidInput { field: &'static str },
    /// The reduced compliance matrix failed its SPD factorization.
    ComplianceNotPositiveDefinite { active_index: usize },
    /// The pressure active set exceeded the caller's iteration budget.
    ActiveSetDidNotConverge { iterations: usize, residual_m: f64 },
    /// Positive pressure reaches the forbidden boundary ring, so truncation
    /// error is not bounded by this request.
    ContactReachedDomainBoundary {
        cell_x: usize,
        cell_y: usize,
        required_clearance_cells: usize,
    },
    /// Cancellation was observed between deterministic active-set iterations.
    Cancelled,
}

impl fmt::Display for FiniteGapContactError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { field } => {
                write!(f, "finite-gap half-space contact received invalid {field}")
            }
            Self::ComplianceNotPositiveDefinite { active_index } => write!(
                f,
                "finite-gap half-space compliance lost positive definiteness at active index {active_index}"
            ),
            Self::ActiveSetDidNotConverge {
                iterations,
                residual_m,
            } => write!(
                f,
                "finite-gap contact active set did not converge in {iterations} changes; residual {residual_m:.9e} m"
            ),
            Self::ContactReachedDomainBoundary {
                cell_x,
                cell_y,
                required_clearance_cells,
            } => write!(
                f,
                "finite-gap contact reaches grid cell ({cell_x}, {cell_y}); at least {required_clearance_cells} inactive boundary cells are required"
            ),
            Self::Cancelled => write!(f, "finite-gap contact solve cancelled"),
        }
    }
}

impl core::error::Error for FiniteGapContactError {}

/// Solve frictionless finite-geometry contact for a prescribed approach.
///
/// The raw pressure field is useful both as a constitutive response and for
/// building a profile/material-specific interpolation table.  Repeated hot
/// dynamics steps should consume such an independently convergence-checked
/// table rather than rerun this dense reference solve at audio rates.
pub fn solve_finite_gap_half_space(
    request: &FiniteGapContactRequest,
    cx: &Cx<'_>,
) -> Result<FiniteGapContactReceipt, FiniteGapContactError> {
    let cell_count = request.grid.cell_count()?;
    if request.undeformed_gap_m.len() != cell_count {
        return Err(FiniteGapContactError::InvalidInput {
            field: "undeformed gap shape",
        });
    }
    if request
        .undeformed_gap_m
        .iter()
        .any(|value| !value.is_finite())
        || !(request.approach_m.is_finite() && request.approach_m >= 0.0)
        || !(request.reduced_modulus_pa.is_finite() && request.reduced_modulus_pa > 0.0)
        || request.maximum_active_set_iterations == 0
        || !(request.complementarity_tolerance_m.is_finite()
            && request.complementarity_tolerance_m > 0.0)
        || request.boundary_clearance_cells == 0
        || 2 * request.boundary_clearance_cells >= request.grid.cells_x
        || 2 * request.boundary_clearance_cells >= request.grid.cells_y
    {
        return Err(FiniteGapContactError::InvalidInput {
            field: "finite-gap request",
        });
    }
    cx.checkpoint()
        .map_err(|_| FiniteGapContactError::Cancelled)?;

    let gap_datum = request
        .undeformed_gap_m
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let gap: Vec<f64> = request
        .undeformed_gap_m
        .iter()
        .map(|value| value - gap_datum)
        .collect();
    let mut active: Vec<usize> = gap
        .iter()
        .enumerate()
        .filter_map(|(index, value)| (*value < request.approach_m).then_some(index))
        .collect();
    let mut pressure = vec![0.0; cell_count];
    let mut opening = vec![0.0; cell_count];
    if active.is_empty() {
        return build_receipt(request, gap_datum, pressure, opening, 0, 0.0);
    }

    // The Boussinesq cell kernel is translation invariant on the uniform
    // grid, so the pairwise compliance depends only on the cell index
    // difference. Precomputing that (2 cells_x - 1) x (2 cells_y - 1)
    // difference table once turns every later matrix/opening entry into a
    // lookup instead of a fresh hypot + two logarithms per pair per
    // active-set iteration (the previous per-pair evaluation dominated the
    // solve for large patches).
    let kernel = ComplianceKernel::build(request.grid, request.reduced_modulus_pa);
    let mut last_residual = f64::INFINITY;
    for iteration in 0..request.maximum_active_set_iterations {
        cx.checkpoint()
            .map_err(|_| FiniteGapContactError::Cancelled)?;
        active.sort_unstable();
        active.dedup();
        let n = active.len();
        let mut compliance = vec![0.0; n * n];
        for (row, &receiver) in active.iter().enumerate() {
            for (column, &source) in active.iter().enumerate().take(row + 1) {
                let value = kernel.compliance(receiver, source);
                compliance[row * n + column] = value;
                compliance[column * n + row] = value;
            }
        }
        let factor = cholesky(&compliance, n).map_err(|error| match error {
            FactorError::NotSpd { index } | FactorError::Singular { index } => {
                FiniteGapContactError::ComplianceNotPositiveDefinite {
                    active_index: index,
                }
            }
        })?;
        let mut active_pressure: Vec<f64> = active
            .iter()
            .map(|&index| request.approach_m - gap[index])
            .collect();
        factor.solve(&mut active_pressure);

        let pressure_tolerance_pa = request.reduced_modulus_pa
            * request.complementarity_tolerance_m
            / request.grid.cell_width_m.min(request.grid.cell_depth_m);
        let retained_active: Vec<usize> = active
            .iter()
            .copied()
            .zip(active_pressure.iter().copied())
            .filter_map(|(index, value)| (value >= -pressure_tolerance_pa).then_some(index))
            .collect();
        if retained_active.len() != active.len() {
            // The unconstrained solve can identify many tensile cells at
            // once for a long, anisotropic contact patch. Removing the whole
            // deterministic inadmissible set is the same active-set pivot as
            // one-at-a-time deletion, while avoiding an iteration count that
            // scales with the initial geometric candidate count. Any cell
            // made compressive by the block pivot is still rediscovered by
            // the inactive-separation check below before publication.
            active = retained_active;
            if active.is_empty() {
                pressure.fill(0.0);
                opening.fill(0.0);
                return build_receipt(request, gap_datum, pressure, opening, iteration + 1, 0.0);
            }
            continue;
        }

        pressure.fill(0.0);
        for (&index, value) in active.iter().zip(active_pressure) {
            pressure[index] = value.max(0.0);
        }
        opening.fill(0.0);
        for (receiver, displacement) in opening.iter_mut().enumerate() {
            let mut value = 0.0;
            for &source in &active {
                value = kernel
                    .compliance(receiver, source)
                    .mul_add(pressure[source], value);
            }
            *displacement = value;
        }

        let separation = |index: usize| gap[index] - request.approach_m + opening[index];
        let mut worst_inactive = None;
        let mut active_cursor = 0_usize;
        last_residual = 0.0;
        for index in 0..cell_count {
            let is_active = active_cursor < active.len() && active[active_cursor] == index;
            if is_active {
                active_cursor += 1;
                last_residual = last_residual.max(separation(index).abs());
            } else {
                let value = separation(index);
                if value < 0.0 {
                    last_residual = last_residual.max(-value);
                    if worst_inactive.is_none_or(|(_, worst): (usize, f64)| value < worst) {
                        worst_inactive = Some((index, value));
                    }
                }
            }
        }
        if let Some((index, value)) = worst_inactive
            && value < -request.complementarity_tolerance_m
        {
            active.push(index);
            continue;
        }
        verify_boundary_clearance(request.grid, &active, request.boundary_clearance_cells)?;
        return build_receipt(
            request,
            gap_datum,
            pressure,
            opening,
            iteration + 1,
            last_residual,
        );
    }
    Err(FiniteGapContactError::ActiveSetDidNotConverge {
        iterations: request.maximum_active_set_iterations,
        residual_m: last_residual,
    })
}

/// Translation-invariant cell-pair compliance, tabulated once per solve
/// over every possible index difference. Entries are the exact
/// rectangular-cell Boussinesq integrals evaluated at the exact
/// difference offsets `(di * width, dj * depth)`, so the assembled matrix
/// is exactly (bitwise) Toeplitz-symmetric.
struct ComplianceKernel {
    cells_x: usize,
    cells_y: usize,
    stride: usize,
    values: Vec<f64>,
}

impl ComplianceKernel {
    fn build(grid: FiniteGapGrid, reduced_modulus_pa: f64) -> Self {
        let stride = 2 * grid.cells_x - 1;
        let rows = 2 * grid.cells_y - 1;
        let mut values = Vec::with_capacity(stride * rows);
        let scale = 1.0 / (core::f64::consts::PI * reduced_modulus_pa);
        for row in 0..rows {
            #[allow(clippy::cast_precision_loss)]
            let dj = row as f64 - (grid.cells_y - 1) as f64;
            for column in 0..stride {
                #[allow(clippy::cast_precision_loss)]
                let di = column as f64 - (grid.cells_x - 1) as f64;
                let integral = rectangle_inverse_distance_integral(
                    di * grid.cell_width_m,
                    dj * grid.cell_depth_m,
                    0.5 * grid.cell_width_m,
                    0.5 * grid.cell_depth_m,
                );
                values.push(integral * scale);
            }
        }
        Self {
            cells_x: grid.cells_x,
            cells_y: grid.cells_y,
            stride,
            values,
        }
    }

    fn compliance(&self, receiver: usize, source: usize) -> f64 {
        let di = (source % self.cells_x) + self.cells_x - 1 - (receiver % self.cells_x);
        let dj = (source / self.cells_x) + self.cells_y - 1 - (receiver / self.cells_x);
        self.values[dj * self.stride + di]
    }
}

/// Exact integral of `1/r` over a rectangle centred at `(x, y)`.
fn rectangle_inverse_distance_integral(x: f64, y: f64, hx: f64, hy: f64) -> f64 {
    let primitive = |a: f64, b: f64| {
        let radius = a.hypot(b);
        let mut value = 0.0;
        if a != 0.0 {
            value += a * (b + radius).ln();
        }
        if b != 0.0 {
            value += b * (a + radius).ln();
        }
        value
    };
    primitive(x + hx, y + hy) - primitive(x - hx, y + hy) - primitive(x + hx, y - hy)
        + primitive(x - hx, y - hy)
}

fn verify_boundary_clearance(
    grid: FiniteGapGrid,
    active: &[usize],
    clearance: usize,
) -> Result<(), FiniteGapContactError> {
    for &index in active {
        let x = index % grid.cells_x;
        let y = index / grid.cells_x;
        if x < clearance
            || y < clearance
            || x >= grid.cells_x - clearance
            || y >= grid.cells_y - clearance
        {
            return Err(FiniteGapContactError::ContactReachedDomainBoundary {
                cell_x: x,
                cell_y: y,
                required_clearance_cells: clearance,
            });
        }
    }
    Ok(())
}

fn build_receipt(
    request: &FiniteGapContactRequest,
    gap_datum: f64,
    pressure: Vec<f64>,
    opening: Vec<f64>,
    iterations: usize,
    residual_m: f64,
) -> Result<FiniteGapContactReceipt, FiniteGapContactError> {
    let area = request.grid.cell_width_m * request.grid.cell_depth_m;
    let normal_force_n = pressure.iter().sum::<f64>() * area;
    let peak_pressure_pa = pressure.iter().copied().fold(0.0_f64, f64::max);
    let active_cells = pressure.iter().filter(|value| **value > 0.0).count();
    let mut centroid = [0.0; 2];
    if normal_force_n > 0.0 {
        for (index, &value) in pressure.iter().enumerate() {
            let force = value * area;
            let (x, y) = request.grid.center(index);
            centroid[0] = x.mul_add(force, centroid[0]);
            centroid[1] = y.mul_add(force, centroid[1]);
        }
        centroid[0] /= normal_force_n;
        centroid[1] /= normal_force_n;
    }
    let mut moments = [0.0; 2];
    if normal_force_n > 0.0 {
        let cell_variance = [
            request.grid.cell_width_m * request.grid.cell_width_m / 12.0,
            request.grid.cell_depth_m * request.grid.cell_depth_m / 12.0,
        ];
        for (index, &value) in pressure.iter().enumerate() {
            let force = value * area;
            let (x, y) = request.grid.center(index);
            let dx = x - centroid[0];
            let dy = y - centroid[1];
            moments[0] += force * (dx * dx + cell_variance[0]);
            moments[1] += force * (dy * dy + cell_variance[1]);
        }
        moments[0] /= normal_force_n;
        moments[1] /= normal_force_n;
    }
    let reversible_energy_j = 0.5
        * area
        * pressure
            .iter()
            .zip(&opening)
            .map(|(load, displacement)| load * displacement)
            .sum::<f64>();
    let mut identity_bytes = Vec::with_capacity(192 + 24 * pressure.len());
    identity_bytes.extend_from_slice(&(request.grid.cells_x as u64).to_le_bytes());
    identity_bytes.extend_from_slice(&(request.grid.cells_y as u64).to_le_bytes());
    identity_bytes.extend_from_slice(&request.grid.cell_width_m.to_bits().to_le_bytes());
    identity_bytes.extend_from_slice(&request.grid.cell_depth_m.to_bits().to_le_bytes());
    identity_bytes.extend_from_slice(&request.approach_m.to_bits().to_le_bytes());
    identity_bytes.extend_from_slice(&request.reduced_modulus_pa.to_bits().to_le_bytes());
    identity_bytes.extend_from_slice(&(request.maximum_active_set_iterations as u64).to_le_bytes());
    identity_bytes.extend_from_slice(&request.complementarity_tolerance_m.to_bits().to_le_bytes());
    identity_bytes.extend_from_slice(&(request.boundary_clearance_cells as u64).to_le_bytes());
    identity_bytes.extend_from_slice(&gap_datum.to_bits().to_le_bytes());
    for value in &request.undeformed_gap_m {
        identity_bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    for value in &pressure {
        identity_bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    for value in &opening {
        identity_bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    identity_bytes.extend_from_slice(&(iterations as u64).to_le_bytes());
    identity_bytes.extend_from_slice(&residual_m.to_bits().to_le_bytes());
    let identity = hash_domain(FINITE_GAP_HALF_SPACE_MODEL_ID, &identity_bytes);
    Ok(FiniteGapContactReceipt {
        identity,
        pressure_pa: pressure,
        elastic_opening_m: opening,
        active_cells,
        active_set_iterations: iterations,
        normal_force_n,
        peak_pressure_pa,
        pressure_centroid_m: centroid,
        equivalent_pressure_semiaxes_m: [(5.0 * moments[0]).sqrt(), (5.0 * moments[1]).sqrt()],
        reversible_energy_j,
        complementarity_residual_m: residual_m,
        undeformed_gap_datum_m: gap_datum,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_geom::fixtures::SphereChart;

    struct NominalOnlySphere(SphereChart);

    impl Chart for NominalOnlySphere {
        fn eval(&self, point: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
            self.0.eval(point, cx)
        }

        fn support(&self) -> fs_geom::Aabb {
            self.0.support()
        }

        fn name(&self) -> &'static str {
            "test/nominal-only-sphere"
        }
    }

    fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 11,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn paraboloid_request(cells: usize, cell_m: f64) -> FiniteGapContactRequest {
        let grid = FiniteGapGrid {
            cells_x: cells,
            cells_y: cells,
            cell_width_m: cell_m,
            cell_depth_m: cell_m,
        };
        let radius_m = 0.020;
        let mut gap = Vec::with_capacity(cells * cells);
        for index in 0..cells * cells {
            let (x, y) = grid.center(index);
            gap.push((x * x + y * y) / (2.0 * radius_m));
        }
        FiniteGapContactRequest {
            grid,
            undeformed_gap_m: gap,
            approach_m: 2.0e-6,
            reduced_modulus_pa: 80.0e9,
            maximum_active_set_iterations: 256,
            complementarity_tolerance_m: 1.0e-11,
            boundary_clearance_cells: 2,
        }
    }

    #[test]
    fn rectangle_self_integral_matches_closed_form() {
        let hx = 0.3_f64;
        let hy = 0.2_f64;
        let expected = 4.0 * (hx * (hy / hx).asinh() + hy * (hx / hy).asinh());
        let observed = rectangle_inverse_distance_integral(0.0, 0.0, hx, hy);
        assert!((observed - expected).abs() < 1.0e-14);
    }

    #[test]
    fn certified_sphere_chart_sampling_brackets_exact_surface_height() {
        let radius_m = 0.020_f64;
        let chart = SphereChart {
            center: Point3::new(0.0, 0.0, radius_m),
            radius: radius_m,
        };
        let request = FiniteGapChartSamplingRequest {
            chart: &chart,
            source_geometry_id: "test/sphere/r20mm/v1",
            grid: FiniteGapGrid {
                cells_x: 7,
                cells_y: 7,
                cell_width_m: 1.0e-3,
                cell_depth_m: 1.0e-3,
            },
            frame: FiniteGapContactFrame {
                surface_point_m: Point3::new(0.0, 0.0, 0.0),
                outward_normal: Vec3::new(0.0, 0.0, -1.0),
                tangent_x: Vec3::new(1.0, 0.0, 0.0),
                tangent_y: Vec3::new(0.0, 1.0, 0.0),
            },
            evidence_requirement: FiniteGapChartEvidenceRequirement::RigorousEnclosure,
            outside_probe_m: 1.0e-5,
            maximum_inward_search_m: 2.0e-3,
            root_tolerance_m: 1.0e-9,
            maximum_bisection_steps: 64,
        };
        let receipt = with_cx(|cx| sample_finite_gap_from_chart(&request, cx)).unwrap();
        assert!(receipt.maximum_surface_bracket_m <= request.root_tolerance_m);
        for index in 0..request.grid.cells_x * request.grid.cells_y {
            let (x, y) = request.grid.center(index);
            let expected = radius_m - (radius_m * radius_m - x * x - y * y).sqrt();
            assert!(
                receipt.surface_height_lower_m[index] <= expected
                    && expected <= receipt.surface_height_upper_m[index],
                "cell={index} bracket=[{},{}] expected={expected}",
                receipt.surface_height_lower_m[index],
                receipt.surface_height_upper_m[index]
            );
        }
    }

    #[test]
    fn nominal_only_chart_is_explicitly_downgraded_to_estimate() {
        let radius_m = 0.020_f64;
        let chart = NominalOnlySphere(SphereChart {
            center: Point3::new(0.0, 0.0, radius_m),
            radius: radius_m,
        });
        let grid = FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 1.0e-3,
            cell_depth_m: 1.0e-3,
        };
        let request = FiniteGapChartSamplingRequest {
            chart: &chart,
            source_geometry_id: "test/nominal-sphere/r20mm/v1",
            grid,
            frame: FiniteGapContactFrame {
                surface_point_m: Point3::new(0.0, 0.0, 0.0),
                outward_normal: Vec3::new(0.0, 0.0, -1.0),
                tangent_x: Vec3::new(1.0, 0.0, 0.0),
                tangent_y: Vec3::new(0.0, 1.0, 0.0),
            },
            evidence_requirement: FiniteGapChartEvidenceRequirement::AllowEstimate,
            outside_probe_m: 1.0e-5,
            maximum_inward_search_m: 2.0e-3,
            root_tolerance_m: 1.0e-9,
            maximum_bisection_steps: 64,
        };
        let receipt = with_cx(|cx| sample_finite_gap_from_chart(&request, cx)).unwrap();
        assert_eq!(receipt.authority, FiniteGapChartSamplingAuthority::Estimate);
        for (index, &observed) in receipt.undeformed_gap_m.iter().enumerate() {
            let (x, y) = grid.center(index);
            let expected = radius_m - (radius_m * radius_m - x * x - y * y).sqrt();
            assert!((observed - expected).abs() <= request.root_tolerance_m);
        }
    }

    #[test]
    fn g1_paraboloid_recovers_hertz_force_and_footprint() {
        let request = paraboloid_request(31, 2.5e-5);
        let receipt = with_cx(|cx| solve_finite_gap_half_space(&request, cx)).unwrap();
        let radius_m = 0.020_f64;
        let expected_force =
            4.0 / 3.0 * request.reduced_modulus_pa * radius_m.sqrt() * request.approach_m.powf(1.5);
        let expected_radius = (radius_m * request.approach_m).sqrt();
        assert!(
            (receipt.normal_force_n / expected_force - 1.0).abs() < 0.06,
            "force={} expected={expected_force}",
            receipt.normal_force_n
        );
        for axis in receipt.equivalent_pressure_semiaxes_m {
            assert!(
                (axis / expected_radius - 1.0).abs() < 0.08,
                "axis={axis} expected={expected_radius}"
            );
        }
        assert!(receipt.complementarity_residual_m <= request.complementarity_tolerance_m);
        assert!(receipt.reversible_energy_j > 0.0);
    }

    #[test]
    fn insufficient_domain_refuses_instead_of_truncating_pressure() {
        let mut request = paraboloid_request(9, 5.0e-5);
        request.approach_m = 8.0e-6;
        let error = with_cx(|cx| solve_finite_gap_half_space(&request, cx)).unwrap_err();
        assert!(matches!(
            error,
            FiniteGapContactError::ContactReachedDomainBoundary { .. }
        ));
    }
}
