//! Passive runtime response curves generated from finite-gap contact solves.
//!
//! The dense half-space solve is a convergence reference, not an audio-rate
//! kernel. This module evaluates it at caller-selected approach nodes and
//! constructs a monotone piecewise-linear force law. Integrating that same
//! line segment supplies a force/work-consistent stored-energy surrogate;
//! disagreement with the independently computed BEM energy is an explicit
//! table-refinement residual rather than hidden interpolation error.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;

use super::{
    FiniteGapChartSamplingAuthority, FiniteGapContactError, FiniteGapContactRequest, FiniteGapGrid,
    MAX_FINITE_GAP_CELLS, solve_finite_gap_half_space,
};

/// Stable identity for finite-gap response curves.
pub const FINITE_GAP_RESPONSE_CURVE_MODEL_ID: &str =
    "org.frankensim.fs-contact.finite-gap-response-curve.v1";

/// Stable identity for scalar-configuration families of response curves.
pub const FINITE_GAP_RESPONSE_FAMILY_MODEL_ID: &str =
    "org.frankensim.fs-contact.finite-gap-response-family.v1";

/// Maximum number of approach nodes admitted by one response curve.
pub const MAX_FINITE_GAP_RESPONSE_NODES: usize = 4_096;

/// Numerical authority of an evaluated response curve.
///
/// Dense nodes may originate from enclosed geometry, but the current linear
/// interpolation has no interval remainder bound between nodes. It therefore
/// remains an estimate until a refinement study or certified interpolant
/// supplies such a bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiniteGapResponseCurveAuthority {
    /// Piecewise interpolation backed by retained reference solves.
    Estimate,
}

/// Request to build one passive response curve from a fixed undeformed gap.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapResponseCurveRequest {
    /// Identity of the chart-sampling or equivalent geometry artifact.
    pub gap_source_identity: ContentHash,
    /// Authority of the supplied gap field.
    pub gap_source_authority: FiniteGapChartSamplingAuthority,
    /// Uniform tangent-plane grid.
    pub grid: FiniteGapGrid,
    /// Undeformed relative gap at cell centres [m].
    pub undeformed_gap_m: Vec<f64>,
    /// Classical two-body reduced modulus [Pa].
    pub reduced_modulus_pa: f64,
    /// Strictly increasing approach nodes beginning at exactly zero [m].
    pub approach_nodes_m: Vec<f64>,
    /// Active-set work bound forwarded to every dense solve.
    pub maximum_active_set_iterations: usize,
    /// Complementarity tolerance forwarded to every dense solve [m].
    pub complementarity_tolerance_m: f64,
    /// Required inactive boundary-cell ring forwarded to every dense solve.
    pub boundary_clearance_cells: usize,
    /// Absolute allowed mismatch between integrated curve work and the BEM
    /// reversible energy at every node [J].
    pub absolute_energy_tolerance_j: f64,
    /// Relative allowed work/energy mismatch at every node.
    pub relative_energy_tolerance: f64,
}

/// One retained dense-solve node in a response curve.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapResponseNode {
    /// Rigid approach [m].
    pub approach_m: f64,
    /// Identity of the complete dense contact receipt.
    pub contact_identity: ContentHash,
    /// Integrated normal force [N].
    pub normal_force_n: f64,
    /// Maximum cell pressure [Pa].
    pub peak_pressure_pa: f64,
    /// Pressure centroid in the tangent frame [m].
    pub pressure_centroid_m: [f64; 2],
    /// Equivalent second-moment pressure semiaxes [m].
    pub equivalent_pressure_semiaxes_m: [f64; 2],
    /// Independently computed half-space strain energy [J].
    pub reference_reversible_energy_j: f64,
    /// Integral of the retained piecewise-linear force curve through this node
    /// [J].
    pub curve_reversible_energy_j: f64,
    /// Signed curve-minus-reference energy residual [J].
    pub energy_residual_j: f64,
}

/// Immutable passive runtime response curve.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapResponseCurve {
    /// Stable identity of the source gap, build controls, and every node.
    pub identity: ContentHash,
    /// Identity of the source gap artifact.
    pub gap_source_identity: ContentHash,
    /// Authority of the source geometry sampling, retained without promotion.
    pub gap_source_authority: FiniteGapChartSamplingAuthority,
    /// Authority of values evaluated between retained dense-solve nodes.
    pub response_authority: FiniteGapResponseCurveAuthority,
    /// Retained approach/response nodes.
    pub nodes: Vec<FiniteGapResponseNode>,
    /// Largest absolute curve/reference energy residual [J].
    pub maximum_absolute_energy_residual_j: f64,
}

/// One force/work-consistent interpolation result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteGapInterpolatedResponse {
    /// Requested rigid approach [m].
    pub approach_m: f64,
    /// Piecewise-linear compressive force [N].
    pub normal_force_n: f64,
    /// Segment tangent `dF/d(delta)` [N/m].
    pub normal_tangent_n_per_m: f64,
    /// Integral of the interpolated force from zero approach [J].
    pub reversible_energy_j: f64,
    /// Linearly interpolated peak cell pressure [Pa].
    pub peak_pressure_pa: f64,
    /// Linearly interpolated pressure centroid [m].
    pub pressure_centroid_m: [f64; 2],
    /// Linearly interpolated equivalent pressure semiaxes [m].
    pub equivalent_pressure_semiaxes_m: [f64; 2],
}

/// A common-grid family of passive response curves indexed by one scalar
/// configuration coordinate such as inclination or temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteGapResponseFamily {
    /// Stable identity of the coordinate definition, nodes, and member curves.
    pub identity: ContentHash,
    /// Semantic coordinate identifier, for example `inclination`.
    pub coordinate_id: String,
    /// Explicit SI unit spelling, for example `rad` or `K`.
    pub coordinate_unit: String,
    /// Strictly increasing configuration nodes.
    pub coordinate_nodes: Vec<f64>,
    /// One response curve per coordinate node on an identical approach grid.
    pub curves: Vec<FiniteGapResponseCurve>,
    /// Authority of values interpolated in approach and configuration.
    pub response_authority: FiniteGapResponseCurveAuthority,
}

/// Typed refusal from response-curve construction or evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum FiniteGapResponseCurveError {
    /// A source field, node set, numerical control, or tolerance is invalid.
    InvalidInput { field: &'static str },
    /// A dense contact solve refused at one approach node.
    Contact {
        node: usize,
        source: FiniteGapContactError,
    },
    /// The dense normal force decreased between successive approaches.
    NonMonotoneForce {
        node: usize,
        previous_force_n: f64,
        force_n: f64,
    },
    /// The force-integral and independently computed BEM energy disagree
    /// beyond the caller's declared table error budget.
    EnergyConsistency {
        node: usize,
        residual_j: f64,
        allowed_j: f64,
    },
    /// Runtime evaluation attempted to extrapolate beyond the built domain.
    ApproachOutsideTable {
        approach_m: f64,
        minimum_m: f64,
        maximum_m: f64,
    },
    /// Family evaluation attempted to extrapolate its configuration coordinate.
    CoordinateOutsideFamily {
        coordinate: f64,
        minimum: f64,
        maximum: f64,
    },
}

impl fmt::Display for FiniteGapResponseCurveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { field } => {
                write!(f, "finite-gap response curve received invalid {field}")
            }
            Self::Contact { node, source } => {
                write!(f, "finite-gap response node {node} refused: {source}")
            }
            Self::NonMonotoneForce {
                node,
                previous_force_n,
                force_n,
            } => write!(
                f,
                "finite-gap response force decreased at node {node}: {previous_force_n:.9e} -> {force_n:.9e} N"
            ),
            Self::EnergyConsistency {
                node,
                residual_j,
                allowed_j,
            } => write!(
                f,
                "finite-gap response node {node} has work/energy residual {residual_j:.9e} J beyond {allowed_j:.9e} J"
            ),
            Self::ApproachOutsideTable {
                approach_m,
                minimum_m,
                maximum_m,
            } => write!(
                f,
                "finite-gap response approach {approach_m:.9e} m lies outside [{minimum_m:.9e}, {maximum_m:.9e}] m"
            ),
            Self::CoordinateOutsideFamily {
                coordinate,
                minimum,
                maximum,
            } => write!(
                f,
                "finite-gap response coordinate {coordinate:.9e} lies outside [{minimum:.9e}, {maximum:.9e}]"
            ),
        }
    }
}

impl core::error::Error for FiniteGapResponseCurveError {}

/// Build one response curve by running the finite-gap reference solver at
/// every declared approach node.
pub fn build_finite_gap_response_curve(
    request: &FiniteGapResponseCurveRequest,
    cx: &Cx<'_>,
) -> Result<FiniteGapResponseCurve, FiniteGapResponseCurveError> {
    validate_request(request)?;
    let mut nodes: Vec<FiniteGapResponseNode> = Vec::with_capacity(request.approach_nodes_m.len());
    for (node, &approach_m) in request.approach_nodes_m.iter().enumerate() {
        let contact = solve_finite_gap_half_space(
            &FiniteGapContactRequest {
                grid: request.grid,
                undeformed_gap_m: request.undeformed_gap_m.clone(),
                approach_m,
                reduced_modulus_pa: request.reduced_modulus_pa,
                maximum_active_set_iterations: request.maximum_active_set_iterations,
                complementarity_tolerance_m: request.complementarity_tolerance_m,
                boundary_clearance_cells: request.boundary_clearance_cells,
            },
            cx,
        )
        .map_err(|source| FiniteGapResponseCurveError::Contact { node, source })?;
        if let Some(previous) = nodes.last() {
            // Do not hide even a rounding-scale negative stiffness by calling
            // it tolerance. The retained nodes themselves define the runtime
            // segment, so any decrease would create an actual negative
            // tangent in that segment.
            if contact.normal_force_n < previous.normal_force_n {
                return Err(FiniteGapResponseCurveError::NonMonotoneForce {
                    node,
                    previous_force_n: previous.normal_force_n,
                    force_n: contact.normal_force_n,
                });
            }
        }
        let curve_reversible_energy_j = nodes.last().map_or(0.0, |previous| {
            previous.curve_reversible_energy_j
                + 0.5
                    * (previous.normal_force_n + contact.normal_force_n)
                    * (approach_m - previous.approach_m)
        });
        let energy_residual_j = curve_reversible_energy_j - contact.reversible_energy_j;
        let allowed_j = request.absolute_energy_tolerance_j
            + request.relative_energy_tolerance
                * curve_reversible_energy_j
                    .abs()
                    .max(contact.reversible_energy_j.abs());
        if energy_residual_j.abs() > allowed_j {
            return Err(FiniteGapResponseCurveError::EnergyConsistency {
                node,
                residual_j: energy_residual_j,
                allowed_j,
            });
        }
        nodes.push(FiniteGapResponseNode {
            approach_m,
            contact_identity: contact.identity,
            normal_force_n: contact.normal_force_n,
            peak_pressure_pa: contact.peak_pressure_pa,
            pressure_centroid_m: contact.pressure_centroid_m,
            equivalent_pressure_semiaxes_m: contact.equivalent_pressure_semiaxes_m,
            reference_reversible_energy_j: contact.reversible_energy_j,
            curve_reversible_energy_j,
            energy_residual_j,
        });
    }
    let maximum_absolute_energy_residual_j = nodes
        .iter()
        .map(|node| node.energy_residual_j.abs())
        .fold(0.0_f64, f64::max);
    let identity = response_curve_identity(request, &nodes);
    Ok(FiniteGapResponseCurve {
        identity,
        gap_source_identity: request.gap_source_identity,
        gap_source_authority: request.gap_source_authority,
        response_authority: FiniteGapResponseCurveAuthority::Estimate,
        nodes,
        maximum_absolute_energy_residual_j,
    })
}

impl FiniteGapResponseCurve {
    /// Evaluate without extrapolation. The returned force, tangent, and energy
    /// come from one line segment and therefore satisfy `dU/d(delta) = F`
    /// within that segment up to binary64 rounding.
    pub fn evaluate(
        &self,
        approach_m: f64,
    ) -> Result<FiniteGapInterpolatedResponse, FiniteGapResponseCurveError> {
        let first = self
            .nodes
            .first()
            .ok_or(FiniteGapResponseCurveError::InvalidInput {
                field: "empty response curve",
            })?;
        let last = self.nodes.last().expect("nonempty checked above");
        if !approach_m.is_finite() || approach_m < first.approach_m || approach_m > last.approach_m
        {
            return Err(FiniteGapResponseCurveError::ApproachOutsideTable {
                approach_m,
                minimum_m: first.approach_m,
                maximum_m: last.approach_m,
            });
        }
        let upper = self
            .nodes
            .partition_point(|node| node.approach_m <= approach_m);
        let right_index = upper.clamp(1, self.nodes.len() - 1);
        let left = &self.nodes[right_index - 1];
        let right = &self.nodes[right_index];
        let width = right.approach_m - left.approach_m;
        let fraction = (approach_m - left.approach_m) / width;
        let tangent = (right.normal_force_n - left.normal_force_n) / width;
        let force = tangent.mul_add(approach_m - left.approach_m, left.normal_force_n);
        let segment_approach_m = approach_m - left.approach_m;
        let energy = left.curve_reversible_energy_j
            + left.normal_force_n * segment_approach_m
            + 0.5 * tangent * segment_approach_m * segment_approach_m;
        Ok(FiniteGapInterpolatedResponse {
            approach_m,
            normal_force_n: force,
            normal_tangent_n_per_m: tangent,
            reversible_energy_j: energy,
            peak_pressure_pa: lerp(left.peak_pressure_pa, right.peak_pressure_pa, fraction),
            pressure_centroid_m: [
                lerp(
                    left.pressure_centroid_m[0],
                    right.pressure_centroid_m[0],
                    fraction,
                ),
                lerp(
                    left.pressure_centroid_m[1],
                    right.pressure_centroid_m[1],
                    fraction,
                ),
            ],
            equivalent_pressure_semiaxes_m: [
                lerp(
                    left.equivalent_pressure_semiaxes_m[0],
                    right.equivalent_pressure_semiaxes_m[0],
                    fraction,
                ),
                lerp(
                    left.equivalent_pressure_semiaxes_m[1],
                    right.equivalent_pressure_semiaxes_m[1],
                    fraction,
                ),
            ],
        })
    }
}

impl FiniteGapResponseFamily {
    /// Construct a family whose member curves share exactly the same approach
    /// grid. The common grid is what makes configuration interpolation retain
    /// force/work consistency inside every rectangular parameter cell.
    pub fn try_new(
        coordinate_id: impl Into<String>,
        coordinate_unit: impl Into<String>,
        coordinate_nodes: Vec<f64>,
        curves: Vec<FiniteGapResponseCurve>,
    ) -> Result<Self, FiniteGapResponseCurveError> {
        let coordinate_id = coordinate_id.into();
        let coordinate_unit = coordinate_unit.into();
        if !valid_coordinate_text(&coordinate_id)
            || !valid_coordinate_text(&coordinate_unit)
            || coordinate_nodes.len() < 2
            || coordinate_nodes.len() > MAX_FINITE_GAP_RESPONSE_NODES
            || coordinate_nodes.len() != curves.len()
            || coordinate_nodes.iter().any(|value| !value.is_finite())
            || coordinate_nodes.windows(2).any(|pair| pair[1] <= pair[0])
        {
            return Err(FiniteGapResponseCurveError::InvalidInput {
                field: "response family coordinate",
            });
        }
        let reference_grid = curves
            .first()
            .ok_or(FiniteGapResponseCurveError::InvalidInput {
                field: "response family curves",
            })?
            .nodes
            .iter()
            .map(|node| node.approach_m.to_bits())
            .collect::<Vec<_>>();
        if reference_grid.len() < 2
            || curves.iter().any(|curve| {
                curve.response_authority != FiniteGapResponseCurveAuthority::Estimate
                    || curve.nodes.len() != reference_grid.len()
                    || curve
                        .nodes
                        .iter()
                        .zip(&reference_grid)
                        .any(|(node, expected)| node.approach_m.to_bits() != *expected)
            })
        {
            return Err(FiniteGapResponseCurveError::InvalidInput {
                field: "response family common approach grid",
            });
        }
        let identity =
            response_family_identity(&coordinate_id, &coordinate_unit, &coordinate_nodes, &curves);
        Ok(Self {
            identity,
            coordinate_id,
            coordinate_unit,
            coordinate_nodes,
            curves,
            response_authority: FiniteGapResponseCurveAuthority::Estimate,
        })
    }

    /// Evaluate the bilinear force surface without extrapolating either the
    /// configuration coordinate or approach. At fixed configuration the
    /// interpolated stored energy differentiates to the interpolated force.
    pub fn evaluate(
        &self,
        coordinate: f64,
        approach_m: f64,
    ) -> Result<FiniteGapInterpolatedResponse, FiniteGapResponseCurveError> {
        let first =
            *self
                .coordinate_nodes
                .first()
                .ok_or(FiniteGapResponseCurveError::InvalidInput {
                    field: "empty response family",
                })?;
        let last = *self
            .coordinate_nodes
            .last()
            .expect("nonempty family checked above");
        if !coordinate.is_finite() || coordinate < first || coordinate > last {
            return Err(FiniteGapResponseCurveError::CoordinateOutsideFamily {
                coordinate,
                minimum: first,
                maximum: last,
            });
        }
        let upper = self
            .coordinate_nodes
            .partition_point(|node| *node <= coordinate);
        let right_index = upper.clamp(1, self.coordinate_nodes.len() - 1);
        let left_coordinate = self.coordinate_nodes[right_index - 1];
        let right_coordinate = self.coordinate_nodes[right_index];
        let fraction = (coordinate - left_coordinate) / (right_coordinate - left_coordinate);
        let left = self.curves[right_index - 1].evaluate(approach_m)?;
        let right = self.curves[right_index].evaluate(approach_m)?;
        Ok(FiniteGapInterpolatedResponse {
            approach_m,
            normal_force_n: lerp(left.normal_force_n, right.normal_force_n, fraction),
            normal_tangent_n_per_m: lerp(
                left.normal_tangent_n_per_m,
                right.normal_tangent_n_per_m,
                fraction,
            ),
            reversible_energy_j: lerp(
                left.reversible_energy_j,
                right.reversible_energy_j,
                fraction,
            ),
            peak_pressure_pa: lerp(left.peak_pressure_pa, right.peak_pressure_pa, fraction),
            pressure_centroid_m: [
                lerp(
                    left.pressure_centroid_m[0],
                    right.pressure_centroid_m[0],
                    fraction,
                ),
                lerp(
                    left.pressure_centroid_m[1],
                    right.pressure_centroid_m[1],
                    fraction,
                ),
            ],
            equivalent_pressure_semiaxes_m: [
                lerp(
                    left.equivalent_pressure_semiaxes_m[0],
                    right.equivalent_pressure_semiaxes_m[0],
                    fraction,
                ),
                lerp(
                    left.equivalent_pressure_semiaxes_m[1],
                    right.equivalent_pressure_semiaxes_m[1],
                    fraction,
                ),
            ],
        })
    }
}

fn validate_request(
    request: &FiniteGapResponseCurveRequest,
) -> Result<(), FiniteGapResponseCurveError> {
    let cell_count = request
        .grid
        .cells_x
        .checked_mul(request.grid.cells_y)
        .ok_or(FiniteGapResponseCurveError::InvalidInput {
            field: "grid cell count",
        })?;
    let nodes = request.approach_nodes_m.len();
    if request.grid.cells_x < 3
        || request.grid.cells_y < 3
        || cell_count > MAX_FINITE_GAP_CELLS
        || !(request.grid.cell_width_m.is_finite() && request.grid.cell_width_m > 0.0)
        || !(request.grid.cell_depth_m.is_finite() && request.grid.cell_depth_m > 0.0)
        || request.undeformed_gap_m.len() != cell_count
        || request
            .undeformed_gap_m
            .iter()
            .any(|value| !value.is_finite())
        || !(request.reduced_modulus_pa.is_finite() && request.reduced_modulus_pa > 0.0)
        || nodes < 2
        || nodes > MAX_FINITE_GAP_RESPONSE_NODES
        || request.approach_nodes_m[0].to_bits() != 0.0_f64.to_bits()
        || request
            .approach_nodes_m
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        || request
            .approach_nodes_m
            .windows(2)
            .any(|pair| pair[1] <= pair[0])
        || request.maximum_active_set_iterations == 0
        || !(request.complementarity_tolerance_m.is_finite()
            && request.complementarity_tolerance_m > 0.0)
        || request.boundary_clearance_cells == 0
        || !(request.absolute_energy_tolerance_j.is_finite()
            && request.absolute_energy_tolerance_j >= 0.0)
        || !(request.relative_energy_tolerance.is_finite()
            && request.relative_energy_tolerance >= 0.0)
    {
        return Err(FiniteGapResponseCurveError::InvalidInput {
            field: "response curve request",
        });
    }
    Ok(())
}

fn lerp(left: f64, right: f64, fraction: f64) -> f64 {
    (right - left).mul_add(fraction, left)
}

fn valid_coordinate_text(value: &str) -> bool {
    !value.trim().is_empty() && value.len() <= 128 && value.is_ascii()
}

fn response_family_identity(
    coordinate_id: &str,
    coordinate_unit: &str,
    coordinate_nodes: &[f64],
    curves: &[FiniteGapResponseCurve],
) -> ContentHash {
    let mut bytes = Vec::with_capacity(64 + coordinate_id.len() + coordinate_unit.len());
    bytes.extend_from_slice(&(coordinate_id.len() as u64).to_le_bytes());
    bytes.extend_from_slice(coordinate_id.as_bytes());
    bytes.extend_from_slice(&(coordinate_unit.len() as u64).to_le_bytes());
    bytes.extend_from_slice(coordinate_unit.as_bytes());
    for (&coordinate, curve) in coordinate_nodes.iter().zip(curves) {
        bytes.extend_from_slice(&coordinate.to_bits().to_le_bytes());
        bytes.extend_from_slice(curve.identity.as_bytes());
    }
    hash_domain(FINITE_GAP_RESPONSE_FAMILY_MODEL_ID, &bytes)
}

fn response_curve_identity(
    request: &FiniteGapResponseCurveRequest,
    nodes: &[FiniteGapResponseNode],
) -> ContentHash {
    let mut bytes = Vec::with_capacity(192 + 112 * nodes.len());
    bytes.extend_from_slice(request.gap_source_identity.as_bytes());
    bytes.push(match request.gap_source_authority {
        FiniteGapChartSamplingAuthority::Enclosure => 0,
        FiniteGapChartSamplingAuthority::Estimate => 1,
    });
    bytes.extend_from_slice(&(request.grid.cells_x as u64).to_le_bytes());
    bytes.extend_from_slice(&(request.grid.cells_y as u64).to_le_bytes());
    bytes.extend_from_slice(&request.grid.cell_width_m.to_bits().to_le_bytes());
    bytes.extend_from_slice(&request.grid.cell_depth_m.to_bits().to_le_bytes());
    bytes.extend_from_slice(&request.reduced_modulus_pa.to_bits().to_le_bytes());
    bytes.extend_from_slice(&(request.maximum_active_set_iterations as u64).to_le_bytes());
    bytes.extend_from_slice(&request.complementarity_tolerance_m.to_bits().to_le_bytes());
    bytes.extend_from_slice(&(request.boundary_clearance_cells as u64).to_le_bytes());
    bytes.extend_from_slice(&request.absolute_energy_tolerance_j.to_bits().to_le_bytes());
    bytes.extend_from_slice(&request.relative_energy_tolerance.to_bits().to_le_bytes());
    for value in &request.undeformed_gap_m {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    for node in nodes {
        bytes.extend_from_slice(&node.approach_m.to_bits().to_le_bytes());
        bytes.extend_from_slice(node.contact_identity.as_bytes());
        for value in [
            node.normal_force_n,
            node.peak_pressure_pa,
            node.pressure_centroid_m[0],
            node.pressure_centroid_m[1],
            node.equivalent_pressure_semiaxes_m[0],
            node.equivalent_pressure_semiaxes_m[1],
            node.reference_reversible_energy_j,
            node.curve_reversible_energy_j,
            node.energy_residual_j,
        ] {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    hash_domain(FINITE_GAP_RESPONSE_CURVE_MODEL_ID, &bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 13,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    #[test]
    fn g1_piecewise_force_energy_is_passive_and_derivative_consistent() {
        let grid = FiniteGapGrid {
            cells_x: 17,
            cells_y: 17,
            cell_width_m: 4.0e-5,
            cell_depth_m: 4.0e-5,
        };
        let radius_m = 0.020;
        let mut gap = Vec::with_capacity(grid.cells_x * grid.cells_y);
        for index in 0..grid.cells_x * grid.cells_y {
            let ix = index % grid.cells_x;
            let iy = index / grid.cells_x;
            let x = (ix as f64 + 0.5 - 0.5 * grid.cells_x as f64) * grid.cell_width_m;
            let y = (iy as f64 + 0.5 - 0.5 * grid.cells_y as f64) * grid.cell_depth_m;
            gap.push((x * x + y * y) / (2.0 * radius_m));
        }
        let request = FiniteGapResponseCurveRequest {
            gap_source_identity: hash_domain("test/finite-gap", b"paraboloid"),
            gap_source_authority: FiniteGapChartSamplingAuthority::Enclosure,
            grid,
            undeformed_gap_m: gap,
            reduced_modulus_pa: 80.0e9,
            approach_nodes_m: vec![0.0, 0.75e-6, 1.5e-6, 2.25e-6],
            maximum_active_set_iterations: 256,
            complementarity_tolerance_m: 1.0e-11,
            boundary_clearance_cells: 2,
            absolute_energy_tolerance_j: 1.0e-7,
            relative_energy_tolerance: 0.20,
        };
        let curve = with_cx(|cx| build_finite_gap_response_curve(&request, cx)).unwrap();
        let approach = 1.125e-6;
        let response = curve.evaluate(approach).unwrap();
        assert!(response.normal_force_n > 0.0);
        assert!(response.normal_tangent_n_per_m > 0.0);
        assert!(response.reversible_energy_j > 0.0);
        let epsilon = 1.0e-10;
        let plus = curve.evaluate(approach + epsilon).unwrap();
        let minus = curve.evaluate(approach - epsilon).unwrap();
        let numerical = (plus.reversible_energy_j - minus.reversible_energy_j) / (2.0 * epsilon);
        assert!((numerical / response.normal_force_n - 1.0).abs() < 2.0e-7);
    }

    #[test]
    fn evaluation_refuses_extrapolation() {
        let curve = linear_test_curve(1.0, 1.0);
        assert!(matches!(
            curve.evaluate(1.1),
            Err(FiniteGapResponseCurveError::ApproachOutsideTable { .. })
        ));
    }

    #[test]
    fn g1_response_family_preserves_force_energy_derivative() {
        let family = FiniteGapResponseFamily::try_new(
            "inclination",
            "rad",
            vec![0.05, 0.15],
            vec![linear_test_curve(2.0, 1.0), linear_test_curve(4.0, 1.0)],
        )
        .unwrap();
        let coordinate = 0.10;
        let approach = 0.4;
        let response = family.evaluate(coordinate, approach).unwrap();
        let epsilon = 1.0e-7;
        let plus = family.evaluate(coordinate, approach + epsilon).unwrap();
        let minus = family.evaluate(coordinate, approach - epsilon).unwrap();
        let numerical = (plus.reversible_energy_j - minus.reversible_energy_j) / (2.0 * epsilon);
        assert!((response.normal_force_n - 1.2).abs() < 1.0e-12);
        assert!((response.normal_tangent_n_per_m - 3.0).abs() < 1.0e-12);
        assert!((numerical - response.normal_force_n).abs() < 2.0e-10);
        assert_eq!(
            family.response_authority,
            FiniteGapResponseCurveAuthority::Estimate
        );
    }

    #[test]
    fn response_family_refuses_mismatched_approach_grids_and_extrapolation() {
        let mismatched = FiniteGapResponseFamily::try_new(
            "inclination",
            "rad",
            vec![0.0, 1.0],
            vec![linear_test_curve(1.0, 1.0), linear_test_curve(2.0, 2.0)],
        );
        assert!(matches!(
            mismatched,
            Err(FiniteGapResponseCurveError::InvalidInput {
                field: "response family common approach grid"
            })
        ));
        let family = FiniteGapResponseFamily::try_new(
            "inclination",
            "rad",
            vec![0.0, 1.0],
            vec![linear_test_curve(1.0, 1.0), linear_test_curve(3.0, 1.0)],
        )
        .unwrap();
        assert!(matches!(
            family.evaluate(1.1, 0.5),
            Err(FiniteGapResponseCurveError::CoordinateOutsideFamily { .. })
        ));
    }

    fn linear_test_curve(force_at_maximum: f64, maximum_approach_m: f64) -> FiniteGapResponseCurve {
        let tangent = force_at_maximum / maximum_approach_m;
        FiniteGapResponseCurve {
            identity: ContentHash([0; 32]),
            gap_source_identity: ContentHash([0; 32]),
            gap_source_authority: FiniteGapChartSamplingAuthority::Estimate,
            response_authority: FiniteGapResponseCurveAuthority::Estimate,
            nodes: vec![
                FiniteGapResponseNode {
                    approach_m: 0.0,
                    contact_identity: ContentHash([0; 32]),
                    normal_force_n: 0.0,
                    peak_pressure_pa: 0.0,
                    pressure_centroid_m: [0.0; 2],
                    equivalent_pressure_semiaxes_m: [0.0; 2],
                    reference_reversible_energy_j: 0.0,
                    curve_reversible_energy_j: 0.0,
                    energy_residual_j: 0.0,
                },
                FiniteGapResponseNode {
                    approach_m: maximum_approach_m,
                    contact_identity: ContentHash([0; 32]),
                    normal_force_n: force_at_maximum,
                    peak_pressure_pa: force_at_maximum,
                    pressure_centroid_m: [0.0; 2],
                    equivalent_pressure_semiaxes_m: [maximum_approach_m; 2],
                    reference_reversible_energy_j: 0.5
                        * tangent
                        * maximum_approach_m
                        * maximum_approach_m,
                    curve_reversible_energy_j: 0.5
                        * tangent
                        * maximum_approach_m
                        * maximum_approach_m,
                    energy_residual_j: 0.0,
                },
            ],
            maximum_absolute_energy_residual_j: 0.0,
        }
    }
}
