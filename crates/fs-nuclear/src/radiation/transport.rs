//! Planar multigroup photon transport with vacuum exterior boundaries.
//!
//! Exact characteristic integrals give the first collision and escape fractions
//! for an isotropic, spatially uniform source in each homogeneous layer or a
//! collimated or diffuse field incident on either exterior boundary. Angular
//! integration uses an explicit positive quadrature. Scattered energy is again
//! isotropic and uniform in its receiving layer (the flat-source collision-
//! probability approximation). Repeated same-group scattering is solved with
//! the original prepared LU owner, once factored per group. Lower-energy groups
//! are then solved in descending order; no fixed scattering-order cutoff or
//! stochastic estimator is introduced.
//!
//! The collision-probability formulation follows SCALE/CENTRM, section 7.4.3.3:
//! <https://scale-manual.ornl.gov/6.3.3/CENTRM.html>. Here the transported quantity
//! is energy. Heat and outgoing-energy coefficients must close the total
//! collision coefficient. Attenuation alone is not a local heating law; see
//! <https://pml.nist.gov/PhysRefData/Xcom/Text/intro.html>.
//!
//! This is a steady, infinite-planar, isotropic-scattering model. Spatial and
//! angular refinement are separate physical-discretization obligations. The
//! arithmetic diagnostics below are not bounds on either error. There is no
//! photon flight-time storage, anisotropic phase function, charged-particle
//! transport, upscatter, dose calculation or inferred nuclear/material data.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use fs_la::prepared_lu::LuWorkspace;
use fs_math::det;

use super::grouped::{
    GroupedRadiationResponse, GroupedRadiationResponseLimits, GroupedRadiationSourceResponse,
    RadiationGroupResponse,
};
use super::{
    ErrorBudget, ROW_CLOSURE_TOLERANCE, RadiationDeposit, RadiationError, Sum, finite, nonnegative,
    product, scaled_bound,
};
use crate::endf::decay::spectra::DecayRadiation;

/// Explicit per-mass material data and original-density conversion.
pub mod material;

/// Bounded spatial convergence with original physical cells and sources retained.
pub mod refinement;

/// Energy balance of one incident photon group in a homogeneous material.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotonCollision {
    /// Total collision/removal coefficient of the incident energy [1/m].
    pub total_per_m: f64,
    /// Energy thermalized locally per incident energy and path length [1/m].
    /// Local charged-particle deposition must be justified by the supplied model.
    pub heat_per_m: f64,
    /// Outgoing photon ENERGY coefficients [1/m], in the original group order.
    /// These are not photon-number scattering cross sections. Entries above
    /// the incident energy group must be zero; heat plus all entries must equal
    /// total_per_m. A zero total requires every coefficient to be zero.
    pub redistribution_per_m: Vec<f64>,
}

/// One physical homogeneous layer, in left-to-right slab order.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotonSlabCell {
    pub cell_id: u64,
    pub thickness_m: f64,
    pub groups: Vec<PhotonCollision>,
}

/// Exterior face through which irradiation enters the ordered slab.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhotonSlabBoundary {
    Left,
    Right,
}

/// Spatial and angular distribution of one independently supplied photon source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PhotonSourceGeometry {
    /// Isotropic emission distributed uniformly through an existing layer.
    UniformLayer { cell_id: u64 },
    /// Collimated energy entering an exterior face. `mu` is the positive cosine
    /// relative to that face's inward normal, in `(0, 1]`. Each path length is
    /// the layer thickness divided by `mu`; `mu` need not be a quadrature node.
    Collimated {
        boundary: PhotonSlabBoundary,
        mu: f64,
    },
    /// Constant incident radiance over the inward hemisphere. The supplied
    /// ordinates carry incident-energy fractions `mu * weight / sum(mu * weight)`.
    /// This cosine-weighted rule preserves the supplied face-crossing energy;
    /// its angular approximation must be refined independently of the slab grid.
    Diffuse { boundary: PhotonSlabBoundary },
}

/// A physical source identity and its explicit spatial/angular distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonSlabSource {
    /// Unique source identity; it need not name a receiving material cell.
    pub source_id: u64,
    pub geometry: PhotonSourceGeometry,
}

/// A positive direction cosine and its probability weight on one hemisphere.
/// Internal isotropic emission uses both signs, each with half the source energy.
/// Diffuse boundary irradiation uses the inward sign and normalized cosine weights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonOrdinate {
    pub mu: f64,
    pub weight: f64,
}

/// Four-direction Gauss-Legendre rule (two positive cosines, mirrored).
/// More ordinates may be explicitly supplied for angular refinement.
#[must_use]
pub fn gauss_legendre_4() -> Vec<PhotonOrdinate> {
    let r = det::sqrt(6.0 / 5.0);
    let inner_weight = (18.0 + det::sqrt(30.0)) / 36.0;
    vec![
        PhotonOrdinate {
            mu: det::sqrt((3.0 - 2.0 * r) / 7.0),
            weight: inner_weight,
        },
        PhotonOrdinate {
            mu: det::sqrt((3.0 + 2.0 * r) / 7.0),
            weight: 1.0 - inner_weight,
        },
    ]
}

/// Preparation bounds, in addition to the original response admission limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PhotonSlabLimits {
    /// At most 128 cells, including when this caller limit is larger. This also
    /// bounds the indivisible original LU factorization cancellation interval.
    pub max_cells: usize,
    pub max_groups: usize,
    /// Number of positive-hemisphere ordinates, before mirroring.
    pub max_ordinates: usize,
    /// Conservative count of retained/scratch scalar slots before expansion.
    pub max_values: usize,
    /// Characteristic, matrix assembly, factor and substitution work allowance.
    pub max_work: usize,
    pub response: GroupedRadiationResponseLimits,
}

impl Default for PhotonSlabLimits {
    fn default() -> Self {
        Self {
            max_cells: 128,
            max_groups: 64,
            max_ordinates: 128,
            max_values: 8_000_000,
            max_work: 200_000_000,
            response: GroupedRadiationResponseLimits {
                max_sources: 128,
                max_groups: 64,
                max_receivers: 128,
                max_source_groups: 8192,
                max_transfers: 1_048_576,
            },
        }
    }
}

/// Prepared source-to-receiver map and dimensionless arithmetic diagnostics.
/// The map reuses the original sparse radiation owner in every future interval.
#[derive(Clone, Debug)]
pub struct PhotonSlabResponse {
    pub response: Arc<GroupedRadiationResponse>,
    /// Charged preparation work, available to compose repeated refinement calls.
    pub work: usize,
    /// Maximum absolute first-flight collision-plus-escape closure defect.
    pub max_first_flight_defect: f64,
    /// Maximum sum of absolute original linear-equation residuals for one solve.
    pub max_solve_residual: f64,
    /// Maximum absolute deposited-plus-escaped minus injected unit energy.
    pub max_response_defect: f64,
}

struct GroupPlan {
    // Destination-major first-flight collision probabilities.
    collision: Vec<f64>,
    escape: Vec<f64>,
    first_flight_defect: f64,
    // The original equation is retained for residual evaluation after each RHS.
    equation: Vec<f64>,
    factor: Option<LuWorkspace>,
}

/// Prepare radioactive photon deposition from the supplied material transport
/// model. Each named physical source emits uniformly and isotropically inside
/// its layer. Every layer is retained as a receiver, including transparent or
/// purely scattering layers, so thermal coupling can bind all optical geometry.
///
/// No coefficient normalization, negative-flux clipping or uncounted scattering
/// truncation is performed. A computed unit-energy fraction above one by an
/// admitted roundoff amount is bounded to one and the entire adjustment is
/// charged by the original response's per-joule numerical audit. Numerically
/// unresolved, nonclosing or nonpositive solutions refuse. The caller must refine the
/// material/source spatial discretization and the angular/energy grids as needed.
#[allow(clippy::too_many_arguments)]
pub fn prepare_slab_response(
    data_id: String,
    radiation: DecayRadiation,
    group_edges_j: Vec<f64>,
    cells: Vec<PhotonSlabCell>,
    source_ids: Vec<u64>,
    ordinates: Vec<PhotonOrdinate>,
    limits: PhotonSlabLimits,
    mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
) -> Result<PhotonSlabResponse, RadiationError> {
    checkpoint()?;
    // Bound expansion before allocating the explicit source descriptions.
    admit_slab_dimensions(
        cells.len(),
        group_edges_j.len().checked_sub(1).unwrap_or(0),
        source_ids.len(),
        ordinates.len(),
        limits,
    )?;
    let sources = source_ids
        .into_iter()
        .map(|cell_id| PhotonSlabSource {
            source_id: cell_id,
            geometry: PhotonSourceGeometry::UniformLayer { cell_id },
        })
        .collect();
    prepare_slab_response_with_sources(
        data_id,
        radiation,
        group_edges_j,
        cells,
        sources,
        ordinates,
        limits,
        checkpoint,
    )
}

/// Prepare internal emission and external irradiation with one transport solve.
///
/// Applying the returned response takes energy in joules for each source/group.
/// A boundary source's input is the energy that actually crosses the slab face,
/// not intensity per unit area or solid angle: no area or additional cosine
/// factor is applied. Direct first collisions follow Beer–Lambert attenuation
/// along each incident path. Diffuse irradiation integrates these paths with
/// normalized cosine weights on the supplied positive angular quadrature.
/// Their heat deposits locally; their outgoing photon energy enters the original
/// isotropic flat-source scattering operator. That
/// operator is factored once per group and shared by all source kinds.
///
/// Source identities are independent of receiver identities. Escape includes
/// direct transmission and all photons scattered out through either face.
/// There is no surface reflection/refraction or interface Fresnel law. As in
/// [`prepare_slab_response`], scattered sources require spatial and angular
/// convergence; direct ballistic absorption is analytic for the supplied beam.
/// [`refinement::refine_slab_response_with_sources`] retains these explicit
/// source geometries while refining and reducing to the original receivers.
#[allow(clippy::too_many_arguments)]
pub fn prepare_slab_response_with_sources(
    data_id: String,
    radiation: DecayRadiation,
    group_edges_j: Vec<f64>,
    cells: Vec<PhotonSlabCell>,
    sources: Vec<PhotonSlabSource>,
    ordinates: Vec<PhotonOrdinate>,
    limits: PhotonSlabLimits,
    mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
) -> Result<PhotonSlabResponse, RadiationError> {
    checkpoint()?;
    if !matches!(radiation, DecayRadiation::Gamma | DecayRadiation::XRay) {
        return Err(invalid("a photon radiation kind is required"));
    }
    if data_id.trim().is_empty() {
        return Err(invalid("a photon transport data identity is required"));
    }
    let n = cells.len();
    let groups = group_edges_j.len().checked_sub(1).unwrap_or(0);
    admit_slab_dimensions(n, groups, sources.len(), ordinates.len(), limits)?;
    if !group_edges_j.iter().copied().all(nonnegative)
        || !group_edges_j.windows(2).all(|v| v[0] < v[1])
    {
        return Err(invalid("finite ascending photon energy edges required"));
    }
    let mut work = Work {
        remaining: limits.max_work,
    };
    let mut weights = Sum::default();
    for ordinate in &ordinates {
        work.charge(1, &mut checkpoint)?;
        if !ordinate.mu.is_finite()
            || ordinate.mu <= 0.0
            || ordinate.mu > 1.0
            || !ordinate.weight.is_finite()
            || ordinate.weight <= 0.0
            || ordinate.weight > 1.0
        {
            return Err(invalid("positive direction cosines and weights required"));
        }
        weights.add(ordinate.weight)?;
    }
    if (value(weights)? - 1.0).abs() > ROW_CLOSURE_TOLERANCE {
        return Err(invalid("hemisphere quadrature weights must sum to one"));
    }
    let mut indices = BTreeMap::new();
    let mut heat = vec![0.0; n * groups];
    let mut redistribution = vec![0.0; n * groups * groups];
    let mut kernel_defect = vec![0.0; n * groups];
    for (cell_index, cell) in cells.iter().enumerate() {
        work.charge(1, &mut checkpoint)?;
        if !cell.thickness_m.is_finite() || cell.thickness_m <= 0.0 || cell.groups.len() != groups {
            return Err(invalid(
                "positive layer thickness and exact collision group coverage required",
            ));
        }
        if indices.insert(cell.cell_id, cell_index).is_some() {
            return Err(invalid("duplicate physical photon layer"));
        }
        for (group, row) in cell.groups.iter().enumerate() {
            validate_collision(row, group, groups)?;
            work.charge(groups + 1, &mut checkpoint)?;
            if row.total_per_m != 0.0 {
                heat[cell_index * groups + group] = ratio(row.heat_per_m, row.total_per_m)?;
                let mut closure = Sum::default();
                closure.add(-1.0)?;
                closure.add(heat[cell_index * groups + group])?;
                for (outgoing, &coefficient) in row.redistribution_per_m.iter().enumerate() {
                    let fraction = ratio(coefficient, row.total_per_m)?;
                    redistribution[(cell_index * groups + group) * groups + outgoing] = fraction;
                    closure.add(fraction)?;
                }
                kernel_defect[cell_index * groups + group] = value(closure)?.abs();
            }
        }
    }
    let mut seen = BTreeSet::new();
    for source in &sources {
        work.charge(1, &mut checkpoint)?;
        if !seen.insert(source.source_id) {
            return Err(invalid("unique photon source identities required"));
        }
        match source.geometry {
            PhotonSourceGeometry::UniformLayer { cell_id } if !indices.contains_key(&cell_id) => {
                return Err(invalid("uniform sources must name original photon layers"));
            }
            PhotonSourceGeometry::Collimated { mu, .. }
                if !mu.is_finite() || mu <= 0.0 || mu > 1.0 =>
            {
                return Err(invalid("incident beam direction cosine must be in (0, 1]"));
            }
            _ => {}
        }
    }
    let mut plans = Vec::with_capacity(groups);
    let mut first_flight_defect = 0.0_f64;
    for group in 0..groups {
        let (collision, escape, defect) =
            first_flights(&cells, group, &ordinates, &mut work, &mut checkpoint)?;
        first_flight_defect = first_flight_defect.max(defect);
        let mut equation = vec![0.0; n * n];
        for destination in 0..n {
            work.charge(n, &mut checkpoint)?;
            let retained = redistribution[(destination * groups + group) * groups + group];
            for source in 0..n {
                equation[destination * n + source] =
                    -multiply(collision[destination * n + source], retained)?;
            }
            equation[destination * n + destination] += 1.0;
        }
        // LU is already a bounded, deterministic owner. Charge the full cubic
        // operation before entering it rather than duplicating a solver here.
        let has_self_scatter =
            (0..n).any(|cell| redistribution[(cell * groups + group) * groups + group] != 0.0);
        let factor = if has_self_scatter {
            work.charge(size(&[n, n, n])?, &mut checkpoint)?;
            let mut factor =
                LuWorkspace::new(n).map_err(|_| RadiationError::Budget("photon factor storage"))?;
            factor
                .factor(&equation)
                .map_err(|_| RadiationError::Arithmetic)?;
            Some(factor)
        } else {
            // Without same-group redistribution the characteristic solution is
            // complete. Do not factor or solve an identity system.
            None
        };
        checkpoint()?;
        plans.push(GroupPlan {
            collision,
            escape,
            first_flight_defect: defect,
            equation,
            factor,
        });
    }
    let mut rows = Vec::with_capacity(sources.len());
    let mut preparation_errors = Vec::with_capacity(sources.len());
    let mut max_solve_residual = 0.0_f64;
    let mut max_response_defect = 0.0_f64;
    for source in sources {
        let source_id = source.source_id;
        let mut source_groups = Vec::with_capacity(groups);
        let mut source_errors = Vec::with_capacity(groups);
        for incident in 0..groups {
            let mut rhs = vec![Sum::default(); n * groups];
            let mut deposited = vec![Sum::default(); n];
            let mut escaped = Sum::default();
            let mut preparation_error = ErrorBudget::new(f64::MAX);
            let first_flight = match source.geometry {
                PhotonSourceGeometry::UniformLayer { cell_id } => {
                    rhs[incident * n + indices[&cell_id]].add(1.0)?;
                    None
                }
                PhotonSourceGeometry::Collimated { boundary, mu } => Some(beam_first_flight(
                    &cells,
                    incident,
                    boundary,
                    mu,
                    &mut work,
                    &mut checkpoint,
                )?),
                PhotonSourceGeometry::Diffuse { boundary } => Some(diffuse_first_flight(
                    &cells,
                    incident,
                    boundary,
                    &ordinates,
                    &mut preparation_error,
                    &mut work,
                    &mut checkpoint,
                )?),
            };
            if let Some((collided, transmitted, defect)) = first_flight {
                first_flight_defect = first_flight_defect.max(defect);
                preparation_error.charge(defect)?;
                escaped.add(transmitted)?;
                for (cell, &energy) in collided.iter().enumerate() {
                    work.charge(incident + 2, &mut checkpoint)?;
                    preparation_error.charge(scaled_bound(
                        energy,
                        kernel_defect[cell * groups + incident],
                    )?)?;
                    deposited[cell].add(multiply(energy, heat[cell * groups + incident])?)?;
                    // Include same-group redistribution: these are the
                    // initial isotropic emissions for the repeated-scatter
                    // equation, not another copy of the incident field.
                    for outgoing in 0..=incident {
                        let fraction =
                            redistribution[(cell * groups + incident) * groups + outgoing];
                        rhs[outgoing * n + cell].add(multiply(energy, fraction)?)?;
                    }
                }
            }
            let mut b = vec![0.0; n];
            let mut emitted = vec![0.0; n];
            // Reused across all active groups; included in the existing linear
            // scratch reservation rather than allocated for each correction.
            let mut solve_residual = vec![0.0; n];
            let mut correction = vec![0.0; n];
            let mut residual_total = Sum::default();
            for group in (0..=incident).rev() {
                for cell in 0..n {
                    b[cell] = value(rhs[group * n + cell])?;
                }
                if b.iter().all(|&x| x == 0.0) {
                    continue;
                }
                let plan = &mut plans[group];
                // Only actual, nonnegative group residuals consume the original
                // unit-energy row allowance. Reserving equal shares for all
                // lower groups can reject an uncoupled group solely because
                // of its energy-grid index, even after an optimal scalar solve.
                let used = value(residual_total)?;
                if used > ROW_CLOSURE_TOLERANCE {
                    return Err(RadiationError::Budget("photon equation residual"));
                }
                let residual = solve_emission(
                    plan,
                    &b,
                    &mut emitted,
                    &mut solve_residual,
                    &mut correction,
                    ROW_CLOSURE_TOLERANCE - used,
                    &mut work,
                    &mut checkpoint,
                )?;
                if !emitted.iter().copied().all(nonnegative) {
                    return Err(RadiationError::Arithmetic);
                }
                for destination in 0..n {
                    work.charge(n + group + 1, &mut checkpoint)?;
                    let mut collided = Sum::default();
                    for source in 0..n {
                        collided.add(multiply(
                            plan.collision[destination * n + source],
                            emitted[source],
                        )?)?;
                    }
                    let collided = value(collided)?;
                    preparation_error.charge(scaled_bound(
                        collided,
                        kernel_defect[destination * groups + group],
                    )?)?;
                    deposited[destination]
                        .add(multiply(collided, heat[destination * groups + group])?)?;
                    for outgoing in 0..group {
                        let fraction =
                            redistribution[(destination * groups + group) * groups + outgoing];
                        rhs[outgoing * n + destination].add(multiply(collided, fraction)?)?;
                    }
                }
                max_solve_residual = max_solve_residual.max(residual);
                residual_total.add(residual)?;
                preparation_error.charge(residual)?;
                for cell in 0..n {
                    work.charge(1, &mut checkpoint)?;
                    preparation_error
                        .charge(scaled_bound(emitted[cell], plan.first_flight_defect)?)?;
                    escaped.add(multiply(plan.escape[cell], emitted[cell])?)?;
                }
            }
            // This is an absolute unit-energy equation residual, not a
            // condition-number-scaled residual that can conceal lost joules.
            if value(residual_total)? > ROW_CLOSURE_TOLERANCE {
                return Err(RadiationError::Budget("photon equation residual"));
            }
            let mut escape_fraction = value(escaped)?;
            let mut closure = Sum::default();
            closure.add(-1.0)?;
            closure.add(escape_fraction)?;
            let mut deposits = Vec::with_capacity(n);
            for (cell, energy) in cells.iter().zip(deposited) {
                let fraction = value(energy)?;
                closure.add(fraction)?;
                // Zero entries intentionally retain optical receiver bindings.
                deposits.push(RadiationDeposit {
                    receiver_id: cell.cell_id,
                    fraction,
                });
            }
            let defect = value(closure)?.abs();
            max_response_defect = max_response_defect.max(defect);
            if defect > ROW_CLOSURE_TOLERANCE {
                return Err(RadiationError::Budget(
                    "photon deposited plus escaped energy",
                ));
            }
            escape_fraction = bounded_fraction(escape_fraction, &mut preparation_error)?;
            for deposit in &mut deposits {
                deposit.fraction = bounded_fraction(deposit.fraction, &mut preparation_error)?;
            }
            source_errors.push(preparation_error.used);
            source_groups.push(RadiationGroupResponse {
                deposits,
                escape_fraction,
            });
        }
        rows.push(GroupedRadiationSourceResponse {
            source_id,
            groups: source_groups,
        });
        preparation_errors.push((source_id, source_errors));
    }
    checkpoint()?;
    let response = Arc::new(
        GroupedRadiationResponse::new(data_id, radiation, group_edges_j, rows, limits.response)?
            .with_preparation_roundoff(preparation_errors)?,
    );
    checkpoint()?;
    Ok(PhotonSlabResponse {
        response,
        work: limits.max_work - work.remaining,
        max_first_flight_defect: first_flight_defect,
        max_solve_residual,
        max_response_defect,
    })
}

// Reuse the original LU with at most two defect-correction solves. The
// compensated residual includes signed FMA product tails, so rounded matrix
// products cannot themselves obscure or inflate the unit-energy residual.
#[allow(clippy::too_many_arguments)]
fn solve_emission(
    plan: &mut GroupPlan,
    rhs: &[f64],
    emitted: &mut [f64],
    residual: &mut [f64],
    correction: &mut [f64],
    tolerance: f64,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
) -> Result<f64, RadiationError> {
    let n = rhs.len();
    let Some(factor) = &mut plan.factor else {
        work.charge(n, checkpoint)?;
        emitted.copy_from_slice(rhs);
        return Ok(0.0);
    };
    work.charge(size(&[n, n])?, checkpoint)?;
    factor
        .solve_factored_into(rhs, emitted)
        .map_err(|_| RadiationError::Arithmetic)?;
    for corrections in 0..=2 {
        let mut norm = Sum::default();
        for row in 0..n {
            work.charge(2 * n + 1, checkpoint)?;
            let mut equation = Sum {
                hi: rhs[row],
                lo: 0.0,
            };
            for column in 0..n {
                let coefficient = plan.equation[row * n + column];
                let product = finite(coefficient * emitted[column])?;
                let tail = finite(coefficient.mul_add(emitted[column], -product))?;
                equation.add(-product)?;
                equation.add(-tail)?;
            }
            residual[row] = value(equation)?;
            norm.add(residual[row].abs())?;
        }
        let norm = value(norm)?;
        if norm <= tolerance {
            checkpoint()?;
            return Ok(norm);
        }
        if corrections == 2 {
            break;
        }
        work.charge(size(&[n, n])?, checkpoint)?;
        factor
            .solve_factored_into(residual, correction)
            .map_err(|_| RadiationError::Arithmetic)?;
        work.charge(2 * n, checkpoint)?;
        for (population, &delta) in emitted.iter_mut().zip(correction.iter()) {
            let mut updated = Sum {
                hi: *population,
                lo: 0.0,
            };
            updated.add(delta)?;
            *population = value(updated)?;
        }
    }
    Err(RadiationError::Budget("photon equation residual"))
}

/// Shared admission for density-scaled and directly supplied collision data.
pub(super) fn validate_collision(
    row: &PhotonCollision,
    incident: usize,
    groups: usize,
) -> Result<(), RadiationError> {
    validate_coefficients(
        row.total_per_m,
        row.heat_per_m,
        &row.redistribution_per_m,
        incident,
        groups,
    )
}

pub(super) fn validate_coefficients(
    total: f64,
    heat: f64,
    redistribution: &[f64],
    incident: usize,
    groups: usize,
) -> Result<(), RadiationError> {
    if !nonnegative(total)
        || !nonnegative(heat)
        || redistribution.len() != groups
        || !redistribution.iter().copied().all(nonnegative)
        || redistribution.iter().skip(incident + 1).any(|&x| x != 0.0)
    {
        return Err(invalid(
            "finite nonnegative closed photon coefficients without upscatter required",
        ));
    }
    if total == 0.0 {
        if heat != 0.0 || redistribution.iter().any(|&x| x != 0.0) {
            return Err(invalid(
                "zero collision rate requires zero energy coefficients",
            ));
        }
        return Ok(());
    }
    let mut balance = Sum::default();
    balance.add(-1.0)?;
    balance.add(ratio(heat, total)?)?;
    for &outgoing in redistribution {
        balance.add(ratio(outgoing, total)?)?;
    }
    if value(balance)?.abs() > ROW_CLOSURE_TOLERANCE {
        return Err(invalid(
            "photon heat and outgoing energy must close total collisions",
        ));
    }
    Ok(())
}

fn first_flights(
    cells: &[PhotonSlabCell],
    group: usize,
    ordinates: &[PhotonOrdinate],
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
) -> Result<(Vec<f64>, Vec<f64>, f64), RadiationError> {
    let n = cells.len();
    let mut collisions = vec![Sum::default(); n * n];
    let mut escapes = vec![Sum::default(); n];
    for ordinate in ordinates {
        let mut pass = Vec::with_capacity(n);
        let mut collision = Vec::with_capacity(n);
        let mut born_pass = Vec::with_capacity(n);
        let mut born_collision = Vec::with_capacity(n);
        for cell in cells {
            work.charge(1, checkpoint)?;
            let tau =
                finite(multiply(cell.groups[group].total_per_m, cell.thickness_m)? / ordinate.mu)?;
            let transmission = det::exp(-tau);
            // A positive but unrepresentable transmitted tail is refused; a
            // later arbitrarily large emitted energy must not resurrect it.
            if transmission == 0.0 {
                return Err(RadiationError::Arithmetic);
            }
            let interaction = -det::expm1(-tau);
            let uniform = if tau == 0.0 {
                1.0
            } else {
                ratio(interaction, tau)?
            };
            pass.push(transmission);
            collision.push(interaction);
            born_pass.push(uniform);
            born_collision.push(uniform_collision(tau, uniform)?);
        }
        let half_weight = multiply(ordinate.weight, 0.5)?;
        for source in 0..n {
            work.charge(1, checkpoint)?;
            collisions[source * n + source]
                .add(multiply(ordinate.weight, born_collision[source])?)?;
            let start = multiply(half_weight, born_pass[source])?;
            let mut ray = start;
            for destination in source + 1..n {
                work.charge(1, checkpoint)?;
                collisions[destination * n + source].add(multiply(ray, collision[destination])?)?;
                ray = multiply(ray, pass[destination])?;
            }
            escapes[source].add(ray)?;
            ray = start;
            for destination in (0..source).rev() {
                work.charge(1, checkpoint)?;
                collisions[destination * n + source].add(multiply(ray, collision[destination])?)?;
                ray = multiply(ray, pass[destination])?;
            }
            escapes[source].add(ray)?;
        }
    }
    let collision: Vec<_> = collisions
        .into_iter()
        .map(value)
        .collect::<Result<_, _>>()?;
    let escape: Vec<_> = escapes.into_iter().map(value).collect::<Result<_, _>>()?;
    let mut max_defect = 0.0_f64;
    for source in 0..n {
        work.charge(n, checkpoint)?;
        let mut closure = Sum::default();
        closure.add(-1.0)?;
        closure.add(escape[source])?;
        for destination in 0..n {
            closure.add(collision[destination * n + source])?;
        }
        max_defect = max_defect.max(value(closure)?.abs());
    }
    if max_defect > ROW_CLOSURE_TOLERANCE {
        return Err(RadiationError::Budget("photon first-flight closure"));
    }
    Ok((collision, escape, max_defect))
}

// Exact first collisions of a unit-energy collimated beam. Only the subsequent
// scattered emission uses the flat-source approximation or angular quadrature.
fn beam_first_flight(
    cells: &[PhotonSlabCell],
    group: usize,
    boundary: PhotonSlabBoundary,
    mu: f64,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
) -> Result<(Vec<f64>, f64, f64), RadiationError> {
    let mut collided = vec![0.0; cells.len()];
    let mut ray = 1.0;
    let mut closure = Sum::default();
    closure.add(-1.0)?;
    for offset in 0..cells.len() {
        work.charge(1, checkpoint)?;
        let index = match boundary {
            PhotonSlabBoundary::Left => offset,
            PhotonSlabBoundary::Right => cells.len() - 1 - offset,
        };
        let cell = &cells[index];
        let tau = finite(multiply(cell.groups[group].total_per_m, cell.thickness_m)? / mu)?;
        let transmission = det::exp(-tau);
        if transmission == 0.0 {
            return Err(RadiationError::Arithmetic);
        }
        collided[index] = multiply(ray, -det::expm1(-tau))?;
        closure.add(collided[index])?;
        ray = multiply(ray, transmission)?;
    }
    closure.add(ray)?;
    let defect = value(closure)?.abs();
    if defect > ROW_CLOSURE_TOLERANCE {
        return Err(RadiationError::Budget("photon beam first-flight closure"));
    }
    Ok((collided, ray, defect))
}

// Constant radiance gives incident energy dQ proportional to mu*dmu. Normalize
// the finite quadrature's first angular moment, not the transport coefficients:
// a mirrored full-range rule need not integrate that half-range moment exactly.
fn diffuse_first_flight(
    cells: &[PhotonSlabCell],
    group: usize,
    boundary: PhotonSlabBoundary,
    ordinates: &[PhotonOrdinate],
    errors: &mut ErrorBudget,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
) -> Result<(Vec<f64>, f64, f64), RadiationError> {
    let mut normalization = Sum::default();
    for ordinate in ordinates {
        work.charge(1, checkpoint)?;
        normalization.add(multiply(ordinate.mu, ordinate.weight)?)?;
    }
    let normalization = value(normalization)?;
    let mut collided = vec![Sum::default(); cells.len()];
    let mut transmitted = Sum::default();
    for ordinate in ordinates {
        let weight = ratio(multiply(ordinate.mu, ordinate.weight)?, normalization)?;
        let (ray_collisions, escape, defect) =
            beam_first_flight(cells, group, boundary, ordinate.mu, work, checkpoint)?;
        // Directional defects must not cancel when their angular average closes.
        // The aggregate closure is checked below and charged by the caller.
        errors.charge(scaled_bound(weight, defect)?)?;
        work.charge(cells.len() + 1, checkpoint)?;
        for (total, energy) in collided.iter_mut().zip(ray_collisions) {
            total.add(multiply(weight, energy)?)?;
        }
        transmitted.add(multiply(weight, escape)?)?;
    }
    work.charge(cells.len() + 1, checkpoint)?;
    let collided = collided
        .into_iter()
        .map(value)
        .collect::<Result<Vec<_>, _>>()?;
    let transmitted = value(transmitted)?;
    let mut closure = Sum { hi: -1.0, lo: 0.0 };
    closure.add(transmitted)?;
    for &energy in &collided {
        closure.add(energy)?;
    }
    let defect = value(closure)?.abs();
    if defect > ROW_CLOSURE_TOLERANCE {
        return Err(RadiationError::Budget(
            "diffuse photon first-flight closure",
        ));
    }
    Ok((collided, transmitted, defect))
}

// 1 - (1-exp(-tau))/tau. The power series retains collision energy in optically
// thin emitting cells when subtracting the escape fraction would round to zero.
fn uniform_collision(tau: f64, escape: f64) -> Result<f64, RadiationError> {
    if tau == 0.0 {
        return Ok(0.0);
    }
    let value = if tau < 0.125 {
        let mut term = tau * 0.5;
        let mut sum = Sum::default();
        sum.add(term)?;
        for order in 2..=12 {
            term *= -tau / (order as f64 + 1.0);
            sum.add(term)?;
        }
        value(sum)?
    } else {
        1.0 - escape
    };
    if !nonnegative(value) || value == 0.0 {
        return Err(RadiationError::Arithmetic);
    }
    Ok(value)
}

fn ratio(numerator: f64, denominator: f64) -> Result<f64, RadiationError> {
    let value = finite(numerator / denominator)?;
    if numerator > 0.0 && value == 0.0 {
        return Err(RadiationError::Arithmetic);
    }
    Ok(value)
}

fn multiply(a: f64, b: f64) -> Result<f64, RadiationError> {
    Ok(product(a, b)?.0)
}

fn bounded_fraction(value: f64, error: &mut ErrorBudget) -> Result<f64, RadiationError> {
    if !nonnegative(value) {
        return Err(RadiationError::Arithmetic);
    }
    if value <= 1.0 {
        return Ok(value);
    }
    let excess = value - 1.0;
    if excess > ROW_CLOSURE_TOLERANCE {
        return Err(RadiationError::Budget("photon response fraction"));
    }
    error.charge(excess)?;
    Ok(1.0)
}

fn value(sum: Sum) -> Result<f64, RadiationError> {
    finite(sum.hi + sum.lo)
}

fn size(dimensions: &[usize]) -> Result<usize, RadiationError> {
    dimensions
        .iter()
        .try_fold(1_usize, |n, &d| n.checked_mul(d))
        .ok_or(RadiationError::Budget("photon extent overflow"))
}

/// Shared pre-allocation admission for core transport and physical adapters.
pub(crate) fn admit_slab_dimensions(
    cells: usize,
    groups: usize,
    sources: usize,
    ordinates: usize,
    limits: PhotonSlabLimits,
) -> Result<(), RadiationError> {
    if cells == 0 || groups == 0 || sources == 0 || ordinates == 0 {
        return Err(invalid(
            "nonempty layers, groups, sources and ordinates required",
        ));
    }
    if cells > limits.max_cells
        || cells > 128
        || groups > limits.max_groups
        || ordinates > limits.max_ordinates
        || sources > limits.response.max_sources
        || groups > limits.response.max_groups
        || cells > limits.response.max_receivers
    {
        return Err(RadiationError::Budget("photon transport dimensions"));
    }
    let source_groups = size(&[sources, groups])?;
    let transfers = size(&[source_groups, cells])?;
    if source_groups > limits.response.max_source_groups
        || transfers > limits.response.max_transfers
    {
        return Err(RadiationError::Budget("photon response expansion"));
    }
    if storage_values(cells, groups, sources)? > limits.max_values {
        return Err(RadiationError::Budget("photon scalar storage"));
    }
    Ok(())
}

fn storage_values(cells: usize, groups: usize, sources: usize) -> Result<usize, RadiationError> {
    size(&[7, groups, cells, cells])?
        .checked_add(size(&[3, cells, groups, groups])?)
        .and_then(|x| x.checked_add(size(&[32, cells, groups]).ok()?))
        .and_then(|x| x.checked_add(size(&[4, sources, groups, cells]).ok()?))
        // Explicit source geometry stays live while output rows accumulate.
        .and_then(|x| x.checked_add(size(&[4, sources]).ok()?))
        .ok_or(RadiationError::Budget("photon scalar storage"))
}

fn invalid(detail: &'static str) -> RadiationError {
    RadiationError::InvalidInput(detail)
}

struct Work {
    remaining: usize,
}
impl Work {
    fn charge(
        &mut self,
        count: usize,
        checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
    ) -> Result<(), RadiationError> {
        checkpoint()?;
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or(RadiationError::Budget("photon transport work"))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
