//! Spatial refinement of the original flat-source photon transport owner.
//!
//! Every original layer is divided into equal sublayers. Each original uniform
//! source is shared equally among them; an exterior beam retains its original
//! boundary, direction and full incident energy. All receiving sublayers are reduced
//! to their original physical cell. Angular ordinates, energy groups and material
//! coefficients stay unchanged. The maximum change between successive spatial
//! grids is an observed convergence diagnostic, not a certified error bound,
//! angular/energy convergence result or nuclear-data accuracy statement.

use super::*;
use crate::radiation::{ErrorBudget, scaled_bound};

/// Per-coefficient agreement required between successive spatial grids.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonSlabRefinement {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    /// Maximum number of doublings after the original, level-zero grid.
    pub max_levels: usize,
}

/// Accepted original-cell response and its last observed spatial-grid change.
#[derive(Clone, Debug)]
pub struct RefinedPhotonSlabResponse {
    /// Preparation work is the total over every grid and coefficient reduction.
    pub transport: PhotonSlabResponse,
    pub max_absolute_change: f64,
    pub levels: usize,
    pub subdivisions_per_layer: usize,
}

/// Refine until every original source/group/receiver and escape coefficient
/// satisfies `change <= absolute + relative * max(old, new)`.
///
/// Ballistic kernels need one successful grid comparison. Any positive outgoing
/// redistribution requires two consecutive successful comparisons: a homogeneous
/// scattering layer and its two equal halves can otherwise agree by symmetry
/// while retaining substantial spatial error. A failed comparison resets that
/// agreement count. This rule does not certify the remaining discretization error.
///
/// Exhausted levels, cell/storage limits, work or cancellation return no response;
/// an unconverged coarse grid is never substituted. One work allowance covers
/// all solves, expansion, aggregation and comparisons. Generated-response
/// arithmetic defects follow the original response's energy audit and remain
/// separate from the spatial-grid change.
#[allow(clippy::too_many_arguments)]
pub fn refine_slab_response(
    data_id: String,
    radiation: DecayRadiation,
    group_edges_j: Vec<f64>,
    cells: Vec<PhotonSlabCell>,
    source_ids: Vec<u64>,
    ordinates: Vec<PhotonOrdinate>,
    limits: PhotonSlabLimits,
    refinement: PhotonSlabRefinement,
    mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
) -> Result<RefinedPhotonSlabResponse, RadiationError> {
    checkpoint()?;
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
    refine_slab_response_with_sources(
        data_id,
        radiation,
        group_edges_j,
        cells,
        sources,
        ordinates,
        limits,
        refinement,
        checkpoint,
    )
}

/// Refine internal emission and external irradiation on the same physical slab.
///
/// Only uniform internal emitters are split and averaged. A boundary beam is
/// never divided by the number of spatial cells: its original incident joules
/// enter the same face at the same cosine on every grid. The return map retains
/// all original source and receiving-cell identities. Source IDs need not be
/// cell IDs and distinct sources may illuminate the same boundary or layer.
///
/// The convergence criterion, consecutive-agreement rule, work/storage limits
/// and no-claim boundaries are those of [`refine_slab_response`]. In particular,
/// refinement resolves the scattered flat-source approximation; it does not
/// change the explicitly supplied angular quadrature or scattering phase law.
#[allow(clippy::too_many_arguments)]
pub fn refine_slab_response_with_sources(
    data_id: String,
    radiation: DecayRadiation,
    group_edges_j: Vec<f64>,
    cells: Vec<PhotonSlabCell>,
    sources: Vec<PhotonSlabSource>,
    ordinates: Vec<PhotonOrdinate>,
    limits: PhotonSlabLimits,
    refinement: PhotonSlabRefinement,
    mut checkpoint: impl FnMut() -> Result<(), RadiationError>,
) -> Result<RefinedPhotonSlabResponse, RadiationError> {
    checkpoint()?;
    if !nonnegative(refinement.absolute_tolerance) || !nonnegative(refinement.relative_tolerance) {
        return Err(invalid(
            "finite nonnegative photon refinement tolerances required",
        ));
    }
    if refinement.max_levels == 0 {
        return Err(RadiationError::Budget("photon refinement needs two grids"));
    }
    let groups = group_edges_j.len().checked_sub(1).unwrap_or(0);
    if cells.is_empty() || groups == 0 || sources.is_empty() || ordinates.is_empty() {
        return Err(invalid("nonempty photon refinement inputs required"));
    }
    // Reserve the original input, previous/current reduced maps and reduction
    // scratch alongside the core's conservative allocation bound. Check this
    // before cloning either authored coefficients or an expanded layer grid.
    let retained_values = size(&[3, cells.len(), groups, groups])?
        .checked_add(size(&[8, sources.len(), groups, cells.len()])?)
        .and_then(|n| n.checked_add(size(&[8, sources.len(), groups]).ok()?))
        .and_then(|n| n.checked_add(size(&[4, cells.len()]).ok()?))
        .and_then(|n| n.checked_add(size(&[4, ordinates.len()]).ok()?))
        .and_then(|n| n.checked_add(size(&[2, group_edges_j.len()]).ok()?))
        .ok_or(RadiationError::Budget("photon refinement storage"))?;
    grid_limits(
        cells.len(),
        groups,
        sources.len(),
        ordinates.len(),
        retained_values,
        limits,
    )?;
    let mut work = Work {
        remaining: limits.max_work,
    };
    let mut indices = BTreeMap::new();
    let mut has_redistribution = false;
    for (index, cell) in cells.iter().enumerate() {
        work.charge(1, &mut checkpoint)?;
        if cell.groups.len() != groups
            || !cell.thickness_m.is_finite()
            || cell.thickness_m <= 0.0
            || indices.insert(cell.cell_id, index).is_some()
        {
            return Err(invalid(
                "distinct positive photon layers with exact group coverage required",
            ));
        }
        for (incident, row) in cell.groups.iter().enumerate() {
            work.charge(2 * groups + 1, &mut checkpoint)?;
            validate_collision(row, incident, groups)?;
            has_redistribution |= row.redistribution_per_m.iter().any(|&value| value > 0.0);
        }
    }
    let mut source_indices = Vec::with_capacity(sources.len());
    let mut seen = BTreeSet::new();
    let mut internal_sources = 0_usize;
    for source in &sources {
        work.charge(1, &mut checkpoint)?;
        if !seen.insert(source.source_id) {
            return Err(invalid("duplicate original photon source"));
        }
        let index = match source.geometry {
            PhotonSourceGeometry::UniformLayer { cell_id } => {
                internal_sources += 1;
                Some(
                    *indices
                        .get(&cell_id)
                        .ok_or_else(|| invalid("unknown original photon source layer"))?,
                )
            }
            PhotonSourceGeometry::Collimated { mu, .. } => {
                if !mu.is_finite() || mu <= 0.0 || mu > 1.0 {
                    return Err(invalid("incident beam direction cosine must be in (0, 1]"));
                }
                None
            }
        };
        source_indices.push(index);
    }

    let mut previous: Option<PhotonSlabResponse> = None;
    let mut subdivisions = 1_usize;
    let required_agreements = if has_redistribution { 2 } else { 1 };
    let mut agreement_streak = 0_usize;
    for level in 0..=refinement.max_levels {
        let count = size(&[cells.len(), subdivisions])?;
        let expanded_sources = size(&[internal_sources, subdivisions])?
            .checked_add(sources.len() - internal_sources)
            .ok_or(RadiationError::Budget("photon refinement sources"))?;
        let mut inner = grid_limits(
            count,
            groups,
            expanded_sources,
            ordinates.len(),
            retained_values,
            limits,
        )?;
        let mut expanded = Vec::with_capacity(count);
        let mut emitters = Vec::with_capacity(expanded_sources);
        for cell in &cells {
            let thickness = ratio(cell.thickness_m, subdivisions as f64)?;
            if multiply(thickness, subdivisions as f64)? != cell.thickness_m {
                return Err(RadiationError::Arithmetic);
            }
            for _ in 0..subdivisions {
                work.charge(size(&[groups, groups + 2])? + 1, &mut checkpoint)?;
                expanded.push(PhotonSlabCell {
                    cell_id: expanded.len() as u64,
                    thickness_m: thickness,
                    groups: cell.groups.clone(),
                });
            }
        }
        for (source, &index) in sources.iter().zip(&source_indices) {
            if let Some(index) = index {
                work.charge(subdivisions, &mut checkpoint)?;
                for part in 0..subdivisions {
                    emitters.push(PhotonSlabSource {
                        source_id: emitters.len() as u64,
                        geometry: PhotonSourceGeometry::UniformLayer {
                            cell_id: (index * subdivisions + part) as u64,
                        },
                    });
                }
            } else {
                work.charge(1, &mut checkpoint)?;
                emitters.push(PhotonSlabSource {
                    source_id: emitters.len() as u64,
                    geometry: source.geometry,
                });
            }
        }
        work.charge(group_edges_j.len() + 2 * ordinates.len(), &mut checkpoint)?;
        inner.max_work = work.remaining;
        let solved = prepare_slab_response_with_sources(
            data_id.clone(),
            radiation,
            group_edges_j.clone(),
            expanded,
            emitters,
            ordinates.clone(),
            inner,
            &mut checkpoint,
        )?;
        work.charge(solved.work, &mut checkpoint)?;
        let mut current = aggregate(
            solved,
            &cells,
            &sources,
            subdivisions,
            limits.response,
            &mut work,
            &mut checkpoint,
        )?;
        if let Some(old) = &previous {
            let (converged, change) =
                compare(old, &current, refinement, &mut work, &mut checkpoint)?;
            agreement_streak = if converged { agreement_streak + 1 } else { 0 };
            if agreement_streak >= required_agreements {
                current.work = limits.max_work - work.remaining;
                checkpoint()?;
                return Ok(RefinedPhotonSlabResponse {
                    transport: current,
                    max_absolute_change: change,
                    levels: level,
                    subdivisions_per_layer: subdivisions,
                });
            }
        }
        previous = Some(current);
        if level < refinement.max_levels {
            subdivisions = subdivisions
                .checked_mul(2)
                .ok_or(RadiationError::Budget("photon spatial subdivisions"))?;
        }
    }
    Err(RadiationError::Budget(
        "photon spatial refinement did not converge",
    ))
}

fn grid_limits(
    cells: usize,
    groups: usize,
    sources: usize,
    ordinates: usize,
    retained: usize,
    mut limits: PhotonSlabLimits,
) -> Result<PhotonSlabLimits, RadiationError> {
    limits.max_values = limits
        .max_values
        .checked_sub(retained)
        .ok_or(RadiationError::Budget("photon refinement retained storage"))?;
    admit_slab_dimensions(cells, groups, sources, ordinates, limits)?;
    Ok(limits)
}

#[allow(clippy::too_many_arguments)]
fn aggregate(
    mut solved: PhotonSlabResponse,
    cells: &[PhotonSlabCell],
    sources: &[PhotonSlabSource],
    subdivisions: usize,
    limits: GroupedRadiationResponseLimits,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
) -> Result<PhotonSlabResponse, RadiationError> {
    let groups = solved.response.group_edges_j().len() - 1;
    let mut rows = Vec::with_capacity(sources.len());
    let mut preparation_errors = Vec::with_capacity(sources.len());
    let mut first_emitter = 0_u64;
    for source in sources {
        let source_id = source.source_id;
        let emitters = match source.geometry {
            PhotonSourceGeometry::UniformLayer { .. } => subdivisions,
            PhotonSourceGeometry::Collimated { .. } => 1,
        };
        let weight = 1.0 / emitters as f64;
        let mut group_rows = Vec::with_capacity(groups);
        let mut group_errors = Vec::with_capacity(groups);
        for group in 0..groups {
            let mut received = vec![Sum::default(); cells.len()];
            let mut escaped = Sum::default();
            let mut errors = ErrorBudget::new(f64::MAX);
            for part in 0..emitters {
                let id = first_emitter + part as u64;
                let row = solved
                    .response
                    .source_group_row(id, group)
                    .ok_or(RadiationError::Arithmetic)?;
                work.charge(row.deposits.len() + 2, checkpoint)?;
                let original_error = solved
                    .response
                    .preparation_roundoff(id, group)
                    .ok_or(RadiationError::Arithmetic)?;
                errors.charge(scaled_bound(original_error, weight)?)?;
                let (energy, defect) = product(row.escape_fraction, weight)?;
                errors.charge(defect)?;
                errors.charge(escaped.add(energy)?)?;
                for deposit in &row.deposits {
                    let target = usize::try_from(deposit.receiver_id)
                        .map_err(|_| RadiationError::Arithmetic)?
                        / subdivisions;
                    let sum = received.get_mut(target).ok_or(RadiationError::Arithmetic)?;
                    let (energy, defect) = product(deposit.fraction, weight)?;
                    errors.charge(defect)?;
                    errors.charge(sum.add(energy)?)?;
                }
            }
            work.charge(cells.len() + 1, checkpoint)?;
            let mut escape = escaped.export(&mut errors)?;
            let mut deposits = Vec::with_capacity(cells.len());
            let mut closure = Sum { hi: -1.0, lo: 0.0 };
            errors.charge(closure.add(escape)?)?;
            for (cell, value) in cells.iter().zip(received) {
                let fraction = value.export(&mut errors)?;
                errors.charge(closure.add(fraction)?)?;
                deposits.push(RadiationDeposit {
                    receiver_id: cell.cell_id,
                    fraction,
                });
            }
            let defect = closure.export(&mut errors)?.abs();
            solved.max_response_defect = solved.max_response_defect.max(defect);
            if defect > ROW_CLOSURE_TOLERANCE {
                return Err(RadiationError::Budget("refined photon response closure"));
            }
            escape = bounded_fraction(escape, &mut errors)?;
            for deposit in &mut deposits {
                deposit.fraction = bounded_fraction(deposit.fraction, &mut errors)?;
            }
            group_rows.push(RadiationGroupResponse {
                deposits,
                escape_fraction: escape,
            });
            group_errors.push(errors.used);
        }
        rows.push(GroupedRadiationSourceResponse {
            source_id,
            groups: group_rows,
        });
        preparation_errors.push((source_id, group_errors));
        first_emitter += emitters as u64;
    }
    solved.response = Arc::new(
        GroupedRadiationResponse::new(
            solved.response.data_id().to_owned(),
            solved.response.radiation(),
            solved.response.group_edges_j().to_vec(),
            rows,
            limits,
        )?
        .with_preparation_roundoff(preparation_errors)?,
    );
    checkpoint()?;
    Ok(solved)
}

fn compare(
    old: &PhotonSlabResponse,
    new: &PhotonSlabResponse,
    refinement: PhotonSlabRefinement,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), RadiationError>,
) -> Result<(bool, f64), RadiationError> {
    let mut converged = true;
    let mut maximum = 0.0_f64;
    for source in old.response.source_ids() {
        for group in 0..old.response.group_edges_j().len() - 1 {
            let a = old
                .response
                .source_group_row(source, group)
                .ok_or(RadiationError::Arithmetic)?;
            let b = new
                .response
                .source_group_row(source, group)
                .ok_or(RadiationError::Arithmetic)?;
            work.charge(a.deposits.len() + 1, checkpoint)?;
            if a.deposits.len() != b.deposits.len() {
                return Err(RadiationError::Arithmetic);
            }
            for (left, right) in a.deposits.iter().zip(&b.deposits) {
                if left.receiver_id != right.receiver_id {
                    return Err(RadiationError::Arithmetic);
                }
            }
            let pairs = a
                .deposits
                .iter()
                .zip(&b.deposits)
                .map(|(a, b)| (a.fraction, b.fraction))
                .chain(std::iter::once((a.escape_fraction, b.escape_fraction)));
            for (a, b) in pairs {
                let change = (a - b).abs();
                let allowance = finite(
                    refinement.absolute_tolerance + refinement.relative_tolerance * a.max(b),
                )?;
                maximum = maximum.max(change);
                converged &= change <= allowance;
            }
        }
    }
    Ok((converged, maximum))
}

#[cfg(test)]
mod tests;
