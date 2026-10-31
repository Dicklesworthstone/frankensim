//! Geometry-based diffuse view factors with retained sampling diagnostics.
//!
//! Rays start uniformly on each named mesh-boundary patch and follow a
//! cosine-weighted outward hemisphere. Every boundary triangle is opaque,
//! including triangles not assigned to a radiation patch. Only the first hit
//! contributes; escaped, unassigned, and back-facing hits are retained losses.
//! Independent Owen-scrambled Sobol replicates give empirical standard errors,
//! not confidence intervals, geometric error bounds, or physical validation.
//!
//! A closed radiosity matrix requires an explicit, bounded projection of the
//! sampled exchange areas onto reciprocity and row closure. It cannot repair
//! sampled losses or introduce an exchange absent in both sampled directions.
//! Zero observed losses do not prove geometric closure. Work is a bounded
//! direct triangle scan, without acceleration or a geometric robustness claim.
//! Open projection retains escaped rays as exchange with an explicit black
//! surrounding reservoir. Hits on unnamed opaque faces are never treated as
//! escape: those faces must be modeled before their heat can be accounted for.

use std::collections::BTreeSet;

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_exec::Cx;
use fs_math::det::{sin_cos, sqrt};
use fs_rand::{Stream, StreamKey, qmc::Sobol};

use super::{
    RadiationSurface, ViewFactorEvidence, ViewFactorMatrix, ViewFactorTolerance, radiation_error,
    radiation_mesh_identity_with_poll,
};
use crate::{ConductionError, ConductionMesh};

const ESTIMATE_DOMAIN: &str = "org.frankensim.fs-conduction.diffuse-geometry-qmc.v1";
const STREAM_DOMAIN: &str = "org.frankensim.fs-conduction.diffuse-patch-stream.v1";

/// Explicit sampling and work limits. There are no hidden retries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiffuseViewFactorConfig {
    /// Study seed. Patch names and replicate indices identify separate streams.
    pub seed: u64,
    /// Positive power of two, per patch per independent scramble.
    pub samples_per_replicate: u32,
    /// Independent scrambles, in 2..=32.
    pub replicates: u32,
    /// Maximum number of opaque boundary triangles, at most one million.
    pub max_triangles: usize,
    /// Budget for the conservative patch * replicate * sample * triangle count.
    pub max_ray_triangle_tests: u64,
}

/// Caller-authorized conservation correction of a loss-free sampled matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViewFactorProjection {
    /// Maximum symmetric diagonal-scaling iterations.
    pub max_iterations: usize,
    /// Positive absolute row-sum tolerance, no larger than 1e-8.
    pub row_sum_tolerance: f64,
    /// Maximum change of any factor from its unmodified sample estimate.
    /// This is an admission budget, not a sampling-error or confidence bound.
    pub max_abs_adjustment: f64,
}

/// Unmodified hit statistics and estimates on the exact supplied mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct DiffuseViewFactorEstimate {
    names: Vec<String>,
    areas: Vec<f64>,
    factors: Vec<Vec<f64>>,
    standard_errors: Vec<Vec<f64>>,
    hit_counts: Vec<Vec<u64>>,
    escaped: Vec<u64>,
    unassigned: Vec<u64>,
    backface: Vec<u64>,
    samples_per_surface: u64,
    ray_triangle_tests: u64,
    mesh_identity: ContentHash,
    identity: ContentHash,
}

impl DiffuseViewFactorEstimate {
    /// Patch names in the caller's order.
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }
    /// Exact selected boundary areas, in square metres.
    #[must_use]
    pub fn areas_m2(&self) -> &[f64] {
        &self.areas
    }
    /// Raw hit fractions. Missing rays are never redistributed.
    #[must_use]
    pub fn factors(&self) -> &[Vec<f64>] {
        &self.factors
    }
    /// Standard deviation of independent replicate means divided by sqrt(R).
    /// Zero means the replicates agreed; it does not prove zero error.
    #[must_use]
    pub fn standard_errors(&self) -> &[Vec<f64>] {
        &self.standard_errors
    }
    /// First, front-facing hits on named patches.
    #[must_use]
    pub fn hit_counts(&self) -> &[Vec<u64>] {
        &self.hit_counts
    }
    /// Rays with no forward boundary intersection.
    #[must_use]
    pub fn escaped_counts(&self) -> &[u64] {
        &self.escaped
    }
    /// First hits on opaque faces outside the named patch partition.
    #[must_use]
    pub fn unassigned_counts(&self) -> &[u64] {
        &self.unassigned
    }
    /// First hits on a named face from its non-radiating side.
    #[must_use]
    pub fn backface_counts(&self) -> &[u64] {
        &self.backface
    }
    /// Total emitted rays from each patch, across all replicates.
    #[must_use]
    pub const fn samples_per_surface(&self) -> u64 {
        self.samples_per_surface
    }
    /// Actual triangle intersection tests, excluding the emitting triangle.
    #[must_use]
    pub const fn ray_triangle_tests(&self) -> u64 {
        self.ray_triangle_tests
    }
    /// Physical mesh identity, including coordinates and tetrahedra.
    #[must_use]
    pub const fn mesh_identity(&self) -> ContentHash {
        self.mesh_identity
    }
    /// Identity of geometry, partition, seed, sampling coordinates and counts.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Symmetrize exchange areas, then symmetrically scale their rows to area.
    ///
    /// This preserves reciprocity and the symmetrized sampled zero pattern.
    /// All sampled losses must be zero. Infeasible support, numerical failure,
    /// exhausted iterations and an excessive correction refuse. No partial
    /// matrix is published and no diagonal exchange is fabricated.
    pub fn project_closed(
        &self,
        cx: &Cx<'_>,
        policy: ViewFactorProjection,
    ) -> Result<ProjectedViewFactors, ConductionError> {
        poll(cx, 0)?;
        if self
            .escaped
            .iter()
            .chain(&self.unassigned)
            .chain(&self.backface)
            .any(|&n| n != 0)
        {
            return Err(refusal(
                "sampled rays leave the named front-facing enclosure",
                "close the physical enclosure and name every receiving patch; never renormalize lost rays",
            ));
        }
        project_reciprocal(cx, &self.areas, &self.factors, self.identity, policy)
    }

    /// Preserve sampled escape to a black surrounding reservoir while imposing
    /// reciprocity on the finite surfaces. The reservoir's effective exchange
    /// area is `sum(A_i F_i,escape)`, not a fictitious meshed surface.
    ///
    /// The same bounded symmetric projection as the closed path acts on the
    /// augmented matrix, including its escape column. Thus the adjustment
    /// budget covers environmental exchange too. Unnamed opaque hits and
    /// back-facing hits refuse; they are not radiation lost to the surroundings.
    pub fn project_open(
        &self,
        cx: &Cx<'_>,
        policy: ViewFactorProjection,
    ) -> Result<ProjectedOpenViewFactors, ConductionError> {
        poll(cx, 0)?;
        if self
            .unassigned
            .iter()
            .chain(&self.backface)
            .any(|&n| n != 0)
        {
            return Err(refusal(
                "sampled radiation hits unmodeled opaque surfaces",
                "name every intercepted radiating face; only rays escaping the mesh may reach the surrounding reservoir",
            ));
        }
        let n = self.areas.len();
        let escape: Vec<f64> = self
            .escaped
            .iter()
            .map(|&count| count as f64 / self.samples_per_surface as f64)
            .collect();
        let bath_area: f64 = self.areas.iter().zip(&escape).map(|(a, f)| a * f).sum();
        if !bath_area.is_finite() || bath_area <= 0.0 {
            return Err(refusal(
                "open projection needs positive representable sampled escape area",
                "use the closed projection for a closed enclosure or refine the escape sampling",
            ));
        }
        let mut areas = self.areas.clone();
        areas.push(bath_area);
        let mut factors = self.factors.clone();
        for (row, &fraction) in factors.iter_mut().zip(&escape) {
            row.push(fraction);
        }
        let mut bath: Vec<f64> = self
            .areas
            .iter()
            .zip(&escape)
            .map(|(a, f)| a * f / bath_area)
            .collect();
        bath.push(0.0);
        factors.push(bath);
        let projected = project_reciprocal(cx, &areas, &factors, self.identity, policy)?;
        let surface_factors: Vec<Vec<f64>> = projected.matrix.factors()[..n]
            .iter()
            .map(|row| row[..n].to_vec())
            .collect();
        // The open solver uses the exact complement of the physical row.
        // Include its rounding difference from the projected bath column in
        // the retained adjustment, rather than concealing it in row tolerance.
        let ambient_factors: Vec<f64> = surface_factors
            .iter()
            .map(|row| 1.0 - row.iter().sum::<f64>())
            .collect();
        let adjustment = ambient_factors
            .iter()
            .zip(&escape)
            .map(|(a, b)| (a - b).abs())
            .fold(projected.max_abs_adjustment, f64::max);
        if ambient_factors.iter().any(|f| !f.is_finite() || *f < 0.0)
            || adjustment > policy.max_abs_adjustment
        {
            return Err(refusal(
                "open conservation correction exceeds its admitted factor budget",
                "increase sampling or explicitly justify a different correction budget",
            ));
        }
        Ok(ProjectedOpenViewFactors {
            surface_factors,
            ambient_factors,
            evidence: ViewFactorEvidence::GeometryQmc {
                estimate: self.identity,
                max_abs_adjustment: adjustment,
            },
            tolerance: ViewFactorTolerance {
                row_sum_abs: policy.row_sum_tolerance,
                reciprocity_rel: 1.0e-12,
            },
            max_abs_adjustment: adjustment,
            iterations: projected.iterations,
        })
    }
}

fn project_reciprocal(
    cx: &Cx<'_>,
    physical_areas: &[f64],
    raw_factors: &[Vec<f64>],
    estimate_identity: ContentHash,
    policy: ViewFactorProjection,
) -> Result<ProjectedViewFactors, ConductionError> {
    if policy.max_iterations == 0
        || policy.max_iterations > 100_000
        || !policy.row_sum_tolerance.is_finite()
        || policy.row_sum_tolerance <= 0.0
        || policy.row_sum_tolerance > 1.0e-8
        || !policy.max_abs_adjustment.is_finite()
        || !(0.0..1.0).contains(&policy.max_abs_adjustment)
    {
        return Err(refusal(
            "invalid reciprocity-projection limits",
            "declare 1..=100000 iterations, row tolerance in (0,1e-8], and factor adjustment in [0,1)",
        ));
    }
    let n = physical_areas.len();
    let scale = physical_areas.iter().copied().fold(0.0_f64, f64::max);
    let areas: Vec<f64> = physical_areas.iter().map(|a| a / scale).collect();
    if areas.iter().any(|a| !a.is_finite() || *a <= 0.0) {
        return Err(refusal(
            "relative patch areas underflow",
            "use a representable enclosure area scale",
        ));
    }
    let mut exchange = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i..n {
            let q = 0.5 * (areas[i] * raw_factors[i][j] + areas[j] * raw_factors[j][i]);
            exchange[i][j] = q;
            exchange[j][i] = q;
        }
    }
    let mut weights = vec![1.0; n];
    for iteration in 0..policy.max_iterations {
        poll(cx, iteration)?;
        let mut factors = vec![vec![0.0; n]; n];
        // Compute each exchange once, so both directions share the same
        // rounded exchange area before division by their own patch area.
        for i in 0..n {
            for j in i..n {
                let q = exchange[i][j] * weights[i] * weights[j];
                factors[i][j] = q / areas[i];
                factors[j][i] = q / areas[j];
            }
        }
        let sums: Vec<f64> = factors.iter().map(|row| row.iter().sum()).collect();
        if sums.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || factors.iter().flatten().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(refusal(
                "sampled exchange support cannot be scaled to positive finite rows",
                "increase sampling or repair the enclosure; no invented exchange is admitted",
            ));
        }
        let residual = sums.iter().map(|s| (s - 1.0).abs()).fold(0.0_f64, f64::max);
        if residual <= 0.5 * policy.row_sum_tolerance {
            // A converged row can lie just above one. Since a view factor
            // must never exceed one, apply ONE common exchange scale,
            // not independent row normalization or pairwise clipping.
            // Retain this correction in the same adjustment budget and
            // recheck closure after rounding. Symmetry/zeros are preserved.
            let normalizer = sums.iter().copied().fold(1.0_f64, f64::max);
            for row in &mut factors {
                for factor in row {
                    *factor /= normalizer;
                }
            }
            let final_residual = factors
                .iter()
                .map(|row| (row.iter().sum::<f64>() - 1.0).abs())
                .fold(0.0_f64, f64::max);
            if final_residual <= policy.row_sum_tolerance {
                let mut adjustment = 0.0_f64;
                for (i, row) in factors.iter().enumerate() {
                    for (j, &factor) in row.iter().enumerate() {
                        adjustment = adjustment.max((factor - raw_factors[i][j]).abs());
                    }
                }
                if adjustment > policy.max_abs_adjustment {
                    return Err(refusal(
                        format!(
                            "conservation correction {adjustment:e} exceeds the declared {} factor budget",
                            policy.max_abs_adjustment
                        ),
                        "increase the sample count or explicitly justify a different correction budget",
                    ));
                }
                let matrix = ViewFactorMatrix::admit(
                    physical_areas.to_vec(),
                    factors,
                    ViewFactorEvidence::GeometryQmc {
                        estimate: estimate_identity,
                        max_abs_adjustment: adjustment,
                    },
                    ViewFactorTolerance {
                        row_sum_abs: policy.row_sum_tolerance,
                        reciprocity_rel: 1.0e-12,
                    },
                )?;
                poll(cx, iteration)?;
                return Ok(ProjectedViewFactors {
                    matrix,
                    max_abs_adjustment: adjustment,
                    iterations: iteration + 1,
                });
            }
        }
        // Simultaneous square-root scaling preserves symmetry and damps
        // the alternating row-normalization oscillation.
        for i in 0..n {
            weights[i] /= sqrt(sums[i]);
            if !weights[i].is_finite() || weights[i] <= 0.0 {
                return Err(refusal(
                    "symmetric conservation scaling is not representable",
                    "increase sampling or repair incompatible patch areas and visibility",
                ));
            }
        }
    }
    Err(refusal(
        format!(
            "reciprocity projection exhausted {} iterations",
            policy.max_iterations
        ),
        "increase the explicit work budget or repair infeasible sampled visibility; no partial matrix is returned",
    ))
}

/// Geometry-derived open-cavity factors. Physical rows retain the external
/// reservoir complement and can be passed directly to `OpenDiffuseEnclosure`.
/// The correction is an admitted numerical model, not a confidence bound.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedOpenViewFactors {
    surface_factors: Vec<Vec<f64>>,
    ambient_factors: Vec<f64>,
    evidence: ViewFactorEvidence,
    tolerance: ViewFactorTolerance,
    max_abs_adjustment: f64,
    iterations: usize,
}

impl ProjectedOpenViewFactors {
    /// Reciprocal physical-surface factors, without the external reservoir row.
    #[must_use]
    pub fn surface_factors(&self) -> &[Vec<f64>] {
        &self.surface_factors
    }
    /// Physical-row complements, directed to the black surroundings.
    #[must_use]
    pub fn ambient_factors(&self) -> &[f64] {
        &self.ambient_factors
    }
    /// Original geometric estimate and actual conservation correction.
    #[must_use]
    pub const fn evidence(&self) -> &ViewFactorEvidence {
        &self.evidence
    }
    /// Admission tolerance used by the shared reciprocal projection.
    #[must_use]
    pub const fn tolerance(&self) -> ViewFactorTolerance {
        self.tolerance
    }
    /// Largest correction, including the surrounding-reservoir exchange.
    #[must_use]
    pub const fn max_abs_adjustment(&self) -> f64 {
        self.max_abs_adjustment
    }
    /// Bounded symmetric-scaling iterations, including the accepted check.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

/// Admitted discrete conservation model, with the actual correction retained.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedViewFactors {
    matrix: ViewFactorMatrix,
    max_abs_adjustment: f64,
    iterations: usize,
}

impl ProjectedViewFactors {
    /// Matrix usable by the existing gray-diffuse radiosity solver.
    #[must_use]
    pub const fn matrix(&self) -> &ViewFactorMatrix {
        &self.matrix
    }
    /// Consume this report and retain its admitted matrix/evidence.
    #[must_use]
    pub fn into_matrix(self) -> ViewFactorMatrix {
        self.matrix
    }
    /// Largest absolute difference from an unmodified sampled factor.
    #[must_use]
    pub const fn max_abs_adjustment(&self) -> f64 {
        self.max_abs_adjustment
    }
    /// Iterations performed, including the accepted row check.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

#[derive(Clone, Copy)]
struct Triangle {
    a: [f64; 3],
    ab: [f64; 3],
    ac: [f64; 3],
    normal: [f64; 3],
}

/// Estimate finite-geometry visibility from 1..=64 named boundary patches.
///
/// Surfaces must be disjoint traces on this exact mesh. Area sampling uses
/// canonical boundary-face order. Names, rather than vector positions, key
/// random streams, so reordering named patches permutes the raw matrix.
/// Geometry is evaluated in normalized mesh coordinates with f64 triangle
/// intersections; mesh discretization and grazing-ray error are not bounded.
pub fn estimate_diffuse_view_factors(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    surfaces: &[RadiationSurface],
    config: DiffuseViewFactorConfig,
) -> Result<DiffuseViewFactorEstimate, ConductionError> {
    poll(cx, 0)?;
    let n = surfaces.len();
    let triangle_count = mesh.boundary().len();
    if !(1..=64).contains(&n)
        || !config.samples_per_replicate.is_power_of_two()
        || !(2..=32).contains(&config.replicates)
        || config.max_triangles == 0
        || config.max_triangles > 1_000_000
        || triangle_count > config.max_triangles
    {
        return Err(refusal(
            "invalid geometry sampling shape or triangle budget",
            "use 1..=64 patches, a power-of-two sample count, 2..=32 replicates, and a sufficient triangle cap up to one million",
        ));
    }
    let samples_per_surface =
        u64::from(config.samples_per_replicate) * u64::from(config.replicates);
    let work = u64::try_from(triangle_count)
        .ok()
        .and_then(|t| t.checked_mul(samples_per_surface))
        .and_then(|t| t.checked_mul(n as u64))
        .ok_or_else(|| {
            refusal(
                "ray-triangle work count overflows",
                "reduce the sampling shape",
            )
        })?;
    if work > config.max_ray_triangle_tests {
        return Err(refusal(
            format!(
                "sampling needs up to {work} triangle tests, exceeding {}",
                config.max_ray_triangle_tests
            ),
            "increase the declared work budget or reduce mesh and sampling size",
        ));
    }
    let mesh_identity = radiation_mesh_identity_with_poll(mesh, |at| poll(cx, at))?;
    let mut names = BTreeSet::new();
    let mut owners = vec![None; triangle_count];
    let mut cumulative = Vec::with_capacity(n);
    for (i, surface) in surfaces.iter().enumerate() {
        poll(cx, i)?;
        if surface.mesh_identity != mesh_identity {
            return Err(refusal(
                "radiation trace belongs to a different physical mesh",
                "bind every patch to the exact mesh being sampled",
            ));
        }
        if !names.insert(surface.name()) {
            return Err(refusal(
                "duplicate radiation patch name",
                "supply unique stable names",
            ));
        }
        let mut area = 0.0;
        let mut ends = Vec::with_capacity(surface.face_slots().len());
        for (local, &slot) in surface.face_slots().iter().enumerate() {
            if local % 256 == 0 {
                poll(cx, local)?;
            }
            if owners[slot].replace(i).is_some() {
                return Err(refusal(
                    "overlapping radiation patches",
                    "partition the named boundary faces without overlap",
                ));
            }
            area += mesh.boundary()[slot].area;
            ends.push(area);
        }
        cumulative.push(ends);
    }
    let triangles = prepare_triangles(cx, mesh)?;
    let sobol = Sobol::new(5);
    let mut counts = vec![vec![0_u64; n]; n];
    let mut factors = vec![vec![0.0; n]; n];
    let mut standard_errors = vec![vec![0.0; n]; n];
    let mut escaped = vec![0_u64; n];
    let mut unassigned = vec![0_u64; n];
    let mut backface = vec![0_u64; n];
    let mut tested = 0_u64;
    let mut identity_hasher = DomainHasher::new(ESTIMATE_DOMAIN);
    identity_hasher.update(mesh_identity.as_bytes());
    identity_hasher.update(&config.seed.to_le_bytes());
    identity_hasher.update(&config.samples_per_replicate.to_le_bytes());
    identity_hasher.update(&config.replicates.to_le_bytes());
    identity_hasher.update(&(n as u64).to_le_bytes());
    for surface in surfaces {
        identity_hasher.update(&(surface.name().len() as u64).to_le_bytes());
        identity_hasher.update(surface.name().as_bytes());
        identity_hasher.update(&(surface.face_slots().len() as u64).to_le_bytes());
        for (local, &slot) in surface.face_slots().iter().enumerate() {
            if local % 256 == 0 {
                poll(cx, local)?;
            }
            identity_hasher.update(&(slot as u64).to_le_bytes());
        }
    }
    for (i, surface) in surfaces.iter().enumerate() {
        let name_hash = hash_domain(STREAM_DOMAIN, surface.name().as_bytes());
        let b = name_hash.as_bytes();
        let name_key = u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
        let mut mean = vec![0.0; n];
        let mut m2 = vec![0.0; n];
        for replicate in 0..config.replicates {
            let words = Stream::at(
                StreamKey {
                    seed: config.seed ^ name_key,
                    kernel: 0x5646_5131,
                    tile: replicate,
                },
                0,
            );
            let scramble = u64::from(words[0]) | (u64::from(words[1]) << 32);
            let mut replicate_counts = vec![0_u64; n];
            let mut losses = [0_u64; 3];
            for sample in 0..config.samples_per_replicate {
                if sample % 64 == 0 {
                    poll(cx, sample as usize)?;
                }
                let mut u = [0.0; 5];
                sobol.point_with_owen_scramble(sample, &mut u, scramble);
                // Center each 32-bit Sobol bin: avoid exactly zero weights
                // and directions lying exactly on hemisphere boundaries.
                for value in &mut u {
                    *value += 1.0 / 8_589_934_592.0;
                }
                let choice = u[0] * surface.area_m2();
                let local_slot = cumulative[i]
                    .partition_point(|&end| end <= choice)
                    .min(surface.face_slots().len() - 1);
                let emitting = surface.face_slots()[local_slot];
                let triangle = &triangles[emitting];
                let radius = sqrt(u[1]);
                let origin = add(
                    triangle.a,
                    add(
                        mul(triangle.ab, radius * (1.0 - u[2])),
                        mul(triangle.ac, radius * u[2]),
                    ),
                );
                let direction = diffuse_direction(triangle.normal, u[3], u[4]);
                let mut closest = f64::INFINITY;
                let mut hit = None;
                for (slot, target) in triangles.iter().enumerate() {
                    if slot == emitting {
                        continue;
                    }
                    if tested % 256 == 0 {
                        poll(cx, sample as usize)?;
                    }
                    tested += 1;
                    if let Some(distance) = intersect(origin, direction, target) {
                        if distance < closest {
                            closest = distance;
                            hit = Some(slot);
                        }
                    }
                }
                match hit {
                    None => losses[0] += 1,
                    Some(slot) => match owners[slot] {
                        None => losses[1] += 1,
                        Some(j) if dot(direction, triangles[slot].normal) < 0.0 => {
                            replicate_counts[j] += 1;
                        }
                        Some(_) => losses[2] += 1,
                    },
                }
            }
            escaped[i] += losses[0];
            unassigned[i] += losses[1];
            backface[i] += losses[2];
            for (j, &count) in replicate_counts.iter().enumerate() {
                counts[i][j] += count;
                let value = count as f64 / f64::from(config.samples_per_replicate);
                let delta = value - mean[j];
                mean[j] += delta / f64::from(replicate + 1);
                m2[j] += delta * (value - mean[j]);
                identity_hasher.update(&count.to_le_bytes());
            }
            for count in losses {
                identity_hasher.update(&count.to_le_bytes());
            }
        }
        for j in 0..n {
            factors[i][j] = counts[i][j] as f64 / samples_per_surface as f64;
            let replicates = f64::from(config.replicates);
            standard_errors[i][j] = sqrt(m2[j].max(0.0) / (replicates * (replicates - 1.0)));
        }
    }
    poll(cx, n)?;
    Ok(DiffuseViewFactorEstimate {
        names: surfaces.iter().map(|s| s.name().to_string()).collect(),
        areas: surfaces.iter().map(RadiationSurface::area_m2).collect(),
        factors,
        standard_errors,
        hit_counts: counts,
        escaped,
        unassigned,
        backface,
        samples_per_surface,
        ray_triangle_tests: tested,
        mesh_identity,
        identity: identity_hasher.finalize(),
    })
}

fn prepare_triangles(cx: &Cx<'_>, mesh: &ConductionMesh) -> Result<Vec<Triangle>, ConductionError> {
    let mut lower = [f64::INFINITY; 3];
    let mut upper = [f64::NEG_INFINITY; 3];
    for (i, point) in mesh.positions().iter().enumerate() {
        if i % 256 == 0 {
            poll(cx, i)?;
        }
        for axis in 0..3 {
            lower[axis] = lower[axis].min(point[axis]);
            upper[axis] = upper[axis].max(point[axis]);
        }
    }
    let scale = (0..3)
        .map(|axis| upper[axis] - lower[axis])
        .fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(refusal(
            "mesh extent cannot be normalized for ray intersections",
            "use finite, representable geometry coordinates",
        ));
    }
    let normalized = |vertex: u32| {
        let p = mesh.positions()[vertex as usize];
        [
            (p[0] - lower[0]) / scale,
            (p[1] - lower[1]) / scale,
            (p[2] - lower[2]) / scale,
        ]
    };
    let mut triangles = Vec::with_capacity(mesh.boundary().len());
    for (slot, face) in mesh.boundary().iter().enumerate() {
        if slot % 256 == 0 {
            poll(cx, slot)?;
        }
        let a = normalized(face.vertices[0]);
        let ab = sub(normalized(face.vertices[1]), a);
        let ac = sub(normalized(face.vertices[2]), a);
        let cross = cross(ab, ac);
        let squared_area = dot(cross, cross);
        if !squared_area.is_finite() || squared_area <= 0.0 {
            return Err(refusal(
                "boundary triangle collapses in normalized ray coordinates",
                "repair geometry with an extreme scale ratio before sampling",
            ));
        }
        triangles.push(Triangle {
            a,
            ab,
            ac,
            normal: face.outward_normal,
        });
    }
    Ok(triangles)
}

fn diffuse_direction(normal: [f64; 3], radial: f64, angular: f64) -> [f64; 3] {
    let helper = if normal[2].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let tangent = cross(helper, normal);
    let tangent = mul(tangent, 1.0 / sqrt(dot(tangent, tangent)));
    let bitangent = cross(normal, tangent);
    let (sin, cos) = sin_cos(core::f64::consts::TAU * angular);
    let radius = sqrt(radial);
    add(
        mul(normal, sqrt(1.0 - radial)),
        add(mul(tangent, radius * cos), mul(bitangent, radius * sin)),
    )
}

/// Two-sided Moller-Trumbore hit on a strictly forward ray. The owning mesh
/// supplies radiating-side orientation; the sorted vertex triple does not.
fn intersect(origin: [f64; 3], direction: [f64; 3], triangle: &Triangle) -> Option<f64> {
    let p = cross(direction, triangle.ac);
    let determinant = dot(triangle.ab, p);
    if determinant == 0.0 {
        return None;
    }
    let t = sub(origin, triangle.a);
    let u = dot(t, p) / determinant;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross(t, triangle.ab);
    let v = dot(direction, q) / determinant;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let distance = dot(triangle.ac, q) / determinant;
    (distance.is_finite() && distance > 0.0).then_some(distance)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn mul(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "radiation-view-factors",
        at,
    })
}
fn refusal(what: impl Into<String>, fix: impl Into<String>) -> ConductionError {
    radiation_error("geometry-view-factors", what, fix)
}
