//! Refractive geometry from an immutable native three-dimensional liquid state.
//!
//! The original VOF owner supplies unclamped cell volume fractions. Its cell
//! centers become samples of the explicit piecewise-linear reconstruction
//! `0.5 - fill`; the existing `fs-viz` marching-tetrahedra owner extracts its
//! zero surface, and the existing `TriMesh` builds visibility acceleration.
//! Normals point from liquid into gas. A normal dielectric material can then
//! apply its original spectral Snell/Fresnel, nested-media and Beer laws.
//!
//! This is a resolved visualization model, not a mass-conservative PLIC
//! reconstruction or a signed-distance field. It neither changes the accepted
//! fluid state nor infers optical constants from numerical LBM density. The
//! caller supplies water/glass optics and the physical observation time.
//! Detached liquid with an empty domain collar is admitted; solid contacts,
//! domain-cut surfaces, unlocated carry and unresolved conversion fills refuse.
//! The default keeps the half-fill geometry. An explicit volume-matched entry
//! selects one global isovalue to reproduce the accepted total liquid volume
//! within a caller tolerance, without per-cell or per-component conservation.
//! Neither model recovers subcell shape or smooth normals.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_lbm::d3q19::freesurface3::{FreeSurface3, FreeSurfaceError3};
use fs_lbm::material::free_surface::MaterialFreeSurface3;
use fs_viz::{Grid3, Grid3Error, IsoMesh3, IsoSurfaceError};

use crate::charts::TriMesh;

/// Explicit physical placement and allocation/work limits for one snapshot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceMeshOptions {
    /// Lower corner of lattice cell `(0,0,0)`, in metres.
    pub origin_m: [f64; 3],
    /// Isotropic lattice cell width in metres, independent of optical index.
    pub cell_size_m: f64,
    /// Largest admitted dense cell-center sample count.
    pub max_cells: usize,
    /// Largest admitted output triangle count and associated BVH workload.
    pub max_triangles: usize,
}

/// Numerical volume agreement and total field-pass budget for one snapshot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceVolumeOptions {
    /// Absolute tolerance in cubic metres; nonnegative.
    pub absolute_tolerance_m3: f64,
    /// Relative tolerance in [0,1). At least one tolerance must be positive.
    pub relative_tolerance: f64,
    /// Maximum allocation-free volume queries, each scanning the admitted grid.
    pub max_volume_evaluations: usize,
}

/// Observed agreement of the final physical mesh with the original liquid fill.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceVolumeReport {
    pub target_volume_m3: f64,
    pub enclosed_volume_m3: f64,
    /// Signed enclosed-minus-target residual, evaluated before SI rescaling.
    pub residual_m3: f64,
    /// One global fill threshold in (0,1), rather than the default half fill.
    pub fill_isovalue: f64,
    pub volume_evaluations: usize,
}

/// A native mesh and its numerical total-volume readback.
#[derive(Debug)]
pub struct VolumeMatchedFreeSurfaceMesh {
    pub mesh: TriMesh,
    pub report: FreeSurfaceVolumeReport,
}

/// No mesh or scene mutation is published on refusal.
#[derive(Debug)]
pub enum FreeSurfaceMeshError {
    InvalidOptions,
    InvalidVolumeOptions,
    UnrepresentableVolume,
    /// The requested tolerance was not attained within the pass budget or
    /// the representable isovalue interval. Flat field plateaus can prevent
    /// any single isovalue from representing a requested total volume.
    VolumeMatchExhausted {
        volume_evaluations: usize,
        target_volume_m3: f64,
        closest_volume_m3: f64,
    },
    /// The final closed world-coordinate mesh did not retain the requested
    /// agreement after coordinate conversion and geometric extraction.
    VolumeMismatch {
        target_volume_m3: f64,
        enclosed_volume_m3: f64,
    },
    CellBudgetExceeded {
        required: usize,
        maximum: usize,
    },
    UnlocatedMass,
    UnsupportedSolid {
        coordinates: [usize; 3],
    },
    UnresolvedFill {
        coordinates: [usize; 3],
        fraction: f64,
    },
    WetDomainBoundary {
        coordinates: [usize; 3],
    },
    /// No resolved crossing remains at the admitted grid resolution.
    UnresolvedSurface,
    /// A degenerate facet or an edge lacking two opposite uses prevents a
    /// reliable liquid/gas boundary; no topology repair is inferred.
    InvalidSurface,
    RendererCapacity,
    AllocationRefused,
    Native(FreeSurfaceError3),
    Grid(Grid3Error),
    Isosurface(IsoSurfaceError),
    Cancelled,
}

impl core::fmt::Display for FreeSurfaceMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "free-surface rendering: {self:?}")
    }
}
impl core::error::Error for FreeSurfaceMeshError {}

fn poll(cx: &Cx<'_>) -> Result<(), FreeSurfaceMeshError> {
    cx.checkpoint().map_err(|_| FreeSurfaceMeshError::Cancelled)
}

/// Build the native renderer's outward liquid boundary at one accepted state.
///
/// Sampling is deterministic x-fastest and consumes the actual unclamped
/// interface mass/density ratio. Its dense storage is moved into the shared
/// mesher, then released before BVH construction. No whole fluid clone, new
/// solver, smoothing pass or per-ray reconstruction is introduced.
///
/// Cancellation is polled per sample/cube/facet and before publication. The
/// original BVH build is one synchronous section bounded by `max_triangles`.
/// Edge closure is checked, but global vertex-link manifoldness and physical
/// accuracy below the cell width are not certified. Use the returned mesh as
/// an ordinary `Shape::Mesh` with independently sourced dielectric optics.
pub fn free_surface_mesh(
    state: &FreeSurface3,
    options: FreeSurfaceMeshOptions,
    cx: &Cx<'_>,
) -> Result<TriMesh, FreeSurfaceMeshError> {
    let (grid, _) = sample_grid(state, options, cx)?;
    let surface = extract_surface(&grid, 0.0, options.max_triangles, cx)?;
    drop(grid);
    let (vertices, triangles) = surface.into_parts();
    validate_surface(&vertices, &triangles, cx)?;
    let mesh = TriMesh::new(vertices, triangles);
    poll(cx)?;
    Ok(mesh)
}

fn validate_options(options: FreeSurfaceMeshOptions) -> Result<(), FreeSurfaceMeshError> {
    if options.origin_m.iter().any(|v| !v.is_finite())
        || !options.cell_size_m.is_finite()
        || options.cell_size_m <= 0.0
        || options.max_cells == 0
        || options.max_triangles == 0
    {
        return Err(FreeSurfaceMeshError::InvalidOptions);
    }
    Ok(())
}

fn sample_grid(
    state: &FreeSurface3,
    options: FreeSurfaceMeshOptions,
    cx: &Cx<'_>,
) -> Result<(Grid3, f64), FreeSurfaceMeshError> {
    poll(cx)?;
    validate_options(options)?;
    let dimensions = state.domain_cells();
    let count = dimensions
        .into_iter()
        .try_fold(1usize, usize::checked_mul)
        .ok_or(FreeSurfaceMeshError::RendererCapacity)?;
    if count > options.max_cells {
        return Err(FreeSurfaceMeshError::CellBudgetExceeded {
            required: count,
            maximum: options.max_cells,
        });
    }
    if state.carry() != 0.0 {
        return Err(FreeSurfaceMeshError::UnlocatedMass);
    }
    let (mut fill_sum, mut correction) = (0.0, 0.0);
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| FreeSurfaceMeshError::AllocationRefused)?;
    for z in 0..dimensions[2] {
        for y in 0..dimensions[1] {
            for x in 0..dimensions[0] {
                poll(cx)?;
                let coordinates = [x, y, z];
                let fill = state
                    .volume_fraction_at(coordinates)
                    .map_err(FreeSurfaceMeshError::Native)?
                    .ok_or(FreeSurfaceMeshError::UnsupportedSolid { coordinates })?;
                if !fill.is_finite() || !(0.0..=1.0).contains(&fill) {
                    return Err(FreeSurfaceMeshError::UnresolvedFill {
                        coordinates,
                        fraction: fill,
                    });
                }
                if fill > 0.0
                    && (0..3).any(|axis| {
                        coordinates[axis] == 0 || coordinates[axis] + 1 == dimensions[axis]
                    })
                {
                    return Err(FreeSurfaceMeshError::WetDomainBoundary { coordinates });
                }
                compensated_add(&mut fill_sum, &mut correction, fill);
                values.push(0.5 - fill);
            }
        }
    }
    let lower = options.origin_m.map(|x| x + 0.5 * options.cell_size_m);
    let upper = core::array::from_fn(|axis| {
        options.origin_m[axis] + (dimensions[axis] as f64 - 0.5) * options.cell_size_m
    });
    let grid = Grid3::from_values(dimensions, lower, upper, options.max_cells, values)
        .map_err(FreeSurfaceMeshError::Grid)?;
    Ok((grid, fill_sum + correction))
}

fn isosurface_error(error: IsoSurfaceError) -> FreeSurfaceMeshError {
    match error {
        IsoSurfaceError::Cancelled => FreeSurfaceMeshError::Cancelled,
        error => FreeSurfaceMeshError::Isosurface(error),
    }
}

fn extract_surface(
    grid: &Grid3,
    iso: f64,
    max_triangles: usize,
    cx: &Cx<'_>,
) -> Result<IsoMesh3, FreeSurfaceMeshError> {
    grid.isosurface_with_control(iso, max_triangles, || cx.checkpoint().is_ok())
        .map_err(isosurface_error)
}

fn validate_surface(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    cx: &Cx<'_>,
) -> Result<(), FreeSurfaceMeshError> {
    if triangles.is_empty() {
        return Err(FreeSurfaceMeshError::UnresolvedSurface);
    }
    // A binary BVH has fewer than twice as many nodes as triangles.
    triangles
        .len()
        .checked_mul(2)
        .and_then(|count| u32::try_from(count).ok())
        .ok_or(FreeSurfaceMeshError::RendererCapacity)?;
    let mut edges = BTreeMap::<(u32, u32), (u8, i8)>::new();
    for &[a, b, c] in triangles {
        poll(cx)?;
        let ab: [f64; 3] =
            core::array::from_fn(|i| vertices[b as usize][i] - vertices[a as usize][i]);
        let ac: [f64; 3] =
            core::array::from_fn(|i| vertices[c as usize][i] - vertices[a as usize][i]);
        let n = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        // Retain the original triangle backend's actual normal domain.
        let length = n.iter().map(|x| x * x).sum::<f64>().sqrt();
        if !length.is_finite() || length <= 1e-12 {
            return Err(FreeSurfaceMeshError::InvalidSurface);
        }
        for (u, v) in [(a, b), (b, c), (c, a)] {
            let edge = edges.entry((u.min(v), u.max(v))).or_default();
            if edge.0 == 2 {
                return Err(FreeSurfaceMeshError::InvalidSurface);
            }
            edge.0 += 1;
            edge.1 += if u < v { 1 } else { -1 };
        }
    }
    if edges
        .values()
        .any(|&(count, orientation)| count != 2 || orientation != 0)
    {
        return Err(FreeSurfaceMeshError::InvalidSurface);
    }
    drop(edges);
    poll(cx)?;
    Ok(())
}

fn compensated_add(sum: &mut f64, correction: &mut f64, value: f64) {
    let next = *sum + value;
    *correction += if sum.abs() >= value.abs() {
        (*sum - next) + value
    } else {
        (value - next) + *sum
    };
    *sum = next;
}

/// Use the material flow's original SI cell width, with an explicit placement
/// and the same sampling/output budgets. No optical data are inferred from
/// the mechanical liquid material.
pub fn material_free_surface_mesh(
    state: &MaterialFreeSurface3<'_>,
    origin_m: [f64; 3],
    max_cells: usize,
    max_triangles: usize,
    cx: &Cx<'_>,
) -> Result<TriMesh, FreeSurfaceMeshError> {
    free_surface_mesh(
        state.lattice(),
        FreeSurfaceMeshOptions {
            origin_m,
            cell_size_m: state.scaling().cell_size_m(),
            max_cells,
            max_triangles,
        },
        cx,
    )
}

/// Match one native liquid snapshot's total volume before building its optical mesh.
///
/// The target is the compensated sum of actual volume fractions times the SI
/// cell volume. Bisection uses the original mesher's analytic, allocation-free
/// tetrahedron volume query; only the accepted isovalue constructs a mesh and
/// BVH. One dense sample buffer is retained across the bounded passes.
///
/// Clipping and matching use cell coordinates to avoid large world-origin
/// cancellation. The final transformed mesh must independently meet
/// `max(absolute_tolerance_m3, relative_tolerance * target_volume_m3)`.
/// Coordinate roundoff, an unattainable plateau volume or exhausted work
/// refuses; tolerances are never silently enlarged.
///
/// One global level can change shape and redistribute enclosed volume among
/// disconnected drops. This is not PLIC, per-cell/per-component conservation,
/// a subcell shape recovery, or a bound on physical reconstruction error.
/// Existing solid-contact, empty-collar and topology limitations still apply.
pub fn volume_matched_free_surface_mesh(
    state: &FreeSurface3,
    options: FreeSurfaceMeshOptions,
    volume_options: FreeSurfaceVolumeOptions,
    cx: &Cx<'_>,
) -> Result<VolumeMatchedFreeSurfaceMesh, FreeSurfaceMeshError> {
    poll(cx)?;
    validate_options(options)?;
    if !volume_options.absolute_tolerance_m3.is_finite()
        || volume_options.absolute_tolerance_m3 < 0.0
        || !volume_options.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&volume_options.relative_tolerance)
        || (volume_options.absolute_tolerance_m3 == 0.0 && volume_options.relative_tolerance == 0.0)
        || volume_options.max_volume_evaluations == 0
    {
        return Err(FreeSurfaceMeshError::InvalidVolumeOptions);
    }
    let cell_volume_m3 = options.cell_size_m * options.cell_size_m * options.cell_size_m;
    if !cell_volume_m3.is_finite() || cell_volume_m3 <= 0.0 {
        return Err(FreeSurfaceMeshError::UnrepresentableVolume);
    }
    let (grid, target_cells) = sample_grid(
        state,
        FreeSurfaceMeshOptions {
            origin_m: [0.0; 3],
            cell_size_m: 1.0,
            ..options
        },
        cx,
    )?;
    if target_cells <= 0.0 {
        return Err(FreeSurfaceMeshError::UnresolvedSurface);
    }
    let target_volume_m3 = target_cells * cell_volume_m3;
    if !target_volume_m3.is_finite() || target_volume_m3 <= 0.0 {
        return Err(FreeSurfaceMeshError::UnrepresentableVolume);
    }
    let tolerance_m3 = volume_options
        .absolute_tolerance_m3
        .max(volume_options.relative_tolerance * target_volume_m3);
    if !tolerance_m3.is_finite() || tolerance_m3 <= 0.0 {
        return Err(FreeSurfaceMeshError::UnrepresentableVolume);
    }
    // The empty collar is at +0.5; only interior levels can be published.
    let (mut lower, mut upper, mut iso) = (-0.5, 0.5, 0.0);
    let (mut closest_cells, mut closest_error) = (0.0, f64::INFINITY);
    let mut volume_evaluations = 0;
    let mut accepted = false;
    while volume_evaluations < volume_options.max_volume_evaluations {
        let volume_cells = grid
            .sublevel_volume_with_control(iso, || cx.checkpoint().is_ok())
            .map_err(isosurface_error)?;
        volume_evaluations += 1;
        let difference = volume_cells - target_cells;
        if difference.abs() < closest_error {
            closest_error = difference.abs();
            closest_cells = volume_cells;
        }
        if difference.abs() * cell_volume_m3 <= tolerance_m3 {
            accepted = true;
            break;
        }
        if difference < 0.0 {
            lower = iso;
        } else {
            upper = iso;
        }
        let next = 0.5 * (lower + upper);
        if next <= lower || next >= upper {
            break;
        }
        iso = next;
    }
    if !accepted {
        return Err(FreeSurfaceMeshError::VolumeMatchExhausted {
            volume_evaluations,
            target_volume_m3,
            closest_volume_m3: closest_cells * cell_volume_m3,
        });
    }
    let surface = extract_surface(&grid, iso, options.max_triangles, cx)?;
    drop(grid);
    let (mut vertices, triangles) = surface.into_parts();
    for vertex in &mut vertices {
        poll(cx)?;
        for axis in 0..3 {
            vertex[axis] = options
                .cell_size_m
                .mul_add(vertex[axis], options.origin_m[axis]);
            if !vertex[axis].is_finite() {
                return Err(FreeSurfaceMeshError::InvalidSurface);
            }
        }
    }
    validate_surface(&vertices, &triangles, cx)?;
    let enclosed_cells = enclosed_volume_in_cells(&vertices, &triangles, options.cell_size_m, cx)?;
    let enclosed_volume_m3 = enclosed_cells * cell_volume_m3;
    let residual_m3 = (enclosed_cells - target_cells) * cell_volume_m3;
    if !enclosed_volume_m3.is_finite() || !residual_m3.is_finite() {
        return Err(FreeSurfaceMeshError::UnrepresentableVolume);
    }
    if residual_m3.abs() > tolerance_m3 {
        return Err(FreeSurfaceMeshError::VolumeMismatch {
            target_volume_m3,
            enclosed_volume_m3,
        });
    }
    poll(cx)?;
    let mesh = TriMesh::new(vertices, triangles);
    poll(cx)?;
    Ok(VolumeMatchedFreeSurfaceMesh {
        mesh,
        report: FreeSurfaceVolumeReport {
            target_volume_m3,
            enclosed_volume_m3,
            residual_m3,
            fill_isovalue: 0.5 - iso,
            volume_evaluations,
        },
    })
}

fn enclosed_volume_in_cells(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    cell_size_m: f64,
    cx: &Cx<'_>,
) -> Result<f64, FreeSurfaceMeshError> {
    let reference = vertices[0];
    let (mut sum, mut correction) = (0.0, 0.0);
    for triangle in triangles {
        poll(cx)?;
        let [a, b, c] = triangle.map(|index| {
            core::array::from_fn::<_, 3, _>(|axis| {
                (vertices[index as usize][axis] - reference[axis]) / cell_size_m
            })
        });
        let volume = (a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
        if !volume.is_finite() {
            return Err(FreeSurfaceMeshError::UnrepresentableVolume);
        }
        compensated_add(&mut sum, &mut correction, volume);
    }
    let volume = sum + correction;
    if !volume.is_finite() || volume <= 0.0 {
        return Err(FreeSurfaceMeshError::UnrepresentableVolume);
    }
    Ok(volume)
}

/// Volume-matched geometry using the material trajectory's original SI cell width.
pub fn volume_matched_material_free_surface_mesh(
    state: &MaterialFreeSurface3<'_>,
    origin_m: [f64; 3],
    max_cells: usize,
    max_triangles: usize,
    volume_options: FreeSurfaceVolumeOptions,
    cx: &Cx<'_>,
) -> Result<VolumeMatchedFreeSurfaceMesh, FreeSurfaceMeshError> {
    volume_matched_free_surface_mesh(
        state.lattice(),
        FreeSurfaceMeshOptions {
            origin_m,
            cell_size_m: state.scaling().cell_size_m(),
            max_cells,
            max_triangles,
        },
        volume_options,
        cx,
    )
}
