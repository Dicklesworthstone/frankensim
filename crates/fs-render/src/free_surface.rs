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
//! Subcell droplets, smooth-normal recovery and exact liquid-volume agreement
//! need a different reconstruction and are not claimed by this adapter.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_lbm::d3q19::freesurface3::{FreeSurface3, FreeSurfaceError3};
use fs_lbm::material::free_surface::MaterialFreeSurface3;
use fs_viz::{Grid3, Grid3Error, IsoSurfaceError};

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

/// No mesh or scene mutation is published on refusal.
#[derive(Debug)]
pub enum FreeSurfaceMeshError {
    InvalidOptions,
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
    /// No resolved half-fill crossing remains at the admitted grid resolution.
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
    poll(cx)?;
    if options.origin_m.iter().any(|v| !v.is_finite())
        || !options.cell_size_m.is_finite()
        || options.cell_size_m <= 0.0
        || options.max_cells == 0
        || options.max_triangles == 0
    {
        return Err(FreeSurfaceMeshError::InvalidOptions);
    }
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
    let surface = grid
        .isosurface_with_control(0.0, options.max_triangles, || cx.checkpoint().is_ok())
        .map_err(|error| match error {
            IsoSurfaceError::Cancelled => FreeSurfaceMeshError::Cancelled,
            error => FreeSurfaceMeshError::Isosurface(error),
        })?;
    drop(grid);
    let (vertices, triangles) = surface.into_parts();
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
    for &[a, b, c] in &triangles {
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
    let mesh = TriMesh::new(vertices, triangles);
    poll(cx)?;
    Ok(mesh)
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
