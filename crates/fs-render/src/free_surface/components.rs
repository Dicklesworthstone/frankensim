//! Independent volume targets for disconnected native positive-fill supports.

use super::*;
use fs_viz::SublevelComponents3;

/// Additional aggregate admission limits for support-component reconstruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FreeSurfaceComponentOptions {
    /// Largest admitted number of disjoint positive-fill supports.
    pub max_components: usize,
    /// Largest sum of all component crop node counts, including gas halos.
    /// Overlapping crop boxes each consume their actual node count.
    pub max_crop_cells: usize,
}

/// One support component, ordered by its first x-fastest original cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceComponentReport {
    /// First positive-fill cell; an identity only within this snapshot.
    pub seed_cell: [usize; 3],
    pub occupied_cells: usize,
    pub volume: FreeSurfaceVolumeReport,
}

/// One renderer mesh/BVH with independent component and combined volume audits.
#[derive(Debug)]
pub struct ComponentMatchedFreeSurfaceMesh {
    pub mesh: TriMesh,
    pub components: Vec<FreeSurfaceComponentReport>,
    pub target_volume_m3: f64,
    pub enclosed_volume_m3: f64,
    pub residual_m3: f64,
    pub volume_evaluations: usize,
    pub cropped_cells: usize,
}

#[derive(Clone, Copy)]
struct Support {
    seed_cell: [usize; 3],
    lower: [usize; 3],
    upper: [usize; 3],
    occupied_cells: usize,
    fill_sum: f64,
    fill_correction: f64,
    crop_cells: usize,
    tolerance_m3: f64,
}

impl Support {
    fn target_cells(self) -> f64 {
        self.fill_sum + self.fill_correction
    }
}

/// Match each disconnected positive-fill support before building one optical mesh.
///
/// Supports are the connected strict-sublevel components of the raw negative-fill
/// field under the original six-tetrahedra decomposition. Body/face diagonals
/// follow that owner's actual edges; zero-fill cells never join components.
/// Labels are snapshot identities, not tracked drops. A connected support may
/// itself separate at its selected positive fill level; no physical topology
/// or per-cell PLIC reconstruction is inferred.
///
/// Each support gets its actual compensated sum of accepted cell fills and an
/// isovalue selected by the original allocation-free volume solver. Its crop
/// keeps a full one-node gas halo in original cell coordinates and masks every
/// other label, even when crop boxes overlap. Only the accepted level is meshed.
/// Raw positive fills are labeled before half-minus-fill can round them away; if
/// the shared geometric model cannot resolve them, the entire operation refuses.
///
/// The total tolerance is max(absolute, relative * total_target); each
/// component receives that allowance in proportion to its original target.
/// Underflow refuses. Component audits and a final combined closed-mesh audit
/// must all meet their original allowances. Query and triangle limits are
/// aggregate across supports; max_crop_cells bounds summed crop work, while
/// max_cells still bounds the original dense sample. One final BVH is built,
/// and cancellation/refusal never returns a partial mesh.
pub fn component_volume_matched_free_surface_mesh(
    state: &FreeSurface3,
    options: FreeSurfaceMeshOptions,
    volume_options: FreeSurfaceVolumeOptions,
    component_options: FreeSurfaceComponentOptions,
    cx: &Cx<'_>,
) -> Result<ComponentMatchedFreeSurfaceMesh, FreeSurfaceMeshError> {
    poll(cx)?;
    validate_options(options)?;
    validate_volume_options(volume_options)?;
    if component_options.max_components == 0 || component_options.max_crop_cells == 0 {
        return Err(FreeSurfaceMeshError::InvalidComponentOptions);
    }
    let cell_volume = cell_volume_m3(options.cell_size_m)?;
    let (raw_grid, target_cells) = sample_grid_with_offset(
        state,
        FreeSurfaceMeshOptions {
            origin_m: [0.0; 3],
            cell_size_m: 1.0,
            ..options
        },
        0.0,
        cx,
    )?;
    if target_cells <= 0.0 {
        return Err(FreeSurfaceMeshError::UnresolvedSurface);
    }
    let target_volume_m3 = target_cells * cell_volume;
    let tolerance_m3 = volume_options
        .absolute_tolerance_m3
        .max(volume_options.relative_tolerance * target_volume_m3);
    if !target_volume_m3.is_finite()
        || target_volume_m3 <= 0.0
        || !tolerance_m3.is_finite()
        || tolerance_m3 <= 0.0
    {
        return Err(FreeSurfaceMeshError::UnrepresentableVolume);
    }
    let labels = raw_grid
        .sublevel_components_with_control(0.0, component_options.max_components, || {
            cx.checkpoint().is_ok()
        })
        .map_err(isosurface_error)?;
    let (mut supports, cropped_cells) =
        support_bounds(&raw_grid, &labels, component_options.max_crop_cells, cx)?;
    // Admit every proportional allowance before any component volume passes.
    for support in &mut supports {
        poll(cx)?;
        let share = tolerance_m3 * (support.target_cells() / target_cells);
        if !share.is_finite() || share <= 0.0 || share > tolerance_m3 {
            return Err(FreeSurfaceMeshError::UnrepresentableVolume);
        }
        support.tolerance_m3 = share;
    }
    let mut components = Vec::new();
    components
        .try_reserve_exact(supports.len())
        .map_err(|_| FreeSurfaceMeshError::AllocationRefused)?;
    let (mut vertices, mut triangles) = (Vec::new(), Vec::new());
    let mut volume_evaluations = 0usize;
    for (label, support) in supports.into_iter().enumerate() {
        poll(cx)?;
        let component_target = support.target_cells();
        let remaining_triangles = options.max_triangles - triangles.len();
        if remaining_triangles == 0 {
            return Err(FreeSurfaceMeshError::Isosurface(
                IsoSurfaceError::TriangleBudgetExceeded {
                    limit: options.max_triangles,
                },
            ));
        }
        let grid = crop_support(&raw_grid, &labels, label, support, cx)?;
        let matched = match match_grid_volume(
            grid,
            component_target,
            FreeSurfaceMeshOptions {
                max_triangles: remaining_triangles,
                ..options
            },
            FreeSurfaceVolumeOptions {
                absolute_tolerance_m3: support.tolerance_m3,
                relative_tolerance: 0.0,
                max_volume_evaluations: volume_options.max_volume_evaluations - volume_evaluations,
            },
            cx,
        ) {
            Err(FreeSurfaceMeshError::VolumeMatchExhausted {
                volume_evaluations: used,
                target_volume_m3,
                closest_volume_m3,
            }) => {
                return Err(FreeSurfaceMeshError::VolumeMatchExhausted {
                    volume_evaluations: volume_evaluations + used,
                    target_volume_m3,
                    closest_volume_m3,
                });
            }
            result => result?,
        };
        volume_evaluations += matched.report.volume_evaluations;
        let offset =
            u32::try_from(vertices.len()).map_err(|_| FreeSurfaceMeshError::RendererCapacity)?;
        vertices
            .len()
            .checked_add(matched.vertices.len())
            .and_then(|count| u32::try_from(count).ok())
            .ok_or(FreeSurfaceMeshError::RendererCapacity)?;
        triangles
            .len()
            .checked_add(matched.triangles.len())
            .and_then(|count| count.checked_mul(2))
            .and_then(|count| u32::try_from(count).ok())
            .ok_or(FreeSurfaceMeshError::RendererCapacity)?;
        vertices
            .try_reserve(matched.vertices.len())
            .map_err(|_| FreeSurfaceMeshError::AllocationRefused)?;
        triangles
            .try_reserve(matched.triangles.len())
            .map_err(|_| FreeSurfaceMeshError::AllocationRefused)?;
        for vertex in matched.vertices {
            poll(cx)?;
            vertices.push(vertex);
        }
        for triangle in matched.triangles {
            poll(cx)?;
            let [a, b, c] = triangle;
            triangles.push([
                a.checked_add(offset)
                    .ok_or(FreeSurfaceMeshError::RendererCapacity)?,
                b.checked_add(offset)
                    .ok_or(FreeSurfaceMeshError::RendererCapacity)?,
                c.checked_add(offset)
                    .ok_or(FreeSurfaceMeshError::RendererCapacity)?,
            ]);
        }
        components.push(FreeSurfaceComponentReport {
            seed_cell: support.seed_cell,
            occupied_cells: support.occupied_cells,
            volume: matched.report,
        });
    }
    drop(labels);
    drop(raw_grid);
    validate_surface(&vertices, &triangles, cx)?;
    // Audit the actual combined arrays after remapping, under the original
    // total allowance. Summing rounded component reports is not this audit.
    let enclosed_cells = enclosed_volume_in_cells(&vertices, &triangles, options.cell_size_m, cx)?;
    let enclosed_volume_m3 = enclosed_cells * cell_volume;
    let residual_m3 = (enclosed_cells - target_cells) * cell_volume;
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
    Ok(ComponentMatchedFreeSurfaceMesh {
        mesh,
        components,
        target_volume_m3,
        enclosed_volume_m3,
        residual_m3,
        volume_evaluations,
        cropped_cells,
    })
}

fn support_bounds(
    grid: &Grid3,
    labels: &SublevelComponents3,
    max_crop_cells: usize,
    cx: &Cx<'_>,
) -> Result<(Vec<Support>, usize), FreeSurfaceMeshError> {
    let mut supports = Vec::new();
    supports
        .try_reserve_exact(labels.component_count())
        .map_err(|_| FreeSurfaceMeshError::AllocationRefused)?;
    supports.resize(
        labels.component_count(),
        Support {
            seed_cell: [0; 3],
            lower: [usize::MAX; 3],
            upper: [0; 3],
            occupied_cells: 0,
            fill_sum: 0.0,
            fill_correction: 0.0,
            crop_cells: 0,
            tolerance_m3: 0.0,
        },
    );
    let [nx, ny, _] = grid.dimensions();
    for (node, &negative_fill) in grid.values().iter().enumerate() {
        poll(cx)?;
        let Some(label) = labels.component_of_node(node) else {
            continue;
        };
        let coordinates = [node % nx, (node / nx) % ny, node / (nx * ny)];
        let support = &mut supports[label];
        if support.occupied_cells == 0 {
            support.seed_cell = coordinates;
        }
        support.occupied_cells += 1;
        for axis in 0..3 {
            support.lower[axis] = support.lower[axis].min(coordinates[axis]);
            support.upper[axis] = support.upper[axis].max(coordinates[axis]);
        }
        compensated_add(
            &mut support.fill_sum,
            &mut support.fill_correction,
            -negative_fill,
        );
    }
    let mut cropped_cells = 0usize;
    for support in &mut supports {
        poll(cx)?;
        // The original sampler already required a dry outer domain collar.
        for axis in 0..3 {
            support.lower[axis] -= 1;
            support.upper[axis] += 1;
        }
        support.crop_cells = (0..3)
            .try_fold(1usize, |count, axis| {
                count.checked_mul(support.upper[axis] - support.lower[axis] + 1)
            })
            .ok_or(FreeSurfaceMeshError::RendererCapacity)?;
        cropped_cells = cropped_cells
            .checked_add(support.crop_cells)
            .ok_or(FreeSurfaceMeshError::RendererCapacity)?;
        if cropped_cells > max_crop_cells {
            return Err(FreeSurfaceMeshError::CropBudgetExceeded {
                required: cropped_cells,
                maximum: max_crop_cells,
            });
        }
    }
    Ok((supports, cropped_cells))
}

fn crop_support(
    raw_grid: &Grid3,
    labels: &SublevelComponents3,
    label: usize,
    support: Support,
    cx: &Cx<'_>,
) -> Result<Grid3, FreeSurfaceMeshError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(support.crop_cells)
        .map_err(|_| FreeSurfaceMeshError::AllocationRefused)?;
    let [nx, ny, _] = raw_grid.dimensions();
    for z in support.lower[2]..=support.upper[2] {
        for y in support.lower[1]..=support.upper[1] {
            for x in support.lower[0]..=support.upper[0] {
                poll(cx)?;
                let node = (z * ny + y) * nx + x;
                values.push(if labels.component_of_node(node) == Some(label) {
                    0.5 + raw_grid.values()[node]
                } else {
                    0.5
                });
            }
        }
    }
    let dimensions = core::array::from_fn(|axis| support.upper[axis] - support.lower[axis] + 1);
    let lower = support.lower.map(|coordinate| coordinate as f64 + 0.5);
    let upper = support.upper.map(|coordinate| coordinate as f64 + 0.5);
    Grid3::from_values(dimensions, lower, upper, support.crop_cells, values)
        .map_err(FreeSurfaceMeshError::Grid)
}

/// Component-matched geometry with the native material trajectory's SI width.
pub fn component_volume_matched_material_free_surface_mesh(
    state: &MaterialFreeSurface3<'_>,
    origin_m: [f64; 3],
    max_cells: usize,
    max_triangles: usize,
    volume_options: FreeSurfaceVolumeOptions,
    component_options: FreeSurfaceComponentOptions,
    cx: &Cx<'_>,
) -> Result<ComponentMatchedFreeSurfaceMesh, FreeSurfaceMeshError> {
    component_volume_matched_free_surface_mesh(
        state.lattice(),
        FreeSurfaceMeshOptions {
            origin_m,
            cell_size_m: state.scaling().cell_size_m(),
            max_cells,
            max_triangles,
        },
        volume_options,
        component_options,
        cx,
    )
}
