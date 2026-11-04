//! Same-step spatial forces from one periodic splitter-plate cell.
//!
//! Patches retain three force components and a common sample clock. Their
//! source positions are actual halfway-reflection geometry, unwrapped around
//! the selected solid cell at a periodic seam. Neither integrating a patch nor
//! partitioning a periodic span establishes a physical isolated radiator.

use super::*;
use core::mem::size_of;
use fs_lbm::d3q19::{BoundaryStepError3, WallPatchProbe3};

/// Authored patch extents in the streamwise and periodic spanwise directions.
/// Both plate faces belong to each x/z patch; the enclosing radius includes
/// plate thickness, leading/trailing faces and diagonal halfway links.
#[derive(Clone, Copy, Debug)]
pub struct SlotJet3dPatchLayout {
    pub cells_x: usize,
    pub cells_z: usize,
}

/// Explicit bounds on the additional spatial source record and probe.
/// Existing flow, fringe and SGS workspaces retain their own allocation rules.
#[derive(Clone, Copy, Debug)]
pub struct SlotJet3dPatchBudget {
    pub maximum_patches: usize,
    /// Force history plus retained patch geometry, excluding allocator overhead.
    pub maximum_record_bytes: usize,
    /// Native wall probe's links, IDs, origins and two exchange buffers.
    pub maximum_probe_bytes: usize,
}

/// A bounding sphere for the selected halfway-reflection points, in cells.
/// It bounds geometry only; a point-force radiation approximation needs a
/// separately authored wavelength/range admission at each observation.
#[derive(Clone, Copy, Debug)]
pub struct SlotJet3dPatch {
    pub center_lu: [f64; 3],
    pub radius_lu: f64,
    pub measured_links: usize,
}

/// Resolved force ON THE FLUID from the splitter's stationary reflections.
/// It is the negative of the plate load, in lattice momentum per step.
/// Patches are ordered z-block then x-block, on the same time grid. The native
/// aggregate/checkpoint record remains a distinct, unchanged API.
#[derive(Debug)]
pub struct SlotJet3dPatchRun {
    pub patches: Vec<SlotJet3dPatch>,
    /// Sample-major: `sample * patches.len() + patch`, three Cartesian axes.
    pub force_on_fluid_lu: Vec<[f64; 3]>,
    pub diagnostics: SlotJet3dDiagnostics,
    pub sgs: Option<SlotJet3dSgsDiagnostics>,
}

fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

fn flow(error: BoundaryStepError3) -> SlotJet3dSgsError {
    BoundarySgsError3::from(error).into()
}

fn allocate<T: Clone>(count: usize, initial: T) -> Result<Vec<T>, AeroacError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(count)
        .map_err(|_| invalid("3-D spatial plate record allocation refused"))?;
    output.resize(count, initial);
    Ok(output)
}

fn admit(
    cfg: &SlotJet3dConfig,
    layout: SlotJet3dPatchLayout,
    budget: SlotJet3dPatchBudget,
) -> Result<(usize, usize), AeroacError> {
    if layout.cells_x == 0 || layout.cells_z == 0 {
        return Err(invalid("3-D plate patch extents must be positive"));
    }
    let patches = cfg
        .plate_length
        .div_ceil(layout.cells_x)
        .checked_mul(cfg.nz.div_ceil(layout.cells_z))
        .filter(|&n| n > 0 && n <= budget.maximum_patches)
        .ok_or_else(|| invalid("3-D plate patch count exceeds budget"))?;
    let forces = patches
        .checked_mul(cfg.steps_record)
        .ok_or_else(|| invalid("3-D plate patch record size overflow"))?;
    forces
        .checked_mul(size_of::<[f64; 3]>())
        .and_then(|bytes| {
            patches
                .checked_mul(size_of::<SlotJet3dPatch>())
                .and_then(|geometry| bytes.checked_add(geometry))
        })
        .filter(|&bytes| bytes <= budget.maximum_record_bytes)
        .ok_or_else(|| invalid("3-D plate patch record exceeds byte budget"))?;
    Ok((patches, forces))
}

fn patch_probe(
    cfg: &SlotJet3dConfig,
    rig: &Rig3,
    layout: SlotJet3dPatchLayout,
    count: usize,
    maximum_probe_bytes: usize,
) -> Result<(Vec<SlotJet3dPatch>, WallPatchProbe3), SlotJet3dSgsError> {
    // Bound these per-link IDs before allocating; the native probe accounts
    // all its retained copies and exchange buffers independently below.
    if rig
        .plate_links
        .len()
        .checked_mul(size_of::<usize>())
        .is_none_or(|bytes| bytes > maximum_probe_bytes)
    {
        return Err(invalid("3-D plate patch selection exceeds probe byte budget").into());
    }
    let mut ids = allocate(rig.plate_links.len(), 0usize)?;
    let mut lower = allocate(count, [f64::INFINITY; 3])?;
    let mut upper = allocate(count, [f64::NEG_INFINITY; 3])?;
    let mut links = allocate(count, 0usize)?;
    let x0 = cfg.nozzle_thickness + cfg.edge_distance;
    let x_count = cfg.plate_length.div_ceil(layout.cells_x);
    for (i, link) in rig.plate_links.iter().enumerate() {
        let e = E3[link.direction];
        let e = [e.0, e.1, e.2];
        let dims = [cfg.nx, cfg.ny, cfg.nz];
        let solid: [usize; 3] = core::array::from_fn(|a| {
            (link.cell[a] as isize - e[a] as isize).rem_euclid(dims[a] as isize) as usize
        });
        let patch = (solid[2] / layout.cells_z) * x_count + (solid[0] - x0) / layout.cells_x;
        ids[i] = patch;
        links[patch] += 1;
        for a in 0..3 {
            // The same halfway link, represented in the solid's periodic cell.
            let position = solid[a] as f64 + 0.5 + 0.5 * f64::from(e[a]);
            lower[patch][a] = lower[patch][a].min(position);
            upper[patch][a] = upper[patch][a].max(position);
        }
    }
    let mut patches = allocate(
        count,
        SlotJet3dPatch {
            center_lu: [0.0; 3],
            radius_lu: 0.0,
            measured_links: 0,
        },
    )?;
    let mut origins = allocate(count, [0.0; 3])?;
    for p in 0..count {
        if links[p] == 0 {
            return Err(invalid("3-D plate patch contains no reflection geometry").into());
        }
        let center = core::array::from_fn(|a| 0.5 * (lower[p][a] + upper[p][a]));
        let radius = det::sqrt(
            (0..3)
                .map(|a| {
                    let half_extent = 0.5 * (upper[p][a] - lower[p][a]);
                    half_extent * half_extent
                })
                .sum(),
        );
        origins[p] = center;
        patches[p] = SlotJet3dPatch {
            center_lu: center,
            radius_lu: radius,
            measured_links: links[p],
        };
    }
    let probe = rig
        .grid
        .wall_patch_probe(&rig.plate_links, &ids, &origins, maximum_probe_bytes)
        .map_err(flow)?;
    Ok((patches, probe))
}

/// Run spatially resolved plate forces with an optional explicitly authored SGS
/// model. The original central-moment solver, fringe, seed, molecular viscosity
/// and stability/speed diagnostics remain the owners of the trajectory.
///
/// Every record is measured before fringe forcing from the actual collision
/// that streams that step. Patches preserve resolved spatial phase; they do not
/// supply unresolved quadrupoles, artificial decorrelation or stochastic noise.
/// No patch is presumed compact, and no span replication or area weighting is
/// inferred. Source-record budgets are admitted before constructing the grid.
///
/// # Errors
/// Original source, closure, native step, memory and cancellation refusals
/// return no partial record. Existing aggregate/checkpoint APIs are unaffected.
pub fn run_slot_jet_3d_patches(
    cfg: &SlotJet3dConfig,
    layout: SlotJet3dPatchLayout,
    budget: SlotJet3dPatchBudget,
    closure: Option<SgsConfig>,
    gate: &CancelGate,
) -> Result<SlotJet3dPatchRun, SlotJet3dSgsError> {
    validate_config(cfg)?;
    let (patch_count, force_count) = admit(cfg, layout, budget)?;
    if let Some(config) = closure {
        config.validate().map_err(BoundarySgsError3::from)?;
    }
    if gate.is_requested() {
        return Err(AeroacError::Cancelled {
            what: "3-D spatial plate record",
        }
        .into());
    }
    let mut forces = Vec::new();
    forces
        .try_reserve_exact(force_count)
        .map_err(|_| invalid("3-D spatial plate record allocation refused"))?;
    let mut rig = Rig3::build(cfg);
    let (patches, mut probe) =
        patch_probe(cfg, &rig, layout, patch_count, budget.maximum_probe_bytes)?;
    let mut workspace = closure.map(|_| rig.grid.sgs_workspace()).transpose()?;
    std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
        for _ in 0..cfg.steps_settle {
            if let (Some(config), Some(workspace)) = (closure, &mut workspace) {
                rig.step_sgs(gate, config, workspace, false)?;
            } else {
                if gate.is_requested() {
                    return Err(AeroacError::Cancelled {
                        what: "3-D spatial jet settle",
                    }
                    .into());
                }
                rig.grid.step();
                rig.fringe.apply(&mut rig.grid);
            }
        }
        rig.check_stable()?;
        let mut recorded = SgsReport::default();
        for _ in 0..cfg.steps_record {
            if let (Some(config), Some(workspace)) = (closure, &mut workspace) {
                let report = rig
                    .grid
                    .step_with_sgs_and_wall_patches(gate, config, workspace, &mut probe)?;
                recorded.fluid_cells = report.fluid_cells;
                recorded.resolved_dissipation_lu +=
                    report.resolved_dissipation_lu / cfg.steps_record as f64;
                recorded.sgs_dissipation_lu += report.sgs_dissipation_lu / cfg.steps_record as f64;
                recorded.maximum_eddy_viscosity_lu = recorded
                    .maximum_eddy_viscosity_lu
                    .max(report.maximum_eddy_viscosity_lu);
                recorded.maximum_eddy_to_molecular = recorded
                    .maximum_eddy_to_molecular
                    .max(report.maximum_eddy_to_molecular);
            } else {
                rig.grid
                    .step_with_wall_patches(gate, &mut probe)
                    .map_err(flow)?;
            }
            forces.extend(
                probe
                    .exchanges()
                    .expect("successful native patch step")
                    .iter()
                    .map(|exchange| exchange.fluid_population_impulse),
            );
            rig.fringe.apply(&mut rig.grid);
        }
        if !recorded.resolved_dissipation_lu.is_finite() || !recorded.sgs_dissipation_lu.is_finite()
        {
            return Err(AeroacError::NonFinite {
                what: "3-D spatial jet mean shear dissipation",
            }
            .into());
        }
        rig.check_stable()?;
        let (speed, plate, fringe) = rig.plane_diagnostics();
        let diagnostics = finish_run(cfg, Vec::new(), speed, plate, fringe).diagnostics;
        Ok(SlotJet3dPatchRun {
            patches,
            force_on_fluid_lu: forces,
            diagnostics,
            sgs: closure.map(|config| SlotJet3dSgsDiagnostics { config, recorded }),
        })
    }))
    .unwrap_or_else(|_| {
        Err(AeroacError::NonFinite {
            what: "3-D spatial jet destabilized during settle or record",
        }
        .into())
    })
}

#[cfg(test)]
mod tests;
