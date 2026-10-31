//! Static equilibrium before modal reduction for a tied, heterogeneous shell body.
//!
//! The original CST/DKT preload energy is evaluated in each component's physical
//! coordinates. The existing junction owner's sparse map `T` carries forces by
//! `T^T` and tangents by `T^T K T`, including eccentric and nonmatching seams.
//! Components are never equilibrated independently and then attached afterward.

use super::tethers::{PreparedPoint, PreparedTether, TetherResponse};
use super::{
    ComponentTransform, CoordinateRows, Element, Equilibrium, PreloadOptions, PreloadProblem,
    positive_factor, solve_problem, validate_options,
};
use crate::shell::arching::ReferenceLoadState;
use crate::shell::junction::{
    EdgeTie, FaceTie, ShellComponent, TieBudget, TieKind, TiedShellModel, tie_shells_with_faces,
};
use crate::shell::stiffened::{ShellBeam, assemble_stiffened_shell_materials};
use crate::shell::{ShellMesh, ShellModel, ShellSection, ShellSupport};
use std::cell::Cell;
use std::sync::Arc;

/// One unloaded component with its original sections, supports and dead loads.
#[derive(Clone, Debug)]
pub struct ShellPreloadComponent<'a> {
    /// Reference midsurface, in the common assembly Cartesian frame [m].
    pub mesh: &'a ShellMesh,
    /// One independently specified membrane/bending section per triangle.
    pub sections: &'a [ShellSection],
    /// Essential support nodes in this component's numbering.
    pub boundary: &'a [usize],
    /// Unrestrained components may obtain all their support through the ties.
    pub support: ShellSupport,
    /// Existing reference-frame linear bonded beams.
    pub beams: &'a [ShellBeam],
    /// Full `6 * node_count` world forces [N] and moments [N m].
    /// Loads at eliminated seam nodes are transferred through the same `T^T`.
    pub loads: &'a [f64],
    /// Physical nodal thickness used by the existing junction offset geometry.
    pub nodal_thickness_m: &'a [f64],
    /// An already loaded scan cannot be treated as an unstressed reference.
    pub reference_state: ReferenceLoadState,
}

/// Full assembled equilibrium and the consistent pencil for small vibrations.
#[derive(Debug)]
pub struct PreloadedShellAssembly {
    /// Same tie maps, reference mass and glue as the original assembly; `model.k`
    /// is the stable, consistent tangent at the accepted assembled equilibrium.
    /// Pass this pencil to the existing modal solver, without adding prestress
    /// or glue a second time. Expand modes with `assembly.expand_components`.
    pub assembly: TiedShellModel,
    /// Equilibrium in retained physical coordinates, before solver scaling.
    pub coordinates: Vec<f64>,
    /// Per-component full physical translations [m] and rotations [rad].
    pub displacements: Vec<Vec<f64>>,
    /// Original assembled force imbalance, with moments divided by the span [N].
    pub residual_force_n: f64,
    /// Reference assembly span used only for coordinate/residual norms [m].
    pub norm_length_m: f64,
    /// Shell, bonded-beam and compliant-glue stored energy [J].
    pub stored_energy_j: f64,
    /// Accepted Newton attempts across the original load continuation.
    pub iterations: usize,
    /// Residual/Jacobian evaluations charged by the original preload owner.
    pub evaluations: usize,
}

/// A material bearing on one component, before its original tie transformation.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyTetherPoint {
    /// Canonical component index in the supplied assembly.
    pub component: usize,
    /// Existing facet nodes, in any order.
    pub nodes: [usize; 3],
    /// Nonnegative barycentric weights summing to one.
    pub weights: [f64; 3],
    /// World-frame reference offset from the interpolated midsurface [m].
    pub arm_m: [f64; 3],
}

/// Opposite endpoint of one massless tensile connection.
#[derive(Clone, Copy, Debug)]
pub enum AssemblyTetherAnchor {
    /// Fixed world position [m]; its reaction is opposite the bearing force.
    World([f64; 3]),
    /// A bearing on the same or another moving shell component.
    Shell(AssemblyTetherPoint),
}

/// Original conservative tether law acting between two physical endpoints.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyTether {
    /// First endpoint; response forces are reported on this bearing.
    pub bearing: AssemblyTetherPoint,
    /// Fixed world anchor or a second moving component bearing.
    pub anchor: AssemblyTetherAnchor,
    /// Tension at the unloaded assembly reference geometry [N].
    pub reference_tension_n: f64,
    /// Taut-branch dT/d(length) [N/m]; zero explicitly prescribes tension.
    pub axial_stiffness_n_per_m: f64,
}

/// Joint shell equilibrium and original tensile-connection response ledger.
#[derive(Debug)]
pub struct PreloadedTetheredShellAssembly {
    /// Same assembled equilibrium, with all tether material/geometric stiffness
    /// already included in its tangent. Inertia is unchanged. Do not add these
    /// connections again when using the returned small-vibration pencil.
    pub equilibrium: PreloadedShellAssembly,
    /// Original endpoint forces, tension, length and relative potential, in
    /// input order. `endpoint_forces_n` gives the equal/opposite bearing and
    /// anchor forces. Tether potential changes are separate from the structural
    /// stored energy retained in `equilibrium`.
    pub tether_responses: Vec<TetherResponse>,
}

/// Equilibrate all components together, then linearize that same physical body.
///
/// Uses [`super::equilibrate_shell_materials`]'s Green--Lagrange membrane force,
/// analytic tangent, moderate-rotation limits, continuation and independent
/// physical-force acceptance. The bending, drilling, beams, glue and offsets
/// retain their existing reference-frame linear laws. This is the CST/DKT
/// shell path; it does not silently replace a caller's shear-shell element.
///
/// All component geometries must explicitly be unloaded. This entry point
/// supplies dead loads; [`equilibrate_tethered_shell_assembly`] also includes
/// conservative, direction-updating tensile connections. Nonconservative
/// follower loads and unilateral contact need their own residual and tangent.
/// At most 64 components, 10,000 total nodes and 40,000 total triangles are
/// admitted. Tie work uses the caller's existing `TieBudget`; nonlinear work
/// uses `PreloadOptions`. Cancellation is checked during component/tie assembly
/// and before and after each injected solve. The solver may check it internally.
/// No partially converged body or modified input is published on failure.
///
/// # Errors
/// Loaded reference geometry, invalid component data/budgets, supported rigid
/// slave nodes (uncomposed competing constraints), unconstrained rigid modes,
/// cancellation, excessive deformation, force imbalance or unstable tangent.
pub fn equilibrate_shell_assembly(
    components: &[ShellPreloadComponent<'_>],
    edge_ties: &[EdgeTie],
    face_ties: &[FaceTie],
    tie_budget: TieBudget,
    options: PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize) -> Result<(Vec<f64>, usize), String>,
    cancelled: impl FnMut() -> bool,
) -> Result<PreloadedShellAssembly, String> {
    equilibrate_tethered_shell_assembly(
        components,
        edge_ties,
        face_ties,
        tie_budget,
        &[],
        options,
        solve,
        cancelled,
    )
    .map(|loaded| loaded.equilibrium)
}

/// Solve the same assembled equilibrium with moving, reciprocal tensile ports.
///
/// Each connection uses the original tether potential and its complete material
/// and geometric Hessian. Its relative point motion is composed with both
/// components' exact tie maps before differentiation, retaining seam loads,
/// offset moments and cross-body tangent blocks. Shared retained coordinates
/// are combined once, including two endpoints on the same component. The load
/// fraction scales the entire tether potential together with the dead loads.
///
/// Reference shell stiffness must already be stable through physical supports
/// and ties. Tethers do not regularize an otherwise singular reference body.
/// At most 4,096 connections, 262,144 aggregate retained connection columns and
/// 8,388,608 scalar tangent contributions are admitted. These are cold work caps;
/// no physical spring or missing support is inferred. Other geometry, Newton,
/// cancellation and moderate-deformation limits are the original assembly's.
///
/// # Errors
/// Original assembly refusals, invalid endpoint geometry/laws, coincident
/// endpoints, mapped-work exhaustion, or a tether trial reaching zero length.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_tethered_shell_assembly(
    components: &[ShellPreloadComponent<'_>],
    edge_ties: &[EdgeTie],
    face_ties: &[FaceTie],
    tie_budget: TieBudget,
    tethers: &[AssemblyTether],
    options: PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize) -> Result<(Vec<f64>, usize), String>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreloadedTetheredShellAssembly, String> {
    validate_options(options)?;
    if tethers.len() > 4096 {
        return Err("assembled preload exceeds 4096 tensile connections".into());
    }
    if components.is_empty() || components.len() > 64 {
        return Err("assembled preload requires 1..=64 shell components".into());
    }
    let mut node_count = 0usize;
    let mut triangle_count = 0usize;
    for component in components {
        if cancelled() {
            return Err("assembled shell preload cancelled".into());
        }
        node_count = node_count
            .checked_add(component.mesh.nodes.len())
            .ok_or("assembled preload node budget overflow")?;
        triangle_count = triangle_count
            .checked_add(component.mesh.tris.len())
            .ok_or("assembled preload triangle budget overflow")?;
        if node_count > 10_000 || triangle_count > 40_000 {
            return Err("assembled preload exceeds aggregate shell geometry budget".into());
        }
        if component.reference_state != ReferenceLoadState::Unloaded {
            return Err(
                "assembled preload requires an explicitly unloaded reference for every component"
                    .into(),
            );
        }
        if component.loads.len() != 6 * component.mesh.nodes.len()
            || component.loads.iter().any(|v| !v.is_finite())
        {
            return Err("invalid assembled shell load cardinality or nonfinite force".into());
        }
    }
    if !components
        .iter()
        .any(|c| c.support != ShellSupport::Free && !c.boundary.is_empty())
    {
        return Err(
            "assembled preload needs physical supports; free rigid modes are not regularized"
                .into(),
        );
    }
    check_slave_supports(components, edge_ties, face_ties)?;
    let mut models = Vec::with_capacity(components.len());
    for component in components {
        if cancelled() {
            return Err("assembled shell preload cancelled".into());
        }
        models.push(
            assemble_stiffened_shell_materials(
                component.mesh,
                component.sections,
                component.boundary,
                component.support,
                component.beams,
            )
            .map_err(|e| e.to_string())?,
        );
    }
    let parts: Vec<_> = components
        .iter()
        .zip(&models)
        .enumerate()
        .map(|(id, (c, model))| ShellComponent {
            id,
            mesh: c.mesh,
            model,
            nodal_thickness_m: c.nodal_thickness_m,
        })
        .collect();
    let mut assembly =
        tie_shells_with_faces(&parts, edge_ties, face_ties, tie_budget, &mut cancelled)
            .map_err(|e| e.to_string())?;
    let single = components.len() == 1 && edge_ties.is_empty() && face_ties.is_empty();
    if single {
        // With no junction the original pencil and arithmetic remain exactly
        // the single-shell owner's, including its sparse storage pattern.
        assembly.model = models[0].clone();
    }
    if assembly.model.free == 0 {
        return Err("assembled preload has no unconstrained coordinate".into());
    }
    let base = positive_factor(&assembly.model.k).map_err(|e| {
        format!("assembled reference stiffness has unresolved rigid modes or instability: {e}")
    })?;
    let length = (0..3)
        .map(|axis| {
            let positions = || {
                components
                    .iter()
                    .flat_map(|c| c.mesh.nodes.iter().map(|p| p[axis]))
            };
            positions().fold(f64::NEG_INFINITY, f64::max)
                - positions().fold(f64::INFINITY, f64::min)
        })
        .fold(0.0_f64, f64::max);
    if !length.is_finite() || length <= 0. {
        return Err("invalid assembled shell norm length".into());
    }
    let mut scale = vec![1.; assembly.model.free];
    for (full, free) in assembly.model.dof_map.iter().enumerate() {
        if let Some(i) = free {
            scale[*i] = if full % 6 < 3 { 1. } else { length };
        }
    }
    let load = assembly
        .project_loads(&components.iter().map(|c| c.loads).collect::<Vec<_>>())
        .map_err(|e| e.to_string())?;
    let mut transforms = Vec::with_capacity(components.len());
    let mut elements = Vec::with_capacity(triangle_count);
    for (c, (component, model)) in components.iter().zip(&models).enumerate() {
        // Read the junction owner's exact retained/slave rows. No interpolation
        // or alternate constraint elimination is reconstructed here.
        let mut rows = vec![Vec::new(); model.free];
        for (full, local) in model.dof_map.iter().enumerate() {
            if let (Some(local), Some(global)) = (local, assembly.map[c][full / 6][full % 6]) {
                rows[*local].push((global, 1.));
            }
        }
        for slave in assembly
            .slave_couplings()
            .iter()
            .filter(|s| s.slave_component == c)
        {
            rows[slave.slave_component_free].clone_from(&slave.master_terms);
        }
        let rows: CoordinateRows = Arc::new(rows);
        for (e, section) in component.sections.iter().enumerate() {
            let mut element = Element::new(component.mesh, model, section, e)?;
            if !single {
                element.transform = Some(Arc::clone(&rows));
            }
            elements.push(element);
        }
        transforms.push(ComponentTransform {
            dof_map: model.dof_map.clone(),
            rows,
        });
    }
    let prepared_tethers = prepare_tethers(
        tethers,
        components,
        &models,
        &transforms,
        single,
        &mut cancelled,
    )?;
    let mut problem = PreloadProblem {
        model: &assembly.model,
        elements,
        transforms,
        tethers: prepared_tethers,
        base,
        scale,
        load,
        fraction: 0.,
        options,
        evaluations: Cell::new(0),
    };
    let state = solve_problem(&mut problem, &mut |p, x, limit| {
        if cancelled() {
            return Err("assembled shell preload cancelled".into());
        }
        let solution = solve(p, x, limit)?;
        if cancelled() {
            return Err("assembled shell preload cancelled".into());
        }
        Ok(solution)
    })?;
    let Equilibrium {
        u,
        k,
        energy,
        residual_force_n,
        iterations,
        evaluations,
    } = state;
    let tether_responses = problem
        .tethers
        .iter()
        .map(|t| t.response(&u))
        .collect::<Result<Vec<_>, _>>()?;
    let displacements = assembly.expand_components(&u).map_err(|e| e.to_string())?;
    if cancelled() {
        return Err("assembled shell preload cancelled".into());
    }
    drop(problem);
    assembly.model.k = k;
    let equilibrium = PreloadedShellAssembly {
        assembly,
        coordinates: u,
        displacements,
        residual_force_n,
        norm_length_m: length,
        stored_energy_j: energy,
        iterations,
        evaluations,
    };
    Ok(PreloadedTetheredShellAssembly {
        equilibrium,
        tether_responses,
    })
}

fn prepare_tethers(
    tethers: &[AssemblyTether],
    components: &[ShellPreloadComponent<'_>],
    models: &[ShellModel],
    maps: &[ComponentTransform],
    single: bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<PreparedTether>, String> {
    let point = |port: AssemblyTetherPoint| {
        let component = components
            .get(port.component)
            .ok_or("unknown tether component")?;
        PreparedPoint::new(
            port.nodes,
            port.weights,
            port.arm_m,
            component.mesh,
            &models[port.component],
        )
    };
    let mut prepared = Vec::with_capacity(tethers.len());
    let mut columns = 0usize;
    let mut tangent_entries = 0usize;
    for tether in tethers {
        if cancelled() {
            return Err("assembled shell preload cancelled".into());
        }
        let bearing = point(tether.bearing)?;
        let (anchor_position, anchor) = match tether.anchor {
            AssemblyTetherAnchor::World(position) => (position, None),
            AssemblyTetherAnchor::Shell(port) => {
                let anchor = point(port)?;
                (anchor.position()?, Some((port.component, anchor)))
            }
        };
        let mut connection = PreparedTether::at_point(
            &bearing,
            anchor_position,
            tether.reference_tension_n,
            tether.axial_stiffness_n_per_m,
        )?;
        if !single || anchor.is_some() {
            connection = connection.remap(
                &maps[tether.bearing.component].rows,
                anchor
                    .as_ref()
                    .map(|(component, point)| (point, &maps[*component].rows)),
            )?;
        }
        let count = connection.column_count();
        columns = columns
            .checked_add(count)
            .ok_or("tether column budget overflow")?;
        tangent_entries = tangent_entries
            .checked_add(
                count
                    .checked_mul(count)
                    .ok_or("tether tangent budget overflow")?,
            )
            .ok_or("tether tangent budget overflow")?;
        if columns > 262_144 || tangent_entries > 8_388_608 {
            return Err("assembled tether map or tangent work budget exceeded".into());
        }
        prepared.push(connection);
    }
    Ok(prepared)
}

// Essential supports on a rigid slave cannot simply be skipped during mortar
// elimination: they constrain the master too. Until that constraint composition
// is owned by junction, require the caller to place those supports on the master.
fn check_slave_supports(
    components: &[ShellPreloadComponent<'_>],
    edge_ties: &[EdgeTie],
    face_ties: &[FaceTie],
) -> Result<(), String> {
    let check = |component: usize, node: usize| -> Result<(), String> {
        let part = components
            .get(component)
            .ok_or("unknown shell tie component")?;
        if part.support != ShellSupport::Free && part.boundary.contains(&node) {
            return Err("rigid shell slave supports require composed mortar constraints; place the support on the master".into());
        }
        Ok(())
    };
    for tie in edge_ties {
        if matches!(tie.kind, TieKind::Rigid) {
            for &node in &tie.b.nodes {
                check(tie.b.component, node)?;
            }
        }
    }
    for tie in face_ties {
        if matches!(tie.kind, TieKind::Rigid) {
            let part = components
                .get(tie.slave_component)
                .ok_or("unknown shell tie component")?;
            for &triangle in &tie.slave_tris {
                for &node in part
                    .mesh
                    .tris
                    .get(triangle)
                    .ok_or("unknown shell tie triangle")?
                {
                    check(tie.slave_component, node)?;
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tether_tests;
