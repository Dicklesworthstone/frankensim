//! One material body assembled from independently oriented component surfaces.
//!
//! Each part keeps its own material, thickness, directors and fluid roles.
//! The original mortar owner ties the complete component pencils before one
//! shared eigensolve. Forces, pressure and every viscous cross term observe
//! those same modes. No component modes are truncated before assembly.
//!
//! This is unloaded linear mechanics. Rigid seams do not imply glue compliance,
//! assembled prestress, a closed fluid skin or resolved interior/exterior air.

use super::{
    AcousticRealizeError, ResolvedShearShell, ShellProjectionBudget, invalid, owner, poll,
    shell_mode_budget_admits,
};
use crate::plate_viscosity::ModalPlateViscosity;
use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection};
use crate::thin_plate::CompactBody;
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_plate::shell::{
    junction::{
        EdgeTie, FaceTie, RotationCoupling, ShellComponent, TieBudget, TieKind, TieReceipt,
        tie_shells_with_faces,
    },
    mitc::{
        ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass, assemble_mitc3p,
    },
    modes_shell,
};

/// One independently oriented material midsurface, before ties are eliminated.
#[derive(Clone, Debug)]
pub struct MaterialShellComponent {
    /// Unique nonempty part identity; never used to select a physical law.
    pub name: String,
    /// Original regions and explicit exterior/cavity facet selections.
    pub specimen: ResolvedShearShell,
    /// Directors belong to this smooth component, never averaged across seams.
    pub directors: NodalDirectors,
    /// Homogeneous GLOBAL supports. Local frames await support in the tie owner.
    pub constraints: Vec<ShellConstraint>,
    /// Original translational and rotary mass discretization.
    pub mass: ShellMass,
    /// Positive numerical drilling stabilization, never material damping.
    pub drilling_alpha: f64,
    /// Original per-component assembly extent.
    pub assembly_budget: MitcBudget,
}

/// Complete component sources and work-conjugate seams, in canonical part order.
#[derive(Clone, Debug)]
pub struct MaterialShellAssembly {
    /// Independently sourced top, back, ribs or other generic structural parts.
    pub components: Vec<MaterialShellComponent>,
    /// Original nonmatching edge ties; indices address `components`.
    pub edge_ties: Vec<EdgeTie>,
    /// Original overlapping face ties; indices address `components`.
    pub face_ties: Vec<FaceTie>,
    /// Explicit mortar work and geometric admission limits.
    pub tie_budget: TieBudget,
}

/// One physical force/moment footprint distributed over any number of parts.
#[derive(Clone, Debug)]
pub struct ShellAssemblyForcePort {
    /// One GLOBAL six-DOF row per component, in original node order. Translation
    /// entries are force fractions; rotation entries are moment/N arms [m].
    /// Empty participation is a correctly sized zero row, never a missing part.
    pub components: Vec<Vec<f64>>,
}

/// Whole-body spectral and projection controls, independent of part materials.
#[derive(Clone, Debug)]
pub struct ShellAssemblyReductionOptions {
    /// Original angular-frequency-squared interval (low,high] [rad²/s²].
    pub eigenvalue_window: (f64, f64),
    /// Lowest modes retained AFTER the complete tied body has been assembled.
    pub n_modes: usize,
    /// Complete force, pressure and full-loss projection allowance.
    pub projection_budget: ShellProjectionBudget,
    /// Cumulative physical nodes across every component, before ties.
    pub max_total_nodes: usize,
    /// Cumulative physical triangles across every component.
    pub max_total_triangles: usize,
}

/// One assembled basis and all work-conjugate material and acoustic projections.
#[derive(Clone, Debug)]
pub struct MaterialShellAssemblyImage {
    source: MaterialShellAssembly,
    physical_modes: Vec<Vec<Vec<f64>>>,
    force_ports: Vec<Vec<f64>>,
    volume_shapes: Vec<f64>,
    cavity_volume_shapes: Vec<f64>,
    intervals: Vec<(f64, f64)>,
    radiators: Vec<CompactBody>,
    damping: ModalPlateViscosity,
    interface_columns: Vec<Vec<Vec<f64>>>,
    tie_receipt: TieReceipt,
}

impl MaterialShellAssemblyImage {
    /// Exact original geometry, material records, fluid roles and tie inputs.
    pub const fn source(&self) -> &MaterialShellAssembly {
        &self.source
    }
    /// Component-major, then mode-major GLOBAL six-DOF physical vectors. Every
    /// component uses the SAME mass-normalized assembled eigenvector signs.
    pub fn physical_modes(&self) -> &[Vec<Vec<f64>>] {
        &self.physical_modes
    }
    /// Work-conjugate modal force/displacement rows in caller port order.
    pub fn force_ports(&self) -> &[Vec<f64>] {
        &self.force_ports
    }
    /// Sum of explicitly exterior-facing normal-volume rows [m²/sqrt(kg)].
    pub fn volume_shapes_m2_per_sqrt_kg(&self) -> &[f64] {
        &self.volume_shapes
    }
    /// Sum of outward-from-cavity rows [m²/sqrt(kg)]; one common compression port.
    pub fn cavity_volume_shapes_m2_per_sqrt_kg(&self) -> &[f64] {
        &self.cavity_volume_shapes
    }
    /// Original assembled eigenvalue residual intervals, not continuum bounds.
    pub fn eigenvalue_intervals(&self) -> &[(f64, f64)] {
        &self.intervals
    }
    /// Complete sum of component modal viscosity matrices, including cross terms.
    pub const fn damping(&self) -> &ModalPlateViscosity {
        &self.damping
    }
    /// Component-major material-interface rows, in each component's source order.
    pub fn interface_columns(&self) -> &[Vec<Vec<f64>>] {
        &self.interface_columns
    }
    /// Original mortar geometry/work receipt.
    pub const fn tie_receipt(&self) -> &TieReceipt {
        &self.tie_receipt
    }
    pub(crate) fn radiators(&self) -> &[CompactBody] {
        &self.radiators
    }
    pub(crate) fn interface_damping_connections(&self, component: usize) -> Vec<ModalConnection> {
        self.source
            .components
            .iter()
            .zip(&self.interface_columns)
            .flat_map(|(part, columns)| part.specimen.interfaces.iter().zip(columns))
            .filter(|(source, _)| source.law().damping_n_s_m != 0.0)
            .map(|(source, column)| ModalConnection {
                left: ModalAttachment {
                    component,
                    shapes: column.clone(),
                },
                right: ModalAttachment {
                    component,
                    shapes: vec![0.0; column.len()],
                },
                stiffness_n_m: 0.0,
                damping_n_s_m: source.law().damping_n_s_m,
                rest_extension_m: 0.0,
            })
            .collect()
    }
}

impl MaterialShellAssembly {
    /// Whether a compliant junction supplies positive viscous loss. No
    /// coefficient is inferred from glue stiffness or a material name.
    pub fn has_viscous_glue(&self) -> bool {
        self.edge_ties
            .iter()
            .map(|tie| tie.kind)
            .chain(self.face_ties.iter().map(|tie| tie.kind))
            .any(|kind| {
                matches!(kind, TieKind::Compliant(glue)
                if glue.viscosity_n_s_m3.is_some_and(|eta| eta.iter().any(|v| *v > 0.0)))
            })
    }

    /// Replace a complete regional material source while keeping all geometry,
    /// force/fluid roles and junction ownership. Actual thickness is re-derived
    /// from the source policy during reduction, including fixed-mass changes.
    pub fn with_region_material(
        &self,
        component: &str,
        region: &str,
        material: &ResolvedMaterialStatePoint,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let matches: Vec<_> = self
            .components
            .iter()
            .enumerate()
            .filter(|(_, part)| part.name == component)
            .map(|(i, _)| i)
            .collect();
        if matches.len() != 1 {
            return Err(invalid("assembly material swap needs one named component"));
        }
        let mut result = self.clone();
        let part = &mut result.components[matches[0]];
        part.specimen = part.specimen.with_region_material(region, material, gate)?;
        poll(gate)?;
        Ok(result)
    }

    /// Assemble full source-bound components and ties BEFORE the one spectral
    /// slice. All force, fluid and loss terms then use that same physical basis.
    /// No fitted frequency, diagonal loss replacement or per-part truncation.
    ///
    /// Explicit pressure ownership is required on every part; an internal part
    /// may select neither fluid. Adjacent triangles of a component must give
    /// each shared node one thickness; represent thickness steps as separate
    /// parts so attachment skin offsets cannot be silently averaged.
    ///
    /// # Errors
    /// Invalid/ambiguous sources, local support frames, constrained slave tie
    /// coordinates, physical-port or cumulative-work limits, material bands,
    /// original shell/mortar/spectral refusals or cancellation.
    pub fn reduce(
        &self,
        options: &ShellAssemblyReductionOptions,
        force_ports: &[ShellAssemblyForcePort],
        gate: &CancelGate,
    ) -> Result<MaterialShellAssemblyImage, AcousticRealizeError> {
        poll(gate)?;
        let n = options.n_modes;
        let budget = options.projection_budget;
        if self.components.is_empty()
            || self.components.len() > 64
            || !shell_mode_budget_admits(n, budget.max_modes)
            || !budget.reference_resistance_n_s_m.is_finite()
            || budget.reference_resistance_n_s_m <= 0.0
            || !options.eigenvalue_window.0.is_finite()
            || !options.eigenvalue_window.1.is_finite()
            || options.eigenvalue_window.0 < 0.0
            || options.eigenvalue_window.1 <= options.eigenvalue_window.0
            || force_ports.is_empty()
            || force_ports.len() > 4096
        {
            return Err(invalid(
                "shell assembly requires bounded parts, modes, ports, window and resistance",
            ));
        }
        let mut names = std::collections::BTreeSet::new();
        let mut total_nodes = 0usize;
        let mut total_triangles = 0usize;
        let mut work = Some(4 * n * n * n);
        for part in &self.components {
            poll(gate)?;
            if part.name.trim().is_empty()
                || !names.insert(part.name.as_str())
                || part
                    .constraints
                    .iter()
                    .any(|c| c.frame != ConstraintFrame::Global)
            {
                return Err(invalid(
                    "shell assembly needs unique part names and global support frames",
                ));
            }
            let source = &part.specimen;
            let nodes = source.mesh.nodes.len();
            let triangles = source.mesh.tris.len();
            total_nodes = total_nodes
                .checked_add(nodes)
                .ok_or_else(|| invalid("assembly node extent overflow"))?;
            total_triangles = total_triangles
                .checked_add(triangles)
                .ok_or_else(|| invalid("assembly triangle extent overflow"))?;
            if total_nodes > options.max_total_nodes
                || total_triangles > options.max_total_triangles
            {
                return Err(invalid(
                    "shell assembly exceeds its cumulative geometry allowance",
                ));
            }
            let surfaces = source.pressure_surfaces().ok_or_else(|| {
                invalid(
                    "every assembled shell part requires explicit exterior and cavity ownership",
                )
            })?;
            work = work
                .and_then(|w| source.loss_projection_terms(n)?.checked_add(w))
                .and_then(|w| {
                    6usize
                        .checked_mul(nodes)?
                        .checked_mul(n)?
                        .checked_mul(force_ports.len())?
                        .checked_add(w)
                })
                .and_then(|w| {
                    fs_plate::shell::pressure::projection_terms(nodes, surfaces.exterior.len(), n)?
                        .checked_add(w)
                })
                .and_then(|w| {
                    fs_plate::shell::pressure::projection_terms(nodes, surfaces.cavity.len(), n)?
                        .checked_add(w)
                });
        }
        // A conservative dense upper bound covers recovery and dual force work
        // through the mortar maps, even when their actual rows are sparse.
        work = work.and_then(|w| {
            total_nodes
                .checked_mul(6)?
                .checked_mul(total_nodes.checked_mul(6)?)?
                .checked_mul(n.checked_add(force_ports.len())?)?
                .checked_add(w)
        });
        // Project the original glue C in the same basis. Bound n sparse
        // products densely and include all n^2 modal inner products.
        if self.has_viscous_glue() {
            work = work.and_then(|w| {
                let dofs = total_nodes.checked_mul(6)?;
                dofs.checked_mul(dofs)?
                    .checked_mul(n)?
                    .checked_add(dofs.checked_mul(n)?.checked_mul(n)?)?
                    .checked_add(w)
            });
        }
        if work.is_none_or(|w| w > budget.max_projection_terms) {
            return Err(invalid(
                "complete assembly projection exceeds its declared allowance",
            ));
        }
        for port in force_ports {
            if port.components.len() != self.components.len()
                || port
                    .components
                    .iter()
                    .zip(&self.components)
                    .any(|(row, part)| {
                        row.len() != 6 * part.specimen.mesh.nodes.len()
                            || row.iter().any(|v| !v.is_finite())
                    })
                || !port.components.iter().flatten().any(|v| *v != 0.0)
            {
                return Err(invalid(
                    "assembly force port requires finite physical rows for every part",
                ));
            }
        }
        self.check_slave_supports()?;
        let mut models = Vec::with_capacity(self.components.len());
        let mut thickness = Vec::with_capacity(self.components.len());
        for part in &self.components {
            poll(gate)?;
            let source = &part.specimen;
            thickness.push(nodal_thickness(source)?);
            let model = assemble_mitc3p(
                &source.mesh,
                &part.directors,
                &source.sections,
                &part.constraints,
                part.mass,
                part.assembly_budget,
                part.drilling_alpha,
                || gate.is_requested(),
            )
            .map_err(owner)?;
            let laws: Vec<_> = source.interfaces.iter().map(|i| i.law().clone()).collect();
            models.push(
                model
                    .with_connections(&laws, source.interface_budget, || gate.is_requested())
                    .map_err(owner)?,
            );
        }
        let parts: Vec<_> = self
            .components
            .iter()
            .enumerate()
            .map(|(id, part)| ShellComponent {
                id,
                mesh: &part.specimen.mesh,
                model: models[id].pencil(),
                nodal_thickness_m: &thickness[id],
            })
            .collect();
        let tied = tie_shells_with_faces(
            &parts,
            &self.edge_ties,
            &self.face_ties,
            self.tie_budget,
            || gate.is_requested(),
        )
        .map_err(owner)?;
        let spectral = fs_modal::SliceOptions {
            mass_diagonal_equilibration: true,
            ..Default::default()
        };
        poll(gate)?;
        let report =
            modes_shell(&tied.model, options.eigenvalue_window, &spectral).map_err(owner)?;
        poll(gate)?;
        if report.modes.len() < n {
            // The count is only known after the whole window is solved: name it, so the next request can ask for it.
            return Err(owner(format_args!(
                "assembled shell window has fewer modes than requested: it holds {}, {n} requested",
                report.modes.len()
            )));
        }
        let modes = &report.modes[..n];
        let mut physical_modes = vec![Vec::with_capacity(n); self.components.len()];
        let mut intervals = Vec::with_capacity(n);
        let mut omega = Vec::with_capacity(n);
        for mode in modes {
            poll(gate)?;
            let w = fs_math::det::sqrt(mode.lambda);
            if !w.is_finite()
                || w <= 0.0
                || self.components.iter().any(|part| {
                    part.specimen
                        .regions
                        .iter()
                        .any(|r| w < r.omega_band.0 || w > r.omega_band.1)
                        || part.specimen.interfaces.iter().any(|r| {
                            let band = r.omega_band_rad_s();
                            w < band.0 || w > band.1
                        })
                })
            {
                return Err(invalid(
                    "assembled shell mode lies outside original component material support",
                ));
            }
            let component_free = tied.expand_component_free(&mode.phi).map_err(owner)?;
            for (i, free) in component_free.iter().enumerate() {
                physical_modes[i].push(models[i].expand(free).map_err(owner)?);
            }
            intervals.push(mode.interval);
            omega.push(w);
        }
        let mut projected = vec![vec![0.0; n]; force_ports.len()];
        for (row, input) in projected.iter_mut().zip(force_ports) {
            poll(gate)?;
            let inputs: Vec<_> = input.components.iter().map(Vec::as_slice).collect();
            let reduced = tied.project_loads(&inputs).map_err(owner)?;
            for (j, mode) in modes.iter().enumerate() {
                let mut value = 0.0;
                let mut scale = 0.0;
                for (i, full) in input.components.iter().enumerate() {
                    for (&force, &displacement) in full.iter().zip(&physical_modes[i][j]) {
                        value += force * displacement;
                        scale += (force * displacement).abs();
                    }
                }
                let expected: f64 = reduced.iter().zip(&mode.phi).map(|(a, b)| a * b).sum();
                if !value.is_finite()
                    || !scale.is_finite()
                    || (value - expected).abs() > 1e-11 * scale.max(1.0)
                {
                    return Err(invalid(
                        "assembled physical force projection lost work conjugacy",
                    ));
                }
                row[j] = value;
            }
            if row.iter().all(|v| *v == 0.0) {
                return Err(invalid(
                    "assembly force footprint is detached from retained modes",
                ));
            }
        }
        let mut volume_shapes = vec![0.0; n];
        let mut cavity_volume_shapes = vec![0.0; n];
        let mut matrix = vec![0.0; n * n];
        if self.has_viscous_glue() {
            poll(gate)?;
            let basis: Vec<_> = modes.iter().map(|mode| mode.phi.clone()).collect();
            add(
                &mut matrix,
                &tied.project_glue_damping(&basis).map_err(owner)?,
            )?;
        }
        let mut interface_columns = Vec::with_capacity(self.components.len());
        let mut exterior_area = 0.0;
        for (i, part) in self.components.iter().enumerate() {
            poll(gate)?;
            let source = &part.specimen;
            let pressure = source.pressure_ports(gate)?;
            add(
                &mut volume_shapes,
                &pressure
                    .exterior
                    .project_modes(&physical_modes[i], budget.max_projection_terms, || {
                        gate.is_requested()
                    })
                    .map_err(owner)?,
            )?;
            add(
                &mut cavity_volume_shapes,
                &pressure
                    .cavity
                    .project_modes(&physical_modes[i], budget.max_projection_terms, || {
                        gate.is_requested()
                    })
                    .map_err(owner)?,
            )?;
            exterior_area += pressure.exterior.area_m2();
            add(&mut matrix, &source.project_material_viscosity(&part.directors,
                &physical_modes[i], &omega, part.assembly_budget,
                budget.max_projection_terms, gate)?)?;
            let laws: Vec<_> = source.interfaces.iter().map(|i| i.law().clone()).collect();
            interface_columns.push(
                models[i]
                    .project_connections(&laws, &physical_modes[i], source.interface_budget, || {
                        gate.is_requested()
                    })
                    .map_err(owner)?,
            );
        }
        // Localized adhesive loss can have an undamped subspace. Reuse the
        // existing PSD channel owner without adding artificial diagonal loss.
        let damping = if self.has_viscous_glue() {
            ModalPlateViscosity::from_psd_modal_matrix(
                matrix,
                n,
                budget.reference_resistance_n_s_m,
                gate,
            )?
        } else {
            ModalPlateViscosity::factor(matrix, n, budget.reference_resistance_n_s_m, gate)?
        };
        let radiators = omega
            .iter()
            .enumerate()
            .map(|(i, &w)| {
                CompactBody::from_mass_normalized_projection(
                    w,
                    volume_shapes[i],
                    projected[0][i],
                    exterior_area,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        poll(gate)?;
        Ok(MaterialShellAssemblyImage {
            source: self.clone(),
            physical_modes,
            force_ports: projected,
            volume_shapes,
            cavity_volume_shapes,
            intervals,
            radiators,
            damping,
            interface_columns,
            tie_receipt: tied.receipt,
        })
    }

    // The current mortar owner eliminates free slave coordinates. It cannot
    // propagate a prescribed slave support back to its masters, so such input
    // must refuse rather than silently dropping one of the physical conditions.
    fn check_slave_supports(&self) -> Result<(), AcousticRealizeError> {
        let check = |component: usize, node: usize, rotations: RotationCoupling| {
            let part = self
                .components
                .get(component)
                .ok_or_else(|| invalid("unknown shell tie component"))?;
            let dofs = if rotations == RotationCoupling::Full {
                6
            } else {
                3
            };
            if part
                .constraints
                .iter()
                .any(|c| c.node == node && c.mask[..dofs].iter().any(|v| *v))
            {
                return Err(invalid(
                    "shell tie slave supports require original mortar constraint composition",
                ));
            }
            Ok(())
        };
        for tie in &self.edge_ties {
            if !matches!(tie.kind, TieKind::Rigid) {
                continue;
            }
            for &node in &tie.b.nodes {
                check(tie.b.component, node, tie.rotations)?;
            }
        }
        for tie in &self.face_ties {
            if !matches!(tie.kind, TieKind::Rigid) {
                continue;
            }
            let part = self
                .components
                .get(tie.slave_component)
                .ok_or_else(|| invalid("unknown shell tie component"))?;
            for &triangle in &tie.slave_tris {
                let nodes = part
                    .specimen
                    .mesh
                    .tris
                    .get(triangle)
                    .ok_or_else(|| invalid("unknown shell tie triangle"))?;
                for &node in nodes {
                    check(tie.slave_component, node, tie.rotations)?;
                }
            }
        }
        Ok(())
    }
}

fn nodal_thickness(source: &ResolvedShearShell) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut values = vec![0.0; source.mesh.nodes.len()];
    for (element, tri) in source.mesh.tris.iter().enumerate() {
        let h = source.regions[source.element_regions[element]].thickness_m;
        for &node in tri {
            if values[node] != 0.0 && values[node] != h {
                return Err(invalid(
                    "shell component has ambiguous nodal thickness; split the thickness step into parts",
                ));
            }
            values[node] = h;
        }
    }
    Ok(values)
}

fn add(target: &mut [f64], input: &[f64]) -> Result<(), AcousticRealizeError> {
    if target.len() != input.len() {
        return Err(invalid(
            "assembled physical projections have different extents",
        ));
    }
    for (target, value) in target.iter_mut().zip(input) {
        *target += value;
        if !target.is_finite() {
            return Err(invalid("assembled physical projection overflows"));
        }
    }
    Ok(())
}
