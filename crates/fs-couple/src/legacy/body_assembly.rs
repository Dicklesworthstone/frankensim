//! Joint thermal/dead-load equilibrium and native membrane prestress of MITC shells.
//!
//! The original elastic shell pencils and thermal condensed loads are assembled
//! through the existing junction transformation. Equilibrium is solved once in
//! that shared space; each original thermal owner then recovers the SAME physical
//! displacement and supplies its signed planar membrane initial-stress tangent.
//! The original junction rows assemble that tangent without rebuilding ties or
//! changing mass. Modal windows belong to `fs-modal`.
//!
//! This is the native flat-reference, small-displacement initial-stress model,
//! with fixed linear ties, homogeneous global supports and dead nodal loads.
//! It is not a geometrically nonlinear body solve: curved component prestress,
//! follower loads, support settlement, contact, beam/solid initial stress and
//! post-buckling require their respective native owners. Every source geometry
//! must be explicitly unloaded; a loaded scan cannot receive preload twice.
//! Static displacement establishes a reference, never a time-step work credit.

use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_modal::{MultiSliceReport, SliceOptions, WindowPlan, slice_windows};
use fs_plate::PlateError;
use fs_plate::shell::arching::ReferenceLoadState;
use fs_plate::shell::junction::{
    EdgeTie, FaceTie, RotationCoupling, ShellComponent, TieBudget, TieKind, TiedShellModel,
    tie_shells_with_faces,
};
use fs_plate::shell::mitc::PreparedShellStatic;
use fs_plate::shell::mitc::thermal::PreparedThermalShell;
use fs_plate::shell::mitc::thermal::prestress::{
    MembranePrestressLimits, MembranePrestressedShell,
};
use fs_sparse::Csr;

/// One immutable thermal/material owner in its original unloaded coordinates.
#[derive(Clone, Copy)]
pub struct ThermalShellComponent<'a> {
    pub thermal: &'a PreparedThermalShell<'a>,
    /// Complete world force/moment vector, including entries on supports [N,Nm].
    /// Thermal eigenstrain loads are supplied by `thermal`, exactly once.
    pub dead_load_6n: &'a [f64],
    /// Original physical nodal thickness for junction eccentricity [m].
    pub nodal_thickness_m: &'a [f64],
    pub reference_state: &'a ReferenceLoadState,
}

/// Explicit input limits; the original sparse factor and tie kernel retain
/// their own allocation and cancellation boundaries.
#[derive(Debug, Clone, Copy)]
pub struct BodyPrestressBudget {
    pub maximum_components: usize,
    pub maximum_total_nodes: usize,
    pub maximum_component_free_coordinates: usize,
    pub ties: TieBudget,
}

/// Original-dimensional equilibrium admission in addition to native relative
/// iterative-refinement checks. Moments are divided by the explicit length.
#[derive(Debug, Clone, Copy)]
pub struct BodyEquilibriumTolerance {
    pub maximum_residual_force_n: f64,
    pub moment_length_m: f64,
}

/// Accepted original equilibrium and its signed, work-conjugate loaded pencil.
/// Private fields keep equilibrium, component recovery and junction maps paired.
#[derive(Debug)]
pub struct PrestressedShellBody<'a> {
    assembly: TiedShellModel,
    elastic_stiffness: Csr,
    components: Vec<MembranePrestressedShell<'a>>,
    coordinates: Vec<f64>,
    component_free_residuals: Vec<Vec<f64>>,
    residual_force_n: f64,
}

impl PrestressedShellBody<'_> {
    /// Native `(K_material + K_G + K_glue, M_reference)` and unchanged tie maps.
    /// Negative stiffness is retained for stability analysis, never clamped.
    pub const fn assembly(&self) -> &TiedShellModel {
        &self.assembly
    }
    pub const fn elastic_stiffness(&self) -> &Csr {
        &self.elastic_stiffness
    }
    pub fn components(&self) -> &[MembranePrestressedShell<'_>] {
        &self.components
    }
    pub fn coordinates(&self) -> &[f64] {
        &self.coordinates
    }
    pub const fn residual_force_n(&self) -> f64 {
        self.residual_force_n
    }

    /// Each component's `K_c q_c - f_thermal,c - f_dead,c` in original free
    /// coordinates, before junction forces. Nonzero entries are balanced by
    /// rigid/compliant interfaces; these are not essential-support reactions.
    pub fn component_free_residuals(&self) -> &[Vec<f64>] {
        &self.component_free_residuals
    }

    /// Equilibrium physical displacement/rotation on each original component.
    pub fn displacements(&self) -> Result<Vec<Vec<f64>>, PlateError> {
        self.assembly.expand_components(&self.coordinates)
    }

    /// Original signed spectral windows. The report includes inertia counts
    /// and residual certificates; positive windows alone do not prove stability
    /// outside those windows. No absolute value or real-frequency conversion is
    /// applied to a zero/negative eigenvalue.
    pub fn eigen_windows<'b>(
        &'b self,
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<PrestressedBodyModes<'b, '_>, PlateError> {
        let model = &self.assembly.model;
        let spectrum =
            slice_windows(&model.k, &model.m, plan, options, cx).map_err(PlateError::Modal)?;
        Ok(PrestressedBodyModes {
            body: self,
            spectrum,
        })
    }
}

/// Spectral output borrowing the exact accepted body and physical-coordinate map.
#[derive(Debug)]
pub struct PrestressedBodyModes<'body, 'source> {
    body: &'body PrestressedShellBody<'source>,
    spectrum: MultiSliceReport,
}

impl PrestressedBodyModes<'_, '_> {
    pub const fn spectrum(&self) -> &MultiSliceReport {
        &self.spectrum
    }

    /// Native mass-normalized modal perturbation; it excludes the static offset.
    pub fn expand_mode(&self, window: usize, mode: usize) -> Result<Vec<Vec<f64>>, PlateError> {
        let mode = self
            .spectrum
            .reports
            .get(window)
            .and_then(|r| r.modes.get(mode))
            .ok_or_else(|| bad("assembled modal index out of range"))?;
        self.body.assembly.expand_components(&mode.phi)
    }

    /// Generalized modal force for each report/mode, through the SAME physical
    /// tie map as equilibrium and expansion: `phi^T T^T f`.
    pub fn project_loads(&self, component_loads: &[&[f64]]) -> Result<Vec<Vec<f64>>, PlateError> {
        let force = self.body.assembly.project_loads(component_loads)?;
        self.spectrum
            .reports
            .iter()
            .map(|report| {
                report
                    .modes
                    .iter()
                    .map(|mode| {
                        let value = mode
                            .phi
                            .iter()
                            .zip(&force)
                            .fold(Dd::ZERO, |sum, (&a, &b)| {
                                sum + Dd::from_f64(a) * Dd::from_f64(b)
                            })
                            .to_f64();
                        if !value.is_finite() {
                            return Err(bad("assembled modal force overflow"));
                        }
                        Ok(value)
                    })
                    .collect()
            })
            .collect()
    }
}

/// Solve joint native thermal/dead-load equilibrium, then assemble the original
/// membrane initial-stress tangents at that accepted state. No component is
/// independently equilibrated before tying. Disconnected supported components
/// are allowed; unresolved rigid mechanisms retain the static owner's refusal.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_thermal_shell_body<'a>(
    components: &[ThermalShellComponent<'a>],
    edge_ties: &[EdgeTie],
    face_ties: &[FaceTie],
    budget: BodyPrestressBudget,
    tolerance: BodyEquilibriumTolerance,
    limits: MembranePrestressLimits,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrestressedShellBody<'a>, PlateError> {
    poll(&mut cancelled)?;
    if components.is_empty()
        || !tolerance.maximum_residual_force_n.is_finite()
        || tolerance.maximum_residual_force_n <= 0.0
        || !tolerance.moment_length_m.is_finite()
        || tolerance.moment_length_m <= 0.0
    {
        return Err(bad(
            "nonempty body and finite positive equilibrium force/length tolerances required",
        ));
    }
    cap(
        "body components",
        components.len(),
        budget.maximum_components,
    )?;
    let (mut nodes, mut free) = (0usize, 0usize);
    for part in components {
        poll(&mut cancelled)?;
        let model = part.thermal.model();
        if *part.reference_state != ReferenceLoadState::Unloaded {
            return Err(bad(
                "assembled preload requires an explicitly unloaded reference",
            ));
        }
        if !model.has_global_support_coordinates() {
            return Err(bad(
                "assembled preload requires original global support coordinates",
            ));
        }
        nodes = nodes
            .checked_add(model.original_geometry().0.nodes.len())
            .ok_or_else(|| bad("body node count overflow"))?;
        free = free
            .checked_add(model.pencil().free)
            .ok_or_else(|| bad("body coordinate count overflow"))?;
        cap("body nodes", nodes, budget.maximum_total_nodes)?;
        cap(
            "body component coordinates",
            free,
            budget.maximum_component_free_coordinates,
        )?;
        model.project_load(part.dead_load_6n)?;
    }
    check_slave_supports(components, edge_ties, face_ties)?;
    let parts: Vec<_> = components
        .iter()
        .enumerate()
        .map(|(id, part)| ShellComponent {
            id,
            mesh: part.thermal.model().original_geometry().0,
            model: part.thermal.model().pencil(),
            nodal_thickness_m: part.nodal_thickness_m,
        })
        .collect();
    let mut assembly =
        tie_shells_with_faces(&parts, edge_ties, face_ties, budget.ties, &mut cancelled)?;
    let loads: Vec<Vec<f64>> = components
        .iter()
        .map(|part| {
            part.thermal
                .load_6n()
                .iter()
                .zip(part.dead_load_6n)
                .enumerate()
                .map(|(row, (&thermal, &dead))| {
                    if row % 64 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let sum = (Dd::from_f64(thermal) + Dd::from_f64(dead)).to_f64();
                    if !sum.is_finite() {
                        return Err(bad("assembled thermal/dead load overflow"));
                    }
                    Ok(sum)
                })
                .collect()
        })
        .collect::<Result<_, _>>()?;
    let rhs = assembly.project_loads(&loads.iter().map(Vec::as_slice).collect::<Vec<_>>())?;
    let coordinates = PreparedShellStatic::try_new(&assembly.model, &mut cancelled)?
        .solve_reduced(&rhs, &mut cancelled)?;

    // Check the actually returned f64 coordinates, not refinement's hidden low
    // words, in original units. Each retained coordinate is an original world
    // translation or rotation; eliminated slave rows have no new coordinate.
    let residual = original_residual(&assembly.model.k, &coordinates, &rhs, &mut cancelled)?;
    let mut residual_force_n = scaled_residual(&assembly, &residual, tolerance.moment_length_m)?;
    if residual_force_n > tolerance.maximum_residual_force_n {
        return Err(PlateError::StaticSolve {
            what: "assembled original force/moment residual exceeds tolerance",
        });
    }
    let displacements = assembly.expand_components(&coordinates)?;
    let component_free = assembly.expand_component_free(&coordinates)?;
    let mut recovered = Vec::with_capacity(components.len());
    let mut component_free_residuals = Vec::with_capacity(components.len());
    for (id, part) in components.iter().enumerate() {
        poll(&mut cancelled)?;
        let model = part.thermal.model();
        let force = model.project_load(&loads[id])?;
        component_free_residuals.push(original_residual(
            &model.pencil().k,
            &component_free[id],
            &force,
            &mut cancelled,
        )?);
        recovered.push(part.thermal.membrane_prestress_at(
            &displacements[id],
            limits,
            &mut cancelled,
        )?);
    }
    // Independently balance the original component forces and original glue
    // force through T^T. This also catches assembly/projection roundoff that an
    // already-congruent K*q residual alone could hide.
    let full_residuals: Vec<_> = components
        .iter()
        .zip(&component_free_residuals)
        .map(|(part, force)| part.thermal.model().expand(force))
        .collect::<Result<_, _>>()?;
    let mut joint_residual =
        assembly.project_loads(&full_residuals.iter().map(Vec::as_slice).collect::<Vec<_>>())?;
    let glue_force = original_residual(
        &assembly.glue_stiffness,
        &coordinates,
        &vec![0.0; coordinates.len()],
        &mut cancelled,
    )?;
    for (residual, glue) in joint_residual.iter_mut().zip(glue_force) {
        *residual = (Dd::from_f64(*residual) + Dd::from_f64(glue)).to_f64();
    }
    residual_force_n = residual_force_n.max(scaled_residual(
        &assembly,
        &joint_residual,
        tolerance.moment_length_m,
    )?);
    if residual_force_n > tolerance.maximum_residual_force_n {
        return Err(PlateError::StaticSolve {
            what: "assembled original component/interface force balance exceeds tolerance",
        });
    }
    let stiffnesses: Vec<_> = recovered.iter().map(|part| &part.pencil().k).collect();
    let tangent = assembly.assemble_component_stiffness(&stiffnesses, &mut cancelled)?;
    let elastic_stiffness = core::mem::replace(&mut assembly.model.k, tangent);
    poll(&mut cancelled)?;
    Ok(PrestressedShellBody {
        assembly,
        elastic_stiffness,
        components: recovered,
        coordinates,
        component_free_residuals,
        residual_force_n,
    })
}

fn bad(what: &'static str) -> PlateError {
    PlateError::BadConstraint { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() {
        Err(PlateError::Cancelled)
    } else {
        Ok(())
    }
}
fn cap(what: &'static str, requested: usize, limit: usize) -> Result<(), PlateError> {
    if requested > limit {
        Err(PlateError::Budget {
            what,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}

fn original_residual(
    matrix: &Csr,
    coordinates: &[f64],
    force: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, PlateError> {
    (0..matrix.nrows())
        .map(|row| {
            if row % 64 == 0 {
                poll(cancelled)?;
            }
            let (columns, values) = matrix.row(row);
            let value = columns
                .iter()
                .zip(values)
                .fold(-Dd::from_f64(force[row]), |sum, (&j, &a)| {
                    sum + Dd::from_f64(a) * Dd::from_f64(coordinates[j])
                })
                .to_f64();
            if !value.is_finite() {
                return Err(bad("assembled original residual overflow"));
            }
            Ok(value)
        })
        .collect()
}

fn scaled_residual(
    assembly: &TiedShellModel,
    residual: &[f64],
    length: f64,
) -> Result<f64, PlateError> {
    let mut scale = vec![0.0; assembly.model.free];
    for (physical, retained) in assembly.model.dof_map.iter().enumerate() {
        if let Some(index) = retained {
            scale[*index] = if physical % 6 < 3 { 1.0 } else { length };
        }
    }
    let mut maximum = 0.0_f64;
    for (&residual, scale) in residual.iter().zip(scale) {
        let force = residual.abs() / scale;
        if !force.is_finite() {
            return Err(bad("assembled force/moment norm is not representable"));
        }
        maximum = maximum.max(force);
    }
    Ok(maximum)
}

// The original junction owner eliminates free slave DOFs; it cannot propagate
// a prescribed slave support onto masters. Refuse overlapping conditions before
// they can be silently discarded. A seam-tangent rotation is conservatively
// admitted only when all slave rotations are free.
fn check_slave_supports(
    parts: &[ThermalShellComponent<'_>],
    edges: &[EdgeTie],
    faces: &[FaceTie],
) -> Result<(), PlateError> {
    let check = |component: usize, node: usize, rotation: RotationCoupling| {
        let part = parts
            .get(component)
            .ok_or_else(|| bad("unknown body tie component"))?;
        let map = &part.thermal.model().pencil().dof_map;
        let start = node
            .checked_mul(6)
            .ok_or_else(|| bad("body tie node overflow"))?;
        let dofs = if rotation == RotationCoupling::TranslationsOnly {
            3
        } else {
            6
        };
        let row = map
            .get(start..start.saturating_add(dofs))
            .ok_or_else(|| bad("unknown body tie node"))?;
        if row.iter().any(Option::is_none) {
            return Err(bad(
                "supported rigid-tie slaves require native constraint composition",
            ));
        }
        Ok(())
    };
    for tie in edges {
        if matches!(tie.kind, TieKind::Rigid) {
            for &node in &tie.b.nodes {
                check(tie.b.component, node, tie.rotations)?;
            }
        }
    }
    for tie in faces {
        if matches!(tie.kind, TieKind::Rigid) {
            let part = parts
                .get(tie.slave_component)
                .ok_or_else(|| bad("unknown body tie component"))?;
            for &triangle in &tie.slave_tris {
                let tri = part
                    .thermal
                    .model()
                    .original_geometry()
                    .0
                    .tris
                    .get(triangle)
                    .ok_or_else(|| bad("unknown body tie facet"))?;
                for &node in tri {
                    check(tie.slave_component, node, tie.rotations)?;
                }
            }
        }
    }
    Ok(())
}
