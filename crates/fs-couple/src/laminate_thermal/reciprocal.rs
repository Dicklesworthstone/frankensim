//! Reciprocal finite heat and compatible layered-shell equilibrium.
//!
//! The original caloric cells and native shell remain the two state owners.
//! At each layer, the SAME source plane-stress potential supplies Psi, entropy
//! S=-Psi_T, and U=Psi+T*S. Reversible heat is mechanical work minus delta U;
//! it is never an elastic-energy residual relabeled as dissipation. Mechanical
//! work uses endpoint-average stress times the actual compatible strain change.
//!
//! This lane admits temperature-INDEPENDENT elastic sources over each complete
//! declared interval. Instantaneous CTE may vary affinely and every thermal
//! cell retains its own temperature. Fixed elasticity keeps the native whole-
//! stack Whitney shear correction fixed: its physical work/storage and the
//! numerical drilling gauge are audited separately, without fictitious local
//! shear entropy. General Q(T) potentials exist in fs-material, but their
//! reciprocal use here requires the nonlocal shear-temperature derivatives.
//!
//! A bounded outer iteration reuses the original sparse implicit conduction
//! solve and one immutable native mechanical factorization. Every trial starts
//! from the accepted caloric history. Both owners advance only after local
//! heat, native material work, and combined energy gates pass. This is a
//! quasistatic endpoint scheme, generally first order in time; it is not exact
//! entropy-conserving integration, finite-rotation mechanics, phase evolution,
//! viscoelasticity, or deformation-dependent transport geometry.
//! The supplied caloric chart represents the same unloaded material path;
//! its source preparation's compatibility of population and work constraints
//! remains required. A common card does not by itself establish that fact.

use super::*;
use fs_material::elastic::thermal::ThermoelasticError;
use fs_material::laminate::thermal::thermoelastic::{
    AffineLaminaThermoelastic, ResolvedLaminaThermoelastic,
};
use fs_material::thermal::transport::ThermalTransportError;
use fs_plate::shell::mitc::strain::{ShellStrainPoint, integration_rule};

pub const MAX_RECIPROCAL_LAMINATE_ITERATIONS: usize = 256;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateThermoelasticControls {
    pub maximum_iterations: usize,
    /// Relaxation of the per-cell reversible joule input, in (0,1].
    pub heat_relaxation: f64,
    pub cell_heat_tolerance_j: f64,
    /// Independent material-work, sum-of-absolute-heat-defects and combined
    /// energy tolerance. This controls algebraic error, not timestep accuracy.
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum LaminateThermoelasticError {
    Invalid(&'static str),
    Source(ThermoelasticError),
    Shell(LaminateShellError),
    Native(PlateError),
    Thermal(DepositionError),
    HeatSolve(LayeredThermalSolveError),
    Capacity(ThermalTransportError),
    NotConverged {
        iterations: usize,
        maximum_heat_defect_j: f64,
        combined_residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for LaminateThermoelasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reciprocal layered thermoelasticity: {self:?}")
    }
}
impl std::error::Error for LaminateThermoelasticError {}

/// Integrated plane-stress elastic-field quantities for one ORIGINAL thermal
/// control volume, accumulated over all facets and its physical thickness.
/// Caloric storage, transverse shear and drilling storage are separate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LaminateCellPotential {
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    pub heat_capacity_j_k: f64,
}

pub struct ReciprocalLaminateState {
    shell: LaminateShellSnapshot,
    sources: Vec<ResolvedLaminaThermoelastic>,
    cells: Vec<LaminateCellPotential>,
    loads: Vec<f64>,
    shear_energy_j: f64,
    gauge_energy_j: f64,
    gauge_force: Vec<f64>,
}
impl ReciprocalLaminateState {
    pub const fn shell(&self) -> &LaminateShellSnapshot {
        &self.shell
    }
    pub fn cell_potentials(&self) -> &[LaminateCellPotential] {
        &self.cells
    }
    pub const fn transverse_shear_energy_j(&self) -> f64 {
        self.shear_energy_j
    }
    /// Numerical gauge storage; never deposited as physical heat.
    pub const fn drilling_gauge_energy_j(&self) -> f64 {
        self.gauge_energy_j
    }
    pub fn external_load_6n(&self) -> &[f64] {
        &self.loads
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateThermoelasticCellReport {
    pub mechanical_work_j: f64,
    pub thermal_strain_work_j: f64,
    pub internal_energy_change_j: f64,
    pub elastic_entropy_change_j_k: f64,
    pub reversible_heat_j: f64,
    pub applied_heat_j: f64,
    pub heat_defect_j: f64,
    pub free_energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LaminateThermoelasticReport {
    pub cells: Vec<LaminateThermoelasticCellReport>,
    /// The original thermal account includes the explicit reversible input.
    pub heat: LayeredThermalSolveReport,
    pub iterations: usize,
    /// Heat from authored sources and boundary reservoirs, excluding exchange
    /// with the elastic field. No thermal event is consumed inside the step.
    pub external_heat_j: f64,
    pub external_mechanical_work_j: f64,
    pub reversible_heat_j: f64,
    pub elastic_internal_energy_change_j: f64,
    pub transverse_shear_energy_change_j: f64,
    pub transverse_shear_work_j: f64,
    pub drilling_gauge_energy_change_j: f64,
    pub drilling_gauge_work_j: f64,
    pub mechanical_work_residual_j: f64,
    /// Caloric + elastic internal + shear + numerical gauge changes, minus
    /// external heat and external mechanical work. Gauge is explicitly named.
    pub combined_energy_residual_j: f64,
    pub maximum_heat_defect_j: f64,
    pub total_absolute_heat_defect_j: f64,
}

pub struct ReciprocalLaminateStep<'m> {
    pub thermal: DepositionSimulation<'m>,
    pub report: LaminateThermoelasticReport,
}

/// Retained native equilibrium, original source intervals, and fixed sparse
/// factorization. The caller retains its original caloric timeline; subsequent
/// steps verify its exact accepted temperatures/time before any solve.
pub struct ReciprocalLayeredThermalShell<'s, 'b, 'a> {
    source: &'s LayeredThermalShell<'b, 'a>,
    laws: Vec<AffineLaminaThermoelastic>,
    solver: PreparedMitcStatic,
    state: ReciprocalLaminateState,
    quadrature: Vec<([f64; 3], f64)>,
}

impl<'s, 'b, 'a> ReciprocalLayeredThermalShell<'s, 'b, 'a> {
    /// One interval per original thermal cell / mechanical sublamina. The
    /// declared principal engineering component bound covers both TOTAL and
    /// elastic strains; native Bernstein hulls admit the whole facet/layer.
    /// The original caloric capacity plus a conservative elastic contribution
    /// must be positive over that entire strain/temperature envelope.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        source: &'s LayeredThermalShell<'b, 'a>,
        thermal: &DepositionSimulation<'_>,
        external_load_6n: &[f64],
        temperature_ranges_k: &[[f64; 2]],
        maximum_absolute_engineering_strain: f64,
        controls: LaminateThermoelasticControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateThermoelasticError> {
        check(&mut cancelled)?;
        admit_controls(controls)?;
        if temperature_ranges_k.len() != source.source.layers().len() {
            return Err(invalid(
                "one complete source interval per thermal cell required",
            ));
        }
        let body = source.thermal.prepared_body();
        let area = body.conduction().area_m2();
        let mut laws = Vec::with_capacity(temperature_ranges_k.len());
        let mut cell = 0;
        for (layer, &count) in body.layer_cell_counts().iter().enumerate() {
            let storage = &body.storages()[layer];
            let span = storage.temperature_span_k();
            for _ in 0..count {
                check(&mut cancelled)?;
                let range = temperature_ranges_k[cell];
                let original = &source.source.layers()[cell];
                let law = AffineLaminaThermoelastic::new(
                    original,
                    range,
                    maximum_absolute_engineering_strain,
                )
                .map_err(LaminateThermoelasticError::Source)?;
                if !law.is_temperature_independent_elasticity() {
                    return Err(invalid(
                        "reciprocal native laminate requires fixed elastic sources; variable Whitney shear entropy is not yet owned",
                    ));
                }
                if range[0] < span[0] || range[1] > span[1] {
                    return Err(invalid(
                        "original caloric chart must cover the complete elastic source interval",
                    ));
                }
                // The native caloric capacity is piecewise affine; endpoints
                // and its ORIGINAL breakpoints attain its minimum.
                let mut minimum = f64::INFINITY;
                for t in range.into_iter().chain(
                    storage
                        .temperatures_k()
                        .filter(|t| *t > range[0] && *t < range[1]),
                ) {
                    check(&mut cancelled)?;
                    minimum = minimum.min(
                        storage
                            .heat_capacity_j_k(t)
                            .map_err(LaminateThermoelasticError::Capacity)?,
                    );
                }
                let volume = area * original.input().thickness_m;
                let elastic_minimum = law
                    .minimum_heat_capacity_j_k(volume)
                    .map_err(LaminateThermoelasticError::Source)?;
                if finite(minimum + elastic_minimum)? <= 0.0 {
                    return Err(invalid(
                        "caloric plus elastic capacity is not positive over the complete admitted envelope",
                    ));
                }
                laws.push(law);
                cell += 1;
            }
        }
        let shell = source
            .resolve(thermal, external_load_6n, &mut cancelled)
            .map_err(LaminateThermoelasticError::Shell)?;
        let solver = PreparedMitcStatic::try_new(Arc::clone(&shell.model), &mut cancelled)
            .map_err(LaminateThermoelasticError::Native)?;
        let quadrature = integration_rule();
        let proposal = assemble_state(
            source,
            &laws,
            &quadrature,
            shell,
            external_load_6n,
            None,
            controls.energy_tolerance_j,
            &mut cancelled,
        )?;
        check(&mut cancelled)?;
        Ok(Self {
            source,
            laws,
            solver,
            state: proposal.state,
            quadrature,
        })
    }

    pub const fn state(&self) -> &ReciprocalLaminateState {
        &self.state
    }
    pub fn time_s(&self) -> f64 {
        self.state.shell.time_s()
    }

    /// Advance the ORIGINAL thermal timeline and native equilibrium atomically.
    /// Signed reversible joules are frozen only inside each thermal Newton
    /// solve, then updated by the outer fixed point. Nonconvergence is explicit;
    /// accepted physical history remains unchanged on every failed proposal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'m>(
        &mut self,
        initial: &DepositionSimulation<'m>,
        target_s: f64,
        external_load_6n: &[f64],
        conduction: LayerConductionConfig,
        heat_controls: LayeredThermalSolveControls,
        controls: LaminateThermoelasticControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReciprocalLaminateStep<'m>, LaminateThermoelasticError> {
        check(&mut cancelled)?;
        admit_controls(controls)?;
        let initial_t = self
            .source
            .thermal
            .cell_temperatures(initial, || thermal_poll(&mut cancelled))
            .map_err(LaminateThermoelasticError::Thermal)?;
        if initial.time_s().to_bits() != self.time_s().to_bits()
            || initial_t
                .iter()
                .zip(self.state.shell.temperatures_k())
                .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err(invalid(
                "accepted shell and original thermal history must remain exactly paired",
            ));
        }
        let mut applied = vec![0.0; self.laws.len()];
        let mut last = (f64::INFINITY, f64::INFINITY);
        for iteration in 0..controls.maximum_iterations {
            check(&mut cancelled)?;
            let step =
                LayeredThermalStep::new(self.source.thermal, initial, target_s, conduction, || {
                    thermal_poll(&mut cancelled)
                })
                .map_err(LaminateThermoelasticError::Thermal)?
                .with_cell_energy_j(&applied, || thermal_poll(&mut cancelled))
                .map_err(LaminateThermoelasticError::Thermal)?;
            let (thermal, _, heat) = step
                .solve(heat_controls, workspace, || thermal_poll(&mut cancelled))
                .map_err(LaminateThermoelasticError::HeatSolve)?;
            let shell = self
                .source
                .resolve_with_solver(
                    &thermal,
                    external_load_6n,
                    Some(&self.solver),
                    &mut cancelled,
                )
                .map_err(LaminateThermoelasticError::Shell)?;
            let proposal = assemble_state(
                self.source,
                &self.laws,
                &self.quadrature,
                shell,
                external_load_6n,
                Some(&self.state),
                controls.energy_tolerance_j,
                &mut cancelled,
            )?;
            let report = account(
                &self.state,
                &proposal,
                &applied,
                heat,
                iteration + 1,
                controls.energy_tolerance_j,
                &mut cancelled,
            )?;
            last = (
                report.maximum_heat_defect_j,
                report.combined_energy_residual_j,
            );
            if report.maximum_heat_defect_j <= controls.cell_heat_tolerance_j
                && report.total_absolute_heat_defect_j <= controls.energy_tolerance_j
                && report.combined_energy_residual_j.abs() <= controls.energy_tolerance_j
            {
                check(&mut cancelled)?;
                self.state = proposal.state;
                return Ok(ReciprocalLaminateStep { thermal, report });
            }
            for (input, cell) in applied.iter_mut().zip(&report.cells) {
                check(&mut cancelled)?;
                *input = finite(*input + controls.heat_relaxation * cell.heat_defect_j)?;
            }
        }
        Err(LaminateThermoelasticError::NotConverged {
            iterations: controls.maximum_iterations,
            maximum_heat_defect_j: last.0,
            combined_residual_j: last.1,
        })
    }
}

struct Proposal {
    state: ReciprocalLaminateState,
    mechanical_work: Vec<f64>,
    expansion_work: Vec<f64>,
    shear_work: f64,
}

#[allow(clippy::too_many_arguments)]
fn assemble_state(
    source: &LayeredThermalShell<'_, '_>,
    laws: &[AffineLaminaThermoelastic],
    quadrature: &[([f64; 3], f64)],
    shell: LaminateShellSnapshot,
    loads: &[f64],
    previous: Option<&ReciprocalLaminateState>,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Proposal, LaminateThermoelasticError> {
    let resolved = laws
        .iter()
        .zip(shell.material.layers())
        .map(|(law, material)| law.at_materials(material))
        .collect::<Result<Vec<_>, _>>()
        .map_err(LaminateThermoelasticError::Source)?;
    let mut cells = vec![LaminateCellPotential::default(); laws.len()];
    let mut mechanical_work = vec![0.0; laws.len()];
    let mut expansion_work = vec![0.0; laws.len()];
    let mut shear_energy = 0.0;
    let mut shear_work = 0.0;
    let section = shell.material.prepared().section();
    let shear = section.shear_pa_m();
    if let Some(old) = previous {
        if old.shell.material.prepared().section().shear_pa_m() != shear {
            return Err(invalid(
                "fixed reciprocal shear section changed unexpectedly",
            ));
        }
    }
    let zero = ShellStrainPoint {
        membrane: [0.0; 3],
        gradient_per_m: [0.0; 3],
        transverse_shear: [0.0; 2],
    };
    for (element, field) in shell.strain_fields.iter().enumerate() {
        check(cancelled)?;
        let area = source
            .discretization
            .mesh
            .facet(element)
            .map_err(LaminateThermoelasticError::Native)?
            .area_m2;
        let section_angle = source.discretization.section_angles_rad[element];
        let mut points = [zero; 7];
        let mut old_points = [zero; 7];
        for (i, &(lambda, weight)) in quadrature.iter().enumerate() {
            check(cancelled)?;
            points[i] = field
                .at_barycentric(lambda)
                .map_err(LaminateThermoelasticError::Native)?;
            let gamma = rotate_shear(points[i].transverse_shear, section_angle);
            let stress = shear_stress(shear, gamma);
            shear_energy = finite(shear_energy + 0.5 * area * weight * dot2(gamma, stress))?;
            if let Some(old) = previous {
                old_points[i] = old.shell.strain_fields[element]
                    .at_barycentric(lambda)
                    .map_err(LaminateThermoelasticError::Native)?;
                let old_gamma = rotate_shear(old_points[i].transverse_shear, section_angle);
                let old_stress = shear_stress(shear, old_gamma);
                shear_work = finite(
                    shear_work
                        + 0.5
                            * area
                            * weight
                            * (0..2)
                                .map(|k| (stress[k] + old_stress[k]) * (gamma[k] - old_gamma[k]))
                                .sum::<f64>(),
                )?;
            }
        }
        for (layer, geometry) in section.layer_integrals().iter().enumerate() {
            check(cancelled)?;
            let angle = section_angle + laws[layer].source().input().orientation_rad;
            let rotation = [det::cos(angle), det::sin(angle)];
            let faces = [geometry.z_bottom_m, geometry.z_top_m];
            let free = resolved[layer].free_strain();
            let limit = laws[layer].maximum_absolute_engineering_strain();
            // Bernstein convexity and affine thickness cover the ENTIRE layer,
            // including shear and internal director motion, before quadrature.
            for control in field.bernstein_controls() {
                for z in faces {
                    let strain = principal(control, z, rotation);
                    let gamma = rotate_shear_with(control.transverse_shear, rotation);
                    if strain
                        .iter()
                        .enumerate()
                        .any(|(k, v)| v.abs() > limit || (*v - free[k]).abs() > limit)
                        || gamma.iter().any(|g| g.abs() > limit)
                    {
                        return Err(invalid(
                            "whole native layer leaves the admitted total/elastic principal strain envelope",
                        ));
                    }
                }
            }
            let half = 0.5 * (faces[1] - faces[0]);
            let middle = faces[0] + half;
            let offset = half / det::sqrt(3.0);
            for (i, &(_, weight)) in quadrature.iter().enumerate() {
                for z in [middle - offset, middle + offset] {
                    check(cancelled)?;
                    let volume = area * weight * half;
                    let strain = principal(points[i], z, rotation);
                    let at = resolved[layer]
                        .in_plane_state(strain, volume)
                        .map_err(LaminateThermoelasticError::Source)?;
                    let cell = &mut cells[layer];
                    cell.free_energy_j = finite(cell.free_energy_j + at.free_energy_j)?;
                    cell.entropy_j_k = finite(cell.entropy_j_k + at.entropy_j_k)?;
                    cell.internal_energy_j = finite(cell.internal_energy_j + at.internal_energy_j)?;
                    cell.heat_capacity_j_k = finite(cell.heat_capacity_j_k + at.heat_capacity_j_k)?;
                    if let Some(old) = previous {
                        let old_strain = principal(old_points[i], z, rotation);
                        let old_at = old.sources[layer]
                            .in_plane_state(old_strain, volume)
                            .map_err(LaminateThermoelasticError::Source)?;
                        let old_free = old.sources[layer].free_strain();
                        let work = 0.5
                            * volume
                            * (0..3)
                                .map(|k| {
                                    (at.stress_pa[k] + old_at.stress_pa[k])
                                        * (strain[k] - old_strain[k])
                                })
                                .sum::<f64>();
                        let expansion = -0.5
                            * volume
                            * (0..3)
                                .map(|k| {
                                    (at.stress_pa[k] + old_at.stress_pa[k])
                                        * (free[k] - old_free[k])
                                })
                                .sum::<f64>();
                        mechanical_work[layer] = finite(mechanical_work[layer] + work)?;
                        expansion_work[layer] = finite(expansion_work[layer] + expansion)?;
                    }
                }
            }
        }
    }
    gate(
        cells.iter().map(|c| c.free_energy_j).sum::<f64>() + shear_energy
            - shell.response.strain_energy_j,
        tolerance,
        "source plane potential and original native material storage disagree",
    )?;
    let (gauge_energy, gauge_force) = drilling_storage(&shell, tolerance, cancelled)?;
    check(cancelled)?;
    Ok(Proposal {
        state: ReciprocalLaminateState {
            shell,
            sources: resolved,
            cells,
            loads: loads.to_vec(),
            shear_energy_j: shear_energy,
            gauge_energy_j: gauge_energy,
            gauge_force,
        },
        mechanical_work,
        expansion_work,
        shear_work,
    })
}

fn drilling_storage(
    shell: &LaminateShellSnapshot,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, Vec<f64>), LaminateThermoelasticError> {
    let model = shell.model();
    let q = model
        .project_load(&shell.response.displacement_6n)
        .map_err(LaminateThermoelasticError::Native)?;
    let material = model
        .project_load(&shell.response.material_nodal_force_6n)
        .map_err(LaminateThermoelasticError::Native)?;
    let thermal = model
        .project_load(&shell.thermal_load_6n)
        .map_err(LaminateThermoelasticError::Native)?;
    let mut force = Vec::with_capacity(q.len());
    let mut energy = 0.0;
    for (row, &value) in q.iter().enumerate() {
        check(cancelled)?;
        let (cols, vals) = model.pencil().k.row(row);
        let total = cols
            .iter()
            .zip(vals)
            .map(|(&col, &a)| a * q[col])
            .sum::<f64>();
        // This subtraction isolates the ORIGINAL numerical drilling gauge,
        // independently of equilibrium loads. It is not a heat correction.
        let f = finite(total - thermal[row] - material[row])?;
        energy = finite(energy + 0.5 * value * f)?;
        force.push(f);
    }
    if energy < -tolerance {
        return Err(invalid("native numerical drilling storage became negative"));
    }
    Ok((energy, force))
}

fn account(
    old: &ReciprocalLaminateState,
    proposal: &Proposal,
    applied: &[f64],
    heat: LayeredThermalSolveReport,
    iterations: usize,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<LaminateThermoelasticReport, LaminateThermoelasticError> {
    let next = &proposal.state;
    let mut cells = Vec::with_capacity(applied.len());
    let mut heat_defect_sum = 0.0;
    let mut heat_defect_abs = 0.0;
    let mut maximum = 0.0_f64;
    let mut reversible = 0.0;
    let mut internal_change = 0.0;
    for i in 0..applied.len() {
        check(cancelled)?;
        let a = &old.cells[i];
        let b = &next.cells[i];
        let mechanical = proposal.mechanical_work[i];
        let expansion = proposal.expansion_work[i];
        let residual = finite(b.free_energy_j - a.free_energy_j - mechanical - expansion)?;
        gate(
            residual,
            tolerance,
            "native layer work does not close its original free energy",
        )?;
        let du = finite(b.internal_energy_j - a.internal_energy_j)?;
        let desired = finite(mechanical - du)?;
        let defect = finite(desired - applied[i])?;
        heat_defect_sum = finite(heat_defect_sum + defect)?;
        heat_defect_abs = finite(heat_defect_abs + defect.abs())?;
        maximum = maximum.max(defect.abs());
        reversible = finite(reversible + desired)?;
        internal_change = finite(internal_change + du)?;
        cells.push(LaminateThermoelasticCellReport {
            mechanical_work_j: mechanical,
            thermal_strain_work_j: expansion,
            internal_energy_change_j: du,
            elastic_entropy_change_j_k: b.entropy_j_k - a.entropy_j_k,
            reversible_heat_j: desired,
            applied_heat_j: applied[i],
            heat_defect_j: defect,
            free_energy_residual_j: residual,
        });
    }
    let u0 = &old.shell.response.displacement_6n;
    let u1 = &next.shell.response.displacement_6n;
    let external_work = finite(
        (0..u0.len())
            .map(|i| 0.5 * (old.loads[i] + next.loads[i]) * (u1[i] - u0[i]))
            .sum(),
    )?;
    let material_work = finite(
        (0..u0.len())
            .map(|i| {
                0.5 * (old.shell.response.material_nodal_force_6n[i]
                    + next.shell.response.material_nodal_force_6n[i])
                    * (u1[i] - u0[i])
            })
            .sum(),
    )?;
    let membrane_work = proposal.mechanical_work.iter().sum::<f64>();
    gate(
        material_work - membrane_work - proposal.shear_work,
        tolerance,
        "native nodal material work and integrated layer/shear work disagree",
    )?;
    let q0 = old
        .shell
        .model()
        .project_load(u0)
        .map_err(LaminateThermoelasticError::Native)?;
    let q1 = next
        .shell
        .model()
        .project_load(u1)
        .map_err(LaminateThermoelasticError::Native)?;
    let gauge_work = finite(
        (0..q0.len())
            .map(|i| 0.5 * (old.gauge_force[i] + next.gauge_force[i]) * (q1[i] - q0[i]))
            .sum(),
    )?;
    let shear_change = finite(next.shear_energy_j - old.shear_energy_j)?;
    let gauge_change = finite(next.gauge_energy_j - old.gauge_energy_j)?;
    gate(
        shear_change - proposal.shear_work,
        tolerance,
        "fixed native shear work/storage mismatch",
    )?;
    gate(
        gauge_change - gauge_work,
        tolerance,
        "fixed numerical gauge work/storage mismatch",
    )?;
    let mechanical_residual =
        finite(membrane_work + proposal.shear_work + gauge_work - external_work)?;
    gate(
        mechanical_residual,
        tolerance,
        "native equilibrium does not close the external work port",
    )?;
    // Use the original compensated caloric residual, avoiding subtraction of
    // large enthalpy datums. Desired heat is not substituted for applied heat.
    let combined = finite(
        heat.energy.residual_j - heat_defect_sum + mechanical_residual + shear_change
            - proposal.shear_work
            + gauge_change
            - gauge_work,
    )?;
    Ok(LaminateThermoelasticReport {
        cells,
        external_heat_j: finite(heat.energy.external_heat_j - applied.iter().sum::<f64>())?,
        heat,
        iterations,
        external_mechanical_work_j: external_work,
        reversible_heat_j: reversible,
        elastic_internal_energy_change_j: internal_change,
        transverse_shear_energy_change_j: shear_change,
        transverse_shear_work_j: proposal.shear_work,
        drilling_gauge_energy_change_j: gauge_change,
        drilling_gauge_work_j: gauge_work,
        mechanical_work_residual_j: mechanical_residual,
        combined_energy_residual_j: combined,
        maximum_heat_defect_j: maximum,
        total_absolute_heat_defect_j: heat_defect_abs,
    })
}

fn principal(point: ShellStrainPoint, z: f64, [c, s]: [f64; 2]) -> [f64; 3] {
    let e: [f64; 3] = core::array::from_fn(|i| point.membrane[i] + z * point.gradient_per_m[i]);
    [
        c * c * e[0] + s * s * e[1] + c * s * e[2],
        s * s * e[0] + c * c * e[1] - c * s * e[2],
        -2.0 * c * s * e[0] + 2.0 * c * s * e[1] + (c * c - s * s) * e[2],
    ]
}
fn rotate_shear(gamma: [f64; 2], angle: f64) -> [f64; 2] {
    rotate_shear_with(gamma, [det::cos(angle), det::sin(angle)])
}
fn rotate_shear_with([x, y]: [f64; 2], [c, s]: [f64; 2]) -> [f64; 2] {
    [c * x + s * y, -s * x + c * y]
}
fn shear_stress(s: &[f64; 4], g: [f64; 2]) -> [f64; 2] {
    [s[0] * g[0] + s[1] * g[1], s[2] * g[0] + s[3] * g[1]]
}
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), LaminateThermoelasticError> {
    if cancelled() {
        Err(LaminateThermoelasticError::Cancelled)
    } else {
        Ok(())
    }
}
fn invalid(what: &'static str) -> LaminateThermoelasticError {
    LaminateThermoelasticError::Invalid(what)
}
fn finite(x: f64) -> Result<f64, LaminateThermoelasticError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("coupled work/heat is unrepresentable"))
    }
}
fn gate(x: f64, tolerance: f64, what: &'static str) -> Result<(), LaminateThermoelasticError> {
    if !x.is_finite() || x.abs() > tolerance {
        Err(invalid(what))
    } else {
        Ok(())
    }
}
fn admit_controls(c: LaminateThermoelasticControls) -> Result<(), LaminateThermoelasticError> {
    if c.maximum_iterations == 0
        || c.maximum_iterations > MAX_RECIPROCAL_LAMINATE_ITERATIONS
        || !c.heat_relaxation.is_finite()
        || c.heat_relaxation <= 0.0
        || c.heat_relaxation > 1.0
        || [c.cell_heat_tolerance_j, c.energy_tolerance_j]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
    {
        Err(invalid(
            "bounded iterations, relaxation and positive finite joule tolerances required",
        ))
    } else {
        Ok(())
    }
}
