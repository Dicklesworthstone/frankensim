//! Native eight-point Hex8 plasticity with reciprocal finite-cell heat feedback.
//!
//! Every original Gauss point binds one distinct active caloric cell with the
//! native reference integration volume and sourced reference density. The
//! existing Standard/B-bar equilibrium or consistent-mass backward-Euler owner
//! advances complete native histories at the actual interval. Physical plus
//! reversible material heat enters the original DepositionSimulation once;
//! free-energy storage, mechanical entropy and numerical losses stay separate.
//!
//! The existing Picard/IQN-ILS driver solves the temperature/history feedback.
//! All trials restart from the same accepted displacement, histories, velocity,
//! exact specific enthalpy and observer. Acceptance returns these owners together.
//! This is a fixed-reference, fully solid, small-strain/rotation model. Point
//! cells represent a supplied thermal network, not a geometric dual mesh;
//! transport geometry and its refinement remain the thermal owner's responsibility.

/// Original indentation hardness admitted with the coupled native candidate.
pub mod hardness;

use std::collections::BTreeSet;
use std::sync::Arc;

use crate::iqn_ils::driver::CouplingMethod;
use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::plastic::section::thermal::source::{
    ThermalSectionSourceError, ThermalSectionSourceFamily,
};
use fs_solid::Formulation;
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_solid::plastic3::hex::thermal::{
    ThermalHexEnergy, ThermalHexError, ThermalHexHistory, ThermalHexMaterialEnergy,
    ThermalHexProblem, ThermalHexTolerances,
};
use fs_solid::plastic3::hex::{HexAssemblyBudget, SmallStrainHexProblem, SmallStrainHexSolution};
use fs_solid::plastic3::{IncrementRecord, IncrementSettings, SmallStrainLawField};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionObserver, DepositionSimulation, StepConfig, ThermalLaw,
};

use super::finite_plastic::FiniteThermoplasticReport;
use super::finite_plastic::feedback::{
    self, CaloricMaterial, MechanicalEnergy, MechanicalProposal, ThermalBinding,
};
pub use super::finite_plastic::{
    FiniteThermoplasticConfig as HexThermoplasticConfig,
    FiniteThermoplasticError as HexThermoplasticError,
};

type History<'a, F> = <<F as ThermalSectionSourceFamily<'a>>::Material as SmallStrainLaw>::State;

/// Original mechanical reduction and explicit native resource limits.
#[derive(Clone, Copy, Debug)]
pub struct HexThermoplasticDiscretization {
    pub formulation: Formulation,
    pub strain_limit: f64,
    pub rotation_limit_rad: f64,
    pub budget: HexAssemblyBudget,
}

/// Accepted state and immutable reference geometry/source families.
/// Cell order is element-major, then the native eight-point Gauss order.
#[derive(Clone, Debug)]
pub struct ThermoplasticHexSolid<'a, F: ThermalSectionSourceFamily<'a>> {
    nodes_m: &'a [[f64; 3]],
    hexahedra: &'a [[usize; 8]],
    families: Arc<[F]>,
    cells: Vec<DepositedCell>,
    solution: SmallStrainHexSolution<History<'a, F>>,
    velocity_mps: Option<Vec<[f64; 3]>>,
    discretization: HexThermoplasticDiscretization,
    time_s: f64,
}

struct PointCaloric<'a>(&'a ResolvedHeatCapacityCurve);
impl CaloricMaterial for PointCaloric<'_> {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        self.0
    }
}

// The shared driver owns caloric storage. The native free-energy and T*s
// changes form the remaining mechanical internal energy; its numerical loss
// is an endpoint-work account, never a second physical heat source.
impl MechanicalEnergy for ThermalHexEnergy {
    fn internal_energy_change_j(&self) -> f64 {
        self.material.stored_change_j + self.material.thermoelastic_storage_change_j
    }
    fn endpoint_work_defect_j(&self) -> f64 {
        self.material.numerical_loss_j
    }
    fn dead_load_endpoint_work_j(&self) -> f64 {
        self.applied_work_j
    }
    fn support_endpoint_work_j(&self) -> f64 {
        self.support_work_j
    }
    fn nonthermal_loss_j(&self) -> f64 {
        0.0
    }
}

impl<'a, F: ThermalSectionSourceFamily<'a>> ThermoplasticHexSolid<'a, F> {
    /// Bind virgin stress-free point histories to their original caloric cells.
    /// One immutable family belongs to each element; the eight cells need not
    /// have equal volumes on a warped brick. No material increment is taken.
    ///
    /// Source resolution at unit duration is used only to obtain its native
    /// initial history and inspect zero stress. It consumes no physical time.
    /// Prestrained or thermally stressed initial states need an explicit native
    /// initialization; they are not silently declared stress free here.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        hexahedra: &'a [[usize; 8]],
        families: Vec<F>,
        point_cells: &[CellId],
        discretization: HexThermoplasticDiscretization,
        cx: &Cx<'_>,
    ) -> Result<Self, HexThermoplasticError> {
        poll(cx)?;
        let points = hexahedra
            .len()
            .checked_mul(8)
            .ok_or(HexThermoplasticError::Invalid("Hex8 point count overflow"))?;
        let dofs = nodes_m
            .len()
            .checked_mul(3)
            .ok_or(HexThermoplasticError::Invalid("Hex8 DOF count overflow"))?;
        if points == 0
            || families.len() != hexahedra.len()
            || point_cells.len() != points
            || nodes_m.len() > discretization.budget.maximum_nodes
            || hexahedra.len() > discretization.budget.maximum_elements
            || points > discretization.budget.maximum_quadrature_points
        {
            return Err(HexThermoplasticError::Invalid(
                "bounded nonempty Hex8 mesh, one family per element and eight cells per element required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(points);
        let mut laws = Vec::with_capacity(points);
        let mut states = Vec::with_capacity(points);
        for (point, &id) in point_cells.iter().enumerate() {
            poll(cx)?;
            let source = &families[point / 8];
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(HexThermoplasticError::Invalid("inactive Hex8 caloric cell"))?;
            if !used.insert(id)
                || !thermal.uses_material_law(
                    cell.material_id,
                    source.caloric_source().curve() as &dyn ThermalLaw,
                )
            {
                return Err(HexThermoplasticError::Invalid(
                    "unique Hex8 point cells and original caloric objects required",
                ));
            }
            source
                .admit_section_enthalpy(cell.specific_h_j_kg)
                .map_err(|error| source_error(point, error))?;
            let resolved = source
                .resolve_section_step(cell.specific_h_j_kg, cell.specific_h_j_kg, 1.0)
                .map_err(|error| source_error(point, error))?;
            let state = resolved.law.initial_state();
            if resolved
                .law
                .stress(&[0.0; 6], &state)
                .iter()
                .any(|&value| value != 0.0)
            {
                return Err(HexThermoplasticError::Invalid(
                    "virgin Hex8 sources must be stress free at actual initial enthalpy",
                ));
            }
            states.push(state);
            laws.push(resolved.law);
            cells.push(cell.clone());
        }
        let fields: Vec<_> = laws.chunks_exact(8).collect();
        // Geometry admission is independent of a future free-coordinate budget.
        // It evaluates no constitutive increment, stress or tangent.
        let reference_supports: Vec<_> = (0..dofs).map(|dof| (dof, 0.0)).collect();
        let reference = SmallStrainHexProblem {
            nodes_m,
            hexahedra,
            laws: SmallStrainLawField::PerQuadrature(&fields),
            formulation: discretization.formulation,
            strain_limit: discretization.strain_limit,
            rotation_limit_rad: discretization.rotation_limit_rad,
            prescribed_m: &reference_supports,
            nodal_forces_n: &[],
            load_path: &[1.0],
            budget: discretization.budget,
            settings: IncrementSettings::default(),
        };
        let volumes = reference
            .reference_quadrature_volumes(cx)
            .map_err(|error| HexThermoplasticError::Hex(ThermalHexError::Mechanical(error)))?;
        for (point, (cell, &volume)) in cells.iter().zip(&volumes).enumerate() {
            poll(cx)?;
            if !super::same_measure(cell.geometry.volume_m3, volume)
                || !super::same_measure(
                    cell.mass_kg,
                    families[point / 8].reference_density_kg_m3() * volume,
                )
            {
                return Err(HexThermoplasticError::Invalid(
                    "each caloric cell must match its native quadrature volume and sourced reference mass",
                ));
            }
        }
        let solution = SmallStrainHexSolution {
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            quadrature_strain: vec![[0.0; 6]; points],
            quadrature_stress_pa: vec![[0.0; 6]; points],
            quadrature_states: states,
            quadrature_volume_m3: volumes,
            element_offsets: (0..=hexahedra.len()).map(|element| 8 * element).collect(),
            increments: Vec::new(),
        };
        poll(cx)?;
        Ok(Self {
            nodes_m,
            hexahedra,
            families: families.into(),
            cells,
            solution,
            velocity_mps: None,
            discretization,
            time_s: thermal.time_s(),
        })
    }

    /// Select consistent-mass motion before the first accepted material interval.
    /// This supplies original physical momentum, not an implicit heat source.
    pub fn with_initial_velocity(
        mut self,
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, HexThermoplasticError> {
        poll(cx)?;
        if !self.solution.increments.is_empty()
            || self.velocity_mps.is_some()
            || velocity_mps.len() != self.nodes_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(HexThermoplasticError::Invalid(
                "finite original nodal velocity may be bound once before Hex8 evolution",
            ));
        }
        self.velocity_mps = Some(velocity_mps.to_vec());
        poll(cx)?;
        Ok(self)
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn solution(&self) -> &SmallStrainHexSolution<History<'a, F>> {
        &self.solution
    }
    pub fn cells(&self) -> &[DepositedCell] {
        &self.cells
    }
    pub fn families(&self) -> &[F] {
        &self.families
    }
    pub fn velocity_mps(&self) -> Option<&[[f64; 3]]> {
        self.velocity_mps.as_deref()
    }

    /// Resolve conductivity on these same source cards and exact point cells.
    /// The existing material-network observer owns transport and thermal substeps.
    pub fn thermal_material_cells<'s>(
        &'s self,
        conductivity: &'s [Option<(
            fs_matdb::PropertyKey,
            fs_material::state_point::MaterialPropertySelection,
        )>],
    ) -> Result<Vec<crate::circuit_enthalpy::material::MaterialCell<'s>>, HexThermoplasticError>
    {
        if conductivity.len() != self.cells.len() {
            return Err(HexThermoplasticError::Invalid(
                "one transport selection per Hex8 point required",
            ));
        }
        self.cells
            .iter()
            .zip(conductivity)
            .enumerate()
            .map(|(point, (cell, selection))| {
                let source = &self.families[point / 8];
                let conductivity = selection
                    .as_ref()
                    .map(|(key, policy)| {
                        source
                            .pin_thermal_transport(cell.specific_h_j_kg, key, policy.clone())
                            .map(|pinned| (key, pinned))
                            .map_err(|error| source_error(point, error))
                    })
                    .transpose()?;
                Ok(crate::circuit_enthalpy::material::MaterialCell {
                    cell: cell.id,
                    material_id: cell.material_id,
                    curve: source.caloric_source().curve(),
                    point: source.thermal_material_point(),
                    conductivity,
                })
            })
            .collect()
    }

    /// Solve one common mechanical/thermal interval. Initial velocity selects
    /// backward-Euler motion; otherwise this is quasistatic equilibrium.
    /// Picard and IQN-ILS share the existing physical, thermal and energy gates.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        tolerances: ThermalHexTolerances,
        thermal_config: StepConfig,
        config: HexThermoplasticConfig,
        method: CouplingMethod,
        cx: &Cx<'_>,
    ) -> Result<HexThermoplasticCandidate<'a, 't, F, O>, HexThermoplasticError<O::Error>> {
        poll(cx)?;
        let calorics: Vec<_> = (0..self.cells.len())
            .map(|point| PointCaloric(self.families[point / 8].caloric_source()))
            .collect();
        let old_h: Vec<_> = self.cells.iter().map(|cell| cell.specific_h_j_kg).collect();
        let binding = ThermalBinding {
            cells: &self.cells,
            materials: FinitePlasticLawField::PerElement(&calorics),
        };
        let d = self.discretization;
        let problem = ThermalHexProblem {
            nodes_m: self.nodes_m,
            hexahedra: self.hexahedra,
            families: &self.families,
            formulation: d.formulation,
            strain_limit: d.strain_limit,
            rotation_limit_rad: d.rotation_limit_rad,
            prescribed_m,
            nodal_forces_n,
            budget: d.budget,
            settings,
        };
        let candidate = feedback::advance_with_enthalpy(
            binding,
            self.time_s,
            thermal,
            observer,
            target_s,
            thermal_config,
            config,
            method,
            cx,
            |enthalpies, duration_s| {
                let history = ThermalHexHistory {
                    displacement_m: self.solution.displacement_m.as_flattened(),
                    quadrature_states: &self.solution.quadrature_states,
                    specific_enthalpy_j_kg: &old_h,
                };
                let step = if let Some(velocity) = &self.velocity_mps {
                    problem.advance_backward_euler(
                        history, velocity, enthalpies, duration_s, tolerances, cx,
                    )
                } else {
                    problem.solve_from(history, enthalpies, duration_s, tolerances, cx)
                }
                .map_err(HexThermoplasticError::Hex)?;
                if step.solution.quadrature_volume_m3 != self.solution.quadrature_volume_m3 {
                    return Err(HexThermoplasticError::Invalid(
                        "native Hex8 quadrature inventory changed during continuation",
                    ));
                }
                if let Some(inertia) = step.inertia
                    && !super::same_measure(
                        inertia.reference_mass_kg,
                        self.cells.iter().map(|cell| cell.mass_kg).sum(),
                    )
                {
                    return Err(HexThermoplasticError::Invalid(
                        "native consistent mass disagrees with original caloric inventory",
                    ));
                }
                let heat_j = step
                    .quadrature_energy
                    .iter()
                    .map(|point| finite(point.physical_heat_j + point.reversible_heat_j))
                    .collect::<Result<_, HexThermoplasticError<O::Error>>>()?;
                Ok(MechanicalProposal {
                    equilibrium: step.solution.increments[0].clone(),
                    energy: step.energy,
                    inertia: step.inertia,
                    heat_j,
                    state: step,
                })
            },
        )?;
        if candidate
            .state
            .specific_enthalpy_j_kg
            .iter()
            .zip(&candidate.cells)
            .any(|(h, cell)| h.to_bits() != cell.specific_h_j_kg.to_bits())
        {
            return Err(HexThermoplasticError::Invalid(
                "native Hex8 proposal must retain the exact accepted caloric coordinates",
            ));
        }
        let solid = Self {
            nodes_m: self.nodes_m,
            hexahedra: self.hexahedra,
            families: Arc::clone(&self.families),
            cells: candidate.cells,
            solution: candidate.state.solution,
            velocity_mps: candidate.state.velocity_mps,
            discretization: self.discretization,
            time_s: target_s,
        };
        poll(cx)?;
        Ok(HexThermoplasticCandidate {
            thermal: candidate.thermal,
            solid,
            observer: candidate.observer,
            report: candidate.report,
            quadrature_energy: candidate.state.quadrature_energy,
        })
    }
}

/// The existing coupled balance report with native Hex8 material/support work.
pub type HexThermoplasticReport = FiniteThermoplasticReport<IncrementRecord, ThermalHexEnergy>;

/// One converged thermal, mechanical and observer candidate. No accepted input
/// is mutated; source, late thermal, convergence and cancellation refusals are atomic.
pub struct HexThermoplasticCandidate<'a, 't, F: ThermalSectionSourceFamily<'a>, O> {
    thermal: DepositionSimulation<'t>,
    solid: ThermoplasticHexSolid<'a, F>,
    observer: O,
    report: HexThermoplasticReport,
    quadrature_energy: Vec<ThermalHexMaterialEnergy>,
}
impl<'a, 't, F: ThermalSectionSourceFamily<'a>, O> HexThermoplasticCandidate<'a, 't, F, O> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &ThermoplasticHexSolid<'a, F> {
        &self.solid
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &HexThermoplasticReport {
        &self.report
    }
    pub fn quadrature_energy(&self) -> &[ThermalHexMaterialEnergy] {
        &self.quadrature_energy
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticHexSolid<'a, F>,
        O,
        HexThermoplasticReport,
        Vec<ThermalHexMaterialEnergy>,
    ) {
        (
            self.thermal,
            self.solid,
            self.observer,
            self.report,
            self.quadrature_energy,
        )
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), HexThermoplasticError<E>> {
    cx.checkpoint()
        .map_err(|_| HexThermoplasticError::Cancelled)
}
fn finite<E>(value: f64) -> Result<f64, HexThermoplasticError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(HexThermoplasticError::Invalid(
            "unrepresentable Hex8 physical heat",
        ))
    }
}
fn source_error<E>(point: usize, source: ThermalSectionSourceError) -> HexThermoplasticError<E> {
    HexThermoplasticError::Hex(ThermalHexError::Material { point, source })
}
