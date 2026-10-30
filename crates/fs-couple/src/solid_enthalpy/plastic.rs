//! Spatial thermoplastic feedback through the original tetrahedral and heat owners.
//!
//! Each physical interval is one load increment. The default is quasistatic;
//! [`dynamics`] adds consistent source-density inertia on the SAME Newton and
//! thermal feedback paths. The native P1/J2
//! Newton solve uses locally refreshed yield data, retaining the previous
//! displacement and plastic history. Physical yield dissipation and reversible
//! thermoelastic exchange enter the original enthalpy graph once per interval.
//! All coupled proposals restart from the same accepted mechanical, thermal and
//! observer states. Only a converged, balanced pair is returned.
//!
//! This uses the energetic-hardening reduction of ThermoplasticSource: affine
//! sourced E(T)/H(T), fixed nu/reference density, integrated isotropic expansion,
//! fully solid caloric charts, small strain AND rotation. Hardening storage and backward-Euler loss
//! are not heat. Spatial conduction, boundary exchange and existing observers
//! remain in DepositionSimulation. Refine both physical load intervals and
//! thermal substeps; a fixed-point gate is not a time-accuracy certificate.

pub mod contact;
/// Finite electrical storage, native 3D plasticity and one shared thermal solve.
pub mod electrical;
/// Inertial motion, plastic heating and reciprocal thermal feedback together.
pub mod dynamics;
/// Finite source gas, actual loaded plastic history and reciprocal cell heat.
pub mod gas;

use std::collections::{BTreeMap, BTreeSet};
use std::convert::Infallible;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::plastic::J2State;
use fs_material::plastic::energy::{EnergyError, EnergyIncrement, EnergyTolerances};
use fs_material::plastic::source::ResolvedJ2Plasticity;
use fs_material::plastic::thermal::{
    FeedbackError, MAX_FEEDBACK_ITERATIONS, SourceQuery, ThermoplasticSource,
};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::{LoadError3, TetLoadGeometry3};
use fs_solid::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, HeatTransfer, StepConfig, ThermalLaw,
};

#[derive(Debug)]
pub enum ThermoplasticError<E = Infallible> {
    Invalid(&'static str),
    Material(FeedbackError),
    Energy(EnergyError),
    Geometry(LoadError3),
    Solid(SmallStrainTetError),
    Thermal(DepositionError),
    Observer(E),
    Cancelled,
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        combined_residual_j: f64,
    },
}
impl<E> From<DepositionError> for ThermoplasticError<E> {
    fn from(e: DepositionError) -> Self {
        if matches!(e, DepositionError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl<E: std::fmt::Debug> std::fmt::Display for ThermoplasticError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "spatial thermoplastic coupling refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ThermoplasticError<E> {}

/// One cell in connectivity order; different tets may use different cards.
pub struct ThermoplasticCell<'a> {
    pub cell: CellId,
    pub query: SourceQuery<'a>,
}

/// Independent numerical gates for the coupled physical interval.
#[derive(Clone, Copy, Debug)]
pub struct ThermoplasticConfig {
    pub maximum_iterations: usize,
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub yield_tolerance_pa: f64,
    /// Relaxation of trial specific enthalpy in (0,1].
    pub enthalpy_relaxation: f64,
}

#[derive(Clone, Debug)]
struct BoundCell<'a> {
    source: ThermoplasticSource<'a>,
    /// Immutable birth reference for later source-pinned interface bindings.
    reference_h_j_kg: f64,
    /// Expected inventory/identity, not an independently evolved thermal owner.
    accepted: DepositedCell,
    volume_m3: f64,
}

/// Accepted loaded solid, with native displacement/history and immutable geometry.
#[derive(Clone, Debug)]
pub struct ThermoplasticSolid<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    geometry: Arc<TetLoadGeometry3<'a>>,
    cells: Vec<BoundCell<'a>>,
    materials: Vec<ResolvedJ2Plasticity>,
    solution: SmallStrainTetSolution<J2State>,
    budget: TetAssemblyBudget,
    time_s: f64,
    conductor_geometry: Option<electrical::geometry::ConductorGeometryBinding>,
}

impl<'a> ThermoplasticSolid<'a> {
    /// Bind a virgin solid to the actual caloric objects and delivered cell
    /// inventory. Volumes come from the existing load-geometry owner. Density
    /// changes cannot create mass on this fixed reference mesh.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<ThermoplasticCell<'a>>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ThermoplasticError> {
        poll(cx)?;
        if sources.len() != tetrahedra.len() {
            return Err(ThermoplasticError::Invalid(
                "exactly one material cell per tetrahedron required",
            ));
        }
        let face_entries = tetrahedra
            .len()
            .checked_mul(4)
            .ok_or(ThermoplasticError::Invalid("face count overflow"))?;
        let geometry = Arc::new(
            TetLoadGeometry3::new(nodes_m, tetrahedra, budget, face_entries, cx)
                .map_err(ThermoplasticError::Geometry)?,
        );
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(sources.len());
        let mut materials = Vec::with_capacity(sources.len());
        for (element, entry) in sources.into_iter().enumerate() {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&entry.cell)
                .ok_or(ThermoplasticError::Invalid("inactive structural cell"))?;
            if !used.insert(entry.cell)
                || !thermal.uses_material_law(
                    cell.material_id,
                    entry.query.caloric.curve() as &dyn ThermalLaw,
                )
            {
                return Err(ThermoplasticError::Invalid(
                    "unique cell mapping and original caloric object required",
                ));
            }
            let source = ThermoplasticSource::try_new(entry.query, cell.specific_h_j_kg)
                .map_err(ThermoplasticError::Material)?;
            let volume_m3 = geometry.element_volumes_m3()[element];
            if !super::same_measure(cell.geometry.volume_m3, volume_m3)
                || !super::same_measure(cell.mass_kg, source.density_kg_m3() * volume_m3)
            {
                return Err(ThermoplasticError::Invalid(
                    "cell volume and mass must match the native tet and sourced density",
                ));
            }
            materials.push(source.reference().clone());
            cells.push(BoundCell {
                source,
                reference_h_j_kg: cell.specific_h_j_kg,
                accepted: cell.clone(),
                volume_m3,
            });
        }
        let solution = SmallStrainTetSolution {
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            element_strain: vec![[0.0; 6]; tetrahedra.len()],
            element_stress: vec![[0.0; 6]; tetrahedra.len()],
            element_states: materials
                .iter()
                .map(SmallStrainLaw::initial_state)
                .collect(),
            increments: Vec::new(),
        };
        poll(cx)?;
        Ok(Self {
            nodes_m,
            tetrahedra,
            geometry,
            cells,
            materials,
            solution,
            budget,
            time_s: thermal.time_s(),
            conductor_geometry: None,
        })
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn geometry(&self) -> &TetLoadGeometry3<'a> {
        &self.geometry
    }
    pub const fn solution(&self) -> &SmallStrainTetSolution<J2State> {
        &self.solution
    }
    pub fn materials(&self) -> &[ResolvedJ2Plasticity] {
        &self.materials
    }

    /// Select source-defined linear overstress in connectivity order before
    /// loading. `None` explicitly retains rate independence for that element.
    /// Each selected viscosity is pinned on the SAME mechanical/caloric card.
    /// All existing static, inertial, contact and electrical candidates use the
    /// actual common physical interval in the original return map and heat audit.
    pub fn with_linear_overstress(self,
        selections: &[Option<fs_material::state_point::MaterialPropertySelection>],
        cx: &Cx<'_>) -> Result<Self, ThermoplasticError> {
        self.with_overstress(selections, false, cx)
    }

    /// Select material-defined nonlinear power-law flow per tetrahedron before
    /// loading. Viscosity, stress scale and exponent are pinned and refreshed at
    /// each actual temperature. Existing equilibrium, inertial, contact and
    /// electrical candidates receive its consistent tangent, retained plastic
    /// history and physical heat on their original shared clock.
    pub fn with_power_law_overstress(self,
        selections: &[Option<fs_material::state_point::MaterialPropertySelection>],
        cx: &Cx<'_>) -> Result<Self, ThermoplasticError> {
        self.with_overstress(selections, true, cx)
    }

    fn with_overstress(mut self,
        selections: &[Option<fs_material::state_point::MaterialPropertySelection>],
        power_law: bool, cx: &Cx<'_>) -> Result<Self, ThermoplasticError> {
        poll(cx)?;
        if selections.len() != self.cells.len() || !self.solution.increments.is_empty()
            || self.solution.element_states.iter().any(|s| s.alpha != 0.0 || s.plastic_strain != [0.0; 6])
        {
            return Err(ThermoplasticError::Invalid(
                "select one optional overstress source per element before mechanical evolution"));
        }
        for (i, (cell, selection)) in self.cells.iter_mut().zip(selections).enumerate() {
            poll(cx)?;
            if let Some(selection) = selection {
                cell.source = if power_law {
                    cell.source.clone().with_power_law_overstress(selection.clone())
                } else {
                    cell.source.clone().with_linear_overstress(selection.clone())
                }.map_err(ThermoplasticError::Material)?;
                self.materials[i] = cell.source.resolve_enthalpy(cell.accepted.specific_h_j_kg)
                    .map_err(ThermoplasticError::Material)?.1;
            }
        }
        poll(cx)?;
        Ok(self)
    }

    /// Lower these exact structural source bindings to the existing material
    /// conduction network. Entries remain in tetrahedron order; linked cells
    /// require an explicit conductivity selection. Selected transport sources
    /// are pinned and share the native rate-independent admission; strain,
    /// stress or rate dependent conductivity needs an evolving response.
    /// Additional thermal cells retain their own bindings in the complete network.
    ///
    /// The returned bindings borrow this accepted solid as their immutable
    /// card/context anchor. Keep it alive while advancing candidate solids and
    /// retain the prepared network to preserve transport claim identities.
    pub fn thermal_material_cells<'s>(
        &'s self,
        conductivity: &'s [Option<(
            fs_matdb::PropertyKey,
            fs_material::state_point::MaterialPropertySelection,
        )>],
    ) -> Result<Vec<crate::circuit_enthalpy::material::MaterialCell<'s>>, ThermoplasticError> {
        if conductivity.len() != self.cells.len() {
            return Err(ThermoplasticError::Invalid(
                "one transport selection per structural material cell required",
            ));
        }
        self.cells
            .iter()
            .zip(conductivity)
            .map(|(bound, selection)| {
                let conductivity = selection
                    .as_ref()
                    .map(|(key, policy)| {
                        bound
                            .source
                            .pin_thermal_transport(
                                bound.accepted.specific_h_j_kg,
                                key,
                                policy.clone(),
                            )
                            .map(|pinned| (key, pinned))
                            .map_err(ThermoplasticError::Material)
                    })
                    .transpose()?;
                Ok(crate::circuit_enthalpy::material::MaterialCell {
                    cell: bound.accepted.id,
                    material_id: bound.accepted.material_id,
                    curve: bound.source.caloric_source().curve(),
                    point: bound.source.thermal_material_point(),
                    conductivity,
                })
            })
            .collect()
    }

    fn enthalpies<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, ThermoplasticError<E>> {
        let mut enthalpies = Vec::with_capacity(self.cells.len());
        for bound in &self.cells {
            poll(cx)?;
            let old = &bound.accepted;
            let cell = thermal
                .cells()
                .get(&old.id)
                .ok_or(ThermoplasticError::Invalid(
                    "structural cell became inactive",
                ))?;
            if cell.id != old.id
                || cell.segment_id != old.segment_id
                || cell.body_id != old.body_id
                || cell.material_id != old.material_id
                || cell.geometry != old.geometry
                || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
                || accepted && cell.specific_h_j_kg.to_bits() != old.specific_h_j_kg.to_bits()
                || !thermal.uses_material_law(
                    cell.material_id,
                    bound.source.caloric_source().curve() as &dyn ThermalLaw,
                )
            {
                return Err(ThermoplasticError::Invalid(
                    "structural identity, reference geometry, mass or accepted enthalpy changed",
                ));
            }
            bound
                .source
                .resolve_enthalpy(cell.specific_h_j_kg)
                .map_err(ThermoplasticError::Material)?;
            enthalpies.push(cell.specific_h_j_kg);
        }
        Ok(enthalpies)
    }

    #[allow(clippy::too_many_arguments)]
    fn mechanics_with_motion<E>(
        &self,
        enthalpies: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        config: ThermoplasticConfig,
        duration_s: f64,
        motion: Option<&[[f64; 3]]>,
        cx: &Cx<'_>,
    ) -> Result<MechanicalProposal, ThermoplasticError<E>> {
        let mut materials = Vec::with_capacity(self.cells.len());
        let mut laws = Vec::with_capacity(self.cells.len());
        let mut temperatures_k = Vec::with_capacity(self.cells.len());
        for (bound, &h) in self.cells.iter().zip(enthalpies) {
            poll(cx)?;
            let (state, law) = bound
                .source
                .resolve_expanded_enthalpy_for_step(h, duration_s)
                .map_err(ThermoplasticError::Material)?;
            temperatures_k.push(state.temperature_k());
            materials.push(law.resolved().clone());
            laws.push(law);
        }
        let problem = SmallStrainTetProblem {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.budget,
            settings,
        };
        let old_u = self.solution.displacement_m.as_flattened();
        let (solution, motion) = if let Some(velocity) = motion {
            // The original source caloric density and delivered inventory were
            // admitted together at binding. Endpoint temperature never creates mass.
            let densities: Vec<_> = self.cells.iter()
                .map(|cell| cell.source.density_kg_m3()).collect();
            let native = problem.advance_backward_euler(old_u, velocity,
                &self.solution.element_states, &densities, duration_s, cx)
                .map_err(ThermoplasticError::Solid)?;
            (native.solution, Some(dynamics::Motion {
                velocity_mps: native.velocity_mps, energy: native.energy,
            }))
        } else {
            (problem.solve_from(old_u, &self.solution.element_states, cx)
                .map_err(ThermoplasticError::Solid)?, None)
        };
        let gradients = problem
            .displacement_gradients(solution.displacement_m.as_flattened(), cx)
            .map_err(ThermoplasticError::Solid)?;
        let mut energy = Vec::with_capacity(self.cells.len());
        let mut heat_j = Vec::with_capacity(self.cells.len());
        let mut stored_change_j = 0.0;
        let mut numerical_loss_j = 0.0;
        let mut internal_work_j = 0.0;
        let mut physical_heat_j = 0.0;
        let mut reversible_heat_j = 0.0;
        let mut thermoelastic_storage_change_j = 0.0;
        let mut coefficient_storage_change_j = 0.0;
        for (i, bound) in self.cells.iter().enumerate() {
            poll(cx)?;
            let law = materials[i].law();
            let gradient_norm_squared: f64 = gradients[i].iter().flatten().map(|g| g * g).sum();
            if !gradient_norm_squared.is_finite()
                || gradient_norm_squared > law.elastic.strain_limit * law.elastic.strain_limit
            {
                return Err(ThermoplasticError::Invalid(
                    "displacement gradient exceeds the small-strain/rotation bound; finite deformation required",
                ));
            }
            let increment = bound
                .source
                .advance_with_energy_for_step(
                    bound.accepted.specific_h_j_kg,
                    enthalpies[i],
                    &self.solution.element_strain[i],
                    &self.solution.element_states[i],
                    &solution.element_strain[i],
                    EnergyTolerances {
                        energy_j_m3: config.energy_tolerance_j
                            / self.cells.len() as f64
                            / bound.volume_m3,
                        yield_pa: config.yield_tolerance_pa,
                    },
                    duration_s,
                )
                .map_err(ThermoplasticError::Material)?;
            heat_j.push(finite(bound.volume_m3 * increment.caloric_heat_j_m3())?);
            reversible_heat_j =
                finite(reversible_heat_j + bound.volume_m3 * increment.reversible_heat_j_m3)?;
            thermoelastic_storage_change_j = finite(
                thermoelastic_storage_change_j
                    + bound.volume_m3 * increment.thermoelastic_storage_change_j_m3,
            )?;
            internal_work_j =
                finite(internal_work_j + bound.volume_m3 * increment.mechanical_work_j_m3)?;
            coefficient_storage_change_j = finite(
                coefficient_storage_change_j
                    + bound.volume_m3 * increment.coefficient_storage_change_j_m3,
            )?;
            stored_change_j =
                finite(stored_change_j + bound.volume_m3 * increment.stored_change_j_m3())?;
            let account = increment.mechanical;
            if account.state != solution.element_states[i] {
                return Err(ThermoplasticError::Invalid(
                    "native equilibrium and energy update produced different plastic history",
                ));
            }
            physical_heat_j = finite(
                physical_heat_j + bound.volume_m3 * account.energy.physical_dissipation_j_m3,
            )?;
            numerical_loss_j = finite(
                numerical_loss_j + bound.volume_m3 * account.energy.numerical_dissipation_j_m3,
            )?;
            energy.push(account);
        }
        let new_u = solution.displacement_m.as_flattened();
        let mut external_work_j = 0.0;
        for &(dof, force) in nodal_forces_n {
            external_work_j = finite(external_work_j + force * (new_u[dof] - old_u[dof]))?;
        }
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&solution.increments[0].reactions_n) {
            external_work_j = finite(external_work_j + reaction * (new_u[dof] - old_u[dof]))?;
        }
        let inertia_work_j = motion.as_ref().map_or(0.0,
            |motion| motion.energy.kinetic_change_j + motion.energy.numerical_loss_j);
        let mechanical_residual_j = finite(internal_work_j + inertia_work_j - external_work_j)?;
        if mechanical_residual_j.abs() > config.energy_tolerance_j {
            return Err(ThermoplasticError::Invalid(
                "native force equilibrium misses the independent boundary-work gate",
            ));
        }
        Ok(MechanicalProposal {
            materials,
            temperatures_k,
            solution,
            displacement_gradients: gradients,
            energy,
            heat_j,
            physical_heat_j,
            reversible_heat_j,
            thermoelastic_storage_change_j,
            coefficient_storage_change_j,
            stored_change_j,
            numerical_loss_j,
            external_work_j,
            mechanical_residual_j,
            motion,
        })
    }

    /// Solve loaded-solid response and spatial heating as one candidate.
    /// Prescribed displacements and dead forces are endpoint values. Pressure,
    /// traction and body forces can be assembled using `geometry()`'s native
    /// load owner. No input history is mutated on success, refusal or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<ThermoplasticCandidate<'a, 't, O>, ThermoplasticError<O::Error>> {
        self.advance_with_motion(thermal, observer, target_s, prescribed_m,
            nodal_forces_n, settings, thermal_config, config, None, cx).map(|(candidate, _)| candidate)
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_motion<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        velocity_mps: Option<&[[f64; 3]]>,
        cx: &Cx<'_>,
    ) -> Result<(ThermoplasticCandidate<'a, 't, O>, Option<dynamics::Motion>), ThermoplasticError<O::Error>> {
        poll(cx)?;
        let duration_s = target_s - self.time_s;
        if thermal.time_s() != self.time_s
            || !positive(duration_s)
            || !target_s.is_finite()
            || config.maximum_iterations == 0
            || config.maximum_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k)
            || !positive(config.energy_tolerance_j)
            || !positive(config.yield_tolerance_pa)
            || !positive(config.enthalpy_relaxation)
            || config.enthalpy_relaxation > 1.0
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(ThermoplasticError::Invalid(
                "matching clocks, positive physical interval and bounded coupling tolerances required",
            ));
        }
        let mut enthalpies = self.enthalpies(thermal, true, cx)?;
        for (i, material) in self.materials.iter().enumerate() {
            let old_temperature = self.cells[i]
                .source
                .caloric_source()
                .curve()
                .state_at_specific_enthalpy(self.cells[i].accepted.specific_h_j_kg)
                .map_err(|e| ThermoplasticError::Material(FeedbackError::Caloric(e)))?
                .temperature_k();
            let old_strain = self.cells[i]
                .source
                .mechanical_strain(old_temperature, &self.solution.element_strain[i])
                .map_err(ThermoplasticError::Material)?;
            self.cells[i].source
                .validate_committed_history(
                    material,
                    &old_strain,
                    &self.solution.element_states[i],
                    config.yield_tolerance_pa,
                )
                .map_err(ThermoplasticError::Energy)?;
        }
        let mut defects = [f64::INFINITY; 3];
        for iteration in 1..=config.maximum_iterations {
            poll(cx)?;
            let trial = self.mechanics_with_motion(
                &enthalpies,
                prescribed_m,
                nodal_forces_n,
                settings,
                config,
                duration_s,
                velocity_mps,
                cx,
            )?;
            let heat = PlasticHeat {
                solid: self,
                inner: observer.clone(),
                cx,
                watts: self
                    .cells
                    .iter()
                    .zip(&trial.heat_j)
                    .map(|(cell, q)| Ok((cell.accepted.id, finite(*q / duration_s)?)))
                    .collect::<Result<_, ThermoplasticError<O::Error>>>()?,
            };
            let (next_thermal, heat, thermal_audit) =
                thermal.advance_observed(&heat, target_s, thermal_config, || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })?;
            let final_h = self.enthalpies(&next_thermal, false, cx)?;
            let actual = self.mechanics_with_motion(&final_h, prescribed_m,
                nodal_forces_n, settings, config, duration_s,
                velocity_mps, cx)?;
            let applied_heat_j = finite(trial.heat_j.iter().sum())?;
            let caloric_heat_j = finite(actual.heat_j.iter().sum())?;
            let external_heat_j = finite(thermal_audit.external_heat_j - applied_heat_j)?;
            let thermal_residual_j = finite(
                thermal_audit.stored_energy_change_j
                    - thermal_audit.incoming_energy_j
                    - external_heat_j
                    - caloric_heat_j,
            )?;
            let combined_residual_j = finite(
                actual.stored_change_j
                    + actual.thermoelastic_storage_change_j
                    + thermal_audit.stored_energy_change_j
                    + actual.numerical_loss_j
                    + actual.motion.as_ref().map_or(0.0, |motion|
                        motion.energy.kinetic_change_j + motion.energy.numerical_loss_j)
                    - actual.external_work_j
                    - external_heat_j
                    - thermal_audit.incoming_energy_j,
            )?;
            defects = [
                trial
                    .temperatures_k
                    .iter()
                    .zip(&actual.temperatures_k)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max),
                finite(
                    trial
                        .heat_j
                        .iter()
                        .zip(&actual.heat_j)
                        .map(|(a, b)| (a - b).abs())
                        .sum(),
                )?,
                combined_residual_j,
            ];
            if defects[0] <= config.temperature_tolerance_k
                && defects[1] <= config.energy_tolerance_j
                && thermal_residual_j.abs() <= config.energy_tolerance_j
                && combined_residual_j.abs() <= config.energy_tolerance_j
            {
                let mut solid = self.clone();
                solid.time_s = target_s;
                solid.materials = actual.materials;
                solid.solution = actual.solution;
                for cell in &mut solid.cells {
                    cell.accepted = next_thermal.cells()[&cell.accepted.id].clone();
                }
                let report = ThermoplasticReport {
                    iterations: iteration,
                    thermal_audit,
                    element_energy: actual.energy,
                    physical_heat_j: actual.physical_heat_j,
                    reversible_heat_j: actual.reversible_heat_j,
                    thermoelastic_storage_change_j: actual.thermoelastic_storage_change_j,
                    coefficient_storage_change_j: actual.coefficient_storage_change_j,
                    applied_heat_j,
                    external_heat_j,
                    stored_mechanical_change_j: actual.stored_change_j,
                    numerical_loss_j: actual.numerical_loss_j,
                    mechanical_work_j: actual.external_work_j,
                    mechanical_residual_j: actual.mechanical_residual_j,
                    thermal_residual_j,
                    combined_residual_j,
                    temperature_defect_k: defects[0],
                    heat_defect_j: defects[1],
                };
                poll(cx)?;
                return Ok((ThermoplasticCandidate {
                    thermal: next_thermal,
                    solid,
                    observer: heat.inner,
                    report,
                }, actual.motion));
            }
            for (h, final_h) in enthalpies.iter_mut().zip(final_h) {
                *h = finite(
                    (1.0 - config.enthalpy_relaxation) * *h + config.enthalpy_relaxation * final_h,
                )?;
            }
        }
        poll(cx)?;
        Err(ThermoplasticError::NotConverged {
            iterations: config.maximum_iterations,
            temperature_defect_k: defects[0],
            heat_defect_j: defects[1],
            combined_residual_j: defects[2],
        })
    }
}

struct MechanicalProposal {
    materials: Vec<ResolvedJ2Plasticity>,
    temperatures_k: Vec<f64>,
    solution: SmallStrainTetSolution<J2State>,
    displacement_gradients: Vec<[[f64; 3]; 3]>,
    energy: Vec<EnergyIncrement>,
    heat_j: Vec<f64>,
    physical_heat_j: f64,
    reversible_heat_j: f64,
    thermoelastic_storage_change_j: f64,
    coefficient_storage_change_j: f64,
    stored_change_j: f64,
    numerical_loss_j: f64,
    external_work_j: f64,
    mechanical_residual_j: f64,
    motion: Option<dynamics::Motion>,
}

/// Thermal-only audit includes plastic heat and reversible exchange as inputs.
/// The coupled fields separate both from external heat and close nodal work.
#[derive(Clone, Debug)]
pub struct ThermoplasticReport {
    pub iterations: usize,
    pub thermal_audit: EnergyAudit,
    /// Native endpoint-modulus J2 accounts. Their sum needs the explicitly
    /// reported coefficient correction for the actual two-temperature storage.
    pub element_energy: Vec<EnergyIncrement>,
    pub physical_heat_j: f64,
    /// Signed reversible thermoelastic exchange into caloric storage [J].
    pub reversible_heat_j: f64,
    /// Change in T times mechanical entropy, distinct from free energy [J].
    pub thermoelastic_storage_change_j: f64,
    /// Frozen-history E(T)/H(T) free-energy change [J], already included in
    /// stored_mechanical_change_j; never added again to the total balance.
    pub coefficient_storage_change_j: f64,
    /// Net plastic plus reversible caloric source applied to thermal transport.
    pub applied_heat_j: f64,
    pub external_heat_j: f64,
    /// Actual elastic and hardening free-energy change, excluding T s_mechanical.
    pub stored_mechanical_change_j: f64,
    /// Constitutive backward-Euler loss only; inertial loss is reported by dynamics.
    pub numerical_loss_j: f64,
    pub mechanical_work_j: f64,
    /// Internal work plus any inertial work, minus actual force/support work.
    pub mechanical_residual_j: f64,
    pub thermal_residual_j: f64,
    /// Complete caloric/material balance, also including kinetic change and
    /// inertial numerical loss when supplied by the dynamics owner.
    pub combined_residual_j: f64,
    pub temperature_defect_k: f64,
    pub heat_defect_j: f64,
}

pub struct ThermoplasticCandidate<'m, 't, O> {
    thermal: DepositionSimulation<'t>,
    solid: ThermoplasticSolid<'m>,
    observer: O,
    report: ThermoplasticReport,
}
impl<'m, 't, O> ThermoplasticCandidate<'m, 't, O> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &ThermoplasticSolid<'m> {
        &self.solid
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &ThermoplasticReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticSolid<'m>,
        O,
        ThermoplasticReport,
    ) {
        (self.thermal, self.solid, self.observer, self.report)
    }
}

#[derive(Clone)]
struct PlasticHeat<'s, 'm, 'cx, O> {
    solid: &'s ThermoplasticSolid<'m>,
    inner: O,
    cx: &'s Cx<'cx>,
    watts: BTreeMap<CellId, f64>,
}
impl<O: DepositionObserver> DepositionObserver for PlasticHeat<'_, '_, '_, O> {
    type Error = ThermoplasticError<O::Error>;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner.begin().map_err(ThermoplasticError::Observer)
    }
    fn requires_step_interval(&self) -> bool {
        self.inner.requires_step_interval()
    }
    fn step_sources(
        &mut self,
        thermal_start: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .step_sources(thermal_start, start_s, dt_s, end_s, check)
            .map_err(ThermoplasticError::Observer)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.solid.enthalpies(thermal, false, self.cx)?;
        let mut rates = self
            .inner
            .sources(thermal, checkpoint)
            .map_err(ThermoplasticError::Observer)?;
        for (&cell, &watts) in &self.watts {
            checkpoint()?;
            let rate = rates.entry(cell).or_insert(0.0);
            *rate = finite(*rate + watts)?;
        }
        Ok(rates)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner
            .conductances(thermal, checkpoint)
            .map_err(ThermoplasticError::Observer)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner
            .heat_transfers(thermal, checkpoint)
            .map_err(ThermoplasticError::Observer)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.inner
            .limit_timestep_s(thermal, proposed_s, checkpoint)
            .map_err(ThermoplasticError::Observer)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.inner
            .event(event, thermal)
            .map_err(ThermoplasticError::Observer)?;
        self.solid.enthalpies(thermal, false, self.cx).map(|_| ())
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .step(start_s, duration_s, old, thermal, checkpoint)
            .map_err(ThermoplasticError::Observer)?;
        self.solid.enthalpies(thermal, false, self.cx).map(|_| ())
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), ThermoplasticError<E>> {
    cx.checkpoint().map_err(|_| ThermoplasticError::Cancelled)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite<E>(value: f64) -> Result<f64, ThermoplasticError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ThermoplasticError::Invalid(
            "unrepresentable thermoplastic work or heat arithmetic",
        ))
    }
}
