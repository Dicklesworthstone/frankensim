//! Temperature-dependent tensor mechanics on the original finite heat residual.
//!
//! With constant elastic tensors, psi(e,z) has no explicit T derivative. Only
//! the dissipative branch rates 1/tau(T) change, so the original exact tensor
//! trajectory's dashpot loss is the complete mechanical heat source. Every
//! implicit endpoint trial rebuilds the original source operators, starts from
//! the same accepted physical memory, and enters the existing enthalpy/Newton
//! owner. The affine-stiffness extension derives the original tensor entropy
//! and routes reversible heat through this same residual, with physical U and
//! internal coefficient work kept distinct from externally supplied work.

/// Source-backed affine stiffness, entropy and internal caloric exchange.
pub mod affine;
pub use affine::{AffineThermalBudget, TensorThermalExchange};
/// Source-integrated free strain in the original full-memory trajectory.
pub mod expansion;
pub use expansion::ThermalExpansionCell;

use super::super::{SolidError, SolidStep, TensorSolidBody, TensorSolidHistory};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    cell_driver_candidate,
};
use crate::material_circuit::MaterialCircuitSolve;
use crate::tensor_material::{MaterialChangeBudget, TensorMaterialBudget, TensorMaterialHistory};
use fs_exec::CancelGate;
use fs_material::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_material::visco_tensor::{TensorThermalState, source::thermal::AffineTensorThermal};
use fs_solid::linear3::visco::SolidProjectionBudget;
use fs_solver::{NewtonKrylovConfig, NewtonReport};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{CellId, DepositionSimulation, EnergyAudit, ThermalLaw};
use std::collections::{BTreeMap, BTreeSet};

/// One complete original element and its whole finite caloric receiver.
/// Bindings follow original mesh-element order, never reduced mode order.
#[derive(Clone, Copy)]
pub struct RelaxationThermalCell<'m> {
    pub cell: CellId,
    /// Source-resolved Cp/density chart used by the actual thermal cell.
    pub caloric: &'m ResolvedHeatCapacityCurve,
    /// Same original card, full external query and owned temperature axis.
    pub point: ThermalMaterialPoint<'m>,
    /// Complete admitted solid interval, covered by every selected source.
    pub temperature_interval_k: [f64; 2],
}

/// Physical and numerical controls for the native shared cell-energy solve.
#[derive(Clone, Copy, Debug)]
pub struct RelaxationThermalSolve {
    pub maximum_cells: usize,
    /// Complete mechanical proposals, including finite-difference stencils and
    /// the final fresh independent admission; not merely Newton iterations.
    pub maximum_trials: usize,
    pub maximum_iterations: usize,
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    /// Element heat attribution and the complete mechanical/caloric balance.
    pub total_energy_tolerance_j: f64,
    /// Relative geometry/mass correspondence allowance, in [0, 1e-8].
    pub inventory_relative_tolerance: f64,
    pub newton: NewtonKrylovConfig,
    pub projection: SolidProjectionBudget,
    pub dynamics: TensorMaterialBudget,
}
impl RelaxationThermalSolve {
    fn native(self) -> CircuitEnthalpySolve {
        CircuitEnthalpySolve {
            maximum_cells: self.maximum_cells,
            maximum_circuit_trials: self.maximum_trials,
            maximum_newton_iterations: self.maximum_iterations,
            difference_step_j: self.difference_step_j,
            cell_energy_tolerance_j: self.cell_energy_tolerance_j,
            thermal_energy_tolerance_j: self.thermal_energy_tolerance_j,
            total_energy_tolerance_j: self.total_energy_tolerance_j,
            newton: self.newton,
            // The common cell-driver path never invokes a circuit solver.
            circuit: MaterialCircuitSolve {
                energy_tolerance_j: self.total_energy_tolerance_j,
                power_tolerance_w: 0.0,
            },
        }
    }
}

#[derive(Debug)]
pub enum RelaxationThermalError<E> {
    Solid(SolidError),
    Phase(PhaseStateError),
    Material(ThermalMaterialError),
    Expansion(fs_solid::linear3::source::SolidMaterialError),
    HeatSource(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for RelaxationThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tensor relaxation/thermal source: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RelaxationThermalError<E> {}
/// Retains the existing thermal/Newton errors and concrete physical cause.
pub type Error<E> = CircuitEnthalpyError<RelaxationThermalError<E>>;
fn invalid<E>(s: &'static str) -> Error<E> {
    CircuitEnthalpyError::Invalid(s)
}
fn check<E>(gate: &CancelGate) -> Result<(), Error<E>> {
    if gate.is_requested() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn solid<E>(e: SolidError) -> Error<E> {
    CircuitEnthalpyError::Source(RelaxationThermalError::Solid(e))
}
fn finite<E>(x: f64) -> Result<f64, Error<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("nonfinite tensor/thermal arithmetic"))
    }
}

// A private lossless copy of the native history, never user-authored or reset.
// Keeping the source body beside this data avoids a self-referential owner.
struct Memory {
    state: Vec<f64>,
    time: f64,
    steps: u64,
    work: f64,
    parameter_work: f64,
    thermal_coefficient_work: f64,
    heat: f64,
    initial_energy: f64,
}
impl Memory {
    fn retain(history: &TensorMaterialHistory<'_>) -> Self {
        Self {
            state: history.state.clone(),
            time: history.time,
            steps: history.steps,
            work: history.work,
            parameter_work: history.parameter_work,
            thermal_coefficient_work: history.thermal_coefficient_work,
            heat: history.heat,
            initial_energy: history.initial_energy,
        }
    }
}

/// Rebuilt original body and its full accepted native history. Readout borrows
/// this exact body, so factors and retained memory cannot be mixed by callers.
pub struct RelaxingSolid<'s> {
    body: TensorSolidBody<'s>,
    memory: Memory,
    thermodynamics: Option<Vec<TensorThermalState>>,
}
impl<'s> RelaxingSolid<'s> {
    pub const fn body(&self) -> &TensorSolidBody<'s> {
        &self.body
    }
    /// Mechanical psi/S/U of each original element when affine thermal storage
    /// was solved. Kinetic and caloric energy remain separate native stores.
    pub fn thermodynamics(&self) -> Option<&[TensorThermalState]> {
        self.thermodynamics.as_deref()
    }
    pub fn history(&self) -> TensorSolidHistory<'_, 's> {
        TensorSolidHistory {
            solid: &self.body,
            inner: TensorMaterialHistory {
                body: &self.body.body,
                state: self.memory.state.clone(),
                time: self.memory.time,
                steps: self.memory.steps,
                work: self.memory.work,
                parameter_work: self.memory.parameter_work,
                thermal_coefficient_work: self.memory.thermal_coefficient_work,
                heat: self.memory.heat,
                initial_energy: self.memory.initial_energy,
            },
        }
    }
}

/// Exact original mechanical record, native thermal account and final sources.
pub struct RelaxationThermalReport<S> {
    pub mechanical: SolidStep,
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    pub mechanical_trials: usize,
    /// Internal coefficient work and signed reversible heat for affine storage;
    /// None for the constant-stiffness relaxation-only model.
    pub thermodynamics: Option<TensorThermalExchange>,
    /// Total deposited material heat, including signed reversible exchange when
    /// present. Positive dashpot dissipation remains in mechanical.element_heat_j.
    pub deposited_material_heat_j: f64,
    /// Authored and callback heat, excluding the internal material loss.
    pub external_thermal_heat_j: f64,
    pub total_energy_residual_j: f64,
}
/// One complete endpoint: source queries, memory, inertia, clock and heat.
pub struct RelaxationThermalCandidate<'s, 't, S> {
    model: RelaxingSolid<'s>,
    thermal: DepositionSimulation<'t>,
    report: RelaxationThermalReport<S>,
}
impl<'s, 't, S> RelaxationThermalCandidate<'s, 't, S> {
    pub const fn model(&self) -> &RelaxingSolid<'s> {
        &self.model
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &RelaxationThermalReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        RelaxingSolid<'s>,
        DepositionSimulation<'t>,
        RelaxationThermalReport<S>,
    ) {
        (self.model, self.thermal, self.report)
    }
}

struct Trial<'s> {
    model: RelaxingSolid<'s>,
    step: SolidStep,
    energy_change_j: f64,
    element_thermal_heat_j: Vec<f64>,
    thermal_heat_j: f64,
    exchange: Option<TensorThermalExchange>,
}
impl<'s> ThermalDriver for Trial<'s> {
    type Model = RelaxingSolid<'s>;
    type Step = (SolidStep, Option<TensorThermalExchange>);
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.thermal_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.energy_change_j,
            supplied_j: self.step.dynamics.supplied,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.model, (self.step, self.exchange))
    }
}

impl<'s> TensorSolidHistory<'_, 's> {
    /// Close temperature-dependent relaxation, tensor motion, local physical
    /// dissipation and the original finite heat network on ONE implicit step.
    /// The callback supplies independent heat/conductance data at that SAME
    /// endpoint and must not include mechanical heat or mutate live histories.
    /// Existing prepared material conduction/contact/radiation owners may be
    /// used; passive thermal cells remain in the same residual.
    ///
    /// Initial guesses are cell-energy increments [J] in sorted CellId order.
    /// Process pending start events and split at interior events, as required
    /// by the existing BackwardEulerStep. Target events remain pending.
    /// Every trial starts from the same complete accepted mechanical history.
    ///
    /// Endpoint relaxation rates are held during native exact-ZOH mechanics;
    /// coupled time accuracy is first order in general. Refine the physical
    /// interval independently of nonlinear tolerances. Constant elasticity and
    /// density, fixed geometry/supports, small strain, fully solid caloric cells
    /// and source-supported temperatures are required. Thermal expansion,
    /// stiffness-dependent entropy, phase change and basis adaptation remain
    /// outside this dissipative-kinetics model. No thermal parameter work is
    /// introduced, and prior external parameter-work accounts are preserved.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_relaxation_thermal_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        force_n: &[[f64; 3]],
        cells: &[RelaxationThermalCell<'_>],
        controls: RelaxationThermalSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        gate: &CancelGate,
    ) -> Result<RelaxationThermalCandidate<'s, 't, S>, Error<E>> {
        self.advance_thermal_model(
            thermal, target_s, force_n, cells, controls, guess_j, resolve, None, None, gate,
        )
    }

    /// Close source-dependent stored tensors and relaxation on the SAME native
    /// tensor/enthalpy residual. All original source knots must establish affine
    /// stored tensors; fixed density, geometry, supports and eigenstrain remain
    /// required. Isotropic E(T) uses fixed nu; orthotropic E/G share an affine
    /// scale with fixed Poisson ratios. Branch tensors may have independent
    /// affine entries and preserve complete physical Maxwell memory.
    ///
    /// With psi=.5 e:C(T):e, S=-psi_T and U=psi+T*S, the caloric source is
    /// dashpot heat minus native coefficient work minus delta(T*S). Reversible
    /// exchange may cool a receiver. Source-operator work is internal and never
    /// increments external parameter_work_j. Mechanical physical energy uses U;
    /// the native mechanical-potential audit retains psi and its own distinct
    /// internal coefficient account. Caloric and kinetic stores are counted once.
    ///
    /// This endpoint-coefficient split remains first order. It adds no thermal
    /// expansion, phase transition, moving geometry or nonlinear C(T) entropy
    /// model. Varying singular branch factors retain the native refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_affine_thermal_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        force_n: &[[f64; 3]],
        cells: &[RelaxationThermalCell<'_>],
        controls: RelaxationThermalSolve,
        budget: AffineThermalBudget,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        gate: &CancelGate,
    ) -> Result<RelaxationThermalCandidate<'s, 't, S>, Error<E>> {
        self.advance_thermal_model(
            thermal,
            target_s,
            force_n,
            cells,
            controls,
            guess_j,
            resolve,
            Some(budget),
            None,
            gate,
        )
    }

    /// Solve actual thermal expansion, inertia, Maxwell memory and reciprocal
    /// caloric heat on the SAME implicit residual. Expansion comes from the
    /// existing isotropic/orthotropic source integrals in fixed mesh frames.
    /// Every interval preserves viscous strain when free strain changes; it
    /// never reinitializes branch history or replaces motion by equilibrium.
    ///
    /// The admitted model has constant stored tensors and constant alpha over
    /// the complete source interval, with source-dependent tau(T) allowed.
    /// The native PHS carries .5*(e-e_th):C:(e-e_th) as affine quadratic storage;
    /// S=V*alpha:sigma includes every branch and U=psi+T*S closes caloric energy.
    /// The actual caloric chart must keep mass*dh/dT-T*V*alpha:C_inst:alpha > 0.
    /// Thermal force is internal, never counted as external port work.
    ///
    /// This is the existing first-order endpoint split on fixed reference
    /// geometry and mass, with small total/elastic/free strains. Varying alpha,
    /// varying stored tensors with expansion, Kelvin thermal strain-rate work,
    /// finite deformation and phase transitions require another admitted model.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_expanding_thermal_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        force_n: &[[f64; 3]],
        cells: &[RelaxationThermalCell<'_>],
        expansion: &[ThermalExpansionCell<'_>],
        controls: RelaxationThermalSolve,
        budget: AffineThermalBudget,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        gate: &CancelGate,
    ) -> Result<RelaxationThermalCandidate<'s, 't, S>, Error<E>> {
        self.advance_thermal_model(
            thermal,
            target_s,
            force_n,
            cells,
            controls,
            guess_j,
            resolve,
            Some(budget),
            Some(expansion),
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_thermal_model<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        force_n: &[[f64; 3]],
        cells: &[RelaxationThermalCell<'_>],
        controls: RelaxationThermalSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        affine: Option<AffineThermalBudget>,
        expansion: Option<&[ThermalExpansionCell<'_>]>,
        gate: &CancelGate,
    ) -> Result<RelaxationThermalCandidate<'s, 't, S>, Error<E>> {
        check(gate)?;
        if expansion.is_none()
            && (self.solid.thermal_expansion_identity.is_some()
                || self.solid.body.has_eigenstrain())
        {
            return Err(invalid(
                "an expanding tensor history must retain its original reciprocal expansion model",
            ));
        }
        if affine.is_some_and(|b| {
            !b.max_absolute_coefficient_work_j.is_finite()
                || b.max_absolute_coefficient_work_j < 0.0
                || !b.relative_memory_tolerance.is_finite()
                || b.relative_memory_tolerance <= 0.0
                || b.relative_memory_tolerance > 1e-8
        }) {
            return Err(invalid(
                "affine thermal feedback needs bounded coefficient change and physical-memory reconstruction controls",
            ));
        }
        let sources = self.solid.material_sources().ok_or_else(|| {
            invalid("thermal relaxation requires original tensor source bindings")
        })?;
        if self.inner.time_s().to_bits() != thermal.time_s().to_bits()
            || cells.len() != sources.len()
            || !controls.inventory_relative_tolerance.is_finite()
            || !(0.0..=1e-8).contains(&controls.inventory_relative_tolerance)
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j <= 0.0
        {
            return Err(invalid(
                "one shared clock, one receiver per original element and bounded physical tolerances required",
            ));
        }
        crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls.native(), guess_j)?;
        let mut volumes = vec![0.0; sources.len()];
        for (point, &element) in self
            .solid
            .projection
            .material_points()
            .iter()
            .zip(self.solid.projection.point_elements())
        {
            check(gate)?;
            volumes[element] = finite(volumes[element] + point.volume_m3())?;
        }
        let mut seen = BTreeSet::new();
        let mut derivatives = Vec::with_capacity(sources.len());
        let mut temperatures_before = Vec::with_capacity(sources.len());
        for ((binding, source), &volume) in cells.iter().zip(sources).zip(&volumes) {
            check(gate)?;
            if !seen.insert(binding.cell) {
                return Err(invalid(
                    "each whole-element caloric receiver must be unique",
                ));
            }
            let cell = thermal
                .cells()
                .get(&binding.cell)
                .ok_or_else(|| invalid("inactive element caloric receiver"))?;
            let curve = binding.caloric.curve();
            if !thermal.uses_material_law(cell.material_id, curve as &dyn ThermalLaw)
                || curve.material_card_identity() != source.material.state().card_identity()
                || binding.point.material.content_hash() != source.material.state().card_identity()
            {
                return Err(invalid(
                    "tensor and caloric owners must retain the same original card and exact borrowed caloric law",
                ));
            }
            let state = curve
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(|e| CircuitEnthalpyError::Source(RelaxationThermalError::Phase(e)))?;
            let query = binding
                .point
                .at(state)
                .map_err(|e| CircuitEnthalpyError::Source(RelaxationThermalError::Material(e)))?;
            if query != source.material.state().query() {
                return Err(invalid(
                    "tensor and caloric initial temperatures and complete physical queries must agree",
                ));
            }
            if affine.is_some() {
                let law = AffineTensorThermal::try_new(
                    &source.material,
                    binding.point.material,
                    binding.point.temperature_axis,
                    binding.temperature_interval_k,
                )
                .map_err(|e| solid(SolidError::MaterialSource(e)))?;
                derivatives.push(
                    law.in_frame(source.material_to_mesh)
                        .map_err(|e| solid(SolidError::MaterialSource(e)))?,
                );
            } else {
                source
                    .material
                    .admit_thermal_relaxation(
                        binding.point.material,
                        binding.point.temperature_axis,
                        binding.temperature_interval_k,
                    )
                    .map_err(|e| solid(SolidError::MaterialSource(e)))?;
            }
            temperatures_before.push(state.temperature_k());
            let [low, high] = binding.temperature_interval_k;
            let density = source.material.density_kg_m3();
            let knots = curve.knots();
            if knots[0].temperature_k > low
                || knots[knots.len() - 1].temperature_k < high
                || knots.iter().any(|k| {
                    k.liquid_mass_fraction != 0.0
                        || k.bulk_density_kg_m3.to_bits() != density.to_bits()
                })
                || knots
                    .windows(2)
                    .any(|k| k[1].temperature_k <= k[0].temperature_k)
                || volume <= 0.0
                || ((cell.geometry.volume_m3 / volume) - 1.0).abs()
                    > controls.inventory_relative_tolerance
                || ((cell.mass_kg / (volume * density)) - 1.0).abs()
                    > controls.inventory_relative_tolerance
            {
                return Err(invalid(
                    "caloric law must cover the solid source interval with the original element volume and constant reference mass",
                ));
            }
            let density_claim = source
                .material
                .state()
                .property("density")
                .unwrap()
                .answer()
                .receipt
                .selected;
            if binding.caloric.receipts().is_empty()
                || binding.caloric.receipts().len() % 2 != 0
                || binding
                    .caloric
                    .receipts()
                    .chunks_exact(2)
                    .any(|pair| pair[1].property != "density" || pair[1].selected != density_claim)
                || binding.caloric.receipts().iter().any(|receipt| {
                    receipt.axis_quantities != *source.material.state().axis_quantities()
                        || !receipt
                            .query_point
                            .iter()
                            .filter(|(a, _)| a != binding.point.temperature_axis)
                            .eq(source
                                .material
                                .state()
                                .query_point()
                                .iter()
                                .filter(|(a, _)| a != binding.point.temperature_axis))
                        || receipt.property == "density" && receipt.selected != density_claim
                })
            {
                return Err(invalid(
                    "caloric source receipts must preserve the complete original context and density claim",
                ));
            }
        }
        let expansion_model = expansion
            .map(|bindings| expansion::prepare(self, thermal, cells, bindings, gate))
            .transpose()?;
        let cards: Vec<_> = cells.iter().map(|b| b.point.material).collect();
        let before_states = affine
            .map(|budget| {
                affine::states(
                    self,
                    &derivatives,
                    &temperatures_before,
                    expansion_model.as_ref().map(|m| m.alpha_mandel.as_slice()),
                    budget.relative_memory_tolerance,
                    gate,
                )
            })
            .transpose()?;
        let before = self.inner.stored_energy_j();
        let dt = target_s - thermal.time_s();
        let (thermal, model, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls.native(),
            guess_j,
            |trial, _remaining| {
                check(gate)?;
                let mut points = Vec::with_capacity(cells.len());
                let mut temperatures_after = Vec::with_capacity(cells.len());
                for binding in cells {
                    check(gate)?;
                    let cell = trial.cells().get(&binding.cell).ok_or_else(|| {
                        invalid("element receiver disappeared during thermal proposal")
                    })?;
                    let phase = binding
                        .caloric
                        .curve()
                        .state_at_specific_enthalpy(cell.specific_h_j_kg)
                        .map_err(|e| {
                            CircuitEnthalpyError::Source(RelaxationThermalError::Phase(e))
                        })?;
                    let [low, high] = binding.temperature_interval_k;
                    if !phase.is_fully_solid() || !(low..=high).contains(&phase.temperature_k()) {
                        return Err(invalid(
                            "thermal trial left the admitted fully-solid relaxation interval",
                        ));
                    }
                    points.push(binding.point.at(phase).map_err(|e| {
                        CircuitEnthalpyError::Source(RelaxationThermalError::Material(e))
                    })?);
                    temperatures_after.push(phase.temperature_k());
                }
                let strains = expansion
                    .map(|bindings| expansion::strains(self, cells, bindings, &points, gate))
                    .transpose()?;
                let mut body = self
                    .solid
                    .at_query_points_with_thermal_strains(
                        &cards,
                        &points,
                        strains.as_deref(),
                        controls.projection,
                        controls.dynamics,
                        gate,
                    )
                    .map_err(solid)?;
                if let Some(model) = &expansion_model {
                    body.thermal_expansion_identity = Some(model.identities.clone());
                }
                let rebound = self
                    .rebind_materials(
                        &body,
                        MaterialChangeBudget {
                            max_absolute_parameter_work_j: affine
                                .map_or(0.0, |b| b.max_absolute_coefficient_work_j),
                            relative_memory_tolerance: affine
                                .map_or(1e-12, |b| b.relative_memory_tolerance),
                        },
                        gate,
                    )
                    .map_err(solid)?;
                let coefficient_work = rebound.parameter_work_j;
                let mut history = rebound.history;
                let frozen_states = affine
                    .map(|budget| {
                        affine::states(
                            &history,
                            &derivatives,
                            &temperatures_after,
                            expansion_model.as_ref().map(|m| m.alpha_mandel.as_slice()),
                            budget.relative_memory_tolerance,
                            gate,
                        )
                    })
                    .transpose()?;
                if affine.is_some() {
                    history.inner.parameter_work = self.inner.parameter_work;
                    history.inner.thermal_coefficient_work =
                        finite(self.inner.thermal_coefficient_work + coefficient_work)?;
                }
                let step = history
                    .advance(force_n, dt, controls.total_energy_tolerance_j, gate)
                    .map_err(solid)?;
                if history.inner.time_s().to_bits() != target_s.to_bits() {
                    return Err(invalid(
                        "tensor and thermal endpoints must share the exact accepted clock",
                    ));
                }
                let mut change = finite(history.inner.stored_energy_j() - before)?;
                let after_states = affine
                    .map(|budget| {
                        affine::states(
                            &history,
                            &derivatives,
                            &temperatures_after,
                            expansion_model.as_ref().map(|m| m.alpha_mandel.as_slice()),
                            budget.relative_memory_tolerance,
                            gate,
                        )
                    })
                    .transpose()?;
                let exchange = match (&before_states, &frozen_states, &after_states) {
                    (Some(before), Some(frozen), Some(after)) => Some(affine::exchange(
                        before,
                        frozen,
                        after,
                        &temperatures_before,
                        &temperatures_after,
                        coefficient_work,
                        controls.total_energy_tolerance_j,
                        gate,
                    )?),
                    _ => None,
                };
                let mut element_thermal_heat_j = step.element_heat_j.clone();
                if let Some(exchange) = &exchange {
                    for ((heat, &reversible), &ts) in element_thermal_heat_j
                        .iter_mut()
                        .zip(&exchange.reversible_heat_j)
                        .zip(&exchange.entropy_energy_change_j)
                    {
                        *heat = finite(*heat + reversible)?;
                        change = finite(change + ts)?;
                    }
                }
                let thermal_heat_j = element_thermal_heat_j
                    .iter()
                    .try_fold(0.0, |sum, &v| finite(sum + v))?;
                let memory = Memory::retain(&history.inner);
                let coefficients = resolve(trial).map_err(|e| {
                    CircuitEnthalpyError::Source(RelaxationThermalError::HeatSource(e))
                })?;
                check(gate)?;
                Ok(DrivenCoefficients {
                    circuit: Trial {
                        model: RelaxingSolid {
                            body,
                            memory,
                            thermodynamics: after_states,
                        },
                        step,
                        energy_change_j: change,
                        element_thermal_heat_j,
                        thermal_heat_j,
                        exchange,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.payload,
                    circuit_trials: 1,
                })
            },
            |candidate, duration, checkpoint| {
                let mut powers = BTreeMap::new();
                for (binding, &heat) in cells.iter().zip(&candidate.element_thermal_heat_j) {
                    checkpoint()?;
                    let power = finite(heat / duration)?;
                    if heat != 0.0 && power == 0.0 {
                        return Err(invalid("tensor material heat underflows power conversion"));
                    }
                    powers.insert(binding.cell, power);
                }
                Ok(powers)
            },
            || gate.is_requested(),
        )?;
        check(gate)?;
        Ok(RelaxationThermalCandidate {
            model,
            thermal,
            report: RelaxationThermalReport {
                mechanical: report.electrical.0,
                thermodynamics: report.electrical.1,
                thermal: report.thermal,
                source: report.source,
                nonlinear: report.nonlinear,
                maximum_cell_residual_j: report.maximum_cell_residual_j,
                mechanical_trials: report.circuit_trials,
                deposited_material_heat_j: report.deposited_joule_j,
                external_thermal_heat_j: report.non_electrical_heat_j,
                total_energy_residual_j: report.total_energy_residual_j,
            },
        })
    }
}
