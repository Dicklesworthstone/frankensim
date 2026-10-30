//! Source-bound two-axis layer shear and its original finite thermal cell.
//!
//! The original material owner admits affine G(T), constant reference density,
//! unchanged area/thickness and solid caloric state. The original modal owner
//! retains both shear directions and their complete response. Each endpoint
//! trial starts from accepted elastic/plastic history; no trial changes clocks.
//! Plastic heat and signed reversible thermoelastic exchange reach the same
//! declared film cell. No receiving mass, heat partition or material is inferred.

use super::{
    CoupledModalSystem, ModalAttachment, ModalFrictionConfig, PlanarTangentialModalCandidate,
    PlanarTangentialModalDrive, PlanarTangentialModalError, PlanarTangentialModalFrame,
    PlanarTangentialModalSystem,
};
use crate::render::schedule::force::coupled::contact::friction::multiple::layer::{
    material_layer_law,
    thermal::{LayerThermalCell, LayerThermalError, continue_material_layer, deposit_heat},
};
use fs_exec::CancelGate;
use fs_matdb::MaterialCard;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::{MaterialStatePointError, ResolvedShearLayer};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::layer_contact::{ThirdBodyLayerContact, ThirdBodyLayerState};
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::{TangentialLaw, TriboError};
use std::cell::{Cell, RefCell};
use std::convert::Infallible;

pub use crate::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::LayerThermalSolve;

/// One existing finite film cell with its original source and caloric objects.
#[derive(Clone, Copy, Debug)]
pub struct PlanarLayerThermalBinding<'a> {
    pub cell: CellId,
    pub material: &'a MaterialCard,
    pub caloric: &'a EquilibriumEnthalpyPhaseCurve,
    pub temperature_axis: &'a str,
}

impl<'a> PlanarLayerThermalBinding<'a> {
    fn original(self) -> LayerThermalCell<'a> {
        LayerThermalCell {
            port: 0,
            cell: self.cell,
            material: self.material,
            caloric: self.caloric,
            temperature_axis: self.temperature_axis,
        }
    }
}

/// Source, caloric, numerical or original two-axis mechanical refusal.
#[derive(Debug)]
pub enum PlanarMaterialLayerError<E = Infallible> {
    Physics(LayerThermalError<E>),
    Mechanics(PlanarTangentialModalError),
}
impl<E> From<LayerThermalError<E>> for PlanarMaterialLayerError<E> {
    fn from(error: LayerThermalError<E>) -> Self {
        Self::Physics(error)
    }
}
impl<E> From<PlanarTangentialModalError> for PlanarMaterialLayerError<E> {
    fn from(error: PlanarTangentialModalError) -> Self {
        Self::Mechanics(error)
    }
}
impl<E> From<MaterialStatePointError> for PlanarMaterialLayerError<E> {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Physics(error.into())
    }
}
impl<E> From<DepositionError> for PlanarMaterialLayerError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Physics(error.into())
    }
}
impl<E> From<NewtonError> for PlanarMaterialLayerError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Physics(error.into())
    }
}
impl<E> From<TriboError> for PlanarMaterialLayerError<E> {
    fn from(error: TriboError) -> Self {
        Self::Physics(LayerThermalError::Mechanics(error.into()))
    }
}
impl<E: core::fmt::Display> core::fmt::Display for PlanarMaterialLayerError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Physics(error) => write!(f, "planar material layer: {error}"),
            Self::Mechanics(error) => write!(f, "planar material layer: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for PlanarMaterialLayerError<E> {}

fn invalid<E>(message: &'static str) -> PlanarMaterialLayerError<E> {
    LayerThermalError::Invalid(message).into()
}
fn number<E>(value: f64) -> Result<f64, PlanarMaterialLayerError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "unrepresentable planar thermoelastic coordinate or energy",
        ))
    }
}
fn numerical<E>(
    error: FiniteDifferenceError<PlanarMaterialLayerError<E>>,
) -> PlanarMaterialLayerError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => invalid(message),
    }
}

/// Two-axis retained film on the unchanged modal and finite-cell owners.
///
/// Constitutive history lives only in `mechanics`; the source retains the same
/// selected claims, geometry and validity domain throughout thermal refresh.
/// Read-only access prevents advancing the mechanical clock without its film.
pub struct PlanarMaterialLayerSystem<'a> {
    mechanics: PlanarTangentialModalSystem<ThirdBodyLayerContact>,
    source: ResolvedShearLayer,
    binding: PlanarLayerThermalBinding<'a>,
    stiffness_slope_n_m_k: f64,
}

/// Accepted held-effort mechanics and complete physical thermal accounting.
pub struct PlanarLayerThermalFrame<P> {
    /// Mechanical free-energy change evaluated with endpoint coefficients.
    /// Coefficient work and mechanical entropy below complete internal energy.
    pub mechanical: PlanarTangentialModalFrame,
    pub thermal: EnergyAudit,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub mechanical_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Actual representable increment inserted into the original cell sources.
    pub deposited_heat_j: f64,
    /// Native nonnegative plastic loss, transferred internally exactly once.
    pub plastic_heat_j: f64,
    /// Original Maxwell-arm heat, including sub-yield stress relaxation.
    pub relaxation_heat_j: f64,
    /// Signed caloric exchange; never classified as frictional dissipation.
    pub reversible_heat_j: f64,
    /// psi(T_new,e_old) - psi(T_old,e_old).
    pub coefficient_free_energy_change_j: f64,
    /// T_new*s_new - T_old*s_old, with s = -k' |e|^2 / 2.
    pub mechanical_entropy_energy_change_j: f64,
    /// Actual accepted change of psi + T*s in both tangent directions.
    pub layer_internal_energy_change_j: f64,
    /// Independently evaluated other heat, excluding film exchanges.
    pub external_heat_j: f64,
    /// Modal + layer internal + caloric storage + network loss - supplies.
    pub total_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: PlanarTangentialModalCandidate<ThirdBodyLayerState>,
    layer: ResolvedShearLayer,
    law: ThirdBodyLayerContact,
    source: P,
    deposited_heat_j: f64,
    plastic_heat_j: f64,
    relaxation_heat_j: f64,
    reversible_heat_j: f64,
    coefficient_work_j: f64,
    entropy_energy_change_j: f64,
    internal_energy_change_j: f64,
    external_heat_j: f64,
}

impl<'a> PlanarMaterialLayerSystem<'a> {
    /// Bind the original resolved film to its actual mass rho*A*h and caloric
    /// object. Source coverage must include static and sample Nyquist response.
    /// All caloric knots must admit one affine shear modulus and fixed density;
    /// non-affine elasticity requires a different heat-capacity/entropy owner.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        network: CoupledModalSystem,
        left: [ModalAttachment; 2],
        right: Option<[ModalAttachment; 2]>,
        source: ResolvedShearLayer,
        binding: PlanarLayerThermalBinding<'a>,
        thermal: &DepositionSimulation<'_>,
        config: ModalFrictionConfig,
        matrix_config: MatrixTangentialConfig,
        gate: &CancelGate,
    ) -> Result<Self, PlanarMaterialLayerError> {
        if gate.is_requested() {
            return Err(DepositionError::Cancelled.into());
        }
        let start = number(network.samples_rendered() as f64 * network.sample_period_s())?;
        if thermal.time_s() != start {
            return Err(invalid(
                "planar layer and original thermal clocks must agree",
            ));
        }
        source.require_angular_frequency_rad_s(0.0)?;
        source
            .require_angular_frequency_rad_s(core::f64::consts::PI / network.sample_period_s())?;
        let law = material_layer_law(&source).map_err(LayerThermalError::Mechanics)?;
        let admitted =
            continue_material_layer(&source, &law.initial_state(), binding.original(), thermal)?;
        // The first refresh pins the already selected claims and therefore
        // changes its policy receipt even at the identical physical state.
        // The original helper checks source identity, geometry and caloric
        // ownership; require the actual law (including temperature) to agree.
        if admitted.law != law {
            return Err(invalid(
                "initial planar layer source must match its accepted thermal cell",
            ));
        }
        let source = admitted.source;
        let knots = binding.caloric.knots();
        let span = (
            knots
                .first()
                .expect("admitted original caloric chart")
                .temperature_k,
            knots
                .last()
                .expect("admitted original caloric chart")
                .temperature_k,
        );
        let stiffness_slope_n_m_k = source.thermal_stiffness_slope_n_m_k(
            binding.material,
            binding.temperature_axis,
            span,
        )?;
        let mechanics = PlanarTangentialModalSystem::new(
            network,
            left,
            right,
            law,
            config,
            matrix_config,
            gate,
        )?;
        Ok(Self {
            mechanics,
            source,
            binding,
            stiffness_slope_n_m_k,
        })
    }

    #[must_use]
    pub const fn mechanics(&self) -> &PlanarTangentialModalSystem<ThirdBodyLayerContact> {
        &self.mechanics
    }
    #[must_use]
    pub const fn source(&self) -> &ResolvedShearLayer {
        &self.source
    }
    #[must_use]
    pub const fn stiffness_slope_n_m_k(&self) -> f64 {
        self.stiffness_slope_n_m_k
    }

    /// Advance the same retained film and original finite thermal network.
    ///
    /// The film temperature is solved from its actual receiving cell: caller
    /// temperatures refuse. Pressure, when supplied, must equal N/A; N must stay
    /// compressive. The law remains the original pressure-independent material
    /// shear law, without invented pressure sensitivity or a changed footprint.
    /// `resolve` supplies only OTHER heat and original link conductances at the
    /// proposed endpoint. It must not publish history or reinsert film heat.
    /// Guesses use the thermal owner's ascending original CellId order.
    ///
    /// For endpoint k(T), native work is delta(psi)_fixed_T + plastic loss.
    /// The caloric source is D - coefficient_work - delta(T*s). This retains
    /// both vector strains and closes internal energy rather than treating
    /// coefficient changes as dissipation. Fixed reference density, solid phase,
    /// source claims, film geometry and the entire affine caloric domain retain
    /// their original admission. Endpoint thermal feedback is first order and
    /// needs refinement; convergence alone does not certify a unique root.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_layer_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        external: &[f64],
        drive: PlanarTangentialModalDrive,
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>) -> Result<StepCoefficients<P>, E>,
        gate: &CancelGate,
    ) -> Result<PlanarLayerThermalFrame<P>, PlanarMaterialLayerError<E>> {
        let checkpoint = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let network = self.mechanics.network();
        let next = network
            .samples_rendered()
            .checked_add(1)
            .ok_or_else(|| invalid("planar layer sample clock overflow"))?;
        let start = number(network.samples_rendered() as f64 * network.sample_period_s())?;
        let target = number(next as f64 * network.sample_period_s())?;
        let pressure = number(drive.normal_force_n / self.source.area_m2())?;
        if thermal.time_s() != start
            || target <= start
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || thermal.cells().len() != energy_guess_j.len()
            || energy_guess_j.iter().any(|x| !x.is_finite())
            || controls.maximum_mechanical_trials == 0
            || controls.maximum_newton_iterations == 0
            || !controls.difference_step_j.is_finite()
            || controls.difference_step_j <= 0.0
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
            || drive.contact_temperature_k.is_some()
            || pressure <= 0.0
            || drive
                .nominal_pressure_pa
                .is_some_and(|p| p.to_bits() != pressure.to_bits())
        {
            return Err(invalid(
                "synchronized bounded thermal solve and retained compressive planar film required",
            ));
        }
        let original_source = self.source.clone();
        let original_state = self.mechanics.state().clone();
        let binding = self.binding.original();
        let admitted =
            continue_material_layer(&original_source, &original_state, binding, thermal)?;
        if admitted.law != *self.mechanics.law() {
            return Err(invalid(
                "accepted planar layer source no longer matches its original thermal cell",
            ));
        }
        let slope = self.stiffness_slope_n_m_k;
        let initial = self
            .mechanics
            .law()
            .thermodynamics(&original_state, slope)?;
        let network_before_j = self
            .mechanics
            .network()
            .total_energy_j()
            .map_err(PlanarTangentialModalError::from)?;
        let old_temperature = original_source.temperature_kelvin();
        let step = BackwardEulerStep::new(thermal, target, checkpoint)?;
        let dt = step.duration_s();
        let trials = Cell::new(0usize);
        let owner = RefCell::new(&mut self.mechanics);
        let evaluate = |x: &[f64]| {
            step.evaluate(
                x,
                |endpoint| {
                    checkpoint()?;
                    if trials.get() >= controls.maximum_mechanical_trials {
                        return Err(LayerThermalError::MechanicalBudget.into());
                    }
                    trials.set(trials.get() + 1);
                    let layer = continue_material_layer(
                        &original_source,
                        &original_state,
                        binding,
                        endpoint,
                    )?;
                    let candidate_drive = PlanarTangentialModalDrive {
                        contact_temperature_k: Some(layer.source.temperature_kelvin()),
                        nominal_pressure_pa: Some(pressure),
                        ..drive
                    };
                    let mechanical = owner.borrow_mut().stage_candidate_with_law(
                        external,
                        candidate_drive,
                        Some(gate),
                        Some((&layer.law, &layer.state)),
                    )?;
                    let current = layer.law.thermodynamics(&mechanical.state, slope)?;
                    let entropy_energy_change_j = if slope == 0.0 {
                        0.0
                    } else {
                        number(
                            layer.source.temperature_kelvin() * current.entropy_j_k
                                - old_temperature * initial.entropy_j_k,
                        )?
                    };
                    let coefficient_work_j = layer.coefficient_free_energy_change_j;
                    let reversible_heat_j = number(-coefficient_work_j - entropy_energy_change_j)?;
                    let internal_energy_change_j =
                        number(current.internal_energy_j - initial.internal_energy_j)?;
                    let relaxation_heat_j = layer.law.step_dissipation(&layer.state,
                        &mechanical.state, dt)?.relaxation_j;
                    let plastic_heat_j = number(mechanical.frame.friction_dissipation_j-relaxation_heat_j)?;
                    let physical_heat_j = number(mechanical.frame.friction_dissipation_j + reversible_heat_j)?;
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(endpoint).map_err(LayerThermalError::Source)?;
                    let external_heat_j = number(
                        dt * endpoint.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                    )?;
                    let deposited_heat_j = deposit_heat(
                        &mut sources_w,
                        binding.cell,
                        physical_heat_j,
                        dt,
                        controls.total_energy_tolerance_j,
                    )?;
                    if number(deposited_heat_j - physical_heat_j)?.abs()
                        > controls.total_energy_tolerance_j
                    {
                        return Err(invalid(
                            "planar film exchange and actual cell deposition do not balance",
                        ));
                    }
                    Ok::<_, PlanarMaterialLayerError<E>>(StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            mechanical,
                            layer: layer.source,
                            law: layer.law,
                            source: payload,
                            deposited_heat_j,
                            plastic_heat_j,
                            relaxation_heat_j,
                            reversible_heat_j,
                            coefficient_work_j,
                            entropy_energy_change_j,
                            internal_energy_change_j,
                            external_heat_j,
                        },
                    })
                },
                checkpoint,
            )
        };
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|trial| trial.residual_j().to_vec()),
        )
        .map_err(numerical)?;
        let construction =
            NewtonKrylovState::new(&problem, energy_guess_j.to_vec(), controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                gate.is_requested() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(LayerThermalError::Interrupted(result.progress).into());
        }
        if !result.report.converged {
            return Err(LayerThermalError::NotConverged(result.report).into());
        }
        // Restage exactly the accepted endpoint after every numerical probe.
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (next_thermal, staged, thermal_audit) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        let reconstructed_residual_j = number(
            staged.mechanical.frame.energy_residual_j + thermal_audit.stored_energy_change_j
                - staged.plastic_heat_j
                - staged.relaxation_heat_j
                - staged.external_heat_j
                + staged.coefficient_work_j
                + staged.entropy_energy_change_j,
        )?;
        // Independently use actual accepted endpoint internal energies; the
        // algebraic free-energy reconstruction must close under the same gate.
        let frame = &staged.mechanical.frame;
        let total_energy_residual_j = number(
            (frame.network_energy_j - network_before_j)
                + staged.internal_energy_change_j
                + thermal_audit.stored_energy_change_j
                + frame.network_dissipation_j
                - frame.external_work_j
                - frame.drive_work_j
                - staged.external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j
            || reconstructed_residual_j.abs() > controls.total_energy_tolerance_j
        {
            return Err(invalid(
                "planar modal, layer internal and caloric energy do not close",
            ));
        }
        checkpoint()?;
        drop(problem);
        let mechanics = owner.into_inner();
        mechanics.network.validate_staged_material_points(Some(gate))
            .map_err(PlanarTangentialModalError::from)?;
        let mechanical = staged.mechanical.frame;
        // All fallible work precedes joint source/history/body/thermal publication.
        mechanics.law = staged.law;
        mechanics.publish_candidate(staged.mechanical);
        self.source = staged.layer;
        *thermal = next_thermal;
        Ok(PlanarLayerThermalFrame {
            mechanical,
            thermal: thermal_audit,
            source: staged.source,
            nonlinear: result.report,
            mechanical_trials: trials.get(),
            maximum_cell_residual_j,
            deposited_heat_j: staged.deposited_heat_j,
            plastic_heat_j: staged.plastic_heat_j,
            relaxation_heat_j: staged.relaxation_heat_j,
            reversible_heat_j: staged.reversible_heat_j,
            coefficient_free_energy_change_j: staged.coefficient_work_j,
            mechanical_entropy_energy_change_j: staged.entropy_energy_change_j,
            layer_internal_energy_change_j: staged.internal_energy_change_j,
            external_heat_j: staged.external_heat_j,
            total_energy_residual_j,
        })
    }
}
