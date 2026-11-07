//! Finite spatial solid heat storage on the original resolved NSF wall faces.
//!
//! Original native fluid steps and one thermal backward-Euler residual share
//! the SAME per-face integrated heat. The solid endpoint supplies each wall's
//! temperature and optional material half-layer resistance; the fluid's native
//! molecular wall solve supplies its actual heat. No boundary-layer correlation,
//! additional gas heat capacity or replacement fluid time integrator is used.
//!
//! Fixed cells, fixed material laws and an integer lattice interval are admitted.
//! The wall temperature/resistance is held at the proposed solid endpoint during
//! that fluid interval, so the coupled time approximation remains first order.
//! Molecular, local and dynamic SGS use the same reversible native wall trial.
//! Dynamic channel momentum/heat fits consume the proposed wall temperatures;
//! molecular conduction still supplies the sole interface heat. Original
//! exploratory closure limits remain. Native fluid source callbacks cannot
//! participate in the wall transaction and remain outside this owner.

use std::cell::{Cell, RefCell};
use std::collections::{BTreeMap, BTreeSet};

use fs_exec::Cx;
use fs_lbm::lattice3::Lattice3;
use fs_lbm::nsf3::{
    EnergyLedger, NsfError, NsfGrid3, NsfThermalWallExchange, NsfThermalWallFace,
    NsfThermalWallTrial, NsfWallThermal,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BoundaryStepCoefficients, StepCoefficients,
};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, ThermalLaw,
};

/// Normal resistance from the original solid cell to its fluid interface.
#[derive(Clone, Copy)]
pub enum NsfSolidLayer<'a> {
    /// Explicitly neglect the solid half-cell and contact resistance. The finite
    /// thermal cell still owns all storage; it is not an infinite reservoir.
    Isothermal,
    /// The original source supplies k at this cell's actual endpoint enthalpy.
    /// k is held through the declared half-layer; no integrated k(T) profile or
    /// phase boundary inside that layer is inferred. Zero k insulates.
    Conducting {
        half_distance_m: f64,
        contact_resistance_m2_k_w: f64,
        conductivity: &'a dyn ConductivityLaw,
    },
}

/// One complete native fluid face coupled to an original finite thermal cell.
/// Its area is exactly the original fluid dx². Several disjoint faces may share
/// a solid cell, but their area must fit its originally exposed ambient area.
#[derive(Clone, Copy)]
pub struct NsfSolidContact<'a> {
    pub fluid_cell: [usize; 3],
    pub axis: usize,
    pub high: bool,
    pub solid_cell: CellId,
    /// Must be the SAME borrowed caloric object registered for `solid_cell`.
    pub caloric: &'a dyn ThermalLaw,
    pub layer: NsfSolidLayer<'a>,
}

impl NsfSolidContact<'_> {
    fn key(&self) -> ([usize; 3], usize, bool) {
        (self.fluid_cell, self.axis, self.high)
    }
}

/// Explicit work ceilings and independently enforced physical joule limits.
#[derive(Clone, Copy, Debug)]
pub struct NsfThermalControls {
    pub maximum_fluid_cells: usize,
    pub maximum_solid_cells: usize,
    pub maximum_contacts: usize,
    /// Includes every residual, numerical derivative and final acceptance query.
    pub maximum_endpoint_evaluations: usize,
    pub maximum_newton_iterations: usize,
    /// Unknowns are native cell-energy increments divided by this positive scale.
    pub energy_scale_j: f64,
    /// Dimensionless directional-difference increment in those scaled coordinates.
    pub finite_difference_step: f64,
    pub cell_tolerance_j: f64,
    pub thermal_tolerance_j: f64,
    pub fluid_tolerance_j: f64,
    pub interface_tolerance_j: f64,
    pub total_tolerance_j: f64,
    pub solver: NewtonKrylovConfig,
}

/// Separate native fluid resolution from the finite solid's update interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NsfThermalIntervalControls {
    /// Positive number of original fluid steps in this one solid interval.
    pub fluid_steps: usize,
    /// Whole-call native cell-step cap, including EVERY nonlinear/derivative
    /// probe and final acceptance. Solid cells in the fluid grid also count.
    /// Each complete candidate reserves its entire interval before starting.
    pub maximum_fluid_cell_steps: usize,
}

/// The native face debit and its one finite solid receiver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfSolidExchange {
    pub face: NsfThermalWallExchange,
    pub solid_cell: CellId,
}

/// Accounts of the same accepted fluid/solid interval.
#[derive(Debug)]
pub struct NsfThermalReport<T> {
    pub fluid: EnergyLedger,
    pub thermal: EnergyAudit,
    /// Accepted original noninterface heat/source payload.
    pub heat: T,
    pub exchanges: Vec<NsfSolidExchange>,
    pub coupled_heat_into_fluid_j: f64,
    /// Original solid ambient, reservoir and authored/resolved external heat,
    /// evaluated independently before adding any fluid interface debit.
    pub solid_external_heat_j: f64,
    pub fluid_energy_residual_j: f64,
    pub interface_energy_residual_j: f64,
    pub maximum_cell_interface_residual_j: f64,
    /// Sum of independently gated original fluid/solid storage defects and
    /// the interface transfer defect, in their retained thermal coordinates.
    pub total_energy_residual_j: f64,
    pub maximum_cell_residual_j: f64,
    pub endpoint_evaluations: usize,
    /// Accepted original fluid steps in this one thermal interval.
    pub fluid_steps: usize,
    /// Reserved native cell-steps across all candidates, including probes.
    pub fluid_cell_steps: usize,
    /// Physical interval from the original lattice duration, not rounded clock subtraction.
    pub duration_s: f64,
    pub nonlinear: NewtonReport,
}

#[derive(Debug)]
pub enum NsfThermalError<E> {
    Invalid(&'static str),
    Fluid(NsfError),
    Thermal(DepositionError),
    Resolve(E),
    Solver(NewtonError),
    NotConverged(NewtonReport),
    WorkBudget {
        performed: usize,
        maximum: usize,
    },
    FluidWorkBudget {
        performed: usize,
        requested: usize,
        maximum: usize,
    },
    Energy {
        component: &'static str,
        residual_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}

impl<E> From<DepositionError> for NsfThermalError<E> {
    fn from(error: DepositionError) -> Self {
        if error == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for NsfThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved fluid/finite solid heat coupling: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for NsfThermalError<E> {}

struct Prepared<'a> {
    contacts: Vec<NsfSolidContact<'a>>,
    original_walls: Vec<NsfThermalWallFace>,
    remaining_ambient: BTreeMap<CellId, f64>,
}

struct Sample<T> {
    heat: T,
    exchanges: Vec<NsfSolidExchange>,
    coupled_heat_j: f64,
    external_heat_j: f64,
    interface_residual_j: f64,
    maximum_cell_interface_residual_j: f64,
}

/// Advance both ORIGINAL owners by one resolved fluid lattice step.
///
/// Thermal time must equal `fluid.steps() * fluid.law().dt_s`; the next endpoint
/// is the next representable lattice timestamp. Native pending/interior thermal
/// events still refuse and must be processed explicitly outside this interval.
/// All source and flow proposals start from the same accepted states.
///
/// `resolve` supplies only independent noninterface heat and original link or
/// reservoir coefficients. Coupled cells' ambient conductances belong to this
/// adapter: it removes exactly the mapped native face area from their authored
/// ambient exchange. A callback overriding those same ambient entries refuses.
/// The callback must be pure and preserve its own source/phase applicability.
///
/// The native per-face integrated heat is added once, with opposite sign, to
/// each finite solid cell. The original fluid and thermal solvers retain their
/// calorics, spatial transport, mechanical work, open boundaries and physical
/// admissions. The final original trials publish together only after local,
/// fluid, thermal, interface and combined energy gates; all earlier failures or
/// cancellation leave both owners unchanged. Refine space and the lattice dt
/// independently of the nonlinear tolerances; no temporal error bound is implied.
#[allow(clippy::too_many_arguments)]
pub fn advance_nsf_solid<'m, L: Lattice3<27>, T, E>(
    cx: &Cx<'_>,
    fluid: &mut NsfGrid3<L>,
    thermal: &mut DepositionSimulation<'m>,
    contacts: &[NsfSolidContact<'_>],
    controls: NsfThermalControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&Cx<'_>, &DepositionSimulation<'m>) -> Result<BoundaryStepCoefficients<T>, E>,
) -> Result<NsfThermalReport<T>, NsfThermalError<E>> {
    let maximum_fluid_cell_steps = controls
        .maximum_fluid_cells
        .checked_mul(controls.maximum_endpoint_evaluations)
        .ok_or(NsfThermalError::Invalid("fluid work ceiling overflow"))?;
    advance_nsf_solid_interval(
        cx,
        fluid,
        thermal,
        contacts,
        controls,
        NsfThermalIntervalControls {
            fluid_steps: 1,
            maximum_fluid_cell_steps,
        },
        energy_guess_j,
        resolve,
    )
}

/// Advance a finite solid over several ORIGINAL native fluid steps.
///
/// Every candidate starts from the same accepted fluid/solid histories. The
/// proposed solid endpoint is held while the native fluid advances at its
/// original stable lattice timestep. Actual chronological face heat is summed
/// and debited once in the solid's original backward-Euler residual; no heat
/// extrapolation or larger fluid timestep is introduced. Independent solid
/// sources and conductances are sampled at the same proposed endpoint.
///
/// This separates update rates, not accuracy requirements: increase the number
/// of native steps only when the solid boundary varies slowly enough. The
/// coupling remains first order; nonlinear convergence is not a time-error
/// estimate. Interior thermal events require a separate interval. Shared work
/// and energy gates cover the complete candidate, and publication is atomic.
#[allow(clippy::too_many_arguments)]
pub fn advance_nsf_solid_interval<'m, L: Lattice3<27>, T, E>(
    cx: &Cx<'_>,
    fluid: &mut NsfGrid3<L>,
    thermal: &mut DepositionSimulation<'m>,
    contacts: &[NsfSolidContact<'_>],
    controls: NsfThermalControls,
    interval: NsfThermalIntervalControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&Cx<'_>, &DepositionSimulation<'m>) -> Result<BoundaryStepCoefficients<T>, E>,
) -> Result<NsfThermalReport<T>, NsfThermalError<E>> {
    checkpoint(cx)?;
    controls
        .solver
        .validate()
        .map_err(NsfThermalError::Solver)?;
    if controls.maximum_fluid_cells == 0
        || controls.maximum_solid_cells == 0
        || controls.maximum_contacts == 0
        || controls.maximum_endpoint_evaluations == 0
        || controls.maximum_newton_iterations == 0
        || interval.fluid_steps == 0
        || interval.maximum_fluid_cell_steps == 0
        || controls
            .maximum_fluid_cells
            .checked_mul(controls.maximum_endpoint_evaluations)
            .is_none()
        || controls.maximum_fluid_cells.checked_mul(6).is_none()
        || fluid.thermal_walls().len() > controls.maximum_fluid_cells.saturating_mul(6)
        || thermal.cells().is_empty()
        || thermal.cells().len() > controls.maximum_solid_cells
        || !positive(controls.energy_scale_j)
        || !positive(controls.finite_difference_step)
        || ![
            controls.cell_tolerance_j,
            controls.thermal_tolerance_j,
            controls.fluid_tolerance_j,
            controls.interface_tolerance_j,
            controls.total_tolerance_j,
        ]
        .into_iter()
        .all(nonnegative)
    {
        return Err(NsfThermalError::Invalid(
            "finite physical controls and bounded fluid/solid work required",
        ));
    }
    let next_index = fluid
        .steps()
        .checked_add(interval.fluid_steps)
        .ok_or(NsfThermalError::Invalid("fluid lattice clock overflow"))?;
    if next_index as u128 > (1_u128 << 53) {
        return Err(NsfThermalError::Invalid(
            "fluid lattice clock is no longer exactly representable",
        ));
    }
    let dt = fluid.law().dt_s;
    let start_s = fluid.steps() as f64 * dt;
    let target_s = next_index as f64 * dt;
    let duration_s = interval.fluid_steps as f64 * dt;
    if !positive(dt)
        || !positive(duration_s)
        || thermal.time_s() != start_s
        || !target_s.is_finite()
        || target_s <= start_s
    {
        return Err(NsfThermalError::Invalid(
            "fluid and solid must share the original representable lattice clock",
        ));
    }
    let cell_steps = fluid
        .cell_count()
        .checked_mul(interval.fluid_steps)
        .ok_or(NsfThermalError::Invalid("native cell-step work overflow"))?;
    let prepared = Prepared::new(cx, fluid, thermal, contacts, controls.maximum_contacts)?;
    let step = BackwardEulerStep::new_with_duration(thermal, target_s, duration_s, || {
        thermal_checkpoint(cx)
    })?;
    let count = step.cell_ids().len();
    if energy_guess_j.len() != count || energy_guess_j.iter().any(|value| !value.is_finite()) {
        return Err(NsfThermalError::Invalid(
            "complete finite solid cell-energy guess required",
        ));
    }
    let evaluations = Cell::new(0usize);
    let fluid_work = Cell::new(0usize);
    let fluid = RefCell::new(fluid);
    let solution;
    let nonlinear;
    {
        let evaluate = |x: &[f64]| {
            charge(
                cx,
                &evaluations,
                controls.maximum_endpoint_evaluations,
                &fluid_work,
                cell_steps,
                interval.maximum_fluid_cell_steps,
            )?;
            let energies = physical_increments(&step, x, controls.energy_scale_j)?;
            step.evaluate_with_boundary_conductances(
                &energies,
                |candidate| {
                    let (walls, mut heat) = prepared.endpoint(cx, candidate, &resolve)?;
                    let mut owner = fluid.borrow_mut();
                    let trial = owner
                        .preview_thermal_wall_interval(
                            walls,
                            interval.fluid_steps,
                            controls.maximum_fluid_cells,
                            cell_steps,
                            || cx.checkpoint().is_err(),
                        )
                        .map_err(|error| flow_error(cx, error))?;
                    prepared.debit(cx, &trial, candidate, step.duration_s(), &mut heat)?;
                    // Numerical probes retain only the native residual; the
                    // complete private fluid candidate is dropped unchanged.
                    Ok::<_, NsfThermalError<E>>(heat)
                },
                || thermal_checkpoint(cx),
            )
        };
        let difference =
            FallibleFiniteDifference::new(count, controls.finite_difference_step, |x: &[f64]| {
                let trial = evaluate(x)?;
                Ok(trial
                    .residual_j()
                    .iter()
                    .map(|r| r / controls.energy_scale_j)
                    .collect())
            })
            .map_err(difference_error)?;
        let guess = energy_guess_j
            .iter()
            .map(|value| value / controls.energy_scale_j)
            .collect();
        let initialized = NewtonKrylovState::new(&difference, guess, controls.solver);
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = initialized.map_err(NsfThermalError::Solver)?;
        let run =
            solver.run_controlled(&difference, controls.maximum_newton_iterations, &mut || {
                difference.has_failed() || cx.checkpoint().is_err()
            });
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        checkpoint(cx)?;
        if !run.report.converged {
            return Err(NsfThermalError::NotConverged(run.report));
        }
        solution = solver.x;
        nonlinear = run.report;
    }
    // The final source evaluation and actual fluid candidate also consume the
    // common work budget. This exact pair, not a replay, is later published.
    charge(
        cx,
        &evaluations,
        controls.maximum_endpoint_evaluations,
        &fluid_work,
        cell_steps,
        interval.maximum_fluid_cell_steps,
    )?;
    let energies = physical_increments(&step, &solution, controls.energy_scale_j)?;
    let fluid = fluid.into_inner();
    let final_trial = step.evaluate_with_boundary_conductances(
        &energies,
        |candidate| {
            // Move the original mutable reference into this FnOnce callback;
            // its private native trial may then live in the returned payload.
            let owner = fluid;
            let (walls, mut heat) = prepared.endpoint(cx, candidate, &resolve)?;
            // This independent native supply query is needed only on the final
            // physical acceptance, not on every Newton/Krylov probe.
            let external_heat_j = step.duration_s()
                * candidate.external_step_power_with_boundary_conductances(
                    &heat.bulk.sources_w,
                    &heat.bulk.conductances_w_k,
                    &heat.boundaries,
                    || thermal_checkpoint(cx),
                )?;
            if !external_heat_j.is_finite() {
                return Err(NsfThermalError::Invalid("external thermal heat overflowed"));
            }
            let trial = owner
                .preview_thermal_wall_interval(
                    walls,
                    interval.fluid_steps,
                    controls.maximum_fluid_cells,
                    cell_steps,
                    || cx.checkpoint().is_err(),
                )
                .map_err(|error| flow_error(cx, error))?;
            let sample = prepared.debit(cx, &trial, candidate, step.duration_s(), &mut heat)?;
            Ok::<_, NsfThermalError<E>>(BoundaryStepCoefficients {
                bulk: StepCoefficients {
                    sources_w: heat.bulk.sources_w,
                    conductances_w_k: heat.bulk.conductances_w_k,
                    payload: (
                        trial,
                        Sample {
                            heat: heat.bulk.payload,
                            exchanges: sample.0,
                            coupled_heat_j: sample.1,
                            external_heat_j,
                            interface_residual_j: sample.2,
                            maximum_cell_interface_residual_j: sample.3,
                        },
                    ),
                },
                boundaries: heat.boundaries,
            })
        },
        || thermal_checkpoint(cx),
    )?;
    let maximum_cell_residual_j = final_trial.maximum_cell_residual_j();
    let (next_thermal, (fluid_trial, sample), thermal_audit) = final_trial.accept(
        controls.cell_tolerance_j,
        controls.thermal_tolerance_j,
        || thermal_checkpoint(cx),
    )?;
    let ledger = *fluid_trial.energy_ledger();
    let mut fluid_defect = Sum::default();
    fluid_defect.add(fluid_trial.grid().total_energy_j() - fluid_trial.previous_energy_j())?;
    for input in [
        ledger.wall_heat_j,
        ledger.wall_work_j,
        ledger.open_boundary_energy_j,
        ledger.open_energy_j,
        ledger.source_energy_j,
        ledger.gravity_work_j,
        ledger.body_force_work_j,
    ] {
        fluid_defect.add(-input)?;
    }
    let fluid_energy_residual_j = fluid_defect.total()?;
    gate_energy("fluid", fluid_energy_residual_j, controls.fluid_tolerance_j)?;
    gate_energy(
        "interface",
        sample.interface_residual_j,
        controls.interface_tolerance_j,
    )?;
    gate_energy(
        "cell interface",
        sample.maximum_cell_interface_residual_j,
        controls.interface_tolerance_j,
    )?;
    // Keep the original datum-stable thermal storage increment. Independent
    // external heat was measured before the interface was added, so this gate
    // cannot silently credit that internal exchange as an external source.
    let mut joint = Sum::default();
    joint.add(fluid_energy_residual_j)?;
    joint.add(thermal_audit.stored_energy_change_j)?;
    joint.add(-sample.external_heat_j)?;
    joint.add(sample.coupled_heat_j)?;
    let total_energy_residual_j = joint.total()?;
    gate_energy(
        "combined",
        total_energy_residual_j,
        controls.total_tolerance_j,
    )?;
    checkpoint(cx)?;
    let report = NsfThermalReport {
        fluid: ledger,
        thermal: thermal_audit,
        heat: sample.heat,
        exchanges: sample.exchanges,
        coupled_heat_into_fluid_j: sample.coupled_heat_j,
        solid_external_heat_j: sample.external_heat_j,
        fluid_energy_residual_j,
        interface_energy_residual_j: sample.interface_residual_j,
        maximum_cell_interface_residual_j: sample.maximum_cell_interface_residual_j,
        total_energy_residual_j,
        maximum_cell_residual_j,
        endpoint_evaluations: evaluations.get(),
        fluid_steps: interval.fluid_steps,
        fluid_cell_steps: fluid_work.get(),
        duration_s,
        nonlinear,
    };
    // No fallible work follows either publication.
    fluid_trial.accept();
    *thermal = next_thermal;
    Ok(report)
}

impl<'a> Prepared<'a> {
    fn new<L: Lattice3<27>, E>(
        cx: &Cx<'_>,
        fluid: &NsfGrid3<L>,
        thermal: &DepositionSimulation<'_>,
        contacts: &[NsfSolidContact<'a>],
        maximum: usize,
    ) -> Result<Self, NsfThermalError<E>> {
        if contacts.is_empty() || contacts.len() > maximum {
            return Err(NsfThermalError::Invalid(
                "nonempty bounded native face mapping required",
            ));
        }
        let area = fluid.law().dx_m * fluid.law().dx_m;
        if !positive(area) {
            return Err(NsfThermalError::Invalid(
                "native fluid face area is unrepresentable",
            ));
        }
        let mut faces = BTreeSet::new();
        let mut covered: BTreeMap<CellId, Sum> = BTreeMap::new();
        for contact in contacts {
            checkpoint(cx)?;
            if contact.axis >= 3 || !faces.insert(contact.key()) {
                return Err(NsfThermalError::Invalid(
                    "fluid faces must be valid and uniquely mapped",
                ));
            }
            let cell = thermal
                .cells()
                .get(&contact.solid_cell)
                .ok_or(NsfThermalError::Invalid("coupled solid cell is inactive"))?;
            if !thermal.uses_material_law(cell.material_id, contact.caloric) {
                return Err(NsfThermalError::Invalid(
                    "coupled solid needs its original registered caloric law",
                ));
            }
            if let NsfSolidLayer::Conducting {
                half_distance_m,
                contact_resistance_m2_k_w,
                ..
            } = contact.layer
            {
                if !positive(half_distance_m) || !nonnegative(contact_resistance_m2_k_w) {
                    return Err(NsfThermalError::Invalid(
                        "solid half-layer needs positive distance and nonnegative contact resistance",
                    ));
                }
            }
            covered.entry(contact.solid_cell).or_default().add(area)?;
        }
        let mut remaining_ambient = BTreeMap::new();
        for (cell, area) in covered {
            checkpoint(cx)?;
            let remaining = thermal.exposed_area_m2(cell)? - area.total()?;
            if !nonnegative(remaining) {
                return Err(NsfThermalError::Invalid(
                    "native fluid faces exceed the solid's exposed area",
                ));
            }
            remaining_ambient.insert(cell, remaining);
        }
        let mut contacts = contacts.to_vec();
        contacts.sort_unstable_by_key(NsfSolidContact::key);
        Ok(Self {
            contacts,
            original_walls: fluid.thermal_walls().to_vec(),
            remaining_ambient,
        })
    }

    fn endpoint<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        resolve: &impl Fn(&Cx<'_>, &DepositionSimulation<'m>) -> Result<BoundaryStepCoefficients<T>, E>,
    ) -> Result<(Vec<NsfThermalWallFace>, BoundaryStepCoefficients<T>), NsfThermalError<E>> {
        checkpoint(cx)?;
        let mut heat = resolve(cx, thermal).map_err(NsfThermalError::Resolve)?;
        let mut walls: BTreeMap<_, _> = self
            .original_walls
            .iter()
            .map(|wall| ((wall.cell, wall.axis, wall.high), *wall))
            .collect();
        for (&cell, &remaining) in &self.remaining_ambient {
            checkpoint(cx)?;
            if heat.boundaries.ambient_w_k.contains_key(&cell) {
                return Err(NsfThermalError::Invalid(
                    "coupled face coverage already owns this cell's ambient conductance",
                ));
            }
            let conductance = remaining * thermal.cells()[&cell].boundary.htc_w_m2_k;
            if !nonnegative(conductance) {
                return Err(NsfThermalError::Invalid(
                    "remaining ambient conductance overflowed",
                ));
            }
            heat.boundaries.ambient_w_k.insert(cell, conductance);
        }
        for contact in &self.contacts {
            checkpoint(cx)?;
            let temperature_k = thermal.temperature_k(contact.solid_cell)?;
            let boundary = match contact.layer {
                NsfSolidLayer::Isothermal => NsfWallThermal::Isothermal(temperature_k),
                NsfSolidLayer::Conducting {
                    half_distance_m,
                    contact_resistance_m2_k_w,
                    conductivity,
                } => {
                    let specific_h = thermal.cells()[&contact.solid_cell].specific_h_j_kg;
                    let k = conductivity.conductivity_w_m_k(specific_h)?;
                    if !nonnegative(k) {
                        return Err(NsfThermalError::Invalid(
                            "original solid conductivity must be finite and nonnegative",
                        ));
                    }
                    if k == 0.0 {
                        NsfWallThermal::Adiabatic
                    } else {
                        let resistance = half_distance_m / k + contact_resistance_m2_k_w;
                        let coefficient = 1.0 / resistance;
                        if !positive(resistance) || !positive(coefficient) {
                            return Err(NsfThermalError::Invalid(
                                "solid normal resistance is unrepresentable",
                            ));
                        }
                        NsfWallThermal::Film {
                            coefficient_w_m2_k: coefficient,
                            ambient_temperature_k: temperature_k,
                        }
                    }
                }
            };
            walls.insert(
                contact.key(),
                NsfThermalWallFace {
                    cell: contact.fluid_cell,
                    axis: contact.axis,
                    high: contact.high,
                    thermal: boundary,
                },
            );
        }
        Ok((walls.into_values().collect(), heat))
    }

    fn debit<L: Lattice3<27>, T, E>(
        &self,
        cx: &Cx<'_>,
        trial: &NsfThermalWallTrial<'_, L>,
        thermal: &DepositionSimulation<'_>,
        duration_s: f64,
        heat: &mut BoundaryStepCoefficients<T>,
    ) -> Result<(Vec<NsfSolidExchange>, f64, f64, f64), NsfThermalError<E>> {
        let actual: BTreeMap<_, _> = trial
            .wall_exchange()
            .iter()
            .map(|face| ((face.cell, face.axis, face.high), *face))
            .collect();
        if actual.len() != trial.wall_exchange().len() {
            return Err(NsfThermalError::Invalid(
                "native thermal face exchange identities changed",
            ));
        }
        let mut exchanges = Vec::with_capacity(self.contacts.len());
        let mut heat_by_cell: BTreeMap<CellId, Sum> = BTreeMap::new();
        let mut coupled = Sum::default();
        for contact in &self.contacts {
            checkpoint(cx)?;
            let face = *actual.get(&contact.key()).ok_or(NsfThermalError::Invalid(
                "original coupled face heat is absent",
            ))?;
            coupled.add(face.heat_into_fluid_j)?;
            heat_by_cell
                .entry(contact.solid_cell)
                .or_default()
                .add(face.heat_into_fluid_j)?;
            exchanges.push(NsfSolidExchange {
                face,
                solid_cell: contact.solid_cell,
            });
        }
        let coupled_heat_j = coupled.total()?;
        let mut defect = Sum::default();
        defect.add(coupled_heat_j)?;
        let mut maximum_cell_interface_residual_j = 0.0_f64;
        for (cell, joules) in heat_by_cell {
            checkpoint(cx)?;
            let joules = joules.total()?;
            let original = heat.bulk.sources_w.get(&cell).copied().unwrap_or(0.0);
            let proposed = original - joules / duration_s;
            if !proposed.is_finite() {
                return Err(NsfThermalError::Invalid(
                    "combined solid interface source overflowed",
                ));
            }
            // Audit the ACTUALLY represented additional source, so a small
            // transfer lost beside either an authored or resolved large source
            // cannot hide behind another cell's opposite lost transfer.
            let authored = thermal.cells()[&cell].boundary.source_w;
            let actual_before = authored + original;
            let actual_after = authored + proposed;
            let added_j = duration_s * (actual_after - actual_before);
            defect.add(added_j)?;
            let local_defect = added_j + joules;
            if !local_defect.is_finite() {
                return Err(NsfThermalError::Invalid(
                    "cell interface heat reduction overflowed",
                ));
            }
            maximum_cell_interface_residual_j =
                maximum_cell_interface_residual_j.max(local_defect.abs());
            heat.bulk.sources_w.insert(cell, proposed);
        }
        Ok((
            exchanges,
            coupled_heat_j,
            defect.total()?,
            maximum_cell_interface_residual_j,
        ))
    }
}

fn physical_increments<E>(
    step: &BackwardEulerStep<'_, '_>,
    x: &[f64],
    scale: f64,
) -> Result<Vec<f64>, NsfThermalError<E>> {
    let physical: Vec<_> = x.iter().map(|value| value * scale).collect();
    if physical.iter().any(|value| !value.is_finite()) {
        return Err(NsfThermalError::Invalid(
            "solid energy increments overflowed",
        ));
    }
    Ok(step.resolvable_increments_j(&physical))
}

fn charge<E>(
    cx: &Cx<'_>,
    work: &Cell<usize>,
    maximum: usize,
    fluid_work: &Cell<usize>,
    cell_steps: usize,
    maximum_cell_steps: usize,
) -> Result<(), NsfThermalError<E>> {
    checkpoint(cx)?;
    if work.get() == maximum {
        return Err(NsfThermalError::WorkBudget {
            performed: work.get(),
            maximum,
        });
    }
    let next_work = fluid_work
        .get()
        .checked_add(cell_steps)
        .filter(|&next| next <= maximum_cell_steps)
        .ok_or(NsfThermalError::FluidWorkBudget {
            performed: fluid_work.get(),
            requested: cell_steps,
            maximum: maximum_cell_steps,
        })?;
    work.set(work.get() + 1);
    fluid_work.set(next_work);
    Ok(())
}

fn gate_energy<E>(
    component: &'static str,
    residual_j: f64,
    allowance_j: f64,
) -> Result<(), NsfThermalError<E>> {
    if !residual_j.is_finite() || residual_j.abs() > allowance_j {
        return Err(NsfThermalError::Energy {
            component,
            residual_j,
            allowance_j,
        });
    }
    Ok(())
}

fn thermal_checkpoint(cx: &Cx<'_>) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}
fn checkpoint<E>(cx: &Cx<'_>) -> Result<(), NsfThermalError<E>> {
    thermal_checkpoint(cx).map_err(Into::into)
}
fn flow_error<E>(cx: &Cx<'_>, error: NsfError) -> NsfThermalError<E> {
    if cx.checkpoint().is_err() {
        NsfThermalError::Cancelled
    } else {
        NsfThermalError::Fluid(error)
    }
}
fn difference_error<E>(error: FiniteDifferenceError<NsfThermalError<E>>) -> NsfThermalError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => NsfThermalError::Invalid(message),
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), NsfThermalError<E>> {
        let next = self.value + value;
        if !value.is_finite() || !next.is_finite() {
            return Err(NsfThermalError::Invalid("heat reduction overflowed"));
        }
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
        if !self.correction.is_finite() {
            return Err(NsfThermalError::Invalid(
                "heat reduction correction overflowed",
            ));
        }
        Ok(())
    }
    fn total<E>(&self) -> Result<f64, NsfThermalError<E>> {
        let result = self.value + self.correction;
        if result.is_finite() {
            Ok(result)
        } else {
            Err(NsfThermalError::Invalid("heat reduction is nonfinite"))
        }
    }
}

#[cfg(test)]
mod tests;
