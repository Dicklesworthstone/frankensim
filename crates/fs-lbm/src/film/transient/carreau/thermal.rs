//! Original caloric storage coupled to resolved wall-normal film momentum.
//!
//! The native backward-Euler residual owns storage, heat transfer and energy
//! acceptance. Each thermal endpoint trial advances the SAME accepted velocity
//! history with the original Carreau source at that trial's face temperatures.
//! Actual local viscous loss enters the native caloric residual once.
//! Endpoint thermal coefficients are held through midpoint momentum: this
//! coupled composition is first order when temperature evolves.

use std::cell::{Cell, RefCell};
use std::collections::BTreeSet;
use std::convert::Infallible;

use fs_conduction::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, BoundaryStepCoefficients, StepCoefficients,
};
use fs_conduction::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use fs_math::dd::Dd;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};

use super::{
    FilmDrive, FilmStepLimits, TransientCarreauControls, TransientCarreauError,
    TransientCarreauFilm, TransientCarreauFilmSpec, TransientCarreauStep, filled,
};

/// Original thermal work/accuracy bounds, separate from momentum controls.
#[derive(Clone, Copy, Debug)]
pub struct ThermalFilmControls {
    /// Includes native walls/solids outside the ordered fluid cells.
    pub maximum_thermal_cells: usize,
    pub maximum_newton_iterations: usize,
    /// Residuals, derivative probes and independent final evaluation.
    pub maximum_endpoint_evaluations: usize,
    /// Whole-call budget for native momentum constitutive evaluations.
    /// Temperature preparation and caloric queries are bounded separately by
    /// the endpoint/cell counts. Trials do not renew this aggregate budget.
    pub maximum_momentum_material_evaluations: u64,
    pub energy_scale_j: f64,
    /// Perturbation of energy coordinates divided by energy_scale_j.
    pub finite_difference_step: f64,
    pub cell_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub combined_energy_tolerance_j: f64,
    pub solver: NewtonKrylovConfig,
}

impl Default for ThermalFilmControls {
    fn default() -> Self {
        Self {
            maximum_thermal_cells: 65_536,
            maximum_newton_iterations: 32,
            maximum_endpoint_evaluations: 4096,
            maximum_momentum_material_evaluations: 16_777_216,
            energy_scale_j: 1.0,
            finite_difference_step: 1e-5,
            cell_energy_tolerance_j: 1e-8,
            total_energy_tolerance_j: 1e-7,
            combined_energy_tolerance_j: 1e-7,
            // Directional differences include constitutive and caloric
            // roundoff. Do not demand an analytical-Jacobian inner residual
            // below that stencil accuracy; the outer native joule checks
            // still decide whether the physical endpoint may be accepted.
            solver: NewtonKrylovConfig {
                forcing_minimum: 1e-4,
                ..NewtonKrylovConfig::default()
            },
        }
    }
}

/// Source-owned heat and actual boundary temperatures at one thermal candidate.
pub struct ThermalFilmExchange<T> {
    /// Wall and free-surface fluid-side traces [K], supplied by the SAME
    /// original interface/conduction model as heat. A finite solid's bulk
    /// temperature is not automatically its fluid-side interface temperature.
    pub boundary_temperatures_k: [f64; 2],
    /// Additional sources and endpoint link/boundary conductances.
    /// Do not include film viscous heat or actuator work: this adapter installs
    /// the actual local viscous losses exactly once.
    pub heat: BoundaryStepCoefficients<T>,
}

/// Refused proposals change neither accepted physical state.
#[derive(Debug)]
pub enum ThermalFilmError<E> {
    Invalid(&'static str),
    Film(TransientCarreauError),
    Thermal(DepositionError),
    Source(E),
    Solver(NewtonError),
    NotConverged(NewtonReport),
    EndpointBudget,
    MaterialBudget,
    CombinedEnergy { residual_j: f64, tolerance_j: f64 },
}
impl<E> From<DepositionError> for ThermalFilmError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<TransientCarreauError> for ThermalFilmError<E> {
    fn from(error: TransientCarreauError) -> Self {
        Self::Film(error)
    }
}
impl<E> From<super::TransientFilmError> for ThermalFilmError<E> {
    fn from(error: super::TransientFilmError) -> Self {
        Self::Film(error.into())
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ThermalFilmError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled thermal film: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ThermalFilmError<E> {}

/// Original caloric timeline and resolved velocity at one accepted time.
///
/// Every fluid cell has mass rho*A*dy and volume A*dy. Its start/end segment
/// spans dy along the common wall normal, in supplied wall-to-surface order.
/// Original solid/wall cells may coexist on the same thermal timeline.
/// The caller binds an admitted liquid caloric law to the original Carreau
/// fluid. Viscosity never supplies capacity, conductivity or phase occupancy.
pub struct ThermalCarreauFilm<'m> {
    flow: TransientCarreauFilm,
    thermal: DepositionSimulation<'m>,
    cells: Vec<CellId>,
    area: f64,
    trial: TransientCarreauFilm,
    face_temperatures: Vec<f64>,
    cell_heat: Vec<f64>,
}

/// Physical exchanges and separate algebraic defects for one interval.
#[derive(Debug)]
pub struct ThermalFilmReport<T> {
    pub flow: TransientCarreauStep,
    pub thermal: EnergyAudit,
    pub exchange: T,
    pub nonlinear: NewtonReport,
    pub endpoint_evaluations: usize,
    pub momentum_material_evaluations: u64,
    pub maximum_cell_residual_j: f64,
    pub deposited_viscous_heat_j: f64,
    /// Sum of deposited cell heat minus area times the native flow loss.
    pub heat_partition_residual_j: f64,
    pub kinetic_energy_change_j: f64,
    /// Signed body + moving-wall + free-surface traction work, J.
    pub mechanical_work_j: f64,
    /// Original external thermal supply BEFORE inserting viscous heat.
    pub other_heat_j: f64,
    /// Actual thermal storage + kinetic change - other heat - mechanical work.
    /// Thermal storage is measured by the native owner relative to accepted
    /// enthalpy; no absolute caloric reference datum is subtracted here.
    pub combined_energy_residual_j: f64,
}

struct Endpoint<T> {
    flow: TransientCarreauStep,
    exchange: T,
    deposited_heat_j: f64,
    other_heat_j: f64,
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn same_positive(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

fn face_temperatures(
    thermal: &DepositionSimulation<'_>,
    cells: &[CellId],
    boundary: [f64; 2],
    output: &mut [f64],
    check: &dyn Fn() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    output[0] = boundary[0];
    let mut previous = thermal.temperature_k(cells[0])?;
    for (i, &cell) in cells.iter().enumerate().skip(1) {
        if i % 128 == 0 {
            check()?;
        }
        let next = thermal.temperature_k(cell)?;
        // Original caloric T interpolated on equal contiguous cells. This
        // never averages apparent viscosity or invents a conductivity law.
        output[i] = finite(0.5 * previous + 0.5 * next)?;
        previous = next;
    }
    output[cells.len()] = boundary[1];
    check()
}

impl<'m> ThermalCarreauFilm<'m> {
    /// Bind actual native inventories, without synthesizing another material
    /// or thermal state. Supplied boundary traces admit the initial source field.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mut flow: TransientCarreauFilm,
        thermal: DepositionSimulation<'m>,
        film_cells: &[CellId],
        wetted_area_m2: f64,
        initial_boundary_temperatures_k: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ThermalFilmError<Infallible>> {
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        let n = flow.velocities_m_s().len();
        let volume = wetted_area_m2 * flow.spacing_m();
        let mass = flow.density_kg_m3() * volume;
        if film_cells.len() != n
            || !same_positive(wetted_area_m2, wetted_area_m2)
            || !same_positive(volume, volume)
            || !same_positive(mass, mass)
            || thermal.time_s() != flow.accepted_time_s()
        {
            return Err(ThermalFilmError::Invalid(
                "complete film inventory, finite area and common accepted time required",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut previous_end = None;
        let mut direction: Option<[f64; 3]> = None;
        let mut material = None;
        for &id in film_cells {
            if cancelled() {
                return Err(DepositionError::Cancelled.into());
            }
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(ThermalFilmError::Invalid("film thermal cell missing"))?;
            if !seen.insert(id)
                || !same_positive(cell.mass_kg, mass)
                || !same_positive(cell.geometry.volume_m3, volume)
                || material.is_some_and(|known| known != cell.material_id)
                || previous_end.is_some_and(|point| point != cell.geometry.start_m)
            {
                return Err(ThermalFilmError::Invalid(
                    "distinct contiguous film cells must retain rho*A*dy, A*dy and one original material",
                ));
            }
            let delta: [f64; 3] =
                std::array::from_fn(|axis| cell.geometry.end_m[axis] - cell.geometry.start_m[axis]);
            let length = fs_math::det::hypot(fs_math::det::hypot(delta[0], delta[1]), delta[2]);
            if !same_positive(length, flow.spacing_m()) {
                return Err(ThermalFilmError::Invalid(
                    "native thermal segment must span one wall-normal film cell",
                ));
            }
            let normal = delta.map(|v| v / length);
            if direction.is_some_and(|old| {
                old.iter()
                    .zip(normal)
                    .any(|(a, b)| (a - b).abs() > 128.0 * f64::EPSILON)
            }) {
                return Err(ThermalFilmError::Invalid(
                    "film thermal cells must follow one straight wall normal",
                ));
            }
            direction = Some(normal);
            previous_end = Some(cell.geometry.end_m);
            material = Some(cell.material_id);
        }
        if thermal
            .material_caloric_identity(material.expect("nonempty native film"))?
            .is_none()
        {
            return Err(ThermalFilmError::Invalid(
                "film requires the original identity-bearing caloric material",
            ));
        }
        let mut temperatures = filled(n + 1)?;
        let checkpoint = RefCell::new(&mut cancelled);
        let check = || {
            if (checkpoint.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        face_temperatures(
            &thermal,
            film_cells,
            initial_boundary_temperatures_k,
            &mut temperatures,
            &check,
        )?;
        flow.set_face_temperatures(&temperatures, || (checkpoint.borrow_mut())())?;
        let mut trial = TransientCarreauFilm::from_profile(
            *flow.material().law(),
            flow.material().temperature_k(),
            TransientCarreauFilmSpec {
                density_kg_m3: flow.density_kg_m3(),
                thickness_m: flow.thickness_m(),
                maximum_laminar_shear_rate_per_s: flow.material().maximum_shear_rate_per_s(),
            },
            flow.velocities_m_s().to_vec(),
            flow.accepted_time_s(),
        )?;
        trial.set_face_temperatures(&temperatures, || (checkpoint.borrow_mut())())?;
        check()?;
        Ok(Self {
            flow,
            thermal,
            cells: film_cells.to_vec(),
            area: wetted_area_m2,
            trial,
            face_temperatures: temperatures,
            cell_heat: filled(n)?,
        })
    }

    pub fn flow(&self) -> &TransientCarreauFilm {
        &self.flow
    }
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub fn film_cells(&self) -> &[CellId] {
        &self.cells
    }
    pub const fn wetted_area_m2(&self) -> f64 {
        self.area
    }

    /// Release both accepted states for an explicit event/geometry handoff.
    /// Rebinding checks the complete physical inventory again.
    pub fn into_parts(self) -> (TransientCarreauFilm, DepositionSimulation<'m>) {
        (self.flow, self.thermal)
    }

    /// Close a thermal interval against its own source temperatures, retained
    /// momentum and locally deposited viscous heat.
    ///
    /// Unknowns are joule increments in native stable CellId order, including
    /// any coupled solid cells. The repeatable source callback supplies original
    /// heat/interface laws and actual wall/surface fluid temperatures. Interior
    /// traces interpolate native caloric T on the fixed uniform normal mesh.
    ///
    /// Momentum uses midpoint with endpoint T held over dt. The composition is
    /// first order in time. Fixed thickness, reference mass/density and admitted
    /// liquid thermal domain are required. Streamwise advection, phase occupancy,
    /// expansion and moving interfaces are outside this locally parallel owner.
    /// Interior events retain native refusal; use into_parts for their explicit
    /// handoff. Nonlinear convergence is not a time/space accuracy certificate.
    ///
    /// Trial flow arrays are reused and ALWAYS restart the same accepted
    /// velocity/time. Native thermal candidates and Krylov retain their normal
    /// allocations. Independent final energy checks precede joint publication.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<T, E>(
        &mut self,
        target_s: f64,
        drive: FilmDrive,
        limits: FilmStepLimits,
        momentum: TransientCarreauControls,
        controls: ThermalFilmControls,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>, &[CellId]) -> Result<ThermalFilmExchange<T>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalFilmReport<T>, ThermalFilmError<E>> {
        let (next, report) = {
            let checkpoint = RefCell::new(cancelled);
            let stopped = Cell::new(false);
            let check = || {
                if stopped.get() || (checkpoint.borrow_mut())() {
                    stopped.set(true);
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            };
            check()?;
            controls
                .solver
                .validate()
                .map_err(ThermalFilmError::Solver)?;
            if !(1..=1_000_000).contains(&controls.maximum_thermal_cells)
                || self.thermal.cells().len() > controls.maximum_thermal_cells
                || !(1..=1024).contains(&controls.maximum_newton_iterations)
                || !(1..=1_000_000).contains(&controls.maximum_endpoint_evaluations)
                || controls.maximum_momentum_material_evaluations == 0
                || !same_positive(controls.energy_scale_j, controls.energy_scale_j)
                || !same_positive(
                    controls.finite_difference_step,
                    controls.finite_difference_step,
                )
                || ![
                    controls.cell_energy_tolerance_j,
                    controls.total_energy_tolerance_j,
                    controls.combined_energy_tolerance_j,
                ]
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0)
            {
                return Err(ThermalFilmError::Invalid(
                    "bounded thermal work, scales and physical energy tolerances required",
                ));
            }
            let step = BackwardEulerStep::new(&self.thermal, target_s, || check())?;
            let count = step.cell_ids().len();
            if energy_guess_j.len() != count || energy_guess_j.iter().any(|v| !v.is_finite()) {
                return Err(ThermalFilmError::Invalid(
                    "complete finite native cell-energy guess required",
                ));
            }
            let dt = step.duration_s();
            let accepted = &self.flow;
            let trial = RefCell::new(&mut self.trial);
            let temperatures = RefCell::new(&mut self.face_temperatures);
            let heat = RefCell::new(&mut self.cell_heat);
            let evaluations = Cell::new(0usize);
            let material_evaluations = Cell::new(0u64);
            let cells = &self.cells;
            let area = self.area;
            let evaluate =
                |x: &[f64]| -> Result<BackwardEulerTrial<'m, Endpoint<T>>, ThermalFilmError<E>> {
                    check()?;
                    if evaluations.get() == controls.maximum_endpoint_evaluations {
                        return Err(ThermalFilmError::EndpointBudget);
                    }
                    evaluations.set(evaluations.get() + 1);
                    let energy = x
                        .iter()
                        .map(|v| finite(v * controls.energy_scale_j))
                        .collect::<Result<Vec<_>, _>>()?;
                    let energy = step.resolvable_increments_j(&energy);
                    step.evaluate_with_boundary_conductances(
                        &energy,
                        |candidate| {
                            let exchange =
                                resolve(candidate, cells).map_err(ThermalFilmError::Source)?;
                            let mut temperatures = temperatures.borrow_mut();
                            face_temperatures(
                                candidate,
                                cells,
                                exchange.boundary_temperatures_k,
                                &mut temperatures,
                                &check,
                            )?;
                            let mut flow = trial.borrow_mut();
                            flow.set_face_temperatures(&temperatures, || check().is_err())?;
                            // Reset only accepted history; preceding thermal trials
                            // are never a starting state. Private geometry/source
                            // identities are fixed by new(); scratch is reusable.
                            for i in 0..cells.len() {
                                if i % 128 == 0 {
                                    check()?;
                                }
                                flow.flow.velocity[i] = accepted.flow.velocity[i];
                            }
                            flow.flow.time = accepted.flow.time;
                            let remaining = controls.maximum_momentum_material_evaluations
                                - material_evaluations.get();
                            if remaining == 0 {
                                return Err(ThermalFilmError::MaterialBudget);
                            }
                            let budget = TransientCarreauControls {
                                maximum_material_evaluations: momentum
                                    .maximum_material_evaluations
                                    .min(remaining),
                                ..momentum
                            };
                            let mut heat = heat.borrow_mut();
                            let mechanical =
                                flow.step_with_heat(dt, drive, limits, budget, &mut heat, || {
                                    check().is_err()
                                })?;
                            material_evaluations
                                .set(material_evaluations.get() + mechanical.material_evaluations);
                            if flow.accepted_time_s() != target_s {
                                return Err(ThermalFilmError::Invalid(
                                    "thermal and momentum endpoint clocks differ",
                                ));
                            }
                            let BoundaryStepCoefficients {
                                mut bulk,
                                boundaries,
                            } = exchange.heat;
                            let other = finite(
                                dt * candidate.external_step_power_with_boundary_conductances(
                                    &bulk.sources_w,
                                    &bulk.conductances_w_k,
                                    &boundaries,
                                    || check(),
                                )?,
                            )?;
                            let mut deposited = Dd::from_f64(0.0);
                            for (&cell, &q) in cells.iter().zip(heat.iter()) {
                                check()?;
                                let joules = finite(area * q)?;
                                let watts = finite(joules / dt)?;
                                if joules < 0.0
                                    || (q != 0.0 && joules == 0.0)
                                    || (joules != 0.0 && watts == 0.0)
                                {
                                    return Err(DepositionError::Arithmetic.into());
                                }
                                deposited = deposited + Dd::from_f64(finite(dt * watts)?);
                                let source = bulk.sources_w.entry(cell).or_insert(0.0);
                                *source = finite(*source + watts)?;
                            }
                            Ok::<_, ThermalFilmError<E>>(BoundaryStepCoefficients {
                                bulk: StepCoefficients {
                                    sources_w: bulk.sources_w,
                                    conductances_w_k: bulk.conductances_w_k,
                                    payload: Endpoint {
                                        flow: mechanical,
                                        exchange: bulk.payload,
                                        deposited_heat_j: finite(deposited.to_f64())?,
                                        other_heat_j: other,
                                    },
                                },
                                boundaries,
                            })
                        },
                        || check(),
                    )
                };
            let difference = FallibleFiniteDifference::new(
                count,
                controls.finite_difference_step,
                |x: &[f64]| -> Result<Vec<f64>, ThermalFilmError<E>> {
                    let result = evaluate(x)?;
                    Ok(result
                        .residual_j()
                        .iter()
                        .map(|v| v / controls.energy_scale_j)
                        .collect())
                },
            )
            .map_err(difference_error)?;
            let initial = energy_guess_j
                .iter()
                .map(|v| v / controls.energy_scale_j)
                .collect();
            let initialized = NewtonKrylovState::new(&difference, initial, controls.solver);
            if let Some(error) = difference.take_failure() {
                return Err(difference_error(error));
            }
            let mut solver = initialized.map_err(ThermalFilmError::Solver)?;
            let result =
                solver.run_controlled(&difference, controls.maximum_newton_iterations, &mut || {
                    difference.has_failed() || check().is_err()
                });
            if let Some(error) = difference.take_failure() {
                return Err(difference_error(error));
            }
            check()?;
            if !result.report.converged {
                return Err(ThermalFilmError::NotConverged(result.report));
            }
            let final_trial = evaluate(&solver.x)?;
            let maximum_cell_residual_j = final_trial.maximum_cell_residual_j();
            let (next, endpoint, audit) = final_trial.accept(
                controls.cell_energy_tolerance_j,
                controls.total_energy_tolerance_j,
                || check(),
            )?;
            let native = endpoint.flow.flow;
            let kinetic = finite(area * native.kinetic_energy_change_j_m2)?;
            let work = finite(
                (Dd::from_f64(native.body_work_j_m2)
                    + Dd::from_f64(native.wall_work_j_m2)
                    + Dd::from_f64(native.surface_work_j_m2))
                .to_f64()
                    * area,
            )?;
            let partition = finite(
                (Dd::from_f64(endpoint.deposited_heat_j)
                    - Dd::from_f64(area) * Dd::from_f64(native.viscous_heat_j_m2))
                .to_f64(),
            )?;
            // Independent endpoint storage/work closure, not a sum of the two
            // solvers' residuals that could hide a missing heat transfer.
            let combined = finite(
                (Dd::from_f64(audit.stored_energy_change_j) + Dd::from_f64(kinetic)
                    - Dd::from_f64(endpoint.other_heat_j)
                    - Dd::from_f64(work))
                .to_f64(),
            )?;
            if combined.abs() > controls.combined_energy_tolerance_j {
                return Err(ThermalFilmError::CombinedEnergy {
                    residual_j: combined,
                    tolerance_j: controls.combined_energy_tolerance_j,
                });
            }
            let report = ThermalFilmReport {
                flow: endpoint.flow,
                thermal: audit,
                exchange: endpoint.exchange,
                nonlinear: result.report,
                endpoint_evaluations: evaluations.get(),
                momentum_material_evaluations: material_evaluations.get(),
                maximum_cell_residual_j,
                deposited_viscous_heat_j: endpoint.deposited_heat_j,
                heat_partition_residual_j: partition,
                kinetic_energy_change_j: kinetic,
                mechanical_work_j: work,
                other_heat_j: endpoint.other_heat_j,
                combined_energy_residual_j: combined,
            };
            check()?;
            (next, report)
        };
        core::mem::swap(&mut self.flow, &mut self.trial);
        self.thermal = next;
        Ok(report)
    }
}

fn difference_error<E>(error: FiniteDifferenceError<ThermalFilmError<E>>) -> ThermalFilmError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => ThermalFilmError::Invalid(message),
    }
}

#[cfg(test)]
mod tests;
