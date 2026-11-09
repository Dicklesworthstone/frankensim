//! Original caloric storage coupled to resolved wall-normal film momentum.
//!
//! The native backward-Euler residual owns storage, heat transfer and energy
//! acceptance. Each thermal endpoint trial advances the SAME accepted velocity
//! history with the original Carreau source at that trial's face temperatures.
//! Actual local bulk viscous loss enters the native caloric residual once.
//! A finite wall slip length requires an explicit heat partition: the declared
//! fraction heats the original wall-adjacent fluid cell, and the complement
//! heats an explicitly named existing cell or leaves to surroundings.
//! Interface heat is distinct from actuator work.
//! Endpoint thermal coefficients are held through midpoint momentum: this
//! coupled composition is first order when temperature evolves.

use std::cell::{Cell, RefCell};
use std::collections::{BTreeMap, BTreeSet};
use std::convert::Infallible;

use fs_conduction::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, BoundaryConductances, BoundaryStepCoefficients,
    StepCoefficients,
};
use fs_conduction::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
};
use fs_math::dd::Dd;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};

use super::{
    FilmDrive, FilmStepLimits, TransientCarreauControls, TransientCarreauError,
    TransientCarreauFilm, TransientCarreauFilmSpec, TransientCarreauStep, filled,
};

pub mod adaptive;

/// Explicit receiver of the wall-interface loss not assigned to the fluid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalFilmSlipHeatReceiver {
    /// Leaves the represented thermal inventory; reported as outward heat.
    Surroundings,
    /// Deposited into this active original caloric cell, e.g. a finite wall.
    Cell(CellId),
}

/// Authored heat partition; no interfacial caloric law is inferred from b.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalFilmSlipHeat {
    /// Fraction in [0,1] deposited into film_cells()[0].
    pub fluid_fraction: f64,
    /// Receiver of the complementary fraction. A named cell must exist even
    /// when its current fraction or instantaneous interface loss is zero.
    pub remainder: ThermalFilmSlipHeatReceiver,
}

/// Original thermal work/accuracy bounds, separate from momentum controls.
#[derive(Clone, Copy, Debug)]
pub struct ThermalFilmControls {
    /// Includes native walls/solids outside the ordered fluid cells.
    pub maximum_thermal_cells: usize,
    pub maximum_newton_iterations: usize,
    /// Residuals, derivative probes and independent final evaluation.
    pub maximum_endpoint_evaluations: usize,
    /// Actual mechanical-drive callback queries across this entire call,
    /// including derivative probes and discarded adaptive comparisons. A held
    /// FilmDrive value consumes no query. Zero admits only paths with no query.
    pub maximum_drive_evaluations: usize,
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
            maximum_drive_evaluations: 4096,
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
    /// Do not include film bulk/interface heat or actuator work: this adapter
    /// installs the actual local bulk loss and declared interface fraction once.
    pub heat: BoundaryStepCoefficients<T>,
}

/// Original mechanical source evaluated at one native thermal candidate.
/// The payload may retain the source identity or a prepared source receipt;
/// the numerical owner neither clones nor reconstructs it.
#[derive(Debug)]
pub struct ThermalFilmDrive<D> {
    pub drive: FilmDrive,
    pub payload: D,
}

/// The exact mechanical sample that produced an accepted native interval.
#[derive(Debug)]
pub struct ThermalFilmDriveSample<D = ()> {
    pub query: adaptive::ThermalFilmQuery,
    /// None for a held/prescribed source. Some(endpoint_s) identifies actual
    /// endpoint caloric feedback, with midpoint_s retained separately.
    pub thermal_time_s: Option<f64>,
    pub drive: FilmDrive,
    pub payload: D,
}

/// Refused proposals change neither accepted physical state.
#[derive(Debug)]
pub enum ThermalFilmError<E> {
    Invalid(&'static str),
    Film(TransientCarreauError),
    Thermal(DepositionError),
    Source(E),
    Drive(E),
    Solver(NewtonError),
    NotConverged(NewtonReport),
    EndpointBudget,
    DriveBudget,
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
    wall_slip_heat: ThermalFilmSlipHeat,
    trial: TransientCarreauFilm,
    face_temperatures: Vec<f64>,
    cell_heat: Vec<f64>,
}

/// Physical exchanges and separate algebraic defects for one interval.
#[derive(Debug)]
pub struct ThermalFilmReport<T, D = ()> {
    pub flow: TransientCarreauStep,
    pub thermal: EnergyAudit,
    pub exchange: T,
    pub drive: ThermalFilmDriveSample<D>,
    pub nonlinear: NewtonReport,
    pub endpoint_evaluations: usize,
    pub drive_evaluations: usize,
    pub momentum_material_evaluations: u64,
    pub maximum_cell_residual_j: f64,
    /// Actual original bulk fluid dissipation inserted into the caloric owner.
    pub deposited_viscous_heat_j: f64,
    /// Total physical wall-interface loss, distinct from bulk fluid heating.
    pub wall_slip_heat_j: f64,
    pub wall_slip_heat_routing: ThermalFilmSlipHeat,
    /// Declared interface fraction deposited into film_cells()[0].
    pub deposited_fluid_slip_heat_j: f64,
    /// Complement deposited into the explicitly named original receiver cell.
    pub deposited_receiver_slip_heat_j: f64,
    /// Complement of the interface loss leaving to unmodelled surroundings.
    /// This is an explicit outward exchange, not heat in a represented wall.
    pub outward_wall_slip_heat_j: f64,
    /// Actual bulk plus routed interface heat in film_cells() order.
    pub deposited_cell_heat_j: Vec<f64>,
    /// Deposited bulk heat minus area times native bulk viscous loss.
    pub heat_partition_residual_j: f64,
    /// Deposited plus outward interface heat minus native interface loss.
    pub wall_slip_partition_residual_j: f64,
    pub kinetic_energy_change_j: f64,
    /// Signed body + moving-wall + free-surface traction work, J.
    pub mechanical_work_j: f64,
    /// Original external thermal supply BEFORE inserting bulk/interface heat.
    pub other_heat_j: f64,
    /// Original net cell heat BEFORE inserting bulk/interface heat, including internal
    /// conduction. Uses the same endpoint source and conductance overrides as
    /// the accepted thermal residual. Different mechanisms within one cell may
    /// cancel: this is not a per-boundary or per-mechanism error certificate.
    pub other_cell_heat_j: BTreeMap<CellId, f64>,
    /// Thermal storage + kinetic change - other heat - mechanical work
    /// + explicitly outward interface heat.
    /// Thermal storage is measured by the native owner relative to accepted
    /// enthalpy; no absolute caloric reference datum is subtracted here.
    pub combined_energy_residual_j: f64,
}

struct Endpoint<T, D> {
    flow: TransientCarreauStep,
    exchange: T,
    drive: ThermalFilmDriveSample<D>,
    deposited_heat_j: f64,
    wall_slip_heat_j: f64,
    deposited_fluid_slip_heat_j: f64,
    deposited_receiver_slip_heat_j: f64,
    outward_wall_slip_heat_j: f64,
    other_heat_j: f64,
    other_cell_heat_j: BTreeMap<CellId, f64>,
}

/// One aggregate counter shared by every native comparison in an adaptive call.
#[derive(Default)]
struct ThermalFilmWork {
    endpoints: Cell<usize>,
    drives: Cell<usize>,
    materials: Cell<u64>,
}

impl ThermalFilmWork {
    fn query_drive<E>(&self, maximum: usize) -> Result<(), ThermalFilmError<E>> {
        if self.drives.get() >= maximum {
            return Err(ThermalFilmError::DriveBudget);
        }
        self.drives.set(self.drives.get() + 1);
        Ok(())
    }
}

/// Forward the actual endpoint overrides to the original thermal rate owner.
/// This is only an observation of the same rates, before viscous insertion.
#[derive(Clone)]
struct EndpointRates<'a> {
    sources: &'a BTreeMap<CellId, f64>,
    links: &'a BTreeMap<u64, f64>,
    boundaries: &'a BoundaryConductances,
}
impl DepositionObserver for EndpointRates<'_> {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        check()?;
        Ok(self.sources.clone())
    }
    fn conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, DepositionError> {
        check()?;
        Ok(self.links.clone())
    }
    fn ambient_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        check()?;
        Ok(self.boundaries.ambient_w_k.clone())
    }
    fn reservoir_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        check()?;
        Ok(self.boundaries.reservoir_w_k.clone())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        Ok(())
    }
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
    /// Nonzero slip requires [`Self::new_with_wall_slip_heat`] so its heat
    /// receiver is always explicit.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        flow: TransientCarreauFilm,
        thermal: DepositionSimulation<'m>,
        film_cells: &[CellId],
        wetted_area_m2: f64,
        initial_boundary_temperatures_k: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ThermalFilmError<Infallible>> {
        if flow.wall_slip_length_m() != 0.0 {
            if cancelled() {
                return Err(DepositionError::Cancelled.into());
            }
            return Err(ThermalFilmError::Invalid(
                "wall slip requires an explicit thermal loss partition",
            ));
        }
        Self::new_with_wall_slip_heat(
            flow,
            thermal,
            film_cells,
            wetted_area_m2,
            initial_boundary_temperatures_k,
            ThermalFilmSlipHeat {
                fluid_fraction: 0.0,
                remainder: ThermalFilmSlipHeatReceiver::Surroundings,
            },
            cancelled,
        )
    }

    /// Bind the SAME native owners with an explicit wall-interface heat split.
    /// The authored fraction heats the original wall-adjacent fluid cell; the
    /// complement heats an explicitly named existing cell or is reported
    /// outward to surroundings. A finite wall uses its own original caloric
    /// material and storage. The heat callback must not deposit these losses
    /// again. The prescribed slip length remains on `flow`; neither
    /// the length nor the partition is inferred from rheology or caloric data.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_wall_slip_heat(
        mut flow: TransientCarreauFilm,
        thermal: DepositionSimulation<'m>,
        film_cells: &[CellId],
        wetted_area_m2: f64,
        initial_boundary_temperatures_k: [f64; 2],
        routing: ThermalFilmSlipHeat,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ThermalFilmError<Infallible>> {
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        if !routing.fluid_fraction.is_finite() || !(0.0..=1.0).contains(&routing.fluid_fraction) {
            return Err(ThermalFilmError::Invalid(
                "wall slip heat fraction must lie in [0,1]",
            ));
        }
        if let ThermalFilmSlipHeatReceiver::Cell(id) = routing.remainder {
            if !thermal.cells().contains_key(&id) {
                return Err(ThermalFilmError::Invalid(
                    "wall slip heat receiver must be an active original cell",
                ));
            }
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
        trial.set_wall_slip_length_m(flow.wall_slip_length_m())?;
        trial.set_face_temperatures(&temperatures, || (checkpoint.borrow_mut())())?;
        check()?;
        Ok(Self {
            flow,
            thermal,
            cells: film_cells.to_vec(),
            area: wetted_area_m2,
            wall_slip_heat: routing,
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
    /// Authored wall-interface partition and its explicit original receiver.
    pub const fn wall_slip_heat_routing(&self) -> ThermalFilmSlipHeat {
        self.wall_slip_heat
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
        self.advance_metered(
            adaptive::ThermalFilmQuery {
                start_s: self.thermal.time_s(),
                midpoint_s: self.thermal.time_s() + 0.5 * (target_s - self.thermal.time_s()),
                endpoint_s: target_s,
                segment: 0,
            },
            drive,
            limits,
            momentum,
            controls,
            energy_guess_j,
            resolve,
            cancelled,
            &ThermalFilmWork::default(),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_metered<T, E>(
        &mut self,
        query: adaptive::ThermalFilmQuery,
        drive: FilmDrive,
        limits: FilmStepLimits,
        momentum: TransientCarreauControls,
        controls: ThermalFilmControls,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>, &[CellId]) -> Result<ThermalFilmExchange<T>, E>,
        cancelled: impl FnMut() -> bool,
        work: &ThermalFilmWork,
    ) -> Result<ThermalFilmReport<T>, ThermalFilmError<E>> {
        self.advance_resolved_metered(
            query,
            limits,
            momentum,
            controls,
            energy_guess_j,
            |_, _| Ok(ThermalFilmDrive { drive, payload: () }),
            |candidate, cells, _| resolve(candidate, cells),
            cancelled,
            work,
            false,
        )
    }

    /// Resolve body acceleration, wall speed and surface traction against EACH
    /// original native caloric candidate, including derivative probes and the
    /// independent final acceptance evaluation. Every candidate restarts the
    /// same accepted velocity history. Source/heat callbacks must be repeatable
    /// and free of physical side effects; discarded queries earn no work/heat.
    ///
    /// The drive callback sees the actual endpoint caloric state and the exact
    /// interval query (segment zero for this unscheduled entry). The heat
    /// callback sees that SAME drive and its owned source payload. It must not
    /// deposit mechanical work or duplicate native viscous heat. Endpoint-held
    /// feedback retains the original first-order coupled composition, support
    /// and independent energy gates. Only the final accepted sample is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_thermal_drive<T, D, E>(
        &mut self,
        target_s: f64,
        limits: FilmStepLimits,
        momentum: TransientCarreauControls,
        controls: ThermalFilmControls,
        energy_guess_j: &[f64],
        drive: impl Fn(
            adaptive::ThermalFilmQuery,
            &DepositionSimulation<'m>,
            &[CellId],
        ) -> Result<ThermalFilmDrive<D>, E>,
        resolve: impl Fn(
            &DepositionSimulation<'m>,
            &[CellId],
            &ThermalFilmDriveSample<D>,
        ) -> Result<ThermalFilmExchange<T>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalFilmReport<T, D>, ThermalFilmError<E>> {
        let query = adaptive::ThermalFilmQuery {
            start_s: self.thermal.time_s(),
            midpoint_s: self.thermal.time_s() + 0.5 * (target_s - self.thermal.time_s()),
            endpoint_s: target_s,
            segment: 0,
        };
        self.advance_resolved_metered(
            query,
            limits,
            momentum,
            controls,
            energy_guess_j,
            |candidate, cells| drive(query, candidate, cells),
            resolve,
            cancelled,
            &ThermalFilmWork::default(),
            true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_resolved_metered<T, D, E>(
        &mut self,
        query: adaptive::ThermalFilmQuery,
        limits: FilmStepLimits,
        momentum: TransientCarreauControls,
        controls: ThermalFilmControls,
        energy_guess_j: &[f64],
        drive: impl Fn(&DepositionSimulation<'m>, &[CellId]) -> Result<ThermalFilmDrive<D>, E>,
        resolve: impl Fn(
            &DepositionSimulation<'m>,
            &[CellId],
            &ThermalFilmDriveSample<D>,
        ) -> Result<ThermalFilmExchange<T>, E>,
        cancelled: impl FnMut() -> bool,
        work: &ThermalFilmWork,
        thermal_drive: bool,
    ) -> Result<ThermalFilmReport<T, D>, ThermalFilmError<E>> {
        let target_s = query.endpoint_s;
        let opening_endpoints = work.endpoints.get();
        let opening_drives = work.drives.get();
        let opening_materials = work.materials.get();
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
                || controls.maximum_drive_evaluations > 1_000_000
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
            let evaluations = &work.endpoints;
            let material_evaluations = &work.materials;
            let cells = &self.cells;
            let area = self.area;
            let routing = self.wall_slip_heat;
            let evaluate = |x: &[f64]| -> Result<
                BackwardEulerTrial<'m, Endpoint<T, D>>,
                ThermalFilmError<E>,
            > {
                check()?;
                if evaluations.get() >= controls.maximum_endpoint_evaluations {
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
                        if thermal_drive {
                            work.query_drive(controls.maximum_drive_evaluations)?;
                        }
                        let resolved = drive(candidate, cells).map_err(ThermalFilmError::Drive)?;
                        check()?;
                        let sample = ThermalFilmDriveSample {
                            query,
                            thermal_time_s: thermal_drive.then_some(candidate.time_s()),
                            drive: resolved.drive,
                            payload: resolved.payload,
                        };
                        let exchange =
                            resolve(candidate, cells, &sample).map_err(ThermalFilmError::Source)?;
                        let mut temperatures = temperatures.borrow_mut();
                        face_temperatures(
                            candidate,
                            cells,
                            exchange.boundary_temperatures_k,
                            &mut temperatures,
                            &check,
                        )?;
                        let mut flow = trial.borrow_mut();
                        // Boundary identity follows the accepted owner even
                        // after an adaptive restore or prior rejected endpoint.
                        flow.flow.wall_slip_length = accepted.flow.wall_slip_length;
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
                        let mechanical = flow.step_with_heat(
                            dt,
                            sample.drive,
                            limits,
                            budget,
                            &mut heat,
                            || check().is_err(),
                        )?;
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
                        let rates = candidate.power_balance_observed(
                            &EndpointRates {
                                sources: &bulk.sources_w,
                                links: &bulk.conductances_w_k,
                                boundaries: &boundaries,
                            },
                            || check(),
                        )?;
                        let other = finite(dt * rates.external_power_w)?;
                        let other_cells = rates
                            .net_power_w
                            .into_iter()
                            .map(|(id, rate)| Ok((id, finite(dt * rate)?)))
                            .collect::<Result<BTreeMap<_, _>, DepositionError>>()?;
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
                        let slip = finite(area * mechanical.flow.wall_slip_heat_j_m2)?;
                        let inward = finite(routing.fluid_fraction * slip)?;
                        let remainder = finite((1.0 - routing.fluid_fraction) * slip)?;
                        let slip_watts = finite(inward / dt)?;
                        if slip < 0.0
                            || (mechanical.flow.wall_slip_heat_j_m2 != 0.0 && slip == 0.0)
                            || (slip != 0.0 && routing.fluid_fraction != 0.0 && inward == 0.0)
                            || (slip != 0.0 && routing.fluid_fraction != 1.0 && remainder == 0.0)
                            || (inward != 0.0 && slip_watts == 0.0)
                        {
                            return Err(DepositionError::Arithmetic.into());
                        }
                        // Preserve the original no-slip source map exactly.
                        if slip_watts != 0.0 {
                            let source = bulk.sources_w.entry(cells[0]).or_insert(0.0);
                            *source = finite(*source + slip_watts)?;
                        }
                        let (receiver_heat, outward) = match routing.remainder {
                            ThermalFilmSlipHeatReceiver::Surroundings => (0.0, remainder),
                            ThermalFilmSlipHeatReceiver::Cell(id) => {
                                if !candidate.cells().contains_key(&id) {
                                    return Err(ThermalFilmError::Invalid(
                                        "wall slip heat receiver left the original inventory",
                                    ));
                                }
                                let watts = finite(remainder / dt)?;
                                if remainder != 0.0 && watts == 0.0 {
                                    return Err(DepositionError::Arithmetic.into());
                                }
                                if watts != 0.0 {
                                    let source = bulk.sources_w.entry(id).or_insert(0.0);
                                    *source = finite(*source + watts)?;
                                }
                                (finite(dt * watts)?, 0.0)
                            }
                        };
                        Ok::<_, ThermalFilmError<E>>(BoundaryStepCoefficients {
                            bulk: StepCoefficients {
                                sources_w: bulk.sources_w,
                                conductances_w_k: bulk.conductances_w_k,
                                payload: Endpoint {
                                    flow: mechanical,
                                    exchange: bulk.payload,
                                    drive: sample,
                                    deposited_heat_j: finite(deposited.to_f64())?,
                                    wall_slip_heat_j: slip,
                                    deposited_fluid_slip_heat_j: finite(dt * slip_watts)?,
                                    deposited_receiver_slip_heat_j: receiver_heat,
                                    outward_wall_slip_heat_j: outward,
                                    other_heat_j: other,
                                    other_cell_heat_j: other_cells,
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
            let mechanical_work_j = finite(
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
            let slip_partition = finite(
                (Dd::from_f64(endpoint.deposited_fluid_slip_heat_j)
                    + Dd::from_f64(endpoint.deposited_receiver_slip_heat_j)
                    + Dd::from_f64(endpoint.outward_wall_slip_heat_j)
                    - Dd::from_f64(endpoint.wall_slip_heat_j))
                .to_f64(),
            )?;
            // Independent endpoint storage/work closure, not a sum of the two
            // solvers' residuals that could hide a missing heat transfer.
            let combined = finite(
                (Dd::from_f64(audit.stored_energy_change_j) + Dd::from_f64(kinetic)
                    - Dd::from_f64(endpoint.other_heat_j)
                    - Dd::from_f64(mechanical_work_j)
                    + Dd::from_f64(endpoint.outward_wall_slip_heat_j))
                .to_f64(),
            )?;
            for residual in [combined, partition, slip_partition] {
                if residual.abs() > controls.combined_energy_tolerance_j {
                    return Err(ThermalFilmError::CombinedEnergy {
                        residual_j: residual,
                        tolerance_j: controls.combined_energy_tolerance_j,
                    });
                }
            }
            let report = ThermalFilmReport {
                flow: endpoint.flow,
                thermal: audit,
                exchange: endpoint.exchange,
                drive: endpoint.drive,
                nonlinear: result.report,
                endpoint_evaluations: evaluations.get() - opening_endpoints,
                drive_evaluations: work.drives.get() - opening_drives,
                momentum_material_evaluations: material_evaluations.get() - opening_materials,
                maximum_cell_residual_j,
                deposited_viscous_heat_j: endpoint.deposited_heat_j,
                wall_slip_heat_j: endpoint.wall_slip_heat_j,
                wall_slip_heat_routing: routing,
                deposited_fluid_slip_heat_j: endpoint.deposited_fluid_slip_heat_j,
                deposited_receiver_slip_heat_j: endpoint.deposited_receiver_slip_heat_j,
                outward_wall_slip_heat_j: endpoint.outward_wall_slip_heat_j,
                deposited_cell_heat_j: heat
                    .borrow()
                    .iter()
                    .enumerate()
                    .map(|(i, q)| {
                        finite(
                            dt * (area * q / dt)
                                + if i == 0 {
                                    endpoint.deposited_fluid_slip_heat_j
                                } else {
                                    0.0
                                }
                                + if routing.remainder
                                    == ThermalFilmSlipHeatReceiver::Cell(cells[i])
                                {
                                    endpoint.deposited_receiver_slip_heat_j
                                } else {
                                    0.0
                                },
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?,
                heat_partition_residual_j: partition,
                wall_slip_partition_residual_j: slip_partition,
                kinetic_energy_change_j: kinetic,
                mechanical_work_j,
                other_heat_j: endpoint.other_heat_j,
                other_cell_heat_j: endpoint.other_cell_heat_j,
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
