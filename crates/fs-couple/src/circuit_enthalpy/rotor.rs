//! Source-bound passive rotor damping on the original finite enthalpy residual.
//!
//! Explicit scalar interface channels add c(T) b b^T to the original rotor C.
//! Their physical slip is b.q_dot and their native midpoint heat is h*c*slip².
//! Arbitrary original C and skew K remain independently evaluated signed work;
//! neither positive net damping nor a numerical energy defect identifies heat.
//! The original symmetric K supplies elastic storage. M/K/G, spin and maps stay
//! fixed, so refreshing a purely viscous coefficient introduces no parameter work.
//!
//! Original material cards, caloric laws, conduction and thermal events retain
//! their owners. Endpoint-temperature coupling is first order in general;
//! refinement of the declared timestep and model bandwidth remains necessary.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate, finite,
    material::{
        ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
    },
};
use fs_matdb::InterfaceSystemCard;
use fs_material::state_point::{MaterialStatePointError, ResolvedLinearInterface};
use fs_math::dd::Dd;
use fs_mbd::rotordynamics::{
    RotorBearingModel, RotorError,
    transient::{RotorTransientState, RotorTransientStep},
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

/// Hard bound on explicitly authored physical damping channels.
pub const MAX_MATERIAL_ROTOR_DAMPERS: usize = 256;

/// One original pure-viscous interface and its physical mechanical/thermal maps.
///
/// The two cells identify the original ordered surface materials. The slip row
/// maps generalized rates to m/s: angular entries therefore have lever-arm units.
/// Its transpose maps force back to the same generalized coordinates. For a
/// relative motion include BOTH sides in that row; a fixed support is an authored
/// boundary. Geometry, area, preload and other source context do not evolve here.
#[derive(Clone, Copy)]
pub struct MaterialRotorDamper<'a> {
    pub interface: &'a InterfaceSystemCard,
    pub source: &'a ResolvedLinearInterface,
    pub slip_row: &'a [f64],
    pub cells: [CellId; 2],
    /// Explicit side whose original caloric temperature supplies the source T.
    pub temperature_side: usize,
    /// Side zero receives this fraction; side one receives the exact remainder.
    /// This is authored heat partition, not an inferred effusivity/flash model.
    pub first_heat_fraction: f64,
}

/// Fixed original rotor and the SAME material/caloric network receiving heat.
#[derive(Clone, Copy)]
pub struct RotorThermalNetwork<'a> {
    /// Fixed baseline M/C/K/G. Authored channel damping is ADDITIONAL to this C;
    /// do not include those same channels a second time in the baseline.
    pub baseline: &'a RotorBearingModel,
    pub spin_rad_s: f64,
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub dampers: &'a [MaterialRotorDamper<'a>],
    /// Every source must support [0, this positive frequency] in rad/s. The
    /// caller qualifies this physical response band; no automatic spectral or
    /// broadband validity claim follows from a small numerical residual.
    pub qualified_angular_frequency_rad_s: f64,
}

/// Original nonlinear/energy gates plus an explicit physical time-resolution cap.
#[derive(Clone, Copy, Debug)]
pub struct RotorThermalControls {
    pub solve: CircuitEnthalpySolve,
    pub relative_motion_tolerance: f64,
    pub maximum_step_s: f64,
}

/// Accepted motion, source states and immutable physical binding.
#[derive(Clone)]
pub struct MaterialRotorThermalState<'a> {
    network: &'a RotorThermalNetwork<'a>,
    time_s: f64,
    state: RotorTransientState,
    sources: Arc<[ResolvedLinearInterface]>,
}
impl MaterialRotorThermalState<'_> {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn state(&self) -> &RotorTransientState {
        &self.state
    }
    pub fn sources(&self) -> &[ResolvedLinearInterface] {
        &self.sources
    }
}

/// Independently evaluated native motion, source heat and original operator work.
#[derive(Clone, Debug)]
pub struct RotorThermalStep {
    pub motion: RotorTransientStep,
    /// One h*c(T)*(b.v_mid)^2 in the original channel order [J].
    pub damper_heat_j: Vec<f64>,
    /// Direct original C contraction; active supply may make it negative.
    pub baseline_velocity_work_j: f64,
    /// Direct original K_skew contraction, excluded from elastic storage.
    pub nonconservative_stiffness_work_j: f64,
    pub elastic_energy_change_j: f64,
    /// Direct original Omega G contraction; zero in exact arithmetic.
    pub gyroscopic_work_j: f64,
    /// Applied force work minus original baseline C/K_skew/gyro exchanges.
    pub supplied_to_storage_and_heat_j: f64,
    pub mechanical_energy_residual_j: f64,
}

#[derive(Clone, Debug)]
pub struct RotorThermalSample {
    pub network: MaterialNetworkSample,
    /// Complete original endpoint interface bundles, shared with accepted state.
    pub dampers: Arc<[ResolvedLinearInterface]>,
}

/// The shared account's electrical field is the original rotor step and its
/// deposited_joule_j field is the actual deposited passive-interface heat.
pub type RotorThermalReport = CircuitEnthalpyReport<RotorThermalSample, RotorThermalStep>;

/// Both accepted histories and their material samples publish together.
#[must_use]
pub struct RotorThermalCandidate<'m, 'a> {
    thermal: DepositionSimulation<'m>,
    rotor: MaterialRotorThermalState<'a>,
    report: RotorThermalReport,
}
impl<'m, 'a> RotorThermalCandidate<'m, 'a> {
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub fn rotor(&self) -> &MaterialRotorThermalState<'a> {
        &self.rotor
    }
    pub fn report(&self) -> &RotorThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        MaterialRotorThermalState<'a>,
        RotorThermalReport,
    ) {
        (self.thermal, self.rotor, self.report)
    }
}

#[derive(Debug)]
pub enum RotorThermalError {
    Material(MaterialNetworkError),
    Interface(MaterialStatePointError),
    Motion(RotorError),
    Energy { residual_j: f64, allowance_j: f64 },
}
impl core::fmt::Display for RotorThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material rotor thermal coupling refused: {self:?}")
    }
}
impl std::error::Error for RotorThermalError {}
type Error = CircuitEnthalpyError<RotorThermalError>;
type Result<T> = core::result::Result<T, Error>;

fn invalid(reason: &'static str) -> Error {
    CircuitEnthalpyError::Invalid(reason)
}
fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(RotorThermalError::Material(other)),
    }
}
fn interface(error: MaterialStatePointError) -> Error {
    CircuitEnthalpyError::Source(RotorThermalError::Interface(error))
}
fn motion(error: RotorError) -> Error {
    match error {
        RotorError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(RotorThermalError::Motion(other)),
    }
}
fn check(cancel: &RefCell<impl FnMut() -> bool>) -> Result<()> {
    if (*cancel.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn scalar(value: Dd) -> Result<f64> {
    finite(value.hi)?;
    finite(value.lo)?;
    Ok(finite(value.to_f64())?)
}
fn balance(value: Dd, allowance_j: f64) -> Result<()> {
    let residual_j = scalar(value)?;
    if residual_j.abs() > allowance_j {
        Err(CircuitEnthalpyError::Source(RotorThermalError::Energy {
            residual_j,
            allowance_j,
        }))
    } else {
        Ok(())
    }
}

impl<'a> RotorThermalNetwork<'a> {
    fn material_network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            links: self.links,
            resistors: &[],
        }
    }

    /// Bind original finite motion and material sources without advancing either.
    /// The borrowed immutable network prevents silent changes to M/K/G, spin,
    /// channel geometry, material assignment or heat routing on continuation.
    pub fn bind(
        &'a self,
        thermal: &DepositionSimulation<'_>,
        state: RotorTransientState,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialRotorThermalState<'a>> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let n = self.baseline.dofs();
        if !self.spin_rad_s.is_finite()
            || !self.qualified_angular_frequency_rad_s.is_finite()
            || self.qualified_angular_frequency_rad_s <= 0.0
            || self.dampers.is_empty()
            || self.dampers.len() > MAX_MATERIAL_ROTOR_DAMPERS
            || state.displacement.len() != n
            || state.velocity.len() != n
            || state
                .displacement
                .iter()
                .chain(&state.velocity)
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "finite original rotor state, spin and qualified source band required",
            ));
        }
        for damper in self.dampers {
            check(&cancel)?;
            if damper.slip_row.len() != n
                || damper.slip_row.iter().any(|x| !x.is_finite())
                || damper.slip_row.iter().all(|x| *x == 0.0)
                || damper.temperature_side > 1
                || damper.cells[0] == damper.cells[1]
                || !damper.first_heat_fraction.is_finite()
                || !(0.0..=1.0).contains(&damper.first_heat_fraction)
            {
                return Err(invalid(
                    "finite nonzero physical slip map and explicit two-cell heat binding required",
                ));
            }
        }
        let network = self.material_network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let originals: Vec<_> = self.dampers.iter().map(|d| (*d.source).clone()).collect();
        let sources = self.sources_at(thermal, &originals, &cancel)?;
        check(&cancel)?;
        Ok(MaterialRotorThermalState {
            network: self,
            time_s: thermal.time_s(),
            state,
            sources,
        })
    }

    fn sources_at(
        &self,
        thermal: &DepositionSimulation<'_>,
        originals: &[ResolvedLinearInterface],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Arc<[ResolvedLinearInterface]>> {
        if originals.len() != self.dampers.len() {
            return Err(invalid("one original source per rotor damper required"));
        }
        let mut sources = Vec::with_capacity(originals.len());
        for (damper, original) in self.dampers.iter().zip(originals) {
            check(cancel)?;
            let first = self
                .cells
                .iter()
                .find(|c| c.cell == damper.cells[0])
                .ok_or_else(|| invalid("inactive first rotor damping material receiver"))?;
            let second = self
                .cells
                .iter()
                .find(|c| c.cell == damper.cells[1])
                .ok_or_else(|| invalid("inactive second rotor damping material receiver"))?;
            // MaterialNetwork separately admits actual caloric/card identity.
            original
                .require_surface_materials(first.point.material.id(), second.point.material.id())
                .map_err(interface)?;
            for id in damper.cells {
                if !thermal.cells().contains_key(&id) {
                    return Err(invalid("rotor damping receivers must remain active"));
                }
            }
            let point = original.resolved().query();
            if !point.axes().contains_key("T") {
                return Err(invalid(
                    "rotor interface source requires its original absolute T coordinate",
                ));
            }
            let temperature = thermal.temperature_k(damper.cells[damper.temperature_side])?;
            let point = point
                .with("T", temperature)
                .map_err(|_| invalid("invalid rotor interface temperature coordinate"))?;
            let source = original
                .at_query_point(damper.interface, &point)
                .map_err(interface)?;
            if source.stiffness_n_m() != 0.0
                || source.damping_n_s_m() <= 0.0
                || source.constant_loss_reference_omega_rad_s().is_some()
            {
                return Err(invalid(
                    "rotor heating requires pure sourced viscosity without elastic or constant-loss anchors",
                ));
            }
            source
                .require_angular_frequency_rad_s(0.0)
                .map_err(interface)?;
            source
                .require_angular_frequency_rad_s(self.qualified_angular_frequency_rad_s)
                .map_err(interface)?;
            sources.push(source);
        }
        Ok(sources.into())
    }
}

impl<'a> MaterialRotorThermalState<'a> {
    /// Solve original rotor motion and finite heat at one common thermal endpoint.
    ///
    /// Every source coefficient is resolved at that proposed endpoint, and every
    /// native mechanical trial restarts the SAME accepted state. Existing
    /// MaterialNetwork conduction and backward-Euler/Newton/event ownership are
    /// reused. Initial guesses are cell-energy increments in sorted CellId order.
    /// No lagged-coefficient fallback or post-acceptance heat patch is performed.
    ///
    /// New coefficient tuples cost O(n³ + channels*n²) native preparation/work.
    /// One exact-coefficient motion result is cached in O(channels+n) storage;
    /// unchanged tuples reuse it while source resolution and work gates stay fresh.
    /// Every trial, including reuse, consumes the original circuit-trial budget.
    /// Source-domain, event, energy or cancellation refusal leaves inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        midpoint_force: &[f64],
        controls: RotorThermalControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorThermalCandidate<'m, 'a>> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        admit_driver_unknowns(thermal, controls.solve, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if self.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || !controls.maximum_step_s.is_finite()
            || controls.maximum_step_s <= 0.0
            || dt > controls.maximum_step_s
            || !controls.relative_motion_tolerance.is_finite()
            || controls.relative_motion_tolerance <= 0.0
            || controls.relative_motion_tolerance >= 1.0
            || midpoint_force.len() != self.state.velocity.len()
            || midpoint_force.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "synchronized clocks, bounded positive interval and finite rotor solve/load required",
            ));
        }
        let network = self.network.material_network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let original = self.network.sources_at(thermal, &self.sources, &cancel)?;
        if original.as_ref() != self.sources.as_ref() {
            return Err(invalid(
                "accepted rotor source and original caloric temperatures differ",
            ));
        }
        // All nonthermal inputs are held for this call. Reuse only exact ordered
        // coefficient tuples; material/source resolution and work gates stay fresh.
        let motion_cache: RefCell<Option<(Vec<u64>, RotorTransientStep)>> = RefCell::new(None);
        let (thermal, rotor, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls.solve,
            initial_energy_increments_j,
            |trial, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let sources = self.network.sources_at(trial, &self.sources, &cancel)?;
                let key: Vec<_> = sources
                    .iter()
                    .map(|s| s.damping_n_s_m().to_bits())
                    .collect();
                let reused = {
                    let cached = motion_cache.borrow();
                    cached
                        .as_ref()
                        .filter(|(stored, _)| stored == &key)
                        .map(|(_, step)| step.clone())
                };
                let step = if let Some(step) = reused {
                    step
                } else {
                    let baseline = self.network.baseline;
                    let n = baseline.dofs();
                    let mut damping = Vec::with_capacity(n * n);
                    for i in 0..n {
                        check(&cancel)?;
                        for j in 0..n {
                            let mut entry = Dd::from_f64(baseline.damping()[i * n + j]);
                            for (damper, source) in self.network.dampers.iter().zip(sources.iter())
                            {
                                entry = entry
                                    + Dd::from_f64(source.damping_n_s_m())
                                        * Dd::from_f64(damper.slip_row[i])
                                        * Dd::from_f64(damper.slip_row[j]);
                            }
                            damping.push(scalar(entry)?);
                        }
                    }
                    let model = RotorBearingModel::new(
                        n,
                        baseline.mass(),
                        &damping,
                        baseline.stiffness(),
                        baseline.gyro_per_rad_s(),
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(motion)?;
                    let operating = model
                        .at_speed(self.network.spin_rad_s, || (*cancel.borrow_mut())())
                        .map_err(motion)?;
                    let native = operating
                        .midpoint_stepper(dt, controls.relative_motion_tolerance, || {
                            (*cancel.borrow_mut())()
                        })
                        .map_err(motion)?;
                    let step = native
                        .step(&self.state, midpoint_force, || (*cancel.borrow_mut())())
                        .map_err(motion)?;
                    *motion_cache.borrow_mut() = Some((key, step.clone()));
                    step
                };
                let driver = self.physical_account(
                    step,
                    Arc::clone(&sources),
                    target_s,
                    dt,
                    controls.solve.total_energy_tolerance_j,
                    &cancel,
                )?;
                Ok(DrivenCoefficients {
                    circuit: driver,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: RotorThermalSample {
                        network: coefficients.source,
                        dampers: sources,
                    },
                    circuit_trials: 1,
                })
            },
            |candidate, duration, checkpoint| {
                let mut heat_by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                for (damper, &heat) in self
                    .network
                    .dampers
                    .iter()
                    .zip(&candidate.step.damper_heat_j)
                {
                    checkpoint()?;
                    let first = Dd::from_f64(heat) * Dd::from_f64(damper.first_heat_fraction);
                    for (cell, value) in [
                        (damper.cells[0], first),
                        (damper.cells[1], Dd::from_f64(heat) - first),
                    ] {
                        let sum = heat_by_cell.entry(cell).or_insert(Dd::ZERO);
                        *sum = *sum + value;
                    }
                }
                let mut powers = BTreeMap::new();
                for (cell, heat) in heat_by_cell {
                    checkpoint()?;
                    let heat = scalar(heat)?;
                    let power = finite(heat / duration)?;
                    if heat > 0.0 && power == 0.0 {
                        return Err(invalid(
                            "rotor heat is unrepresentable on the thermal interval",
                        ));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(RotorThermalCandidate {
            thermal,
            rotor,
            report,
        })
    }

    fn physical_account(
        &self,
        step: RotorTransientStep,
        sources: Arc<[ResolvedLinearInterface]>,
        target_s: f64,
        dt: f64,
        tolerance_j: f64,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<RotorDriver<'a>> {
        let baseline = self.network.baseline;
        let n = baseline.dofs();
        let half = Dd::from_f64(0.5);
        let h = Dd::from_f64(dt);
        let velocity: Vec<_> = self
            .state
            .velocity
            .iter()
            .zip(&step.state_after.velocity)
            .map(|(&a, &b)| half * (Dd::from_f64(a) + Dd::from_f64(b)))
            .collect();
        let position: Vec<_> = self
            .state
            .displacement
            .iter()
            .zip(&step.state_after.displacement)
            .map(|(&a, &b)| half * (Dd::from_f64(a) + Dd::from_f64(b)))
            .collect();
        let mut heat = Dd::ZERO;
        let mut damper_heat_j = Vec::with_capacity(sources.len());
        for (damper, source) in self.network.dampers.iter().zip(sources.iter()) {
            check(cancel)?;
            let mut slip = Dd::ZERO;
            for (&b, &v) in damper.slip_row.iter().zip(&velocity) {
                slip = slip + Dd::from_f64(b) * v;
            }
            let physical = h * Dd::from_f64(source.damping_n_s_m()) * slip * slip;
            let physical_j = scalar(physical)?;
            if physical_j < 0.0 {
                return Err(invalid("passive rotor source produced negative heat"));
            }
            heat = heat + physical;
            damper_heat_j.push(physical_j);
        }
        let mut original_c = Dd::ZERO;
        let mut skew_k = Dd::ZERO;
        let mut elastic = Dd::ZERO;
        let mut gyro = Dd::ZERO;
        for i in 0..n {
            check(cancel)?;
            let dq = Dd::from_f64(step.state_after.displacement[i])
                - Dd::from_f64(self.state.displacement[i]);
            for j in 0..n {
                let k = Dd::from_f64(baseline.stiffness()[i * n + j]);
                let kt = Dd::from_f64(baseline.stiffness()[j * n + i]);
                original_c = original_c
                    + h * velocity[i] * Dd::from_f64(baseline.damping()[i * n + j]) * velocity[j];
                skew_k = skew_k + h * velocity[i] * half * (k - kt) * position[j];
                elastic = elastic + dq * half * (k + kt) * position[j];
                gyro = gyro
                    + h * Dd::from_f64(self.network.spin_rad_s)
                        * velocity[i]
                        * Dd::from_f64(baseline.gyro_per_rad_s()[i * n + j])
                        * velocity[j];
            }
        }
        let stored = Dd::from_f64(step.kinetic_energy_change_j) + elastic;
        let supplied = Dd::from_f64(step.input_work_j) - original_c - skew_k - gyro;
        let residual = stored + heat - supplied;
        // These compare independent physical contractions, never allocate heat
        // by subtracting an energy residual or an arbitrary total damping loss.
        for defect in [
            Dd::from_f64(step.velocity_force_work_j) - original_c - heat,
            Dd::from_f64(step.stiffness_force_work_j) - elastic - skew_k,
            Dd::from_f64(step.gyroscopic_work_j) - gyro,
            Dd::from_f64(step.work_energy_residual_j),
            residual,
        ] {
            balance(defect, tolerance_j)?;
        }
        let next = MaterialRotorThermalState {
            network: self.network,
            time_s: target_s,
            state: step.state_after.clone(),
            sources,
        };
        let work = RotorThermalStep {
            motion: step,
            damper_heat_j,
            baseline_velocity_work_j: scalar(original_c)?,
            nonconservative_stiffness_work_j: scalar(skew_k)?,
            elastic_energy_change_j: scalar(elastic)?,
            gyroscopic_work_j: scalar(gyro)?,
            supplied_to_storage_and_heat_j: scalar(supplied)?,
            mechanical_energy_residual_j: scalar(residual)?,
        };
        Ok(RotorDriver {
            rotor: next,
            step: work,
            heat_j: scalar(heat)?,
            stored,
            supplied,
        })
    }
}

struct RotorDriver<'a> {
    rotor: MaterialRotorThermalState<'a>,
    step: RotorThermalStep,
    heat_j: f64,
    stored: Dd,
    supplied: Dd,
}
impl<'a> ThermalDriver for RotorDriver<'a> {
    type Model = MaterialRotorThermalState<'a>;
    type Step = RotorThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.stored.to_f64(),
            supplied_j: self.supplied.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        self.stored - self.supplied
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.rotor, self.step)
    }
}
