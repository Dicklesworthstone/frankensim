//! Quasi-steady mixed streams coupled to one finite original material wall.
//!
//! Fluid inventories are negligible on this time scale. Wall storage advances
//! by backward Euler; every endpoint port retains its original phase PH owner.
//! The stored wall coordinate and its mechanical interpretation remain explicit.

use super::{
    MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS, MixtureHeatExchangerError,
    MixtureHeatExchangerInlet, positive,
    wall_stream::{MixtureWallStreamControl, MixtureWallStreamStep, PreparedMixtureWallStream},
};
use crate::chemistry::ThermalConductance;
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, SharedWallStorage, WallContact, WallDrag, WallExchangeError, WallThermal,
};
use fs_material::thermochemical_gas::nonideal::flash::enthalpy::MAX_ENTHALPY_FLASH_EVALUATIONS;
use fs_math::eft::two_sum;
use fs_qty::{Pressure, Temperature, Time};
use fs_thermal::enthalpy::ThermalStorageBasis;

/// Maximum original wall-state trials in one advance.
pub const MAX_FINITE_WALL_EVALUATIONS: usize = 132;
/// Maximum scalar stream-duty trials across both ports and all wall trials.
pub const MAX_FINITE_WALL_DUTY_EVALUATIONS: usize =
    2 * MAX_FINITE_WALL_EVALUATIONS * MAX_HEAT_EXCHANGER_DUTY_EVALUATIONS;
/// Maximum original PT selections, including the two inlet preparations.
pub const MAX_FINITE_WALL_PT_EVALUATIONS: usize =
    2 + MAX_FINITE_WALL_DUTY_EVALUATIONS * MAX_ENTHALPY_FLASH_EVALUATIONS;

/// Explicit physical work convention for the SAME original wall storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WallStreamMechanicalConstraint {
    /// The original source stores U at fixed volume; wall pressure work is zero.
    IsochoricInternalEnergy,
    /// Caller-declared stationary-solid approximation: H is retained as H,
    /// while expansion, p dV and the H/U increment difference are neglected.
    NegligibleExpansionEnthalpy,
    /// The source stores H at this caller-declared constant pressure and must
    /// supply its original density at both energies. Report p Delta V and
    /// Delta U = Delta H - p Delta V; no reference geometry is mutated here.
    /// The caller must use the pressure at which that source law was authored.
    IsobaricEnthalpy { pressure: Pressure },
}

/// A supplied steady feed and independent wall interface; no wall mass is copied.
#[derive(Clone, Copy, Debug)]
pub struct MixtureWallPort<'m> {
    /// Immutable composition, original calorics, fixed pressure and mass flow.
    pub inlet: MixtureHeatExchangerInlet<'m>,
    /// Complete nonnegative wall-to-chamber conductance.
    pub conductance: ThermalConductance,
    /// Independent source policy, duty bracket and per-port work allowances.
    pub controls: MixtureWallStreamControl,
}

/// Physical tolerances and whole-advance work allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteWallStreamControl {
    /// Increasing finite signed heat into the wall [J], containing zero.
    pub wall_heat_bracket_j: [f64; 2],
    /// Absolute wall balance and integrated-contact allowance [J].
    pub absolute_energy_tolerance_j: f64,
    /// Relative allowance scaled by wall and individual stream heat changes.
    pub relative_energy_tolerance: f64,
    /// Numerical allowance for negative total wall-plus-stream entropy [J/K].
    pub entropy_tolerance_j_k: f64,
    /// Original wall-state trials, including the initial state.
    pub maximum_wall_evaluations: usize,
    /// Total scalar duty queries across all ports and wall trials.
    pub maximum_duty_evaluations: usize,
    /// Total native PT work, charging each immutable inlet preparation once.
    pub maximum_pt_evaluations: usize,
}
impl FiniteWallStreamControl {
    /// Declare the supported wall heat interval and finite default allowances.
    #[must_use]
    pub fn new(wall_heat_bracket_j: [f64; 2]) -> Self {
        Self {
            wall_heat_bracket_j,
            absolute_energy_tolerance_j: 1e-5,
            relative_energy_tolerance: 1e-9,
            entropy_tolerance_j_k: 1e-7,
            maximum_wall_evaluations: MAX_FINITE_WALL_EVALUATIONS,
            maximum_duty_evaluations: 8192,
            maximum_pt_evaluations: 65_536,
        }
    }

    /// Validate controls without changing any physical owner.
    pub fn validate(self) -> Result<(), FiniteWallStreamError> {
        let [lo, hi] = self.wall_heat_bracket_j;
        if !lo.is_finite() || !hi.is_finite() || lo > 0.0 || hi < 0.0 || lo >= hi
            || !positive(self.absolute_energy_tolerance_j)
            || !self.relative_energy_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_energy_tolerance)
            || !self.entropy_tolerance_j_k.is_finite() || self.entropy_tolerance_j_k < 0.0
            || !(1..=MAX_FINITE_WALL_EVALUATIONS).contains(&self.maximum_wall_evaluations)
            || !(1..=MAX_FINITE_WALL_DUTY_EVALUATIONS).contains(&self.maximum_duty_evaluations)
            || !(1..=MAX_FINITE_WALL_PT_EVALUATIONS).contains(&self.maximum_pt_evaluations)
        {
            return Err(FiniteWallStreamError::Invalid("finite wall-stream controls"));
        }
        Ok(())
    }

    fn tolerance(self, scale: f64) -> f64 {
        self.absolute_energy_tolerance_j + self.relative_energy_tolerance * scale
    }
}

/// A refusal publishes neither a changed wall nor a partial outlet collection.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteWallStreamError {
    /// Invalid input or unrepresentable physical coordinate.
    Invalid(&'static str),
    /// A required original entropy, volume or physical storage convention is absent.
    Unsupported(&'static str),
    /// Original material-wall admission or caloric refusal.
    Wall(WallExchangeError),
    /// Original stream refusal. Completed counts precede this port evaluation;
    /// its nested source error retains local attempted work where available.
    Stream {
        port: usize,
        completed_duty_evaluations: usize,
        completed_pt_evaluations: usize,
        source: MixtureHeatExchangerError,
    },
    /// Whole-advance allowance exhausted before another bounded query.
    WorkLimit { wall_evaluations: usize, duty_evaluations: usize, pt_evaluations: usize },
    /// The admitted original wall heat interval does not bracket the balance.
    NotBracketed { wall_heat_bracket_j: [f64; 2], residual_j: [f64; 2] },
    /// Floating-point heat resolution exhausted before physical closure.
    NoConvergence,
    /// Actual original wall storage or native stream heat failed closure.
    Balance {
        energy_residual_j: f64,
        wall_target_residual_j: f64,
        maximum_contact_residual_j: f64,
        mechanical_energy_residual_j: f64,
    },
    /// Original combined entropy change violated the declared allowance.
    Entropy { production_j_k: f64, tolerance_j_k: f64 },
    /// Cancellation before atomic wall publication.
    Cancelled,
}
impl core::fmt::Display for FiniteWallStreamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite wall with original mixture streams: {self:?}")
    }
}
impl std::error::Error for FiniteWallStreamError {}
type Error = FiniteWallStreamError;

/// One accepted wall endpoint and instantaneous outlet collection.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWallStreamStep {
    /// Backward-Euler source interval; zero retains instantaneous flow semantics.
    pub elapsed: Time,
    /// End-of-interval quasi-steady outlet states in the supplied port order.
    /// Their counters cover only the accepted wall trial, excluding preparation.
    pub streams: Vec<MixtureWallStreamStep>,
    /// Original storage coordinate; H is never relabelled as rigid-vessel U.
    pub wall_storage_basis: ThermalStorageBasis,
    /// Explicit mechanical constraint used to interpret that coordinate.
    pub mechanical_constraint: WallStreamMechanicalConstraint,
    /// Initial and accepted original wall temperatures [K].
    pub wall_temperature_k: [f64; 2],
    /// Actual original H/U storage change, counted once [J].
    pub wall_storage_change_j: f64,
    /// Signed outward wall pressure-volume work [J].
    pub wall_pressure_volume_work_j: f64,
    /// Wall Delta U; approximate only in the explicitly selected negligible-
    /// expansion mode. Exactly Delta H - p Delta V in the isobaric mode.
    pub wall_internal_energy_change_j: f64,
    /// Represented Delta U plus pressure-volume work minus stored Delta E [J].
    pub mechanical_energy_residual_j: f64,
    /// Original density-derived initial/final volumes in isobaric mode [m3].
    pub wall_volume_m3: Option<[f64; 2]>,
    /// Integrated original stream enthalpy changes in port order [J].
    pub stream_heat_j: Vec<f64>,
    /// Actual wall storage change plus all stream heat [J].
    pub energy_residual_j: f64,
    /// Actual storage change minus the proposed scalar heat coordinate [J].
    pub wall_target_residual_j: f64,
    /// Largest absolute integrated stream contact-law residual [J].
    pub maximum_contact_residual_j: f64,
    /// Original material entropy integral [J/K].
    pub wall_entropy_change_j_k: f64,
    /// Integrated original stream entropy changes in port order [J/K].
    pub stream_entropy_change_j_k: Vec<f64>,
    /// Discrete wall-plus-stream entropy account, never clipped [J/K].
    /// Includes backward-Euler temporal dissipation; this does not isolate the
    /// exact continuous-time physical entropy production.
    pub entropy_production_j_k: f64,
    /// Original wall entropy minus Delta E_wall / T_wall_end [J/K].
    /// This nonnegative concavity defect quantifies the temporal entropy
    /// contribution of the backward-Euler wall update, subject to tolerances.
    pub wall_temporal_entropy_defect_j_k: f64,
    /// Original wall-state trial count, including initial and bracket trials.
    pub wall_evaluations: usize,
    /// Stream solves started across all wall trials.
    pub stream_evaluations: usize,
    /// Scalar duty queries across all stream solves.
    pub duty_evaluations: usize,
    /// Original PT calls across preparations and every nested trial.
    pub pt_evaluations: usize,
    /// Sum of original EOS work ceilings, not a measured EOS count.
    pub eos_evaluation_upper_bound: usize,
}

#[derive(Default)]
struct Work { wall: usize, streams: usize, duty: usize, pt: usize, eos: usize }
impl Work {
    fn exhausted(&self) -> Error {
        Error::WorkLimit { wall_evaluations: self.wall, duty_evaluations: self.duty, pt_evaluations: self.pt }
    }
    fn stream_error(&self, port: usize, source: MixtureHeatExchangerError) -> Error {
        if source == MixtureHeatExchangerError::Cancelled { Error::Cancelled } else {
            Error::Stream { port, completed_duty_evaluations: self.duty,
                completed_pt_evaluations: self.pt, source }
        }
    }
}
struct Trial {
    proposed_heat: f64,
    contact: WallContact,
    temperature: f64,
    change: f64,
    streams: Vec<MixtureWallStreamStep>,
    residual: f64,
    scale: f64,
}

/// Advance one finite original wall and one or two perfectly mixed steady ports.
///
/// This is an explicit model reduction: fluids have negligible holdup, each
/// chamber has a uniform outlet temperature, pressure and feed stay fixed, and
/// the wall is spatially uniform. For each wall trial the original PH ports
/// resolve their instantaneous outlet states. Backward Euler closes
/// `Delta E_wall + dt * sum(mass_flow * Delta h_stream) = 0` with actual source
/// increments. There is no additional latent term or stream pressure-work term.
///
/// The wall's own conductance must be zero: interface conductances belong only
/// to ports. Original material entropy is required. The final combined entropy
/// uses the wall integral and stream increments, not another boundary entropy
/// debit. This discrete entropy account includes temporal dissipation from
/// backward Euler; it is not an isolated continuous-time physical production.
/// The method is first-order in wall time; refine the interval to resolve
/// transients. Zero duration leaves the wall unchanged while returning the
/// instantaneous steady outlets. Nonzero throughflow can coexist with zero
/// net wall heat. Any refusal, including final cancellation, leaves the wall
/// unchanged. The source-domain brackets must support the admitted trials.
pub fn advance_wall_with_mixture_streams(
    wall: &mut SharedWallStorage<'_>,
    ports: &[MixtureWallPort<'_>],
    duration: Time,
    mechanical_constraint: WallStreamMechanicalConstraint,
    controls: FiniteWallStreamControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FiniteWallStreamStep, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let dt = duration.value();
    if !dt.is_finite() || dt < 0.0 || !(1..=2).contains(&ports.len()) {
        return Err(Error::Invalid("finite nonnegative duration and one or two flowing ports required"));
    }
    let conductance = match wall.thermal {
        WallThermal::Lumped { conductance_w_per_k, .. }
        | WallThermal::Caloric { conductance_w_per_k, .. } => conductance_w_per_k,
        _ => return Err(Error::Unsupported("finite original material wall required")),
    };
    if conductance != 0.0 {
        return Err(Error::Invalid("shared wall conductance belongs exclusively to ports"));
    }
    let source = wall.heat_storage.as_ref().ok_or(Error::Unsupported("original material storage required"))?;
    let initial = source.bind_contact(WallContact { drag: WallDrag::None, thermal: wall.thermal }).map_err(wall_error)?;
    let initial_temperature = source.temperature_k(initial.thermal).map_err(wall_error)?;
    let basis = source.storage_basis();
    let _ = mechanical_account(source, initial, initial, mechanical_constraint, 0.0)?;
    let _ = wall_entropy(source, initial, initial)?;
    for (index, port) in ports.iter().enumerate() {
        if !port.conductance.value().is_finite() || port.conductance.value() < 0.0 {
            return Err(Error::Invalid("nonnegative finite port conductance required"));
        }
        port.controls.validate().map_err(|source| Error::Stream {
            port: index, completed_duty_evaluations: 0, completed_pt_evaluations: 0, source,
        })?;
    }
    let mut work = Work::default();
    let mut prepared = Vec::with_capacity(ports.len());
    for (index, port) in ports.iter().enumerate() {
        poll(&mut cancelled)?;
        let remaining = controls.maximum_pt_evaluations - work.pt;
        if remaining == 0 { return Err(work.exhausted()); }
        let mut control = port.controls;
        if dt != 0.0 {
            // Each port's native root/contact gets a share of the energy budget.
            control.absolute_power_tolerance_w = control.absolute_power_tolerance_w
                .min(0.125 * controls.absolute_energy_tolerance_j / dt / ports.len() as f64);
            control.relative_power_tolerance = control.relative_power_tolerance
                .min(0.125 * controls.relative_energy_tolerance);
        }
        control.maximum_pt_evaluations = control.maximum_pt_evaluations.min(remaining);
        let stream = PreparedMixtureWallStream::new(port.inlet, control, &mut cancelled)
            .map_err(|e| work.stream_error(index, e))?;
        work.pt += stream.preparation_pt_evaluations();
        work.eos += stream.preparation_eos_evaluation_upper_bound();
        prepared.push(stream);
    }
    let mut evaluate = |heat: f64, work: &mut Work| -> Result<Trial, Error> {
        poll(&mut cancelled)?;
        if work.wall == controls.maximum_wall_evaluations { return Err(work.exhausted()); }
        work.wall += 1;
        let point = source.after(initial, heat).map_err(wall_error)?;
        let mut streams = Vec::with_capacity(ports.len());
        let mut integrated_heat = [0.0; 2];
        let mut scale = point.energy_change_j.abs();
        for (index, (port, prepared)) in ports.iter().zip(&prepared).enumerate() {
            if work.duty == controls.maximum_duty_evaluations { return Err(work.exhausted()); }
            work.streams += 1;
            let stream = prepared.evaluate(
                Temperature::new(point.temperature_k), port.conductance,
                controls.maximum_duty_evaluations - work.duty,
                controls.maximum_pt_evaluations - work.pt, &mut cancelled,
            ).map_err(|e| work.stream_error(index, e))?;
            work.duty += stream.duty_evaluations;
            work.pt += stream.pt_evaluations;
            work.eos += stream.eos_evaluation_upper_bound;
            let integrated = dt * stream.enthalpy_flow_change_w;
            if !integrated.is_finite()
                || dt != 0.0 && stream.enthalpy_flow_change_w != 0.0 && integrated == 0.0
            {
                return Err(Error::Invalid("unrepresentable integrated stream heat"));
            }
            integrated_heat[index] = integrated;
            scale = scale.max(integrated.abs());
            streams.push(stream);
        }
        let residual = sum_three(point.energy_change_j, integrated_heat[0], integrated_heat[1]);
        if !residual.is_finite() { return Err(Error::Invalid("unrepresentable wall energy residual")); }
        Ok(Trial { proposed_heat: heat, contact: point.contact, temperature: point.temperature_k,
            change: point.energy_change_j, streams, residual, scale })
    };
    let zero = evaluate(0.0, &mut work)?;
    let accepted = if dt == 0.0 || zero.residual.abs() <= 0.5 * controls.tolerance(zero.scale) {
        zero
    } else {
        let passive_bound = -zero.residual;
        let mut low = controls.wall_heat_bracket_j[0].max(passive_bound.min(0.0));
        let mut high = controls.wall_heat_bracket_j[1].min(passive_bound.max(0.0));
        let mut source_low: f64 = 0.0;
        let mut source_high: f64 = 0.0;
        for port in ports {
            let limits = source.heat_limits(initial, port.inlet.temperature.value()).map_err(wall_error)?;
            source_low = source_low.min(-limits[1]);
            source_high = source_high.max(-limits[0]);
        }
        low = low.max(source_low);
        high = high.min(source_high);
        if !low.is_finite() || !high.is_finite() || low >= high {
            return Err(Error::Invalid("empty supported finite-wall heat bracket"));
        }
        let endpoint = evaluate(if passive_bound > 0.0 { high } else { low }, &mut work)?;
        let (mut lower, mut upper) = if passive_bound > 0.0 {
            ([low, zero.residual], [high, endpoint.residual])
        } else {
            ([low, endpoint.residual], [high, zero.residual])
        };
        let mut candidate = endpoint;
        loop {
            let tolerance = 0.5 * controls.tolerance(candidate.scale);
            if !positive(tolerance) { return Err(Error::Invalid("unrepresentable wall energy tolerance")); }
            if candidate.residual.abs() <= tolerance { break candidate; }
            if lower[1] > 0.0 || upper[1] < 0.0 {
                return Err(Error::NotBracketed {
                    wall_heat_bracket_j: [lower[0], upper[0]], residual_j: [lower[1], upper[1]],
                });
            }
            let fraction = -lower[1] / (upper[1] - lower[1]);
            let fraction = if (0.05..0.95).contains(&fraction) { fraction } else { 0.5 };
            let heat = lower[0] + fraction * (upper[0] - lower[0]);
            if heat == lower[0] || heat == upper[0] { return Err(Error::NoConvergence); }
            candidate = evaluate(heat, &mut work)?;
            if candidate.residual < 0.0 { lower = [heat, candidate.residual]; }
            else { upper = [heat, candidate.residual]; }
        }
    };
    let target_residual = accepted.change - accepted.proposed_heat;
    let maximum_contact_residual = accepted.streams.iter()
        .map(|s| dt * s.heat_transfer_residual_w.abs()).fold(0.0, f64::max);
    let tolerance = controls.tolerance(accepted.scale);
    let (pressure_work, internal_change, volume) =
        mechanical_account(source, initial, accepted.contact, mechanical_constraint, accepted.change)?;
    let mechanical_residual = sum_three(internal_change, pressure_work, -accepted.change);
    if !positive(tolerance)
        || [accepted.residual, target_residual, maximum_contact_residual, mechanical_residual].iter()
            .any(|r| !r.is_finite() || r.abs() > tolerance)
    {
        return Err(Error::Balance { energy_residual_j: accepted.residual,
            wall_target_residual_j: target_residual,
            maximum_contact_residual_j: maximum_contact_residual,
            mechanical_energy_residual_j: mechanical_residual });
    }
    let entropy_wall = wall_entropy(source, initial, accepted.contact)?;
    let stream_heat: Vec<_> = accepted.streams.iter().map(|s| dt * s.enthalpy_flow_change_w).collect();
    let stream_entropy: Vec<_> = accepted.streams.iter().map(|s| dt * s.entropy_flow_change_w_k).collect();
    let production = sum_three(entropy_wall, stream_entropy[0], stream_entropy.get(1).copied().unwrap_or(0.0));
    let temporal_defect = entropy_wall - accepted.change / accepted.temperature;
    if !temporal_defect.is_finite() {
        return Err(Error::Invalid("unrepresentable original wall temporal entropy defect"));
    }
    if !production.is_finite() || production < -controls.entropy_tolerance_j_k {
        return Err(Error::Entropy { production_j_k: production, tolerance_j_k: controls.entropy_tolerance_j_k });
    }
    poll(&mut cancelled)?;
    // This is the only mutation. A zero-time query preserves even cache fields.
    if dt != 0.0 { wall.thermal = accepted.contact.thermal; }
    Ok(FiniteWallStreamStep {
        elapsed: duration, streams: accepted.streams,
        wall_storage_basis: basis, mechanical_constraint,
        wall_temperature_k: [initial_temperature, accepted.temperature],
        wall_storage_change_j: accepted.change,
        wall_pressure_volume_work_j: pressure_work,
        wall_internal_energy_change_j: internal_change,
        mechanical_energy_residual_j: mechanical_residual,
        wall_volume_m3: volume,
        stream_heat_j: stream_heat,
        energy_residual_j: accepted.residual,
        wall_target_residual_j: target_residual,
        maximum_contact_residual_j: maximum_contact_residual,
        wall_entropy_change_j_k: entropy_wall,
        stream_entropy_change_j_k: stream_entropy,
        entropy_production_j_k: production,
        wall_temporal_entropy_defect_j_k: temporal_defect,
        wall_evaluations: work.wall, stream_evaluations: work.streams,
        duty_evaluations: work.duty, pt_evaluations: work.pt,
        eos_evaluation_upper_bound: work.eos,
    })
}

fn wall_entropy(
    source: &MaterialWallCaloric<'_>, before: WallContact, after: WallContact,
) -> Result<f64, Error> {
    let entropy = if let Some(law) = source.thermal_law() {
        let (mass, initial) = energy_state(before)?;
        let (_, final_energy) = energy_state(after)?;
        mass * law.specific_entropy_change_j_kg_k(initial, final_energy)
            .map_err(|e| wall_error(e.into()))?
            .ok_or(Error::Unsupported("original wall entropy integral required"))?
    } else if let Some(capacity) = source.storage() {
        capacity.sensible_entropy_change_j_k(
            source.temperature_k(before.thermal).map_err(wall_error)?,
            source.temperature_k(after.thermal).map_err(wall_error)?,
        ).map_err(|e| wall_error(e.into()))?
    } else { return Err(Error::Unsupported("original wall entropy integral required")); };
    if !entropy.is_finite() { return Err(Error::Invalid("nonfinite original wall entropy")); }
    Ok(entropy)
}

fn mechanical_account(
    source: &MaterialWallCaloric<'_>, before: WallContact, after: WallContact,
    constraint: WallStreamMechanicalConstraint, change: f64,
) -> Result<(f64, f64, Option<[f64; 2]>), Error> {
    match (source.storage_basis(), constraint) {
        (ThermalStorageBasis::InternalEnergy, WallStreamMechanicalConstraint::IsochoricInternalEnergy)
        | (ThermalStorageBasis::Enthalpy, WallStreamMechanicalConstraint::NegligibleExpansionEnthalpy) =>
            Ok((0.0, change, None)),
        (ThermalStorageBasis::Enthalpy, WallStreamMechanicalConstraint::IsobaricEnthalpy { pressure }) => {
            if !positive(pressure.value()) { return Err(Error::Invalid("positive authored wall pressure required")); }
            let law = source.thermal_law().ok_or(Error::Unsupported("isobaric wall requires original density law"))?;
            let (mass, initial) = energy_state(before)?;
            let (_, final_energy) = energy_state(after)?;
            let density = |energy| -> Result<f64, Error> {
                let rho = law.bulk_density_kg_m3(energy).map_err(|e| wall_error(e.into()))?
                    .ok_or(Error::Unsupported("original wall density required for pressure-volume work"))?;
                if !positive(rho) { return Err(Error::Invalid("positive original wall density required")); }
                Ok(rho)
            };
            let rho0 = density(initial)?;
            let rho1 = density(final_energy)?;
            let volumes = [mass / rho0, mass / rho1];
            // Form the density difference before dividing to retain small changes.
            let delta_volume = (rho0 - rho1) / rho1 * volumes[0];
            let work = pressure.value() * delta_volume;
            let internal = change - work;
            if volumes.iter().any(|v| !positive(*v)) || !delta_volume.is_finite()
                || !work.is_finite() || !internal.is_finite()
                || rho0 != rho1 && delta_volume == 0.0
                || delta_volume != 0.0 && work == 0.0
            {
                return Err(Error::Invalid("unrepresentable original wall pressure-volume account"));
            }
            Ok((work, internal, Some(volumes)))
        }
        _ => Err(Error::Unsupported("wall storage and explicit mechanical constraint are incompatible")),
    }
}
fn energy_state(contact: WallContact) -> Result<(f64, f64), Error> {
    match contact.thermal {
        WallThermal::Caloric { mass_kg, specific_energy_j_kg, .. } => Ok((mass_kg, specific_energy_j_kg)),
        _ => Err(Error::Unsupported("original mass and specific-energy wall state required")),
    }
}
fn wall_error(error: WallExchangeError) -> Error {
    if error == WallExchangeError::Cancelled { Error::Cancelled } else { Error::Wall(error) }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

// Retain the rounding error from both additions in a three-term balance.
fn sum_three(a: f64, b: f64, c: f64) -> f64 {
    let (ab, error_ab) = two_sum(a, b);
    let (sum, error_c) = two_sum(ab, c);
    sum + (error_ab + error_c)
}
