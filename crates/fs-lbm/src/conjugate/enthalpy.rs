//! Stationary Ohmic heating through the existing conservative enthalpy timeline.
//!
//! Current is recomputed at every thermal substep's accepted start, using both
//! temperature and specific enthalpy. Thus a phase fraction can change electrical
//! response at constant temperature without an apparent infinite heat capacity.
//! The original electrical contact deposition and thermal integration own all
//! numerical work. This adapter neither solves melt flow nor changes cell mass.

/// Solved material-driven face fluxes on the same conservative heat timeline.
pub mod advection;
pub mod transport;

use std::collections::{BTreeMap, BTreeSet};
use fs_conduction::transient::deposition::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation,
    EnergyAudit, Event, StepConfig,
};
use fs_exec::CancelGate;
use super::{ChtError, ElectricSetup, ElectricSolution, ElectricalContact, VoxelDomain};

/// Fixed spatial/current inputs during one advance. Drive changes must coincide
/// with separate advance endpoints; zero voltage is a short, not an open circuit.
#[derive(Clone, Copy)]
pub struct JouleEnthalpyInput<'a> {
    /// Existing electrical geometry. Exterior faces remain insulating.
    pub domain: &'a VoxelDomain,
    /// Original voltage/current electrodes; solid presets are ignored.
    pub setup: &'a ElectricSetup,
    /// Original area-specific resistance and explicit heat partition.
    pub contacts: &'a [ElectricalContact],
    /// One unique active thermal cell per voxel, in x-fastest electrical order.
    /// Each volume and centreline midpoint must match the receiving voxel.
    /// Unmapped thermal cells may exchange heat but receive no electrical source.
    pub cell_ids: &'a [CellId],
    /// Electrical cell allocation ceiling, separate from thermal work limits.
    pub max_cells: usize,
    /// Current-solve ceiling for the WHOLE advance, including final diagnostics.
    /// There is one solve per actual thermal substep plus one endpoint solve.
    pub max_current_solves: usize,
    /// Absolute tolerance for each electrode/Joule/mapped-power check [W].
    pub power_tolerance_w: f64,
    /// Absolute accumulated electrical and transfer discrepancy limit [J].
    /// Independent of the thermal owner's `StepConfig::energy_tolerance_j`.
    pub energy_tolerance_j: f64,
}

/// Actual electrical evaluation coordinates and the corresponding current solve.
#[derive(Clone, Debug, PartialEq)]
pub struct EnthalpyCurrentState {
    /// Physical evaluation time [s].
    pub time_s: f64,
    /// Actual receiving-cell temperatures [K], including electrode cells.
    pub temperature_k: Vec<f64>,
    /// Actual receiving-cell specific enthalpies [J/kg], not inferred from T.
    pub specific_enthalpy_j_kg: Vec<f64>,
    /// Exact callback resistivities [ohm m]. None means explicit insulation OR
    /// an ideal electrode, for which no bulk electrical property is requested.
    pub resistivity_ohm_m: Vec<Option<f64>>,
    /// Original source-producing solution, with bulk/contact heat counted once.
    pub electric: ElectricSolution,
}

/// Accounts from one accepted advance; these are numerical checks, not source,
/// phase-kinetic, space/time accuracy or experimental certificates.
#[derive(Clone, Debug, PartialEq)]
pub struct JouleEnthalpyReport {
    /// Original thermal audit. Its external heat ALREADY includes generated
    /// Joule heat; never add `deposited_joule_j` to it a second time.
    pub thermal: EnergyAudit,
    /// Sum of applied electrode work over the actual thermal substeps [J].
    pub electrode_energy_j: f64,
    /// Integrated electrical solver's total irreversible heat [J].
    pub joule_energy_j: f64,
    /// Sum of mapped cell powers integrated over the SAME substeps [J].
    pub deposited_joule_j: f64,
    /// Electrode work minus electrical Joule heat [J].
    pub electrical_defect_j: f64,
    /// Electrical Joule heat minus mapped thermal source energy [J].
    pub transfer_defect_j: f64,
    /// Actual current solves, including the endpoint support/diagnostic solve.
    pub current_solves: usize,
    /// Last source actually applied; None for a zero-duration advance.
    pub last_applied: Option<EnthalpyCurrentState>,
    /// Duration over which the last source was applied, in seconds.
    pub last_duration_s: Option<f64>,
    /// Fresh electrical/material evaluation at the accepted endpoint. This did
    /// NOT supply heat during the preceding explicit step and is not substituted
    /// for the actual applied fields in either energy account.
    pub endpoint: EnthalpyCurrentState,
}

/// Original electrical/thermal failures and coupling admission refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum JouleEnthalpyError {
    /// Bad map, limits, power handoff or accumulated electrical discrepancy.
    Invalid(&'static str),
    /// Whole-advance current work budget exhausted before the next evaluation.
    CurrentBudget,
    /// Original electrical or caller material-law failure.
    Electric(ChtError),
    /// Original enthalpy-domain, thermal, observer or cancellation failure.
    Thermal(DepositionError),
}
impl core::fmt::Display for JouleEnthalpyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "Joule/enthalpy coupling: {s}"),
            Self::CurrentBudget => write!(f, "Joule/enthalpy current-solve budget exhausted"),
            Self::Electric(e) => write!(f, "Joule/enthalpy coupling: {e}"),
            Self::Thermal(e) => write!(f, "Joule/enthalpy coupling: {e}"),
        }
    }
}
impl std::error::Error for JouleEnthalpyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Electric(e) => Some(e), Self::Thermal(e) => Some(e), _ => None }
    }
}
impl From<ChtError> for JouleEnthalpyError {
    fn from(e: ChtError) -> Self { Self::Electric(e) }
}
impl From<DepositionError> for JouleEnthalpyError {
    fn from(e: DepositionError) -> Self { Self::Thermal(e) }
}
fn invalid(s: &'static str) -> JouleEnthalpyError { JouleEnthalpyError::Invalid(s) }
fn finite(x: f64) -> Result<f64, JouleEnthalpyError> {
    if x.is_finite() { Ok(x) } else { Err(invalid("nonfinite source or accumulated energy")) }
}
fn electrode_power(e: &super::ElectrodeSolution) -> Result<f64, JouleEnthalpyError> {
    // An isolated ideal zero-current terminal has no defined voltage, but does
    // exactly zero work. Do not invent a potential or accept any other NaN/Inf.
    if e.voltage_v.is_nan() && e.current_a == 0.0 { return Ok(0.0); }
    finite(e.voltage_v)?;
    finite(e.current_a)?;
    finite(e.voltage_v * e.current_a)
}
fn checkpoint(gate: &CancelGate) -> Result<(), DepositionError> {
    if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) }
}

fn admit_map(thermal: &DepositionSimulation<'_>, input: JouleEnthalpyInput<'_>, gate: &CancelGate)
    -> Result<(), JouleEnthalpyError> {
    checkpoint(gate)?;
    let n = input.domain.cell_count();
    if n > input.max_cells || input.cell_ids.len() != n || input.max_current_solves == 0
        || !input.power_tolerance_w.is_finite() || input.power_tolerance_w < 0.0
        || !input.energy_tolerance_j.is_finite() || input.energy_tolerance_j < 0.0 {
        return Err(invalid("invalid electrical cell map or solve/power/energy limits"));
    }
    let mut seen = BTreeSet::new();
    for (c, &id) in input.cell_ids.iter().enumerate() {
        if c % 256 == 0 { checkpoint(gate)?; }
        let cell = thermal.cells().get(&id).ok_or(invalid("electrical map needs an active thermal cell"))?;
        if !seen.insert(id) { return Err(invalid("electrical cells cannot share a thermal receiving cell")); }
        let volume = input.domain.volume(c);
        let actual = cell.geometry.volume_m3;
        if !volume.is_finite() || volume <= 0.0
            || (actual-volume).abs() > 64.0*f64::EPSILON*actual.max(volume) {
            return Err(invalid("electrical and thermal cell volumes do not match"));
        }
        let [x,y,z] = input.domain.coords(c);
        let centre = input.domain.center(x,y,z);
        let widths = input.domain.widths(c);
        for a in 0..3 {
            let mid = 0.5*cell.geometry.start_m[a] + 0.5*cell.geometry.end_m[a];
            if !centre[a].is_finite() || !widths[a].is_finite() || widths[a] <= 0.0
                || (mid-centre[a]).abs() > 64.0*f64::EPSILON*centre[a].abs().max(widths[a]) {
                return Err(invalid("electrical and thermal cell centres do not match"));
            }
        }
    }
    Ok(())
}

type ResistivityLaw<'a> = dyn Fn(usize, &DepositedCell, f64) -> Result<Option<f64>, ChtError> + 'a;

#[derive(Clone)]
struct Account<'a> {
    input: JouleEnthalpyInput<'a>,
    gate: &'a CancelGate,
    solves: usize,
    electrode_j: f64,
    joule_j: f64,
    mapped_j: f64,
    pending_mapped_w: f64,
    pending_electrode_w: f64,
    last: Option<EnthalpyCurrentState>,
    last_dt: Option<f64>,
}
impl Account<'_> {
    fn evaluate(&mut self, thermal: &DepositionSimulation<'_>, law: &ResistivityLaw<'_>) -> Result<EnthalpyCurrentState, JouleEnthalpyError> {
        self.evaluate_with_contacts(thermal, law, self.input.contacts)
    }

    // Keep a single current solve, power transfer and work-accounting owner.
    // Material-interface consumers pass contacts resolved at this SAME epoch;
    // the original entry point passes its unchanged authored contact slice.
    fn evaluate_with_contacts(&mut self, thermal: &DepositionSimulation<'_>,
        law: &ResistivityLaw<'_>, contacts: &[ElectricalContact],
    ) -> Result<EnthalpyCurrentState, JouleEnthalpyError> {
        self.evaluate_with_contact_topology(thermal, law, contacts, &[])
    }

    // Open and closed interfaces use the same current and energy accounting.
    // A missing material law still refuses; openings are an explicit state.
    fn evaluate_with_contact_topology(&mut self, thermal: &DepositionSimulation<'_>,
        law: &ResistivityLaw<'_>, contacts: &[ElectricalContact], open_faces: &[[usize; 2]],
    ) -> Result<EnthalpyCurrentState, JouleEnthalpyError> {
        checkpoint(self.gate)?;
        if self.solves >= self.input.max_current_solves { return Err(JouleEnthalpyError::CurrentBudget); }
        self.solves += 1;
        let mut temperature_k = Vec::with_capacity(self.input.cell_ids.len());
        let mut specific_enthalpy_j_kg = Vec::with_capacity(self.input.cell_ids.len());
        for (c,&id) in self.input.cell_ids.iter().enumerate() {
            if c % 256 == 0 { checkpoint(self.gate)?; }
            temperature_k.push(thermal.temperature_k(id)?);
            specific_enthalpy_j_kg.push(thermal.cells()[&id].specific_h_j_kg);
        }
        let mut resistivity_ohm_m = vec![None; self.input.cell_ids.len()];
        let electric = self.input.setup.solve_with_cell_contact_topology(
            self.input.domain, &temperature_k, contacts, open_faces,
            |c,_,t| {
                let value = (law)(c, &thermal.cells()[&self.input.cell_ids[c]], t)?;
                resistivity_ohm_m[c] = value;
                Ok(value)
            }, self.gate,
        )?;
        let summed = electric.joule_w.iter().try_fold(0.0, |sum,&q| {
            if !q.is_finite() || q < 0.0 { return Err(invalid("invalid electrical heat")); }
            finite(sum+q)
        })?;
        let electrode = electric.electrodes.iter().try_fold(0.0, |sum,e| finite(sum+electrode_power(e)?))?;
        if electric.joule_w.len() != self.input.cell_ids.len()
            || !electric.total_joule_w.is_finite() || !electric.delivered_w.is_finite()
            || (summed-electric.total_joule_w).abs() > self.input.power_tolerance_w
            || (electrode-electric.delivered_w).abs() > self.input.power_tolerance_w
            || (electrode-summed).abs() > self.input.power_tolerance_w {
            return Err(invalid("electrode work, electrical heat and cell deposition do not close"));
        }
        checkpoint(self.gate)?;
        Ok(EnthalpyCurrentState { time_s: thermal.time_s(), temperature_k, specific_enthalpy_j_kg,
            resistivity_ohm_m, electric })
    }
}
impl<'a> Account<'a> {
    fn new(input: JouleEnthalpyInput<'a>, gate: &'a CancelGate) -> Self {
        Self { input, gate, solves:0, electrode_j:0.0, joule_j:0.0, mapped_j:0.0,
            pending_mapped_w:0.0, pending_electrode_w:0.0, last:None, last_dt:None }
    }
    fn apply(&mut self, state: EnthalpyCurrentState) -> Result<BTreeMap<CellId,f64>,JouleEnthalpyError> {
        let mapped: BTreeMap<_,_> = self.input.cell_ids.iter().copied()
            .zip(state.electric.joule_w.iter().copied()).collect();
        self.pending_mapped_w = mapped.values().try_fold(0.0, |sum,&q| finite(sum+q))?;
        self.pending_electrode_w = state.electric.electrodes.iter()
            .try_fold(0.0, |sum,e| finite(sum+electrode_power(e)?))?;
        self.last = Some(state);
        self.last_dt = None;
        Ok(mapped)
    }
    fn record(&mut self, start_s:f64, dt:f64) -> Result<(),JouleEnthalpyError> {
        let state = self.last.as_ref().ok_or(invalid("thermal step has no electrical source"))?;
        if state.time_s != start_s || self.last_dt.is_some() || !dt.is_finite() || dt <= 0.0 {
            return Err(invalid("electrical source epoch does not match the thermal substep"));
        }
        self.electrode_j = finite(self.electrode_j+finite(dt*self.pending_electrode_w)?)?;
        self.joule_j = finite(self.joule_j+finite(dt*state.electric.total_joule_w)?)?;
        self.mapped_j = finite(self.mapped_j+finite(dt*self.pending_mapped_w)?)?;
        self.last_dt = Some(dt);
        Ok(())
    }
}

impl Account<'_> {
    fn finish(self, audit: EnergyAudit, endpoint: EnthalpyCurrentState) -> Result<JouleEnthalpyReport,JouleEnthalpyError> {
        let electrical_defect_j = finite(self.electrode_j-self.joule_j)?;
        let transfer_defect_j = finite(self.joule_j-self.mapped_j)?;
        if electrical_defect_j.abs() > self.input.energy_tolerance_j || transfer_defect_j.abs() > self.input.energy_tolerance_j {
            return Err(invalid("accumulated electrical/thermal source energy does not close"));
        }
        checkpoint(self.gate)?;
        Ok(JouleEnthalpyReport { thermal:audit, electrode_energy_j:self.electrode_j,
            joule_energy_j:self.joule_j, deposited_joule_j:self.mapped_j,
            electrical_defect_j, transfer_defect_j, current_solves:self.solves,
            last_applied:self.last, last_duration_s:self.last_dt, endpoint })
    }
}

#[derive(Clone)]
struct Source<'a> {
    account: Account<'a>,
    law: &'a ResistivityLaw<'a>,
}
impl DepositionObserver for Source<'_> {
    type Error = JouleEnthalpyError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(),Self::Error> { Ok(()) }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>, check: &mut dyn FnMut()->Result<(),DepositionError>)
        -> Result<BTreeMap<CellId,f64>,Self::Error> {
        check()?;
        let state = self.account.evaluate(thermal,self.law)?;
        self.account.apply(state)
    }
    fn step(&mut self, start_s:f64, dt:f64, _: &BTreeMap<CellId,f64>, _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut()->Result<(),DepositionError>) -> Result<(),Self::Error> {
        check()?;
        self.account.record(start_s,dt)
    }
}

/// Advance the existing enthalpy network with generated Ohmic heat, once per
/// actual substep. `resolve(voxel_index, thermal_cell, kelvin)` sees invariant
/// mass, material ID and specific enthalpy as well as T. A caller can resolve
/// phase fraction from its EXISTING enthalpy chart, including at fixed T, and
/// evaluate its admitted electrical law. `None` explicitly declares insulation;
/// missing data must be an error. The callback is not called for ideal electrodes.
/// Source receipts and any interior callback side effects remain caller-owned.
///
/// Thermal and electrical models share this explicit cell correspondence, not
/// a material-name rule. Electrical contacts remain distinct from thermal links.
/// Thermal links/geometry/laws retain their existing owners; they are not inferred
/// from electrical resistance. Initialize mapped cells before calling this API.
/// No accepted state is changed until the caller adopts the returned candidate.
/// Final electrical support, all original heat checks, the whole current-work
/// budget, electrical transfer and cancellation must pass before returning it.
///
/// This is explicit stationary Ohmic/caloric coupling, not a fluid solver,
/// implicit source iteration or moving phase front. The conduction CFL does
/// NOT bound source stiffness; refine `max_dt_s`. The caller must justify each
/// thermal control volume and its fixed-volume caloric approximation. No density
/// jump, pressure/expansion work, latent-heat double counting, geometry change,
/// magnetic induction, kinetic phase law or experimental validation is supplied.
#[allow(clippy::too_many_arguments)]
pub fn advance_joule_enthalpy<'t>(
    thermal: &DepositionSimulation<'t>, target_s: f64, config: StepConfig,
    input: JouleEnthalpyInput<'_>, gate: &CancelGate,
    resolve: impl Fn(usize, &DepositedCell, f64) -> Result<Option<f64>, ChtError>,
) -> Result<(DepositionSimulation<'t>, JouleEnthalpyReport), JouleEnthalpyError> {
    admit_map(thermal,input,gate)?;
    let source = Source { account:Account::new(input,gate), law:&resolve };
    let (next,mut source,audit) = thermal.advance_observed(&source,target_s,config,||checkpoint(gate))?;
    // Endpoint support contributes no additional heat.
    let endpoint = source.account.evaluate(&next,source.law)?;
    let report = source.account.finish(audit,endpoint)?;
    Ok((next,report))
}

