//! Finite-rate chemistry during prescribed adiabatic compression or expansion.
//!
//! Each time interval applies the original material's isentropic half-volume
//! map, the original finite-rate rigid reactor, then another half-volume map.
//! This symmetric splitting resolves pressure work and reaction heat together
//! without replacing the source calorics with a frozen heat-capacity ratio.
//! Volume follows a constant logarithmic rate. The caller controls temporal
//! splitting error through `substeps`; the chemistry owner still controls its
//! own integration error. These are homogeneous gas states, without piston
//! inertia, shocks, wall heat transfer or phase change.

use super::{
    ClosedGasParcel, MaterialKineticsError, MaterialKineticsSolve, ThermochemicalGasModel,
    advance_closed_parcel_adiabatic, poll, positive,
};
use fs_material::thermochemical_gas::energy::entropy::{GasEntropyError, GasEntropySolve};
use fs_qty::{Time, Volume};
use fs_thermochem::kinetics::IdealGasKineticsV1;

/// Maximum prescribed trajectory intervals in one atomic update.
pub const MAX_COMPRESSION_SUBSTEPS: usize = 4096;

/// Explicit controls for chemistry, reversible volume work and time splitting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialCompressionSolve {
    /// Native finite-rate integration and material endpoint tolerances.
    pub kinetics: MaterialKineticsSolve,
    /// Original caloric/entropy inversion with an explicit temperature bracket.
    pub volume: GasEntropySolve,
    /// Equal time intervals along the logarithmic volume trajectory, 1..=4096.
    /// A fixed-volume request uses one native reactor call regardless of this.
    pub substeps: usize,
}

/// Actual work and balances from the accepted original-owner operations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialCompressionStep {
    /// Requested physical elapsed time.
    pub elapsed_time: Time,
    /// Volume before the stroke [m³].
    pub initial_volume_m3: f64,
    /// Prescribed volume at the end of the stroke [m³].
    pub final_volume_m3: f64,
    /// Pressure work delivered to the gas, positive for compression [J].
    pub boundary_work_on_gas_j: f64,
    /// Sum of original-owner stored energy changes, including chemistry's
    /// material reconciliation residual [J]. Volume increments use the source's
    /// reference-invariant caloric difference, not subtraction of large U values.
    pub internal_energy_change_j: f64,
    /// Stored energy change minus boundary work; there is no external heat [J].
    /// Accumulated directly from chemical material-reconciliation residuals,
    /// without subtracting potentially much larger work/energy totals.
    pub energy_residual_j: f64,
    /// Sum of original material volume-step and native reactor entropy changes
    /// [J/K]. Native chemical entries precede material endpoint reconciliation;
    /// its entropy correction is excluded. This is a stage diagnostic, not an
    /// independently bounded entropy difference of the published parcel.
    pub entropy_change_j_k: f64,
    /// Actual splitting intervals (one for a fixed-volume reactor).
    pub substeps: usize,
    /// Accepted native chemical integration steps across all intervals.
    pub chemical_accepted_steps: usize,
    /// Rejected native chemical integration steps across all intervals.
    pub chemical_rejected_steps: usize,
}

/// Refusal leaves the caller's parcel unchanged, including all source receipts.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialCompressionError {
    /// Invalid time, volume, work budget or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original chemistry or material reconstruction refused a state.
    Kinetics(MaterialKineticsError),
    /// Original material entropy closure refused a volume step.
    Volume(GasEntropyError),
    /// Caller cancelled before publication.
    Cancelled,
}
impl core::fmt::Display for MaterialCompressionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting material compression: {self:?}")
    }
}
impl std::error::Error for MaterialCompressionError {}
impl From<MaterialKineticsError> for MaterialCompressionError {
    fn from(error: MaterialKineticsError) -> Self {
        match error {
            MaterialKineticsError::Cancelled => Self::Cancelled,
            other => Self::Kinetics(other),
        }
    }
}
impl From<GasEntropyError> for MaterialCompressionError {
    fn from(error: GasEntropyError) -> Self {
        match error {
            GasEntropyError::Cancelled => Self::Cancelled,
            other => Self::Volume(other),
        }
    }
}

/// Advance a sealed, adiabatic, reacting gas through a prescribed volume stroke.
///
/// For positive duration, `V(t) = V_initial exp(t ln(V_final/V_initial)/duration)`.
/// Symmetric composition of exact subflows is second order in the splitting
/// interval; the chosen native chemistry integrator may impose a lower order.
/// There is no automatic estimate of splitting error: refine `substeps` when
/// compression and chemical time scales compete. Every intermediate state must
/// remain inside the original species/transport domains and supplied brackets.
/// Formation energies are already in the native calorics; no reaction heat is
/// added separately. All changes publish together after the last cancel poll.
pub fn advance_closed_parcel_compression(
    parcel: &mut ClosedGasParcel,
    bank: &ThermochemicalGasModel,
    mechanism: &IdealGasKineticsV1,
    final_volume: Volume,
    duration: Time,
    controls: MaterialCompressionSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MaterialCompressionStep, MaterialCompressionError> {
    poll(&mut cancelled)?;
    let initial_volume = parcel.volume_m3();
    let final_volume = final_volume.value();
    let duration_s = duration.value();
    if !positive(final_volume)
        || !duration_s.is_finite()
        || duration_s < 0.0
        || !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&controls.substeps)
    {
        return Err(MaterialCompressionError::Invalid(
            "finite positive volume, nonnegative duration and 1..=4096 substeps required",
        ));
    }
    let moving = final_volume != initial_volume;
    if moving && duration_s == 0.0 {
        return Err(MaterialCompressionError::Invalid(
            "a finite compression or expansion stroke requires positive elapsed time",
        ));
    }
    let ratio = final_volume / initial_volume;
    if !positive(ratio) {
        return Err(MaterialCompressionError::Invalid(
            "volume ratio must remain finite and positive",
        ));
    }
    let log_ratio = fs_math::det::ln(ratio);
    let steps = if moving { controls.substeps } else { 1 };
    let mut candidate = parcel.clone();
    // Admit original volume controls even for the fixed-volume/no-op path.
    candidate.set_volume_isentropic(initial_volume, controls.volume, &mut cancelled)?;
    let mut work = (0.0, 0.0);
    let mut energy = (0.0, 0.0);
    let mut residual = (0.0, 0.0);
    let mut entropy = (0.0, 0.0);
    let mut accepted = 0;
    let mut rejected = 0;
    let mut elapsed = 0.0;
    for index in 0..steps {
        poll(&mut cancelled)?;
        let end_time = if index + 1 == steps {
            duration_s
        } else {
            duration_s * ((index + 1) as f64 / steps as f64)
        };
        let dt = end_time - elapsed;
        if duration_s > 0.0 && !positive(dt) {
            return Err(MaterialCompressionError::Invalid(
                "trajectory interval does not advance representable physical time",
            ));
        }
        let midpoint =
            initial_volume * fs_math::det::exp(log_ratio * ((index as f64 + 0.5) / steps as f64));
        let endpoint = if index + 1 == steps {
            final_volume
        } else {
            initial_volume * fs_math::det::exp(log_ratio * ((index + 1) as f64 / steps as f64))
        };
        let first = candidate.set_volume_isentropic(midpoint, controls.volume, &mut cancelled)?;
        let chemical = advance_closed_parcel_adiabatic(
            &mut candidate,
            bank,
            mechanism,
            Time::new(dt),
            controls.kinetics,
            &mut cancelled,
        )?;
        let second = candidate.set_volume_isentropic(endpoint, controls.volume, &mut cancelled)?;
        accumulate(&mut work, -first.work_by_gas_j);
        accumulate(&mut work, -second.work_by_gas_j);
        accumulate(&mut energy, first.stored_energy_change_j);
        accumulate(&mut energy, chemical.balance.energy_residual_j);
        accumulate(&mut energy, second.stored_energy_change_j);
        accumulate(&mut residual, chemical.balance.energy_residual_j);
        accumulate(&mut entropy, first.entropy_change_j_k);
        accumulate(
            &mut entropy,
            chemical.chemistry.balance().entropy_change_j_k,
        );
        accumulate(&mut entropy, second.entropy_change_j_k);
        accepted += chemical.chemistry.balance().accepted_steps;
        rejected += chemical.chemistry.balance().rejected_steps;
        elapsed = end_time;
    }
    let report = MaterialCompressionStep {
        elapsed_time: duration,
        initial_volume_m3: initial_volume,
        final_volume_m3: final_volume,
        boundary_work_on_gas_j: work.0,
        internal_energy_change_j: energy.0,
        energy_residual_j: residual.0,
        entropy_change_j_k: entropy.0,
        substeps: steps,
        chemical_accepted_steps: accepted,
        chemical_rejected_steps: rejected,
    };
    if !report.boundary_work_on_gas_j.is_finite()
        || !report.internal_energy_change_j.is_finite()
        || !report.energy_residual_j.is_finite()
        || !report.entropy_change_j_k.is_finite()
    {
        return Err(MaterialCompressionError::Invalid(
            "accumulated work, energy or entropy is unrepresentable",
        ));
    }
    poll(&mut cancelled)?;
    *parcel = candidate;
    Ok(report)
}

fn accumulate((sum, correction): &mut (f64, f64), value: f64) {
    let adjusted = value - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}
