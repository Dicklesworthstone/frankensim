//! Transient first-order transfer of one conserved constituent in equal full CSTRs.
//!
//! Both undissolved and dissolved constituent remain in each tank and advect:
//! `dU_i/dt = a*(U_(i-1)-U_i) - k_i*U_i`,
//! `dD_i/dt = a*(D_(i-1)-D_i) + k_i*U_i`, with `a = Q/(V/N)`.
//! U and D denote masses of the SAME constituent (e.g. copper), not masses of
//! different compounds. Inlet concentrations supply the first tank. Conversion
//! is an internal transfer; dissolved effluent is not final refining recovery.
//!
//! Basis: Fogler, *Elements of Chemical Reaction Engineering*, chapter 1 general
//! mole balance and section 17.4 tanks in series / first-order conversion:
//! <https://websites.umich.edu/~elements/5e/01chap/summary.html>,
//! <https://websites.umich.edu/~elements/5e/17chap/Ch17_Web-Additional%20Material.pdf>.
//! This is an explicitly declared dilute, irreversible local law. No reagent
//! limitation, passivation, adsorption, solid/liquid slip, heat or changing
//! carrier volume is inferred. Both forms have the same hydraulic residence.
//!
//! Positive uniformization propagates the original state and independent
//! outlet/transfer integrals together. Numerical residuals are checked, never
//! plugged. The bounded floating-point calculation is not an interval proof.

use super::residence::{MAX_INTERVAL_EXPOSURE, MAX_TANKS};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::det;
use fs_qty::{Density, Frequency, Mass, Time, Volume, VolumetricFlowRate};

/// A local transient rate declaration, distinct from batch or stationary cards.
pub const TRANSIENT_TRANSFER_LAW: &str = "transient-first-order-constituent-transfer";
/// Nonnegative inverse-second transfer coefficient, with no ultimate-yield factor.
pub const TRANSFER_RATE_PARAMETER: &str = "transfer-rate";
const STEP_EXPOSURE: f64 = 8.0;
const POISSON_TERMS: usize = 256;
const BALANCE_TOLERANCE: f64 = 2.0e-11;

/// Unsupported declaration, operating point or numerical range; no state published.
#[derive(Clone, Debug, PartialEq)]
pub enum ReactiveError {
    /// Original material-card admission failure.
    Source(MatDbError),
    /// Invalid dimensions, geometry, identity, state or source declaration.
    Invalid(&'static str),
    /// At least one original local law does not admit the supplied conditions.
    OutsideValidity,
    /// A physical result or coefficient cannot be represented.
    Arithmetic(&'static str),
    /// Split the interval; combined transport/reaction exposure exceeds the bound.
    WorkLimit,
    /// Independently propagated phase balances did not close.
    Conservation,
}
impl core::fmt::Display for ReactiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient constituent transfer: {self:?}")
    }
}
impl std::error::Error for ReactiveError {}
impl From<MatDbError> for ReactiveError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
type Result<T> = core::result::Result<T, ReactiveError>;

/// Immutable source declaration for a stationary coefficient in a transient ODE.
/// The caller holds the admitted source conditions fixed over each step.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalTransferCard {
    source: ConstitutiveModelCard,
    rate: Frequency,
}
impl LocalTransferCard {
    /// Admit only this explicit law, one inverse-second coefficient, original
    /// source identities, and finite, typed applicability axes. Finite batch-age
    /// windows and stationary residence domains cannot authorize this local law.
    pub fn from_card(source: ConstitutiveModelCard) -> Result<Self> {
        source.validate()?;
        if source.law.0 != TRANSIENT_TRANSFER_LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.parameters.len() != 1
            || source.sources.is_empty()
            || source.provenance.artifact.is_none()
            || source.validity.bounds().iter().any(|(axis, (lo, hi))| {
                !lo.is_finite()
                    || !hi.is_finite()
                    || lo > hi
                    || !source.validity.axis_quantities().contains_key(axis)
                    || matches!(
                        axis.as_str(),
                        "elapsed-time" | "residence-time" | "mean-residence-time"
                    )
            })
        {
            return Err(ReactiveError::Invalid(
                "requires an explicit transient local law and finite typed source conditions",
            ));
        }
        let parameter = source
            .parameters
            .get(TRANSFER_RATE_PARAMETER)
            .ok_or(ReactiveError::Invalid("missing transfer-rate"))?;
        if parameter.dims != Frequency::DIMS
            || !parameter.value.is_finite()
            || parameter.value < 0.0
        {
            return Err(ReactiveError::Invalid(
                "transfer-rate must be finite nonnegative inverse seconds",
            ));
        }
        let rate = Frequency::new(parameter.value);
        Ok(Self { source, rate })
    }
    /// Complete original law; state propagation adds no calibration authority.
    #[must_use]
    pub const fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Source-declared coefficient, not an inferred ore or mineral constant.
    #[must_use]
    pub const fn rate(&self) -> Frequency {
        self.rate
    }
}

/// Mass of the same named constituent in two phases, in kilograms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConstituentInventory {
    /// Constituent still in the undissolved feed population.
    pub undissolved: Mass,
    /// Constituent transferred into the dissolved population, still present.
    pub dissolved: Mass,
}

/// Independent interval accounts; conversion transfers between phases internally.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveInterval {
    /// Elapsed physical time.
    pub duration: Time,
    /// Actual carrier flow, including zero during a stop.
    pub train_flow: VolumetricFlowRate,
    /// Per-tank fixed conditions admitted against the original source cards.
    pub conditions: Vec<QueryPoint>,
    /// Total phase inventories before this step.
    pub initial: ConstituentInventory,
    /// Total phase inventories after this step.
    pub final_inventory: ConstituentInventory,
    /// Constituents entering during this interval, independently calculated.
    pub inlet: ConstituentInventory,
    /// Constituents leaving the final tank, independently integrated.
    pub outlet: ConstituentInventory,
    /// Undissolved-to-dissolved transfer inside all tanks, independently integrated.
    pub converted: Mass,
    /// Undissolved: final - initial + outlet - inlet + converted.
    pub undissolved_residual: Mass,
    /// Dissolved: final - initial + outlet - inlet - converted.
    pub dissolved_residual: Mass,
    /// Bounded positive propagation substeps (zero for no physical evolution).
    pub substeps: usize,
}

/// Equal full ideal tanks with explicit local laws and retained phase inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveTankCascade {
    constituent: String,
    active_volume: Volume,
    laws: Vec<LocalTransferCard>,
    tanks: Vec<ConstituentInventory>,
}
impl ReactiveTankCascade {
    /// Initialize 1..=64 equal tanks. Laws and inventories are in upstream-to-
    /// downstream order. No stationary/batch card is automatically promoted.
    pub fn new(
        constituent: impl Into<String>,
        active_volume: Volume,
        laws: Vec<LocalTransferCard>,
        tanks: Vec<ConstituentInventory>,
    ) -> Result<Self> {
        let constituent = constituent.into();
        if constituent.trim().is_empty()
            || constituent.len() > 256
            || laws.is_empty()
            || laws.len() > MAX_TANKS
            || laws.len() != tanks.len()
            || !active_volume.value().is_finite()
            || active_volume.value() <= 0.0
        {
            return Err(ReactiveError::Invalid(
                "requires named constituent, positive volume, and matching 1..=64 laws/inventories",
            ));
        }
        for tank in &tanks {
            nonnegative(tank.undissolved.value(), "undissolved initial inventory")?;
            nonnegative(tank.dissolved.value(), "dissolved initial inventory")?;
        }
        if active_volume.value() / tanks.len() as f64 == 0.0 {
            return Err(ReactiveError::Arithmetic("individual tank volume"));
        }
        let result = Self {
            constituent,
            active_volume,
            laws,
            tanks,
        };
        result.inventory()?;
        Ok(result)
    }
    /// Caller-declared common constituent identity, not inferred mineral chemistry.
    #[must_use]
    pub fn constituent(&self) -> &str {
        &self.constituent
    }
    /// Active mixed volume of the complete train.
    #[must_use]
    pub const fn active_volume(&self) -> Volume {
        self.active_volume
    }
    /// Immutable upstream-to-downstream local law sources.
    #[must_use]
    pub fn laws(&self) -> &[LocalTransferCard] {
        &self.laws
    }
    /// Current retained masses, in upstream-to-downstream order.
    #[must_use]
    pub fn tanks(&self) -> &[ConstituentInventory] {
        &self.tanks
    }
    /// Total retained constituent, separated by phase.
    pub fn inventory(&self) -> Result<ConstituentInventory> {
        inventory(&self.tanks)
    }

    /// Advance fixed flow/inlet/conditions over one physical interval. A stopped
    /// flow admits continued transfer; neither dissolved product nor carrier
    /// disappears. All validation and balances complete before state publication.
    #[allow(clippy::too_many_lines)] // Keep admission, independent accounts and atomic publication together.
    pub fn advance(
        &mut self,
        duration: Time,
        actual_train_flow: VolumetricFlowRate,
        inlet_undissolved: Density,
        inlet_dissolved: Density,
        conditions: &[QueryPoint],
    ) -> Result<ReactiveInterval> {
        let (dt, flow) = (duration.value(), actual_train_flow.value());
        for (value, label) in [
            (dt, "duration"),
            (flow, "actual train flow"),
            (inlet_undissolved.value(), "undissolved inlet concentration"),
            (inlet_dissolved.value(), "dissolved inlet concentration"),
        ] {
            nonnegative(value, label)?;
        }
        if conditions.len() != self.laws.len() {
            return Err(ReactiveError::Invalid(
                "one operating point per local law is required",
            ));
        }
        for (law, point) in self.laws.iter().zip(conditions) {
            if !law
                .source
                .validity
                .contains_typed(point.axes(), point.axis_quantities())
            {
                return Err(ReactiveError::OutsideValidity);
            }
        }
        let initial = self.inventory()?;
        let zero = ConstituentInventory {
            undissolved: Mass::new(0.0),
            dissolved: Mass::new(0.0),
        };
        let mut result = ReactiveInterval {
            duration,
            train_flow: actual_train_flow,
            conditions: conditions.to_vec(),
            initial,
            final_inventory: initial,
            inlet: zero,
            outlet: zero,
            converted: Mass::new(0.0),
            undissolved_residual: Mass::new(0.0),
            dissolved_residual: Mass::new(0.0),
            substeps: 0,
        };
        if dt == 0.0 {
            return Ok(result);
        }
        let n = self.tanks.len();
        let cell_volume = self.active_volume.value() / n as f64;
        let advection = flow / cell_volume;
        if !advection.is_finite() || (flow > 0.0 && advection == 0.0) {
            return Err(ReactiveError::Arithmetic("hydraulic transfer rate"));
        }
        let lambda = self
            .laws
            .iter()
            .map(|law| advection + law.rate.value())
            .fold(0.0, f64::max);
        if !lambda.is_finite() {
            return Err(ReactiveError::Arithmetic("combined transfer rate"));
        }
        if lambda == 0.0 {
            return Ok(result);
        }
        let exposure = product(lambda, dt, "combined interval exposure")?;
        if exposure > MAX_INTERVAL_EXPOSURE {
            return Err(ReactiveError::WorkLimit);
        }
        let steps = (exposure / STEP_EXPOSURE).ceil().max(1.0) as usize;
        let h = exposure / steps as f64;
        let exchanged = product(flow, dt, "exchanged volume")?;
        result.inlet = ConstituentInventory {
            undissolved: Mass::new(product(
                exchanged,
                inlet_undissolved.value(),
                "undissolved inlet mass",
            )?),
            dissolved: Mass::new(product(
                exchanged,
                inlet_dissolved.value(),
                "dissolved inlet mass",
            )?),
        };
        let input = [
            result.inlet.undissolved.value() / steps as f64,
            result.inlet.dissolved.value() / steps as f64,
        ];
        if (result.inlet.undissolved.value() > 0.0 && input[0] == 0.0)
            || (result.inlet.dissolved.value() > 0.0 && input[1] == 0.0)
        {
            return Err(ReactiveError::Arithmetic("positive substep inlet mass"));
        }
        let advect = advection / lambda;
        let rates: Vec<_> = self
            .laws
            .iter()
            .map(|law| law.rate.value() / lambda)
            .collect();
        // Subtract the original rounded total hazard, before normalization;
        // separately rounded a/lambda + k/lambda may otherwise exceed one.
        let holds: Vec<_> = self
            .laws
            .iter()
            .map(|law| (lambda - (advection + law.rate.value())) / lambda)
            .collect();
        if (advection > 0.0 && advect == 0.0)
            || rates
                .iter()
                .zip(&self.laws)
                .any(|(r, l)| *r == 0.0 && l.rate.value() > 0.0)
        {
            return Err(ReactiveError::Arithmetic("positive transition coefficient"));
        }
        let mut state = self.tanks.clone();
        let mut outputs = [0.0; 3];
        for _ in 0..steps {
            let propagated = propagate(&state, input, h, advect, &rates, &holds)?;
            for (total, increment) in outputs.iter_mut().zip(&propagated[2 * n..]) {
                *total = finite(*total + increment, "integrated outlet/transfer")?;
            }
            state = (0..n)
                .map(|i| ConstituentInventory {
                    undissolved: Mass::new(propagated[2 * i]),
                    dissolved: Mass::new(propagated[2 * i + 1]),
                })
                .collect();
        }
        result.final_inventory = inventory(&state)?;
        result.outlet = ConstituentInventory {
            undissolved: Mass::new(outputs[0]),
            dissolved: Mass::new(outputs[1]),
        };
        result.converted = Mass::new(outputs[2]);
        let phase_balance =
            |start: f64, end: f64, inlet: f64, outlet: f64, transfer: f64| -> Result<Mass> {
                let residual = (end - start) + (outlet - inlet) + transfer;
                let scale = start.max(end).max(inlet).max(outlet).max(transfer.abs());
                if !residual.is_finite() || residual.abs() > BALANCE_TOLERANCE * scale {
                    return Err(ReactiveError::Conservation);
                }
                Ok(Mass::new(residual))
            };
        result.undissolved_residual = phase_balance(
            initial.undissolved.value(),
            result.final_inventory.undissolved.value(),
            result.inlet.undissolved.value(),
            outputs[0],
            outputs[2],
        )?;
        result.dissolved_residual = phase_balance(
            initial.dissolved.value(),
            result.final_inventory.dissolved.value(),
            result.inlet.dissolved.value(),
            outputs[1],
            -outputs[2],
        )?;
        result.substeps = steps;
        self.tanks = state;
        Ok(result)
    }
}

// Uniformization: exp(A*t) = exp(-h) sum h^j P^j/j!, P=I+A/lambda.
// Two absorbing outlet states and one transfer reward integrate their original
// fluxes. Source-scaled recurrence avoids forming tiny probabilities before
// multiplying them by physical masses. The affine source contributes input*p_j/(j+1).
fn propagate(
    state: &[ConstituentInventory],
    input: [f64; 2],
    h: f64,
    a: f64,
    rates: &[f64],
    holds: &[f64],
) -> Result<Vec<f64>> {
    let n = state.len();
    let mut weight = det::exp(-h);
    let mut term = vec![0.0; 2 * n + 3];
    for (i, tank) in state.iter().enumerate() {
        term[2 * i] = product(
            weight,
            tank.undissolved.value(),
            "scaled undissolved inventory",
        )?;
        term[2 * i + 1] = product(weight, tank.dissolved.value(), "scaled dissolved inventory")?;
    }
    let mut sum = term.clone();
    let mut next = vec![0.0; term.len()];
    for j in 1..POISSON_TERMS {
        let factor = h / j as f64;
        next.fill(0.0);
        for i in 0..n {
            let u = term[2 * i];
            let d = term[2 * i + 1];
            next[2 * i] += u * holds[i];
            next[2 * i + 1] += d * (1.0 - a) + u * rates[i];
            let destination = if i + 1 == n { 2 * n } else { 2 * (i + 1) };
            next[destination] += u * a;
            next[destination + 1] += d * a;
            next[2 * n + 2] += u * rates[i];
        }
        for i in 2 * n..2 * n + 3 {
            next[i] += term[i];
        }
        for value in &mut next {
            *value *= factor;
        }
        next[0] += input[0] * (weight / j as f64);
        next[1] += input[1] * (weight / j as f64);
        for (total, value) in sum.iter_mut().zip(&next) {
            *total = finite(*total + value, "positive propagation")?;
        }
        std::mem::swap(&mut term, &mut next);
        weight *= factor;
    }
    let mut u_path = a > 0.0 && input[0] > 0.0;
    let mut d_path = a > 0.0 && input[1] > 0.0;
    let mut conversion = false;
    for (i, tank) in state.iter().enumerate() {
        u_path = tank.undissolved.value() > 0.0 || (a > 0.0 && u_path);
        let reacts = u_path && rates[i] > 0.0;
        conversion |= reacts;
        d_path = tank.dissolved.value() > 0.0 || (a > 0.0 && d_path) || reacts;
        if (u_path && sum[2 * i] == 0.0) || (d_path && sum[2 * i + 1] == 0.0) {
            return Err(ReactiveError::Arithmetic(
                "positive retained phase underflow",
            ));
        }
    }
    if (a > 0.0 && u_path && sum[2 * n] == 0.0)
        || (a > 0.0 && d_path && sum[2 * n + 1] == 0.0)
        || (conversion && sum[2 * n + 2] == 0.0)
    {
        return Err(ReactiveError::Arithmetic(
            "positive outlet/transfer underflow",
        ));
    }
    Ok(sum)
}
fn inventory(tanks: &[ConstituentInventory]) -> Result<ConstituentInventory> {
    let mut totals = [0.0; 2];
    for tank in tanks {
        totals[0] = finite(
            totals[0] + tank.undissolved.value(),
            "total undissolved inventory",
        )?;
        totals[1] = finite(
            totals[1] + tank.dissolved.value(),
            "total dissolved inventory",
        )?;
    }
    finite(totals[0] + totals[1], "total constituent inventory")?;
    Ok(ConstituentInventory {
        undissolved: Mass::new(totals[0]),
        dissolved: Mass::new(totals[1]),
    })
}
fn finite(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(ReactiveError::Arithmetic(what))
    }
}
fn nonnegative(value: f64, what: &'static str) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(ReactiveError::Invalid(what))
    }
}
fn product(a: f64, b: f64, what: &'static str) -> Result<f64> {
    let value = finite(a * b, what)?;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        Err(ReactiveError::Arithmetic(what))
    } else {
        Ok(value)
    }
}

#[cfg(test)]
#[path = "reactive_tests.rs"]
mod tests;
