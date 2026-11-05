//! Nonideal hydraulic hold-up with unequal tanks, bypass and stagnant exchange.
//!
//! For each tank, active and stagnant constituent masses obey
//! `dMa/dt = Qin*Cin - Qtrain*Ca - qex*(Ca-Cs)` and
//! `dMs/dt = qex*(Ca-Cs)`. Downstream tanks receive the preceding active outlet.
//! Total feed splits into `Qbypass = alpha*Qtotal` and
//! `Qtrain = (1-alpha)*Qtotal`; the bypass rejoins the final outlet instantly.
//! Active residence therefore uses Qtrain, not Qtotal. Zero exchange means an
//! isolated dead volume; positive exchange retains and later releases material.
//!
//! Model families: Fogler, *Elements of Chemical Reaction Engineering*, chapter
//! 18, sections 18.7-18.8, <https://public.websites.umich.edu/~elements/5e/18chap/Fogler_Web_Ch18_final.pdf>;
//! Dagadu et al., *Applied Radiation and Isotopes* 70 (2012), 156-161,
//! <https://doi.org/10.1016/j.apradiso.2011.09.003>. These citations establish
//! model families, not numerical parameters or a calibration for any plant.
//!
//! The caller declares geometry, inventories, provenance and all interval flows.
//! This is passive transport of one conserved constituent, not leach recovery,
//! tracer deconvolution, a two-phase slip model or a fit to operating data.
//! Exchange may continue while external flow stops. A throughput-volume clock
//! is equivalent only when every internal flow scales with the external flow.
//!
//! Positive exponential propagation integrates outlet and both exchange fluxes
//! independently. Conservation is checked without a balancing plug. A refusal
//! or cancellation never publishes a partial state. Binary64 results are not
//! interval enclosures.

use super::{MAX_INTERVAL_EXPOSURE, MAX_TANKS};
use fs_matdb::Provenance;
use fs_math::det;
use fs_qty::{Density, Dimensionless, Mass, Time, Volume, VolumetricFlowRate};

const STEP_EXPOSURE: f64 = 8.0;
const POISSON_TERMS: usize = 256;
const BALANCE_TOLERANCE: f64 = 2.0e-11;

/// Invalid declarations, unsupported numeric range, work bound or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum NonidealError {
    /// Geometry, inventory, flow, source or array dimensions are invalid.
    Invalid(&'static str),
    /// A physical quantity overflows or a positive path is unrepresentable.
    Arithmetic(&'static str),
    /// Split the interval; the maximum transition-rate exposure exceeds 4096.
    WorkLimit,
    /// Independently computed physical accounts did not close.
    Conservation,
    /// The caller cancelled before a complete state was admitted.
    Cancelled,
}
impl core::fmt::Display for NonidealError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonideal tank residence: {self:?}")
    }
}
impl std::error::Error for NonidealError {}
type Result<T> = core::result::Result<T, NonidealError>;

/// Fixed filled volumes of one tank; the stagnant region has no direct outlet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TankGeometry {
    active: Volume,
    stagnant: Volume,
}
impl TankGeometry {
    /// Positive active volume and nonnegative stagnant volume, in cubic metres.
    pub fn new(active: Volume, stagnant: Volume) -> Result<Self> {
        nonnegative(active.value(), "active volume")?;
        nonnegative(stagnant.value(), "stagnant volume")?;
        if active.value() == 0.0 {
            return Err(NonidealError::Invalid("active volume must be positive"));
        }
        finite(active.value() + stagnant.value(), "total tank volume")?;
        Ok(Self { active, stagnant })
    }
    /// Filled mobile volume.
    pub const fn active_volume(self) -> Volume {
        self.active
    }
    /// Filled stagnant volume, including an isolated dead region at zero exchange.
    pub const fn stagnant_volume(self) -> Volume {
        self.stagnant
    }
}

/// Mass of the same named passive constituent in the two physical regions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZoneInventory {
    /// Mass in the region connected to the main flow, kg.
    pub active: Mass,
    /// Mass in the stagnant region, kg.
    pub stagnant: Mass,
}

/// Piecewise-constant externally declared flow conditions, retained in receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidenceFlows {
    total: VolumetricFlowRate,
    bypass: Dimensionless,
    exchange: Vec<VolumetricFlowRate>,
}
impl ResidenceFlows {
    /// Exchange is the equal volumetric counterflow in each direction, not its
    /// signed net constituent flux. Supply one value per tank, including zeros.
    pub fn new(
        total: VolumetricFlowRate,
        bypass: Dimensionless,
        exchange: Vec<VolumetricFlowRate>,
    ) -> Result<Self> {
        nonnegative(total.value(), "total feed flow")?;
        if !bypass.value().is_finite() || !(0.0..=1.0).contains(&bypass.value()) {
            return Err(NonidealError::Invalid("bypass fraction must be in [0,1]"));
        }
        if exchange.is_empty() || exchange.len() > MAX_TANKS {
            return Err(NonidealError::Invalid("requires 1..=64 exchange flows"));
        }
        for flow in &exchange {
            nonnegative(flow.value(), "stagnant exchange flow")?;
        }
        Ok(Self {
            total,
            bypass,
            exchange,
        })
    }
    /// Total feed, before bypass splitting, in cubic metres per second.
    pub const fn total_flow(&self) -> VolumetricFlowRate {
        self.total
    }
    /// Instantaneously bypassed fraction of total feed.
    pub const fn bypass_fraction(&self) -> Dimensionless {
        self.bypass
    }
    /// Original per-tank exchange flows, in physical order.
    pub fn exchange_flows(&self) -> &[VolumetricFlowRate] {
        &self.exchange
    }
}

/// Complete interval accounts. Outlet and exchange masses are integrated fluxes.
#[derive(Clone, Debug, PartialEq)]
pub struct NonidealInterval {
    pub duration: Time,
    pub flows: ResidenceFlows,
    pub inlet_concentration: Density,
    pub total_inlet_volume: Volume,
    pub train_volume: Volume,
    pub bypass_volume: Volume,
    pub inlet_mass: Mass,
    pub train_inlet_mass: Mass,
    pub train_outlet_mass: Mass,
    pub bypass_mass: Mass,
    pub outlet_mass: Mass,
    /// Flow-weighted combined outlet concentration; absent at zero output volume.
    pub outlet_mean: Option<Density>,
    pub initial_inventory: ZoneInventory,
    pub final_inventory: ZoneInventory,
    pub transferred_to_stagnant: Mass,
    pub transferred_from_stagnant: Mass,
    pub active_residual: Mass,
    pub stagnant_residual: Mass,
    pub conservation_residual: Mass,
    pub substeps: usize,
}

/// A retained passive constituent in up to 64 unequal full tanks.
#[derive(Clone, Debug, PartialEq)]
pub struct NonidealTankCascade {
    constituent: String,
    geometry: Vec<TankGeometry>,
    inventory: Vec<ZoneInventory>,
    provenance: Provenance,
}
impl NonidealTankCascade {
    /// Declare the original geometry and both initial inventories explicitly.
    pub fn new(
        constituent: impl Into<String>,
        geometry: Vec<TankGeometry>,
        inventory: Vec<ZoneInventory>,
        provenance: Provenance,
    ) -> Result<Self> {
        let constituent = constituent.into();
        if constituent.trim().is_empty()
            || constituent.len() > 256
            || geometry.is_empty()
            || geometry.len() > MAX_TANKS
            || geometry.len() != inventory.len()
            || provenance.source.trim().is_empty()
            || provenance.license.trim().is_empty()
        {
            return Err(NonidealError::Invalid(
                "requires named constituent, matching 1..=64 tanks and original provenance",
            ));
        }
        for (g, m) in geometry.iter().zip(&inventory) {
            nonnegative(m.active.value(), "initial active mass")?;
            nonnegative(m.stagnant.value(), "initial stagnant mass")?;
            if g.stagnant.value() == 0.0 && m.stagnant.value() != 0.0 {
                return Err(NonidealError::Invalid(
                    "mass cannot occupy zero stagnant volume",
                ));
            }
            ratio(
                m.active.value(),
                g.active.value(),
                "initial active concentration",
            )?;
            if g.stagnant.value() > 0.0 {
                ratio(
                    m.stagnant.value(),
                    g.stagnant.value(),
                    "initial stagnant concentration",
                )?;
            }
        }
        totals(&inventory)?;
        Ok(Self {
            constituent,
            geometry,
            inventory,
            provenance,
        })
    }
    pub fn constituent(&self) -> &str {
        &self.constituent
    }
    pub fn geometry(&self) -> &[TankGeometry] {
        &self.geometry
    }
    pub fn inventories(&self) -> &[ZoneInventory] {
        &self.inventory
    }
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }
    pub fn total_inventory(&self) -> Result<ZoneInventory> {
        totals(&self.inventory)
    }

    /// Advance with declared constant flows and feed over this interval.
    pub fn advance(
        &mut self,
        duration: Time,
        flows: &ResidenceFlows,
        feed: Density,
    ) -> Result<NonidealInterval> {
        self.advance_with_cancel(duration, flows, feed, || false)
    }

    /// Cancellation is checked before work, per substep and every 16 series
    /// terms. All state changes are committed together after physical closure.
    pub fn advance_with_cancel(
        &mut self,
        duration: Time,
        flows: &ResidenceFlows,
        feed: Density,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NonidealInterval> {
        nonnegative(duration.value(), "interval duration")?;
        nonnegative(feed.value(), "feed concentration")?;
        if flows.exchange.len() != self.geometry.len() {
            return Err(NonidealError::Invalid(
                "one exchange flow is required per tank",
            ));
        }
        if cancelled() {
            return Err(NonidealError::Cancelled);
        }
        let qtrain = product(
            flows.total.value(),
            1.0 - flows.bypass.value(),
            "train flow",
        )?;
        let qbypass = product(flows.total.value(), flows.bypass.value(), "bypass flow")?;
        let total_volume = product(flows.total.value(), duration.value(), "total inlet volume")?;
        let train_volume = product(qtrain, duration.value(), "train volume")?;
        let bypass_volume = product(qbypass, duration.value(), "bypass volume")?;
        let inlet = product(total_volume, feed.value(), "inlet mass")?;
        let train_input = product(train_volume, feed.value(), "train inlet mass")?;
        let bypass = product(bypass_volume, feed.value(), "bypass mass")?;
        let initial = totals(&self.inventory)?;
        let mut state = self.inventory.clone();
        let (lambda, transitions) = transitions(&self.geometry, flows, qtrain)?;
        let mut outputs = [0.0; 3];
        let mut steps = 0;
        if duration.value() > 0.0 && lambda > 0.0 {
            let exposure = product(lambda, duration.value(), "transition exposure")?;
            if exposure > MAX_INTERVAL_EXPOSURE {
                return Err(NonidealError::WorkLimit);
            }
            steps = (exposure / STEP_EXPOSURE).ceil() as usize;
            let h = exposure / steps as f64;
            let input = ratio(train_input, steps as f64, "substep inlet mass")?;
            for _ in 0..steps {
                if cancelled() {
                    return Err(NonidealError::Cancelled);
                }
                let (next, fluxes) = propagate(&state, input, h, &transitions, &mut cancelled)?;
                for (out, increment) in outputs.iter_mut().zip(fluxes) {
                    *out = finite(*out + increment, "integrated outlet/exchange")?;
                }
                state = next;
            }
        }
        let final_inventory = totals(&state)?;
        // Concentrations, as well as masses, must remain representable.
        for (g, m) in self.geometry.iter().zip(&state) {
            ratio(
                m.active.value(),
                g.active.value(),
                "final active concentration",
            )?;
            if g.stagnant.value() > 0.0 {
                ratio(
                    m.stagnant.value(),
                    g.stagnant.value(),
                    "final stagnant concentration",
                )?;
            }
        }
        let outlet = finite(outputs[0] + bypass, "combined outlet mass")?;
        let mean = if total_volume > 0.0 {
            Some(Density::new(ratio(
                outlet,
                total_volume,
                "combined outlet concentration",
            )?))
        } else {
            None
        };
        let active_residual = balance(
            (final_inventory.active.value() - initial.active.value())
                + (outputs[0] - train_input)
                + (outputs[1] - outputs[2]),
            &[
                initial.active.value(),
                final_inventory.active.value(),
                train_input,
                outputs[0],
                outputs[1],
                outputs[2],
            ],
        )?;
        let stagnant_residual = balance(
            (final_inventory.stagnant.value() - initial.stagnant.value())
                + (outputs[2] - outputs[1]),
            &[
                initial.stagnant.value(),
                final_inventory.stagnant.value(),
                outputs[1],
                outputs[2],
            ],
        )?;
        let start = initial.active.value() + initial.stagnant.value();
        let end = final_inventory.active.value() + final_inventory.stagnant.value();
        let conservation_residual = balance(
            (end - start) + (outlet - inlet),
            &[start, end, inlet, outlet],
        )?;
        if cancelled() {
            return Err(NonidealError::Cancelled);
        }
        let result = NonidealInterval {
            duration,
            flows: flows.clone(),
            inlet_concentration: feed,
            total_inlet_volume: Volume::new(total_volume),
            train_volume: Volume::new(train_volume),
            bypass_volume: Volume::new(bypass_volume),
            inlet_mass: Mass::new(inlet),
            train_inlet_mass: Mass::new(train_input),
            train_outlet_mass: Mass::new(outputs[0]),
            bypass_mass: Mass::new(bypass),
            outlet_mass: Mass::new(outlet),
            outlet_mean: mean,
            initial_inventory: initial,
            final_inventory,
            transferred_to_stagnant: Mass::new(outputs[1]),
            transferred_from_stagnant: Mass::new(outputs[2]),
            active_residual,
            stagnant_residual,
            conservation_residual,
            substeps: steps,
        };
        self.inventory = state;
        Ok(result)
    }
}

#[derive(Clone, Copy)]
struct Transition {
    through: f64,
    into: f64,
    from: f64,
    active_hold: f64,
    stagnant_hold: f64,
}

fn transitions(
    geometry: &[TankGeometry],
    flows: &ResidenceFlows,
    qtrain: f64,
) -> Result<(f64, Vec<Transition>)> {
    let mut raw = Vec::with_capacity(geometry.len());
    let mut lambda: f64 = 0.0;
    for (g, exchange) in geometry.iter().zip(&flows.exchange) {
        if g.stagnant.value() == 0.0 && exchange.value() != 0.0 {
            return Err(NonidealError::Invalid(
                "positive exchange requires positive stagnant volume",
            ));
        }
        let through = ratio(qtrain, g.active.value(), "active discharge rate")?;
        let into = ratio(exchange.value(), g.active.value(), "active exchange rate")?;
        let from = if g.stagnant.value() > 0.0 {
            ratio(
                exchange.value(),
                g.stagnant.value(),
                "stagnant exchange rate",
            )?
        } else {
            0.0
        };
        let active_exit = finite(through + into, "total active transition rate")?;
        lambda = lambda.max(active_exit).max(from);
        raw.push((through, into, from, active_exit));
    }
    if lambda == 0.0 {
        return Ok((lambda, Vec::new()));
    }
    let transitions = raw
        .into_iter()
        .map(|(through, into, from, active_exit)| {
            Ok(Transition {
                through: ratio(through, lambda, "normalized discharge rate")?,
                into: ratio(into, lambda, "normalized active exchange rate")?,
                from: ratio(from, lambda, "normalized stagnant exchange rate")?,
                active_hold: (lambda - active_exit) / lambda,
                stagnant_hold: (lambda - from) / lambda,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((lambda, transitions))
}

// exp(A*dt) = exp(-h) sum_j h^j P^j/j!, P=I+A/lambda.
// Three absorbing rewards integrate the train outlet and directional exchange.
// The affine input enters as input*p_(j-1)/j. Scale by physical mass before
// recurrence so a tiny bare probability does not erase a representable flux.
fn propagate(
    state: &[ZoneInventory],
    input: f64,
    h: f64,
    rates: &[Transition],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<ZoneInventory>, [f64; 3])> {
    let n = state.len();
    let mut weight = det::exp(-h);
    let mut term = vec![0.0; 2 * n + 3];
    for (i, (m, r)) in state.iter().zip(rates).enumerate() {
        // Disconnected compartments have the exact identity solution. Exclude
        // them from the series so an unrelated fast tank cannot perturb or
        // underflow their permanently retained inventory.
        if r.through > 0.0 || r.into > 0.0 || r.from > 0.0 {
            term[2 * i] = product(weight, m.active.value(), "scaled active mass")?;
        }
        if r.into > 0.0 || r.from > 0.0 {
            term[2 * i + 1] = product(weight, m.stagnant.value(), "scaled stagnant mass")?;
        }
    }
    let mut sum = term.clone();
    let mut next = vec![0.0; term.len()];
    for j in 1..POISSON_TERMS {
        if j % 16 == 0 && cancelled() {
            return Err(NonidealError::Cancelled);
        }
        next.fill(0.0);
        for (i, r) in rates.iter().enumerate() {
            let active = term[2 * i];
            let stagnant = term[2 * i + 1];
            next[2 * i] += active * r.active_hold + stagnant * r.from;
            next[2 * i + 1] += stagnant * r.stagnant_hold + active * r.into;
            let dest = if i + 1 == n { 2 * n } else { 2 * (i + 1) };
            next[dest] += active * r.through;
            next[2 * n + 1] += active * r.into;
            next[2 * n + 2] += stagnant * r.from;
        }
        for i in 2 * n..2 * n + 3 {
            next[i] += term[i];
        }
        let factor = h / j as f64;
        for value in &mut next {
            *value *= factor;
        }
        next[0] += input * (weight / j as f64);
        for (total, value) in sum.iter_mut().zip(&next) {
            *total = finite(*total + value, "positive propagation")?;
        }
        std::mem::swap(&mut term, &mut next);
        weight *= factor;
    }
    let mut upstream = input > 0.0;
    let mut into_path = false;
    let mut from_path = false;
    for (i, (m, r)) in state.iter().zip(rates).enumerate() {
        if r.through == 0.0 && r.into == 0.0 && r.from == 0.0 {
            sum[2 * i] = m.active.value();
        }
        if r.into == 0.0 && r.from == 0.0 {
            sum[2 * i + 1] = m.stagnant.value();
        }
        let active =
            upstream || m.active.value() > 0.0 || (m.stagnant.value() > 0.0 && r.from > 0.0);
        let stagnant = m.stagnant.value() > 0.0 || (active && r.into > 0.0);
        if (active && sum[2 * i] == 0.0) || (stagnant && sum[2 * i + 1] == 0.0) {
            return Err(NonidealError::Arithmetic(
                "positive retained inventory underflow",
            ));
        }
        into_path |= active && r.into > 0.0;
        from_path |= stagnant && r.from > 0.0;
        upstream = active && r.through > 0.0;
    }
    if (upstream && sum[2 * n] == 0.0)
        || (into_path && sum[2 * n + 1] == 0.0)
        || (from_path && sum[2 * n + 2] == 0.0)
    {
        return Err(NonidealError::Arithmetic(
            "positive outlet/exchange underflow",
        ));
    }
    let inventories = (0..n)
        .map(|i| ZoneInventory {
            active: Mass::new(sum[2 * i]),
            stagnant: Mass::new(sum[2 * i + 1]),
        })
        .collect();
    Ok((inventories, [sum[2 * n], sum[2 * n + 1], sum[2 * n + 2]]))
}

fn totals(state: &[ZoneInventory]) -> Result<ZoneInventory> {
    let (mut active, mut stagnant) = (0.0, 0.0);
    for m in state {
        active = finite(active + m.active.value(), "total active mass")?;
        stagnant = finite(stagnant + m.stagnant.value(), "total stagnant mass")?;
    }
    finite(active + stagnant, "total retained mass")?;
    Ok(ZoneInventory {
        active: Mass::new(active),
        stagnant: Mass::new(stagnant),
    })
}
fn nonnegative(value: f64, context: &'static str) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(NonidealError::Invalid(context))
    }
}
fn finite(value: f64, context: &'static str) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(NonidealError::Arithmetic(context))
    }
}
fn product(a: f64, b: f64, context: &'static str) -> Result<f64> {
    let value = finite(a * b, context)?;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        Err(NonidealError::Arithmetic(context))
    } else {
        Ok(value)
    }
}
fn ratio(a: f64, b: f64, context: &'static str) -> Result<f64> {
    let value = finite(a / b, context)?;
    if a > 0.0 && value == 0.0 {
        Err(NonidealError::Arithmetic(context))
    } else {
        Ok(value)
    }
}
fn balance(residual: f64, terms: &[f64]) -> Result<Mass> {
    let scale = terms.iter().copied().fold(0.0, f64::max);
    if !residual.is_finite() || residual.abs() > BALANCE_TOLERANCE * scale {
        Err(NonidealError::Conservation)
    } else {
        Ok(Mass::new(residual))
    }
}

#[cfg(test)]
#[path = "nonideal/tests.rs"]
mod tests;
