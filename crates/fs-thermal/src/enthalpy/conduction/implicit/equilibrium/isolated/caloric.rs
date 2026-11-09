//! Common-temperature equilibrium of a closed original caloric inventory.
//!
//! The conserved coordinate is the caller's original stored specific energy,
//! not an inferred internal energy, enthalpy, heat capacity or phase fraction.
//! Each immutable T(e) must be continuous and nondecreasing on its bounds.
//! Set-valued inverses preserve flat latent intervals. Their allocation is
//! accepted only when conservation determines it, or the seed already balances.

use fs_math::dd::Dd;

/// One fixed parcel and its original continuous caloric domain, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricInventory {
    pub mass_kg: f64,
    pub initial_specific_energy_j_kg: f64,
    pub specific_energy_bounds_j_kg: [f64; 2],
}

/// Temperature and the original one-sided caloric derivative at the same e.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricState {
    pub temperature_k: f64,
    pub temperature_slope_kg_k_j: f64,
}

/// Physical acceptance and shared work limits, including refused queries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricEquilibriumControls {
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub maximum_temperature_iterations: usize,
    pub maximum_inverse_iterations: usize,
    pub maximum_caloric_evaluations: usize,
}

impl Default for CaloricEquilibriumControls {
    fn default() -> Self {
        Self {
            temperature_tolerance_k: 1e-9,
            energy_tolerance_j: 1e-9,
            maximum_temperature_iterations: 96,
            maximum_inverse_iterations: 128,
            maximum_caloric_evaluations: 1_000_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum CaloricEquilibriumError<E> {
    Invalid(&'static str),
    /// Includes original material-domain errors and checkpoint cancellation.
    Source(E),
    EvaluationBudget {
        used: usize,
        maximum: usize,
    },
    InverseNotResolved {
        cell: usize,
    },
    NoCommonTemperature,
    NoConservedEnergyRoot,
    /// Several independent latent allocations remain after the one energy row.
    UnderdeterminedPhaseAllocation {
        cells: usize,
    },
    NotConverged,
    Unrepresentable(&'static str),
}
impl<E: core::fmt::Debug> core::fmt::Display for CaloricEquilibriumError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "closed caloric equilibrium: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for CaloricEquilibriumError<E> {}

#[derive(Clone, Debug, PartialEq)]
pub struct CaloricEquilibrium {
    pub specific_energies_j_kg: Vec<f64>,
    pub temperatures_k: Vec<f64>,
    /// Sum m_i * (e_i - e_i,initial), evaluated in the original coordinates.
    pub energy_residual_j: f64,
    pub temperature_spread_k: f64,
    pub temperature_iterations: usize,
    pub caloric_evaluations: usize,
    pub retained_initial_state: bool,
}

#[derive(Clone, Copy)]
struct Sample {
    energy: f64,
    state: CaloricState,
}

struct Queries<F, C> {
    source: F,
    checkpoint: C,
    used: usize,
    maximum: usize,
}
impl<E, F, C> Queries<F, C>
where
    F: FnMut(usize, f64) -> Result<CaloricState, E>,
    C: FnMut() -> Result<(), E>,
{
    fn poll(&mut self) -> Result<(), CaloricEquilibriumError<E>> {
        (self.checkpoint)().map_err(CaloricEquilibriumError::Source)
    }
    fn get(&mut self, cell: usize, energy: f64) -> Result<Sample, CaloricEquilibriumError<E>> {
        self.poll()?;
        if self.used >= self.maximum {
            return Err(CaloricEquilibriumError::EvaluationBudget {
                used: self.used,
                maximum: self.maximum,
            });
        }
        self.used += 1;
        let state = (self.source)(cell, energy).map_err(CaloricEquilibriumError::Source)?;
        if !state.temperature_k.is_finite()
            || state.temperature_k <= 0.0
            || !state.temperature_slope_kg_k_j.is_finite()
            || state.temperature_slope_kg_k_j < 0.0
        {
            return Err(CaloricEquilibriumError::Invalid(
                "invalid original caloric state",
            ));
        }
        Ok(Sample { energy, state })
    }
}

#[derive(Clone, Copy)]
struct Inverse {
    lower: f64,
    upper: f64,
    slope: f64,
    latent: bool,
    hint: Option<f64>,
}

// Find BOTH generalized inverses. A positive derivative at a plateau endpoint
// may be a legitimate one-sided choice; it cannot erase the adjacent latent
// interval. Binary searches stop at adjacent representable source coordinates.
fn invert<E, F, C>(
    cell: usize,
    temperature: f64,
    bounds: [Sample; 2],
    controls: CaloricEquilibriumControls,
    queries: &mut Queries<F, C>,
) -> Result<Inverse, CaloricEquilibriumError<E>>
where
    F: FnMut(usize, f64) -> Result<CaloricState, E>,
    C: FnMut() -> Result<(), E>,
{
    let mut hint = None;
    let mut search = |lower: bool| -> Result<Sample, CaloricEquilibriumError<E>> {
        let [mut lo, mut hi] = bounds;
        if lower && lo.state.temperature_k == temperature {
            return Ok(lo);
        }
        if !lower && hi.state.temperature_k == temperature {
            return Ok(hi);
        }
        for _ in 0..controls.maximum_inverse_iterations {
            let energy = lo.energy.midpoint(hi.energy);
            if energy == lo.energy || energy == hi.energy {
                return Ok(if lower { hi } else { lo });
            }
            let sample = queries.get(cell, energy)?;
            let t = sample.state.temperature_k;
            if t < lo.state.temperature_k || t > hi.state.temperature_k {
                return Err(CaloricEquilibriumError::Invalid(
                    "original T(e) is not monotone",
                ));
            }
            if sample.state.temperature_slope_kg_k_j == 0.0 && t != temperature {
                hint = Some(t);
            }
            if t < temperature || (!lower && t == temperature) {
                lo = sample;
            } else {
                hi = sample;
            }
        }
        Err(CaloricEquilibriumError::InverseNotResolved { cell })
    };
    let lower = search(true)?;
    let upper = search(false)?;
    if lower.energy <= upper.energy {
        let middle = queries.get(cell, lower.energy.midpoint(upper.energy))?;
        if middle.state.temperature_k != temperature {
            return Err(CaloricEquilibriumError::Invalid(
                "inconsistent caloric inverse interval",
            ));
        }
        let latent = lower.energy < upper.energy && middle.state.temperature_slope_kg_k_j == 0.0;
        // A merely rounded sensible inverse is not a latent phase interval.
        return Ok(Inverse {
            lower: if latent { lower.energy } else { middle.energy },
            upper: if latent { upper.energy } else { middle.energy },
            slope: middle.state.temperature_slope_kg_k_j,
            latent,
            hint,
        });
    }
    let closest = if (lower.state.temperature_k - temperature).abs()
        <= (upper.state.temperature_k - temperature).abs()
    {
        lower
    } else {
        upper
    };
    if (closest.state.temperature_k - temperature).abs() > 0.5 * controls.temperature_tolerance_k {
        return Err(CaloricEquilibriumError::InverseNotResolved { cell });
    }
    Ok(Inverse {
        lower: closest.energy,
        upper: closest.energy,
        slope: closest.state.temperature_slope_kg_k_j,
        latent: false,
        hint,
    })
}

fn finite_dd<E>(value: Dd) -> Result<Dd, CaloricEquilibriumError<E>> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(CaloricEquilibriumError::Unrepresentable(
            "original-coordinate energy arithmetic",
        ))
    }
}

fn energy_change<E>(
    inventory: &[CaloricInventory],
    energies: &[f64],
) -> Result<Dd, CaloricEquilibriumError<E>> {
    let mut sum = Dd::ZERO;
    for (cell, &energy) in inventory.iter().zip(energies) {
        let difference =
            finite_dd(Dd::from_f64(energy) - Dd::from_f64(cell.initial_specific_energy_j_kg))?;
        let term = finite_dd(Dd::from_f64(cell.mass_kg) * difference)?;
        if difference != Dd::ZERO && term == Dd::ZERO {
            return Err(CaloricEquilibriumError::Unrepresentable(
                "nonzero energy change underflows",
            ));
        }
        sum = finite_dd(sum + term)?;
    }
    Ok(sum)
}

fn finish<E, F, C>(
    inventory: &[CaloricInventory],
    energies: Vec<f64>,
    controls: CaloricEquilibriumControls,
    iterations: usize,
    retained: bool,
    queries: &mut Queries<F, C>,
) -> Result<CaloricEquilibrium, CaloricEquilibriumError<E>>
where
    F: FnMut(usize, f64) -> Result<CaloricState, E>,
    C: FnMut() -> Result<(), E>,
{
    let residual = energy_change::<E>(inventory, &energies)?.to_f64();
    let mut temperatures = Vec::with_capacity(inventory.len());
    let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
    for (i, (&energy, cell)) in energies.iter().zip(inventory).enumerate() {
        let [lo, hi] = cell.specific_energy_bounds_j_kg;
        if !energy.is_finite() || !(lo..=hi).contains(&energy) {
            return Err(CaloricEquilibriumError::Unrepresentable(
                "final original energy coordinate",
            ));
        }
        let temperature = queries.get(i, energy)?.state.temperature_k;
        minimum = minimum.min(temperature);
        maximum = maximum.max(temperature);
        temperatures.push(temperature);
    }
    let spread = maximum - minimum;
    if residual.abs() > controls.energy_tolerance_j || spread > controls.temperature_tolerance_k {
        return Err(CaloricEquilibriumError::NotConverged);
    }
    queries.poll()?;
    Ok(CaloricEquilibrium {
        specific_energies_j_kg: energies,
        temperatures_k: temperatures,
        energy_residual_j: residual,
        temperature_spread_k: spread,
        temperature_iterations: iterations,
        caloric_evaluations: queries.used,
        retained_initial_state: retained,
    })
}

struct Probe {
    inverse: Vec<Inverse>,
    lower_energy: Vec<f64>,
    upper_energy: Vec<f64>,
    lower_residual: Dd,
    upper_residual: Dd,
    hint: Option<f64>,
}

fn probe<E, F, C>(
    inventory: &[CaloricInventory],
    bounds: &[[Sample; 2]],
    temperature: f64,
    controls: CaloricEquilibriumControls,
    queries: &mut Queries<F, C>,
) -> Result<Probe, CaloricEquilibriumError<E>>
where
    F: FnMut(usize, f64) -> Result<CaloricState, E>,
    C: FnMut() -> Result<(), E>,
{
    let mut inverse = Vec::with_capacity(inventory.len());
    let mut hint = None;
    for (i, &range) in bounds.iter().enumerate() {
        let value = invert(i, temperature, range, controls, queries)?;
        hint = hint.or(value.hint);
        inverse.push(value);
    }
    let lower_energy: Vec<_> = inverse.iter().map(|v| v.lower).collect();
    let upper_energy: Vec<_> = inverse.iter().map(|v| v.upper).collect();
    Ok(Probe {
        lower_residual: energy_change(inventory, &lower_energy)?,
        upper_residual: energy_change(inventory, &upper_energy)?,
        inverse,
        lower_energy,
        upper_energy,
        hint,
    })
}

fn allocation<E>(
    inventory: &[CaloricInventory],
    probe: &Probe,
    tolerance: f64,
) -> Result<Option<Vec<f64>>, CaloricEquilibriumError<E>> {
    if !probe.lower_residual.lt(Dd::ZERO) || !Dd::ZERO.lt(probe.upper_residual) {
        if probe.lower_residual.abs().to_f64() <= tolerance {
            return Ok(Some(probe.lower_energy.clone()));
        }
        if probe.upper_residual.abs().to_f64() <= tolerance {
            return Ok(Some(probe.upper_energy.clone()));
        }
        return Ok(None);
    }
    // An energy tolerance cannot choose a phase allocation. A strictly
    // interior latent inventory still needs its conservation-determined value
    // or an explicit refusal, even when an endpoint lies within that tolerance.
    let latent: Vec<_> = probe
        .inverse
        .iter()
        .enumerate()
        .filter_map(|(i, v)| v.latent.then_some(i))
        .collect();
    if latent.len() != 1 {
        return Err(CaloricEquilibriumError::UnderdeterminedPhaseAllocation {
            cells: latent.len(),
        });
    }
    let index = latent[0];
    let mut energies = probe.lower_energy.clone();
    // Preserve every other inverse and solve the remaining ORIGINAL energy row.
    // No material-independent phase fraction or distribution rule is selected.
    energies[index] = inventory[index].initial_specific_energy_j_kg;
    let others = energy_change::<E>(inventory, &energies)?;
    let required = finite_dd(
        Dd::from_f64(inventory[index].initial_specific_energy_j_kg)
            - others / Dd::from_f64(inventory[index].mass_kg),
    )?
    .to_f64();
    if !(probe.inverse[index].lower..=probe.inverse[index].upper).contains(&required) {
        return Err(CaloricEquilibriumError::Unrepresentable(
            "conserved latent allocation",
        ));
    }
    energies[index] = required;
    Ok(Some(energies))
}

/// Solve one thermally connected closed component using its actual T(e) laws.
///
/// Source calls are immutable and may be repeated in arbitrary order. Both
/// generalized inverses are bracketed, so sharp plateaus need no artificial
/// heat capacity. The outer scalar solve uses original sensible slopes when
/// available and otherwise bisects. Every original query shares one budget.
/// No trajectory, pressure-work conversion, entropy amount or reachability is
/// inferred. A law's own domain/refusal and each checkpoint error are retained.
pub fn equilibrate_caloric_inventory<E>(
    inventory: &[CaloricInventory],
    controls: CaloricEquilibriumControls,
    source: impl FnMut(usize, f64) -> Result<CaloricState, E>,
    checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<CaloricEquilibrium, CaloricEquilibriumError<E>> {
    if inventory.is_empty()
        || inventory.len() > 1_048_576
        || !controls.temperature_tolerance_k.is_finite()
        || controls.temperature_tolerance_k < 0.0
        || !controls.energy_tolerance_j.is_finite()
        || controls.energy_tolerance_j < 0.0
        || !(1..=256).contains(&controls.maximum_temperature_iterations)
        || !(1..=2048).contains(&controls.maximum_inverse_iterations)
        || controls.maximum_caloric_evaluations == 0
    {
        return Err(CaloricEquilibriumError::Invalid(
            "invalid inventory size or solve controls",
        ));
    }
    let mut queries = Queries {
        source,
        checkpoint,
        used: 0,
        maximum: controls.maximum_caloric_evaluations,
    };
    let mut initial = Vec::with_capacity(inventory.len());
    let mut energies = Vec::with_capacity(inventory.len());
    let (mut low_temperature, mut high_temperature) = (f64::INFINITY, 0.0_f64);
    for (i, cell) in inventory.iter().enumerate() {
        queries.poll()?;
        let [lo, hi] = cell.specific_energy_bounds_j_kg;
        let energy = cell.initial_specific_energy_j_kg;
        if !cell.mass_kg.is_finite()
            || cell.mass_kg <= 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || !energy.is_finite()
            || !(lo..=hi).contains(&energy)
        {
            return Err(CaloricEquilibriumError::Invalid(
                "invalid original caloric parcel",
            ));
        }
        let state = queries.get(i, energy)?;
        low_temperature = low_temperature.min(state.state.temperature_k);
        high_temperature = high_temperature.max(state.state.temperature_k);
        initial.push(state);
        energies.push(energy);
    }
    if high_temperature - low_temperature <= controls.temperature_tolerance_k {
        return finish(inventory, energies, controls, 0, true, &mut queries);
    }
    let mut bounds = Vec::with_capacity(inventory.len());
    for (i, cell) in inventory.iter().enumerate() {
        let [lo, hi] = cell.specific_energy_bounds_j_kg;
        let lower = queries.get(i, lo)?;
        let upper = queries.get(i, hi)?;
        if lower.state.temperature_k > initial[i].state.temperature_k
            || upper.state.temperature_k < initial[i].state.temperature_k
        {
            return Err(CaloricEquilibriumError::Invalid(
                "original caloric domain is not monotone",
            ));
        }
        low_temperature = low_temperature.max(lower.state.temperature_k);
        high_temperature = high_temperature.min(upper.state.temperature_k);
        bounds.push([lower, upper]);
    }
    if low_temperature > high_temperature {
        return Err(CaloricEquilibriumError::NoCommonTemperature);
    }
    let low = probe(inventory, &bounds, low_temperature, controls, &mut queries)?;
    if let Some(energies) = allocation(inventory, &low, controls.energy_tolerance_j)? {
        return finish(inventory, energies, controls, 0, false, &mut queries);
    }
    let high = probe(inventory, &bounds, high_temperature, controls, &mut queries)?;
    if let Some(energies) = allocation(inventory, &high, controls.energy_tolerance_j)? {
        return finish(inventory, energies, controls, 0, false, &mut queries);
    }
    if !low.upper_residual.lt(Dd::ZERO) || !Dd::ZERO.lt(high.lower_residual) {
        return Err(CaloricEquilibriumError::NoConservedEnergyRoot);
    }
    let mut temperature = low_temperature.midpoint(high_temperature);
    let mut hint = low.hint.or(high.hint);
    for iteration in 1..=controls.maximum_temperature_iterations {
        queries.poll()?;
        if let Some(t) = hint.filter(|t| *t > low_temperature && *t < high_temperature) {
            temperature = t;
        }
        if temperature <= low_temperature || temperature >= high_temperature {
            return Err(CaloricEquilibriumError::NotConverged);
        }
        let current = probe(inventory, &bounds, temperature, controls, &mut queries)?;
        if let Some(energies) = allocation(inventory, &current, controls.energy_tolerance_j)? {
            return finish(
                inventory,
                energies,
                controls,
                iteration,
                false,
                &mut queries,
            );
        }
        if current.upper_residual.lt(Dd::ZERO) {
            low_temperature = temperature;
        } else {
            high_temperature = temperature;
        }
        let mut capacity = Dd::ZERO;
        for (cell, inverse) in inventory.iter().zip(&current.inverse) {
            if inverse.slope == 0.0 {
                capacity = Dd::ZERO;
                break;
            }
            capacity = capacity + Dd::from_f64(cell.mass_kg) / Dd::from_f64(inverse.slope);
        }
        let correction = current.lower_residual.to_f64() / capacity.to_f64();
        let newton = temperature - correction;
        // Periodic bisection guarantees bracket progress without rejecting
        // useful Newton corrections merely because the root is near an end.
        temperature = if iteration % 4 != 0
            && newton.is_finite()
            && newton > low_temperature
            && newton < high_temperature
        {
            newton
        } else {
            low_temperature.midpoint(high_temperature)
        };
        hint = current.hint;
    }
    Err(CaloricEquilibriumError::NotConverged)
}

#[cfg(test)]
#[path = "caloric_tests.rs"]
mod tests;
