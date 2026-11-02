//! Forward hold-up in a train of equal, ideal continuously stirred tanks.
//!
//! Each tank obeys `dC_i/dt = Q/(V/N) * (C_(i-1)-C_i)`. The caller supplies
//! active mixed volume and actual flow through the train, including its own
//! treatment of bypasses. All tanks stay full; inlet and outlet flow are equal.
//! Concentration is conserved tracer/metal mass per slurry volume, not ore grade
//! or bulk slurry density. Reaction, settling, segregation and deconvolution are
//! outside this hydraulic model.
//!
//! The equal-tank model and its Erlang residence-time distribution follow
//! Fogler, *Elements of Chemical Reaction Engineering*, §17.4, eqs.17-4/17-5:
//! <https://websites.umich.edu/~elements/5e/17chap/Ch17_Web-Additional%20Material.pdf>.
//! Model geometry and initial conditions still need the caller's provenance;
//! mathematical agreement is not a calibration claim for a real plant.

use fs_matdb::Provenance;
use fs_math::det;
use fs_qty::{Density, Mass, Time, Volume, VolumetricFlowRate};

/// Maximum explicitly represented tanks in this bounded synchronous evaluator.
pub const MAX_TANKS: usize = 64;
/// Maximum `Q * duration / (active_volume / tank_count)` per call.
/// Longer intervals can be advanced through several calls with the same input.
pub const MAX_INTERVAL_EXPOSURE: f64 = 4096.0;
const STEP_EXPOSURE: f64 = 8.0;
const POISSON_TERMS: usize = 256;
const CONSERVATION_TOLERANCE: f64 = 2.0e-11;

/// Invalid physical inputs, excessive interval work, or numerical range loss.
#[derive(Clone, Debug, PartialEq)]
pub enum ResidenceError {
    /// A physical input or its provenance is absent or invalid.
    Invalid(&'static str),
    /// A requested positive quantity cannot be represented by this evaluator.
    Arithmetic(&'static str),
    /// Split this interval into smaller calls; the state was not advanced.
    WorkLimit,
    /// Independently integrated outlet and physical inventory failed closure.
    Conservation,
}
impl core::fmt::Display for ResidenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "stirred-tank residence: {self:?}")
    }
}
impl std::error::Error for ResidenceError {}

/// One piecewise-constant inlet/flow interval, including physical mass hold-up.
#[derive(Clone, Debug, PartialEq)]
pub struct TankInterval {
    /// Elapsed interval, in seconds.
    pub duration: Time,
    /// Actual train flow, in cubic metres per second.
    pub train_flow: VolumetricFlowRate,
    /// Conserved constituent mass per cubic metre entering the train.
    pub inlet_concentration: Density,
    /// Slurry volume crossing each train boundary during this interval.
    pub exchanged_volume: Volume,
    /// Constituent mass that entered during this interval.
    pub inlet_mass: Mass,
    /// Constituent mass independently integrated at the final tank outlet.
    pub outlet_mass: Mass,
    /// Flow-weighted interval outlet concentration; absent when no volume exits.
    pub outlet_mean: Option<Density>,
    /// Constituent mass held in the tanks before advancing.
    pub initial_inventory: Mass,
    /// Constituent mass held in the tanks after advancing.
    pub final_inventory: Mass,
    /// `(final - initial) + (outlet - inlet)`, without a balance plug.
    pub conservation_residual: Mass,
    /// Bounded exponential propagation steps used (zero for no exchanged flow).
    pub substeps: usize,
}

/// Physical, retained constituent inventories in equal, perfectly mixed tanks.
#[derive(Clone, Debug, PartialEq)]
pub struct StirredTankCascade {
    active_volume: Volume,
    concentrations: Vec<Density>,
    inventory: Mass,
    provenance: Provenance,
}

impl StirredTankCascade {
    /// Initialize one concentration per tank; the vector length is the tank count.
    /// `active_volume` is the combined wet, mixed volume of these tanks.
    pub fn new(
        active_volume: Volume,
        initial_concentrations: Vec<Density>,
        provenance: Provenance,
    ) -> Result<Self, ResidenceError> {
        let n = initial_concentrations.len();
        if n == 0
            || n > MAX_TANKS
            || !active_volume.value().is_finite()
            || active_volume.value() <= 0.0
            || initial_concentrations
                .iter()
                .any(|c| !c.value().is_finite() || c.value() < 0.0)
            || provenance.source.trim().is_empty()
            || provenance.license.trim().is_empty()
        {
            return Err(ResidenceError::Invalid(
                "requires 1..=64 tanks, positive finite active volume, nonnegative finite concentrations and provenance",
            ));
        }
        let cell_volume = active_volume.value() / n as f64;
        if cell_volume == 0.0 {
            return Err(ResidenceError::Arithmetic("individual tank volume"));
        }
        let inventory = inventory(cell_volume, &initial_concentrations)?;
        Ok(Self {
            active_volume,
            concentrations: initial_concentrations,
            inventory: Mass::new(inventory),
            provenance,
        })
    }

    /// Initialize every tank at the same conserved constituent concentration.
    pub fn uniform(
        tank_count: usize,
        active_volume: Volume,
        initial_concentration: Density,
        provenance: Provenance,
    ) -> Result<Self, ResidenceError> {
        if tank_count == 0 || tank_count > MAX_TANKS {
            return Err(ResidenceError::Invalid("tank count must be 1..=64"));
        }
        Self::new(
            active_volume,
            vec![initial_concentration; tank_count],
            provenance,
        )
    }

    /// Total active mixed volume, excluding any caller-declared bypass/dead space.
    pub const fn active_volume(&self) -> Volume {
        self.active_volume
    }
    /// Current concentrations in physical upstream-to-downstream order.
    pub fn tank_concentrations(&self) -> &[Density] {
        &self.concentrations
    }
    /// Total retained conserved constituent mass, in kilograms.
    pub const fn inventory(&self) -> Mass {
        self.inventory
    }
    /// Original source of geometry, initial state and model use.
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }

    /// Advance with constant feed concentration and actual train flow over this
    /// interval. Flow may change between calls, or stop at zero without losing
    /// inventory. No change is published if the interval refuses.
    ///
    /// Positive Poisson weights propagate the linear tank balances. Exposure is
    /// split into at most 512 steps of at most 8; each uses 256 terms. Outlet mass
    /// is an independently integrated output of these balances, not a residual
    /// chosen to force closure. The returned residual must close to 2e-11 of the
    /// largest inlet/outlet/initial/final mass. This is floating-point numerical
    /// evaluation, not a certified interval enclosure.
    pub fn advance(
        &mut self,
        duration: Time,
        actual_train_flow: VolumetricFlowRate,
        inlet_concentration: Density,
    ) -> Result<TankInterval, ResidenceError> {
        let (dt, q, feed) = (
            duration.value(),
            actual_train_flow.value(),
            inlet_concentration.value(),
        );
        if [dt, q, feed].iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(ResidenceError::Invalid(
                "duration, actual train flow and concentration must be finite and nonnegative",
            ));
        }
        let initial = self.inventory.value();
        let mut interval = TankInterval {
            duration,
            train_flow: actual_train_flow,
            inlet_concentration,
            exchanged_volume: Volume::new(0.0),
            inlet_mass: Mass::new(0.0),
            outlet_mass: Mass::new(0.0),
            outlet_mean: None,
            initial_inventory: self.inventory,
            final_inventory: self.inventory,
            conservation_residual: Mass::new(0.0),
            substeps: 0,
        };
        if dt == 0.0 || q == 0.0 {
            return Ok(interval);
        }
        let exchanged = product(q, dt, "exchanged volume")?;
        let n = self.concentrations.len();
        let cell_volume = self.active_volume.value() / n as f64;
        let exposure = exchanged / cell_volume;
        if !exposure.is_finite() || exposure > MAX_INTERVAL_EXPOSURE {
            return Err(ResidenceError::WorkLimit);
        }
        if exposure == 0.0 {
            return Err(ResidenceError::Arithmetic(
                "dimensionless residence exposure",
            ));
        }
        let steps = (exposure / STEP_EXPOSURE).ceil().max(1.0) as usize;
        let h = exposure / steps as f64;
        let step_volume = exchanged / steps as f64;
        let inlet_mass = product(exchanged, feed, "inlet mass")?;
        let (_, feed_integrals, feed_mean) = scaled_coefficients(h, n, feed);
        let mut state = self.concentrations.clone();
        let mut outlet_mass = 0.0;
        for _ in 0..steps {
            // ∫C_out du = Cin*H_N + Σ C_(N-1-j)(0)*T_j, where
            // T_j=Σ_(k>j)p_k and H_N=Σ_(k>N)(k-N)*p_k. Dividing by h
            // gives an interval mean without subtracting nearly equal states.
            let mut outlet_mean = feed_mean;
            let mut next: Vec<_> = feed_integrals[..n].iter().map(|value| value * h).collect();
            for (i, concentration) in state.iter().enumerate() {
                let (weights, integrals, _) = scaled_coefficients(h, n - i, concentration.value());
                outlet_mean += integrals[n - i - 1];
                for j in 0..n - i {
                    next[i + j] += weights[j];
                }
            }
            // Every positive upstream source reaches each downstream tank over
            // positive time in this ideal model. A zero on such a path denotes
            // unrepresentable underflow, not a physically empty output.
            let mut positive_path = feed > 0.0;
            for (i, concentration) in state.iter().enumerate() {
                positive_path |= concentration.value() > 0.0;
                if positive_path && next[i] == 0.0 {
                    return Err(ResidenceError::Arithmetic(
                        "positive tank response underflow",
                    ));
                }
            }
            if positive_path && outlet_mean == 0.0 {
                return Err(ResidenceError::Arithmetic(
                    "positive outlet response underflow",
                ));
            }
            outlet_mass += product(step_volume, outlet_mean, "integrated outlet mass")?;
            if next
                .iter()
                .any(|concentration| !concentration.is_finite() || *concentration < 0.0)
            {
                return Err(ResidenceError::Arithmetic("tank concentration"));
            }
            state = next.into_iter().map(Density::new).collect();
        }
        let final_inventory = inventory(cell_volume, &state)?;
        let mean = outlet_mass / exchanged;
        let residual = (final_inventory - initial) + (outlet_mass - inlet_mass);
        let scale = initial
            .max(final_inventory)
            .max(inlet_mass)
            .max(outlet_mass);
        if !outlet_mass.is_finite() || !mean.is_finite() || (outlet_mass > 0.0 && mean == 0.0) {
            return Err(ResidenceError::Arithmetic("interval outlet"));
        }
        if !residual.is_finite() || residual.abs() > CONSERVATION_TOLERANCE * scale {
            return Err(ResidenceError::Conservation);
        }
        interval.exchanged_volume = Volume::new(exchanged);
        interval.inlet_mass = Mass::new(inlet_mass);
        interval.outlet_mass = Mass::new(outlet_mass);
        interval.outlet_mean = Some(Density::new(mean));
        interval.final_inventory = Mass::new(final_inventory);
        interval.conservation_residual = Mass::new(residual);
        interval.substeps = steps;
        self.concentrations = state;
        self.inventory = interval.final_inventory;
        Ok(interval)
    }
}

fn product(a: f64, b: f64, context: &'static str) -> Result<f64, ResidenceError> {
    let value = a * b;
    if !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        Err(ResidenceError::Arithmetic(context))
    } else {
        Ok(value)
    }
}

fn inventory(cell_volume: f64, state: &[Density]) -> Result<f64, ResidenceError> {
    let mut mass = 0.0;
    for concentration in state {
        mass += product(cell_volume, concentration.value(), "tank inventory")?;
    }
    if mass.is_finite() {
        Ok(mass)
    } else {
        Err(ResidenceError::Arithmetic("total tank inventory"))
    }
}

// For h<=8, omitted Poisson terms after k=256 are far below binary64 relative
// precision. Summing tails from small to large avoids 1-CDF cancellation.
fn scaled_coefficients(
    h: f64,
    n: usize,
    concentration: f64,
) -> ([f64; MAX_TANKS], [f64; MAX_TANKS], f64) {
    // q[k-1] = C*p_k/h, initialized without dividing by h. Scaling before
    // recurrence keeps representable responses when a bare probability such
    // as h²/2 underflows (e.g. h=1e-200, C=1e200). Never form that bare term.
    let mut q = [0.0; POISSON_TERMS];
    q[0] = concentration * det::exp(-h);
    for k in 1..POISSON_TERMS {
        q[k] = q[k - 1] * (h / (k + 1) as f64);
    }
    let mut tails = [0.0; MAX_TANKS];
    let mut sum = 0.0;
    let mut integrated_feed = 0.0;
    for k in (0..POISSON_TERMS).rev() {
        sum += q[k];
        if k < n {
            tails[k] = sum;
        }
        if k >= n {
            integrated_feed += (k + 1 - n) as f64 * q[k];
        }
    }
    let mut weights = [0.0; MAX_TANKS];
    weights[0] = q[0];
    for k in 1..n {
        weights[k] = q[k - 1] * h;
    }
    (weights, tails, integrated_feed)
}

#[cfg(test)]
#[path = "residence_tests.rs"]
mod tests;
