//! Finite-bath, spherical shrinking-core leaching with resistances in series.
//!
//! The outer particle radius R remains fixed while the reactive core recedes.
//! For the declared equation a A(aq) + b B(s) -> products, reagent flux is
//! first order in A. With u = 1-r/R, X = u(3-3u+u²), the governing equation is
//!
//! dt/du = (a/b) rho_B R [r_f (1-u)² + r_s + R r_d u(1-u)] / C_A(u),
//! C_A(u) = [n_A0 - (a/b) n_B0 X(u)] / V.
//!
//! rho_B is reactive-solid molar density, r_f=1/k_f and r_s=1/k_s are in
//! s/m, and r_d=1/D_e is in s/m². A zero resistance explicitly neglects that
//! mechanism; it does not mean zero transfer. At least one must be positive.
//! Equation (3) of Shiea, Crema & Macchi, Chemical Engineering Science 292
//! (2024) 119987, <https://doi.org/10.1016/j.ces.2024.119987>, gives the
//! diameter form. The equation above follows by l=2r, L=2R and the finite
//! well-mixed reagent balance, without assuming a rate-limiting mechanism.
//!
//! One initially unreacted monodisperse mineral population, fixed liquid
//! volume/temperature/transport coefficients, no inflow, product inhibition,
//! equilibrium, pore evolution, fragmentation, or activity/speciation solve.
//! The source supplies coefficients and applicability, never the mineral name.
//! In particular, this is not an elementary mechanism for every balanced
//! mineral equation. Constant-concentration calculations are exposed separately.
//! Numerical quadrature errors are estimates, not interval certificates or
//! evidence of ore calibration. Finite reagent depletion is integrated, not
//! applied as a clip to an excess-reagent recovery curve.

/// Persistent mineral/size populations sharing one finite reagent bath.
#[path = "shrinking_core/population.rs"]
pub mod population;

use super::{ReactantAccount, stock};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Length, Qty, QuantitySpec, Temperature, Time,
    Volume,
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{ConservedReaction, ReactionAmount, ReactionBalanceError},
};

/// Source model identifier; version 1 fixes spherical, first-order transport.
pub const LAW: &str = "spherical-shrinking-core-reaction";
/// Reciprocal film mass-transfer coefficient, in s/m.
pub const FILM_RESISTANCE: &str = "film-transfer-resistance";
/// Reciprocal surface-reaction velocity coefficient, in s/m.
pub const SURFACE_RESISTANCE: &str = "surface-reaction-resistance";
/// Reciprocal effective product-layer diffusivity, in s/m².
pub const DIFFUSION_RESISTANCE: &str = "product-layer-diffusion-resistance";
/// Moles of reactive mineral per cubic metre of the original solid core.
pub const SOLID_MOLAR_DENSITY: &str = "solid-molar-density";
/// Initial outer radius of each particle, in metres (not diameter).
pub const RADIUS_AXIS: &str = "initial-particle-radius";
/// Actual well-mixed aqueous reagent concentration, in mol/m³.
pub const CONCENTRATION_AXIS: &str = "reagent-concentration";
/// Original solid inventory divided by liquid volume; remains a batch descriptor.
pub const LOADING_AXIS: &str = "initial-solid-loading";

use super::{ELAPSED_TIME_AXIS, TEMPERATURE_AXIS, VOLUME_AXIS};
/// Reciprocal transfer velocity in coherent SI s/m.
pub type TransferResistance = Qty<-1, 0, 1, 0, 0, 0>;
/// Reciprocal diffusivity in coherent SI s/m².
pub type DiffusionResistance = Qty<-2, 0, 1, 0, 0, 0>;

const AXES: [(&str, Dims); 6] = [
    (ELAPSED_TIME_AXIS, Time::DIMS),
    (TEMPERATURE_AXIS, Temperature::DIMS),
    (VOLUME_AXIS, Volume::DIMS),
    (RADIUS_AXIS, Length::DIMS),
    (CONCENTRATION_AXIS, AmountConcentration::DIMS),
    (LOADING_AXIS, AmountConcentration::DIMS),
];

/// Source, physical-input or bounded-numerical refusal; inputs are never mutated.
#[derive(Clone, Debug, PartialEq)]
pub enum ShrinkingCoreError {
    /// Original material-card or typed-query admission failure.
    Source(MatDbError),
    /// Original exact-reaction or finite-inventory admission failure.
    Reaction(ReactionBalanceError),
    /// Invalid model, physical quantity or species-role association.
    Invalid(&'static str),
    /// A required source coordinate is absent, mistyped or outside its support.
    OutsideValidity,
    /// A positive physical quantity or progress cannot be represented.
    Arithmetic(&'static str),
    /// Requested tolerance exceeded the caller's work or subdivision budget.
    Budget,
    /// Caller cancelled before any result was published.
    Cancelled,
}
impl core::fmt::Display for ShrinkingCoreError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shrinking-core leach: {self:?}")
    }
}
impl std::error::Error for ShrinkingCoreError {}
impl From<MatDbError> for ShrinkingCoreError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e)
    }
}
impl From<ReactionBalanceError> for ShrinkingCoreError {
    fn from(e: ReactionBalanceError) -> Self {
        Self::Reaction(e)
    }
}
type Result<T> = core::result::Result<T, ShrinkingCoreError>;

/// Explicit accuracy and work limits for scalar quadrature and inversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LeachNumerics {
    /// Positive absolute time residual allowance, in seconds.
    pub absolute_time_tolerance: Time,
    /// Relative time residual allowance, in (0, 0.01].
    pub relative_time_tolerance: f64,
    /// Shared integrand-evaluation limit across all root iterations, at most 1e6.
    pub max_evaluations: usize,
    /// Bisection limit, between one and 128.
    pub max_iterations: usize,
}
impl Default for LeachNumerics {
    fn default() -> Self {
        Self {
            absolute_time_tolerance: Time::new(1e-9),
            relative_time_tolerance: 1e-9,
            max_evaluations: 65_536,
            max_iterations: 96,
        }
    }
}
impl LeachNumerics {
    fn tolerance(self, t: f64) -> Result<f64> {
        if !(self.absolute_time_tolerance.value().is_finite()
            && self.absolute_time_tolerance.value() > 0.0
            && self.relative_time_tolerance.is_finite()
            && self.relative_time_tolerance > 0.0
            && self.relative_time_tolerance <= 0.01
            && (5..=1_000_000).contains(&self.max_evaluations)
            && (1..=128).contains(&self.max_iterations))
        {
            return Err(ShrinkingCoreError::Invalid(
                "invalid time tolerance or numerical work limit",
            ));
        }
        positive(
            self.absolute_time_tolerance.value() + self.relative_time_tolerance * t,
            "time tolerance range loss",
        )
    }
}

/// Immutable source-bound mineral/reagent association and series resistances.
#[derive(Clone, Debug)]
pub struct ShrinkingCoreReaction {
    reaction: ConservedReaction,
    source: ConstitutiveModelCard,
    solid: SpeciesId,
    reagent: SpeciesId,
    coefficients: [i128; 2], // solid, reagent
    resistances: [f64; 3],   // film, surface, product-layer
    molar_density: AmountConcentration,
}

/// Completed finite-bath calculation, with original source and species accounts.
#[derive(Clone, Debug)]
pub struct ShrinkingCoreAdvance {
    model: ShrinkingCoreReaction,
    /// Original complete canonical species inventory, in mol.
    pub initial_inventory: Vec<(SpeciesId, Amount)>,
    /// Actual complete final species inventory from the conserved-reaction owner.
    pub final_inventory: Vec<(SpeciesId, Amount)>,
    /// Extent in moles of the exact integer equation as written.
    pub extent: Amount,
    /// Mineral conversion fraction, distinct from metal recovery or payability.
    pub conversion: Dimensionless,
    /// Radius of the remaining reactive core; the outer residue radius stays fixed.
    pub core_radius: Length,
    /// Gross consumption/formation from the original stoichiometry.
    pub amounts: Vec<ReactionAmount>,
    /// Solid and aqueous-reagent stock accounts, in that order.
    pub reactants: [ReactantAccount; 2],
    /// Initial source conditions, including original solid loading.
    pub initial_point: QueryPoint,
    /// Final source conditions with actual depleted reagent concentration.
    pub final_point: QueryPoint,
    /// Requested elapsed time in seconds, including any time after completion.
    pub elapsed: Time,
    /// Absolute inversion residual plus quadrature-error estimate, in seconds.
    /// Zero after admitted complete dissolution; not a certified error bound.
    pub estimated_time_error: Time,
    /// Actual integrand calls, bounded by the caller's shared work limit.
    pub evaluations: usize,
}
impl ShrinkingCoreAdvance {
    /// Original immutable equation, source hashes, coefficients and applicability.
    pub fn model(&self) -> &ShrinkingCoreReaction {
        &self.model
    }
}

impl ShrinkingCoreReaction {
    /// Admit exactly two consumed species with caller-declared solid/reagent roles.
    /// The source has four parameters: three nonnegative reciprocal transport
    /// coefficients and positive reactive-solid molar density. Source provenance
    /// and all six typed applicability axes are mandatory. No rate is inferred.
    pub fn from_card(
        reaction: ConservedReaction,
        solid: SpeciesId,
        reagent: SpeciesId,
        source: ConstitutiveModelCard,
    ) -> Result<Self> {
        source.validate()?;
        if source.law.0 != LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.parameters.len() != 4
            || source.sources.is_empty()
            || source.provenance.source.trim().is_empty()
            || source.provenance.license.trim().is_empty()
        {
            return Err(ShrinkingCoreError::Invalid(
                "requires version 1 source, zero hidden state and four sourced parameters",
            ));
        }
        let mut parameters = [0.0; 4];
        for (i, (name, dims)) in [
            (FILM_RESISTANCE, TransferResistance::DIMS),
            (SURFACE_RESISTANCE, TransferResistance::DIMS),
            (DIFFUSION_RESISTANCE, DiffusionResistance::DIMS),
            (SOLID_MOLAR_DENSITY, AmountConcentration::DIMS),
        ]
        .iter()
        .enumerate()
        {
            let p = source
                .parameters
                .get(*name)
                .ok_or(ShrinkingCoreError::Invalid(
                    "missing shrinking-core parameter",
                ))?;
            if p.dims != *dims || !p.value.is_finite() || p.value < 0.0 {
                return Err(ShrinkingCoreError::Invalid(
                    "parameter dimensions or sign are invalid",
                ));
            }
            parameters[i] = p.value;
        }
        if parameters[3] <= 0.0 || parameters[..3].iter().all(|r| *r == 0.0) {
            return Err(ShrinkingCoreError::Invalid(
                "solid molar density and at least one resistance must be positive",
            ));
        }
        for (axis, (lo, hi)) in source.validity.bounds() {
            if !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
                || !source.validity.axis_quantities().contains_key(axis)
            {
                return Err(ShrinkingCoreError::Invalid(
                    "all source axes require finite typed bounds",
                ));
            }
        }
        for (axis, dims) in AXES {
            let (lo, hi) = source
                .validity
                .bound(axis)
                .ok_or(ShrinkingCoreError::Invalid("missing owned source axis"))?;
            if source.validity.axis_quantities().get(axis) != Some(&QuantitySpec::dimensional(dims))
                || lo < 0.0
                || (axis == ELAPSED_TIME_AXIS && (lo != 0.0 || hi <= 0.0))
                || ([TEMPERATURE_AXIS, VOLUME_AXIS, RADIUS_AXIS].contains(&axis) && lo <= 0.0)
            {
                return Err(ShrinkingCoreError::Invalid(
                    "source axes need canonical units and physically valid bounds",
                ));
            }
        }
        if solid == reagent
            || reaction
                .terms()
                .iter()
                .filter(|t| t.coefficient < 0)
                .count()
                != 2
        {
            return Err(ShrinkingCoreError::Invalid(
                "declare two distinct consumed solid and reagent species",
            ));
        }
        let mut coefficients = [0; 2];
        for (i, id) in [&solid, &reagent].iter().enumerate() {
            coefficients[i] = -reaction
                .terms()
                .iter()
                .find(|t| &t.species == *id && t.coefficient < 0)
                .ok_or(ShrinkingCoreError::Invalid(
                    "mapped species is not a reactant",
                ))?
                .coefficient;
        }
        Ok(Self {
            reaction,
            source,
            solid,
            reagent,
            coefficients,
            resistances: [parameters[0], parameters[1], parameters[2]],
            molar_density: AmountConcentration::new(parameters[3]),
        })
    }

    /// Original exact equation and elemental/charge conservation certificate.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Complete original source card, with fit evidence and validity unchanged.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Caller-declared solid then reagent roles, without parsing chemical labels.
    pub fn reactants(&self) -> (&SpeciesId, &SpeciesId) {
        (&self.solid, &self.reagent)
    }

    /// Analytical elapsed time at a prescribed penetration 1-r/R under a
    /// separately maintained constant bulk concentration. This is NOT the
    /// finite-bath solution. Its full time/source domain is checked. The loading
    /// and volume remain explicit context, not inferred from particle radius.
    pub fn constant_concentration_time(
        &self,
        penetration: Dimensionless,
        radius: Length,
        concentration: AmountConcentration,
        volume: Volume,
        initial_solid_loading: AmountConcentration,
        temperature: Temperature,
        conditions: &QueryPoint,
    ) -> Result<Time> {
        let u = penetration.value();
        if !u.is_finite() || !(0.0..=1.0).contains(&u) {
            return Err(ShrinkingCoreError::Invalid("penetration must lie in [0,1]"));
        }
        positive(
            concentration.value(),
            "constant concentration must be positive",
        )?;
        let scale = self.scale(radius)?;
        let [rf, rs, rd] = self.resistances;
        let time = scale
            * (rf * conversion(u) / 3.0
                + rs * u
                + radius.value() * rd * u * u * (3.0 - 2.0 * u) / 6.0)
            / concentration.value();
        nonnegative(time, "constant-concentration time range loss")?;
        if u > 0.0 && time == 0.0 {
            return Err(ShrinkingCoreError::Arithmetic(
                "positive conversion time underflow",
            ));
        }
        for elapsed in [Time::new(0.0), Time::new(time)] {
            self.point(
                volume,
                radius,
                concentration,
                initial_solid_loading,
                elapsed,
                temperature,
                conditions,
            )?;
        }
        Ok(Time::new(time))
    }

    /// Advance an initially unreacted, monodisperse population in a closed bath.
    /// Both reactant stocks deplete according to the original equation. The
    /// liquid volume, outer radius and temperature stay fixed. Additional source
    /// conditions remain mandatory; callers cannot replace the six owned axes.
    /// Adaptive Simpson quadrature and bisection share an explicit work budget.
    /// Cancellation is polled on every numerical evaluation and at publication.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        inventory: &[(SpeciesId, Amount)],
        volume: Volume,
        radius: Length,
        elapsed: Time,
        temperature: Temperature,
        conditions: &QueryPoint,
        numerics: LeachNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShrinkingCoreAdvance> {
        if cancelled() {
            return Err(ShrinkingCoreError::Cancelled);
        }
        nonnegative(
            elapsed.value(),
            "elapsed time must be finite and nonnegative",
        )?;
        let tolerance = numerics.tolerance(elapsed.value())?;
        positive(volume.value(), "liquid volume must be positive")?;
        let scale = self.scale(radius)?;
        let initial = self
            .reaction
            .advance_inventory(inventory, Amount::new(0.0))?;
        let solid0 = stock(&initial, &self.solid);
        let reagent0 = stock(&initial, &self.reagent);
        let loading = quotient(solid0, volume.value())?;
        let concentration = quotient(reagent0, volume.value())?;
        let initial_point = self.point(
            volume,
            radius,
            AmountConcentration::new(concentration),
            AmountConcentration::new(loading),
            Time::new(0.0),
            temperature,
            conditions,
        )?;
        // The whole requested time is admitted, including a possible post-dissolution hold.
        self.point(
            volume,
            radius,
            AmountConcentration::new(concentration),
            AmountConcentration::new(loading),
            elapsed,
            temperature,
            conditions,
        )?;
        let mut evaluations = 0;
        let (u, time_error) = if elapsed.value() == 0.0 || solid0 == 0.0 || reagent0 == 0.0 {
            (0.0, 0.0)
        } else {
            let demand = positive(
                (self.coefficients[1] as f64 / self.coefficients[0] as f64) * loading,
                "reagent concentration demand range loss",
            )?;
            let kernel = Kernel {
                scale,
                film: self.resistances[0],
                surface: self.resistances[1],
                diffusion: radius.value() * self.resistances[2],
                concentration,
                demand,
            };
            nonnegative(kernel.diffusion, "diffusion resistance scale range loss")?;
            solve(
                &kernel,
                elapsed.value(),
                tolerance,
                numerics,
                &mut evaluations,
                &mut cancelled,
            )?
        };
        let x = conversion(u);
        let extent_value = solid0 * x / self.coefficients[0] as f64;
        nonnegative(extent_value, "extent range loss")?;
        if u > 0.0 && extent_value == 0.0 {
            return Err(ShrinkingCoreError::Arithmetic("positive extent underflow"));
        }
        let extent = Amount::new(extent_value);
        let final_inventory = self.reaction.advance_inventory(&initial, extent)?;
        let remaining = [
            stock(&final_inventory, &self.solid),
            stock(&final_inventory, &self.reagent),
        ];
        if u > 0.0 && ((u < 1.0 && remaining[0] == 0.0) || remaining[1] == 0.0) {
            return Err(ShrinkingCoreError::Arithmetic(
                "finite-time progress erased positive remaining stock",
            ));
        }
        if u == 1.0 && remaining[0] != 0.0 {
            return Err(ShrinkingCoreError::Arithmetic(
                "complete dissolution extent is not exactly representable",
            ));
        }
        let final_point = self.point(
            volume,
            radius,
            AmountConcentration::new(quotient(remaining[1], volume.value())?),
            AmountConcentration::new(loading),
            elapsed,
            temperature,
            conditions,
        )?;
        let core_radius = radius.value() * (1.0 - u);
        if u < 1.0 && core_radius == 0.0 {
            return Err(ShrinkingCoreError::Arithmetic(
                "remaining core radius underflow",
            ));
        }
        let reactants = core::array::from_fn(|i| ReactantAccount {
            species: if i == 0 {
                self.solid.clone()
            } else {
                self.reagent.clone()
            },
            coefficient: self.coefficients[i],
            initial: Amount::new(if i == 0 { solid0 } else { reagent0 }),
            consumed: Amount::new(self.coefficients[i] as f64 * extent_value),
            remaining: Amount::new(remaining[i]),
        });
        if cancelled() {
            return Err(ShrinkingCoreError::Cancelled);
        }
        Ok(ShrinkingCoreAdvance {
            model: self.clone(),
            initial_inventory: initial,
            final_inventory,
            extent,
            conversion: Dimensionless::new(x),
            core_radius: Length::new(core_radius),
            amounts: self.reaction.amounts(extent)?,
            reactants,
            initial_point,
            final_point,
            elapsed,
            estimated_time_error: Time::new(time_error),
            evaluations,
        })
    }

    fn scale(&self, radius: Length) -> Result<f64> {
        positive(radius.value(), "initial radius must be finite and positive")?;
        positive(
            (self.coefficients[1] as f64 / self.coefficients[0] as f64)
                * self.molar_density.value()
                * radius.value(),
            "particle time scale range loss",
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn point(
        &self,
        volume: Volume,
        radius: Length,
        concentration: AmountConcentration,
        loading: AmountConcentration,
        elapsed: Time,
        temperature: Temperature,
        conditions: &QueryPoint,
    ) -> Result<QueryPoint> {
        if conditions
            .axes()
            .keys()
            .any(|axis| AXES.iter().any(|(own, _)| axis == own))
        {
            return Err(ShrinkingCoreError::Invalid(
                "caller conditions must not replace owned physical axes",
            ));
        }
        let values = [
            elapsed.value(),
            temperature.value(),
            volume.value(),
            radius.value(),
            concentration.value(),
            loading.value(),
        ];
        let mut point = conditions.clone();
        for ((axis, dims), value) in AXES.into_iter().zip(values) {
            point = point.with_quantity(axis, QuantitySpec::dimensional(dims), value)?;
        }
        if !self
            .source
            .validity
            .contains_typed(point.axes(), point.axis_quantities())
        {
            return Err(ShrinkingCoreError::OutsideValidity);
        }
        Ok(point)
    }
}

// Penetration form retains tiny positive conversions without 1-(1-u)^3 cancellation.
fn conversion(u: f64) -> f64 {
    u * (3.0 - u * (3.0 - u))
}
fn nonnegative(x: f64, what: &'static str) -> Result<f64> {
    if x.is_finite() && x >= 0.0 {
        Ok(x)
    } else {
        Err(ShrinkingCoreError::Arithmetic(what))
    }
}
fn positive(x: f64, what: &'static str) -> Result<f64> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(ShrinkingCoreError::Arithmetic(what))
    }
}
fn quotient(n: f64, v: f64) -> Result<f64> {
    let x = nonnegative(n / v, "concentration range loss")?;
    if n > 0.0 && x == 0.0 {
        return Err(ShrinkingCoreError::Arithmetic(
            "positive concentration underflow",
        ));
    }
    Ok(x)
}

struct Kernel {
    scale: f64,
    film: f64,
    surface: f64,
    diffusion: f64,
    concentration: f64,
    demand: f64,
}
impl Kernel {
    fn concentration(&self, u: f64) -> f64 {
        (-self.demand).mul_add(conversion(u), self.concentration)
    }
    fn value(&self, u: f64) -> Result<f64> {
        let c = positive(
            self.concentration(u),
            "depleted reagent is an asymptotic boundary",
        )?;
        let y = 1.0 - u;
        nonnegative(
            self.scale * (self.film * y * y + self.surface + self.diffusion * u * y) / c,
            "quadrature integrand range loss",
        )
    }
}

fn solve(
    kernel: &Kernel,
    target: f64,
    tolerance: f64,
    settings: LeachNumerics,
    calls: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64)> {
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..settings.max_iterations {
        if cancelled() {
            return Err(ShrinkingCoreError::Cancelled);
        }
        let mid = lo + (hi - lo) * 0.5;
        if mid == lo || mid == hi {
            return Err(ShrinkingCoreError::Arithmetic(
                "penetration resolution exhausted",
            ));
        }
        if kernel.concentration(mid) <= 0.0 {
            hi = mid;
            continue;
        }
        let (time, error) = integrate(
            kernel,
            mid,
            tolerance / 8.0,
            settings.max_evaluations,
            calls,
            cancelled,
        )?;
        if (time - target).abs() + error <= tolerance {
            return Ok((mid, (time - target).abs() + error));
        }
        if time > target {
            hi = mid;
        } else {
            // Do not integrate the near-depletion endpoint just to answer an
            // early-time request. With excess reagent, the remaining integral
            // is bounded above by holding concentration at its final minimum.
            // Express the polynomial in remaining radius to avoid cancellation.
            let final_c = kernel.concentration(1.0);
            if final_c > 0.0 {
                let y = 1.0 - mid;
                let tail = kernel.scale
                    * (kernel.film * y * y * y / 3.0
                        + kernel.surface * y
                        + kernel.diffusion * y * y * (0.5 - y / 3.0))
                    / final_c;
                if time + error + tail <= target {
                    return Ok((1.0, 0.0));
                }
            }
            lo = mid;
        }
    }
    Err(ShrinkingCoreError::Budget)
}

// Adaptive Simpson estimates are used for acceptance, not advertised as rigorous
// enclosures. Shared call/depth limits make refusal/cancellation bounded.
fn integrate(
    kernel: &Kernel,
    upper: f64,
    tolerance: f64,
    max_calls: usize,
    calls: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64)> {
    let mut eval = |u| {
        if cancelled() {
            return Err(ShrinkingCoreError::Cancelled);
        }
        if *calls >= max_calls {
            return Err(ShrinkingCoreError::Budget);
        }
        *calls += 1;
        kernel.value(u)
    };
    let f0 = eval(0.0)?;
    let fm = eval(upper * 0.5)?;
    let f1 = eval(upper)?;
    let whole = upper * (f0 / 6.0 + 2.0 * fm / 3.0 + f1 / 6.0);
    let mut stack = vec![(0.0, upper, f0, fm, f1, whole, tolerance, 0u32)];
    let (mut total, mut error_sum) = (0.0, 0.0);
    while let Some((a, b, fa, fm, fb, whole, tol, depth)) = stack.pop() {
        let m = a + (b - a) * 0.5;
        let left_mid = a + (m - a) * 0.5;
        let right_mid = m + (b - m) * 0.5;
        if left_mid == a || right_mid == b {
            return Err(ShrinkingCoreError::Budget);
        }
        let fl = eval(left_mid)?;
        let fr = eval(right_mid)?;
        let left = (m - a) * (fa / 6.0 + 2.0 * fl / 3.0 + fm / 6.0);
        let right = (b - m) * (fm / 6.0 + 2.0 * fr / 3.0 + fb / 6.0);
        let delta = left + right - whole;
        nonnegative(left + right, "integrated time overflow")?;
        let error = delta.abs() / 15.0;
        if error <= tol {
            total += left + right + delta / 15.0;
            error_sum += error;
        } else {
            if depth == 48 {
                return Err(ShrinkingCoreError::Budget);
            }
            stack.push((m, b, fm, fr, fb, right, tol * 0.5, depth + 1));
            stack.push((a, m, fa, fl, fm, left, tol * 0.5, depth + 1));
        }
    }
    Ok((
        positive(total, "integrated positive time range loss")?,
        nonnegative(error_sum, "quadrature error range loss")?,
    ))
}

#[cfg(test)]
#[path = "shrinking_core/tests.rs"]
mod tests;
