//! Finite-bath aqueous Langmuir adsorption and desorption of one analytical metal.
//!
//! The source-declared isotherm is `q = q_max c / (c_half + c)`, where `q`
//! is mol metal/kg dry adsorbent and `c` is mol metal/m³ solution. The form
//! follows USGS PHREEQC example 19, equation 35:
//! <https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-81.htm>.
//! Parameters and chemical support must come from the supplied original card;
//! an element or carbon name supplies neither a capacity nor an affinity.
//!
//! A contact solves this isotherm together with the finite metal balance,
//! `V c + m q = n_dissolved_initial + n_adsorbed_initial`. Initially loaded
//! adsorbent may release metal. Concentration support applies to the solved
//! equilibrium, not to the nonequilibrium starting liquid. The carrier volume
//! and dry adsorbent mass are fixed; amounts refer to the same analytical
//! element in both phases. This is a dilute empirical partition model, not
//! aqueous speciation, charge balance, solvent adsorption, competitive uptake,
//! leaching, elution chemistry, a kinetic rate or a prediction of plant recovery.

/// Persistent finite-time uptake and desorption in one shared solution bath.
#[path = "aqueous_adsorption/kinetics.rs"]
pub mod kinetics;

/// Finite-time sorption into retained populations with actual solution through-flow.
#[path = "aqueous_adsorption/flow.rs"]
pub mod flow;

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::det;
use fs_qty::{
    Amount, AmountConcentration, Dims, Mass, Qty, QuantitySpec, Temperature, Volume,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

/// Moles of the declared analytical element per kilogram of dry adsorbent.
pub type AdsorbentLoading = Qty<0, -1, 0, 0, 0, 1>;
/// Version-one, stateless equilibrium law; not a finite-time uptake curve.
pub const AQUEOUS_LANGMUIR_LAW: &str = "aqueous-langmuir-equilibrium";
/// Positive saturation loading, mol/kg dry adsorbent.
pub const CAPACITY_PARAMETER: &str = "saturation-loading";
/// Positive solution concentration giving half saturation, mol/m³.
pub const HALF_SATURATION_PARAMETER: &str = "half-saturation-concentration";
/// Actual solved solution concentration; callers cannot override it.
pub const CONCENTRATION_AXIS: &str = "aqueous-equilibrium-concentration";
/// Actual absolute contact temperature; no temperature correction is inferred.
pub const TEMPERATURE_AXIS: &str = "temperature";

/// Refusal leaves the supplied batch and original source unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum AqueousAdsorptionError {
    /// Original card or typed-coordinate refusal.
    Material(MatDbError),
    /// Missing, inconsistent or nonphysical model/input declaration.
    Invalid(&'static str),
    /// A solved equilibrium or other declared condition lacks source support.
    OutsideValidity {
        /// Source coordinate in its original convention.
        axis: String,
        /// Requested value in coherent SI.
        value: f64,
        /// Inclusive source lower bound.
        minimum: f64,
        /// Inclusive source upper bound.
        maximum: f64,
    },
    /// Arithmetic cannot retain a required positive amount or close its balance.
    Arithmetic(&'static str),
}

impl core::fmt::Display for AqueousAdsorptionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "aqueous adsorption: {self:?}")
    }
}
impl std::error::Error for AqueousAdsorptionError {}
impl From<MatDbError> for AqueousAdsorptionError {
    fn from(value: MatDbError) -> Self {
        Self::Material(value)
    }
}
type Result<T> = core::result::Result<T, AqueousAdsorptionError>;

/// Explicit analytical-metal stocks; adsorbed amount is not inferred from a grade.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdsorptionInventory {
    /// Metal in the fixed-volume solution [mol].
    pub dissolved: Amount,
    /// Metal already on the supplied dry adsorbent [mol].
    pub adsorbed: Amount,
}

/// Immutable original equilibrium card and its explicit analytical-element basis.
#[derive(Clone, Debug, PartialEq)]
pub struct AqueousLangmuir {
    element: ElementId,
    source: ConstitutiveModelCard,
    capacity: AdsorbentLoading,
    half_saturation: AmountConcentration,
}

/// Source-admitted isotherm response at instantaneous solution conditions.
/// This is a local loading target, not a finite-contact inventory solution.
#[derive(Clone, Debug, PartialEq)]
pub struct LangmuirLoading {
    model: AqueousLangmuir,
    /// Actual aqueous analytical-element concentration [mol/m³].
    pub concentration: AmountConcentration,
    /// Actual absolute temperature with the original semantic source coordinate.
    pub temperature: Temperature,
    /// Equilibrium loading per unit dry adsorbent mass [mol/kg].
    pub loading: AdsorbentLoading,
    /// Equilibrium vacant capacity [mol/kg], evaluated independently of loading.
    pub vacant_loading: AdsorbentLoading,
    /// Complete actual source point, including explicitly supplied extra conditions.
    pub conditions: QueryPoint,
}
impl LangmuirLoading {
    /// Original isotherm, analytical-element identity and immutable source card.
    pub fn model(&self) -> &AqueousLangmuir {
        &self.model
    }
}

/// Equilibrium of one complete, finite contact; binary64 balances are not certificates.
#[derive(Clone, Debug)]
pub struct AdsorptionContact {
    model: AqueousLangmuir,
    /// Original and solved analytical-metal inventories.
    pub initial: AdsorptionInventory,
    /// Final solution and adsorbent stocks, evaluated without clipping.
    pub final_inventory: AdsorptionInventory,
    /// Actual constant solution volume.
    pub solution_volume: Volume,
    /// Actual dry adsorbent mass.
    pub adsorbent_mass: Mass,
    /// Solved concentration and all admitted original source conditions.
    pub conditions: QueryPoint,
    /// Final dissolved concentration [mol/m³].
    pub equilibrium_concentration: AmountConcentration,
    /// Final loading [mol/kg dry adsorbent]; absent when no adsorbent exists.
    pub loading: Option<AdsorbentLoading>,
    /// Available sites [mol metal], evaluated independently of occupied sites.
    pub vacant_capacity: Amount,
    /// Signed transfer to the adsorbent [mol]; negative denotes desorption.
    pub transferred_to_adsorbent: Amount,
    /// Independently recomputed final-minus-initial metal balance [mol].
    pub conservation_residual: Amount,
}

impl AdsorptionContact {
    /// Complete source owner; no new fit or stronger source authority is minted.
    pub fn model(&self) -> &AqueousLangmuir {
        &self.model
    }
}

impl AqueousLangmuir {
    /// Admit exactly two dimensioned parameters and finite typed source support.
    /// The card must declare a stateless version-one equilibrium law, original
    /// source hashes, concentration support, and absolute-temperature support.
    pub fn from_card(element: ElementId, source: ConstitutiveModelCard) -> Result<Self> {
        source.validate()?;
        if source.law.0 != AQUEOUS_LANGMUIR_LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || source.sources.is_empty()
            || source.parameters.len() != 2
            || source.validity.bound("elapsed-time").is_some()
        {
            return Err(AqueousAdsorptionError::Invalid(
                "require a sourced stateless version-one aqueous equilibrium card",
            ));
        }
        let capacity = parameter(&source, CAPACITY_PARAMETER, AdsorbentLoading::DIMS)?;
        let half_saturation = parameter(
            &source,
            HALF_SATURATION_PARAMETER,
            AmountConcentration::DIMS,
        )?;
        for (axis, &(lo, hi)) in source.validity.bounds() {
            if !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
                || !source.validity.axis_quantities().contains_key(axis)
            {
                return Err(AqueousAdsorptionError::Invalid(
                    "every source condition needs finite ordered typed bounds",
                ));
            }
        }
        for (axis, quantity, positive_lower) in [
            (CONCENTRATION_AXIS, concentration_spec(), false),
            (TEMPERATURE_AXIS, temperature_spec(), true),
        ] {
            let (lo, hi) = source
                .validity
                .bound(axis)
                .ok_or(AqueousAdsorptionError::Invalid(
                    "missing required source domain",
                ))?;
            if source.validity.axis_quantities().get(axis) != Some(&quantity)
                || lo < 0.0
                || (positive_lower && lo == 0.0)
                || hi <= 0.0
            {
                return Err(AqueousAdsorptionError::Invalid(
                    "invalid concentration or absolute-temperature domain",
                ));
            }
        }
        Ok(Self {
            element,
            source,
            capacity: AdsorbentLoading::new(capacity),
            half_saturation: AmountConcentration::new(half_saturation),
        })
    }

    /// Analytical element conserved between solution and adsorbent.
    pub fn element(&self) -> &ElementId {
        &self.element
    }
    /// Original law, parameters, typed support, source hashes and provenance.
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Original saturation loading, mol/kg dry adsorbent.
    pub fn capacity(&self) -> AdsorbentLoading {
        self.capacity
    }
    /// Original half-saturation solution concentration, mol/m³.
    pub fn half_saturation(&self) -> AmountConcentration {
        self.half_saturation
    }

    /// Evaluate the original isotherm at the current solution concentration.
    ///
    /// This owner supplies actual concentration and semantic absolute temperature;
    /// caller overrides refuse. Every additional original source condition must
    /// be supplied explicitly. Occupied and vacant loadings use the original
    /// stable fractions, retaining a representable vacant tail when the occupied
    /// loading rounds to capacity. No liquid inventory, finite-contact balance,
    /// kinetic coefficient or changed source authority is inferred.
    pub fn loading_at(
        &self,
        concentration: AmountConcentration,
        temperature: Temperature,
        conditions: &QueryPoint,
    ) -> Result<LangmuirLoading> {
        nonnegative(
            concentration.value(),
            "instantaneous concentration must be nonnegative",
        )?;
        positive(temperature.value(), "absolute temperature must be positive")?;
        if [CONCENTRATION_AXIS, TEMPERATURE_AXIS]
            .iter()
            .any(|axis| conditions.axes().contains_key(*axis))
        {
            return Err(AqueousAdsorptionError::Invalid(
                "caller cannot override actual concentration or temperature",
            ));
        }
        let query = conditions
            .clone()
            .with_quantity(
                CONCENTRATION_AXIS,
                concentration_spec(),
                concentration.value(),
            )?
            .with_quantity(TEMPERATURE_AXIS, temperature_spec(), temperature.value())?;
        self.admit(&query)?;
        let (occupied, vacant) = fractions(concentration.value(), self.half_saturation.value())?;
        let loading = multiply(
            self.capacity.value(),
            occupied,
            "instantaneous equilibrium loading",
        )?;
        let vacant_loading = multiply(
            self.capacity.value(),
            vacant,
            "instantaneous equilibrium vacant loading",
        )?;
        Ok(LangmuirLoading {
            model: self.clone(),
            concentration,
            temperature,
            loading: AdsorbentLoading::new(loading),
            vacant_loading: AdsorbentLoading::new(vacant_loading),
            conditions: query,
        })
    }

    /// Solve one ideal equilibrium contact at fixed volume and dry adsorbent mass.
    ///
    /// Supply every additional source condition (for example pH or competing
    /// ion concentration) explicitly. Those coordinates stay fixed. The owner
    /// supplies temperature and solved concentration; caller overrides refuse.
    /// Nonzero initial loading is allowed only up to the declared site capacity.
    /// No contact time, removal percentage or rate is inferred.
    pub fn equilibrate(
        &self,
        initial: AdsorptionInventory,
        solution_volume: Volume,
        adsorbent_mass: Mass,
        temperature: Temperature,
        conditions: &QueryPoint,
    ) -> Result<AdsorptionContact> {
        nonnegative(initial.dissolved.value(), "invalid dissolved inventory")?;
        nonnegative(initial.adsorbed.value(), "invalid adsorbed inventory")?;
        positive(solution_volume.value(), "solution volume must be positive")?;
        nonnegative(
            adsorbent_mass.value(),
            "dry adsorbent mass must be nonnegative",
        )?;
        positive(temperature.value(), "absolute temperature must be positive")?;
        if conditions.axes().contains_key(CONCENTRATION_AXIS)
            || conditions.axes().contains_key(TEMPERATURE_AXIS)
        {
            return Err(AqueousAdsorptionError::Invalid(
                "caller cannot override actual temperature or solved concentration",
            ));
        }
        let capacity = multiply(
            self.capacity.value(),
            adsorbent_mass.value(),
            "total site capacity",
        )?;
        if initial.adsorbed.value() > capacity {
            return Err(AqueousAdsorptionError::Invalid(
                "initial adsorbed stock exceeds the supplied adsorbent capacity",
            ));
        }
        let total = initial.dissolved.value() + initial.adsorbed.value();
        arithmetic(total, false, "total metal inventory")?;
        let half_amount = multiply(
            self.half_saturation.value(),
            solution_volume.value(),
            "half-saturation amount",
        )?;
        // x is final dissolved amount. Its quadratic is
        // x² + (B + C - N)x - NB = 0, B=V*c_half, C=m*q_max.
        // Scale coefficients and use the opposite-root identity when the
        // positive root would otherwise subtract two nearly equal numbers.
        let dissolved = if total == 0.0 || capacity == 0.0 {
            total
        } else {
            let scale = total.max(capacity).max(half_amount);
            let n = divide(total, scale, "scaled inventory")?;
            let c = divide(capacity, scale, "scaled site capacity")?;
            let b = divide(half_amount, scale, "scaled half saturation")?;
            let linear = sum3(b, c, -n);
            let cross = multiply(4.0 * n, b, "quadratic discriminant")?;
            let root = det::sqrt(linear.mul_add(linear, cross));
            let scaled = if linear >= 0.0 {
                let ratio = divide(2.0 * b, root + linear, "dissolved root ratio")?;
                multiply(n, ratio, "dissolved root")?
            } else {
                0.5 * (root - linear)
            };
            multiply(scaled, scale, "dissolved equilibrium inventory")?
        };
        arithmetic(dissolved, total > 0.0, "positive dissolved inventory")?;
        if dissolved > total {
            return Err(AqueousAdsorptionError::Arithmetic(
                "dissolved equilibrium exceeds total inventory",
            ));
        }
        let concentration = divide(
            dissolved,
            solution_volume.value(),
            "equilibrium concentration",
        )?;
        let (occupied, vacant) = fractions(concentration, self.half_saturation.value())?;
        let adsorbed = multiply(capacity, occupied, "adsorbed equilibrium inventory")?;
        let vacant_capacity = multiply(capacity, vacant, "vacant site inventory")?;
        let residual = sum3(dissolved, adsorbed, -total);
        let sites_residual = sum3(adsorbed, vacant_capacity, -capacity);
        if !residual.is_finite()
            || !sites_residual.is_finite()
            || (total > 0.0 && residual.abs() / total > 128.0 * f64::EPSILON)
            || (capacity > 0.0 && sites_residual.abs() / capacity > 128.0 * f64::EPSILON)
        {
            return Err(AqueousAdsorptionError::Arithmetic(
                "equilibrium metal or site balance does not close",
            ));
        }
        let query = conditions
            .clone()
            .with_quantity(TEMPERATURE_AXIS, temperature_spec(), temperature.value())?
            .with_quantity(CONCENTRATION_AXIS, concentration_spec(), concentration)?;
        self.admit(&query)?;
        let loading = if adsorbent_mass.value() == 0.0 {
            None
        } else {
            Some(AdsorbentLoading::new(divide(
                adsorbed,
                adsorbent_mass.value(),
                "equilibrium loading",
            )?))
        };
        Ok(AdsorptionContact {
            model: self.clone(),
            initial,
            final_inventory: AdsorptionInventory {
                dissolved: Amount::new(dissolved),
                adsorbed: Amount::new(adsorbed),
            },
            solution_volume,
            adsorbent_mass,
            conditions: query,
            equilibrium_concentration: AmountConcentration::new(concentration),
            loading,
            vacant_capacity: Amount::new(vacant_capacity),
            transferred_to_adsorbent: Amount::new(adsorbed - initial.adsorbed.value()),
            conservation_residual: Amount::new(residual),
        })
    }

    fn admit(&self, query: &QueryPoint) -> Result<()> {
        for (axis, &(minimum, maximum)) in self.source.validity.bounds() {
            let value = *query
                .axes()
                .get(axis)
                .ok_or(AqueousAdsorptionError::Invalid(
                    "missing original source condition",
                ))?;
            if query.axis_quantities().get(axis) != self.source.validity.axis_quantities().get(axis)
            {
                return Err(AqueousAdsorptionError::Invalid(
                    "source condition has wrong quantity",
                ));
            }
            if value < minimum || value > maximum {
                return Err(AqueousAdsorptionError::OutsideValidity {
                    axis: axis.clone(),
                    value,
                    minimum,
                    maximum,
                });
            }
        }
        Ok(())
    }
}

fn parameter(source: &ConstitutiveModelCard, name: &str, dims: Dims) -> Result<f64> {
    let value = source
        .parameters
        .get(name)
        .ok_or(AqueousAdsorptionError::Invalid(
            "missing Langmuir parameter",
        ))?;
    if value.dims != dims {
        return Err(AqueousAdsorptionError::Invalid(
            "Langmuir parameter has wrong dimensions",
        ));
    }
    positive(value.value, "Langmuir parameters must be positive")?;
    Ok(value.value)
}

fn concentration_spec() -> QuantitySpec {
    QuantitySpec::dimensional(AmountConcentration::DIMS)
}
fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn nonnegative(value: f64, reason: &'static str) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        Err(AqueousAdsorptionError::Invalid(reason))
    } else {
        Ok(())
    }
}
fn positive(value: f64, reason: &'static str) -> Result<()> {
    if !value.is_finite() || value <= 0.0 {
        Err(AqueousAdsorptionError::Invalid(reason))
    } else {
        Ok(())
    }
}
fn arithmetic(value: f64, positive: bool, reason: &'static str) -> Result<f64> {
    if !value.is_finite() || value < 0.0 || (positive && value == 0.0) {
        Err(AqueousAdsorptionError::Arithmetic(reason))
    } else {
        Ok(value)
    }
}
fn multiply(a: f64, b: f64, reason: &'static str) -> Result<f64> {
    arithmetic(a * b, a > 0.0 && b > 0.0, reason)
}
fn divide(a: f64, b: f64, reason: &'static str) -> Result<f64> {
    arithmetic(a / b, a > 0.0, reason)
}
fn fractions(c: f64, half: f64) -> Result<(f64, f64)> {
    if c <= half {
        let ratio = divide(c, half, "occupied fraction ratio")?;
        Ok((ratio / (1.0 + ratio), 1.0 / (1.0 + ratio)))
    } else {
        let ratio = divide(half, c, "vacant fraction ratio")?;
        Ok((1.0 / (1.0 + ratio), ratio / (1.0 + ratio)))
    }
}
fn sum3(a: f64, b: f64, c: f64) -> f64 {
    let mut total: f64 = 0.0;
    let mut correction: f64 = 0.0;
    for value in [a, b, c] {
        let next = total + value;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    total + correction
}

#[cfg(test)]
#[path = "aqueous_adsorption/tests.rs"]
mod tests;
