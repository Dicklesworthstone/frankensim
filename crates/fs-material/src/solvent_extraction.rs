//! Finite-inventory equilibrium contact between aqueous and organic phases.
//!
//! A source-declared table gives organic total analytical concentration `y`
//! versus aqueous concentration `x`, both in mol/m³ of the named element.
//! This axis convention is unchanged for stripping. The IUPAC concentration
//! distribution convention is organic/aqueous:
//! <https://doi.org/10.1351/goldbook.D01817>. A primary extraction/stripping
//! study explicitly distinguishes equal-volume formulas from balances using
//! both phase volumes: <https://doi.org/10.3390/met12061053>, equations 4–7.
//!
//! For fixed positive phase volumes, an ideal contact solves
//! `Vaq*x + Vorg*y(x) = naq_initial + norg_initial`. Strictly increasing x
//! and nondecreasing y make this balance strictly increasing, including flat
//! organic-loading segments. The intersection is solved algebraically on the
//! admitted piecewise-linear segment; there is no iterative stopping tolerance
//! and no extrapolation beyond the original table. Double-double intermediates
//! retain small balance differences; binary64 outputs are not a certificate.
//!
//! This is one analytical-element balance, not a species/charge balance or a
//! prediction of pH, extractant consumption, speciation, heat, rate, entrainment,
//! stage efficiency or phase-volume changes. Different acid/extractant systems
//! need their own source cards. Equilibrium concentration support applies to
//! the solved endpoint, not to initially nonequilibrated feed concentrations.

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::dd::Dd;
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, QuantitySpec, Temperature, Volume,
    chemistry::ElementId,
};

/// Version-1, stateless, linear interpolation of an original equilibrium table.
pub const ISOTHERM_LAW: &str = "tabulated-solvent-extraction-isotherm";
/// Parameters `aqueous-knot-0`, `aqueous-knot-1`, ... in mol/m³.
pub const AQUEOUS_KNOT_PREFIX: &str = "aqueous-knot-";
/// Matching parameters `organic-knot-0`, `organic-knot-1`, ... in mol/m³.
pub const ORGANIC_KNOT_PREFIX: &str = "organic-knot-";
pub const AQUEOUS_CONCENTRATION_AXIS: &str = "aqueous-equilibrium-concentration";
pub const ORGANIC_CONCENTRATION_AXIS: &str = "organic-equilibrium-concentration";
pub const VOLUME_RATIO_AXIS: &str = "organic-aqueous-volume-ratio";
pub const TEMPERATURE_AXIS: &str = "temperature";
/// Bounded admission and segment search, without a nonlinear plant solver.
pub const MAX_KNOTS: usize = 256;

const OWNED_AXES: [&str; 4] = [
    AQUEOUS_CONCENTRATION_AXIS,
    ORGANIC_CONCENTRATION_AXIS,
    VOLUME_RATIO_AXIS,
    TEMPERATURE_AXIS,
];
type Result<T> = core::result::Result<T, SolventExtractionError>;

/// Invalid declaration, unsupported equilibrium or unrepresentable arithmetic.
#[derive(Clone, Debug, PartialEq)]
pub enum SolventExtractionError {
    Source(MatDbError),
    Invalid(&'static str),
    /// The closed-system balance has no intersection on the supplied table.
    OutsideIsotherm,
    /// Actual equilibrium or other conditions leave the original source domain.
    OutsideValidity,
    Arithmetic(&'static str),
}
impl core::fmt::Display for SolventExtractionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "solvent extraction: {self:?}")
    }
}
impl std::error::Error for SolventExtractionError {}
impl From<MatDbError> for SolventExtractionError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}

/// A phase's finite total analytical inventory; chemical species are unspecified.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseInventory {
    solute: ElementId,
    volume: Volume,
    amount: Amount,
}
impl PhaseInventory {
    pub fn new(solute: ElementId, volume: Volume, amount: Amount) -> Result<Self> {
        if !volume.value().is_finite() || volume.value() <= 0.0 {
            return Err(SolventExtractionError::Invalid(
                "phase volume must be finite and positive",
            ));
        }
        if !amount.value().is_finite() || amount.value() < 0.0 {
            return Err(SolventExtractionError::Invalid(
                "solute amount must be finite and nonnegative",
            ));
        }
        represented(amount.value() / volume.value(), amount.value() > 0.0)?;
        Ok(Self {
            solute,
            volume,
            amount,
        })
    }
    pub fn solute(&self) -> &ElementId {
        &self.solute
    }
    pub const fn volume(&self) -> Volume {
        self.volume
    }
    pub const fn amount(&self) -> Amount {
        self.amount
    }
    pub fn concentration(&self) -> AmountConcentration {
        AmountConcentration::new(self.amount.value() / self.volume.value())
    }
}

/// One original table point; organic loading is always the ordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumPoint {
    pub aqueous: AmountConcentration,
    pub organic: AmountConcentration,
}

/// An explicitly identified analytical element and its complete original card.
#[derive(Clone, Debug, PartialEq)]
pub struct SolventExtractionIsotherm {
    solute: ElementId,
    source: ConstitutiveModelCard,
    points: Vec<EquilibriumPoint>,
}

/// Both final phases, their original source, and the finite material account.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumStage {
    model: SolventExtractionIsotherm,
    pub initial_aqueous: PhaseInventory,
    pub initial_organic: PhaseInventory,
    pub aqueous: PhaseInventory,
    pub organic: PhaseInventory,
    /// Positive for extraction, negative for stripping, in mol of the element.
    pub transfer_to_organic: Amount,
    /// Final minus initial total, calculated independently, without a balance plug.
    pub conservation_residual: Amount,
    /// Index of the lower knot; an internal knot chooses the segment to its left.
    pub segment_index: usize,
    /// Actual final concentrations, volume ratio, temperature and source conditions.
    pub conditions: QueryPoint,
}
impl EquilibriumStage {
    pub fn model(&self) -> &SolventExtractionIsotherm {
        &self.model
    }
}

impl SolventExtractionIsotherm {
    /// Read all knots from the original card, keeping their provenance attached.
    ///
    /// Exactly 2..=MAX_KNOTS zero-indexed pairs of concentration parameters are
    /// required. No sorting, smoothing, inferred origin or extra parameters are
    /// allowed. All source axes need finite typed bounds; the four owned axes
    /// must use the canonical quantities declared by this module. The supplied
    /// element is an explicit caller binding, not inferred from a material name.
    pub fn from_card(solute: ElementId, source: ConstitutiveModelCard) -> Result<Self> {
        source.validate()?;
        let count = source.parameters.len() / 2;
        if source.law.0 != ISOTHERM_LAW
            || source.law_version != 1
            || source.state_schema_version != 0
            || source.initial_state != InitialStatePolicy::ZeroInternalState
            || !(2..=MAX_KNOTS).contains(&count)
            || source.parameters.len() != count * 2
            || source.sources.is_empty()
            || source.provenance.source.trim().is_empty()
            || source.provenance.license.trim().is_empty()
        {
            return Err(SolventExtractionError::Invalid(
                "requires version 1, zero state, paired knots and original sources",
            ));
        }
        for (axis, &(lo, hi)) in source.validity.bounds() {
            if !lo.is_finite()
                || !hi.is_finite()
                || lo > hi
                || !source.validity.axis_quantities().contains_key(axis)
            {
                return Err(SolventExtractionError::Invalid(
                    "every source axis needs finite typed bounds",
                ));
            }
        }
        for (axis, dims) in [
            (AQUEOUS_CONCENTRATION_AXIS, AmountConcentration::DIMS),
            (ORGANIC_CONCENTRATION_AXIS, AmountConcentration::DIMS),
            (VOLUME_RATIO_AXIS, Dimensionless::DIMS),
            (TEMPERATURE_AXIS, Temperature::DIMS),
        ] {
            let (lo, _) = source
                .validity
                .bound(axis)
                .ok_or(SolventExtractionError::Invalid("missing owned source axis"))?;
            if source.validity.axis_quantities().get(axis) != Some(&QuantitySpec::dimensional(dims))
                || lo < 0.0
                || ((axis == VOLUME_RATIO_AXIS || axis == TEMPERATURE_AXIS) && lo == 0.0)
            {
                return Err(SolventExtractionError::Invalid(
                    "owned axes need canonical quantities and physical bounds",
                ));
            }
        }
        let parameter = |prefix: &str, index: usize| -> Result<f64> {
            let p = source.parameters.get(&format!("{prefix}{index}")).ok_or(
                SolventExtractionError::Invalid("missing indexed concentration knot"),
            )?;
            if p.dims != AmountConcentration::DIMS || p.value < 0.0 {
                return Err(SolventExtractionError::Invalid(
                    "knots must be nonnegative mol/m3",
                ));
            }
            Ok(p.value)
        };
        let mut points: Vec<EquilibriumPoint> = Vec::with_capacity(count);
        for index in 0..count {
            let x = parameter(AQUEOUS_KNOT_PREFIX, index)?;
            let y = parameter(ORGANIC_KNOT_PREFIX, index)?;
            if points
                .last()
                .is_some_and(|p| x <= p.aqueous.value() || y < p.organic.value())
            {
                return Err(SolventExtractionError::Invalid(
                    "aqueous knots must increase strictly; organic knots must not decrease",
                ));
            }
            for (axis, value) in [
                (AQUEOUS_CONCENTRATION_AXIS, x),
                (ORGANIC_CONCENTRATION_AXIS, y),
            ] {
                let (lo, hi) = source.validity.bound(axis).expect("validated axis");
                if value < lo || value > hi {
                    return Err(SolventExtractionError::Invalid(
                        "table exceeds original concentration support",
                    ));
                }
            }
            points.push(EquilibriumPoint {
                aqueous: AmountConcentration::new(x),
                organic: AmountConcentration::new(y),
            });
        }
        Ok(Self {
            solute,
            source,
            points,
        })
    }
    pub fn solute(&self) -> &ElementId {
        &self.solute
    }
    pub fn source(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    pub fn points(&self) -> &[EquilibriumPoint] {
        &self.points
    }

    /// Equilibrate a closed contact without mutating either initial inventory.
    ///
    /// Extraction and stripping use this same balance and axis convention;
    /// the source chemistry and starting inventories determine the direction.
    /// All other declared coordinates must be supplied in `conditions`. The
    /// caller cannot replace any owned coordinate, even for zero transfer.
    pub fn equilibrate(
        &self,
        aqueous: &PhaseInventory,
        organic: &PhaseInventory,
        temperature: Temperature,
        conditions: &QueryPoint,
    ) -> Result<EquilibriumStage> {
        if aqueous.solute != self.solute || organic.solute != self.solute {
            return Err(SolventExtractionError::Invalid(
                "both inventories must name the isotherm's analytical element",
            ));
        }
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(SolventExtractionError::Invalid(
                "temperature must be finite positive kelvin",
            ));
        }
        if OWNED_AXES
            .iter()
            .any(|axis| conditions.axes().contains_key(*axis))
        {
            return Err(SolventExtractionError::Invalid(
                "caller conditions override owned equilibrium axes",
            ));
        }
        let va = aqueous.volume.value();
        let vo = organic.volume.value();
        let ratio = represented(vo / va, true)?;
        let total = sum(
            Dd::from_f64(aqueous.amount.value()),
            Dd::from_f64(organic.amount.value()),
        )?;
        let totals = self
            .points
            .iter()
            .map(|p| {
                sum(
                    product(va, p.aqueous.value())?,
                    product(vo, p.organic.value())?,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        if total.lt(totals[0]) || totals[totals.len() - 1].lt(total) {
            return Err(SolventExtractionError::OutsideIsotherm);
        }
        let index = totals
            .windows(2)
            .position(|t| !t[1].lt(total))
            .ok_or(SolventExtractionError::OutsideIsotherm)?;
        let lo = self.points[index];
        let hi = self.points[index + 1];
        let width = totals[index + 1] - totals[index];
        if !Dd::ZERO.lt(width) {
            return Err(SolventExtractionError::Arithmetic(
                "unresolved segment balance width",
            ));
        }
        let fraction = (total - totals[index]) / width;
        let interpolate = |a: f64, b: f64| -> Result<f64> {
            let value = Dd::from_f64(a) + fraction * (Dd::from_f64(b) - Dd::from_f64(a));
            let value = represented(
                value.to_f64(),
                a > 0.0 || (b > 0.0 && Dd::ZERO.lt(fraction)),
            )?;
            if value < a || value > b {
                return Err(SolventExtractionError::Arithmetic(
                    "interpolated concentration escaped its segment",
                ));
            }
            Ok(value)
        };
        let x = interpolate(lo.aqueous.value(), hi.aqueous.value())?;
        let y = interpolate(lo.organic.value(), hi.organic.value())?;
        let final_aqueous = PhaseInventory::new(
            self.solute.clone(),
            aqueous.volume,
            Amount::new(product(va, x)?.to_f64()),
        )?;
        let final_organic = PhaseInventory::new(
            self.solute.clone(),
            organic.volume,
            Amount::new(product(vo, y)?.to_f64()),
        )?;
        let mut point = conditions.clone();
        for (axis, dims, value) in [
            (
                AQUEOUS_CONCENTRATION_AXIS,
                AmountConcentration::DIMS,
                final_aqueous.concentration().value(),
            ),
            (
                ORGANIC_CONCENTRATION_AXIS,
                AmountConcentration::DIMS,
                final_organic.concentration().value(),
            ),
            (VOLUME_RATIO_AXIS, Dimensionless::DIMS, ratio),
            (TEMPERATURE_AXIS, Temperature::DIMS, temperature.value()),
        ] {
            point = point.with_quantity(axis, QuantitySpec::dimensional(dims), value)?;
        }
        if !self
            .source
            .validity
            .contains_typed(point.axes(), point.axis_quantities())
        {
            return Err(SolventExtractionError::OutsideValidity);
        }
        let final_total = sum(
            Dd::from_f64(final_aqueous.amount.value()),
            Dd::from_f64(final_organic.amount.value()),
        )?;
        let residual = (final_total - total).to_f64();
        if !residual.is_finite()
            || (total.to_f64() > 0.0 && residual.abs() / total.to_f64() > 32.0 * f64::EPSILON)
        {
            return Err(SolventExtractionError::Arithmetic(
                "equilibrium material balance failed",
            ));
        }
        Ok(EquilibriumStage {
            model: self.clone(),
            initial_aqueous: aqueous.clone(),
            initial_organic: organic.clone(),
            transfer_to_organic: Amount::new(final_organic.amount.value() - organic.amount.value()),
            aqueous: final_aqueous,
            organic: final_organic,
            conservation_residual: Amount::new(residual),
            segment_index: index,
            conditions: point,
        })
    }
}

fn represented(value: f64, positive: bool) -> Result<f64> {
    if !value.is_finite() || value < 0.0 || (positive && value == 0.0) {
        return Err(SolventExtractionError::Arithmetic(
            "positive quantity overflowed or underflowed",
        ));
    }
    Ok(value)
}
fn product(a: f64, b: f64) -> Result<Dd> {
    represented(a * b, a > 0.0 && b > 0.0)?;
    Ok(Dd::from_f64(a) * Dd::from_f64(b))
}
fn sum(a: Dd, b: Dd) -> Result<Dd> {
    represented(a.to_f64() + b.to_f64(), false)?;
    Ok(a + b)
}

#[cfg(test)]
#[path = "solvent_extraction/tests.rs"]
mod tests;
