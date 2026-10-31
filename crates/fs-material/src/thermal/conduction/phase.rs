//! Compiled material conductivity on the original sensible/latent enthalpy path.
//!
//! A scalar source must describe the admitted heat-flow direction. No isotropy,
//! phase mixing rule, fluid motion or changing reference geometry is inferred.
//! Pure-phase temperature segments and isothermal liquid-mass-fraction segments
//! reuse the canonical one-coordinate source-line owner. Every source crossing
//! and caloric knot survives compilation. Simultaneously changing temperature
//! and fraction require a multidimensional path owner and are refused here.
//!
//! The native conduction observer and implicit chain/graph retain all thermal
//! history, storage, fluxes, nonlinear solves and energy acceptance. Evaluation
//! uses the original caloric interpolation plus a binary search of this immutable
//! conductivity path; it performs no source queries, allocations or hashing.

use crate::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError, SolidLiquidPhase};
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use crate::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimId, ClaimSelection, PropertyKey, PropertyLineSlice, PropertySupportError};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Combined source/caloric knot ceiling; each individual source line also has
/// the canonical `MAX_PROPERTY_LINE_KNOTS` bound. No adaptive resampling occurs.
pub const MAX_PHASE_CONDUCTIVITY_KNOTS: usize = 65_536;

/// Original source/caloric refusals remain inspectable during preparation and
/// direct evaluation. Native thermal traits lower these to `MaterialDomain`.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseConductivityError {
    Invalid(&'static str),
    Phase(PhaseStateError),
    Material(ThermalMaterialError),
    Source(PropertySupportError),
    Thermal(DepositionError),
    OutsideEnthalpyDomain {
        specific_enthalpy_j_kg: f64,
        bounds_j_kg: [f64; 2],
    },
    Unrepresentable(&'static str),
}
impl core::fmt::Display for PhaseConductivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "phase conductivity: {self:?}")
    }
}
impl std::error::Error for PhaseConductivityError {}
impl From<PhaseStateError> for PhaseConductivityError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<ThermalMaterialError> for PhaseConductivityError {
    fn from(error: ThermalMaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<PropertySupportError> for PhaseConductivityError {
    fn from(error: PropertySupportError) -> Self {
        Self::Source(error)
    }
}
impl From<DepositionError> for PhaseConductivityError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

#[derive(Clone, Copy, Debug)]
struct Interval {
    temperature_slope: f64,
    conductivity_slope: f64,
    source_line: usize,
}

/// Exact compiled scalar response, distinct from a fresh source-query receipt.
/// Source authority remains on the retained lines, not on this numerical sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseConductivitySample {
    specific_enthalpy: f64,
    temperature: f64,
    liquid_fraction: f64,
    conductivity: f64,
    temperature_slope: f64,
    conductivity_slope: f64,
    source_line: usize,
}
impl PhaseConductivitySample {
    pub const fn specific_enthalpy_j_kg(self) -> f64 {
        self.specific_enthalpy
    }
    pub const fn temperature_k(self) -> f64 {
        self.temperature
    }
    pub const fn liquid_mass_fraction(self) -> f64 {
        self.liquid_fraction
    }
    pub fn phase(self) -> SolidLiquidPhase {
        if self.liquid_fraction == 0.0 {
            SolidLiquidPhase::Solid
        } else if self.liquid_fraction == 1.0 {
            SolidLiquidPhase::Liquid
        } else {
            SolidLiquidPhase::SolidLiquid
        }
    }
    pub const fn conductivity_w_m_k(self) -> f64 {
        self.conductivity
    }
    pub const fn temperature_enthalpy_derivative_kg_k_j(self) -> f64 {
        self.temperature_slope
    }
    pub const fn conductivity_enthalpy_derivative_kg_w_m_k_j(self) -> f64 {
        self.conductivity_slope
    }
    /// Retained source line used by the selected one-sided interval.
    pub const fn source_line_index(self) -> usize {
        self.source_line
    }
}

/// One pinned source compiled on an explicitly bounded original caloric path.
///
/// Registration in both native conduction traits uses this SAME object; its
/// caloric binding is the original borrowed curve, not a cloned replacement.
/// The declared source context supplies fixed pressure/composition/etc. Only
/// the named temperature and optional liquid mass fraction coordinates evolve.
#[derive(Clone, Debug)]
pub struct PreparedPhaseConductivity<'a> {
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    knots: Vec<[f64; 2]>,
    intervals: Vec<Interval>,
    source_lines: Vec<PropertyLineSlice>,
    selected: ClaimId,
    identity: ContentHash,
}

impl<'a> PreparedPhaseConductivity<'a> {
    /// Compile continuous source support over the complete requested h interval.
    ///
    /// The initial claim selection is pinned across every caloric segment.
    /// An evolving phase fraction requires its declared existing source axis,
    /// including when the selected conductivity happens to be phase independent.
    /// Source holes, discrete-only data, unsupported phase paths, negative k,
    /// source switching and unrepresentable breakpoints/tangents refuse without
    /// publishing a partial law. Explicit constant-zero insulation is admitted;
    /// a zero k with nonzero selected tangent is outside the native implicit law.
    /// Source-line admission is synchronous and individually bounded; cancellation
    /// is checked around every query and at most every 256 retained source knots.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        caloric: &'a EquilibriumEnthalpyPhaseCurve,
        point: ThermalMaterialPoint<'_>,
        key: &PropertyKey,
        selection: ClaimSelection,
        enthalpy_range_j_kg: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PhaseConductivityError> {
        checkpoint()?;
        let invalid = PhaseConductivityError::Invalid;
        if key.quantity() != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && key.quantity()
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                ))
        {
            return Err(invalid(
                "a scalar static thermal conductivity key is required",
            ));
        }
        let [lo, hi] = enthalpy_range_j_kg;
        if !lo.is_finite() || !hi.is_finite() || lo >= hi || !(hi - lo).is_finite() {
            return Err(invalid("a finite increasing enthalpy interval is required"));
        }
        // These calls enforce the original card identity and typed coordinates.
        point.at(caloric.state_at_specific_enthalpy(lo)?)?;
        point.at(caloric.state_at_specific_enthalpy(hi)?)?;
        for receipt in caloric.source_receipts() {
            checkpoint()?;
            for (axis, value) in &receipt.query_point {
                if axis != point.temperature_axis
                    && Some(axis.as_str()) != point.liquid_fraction_axis
                    && (point.context.axes().get(axis) != Some(value)
                        || point.context.axis_quantities().get(axis)
                            != receipt.axis_quantities.get(axis))
                {
                    return Err(invalid(
                        "conductivity and caloric sources require the same fixed external conditions",
                    ));
                }
            }
        }
        let original = caloric.knots();
        let first = original
            .partition_point(|k| k.specific_enthalpy_j_kg <= lo)
            .saturating_sub(1)
            .min(original.len() - 2);
        let last = original
            .partition_point(|k| k.specific_enthalpy_j_kg < hi)
            .min(original.len() - 1);
        if last - first >= MAX_PHASE_CONDUCTIVITY_KNOTS {
            return Err(DepositionError::Budget("phase conductivity caloric knots").into());
        }
        let mut knots = Vec::new();
        let mut intervals = Vec::new();
        let mut source_lines = Vec::new();
        let mut selected = None;
        let mut hash = DomainHasher::new("org.frankensim.fs-material.phase-conductivity.v1");
        hash.update(caloric.identity().as_bytes());
        for h in enthalpy_range_j_kg {
            hash.update(&h.to_bits().to_le_bytes());
        }
        for pair in original[first..=last].windows(2) {
            checkpoint()?;
            let a = lo.max(pair[0].specific_enthalpy_j_kg);
            let b = hi.min(pair[1].specific_enthalpy_j_kg);
            let changes_t = pair[0].temperature_k != pair[1].temperature_k;
            let changes_fraction = pair[0].liquid_mass_fraction != pair[1].liquid_mass_fraction;
            if changes_t
                && (changes_fraction || ![0.0, 1.0].contains(&pair[0].liquid_mass_fraction))
            {
                return Err(invalid(
                    "temperature may vary only on a pure-phase segment; varying T and fraction needs a multidimensional source path",
                ));
            }
            if changes_fraction && point.liquid_fraction_axis.is_none() {
                return Err(invalid(
                    "latent transport requires an explicitly evolving liquid mass fraction axis",
                ));
            }
            let lower = point.at(caloric.state_at_specific_enthalpy(a)?)?;
            let upper = point.at(caloric.state_at_specific_enthalpy(b)?)?;
            let axis = if changes_fraction {
                point.liquid_fraction_axis.unwrap_or(point.temperature_axis)
            } else {
                point.temperature_axis
            };
            let line = point.material.claims().query_line_typed(
                key,
                axis,
                &lower,
                &upper,
                selected.map(ClaimSelection::Pinned).unwrap_or(selection),
            )?;
            checkpoint()?;
            let answer = &line
                .knots
                .first()
                .ok_or(invalid("empty conductivity source line"))?
                .1;
            if selected.is_none() {
                selected = Some(answer.receipt.selected);
                let requirement =
                    ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative)
                        .map_err(|error| ThermalMaterialError::Material(error))?;
                let bundle = resolve_material_state_point(
                    point.material,
                    &lower,
                    &[requirement],
                    MaterialPropertySelection::PinnedByProperty(vec![(
                        key.name().to_owned(),
                        answer.receipt.selected,
                    )]),
                )
                .map_err(|error| ThermalMaterialError::Material(error))?;
                point.admit_static(
                    bundle
                        .properties()
                        .first()
                        .ok_or(invalid("empty conductivity bundle"))?,
                )?;
            }
            let temperature_slope = caloric.temperature_derivative_at_specific_enthalpy(a)?;
            let x0 = lower.axes()[axis];
            let x1 = upper.axes()[axis];
            let line_index = source_lines.len();
            hash.update(&(line_index as u64).to_le_bytes());
            hash.update(&(line.knots.len() as u64).to_le_bytes());
            let mut previous = None;
            for (index, (x, answer)) in line.knots.iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint()?;
                }
                let k = answer.evidence.value.value;
                if !k.is_finite() || k < 0.0 {
                    return Err(invalid(
                        "source conductivity must be finite and nonnegative throughout the path",
                    ));
                }
                hash.update(answer.receipt.content_hash().as_bytes());
                let h = if index == 0 {
                    a
                } else if index + 1 == line.knots.len() {
                    b
                } else {
                    (b - a).mul_add((*x - x0) / (x1 - x0), a)
                };
                Self::append(
                    &mut knots,
                    &mut intervals,
                    &mut previous,
                    [h, k],
                    temperature_slope,
                    line_index,
                )?;
            }
            // A phase-independent source may be constant on the entire latent
            // interval. Retain both caloric endpoints even for one source point.
            if line.knots.len() == 1 {
                Self::append(
                    &mut knots,
                    &mut intervals,
                    &mut previous,
                    [b, line.knots[0].1.evidence.value.value],
                    temperature_slope,
                    line_index,
                )?;
            }
            source_lines.push(line);
        }
        hash.update(&(knots.len() as u64).to_le_bytes());
        for (index, knot) in knots.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let interval = &intervals[index.min(intervals.len() - 1)];
            if knot[1] == 0.0 && interval.conductivity_slope != 0.0 {
                return Err(invalid(
                    "zero conductivity requires a zero one-sided tangent for native implicit transport",
                ));
            }
            for value in knot {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        checkpoint()?;
        Ok(Self {
            caloric,
            knots,
            intervals,
            source_lines,
            selected: selected.ok_or(invalid("empty caloric conductivity path"))?,
            identity: hash.finalize(),
        })
    }

    fn append(
        knots: &mut Vec<[f64; 2]>,
        intervals: &mut Vec<Interval>,
        previous: &mut Option<[f64; 2]>,
        knot: [f64; 2],
        temperature_slope: f64,
        source_line: usize,
    ) -> Result<(), PhaseConductivityError> {
        if let Some(before) = *previous {
            let width = knot[0] - before[0];
            let change = knot[1] - before[1];
            let slope = change / width;
            if !width.is_finite()
                || width <= 0.0
                || !slope.is_finite()
                || (change != 0.0 && slope == 0.0)
            {
                return Err(PhaseConductivityError::Unrepresentable(
                    "source knots or conductivity tangent cannot be represented on the caloric interval",
                ));
            }
            intervals.push(Interval {
                temperature_slope,
                conductivity_slope: slope,
                source_line,
            });
        } else if let Some(last) = knots.last() {
            if last != &knot {
                return Err(PhaseConductivityError::Invalid(
                    "conductivity source must be continuous at a shared caloric knot",
                ));
            }
            *previous = Some(knot);
            return Ok(());
        }
        if knots.len() >= MAX_PHASE_CONDUCTIVITY_KNOTS {
            return Err(DepositionError::Budget("combined phase conductivity knots").into());
        }
        knots.push(knot);
        *previous = Some(knot);
        Ok(())
    }

    pub const fn caloric_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.caloric
    }
    pub const fn selected_claim(&self) -> ClaimId {
        self.selected
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Inclusive bounds in the original caloric enthalpy datum.
    pub fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [self.knots[0][0], self.knots[self.knots.len() - 1][0]]
    }
    /// Combined original caloric/source breakpoints, as [h in J/kg, k in W/(m K)].
    pub fn knots(&self) -> &[[f64; 2]] {
        &self.knots
    }
    /// Exact source lines and point receipts, in increasing caloric-interval order.
    pub fn source_lines(&self) -> &[PropertyLineSlice] {
        &self.source_lines
    }

    /// Evaluate without allocation, hashing, source selection or extrapolation.
    /// Interior knots select the right interval; the upper bound selects left.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<PhaseConductivitySample, PhaseConductivityError> {
        let bounds = self.enthalpy_bounds_j_kg();
        if !h.is_finite() || !(bounds[0]..=bounds[1]).contains(&h) {
            return Err(PhaseConductivityError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg: h,
                bounds_j_kg: bounds,
            });
        }
        let index = self
            .knots
            .partition_point(|knot| knot[0] <= h)
            .min(self.knots.len() - 1)
            - 1;
        let [a, ka] = self.knots[index];
        let [b, kb] = self.knots[index + 1];
        let fraction = (h - a) / (b - a);
        let k = if h == a {
            ka
        } else if h == b {
            kb
        } else if fraction <= 0.5 {
            (kb - ka).mul_add(fraction, ka)
        } else {
            (ka - kb).mul_add(1.0 - fraction, kb)
        };
        let interval = self.intervals[index];
        if !k.is_finite() || k < 0.0 || (k == 0.0 && interval.conductivity_slope != 0.0) {
            return Err(PhaseConductivityError::Unrepresentable(
                "conductivity and its native implicit tangent are not representable",
            ));
        }
        let caloric = self.caloric.values_at_specific_enthalpy(h)?;
        Ok(PhaseConductivitySample {
            specific_enthalpy: h,
            temperature: caloric.temperature_k,
            liquid_fraction: caloric.liquid_mass_fraction,
            conductivity: k,
            temperature_slope: interval.temperature_slope,
            conductivity_slope: interval.conductivity_slope,
            source_line: interval.source_line,
        })
    }
}

impl ConductivityLaw for PreparedPhaseConductivity<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.state_at_specific_enthalpy(h)
            .map(|sample| sample.conductivity_w_m_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }
}
impl ImplicitConductionLaw for PreparedPhaseConductivity<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }
    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.enthalpy_bounds_j_kg()
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let sample = self
            .state_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(ImplicitConductionState {
            temperature_k: sample.temperature_k(),
            temperature_slope_kg_k_j: sample.temperature_enthalpy_derivative_kg_k_j(),
            conductivity_w_m_k: sample.conductivity_w_m_k(),
            conductivity_slope_kg_w_m_k_j: sample.conductivity_enthalpy_derivative_kg_w_m_k_j(),
        })
    }
    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        let sample = self
            .state_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(ImplicitConductionObservables {
            temperature_k: sample.temperature_k(),
            phase_mass_fraction: Some(sample.liquid_mass_fraction()),
        })
    }
}
