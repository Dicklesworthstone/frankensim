//! Explicit, bounded reconciliation of independent NASA-9 caloric references.
//!
//! Cp remains source-authored. Constant h/s offsets make contiguous charts
//! agree at their shared endpoints; u=h-RT and g=h-Ts follow the same owner.
//! The chosen anchor region is unchanged. This is an opt-in model adjustment,
//! not a claim that independently fitted source cards were already continuous.

use crate::{
    MolarEnergyQuantityV1, MolarEnthalpyV1, MolarEntropyV1, MolarThermalQuantityV1,
    Nasa9PropertyV1, Nasa9StandardStateModelV1, StandardStateMolarPropertiesV1, ThermochemErrorV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};
use fs_qty::Temperature;
use std::sync::Arc;

/// Operation-tree version for explicit continuous caloric references.
pub const NASA9_CONTINUOUS_CALORIC_EVALUATOR_VERSION_V1: u32 = 2;

/// One original source identity and its explicit constant caloric offsets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Nasa9CaloricRegionCorrectionV1 {
    card_identity: [u8; 32],
    enthalpy_offset_bits: u64,
    entropy_offset_bits: u64,
}
impl Nasa9CaloricRegionCorrectionV1 {
    /// Original immutable NASA-9 source card; its coefficients are unchanged.
    pub const fn card_identity(&self) -> &[u8; 32] {
        &self.card_identity
    }
    /// Constant added to source enthalpy throughout this region, J/mol.
    pub fn enthalpy_offset(&self) -> MolarEnthalpyV1 {
        MolarEnthalpyV1::new(f64::from_bits(self.enthalpy_offset_bits))
    }
    /// Constant added to source entropy throughout this region, J/(mol K).
    pub fn entropy_offset(&self) -> MolarEntropyV1 {
        MolarEntropyV1::new(f64::from_bits(self.entropy_offset_bits))
    }
}

/// Declared anchor, correction budgets, original identities, and exact offsets.
///
/// The model and accepted receipts share this immutable bounded metadata.
/// Budgets bound the adjustment to each raw source region, not source-data error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Nasa9CaloricReferenceV1 {
    anchor_temperature_bits: u64,
    anchor_region: usize,
    maximum_enthalpy_correction_bits: u64,
    maximum_entropy_correction_bits: u64,
    regions: Vec<Nasa9CaloricRegionCorrectionV1>,
}
impl Nasa9CaloricReferenceV1 {
    /// Version of the complete corrected evaluation operation tree.
    pub const fn evaluator_version(&self) -> u32 {
        NASA9_CONTINUOUS_CALORIC_EVALUATOR_VERSION_V1
    }
    /// Declared temperature selecting the region whose h/s stay unchanged.
    pub fn anchor_temperature(&self) -> Temperature {
        Temperature::new(f64::from_bits(self.anchor_temperature_bits))
    }
    /// Original source region whose caloric reference is preserved.
    pub const fn anchor_region(&self) -> usize {
        self.anchor_region
    }
    /// Maximum absolute adjustment to any source enthalpy, J/mol.
    pub fn maximum_enthalpy_correction(&self) -> MolarEnergyQuantityV1 {
        MolarEnergyQuantityV1::new(f64::from_bits(self.maximum_enthalpy_correction_bits))
    }
    /// Maximum absolute adjustment to any source entropy, J/(mol K).
    pub fn maximum_entropy_correction(&self) -> MolarThermalQuantityV1 {
        MolarThermalQuantityV1::new(f64::from_bits(self.maximum_entropy_correction_bits))
    }
    /// Corrections in the exact original source-region order.
    pub fn regions(&self) -> &[Nasa9CaloricRegionCorrectionV1] {
        &self.regions
    }
}

fn refused(region: usize, reason: &'static str) -> ThermochemErrorV1 {
    ThermochemErrorV1::CaloricReferenceRefused { region, reason }
}

fn shifted(value: f64, offset: f64) -> f64 {
    if offset == 0.0 { value } else { value + offset }
}

impl Nasa9StandardStateModelV1 {
    /// Explicitly reconcile contiguous source h/s charts to one anchor region.
    ///
    /// The original cards and Cp polynomials are retained byte for byte.
    /// All h/s values in the anchor region remain unchanged; other regions
    /// receive constant offsets. Usually the anchor region contains the
    /// formation/reference temperature. Anchoring a different region changes
    /// low-temperature formation values within the caller's explicit budget.
    ///
    /// Every cumulative offset is checked against its declared molar bound.
    /// Gaps and unrepresentable one-sided endpoint agreement refuse. This does
    /// not prove positive Cp/Cv; consumers still validate caloric stability.
    /// Repeated calls always derive from original cards, never compound offsets.
    /// Construction uses at most two source evaluations per seam (16 regions).
    pub fn with_caloric_continuity(
        &self,
        anchor_temperature: Temperature,
        maximum_enthalpy_correction: MolarEnergyQuantityV1,
        maximum_entropy_correction: MolarThermalQuantityV1,
    ) -> Result<Self, ThermochemErrorV1> {
        let t = anchor_temperature.value();
        if !t.is_finite() || t <= 0.0 {
            return Err(ThermochemErrorV1::InvalidEvaluationTemperature { bits: t.to_bits() });
        }
        let anchor = self
            .region_for_temperature(t)
            .ok_or(ThermochemErrorV1::TemperatureOutOfRange { bits: t.to_bits() })?
            .0;
        let maximum_h = maximum_enthalpy_correction.value();
        let maximum_s = maximum_entropy_correction.value();
        if !maximum_h.is_finite() || maximum_h < 0.0 || !maximum_s.is_finite() || maximum_s < 0.0 {
            return Err(refused(
                anchor,
                "caloric correction budgets must be finite and nonnegative",
            ));
        }
        for (index, pair) in self.regions.windows(2).enumerate() {
            if pair[0].temperature_max().value().to_bits()
                != pair[1].temperature_min().value().to_bits()
            {
                return Err(refused(
                    index + 1,
                    "continuous caloric references require contiguous regions",
                ));
            }
        }
        let mut corrections: Vec<_> = self
            .regions
            .iter()
            .map(|region| Nasa9CaloricRegionCorrectionV1 {
                card_identity: *region.card_identity(),
                enthalpy_offset_bits: 0.0_f64.to_bits(),
                entropy_offset_bits: 0.0_f64.to_bits(),
            })
            .collect();
        // Walk outward from the unchanged anchor. Each seam is visited once.
        for index in (anchor + 1..self.regions.len()).chain((0..anchor).rev()) {
            let neighbor = if index > anchor { index - 1 } else { index + 1 };
            let boundary = if index > anchor {
                self.regions[index].temperature_min().value()
            } else {
                self.regions[index].temperature_max().value()
            };
            let original = self.evaluate_source_region(index, boundary)?;
            let adjacent = self.evaluate_source_region(neighbor, boundary)?;
            let h = shifted(
                adjacent.h().value(),
                corrections[neighbor].enthalpy_offset().value(),
            );
            let s = shifted(
                adjacent.s().value(),
                corrections[neighbor].entropy_offset().value(),
            );
            let dh = h - original.h().value();
            let ds = s - original.s().value();
            if !dh.is_finite() || dh.abs() > maximum_h {
                return Err(refused(
                    index,
                    "enthalpy correction exceeds the declared molar bound",
                ));
            }
            if !ds.is_finite() || ds.abs() > maximum_s {
                return Err(refused(
                    index,
                    "entropy correction exceeds the declared molar bound",
                ));
            }
            if shifted(original.h().value(), dh).to_bits() != h.to_bits()
                || shifted(original.s().value(), ds).to_bits() != s.to_bits()
            {
                return Err(refused(
                    index,
                    "caloric seam agreement is not representable",
                ));
            }
            corrections[index].enthalpy_offset_bits = dh.to_bits();
            corrections[index].entropy_offset_bits = ds.to_bits();
        }
        let mut reconciled = self.clone();
        reconciled.caloric_reference = Some(Arc::new(Nasa9CaloricReferenceV1 {
            anchor_temperature_bits: t.to_bits(),
            anchor_region: anchor,
            maximum_enthalpy_correction_bits: maximum_h.to_bits(),
            maximum_entropy_correction_bits: maximum_s.to_bits(),
            regions: corrections,
        }));
        Ok(reconciled)
    }

    /// Explicit reconciliation policy, if this model was opted into continuity.
    pub fn caloric_reference(&self) -> Option<&Nasa9CaloricReferenceV1> {
        self.caloric_reference.as_deref()
    }

    pub(crate) fn apply_caloric_reference(
        &self,
        region: usize,
        t: f64,
        mut properties: StandardStateMolarPropertiesV1,
    ) -> Result<StandardStateMolarPropertiesV1, ThermochemErrorV1> {
        let Some(reference) = self.caloric_reference() else {
            return Ok(properties);
        };
        let correction = &reference.regions[region];
        let dh = correction.enthalpy_offset().value();
        let ds = correction.entropy_offset().value();
        if dh == 0.0 && ds == 0.0 {
            return Ok(properties);
        }
        let h = shifted(properties.h().value(), dh);
        let s = shifted(properties.s().value(), ds);
        let u = h - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t;
        let g = h - t * s;
        for (property, value) in [
            (Nasa9PropertyV1::Enthalpy, h),
            (Nasa9PropertyV1::Entropy, s),
            (Nasa9PropertyV1::InternalEnergy, u),
            (Nasa9PropertyV1::GibbsEnergy, g),
        ] {
            if !value.is_finite() {
                return Err(ThermochemErrorV1::NonFiniteEvaluation {
                    property,
                    bits: value.to_bits(),
                });
            }
        }
        properties.h = crate::MolarEnthalpyV1::new(h);
        properties.s = crate::MolarEntropyV1::new(s);
        properties.u = crate::MolarInternalEnergyV1::new(u);
        properties.g = crate::MolarGibbsEnergyV1::new(g);
        Ok(properties)
    }

    // Relative caloric changes integrate the retained Cp chart on each segment.
    // This preserves small changes across a seam instead of subtracting its
    // independently rounded offsets. Work is bounded by the admitted 16 regions.
    pub(crate) fn integrate_caloric_regions(
        &self,
        first: usize,
        last: usize,
        lo: f64,
        hi: f64,
        mut integrate: impl FnMut([f64; 9], f64, f64) -> f64,
    ) -> f64 {
        let mut sum = 0.0;
        let mut correction = 0.0;
        for index in first..=last {
            let region = &self.regions[index];
            let start = lo.max(region.temperature_min().value());
            let end = hi.min(region.temperature_max().value());
            if start < end {
                crate::entropy::accumulate(
                    &mut sum,
                    &mut correction,
                    integrate(region.coefficients, start, end),
                );
            }
        }
        sum
    }
}
