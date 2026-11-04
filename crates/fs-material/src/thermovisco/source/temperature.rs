//! Original-card calibration for exact thermoviscoelastic temperature kinetics.
//!
//! Parameter roles and SI quantities are explicit. The selected original claims
//! must cover the whole mechanical interval at the same fixed context, and their
//! complete source-knot union must give constant calibration coefficients. This
//! preserves the calibration's reference temperature and source receipts instead
//! of synthesizing new material claims or fitting a sampled shift curve.

use super::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ResolvedThermoMaxwell,
    ScalarAdmissibility, ScalarPropertyRequirement, SourceError, poll,
    resolve_material_state_point,
};
use crate::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use std::collections::BTreeSet;

/// Explicit original property keys; no material-name calibration is inferred.
#[derive(Clone, Copy, Debug)]
pub enum TemperatureShiftSourceQuery<'a> {
    /// `log10(a_T) = -C1 (T-T_ref)/(C2+T-T_ref)`.
    WilliamsLandelFerry {
        /// Absolute reference temperature [K].
        reference_temperature: &'a PropertyKey,
        /// Positive, dimensionless base-ten coefficient.
        c1: &'a PropertyKey,
        /// Positive temperature DIFFERENCE [K], not an absolute temperature.
        c2: &'a PropertyKey,
    },
    /// `ln(a_T) = (E/R) (1/T-1/T_ref)`.
    Arrhenius {
        /// Absolute reference temperature [K].
        reference_temperature: &'a PropertyKey,
        /// Positive molar activation energy [J/mol], not energy per particle.
        activation_energy: &'a PropertyKey,
    },
}

/// Immutable calibration and complete original parameter-source receipts.
#[derive(Clone, Debug)]
pub struct ResolvedTemperatureShift {
    declaration: TemperatureShift,
    samples: Vec<ResolvedMaterialStatePoint>,
    // Role order: reference temperature, then C1/C2 or molar activation energy.
    keys: Vec<PropertyKey>,
}
impl ResolvedTemperatureShift {
    #[must_use]
    pub const fn declaration(&self) -> &TemperatureShift {
        &self.declaration
    }

    /// Complete calibration source-knot union, including both interval ends and
    /// the reference temperature, at the original fixed material context.
    #[must_use]
    pub fn samples(&self) -> &[ResolvedMaterialStatePoint] {
        &self.samples
    }

    /// Ordered original role keys: reference temperature, then C1/C2 or E.
    #[must_use]
    pub fn keys(&self) -> &[PropertyKey] {
        &self.keys
    }
}

impl ResolvedThermoMaxwell {
    /// Attach calibrated temperature kinetics from this SAME original material
    /// card, without changing the mechanical or caloric source owner.
    ///
    /// All parameter claims must continuously cover the existing mechanical
    /// temperature envelope, at exactly its nonthermal coordinates and quantity
    /// conventions. Every original source breakpoint is retained and checked;
    /// a coefficient varying even at an interior breakpoint refuses. Explicit
    /// semantic absolute/difference temperatures are accepted in their own roles.
    ///
    /// The existing relaxation multiplier must be constant; its original value
    /// remains the reference time scale. A previously tabulated temperature
    /// response or a second analytic shift refuses rather than being multiplied
    /// silently. No claim, branch history, modulus or heat-capacity law is fitted,
    /// reset or replaced.
    pub fn with_temperature_shift_from_source(
        mut self,
        card: &MaterialCard,
        query: TemperatureShiftSourceQuery<'_>,
        selection: SelectionPolicy,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SourceError> {
        poll(&mut cancelled)?;
        if card.content_hash() != self.samples[0].card_identity() {
            return Err(SourceError::Invalid(
                "temperature-shift calibration requires the original mechanical material card",
            ));
        }
        if self.law.temperature_shift().is_some() {
            return Err(SourceError::Invalid(
                "temperature kinetics are already attached",
            ));
        }
        let (keys, quantities): (Vec<PropertyKey>, Vec<(Dims, Option<QuantityKind>)>) = match query
        {
            TemperatureShiftSourceQuery::WilliamsLandelFerry {
                reference_temperature,
                c1,
                c2,
            } => (
                vec![reference_temperature.clone(), c1.clone(), c2.clone()],
                vec![
                    (Temperature::DIMS, Some(QuantityKind::AbsoluteTemperature)),
                    (Dims([0; 6]), None),
                    (Temperature::DIMS, Some(QuantityKind::TemperatureDifference)),
                ],
            ),
            TemperatureShiftSourceQuery::Arrhenius {
                reference_temperature,
                activation_energy,
            } => (
                vec![reference_temperature.clone(), activation_energy.clone()],
                vec![
                    (Temperature::DIMS, Some(QuantityKind::AbsoluteTemperature)),
                    (Dims([2, 1, -2, 0, 0, -1]), None),
                ],
            ),
        };
        let mut names = BTreeSet::new();
        let mut requirements = Vec::with_capacity(keys.len());
        for (key, (dims, kind)) in keys.iter().zip(quantities) {
            let quantity = key.quantity();
            let semantic_match = kind.is_some_and(|kind| {
                quantity == QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
            });
            if !names.insert(key.name())
                || (quantity != QuantitySpec::dimensional(dims) && !semantic_match)
            {
                return Err(SourceError::Invalid(
                    "distinct calibration roles with explicit reference/difference temperature or molar-energy quantities are required",
                ));
            }
            requirements.push(
                ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
                    .map_err(SourceError::Material)?,
            );
        }
        let axis = self.temperature_axis.as_str();
        let lower = self.samples[0].query();
        let upper = self.samples[self.samples.len() - 1].query();
        let low = self.law.knots[0].temperature_k;
        let high = self.law.knots[self.law.knots.len() - 1].temperature_k;
        let mut grid = vec![low, high];
        let mut pins = Vec::with_capacity(keys.len());
        for key in &keys {
            poll(&mut cancelled)?;
            let support = card
                .claims()
                .query_line_typed(key, axis, &lower, &upper, ClaimSelection::Policy(selection))
                .map_err(|e| SourceError::Support {
                    property: key.name().to_owned(),
                    reason: e.to_string(),
                })?;
            let selected = support.knots[0].1.receipt.selected;
            let claim = card.claims().claim(selected).ok_or(SourceError::Invalid(
                "selected temperature-shift parameter claim is absent",
            ))?;
            if claim
                .value
                .axes()
                .any(|(name, dims)| name == axis && dims != Temperature::DIMS)
            {
                return Err(SourceError::Invalid(
                    "calibration temperature coordinate must have temperature dimensions",
                ));
            }
            for (index, (temperature_k, _)) in support.knots.into_iter().enumerate() {
                if index % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                if temperature_k > low && temperature_k < high {
                    insert_temperature(&mut grid, temperature_k)?;
                }
            }
            pins.push((key.name().to_owned(), selected));
        }
        let resolve = |temperature_k: f64| -> Result<ResolvedMaterialStatePoint, SourceError> {
            let point =
                lower
                    .clone()
                    .with(axis, temperature_k)
                    .map_err(|e| SourceError::Support {
                        property: axis.to_owned(),
                        reason: e.to_string(),
                    })?;
            let sample = resolve_material_state_point(
                card,
                &point,
                &requirements,
                MaterialPropertySelection::PinnedByProperty(pins.clone()),
            )
            .map_err(SourceError::Material)?;
            for property in sample.properties() {
                if !property.has_only_fixed_external_conditions(Some(axis), true) {
                    return Err(SourceError::Invalid(
                        "temperature-shift calibration may depend only on temperature and fixed external conditions",
                    ));
                }
                if property
                    .angular_frequency_band_rad_s()
                    .map_err(|e| SourceError::Support {
                        property: "temperature-shift calibration".into(),
                        reason: e.to_string(),
                    })?
                    .is_some_and(|(lo, hi)| lo > 0.0 || hi < 0.0)
                {
                    return Err(SourceError::Invalid(
                        "temperature-shift calibration must support the static/relaxed endpoint",
                    ));
                }
            }
            Ok(sample)
        };
        let first = resolve(low)?;
        let values = calibration_values(&first, &keys)?;
        let model = match query {
            TemperatureShiftSourceQuery::WilliamsLandelFerry { .. } => {
                TemperatureShiftModel::WilliamsLandelFerry {
                    c1: values[1],
                    c2_k: values[2],
                }
            }
            TemperatureShiftSourceQuery::Arrhenius { .. } => TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: values[1],
            },
        };
        let declaration = TemperatureShift::new(model, values[0], [low, high]).map_err(|e| {
            SourceError::Support {
                property: "temperature-shift calibration".into(),
                reason: e.to_string(),
            }
        })?;
        insert_temperature(&mut grid, declaration.reference_temperature_k())?;
        let mut samples = Vec::with_capacity(grid.len());
        samples.push(first);
        for &temperature_k in &grid[1..] {
            poll(&mut cancelled)?;
            let sample = resolve(temperature_k)?;
            if calibration_values(&sample, &keys)? != values {
                return Err(SourceError::Invalid(
                    "temperature-shift calibration parameters vary over the original source envelope",
                ));
            }
            samples.push(sample);
        }
        poll(&mut cancelled)?;
        self.law = self
            .law
            .with_temperature_shift(declaration)
            .map_err(SourceError::Mechanical)?;
        self.temperature_shift_source = Some(ResolvedTemperatureShift {
            declaration,
            samples,
            keys,
        });
        Ok(self)
    }

    /// Original shift calibration, separate from the unchanged mechanical
    /// coefficient bundles returned by `samples()`.
    #[must_use]
    pub const fn temperature_shift_source(&self) -> Option<&ResolvedTemperatureShift> {
        self.temperature_shift_source.as_ref()
    }
}

fn calibration_values(
    sample: &ResolvedMaterialStatePoint,
    keys: &[PropertyKey],
) -> Result<Vec<f64>, SourceError> {
    keys.iter()
        .map(|key| {
            sample
                .property_by_key(key)
                .map(|property| property.value_si())
                .ok_or(SourceError::Invalid(
                    "original calibration property is absent",
                ))
        })
        .collect()
}

fn insert_temperature(grid: &mut Vec<f64>, temperature_k: f64) -> Result<(), SourceError> {
    if let Err(index) = grid.binary_search_by(|value| value.total_cmp(&temperature_k)) {
        if grid.len() >= 1024 {
            return Err(SourceError::Invalid(
                "complete calibration source grid exceeds 1024 knots",
            ));
        }
        grid.insert(index, temperature_k);
    }
    Ok(())
}
