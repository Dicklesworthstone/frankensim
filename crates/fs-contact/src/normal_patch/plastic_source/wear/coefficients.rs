//! Original wear claims sampled at actual local pressure and irreversible speed.
//!
//! These coordinates only select an existing source value. They do not supply
//! a friction closure, fabricate material calibration, or infer a wear regime.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::ArchardLaw;

use super::{SourcePlasticArchardWear, SourcePlasticWearError, checkpoint, invalid};

/// Coordinates evolved by the original local contact/sliding owner.
/// At least one coordinate is required; all other source conditions stay fixed
/// except the existing explicit interface-temperature continuation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalWearAxes {
    /// Actual beginning cell pressure [Pa], with static pressure semantics.
    pub pressure: Option<String>,
    /// Actual unsigned irreversible sliding distance / duration [m/s].
    /// This requires the dimensional length/time quantity, not a strain rate.
    pub sliding_speed: Option<String>,
}

impl LocalWearAxes {
    pub(super) fn owned(&self) -> Vec<(&str, QuantitySpec)> {
        self.pressure
            .as_deref()
            .map(|name| {
                (
                    name,
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Pressure,
                        ValueForm::Static,
                    )),
                )
            })
            .into_iter()
            .chain(
                self.sliding_speed
                    .as_deref()
                    .map(|name| (name, QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0])))),
            )
            .collect()
    }

    pub(super) fn admit(
        &self,
        source: &ResolvedInterfaceStatePoint,
        fixed: &ResolvedInterfaceStatePoint,
    ) -> Result<(), SourcePlasticWearError> {
        let owned = self.owned();
        if owned.is_empty()
            || owned.iter().enumerate().any(|(index, (name, quantity))| {
                name.is_empty()
                    || fixed.state_coordinate(name).is_some()
                    || owned[..index].iter().any(|(other, _)| other == name)
                    || source.axis_quantities().get(*name) != Some(quantity)
                    || !source
                        .state_coordinate(name)
                        .is_some_and(|value| value.is_finite() && value >= 0.0)
            })
        {
            return Err(invalid(
                "local wear requires distinct typed pressure/speed coordinates outside the fixed normal query",
            ));
        }
        let is_owned = |name: &str| owned.iter().any(|(axis, _)| *axis == name);
        if !source
            .query_point()
            .iter()
            .filter(|(name, _)| !is_owned(name))
            .eq(fixed.query_point().iter())
            || !source
                .axis_quantities()
                .iter()
                .filter(|(name, _)| !is_owned(name))
                .eq(fixed.axis_quantities().iter())
        {
            return Err(invalid(
                "local wear must retain every original fixed interface coordinate and quantity",
            ));
        }
        Ok(())
    }
}

/// Actual coefficient queries used by one accepted removal step.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalWearCoefficientSamples {
    /// Exact original pinned A/B claims at each distinct local query.
    pub sources: Vec<ResolvedInterfaceStatePoint>,
    /// Original-cell indices into `sources`. Unloaded or non-sliding cells
    /// have no sample and make no source-coverage claim at zero pressure/speed.
    pub cell_source_indices: Vec<Option<usize>>,
}

pub(super) struct EvaluatedLocalCoefficients {
    pub receipt: LocalWearCoefficientSamples,
    pub laws: Vec<[ArchardLaw; 2]>,
}

impl EvaluatedLocalCoefficients {
    pub(super) fn coefficient(&self, side: usize, cell: usize) -> f64 {
        self.receipt.cell_source_indices[cell]
            .map_or(0.0, |slot| self.laws[slot][side].wear_coefficient)
    }
}

impl SourcePlasticArchardWear {
    pub(super) fn sample_local_coefficients(
        &self,
        pressure: &[f64],
        slip: &[f64],
        duration_s: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<Option<EvaluatedLocalCoefficients>, SourcePlasticWearError> {
        if duration_s.is_some_and(|duration| !duration.is_finite() || duration <= 0.0) {
            return Err(invalid("positive finite wear-step duration"));
        }
        let Some(axes) = &self.local_axes else {
            return Ok(None);
        };
        let mut sources = Vec::new();
        let mut laws = Vec::new();
        let mut cell_source_indices = Vec::with_capacity(pressure.len());
        let mut indices = BTreeMap::new();
        for (&p, &distance) in pressure.iter().zip(slip) {
            checkpoint(cx)?;
            if p == 0.0 || distance == 0.0 {
                cell_source_indices.push(None);
                continue;
            }
            let speed = if axes.sliding_speed.is_some() {
                let duration = duration_s.ok_or_else(|| {
                    invalid("speed-dependent wear requires the actual sliding interval")
                })?;
                let speed = distance / duration;
                if !speed.is_finite() || speed <= 0.0 {
                    return Err(invalid("representable positive actual sliding speed"));
                }
                speed
            } else {
                0.0
            };
            let key = (
                if axes.pressure.is_some() {
                    p.to_bits()
                } else {
                    0
                },
                speed.to_bits(),
            );
            let slot = if let Some(&slot) = indices.get(&key) {
                slot
            } else {
                let mut query = self.coefficients.query();
                for (axis, value) in [
                    (axes.pressure.as_deref(), p),
                    (axes.sliding_speed.as_deref(), speed),
                ] {
                    if let Some(axis) = axis {
                        query = query
                            .with(axis, value)
                            .map_err(|_| invalid("finite original local wear query"))?;
                    }
                }
                let source = self.coefficients.at_query_point(
                    self.contact.normal_source().elastic().interface().card(),
                    &query,
                )?;
                let mut local = self.laws;
                for side in 0..2 {
                    let k = source
                        .property_by_key(&self.coefficient_keys[side])
                        .ok_or_else(|| invalid("original local wear coefficient disappeared"))?
                        .value_si();
                    if !k.is_finite() || k < 0.0 {
                        return Err(invalid("nonnegative finite local Archard coefficient"));
                    }
                    local[side].wear_coefficient = k;
                }
                let slot = sources.len();
                sources.push(source);
                laws.push(local);
                indices.insert(key, slot);
                slot
            };
            cell_source_indices.push(Some(slot));
        }
        Ok(Some(EvaluatedLocalCoefficients {
            receipt: LocalWearCoefficientSamples {
                sources,
                cell_source_indices,
            },
            laws,
        }))
    }
}
