//! Source-backed generalized-Maxwell coefficients for the existing thermal law.
//!
//! The selected scalar, curve or sampled-surface claims must continuously cover
//! the same fixed material context. Their complete temperature-line knot union
//! is lowered without resampling, including every crossed source-cell edge.
//! Reference branch coefficients are constant; the three explicit thermal
//! responses are the common modulus multiplier, time multiplier and free strain.
//! Neither a material name nor an absent branch supplies an implicit property.

use super::thermal::{CaloricPoint, FeedbackError};
use super::{Limits, ThermalKnot, ThermoMaxwell};
use crate::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use crate::state_point::{
    MAX_MATERIAL_STATE_PROPERTIES, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature, Time};
use std::collections::BTreeSet;

/// Explicit roles and source context; branch order is never sorted or inferred.
#[derive(Clone, Copy, Debug)]
pub struct ThermoMaxwellSourceQuery<'a> {
    /// Equilibrium axial modulus [Pa].
    pub equilibrium_modulus: &'a PropertyKey,
    /// Reference branch (axial modulus [Pa], relaxation time [s]) keys.
    /// Up to 30 pairs fit the shared 64-property atomic-resolution budget.
    pub branches: &'a [(PropertyKey, PropertyKey)],
    /// Common multiplier of all spring moduli [1], explicitly sourced.
    pub modulus_scale: &'a PropertyKey,
    /// Common multiplier of all relaxation times [1], explicitly sourced.
    pub relaxation_scale: &'a PropertyKey,
    /// Free engineering strain [1], not an expansion coefficient [1/K].
    pub free_strain: &'a PropertyKey,
    /// Absolute-temperature coordinate; all other coordinates stay fixed.
    pub temperature_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    pub selection: SelectionPolicy,
}

/// Source or physical-law failure; construction never publishes a partial law.
#[derive(Debug)]
pub enum SourceError {
    Invalid(&'static str),
    Support { property: String, reason: String },
    Material(MaterialStatePointError),
    Mechanical(super::Error),
    Caloric(PhaseStateError),
    Feedback(FeedbackError),
    Cancelled,
}
impl core::fmt::Display for SourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "thermovisco source: {s}"),
            Self::Support { property, reason } => {
                write!(f, "thermovisco source {property}: {reason}")
            }
            Self::Material(e) => write!(f, "thermovisco source: {e}"),
            Self::Mechanical(e) => write!(f, "thermovisco source: {e}"),
            Self::Caloric(e) => write!(f, "thermovisco source: {e}"),
            Self::Feedback(e) => write!(f, "thermovisco source: {e}"),
            Self::Cancelled => f.write_str("thermovisco source: cancelled"),
        }
    }
}
impl std::error::Error for SourceError {}

/// An executable law and the exact original source bundles that produced it.
/// No new material registry or constitutive update is introduced. Fitted branch
/// data remain fitted data; source receipts do not imply empirical validation.
#[derive(Clone, Debug)]
pub struct ResolvedThermoMaxwell {
    law: ThermoMaxwell,
    samples: Vec<ResolvedMaterialStatePoint>,
    // Role order: equilibrium, each branch's E/tau, then scale/time/free strain.
    keys: Vec<PropertyKey>,
    temperature_axis: String,
}
impl ResolvedThermoMaxwell {
    /// Resolve one continuously supported temperature line at fixed conditions.
    ///
    /// Every source breakpoint inside the envelope is retained, including those
    /// absent from other properties' grids. All base E/tau coefficients must be
    /// constant across that complete union; temperature dependence must use the
    /// declared shared scale laws. This refuses independent branch-specific
    /// shifts rather than approximating them by one invented common multiplier.
    pub fn try_new(
        card: &MaterialCard,
        query: ThermoMaxwellSourceQuery<'_>,
        limits: Limits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SourceError> {
        poll(&mut cancelled)?;
        let [low, high] = span(query.lower, query.upper, query.temperature_axis)?;
        if query.branches.len() > (MAX_MATERIAL_STATE_PROPERTIES - 4) / 2 {
            return Err(SourceError::Invalid(
                "source branches exceed the atomic material-property budget",
            ));
        }
        let mut keys: Vec<PropertyKey> = vec![query.equilibrium_modulus.clone()];
        let mut dims = vec![Pressure::DIMS];
        for (modulus, time) in query.branches {
            keys.extend([modulus.clone(), time.clone()]);
            dims.extend([Pressure::DIMS, Time::DIMS]);
        }
        let base_count = keys.len();
        keys.extend([
            query.modulus_scale.clone(),
            query.relaxation_scale.clone(),
            query.free_strain.clone(),
        ]);
        dims.extend([Dims([0; 6]); 3]);
        let mut names = BTreeSet::new();
        let mut requirements = Vec::with_capacity(keys.len());
        for (index, (key, dimension)) in keys.iter().zip(dims).enumerate() {
            if !names.insert(key.name()) || key.quantity() != QuantitySpec::dimensional(dimension) {
                return Err(SourceError::Invalid(
                    "distinct role keys with the declared SI dimensions are required",
                ));
            }
            let domain = if index == keys.len() - 1 {
                ScalarAdmissibility::Finite
            } else if index > 0 && index < base_count && index % 2 == 1 {
                ScalarAdmissibility::NonNegative
            } else {
                ScalarAdmissibility::StrictlyPositive
            };
            requirements.push(
                ScalarPropertyRequirement::try_with_key(key, domain)
                    .map_err(SourceError::Material)?,
            );
        }
        let mut grid = vec![low, high];
        let mut pins = Vec::with_capacity(keys.len());
        for key in &keys {
            poll(&mut cancelled)?;
            let support = card
                .claims()
                .query_line_typed(
                    key,
                    query.temperature_axis,
                    query.lower,
                    query.upper,
                    ClaimSelection::Policy(query.selection),
                )
                .map_err(|e| SourceError::Support {
                    property: key.name().to_owned(),
                    reason: e.to_string(),
                })?;
            let selected = support.knots[0].1.receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or(SourceError::Invalid("selected source claim is absent"))?;
            if claim
                .value
                .axes()
                .any(|(axis, dims)| axis == query.temperature_axis && dims != Temperature::DIMS)
            {
                return Err(SourceError::Invalid(
                    "source temperature coordinate must have temperature dimensions",
                ));
            }
            for (index, (t, _)) in support.knots.into_iter().enumerate() {
                if index % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                if t <= low || t >= high {
                    continue;
                }
                if let Err(at) = grid.binary_search_by(|value| value.total_cmp(&t)) {
                    if grid.len() >= 1024 {
                        return Err(SourceError::Invalid(
                            "source-complete thermal grid exceeds 1024 knots",
                        ));
                    }
                    grid.insert(at, t);
                }
            }
            pins.push((key.name().to_owned(), selected));
        }
        let mut samples = Vec::with_capacity(grid.len());
        let mut knots = Vec::with_capacity(grid.len());
        let mut reference: Option<Vec<f64>> = None;
        for temperature_k in grid {
            poll(&mut cancelled)?;
            // Updating an existing coordinate retains its original quantity kind.
            let point = query
                .lower
                .clone()
                .with(query.temperature_axis, temperature_k)
                .map_err(|e| SourceError::Support {
                    property: query.temperature_axis.into(),
                    reason: e.to_string(),
                })?;
            let sample = resolve_material_state_point(
                card,
                &point,
                &requirements,
                MaterialPropertySelection::PinnedByProperty(pins.clone()),
            )
            .map_err(SourceError::Material)?;
            let mut values = Vec::with_capacity(keys.len());
            for key in &keys {
                let property = sample
                    .property(key.name())
                    .ok_or(SourceError::Invalid("required source value is absent"))?;
                if !property.has_only_fixed_external_conditions(Some(query.temperature_axis), true)
                {
                    return Err(SourceError::Invalid(
                        "coefficients may depend only on the owned temperature and fixed external conditions; mechanical or frequency value axes require an evolving owner",
                    ));
                }
                if property
                    .angular_frequency_band_rad_s()
                    .map_err(|e| SourceError::Support {
                        property: key.name().into(),
                        reason: format!("{e:?}"),
                    })?
                    .is_some_and(|(lo, hi)| lo > 0.0 || hi < 0.0)
                {
                    return Err(SourceError::Invalid(
                        "material coefficients do not support the static/relaxed endpoint",
                    ));
                }
                values.push(property.value_si());
            }
            if let Some(first) = &reference {
                if values[..base_count] != first[..] {
                    return Err(SourceError::Invalid(
                        "reference branch coefficients vary with temperature; declare a supported common scale law",
                    ));
                }
            } else {
                reference = Some(values[..base_count].to_vec());
            }
            knots.push(ThermalKnot {
                temperature_k,
                modulus_scale: values[base_count],
                relaxation_scale: values[base_count + 1],
                free_strain: values[base_count + 2],
            });
            samples.push(sample);
        }
        let reference = reference.ok_or(SourceError::Invalid("empty source envelope"))?;
        let branches = reference[1..]
            .chunks_exact(2)
            .map(|pair| (pair[0], pair[1]))
            .collect();
        let law = ThermoMaxwell::new(reference[0], branches, knots, limits)
            .map_err(SourceError::Mechanical)?;
        poll(&mut cancelled)?;
        Ok(Self {
            law,
            samples,
            keys,
            temperature_axis: query.temperature_axis.into(),
        })
    }

    /// Same law consumed by existing deposition, bridge and material-point solvers.
    #[must_use]
    pub const fn law(&self) -> &ThermoMaxwell {
        &self.law
    }
    /// Select the existing law's history integrator while retaining every
    /// source sample, original role key and material context. Select before
    /// creating points; the caloric prescribed-strain path uses this same law.
    #[must_use]
    pub fn with_history_integration(mut self, integration: super::HistoryIntegration) -> Self {
        self.law = self.law.with_history_integration(integration);
        self
    }
    /// Original card/claim/query receipts in increasing temperature order.
    #[must_use]
    pub fn samples(&self) -> &[ResolvedMaterialStatePoint] {
        &self.samples
    }
    /// Ordered role binding: equilibrium, branch E/tau pairs, modulus/time/free strain.
    #[must_use]
    pub fn keys(&self) -> &[PropertyKey] {
        &self.keys
    }

    /// Admit another resolved property bundle as the same physical material state.
    /// Only the declared temperature coordinate may change, within this law's
    /// supported interval. Card identity, other coordinates and their quantity
    /// conventions remain exact. This does not reselect sources or certify the
    /// other bundle's constitutive law; its original resolver remains responsible.
    pub fn admit_state_point(&self, state: &ResolvedMaterialStatePoint) -> Result<(), SourceError> {
        let reference = &self.samples[0];
        if state.card_identity() != reference.card_identity()
            || state.axis_quantities() != reference.axis_quantities()
            || !same_context(
                state.query_point(),
                reference.query_point(),
                &self.temperature_axis,
            )
        {
            return Err(SourceError::Invalid(
                "property bundle and mechanical law must share the same card and fixed context",
            ));
        }
        let temperature = state
            .query_point()
            .iter()
            .find(|(key, _)| key == &self.temperature_axis)
            .map(|(_, value)| *value)
            .ok_or(SourceError::Invalid("property bundle lacks temperature"))?;
        let low = self.law.knots[0].temperature_k;
        let high = self.law.knots[self.law.knots.len() - 1].temperature_k;
        if !temperature.is_finite() || !(low..=high).contains(&temperature) {
            return Err(SourceError::Invalid(
                "property bundle temperature is outside the mechanical source envelope",
            ));
        }
        Ok(())
    }

    fn validate_caloric_source(
        &self,
        caloric: &ResolvedHeatCapacityCurve,
    ) -> Result<(), SourceError> {
        let reference = &self.samples[0];
        let low = self.law.knots[0].temperature_k;
        let high = self.law.knots[self.law.knots.len() - 1].temperature_k;
        if caloric.curve().material_card_identity() != reference.card_identity()
            || caloric
                .curve()
                .knots()
                .iter()
                .any(|k| k.temperature_k < low || k.temperature_k > high)
            || caloric.receipts().iter().any(|r| {
                r.axis_quantities != *reference.axis_quantities()
                    || !same_context(
                        &r.query_point,
                        reference.query_point(),
                        &self.temperature_axis,
                    )
            })
        {
            return Err(SourceError::Invalid(
                "caloric and mechanical models must share a card, fixed context and temperature support",
            ));
        }
        Ok(())
    }

    /// Reference density [kg/m3] from a caloric bundle proven to belong to this
    /// exact mechanical source context. This is the mass authority for closed
    /// deposited specimens; callers must not pair the law with an unrelated
    /// hand-authored density when this shared bundle is available.
    pub fn reference_density_kg_m3(
        &self,
        caloric: &ResolvedHeatCapacityCurve,
    ) -> Result<f64, SourceError> {
        self.validate_caloric_source(caloric)?;
        let density = caloric.curve().knots()[0].bulk_density_kg_m3;
        if !density.is_finite()
            || density <= 0.0
            || caloric
                .curve()
                .knots()
                .iter()
                .any(|k| k.bulk_density_kg_m3 != density)
        {
            return Err(SourceError::Invalid(
                "source-bound reference density must be finite, positive and constant over the admitted thermal span",
            ));
        }
        Ok(density)
    }

    /// Initialize heat feedback from the SAME material card and fixed context.
    /// Caloric interpolation may introduce extra knots but must stay inside the
    /// mechanical span and retain identical nonthermal coordinates and conventions.
    /// The caloric owner admits globally affine common spring and free-strain
    /// laws with mechanical entropy and reversible heat. Expansion requires
    /// adequate total caloric capacity; density changes, derivative jumps and
    /// phase change still refuse. Every original source knot is retained.
    pub fn caloric_point<'a>(
        &'a self,
        caloric: &'a ResolvedHeatCapacityCurve,
        temperature_k: f64,
        total_strain: f64,
        volume_m3: f64,
        surface_area_m2: f64,
    ) -> Result<CaloricPoint<'a>, SourceError> {
        let density = self.reference_density_kg_m3(caloric)?;
        let state = caloric
            .state_at_temperature(temperature_k)
            .map_err(SourceError::Caloric)?;
        let point = self
            .law
            .birth(state.temperature_k(), total_strain)
            .map_err(SourceError::Mechanical)?;
        let mass = volume_m3 * density;
        CaloricPoint::try_new(
            point,
            caloric.curve(),
            state.specific_enthalpy_j_kg(),
            mass,
            volume_m3,
            surface_area_m2,
        )
        .map_err(SourceError::Feedback)
    }
}

fn same_context(a: &[(String, f64)], b: &[(String, f64)], axis: &str) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|((ka, va), (kb, vb))| ka == kb && (ka == axis || va.to_bits() == vb.to_bits()))
}
fn span(lower: &QueryPoint, upper: &QueryPoint, axis: &str) -> Result<[f64; 2], SourceError> {
    let low = lower
        .axes()
        .get(axis)
        .copied()
        .ok_or(SourceError::Invalid("lower query lacks temperature"))?;
    let high = upper
        .axes()
        .get(axis)
        .copied()
        .ok_or(SourceError::Invalid("upper query lacks temperature"))?;
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if !low.is_finite()
        || !high.is_finite()
        || low <= 0.0
        || low >= high
        || lower.axis_quantities() != upper.axis_quantities()
        || lower.axes().len() != upper.axes().len()
        || lower.axes().iter().any(|(name, value)| {
            upper
                .axes()
                .get(name)
                .is_none_or(|other| name != axis && value.to_bits() != other.to_bits())
        })
        || lower
            .axis_quantities()
            .get(axis)
            .is_some_and(|q| *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS))
    {
        return Err(SourceError::Invalid(
            "one increasing positive absolute-temperature axis and unchanged source context are required",
        ));
    }
    Ok([low, high])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SourceError> {
    if cancelled() {
        Err(SourceError::Cancelled)
    } else {
        Ok(())
    }
}
