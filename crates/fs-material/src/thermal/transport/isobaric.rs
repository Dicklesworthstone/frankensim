//! Source Cp(T) on the existing fixed-mass, constant-pressure heat timeline.
//!
//! This caloric adapter needs neither conductivity nor a density-temperature law.
//! The original thermal owner evolves m*h and exchanges heat. Geometry remains
//! caller-authored heat-transfer geometry; this law supplies no current volume,
//! pressure evolution, phase change, mechanical storage or separate p*dV account.
//! The explicit [`expansion`] composition adds free geometry and pressure work
//! when compatible original expansion and density sources support this isobar.
//! The caller selects true single-phase isobaric Cp at fixed composition. A
//! same-dimension Cv, apparent capacity or strained-body capacity cannot substitute.
//!
//! The shared capacity kernel owns integration, entropy and inversion. Original
//! source knots, receipts and uncertainty survive; no resampled caloric chart,
//! second integrator or material-name branch is introduced.

/// Source-exact free expansion with explicit constant-pressure boundary work.
#[path = "isobaric/expansion.rs"]
pub mod expansion;

/// Source-exact Cp(T) and conductivity on native implicit conduction networks.
#[path = "isobaric/conduction.rs"]
pub mod conduction;

/// Exact original phase heat capacities on the finite-rate crystallization owner.
#[path = "isobaric/crystallization.rs"]
pub mod crystallization;

use super::{
    SPECIFIC_HEAT, ThermalTransportError, admit_temperature_point,
    capacity::{CapacityCurve, CapacityInverse},
    positive, source_line,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::{Pressure, QuantitySpec};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};

/// Explicit Cp role and one source-supported isobar at fixed composition.
#[derive(Clone, Copy, Debug)]
pub struct IsobaricHeatCapacityQuery<'a> {
    pub specific_heat_cp: &'a PropertyKey,
    pub temperature_axis: &'a str,
    /// Positive absolute pressure, typed in pascals and explicitly bounded by
    /// the selected source. A free extra query coordinate proves no support.
    pub pressure_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    pub selection: ClaimSelection,
}

/// Immutable single-phase caloric law. Zero h is at the caller's datum T.
/// Source Cp alone does not identify an absolute formation-energy datum.
/// Materials exchanging mass must use compatible externally established datums.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricHeatCapacityLaw {
    material_card: ContentHash,
    endpoints: [QueryPoint; 2],
    temperature_axis: String,
    pressure_axis: String,
    source: PropertyLineSlice,
    knots: Vec<[f64; 2]>,
    pressure_pa: f64,
    datum_temperature_k: f64,
    enthalpy_span: [f64; 2],
    inverse: Option<CapacityInverse>,
    maximum_slope: f64,
    identity: ContentHash,
}

/// One source-exact passive exchange with a fixed-temperature reservoir.
/// At a material-support endpoint, derivatives are admitted-side tangents;
/// they do not authorize perturbations beyond the original caloric source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearHeatExchange {
    /// Final absolute material temperature [K], rounded only after inversion.
    pub temperature_k: f64,
    /// Signed integral Cp dT [J/kg], retaining the solved relative temperature
    /// increment before rounding its absolute endpoint. Positive means heating.
    pub specific_enthalpy_change_j_kg: f64,
    /// Conductance-weighted mean material temperature [K]. With positive
    /// exposure X, delta_h = X * (reservoir - mean_temperature) up to rounding.
    /// Computed directly for weak exchange, not by subtracting delta_h/X from
    /// the reservoir. At zero exchange this is the initial temperature.
    pub mean_temperature_k: f64,
    /// Analytic d h_out / d h_in at fixed reservoir temperature, exposure and
    /// original Cp source. Nonnegative and at most one; a sufficiently strong
    /// exchange may underflow to zero. At equilibrium this is exp(-exposure/Cp),
    /// including when the absolute outlet rounds to the reservoir temperature.
    pub inlet_enthalpy_gain: f64,
    /// Stable complement of `inlet_enthalpy_gain`. Retained independently so
    /// the derivative of transferred heat survives very weak exchange.
    pub inlet_enthalpy_loss: f64,
    /// d h_out / d T_reservoir [J/(kg K)], with initial state, exposure and
    /// source fixed. Integrates original Cp along the actual fluid path; the
    /// reservoir itself need not lie inside material temperature support.
    pub reservoir_enthalpy_derivative_j_kg_k: f64,
    /// d T_mean / d T_reservoir, with initial state and exposure fixed.
    /// A direct positive moment retains weak exchange without subtracting
    /// nearly equal heat and exposure. Zero exposure has zero wall response.
    pub mean_reservoir_derivative: f64,
}

impl IsobaricHeatCapacityLaw {
    /// Admit the complete typed Cp path and source pressure, preserving all
    /// breakpoints. Cp must be positive throughout every affine segment.
    /// Pressure, composition and all other source coordinates stay fixed.
    pub fn try_new(
        card: &MaterialCard,
        query: IsobaricHeatCapacityQuery<'_>,
        datum_temperature_k: f64,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let key = query.specific_heat_cp;
        if key.quantity() != QuantitySpec::dimensional(SPECIFIC_HEAT)
            || key.hardness_test().is_some()
            || key.elastic_component().is_some()
            || key.strain_component().is_some()
            || key.stress_component().is_some()
            || query.temperature_axis == query.pressure_axis
        {
            return Err(invalid(
                "isobaric calorimetry requires a scalar static specific Cp key",
            ));
        }
        for point in [query.lower, query.upper] {
            admit_temperature_point(point, query.temperature_axis)?;
            if !point.axis_quantities().get(query.pressure_axis).is_some_and(|quantity| {
                *quantity == QuantitySpec::dimensional(Pressure::DIMS)
                    || *quantity == QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Pressure, ValueForm::Static))
            })
                || !point
                    .axes()
                    .get(query.pressure_axis)
                    .is_some_and(|&p| positive(p))
            {
                return Err(invalid(
                    "isobaric Cp requires positive static pressure in pascals, not a stress or dynamic alias",
                ));
            }
        }
        if query.lower.axis_quantities() != query.upper.axis_quantities()
            || !query.lower.axes().keys().eq(query.upper.axes().keys())
            || query.lower.axes().iter().any(|(axis, value)| {
                axis != query.temperature_axis && query.upper.axes()[axis] != *value
            })
            || query.lower.axes()[query.temperature_axis]
                >= query.upper.axes()[query.temperature_axis]
        {
            return Err(invalid(
                "the Cp interval must retain its complete fixed source context",
            ));
        }
        let source = source_line(
            card,
            key,
            query.temperature_axis,
            query.lower,
            query.upper,
            query.selection,
        )?;
        let claim = card
            .claims()
            .claim(source.knots[0].1.receipt.selected)
            .ok_or(invalid("selected Cp source is missing"))?;
        if !claim.validity.bounds().contains_key(query.pressure_axis) {
            return Err(invalid(
                "selected Cp source must explicitly support the declared pressure",
            ));
        }
        let knots: Vec<_> = source
            .knots
            .iter()
            .map(|(t, answer)| [*t, answer.evidence.value.value])
            .collect();
        if knots.len() < 2
            || knots.len() > fs_matdb::MAX_PROPERTY_LINE_KNOTS
            || knots.iter().any(|[t, cp]| !positive(*t) || !positive(*cp))
            || knots.windows(2).any(|pair| pair[0][0] >= pair[1][0])
        {
            return Err(invalid(
                "source Cp requires bounded ordered temperatures and positive capacities",
            ));
        }
        let curve = CapacityCurve::new(&knots);
        curve.capacity(datum_temperature_k)?;
        let enthalpy_span = [
            curve.energy_change(datum_temperature_k, knots[0][0])?,
            curve.energy_change(datum_temperature_k, knots[knots.len() - 1][0])?,
        ];
        let minimum = knots.iter().map(|k| k[1]).fold(f64::INFINITY, f64::min);
        let maximum_slope = (1.0 / minimum).next_up();
        if !positive(maximum_slope)
            || !enthalpy_span.iter().all(|h| h.is_finite())
            || enthalpy_span[0] >= enthalpy_span[1]
        {
            return Err(invalid(
                "isobaric caloric range or inverse slope is unrepresentable",
            ));
        }
        let pressure_pa = query.lower.axes()[query.pressure_axis];
        let inverse = CapacityInverse::new(&curve, datum_temperature_k).ok();
        let mut hash = DomainHasher::new("org.frankensim.fs-material.isobaric-cp-law.v3");
        hash.update(card.content_hash().as_bytes());
        for axis in [query.temperature_axis, query.pressure_axis] {
            hash.update(&(axis.len() as u64).to_le_bytes());
            hash.update(axis.as_bytes());
        }
        hash.update(&(knots.len() as u64).to_le_bytes());
        for (_, answer) in &source.knots {
            hash.update(answer.receipt.content_hash().as_bytes());
        }
        hash.update(&pressure_pa.to_bits().to_le_bytes());
        hash.update(&datum_temperature_k.to_bits().to_le_bytes());
        Ok(Self {
            material_card: card.content_hash(),
            endpoints: [query.lower.clone(), query.upper.clone()],
            temperature_axis: query.temperature_axis.to_owned(),
            pressure_axis: query.pressure_axis.to_owned(),
            source,
            knots,
            pressure_pa,
            datum_temperature_k,
            enthalpy_span,
            inverse,
            maximum_slope,
            identity: hash.finalize(),
        })
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    /// Unmodified original answers, receipts and uncertainty at every source knot.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        &self.source
    }
    #[must_use]
    pub const fn query_endpoints(&self) -> &[QueryPoint; 2] {
        &self.endpoints
    }
    /// Original typed temperature coordinate of this source path.
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    /// Original source-supported static pressure coordinate of this isobar.
    #[must_use]
    pub fn pressure_axis(&self) -> &str {
        &self.pressure_axis
    }
    #[must_use]
    pub const fn pressure_pa(&self) -> f64 {
        self.pressure_pa
    }
    #[must_use]
    pub const fn datum_temperature_k(&self) -> f64 {
        self.datum_temperature_k
    }
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        [self.knots[0][0], self.knots[self.knots.len() - 1][0]]
    }
    #[must_use]
    pub const fn specific_enthalpy_span_j_kg(&self) -> [f64; 2] {
        self.enthalpy_span
    }

    pub fn specific_heat_cp_j_kg_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).capacity(temperature_k)
    }

    /// Integral Cp dT from the selected datum; no density or volume enters.
    /// Reuse the caloric index and one original affine segment, retaining the
    /// canonical support endpoints and direct integration at numerical edges.
    pub fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, ThermalTransportError> {
        let [cold, hot] = self.temperature_span_k();
        let [low, high] = self.enthalpy_span;
        if temperature_k == cold {
            return Ok(low);
        }
        if temperature_k == hot {
            return Ok(high);
        }
        let curve = CapacityCurve::new(&self.knots);
        let direct = || curve.energy_change(self.datum_temperature_k, temperature_k);
        let enthalpy = self.inverse.as_ref().map_or_else(
            direct,
            |inverse| inverse.energy_change_from_datum(&curve, temperature_k),
        )?;
        if enthalpy < low || enthalpy > high {
            return direct();
        }
        Ok(enthalpy)
    }

    /// Datum-aligned changes reuse the caloric index. Other changes integrate
    /// crossed source segments directly; no large reference energies are subtracted.
    pub fn specific_enthalpy_change_j_kg(
        &self,
        initial_k: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        if initial_k == self.datum_temperature_k {
            return self.specific_enthalpy_j_kg(final_k);
        }
        if final_k == self.datum_temperature_k {
            return self.specific_enthalpy_j_kg(initial_k).map(|h| -h);
        }
        CapacityCurve::new(&self.knots).energy_change(initial_k, final_k)
    }

    /// Reversible integral Cp/T dT on this same isobar, not absolute source entropy.
    pub fn specific_entropy_change_at_temperatures_j_kg_k(
        &self,
        initial_k: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).entropy_change(initial_k, final_k)
    }

    /// Same analytic capacity inverse for heating or cooling, without extrapolation.
    pub fn temperature_after_specific_heat(
        &self,
        initial_k: f64,
        heat_j_kg: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).temperature_after_energy_change(initial_k, heat_j_kg)
    }

    /// Exact source-Cp law for passive exchange with a fixed-temperature reservoir.
    /// Solves `integral(initial..final) Cp(T)/(reservoir_k-T) dT = exposure`,
    /// where `exposure = conductance_time_over_mass_j_kg_k` is nonnegative.
    /// For a lumped body it is `G*dt/m`; for a steady plug-flow heat exchanger
    /// it is `UA/m_dot`. G and UA have units W/K. The caller owns the validity
    /// of that transfer model and the fixed-pressure, fixed-composition path.
    ///
    /// Source breakpoints and affine Cp laws are integrated without resampling
    /// or freezing capacity. Constant-Cp intervals use the analytic exponential;
    /// the last varying-Cp interval solves relative temperature travel to
    /// floating-point resolution. Zero exposure or equal temperatures preserve the
    /// initial temperature. The reservoir may lie outside source support, but
    /// the entire material temperature path must remain inside it.
    ///
    /// Use [`Self::linear_heat_exchange`] to retain the integrated heat and
    /// exchange mean without reconstructing them from rounded endpoints.
    /// There is no source-uncertainty or
    /// heat-transfer-model error bound implied by this numerical inversion.
    ///
    /// # Errors
    /// Refuses invalid temperatures/exposure, a path beyond the original source
    /// interval, unrepresentable nonzero motion or exhausted inversion work.
    pub fn temperature_after_linear_exchange(
        &self,
        initial_k: f64,
        reservoir_k: f64,
        conductance_time_over_mass_j_kg_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).temperature_after_linear_exchange(
            initial_k,
            reservoir_k,
            conductance_time_over_mass_j_kg_k,
        )
    }

    /// The same exchange with its signed source-Cp enthalpy integral and direct
    /// conductance-weighted mean temperature. These use the retained relative
    /// travel through every original source interval, preserving weak exchange
    /// that subtracting two rounded absolute temperatures would distort.
    /// Admission and refusal semantics are identical to
    /// [`Self::temperature_after_linear_exchange`], including refusal of nonzero
    /// exchange when the absolute outlet temperature cannot move at this scale.
    pub fn linear_heat_exchange(
        &self,
        initial_k: f64,
        reservoir_k: f64,
        conductance_time_over_mass_j_kg_k: f64,
    ) -> Result<LinearHeatExchange, ThermalTransportError> {
        source_linear_heat_exchange(
            &self.knots,
            initial_k,
            reservoir_k,
            conductance_time_over_mass_j_kg_k,
        )
    }
}

/// Shared numerical owner for already admitted original specific-Cp knots.
/// Source identity and pressure/composition admission remain with each adapter.
pub(crate) fn source_linear_heat_exchange(
    knots: &[[f64; 2]],
    initial_k: f64,
    reservoir_k: f64,
    exposure_j_kg_k: f64,
) -> Result<LinearHeatExchange, ThermalTransportError> {
    let answer = CapacityCurve::new(knots).linear_heat_exchange(
        initial_k, reservoir_k, exposure_j_kg_k,
    )?;
    Ok(LinearHeatExchange {
        temperature_k: answer.temperature,
        specific_enthalpy_change_j_kg: answer.energy_change,
        mean_temperature_k: answer.mean_temperature,
        inlet_enthalpy_gain: answer.inlet_energy_gain,
        inlet_enthalpy_loss: answer.inlet_energy_loss,
        reservoir_enthalpy_derivative_j_kg_k: answer.reservoir_energy_derivative,
        mean_reservoir_derivative: answer.mean_reservoir_derivative,
    })
}

/// Exact exposure to a boundary of one already admitted original Cp interval.
/// The phase owner uses this to split sensible and latent residence without
/// copying the logarithmic exchange kernel or approximating a transition time.
pub(crate) fn source_linear_heat_exchange_exposure(
    knots: &[[f64; 2]; 2],
    initial_k: f64,
    final_k: f64,
    reservoir_k: f64,
) -> Result<f64, ThermalTransportError> {
    CapacityCurve::new(knots).exchange_integral(0, initial_k, final_k, reservoir_k)
}

impl ThermalLaw for IsobaricHeatCapacityLaw {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity.as_bytes())
    }

    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }

    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        let [low, high] = self.enthalpy_span;
        if !specific_h_j_kg.is_finite() || specific_h_j_kg < low || specific_h_j_kg > high {
            return Err(DepositionError::MaterialDomain);
        }
        let [cold, hot] = self.temperature_span_k();
        if specific_h_j_kg == low {
            return Ok(cold);
        }
        if specific_h_j_kg == high {
            return Ok(hot);
        }
        self.inverse.as_ref().map_or_else(
            || self.temperature_after_specific_heat(self.datum_temperature_k, specific_h_j_kg),
            |inverse| inverse.temperature(&CapacityCurve::new(&self.knots), specific_h_j_kg),
        )
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let initial_k = self.temperature_k(initial_h_j_kg)?;
        let final_k = self.temperature_k(final_h_j_kg)?;
        if initial_h_j_kg == final_h_j_kg {
            return Ok(Some(0.0));
        }
        let dh = final_h_j_kg - initial_h_j_kg;
        let change = if initial_k == final_k {
            dh / initial_k
        } else {
            let ds = self
                .specific_entropy_change_at_temperatures_j_kg_k(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            let represented_dh = self
                .specific_enthalpy_change_j_kg(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            // Preserve actual stored dh when inversion rounds the endpoints.
            dh * (ds / represented_dh)
        };
        if !change.is_finite()
            || change == 0.0
            || change.is_sign_positive() != (final_h_j_kg > initial_h_j_kg)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(change))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue,
        Provenance, UncertaintyModel,
    };
    use fs_qty::Temperature;

    #[test]
    fn static_pressure_quantity_is_preserved_without_accepting_stress_or_rms_aliases() {
        let typed = |kind, form| QuantitySpec::semantic(SemanticType::new(kind, form));
        for (pressure, accepted) in [
            (QuantitySpec::dimensional(Pressure::DIMS), true),
            (typed(QuantityKind::Pressure, ValueForm::Static), true),
            (typed(QuantityKind::Stress, ValueForm::Static), false),
            (typed(QuantityKind::Pressure, ValueForm::Rms), false),
        ] {
            let temperature = typed(QuantityKind::AbsoluteTemperature, ValueForm::Static);
            let key = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT);
            let mut claims = ClaimSet::new();
            claims.insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(290.0, 2.0), (300.0, 2.0), (310.0, 4.0)],
                    dims: SPECIFIC_HEAT,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature, 290.0, 310.0)
                    .with_quantity("p", pressure, 100000.0, 100000.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic isobaric capacity".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            }).unwrap();
            let card = MaterialCard::assemble(MaterialStateId {
                chemistry: "synthetic caloric fluid".into(),
                phase: "liquid".into(), process: "numerical fixture".into(), revision: 0,
            }, claims, vec![]).unwrap();
            let point = |t| QueryPoint::new()
                .with_quantity("T", temperature, t).unwrap()
                .with_quantity("p", pressure, 100000.0).unwrap();
            let lower = point(290.0);
            let upper = point(310.0);
            let law = IsobaricHeatCapacityLaw::try_new(&card, IsobaricHeatCapacityQuery {
                specific_heat_cp: &key, temperature_axis: "T", pressure_axis: "p",
                lower: &lower, upper: &upper, selection: ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
            }, 300.0);
            assert_eq!(law.is_ok(), accepted, "{pressure:?}");
            if let Ok(law) = law {
                assert_eq!(law.material_card_identity(), card.content_hash());
                assert_eq!(law.query_endpoints(), &[lower.clone(), upper.clone()]);
                assert_eq!(law.query_endpoints()[0].axis_quantities()["p"], pressure);
                assert_eq!(law.caloric_identity(), Some(*law.identity().as_bytes()));
                for target in [290.0, 293.25, 300.0, 306.75, 310.0] {
                    let h = law.specific_enthalpy_j_kg(target).unwrap();
                    assert!((law.temperature_k(h).unwrap() - target).abs()
                        <= 4.0 * f64::EPSILON * target);
                }
                assert!((law.specific_enthalpy_change_j_kg(300.0, 310.0).unwrap() - 30.0).abs() < 1e-13);
                assert!((law.specific_enthalpy_change_j_kg(305.0, 300.0).unwrap() + 12.5).abs() < 1e-13);
                for datum in [290.0_f64, 295.25, 300.0, 305.25, 310.0] {
                    let indexed = IsobaricHeatCapacityLaw::try_new(&card, IsobaricHeatCapacityQuery {
                        specific_heat_cp: &key, temperature_axis: "T", pressure_axis: "p",
                        lower: &lower, upper: &upper,
                        selection: ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
                    }, datum).unwrap();
                    let mut direct = indexed.clone();
                    direct.inverse = None;
                    for target in [290.0, 293.25, 300.0, 306.75, 310.0,
                        datum.next_down(), datum, datum.next_up()]
                    {
                        if !(290.0..=310.0).contains(&target) { continue; }
                        let expected = direct.specific_enthalpy_j_kg(target).unwrap();
                        let h = indexed.specific_enthalpy_j_kg(target).unwrap();
                        assert!((h - expected).abs() <= 8.0 * f64::EPSILON * expected.abs());
                        assert_eq!(indexed.specific_enthalpy_change_j_kg(datum, target).unwrap(), h);
                        assert_eq!(indexed.specific_enthalpy_change_j_kg(target, datum).unwrap(), -h);
                        assert!((indexed.temperature_k(h).unwrap() - target).abs()
                            <= 8.0 * f64::EPSILON * target);
                        if target == 290.0 || target == 310.0 || (target - datum).abs() < 1e-10 {
                            assert_eq!(h, expected);
                        }
                    }
                    // Unaligned relative queries retain their direct small increment.
                    let a = 306.75_f64;
                    let dh = indexed.specific_enthalpy_change_j_kg(a, a.next_up()).unwrap();
                    assert!(dh > 0.0);
                    assert_eq!(dh, direct.specific_enthalpy_change_j_kg(a, a.next_up()).unwrap());
                    for outside in [289.0, 311.0, f64::NAN, f64::INFINITY] {
                        assert!(indexed.specific_enthalpy_j_kg(outside).is_err());
                        assert!(indexed.specific_enthalpy_change_j_kg(outside, datum).is_err());
                    }
                }
            }
        }
    }
}
