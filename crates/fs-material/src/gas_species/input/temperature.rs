//! Isobaric property feedback from one retained gas composition and source law.
//!
//! Density, viscosity, conductivity and heat capacity are evaluated together.
//! This is a property law for low-Mach/Boussinesq consumers, not a compressible
//! mass or energy integrator. In particular, warming humid air preserves its
//! water mole fraction, not its relative humidity; no water enters from an
//! unspecified reservoir. Every endpoint retains the original phase admission.

use super::{GasInputError, GasRecord};
use crate::MaterialError;
use crate::gas::{GasSpec, GasState, saturation_pressure_water_pa};
use crate::gas_species::{GasSpecies, resolve_gas_mixture};
use crate::thermochemical_gas::ThermochemicalGasModel;
use fs_matdb::QueryPoint;
use fs_qty::{
    Temperature,
    semantic::{Composition, CompositionBasis},
};

#[derive(Clone, Debug, PartialEq)]
enum Source {
    MoistAir {
        composition: Composition,
    },
    Mixture {
        species: Vec<GasSpecies>,
        composition: Composition,
    },
    Thermochemical {
        model: ThermochemicalGasModel,
        point: QueryPoint,
    },
}

/// Original admitted reference and its temperature-dependent gas law.
/// Pressure, species fractions, source choices and all other conditions stay
/// fixed. Buoyancy consumers use the local density departure from `reference`.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricGasLaw {
    reference: GasRecord,
    source: Source,
}

impl IsobaricGasLaw {
    /// Bind the shared legacy moist-air or named-species input grammar.
    /// Thermochemical snapshots need their original model; use
    /// [`super::packs::read_source_temperature_law`] for canonical source packs.
    pub fn try_from_record(reference: GasRecord) -> Result<Self, GasInputError> {
        let source = match &reference {
            GasRecord::Legacy(state) => Source::MoistAir {
                composition: Composition::new(
                    CompositionBasis::MoleFraction,
                    vec![1.0 - state.water_mole_fraction, state.water_mole_fraction],
                )
                .map_err(|_| GasInputError::Input("invalid retained humid-air composition"))?,
            },
            GasRecord::Mixture(_, components) => Source::Mixture {
                species: components.iter().map(|(species, _)| *species).collect(),
                composition: Composition::new(
                    CompositionBasis::MoleFraction,
                    components.iter().map(|(_, fraction)| *fraction).collect(),
                )
                .map_err(|_| GasInputError::Input("invalid retained gas composition"))?,
            },
            GasRecord::Thermochemical(_) => {
                return Err(GasInputError::Input(
                    "temperature feedback requires the original thermochemical model, not a snapshot",
                ));
            }
        };
        reference
            .state()
            .validate_coherence()
            .map_err(GasInputError::Material)?;
        let law = Self { reference, source };
        // GasRecord variants are public. Do not turn an unrelated, numerically
        // coherent snapshot into a different gas when feedback begins.
        let reproduced = law
            .at_temperature(Temperature::new(law.reference.state().temperature))
            .map_err(GasInputError::Material)?;
        let a = law.reference.state();
        let b = reproduced.state();
        for (x, y) in [
            (a.density, b.density),
            (a.specific_heat_cp, b.specific_heat_cp),
            (a.dynamic_viscosity, b.dynamic_viscosity),
            (a.thermal_conductivity, b.thermal_conductivity),
            (a.specific_gas_constant, b.specific_gas_constant),
            (a.water_mole_fraction, b.water_mole_fraction),
        ] {
            if (x - y).abs() > 128.0 * f64::EPSILON * x.abs().max(y.abs()) {
                return Err(GasInputError::Input(
                    "retained gas does not match its declared source law",
                ));
            }
        }
        Ok(law)
    }

    /// Retain an original NASA-9/transport model and typed reference query.
    /// This path is explicitly frozen gas: it makes no condensation or
    /// reaction claim. The model enforces all source domains at every query.
    pub fn try_from_model(
        model: ThermochemicalGasModel,
        point: QueryPoint,
    ) -> Result<Self, MaterialError> {
        let state = model.evaluate(&point)?;
        let model = model.pinned_to_snapshot(&state)?;
        let reference = GasRecord::Thermochemical(state);
        Ok(Self {
            reference,
            source: Source::Thermochemical { model, point },
        })
    }

    /// Complete original reference, including source receipts where present.
    #[must_use]
    pub const fn reference(&self) -> &GasRecord {
        &self.reference
    }

    /// Admit a continuous, stable caloric interval on this retained gas law.
    ///
    /// This initialization check retains the original pressure, composition and
    /// transport/phase admission at both ends. Source-backed gases reuse the
    /// existing bounded NASA-9 coefficient analysis to establish positive Cv
    /// (and therefore Cp) throughout the interval. Raw sources must remain in
    /// one region per species; crossing seams requires an already explicit
    /// continuous caloric reference. No source is fitted or reconciled here.
    ///
    /// This supplies interval admission, not an energy integrator. Consumers
    /// still evaluate the original Cp and enthalpy change at each trial.
    pub fn admit_caloric_interval(
        &self,
        lower: Temperature,
        upper: Temperature,
    ) -> Result<(), MaterialError> {
        if !lower.value().is_finite() || !upper.value().is_finite()
            || lower.value() <= 0.0 || lower.value() >= upper.value()
        {
            return Err(MaterialError::State {
                what: "caloric interval requires ordered positive absolute temperatures".into(),
            });
        }
        self.at_temperature(lower)?;
        self.at_temperature(upper)?;
        if let Source::Thermochemical { model, .. } = &self.source {
            // Only the original coefficient-bound initialization is used.
            // The root-solver fields are validated but no root solve is run.
            let controls = crate::thermochemical_gas::energy::GasEnergySolve {
                minimum_temperature_k: lower.value(),
                maximum_temperature_k: upper.value(),
                energy_tolerance_j_kg: 1.0,
                maximum_iterations: 1,
            };
            model.maximum_isochoric_temperature_slope(controls).map_err(|error| {
                MaterialError::State {
                    what: format!("retained gas continuous caloric interval: {error}"),
                }
            })?;
        }
        Ok(())
    }

    /// Isobaric specific enthalpy change `h(to) - h(from)` [J/kg] of this
    /// same retained gas. Both endpoints retain the original pressure,
    /// composition, source claims and phase/domain admission.
    ///
    /// Thermochemical gases evaluate their original NASA-9 enthalpy changes
    /// without subtracting large absolute formation enthalpies. Crossed
    /// temperatures must have contiguous source coverage. Raw region-reference
    /// jumps remain as authored; this check neither reconciles those jumps nor
    /// certifies a continuous caloric branch. An explicitly reconciled model
    /// integrates its retained piecewise Cp chart across each seam.
    /// Legacy and named-species engineering models are
    /// calorically perfect, so their own constant mixture Cp gives the exact
    /// change within that declared model. No gas is replaced by dry air.
    ///
    /// # Errors
    /// Refuses unsupported endpoints, supersaturation, caloric source gaps or
    /// a nonzero enthalpy change that is not finite, representable, or directed
    /// consistently with the positive-Cp temperature change.
    pub fn specific_enthalpy_change(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, MaterialError> {
        if from.value() != self.reference.state().temperature {
            self.at_temperature(from)?;
        }
        let endpoint = self.at_temperature(to)?;
        let delta = to.value() - from.value();
        let change = match &self.source {
            Source::Thermochemical { model, .. } => {
                let lower = from.value().min(to.value());
                let upper = from.value().max(to.value());
                // Raw NASA-9 endpoint differences intentionally allow gaps.
                // A thermal increment needs support along the crossed path,
                // without altering the retained source-reference convention.
                for component in model.thermodynamics().components() {
                    for pair in component.regions().windows(2) {
                        let end = pair[0].temperature_max().value();
                        let start = pair[1].temperature_min().value();
                        if end < start && lower < start && upper > end {
                            return Err(MaterialError::State {
                                what: format!(
                                    "retained gas caloric source {} has an unsupported temperature gap ({end}, {start}) K",
                                    component.species().as_str()
                                ),
                            });
                        }
                    }
                }
                model
                    .thermodynamics()
                    .specific_enthalpy_change(from, to)
                    .map_err(|error| MaterialError::State {
                        what: format!("retained gas caloric source: {error}"),
                    })?
                    .value()
            }
            Source::MoistAir { .. } | Source::Mixture { .. } => {
                endpoint.state().specific_heat_cp * delta
            }
        };
        if !change.is_finite()
            || (delta != 0.0
                && (change == 0.0 || change.is_sign_positive() != delta.is_sign_positive()))
        {
            return Err(MaterialError::State {
                what: "retained gas enthalpy change must be finite, representable and follow positive Cp".into(),
            });
        }
        Ok(change)
    }

    /// Resolve all properties at one local absolute temperature and the
    /// original pressure/composition. Does not mutate the reference or model.
    /// Supersaturation in legacy/named humid gases and source-domain failures
    /// refuse; no extrapolation or air-property substitution is performed.
    pub fn at_temperature(&self, temperature: Temperature) -> Result<GasRecord, MaterialError> {
        let t = temperature.value();
        if !t.is_finite() || t <= 0.0 {
            return Err(MaterialError::State {
                what: "gas temperature must be finite and positive".into(),
            });
        }
        let reference = self.reference.state();
        match &self.source {
            Source::MoistAir { composition } => {
                let water = reference.water_mole_fraction;
                if water == 0.0 {
                    return Ok(GasRecord::Legacy(GasState::try_new(
                        &GasSpec::dry_air_ussa1976(),
                        t,
                        reference.pressure,
                    )?));
                }
                // Compare mole fractions directly. A round trip through RH
                // can put an admitted saturated reference one ULP above 1.
                if water > 0.15 || water > saturation_pressure_water_pa(t)? / reference.pressure {
                    return Err(MaterialError::Parameters {
                        what: "retained water mole fraction exceeds saturation or the dilute-mixture ceiling".into(),
                    });
                }
                Ok(GasRecord::Legacy(GasState::try_new_mixture(
                    &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()],
                    composition,
                    t,
                    reference.pressure,
                    Some(1),
                )?))
            }
            Source::Mixture {
                species,
                composition,
            } => {
                let state = resolve_gas_mixture(species, composition, t, reference.pressure)?;
                Ok(GasRecord::Mixture(
                    state,
                    species
                        .iter()
                        .copied()
                        .zip(composition.fractions().iter().copied())
                        .collect(),
                ))
            }
            Source::Thermochemical { model, point } => {
                let point =
                    point
                        .clone()
                        .with("temperature", t)
                        .map_err(|error| MaterialError::State {
                            what: error.to_string(),
                        })?;
                Ok(GasRecord::Thermochemical(model.evaluate(&point)?))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_transport::FluidTransport;
    use crate::gas_species::input::read_gas_record;

    fn law(text: &str) -> IsobaricGasLaw {
        IsobaricGasLaw::try_from_record(
            read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn local_temperature_drives_density_and_transport_at_fixed_pressure_and_composition() {
        for text in [
            "300 100000 0",
            "300 100000 0.5",
            "mixture 300 100000 mole 2 He 0.4 N2 0.6",
        ] {
            let law = law(text);
            let hot = law.at_temperature(Temperature::new(310.0)).unwrap();
            let (a, b) = (law.reference().state(), hot.state());
            assert!((b.density / a.density - 300.0 / 310.0).abs() < 2e-14);
            assert_eq!(a.pressure, b.pressure);
            assert!((a.water_mole_fraction - b.water_mole_fraction).abs() < 1e-16);
            assert!((a.specific_gas_constant - b.specific_gas_constant).abs() < 1e-12);
            assert!(b.dynamic_viscosity > a.dynamic_viscosity);
            assert!(b.thermal_conductivity > a.thermal_conductivity);
            assert_eq!(hot.transport_state().unwrap().density_kg_m3(), b.density);
        }
    }

    #[test]
    fn cooling_cannot_create_condensation_or_change_the_retained_source() {
        let humid = law("300 100000 0.9");
        let original = humid.clone();
        assert!(humid.at_temperature(Temperature::new(280.0)).is_err());
        assert!(humid.at_temperature(Temperature::new(f64::NAN)).is_err());
        assert_eq!(humid, original);
        let helium = law("mixture 300 100000 mole 1 He 1");
        assert!(helium.at_temperature(Temperature::new(400.0)).is_err());
        assert!(
            IsobaricGasLaw::try_from_record(GasRecord::Legacy(*helium.reference().state()))
                .is_err()
        );
        let GasRecord::Mixture(mut state, composition) = helium.reference().clone() else {
            unreachable!()
        };
        state.water_mole_fraction = 0.05;
        assert!(IsobaricGasLaw::try_from_record(GasRecord::Mixture(state, composition)).is_err());
    }

    #[test]
    fn saturated_reference_is_admitted_without_a_humidity_round_trip() {
        let saturated = law("300 50001 1");
        assert_eq!(
            saturated.at_temperature(Temperature::new(300.0)).unwrap(),
            *saturated.reference(),
        );
        let warm = saturated.at_temperature(Temperature::new(301.0)).unwrap();
        assert_eq!(
            warm.state().water_mole_fraction,
            saturated.reference().state().water_mole_fraction
        );
        assert!(saturated.at_temperature(Temperature::new(299.0)).is_err());
    }

    #[test]
    fn retained_gas_caloric_change_uses_its_own_composition_and_reverses() {
        for text in [
            "300 100000 0",
            "300 100000 0.5",
            "mixture 300 100000 mole 2 He 0.4 N2 0.6",
        ] {
            let law = law(text);
            let (a, b) = (Temperature::new(300.0), Temperature::new(310.0));
            law.admit_caloric_interval(a, b).unwrap();
            assert!(law.admit_caloric_interval(b, a).is_err());
            assert!(law.admit_caloric_interval(a, a).is_err());
            let change = law.specific_enthalpy_change(a, b).unwrap();
            assert_eq!(change, law.reference().state().specific_heat_cp * 10.0);
            assert_eq!(law.specific_enthalpy_change(b, a).unwrap(), -change);
            assert_eq!(law.specific_enthalpy_change(a, a).unwrap(), 0.0);
        }
        let humid = law("300 100000 0.9");
        assert!(
            humid
                .specific_enthalpy_change(Temperature::new(300.0), Temperature::new(280.0))
                .is_err()
        );
        assert!(
            humid
                .specific_enthalpy_change(Temperature::new(280.0), Temperature::new(300.0))
                .is_err()
        );
    }

    #[test]
    fn retained_caloric_path_refuses_gaps_without_rewriting_raw_region_references() {
        use crate::thermochemical_gas::tests::{point, species_in_domain, transport_in_domain};
        use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
        use fs_thermochem::{Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

        for gap in [0.0, 10.0] {
            let low = species_in_domain("A", 0.028, 3.5, 0.001, 1e18, [250.0, 300.0]);
            let high =
                species_in_domain("A", 0.028, 3.5, 0.001, 1e18 + 1024.0, [300.0 + gap, 400.0]);
            let species = Nasa9StandardStateModelV1::new(
                low.species().clone(),
                low.molar_mass(),
                low.convention().clone(),
                vec![low.regions()[0].clone(), high.regions()[0].clone()],
            )
            .unwrap();
            let model = ThermochemicalGasModel::new(
                FrozenIdealGasMixtureModelV1::new(
                    vec![species],
                    Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
                )
                .unwrap(),
                vec![transport_in_domain(
                    "A",
                    1.7e-5,
                    0.025,
                    [250.0, 400.0],
                    [50000.0, 150000.0],
                )],
                None,
            )
            .unwrap();
            let law =
                IsobaricGasLaw::try_from_model(model.clone(), point(280.0, 100000.0)).unwrap();
            for (a, b) in [(280.0, 290.0), (320.0, 330.0)] {
                let (a, b) = (Temperature::new(a), Temperature::new(b));
                law.admit_caloric_interval(a, b).unwrap();
                let change = law.specific_enthalpy_change(a, b).unwrap();
                assert_eq!(
                    change,
                    model
                        .thermodynamics()
                        .specific_enthalpy_change(a, b)
                        .unwrap()
                        .value()
                );
                assert_eq!(law.specific_enthalpy_change(b, a).unwrap(), -change);
            }
            // Stable relative arithmetic survives a formation reference that
            // would erase this small sensible increment by absolute subtraction.
            let (a, b) = (289.0, 289.0 + 1e-6);
            let change = law
                .specific_enthalpy_change(Temperature::new(a), Temperature::new(b))
                .unwrap();
            let expected = (R / 0.028) * (3.5 + 0.0005 * (a + b)) * (b - a);
            assert!((change / expected - 1.0).abs() < 1e-13);

            let (a, b) = (Temperature::new(290.0), Temperature::new(320.0));
            assert!(law.admit_caloric_interval(a, b).is_err());
            if gap == 0.0 {
                let reconciled = model.with_caloric_continuity(
                    Temperature::new(280.0),
                    fs_thermochem::MolarEnergyQuantityV1::new(1e6),
                    fs_thermochem::MolarThermalQuantityV1::new(1.0),
                ).unwrap();
                let continuous = IsobaricGasLaw::try_from_model(
                    reconciled, point(280.0, 100000.0),
                ).unwrap();
                continuous.admit_caloric_interval(a, b).unwrap();
                let expected = (R / 0.028) * (3.5 + 0.0005 * 610.0) * 30.0;
                assert!((continuous.specific_enthalpy_change(a, b).unwrap()
                    / expected - 1.0).abs() < 1e-13);
            }
            // Both original endpoints, and the raw endpoint-difference API,
            // support this query even when the intervening 300..310 K does not.
            law.at_temperature(a).unwrap();
            law.at_temperature(b).unwrap();
            let original = model
                .thermodynamics()
                .specific_enthalpy_change(a, b)
                .unwrap()
                .value();
            if gap == 0.0 {
                assert_eq!(law.specific_enthalpy_change(a, b).unwrap(), original);
                assert_eq!(law.specific_enthalpy_change(b, a).unwrap(), -original);
                let continuous_cp_increment = (R / 0.028) * (3.5 + 0.0005 * 610.0) * 30.0;
                assert!(
                    (original - continuous_cp_increment).abs() > 500.0,
                    "the original region-reference jump must not be silently removed"
                );
            } else {
                for (from, to) in [
                    (a, b),
                    (b, a),
                    (Temperature::new(300.0), Temperature::new(310.0)),
                ] {
                    let error = law.specific_enthalpy_change(from, to).unwrap_err();
                    assert!(error.to_string().contains("unsupported temperature gap"));
                }
            }
        }
    }

    #[test]
    fn thermochemical_feedback_retains_the_original_models_and_query_conditions() {
        use crate::thermochemical_gas::tests::{model_in_domains, point};
        let model = model_in_domains(
            vec![0.8, 0.2],
            0.001,
            [250.0, 400.0],
            [290.0, 330.0],
            [50000.0, 150000.0],
        );
        let law = IsobaricGasLaw::try_from_model(model.clone(), point(300.0, 100000.0)).unwrap();
        let GasRecord::Thermochemical(reference) = law.reference() else {
            unreachable!()
        };
        let expected = model
            .evaluate_temperature_from_snapshot(reference, Temperature::new(315.0))
            .unwrap();
        assert_eq!(
            law.at_temperature(Temperature::new(315.0)).unwrap(),
            GasRecord::Thermochemical(expected)
        );
        assert!(law.at_temperature(Temperature::new(340.0)).is_err());
        assert!(IsobaricGasLaw::try_from_record(law.reference().clone()).is_err());
        assert_eq!(law.reference().state().temperature, 300.0);
        let change = law
            .specific_enthalpy_change(Temperature::new(300.0), Temperature::new(315.0))
            .unwrap();
        let expected = model
            .thermodynamics()
            .specific_enthalpy_change(Temperature::new(300.0), Temperature::new(315.0))
            .unwrap()
            .value();
        assert_eq!(change, expected);
        assert_eq!(
            law.specific_enthalpy_change(Temperature::new(315.0), Temperature::new(300.0))
                .unwrap(),
            -change
        );
        assert!(
            law.specific_enthalpy_change(Temperature::new(300.0), Temperature::new(340.0))
                .is_err()
        );
    }
}
