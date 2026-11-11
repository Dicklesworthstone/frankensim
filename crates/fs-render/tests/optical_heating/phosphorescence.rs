//! G1/G3/G4: intersystem crossing, afterglow and native optical/thermal closure.

use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::fluorescence::FluorescenceSource;
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::luminescence::{
    LuminescenceSource, LuminescentOpticalHeating, SingletTripletLimits, SingletTripletSource,
};
use fs_render::volumes::transport::luminescence::LuminescenceGrid;
use fs_thermal::enthalpy::DepositionObserver;
use fs_thermal::luminescence::ArrheniusQuenching;

type Heating = LuminescentOpticalHeating<1, 2>;

fn coupled(
    heating: Heating,
    thermal: &DepositionSimulation<'_>,
    triplet: LuminescenceSource<2>,
    crossing_rates_s: [f64; 2],
) -> Heating {
    heating
        .with_singlet_triplet(
            SingletTripletSource::try_new(
                triplet,
                crossing_rates_s[0],
                ArrheniusQuenching::try_new(crossing_rates_s[1], 0.0).unwrap(),
                SingletTripletLimits::default(),
            )
            .unwrap(),
            thermal,
            || Ok(()),
        )
        .unwrap()
}

fn source(
    level_nm: f64,
    lifetime_s: f64,
    yield_: f64,
    weights: [f64; 2],
    id: u8,
) -> LuminescenceSource<2> {
    LuminescenceSource::try_new(
        FluorescenceSource::try_new(
            ContentHash([id; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            yield_,
            [650.0, 800.0],
            weights,
        )
        .unwrap(),
        level_nm,
        lifetime_s,
    )
    .unwrap()
}
fn singlet() -> LuminescenceSource<2> {
    source(500.0, 1.0, 1.0, [1.0, 0.0], 181)
}
fn triplet() -> LuminescenceSource<2> {
    source(600.0, 5.0, 0.5, [0.0, 1.0], 182)
}

fn absorption() -> MassAbsorptionSpectrum {
    MassAbsorptionSpectrum::source_table(
        &[(380.0, 1e4), (500.0, 1e4), (650.0, 0.0), (800.0, 0.0)],
        ContentHash([183; 32]),
    )
    .unwrap()
}

fn setup(cx: &Cx<'_>, singlet: LuminescenceSource<2>) -> (DepositionSimulation<'static>, Heating) {
    let density = [1.0];
    let medium = grid(cx, &density, absorption());
    let deposition = medium
        .deposit_beam_power(
            cx,
            Point3::new(0.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            &[400.0],
            &[10.0],
            1,
        )
        .unwrap();
    let thermal = thermal(vec![cell(1, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let heating = Heating::try_new(
        &medium,
        singlet,
        deposition,
        &[(0, CellId(1))],
        &thermal,
        1,
        1e-12,
        1e-8,
        || Ok(()),
    )
    .unwrap();
    (thermal, heating)
}

fn pulse(heating: Heating, end_s: f64) -> Heating {
    heating
        .with_pump_schedule(
            SpectralPowerSchedule::try_new(
                &[SpectralPowerSegment {
                    start_s: 0.0,
                    end_s,
                    start_power_w: [10.0],
                    end_power_w: [10.0],
                }],
                1,
            )
            .unwrap(),
        )
        .unwrap()
}

// Independent Bateman solution: kST=3/s, intrinsic S removal=1/s,
// T removal=0.2/s, an 8 W singlet pump and Et/Es=500/600.
fn charged_stores(time_s: f64) -> [f64; 2] {
    [
        2.0 * (1.0 - (-4.0 * time_s).exp()),
        (500.0 / 600.0)
            * 6.0
            * ((1.0 - (-0.2 * time_s).exp()) / 0.2
                - ((-0.2 * time_s).exp() - (-4.0 * time_s).exp()) / 3.8),
    ]
}
fn dark_stores(initial: [f64; 2], duration_s: f64) -> [f64; 2] {
    let fast = (-4.0 * duration_s).exp();
    let slow = (-0.2 * duration_s).exp();
    [
        initial[0] * fast,
        initial[1] * slow + (500.0 / 600.0) * 3.0 * initial[0] * (slow - fast) / 3.8,
    ]
}

fn partition(heating: &Heating) {
    let e = heating.energy();
    let heat = e.prompt_heat_j + e.nonradiative_heat_j + e.relaxation_heat_j + e.crossing_heat_j
        - e.reverse_crossing_heat_absorbed_j
        + e.secondary_absorption_heat_j;
    close(heat, heating.pump_energy().thermal_input_j);
    close(
        heating.pump_energy().absorbed_j[0],
        e.escaped_j.iter().sum::<f64>()
            + e.unresolved_radiation_j.iter().sum::<f64>()
            + heat
            + e.stored_excitation_j,
    );
    assert!(e.residual_j.abs() < 1e-8);
    assert!(e.absolute_roundoff_j < 1e-8);
}

#[test]
fn g0_native_rates_preserve_quenching_below_the_radiative_yield_ulp() {
    let quenching = ArrheniusQuenching::try_new(1e-18, 0.0).unwrap();
    let singlet = singlet().with_quenching(quenching);
    let triplet = source(600.0, 1.0, 1.0, [0.0, 1.0], 184).with_quenching(quenching);
    for emitter in [&singlet, &triplet] {
        let decay = emitter.decay_at_temperature(300.0).unwrap();
        assert_eq!(decay.quantum_yield(), 1.0);
        assert_eq!(decay.nonradiative_yield(), 1e-18);
    }
    let coupling = SingletTripletSource::try_new(
        triplet,
        0.0,
        ArrheniusQuenching::try_new(0.0, 0.0).unwrap(),
        SingletTripletLimits::default(),
    )
    .unwrap();
    let native = coupling.model_at_temperature(&singlet, 300.0).unwrap();
    assert_eq!(native.radiative_rates_s(), [1.0; 2]);
    // Reconstructing these rates as 1 - quantum_yield would silently turn
    // both positive heat channels into zero despite explicit quenching data.
    assert_eq!(native.nonradiative_rates_s(), [1e-18; 2]);
}

#[test]
fn g1_intersystem_crossing_uses_two_populations_and_the_physical_energy_gap() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, singlet());
        let heating = coupled(heating, &thermal, triplet(), [3.0, 0.0]);
        let (next, accepted, audit) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        let stores = charged_stores(1.0);
        close(accepted.energy().stored_excitation_j, stores.iter().sum());
        close(accepted.energy().prompt_heat_j, 2.0);
        let singlet_integral = (8.0 - stores[0]) / 4.0;
        let triplet_input = singlet_integral * 3.0 * 500.0 / 600.0;
        let triplet_integral = (triplet_input - stores[1]) / 0.2;
        close(
            accepted.energy().escaped_j[0],
            singlet_integral * 500.0 / 650.0,
        );
        close(
            accepted.energy().escaped_j[1],
            triplet_integral * 0.1 * 600.0 / 800.0,
        );
        close(
            accepted.energy().crossing_heat_j,
            singlet_integral * 3.0 / 6.0,
        );
        close(
            accepted.energy().nonradiative_heat_j,
            triplet_integral * 0.1,
        );
        close(
            audit.external_heat_j,
            accepted.pump_energy().thermal_input_j,
        );
        close(next.energy_j(), accepted.pump_energy().thermal_input_j);
        assert_eq!(audit.substeps, 1);
        partition(&accepted);
    });
}

#[test]
fn g1_pump_off_phosphorescence_reaches_the_existing_volume_radiance_path() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, singlet());
        let heating = pulse(coupled(heating, &thermal, triplet(), [3.0, 0.0]), 1.0);
        let (_, accepted, audit) = thermal
            .advance_observed(&heating, 10.0, config(10.0), || Ok(()))
            .unwrap();
        let stores = dark_stores(charged_stores(1.0), 9.0);
        assert_eq!(audit.substeps, 2);
        close(accepted.pump_energy().absorbed_j[0], 10.0);
        close(accepted.energy().stored_excitation_j, stores.iter().sum());
        let emitted = accepted.cell_emission_w()[0].1;
        close(emitted[0], stores[0] * 500.0 / 650.0);
        close(emitted[1], stores[1] * 0.1 * 600.0 / 800.0);
        assert!(emitted[1] > 1e12 * emitted[0]);
        assert!(accepted.last_mean_cell_emission_w()[0].1[1] > emitted[1]);
        let density = [1.0];
        let medium = grid(cx, &density, absorption());
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            accepted.emission_source_identity(),
            [650.0, 800.0],
            accepted.cell_emission_w(),
            1,
        )
        .unwrap();
        let ray = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                3.0,
            )
            .unwrap();
        close(
            ray.emitted_radiance_w_m2_sr[1],
            emitted[1] / (4.0 * core::f64::consts::PI),
        );
        partition(&accepted);
    });
}

#[test]
fn g3_split_steps_and_restart_preserve_afterglow_and_thermal_energy() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, singlet());
        let heating = pulse(coupled(heating, &thermal, triplet(), [3.0, 0.0]), 0.75);
        let (whole, one, _) = thermal
            .advance_observed(&heating, 3.0, config(3.0), || Ok(()))
            .unwrap();
        let (fine, many, _) = thermal
            .advance_observed(&heating, 3.0, config(0.125), || Ok(()))
            .unwrap();
        let (middle, saved, _) = thermal
            .advance_observed(&heating, 0.5, config(0.5), || Ok(()))
            .unwrap();
        let (resumed, replay, _) = middle
            .advance_observed(&saved, 3.0, config(3.0), || Ok(()))
            .unwrap();
        for other in [&many, &replay] {
            close(
                other.energy().stored_excitation_j,
                one.energy().stored_excitation_j,
            );
            close(
                other.pump_energy().thermal_input_j,
                one.pump_energy().thermal_input_j,
            );
            for (a, b) in other.energy().escaped_j.iter().zip(one.energy().escaped_j) {
                close(*a, b);
            }
            partition(other);
        }
        close(whole.energy_j(), fine.energy_j());
        close(whole.energy_j(), resumed.energy_j());
        let (_, duplicate, _) = middle
            .advance_observed(&saved, 3.0, config(3.0), || Ok(()))
            .unwrap();
        assert_eq!(duplicate.energy(), replay.energy());
        assert_eq!(duplicate.cell_emission_w(), replay.cell_emission_w());
    });
}

#[test]
fn g1_reverse_crossing_draws_its_energy_from_the_native_thermal_bath() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, singlet());
        let heating = pulse(coupled(heating, &thermal, triplet(), [3.0, 2.0]), 0.5);
        let (next, accepted, _) = thermal
            .advance_observed(&heating, 2.0, config(2.0), || Ok(()))
            .unwrap();
        assert!(accepted.energy().crossing_heat_j > 0.0);
        assert!(accepted.energy().reverse_crossing_heat_absorbed_j > 0.0);
        assert!(accepted.cell_emission_w()[0].1[0] > 0.0);
        close(next.energy_j(), accepted.pump_energy().thermal_input_j);
        partition(&accepted);
    });
}

#[test]
fn g4_implicit_retry_discards_both_tentative_populations_and_all_channels() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, singlet());
        let heating = coupled(heating, &thermal, triplet(), [3.0, 0.5]);
        let mut candidate = heating.clone();
        candidate
            .step_sources(&thermal, 0.0, 1.0, 1.0, &mut || Ok(()))
            .unwrap();
        let retried = candidate
            .step_sources(&thermal, 0.0, 0.25, 0.25, &mut || Ok(()))
            .unwrap();
        let mut direct = heating.clone();
        assert_eq!(
            retried,
            direct
                .step_sources(&thermal, 0.0, 0.25, 0.25, &mut || Ok(()))
                .unwrap()
        );
        let (next, accepted, _) = thermal
            .advance_observed(&heating, 0.25, config(0.25), || Ok(()))
            .unwrap();
        candidate
            .step(
                0.0,
                0.25,
                &std::collections::BTreeMap::from([(CellId(1), 300.0)]),
                &next,
                &mut || Ok(()),
            )
            .unwrap();
        assert_eq!(candidate.energy(), accepted.energy());
        assert_eq!(candidate.cell_emission_w(), accepted.cell_emission_w());
        assert_eq!(heating.energy().stored_excitation_j, 0.0);
    });
}

#[test]
fn g4_cancelled_phosphorescence_does_not_publish_heat_or_excitation() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, singlet());
        let heating = coupled(heating, &thermal, triplet(), [3.0, 0.0]);
        let mut calls = 0;
        assert!(
            thermal
                .advance_observed(&heating, 1.0, config(0.01), || {
                    calls += 1;
                    if calls == 150 {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert_eq!(calls, 150);
        assert_eq!(heating.energy().stored_excitation_j, 0.0);
        assert_eq!(heating.energy().escaped_j, [0.0; 2]);
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(thermal.energy_j(), 0.0);
    });
}

#[test]
fn g0_triplet_zero_weight_blue_node_does_not_imply_an_anti_stokes_photon() {
    let spectrum = |weights| {
        FluorescenceSource::try_new(
            ContentHash([185; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            0.5,
            [550.0, 800.0],
            weights,
        )
        .unwrap()
    };
    let source = LuminescenceSource::try_new(spectrum([0.0, 1.0]), 600.0, 5.0).unwrap();
    assert_eq!(source.fluorescence().photon_fractions(), &[0.0, 1.0]);
    assert!(LuminescenceSource::try_new(spectrum([1.0, 1.0]), 600.0, 5.0).is_err());
}

#[test]
fn g1_no_crossing_recovers_the_original_single_state_history() {
    with_cx(|cx| {
        let (thermal, heating) = setup(cx, singlet());
        let scalar = pulse(heating.clone(), 0.5);
        let paired = pulse(coupled(heating, &thermal, triplet(), [0.0, 0.0]), 0.5);
        let (_, scalar, _) = thermal
            .advance_observed(&scalar, 2.0, config(2.0), || Ok(()))
            .unwrap();
        let (_, paired, _) = thermal
            .advance_observed(&paired, 2.0, config(2.0), || Ok(()))
            .unwrap();
        close(
            scalar.energy().stored_excitation_j,
            paired.energy().stored_excitation_j,
        );
        close(
            scalar.pump_energy().thermal_input_j,
            paired.pump_energy().thermal_input_j,
        );
        close(scalar.energy().escaped_j[0], paired.energy().escaped_j[0]);
        assert_eq!(paired.energy().escaped_j[1], 0.0);
        assert_eq!(paired.energy().crossing_heat_j, 0.0);
        assert_eq!(paired.energy().reverse_crossing_heat_absorbed_j, 0.0);
        partition(&paired);
    });
}

#[cfg(feature = "material-state")]
mod material_rates {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::phosphorescence::*;
    use fs_material::state_point::MaterialPropertySelection;
    use fs_qty::{
        QuantitySpec, Temperature,
        semantic::{QuantityKind, SemanticType, ValueForm},
    };
    use std::sync::Arc;

    const HC_J_NM: f64 = 6.626_070_15e-34 * 299_792_458.0 * 1e9;

    include!("saturation_material.rs");

    fn claims() -> Vec<PropertyClaim> {
        [
            (
                PHOSPHORESCENCE_ENERGY_PROPERTIES,
                PHOSPHORESCENCE_ENERGY_DIMS,
                [HC_J_NM / 500.0, HC_J_NM / 600.0],
            ),
            (
                PHOSPHORESCENCE_RADIATIVE_RATE_PROPERTIES,
                PHOSPHORESCENCE_RATE_DIMS,
                [1.0, 0.1],
            ),
            (
                PHOSPHORESCENCE_NONRADIATIVE_RATE_PROPERTIES,
                PHOSPHORESCENCE_RATE_DIMS,
                [0.0, 0.1],
            ),
            (
                PHOSPHORESCENCE_CROSSING_RATE_PROPERTIES,
                PHOSPHORESCENCE_RATE_DIMS,
                [3.0, 0.0],
            ),
        ]
        .into_iter()
        .flat_map(|(names, dims, values)| {
            names
                .into_iter()
                .zip(values)
                .map(move |(name, value)| PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained(),
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: format!("manufactured emitting material: {name}"),
                        license: "internal-test-use".into(),
                        artifact: None,
                    },
                })
        })
        .collect()
    }

    fn card(claims: Vec<PropertyClaim>) -> Arc<MaterialCard> {
        let mut set = ClaimSet::new();
        for claim in claims {
            set.insert_claim(claim).unwrap();
        }
        Arc::new(
            MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "synthetic coupled phosphorescence".into(),
                    phase: "solid".into(),
                    process: "manufactured integration".into(),
                    revision: 0,
                },
                set,
                Vec::new(),
            )
            .unwrap(),
        )
    }

    /// The absolute temperature the native thermal owner queries. A claim's
    /// temperature validity declares the same quantity.
    fn kelvin() -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ))
    }

    fn point() -> QueryPoint {
        QueryPoint::new()
            .with_quantity("T", kelvin(), 300.0)
            .unwrap()
    }

    fn native_sources(
        card: Arc<MaterialCard>,
        point: &QueryPoint,
    ) -> Result<
        (LuminescenceSource<2>, SingletTripletSource<2>),
        fs_render::volumes::heating::OpticalHeatingError,
    > {
        let material = PhosphorescenceStatePoint::resolve(
            &card,
            point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        // Deliberately inconsistent scalar yields prove the original native
        // radiative/dark rates, not a spectrum's legacy yield, own the kinetics.
        let singlet = FluorescenceSource::try_new(
            ContentHash([201; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            0.0,
            [650.0, 800.0],
            [1.0, 0.0],
        )
        .unwrap();
        let triplet = FluorescenceSource::try_new(
            ContentHash([202; 32]),
            (380.0, 500.0),
            (290.0, 800.0),
            1.0,
            [650.0, 800.0],
            [0.0, 1.0],
        )
        .unwrap();
        SingletTripletSource::try_material_temperature(
            card,
            material,
            "T",
            singlet,
            triplet,
            SingletTripletLimits::default(),
        )
    }

    fn bind(
        heating: Heating,
        coupling: SingletTripletSource<2>,
        thermal: &DepositionSimulation<'_>,
    ) -> Heating {
        heating
            .with_singlet_triplet(coupling, thermal, || Ok(()))
            .unwrap()
    }

    fn charged_for(beta: f64, t: f64) -> [f64; 2] {
        [
            2.0 * (1.0 - (-4.0 * t).exp()),
            (500.0 / 600.0)
                * 6.0
                * ((1.0 - (-beta * t).exp()) / beta
                    - ((-beta * t).exp() - (-4.0 * t).exp()) / (4.0 - beta)),
        ]
    }

    fn dark_for(initial: [f64; 2], beta: f64, t: f64) -> [f64; 2] {
        let fast = (-4.0 * t).exp();
        let slow = (-beta * t).exp();
        [
            initial[0] * fast,
            initial[1] * slow + (500.0 / 600.0) * 3.0 * initial[0] * (slow - fast) / (4.0 - beta),
        ]
    }

    #[test]
    fn g1_source_backed_stable_triplet_stores_energy_without_a_fake_lifetime() {
        with_cx(|cx| {
            let mut data = claims();
            for index in [3, 5] {
                data[index].value = PropertyValue::Scalar {
                    value: 0.0,
                    dims: PHOSPHORESCENCE_RATE_DIMS,
                };
            }
            let (singlet, coupling) = native_sources(card(data), &point()).unwrap();
            assert!(singlet.decay_at_temperature(300.0).is_err());
            assert!(coupling.triplet().decay_at_temperature(300.0).is_err());
            let (thermal, heating) = setup(cx, singlet);
            let heating = pulse(bind(heating, coupling, &thermal), 1.0);
            let (_, accepted, _) = thermal
                .advance_observed(&heating, 20.0, config(20.0), || Ok(()))
                .unwrap();
            let stores = accepted.state_stored_energy_j(0).unwrap();
            close(stores[0], 0.0);
            // 8 J into S, 3/4 of S transitions populate T, Et/Es=5/6.
            close(stores[1], 5.0);
            close(accepted.cell_emission_w()[0].1[1], 0.0);
            close(accepted.energy().emitted_j[1], 0.0);
            assert!(accepted.energy().emitted_j[0] > 0.0);
            partition(&accepted);
        });
    }

    #[test]
    fn g1_material_rates_refresh_at_the_actual_thermal_temperature() {
        with_cx(|cx| {
            let mut data = claims();
            data[3].value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.1), (400.0, 1.1)],
                dims: PHOSPHORESCENCE_RATE_DIMS,
            };
            data[3].validity =
                ValidityDomain::unconstrained().with_quantity("T", kelvin(), 300.0, 400.0);
            data[3].interpolation = InterpolationPolicy::LinearInside;
            let (singlet, coupling) = native_sources(card(data), &point()).unwrap();
            let (thermal, heating) = setup(cx, singlet);
            let heating = pulse(bind(heating, coupling, &thermal), 1.0);
            let (warm, charged, _) = thermal
                .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
                .unwrap();
            let expected = charged_for(0.2, 1.0);
            let stores = charged.state_stored_energy_j(0).unwrap();
            close(stores[0], expected[0]);
            close(stores[1], expected[1]);
            let temperature = warm.temperature_k(CellId(1)).unwrap();
            assert!(temperature > 300.0);
            let kr = 0.1 + 0.01 * (temperature - 300.0);
            close(
                charged.cell_emission_w()[0].1[1],
                stores[1] * kr * 600.0 / 800.0,
            );
            let (_, dark, _) = warm
                .advance_observed(&charged, 2.0, config(1.0), || Ok(()))
                .unwrap();
            let expected = dark_for(stores, kr + 0.1, 1.0);
            let actual = dark.state_stored_energy_j(0).unwrap();
            close(actual[0], expected[0]);
            close(actual[1], expected[1]);
            partition(&charged);
            partition(&dark);
        });
    }

    #[test]
    fn g1_sourced_quencher_conditions_change_afterglow_and_native_heat() {
        with_cx(|cx| {
            let oxygen = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AmountConcentration,
                ValueForm::Static,
            ));
            let mut data = claims();
            data[5].value = PropertyValue::Curve {
                abscissa: "dissolved_O2".into(),
                abscissa_dims: oxygen.dims(),
                knots: vec![(0.0, 0.1), (4.0, 2.1)],
                dims: PHOSPHORESCENCE_RATE_DIMS,
            };
            data[5].validity =
                ValidityDomain::unconstrained().with_quantity("dissolved_O2", oxygen, 0.0, 4.0);
            data[5].interpolation = InterpolationPolicy::LinearInside;
            let card = card(data);
            let mut results = Vec::new();
            for concentration in [0.0, 4.0] {
                let point = point()
                    .with_quantity("dissolved_O2", oxygen, concentration)
                    .unwrap();
                let (singlet, coupling) = native_sources(card.clone(), &point).unwrap();
                let (thermal, heating) = setup(cx, singlet);
                let heating = pulse(bind(heating, coupling, &thermal), 1.0);
                let (_, accepted, _) = thermal
                    .advance_observed(&heating, 5.0, config(1.0), || Ok(()))
                    .unwrap();
                let beta = 0.2 + 0.5 * concentration;
                let expected = dark_for(charged_for(beta, 1.0), beta, 4.0);
                let actual = accepted.state_stored_energy_j(0).unwrap();
                close(actual[0], expected[0]);
                close(actual[1], expected[1]);
                close(
                    accepted.cell_emission_w()[0].1[1],
                    expected[1] * 0.1 * 600.0 / 800.0,
                );
                partition(&accepted);
                results.push(accepted);
            }
            assert!(
                results[1].energy().stored_excitation_j < results[0].energy().stored_excitation_j
            );
            assert!(
                results[1].energy().nonradiative_heat_j > results[0].energy().nonradiative_heat_j
            );
            assert!(results[1].energy().emitted_j[1] < results[0].energy().emitted_j[1]);
            assert_ne!(
                results[0].emission_source_identity(),
                results[1].emission_source_identity()
            );
        });
    }

    #[test]
    fn g0_material_binding_refuses_unaccounted_level_work_and_wrong_singlet() {
        let mut data = claims();
        data[0].value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, HC_J_NM / 500.0), (400.0, 0.9 * HC_J_NM / 500.0)],
            dims: PHOSPHORESCENCE_ENERGY_DIMS,
        };
        data[0].validity =
            ValidityDomain::unconstrained().with_quantity("T", kelvin(), 300.0, 400.0);
        data[0].interpolation = InterpolationPolicy::LinearInside;
        assert!(native_sources(card(data), &point()).is_err());

        // A second temperature axis cannot silently remain frozen while the
        // native thermal owner advances an unrelated coordinate named T.
        let mut data = claims();
        data[3].value = PropertyValue::Curve {
            abscissa: "T_material".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 0.1), (400.0, 1.1)],
            dims: PHOSPHORESCENCE_RATE_DIMS,
        };
        data[3].validity =
            ValidityDomain::unconstrained().with_quantity("T_material", kelvin(), 300.0, 400.0);
        data[3].interpolation = InterpolationPolicy::LinearInside;
        let two_temperatures = point()
            .with_quantity("T_material", kelvin(), 300.0)
            .unwrap();
        assert!(native_sources(card(data), &two_temperatures).is_err());
        with_cx(|cx| {
            let (_, coupling) = native_sources(card(claims()), &point()).unwrap();
            let (thermal, wrong) = setup(cx, singlet());
            assert!(
                wrong
                    .with_singlet_triplet(coupling, &thermal, || Ok(()))
                    .is_err()
            );
        });
    }
}

/// A dark branch below one radiative-yield ulp must still heat the native cell.
#[test]
fn g1_sub_ulp_quenching_reaches_the_phosphorescence_heat_account() {
    with_cx(|cx| {
        let source = singlet().with_quenching(ArrheniusQuenching::try_new(1e-18, 0.0).unwrap());
        assert_eq!(
            source.decay_at_temperature(300.0).unwrap().quantum_yield(),
            1.0
        );
        let (thermal, heating) = setup(cx, source);
        let heating = coupled(heating, &thermal, triplet(), [0.0; 2]);
        let (_, accepted, _) = thermal
            .advance_observed(&heating, 1.0, config(1.0), || Ok(()))
            .unwrap();
        // The original 400-nm 10-W beam feeds an 8-W, 500-nm store.
        // At tau=1 and t=1 its integrated release is 8/e joules.
        let expected_heat = 8.0 * (-1.0_f64).exp() * 1e-18;
        let heat = accepted.energy().nonradiative_heat_j;
        assert!(heat > 0.0);
        assert!((heat / expected_heat - 1.0).abs() < 1e-10);
    });
}
