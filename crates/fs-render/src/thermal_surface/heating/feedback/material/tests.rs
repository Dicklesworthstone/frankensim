use super::*;
use crate::conductor::{
    CoatingIorSample, CoherentAbsorbingCoating, ConductorDataStatus, ConductorError,
    ConductorSource,
};
use crate::thermal_surface::CoatingExchangeQuadrature;
use crate::thermal_surface::heating::coating_exchange::{
    CoatingExchangeBath, CoatingExchangeSource, IsothermalCoatingStack, NativeCoatingExchange,
    NativeCoatingExchangeBudget,
};
use crate::thermal_surface::heating::{
    CoatingOpticsModel, NativeCoatingEmission, StatefulCoatedSurfaceIlluminationSource,
    StatefulCoatingLayer, StatefulSurfaceIlluminationSource, SubstrateOpticsModel,
    SurfaceHeatingBudget, SurfaceHeatingReceiver, SurfaceHeatingTolerance, SurfaceIncidence,
    SurfaceOpticalHeating,
};
use crate::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECTRAL_CONDUCTOR_ETA_PROPERTIES, SPECTRAL_CONDUCTOR_K_PROPERTIES,
    SPECTRAL_CONDUCTOR_K_QUANTITY, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use std::collections::BTreeSet;

const ID: ContentHash = ContentHash([192; 32]);
const CELL: CellId = CellId(1);
const MATERIAL: MaterialId = MaterialId(7);
const WAVELENGTH: f64 = 10_000.0;

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("thermal", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
}

fn claim(name: &str, q: QuantitySpec, values: [f64; 3]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, q),
        value: PropertyValue::Curve {
            abscissa: "thermal".into(),
            abscissa_dims: Temperature::DIMS,
            knots: [300.0, 400.0, 600.0].into_iter().zip(values).collect(),
            dims: q.dims(),
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "thermal",
                quantity(QuantityKind::AbsoluteTemperature),
                300.0,
                600.0,
            )
            .with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                90_000.0,
                110_000.0,
            ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("synthetic native photothermal claim: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}

fn fixture(
    change: impl FnOnce(&mut Vec<PropertyClaim>),
) -> (Arc<MaterialCard>, SpectralConductorStatePoint) {
    let mut claims = Vec::new();
    for i in 0..9 {
        claims.extend([
            claim(
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                [(8000.0 + 500.0 * i as f64) * 1e-9; 3],
            ),
            claim(
                SPECTRAL_CONDUCTOR_ETA_PROPERTIES[i],
                QuantitySpec::dimensional(Dims::NONE),
                [1.5, 4.0, 2.0],
            ),
            claim(
                SPECTRAL_CONDUCTOR_K_PROPERTIES[i],
                SPECTRAL_CONDUCTOR_K_QUANTITY,
                [3.0, 1.0, 2.0],
            ),
        ]);
    }
    claims.push(claim(
        "density",
        QuantitySpec::dimensional(Density::DIMS),
        [2000.0, 1900.0, 1800.0],
    ));
    change(&mut claims);
    let mut names = BTreeSet::new();
    let requirements: Vec<_> = claims
        .iter()
        .filter(|source| names.insert(source.key.name().to_owned()))
        .map(|source| {
            ScalarPropertyRequirement::try_with_key(&source.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect();
    let mut set = ClaimSet::new();
    for source in claims {
        set.insert_claim(source).unwrap();
    }
    let card = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic source-bound opaque material".into(),
                phase: "solid".into(),
                process: "photothermal regression".into(),
                revision: 0,
            },
            set,
            Vec::new(),
        )
        .unwrap(),
    );
    let resolved = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = SpectralConductorStatePoint::try_from_material(resolved).unwrap();
    (card, source)
}

fn law(card: Arc<MaterialCard>, source: SpectralConductorStatePoint) -> SurfaceConductorLaw {
    SurfaceConductorLaw::try_material_temperature(
        card,
        source,
        "thermal",
        (300.0, 600.0),
        501.0,
        28,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-5),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn normal_absorptance(temperature: f64) -> f64 {
    let (n, k) = if temperature <= 400.0 {
        let t = (temperature - 300.0) / 100.0;
        (1.5 + 2.5 * t, 3.0 - 2.0 * t)
    } else {
        let t = (temperature - 400.0) / 200.0;
        (4.0 - 2.0 * t, 1.0 + t)
    };
    4.0 * n / ((n + 1.0).powi(2) + k * k)
}

fn planck(t: f64, wavelength_nm: f64) -> f64 {
    let lambda = wavelength_nm * 1e-9;
    let h = 6.62607015e-34;
    let c = 299_792_458.0;
    let k = 1.380649e-23;
    2.0 * h * c * c / lambda.powi(5) / (h * c / (lambda * k * t)).exp_m1() * 1e-9
}

#[test]
fn original_claims_and_all_properties_refresh_without_a_new_temperature_fit() {
    let (card, original) = fixture(|claims| {
        let mut replacement = claim(
            SPECTRAL_CONDUCTOR_ETA_PROPERTIES[0],
            QuantitySpec::dimensional(Dims::NONE),
            [50.0; 3],
        );
        replacement.validity = replacement.validity.with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            400.0,
            600.0,
        );
        claims.push(replacement);
    });
    let law = law(Arc::clone(&card), original.clone());
    assert_eq!(law.row_count(), 0);
    assert_eq!(law.budget_size(), 28);
    assert_eq!(law.material_source(), Some(&original));
    let Model::Material(retained) = &law.model else {
        panic!("original material model")
    };
    let refreshed = retained.source_at_temperature(450.0).unwrap();
    assert_eq!(
        refreshed,
        original.at_query_point(&card, &point(450.0)).unwrap()
    );
    assert_eq!(refreshed.resolved().properties().len(), 28);
    for property in original.resolved().properties() {
        let after = refreshed
            .resolved()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            after.answer().receipt.selected,
            property.answer().receipt.selected
        );
    }
    near(
        refreshed.resolved().property("density").unwrap().value_si(),
        1875.0,
    );
    assert_eq!(refreshed.resolved().query().axes()["pressure"], 100_000.0);
    for t in [300.0, 350.0, 400.0, 450.0, 600.0] {
        let optics = law.at_temperature(t).unwrap();
        near(
            optics.power_balance(WAVELENGTH, 1.0, 1.0).unwrap().absorbed[0],
            normal_absorptance(t),
        );
        assert_eq!(
            optics.source().status(),
            ConductorDataStatus::MaterialStateClaims
        );
        assert!(optics.complex_ior(550.0).is_err());
    }
    near(
        law.at_temperature(350.0)
            .unwrap()
            .complex_ior(WAVELENGTH)
            .unwrap()
            .eta(),
        2.75,
    );
    assert!(law.at_temperature(299.0).is_err());
    assert!(law.at_temperature(601.0).is_err());
    assert!(law.at_specific_enthalpy(400.0, 400.0).is_err());
    assert!(
        SurfaceConductorLaw::try_material_temperature(
            Arc::clone(&card),
            original.clone(),
            "thermal",
            (300.0, 600.0),
            501.0,
            27
        )
        .is_err()
    );
    assert!(
        SurfaceConductorLaw::try_material_temperature(
            card,
            original,
            "pressure",
            (300.0, 600.0),
            501.0,
            28
        )
        .is_err()
    );
}

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (200.0..=800.0).contains(&h) {
            Ok(h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static SENSIBLE: Sensible = Sensible;

fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}

fn native(t: f64, external_w: f64) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 1,
                body_id: 1,
                material_id: MATERIAL,
                mass_kg: 1.0,
                specific_h_j_kg: t,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: external_w,
                },
            }),
        }],
        &[(MATERIAL, &SENSIBLE)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(0.0, config(0.1), || Ok(()))
    .unwrap()
    .0
}

fn receiver() -> SurfaceHeatingReceiver {
    SurfaceHeatingReceiver {
        thermal_cell: CELL,
        material_id: MATERIAL,
        temperature_bounds_k: (300.0, 600.0),
    }
}

fn illumination(
    native: &DepositionSimulation<'_>,
    optics: Arc<SurfaceConductorLaw>,
    power_w: f64,
) -> SurfaceOpticalHeating<1> {
    SurfaceOpticalHeating::try_new_with_stateful(
        &[],
        &[StatefulSurfaceIlluminationSource {
            source_id: 1,
            thermal_cell: CELL,
            material_id: MATERIAL,
            temperature_bounds_k: (300.0, 600.0),
            incidence: SurfaceIncidence {
                wavelengths_nm: [WAVELENGTH],
                incident_eta: [1.0],
                incidence_cosine: 1.0,
                s_power_fraction: [0.5],
            },
            optics,
            power: SpectralPowerSchedule::try_new(
                &[SpectralPowerSegment {
                    start_s: 0.0,
                    end_s: 1.0,
                    start_power_w: [power_w],
                    end_power_w: [power_w],
                }],
                1,
            )
            .unwrap(),
        }],
        native,
        SurfaceHeatingBudget {
            max_sources: 1,
            max_total_power_segments: 1,
            max_source_evaluations: 128,
        },
        28,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn native_absorption_heats_once_and_emission_reads_the_same_updated_source() {
    let (card, original) = fixture(|_| {});
    let law = Arc::new(law(card, original));
    let cold = native(300.0, 0.0);
    let light = illumination(&cold, Arc::clone(&law), 100.0);
    let (first, light, _) = cold
        .advance_observed(&light, 0.1, config(0.1), || Ok(()))
        .unwrap();
    let t1 = 300.0 + 10.0 * normal_absorptance(300.0);
    near(first.temperature_k(CELL).unwrap(), t1);
    let (second, light, audit) = first
        .advance_observed(&light, 0.2, config(0.1), || Ok(()))
        .unwrap();
    let t2 = t1 + 10.0 * normal_absorptance(t1);
    near(second.temperature_k(CELL).unwrap(), t2);
    near(light.energy().absorbed_j[0], t2 - 300.0);
    near(
        light.energy().thermal_input_j,
        second.energy_j() - cold.energy_j(),
    );
    near(audit.external_heat_j, t2 - t1);
    assert!(light.response_evaluations() >= 1);
    let emitter = NativeCoatingEmission::try_opaque(
        &second,
        &[],
        receiver(),
        &SubstrateOpticsModel::Law(law),
        ID,
        28,
        || Ok(()),
    )
    .unwrap();
    near(
        emitter
            .emitter()
            .radiance_w_m2_sr_nm(WAVELENGTH, 1.0)
            .unwrap(),
        normal_absorptance(t2) * planck(t2, WAVELENGTH),
    );
    assert!(emitter.emitter().radiance_w_m2_sr_nm(550.0, 1.0).is_err());
}

#[test]
fn native_infrared_exchange_uses_source_feedback_for_heating_cooling_and_equilibrium() {
    let (card, original) = fixture(|_| {});
    let law = Arc::new(law(card, original));
    for bath_temperature in [300.0, 400.0, 500.0] {
        let initial = native(400.0, 0.0);
        let physical = CoatingExchangeSource {
            source_id: 1,
            source_identity: ID,
            receiver: receiver(),
            area_m2: 0.01,
            stack: IsothermalCoatingStack::Opaque {
                films: vec![],
                substrate: SubstrateOpticsModel::Law(Arc::clone(&law)),
            },
            bath: CoatingExchangeBath::Blackbody {
                temperature_k: bath_temperature,
            },
            quadrature: CoatingExchangeQuadrature {
                wavelength_bins: 8,
                cosine_bins: 4,
            },
        };
        let budget = NativeCoatingExchangeBudget {
            max_sources: 1,
            max_total_optical_knots: 28,
            max_source_evaluations: 128,
            max_optical_evaluations: 8192,
        };
        let source =
            NativeCoatingExchange::try_new(&initial, &[physical.clone()], budget, || Ok(()))
                .unwrap();
        assert!(
            NativeCoatingExchange::try_new(
                &initial,
                &[physical],
                NativeCoatingExchangeBudget {
                    max_total_optical_knots: 27,
                    ..budget
                },
                || Ok(())
            )
            .is_err()
        );
        let (changed, accepted, audit) = initial
            .advance_observed(&source, 0.1, config(0.05), || Ok(()))
            .unwrap();
        let delta = changed.temperature_k(CELL).unwrap() - 400.0;
        assert_eq!(delta.signum(), (bath_temperature - 400.0).signum());
        near(
            changed.energy_j() - initial.energy_j(),
            -accepted.energy().outward_j,
        );
        near(audit.external_heat_j, accepted.energy().thermal_external_j);
        if bath_temperature == 400.0 {
            assert_eq!(delta, 0.0);
            assert_eq!(accepted.energy().outward_j, 0.0);
        } else {
            let t = changed.temperature_k(CELL).unwrap();
            near(
                accepted
                    .emitter(1)
                    .unwrap()
                    .radiance_w_m2_sr_nm(WAVELENGTH, 1.0)
                    .unwrap(),
                normal_absorptance(t) * planck(t, WAVELENGTH),
            );
        }
    }
}

#[test]
fn exact_transparency_is_refused_by_the_opaque_owner_initially_and_after_heating() {
    let (card, original) = fixture(|claims| {
        let target = claims
            .iter_mut()
            .find(|c| c.key.name() == SPECTRAL_CONDUCTOR_K_PROPERTIES[0])
            .unwrap();
        if let PropertyValue::Curve { knots, .. } = &mut target.value {
            for (_, value) in knots {
                *value = 0.0;
            }
        }
    });
    assert!(matches!(
        SurfaceConductorLaw::try_material_temperature(
            card,
            original,
            "thermal",
            (300.0, 600.0),
            501.0,
            28
        ),
        Err(SurfaceHeatingError::MaterialOptics(_))
    ));
    let (card, original) = fixture(|claims| {
        let target = claims
            .iter_mut()
            .find(|c| c.key.name() == SPECTRAL_CONDUCTOR_K_PROPERTIES[0])
            .unwrap();
        if let PropertyValue::Curve { knots, .. } = &mut target.value {
            knots[1].1 = 0.0;
        }
    });
    let law = Arc::new(law(card, original));
    assert!(matches!(
        law.at_temperature(400.0),
        Err(SurfaceHeatingError::MaterialOptics(_))
    ));
    let initial = native(300.0, 100.0);
    let dark = illumination(&initial, law, 0.0);
    assert!(
        initial
            .advance_observed(&dark, 1.0, config(1.0), || Ok(()))
            .is_err()
    );
    assert_eq!(initial.temperature_k(CELL).unwrap(), 300.0);
    assert_eq!(dark.time_s(), 0.0);
    assert_eq!(dark.energy().thermal_input_j, 0.0);
}

#[test]
fn source_expiry_or_moving_band_in_darkness_refuses_the_native_candidate() {
    let (card, original) = fixture(|claims| {
        let density = claims
            .iter_mut()
            .find(|c| c.key.name() == "density")
            .unwrap();
        density.validity = density.validity.clone().with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
            325.0,
        );
        let mut replacement = claim(
            "density",
            QuantitySpec::dimensional(Density::DIMS),
            [1600.0; 3],
        );
        replacement.validity = replacement.validity.with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            320.0,
            600.0,
        );
        claims.push(replacement);
    });
    let law = Arc::new(law(card, original));
    let initial = native(300.0, 100.0);
    let dark = illumination(&initial, Arc::clone(&law), 0.0);
    assert!(matches!(
        law.at_temperature(350.0),
        Err(SurfaceHeatingError::MaterialSource(_))
    ));
    assert!(
        initial
            .advance_observed(&dark, 0.5, config(0.25), || Ok(()))
            .is_err()
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CELL).unwrap(), 300.0);
    assert_eq!(dark.time_s(), 0.0);
    assert_eq!(dark.energy().incident_j, [0.0]);
    assert_eq!(dark.response_evaluations(), 0);

    // The complete source and opaque n/k remain valid as the spectral band
    // moves from 8–12 to 18–22 μm. A retained 10 μm beam must nevertheless
    // refuse at the first unsupported endpoint, even when it carries no power.
    let (card, original) = fixture(|claims| {
        for name in SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES {
            let source = claims.iter_mut().find(|c| c.key.name() == name).unwrap();
            if let PropertyValue::Curve { knots, .. } = &mut source.value {
                knots[1].1 += 10e-6;
                knots[2].1 += 10e-6;
            }
        }
    });
    let moving = Arc::new(
        SurfaceConductorLaw::try_material_temperature(
            card,
            original,
            "thermal",
            (300.0, 600.0),
            501.0,
            28,
        )
        .unwrap(),
    );
    assert!(matches!(
        moving
            .at_temperature(350.0)
            .unwrap()
            .complex_ior(WAVELENGTH),
        Err(ConductorError::InvalidWavelength)
    ));
    let initial = native(300.0, 100.0);
    let bare = illumination(&initial, Arc::clone(&moving), 0.0);
    let film = CoherentAbsorbingCoating::try_new_in_band(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(4000.0 + 2500.0 * i as f64, 1.8, 0.1).unwrap()
        }),
        1e-6,
        2501.0,
        ConductorSource::try_new(ID, ConductorDataStatus::Representative).unwrap(),
    )
    .unwrap();
    let coated = SurfaceOpticalHeating::try_new_coated_with_stateful(
        &[],
        &[StatefulCoatedSurfaceIlluminationSource {
            source_id: 2,
            films: vec![StatefulCoatingLayer {
                receiver: receiver(),
                optics: CoatingOpticsModel::Fixed(film),
            }],
            substrate_receiver: receiver(),
            substrate: SubstrateOpticsModel::Law(moving),
            incidence: SurfaceIncidence {
                wavelengths_nm: [WAVELENGTH],
                incident_eta: [1.0],
                incidence_cosine: 1.0,
                s_power_fraction: [0.5],
            },
            power: SpectralPowerSchedule::try_new(&[], 1).unwrap(),
        }],
        &initial,
        SurfaceHeatingBudget {
            max_sources: 1,
            max_total_power_segments: 1,
            max_source_evaluations: 128,
        },
        28,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
    .unwrap();
    for observer in [bare, coated] {
        assert!(matches!(
            initial.advance_observed(&observer, 0.5, config(0.1), || Ok(())),
            Err(SurfaceHeatingError::Optics(
                ConductorError::InvalidWavelength
            ))
        ));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.temperature_k(CELL).unwrap(), 300.0);
        assert_eq!(observer.time_s(), 0.0);
        assert_eq!(observer.energy().incident_j, [0.0]);
        assert_eq!(observer.energy().thermal_input_j, 0.0);
        assert_eq!(observer.response_evaluations(), 0);
    }
}
