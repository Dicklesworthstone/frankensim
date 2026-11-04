use super::*;
use crate::conductor::{CoatingLayer, ConductorError, ConductorIorSample, ConductorOptics};
use crate::dielectric::{DielectricIor, DielectricIorSample, TabulatedIor};
use crate::thermal_surface::CoatingExchangeQuadrature;
use crate::thermal_surface::heating::coating_exchange::{
    CoatingExchangeBath, CoatingExchangeSource, IsothermalCoatingStack, NativeCoatingExchange,
    NativeCoatingExchangeBudget,
};
use crate::thermal_surface::heating::{
    CoatingOpticsModel, CoatingSurfaceSource, NativeCoatingEmission,
    StatefulCoatedSurfaceIlluminationSource, StatefulCoatingLayer,
    StatefulTransmittingSurfaceIlluminationSource, SubstrateOpticsModel, SurfaceHeatingBudget,
    SurfaceHeatingReceiver, SurfaceHeatingTolerance, SurfaceIncidence, SurfaceOpticalHeating,
    TransmittingSurfacePowerResponse,
};
use crate::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_blake3::ContentHash;
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
use fs_math::c64::C64;
use fs_qty::{Density, Dims, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use std::collections::BTreeSet;

const ID: ContentHash = ContentHash([211; 32]);
const CELL: CellId = CellId(1);
const MATERIAL: MaterialId = MaterialId(7);
const WAVELENGTH: f64 = 10_000.0;
const DEPTH: f64 = 2e-6;

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
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
            source: format!("synthetic finite-film temperature claim: {name}"),
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
                [1.5, 2.5, 1.75],
            ),
            claim(
                SPECTRAL_CONDUCTOR_K_PROPERTIES[i],
                SPECTRAL_CONDUCTOR_K_QUANTITY,
                [0.0, 0.4, 0.1],
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
                chemistry: "synthetic finite absorbing film".into(),
                phase: "solid".into(),
                process: "native photothermal regression".into(),
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
    (
        card,
        SpectralConductorStatePoint::try_from_material(resolved).unwrap(),
    )
}
fn law(card: Arc<MaterialCard>, original: SpectralConductorStatePoint) -> SurfaceCoatingLaw {
    SurfaceCoatingLaw::try_material_temperature(
        card,
        original,
        "thermal",
        (300.0, 600.0),
        DEPTH,
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
fn constants(t: f64) -> (f64, f64) {
    if t <= 400.0 {
        let f = (t - 300.0) / 100.0;
        (1.5 + f, 0.4 * f)
    } else {
        let f = (t - 400.0) / 200.0;
        (2.5 - 0.75 * f, 0.4 - 0.3 * f)
    }
}

// Independent single-slab amplitudes from the two-interface geometric series,
// not the native positive-loss modes. Byrnes (2016), sections 3, 4.1 and 4.4:
// https://arxiv.org/abs/1603.02720. Vacuum ports and normal incidence.
fn slab_power(t: f64) -> [f64; 3] {
    let (n, k) = constants(t);
    let index = C64::new(n, k);
    let phase = core::f64::consts::TAU * DEPTH / (WAVELENGTH * 1e-9);
    let z = C64::new((phase * n).cos(), (phase * n).sin()).scale((-phase * k).exp());
    let r01 = (C64::ONE - index) / (C64::ONE + index);
    let denominator = C64::ONE - r01 * r01 * z * z;
    let r = r01 * (C64::ONE - z * z) / denominator;
    let transmission =
        index.scale(4.0) * z / ((C64::ONE + index) * (C64::ONE + index) * denominator);
    let reflected = r.norm_sq();
    let transmitted = transmission.norm_sq();
    [
        reflected,
        transmitted,
        if k == 0.0 {
            0.0
        } else {
            1.0 - reflected - transmitted
        },
    ]
}
fn planck(t: f64) -> f64 {
    let lambda = WAVELENGTH * 1e-9;
    let (h, c, k) = (6.62607015e-34, 299_792_458.0, 1.380649e-23);
    2.0 * h * c * c / lambda.powi(5) / (h * c / (lambda * k * t)).exp_m1() * 1e-9
}
fn vacuum() -> DielectricIor {
    TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(4000.0 + 2500.0 * i as f64, 1.0).unwrap()
        }),
        2501.0,
    )
    .unwrap()
    .into()
}
fn substrate() -> ConductorOptics {
    ConductorOptics::try_new_in_band(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new_in_band(4000.0 + 2500.0 * i as f64, 0.3, 3.0).unwrap()
        }),
        2501.0,
        ConductorSource::try_new(ID, ConductorDataStatus::Representative).unwrap(),
    )
    .unwrap()
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
fn film(optics: Arc<SurfaceCoatingLaw>) -> StatefulCoatingLayer {
    StatefulCoatingLayer {
        receiver: receiver(),
        optics: CoatingOpticsModel::Law(optics),
    }
}
fn illumination(
    native: &DepositionSimulation<'_>,
    optics: Arc<SurfaceCoatingLaw>,
    power_w: f64,
    opaque: bool,
    maximum: usize,
) -> Result<SurfaceOpticalHeating<1>, SurfaceHeatingError> {
    let films = vec![film(optics)];
    let incidence = SurfaceIncidence {
        wavelengths_nm: [WAVELENGTH],
        incident_eta: [1.0],
        incidence_cosine: 1.0,
        s_power_fraction: [0.5],
    };
    let power = SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_power_w: [power_w],
            end_power_w: [power_w],
        }],
        1,
    )
    .unwrap();
    let transmitting = StatefulTransmittingSurfaceIlluminationSource {
        source_id: 1,
        films: films.clone(),
        incidence,
        transmitted_eta: [1.0],
        reverse: false,
        source_identity: ID,
        power: power.clone(),
    };
    let coated = StatefulCoatedSurfaceIlluminationSource {
        source_id: 1,
        films,
        substrate_receiver: receiver(),
        substrate: SubstrateOpticsModel::Fixed(substrate()),
        incidence,
        power,
    };
    let source = if opaque {
        CoatingSurfaceSource::StatefulOpaque(&coated)
    } else {
        CoatingSurfaceSource::StatefulTransmitting(&transmitting)
    };
    SurfaceOpticalHeating::try_new_coating_stacks(
        &[source],
        native,
        SurfaceHeatingBudget {
            max_sources: 1,
            max_total_power_segments: 1,
            max_source_evaluations: 128,
        },
        maximum,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
}

#[test]
fn complete_original_claims_refresh_without_fitted_temperature_rows() {
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
    assert_eq!(law.caloric_identity(), None);
    assert_eq!(law.specific_enthalpy_bounds_j_kg(), None);
    let Model::Material(retained) = &law.model else {
        panic!("material film")
    };
    let cloned = law.clone();
    let Model::Material(shared) = &cloned.model else {
        panic!("shared material film")
    };
    assert!(Arc::ptr_eq(retained, shared));
    let refreshed = retained.source_at_temperature(450.0).unwrap();
    assert_eq!(
        refreshed,
        original.at_query_point(&card, &point(450.0)).unwrap()
    );
    for property in original.resolved().properties() {
        assert_eq!(
            refreshed
                .resolved()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected,
            property.answer().receipt.selected
        );
    }
    near(
        refreshed.resolved().property("density").unwrap().value_si(),
        1875.0,
    );
    assert_eq!(refreshed.resolved().query().axes()["pressure"], 100_000.0);
    for t in [300.0, 350.0, 400.0, 450.0, 600.0] {
        let coating = law.at_temperature(t).unwrap();
        let (n, k) = constants(t);
        for sample in coating.samples() {
            near(sample.eta(), n);
            near(sample.k(), k);
        }
        assert_eq!(coating.thickness_m(), DEPTH);
        assert_eq!(
            coating.source().status(),
            ConductorDataStatus::MaterialStateClaims
        );
    }
    for (axis, maximum) in [("thermal", 27), ("pressure", 28), ("absent", 28)] {
        assert!(
            SurfaceCoatingLaw::try_material_temperature(
                Arc::clone(&card),
                original.clone(),
                axis,
                (300.0, 600.0),
                DEPTH,
                501.0,
                maximum
            )
            .is_err()
        );
    }
    assert!(law.at_temperature(299.0).is_err());
    assert!(law.at_specific_enthalpy(350.0, 350.0).is_err());
}

#[test]
fn transparent_film_becomes_absorbing_and_heat_matches_independent_slab_and_emission() {
    let (card, original) = fixture(|_| {});
    let law = Arc::new(law(card, original));
    let cold = native(300.0, 100.0);
    assert!(matches!(
        illumination(&cold, Arc::clone(&law), 100.0, false, 27),
        Err(SurfaceHeatingError::Budget)
    ));
    let light = illumination(&cold, Arc::clone(&law), 100.0, false, 28).unwrap();
    let (first, light, _) = cold
        .advance_observed(&light, 0.1, config(0.1), || Ok(()))
        .unwrap();
    // A perfectly transparent film receives no optical heat; the separate
    // authored thermal source first brings it into an absorbing source state.
    assert_eq!(light.energy().absorbed_j, [0.0]);
    near(first.temperature_k(CELL).unwrap(), 310.0);
    let [r, t, a] = slab_power(310.0);
    assert!(a > 0.0 && t > 0.0);
    let prepared = law.at_temperature(310.0).unwrap();
    let response = TransmittingSurfacePowerResponse::from_dielectric_stack(
        [WAVELENGTH],
        &[CoatingLayer::Absorbing(&prepared)],
        [1.0],
        [1.0],
        1.0,
        false,
        [0.5],
        ID,
    )
    .unwrap();
    near(response.total().reflected_fraction()[0], r);
    near(response.total().transmitted_fraction()[0], t);
    near(response.film_absorbed_fractions()[0][0], a);
    let (second, light, audit) = first
        .advance_observed(&light, 0.2, config(0.1), || Ok(()))
        .unwrap();
    let final_t = 320.0 + 10.0 * a;
    near(second.temperature_k(CELL).unwrap(), final_t);
    near(light.energy().absorbed_j[0], 10.0 * a);
    near(light.energy().thermal_input_j, 10.0 * a);
    near(audit.external_heat_j, 10.0 * (1.0 + a));
    near(
        light.energy().incident_j[0],
        light.energy().reflected_j[0]
            + light.energy().transmitted_j[0]
            + light.energy().absorbed_j[0],
    );
    near(
        second.energy_j() - cold.energy_j(),
        20.0 + light.energy().thermal_input_j,
    );
    let layers = [film(law)];
    assert!(
        NativeCoatingEmission::try_transmitting(&second, &layers, vacuum(), false, ID, 27, || Ok(
            ()
        ))
        .is_err()
    );
    for reverse in [false, true] {
        let emitter = NativeCoatingEmission::try_transmitting(
            &second,
            &layers,
            vacuum(),
            reverse,
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
            slab_power(final_t)[2] * planck(final_t),
        );
        assert!(emitter.emitter().radiance_w_m2_sr_nm(550.0, 1.0).is_err());
    }
}

#[test]
fn finite_film_infrared_exchange_retains_transparency_and_thermal_equilibrium() {
    let (card, original) = fixture(|_| {});
    let law = Arc::new(law(card, original));
    for (initial_t, bath_t) in [
        (300.0, 500.0),
        (400.0, 300.0),
        (400.0, 400.0),
        (400.0, 500.0),
    ] {
        let initial = native(initial_t, 0.0);
        let patch = CoatingExchangeSource {
            source_id: 1,
            source_identity: ID,
            receiver: receiver(),
            area_m2: 0.01,
            stack: IsothermalCoatingStack::Transmitting {
                films: vec![CoatingOpticsModel::Law(Arc::clone(&law))],
                opposite_ior: vacuum(),
                reverse: false,
            },
            bath: CoatingExchangeBath::Blackbody {
                temperature_k: bath_t,
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
        assert!(
            NativeCoatingExchange::try_new(
                &initial,
                &[patch.clone()],
                NativeCoatingExchangeBudget {
                    max_total_optical_knots: 27,
                    ..budget
                },
                || Ok(())
            )
            .is_err()
        );
        let exchange =
            NativeCoatingExchange::try_new(&initial, &[patch], budget, || Ok(())).unwrap();
        let (changed, accepted, audit) = initial
            .advance_observed(&exchange, 0.1, config(0.05), || Ok(()))
            .unwrap();
        let delta = changed.temperature_k(CELL).unwrap() - initial_t;
        near(
            changed.energy_j() - initial.energy_j(),
            -accepted.energy().outward_j,
        );
        near(audit.external_heat_j, accepted.energy().thermal_external_j);
        if initial_t == 300.0 || initial_t == bath_t {
            assert_eq!(delta, 0.0);
            assert_eq!(accepted.energy().outward_j, 0.0);
        } else {
            assert_eq!(delta.signum(), (bath_t - initial_t).signum());
            let t = changed.temperature_k(CELL).unwrap();
            near(
                accepted
                    .emitter(1)
                    .unwrap()
                    .radiance_w_m2_sr_nm(WAVELENGTH, 1.0)
                    .unwrap(),
                slab_power(t)[2] * planck(t),
            );
        }
    }
}

#[test]
fn expired_auxiliary_claims_refuse_dark_opaque_and_transmitting_candidates() {
    let (card, original) = fixture(|claims| {
        let density = claims
            .iter_mut()
            .find(|c| c.key.name() == "density")
            .unwrap();
        density.validity = density.validity.clone().with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
            315.0,
        );
        let mut replacement = claim(
            "density",
            QuantitySpec::dimensional(Density::DIMS),
            [1600.0; 3],
        );
        replacement.validity = replacement.validity.with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            310.0,
            600.0,
        );
        claims.push(replacement);
    });
    let expired = Arc::new(law(card, original));
    assert!(matches!(
        expired.at_temperature(320.0),
        Err(SurfaceHeatingError::MaterialSource(_))
    ));
    let initial = native(300.0, 100.0);
    for opaque in [false, true] {
        let dark = illumination(&initial, Arc::clone(&expired), 0.0, opaque, 28).unwrap();
        assert!(matches!(
            initial.advance_observed(&dark, 0.5, config(0.1), || Ok(())),
            Err(SurfaceHeatingError::MaterialSource(_))
        ));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.temperature_k(CELL).unwrap(), 300.0);
        assert_eq!(dark.time_s(), 0.0);
        assert_eq!(dark.energy().thermal_input_j, 0.0);
        assert_eq!(dark.response_evaluations(), 0);
    }
    let (card, original) = fixture(|_| {});
    let light = illumination(&initial, Arc::new(law(card, original)), 100.0, false, 28).unwrap();
    assert!(
        initial
            .advance_observed(&light, 0.2, config(0.1), || Err(DepositionError::Cancelled))
            .is_err()
    );
    assert_eq!(light.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CELL).unwrap(), 300.0);
}

#[test]
fn moving_film_wavelengths_refuse_unsupported_dark_endpoints() {
    let (card, original) = fixture(|claims| {
        for name in SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES {
            let source = claims.iter_mut().find(|c| c.key.name() == name).unwrap();
            if let PropertyValue::Curve { knots, .. } = &mut source.value {
                knots[1].1 += 10e-6;
                knots[2].1 += 10e-6;
            }
        }
    });
    let law = Arc::new(law(card, original));
    assert!(matches!(
        law.at_temperature(350.0).unwrap().complex_ior(WAVELENGTH),
        Err(ConductorError::InvalidWavelength)
    ));
    let initial = native(300.0, 100.0);
    for opaque in [false, true] {
        let dark = illumination(&initial, Arc::clone(&law), 0.0, opaque, 28).unwrap();
        assert!(matches!(
            initial.advance_observed(&dark, 0.5, config(0.1), || Ok(())),
            Err(SurfaceHeatingError::Optics(
                ConductorError::InvalidWavelength
            ))
        ));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.temperature_k(CELL).unwrap(), 300.0);
        assert_eq!(dark.time_s(), 0.0);
        assert_eq!(dark.energy().incident_j, [0.0]);
        assert_eq!(dark.energy().thermal_input_j, 0.0);
        assert_eq!(dark.response_evaluations(), 0);
    }
}
