use super::*;
use crate::dielectric::{DielectricGlass, DielectricIor, SellmeierIor};
use crate::thermal_surface::CoatingExchangeQuadrature;
use crate::thermal_surface::heating::coating_exchange::{
    CoatingExchangeBath, CoatingExchangeSource, IsothermalCoatingStack, NativeCoatingExchange,
    NativeCoatingExchangeBudget, NativeCoatingExchangeError,
};
use crate::thermal_surface::heating::{
    NativeCoatingEmission, NativeCoatingEmissionError, SurfaceHeatingReceiver,
};
use crate::tracer::Material;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
    SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, ScalarAdmissibility, ScalarPropertyRequirement,
    SpectralSellmeierStatePoint, SpectralTabulatedDielectricStatePoint,
    VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use std::collections::BTreeSet;

mod illumination;
mod schott;

const ID: ContentHash = ContentHash([213; 32]);
const CELL: CellId = CellId(1);
const MATERIAL: MaterialId = MaterialId(7);
const DEPTH: f64 = 0.02;
const AREA: f64 = 0.02;
const Q: CoatingExchangeQuadrature = CoatingExchangeQuadrature {
    wavelength_bins: 8,
    cosine_bins: 4,
};

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
            knots: [300.0, 450.0, 600.0].into_iter().zip(values).collect(),
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
            source: format!("synthetic native slab temperature claim: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}
fn index(wavelength: f64, t: f64) -> f64 {
    1.4 + 0.02 * (wavelength - 8000.0) / 6000.0 + 0.0001 * (t - 300.0)
}
fn absorption(wavelength: f64, t: f64) -> f64 {
    let baseline = if t <= 450.0 {
        20.0 + 60.0 * (t - 300.0) / 150.0
    } else {
        80.0 - 40.0 * (t - 450.0) / 150.0
    };
    let feature = (1.0 - (wavelength - 11_000.0).abs() / 750.0).max(0.0);
    baseline * (1.0 + 2.0 * feature)
}
fn fixture(
    change: impl FnOnce(&mut Vec<PropertyClaim>),
) -> (Arc<MaterialCard>, SpectralOpticalStatePoint) {
    let mut claims = Vec::new();
    for i in 0..9 {
        let wavelength = 8000.0 + 750.0 * i as f64;
        claims.extend([
            claim(
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                [wavelength * 1e-9; 3],
            ),
            claim(
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[i],
                QuantitySpec::dimensional(Dims::NONE),
                [300.0, 450.0, 600.0].map(|t| index(wavelength, t)),
            ),
            claim(
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                [300.0, 450.0, 600.0].map(|t| absorption(wavelength, t)),
            ),
        ]);
    }
    claims.push(claim(
        "density",
        QuantitySpec::dimensional(Density::DIMS),
        [2500.0, 2400.0, 2300.0],
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
                chemistry: "synthetic dielectric slab".into(),
                phase: "solid".into(),
                process: "native opacity feedback regression".into(),
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
        SpectralOpticalStatePoint::TabulatedDielectric(
            SpectralTabulatedDielectricStatePoint::try_from_material(resolved).unwrap(),
        ),
    )
}
fn law(card: Arc<MaterialCard>, source: SpectralOpticalStatePoint) -> Arc<MaterialSlabLaw> {
    Arc::new(
        MaterialSlabLaw::try_new(card, source, "thermal", (300.0, 600.0), Some(751.0), 28).unwrap(),
    )
}
fn glass(binding: MaterialOpticalBinding) -> DielectricGlass {
    let Material::Dielectric { glass, surface } = binding.material() else {
        panic!("dielectric slab")
    };
    assert!(surface.is_delta());
    assert!(surface.coating().is_none());
    glass
}
fn vacuum() -> DielectricIor {
    SellmeierIor::try_new_in_band([0.0; 3], [0.0; 3], (7900.0, 14_100.0))
        .unwrap()
        .into()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * b.abs().max(1e-5),
        "{a:.17e} != {b:.17e}"
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
fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.1,
        energy_tolerance_j: 1e-9,
    }
}
fn native(t: f64) -> DepositionSimulation<'static> {
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
                    source_w: 0.0,
                },
            }),
        }],
        &[(MATERIAL, &SENSIBLE)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
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
fn source(optics: Arc<MaterialSlabLaw>, bath: f64) -> CoatingExchangeSource {
    CoatingExchangeSource {
        source_id: 10,
        source_identity: ID,
        receiver: receiver(),
        area_m2: AREA,
        stack: IsothermalCoatingStack::MaterialIncoherentSlab {
            law: optics,
            thickness_m: DEPTH,
            opposite_ior: vacuum(),
        },
        bath: CoatingExchangeBath::Blackbody {
            temperature_k: bath,
        },
        quadrature: Q,
    }
}
fn budget() -> NativeCoatingExchangeBudget {
    NativeCoatingExchangeBudget {
        max_sources: 2,
        max_total_optical_knots: 28,
        max_source_evaluations: 128,
        max_optical_evaluations: 4096,
    }
}

// Independent incoherent geometric series for a slab between two vacuum ports.
// Ordinary real-index Fresnel powers and the refracted bulk path are evaluated
// directly, without calling the production slab/absorption/thermal owners.
fn absorptance(wavelength: f64, cosine: f64, optical_t: f64) -> [f64; 2] {
    let n = index(wavelength, optical_t);
    let ct = (1.0 - (1.0 - cosine * cosine) / (n * n)).sqrt();
    let s = ((cosine - n * ct) / (cosine + n * ct)).powi(2);
    let p = ((n * cosine - ct) / (n * cosine + ct)).powi(2);
    let a = (-absorption(wavelength, optical_t) * DEPTH / ct).exp();
    [s, p].map(|r| (1.0 - r) * (1.0 - a) / (1.0 - r * a))
}
fn planck(wavelength: f64, t: f64) -> f64 {
    let lambda = wavelength * 1e-9;
    let (h, c, k) = (6.62607015e-34, 299_792_458.0, 1.380649e-23);
    2.0 * h * c * c / lambda.powi(5) / (h * c / (lambda * k * t)).exp_m1() * 1e-9
}
fn outward(body_t: f64, optical_t: f64, bath_t: f64) -> f64 {
    let width = 6000.0 / Q.wavelength_bins as f64;
    let mut watts = 0.0;
    for i in 0..Q.wavelength_bins {
        let wavelength = 8000.0 + (i as f64 + 0.5) * width;
        for j in 0..Q.cosine_bins {
            let cosine = (j as f64 + 0.5) / Q.cosine_bins as f64;
            let [s, p] = absorptance(wavelength, cosine, optical_t);
            watts += AREA * core::f64::consts::PI * width / Q.cosine_bins as f64
                * cosine
                * (s + p)
                * (planck(wavelength, body_t) - planck(wavelength, bath_t));
        }
    }
    watts
}

#[test]
fn original_claims_and_spectral_feature_refresh_without_temperature_fitting() {
    let (card, original) = fixture(|claims| {
        let mut replacement = claim(
            SPECTRAL_DIELECTRIC_IOR_PROPERTIES[0],
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
    assert_eq!(law.material_source(), &original);
    assert_eq!(law.budget_size(), 28);
    let refreshed = law.source_at_temperature(525.0).unwrap();
    assert_eq!(
        refreshed,
        original.at_query_point(&card, &point(525.0)).unwrap()
    );
    assert_eq!(refreshed.resolved().query().axes()["pressure"], 100_000.0);
    near(
        refreshed.resolved().property("density").unwrap().value_si(),
        2350.0,
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
    for t in [300.0, 375.0, 450.0, 525.0, 600.0] {
        let glass = glass(law.at_temperature(t).unwrap());
        for wavelength in [8001.0, 10_625.0, 11_000.0, 13_999.0] {
            near(glass.ior().eval(wavelength).unwrap(), index(wavelength, t));
            near(
                glass.absorption_ref().extinction_per_m(wavelength).unwrap(),
                absorption(wavelength, t),
            );
        }
        assert!(glass.ior().eval(550.0).is_err());
        assert!(glass.absorption_ref().extinction_per_m(550.0).is_err());
    }
    for (axis, limit, gap) in [
        ("thermal", 27, Some(751.0)),
        ("pressure", 28, Some(751.0)),
        ("thermal", 28, None),
    ] {
        assert!(
            MaterialSlabLaw::try_new(
                Arc::clone(&card),
                original.clone(),
                axis,
                (300.0, 600.0),
                gap,
                limit
            )
            .is_err()
        );
    }
    assert!(law.at_temperature(601.0).is_err());
    assert_eq!(law.material_source(), &original);
}

#[test]
fn analytic_sellmeier_temperature_and_source_domain_stay_analytic() {
    let (card, original) = fixture(|claims| {
        for i in 0..3 {
            claims.push(claim(
                VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[i],
                QuantitySpec::dimensional(Dims::NONE),
                if i == 0 { [1.1, 1.2, 1.3] } else { [0.0; 3] },
            ));
            claims.push(claim(
                VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
                if i == 0 {
                    [1e-12, 2e-12, 3e-12]
                } else {
                    [0.0; 3]
                },
            ));
        }
        for (i, values) in [[8e-6, 8.5e-6, 9e-6], [14e-6; 3]].into_iter().enumerate() {
            claims.push(claim(
                SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[i],
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                values,
            ));
        }
    });
    let original = SpectralOpticalStatePoint::DielectricSellmeier(
        SpectralSellmeierStatePoint::try_from_material(original.resolved().clone()).unwrap(),
    );
    let law =
        MaterialSlabLaw::try_new(card, original, "thermal", (300.0, 600.0), None, 36).unwrap();
    let warm = glass(law.at_temperature(375.0).unwrap());
    assert!(matches!(warm.ior(), DielectricIor::Sellmeier(_)));
    let wavelength: f64 = 10_625.0;
    let square = (wavelength * 1e-9).powi(2);
    near(
        warm.ior().eval(wavelength).unwrap(),
        (1.0 + 1.15 * square / (square - 1.5e-12)).sqrt(),
    );
    near(
        warm.absorption_ref().extinction_per_m(wavelength).unwrap(),
        absorption(wavelength, 375.0),
    );
    assert!(warm.ior().eval(8100.0).is_err());
    assert!(
        glass(law.at_temperature(300.0).unwrap())
            .ior()
            .eval(8100.0)
            .is_ok()
    );
    near(warm.wavelength_bounds_nm().unwrap().0, 8250.0);
}

#[test]
fn native_opacity_emission_and_accepted_heat_follow_original_temperature_claims() {
    let (card, original) = fixture(|_| {});
    let law = law(card, original);
    let thermal = native(400.0);
    let initial = NativeCoatingExchange::try_new(
        &thermal,
        &[source(Arc::clone(&law), 600.0)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    let emitter = initial.emitter(10).unwrap();
    let old_optics = emitter.optical_identity();
    let sample = emitter
        .evaluate(10_625.0, 0.6)
        .unwrap()
        .polarized_radiance_w_m2_sr_nm();
    let a = absorptance(10_625.0, 0.6, 400.0);
    for pol in 0..2 {
        near(sample[pol], a[pol] * planck(10_625.0, 400.0) * 0.5);
    }
    near(initial.sources_w()[&CELL], -outward(400.0, 400.0, 600.0));

    // Explicit native caloric storage is 1 J/K here. The native monotonicity
    // owner combines max_dt_s harmonically with the positive radiative secant,
    // even when max_dt_s is the smaller limit. Independently integrate the
    // physical reference on those accepted intervals, including the final cut.
    let reference = |feedback: bool| {
        let (mut time, mut temperature, mut steps) = (0.0, 400.0, 0);
        while time < 0.2 {
            let optical_t = if feedback { temperature } else { 400.0 };
            let power = outward(temperature, optical_t, 600.0);
            let conductance = power / (temperature - 600.0);
            let maximum = config().max_dt_s;
            let limit = maximum / (1.0 + maximum * conductance / 0.8);
            let dt = limit.min(0.2 - time);
            time = if dt == 0.2 - time { 0.2 } else { time + dt };
            temperature -= dt * power;
            steps += 1;
            assert!(steps <= 16);
        }
        (temperature, steps)
    };
    let (expected_t, expected_steps) = reference(true);
    let (frozen_t, _) = reference(false);
    assert!((expected_t - frozen_t).abs() > 1e-5);
    let (advanced, observed, audit) = thermal
        .advance_observed(&initial, 0.2, config(), || Ok(()))
        .unwrap();
    let actual_t = advanced.temperature_k(CELL).unwrap();
    near(actual_t, expected_t);
    near(advanced.energy_j() - thermal.energy_j(), expected_t - 400.0);
    near(observed.energy().thermal_external_j, expected_t - 400.0);
    assert_eq!(audit.substeps, expected_steps);
    near(audit.external_heat_j, -observed.energy().outward_j);
    near(
        observed.sources_w()[&CELL],
        -outward(actual_t, actual_t, 600.0),
    );
    assert_ne!(observed.emitter(10).unwrap().optical_identity(), old_optics);
    assert_eq!(
        observed.emitter(10).unwrap().film_temperature_k(0),
        Some(actual_t)
    );
    assert_eq!(
        observed.work().optical_evaluations,
        expected_steps * Q.wavelength_bins * Q.cosine_bins
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(initial.emitter(10).unwrap().optical_identity(), old_optics);
    assert_eq!(initial.energy().outward_j, 0.0);
}

#[test]
fn transparent_and_isothermal_material_slabs_exchange_no_heat() {
    for (transparent, bath) in [(true, 600.0), (false, 400.0)] {
        let (card, original) = fixture(|claims| {
            if transparent {
                for source in claims {
                    if SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES.contains(&source.key.name())
                    {
                        if let PropertyValue::Curve { knots, .. } = &mut source.value {
                            for (_, value) in knots {
                                *value = 0.0;
                            }
                        }
                    }
                }
            }
        });
        let law = law(card, original);
        let thermal = native(400.0);
        // Two exposed faces share one law and one actual caloric store; the
        // complete material property bundle is charged once, not per patch.
        let front = source(Arc::clone(&law), bath);
        let mut rear = front.clone();
        rear.source_id = 20;
        let initial =
            NativeCoatingExchange::try_new(&thermal, &[front, rear], budget(), || Ok(())).unwrap();
        let (advanced, observed, audit) = thermal
            .advance_observed(&initial, 0.2, config(), || Ok(()))
            .unwrap();
        assert_eq!(advanced.temperature_k(CELL).unwrap(), 400.0);
        assert_eq!(observed.energy().outward_j, 0.0);
        assert_eq!(audit.external_heat_j, 0.0);
        assert_eq!(observed.work().optical_evaluations, 0);
        if transparent {
            assert_eq!(
                observed
                    .emitter(10)
                    .unwrap()
                    .evaluate(10_625.0, 0.6)
                    .unwrap()
                    .polarized_radiance_w_m2_sr_nm(),
                [0.0; 2]
            );
        }
    }
}

#[test]
fn unsupported_endpoint_keeps_accepted_heat_emission_and_cache_unchanged() {
    let (card, original) = fixture(|claims| {
        let density = claims
            .iter_mut()
            .find(|source| source.key.name() == "density")
            .unwrap();
        density.validity = density.validity.clone().with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
            400.0,
        );
    });
    let law = law(card, original);
    assert!(matches!(
        law.at_temperature(400.1),
        Err(SurfaceHeatingError::MaterialSource(_))
    ));
    let thermal = native(400.0);
    let patch = source(Arc::clone(&law), 600.0);
    let initial =
        NativeCoatingExchange::try_new(&thermal, &[patch.clone()], budget(), || Ok(())).unwrap();
    let identity = initial.emitter(10).unwrap().identity();
    let optics = initial.emitter(10).unwrap().optical_identity();
    let watts = initial.sources_w().clone();
    let work = initial.work();
    assert!(matches!(
        thermal.advance_observed(&initial, 0.1, config(), || Ok(())),
        Err(NativeCoatingExchangeError::Binding(
            NativeCoatingEmissionError::Binding(SurfaceHeatingError::MaterialSource(_))
        ))
    ));
    assert_eq!(thermal.temperature_k(CELL).unwrap(), 400.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.emitter(10).unwrap().identity(), identity);
    assert_eq!(initial.emitter(10).unwrap().optical_identity(), optics);
    assert_eq!(initial.sources_w(), &watts);
    assert_eq!(initial.work(), work);
    assert_eq!(initial.energy().outward_j, 0.0);
    let mut insufficient = budget();
    insufficient.max_total_optical_knots = 27;
    assert!(matches!(
        NativeCoatingExchange::try_new(&thermal, &[patch], insufficient, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    assert!(
        NativeCoatingEmission::try_material_incoherent_slab(
            &thermal,
            receiver(),
            &law,
            DEPTH,
            vacuum(),
            ContentHash([0; 32]),
            28,
            || Ok(())
        )
        .is_err()
    );
}
