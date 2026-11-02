//! G1/G3: native storage, coherent cooling, bath events and independent refinement.
use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack,
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::dielectric::{CauchyIor, DielectricIor};
use fs_render::thermal_surface::CoatingExchangeQuadrature;
use fs_render::thermal_surface::heating::coating_exchange::*;
use fs_render::thermal_surface::heating::{
    CoatingOpticsModel, SubstrateOpticsModel, SurfaceCoatingLaw, SurfaceHeatingReceiver,
};
use fs_thermal::enthalpy::{Reservoir, ThermalLink};
use std::sync::Arc;

mod infrared;

struct Hot;
impl ThermalLaw for Hot {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (-9000.0..=100_000.0).contains(&h) {
            Ok(1000.0 + h / 10.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.1
    }
}
static HOT: Hot = Hot;
struct HotLatent;
impl ThermalLaw for HotLatent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (10.0..=110.0).contains(&h) {
            Ok(2500.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([231; 32])
    }
}
static HOT_LATENT: HotLatent = HotLatent;
fn authority() -> ConductorSource {
    ConductorSource::try_new(ContentHash([232; 32]), ConductorDataStatus::Representative).unwrap()
}
fn film(n: f64, k: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        std::array::from_fn(|i| CoatingIorSample::try_new(380.0 + 50.0 * i as f64, n, k).unwrap()),
        100e-9,
        authority(),
    )
    .unwrap()
}
fn metal() -> ConductorOptics {
    ConductorOptics::try_new(
        std::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.5, 3.0).unwrap()
        }),
        authority(),
    )
    .unwrap()
}
fn temperature_film(t: f64) -> CoherentAbsorbingCoating {
    film(1.5 + 0.0001 * (t - 1000.0), 0.05 + 0.00004 * (t - 1000.0))
}
fn temperature_model() -> CoatingOpticsModel {
    CoatingOpticsModel::Law(Arc::new(
        SurfaceCoatingLaw::try_temperature(
            authority(),
            &[
                (1000.0, temperature_film(1000.0)),
                (5000.0, temperature_film(5000.0)),
            ],
            2,
        )
        .unwrap(),
    ))
}
fn specimen(t: f64, mass: f64) -> DepositedCell {
    let mut c = cell(11, mass, 1.0, 0.0);
    c.specific_h_j_kg = 10.0 * (t - 1000.0);
    c.boundary.ambient_k = 1800.0;
    c
}
fn source(area: f64, bath: CoatingExchangeBath, dynamic: bool) -> CoatingExchangeSource {
    CoatingExchangeSource {
        source_id: 10,
        source_identity: ContentHash([233; 32]),
        receiver: SurfaceHeatingReceiver {
            thermal_cell: CellId(11),
            material_id: MaterialId(7),
            temperature_bounds_k: (100.0, 5000.0),
        },
        area_m2: area,
        bath,
        stack: IsothermalCoatingStack::Opaque {
            films: vec![if dynamic {
                temperature_model()
            } else {
                CoatingOpticsModel::Fixed(film(1.8, 0.2))
            }],
            substrate: SubstrateOpticsModel::Fixed(metal()),
        },
        quadrature: CoatingExchangeQuadrature {
            wavelength_bins: 8,
            cosine_bins: 6,
        },
    }
}
fn budget() -> NativeCoatingExchangeBudget {
    NativeCoatingExchangeBudget {
        max_sources: 8,
        max_total_optical_knots: 16,
        max_source_evaluations: 10_000,
        max_optical_evaluations: 1_000_000,
    }
}
fn exchange(
    native: &DepositionSimulation<'_>,
    sources: &[CoatingExchangeSource],
) -> NativeCoatingExchange {
    NativeCoatingExchange::try_new(native, sources, budget(), || Ok(())).unwrap()
}
fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn cooling_warming_equilibrium_and_vacuum_use_native_storage_without_bath_overshoot() {
    for (body_t, bath_t) in [
        (3100.0, Some(1800.0)),
        (1800.0, Some(3100.0)),
        (2500.0, Some(2500.0)),
        (3100.0, None),
    ] {
        let native = thermal(vec![specimen(body_t, 0.001)], &HOT, vec![]);
        let bath = bath_t.map_or(CoatingExchangeBath::Vacuum, |temperature_k| {
            CoatingExchangeBath::Blackbody { temperature_k }
        });
        let source = source(0.001, bath, false);
        let initial = exchange(&native, &[source]);
        let initial_optics = initial.emitter(10).unwrap().optical_identity();
        let (advanced, observed, audit) = native
            .advance_observed(&initial, 0.2, config(1e6), || Ok(()))
            .unwrap();
        let t = advanced.temperature_k(CellId(11)).unwrap();
        let bath = bath_t.unwrap_or(0.0);
        assert!(t >= body_t.min(bath) && t <= body_t.max(bath));
        if body_t == bath {
            assert_eq!(t, body_t);
            assert_eq!(observed.energy().outward_j, 0.0);
        } else {
            assert_ne!(t, body_t);
        }
        near(
            advanced.energy_j() - native.energy_j(),
            observed.energy().thermal_external_j,
            1e-11,
        );
        near(
            audit.external_heat_j,
            observed.energy().thermal_external_j,
            1e-11,
        );
        near(
            observed.energy().outward_j,
            -observed.energy().thermal_external_j,
            1e-12,
        );
        near(observed.energy().duration_s, 0.2, 1e-14);
        assert_eq!(
            observed.work().optical_evaluations,
            0,
            "fixed optics must reuse its quadrature despite changing T"
        );
        assert_eq!(
            observed.emitter(10).unwrap().optical_identity(),
            initial_optics
        );
        assert_eq!(observed.emitter(10).unwrap().film_temperature_k(0), Some(t));
        assert_eq!(advanced.cells()[&CellId(11)].mass_kg, 0.001);
        assert_eq!(native.time_s(), 0.0);
    }
}

#[test]
fn dynamic_optics_ambient_events_canonical_patches_and_continuation_share_one_native_clock() {
    let mut body = specimen(3100.0, 1.0);
    body.boundary.ambient_k = 3000.0;
    body.boundary.source_w = 2.0;
    let mut later = body.boundary;
    later.ambient_k = 3500.0;
    let native = thermal(
        vec![body],
        &HOT,
        vec![TimedEvent {
            time_s: 0.05,
            event: Event::Boundary(CellId(11), later),
        }],
    );
    let first = source(0.0001, CoatingExchangeBath::NativeAmbient, true);
    let mut second = first.clone();
    second.source_id = 20;
    let original = exchange(&native, &[second.clone(), first.clone()]);
    let ordered = exchange(&native, &[first, second]);
    assert_eq!(original.identity(), ordered.identity());
    let (whole, observation, audit) = native
        .advance_observed(&original, 0.1, config(0.025), || Ok(()))
        .unwrap();
    let (middle, partial, _) = native
        .advance_observed(&ordered, 0.05, config(0.025), || Ok(()))
        .unwrap();
    let (split, resumed, _) = middle
        .advance_observed(&partial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    assert!(middle.temperature_k(CellId(11)).unwrap() < 3100.0);
    assert!(whole.temperature_k(CellId(11)).unwrap() > middle.temperature_k(CellId(11)).unwrap());
    near(whole.energy_j(), split.energy_j(), 2e-13);
    near(
        observation.energy().outward_j,
        resumed.energy().outward_j,
        2e-10,
    );
    near(
        audit.external_heat_j,
        observation.energy().thermal_external_j + 0.2,
        2e-10,
    );
    assert!(observation.work().optical_evaluations > 0);
    assert_eq!(observation.energy().by_source_outward_j.len(), 2);
    near(
        observation.energy().by_source_outward_j[&10],
        observation.energy().by_source_outward_j[&20],
        1e-14,
    );
    assert_eq!(whole.cells()[&CellId(11)].boundary, later);
    assert_eq!(original.time_s(), 0.0);
}

#[test]
fn latent_enthalpy_changes_coherent_cooling_at_constant_temperature() {
    let mut body = cell(11, 1.0, 1.0, 0.0);
    body.specific_h_j_kg = 70.0;
    let native = thermal(vec![body], &HOT_LATENT, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([231; 32]),
            (2500.0, 2500.0),
            &[(10.0, film(1.4, 0.02)), (110.0, film(2.2, 0.6))],
            2,
        )
        .unwrap(),
    );
    let mut source = source(0.00001, CoatingExchangeBath::Vacuum, false);
    source.stack = IsothermalCoatingStack::Transmitting {
        films: vec![
            CoatingOpticsModel::Law(Arc::clone(&law)),
            CoatingOpticsModel::Law(law),
        ],
        opposite_ior: DielectricIor::Cauchy(CauchyIor::try_constant(1.0).unwrap()),
        reverse: true,
    };
    source.receiver.temperature_bounds_k = (2500.0, 2500.0);
    let mut limits = budget();
    limits.max_total_optical_knots = 2;
    let initial =
        NativeCoatingExchange::try_new(&native, &[source.clone()], limits, || Ok(())).unwrap();
    let old_optics = initial.emitter(10).unwrap().optical_identity();
    let (after, evolved, _) = native
        .advance_observed(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(after.temperature_k(CellId(11)).unwrap(), 2500.0);
    assert!(after.cells()[&CellId(11)].specific_h_j_kg < 70.0);
    assert_ne!(evolved.emitter(10).unwrap().optical_identity(), old_optics);
    assert!(evolved.work().optical_evaluations >= 4 * 48);
    near(
        after.energy_j() - native.energy_j(),
        evolved.energy().thermal_external_j,
        1e-11,
    );
    limits.max_total_optical_knots = 1;
    assert!(matches!(
        NativeCoatingExchange::try_new(&native, &[source], limits, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
}

#[test]
fn coverage_events_domains_budgets_and_late_cancellation_refuse_atomically() {
    let body = specimen(3100.0, 1.0);
    let native = thermal(vec![body.clone()], &HOT, vec![]);
    let valid = source(0.001, CoatingExchangeBath::Vacuum, false);
    let initial = exchange(&native, &[valid.clone()]);
    let saved = initial.emitter(10).unwrap().identity();
    let mut excessive = valid.clone();
    excessive.area_m2 = 6.0001;
    assert!(matches!(
        NativeCoatingExchange::try_new(&native, &[excessive], budget(), || Ok(())),
        Err(NativeCoatingExchangeError::Invalid(_))
    ));
    assert!(
        NativeCoatingExchange::try_new(&native, &[valid.clone(), valid.clone()], budget(), || Ok(
            ()
        ))
        .is_err()
    );
    let mut limits = budget();
    limits.max_optical_evaluations = 47;
    assert!(matches!(
        NativeCoatingExchange::try_new(&native, &[valid.clone()], limits, || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    limits = budget();
    limits.max_source_evaluations = 1;
    let limited =
        NativeCoatingExchange::try_new(&native, &[valid.clone()], limits, || Ok(())).unwrap();
    assert!(matches!(
        native.advance_observed(&limited, 0.01, config(0.01), || Ok(())),
        Err(NativeCoatingExchangeError::Budget)
    ));
    let mut calls = 0;
    assert!(matches!(
        native.advance_observed(&initial, 0.1, config(0.025), || {
            calls += 1;
            if calls == 30 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(NativeCoatingExchangeError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(native.cells()[&CellId(11)], body);
    assert_eq!(initial.emitter(10).unwrap().identity(), saved);
    let retry = native
        .advance_observed(&initial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    let repeat = native
        .advance_observed(&initial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    assert_eq!(retry.0.energy_j().to_bits(), repeat.0.energy_j().to_bits());
    assert_eq!(retry.1.energy(), repeat.1.energy());

    let mut other = specimen(3100.0, 1.0);
    other.id = CellId(12);
    let contacted = thermal(
        vec![body, other],
        &HOT,
        vec![TimedEvent {
            time_s: 0.001,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CellId(11),
                b: CellId(12),
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        }],
    );
    let exposed = source(
        5.5,
        CoatingExchangeBath::Blackbody {
            temperature_k: 3100.0,
        },
        false,
    );
    let radiating = exchange(&contacted, &[exposed]);
    assert!(matches!(
        contacted.advance_observed(&radiating, 0.001, config(1.0), || Ok(())),
        Err(NativeCoatingExchangeError::Invalid(_))
    ));
    assert!(contacted.links().is_empty());
    let mut bounded = valid;
    bounded.receiver.temperature_bounds_k = (3100.0, 3100.0);
    let at_limit = exchange(&native, &[bounded]);
    assert!(matches!(
        native.advance_observed(&at_limit, 0.01, config(0.01), || Ok(())),
        Err(NativeCoatingExchangeError::Binding(_))
    ));
}

fn planck(t: f64, wavelength_nm: f64) -> f64 {
    let wavelength = wavelength_nm * 1e-9;
    let (h, c, k) = (6.62607015e-34, 299792458.0, 1.380649e-23);
    2.0 * h * c * c / wavelength.powi(5) / (h * c / (wavelength * k * t)).exp_m1() * 1e-9
}
// Independent original coherent owner + direct SI Planck, at the declared rule.
fn outward(t: f64) -> f64 {
    let optics = temperature_film(t);
    let layers = [CoatingLayer::Absorbing(&optics)];
    let stack = CoherentCoatingStack::try_new(&layers).unwrap();
    let mut total = 0.0;
    for w in 0..8 {
        let lambda = 380.0 + (w as f64 + 0.5) * 50.0;
        for a in 0..6 {
            let mu = (a as f64 + 0.5) / 6.0;
            let work = stack
                .layer_power_balance(&metal(), lambda, 1.0, mu)
                .unwrap();
            total += std::f64::consts::PI * 0.001 * 50.0 / 6.0
                * mu
                * (work.total.absorbed[0] + work.total.absorbed[1])
                * (planck(t, lambda) - planck(1700.0, lambda));
        }
    }
    total
}

#[test]
fn temperature_feedback_with_native_convection_and_reservoir_converges_to_independent_rk4() {
    let mut body = specimen(3100.0, 0.01);
    body.boundary.ambient_k = 1700.0;
    body.boundary.htc_w_m2_k = 0.02;
    body.boundary.reservoir = Some(Reservoir {
        area_m2: 1.0,
        conductance_w_k: 0.02,
        temperature_k: 1800.0,
    });
    body.boundary.source_w = 5.0;
    let native = thermal(vec![body], &HOT, vec![]);
    let source = source(0.001, CoatingExchangeBath::NativeAmbient, true);
    let initial = exchange(&native, &[source]);
    near(-initial.sources_w()[&CellId(11)], outward(3100.0), 5e-12);
    let derivative =
        |t: f64| (5.0 - outward(t) - 0.02 * 5.0 * (t - 1700.0) - 0.02 * (t - 1800.0)) / 0.1;
    let duration = 0.15;
    let dt = duration / 1024.0;
    let mut reference = 3100.0;
    let mut emitted = 0.0;
    for _ in 0..1024 {
        let k1 = derivative(reference);
        let p1 = outward(reference);
        let k2 = derivative(reference + dt * k1 / 2.0);
        let p2 = outward(reference + dt * k1 / 2.0);
        let k3 = derivative(reference + dt * k2 / 2.0);
        let p3 = outward(reference + dt * k2 / 2.0);
        let k4 = derivative(reference + dt * k3);
        let p4 = outward(reference + dt * k3);
        reference += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        emitted += dt * (p1 + 2.0 * p2 + 2.0 * p3 + p4) / 6.0;
    }
    let mut errors = Vec::new();
    let mut energy_errors = Vec::new();
    for dt in [0.025, 0.0125, 0.00625] {
        let (advanced, observed, audit) = native
            .advance_observed(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        errors.push((advanced.temperature_k(CellId(11)).unwrap() - reference).abs());
        energy_errors.push((observed.energy().outward_j - emitted).abs());
        near(
            advanced.energy_j() - native.energy_j(),
            audit.external_heat_j,
            1e-10,
        );
        assert!(observed.work().optical_evaluations > 48);
        near(observed.energy().balance_residual_j, 0.0, 1e-12);
    }
    assert!(
        errors[1] < 0.8 * errors[0] && errors[2] < 0.8 * errors[1],
        "{errors:?}"
    );
    assert!(
        energy_errors[1] < 0.8 * energy_errors[0] && energy_errors[2] < 0.8 * energy_errors[1],
        "{energy_errors:?}"
    );
}

#[test]
fn illumination_and_coating_cooling_balance_on_the_same_native_interval() {
    use fs_render::thermal_surface::heating::{
        SurfaceHeatingBudget, SurfaceHeatingTolerance, SurfaceIlluminationSource,
        SurfaceOpticalHeating, SurfacePowerResponse,
    };
    use fs_render::volumes::heating::illumination::schedule::{
        SpectralPowerSchedule, SpectralPowerSegment,
    };
    use fs_thermal::observers::CoupledObservers;
    let native = thermal(vec![specimen(2500.0, 1.0)], &HOT, vec![]);
    let mut emitting = source(
        0.001,
        CoatingExchangeBath::Blackbody {
            temperature_k: 1800.0,
        },
        false,
    );
    emitting.stack = IsothermalCoatingStack::Opaque {
        films: vec![],
        substrate: SubstrateOpticsModel::Fixed(metal()),
    };
    let cooling = exchange(&native, &[emitting]);
    let loss = -cooling.sources_w()[&CellId(11)];
    assert!(loss > 0.0);
    let response =
        SurfacePowerResponse::from_conductor([580.0], &metal(), [1.0], 1.0, [0.5]).unwrap();
    let incident = loss / response.absorbed_fraction()[0];
    let heating = SurfaceOpticalHeating::try_new(
        &[SurfaceIlluminationSource {
            source_id: 90,
            thermal_cell: CellId(11),
            material_id: MaterialId(7),
            temperature_bounds_k: (100.0, 5000.0),
            response,
            power: SpectralPowerSchedule::try_new(
                &[SpectralPowerSegment {
                    start_s: 0.0,
                    end_s: 1.0,
                    start_power_w: [incident],
                    end_power_w: [incident],
                }],
                1,
            )
            .unwrap(),
        }],
        &native,
        SurfaceHeatingBudget {
            max_sources: 1,
            max_total_power_segments: 1,
            max_source_evaluations: 128,
        },
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
    .unwrap();
    let pair = CoupledObservers::new(heating, cooling);
    let (after, accepted, audit) = native
        .advance_observed(&pair, 0.1, config(0.025), || Ok(()))
        .unwrap();
    near(after.temperature_k(CellId(11)).unwrap(), 2500.0, 1e-14);
    near(after.energy_j(), native.energy_j(), 1e-14);
    near(
        accepted.first.energy().thermal_input_j,
        accepted.second.energy().outward_j,
        1e-12,
    );
    near(audit.external_heat_j, 0.0, 1e-12);
    assert_eq!(accepted.first.time_s(), accepted.second.time_s());
}
