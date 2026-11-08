//! Native storage and observed light use one accepted spectral emissivity law.
use super::*;
use fs_blake3::ContentHash;
use fs_render::thermal_surface::SpectralEmissivity;
use fs_render::thermal_surface::diffuse_exchange::DiffuseExchangeQuadrature;
use fs_render::thermal_surface::heating::SurfaceHeatingReceiver;
use fs_render::thermal_surface::heating::diffuse_exchange::*;
use fs_render::thermal_surface::heating::diffuse_feedback::SurfaceDiffuseLaw;
use std::sync::Arc;

struct Hot(bool);
impl ThermalLaw for Hot {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if self.0 {
            if h.is_finite() && (10.0..=110.0).contains(&h) {
                Ok(2500.0)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        } else if h.is_finite() && (-9000.0..=100_000.0).contains(&h) {
            Ok(1000.0 + h / 10.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        if self.0 { 0.0 } else { 0.1 }
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.0.then_some([237; 32])
    }
}
static HOT: Hot = Hot(false);
static HOT_LATENT: Hot = Hot(true);
fn specimen(t: f64) -> DepositedCell {
    let mut body = cell(11, 1.0, 1.0, 0.0);
    body.specific_h_j_kg = 10.0 * (t - 1000.0);
    body
}
fn budget() -> NativeDiffuseExchangeBudget {
    NativeDiffuseExchangeBudget {
        max_sources: 8,
        max_total_law_rows: 16,
        max_total_spectral_nodes: 4096,
        max_source_evaluations: 10_000,
        max_spectral_evaluations: 1_000_000,
    }
}
fn source(bath: DiffuseExchangeBath) -> DiffuseExchangeSource {
    let epsilon = SpectralEmissivity::constant_in_band(0.5, (100.0, 100_000.0)).unwrap();
    DiffuseExchangeSource {
        source_id: 10,
        source_identity: ContentHash([238; 32]),
        receiver: SurfaceHeatingReceiver {
            thermal_cell: CellId(11),
            material_id: MaterialId(7),
            temperature_bounds_k: (100.0, 5000.0),
        },
        area_m2: 0.0001,
        law: Arc::new(
            SurfaceDiffuseLaw::try_temperature(
                ContentHash([239; 32]),
                &[(100.0, epsilon), (5000.0, epsilon)],
                2,
            )
            .unwrap(),
        ),
        bath,
        quadrature: DiffuseExchangeQuadrature {
            subdivisions_per_interval: 16,
        },
        maximum_omitted_exchange_w: Some(1.0),
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn heating_cooling_equilibrium_and_vacuum_share_native_energy_and_observed_surface() {
    for (body_t, bath) in [
        (3100.0, Some(1800.0)),
        (1800.0, Some(3100.0)),
        (2500.0, Some(2500.0)),
        (3100.0, None),
    ] {
        let native = thermal(vec![specimen(body_t)], &HOT, vec![]);
        let source = source(bath.map_or(DiffuseExchangeBath::Vacuum, |temperature_k| {
            DiffuseExchangeBath::Blackbody { temperature_k }
        }));
        let initial =
            NativeDiffuseExchange::try_new(&native, &[source.clone()], budget(), || Ok(()))
                .unwrap();
        let (advanced, observed, audit) = native
            .advance_observed(&initial, 0.02, config(0.005), || Ok(()))
            .unwrap();
        let t = advanced.temperature_k(CellId(11)).unwrap();
        assert!(t >= body_t.min(bath.unwrap_or(0.0)) && t <= body_t.max(bath.unwrap_or(0.0)));
        if bath == Some(body_t) {
            assert_eq!(t, body_t);
            assert_eq!(observed.energy().outward_j, 0.0);
        } else {
            assert_ne!(t, body_t);
        }
        near(
            advanced.energy_j() - native.energy_j(),
            observed.energy().thermal_external_j,
            2.0e-10,
        );
        near(
            audit.external_heat_j,
            observed.energy().thermal_external_j,
            2.0e-12,
        );
        near(
            observed.energy().outward_j,
            -observed.energy().thermal_external_j,
            2.0e-12,
        );
        assert_eq!(observed.work().response_preparations, 0);
        let emitted = observed.emitter(10).unwrap();
        assert_eq!(emitted.temperature_k(), t);
        assert_eq!(
            *emitted,
            source
                .law
                .resolve_surface(&advanced, CellId(11), MaterialId(7))
                .unwrap()
        );
        assert_eq!(native.time_s(), 0.0);
    }
}

#[test]
fn latent_storage_changes_spectral_cooling_and_light_at_exactly_constant_temperature() {
    let mut body = cell(11, 1.0, 1.0, 0.0);
    body.specific_h_j_kg = 60.0;
    let native = thermal(vec![body], &HOT_LATENT, vec![]);
    let mut source = source(DiffuseExchangeBath::Vacuum);
    source.law = Arc::new(
        SurfaceDiffuseLaw::try_enthalpy(
            ContentHash([240; 32]),
            ContentHash([237; 32]),
            (2500.0, 2500.0),
            &[
                (
                    10.0,
                    SpectralEmissivity::table_in_band(&[
                        (100.0, 0.1),
                        (1000.0, 0.4),
                        (100_000.0, 0.1),
                    ])
                    .unwrap(),
                ),
                (
                    110.0,
                    SpectralEmissivity::table_in_band(&[
                        (100.0, 0.8),
                        (1000.0, 0.9),
                        (100_000.0, 0.8),
                    ])
                    .unwrap(),
                ),
            ],
            2,
        )
        .unwrap(),
    );
    let initial =
        NativeDiffuseExchange::try_new(&native, &[source.clone()], budget(), || Ok(())).unwrap();
    let (advanced, observed, audit) = native
        .advance_observed(&initial, 0.02, config(0.005), || Ok(()))
        .unwrap();
    assert_eq!(advanced.temperature_k(CellId(11)).unwrap(), 2500.0);
    assert!(advanced.cells()[&CellId(11)].specific_h_j_kg < 60.0);
    assert!(
        observed.rates(10).unwrap().outward_power_w < initial.rates(10).unwrap().outward_power_w
    );
    assert!(
        observed
            .emitter(10)
            .unwrap()
            .radiance_w_m2_sr_nm(1000.0)
            .unwrap()
            < initial
                .emitter(10)
                .unwrap()
                .radiance_w_m2_sr_nm(1000.0)
                .unwrap()
    );
    assert_eq!(
        *observed.emitter(10).unwrap(),
        source
            .law
            .resolve_surface(&advanced, CellId(11), MaterialId(7))
            .unwrap()
    );
    near(
        advanced.energy_j() - native.energy_j(),
        audit.external_heat_j,
        2.0e-12,
    );
    near(audit.external_heat_j, -observed.energy().outward_j, 2.0e-12);
    assert_eq!(observed.work().response_preparations, 0);
}

#[test]
fn visible_only_sources_fail_coverage_and_budget_cancellation_leave_native_history_unchanged() {
    let native = thermal(vec![specimen(2500.0)], &HOT, vec![]);
    let mut visible = source(DiffuseExchangeBath::Vacuum);
    let epsilon = SpectralEmissivity::constant(0.5).unwrap();
    visible.law = Arc::new(
        SurfaceDiffuseLaw::try_temperature(
            ContentHash([241; 32]),
            &[(100.0, epsilon), (5000.0, epsilon)],
            2,
        )
        .unwrap(),
    );
    assert!(matches!(
        NativeDiffuseExchange::try_new(&native, &[visible.clone()], budget(), || Ok(())),
        Err(NativeDiffuseExchangeError::Coverage { .. })
    ));
    visible.maximum_omitted_exchange_w = None;
    let partial = NativeDiffuseExchange::try_new(&native, &[visible], budget(), || Ok(())).unwrap();
    assert!(partial.rates(10).unwrap().omitted_exchange_upper_bound_w > 200.0);
    let mut tight = budget();
    tight.max_source_evaluations = 1;
    let bounded = NativeDiffuseExchange::try_new(
        &native,
        &[source(DiffuseExchangeBath::Vacuum)],
        tight,
        || Ok(()),
    )
    .unwrap();
    assert!(
        native
            .advance_observed(&bounded, 0.02, config(0.005), || Ok(()))
            .is_err()
    );
    let mut calls = 0;
    assert!(
        native
            .advance_observed(&partial, 0.02, config(0.005), || {
                calls += 1;
                if calls > 20 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(partial.time_s(), 0.0);
    assert_eq!(bounded.energy().outward_j, 0.0);
    let mut covered = source(DiffuseExchangeBath::Vacuum);
    covered.area_m2 = 7.0;
    covered.maximum_omitted_exchange_w = None;
    assert!(NativeDiffuseExchange::try_new(&native, &[covered], budget(), || Ok(())).is_err());
}

#[test]
fn native_ambient_events_and_split_advance_keep_one_radiation_clock() {
    let mut body = specimen(2500.0);
    body.boundary.ambient_k = 1800.0;
    let mut later = body.boundary;
    later.ambient_k = 3100.0;
    let native = thermal(
        vec![body],
        &HOT,
        vec![TimedEvent {
            time_s: 0.01,
            event: Event::Boundary(CellId(11), later),
        }],
    );
    let source = source(DiffuseExchangeBath::NativeAmbient);
    let initial = NativeDiffuseExchange::try_new(&native, &[source], budget(), || Ok(())).unwrap();
    let (whole, observed, _) = native
        .advance_observed(&initial, 0.02, config(0.005), || Ok(()))
        .unwrap();
    let (middle, partial, _) = native
        .advance_observed(&initial, 0.01, config(0.005), || Ok(()))
        .unwrap();
    let (split, resumed, _) = middle
        .advance_observed(&partial, 0.02, config(0.005), || Ok(()))
        .unwrap();
    assert!(middle.temperature_k(CellId(11)).unwrap() < 2500.0);
    assert!(whole.temperature_k(CellId(11)).unwrap() > middle.temperature_k(CellId(11)).unwrap());
    near(whole.energy_j(), split.energy_j(), 2.0e-12);
    near(
        observed.energy().outward_j,
        resumed.energy().outward_j,
        2.0e-10,
    );
    assert_eq!(observed.emitter(10), resumed.emitter(10));
}
