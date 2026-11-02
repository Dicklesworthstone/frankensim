use super::*;
use fs_render::dielectric::{DielectricIorSample, TabulatedIor};
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;
use fs_render::thermal_surface::heating::{SurfacePowerResponse, TransmittingSurfacePowerResponse};

fn ir_film(k: f64, gap: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        std::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, 1.5, k).unwrap()
        }),
        2e-6,
        gap,
        authority(),
    )
    .unwrap()
}
fn ir_metal(eta: f64, k: f64, gap: f64) -> ConductorOptics {
    ConductorOptics::try_new_in_band(
        std::array::from_fn(|i| {
            ConductorIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, eta, k).unwrap()
        }),
        gap,
        authority(),
    )
    .unwrap()
}
fn ir_vacuum() -> DielectricIor {
    TabulatedIor::try_new_in_band(
        std::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, 1.0).unwrap()
        }),
        2000.0,
    )
    .unwrap()
    .into()
}

#[test]
fn infrared_state_interpolation_preserves_source_band_gap_policy_and_packet_work() {
    let coating = SurfaceCoatingLaw::try_temperature(
        authority(),
        &[(250.0, ir_film(0.1, 2500.0)), (350.0, ir_film(0.5, 2500.0))],
        2,
    )
    .unwrap();
    let conductor = SurfaceConductorLaw::try_temperature(
        authority(),
        &[
            (250.0, ir_metal(1.4, 0.2, 2500.0)),
            (350.0, ir_metal(1.8, 0.6, 2500.0)),
        ],
        2,
    )
    .unwrap();
    let film = coating.at_temperature(300.0).unwrap();
    let metal = conductor.at_temperature(300.0).unwrap();
    assert_eq!(film.wavelength_bounds_nm(), (4000.0, 20_000.0));
    assert_eq!(metal.wavelength_bounds_nm(), (4000.0, 20_000.0));
    assert_eq!(film.max_knot_span_nm(), 2500.0);
    assert_eq!(metal.max_knot_span_nm(), 2500.0);
    near(film.complex_ior(10_000.0).unwrap().k(), 0.3, 1e-14);
    near(metal.complex_ior(10_000.0).unwrap().eta(), 1.6, 1e-14);
    assert!(film.complex_ior(550.0).is_err());
    assert!(
        SurfaceCoatingLaw::try_temperature(
            authority(),
            &[(250.0, ir_film(0.1, 2000.0)), (350.0, ir_film(0.5, 2500.0))],
            2
        )
        .is_err()
    );
    assert!(
        SurfaceConductorLaw::try_temperature(
            authority(),
            &[
                (250.0, ir_metal(1.4, 0.2, 2000.0)),
                (350.0, ir_metal(1.8, 0.6, 2500.0))
            ],
            2
        )
        .is_err()
    );
    let layers = [CoatingLayer::Absorbing(&film)];
    let source = ContentHash([239; 32]);
    for reverse in [false, true] {
        let response = TransmittingSurfacePowerResponse::from_dielectric_stack(
            [10_000.0],
            &layers,
            [1.0],
            [1.0],
            0.7,
            reverse,
            [0.4],
            source,
        )
        .unwrap();
        let total = SurfacePowerResponse::from_dielectric_stack(
            [10_000.0],
            &layers,
            [1.0],
            [1.0],
            0.7,
            reverse,
            [0.4],
            source,
        )
        .unwrap();
        assert_eq!(response.total(), &total);
        assert!(response.film_absorbed_fractions()[0][0] > 0.0);
        assert!(
            SurfacePowerResponse::from_dielectric_stack(
                [550.0],
                &layers,
                [1.0],
                [1.0],
                0.7,
                reverse,
                [0.4],
                source
            )
            .is_err()
        );
    }
}

#[test]
fn room_temperature_infrared_cooling_uses_native_storage_and_converges_to_si_planck_rk4() {
    let native = thermal(vec![specimen(300.0, 0.001)], &HOT, vec![]);
    let mut physical = source(
        0.001,
        CoatingExchangeBath::Blackbody {
            temperature_k: 250.0,
        },
        false,
    );
    physical.stack = IsothermalCoatingStack::Opaque {
        films: vec![],
        substrate: SubstrateOpticsModel::Fixed(ir_metal(1.0, 1e-20, 2000.0)),
    };
    physical.quadrature = CoatingExchangeQuadrature {
        wavelength_bins: 128,
        cosine_bins: 1,
    };
    let initial = exchange(&native, &[physical]);
    assert_eq!(
        initial.emitter(10).unwrap().wavelength_bounds_nm(),
        (4000.0, 20_000.0)
    );
    let outward = |t: f64| {
        (0..128)
            .map(|i| {
                let lambda = 4000.0 + (i as f64 + 0.5) * 125.0;
                std::f64::consts::PI * 0.001 * 125.0 * (planck(t, lambda) - planck(250.0, lambda))
            })
            .sum::<f64>()
    };
    near(-initial.sources_w()[&CellId(11)], outward(300.0), 3e-12);
    assert!(outward(300.0) > 0.1);
    let duration = 0.1;
    let step = duration / 256.0;
    let derivative = |t| -outward(t) / 0.01;
    let mut reference = 300.0;
    for _ in 0..256 {
        let a = derivative(reference);
        let b = derivative(reference + 0.5 * step * a);
        let c = derivative(reference + 0.5 * step * b);
        let d = derivative(reference + step * c);
        reference += step * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let mut errors = Vec::new();
    for dt in [0.025, 0.0125, 0.00625] {
        let (after, observed, audit) = native
            .advance_observed(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        let temperature = after.temperature_k(CellId(11)).unwrap();
        assert!(temperature < 299.0 && temperature > 250.0);
        errors.push((temperature - reference).abs());
        near(
            after.energy_j() - native.energy_j(),
            -observed.energy().outward_j,
            1e-11,
        );
        near(
            audit.external_heat_j,
            observed.energy().thermal_external_j,
            1e-11,
        );
        assert_eq!(
            observed.work().optical_evaluations,
            0,
            "fixed IR optics reuse the prepared response"
        );
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
}

struct InfraredLatent;
impl ThermalLaw for InfraredLatent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (10.0..=110.0).contains(&h) {
            Ok(310.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([240; 32])
    }
}

#[test]
fn infrared_latent_enthalpy_feedback_changes_optics_at_fixed_temperature() {
    let mut body = cell(11, 0.001, 1.0, 0.0);
    body.specific_h_j_kg = 70.0;
    let native = thermal(vec![body], &InfraredLatent, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([240; 32]),
            (310.0, 310.0),
            &[(10.0, ir_film(0.01, 2500.0)), (110.0, ir_film(0.6, 2500.0))],
            2,
        )
        .unwrap(),
    );
    let mut physical = source(0.0001, CoatingExchangeBath::Vacuum, false);
    physical.receiver.temperature_bounds_k = (310.0, 310.0);
    physical.stack = IsothermalCoatingStack::Transmitting {
        films: vec![CoatingOpticsModel::Law(law)],
        opposite_ior: ir_vacuum(),
        reverse: false,
    };
    let initial = exchange(&native, &[physical]);
    let identity = initial.emitter(10).unwrap().optical_identity();
    let (after, observed, audit) = native
        .advance_observed(&initial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    assert_eq!(after.temperature_k(CellId(11)).unwrap(), 310.0);
    assert!(after.cells()[&CellId(11)].specific_h_j_kg < 70.0);
    assert_ne!(observed.emitter(10).unwrap().optical_identity(), identity);
    assert_eq!(
        observed.emitter(10).unwrap().wavelength_bounds_nm(),
        (4000.0, 20_000.0)
    );
    assert!(observed.work().optical_evaluations >= 4 * 48);
    near(audit.external_heat_j, -observed.energy().outward_j, 1e-12);
    near(
        after.energy_j() - native.energy_j(),
        audit.external_heat_j,
        1e-12,
    );
}
