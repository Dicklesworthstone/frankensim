//! Surface R/T/A becomes heat exactly once on the native enthalpy timeline.

use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::dielectric::CauchyIor;
use fs_render::thermal_surface::heating::{
    SurfaceHeatingBudget, SurfaceHeatingError, SurfaceHeatingTolerance, SurfaceIlluminationSource,
    SurfaceOpticalHeating, SurfacePowerResponse,
};
use fs_render::thermal_surface::{
    GrayDiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity,
    ThermalSurfaceError,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_thermal::enthalpy::ThermalLink;

const KNOTS: [f64; 9] = [
    380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0,
];

fn optical_source(identity: u8) -> ConductorSource {
    ConductorSource::try_new(
        ContentHash([identity; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap()
}

fn metal(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        KNOTS.map(|wavelength| ConductorIorSample::try_new(wavelength, eta, k).unwrap()),
        optical_source(171),
    )
    .unwrap()
}

fn segment<const N: usize>(
    start_s: f64,
    end_s: f64,
    start_power_w: [f64; N],
    end_power_w: [f64; N],
) -> SpectralPowerSegment<N> {
    SpectralPowerSegment {
        start_s,
        end_s,
        start_power_w,
        end_power_w,
    }
}

fn source<const N: usize>(
    source_id: u64,
    thermal_cell: CellId,
    response: SurfacePowerResponse<N>,
    segments: &[SpectralPowerSegment<N>],
) -> SurfaceIlluminationSource<N> {
    SurfaceIlluminationSource {
        source_id,
        thermal_cell,
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
        response,
        power: SpectralPowerSchedule::try_new(segments, 16).unwrap(),
    }
}

fn limits() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 8,
        max_total_power_segments: 32,
        max_source_evaluations: 4096,
    }
}

fn heating<const N: usize>(
    thermal: &DepositionSimulation<'_>,
    sources: &[SurfaceIlluminationSource<N>],
) -> SurfaceOpticalHeating<N> {
    SurfaceOpticalHeating::try_new(
        sources,
        thermal,
        limits(),
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn clear_interface_transmits_power_without_heating_the_receiver() {
    let response = SurfacePowerResponse::from_dielectric(
        [480.0, 680.0],
        [1.0; 2],
        [1.5; 2],
        1.0,
        [0.0, 1.0],
        ContentHash([172; 32]),
    )
    .unwrap();
    let brewster = SurfacePowerResponse::from_dielectric(
        [580.0],
        [1.0],
        [1.5],
        1.0 / 3.25_f64.sqrt(),
        [0.0],
        ContentHash([173; 32]),
    )
    .unwrap();
    assert!(brewster.reflected_fraction()[0] < 1e-28);
    assert_eq!(brewster.absorbed_fraction(), &[0.0]);
    let original = thermal(vec![cell(11, 2.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let input = heating(
        &original,
        &[source(
            1,
            CellId(11),
            response,
            &[segment(0.25, 1.25, [8.0, 12.0], [8.0, 12.0])],
        )],
    );
    let (clear, accepted, audit) = original
        .advance_observed(&input, 2.0, config(3.0), || Ok(()))
        .unwrap();
    for (lane, fluence) in [8.0, 12.0].into_iter().enumerate() {
        close(accepted.energy().incident_j[lane], fluence);
        close(accepted.energy().reflected_j[lane], 0.04 * fluence);
        close(accepted.energy().transmitted_j[lane], 0.96 * fluence);
        assert_eq!(accepted.energy().absorbed_j[lane], 0.0);
    }
    assert_eq!(clear.cells(), original.cells());
    assert_eq!(accepted.energy().thermal_input_j, 0.0);
    assert_eq!(audit.external_heat_j, 0.0);
}

#[test]
fn gray_and_spectral_kirchhoff_absorption_retain_their_distinct_source_support() {
    let gray = GrayDiffuseThermalSurface::try_new(350.0, 0.25, ContentHash([191; 32])).unwrap();
    let emissivity = SpectralEmissivity::table(&[(380.0, 0.0), (780.0, 1.0)]).unwrap();
    let spectral_surface =
        SpectralDiffuseThermalSurface::try_new(350.0, emissivity, ContentHash([192; 32])).unwrap();
    let blackbody = GrayDiffuseThermalSurface::try_new(350.0, 1.0, ContentHash([193; 32])).unwrap();
    let wavelengths = [480.0, 680.0];
    let gray_response = SurfacePowerResponse::from_gray_diffuse(wavelengths, &gray).unwrap();
    let spectral_response = SurfacePowerResponse::from_spectral_diffuse(
        wavelengths,
        &emissivity,
        ContentHash([192; 32]),
    )
    .unwrap();
    for (lane, wavelength) in wavelengths.into_iter().enumerate() {
        let black = blackbody.radiance_w_m2_sr_nm(wavelength).unwrap();
        close(
            gray_response.absorbed_fraction()[lane],
            gray.radiance_w_m2_sr_nm(wavelength).unwrap() / black,
        );
        close(
            spectral_response.absorbed_fraction()[lane],
            spectral_surface.radiance_w_m2_sr_nm(wavelength).unwrap() / black,
        );
    }
    let original = thermal(
        vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let pulse = [segment(0.0, 1.0, [10.0, 20.0], [10.0, 20.0])];
    let input = heating(
        &original,
        &[
            source(1, CellId(11), gray_response, &pulse),
            source(2, CellId(12), spectral_response, &pulse),
        ],
    );
    let (hot, accepted, audit) = original
        .advance_observed(&input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    // The original emission law fixes absorptance via Kirchhoff. Its emitted
    // radiance is an observation, not an automatic second thermal source.
    close(
        hot.cells()[&CellId(11)].specific_h_j_kg,
        0.25 * (10.0 + 20.0),
    );
    close(
        hot.cells()[&CellId(12)].specific_h_j_kg,
        0.25 * 10.0 + 0.75 * 20.0,
    );
    close(accepted.energy().thermal_input_j, 25.0);
    close(audit.external_heat_j, 25.0);
    assert_eq!(accepted.energy().transmitted_j, [0.0; 2]);
    assert!(matches!(
        SurfacePowerResponse::from_spectral_diffuse(
            [10_000.0],
            &emissivity,
            ContentHash([192; 32])
        ),
        Err(SurfaceHeatingError::Surface(
            ThermalSurfaceError::Wavelength
        ))
    ));
    let infrared = SurfacePowerResponse::from_gray_diffuse([10_000.0], &gray).unwrap();
    let infrared_input = heating(
        &original,
        &[source(
            3,
            CellId(11),
            infrared,
            &[segment(0.0, 1.0, [8.0], [8.0])],
        )],
    );
    let (infrared_hot, infrared_accepted, _) = original
        .advance_observed(&infrared_input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    close(infrared_hot.energy_j(), 2.0);
    close(infrared_accepted.energy().absorbed_j[0], 2.0);
}

#[test]
fn lossless_stacks_add_no_heat_and_unresolved_weak_film_loss_refuses() {
    let films = [
        CoherentDielectricCoating::try_new(
            CauchyIor::try_constant(1.4).unwrap(),
            110e-9,
            optical_source(194),
        )
        .unwrap(),
        CoherentDielectricCoating::try_new(
            CauchyIor::try_constant(1.8).unwrap(),
            70e-9,
            optical_source(195),
        )
        .unwrap(),
    ];
    let response = SurfacePowerResponse::from_dielectric_stack(
        [480.0, 680.0],
        &[
            CoatingLayer::Dielectric(&films[0]),
            CoatingLayer::Dielectric(&films[1]),
        ],
        [1.0; 2],
        [1.5; 2],
        0.65,
        false,
        [0.2, 0.9],
        ContentHash([196; 32]),
    )
    .unwrap();
    assert_eq!(response.absorbed_fraction(), &[0.0; 2]);
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let input = SurfaceOpticalHeating::try_new(
        &[source(
            1,
            CellId(11),
            response,
            &[segment(0.0, 0.5, [1e6; 2], [1e6; 2])],
        )],
        &original,
        limits(),
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-7,
        },
        || Ok(()),
    )
    .unwrap();
    let (clear, accepted, audit) = original
        .advance_observed(&input, 0.5, config(0.25), || Ok(()))
        .unwrap();
    assert_eq!(accepted.energy().incident_j, [500_000.0; 2]);
    assert_eq!(accepted.energy().absorbed_j, [0.0; 2]);
    assert_eq!(clear.cells(), original.cells());
    assert_eq!(audit.external_heat_j, 0.0);
    for lane in 0..2 {
        let outgoing = accepted.energy().reflected_j[lane] + accepted.energy().transmitted_j[lane];
        assert!((outgoing - 500_000.0).abs() < 1e-7);
        assert!(accepted.energy().reflected_j[lane] > 0.0);
        assert!(accepted.energy().transmitted_j[lane] > 0.0);
    }
    let weak = CoherentAbsorbingCoating::try_new(
        KNOTS.map(|wavelength| CoatingIorSample::try_new(wavelength, 1.8, 1e-20).unwrap()),
        100e-9,
        optical_source(197),
    )
    .unwrap();
    assert!(matches!(
        SurfacePowerResponse::from_dielectric_stack(
            [580.0],
            &[CoatingLayer::Absorbing(&weak)],
            [1.0],
            [1.5],
            1.0,
            false,
            [0.5],
            ContentHash([198; 32]),
        ),
        Err(SurfaceHeatingError::Invalid(_))
    ));
}

#[test]
fn spectral_metal_pulses_store_analytic_absorption_in_native_latent_enthalpy() {
    let indices = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
    let optics = ConductorOptics::try_new(
        core::array::from_fn(|i| ConductorIorSample::try_new(KNOTS[i], 1.0, indices[i]).unwrap()),
        optical_source(174),
    )
    .unwrap();
    let response =
        SurfacePowerResponse::from_conductor([480.0, 680.0], &optics, [1.0; 2], 1.0, [0.0, 1.0])
            .unwrap();
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &LATENT, vec![]);
    let input = heating(
        &original,
        &[source(
            3,
            CellId(11),
            response,
            &[
                segment(0.125, 0.25, [80.0, 0.0], [80.0, 0.0]),
                segment(0.5, 0.75, [0.0; 2], [80.0, 120.0]),
                segment(0.75, 1.0, [80.0, 120.0], [0.0; 2]),
            ],
        )],
    );
    // Even a requested step longer than the source program must retain its
    // exact box and triangle areas. At normal incidence A=4n/((1+n)^2+k^2).
    let (hot, accepted, audit) = original
        .advance_observed(&input, 1.5, config(2.0), || Ok(()))
        .unwrap();
    let absorption = [4.0 / 5.0, 4.0 / 13.0];
    let mut stored = 0.0;
    for (lane, a) in absorption.into_iter().enumerate() {
        let joules = 30.0 * a;
        close(accepted.energy().incident_j[lane], 30.0);
        close(accepted.energy().absorbed_j[lane], joules);
        close(accepted.energy().reflected_j[lane], 30.0 - joules);
        assert_eq!(accepted.energy().transmitted_j[lane], 0.0);
        stored += joules;
    }
    close(hot.cells()[&CellId(11)].specific_h_j_kg, stored);
    close(hot.temperature_k(CellId(11)).unwrap(), 310.0);
    close(audit.external_heat_j, stored);
    close(accepted.energy().thermal_input_j, stored);
    assert!(stored > 10.0 && stored < 110.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(original.energy_j(), 0.0);
}

#[test]
fn absorbing_film_on_glass_heats_only_from_absorption_not_transmission() {
    let film = CoherentAbsorbingCoating::try_new(
        KNOTS.map(|wavelength| CoatingIorSample::try_new(wavelength, 1.8, 0.2).unwrap()),
        100e-9,
        optical_source(175),
    )
    .unwrap();
    let response = SurfacePowerResponse::from_dielectric_stack(
        [480.0, 680.0],
        &[CoatingLayer::Absorbing(&film)],
        [1.0; 2],
        [1.5; 2],
        1.0,
        false,
        [0.5; 2],
        ContentHash([176; 32]),
    )
    .unwrap();
    // Independent normal-incidence Airy amplitudes for n_f=1.8+0.2i,
    // d=100 nm, n_i=1, n_t=1.5: R=|r|^2 and T=(n_t/n_i)|t|^2.
    let expected = [
        [
            0.079_229_034_287_110_24,
            0.554_815_472_007_096_1,
            0.365_955_493_705_793_6,
        ],
        [
            0.126_379_946_617_055_43,
            0.614_073_948_238_498_9,
            0.259_546_105_144_445_65,
        ],
    ];
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let input = heating(
        &original,
        &[source(
            1,
            CellId(11),
            response,
            &[segment(0.0, 1.0, [10.0; 2], [10.0; 2])],
        )],
    );
    let (hot, accepted, audit) = original
        .advance_observed(&input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    let mut absorbed = 0.0;
    for (lane, [r, t, a]) in expected.into_iter().enumerate() {
        close(response.reflected_fraction()[lane], r);
        close(response.transmitted_fraction()[lane], t);
        close(response.absorbed_fraction()[lane], a);
        close(accepted.energy().reflected_j[lane], 10.0 * r);
        close(accepted.energy().transmitted_j[lane], 10.0 * t);
        close(accepted.energy().absorbed_j[lane], 10.0 * a);
        absorbed += 10.0 * a;
    }
    close(hot.energy_j(), absorbed);
    close(audit.external_heat_j, absorbed);
    assert!(accepted.energy().transmitted_j.iter().sum::<f64>() > absorbed);
}

#[test]
fn coated_polarized_sources_aggregate_with_conduction_and_dark_native_cooling() {
    let coating = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.4).unwrap(),
        100e-9,
        optical_source(177),
    )
    .unwrap();
    let optics = metal(1.0, 2.0).with_coating(coating);
    let s = SurfacePowerResponse::from_conductor([480.0, 680.0], &optics, [1.0; 2], 0.5, [1.0; 2])
        .unwrap();
    let p = SurfacePowerResponse::from_conductor([480.0, 680.0], &optics, [1.0; 2], 0.5, [0.0; 2])
        .unwrap();
    assert_ne!(s.absorbed_fraction(), p.absorbed_fraction());
    let clear = SurfacePowerResponse::from_dielectric(
        [480.0, 680.0],
        [1.0; 2],
        [1.5; 2],
        0.5,
        [0.5; 2],
        ContentHash([178; 32]),
    )
    .unwrap();
    let mut cooling = boundary(0.0);
    cooling.htc_w_m2_k = 0.1;
    let original = thermal(
        vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(11),
                    b: CellId(12),
                    area_m2: 1.0,
                    conductance_w_k: 0.5,
                }),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(11), cooling),
            },
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(12), cooling),
            },
        ],
    );
    let mut sources = vec![
        source(
            9,
            CellId(11),
            s,
            &[segment(0.0, 1.0, [4.0, 8.0], [4.0, 8.0])],
        ),
        source(
            1,
            CellId(11),
            p,
            &[segment(0.0, 1.0, [6.0, 2.0], [6.0, 2.0])],
        ),
        source(
            5,
            CellId(12),
            clear,
            &[segment(0.0, 1.0, [8.0, 6.0], [8.0, 6.0])],
        ),
    ];
    let input = heating(&original, &sources);
    sources.reverse();
    let reversed = heating(&original, &sources);
    assert_eq!(input.model_identity(), reversed.model_identity());
    let (hot, lit, audit) = original
        .advance_observed(&input, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let absorbed = s.absorbed_fraction()[0] * 4.0
        + s.absorbed_fraction()[1] * 8.0
        + p.absorbed_fraction()[0] * 6.0
        + p.absorbed_fraction()[1] * 2.0;
    close(hot.energy_j(), absorbed);
    close(audit.external_heat_j, absorbed);
    // The transparent receiver warms through the native contact, not through
    // its transmitted optical power. Intercepted watts need no area/cos factor.
    assert!(hot.temperature_k(CellId(12)).unwrap() > 300.0);
    assert!(hot.temperature_k(CellId(11)).unwrap() > hot.temperature_k(CellId(12)).unwrap());
    let (cold, dark, dark_audit) = hot
        .advance_observed(&lit, 3.0, config(0.125), || Ok(()))
        .unwrap();
    assert!(cold.energy_j() < hot.energy_j());
    assert!(dark_audit.external_heat_j < 0.0);
    close(
        dark_audit.stored_energy_change_j,
        dark_audit.external_heat_j,
    );
    assert_eq!(dark.energy().incident_j, lit.energy().incident_j);
    assert_eq!(dark.energy().thermal_input_j, lit.energy().thermal_input_j);
    let (whole, whole_input, _) = original
        .advance_observed(&reversed, 3.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(whole.cells(), cold.cells());
    assert_eq!(whole_input.energy(), dark.energy());
}

#[test]
fn weak_conductor_absorption_survives_reflectance_rounding_to_one() {
    let optics = metal(1e-20, 1.0);
    let response =
        SurfacePowerResponse::from_conductor([580.0], &optics, [1.0], 1.0, [0.5]).unwrap();
    assert_eq!(response.reflected_fraction(), &[1.0]);
    let a = response.absorbed_fraction()[0];
    assert!(a > 0.0 && (a / 2e-20 - 1.0).abs() < 1e-13);
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let input = SurfaceOpticalHeating::try_new(
        &[source(
            1,
            CellId(11),
            response,
            &[segment(0.0, 1.0, [1e20], [1e20])],
        )],
        &original,
        limits(),
        // Reflected watts are resolved only to their f64 spacing at 1e20 W;
        // direct absorption is separately representable and must reach heat.
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e5,
        },
        || Ok(()),
    )
    .unwrap();
    let (hot, accepted, audit) = original
        .advance_observed(&input, 1.0, config(0.25), || Ok(()))
        .unwrap();
    close(accepted.energy().absorbed_j[0], 2.0);
    close(accepted.energy().residual_j[0], -2.0);
    close(hot.energy_j(), 2.0);
    close(audit.external_heat_j, 2.0);
    close(hot.temperature_k(CellId(11)).unwrap(), 300.5);
}

#[test]
fn admission_bounds_budgets_and_late_cancellation_preserve_input_state() {
    let response =
        SurfacePowerResponse::from_conductor([580.0], &metal(1.0, 2.0), [1.0], 1.0, [0.5]).unwrap();
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let valid = source(1, CellId(11), response, &[segment(0.0, 2.0, [8.0], [8.0])]);
    let tolerance = SurfaceHeatingTolerance {
        power_roundoff_w: 1e-10,
    };
    let mut wrong_material = valid.clone();
    wrong_material.material_id = MaterialId(99);
    assert!(
        SurfaceOpticalHeating::try_new(
            &[wrong_material],
            &original,
            limits(),
            tolerance,
            || Ok(())
        )
        .is_err()
    );
    assert!(
        SurfaceOpticalHeating::try_new(
            &[valid.clone(), valid.clone()],
            &original,
            limits(),
            tolerance,
            || Ok(())
        )
        .is_err()
    );
    let mut exhausted = limits();
    exhausted.max_total_power_segments = 0;
    assert!(matches!(
        SurfaceOpticalHeating::try_new(
            &[valid.clone()],
            &original,
            exhausted,
            tolerance,
            || Ok(())
        ),
        Err(SurfaceHeatingError::Budget)
    ));
    exhausted = limits();
    exhausted.max_source_evaluations = 1;
    let bounded = SurfaceOpticalHeating::try_new(
        &[valid.clone()],
        &original,
        exhausted,
        tolerance,
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        original.advance_observed(&bounded, 2.0, config(0.5), || Ok(())),
        Err(SurfaceHeatingError::Budget)
    ));
    let mut limited = valid.clone();
    limited.temperature_bounds_k = (300.0, 300.75);
    let domain_limited = heating(&original, &[limited]);
    assert!(matches!(
        original.advance_observed(&domain_limited, 2.0, config(0.5), || Ok(())),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(11),
            ..
        })
    ));
    let input = heating(&original, &[valid]);
    let mut checkpoints = 0usize;
    original
        .advance_observed(&input, 2.0, config(0.125), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert!(checkpoints > 8);
    let mut calls = 0usize;
    assert!(matches!(
        original.advance_observed(&input, 2.0, config(0.125), || {
            calls += 1;
            if calls == checkpoints - 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    // Positive absorption cannot silently disappear in a power product, even
    // when integrating that vanished power would yield representable energy.
    // Nor may a positive native interval silently round its heat to zero.
    for (emissivity, power_w, duration_s) in [(1e-200, 1e-200, 1e200), (1.0, 1e-200, 1e-200)] {
        let gray =
            GrayDiffuseThermalSurface::try_new(300.0, emissivity, ContentHash([199; 32])).unwrap();
        let response = SurfacePowerResponse::from_gray_diffuse([580.0], &gray).unwrap();
        let underflow = heating(
            &original,
            &[source(
                8,
                CellId(11),
                response,
                &[segment(0.0, duration_s, [power_w], [power_w])],
            )],
        );
        assert!(matches!(
            original.advance_observed(&underflow, duration_s, config(duration_s), || Ok(())),
            Err(SurfaceHeatingError::Invalid(_))
        ));
        assert_eq!(underflow.time_s(), 0.0);
        assert_eq!(underflow.energy().incident_j, [0.0]);
        assert_eq!(underflow.energy().absorbed_j, [0.0]);
        assert_eq!(underflow.energy().thermal_input_j, 0.0);
        assert_eq!(underflow.source_evaluations(), 0);
    }
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.energy_j(), 0.0);
    for unchanged in [&input, &bounded, &domain_limited] {
        assert_eq!(unchanged.time_s(), 0.0);
        assert_eq!(unchanged.energy().thermal_input_j, 0.0);
        assert_eq!(unchanged.energy().incident_j, [0.0]);
    }
}
