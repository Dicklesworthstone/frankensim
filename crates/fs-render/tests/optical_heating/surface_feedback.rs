//! Constitutive surface optics follow the original native temperature/enthalpy.

use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorOptics,
    ConductorSource,
};
use fs_render::dielectric::CauchyIor;
use fs_render::thermal_surface::heating::diffuse_feedback::SurfaceDiffuseLaw;
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;
use fs_render::thermal_surface::SpectralEmissivity;
use fs_render::thermal_surface::heating::{
    StatefulDiffuseSurfaceIlluminationSource, StatefulSurfaceIlluminationSource,
    SurfaceHeatingBudget, SurfaceHeatingError,
    SurfaceHeatingTolerance, SurfaceIlluminationSource, SurfaceIncidence, SurfaceOpticalHeating,
    SurfacePowerResponse,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use std::sync::Arc;

const KNOTS: [f64; 9] = [
    380.0, 430.0, 480.0, 530.0, 580.0, 630.0, 680.0, 730.0, 780.0,
];

fn declared_source(identity: u8) -> ConductorSource {
    ConductorSource::try_new(
        ContentHash([identity; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap()
}

fn metal(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        KNOTS.map(|wavelength| ConductorIorSample::try_new(wavelength, eta, k).unwrap()),
        declared_source(201),
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
    id: u64,
    cell: CellId,
    wavelengths_nm: [f64; N],
    optics: Arc<SurfaceConductorLaw>,
    segments: &[SpectralPowerSegment<N>],
) -> StatefulSurfaceIlluminationSource<N> {
    StatefulSurfaceIlluminationSource {
        source_id: id,
        thermal_cell: cell,
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
        incidence: SurfaceIncidence {
            wavelengths_nm,
            incident_eta: [1.0; N],
            incidence_cosine: 1.0,
            s_power_fraction: [0.5; N],
        },
        optics,
        power: SpectralPowerSchedule::try_new(segments, 16).unwrap(),
    }
}

fn budget() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 8,
        max_total_power_segments: 32,
        max_source_evaluations: 4096,
    }
}

fn tolerance() -> SurfaceHeatingTolerance {
    SurfaceHeatingTolerance {
        power_roundoff_w: 1e-10,
    }
}

fn heating<const N: usize>(
    thermal: &DepositionSimulation<'_>,
    sources: &[StatefulSurfaceIlluminationSource<N>],
) -> SurfaceOpticalHeating<N> {
    SurfaceOpticalHeating::try_new_with_stateful(
        &[],
        sources,
        thermal,
        budget(),
        32,
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}

struct LatentDatum {
    shift: f64,
    identity: u8,
}
impl ThermalLaw for LatentDatum {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        LATENT.temperature_k(h + self.shift)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        LATENT.maximum_temperature_slope()
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([self.identity; 32])
    }
}
static IDENTIFIED_LATENT: LatentDatum = LatentDatum {
    shift: 0.0,
    identity: 202,
};
static NEGATIVE_DATUM: LatentDatum = LatentDatum {
    shift: 100.0,
    identity: 203,
};

#[test]
fn interpolate_complex_indices_before_fresnel_and_retain_the_fixed_coating() {
    let mut native_cell = cell(11, 1.0, 1.0, 0.0);
    native_cell.specific_h_j_kg = 4.0; // Native cp=4: T=301 K.
    let native = thermal(vec![native_cell], &SENSIBLE, vec![]);
    let law = SurfaceConductorLaw::try_temperature(
        declared_source(204),
        &[(300.0, metal(1.0, 1.0)), (302.0, metal(3.0, 3.0))],
        2,
    )
    .unwrap();
    let resolved = law.resolve(&native, CellId(11), MaterialId(7)).unwrap();
    let ior = resolved.complex_ior(580.0).unwrap();
    close(ior.eta(), 2.0);
    close(ior.k(), 2.0);
    let response =
        SurfacePowerResponse::from_conductor([580.0], &resolved, [1.0], 1.0, [0.5]).unwrap();
    close(response.absorbed_fraction()[0], 8.0 / 13.0);
    let blended_endpoint_absorption = 0.5 * (4.0 / 5.0 + 12.0 / 25.0);
    assert!((response.absorbed_fraction()[0] - blended_endpoint_absorption).abs() > 0.01);

    let film = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.4).unwrap(),
        100e-9,
        declared_source(205),
    )
    .unwrap();
    let coated = SurfaceConductorLaw::try_temperature(
        declared_source(206),
        &[
            (300.0, metal(1.0, 1.0).with_coating(film)),
            (302.0, metal(3.0, 3.0).with_coating(film)),
        ],
        2,
    )
    .unwrap();
    let actual = coated.resolve(&native, CellId(11), MaterialId(7)).unwrap();
    let expected = metal(2.0, 2.0).with_coating(film);
    for wavelength in [480.0, 680.0] {
        let actual = actual.power_balance(wavelength, 1.0, 0.6).unwrap();
        let expected = expected.power_balance(wavelength, 1.0, 0.6).unwrap();
        for polarization in 0..2 {
            close(
                actual.reflected[polarization],
                expected.reflected[polarization],
            );
            close(
                actual.absorbed[polarization],
                expected.absorbed[polarization],
            );
        }
    }
    assert!(
        SurfaceConductorLaw::try_temperature(
            declared_source(207),
            &[
                (300.0, metal(1.0, 1.0)),
                (302.0, metal(3.0, 3.0).with_coating(film))
            ],
            2,
        )
        .is_err()
    );
}

#[test]
fn temperature_feedback_converges_to_an_independent_separated_ode() {
    let law = Arc::new(
        SurfaceConductorLaw::try_temperature(
            declared_source(208),
            &[(300.0, metal(1.0, 1.0)), (310.0, metal(1.0, 11.0))],
            2,
        )
        .unwrap(),
    );
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let input = heating(
        &original,
        &[source(
            1,
            CellId(11),
            [580.0],
            law,
            &[segment(0.0, 2.0, [16.0], [16.0])],
        )],
    );
    // Normal incidence, n=1, k=1+x, x=T-300, native heat capacity=4 J/K:
    // dx/dt=16/(4+(1+x)^2), hence 16t=5x+x^2+x^3/3.
    // Invert that strictly increasing polynomial, not a second time integrator.
    let (mut lo, mut hi) = (0.0, 10.0);
    for _ in 0..80 {
        let x = 0.5 * (lo + hi);
        if 5.0 * x + x * x + x * x * x / 3.0 < 32.0 {
            lo = x;
        } else {
            hi = x;
        }
    }
    let exact_temperature = 300.0 + 0.5 * (lo + hi);
    let mut errors = Vec::new();
    for dt in [0.5, 0.25, 0.125] {
        let (hot, accepted, audit) = original
            .advance_observed(&input, 2.0, config(dt), || Ok(()))
            .unwrap();
        let temperature = hot.temperature_k(CellId(11)).unwrap();
        errors.push((temperature - exact_temperature).abs());
        assert!(accepted.response_evaluations() > 0);
        close(accepted.energy().incident_j[0], 32.0);
        close(
            accepted.energy().thermal_input_j,
            4.0 * (temperature - 300.0),
        );
        close(audit.external_heat_j, accepted.energy().thermal_input_j);
        close(
            accepted.energy().reflected_j[0] + accepted.energy().absorbed_j[0],
            32.0,
        );
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 0.08);
    assert_eq!(input.response_evaluations(), 0);
    assert_eq!(original.energy_j(), 0.0);
}

#[test]
fn equal_temperature_latent_states_have_different_optics_and_heat() {
    let law = Arc::new(
        SurfaceConductorLaw::try_enthalpy(
            declared_source(209),
            ContentHash([202; 32]),
            (300.0, 500.0),
            &[
                (0.0, metal(1.0, 1.0)),
                (10.0, metal(1.0, 1.0)),
                (110.0, metal(1.0, 3.0)),
                (200.0, metal(1.0, 3.0)),
            ],
            4,
        )
        .unwrap(),
    );
    let mut cells = vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)];
    cells[0].specific_h_j_kg = 30.0;
    cells[1].specific_h_j_kg = 90.0;
    let original = thermal(cells, &IDENTIFIED_LATENT, vec![]);
    assert_eq!(original.temperature_k(CellId(11)).unwrap(), 310.0);
    assert_eq!(original.temperature_k(CellId(12)).unwrap(), 310.0);
    let sources: Vec<_> = [CellId(11), CellId(12)]
        .into_iter()
        .map(|id| {
            source(
                id.0,
                id,
                [580.0],
                law.clone(),
                &[segment(0.0, 1.0, [5.0], [5.0])],
            )
        })
        .collect();
    let input = heating(&original, &sources);
    let absorption = |h: f64| {
        let k = 1.0 + 2.0 * (h - 10.0) / 100.0;
        4.0 / (4.0 + k * k)
    };
    let (first, accepted, _) = original
        .advance_observed(&input, 0.125, config(0.125), || Ok(()))
        .unwrap();
    for (id, initial) in [(CellId(11), 30.0), (CellId(12), 90.0)] {
        close(
            first.cells()[&id].specific_h_j_kg,
            initial + 0.625 * absorption(initial),
        );
        close(first.temperature_k(id).unwrap(), 310.0);
    }
    let (second, updated, audit) = first
        .advance_observed(&accepted, 0.25, config(0.125), || Ok(()))
        .unwrap();
    for id in [CellId(11), CellId(12)] {
        let h = first.cells()[&id].specific_h_j_kg;
        close(
            second.cells()[&id].specific_h_j_kg,
            h + 0.625 * absorption(h),
        );
        close(second.temperature_k(id).unwrap(), 310.0);
    }
    assert!(updated.response_evaluations() > 0);
    close(
        second.energy_j() - original.energy_j(),
        updated.energy().thermal_input_j,
    );
    close(
        audit.external_heat_j,
        updated.energy().thermal_input_j - accepted.energy().thermal_input_j,
    );
}

#[test]
fn negative_enthalpy_datums_are_bound_and_dark_domain_failures_are_transactional() {
    let law = Arc::new(
        SurfaceConductorLaw::try_enthalpy(
            declared_source(210),
            ContentHash([203; 32]),
            (300.0, 500.0),
            &[(-100.0, metal(1.0, 1.0)), (100.0, metal(1.0, 3.0))],
            2,
        )
        .unwrap(),
    );
    let mut native_cell = cell(11, 1.0, 1.0, 0.0);
    native_cell.specific_h_j_kg = -100.0;
    let original = thermal(vec![native_cell], &NEGATIVE_DATUM, vec![]);
    let input = heating(
        &original,
        &[source(
            1,
            CellId(11),
            [580.0],
            law,
            &[segment(0.0, 1.0, [8.0], [8.0])],
        )],
    );
    let (hot, accepted, audit) = original
        .advance_observed(&input, 1.0, config(1.0), || Ok(()))
        .unwrap();
    close(hot.cells()[&CellId(11)].specific_h_j_kg, -93.6);
    close(audit.stored_energy_change_j, 6.4);
    close(accepted.energy().thermal_input_j, 6.4);
    let foreign = SurfaceConductorLaw::try_enthalpy(
        declared_source(211),
        ContentHash([204; 32]),
        (300.0, 500.0),
        &[(-100.0, metal(1.0, 1.0)), (100.0, metal(1.0, 3.0))],
        2,
    )
    .unwrap();
    assert!(
        foreign
            .resolve(&original, CellId(11), MaterialId(7))
            .is_err()
    );

    let dark_law = Arc::new(
        SurfaceConductorLaw::try_temperature(
            declared_source(212),
            &[(300.0, metal(1.0, 1.0)), (300.5, metal(1.0, 2.0))],
            2,
        )
        .unwrap(),
    );
    let externally_heated = thermal(vec![cell(11, 1.0, 1.0, 8.0)], &SENSIBLE, vec![]);
    let dark = heating(
        &externally_heated,
        &[source(2, CellId(11), [580.0], dark_law, &[])],
    );
    assert!(
        externally_heated
            .advance_observed(&dark, 1.0, config(0.25), || Ok(()))
            .is_err()
    );
    assert_eq!(externally_heated.time_s(), 0.0);
    assert_eq!(externally_heated.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    assert_eq!(dark.time_s(), 0.0);
    assert_eq!(dark.energy().incident_j, [0.0]);
    assert_eq!(dark.energy().thermal_input_j, 0.0);
    assert_eq!(dark.response_evaluations(), 0);
}

#[test]
fn constant_constitutive_rows_preserve_fixed_packet_bits_and_skip_fresnel_rebuilds() {
    let optics = metal(1.0, 2.0);
    let law = Arc::new(
        SurfaceConductorLaw::try_temperature(
            declared_source(213),
            &[(300.0, optics), (600.0, optics)],
            2,
        )
        .unwrap(),
    );
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let stateful = source(
        1,
        CellId(11),
        [480.0, 680.0],
        law,
        &[
            segment(0.125, 0.5, [8.0, 4.0], [8.0, 4.0]),
            segment(0.5, 0.75, [8.0, 4.0], [0.0; 2]),
        ],
    );
    let fixed = SurfaceIlluminationSource {
        source_id: stateful.source_id,
        thermal_cell: stateful.thermal_cell,
        material_id: stateful.material_id,
        temperature_bounds_k: stateful.temperature_bounds_k,
        response: SurfacePowerResponse::from_conductor(
            [480.0, 680.0],
            &optics,
            [1.0; 2],
            1.0,
            [0.5; 2],
        )
        .unwrap(),
        power: stateful.power.clone(),
    };
    let frozen =
        SurfaceOpticalHeating::try_new(&[fixed], &original, budget(), tolerance(), || Ok(()))
            .unwrap();
    let dynamic = heating(&original, &[stateful.clone()]);
    let (frozen_hot, frozen_account, _) = original
        .advance_observed(&frozen, 1.5, config(0.125), || Ok(()))
        .unwrap();
    let (dynamic_hot, dynamic_account, _) = original
        .advance_observed(&dynamic, 1.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(frozen_hot.cells(), dynamic_hot.cells());
    assert_eq!(frozen_account.energy(), dynamic_account.energy());
    assert_eq!(dynamic_account.response_evaluations(), 0);
    close(dynamic_account.energy().thermal_input_j, 3.0);
    assert!(matches!(
        SurfaceOpticalHeating::try_new_with_stateful(
            &[],
            &[stateful],
            &original,
            budget(),
            1,
            tolerance(),
            || Ok(())
        ),
        Err(SurfaceHeatingError::Budget)
    ));
}

#[test]
fn shared_stateful_sources_are_canonical_across_split_advances_and_cancellation() {
    let law = Arc::new(
        SurfaceConductorLaw::try_temperature(
            declared_source(214),
            &[(300.0, metal(1.0, 1.0)), (400.0, metal(1.0, 11.0))],
            2,
        )
        .unwrap(),
    );
    let original = thermal(
        vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let mut sources = vec![
        source(
            3,
            CellId(11),
            [580.0],
            law.clone(),
            &[segment(0.0, 1.5, [3.0], [3.0])],
        ),
        source(
            1,
            CellId(11),
            [580.0],
            law.clone(),
            &[segment(0.25, 1.25, [4.0], [0.0])],
        ),
        source(
            2,
            CellId(12),
            [580.0],
            law,
            &[segment(0.5, 1.5, [6.0], [6.0])],
        ),
    ];
    let fixed = SurfaceIlluminationSource {
        source_id: 4,
        thermal_cell: CellId(12),
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
        response: SurfacePowerResponse::from_conductor(
            [580.0],
            &metal(1.0, 2.0),
            [1.0],
            1.0,
            [0.5],
        )
        .unwrap(),
        power: SpectralPowerSchedule::try_new(&[segment(0.125, 0.75, [2.0], [2.0])], 1).unwrap(),
    };
    // Three source bindings share just two native-state rows in this budget.
    let input = SurfaceOpticalHeating::try_new_with_stateful(
        &[fixed.clone()],
        &sources,
        &original,
        budget(),
        2,
        tolerance(),
        || Ok(()),
    )
    .unwrap();
    sources.reverse();
    let reversed = SurfaceOpticalHeating::try_new_with_stateful(
        &[fixed.clone()],
        &sources,
        &original,
        budget(),
        2,
        tolerance(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(input.model_identity(), reversed.model_identity());
    let (whole, complete, _) = original
        .advance_observed(&input, 2.0, config(0.125), || Ok(()))
        .unwrap();
    let (middle, partial, _) = original
        .advance_observed(&reversed, 0.75, config(0.125), || Ok(()))
        .unwrap();
    let (split, continued, _) = middle
        .advance_observed(&partial, 2.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(whole.cells(), split.cells());
    assert_eq!(complete.energy(), continued.energy());
    close(complete.energy().incident_j[0], 13.75);
    close(whole.energy_j(), complete.energy().thermal_input_j);
    let mut collision = fixed;
    collision.source_id = sources[0].source_id;
    assert!(
        SurfaceOpticalHeating::try_new_with_stateful(
            &[collision],
            &sources,
            &original,
            budget(),
            2,
            tolerance(),
            || Ok(())
        )
        .is_err()
    );
    let mut checkpoints = 0usize;
    original
        .advance_observed(&input, 2.0, config(0.125), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0usize;
    assert!(
        original
            .advance_observed(&input, 2.0, config(0.125), || {
                calls += 1;
                if calls == checkpoints - 2 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(original.energy_j(), 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(input.response_evaluations(), 0);
}


#[test]
fn diffuse_enthalpy_feedback_changes_absorption_inside_one_isothermal_melt_plateau() {
    let law = Arc::new(
        SurfaceDiffuseLaw::try_enthalpy(
            ContentHash([215; 32]),
            ContentHash([202; 32]),
            (300.0, 500.0),
            &[
                (0.0, SpectralEmissivity::constant(0.2).unwrap()),
                (10.0, SpectralEmissivity::constant(0.2).unwrap()),
                (110.0, SpectralEmissivity::constant(0.8).unwrap()),
                (200.0, SpectralEmissivity::constant(0.8).unwrap()),
            ],
            4,
        )
        .unwrap(),
    );
    let mut cells = vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)];
    cells[0].specific_h_j_kg = 30.0;
    cells[1].specific_h_j_kg = 90.0;
    let original = thermal(cells, &IDENTIFIED_LATENT, vec![]);
    assert_eq!(original.temperature_k(CellId(11)).unwrap(), 310.0);
    assert_eq!(original.temperature_k(CellId(12)).unwrap(), 310.0);

    let sources: Vec<_> = [CellId(11), CellId(12)]
        .into_iter()
        .map(|id| StatefulDiffuseSurfaceIlluminationSource {
            source_id: id.0,
            thermal_cell: id,
            material_id: MaterialId(7),
            temperature_bounds_k: (300.0, 500.0),
            wavelengths_nm: [580.0],
            optics: Arc::clone(&law),
            power: SpectralPowerSchedule::try_new(
                &[segment(0.0, 0.25, [5.0], [5.0])],
                1,
            )
            .unwrap(),
        })
        .collect();
    let input = SurfaceOpticalHeating::try_new_with_stateful_diffuse(
        &[],
        &sources,
        &original,
        budget(),
        4,
        tolerance(),
        || Ok(()),
    )
    .unwrap();

    let emissivity = |h: f64| 0.2 + 0.6 * (h - 10.0) / 100.0;
    let (first, accepted, _) = original
        .advance_observed(&input, 0.125, config(0.125), || Ok(()))
        .unwrap();
    for (id, initial) in [(CellId(11), 30.0), (CellId(12), 90.0)] {
        close(
            first.cells()[&id].specific_h_j_kg,
            initial + 0.625 * emissivity(initial),
        );
        close(first.temperature_k(id).unwrap(), 310.0);
    }
    let (second, updated, audit) = first
        .advance_observed(&accepted, 0.25, config(0.125), || Ok(()))
        .unwrap();
    for id in [CellId(11), CellId(12)] {
        let h = first.cells()[&id].specific_h_j_kg;
        close(
            second.cells()[&id].specific_h_j_kg,
            h + 0.625 * emissivity(h),
        );
        close(second.temperature_k(id).unwrap(), 310.0);
    }
    assert!(updated.response_evaluations() > 0);
    close(
        second.energy_j() - original.energy_j(),
        updated.energy().thermal_input_j,
    );
    close(
        audit.external_heat_j,
        updated.energy().thermal_input_j - accepted.energy().thermal_input_j,
    );
    close(
        updated.energy().reflected_j[0] + updated.energy().absorbed_j[0],
        updated.energy().incident_j[0],
    );
}

#[test]
fn constant_diffuse_state_rows_match_the_fixed_spectral_response_and_cache_it() {
    let emissivity = SpectralEmissivity::table(&[
        (380.0, 0.25),
        (580.0, 0.6),
        (780.0, 0.4),
    ])
    .unwrap();
    let law = Arc::new(
        SurfaceDiffuseLaw::try_temperature(
            ContentHash([216; 32]),
            &[(300.0, emissivity), (600.0, emissivity)],
            2,
        )
        .unwrap(),
    );
    let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
    let stateful = StatefulDiffuseSurfaceIlluminationSource {
        source_id: 1,
        thermal_cell: CellId(11),
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
        wavelengths_nm: [480.0, 680.0],
        optics: law,
        power: SpectralPowerSchedule::try_new(
            &[segment(0.0, 1.0, [8.0, 4.0], [8.0, 4.0])],
            1,
        )
        .unwrap(),
    };
    let fixed = SurfaceIlluminationSource {
        source_id: 1,
        thermal_cell: CellId(11),
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
        response: SurfacePowerResponse::from_spectral_diffuse(
            [480.0, 680.0],
            &emissivity,
            ContentHash([216; 32]),
        )
        .unwrap(),
        power: stateful.power.clone(),
    };
    let frozen =
        SurfaceOpticalHeating::try_new(&[fixed], &original, budget(), tolerance(), || Ok(()))
            .unwrap();
    let dynamic = SurfaceOpticalHeating::try_new_with_stateful_diffuse(
        &[],
        &[stateful],
        &original,
        budget(),
        2,
        tolerance(),
        || Ok(()),
    )
    .unwrap();
    let (fixed_state, fixed_history, _) = original
        .advance_observed(&frozen, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let (dynamic_state, dynamic_history, _) = original
        .advance_observed(&dynamic, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(fixed_state.cells(), dynamic_state.cells());
    assert_eq!(fixed_history.energy(), dynamic_history.energy());
    assert_eq!(dynamic_history.response_evaluations(), 0);
}


#[test]
fn diffuse_state_law_uses_the_same_latent_emissivity_for_absorption_and_emission() {
    let law = SurfaceDiffuseLaw::try_enthalpy(
        ContentHash([217; 32]),
        ContentHash([202; 32]),
        (300.0, 500.0),
        &[
            (0.0, SpectralEmissivity::constant(0.2).unwrap()),
            (10.0, SpectralEmissivity::constant(0.2).unwrap()),
            (110.0, SpectralEmissivity::constant(0.8).unwrap()),
            (200.0, SpectralEmissivity::constant(0.8).unwrap()),
        ],
        4,
    )
    .unwrap();
    let mut cells = vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)];
    cells[0].specific_h_j_kg = 30.0;
    cells[1].specific_h_j_kg = 90.0;
    let native = thermal(cells, &IDENTIFIED_LATENT, vec![]);
    assert_eq!(native.temperature_k(CellId(11)).unwrap(), 310.0);
    assert_eq!(native.temperature_k(CellId(12)).unwrap(), 310.0);

    let low = law.resolve_surface(&native, CellId(11), MaterialId(7)).unwrap();
    let high = law.resolve_surface(&native, CellId(12), MaterialId(7)).unwrap();
    let low_e = 0.2 + 0.6 * 20.0 / 100.0;
    let high_e = 0.2 + 0.6 * 80.0 / 100.0;
    close(low.emissivity().at(580.0).unwrap(), low_e);
    close(high.emissivity().at(580.0).unwrap(), high_e);
    let blackbody_factor = low.radiance_w_m2_sr_nm(580.0).unwrap() / low_e;
    close(
        high.radiance_w_m2_sr_nm(580.0).unwrap(),
        blackbody_factor * high_e,
    );
    let low_response =
        law.resolve_response(&native, CellId(11), MaterialId(7), [580.0]).unwrap();
    let high_response =
        law.resolve_response(&native, CellId(12), MaterialId(7), [580.0]).unwrap();
    close(low_response.absorbed_fraction()[0], low_e);
    close(high_response.absorbed_fraction()[0], high_e);
}
