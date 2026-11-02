//! G1/G3: independent native film states drive one coherent heat partition.

use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack,
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorOptics,
    ConductorSource,
};
use fs_render::dielectric::CauchyIor;
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;
use fs_render::thermal_surface::heating::{
    CoatedSurfaceIlluminationSource, CoatedSurfacePowerResponse, CoatingOpticsModel,
    StatefulCoatedSurfaceIlluminationSource, StatefulCoatingLayer, SubstrateOpticsModel,
    SurfaceCoatingLaw, SurfaceHeatingBudget, SurfaceHeatingError, SurfaceHeatingReceiver,
    SurfaceHeatingTolerance, SurfaceIncidence, SurfaceOpticalHeating,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use std::sync::Arc;

const WAVELENGTHS: [f64; 2] = [480.0, 680.0];
const S_FRACTION: [f64; 2] = [0.25, 0.75];
const COSINE: f64 = 0.8;

fn declared_source(id: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([id; 32]), ConductorDataStatus::Representative).unwrap()
}

fn film(eta: f64, k: f64, thickness_nm: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        thickness_nm * 1e-9,
        declared_source(221),
    )
    .unwrap()
}

fn metal(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        declared_source(222),
    )
    .unwrap()
}

fn receiver(id: u64) -> SurfaceHeatingReceiver {
    SurfaceHeatingReceiver {
        thermal_cell: CellId(id),
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
    }
}

fn incidence() -> SurfaceIncidence<2> {
    SurfaceIncidence {
        wavelengths_nm: WAVELENGTHS,
        incident_eta: [1.0; 2],
        incidence_cosine: COSINE,
        s_power_fraction: S_FRACTION,
    }
}

fn power(watts: [f64; 2]) -> SpectralPowerSchedule<2> {
    SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_power_w: watts,
            end_power_w: watts,
        }],
        1,
    )
    .unwrap()
}

fn budget() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 8,
        max_total_power_segments: 16,
        max_source_evaluations: 4096,
    }
}

fn tolerance() -> SurfaceHeatingTolerance {
    SurfaceHeatingTolerance {
        power_roundoff_w: 1e-10,
    }
}

fn heating(
    thermal: &DepositionSimulation<'_>,
    sources: &[StatefulCoatedSurfaceIlluminationSource<2>],
    max_rows: usize,
) -> SurfaceOpticalHeating<2> {
    SurfaceOpticalHeating::try_new_coated_with_stateful(
        &[],
        sources,
        thermal,
        budget(),
        max_rows,
        tolerance(),
        || Ok(()),
    )
    .unwrap()
}

// The reference calls the original coherent field-work owner directly. It
// neither resolves a constitutive table nor invokes a heating observer.
fn optical_watts(
    films: &[CoherentAbsorbingCoating],
    substrate: ConductorOptics,
    incident: [f64; 2],
) -> Vec<f64> {
    let layers: Vec<_> = films.iter().map(CoatingLayer::Absorbing).collect();
    let stack = CoherentCoatingStack::try_new(&layers).unwrap();
    let mut watts = vec![0.0; films.len() + 1];
    for lane in 0..2 {
        let partition = stack
            .layer_power_balance(&substrate, WAVELENGTHS[lane], 1.0, COSINE)
            .unwrap();
        for (target, fractions) in watts.iter_mut().zip(
            partition
                .layer_absorbed()
                .iter()
                .chain(core::iter::once(&partition.substrate_absorbed)),
        ) {
            *target += incident[lane]
                * (S_FRACTION[lane] * fractions[0] + (1.0 - S_FRACTION[lane]) * fractions[1]);
        }
    }
    watts
}

#[test]
fn constant_laws_retain_static_per_film_heat_and_skip_packet_rebuilds() {
    let absorbing = film(1.8, 0.2, 75.0);
    let clear = CoherentDielectricCoating::try_new(
        CauchyIor::try_new(1.45, 0.004, 0.0).unwrap(),
        110e-9,
        declared_source(223),
    )
    .unwrap();
    let substrate = metal(0.5, 3.0);
    let film_law = Arc::new(
        SurfaceCoatingLaw::try_temperature(
            declared_source(224),
            &[(300.0, absorbing), (600.0, absorbing)],
            2,
        )
        .unwrap(),
    );
    let substrate_law = Arc::new(
        SurfaceConductorLaw::try_temperature(
            declared_source(225),
            &[(300.0, substrate), (600.0, substrate)],
            2,
        )
        .unwrap(),
    );
    let source = StatefulCoatedSurfaceIlluminationSource {
        source_id: 1,
        films: vec![
            StatefulCoatingLayer {
                receiver: receiver(11),
                optics: CoatingOpticsModel::Law(film_law),
            },
            StatefulCoatingLayer {
                receiver: receiver(12),
                optics: CoatingOpticsModel::Lossless(clear),
            },
        ],
        substrate_receiver: receiver(13),
        substrate: SubstrateOpticsModel::Law(substrate_law),
        incidence: incidence(),
        power: SpectralPowerSchedule::try_new(
            &[
                SpectralPowerSegment {
                    start_s: 0.125,
                    end_s: 0.5,
                    start_power_w: [8.0, 4.0],
                    end_power_w: [8.0, 4.0],
                },
                SpectralPowerSegment {
                    start_s: 0.5,
                    end_s: 0.75,
                    start_power_w: [8.0, 4.0],
                    end_power_w: [0.0; 2],
                },
            ],
            2,
        )
        .unwrap(),
    };
    let layers = [
        CoatingLayer::Absorbing(&absorbing),
        CoatingLayer::Dielectric(&clear),
    ];
    let fixed = CoatedSurfaceIlluminationSource {
        source_id: 1,
        film_receivers: vec![receiver(11), receiver(12)],
        substrate_receiver: receiver(13),
        response: CoatedSurfacePowerResponse::from_coated_conductor(
            WAVELENGTHS,
            &substrate,
            CoherentCoatingStack::try_new(&layers).unwrap(),
            [1.0; 2],
            COSINE,
            S_FRACTION,
        )
        .unwrap(),
        power: source.power.clone(),
    };
    let mut steady = fixed.clone();
    steady.source_id = 0;
    steady.power = power([1.0, 0.5]);
    let original = thermal(
        vec![
            cell(11, 1.0, 1.0, 1.0),
            cell(12, 2.0, 1.0, 0.0),
            cell(13, 3.0, 1.0, 2.0),
        ],
        &SENSIBLE,
        vec![],
    );
    let frozen = SurfaceOpticalHeating::try_new_coated(
        &[fixed, steady.clone()],
        &original,
        budget(),
        tolerance(),
        || Ok(()),
    )
    .unwrap();
    let dynamic = SurfaceOpticalHeating::try_new_coated_with_stateful(
        &[steady],
        &[source.clone()],
        &original,
        budget(),
        4,
        tolerance(),
        || Ok(()),
    )
    .unwrap();
    let (frozen_hot, frozen_heat, _) = original
        .advance_observed(&frozen, 1.5, config(0.125), || Ok(()))
        .unwrap();
    let (dynamic_hot, dynamic_heat, _) = original
        .advance_observed(&dynamic, 1.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(frozen_hot.cells(), dynamic_hot.cells());
    assert_eq!(frozen_heat.energy(), dynamic_heat.energy());
    assert_eq!(dynamic_heat.response_evaluations(), 0);
    assert_eq!(dynamic_hot.cells()[&CellId(12)].specific_h_j_kg, 0.0);
    close(
        dynamic_hot.energy_j(),
        dynamic_heat.energy().thermal_input_j + 4.5,
    );
    assert!(matches!(
        SurfaceOpticalHeating::try_new_coated_with_stateful(
            &[],
            &[source],
            &original,
            budget(),
            3,
            tolerance(),
            || Ok(()),
        ),
        Err(SurfaceHeatingError::Budget)
    ));
}

const MASSES: [f64; 3] = [0.5, 0.75, 1.25];
const INITIAL_H: [f64; 3] = [0.0, 8.0, 16.0];
const EXTERNAL_W: [f64; 3] = [1.0, 0.0, 0.5];
const INCIDENT_W: [f64; 2] = [12.0, 8.0];

fn first_film(t: f64) -> CoherentAbsorbingCoating {
    film(1.4 + 0.08 * (t - 300.0), 0.05 + 0.05 * (t - 300.0), 75.0)
}
fn second_film(t: f64) -> CoherentAbsorbingCoating {
    film(2.3 - 0.06 * (t - 300.0), 0.2 + 0.04 * (t - 300.0), 120.0)
}
fn varying_substrate(t: f64) -> ConductorOptics {
    metal(0.3 + 0.03 * (t - 300.0), 3.0 - 0.1 * (t - 300.0))
}
fn derivative(h: [f64; 3]) -> [f64; 3] {
    let t = h.map(|h| 300.0 + h / 4.0);
    let watts = optical_watts(
        &[first_film(t[0]), second_film(t[1])],
        varying_substrate(t[2]),
        INCIDENT_W,
    );
    core::array::from_fn(|i| (watts[i] + EXTERNAL_W[i]) / MASSES[i])
}

fn reference_rk4() -> [f64; 3] {
    // Independent high-resolution integration of dh_i/dt=P_i(T_1,T_2,T_s)/m_i.
    // The native observer uses accepted-start coefficients; this fourth-order
    // reference does not use its time stepping or optical response/cache path.
    let mut h = INITIAL_H;
    let dt = 1.0 / 256.0;
    for _ in 0..256 {
        let k1 = derivative(h);
        let k2 = derivative(core::array::from_fn(|i| h[i] + 0.5 * dt * k1[i]));
        let k3 = derivative(core::array::from_fn(|i| h[i] + 0.5 * dt * k2[i]));
        let k4 = derivative(core::array::from_fn(|i| h[i] + dt * k3[i]));
        h = core::array::from_fn(|i| h[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    h
}

#[test]
fn separate_film_and_substrate_states_match_independent_ode_under_refinement() {
    let mut cells = Vec::new();
    for i in 0..3 {
        let mut native = cell(11 + i as u64, MASSES[i], 1.0, EXTERNAL_W[i]);
        native.specific_h_j_kg = INITIAL_H[i];
        cells.push(native);
    }
    let original = thermal(cells, &SENSIBLE, vec![]);
    let source = StatefulCoatedSurfaceIlluminationSource {
        source_id: 1,
        films: vec![
            StatefulCoatingLayer {
                receiver: receiver(11),
                optics: CoatingOpticsModel::Law(Arc::new(
                    SurfaceCoatingLaw::try_temperature(
                        declared_source(226),
                        &[(300.0, first_film(300.0)), (320.0, first_film(320.0))],
                        2,
                    )
                    .unwrap(),
                )),
            },
            StatefulCoatingLayer {
                receiver: receiver(12),
                optics: CoatingOpticsModel::Law(Arc::new(
                    SurfaceCoatingLaw::try_temperature(
                        declared_source(227),
                        &[(300.0, second_film(300.0)), (320.0, second_film(320.0))],
                        2,
                    )
                    .unwrap(),
                )),
            },
        ],
        substrate_receiver: receiver(13),
        substrate: SubstrateOpticsModel::Law(Arc::new(
            SurfaceConductorLaw::try_temperature(
                declared_source(228),
                &[
                    (300.0, varying_substrate(300.0)),
                    (320.0, varying_substrate(320.0)),
                ],
                2,
            )
            .unwrap(),
        )),
        incidence: incidence(),
        power: power(INCIDENT_W),
    };
    let input = heating(&original, &[source], 6);
    let (first, accepted, _) = original
        .advance_observed(&input, 0.125, config(0.125), || Ok(()))
        .unwrap();
    let initial_rate = derivative(INITIAL_H);
    let first_h = core::array::from_fn(|i| first.cells()[&CellId(11 + i as u64)].specific_h_j_kg);
    for i in 0..3 {
        close(first_h[i], INITIAL_H[i] + 0.125 * initial_rate[i]);
    }
    let (second, refreshed, _) = first
        .advance_observed(&accepted, 0.25, config(0.125), || Ok(()))
        .unwrap();
    let second_rate = derivative(first_h);
    for i in 0..3 {
        close(
            second.cells()[&CellId(11 + i as u64)].specific_h_j_kg,
            first_h[i] + 0.125 * second_rate[i],
        );
    }
    // Both films and the substrate change, but the coupled stack rebuilds once.
    assert_eq!(refreshed.response_evaluations(), 1);
    let reference = reference_rk4();
    assert!(
        (0..3)
            .map(|i| (reference[i] - INITIAL_H[i] - initial_rate[i]).abs())
            .sum::<f64>()
            > 0.01
    );
    let mut errors = Vec::new();
    for dt in [0.0625, 0.03125, 0.015625] {
        let (hot, heat, audit) = original
            .advance_observed(&input, 1.0, config(dt), || Ok(()))
            .unwrap();
        errors.push(
            (0..3)
                .map(|i| (hot.cells()[&CellId(11 + i as u64)].specific_h_j_kg - reference[i]).abs())
                .sum::<f64>(),
        );
        assert_eq!(heat.response_evaluations(), (1.0 / dt) as usize - 1);
        close(
            hot.energy_j() - original.energy_j(),
            heat.energy().thermal_input_j + 1.5,
        );
        close(audit.external_heat_j, heat.energy().thermal_input_j + 1.5);
        for lane in 0..2 {
            close(heat.energy().incident_j[lane], INCIDENT_W[lane]);
            close(
                heat.energy().reflected_j[lane] + heat.energy().absorbed_j[lane],
                INCIDENT_W[lane],
            );
            assert_eq!(heat.energy().transmitted_j[lane], 0.0);
        }
    }
    assert!(
        errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 0.2, "{errors:?}");
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
}

struct IdentifiedLatent;
impl ThermalLaw for IdentifiedLatent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        LATENT.temperature_k(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        LATENT.maximum_temperature_slope()
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([229; 32])
    }
}
static IDENTIFIED_LATENT: IdentifiedLatent = IdentifiedLatent;

fn latent_film(h: f64) -> CoherentAbsorbingCoating {
    film(1.4 + 0.01 * (h - 10.0), 0.02 + 0.0078 * (h - 10.0), 90.0)
}

#[test]
fn latent_films_share_only_identical_receiver_packets_and_survive_split_and_cancel() {
    let mut cells = vec![
        cell(11, 1.0, 1.0, 0.0),
        cell(12, 1.0, 1.0, 0.0),
        cell(13, 1.0, 1.0, 0.0),
        cell(14, 1.0, 1.0, 0.0),
    ];
    for (cell, h) in cells.iter_mut().zip([30.0, 90.0, 50.0, 50.0]) {
        cell.specific_h_j_kg = h;
    }
    let original = thermal(cells, &IDENTIFIED_LATENT, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            declared_source(230),
            ContentHash([229; 32]),
            (310.0, 310.0),
            &[(10.0, latent_film(10.0)), (110.0, latent_film(110.0))],
            2,
        )
        .unwrap(),
    );
    let make = |id, film_cell, substrate_cell, watts| StatefulCoatedSurfaceIlluminationSource {
        source_id: id,
        films: vec![StatefulCoatingLayer {
            receiver: receiver(film_cell),
            optics: CoatingOpticsModel::Law(Arc::clone(&law)),
        }],
        substrate_receiver: receiver(substrate_cell),
        substrate: SubstrateOpticsModel::Fixed(metal(0.5, 3.0)),
        incidence: incidence(),
        power: power(watts),
    };
    let mut sources = vec![
        make(3, 11, 13, [3.0, 2.0]),
        make(1, 11, 13, [2.0, 1.0]),
        make(2, 12, 14, [5.0, 3.0]),
    ];
    // Two native latent coordinates and three illuminations share just two
    // constitutive rows, but only the first two sources share a response.
    let input = heating(&original, &sources, 2);
    sources.reverse();
    let reversed = heating(&original, &sources, 2);
    assert_eq!(input.model_identity(), reversed.model_identity());
    let (first, accepted, _) = original
        .advance_observed(&input, 0.125, config(0.125), || Ok(()))
        .unwrap();
    let mut film_heat = Vec::new();
    for (film_cell, substrate_cell, initial_h) in [(11, 13, 30.0), (12, 14, 90.0)] {
        let watts = optical_watts(&[latent_film(initial_h)], metal(0.5, 3.0), [5.0, 3.0]);
        film_heat.push(watts[0]);
        close(
            first.cells()[&CellId(film_cell)].specific_h_j_kg,
            initial_h + 0.125 * watts[0],
        );
        close(
            first.cells()[&CellId(substrate_cell)].specific_h_j_kg,
            50.0 + 0.125 * watts[1],
        );
        close(first.temperature_k(CellId(film_cell)).unwrap(), 310.0);
    }
    assert!((film_heat[0] - film_heat[1]).abs() > 0.1);
    let (second, refreshed, _) = first
        .advance_observed(&accepted, 0.25, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(refreshed.response_evaluations(), 2);
    for (film_cell, substrate_cell) in [(11, 13), (12, 14)] {
        let h = first.cells()[&CellId(film_cell)].specific_h_j_kg;
        let watts = optical_watts(&[latent_film(h)], metal(0.5, 3.0), [5.0, 3.0]);
        close(
            second.cells()[&CellId(film_cell)].specific_h_j_kg,
            h + 0.125 * watts[0],
        );
        close(
            second.cells()[&CellId(substrate_cell)].specific_h_j_kg,
            first.cells()[&CellId(substrate_cell)].specific_h_j_kg + 0.125 * watts[1],
        );
        close(second.temperature_k(CellId(film_cell)).unwrap(), 310.0);
    }
    let mut checkpoints = 0;
    let (whole, complete, _) = original
        .advance_observed(&input, 1.0, config(0.125), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(complete.response_evaluations(), 14);
    let (middle, partial, _) = original
        .advance_observed(&reversed, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let (split, continued, _) = middle
        .advance_observed(&partial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(whole.cells(), split.cells());
    assert_eq!(complete.energy(), continued.energy());
    close(complete.energy().incident_j.iter().sum(), 16.0);
    close(
        whole.energy_j() - original.energy_j(),
        complete.energy().thermal_input_j,
    );
    let mut calls = 0;
    assert!(matches!(
        original.advance_observed(&input, 1.0, config(0.125), || {
            calls += 1;
            if calls == checkpoints - 2 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(input.response_evaluations(), 0);
    let (retry, retried, _) = original
        .advance_observed(&input, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), whole.cells());
    assert_eq!(retried.energy(), complete.energy());
    assert_eq!(retried.response_evaluations(), 14);
}

#[test]
fn dark_film_domains_bind_at_accepted_endpoints_without_partial_state_or_heat() {
    let law = Arc::new(
        SurfaceCoatingLaw::try_temperature(
            declared_source(231),
            &[(300.0, film(1.5, 0.1, 90.0)), (300.5, film(2.0, 0.5, 90.0))],
            2,
        )
        .unwrap(),
    );
    let original = thermal(
        vec![cell(11, 1.0, 1.0, 8.0), cell(12, 1.0, 1.0, 0.0)],
        &SENSIBLE,
        vec![],
    );
    let source = StatefulCoatedSurfaceIlluminationSource {
        source_id: 1,
        films: vec![StatefulCoatingLayer {
            receiver: receiver(11),
            optics: CoatingOpticsModel::Law(law),
        }],
        substrate_receiver: receiver(12),
        substrate: SubstrateOpticsModel::Fixed(metal(0.5, 3.0)),
        incidence: incidence(),
        power: SpectralPowerSchedule::try_new(&[], 0).unwrap(),
    };
    let input = heating(&original, &[source], 2);
    let (partial, accepted, audit) = original
        .advance_observed(&input, 0.25, config(0.125), || Ok(()))
        .unwrap();
    close(partial.temperature_k(CellId(11)).unwrap(), 300.5);
    close(audit.external_heat_j, 2.0);
    assert_eq!(accepted.response_evaluations(), 0);
    assert_eq!(accepted.energy().incident_j, [0.0; 2]);
    assert_eq!(accepted.energy().thermal_input_j, 0.0);
    assert!(matches!(
        partial.advance_observed(&accepted, 0.5, config(0.125), || Ok(())),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(11),
            ..
        })
    ));
    // A failure after earlier accepted internal intervals still returns no
    // replacement thermal state, energy account or updated optical cache.
    assert!(matches!(
        original.advance_observed(&input, 1.0, config(0.125), || Ok(())),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(11),
            ..
        })
    ));
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.response_evaluations(), 0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(partial.time_s(), 0.25);
    close(partial.cells()[&CellId(11)].specific_h_j_kg, 2.0);
    assert_eq!(accepted.time_s(), 0.25);
    assert_eq!(accepted.response_evaluations(), 0);
}
