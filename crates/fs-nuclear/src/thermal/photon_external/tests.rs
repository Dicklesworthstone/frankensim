//! G1/G3/G4: prescribed beam energy on actual material, transport and heat owners.

use super::*;
use crate::radiation::transport::material::{PhotonMassCollision, PhotonMassMaterial};
use crate::radiation::transport::{PhotonSlabBoundary, PhotonSourceGeometry};
use fs_thermal::enthalpy::{Boundary, DepositedCell, DepositionBudget, ThermalLaw, TimedEvent};

struct Calorics {
    maximum_k: f64,
}
impl ThermalLaw for Calorics {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([19; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = h / 2.0;
        if t.is_finite() && t > 0.0 && t <= self.maximum_k {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.5
    }
}
static LAW: Calorics = Calorics { maximum_k: 1e6 };
static LIMITED: Calorics = Calorics { maximum_k: 300.01 };

fn controls(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-8,
    }
}
fn allowance() -> PhotonHeatingBudget {
    PhotonHeatingBudget {
        max_schedule_values: 64,
        energy_tolerance_j: 1e-8,
    }
}
fn edges() -> Vec<f64> {
    vec![0.0, 1e-14, 2e-14]
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a:?} != {b:?}");
}

fn cell(id: u64, material_id: u64, start: f64, mass: f64, source_w: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id + 100,
        body_id: id + 200,
        material_id: MaterialId(material_id),
        mass_kg: mass,
        specific_h_j_kg: 600.0,
        geometry: CellGeometry {
            start_m: [start, 0.0, 0.0],
            end_m: [start + 1.0, 0.0, 0.0],
            volume_m3: 1.0,
            surface_area_m2: 6.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w,
        },
    }
}
fn thermal(
    law: &'static dyn ThermalLaw,
    first_mass: f64,
    boundary_heat: f64,
) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        vec![
            cell(11, 1, 0.0, first_mass, boundary_heat),
            cell(29, 2, 1.0, 3.0, 0.0),
        ]
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect(),
        &[(MaterialId(1), law), (MaterialId(2), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0
}
fn layers(scattering: bool) -> Vec<PhotonSlabLayer> {
    [(11, 1, [0.25, 0.5]), (29, 2, [1.0 / 3.0, 1.0 / 6.0])]
        .into_iter()
        .map(|(id, material, coeff)| PhotonSlabLayer {
            cell_id: CellId(id),
            material: Arc::new(
                PhotonMassMaterial::new(
                    format!("manufactured slab material {material}"),
                    material,
                    DecayRadiation::Gamma,
                    edges(),
                    coeff
                        .into_iter()
                        .enumerate()
                        .map(|(g, coefficient)| {
                            let mut redistribution = vec![0.0; 2];
                            if scattering {
                                redistribution[g] = coefficient;
                            }
                            PhotonMassCollision {
                                total_m2_kg: coefficient,
                                heat_m2_kg: if scattering { 0.0 } else { coefficient },
                                redistribution_m2_kg: redistribution,
                            }
                        })
                        .collect(),
                    PhotonSlabLimits::default(),
                    || Ok(()),
                )
                .unwrap(),
            ),
        })
        .collect()
}
fn schedule(knots: &[(f64, [f64; 2])]) -> PhotonPowerSchedule {
    PhotonPowerSchedule::new(
        "authored beam history".into(),
        DecayRadiation::Gamma,
        edges(),
        knots
            .iter()
            .map(|(time, power)| PhotonPowerKnot {
                time_s: *time,
                group_power_w: power.to_vec(),
            })
            .collect(),
        64,
        || Ok(()),
    )
    .unwrap()
}
fn beam(
    id: u64,
    boundary: PhotonSlabBoundary,
    mu: f64,
    power: PhotonPowerSchedule,
) -> PhotonIllumination {
    PhotonIllumination {
        source: PhotonSlabSource {
            source_id: id,
            geometry: PhotonSourceGeometry::Collimated { boundary, mu },
        },
        power,
    }
}
fn pulse() -> PhotonIllumination {
    beam(
        9001,
        PhotonSlabBoundary::Left,
        0.5,
        schedule(&[(0.0, [0.0, 0.0]), (1.0, [20.0, 10.0]), (3.0, [0.0, 0.0])]),
    )
}
fn geometry() -> PhotonSlabGeometry {
    PhotonSlabGeometry {
        area_m2: 1.0,
        position_tolerance_m: 1e-12,
    }
}
fn bind(
    thermal: &DepositionSimulation<'_>,
    sources: Vec<PhotonIllumination>,
    scattering: bool,
) -> ExternalPhotonHeating {
    ExternalPhotonHeating::from_slab(
        thermal,
        "manufactured planar response".into(),
        geometry(),
        layers(scattering),
        sources,
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        PhotonSlabLimits::default(),
        allowance(),
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn g1_external_spectral_pulse_heats_native_cells_and_clips_at_true_knots() {
    let original = thermal(&LAW, 2.0, 1.0);
    let heating = bind(&original, vec![pulse()], false);
    assert_eq!(heating.sources()[0].0, 9001); // No source cell with this identity exists.
    assert_eq!(heating.sources()[0].1.data_id(), "authored beam history");
    let (heated, accepted, thermal_audit) = heating
        .advance_to(&original, 3.0, controls(100.0), || Ok(()))
        .unwrap();
    let first = 30.0 * (1.0 - (-1.0_f64).exp()) + 15.0 * (1.0 - (-2.0_f64).exp());
    let second = 30.0 * (-1.0_f64).exp() * (1.0 - (-2.0_f64).exp())
        + 15.0 * (-2.0_f64).exp() * (1.0 - (-1.0_f64).exp());
    let escaped = 45.0 * (-3.0_f64).exp();
    close(
        heated.temperature_k(CellId(11)).unwrap(),
        300.0 + (first + 3.0) / 4.0,
    );
    close(
        heated.temperature_k(CellId(29)).unwrap(),
        300.0 + second / 6.0,
    );
    let audit = accepted.audit();
    assert_eq!(audit.substeps, 2);
    close(audit.incident_j, 45.0);
    close(audit.absorbed_j, first + second);
    close(audit.escaped_j, escaped);
    close(audit.thermal_credit_j, audit.absorbed_j);
    close(audit.residual_j, 0.0);
    close(thermal_audit.external_heat_j, first + second + 3.0);
    close(
        thermal_audit.stored_energy_change_j,
        thermal_audit.external_heat_j,
    );
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(heating.time_s(), 0.0);
    let (after, dark, _) = accepted
        .advance_to(&heated, 4.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(dark.audit().incident_j, 0.0);
    close(
        after.temperature_k(CellId(29)).unwrap(),
        heated.temperature_k(CellId(29)).unwrap(),
    );
}

#[test]
fn g3_external_pulse_restart_and_response_refinement_preserve_source_energy() {
    let original = thermal(&LAW, 2.0, 0.0);
    let heating = bind(&original, vec![pulse()], false);
    let (whole, _, _) = heating
        .advance_to(&original, 3.0, controls(0.25), || Ok(()))
        .unwrap();
    let (middle, partial, _) = heating
        .advance_to(&original, 1.5, controls(0.25), || Ok(()))
        .unwrap();
    let (split, continued, _) = partial
        .advance_to(&middle, 3.0, controls(0.25), || Ok(()))
        .unwrap();
    for id in [CellId(11), CellId(29)] {
        assert_eq!(
            whole.cells()[&id].specific_h_j_kg,
            split.cells()[&id].specific_h_j_kg
        );
    }
    close(
        partial.audit().incident_j + continued.audit().incident_j,
        45.0,
    );
    let (refined, response) = ExternalPhotonHeating::from_refined_slab(
        &original,
        "refined beam".into(),
        geometry(),
        layers(false),
        vec![pulse()],
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        PhotonSlabLimits::default(),
        PhotonSlabRefinement {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            max_levels: 3,
        },
        allowance(),
        || Ok(()),
    )
    .unwrap();
    assert!(response.subdivisions_per_layer > 1);
    let (heated, refined, _) = refined
        .advance_to(&original, 3.0, controls(0.25), || Ok(()))
        .unwrap();
    close(refined.audit().incident_j, 45.0);
    for id in [CellId(11), CellId(29)] {
        close(
            heated.temperature_k(id).unwrap(),
            whole.temperature_k(id).unwrap(),
        );
    }
}

#[test]
fn g1_opposed_beams_keep_independent_spectra_and_original_receivers() {
    let original = thermal(&LAW, 2.0, 0.0);
    let sources = vec![
        beam(
            9001,
            PhotonSlabBoundary::Left,
            1.0,
            schedule(&[(0.0, [4.0, 0.0])]),
        ),
        beam(
            77,
            PhotonSlabBoundary::Right,
            1.0,
            schedule(&[(0.0, [0.0, 6.0])]),
        ),
    ];
    let heating = bind(&original, sources, false);
    let (heated, accepted, _) = heating
        .advance_to(&original, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let a = 1.0 - (-0.5_f64).exp();
    let b = (-0.5_f64).exp() * (1.0 - (-1.0_f64).exp());
    close(
        heated.temperature_k(CellId(11)).unwrap(),
        300.0 + (4.0 * a + 6.0 * b) / 4.0,
    );
    close(
        heated.temperature_k(CellId(29)).unwrap(),
        300.0 + (6.0 * a + 4.0 * b) / 6.0,
    );
    close(accepted.audit().incident_j, 10.0);
    close(accepted.audit().escaped_j, 10.0 * (-1.5_f64).exp());
    assert_eq!(
        accepted.sources().iter().map(|s| s.0).collect::<Vec<_>>(),
        vec![77, 9001]
    );
}

#[test]
fn g1_g4_scattering_only_receivers_remain_bound_without_artificial_heat() {
    let original = thermal(&LAW, 2.0, 0.0);
    let heating = bind(&original, vec![pulse()], true);
    let (same, accepted, _) = heating
        .advance_to(&original, 3.0, controls(1.0), || Ok(()))
        .unwrap();
    close(accepted.audit().incident_j, 45.0);
    close(accepted.audit().escaped_j, 45.0);
    assert_eq!(accepted.audit().absorbed_j, 0.0);
    assert_eq!(same.temperature_k(CellId(11)).unwrap(), 300.0);
    let changed_mass = thermal(&LAW, 2.1, 0.0);
    assert!(matches!(
        heating.advance_to(&changed_mass, 0.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::BindingChanged(CellId(11)))
    ));
}

#[test]
fn g4_photon_cancellation_and_caloric_refusal_discard_both_histories() {
    let original = thermal(&LAW, 2.0, 0.0);
    let heating = bind(&original, vec![pulse()], false);
    let mut calls = 0;
    let cancelled = heating.advance_to(&original, 3.0, controls(0.1), || {
        calls += 1;
        if calls == 200 {
            Err(DepositionError::InvalidInput("exact caller cancellation"))
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(FissionThermalError::Thermal(DepositionError::InvalidInput(
            "exact caller cancellation"
        )))
    ));
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(heating.time_s(), 0.0);
    assert_eq!(heating.audit(), ExternalPhotonHeatingAudit::default());
    let limited = thermal(&LIMITED, 2.0, 0.0);
    let source = bind(&limited, vec![pulse()], false);
    assert!(
        source
            .advance_to(&limited, 3.0, controls(1.0), || Ok(()))
            .is_err()
    );
    assert_eq!(limited.temperature_k(CellId(11)).unwrap(), 300.0);
    assert_eq!(source.time_s(), 0.0);
    let wrong_clock = original
        .advance_to(1.0, controls(1.0), || Ok(()))
        .unwrap()
        .0;
    assert!(matches!(
        heating.advance_to(&wrong_clock, 1.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::ClockMismatch { .. })
    ));
}

#[test]
fn g0_g4_source_coverage_grid_budget_and_knot_refusals_are_early() {
    let original = thermal(&LAW, 2.0, 0.0);
    let heating = bind(&original, vec![pulse()], false);
    let receivers = vec![
        RadiationReceiver {
            cell_id: CellId(11),
            material_id: MaterialId(1),
        },
        RadiationReceiver {
            cell_id: CellId(29),
            material_id: MaterialId(2),
        },
    ];
    let mut wrong = schedule(&[(0.0, [1.0, 2.0])]);
    wrong.edges_j = vec![0.0, 0.5e-14, 2e-14].into();
    assert!(
        ExternalPhotonHeating::new(
            &original,
            Arc::clone(heating.response()),
            vec![(9001, wrong)],
            receivers.clone(),
            allowance(),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        ExternalPhotonHeating::new(
            &original,
            Arc::clone(heating.response()),
            vec![(11, pulse().power)],
            receivers.clone(),
            allowance(),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        ExternalPhotonHeating::new(
            &original,
            Arc::clone(heating.response()),
            vec![(9001, pulse().power)],
            receivers,
            PhotonHeatingBudget {
                max_schedule_values: 5,
                ..allowance()
            },
            || Ok(())
        )
        .is_err()
    );
    assert!(
        PhotonPowerSchedule::new(
            "bad".into(),
            DecayRadiation::Gamma,
            edges(),
            vec![PhotonPowerKnot {
                time_s: 1.0,
                group_power_w: vec![1.0, 2.0]
            }],
            8,
            || Ok(())
        )
        .is_err()
    );
    let mut direct = heating.clone();
    assert!(
        direct
            .step_sources(&original, 0.0, 2.0, 2.0, &mut || Ok(()))
            .is_err()
    );
    assert!(direct.pending.is_none());
    let mut tiny = PhotonSlabLimits::default();
    tiny.max_cells = 1;
    assert!(
        ExternalPhotonHeating::from_slab(
            &original,
            "bounded".into(),
            geometry(),
            layers(false),
            vec![pulse()],
            vec![PhotonOrdinate {
                mu: 1.0,
                weight: 1.0
            }],
            tiny,
            allowance(),
            || Ok(())
        )
        .is_err()
    );
}
