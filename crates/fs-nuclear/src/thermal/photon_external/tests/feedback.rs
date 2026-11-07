//! G1/G3/G4: temperature feedback on actual photon and native thermal owners.

use super::*;
use crate::radiation::transport::material::temperature::{
    PhotonMassTemperatureKnot, TemperaturePhotonMassMaterial,
};

fn temperature_layers(growing: bool, maximum_k: f64) -> Vec<TemperaturePhotonSlabLayer> {
    layers(false)
        .into_iter()
        .map(|layer| {
            let original = layer.material.as_ref().clone();
            let mut hot = original.groups().to_vec();
            if growing && layer.cell_id == CellId(11) {
                for row in &mut hot {
                    // tau_0(T)=0.5+0.05*(T-300) at the original density 2.
                    row.total_m2_kg += 0.025 * (maximum_k - 300.0);
                    row.heat_m2_kg = row.total_m2_kg;
                }
            }
            let hot = PhotonMassMaterial::new(
                format!("manufactured hot material {}", original.material_id()),
                original.material_id(),
                original.radiation(),
                original.group_edges_j().to_vec(),
                hot,
                PhotonSlabLimits::default(),
                || Ok(()),
            )
            .unwrap();
            TemperaturePhotonSlabLayer {
                cell_id: layer.cell_id,
                material: Arc::new(
                    TemperaturePhotonMassMaterial::new(
                        format!("manufactured temperature data {}", original.material_id()),
                        vec![
                            PhotonMassTemperatureKnot {
                                temperature_k: 300.0,
                                material: original,
                            },
                            PhotonMassTemperatureKnot {
                                temperature_k: maximum_k,
                                material: hot,
                            },
                        ],
                        PhotonSlabLimits::default(),
                        || Ok(()),
                    )
                    .unwrap(),
                ),
            }
        })
        .collect()
}

fn feedback(
    thermal: &DepositionSimulation<'_>,
    layers: Vec<TemperaturePhotonSlabLayer>,
    sources: Vec<PhotonIllumination>,
    refinement: Option<PhotonSlabRefinement>,
    feedback_budget: PhotonFeedbackBudget,
) -> Result<ExternalPhotonHeating, FissionThermalError> {
    ExternalPhotonHeating::from_temperature_slab(
        thermal,
        "manufactured temperature-dependent response".into(),
        geometry(),
        layers,
        sources,
        vec![PhotonOrdinate {
            mu: 1.0,
            weight: 1.0,
        }],
        PhotonSlabLimits::default(),
        refinement,
        allowance(),
        feedback_budget,
        || Ok(()),
    )
}

fn steady() -> PhotonIllumination {
    beam(
        9001,
        PhotonSlabBoundary::Left,
        1.0,
        schedule(&[(0.0, [8.0, 0.0])]),
    )
}

#[test]
fn g1_heated_shielding_converges_to_an_independent_nonlinear_solution() {
    let original = thermal(&LAW, 2.0, 0.0);
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let target = 4.0;
    let a = 0.5_f64;
    let b = 0.05_f64;
    // dtheta/dt=P/C*(1-exp(-a-b*theta)); exp(b*theta) obeys
    // a linear ODE. This reference does not call the transport owner.
    let theta = ((-a).exp() + (1.0 - (-a).exp()) * (8.0 * b * target / 4.0).exp()).ln() / b;
    let exact_second = 300.0 + (8.0 * target - 4.0 * theta) * (1.0 - (-1.0_f64).exp()) / 6.0;
    let mut errors = Vec::new();
    for dt in [0.5, 0.25, 0.125] {
        let (heated, accepted, thermal_audit) = owner
            .advance_to(&original, target, controls(dt), || Ok(()))
            .unwrap();
        let t1 = heated.temperature_k(CellId(11)).unwrap();
        let t2 = heated.temperature_k(CellId(29)).unwrap();
        errors.push((t1 - (300.0 + theta)).abs());
        assert!((t2 - exact_second).abs() < errors.last().unwrap() * 4.0 / 6.0);
        let heat = 4.0 * (t1 - 300.0) + 6.0 * (t2 - 300.0);
        close(accepted.audit().incident_j, 32.0);
        close(accepted.audit().absorbed_j, heat);
        close(accepted.audit().escaped_j + heat, 32.0);
        close(thermal_audit.external_heat_j, heat);
        assert_eq!(
            accepted.feedback_work().unwrap().preparations + 1,
            accepted.audit().substeps
        );
        assert!(accepted.prepared_temperatures_k().unwrap()[0].1 < t1);
        assert_eq!(heated.cells()[&CellId(11)].mass_kg, 2.0);
        assert_eq!(
            heated.cells()[&CellId(11)].geometry,
            original.cells()[&CellId(11)].geometry
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.2,
            "{errors:?}"
        );
    }
    let (fixed, _, _) = bind(&original, vec![steady()], false)
        .advance_to(&original, target, controls(0.125), || Ok(()))
        .unwrap();
    assert!(300.0 + theta > fixed.temperature_k(CellId(11)).unwrap() + 0.2);
    assert!(exact_second < fixed.temperature_k(CellId(29)).unwrap());
}

#[test]
fn g3_constant_tables_reproduce_fixed_transport_and_native_restart() {
    let original = thermal(&LAW, 2.0, 0.0);
    let fixed = bind(&original, vec![pulse()], false);
    let owner = feedback(
        &original,
        temperature_layers(false, 350.0),
        vec![pulse()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let (expected, _, _) = fixed
        .advance_to(&original, 3.0, controls(0.25), || Ok(()))
        .unwrap();
    let (whole, accepted, _) = owner
        .advance_to(&original, 3.0, controls(0.25), || Ok(()))
        .unwrap();
    let (middle, prefix, _) = owner
        .advance_to(&original, 1.5, controls(0.25), || Ok(()))
        .unwrap();
    let (split, suffix, _) = prefix
        .advance_to(&middle, 3.0, controls(0.25), || Ok(()))
        .unwrap();
    for id in [CellId(11), CellId(29)] {
        assert_eq!(
            whole.cells()[&id].specific_h_j_kg,
            expected.cells()[&id].specific_h_j_kg
        );
        assert_eq!(
            whole.cells()[&id].specific_h_j_kg,
            split.cells()[&id].specific_h_j_kg
        );
    }
    close(
        prefix.audit().absorbed_j + suffix.audit().absorbed_j,
        accepted.audit().absorbed_j,
    );
    close(
        prefix.audit().escaped_j + suffix.audit().escaped_j,
        accepted.audit().escaped_j,
    );
    let refined = feedback(
        &original,
        temperature_layers(false, 350.0),
        vec![pulse()],
        Some(PhotonSlabRefinement {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-10,
            max_levels: 2,
        }),
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let (fine, _, _) = refined
        .advance_to(&original, 3.0, controls(0.25), || Ok(()))
        .unwrap();
    for id in [CellId(11), CellId(29)] {
        close(
            fine.temperature_k(id).unwrap(),
            whole.temperature_k(id).unwrap(),
        );
    }
}

#[test]
fn g1_native_temperature_changes_heat_and_multigroup_scattering_partition() {
    let original = thermal(&LAW, 2.0, 1.0);
    let mut materials = temperature_layers(false, 320.0);
    let make = |hot: bool| {
        PhotonMassMaterial::new(
            if hot {
                "manufactured warm scattering"
            } else {
                "manufactured cold scattering"
            }
            .into(),
            1,
            DecayRadiation::Gamma,
            edges(),
            vec![
                PhotonMassCollision {
                    total_m2_kg: 0.25,
                    heat_m2_kg: if hot { 0.25 } else { 0.0 },
                    redistribution_m2_kg: if hot { vec![0.0, 0.0] } else { vec![0.25, 0.0] },
                },
                PhotonMassCollision {
                    total_m2_kg: 0.5,
                    heat_m2_kg: if hot { 0.25 } else { 0.0 },
                    redistribution_m2_kg: if hot {
                        vec![0.25, 0.0]
                    } else {
                        vec![0.25, 0.25]
                    },
                },
            ],
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap()
    };
    materials[0].material = Arc::new(
        TemperaturePhotonMassMaterial::new(
            "manufactured heat/scatter transition".into(),
            vec![
                PhotonMassTemperatureKnot {
                    temperature_k: 300.0,
                    material: make(false),
                },
                PhotonMassTemperatureKnot {
                    temperature_k: 302.0,
                    material: make(true),
                },
            ],
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    let owner = feedback(
        &original,
        materials,
        vec![beam(
            9001,
            PhotonSlabBoundary::Left,
            1.0,
            schedule(&[(0.0, [0.0, 1.0])]),
        )],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let cold = owner
        .response()
        .apply(&[(9001, &[0.0, 1.0])], 1e-8, || Ok(()))
        .unwrap();
    assert_eq!(cold.receiver_energy_j[&11], 0.0);
    let (heated, accepted, audit) = owner
        .advance_to(&original, 1.0, controls(0.25), || Ok(()))
        .unwrap();
    let warm = accepted
        .response()
        .apply(&[(9001, &[0.0, 1.0])], 1e-8, || Ok(()))
        .unwrap();
    assert!(warm.receiver_energy_j[&11] > 0.0);
    assert!(heated.temperature_k(CellId(11)).unwrap() > 300.25);
    close(
        accepted.audit().absorbed_j + accepted.audit().escaped_j,
        1.0,
    );
    close(audit.external_heat_j, accepted.audit().absorbed_j + 1.0);
    assert_eq!(accepted.feedback_work().unwrap().preparations, 3);
}

struct Latent;
impl ThermalLaw for Latent {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([41; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || h <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 600.0 {
            h / 2.0
        } else if h <= 620.0 {
            300.0
        } else {
            (h - 20.0) / 2.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.5
    }
}
static LATENT: Latent = Latent;

#[test]
fn g3_latent_heat_keeps_the_original_calorics_and_reuses_unchanged_response() {
    let original = thermal(&LATENT, 2.0, 0.0);
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let (heated, accepted, audit) = owner
        .advance_to(&original, 2.0, controls(0.25), || Ok(()))
        .unwrap();
    assert_eq!(
        accepted.feedback_work().unwrap(),
        PhotonFeedbackWork::default()
    );
    assert!(Arc::ptr_eq(owner.response(), accepted.response()));
    for id in [CellId(11), CellId(29)] {
        assert_eq!(heated.temperature_k(id).unwrap(), 300.0);
        assert!(heated.cells()[&id].specific_h_j_kg > 600.0);
    }
    close(audit.stored_energy_change_j, accepted.audit().absorbed_j);
    close(audit.external_heat_j, accepted.audit().absorbed_j);

    let original = thermal(&LAW, 2.0, 1.0);
    let dark = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![beam(
            9001,
            PhotonSlabBoundary::Left,
            1.0,
            schedule(&[(0.0, [0.0, 0.0])]),
        )],
        None,
        PhotonFeedbackBudget {
            max_transport_preparations: 1,
            ..PhotonFeedbackBudget::default()
        },
    )
    .unwrap();
    let (heated, accepted, _) = dark
        .advance_to(&original, 2.0, controls(0.25), || Ok(()))
        .unwrap();
    assert_eq!(
        accepted.feedback_work().unwrap(),
        PhotonFeedbackWork::default()
    );
    assert_eq!(accepted.audit().incident_j, 0.0);
    assert!(Arc::ptr_eq(dark.response(), accepted.response()));
    close(heated.temperature_k(CellId(11)).unwrap(), 300.5);
}

#[test]
fn g4_final_temperature_domain_and_preparation_budgets_refuse_atomically() {
    let original = thermal(&LAW, 2.0, 0.0);
    let narrow = feedback(
        &original,
        temperature_layers(true, 300.01),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    assert!(matches!(
        narrow.advance_to(&original, 1.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::PhotonTemperatureDomain {
            cell_id: CellId(11),
            ..
        })
    ));
    assert_eq!(narrow.time_s(), 0.0);
    assert_eq!(narrow.audit(), ExternalPhotonHeatingAudit::default());
    assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
    let tiny = PhotonFeedbackBudget {
        max_transport_preparations: 1,
        ..PhotonFeedbackBudget::default()
    };
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        tiny,
    )
    .unwrap();
    let initial_work = owner.feedback_work().unwrap();
    assert!(matches!(
        owner.advance_to(&original, 1.0, controls(0.25), || Ok(())),
        Err(FissionThermalError::Radiation(RadiationError::Budget(_)))
    ));
    assert_eq!(owner.feedback_work().unwrap(), initial_work);
    assert_eq!(
        owner.prepared_temperatures_k().unwrap(),
        &[(CellId(11), 300.0), (CellId(29), 300.0)]
    );
    let work_limited = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget {
            max_transport_work: initial_work.transport_work,
            ..PhotonFeedbackBudget::default()
        },
    )
    .unwrap();
    assert!(matches!(
        work_limited.advance_to(&original, 1.0, controls(0.25), || Ok(())),
        Err(FissionThermalError::Radiation(RadiationError::Budget(_)))
    ));
    assert!(matches!(
        feedback(
            &original,
            temperature_layers(true, 320.0),
            vec![steady()],
            None,
            PhotonFeedbackBudget {
                max_material_values: 47,
                ..PhotonFeedbackBudget::default()
            },
        ),
        Err(FissionThermalError::Radiation(RadiationError::Budget(_)))
    ));
}

#[test]
fn g4_feedback_cancellation_preserves_prepared_state_and_original_error() {
    let original = thermal(&LAW, 2.0, 0.0);
    let owner = feedback(
        &original,
        temperature_layers(true, 320.0),
        vec![steady()],
        None,
        PhotonFeedbackBudget::default(),
    )
    .unwrap();
    let mut total_calls = 0;
    owner
        .advance_to(&original, 1.0, controls(0.25), || {
            total_calls += 1;
            Ok(())
        })
        .unwrap();
    let expected = DepositionError::InvalidInput("original photon feedback callback");
    let original_response = Arc::clone(owner.response());
    let original_work = owner.feedback_work();
    for stop in [1, total_calls / 3, 2 * total_calls / 3, total_calls] {
        let mut calls = 0;
        let result = owner.advance_to(&original, 1.0, controls(0.25), || {
            calls += 1;
            if calls == stop {
                Err(expected.clone())
            } else {
                Ok(())
            }
        });
        assert!(matches!(result, Err(FissionThermalError::Thermal(error)) if error == expected));
        assert!(Arc::ptr_eq(&original_response, owner.response()));
        assert_eq!(owner.feedback_work(), original_work);
        assert_eq!(owner.time_s(), 0.0);
        assert_eq!(owner.audit(), ExternalPhotonHeatingAudit::default());
    }
}
