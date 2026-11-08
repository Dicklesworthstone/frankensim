//! G1/G3/G4: native radioactive photons with temperature-dependent shielding.

use super::*;
use crate::radiation::transport::material::temperature::{
    PhotonMassTemperatureKnot, TemperaturePhotonMassMaterial,
};

fn fixed_layers() -> Vec<PhotonSlabLayer> {
    [(1, [0.0, 0.0]), (2, [0.5, 0.5]), (3, [0.25, 0.25])]
        .into_iter()
        .map(|(id, coefficients)| PhotonSlabLayer {
            cell_id: CellId(id),
            material: material(id, coefficients, false),
        })
        .collect()
}

fn temperature_layers(growing: bool, shield_maximum_k: f64) -> Vec<TemperaturePhotonSlabLayer> {
    fixed_layers()
        .into_iter()
        .map(|layer| {
            let maximum_k = if layer.cell_id == CellId(2) {
                shield_maximum_k
            } else {
                350.0
            };
            let original = layer.material.as_ref().clone();
            let mut hot = original.groups().to_vec();
            if growing && layer.cell_id == CellId(2) {
                for row in &mut hot {
                    // At the original density 1 kg/m³ and depth 1 m:
                    // tau(T) = 0.5 + 0.5*(T - 300 K), in both energy groups.
                    row.total_m2_kg += 0.5 * (maximum_k - 300.0);
                    row.heat_m2_kg = row.total_m2_kg;
                }
            }
            let hot = PhotonMassMaterial::new(
                format!(
                    "manufactured hot internal-photon material {}",
                    original.material_id()
                ),
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
                        format!(
                            "manufactured internal-photon temperature law {}",
                            original.material_id()
                        ),
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

fn bind_temperature(
    thermal: &DepositionSimulation<'_>,
    partition: Arc<PreparedDecayRadiation>,
    growing: bool,
    maximum_k: f64,
) -> NuclearHeating<DecayCellSource> {
    plain(thermal, &partition)
        .with_temperature_photon_slab(
            thermal,
            "manufactured radioactive thermal shielding".into(),
            geometry(),
            temperature_layers(growing, maximum_k),
            vec![DecayRadiationEmitter {
                cell_id: CellId(1),
                partition,
            }],
            vec![PhotonOrdinate {
                mu: 1.0,
                weight: 1.0,
            }],
            PhotonSlabLimits::default(),
            None,
            PhotonFeedbackBudget::default(),
            || Ok(()),
        )
        .unwrap()
        .0
}

fn emitted_j(time_s: f64) -> f64 {
    // Two equal-rate serial decays with 1 s half-lives. This closed Bateman
    // energy integral does not call the nuclear or photon implementation.
    let x = std::f64::consts::LN_2 * time_s;
    INITIAL * EV_J * (1e6 * -(-x).exp_m1() + 2e6 * (1.0 - (-x).exp() * (1.0 + x)))
}

#[test]
fn g1_internal_photon_temperature_feedback_follows_decay_energy_and_nonlinear_absorption() {
    let original = thermal_grid(cells(1.0));
    let source = bind_temperature(&original, partition(), true, 350.0);
    let photons = emitted_j(1.0);
    let capacity = 2.0;
    let a = 0.5_f64;
    let b = 0.5_f64;
    // Half the isotropic emission reaches the shield. Reparameterizing time
    // by total emitted joules Q gives C*dtheta/dQ=(1-exp(-a-b*theta))/2.
    // The resulting exponential solution includes the changing daughter spectrum.
    let exact_theta =
        ((-a).exp() - (-a).exp_m1() * (b * photons / (2.0 * capacity)).exp()).ln() / b;
    let mut errors = Vec::new();
    for dt in [0.125, 0.0625, 0.03125] {
        let (heated, evolved, thermal_audit) = source
            .advance_to(&original, 1.0, controls(dt), || Ok(()))
            .unwrap();
        let shield = heated.temperature_k(CellId(2)).unwrap();
        errors.push((shield - (300.0 + exact_theta)).abs());
        let shield_heat = capacity * (shield - 300.0);
        let receiver_heat = (0.5 * photons - shield_heat) * -(-1.0_f64).exp_m1();
        close(heated.temperature_k(CellId(1)).unwrap(), 300.0, 1e-12);
        close(
            heated.temperature_k(CellId(3)).unwrap(),
            300.0 + receiver_heat / 16.0,
            2e-12,
        );
        let audit = evolved.audit();
        close(audit.released.radiation_j, photons, 2e-11);
        close(
            audit.radiation_deposited_j,
            shield_heat + receiver_heat,
            2e-11,
        );
        close(
            audit.radiation_escaped_j,
            photons - shield_heat - receiver_heat,
            2e-11,
        );
        close(audit.thermal_credit_j, shield_heat + receiver_heat, 2e-11);
        close(
            thermal_audit.external_heat_j,
            shield_heat + receiver_heat,
            2e-11,
        );
        assert_eq!(audit.released.deposited_j, 0.0);
        assert_eq!(audit.radiation_ungrouped_j, 0.0);
        let populations = evolved.source(CellId(1)).unwrap().inventory.populations();
        close(populations[0], INITIAL * 0.5, 100.0);
        close(
            populations[1],
            INITIAL * 0.5 * std::f64::consts::LN_2,
            100.0,
        );
        assert_eq!(
            evolved.photon_feedback_work().unwrap().preparations + 1,
            audit.substeps
        );
        assert!(evolved.prepared_photon_temperatures_k().unwrap()[1].1 < shield);
        assert!(!Arc::ptr_eq(
            source.grouped_radiation_response().unwrap(),
            evolved.grouped_radiation_response().unwrap()
        ));
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.2,
            "{errors:?}"
        );
    }
    assert!(errors[2] < 0.03, "{errors:?}");
    let constant_theta = 0.5 * photons * -(-a).exp_m1() / capacity;
    assert!(exact_theta - constant_theta > 0.6);
}

#[test]
fn g3_internal_photon_temperature_constants_reproduce_fixed_transport_with_local_heat_once() {
    let original = thermal_grid(cells(1.0));
    // Route half of each original gamma decay energy locally and half outward.
    let assembly = crate::radiation::decay::tests::gamma_chain(1e6, 2e6, 0.5, 0.5);
    let partition = Arc::new(
        PreparedDecayRadiation::new(
            Arc::new(
                PreparedDecayEmission::new(
                    assembly,
                    DecayRadiation::Gamma,
                    edges(),
                    DecayEmissionLimits::default(),
                    || Ok(()),
                )
                .unwrap(),
            ),
            DecayRadiationClosure {
                absolute_j_per_decay: 0.0,
                relative: 1e-12,
            },
            DecayEmissionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    let source = bind_temperature(&original, Arc::clone(&partition), false, 350.0);
    // Replacing the temperature binding must also remove its future refreshes.
    let (fixed, _) = source
        .clone()
        .with_photon_slab(
            &original,
            "manufactured fixed replacement".into(),
            geometry(),
            fixed_layers(),
            vec![DecayRadiationEmitter {
                cell_id: CellId(1),
                partition,
            }],
            vec![PhotonOrdinate {
                mu: 1.0,
                weight: 1.0,
            }],
            PhotonSlabLimits::default(),
            || Ok(()),
        )
        .unwrap();
    assert!(fixed.photon_feedback_work().is_none());
    assert!(fixed.prepared_photon_temperatures_k().is_none());
    let (heated, evolved, thermal_audit) = source
        .advance_to(&original, 1.0, controls(0.125), || Ok(()))
        .unwrap();
    let (reference, fixed, _) = fixed
        .advance_to(&original, 1.0, controls(0.125), || Ok(()))
        .unwrap();
    for id in [CellId(1), CellId(2), CellId(3)] {
        close(
            heated.temperature_k(id).unwrap(),
            reference.temperature_k(id).unwrap(),
            2e-12,
        );
    }
    assert_eq!(
        evolved.source(CellId(1)).unwrap().inventory.populations(),
        fixed.source(CellId(1)).unwrap().inventory.populations()
    );
    let total = emitted_j(1.0);
    close(evolved.audit().released.deposited_j, 0.5 * total, 2e-11);
    close(evolved.audit().released.radiation_j, 0.5 * total, 2e-11);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + 0.5 * total / 8.0,
        2e-12,
    );
    close(
        evolved.audit().thermal_credit_j,
        0.5 * total + evolved.audit().radiation_deposited_j,
        2e-11,
    );
    close(
        thermal_audit.external_heat_j,
        evolved.audit().thermal_credit_j,
        2e-11,
    );
}

#[test]
fn g4_internal_photon_temperature_refusal_keeps_inventory_response_and_cache_together() {
    let original = thermal_grid(cells(1.0));
    let source = bind_temperature(&original, partition(), true, 300.5);
    let (middle, accepted, _) = source
        .advance_to(&original, 0.25, controls(0.25), || Ok(()))
        .unwrap();
    let middle_temperature = middle.temperature_k(CellId(2)).unwrap();
    assert!(middle_temperature > 300.2 && middle_temperature < 300.5);
    let response = Arc::clone(accepted.grouped_radiation_response().unwrap());
    let temperatures = accepted.prepared_photon_temperatures_k().unwrap().to_vec();
    let work = accepted.photon_feedback_work();
    let populations = accepted
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .to_vec();
    let audit = accepted.audit();
    // The next interval must rebuild at the newly heated start temperature.
    // Its thermal endpoint exceeds the opacity domain and rejects that response
    // together with the tentative parent/daughter histories and heat credit.
    assert!(temperatures[1].1 < middle_temperature);
    let failure = accepted
        .advance_to(&middle, 0.5, controls(0.25), || Ok(()))
        .err()
        .unwrap();
    assert!(matches!(
        failure,
        FissionThermalError::PhotonTemperatureDomain {
            cell_id: CellId(2),
            ..
        }
    ));
    assert_eq!(accepted.time_s(), 0.25);
    assert_eq!(middle.time_s(), 0.25);
    assert_eq!(middle.temperature_k(CellId(2)).unwrap(), middle_temperature);
    assert_eq!(
        accepted.source(CellId(1)).unwrap().inventory.populations(),
        populations.as_slice()
    );
    assert_eq!(accepted.audit(), audit);
    assert_eq!(accepted.photon_feedback_work(), work);
    assert_eq!(
        accepted.prepared_photon_temperatures_k().unwrap(),
        temperatures.as_slice()
    );
    assert!(Arc::ptr_eq(
        &response,
        accepted.grouped_radiation_response().unwrap()
    ));
}
