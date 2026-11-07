//! Local source temperature controls molecular friction without creating heat.

use super::*;
use fs_lbm::material::free_surface::wetting::{
    MolecularKineticControls, MolecularKineticError, MolecularKineticLaw,
};

#[test]
fn g1_g4_heated_molecular_wetting_uses_source_temperature_and_refuses_atomically() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    // The only supported source speed is zero. Predicted moving angles must
    // retain these static claims, even while their actual temperature changes.
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 0.0);
    let source = liquid_source(&liquid, false);
    let law = LiquidTemperatureLaw::try_from_material(
        &liquid, source.material().clone(), "T", LiquidViscosityBasis::Dynamic,
    ).unwrap();
    let chart = caloric(&liquid);
    let gate = CancelGate::new_clock_free();

    for (maximum_temperature_k, activation_temperature_k, reference_frequency_per_s) in [
        (500.0, 0.0, 1e12),
        (410.0, 0.0, 1e12),
        (500.0, 2000.0, 1e12),
        (410.0, 2000.0, 1e12),
        (500.0, 0.0, 1e9),
    ]
    {
        with_cx(&gate, |cx| {
            let initial = cells().into_iter().map(|mut cell| {
                cell.velocity_m_s = [0.0; 3];
                cell
            }).collect::<Vec<_>>();
            let flow = plan_free_surface_from_material(&source, 0.001, 0.00002, 0.1)
                .unwrap()
                .free_surface3_from_cells(domain(), &initial, cx)
                .unwrap();
            let controls = MolecularKineticControls {
                // Synthetic model calibration, deliberately separate from the
                // material card and the static contact-angle measurements.
                // At 1 THz, 2*lambda*K_ref = 1000 m/s. With sigma=0.004 N/m
                // and theta_0=1.2, both signs remain angle-admissible below
                // 46 m/s throughout 400..500 K, including activated kinetics.
                // The original 1 GHz calibration remains a refusal case.
                jump_length_m: 5e-10,
                equilibrium_jump_frequency_per_s: reference_frequency_per_s,
                minimum_temperature_k: 300.0,
                maximum_temperature_k,
                maximum_speed_m_s: 100.0,
                minimum_angle_radians: 0.01,
                maximum_angle_radians: 3.0,
            };
            let model = if activation_temperature_k == 0.0 {
                MolecularKineticLaw::try_new(controls)
            } else {
                MolecularKineticLaw::try_new_arrhenius(
                    controls, 400.0, activation_temperature_k * 1.380_649e-23,
                )
            }.unwrap();
            let mut thermal = flow.with_enthalpy(
                &law, &chart, FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0, &gate,
            ).unwrap()
                .with_dynamic_planar_wetting(
                    PlanarWall3::XMin, &interface,
                    MaterialPropertySelection::SingleClaimOnly, wetting_controls(), &gate,
                ).unwrap()
                .with_molecular_kinetic_wetting(model, &gate).unwrap();
            let before = thermal.total_enthalpy_j().unwrap();
            let mass = thermal.flow().ledger_mass_kg().unwrap();
            let bits = thermal.flow().lattice().grid().state_bits();
            let samples = thermal.wetting_samples().unwrap().clone();
            let material = thermal.material_at(WALL_CELL).unwrap().clone();
            let heat = initial.iter().map(|cell| FreeSurfaceCellHeat {
                coordinates: cell.coordinates, energy_j: 0.04,
            }).collect::<Vec<_>>();
            let supplied = heat.iter().map(|entry| entry.energy_j).sum::<f64>();
            let result = thermal.step_with_heat(&heat, &gate);

            let expected_refusal = if maximum_temperature_k == 410.0 {
                Some("temperature")
            } else if reference_frequency_per_s == 1e9 {
                Some("requested speed has no finite contact angle")
            } else {
                None
            };
            if let Some(expected_reason) = expected_refusal {
                assert!(matches!(
                    result,
                    Err(FreeSurfaceEnthalpyError::Wetting(
                        FreeSurfaceWettingError::MolecularKinetic(
                            MolecularKineticError::OutsideDomain(reason)
                        )
                    )) if reason == expected_reason
                ));
                assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
                assert_eq!(thermal.flow().ledger_mass_kg().unwrap(), mass);
                assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
                assert_eq!(thermal.wetting_samples().unwrap(), &samples);
                assert_eq!(thermal.material_at(WALL_CELL).unwrap(), &material);
                return;
            }

            let report = result.unwrap();
            near(thermal.total_enthalpy_j().unwrap() - before, supplied, 2e-9);
            near(report.external_heat_j, supplied, 2e-12);
            assert_eq!(report.viscous_heat_j, 0.0);
            let mut heated_moving_samples = 0;
            for (&position, sample) in thermal.wetting_samples().unwrap() {
                let prediction = sample.molecular_kinetic.unwrap();
                assert!(sample.hydrodynamic.is_none());
                assert_eq!(
                    sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS), Some(0.0)
                );
                assert_eq!(
                    sample.material.state_coordinate("T"), Some(prediction.temperature_k)
                );
                assert_eq!(prediction.surface_tension_n_m, 0.004);
                if let Some(state) = thermal.material_at(position) {
                    near(prediction.temperature_k, temperature(state), 2e-12);
                    if prediction.signed_speed_m_s != 0.0 && temperature(state) > 401.0 {
                        heated_moving_samples += 1;
                    }
                } else {
                    assert_eq!(prediction.temperature_k, 400.0);
                }
                let frequency_ratio = (activation_temperature_k
                    * (1.0 / 400.0 - 1.0 / prediction.temperature_k)).exp();
                let expected_frequency = reference_frequency_per_s * frequency_ratio;
                near(prediction.equilibrium_jump_frequency_per_s, expected_frequency, 2e-12);
                let velocity_scale = 2.0 * controls.jump_length_m * expected_frequency;
                let expected_force = 2.0 * 1.380_649e-23 * prediction.temperature_k
                    / controls.jump_length_m.powi(2)
                    * (prediction.signed_speed_m_s / velocity_scale).asinh();
                near(prediction.driving_force_per_length_n_m, expected_force, 2e-12);
                assert!(prediction.signed_speed_m_s.abs() < 46.0);
                assert_eq!(prediction.applied_angle_radians, sample.applied_angle_radians);
                assert!(prediction.dissipation_per_length_w_m >= 0.0);
                if sample.regime == ContactLineRegime2::Pinned {
                    assert_eq!(prediction.signed_speed_m_s, 0.0);
                    assert_eq!(prediction.dissipation_per_length_w_m, 0.0);
                }
            }
            assert!(heated_moving_samples > 0);
        });
    }
}
