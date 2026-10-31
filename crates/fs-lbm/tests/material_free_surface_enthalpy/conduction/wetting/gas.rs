//! Original atmosphere transport drives the native two-phase wetting boundary.
//! Manufactured source curves exercise coupling, not empirical gas accuracy.

use super::*;
use fs_lbm::material::free_surface::wetting::{CoxVoinovControls, CoxVoinovLaw};
use fs_material::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};

fn cox_controls() -> CoxVoinovControls {
    CoxVoinovControls {
        microscopic_cutoff_m: 1e-9,
        observation_length_m: 0.001,
        maximum_capillary_number: 0.01,
        minimum_angle_radians: 0.01,
        maximum_angle_radians: 3.0,
    }
}

fn liquid_card(changing: bool) -> MaterialCard {
    let base = source_card(0.004, false, 500.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == "dynamic_viscosity" {
            claim.value = if changing {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, 0.12), (500.0, 0.08)],
                    dims: claim.key.dims(),
                }
            } else {
                PropertyValue::Scalar {
                    value: 0.1,
                    dims: claim.key.dims(),
                }
            };
            claim.interpolation = if changing {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn gas_model(viscosity_scale: f64) -> ThermochemicalGasModel {
    let base = gas_source::model(100000.0, 250000.0, 100.0);
    let mut transport = base.transport_materials().to_vec();
    for binding in &mut transport {
        let mut claims = ClaimSet::new();
        for id in binding.material.claims().observation_ids() {
            assert_eq!(
                claims
                    .register_observation(
                        binding.material.claims().observation(id).unwrap().clone(),
                    )
                    .unwrap(),
                id
            );
        }
        for (_, original) in binding.material.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == "dynamic_viscosity" {
                let PropertyValue::Curve { knots, .. } = &mut claim.value else {
                    panic!("the shared synthetic gas fixture must retain its pressure curve");
                };
                for (_, value) in knots {
                    *value *= viscosity_scale;
                }
            }
            claims.insert_claim(claim).unwrap();
        }
        binding.material = MaterialCard::assemble(
            binding.material.id().clone(),
            claims,
            binding.material.models().to_vec(),
        )
        .unwrap();
    }
    ThermochemicalGasModel::new(base.thermodynamics().clone(), transport, None).unwrap()
}

fn with_wetting<R>(
    liquid: &MaterialCard,
    interface: &InterfaceSystemCard,
    outer: &ResolvedThermochemicalGasState,
    sourced: bool,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(liquid, false);
    let law = LiquidTemperatureLaw::try_from_material(
        liquid,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(liquid);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let initial = cells()
            .into_iter()
            .map(|mut cell| {
                cell.velocity_m_s = [0.0; 3];
                cell
            })
            .collect::<Vec<_>>();
        let flow = plan_free_surface_from_material(&source, 0.001, 0.00002, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial, cx)
            .unwrap();
        let thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap()
            .with_dynamic_planar_wetting(
                PlanarWall3::XMin,
                interface,
                MaterialPropertySelection::SingleClaimOnly,
                wetting_controls(),
                &gate,
            )
            .unwrap();
        let mut thermal = if sourced {
            thermal
                .with_cox_voinov_gas_wetting(cox_controls(), outer, &gate)
                .unwrap()
        } else {
            thermal
                .with_cox_voinov_wetting(
                    CoxVoinovLaw::try_new_with_viscosity_ratio(
                        cox_controls(),
                        outer.state().dynamic_viscosity / 0.1,
                        &gate,
                    )
                    .unwrap(),
                    &gate,
                )
                .unwrap()
        };
        run(&mut thermal, &gate)
    })
}

#[test]
fn constant_gas_and_liquid_sources_preserve_the_original_fixed_ratio_trajectory() {
    let liquid = liquid_card(false);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let outer = gas_model(5000.0)
        .evaluate(&gas_source::point(450.0, 100000.0))
        .unwrap();
    let run = |sourced| {
        with_wetting(&liquid, &interface, &outer, sourced, |thermal, gate| {
            let energy = thermal.total_enthalpy_j().unwrap();
            let mut states = Vec::new();
            for _ in 0..2 {
                let report = thermal.step(gate).unwrap();
                assert_eq!(report.external_heat_j, 0.0);
                assert_eq!(report.viscous_heat_j, 0.0);
                near(thermal.total_enthalpy_j().unwrap(), energy, 2e-11);
                states.push(thermal.flow().lattice().grid().state_bits());
            }
            states
        })
    };
    assert_eq!(run(true), run(false));
}

#[test]
fn changing_only_gas_transport_data_changes_contact_angles_and_following_motion() {
    let liquid = liquid_card(false);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let run = |scale| {
        let outer = gas_model(scale)
            .evaluate(&gas_source::point(450.0, 100000.0))
            .unwrap();
        with_wetting(&liquid, &interface, &outer, true, |thermal, gate| {
            thermal.step(gate).unwrap();
            let first = thermal.flow().lattice().grid().state_bits();
            let mut moving = 0;
            let angles = thermal
                .wetting_samples()
                .unwrap()
                .iter()
                .map(|(&p, sample)| {
                    assert_eq!(sample.outer_gas.as_deref(), Some(&outer));
                    let hydrodynamic = sample.hydrodynamic.unwrap();
                    // Independent pure-source pressure interpolation, not a Cox helper.
                    let mu_g = scale * 1.8e-5 * (0.5 + 99000.0 / 249000.0);
                    near(
                        hydrodynamic.outer_to_liquid_viscosity_ratio,
                        mu_g / 0.1,
                        2e-13,
                    );
                    moving += usize::from(sample.signed_speed_m_s != 0.0);
                    (p, sample.applied_angle_radians.to_bits())
                })
                .collect::<Vec<_>>();
            assert!(moving > 0);
            thermal.step(gate).unwrap();
            (first, angles, thermal.flow().lattice().grid().state_bits())
        })
    };
    let thin_outer = run(1.0);
    let viscous_outer = run(5000.0);
    assert_eq!(
        thin_outer.0, viscous_outer.0,
        "initial pinned angles coincide"
    );
    assert_ne!(thin_outer.1, viscous_outer.1);
    assert_ne!(
        thin_outer.2, viscous_outer.2,
        "gas data must reach native momentum"
    );
}

#[test]
fn liquid_heating_updates_the_ratio_while_preserving_the_independent_gas_snapshot() {
    let liquid = liquid_card(true);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let outer = gas_model(5000.0)
        .evaluate(&gas_source::point(450.0, 100000.0))
        .unwrap();
    with_wetting(&liquid, &interface, &outer, true, |thermal, gate| {
        let before = thermal.total_enthalpy_j().unwrap();
        let heat = cells()
            .iter()
            .map(|cell| FreeSurfaceCellHeat {
                coordinates: cell.coordinates,
                energy_j: 0.04,
            })
            .collect::<Vec<_>>();
        let supplied = heat.iter().map(|source| source.energy_j).sum::<f64>();
        let report = thermal.step_with_heat(&heat, gate).unwrap();
        near(thermal.total_enthalpy_j().unwrap() - before, supplied, 2e-9);
        near(report.external_heat_j, supplied, 2e-12);
        assert_eq!(report.viscous_heat_j, 0.0);
        assert_eq!(thermal.wetting_outer_gas(), Some(&outer));
        let initial_ratio = outer.state().dynamic_viscosity / 0.1;
        let mut heated_moving = 0;
        for (&position, sample) in thermal.wetting_samples().unwrap() {
            let hydrodynamic = sample.hydrodynamic.unwrap();
            assert_eq!(sample.outer_gas.as_deref(), Some(&outer));
            let mu_l = if let Some(state) = thermal.material_at(position) {
                let mu = 0.1 - 0.0002 * (temperature(state) - 400.0);
                near(state.dynamic_viscosity_pa_s(), mu, 2e-13);
                if sample.signed_speed_m_s != 0.0 && temperature(state) > 410.0 {
                    heated_moving += 1;
                    assert!(hydrodynamic.outer_to_liquid_viscosity_ratio > initial_ratio);
                }
                mu
            } else {
                0.1
            };
            near(
                hydrodynamic.outer_to_liquid_viscosity_ratio,
                outer.state().dynamic_viscosity / mu_l,
                2e-13,
            );
            near(
                hydrodynamic.signed_capillary_number,
                mu_l * sample.signed_speed_m_s / 0.004,
                2e-13,
            );
            assert_eq!(
                sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                Some(0.0)
            );
            if sample.regime != ContactLineRegime2::Pinned {
                let current = CoxVoinovLaw::try_new_with_viscosity_ratio(
                    cox_controls(),
                    hydrodynamic.outer_to_liquid_viscosity_ratio,
                    gate,
                )
                .unwrap();
                near(
                    sample.applied_angle_radians,
                    current
                        .apparent_angle_radians(1.2, hydrodynamic.signed_capillary_number)
                        .unwrap(),
                    2e-12,
                );
            }
        }
        assert!(heated_moving > 0);
    });
}

#[test]
fn finite_atmosphere_updates_wetting_from_its_actual_accepted_temperature_and_pressure() {
    let liquid = liquid_card(true);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let mut atmosphere =
        IsochoricGasParcel::new(gas_model(5000.0), gas_source::point(450.0, 100000.0), 1e-8)
            .unwrap();
    let initial_gas = atmosphere.state().clone();
    with_wetting(&liquid, &interface, &initial_gas, true, |thermal, gate| {
        let energy = thermal.total_enthalpy_j().unwrap()
            + atmosphere.mass_kg() * atmosphere.state().specific_internal_energy_j_kg();
        let report = thermal
            .step_with_finite_gas_environment(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &mut atmosphere,
                &[gas_face()],
                gas_controls(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        assert!(atmosphere.state().state().temperature < 450.0);
        assert_ne!(
            atmosphere.state().state().pressure,
            initial_gas.state().pressure
        );
        assert_eq!(thermal.wetting_outer_gas(), Some(atmosphere.state()));
        let mu_g =
            5000.0 * 1.8e-5 * (0.5 + (atmosphere.state().state().pressure - 1000.0) / 249000.0);
        near(atmosphere.state().state().dynamic_viscosity, mu_g, 2e-13);
        for sample in thermal.wetting_samples().unwrap().values() {
            assert_eq!(sample.outer_gas.as_deref(), Some(atmosphere.state()));
            let prediction = sample.hydrodynamic.unwrap();
            near(
                prediction.outer_to_liquid_viscosity_ratio,
                mu_g / prediction.dynamic_viscosity_pa_s,
                2e-13,
            );
        }
        let after = thermal.total_enthalpy_j().unwrap()
            + atmosphere.mass_kg() * atmosphere.state().specific_internal_energy_j_kg();
        near(after, energy, 2e-10);
        assert!(
            report
                .finite_gas
                .unwrap()
                .combined_thermal_balance_defect_j
                .abs()
                < 1e-9
        );
    });
}

#[test]
fn a_different_gas_source_and_cancellation_preserve_every_accepted_participant() {
    let liquid = liquid_card(false);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let outer = gas_model(1.0)
        .evaluate(&gas_source::point(450.0, 100000.0))
        .unwrap();
    with_wetting(&liquid, &interface, &outer, true, |thermal, gate| {
        let mut other =
            IsochoricGasParcel::new(gas_model(5000.0), gas_source::point(450.0, 100000.0), 1e-8)
                .unwrap();
        let flow = thermal.flow().lattice().grid().state_bits();
        let samples = thermal.wetting_samples().unwrap().clone();
        let energy = thermal.total_enthalpy_j().unwrap();
        let gas_before = other.state().clone();
        assert!(
            thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &mut other,
                    &[gas_face()],
                    gas_controls(),
                    false,
                    gate,
                )
                .is_err()
        );
        assert_eq!(other.state(), &gas_before);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), flow);
        assert_eq!(thermal.wetting_samples().unwrap(), &samples);
        assert_eq!(thermal.wetting_outer_gas(), Some(&outer));
        assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(thermal.step(&cancelled).is_err());
        assert_eq!(thermal.flow().lattice().grid().state_bits(), flow);
        assert_eq!(thermal.wetting_samples().unwrap(), &samples);
    });
}
