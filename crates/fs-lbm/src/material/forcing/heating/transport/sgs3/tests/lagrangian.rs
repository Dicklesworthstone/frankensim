//! G0/G1/G4: actual material epochs and physical publication with joint memory.
use super::*;
use crate::sgs3::dynamic::lagrangian::scalar::LagrangianScalarConfig;
use crate::sgs3::dynamic::lagrangian::{AdaptiveLagrangianMemory, LagrangianSgsConfig};
use crate::thermal::sources::sgs3::lagrangian::{
    ThermalLagrangianSgsConfig, ThermalLagrangianSgsReport, ThermalLagrangianSgsWorkspace3,
};

fn config() -> ThermalLagrangianSgsConfig {
    ThermalLagrangianSgsConfig {
        momentum: LagrangianSgsConfig {
            memory_time: 4.0,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 10.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        heat: LagrangianScalarConfig {
            memory_time: 2.0,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 10.0,
            maximum_diffusivity_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_mach: 0.2,
    }
}

#[allow(clippy::too_many_arguments)]
fn run(
    state: &mut MaterialSlab<ThermalLbm3>,
    workspace: &mut ThermalLagrangianSgsWorkspace3,
    time: f64,
    mode: usize,
    gate: &CancelGate,
    source: impl FnMut(usize, f64, &FluidTransportState) -> Result<FluidTransportState, MaterialError>,
) -> Result<(HeatBalance, ThermalLagrangianSgsReport), MaterialTransportError> {
    let (power, forces) = loads();
    match mode {
        0 => state
            .step_with_material_transport_and_lagrangian_sgs(
                &power,
                &forces,
                limits(),
                time,
                config(),
                workspace,
                gate,
                source,
            )
            .map(|r| (r.material.heat, r.sgs)),
        1 => state
            .step_with_boussinesq_material_transport_and_lagrangian_sgs(
                &power,
                &forces,
                limits(),
                buoyancy(),
                time,
                config(),
                workspace,
                gate,
                source,
            )
            .map(|r| (r.material.heat, r.sgs)),
        _ => state
            .step_with_boussinesq_material_transport_and_lagrangian_sgs_endpoint_averaged(
                &power,
                &forces,
                limits(),
                buoyancy(),
                time,
                config(),
                workspace,
                gate,
                source,
            )
            .map(|r| (r.material.transport.heat, r.sgs)),
    }
}

#[test]
fn g0_zero_histories_preserve_all_original_material_source_transactions() {
    let gate = CancelGate::new_clock_free();
    let (power, forces) = loads();
    for mode in 0..3 {
        let mut actual = endpoint_slab();
        let mut expected = endpoint_slab();
        for state in [&mut actual, &mut expected] {
            for cell in PLANE..CELLS - PLANE {
                state.lattice.grid.f[cell] = equilibrium3(1.0, [0.0; 3]);
                state.lattice.grid.fext[cell] = [0.0; 3];
            }
        }
        let mut work = actual
            .lattice
            .thermal_lagrangian_sgs_workspace(usize::MAX)
            .unwrap();
        let (heat, report) = run(&mut actual, &mut work, 0.0, mode, &gate, |_, t, _| {
            Ok(varying(t))
        })
        .unwrap();
        let baseline = match mode {
            0 => {
                expected
                    .step_with_material_transport(&power, &forces, limits(), &gate, |_, t, _| {
                        Ok(varying(t))
                    })
                    .unwrap()
                    .heat
            }
            1 => {
                expected
                    .step_with_boussinesq_material_transport(
                        &power,
                        &forces,
                        limits(),
                        buoyancy(),
                        &gate,
                        |_, t, _| Ok(varying(t)),
                    )
                    .unwrap()
                    .heat
            }
            _ => {
                expected
                    .step_with_boussinesq_material_transport_endpoint_averaged(
                        &power,
                        &forces,
                        limits(),
                        buoyancy(),
                        &gate,
                        |_, t, _| Ok(varying(t)),
                    )
                    .unwrap()
                    .transport
                    .heat
            }
        };
        assert_eq!(heat, baseline);
        assert_eq!(report.transport.momentum.maximum_eddy_viscosity_lu, 0.0);
        assert_eq!(report.transport.maximum_eddy_diffusivity_lu, 0.0);
        assert_eq!(bits(&actual), bits(&expected));
        assert_eq!(
            actual.lattice.last_wall_exchange(),
            expected.lattice.last_wall_exchange()
        );
        assert_eq!(work.accepted_time(), Some(0.0));
        for cell in PLANE..CELLS - PLANE {
            assert_eq!(
                actual.lattice.grid.collision_model_at(cell),
                expected.lattice.grid.collision_model_at(cell)
            );
        }
    }
}

#[test]
fn g1_joint_history_heat_drives_original_endpoint_density_and_source_epochs() {
    let gate = CancelGate::new_clock_free();
    let (power, external) = loads();
    let mut actual = endpoint_slab();
    let original_molecular = actual.lattice.grid.tau.clone();
    let original_thermal = actual.lattice.tau_g;
    let mut work = actual
        .lattice
        .thermal_lagrangian_sgs_workspace(usize::MAX)
        .unwrap()
        .with_adaptive_momentum_memory(AdaptiveLagrangianMemory {
            coefficient: 1.5,
            minimum_time: 0.2,
        })
        .unwrap();
    let pointer = work.relaxation().thermal.as_ptr();
    for step in 0..3 {
        let reads = Rc::new(Cell::new(0));
        let result = actual
            .step_with_boussinesq_material_transport_and_lagrangian_sgs_endpoint_averaged(
                &power,
                &external,
                limits(),
                buoyancy(),
                step as f64,
                config(),
                &mut work,
                &gate,
                |cell, t, _| {
                    Ok(Ticket {
                        value: nonlinear_density(t),
                        query_k: t,
                        cell,
                        reads: reads.clone(),
                    })
                },
            )
            .unwrap();
        assert_eq!(reads.get(), 3 * FLUID);
        assert_eq!(work.accepted_time(), Some(step as f64));
        assert_eq!(
            work.accepted_momentum_memory_time_range(),
            Some(result.sgs.momentum_memory_time_range)
        );
        assert!(result.sgs.transport.maximum_eddy_diffusivity_lu > 0.0);
        let material = &result.material;
        for sample in material
            .transport
            .applied
            .iter()
            .chain(&material.force_endpoint)
            .chain(&material.transport.endpoint)
        {
            assert_eq!(sample.material.material().query_k, sample.temperature_k);
            assert_eq!(sample.material.material().cell, sample.cell);
            assert_eq!(
                *sample.material.resolved_transport(),
                nonlinear_density(sample.temperature_k)
            );
            assert_eq!(
                sample.transport.density_kg_m3(),
                reference().density_kg_m3()
            );
            assert_eq!(
                sample.transport.specific_heat_j_kg_k(),
                reference().specific_heat_j_kg_k()
            );
        }
        for (opening, predicted) in material
            .transport
            .applied
            .iter()
            .zip(&material.force_endpoint)
        {
            let d0 = opening.material.resolved_transport().density_kg_m3() - 2.0;
            let d1 = predicted.material.resolved_transport().density_kg_m3() - 2.0;
            for axis in 0..3 {
                let force = external[predicted.cell][axis]
                    + (0.5 * d0 + 0.5 * d1) * buoyancy().gravity_m_s2[axis];
                assert_eq!(
                    actual.lattice.grid.fext[predicted.cell][axis],
                    actual.force_units.component(force).unwrap()
                );
            }
        }
        for sample in &material.transport.endpoint {
            assert_eq!(
                sample.temperature_k,
                actual
                    .temperature_si(actual.lattice.gpop[sample.cell].iter().sum())
                    .unwrap()
            );
        }
        assert!(material.transport.heat.defect_j.abs() <= limits().balance_tolerance_j);
        assert_eq!(actual.lattice.grid.tau, original_molecular);
        assert_eq!(actual.lattice.tau_g, original_thermal);
        assert_eq!(pointer, work.relaxation().thermal.as_ptr());
    }
}

#[test]
fn g4_late_material_refusal_and_cancellation_preserve_warmed_joint_history() {
    for mode in 0..3 {
        for refusal in 0..3 {
            let mut actual = endpoint_slab();
            let mut expected = endpoint_slab();
            let mut work = actual
                .lattice
                .thermal_lagrangian_sgs_workspace(usize::MAX)
                .unwrap();
            let mut clean = expected
                .lattice
                .thermal_lagrangian_sgs_workspace(usize::MAX)
                .unwrap();
            let gate = CancelGate::new_clock_free();
            assert_eq!(
                run(&mut actual, &mut work, 0.0, mode, &gate, |_, t, _| Ok(
                    varying(t)
                ))
                .unwrap(),
                run(&mut expected, &mut clean, 0.0, mode, &gate, |_, t, _| Ok(
                    varying(t)
                ))
                .unwrap()
            );
            let before = bits(&actual);
            let old_heat = actual.lattice.last_wall_exchange();
            let old_time_range = work.accepted_momentum_memory_time_range();
            let epochs = if mode == 2 { 3 } else { 2 };
            let mut calls = 0;
            let result = run(&mut actual, &mut work, 1.0, mode, &gate, |_, t, _| {
                calls += 1;
                if calls == (epochs - 1) * FLUID + 1 {
                    if refusal == 0 {
                        return Err(MaterialError::State {
                            what: "late original source refusal".into(),
                        });
                    }
                    if refusal == 1 {
                        return Ok(FluidTransportState::try_new(2.3, 3.0, 0.02, 0.01).unwrap());
                    }
                }
                if calls == epochs * FLUID && refusal == 2 {
                    gate.request();
                }
                Ok(varying(t))
            });
            assert!(result.is_err());
            assert!(calls > (epochs - 1) * FLUID);
            assert_eq!(bits(&actual), before);
            assert_eq!(actual.lattice.last_wall_exchange(), old_heat);
            assert_eq!(work.accepted_time(), Some(0.0));
            assert_eq!(work.accepted_momentum_memory_time_range(), old_time_range);
            let retry_gate = CancelGate::new_clock_free();
            for step in 1..=2 {
                let replay = run(
                    &mut actual,
                    &mut work,
                    step as f64,
                    mode,
                    &retry_gate,
                    |_, t, _| Ok(varying(t)),
                )
                .unwrap();
                let clean_result = run(
                    &mut expected,
                    &mut clean,
                    step as f64,
                    mode,
                    &retry_gate,
                    |_, t, _| Ok(varying(t)),
                )
                .unwrap();
                assert_eq!(replay, clean_result);
                assert_eq!(bits(&actual), bits(&expected));
                assert_eq!(work.accepted_time(), Some(step as f64));
            }
        }
    }
}
