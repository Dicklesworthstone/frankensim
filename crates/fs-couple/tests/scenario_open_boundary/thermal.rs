//! Original fluid and thermal kernels with scenario-controlled reservoirs.
use super::*;
use fs_couple::scenario_open_boundary::thermal as coupled;
use fs_lbm::lattice3::open_boundaries::thermal::{OpenThermalDiffusion3, OpenThermalError3};
use fs_lbm::thermal::boundary::ThermalWall;
use fs_lbm::thermal::open::{
    ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenError3, ThermalOpenFace3,
    ThermalOpenField3, ThermalOpenLimits3,
};
use fs_lbm::thermal::tensor::{DiffusivityField, DiffusivityTensor};

fn heat_boundaries() -> ThermalOpenBoundaries3 {
    ThermalOpenBoundaries3 {
        faces: [
            ThermalOpenFace3::Periodic,
            ThermalOpenFace3::Periodic,
            ThermalOpenFace3::Periodic,
            ThermalOpenFace3::Periodic,
            ThermalOpenFace3::Reservoir { temperature: 0.65 },
            ThermalOpenFace3::Reservoir { temperature: 0.9 },
        ],
        internal_wall: ThermalWall::Adiabatic,
    }
}

fn limits() -> ThermalOpenLimits3 {
    ThermalOpenLimits3 {
        maximum_cells: 80,
        maximum_storage_bytes: 1 << 20,
        temperature_range: [0.0, 1.0],
        balance_tolerance: 1e-10,
    }
}

fn tensor() -> DiffusivityTensor<3> {
    DiffusivityTensor::from_lattice([[0.07, 0.005, 0.0], [0.005, 0.05, 0.002], [0.0, 0.002, 0.06]])
        .unwrap()
}

fn heat_total(field: &ThermalOpenField3) -> f64 {
    field.populations().iter().flat_map(|g| g.iter()).sum()
}

fn heat_balance(
    before: f64,
    field: &ThermalOpenField3,
    source: &[f64],
    report: ThermalOpenBalance3,
) {
    let change = heat_total(field) - before;
    let supplied: f64 =
        source.iter().sum::<f64>() + report.reservoirs.iter().map(|r| r.heat_in).sum::<f64>();
    assert!((change - supplied).abs() < 2e-12);
    assert!((change - report.stored_change).abs() < 2e-12);
    assert!(report.combined_defect.abs() < 2e-12);
}

fn uniform<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut actual = grid_with_collision::<L, Q>(cx, collision);
        let mut reference = grid_with_collision::<L, Q>(cx, collision);
        let initial: Vec<_> = (0..80).map(|i| 0.4 + 0.002 * (i % 7) as f64).collect();
        let mut heat = actual
            .prepare_thermal_field(cx, &initial, 0.83, heat_boundaries(), limits())
            .unwrap();
        let mut reference_heat = reference
            .prepare_thermal_field(cx, &initial, 0.83, heat_boundaries(), limits())
            .unwrap();
        let tensor = tensor();
        let tau = vec![0.87; 80];
        let source = vec![0.001; 80];
        for step in 0..4 {
            let time = 0.1 * (step + 1) as f64;
            let inlet = velocity([0.1, -0.05, 0.2 + time]);
            let outlet = pressure(100.0 + 200.0 * time);
            let frames = FrameTree::new();
            let targets = forced_targets(&inlet, &outlet, &frames, time, &scaling);
            reference.rebind_open_boundaries(targets, &gate).unwrap();
            let diffusion = if step % 2 == 0 {
                OpenThermalDiffusion3::Relaxation(Some(&tau))
            } else {
                OpenThermalDiffusion3::Tensor(DiffusivityField::Uniform(&tensor))
            };
            let expected = match diffusion {
                OpenThermalDiffusion3::Relaxation(tau) => reference.step_serial_with_thermal(
                    &mut reference_heat,
                    cx,
                    None,
                    &source,
                    tau,
                    [0.0; 3],
                ),
                OpenThermalDiffusion3::Tensor(tensor) => reference
                    .step_serial_with_thermal_tensor_diffusion(
                        &mut reference_heat,
                        cx,
                        None,
                        &source,
                        tensor,
                        [0.0; 3],
                    ),
            }
            .unwrap();
            let before = heat_total(&heat);
            let report = coupled::step_serial(
                &mut actual,
                &mut heat,
                controls(&inlet, &outlet),
                &frames,
                time,
                &scaling,
                cx,
                None,
                &source,
                diffusion,
                [0.0; 3],
            )
            .unwrap();
            assert_eq!(actual.populations(), reference.populations());
            assert_eq!(heat.populations(), reference_heat.populations());
            assert_eq!(report.transport.momentum, expected.momentum);
            assert_eq!(report.transport.heat, expected.heat);
            assert_eq!(actual.boundaries(), targets);
            assert_eq!(actual.completed_steps(), step + 1);
            assert_eq!(heat.completed_steps(), step + 1);
            assert!(report.samples[4].is_some() && report.samples[5].is_some());
            heat_balance(before, &heat, &source, report.transport.heat);
        }
    });
}

fn forced<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut actual = forced_grid::<L, Q>(cx, collision);
        let mut reference = forced_grid::<L, Q>(cx, collision);
        let initial: Vec<_> = (0..80).map(|i| 0.4 + 0.002 * (i % 7) as f64).collect();
        let mut heat = actual
            .prepare_thermal_field(cx, &initial, 0.83, heat_boundaries(), limits())
            .unwrap();
        let mut reference_heat = reference
            .prepare_thermal_field(cx, &initial, 0.83, heat_boundaries(), limits())
            .unwrap();
        let tensor = tensor();
        let source = vec![0.001; 80];
        for step in 0..4 {
            let time = 0.1 * (step + 1) as f64;
            let inlet = velocity([0.1, -0.05, 0.2 + time]);
            let outlet = pressure(100.0 + 200.0 * time);
            let frames = FrameTree::new();
            let targets = forced_targets(&inlet, &outlet, &frames, time, &scaling);
            reference.rebind_open_boundaries(targets, &gate).unwrap();
            let changing = step >= 2;
            let diffusion = if step % 2 == 0 {
                OpenThermalDiffusion3::Relaxation(None)
            } else {
                OpenThermalDiffusion3::Tensor(DiffusivityField::Uniform(&tensor))
            };
            let previous = actual.forces().to_vec();
            let endpoint = if changing {
                force_epoch(step + 1)
            } else {
                previous.clone()
            };
            let mut expected_buffer = endpoint.clone();
            let mut actual_buffer = endpoint.clone();
            let expected = match (changing, diffusion) {
                (false, OpenThermalDiffusion3::Relaxation(tau)) => reference
                    .step_serial_with_thermal(
                        &mut reference_heat,
                        cx,
                        None,
                        &source,
                        tau,
                        [0.0; 3],
                    ),
                (true, OpenThermalDiffusion3::Relaxation(tau)) => reference
                    .advance_serial_with_thermal(
                        &mut reference_heat,
                        cx,
                        None,
                        &mut expected_buffer,
                        &source,
                        tau,
                        [0.0; 3],
                    ),
                (false, OpenThermalDiffusion3::Tensor(tensor)) => reference
                    .step_serial_with_thermal_tensor_diffusion(
                        &mut reference_heat,
                        cx,
                        None,
                        &source,
                        tensor,
                        [0.0; 3],
                    ),
                (true, OpenThermalDiffusion3::Tensor(tensor)) => reference
                    .advance_serial_with_thermal_tensor_diffusion(
                        &mut reference_heat,
                        cx,
                        None,
                        &mut expected_buffer,
                        &source,
                        tensor,
                        [0.0; 3],
                    ),
            }
            .unwrap();
            let before = forced_totals(&actual);
            let thermal_before = heat_total(&heat);
            let report = coupled::step_forced_serial(
                &mut actual,
                &mut heat,
                controls(&inlet, &outlet),
                &frames,
                time,
                &scaling,
                cx,
                None,
                if changing {
                    Some(&mut actual_buffer)
                } else {
                    None
                },
                &source,
                diffusion,
                [0.0; 3],
            )
            .unwrap();
            assert_eq!(actual.populations(), reference.populations());
            assert_eq!(heat.populations(), reference_heat.populations());
            assert_eq!(report.transport.momentum, expected.momentum);
            assert_eq!(report.transport.heat, expected.heat);
            assert_eq!(actual.forces(), endpoint);
            assert_eq!(actual_buffer, expected_buffer);
            assert_eq!(actual.boundaries(), targets);
            assert_eq!(actual.completed_steps(), step + 1);
            assert_eq!(heat.completed_steps(), step + 1);
            let after = forced_totals(&actual);
            assert!((after.0 - before.0 - report.transport.momentum.mass_change).abs() < 3e-12);
            for axis in 0..3 {
                let impulse: f64 = previous
                    .iter()
                    .zip(&endpoint)
                    .map(|(a, b)| 0.5 * (a[axis] + b[axis]))
                    .sum();
                assert!(
                    (after.1[axis]
                        - before.1[axis]
                        - impulse
                        - report.transport.momentum.impulse[axis])
                        .abs()
                        < 3e-12
                );
            }
            heat_balance(thermal_before, &heat, &source, report.transport.heat);
        }
    });
}

#[test]
fn reservoirs_reach_original_scalar_and_tensor_heat_on_both_native_lattices() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let tau = scaling.flow().tau;
    uniform::<D3q19, 19>(CollisionModel3::Bgk { tau });
    uniform::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
    forced::<D3q19, 19>(CollisionModel3::Bgk { tau });
    forced::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn closing_heat_refusal_keeps_scenario_targets_and_force_buffer_for_exact_retry() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = forced_grid::<D3q19, 19>(
            cx,
            CollisionModel3::Bgk {
                tau: scaling.flow().tau,
            },
        );
        let mut thermal = flow
            .prepare_thermal_field(cx, &[0.5; 80], 0.83, heat_boundaries(), limits())
            .unwrap();
        let before = (
            flow.populations().to_vec(),
            flow.forces().to_vec(),
            flow.boundaries(),
            flow.last_wall_exchange(),
            thermal.populations().to_vec(),
            thermal.last_balance(),
        );
        let mut endpoint = force_epoch(1);
        let saved_endpoint = endpoint.clone();
        let inlet = velocity([0.1, 0.0, 0.4]);
        let outlet = pressure(200.0);
        let frames = FrameTree::new();
        assert!(matches!(
            coupled::step_forced_serial(
                &mut flow,
                &mut thermal,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                cx,
                None,
                Some(&mut endpoint),
                &[0.8; 80],
                OpenThermalDiffusion3::Relaxation(None),
                [0.0; 3]
            ),
            Err(coupled::Error::Native(OpenThermalError3::Thermal(
                ThermalOpenError3::Temperature { .. }
            )))
        ));
        assert_eq!(flow.populations(), before.0);
        assert_eq!(flow.forces(), before.1);
        assert_eq!(flow.boundaries(), before.2);
        assert_eq!(flow.last_wall_exchange(), before.3);
        assert_eq!(thermal.populations(), before.4);
        assert_eq!(thermal.last_balance(), before.5);
        assert_eq!(endpoint, saved_endpoint);
        assert_eq!(flow.completed_steps(), 0);
        assert_eq!(thermal.completed_steps(), 0);
        let mut reference = forced_grid::<D3q19, 19>(
            cx,
            CollisionModel3::Bgk {
                tau: scaling.flow().tau,
            },
        );
        let mut reference_heat = reference
            .prepare_thermal_field(cx, &[0.5; 80], 0.83, heat_boundaries(), limits())
            .unwrap();
        reference
            .rebind_open_boundaries(
                forced_targets(&inlet, &outlet, &frames, 0.5, &scaling),
                &gate,
            )
            .unwrap();
        let mut expected_endpoint = saved_endpoint;
        reference
            .advance_serial_with_thermal(
                &mut reference_heat,
                cx,
                None,
                &mut expected_endpoint,
                &[0.001; 80],
                None,
                [0.0; 3],
            )
            .unwrap();
        coupled::step_forced_serial(
            &mut flow,
            &mut thermal,
            controls(&inlet, &outlet),
            &frames,
            0.5,
            &scaling,
            cx,
            None,
            Some(&mut endpoint),
            &[0.001; 80],
            OpenThermalDiffusion3::Relaxation(None),
            [0.0; 3],
        )
        .unwrap();
        assert_eq!(flow.populations(), reference.populations());
        assert_eq!(thermal.populations(), reference_heat.populations());
        assert_eq!(flow.forces(), reference.forces());
        assert_eq!(flow.boundaries(), reference.boundaries());
        assert_eq!(endpoint, expected_endpoint);
    });
}

#[test]
fn heat_coupling_retains_selected_profile_refusal_and_unselected_trace() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = forced_grid::<D3q19, 19>(
            cx,
            CollisionModel3::Bgk {
                tau: scaling.flow().tau,
            },
        );
        flow.set_velocity_profile(Face3::ZMin, InletFlux3::Unconstrained, &gate, |p| {
            Ok([0.002, -0.001, 0.01 + 0.001 * p[0] as f64])
        })
        .unwrap();
        let profile = flow.velocity_profile(Face3::ZMin).unwrap().to_vec();
        let mut thermal = flow
            .prepare_thermal_field(cx, &[0.5; 80], 0.83, heat_boundaries(), limits())
            .unwrap();
        let before = (
            flow.populations().to_vec(),
            flow.forces().to_vec(),
            flow.boundaries(),
            thermal.populations().to_vec(),
        );
        let inlet = velocity([0.0; 3]);
        let outlet = pressure(100.0);
        let frames = FrameTree::new();
        let mut endpoint = force_epoch(1);
        let saved_endpoint = endpoint.clone();
        assert!(matches!(
            coupled::step_forced_serial(
                &mut flow,
                &mut thermal,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                cx,
                None,
                Some(&mut endpoint),
                &[0.0; 80],
                OpenThermalDiffusion3::Relaxation(None),
                [0.0; 3]
            ),
            Err(coupled::Error::Scenario(_))
        ));
        assert_eq!(flow.populations(), before.0);
        assert_eq!(flow.forces(), before.1);
        assert_eq!(flow.boundaries(), before.2);
        assert_eq!(thermal.populations(), before.3);
        assert_eq!(endpoint, saved_endpoint);
        assert_eq!(thermal.completed_steps(), 0);
        let mut only_outlet = controls(&inlet, &outlet);
        only_outlet[4] = None;
        let report = coupled::step_forced_serial(
            &mut flow,
            &mut thermal,
            only_outlet,
            &frames,
            0.5,
            &scaling,
            cx,
            None,
            None,
            &[0.0; 80],
            OpenThermalDiffusion3::Relaxation(None),
            [0.0; 3],
        )
        .unwrap();
        assert_eq!(report.samples[4], None);
        assert_eq!(flow.velocity_profile(Face3::ZMin), Some(profile.as_slice()));
        assert_eq!(flow.forces(), before.1);
        for (cell, target) in profile.iter().enumerate() {
            for (value, expected) in flow.macroscopic(cell).unwrap().1.iter().zip(target) {
                near(*value, *expected);
            }
        }
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(thermal.completed_steps(), 1);
    });
}
