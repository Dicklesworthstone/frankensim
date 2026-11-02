//! G0/G1 transport identities and physical balance checks on the actual slab;
//! these small fixtures are not turbulent convection/DNS validation.

use super::*;
use crate::d3q19::equilibrium3;
use crate::sgs3::{SgsModel, SgsUse};

fn config(model: SgsModel) -> ThermalSgsConfig {
    ThermalSgsConfig {
        closure: SgsConfig {
            model,
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.7,
        maximum_mach: 0.2,
    }
}

fn limits(cells: usize) -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: cells,
        temperature_range: [-0.2, 1.2],
        balance_tolerance: 2e-12,
    }
}

fn bits(s: &ThermalLbm3) -> Vec<u64> {
    s.grid
        .f
        .iter()
        .flatten()
        .chain(s.gpop.iter().flatten())
        .chain(s.grid.fext.iter().flatten())
        .chain(&s.grid.tau)
        .chain([&s.tau_g, &s.t_bottom, &s.t_top, &s.gbeta, &s.t_ref])
        .map(|v| v.to_bits())
        .collect()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance}"
    );
}

fn physical_velocity(s: &mut ThermalLbm3, cell: usize, u: [f64; 3]) {
    let force = [1e-5 * (cell % 3) as f64, -2e-5 * (cell % 5) as f64, 3e-5];
    s.grid.f[cell] = equilibrium3(1.0, u);
    crate::thermal::forcing::reframe3(&mut s.grid.f[cell], u, [0.0; 3], force);
    s.grid.fext[cell] = force;
}

fn seed(s: &mut ThermalLbm3) {
    let (nx, ny, height) = (s.grid.nx, s.grid.ny, s.grid.nz - 2);
    for z in 1..s.grid.nz - 1 {
        let h = z as f64 - 0.5;
        let wall = 4.0 * h * (height as f64 - h) / (height * height) as f64;
        for y in 0..ny {
            for x in 0..nx {
                let sx = fs_math::det::sin(core::f64::consts::TAU * x as f64 / nx as f64);
                let sy = fs_math::det::sin(core::f64::consts::TAU * y as f64 / ny as f64);
                let cell = s.grid.idx(x, y, z);
                physical_velocity(
                    s,
                    cell,
                    [
                        wall * (0.01 + 0.003 * sy),
                        wall * 0.004 * sx,
                        wall * 0.002 * (sx - sy),
                    ],
                );
                s.gpop[cell][0] += 0.005 * sx;
            }
        }
    }
}

#[test]
fn zero_closures_preserve_sources_forces_and_local_molecular_fields_bitwise() {
    let gate = CancelGate::new_clock_free();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        for local in [false, true] {
            let mut reference = ThermalLbm3::slab(4, 3, 6, 0.8, 0.9, 1e-5);
            seed(&mut reference);
            let mut actual = reference.snapshot().0;
            let n = actual.gpop.len();
            let mut workspace = actual.thermal_sgs_workspace(n).unwrap();
            let pointer = workspace.relaxation().momentum.as_ptr();
            let momentum: Vec<_> = (0..n).map(|i| 0.7 + 0.01 * (i % 5) as f64).collect();
            let thermal: Vec<_> = (0..n).map(|i| 0.8 + 0.01 * (i % 7) as f64).collect();
            let molecular = local.then_some(RelaxationField {
                momentum: &momentum,
                thermal: &thermal,
            });
            let mut heat = vec![0.0001; n];
            heat[..12].fill(0.0);
            heat[n - 12..].fill(0.0);
            for force in [[1e-5, 2e-5, -1e-5], [-2e-5, 1e-5, 0.0]] {
                let forces = vec![force; n];
                let expected = match molecular {
                    Some(field) => {
                        reference.step_with_sources_and_transport(&heat, &forces, field, limits(n))
                    }
                    None => reference.step_with_sources(&heat, &forces, limits(n)),
                }
                .unwrap();
                let observed = actual
                    .step_with_sources_and_sgs(
                        &heat,
                        &forces,
                        config(model),
                        molecular,
                        &mut workspace,
                        limits(n),
                        &gate,
                    )
                    .unwrap();
                assert_eq!(observed.heat, expected);
                assert_eq!(observed.sgs.momentum.sgs_dissipation_lu, 0.0);
                assert_eq!(observed.sgs.maximum_eddy_diffusivity_lu, 0.0);
                assert_eq!(bits(&actual), bits(&reference));
                assert_eq!(pointer, workspace.relaxation().momentum.as_ptr());
            }
        }
    }
}

#[test]
fn quadratic_halfway_wall_shear_drives_the_independent_momentum_and_heat_rates() {
    let gate = CancelGate::new_clock_free();
    let mut state = ThermalLbm3::slab(3, 2, 8, 0.8, 0.9, 0.0);
    let n = state.gpop.len();
    let mut molecular_momentum = vec![0.8; n];
    let mut molecular_thermal = vec![0.9; n];
    let mut total_momentum = molecular_momentum.clone();
    let mut total_thermal = molecular_thermal.clone();
    let settings = config(SgsModel::Smagorinsky { cs: 0.23 });
    let mut expected_dissipation = 0.0;
    let mut expected_resolved = 0.0;
    for z in 1..state.grid.nz - 1 {
        let h = z as f64 - 0.5;
        let shear = 0.0005 * (8.0 - 2.0 * h);
        let eddy = 0.23_f64.powi(2) * shear.abs();
        for y in 0..2 {
            for x in 0..3 {
                let cell = state.grid.idx(x, y, z);
                physical_velocity(&mut state, cell, [0.0005 * h * (8.0 - h), 0.0, 0.0]);
                state.gpop[cell][0] += 0.01 * ((x + 2 * y) % 3) as f64;
                let nu = 0.05 + 0.001 * z as f64;
                let alpha = 0.06 + 0.002 * z as f64;
                molecular_momentum[cell] = 0.5 + 3.0 * nu;
                molecular_thermal[cell] = 0.5 + 4.0 * alpha;
                total_momentum[cell] = 0.5 + 3.0 * (nu + eddy);
                total_thermal[cell] = 0.5 + 4.0 * (alpha + eddy / settings.turbulent_prandtl);
                expected_dissipation += eddy * shear * shear;
                expected_resolved += nu * shear * shear;
            }
        }
    }
    let original_tau = state.grid.tau.clone();
    let mut oracle = state.snapshot().0;
    let mut molecular_only = state.snapshot().0;
    let mut workspace = state.thermal_sgs_workspace(n).unwrap();
    let forces = vec![[2e-5, 0.0, 0.0]; n];
    let heat = vec![0.0; n];
    let expected = oracle
        .step_with_sources_and_transport(
            &heat,
            &forces,
            RelaxationField {
                momentum: &total_momentum,
                thermal: &total_thermal,
            },
            limits(n),
        )
        .unwrap();
    molecular_only
        .step_with_sources_and_transport(
            &heat,
            &forces,
            RelaxationField {
                momentum: &molecular_momentum,
                thermal: &molecular_thermal,
            },
            limits(n),
        )
        .unwrap();
    let actual = state
        .step_with_sources_and_sgs(
            &heat,
            &forces,
            settings,
            Some(RelaxationField {
                momentum: &molecular_momentum,
                thermal: &molecular_thermal,
            }),
            &mut workspace,
            limits(n),
            &gate,
        )
        .unwrap();
    for ((&tau, &expected), (&thermal, &expected_heat)) in workspace
        .relaxation()
        .momentum
        .iter()
        .zip(&total_momentum)
        .zip(workspace.relaxation().thermal.iter().zip(&total_thermal))
    {
        near(tau, expected, 3e-15);
        near(thermal, expected_heat, 3e-15);
    }
    near(
        actual.sgs.momentum.sgs_dissipation_lu,
        expected_dissipation,
        1e-18,
    );
    near(
        actual.sgs.momentum.resolved_dissipation_lu,
        expected_resolved,
        1e-16,
    );
    for (a, b) in state
        .grid
        .f
        .iter()
        .flatten()
        .zip(oracle.grid.f.iter().flatten())
    {
        near(*a, *b, 2e-14);
    }
    for (a, b) in state
        .gpop
        .iter()
        .flatten()
        .zip(oracle.gpop.iter().flatten())
    {
        near(*a, *b, 2e-14);
    }
    near(actual.heat.stored_change, expected.stored_change, 1e-12);
    near(
        actual.heat.stored_change,
        actual.heat.bottom_wall + actual.heat.top_wall,
        2e-12,
    );
    assert_eq!(actual.heat.source, 0.0); // No modeled dissipation relabeled as heat.
    assert_ne!(state.grid.f, molecular_only.grid.f);
    assert!(
        state
            .gpop
            .iter()
            .zip(&molecular_only.gpop)
            .any(|(a, b)| (a.iter().sum::<f64>() - b.iter().sum::<f64>()).abs() > 1e-9)
    );
    assert_eq!(state.grid.tau, original_tau);
    assert_eq!(state.tau_g, 0.9);
}

#[test]
fn all_models_use_one_eddy_field_for_flow_and_declared_turbulent_prandtl() {
    let gate = CancelGate::new_clock_free();
    for model in [
        SgsModel::Smagorinsky { cs: 0.2 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let mut a = ThermalLbm3::slab(4, 3, 6, 0.8, 0.9, 1e-5);
        seed(&mut a);
        let mut b = a.snapshot().0;
        let n = a.gpop.len();
        let mut wa = a.thermal_sgs_workspace(n).unwrap();
        let mut wb = b.thermal_sgs_workspace(n).unwrap();
        let mut ca = config(model);
        ca.turbulent_prandtl = 0.5;
        let mut cb = config(model);
        cb.turbulent_prandtl = 2.0;
        let ra = a
            .step_with_sources_and_sgs(
                &vec![0.0; n],
                &vec![[0.0; 3]; n],
                ca,
                None,
                &mut wa,
                limits(n),
                &gate,
            )
            .unwrap();
        let rb = b
            .step_with_sources_and_sgs(
                &vec![0.0; n],
                &vec![[0.0; 3]; n],
                cb,
                None,
                &mut wb,
                limits(n),
                &gate,
            )
            .unwrap();
        assert_eq!(a.grid.f, b.grid.f);
        assert_eq!(wa.relaxation().momentum, wb.relaxation().momentum);
        assert!(ra.sgs.momentum.maximum_eddy_viscosity_lu > 0.0);
        near(
            ra.sgs.maximum_eddy_diffusivity_lu,
            4.0 * rb.sgs.maximum_eddy_diffusivity_lu,
            1e-18,
        );
        for cell in 12..n - 12 {
            let nu_t = (wa.momentum[cell] - 0.8) / 3.0;
            near((wa.thermal[cell] - 0.9) / 4.0, nu_t / 0.5, 2e-16);
            near((wb.thermal[cell] - 0.9) / 4.0, nu_t / 2.0, 2e-16);
        }
        assert_ne!(a.gpop, b.gpop);
        for balance in [ra.heat, rb.heat] {
            near(
                balance.stored_change,
                balance.bottom_wall + balance.top_wall,
                2e-12,
            );
            assert_eq!(balance.source, 0.0);
        }
    }
}

#[test]
fn failed_preparation_source_and_endpoint_admission_preserve_all_physical_state() {
    let gate = CancelGate::new_clock_free();
    let mut state = ThermalLbm3::slab(4, 3, 6, 0.8, 0.9, 0.0);
    seed(&mut state);
    let n = state.gpop.len();
    let before = bits(&state);
    let mut workspace = state.thermal_sgs_workspace(n).unwrap();
    let settings = config(SgsModel::Smagorinsky { cs: 0.2 });
    let forces = vec![[0.0; 3]; n];
    let mut heat = vec![0.0; n];
    for pr in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut bad = settings;
        bad.turbulent_prandtl = pr;
        assert!(
            state
                .step_with_sources_and_sgs(
                    &heat,
                    &forces,
                    bad,
                    None,
                    &mut workspace,
                    limits(n),
                    &gate
                )
                .is_err()
        );
        assert_eq!(bits(&state), before);
    }
    let mut restricted = settings;
    restricted.closure.maximum_eddy_to_molecular = 0.0;
    assert!(matches!(
        state.step_with_sources_and_sgs(
            &heat,
            &forces,
            restricted,
            None,
            &mut workspace,
            limits(n),
            &gate
        ),
        Err(ThermalSgsError::Closure(
            SgsError::ViscosityRatioExceeded { .. }
        ))
    ));
    assert_eq!(bits(&state), before);
    heat[20] = 100.0;
    assert!(matches!(
        state.step_with_sources_and_sgs(
            &heat,
            &forces,
            settings,
            None,
            &mut workspace,
            limits(n),
            &gate
        ),
        Err(ThermalSgsError::Source(
            ThermalSourceError::Temperature { .. }
        ))
    ));
    assert_eq!(bits(&state), before);
    gate.request();
    assert_eq!(
        state.prepare_thermal_sgs(settings, None, &mut workspace, &gate),
        Err(ThermalSgsError::Cancelled)
    );
    assert_eq!(bits(&state), before);

    // Incoming state passes the Mach gate; the actual forced endpoint fails.
    let gate = CancelGate::new_clock_free();
    let mut state = ThermalLbm3::slab(2, 2, 4, 0.8, 0.9, 0.0);
    let n = state.gpop.len();
    let before = bits(&state);
    let mut workspace = state.thermal_sgs_workspace(n).unwrap();
    let mut strict = settings;
    strict.maximum_mach = 0.001;
    assert!(matches!(
        state.step_with_sources_and_sgs(
            &vec![0.0; n],
            &vec![[0.002, 0.0, 0.0]; n],
            strict,
            None,
            &mut workspace,
            limits(n),
            &gate
        ),
        Err(ThermalSgsError::Mach { .. })
    ));
    assert_eq!(bits(&state), before);
}

#[test]
fn retained_material_adapter_preserves_si_heat_with_both_temperature_orientations() {
    use crate::material::{MaterialThermalSlab3, plan_scaling_from_material};
    use fs_material::fluid_transport::FluidTransportState;
    let gate = CancelGate::new_clock_free();
    let material = FluidTransportState::try_new(2.0, 3.0, 0.02, 0.01).unwrap();
    let scale = plan_scaling_from_material(&material, 0.01, 1e-4, 0.1).unwrap();
    for (bottom, top) in [(310.0, 300.0), (300.0, 310.0)] {
        for cs in [0.0, 0.2] {
            let boundary = MaterialThermalSlab3 {
                nx: 4,
                ny: 3,
                nz: 6,
                bottom_temperature_k: bottom,
                top_temperature_k: top,
                gravity_m_s2: 0.0,
                expansion_per_k: 0.001,
                max_cells: 96,
            };
            let mut slab = scale.thermal_slab3(boundary).unwrap();
            let mut reference = scale.thermal_slab3(boundary).unwrap();
            seed(&mut slab.lattice);
            seed(&mut reference.lattice);
            let n = slab.lattice.gpop.len();
            let mut workspace = slab.lattice.thermal_sgs_workspace(n).unwrap();
            let mut power = vec![100.0; n];
            power[..12].fill(0.0);
            power[n - 12..].fill(0.0);
            let forces = vec![[0.1, 0.0, 0.0]; n];
            let limits = super::super::HeatSourceLimits {
                max_cells: n,
                temperature_range_k: [290.0, 320.0],
                balance_tolerance_j: 1e-12,
            };
            let expected = reference
                .step_with_power_density(&power, &forces, limits)
                .unwrap();
            let result = slab
                .step_with_power_density_and_sgs(
                    &power,
                    &forces,
                    limits,
                    config(SgsModel::Smagorinsky { cs }),
                    &mut workspace,
                    &gate,
                )
                .unwrap();
            near(
                result.heat.source_j,
                72.0 * 100.0 * 0.01_f64.powi(3) * 1e-4,
                1e-20,
            );
            near(
                result.heat.stored_change_j - result.heat.source_j,
                result.heat.bottom_wall_j + result.heat.top_wall_j,
                1e-12,
            );
            if cs == 0.0 {
                assert_eq!(result.heat, expected);
                assert_eq!(bits(&slab.lattice), bits(&reference.lattice));
            } else {
                assert!(result.sgs.momentum.sgs_dissipation_lu > 0.0);
                assert_ne!(slab.lattice.gpop, reference.lattice.gpop);
            }
        }
    }
}
