//! G1/G3: prescribed wall work drives native Couette flow and source-owned heat.
//! All material parameters here are declared synthetic calibration inputs.

use fs_exec::CancelGate;
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::material::transport::rheology::RheologySolve;
use fs_lbm::material::{
    MaterialSlab3, MaterialThermalSlab3, MaterialThermalWall, plan_scaling_from_material,
};
use fs_lbm::thermal::ThermalForcingError;
use fs_lbm::thermal::sources::{HeatSourceLimits, ThermalSourceError};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use std::cell::Cell;

const DX: f64 = 0.01;
const DT: f64 = 0.02;
const RHO: f64 = 1000.0;
const MU: f64 = 0.5;
const SPEED: f64 = 0.03;
const DRIVE: [[f64; 3]; 2] = [[0.0; 3], [0.024, 0.018, 0.0]];

fn material(cp: f64) -> FluidTransportState {
    // Both diffusivities are 0.1 in lattice units; varying Cp also varies k
    // explicitly, so a thermal feedback comparison does not change diffusion.
    FluidTransportState::try_new(RHO, cp, 0.5 * cp, MU).unwrap()
}

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 100,
        relaxation: 0.75,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-12,
    }
}

fn limits(height: usize) -> HeatSourceLimits {
    HeatSourceLimits {
        max_cells: height + 2,
        temperature_range_k: [299.0, 302.0],
        balance_tolerance_j: 1e-14,
    }
}

fn slab(height: usize, source: FluidTransportState) -> MaterialSlab3 {
    let mut slab = plan_scaling_from_material(&source, DX, DT, SPEED)
        .unwrap()
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 1,
            ny: 1,
            nz: height,
            bottom_temperature_k: 301.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: height + 2,
        })
        .unwrap();
    // Keep the original signed temperature map, but prescribe two equal wall
    // temperatures. Starting from uniform 300 K isolates viscous self-heating.
    slab.set_thermal_walls([MaterialThermalWall::Temperature(300.0); 2])
        .unwrap();
    slab.lattice.gpop.fill([0.0; 7]);
    slab.lattice.reset_advection_history();
    slab.lattice.grid.set_trt_magic(3.0 / 16.0).unwrap();
    slab
}

fn drive(slab: &mut MaterialSlab3) {
    let center = 0.5 * (slab.lattice.grid.nz as f64 - 1.0) * DX;
    slab.set_plate_motion_si(DRIVE, [0.0, 0.0, center]).unwrap();
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}

fn bits<const Q: usize>(field: &[[f64; Q]]) -> Vec<u64> {
    field.iter().flatten().map(|v| v.to_bits()).collect()
}

fn steady_couette(height: usize) -> f64 {
    let source = material(1.0);
    let mut slab = slab(height, source);
    drive(&mut slab);
    let heat = vec![0.0; height + 2];
    let forces = vec![[0.0; 3]; height + 2];
    let gate = CancelGate::new_clock_free();
    let mass_before = slab.lattice.grid.total_mass();
    let mut last = None;
    // 50 H² steps exceed 49 slowest-mode diffusion e-foldings for nu=alpha=0.1.
    for _ in 0..50 * height * height {
        last = Some(
            slab.step_with_material_rheology_and_viscous_heating(
                &heat,
                &forces,
                limits(height),
                solve(),
                &gate,
                |_, _| Ok(source),
            )
            .unwrap(),
        );
    }
    let step = last.unwrap();
    near(
        slab.lattice.grid.total_mass(),
        mass_before,
        2e-10 * mass_before,
    );
    let h = height as f64 * DX;
    let rate = SPEED / h;
    let power_density = MU * rate * rate;
    let expected_heat = power_density * (height as f64 * DX.powi(3)) * DT;
    // The kinetic stress estimate has finite-grid corrections; the wall work
    // and analytic shear are independent of that estimator.
    near(
        step.viscous.viscous_heat_j,
        expected_heat,
        2e-3 * expected_heat,
    );
    near(
        step.transport.heat.source_j,
        step.viscous.viscous_heat_j,
        1e-20,
    );
    assert_eq!(step.viscous.authored_heat_j, 0.0);
    near(
        step.transport.heat.stored_change_j,
        0.0,
        1e-8 * expected_heat,
    );
    near(
        step.transport.heat.bottom_wall_j + step.transport.heat.top_wall_j,
        -step.viscous.viscous_heat_j,
        1e-8 * expected_heat,
    );
    let native = slab.lattice.grid.last_plate_exchange().unwrap();
    let physical = slab.last_plate_exchange_si().unwrap().unwrap();
    let mass_unit = RHO * DX.powi(3);
    let impulse_unit = mass_unit * DX / DT;
    let work_unit = impulse_unit * DX / DT;
    for side in 0..2 {
        assert_eq!(physical[side].measured_links, 5);
        near(
            physical[side].wall_work_j,
            native[side].wall_work * work_unit,
            1e-22,
        );
        near(physical[side].fluid_mass_change_kg, 0.0, 1e-15);
        for axis in 0..3 {
            near(
                physical[side].wall_impulse_n_s[axis],
                native[side].wall_impulse[axis] * impulse_unit,
                1e-18,
            );
            near(
                physical[side].wall_impulse_n_s[axis]
                    + physical[side].fluid_population_impulse_n_s[axis],
                physical[side].wall_velocity_mass_impulse_n_s[axis],
                1e-16,
            );
            near(
                physical[side].wall_angular_impulse_n_m_s[axis]
                    + physical[side].fluid_population_angular_impulse_n_m_s[axis],
                physical[side].wall_velocity_mass_angular_impulse_n_m_s[axis],
                1e-18,
            );
        }
    }
    near(physical[0].wall_work_j, 0.0, 1e-22);
    near(
        physical[1].wall_work_j,
        -expected_heat,
        2e-3 * expected_heat,
    );
    let peak = MU * SPEED * SPEED / (8.0 * source.thermal_conductivity_w_m_k());
    let mut thermal_error = 0.0_f64;
    for z in 1..=height {
        let fraction = (z as f64 - 0.5) / height as f64;
        let velocity = slab.lattice.grid.moments3(slab.lattice.grid.idx(0, 0, z)).1;
        for axis in 0..3 {
            near(
                velocity[axis] * DX / DT,
                DRIVE[1][axis] * fraction,
                5e-7 * SPEED,
            );
        }
        let expected = 300.0 + 4.0 * peak * fraction * (1.0 - fraction);
        let actual = slab
            .temperature_si(slab.lattice.temperature(0, 0, z))
            .unwrap();
        thermal_error = thermal_error.max((actual - expected).abs() / peak);
    }
    thermal_error
}

#[test]
fn couette_material_heat_and_wall_work_match_the_analytic_solution_under_refinement() {
    let coarse = steady_couette(6);
    let fine = steady_couette(12);
    assert!(
        coarse < 0.06 && fine < 0.02,
        "thermal errors {coarse}, {fine}"
    );
    assert!(
        fine < 0.35 * coarse,
        "thermal refinement errors {coarse}, {fine}"
    );
}

#[test]
fn wall_driven_heat_reaches_the_original_temperature_and_shear_viscosity_source() {
    let height = 4;
    let source = material(0.001);
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(MU, 0.05, 3.0, 0.5, 2.0, [0.0, 100.0]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 40_000.0,
            },
            300.0,
            [295.0, 310.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    );
    let mut heated = slab(height, source);
    let mut control = slab(height, source);
    drive(&mut heated);
    drive(&mut control);
    let heat = vec![0.0; height + 2];
    let forces = vec![[0.0; 3]; height + 2];
    let gate = CancelGate::new_clock_free();
    let mut last_hot = None;
    let mut last_cold = None;
    for _ in 0..160 {
        last_hot = Some(
            heated
                .step_with_material_rheology_and_viscous_heating(
                    &heat,
                    &forces,
                    limits(height),
                    solve(),
                    &gate,
                    |point, _| {
                        law.evaluate(point.temperature_k, point.shear_rate_per_s)?
                            .bind(source)
                    },
                )
                .unwrap(),
        );
        last_cold = Some(
            control
                .step_with_material_rheology(
                    &heat,
                    &forces,
                    limits(height),
                    solve(),
                    &gate,
                    |point, _| {
                        law.evaluate(point.temperature_k, point.shear_rate_per_s)?
                            .bind(source)
                    },
                )
                .unwrap(),
        );
    }
    let hot = last_hot.unwrap();
    let cold = last_cold.unwrap();
    assert_ne!(bits(&heated.lattice.grid.f), bits(&control.lattice.grid.f));
    assert_ne!(bits(&heated.lattice.gpop), bits(&control.lattice.gpop));
    let hot_center = hot.transport.applied.iter().find(|s| s.cell == 2).unwrap();
    let cold_center = cold.applied.iter().find(|s| s.cell == 2).unwrap();
    assert!(hot_center.temperature_k > cold_center.temperature_k + 0.01);
    assert!(
        hot_center.transport.dynamic_viscosity_pa_s()
            < cold_center.transport.dynamic_viscosity_pa_s()
    );
    for sample in hot.transport.applied.iter().chain(&hot.transport.endpoint) {
        let retained = sample.material.source();
        assert_eq!(retained.rheology().law(), &law);
        assert_eq!(retained.original_transport(), &source);
        assert_eq!(
            retained.rheology().temperature().temperature_k(),
            sample.temperature_k
        );
        assert_eq!(
            retained.rheology().shear_rate_per_s(),
            sample.material.point().shear_rate_per_s
        );
        assert_eq!(
            retained.rheology().dynamic_viscosity_pa_s(),
            sample.transport.dynamic_viscosity_pa_s()
        );
    }
    assert!(hot.viscous.viscous_heat_j > 0.0);
    assert_eq!(hot.viscous.authored_heat_j, 0.0);
    assert!(heated.last_plate_exchange_si().unwrap().unwrap()[1].wall_work_j < 0.0);
}

#[test]
fn stationary_parity_and_native_or_cancelled_trials_preserve_the_driven_history() {
    let height = 4;
    let source = material(1.0);
    let mut actual = slab(height, source);
    let mut control = slab(height, source);
    actual.set_plate_motion_si([[0.0; 3]; 2], [0.0; 3]).unwrap();
    let heat = vec![0.0; height + 2];
    let forces = vec![[0.0; 3]; height + 2];
    let gate = CancelGate::new_clock_free();
    for _ in 0..3 {
        actual
            .step_with_material_rheology_and_viscous_heating(
                &heat,
                &forces,
                limits(height),
                solve(),
                &gate,
                |_, _| Ok(source),
            )
            .unwrap();
        control
            .step_with_material_rheology_and_viscous_heating(
                &heat,
                &forces,
                limits(height),
                solve(),
                &gate,
                |_, _| Ok(source),
            )
            .unwrap();
        assert_eq!(bits(&actual.lattice.grid.f), bits(&control.lattice.grid.f));
        assert_eq!(bits(&actual.lattice.gpop), bits(&control.lattice.gpop));
    }
    drive(&mut actual);
    drive(&mut control);
    for slab in [&mut actual, &mut control] {
        slab.step_with_material_rheology_and_viscous_heating(
            &heat,
            &forces,
            limits(height),
            solve(),
            &gate,
            |_, _| Ok(source),
        )
        .unwrap();
    }
    let f_before = bits(&actual.lattice.grid.f);
    let g_before = bits(&actual.lattice.gpop);
    let tau_before = actual.lattice.grid.tau.clone();
    let force_before = actual.lattice.grid.fext.clone();
    let velocity_before = actual.lattice.grid.plate_velocities();
    let wall_before = actual.last_plate_exchange_si().unwrap();
    let thermal_before = actual.last_wall_exchange().unwrap();
    let unchanged = |slab: &MaterialSlab3| {
        assert_eq!(bits(&slab.lattice.grid.f), f_before);
        assert_eq!(bits(&slab.lattice.gpop), g_before);
        assert_eq!(slab.lattice.grid.tau, tau_before);
        assert_eq!(slab.lattice.grid.fext, force_before);
        assert_eq!(slab.lattice.grid.plate_velocities(), velocity_before);
        assert_eq!(slab.last_plate_exchange_si().unwrap(), wall_before);
        assert_eq!(slab.last_wall_exchange().unwrap(), thermal_before);
        assert_eq!(slab.lattice.grid.trt_magic(), Some(3.0 / 16.0));
    };
    for bad in [[[0.0; 3], [0.0, 0.0, 0.01]], [[0.0; 3], [1.0, 0.0, 0.0]]] {
        assert!(actual.set_plate_motion_si(bad, [0.0; 3]).is_err());
        unchanged(&actual);
    }
    // This finite actuator force is too large for the resulting wall-adjacent
    // flow. The checked thermal owner must return the native error, not panic
    // after advancing temperature or replacing the accepted wall receipt.
    let large_force = vec![[20_000.0, 0.0, 0.0]; height + 2];
    let failed = actual.step_with_material_rheology_and_viscous_heating(
        &heat,
        &large_force,
        limits(height),
        solve(),
        &gate,
        |_, _| Ok(source),
    );
    assert!(matches!(
        failed,
        Err(MaterialTransportError::Thermal(ThermalSourceError::Force(
            ThermalForcingError::Flow3(_)
        )))
    ));
    unchanged(&actual);
    let cancel = CancelGate::new_clock_free();
    let queries = Cell::new(0);
    let failed = actual.step_with_material_rheology_and_viscous_heating(
        &heat,
        &forces,
        limits(height),
        solve(),
        &cancel,
        |_, _| {
            queries.set(queries.get() + 1);
            // A constant source needs one query per cell; this is the first
            // endpoint query, after native flow and thermal transport.
            if queries.get() == height + 1 {
                cancel.request();
            }
            Ok(source)
        },
    );
    assert!(matches!(failed, Err(MaterialTransportError::Cancelled)));
    assert_eq!(queries.get(), height + 1);
    unchanged(&actual);
    for slab in [&mut actual, &mut control] {
        slab.step_with_material_rheology_and_viscous_heating(
            &heat,
            &forces,
            limits(height),
            solve(),
            &gate,
            |_, _| Ok(source),
        )
        .unwrap();
    }
    assert_eq!(bits(&actual.lattice.grid.f), bits(&control.lattice.grid.f));
    assert_eq!(bits(&actual.lattice.gpop), bits(&control.lattice.gpop));
    assert_eq!(
        actual.last_plate_exchange_si().unwrap(),
        control.last_plate_exchange_si().unwrap()
    );
    assert_eq!(
        actual.last_wall_exchange().unwrap(),
        control.last_wall_exchange().unwrap()
    );
}
