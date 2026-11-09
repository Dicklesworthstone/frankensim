use super::*;
use fs_material::fluid_acoustics::{FluidAcousticState, thermoviscous::ThermoviscousFluidState};
use fs_material::fluid_transport::{FluidTransportState, ThermalFluidTransportState};
use fs_math::det;

fn medium(expansion: f64, viscosity: f64, bulk: f64) -> BulkThermoviscousEquilibrium {
    let transport = FluidTransportState::try_new(2.0, 8.0, 0.6, viscosity).unwrap();
    let beta_s = 1.0 / 200.0;
    let beta_t = beta_s + expansion * expansion * 12.5 / 16.0;
    BulkThermoviscousEquilibrium::from_states(
        ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(2.0, 10.0).unwrap(),
            transport,
            8.0 * beta_s / beta_t,
        )
        .unwrap(),
        ThermalFluidTransportState::try_new(transport, 12.5, Some(expansion)).unwrap(),
        bulk,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-300),
        "{a:e} != {b:e}"
    );
}
fn near_complex(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-300),
        "{a:?} != {b:?}"
    );
}

#[test]
fn g1_loaded_duct_matching_delegates_physical_fields_and_actual_area_once() {
    let material = medium(-0.04, 0.003, 0.002);
    let (length, omega, velocity) = (0.04, 100.0, C64::new(0.2, -0.3));
    let native = material.prepare_wave_absorption(omega).unwrap();
    let mut heats = Vec::new();
    for radius in [0.01, 0.02] {
        let cylinder =
            BulkTractionCylinder::new(material.absorption(), radius, length, omega).unwrap();
        let area = core::f64::consts::PI * radius * radius;
        let load = cylinder
            .wave()
            .traction_impedance_pa_s_m()
            .scale(1.0 / area);
        let drive = velocity.scale(area);
        let response = cylinder
            .absorption_response(&material, Some(load), drive, ComposeOptions::default())
            .unwrap();
        assert_eq!(response.equilibrium(), &material);
        assert_eq!(response.area_m2(), area);
        assert_eq!(response.length_m(), length);
        assert_eq!(response.angular_frequency_rad_s(), omega);
        assert_eq!(response.inlet_volume_flow_m3_s(), drive);
        let original = cylinder
            .response(Some(load), ComposeOptions::default())
            .unwrap();
        assert_eq!(
            response.response().inlet_traction_impedance_pa_s_m3,
            original.inlet_traction_impedance_pa_s_m3
        );
        for x in [0.0, 0.013, length] {
            let a = response.at_distance(x).unwrap();
            let b = native.at_distance(velocity, x).unwrap();
            for (actual, expected) in [
                (a.velocity_m_s, b.velocity_m_s),
                (a.pressure_pa, b.pressure_pa),
                (a.normal_traction_pa, b.normal_traction_pa),
                (a.temperature_perturbation_k, b.temperature_perturbation_k),
                (a.conductive_heat_flux_w_m2, b.conductive_heat_flux_w_m2),
            ] {
                near_complex(actual, expected, 3e-14);
            }
        }
        let a = response.absorption_between(0.0, length).unwrap();
        let b = native
            .absorption_between(velocity, 0.0, length, area)
            .unwrap();
        near(a.viscous_heat_w, b.viscous_heat_w, 3e-14);
        near(a.conductive_entropy_w_k, b.conductive_entropy_w_k, 3e-14);
        near(a.acoustic_power_loss_w, b.acoustic_power_loss_w, 3e-14);
        heats.push(a.acoustic_power_loss_w);
    }
    near(heats[1], 4.0 * heats[0], 3e-14);
}

#[test]
fn g1_reflective_loads_recover_original_endpoints_and_viscous_work_balance() {
    let material = medium(0.0, 0.003, 0.002);
    let (radius, length, omega, drive) = (0.02, 0.04, 1000.0, C64::new(2e-5, -1e-5));
    let cylinder = BulkTractionCylinder::new(material.absorption(), radius, length, omega).unwrap();
    let area = core::f64::consts::PI * radius * radius;
    let scale = drive.abs() / area;
    let reference = cylinder.traction_reference_pa_s_m3();
    for load in [
        None,
        Some(C64::ZERO),
        Some(C64::new(1.4, 0.3).scale(reference)),
    ] {
        let result = cylinder
            .absorption_response(&material, load, drive, ComposeOptions::default())
            .unwrap();
        let response = result.response();
        let front = result.at_distance(0.0).unwrap();
        let back = result.at_distance(length).unwrap();
        near_complex(front.velocity_m_s.scale(area), drive, 2e-13);
        near_complex(
            front.normal_traction_pa,
            response.inlet_traction_impedance_pa_s_m3 * drive,
            2e-13,
        );
        near_complex(
            front.pressure_pa,
            response.inlet_pressure_over_flow_pa_s_m3 * drive,
            2e-13,
        );
        assert!(
            (back.velocity_m_s - response.outlet_over_inlet_flow * drive / C64::from_re(area))
                .abs()
                < 3e-13 * scale
        );
        assert!(
            (back.normal_traction_pa - response.outlet_traction_over_inlet_flow_pa_s_m3 * drive)
                .abs()
                < 3e-13 * reference * drive.abs()
        );
        let heat = result.absorption_between(0.0, length).unwrap();
        let outlet_work = (response.outlet_traction_over_inlet_flow_pa_s_m3
            * response.outlet_over_inlet_flow.conj())
        .re;
        let original_work =
            0.5 * (response.inlet_traction_impedance_pa_s_m3.re - outlet_work) * drive.norm_sq();
        assert!(heat.viscous_heat_w > 0.0);
        assert_eq!(heat.conductive_entropy_w_k, 0.0);
        near(heat.viscous_heat_w, original_work, 3e-12);
        assert!(heat.balance_residual_w.abs() < 3e-13 * heat.left_acoustic_power_w.abs());
    }
}

#[test]
fn g1_native_duct_thermal_absorption_retains_more_than_mechanical_work() {
    let material = medium(0.08, 0.003, 0.002);
    let (radius, length, omega, drive) = (0.02, 0.01, 10_000.0, C64::new(2e-5, -1e-5));
    let cylinder = BulkTractionCylinder::new(material.absorption(), radius, length, omega).unwrap();
    let load = cylinder.wave().traction_impedance_pa_s_m() / C64::from_re(cylinder.area_m2);
    let result = cylinder
        .absorption_response(&material, Some(load), drive, ComposeOptions::default())
        .unwrap();
    let heat = result.absorption_between(0.0, length).unwrap();
    let front = result.at_distance(0.0).unwrap();
    let back = result.at_distance(length).unwrap();
    let mechanical_work =
        (front.mechanical_power_flux_w_m2 - back.mechanical_power_flux_w_m2) * result.area_m2();
    assert!(heat.acoustic_power_loss_w > 10.0 * mechanical_work);
    assert!(heat.thermal_availability_loss_w > 100.0 * heat.viscous_heat_w);
    near(
        heat.acoustic_power_loss_w,
        heat.viscous_heat_w + 12.5 * heat.conductive_entropy_w_k,
        3e-15,
    );
    assert!(heat.balance_residual_w.abs() < 3e-14 * heat.left_acoustic_power_w);
    let doubled = cylinder
        .absorption_response(
            &material,
            Some(load),
            drive.scale(2.0),
            ComposeOptions::default(),
        )
        .unwrap();
    near(
        doubled
            .absorption_between(0.0, length)
            .unwrap()
            .acoustic_power_loss_w,
        4.0 * heat.acoustic_power_loss_w,
        3e-14,
    );
}

#[test]
fn g1_tiny_reflection_comes_from_original_load_not_cancelled_endpoint_fields() {
    let material = medium(0.0, 0.003, 0.002);
    let (radius, length, omega) = (0.02, 0.04, 1e-8);
    let cylinder = BulkTractionCylinder::new(material.absorption(), radius, length, omega).unwrap();
    let z = cylinder
        .wave()
        .traction_impedance_pa_s_m()
        .scale(1.0 / cylinder.area_m2);
    let drive = C64::new(2e-5, -1e-5);
    let matched = cylinder
        .response(Some(z), ComposeOptions::default())
        .unwrap();
    assert_eq!(
        cylinder.backward_at_outlet(&matched, Some(z), drive),
        C64::ZERO
    );
    let load = C64::new(z.re, z.im + 1e-24);
    assert_ne!(load, z);
    let response = cylinder
        .response(Some(load), ComposeOptions::default())
        .unwrap();
    let reflected = cylinder.backward_at_outlet(&response, Some(load), drive);
    let reflection = (z - load) / (z + load);
    let k = cylinder.wave().wavenumber_rad_m();
    let (s, c) = det::sin_cos(k.re * length);
    let propagation = C64::new(c, s).scale(det::exp(-k.im * length));
    let forward = (drive / C64::from_re(cylinder.area_m2))
        / (C64::ONE + reflection * propagation * propagation);
    let expected = reflection * forward * propagation;
    assert!(reflected.abs() > 0.0 && reflected.abs() < 1e-25 * forward.abs());
    near_complex(reflected, expected, 3e-13);
}

#[test]
fn g1_opaque_loaded_bulk_cylinder_and_silence_remain_finite() {
    let material = medium(0.08, 0.003, 0.002);
    let (radius, length, omega) = (0.02, 10_000.0, 1000.0);
    let cylinder = BulkTractionCylinder::new(material.absorption(), radius, length, omega).unwrap();
    for drive in [C64::new(2e-5, -1e-5), C64::ZERO] {
        let result = cylinder
            .absorption_response(&material, None, drive, ComposeOptions::default())
            .unwrap();
        assert_eq!(result.at_distance(length).unwrap().velocity_m_s, C64::ZERO);
        let heat = result.absorption_between(0.0, length).unwrap();
        assert_eq!(heat.right_acoustic_power_w, 0.0);
        if drive == C64::ZERO {
            assert_eq!(heat.acoustic_power_loss_w, 0.0);
            assert_eq!(heat.balance_residual_w, 0.0);
        } else {
            near(
                heat.acoustic_power_loss_w,
                heat.left_acoustic_power_w,
                3e-14,
            );
        }
    }
}

#[test]
fn g0_bound_bulk_absorption_refuses_other_material_invalid_drive_and_geometry() {
    let material = medium(0.04, 0.003, 0.002);
    let cylinder = BulkTractionCylinder::new(material.absorption(), 0.02, 0.04, 1000.0).unwrap();
    for wrong in [
        medium(0.08, 0.003, 0.002),
        medium(0.04, 0.006, 0.002),
        medium(0.04, 0.003, 0.004),
    ] {
        assert!(
            cylinder
                .absorption_response(&wrong, None, C64::ONE, ComposeOptions::default())
                .is_err()
        );
    }
    for bad in [C64::new(f64::NAN, 0.0), C64::new(0.0, f64::INFINITY)] {
        assert!(
            cylinder
                .absorption_response(&material, None, bad, ComposeOptions::default())
                .is_err()
        );
        assert!(
            cylinder
                .absorption_response(&material, Some(bad), C64::ONE, ComposeOptions::default())
                .is_err()
        );
    }
    let response = cylinder
        .absorption_response(
            &material,
            None,
            C64::from_re(1e-5),
            ComposeOptions::default(),
        )
        .unwrap();
    for bad in [-1.0, 0.05, f64::INFINITY] {
        assert!(response.at_distance(bad).is_err());
    }
    assert!(response.absorption_between(0.03, 0.02).is_err());
    assert!(response.absorption_between(0.0, 0.05).is_err());
    assert_eq!(
        response
            .absorption_between(0.02, 0.02)
            .unwrap()
            .acoustic_power_loss_w,
        0.0
    );
}
