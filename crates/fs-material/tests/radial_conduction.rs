//! Radial geometry and finite source-backed storage use one shared heat owner.
use fs_material::thermal::conduction::{ConductivityCurve, LayerConductionConfig, LayerGeometry};
use fs_thermal::enthalpy::DepositionError;

fn config() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12,
    }
}
fn constant(k: f64) -> ConductivityCurve {
    ConductivityCurve::new(vec![[200.0, k], [1000.0, k]]).unwrap()
}
fn close(a: f64, b: f64, allowance: f64) {
    assert!(
        (a - b).abs() <= allowance,
        "{a:.16e} != {b:.16e}, allowance={allowance}"
    );
}

#[test]
fn g1_radial_metrics_recover_independent_shell_resistance_and_interface_temperatures() {
    let a = constant(20.0);
    let b = constant(0.6);
    let pi = core::f64::consts::PI;
    for cylinder in [false, true] {
        let geometry = if cylinder {
            LayerGeometry::cylindrical([0.01, 0.04], 0.02, 0.7).unwrap()
        } else {
            LayerGeometry::spherical([0.01, 0.04], 0.02).unwrap()
        };
        let [ra, rb] = if cylinder {
            [
                2.0_f64.ln() / (2.0 * pi * 0.7 * 20.0),
                2.0_f64.ln() / (2.0 * pi * 0.7 * 0.6),
            ]
        } else {
            [
                (100.0 - 50.0) / (4.0 * pi * 20.0),
                (50.0 - 25.0) / (4.0 * pi * 0.6),
            ]
        };
        let contact = 0.03;
        let rc = contact / geometry.contact_area_m2();
        for temperatures in [[600.0, 300.0], [300.0, 600.0], [450.0; 2]] {
            let result = geometry
                .conduction([&a, &b], contact)
                .unwrap()
                .evaluate(temperatures, config(), || Ok(()))
                .unwrap();
            let g = 1.0 / (ra + rc + rb);
            let q = g * (temperatures[0] - temperatures[1]);
            close(result.power_a_to_b_w, q, 1e-9);
            close(result.conductance_w_k, g, 1e-11);
            close(
                result.interface_temperature_k[0],
                temperatures[0] - q * ra,
                1e-9,
            );
            close(
                result.interface_temperature_k[1],
                temperatures[1] + q * rb,
                1e-9,
            );
            close(result.endpoint_derivative_w_k[0], g, 1e-10);
            close(result.endpoint_derivative_w_k[1], -g, 1e-10);
        }
        let reversed = if cylinder {
            LayerGeometry::cylindrical([0.04, 0.01], 0.02, 0.7).unwrap()
        } else {
            LayerGeometry::spherical([0.04, 0.01], 0.02).unwrap()
        };
        let forward = geometry
            .conduction([&a, &b], contact)
            .unwrap()
            .evaluate([600.0, 300.0], config(), || Ok(()))
            .unwrap();
        let reverse = reversed
            .conduction([&b, &a], contact)
            .unwrap()
            .evaluate([300.0, 600.0], config(), || Ok(()))
            .unwrap();
        assert_eq!(reverse.power_a_to_b_w, -forward.power_a_to_b_w);
        assert_eq!(
            reverse.interface_temperature_k,
            [
                forward.interface_temperature_k[1],
                forward.interface_temperature_k[0]
            ]
        );
    }
}

#[test]
fn g1_nonlinear_radial_heat_integrates_original_k_and_retains_exact_endpoint_tangents() {
    let k = ConductivityCurve::new(vec![[200.0, 4.0], [450.0, 6.5], [1000.0, 12.0]]).unwrap();
    let pi = core::f64::consts::PI;
    for cylinder in [false, true] {
        let geometry = if cylinder {
            LayerGeometry::cylindrical([0.01, 0.04], 0.02, 0.7).unwrap()
        } else {
            LayerGeometry::spherical([0.01, 0.04], 0.02).unwrap()
        };
        let metric = if cylinder {
            4.0_f64.ln() / (2.0 * pi * 0.7)
        } else {
            (100.0 - 25.0) / (4.0 * pi)
        };
        let result = geometry
            .conduction([&k; 2], 0.0)
            .unwrap()
            .evaluate([600.0, 300.0], config(), || Ok(()))
            .unwrap();
        // Independent integral of k=2+0.01*T across the complete radial domain.
        let integral = 2.0 * 300.0 + 0.005 * (600.0 * 600.0 - 300.0 * 300.0);
        close(result.power_a_to_b_w, integral / metric, 1e-8);
        close(result.endpoint_derivative_w_k[0], 8.0 / metric, 1e-10);
        close(result.endpoint_derivative_w_k[1], -5.0 / metric, 1e-10);
        let mut calls = 0;
        geometry
            .conduction([&k; 2], 0.0)
            .unwrap()
            .evaluate([600.0, 300.0], config(), || {
                calls += 1;
                Ok(())
            })
            .unwrap();
        let mut retry = 0;
        assert_eq!(
            geometry
                .conduction([&k; 2], 0.0)
                .unwrap()
                .evaluate([600.0, 300.0], config(), || {
                    retry += 1;
                    if retry == calls {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                }),
            Err(DepositionError::Cancelled)
        );
    }
}

#[test]
fn g3_thin_shells_preserve_radial_differences_and_reject_invalid_geometry() {
    let ri = 1.0_f64.next_down();
    let rc = 1.0;
    let ro = 1.0_f64.next_up();
    let spherical = LayerGeometry::spherical([ri, ro], rc).unwrap();
    let cylinder = LayerGeometry::cylindrical([ri, ro], rc, 1.0).unwrap();
    for geometry in [spherical, cylinder] {
        let distances = geometry.equivalent_distances_m();
        assert!(distances.into_iter().all(|d| d > 0.0));
        close(distances[0] / (rc - ri), 1.0, 5e-16);
        close(distances[1] / (ro - rc), 1.0, 5e-16);
        assert_eq!(geometry.radial_coordinates_m(), Some(([ri, ro], rc)));
    }
    for radii in [
        [0.0, 2.0],
        [1.0, 2.0],
        [2.0, 3.0],
        [f64::NAN, 2.0],
        [0.5, f64::INFINITY],
    ] {
        assert!(LayerGeometry::spherical(radii, 1.0).is_err());
        assert!(LayerGeometry::cylindrical(radii, 1.0, 1.0).is_err());
    }
    assert!(LayerGeometry::spherical([1e200, 3e200], 2e200).is_err());
    assert!(LayerGeometry::cylindrical([0.5, 2.0], 1.0, 0.0).is_err());
    let k = constant(1.0);
    let finite = LayerGeometry::spherical([1.0, 2.0], 1.5)
        .unwrap()
        .conduction([&k; 2], 0.0)
        .unwrap()
        .evaluate([600.0, 300.0], config(), || Ok(()))
        .unwrap();
    let far = LayerGeometry::spherical([1.0, 1e6], 2.0)
        .unwrap()
        .conduction([&k; 2], 0.0)
        .unwrap()
        .evaluate([600.0, 300.0], config(), || Ok(()))
        .unwrap();
    close(finite.conductance_w_k, 8.0 * core::f64::consts::PI, 1e-10);
    close(
        far.conductance_w_k / (4.0 * core::f64::consts::PI),
        1.0,
        1.1e-6,
    );
}

#[path = "radial_conduction/finite.rs"]
mod finite;
