use super::*;
use crate::multicomponent::characteristic_tests::caloric_gas;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-10 * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn axial(
    gas: &MixtureGas,
    densities: &[f64],
    velocity: f64,
    normal: f64,
    pressure: f64,
) -> ([f64; 4], FaceFlux) {
    let state = gas.conserved(densities, velocity, 1.0).unwrap();
    let mut flux = [0.0; 4];
    let report = gas
        .pressure_outlet_flux(state.as_cell(), normal, pressure, &mut flux, || false)
        .unwrap();
    (flux, report)
}

#[test]
fn g1_pressure_outlet_matches_mixture_shock_fan_and_sonic_fluxes() {
    let gas = caloric_gas(false, 0.0);
    for densities in [[1.0_f64, 0.0], [0.3, 0.7]] {
        // Manufactured species have R_A=1, Cp_A=3.5 and R_B=.5, Cp_B=2.5.
        // These independent constant-Cp Riemann relations include the actual
        // source mixture rather than either component's gamma.
        let r = densities[0] + 0.5 * densities[1];
        let cp = 3.5 * densities[0] + 2.5 * densities[1];
        let cv = cp - r;
        let gamma = cp / cv;
        let ci = gamma.sqrt();
        for (ui, backpressure) in [
            (0.2_f64, 0.8_f64),
            (0.6, 1.2),
            (0.0, 0.02),
            (0.0, 1e-12),
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 0.8),
            (2.0, 1.1),
        ] {
            let (density, velocity, pressure) = if backpressure > 1.0 {
                let shock_speed =
                    ui - ci * (((gamma + 1.0) * backpressure + gamma - 1.0) / (2.0 * gamma)).sqrt();
                if shock_speed >= 0.0 {
                    (1.0, ui, 1.0)
                } else {
                    let ratio = (gamma - 1.0) / (gamma + 1.0);
                    let rho = (backpressure + ratio) / (1.0 + ratio * backpressure);
                    let u = ui
                        - (backpressure - 1.0)
                            * (2.0 / ((gamma + 1.0) * (backpressure + ratio))).sqrt();
                    (rho, u, backpressure)
                }
            } else if ui >= ci {
                (1.0, ui, 1.0)
            } else {
                let acoustic = ci * backpressure.powf((gamma - 1.0) / (2.0 * gamma));
                let u = ui + 2.0 * (ci - acoustic) / (gamma - 1.0);
                if u <= acoustic {
                    (backpressure.powf(1.0 / gamma), u, backpressure)
                } else {
                    let sonic = (2.0 * ci + (gamma - 1.0) * ui) / (gamma + 1.0);
                    let ratio = sonic / ci;
                    (
                        ratio.powf(2.0 / (gamma - 1.0)),
                        sonic,
                        ratio.powf(2.0 * gamma / (gamma - 1.0)),
                    )
                }
            };
            let temperature = pressure / (density * r);
            let mass = density * velocity;
            let enthalpy =
                cv * (temperature - 0.001) + 0.5 * velocity * velocity + pressure / density;
            for normal in [-1.0, 1.0] {
                let (flux, report) = axial(&gas, &densities, normal * ui, normal, backpressure);
                for k in 0..2 {
                    close(flux[k], normal * densities[k] * mass);
                }
                if densities[1] == 0.0 {
                    assert_eq!(flux[1], 0.0);
                }
                close(flux[2], mass * velocity + pressure);
                close(flux[3], normal * mass * enthalpy);
                assert!(!report.used_hll_fallback);
                assert!(
                    report.signal_speed >= velocity + (gamma * pressure / density).sqrt() - 1e-12
                );
            }
        }
    }
}

#[test]
fn g1_variable_calorics_keep_the_original_entropy_and_energy_references() {
    let densities = [0.3, 0.7];
    let r = densities[0] + 0.5 * densities[1];
    let a = 3.5 * densities[0] + 2.5 * densities[1];
    let b = 0.25 * densities[0] + 0.35 * densities[1];
    let pressure = 0.8;
    let initial_temperature = 1.0 / r;
    let ordinary = caloric_gas(true, 0.0);
    let (expected, _) = axial(&ordinary, &densities, 0.2, 1.0, pressure);
    let mass = expected[0] + expected[1];
    let velocity = (expected[2] - pressure) / mass;
    let density = mass / velocity;
    let temperature = pressure / (density * r);
    close(
        a * (temperature / initial_temperature).ln() + b * (temperature - initial_temperature)
            - r * pressure.ln(),
        0.0,
    );
    let thermal =
        (a - r) * (temperature - 0.001) + 0.5 * b * (temperature * temperature - 0.001_f64.powi(2));
    close(
        expected[3] / mass,
        thermal + r * temperature + 0.5 * velocity * velocity,
    );
    for reference in [-1e18, 1e18] {
        let shifted = caloric_gas(true, reference);
        let (actual, _) = axial(&shifted, &densities, 0.2, 1.0, pressure);
        assert_eq!(actual, expected);
    }
}

#[test]
fn g1_near_vacuum_mixture_outlets_admit_the_sonic_face_without_a_cold_conservative_tail() {
    use crate::multicomponent::tests::{gas as mixture, source};
    use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

    let gas = mixture(vec![
        source("A", R, 3.5, 0.0, 30.0),
        source("B", 2.0 * R, 3.5, 0.4, -60.0),
        source("absent-C", 0.5 * R, 4.5, 0.2, 7.0),
    ]);
    let densities = [0.35, 0.65, 0.0];
    // Independent variable-Cp sonic oracle for this exact three-species fixture.
    // Cp(T)=a+b*T, Rmix=r. Integrating in sqrt(T) removes the cold singularity.
    let r = 0.675_f64;
    let a = 3.5 * r;
    let b = 0.13;
    let initial_temperature = 1.0 / r;
    let acoustic = |t: f64| (r * t * (a + b * t) / (a - r + b * t)).sqrt();
    let invariant = |t: f64| {
        let panels = 512;
        let start = t.sqrt();
        let width = (initial_temperature.sqrt() - start) / panels as f64;
        let integrand = |z: f64| 2.0 * ((a + b * z * z) * (a - r + b * z * z) / r).sqrt();
        let mut value = integrand(start) + integrand(initial_temperature.sqrt());
        for panel in 1..panels {
            value +=
                if panel % 2 == 0 { 2.0 } else { 4.0 } * integrand(start + panel as f64 * width);
        }
        value * width / 3.0
    };
    let (mut low, mut high) = (0.001_f64, initial_temperature);
    for _ in 0..64 {
        let t = low.midpoint(high);
        if invariant(t) > acoustic(t) {
            low = t;
        } else {
            high = t;
        }
    }
    let temperature = low.midpoint(high);
    let speed = acoustic(temperature);
    let pressure = ((a * (temperature / initial_temperature).ln()
        + b * (temperature - initial_temperature))
        / r)
        .exp();
    let density = pressure / (r * temperature);
    let mass = density * speed;
    let thermal =
        (a - r) * (temperature - 0.001) + 0.5 * b * (temperature * temperature - 0.001_f64.powi(2));
    let enthalpy = thermal + r * temperature + 0.5 * speed * speed;
    let inside = gas.conserved(&densities, 0.0, 1.0).unwrap();
    let bath = gas.conserved(&[0.0, 2e-12, 0.0], 0.0, 1e-12).unwrap();
    for normal in [-1.0, 1.0] {
        for backpressure in [0.01, 1e-12, 0.0] {
            let mut flux = [f64::NAN; 5];
            let receipt = gas
                .pressure_outlet_flux(inside.as_cell(), normal, backpressure, &mut flux, || false)
                .unwrap();
            close(flux[0], normal * 0.35 * mass);
            close(flux[1], normal * 0.65 * mass);
            assert_eq!(flux[2], 0.0);
            close(flux[3], mass * speed + pressure);
            close(flux[4], normal * mass * enthalpy);
            assert!(receipt.signal_speed >= 2.0 * speed - 1e-9);
        }
        let mut reservoir = [f64::NAN; 5];
        gas.reservoir_flux(
            inside.as_cell(),
            bath.as_cell(),
            normal,
            &mut reservoir,
            || false,
        )
        .unwrap();
        close(reservoir[0], normal * 0.35 * mass);
        close(reservoir[1], normal * 0.65 * mass);
        assert_eq!(reservoir[2], 0.0);
        close(reservoir[3], mass * speed + pressure);
        close(reservoir[4], normal * mass * enthalpy);

        let shear = [0.2, -0.3];
        let mut spatial = [f64::NAN; 7];
        gas.pressure_outlet_spatial_flux(
            SpatialCell {
                partial_densities: &densities,
                normal_momentum: 0.0,
                transverse_momentum: shear,
                energy: inside.energy + 0.5 * (shear[0] * shear[0] + shear[1] * shear[1]),
            },
            normal,
            1e-12,
            &mut spatial,
            || false,
        )
        .unwrap();
        for k in 0..4 {
            close(spatial[k], reservoir[k]);
        }
        close(spatial[4], normal * mass * shear[0]);
        close(spatial[5], normal * mass * shear[1]);
        close(
            spatial[6],
            normal * mass * (enthalpy + 0.5 * (shear[0] * shear[0] + shear[1] * shear[1])),
        );
    }

    // A source window that excludes the sonic state still refuses. Only an
    // unused remote tail may be omitted; the actual face is never extrapolated.
    let mut domain = gas.domain();
    domain.minimum_temperature = 1.4;
    let restricted = MixtureGas::new(gas.source().clone(), domain, || false).unwrap();
    let inside = restricted.conserved(&densities, 0.0, 1.0).unwrap();
    let mut refused = [37.0; 5];
    for backpressure in [1e-12, 0.0] {
        assert!(matches!(
            restricted.pressure_outlet_flux(
                inside.as_cell(),
                1.0,
                backpressure,
                &mut refused,
                || false,
            ),
            Err(Error::CharacteristicPath(_))
        ));
        assert_eq!(refused, [37.0; 5]);
    }
}

#[test]
fn g3_planar_and_spatial_outlets_transport_all_shear_work_once() {
    let gas = caloric_gas(true, 0.0);
    let densities = [0.3, 0.7];
    for normal in [-1.0, 1.0] {
        let (reference, report) = axial(&gas, &densities, normal * 0.2, normal, 0.8);
        let mass = reference[0] + reference[1];
        for transverse in [[0.0, 0.0], [0.3, -0.4], [-2.0, 3.0]] {
            let planar = NormalPrimitive {
                partial_densities: &densities,
                normal_velocity: normal * 0.2,
                transverse_velocity: transverse[0],
                pressure: 1.0,
            };
            let spatial = SpatialPrimitive {
                partial_densities: &densities,
                normal_velocity: planar.normal_velocity,
                transverse_velocity: transverse,
                pressure: 1.0,
            };
            let mut q2 = [0.0; 5];
            let mut q3 = [0.0; 6];
            gas.normal_conserved_into(planar, &mut q2).unwrap();
            gas.spatial_conserved_into(spatial, &mut q3).unwrap();
            let mut f2 = [0.0; 5];
            let mut f3 = [0.0; 6];
            let r2 = gas
                .pressure_outlet_normal_flux(
                    NormalCell {
                        partial_densities: &q2[..2],
                        normal_momentum: q2[2],
                        transverse_momentum: q2[3],
                        energy: q2[4],
                    },
                    normal,
                    0.8,
                    &mut f2,
                    || false,
                )
                .unwrap();
            let r3 = gas
                .pressure_outlet_spatial_flux(
                    SpatialCell {
                        partial_densities: &q3[..2],
                        normal_momentum: q3[2],
                        transverse_momentum: [q3[3], q3[4]],
                        energy: q3[5],
                    },
                    normal,
                    0.8,
                    &mut f3,
                    || false,
                )
                .unwrap();
            let mut traced2 = [0.0; 5];
            let mut traced3 = [0.0; 6];
            gas.pressure_outlet_normal_flux_from_primitive(
                planar,
                normal,
                0.8,
                &mut traced2,
                || false,
            )
            .unwrap();
            gas.pressure_outlet_spatial_flux_from_primitive(
                spatial,
                normal,
                0.8,
                &mut traced3,
                || false,
            )
            .unwrap();
            for k in 0..3 {
                close(f2[k], reference[k]);
                close(f3[k], reference[k]);
            }
            close(f2[3], mass * transverse[0]);
            close(f3[3], mass * transverse[0]);
            close(f3[4], mass * transverse[1]);
            close(f2[4], reference[3] + 0.5 * mass * transverse[0].powi(2));
            close(
                f3[5],
                reference[3] + 0.5 * mass * (transverse[0].powi(2) + transverse[1].powi(2)),
            );
            for (actual, expected) in f2.into_iter().zip(traced2) {
                close(actual, expected);
            }
            for (actual, expected) in f3.into_iter().zip(traced3) {
                close(actual, expected);
            }
            close(r2.signal_speed, report.signal_speed);
            close(r3.signal_speed, report.signal_speed);
            assert!(!r2.used_hll_fallback && !r3.used_hll_fallback);
        }
    }
}

#[test]
fn g0_backflow_invalid_axes_and_unresolved_shear_leave_outputs_untouched() {
    let gas = caloric_gas(true, 0.0);
    let state = gas.conserved(&[0.3, 0.7], 0.0, 1.0).unwrap();
    for (normal, pressure) in [
        (1.0, 2.0),
        (-1.0, 2.0),
        (0.0, 0.8),
        (1.0, -1.0),
        (1.0, f64::NAN),
        (1.0, f64::INFINITY),
    ] {
        let mut output = [7.0; 4];
        assert!(
            gas.pressure_outlet_flux(state.as_cell(), normal, pressure, &mut output, || false)
                .is_err()
        );
        assert_eq!(output, [7.0; 4]);
    }
    let mut short = [7.0; 3];
    assert!(
        gas.pressure_outlet_flux(state.as_cell(), 1.0, 0.8, &mut short, || false)
            .is_err()
    );
    assert_eq!(short, [7.0; 3]);
    let planar = NormalPrimitive {
        partial_densities: &state.partial_densities,
        normal_velocity: 0.0,
        transverse_velocity: 1e10,
        pressure: 1.0,
    };
    let spatial = SpatialPrimitive {
        partial_densities: &state.partial_densities,
        normal_velocity: 0.0,
        transverse_velocity: [0.0, 1e10],
        pressure: 1.0,
    };
    let mut p = [7.0; 5];
    let mut s = [7.0; 6];
    assert!(
        gas.pressure_outlet_normal_flux_from_primitive(planar, 1.0, 0.8, &mut p, || false)
            .is_err()
    );
    assert!(
        gas.pressure_outlet_spatial_flux_from_primitive(spatial, 1.0, 0.8, &mut s, || false)
            .is_err()
    );
    assert_eq!(p, [7.0; 5]);
    assert_eq!(s, [7.0; 6]);
}

#[test]
fn g4_all_outlet_paths_publish_only_after_the_final_cancellation_check() {
    let gas = caloric_gas(true, 0.0);
    let densities = [0.3, 0.7];
    let axial = gas.conserved(&densities, 0.2, 1.0).unwrap();
    let planar = NormalPrimitive {
        partial_densities: &densities,
        normal_velocity: 0.2,
        transverse_velocity: 0.3,
        pressure: 1.0,
    };
    let spatial = SpatialPrimitive {
        partial_densities: &densities,
        normal_velocity: 0.2,
        transverse_velocity: [0.3, -0.4],
        pressure: 1.0,
    };
    let mut q2 = [0.0; 5];
    let mut q3 = [0.0; 6];
    gas.normal_conserved_into(planar, &mut q2).unwrap();
    gas.spatial_conserved_into(spatial, &mut q3).unwrap();
    let planar_cell = NormalCell {
        partial_densities: &q2[..2],
        normal_momentum: q2[2],
        transverse_momentum: q2[3],
        energy: q2[4],
    };
    let spatial_cell = SpatialCell {
        partial_densities: &q3[..2],
        normal_momentum: q3[2],
        transverse_momentum: [q3[3], q3[4]],
        energy: q3[5],
    };
    for path in 0..5 {
        let invoke = |output: &mut [f64], cancelled: &mut dyn FnMut() -> bool| match path {
            0 => gas.pressure_outlet_flux(axial.as_cell(), 1.0, 0.8, output, cancelled),
            1 => gas.pressure_outlet_normal_flux(planar_cell, 1.0, 0.8, output, cancelled),
            2 => {
                gas.pressure_outlet_normal_flux_from_primitive(planar, 1.0, 0.8, output, cancelled)
            }
            3 => gas.pressure_outlet_spatial_flux(spatial_cell, 1.0, 0.8, output, cancelled),
            _ => gas
                .pressure_outlet_spatial_flux_from_primitive(spatial, 1.0, 0.8, output, cancelled),
        };
        let size = [4, 5, 5, 6, 6][path];
        let mut output = [7.0; 6];
        let mut count = 0;
        invoke(&mut output[..size], &mut || {
            count += 1;
            false
        })
        .unwrap();
        assert!(count > 4);
        for stop in [1, count / 2, count] {
            output.fill(7.0);
            let mut visits = 0;
            assert_eq!(
                invoke(&mut output[..size], &mut || {
                    visits += 1;
                    visits >= stop
                }),
                Err(Error::Cancelled),
                "path={path}, stop={stop}"
            );
            assert_eq!(output, [7.0; 6]);
        }
    }
}
