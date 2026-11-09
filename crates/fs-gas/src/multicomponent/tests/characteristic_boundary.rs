//! Original-species acoustic characteristics through the actual mixture owner.

use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-9 * (1.0 + expected.abs()),
        "actual={actual}, expected={expected}"
    );
}

fn physical_flux(state: &Conserved, velocity: f64, pressure: f64) -> [f64; 4] {
    [
        state.partial_densities[0] * velocity,
        state.partial_densities[1] * velocity,
        state.momentum * velocity + pressure,
        (state.energy + pressure) * velocity,
    ]
}

fn close_flux(actual: &[f64; 4], expected: [f64; 4]) {
    for (&value, reference) in actual.iter().zip(expected) {
        close(value, reference);
    }
}

fn variable_source(formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.2, formation),
        source("B", 2.0 * R, 4.2, 0.4, -2.0 * formation),
    ])
}

#[test]
fn g3_multidimensional_characteristic_outlets_retain_species_shear_and_full_energy() {
    let gas = variable_source(1e18);
    for densities in [[1.0, 0.0], [0.0, 2.0], [0.2, 0.8]] {
        for normal in [-1.0, 1.0] {
            for (speed, reference_pressure) in [(0.3, 0.9), (0.3, 0.001), (3.0, 1e40)] {
                let velocity = normal * speed;
                let axial = gas.conserved(&densities, velocity, 1.0).unwrap();
                let mut expected = [0.0; 4];
                let axial_report = gas
                    .characteristic_outlet_flux(
                        axial.as_cell(),
                        normal,
                        reference_pressure,
                        normal * 0.3,
                        &mut expected,
                        || false,
                    )
                    .unwrap();
                let mass_flux = expected[0] + expected[1];
                for shear in [[0.0, 0.0], [2.0, -3.0], [-3.0, 2.0]] {
                    let planar = NormalPrimitive {
                        partial_densities: &densities,
                        normal_velocity: velocity,
                        transverse_velocity: shear[0],
                        pressure: 1.0,
                    };
                    let spatial = SpatialPrimitive {
                        partial_densities: &densities,
                        normal_velocity: velocity,
                        transverse_velocity: shear,
                        pressure: 1.0,
                    };
                    let mut q2 = [0.0; 5];
                    let mut q3 = [0.0; 6];
                    gas.normal_conserved_into(planar, &mut q2).unwrap();
                    gas.spatial_conserved_into(spatial, &mut q3).unwrap();
                    let mut f2 = [0.0; 5];
                    let mut f3 = [0.0; 6];
                    for reconstructed in [false, true] {
                        let report2 = if reconstructed {
                            gas.characteristic_outlet_normal_flux_from_primitive(
                                planar,
                                normal,
                                reference_pressure,
                                normal * 0.3,
                                &mut f2,
                                || false,
                            )
                        } else {
                            gas.characteristic_outlet_normal_flux(
                                NormalCell {
                                    partial_densities: &q2[..2],
                                    normal_momentum: q2[2],
                                    transverse_momentum: q2[3],
                                    energy: q2[4],
                                },
                                normal,
                                reference_pressure,
                                normal * 0.3,
                                &mut f2,
                                || false,
                            )
                        }
                        .unwrap();
                        let report3 = if reconstructed {
                            gas.characteristic_outlet_spatial_flux_from_primitive(
                                spatial,
                                normal,
                                reference_pressure,
                                normal * 0.3,
                                &mut f3,
                                || false,
                            )
                        } else {
                            gas.characteristic_outlet_spatial_flux(
                                SpatialCell {
                                    partial_densities: &q3[..2],
                                    normal_momentum: q3[2],
                                    transverse_momentum: [q3[3], q3[4]],
                                    energy: q3[5],
                                },
                                normal,
                                reference_pressure,
                                normal * 0.3,
                                &mut f3,
                                || false,
                            )
                        }
                        .unwrap();
                        for i in 0..3 {
                            close(f2[i], expected[i]);
                            close(f3[i], expected[i]);
                        }
                        close(f2[3], mass_flux * shear[0]);
                        close(f3[3], mass_flux * shear[0]);
                        close(f3[4], mass_flux * shear[1]);
                        close(f2[4], expected[3] + 0.5 * mass_flux * shear[0].powi(2));
                        close(
                            f3[5],
                            expected[3] + 0.5 * mass_flux * (shear[0].powi(2) + shear[1].powi(2)),
                        );
                        close(report2.signal_speed, axial_report.signal_speed);
                        close(report3.signal_speed, axial_report.signal_speed);
                        for i in 0..2 {
                            if densities[i] == 0.0 {
                                assert_eq!(f2[i], 0.0);
                                assert_eq!(f3[i], 0.0);
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn g4_multidimensional_characteristic_refusals_preserve_complete_output() {
    let gas = variable_source(0.0);
    let original = [17.0; 6];
    let mut output = original;
    let primitive = SpatialPrimitive {
        partial_densities: &[0.4, 0.6],
        normal_velocity: 0.3,
        transverse_velocity: [2.0, -3.0],
        pressure: 1.0,
    };
    for (normal, pressure, velocity) in [(0.0, 1.0, 0.3), (1.0, 0.0, 0.3), (1.0, 1.0, -0.3)] {
        assert!(
            gas.characteristic_outlet_spatial_flux_from_primitive(
                primitive,
                normal,
                pressure,
                velocity,
                &mut output,
                || false,
            )
            .is_err()
        );
        assert_eq!(output, original);
    }
    for invalid in [
        SpatialPrimitive {
            normal_velocity: -0.3,
            ..primitive
        },
        SpatialPrimitive {
            transverse_velocity: [f64::NAN, 0.0],
            ..primitive
        },
        SpatialPrimitive {
            transverse_velocity: [1e20, 0.0],
            ..primitive
        },
    ] {
        assert!(
            gas.characteristic_outlet_spatial_flux_from_primitive(
                invalid,
                1.0,
                1.0,
                0.3,
                &mut output,
                || false,
            )
            .is_err()
        );
        assert_eq!(output, original);
    }
    // Inject cancellation after the complete characteristic solve reaches
    // its final publication checkpoint, not just at API entry.
    let mut calls = 0;
    gas.characteristic_outlet_spatial_flux_from_primitive(
        primitive,
        1.0,
        0.9,
        0.3,
        &mut output,
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    output = original;
    let mut actual = 0;
    assert_eq!(
        gas.characteristic_outlet_spatial_flux_from_primitive(
            primitive,
            1.0,
            0.9,
            0.3,
            &mut output,
            || {
                actual += 1;
                actual == calls
            },
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(output, original);
}

#[test]
fn g3_compact_multidimensional_ports_preserve_donor_species_and_stagnation_energy() {
    let gas = variable_source(1e18);
    for (densities, bath_densities) in [([1.0, 0.0], [0.0, 2.0]), ([0.25, 0.75], [0.7, 0.3])] {
        for normal in [-1.0, 1.0] {
            for (speed, bath_pressure) in [(0.4, 0.8), (0.4, 3.0), (0.0, 20.0), (0.4, 0.01)] {
                let velocity = normal * speed;
                let axial = gas.conserved(&densities, velocity, 1.0).unwrap();
                let bath = gas.conserved(&bath_densities, 0.0, bath_pressure).unwrap();
                let mut expected = [0.0; 4];
                let reference = gas
                    .reservoir_flux(
                        axial.as_cell(),
                        bath.as_cell(),
                        normal,
                        &mut expected,
                        || false,
                    )
                    .unwrap();
                let mass = expected[0] + expected[1];
                let incoming = normal * mass < 0.0;
                let donor = if incoming { bath_densities } else { densities };
                let donor_density = donor[0] + donor[1];
                for shear in [[2.0, -3.0], [-3.0, 2.0]] {
                    let planar = NormalPrimitive {
                        partial_densities: &densities,
                        normal_velocity: velocity,
                        transverse_velocity: shear[0],
                        pressure: 1.0,
                    };
                    let spatial = SpatialPrimitive {
                        partial_densities: &densities,
                        normal_velocity: velocity,
                        transverse_velocity: shear,
                        pressure: 1.0,
                    };
                    let mut q2 = [0.0; 5];
                    let mut q3 = [0.0; 6];
                    gas.normal_conserved_into(planar, &mut q2).unwrap();
                    gas.spatial_conserved_into(spatial, &mut q3).unwrap();
                    let mut f2 = [0.0; 5];
                    let mut f3 = [0.0; 6];
                    for reconstructed in [false, true] {
                        let report2 = if reconstructed {
                            gas.reservoir_normal_flux_from_primitive(
                                planar,
                                bath.as_cell(),
                                normal,
                                &mut f2,
                                || false,
                            )
                        } else {
                            gas.reservoir_normal_flux(
                                NormalCell {
                                    partial_densities: &q2[..2],
                                    normal_momentum: q2[2],
                                    transverse_momentum: q2[3],
                                    energy: q2[4],
                                },
                                bath.as_cell(),
                                normal,
                                &mut f2,
                                || false,
                            )
                        }
                        .unwrap();
                        let report3 = if reconstructed {
                            gas.reservoir_spatial_flux_from_primitive(
                                spatial,
                                bath.as_cell(),
                                normal,
                                &mut f3,
                                || false,
                            )
                        } else {
                            gas.reservoir_spatial_flux(
                                SpatialCell {
                                    partial_densities: &q3[..2],
                                    normal_momentum: q3[2],
                                    transverse_momentum: [q3[3], q3[4]],
                                    energy: q3[5],
                                },
                                bath.as_cell(),
                                normal,
                                &mut f3,
                                || false,
                            )
                        }
                        .unwrap();
                        for i in 0..3 {
                            close(f2[i], expected[i]);
                            close(f3[i], expected[i]);
                        }
                        for i in 0..2 {
                            close(f3[i], mass * donor[i] / donor_density);
                            if donor[i] == 0.0 {
                                assert_eq!(f2[i], 0.0);
                                assert_eq!(f3[i], 0.0);
                            }
                        }
                        let tangents = if incoming { [0.0; 2] } else { shear };
                        close(f2[3], mass * tangents[0]);
                        close(f3[3], mass * tangents[0]);
                        close(f3[4], mass * tangents[1]);
                        close(f2[4], expected[3] + 0.5 * mass * tangents[0].powi(2));
                        close(
                            f3[5],
                            expected[3] + 0.5 * mass * (tangents[0].powi(2) + tangents[1].powi(2)),
                        );
                        if incoming {
                            assert_eq!(f2[3], 0.0);
                            assert_eq!([f3[3], f3[4]], [0.0; 2]);
                            close(f3[5], mass * (bath.energy + bath_pressure) / donor_density);
                        }
                        close(report2.signal_speed, reference.signal_speed);
                        close(report3.signal_speed, reference.signal_speed);
                    }
                }
            }
        }
    }
}

#[test]
fn g4_compact_mixture_ports_publish_nothing_after_refusal_or_late_cancellation() {
    let gas = variable_source(0.0);
    let densities = [0.4, 0.6];
    let interior = SpatialPrimitive {
        partial_densities: &densities,
        normal_velocity: -0.3,
        transverse_velocity: [2.0, -3.0],
        pressure: 1.0,
    };
    let bath = gas.conserved(&[0.0, 2.0], 0.0, 3.0).unwrap();
    let original = [17.0; 6];
    let mut output = original;
    let moving = Cell {
        momentum: 0.1,
        ..bath.as_cell()
    };
    assert!(
        gas.reservoir_spatial_flux_from_primitive(interior, moving, 1.0, &mut output, || false,)
            .is_err()
    );
    assert_eq!(output, original);
    assert!(
        gas.reservoir_spatial_flux_from_primitive(
            interior,
            bath.as_cell(),
            1.0,
            &mut output[..5],
            || false,
        )
        .is_err()
    );
    assert_eq!(output, original);
    let unresolved = SpatialPrimitive {
        transverse_velocity: [1e20, 0.0],
        ..interior
    };
    assert!(
        gas.reservoir_spatial_flux_from_primitive(
            unresolved,
            bath.as_cell(),
            1.0,
            &mut output,
            || false,
        )
        .is_err()
    );
    assert_eq!(output, original);
    let planar = NormalPrimitive {
        partial_densities: &densities,
        normal_velocity: -0.3,
        transverse_velocity: 2.0,
        pressure: 1.0,
    };
    for spatial in [false, true] {
        let mut calls = 0;
        if spatial {
            gas.reservoir_spatial_flux_from_primitive(
                interior,
                bath.as_cell(),
                1.0,
                &mut output,
                || {
                    calls += 1;
                    false
                },
            )
            .unwrap();
        } else {
            gas.reservoir_normal_flux_from_primitive(
                planar,
                bath.as_cell(),
                1.0,
                &mut output[..5],
                || {
                    calls += 1;
                    false
                },
            )
            .unwrap();
        }
        output = original;
        let mut seen = 0;
        let result = if spatial {
            gas.reservoir_spatial_flux_from_primitive(
                interior,
                bath.as_cell(),
                1.0,
                &mut output,
                || {
                    seen += 1;
                    seen == calls
                },
            )
        } else {
            gas.reservoir_normal_flux_from_primitive(
                planar,
                bath.as_cell(),
                1.0,
                &mut output[..5],
                || {
                    seen += 1;
                    seen == calls
                },
            )
        };
        assert_eq!(result, Err(Error::Cancelled));
        assert_eq!(output, original);
    }
}

// Independent integration for Cp(T)=a+b*T with gas constant r, at fixed Y.
fn acoustic_integral(temperature: f64, a: f64, b: f64, r: f64) -> f64 {
    let panels = 512;
    let width = (temperature - 1.0) / f64::from(panels);
    let mut sum = 0.0;
    for i in 0..=panels {
        let t = 1.0 + f64::from(i) * width;
        let cp = a + b * t;
        let sound_speed = (cp / (cp - r) * r * t).sqrt();
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * cp / sound_speed;
    }
    width * sum / 3.0
}

#[test]
fn mixture_characteristic_outlet_passes_original_variable_cp_simple_waves_and_energy_datums() {
    let gas = variable_source(0.0);
    let shifted = variable_source(1e18);
    for fraction in [0.0, 0.25, 0.8, 1.0] {
        let r: f64 = fraction + (1.0 - fraction) * 0.5;
        let a = fraction * 3.5 + (1.0 - fraction) * 2.1;
        let b = 0.2;
        let c0 = ((a + b) / (a + b - r) * r).sqrt();
        let reference_velocity = 0.4 * c0;
        for normal in [-1.0, 1.0] {
            for temperature in [0.95_f64, 1.0, 1.05] {
                // Fixed-composition entropy: a ln(T)+b(T-1)-r ln(p)=0.
                let pressure = ((a * temperature.ln() + b * (temperature - 1.0)) / r).exp();
                let density = pressure / (r * temperature);
                let velocity =
                    normal * (reference_velocity + acoustic_integral(temperature, a, b, r));
                let state = gas
                    .conserved(
                        &[fraction * density, (1.0 - fraction) * density],
                        velocity,
                        pressure,
                    )
                    .unwrap();
                let mut flux = [0.0; 4];
                gas.characteristic_outlet_flux(
                    state.as_cell(),
                    normal,
                    1.0,
                    normal * reference_velocity,
                    &mut flux,
                    || false,
                )
                .unwrap();
                close_flux(&flux, physical_flux(&state, velocity, pressure));
                let mut shifted_flux = [0.0; 4];
                shifted
                    .characteristic_outlet_flux(
                        state.as_cell(),
                        normal,
                        1.0,
                        normal * reference_velocity,
                        &mut shifted_flux,
                        || false,
                    )
                    .unwrap();
                assert_eq!(flux, shifted_flux);
                if fraction == 0.0 {
                    assert_eq!(flux[0], 0.0);
                } else if fraction == 1.0 {
                    assert_eq!(flux[1], 0.0);
                }
            }
        }
    }
}

#[test]
fn mixture_characteristic_outlet_changes_only_the_incoming_acoustic_mode() {
    // Different molecular masses, same gamma=1.4. The independent invariant
    // is u +/- 5c; Cp/Cv and r depend on the CURRENT partial densities.
    let gas = standard(1e18);
    let fractions = [0.25, 0.75];
    let r = 0.625;
    let c0 = (1.4_f64 * r).sqrt();
    for normal in [-1.0, 1.0] {
        let reference_velocity = 0.2 * c0;
        let excess_velocity = 0.2 * c0;
        let state = gas
            .conserved(
                &[fractions[0] / r, fractions[1] / r],
                normal * (reference_velocity + excess_velocity),
                1.0,
            )
            .unwrap();
        let face_velocity = normal * (reference_velocity + 0.5 * excess_velocity);
        let face_sound = c0 + excess_velocity / 10.0;
        let temperature = (face_sound / c0).powi(2);
        let pressure = temperature.powf(3.5);
        let density = temperature.powf(2.5) / r;
        let mass = density * face_velocity;
        // Independent sensible-energy integral from this source's Tmin=.001.
        let energy = density * (2.5 * r * (temperature - 0.001) + 0.5 * face_velocity.powi(2));
        let expected = [
            fractions[0] * mass,
            fractions[1] * mass,
            mass * face_velocity + pressure,
            (energy + pressure) * face_velocity,
        ];
        let mut flux = [0.0; 4];
        gas.characteristic_outlet_flux(
            state.as_cell(),
            normal,
            1.0,
            normal * reference_velocity,
            &mut flux,
            || false,
        )
        .unwrap();
        close_flux(&flux, expected);
    }
}

#[test]
fn mixture_characteristic_outlet_preserves_entropy_and_composition_contacts() {
    let gas = variable_source(0.0);
    for normal in [-1.0, 1.0] {
        for densities in [[1.0, 0.0], [0.0, 2.0], [0.2, 0.8], [1.3, 2.7]] {
            let velocity = normal * 0.3;
            let state = gas.conserved(&densities, velocity, 1.0).unwrap();
            let mut flux = [0.0; 4];
            gas.characteristic_outlet_flux(
                state.as_cell(),
                normal,
                1.0,
                velocity,
                &mut flux,
                || false,
            )
            .unwrap();
            close_flux(&flux, physical_flux(&state, velocity, 1.0));
        }
    }
}

#[test]
fn mixture_characteristic_outlet_has_sonic_discharge_and_supersonic_causality() {
    let gas = standard(0.0);
    let r = 0.7;
    let c0 = (1.4_f64 * r).sqrt();
    for normal in [-1.0, 1.0] {
        let state = gas.conserved(&[0.4, 0.6], normal * 0.2 * c0, r).unwrap();
        let sonic_sound = (0.2 * c0 + 5.0 * c0) / 6.0;
        let temperature = (sonic_sound / c0).powi(2);
        let pressure = r * temperature.powf(3.5);
        let density = temperature.powf(2.5);
        let sonic = gas
            .conserved(
                &[0.4 * density, 0.6 * density],
                normal * sonic_sound,
                pressure,
            )
            .unwrap();
        let mut flux = [0.0; 4];
        gas.characteristic_outlet_flux(
            state.as_cell(),
            normal,
            0.001 * r,
            normal * 0.2 * c0,
            &mut flux,
            || false,
        )
        .unwrap();
        close_flux(&flux, physical_flux(&sonic, normal * sonic_sound, pressure));
        let supersonic = gas.conserved(&[0.4, 0.6], normal * 2.0 * c0, r).unwrap();
        for remote_pressure in [1e-100, 1e100] {
            gas.characteristic_outlet_flux(
                supersonic.as_cell(),
                normal,
                remote_pressure,
                normal * 0.2 * c0,
                &mut flux,
                || false,
            )
            .unwrap();
            close_flux(&flux, physical_flux(&supersonic, normal * 2.0 * c0, r));
        }
    }
}

#[test]
fn mixture_characteristic_outlet_refuses_backflow_and_cancellation_without_output_changes() {
    let gas = variable_source(0.0);
    let state = gas.conserved(&[0.4, 0.6], 0.3, 1.0).unwrap();
    let original = [17.0, -19.0, 23.0, -29.0];
    let mut flux = original;
    for (normal, pressure, velocity) in [
        (0.0, 1.0, 0.3),
        (1.0, 0.0, 0.3),
        (1.0, -1.0, 0.3),
        (1.0, f64::NAN, 0.3),
        (1.0, f64::INFINITY, 0.3),
        (1.0, 1.0, f64::NAN),
        (1.0, 1.0, f64::INFINITY),
        (1.0, 1.0, -0.3),
        (-1.0, 1.0, -0.3),
        (1.0, 100.0, 0.3),
    ] {
        assert!(
            gas.characteristic_outlet_flux(
                state.as_cell(),
                normal,
                pressure,
                velocity,
                &mut flux,
                || false,
            )
            .is_err()
        );
        assert_eq!(flux, original);
    }
    assert!(matches!(
        gas.characteristic_outlet_flux(state.as_cell(), 1.0, 1e-40, 0.3, &mut flux, || false,),
        Err(Error::CharacteristicPath(_))
    ));
    assert_eq!(flux, original);
    for after in [0, 12] {
        let mut calls = 0;
        let result =
            gas.characteristic_outlet_flux(state.as_cell(), 1.0, 0.8, 0.3, &mut flux, || {
                calls += 1;
                calls > after
            });
        assert_eq!(result, Err(Error::Cancelled));
        assert_eq!(flux, original);
    }
    assert!(
        gas.characteristic_outlet_flux(state.as_cell(), 1.0, 1.0, 0.3, &mut flux[..3], || false,)
            .is_err()
    );
    assert_eq!(flux, original);
}

fn acoustic_end(pressure_pa: f64, velocity_m_per_s: f64) -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Transmissive,
        right: Boundary::CharacteristicOutlet {
            pressure_pa,
            velocity_m_per_s,
        },
    }
}

#[test]
fn native_mixture_rk_outlet_tracks_advected_composition_without_pressure_reflection() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        // Equal constant Cp/R and different molecular masses give a linear,
        // pressure-compatible conservative contact, including reconstruction.
        // Variable-Cp local entropy/composition contacts are checked above.
        let gas = gas(vec![
            source("A", R, 3.5, 0.0, 1e18),
            source("B", 2.0 * R, 3.5, 0.0, -2e18),
        ]);
        let mut field = flow(gas, 48, acoustic_end(1.0, 0.6), |x| {
            let fraction = if x < 0.75 { 0.8 } else { 0.2 };
            let density = 1.0 / (fraction + 0.5 * (1.0 - fraction));
            ([fraction * density, (1.0 - fraction) * density], 0.6, 1.0)
        });
        let mut final_fraction = 0.0;
        while field.time().value() < 0.6 {
            let report = field
                .step(
                    Time::new(0.6 - field.time().value()),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    || false,
                )
                .unwrap();
            final_fraction = report.right_boundary_change.species_masses[0]
                / report.right_boundary_change.mass();
            for residual in &report.conservation_residual.species_masses {
                close(*residual, 0.0);
            }
            close(report.conservation_residual.momentum, 0.0);
            close(report.conservation_residual.energy, 0.0);
            for i in 0..field.cell_count() {
                let primitive = field
                    .gas()
                    .primitive(field.cell(i).unwrap(), || false)
                    .unwrap();
                close(primitive.pressure, 1.0);
                close(primitive.velocity, 0.6);
            }
        }
        assert!(
            final_fraction > 0.6,
            "outgoing species fraction={final_fraction}"
        );
    }
}

fn outgoing_pulse_remaining_pressure(right: Boundary) -> f64 {
    let mut field = flow(
        standard(0.0),
        80,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right,
        },
        |x| {
            let pressure = 1.0 + 1e-3 * (-((x - 0.65) / 0.055).powi(2)).exp();
            let density = pressure.powf(1.0 / 1.4);
            let sound_speed = (1.4 * pressure / density).sqrt();
            let velocity = 0.3 + 5.0 * (sound_speed - 1.4_f64.sqrt());
            ([0.4 * density, 0.6 * density], velocity, pressure)
        },
    );
    run(&mut field, 0.6);
    (0..field.cell_count())
        .map(|i| {
            (field
                .gas()
                .primitive(field.cell(i).unwrap(), || false)
                .unwrap()
                .pressure
                - 1.0)
                .powi(2)
        })
        .sum::<f64>()
        .sqrt()
}

#[test]
fn native_mixture_rk_acoustic_outlet_releases_pressure_pulse_without_fixed_pressure_reflection() {
    let acoustic = outgoing_pulse_remaining_pressure(Boundary::CharacteristicOutlet {
        pressure_pa: 1.0,
        velocity_m_per_s: 0.3,
    });
    let fixed = outgoing_pulse_remaining_pressure(Boundary::PressureOutlet { pressure_pa: 1.0 });
    assert!(fixed > 1e-4, "fixed-pressure reflection={fixed}");
    assert!(
        acoustic < 0.02 * fixed,
        "acoustic remainder={acoustic}, fixed reflection={fixed}"
    );
}

#[test]
fn native_mixture_rk_rejects_outlet_backflow_without_changing_the_accepted_field() {
    let mut field = flow(variable_source(0.0), 8, acoustic_end(100.0, 0.3), |_| {
        ([0.4, 0.6], 0.3, 1.0)
    });
    let original = field.states.clone();
    let time = field.time();
    assert!(
        field
            .step(
                Time::new(0.001),
                StepControl {
                    maximum_halvings: 0,
                    ..StepControl::default()
                },
                || false,
            )
            .is_err()
    );
    assert_eq!(field.states, original);
    assert_eq!(field.time(), time);
    field
        .set_boundaries(acoustic_end(1.0, 0.3), || false)
        .unwrap();
    assert_eq!(
        field.step(Time::new(0.001), StepControl::default(), || true),
        Err(Error::Cancelled)
    );
    assert_eq!(field.states, original);
    assert_eq!(field.time(), time);
    field
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
}
