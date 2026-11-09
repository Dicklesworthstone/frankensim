//! Public full-FE regressions. Native execution requires the complete crate
//! dependency closure; the physical response contractions run independently.

use super::super::resonance::{
    ResonanceSurveyControls, StructuralInteriorReceiver, StructuralInteriorSources,
};
use super::*;
use fs_modal::harmonic::HarmonicFrequencyControls;

fn frequency_controls() -> HarmonicFrequencyControls {
    HarmonicFrequencyControls {
        maximum_additional_payload_bytes: 4 * 1024 * 1024,
        maximum_work: 4_000_000_000,
    }
}

fn close_complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

fn quarter_energy(matrix: &Csr, field: &[C64]) -> f64 {
    0.25 * field
        .iter()
        .enumerate()
        .map(|(row, z)| {
            let mz = product(matrix, field, row);
            z.re * mz.re + z.im * mz.im
        })
        .sum::<f64>()
}

#[test]
fn original_nonuniform_spatial_tangent_closes_every_differentiated_row_and_work_channel() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = HeterogeneousHelmholtzFe::new(
        &geometry,
        5,
        &mut |_, p: [f64; 3]| {
            Ok::<_, AcousticFeError>(Medium {
                density_kg_m3: 1.1 + 0.2 * p[0],
                sound_speed_m_s: 3.0 + 0.4 * p[1],
            })
        },
        &mut || false,
    )
    .unwrap();
    let mass = dense(&[2.0, 0.4, 0.4, 1.1], 2);
    let stiffness = dense(&[70.0, -13.0, -13.0, 40.0], 2);
    let damping = dense(&[0.7, 0.3, -0.1, 0.5], 2);
    let traces = [
        vec![([0, 1, 2], [0.7, -0.3, 0.9])],
        vec![([2, 0, 1], [1.1, -0.2, 0.8])],
    ];
    let wall = ResistiveBoundary {
        faces: vec![[0, 1, 3]],
        admittance_m_pa_s: 0.08,
    };
    let surface = geometry.boundary_mass(&wall.faces).unwrap();
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &traces,
            &[wall],
            controls(),
            cx,
        )
        .unwrap();
        let force = [C64::new(0.7, -0.2), C64::new(-0.3, 0.4)];
        let force_slope = [C64::new(-0.02, 0.04), C64::new(0.03, -0.01)];
        let flow = model
            .project_inward_velocity(&[(
                [1, 2, 3],
                [
                    C64::new(0.1, 0.03),
                    C64::new(-0.04, 0.07),
                    C64::new(0.02, -0.01),
                ],
            )])
            .unwrap();
        let flow_slope = model
            .project_inward_velocity(&[(
                [1, 2, 3],
                [
                    C64::new(-0.01, 0.002),
                    C64::new(0.006, -0.003),
                    C64::new(0.004, 0.007),
                ],
            )])
            .unwrap();
        let omega = 3.7;
        let prepared = model.prepare_angular(omega, cx).unwrap();
        let result = prepared
            .solve_with_frequency_derivative(
                &force,
                &flow,
                Some(StructuralInteriorSources {
                    structural_force: &force_slope,
                    inward_flow_load: &flow_slope,
                }),
                frequency_controls(),
                cx,
            )
            .unwrap();
        let r = &result.response;
        let du = &result.structural_displacement_derivative_per_rad_s;
        let dv = &result.structural_velocity_derivative_per_rad_s;
        let dpsi = &result.pressure_potential_derivative_pa_s_per_rad_s;
        let dp = &result.pressure_derivative_pa_per_rad_s;
        for row in 0..2 {
            let residual = product(&stiffness, du, row)
                - product(&mass, du, row).scale(omega * omega)
                - product(&mass, &r.structural_displacement, row).scale(2.0 * omega)
                + product(&damping, dv, row)
                - result.fluid_force_derivative_per_rad_s[row]
                - force_slope[row];
            assert!(
                residual.abs() < 2e-10,
                "structural derivative row {row}: {residual:?}"
            );
            close_complex(
                dv[row],
                C64::new(0.0, -1.0) * (r.structural_displacement[row] + du[row].scale(omega)),
                1e-13,
            );
        }
        for row in 0..geometry.ndof() {
            let residual = product(interior.stiffness(), dpsi, row)
                - product(interior.mass(), dpsi, row).scale(omega * omega)
                - product(interior.mass(), &r.pressure_potential_pa_s, row).scale(2.0 * omega)
                + product(&surface, dp, row).scale(0.08)
                + result.displaced_flow_derivative_m3_s_per_rad_s[row]
                - flow_slope[row];
            assert!(
                residual.abs() < 2e-10,
                "acoustic derivative row {row}: {residual:?}"
            );
            close_complex(
                dp[row],
                C64::new(0.0, -1.0) * (r.pressure_potential_pa_s[row] + dpsi[row].scale(omega)),
                1e-13,
            );
        }
        for (value, expected) in [
            (
                r.energy.structural_kinetic_j,
                quarter_energy(&mass, &r.structural_velocity),
            ),
            (
                r.energy.structural_strain_j,
                quarter_energy(&stiffness, &r.structural_displacement),
            ),
            (
                r.energy.fluid_kinetic_j,
                quarter_energy(interior.stiffness(), &r.pressure_potential_pa_s),
            ),
            (
                r.energy.fluid_compression_j,
                quarter_energy(interior.mass(), &r.pressure_pa),
            ),
        ] {
            assert!((value - expected).abs() < 2e-13 * expected.abs().max(1.0));
        }
        assert!(r.structural_power_balance_residual_w().abs() < 1e-10);
        assert!(r.acoustic_power_balance_residual_w().abs() < 1e-10);
        assert!(
            result
                .power_derivative
                .structural_balance_residual_w_per_rad_s()
                .abs()
                < 1e-10
        );
        assert!(
            result
                .power_derivative
                .acoustic_balance_residual_w_per_rad_s()
                .abs()
                < 1e-10
        );
        assert!(result.power_balance_residual_derivative_w_per_rad_s.abs() < 1e-10);

        // Independent nearby whole-field solves include the authored source change.
        let h = 1e-5;
        let at = |delta: f64| {
            let f: Vec<_> = force
                .iter()
                .zip(force_slope)
                .map(|(&f, df)| f + df.scale(delta))
                .collect();
            let s: Vec<_> = flow
                .iter()
                .zip(&flow_slope)
                .map(|(&s, &ds)| s + ds.scale(delta))
                .collect();
            model
                .prepare_angular(omega + delta, cx)
                .unwrap()
                .solve(&f, &s, cx)
                .unwrap()
        };
        let low = at(-h);
        let high = at(h);
        for (tangent, (a, b)) in du
            .iter()
            .zip(
                low.structural_displacement
                    .iter()
                    .zip(&high.structural_displacement),
            )
            .chain(dp.iter().zip(low.pressure_pa.iter().zip(&high.pressure_pa)))
            .chain(
                result
                    .fluid_force_derivative_per_rad_s
                    .iter()
                    .zip(low.fluid_force.iter().zip(&high.fluid_force)),
            )
        {
            close_complex(*tangent, (*b - *a).scale(0.5 / h), 2e-7);
        }
        for (actual, a, b) in [
            (
                result.power_derivative.structural_input_w_per_rad_s,
                low.average_structural_input_power_w,
                high.average_structural_input_power_w,
            ),
            (
                result.power_derivative.acoustic_input_w_per_rad_s,
                low.average_acoustic_input_power_w,
                high.average_acoustic_input_power_w,
            ),
            (
                result.power_derivative.fluid_to_structure_w_per_rad_s,
                low.average_fluid_to_structure_power_w,
                high.average_fluid_to_structure_power_w,
            ),
            (
                result.power_derivative.structural_dissipation_w_per_rad_s,
                low.average_structural_dissipation_w,
                high.average_structural_dissipation_w,
            ),
            (
                result.power_derivative.wall_dissipation_w_per_rad_s,
                low.average_wall_dissipation_w,
                high.average_wall_dissipation_w,
            ),
        ] {
            assert!((actual - (b - a) / (2.0 * h)).abs() < 2e-7 * actual.abs().max(1.0));
        }
    });
}

#[test]
fn full_spatial_coupling_moves_resonance_and_each_physical_receiver_keeps_its_phase() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[0.8], 1);
    let stiffness = dense(&[4.0], 1);
    let damping = dense(&[0.3], 1);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[vec![([0, 1, 2], [1.0; 3])]],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let grid: Vec<_> = (0..=32).map(|i| 2.0 + 0.2 * i as f64).collect();
        let force = [C64::new(1.0, 0.3)];
        let zero = vec![C64::ZERO; geometry.ndof()];
        let search = ResonanceSurveyControls {
            absolute_frequency_tolerance_rad_s: 1e-7,
            relative_frequency_tolerance: 1e-8,
            maximum_evaluations: 384,
            maximum_peaks: 8,
        };
        let receivers = [
            StructuralInteriorReceiver::StructuralDisplacement(&[1.0]),
            StructuralInteriorReceiver::StructuralVelocity(&[1.0]),
            StructuralInteriorReceiver::PressurePoint {
                tetrahedron: 0,
                barycentric: [1.0, 0.0, 0.0, 0.0],
            },
        ];
        let mut peaks = Vec::new();
        for (kind, receiver) in receivers.into_iter().enumerate() {
            let survey = model
                .survey_resonances(
                    &grid,
                    &force,
                    &zero,
                    receiver,
                    frequency_controls(),
                    search,
                    cx,
                )
                .unwrap();
            assert!(!survey.peaks.is_empty());
            let peak = &survey.peaks[0];
            assert!(peak.half_power.is_none());
            let omega = peak.peak.angular_frequency_rad_s;
            assert!(omega > 3.0 && omega < 7.0);
            let observed = |omega| {
                let r = model
                    .prepare_angular(omega, cx)
                    .unwrap()
                    .solve(&force, &zero, cx)
                    .unwrap();
                match kind {
                    0 => r.structural_displacement[0],
                    1 => r.structural_velocity[0],
                    _ => interior
                        .pressure_at(&r.pressure_pa, 0, [1.0, 0.0, 0.0, 0.0])
                        .unwrap(),
                }
            };
            close_complex(peak.observation, observed(omega), 1e-10);
            for offset in [-0.02, 0.02] {
                assert!(observed(omega + offset).abs() < peak.observation.abs());
            }
            peaks.push(omega);
        }
        // The original closed cavity changes the response; velocity has its
        // own stationary amplitude, not the displacement peak relabelled.
        let dry_omega = fs_math::det::sqrt(4.0 / 0.8);
        assert!(peaks[0] > dry_omega + 0.5);
        assert!(peaks[1] > peaks[0] + 1e-5);
    });
}

#[test]
fn zero_primary_drive_keeps_source_tangent_and_refusals_preserve_the_factor() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 2.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[10.0], 1);
    let damping = dense(&[0.8], 1);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[vec![([0, 1, 2], [1.0; 3])]],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let prepared = model.prepare_angular(2.1, cx).unwrap();
        let force = [C64::new(0.8, 0.3)];
        let zero = [C64::ZERO; 4];
        let ordinary = prepared.solve(&force, &zero, cx).unwrap();
        let derivative = prepared
            .solve_with_frequency_derivative(
                &[C64::ZERO],
                &zero,
                Some(StructuralInteriorSources {
                    structural_force: &force,
                    inward_flow_load: &zero,
                }),
                frequency_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(derivative.response.structural_displacement, [C64::ZERO]);
        assert_eq!(derivative.response.energy.total_j(), 0.0);
        for (actual, expected) in derivative
            .structural_displacement_derivative_per_rad_s
            .iter()
            .zip(&ordinary.structural_displacement)
            .chain(
                derivative
                    .pressure_derivative_pa_per_rad_s
                    .iter()
                    .zip(&ordinary.pressure_pa),
            )
        {
            close_complex(*actual, *expected, 1e-12);
        }
        let too_small = HarmonicFrequencyControls {
            maximum_additional_payload_bytes: 0,
            ..frequency_controls()
        };
        assert!(matches!(
            prepared.solve_with_frequency_derivative(&force, &zero, None, too_small, cx),
            Err(StructuralInteriorError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        assert!(
            prepared
                .solve_with_frequency_derivative(
                    &force,
                    &zero,
                    Some(StructuralInteriorSources {
                        structural_force: &[],
                        inward_flow_load: &zero
                    }),
                    frequency_controls(),
                    cx
                )
                .is_err()
        );
        let replay = prepared.solve(&force, &zero, cx).unwrap();
        assert_eq!(
            ordinary.structural_displacement,
            replay.structural_displacement
        );
        assert_eq!(ordinary.pressure_pa, replay.pressure_pa);
        assert_eq!(ordinary.energy, replay.energy);
    });
}

#[test]
fn original_fe_interface_direction_reaches_forced_fields_and_mixed_work() {
    use super::super::resonance::{
        StructuralInteriorParameterDirection, StructuralInteriorResponseField as Field,
        StructuralInteriorSourceParameterJet,
    };
    use fs_modal::harmonic::HarmonicFrequencyParameterControls;

    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 2.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[10.0], 1);
    let damping = dense(&[0.8], 1);
    let limits = HarmonicFrequencyParameterControls {
        maximum_additional_payload_bytes: 16 * 1024 * 1024,
        maximum_work: 4_000_000_000,
    };
    with_cx(|cx| {
        let make = |scale: f64| {
            StructuralInterior::new(
                &interior,
                StructuralMatrices {
                    mass: &mass,
                    stiffness: &stiffness,
                    damping: Some(&damping),
                },
                &[vec![([0, 1, 2], [scale, -0.4 * scale, 0.7 * scale])]],
                &[ResistiveBoundary {
                    faces: vec![[0, 1, 3]],
                    admittance_m_pa_s: 0.08,
                }],
                controls(),
                cx,
            )
            .unwrap()
        };
        let model = make(1.0);
        let omega = 2.1;
        let prepared = model.prepare_angular(omega, cx).unwrap();
        let force = [C64::new(0.8, 0.3)];
        let flow = [C64::new(0.02, -0.01), C64::ZERO, C64::ZERO, C64::ZERO];
        let source = StructuralInteriorSourceParameterJet::fixed(StructuralInteriorSources {
            structural_force: &force,
            inward_flow_load: &flow,
        });
        let zero_structure = [0.0];
        let zero_fluid = vec![0.0; model.pressure_dofs().pow(2)];
        let direction = StructuralInteriorParameterDirection {
            structural_mass: &zero_structure,
            structural_damping: &zero_structure,
            structural_stiffness: &zero_structure,
            acoustic_mass: &zero_fluid,
            acoustic_stiffness: &zero_fluid,
            wall_damping: &zero_fluid,
            coupling: model.coupling(),
        };
        let result = prepared
            .solve_with_frequency_parameter_derivative(source, direction, limits, cx)
            .unwrap();
        let original = prepared.solve(&force, &flow, cx).unwrap();
        let step = 1e-5;
        let low_model = make(1.0 - step);
        let high_model = make(1.0 + step);
        let low = low_model
            .prepare_angular(omega, cx)
            .unwrap()
            .solve_with_frequency_derivative(&force, &flow, None, frequency_controls(), cx)
            .unwrap();
        let high = high_model
            .prepare_angular(omega, cx)
            .unwrap()
            .solve_with_frequency_derivative(&force, &flow, None, frequency_controls(), cx)
            .unwrap();
        for (kind, original, low, high, low_w, high_w) in [
            (
                Field::StructuralDisplacement,
                original.structural_displacement.as_slice(),
                low.response.structural_displacement.as_slice(),
                high.response.structural_displacement.as_slice(),
                low.structural_displacement_derivative_per_rad_s.as_slice(),
                high.structural_displacement_derivative_per_rad_s.as_slice(),
            ),
            (
                Field::StructuralVelocity,
                original.structural_velocity.as_slice(),
                low.response.structural_velocity.as_slice(),
                high.response.structural_velocity.as_slice(),
                low.structural_velocity_derivative_per_rad_s.as_slice(),
                high.structural_velocity_derivative_per_rad_s.as_slice(),
            ),
            (
                Field::PressurePotential,
                original.pressure_potential_pa_s.as_slice(),
                low.response.pressure_potential_pa_s.as_slice(),
                high.response.pressure_potential_pa_s.as_slice(),
                low.pressure_potential_derivative_pa_s_per_rad_s.as_slice(),
                high.pressure_potential_derivative_pa_s_per_rad_s.as_slice(),
            ),
            (
                Field::Pressure,
                original.pressure_pa.as_slice(),
                low.response.pressure_pa.as_slice(),
                high.response.pressure_pa.as_slice(),
                low.pressure_derivative_pa_per_rad_s.as_slice(),
                high.pressure_derivative_pa_per_rad_s.as_slice(),
            ),
            (
                Field::FluidForce,
                original.fluid_force.as_slice(),
                low.response.fluid_force.as_slice(),
                high.response.fluid_force.as_slice(),
                low.fluid_force_derivative_per_rad_s.as_slice(),
                high.fluid_force_derivative_per_rad_s.as_slice(),
            ),
            (
                Field::DisplacedFlow,
                original.displaced_flow_m3_s.as_slice(),
                low.response.displaced_flow_m3_s.as_slice(),
                high.response.displaced_flow_m3_s.as_slice(),
                low.displaced_flow_derivative_m3_s_per_rad_s.as_slice(),
                high.displaced_flow_derivative_m3_s_per_rad_s.as_slice(),
            ),
        ] {
            let jet = result.field_parameter_jet(kind);
            for i in 0..original.len() {
                close_complex(jet.values[i], original[i], 2e-11);
                close_complex(
                    jet.parameter_derivative[i],
                    (high[i] - low[i]).scale(0.5 / step),
                    2e-7,
                );
                close_complex(
                    jet.angular_parameter_derivative[i],
                    (high_w[i] - low_w[i]).scale(0.5 / step),
                    3e-7,
                );
            }
        }
        for (jet, low, high, low_w, high_w) in [
            (
                result.power.structural_input,
                low.response.average_structural_input_power_w,
                high.response.average_structural_input_power_w,
                low.power_derivative.structural_input_w_per_rad_s,
                high.power_derivative.structural_input_w_per_rad_s,
            ),
            (
                result.power.acoustic_input,
                low.response.average_acoustic_input_power_w,
                high.response.average_acoustic_input_power_w,
                low.power_derivative.acoustic_input_w_per_rad_s,
                high.power_derivative.acoustic_input_w_per_rad_s,
            ),
            (
                result.power.fluid_to_structure,
                low.response.average_fluid_to_structure_power_w,
                high.response.average_fluid_to_structure_power_w,
                low.power_derivative.fluid_to_structure_w_per_rad_s,
                high.power_derivative.fluid_to_structure_w_per_rad_s,
            ),
            (
                result.power.structural_dissipation,
                low.response.average_structural_dissipation_w,
                high.response.average_structural_dissipation_w,
                low.power_derivative.structural_dissipation_w_per_rad_s,
                high.power_derivative.structural_dissipation_w_per_rad_s,
            ),
            (
                result.power.wall_dissipation,
                low.response.average_wall_dissipation_w,
                high.response.average_wall_dissipation_w,
                low.power_derivative.wall_dissipation_w_per_rad_s,
                high.power_derivative.wall_dissipation_w_per_rad_s,
            ),
        ] {
            assert!((jet.parameter_derivative - (high - low) * (0.5 / step)).abs() < 2e-7);
            assert!(
                (jet.angular_parameter_derivative - (high_w - low_w) * (0.5 / step)).abs() < 3e-7
            );
        }
        for balance in [
            result.power.structural_balance().unwrap(),
            result.power.acoustic_balance().unwrap(),
        ] {
            assert!(balance.sample.average_power_w.abs() < 2e-11);
            assert!(balance.sample.angular_derivative.abs() < 2e-11);
            assert!(balance.angular_second_derivative.abs() < 2e-10);
            assert!(balance.parameter_derivative.abs() < 2e-11);
            assert!(balance.angular_parameter_derivative.abs() < 2e-10);
        }
        let too_small = HarmonicFrequencyParameterControls {
            maximum_additional_payload_bytes: 0,
            ..limits
        };
        assert!(matches!(
            prepared.solve_with_frequency_parameter_derivative(source, direction, too_small, cx),
            Err(StructuralInteriorError::Harmonic(
                HarmonicResponseError::Budget { .. }
            ))
        ));
        let replay = prepared.solve(&force, &flow, cx).unwrap();
        assert_eq!(
            original.structural_displacement,
            replay.structural_displacement
        );
        assert_eq!(original.pressure_pa, replay.pressure_pa);
    });
}
