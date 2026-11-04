//! G1/G3: original material waves, physical sources, cutoff limits and work.

use super::*;
use crate::modal::ModalSourceResponse;

fn pressure(open: C64, impedance: C64) -> PortExcitation {
    PortExcitation::Pressure {
        open_circuit_pressure_pa: open,
        impedance_pa_s_m3: impedance,
    }
}

fn balanced(response: &ModalSourceResponse) {
    let scale = response
        .generator_powers_w
        .iter()
        .map(|p| p.abs())
        .sum::<f64>()
        + response.source_dissipation_w.abs()
        + response.inlet_power_w.abs();
    assert!(
        (response.generator_power_w - response.source_dissipation_w - response.inlet_power_w).abs()
            <= 2e-10 * scale.max(1e-30)
    );
    assert!(
        response
            .source_dissipations_w
            .iter()
            .all(|&power| power >= 0.0)
    );
}

#[test]
fn material_native_waves_and_source_drives_retain_original_gas_results() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let source = duct(&[0.014, 0.02, 0.017], 0.011);
    let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
    let reference = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.014 * 0.014);
    let sources = [
        pressure(C64::new(1.0, 0.3), C64::new(0.5, 0.2).scale(reference)),
        PortExcitation::VolumeFlow(C64::new(0.2, -0.1).scale(1.0 / reference)),
        pressure(C64::new(-0.2, 0.7), C64::new(0.3, -0.1).scale(reference)),
    ];
    for omega in [5000.0, 60_000.0] {
        for boundary in [
            MouthBoundary::Rigid,
            MouthBoundary::PressureRelease,
            MouthBoundary::Flanged,
        ] {
            for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
                let native = if boundary == MouthBoundary::Flanged {
                    geometry
                        .mouth_waves_flanged(&gas, &gas, omega, loss)
                        .unwrap()
                } else {
                    let termination = if boundary == MouthBoundary::Rigid {
                        Termination::Closed
                    } else {
                        Termination::IdealOpen
                    };
                    geometry
                        .mouth_waves(&gas, omega, loss, termination)
                        .unwrap()
                };
                let shared = geometry
                    .mouth_waves_from_material_with_loss(&gas, omega, loss, boundary)
                    .unwrap();
                let reflection = geometry
                    .input_reflection_from_material_with_loss(&gas, omega, loss, boundary)
                    .unwrap();
                assert_eq!(
                    shared.input.reflection_matrix,
                    native.input.reflection_matrix
                );
                assert_eq!(reflection.reflection_matrix, native.input.reflection_matrix);
                assert_eq!(
                    shared.input.characteristic_impedances,
                    native.input.characteristic_impedances
                );
                assert_eq!(shared.native_flow_map, native.native_flow_map);
                let expected = native.drive_sources(&sources).unwrap();
                let actual = geometry
                    .drive_sources_from_material_with_loss(&gas, omega, loss, boundary, &sources)
                    .unwrap();
                assert_eq!(actual.inlet.inlet_pressure, expected.inlet.inlet_pressure);
                assert_eq!(actual.inlet.inlet_flow, expected.inlet.inlet_flow);
                assert_eq!(actual.mouth_flow, expected.mouth_flow);
                assert_eq!(
                    actual.inlet.generator_power_w,
                    expected.inlet.generator_power_w
                );
                balanced(&actual.inlet);
                if loss == LossModel::Lossless {
                    let acoustic = geometry
                        .mouth_waves_from_material(&gas, omega, boundary)
                        .unwrap();
                    let reflection = geometry
                        .input_reflection_from_material(&gas, omega, boundary)
                        .unwrap();
                    let driven = geometry
                        .drive_sources_from_material(&gas, omega, boundary, &sources)
                        .unwrap();
                    assert_eq!(acoustic.native_flow_map, native.native_flow_map);
                    assert_eq!(reflection.reflection_matrix, native.input.reflection_matrix);
                    assert_eq!(driven.inlet.inlet_pressure, expected.inlet.inlet_pressure);
                    assert_eq!(driven.inlet.inlet_flow, expected.inlet.inlet_flow);
                    assert_eq!(driven.mouth_flow, expected.mouth_flow);
                }
            }
        }
    }
}

#[test]
fn physical_liquid_sources_survive_exact_cutoff_poles_without_exposing_native_references() {
    let state = liquid();
    let radius = j1_roots(1).unwrap()[0] / 128.0;
    let length = 0.009;
    let omega = state.sound_speed_m_s() * 128.0;
    let source = duct(&[radius], length);
    let geometry = PreparedModalDuct::new(&source, 2, 1).unwrap();
    let area = core::f64::consts::PI * radius * radius;
    let reference = state.characteristic_impedance_pa_s_m() / area;
    let open = C64::new(0.3, -0.2);
    let impedance = C64::new(0.4, 0.1).scale(reference);
    let sources = [
        PortExcitation::VolumeFlow(C64::ZERO),
        pressure(open, impedance),
    ];

    // The higher-mode input impedance really is a pole: its rigid-mouth
    // cutoff field has U(x)=0. A finite pressure drive still has P=p_open.
    assert!(matches!(
        geometry.input_impedance_from_material(&state, omega, MouthBoundary::Rigid),
        Err(FluidDuctError::Duct(DuctError::Singular))
    ));
    let closed = geometry
        .drive_sources_from_material(&state, omega, MouthBoundary::Rigid, &sources)
        .unwrap();
    assert_eq!(closed.inlet.inlet_flow, vec![C64::ZERO; 2]);
    assert_eq!(closed.mouth_flow, vec![C64::ZERO; 2]);
    assert_eq!(closed.inlet.inlet_pressure[0], C64::ZERO);
    near(closed.inlet.inlet_pressure[1], open, open.abs());
    assert_eq!(closed.inlet.generator_power_w, 0.0);
    assert_eq!(closed.inlet.source_dissipation_w, 0.0);
    let quiet = geometry
        .drive_sources_from_material(
            &state,
            omega,
            MouthBoundary::Rigid,
            &[PortExcitation::VolumeFlow(C64::ZERO); 2],
        )
        .unwrap();
    assert_eq!(quiet.inlet.inlet_pressure, vec![C64::ZERO; 2]);

    // Prescribing nonzero flow into that exact closed cutoff is inconsistent;
    // no resistance or finite native characteristic impedance is invented.
    assert!(matches!(
        geometry.drive_sources_from_material(
            &state,
            omega,
            MouthBoundary::Rigid,
            &[
                PortExcitation::VolumeFlow(C64::ZERO),
                PortExcitation::VolumeFlow(C64::from_re(1.0 / reference))
            ]
        ),
        Err(FluidDuctError::Duct(DuctError::Singular))
    ));
    for boundary in [
        MouthBoundary::Rigid,
        MouthBoundary::PressureRelease,
        MouthBoundary::Flanged,
    ] {
        assert!(matches!(
            geometry.input_reflection_from_material(&state, omega, boundary),
            Err(FluidDuctError::Duct(DuctError::Singular))
        ));
        assert!(matches!(
            geometry.mouth_waves_from_material(&state, omega, boundary),
            Err(FluidDuctError::Duct(DuctError::Singular))
        ));
    }

    // Pressure release gives the independent constant-flow cutoff solution:
    // P_in=-i omega rho L U/S, including its source resistance work.
    let inertance = C64::new(0.0, -omega * state.density_kg_m3() * length / area);
    let expected_flow = open / (impedance + inertance);
    let response = geometry
        .drive_sources_from_material(&state, omega, MouthBoundary::PressureRelease, &sources)
        .unwrap();
    near(
        response.inlet.inlet_pressure[1],
        inertance * expected_flow,
        open.abs(),
    );
    near(
        response.inlet.inlet_flow[1],
        expected_flow,
        expected_flow.abs(),
    );
    near(response.mouth_flow[1], expected_flow, expected_flow.abs());
    near(
        C64::from_re(response.inlet.source_dissipation_w),
        C64::from_re(0.5 * impedance.re * expected_flow.norm_sq()),
        response.inlet.generator_power_w.abs(),
    );
    balanced(&response.inlet);
    for factor in [1.0 - 1e-7, 1.0 + 1e-7] {
        let nearby = geometry
            .drive_sources_from_material(
                &state,
                omega * factor,
                MouthBoundary::PressureRelease,
                &sources,
            )
            .unwrap();
        assert!((nearby.inlet.inlet_flow[1] / expected_flow - C64::ONE).abs() < 1e-6);
    }
}

// Independent physical source closure applied to the existing global
// pressure/flow ODE oracle; no native-wave or production source helper.
fn source_oracle(
    z: &[C64],
    map: &[C64],
    sources: &[PortExcitation],
) -> (Vec<C64>, Vec<C64>, Vec<C64>) {
    let n = sources.len();
    let mut boundary = z.to_vec();
    let mut rhs = Vec::with_capacity(n);
    for (i, &source) in sources.iter().enumerate() {
        match source {
            PortExcitation::VolumeFlow(flow) => {
                boundary[i * n..(i + 1) * n].fill(C64::ZERO);
                boundary[i * n + i] = C64::ONE;
                rhs.push(flow);
            }
            PortExcitation::Pressure {
                open_circuit_pressure_pa,
                impedance_pa_s_m3,
            } => {
                boundary[i * n + i] = boundary[i * n + i] + impedance_pa_s_m3;
                rhs.push(open_circuit_pressure_pa);
            }
        }
        let scale = boundary[i * n..(i + 1) * n]
            .iter()
            .map(|z| z.abs())
            .fold(0.0_f64, f64::max);
        for z in &mut boundary[i * n..(i + 1) * n] {
            *z = z.scale(1.0 / scale);
        }
        rhs[i] = rhs[i].scale(1.0 / scale);
    }
    lu_complex(&boundary, n).unwrap().solve(&mut rhs);
    (apply(z, &rhs), apply(map, &rhs), rhs)
}

#[test]
fn coherent_material_sources_match_global_fields_and_full_radiated_power() {
    let state = liquid();
    let cutoff_radius = j1_roots(1).unwrap()[0] / 128.0;
    let source = duct(&[cutoff_radius, 0.045, 0.038], 0.009);
    let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
    let omega = 128.0 * state.sound_speed_m_s();
    let reference = state.characteristic_impedance_pa_s_m()
        / (core::f64::consts::PI * cutoff_radius * cutoff_radius);
    let sources = [
        pressure(C64::new(2.0, 0.3), C64::new(0.6, 0.2).scale(reference)),
        PortExcitation::VolumeFlow(C64::new(0.25, -0.1).scale(1.0 / reference)),
        pressure(C64::new(0.4, 0.8), C64::new(0.3, -0.2).scale(reference)),
    ];
    let radiation = crate::flanged::flanged_modal_impedance_from_material(
        &state,
        geometry.mouth_radius_m(),
        omega,
        3,
    )
    .unwrap();
    let (z, map) =
        direct_pressure_flow(&geometry, &state, omega, MouthBoundary::Flanged, &radiation);
    let (expected_pressure, expected_mouth, expected_flow) = source_oracle(&z, &map, &sources);
    let actual = geometry
        .drive_sources_from_material(&state, omega, MouthBoundary::Flanged, &sources)
        .unwrap();
    for i in 0..3 {
        near(actual.inlet.inlet_pressure[i], expected_pressure[i], 2.0);
        near(
            actual.inlet.inlet_flow[i],
            expected_flow[i],
            2.0 / reference,
        );
        near(actual.mouth_flow[i], expected_mouth[i], 2.0 / reference);
    }
    let radiated = power(&apply(&radiation, &actual.mouth_flow), &actual.mouth_flow);
    assert!(radiated > 0.0);
    near(
        C64::from_re(actual.inlet.inlet_power_w),
        C64::from_re(radiated),
        radiated,
    );
    balanced(&actual.inlet);

    // Preserve all phases and zero-drive impedances in a superposition check.
    let quiet = |source| match source {
        PortExcitation::VolumeFlow(_) => PortExcitation::VolumeFlow(C64::ZERO),
        PortExcitation::Pressure {
            impedance_pa_s_m3, ..
        } => pressure(C64::ZERO, impedance_pa_s_m3),
    };
    let a = geometry
        .drive_sources_from_material(
            &state,
            omega,
            MouthBoundary::Flanged,
            &[sources[0], quiet(sources[1]), quiet(sources[2])],
        )
        .unwrap();
    let b = geometry
        .drive_sources_from_material(
            &state,
            omega,
            MouthBoundary::Flanged,
            &[quiet(sources[0]), sources[1], sources[2]],
        )
        .unwrap();
    for i in 0..3 {
        near(
            actual.inlet.inlet_pressure[i],
            a.inlet.inlet_pressure[i] + b.inlet.inlet_pressure[i],
            2.0,
        );
        near(
            actual.mouth_flow[i],
            a.mouth_flow[i] + b.mouth_flow[i],
            2.0 / reference,
        );
    }
    assert!(
        (actual.inlet.inlet_power_w - a.inlet.inlet_power_w - b.inlet.inlet_power_w).abs()
            > 1e-4 * actual.inlet.inlet_power_w.abs()
    );

    let transport =
        FluidTransportState::try_new(state.density_kg_m3(), 4200.0, 0.6, 0.001).unwrap();
    let material = ThermoviscousFluidState::try_new(state, transport, 4100.0).unwrap();
    let lossy = geometry
        .drive_sources_from_material_with_loss(
            &material,
            omega,
            LossModel::Bessel,
            MouthBoundary::Flanged,
            &sources,
        )
        .unwrap();
    let radiated = power(&apply(&radiation, &lossy.mouth_flow), &lossy.mouth_flow);
    assert!(lossy.inlet.inlet_power_w > radiated);
    balanced(&lossy.inlet);
}

#[test]
fn source_queries_keep_one_material_snapshot_and_original_refusals() {
    let state = liquid();
    let transport =
        FluidTransportState::try_new(state.density_kg_m3(), 4200.0, 0.6, 0.001).unwrap();
    let material = CountedFluid {
        state: ThermoviscousFluidState::try_new(state, transport, 4100.0).unwrap(),
        acoustic_samples: Cell::new(0),
        complete_samples: Cell::new(0),
    };
    let source = duct(&[0.025, 0.038, 0.028], 0.009);
    let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
    let omega = state.sound_speed_m_s() * 110.0;
    let boundary = MouthBoundary::Flanged;
    let sources = [
        pressure(C64::ONE, C64::from_re(1e9)),
        PortExcitation::VolumeFlow(C64::ZERO),
        PortExcitation::VolumeFlow(C64::ZERO),
    ];
    geometry
        .input_reflection_from_material(&material, omega, boundary)
        .unwrap();
    geometry
        .mouth_waves_from_material(&material, omega, boundary)
        .unwrap();
    geometry
        .drive_sources_from_material(&material, omega, boundary, &sources)
        .unwrap();
    geometry
        .input_reflection_from_material_with_loss(&material, omega, LossModel::Bessel, boundary)
        .unwrap();
    geometry
        .mouth_waves_from_material_with_loss(&material, omega, LossModel::Bessel, boundary)
        .unwrap();
    geometry
        .drive_sources_from_material_with_loss(
            &material,
            omega,
            LossModel::Bessel,
            boundary,
            &sources,
        )
        .unwrap();
    assert_eq!(material.acoustic_samples.get(), 3);
    assert_eq!(material.complete_samples.get(), 3);
    assert!(
        geometry
            .drive_sources_from_material(&material, 0.0, boundary, &sources)
            .is_err()
    );
    assert_eq!(material.acoustic_samples.get(), 3);
    assert!(
        geometry
            .drive_sources_from_material(&state, omega, boundary, &[])
            .is_err()
    );
    let mut bad = sources;
    bad[0] = pressure(C64::ONE, C64::from_re(-1.0));
    assert!(
        geometry
            .drive_sources_from_material(&state, omega, boundary, &bad)
            .is_err()
    );
    bad[0] = PortExcitation::VolumeFlow(C64::new(f64::NAN, 0.0));
    assert!(
        geometry
            .drive_sources_from_material(&state, omega, boundary, &bad)
            .is_err()
    );
    assert!(matches!(
        geometry.input_reflection_from_material(&state, 1e9, boundary),
        Err(FluidDuctError::Duct(DuctError::RadiationKaTooLarge { .. }))
    ));
}

#[test]
fn original_water_states_move_cutoffs_and_drive_their_own_pressure_flow_fields() {
    use fs_qty::{Pressure, Temperature};
    use fs_thermochem::water::if97::{If97Region, If97State};
    let water = [280.0, 350.0].map(|temperature| {
        If97State::at_pressure_temperature(
            Pressure::new(3_000_000.0),
            Temperature::new(temperature),
            If97Region::Region1,
        )
        .unwrap()
    });
    let states = water
        .each_ref()
        .map(|state| state.acoustic_state().unwrap());
    assert!(states[0].sound_speed_m_s() < states[1].sound_speed_m_s());
    assert_ne!(states[0].density_kg_m3(), states[1].density_kg_m3());
    let (radius, length) = (0.03, 0.009);
    let transverse = j1_roots(1).unwrap()[0] / radius;
    let omega =
        f64::midpoint(states[0].sound_speed_m_s(), states[1].sound_speed_m_s()) * transverse;
    let source = duct(&[radius], length);
    let geometry = PreparedModalDuct::new(&source, 2, 1).unwrap();
    for (index, material) in water.iter().enumerate() {
        let state = states[index];
        let waves = geometry
            .input_reflection_from_material(material, omega, MouthBoundary::PressureRelease)
            .unwrap();
        near(
            C64::from_re(waves.cutoffs_at_throat_hz[0]),
            C64::from_re(state.sound_speed_m_s() * transverse / core::f64::consts::TAU),
            1.0,
        );
        if index == 0 {
            assert!(waves.characteristic_impedances[1].re > 0.0);
            assert_eq!(waves.characteristic_impedances[1].im, 0.0);
        } else {
            assert_eq!(waves.characteristic_impedances[1].re, 0.0);
            assert!(waves.characteristic_impedances[1].im < 0.0);
        }
        let area = core::f64::consts::PI * radius * radius;
        let k = omega / state.sound_speed_m_s();
        let (c, s_over_q) = axial_functions(k * k - transverse * transverse, length);
        let load = C64::new(0.0, -omega * state.density_kg_m3() * s_over_q / (area * c));
        let resistance = C64::from_re(0.4 * state.characteristic_impedance_pa_s_m() / area);
        let open = C64::new(1.0, -0.3);
        let expected_flow = open / (resistance + load);
        let result = geometry
            .drive_sources_from_material(
                material,
                omega,
                MouthBoundary::PressureRelease,
                &[
                    PortExcitation::VolumeFlow(C64::ZERO),
                    pressure(open, resistance),
                ],
            )
            .unwrap();
        near(
            result.inlet.inlet_pressure[1],
            load * expected_flow,
            open.abs(),
        );
        near(
            result.inlet.inlet_flow[1],
            expected_flow,
            expected_flow.abs(),
        );
        near(
            result.mouth_flow[1],
            expected_flow.scale(1.0 / c),
            expected_flow.abs(),
        );
        balanced(&result.inlet);
    }
}
