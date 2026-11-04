//! Independent pressure/flow controls for material-driven transverse waves.

use super::*;
use crate::modal::{MAX_MODES, j1_roots, junction_projection, outlet};
use crate::{Duct, Segment, Termination};
use fs_la::eigen_complex::lu_complex;
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use std::cell::Cell;

mod sources;

fn liquid() -> FluidAcousticState {
    // Authored coefficients for analytic controls, not measured water data.
    FluidAcousticState::try_new(1000.0, 1500.0).unwrap()
}

fn duct(radii: &[f64], length: f64) -> Duct {
    Duct {
        segments: radii
            .iter()
            .map(|&radius| Segment::Cylinder { radius, length })
            .collect(),
    }
}

fn near(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * scale.max(expected.abs()).max(1e-30),
        "actual {actual:?}, expected {expected:?}, scale {scale:e}"
    );
}

// Independent real functions for the exact pressure/flow ODE, analytic at q=0.
fn axial_functions(q_squared: f64, length: f64) -> (f64, f64) {
    if q_squared == 0.0 {
        (1.0, length)
    } else if q_squared > 0.0 {
        let q = q_squared.sqrt();
        ((q * length).cos(), (q * length).sin() / q)
    } else {
        let q = (-q_squared).sqrt();
        ((q * length).cosh(), (q * length).sinh() / q)
    }
}

#[test]
fn uniform_liquid_modes_resolve_propagating_and_evanescent_wavelengths() {
    let state = liquid();
    let (radius, length) = (0.03, 0.012);
    let source = duct(&[radius], length);
    let geometry = PreparedModalDuct::new(&source, 4, 1).unwrap();
    let roots = j1_roots(3).unwrap();
    let area = core::f64::consts::PI * radius * radius;
    let reference = state.characteristic_impedance_pa_s_m() / area;
    for k in [20.0, 100.0, 150.0, 300.0, 450.0] {
        let omega = state.sound_speed_m_s() * k;
        for boundary in [MouthBoundary::Rigid, MouthBoundary::PressureRelease] {
            let response = geometry
                .mouth_response_from_material(&state, omega, boundary)
                .unwrap();
            let impedance = geometry
                .input_impedance_from_material(&state, omega, boundary)
                .unwrap();
            assert_eq!(response.input.impedance_matrix, impedance.impedance_matrix);
            for mode in 0..4 {
                let transverse = if mode == 0 {
                    0.0
                } else {
                    roots[mode - 1] / radius
                };
                let q2 = k * k - transverse * transverse;
                let (c, s_over_q) = axial_functions(q2, length);
                let inertia = omega * state.density_kg_m3() / area;
                let (z, flow) = if boundary == MouthBoundary::PressureRelease {
                    (
                        C64::new(0.0, -inertia * s_over_q / c),
                        C64::from_re(1.0 / c),
                    )
                } else {
                    (C64::new(0.0, inertia * c / (q2 * s_over_q)), C64::ZERO)
                };
                for other in 0..4 {
                    near(
                        response.input.impedance_matrix[mode * 4 + other],
                        if mode == other { z } else { C64::ZERO },
                        reference,
                    );
                    near(
                        response.flow_map[mode * 4 + other],
                        if mode == other { flow } else { C64::ZERO },
                        1.0,
                    );
                }
            }
        }
    }
}

// Solve every station's endpoint pressure/flow equations simultaneously. The
// unknowns are P_left and Z0 U_left, never production reflection coordinates.
// Small test lengths keep this independent global oracle well conditioned.
fn direct_pressure_flow(
    geometry: &PreparedModalDuct<'_>,
    state: &FluidAcousticState,
    omega: f64,
    boundary: MouthBoundary,
    load: &[C64],
) -> (Vec<C64>, Vec<C64>) {
    let n = geometry.n_modes;
    let dimension = 2 * n * geometry.stations.len();
    let k = omega / state.sound_speed_m_s();
    let endpoint = |station: usize, right: bool| {
        let span = &geometry.stations[station];
        let area = core::f64::consts::PI * span.radius * span.radius;
        let z0 = state.characteristic_impedance_pa_s_m() / area;
        let inertia = omega * state.density_kg_m3() / area;
        let mut p = vec![C64::ZERO; n * dimension];
        let mut u = p.clone();
        for mode in 0..n {
            let transverse = if mode == 0 {
                0.0
            } else {
                geometry.gammas[mode - 1] / span.radius
            };
            let q2 = k * k - transverse * transverse;
            let (c, s_over_q) = axial_functions(q2, if right { span.length } else { 0.0 });
            let a = 2 * n * station + mode;
            let b = a + n;
            p[mode * dimension + a] = C64::from_re(c);
            p[mode * dimension + b] = C64::new(0.0, inertia * s_over_q / z0);
            u[mode * dimension + a] = C64::new(0.0, q2 * s_over_q / inertia);
            u[mode * dimension + b] = C64::from_re(c / z0);
        }
        (p, u)
    };
    let (pin, uin) = endpoint(0, false);
    let (pout, uout) = endpoint(geometry.stations.len() - 1, true);
    let mut matrix = vec![C64::ZERO; dimension * dimension];
    matrix[..n * dimension].copy_from_slice(&uin);
    let mut row = n;
    for station in 0..geometry.stations.len() - 1 {
        let (p, u) = endpoint(station, true);
        let (q, v) = endpoint(station + 1, false);
        let (a, b) = (
            geometry.stations[station].radius,
            geometry.stations[station + 1].radius,
        );
        let projection = junction_projection(a.min(b), a.max(b), &geometry.gammas, n);
        for mode in 0..n {
            for col in 0..dimension {
                let mut pressure = if a < b {
                    p[mode * dimension + col]
                } else {
                    q[mode * dimension + col]
                };
                let mut flow = if a < b {
                    v[mode * dimension + col]
                } else {
                    u[mode * dimension + col]
                };
                for other in 0..n {
                    pressure = pressure
                        - projection[mode * n + other]
                            * if a < b {
                                q[other * dimension + col]
                            } else {
                                p[other * dimension + col]
                            };
                    flow = flow
                        - projection[other * n + mode]
                            * if a < b {
                                u[other * dimension + col]
                            } else {
                                v[other * dimension + col]
                            };
                }
                matrix[(row + mode) * dimension + col] = pressure;
                matrix[(row + n + mode) * dimension + col] = flow;
            }
        }
        row += 2 * n;
    }
    for mode in 0..n {
        for col in 0..dimension {
            matrix[(row + mode) * dimension + col] = match boundary {
                MouthBoundary::Rigid => uout[mode * dimension + col],
                MouthBoundary::PressureRelease => pout[mode * dimension + col],
                MouthBoundary::Flanged => (0..n).fold(pout[mode * dimension + col], |p, other| {
                    p - load[mode * n + other] * uout[other * dimension + col]
                }),
            };
        }
    }
    let scales: Vec<_> = matrix
        .chunks_exact(dimension)
        .map(|row| row.iter().map(|value| value.abs()).fold(0.0_f64, f64::max))
        .collect();
    for (row, &scale) in matrix.chunks_exact_mut(dimension).zip(&scales) {
        for value in row {
            *value = value.scale(1.0 / scale);
        }
    }
    let lu = lu_complex(&matrix, dimension).unwrap();
    let mut impedance = vec![C64::ZERO; n * n];
    let mut flow_map = impedance.clone();
    for col in 0..n {
        let mut solution = vec![C64::ZERO; dimension];
        solution[col] = C64::from_re(1.0 / scales[col]);
        lu.solve(&mut solution);
        for mode in 0..n {
            for i in 0..dimension {
                impedance[mode * n + col] =
                    impedance[mode * n + col] + pin[mode * dimension + i] * solution[i];
                flow_map[mode * n + col] =
                    flow_map[mode * n + col] + uout[mode * dimension + i] * solution[i];
            }
        }
    }
    (impedance, flow_map)
}

#[test]
fn liquid_cutoffs_and_radius_junctions_match_global_pressure_flow_equations() {
    let state = liquid();
    let radius = j1_roots(1).unwrap()[0] / 128.0;
    let omega = 128.0 * state.sound_speed_m_s();
    for radii in [
        vec![radius],
        vec![radius, 0.045, 0.038],
        vec![0.038, 0.045, radius],
    ] {
        let source = duct(&radii, 0.009);
        let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
        for boundary in [MouthBoundary::PressureRelease, MouthBoundary::Flanged] {
            let response = geometry
                .mouth_response_from_material(&state, omega, boundary)
                .unwrap();
            let load = crate::flanged::flanged_modal_impedance_from_material(
                &state,
                geometry.mouth_radius_m(),
                omega,
                3,
            )
            .unwrap();
            let (expected_z, expected_flow) =
                direct_pressure_flow(&geometry, &state, omega, boundary, &load);
            for (&z, &expected) in response.input.impedance_matrix.iter().zip(&expected_z) {
                near(z, expected, 1e9);
            }
            for (&flow, &expected) in response.flow_map.iter().zip(&expected_flow) {
                near(flow, expected, 1.0);
            }
            if radii.len() == 1 && boundary == MouthBoundary::PressureRelease {
                let inertia = omega * state.density_kg_m3() * 0.009
                    / (core::f64::consts::PI * radius * radius);
                near(
                    response.input.impedance_matrix[4],
                    C64::new(0.0, -inertia),
                    inertia,
                );
                near(response.flow_map[4], C64::ONE, 1.0);
            }
        }
    }
    let source = duct(&[radius], 0.009);
    let geometry = PreparedModalDuct::new(&source, 2, 1).unwrap();
    assert!(matches!(
        geometry.input_impedance_from_material(&state, omega, MouthBoundary::Rigid),
        Err(FluidDuctError::Duct(DuctError::Singular))
    ));
}

#[test]
fn gas_material_paths_retain_native_boundary_and_wall_loss_results() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let source = duct(&[0.014, 0.02, 0.017], 0.011);
    let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
    for omega in [5000.0, 60_000.0] {
        for boundary in [
            MouthBoundary::Rigid,
            MouthBoundary::PressureRelease,
            MouthBoundary::Flanged,
        ] {
            for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
                let native = if boundary == MouthBoundary::Flanged {
                    outlet::mouth_response_flanged(&source, &gas, &gas, omega, loss, 3, 1).unwrap()
                } else {
                    let termination = if boundary == MouthBoundary::Rigid {
                        Termination::Closed
                    } else {
                        Termination::IdealOpen
                    };
                    outlet::mouth_response(&source, &gas, omega, loss, termination, 3, 1).unwrap()
                };
                let shared = geometry
                    .mouth_response_from_material_with_loss(&gas, omega, loss, boundary)
                    .unwrap();
                assert_eq!(shared.input.impedance_matrix, native.input.impedance_matrix);
                assert_eq!(shared.flow_map, native.flow_map);
                assert_eq!(
                    shared.input.cutoffs_at_throat_hz,
                    native.input.cutoffs_at_throat_hz
                );
                if loss == LossModel::Lossless {
                    let acoustic = geometry
                        .mouth_response_from_material(&gas, omega, boundary)
                        .unwrap();
                    assert_eq!(
                        acoustic.input.impedance_matrix,
                        native.input.impedance_matrix
                    );
                    assert_eq!(acoustic.flow_map, native.flow_map);
                }
            }
        }
    }
}

fn apply(matrix: &[C64], flow: &[C64]) -> Vec<C64> {
    matrix
        .chunks_exact(flow.len())
        .map(|row| {
            row.iter()
                .zip(flow)
                .fold(C64::ZERO, |sum, (&z, &u)| sum + z * u)
        })
        .collect()
}

fn power(pressure: &[C64], flow: &[C64]) -> f64 {
    pressure
        .iter()
        .zip(flow)
        .map(|(&p, &u)| 0.5 * (u.conj() * p).re)
        .sum()
}

#[test]
fn sourced_water_drives_mode_conversion_and_complete_radiation_power() {
    use fs_qty::{Pressure, Temperature};
    use fs_thermochem::water::if97::{If97Region, If97State};
    let water = If97State::at_pressure_temperature(
        Pressure::new(3_000_000.0),
        Temperature::new(300.0),
        If97Region::Region1,
    )
    .unwrap();
    let state = water.acoustic_state().unwrap();
    let source = duct(&[0.025, 0.038, 0.028], 0.009);
    let geometry = PreparedModalDuct::new(&source, 4, 1).unwrap();
    let flows = [
        C64::new(1.0, 0.2),
        C64::new(-0.3, 0.7),
        C64::new(0.1, -0.4),
        C64::new(0.2, 0.1),
    ];
    for ka in [1.0, 4.3, 8.1] {
        let omega = ka * state.sound_speed_m_s() / geometry.mouth_radius_m();
        let response = geometry
            .mouth_response_from_material(&water, omega, MouthBoundary::Flanged)
            .unwrap();
        let radiation = crate::flanged::flanged_modal_impedance_from_material(
            &water,
            geometry.mouth_radius_m(),
            omega,
            4,
        )
        .unwrap();
        let mouth_flow = apply(&response.flow_map, &flows);
        let input_power = power(&apply(&response.input.impedance_matrix, &flows), &flows);
        let output_power = power(&apply(&radiation, &mouth_flow), &mouth_flow);
        assert!(output_power > 0.0);
        assert!((input_power / output_power - 1.0).abs() < 2e-8);
        for row in 0..4 {
            for col in 0..4 {
                near(
                    response.input.impedance_matrix[row * 4 + col],
                    response.input.impedance_matrix[col * 4 + row],
                    1e9,
                );
            }
        }
        assert!(
            response.flow_map[4].abs() > 1e-5,
            "the radius changes must activate higher modes"
        );
        near(
            C64::from_re(response.input.cutoffs_at_throat_hz[0]),
            C64::from_re(
                state.sound_speed_m_s() * j1_roots(1).unwrap()[0]
                    / (core::f64::consts::TAU * 0.025),
            ),
            1.0,
        );
    }
}

struct CountedFluid {
    state: ThermoviscousFluidState,
    acoustic_samples: Cell<usize>,
    complete_samples: Cell<usize>,
}
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.acoustic_samples.set(self.acoustic_samples.get() + 1);
        Ok(self.state.acoustic())
    }
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("the modal solver must reuse its single complete snapshot")
    }
}
impl ThermoviscousFluid for CountedFluid {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.complete_samples.set(self.complete_samples.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn one_source_snapshot_drives_every_radius_boundary_and_loss_term() {
    let acoustic = liquid();
    let transport =
        FluidTransportState::try_new(acoustic.density_kg_m3(), 4200.0, 0.6, 0.001).unwrap();
    let material = CountedFluid {
        state: ThermoviscousFluidState::try_new(acoustic, transport, 4100.0).unwrap(),
        acoustic_samples: Cell::new(0),
        complete_samples: Cell::new(0),
    };
    let source = duct(&[0.025, 0.038, 0.028], 0.009);
    let geometry = PreparedModalDuct::new(&source, 3, 1).unwrap();
    let omega = 1.2 * acoustic.sound_speed_m_s() / geometry.mouth_radius_m();
    let incident_flow = [C64::ONE, C64::ZERO, C64::ZERO];
    for boundary in [
        MouthBoundary::Rigid,
        MouthBoundary::PressureRelease,
        MouthBoundary::Flanged,
    ] {
        geometry
            .mouth_response_from_material(&material, omega, boundary)
            .unwrap();
        for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
            let response = geometry
                .mouth_response_from_material_with_loss(&material, omega, loss, boundary)
                .unwrap();
            let impedance = geometry
                .input_impedance_from_material_with_loss(&material.state, omega, loss, boundary)
                .unwrap();
            assert_eq!(response.input.impedance_matrix, impedance.impedance_matrix);
            let incoming = power(
                &apply(&response.input.impedance_matrix, &incident_flow),
                &incident_flow,
            );
            let outgoing = if boundary == MouthBoundary::Flanged {
                let radiation = crate::flanged::flanged_modal_impedance_from_material(
                    &acoustic,
                    geometry.mouth_radius_m(),
                    omega,
                    3,
                )
                .unwrap();
                let flow = apply(&response.flow_map, &incident_flow);
                power(&apply(&radiation, &flow), &flow)
            } else {
                0.0
            };
            if loss == LossModel::Lossless {
                // Ideal boundaries carry zero real power; normalize roundoff
                // by the physical pressure/flow scale, not that zero result.
                let reference_power = 0.5 * acoustic.characteristic_impedance_pa_s_m()
                    / (core::f64::consts::PI * 0.025 * 0.025);
                assert!((incoming - outgoing).abs() < 2e-10 * reference_power);
            } else {
                assert!(
                    incoming > outgoing,
                    "the material's wall diffusion must dissipate power"
                );
            }
        }
    }
    assert_eq!(material.acoustic_samples.get(), 3);
    assert_eq!(material.complete_samples.get(), 9);
    for omega in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            geometry.input_impedance_from_material(&material, omega, MouthBoundary::Flanged),
            Err(FluidDuctError::Duct(DuctError::BadParameter { .. }))
        ));
    }
    assert_eq!(material.acoustic_samples.get(), 3);
    assert!(PreparedModalDuct::new(&source, MAX_MODES + 1, 1).is_err());
    assert!(matches!(
        geometry.input_impedance_from_material(&material, 1e9, MouthBoundary::Flanged),
        Err(FluidDuctError::Duct(DuctError::RadiationKaTooLarge { .. }))
    ));

    struct FailedSource;
    impl FluidAcoustics for FailedSource {
        fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
            Err(MaterialError::Parameters {
                what: "original source refused its phase".into(),
            })
        }
    }
    assert!(
        matches!(geometry.input_impedance_from_material(&FailedSource, 1000.0, MouthBoundary::Rigid),
        Err(FluidDuctError::Material(MaterialError::Parameters { what })) if what == "original source refused its phase")
    );
}
