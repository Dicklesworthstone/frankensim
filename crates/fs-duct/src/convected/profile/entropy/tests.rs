//! G1/G3: physical entropy transport, independent primitive Euler and compact
//! invariants, source/flow reversal, and boundary work. G4: bounded refusal.

use super::*;
use crate::convected::phase;
use crate::convected::profile::IsentropicMeanFlow;
use crate::convected::termination::FlowImpedance;
use crate::profile::{RadiusProfile, RiccatiSteps};
use fs_material::gas::{GasSpec, GasState};

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn close(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(1e-30),
        "{a:?} != {b:?}; scale={scale:e}, tolerance={tolerance:e}"
    );
}

fn physical_incoming(field: EntropyProfileField, side: usize) -> C64 {
    let direction = if side == 0 { 1.0 } else { -1.0 };
    (field.pressure_pa
        + field
            .particle_velocity_m_s
            .scale(direction * field.mean.density_kg_m3 * field.mean.sound_speed_m_s))
    .scale(0.5)
}

#[test]
fn uniform_flow_transports_original_entropy_without_generating_sound() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 0.0, 0.6).unwrap();
    let sigma = C64::new(0.002, -0.001);
    let omega = 2400.0;
    for sign in [-1.0, 1.0] {
        let velocity = sign * 0.35 * gas.sound_speed;
        let flow = IsentropicMeanFlow::new(&gas, &shape, velocity).unwrap();
        let acoustic = flow.prepare(omega, RiccatiSteps(32)).unwrap();
        let transport = acoustic.prepare_entropy().unwrap();
        assert_eq!(
            transport.entropy_inlet(),
            if sign > 0.0 {
                EntropyInlet::Left
            } else {
                EntropyInlet::Right
            }
        );
        let field = transport.prepare_excitation([C64::ZERO; 2], sigma).unwrap();
        assert_eq!(field.outgoing_power_waves(), [C64::ZERO; 2]);
        for index in [0, 7, 16, 32] {
            let actual = field.at_node(index).unwrap();
            let upstream = if sign > 0.0 { 0.0 } else { shape.length() };
            let expected = sigma * phase(omega * (actual.x_m - upstream) / velocity);
            close(actual.relative_entropy, expected, sigma.abs(), 3e-13);
            assert_eq!(actual.pressure_pa, C64::ZERO);
            assert_eq!(actual.particle_velocity_m_s, C64::ZERO);
            assert_eq!(actual.acoustic_power_w, 0.0);
            assert_eq!(actual.entropy_conversion_power_per_length_w_m, 0.0);
            close(
                actual.density_perturbation_kg_m3,
                expected.scale(-gas.density),
                gas.density * sigma.abs(),
                3e-13,
            );
            close(
                actual.temperature_perturbation_k,
                expected.scale(gas.temperature),
                gas.temperature * sigma.abs(),
                3e-13,
            );
            close(
                actual.mass_flow_perturbation_kg_s,
                expected.scale(-flow.mass_flow_kg_s()),
                flow.mass_flow_kg_s().abs() * sigma.abs(),
                3e-13,
            );
            close(
                actual.stagnation_enthalpy_perturbation_j_kg,
                expected.scale(gas.specific_heat_cp * gas.temperature),
                gas.specific_heat_cp * gas.temperature * sigma.abs(),
                3e-13,
            );
        }
        close(
            transport.entropy_transfer(),
            phase(omega * shape.length() / velocity.abs()),
            1.0,
            3e-13,
        );
    }
}

#[test]
fn zero_entropy_reuses_original_acoustics_and_coherent_sources_superpose() {
    let gas = gas();
    let shape =
        RadiusProfile::spline(vec![(0.0, 0.012), (0.2, 0.014), (0.6, 0.017)], None).unwrap();
    let flow = IsentropicMeanFlow::new(&gas, &shape, 0.4 * gas.sound_speed).unwrap();
    let acoustic = flow.prepare(1800.0, RiccatiSteps(32)).unwrap();
    let transport = acoustic.prepare_entropy().unwrap();
    assert!(core::ptr::eq(transport.acoustic_profile(), &acoustic));
    assert_eq!(transport.node_count(), acoustic.node_count());
    assert_eq!(acoustic.node_position_m(32), Some(0.2));
    let incoming = [C64::new(2.0, -0.3), C64::new(-0.4, 0.7)];
    let original = acoustic.prepare_excitation(incoming).unwrap();
    let a = transport.prepare_excitation(incoming, C64::ZERO).unwrap();
    let e = transport
        .prepare_excitation([C64::ZERO; 2], C64::new(0.001, 0.0003))
        .unwrap();
    let combined = transport
        .prepare_excitation(incoming, e.upstream_relative_entropy())
        .unwrap();
    for index in [0, 13, 32, 49, 64] {
        let old = original.at_node(index).unwrap();
        let zero = a.at_node(index).unwrap();
        assert_eq!(zero.pressure_pa, old.pressure_pa);
        assert_eq!(zero.particle_velocity_m_s, old.particle_velocity_m_s);
        assert_eq!(
            zero.mass_flow_perturbation_kg_s,
            old.mass_flow_perturbation_kg_s
        );
        assert_eq!(
            zero.stagnation_enthalpy_perturbation_j_kg,
            old.stagnation_enthalpy_perturbation_j_kg
        );
        let entropy = e.at_node(index).unwrap();
        let total = combined.at_node(index).unwrap();
        for (actual, expected) in [
            (total.pressure_pa, zero.pressure_pa + entropy.pressure_pa),
            (
                total.particle_velocity_m_s,
                zero.particle_velocity_m_s + entropy.particle_velocity_m_s,
            ),
            (
                total.mass_flow_perturbation_kg_s,
                zero.mass_flow_perturbation_kg_s + entropy.mass_flow_perturbation_kg_s,
            ),
            (
                total.stagnation_enthalpy_perturbation_j_kg,
                zero.stagnation_enthalpy_perturbation_j_kg
                    + entropy.stagnation_enthalpy_perturbation_j_kg,
            ),
        ] {
            close(actual, expected, expected.abs(), 3e-13);
        }
    }
    for (side, index) in [0, transport.node_count() - 1].into_iter().enumerate() {
        close(
            physical_incoming(combined.at_node(index).unwrap(), side),
            incoming[side],
            3.0,
            1e-12,
        );
    }
}

// Independent dimensional primitive-variable LEE. This integrates pressure,
// velocity and entropy; it never calls the production B/j generator, block
// exponential or source response. Original geometry derivatives determine the
// steady gradients through mass conservation and steady momentum.
fn primitive_rhs(flow: &IsentropicMeanFlow<'_>, omega: f64, x: f64, y: [C64; 3]) -> [C64; 3] {
    let s = flow.at_position(x).unwrap();
    let geometry = flow.radius_profile();
    let alpha = 2.0 * geometry.slope(x).unwrap() / geometry.radius(x).unwrap();
    let ux = -s.velocity_m_s * alpha / (1.0 - s.mach * s.mach);
    let rho_log_gradient = s.mach * s.mach * alpha / (1.0 - s.mach * s.mach);
    let c2 = s.sound_speed_m_s * s.sound_speed_m_s;
    let momentum = C64::new(0.0, omega * s.density_kg_m3) * y[1] - y[1].scale(s.density_kg_m3 * ux)
        + y[0].scale(rho_log_gradient)
        + y[2].scale(s.density_kg_m3 * s.velocity_m_s * ux);
    let pressure_transport = C64::new(0.0, omega) * y[0]
        - y[1].scale(s.density_kg_m3 * c2 * (alpha + rho_log_gradient))
        + y[0].scale(flow.source().gamma * s.velocity_m_s * rho_log_gradient);
    let velocity_derivative = (pressure_transport - momentum.scale(s.velocity_m_s))
        .scale(1.0 / (s.density_kg_m3 * (c2 - s.velocity_m_s * s.velocity_m_s)));
    [
        momentum - velocity_derivative.scale(s.density_kg_m3 * s.velocity_m_s),
        velocity_derivative,
        C64::new(0.0, omega / s.velocity_m_s) * y[2],
    ]
}

fn primitive_end(flow: &IsentropicMeanFlow<'_>, omega: f64, mut y: [C64; 3], n: usize) -> [C64; 3] {
    let length = flow.radius_profile().length();
    let h = length / n as f64;
    for i in 0..n {
        let x = i as f64 * h;
        let k1 = primitive_rhs(flow, omega, x, y);
        let k2 = primitive_rhs(
            flow,
            omega,
            x + h * 0.5,
            core::array::from_fn(|j| y[j] + k1[j].scale(h * 0.5)),
        );
        let k3 = primitive_rhs(
            flow,
            omega,
            x + h * 0.5,
            core::array::from_fn(|j| y[j] + k2[j].scale(h * 0.5)),
        );
        let end = if i + 1 == n { length } else { x + h };
        let k4 = primitive_rhs(
            flow,
            omega,
            end,
            core::array::from_fn(|j| y[j] + k3[j].scale(h)),
        );
        y = core::array::from_fn(|j| {
            y[j] + (k1[j] + k2[j].scale(2.0) + k3[j].scale(2.0) + k4[j]).scale(h / 6.0)
        });
    }
    y
}

#[test]
fn full_entropy_acoustic_fields_refine_against_independent_primitive_euler() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, -1.2, 0.4).unwrap();
    let omega = 2100.0;
    let incoming = [C64::new(2.0, -0.3), C64::new(-0.4, 0.7)];
    let sigma = C64::new(0.001, 0.0003);
    for sign in [-1.0, 1.0] {
        let flow = IsentropicMeanFlow::new(&gas, &shape, sign * 0.25 * gas.sound_speed).unwrap();
        let mut errors = Vec::new();
        for panels in [8, 16, 32, 64] {
            let acoustic = flow.prepare(omega, RiccatiSteps(panels)).unwrap();
            let prepared = acoustic.prepare_entropy().unwrap();
            let field = prepared.prepare_excitation(incoming, sigma).unwrap();
            let left = field.at_node(0).unwrap();
            let right = field.at_node(panels).unwrap();
            let initial = [
                left.pressure_pa,
                left.particle_velocity_m_s,
                left.relative_entropy,
            ];
            let reference = primitive_end(&flow, omega, initial, 4096);
            let finer = primitive_end(&flow, omega, initial, 8192);
            let scales = [gas.pressure * gas.gamma, gas.sound_speed, 1.0];
            let actual = [
                right.pressure_pa,
                right.particle_velocity_m_s,
                right.relative_entropy,
            ];
            let mut error: f64 = 0.0;
            for j in 0..3 {
                close(reference[j], finer[j], scales[j] * sigma.abs(), 3e-9);
                error = error.max((actual[j] - finer[j]).abs() / (scales[j] * sigma.abs()));
            }
            errors.push(error);
            close(physical_incoming(left, 0), incoming[0], 3.0, 2e-12);
            close(physical_incoming(right, 1), incoming[1], 3.0, 2e-12);
            close(
                if sign > 0.0 {
                    left.relative_entropy
                } else {
                    right.relative_entropy
                },
                sigma,
                sigma.abs(),
                3e-14,
            );
        }
        assert!(
            errors[0] > 1e-7,
            "nonvacuous entropy conversion ladder: {errors:?}"
        );
        for pair in errors.windows(2) {
            assert!(
                pair[1] < 0.13 * pair[0],
                "fourth-order entropy/LEE refinement: {errors:?}"
            );
        }
        assert!(errors[3] < 2e-7, "resolved full LEE field: {errors:?}");
    }
}

#[test]
fn compact_entropy_noise_matches_constant_total_mass_and_stagnation_enthalpy() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, -1.2, 0.4).unwrap();
    let sigma = 0.001;
    for sign in [-1.0, 1.0] {
        let flow = IsentropicMeanFlow::new(&gas, &shape, sign * 0.25 * gas.sound_speed).unwrap();
        let acoustic = flow.prepare(1e-4, RiccatiSteps(64)).unwrap();
        let prepared = acoustic.prepare_entropy().unwrap();
        let field = prepared
            .prepare_excitation([C64::ZERO; 2], C64::from_re(sigma))
            .unwrap();
        let left = field.at_node(0).unwrap();
        let right = field.at_node(64).unwrap();
        let l = left.mean;
        let r = right.mean;
        // Only outgoing sound is present. At zero frequency total m', H' and
        // sigma are invariant. This independent 2x2 algebra uses only end states.
        let a = (1.0 - l.mach) / l.density_kg_m3;
        let b = (1.0 + r.mach) / r.density_kg_m3;
        let ml = -l.area_m2 / l.sound_speed_m_s * (1.0 - l.mach);
        let mr = r.area_m2 / r.sound_speed_m_s * (1.0 + r.mach);
        let thermal = gas.specific_heat_cp * (r.temperature_k - l.temperature_k) * sigma;
        let expected_left = thermal / (a - b * ml / mr);
        let expected_right = expected_left * ml / mr;
        assert!(expected_left.abs() > 1.0 && expected_right.abs() > 1.0);
        close(
            left.pressure_pa,
            C64::from_re(expected_left),
            expected_left.abs(),
            2e-6,
        );
        close(
            right.pressure_pa,
            C64::from_re(expected_right),
            expected_right.abs(),
            2e-6,
        );
        close(
            left.mass_flow_perturbation_kg_s,
            right.mass_flow_perturbation_kg_s,
            left.mass_flow_perturbation_kg_s.abs(),
            2e-6,
        );
        close(
            left.stagnation_enthalpy_perturbation_j_kg,
            right.stagnation_enthalpy_perturbation_j_kg,
            left.stagnation_enthalpy_perturbation_j_kg.abs(),
            2e-6,
        );
    }
}

#[test]
fn mirrored_geometry_and_flow_reverse_entropy_inlet_and_vector_fields() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 1.2, 0.4).unwrap();
    let flow = IsentropicMeanFlow::new(&gas, &shape, 0.4 * gas.sound_speed).unwrap();
    let end = flow.at_position(shape.length()).unwrap();
    let reverse_gas = GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        end.temperature_k,
        end.pressure_pa,
    )
    .unwrap();
    let reverse_shape =
        RadiusProfile::exponential(shape.outlet_radius(), -1.2, shape.length()).unwrap();
    let reverse = IsentropicMeanFlow::new(&reverse_gas, &reverse_shape, -end.velocity_m_s).unwrap();
    let a = flow.prepare(1800.0, RiccatiSteps(64)).unwrap();
    let b = reverse.prepare(1800.0, RiccatiSteps(64)).unwrap();
    let a = a.prepare_entropy().unwrap();
    let b = b.prepare_entropy().unwrap();
    let incoming = [C64::new(2.0, -0.3), C64::new(-0.4, 0.7)];
    let sigma = C64::new(0.001, 0.0003);
    let a = a.prepare_excitation(incoming, sigma).unwrap();
    let b = b
        .prepare_excitation([incoming[1], incoming[0]], sigma)
        .unwrap();
    for index in [0, 13, 32, 49, 64] {
        let left = a.at_node(index).unwrap();
        let right = b.at_node(64 - index).unwrap();
        for (actual, expected) in [
            (left.pressure_pa, right.pressure_pa),
            (left.particle_velocity_m_s, -right.particle_velocity_m_s),
            (left.relative_entropy, right.relative_entropy),
            (
                left.mass_flow_perturbation_kg_s,
                -right.mass_flow_perturbation_kg_s,
            ),
            (
                left.stagnation_enthalpy_perturbation_j_kg,
                right.stagnation_enthalpy_perturbation_j_kg,
            ),
        ] {
            close(actual, expected, actual.abs().max(expected.abs()), 3e-9);
        }
    }
}

#[test]
fn acoustic_generation_matches_signed_local_entropy_conversion() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 1.2, 0.4).unwrap();
    for sign in [-1.0, 1.0] {
        let flow = IsentropicMeanFlow::new(&gas, &shape, sign * 0.4 * gas.sound_speed).unwrap();
        let mut errors = Vec::new();
        for n in [32, 64, 128] {
            let acoustic = flow.prepare(1800.0, RiccatiSteps(n)).unwrap();
            let transport = acoustic.prepare_entropy().unwrap();
            let field = transport
                .prepare_excitation([C64::ZERO; 2], C64::from_re(0.001))
                .unwrap();
            let mut integral = 0.0;
            for i in 0..=n {
                let weight = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                integral += weight
                    * field
                        .at_node(i)
                        .unwrap()
                        .entropy_conversion_power_per_length_w_m;
            }
            integral *= shape.length() / (3.0 * n as f64);
            let generated = field.net_acoustic_generation_w();
            assert!(generated > 1e-8);
            let boundary = field.at_node(n).unwrap().acoustic_power_w
                - field.at_node(0).unwrap().acoustic_power_w;
            assert!((boundary - generated).abs() < 2e-12 * generated);
            errors.push((integral - generated).abs() / generated);
        }
        for pair in errors.windows(2) {
            assert!(
                pair[1] < 0.15 * pair[0],
                "conversion-work refinement: {errors:?}"
            );
        }
        assert!(errors[2] < 1e-6, "conversion work: {errors:?}");
    }
}

#[test]
fn loaded_entropy_boundaries_use_total_mass_and_enthalpy_with_coherent_drives() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 1.2, 0.4).unwrap();
    let drives = [C64::new(0.8, -0.2), C64::new(-0.3, 0.4)];
    let sigma = C64::new(0.001, 0.0003);
    for sign in [-1.0, 1.0] {
        let flow = IsentropicMeanFlow::new(&gas, &shape, sign * 0.4 * gas.sound_speed).unwrap();
        let acoustic = flow.prepare(1800.0, RiccatiSteps(64)).unwrap();
        let r = acoustic.reference_impedances_j_s_kg2();
        let z = [
            FlowImpedance::try_new(C64::new(0.7, 0.3).scale(r[0])).unwrap(),
            FlowImpedance::try_new(C64::new(1.3, -0.4).scale(r[1])).unwrap(),
        ];
        let transport = acoustic.prepare_entropy().unwrap();
        let loads = transport.prepare_terminations(z, 1e8).unwrap();
        let actual = loads.solve(drives, sigma).unwrap();
        let only_entropy = loads.solve([C64::ZERO; 2], sigma).unwrap();
        let no_entropy = loads.solve(drives, C64::ZERO).unwrap();
        let original_loads = acoustic.prepare_terminations(z, 1e8).unwrap();
        let original = original_loads.solve(drives).unwrap();
        for index in [0, 23, 64] {
            let old = original.field().at_node(index).unwrap();
            let zero = no_entropy.field().at_node(index).unwrap();
            close(
                zero.pressure_pa,
                old.pressure_pa,
                old.pressure_pa.abs(),
                3e-13,
            );
            let total = actual.field().at_node(index).unwrap();
            let entropy = only_entropy.field().at_node(index).unwrap();
            close(
                total.pressure_pa,
                zero.pressure_pa + entropy.pressure_pa,
                total.pressure_pa.abs(),
                3e-13,
            );
        }
        for side in 0..2 {
            let end = actual.terminals()[side];
            let field = actual
                .field()
                .at_node(if side == 0 { 0 } else { 64 })
                .unwrap();
            let mass = (field.particle_velocity_m_s.scale(field.mean.density_kg_m3)
                + field
                    .density_perturbation_kg_m3
                    .scale(field.mean.velocity_m_s))
            .scale(field.mean.area_m2 * if side == 0 { -1.0 } else { 1.0 });
            let enthalpy = field.temperature_perturbation_k.scale(gas.specific_heat_cp)
                + field.particle_velocity_m_s.scale(field.mean.velocity_m_s);
            close(
                end.outward_mass_flow_perturbation_kg_s,
                mass,
                mass.abs(),
                3e-13,
            );
            close(
                end.stagnation_enthalpy_perturbation_j_kg,
                enthalpy,
                enthalpy.abs(),
                3e-13,
            );
            close(
                enthalpy - z[side].value_j_s_kg2() * mass,
                drives[side],
                1000.0 * sigma.abs(),
                2e-10,
            );
            close(
                end.boundary_residual_j_kg,
                C64::ZERO,
                1000.0 * sigma.abs(),
                2e-10,
            );
            let balance = end.enthalpy_mass_pair_power_w - end.resistance_pair_power_w
                + end.source_pair_work_into_duct_w;
            assert!(balance.abs() < 2e-10 * end.resistance_pair_power_w.max(1e-10));
        }
    }
}

#[test]
fn advected_entropy_refuses_zero_flow_unresolved_work_and_invalid_sources() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 0.0, 0.4).unwrap();
    for velocity in [0.0, 1e-100] {
        let flow = IsentropicMeanFlow::new(&gas, &shape, velocity).unwrap();
        let acoustic = flow.prepare(1800.0, RiccatiSteps(1)).unwrap();
        assert!(acoustic.prepare_entropy().is_err());
    }
    let flow = IsentropicMeanFlow::new(&gas, &shape, 0.3 * gas.sound_speed).unwrap();
    let acoustic = flow.prepare(1800.0, RiccatiSteps(16)).unwrap();
    let mut polls = 0;
    assert!(matches!(
        acoustic.prepare_entropy_controlled(&mut || {
            polls += 1;
            polls == 4
        }),
        Err(DuctError::Cancelled { .. })
    ));
    let transport = acoustic.prepare_entropy().unwrap();
    assert!(
        transport
            .prepare_excitation([C64::ZERO; 2], C64::new(f64::NAN, 0.0))
            .is_err()
    );
    assert!(
        transport
            .excite_power_waves([C64::new(f64::INFINITY, 0.0), C64::ZERO], C64::ZERO)
            .is_err()
    );
    let field = transport
        .prepare_excitation([C64::ZERO; 2], C64::ONE)
        .unwrap();
    assert!(field.at_node(transport.node_count()).is_err());
    let z = [FlowImpedance::try_new(C64::ZERO).unwrap(); 2];
    assert!(transport.prepare_terminations(z, 0.0).is_err());
    let loads = transport.prepare_terminations(z, 1e8).unwrap();
    assert!(
        loads
            .solve([C64::ZERO; 2], C64::new(0.0, f64::NAN))
            .is_err()
    );
}
