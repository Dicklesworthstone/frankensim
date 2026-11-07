use super::*;
use crate::convected::UniformFlowCylinder;
use crate::convected::termination::FlowImpedance;
use crate::profile::{RadiusProfile, profile_matrix};
use crate::{LossModel, Segment};
use fs_material::gas::{GasSpec, GasState};

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn omega() -> f64 {
    core::f64::consts::TAU * 723.0
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(1e-30),
        "{a:e} != {b:e}, scale {scale:e}"
    );
}
fn complex_close(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale.max(1e-30),
        "{a:?} != {b:?}, scale {scale:e}"
    );
}
fn loads(references: [f64; 2]) -> [FlowImpedance; 2] {
    [
        FlowImpedance::try_new(C64::new(0.7, 0.3).scale(references[0])).unwrap(),
        FlowImpedance::try_new(C64::new(1.3, -0.4).scale(references[1])).unwrap(),
    ]
}

#[test]
fn uniform_profile_recovers_exact_both_direction_fields_and_loaded_cylinder() {
    let gas = gas();
    let radius = 0.012;
    let length = 0.83;
    let shape = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    let geometry = Segment::Cylinder { radius, length };
    let incoming = [C64::new(1.3, -0.2), C64::new(-0.4, 0.7)];
    for mach in [-0.55, 0.0, 0.4, 0.75] {
        let velocity = mach * gas.sound_speed;
        let flow = IsentropicMeanFlow::new(&gas, &shape, velocity).unwrap();
        let profile = flow.prepare(omega(), RiccatiSteps(64)).unwrap();
        let cylinder = UniformFlowCylinder::try_new(&gas, &geometry, omega(), velocity).unwrap();
        let actual = profile.prepare_excitation(incoming).unwrap();
        let expected = cylinder.prepare_excitation(incoming).unwrap();
        let s = profile.scattering();
        complex_close(s[0], C64::ZERO, 1.0, 3e-12);
        complex_close(s[3], C64::ZERO, 1.0, 3e-12);
        complex_close(
            s[1],
            phase(omega() * length / (gas.sound_speed - velocity)),
            1.0,
            3e-12,
        );
        complex_close(
            s[2],
            phase(omega() * length / (gas.sound_speed + velocity)),
            1.0,
            3e-12,
        );
        for index in [0, 11, 32, 64] {
            let a = actual.at_node(index).unwrap();
            let b = expected.at_position(a.x_m).unwrap();
            for (a, b) in [
                (a.pressure_pa, b.pressure_pa),
                (a.particle_velocity_m_s, b.particle_velocity_m_s),
                (a.mass_flow_perturbation_kg_s, b.mass_flow_perturbation_kg_s),
                (
                    a.stagnation_enthalpy_perturbation_j_kg,
                    b.stagnation_enthalpy_perturbation_j_kg,
                ),
            ] {
                complex_close(a, b, a.abs().max(b.abs()), 3e-11);
            }
            close(
                a.disturbance_energy_density_j_m3,
                b.disturbance_energy_density_j_m3,
                b.disturbance_energy_density_j_m3,
                3e-11,
            );
            close(
                a.disturbance_power_w,
                b.disturbance_power_w,
                expected.boundary_power().incident_disturbance_power_w,
                3e-11,
            );
        }
        let z = loads(profile.reference_impedances_j_s_kg2());
        let profile_loads = profile.prepare_terminations(z, 1e8).unwrap();
        let cylinder_loads = cylinder.prepare_terminations(z, 1e8).unwrap();
        let drives = [C64::new(1.0, 0.2), C64::new(-0.3, 0.5)];
        let a = profile_loads.solve(drives).unwrap();
        let b = cylinder_loads.solve(drives).unwrap();
        for index in [0, 32, 64] {
            let p = a.field().at_node(index).unwrap();
            let q = b.at_position(p.x_m).unwrap();
            complex_close(p.pressure_pa, q.pressure_pa, q.pressure_pa.abs(), 3e-11);
        }
        close(
            a.power().absorbed_disturbance_power_w,
            b.power().absorbed_disturbance_power_w,
            b.power().absorbed_disturbance_power_w,
            3e-11,
        );
    }
}

#[test]
fn original_mean_state_conserves_mass_entropy_and_stagnation_enthalpy() {
    let shape = RadiusProfile::exponential(0.012, 1.1, 0.8).unwrap();
    for gamma in [1.0 + 1e-10, 1.4, 5.0 / 3.0] {
        let mut spec = GasSpec::dry_air_ussa1976();
        spec.gamma = gamma;
        let gas = GasState::try_new(&spec, 293.15, 101_325.0).unwrap();
        for sign in [-1.0, 1.0] {
            let u = sign * 0.45 * gas.sound_speed;
            let flow = IsentropicMeanFlow::new(&gas, &shape, u).unwrap();
            assert!(core::ptr::eq(flow.source(), &gas));
            assert!(core::ptr::eq(flow.radius_profile(), &shape));
            let h0 = gas.specific_heat_cp * gas.temperature + 0.5 * u * u;
            for i in 0..=20 {
                let s = flow
                    .at_position(shape.length() * f64::from(i) / 20.0)
                    .unwrap();
                close(
                    s.density_kg_m3 * s.area_m2 * s.velocity_m_s,
                    flow.mass_flow_kg_s(),
                    flow.mass_flow_kg_s().abs(),
                    5e-13,
                );
                close(
                    gas.specific_heat_cp * s.temperature_k + 0.5 * s.velocity_m_s * s.velocity_m_s,
                    h0,
                    h0,
                    5e-13,
                );
                close(
                    det::ln(s.pressure_pa / gas.pressure)
                        - gamma * det::ln(s.density_kg_m3 / gas.density),
                    0.0,
                    1.0,
                    5e-13,
                );
                assert_eq!(s.mach.is_sign_positive(), sign > 0.0);
            }
        }
    }
}

#[test]
fn no_flow_recovers_original_variable_area_webster_propagation() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 1.2, 0.8).unwrap();
    let flow = IsentropicMeanFlow::new(&gas, &shape, 0.0).unwrap();
    let prepared = flow.prepare(omega(), RiccatiSteps(64)).unwrap();
    let field = prepared
        .prepare_excitation([C64::new(1.2, -0.3), C64::new(-0.2, 0.5)])
        .unwrap();
    let left = field.at_node(0).unwrap();
    let right = field.at_node(64).unwrap();
    let t = profile_matrix(&shape, &gas, omega(), LossModel::Lossless, RiccatiSteps(64))
        .unwrap()
        .matrix;
    complex_close(
        left.pressure_pa,
        t[0] * right.pressure_pa + t[1] * right.volume_flow_m3_s,
        left.pressure_pa.abs(),
        2e-11,
    );
    complex_close(
        left.volume_flow_m3_s,
        t[2] * right.pressure_pa + t[3] * right.volume_flow_m3_s,
        left.volume_flow_m3_s.abs(),
        2e-11,
    );
    close(
        left.pressure_work_power_w,
        left.disturbance_power_w,
        field.incident_disturbance_power_w(),
        2e-12,
    );
}

// Independent primitive-variable linearized Euler oracle. It uses the actual
// geometry derivative and background gradients, NOT the H/m Magnus generator.
fn primitive_rhs(flow: &IsentropicMeanFlow<'_>, omega: f64, x: f64, y: [C64; 2]) -> [C64; 2] {
    let s = flow.at_position(x).unwrap();
    let radius = flow.radius_profile().radius(x).unwrap();
    let alpha = 2.0 * flow.radius_profile().slope(x).unwrap() / radius;
    let velocity_gradient = -s.velocity_m_s * alpha / (1.0 - s.mach * s.mach);
    let density_log_gradient = s.mach * s.mach * alpha / (1.0 - s.mach * s.mach);
    let c2 = s.sound_speed_m_s * s.sound_speed_m_s;
    let momentum = C64::new(0.0, omega * s.density_kg_m3) * y[1]
        - y[1].scale(s.density_kg_m3 * velocity_gradient)
        + y[0].scale(density_log_gradient);
    let continuity = C64::new(0.0, omega) * y[0]
        - y[1].scale(s.density_kg_m3 * c2 * (alpha + density_log_gradient))
        + y[0].scale(flow.source().gamma * s.velocity_m_s * density_log_gradient);
    let ux = (continuity - momentum.scale(s.velocity_m_s))
        .scale(1.0 / (s.density_kg_m3 * (c2 - s.velocity_m_s * s.velocity_m_s)));
    [momentum - ux.scale(s.density_kg_m3 * s.velocity_m_s), ux]
}
fn integrate_primitive(flow: &IsentropicMeanFlow<'_>, omega: f64, mut y: [C64; 2]) -> [C64; 2] {
    let n = 4096;
    let h = flow.radius_profile().length() / n as f64;
    for i in 0..n {
        let x = i as f64 * h;
        let k1 = primitive_rhs(flow, omega, x, y);
        let k2 = primitive_rhs(
            flow,
            omega,
            x + 0.5 * h,
            core::array::from_fn(|j| y[j] + k1[j].scale(0.5 * h)),
        );
        let k3 = primitive_rhs(
            flow,
            omega,
            x + 0.5 * h,
            core::array::from_fn(|j| y[j] + k2[j].scale(0.5 * h)),
        );
        let end = if i + 1 == n {
            flow.radius_profile().length()
        } else {
            x + h
        };
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
fn nonuniform_flow_converges_against_independent_primitive_euler() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 1.4, 0.8).unwrap();
    for sign in [-1.0, 1.0] {
        let flow = IsentropicMeanFlow::new(&gas, &shape, sign * 0.5 * gas.sound_speed).unwrap();
        let mut errors = Vec::new();
        for steps in [16, 32, 64, 128] {
            let prepared = flow.prepare(omega(), RiccatiSteps(steps)).unwrap();
            let field = prepared
                .prepare_excitation([C64::new(1.1, -0.2), C64::new(-0.4, 0.3)])
                .unwrap();
            let left = field.at_node(0).unwrap();
            let right = field.at_node(steps).unwrap();
            let oracle = integrate_primitive(
                &flow,
                omega(),
                [left.pressure_pa, left.particle_velocity_m_s],
            );
            let error = (right.pressure_pa - oracle[0]).abs()
                + (right.particle_velocity_m_s - oracle[1]).abs() * gas.density * gas.sound_speed;
            errors.push(error);
            for node in [0, steps / 2, steps] {
                let s = field.at_node(node).unwrap();
                close(
                    s.disturbance_power_w,
                    left.disturbance_power_w,
                    field.incident_disturbance_power_w(),
                    3e-11,
                );
            }
        }
        for pair in errors.windows(2) {
            assert!(pair[1] < 0.12 * pair[0], "fourth-order ladder {errors:?}");
        }
        assert!(errors[3] < 3e-7, "primitive Euler residual {errors:?}");
    }
}

#[test]
fn variable_flow_loads_obey_physical_laws_and_coherent_source_work() {
    let gas = gas();
    let shape =
        RadiusProfile::spline(vec![(0.0, 0.012), (0.4, 0.0135), (0.8, 0.017)], None).unwrap();
    let flow = IsentropicMeanFlow::new(&gas, &shape, 0.35 * gas.sound_speed).unwrap();
    let profile = flow.prepare(omega(), RiccatiSteps(32)).unwrap();
    assert_eq!(profile.steps_used(), 64);
    assert_eq!(profile.node_position_m(32), Some(0.4));
    let z = loads(profile.reference_impedances_j_s_kg2());
    let prepared = profile.prepare_terminations(z, 1e8).unwrap();
    let drives = [C64::new(1.2, -0.3), C64::new(-0.2, 0.8)];
    let actual = prepared.solve(drives).unwrap();
    let only_left = prepared.solve([drives[0], C64::ZERO]).unwrap();
    let only_right = prepared.solve([C64::ZERO, drives[1]]).unwrap();
    let scale = actual.power().absorbed_disturbance_power_w;
    assert!(scale > 0.0);
    close(actual.power().balance_residual_w, 0.0, scale, 3e-11);
    close(actual.power().interior_flux_residual_w, 0.0, scale, 3e-11);
    for i in 0..2 {
        let end = actual.terminals()[i];
        complex_close(
            end.enthalpy_j_kg,
            z[i].value_j_s_kg2() * end.outward_mass_flow_kg_s + drives[i],
            drives[i].abs(),
            3e-11,
        );
        let local = actual.field().at_node(if i == 0 { 0 } else { 64 }).unwrap();
        let rho = local.mean.density_kg_m3;
        let u = local.mean.velocity_m_s;
        let c = local.mean.sound_speed_m_s;
        let physical_h = local.pressure_pa.scale(1.0 / rho) + local.particle_velocity_m_s.scale(u);
        let physical_m = (local.particle_velocity_m_s.scale(rho)
            + local.pressure_pa.scale(u / (c * c)))
        .scale(local.mean.area_m2);
        close(
            local.disturbance_power_w,
            0.5 * (physical_h * physical_m.conj()).re,
            scale,
            3e-12,
        );
    }
    for index in [0, 13, 32, 64] {
        let p = actual.field().at_node(index).unwrap().pressure_pa;
        let sum = only_left.field().at_node(index).unwrap().pressure_pa
            + only_right.field().at_node(index).unwrap().pressure_pa;
        complex_close(p, sum, p.abs(), 3e-12);
    }
    assert!(
        (actual.power().absorbed_disturbance_power_w
            - only_left.power().absorbed_disturbance_power_w
            - only_right.power().absorbed_disturbance_power_w)
            .abs()
            > 1e-4 * scale
    );
    assert_eq!(
        prepared
            .solve([C64::ZERO; 2])
            .unwrap()
            .power()
            .absorbed_disturbance_power_w,
        0.0
    );
}

#[test]
fn mirrored_geometry_and_reversed_flow_exchange_the_physical_ends() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 1.1, 0.8).unwrap();
    let flow = IsentropicMeanFlow::new(&gas, &shape, 0.45 * gas.sound_speed).unwrap();
    let end = flow.at_position(shape.length()).unwrap();
    let reversed_gas = GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        end.temperature_k,
        end.pressure_pa,
    )
    .unwrap();
    let reversed_shape =
        RadiusProfile::exponential(shape.outlet_radius(), -1.1, shape.length()).unwrap();
    let reverse =
        IsentropicMeanFlow::new(&reversed_gas, &reversed_shape, -end.velocity_m_s).unwrap();
    let forward = flow.prepare(omega(), RiccatiSteps(64)).unwrap();
    let backward = reverse.prepare(omega(), RiccatiSteps(64)).unwrap();
    let a = forward.scattering();
    let b = backward.scattering();
    for (x, y) in [(a[0], b[3]), (a[1], b[2]), (a[2], b[1]), (a[3], b[0])] {
        complex_close(x, y, 1.0, 3e-11);
    }
    let a = forward
        .prepare_excitation([C64::new(1.2, -0.3), C64::new(-0.4, 0.2)])
        .unwrap();
    let b = backward
        .prepare_excitation([C64::new(-0.4, 0.2), C64::new(1.2, -0.3)])
        .unwrap();
    for i in [0, 13, 32, 64] {
        let x = a.at_node(i).unwrap();
        let y = b.at_node(64 - i).unwrap();
        complex_close(x.pressure_pa, y.pressure_pa, x.pressure_pa.abs(), 3e-11);
        complex_close(
            x.particle_velocity_m_s,
            -y.particle_velocity_m_s,
            x.particle_velocity_m_s.abs(),
            3e-11,
        );
        close(
            x.disturbance_power_w,
            -y.disturbance_power_w,
            a.incident_disturbance_power_w(),
            3e-11,
        );
    }
}

#[test]
fn whole_profile_admission_sees_hidden_throats_and_bulges_and_bounds_work() {
    let gas = gas();
    let hidden_throat =
        RadiusProfile::spline(vec![(0.0, 0.012), (1.0, 0.012)], Some([-0.03, 0.03])).unwrap();
    assert!(IsentropicMeanFlow::new(&gas, &hidden_throat, 0.4 * gas.sound_speed).is_err());
    let bulge =
        RadiusProfile::spline(vec![(0.0, 0.012), (1.0, 0.012)], Some([0.06, -0.06])).unwrap();
    let flow = IsentropicMeanFlow::new(&gas, &bulge, 0.25 * gas.sound_speed).unwrap();
    let bound = flow.transverse_cutoff_bound_rad_s();
    assert!(flow.prepare(bound, RiccatiSteps(1)).is_err());
    assert!(bound < super::super::FIRST_TRANSVERSE_ROOT * gas.sound_speed / 0.012);
    for velocity in [gas.sound_speed, -gas.sound_speed, f64::NAN] {
        assert!(IsentropicMeanFlow::new(&gas, &bulge, velocity).is_err());
    }
    let near_sonic = gas.sound_speed * det::sqrt(1.0 - 0.5 * MINIMUM_SUBSONIC_MARGIN);
    assert!(IsentropicMeanFlow::new(&gas, &bulge, near_sonic).is_err());
    for steps in [0, MAX_RICCATI_STEPS + 1] {
        assert!(flow.prepare(omega(), RiccatiSteps(steps)).is_err());
    }
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(flow.prepare(frequency, RiccatiSteps(2)).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        flow.prepare_controlled(omega(), RiccatiSteps(20), &mut || {
            polls += 1;
            polls == 4
        }),
        Err(DuctError::Cancelled { completed_steps: 3 })
    ));
    let profile = flow.prepare(omega(), RiccatiSteps(2)).unwrap();
    let silence = profile.prepare_excitation([C64::ZERO; 2]).unwrap();
    assert!(silence.at_node(profile.node_count()).is_err());
    assert!(
        profile
            .prepare_excitation([C64::new(f64::NAN, 0.0), C64::ZERO])
            .is_err()
    );
    assert!(
        profile
            .excite_power_waves([C64::from_re(f64::MAX), C64::ZERO])
            .is_err()
    );
    let mut incoherent = gas;
    incoherent.density *= 1.1;
    assert!(IsentropicMeanFlow::new(&incoherent, &bulge, 0.0).is_err());
    let many_spans =
        RadiusProfile::spline(vec![(0.0, 0.012), (0.5, 0.013), (1.0, 0.014)], None).unwrap();
    let mean = IsentropicMeanFlow::new(&gas, &many_spans, 0.0).unwrap();
    assert!(
        mean.prepare(omega(), RiccatiSteps(MAX_RICCATI_STEPS))
            .is_err()
    );
}

#[test]
fn undamped_round_trip_resonance_refuses_and_profile_identity_is_retained() {
    let gas = gas();
    let shape = RadiusProfile::exponential(0.012, 0.0, 0.83).unwrap();
    let flow = IsentropicMeanFlow::new(&gas, &shape, 0.4 * gas.sound_speed).unwrap();
    let resonance = core::f64::consts::PI * gas.sound_speed * (1.0 - 0.4 * 0.4) / shape.length();
    let prepared = flow.prepare(resonance, RiccatiSteps(32)).unwrap();
    let zero = FlowImpedance::try_new(C64::ZERO).unwrap();
    assert!(matches!(
        prepared.prepare_terminations([zero; 2], 1e8),
        Err(DuctError::Singular)
    ));
    assert!(prepared.prepare_terminations([zero; 2], 0.9).is_err());
    let matched = prepared
        .reference_impedances_j_s_kg2()
        .map(|r| FlowImpedance::try_new(C64::from_re(r)).unwrap());
    let loads = prepared.prepare_terminations(matched, 1.000_001).unwrap();
    assert!(core::ptr::eq(loads.profile().mean_flow(), &flow));
    assert!(core::ptr::eq(loads.profile().mean_flow().source(), &gas));
    assert!(loads.solve([C64::new(f64::NAN, 0.0), C64::ZERO]).is_err());
    assert!(loads.solve([C64::from_re(f64::MAX); 2]).is_err());
}
