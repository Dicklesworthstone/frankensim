use super::*;
use crate::dg_port::DgBranch;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::MimoLoadLimits;
use fs_vfit::model::PoleTerm;

const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const SIDE: PortId = PortId(19);
const SECOND: PortId = PortId(11);
const DT: f64 = 1.0 / 48_000.0;
const WORK: usize = 100_000_000;
const RADIUS: f64 = 0.007;
const LENGTH: f64 = 0.16;
const BRANCH_RADIUS: f64 = 0.004;
const BRANCH_LENGTH: f64 = 0.35;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn zc(radius: f64) -> f64 {
    let gas = gas();
    gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius)
}
fn column() -> f64 {
    det::sqrt(100.0 / zc(RADIUS))
}
fn oscillator() -> CoupledModalSystem {
    let component = ModalAcousticTimeModel::try_new(
        48000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 1800.0,
            damping_ratio: 3.0 / 3600.0,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    CoupledModalSystem::new(
        vec![component],
        vec![],
        ModalCouplingConfig {
            max_modes: 4,
            max_connections: 4,
            max_setup_terms: WORK,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn cylinder(radius: f64, length: f64, elements: usize, order: usize) -> DgDuct {
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        &gas(),
        DgMesh { elements, order },
        DgBoundary::Port(INLET),
        DgBoundary::Port(OUTLET),
    )
    .unwrap()
}
fn source(mutual: bool, strictly_proper: bool) -> MimoModel {
    let z = zc(RADIUS);
    let zb = zc(BRANCH_RADIUS);
    let d = vec![
        0.04 * zb,
        0.003 * det::sqrt(z * zb),
        0.003 * det::sqrt(z * zb),
        0.02 * z,
    ];
    MimoModel {
        ports: 2,
        poles: if mutual {
            vec![PoleTerm::Real {
                pole: -600.0,
                residue: 0.0,
            }]
        } else {
            vec![]
        },
        residues: if mutual {
            vec![d.iter().map(|&v| C64::from_re(300.0 * v)).collect()]
        } else {
            vec![]
        },
        d: if strictly_proper {
            vec![0.0; 4]
        } else if mutual {
            d
        } else {
            vec![d[0], 0.0, 0.0, d[3]]
        },
        e: vec![0.0; 4],
        reciprocal: true,
    }
}
fn terminal(source: &MimoModel) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        source,
        DT,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 65,
            max_states: 16,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}
fn loaded(rung: usize, mutual: bool, strictly_proper: bool, reverse: bool) -> DgPortLoad {
    let main = cylinder(RADIUS, LENGTH, 2 * rung, 2)
        .with_side_ports(&[DgSidePort {
            id: SIDE,
            downstream_element: rung,
        }])
        .unwrap();
    let branch = cylinder(BRANCH_RADIUS, BRANCH_LENGTH, rung, 2);
    let mut model = source(mutual, strictly_proper);
    let ports = if reverse {
        model.d.reverse();
        for residue in &mut model.residues {
            residue.reverse();
        }
        vec![OUTLET, SIDE]
    } else {
        vec![SIDE, OUTLET]
    };
    DgPortLoad::new_with_terminal_branches(
        main,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        ports,
        terminal(&model),
        vec![DgBranch {
            port: SIDE,
            duct: branch,
        }],
    )
    .unwrap()
}

// Independent physical cylinder matrices in exp(-i omega t), with both flow
// coordinates positive downstream. No DG matrix, coupled LU or reduced fitted
// impedance supplies the reference. The 2x2 inverse solves only the two actual
// external flows at unit junction pressure, including mutual exterior reaction.
fn transfer(omega: f64, length: f64, radius: f64) -> (f64, C64, C64) {
    let (sin, cos) = det::sin_cos(omega * length / gas().sound_speed);
    (
        cos,
        C64::new(0.0, -zc(radius) * sin),
        C64::new(0.0, -sin / zc(radius)),
    )
}
fn reference(
    omega: f64,
    mutual: bool,
    strictly_proper: bool,
    force: C64,
) -> (C64, C64, C64, C64, [C64; 2], [C64; 2]) {
    let s = C64::new(0.0, -omega);
    let original = source(mutual, strictly_proper);
    let z: Vec<_> = (0..4)
        .map(|i| {
            C64::from_re(original.d[i])
                + if mutual {
                    original.residues[0][i] / (s + C64::from_re(600.0))
                } else {
                    C64::ZERO
                }
        })
        .collect();
    let (a, b, c) = transfer(omega, LENGTH / 2.0, RADIUS);
    let (ab, bb, cb) = transfer(omega, BRANCH_LENGTH, BRANCH_RADIUS);
    let m = [
        z[0].scale(ab) + bb,
        z[1].scale(ab),
        z[2].scale(a),
        z[3].scale(a) + b,
    ];
    let determinant = m[0] * m[3] - m[1] * m[2];
    let u = [(m[3] - m[1]) / determinant, (m[0] - m[2]) / determinant];
    let p = [z[0] * u[0] + z[1] * u[1], z[2] * u[0] + z[3] * u[1]];
    let branch_flow = cb * p[0] + u[0].scale(ab);
    let y = branch_flow + c * p[1] + u[1].scale(a);
    let input_per_pressure = c + y.scale(a);
    let impedance = (C64::from_re(a) + b * y) / input_per_pressure;
    let displacement = force
        / (C64::from_re(1800.0_f64.powi(2) - omega * omega)
            + s.scale(3.0)
            + s.scale(column().powi(2)) * impedance);
    let input = (s * displacement).scale(column());
    let junction_pressure = input / input_per_pressure;
    (
        displacement,
        input,
        junction_pressure,
        junction_pressure * branch_flow,
        [u[0] * junction_pressure, u[1] * junction_pressure],
        [p[0] * junction_pressure, p[1] * junction_pressure],
    )
}
#[track_caller]
fn close(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + 1e-24,
        "actual={actual:?}, expected={expected:?}, relative={relative}"
    );
}
fn peak_width(mut norm: impl FnMut(f64) -> f64) -> (f64, f64) {
    let (mut lo, mut hi) = (1700.0, 1950.0);
    for _ in 0..60 {
        let a = lo + (hi - lo) / 3.0;
        let b = hi - (hi - lo) / 3.0;
        if norm(a) < norm(b) {
            lo = a;
        } else {
            hi = b;
        }
    }
    let peak = 0.5 * (lo + hi);
    let half = 0.5 * norm(peak);
    assert!(norm(1600.0) < half && norm(2100.0) < half);
    let (mut lo, mut hi) = (1600.0, peak);
    for _ in 0..42 {
        let mid = 0.5 * (lo + hi);
        if norm(mid) < half {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let lower = 0.5 * (lo + hi);
    let (mut lo, mut hi) = (peak, 2100.0);
    for _ in 0..42 {
        let mid = 0.5 * (lo + hi);
        if norm(mid) > half {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (peak, 0.5 * (lo + hi) - lower)
}
fn check_power(response: &DgLoadedModalHarmonicResponse) {
    let scale = response.authored_external_power_w.abs()
        + response.total_dissipation_w()
        + response.terminal_power_w.abs();
    assert!(response.power_residual_w.abs() <= 2e-8 * scale + 1e-25);
    assert!(response.interface_power_residual_w.abs() <= 2e-8 * scale + 1e-25);
    assert!(response.solve_relative_residual <= 1e-11);
    let storage = response.branches.iter().fold(
        response.modal.mean_stored_energy_j + response.duct.mean_stored_energy_j,
        |sum, branch| sum + branch.duct.mean_stored_energy_j,
    );
    assert_eq!(response.mean_stored_energy_j, storage);
    for branch in &response.branches {
        close(branch.bore_pressure_pa, branch.inlet_pressure_pa, 1e-9);
        assert!(
            branch.junction_power_residual_w.abs()
                <= 2e-8 * (branch.bore_power_w.abs() + branch.branch_power_w.abs()) + 1e-25
        );
    }
}

#[test]
fn finite_branch_peak_and_linewidth_converge_to_analytic_cylinder_network() {
    let system = oscillator();
    let force = C64::new(1e-5, -4e-6);
    let gate = CancelGate::new_clock_free();
    let exact = peak_width(|w| reference(w, false, false, force).0.norm_sq());
    let mut errors = Vec::new();
    for rung in [1, 2, 4] {
        let load = loaded(rung, false, false, false);
        let actual = peak_width(|omega| {
            system
                .harmonic_response_with_load(omega, &[force], &load, &[column()], WORK, &gate)
                .unwrap()
                .modal
                .displacement_m_sqrt_kg[0]
                .norm_sq()
        });
        errors.push(((actual.0 - exact.0).abs(), (actual.1 - exact.1).abs()));
        let response = system
            .harmonic_response_with_load(actual.0, &[force], &load, &[column()], WORK, &gate)
            .unwrap();
        assert_eq!(response.branches.len(), 1);
        assert_eq!(
            response.branches[0].duct.energy_coordinates.len(),
            load.branches()[0].duct().total_dofs()
        );
        check_power(&response);
        assert_eq!(load.branches()[0].state().step_count(), 0);
    }
    for i in 1..3 {
        assert!(
            errors[i].0 < 0.5 * errors[i - 1].0,
            "peak/width errors={errors:?}"
        );
        assert!(
            errors[i].1 < 0.5 * errors[i - 1].1,
            "peak/width errors={errors:?}"
        );
    }
    assert!(
        1200.0 * errors[2].0 / (exact.0 * core::f64::consts::LN_2) < 0.1,
        "peak/width errors={errors:?}, analytic={exact:?}"
    );
    assert!(
        errors[2].1 < 0.001 * exact.1,
        "peak/width errors={errors:?}, analytic={exact:?}"
    );
}

#[test]
fn original_mutual_terminal_uses_branch_outlet_and_preserves_live_histories() {
    let system = oscillator();
    let force = C64::new(1e-5, 3e-6);
    let gate = CancelGate::new_clock_free();
    for strictly_proper in [false, true] {
        let load = loaded(12, true, strictly_proper, false);
        let reverse = loaded(12, true, strictly_proper, true);
        // Local junction pressure is sensitive near the branch half-wave pole.
        // Resolve that observable independently of the structural-peak ladder;
        // do not loosen its physical tolerance to cover an underresolved field.
        let expected_pressure = reference(3100.0, true, strictly_proper, force).2;
        let coarse = loaded(4, true, strictly_proper, false);
        let middle = loaded(8, true, strictly_proper, false);
        let mut pressure_errors = Vec::new();
        for original in [&coarse, &middle, &load] {
            let response = system
                .harmonic_response_with_load(3100.0, &[force], original, &[column()], WORK, &gate)
                .unwrap();
            pressure_errors.push(
                (response.branches[0].inlet_pressure_pa - expected_pressure).abs()
                    / expected_pressure.abs(),
            );
        }
        for i in 1..3 {
            assert!(
                pressure_errors[i] < 0.5 * pressure_errors[i - 1],
                "proper={strictly_proper}, branch half-wave pressure errors={pressure_errors:?}"
            );
        }
        assert!(
            pressure_errors[2] < 3e-4,
            "proper={strictly_proper}, branch half-wave pressure errors={pressure_errors:?}"
        );
        let main_state = load.state().clone();
        let branch_state = load.branches()[0].state().clone();
        let terminal_state = load.terminal().unwrap().pole_states().to_vec();
        let frames = load.last_branch_frames().to_vec();
        for omega in [400.0, 1200.0, 1800.0, 3100.0] {
            let actual = system
                .harmonic_response_with_load(omega, &[force], &load, &[column()], WORK, &gate)
                .unwrap();
            let reordered = system
                .harmonic_response_with_load(omega, &[force], &reverse, &[column()], WORK, &gate)
                .unwrap();
            let expected = reference(omega, true, strictly_proper, force);
            close(actual.modal.displacement_m_sqrt_kg[0], expected.0, 2e-4);
            close(actual.inlet_flow_m3_s, expected.1, 2e-4);
            close(actual.branches[0].inlet_pressure_pa, expected.2, 3e-4);
            close(actual.branches[0].junction_flow_m3_s, expected.3, 3e-4);
            close(actual.inlet_pressure_pa, reordered.inlet_pressure_pa, 1e-9);
            let load_response = actual.terminal.as_ref().unwrap();
            let reversed = reordered.terminal.as_ref().unwrap();
            assert_eq!(load_response.ports, [SIDE, OUTLET]);
            assert_eq!(reversed.ports, [OUTLET, SIDE]);
            for i in 0..2 {
                close(load_response.outward_flow_m3_s[i], expected.4[i], 3e-4);
                close(load_response.pressure_pa[i], expected.5[i], 3e-4);
                close(
                    load_response.outward_flow_m3_s[i],
                    reversed.outward_flow_m3_s[1 - i],
                    1e-9,
                );
                close(load_response.state[i], reversed.state[1 - i], 1e-9);
            }
            let outlet = actual.branches[0]
                .duct
                .ports
                .iter()
                .find(|p| p.id == OUTLET)
                .unwrap();
            close(load_response.pressure_pa[0], outlet.pressure_pa, 1e-9);
            close(
                load_response.outward_flow_m3_s[0],
                -outlet.inflow_m3_s,
                1e-9,
            );
            assert!(
                (actual.branches[0].junction_flow_m3_s - load_response.outward_flow_m3_s[0]).abs()
                    > 1e-4 * actual.branches[0].junction_flow_m3_s.abs()
            );
            check_power(&actual);
        }
        assert_eq!(load.state(), &main_state);
        assert_eq!(load.branches()[0].state(), &branch_state);
        assert_eq!(load.terminal().unwrap().pole_states(), terminal_state);
        assert_eq!(load.last_branch_frames(), frames);
        // Use the existing mean-port owner's admitted coarse fixture to create
        // physical histories. This checks history independence separately from
        // the finer spatial accuracy oracle; it does not change that oracle.
        let mut warm = loaded(1, true, strictly_proper, false);
        let before = system
            .harmonic_response_with_load(1800.0, &[force], &warm, &[column()], WORK, &gate)
            .unwrap();
        for sample in 0..2 {
            warm.preview_step(1e-9).unwrap_or_else(|error|
                panic!("original mean history: proper={strictly_proper}, sample={sample}, {error:?}"))
                .accept();
        }
        let state = warm.state().clone();
        let branch = warm.branches()[0].state().clone();
        let memory = warm.terminal().unwrap().pole_states().to_vec();
        let frames = warm.last_branch_frames().to_vec();
        assert!(branch.stored_energy_j() > 0.0 && memory.iter().any(|z| z.abs() > 0.0));
        assert_eq!(
            before,
            system
                .harmonic_response_with_load(1800.0, &[force], &warm, &[column()], WORK, &gate,)
                .unwrap()
        );
        assert_eq!(warm.state(), &state);
        assert_eq!(warm.branches()[0].state(), &branch);
        assert_eq!(warm.terminal().unwrap().pole_states(), memory);
        assert_eq!(warm.last_branch_frames(), frames);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            system.harmonic_response_with_load(
                1800.0,
                &[force],
                &load,
                &[column()],
                WORK,
                &cancelled,
            ),
            Err(DgHarmonicCouplingError::Modal(
                ModalCouplingError::Cancelled
            ))
        ));
        assert!(
            system
                .harmonic_response_with_load(1800.0, &[force], &load, &[column()], 1, &gate)
                .is_err()
        );
    }
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn bare_original_branch_pole_remains_a_finite_joint_frequency_problem() {
    let main = cylinder(RADIUS, LENGTH, 2, 2)
        .with_side_ports(&[DgSidePort {
            id: SIDE,
            downstream_element: 1,
        }])
        .unwrap();
    let branch = cylinder(BRANCH_RADIUS, BRANCH_LENGTH, 1, 1);
    // Its original four-state skew generator has one nonzero singular value
    // in the 2x2 pressure/flow block. Obtain that exact spatial pole without
    // ever asking an isolated harmonic solve to invert it.
    let original = branch.continuous_operator();
    assert_eq!(original.dimension(), 4);
    let l = original.generator_per_s();
    let g = [l[2], l[3], l[6], l[7]];
    let a = g[0] * g[0] + g[1] * g[1];
    let d = g[2] * g[2] + g[3] * g[3];
    let c = g[0] * g[2] + g[1] * g[3];
    let omega = det::sqrt(0.5 * (a + d + det::sqrt((a - d).powi(2) + 4.0 * c * c)));
    let load = DgPortLoad::new_with_terminal_branches(
        main,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        vec![SIDE, OUTLET],
        terminal(&source(true, false)),
        vec![DgBranch {
            port: SIDE,
            duct: branch,
        }],
    )
    .unwrap();
    let system = oscillator();
    let force = C64::new(1e-5, -3e-6);
    let gate = CancelGate::new_clock_free();
    let exact = system
        .harmonic_response_with_load(omega, &[force], &load, &[column()], WORK, &gate)
        .unwrap();
    check_power(&exact);
    assert!(exact.branches[0].duct.mean_stored_energy_j > 0.0);
    for scale in [1.0 - 1e-7, 1.0 + 1e-7] {
        let nearby = system
            .harmonic_response_with_load(omega * scale, &[force], &load, &[column()], WORK, &gate)
            .unwrap();
        close(exact.inlet_pressure_pa, nearby.inlet_pressure_pa, 1e-4);
        close(
            exact.branches[0].junction_flow_m3_s,
            nearby.branches[0].junction_flow_m3_s,
            1e-4,
        );
        check_power(&nearby);
    }
}

#[test]
fn complete_two_branch_network_retains_local_walls_bands_and_both_orderings() {
    let radius = 0.001;
    let lengths = [0.05, 0.08];
    let resistances = [0.04 * zc(radius), 0.02 * zc(RADIUS), 0.07 * zc(radius)];
    let create = |walls: bool, reverse: bool| {
        let main = cylinder(RADIUS, LENGTH, 6, 2)
            .with_side_ports(&[
                DgSidePort {
                    id: SIDE,
                    downstream_element: 2,
                },
                DgSidePort {
                    id: SECOND,
                    downstream_element: 4,
                },
            ])
            .unwrap();
        let first = cylinder(radius, lengths[0], 2, 2);
        let second = cylinder(radius, lengths[1], 2, 2);
        let second = if walls {
            second.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
        } else {
            second
        };
        let mut branches = vec![
            DgBranch {
                port: SIDE,
                duct: first,
            },
            DgBranch {
                port: SECOND,
                duct: second,
            },
        ];
        let order = if reverse { [1, 2, 0] } else { [0, 1, 2] };
        if reverse {
            branches.reverse();
        }
        let ids = [SIDE, OUTLET, SECOND];
        let mut d = vec![0.0; 9];
        for i in 0..3 {
            d[3 * i + i] = resistances[order[i]];
        }
        let source = MimoModel {
            ports: 3,
            poles: vec![],
            residues: vec![],
            d,
            e: vec![0.0; 9],
            reciprocal: true,
        };
        DgPortLoad::new_with_terminal_branches(
            main,
            DT,
            1e-11,
            MAX_DG_PROPAGATOR_DOFS,
            order.iter().map(|&i| ids[i]).collect(),
            terminal(&source),
            branches,
        )
        .unwrap()
    };
    let system = oscillator();
    let force = C64::new(1e-5, -4e-6);
    let gate = CancelGate::new_clock_free();
    for walls in [false, true] {
        let load = create(walls, false);
        let reverse = create(walls, true);
        for omega in [700.0, 1700.0] {
            let actual = system
                .harmonic_response_with_load(omega, &[force], &load, &[column()], WORK, &gate)
                .unwrap();
            let reordered = system
                .harmonic_response_with_load(omega, &[force], &reverse, &[column()], WORK, &gate)
                .unwrap();
            close(actual.inlet_pressure_pa, reordered.inlet_pressure_pa, 2e-9);
            assert_eq!(
                actual.branches.iter().map(|b| b.port).collect::<Vec<_>>(),
                [SIDE, SECOND]
            );
            assert_eq!(
                reordered
                    .branches
                    .iter()
                    .map(|b| b.port)
                    .collect::<Vec<_>>(),
                [SECOND, SIDE]
            );
            for (i, branch) in actual.branches.iter().enumerate() {
                close(
                    branch.junction_flow_m3_s,
                    reordered.branches[1 - i].junction_flow_m3_s,
                    2e-9,
                );
                assert_eq!(
                    branch.duct.energy_coordinates.len(),
                    load.branches()[i].duct().total_dofs()
                );
                let terminal = actual.terminal.as_ref().unwrap();
                let row = terminal
                    .ports
                    .iter()
                    .position(|id| *id == branch.port)
                    .unwrap();
                let outlet = branch.duct.ports.iter().find(|p| p.id == OUTLET).unwrap();
                close(outlet.pressure_pa, terminal.pressure_pa[row], 2e-9);
                close(-outlet.inflow_m3_s, terminal.outward_flow_m3_s[row], 2e-9);
            }
            check_power(&actual);
            if walls {
                assert!(actual.branches[1].duct.mean_wall_energy_j > 0.0);
                assert!(actual.branches[1].duct.wall_dissipation_w > 0.0);
                assert_eq!(actual.branches[0].duct.wall_dissipation_w, 0.0);
                assert_eq!(actual.duct.wall_dissipation_w, 0.0);
            } else {
                // Independent scalar recursion through TWO distinct junctions.
                // This oracle is valid only for this diagonal terminal fixture;
                // the mutual test above uses the full physical 2x2 reference.
                let propagate = |length, radius, z: C64| {
                    let (a, b, c) = transfer(omega, length, radius);
                    (z.scale(a) + b) / (c * z + C64::from_re(a))
                };
                let zb0 = propagate(lengths[0], radius, C64::from_re(resistances[0]));
                let zb1 = propagate(lengths[1], radius, C64::from_re(resistances[2]));
                let distal = propagate(LENGTH / 3.0, RADIUS, C64::from_re(resistances[1]));
                let second_junction = (distal.recip() + zb1.recip()).recip();
                let middle = propagate(LENGTH / 3.0, RADIUS, second_junction);
                let first_junction = (middle.recip() + zb0.recip()).recip();
                let impedance = propagate(LENGTH / 3.0, RADIUS, first_junction);
                close(
                    actual.inlet_pressure_pa / actual.inlet_flow_m3_s,
                    impedance,
                    2e-4,
                );
            }
        }
        if walls {
            // The main duct and rational source admit this query. The actual
            // second branch is the owner whose thermal band must refuse it.
            assert!(
                system
                    .harmonic_response_with_load(
                        13_000.0,
                        &[force],
                        &load,
                        &[column()],
                        WORK,
                        &gate,
                    )
                    .is_err()
            );
        }
    }
}
