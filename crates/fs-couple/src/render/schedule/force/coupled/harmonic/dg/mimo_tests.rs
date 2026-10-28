use super::*;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_duct::dg::{DgBoundary, DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_vfit::mimo::{MimoModel, load::MimoLoadLimits};
use fs_vfit::model::PoleTerm;

fn close(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 3e-9 * a.abs().max(b.abs()) + 1e-21,
        "{a:?} != {b:?}"
    );
}
fn fixture(reverse: bool) -> (CoupledModalSystem, DgPortLoad) {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let radius = 0.003;
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: 0.08,
            }],
        },
        &gas,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Port(PortId(8)),
    )
    .unwrap()
    .with_side_ports(&[DgSidePort {
        id: PortId(9),
        downstream_element: 1,
    }])
    .unwrap();
    let mut source = MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Real {
            pole: -500.0,
            residue: 0.0,
        }],
        residues: vec![
            [0.25, 0.06, 0.06, 0.2]
                .into_iter()
                .map(|v| C64::from_re(v * zc * 300.0))
                .collect(),
        ],
        d: vec![0.15 * zc, 0.03 * zc, 0.03 * zc, 0.2 * zc],
        e: vec![0.0; 4],
        reciprocal: true,
    };
    let ports = if reverse {
        source.d.reverse();
        source.residues[0].reverse();
        vec![PortId(8), PortId(9)]
    } else {
        vec![PortId(9), PortId(8)]
    };
    let dt = 1.0 / 48000.0;
    let terminal = PreparedMimoLoad::new(
        &source,
        dt,
        MimoLoadLimits {
            band_rad_s: [100.0, core::f64::consts::PI / dt],
            sampled_passivity_points: 33,
            max_states: 8,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    let load = DgPortLoad::new_with_terminal_ports(
        duct,
        dt,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        ports,
        terminal,
    )
    .unwrap();
    let component = ModalAcousticTimeModel::try_new(
        48000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 1400.0,
            damping_ratio: 0.015,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let system = CoupledModalSystem::new(
        vec![component],
        vec![],
        ModalCouplingConfig {
            max_modes: 4,
            max_connections: 4,
            max_setup_terms: 1_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    (system, load)
}

#[test]
fn original_mutual_terminal_matrix_and_port_permutation_match_physical_condensation() {
    // The reference eliminates TWO physical outward flows using a closed-form
    // 2x2 inverse. Production retains original terminal states and all three DG
    // ports. Reversed authored channel order exercises every original mapping.
    let (system, load) = fixture(false);
    let (reordered, other) = fixture(true);
    let gate = CancelGate::new_clock_free();
    let force = C64::new(0.002, -0.001);
    let b = 0.002;
    let ids = [PortId(7), PortId(9), PortId(8)];
    let initial = load.state().clone();
    let history = load.terminal().unwrap().pole_states().to_vec();
    for omega in [400.0, 1200.0, 1800.0, 3500.0] {
        let s = C64::new(0.0, -omega);
        let columns: Vec<_> = ids
            .iter()
            .map(|&id| {
                load.duct()
                    .harmonic_response(omega, &[(id, C64::ONE)], 10_000_000, 1e-10, &mut || false)
                    .unwrap()
            })
            .collect();
        let z = |i: usize, j: usize| {
            columns[j]
                .ports
                .iter()
                .find(|p| p.id == ids[i])
                .unwrap()
                .pressure_pa
        };
        let terminal = load.terminal().unwrap().source_model().eval(s);
        let a = z(1, 1) + terminal[0];
        let c = z(1, 2) + terminal[1];
        let d = z(2, 1) + terminal[2];
        let e = z(2, 2) + terminal[3];
        let determinant = a * e - c * d;
        let alpha = [
            (e * z(1, 0) - c * z(2, 0)) / determinant,
            (a * z(2, 0) - d * z(1, 0)) / determinant,
        ];
        let zin = z(0, 0) - z(0, 1) * alpha[0] - z(0, 2) * alpha[1];
        let displacement = force
            / (C64::new(1400.0_f64.powi(2) - omega * omega, -42.0 * omega) + s.scale(b * b) * zin);
        let inflow = s.scale(b) * displacement;
        let expected = [alpha[0] * inflow, alpha[1] * inflow];
        let actual = system
            .harmonic_response_with_load(omega, &[force], &load, &[b], 10_000_000, &gate)
            .unwrap();
        let permuted = reordered
            .harmonic_response_with_load(omega, &[force], &other, &[b], 10_000_000, &gate)
            .unwrap();
        close(actual.modal.displacement_m_sqrt_kg[0], displacement);
        close(actual.inlet_pressure_pa, zin * inflow);
        close(actual.inlet_flow_m3_s, inflow);
        close(actual.inlet_pressure_pa, permuted.inlet_pressure_pa);
        let t = actual.terminal.as_ref().unwrap();
        let reversed = permuted.terminal.as_ref().unwrap();
        assert_eq!(t.ports, [PortId(9), PortId(8)]);
        assert_eq!(t.state.len(), 2);
        let mut power = 0.0;
        for i in 0..2 {
            let pressure = terminal[i * 2] * expected[0] + terminal[i * 2 + 1] * expected[1];
            close(t.outward_flow_m3_s[i], expected[i]);
            close(t.pressure_pa[i], pressure);
            close(t.state[i], expected[i] / (s + C64::from_re(500.0)));
            close(t.outward_flow_m3_s[i], reversed.outward_flow_m3_s[1 - i]);
            close(t.state[i], reversed.state[1 - i]);
            power += 0.5 * (pressure * expected[i].conj()).re;
        }
        for i in 0..actual.duct.pressure_pa.len() {
            let p = columns[0].pressure_pa[i] * inflow
                - columns[1].pressure_pa[i] * expected[0]
                - columns[2].pressure_pa[i] * expected[1];
            close(actual.duct.pressure_pa[i], p);
        }
        assert!((actual.terminal_power_w - power).abs() <= 3e-9 * power.abs() + 1e-24);
        assert!(
            (actual.mean_stored_energy_j
                - actual.modal.mean_stored_energy_j
                - actual.duct.mean_stored_energy_j)
                .abs()
                < 1e-20
        );
        assert!(
            actual.power_residual_w.abs() < 3e-9 * actual.authored_external_power_w.abs() + 1e-24
        );
        assert!(
            actual.interface_power_residual_w.abs()
                < 3e-9 * actual.authored_external_power_w.abs() + 1e-24
        );
    }
    assert_eq!(load.state(), &initial);
    assert_eq!(load.terminal().unwrap().pole_states(), history);
    assert_eq!(system.samples_rendered(), 0);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.harmonic_response_with_load(1000.0, &[force], &load, &[b], 10_000_000, &cancelled),
        Err(DgHarmonicCouplingError::Modal(
            ModalCouplingError::Cancelled
        ))
    ));
    for omega in [99.0, 200000.0] {
        assert!(
            system
                .harmonic_response_with_load(omega, &[force], &load, &[b], 10_000_000, &gate)
                .is_err()
        );
    }
    assert!(
        system
            .harmonic_response_with_load(1000.0, &[force], &load, &[b], 1, &gate)
            .is_err()
    );
}
