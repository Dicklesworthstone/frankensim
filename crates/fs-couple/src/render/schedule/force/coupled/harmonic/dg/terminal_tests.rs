use super::*;
use crate::dg_port::DgPortLoad;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_duct::dg::{DgBoundary, DgMesh, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::model::PoleTerm;

const INLET: PortId = PortId(7);
const TERMINAL: PortId = PortId(3);
const DT: f64 = 1.0 / 48_000.0;
const WORK: usize = 10_000_000;
const RADIUS: f64 = 0.007;
const LENGTH: f64 = 0.08;
const MODAL_OMEGA: f64 = 1800.0;
const MODAL_DRAG: f64 = 3.0;
const TERMINAL_OMEGA: f64 = 1300.0;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn coefficients() -> (f64, f64, f64, f64) {
    let gas = gas();
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * RADIUS * RADIUS);
    (zc, 0.02 * zc, 0.2 * zc * TERMINAL_OMEGA, det::sqrt(100.0 / zc))
}

fn oscillator() -> CoupledModalSystem {
    let model = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: MODAL_OMEGA,
        damping_ratio: MODAL_DRAG / (2.0 * MODAL_OMEGA),
        pressure_per_modal_velocity: C64::ZERO,
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    CoupledModalSystem::new(vec![model], vec![], ModalCouplingConfig {
        max_modes: 4, max_connections: 4, max_setup_terms: WORK,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-9,
    }, &CancelGate::new_clock_free()).unwrap()
}

fn loaded(elements: usize, decay: f64) -> DgPortLoad {
    let (_, resistance, gain, _) = coefficients();
    let imaginary = det::sqrt(TERMINAL_OMEGA * TERMINAL_OMEGA - decay * decay);
    // A physical damped oscillator: Zt(s)=R+K*s/(s²+2*a*s+w0²).
    // Both the proper rational source and its actual mean-output admission
    // retain the original conjugate pair; this is not a replacement integrator.
    let source = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Pair {
            pole: C64::new(-decay, imaginary), residue: C64::ZERO,
        }],
        residues: vec![vec![C64::new(0.5 * gain, 0.5 * gain * decay / imaginary)]],
        d: vec![resistance], e: vec![0.0], reciprocal: true,
    };
    let terminal = PreparedMimoLoad::new(&source, DT, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / DT],
        sampled_passivity_points: 33, max_states: 8, max_memory_bytes: 1 << 20,
    }).unwrap();
    let duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: RADIUS, length: LENGTH }] },
        &gas(), DgMesh { elements, order: 2 },
        DgBoundary::Port(INLET), DgBoundary::Port(TERMINAL),
    ).unwrap();
    DgPortLoad::new_with_terminal(
        duct, DT, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal,
    ).unwrap()
}

// Independent physical transfer in exp(-i*w*t). Its cylinder is analytical,
// and its terminal is the scalar oscillator law above. Neither the DG source
// matrix, rational state-space evaluator nor the coupled LU supplies this oracle.
fn reference(omega: f64, decay: f64, gain_scale: f64, force: C64) -> (C64, C64, C64, C64) {
    let (zc, resistance, gain, column) = coefficients();
    let s = C64::new(0.0, -omega);
    let terminal = C64::from_re(resistance)
        + s.scale(gain * gain_scale)
            / (s * s + s.scale(2.0 * decay) + C64::from_re(TERMINAL_OMEGA.powi(2)));
    let (sin, cos) = det::sin_cos(omega * LENGTH / gas().sound_speed);
    let transfer_flow = C64::from_re(cos) + C64::new(0.0, -sin / zc) * terminal;
    let inlet = (terminal.scale(cos) + C64::new(0.0, -zc * sin)) / transfer_flow;
    let q = force / (C64::from_re(MODAL_OMEGA.powi(2) - omega * omega)
        + s.scale(MODAL_DRAG) + s.scale(column * column) * inlet);
    let inflow = (s * q).scale(column);
    let outward = inflow / transfer_flow;
    (q, inflow, outward, terminal * outward)
}

fn peak_and_linewidth(mut magnitude_squared: impl FnMut(f64) -> f64) -> (f64, f64) {
    // This bracket contains one isolated structural peak. Resolve its actual
    // maximum and both half-power crossings; do not infer linewidth from a
    // lossless pole or replace its frequency-dependent radiation by a constant.
    let (mut lo, mut hi) = (1720.0, 1870.0);
    for _ in 0..64 {
        let a = lo + (hi - lo) / 3.0;
        let b = hi - (hi - lo) / 3.0;
        if magnitude_squared(a) < magnitude_squared(b) { lo = a; } else { hi = b; }
    }
    let peak = 0.5 * (lo + hi);
    let half = 0.5 * magnitude_squared(peak);
    assert!(magnitude_squared(1650.0) < half && magnitude_squared(1950.0) < half);
    let (mut lo, mut hi) = (1650.0, peak);
    for _ in 0..46 {
        let mid = 0.5 * (lo + hi);
        if magnitude_squared(mid) < half { lo = mid; } else { hi = mid; }
    }
    let lower = 0.5 * (lo + hi);
    let (mut lo, mut hi) = (peak, 1950.0);
    for _ in 0..46 {
        let mid = 0.5 * (lo + hi);
        if magnitude_squared(mid) > half { lo = mid; } else { hi = mid; }
    }
    (peak, 0.5 * (lo + hi) - lower)
}

#[test]
fn retained_terminal_resonance_and_linewidth_converge_to_the_physical_cylinder() {
    let system = oscillator();
    let (_, _, _, column) = coefficients();
    let force = C64::new(1e-5, -4e-6);
    let gate = CancelGate::new_clock_free();
    let exact = peak_and_linewidth(|w| reference(w, 25.0, 1.0, force).0.norm_sq());
    let resistive = peak_and_linewidth(|w| reference(w, 25.0, 0.0, force).0.norm_sq());
    assert!((exact.0 - resistive.0).abs() > 10.0);
    assert!((exact.1 - resistive.1).abs() > 0.1);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let load = loaded(elements, 25.0);
        let actual = peak_and_linewidth(|omega| {
            system.harmonic_response_with_load(
                omega, &[force], &load, &[column], WORK, &gate,
            ).unwrap().modal.displacement_m_sqrt_kg[0].norm_sq()
        });
        errors.push(((actual.0 - exact.0).abs(), (actual.1 - exact.1).abs()));
        assert_eq!(load.state().step_count(), 0);
        assert!(load.terminal().unwrap().pole_states().iter().all(|z| *z == C64::ZERO));
    }
    for i in 1..3 {
        assert!(errors[i].0 < 0.5 * errors[i - 1].0, "peak errors: {errors:?}");
        assert!(errors[i].1 < 0.5 * errors[i - 1].1, "linewidth errors: {errors:?}");
    }
    assert!(1200.0 * errors[2].0 / (exact.0 * core::f64::consts::LN_2) < 0.1);
    assert!(errors[2].1 < 1e-3 * exact.1, "linewidth errors: {errors:?}");
    assert_eq!(system.samples_rendered(), 0);
}

fn complex_close(actual: C64, expected: C64, relative: f64) {
    assert!((actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + 1e-27,
        "actual={actual:?}, expected={expected:?}, relative={relative}");
}

#[test]
fn lightly_damped_bare_terminal_resonance_retains_finite_original_states_and_physical_power() {
    // Strictly stable original admission excludes imaginary-axis poles. Approach
    // the bare resonator with three positive linewidths instead, while solving
    // the complete joint equations at its natural frequency without pole division.
    let system = oscillator();
    let (_, resistance, gain, column) = coefficients();
    let omega = TERMINAL_OMEGA;
    let s = C64::new(0.0, -omega);
    let force = C64::new(1e-5, -4e-6);
    let gate = CancelGate::new_clock_free();
    for decay in [25.0, 0.25, 0.0025] {
        let load = loaded(4, decay);
        let response = system.harmonic_response_with_load(
            omega, &[force], &load, &[column], WORK, &gate,
        ).unwrap();
        let terminal = response.terminal.as_ref().unwrap();
        assert_eq!(terminal.ports, [TERMINAL]);
        assert_eq!(terminal.state.len(), 2);
        let flow = terminal.outward_flow_m3_s[0];
        let q = flow / (s * s + s.scale(2.0 * decay) + C64::from_re(omega * omega));
        let velocity = s * q;
        let imaginary = det::sqrt(omega * omega - decay * decay);
        // The original real chart is [Re(z), -Im(z)] for each pair.
        complex_close(terminal.state[0], (s + C64::from_re(decay)) * q, 2e-8);
        complex_close(terminal.state[1], q.scale(-imaginary), 2e-8);
        let pressure = flow.scale(resistance) + velocity.scale(gain);
        complex_close(terminal.pressure_pa[0], pressure, 2e-8);
        let exact = reference(omega, decay, 1.0, force);
        complex_close(response.modal.displacement_m_sqrt_kg[0], exact.0, 2e-4);
        complex_close(response.inlet_flow_m3_s, exact.1, 2e-4);
        complex_close(flow, exact.2, 2e-4);
        complex_close(pressure, exact.3, 2e-4);

        let power = 0.5 * resistance * flow.norm_sq() + decay * gain * velocity.norm_sq();
        let stored = 0.25 * gain * (velocity.norm_sq() + omega * omega * q.norm_sq());
        assert!(power > 0.0 && power.is_finite() && stored > 0.0 && stored.is_finite());
        for actual in [terminal.mean_input_power_w[0], terminal.input_power_w, response.terminal_power_w] {
            assert!((actual - power).abs() <= 2e-8 * power + 1e-28,
                "terminal power={actual:e}, physical loss={power:e}");
        }
        let full_power_scale = response.authored_external_power_w.abs() + power;
        assert!(response.power_residual_w.abs() <= 2e-8 * full_power_scale + 1e-27);
        assert!(response.interface_power_residual_w.abs() <= 2e-8 * full_power_scale + 1e-27);
        assert!(response.solve_relative_residual <= 1e-11);
        // The generic owner does not assert a storage certificate for an
        // arbitrary rational fit. This fixture's extra storage was derived
        // independently from its physical oscillator, not added to that report.
        assert_eq!(response.mean_stored_energy_j,
            response.modal.mean_stored_energy_j + response.duct.mean_stored_energy_j);
        assert_eq!(load.terminal().unwrap().pole_states(), [C64::ZERO]);
        assert_eq!(load.state().step_count(), 0);
    }
    assert_eq!(system.samples_rendered(), 0);
}
