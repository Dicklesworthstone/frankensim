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
use crate::dg_port::DgTerminalChannelInput;

const INLET: PortId = PortId(7);
const TERMINAL: PortId = PortId(3);
const DT: f64 = 1.0 / 48_000.0;
const WORK: usize = 10_000_000;
const RADIUS: f64 = 0.007;
const LENGTH: f64 = 0.08;
const MODAL_OMEGA: f64 = 1800.0;
const MODAL_DRAG: f64 = 3.0;
const TERMINAL_OMEGA: f64 = 1300.0;

pub(super) fn pole_preparation() -> super::poles::DgPolePreparation {
    super::poles::DgPolePreparation { reference_rate_per_s: 1800.0,
        maximum_dimension: 256, maximum_payload_bytes: 64 << 20,
        maximum_work: 100_000_000 }
}
pub(super) fn pole_options() -> super::poles::PoleOptions {
    super::poles::PoleOptions { spectral_scale: 1800.0, step_tolerance: 1e-12,
        residual_tolerance: 1e-11, max_iterations: 32, max_dimension: 256,
        max_dense_bytes: 64 << 20, max_work: 1_000_000_000 }
}

#[test]
fn matched_terminal_pole_converges_to_the_physical_loaded_oscillator() {
    // G1: a uniform duct terminated by its characteristic impedance presents
    // Zc at the inlet. The physical oscillator has c=MODAL_DRAG+b^2*Zc;
    // solve that quadratic directly, with no peak/linewidth-to-pole conversion.
    let gate = CancelGate::new_clock_free();
    let system = oscillator();
    let (zc, _, _, column) = coefficients();
    let rate = 0.5 * (MODAL_DRAG + column * column * zc);
    let expected = C64::new(-rate, -det::sqrt(MODAL_OMEGA.powi(2) - rate * rate));
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let source = MimoModel { ports: 1, poles: vec![], residues: vec![],
            d: vec![zc], e: vec![0.0], reciprocal: true };
        let terminal = PreparedMimoLoad::new(&source, DT, MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 33, max_states: 8, max_memory_bytes: 1 << 20,
        }).unwrap();
        let duct = DgDuct::new(&Duct { segments: vec![Segment::Cylinder {
            radius: RADIUS, length: LENGTH }] }, &gas(), DgMesh { elements, order: 2 },
            DgBoundary::Port(INLET), DgBoundary::Port(TERMINAL)).unwrap();
        let load = DgPortLoad::new_with_terminal(duct, DT, 1e-11,
            MAX_DG_PROPAGATOR_DOFS, terminal).unwrap();
        let original = load.state().clone();
        let prepared = system.prepare_resonance_poles_with_load(
            &load, &[column], pole_preparation(), &gate).unwrap();
        let response = system.harmonic_response_with_load(
            -expected.im, &[C64::from_re(1e-4)], &load, &[column], WORK, &gate).unwrap();
        let seed = prepared.seed_from_harmonic_response(&response).unwrap();
        let pole = prepared.resonance_pole(C64::new(-40.0, expected.im),
            &seed, &pole_options(), &gate).unwrap();
        let error = (pole.laplace_pole_per_s - expected).abs() / expected.abs();
        let impedance_error = (pole.inlet_pressure_pa - pole.inlet_flow_m3_s.scale(zc)).abs()
            / pole.inlet_pressure_pa.abs();
        errors.push((error, impedance_error));
        assert!(pole.relative_residual < 1e-11);
        assert_eq!(load.state(), &original);
        assert_eq!(system.samples_rendered(), 0);
    }
    // The target is the continuum characteristic load. Its finite DG image
    // keeps a separate spatial error; no numerical eigen-residual certifies it.
    for pair in errors.windows(2) {
        assert!(pair[1].0 < 0.3 * pair[0].0 + 1e-10);
        assert!(pair[1].1 < 0.3 * pair[0].1 + 1e-9);
    }
    assert!(errors[2].0 < 5e-7);
    assert!(errors[2].1 < 1e-5);
}

#[test]
fn a_dark_real_material_relaxation_survives_the_complete_loaded_pencil() {
    // G1/G3: this positive-storage Maxwell coordinate has zero mechanical
    // projection. Its exact pole -1/tau must not disappear during elimination.
    let gate = CancelGate::new_clock_free();
    let system = oscillator().with_material_relaxation(vec![fs_phs::RelaxationBranch {
        projection: vec![0.0, 0.0], stiffness: 3.0, relaxation_time_s: 0.002,
    }], None, &gate).unwrap();
    let load = loaded(1, 25.0);
    let prepared = system.prepare_resonance_poles_with_load(&load, &[coefficients().3],
        pole_preparation(), &gate).unwrap();
    let memory = prepared.coordinates().material_memory.start;
    let mut seed = vec![C64::ZERO; prepared.dimension()];
    seed[memory] = C64::ONE;
    let pole = prepared.resonance_pole(C64::from_re(-490.0), &seed,
        &pole_options(), &gate).unwrap();
    assert!((pole.laplace_pole_per_s - C64::from_re(-500.0)).abs() < 1e-10);
    assert_eq!(pole.laplace_pole_per_s.im, 0.0);
    assert_eq!(pole.inlet_flow_m3_s, C64::ZERO);
    assert_eq!(pole.inlet_pressure_pa, C64::ZERO);
    assert_eq!(pole.state[memory], C64::ONE);
    assert!(pole.state.iter().enumerate().all(|(i, value)|
        i == memory || *value == C64::ZERO));

    let mut small = pole_preparation();
    small.maximum_payload_bytes = 1;
    assert!(system.prepare_resonance_poles_with_load(&load, &[coefficients().3], small, &gate).is_err());
    small = pole_preparation(); small.maximum_dimension = prepared.dimension() - 1;
    assert!(system.prepare_resonance_poles_with_load(&load, &[coefficients().3], small, &gate).is_err());
    let mut limited = pole_options(); limited.max_work = 1;
    assert!(prepared.resonance_pole(C64::from_re(-490.0), &seed, &limited, &gate).is_err());
    assert!(prepared.resonance_pole(C64::new(-1.0, -1e6), &seed, &pole_options(), &gate).is_err());
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(prepared.resonance_pole(C64::from_re(-490.0), &seed, &pole_options(), &cancelled).is_err());
    assert!(system.prepare_resonance_poles_with_load(&load, &[coefficients().3],
        pole_preparation(), &cancelled).is_err());
    assert_eq!(system.samples_rendered(), 0);
    assert_eq!(load.state().step_count(), 0);
    assert!(load.terminal().unwrap().pole_states().iter().all(|z| *z == C64::ZERO));
}

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

pub(super) fn gap_kind(circular: bool) -> DgTerminalChannelKind {
    if circular {
        DgTerminalChannelKind::Circular {
            radius_m: 0.001, length_m: 0.02, radial_cells: 2, maximum_mach: 0.1,
        }
    } else {
        DgTerminalChannelKind::Slit {
            width_m: 0.01, height_m: 0.0005, length_m: 0.02,
            transverse_cells: 4, maximum_mach: 0.1,
        }
    }
}

// Independent 2x2 solutions derived from geometry and Newton's momentum law.
// The slit uses reflection symmetry of FOUR original cells. The circular
// case retains both unequal annular masses and the zero-shear axis. Neither
// this oracle nor its energy integrals read the channel owner's coefficients.
pub(super) fn gap_oracle(omega: f64, circular: bool) -> (C64, Vec<C64>, Vec<f64>, Vec<f64>) {
    let gas = gas();
    let s = C64::new(0.0, -omega);
    let length = 0.02;
    if circular {
        let area = core::f64::consts::PI * 0.0005_f64.powi(2);
        let mass = gas.density * length * area;
        let shear = 2.0 * core::f64::consts::PI * gas.dynamic_viscosity * length;
        let sm = s.scale(mass);
        let determinant = (C64::from_re(shear) + sm)
            * (C64::from_re(5.0 * shear) + sm.scale(3.0))
            - C64::from_re(shear * shear);
        let inner = (C64::from_re(8.0 * shear) + sm.scale(3.0)).scale(area) / determinant;
        let outer = (C64::from_re(4.0 * shear) + sm.scale(3.0)).scale(area) / determinant;
        let impedance = C64::from_re(1.0) / (inner + outer.scale(3.0)).scale(area);
        (impedance, vec![inner, outer], vec![mass, 3.0 * mass], vec![0.0, shear, 4.0 * shear])
    } else {
        let area = 0.01 * 0.0005 / 4.0;
        let mass = gas.density * length * area;
        let shear = gas.dynamic_viscosity * length * 0.01 / (0.0005 / 4.0);
        let sm = s.scale(mass);
        let determinant = (C64::from_re(3.0 * shear) + sm)
            * (C64::from_re(shear) + sm) - C64::from_re(shear * shear);
        let edge = (C64::from_re(2.0 * shear) + sm).scale(area) / determinant;
        let centre = (C64::from_re(4.0 * shear) + sm).scale(area) / determinant;
        let impedance = C64::from_re(1.0) / (edge + centre).scale(2.0 * area);
        (impedance, vec![edge, centre, centre, edge], vec![mass; 4],
            vec![2.0 * shear, shear, shear, shear, 2.0 * shear])
    }
}

fn reference_with_gap(omega: f64, gap_impedance: C64, force: C64) -> (C64, C64, C64, C64) {
    let (zc, resistance, gain, column) = coefficients();
    let s = C64::new(0.0, -omega);
    let terminal = C64::from_re(resistance)
        + s.scale(gain) / (s * s + s.scale(50.0) + C64::from_re(TERMINAL_OMEGA.powi(2)));
    let series = terminal + gap_impedance;
    let (sin, cos) = det::sin_cos(omega * LENGTH / gas().sound_speed);
    let transfer = C64::from_re(cos) + C64::new(0.0, -sin / zc) * series;
    let inlet = (series.scale(cos) + C64::new(0.0, -zc * sin)) / transfer;
    let q = force / (C64::from_re(MODAL_OMEGA.powi(2) - omega * omega)
        + s.scale(MODAL_DRAG) + s.scale(column * column) * inlet);
    let inflow = (s * q).scale(column);
    let outward = inflow / transfer;
    (q, inflow, outward, terminal * outward)
}

#[test]
fn slit_complex_impedance_converges_to_continuous_unsteady_poiseuille_flow() {
    let omega = 1000.0;
    let material = gas();
    let s = C64::new(0.0, -omega);
    let (width, height, length) = (0.01, 0.0005, 0.02);
    // Exact no-slip solution of rho*s*v=dp/L+mu*v'' on [-h/2,h/2]:
    // Q/dp = w*h/(rho*L*s) * (1-tanh(chi)/chi), chi=h/2*sqrt(rho*s/mu).
    // Evaluate complex tanh directly with real exponential/trigonometric
    // identities, independently of the FV operator and coupled LU.
    let chi = s.scale(material.density / material.dynamic_viscosity).sqrt().scale(0.5 * height);
    let exponential = det::exp(2.0 * chi.re);
    let sinh = 0.5 * (exponential - 1.0 / exponential);
    let cosh = 0.5 * (exponential + 1.0 / exponential);
    let (sin, cos) = det::sin_cos(2.0 * chi.im);
    let tanh = C64::new(sinh, sin).scale(1.0 / (cosh + cos));
    let admittance = (C64::ONE - tanh / chi)
        .scale(width * height / (material.density * length)) / s;
    let expected = C64::ONE / admittance;
    assert!(expected.re > 0.0 && expected.im < 0.0);
    let system = oscillator();
    let (_, _, _, column) = coefficients();
    let gate = CancelGate::new_clock_free();
    let mut errors = Vec::new();
    for transverse_cells in [8, 16, 32] {
        let load = loaded(2, 25.0).with_terminal_channels(vec![DgTerminalChannelInput {
            port: TERMINAL, kind: DgTerminalChannelKind::Slit {
                width_m: width, height_m: height, length_m: length,
                transverse_cells, maximum_mach: 0.1,
            },
        }], 32).unwrap();
        let response = system.harmonic_response_with_load(
            omega, &[C64::new(1e-5, -4e-6)], &load, &[column], WORK, &gate,
        ).unwrap();
        let channel = response.channels[0].channel.as_ref().unwrap();
        assert_eq!(channel.velocities_m_s.len(), transverse_cells);
        let impedance = channel.pressure_drop_pa / channel.flow_m3_s;
        errors.push((impedance - expected).abs() / expected.abs());
    }
    assert!(errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1], "{errors:?}");
    assert!(errors[2] < 3e-3, "{errors:?}");
}

#[test]
fn continuous_slit_and_circular_fields_load_resonance_with_original_storage_and_heat() {
    let system = oscillator();
    let (_, _, _, column) = coefficients();
    let force = C64::new(1e-5, -4e-6);
    let gate = CancelGate::new_clock_free();
    for circular in [false, true] {
        let load = loaded(4, 25.0).with_terminal_channels(
            vec![DgTerminalChannelInput { port: TERMINAL, kind: gap_kind(circular) }], 8,
        ).unwrap();
        for omega in [600.0, 1300.0, 1800.0] {
            let (impedance, unit_velocity, masses, faces) = gap_oracle(omega, circular);
            let expected = reference_with_gap(omega, impedance, force);
            let response = system.harmonic_response_with_load(
                omega, &[force], &load, &[column], WORK, &gate,
            ).unwrap();
            complex_close(response.modal.displacement_m_sqrt_kg[0], expected.0, 3e-4);
            complex_close(response.inlet_flow_m3_s, expected.1, 3e-4);
            let exterior = response.terminal.as_ref().unwrap();
            complex_close(exterior.outward_flow_m3_s[0], expected.2, 3e-4);
            complex_close(exterior.pressure_pa[0], expected.3, 3e-4);
            assert_eq!(response.channels.len(), 1);
            let channel = &response.channels[0];
            assert_eq!(channel.port, TERMINAL);
            complex_close(channel.exterior_pressure_pa, exterior.pressure_pa[0], 1e-12);
            let field = channel.channel.as_ref().unwrap();
            assert_eq!(field.velocities_m_s.len(), unit_velocity.len());
            let drop = impedance * exterior.outward_flow_m3_s[0];
            complex_close(channel.pressure_drop_pa, drop, 2e-9);
            complex_close(channel.duct_pressure_pa - channel.exterior_pressure_pa, drop, 2e-9);
            for (&velocity, unit) in field.velocities_m_s.iter().zip(unit_velocity) {
                complex_close(velocity, unit * drop, 2e-9);
            }
            let stored = 0.25 * masses.iter().zip(&field.velocities_m_s)
                .map(|(mass, velocity)| mass * velocity.norm_sq()).sum::<f64>();
            let mut heat = 0.5 * faces[0] * field.velocities_m_s[0].norm_sq();
            for i in 1..field.velocities_m_s.len() {
                heat += 0.5 * faces[i]
                    * (field.velocities_m_s[i] - field.velocities_m_s[i - 1]).norm_sq();
            }
            heat += 0.5 * faces[field.velocities_m_s.len()]
                * field.velocities_m_s.last().unwrap().norm_sq();
            assert!(stored > 0.0 && heat > 0.0);
            assert!((field.mean_stored_energy_j - stored).abs() < 1e-12 * stored);
            assert!((field.viscous_dissipation_w - heat).abs() < 1e-12 * heat);
            assert!((field.input_power_w - heat).abs() < 2e-8 * heat);
            let complex_work = (drop * field.flow_m3_s.conj()).scale(0.5);
            assert!((complex_work.im + 2.0 * omega * stored).abs()
                < 2e-8 * complex_work.abs());
            assert!((response.mean_stored_energy_j - response.modal.mean_stored_energy_j
                - response.duct.mean_stored_energy_j - stored).abs()
                < 2e-12 * response.mean_stored_energy_j);
            let power_scale = response.authored_external_power_w.abs() + heat
                + response.terminal_power_w.abs();
            assert!(response.power_residual_w.abs() < 2e-8 * power_scale + 1e-27);
            assert!(response.interface_power_residual_w.abs() < 2e-8 * power_scale + 1e-27);
            assert!(channel.flow_residual_m3_s.abs() < 2e-9 * field.flow_m3_s.abs());
            assert!(channel.interface_power_residual_w.abs() < 2e-8 * heat);
            assert!(field.solve_relative_residual <= 1e-11);
        }
        assert_eq!(load.state().step_count(), 0);
        assert!(load.terminal_channels()[0].channel().unwrap().velocities_m_s()
            .iter().all(|value| *value == 0.0));
        assert_eq!(load.terminal_channels()[0].stored_energy_j(), 0.0);
    }
    assert_eq!(system.samples_rendered(), 0);
}

#[test]
fn original_channel_continuous_equations_are_independent_of_step_history_and_clock() {
    use crate::bernoulli_aperture::slit_channel::{SlitChannel, SlitChannelSpec, SlitExit};
    let omega = 1800.0;
    let drop = C64::new(0.003, -0.002);
    let (_, unit_velocity, _, _) = gap_oracle(omega, false);
    let velocity: Vec<_> = unit_velocity.iter().map(|&v| v * drop).collect();
    let spec = SlitChannelSpec {
        length_m: 0.02, width_m: 0.01, height_m: 0.0005, transverse_cells: 4,
        sample_period_s: DT, exit: SlitExit::StaticPressure,
    };
    let mut previous = None;
    for period in [DT, 17.0 * DT] {
        let mut channel = SlitChannel::from_gas(&gas(), SlitChannelSpec {
            sample_period_s: period, ..spec
        }).unwrap();
        channel.step(0.01).unwrap();
        let accepted = channel.velocities_m_s().to_vec();
        let energy = channel.stored_energy_j();
        let view = channel.continuous_operator().unwrap();
        let response = view.observe_harmonic(omega, &velocity, drop, WORK, 1e-11, &mut || false).unwrap();
        if let Some(previous) = previous { assert_eq!(response, previous); }
        previous = Some(response);
        assert_eq!(channel.velocities_m_s(), accepted);
        assert_eq!(channel.stored_energy_j(), energy);
        assert!(matches!(view.observe_harmonic(omega, &velocity, drop, WORK, 1e-11, &mut || true),
            Err(AcousticRealizeError::Cancelled)));
        assert!(view.observe_harmonic(omega, &velocity, drop, 0, 1e-11, &mut || false).is_err());
        assert!(view.observe_harmonic(0.0, &velocity, drop, WORK, 1e-11, &mut || false).is_err());
        let mut wrong = velocity.clone();
        wrong[1] = wrong[1].scale(1.01);
        assert!(view.observe_harmonic(omega, &wrong, drop, WORK, 1e-11, &mut || false).is_err());
    }
    let nonlinear = SlitChannel::from_gas(&gas(), SlitChannelSpec {
        exit: SlitExit::UniformKineticHead, ..spec
    }).unwrap();
    assert!(nonlinear.continuous_operator().is_err());
}

#[test]
fn rigid_terminal_retains_exterior_states_and_finite_gap_admission_is_bounded() {
    let system = oscillator();
    let (_, _, _, column) = coefficients();
    let force = C64::new(1e-5, -4e-6);
    let gate = CancelGate::new_clock_free();
    let closed = loaded(4, 25.0).with_terminal_channels(vec![DgTerminalChannelInput {
        port: TERMINAL, kind: DgTerminalChannelKind::Closed,
    }], 0).unwrap();
    let response = system.harmonic_response_with_load(
        1800.0, &[force], &closed, &[column], WORK, &gate,
    ).unwrap();
    assert_eq!(response.terminal.as_ref().unwrap().state.len(), 2);
    let closure = &response.channels[0];
    assert!(closure.channel.is_none());
    assert!(closure.outward_flow_m3_s.abs() < 1e-10 * response.inlet_flow_m3_s.abs());
    assert!(closure.pressure_drop_pa.abs() > 0.0);
    let (zc, _, _, _) = coefficients();
    let (sin, cos) = det::sin_cos(1800.0 * LENGTH / gas().sound_speed);
    let impedance = C64::new(0.0, zc * cos / sin);
    complex_close(response.inlet_pressure_pa, impedance * response.inlet_flow_m3_s, 3e-4);

    let finite = loaded(1, 25.0).with_terminal_channels(vec![DgTerminalChannelInput {
        port: TERMINAL, kind: gap_kind(false),
    }], 8).unwrap();
    assert!(system.harmonic_response_with_load(1800.0, &[force], &finite, &[column], 1, &gate).is_err());
    assert!(system.harmonic_response_with_load(1800.0, &[C64::from_re(1e10)], &finite, &[column], WORK, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(system.harmonic_response_with_load(
        1800.0, &[force], &finite, &[column], WORK, &cancelled,
    ), Err(DgHarmonicCouplingError::Modal(ModalCouplingError::Cancelled))));
    let large = loaded(1, 25.0).with_terminal_channels(vec![DgTerminalChannelInput {
        port: TERMINAL, kind: DgTerminalChannelKind::Slit {
            width_m: 0.01, height_m: 0.0005, length_m: 0.02,
            transverse_cells: 1024, maximum_mach: 0.1,
        },
    }], 1024).unwrap();
    assert!(system.harmonic_response_with_load(
        1800.0, &[force], &large, &[column], usize::MAX, &gate,
    ).is_err());
}

