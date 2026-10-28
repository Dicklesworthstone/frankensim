use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_duct::dg::{DgBoundary, DgMesh};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;

const WORK: usize = 1_000_000;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn cylinder(elements: usize, order: usize, radius: f64, length: f64) -> DgDuct {
    DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius, length }] },
        &air(), DgMesh { elements, order },
        DgBoundary::Port(PortId(7)), DgBoundary::IdealOpen,
    ).unwrap()
}

fn oscillator(omega: f64) -> CoupledModalSystem {
    let model = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO,
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    CoupledModalSystem::new(vec![model], vec![], ModalCouplingConfig {
        max_modes: 4, max_connections: 4, max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-9,
    }, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn loaded_cavity_neck_resonance_converges_to_the_analytic_cylinder() {
    // G1/G2: the cavity-neck oscillator uses q=sqrt(L_n)*displaced_volume.
    // The independent ideal-open cylinder gives Z=-i Zc tan(omega L/c).
    // Its joint characteristic equation is 1/C-omega²L_n-omega Zc tan(kL)=0.
    // Only DG spatial resolution changes; no isolated DG impedance is inverted.
    let gas = air();
    let radius = 0.003;
    let length = 0.08;
    let area = core::f64::consts::PI * radius * radius;
    let inertance = gas.density * 0.02 / area;
    let compliance = 3e-6 / (gas.density * gas.sound_speed * gas.sound_speed);
    let zc = gas.density * gas.sound_speed / area;
    let isolated = 1.0 / det::sqrt(compliance * inertance);
    let pole = core::f64::consts::PI * gas.sound_speed / (2.0 * length);
    let characteristic = |omega: f64| {
        let (sin, cos) = det::sin_cos(omega * length / gas.sound_speed);
        1.0 / compliance - omega * omega * inertance - omega * zc * sin / cos
    };
    let (mut lo, mut hi) = (0.01 * isolated.min(pole), 0.95 * isolated.min(pole));
    assert!(characteristic(lo) > 0.0 && characteristic(hi) < 0.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if characteristic(mid) > 0.0 { lo = mid; } else { hi = mid; }
    }
    let exact = 0.5 * (lo + hi);
    let system = oscillator(isolated);
    let column = [1.0 / det::sqrt(inertance)];
    let force = C64::from_re(1e-8);
    let gate = CancelGate::new_clock_free();
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let duct = cylinder(elements, 2, radius, length);
        let stiffness = |omega| {
            let response = system.harmonic_response_with_duct(
                omega, &[force], &duct, PortId(7), &column, WORK, &gate,
            ).unwrap();
            (force / response.modal.displacement_m_sqrt_kg[0]).re
        };
        let (mut lo, mut hi) = (0.8 * exact, 1.2 * exact);
        assert!(stiffness(lo) > 0.0 && stiffness(hi) < 0.0);
        for _ in 0..28 {
            let mid = 0.5 * (lo + hi);
            if stiffness(mid) > 0.0 { lo = mid; } else { hi = mid; }
        }
        let numerical = 0.5 * (lo + hi);
        errors.push((1200.0 * det::ln(numerical / exact) / core::f64::consts::LN_2).abs());
    }
    assert!(errors[1] < errors[0], "loaded resonance errors in cents: {errors:?}");
    assert!(errors[2] < errors[1], "loaded resonance errors in cents: {errors:?}");
    assert!(errors[2] < 0.1, "loaded resonance errors in cents: {errors:?}");
}

#[test]
fn bare_lossless_duct_pole_is_a_finite_joint_antiresonance() {
    // G0/G2: one order-one element has L=[0,G;-G^T,0] with a 2x2 G.
    // The two eigenvalues of GG^T determine its bare poles analytically;
    // this reference never uses the harmonic LU or its assembled joint pencil.
    let duct = cylinder(1, 1, 0.003, 0.08);
    let view = duct.continuous_operator();
    assert_eq!(view.dimension(), 4);
    let l = view.generator_per_s();
    let g = [l[2], l[3], l[6], l[7]];
    let a = g[0] * g[0] + g[1] * g[1];
    let d = g[2] * g[2] + g[3] * g[3];
    let c = g[0] * g[2] + g[1] * g[3];
    let discriminant = det::sqrt((a - d) * (a - d) + 4.0 * c * c);
    let omega = det::sqrt(0.5 * (a + d - discriminant));
    let system = oscillator(0.6 * omega);
    let b = 0.02;
    let force = C64::new(0.002, -0.001);
    let response = system.harmonic_response_with_duct(
        omega, &[force], &duct, PortId(7), &[b], WORK,
        &CancelGate::new_clock_free(),
    ).unwrap();
    let free_scale = force.abs() / (omega * omega);
    assert!(response.modal.displacement_m_sqrt_kg[0].abs() < 1e-10 * free_scale);
    assert!((response.inlet_pressure_pa - force.scale(1.0 / b)).abs()
        < 1e-10 * force.abs() / b);
    assert!(response.inlet_flow_m3_s.abs() < 1e-10 * omega * b * free_scale);
    assert!(response.mean_stored_energy_j > 0.0);
    assert!(response.solve_relative_residual <= 1e-11);
    assert!(response.power_residual_w.abs() < 1e-18);
    assert_eq!(system.samples_rendered(), 0);
}
