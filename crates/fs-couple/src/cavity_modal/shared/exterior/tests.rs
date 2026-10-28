use super::*;
use crate::cavity_modal::{HelmholtzModalCavity, shared::CavitySurface};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::force::coupled::ModalCouplingConfig;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_scenario::{HelmholtzCavity, NeckEnd};

fn cells(sign: f64) -> Vec<AreaSample> {
    vec![
        AreaSample { position_m: [-0.025, 0.0], area_m2: 0.0002, mode_shape: vec![100.0 * sign] },
        AreaSample { position_m: [0.025, 0.0], area_m2: 0.0002, mode_shape: vec![-40.0 * sign] },
    ]
}
fn cavity(rate: u32, sign: f64) -> SharedCavity {
    let source = HelmholtzModalCavity::new(HelmholtzCavity {
        volume_m3: 0.002, neck_radius_m: 0.01, neck_length_m: 0.02,
        inner_end: NeckEnd::Flanged, outer_end: NeckEnd::Flanged,
    }, GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
        CavityNeckLoss::ExteriorMouth).unwrap();
    let budget = ModalAcousticTimeBudget::audible_reference();
    let surface = ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: 700.0, damping_ratio: 0.0, pressure_per_modal_velocity: C64::ZERO,
    }], budget).unwrap();
    SharedCavity::new(source, vec![CavitySurface {
        model: surface, outward_volume_shapes: vec![0.012 * sign],
    }], rate, ModalAcousticState { displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: 0.001 },
        budget, ModalCouplingConfig {
            max_modes: 8, max_connections: 4, max_setup_terms: 4096, nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1.0, maximum_abs_pressure_pa: 1e5,
            maximum_abs_connection_force_n: 1000.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8,
        }, &CancelGate::new_clock_free()).unwrap()
}
fn mouth(rings: usize) -> CircularBaffledMouth {
    CircularBaffledMouth { centre_m: [0.0, 0.04], radial_rings: rings, first_ring_sectors: 8 }
}
fn options() -> BaffledExteriorOptions {
    BaffledExteriorOptions { band_rad_s: [200.0, 1800.0], frequency_samples: 65,
        maximum_cells: 1024, maximum_sweep_terms: 200_000_000, maximum_poles: 64 }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-13 + 2e-9 * a.abs().max(b.abs()), "{a:e} != {b:e}");
}

#[test]
fn actual_circular_mouth_quadrature_converges_to_flanged_inertance_and_power() {
    // G1/G3: independent low-ka circular-piston limits, not another evaluation
    // of the discrete Rayleigh self-cell formula. Three geometric rungs retain
    // the original mouth radius, area, gas and mass-normalized neck coordinate.
    let cavity = cavity(192000, 0.0);
    let gas = cavity.source.gas();
    let radius = cavity.source.geometry().neck_radius_m;
    let projection = cavity.source.neck_volume_shape_m2_per_sqrt_kg();
    let omega = 100.0;
    let exact_inertance = 8.0 * gas.density / (3.0 * core::f64::consts::PI.powi(2) * radius);
    let compact_resistance = gas.density * omega * omega
        / (2.0 * core::f64::consts::PI * gas.sound_speed);
    let mut previous_error = f64::INFINITY;
    for rings in [2, 4, 8] {
        let (samples, _) = cavity.baffled_samples(&cells(0.0), mouth(rings), options(),
            &CancelGate::new_clock_free()).unwrap();
        let mouth_area: f64 = samples[2..].iter().map(|cell| cell.area_m2).sum();
        close(mouth_area, cavity.source.neck_area_m2());
        let surface = PreparedSurface::new(&samples).unwrap();
        let z = surface.power_consistent_impedance(omega, Medium {
            density: gas.density, sound_speed: gas.sound_speed,
        }).unwrap()[3].scale(1.0 / (projection * projection));
        let error = (-z.im / (omega * exact_inertance) - 1.0).abs();
        assert!(error < previous_error, "rings={rings}: error={error} previous={previous_error}");
        previous_error = error;
        assert!((z.re / compact_resistance - 1.0).abs() < 2e-5);
    }
    assert!(previous_error < 0.008);
}

#[test]
fn geometry_fit_retains_signed_mutual_load_and_passive_joint_motion() {
    let gate = CancelGate::new_clock_free();
    let original = cavity(192000, 1.0);
    let initial = original.network().total_energy_j().unwrap();
    let (mut whole, report) = original.with_baffled_exterior(&cells(1.0), mouth(2), options(), &gate).unwrap();
    assert_eq!(report.physical_fields, 2);
    assert_eq!(report.surface_cells, 34);
    assert!(report.peak_complex_error <= 0.05 && report.rms_complex_error <= 0.02);
    assert!(report.peak_resistance_error <= 0.10 && report.rms_resistance_error <= 0.05);
    // The accepted model is in original modal coordinates. Pull its neck column
    // back to physical Q for the second clock and the independent power check.
    let mut model = whole.network().radiation_model().unwrap().clone();
    let projection = whole.source.neck_volume_shape_m2_per_sqrt_kg();
    for pole in &mut model.poles { pole.coupling[1] /= projection; }
    let mut half = cavity(384000, 1.0).with_exterior_radiation(&model, &[0, 1], 64, &gate).unwrap();
    let (opposite, _) = cavity(192000, -1.0).prepare_baffled_exterior(
        &cells(-1.0), mouth(2), options(), &gate).unwrap();
    let omega = 975.0; // held-out odd frequency, never used by the fit.
    let z = model.impedance(omega).unwrap();
    let flipped = opposite.impedance(omega).unwrap();
    close(z[0].re, flipped[0].re); close(z[0].im, flipped[0].im);
    close(z[3].re, flipped[3].re); close(z[3].im, flipped[3].im);
    close(z[1].re, -flipped[1].re); close(z[1].im, -flipped[1].im);
    assert!(z[1].abs() > 0.0);
    // Independent outgoing-pressure Gram expression for the mutual resistance.
    // This catches a missing mouth-area factor or a lost signed body trace.
    let (samples, _) = whole.baffled_samples(&cells(1.0), mouth(2), options(), &gate).unwrap();
    let gas = whole.source.gas();
    let mut mutual = 0.0;
    for body in &samples[..2] {
        for neck in &samples[2..] {
            let r = (body.position_m[0] - neck.position_m[0]).hypot(body.position_m[1] - neck.position_m[1]);
            let kr = omega * r / gas.sound_speed;
            mutual += body.area_m2 * body.mode_shape[0]
                * neck.area_m2 / whole.source.neck_area_m2() * det::sin(kr) / kr;
        }
    }
    mutual *= gas.density * omega * omega / (2.0 * core::f64::consts::PI * gas.sound_speed);
    assert!((z[1].re / mutual - 1.0).abs() < 0.10, "fit={} physical={mutual}", z[1].re);
    let mut loss = 0.0;
    for _ in 0..8 {
        let a = whole.step(&[0.0], 0.0, &gate).unwrap().clone();
        let b = half.step(&[0.0], 0.0, &gate).unwrap().clone();
        let c = half.step(&[0.0], 0.0, &gate).unwrap().clone();
        loss += a.total_dissipation_j();
        close(a.total_dissipation_j(), b.total_dissipation_j() + c.total_dissipation_j());
        assert!(a.total_dissipation_j() >= 0.0);
        assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
        for (a, b) in whole.network().components().iter().zip(half.network().components()) {
            for (a, b) in a.states().iter().zip(b.states()) {
                close(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
                close(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
            }
        }
        close(whole.neck_exterior_pressure_pa().unwrap(), half.neck_exterior_pressure_pa().unwrap());
    }
    close(whole.network().total_energy_j().unwrap() + loss, initial);
    assert!(whole.neck_exterior_pressure_pa().unwrap().abs() > 0.0);
}

#[test]
fn geometry_attachment_refuses_changed_basis_double_loading_and_unresolved_work() {
    let gate = CancelGate::new_clock_free();
    let mut cavity = cavity(192000, 1.0);
    let before: Vec<_> = cavity.network().components().iter().map(|m| m.states().to_vec()).collect();
    assert!(cavity.prepare_baffled_exterior(&cells(2.0), mouth(2), options(), &gate).is_err());
    assert!(cavity.prepare_baffled_exterior(&cells(1.0), mouth(2),
        BaffledExteriorOptions { maximum_cells: 33, ..options() }, &gate).is_err());
    assert!(cavity.prepare_baffled_exterior(&cells(1.0), mouth(2),
        BaffledExteriorOptions { maximum_sweep_terms: 1, ..options() }, &gate).is_err());
    assert!(cavity.prepare_baffled_exterior(&cells(1.0), mouth(2),
        BaffledExteriorOptions { band_rad_s: [200.0, 500.0], ..options() }, &gate).is_err());
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(cavity.prepare_baffled_exterior(&cells(1.0), mouth(2), options(), &cancelled),
        Err(BaffledExteriorError::Cancelled)));
    for (model, expected) in cavity.network().components().iter().zip(before) { assert_eq!(model.states(), expected); }
    cavity.source.loss = CavityNeckLoss::RadiationAtResonance;
    assert!(cavity.prepare_baffled_exterior(&cells(1.0), mouth(2), options(), &gate).is_err());
    cavity.source.loss = CavityNeckLoss::ExteriorMouth;
    cavity.source.geometry.outer_end = NeckEnd::Unflanged;
    assert!(cavity.prepare_baffled_exterior(&cells(1.0), mouth(2), options(), &gate).is_err());
    cavity.source.geometry.outer_end = NeckEnd::Flanged;
    cavity.step(&[0.0], 0.0, &gate).unwrap();
    assert!(cavity.prepare_baffled_exterior(&cells(1.0), mouth(2), options(), &gate).is_err());
}
