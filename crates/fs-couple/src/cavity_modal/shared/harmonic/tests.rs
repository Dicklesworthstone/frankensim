use super::*;
use crate::render::plate::impact::radiation::Pole;
use fs_material::gas::GasSpec;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 8, max_connections: 6, max_setup_terms: 20000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5, maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-14,
        energy_relative_tolerance: 1e-8 }
}

fn create(sign: f64) -> (SharedCavity, RadiationModel) {
    let gate = CancelGate::new_clock_free();
    let source = HelmholtzModalCavity::new(HelmholtzCavity {
        volume_m3: 1e-7, neck_radius_m: 1e-4, neck_length_m: 0.1,
        inner_end: fs_scenario::NeckEnd::Flanged, outer_end: fs_scenario::NeckEnd::Flanged,
    }, GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
        CavityNeckLoss::ExteriorMouth).unwrap();
    let surface = CavitySurface { model: ModalAcousticTimeModel::try_new(48000,
        vec![ModalAcousticMode { angular_frequency_rad_s: 500.0, damping_ratio: 0.02,
            pressure_per_modal_velocity: C64::ZERO }], ModalAcousticTimeBudget::audible_reference()).unwrap(),
        outward_volume_shapes: vec![sign * 1e-4] };
    let cavity = SharedCavity::new_with_cylindrical_thermal_wall(source, vec![surface], 48000,
        ModalAcousticState { displacement_m_sqrt_kg: 2e-7, velocity_m_sqrt_kg_per_s: 0.0 },
        ModalAcousticTimeBudget::audible_reference(), config(),
        Some(ViscousNeckOptions { poles: 4, band_rad_s: [10.0, 1000.0], maximum_relative_error: 1e-3 }),
        CylindricalCavityThermalOptions { radius_m: 0.001, poles: 4,
            band_rad_s: [10.0, 1000.0], maximum_relative_error: 1e-3 },
        Some(vec![0.0; 4]), &gate).unwrap();
    let radiation = RadiationModel { ports: 2, poles: vec![Pole {
        omega: 700.0, zeta: 0.12, coupling: vec![sign * 80.0, -3e5],
    }] };
    (cavity.with_exterior_radiation(&radiation, &[0, 1], 1, &gate).unwrap(), radiation)
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-9 * a.abs().max(b.abs()) + 1e-25, "{a:e} != {b:e}");
}
fn complex_close(a: C64, b: C64) {
    assert!((a - b).abs() <= 2e-9 * a.abs().max(b.abs()) + 1e-22, "{a:?} != {b:?}");
}

#[test]
fn loaded_harmonic_cavity_matches_physical_pressure_flow_and_independent_power() {
    // Independent physical two-variable reduction in surface displacement and
    // neck VOLUME FLOW, not the production modal/Maxwell/radiation matrix. The
    // thermal reference retains the original pressure-diffusion compliance;
    // production instead advances its inverted positive stiffness spectrum.
    let (cavity, radiation) = create(1.0);
    let source = cavity.source();
    let area = source.neck_area_m2();
    let cad = area * area / source.compression_stiffness_n_m();
    let neck_inertia = source.neck_mass_kg() / (area * area);
    let bore_inertia = source.gas().density * source.geometry().neck_length_m / area;
    let thermal = cavity.cavity_thermal_compression().unwrap().diffusion();
    let viscous = cavity.neck_viscous_diffusion().unwrap();
    let a = 1e-4;
    let surface_force = C64::new(0.001, -0.0003);
    let incident = C64::new(0.02, 0.01);
    let gate = CancelGate::new_clock_free();
    for omega in [10.0, 40.0, 300.0, 500.0, 700.0, 1000.0] {
        let s = C64::new(0.0, -omega);
        let compliance = thermal.shunt_terms().iter().fold(C64::from_re(1.0), |sum, term|
            sum + C64::from_re(term.gain_per_s) / (s + C64::from_re(term.rate_per_s))).scale(cad);
        let neck = s.scale(neck_inertia)
            + viscous.series_terms().iter().fold(C64::from_re(viscous.series_drag_per_s()), |sum, term|
                sum + (s / (s + C64::from_re(term.rate_per_s))).scale(term.gain_per_s)).scale(bore_inertia);
        let z = radiation.impedance(omega).unwrap();
        let d = C64::new(500.0_f64.powi(2) - omega * omega, -2.0 * 0.02 * 500.0 * omega);
        let aa = d + C64::from_re(a * a) / compliance + s * z[0];
        let ab = C64::from_re(a) / (s * compliance) + z[1];
        let ba = C64::from_re(a) / compliance + s * z[2];
        let bb = neck + (s * compliance).recip() + z[3];
        let determinant = aa * bb - ab * ba;
        let q = (surface_force * bb + ab * incident) / determinant;
        let flow = (-aa * incident - ba * surface_force) / determinant;
        let v = s * q;
        let pressure = -(v.scale(a) + flow) / (s * compliance);
        let exterior = z[2] * v + z[3] * flow;
        let actual = cavity.harmonic_response(omega, &[surface_force], incident, 100_000, &gate).unwrap();
        complex_close(actual.network.displacement_m_sqrt_kg[0], q);
        complex_close(actual.neck_outflow_m3_s, flow);
        complex_close(actual.pressure_increment_pa, pressure);
        complex_close(actual.neck_exterior_pressure_pa, exterior);
        complex_close(actual.volume_change_m3, -compliance * pressure);
        complex_close(actual.surface_volume_rate_m3_s, v.scale(a));
        complex_close(s * actual.volume_change_m3,
            actual.surface_volume_rate_m3_s + actual.neck_outflow_m3_s);
        // The original physical pressure diffusion has availability and heat
        // different in form from the transformed Maxwell coordinates.
        let mut availability = 0.25 * cad * pressure.norm_sq();
        let mut thermal_power = 0.0;
        for term in thermal.shunt_terms() {
            let temperature_pressure = pressure.scale(term.rate_per_s) / (s + C64::from_re(term.rate_per_s));
            availability += 0.25 * cad * term.gain_per_s / term.rate_per_s * temperature_pressure.norm_sq();
            thermal_power += 0.5 * cad * term.gain_per_s * (pressure - temperature_pressure).norm_sq();
        }
        let mut viscous_power = 0.5 * bore_inertia * viscous.series_drag_per_s() * flow.norm_sq();
        for term in viscous.series_terms() {
            let slip = s * flow / (s + C64::from_re(term.rate_per_s));
            viscous_power += 0.5 * bore_inertia * term.gain_per_s * slip.norm_sq();
        }
        let radiation_power = 0.5 * ((z[0] * v + z[1] * flow) * v.conj() + exterior * flow.conj()).re;
        close(actual.mean_compression_energy_j, availability);
        close(actual.mean_neck_kinetic_energy_j, 0.25 * neck_inertia * flow.norm_sq());
        close(actual.network.material_dissipation_w, thermal_power);
        close(actual.network.connection_dissipation_w, viscous_power);
        close(actual.network.radiation_dissipation_w, radiation_power);
        close(actual.network.component_dissipation_w, 0.5 * 2.0 * 0.02 * 500.0 * v.norm_sq());
        close(actual.surface_input_power_w, 0.5 * (surface_force * v.conj()).re);
        close(actual.incident_pressure_input_power_w, -0.5 * (incident * flow.conj()).re);
        close(actual.network.external_power_w, actual.surface_input_power_w + actual.incident_pressure_input_power_w);
        close(actual.network.external_power_w, actual.network.total_dissipation_w());
        assert!(actual.network.solve_relative_residual <= config().solve_relative_tolerance);
    }
}

#[test]
fn harmonic_response_preserves_signed_basis_and_ignores_accepted_transients() {
    let (mut forward, _) = create(1.0);
    let (reverse, _) = create(-1.0);
    let gate = CancelGate::new_clock_free();
    let force = C64::new(0.002, 0.0007);
    let incident = C64::new(0.03, -0.01);
    let before = forward.harmonic_response(480.0, &[force], incident, 100_000, &gate).unwrap();
    let signed = reverse.harmonic_response(480.0, &[-force], incident, 100_000, &gate).unwrap();
    complex_close(before.pressure_increment_pa, signed.pressure_increment_pa);
    complex_close(before.neck_exterior_pressure_pa, signed.neck_exterior_pressure_pa);
    complex_close(before.neck_outflow_m3_s, signed.neck_outflow_m3_s);
    complex_close(before.network.displacement_m_sqrt_kg[0], -signed.network.displacement_m_sqrt_kg[0]);
    close(before.network.external_power_w, signed.network.external_power_w);
    for _ in 0..8 { forward.step(&[0.001], 0.02, &gate).unwrap(); }
    assert_eq!(before, forward.harmonic_response(480.0, &[force], incident, 100_000, &gate).unwrap());
    assert_eq!(forward.network().samples_rendered(), 8);
}

#[test]
fn harmonic_cavity_band_budget_and_cancel_refusals_leave_histories_unchanged() {
    let (mut cavity, _) = create(1.0);
    let gate = CancelGate::new_clock_free();
    cavity.step(&[0.001], 0.02, &gate).unwrap();
    let states: Vec<_> = cavity.network().components().iter().map(|m| m.states().to_vec()).collect();
    let memory = cavity.network().material_memory().to_vec();
    let radiation = cavity.network().radiation_state().to_vec();
    let observation = cavity.observe().unwrap();
    let frame = cavity.network().last_frame().unwrap().clone();
    for omega in [0.0, 1.0, 1001.0, f64::NAN, f64::INFINITY] {
        assert!(cavity.harmonic_response(omega, &[C64::ZERO], C64::ZERO, 100_000, &gate).is_err());
    }
    assert!(cavity.harmonic_response(500.0, &[], C64::ZERO, 100_000, &gate).is_err());
    assert!(cavity.harmonic_response(500.0, &[C64::ZERO], C64::new(0.0, f64::NAN), 100_000, &gate).is_err());
    assert!(cavity.harmonic_response(500.0, &[C64::ZERO], C64::ZERO, 1, &gate).is_err());
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(cavity.harmonic_response(500.0, &[C64::ZERO], C64::ZERO, 100_000, &cancelled),
        Err(ModalCouplingError::Cancelled)));
    assert_eq!(cavity.network().material_memory(), memory);
    assert_eq!(cavity.network().radiation_state(), radiation);
    assert_eq!(cavity.observe().unwrap(), observation);
    assert_eq!(cavity.network().last_frame().unwrap(), &frame);
    for (m, state) in cavity.network().components().iter().zip(states) { assert_eq!(m.states(), state); }
}
