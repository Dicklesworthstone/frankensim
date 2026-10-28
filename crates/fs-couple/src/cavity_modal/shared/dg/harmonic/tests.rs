use super::*;
use crate::cavity_modal::shared::dg::tests::{cavity, load};
use fs_duct::dg::PortId;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-9 * a.abs().max(b.abs()) + 1e-25,
        "{a:e} != {b:e}"
    );
}
fn complex_close(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 2e-9 * a.abs().max(b.abs()) + 1e-22,
        "{a:?} != {b:?}"
    );
}
fn create(rate: u32) -> SharedCavityDg {
    cavity(rate, true, 1e6)
        .with_downstream_duct(load(rate, true, 0.001), &CancelGate::new_clock_free())
        .unwrap()
}

#[test]
fn joint_cavity_duct_matches_physical_flow_pressure_storage_and_each_loss() {
    // Independent physical reduction in surface displacement and neck VOLUME
    // FLOW. Thermal compression uses the original pressure-diffusion law, not
    // production's inverted Maxwell spectrum. The isolated DG field supplies
    // an impedance only to this oracle; production retains every DG state.
    let coupled = create(48000);
    let source = coupled.cavity().source();
    let area = source.neck_area_m2();
    let cad = area * area / source.compression_stiffness_n_m();
    let neck_inertia = source.neck_mass_kg() / (area * area);
    let bore_inertia = source.gas().density * source.geometry().neck_length_m / area;
    let thermal = coupled
        .cavity()
        .cavity_thermal_compression()
        .unwrap()
        .diffusion();
    let viscous = coupled.cavity().neck_viscous_diffusion().unwrap();
    let a = 1e-4;
    let force = C64::new(0.001, -0.0003);
    let incident = C64::new(0.02, 0.01);
    let gate = CancelGate::new_clock_free();
    for omega in [700.0, 1400.0, 2000.0, 3500.0] {
        let s = C64::new(0.0, -omega);
        let compliance = thermal
            .shunt_terms()
            .iter()
            .fold(C64::from_re(1.0), |sum, t| {
                sum + C64::from_re(t.gain_per_s) / (s + C64::from_re(t.rate_per_s))
            })
            .scale(cad);
        let neck = s.scale(neck_inertia)
            + viscous
                .series_terms()
                .iter()
                .fold(C64::from_re(viscous.series_drag_per_s()), |sum, t| {
                    sum + (s / (s + C64::from_re(t.rate_per_s))).scale(t.gain_per_s)
                })
                .scale(bore_inertia);
        let unit = coupled
            .load()
            .duct()
            .harmonic_response(
                omega,
                &[(PortId(7), C64::from_re(1.0))],
                100_000_000,
                1e-10,
                &mut || false,
            )
            .unwrap();
        let impedance = unit.ports[0].pressure_pa;
        let aa =
            C64::from_re(1400.0_f64.powi(2) - omega * omega) + C64::from_re(a * a) / compliance;
        let ab = C64::from_re(a) / (s * compliance);
        let ba = C64::from_re(a) / compliance;
        let bb = neck + impedance + (s * compliance).recip();
        let determinant = aa * bb - ab * ba;
        let q = (force * bb + ab * incident) / determinant;
        let flow = (-aa * incident - ba * force) / determinant;
        let v = s * q;
        let pressure = -(v.scale(a) + flow) / (s * compliance);
        let duct_pressure = impedance * flow;
        let actual = coupled
            .harmonic_response(omega, &[force], incident, 100_000_000, &gate)
            .unwrap();
        complex_close(actual.cavity.network.displacement_m_sqrt_kg[0], q);
        complex_close(actual.cavity.neck_outflow_m3_s, flow);
        complex_close(actual.duct_inlet_flow_m3_s, flow);
        complex_close(actual.duct_inlet_pressure_pa, duct_pressure);
        complex_close(actual.cavity.pressure_increment_pa, pressure);
        complex_close(actual.cavity.volume_change_m3, -compliance * pressure);
        complex_close(
            s * actual.cavity.volume_change_m3,
            actual.cavity.surface_volume_rate_m3_s + flow,
        );
        for (&value, &reference) in actual.duct.pressure_pa.iter().zip(&unit.pressure_pa) {
            complex_close(value, reference * flow);
        }
        for (&value, &reference) in actual
            .duct
            .volume_flow_m3_s
            .iter()
            .zip(&unit.volume_flow_m3_s)
        {
            complex_close(value, reference * flow);
        }
        let mut availability = 0.25 * cad * pressure.norm_sq();
        let mut thermal_power = 0.0;
        for t in thermal.shunt_terms() {
            let temperature_pressure =
                pressure.scale(t.rate_per_s) / (s + C64::from_re(t.rate_per_s));
            availability +=
                0.25 * cad * t.gain_per_s / t.rate_per_s * temperature_pressure.norm_sq();
            thermal_power += 0.5 * cad * t.gain_per_s * (pressure - temperature_pressure).norm_sq();
        }
        let mut viscous_power = 0.5 * bore_inertia * viscous.series_drag_per_s() * flow.norm_sq();
        for t in viscous.series_terms() {
            viscous_power += 0.5
                * bore_inertia
                * t.gain_per_s
                * (s * flow / (s + C64::from_re(t.rate_per_s))).norm_sq();
        }
        close(actual.cavity.mean_compression_energy_j, availability);
        close(
            actual.cavity.mean_neck_kinetic_energy_j,
            0.25 * neck_inertia * flow.norm_sq(),
        );
        close(actual.cavity.network.material_dissipation_w, thermal_power);
        close(
            actual.cavity.network.connection_dissipation_w,
            viscous_power,
        );
        close(
            actual.duct.mean_stored_energy_j,
            unit.mean_stored_energy_j * flow.norm_sq(),
        );
        close(
            actual.duct.wall_dissipation_w,
            unit.wall_dissipation_w * flow.norm_sq(),
        );
        close(
            actual.duct.input_power_w,
            0.5 * (duct_pressure * flow.conj()).re,
        );
        close(actual.duct.input_power_w, actual.duct.wall_dissipation_w);
        close(
            actual.authored_external_power_w,
            0.5 * (force * v.conj() - incident * flow.conj()).re,
        );
        close(
            actual.authored_external_power_w,
            actual.total_dissipation_w(),
        );
        close(
            actual.mean_stored_energy_j,
            actual.cavity.network.mean_stored_energy_j + actual.duct.mean_stored_energy_j,
        );
        assert!(
            actual.interface_power_residual_w.abs()
                <= 2e-9 * actual.authored_external_power_w.abs() + 1e-25
        );
        assert!(
            actual.power_residual_w.abs() <= 2e-9 * actual.authored_external_power_w.abs() + 1e-25
        );
        assert!(actual.solve_relative_residual <= 1e-10);
    }
}

#[test]
fn continuous_joint_query_ignores_clock_and_transients_and_preserves_all_histories() {
    let mut coupled = create(48000);
    let other_clock = create(96000);
    let gate = CancelGate::new_clock_free();
    let force = C64::new(0.001, 0.0007);
    let incident = C64::new(0.02, -0.01);
    let before = coupled
        .harmonic_response(1400.0, &[force], incident, 100_000_000, &gate)
        .unwrap();
    assert_eq!(
        before,
        other_clock
            .harmonic_response(1400.0, &[force], incident, 100_000_000, &gate)
            .unwrap()
    );
    for _ in 0..8 {
        coupled.step(&[0.001], 0.02, &gate).unwrap();
    }
    assert_eq!(
        before,
        coupled
            .harmonic_response(1400.0, &[force], incident, 100_000_000, &gate)
            .unwrap()
    );
    let states: Vec<_> = coupled
        .cavity()
        .network()
        .components()
        .iter()
        .map(|m| m.states().to_vec())
        .collect();
    let memory = coupled.cavity().network().material_memory().to_vec();
    let duct = coupled.load().state().clone();
    let frame = coupled.last_frame().unwrap().clone();
    let observation = coupled.cavity().observe().unwrap();
    for omega in [0.0, 200.0, 10001.0, f64::NAN, f64::INFINITY] {
        assert!(
            coupled
                .harmonic_response(omega, &[force], incident, 100_000_000, &gate)
                .is_err()
        );
    }
    assert!(
        coupled
            .harmonic_response(1400.0, &[], incident, 100_000_000, &gate)
            .is_err()
    );
    assert!(
        coupled
            .harmonic_response(
                1400.0,
                &[force],
                C64::new(0.0, f64::NAN),
                100_000_000,
                &gate
            )
            .is_err()
    );
    assert!(
        coupled
            .harmonic_response(1400.0, &[force], incident, 1, &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        coupled.harmonic_response(1400.0, &[force], incident, 100_000_000, &cancelled),
        Err(SharedCavityDgError::Modal(ModalCouplingError::Cancelled))
    ));
    for (model, state) in coupled.cavity().network().components().iter().zip(states) {
        assert_eq!(model.states(), state);
    }
    assert_eq!(coupled.cavity().network().material_memory(), memory);
    assert_eq!(coupled.load().state(), &duct);
    assert_eq!(coupled.last_frame(), Some(&frame));
    assert_eq!(coupled.cavity().observe().unwrap(), observation);
    assert_eq!(coupled.cavity().network().samples_rendered(), 8);
    assert_eq!(coupled.load().state().step_count(), 8);
}
