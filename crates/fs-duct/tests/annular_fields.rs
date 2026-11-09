//! G0/G1/G3: original annular fields, physical heat and independent limits.
use fs_duct::DuctError;
use fs_duct::fluid::FluidDuctError;
use fs_duct::fluid::annular::{AnnularAdmission, AnnularSection};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use std::cell::Cell;
use std::f64::consts::PI;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn admission() -> AnnularAdmission {
    AnnularAdmission {
        maximum_transverse_phase: 0.8,
        maximum_cell_diffusion_phase: 0.75,
        maximum_relative_residual: 1e-8,
        maximum_payload_bytes: 2 << 20,
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?} at {tolerance}"
    );
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:.17e} != {expected:.17e} at {tolerance}"
    );
}

#[test]
fn no_slip_velocity_and_temperature_recover_independent_annular_poiseuille_profiles() {
    let gas = air();
    let (a, b, omega) = (0.0004_f64, 0.0012_f64, 0.00001);
    let pressure = C64::new(3.0, -2.0);
    let flow = C64::new(2e-8, 5e-9);
    let geometry = b.powi(4) - a.powi(4) - (b * b - a * a).powi(2) / (b / a).ln();
    let mut previous = f64::INFINITY;
    for cells in [32, 64, 128, 256] {
        let section = AnnularSection::new(a, b, cells).unwrap();
        let prepared = section.prepare_fields(&gas, omega, admission()).unwrap();
        let wave = section.evaluate(&gas, omega, admission()).unwrap();
        assert_eq!(prepared.wave().wave.wavenumber, wave.wave.wavenumber);
        assert_eq!(
            prepared.wave().wave.characteristic_impedance,
            wave.wave.characteristic_impedance
        );
        assert_eq!(prepared.wave().viscous.mean, wave.viscous.mean);
        assert_eq!(
            prepared.wave().thermal.unwrap().mean,
            wave.thermal.unwrap().mean
        );
        assert_eq!(prepared.angular_frequency_rad_s(), omega);
        let field = prepared
            .field_in_cell(cells / 2, 0.0, pressure, flow)
            .unwrap();
        let r = field.radius_m;
        let radial = b * b - r * r - (b * b - a * a) * (b / r).ln() / (b / a).ln();
        let expected_velocity = flow.scale(2.0 * radial / (PI * geometry));
        // The steady annular Poiseuille solution and diffusion Poisson limit
        // use physical radii and logarithms, not the production P1 matrices.
        let expected_temperature =
            C64::new(0.0, -omega * radial / (4.0 * gas.thermal_conductivity)) * pressure;
        near(field.axial_velocity_m_s, expected_velocity, 0.001);
        near(
            prepared
                .gas_temperature_in_cell(cells / 2, 0.0, pressure)
                .unwrap(),
            expected_temperature,
            0.001,
        );
        let error = (field.axial_velocity_m_s - expected_velocity).abs();
        if previous.is_finite() {
            assert!(
                (3.8..4.2).contains(&(previous / error)),
                "N={cells}, ratio={}",
                previous / error
            );
        }
        previous = error;
        for (cell, fraction) in [(0, 0.0), (cells - 1, 1.0)] {
            let wall = prepared
                .field_in_cell(cell, fraction, pressure, flow)
                .unwrap();
            assert_eq!(wall.axial_velocity_m_s, C64::ZERO);
            assert_eq!(wall.thermal_pressure_response_pa, Some(C64::ZERO));
            assert_eq!(
                prepared
                    .gas_temperature_in_cell(cell, fraction, pressure)
                    .unwrap(),
                C64::ZERO
            );
        }
    }
}

#[test]
fn radial_fields_integrate_to_the_original_flow_and_independent_passive_work() {
    let gas = air();
    let section = AnnularSection::new(0.0004, 0.0012, 256).unwrap();
    let pressure = C64::new(30.0, 7.0);
    let flow = C64::new(4e-8, -3e-8);
    for omega in [10.0, 1000.0, 50000.0] {
        let prepared = section.prepare_fields(&gas, omega, admission()).unwrap();
        let mut integrated_flow = C64::ZERO;
        let mut area = 0.0;
        let mut heat = [0.0; 2];
        for cell in 0..section.elements() {
            let r0 = section.radius_m(cell).unwrap();
            let r1 = section.radius_m(cell + 1).unwrap();
            relative(
                prepared.cell_area_m2(cell).unwrap(),
                PI * (r1 - r0) * (r1 + r0),
                2e-13,
            );
            area += prepared.cell_area_m2(cell).unwrap();
            // Independent two-point Gauss quadrature of 2*pi*r*u(r); exact
            // for the P1 velocity times the original physical radius.
            for fraction in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
                let field = prepared
                    .field_in_cell(cell, fraction, pressure, flow)
                    .unwrap();
                integrated_flow = integrated_flow
                    + field
                        .axial_velocity_m_s
                        .scale(PI * field.radius_m * (r1 - r0));
                relative(
                    field.viscous_dissipation_w_m3,
                    0.5 * gas.dynamic_viscosity * field.radial_velocity_gradient_per_s.norm_sq(),
                    3e-14,
                );
            }
            let field = prepared.field_in_cell(cell, 0.5, pressure, flow).unwrap();
            let cell_heat = prepared
                .cell_dissipation_per_axial_length(cell, pressure, flow)
                .unwrap();
            relative(
                cell_heat.viscous_w_m,
                field.viscous_dissipation_w_m3 * prepared.cell_area_m2(cell).unwrap(),
                2e-14,
            );
            relative(
                cell_heat.thermal_w_m,
                field.thermal_dissipation_w_m3 * prepared.cell_area_m2(cell).unwrap(),
                2e-14,
            );
            heat[0] += cell_heat.viscous_w_m;
            heat[1] += cell_heat.thermal_w_m;
        }
        near(integrated_flow, flow, 2e-13);
        relative(area, section.area_m2(), 2e-14);
        let total = prepared
            .dissipation_per_axial_length(pressure, flow)
            .unwrap();
        relative(total.viscous_w_m, heat[0], 2e-14);
        relative(total.thermal_w_m, heat[1], 2e-14);
        let wave = prepared.wave().wave;
        // Compare positive ORIGINAL radial work with the independently formed
        // pressure/flow telegraph work; never use that difference as heat.
        relative(
            total.viscous_w_m,
            0.5 * (wave.wavenumber * wave.characteristic_impedance).im * flow.norm_sq(),
            2e-8,
        );
        relative(
            total.thermal_w_m,
            0.5 * (wave.wavenumber / wave.characteristic_impedance).im * pressure.norm_sq(),
            2e-8,
        );
        assert_eq!(
            prepared
                .dissipation_per_axial_length(C64::ZERO, flow)
                .unwrap()
                .thermal_w_m,
            0.0
        );
        assert_eq!(
            prepared
                .dissipation_per_axial_length(pressure, C64::ZERO)
                .unwrap()
                .viscous_w_m,
            0.0
        );
    }
}

fn travel(k: C64, x: f64) -> C64 {
    C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp())
}

#[test]
fn coherent_counterwaves_lose_exactly_the_observed_local_irreversible_heat() {
    let gas = air();
    let section = AnnularSection::new(0.001, 0.002, 512).unwrap();
    let prepared = section.prepare_fields(&gas, 4000.0, admission()).unwrap();
    let wave = prepared.wave().wave;
    let fields = |x: f64| {
        let forward = C64::new(20.0, 3.0) * travel(wave.wavenumber, x);
        let backward = C64::new(-6.0, 7.0) * travel(wave.wavenumber, 0.2 - x);
        (
            forward + backward,
            (forward - backward) / wave.characteristic_impedance,
        )
    };
    let power = |x: f64| {
        let (p, u) = fields(x);
        0.5 * (p * u.conj()).re
    };
    for x in [0.03, 0.07, 0.13, 0.18] {
        let (p, u) = fields(x);
        let heat = prepared.dissipation_per_axial_length(p, u).unwrap();
        let h = 1e-6;
        let power_gradient = (power(x - h) - power(x + h)) / (2.0 * h);
        relative(heat.total_w_m, power_gradient, 3e-8);
        assert!(heat.viscous_w_m > 0.0 && heat.thermal_w_m > 0.0);
    }
}

#[test]
fn liquid_without_thermal_expansion_has_only_its_original_shear_loss() {
    let state = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap(),
        4200.0,
    )
    .unwrap();
    let section = AnnularSection::new(0.001, 0.002, 256).unwrap();
    let prepared = section.prepare_fields(&state, 1000.0, admission()).unwrap();
    let pressure = C64::new(100.0, 20.0);
    let flow = C64::new(2e-9, -7e-9);
    let field = prepared.field_in_cell(20, 0.25, pressure, flow).unwrap();
    assert!(std::ptr::eq(prepared.source(), &state));
    assert_eq!(*prepared.state(), state);
    assert_eq!(field.thermal_pressure_response_pa, None);
    assert_eq!(field.thermal_dissipation_w_m3, 0.0);
    let heat = prepared
        .dissipation_per_axial_length(pressure, flow)
        .unwrap();
    assert_eq!(heat.thermal_w_m, 0.0);
    assert!(heat.viscous_w_m > 0.0);
}

struct OneSnapshot {
    state: ThermoviscousFluidState,
    queries: Cell<usize>,
}
impl FluidAcoustics for OneSnapshot {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        panic!("query the complete state once")
    }
}
impl FluidTransport for OneSnapshot {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("query the complete state once")
    }
}
impl ThermoviscousFluid for OneSnapshot {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn original_source_once_and_transactional_admission_remain_bound_to_every_field() {
    let source = OneSnapshot {
        state: air().thermoviscous_state().unwrap(),
        queries: Cell::new(0),
    };
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    assert!(matches!(
        section.prepare_fields_controlled(&source, 1000.0, admission(), || true),
        Err(FluidDuctError::Duct(DuctError::Cancelled {
            completed_steps: 0
        }))
    ));
    assert_eq!(source.queries.get(), 0);
    let prepared = section
        .prepare_fields(&source, 1000.0, admission())
        .unwrap();
    assert!(std::ptr::eq(prepared.source(), &source));
    assert!(std::ptr::eq(prepared.section(), &section));
    for _ in 0..10 {
        prepared
            .field_in_cell(14, 0.4, C64::ONE, C64::from_re(1e-8))
            .unwrap();
        prepared
            .dissipation_per_axial_length(C64::ONE, C64::from_re(1e-8))
            .unwrap();
    }
    assert_eq!(source.queries.get(), 1);
    for (cell, fraction, p, q) in [
        (128, 0.5, C64::ONE, C64::ONE),
        (0, -0.1, C64::ONE, C64::ONE),
        (0, 1.1, C64::ONE, C64::ONE),
        (0, f64::NAN, C64::ONE, C64::ONE),
        (0, 0.5, C64::new(f64::NAN, 0.0), C64::ONE),
        (0, 0.5, C64::ONE, C64::new(0.0, f64::INFINITY)),
    ] {
        assert!(prepared.field_in_cell(cell, fraction, p, q).is_err());
    }
    assert!(
        prepared
            .cell_dissipation_per_axial_length(128, C64::ONE, C64::ONE)
            .is_err()
    );
    assert!(
        prepared
            .dissipation_per_axial_length(C64::new(f64::INFINITY, 0.0), C64::ONE)
            .is_err()
    );
    let zero = prepared
        .field_in_cell(0, 0.5, C64::ZERO, C64::ZERO)
        .unwrap();
    assert_eq!(zero.axial_velocity_m_s, C64::ZERO);
    assert_eq!(zero.viscous_dissipation_w_m3, 0.0);
    assert_eq!(zero.thermal_dissipation_w_m3, 0.0);
    let mut polls = 0;
    assert!(matches!(
        section.prepare_fields_controlled(&source, 1000.0, admission(), || {
            polls += 1;
            polls == 6
        }),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert_eq!(polls, 6);
    let mut small = admission();
    small.maximum_payload_bytes = 192 * (section.elements() + 1) + 4095;
    let queries = source.queries.get();
    assert!(section.prepare_fields(&source, 1000.0, small).is_err());
    assert_eq!(source.queries.get(), queries);
}
