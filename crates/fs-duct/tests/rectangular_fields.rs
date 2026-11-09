//! G1/G3: original rectangular fields, independent transverse physics and work.
use fs_duct::DuctError;
use fs_duct::fluid::FluidDuctError;
use fs_duct::fluid::rectangular::{RectangularAdmission, RectangularSection};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::c64::C64;
use std::cell::Cell;
use std::f64::consts::PI;

// Numerical fixtures exercise the constitutive model; no measured-air or
// measured-liquid accuracy is inferred from these independently authored data.
fn fluid() -> ThermoviscousFluidState {
    ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        FluidTransportState::try_new(1.2, 1005.0, 0.026, 1.8e-5).unwrap(),
        718.0,
    ).unwrap()
}
fn admission() -> RectangularAdmission {
    RectangularAdmission {
        maximum_relative_tail: 1e-3,
        maximum_transverse_phase: 0.8,
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?} at {tolerance}");
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:.17e} != {expected:.17e} at {tolerance}");
}

#[test]
fn no_slip_and_isothermal_fields_recover_the_independent_slit_poiseuille_limit() {
    let material = fluid();
    let (a, b, omega) = (0.0003, 0.03, 1e-5);
    let section = RectangularSection::new(a, b, 512).unwrap();
    let prepared = section.prepare_fields(&material, omega, admission()).unwrap();
    let pressure = C64::new(30.0, -11.0);
    let flow = C64::new(4e-9, 7e-9);
    let sidewall_sum: f64 = (0..1024).map(|i| {
        let m = (2 * i + 1) as f64;
        (m * PI * b / (2.0 * a)).tanh() / m.powi(5)
    }).sum();
    let flow_correction = 1.0 - 192.0 * a * sidewall_sum / (PI.powi(5) * b);
    let thermal_diffusion = material.transport().thermal_diffusivity_m2_s();
    for fraction in [0.1, 0.25, 0.5, 0.8] {
        let y = a * fraction;
        let point = prepared.field_at(y, 0.5 * b, pressure, flow).unwrap();
        // At the long-side midpoint its distant end-wall correction is
        // exponentially tiny. These are the direct slit Poisson solutions.
        let velocity = flow.scale(6.0 * fraction * (1.0 - fraction) / (a * b * flow_correction));
        let thermal = pressure * C64::new(0.0, -omega * y * (a - y) / (2.0 * thermal_diffusion));
        near(point.axial_velocity_m_s, velocity, 3e-8);
        near(point.thermal_pressure_response_pa.unwrap(), thermal, 3e-8);
        if fraction != 0.5 {
            near(point.transverse_velocity_gradient_per_s[0],
                flow.scale(6.0 * (1.0 - 2.0 * fraction) / (a * a * b * flow_correction)), 3e-6);
        }
    }
    for (y, z) in [(0.0, 0.2 * b), (a, 0.7 * b), (0.3 * a, 0.0), (0.8 * a, b)] {
        let wall = prepared.field_at(y, z, pressure, flow).unwrap();
        assert_eq!(wall.axial_velocity_m_s, C64::ZERO);
        assert_eq!(wall.thermal_pressure_response_pa, Some(C64::ZERO));
    }
}

#[test]
fn original_fields_recover_volume_flow_and_both_positive_gradient_integrals() {
    let material = fluid();
    let (a, b, omega) = (0.001, 0.0015, 800.0);
    let section = RectangularSection::new(a, b, 24).unwrap();
    let prepared = section.prepare_fields(&material, omega, admission()).unwrap();
    let mean = section.evaluate(&material, omega, admission()).unwrap();
    assert_eq!(prepared.wave().wave.wavenumber, mean.wave.wavenumber);
    assert_eq!(prepared.wave().wave.characteristic_impedance, mean.wave.characteristic_impedance);
    assert_eq!(prepared.wave().viscous.mean, mean.viscous.mean);
    assert_eq!(prepared.wave().thermal.unwrap().mean, mean.thermal.unwrap().mean);
    assert!(std::ptr::eq(prepared.section(), &section));
    assert_eq!(prepared.angular_frequency_rad_s(), omega);
    let (pressure, flow) = (C64::new(30.0, 7.0), C64::new(4e-8, -3e-8));
    let cells = 48;
    let points = [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()];
    let weight = a * b / (4.0 * (cells * cells) as f64);
    let mut integrated_flow = C64::ZERO;
    let mut heat = [0.0; 2];
    for i in 0..cells {
        for j in 0..cells {
            for fi in points {
                for fj in points {
                    let field = prepared.field_at(
                        a * (i as f64 + fi) / cells as f64,
                        b * (j as f64 + fj) / cells as f64,
                        pressure, flow,
                    ).unwrap();
                    integrated_flow = integrated_flow + field.axial_velocity_m_s.scale(weight);
                    heat[0] += weight * field.viscous_dissipation_w_m3;
                    heat[1] += weight * field.thermal_dissipation_w_m3;
                    relative(field.viscous_dissipation_w_m3,
                        0.5 * material.transport().dynamic_viscosity_pa_s()
                            * field.transverse_velocity_gradient_per_s.iter().map(|g| g.norm_sq()).sum::<f64>(),
                        2e-14);
                }
            }
        }
    }
    near(integrated_flow, flow, 3e-6);
    let work = prepared.dissipation_per_axial_length(pressure, flow).unwrap();
    relative(work.viscous_w_m, heat[0], 2e-5);
    relative(work.thermal_w_m, heat[1], 2e-5);
    let wave = prepared.wave().wave;
    relative(work.viscous_w_m, 0.5 * (wave.wavenumber * wave.characteristic_impedance).im * flow.norm_sq(), 3e-13);
    relative(work.thermal_w_m, 0.5 * (wave.wavenumber / wave.characteristic_impedance).im * pressure.norm_sq(), 3e-13);
    assert_eq!(prepared.dissipation_per_axial_length(C64::ZERO, flow).unwrap().thermal_w_m, 0.0);
    assert_eq!(prepared.dissipation_per_axial_length(pressure, C64::ZERO).unwrap().viscous_w_m, 0.0);
}

// Independent TWO-dimensional separation: no analytic cosh elimination or
// mean/gradient-integral helper from the production field owner is used.
fn double_sine(a: f64, b: f64, y: f64, z: f64, omega: f64, diffusivity: f64) -> C64 {
    let mut sum = C64::ZERO;
    for i in 0..384 {
        let m = (2 * i + 1) as f64;
        for j in 0..384 {
            let n = (2 * j + 1) as f64;
            let rate = diffusivity * PI * PI * (m * m / (a * a) + n * n / (b * b));
            let weight = 16.0 * (m * PI * y / a).sin() * (n * PI * z / b).sin() / (PI * PI * m * n);
            sum = sum + (C64::new(0.0, -omega) / C64::new(rate, -omega)).scale(weight);
        }
    }
    sum
}

#[test]
fn transverse_fields_converge_to_an_independent_two_dimensional_diffusion_solution() {
    let material = fluid();
    let (a, b, omega) = (0.001, 0.0017, 2000.0);
    let points = [(0.17, 0.23), (0.37, 0.41), (0.68, 0.76)];
    let references = points.map(|(y, z)| double_sine(a, b, y * a, z * b, omega,
        material.transport().thermal_diffusivity_m2_s()));
    let mut errors = Vec::new();
    for terms in [8, 32, 256] {
        let section = RectangularSection::new(a, b, terms).unwrap();
        let mut controls = admission();
        controls.maximum_relative_tail = 0.2;
        let prepared = section.prepare_fields(&material, omega, controls).unwrap();
        let error = points.into_iter().zip(references).map(|((y, z), reference)| {
            let field = prepared.field_at(y * a, z * b, C64::ONE, C64::ZERO).unwrap();
            (field.thermal_pressure_response_pa.unwrap() - reference).norm_sq()
        }).sum::<f64>().sqrt();
        errors.push(error);
    }
    assert!(errors[1] < errors[0] * 0.1 && errors[2] < errors[1] * 0.1, "{errors:?}");
    assert!(errors[2] < 2e-7, "{errors:?}");
}

fn travel(k: C64, x: f64) -> C64 {
    C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp())
}

#[test]
fn coherent_counterwaves_lose_the_actual_local_irreversible_work() {
    let material = fluid();
    let section = RectangularSection::new(0.001, 0.002, 128).unwrap();
    let prepared = section.prepare_fields(&material, 4000.0, admission()).unwrap();
    let wave = prepared.wave().wave;
    let fields = |x: f64| {
        let forward = C64::new(20.0, 3.0) * travel(wave.wavenumber, x);
        let backward = C64::new(-6.0, 7.0) * travel(wave.wavenumber, 0.2 - x);
        (forward + backward, (forward - backward) / wave.characteristic_impedance)
    };
    let power = |x: f64| {
        let (p, q) = fields(x);
        0.5 * (p * q.conj()).re
    };
    for x in [0.03, 0.07, 0.13, 0.18] {
        let (p, q) = fields(x);
        let work = prepared.dissipation_per_axial_length(p, q).unwrap();
        let h = 1e-6;
        relative(work.total_w_m, (power(x - h) - power(x + h)) / (2.0 * h), 3e-8);
        assert!(work.viscous_w_m > 0.0 && work.thermal_w_m > 0.0);
    }
}

#[test]
fn cp_equal_cv_retains_only_original_liquid_shear_and_side_interchange() {
    let material = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1500.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap(),
        4200.0,
    ).unwrap();
    let section = RectangularSection::new(0.001, 0.002, 128).unwrap();
    let swapped = RectangularSection::new(0.002, 0.001, 128).unwrap();
    let prepared = section.prepare_fields(&material, 1000.0, admission()).unwrap();
    let reversed = swapped.prepare_fields(&material, 1000.0, admission()).unwrap();
    let (p, q) = (C64::new(100.0, 20.0), C64::new(2e-9, -7e-9));
    let field = prepared.field_at(0.0003, 0.0007, p, q).unwrap();
    assert_eq!(field, reversed.field_at(0.0003, 0.0007, p, q).unwrap());
    assert_eq!(field.thermal_pressure_response_pa, None);
    assert_eq!(field.thermal_pressure_gradient_pa_m, None);
    assert_eq!(field.thermal_dissipation_w_m3, 0.0);
    let work = prepared.dissipation_per_axial_length(p, q).unwrap();
    assert_eq!(work.thermal_w_m, 0.0);
    assert!(work.viscous_w_m > 0.0);
}

struct OneSnapshot {
    state: ThermoviscousFluidState,
    queries: Cell<usize>,
}
impl FluidAcoustics for OneSnapshot {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        panic!("the original owner must only be queried for its complete state")
    }
}
impl FluidTransport for OneSnapshot {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("the original owner must only be queried for its complete state")
    }
}
impl ThermoviscousFluid for OneSnapshot {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.queries.set(self.queries.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn one_source_query_bounded_fields_and_atomic_refusals() {
    let source = OneSnapshot { state: fluid(), queries: Cell::new(0) };
    let section = RectangularSection::new(0.001, 0.002, 128).unwrap();
    let prepared = section.prepare_fields(&source, 1000.0, admission()).unwrap();
    assert!(std::ptr::eq(prepared.source(), &source));
    assert_eq!(*prepared.state(), source.state);
    let (p, q) = (C64::new(10.0, 2.0), C64::new(2e-9, -7e-9));
    for _ in 0..5 {
        prepared.field_at(0.0003, 0.0007, p, q).unwrap();
        prepared.dissipation_per_axial_length(p, q).unwrap();
    }
    assert_eq!(source.queries.get(), 1);
    for (y, z) in [(-1e-6, 0.001), (0.002, 0.001), (0.0005, -1e-6), (0.0005, 0.003),
        (f64::NAN, 0.001), (0.0005, f64::INFINITY)]
    {
        assert!(prepared.field_at(y, z, p, q).is_err());
    }
    assert!(prepared.field_at(0.0003, 0.0007, C64::new(f64::NAN, 0.0), q).is_err());
    assert!(prepared.dissipation_per_axial_length(p, C64::from_re(f64::INFINITY)).is_err());
    assert!(matches!(prepared.field_at_controlled(0.0003, 0.0007, p, q, || true),
        Err(DuctError::Cancelled { .. })));
    assert!(matches!(section.prepare_fields_controlled(&source, 1000.0, admission(), || true),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))));
    assert_eq!(source.queries.get(), 1);
    let mut polls = 0;
    assert!(matches!(section.prepare_fields_controlled(&source, 1000.0, admission(), || {
        polls += 1;
        polls == 12
    }), Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))));
    let unresolved = RectangularSection::new(0.001, 0.002, 1).unwrap();
    assert!(unresolved.prepare_fields(&source, 1000.0, admission()).is_err());
    assert!(section.prepare_fields(&source, 1e8, admission()).is_err());
}
