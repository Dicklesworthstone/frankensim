//! Local heat diffusion, internal-reservoir temperatures and positive loss in
//! the original rectangular/slit molecular wave.

use fs_duct::DuctError;
use fs_duct::fluid::FluidDuctError;
use fs_duct::fluid::rectangular::{RectangularAdmission, RectangularSection};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;
use std::f64::consts::PI;

// Manufactured reservoirs test the equations, not measured air relaxation.
// The production pressure-dependent material law evaluates every source process.
fn material(count: usize, time_scale: f64) -> RelaxationState {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let species = SpeciesId::new("synthetic").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: (0..count)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * i as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (0.002 + 0.01 * i as f64) * time_scale,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap();
    law.evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
        .unwrap()
}

fn admission() -> RectangularAdmission {
    RectangularAdmission {
        maximum_relative_tail: 1e-5,
        maximum_transverse_phase: 0.8,
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:?} != {expected:?}"
    );
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn absent_reservoirs_preserve_the_original_wave_and_resolved_equilibrium_fields() {
    let state = material(0, 1.0);
    let section = RectangularSection::new(0.0004, 0.0015, 128).unwrap();
    for omega in [20.0, 1000.0, 5000.0] {
        let original = section
            .prepare_fields(state.gas(), omega, admission())
            .unwrap();
        let molecular = section
            .prepare_molecular_fields(&state, omega, admission())
            .unwrap();
        assert!(std::ptr::eq(molecular.source(), &state));
        assert!(std::ptr::eq(molecular.section(), &section));
        assert_eq!(
            molecular.wave().wave.wavenumber,
            original.wave().wave.wavenumber
        );
        assert_eq!(
            molecular.wave().wave.characteristic_impedance,
            original.wave().wave.characteristic_impedance
        );
        let (p, q) = (C64::new(3.0, -2.0), C64::new(2e-8, 1e-8));
        let e = original.field_at(0.00013, 0.00062, p, q).unwrap();
        let m = molecular.field_at(0.00013, 0.00062, p, q).unwrap();
        assert_eq!(m.axial_velocity_m_s, e.axial_velocity_m_s);
        assert_eq!(
            m.transverse_velocity_gradient_per_s,
            e.transverse_velocity_gradient_per_s
        );
        near(
            m.gas_temperature_k,
            original.gas_temperature_at(0.00013, 0.00062, p).unwrap(),
            3e-15,
        );
        relative(
            m.conduction_dissipation_w_m3,
            e.thermal_dissipation_w_m3,
            3e-14,
        );
        assert_eq!(m.relaxation_dissipation_w_m3, 0.0);
        let work = molecular.dissipation_per_axial_length(p, q).unwrap();
        let old_work = original.dissipation_per_axial_length(p, q).unwrap();
        relative(work.total_w_m, old_work.total_w_m, 3e-14);
        assert_eq!(work.relaxation_w_m, 0.0);
        assert!(molecular.reservoir_field_at(0, 0.0001, 0.0005, p).is_err());
    }
}

// Independent TWO-dimensional eigenmode energy equation in the native sign.
// Production eliminates the long direction analytically and sums only one axis.
fn double_sine_temperature(
    state: &RelaxationState,
    omega: f64,
    a: f64,
    b: f64,
    y: f64,
    z: f64,
    p: C64,
) -> C64 {
    let gas = state.gas();
    let s = C64::new(0.0, -omega);
    let cp = state.processes().iter().fold(
        C64::from_re(state.frozen_cv_j_kg_k() + gas.specific_gas_constant),
        |sum, r| {
            sum + C64::from_re(r.heat_capacity_j_kg_k) / C64::new(1.0, -omega * r.relaxation_time_s)
        },
    );
    let mut temperature = C64::ZERO;
    for i in 0..384 {
        let m = (2 * i + 1) as f64;
        let sy = (m * PI * y / a).sin();
        for j in 0..384 {
            let n = (2 * j + 1) as f64;
            let rate = gas.thermal_conductivity / gas.density
                * PI
                * PI
                * (m * m / (a * a) + n * n / (b * b));
            let weight = 16.0 * sy * (n * PI * z / b).sin() / (PI * PI * m * n);
            temperature = temperature + (s / (s * cp + C64::from_re(rate))).scale(weight);
        }
    }
    temperature * p / C64::from_re(gas.density)
}

#[test]
fn molecular_wall_temperature_solves_the_independent_two_dimensional_energy_equation() {
    let state = material(2, 1.0);
    let (a, b) = (0.0007, 0.0015);
    let p = C64::new(2.0, -0.4);
    for omega in [50.0, 1200.0, 5000.0] {
        let y = 0.31 * a;
        let z = 0.43 * b;
        let expected = double_sine_temperature(&state, omega, a, b, y, z, p);
        let mut previous = f64::INFINITY;
        for terms in [16, 64, 256] {
            let section = RectangularSection::new(a, b, terms).unwrap();
            let fields = section
                .prepare_molecular_fields(
                    &state,
                    omega,
                    RectangularAdmission {
                        maximum_relative_tail: 0.1,
                        maximum_transverse_phase: 0.8,
                    },
                )
                .unwrap();
            let temperature = fields
                .field_at(y, z, p, C64::ZERO)
                .unwrap()
                .gas_temperature_k;
            let error = (temperature - expected).abs();
            assert!(
                error < previous,
                "refinement did not improve the coupled field"
            );
            previous = error;
            if terms == 256 {
                near(temperature, expected, 4e-7);
            }
        }
    }
}

#[test]
fn original_reservoir_identity_thermal_lag_and_all_four_isothermal_walls_are_retained() {
    let state = material(2, 1.0);
    let section = RectangularSection::new(0.0006, 0.0018, 128).unwrap();
    let omega = 900.0;
    let fields = section
        .prepare_molecular_fields(&state, omega, admission())
        .unwrap();
    let p = C64::new(3.0, 1.0);
    let local = fields
        .field_at(0.0002, 0.0011, p, C64::new(2e-8, 1e-8))
        .unwrap();
    let mut work = 0.0;
    for (index, process) in state.processes().iter().enumerate() {
        let r = fields.reservoir_field_at(index, 0.0002, 0.0011, p).unwrap();
        assert!(std::ptr::eq(r.process, process));
        near(r.gas_temperature_k, local.gas_temperature_k, 2e-15);
        near(
            r.reservoir_temperature_k,
            local.gas_temperature_k / C64::new(1.0, -omega * process.relaxation_time_s),
            2e-15,
        );
        near(
            C64::new(0.0, -omega) * r.reservoir_temperature_k,
            r.temperature_difference_k / C64::from_re(process.relaxation_time_s),
            3e-15,
        );
        relative(
            r.relaxation_dissipation_w_m3,
            state.gas().density * process.heat_capacity_j_kg_k
                / (2.0 * state.gas().temperature * process.relaxation_time_s)
                * r.temperature_difference_k.norm_sq(),
            3e-15,
        );
        work += r.relaxation_dissipation_w_m3;
    }
    relative(work, local.relaxation_dissipation_w_m3, 3e-15);
    for (y, z) in [
        (0.0, 0.0009),
        (0.0006, 0.0009),
        (0.0003, 0.0),
        (0.0003, 0.0018),
    ] {
        let wall = fields.field_at(y, z, p, C64::from_re(1e-8)).unwrap();
        assert_eq!(wall.gas_temperature_k, C64::ZERO);
        assert_eq!(wall.axial_velocity_m_s, C64::ZERO);
        assert_eq!(wall.relaxation_dissipation_w_m3, 0.0);
        assert_eq!(
            fields
                .reservoir_field_at(0, y, z, p)
                .unwrap()
                .reservoir_temperature_k,
            C64::ZERO
        );
    }
}

#[test]
fn positive_gradient_and_reservoir_work_close_the_original_native_wave_power_balance() {
    let state = material(2, 1.0);
    let section = RectangularSection::new(0.0005, 0.0015, 128).unwrap();
    let (p, q) = (C64::new(3.0, -2.0), C64::new(1e-8, 3e-8));
    for omega in [1.0, 50.0, 500.0, 3000.0, 10000.0] {
        let fields = section
            .prepare_molecular_fields(&state, omega, admission())
            .unwrap();
        let wave = fields.wave().wave;
        let work = fields.dissipation_per_axial_length(p, q).unwrap();
        let shear = 0.5 * (wave.wavenumber * wave.characteristic_impedance).im * q.norm_sq();
        let heat = 0.5 * (wave.wavenumber / wave.characteristic_impedance).im * p.norm_sq();
        relative(work.viscous_w_m, shear, 3e-12);
        relative(work.conduction_w_m + work.relaxation_w_m, heat, 3e-12);
        relative(work.total_w_m, shear + heat, 3e-12);
        assert!(work.viscous_w_m > 0.0 && work.conduction_w_m > 0.0 && work.relaxation_w_m > 0.0);
        let individual = (0..state.processes().len())
            .map(|i| fields.reservoir_dissipation_per_axial_length(i, p).unwrap())
            .sum::<f64>();
        relative(individual, work.relaxation_w_m, 3e-15);
    }
}

#[test]
fn local_work_quadrature_recovers_each_analytic_positive_channel() {
    let state = material(2, 1.0);
    let (a, b) = (0.0005, 0.0015);
    let section = RectangularSection::new(a, b, 128).unwrap();
    let fields = section
        .prepare_molecular_fields(&state, 1500.0, admission())
        .unwrap();
    let (p, q) = (C64::new(2.0, 0.7), C64::new(1.5e-8, -1e-8));
    let (ny, nz) = (96, 144);
    let weight = |i: usize, n: usize| {
        if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    let mut integral = [0.0; 3];
    for i in 0..=ny {
        for j in 0..=nz {
            let field = fields
                .field_at(a * i as f64 / ny as f64, b * j as f64 / nz as f64, p, q)
                .unwrap();
            let w = weight(i, ny) * weight(j, nz) * a * b / (9.0 * ny as f64 * nz as f64);
            for (sum, value) in integral.iter_mut().zip([
                field.viscous_dissipation_w_m3,
                field.conduction_dissipation_w_m3,
                field.relaxation_dissipation_w_m3,
            ]) {
                *sum += w * value;
            }
        }
    }
    let work = fields.dissipation_per_axial_length(p, q).unwrap();
    for (measured, exact) in
        integral
            .into_iter()
            .zip([work.viscous_w_m, work.conduction_w_m, work.relaxation_w_m])
    {
        relative(measured, exact, 2e-4);
    }
}

#[test]
fn equilibrated_and_frozen_reservoir_lags_remain_representable_without_temperature_subtraction() {
    let section = RectangularSection::new(0.0006, 0.0015, 128).unwrap();
    for time_scale in [1e-14, 1e12] {
        let state = material(1, time_scale);
        let fields = section
            .prepare_molecular_fields(&state, 1000.0, admission())
            .unwrap();
        let r = fields
            .reservoir_field_at(0, 0.0003, 0.00075, C64::ONE)
            .unwrap();
        near(
            C64::new(0.0, -1000.0) * r.reservoir_temperature_k,
            r.temperature_difference_k / C64::from_re(r.process.relaxation_time_s),
            4e-15,
        );
        assert!(r.relaxation_dissipation_w_m3.is_finite() && r.relaxation_dissipation_w_m3 > 0.0);
        let ratio = (r.reservoir_temperature_k / r.gas_temperature_k).abs();
        if time_scale < 1.0 {
            assert!((ratio - 1.0).abs() < 1e-14);
        } else {
            assert!(ratio < 1e-10);
        }
    }
}

#[test]
fn molecular_field_domain_source_index_and_bounded_cancellation_are_explicit() {
    let state = material(2, 1.0);
    let section = RectangularSection::new(0.0005, 0.0015, 128).unwrap();
    assert!(matches!(
        section.prepare_molecular_fields_controlled(&state, 1000.0, admission(), || true),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let fields = section
        .prepare_molecular_fields(&state, 1000.0, admission())
        .unwrap();
    for (y, z) in [
        (-1e-8, 0.0),
        (0.00050001, 0.0),
        (0.0, f64::NAN),
        (0.0, 0.00150001),
    ] {
        assert!(fields.field_at(y, z, C64::ONE, C64::ZERO).is_err());
    }
    assert!(
        fields
            .reservoir_field_at(2, 0.0002, 0.0008, C64::ONE)
            .is_err()
    );
    assert!(
        fields
            .reservoir_dissipation_per_axial_length(2, C64::ONE)
            .is_err()
    );
    assert!(
        fields
            .dissipation_per_axial_length(C64::new(f64::NAN, 0.0), C64::ZERO)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        fields.field_at_controlled(0.0002, 0.0008, C64::ONE, C64::ZERO, || {
            polls += 1;
            polls == 2
        }),
        Err(DuctError::Cancelled {
            completed_steps: 32
        })
    ));
    assert!(matches!(
        fields.reservoir_field_at_controlled(0, 0.0002, 0.0008, C64::ONE, || true),
        Err(DuctError::Cancelled { completed_steps: 0 })
    ));
}
