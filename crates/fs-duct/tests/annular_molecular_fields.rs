//! Original molecular temperature fields and independent positive loss forms.
//! Authored reservoirs exercise the equations, not experimental calibration.

use fs_duct::DuctError;
use fs_duct::fluid::FluidDuctError;
use fs_duct::fluid::annular::{AnnularAdmission, AnnularSection};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;

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

fn admission() -> AnnularAdmission {
    AnnularAdmission {
        maximum_transverse_phase: 0.8,
        maximum_cell_diffusion_phase: 0.75,
        maximum_relative_residual: 1e-8,
        maximum_payload_bytes: 2 << 20,
    }
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e} at {tolerance}"
    );
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?} at {tolerance}"
    );
}
fn cp(state: &RelaxationState, omega: f64) -> C64 {
    state.processes().iter().fold(
        C64::from_re(state.frozen_cv_j_kg_k() + state.gas().specific_gas_constant),
        |value, process| {
            value
                + C64::from_re(process.heat_capacity_j_kg_k)
                    / C64::new(1.0, -omega * process.relaxation_time_s)
        },
    )
}

// Independent original reservoir temperatures obey (1-i*omega*tau)*Tj=T.
// The positive transfer form uses the actual temperature difference; neither
// Im(Cp), the wave attenuation nor a retained mass/loss matrix is the oracle.
fn reservoir_loss(state: &RelaxationState, omega: f64, temperature: C64) -> f64 {
    let work = state
        .processes()
        .iter()
        .map(|process| {
            let internal = temperature / C64::new(1.0, -omega * process.relaxation_time_s);
            process.heat_capacity_j_kg_k / process.relaxation_time_s
                * (temperature - internal).norm_sq()
        })
        .sum::<f64>();
    0.5 * state.gas().density / state.gas().temperature * work
}

#[test]
fn retained_temperature_and_native_mass_integrals_recover_each_original_loss() {
    let state = material(2, 1.0);
    let gas = state.gas();
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    let pressure = C64::new(30.0, 7.0);
    let flow = C64::new(4e-8, -3e-8);
    for omega in [50.0, 500.0, 3000.0] {
        let fields = section
            .prepare_molecular_fields(&state, omega, admission())
            .unwrap();
        let native = section
            .evaluate_molecular(&state, omega, admission())
            .unwrap();
        assert!(core::ptr::eq(fields.source(), &state));
        assert!(core::ptr::eq(fields.section(), &section));
        assert_eq!(fields.wave().wave.wavenumber, native.wave.wavenumber);
        assert_eq!(
            fields.wave().wave.characteristic_impedance,
            native.wave.characteristic_impedance
        );
        near(fields.complex_cp_j_kg_k(), cp(&state, omega), 2e-15);
        let mut total = [0.0; 3];
        for cell in 0..section.elements() {
            let left = section.radius_m(cell).unwrap();
            let right = section.radius_m(cell + 1).unwrap();
            let temperatures = [0.0, 1.0].map(|fraction| {
                fields
                    .gas_temperature_in_cell(cell, fraction, pressure)
                    .unwrap()
            });
            let gradient = (temperatures[1] - temperatures[0]).scale(1.0 / (right - left));
            let mut oracle = [0.0; 3];
            // Two physical Gauss points exactly integrate r times the square
            // of a P1 temperature; no production mass weights are consulted.
            for fraction in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
                let field = fields
                    .field_in_cell(cell, fraction, pressure, flow)
                    .unwrap();
                let temperature = fields
                    .gas_temperature_in_cell(cell, fraction, pressure)
                    .unwrap();
                near(
                    temperature,
                    field.thermal_pressure_response_pa.unwrap()
                        / C64::from_re(gas.density)
                        / cp(&state, omega),
                    3e-15,
                );
                let independent = [
                    0.5 * gas.dynamic_viscosity * field.radial_velocity_gradient_per_s.norm_sq(),
                    0.5 * gas.thermal_conductivity / gas.temperature * gradient.norm_sq(),
                    reservoir_loss(&state, omega, temperature),
                ];
                for (actual, expected) in [
                    field.viscous_dissipation_w_m3,
                    field.thermal_dissipation_w_m3,
                    field.molecular_dissipation_w_m3,
                ]
                .into_iter()
                .zip(independent)
                {
                    relative(actual, expected, 3e-11);
                }
                for (sum, value) in oracle.iter_mut().zip(independent) {
                    *sum += core::f64::consts::PI * field.radius_m * (right - left) * value;
                }
            }
            let ring = fields
                .cell_dissipation_per_axial_length(cell, pressure, flow)
                .unwrap();
            for ((sum, actual), expected) in total
                .iter_mut()
                .zip([ring.viscous_w_m, ring.thermal_w_m, ring.molecular_w_m])
                .zip(oracle)
            {
                relative(actual, expected, 3e-11);
                *sum += actual;
            }
        }
        let whole = fields.dissipation_per_axial_length(pressure, flow).unwrap();
        for (actual, expected) in [whole.viscous_w_m, whole.thermal_w_m, whole.molecular_w_m]
            .into_iter()
            .zip(total)
        {
            assert!(actual > 0.0);
            relative(actual, expected, 3e-14);
        }
        let wave = native.wave;
        relative(
            whole.thermal_w_m + whole.molecular_w_m,
            0.5 * (wave.wavenumber / wave.characteristic_impedance).im * pressure.norm_sq(),
            2e-8,
        );
        for (cell, fraction) in [(0, 0.0), (section.elements() - 1, 1.0)] {
            let wall = fields
                .field_in_cell(cell, fraction, pressure, flow)
                .unwrap();
            assert_eq!(wall.axial_velocity_m_s, C64::ZERO);
            assert_eq!(wall.molecular_dissipation_w_m3, 0.0);
            assert_eq!(
                fields
                    .gas_temperature_in_cell(cell, fraction, pressure)
                    .unwrap(),
                C64::ZERO
            );
        }
    }
}

#[test]
fn no_reservoirs_preserve_original_fields_temperature_and_finite_integrals_exactly() {
    let state = material(0, 1.0);
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    let pressure = C64::new(2.0, -3.0);
    let flow = C64::new(3e-8, 7e-9);
    for omega in [0.001, 100.0, 3000.0] {
        let original = section
            .prepare_fields(state.gas(), omega, admission())
            .unwrap();
        let molecular = section
            .prepare_molecular_fields(&state, omega, admission())
            .unwrap();
        assert!(core::ptr::eq(molecular.source(), &state));
        for (cell, fraction) in [(0, 0.0), (30, 0.37), (127, 1.0)] {
            assert_eq!(
                molecular
                    .field_in_cell(cell, fraction, pressure, flow)
                    .unwrap(),
                original
                    .field_in_cell(cell, fraction, pressure, flow)
                    .unwrap()
            );
            assert_eq!(
                molecular
                    .gas_temperature_in_cell(cell, fraction, pressure)
                    .unwrap(),
                original
                    .gas_temperature_in_cell(cell, fraction, pressure)
                    .unwrap()
            );
        }
        let reference = state.gas().characteristic_impedance / section.area_m2();
        let a = original
            .prepare_excitation(0.2, reference, [pressure, -pressure])
            .unwrap();
        let b = molecular
            .prepare_excitation(0.2, reference, [pressure, -pressure])
            .unwrap();
        assert_eq!(a.two_port().scattering, b.two_port().scattering);
        assert_eq!(
            a.dissipation_between(0.02, 0.17).unwrap(),
            b.dissipation_between(0.02, 0.17).unwrap()
        );
        assert_eq!(
            b.dissipation_between(0.0, 0.2).unwrap().molecular_power_w,
            0.0
        );
    }
}

#[test]
fn coherent_molecular_segments_close_all_channels_partition_and_retain_node_work() {
    let state = material(2, 1.0);
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    let fields = section
        .prepare_molecular_fields(&state, 1000.0, admission())
        .unwrap();
    let reference = state.gas().characteristic_impedance / section.area_m2();
    let length = 0.23;
    let drive = [C64::new(20.0, 3.0), C64::new(-7.0, 9.0)];
    let excitation = fields.prepare_excitation(length, reference, drive).unwrap();
    assert!(core::ptr::eq(excitation.fields(), &fields));
    let whole = excitation.dissipation_between(0.0, length).unwrap();
    relative(
        whole.absorbed_power_w,
        excitation.boundary_power().absorbed_power_w,
        4e-8,
    );
    let mut partition = [0.0; 3];
    for interval in 0..7 {
        let start = length * interval as f64 / 7.0;
        let end = length * (interval + 1) as f64 / 7.0;
        for cell in 0..section.elements() {
            let power = excitation
                .cell_dissipation_between(cell, start, end)
                .unwrap();
            for (sum, value) in partition.iter_mut().zip([
                power.viscous_power_w,
                power.thermal_power_w,
                power.molecular_power_w,
            ]) {
                *sum += value;
            }
        }
    }
    for (actual, expected) in partition.into_iter().zip([
        whole.viscous_power_w,
        whole.thermal_power_w,
        whole.molecular_power_w,
    ]) {
        relative(actual, expected, 4e-13);
    }
    // Independent axial Simpson integration of the observed physical molecular
    // temperature. Its radial P1 square is integrated at physical Gauss points.
    let mut radial = 0.0;
    for cell in 0..section.elements() {
        let left = section.radius_m(cell).unwrap();
        let right = section.radius_m(cell + 1).unwrap();
        for fraction in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
            let radius = left + fraction * (right - left);
            let temperature = fields
                .gas_temperature_in_cell(cell, fraction, C64::ONE)
                .unwrap();
            radial += core::f64::consts::PI
                * radius
                * (right - left)
                * reservoir_loss(&state, 1000.0, temperature);
        }
    }
    let mut integral = 0.0;
    let panels = 1024;
    for i in 0..=panels {
        let field = excitation
            .at_depth(length * i as f64 / panels as f64)
            .unwrap();
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * radial * field.pressure_pa.norm_sq();
    }
    relative(
        whole.molecular_power_w,
        integral * length / (3.0 * panels as f64),
        2e-10,
    );
    let axial = excitation.at_depth(0.07).unwrap();
    assert_eq!(
        excitation.gas_temperature_in_cell(30, 0.4, 0.07).unwrap(),
        fields
            .gas_temperature_in_cell(30, 0.4, axial.pressure_pa)
            .unwrap()
    );

    let standing = fields
        .prepare_excitation(length, reference, [C64::ONE, -C64::ONE])
        .unwrap();
    let center = 0.5 * length;
    let start = center - 0.5e-9;
    let end = center + 0.5e-9;
    let middle = standing.at_depth(center).unwrap();
    assert_eq!(middle.pressure_pa, C64::ZERO);
    let wave = fields.wave().wave;
    let gradient = C64::new(0.0, 1.0)
        * wave.wavenumber
        * wave.characteristic_impedance
        * middle.volume_flow_m3_s;
    let node = standing.dissipation_between(start, end).unwrap();
    assert!(node.molecular_power_w > 0.0);
    relative(
        node.molecular_power_w,
        radial * gradient.norm_sq() * (end - start).powi(3) / 12.0,
        5e-8,
    );
    let opaque_length = 1000.0 / wave.wavenumber.im;
    let opaque = fields
        .prepare_excitation(opaque_length, reference, drive)
        .unwrap();
    let power = opaque.dissipation_between(0.0, opaque_length).unwrap();
    assert!(power.molecular_power_w.is_finite() && power.molecular_power_w > 0.0);
    relative(
        power.absorbed_power_w,
        opaque.boundary_power().absorbed_power_w,
        4e-8,
    );
}

#[test]
fn isothermal_wall_suppression_extreme_rates_and_original_refusals_are_retained() {
    let state = material(2, 1.0);
    let narrow = AnnularSection::new(0.0004, 0.00042, 128).unwrap();
    let wide = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    let molecular_per_volume = |section: &AnnularSection| {
        section
            .prepare_molecular_fields(&state, 500.0, admission())
            .unwrap()
            .dissipation_per_axial_length(C64::ONE, C64::ZERO)
            .unwrap()
            .molecular_w_m
            / section.area_m2()
    };
    assert!(
        molecular_per_volume(&narrow) < 1e-3 * molecular_per_volume(&wide),
        "the isothermal walls must suppress the original reservoir temperature, not add a bulk attenuation"
    );
    for scale in [1e-9, 1e9] {
        let limiting = material(2, scale);
        let fields = wide
            .prepare_molecular_fields(&limiting, 500.0, admission())
            .unwrap();
        let native = wide
            .evaluate_molecular(&limiting, 500.0, admission())
            .unwrap();
        assert_eq!(fields.wave().wave.wavenumber, native.wave.wavenumber);
        let work = fields
            .dissipation_per_axial_length(C64::ONE, C64::ZERO)
            .unwrap();
        assert!(work.molecular_w_m.is_finite() && work.molecular_w_m > 0.0);
        relative(
            work.thermal_w_m + work.molecular_w_m,
            0.5 * (native.wave.wavenumber / native.wave.characteristic_impedance).im,
            2e-8,
        );
    }
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            wide.prepare_molecular_fields(&state, omega, admission())
                .is_err()
        );
    }
    let mut polls = 0;
    wide.prepare_molecular_fields_controlled(&state, 500.0, admission(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut calls = 0;
    assert!(matches!(
        wide.prepare_molecular_fields_controlled(&state, 500.0, admission(), || {
            calls += 1;
            calls >= polls - 1
        }),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let mut insufficient = admission();
    insufficient.maximum_payload_bytes = 1;
    assert!(
        wide.prepare_molecular_fields(&state, 500.0, insufficient)
            .is_err()
    );
    let fields = wide
        .prepare_molecular_fields(&state, 500.0, admission())
        .unwrap();
    assert!(
        fields
            .gas_temperature_in_cell(wide.elements(), 0.5, C64::ONE)
            .is_err()
    );
    assert!(fields.field_in_cell(0, -0.1, C64::ONE, C64::ZERO).is_err());
    assert!(
        fields
            .field_in_cell(0, 0.5, C64::new(f64::NAN, 0.0), C64::ZERO)
            .is_err()
    );
}
