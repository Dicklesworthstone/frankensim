//! G1/G3: the original molecular energy equation on both actual annular walls.
//! Manufactured reservoir data test the equations, not empirical rate validity.

use fs_duct::DuctError;
use fs_duct::fluid::FluidDuctError;
use fs_duct::fluid::annular::{AnnularAdmission, AnnularSection, AnnularWave};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;

fn material(count: usize, time_scale: f64, pressure: f64) -> RelaxationState {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
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
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:?} != {expected:?} at tolerance {tolerance}"
    );
}
fn compliance(response: AnnularWave, omega: f64) -> C64 {
    response.wave.wavenumber
        / response.wave.characteristic_impedance
        / C64::from_re(omega * response.area_m2())
}

// Direct +i omega reservoir sum; independent of caloric_response's Debye arithmetic.
fn heat_capacity(state: &RelaxationState, omega: f64) -> C64 {
    state.processes().iter().fold(
        C64::from_re(state.frozen_cv_j_kg_k() + state.gas().specific_gas_constant),
        |cp, p| {
            cp + C64::from_re(p.heat_capacity_j_kg_k) / C64::new(1.0, omega * p.relaxation_time_s)
        },
    )
}

#[test]
fn absent_reservoirs_keep_the_original_annular_fields_and_scattering_exactly() {
    let state = material(0, 1.0, 101325.0);
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    for omega in [0.001, 100.0, 3000.0] {
        let original = section.evaluate(state.gas(), omega, admission()).unwrap();
        let molecular = section
            .evaluate_molecular(&state, omega, admission())
            .unwrap();
        assert_eq!(molecular.wave.wavenumber, original.wave.wavenumber);
        assert_eq!(
            molecular.wave.characteristic_impedance,
            original.wave.characteristic_impedance
        );
        for (actual, expected) in [molecular.viscous, molecular.thermal.unwrap()]
            .into_iter()
            .zip([original.viscous, original.thermal.unwrap()])
        {
            assert_eq!(actual.mean, expected.mean);
            assert_eq!(actual.wall, expected.wall);
            assert_eq!(actual.gradient_work, expected.gradient_work);
            assert_eq!(actual.relative_work_defect, expected.relative_work_defect);
        }
        let reference = state.gas().characteristic_impedance / section.area_m2();
        assert_eq!(
            molecular.two_port(0.1, reference).unwrap().scattering,
            original.two_port(0.1, reference).unwrap().scattering
        );
    }
}

// Independent primitive-temperature shooting, theta = rho*T'/p':
// theta_xx + gap/r theta_x = gap²*rho/κ * (s*cp(s)*theta - s).
// Uniform-radius RK4 integrates a forced and a homogeneous solution plus
// their area means; the second wall fixes the unknown initial derivative.
// No P1 matrices, logarithmic mesh, tridiagonal solve or nodal quadrature.
fn temperature_reference(
    a: f64,
    b: f64,
    omega: f64,
    state: &RelaxationState,
    panels: usize,
) -> C64 {
    let gap = b - a;
    let gas = state.gas();
    let s = C64::new(0.0, omega);
    let cp = heat_capacity(state, omega);
    let scale = gap * gap * gas.density / gas.thermal_conductivity;
    let derivative = |x: f64, y: [C64; 6]| {
        let r = a + gap * x;
        let weight = 2.0 * r / (a + b);
        [
            y[1],
            (s * cp * y[0] - s).scale(scale) - y[1].scale(gap / r),
            y[0].scale(weight),
            y[4],
            (s * cp * y[3]).scale(scale) - y[4].scale(gap / r),
            y[3].scale(weight),
        ]
    };
    let mut y = [C64::ZERO; 6];
    y[4] = C64::ONE;
    let step = 1.0 / panels as f64;
    for i in 0..panels {
        let x = i as f64 * step;
        let k1 = derivative(x, y);
        let k2 = derivative(
            x + 0.5 * step,
            core::array::from_fn(|j| y[j] + k1[j].scale(0.5 * step)),
        );
        let k3 = derivative(
            x + 0.5 * step,
            core::array::from_fn(|j| y[j] + k2[j].scale(0.5 * step)),
        );
        let k4 = derivative(x + step, core::array::from_fn(|j| y[j] + k3[j].scale(step)));
        for j in 0..6 {
            y[j] = y[j] + (k1[j] + k2[j].scale(2.0) + k3[j].scale(2.0) + k4[j]).scale(step / 6.0);
        }
    }
    y[2] - (y[0] / y[3]) * y[5]
}

#[test]
fn radial_refinement_converges_to_the_independent_molecular_energy_equation() {
    let state = material(2, 1.0, 101325.0);
    let gas = state.gas();
    let (a, b) = (0.0004, 0.0012);
    for omega in [50.0, 500.0, 3000.0] {
        let temperature = temperature_reference(a, b, omega, &state, 4096);
        near(
            temperature,
            temperature_reference(a, b, omega, &state, 2048),
            2e-9,
        );
        let expected_h = temperature * heat_capacity(&state, omega);
        let expected_compliance = (C64::ONE - temperature.scale(gas.specific_gas_constant))
            .scale(1.0 / gas.pressure)
            .conj();
        let mut previous = f64::INFINITY;
        for elements in [64, 128, 256] {
            let section = AnnularSection::new(a, b, elements).unwrap();
            let wave = section
                .evaluate_molecular(&state, omega, admission())
                .unwrap();
            let thermal = wave.thermal.unwrap();
            let error = (thermal.mean - expected_h).abs();
            if previous.is_finite() {
                assert!(
                    (3.8..4.2).contains(&(previous / error)),
                    "omega={omega}, cells={elements}, error={error}, previous={previous}"
                );
            }
            previous = error;
            near(thermal.mean, expected_h, 0.003);
            near(compliance(wave, omega), expected_compliance, 3e-4);
            near(thermal.mean + thermal.wall, C64::ONE, 3e-16);
            assert!(thermal.gradient_work > 0.0 && thermal.relative_work_defect < 1e-8);
        }
    }
    let section = AnnularSection::new(a, b, 128).unwrap();
    let molecular = section
        .evaluate_molecular(&state, 500.0, admission())
        .unwrap();
    let equilibrium = section.evaluate(gas, 500.0, admission()).unwrap();
    assert!(
        (molecular.wave.wavenumber - equilibrium.wave.wavenumber).abs()
            > 1e-3 * molecular.wave.wavenumber.abs()
    );
}

#[test]
fn narrow_annuli_match_the_complex_caloric_slab_solution() {
    let state = material(2, 0.01, 101325.0);
    let gas = state.gas();
    // Small outer radius keeps the fully viscous wave below the declared
    // transverse-phase ceiling even at the highest frequency.
    let (a, b) = (0.001, 0.00102);
    let section = AnnularSection::new(a, b, 256).unwrap();
    for omega in [1000.0, 10000.0, 50000.0] {
        let cp = heat_capacity(&state, omega);
        let q = (C64::new(0.0, omega) * cp)
            .scale(gas.density / gas.thermal_conductivity)
            .sqrt()
            .scale(0.5 * (b - a));
        let tanh = C64::new((2.0 * q.re).sinh(), (2.0 * q.im).sin())
            .scale(1.0 / ((2.0 * q.re).cosh() + (2.0 * q.im).cos()));
        let h = C64::ONE - tanh / q;
        let expected_compliance = (C64::ONE - h / cp.scale(1.0 / gas.specific_gas_constant))
            .scale(1.0 / gas.pressure)
            .conj();
        let wave = section
            .evaluate_molecular(&state, omega, admission())
            .unwrap();
        near(wave.thermal.unwrap().mean, h, 2e-5);
        near(compliance(wave, omega), expected_compliance, 1e-6);
    }
}

#[test]
fn isothermal_equilibrated_and_frozen_limits_keep_the_original_storage() {
    let fast = material(2, 1e-9, 101325.0);
    let slow = material(2, 1e9, 101325.0);
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    for state in [&fast, &slow] {
        near(
            compliance(
                section
                    .evaluate_molecular(state, 0.001, admission())
                    .unwrap(),
                0.001,
            ),
            C64::from_re(1.0 / state.gas().pressure),
            2e-6,
        );
    }
    let gas = slow.gas();
    let transport = gas.thermoviscous_state().unwrap().transport();
    let cv = slow.frozen_cv_j_kg_k();
    let frozen = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(gas.density, slow.frozen_sound_speed_m_s()).unwrap(),
        FluidTransportState::try_new(
            gas.density,
            cv + gas.specific_gas_constant,
            transport.thermal_conductivity_w_m_k(),
            transport.dynamic_viscosity_pa_s(),
        )
        .unwrap(),
        cv,
    )
    .unwrap();
    for omega in [500.0, 3000.0] {
        let equilibrated = section
            .evaluate_molecular(&fast, omega, admission())
            .unwrap();
        let expected_eq = section.evaluate(fast.gas(), omega, admission()).unwrap();
        let unrelaxed = section
            .evaluate_molecular(&slow, omega, admission())
            .unwrap();
        let expected_frozen = section.evaluate(&frozen, omega, admission()).unwrap();
        near(
            equilibrated.wave.wavenumber,
            expected_eq.wave.wavenumber,
            2e-8,
        );
        near(
            unrelaxed.wave.wavenumber,
            expected_frozen.wave.wavenumber,
            2e-8,
        );
        near(
            unrelaxed.wave.characteristic_impedance,
            expected_frozen.wave.characteristic_impedance,
            2e-8,
        );
        assert!(
            (equilibrated.wave.wavenumber - unrelaxed.wave.wavenumber).abs()
                > 1e-3 * unrelaxed.wave.wavenumber.abs()
        );
    }
}

#[test]
fn original_pressure_dependent_rates_preserve_similarity_and_passive_scattering() {
    let mut reference = None;
    for scale in [1.0_f64, 0.5, 2.0] {
        let state = material(2, 1.0, 101325.0 * scale);
        let section = AnnularSection::new(0.0004 / scale, 0.0012 / scale, 128).unwrap();
        let wave = section
            .evaluate_molecular(&state, 500.0 * scale, admission())
            .unwrap();
        let block = wave
            .two_port(
                0.2 / scale,
                state.gas().characteristic_impedance / section.area_m2(),
            )
            .unwrap();
        let k = wave.wave.wavenumber.scale(1.0 / scale);
        let z = wave
            .wave
            .characteristic_impedance
            .scale(1.0 / scale.powi(3));
        if let Some((old_k, old_z, old_s)) = reference {
            near(k, old_k, 3e-12);
            near(z, old_z, 3e-12);
            for (actual, expected) in block.scattering.into_iter().zip(old_s) {
                near(actual, expected, 3e-12);
            }
        } else {
            reference = Some((k, z, block.scattering));
        }
        let [r, t, reverse, opposite] = block.scattering;
        assert_eq!(r, opposite);
        assert_eq!(t, reverse);
        for phase in [0.0_f64, 0.7, 1.9, core::f64::consts::PI] {
            let other = C64::new(phase.cos(), phase.sin()).scale(0.7);
            let outgoing = (r + t * other).norm_sq() + (reverse + opposite * other).norm_sq();
            assert!(outgoing < 1.0 + other.norm_sq());
        }
        // Sixty physical amplitude e-foldings, for every pressure-scaled gas.
        let opaque = wave
            .two_port(
                60.0 / wave.wave.wavenumber.im,
                state.gas().characteristic_impedance / section.area_m2(),
            )
            .unwrap();
        assert!(
            opaque
                .scattering
                .iter()
                .all(|v| v.re.is_finite() && v.im.is_finite())
        );
        assert!(opaque.scattering[1].abs() < 1e-20);
    }
}

#[test]
fn complex_thermal_resolution_payload_and_cancellation_refuse_without_mutation() {
    let state = material(2, 1.0, 101325.0);
    let before = state.clone();
    let section = AnnularSection::new(0.0004, 0.0012, 128).unwrap();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e6] {
        assert!(
            section
                .evaluate_molecular(&state, omega, admission())
                .is_err()
        );
    }
    let exact_bytes = 192 * (section.elements() + 1) + 4096;
    let controls = AnnularAdmission {
        maximum_payload_bytes: exact_bytes,
        ..admission()
    };
    let original = section.evaluate_molecular(&state, 500.0, controls).unwrap();
    assert!(
        section
            .evaluate_molecular(
                &state,
                500.0,
                AnnularAdmission {
                    maximum_payload_bytes: exact_bytes - 1,
                    ..controls
                }
            )
            .is_err()
    );
    assert!(
        section
            .evaluate_molecular(
                &state,
                500.0,
                AnnularAdmission {
                    maximum_cell_diffusion_phase: 1e-9,
                    ..controls
                }
            )
            .is_err()
    );
    assert!(
        section
            .evaluate_molecular(
                &state,
                500.0,
                AnnularAdmission {
                    maximum_relative_residual: f64::NAN,
                    ..controls
                }
            )
            .is_err()
    );
    let mut polls = 0;
    let refused = section.evaluate_molecular_controlled(&state, 500.0, controls, || {
        polls += 1;
        polls == 16 // after the complete viscous field, during molecular thermal work
    });
    assert!(matches!(
        refused,
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert_eq!(state, before);
    let replay = section.evaluate_molecular(&state, 500.0, controls).unwrap();
    assert_eq!(replay.wave.wavenumber, original.wave.wavenumber);
    assert_eq!(
        replay.wave.characteristic_impedance,
        original.wave.characteristic_impedance
    );
}
