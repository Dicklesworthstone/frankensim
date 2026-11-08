//! G0/G1/G3: finite thermal strain, entropy and adiabatic waves share one potential.

#[path = "hyper_thermoelastic/source.rs"]
mod source;

use fs_material::hyper::acoustic::HyperelasticAcousticState;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::thermoelastic::{FiniteThermoelasticState, ThermalReferenceJet};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};

const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const G: [f64; 9] = [1.09, 0.07, -0.025, 0.03, 0.97, 0.04, -0.015, 0.02, 1.04];
const G1: [f64; 9] = [
    0.002, 0.001, -0.0004, -0.0003, -0.001, 0.0006, 0.0002, 0.0008, 0.0015,
];
const G2: [f64; 9] = [2e-5, -3e-5, 1e-5, 1e-5, 4e-5, 0.0, -2e-5, 1e-5, -1e-5];
const F: [f64; 9] = [1.22, 0.16, 0.04, -0.05, 0.91, 0.12, 0.03, -0.06, 1.13];

fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 4.0,
            c01: 2.0,
            kappa: 40.0,
        },
        4.0,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let tolerance = absolute + relative * expected.abs();
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} vs {expected:e}, allowance {tolerance:e}"
    );
}

fn diagonal(value: f64) -> [f64; 9] {
    [value, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0, value]
}

fn multiply(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn cof(a: [f64; 9]) -> [f64; 9] {
    [
        a[4] * a[8] - a[5] * a[7],
        a[5] * a[6] - a[3] * a[8],
        a[3] * a[7] - a[4] * a[6],
        a[2] * a[7] - a[1] * a[8],
        a[0] * a[8] - a[2] * a[6],
        a[1] * a[6] - a[0] * a[7],
        a[1] * a[5] - a[2] * a[4],
        a[2] * a[3] - a[0] * a[5],
        a[0] * a[4] - a[1] * a[3],
    ]
}

fn determinant(a: [f64; 9]) -> f64 {
    let c = cof(a);
    a[0] * c[0] + a[1] * c[1] + a[2] * c[2]
}

// Independent invariant expression and explicit matrix inverse, without the
// native thermal chain rule or reference-deformation owner.
fn energy(f: [f64; 9], g: [f64; 9]) -> f64 {
    let jg = determinant(g);
    let cofactor = cof(g);
    let inverse = core::array::from_fn(|ij| cofactor[3 * (ij % 3) + ij / 3] / jg);
    let fe = multiply(f, inverse);
    let je = determinant(fe);
    let i1: f64 = fe.iter().map(|x| x * x).sum();
    let i2: f64 = cof(fe).iter().map(|x| x * x).sum();
    let iso = je.powf(-2.0 / 3.0);
    jg * (4.0 * (iso * i1 - 3.0) + 2.0 * (iso * iso * i2 - 3.0) + 20.0 * (je - 1.0).powi(2))
}

fn reference_at(delta_k: f64) -> ReferenceDeformation {
    ReferenceDeformation::try_new(core::array::from_fn(|i| {
        G[i] + delta_k * G1[i] + 0.5 * delta_k * delta_k * G2[i]
    }))
    .unwrap()
}

fn jet() -> ThermalReferenceJet {
    ThermalReferenceJet::try_new(reference_at(0.0), G1, G2).unwrap()
}

fn derivative(values: [f64; 4], step: f64) -> f64 {
    (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * step)
}

#[test]
fn g1_isotropic_reference_entropy_capacity_and_waves_match_closed_forms() {
    let material = law();
    let (s, a, b, temperature, rho0, caloric) = (1.1_f64, 0.01, -0.0003, 300.0, 6.0, 100.0);
    let reference = ReferenceDeformation::try_new(diagonal(s)).unwrap();
    let thermal = ThermalReferenceJet::try_new(reference, diagonal(a), diagonal(b)).unwrap();
    for f in [diagonal(s), F] {
        let state = FiniteThermoelasticState::try_new(&material, &f, temperature, thermal).unwrap();
        let q = determinant(f);
        let iso = q.powf(-2.0 / 3.0);
        let i1: f64 = f.iter().map(|x| x * x).sum();
        let i2: f64 = cof(f).iter().map(|x| x * x).sum();
        let dev = 4.0 * (iso * i1 - 3.0) + 2.0 * (iso * iso * i2 - 3.0);
        // Psi(s)=B*s^3+C/s^3-40*q, at fixed TOTAL F. This also
        // checks G'' through the psi_s*b term away from the natural state.
        let (big_b, big_c) = (dev + 20.0, 20.0 * q * q);
        let psi = big_b * s.powi(3) + big_c / s.powi(3) - 40.0 * q;
        let psi_s = 3.0 * big_b * s * s - 3.0 * big_c / s.powi(4);
        let psi_ss = 6.0 * big_b * s + 12.0 * big_c / s.powi(5);
        let entropy = -psi_s * a;
        near(state.free_energy_j_m3(), psi, 8e-14, 2e-13);
        near(state.entropy_j_m3_k(), entropy, 2e-14, 2e-12);
        near(
            state.internal_energy_j_m3(),
            psi + temperature * entropy,
            4e-12,
            2e-12,
        );
        near(
            state.mechanical_heat_capacity_j_m3_k(),
            -temperature * (psi_ss * a * a + psi_s * b),
            2e-12,
            2e-12,
        );
    }
    let state =
        FiniteThermoelasticState::try_new(&material, &diagonal(s), temperature, thermal).unwrap();
    let mechanical_capacity = -9.0 * temperature * 40.0 * s * a * a;
    let total_capacity = caloric + mechanical_capacity;
    let beta = -3.0 * 40.0 * s * a;
    for i in 0..9 {
        near(
            state.stress_temperature_derivative_pa_k()[i],
            if i % 4 == 0 { beta } else { 0.0 },
            5e-13,
            2e-12,
        );
    }
    let adiabatic = state.adiabatic_acoustics(rho0, caloric).unwrap();
    near(
        adiabatic.fixed_deformation_heat_capacity_j_m3_k(),
        total_capacity,
        3e-12,
        2e-12,
    );
    // The spatial shear modulus is 2*(c10+c01)=12. Thermal coupling
    // adds only to the longitudinal modulus at this isotropic natural state.
    let rho = rho0 / s.powi(3);
    let correction = 9.0 * temperature * 40.0_f64.powi(2) * s * a * a / total_capacity;
    for n in [[1.0, 0.0, 0.0], [0.3, -0.4, 0.5]] {
        let waves = adiabatic.acoustic().directional_waves(n).unwrap();
        near(waves.current_density_kg_m3, rho, 2e-14, 2e-14);
        for (mode, modulus) in waves.modes.iter().zip([12.0, 12.0, 56.0 + correction]) {
            near(mode.directional_modulus_pa, modulus, 2e-11, 2e-12);
            near(mode.squared_phase_speed_m2_s2, modulus / rho, 2e-11, 2e-12);
            let velocity = mode.energy_velocity_m_s.unwrap();
            near(
                (0..3).map(|i| waves.direction[i] * velocity[i]).sum(),
                mode.phase_speed_m_s.unwrap(),
                2e-11,
                2e-12,
            );
        }
    }
}

#[test]
fn g1_noncommuting_reference_jets_match_energy_stress_and_maxwell_derivatives() {
    let material = law();
    let temperature = 310.0;
    let state = FiniteThermoelasticState::try_new(&material, &F, temperature, jet()).unwrap();
    assert!(
        multiply(G, G1)
            .iter()
            .zip(multiply(G1, G))
            .any(|(a, b)| (a - b).abs() > 1e-5)
    );
    assert_eq!(state.deformation_gradient(), &F);
    near(state.free_energy_j_m3(), energy(F, G), 5e-14, 2e-13);
    let h = 0.05;
    let shifts = [-2.0 * h, -h, h, 2.0 * h];
    let energies = shifts.map(|delta| energy(F, *reference_at(delta).deformation_gradient()));
    let psi_t = derivative(energies, h);
    let psi_tt = (-energies[3] + 16.0 * energies[2] - 30.0 * energy(F, G) + 16.0 * energies[1]
        - energies[0])
        / (12.0 * h * h);
    near(state.entropy_j_m3_k(), -psi_t, 2e-10, 3e-8);
    near(
        state.mechanical_heat_capacity_j_m3_k(),
        -temperature * psi_tt,
        8e-8,
        3e-7,
    );
    let stresses = shifts.map(|delta| {
        ReferencedHyperelastic::new(&material, reference_at(delta))
            .first_piola(&F)
            .unwrap()
    });
    for i in 0..9 {
        let beta = state.stress_temperature_derivative_pa_k()[i];
        near(beta, derivative(stresses.map(|p| p[i]), h), 5e-10, 5e-8);
        // Maxwell symmetry dP/dT = -dS/dF is checked through independent
        // F perturbations, not by differentiating the returned Piola array.
        let epsilon = 2e-5;
        let entropy = [-2.0 * epsilon, -epsilon, epsilon, 2.0 * epsilon].map(|delta| {
            let mut f = F;
            f[i] += delta;
            FiniteThermoelasticState::try_new(&material, &f, temperature, jet())
                .unwrap()
                .entropy_j_m3_k()
        });
        near(beta, -derivative(entropy, epsilon), 2e-9, 5e-8);
    }
    let affine = ThermalReferenceJet::try_new(reference_at(0.0), G1, [0.0; 9]).unwrap();
    let affine_state =
        FiniteThermoelasticState::try_new(&material, &F, temperature, affine).unwrap();
    assert!(
        (state.mechanical_heat_capacity_j_m3_k() - affine_state.mechanical_heat_capacity_j_m3_k())
            .abs()
            > 1e-3
    );
}

#[test]
fn g3_strained_adiabatic_waves_keep_full_coupling_objectivity_and_current_density() {
    let material = law();
    let (temperature, rho0, caloric) = (310.0, 6.0, 50.0);
    let state = FiniteThermoelasticState::try_new(&material, &F, temperature, jet()).unwrap();
    let adiabatic = state.adiabatic_acoustics(rho0, caloric).unwrap();
    let iso = HyperelasticAcousticState::from_material(
        &ReferencedHyperelastic::new(&material, reference_at(0.0)),
        &F,
        rho0,
    )
    .unwrap();
    let beta = state.stress_temperature_derivative_pa_k();
    assert!([1, 2, 3, 5, 6, 7].iter().any(|&i| beta[i].abs() > 1e-4));
    let scale = temperature / adiabatic.fixed_deformation_heat_capacity_j_m3_k();
    for i in 0..9 {
        near(
            adiabatic.temperature_gradient_k()[i],
            scale * beta[i],
            3e-14,
            2e-13,
        );
        for j in 0..9 {
            near(
                adiabatic.isentropic_tangent_pa()[i][j] - state.isothermal_tangent_pa()[i][j],
                scale * beta[i] * beta[j],
                3e-13,
                2e-11,
            );
        }
    }
    let direction = [0.3, -0.4, 0.5];
    let waves = adiabatic.acoustic().directional_waves(direction).unwrap();
    let iso_waves = iso.directional_waves(direction).unwrap();
    let q: [f64; 3] =
        core::array::from_fn(|j| (0..3).map(|i| F[3 * i + j] * waves.direction[i]).sum());
    let thermal_force: [f64; 3] =
        core::array::from_fn(|i| (0..3).map(|j| beta[3 * i + j] * q[j]).sum());
    for i in 0..3 {
        for j in 0..3 {
            near(
                waves.acoustic_tensor_pa[i][j] - iso_waves.acoustic_tensor_pa[i][j],
                scale * thermal_force[i] * thermal_force[j] / determinant(F),
                8e-13,
                2e-10,
            );
        }
    }
    let rotation = [0.6, -0.8, 0.0, 0.8, 0.6, 0.0, 0.0, 0.0, 1.0];
    let rotated_f = multiply(rotation, F);
    let rotated =
        FiniteThermoelasticState::try_new(&material, &rotated_f, temperature, jet()).unwrap();
    for (a, b) in [
        (rotated.free_energy_j_m3(), state.free_energy_j_m3()),
        (rotated.entropy_j_m3_k(), state.entropy_j_m3_k()),
        (rotated.internal_energy_j_m3(), state.internal_energy_j_m3()),
        (
            rotated.mechanical_heat_capacity_j_m3_k(),
            state.mechanical_heat_capacity_j_m3_k(),
        ),
    ] {
        near(a, b, 3e-12, 3e-12);
    }
    for (actual, expected) in rotated
        .stress_temperature_derivative_pa_k()
        .iter()
        .zip(multiply(rotation, *beta))
    {
        near(*actual, expected, 3e-13, 3e-12);
    }
    let rotated_direction =
        core::array::from_fn(|i| (0..3).map(|j| rotation[3 * i + j] * direction[j]).sum());
    let next_waves = rotated
        .adiabatic_acoustics(rho0, caloric)
        .unwrap()
        .acoustic()
        .directional_waves(rotated_direction)
        .unwrap();
    near(
        next_waves.current_density_kg_m3,
        rho0 / determinant(F),
        2e-14,
        2e-14,
    );
    for i in 0..3 {
        for j in 0..3 {
            let expected: f64 = (0..3)
                .flat_map(|a| (0..3).map(move |b| (a, b)))
                .map(|(a, b)| {
                    rotation[3 * i + a] * waves.acoustic_tensor_pa[a][b] * rotation[3 * j + b]
                })
                .sum();
            near(next_waves.acoustic_tensor_pa[i][j], expected, 2e-11, 2e-12);
        }
        near(
            next_waves.modes[i].squared_phase_speed_m2_s2,
            waves.modes[i].squared_phase_speed_m2_s2,
            2e-11,
            2e-12,
        );
        let old_velocity = waves.modes[i].energy_velocity_m_s.unwrap();
        let new_velocity = next_waves.modes[i].energy_velocity_m_s.unwrap();
        for a in 0..3 {
            near(
                new_velocity[a],
                (0..3).map(|b| rotation[3 * a + b] * old_velocity[b]).sum(),
                3e-11,
                3e-12,
            );
        }
    }
}

#[test]
fn g0_zero_thermal_jet_preserves_exact_isothermal_response_and_signed_instability() {
    let material = law();
    for (g, f) in [(G, F), (diagonal(5.0), diagonal(5.0))] {
        let reference = ReferenceDeformation::try_new(g).unwrap();
        let thermal = ThermalReferenceJet::try_new(reference, [0.0; 9], [0.0; 9]).unwrap();
        let state = FiniteThermoelasticState::try_new(&material, &f, 300.0, thermal).unwrap();
        assert_eq!(state.entropy_j_m3_k(), 0.0);
        assert_eq!(state.mechanical_heat_capacity_j_m3_k(), 0.0);
        assert_eq!(state.stress_temperature_derivative_pa_k(), &[0.0; 9]);
        let adiabatic = state.adiabatic_acoustics(6.0, 50.0).unwrap();
        let expected = HyperelasticAcousticState::from_material(
            &ReferencedHyperelastic::new(&material, reference),
            &f,
            6.0,
        )
        .unwrap();
        assert_eq!(
            adiabatic.isentropic_tangent_pa(),
            state.isothermal_tangent_pa()
        );
        assert_eq!(adiabatic.acoustic(), &expected);
        assert_eq!(
            adiabatic
                .acoustic()
                .directional_waves([0.3, 0.4, -0.5])
                .unwrap(),
            expected.directional_waves([0.3, 0.4, -0.5]).unwrap()
        );
    }
    let unstable_law = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 1.0,
            lambda: 8.0,
        },
        4.0,
    )
    .unwrap();
    let zero =
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), [0.0; 9], [0.0; 9]).unwrap();
    let state =
        FiniteThermoelasticState::try_new(&unstable_law, &diagonal(3.0), 300.0, zero).unwrap();
    let waves = state
        .adiabatic_acoustics(6.0, 50.0)
        .unwrap()
        .acoustic()
        .directional_waves([1.0, 0.0, 0.0])
        .unwrap();
    assert!(waves.modes[0].squared_phase_speed_m2_s2 < 0.0);
    assert_eq!(waves.modes[0].phase_speed_m_s, None);
    assert_eq!(waves.modes[0].energy_velocity_m_s, None);
}

#[test]
fn g0_thermal_and_caloric_domains_refuse_without_mutating_the_snapshot() {
    let material = law();
    assert!(
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), [f64::NAN; 9], [0.0; 9])
            .is_err()
    );
    assert!(
        ThermalReferenceJet::try_new(
            ReferenceDeformation::identity(),
            [0.0; 9],
            [f64::INFINITY; 9]
        )
        .is_err()
    );
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(FiniteThermoelasticState::try_new(&material, &F, temperature, jet()).is_err());
    }
    for f in [
        [0.0; 9],
        [f64::NAN; 9],
        [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        diagonal(10.0),
    ] {
        assert!(FiniteThermoelasticState::try_new(&material, &f, 300.0, jet()).is_err());
    }
    let thermal =
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), diagonal(0.02), [0.0; 9])
            .unwrap();
    let state = FiniteThermoelasticState::try_new(&material, &I, 300.0, thermal).unwrap();
    let before = state.clone();
    let cancellation_capacity = -state.mechanical_heat_capacity_j_m3_k();
    assert!(cancellation_capacity > 0.0);
    for capacity in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        cancellation_capacity,
        0.5 * cancellation_capacity,
    ] {
        assert!(state.adiabatic_acoustics(6.0, capacity).is_err());
    }
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            state
                .adiabatic_acoustics(density, 2.0 * cancellation_capacity)
                .is_err()
        );
    }
    assert!(
        state
            .adiabatic_acoustics(6.0, 2.0 * cancellation_capacity)
            .is_ok()
    );
    assert_eq!(state, before);
    let mut invalid_law = material.clone();
    invalid_law.model = HyperelasticModel::NeoHookean {
        mu: -1.0,
        lambda: 1.0,
    };
    assert!(FiniteThermoelasticState::try_new(&invalid_law, &I, 300.0, thermal).is_err());
    assert!(
        ThermalReferenceJet::try_new(ReferenceDeformation::identity(), [1e300; 9], [0.0; 9])
            .and_then(|jet| FiniteThermoelasticState::try_new(&material, &I, 300.0, jet))
            .is_err()
    );
}

#[test]
fn simultaneous_deformation_and_temperature_rates_match_full_path_differences() {
    let material = law();
    let temperature = 310.0;
    let fdot = [0.03, -0.04, 0.02, -0.01, 0.05, 0.03, -0.02, 0.01, -0.04];
    let tdot = 1.7;
    let state = FiniteThermoelasticState::try_new(&material, &F, temperature, jet()).unwrap();
    let rate = state.rate(fdot, tdot).unwrap();
    assert_eq!(rate.deformation_gradient_rate_per_s(), &fdot);
    assert_eq!(rate.temperature_rate_k_s(), tdot);
    let h = 0.002;
    let samples = [-2.0 * h, -h, h, 2.0 * h].map(|time| {
        let delta = tdot * time;
        let f = core::array::from_fn(|i| F[i] + time * fdot[i]);
        let thermal = ThermalReferenceJet::try_new(
            reference_at(delta),
            core::array::from_fn(|i| G1[i] + delta * G2[i]),
            G2,
        )
        .unwrap();
        FiniteThermoelasticState::try_new(&material, &f, temperature + delta, thermal).unwrap()
    });
    for i in 0..9 {
        near(
            rate.first_piola_rate_pa_s()[i],
            derivative(samples.each_ref().map(|s| s.first_piola_pa()[i]), h),
            3e-8,
            2e-7,
        );
        near(
            rate.elastic_deformation_gradient_rate_per_s()[i],
            derivative(
                samples
                    .each_ref()
                    .map(|s| s.elastic_deformation_gradient()[i]),
                h,
            ),
            2e-11,
            2e-8,
        );
    }
    for (actual, values) in [
        (
            rate.free_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.free_energy_j_m3()),
        ),
        (
            rate.entropy_rate_w_m3_k(),
            samples.each_ref().map(|s| s.entropy_j_m3_k()),
        ),
        (
            rate.internal_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.internal_energy_j_m3()),
        ),
    ] {
        near(actual, derivative(values, h), 3e-8, 2e-7);
    }
    near(
        rate.internal_energy_rate_w_m3(),
        rate.free_energy_rate_w_m3()
            + temperature * rate.entropy_rate_w_m3_k()
            + state.entropy_j_m3_k() * tdot,
        2e-13,
        2e-13,
    );
}

#[test]
fn thermal_rate_tracks_changing_observers_without_objectivizing_piola_by_fiat() {
    let material = law();
    let state = FiniteThermoelasticState::try_new(&material, &F, 310.0, jet()).unwrap();
    let fdot = [0.03, -0.04, 0.02, -0.01, 0.05, 0.03, -0.02, 0.01, -0.04];
    let tdot = -2.1;
    let rate = state.rate(fdot, tdot).unwrap();
    let q = [0.6, -0.8, 0.0, 0.8, 0.6, 0.0, 0.0, 0.0, 1.0];
    let spin = [0.0, -0.7, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0];
    let rotated_f = multiply(q, F);
    let q_fdot = multiply(q, fdot);
    let spin_f = multiply(spin, rotated_f);
    let rotated_fdot = core::array::from_fn(|i| q_fdot[i] + spin_f[i]);
    let rotated = FiniteThermoelasticState::try_new(&material, &rotated_f, 310.0, jet()).unwrap();
    let observed = rotated.rate(rotated_fdot, tdot).unwrap();
    let q_pdot = multiply(q, *rate.first_piola_rate_pa_s());
    let spin_p = multiply(spin, multiply(q, *state.first_piola_pa()));
    let q_fedot = multiply(q, *rate.elastic_deformation_gradient_rate_per_s());
    let spin_fe = multiply(spin, multiply(q, *state.elastic_deformation_gradient()));
    for i in 0..9 {
        near(
            observed.first_piola_rate_pa_s()[i],
            q_pdot[i] + spin_p[i],
            2e-11,
            2e-11,
        );
        near(
            observed.elastic_deformation_gradient_rate_per_s()[i],
            q_fedot[i] + spin_fe[i],
            2e-13,
            2e-12,
        );
    }
    for (actual, expected) in [
        (
            observed.free_energy_rate_w_m3(),
            rate.free_energy_rate_w_m3(),
        ),
        (observed.entropy_rate_w_m3_k(), rate.entropy_rate_w_m3_k()),
        (
            observed.internal_energy_rate_w_m3(),
            rate.internal_energy_rate_w_m3(),
        ),
    ] {
        near(actual, expected, 2e-11, 2e-11);
    }
}

#[test]
fn thermal_free_expansion_has_no_elastic_rate_and_rate_failures_preserve_state() {
    let material = law();
    let (stretch, expansion, tdot) = (1.1, 0.002, 3.0);
    let thermal = ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(diagonal(stretch)).unwrap(),
        diagonal(expansion),
        [0.0; 9],
    )
    .unwrap();
    let state =
        FiniteThermoelasticState::try_new(&material, &diagonal(stretch), 300.0, thermal).unwrap();
    let rate = state.rate(diagonal(expansion * tdot), tdot).unwrap();
    for actual in rate
        .first_piola_rate_pa_s()
        .iter()
        .chain(rate.elastic_deformation_gradient_rate_per_s())
    {
        near(*actual, 0.0, 2e-12, 0.0);
    }
    for actual in [
        rate.free_energy_rate_w_m3(),
        rate.entropy_rate_w_m3_k(),
        rate.internal_energy_rate_w_m3(),
    ] {
        near(actual, 0.0, 2e-11, 0.0);
    }
    let before = state.clone();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(state.rate([bad; 9], 0.0).is_err());
        assert!(state.rate([0.0; 9], bad).is_err());
    }
    assert!(state.rate([f64::MAX; 9], 0.0).is_err());
    let zero = state.rate([0.0; 9], 0.0).unwrap();
    assert_eq!(zero.first_piola_rate_pa_s(), &[0.0; 9]);
    assert_eq!(zero.elastic_deformation_gradient_rate_per_s(), &[0.0; 9]);
    assert_eq!(zero.free_energy_rate_w_m3(), 0.0);
    assert_eq!(zero.entropy_rate_w_m3_k(), 0.0);
    assert_eq!(zero.internal_energy_rate_w_m3(), 0.0);
    assert_eq!(state, before);
}
