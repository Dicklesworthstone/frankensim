//! Reversible temperature derivatives retain the native finite-strain memory.

use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::incompressible_plane_stress::{
    IncompressiblePlaneStressPotential, incompressible_plane_stress_with_kinematics,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::thermoelastic::incompressible_plane_stress::IncompressibleThermoelasticPlaneStressState;
use fs_material::hyper::thermoelastic::{FiniteThermoelasticState, ThermalReferenceJet};
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm};
use fs_material::visco_hyper::{
    BranchEnergy, FiniteStrainBranch, ViscoHyperelastic, ViscoRemainderState, ViscoState,
};

const G: [f64; 9] = [1.09, 0.07, -0.025, 0.03, 0.97, 0.04, -0.015, 0.02, 1.04];
const G1: [f64; 9] = [
    0.002, 0.001, -0.0004, -0.0003, -0.001, 0.0006, 0.0002, 0.0008, 0.0015,
];
const G2: [f64; 9] = [2e-5, -3e-5, 1e-5, 1e-5, 4e-5, 0.0, -2e-5, 1e-5, -1e-5];
const F: [f64; 9] = [1.22, 0.16, 0.04, -0.05, 0.91, 0.12, 0.03, -0.06, 1.13];
const T: f64 = 310.0;

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let tolerance = absolute + relative * expected.abs();
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} vs {expected:e}, allowance {tolerance:e}"
    );
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn cofactor(a: &[f64; 9]) -> [f64; 9] {
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

fn determinant(a: &[f64; 9]) -> f64 {
    let c = cofactor(a);
    a[0] * c[0] + a[1] * c[1] + a[2] * c[2]
}

fn jet(delta: f64) -> ThermalReferenceJet {
    let g = core::array::from_fn(|i| G[i] + delta * G1[i] + 0.5 * delta * delta * G2[i]);
    let gt = core::array::from_fn(|i| G1[i] + delta * G2[i]);
    ThermalReferenceJet::try_new(ReferenceDeformation::try_new(g).unwrap(), gt, G2).unwrap()
}

fn law(ogden: bool) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 8.0,
                lambda: 12.0,
            },
            4.0,
        )
        .unwrap(),
        vec![
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 14.0 },
                tau_s: 0.4,
            },
            FiniteStrainBranch {
                energy: if ogden {
                    BranchEnergy::Ogden {
                        terms: vec![
                            OgdenTerm {
                                mu_pa: 5.0,
                                alpha: 1.3,
                            },
                            OgdenTerm {
                                mu_pa: -1.0,
                                alpha: -2.0,
                            },
                        ],
                    }
                } else {
                    BranchEnergy::NeoHookean { mu_pa: 5.0 }
                },
                tau_s: 1.3,
            },
        ],
    )
    .unwrap()
    .with_reference_deformation(*jet(0.0).reference())
}

fn history(law: &ViscoHyperelastic) -> ViscoState {
    let f1 = [1.3, 0.14, 0.02, -0.03, 0.91, 0.08, 0.01, -0.04, 1.08];
    let f2 = [1.16, -0.08, 0.11, 0.07, 1.11, -0.04, 0.025, 0.05, 0.94];
    let (first, _) = law.step(&law.initial_state(), &f1, 0.07).unwrap();
    law.step(&first, &f2, 0.11).unwrap().0
}

fn shifted(
    original: &ViscoHyperelastic,
    accepted: &ViscoState,
    delta: f64,
) -> (ViscoHyperelastic, ViscoState) {
    let changed = original
        .clone()
        .with_reference_deformation(*jet(delta).reference());
    let state = changed
        .change_reference_from(original, accepted)
        .unwrap()
        .state;
    assert_eq!(
        state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    assert_eq!(state.dissipated, accepted.dissipated);
    (changed, state)
}

// Independent invariant energy, without native thermal derivatives, AD,
// eigensystems or reference pullbacks. Both branches here are Neo-Hookean.
fn independent_energy(f: &[f64; 9], accepted: &ViscoState, delta: f64) -> f64 {
    let g = core::array::from_fn(|i| G[i] + delta * G1[i] + 0.5 * delta * delta * G2[i]);
    let jg = determinant(&g);
    let cof = cofactor(&g);
    let inverse = core::array::from_fn(|ij| cof[3 * (ij % 3) + ij / 3] / jg);
    let fe = multiply(f, &inverse);
    let je = determinant(&fe);
    let log_j = je.ln();
    let i1: f64 = fe.iter().map(|entry| entry * entry).sum();
    let mut energy = 4.0 * (i1 - 3.0) - 8.0 * log_j + 6.0 * log_j * log_j;
    for (cv, mu) in accepted
        .inverse_viscous_cauchy_green
        .iter()
        .zip([14.0, 5.0])
    {
        let metric = [
            cv[0], cv[3], cv[4], cv[3], cv[1], cv[5], cv[4], cv[5], cv[2],
        ];
        let fe_cv = multiply(&fe, &metric);
        let trace: f64 = fe_cv.iter().zip(fe).map(|(left, right)| left * right).sum();
        energy += 0.5 * mu * (je.powf(-2.0 / 3.0) * trace - 3.0);
    }
    jg * energy
}

#[test]
fn finite_history_entropy_internal_energy_and_capacity_follow_original_potential() {
    let law = law(false);
    let accepted = history(&law);
    let before = accepted.clone();
    assert!(
        accepted
            .inverse_viscous_cauchy_green
            .iter()
            .any(|cv| cv[3..].iter().any(|x| x.abs() > 1e-4))
    );
    assert!(accepted.dissipated > 0.0);
    let state = law
        .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
        .unwrap();
    let h = 0.04;
    let sample: [f64; 5] =
        core::array::from_fn(|i| independent_energy(&F, &accepted, (i as f64 - 2.0) * h));
    let psi_t = (sample[0] - 8.0 * sample[1] + 8.0 * sample[3] - sample[4]) / (12.0 * h);
    let psi_tt = (-sample[4] + 16.0 * sample[3] - 30.0 * sample[2] + 16.0 * sample[1] - sample[0])
        / (12.0 * h * h);
    near(state.free_energy_j_m3(), sample[2], 2e-12, 2e-12);
    near(state.entropy_j_m3_k(), -psi_t, 2e-10, 2e-8);
    near(
        state.internal_energy_j_m3(),
        sample[2] - T * psi_t,
        8e-8,
        2e-8,
    );
    near(
        state.mechanical_heat_capacity_j_m3_k(),
        -T * psi_tt,
        2e-6,
        2e-6,
    );
    let reset = law
        .frozen_thermoelastic_state(&law.initial_state(), &F, T, jet(0.0))
        .unwrap();
    assert!((state.entropy_j_m3_k() - reset.entropy_j_m3_k()).abs() > 1e-5);
    assert_eq!(accepted, before);
}

#[test]
fn native_ogden_response_is_instantaneous_and_independent_of_relaxation_times() {
    let law = law(true);
    let accepted = history(&law);
    let before = accepted.clone();
    let state = law
        .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
        .unwrap();
    let frozen = law.freeze(&accepted, 0.0).unwrap();
    let (p, energy) = frozen.stress_and_potential(&F).unwrap();
    assert_eq!(state.first_piola_pa(), &p);
    assert_eq!(state.free_energy_j_m3(), energy);
    assert_eq!(
        state.isothermal_tangent_pa(),
        &frozen.consistent_tangent(&F).unwrap()
    );
    near(
        energy,
        law.free_energy(&accepted, &F).unwrap(),
        2e-12,
        2e-12,
    );
    let mut shifted_times = law.clone();
    shifted_times.branches[0].tau_s *= 1e-8;
    shifted_times.branches[1].tau_s *= 1e8;
    assert_eq!(
        state,
        shifted_times
            .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
            .unwrap()
    );
    let interval = law.freeze(&accepted, 0.2).unwrap().energy_density(&F);
    assert!((interval - energy).abs() > 1e-3);
    assert_eq!(accepted, before);
}

#[test]
fn all_mixed_thermal_stress_components_obey_maxwell_relation_and_objectivity() {
    let law = law(true);
    let accepted = history(&law);
    let state = law
        .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
        .unwrap();
    let h = 0.02;
    let stresses: [[f64; 9]; 4] = core::array::from_fn(|i| {
        let delta = [-2.0, -1.0, 1.0, 2.0][i] * h;
        let (changed, retained) = shifted(&law, &accepted, delta);
        changed
            .freeze(&retained, 0.0)
            .unwrap()
            .first_piola(&F)
            .unwrap()
    });
    for i in 0..9 {
        let derivative = (stresses[0][i] - 8.0 * stresses[1][i] + 8.0 * stresses[2][i]
            - stresses[3][i])
            / (12.0 * h);
        near(
            state.stress_temperature_derivative_pa_k()[i],
            derivative,
            5e-9,
            3e-7,
        );
        let mut plus = F;
        let mut minus = F;
        plus[i] += 1e-5;
        minus[i] -= 1e-5;
        let entropy_plus = law
            .frozen_thermoelastic_state(&accepted, &plus, T, jet(0.0))
            .unwrap()
            .entropy_j_m3_k();
        let entropy_minus = law
            .frozen_thermoelastic_state(&accepted, &minus, T, jet(0.0))
            .unwrap()
            .entropy_j_m3_k();
        near(
            state.stress_temperature_derivative_pa_k()[i],
            -(entropy_plus - entropy_minus) / 2e-5,
            3e-8,
            2e-6,
        );
    }
    let rotation = [0.8, -0.6, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 1.0];
    let rotated = law
        .frozen_thermoelastic_state(&accepted, &multiply(&rotation, &F), T, jet(0.0))
        .unwrap();
    near(
        rotated.free_energy_j_m3(),
        state.free_energy_j_m3(),
        2e-11,
        2e-11,
    );
    near(
        rotated.entropy_j_m3_k(),
        state.entropy_j_m3_k(),
        2e-12,
        2e-10,
    );
    near(
        rotated.mechanical_heat_capacity_j_m3_k(),
        state.mechanical_heat_capacity_j_m3_k(),
        3e-9,
        3e-8,
    );
    let rotated_beta = multiply(&rotation, state.stress_temperature_derivative_pa_k());
    for (actual, expected) in rotated
        .stress_temperature_derivative_pa_k()
        .iter()
        .zip(rotated_beta)
    {
        near(*actual, expected, 3e-10, 3e-8);
    }
}

#[test]
fn no_branches_recovers_elastic_owner_and_a_constant_jet_has_no_mechanical_entropy() {
    let mut law = law(false);
    law.branches.clear();
    let accepted = law.initial_state();
    let elastic = FiniteThermoelasticState::try_new(&law.equilibrium, &F, T, jet(0.0)).unwrap();
    let frozen = law
        .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
        .unwrap();
    // The numerical snapshot coincides, but only the elastic entry retains
    // its separate direct stable finite-increment owner.
    assert_eq!(frozen.free_energy_j_m3(), elastic.free_energy_j_m3());
    assert_eq!(frozen.entropy_j_m3_k(), elastic.entropy_j_m3_k());
    assert_eq!(
        frozen.internal_energy_j_m3(),
        elastic.internal_energy_j_m3()
    );
    assert_eq!(
        frozen.mechanical_heat_capacity_j_m3_k(),
        elastic.mechanical_heat_capacity_j_m3_k()
    );
    assert_eq!(frozen.first_piola_pa(), elastic.first_piola_pa());
    assert_eq!(
        frozen.isothermal_tangent_pa(),
        elastic.isothermal_tangent_pa()
    );
    assert_eq!(
        frozen.stress_temperature_derivative_pa_k(),
        elastic.stress_temperature_derivative_pa_k()
    );
    let with_memory = self::law(true);
    let history = history(&with_memory);
    let constant =
        ThermalReferenceJet::try_new(*with_memory.reference_deformation(), [0.0; 9], [0.0; 9])
            .unwrap();
    let state = with_memory
        .frozen_thermoelastic_state(&history, &F, T, constant)
        .unwrap();
    assert_eq!(state.entropy_j_m3_k(), 0.0);
    assert_eq!(state.mechanical_heat_capacity_j_m3_k(), 0.0);
    assert_eq!(state.stress_temperature_derivative_pa_k(), &[0.0; 9]);
    assert_eq!(state.internal_energy_j_m3(), state.free_energy_j_m3());
}

#[test]
fn binding_history_domain_and_public_parameter_errors_refuse_without_mutation() {
    let law = law(false);
    let accepted = history(&law);
    let before = accepted.clone();
    assert!(
        law.frozen_thermoelastic_state(&accepted, &F, T, jet(0.1))
            .is_err()
    );
    let wrong_law = law
        .clone()
        .with_reference_deformation(*jet(0.1).reference());
    assert!(
        wrong_law
            .frozen_thermoelastic_state(&accepted, &F, T, jet(0.1))
            .is_err()
    );
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.frozen_thermoelastic_state(&accepted, &F, temperature, jet(0.0))
                .is_err()
        );
    }
    let mut invalid = accepted.clone();
    invalid.inverse_viscous_cauchy_green[0] = [1.0, -1.0, -1.0, 0.0, 0.0, 0.0];
    assert!(
        law.frozen_thermoelastic_state(&invalid, &F, T, jet(0.0))
            .is_err()
    );
    invalid = accepted.clone();
    invalid.dissipated = -1.0;
    assert!(
        law.frozen_thermoelastic_state(&invalid, &F, T, jet(0.0))
            .is_err()
    );
    let mut invalid_law = law.clone();
    invalid_law.branches[0].tau_s = 0.0;
    assert!(
        invalid_law
            .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
            .is_err()
    );
    invalid_law = law.clone();
    invalid_law.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: -1.0 };
    assert!(
        invalid_law
            .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
            .is_err()
    );
    invalid_law = law.clone();
    invalid_law.equilibrium.stretch_limit = f64::NAN;
    assert!(
        invalid_law
            .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
            .is_err()
    );
    assert!(
        law.frozen_thermoelastic_state(&accepted, &[0.0; 9], T, jet(0.0))
            .is_err()
    );
    assert!(
        law.frozen_thermoelastic_state(&accepted, &F.map(|entry| entry * 10.0), T, jet(0.0))
            .is_err()
    );
    assert_eq!(accepted, before);
}

#[test]
fn elastic_finite_increments_survive_the_shared_kernel_and_visco_endpoints_refuse() {
    let law = law(false);
    let before = FiniteThermoelasticState::try_new(&law.equilibrium, &F, T, jet(0.0)).unwrap();
    let mut after_f = F;
    after_f[0] += 0.025;
    after_f[4] -= 0.015;
    let after =
        FiniteThermoelasticState::try_new(&law.equilibrium, &after_f, T + 0.3, jet(0.3)).unwrap();
    near(
        before.free_energy_increment_to(&after).unwrap(),
        after.free_energy_j_m3() - before.free_energy_j_m3(),
        2e-12,
        2e-12,
    );
    near(
        before.internal_energy_increment_to(&after).unwrap(),
        after.internal_energy_j_m3() - before.internal_energy_j_m3(),
        2e-12,
        2e-12,
    );
    assert_eq!(before.free_energy_increment_to(&before).unwrap(), 0.0);
    assert_eq!(before.internal_energy_increment_to(&before).unwrap(), 0.0);
    let accepted = history(&law);
    let frozen = law
        .frozen_thermoelastic_state(&accepted, &F, T, jet(0.0))
        .unwrap();
    assert!(before.free_energy_increment_to(&frozen).is_err());
    assert!(frozen.free_energy_increment_to(&before).is_err());
    assert!(frozen.free_energy_increment_to(&frozen).is_err());
    assert!(before.internal_energy_increment_to(&frozen).is_err());
    assert!(frozen.internal_energy_increment_to(&before).is_err());
    assert!(frozen.internal_energy_increment_to(&frozen).is_err());
}

#[test]
fn plane_stress_thermal_reduction_retains_memory_and_differentiates_transverse_balance() {
    let law = law(true);
    let accepted = history(&law);
    let before = accepted.clone();
    let in_plane = [1.18, 0.14, -0.03, 0.96];
    let settings = PlaneStressSettings::default();
    let state = law
        .frozen_thermoelastic_plane_stress_state(&accepted, &in_plane, T, jet(0.0), None, settings)
        .unwrap();
    let full = state.full_state();
    for i in [2, 5, 8] {
        assert!(full.first_piola_pa()[i].abs() < 2e-8);
    }
    assert!(state.mechanical_heat_capacity_j_m3_k() > full.mechanical_heat_capacity_j_m3_k());
    let f = full.deformation_gradient();
    let seed = Some([f[2], f[5], f[8]]);
    let h = 0.04;
    let samples = [-2.0, -1.0, 0.0, 1.0, 2.0].map(|offset| {
        let delta = offset * h;
        let (changed, retained) = shifted(&law, &accepted, delta);
        changed
            .frozen_thermoelastic_plane_stress_state(
                &retained,
                &in_plane,
                T + delta,
                jet(delta),
                seed,
                settings,
            )
            .unwrap()
    });
    let energies = samples.each_ref().map(|sample| sample.free_energy_j_m3());
    let psi_t = (energies[0] - 8.0 * energies[1] + 8.0 * energies[3] - energies[4]) / (12.0 * h);
    let psi_tt = (-energies[4] + 16.0 * energies[3] - 30.0 * energies[2] + 16.0 * energies[1]
        - energies[0])
        / (12.0 * h * h);
    near(state.entropy_j_m3_k(), -psi_t, 2e-9, 3e-7);
    near(
        state.mechanical_heat_capacity_j_m3_k(),
        -T * psi_tt,
        1e-5,
        3e-5,
    );
    for i in 0..4 {
        let p: [f64; 5] = samples.each_ref().map(|sample| sample.response().piola[i]);
        let derivative = (p[0] - 8.0 * p[1] + 8.0 * p[3] - p[4]) / (12.0 * h);
        near(
            state.stress_temperature_derivative_pa_k()[i],
            derivative,
            2e-8,
            2e-6,
        );
    }
    for (i, index) in [2, 5, 8].into_iter().enumerate() {
        let z: [f64; 5] = samples
            .each_ref()
            .map(|sample| sample.full_state().deformation_gradient()[index]);
        let derivative = (z[0] - 8.0 * z[1] + 8.0 * z[3] - z[4]) / (12.0 * h);
        near(
            state.transverse_temperature_derivative_per_k()[i],
            derivative,
            2e-9,
            2e-6,
        );
    }
    assert_eq!(accepted, before);
}

fn remainder_law(ogden: bool) -> ViscoHyperelastic {
    let mut material = law(ogden);
    material.equilibrium = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 4.0,
            c01: 0.0,
            kappa: 40.0,
        },
        4.0,
    )
    .unwrap();
    material
}

fn remainder_history(material: &ViscoHyperelastic) -> ViscoRemainderState {
    let first = material
        .step_incompressible_plane_stress_response(
            &material.initial_remainder_state(),
            &[1.25, 0.11, -0.02, 0.95],
            0.07,
            None,
            PlaneStressSettings::default(),
        )
        .unwrap()
        .material
        .state;
    material
        .step_incompressible_plane_stress_response(
            &first,
            &[1.13, -0.08, 0.04, 1.08],
            0.11,
            None,
            PlaneStressSettings::default(),
        )
        .unwrap()
        .material
        .state
}

struct ExactFrozenOracle {
    deformation: [f64; 9],
    elastic_deformation: [f64; 9],
    piola: [f64; 9],
    pressure: f64,
    energy: f64,
    entropy: f64,
    internal_energy: f64,
}

// With c01=0 and NH branches, the exact-volume stored energy is a quadratic
// form in F with M=8I+14Cv0^-1+5Cv1^-1. Its transverse minimum is analytic.
// No native thermal derivative, spectral lift, pressure solve or AD is used.
fn exact_frozen_oracle(
    accepted: &ViscoRemainderState,
    a: &[f64; 4],
    delta: f64,
) -> ExactFrozenOracle {
    let g: [f64; 9] = core::array::from_fn(|i| G[i] + delta * G1[i] + 0.5 * delta * delta * G2[i]);
    let gt: [f64; 9] = core::array::from_fn(|i| G1[i] + delta * G2[i]);
    let j = determinant(&g);
    let cg = cofactor(&g);
    let h = core::array::from_fn(|ij| cg[3 * (ij % 3) + ij / 3] / j);
    let ht = core::array::from_fn(|ij| h[3 * (ij % 3) + ij / 3]);
    let mut metric = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0];
    for (cv, mu) in accepted
        .inverse_viscous_cauchy_green
        .iter()
        .zip([14.0, 5.0])
    {
        let full = [
            cv[0], cv[3], cv[4], cv[3], cv[1], cv[5], cv[4], cv[5], cv[2],
        ];
        near(determinant(&full), 1.0, 3e-13, 0.0);
        for i in 0..9 {
            metric[i] += mu * full[i];
        }
    }
    let q = multiply(&multiply(&h, &metric), &ht);
    let direction = multiply(&h, &gt);
    let direction_t = core::array::from_fn(|ij| direction[3 * (ij % 3) + ij / 3]);
    let left = multiply(&direction, &q);
    let right = multiply(&q, &direction_t);
    let qt = core::array::from_fn(|i| -left[i] - right[i]);
    let jt: f64 = cg.iter().zip(gt).map(|(c, rate)| c * rate).sum();
    let d = a[0] * a[3] - a[1] * a[2];
    let mut f = [a[0], a[1], 0.0, a[2], a[3], 0.0, 0.0, 0.0, j / d];
    let mut ft = [0.0; 9];
    ft[8] = jt / d;
    for row in 0..2 {
        let numerator = f[3 * row] * q[2] + f[3 * row + 1] * q[5];
        let numerator_t = f[3 * row] * qt[2] + f[3 * row + 1] * qt[5];
        f[3 * row + 2] = -numerator / q[8];
        ft[3 * row + 2] = -(numerator_t + f[3 * row + 2] * qt[8]) / q[8];
    }
    let fq = multiply(&f, &q);
    let fqt = multiply(&f, &qt);
    let quadratic: f64 = f.iter().zip(fq).map(|(x, y)| x * y).sum();
    let quadratic_t: f64 = (0..9).map(|i| 2.0 * ft[i] * fq[i] + f[i] * fqt[i]).sum();
    let energy = 0.5 * j * (quadratic - 3.0 * 27.0);
    let entropy = -0.5 * jt * (quadratic - 3.0 * 27.0) - 0.5 * j * quadratic_t;
    let cof = cofactor(&f);
    let normal = f[8] * f[8] * q[8];
    let piola = core::array::from_fn(|i| j * fq[i] - normal * cof[i]);
    ExactFrozenOracle {
        deformation: f,
        elastic_deformation: multiply(&f, &h),
        piola,
        pressure: normal - quadratic / 3.0,
        energy,
        entropy,
        internal_energy: energy + (T + delta) * entropy,
    }
}

fn four_point(values: [f64; 4], h: f64) -> f64 {
    (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * h)
}

#[test]
fn exact_thermal_observer_uses_physical_aged_remainder_and_never_advances_memory() {
    for ogden in [false, true] {
        let material = remainder_law(ogden);
        let accepted = remainder_history(&material);
        let before = accepted.clone();
        assert!(accepted.dissipated > 0.0);
        let a = [1.18, 0.14, -0.03, 0.96];
        let settings = PlaneStressSettings::default();
        let state = material
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &accepted,
                &a,
                T,
                jet(0.0),
                None,
                settings,
            )
            .unwrap();
        let frozen = material.freeze_remainder(&accepted, 0.0).unwrap();
        let native = incompressible_plane_stress_with_kinematics(
            &frozen,
            &a,
            Some([
                accepted.deformation_gradient[2],
                accepted.deformation_gradient[5],
            ]),
            settings,
        )
        .unwrap();
        assert_eq!(state.response().full_piola, native.response.full_piola);
        assert_eq!(state.response().tangent, native.response.tangent);
        assert_eq!(state.response().pressure_pa, native.response.pressure_pa);
        near(
            state.free_energy_j_m3(),
            material
                .remainder_free_energy(&accepted, &state.response().deformation_gradient)
                .unwrap(),
            3e-11,
            3e-11,
        );
        let interval = material
            .freeze_remainder(&accepted, 0.2)
            .unwrap()
            .incompressible_remainder_stress_and_potential(&state.response().deformation_gradient)
            .unwrap()
            .1;
        assert!((interval - state.free_energy_j_m3()).abs() > 1e-6);
        let reset = material
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &material.initial_remainder_state(),
                &a,
                T,
                jet(0.0),
                None,
                settings,
            )
            .unwrap();
        assert!((reset.free_energy_j_m3() - state.free_energy_j_m3()).abs() > 1e-4);
        let mut shifted_times = material.clone();
        shifted_times.branches[0].tau_s *= 1e-8;
        shifted_times.branches[1].tau_s *= 1e8;
        let slow = shifted_times
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &accepted,
                &a,
                T,
                jet(0.0),
                None,
                settings,
            )
            .unwrap();
        assert_eq!(slow.free_energy_j_m3(), state.free_energy_j_m3());
        assert_eq!(slow.entropy_j_m3_k(), state.entropy_j_m3_k());
        assert_eq!(
            slow.mechanical_heat_capacity_j_m3_k(),
            state.mechanical_heat_capacity_j_m3_k()
        );
        assert_eq!(
            slow.rate([0.08, -0.03, 0.01, -0.05], 2.3).unwrap(),
            state.rate([0.08, -0.03, 0.01, -0.05], 2.3).unwrap()
        );
        assert_eq!(accepted, before);
    }
}

#[test]
fn exact_thermal_memory_matches_independent_relaxed_potential_and_motion() {
    let material = remainder_law(false);
    let accepted = remainder_history(&material);
    let a = [1.18, 0.14, -0.03, 0.96];
    let state = material
        .frozen_incompressible_thermoelastic_plane_stress_state(
            &accepted,
            &a,
            T,
            jet(0.0),
            None,
            PlaneStressSettings::default(),
        )
        .unwrap();
    let oracle = exact_frozen_oracle(&accepted, &a, 0.0);
    for (actual, expected) in [
        (state.free_energy_j_m3(), oracle.energy),
        (state.entropy_j_m3_k(), oracle.entropy),
        (state.internal_energy_j_m3(), oracle.internal_energy),
        (state.response().pressure_pa, oracle.pressure),
    ] {
        near(actual, expected, 3e-10, 3e-10);
    }
    for i in 0..9 {
        near(
            state.response().deformation_gradient[i],
            oracle.deformation[i],
            3e-11,
            3e-11,
        );
        near(
            state.response().full_piola[i],
            oracle.piola[i],
            3e-10,
            3e-10,
        );
    }
    let h = 0.02;
    let thermal = [-2.0 * h, -h, h, 2.0 * h].map(|delta| exact_frozen_oracle(&accepted, &a, delta));
    near(
        state.mechanical_heat_capacity_j_m3_k(),
        T * four_point(thermal.each_ref().map(|s| s.entropy), h),
        3e-8,
        3e-7,
    );
    let adot = [0.08, -0.03, 0.01, -0.05];
    let tdot = 2.3;
    let rate = state.rate(adot, tdot).unwrap();
    let samples = [-2.0 * h, -h, h, 2.0 * h].map(|time| {
        exact_frozen_oracle(
            &accepted,
            &core::array::from_fn(|i| a[i] + time * adot[i]),
            tdot * time,
        )
    });
    for i in 0..9 {
        near(
            rate.material().deformation_gradient_rate_per_s()[i],
            four_point(samples.each_ref().map(|s| s.deformation[i]), h),
            4e-8,
            3e-7,
        );
        near(
            rate.material().elastic_deformation_gradient_rate_per_s()[i],
            four_point(samples.each_ref().map(|s| s.elastic_deformation[i]), h),
            4e-8,
            3e-7,
        );
        near(
            rate.material().first_piola_rate_pa_s()[i],
            four_point(samples.each_ref().map(|s| s.piola[i]), h),
            3e-7,
            3e-7,
        );
    }
    for (actual, values) in [
        (
            rate.pressure_rate_pa_s(),
            samples.each_ref().map(|s| s.pressure),
        ),
        (
            rate.material().free_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.energy),
        ),
        (
            rate.material().entropy_rate_w_m3_k(),
            samples.each_ref().map(|s| s.entropy),
        ),
        (
            rate.material().internal_energy_rate_w_m3(),
            samples.each_ref().map(|s| s.internal_energy),
        ),
    ] {
        near(actual, four_point(values, h), 3e-7, 3e-7);
    }
    for i in [2, 5, 8] {
        near(rate.material().first_piola_rate_pa_s()[i], 0.0, 2e-10, 0.0);
    }
}

#[test]
fn exact_thermal_entropy_matches_retained_reference_work_without_relaxation_heat() {
    let material = remainder_law(true);
    let aged = remainder_history(&material);
    let a = [1.18, 0.14, -0.03, 0.96];
    let settings = PlaneStressSettings::default();
    let accepted = material
        .step_incompressible_plane_stress_response(&aged, &a, 0.0, None, settings)
        .unwrap()
        .material
        .state;
    let before = accepted.clone();
    let center = material
        .frozen_incompressible_thermoelastic_plane_stress_state(
            &accepted,
            &a,
            T,
            jet(0.0),
            None,
            settings,
        )
        .unwrap();
    let h = 0.02;
    let work = [-2.0 * h, -h, h, 2.0 * h].map(|delta| {
        let target = material
            .clone()
            .with_reference_deformation(*jet(delta).reference());
        let event = target
            .change_incompressible_plane_stress_reference_from(&material, &accepted, None, settings)
            .unwrap();
        assert_eq!(
            event.material.state.inverse_viscous_cauchy_green,
            accepted.inverse_viscous_cauchy_green
        );
        assert_eq!(event.material.state.dissipated, accepted.dissipated);
        assert_eq!(event.material.dissipated_pa, 0.0);
        let endpoint = target
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &event.material.state,
                &a,
                T + delta,
                jet(delta),
                None,
                settings,
            )
            .unwrap();
        near(
            endpoint.free_energy_j_m3(),
            event.material.state.stored,
            3e-10,
            3e-10,
        );
        [
            event.reference_work_pa + event.constraint_reference_work_pa,
            event.constraint_reference_work_pa,
            event.endpoint_work_defect_pa,
        ]
    });
    near(
        four_point(work.map(|w| w[0]), h),
        -center.entropy_j_m3_k(),
        3e-8,
        3e-6,
    );
    near(
        four_point(work.map(|w| w[1]), h),
        center.response().pressure_pa * center.natural_volume_derivative_per_k(),
        3e-8,
        3e-6,
    );
    near(four_point(work.map(|w| w[2]), h), 0.0, 3e-8, 0.0);
    assert_eq!(accepted, before);
}

#[test]
fn exact_thermal_visco_native_limit_binding_refusals_and_increment_ownership() {
    let mut material = remainder_law(false);
    material.branches.clear();
    let accepted = material.initial_remainder_state();
    let a = [1.18, 0.14, -0.03, 0.96];
    let settings = PlaneStressSettings::default();
    let frozen = material
        .frozen_incompressible_thermoelastic_plane_stress_state(
            &accepted,
            &a,
            T,
            jet(0.0),
            None,
            settings,
        )
        .unwrap();
    let native = IncompressibleThermoelasticPlaneStressState::try_new(
        &material.equilibrium,
        &a,
        T,
        jet(0.0),
        None,
        settings,
    )
    .unwrap();
    for (actual, expected) in [
        (frozen.free_energy_j_m3(), native.free_energy_j_m3()),
        (frozen.entropy_j_m3_k(), native.entropy_j_m3_k()),
        (frozen.internal_energy_j_m3(), native.internal_energy_j_m3()),
        (
            frozen.mechanical_heat_capacity_j_m3_k(),
            native.mechanical_heat_capacity_j_m3_k(),
        ),
        (frozen.response().pressure_pa, native.response().pressure_pa),
    ] {
        near(actual, expected, 3e-11, 3e-11);
    }
    assert_eq!(native.free_energy_increment_to(&native).unwrap(), 0.0);
    assert!(frozen.free_energy_increment_to(&frozen).is_err());
    assert!(frozen.internal_energy_increment_to(&frozen).is_err());
    assert!(native.free_energy_increment_to(&frozen).is_err());
    assert!(frozen.internal_energy_increment_to(&native).is_err());
    let material = remainder_law(false);
    let accepted = remainder_history(&material);
    let saved = accepted.clone();
    assert!(
        material
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &accepted,
                &a,
                T,
                jet(0.1),
                None,
                settings,
            )
            .is_err()
    );
    let foreign = material
        .clone()
        .with_reference_deformation(*jet(0.1).reference());
    assert!(
        foreign
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &accepted,
                &a,
                T,
                jet(0.1),
                None,
                settings,
            )
            .is_err()
    );
    let mut invalid = accepted.clone();
    invalid.inverse_viscous_cauchy_green[0] = [1.0, -1.0, -1.0, 0.0, 0.0, 0.0];
    assert!(
        material
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &invalid,
                &a,
                T,
                jet(0.0),
                None,
                settings,
            )
            .is_err()
    );
    invalid = accepted.clone();
    invalid.dissipated = -1.0;
    assert!(
        material
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &invalid,
                &a,
                T,
                jet(0.0),
                None,
                settings,
            )
            .is_err()
    );
    let unsupported = law(false); // Original visco mixed schema refuses coupled NH.
    assert!(
        unsupported
            .frozen_incompressible_thermoelastic_plane_stress_state(
                &unsupported.initial_remainder_state(),
                &a,
                T,
                jet(0.0),
                None,
                settings,
            )
            .is_err()
    );
    assert_eq!(accepted, saved);
}
