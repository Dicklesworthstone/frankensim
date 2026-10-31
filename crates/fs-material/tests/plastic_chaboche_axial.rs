//! Independent axial equations, reversal and work for nonlinear backstress.
use fs_material::elastic::IsotropicElastic;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::UniaxialVoce;
use fs_material::plastic::voce::VoceJ2;

const E: f64 = 1200.0;
const Y: f64 = 3.0;
const HI: f64 = 40.0;
const HK: f64 = 30.0;
const Q: f64 = 2.0;
const B: f64 = 80.0;
const BRANCHES: [(f64, f64); 2] = [(160.0, 120.0), (60.0, 25.0)];

fn voce(kinematic: f64) -> VoceJ2 {
    VoceJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(E, 0.2, 0.05).unwrap(),
            Y,
            HI,
            kinematic,
        )
        .unwrap(),
        Q,
        B,
    )
    .unwrap()
}
fn law() -> UniaxialChaboche<2> {
    UniaxialChaboche::new(
        ChabocheJ2::new(
            voce(HK),
            BRANCHES.map(|(c, g)| BackstressBranch::new(c, g).unwrap()),
        )
        .unwrap(),
    )
    .unwrap()
}
fn radius(alpha: f64) -> f64 {
    Y + HI * alpha + Q * (1.0 - (-B * alpha).exp())
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.15e} != {b:.15e}, allowance {tolerance}"
    );
}

#[test]
fn axial_return_matches_independent_prescribed_plastic_increment_and_full_tensor() {
    let law = law();
    let virgin = law.initial_state();
    for p in [0.001_f64, 0.006, 0.015] {
        let translations = BRANCHES.map(|(c, g)| c * p / (1.0 + g * p));
        let stress = radius(p) + HK * p + translations.iter().sum::<f64>();
        let strain = p + stress / E;
        let actual = law.response(strain, &virgin).unwrap();
        close(actual.stress_pa, stress, 2e-11);
        close(actual.state.plastic_strain(), p, 2e-14);
        close(actual.state.equivalent_plastic_strain(), p, 2e-14);
        for (a, b) in actual.state.backstress_pa().iter().zip(translations) {
            close(*a, b, 3e-12);
        }
        let tensor = law
            .backbone()
            .response(&actual.strain, virgin.tensor_history())
            .unwrap();
        close(tensor.stress_pa[0], stress, 3e-11);
        for transverse in &tensor.stress_pa[1..] {
            close(*transverse, 0.0, 3e-11);
        }
        let h = 1e-7;
        let fd = (law.stress(strain + h, &virgin) - law.stress(strain - h, &virgin)) / (2.0 * h);
        close(actual.tangent_pa, fd, 2e-6);
        assert_eq!(virgin, law.initial_state());
    }
}

#[test]
fn reversal_retains_independent_recovering_backstress_and_exact_axial_tangent() {
    let law = law();
    let mut old = law.response(0.015, &law.initial_state()).unwrap().state;
    for (sign, increment) in [(-1.0, 0.004), (1.0, 0.008), (-1.0, 0.013)] {
        let p = old.plastic_strain() + sign * increment;
        let alpha = old.equivalent_plastic_strain() + increment;
        let previous = old.backstress_pa();
        let translation: [f64; 2] = core::array::from_fn(|i| {
            let (c, g) = BRANCHES[i];
            (previous[i] + c * sign * increment) / (1.0 + g * increment)
        });
        let stress = HK * p + translation.iter().sum::<f64>() + sign * radius(alpha);
        let strain = p + stress / E;
        let actual = law.response(strain, &old).unwrap();
        close(actual.stress_pa, stress, 2e-11);
        close(actual.state.plastic_strain(), p, 2e-14);
        close(actual.state.equivalent_plastic_strain(), alpha, 2e-14);
        for (a, b) in actual.state.backstress_pa().iter().zip(translation) {
            close(*a, b, 3e-12);
        }
        let slope = HI
            + HK
            + Q * B * (-B * alpha).exp()
            + BRANCHES
                .iter()
                .enumerate()
                .map(|(i, &(c, g))| (c - g * sign * previous[i]) / (1.0 + g * increment).powi(2))
                .sum::<f64>();
        close(actual.tangent_pa, E * slope / (E + slope), 2e-9);
        let h = 1e-7;
        close(
            actual.tangent_pa,
            (law.stress(strain + h, &old) - law.stress(strain - h, &old)) / (2.0 * h),
            2e-6,
        );
        old = actual.state;
    }
}

#[test]
fn zero_recovery_matches_additional_linear_prager_hardening_on_reversal() {
    let law = UniaxialChaboche::new(
        ChabocheJ2::new(voce(HK), [BackstressBranch::new(160.0, 0.0).unwrap()]).unwrap(),
    )
    .unwrap();
    let reference = UniaxialVoce::new(voce(HK + 160.0)).unwrap();
    let (mut state, mut old_reference) = (law.initial_state(), reference.initial_state());
    for strain in [0.001, 0.012, 0.007, -0.011, 0.013, -0.008, 0.0] {
        let actual = law.response(strain, &state).unwrap();
        let expected = reference.response(strain, &old_reference);
        close(actual.stress_pa, expected.stress_pa, 3e-11);
        close(actual.tangent_pa, expected.tangent_pa, 2e-8);
        close(
            actual.state.plastic_strain(),
            expected.state.plastic_strain(),
            2e-14,
        );
        close(
            actual.state.backstress_pa()[0],
            160.0 * actual.state.plastic_strain(),
            2e-12,
        );
        state = actual.state;
        old_reference = expected.state;
    }
}

#[test]
fn cyclic_axial_work_distinguishes_recovery_heat_storage_and_numerical_loss() {
    let law = law();
    let mut state = law.initial_state();
    let mut old_strain = 0.0;
    let mut recovery = 0.0;
    for target in [0.014, -0.014, 0.014, -0.014, 0.0] {
        let start = old_strain;
        for step in 1..=12 {
            let strain = start + (target - start) * step as f64 / 12.0;
            let actual = law.response(strain, &state).unwrap();
            let audit = law
                .advance_with_energy(
                    old_strain,
                    &state,
                    strain,
                    EnergyTolerances {
                        yield_pa: 2e-10,
                        energy_j_m3: 2e-10,
                    },
                )
                .unwrap();
            let delta =
                actual.state.equivalent_plastic_strain() - state.equivalent_plastic_strain();
            close(
                audit.energy.physical_dissipation_j_m3,
                Y * delta + audit.recovery_dissipation_j_m3,
                2e-12,
            );
            close(
                audit.energy.mechanical_work_j_m3,
                actual.stress_pa * (strain - old_strain),
                2e-11,
            );
            assert!(audit.energy.numerical_dissipation_j_m3 >= 0.0);
            assert!(audit.energy.residual_j_m3.abs() < 2e-10);
            recovery += audit.recovery_dissipation_j_m3;
            state = actual.state;
            old_strain = strain;
        }
    }
    assert!(recovery > 0.01);
    let snapshot = state.clone();
    assert!(law.response(0.1, &state).is_err());
    assert!(law.response(f64::NAN, &state).is_err());
    assert_eq!(state, snapshot);
}
