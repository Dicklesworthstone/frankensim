//! Independent invariant, return-map and work checks for the authored
//! pressure-sensitive small-strain law. These tests do not fit foam data or
//! qualify finite crushing, evolving porosity or a physical heat partition.

use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::porous::{PorousPlasticity, PorousState};
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt};

const YOUNG: f64 = 1.0e9;
const POISSON: f64 = 0.22;
const YIELD: f64 = 1.0e6;
const HARDENING: f64 = 1.0e8;
const SENSITIVITY: f64 = 1.5;
const STRAIN_LIMIT: f64 = 0.05;
const WEIGHT: Voigt = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];

fn material(hardening: f64, sensitivity: f64) -> PorousPlasticity {
    PorousPlasticity::new(
        IsotropicElastic::new(YOUNG, POISSON, STRAIN_LIMIT).unwrap(),
        YIELD,
        hardening,
        sensitivity,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}",
    );
}

fn contract(a: &Voigt, b: &Voigt) -> f64 {
    (0..6).map(|i| WEIGHT[i] * a[i] * b[i]).sum()
}

fn elastic_energy(strain: &Voigt) -> f64 {
    let bulk = YOUNG / (3.0 * (1.0 - 2.0 * POISSON));
    let shear = YOUNG / (2.0 * (1.0 + POISSON));
    let trace: f64 = strain[..3].iter().sum();
    let dev: Voigt = std::array::from_fn(|i| strain[i] - if i < 3 { trace / 3.0 } else { 0.0 });
    shear * contract(&dev, &dev) + 0.5 * bulk * trace * trace
}

// Proper orthogonal rotation from the unit quaternion (1, 2, 2, 0)/3.
// Every physical tensor component participates; this is not an axis swap.
fn rotate(tensor: &Voigt) -> Voigt {
    let r = [
        [1.0 / 9.0, 8.0 / 9.0, 4.0 / 9.0],
        [8.0 / 9.0, 1.0 / 9.0, -4.0 / 9.0],
        [-4.0 / 9.0, 4.0 / 9.0, -7.0 / 9.0],
    ];
    let matrix = [
        [tensor[0], tensor[3], tensor[5]],
        [tensor[3], tensor[1], tensor[4]],
        [tensor[5], tensor[4], tensor[2]],
    ];
    let mut rotated = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for a in 0..3 {
                for b in 0..3 {
                    rotated[i][j] += r[i][a] * matrix[a][b] * r[j][b];
                }
            }
        }
    }
    [
        rotated[0][0],
        rotated[1][1],
        rotated[2][2],
        rotated[0][1],
        rotated[1][2],
        rotated[2][0],
    ]
}

#[test]
fn hydrostatic_and_shear_returns_match_independent_closed_forms() {
    let bulk = YOUNG / (3.0 * (1.0 - 2.0 * POISSON));
    let shear = YOUNG / (2.0 * (1.0 + POISSON));
    let c = 1.0 + SENSITIVITY * SENSITIVITY / 9.0;
    let root_c = c.sqrt();
    for hardening in [0.0, HARDENING] {
        let law = material(hardening, SENSITIVITY);
        let state = law.initial_state();

        let shortening = 0.01;
        let hydro = law
            .response(
                &[-shortening, -shortening, -shortening, 0.0, 0.0, 0.0],
                &state,
            )
            .unwrap();
        let hydro_modulus = bulk * SENSITIVITY * SENSITIVITY / c;
        let increment =
            (3.0 * bulk * shortening * SENSITIVITY / root_c - YIELD) / (hardening + hydro_modulus);
        let pressure = root_c / SENSITIVITY * (YIELD + hardening * increment);
        close(hydro.plastic_multiplier, increment, 1.0e-12);
        close(hydro.state.alpha, increment, 1.0e-12);
        close(hydro.yield_residual_pa, 0.0, 1.0e-4);
        for i in 0..3 {
            close(hydro.stress_pa[i], -pressure, 0.1);
            close(
                hydro.state.plastic_strain[i],
                -SENSITIVITY * increment / (3.0 * root_c),
                1.0e-12,
            );
            let isotropic_tangent: f64 = hydro.tangent_pa[i][..3].iter().sum();
            close(
                isotropic_tangent,
                3.0 * bulk * hardening / (hardening + hydro_modulus),
                2.0,
            );
        }
        for i in 3..6 {
            close(hydro.stress_pa[i], 0.0, 1.0e-8);
            close(hydro.state.plastic_strain[i], 0.0, 1.0e-14);
        }

        // Tensor shear epsilon_xy=t corresponds to engineering shear 2*t.
        let tensor_shear = 0.012;
        let shear_response = law
            .response(&[0.0, 0.0, 0.0, tensor_shear, 0.0, 0.0], &state)
            .unwrap();
        let root_three = 3.0_f64.sqrt();
        let increment = (2.0 * root_three * shear * tensor_shear / root_c - YIELD)
            / (hardening + 3.0 * shear / c);
        let expected_shear = root_c / root_three * (YIELD + hardening * increment);
        close(shear_response.plastic_multiplier, increment, 1.0e-12);
        close(shear_response.state.alpha, increment, 1.0e-12);
        close(shear_response.stress_pa[3], expected_shear, 0.1);
        close(
            shear_response.state.plastic_strain[3],
            root_three * increment / (2.0 * root_c),
            1.0e-12,
        );
        close(
            shear_response.tangent_pa[3][3],
            2.0 * shear * hardening / (hardening + 3.0 * shear / c),
            2.0,
        );
        for i in [0, 1, 2, 4, 5] {
            close(shear_response.stress_pa[i], 0.0, 1.0e-8);
            close(shear_response.state.plastic_strain[i], 0.0, 1.0e-14);
        }
    }
}

#[test]
fn mixed_pressure_and_shear_have_the_consistent_weighted_symmetric_tangent() {
    let law = material(HARDENING, SENSITIVITY);
    let old = law
        .response(
            &[-0.006, -0.004, -0.003, 0.001, 0.0, 0.0],
            &law.initial_state(),
        )
        .unwrap()
        .state;
    let strain = [-0.012, -0.005, 0.002, 0.004, -0.003, 0.002];
    let actual = law.response(&strain, &old).unwrap();
    assert!(actual.plastic_multiplier > 0.0);
    let scale = actual
        .tangent_pa
        .iter()
        .flatten()
        .fold(0.0_f64, |s, v| s.max(v.abs()));
    let h = 1.0e-7;
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += h;
        minus[j] -= h;
        let upper = law.response(&plus, &old).unwrap();
        let lower = law.response(&minus, &old).unwrap();
        for i in 0..6 {
            close(
                actual.tangent_pa[i][j],
                (upper.stress_pa[i] - lower.stress_pa[i]) / (2.0 * h),
                2.0e-7 * scale,
            );
            close(
                WEIGHT[i] * actual.tangent_pa[i][j],
                WEIGHT[j] * actual.tangent_pa[j][i],
                1.0e-12 * scale,
            );
        }
    }
}

#[test]
fn zero_pressure_sensitivity_matches_native_j2_history_through_reversal() {
    let porous = material(HARDENING, 0.0);
    let j2 = J2Plasticity::new(
        IsotropicElastic::new(YOUNG, POISSON, STRAIN_LIMIT).unwrap(),
        YIELD,
        HARDENING,
    )
    .unwrap();
    let mut porous_state = porous.initial_state();
    let mut j2_state = j2.initial_state();
    let history: [Voigt; 4] = [
        [0.005, -0.001, -0.002, 0.003, 0.0, 0.0],
        [0.001, 0.0002, -0.0001, 0.0002, 0.0, 0.0],
        [-0.006, 0.002, 0.001, -0.004, 0.001, 0.002],
        [0.004, 0.004, 0.004, 0.0, 0.0, 0.0],
    ];
    for strain in history {
        let actual = porous.response(&strain, &porous_state).unwrap();
        let expected_stress = j2.stress(&strain, &j2_state);
        let expected_tangent = j2.tangent(&strain, &j2_state);
        let next = j2.update_state(&strain, &j2_state);
        for i in 0..6 {
            close(actual.stress_pa[i], expected_stress[i], 0.1);
            close(
                actual.state.plastic_strain[i],
                next.plastic_strain[i],
                1.0e-12,
            );
            for j in 0..6 {
                close(actual.tangent_pa[i][j], expected_tangent[i][j], 2.0);
            }
        }
        close(actual.state.alpha, next.alpha, 1.0e-12);
        close(actual.state.plastic_strain[..3].iter().sum(), 0.0, 1.0e-13);
        porous_state = actual.state;
        j2_state = next;
    }
}

#[test]
fn response_rotates_with_the_history_and_preserves_stress_unit_rescaling() {
    let law = material(HARDENING, SENSITIVITY);
    let old = law
        .response(
            &[-0.006, -0.004, -0.003, 0.001, 0.0, 0.0],
            &law.initial_state(),
        )
        .unwrap()
        .state;
    let strain = [-0.012, -0.005, 0.002, 0.004, -0.003, 0.002];
    let original = law.response(&strain, &old).unwrap();
    let rotated_history = PorousState {
        plastic_strain: rotate(&old.plastic_strain),
        alpha: old.alpha,
    };
    let rotated = law.response(&rotate(&strain), &rotated_history).unwrap();
    let expected_stress = rotate(&original.stress_pa);
    let expected_plastic = rotate(&original.state.plastic_strain);
    for i in 0..6 {
        close(rotated.stress_pa[i], expected_stress[i], 0.1);
        close(
            rotated.state.plastic_strain[i],
            expected_plastic[i],
            1.0e-12,
        );
    }
    close(rotated.state.alpha, original.state.alpha, 1.0e-12);

    let direction = [0.7, -0.3, 0.2, 0.4, -0.2, 0.1];
    let original_directional_stress: Voigt = std::array::from_fn(|i| {
        (0..6)
            .map(|j| original.tangent_pa[i][j] * direction[j])
            .sum()
    });
    let rotated_direction = rotate(&direction);
    let expected_directional_stress = rotate(&original_directional_stress);
    for i in 0..6 {
        let actual: f64 = (0..6)
            .map(|j| rotated.tangent_pa[i][j] * rotated_direction[j])
            .sum();
        close(actual, expected_directional_stress[i], 2.0);
    }

    for scale in [1.0e-6, 1.0e6] {
        let scaled = PorousPlasticity::new(
            IsotropicElastic::new(YOUNG * scale, POISSON, STRAIN_LIMIT).unwrap(),
            YIELD * scale,
            HARDENING * scale,
            SENSITIVITY,
        )
        .unwrap()
        .response(&strain, &old)
        .unwrap();
        close(scaled.state.alpha, original.state.alpha, 1.0e-12);
        close(
            scaled.plastic_multiplier,
            original.plastic_multiplier,
            1.0e-12,
        );
        for i in 0..6 {
            close(scaled.stress_pa[i] / scale, original.stress_pa[i], 0.1);
            close(
                scaled.state.plastic_strain[i],
                original.state.plastic_strain[i],
                1.0e-12,
            );
            for j in 0..6 {
                close(
                    scaled.tangent_pa[i][j] / scale,
                    original.tangent_pa[i][j],
                    2.0,
                );
            }
        }
    }
}

#[test]
fn invalid_parameters_histories_and_tensor_strain_domains_refuse_without_mutation() {
    for (yield_stress, hardening, sensitivity) in [
        (0.0, HARDENING, SENSITIVITY),
        (f64::NAN, HARDENING, SENSITIVITY),
        (YIELD, -1.0, SENSITIVITY),
        (YIELD, f64::INFINITY, SENSITIVITY),
        (YIELD, HARDENING, -0.1),
        (YIELD, HARDENING, f64::NAN),
    ] {
        assert!(
            PorousPlasticity::new(
                IsotropicElastic::new(YOUNG, POISSON, STRAIN_LIMIT).unwrap(),
                yield_stress,
                hardening,
                sensitivity,
            )
            .is_err()
        );
    }
    let invalid_elastic = IsotropicElastic {
        youngs: YOUNG,
        poisson: 0.5,
        strain_limit: STRAIN_LIMIT,
    };
    assert!(PorousPlasticity::new(invalid_elastic, YIELD, HARDENING, SENSITIVITY).is_err());

    let law = material(HARDENING, SENSITIVITY);
    let invalid_histories = [
        PorousState {
            plastic_strain: [0.0; 6],
            alpha: -0.001,
        },
        PorousState {
            plastic_strain: [0.0; 6],
            alpha: f64::NAN,
        },
        PorousState {
            plastic_strain: [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0],
            alpha: 0.01,
        },
        PorousState {
            plastic_strain: [0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
            alpha: 0.0,
        },
        PorousState {
            plastic_strain: [-0.003, -0.003, -0.003, 0.0, 0.0, 0.0],
            alpha: 0.001,
        },
    ];
    for state in invalid_histories {
        let bits = (
            state.alpha.to_bits(),
            state.plastic_strain.map(f64::to_bits),
        );
        assert!(law.response(&[0.0; 6], &state).is_err());
        assert_eq!(
            bits,
            (
                state.alpha.to_bits(),
                state.plastic_strain.map(f64::to_bits)
            )
        );
    }

    let accepted = law
        .response(
            &[-0.006, -0.004, -0.003, 0.001, 0.0, 0.0],
            &law.initial_state(),
        )
        .unwrap()
        .state;
    let saved = accepted.clone();
    for strain in [
        [0.03, 0.03, 0.03, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.04, 0.0, 0.0],
        [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0],
    ] {
        assert!(law.response(&strain, &accepted).is_err());
        assert_eq!(accepted, saved);
    }
    let mut changed = material(HARDENING, SENSITIVITY);
    changed.pressure_sensitivity = -1.0;
    assert!(changed.response(&[0.0; 6], &accepted).is_err());
    assert_eq!(accepted, saved);

    // Finite, norm-consistent internal variables alone do not prove that an
    // old total strain was an accepted state of this fixed constitutive law.
    let uncommitted = [-0.01, -0.01, -0.01, 0.0, 0.0, 0.0];
    assert!(
        law.advance_with_energy(
            &uncommitted,
            &law.initial_state(),
            &uncommitted,
            EnergyTolerances {
                energy_j_m3: 1.0e-6,
                yield_pa: 1.0e-4
            },
        )
        .is_err()
    );
}

#[test]
fn work_account_separates_heat_storage_and_numerical_loss_during_load_and_unload() {
    let law = material(HARDENING, SENSITIVITY);
    let strain = [-0.012, -0.006, -0.003, 0.004, -0.002, 0.001];
    let tolerances = EnergyTolerances {
        energy_j_m3: 1.0e-6,
        yield_pa: 1.0e-4,
    };
    let increment = law
        .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, tolerances)
        .unwrap();
    let response = &increment.response;
    let elastic_strain: Voigt =
        std::array::from_fn(|i| strain[i] - response.state.plastic_strain[i]);
    let expected_elastic = elastic_energy(&elastic_strain);
    let expected_hardening = 0.5 * HARDENING * response.state.alpha * response.state.alpha;
    let expected_heat = YIELD * response.state.alpha;
    let expected_numerical = expected_elastic + expected_hardening;
    let expected_work = contract(&response.stress_pa, &strain);
    assert!(expected_heat > 0.0 && expected_numerical > 0.0);
    assert_eq!(increment.strain, strain);
    close(increment.energy.initial.total_j_m3(), 0.0, 1.0e-12);
    close(
        increment.energy.final_state.elastic_j_m3,
        expected_elastic,
        1.0e-6,
    );
    close(
        increment.energy.final_state.hardening_j_m3,
        expected_hardening,
        1.0e-6,
    );
    close(
        increment.energy.final_state.kinematic_hardening_j_m3,
        0.0,
        1.0e-12,
    );
    close(increment.energy.mechanical_work_j_m3, expected_work, 1.0e-6);
    close(
        increment.energy.physical_dissipation_j_m3,
        expected_heat,
        1.0e-6,
    );
    close(
        increment.energy.numerical_dissipation_j_m3,
        expected_numerical,
        1.0e-6,
    );
    close(
        expected_work,
        expected_elastic + expected_hardening + expected_heat + expected_numerical,
        1.0e-6,
    );
    close(increment.energy.residual_j_m3, 0.0, 1.0e-6);
    close(
        contract(&response.stress_pa, &response.state.plastic_strain),
        (YIELD + HARDENING * response.state.alpha) * response.state.alpha,
        1.0e-6,
    );
    let stored = law.stored_energy(&strain, &response.state).unwrap();
    close(stored.elastic_j_m3, expected_elastic, 1.0e-6);
    close(stored.hardening_j_m3, expected_hardening, 1.0e-6);

    // At zero elastic strain the stress vanishes. The endpoint-work convention
    // leaves the removed elastic storage as numerical loss, not physical heat.
    let released_strain = response.state.plastic_strain;
    let unloaded = law
        .advance_with_energy(&strain, &response.state, &released_strain, tolerances)
        .unwrap();
    assert_eq!(unloaded.response.state, response.state);
    for stress in unloaded.response.stress_pa {
        close(stress, 0.0, 1.0e-8);
    }
    close(unloaded.energy.final_state.elastic_j_m3, 0.0, 1.0e-12);
    close(
        unloaded.energy.final_state.hardening_j_m3,
        expected_hardening,
        1.0e-6,
    );
    close(unloaded.energy.physical_dissipation_j_m3, 0.0, 1.0e-12);
    close(unloaded.energy.mechanical_work_j_m3, 0.0, 1.0e-8);
    close(
        unloaded.energy.numerical_dissipation_j_m3,
        expected_elastic,
        1.0e-6,
    );
    close(unloaded.energy.residual_j_m3, 0.0, 1.0e-6);
}
