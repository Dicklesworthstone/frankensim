//! The original finite-strain branch potential is usable by spatial Newton
//! solves without confusing its value with physical storage or heat.

use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

fn law() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 12_000.0,
                lambda: 24_000.0,
            },
            4.0,
        )
        .unwrap(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 30_000.0 },
            tau_s: 0.2,
        }],
    )
    .unwrap()
}

#[test]
fn frozen_incremental_potential_generates_native_stress_and_tangent() {
    let law = law();
    let previous_f = [1.2, 0.1, 0.0, 0.02, 0.94, 0.0, 0.0, 0.0, 0.91];
    let (accepted, _) = law.step(&law.initial_state(), &previous_f, 0.04).unwrap();
    let snapshot = accepted.clone();
    let f = [1.16, 0.07, 0.02, -0.01, 0.97, 0.03, 0.01, 0.0, 0.93];
    let direction = [0.3, -0.2, 0.1, 0.07, -0.13, 0.17, -0.08, 0.11, 0.09];
    let h = 2e-6;
    let plus = core::array::from_fn(|i| f[i] + h * direction[i]);
    let minus = core::array::from_fn(|i| f[i] - h * direction[i]);
    for dt in [0.0, 0.025, 0.3] {
        let frozen = law.freeze(&accepted, dt).unwrap();
        let (p, potential) = frozen.stress_and_potential(&f).unwrap();
        let tangent = frozen.consistent_tangent(&f).unwrap();
        let response = law.step_response(&accepted, &f, dt).unwrap();
        assert_eq!(p, response.first_piola_pa);
        assert_eq!(potential, response.incremental_potential_pa);
        let fd_work = (frozen.energy_density(&plus) - frozen.energy_density(&minus)) / (2.0 * h);
        let expected_work: f64 = p.iter().zip(direction).map(|(a, b)| a * b).sum();
        assert!(
            (fd_work - expected_work).abs() < 2e-5,
            "dt={dt}: {fd_work} versus {expected_work}"
        );
        let p_plus = frozen.first_piola(&plus).unwrap();
        let p_minus = frozen.first_piola(&minus).unwrap();
        for row in 0..9 {
            let fd = (p_plus[row] - p_minus[row]) / (2.0 * h);
            let exact: f64 = tangent[row].iter().zip(direction).map(|(a, b)| a * b).sum();
            assert!(
                (fd - exact).abs() < 1e-4,
                "dt={dt}, row={row}: {fd} versus {exact}"
            );
        }
        assert_eq!(accepted, snapshot);
    }
}

#[test]
fn material_response_separates_incremental_potential_storage_and_relaxation_heat() {
    let law = law();
    let f = [1.25, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let jump = law.step_response(&law.initial_state(), &f, 0.0).unwrap();
    assert_eq!(jump.dissipated_pa, 0.0);
    assert!((jump.incremental_potential_pa - jump.state.stored).abs() < 1e-10);
    let held = law.step_response(&jump.state, &f, 0.1).unwrap();
    assert!(held.dissipated_pa > 0.0);
    assert!(held.state.stored < jump.state.stored);
    assert!(held.incremental_potential_pa > held.state.stored);
    assert!(held.loading_work_pa.abs() < 1e-9);
    assert!(
        (held.loading_work_pa - (held.state.stored - jump.state.stored) - held.dissipated_pa).abs()
            < 1e-8
    );
    assert!(law.freeze(&held.state, -0.1).is_err());
    let mut invalid_history = held.state.clone();
    invalid_history.work = f64::NAN;
    assert!(law.freeze(&invalid_history, 0.1).is_err());
}
