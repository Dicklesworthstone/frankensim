//! G0/G1/G3: original finite-strain histories, objective sheet work and tangents.

use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::sheet::SheetLaw;
use fs_material::visco_hyper::membrane::{FrozenViscoMembrane, ViscoMembraneState};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const THICKNESS: f64 = 0.0025;
const TREATMENTS: [VolumetricTreatment; 2] = [
    VolumetricTreatment::OriginalMaterial,
    VolumetricTreatment::ExactIncompressible,
];

fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} versus {b:e}; scale {scale:e}"
    );
}

fn law(kappa: f64) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa,
            },
            4.0,
        )
        .unwrap(),
        vec![
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 1.6e5 },
                tau_s: 0.12,
            },
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 8e4 },
                tau_s: 0.8,
            },
        ],
    )
    .unwrap()
}

fn frozen<'a>(
    law: &'a ViscoHyperelastic,
    state: &'a ViscoMembraneState,
    dt: f64,
) -> FrozenViscoMembrane<'a> {
    FrozenViscoMembrane::new(law, state, dt, THICKNESS, PlaneStressSettings::default()).unwrap()
}

#[test]
fn g0_loading_hold_and_unload_retain_native_memory_and_physical_work() {
    let law = law(8e5);
    let loaded = [0.3, 0.0, 0.04];
    for treatment in TREATMENTS {
        let mut state = ViscoMembraneState::initial(&law, treatment);
        let mut held_stress = Vec::new();
        for (strain, dt) in [
            ([0.0; 3], 0.0),
            (loaded, 0.03),
            (loaded, 0.08),
            (loaded, 0.08),
            ([0.02, 0.15, -0.13], 0.06),
            ([0.0; 3], 0.03),
        ] {
            let before = state.clone();
            let adapter = frozen(&law, &state, dt);
            let trial = adapter.response(strain).unwrap();
            let endpoint = adapter.step_response(strain).unwrap();
            assert_eq!(state, before);
            assert_eq!(trial.energy_j_m2, endpoint.sheet.energy_j_m2);
            assert_eq!(trial.resultants_n_m, endpoint.sheet.resultants_n_m);
            assert_eq!(trial.tangent_n_m, endpoint.sheet.tangent_n_m);
            close(
                endpoint.loading_work_j_m2,
                endpoint.stored_energy_j_m2 - THICKNESS * before.stored_pa()
                    + endpoint.dissipated_j_m2,
                1e3,
                5e-12,
            );
            assert!(endpoint.dissipated_j_m2 >= 0.0);
            assert!(endpoint.state.dissipated_pa() >= before.dissipated_pa());
            for i in [2, 5, 8] {
                close(endpoint.first_piola_pa[i], 0.0, 1e6, 2e-10);
            }
            if strain == loaded {
                held_stress.push(endpoint.sheet.resultants_n_m[0]);
                assert!(endpoint.dissipated_j_m2 > 0.0);
                assert!((endpoint.sheet.energy_j_m2 - endpoint.stored_energy_j_m2).abs() > 0.01);
            }
            state = endpoint.state;
        }
        assert!(held_stress[2] < held_stress[1] && held_stress[1] < held_stress[0]);
        close(
            state.work_pa() + state.reference_work_pa(),
            state.stored_pa() + state.dissipated_pa(),
            1e6,
            5e-12,
        );
        assert!(state.dissipated_pa() > 0.0);
    }
}

#[test]
fn g1_green_strain_tangent_after_nonproportional_loading_is_the_original_potential_hessian() {
    let law = law(8e5);
    for treatment in TREATMENTS {
        let mut state = ViscoMembraneState::initial(&law, treatment);
        for strain in [[0.2, -0.03, 0.1], [0.03, 0.18, -0.12]] {
            state = frozen(&law, &state, 0.04)
                .step_response(strain)
                .unwrap()
                .state;
        }
        let before = state.clone();
        let strain = [0.12, 0.05, 0.22];
        let adapter = frozen(&law, &state, 0.08);
        let response = adapter.response(strain).unwrap();
        let h = 2e-6;
        for col in 0..3 {
            let mut plus = strain;
            let mut minus = strain;
            plus[col] += h;
            minus[col] -= h;
            let plus = adapter.response(plus).unwrap();
            let minus = adapter.response(minus).unwrap();
            close(
                (plus.energy_j_m2 - minus.energy_j_m2) / (2.0 * h),
                response.resultants_n_m[col],
                1e3,
                2e-7,
            );
            for row in 0..3 {
                close(
                    (plus.resultants_n_m[row] - minus.resultants_n_m[row]) / (2.0 * h),
                    response.tangent_n_m[3 * row + col],
                    1e4,
                    2e-7,
                );
                close(
                    response.tangent_n_m[3 * row + col],
                    response.tangent_n_m[3 * col + row],
                    1e4,
                    2e-12,
                );
            }
        }
        assert_eq!(state, before);
        let virgin = ViscoMembraneState::initial(&law, treatment);
        let no_history = frozen(&law, &virgin, 0.08).response(strain).unwrap();
        assert!((response.resultants_n_m[2] - no_history.resultants_n_m[2]).abs() > 1.0);
    }
}

#[test]
fn g1_exact_elastic_limit_matches_an_independent_green_metric_energy() {
    let mu = 2e5;
    let law = ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 0.5 * mu,
                c01: 0.0,
                kappa: 1e300,
            },
            4.0,
        )
        .unwrap(),
        vec![],
    )
    .unwrap();
    let state = ViscoMembraneState::initial(&law, VolumetricTreatment::ExactIncompressible);
    let strain = [0.22, -0.06, 0.17];
    let response = frozen(&law, &state, 3.0).step_response(strain).unwrap();
    let a = 1.0 + 2.0 * strain[0];
    let b = 1.0 + 2.0 * strain[1];
    let g = strain[2];
    let determinant = a * b - g * g;
    close(
        response.stored_energy_j_m2,
        THICKNESS * mu * 0.5 * (a + b + 1.0 / determinant - 3.0),
        1e3,
        1e-12,
    );
    let expected = [
        THICKNESS * mu * (1.0 - b / determinant.powi(2)),
        THICKNESS * mu * (1.0 - a / determinant.powi(2)),
        THICKNESS * mu * g / determinant.powi(2),
    ];
    for i in 0..3 {
        close(response.sheet.resultants_n_m[i], expected[i], 1e3, 1e-12);
    }
    close(
        response.state.deformation_gradient()[8],
        1.0 / determinant.sqrt(),
        1.0,
        1e-14,
    );
    assert_eq!(response.dissipated_j_m2, 0.0);
}

#[test]
fn g3_exact_volume_retains_natural_reference_and_is_independent_of_unused_bulk() {
    let reference =
        ReferenceDeformation::try_new([1.08, 0.04, 0.16, 0.0, 0.95, -0.09, 0.0, 0.0, 1.03])
            .unwrap();
    let stiff = law(1e300).with_reference_deformation(reference.clone());
    let law = law(8e5).with_reference_deformation(reference);
    let mut state = ViscoMembraneState::initial(&law, VolumetricTreatment::ExactIncompressible);
    let mut stiff_state =
        ViscoMembraneState::initial(&stiff, VolumetricTreatment::ExactIncompressible);
    for strain in [[0.2, -0.03, 0.12], [0.01, 0.19, -0.09], [0.0; 3]] {
        let endpoint = frozen(&law, &state, 0.08).step_response(strain).unwrap();
        let other = frozen(&stiff, &stiff_state, 0.08)
            .step_response(strain)
            .unwrap();
        assert_eq!(endpoint.state, other.state);
        assert_eq!(endpoint.sheet.energy_j_m2, other.sheet.energy_j_m2);
        assert_eq!(endpoint.sheet.resultants_n_m, other.sheet.resultants_n_m);
        assert_eq!(endpoint.sheet.tangent_n_m, other.sheet.tangent_n_m);
        assert_eq!(endpoint.first_piola_pa, other.first_piola_pa);
        let f = endpoint.state.deformation_gradient();
        close(
            f[0] * f[4] * f[8],
            law.reference_deformation().volume_ratio(),
            1.0,
            1e-14,
        );
        assert!(endpoint.pressure_pa.is_some());
        state = endpoint.state;
        stiff_state = other.state;
    }
}

#[test]
fn g0_refused_metric_interval_and_thickness_leave_accepted_history_unchanged() {
    let law = law(8e5);
    for treatment in TREATMENTS {
        let state = ViscoMembraneState::initial(&law, treatment);
        let before = state.clone();
        let settings = PlaneStressSettings::default();
        for dt in [-0.1, f64::NAN, f64::INFINITY] {
            assert!(FrozenViscoMembrane::new(&law, &state, dt, THICKNESS, settings).is_err());
        }
        for thickness in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(FrozenViscoMembrane::new(&law, &state, 0.1, thickness, settings).is_err());
        }
        let adapter = frozen(&law, &state, 0.1);
        for strain in [[-0.5, 0.0, 0.0], [0.0, 0.0, 2.0], [f64::NAN, 0.0, 0.0]] {
            assert!(adapter.response(strain).is_err());
            assert!(adapter.step_response(strain).is_err());
        }
        assert_eq!(state, before);
    }
}
