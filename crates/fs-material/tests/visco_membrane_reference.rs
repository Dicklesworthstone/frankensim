//! Native natural-strain events preserve membrane memory and physical work.

use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::sheet::SheetLaw;
use fs_material::visco_hyper::membrane::{
    FrozenViscoMembrane, ViscoMembraneReferenceChange, ViscoMembraneState,
};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const THICKNESS: f64 = 0.0025;
const STRAIN: [f64; 3] = [0.12, 0.05, 0.22];
const TREATMENTS: [VolumetricTreatment; 2] = [
    VolumetricTreatment::OriginalMaterial,
    VolumetricTreatment::ExactIncompressible,
];

fn settings() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-13,
        ..Default::default()
    }
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

fn close(actual: f64, expected: f64, scale: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= scale * relative,
        "{actual:e} != {expected:e} at scale {scale:e}"
    );
}

fn loaded(law: &ViscoHyperelastic, treatment: VolumetricTreatment) -> ViscoMembraneState {
    let mut state = ViscoMembraneState::initial(law, treatment);
    for strain in [[0.2, -0.03, 0.1], [0.03, 0.18, -0.12], STRAIN] {
        state = FrozenViscoMembrane::new(law, &state, 0.04, THICKNESS, settings())
            .unwrap()
            .step_response(strain)
            .unwrap()
            .state;
    }
    state
}

fn metrics(state: &ViscoMembraneState) -> &[[f64; 6]] {
    match state {
        ViscoMembraneState::Full(state) => &state.inverse_viscous_cauchy_green,
        ViscoMembraneState::Incompressible(state) => &state.inverse_viscous_cauchy_green,
    }
}

fn reference(law: ViscoHyperelastic, stretches: [f64; 3]) -> ViscoHyperelastic {
    law.with_reference_deformation(
        ReferenceDeformation::from_principal_stretches(
            stretches,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
        .unwrap(),
    )
}

fn check_work(change: &ViscoMembraneReferenceChange, before: &ViscoMembraneState) {
    close(
        change.stored_before_j_m2,
        THICKNESS * before.stored_pa(),
        1e3,
        1e-13,
    );
    close(
        change.stored_after_j_m2,
        THICKNESS * change.state.stored_pa(),
        1e3,
        1e-13,
    );
    close(
        change.reference_work_j_m2 + change.native_transverse_loading_work_j_m2,
        change.stored_after_j_m2 - change.stored_before_j_m2,
        1e3,
        2e-11,
    );
    close(
        change.transverse_endpoint_work_j_m2 + change.constraint_reference_work_j_m2
            - change.native_transverse_loading_work_j_m2,
        change.endpoint_work_defect_j_m2,
        1e3,
        1e-13,
    );
    close(
        THICKNESS * (change.state.reference_work_pa() - before.reference_work_pa()),
        change.reference_work_j_m2,
        1e3,
        2e-11,
    );
    close(
        THICKNESS * (change.state.work_pa() - before.work_pa()),
        change.native_transverse_loading_work_j_m2,
        1e3,
        2e-11,
    );
    assert_eq!(change.state.dissipated_pa(), before.dissipated_pa());
    assert_eq!(metrics(&change.state), metrics(before));
    for i in [0, 1, 3, 4] {
        assert_eq!(
            change.state.deformation_gradient()[i],
            before.deformation_gradient()[i]
        );
    }
    for i in [2, 5, 8] {
        close(change.first_piola_pa[i], 0.0, 1e6, 2e-11);
    }
    let endpoint: f64 = [2, 5, 8]
        .iter()
        .map(|&i| {
            THICKNESS
                * change.first_piola_pa[i]
                * (change.state.deformation_gradient()[i] - before.deformation_gradient()[i])
        })
        .sum();
    close(change.transverse_endpoint_work_j_m2, endpoint, 1e3, 1e-13);
}

#[test]
fn changing_natural_strain_preserves_nonproportional_memory_and_original_work() {
    let previous = law(8e5);
    let target = reference(previous.clone(), [1.04, 0.98, 1.03]);
    for treatment in TREATMENTS {
        let accepted = loaded(&previous, treatment);
        let before = accepted.clone();
        let change = target
            .change_membrane_reference_from(&previous, &accepted, STRAIN, THICKNESS, settings())
            .unwrap();
        assert_eq!(accepted, before);
        check_work(&change, &before);
        assert!(change.reference_work_j_m2.abs() > 0.01);
        if let Some(pressure) = change.pressure_pa {
            let f = change.state.deformation_gradient();
            let target_volume = target.reference_deformation().volume_ratio();
            close(f[0] * f[4] * f[8], target_volume, 2.0, 1e-13);
            close(
                change.constraint_reference_work_j_m2,
                THICKNESS
                    * pressure
                    * (target_volume - previous.reference_deformation().volume_ratio()),
                1e3,
                1e-13,
            );
            assert!(change.constraint_reference_work_j_m2.abs() > 0.01);
        } else {
            assert_eq!(change.constraint_reference_work_j_m2, 0.0);
        }
        // Subsequent physical time continues the old metrics under the new G.
        let held = FrozenViscoMembrane::new(&target, &change.state, 0.05, THICKNESS, settings())
            .unwrap()
            .step_response(STRAIN)
            .unwrap();
        assert!(held.dissipated_j_m2 > 0.0);
        assert!(held.state.dissipated_pa() > accepted.dissipated_pa());
    }
}

#[test]
fn changed_reference_endpoint_has_the_original_consistent_green_strain_tangent() {
    let previous = law(8e5);
    let target = reference(previous.clone(), [1.04, 0.98, 1.03]);
    for treatment in TREATMENTS {
        let accepted = loaded(&previous, treatment);
        let change = target
            .change_membrane_reference_from(&previous, &accepted, STRAIN, THICKNESS, settings())
            .unwrap();
        let frozen =
            FrozenViscoMembrane::new(&target, &change.state, 0.0, THICKNESS, settings()).unwrap();
        let h = 2e-6;
        for col in 0..3 {
            let mut plus = STRAIN;
            let mut minus = STRAIN;
            plus[col] += h;
            minus[col] -= h;
            let plus = frozen.response(plus).unwrap();
            let minus = frozen.response(minus).unwrap();
            close(
                (plus.energy_j_m2 - minus.energy_j_m2) / (2.0 * h),
                change.sheet.resultants_n_m[col],
                1e3,
                2e-7,
            );
            for row in 0..3 {
                close(
                    (plus.resultants_n_m[row] - minus.resultants_n_m[row]) / (2.0 * h),
                    change.sheet.tangent_n_m[3 * row + col],
                    1e4,
                    2e-7,
                );
            }
        }
        close(
            change.sheet.energy_j_m2,
            change.stored_after_j_m2,
            1e3,
            2e-11,
        );
    }
}

#[test]
fn kinetic_shift_alone_has_no_reference_work_and_keeps_the_accepted_metrics() {
    let previous = law(8e5);
    let mut target = previous.clone();
    for branch in &mut target.branches {
        branch.tau_s *= 0.5;
    }
    for treatment in TREATMENTS {
        let accepted = loaded(&previous, treatment);
        let change = target
            .change_membrane_reference_from(&previous, &accepted, STRAIN, THICKNESS, settings())
            .unwrap();
        check_work(&change, &accepted);
        assert_eq!(change.reference_work_j_m2, 0.0);
        assert_eq!(change.constraint_reference_work_j_m2, 0.0);
        close(change.native_transverse_loading_work_j_m2, 0.0, 1e3, 1e-13);
        let slow = FrozenViscoMembrane::new(&previous, &accepted, 0.08, THICKNESS, settings())
            .unwrap()
            .step_response(STRAIN)
            .unwrap();
        let fast = FrozenViscoMembrane::new(&target, &change.state, 0.08, THICKNESS, settings())
            .unwrap()
            .step_response(STRAIN)
            .unwrap();
        assert!(fast.dissipated_j_m2 > slow.dissipated_j_m2);
        assert!(fast.stored_energy_j_m2 < slow.stored_energy_j_m2);
    }
}

#[test]
fn exact_volume_reference_change_is_independent_of_unused_bulk_stiffness() {
    let mut responses = Vec::new();
    for kappa in [8e5, 1e300] {
        let previous = law(kappa);
        let target = reference(previous.clone(), [1.2, 0.9, 1.1]);
        let accepted = loaded(&previous, VolumetricTreatment::ExactIncompressible);
        let change = target
            .change_membrane_reference_from(&previous, &accepted, STRAIN, THICKNESS, settings())
            .unwrap();
        check_work(&change, &accepted);
        responses.push(change);
    }
    let [a, b] = [&responses[0], &responses[1]];
    assert_eq!(
        a.state.deformation_gradient(),
        b.state.deformation_gradient()
    );
    assert_eq!(metrics(&a.state), metrics(&b.state));
    assert_eq!(a.sheet, b.sheet);
    assert_eq!(a.first_piola_pa, b.first_piola_pa);
    assert_eq!(a.pressure_pa, b.pressure_pa);
    assert_eq!(a.reference_work_j_m2, b.reference_work_j_m2);
    assert_eq!(
        a.constraint_reference_work_j_m2,
        b.constraint_reference_work_j_m2
    );
    assert_eq!(
        a.native_transverse_loading_work_j_m2,
        b.native_transverse_loading_work_j_m2
    );
    assert_eq!(a.endpoint_work_defect_j_m2, b.endpoint_work_defect_j_m2);
}

#[test]
fn inadmissible_reference_events_leave_installed_strain_and_memory_unchanged() {
    let previous = law(8e5);
    for treatment in TREATMENTS {
        let accepted = loaded(&previous, treatment);
        let before = accepted.clone();
        let target = reference(previous.clone(), [1.04, 0.98, 1.03]);
        let mut changed_strain = STRAIN;
        changed_strain[2] += 1e-4;
        assert!(
            target
                .change_membrane_reference_from(
                    &previous,
                    &accepted,
                    changed_strain,
                    THICKNESS,
                    settings(),
                )
                .is_err()
        );
        for thickness in [0.0, -0.001, f64::NAN, f64::INFINITY, 1e308] {
            assert!(
                target
                    .change_membrane_reference_from(
                        &previous,
                        &accepted,
                        STRAIN,
                        thickness,
                        settings(),
                    )
                    .is_err()
            );
        }
        let mut replacement = target.clone();
        replacement.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 2e5 };
        assert!(replacement.change_membrane_reference_from(
            &previous, &accepted, STRAIN, THICKNESS, settings(),
        ).is_err());
        let outside = reference(previous.clone(), [10.0; 3]);
        assert!(outside.change_membrane_reference_from(
            &previous, &accepted, STRAIN, THICKNESS, settings(),
        ).is_err());
        assert_eq!(accepted, before);
    }
}
