//! G0/G3: reference-aware thin reductions retain strain, work and history.

use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::incompressible_plane_stress::{
    IncompressiblePlaneStressPotential, incompressible_plane_stress,
};
use fs_material::hyper::plane_stress::{PlaneStressSettings, plane_stress, plane_stress_with_seed};
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, Volumetric};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| (0..3).map(|k| a[3 * (i / 3) + k] * b[3 * k + i % 3]).sum())
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} versus {expected:e}; scale={scale:e}, tolerance={tolerance:e}"
    );
}

fn settings() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-12,
        ..PlaneStressSettings::default()
    }
}

fn fibre() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::HolzapfelGasserOgden {
            c10: 5e4,
            k1: 1.2e5,
            k2: 1.1,
            dispersion: 0.05,
            fibres: vec![[0.8, 0.0, 0.6], [0.0, 0.6, 0.8]],
            invariants: FibreInvariants::Full,
            volumetric: Volumetric::Quadratic { kappa: 8e5 },
        },
        4.0,
    )
    .unwrap()
}

fn mooney(kappa: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 6e4,
            c01: 4e4,
            kappa,
        },
        4.0,
    )
    .unwrap()
}

// G = Q^T U has a nonzero lower row and all three shear directions. F=U
// gives the rigid elastic rotation Q, hence zero elastic strain and stress.
fn reference() -> (ReferenceDeformation, [f64; 9]) {
    let (s, c) = (0.6, 0.8);
    let qt = [c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c];
    let u = [1.08, 0.05, 0.08, 0.0, 0.94, -0.06, 0.0, 0.0, 1.03];
    (ReferenceDeformation::try_new(multiply(&qt, &u)).unwrap(), u)
}

#[test]
fn g0_generic_reductions_preserve_native_limits_including_coupled_neo_hookean() {
    let neo = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2e5,
            lambda: 0.0,
        },
        4.0,
    )
    .unwrap();
    for law in [neo, mooney(8e5), fibre()] {
        let a = [1.18, 0.13, -0.04, 0.96];
        let native = law.plane_stress(&a, settings()).unwrap();
        let generic = plane_stress(&law, &a, settings()).unwrap();
        let native_exact = law.incompressible_plane_stress(&a, settings()).unwrap();
        let generic_exact = incompressible_plane_stress(&law, &a, settings()).unwrap();
        for (n, g) in [
            (native.energy_density_pa, generic.energy_density_pa),
            (
                native_exact.energy_density_pa,
                generic_exact.energy_density_pa,
            ),
            (native_exact.pressure_pa, generic_exact.pressure_pa),
        ] {
            close(g, n, 1e6, 2e-11);
        }
        for i in 0..4 {
            close(generic.piola[i], native.piola[i], 1e6, 2e-11);
            close(generic_exact.piola[i], native_exact.piola[i], 1e6, 2e-11);
            for j in 0..4 {
                close(generic.tangent[i][j], native.tangent[i][j], 1e6, 2e-10);
                close(
                    generic_exact.tangent[i][j],
                    native_exact.tangent[i][j],
                    1e6,
                    2e-10,
                );
            }
        }
    }
}

#[test]
fn g3_general_reference_rotation_recovers_stress_free_thickness_and_both_shears() {
    let law = fibre();
    let (reference, u) = reference();
    assert!(reference.deformation_gradient()[6].abs() > 0.1);
    assert!(reference.deformation_gradient()[7].abs() > 0.01);
    let material = ReferencedHyperelastic::new(&law, reference);
    let a = [u[0], u[1], u[3], u[4]];
    let finite = plane_stress(&material, &a, settings()).unwrap();
    let exact = incompressible_plane_stress(&material, &a, settings()).unwrap();
    for (f, energy, p) in [
        (
            finite.deformation_gradient,
            finite.energy_density_pa,
            finite.piola,
        ),
        (
            exact.deformation_gradient,
            exact.energy_density_pa,
            exact.piola,
        ),
    ] {
        for i in 0..9 {
            close(f[i], u[i], 1.0, 2e-9);
        }
        close(energy, 0.0, 1e6, 2e-12);
        for value in p {
            close(value, 0.0, 1e6, 2e-10);
        }
    }
    close(
        determinant(&exact.deformation_gradient),
        reference.volume_ratio(),
        1.0,
        2e-14,
    );
    assert!((reference.volume_ratio() - 1.0).abs() > 0.01);
    let finite_full = material.first_piola(&finite.deformation_gradient).unwrap();
    for i in [2, 5, 8] {
        close(finite_full[i], 0.0, 1e6, 2e-11);
        close(exact.full_piola[i], 0.0, 1e6, 2e-11);
    }
}

#[test]
fn g0_general_reference_reduced_stress_and_tangent_match_independent_work_differences() {
    let law = fibre();
    let material = ReferencedHyperelastic::new(&law, reference().0);
    let a = [1.22, 0.17, -0.04, 0.97];
    let finite = plane_stress(&material, &a, settings()).unwrap();
    let exact = incompressible_plane_stress(&material, &a, settings()).unwrap();
    let h = 2e-6;
    for col in 0..4 {
        let mut plus = a;
        let mut minus = a;
        plus[col] += h;
        minus[col] -= h;
        let fp = plane_stress(&material, &plus, settings()).unwrap();
        let fm = plane_stress(&material, &minus, settings()).unwrap();
        let ep = incompressible_plane_stress(&material, &plus, settings()).unwrap();
        let em = incompressible_plane_stress(&material, &minus, settings()).unwrap();
        close(
            (fp.energy_density_pa - fm.energy_density_pa) / (2.0 * h),
            finite.piola[col],
            1e6,
            3e-8,
        );
        close(
            (ep.energy_density_pa - em.energy_density_pa) / (2.0 * h),
            exact.piola[col],
            1e6,
            3e-8,
        );
        for row in 0..4 {
            close(
                (fp.piola[row] - fm.piola[row]) / (2.0 * h),
                finite.tangent[row][col],
                1e6,
                3e-7,
            );
            close(
                (ep.piola[row] - em.piola[row]) / (2.0 * h),
                exact.tangent[row][col],
                1e6,
                3e-7,
            );
            close(
                finite.tangent[row][col],
                finite.tangent[col][row],
                1e6,
                2e-11,
            );
            close(exact.tangent[row][col], exact.tangent[col][row], 1e6, 2e-11);
        }
    }
}

#[test]
fn g3_referenced_exact_volume_response_is_independent_of_the_bulk_penalty() {
    let a = [1.22, 0.17, -0.04, 0.97];
    let g = reference().0;
    let soft = mooney(3e5);
    let extreme = mooney(1e300);
    let expected =
        incompressible_plane_stress(&ReferencedHyperelastic::new(&soft, g), &a, settings())
            .unwrap();
    let actual =
        incompressible_plane_stress(&ReferencedHyperelastic::new(&extreme, g), &a, settings())
            .unwrap();
    assert_eq!(actual.deformation_gradient, expected.deformation_gradient);
    assert_eq!(actual.energy_density_pa, expected.energy_density_pa);
    assert_eq!(actual.full_piola, expected.full_piola);
    assert_eq!(actual.tangent, expected.tangent);
    close(
        determinant(&actual.deformation_gradient),
        g.volume_ratio(),
        1.0,
        2e-14,
    );
}

#[test]
fn g3_spatial_rotation_preserves_referenced_thin_energy_and_rotates_piola() {
    let law = fibre();
    let material = ReferencedHyperelastic::new(&law, reference().0);
    let a = [1.22, 0.17, -0.04, 0.97];
    let r = [0.8, -0.6, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 1.0];
    let rotated = [
        0.8 * a[0] - 0.6 * a[2],
        0.8 * a[1] - 0.6 * a[3],
        0.6 * a[0] + 0.8 * a[2],
        0.6 * a[1] + 0.8 * a[3],
    ];
    let finite = plane_stress(&material, &a, settings()).unwrap();
    let finite_rotated = plane_stress(&material, &rotated, settings()).unwrap();
    let exact = incompressible_plane_stress(&material, &a, settings()).unwrap();
    let exact_rotated = incompressible_plane_stress(&material, &rotated, settings()).unwrap();
    for (f, rotated_f, p, rotated_p, e, rotated_e) in [
        (
            finite.deformation_gradient,
            finite_rotated.deformation_gradient,
            finite.piola,
            finite_rotated.piola,
            finite.energy_density_pa,
            finite_rotated.energy_density_pa,
        ),
        (
            exact.deformation_gradient,
            exact_rotated.deformation_gradient,
            exact.piola,
            exact_rotated.piola,
            exact.energy_density_pa,
            exact_rotated.energy_density_pa,
        ),
    ] {
        close(rotated_e, e, 1e6, 2e-11);
        let expected_f = multiply(&r, &f);
        for i in 0..9 {
            close(rotated_f[i], expected_f[i], 1.0, 2e-10);
        }
        for j in 0..2 {
            close(rotated_p[j], 0.8 * p[j] - 0.6 * p[2 + j], 1e6, 2e-10);
            close(rotated_p[2 + j], 0.6 * p[j] + 0.8 * p[2 + j], 1e6, 2e-10);
        }
    }
}

#[test]
fn g0_frozen_viscous_thin_reductions_retain_memory_and_native_algorithmic_tangents() {
    let g = reference().0;
    let law = ViscoHyperelastic::new(
        mooney(8e5),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 1.4e5 },
            tau_s: 0.12,
        }],
    )
    .unwrap()
    .with_reference_deformation(g);
    let loading = [1.2, 0.1, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 1.0 / 1.08];
    let old_f = multiply(&loading, g.deformation_gradient());
    let accepted = law.step(&law.initial_state(), &old_f, 0.07).unwrap().0;
    let accepted_remainder = law
        .step_remainder_response(&law.initial_remainder_state(), &old_f, 0.07)
        .unwrap()
        .state;
    let full_snapshot = accepted.clone();
    let remainder_snapshot = accepted_remainder.clone();
    let dt = 0.04;
    let frozen = law.freeze(&accepted, dt).unwrap();
    let frozen_remainder = law.freeze_remainder(&accepted_remainder, dt).unwrap();
    let a = [1.17, 0.1, -0.02, 0.97];
    let finite = plane_stress(&frozen, &a, settings()).unwrap();
    let exact = incompressible_plane_stress(&frozen, &a, settings()).unwrap();
    let remainder = incompressible_plane_stress(&frozen_remainder, &a, settings()).unwrap();
    for i in 0..9 {
        close(exact.full_piola[i], remainder.full_piola[i], 1e6, 1e-11);
    }
    for i in [2, 5, 8] {
        close(
            frozen.first_piola(&finite.deformation_gradient).unwrap()[i],
            0.0,
            1e6,
            2e-11,
        );
        close(exact.full_piola[i], 0.0, 1e6, 2e-11);
    }
    close(
        determinant(&exact.deformation_gradient),
        g.volume_ratio(),
        1.0,
        2e-14,
    );
    let h = 2e-6;
    for col in 0..4 {
        let mut plus = a;
        let mut minus = a;
        plus[col] += h;
        minus[col] -= h;
        let fp = plane_stress(&frozen, &plus, settings()).unwrap();
        let fm = plane_stress(&frozen, &minus, settings()).unwrap();
        let ep = incompressible_plane_stress(&frozen_remainder, &plus, settings()).unwrap();
        let em = incompressible_plane_stress(&frozen_remainder, &minus, settings()).unwrap();
        for row in 0..4 {
            close(
                (fp.piola[row] - fm.piola[row]) / (2.0 * h),
                finite.tangent[row][col],
                1e6,
                4e-7,
            );
            close(
                (ep.piola[row] - em.piola[row]) / (2.0 * h),
                exact.tangent[row][col],
                1e6,
                4e-7,
            );
        }
    }
    assert_eq!(accepted, full_snapshot);
    assert_eq!(accepted_remainder, remainder_snapshot);
    let next = law
        .step_response(&accepted, &finite.deformation_gradient, dt)
        .unwrap();
    for (i, full) in [0, 1, 3, 4].into_iter().enumerate() {
        close(next.first_piola_pa[full], finite.piola[i], 1e6, 2e-11);
    }
    assert!(next.dissipated_pa > 0.0);
    assert_ne!(
        next.state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    close(
        frozen.incompressible_reference_volume_ratio(),
        g.volume_ratio(),
        1.0,
        0.0,
    );
}

#[test]
fn g0_generic_reductions_refuse_inversion_domain_exit_and_exhausted_budgets() {
    let law = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2e5,
            lambda: 3e5,
        },
        1.4,
    )
    .unwrap();
    let material = ReferencedHyperelastic::new(&law, reference().0);
    for a in [
        [1.0, 0.0, 0.0, -1.0],
        [2.0, 0.0, 0.0, 1.0],
        [f64::NAN, 0.0, 0.0, 1.0],
    ] {
        assert!(plane_stress(&material, &a, settings()).is_err());
        assert!(incompressible_plane_stress(&material, &a, settings()).is_err());
    }
    let empty = PlaneStressSettings {
        max_iterations: 0,
        ..settings()
    };
    assert!(plane_stress(&material, &[1.0, 0.0, 0.0, 1.0], empty).is_err());
    assert!(incompressible_plane_stress(&material, &[1.0, 0.0, 0.0, 1.0], empty).is_err());
}

#[test]
fn g0_prescribed_reference_seed_admits_large_natural_expansion_without_widening_elastic_domain() {
    let law = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2e5,
            lambda: 3e5,
        },
        1.15,
    )
    .unwrap();
    let g = ReferenceDeformation::try_new([4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 4.0]).unwrap();
    let material = ReferencedHyperelastic::new(&law, g);
    let a = [4.0, 0.0, 0.0, 4.0];
    assert!(plane_stress(&material, &a, settings()).is_err());
    let response = plane_stress_with_seed(&material, &a, [0.0, 0.0, 4.0], settings()).unwrap();
    assert_eq!(response.deformation_gradient, *g.deformation_gradient());
    assert_eq!(response.piola, [0.0; 4]);
    assert_eq!(response.energy_density_pa, 0.0);
    let exact = incompressible_plane_stress(&material, &a, settings()).unwrap();
    assert_eq!(exact.deformation_gradient, *g.deformation_gradient());
    for seed in [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -4.0],
        [0.0, f64::INFINITY, 4.0],
        [0.0, 0.0, 1.0],
    ] {
        assert!(plane_stress_with_seed(&material, &a, seed, settings()).is_err());
    }

    // The natural-volume guess must not impose incompressibility on a soft
    // bulk law: it is outside this domain, while F22=1 is admissible and exact.
    let soft_bulk = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2e5,
            lambda: 0.0,
        },
        1.8,
    )
    .unwrap();
    let soft = ReferencedHyperelastic::new(&soft_bulk, ReferenceDeformation::identity());
    let a = [1.5, 0.0, 0.0, 1.5];
    let response = plane_stress_with_seed(&soft, &a, [0.0, 0.0, 1.0 / 2.25], settings()).unwrap();
    assert_eq!(response.deformation_gradient[8], 1.0);
}
