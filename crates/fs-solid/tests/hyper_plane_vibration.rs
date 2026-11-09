//! Independent discrete spectra about loaded native plane equilibria.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper2d::vibration::{
    PlaneReferenceDensity, PlaneVibration, PlaneVibrationError, PlaneVibrationFormulation,
    PlaneVibrationSettings, VibrationSettings,
};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 73,
                    kernel_id: 59,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn settings() -> PlaneVibrationSettings {
    PlaneVibrationSettings {
        spectrum: VibrationSettings {
            maximum_free_dofs: 32,
            relative_residual_tolerance: 1e-8,
        },
        maximum_mass_assembly_bytes: 1024 * 1024,
    }
}

fn problem<'a>(mesh: &'a Mesh2, material: &'a Hyperelastic) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings {
            tol: 2e-9,
            max_iters: 40,
            load_steps: 1,
            max_backtracks: 30,
        },
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-8 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn check_modes(spectrum: &PlaneVibration<'_, '_>, expected: &[f64], diagonal_mass: f64) {
    assert_eq!(spectrum.modes().len(), expected.len());
    let n = expected.len();
    for (mode, &lambda) in spectrum.modes().iter().zip(expected) {
        near(mode.lambda, lambda);
        assert!(mode.residual <= 1e-8 * expected.iter().map(|v| v.abs()).fold(0.0, f64::max));
        near(mode.phi.iter().map(|v| v * v * diagonal_mass).sum(), 1.0);
        for i in 0..n {
            let residual: f64 = (0..n)
                .map(|j| {
                    (spectrum.stiffness_n_m()[i * n + j]
                        - mode.lambda * spectrum.mass_kg()[i * n + j])
                        * mode.phi[j]
                })
                .sum();
            near(residual, 0.0);
        }
    }
    for i in 0..n {
        for j in 0..n {
            near(
                spectrum.mass_kg()[i * n + j],
                if i == j { diagonal_mass } else { 0.0 },
            );
        }
    }
}

#[test]
fn finite_and_exact_plane_stress_modes_use_relaxed_thickness_and_original_mass() {
    with_cx(|cx, _| {
        let (mu, rho) = (2.0, 3.0);
        for mesh in [
            Mesh2::quads(1.0, 1.0, 2, 2),
            Mesh2::triangles(1.0, 1.0, 2, 2),
        ] {
            let p1 = mesh.elems[0].len() == 3;
            for stretch in [1.0_f64, 1.3] {
                let boundary = |x: f64, y: f64| [(stretch - 1.0) * x, (stretch - 1.0) * y];
                for lambda in [40.0, 4000.0] {
                    let material =
                        Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0)
                            .unwrap();
                    let mut p = problem(&mesh, &material);
                    p.dirichlet = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top]
                        .into_iter()
                        .map(|patch| (patch, &boundary as _))
                        .collect();
                    let density = vec![rho; mesh.elems.len()];
                    let local = PlaneStressSettings::default();
                    // Independent traction equation: mu(t²-1)+lambda ln(s²t)=0.
                    let (mut lo, mut hi) = (1.0 / (stretch * stretch), 1.0);
                    for _ in 0..80 {
                        let t = 0.5 * (lo + hi);
                        if mu * (t * t - 1.0) + lambda * (stretch * stretch * t).ln() > 0.0 {
                            hi = t;
                        } else {
                            lo = t;
                        }
                    }
                    let t = 0.5 * (lo + hi);
                    for (form, coupling) in [
                        (
                            PlaneVibrationFormulation::PlaneStress(local),
                            mu * t * t / (stretch * stretch)
                                * (1.0 + 2.0 * lambda / (lambda + 2.0 * mu * t * t)),
                        ),
                        (
                            PlaneVibrationFormulation::IncompressiblePlaneStress(local),
                            3.0 * mu / stretch.powi(6),
                        ),
                    ] {
                        // Native shape integrals give two coupled P1 eigenvalues
                        // or two repeated Q1 eigenvalues, without an eigensolver oracle.
                        let expected = if p1 {
                            [
                                8.0 * (4.0 * mu + coupling) / rho,
                                8.0 * (4.0 * mu + 3.0 * coupling) / rho,
                            ]
                        } else {
                            [12.0 * (2.0 * mu + coupling) / rho; 2]
                        };
                        let h = 0.02;
                        let spectrum = p
                            .solve_vibration(
                                form,
                                h,
                                PlaneReferenceDensity::PerElement(&density),
                                settings(),
                                cx,
                            )
                            .unwrap();
                        assert_eq!(spectrum.free_dofs(), &[8, 9]);
                        check_modes(&spectrum, &expected, rho * h / if p1 { 8.0 } else { 9.0 });
                        near(spectrum.reference_mass_kg(), rho * h);
                        for (point, u) in mesh.nodes.iter().zip(spectrum.displacement_m()) {
                            near(u[0], (stretch - 1.0) * point[0]);
                            near(u[1], (stretch - 1.0) * point[1]);
                        }
                    }
                }
            }
        }
    });
}

#[test]
fn mean_dilatation_keeps_its_condensed_hessian_and_signed_unstable_modes() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
        let (mu, lambda, rho, h) = (2.0, 100.0, 3.0, 0.02);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let density = vec![rho; mesh.elems.len()];
        for stretch in [1.0_f64, 1.3, 2.0] {
            let boundary = |x: f64, y: f64| [(stretch - 1.0) * x, (stretch - 1.0) * y];
            let mut p = problem(&mesh, &material);
            p.dirichlet = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top]
                .into_iter()
                .map(|patch| (patch, &boundary as _))
                .collect();
            let initial: Vec<_> = mesh
                .nodes
                .iter()
                .map(|&[x, y]| [0.05 * x, 0.05 * y])
                .collect();
            let spectrum = p
                .solve_vibration_from(
                    &initial,
                    PlaneVibrationFormulation::MeanDilatation,
                    h,
                    PlaneReferenceDensity::PerElement(&density),
                    settings(),
                    cx,
                )
                .unwrap();
            // Integral grad(N)^2=4/3 and integral N²=1/9 for the interior
            // Q1 hat. The condensed bulk term uses each CELL'S mean gradient.
            let expected = (mu * (24.0 + 12.0 / stretch.powi(2))
                + 9.0 * lambda * (1.0 - 2.0 * stretch.ln()) / stretch.powi(2))
                / rho;
            check_modes(&spectrum, &[expected; 2], rho * h / 9.0);
            if stretch == 2.0 {
                assert!(spectrum.modes().iter().all(|mode| mode.interval.1 < 0.0));
            }
        }
    });
}

#[test]
fn reference_inventory_scaling_and_free_rigid_modes_remain_physical() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 40.0,
            },
            4.0,
        )
        .unwrap();
        let p = problem(&mesh, &material);
        let form =
            PlaneVibrationFormulation::IncompressiblePlaneStress(PlaneStressSettings::default());
        let cells = [3.0];
        let nodes = [3.0; 4];
        let base = p
            .solve_vibration(
                form,
                0.02,
                PlaneReferenceDensity::PerElement(&cells),
                settings(),
                cx,
            )
            .unwrap();
        let nodal = p
            .solve_vibration(
                form,
                0.02,
                PlaneReferenceDensity::Nodal(&nodes),
                settings(),
                cx,
            )
            .unwrap();
        assert_eq!(base.mass_kg(), nodal.mass_kg());
        let changed_density = [12.0; 4];
        let scaled = p
            .solve_vibration(
                form,
                0.2,
                PlaneReferenceDensity::Nodal(&changed_density),
                settings(),
                cx,
            )
            .unwrap();
        assert_eq!(base.modes().len(), 8);
        let largest = base
            .modes()
            .iter()
            .map(|mode| mode.lambda.abs())
            .fold(0.0, f64::max);
        assert_eq!(
            base.modes()
                .iter()
                .filter(|mode| mode.lambda.abs() < 1e-10 * largest)
                .count(),
            3
        );
        for (before, after) in base.modes().iter().zip(scaled.modes()) {
            near(after.lambda, before.lambda / 4.0);
        }
        for (&before, &after) in base.stiffness_n_m().iter().zip(scaled.stiffness_n_m()) {
            near(after, 10.0 * before);
        }
        for (&before, &after) in base.mass_kg().iter().zip(scaled.mass_kg()) {
            near(after, 40.0 * before);
        }
        let acceleration = vec![[2.0, -3.0]; 4];
        let forces = scaled
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        near(
            forces.iter().step_by(2).sum(),
            scaled.reference_mass_kg() * 2.0,
        );
        near(
            forces.iter().skip(1).step_by(2).sum(),
            -scaled.reference_mass_kg() * 3.0,
        );
        let mut limited = settings();
        limited.spectrum.maximum_free_dofs = 7;
        assert!(
            p.solve_vibration(
                form,
                0.02,
                PlaneReferenceDensity::PerElement(&cells),
                limited,
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn incremental_materials_and_cancelled_work_are_refused_before_publication() {
    use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
    with_cx(|cx, gate| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let law = ViscoHyperelastic::new(
            Hyperelastic::new(
                HyperelasticModel::NeoHookean {
                    mu: 2.0,
                    lambda: 40.0,
                },
                4.0,
            )
            .unwrap(),
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 1.0 },
                tau_s: 0.2,
            }],
        )
        .unwrap();
        let accepted = law.initial_state();
        let frozen = law.freeze(&accepted, 0.1).unwrap();
        let mut p = HyperProblem {
            mesh: &mesh,
            material: &frozen,
            dirichlet: vec![],
            traction: vec![],
            settings: NewtonSettings::default(),
        };
        let density = [3.0];
        let form = PlaneVibrationFormulation::PlaneStress(PlaneStressSettings::default());
        assert!(matches!(
            p.solve_vibration(
                form,
                0.02,
                PlaneReferenceDensity::PerElement(&density),
                settings(),
                cx
            ),
            Err(PlaneVibrationError::Solid(SolidError::InvalidInput { .. }))
        ));
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 40.0,
            },
            4.0,
        )
        .unwrap();
        p.material = &material;
        gate.request();
        assert!(matches!(
            p.solve_vibration(
                form,
                0.02,
                PlaneReferenceDensity::PerElement(&density),
                settings(),
                cx
            ),
            Err(PlaneVibrationError::Solid(SolidError::Cancelled))
        ));
    });
}

#[test]
fn held_point_references_keep_their_stiffness_and_original_mass() {
    use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
    use fs_solid::hyper2d::{HyperelasticMaterials, PlaneMaterial};

    struct PointLaws<'a>(Vec<ReferencedHyperelastic<'a, Hyperelastic>>);
    impl HyperelasticMaterials for PointLaws<'_> {
        fn material_count(&self) -> Option<usize> {
            Some(self.0.len() / 4)
        }
        fn material_at(&self, _: usize) -> Option<&dyn PlaneMaterial> {
            None
        }
        fn material_point_count(&self, _: usize) -> Option<usize> {
            Some(4)
        }
        fn material_at_point(&self, element: usize, point: usize) -> Option<&dyn PlaneMaterial> {
            self.0
                .get(4 * element + point)
                .map(|law| law as &dyn PlaneMaterial)
        }
    }
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
        let (mu, lambda, rho, h, stretch) = (2.0, 40.0, 3.0, 0.02, 1.2);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        // Each point retains its own reference owner. Element-only lookup
        // intentionally fails; a consumer cannot silently use a bare/first law.
        let field = PointLaws(
            (0..16)
                .map(|point| {
                    let thickness = if (point / 4) % 2 == 0 { 1.0 } else { 2.0 };
                    ReferencedHyperelastic::new(
                        &material,
                        ReferenceDeformation::try_new([
                            stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, thickness,
                        ])
                        .unwrap(),
                    )
                })
                .collect(),
        );
        let boundary = |x: f64, y: f64| [(stretch - 1.0) * x, (stretch - 1.0) * y];
        let mut p = problem(&mesh, &material);
        p.material = &field;
        p.dirichlet = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top]
            .into_iter()
            .map(|patch| (patch, &boundary as _))
            .collect();
        let density = vec![rho; mesh.elems.len()];
        let local = PlaneStressSettings::default();
        let spectrum = p
            .solve_vibration(
                PlaneVibrationFormulation::PlaneStress(local),
                h,
                PlaneReferenceDensity::PerElement(&density),
                settings(),
                cx,
            )
            .unwrap();
        // W0=det(G) We(F G^-1): at F=G the in-plane Hessian
        // scales by Gzz, not det(G). These four cells have mean Gzz=1.5;
        // their opposite cross terms cancel. Original inertia does not scale.
        let lambda_plane = 2.0 * mu * lambda / (lambda + 2.0 * mu);
        let expected = 1.5 * 12.0 * (lambda_plane + 3.0 * mu) / rho;
        check_modes(&spectrum, &[expected; 2], rho * h / 9.0);
        near(spectrum.reference_mass_kg(), rho * h);
        for (u, &[x, y]) in spectrum.displacement_m().iter().zip(&mesh.nodes) {
            near(u[0], (stretch - 1.0) * x);
            near(u[1], (stretch - 1.0) * y);
        }
        let flat: Vec<_> = spectrum
            .displacement_m()
            .iter()
            .flatten()
            .copied()
            .collect();
        for (point, state) in p
            .plane_stress_states(&flat, local, cx)
            .unwrap()
            .iter()
            .enumerate()
        {
            near(
                state.deformation_gradient[8],
                if (point / 4) % 2 == 0 { 1.0 } else { 2.0 },
            );
        }
    });
}
