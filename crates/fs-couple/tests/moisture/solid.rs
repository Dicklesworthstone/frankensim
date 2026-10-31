//! G1 compatible swelling and G3 data-only changes on the native solid owner.
//! Every coefficient here is authored for verification, not wood validation.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::solid::*;
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig, source::SolidStiffnessModel};

fn axes() -> HygroscopicAxes {
    HygroscopicAxes {
        principal_to_mesh: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        principal_frame: ContentHash([3; 32]),
        mesh_frame: ContentHash([4; 32]),
    }
}
fn controls() -> MoistureSolidControls {
    MoistureSolidControls {
        assembly: TetAssemblyBudget {
            maximum_nodes: 256,
            maximum_tetrahedra: 512,
            maximum_free_dofs: 768,
            minimum_scaled_jacobian: 1e-8,
        },
        equilibrium: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-10,
            maximum_iterations: 4096,
        },
        maximum_absolute_mandel_strain: 0.02,
        maximum_nodal_displacement_m: 0.02,
        mass_tolerance_kg: 1e-12,
    }
}
fn bind(
    slab: &MoistureSlab,
    model: SolidStiffnessModel,
    axes: HygroscopicAxes,
    fixed: Vec<usize>,
    cx: &Cx<'_>,
) -> MoistureSolid {
    MoistureSolid::new(
        slab,
        0.01,
        model,
        MaterialPropertySelection::SingleClaimOnly,
        MaterialPropertySelection::SingleClaimOnly,
        axes,
        fixed,
        controls(),
        cx,
    )
    .unwrap()
}
fn free_supports() -> Vec<usize> {
    vec![0, 1, 2, 3, 5, 9]
}

#[test]
fn finite_solution_activity_drives_native_swelling_and_late_refusal_is_atomic() {
    with_cx(|cx| {
        let volume = 1e-4;
        let original = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let binding = bind(
            &original,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        );
        let faces = [
            Some(atmosphere::sorption(0.2)),
            Some(atmosphere::sorption(0.2)),
        ];
        let pure = solution::reservoir(volume, vec![0.0]);
        let saline = solution::reservoir(volume, vec![0.03]);
        let (pure_water, pure_solid) = binding
            .advance_finite_solution(
                &pure,
                &original,
                0.25,
                &faces,
                config(),
                solution::controls(volume),
                cx,
            )
            .unwrap();
        let (salt_water, salt_solid) = binding
            .advance_finite_solution(
                &saline,
                &original,
                0.25,
                &faces,
                config(),
                solution::controls(volume),
                cx,
            )
            .unwrap();
        let activity = salt_water.reservoir().solution().water_activity().unwrap();
        // Two half-cell faces give alpha = 4*dt*D/L² = 1. Source-integrated
        // beta is constant, so compatible free strain is beta*(M_end-M_start).
        let expected_moisture = (0.1 + 0.2 * activity) / 2.0;
        near(
            salt_water.moisture().moisture()[0],
            expected_moisture,
            2e-13,
        );
        let strain = [0.01, 0.04, 0.02].map(|beta| beta * (expected_moisture - 0.1));
        for (x, u) in binding
            .reference_nodes_m()
            .iter()
            .zip(salt_solid.solution().displacement_m())
        {
            for k in 0..3 {
                near(u[k], strain[k] * x[k], 1e-10);
            }
        }
        for element in salt_solid.stress().elements() {
            for stress in element.stress_mandel_pa() {
                near(*stress, 0.0, 0.1);
            }
        }
        assert!(salt_water.step().water_to_solid_kg < pure_water.step().water_to_solid_kg);
        assert!(salt_solid.wet_mass_kg() < pure_solid.wet_mass_kg());
        assert!(
            salt_solid
                .geometry()
                .element_volumes_m3()
                .iter()
                .sum::<f64>()
                < pure_solid
                    .geometry()
                    .element_volumes_m3()
                    .iter()
                    .sum::<f64>()
        );
        near(
            salt_solid.wet_mass_kg(),
            salt_water.moisture().wet_mass_kg(),
            1e-13,
        );
        assert_eq!(salt_solid.time_s(), salt_water.moisture().time_s());

        let limited = slab(card([0.01; 2], 0.12), 0.1, vec![0.1], cx);
        let limited_solid = bind(
            &limited,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        );
        // Native diffusion accepts the proposed water exchange, but the actual
        // endpoint elastic source is outside its admitted moisture range.
        let water_only = saline
            .advance(
                &limited,
                0.25,
                &faces,
                config(),
                solution::controls(volume),
                cx,
            )
            .unwrap();
        assert!(water_only.moisture().moisture()[0] > 0.12);
        assert!(
            limited_solid
                .advance_finite_solution(
                    &saline,
                    &limited,
                    0.25,
                    &faces,
                    config(),
                    solution::controls(volume),
                    cx,
                )
                .is_err()
        );
        assert_eq!(saline.solution().solvent_volume_m3(), volume);
        assert_eq!(saline.solution().solute_moles(), &[0.03]);
        assert_eq!(limited.time_s(), 0.0);
        assert_eq!(limited.moisture(), &[0.1]);
    });
}

#[test]
fn transported_moisture_drives_free_shape_wet_mass_and_current_density() {
    with_cx(|cx| {
        let original = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let binding = bind(
            &original,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        );
        let step = binding
            .advance(&original, 0.25, &faces(0.2), config(), cx)
            .unwrap();
        // One-cell backward Euler gives .15; beta*(.15-.1) is directional strain.
        near(step.moisture.moisture()[0], 0.15, 1e-13);
        let strain = [0.0005, 0.002, 0.001];
        for (x, u) in binding
            .reference_nodes_m()
            .iter()
            .zip(step.solid.solution().displacement_m())
        {
            for k in 0..3 {
                near(u[k], strain[k] * x[k], 1e-10);
            }
        }
        for element in step.solid.stress().elements() {
            for stress in element.stress_mandel_pa() {
                near(*stress, 0., 0.1);
            }
        }
        let volume = 0.0001 * (1. + strain[0]) * (1. + strain[1]) * (1. + strain[2]);
        near(
            step.solid.geometry().element_volumes_m3().iter().sum(),
            volume,
            1e-13,
        );
        near(step.solid.wet_mass_kg(), 0.0575, 1e-13);
        for rho in step.solid.current_density_kg_m3() {
            near(*rho, 0.0575 / volume, 1e-5);
        }
        near(
            step.solid.tangent().unwrap().pencil().total_mass_kg,
            0.0575,
            1e-13,
        );
        near(
            step.solid
                .current_density_kg_m3()
                .iter()
                .zip(step.solid.geometry().element_volumes_m3())
                .map(|(rho, v)| rho * v)
                .sum(),
            step.moisture.wet_mass_kg(),
            1e-13,
        );
        assert_eq!(step.solid.time_s(), step.moisture.time_s());
        assert_eq!(original.time_s(), 0.);
        assert_eq!(original.moisture(), &[0.1]);
        let again = binding
            .advance(&original, 0.25, &faces(0.2), config(), cx)
            .unwrap();
        assert_eq!(again.solid.solution(), step.solid.solution());
    });
}

#[test]
fn constrained_body_has_analytical_compression_and_no_invented_free_dofs() {
    with_cx(|cx| {
        let original = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let binding = bind(
            &original,
            SolidStiffnessModel::Isotropic,
            axes(),
            (0..24).collect(),
            cx,
        );
        let step = binding
            .advance(&original, 0.25, &faces(0.2), config(), cx)
            .unwrap();
        let e = step.solid.elastic_states()[0]
            .property("young_modulus")
            .unwrap()
            .value_si();
        let mu = e / 2.5;
        let lambda = 0.4 * e;
        let eps = [0.0005, 0.002, 0.001];
        let trace = eps.iter().sum::<f64>();
        let expected = eps.map(|strain| -(lambda * trace + 2. * mu * strain));
        for state in step.solid.stress().elements() {
            for i in 0..3 {
                near(state.stress_mandel_pa()[i], expected[i], 1e-5);
            }
        }
        let energy = 0.5 * 0.0001 * (0..3).map(|i| -eps[i] * expected[i]).sum::<f64>();
        near(step.solid.stress().elastic_energy_j(), energy, 1e-9);
        assert!(step.solid.tangent().is_none());
        for u in step.solid.solution().displacement_m() {
            assert_eq!(u, &[0.; 3]);
        }
        for k in 0..3 {
            near(
                step.solid
                    .stress()
                    .nodal_internal_force_n()
                    .iter()
                    .map(|f| f[k])
                    .sum(),
                0.,
                1e-7,
            );
        }
    });
}

#[test]
fn data_only_direction_and_modulus_changes_reach_geometry_and_native_pencil() {
    with_cx(|cx| {
        let wet = slab(card([0.01; 2], 0.3), 0.1, vec![0.15], cx);
        let normal = bind(
            &wet,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        )
        .resolve(&wet, cx)
        .unwrap();
        let mut rotated = axes();
        rotated.principal_to_mesh = [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]];
        let rotated = bind(
            &wet,
            SolidStiffnessModel::Isotropic,
            rotated,
            free_supports(),
            cx,
        )
        .resolve(&wet, cx)
        .unwrap();
        near(normal.solution().displacement_m()[4][0], 0.00005, 1e-10);
        near(rotated.solution().displacement_m()[4][0], 0.0002, 1e-10);
        let changed = slab(
            card_with_swelling([0.01; 2], 0.3, [0.02, 0.04, 0.02]),
            0.1,
            vec![0.15],
            cx,
        );
        let other = bind(
            &changed,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        )
        .resolve(&changed, cx)
        .unwrap();
        near(
            other.solution().displacement_m()[4][0],
            2. * normal.solution().displacement_m()[4][0],
            1e-10,
        );
        let model = SolidStiffnessModel::OrientedOrthotropic {
            principal_to_mesh: axes().principal_to_mesh,
            orientation_identity: ContentHash([9; 32]),
        };
        let orthotropic = bind(&wet, model, axes(), free_supports(), cx)
            .resolve(&wet, cx)
            .unwrap();
        assert_ne!(
            normal.tangent().unwrap().elastic_stiffness().row(0).1,
            orthotropic.tangent().unwrap().elastic_stiffness().row(0).1
        );
        near(orthotropic.wet_mass_kg(), normal.wet_mass_kg(), 1e-14);
    });
}

#[test]
fn nonuniform_water_keeps_local_source_states_and_refusals_keep_old_inventory() {
    with_cx(|cx| {
        let field = slab(card([0.01; 2], 0.3), 0.1, vec![0.2, 0.1], cx);
        let binding = bind(
            &field,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        );
        let state = binding.resolve(&field, cx).unwrap();
        assert_eq!(
            state.expansion_paths()[1].free_principal_strains(),
            &[0.; 3]
        );
        assert!(state.expansion_paths()[0].free_principal_strains()[0] > 0.);
        assert_ne!(
            state.materials()[0].material_state_identity,
            state.materials()[6].material_state_identity
        );
        let restricted = slab(card([0.01; 2], 0.11), 0.1, vec![0.1], cx);
        let binding = bind(
            &restricted,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            cx,
        );
        assert!(
            binding
                .advance(&restricted, 0.25, &faces(0.2), config(), cx)
                .is_err()
        );
        assert_eq!(restricted.time_s(), 0.);
        assert_eq!(restricted.moisture(), &[0.1]);
        assert!(binding.resolve(&field, cx).is_err());
    });
}

#[test]
fn shared_one_dof_modal_pencil_changes_with_transported_water_and_material() {
    with_cx(|cx| {
        let original = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let fixed = (0..24).filter(|d| *d != 12).collect();
        let binding = bind(&original, SolidStiffnessModel::Isotropic, axes(), fixed, cx);
        let before = binding.resolve(&original, cx).unwrap();
        let after = binding
            .advance(&original, 0.25, &faces(0.2), config(), cx)
            .unwrap();
        // One generalized DOF: lambda=K00/M00 exactly, with native matrices.
        let lambda = |s: &MoistureSolidState| {
            let p = s.tangent().unwrap().pencil();
            assert_eq!(p.free_dofs.as_slice(), &[12]);
            p.stiffness.row(0).1[0] / p.mass.row(0).1[0]
        };
        let dry = lambda(&before);
        let wet = lambda(&after.solid);
        assert!(wet > 0. && wet < dry);
        assert_ne!(
            before.materials()[0].material_state_identity,
            after.solid.materials()[0].material_state_identity
        );
        assert!(after.solid.stress().elastic_energy_j() > 0.);
    });
}

#[test]
fn finite_vapor_uptake_reaches_compatible_shape_prestress_and_native_modes_atomically() {
    use super::atmosphere::{atmosphere, controls as atmosphere_controls, sorption};
    with_cx(|cx| {
        let gas = atmosphere(0.8, 0.2);
        let accepted_gas = gas.clone();
        let initial = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], cx);
        let fixed = || (0..24).filter(|d| *d != 12).collect();
        let binding = bind(
            &initial,
            SolidStiffnessModel::Isotropic,
            axes(),
            fixed(),
            cx,
        );
        let before = binding.resolve(&initial, cx).unwrap();
        let boundaries = [Some(sorption(0.25)), Some(sorption(0.25))];
        let (candidate, solid) = binding
            .advance_finite_atmosphere(
                &gas,
                &initial,
                0.25,
                &boundaries,
                config(),
                atmosphere_controls(),
                cx,
            )
            .unwrap();
        let water = candidate.step().water_to_solid_kg;
        assert!(water > 0.0);
        assert!(
            candidate.atmosphere().relative_humidity().unwrap() < gas.relative_humidity().unwrap()
        );
        near(
            candidate.atmosphere().gas().mass_kg() + solid.wet_mass_kg(),
            gas.gas().mass_kg() + before.wet_mass_kg(),
            2e-13,
        );
        near(solid.wet_mass_kg() - before.wet_mass_kg(), water, 1e-13);
        assert_eq!(solid.time_s(), candidate.moisture().time_s());
        let increment = candidate.moisture().moisture()[0] - 0.1;
        for (actual, beta) in solid.expansion_paths()[0]
            .free_principal_strains()
            .iter()
            .zip([0.01, 0.04, 0.02])
        {
            near(*actual, beta * increment, 1e-13);
        }
        assert!(solid.solution().displacement_m()[4][0] > 0.0);
        for ((reference, displacement), deformed) in binding
            .reference_nodes_m()
            .iter()
            .zip(solid.solution().displacement_m())
            .zip(solid.geometry().nodes_m())
        {
            for axis in 0..3 {
                near(deformed[axis], reference[axis] + displacement[axis], 1e-13);
            }
        }
        assert!(solid.stress().elastic_energy_j() > 0.0);
        let tangent = solid.tangent().unwrap();
        assert!(
            tangent
                .geometric_stiffness()
                .row(0)
                .1
                .iter()
                .any(|value| value.abs() > 1e-9)
        );
        near(
            tangent.pencil().total_mass_kg,
            candidate.moisture().wet_mass_kg(),
            1e-13,
        );
        // One native generalized DOF supplies its exact eigenvalue; this uses
        // the actual signed K_e+K_g and transported wet M, without a new solver.
        let eigenvalue = |state: &MoistureSolidState| {
            let pencil = state.tangent().unwrap().pencil();
            assert_eq!(pencil.free_dofs.as_slice(), &[12]);
            pencil.stiffness.row(0).1[0] / pencil.mass.row(0).1[0]
        };
        assert!(eigenvalue(&solid) > 0.0 && eigenvalue(&solid) < eigenvalue(&before));

        let restricted = slab(card([0.01; 2], 0.105), 0.1, vec![0.1], cx);
        let restricted_binding = bind(
            &restricted,
            SolidStiffnessModel::Isotropic,
            axes(),
            fixed(),
            cx,
        );
        let accepted_solid = restricted_binding.resolve(&restricted, cx).unwrap();
        let transport_only = gas
            .advance(
                &restricted,
                0.25,
                &boundaries,
                config(),
                atmosphere_controls(),
                cx,
            )
            .unwrap();
        assert!(transport_only.moisture().moisture()[0] > 0.105);
        assert!(
            restricted_binding
                .advance_finite_atmosphere(
                    &gas,
                    &restricted,
                    0.25,
                    &boundaries,
                    config(),
                    atmosphere_controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(gas, accepted_gas);
        assert_eq!(restricted.moisture(), &[0.1]);
        assert_eq!(restricted.time_s(), 0.0);
        let retained = restricted_binding.resolve(&restricted, cx).unwrap();
        assert_eq!(retained.solution(), accepted_solid.solution());
        assert_eq!(
            retained.tangent().unwrap().identity(),
            accepted_solid.tangent().unwrap().identity()
        );
    });
}

#[test]
fn cancelled_swell_candidate_never_debits_water_or_publishes_geometry() {
    // Cancellation at a native solid checkpoint retains the same public signal
    // as cancellation before diffusion or at the final publication boundary.
    assert!(matches!(
        MoistureError::from(fs_solid::linear3::TetElasticError::Cancelled),
        MoistureError::Cancelled
    ));
    assert!(matches!(
        MoistureError::from(fs_solid::linear3::source::SolidMaterialError::Solid(
            fs_solid::linear3::TetElasticError::Cancelled
        )),
        MoistureError::Cancelled
    ));
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let original = slab(card([0.01; 2], 0.3), 0.1, vec![0.1], &cx);
        let binding = bind(
            &original,
            SolidStiffnessModel::Isotropic,
            axes(),
            free_supports(),
            &cx,
        );
        gate.request();
        assert!(matches!(
            binding.advance(&original, 0.25, &faces(0.2), config(), &cx),
            Err(MoistureError::Cancelled)
        ));
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.time_s(), 0.);
    });
}
