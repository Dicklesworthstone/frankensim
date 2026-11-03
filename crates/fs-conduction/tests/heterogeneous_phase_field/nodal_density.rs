//! Source-backed phase transport with one immutable graded reference inventory.
use super::*;

const PROFILE: [f64; 8] = [3.0, 15.0, 3.0, 3.0, 6.0, 18.0, 6.0, 6.0];
// Exact integrals of (3+12x) N_i and (6+12x) N_i over unit right tets.
const MASSES: [f64; 8] = [0.225, 0.325, 0.225, 0.225, 0.35, 0.45, 0.35, 0.35];

#[test]
fn original_sources_and_graded_masses_drive_independent_mixed_phase_heat_balance() {
    let mesh = mesh(false);
    let boundary = boundary(&mesh);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    let target = [20.0, 95.0, 60.0, 160.0, 1020.0, 1095.0, 1060.0, 1160.0];
    let mut q = [0.0; 8];
    for material in 0..2 {
        let first = 4 * material;
        let t: [f64; 4] = std::array::from_fn(|i| temperature(target[first + i], material));
        let k = average_k(&target[first..first + 4], material);
        let laplace = [
            3.0 * t[0] - t[1] - t[2] - t[3],
            t[1] - t[0],
            t[2] - t[0],
            t[3] - t[0],
        ];
        // Independent right-tet stiffness k*V*grad(N_i).grad(N_j).
        // Invert the source integral V*(I+11^T)/20 analytically.
        let load: [f64; 4] = std::array::from_fn(|i| {
            MASSES[first + i] * (target[first + i] - OLD[first + i]) / DT + k * laplace[i] / 6.0
        });
        let sum: f64 = load.iter().sum();
        for i in 0..4 {
            q[first + i] = 120.0 * (load[i] - sum / 5.0);
        }
    }
    let heat = ScalarField::Nodal(q.to_vec());
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        let storage =
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials_with_nodal_reference_density(
                cx,
                &mesh,
                &materials,
                &[0, 1],
                &PROFILE,
                budget(),
            )
            .unwrap();
        let masses = storage.reference_nodal_masses_kg().to_vec();
        for (got, expected) in masses.iter().zip(MASSES) {
            near(*got, expected, 1e-15);
        }
        assert_eq!(
            storage.nodal_reference_density_kg_m3(),
            Some(PROFILE.as_slice())
        );
        for material in 0..2 {
            assert!(storage.phase_transport_material(material).is_none());
            assert!(storage.reference_material(material).is_none());
            assert!(std::ptr::eq(
                storage.phase_conductivity_source(material).unwrap(),
                &sources[material]
            ));
            assert!(std::ptr::eq(
                storage.phase_curve_at(4 * material).unwrap(),
                &curves[material]
            ));
        }
        assert!(storage.phase_conductivity_source(2).is_none());
        let result = storage
            .advance(
                cx,
                storage.phase_transport_problem(&boundary, &heat).unwrap(),
                None,
                &OLD,
                DT,
                config(),
            )
            .unwrap();
        for i in 0..8 {
            near(result.specific_enthalpy_j_kg[i], target[i], 2e-8);
            near(result.temperature[i], temperature(target[i], i / 4), 2e-8);
            near(
                result.liquid_mass_fraction[i],
                [0.2, 0.95, 0.6, 1.0][i % 4],
                2e-10,
            );
        }
        let expected_energy: f64 = (0..8).map(|i| MASSES[i] * (target[i] - OLD[i])).sum();
        near(result.stored_energy_change_j, expected_energy, 2e-9);
        near(result.energy_residual_j, 0.0, 1e-10);
        let flux = storage
            .phase_transport_heat_flux(cx, &result.specific_enthalpy_j_kg)
            .unwrap();
        for material in 0..2 {
            let h = &target[4 * material..4 * material + 4];
            for axis in 0..3 {
                near(
                    flux[material][axis],
                    -average_k(h, material)
                        * (temperature(h[axis + 1], material) - temperature(h[0], material)),
                    2e-7,
                );
            }
        }
        // Equal volumetric deposition on the latent plateau: exact local
        // increment 0.25/m_i, so the denser corner melts more slowly.
        let latent = [20.0, 20.0, 20.0, 20.0, 1020.0, 1020.0, 1020.0, 1020.0];
        let heat = ScalarField::Uniform(12.0);
        let result = storage
            .advance(
                cx,
                storage.phase_transport_problem(&boundary, &heat).unwrap(),
                None,
                &latent,
                0.5,
                config(),
            )
            .unwrap();
        for i in 0..8 {
            near(
                result.specific_enthalpy_j_kg[i],
                latent[i] + 0.25 / MASSES[i],
                2e-9,
            );
        }
        assert!(result.liquid_mass_fraction[1] < result.liquid_mass_fraction[0]);
        assert_eq!(storage.reference_nodal_masses_kg(), masses);
    });
}

#[test]
fn profiled_source_contact_and_density_adjoint_share_the_actual_heat_inventory() {
    let mesh = mesh(true);
    let boundary = boundary(&mesh);
    let interfaces = contact(&mesh, &boundary);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    let zero = ScalarField::Uniform(0.0);
    let seed = [0.7, -0.4, 0.2, 1.1, -0.3, 0.6, -0.2, 0.4];
    let direction = [0.4, -0.7, 0.2, 0.5, -0.3, 0.8, -0.6, 0.1];
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        let storage =
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials_with_nodal_reference_density(
                cx,
                &mesh,
                &materials,
                &[0, 1],
                &PROFILE,
                budget(),
            )
            .unwrap();
        let problem = storage.phase_transport_problem(&boundary, &zero).unwrap();
        let response = storage
            .linearize_step(cx, problem, Some(&interfaces), &OLD, DT, config())
            .unwrap();
        let h = &response.primal().specific_enthalpy_j_kg;
        // Independent original k(h) moments and consistent interface matrix,
        // replacing only the analytic reference-storage term in that oracle.
        for (i, residual) in contact_residual(h, &OLD).into_iter().enumerate() {
            near(
                residual + (MASSES[i] - DENSITY[i / 4] / 24.0) * (h[i] - OLD[i]),
                0.0,
                2e-10,
            );
        }
        let first_energy: f64 = (0..4).map(|i| MASSES[i] * (h[i] - OLD[i])).sum();
        near(
            first_energy,
            -DT * response.primal().contact_fluxes[0].heat_rate_a_to_b_w,
            1e-10,
        );
        near(response.primal().stored_energy_change_j, 0.0, 1e-10);
        let gradient = response
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-12,
                    max_iterations: 128,
                    restart: 8,
                },
            )
            .unwrap();
        let nodal = gradient.nodal_reference_density.as_ref().unwrap();
        let action = response
            .apply_nodal_reference_density_jacobian(cx, &direction)
            .unwrap();
        for i in 0..8 {
            let first = 4 * (i / 4);
            let sum: f64 = direction[first..first + 4].iter().sum();
            near(
                action[i],
                (sum + direction[i]) * (h[i] - OLD[i]) / 120.0,
                1e-12,
            );
        }
        near(
            dot(nodal, &direction),
            -dot(&gradient.transport.adjoint, &action),
            1e-11,
        );
        let objective = |delta| {
            let profile: Vec<_> = PROFILE
                .iter()
                .zip(direction)
                .map(|(rho, d)| rho + delta * d)
                .collect();
            let perturbed = HeterogeneousEnthalpyBackwardEuler::from_phase_materials_with_nodal_reference_density(
                cx, &mesh, &materials, &[0, 1], &profile, budget(),
            ).unwrap();
            let next = perturbed
                .advance(
                    cx,
                    perturbed.phase_transport_problem(&boundary, &zero).unwrap(),
                    Some(&interfaces),
                    &OLD,
                    DT,
                    config(),
                )
                .unwrap();
            (0..8)
                .map(|i| seed[i] * (next.specific_enthalpy_j_kg[i] - OLD[i]))
                .sum::<f64>()
        };
        let e = 0.001;
        let finite_difference = (-objective(2.0 * e) + 8.0 * objective(e) - 8.0 * objective(-e)
            + objective(-2.0 * e))
            / (12.0 * e);
        near(dot(nodal, &direction), finite_difference, 2e-7);
    });
}

#[test]
fn constant_profile_preserves_source_transport_bits_without_exposing_scalar_inventory() {
    let mesh = mesh(false);
    let boundary = boundary(&mesh);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    let zero = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        let scalar = HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
            cx,
            &mesh,
            &materials,
            &[0, 1],
            budget(),
        )
        .unwrap();
        let profiled =
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials_with_nodal_reference_density(
                cx,
                &mesh,
                &materials,
                &[0, 1],
                &[2.0, 2.0, 2.0, 2.0, 6.0, 6.0, 6.0, 6.0],
                budget(),
            )
            .unwrap();
        assert_eq!(
            profiled.reference_nodal_masses_kg(),
            scalar.reference_nodal_masses_kg()
        );
        let advance = |storage: &HeterogeneousEnthalpyBackwardEuler<'_, '_>| {
            storage
                .advance(
                    cx,
                    storage.phase_transport_problem(&boundary, &zero).unwrap(),
                    None,
                    &OLD,
                    DT,
                    config(),
                )
                .unwrap()
        };
        let expected = advance(&scalar);
        let actual = advance(&profiled);
        assert_eq!(
            actual.specific_enthalpy_j_kg,
            expected.specific_enthalpy_j_kg
        );
        assert_eq!(
            actual.stored_energy_change_j.to_bits(),
            expected.stored_energy_change_j.to_bits()
        );
        for material in 0..2 {
            assert!(scalar.phase_transport_material(material).is_some());
            assert!(scalar.reference_material(material).is_some());
            assert!(profiled.phase_transport_material(material).is_none());
            assert!(profiled.reference_material(material).is_none());
            assert!(std::ptr::eq(
                profiled.phase_conductivity_source(material).unwrap(),
                scalar.phase_conductivity_source(material).unwrap()
            ));
        }
        assert!(profiled.with_nodal_reference_density(cx, &PROFILE).is_err());
        assert!(scalar.with_nodal_reference_density(cx, &PROFILE).is_err());
    });
}

#[test]
fn profiled_source_preparation_rejects_invalid_inventory_and_never_assembles_superseded_mass() {
    let mesh = mesh(false);
    let cards = cards();
    let curves = charts(&cards);
    let sources = sources(&cards, &curves);
    with_cx(|cx| {
        let materials = materials(cx, &sources, DENSITY);
        let prepare = |profile: &[f64], ids: &[usize], budget| {
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials_with_nodal_reference_density(
                cx, &mesh, &materials, ids, profile, budget,
            )
        };
        assert!(prepare(&PROFILE[..7], &[0, 1], budget()).is_err());
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut bad = PROFILE;
            bad[6] = value;
            assert!(prepare(&bad, &[0, 1], budget()).is_err());
        }
        assert!(prepare(&PROFILE, &[0, 2], budget()).is_err());
        assert!(
            prepare(
                &PROFILE,
                &[0, 1],
                EnthalpyBudget {
                    max_vertices: 7,
                    max_elements: 2
                }
            )
            .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(HeterogeneousEnthalpyBackwardEuler::from_phase_materials_with_nodal_reference_density(
                cancelled, &mesh, &materials, &[0, 1], &PROFILE, budget(),
            ).is_err());
        });
        let storage = prepare(&PROFILE, &[0, 1], budget()).unwrap();
        let boundary = boundary(&mesh);
        let zero = ScalarField::Uniform(0.0);
        let mut bad_h = OLD;
        bad_h[4] = 800.0;
        assert!(
            storage
                .advance(
                    cx,
                    storage.phase_transport_problem(&boundary, &zero).unwrap(),
                    None,
                    &bad_h,
                    DT,
                    config()
                )
                .is_err()
        );
        let mut problem = storage.phase_transport_problem(&boundary, &zero).unwrap();
        let unrelated = ConductivityModel::isotropic_declared(1.0).unwrap();
        problem.material = &unrelated;
        assert!(
            storage
                .advance(cx, problem, None, &OLD, DT, config())
                .is_err()
        );

        // A superseded scalar inventory would overflow on this larger tet.
        // The authored profile is valid and must be assembled directly.
        let large = ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
            vec![[0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.]],
        )
        .unwrap();
        let huge = self::materials(cx, &sources, [f64::MAX; 2]);
        assert!(
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials(
                cx,
                &large,
                &huge,
                &[0],
                budget()
            )
            .is_err()
        );
        let storage =
            HeterogeneousEnthalpyBackwardEuler::from_phase_materials_with_nodal_reference_density(
                cx,
                &large,
                &huge,
                &[0],
                &PROFILE[..4],
                budget(),
            )
            .unwrap();
        for (mass, unit) in storage.reference_nodal_masses_kg().iter().zip(MASSES) {
            near(*mass, 8.0 * unit, 2e-14);
        }
    });
}
