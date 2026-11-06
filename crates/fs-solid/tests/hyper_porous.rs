//! Force-driven finite porous compaction through the existing tetrahedral solver.
//! The original mesh and caller-owned history survive loading/unloading.
//! Synthetic homogeneous coefficients do not establish measured foam accuracy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::porous::source::{FINITE_POROUS_LAW_ID, ResolvedFinitePorous};
use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousLimits, FinitePorousParameters,
};
use fs_qty::Dims;
use fs_solid::hyper3::plastic::porous::{FinitePorousTetError, FinitePorousTetSolid};
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetError,
};
use fs_solid::hyper3::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;

const MU: f64 = 80.0;
const K: f64 = 160.0;
const Y: f64 = 2.0;
const H: f64 = 8.0;
const A: f64 = 1.2;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.3e}"
    );
}

#[test]
fn native_force_equilibrium_retains_compaction_and_elastic_reload_without_a_new_solver() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 98,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let nodes = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let tets = [[0, 1, 2, 3]];
        // Remove six rigid motions; all three axial stretches and all three
        // independent shears remain unknowns in the native equilibrium solve.
        let prescribed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let material = FinitePorous::new(
            FinitePorousParameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: K,
                yield_stress_pa: Y,
                hardening_modulus_pa: H,
                pressure_sensitivity: A,
            },
            FinitePorousLimits {
                minimum_jacobian: 0.2,
                maximum_elastic_log_strain: 0.5,
                maximum_plastic_log_strain: 2.0,
                maximum_accumulated_plastic_strain: 3.0,
            },
        )
        .unwrap();
        let beta = A / (1.0 + A * A / 9.0).sqrt();
        let peak_stretch = 0.8_f64;
        let alpha = (-3.0 * K * beta * peak_stretch.ln() - Y) / (H + K * beta * beta);
        let plastic_j = (-beta * alpha).exp();
        let unloaded_stretch = (-beta * alpha / 3.0).exp();
        let reload_stretch = (peak_stretch * unloaded_stretch).sqrt();
        let peak_kirchhoff = -(Y + H * alpha) / beta;
        let mut accepted_state = material.initial_state();
        let mut accepted_u = vec![0.0; 12];
        let mut accumulated_heat = 0.0;

        for (step, stretch) in [peak_stretch, unloaded_stretch, reload_stretch]
            .into_iter()
            .enumerate()
        {
            let expected_kirchhoff = match step {
                0 => peak_kirchhoff,
                1 => 0.0,
                _ => 0.5 * peak_kirchhoff,
            };
            let nominal = expected_kirchhoff / stretch;
            // V_ref=1/6 and gradients of the three axial vertices are unit
            // Cartesian vectors. These are dead referential force resultants.
            let forces = [(3, nominal / 6.0), (7, nominal / 6.0), (11, nominal / 6.0)];
            let saved_state = accepted_state.clone();
            let frozen = [material.frozen(&accepted_state).unwrap()];
            let problem = HyperTetProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                material: HyperelasticLawField::PerElement(&frozen),
                prescribed_m: &prescribed,
                nodal_forces_n: &forces,
                body_force_n_m3: [0.0; 3],
                budget: TetAssemblyBudget::standard(),
                settings: HyperTetSettings::default(),
            };
            let equilibrium = problem.solve_from(&accepted_u, &cx).unwrap();
            assert!(equilibrium.steps[0].linear_iterations > 0);
            let candidate_u: Vec<_> = equilibrium
                .displacement_m
                .iter()
                .flatten()
                .copied()
                .collect();
            let fields = problem.element_responses(&candidate_u, &cx).unwrap();
            assert_eq!(fields.len(), 1);
            let field = &fields[0];
            let response = material
                .response(&field.deformation_gradient, &accepted_state)
                .unwrap();
            material
                .validate_committed(&field.deformation_gradient, &response.state, 2e-7)
                .unwrap();
            // Native Newton iterations and field observation are frozen-history
            // trials. The caller publishes only after every candidate is admitted.
            assert_eq!(accepted_state, saved_state);
            let expected_f = [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch];
            for i in 0..9 {
                near(field.deformation_gradient[i], expected_f[i], 3e-8);
                let diagonal = i % 4 == 0;
                near(
                    field.first_piola_pa[i],
                    if diagonal { nominal } else { 0.0 },
                    3e-7,
                );
                near(
                    field.cauchy_pa[i],
                    if diagonal {
                        expected_kirchhoff / stretch.powi(3)
                    } else {
                        0.0
                    },
                    5e-7,
                );
            }
            near(field.reference_volume_m3, 1.0 / 6.0, 1e-15);
            near(field.current_volume_m3, stretch.powi(3) / 6.0, 3e-8);
            near(response.state.alpha, alpha, 3e-8);
            near(response.plastic_jacobian, plastic_j, 3e-8);
            let storage = 0.5 * expected_kirchhoff.powi(2) / K + 0.5 * H * alpha * alpha;
            let heat = if step == 0 { Y * alpha } else { 0.0 };
            near(response.stored_energy_j_m3, storage, 4e-7);
            near(response.plastic_dissipation_j_m3, heat, 3e-7);
            near(equilibrium.strain_energy_j, (storage + heat) / 6.0, 4e-7);
            near(field.potential_density_j_m3, storage + heat, 4e-7);
            if step > 0 {
                assert_eq!(
                    response.state, accepted_state,
                    "unload/reload must remain elastic"
                );
                near(response.plastic_increment, 0.0, 0.0);
            }
            accumulated_heat += response.plastic_dissipation_j_m3 / 6.0;
            accepted_state = response.state;
            accepted_u = candidate_u;
        }
        assert!(accepted_state.alpha > 0.5);
        assert!(material.plastic_jacobian(&accepted_state).unwrap() < 0.55);
        near(accumulated_heat, Y * alpha / 6.0, 4e-7);
    });
}

fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 98,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate)
    })
}

fn sourced_material(pressure_sensitivity: f64) -> ResolvedFinitePorous {
    const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
    let source = fs_blake3::hash_bytes(b"authored retained porous compression fixture");
    let card = ConstitutiveModelCard {
        law: LawId(FINITE_POROUS_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("shear_modulus_pa", MU, PA),
            ("bulk_modulus_pa", K, PA),
            ("yield_stress_pa", Y, PA),
            ("hardening_modulus_pa", H, PA),
            ("pressure_sensitivity", pressure_sensitivity, Dims::NONE),
            ("minimum_jacobian", 0.2, Dims::NONE),
            ("maximum_elastic_log_strain", 0.5, Dims::NONE),
            ("maximum_plastic_log_strain", 2.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 3.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
        sources: vec![source],
        provenance: Provenance {
            source: "authored retained porous compression fixture".into(),
            license: "test fixture".into(),
            artifact: Some(source),
        },
    };
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    ResolvedFinitePorous::try_from_model_card_at(&card, &context).unwrap()
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-7,
        maximum_endpoint_work_defect_j: 10.0,
        yield_tolerance_pa: 2e-7,
    }
}

#[test]
fn retained_source_drives_compaction_unloading_density_and_atomic_retry() {
    with_cx(|cx, _| {
        let nodes = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let tets = [[0, 1, 2, 3]];
        let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let material = sourced_material(A);
        let mut solid = FinitePorousTetSolid::try_new(
            &nodes,
            &tets,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let beta = A / (1.0 + A * A / 9.0).sqrt();
        let peak = 0.8_f64;
        let alpha = (-3.0 * K * beta * peak.ln() - Y) / (H + K * beta * beta);
        let plastic_j = (-beta * alpha).exp();
        let release = (-beta * alpha / 3.0).exp();
        let peak_tau = -(Y + H * alpha) / beta;
        let mut previous_stretch = 1.0;
        let mut previous_storage = 0.0;
        let mut total_heat = 0.0;

        for (index, stretch) in [peak, release, (peak * release).sqrt()]
            .into_iter()
            .enumerate()
        {
            let tau = match index {
                0 => peak_tau,
                1 => 0.0,
                _ => 0.5 * peak_tau,
            };
            let nominal = tau / stretch;
            let forces = [(3, nominal / 6.0), (7, nominal / 6.0), (11, nominal / 6.0)];
            let accepted_states = solid.element_states().to_vec();
            let accepted_u = solid.displacement_m().to_vec();
            let accepted_f = solid.element_deformation().to_vec();
            let accepted_j = solid.element_jacobians().to_vec();
            let step = solid
                .advance(
                    &supports,
                    &forces,
                    [0.0; 3],
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert!(step.equilibrium.steps[0].linear_iterations > 0);
            assert_eq!(solid.element_states(), accepted_states);
            assert_eq!(solid.displacement_m(), accepted_u);
            assert_eq!(solid.element_deformation(), accepted_f);
            assert_eq!(solid.element_jacobians(), accepted_j);

            if index == 0 {
                assert!(step.energy.endpoint_work_defect_j > 0.1);
                let refused = solid
                    .advance(
                        &supports,
                        &forces,
                        [0.0; 3],
                        HyperTetSettings::default(),
                        FinitePlasticStepControls {
                            maximum_endpoint_work_defect_j: 0.0,
                            ..controls()
                        },
                        cx,
                    )
                    .unwrap_err();
                assert!(matches!(
                    refused,
                    FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate {
                        quantity: "porous signed endpoint integration defect",
                        ..
                    })
                ));
                // This gate runs after converged equilibrium and every return map.
                // A refused body must leave the original histories ready for retry.
                assert_eq!(solid.element_states(), accepted_states);
                assert_eq!(solid.displacement_m(), accepted_u);
                assert_eq!(solid.element_jacobians(), accepted_j);
                let retry = solid
                    .advance(
                        &supports,
                        &forces,
                        [0.0; 3],
                        HyperTetSettings::default(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                assert_eq!(retry.equilibrium, step.equilibrium);
                assert_eq!(retry.energy, step.energy);
                assert_eq!(retry.solid.element_states(), step.solid.element_states());
                assert_eq!(
                    retry.solid.element_plastic_jacobians(),
                    step.solid.element_plastic_jacobians()
                );
            } else {
                assert_eq!(
                    step.solid.element_states(),
                    accepted_states,
                    "unloading and subyield reloading must be elastic"
                );
            }

            let stored = (0.5 * tau * tau / K + 0.5 * H * alpha * alpha) / 6.0;
            let heat = if index == 0 { Y * alpha / 6.0 } else { 0.0 };
            let endpoint_work = 0.5 * nominal * (stretch - previous_stretch);
            near(step.energy.stored_before_j, previous_storage, 4e-7);
            near(step.energy.stored_after_j, stored, 4e-7);
            near(step.energy.plastic_dissipation_j, heat, 4e-7);
            near(step.energy.incremental_potential_j, stored + heat, 4e-7);
            near(step.equilibrium.strain_energy_j, stored + heat, 4e-7);
            near(step.energy.internal_endpoint_work_j, endpoint_work, 4e-7);
            near(step.energy.dead_load_endpoint_work_j, endpoint_work, 4e-7);
            near(step.energy.support_endpoint_work_j, 0.0, 1e-15);
            near(
                step.energy.endpoint_work_defect_j,
                endpoint_work - stored + previous_storage - heat,
                8e-7,
            );
            near(step.energy.equilibrium_work_residual_j, 0.0, 1e-7);
            near(step.solid.element_states()[0].alpha, alpha, 3e-8);
            near(step.solid.element_jacobians()[0], stretch.powi(3), 3e-8);
            near(step.solid.element_plastic_jacobians()[0], plastic_j, 3e-8);
            let density = step.solid.current_density_kg_m3(&[600.0], cx).unwrap()[0];
            near(density, 600.0 / stretch.powi(3), 6e-5);
            near(
                density * step.solid.element_jacobians()[0] / 6.0,
                600.0 / 6.0,
                1e-12,
            );
            for component in 0..9 {
                let diagonal = component % 4 == 0;
                near(
                    step.solid.element_deformation()[0][component],
                    if diagonal { stretch } else { 0.0 },
                    3e-8,
                );
                near(
                    step.solid.element_first_piola_pa()[0][component],
                    if diagonal { nominal } else { 0.0 },
                    3e-7,
                );
            }
            assert!(std::ptr::eq(step.solid.reference_nodes_m(), &nodes[..]));
            assert!(std::ptr::eq(step.solid.reference_tetrahedra(), &tets[..]));
            match step.solid.material() {
                FinitePlasticLawField::Uniform(source) => {
                    assert!(std::ptr::eq(source, &material));
                    assert_eq!(
                        source.context(),
                        &QueryPoint::new().with("T", 300.0).unwrap()
                    );
                    assert_eq!(source.model_card().law.0, FINITE_POROUS_LAW_ID);
                }
                FinitePlasticLawField::PerElement(_) => panic!("uniform binding changed"),
            }
            total_heat += step.energy.plastic_dissipation_j;
            previous_stretch = stretch;
            previous_storage = stored;
            solid = step.solid;
        }
        near(total_heat, Y * alpha / 6.0, 4e-7);
        assert!(solid.element_plastic_jacobians()[0] < 0.55);
        assert!(solid.current_density_kg_m3(&[], cx).is_err());
        for invalid_density in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(solid.current_density_kg_m3(&[invalid_density], cx).is_err());
        }
    });
}

#[test]
fn heterogeneous_retained_sources_preserve_pressure_sensitivity_and_cancellation() {
    with_cx(|cx, gate| {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let materials = [sourced_material(A), sourced_material(0.0)];
        assert!(matches!(
            FinitePorousTetSolid::try_new(
                &nodes,
                &tets,
                FinitePlasticLawField::PerElement(&materials[..1]),
                TetAssemblyBudget::standard(),
                cx,
            ),
            Err(FinitePorousTetError::Solid(FinitePlasticTetError::Native(
                HyperTetError::LawCountMismatch { .. }
            )))
        ));
        let original = FinitePorousTetSolid::try_new(
            &nodes,
            &tets,
            FinitePlasticLawField::PerElement(&materials),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let stretch = 0.9_f64;
        let mut prescribed = Vec::new();
        for (index, node) in nodes.iter().enumerate() {
            let origin = if index < 4 { [0.0; 3] } else { [2.0, 0.0, 0.0] };
            for axis in 0..3 {
                prescribed.push((
                    3 * index + axis,
                    (stretch - 1.0) * (node[axis] - origin[axis]),
                ));
            }
        }
        let step = original
            .advance(
                &prescribed,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let beta = A / (1.0 + A * A / 9.0).sqrt();
        let alpha = (-3.0 * K * beta * stretch.ln() - Y) / (H + K * beta * beta);
        near(step.solid.element_states()[0].alpha, alpha, 1e-12);
        assert_eq!(
            step.solid.element_states()[1],
            materials[1].law().initial_state()
        );
        near(
            step.solid.element_plastic_jacobians()[0],
            (-beta * alpha).exp(),
            1e-12,
        );
        near(step.solid.element_plastic_jacobians()[1], 1.0, 0.0);
        near(
            step.solid.element_first_piola_pa()[0][0],
            -(Y + H * alpha) / (beta * stretch),
            1e-10,
        );
        near(
            step.solid.element_first_piola_pa()[1][0],
            3.0 * K * stretch.ln() / stretch,
            1e-10,
        );
        near(step.energy.plastic_dissipation_j, Y * alpha / 6.0, 1e-12);
        for (index, density) in step
            .solid
            .current_density_kg_m3(&[500.0, 700.0], cx)
            .unwrap()
            .iter()
            .enumerate()
        {
            near(*density, [500.0, 700.0][index] / stretch.powi(3), 1e-10);
        }
        match step.solid.material() {
            FinitePlasticLawField::PerElement(source) => {
                assert!(std::ptr::eq(source, &materials[..]))
            }
            FinitePlasticLawField::Uniform(_) => panic!("heterogeneous bindings were averaged"),
        }
        let accepted = step.solid.element_states().to_vec();
        gate.request();
        assert!(matches!(
            step.solid.advance(
                &prescribed,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx
            ),
            Err(FinitePorousTetError::Solid(FinitePlasticTetError::Native(
                HyperTetError::Cancelled
            )))
        ));
        assert_eq!(step.solid.element_states(), accepted);
    });
}
