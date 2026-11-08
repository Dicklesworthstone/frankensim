//! Table-driven forming reaches the native P1/Q1 equilibrium and layer histories.
//! Authored coefficients test the coupling, not calibration of a particular metal.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::{ReductionConfig, UniaxialTabulated};
use fs_material::plastic::section::{
    LayeredPlaneStressSection, PlaneStressLayer, SectionBudget, ThicknessQuadrature,
};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening, TabulatedJ2};
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, FiberLaw, FiberState, Section};
use fs_solid::{
    ForceBasedElement, Formulation, Frame2d, Geometry, IncrementSettings, Mesh2,
    PlaneAssemblyBudget, PlaneKind, SectionModel, SmallStrainLawField, SmallStrainPlaneProblem,
    SmallStrainTetError, SmallStrainTetProblem, TetAssemblyBudget,
};

fn metal() -> TabulatedJ2 {
    TabulatedJ2::new(
        J2Plasticity::new(IsotropicElastic::new(1000.0, 0.3, 0.05).unwrap(), 1.0, 0.0).unwrap(),
        TabulatedHardening::new(
            [(0.0, 1.0), (0.002, 1.2), (0.006, 3.2), (0.02, 4.6)]
                .map(|(plastic_strain, yield_stress_pa)| HardeningPoint {
                    plastic_strain,
                    yield_stress_pa,
                })
                .to_vec(),
        )
        .unwrap(),
    )
    .unwrap()
}
fn controls() -> ReductionConfig {
    ReductionConfig {
        max_iterations: 48,
        stress_tolerance_pa: 1e-11,
        strain_limit: 0.05,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 92,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn supports(mesh: &Mesh2, axial_strain: Option<f64>) -> Vec<(usize, f64)> {
    let mut prescribed = Vec::new();
    for (i, x) in mesh.nodes.iter().enumerate() {
        for j in 0..2 {
            if x[j] == 0.0 {
                prescribed.push((2 * i + j, 0.0));
            }
        }
        if x[0] == 1.0
            && let Some(strain) = axial_strain
        {
            prescribed.push((2 * i, strain));
        }
    }
    prescribed
}
fn problem<'a>(
    mesh: &'a Mesh2,
    law: &'a TabulatedJ2,
    fixed: &'a [(usize, f64)],
) -> SmallStrainPlaneProblem<'a, TabulatedJ2> {
    SmallStrainPlaneProblem {
        mesh,
        laws: SmallStrainLawField::Uniform(law),
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        thickness_m: 0.2,
        reduction: controls(),
        rotation_limit_rad: 0.02,
        prescribed_m: fixed,
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget: PlaneAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn p1_and_q1_form_and_release_at_the_measured_curve_without_losing_thickness_history() {
    let law = metal();
    // The middle segment is Y(alpha)=0.2+500 alpha. At alpha=.004,
    // sigma=2.2 and epsilon=alpha+sigma/E=.0062, independently of the mesh.
    for mesh in [
        Mesh2::triangles(1.0, 1.0, 2, 2),
        Mesh2::quads(1.0, 1.0, 2, 2),
    ] {
        let fixed = supports(&mesh, Some(0.0062));
        let loaded = with_cx(|cx| problem(&mesh, &law, &fixed).solve(cx)).unwrap();
        for (state, stress) in loaded
            .quadrature_states
            .iter()
            .zip(&loaded.quadrature_stress_pa)
        {
            near(state.history.alpha, 0.004, 2e-10);
            near(stress[0], 2.2, 2e-8);
            near(state.strain[0], 0.0062, 2e-10);
            for j in 1..3 {
                near(state.history.plastic_strain[j], -0.002, 2e-10);
                near(state.strain[j], -0.00266, 2e-10);
            }
            for value in &stress[1..] {
                near(*value, 0.0, 2e-8);
            }
        }
        let reaction: f64 = fixed
            .iter()
            .zip(&loaded.increments[0].reactions_n)
            .filter(|((dof, _), _)| dof % 2 == 0 && mesh.nodes[dof / 2][0] == 1.0)
            .map(|(_, force)| force)
            .sum();
        near(reaction, 0.44, 3e-9);
        let released = supports(&mesh, None);
        let unloaded = with_cx(|cx| {
            problem(&mesh, &law, &released).solve_from(
                loaded.displacement_m.as_flattened(),
                &loaded.quadrature_states,
                cx,
            )
        })
        .unwrap();
        for (x, u) in mesh.nodes.iter().zip(&unloaded.displacement_m) {
            near(u[0], 0.004 * x[0], 2e-10);
            near(u[1], -0.002 * x[1], 2e-10);
        }
        for (state, stress) in unloaded
            .quadrature_states
            .iter()
            .zip(&unloaded.quadrature_stress_pa)
        {
            near(state.history.alpha, 0.004, 2e-10);
            near(state.strain[2], -0.002, 2e-10);
            for value in stress {
                near(*value, 0.0, 2e-8);
            }
        }
        let saved = unloaded.clone();
        let unsupported = supports(&mesh, Some(0.035));
        assert!(
            with_cx(|cx| problem(&mesh, &law, &unsupported).solve_from(
                unloaded.displacement_m.as_flattened(),
                &unloaded.quadrature_states,
                cx,
            ))
            .is_err()
        );
        assert_eq!(
            unloaded, saved,
            "refused extrapolation must retain the accepted state"
        );
    }
}

#[test]
fn tetrahedral_force_control_recovers_calibrated_equilibrium_after_domain_overshoot() {
    // A nearly flat measured segment followed by strong hardening makes the
    // full Newton proposal overshoot the final measurement. The physical
    // equilibrium is in the interior: alpha=.004, sigma=3.002, eps=.007002.
    let law = TabulatedJ2::new(
        J2Plasticity::new(IsotropicElastic::new(1000.0, 0.3, 0.05).unwrap(), 1.0, 0.0).unwrap(),
        TabulatedHardening::new(
            [(0.0, 1.0), (0.002, 1.002), (0.006, 5.002), (0.008, 7.002)]
                .map(|(plastic_strain, yield_stress_pa)| HardeningPoint {
                    plastic_strain,
                    yield_stress_pa,
                })
                .to_vec(),
        )
        .unwrap(),
    )
    .unwrap();
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tets = [[0, 1, 2, 3]];
    // Three free diagonal strains, with zero transverse traction. The volume
    // is 1/6, hence the nodal axial force is sigma/6 for a unit axial gradient.
    let prescribed = [0, 1, 2, 4, 5, 6, 8, 9, 10].map(|dof| (dof, 0.0));
    let forces = [(3, 3.002 / 6.0)];
    let mut problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(&law),
        prescribed_m: &prescribed,
        nodal_forces_n: &forces,
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings {
            max_backtracks: 0,
            ..IncrementSettings::default()
        },
    };
    assert!(matches!(
        with_cx(|cx| problem.solve(cx)),
        Err(SmallStrainTetError::NewtonStalled { .. })
    ));
    problem.settings = IncrementSettings::default();
    let loaded = with_cx(|cx| problem.solve(cx)).unwrap();
    near(loaded.element_states[0].alpha, 0.004, 2e-10);
    near(loaded.element_stress[0][0], 3.002, 2e-8);
    near(loaded.element_strain[0][0], 0.007002, 2e-10);
    for j in 1..3 {
        near(loaded.element_states[0].plastic_strain[j], -0.002, 2e-10);
        near(loaded.element_strain[0][j], -0.0029006, 2e-10);
    }
    for value in &loaded.element_stress[0][1..] {
        near(*value, 0.0, 2e-8);
    }

    // Rejected trial steps must never accumulate plastic strain. Elastic
    // release from the actual accepted state leaves exactly the permanent set.
    problem.nodal_forces_n = &[];
    let released = with_cx(|cx| {
        problem.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.element_states,
            cx,
        )
    })
    .unwrap();
    near(released.element_states[0].alpha, 0.004, 2e-10);
    near(released.element_strain[0][0], 0.004, 2e-10);
    for j in 1..3 {
        near(released.element_strain[0][j], -0.002, 2e-10);
    }
    for value in &released.element_stress[0] {
        near(*value, 0.0, 2e-8);
    }

    // The recovery does not extrapolate a material table or mutate history
    // when no admissible equilibrium exists beyond its final measured yield.
    let unsupported = [(3, 8.0 / 6.0)];
    problem.nodal_forces_n = &unsupported;
    let saved = released.clone();
    assert!(
        with_cx(|cx| problem.solve_from(
            released.displacement_m.as_flattened(),
            &released.element_states,
            cx,
        ))
        .is_err()
    );
    assert_eq!(released, saved);
}

#[test]
fn layered_bending_crosses_table_knots_and_springs_back_with_native_energy() {
    let section = LayeredPlaneStressSection::new(
        vec![PlaneStressLayer {
            bottom_z_m: -0.1,
            top_z_m: 0.1,
            subdivisions: 6,
            law: metal(),
        }],
        ThicknessQuadrature::Gauss3,
        controls(),
        SectionBudget {
            maximum_layers: 1,
            maximum_quadrature_points: 32,
        },
    )
    .unwrap();
    let virgin = section.initial_state();
    let peak = [0.0, 0.0, 0.0, 0.06, 0.0, 0.0];
    let loaded = section.evaluate(peak, &virgin).unwrap();
    assert!(loaded.state.points.iter().any(|p| p.history.alpha == 0.0));
    assert!(loaded.state.points.iter().any(|p| p.history.alpha > 0.002));
    let gates = EnergyTolerances {
        energy_j_m3: 1e-9,
        yield_pa: 1e-8,
    };
    let loading_energy = section
        .audit_energy(&virgin, &loaded.state, gates, 1e-10)
        .unwrap();
    let heat: f64 = section
        .points()
        .iter()
        .zip(&loaded.state.points)
        .map(|(point, state)| point.weight_m * state.history.alpha)
        .sum();
    near(loading_energy.physical_dissipation_j_m2, heat, 1e-11);
    assert!(loading_energy.final_state.hardening_j_m2 > 0.0);
    assert!(loading_energy.numerical_dissipation_j_m2 > 0.0);
    // Independent elastic plate compliance determines the released curvature.
    let ei = 1000.0 * 0.2_f64.powi(3) / 12.0;
    let mut residual = peak;
    residual[3] -= (loaded.resultants[3] - 0.3 * loaded.resultants[4]) / ei;
    residual[4] -= (loaded.resultants[4] - 0.3 * loaded.resultants[3]) / ei;
    residual[5] -= loaded.resultants[5] / (ei / 2.6);
    assert!(residual[3] > 0.0 && residual[3] < peak[3]);
    let unloaded = section.evaluate(residual, &loaded.state).unwrap();
    for resultant in unloaded.resultants {
        near(resultant, 0.0, 2e-11);
    }
    for (before, after) in loaded.state.points.iter().zip(&unloaded.state.points) {
        assert_eq!(before.history, after.history);
    }
    let unloading_energy = section
        .audit_energy(&loaded.state, &unloaded.state, gates, 1e-10)
        .unwrap();
    near(unloading_energy.physical_dissipation_j_m2, 0.0, 1e-13);
    near(unloading_energy.balance_residual_j_m2, 0.0, 1e-10);
    assert_eq!(virgin, section.initial_state());
}

fn fiber_section() -> Section {
    let axial = UniaxialTabulated::new(metal()).unwrap();
    Section::new(
        [-0.1, 0.1]
            .map(|y| Fiber {
                y,
                area: 1.0,
                law: axial.clone().into(),
            })
            .to_vec(),
    )
}

fn fiber_histories(section: &Section) -> Vec<fs_material::plastic::J2State> {
    section
        .states
        .iter()
        .map(|state| {
            let FiberState::Tabulated(state) = state else {
                panic!("expected tabulated history")
            };
            state.tensor_history()
        })
        .collect()
}

#[test]
fn tabulated_fibers_drive_beam_and_frame_permanent_deflection() {
    let mut beam = ForceBasedElement::new(1.0, &fiber_section);
    let force = 0.44;
    let (mut loaded_reference, mut residual_reference) = (0.0, 0.0);
    // Independent force equilibrium fixes sigma=M/(2 A y). Invert each
    // measured stress segment, then integrate curvature by virtual work.
    for (&xi, &weight) in beam.xi.iter().zip(&beam.wi) {
        let arm = 1.0 - xi;
        let sigma = force * arm / 0.2;
        let plastic = if sigma <= 1.0 {
            0.0
        } else if sigma <= 1.2 {
            (sigma - 1.0) / 100.0
        } else {
            0.002 + (sigma - 1.2) / 500.0
        };
        loaded_reference += weight * arm * (plastic + sigma / 1000.0) / 0.1;
        residual_reference += weight * arm * plastic / 0.1;
    }
    near(
        beam.tip_deflection_under_shear(force).unwrap(),
        loaded_reference,
        2e-8,
    );
    let loaded: Vec<_> = beam.sections.iter().map(fiber_histories).collect();
    near(loaded[0][0].alpha, 0.004, 2e-8);
    near(
        beam.tip_deflection_under_shear(0.0).unwrap(),
        residual_reference,
        2e-8,
    );
    assert!(residual_reference > 0.0);
    for (section, history) in beam.sections.iter().zip(loaded) {
        assert_eq!(fiber_histories(section), history);
    }

    let mut frame = Frame2d::new(Geometry::Linear);
    let base = frame.add_node(0.0, 0.0);
    let tip = frame.add_node(1.0, 0.0);
    frame.fix(base, [true; 3]);
    frame
        .add_element(base, tip, 5, 0.0, &|| SectionModel::Fiber(fiber_section()))
        .unwrap();
    let mut loads = vec![0.0; frame.ndof()];
    loads[3 * tip + 1] = force;
    frame.static_load(&loads, 8).unwrap();
    near(frame.displacements()[3 * tip + 1], loaded_reference, 2e-7);
    frame.static_load(&vec![0.0; frame.ndof()], 8).unwrap();
    near(frame.displacements()[3 * tip + 1], residual_reference, 2e-7);
}

#[test]
fn a_later_fiber_refusing_table_extrapolation_cannot_partially_commit_the_section() {
    let long = UniaxialTabulated::new(metal()).unwrap();
    let short = UniaxialTabulated::new(
        TabulatedJ2::new(
            metal().backbone().clone(),
            TabulatedHardening::new(vec![
                HardeningPoint {
                    plastic_strain: 0.0,
                    yield_stress_pa: 1.0,
                },
                HardeningPoint {
                    plastic_strain: 0.002,
                    yield_stress_pa: 1.2,
                },
            ])
            .unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let material = FiberMaterial::authored(long.clone().into(), 7850.0).unwrap();
    assert_eq!(material.density_kg_m3(), 7850.0);
    assert!(FiberMaterial::authored(long.clone().into(), f64::NAN).is_err());
    let mut section = Section::new(vec![
        Fiber {
            y: 0.0,
            area: 1.0,
            law: long.into(),
        },
        Fiber {
            y: 0.0,
            area: 1.0,
            law: short.into(),
        },
    ]);
    section.try_commit(0.002, 0.0).unwrap();
    let accepted = fiber_histories(&section);
    let (FiberLaw::Tabulated(law), FiberState::Tabulated(state)) =
        (&section.fibers[0].law, &section.states[0])
    else {
        unreachable!()
    };
    assert!(
        law.response(0.0062, state)
            .unwrap()
            .state
            .equivalent_plastic_strain()
            > state.equivalent_plastic_strain()
    );
    let error = section.try_commit(0.0062, 0.0).unwrap_err();
    assert!(error.to_string().contains("tabulated"));
    assert_eq!(fiber_histories(&section), accepted);
}
