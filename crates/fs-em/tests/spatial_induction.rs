//! Independent A-phi induction oracles on actual 3D Whitney fields.
//! Synthetic materials qualify equations/discretization, not empirical devices.
#[path = "spatial_induction/nonlinear.rs"]
mod nonlinear;
use fs_em::spatial::{
    SpatialMagneticError, SpatialMagneticMaterial, SpatialRecoilMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        induction::{SpatialInduction, SpatialInductionControls, SpatialInductionLimits},
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY, ResolvedConductorState,
    },
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_solver::CheckedCgConfig;
use std::collections::BTreeSet;

fn source(
    mu: f64,
    sigma: Option<f64>,
    br: Option<f64>,
) -> (SpatialMagneticMaterial, Option<ResolvedConductorState>) {
    let mut values = vec![(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        mu,
    )];
    if let Some(sigma) = sigma {
        values.push((
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            sigma,
        ));
    }
    if let Some(br) = br {
        values.extend([
            (RECOIL_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu),
            (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, br),
            (
                RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                -10.0,
            ),
            (
                RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                10.0,
            ),
            (
                RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
                FIELD_STRENGTH_DIMS,
                10.0,
            ),
        ]);
    }
    let mut claims = ClaimSet::new();
    let mut requirements = vec![];
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured 3D induction material".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let bundle = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor =
        sigma.map(|_| ResolvedConductorState::try_from_material(bundle.clone()).unwrap());
    let magnetic = if br.is_some() {
        SpatialMagneticMaterial::Recoil(
            SpatialRecoilMaterial::new(
                ResolvedRecoilMaterial::try_from_material(bundle).unwrap(),
                [0.0, 0.0, 1.0],
            )
            .unwrap(),
        )
    } else {
        SpatialMagneticMaterial::Linear(bundle)
    };
    (magnetic, conductor)
}
fn grid(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let mut positions = vec![];
    for k in 0..=n {
        for j in 0..=n {
            for i in 0..=n {
                positions.push([
                    i as f64 / n as f64,
                    j as f64 / n as f64,
                    k as f64 / n as f64,
                ]);
            }
        }
    }
    let id = |i, j, k| ((k * (n + 1) + j) * (n + 1) + i) as u32;
    let mut tets = vec![];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let (a, b, c, d) = (
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i, j + 1, k),
                    id(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i, j + 1, k + 1),
                    id(i + 1, j + 1, k + 1),
                );
                tets.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    (positions, tets)
}
fn assembly() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 4000,
        max_dofs: 20000,
        max_triplets: 200000,
    }
}
fn preparation() -> SpatialInductionLimits {
    SpatialInductionLimits {
        maximum_unknowns: 20000,
        maximum_conducting_cells: 4000,
        maximum_assembly_terms: 1_000_000,
    }
}
fn static_controls() -> VectorMagneticControls {
    let solve = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 3000,
        max_corrections: 2,
    };
    VectorMagneticControls {
        field: solve,
        gauge: solve,
        maximum_gauge_product_terms: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 2e-8,
    }
}
fn controls() -> SpatialInductionControls {
    SpatialInductionControls {
        relative_tolerance: 2e-11,
        absolute_residual_a: 1e-11,
        restart: 64,
        maximum_cycles: 100,
        maximum_krylov_scalars: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 2e-8,
        maximum_cell_magnetic_peclet: 100.0,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a:.16e} != {b:.16e}, error {:.6e}",
        (a - b).abs()
    );
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn trace(
    problem: &SpatialVectorMagnetostatics,
    b: [f64; 3],
    gradient: bool,
) -> VectorMagneticBoundary {
    let psi = |v: usize| {
        if gradient {
            0.13 * ((v * 17) % 11) as f64
        } else {
            0.0
        }
    };
    let a: Vec<_> = problem
        .mesh()
        .edges
        .iter()
        .map(|&[u, v]| {
            let p = problem.positions_m()[u as usize];
            let q = problem.positions_m()[v as usize];
            let mid = std::array::from_fn(|j| 0.5 * (p[j] + q[j]));
            let vector = cross(b, mid).map(|v| 0.5 * v);
            (0..3).map(|j| vector[j] * (q[j] - p[j])).sum::<f64>() + psi(v as usize)
                - psi(u as usize)
        })
        .collect();
    VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: problem.face_edges(face).unwrap().map(|e| a[e]),
            })
            .collect(),
    }
}
fn cube(
    n: usize,
    br: Option<f64>,
) -> (
    SpatialVectorMagnetostatics,
    Vec<Option<ResolvedConductorState>>,
) {
    let (positions, tets) = grid(n);
    let (m, c) = source(2.0, Some(5.0), br);
    let count = tets.len();
    (
        SpatialVectorMagnetostatics::new(positions, tets, vec![m; count], assembly(), &mut || {
            false
        })
        .unwrap(),
        vec![c; count],
    )
}
fn diagonal(problem: &SpatialVectorMagnetostatics) -> usize {
    problem.mesh().edges.binary_search(&[0, 7]).unwrap()
}

#[test]
fn uniform_translation_is_cancelled_by_real_floating_charge_potential_and_is_gauge_invariant() {
    let (problem, sources) = cube(2, None);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.7, -0.3, 0.2]; problem.mesh().vertex_count];
    let b = [0.2, -0.1, 0.8];
    let expected_gradient = cross(velocity[0], b);
    let mut outcomes = vec![];
    for gauge in [false, true] {
        let initial = problem
            .solve(
                &zero,
                &trace(&problem, b, gauge),
                static_controls(),
                &mut || false,
            )
            .unwrap();
        let state =
            SpatialInduction::new(&initial, sources.clone(), preparation(), 2e-9, &mut || {
                false
            })
            .unwrap();
        let candidate = state
            .step_candidate(&zero, &velocity, 0.1, controls(), &mut || false)
            .unwrap();
        let step = candidate.step();
        near(step.heat_j, 0.0, 1e-17);
        near(step.mechanical_work_j, 0.0, 2e-10);
        near(step.stored_before_j, step.stored_after_j, 1e-10);
        for field in &step.endpoint_flux_density_t {
            for j in 0..3 {
                near(field[j], b[j], 2e-9);
            }
        }
        for cell in &step.cells {
            let vertices = problem.mesh().tets[cell.cell];
            for k in 0..4 {
                for j in 0..3 {
                    near(cell.current_density_a_m2[k][j], 0.0, 2e-8);
                }
                let p = problem.positions_m()[vertices[k] as usize];
                let q = problem.positions_m()[vertices[0] as usize];
                let expected = (0..3)
                    .map(|j| expected_gradient[j] * (p[j] - q[j]))
                    .sum::<f64>();
                near(cell.potential_v[k] - cell.potential_v[0], expected, 2e-9);
            }
        }
        assert!(step.maximum_ampere_residual_a < 2e-8 && step.maximum_charge_residual_a < 2e-8);
        outcomes.push(step.endpoint_flux_density_t.clone());
        assert_eq!(state.time_s(), 0.0);
        near(candidate.state().time_s(), 0.1, 1e-15);
    }
    for (a, b) in outcomes[0].iter().zip(&outcomes[1]) {
        for j in 0..3 {
            near(a[j], b[j], 2e-9);
        }
    }
}

#[test]
fn charge_eliminated_cube_mode_matches_independent_decay_and_second_order_time_refinement() {
    // Unit cube's six ordered tetrahedra have ONE free edge (the body diagonal).
    // Exact integration gives K=4/mu. Projecting W onto nodal gradients gives
    // ||W||²=1/5, ||grad q||²=1/8, hence M=3*sigma/40.
    // q at cube vertices is 0,1/4,...,1/4,1/2. These constants are analytic,
    // not extracted from the production matrices.
    let (problem, sources) = cube(1, None);
    let edge = diagonal(&problem);
    let old = vec![[4.0; 3]; 6];
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let initial = problem
        .solve(
            &old,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    near(initial.potential_wb()[edge], 1.0, 2e-10);
    let state = SpatialInduction::new(&initial, sources.clone(), preparation(), 2e-9, &mut || {
        false
    })
    .unwrap();
    let (m, k, dt): (f64, f64, f64) = (3.0 * 5.0 / 40.0, 4.0 / 2.0, 0.04);
    let first = state
        .step_candidate(&zero, &velocity, dt, controls(), &mut || false)
        .unwrap();
    let a1 = m / (m + 0.5 * k * dt);
    near(first.state().potential_wb()[edge], a1, 2e-10);
    near(first.step().heat_j, m * (a1 - 1.0).powi(2) / dt, 2e-10);
    near(first.step().impressed_work_j, 0.5 * k * (a1 - 1.0), 2e-10);
    let (state, _) = first.into_parts();
    let second = state
        .step_candidate(&zero, &velocity, dt, controls(), &mut || false)
        .unwrap();
    let a2 = (m - 0.5 * k * dt) / (m + 0.5 * k * dt) * a1;
    near(second.state().potential_wb()[edge], a2, 2e-10);
    near(second.step().heat_j, 0.5 * k * (a1 * a1 - a2 * a2), 2e-10);
    near(second.step().impressed_work_j, 0.0, 1e-14);
    let rate = (a2 - a1) / dt;
    for c in &second.step().cells {
        let tet = problem.mesh().tets[c.cell];
        let q = |v: u32| {
            if v == 0 {
                0.0
            } else if v == 7 {
                0.5
            } else {
                0.25
            }
        };
        for i in 0..4 {
            near(
                c.potential_v[i] - c.potential_v[0],
                -rate * (q(tet[i]) - q(tet[0])),
                2e-9,
            );
        }
    }
    // Independently solve m a'+k a=k*(1-t/T), a(0)=1.
    let total = 0.2;
    let lambda = k / m;
    let exact = (1.0 - (-lambda * total).exp()) / (lambda * total);
    let mut errors = vec![];
    for steps in [4, 8, 16] {
        let mut state =
            SpatialInduction::new(&initial, sources.clone(), preparation(), 2e-9, &mut || {
                false
            })
            .unwrap();
        for i in 1..=steps {
            let amplitude = 4.0 * (1.0 - i as f64 / steps as f64);
            let source = vec![[amplitude; 3]; 6];
            state = state
                .step_candidate(
                    &source,
                    &velocity,
                    total / steps as f64,
                    controls(),
                    &mut || false,
                )
                .unwrap()
                .into_parts()
                .0;
        }
        errors.push((state.potential_wb()[edge] - exact).abs());
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!((1.8..2.2).contains(&order), "{errors:?}, order {order}");
    }
}

#[test]
fn three_dimensional_radial_motion_changes_the_field_and_closes_independent_lorentz_work() {
    let (problem, sources) = cube(1, None);
    let edge = diagonal(&problem);
    let source = vec![[4.0; 3]; 6];
    let initial = problem
        .solve(
            &source,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let alpha = 0.6;
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| alpha * x))
        .collect();
    let (dt, m, k): (f64, f64, f64) = (0.03, 3.0 * 5.0 / 40.0, 4.0 / 2.0);
    // On this exact mode, projection of v cross curl(W) is -2*alpha*(W-grad q).
    // Thus m a'+(k+2*alpha*m)a=j. Heat and mechanical work are independent
    // closed-form quadratic integrals in that same projected field.
    let effective = k + 2.0 * alpha * m;
    let a1 = (m - 0.5 * dt * effective + dt * k) / (m + 0.5 * dt * effective);
    let amid = 0.5 * (1.0 + a1);
    let rate = (a1 - 1.0) / dt;
    let conductor_rate = rate + 2.0 * alpha * amid;
    let candidate = state
        .step_candidate(&source, &velocity, dt, controls(), &mut || false)
        .unwrap();
    let s = candidate.step();
    near(candidate.state().potential_wb()[edge], a1, 2e-10);
    assert!((a1 - 1.0).abs() > 0.01);
    near(s.heat_j, m * conductor_rate.powi(2) * dt, 2e-10);
    near(
        s.mechanical_work_j,
        -2.0 * alpha * amid * m * conductor_rate * dt,
        2e-10,
    );
    near(s.impressed_work_j, k * (a1 - 1.0), 2e-10);
    near(
        s.stored_after_j - s.stored_before_j,
        0.5 * k * (a1 * a1 - 1.0),
        2e-10,
    );
    near(s.energy_residual_j, 0.0, 2e-10);
    assert!(s.heat_j > 0.0 && s.mechanical_work_j.abs() > 1e-5);
    for c in &s.cells {
        let tet = problem.mesh().tets[c.cell];
        let q = |v: u32| {
            if v == 0 {
                0.0
            } else if v == 7 {
                0.5
            } else {
                0.25
            }
        };
        let r = |v: u32| {
            if v == 0 {
                0.0
            } else if v == 7 {
                1.0
            } else {
                -0.5
            }
        };
        for i in 0..4 {
            let phi = -rate * (q(tet[i]) - q(tet[0])) + alpha * amid * (r(tet[i]) - r(tet[0]));
            near(c.potential_v[i] - c.potential_v[0], phi, 2e-9);
        }
    }
}

#[test]
fn insulating_separation_preserves_independent_floating_conductors() {
    let (positions, tets) = grid(3);
    let (metal, conductor) = source(2.0, Some(5.0), None);
    let (air, _) = source(2.0, None, None);
    let mut materials = vec![];
    let mut sources = vec![];
    for tet in &tets {
        let x = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
        let conducting = x < 1.0 / 3.0 || x > 2.0 / 3.0;
        materials.push(if conducting {
            metal.clone()
        } else {
            air.clone()
        });
        sources.push(if conducting { conductor.clone() } else { None });
    }
    let problem =
        SpatialVectorMagnetostatics::new(positions, tets, materials, assembly(), &mut || false)
            .unwrap();
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let b = [0.0, 0.0, 0.8];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, b, false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let velocity = vec![[0.7, 0.0, 0.0]; problem.mesh().vertex_count];
    let candidate = state
        .step_candidate(&zero, &velocity, 0.1, controls(), &mut || false)
        .unwrap();
    let components: BTreeSet<_> = candidate.step().cells.iter().map(|c| c.component).collect();
    assert_eq!(components.len(), 2);
    near(candidate.step().heat_j, 0.0, 2e-16);
    for c in &candidate.step().cells {
        for j in c.current_density_a_m2 {
            for x in j {
                near(x, 0.0, 2e-8);
            }
        }
    }
    assert!(candidate.step().maximum_charge_residual_a < 2e-8);
}

#[test]
fn permanent_magnet_recoil_storage_survives_induction_and_source_limits_remain_active() {
    let (problem, sources) = cube(1, Some(0.4));
    let edge = diagonal(&problem);
    let old = vec![[4.0; 3]; 6];
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let initial = problem
        .solve(
            &old,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let candidate = state
        .step_candidate(&zero, &velocity, 0.04, controls(), &mut || false)
        .unwrap();
    let a0 = state.potential_wb()[edge];
    let a1 = candidate.state().potential_wb()[edge];
    let k = 4.0 / 2.0;
    // Spatially constant Br has zero pairing with the internal curl mode.
    // Its fixed-remanence storage reference contributes Br²/(2mu) separately.
    near(
        candidate.step().stored_before_j,
        0.5 * k * a0 * a0 + 0.4 * 0.4 / (2.0 * 2.0),
        2e-10,
    );
    near(
        candidate.step().stored_after_j,
        0.5 * k * a1 * a1 + 0.4 * 0.4 / (2.0 * 2.0),
        2e-10,
    );
    near(candidate.step().energy_residual_j, 0.0, 2e-10);
    let too_large = vec![[2000.0; 3]; 6];
    assert!(matches!(
        state.step_candidate(&too_large, &velocity, 0.1, controls(), &mut || false),
        Err(SpatialMagneticError::Recoil { .. })
    ));
    assert_eq!(state.time_s(), 0.0);
}

#[test]
fn source_budget_current_and_late_cancellation_refuse_without_partial_history() {
    let (problem, sources) = cube(1, None);
    let source = vec![[4.0; 3]; 6];
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let initial = problem
        .solve(
            &source,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let mut changed = sources.clone();
    changed[0] = source_for_refusal();
    assert!(SpatialInduction::new(&initial, changed, preparation(), 2e-9, &mut || false).is_err());
    let mut too_small = preparation();
    too_small.maximum_assembly_terms = 1;
    assert!(
        SpatialInduction::new(&initial, sources.clone(), too_small, 2e-9, &mut || false).is_err()
    );
    let state =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let before = state.potential_wb().to_vec();
    let mut bad = controls();
    bad.maximum_krylov_scalars = 1;
    assert!(
        state
            .step_candidate(&zero, &velocity, 0.1, bad, &mut || false)
            .is_err()
    );
    let mut exhausted = controls();
    exhausted.restart = 1;
    exhausted.maximum_cycles = 1;
    assert!(matches!(
        state.step_candidate(&zero, &velocity, 0.1, exhausted, &mut || false),
        Err(SpatialMagneticError::Nonconvergence(_))
    ));
    let mut broken = zero.clone();
    broken[0] = [1.0, 2.0, 3.0];
    assert!(
        state
            .step_candidate(&broken, &velocity, 0.1, controls(), &mut || false)
            .is_err()
    );
    let moving = vec![[10.0; 3]; 8];
    let mut peclet = controls();
    peclet.maximum_cell_magnetic_peclet = 0.01;
    assert!(
        state
            .step_candidate(&zero, &moving, 0.1, peclet, &mut || false)
            .is_err()
    );
    let mut calls = 0;
    let successful = state
        .step_candidate(&zero, &velocity, 0.1, controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 30);
    let mut counter = 0;
    assert!(matches!(
        state.step_candidate(&zero, &velocity, 0.1, controls(), &mut || {
            counter += 1;
            counter == calls - 1
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    assert_eq!(state.potential_wb(), before);
    assert_eq!(state.time_s(), 0.0);
    let repeated = state
        .step_candidate(&zero, &velocity, 0.1, controls(), &mut || false)
        .unwrap();
    assert_eq!(
        successful.state().potential_wb(),
        repeated.state().potential_wb()
    );
    assert_eq!(
        successful.step().heat_j.to_bits(),
        repeated.step().heat_j.to_bits()
    );
}
fn source_for_refusal() -> Option<ResolvedConductorState> {
    source(2.0, Some(6.0), None).1
}
