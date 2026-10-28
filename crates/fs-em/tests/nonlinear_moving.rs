//! Independent P1, quadrature and time-domain oracles for nonlinear motion.
//! Synthetic source laws establish numerical behavior, not material validation.
#[path = "winding_ohmic/common.rs"]
mod linear_material;
#[path = "nonlinear_eddy/common.rs"]
mod material;

use fs_em::planar::{
    eddy::{
        EddyConductor,
        moving::MovingEddyControls,
        nonlinear::{NonlinearEddyError, PlanarNonlinearEddy, moving::NonlinearMovingEddyControls},
    },
    nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
    permanent::PlanarRecoilMaterial,
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
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
use fs_rep_mesh::{Metric2, TriComplex2};

fn controls() -> NonlinearMovingEddyControls {
    NonlinearMovingEddyControls {
        moving: MovingEddyControls {
            relative_tolerance: 1e-12,
            absolute_residual_a_m: 1e-13,
            restart: 8,
            max_cycles: 30,
            max_krylov_scalars: 10000,
            max_cell_magnetic_peclet: 100.0,
            energy_tolerance_j: 1e-9,
            current_balance_tolerance_a: 1e-8,
        },
        absolute_residual_tolerance_a_m: 1e-11,
        relative_residual_tolerance: 1e-11,
        max_newton_iterations: 30,
        max_backtracks: 30,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn rectangle() -> TriComplex2 {
    let vertices: Vec<_> = (0..3)
        .flat_map(|j| (0..4).map(move |i| [f64::from(i), f64::from(j) / 2.0]))
        .collect();
    let mut faces = Vec::new();
    for j in 0..2 {
        for i in 0..3 {
            let a = 4 * j + i;
            faces.extend([[a, a + 1, a + 5], [a, a + 5, a + 4]]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "nonlinear-moving-two-modes",
        vertices,
        faces,
        Metric2::planar(0.3).unwrap(),
    )
    .unwrap()
}
fn boundary(mesh: &TriComplex2, affine: f64) -> Vec<(usize, f64)> {
    mesh.vertices()
        .iter()
        .enumerate()
        .filter(|(v, _)| *v != 5 && *v != 6)
        .map(|(v, p)| (v, affine * p[0]))
        .collect()
}
fn state<'a>(
    problem: &'a PlanarNonlinearMagnetostatics,
    source: &MaterialCard,
    current: f64,
    boundary: &[(usize, f64)],
    guess: &[f64],
) -> PlanarNonlinearEddy<'a> {
    let initial = problem
        .solve(
            &vec![current; problem.mesh().faces().len()],
            boundary,
            guess,
            material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    PlanarNonlinearEddy::new(
        &initial,
        material::conductors(
            source,
            problem.mesh().faces().len(),
            300.0,
            ElectricalPropertyBasis::Conductivity,
        ),
        &[],
        boundary,
        0.0,
        material::controls(),
        &mut || false,
    )
    .unwrap()
}

// Elementary area/shape formulas, independent of production Hodge/gradient code.
fn geometry(points: [[f64; 2]; 3], depth: f64) -> (f64, [[f64; 2]; 3]) {
    let det = (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
        - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1]);
    (
        det.abs() * depth / 2.0,
        std::array::from_fn(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            [
                (points[j][1] - points[k][1]) / det,
                (points[k][0] - points[j][0]) / det,
            ]
        }),
    )
}
fn law(b: f64, saturating: bool) -> (f64, f64) {
    if saturating && b > 1.0 {
        (4.0 * b - 3.0, 0.5 + (b - 1.0) + 2.0 * (b - 1.0).powi(2))
    } else {
        (b, 0.5 * b * b)
    }
}
fn gradient(potential: &[f64], face: [usize; 3], g: [[f64; 2]; 3]) -> [f64; 2] {
    std::array::from_fn(|axis| (0..3).map(|i| potential[face[i]] * g[i][axis]).sum())
}

fn dense_reference(
    mesh: &TriComplex2,
    a0: &[f64],
    velocity: &[[f64; 2]],
    dt: f64,
    saturating: bool,
) -> Vec<f64> {
    let n = a0.len();
    let sigma = 18.0;
    let mut c = vec![vec![0.0; n]; n];
    let mut d = c.clone();
    let mut moments = vec![0.0; n];
    let mut advected = vec![0.0; n];
    let mut weight = 0.0;
    for face in mesh.faces() {
        let f = face.map(|v| v as usize);
        let (volume, g) = geometry(f.map(|v| mesh.vertices()[v]), 0.3);
        weight += sigma * volume;
        for i in 0..3 {
            moments[f[i]] += sigma * volume / 3.0;
            advected[f[i]] += sigma * volume / 3.0
                * f.iter()
                    .map(|&v| velocity[v][0] * g[i][0] + velocity[v][1] * g[i][1])
                    .sum::<f64>();
            for j in 0..3 {
                c[f[i]][f[j]] += sigma * volume / 12.0 * if i == j { 2.0 } else { 1.0 };
                for q in 0..3 {
                    d[f[i]][f[j]] += sigma * volume / 12.0
                        * if i == q { 2.0 } else { 1.0 }
                        * (velocity[f[q]][0] * g[j][0] + velocity[f[q]][1] * g[j][1]);
                }
            }
        }
    }
    for i in 0..n {
        for j in 0..n {
            c[i][j] -= moments[i] * moments[j] / weight;
            d[i][j] -= moments[i] * advected[j] / weight;
        }
    }
    let residual = |a: &[f64]| {
        let mut r: Vec<_> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| c[i][j] * (a[j] - a0[j]) / dt + d[i][j] * a[j])
                    .sum::<f64>()
            })
            .collect();
        for face in mesh.faces() {
            let f = face.map(|v| v as usize);
            let (volume, g) = geometry(f.map(|v| mesh.vertices()[v]), 0.3);
            let grad = gradient(a, f, g);
            let b = grad[0].hypot(grad[1]);
            let nu = if b == 0.0 {
                1.0
            } else {
                law(b, saturating).0 / b
            };
            for i in 0..3 {
                r[f[i]] += volume * nu * (grad[0] * g[i][0] + grad[1] * g[i][1]);
            }
        }
        [r[5], r[6]]
    };
    let mut a = a0.to_vec();
    for _ in 0..20 {
        let r = residual(&a);
        if r[0].hypot(r[1]) < 1e-12 {
            return a;
        }
        // Independent finite-difference 2x2 Newton oracle: no production
        // constitutive tangent, FEEC assembly, or Krylov implementation.
        let mut jac = [[0.0; 2]; 2];
        for (col, v) in [5, 6].into_iter().enumerate() {
            let mut lo = a.clone();
            let mut hi = a.clone();
            lo[v] -= 1e-6;
            hi[v] += 1e-6;
            let lo = residual(&lo);
            let hi = residual(&hi);
            for row in 0..2 {
                jac[row][col] = (hi[row] - lo[row]) / 2e-6;
            }
        }
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        a[5] -= (r[0] * jac[1][1] - jac[0][1] * r[1]) / det;
        a[6] -= (jac[0][0] * r[1] - r[0] * jac[1][0]) / det;
    }
    panic!("independent two-mode reference did not converge")
}

#[test]
fn saturating_motion_matches_independent_nonlinear_field_and_three_way_quadrature() {
    let mut results = Vec::new();
    for saturating in [false, true] {
        let knots = if saturating {
            vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]
        } else {
            vec![(0.0, 0.0), (4.0, 4.0)]
        };
        let source = material::source(&knots, 18.0);
        let problem = material::problem(&source, rectangle());
        let initial = state(
            &problem,
            &source,
            12.0,
            &boundary(problem.mesh(), 0.0),
            &[0.0; 12],
        );
        let velocity: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|p| [0.2 + 0.1 * p[1], 0.0])
            .collect();
        let dt = 0.08;
        let candidate = initial
            .step_with_eulerian_velocity_candidate(
                dt,
                &velocity,
                &[0.0; 12],
                controls(),
                &mut || false,
            )
            .unwrap();
        let reference = dense_reference(
            problem.mesh(),
            initial.potential_wb_m(),
            &velocity,
            dt,
            saturating,
        );
        for (&a, &b) in candidate.state().potential_wb_m().iter().zip(&reference) {
            near(a, b);
        }
        let rate: Vec<_> = reference
            .iter()
            .zip(initial.potential_wb_m())
            .map(|(a, b)| (a - b) / dt)
            .collect();
        let mut mean = 0.0;
        let mut total_volume = 0.0;
        let mut cells = Vec::new();
        let mut before_energy = 0.0;
        let mut after_energy = 0.0;
        let mut numerical = 0.0;
        for face in problem.mesh().faces() {
            let f = face.map(|v| v as usize);
            let (volume, g) = geometry(f.map(|v| problem.mesh().vertices()[v]), 0.3);
            let before = gradient(initial.potential_wb_m(), f, g);
            let after = gradient(&reference, f, g);
            let (h, psi) = law(after[0].hypot(after[1]), saturating);
            let old_psi = law(before[0].hypot(before[1]), saturating).1;
            before_energy += volume * old_psi;
            after_energy += volume * psi;
            let magnitude = after[0].hypot(after[1]);
            let nu = if magnitude == 0.0 { 1.0 } else { h / magnitude };
            numerical += volume
                * (nu * (after[0] * (after[0] - before[0]) + after[1] * (after[1] - before[1]))
                    - (psi - old_psi));
            let u = f.map(|v| velocity[v][0] * after[0] + velocity[v][1] * after[1]);
            let r: [f64; 3] = std::array::from_fn(|i| rate[f[i]] + u[i]);
            mean += volume * r.iter().sum::<f64>() / 3.0;
            total_volume += volume;
            cells.push((volume, after, u, r));
        }
        mean /= total_volume;
        let mut heat = 0.0;
        let mut work = 0.0;
        for (actual, (volume, grad, u, r)) in candidate.step().cells.iter().zip(cells) {
            let mut q = 0.0;
            let mut w = 0.0;
            for k in 0..3 {
                let bary: [f64; 3] =
                    std::array::from_fn(|i| if i == k { 2.0 / 3.0 } else { 1.0 / 6.0 });
                let e = mean - (0..3).map(|i| bary[i] * r[i]).sum::<f64>();
                let motion = (0..3).map(|i| bary[i] * u[i]).sum::<f64>();
                q += dt * volume / 3.0 * 18.0 * e * e;
                w += dt * volume / 3.0 * 18.0 * e * motion;
            }
            near(actual.joule.heat_j, q);
            near(actual.mechanical_work_j, w);
            for i in 0..3 {
                near(actual.joule.current_density_a_m2[i], 18.0 * (mean - r[i]));
                near(
                    actual.laboratory_electric_field_v_m[i] + actual.motional_electric_field_v_m[i],
                    mean - r[i],
                );
                for axis in 0..2 {
                    near(
                        actual.lorentz_force_density_n_m3[i][axis],
                        18.0 * (mean - r[i]) * grad[axis],
                    );
                }
            }
            heat += q;
            work += w;
        }
        let step = candidate.step();
        near(step.stored_before_j, before_energy);
        near(step.stored_after_j, after_energy);
        near(step.numerical_dissipation_j, numerical);
        near(step.heat_j, heat);
        near(step.mechanical_work_j, work);
        near(step.energy_residual_j, 0.0);
        near(after_energy - before_energy + heat + work + numerical, 0.0);
        near(step.region_currents_a[0].1, 0.0);
        assert!(heat > 0.0 && numerical > 0.0 && work < -1e-5);
        assert!(step.report.newton_iterations > 0 && step.report.linear_iterations > 0);
        assert!(step.report.residual_norm_a_m <= step.report.residual_target_a_m);
        assert_eq!(initial.time_s(), 0.0);
        results.push(reference[5]);
    }
    assert!(
        (results[0] - results[1]).abs() > 0.1,
        "source H(B) must change actual moving induction"
    );
}

#[test]
fn zero_motion_continues_the_original_saturating_induction_and_its_heat() {
    let source = material::source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], 18.0);
    let problem = material::problem(&source, rectangle());
    let mut moving = state(
        &problem,
        &source,
        12.0,
        &boundary(problem.mesh(), 0.0),
        &[0.0; 12],
    );
    let mut stationary = moving.clone();
    for tick in 1..=4 {
        let t = f64::from(tick) * 0.04;
        let a = moving
            .step_with_eulerian_velocity_candidate(
                t,
                &[[0.0; 2]; 12],
                &[0.0; 12],
                controls(),
                &mut || false,
            )
            .unwrap();
        let b = stationary
            .step_to_candidate(
                t,
                stationary.conductors().to_vec(),
                &[0.0; 12],
                material::controls(),
                &mut || false,
            )
            .unwrap();
        for (&x, &y) in a
            .state()
            .potential_wb_m()
            .iter()
            .zip(b.state().potential_wb_m())
        {
            near(x, y);
        }
        near(a.step().heat_j, b.step().heat_j);
        near(
            a.step().numerical_dissipation_j,
            b.step().numerical_dissipation_j,
        );
        near(a.step().mechanical_work_j, 0.0);
        moving = a.into_parts().0;
        stationary = b.into_parts().0;
    }
}

#[test]
fn saturated_uniform_field_retains_floating_charge_cancellation() {
    let source = material::source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], 18.0);
    let problem = material::problem(&source, rectangle());
    let potential: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .map(|p| 1.5 * p[0])
        .collect();
    let initial = state(
        &problem,
        &source,
        0.0,
        &boundary(problem.mesh(), 1.5),
        &potential,
    );
    let next = initial
        .step_with_eulerian_velocity_candidate(
            0.1,
            &[[0.4, 0.0]; 12],
            &[0.0; 12],
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in next.state().potential_wb_m().iter().zip(&potential) {
        near(a, b);
    }
    near(next.step().heat_j, 0.0);
    near(next.step().mechanical_work_j, 0.0);
    for cell in &next.step().cells {
        for i in 0..3 {
            near(cell.joule.current_density_a_m2[i], 0.0);
            near(cell.laboratory_electric_field_v_m[i], 0.6);
            near(cell.motional_electric_field_v_m[i], -0.6);
        }
    }
}

fn recoil_case(sign: f64, limit: f64) -> (PlanarNonlinearMagnetostatics, ResolvedConductorState) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let values = [
        (
            RECOIL_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1.0,
        ),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, 0.8),
        (
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            -limit,
        ),
        (
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            limit,
        ),
        (
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            limit,
        ),
        (basis.property_name(), basis.dimensions(), 18.0),
    ];
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
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
                    source: "manufactured stationary recoil bias".into(),
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
            chemistry: "synthetic conductive magnet".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let magnet = PlanarRecoilMaterial::new(
        ResolvedRecoilMaterial::try_from_material(material).unwrap(),
        [sign, 0.0],
    )
    .unwrap();
    let metal_source = linear_material::card(1.0 / 18.0, 0.0, 1.0, basis);
    let metal = linear_material::bundle(&metal_source, 300.0, basis);
    let conductor =
        ResolvedConductorState::try_from_material_with_basis(metal.clone(), basis).unwrap();
    let mesh = TriComplex2::from_indexed_triangles(
        "moving-fixed-recoil",
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        vec![[0, 1, 2], [1, 3, 2]],
        Metric2::planar(0.2).unwrap(),
    )
    .unwrap();
    (
        PlanarNonlinearMagnetostatics::new(
            mesh,
            vec![
                IsotropicMagneticMaterial::Recoil(magnet),
                IsotropicMagneticMaterial::Linear(metal),
            ],
            WeightedAssemblyLimits {
                max_cells: 10,
                max_dofs: 10,
                max_triplets: 100,
            },
            &mut || false,
        )
        .unwrap(),
        conductor,
    )
}

fn recoil_state<'a>(
    problem: &'a PlanarNonlinearMagnetostatics,
    conductor: ResolvedConductorState,
    sign: f64,
) -> PlanarNonlinearEddy<'a> {
    let boundary = [(0, 0.0), (1, 0.0)];
    let field = problem
        .solve(
            &[0.0; 2],
            &boundary,
            &[0.0, 0.0, sign * 1.6 / 3.0, sign * 0.8 / 3.0],
            material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    PlanarNonlinearEddy::new(
        &field,
        vec![EddyConductor {
            triangle: 1,
            region: 7,
            conductor,
        }],
        &[0],
        &boundary,
        0.0,
        material::controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn stationary_magnet_and_moving_conductor_match_exact_p1_ode_and_time_refinement() {
    for sign in [-1.0, 1.0] {
        let (problem, conductor) = recoil_case(sign, 5.0);
        let initial = recoil_state(&problem, conductor, sign);
        // Magnet vertices 0,1,2 are stationary. In the metal,
        // v=0.2*N3*(1,-1) has zero divergence and is tangent to the magnet interface.
        let velocity = [[0.0; 2], [0.0; 2], [0.0; 2], [0.2, -0.2]];
        let dt = 0.1;
        let next = initial
            .step_with_eulerian_velocity_candidate(
                dt,
                &velocity,
                &[0.0; 2],
                controls(),
                &mut || false,
            )
            .unwrap();
        // Independent two-node P1 system, after dividing by cell volume:
        // C=[[1,-1/2],[-1/2,1]], K=2*C, D=[[0.1,0],[-0.2,0]], f=[signed Br,0].
        // Therefore u=2*Br/3 stays fixed and wdot+2*w=2*Br/3+0.2*u.
        let u = sign * 1.6 / 3.0;
        let w0 = sign * 0.8 / 3.0;
        let exact = (w0 + dt * (2.0 * w0 + 0.2 * u)) / (1.0 + 2.0 * dt);
        near(next.state().potential_wb_m()[2], u);
        near(next.state().potential_wb_m()[3], exact);
        near(
            next.step().stored_after_j,
            next.step().stored_before_j + 0.1 * (exact - w0).powi(2),
        );
        near(
            next.step().numerical_dissipation_j,
            0.1 * (exact - w0).powi(2),
        );
        near(next.step().energy_residual_j, 0.0);
        assert!(next.step().mechanical_work_j < 0.0 && next.step().heat_j > 0.0);
        let final_time = 0.2;
        let equilibrium = w0 + 0.1 * u;
        let continuous = equilibrium + (w0 - equilibrium) * f64::exp(-2.0 * final_time);
        let mut previous = None;
        for count in [2, 4, 8] {
            let mut state = initial.clone();
            let mut numerical = 0.0;
            for tick in 1..=count {
                let next = state
                    .step_with_eulerian_velocity_candidate(
                        f64::from(tick) * final_time / f64::from(count),
                        &velocity,
                        &[0.0; 2],
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                numerical += next.step().numerical_dissipation_j;
                near(next.step().energy_residual_j, 0.0);
                state = next.into_parts().0;
            }
            near(state.potential_wb_m()[2], u);
            let error = (state.potential_wb_m()[3] - continuous).abs();
            if let Some((old_error, old_loss)) = previous {
                assert!(old_error / error > 1.8 && old_error / error < 2.1);
                assert!(numerical < old_loss);
            }
            previous = Some((error, numerical));
        }
    }
}

#[test]
fn nonlinear_source_velocity_budget_and_cancellation_refusals_preserve_history() {
    let source = material::source(&[(0.0, 0.0), (1.0, 0.1), (4.0, 12.1)], 18.0);
    let problem = material::problem(&source, rectangle());
    let original = state(
        &problem,
        &source,
        12.0,
        &boundary(problem.mesh(), 0.0),
        &[0.0; 12],
    );
    let velocity = [[0.2, 0.0]; 12];
    let saved = original.potential_wb_m().to_vec();
    let mut bad = controls();
    bad.moving.max_cell_magnetic_peclet = 5.0;
    // The low-field slope controls the whole-source Peclet bound even when
    // the current operating point lies on a much stiffer saturated segment.
    assert!(
        original
            .step_with_eulerian_velocity_candidate(0.08, &velocity, &[0.0; 12], bad, &mut || false)
            .is_err()
    );
    bad = controls();
    bad.max_newton_iterations = 0;
    assert!(
        original
            .step_with_eulerian_velocity_candidate(0.08, &velocity, &[0.0; 12], bad, &mut || false)
            .is_err()
    );
    bad = controls();
    bad.moving.max_krylov_scalars = 0;
    assert!(
        original
            .step_with_eulerian_velocity_candidate(0.08, &velocity, &[0.0; 12], bad, &mut || false)
            .is_err()
    );
    assert!(
        original
            .step_with_eulerian_velocity_candidate(
                0.08,
                &[[f64::NAN, 0.0]; 12],
                &[0.0; 12],
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        original
            .step_with_eulerian_velocity_candidate(
                0.08,
                &velocity[..11],
                &[0.0; 12],
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        original
            .step_with_eulerian_velocity_candidate(
                0.08,
                &velocity,
                &[1e6; 12],
                controls(),
                &mut || false
            )
            .is_err()
    );
    let mut polls = 0;
    original
        .step_with_eulerian_velocity_candidate(0.08, &velocity, &[0.0; 12], controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut called = 0;
    let interrupted = original.step_with_eulerian_velocity_candidate(
        0.08,
        &velocity,
        &[0.0; 12],
        controls(),
        &mut || {
            called += 1;
            called > polls / 2
        },
    );
    assert!(matches!(interrupted, Err(NonlinearEddyError::Cancelled)));
    assert_eq!(original.potential_wb_m(), saved);
    assert_eq!(original.time_s(), 0.0);
    let (magnet, conductor) = recoil_case(1.0, 0.3);
    let original = recoil_state(&magnet, conductor, 1.0);
    assert!(
        original
            .step_with_eulerian_velocity_candidate(
                0.1,
                &[[0.0; 2], [0.0; 2], [0.0; 2], [0.2, -0.2]],
                &[0.0, -1e6],
                controls(),
                &mut || false
            )
            .is_err()
    );
    near(original.potential_wb_m()[2], 1.6 / 3.0);
    assert_eq!(original.time_s(), 0.0);
    let field = original.field(&mut || false).unwrap();
    let mut conductors = original.conductors().to_vec();
    conductors.push(EddyConductor {
        triangle: 0,
        region: 7,
        conductor: ResolvedConductorState::try_from_material_with_basis(
            magnet.materials()[0].material().clone(),
            ElectricalPropertyBasis::Conductivity,
        )
        .unwrap(),
    });
    let stationary_magnet = PlanarNonlinearEddy::new(
        &field,
        conductors,
        &[],
        &[(0, 0.0), (1, 0.0)],
        0.0,
        material::controls(),
        &mut || false,
    )
    .unwrap();
    assert!(
        stationary_magnet
            .step_with_eulerian_velocity_candidate(
                0.1,
                &[[0.2, 0.0]; 4],
                &[0.0; 2],
                controls(),
                &mut || false,
            )
            .is_err(),
        "a fixed recoil source cannot silently acquire material velocity"
    );
}
