//! G1/G3 implicit derivatives of manufactured coupled convection laws.

use super::*;
use fs_conduction::ThermalBc;
use fs_conduction::transient::enthalpy::EnthalpyStepSolution;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyStepLinearization;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};

const NAMES: [&str; 2] = ["natural", "air"];
const H_SLOPES: [f64; 2] = [0.04, -0.015];
const CARRY: [f64; 4] = [0.2, -0.1, 0.3, 0.4];
const WALL: [f64; 2] = [0.15, -0.1];
const HEAT: [f64; 2] = [0.03, -0.02];

fn references() -> Vec<Vec<f64>> {
    vec![vec![0.0, 0.15], vec![0.2, 0.35]]
}

fn means(mesh: &ConductionMesh, temperature: &[f64]) -> [f64; 2] {
    let mut areas = [0.0; 2];
    let mut sums = [0.0; 2];
    for face in mesh.boundary() {
        let i = usize::from(face.centroid[2] != 0.0);
        areas[i] += face.area;
        for &vertex in &face.vertices {
            sums[i] += face.area / 3.0 * temperature[vertex as usize];
        }
    }
    [sums[0] / areas[0], sums[1] / areas[1]]
}

fn boundary(fixture: &Fixture, temperature: &[f64]) -> ThermalBoundary {
    let mean = means(&fixture.mesh, temperature);
    let h = [
        3.0 + H_SLOPES[0] * (mean[0] - 350.0),
        4.0 + H_SLOPES[1] * (mean[1] - 350.0),
    ];
    let r = [
        325.0 + 0.15 * (mean[1] - 350.0),
        335.0 + 0.2 * (mean[0] - 350.0) + 0.35 * (mean[1] - 350.0),
    ];
    ThermalBoundaryBuilder::new(&fixture.mesh)
        .region(
            NAMES[0],
            |face| face.centroid[2] == 0.0,
            ThermalBc::robin(h[0], r[0]).unwrap(),
        )
        .unwrap()
        .remainder(NAMES[1], ThermalBc::robin(h[1], r[1]).unwrap())
        .unwrap()
        .finish()
        .unwrap()
}

fn temperature(fixture: &Fixture, h: &[f64]) -> Vec<f64> {
    h.iter()
        .map(|&h| {
            fixture
                .curve
                .state_at_specific_enthalpy(h)
                .unwrap()
                .temperature_k()
        })
        .collect()
}

fn solve(
    fixture: &Fixture,
    cx: &Cx<'_>,
    old: &[f64],
    source: &ScalarField,
) -> (EnthalpyStepSolution, ThermalBoundary) {
    let mut wall = temperature(fixture, old);
    for _ in 0..64 {
        let bc = boundary(fixture, &wall);
        let problem = ConductionProblem {
            boundary: &bc,
            ..fixture.problem(source)
        };
        let accepted = fixture
            .stepper(cx)
            .advance(cx, problem, None, old, DT, config())
            .unwrap();
        let change = wall
            .iter()
            .zip(&accepted.temperature)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        if change < 1e-11 {
            return (accepted, bc);
        }
        wall = accepted.temperature;
    }
    panic!("manufactured convection fixed point did not converge");
}

fn linearize<'m>(fixture: &'m Fixture, cx: &Cx<'_>, old: &[f64]) -> EnthalpyStepLinearization<'m> {
    let source = ScalarField::Nodal(SOURCE.to_vec());
    let (accepted, boundary) = solve(fixture, cx, old, &source);
    fixture
        .stepper(cx)
        .linearize_accepted(
            cx,
            ConductionProblem {
                boundary: &boundary,
                ..fixture.problem(&source)
            },
            None,
            old,
            DT,
            config(),
            accepted,
        )
        .unwrap()
}

fn residual(fixture: &Fixture, cx: &Cx<'_>, h: &[f64]) -> Vec<f64> {
    let t = temperature(fixture, h);
    let bc = boundary(fixture, &t);
    let source = ScalarField::Nodal(SOURCE.to_vec());
    let system =
        fs_conduction::assemble_operator(cx, &fixture.mesh, &bc, &fixture.material, &source, &t)
            .unwrap();
    let mut residual = vec![0.0; 4];
    system.operator.spmv(&t, &mut residual);
    let stepper = fixture.stepper(cx);
    for i in 0..4 {
        residual[i] = stepper.reference_nodal_masses_kg()[i] * (h[i] - OLD[i])
            + DT * (residual[i] - system.load[i]);
    }
    residual
}

fn objective(fixture: &Fixture, cx: &Cx<'_>, old: &[f64], q: &[f64]) -> f64 {
    let source = ScalarField::Nodal(q.to_vec());
    let (accepted, bc) = solve(fixture, cx, old, &source);
    let mean = means(&fixture.mesh, &accepted.temperature);
    let mut value = dot(&CARRY, &accepted.specific_enthalpy_j_kg)
        + accepted
            .temperature
            .iter()
            .zip(WEIGHTS)
            .map(|(t, w)| w * (t - 350.0))
            .sum::<f64>()
        + WALL
            .iter()
            .zip(mean)
            .map(|(w, t)| w * (t - 350.0))
            .sum::<f64>();
    for (slot, face) in fixture.mesh.boundary().iter().enumerate() {
        let region = bc.region_for(slot).unwrap();
        let i = usize::from(bc.region_names()[region] != NAMES[0]);
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(h),
            t_ref: ScalarField::Uniform(r),
        } = bc.conditions()[region]
        else {
            panic!("manufactured uniform Robin row");
        };
        for &vertex in &face.vertices {
            value += HEAT[i] * h * (face.area / 3.0) * (accepted.temperature[vertex as usize] - r);
        }
    }
    value
}

#[test]
fn coupled_robin_enthalpy_tangent_and_physical_heat_goals_match_resolves() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        let base = linearize(&fixture, cx, &OLD);
        let retained = base.primal().clone();
        let full = base
            .with_robin_feedback(cx, &NAMES, &references(), &H_SLOPES, 32)
            .unwrap();
        assert_eq!(
            full.primal().specific_enthalpy_j_kg,
            retained.specific_enthalpy_j_kg
        );
        assert_eq!(full.primal().temperature, retained.temperature);
        let u = [0.7, -1.2, 0.6, 0.3];
        let v = [-0.5, 0.8, 1.4, -0.2];
        let jv = full.apply_jacobian(cx, &v).unwrap();
        let jtu = full.apply_jacobian_transpose(cx, &u).unwrap();
        assert!((dot(&u, &jv) - dot(&jtu, &v)).abs() < 1e-13);
        let epsilon = 1e-4;
        let shifted = |sign| {
            retained
                .specific_enthalpy_j_kg
                .iter()
                .zip(v)
                .map(|(h, v)| h + sign * epsilon * v)
                .collect::<Vec<_>>()
        };
        let high = residual(&fixture, cx, &shifted(1.0));
        let low = residual(&fixture, cx, &shifted(-1.0));
        for i in 0..4 {
            assert!((jv[i] - (high[i] - low[i]) / (2.0 * epsilon)).abs() < 2e-8);
        }
        let response = full.robin_response(cx, &NAMES, linear()).unwrap();
        let gradient = response
            .pullback(cx, &CARRY, &WEIGHTS, &WALL, &HEAT)
            .unwrap();
        assert!(gradient.transport.relative_residual <= linear().tolerance);
        for i in 0..4 {
            let history_fd = five_point(|delta| {
                let mut old = OLD;
                old[i] += delta;
                objective(&fixture, cx, &old, &SOURCE)
            });
            let source_fd = five_point(|delta| {
                let mut source = SOURCE;
                source[i] += delta;
                objective(&fixture, cx, &OLD, &source)
            });
            assert!(
                (gradient.transport.previous_specific_enthalpy[i] - history_fd).abs() < 2e-8,
                "history {i}: {} vs {history_fd}",
                gradient.transport.previous_specific_enthalpy[i]
            );
            assert!(
                (gradient.transport.source_density[i] - source_fd).abs() < 2e-8,
                "source {i}: {} vs {source_fd}",
                gradient.transport.source_density[i]
            );
        }
    });
}

#[test]
fn coupled_attachment_preserves_plateau_carry_and_bound_heterogeneous_density() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        let base = linearize(&fixture, cx, &[200.0; 4]);
        let full = base
            .with_robin_feedback(cx, &NAMES, &references(), &H_SLOPES, 32)
            .unwrap();
        assert_eq!(full.primal().temperature, [350.0; 4]);
        let seed = full.temperature_pullback(cx, &WEIGHTS).unwrap();
        assert_eq!(seed, [0.0; 4]);
        let gradient = full.pullback(cx, &[1.0; 4], linear()).unwrap();
        for &value in &gradient.previous_specific_enthalpy {
            assert!((value - 1.0).abs() < 1e-13);
        }
        for &value in &gradient.source_density {
            assert!((value - DT / 2.0).abs() < 1e-13);
        }

        let source = ScalarField::Nodal(SOURCE.to_vec());
        let (accepted, bc) = solve(&fixture, cx, &OLD, &source);
        let materials = [ReferenceEnthalpyMaterial {
            curve: &fixture.curve,
            reference_density_kg_m3: 2.0,
        }];
        let assigned = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &fixture.mesh,
            &materials,
            &[0],
            EnthalpyBudget {
                max_vertices: 4,
                max_elements: 1,
            },
        )
        .unwrap();
        let linearization = assigned
            .linearize_accepted(
                cx,
                ConductionProblem {
                    boundary: &bc,
                    ..fixture.problem(&source)
                },
                None,
                &OLD,
                DT,
                config(),
                accepted,
            )
            .unwrap()
            .with_robin_feedback(cx, &NAMES, &references(), &H_SLOPES, 32)
            .unwrap();
        let seed = linearization
            .transport()
            .temperature_pullback(cx, &WEIGHTS)
            .unwrap();
        let gradient = linearization.pullback(cx, &seed, linear()).unwrap();
        let difference: Vec<_> = linearization
            .primal()
            .specific_enthalpy_j_kg
            .iter()
            .zip(OLD)
            .map(|(new, old)| new - old)
            .collect();
        let expected =
            -fixture.mesh.element_volume(0) / 4.0 * dot(&gradient.transport.adjoint, &difference);
        assert!((gradient.element_reference_density[0] - expected).abs() < 1e-13);
        assert_eq!(
            gradient.material_reference_density,
            gradient.element_reference_density
        );
    });
}

#[test]
fn coupled_attachment_refuses_bad_laws_repeated_attachment_and_work_overruns() {
    let fixture = Fixture::new();
    with_cx(|cx| {
        for (names, matrix, slopes, cap) in [
            (NAMES.to_vec(), references(), H_SLOPES.to_vec(), 31),
            (
                vec!["natural", "natural"],
                references(),
                H_SLOPES.to_vec(),
                32,
            ),
            (
                vec!["natural", "missing"],
                references(),
                H_SLOPES.to_vec(),
                32,
            ),
            (NAMES.to_vec(), vec![vec![0.0]; 2], H_SLOPES.to_vec(), 32),
            (NAMES.to_vec(), references(), vec![0.0], 32),
            (
                NAMES.to_vec(),
                vec![vec![f64::NAN; 2]; 2],
                H_SLOPES.to_vec(),
                32,
            ),
            (NAMES.to_vec(), references(), vec![f64::INFINITY; 2], 32),
        ] {
            assert!(
                linearize(&fixture, cx, &OLD)
                    .with_robin_feedback(cx, &names, &matrix, &slopes, cap)
                    .is_err()
            );
        }
        let full = linearize(&fixture, cx, &OLD)
            .with_robin_feedback(cx, &NAMES, &references(), &H_SLOPES, 32)
            .unwrap();
        assert!(
            full.with_robin_feedback(cx, &NAMES, &references(), &H_SLOPES, 64)
                .is_err()
        );
        let base = linearize(&fixture, cx, &OLD);
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                base.with_robin_feedback(cancelled, &NAMES, &references(), &H_SLOPES, 32),
                Err(EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
                    ConductionError::Cancelled { .. }
                )))
            ));
        });
    });
}
