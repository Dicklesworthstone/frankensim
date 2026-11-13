//! G1/G3: canonical vector potential reaches actual oriented edge cochains.
#[path = "vector/induction.rs"]
mod induction;
use super::*;
use fs_couple::scenario_magnetic::{
    MagneticPotentialBoundary3, solve_vector_potential3, vector_potential_boundary_at,
};
use fs_em::spatial::vector_potential::{SpatialVectorMagnetostatics, VectorMagneticControls};
use fs_ga::{Quat, Vec3};
use fs_scenario::payload::VectorPayload;
use fs_scenario::{Frame, FrameId, FrameMotion};

const A_DIMS: Dims = Dims([1, 1, -2, 0, -1, 0]);

fn vector_problem() -> SpatialVectorMagnetostatics {
    let scalar = super::problem();
    SpatialVectorMagnetostatics::new(
        scalar.positions_m().to_vec(),
        scalar.mesh().tets.clone(),
        scalar.materials().to_vec(),
        common::limits(),
        &mut || false,
    )
    .unwrap()
}

fn vector_controls() -> VectorMagneticControls {
    VectorMagneticControls {
        field: common::config(),
        gauge: common::config(),
        maximum_gauge_product_terms: 10000,
        current_tolerance_a: 1e-10,
        flux_tolerance_wb: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}

fn potential_row(value: [f64; 3], frame: FrameId, parity: OrientationParity) -> BoundaryCondition {
    let meta = PayloadMeta::new(
        QuantityContract::Dimensions(A_DIMS),
        PayloadId::new("basis/frame-cartesian").unwrap(),
        frame,
        parity,
        ReferenceSemantics::Continuous,
    )
    .unwrap();
    let samples = [[0.0; 3], value]
        .map(|value| value.map(|v| QtyAny::new(v, A_DIMS)).to_vec())
        .to_vec();
    BoundaryCondition {
        region: "vector-potential".into(),
        physics: Physics::Magnetics,
        kind: BcKind::MagneticVectorPotential,
        value: Some(BcValue::Typed(Payload::Vector(
            VectorPayload::new(
                meta,
                SampleSource::table(
                    [0.0, 1.0].map(|v| QtyAny::new(v, TIME)).to_vec(),
                    samples,
                    TableInterpolation::Linear,
                    OutsideDomainPolicy::Refuse,
                )
                .unwrap(),
            )
            .unwrap(),
        ))),
        compatibility: None,
        frame: frame.0,
    }
}

fn edge_history(problem: &SpatialVectorMagnetostatics, b: [f64; 3]) -> Vec<BoundaryCondition> {
    problem
        .mesh()
        .edges
        .iter()
        .enumerate()
        .map(|(index, &[a, c])| {
            let p = problem.positions_m()[a as usize];
            let q = problem.positions_m()[c as usize];
            let mid: [f64; 3] = std::array::from_fn(|j| 0.5 * (p[j] + q[j]));
            // A(x)=B cross x / 2 is affine along each edge, so these are
            // its exact line averages, independently of the adapter's projection.
            let average = [
                0.5 * (b[1] * mid[2] - b[2] * mid[1]),
                0.5 * (b[2] * mid[0] - b[0] * mid[2]),
                0.5 * (b[0] * mid[1] - b[1] * mid[0]),
            ];
            let mut row = potential_row(average, WORLD, OrientationParity::Even);
            row.region = format!("original-edge-{index}");
            row
        })
        .collect()
}

fn edge_bindings<'a>(
    problem: &SpatialVectorMagnetostatics,
    rows: &'a [BoundaryCondition],
) -> Vec<MagneticPotentialBoundary3<'a>> {
    problem
        .exterior_faces()
        .into_iter()
        .map(|face| MagneticPotentialBoundary3::EdgeAverages {
            face,
            values: problem.face_edges(face).unwrap().map(|edge| &rows[edge]),
        })
        .collect()
}

#[test]
fn rotated_uniform_potential_uses_edge_integrals_without_creating_magnetic_flux() {
    let problem = vector_problem();
    let mut scenario = Scenario::new("rotating-axes-potential", 23, Environment::earth_lab());
    scenario.frames.add(Frame {
        id: FrameId(1),
        name: "rotated".into(),
        parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(
                Vec3::new(0.0, 0.0, 1.0),
                std::f64::consts::FRAC_PI_2,
            ),
            translation: Vec3::new(1e12, -1e12, 1e12),
        },
    });
    scenario.base_bcs = vec![potential_row(
        [1.0, -2.0, 3.0],
        FrameId(1),
        OrientationParity::Even,
    )];
    assert!(scenario.validate().is_empty());
    let decoded = parse_ir(&write_ir(&scenario)).unwrap();
    let restored = decoded.scenario();
    let faces = problem.exterior_faces();
    let rows = [MagneticPotentialBoundary3::Uniform {
        faces: &faces,
        row: &restored.base_bcs[0],
    }];
    for time in [0.0, 0.5, 1.0] {
        let result = solve_vector_potential3(
            &problem,
            &rows,
            &restored.frames,
            time,
            &[[0.0; 3]],
            vector_controls(),
            &mut || false,
        )
        .unwrap();
        let world = [2.0 * time, time, 3.0 * time];
        let mut seen = std::collections::BTreeMap::new();
        for trace in &result.boundary.potential {
            for (j, edge) in problem
                .face_edges(trace.face)
                .unwrap()
                .into_iter()
                .enumerate()
            {
                let [a, b] = problem.mesh().edges[edge];
                let p = problem.positions_m()[a as usize];
                let q = problem.positions_m()[b as usize];
                let expected: f64 = (0..3).map(|axis| world[axis] * (q[axis] - p[axis])).sum();
                near(trace.edge_integrals_wb[j], expected);
                if let Some(prior) = seen.insert(edge, trace.edge_integrals_wb[j]) {
                    assert_eq!(prior, trace.edge_integrals_wb[j]);
                }
            }
        }
        assert_eq!(result.time_s, time);
        for value in result.field.flux_density_t()[0] {
            near(value, 0.0);
        }
        near(result.field.stored_energy_j(), 0.0);
    }
}

#[test]
fn explicit_edge_average_histories_produce_nonzero_analytical_field() {
    let problem = vector_problem();
    let b = [0.3, -0.2, 0.4];
    let mut scenario = Scenario::new("edge-average-flux", 79, Environment::earth_lab());
    scenario.base_bcs = edge_history(&problem, b);
    let decoded = parse_ir(&write_ir(&scenario)).unwrap();
    let restored = decoded.scenario();
    let rows = edge_bindings(&problem, &restored.base_bcs);
    for time in [0.0, 0.25, 0.5, 1.0] {
        let result = solve_vector_potential3(
            &problem,
            &rows,
            &restored.frames,
            time,
            &[[0.0; 3]],
            vector_controls(),
            &mut || false,
        )
        .unwrap();
        for axis in 0..3 {
            near(result.field.flux_density_t()[0][axis], time * b[axis]);
            near(
                result.field.field_strength_a_m()[0][axis],
                time * b[axis] / 2.0,
            );
        }
        near(
            result.field.stored_energy_j(),
            time * time * b.iter().map(|v| v * v).sum::<f64>() / 24.0,
        );
        assert!(result.field.maximum_divergence_flux_wb <= vector_controls().flux_tolerance_wb);
    }
}

#[test]
fn vector_metadata_geometry_and_shared_edge_conflicts_refuse() {
    let problem = vector_problem();
    let frames = FrameTree::new();
    let faces = problem.exterior_faces();
    let valid = potential_row([0.0; 3], WORLD, OrientationParity::Even);
    let odd = potential_row([0.0; 3], WORLD, OrientationParity::Odd);
    let wrong_quantity = flux_row(0.0, OrientationParity::Odd);
    let mut wrong_frame = valid.clone();
    wrong_frame.frame = 99;
    for row in [&odd, &wrong_quantity, &wrong_frame] {
        let bindings = [MagneticPotentialBoundary3::Uniform { faces: &faces, row }];
        assert!(
            vector_potential_boundary_at(&problem, &bindings, &frames, 0.5, &mut || false).is_err()
        );
    }
    for selected in [&[faces[0], faces[0]][..], &[usize::MAX][..], &[][..]] {
        let bindings = [MagneticPotentialBoundary3::Uniform {
            faces: selected,
            row: &valid,
        }];
        assert!(
            vector_potential_boundary_at(&problem, &bindings, &frames, 0.5, &mut || false).is_err()
        );
    }
    let authored = edge_history(&problem, [0.3, -0.2, 0.4]);
    let mut bindings = edge_bindings(&problem, &authored);
    let inconsistent = potential_row([1.0, 0.0, 0.0], WORLD, OrientationParity::Even);
    if let MagneticPotentialBoundary3::EdgeAverages { values, .. } = &mut bindings[0] {
        values[0] = &inconsistent;
    }
    assert!(matches!(
        solve_vector_potential3(
            &problem,
            &bindings,
            &frames,
            0.5,
            &[[0.0; 3]],
            vector_controls(),
            &mut || false
        ),
        Err(ScenarioMagneticError::Native(_))
    ));
    let bindings = edge_bindings(&problem, &authored);
    for stop in [1, 20, 50] {
        let mut polls = 0;
        assert!(matches!(
            solve_vector_potential3(
                &problem,
                &bindings,
                &frames,
                0.5,
                &[[0.0; 3]],
                vector_controls(),
                &mut || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(ScenarioMagneticError::Cancelled)
        ));
    }
    let result = solve_vector_potential3(
        &problem,
        &bindings,
        &frames,
        0.5,
        &[[0.0; 3]],
        vector_controls(),
        &mut || false,
    )
    .unwrap();
    near(result.field.flux_density_t()[0][0], 0.15);
}

#[test]
fn interior_faces_cannot_become_magnetic_boundary_traces() {
    let problem = SpatialVectorMagnetostatics::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
        vec![[0, 1, 2, 3], [0, 2, 1, 4]],
        vec![SpatialMagneticMaterial::Linear(material(2.0)); 2],
        common::limits(),
        &mut || false,
    )
    .unwrap();
    let face = problem.mesh().faces.binary_search(&[0, 1, 2]).unwrap();
    assert_eq!(problem.exterior_normal(face), None);
    let row = potential_row([0.0; 3], WORLD, OrientationParity::Even);
    let bindings = [MagneticPotentialBoundary3::Uniform {
        faces: &[face],
        row: &row,
    }];
    assert!(
        vector_potential_boundary_at(&problem, &bindings, &FrameTree::new(), 0.5, &mut || false)
            .is_err()
    );
}
