//! G0/G3/G4: resolved 3D induction, finite rotor inertia and transactional history.
//! Synthetic coefficients and a faceted cylinder are numerical fixtures, not
//! validation against a measured machine or a continuum braking coefficient.

#[path = "spatial_induction_rotor/thermal.rs"]
mod thermal;

use fs_couple::{
    induction_rotor::RotorControls,
    spatial_induction::{
        CylindricalConductor, SpatialInductionRotor, SpatialRotorCandidate, SpatialRotorError,
    },
};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        induction::{
            SpatialInduction, SpatialInductionControls, SpatialInductionLimits,
            SpatialInductionStep, mechanical::SpatialLorentzLoadLimits,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_ga::{Mat3, Se3, So3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_mbd::articulated::{ArticulatedModel, JointModel, Link, SpatialInertia};
use fs_qty::{QuantitySpec, Temperature};
use fs_solver::CheckedCgConfig;

const SIGMA: f64 = 6.0;
const MU: f64 = 0.2;
const INERTIA: f64 = 0.2;
const HALF_LENGTH: f64 = 0.3;
const DT: f64 = 0.02;

fn source() -> ResolvedMaterialStatePoint {
    let basis = ElectricalPropertyBasis::Conductivity;
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (basis.property_name(), basis.dimensions(), SIGMA),
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            MU,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    1000.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic spatial induction rotor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "spatial-rotor-fixture".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap(),
        &[
            basis.requirement().unwrap(),
            ScalarPropertyRequirement::try_new(
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn six_volume(p: [[f64; 3]; 4]) -> f64 {
    dot(
        difference(p[1], p[0]),
        cross(difference(p[2], p[0]), difference(p[3], p[0])),
    )
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a:.16e} != {b:.16e}, error {:.6e}",
        (a - b).abs()
    );
}

fn problem(
    sectors: usize,
) -> (
    SpatialVectorMagnetostatics,
    Vec<Option<ResolvedConductorState>>,
) {
    let mut section = vec![[0.0; 2]];
    for radius in [0.5, 1.0] {
        for sector in 0..sectors {
            let angle = std::f64::consts::TAU * sector as f64 / sectors as f64;
            section.push([radius * angle.cos(), radius * angle.sin()]);
        }
    }
    let mut triangles = Vec::new();
    for i in 0..sectors {
        let a = 1 + i;
        let b = 1 + (i + 1) % sectors;
        triangles.extend([
            [0, a, b],
            [a, a + sectors, b + sectors],
            [a, b + sectors, b],
        ]);
    }
    let positions: Vec<_> = [-HALF_LENGTH, 0.0, HALF_LENGTH]
        .into_iter()
        .flat_map(|z| section.iter().map(move |p| [p[0], p[1], z]))
        .collect();
    let mut tets = Vec::new();
    for layer in 0..2 {
        for triangle in &triangles {
            // Global ordering chooses the same diagonal on each shared prism
            // face. The middle layer/ring supplies genuine free interior dofs.
            let mut sorted = *triangle;
            sorted.sort_unstable();
            let [a, b, c] = sorted.map(|v| (v + layer * section.len()) as u32);
            let [aa, bb, cc] = [a, b, c].map(|v| v + section.len() as u32);
            for mut tet in [[a, b, c, cc], [a, b, bb, cc], [a, aa, bb, cc]] {
                if six_volume(tet.map(|v| positions[v as usize])) < 0.0 {
                    tet.swap(0, 1);
                }
                assert!(six_volume(tet.map(|v| positions[v as usize])) > 0.0);
                tets.push(tet);
            }
        }
    }
    let bundle = source();
    let conductor = ResolvedConductorState::try_from_material(bundle.clone()).unwrap();
    let count = tets.len();
    (
        SpatialVectorMagnetostatics::new(
            positions,
            tets,
            vec![SpatialMagneticMaterial::Linear(bundle); count],
            WeightedAssemblyLimits {
                max_cells: 512,
                max_dofs: 4096,
                max_triplets: 100000,
            },
            &mut || false,
        )
        .unwrap(),
        vec![Some(conductor); count],
    )
}

fn field_controls() -> SpatialInductionControls {
    SpatialInductionControls {
        relative_tolerance: 1e-11,
        absolute_residual_a: 1e-11,
        restart: 64,
        maximum_cycles: 100,
        maximum_krylov_scalars: 2_000_000,
        current_tolerance_a: 2e-9,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 2e-9,
        maximum_cell_magnetic_peclet: 30.0,
    }
}
fn controls() -> RotorControls {
    RotorControls {
        maximum_iterations: 64,
        relaxation: 0.7,
        angular_velocity_tolerance_rad_s: 1e-10,
        energy_tolerance_j: 5e-9,
        angular_impulse_tolerance_nms: 1e-10,
        linear_impulse_tolerance_ns: 1e-10,
    }
}
fn limits() -> SpatialLorentzLoadLimits {
    SpatialLorentzLoadLimits {
        max_vertices: 256,
        max_cells: 512,
    }
}
fn shape() -> CylindricalConductor {
    CylindricalConductor {
        inner_radius_m: 0.0,
        outer_radius_m: 1.0,
        half_length_m: HALF_LENGTH,
        maximum_boundary_error_m: 0.025,
    }
}
fn model(axis_x: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "native spatial electromagnetic rotor",
        None,
        Se3::from_parts(So3::identity(), Vec3::new(axis_x, 0.0, 0.0)).unwrap(),
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
        SpatialInertia::new(
            2.0,
            Vec3::default(),
            Mat3 {
                m: [INERTIA, 0.0, 0.0, 0.0, INERTIA, 0.0, 0.0, 0.0, INERTIA],
            },
        )
        .unwrap(),
    )])
    .unwrap()
}
fn uniform_trace(problem: &SpatialVectorMagnetostatics, b: [f64; 3]) -> VectorMagneticBoundary {
    // Exact Whitney integrals of A = (B cross r)/2 on the original edges.
    let a: Vec<_> = problem
        .mesh()
        .edges
        .iter()
        .map(|&[u, v]| {
            let p = problem.positions_m()[u as usize];
            let q = problem.positions_m()[v as usize];
            let mid = std::array::from_fn(|i| 0.5 * (p[i] + q[i]));
            dot(cross(b, mid).map(|v| 0.5 * v), difference(q, p))
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
fn initial<'m>(
    problem: &'m SpatialVectorMagnetostatics,
    sources: Vec<Option<ResolvedConductorState>>,
) -> SpatialInduction<'m> {
    // The initial exterior A trace imposes transverse B=(1,0,0); all original
    // interior field dofs remain free during the subsequent induction solve.
    let boundary = uniform_trace(problem, [1.0, 0.0, 0.0]);
    let cg = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 3000,
        max_corrections: 2,
    };
    let field = problem
        .solve(
            &vec![[0.0; 3]; problem.mesh().tets.len()],
            &boundary,
            VectorMagneticControls {
                field: cg,
                gauge: cg,
                maximum_gauge_product_terms: 2_000_000,
                current_tolerance_a: 2e-9,
                flux_tolerance_wb: 2e-9,
                energy_tolerance_j: 2e-9,
            },
            &mut || false,
        )
        .unwrap();
    SpatialInduction::new(
        &field,
        sources,
        SpatialInductionLimits {
            maximum_unknowns: 4096,
            maximum_conducting_cells: 512,
            maximum_assembly_terms: 1_000_000,
        },
        2e-9,
        &mut || false,
    )
    .unwrap()
}

// Independent polynomial quadrature: integral(lambda_i lambda_j) is V/20
// off-diagonal and V/10 on-diagonal. Using V/4 at the four force samples
// would give the wrong torque and mechanical work even if net force agreed.
fn integrate(
    problem: &SpatialVectorMagnetostatics,
    step: &SpatialInductionStep,
    dt: f64,
) -> ([f64; 3], [f64; 3], f64) {
    let (mut force, mut torque, mut heat) = ([0.0; 3], [0.0; 3], 0.0);
    for cell in &step.cells {
        let p = problem.mesh().tets[cell.cell].map(|v| problem.positions_m()[v as usize]);
        let volume = six_volume(p) / 6.0;
        let j = cell.current_density_a_m2;
        let f = j.map(|current| cross(current, step.stage_flux_density_t[cell.cell]));
        for i in 0..4 {
            for axis in 0..3 {
                force[axis] += volume * f[i][axis] / 4.0;
            }
            for k in 0..4 {
                let weight = volume * if i == k { 0.1 } else { 0.05 };
                let moment = cross(p[i], f[k]);
                for axis in 0..3 {
                    torque[axis] += weight * moment[axis];
                }
                heat += dt * weight * dot(j[i], j[k]) / SIGMA;
            }
        }
    }
    (force, torque, heat)
}
fn stored(problem: &SpatialVectorMagnetostatics, b: &[[f64; 3]]) -> f64 {
    problem
        .mesh()
        .tets
        .iter()
        .zip(b)
        .map(|(tet, &flux)| {
            let volume = six_volume(tet.map(|v| problem.positions_m()[v as usize])) / 6.0;
            volume * dot(flux, flux) / (2.0 * MU)
        })
        .sum()
}
fn assert_balance(
    original: &SpatialInductionRotor<'_>,
    candidate: &SpatialRotorCandidate<'_>,
    dt: f64,
    shaft_torque: f64,
) {
    let before = original.angular_velocity_rad_s();
    let after = candidate.state().angular_velocity_rad_s();
    let midpoint = 0.5 * (before + after);
    let (force, torque, heat) = integrate(original.field().problem(), &candidate.induction, dt);
    let b0 = original.field().flux_density_t(&mut || false).unwrap();
    let u0 = stored(original.field().problem(), &b0);
    let u1 = stored(
        original.field().problem(),
        &candidate.induction.endpoint_flux_density_t,
    );
    let kinetic = 0.5 * INERTIA * (after - before) * (after + before);
    near(
        candidate.balance.field_angular_velocity_rad_s,
        midpoint,
        1e-10,
    );
    near(
        candidate.balance.mechanical_midpoint_velocity_rad_s,
        midpoint,
        1e-12,
    );
    near(
        candidate.state().angle_rad(),
        original.angle_rad() + dt * midpoint,
        1e-12,
    );
    near(
        candidate.balance.electromagnetic_torque_nm,
        torque[2],
        1e-10,
    );
    near(candidate.induction.heat_j, heat, 1e-10);
    near(candidate.induction.stored_before_j, u0, 2e-9);
    near(candidate.induction.stored_after_j, u1, 2e-9);
    near(
        candidate.induction.mechanical_work_j,
        dt * midpoint * torque[2],
        2e-9,
    );
    near(candidate.mechanics.kinetic_energy_change, kinetic, 1e-12);
    near(
        INERTIA * (after - before),
        dt * (torque[2] + shaft_torque),
        1e-10,
    );
    near(
        candidate.balance.shaft_work_j,
        dt * midpoint * shaft_torque,
        1e-12,
    );
    near(
        kinetic + u1 - u0 + heat,
        dt * midpoint * shaft_torque + candidate.induction.magnetic_boundary_work_j,
        2e-8,
    );
    near(candidate.balance.total_energy_residual_j, 0.0, 5e-9);
    near(candidate.mechanics.support_work, 0.0, 1e-15);
    near(candidate.induction.impressed_work_j, 0.0, 1e-15);
    near(candidate.induction.terminal_work_j, 0.0, 1e-15);
    let body = candidate.mechanics.external_impulse_world.to_array();
    let support = candidate.mechanics.support_impulse_world.to_array();
    for axis in 0..3 {
        near(body[axis], dt * torque[axis], 1e-10);
        near(body[axis + 3], dt * force[axis], 1e-10);
        near(body[axis + 3] + support[axis + 3], 0.0, 1e-10);
        near(
            body[axis] + support[axis],
            if axis == 2 {
                INERTIA * (after - before)
            } else {
                0.0
            },
            1e-10,
        );
    }
}

#[test]
fn applied_rotating_trace_drives_native_rotor_and_retains_endpoint() {
    use fs_couple::scenario_magnetic::{
        MagneticPotentialBoundary3,
        induction::{ScenarioMagneticRotorError, step_rotor3},
    };
    use fs_qty::{Dims, QtyAny};
    use fs_scenario::payload::{
        OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta, QuantityContract,
        ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
    };
    use fs_scenario::{BcKind, BcValue, BoundaryCondition, FrameTree, Physics, WORLD};

    let (problem, sources) = problem(16);
    let original = SpatialInductionRotor::new(
        initial(&problem, sources),
        model(0.0),
        shape(),
        0.0,
        0.0,
        limits(),
        &mut || false,
    )
    .unwrap();
    let original_a = original.field().potential_wb().to_vec();
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let boundary = uniform_trace(&problem, [1.0, 0.1, 0.0]);
    let mut checkpoints = 0;
    let candidate = original
        .step_with_boundary_candidate(
            DT,
            &zero,
            &boundary,
            0.0,
            field_controls(),
            controls(),
            &mut || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    assert!(candidate.state().angular_velocity_rad_s() > 0.0);
    assert!(candidate.balance.electromagnetic_torque_nm > 0.0);
    assert!(candidate.induction.heat_j > 0.0);
    assert!(candidate.induction.magnetic_boundary_work_j > 0.0);
    assert_balance(&original, &candidate, DT, 0.0);

    // Authored line-average A histories drive this same original rotor. Their
    // initial values match the retained transverse field; no field is reset.
    let dims = Dims([1, 1, -2, 0, -1, 0]);
    let authored: Vec<_> = problem
        .mesh()
        .edges
        .iter()
        .enumerate()
        .map(|(edge, &[a, b])| {
            let p = problem.positions_m()[a as usize];
            let q = problem.positions_m()[b as usize];
            let mid = std::array::from_fn(|axis| 0.5 * (p[axis] + q[axis]));
            let values = [[1.0, 0.0, 0.0], [1.0, 0.1, 0.0]]
                .map(|b| cross(b, mid).map(|a| QtyAny::new(0.5 * a, dims)).to_vec())
                .to_vec();
            let meta = PayloadMeta::new(
                QuantityContract::Dimensions(dims),
                PayloadId::new("basis/frame-cartesian").unwrap(),
                WORLD,
                OrientationParity::Even,
                ReferenceSemantics::Continuous,
            )
            .unwrap();
            BoundaryCondition {
                region: format!("rotor-original-edge-{edge}"),
                physics: Physics::Magnetics,
                kind: BcKind::MagneticVectorPotential,
                value: Some(BcValue::Typed(Payload::Vector(
                    VectorPayload::new(
                        meta,
                        SampleSource::table(
                            [0.0, DT]
                                .map(|time| QtyAny::new(time, Dims([0, 0, 1, 0, 0, 0])))
                                .to_vec(),
                            values,
                            TableInterpolation::Linear,
                            OutsideDomainPolicy::Refuse,
                        )
                        .unwrap(),
                    )
                    .unwrap(),
                ))),
                compatibility: None,
                frame: WORLD.0,
            }
        })
        .collect();
    let rows: Vec<_> = problem
        .exterior_faces()
        .into_iter()
        .map(|face| MagneticPotentialBoundary3::EdgeAverages {
            face,
            values: problem
                .face_edges(face)
                .unwrap()
                .map(|edge| &authored[edge]),
        })
        .collect();
    let frames = FrameTree::new();
    let scenario = step_rotor3(
        &original,
        &rows,
        &frames,
        &zero,
        0.0,
        DT,
        field_controls(),
        controls(),
        &mut || false,
    )
    .unwrap();
    near(
        scenario.state().angular_velocity_rad_s(),
        candidate.state().angular_velocity_rad_s(),
        2e-9,
    );
    near(scenario.induction.heat_j, candidate.induction.heat_j, 2e-9);
    near(
        scenario.induction.magnetic_boundary_work_j,
        candidate.induction.magnetic_boundary_work_j,
        2e-9,
    );
    for (&a, &b) in scenario
        .state()
        .field()
        .potential_wb()
        .iter()
        .zip(candidate.state().field().potential_wb())
    {
        near(a, b, 2e-9);
    }
    assert!(matches!(
        step_rotor3(
            &original,
            &rows,
            &frames,
            &[],
            0.0,
            DT,
            field_controls(),
            controls(),
            &mut || false
        ),
        Err(ScenarioMagneticRotorError::NativeRotor(
            SpatialRotorError::Invalid(_)
        ))
    ));
    assert!(matches!(
        step_rotor3(
            &original,
            &rows,
            &frames,
            &zero,
            0.0,
            0.0,
            field_controls(),
            controls(),
            &mut || false
        ),
        Err(ScenarioMagneticRotorError::Scenario(_))
    ));

    // Recompute from the original accepted field at the converged speed.
    // This independently checks that speed trials did not advance the trace.
    let midpoint = 0.5 * candidate.state().angular_velocity_rad_s();
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| [-midpoint * p[1], midpoint * p[0], 0.0])
        .collect();
    let direct = original
        .field()
        .step_with_boundary_candidate(
            &zero,
            &velocity,
            &boundary,
            DT,
            field_controls(),
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in candidate
        .state()
        .field()
        .potential_wb()
        .iter()
        .zip(direct.state().potential_wb())
    {
        near(a, b, 2e-9);
    }
    near(candidate.induction.heat_j, direct.step().heat_j, 2e-9);
    near(
        candidate.induction.magnetic_boundary_work_j,
        direct.step().magnetic_boundary_work_j,
        2e-9,
    );

    let mut polled = 0;
    let refused = original.step_with_boundary_candidate(
        DT,
        &zero,
        &boundary,
        0.0,
        field_controls(),
        controls(),
        &mut || {
            polled += 1;
            polled == checkpoints
        },
    );
    assert!(matches!(refused, Err(SpatialRotorError::Cancelled)));
    let incomplete = VectorMagneticBoundary { potential: vec![] };
    assert!(
        original
            .step_with_boundary_candidate(
                DT,
                &zero,
                &incomplete,
                0.0,
                field_controls(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(original.field().potential_wb(), original_a);
    assert_eq!(original.field().time_s(), 0.0);
    assert_eq!(original.angular_velocity_rad_s(), 0.0);

    let accepted = candidate.into_state();
    let hold = accepted
        .step_candidate(
            2.0 * DT,
            &zero,
            0.0,
            field_controls(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(hold.induction.magnetic_boundary_work_j, 0.0);
    assert_balance(&accepted, &hold, DT, 0.0);
    for trace in &boundary.potential {
        for (edge, expected) in problem
            .face_edges(trace.face)
            .unwrap()
            .into_iter()
            .zip(trace.edge_integrals_wb)
        {
            assert_eq!(hold.state().field().potential_wb()[edge], expected);
        }
    }
}

#[test]
fn transverse_braking_uses_finite_inertia_and_the_same_accepted_field_stage() {
    let (problem, sources) = problem(16);
    let field = initial(&problem, sources);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    for speed in [2.0, -2.0] {
        let original = SpatialInductionRotor::new(
            field.clone(),
            model(0.0),
            shape(),
            0.37,
            speed,
            limits(),
            &mut || false,
        )
        .unwrap();
        near(
            original.boundary_error_m(),
            1.0 - (std::f64::consts::PI / 16.0).cos(),
            1e-12,
        );
        let potential = original.field().potential_wb().to_vec();
        let candidate = original
            .step_candidate(DT, &zero, 0.0, field_controls(), controls(), &mut || false)
            .unwrap();
        let after = candidate.state().angular_velocity_rad_s();
        assert!(candidate.balance.iterations > 1);
        assert!(after * speed > 0.0 && after.abs() < speed.abs());
        assert!(candidate.balance.electromagnetic_torque_nm * speed < 0.0);
        assert!(candidate.induction.heat_j > 0.0);
        assert!(candidate.induction.stored_after_j > candidate.induction.stored_before_j);
        assert_balance(&original, &candidate, DT, 0.0);

        // Re-solving from the ORIGINAL accepted field at the mechanical
        // midpoint must reproduce the coupled update. Repeatedly advancing
        // field history during mechanical iterations would fail this check.
        let midpoint = 0.5 * (speed + after);
        let velocity: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| [-midpoint * p[1], midpoint * p[0], 0.0])
            .collect();
        let direct = original
            .field()
            .step_candidate(&zero, &velocity, DT, field_controls(), &mut || false)
            .unwrap();
        for (&coupled, &native) in candidate
            .state()
            .field()
            .potential_wb()
            .iter()
            .zip(direct.state().potential_wb())
        {
            near(coupled, native, 2e-9);
        }
        near(candidate.induction.heat_j, direct.step().heat_j, 2e-9);
        near(
            candidate.induction.mechanical_work_j,
            direct.step().mechanical_work_j,
            2e-9,
        );
        assert!(
            candidate
                .state()
                .field()
                .potential_wb()
                .iter()
                .zip(&potential)
                .any(|(a, b)| (a - b).abs() > 1e-8)
        );
        assert_eq!(original.field().potential_wb(), potential);
        assert_eq!(original.field().time_s(), 0.0);
        assert_eq!(original.angular_velocity_rad_s(), speed);
    }
}

#[test]
fn shaft_driven_cylinder_retains_field_history_and_full_bearing_reaction() {
    let (problem, sources) = problem(16);
    let original = SpatialInductionRotor::new(
        initial(&problem, sources),
        model(0.0),
        shape(),
        -0.2,
        0.0,
        limits(),
        &mut || false,
    )
    .unwrap();
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let driven = original
        .step_candidate(DT, &zero, 1.0, field_controls(), controls(), &mut || false)
        .unwrap();
    assert!(driven.state().angular_velocity_rad_s() > 0.0);
    assert!(driven.state().angular_velocity_rad_s() < DT / INERTIA);
    assert!(driven.induction.heat_j > 0.0 && driven.balance.shaft_work_j > 0.0);
    near(driven.mechanics.support_impulse_world.torque.z, DT, 1e-10);
    assert_balance(&original, &driven, DT, 1.0);
    let accepted = driven.into_state();
    let potential = accepted.field().potential_wb().to_vec();
    let coast = accepted
        .step_candidate(
            2.0 * DT,
            &zero,
            0.0,
            field_controls(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(coast.induction.heat_j > 0.0);
    assert_balance(&accepted, &coast, DT, 0.0);
    assert_eq!(coast.state().field().time_s(), 2.0 * DT);
    assert_eq!(accepted.field().potential_wb(), potential);
    assert_eq!(accepted.field().time_s(), DT);
    assert_eq!(original.field().time_s(), 0.0);
    assert_eq!(original.angular_velocity_rad_s(), 0.0);
}

#[test]
fn geometry_refusal_and_failed_feedback_publish_neither_native_history() {
    let (square, sources) = problem(4);
    assert!(matches!(
        SpatialInductionRotor::new(
            initial(&square, sources),
            model(0.0),
            shape(),
            0.0,
            2.0,
            limits(),
            &mut || false
        ),
        Err(SpatialRotorError::Invalid(_))
    ));
    let (problem, sources) = problem(16);
    let field = initial(&problem, sources);
    assert!(matches!(
        SpatialInductionRotor::new(
            field.clone(),
            model(0.1),
            shape(),
            0.0,
            2.0,
            limits(),
            &mut || false
        ),
        Err(SpatialRotorError::Invalid(_))
    ));
    let original =
        SpatialInductionRotor::new(field, model(0.0), shape(), 0.25, 2.0, limits(), &mut || {
            false
        })
        .unwrap();
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let potential = original.field().potential_wb().to_vec();
    assert!(matches!(
        original.step_candidate(
            DT,
            &zero,
            0.0,
            field_controls(),
            RotorControls {
                maximum_iterations: 1,
                ..controls()
            },
            &mut || false
        ),
        Err(SpatialRotorError::Nonconvergence { iterations: 1, .. })
    ));
    let mut internal_source = zero.clone();
    internal_source[0] = [1.0, 0.0, 0.0];
    assert!(matches!(
        original.step_candidate(
            DT,
            &internal_source,
            0.0,
            field_controls(),
            controls(),
            &mut || false
        ),
        Err(SpatialRotorError::Invalid(_))
    ));

    let mut successful_polls = 0;
    original
        .step_candidate(DT, &zero, 0.0, field_controls(), controls(), &mut || {
            successful_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    let cancelled =
        original.step_candidate(DT, &zero, 0.0, field_controls(), controls(), &mut || {
            polls += 1;
            polls == successful_polls
        });
    assert!(matches!(cancelled, Err(SpatialRotorError::Cancelled)));
    assert_eq!(polls, successful_polls);
    assert_eq!(original.field().potential_wb(), potential);
    assert_eq!(original.field().time_s(), 0.0);
    assert_eq!(original.angular_velocity_rad_s(), 2.0);
    assert_eq!(original.angle_rad(), 0.25);
}
