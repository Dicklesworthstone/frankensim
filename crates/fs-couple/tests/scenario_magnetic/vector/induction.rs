//! G1/G3/G4: canonical A histories produce native Faraday heat and boundary work.
use super::*;
use fs_couple::scenario_magnetic::induction::{induction_boundary3, step_induction3};
use fs_em::spatial::vector_potential::induction::{
    SpatialInduction, SpatialInductionControls, SpatialInductionLimits,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY, ResolvedConductorState,
    },
    field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};

fn preparation() -> SpatialInductionLimits {
    SpatialInductionLimits {
        maximum_unknowns: 20000,
        maximum_conducting_cells: 4000,
        maximum_assembly_terms: 1_000_000,
    }
}
fn induction_controls() -> SpatialInductionControls {
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
fn conducting_problem() -> (SpatialVectorMagnetostatics, ResolvedConductorState) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            2.0,
        ),
        (
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            5.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured magnetic scenario".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
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
    let conductor = ResolvedConductorState::try_from_material(bundle.clone()).unwrap();
    let geometry = vector_problem();
    let problem = SpatialVectorMagnetostatics::new(
        geometry.positions_m().to_vec(),
        geometry.mesh().tets.clone(),
        vec![SpatialMagneticMaterial::Linear(bundle)],
        common::limits(),
        &mut || false,
    )
    .unwrap();
    (problem, conductor)
}

fn retime(
    rows: &[BoundaryCondition],
    times: &[f64],
    factors: &[f64],
    interpolation: TableInterpolation,
    outside: OutsideDomainPolicy,
) -> Vec<BoundaryCondition> {
    rows.iter()
        .map(|row| {
            let Some(BcValue::Typed(Payload::Vector(value))) = &row.value else {
                panic!("vector fixture")
            };
            let SampleSource::Table { values, .. } = value.source() else {
                panic!("table fixture")
            };
            let target = values.last().unwrap();
            let samples = factors
                .iter()
                .map(|factor| {
                    target
                        .iter()
                        .map(|q| QtyAny::new(factor * q.value, A_DIMS))
                        .collect()
                })
                .collect();
            let mut output = row.clone();
            output.value = Some(BcValue::Typed(Payload::Vector(
                VectorPayload::new(
                    value.meta().clone(),
                    SampleSource::table(
                        times.iter().map(|time| QtyAny::new(*time, TIME)).collect(),
                        samples,
                        interpolation,
                        outside,
                    )
                    .unwrap(),
                )
                .unwrap(),
            )));
            output
        })
        .collect()
}

#[test]
fn scenario_ramp_drives_analytical_faraday_current_heat_and_signed_boundary_work() {
    let (problem, conductor) = conducting_problem();
    // The authored one-second endpoint is Bz=2 T, so Bz=.2 T at .1 s.
    let authored = edge_history(&problem, [0.0, 0.0, 2.0]);
    let rows = edge_bindings(&problem, &authored);
    let frames = FrameTree::new();
    let current = [[0.0; 3]];
    let velocity = [[0.0; 3]; 4];
    let initial = solve_vector_potential3(
        &problem,
        &rows,
        &frames,
        0.0,
        &current,
        vector_controls(),
        &mut || false,
    )
    .unwrap();
    let state = SpatialInduction::new(
        &initial.field,
        vec![Some(conductor)],
        preparation(),
        1e-9,
        &mut || false,
    )
    .unwrap();
    let result = step_induction3(
        &state,
        &rows,
        &frames,
        &current,
        &velocity,
        0.1,
        induction_controls(),
        &mut || false,
    )
    .unwrap();
    let step = result.step();
    near(result.state().time_s(), 0.1);
    near(step.endpoint_flux_density_t[0][2], 0.2);
    near(step.heat_j, 1.0 / 160.0);
    near(step.stored_after_j, 0.2_f64.powi(2) / 24.0);
    near(
        step.magnetic_boundary_work_j,
        step.stored_after_j + step.heat_j,
    );
    near(step.energy_residual_j, 0.0);
    assert_eq!(step.impressed_work_j, 0.0);
    assert_eq!(step.mechanical_work_j, 0.0);
    for (local, &vertex) in problem.mesh().tets[0].iter().enumerate() {
        let p = problem.positions_m()[vertex as usize];
        let expected = [p[1] - 0.25, 0.25 - p[0], 0.0];
        for (axis, value) in expected.into_iter().enumerate() {
            near(
                step.cells[0].laboratory_electric_field_v_m[local][axis],
                value,
            );
            near(step.cells[0].current_density_a_m2[local][axis], 5.0 * value);
        }
    }
    let sampled =
        vector_potential_boundary_at(&problem, &rows, &frames, 0.1, &mut || false).unwrap();
    let endpoint = induction_boundary3(
        &state,
        &rows,
        &frames,
        0.1,
        induction_controls().flux_tolerance_wb,
        &mut || false,
    )
    .unwrap();
    assert_eq!(endpoint.potential.len(), sampled.potential.len());
    for (resolved, sampled) in endpoint.potential.iter().zip(&sampled.potential) {
        assert_eq!(resolved.face, sampled.face);
        assert_eq!(resolved.edge_integrals_wb, sampled.edge_integrals_wb);
    }
    let native = state
        .step_with_boundary_candidate(
            &current,
            &velocity,
            &endpoint,
            0.1,
            induction_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.state().potential_wb(), native.state().potential_wb());
    assert_eq!(step.heat_j.to_bits(), native.step().heat_j.to_bits());
    let half = step_induction3(
        &state,
        &rows,
        &frames,
        &current,
        &velocity,
        0.05,
        induction_controls(),
        &mut || false,
    )
    .unwrap();
    let second = step_induction3(
        half.state(),
        &rows,
        &frames,
        &current,
        &velocity,
        0.05,
        induction_controls(),
        &mut || false,
    )
    .unwrap();
    near(half.step().heat_j + second.step().heat_j, step.heat_j);
    near(second.step().stored_after_j, step.stored_after_j);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.field.potential_wb());
}

#[test]
fn incompatible_history_intervals_and_late_cancellation_do_not_advance_state() {
    let problem = vector_problem();
    let authored = edge_history(&problem, [0.0, 0.0, 0.2]);
    let rows = edge_bindings(&problem, &authored);
    let frames = FrameTree::new();
    let current = [[0.0; 3]];
    let velocity = [[0.0; 3]; 4];
    let initial = solve_vector_potential3(
        &problem,
        &rows,
        &frames,
        0.0,
        &current,
        vector_controls(),
        &mut || false,
    )
    .unwrap();
    let state = SpatialInduction::new(&initial.field, vec![None], preparation(), 1e-9, &mut || {
        false
    })
    .unwrap();
    for (times, factors, interpolation, outside) in [
        (
            vec![0.0, 1.0],
            vec![1.0, 2.0],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        ),
        (
            vec![0.0, 0.05, 1.0],
            vec![0.0, 0.5, 1.0],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        ),
        (
            vec![0.0, 0.1],
            vec![0.0, 1.0],
            TableInterpolation::StepLeft,
            OutsideDomainPolicy::Refuse,
        ),
        (
            vec![0.0, 1.0],
            vec![0.0, 1.0],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Periodic,
        ),
        (
            vec![0.05, 1.0],
            vec![0.0, 1.0],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Clamp,
        ),
    ] {
        let changed = retime(&authored, &times, &factors, interpolation, outside);
        let changed_rows = edge_bindings(&problem, &changed);
        assert!(matches!(
            step_induction3(
                &state,
                &changed_rows,
                &frames,
                &current,
                &velocity,
                0.1,
                induction_controls(),
                &mut || false
            ),
            Err(ScenarioMagneticError::Invalid(_))
        ));
    }
    for duration in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(
            step_induction3(
                &state,
                &rows,
                &frames,
                &current,
                &velocity,
                duration,
                induction_controls(),
                &mut || false
            )
            .is_err()
        );
    }
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            induction_boundary3(&state, &rows, &frames, 0.1, tolerance, &mut || false),
            Err(ScenarioMagneticError::Invalid(_))
        ));
    }
    assert!(matches!(
        step_induction3(
            &state,
            &rows[..1],
            &frames,
            &current,
            &velocity,
            0.1,
            induction_controls(),
            &mut || false
        ),
        Err(ScenarioMagneticError::Native(_))
    ));
    let mut calls = 0;
    step_induction3(
        &state,
        &rows,
        &frames,
        &current,
        &velocity,
        0.1,
        induction_controls(),
        &mut || {
            calls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, calls] {
        let mut polls = 0;
        assert!(matches!(
            step_induction3(
                &state,
                &rows,
                &frames,
                &current,
                &velocity,
                0.1,
                induction_controls(),
                &mut || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(ScenarioMagneticError::Cancelled)
        ));
    }
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.field.potential_wb());
    // Endpoint knots and the subsequent constant clamp interval remain usable.
    let clamped = retime(
        &authored,
        &[0.0, 0.1],
        &[0.0, 1.0],
        TableInterpolation::Linear,
        OutsideDomainPolicy::Clamp,
    );
    let clamped_rows = edge_bindings(&problem, &clamped);
    let ramp = step_induction3(
        &state,
        &clamped_rows,
        &frames,
        &current,
        &velocity,
        0.1,
        induction_controls(),
        &mut || false,
    )
    .unwrap();
    let hold = step_induction3(
        ramp.state(),
        &clamped_rows,
        &frames,
        &current,
        &velocity,
        0.1,
        induction_controls(),
        &mut || false,
    )
    .unwrap();
    near(hold.state().time_s(), 0.2);
    near(hold.step().magnetic_boundary_work_j, 0.0);
}

#[test]
fn only_the_selected_frame_ancestry_must_be_fixed() {
    let problem = vector_problem();
    let authored = edge_history(&problem, [0.0, 0.0, 0.2]);
    let rows = edge_bindings(&problem, &authored);
    let current = [[0.0; 3]];
    let velocity = [[0.0; 3]; 4];
    let mut frames = FrameTree::new();
    frames.add(Frame {
        id: FrameId(2),
        name: "rotating-parent".into(),
        parent: WORLD,
        motion: FrameMotion::Rotating {
            axis: [0.0, 0.0, 1.0],
            center: Vec3::new(0.0, 0.0, 0.0),
            rate: QtyAny::new(1.0, Dims([0, 0, -1, 0, 0, 0])),
        },
    });
    frames.add(Frame {
        id: FrameId(1),
        name: "fixed-child".into(),
        parent: FrameId(2),
        motion: FrameMotion::Fixed {
            orientation: Quat::identity(),
            translation: Vec3::new(0.0, 0.0, 0.0),
        },
    });
    let initial = solve_vector_potential3(
        &problem,
        &rows,
        &frames,
        0.0,
        &current,
        vector_controls(),
        &mut || false,
    )
    .unwrap();
    let state = SpatialInduction::new(&initial.field, vec![None], preparation(), 1e-9, &mut || {
        false
    })
    .unwrap();
    // A separate moving scenario frame cannot disable a WORLD-authored drive.
    step_induction3(
        &state,
        &rows,
        &frames,
        &current,
        &velocity,
        0.1,
        induction_controls(),
        &mut || false,
    )
    .unwrap();
    let moving: Vec<_> = authored
        .iter()
        .map(|row| {
            let Some(BcValue::Typed(Payload::Vector(value))) = &row.value else {
                panic!("vector fixture")
            };
            let meta = PayloadMeta::new(
                QuantityContract::Dimensions(A_DIMS),
                PayloadId::new("basis/frame-cartesian").unwrap(),
                FrameId(1),
                OrientationParity::Even,
                ReferenceSemantics::Continuous,
            )
            .unwrap();
            let mut row = row.clone();
            row.frame = 1;
            row.value = Some(BcValue::Typed(Payload::Vector(
                VectorPayload::new(meta, value.source().clone()).unwrap(),
            )));
            row
        })
        .collect();
    let moving_rows = edge_bindings(&problem, &moving);
    assert!(matches!(
        step_induction3(
            &state,
            &moving_rows,
            &frames,
            &current,
            &velocity,
            0.1,
            induction_controls(),
            &mut || false
        ),
        Err(ScenarioMagneticError::Invalid(_))
    ));
    assert_eq!(state.time_s(), 0.0);
}
