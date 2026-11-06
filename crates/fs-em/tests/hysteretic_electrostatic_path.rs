//! G1/G3/G4: load-path refinement on the original Gauss and vector-play owners.
//! Synthetic SI fixtures; the error estimates are not continuum certificates.
#![cfg(feature = "dielectric-hysteresis")]

use fs_em::electrostatic::nonlinear::{
    ElectrodeControl, NonlinearElectrostaticControls,
    hysteretic::{
        PlanarHystereticElectrostatics,
        path::{HystereticPathControls, HystereticPathError},
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{
            hysteresis::{
                source::{
                    DielectricHysteresisBranchSources, DielectricHysteresisSources,
                    ResolvedDielectricHysteresis,
                },
                vector::{ResolvedVectorDielectricHysteresis, VectorDielectricHysteresisControls},
            },
            reversible::{ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS},
        },
    },
    state_point::MaterialPropertySelection,
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;
use std::sync::Arc;

const GAP: f64 = 1e-3;
const WIDTH: f64 = 2e-2;
const DEPTH: f64 = 5e-3;
const AREA: f64 = WIDTH * DEPTH;
const VOLUME: f64 = AREA * GAP;
const EPS: f64 = 1e-10;
const CHI: f64 = 2e-10;
const EC: f64 = 1e5;
const PS: f64 = 4e-5;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}
fn source(switching: bool) -> ResolvedVectorDielectricHysteresis {
    let branch = DielectricHysteresisBranchSources {
        susceptibility: PropertyKey::new("chi", ELECTRIC_PERMITTIVITY_DIMS),
        pinning_field: PropertyKey::new("Ec", ELECTRIC_FIELD_DIMS),
        saturation_polarization: PropertyKey::new("Ps", ELECTRIC_DISPLACEMENT_DIMS),
    };
    let roles = DielectricHysteresisSources {
        background_permittivity: PropertyKey::new("epsilon", ELECTRIC_PERMITTIVITY_DIMS),
        field_limit: PropertyKey::new("field_limit", ELECTRIC_FIELD_DIMS),
        branches: if switching {
            vec![branch.clone()]
        } else {
            vec![]
        },
    };
    let mut values = vec![
        (roles.background_permittivity.clone(), EPS),
        (roles.field_limit.clone(), 2e6),
    ];
    if switching {
        values.extend([
            (branch.susceptibility, CHI),
            (branch.pinning_field, EC),
            (branch.saturation_polarization, PS),
        ]);
    }
    let mut claims = ClaimSet::new();
    for (key, value) in values {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "analytical electric loading path".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "quadratic isotropic play".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedVectorDielectricHysteresis::from_isotropic_source(
        ResolvedDielectricHysteresis::resolve(
            Arc::new(card),
            &QueryPoint::new(),
            roles,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
    .unwrap()
}
fn native() -> NonlinearElectrostaticControls {
    NonlinearElectrostaticControls {
        absolute_residual_tolerance_c: 1e-22,
        relative_residual_tolerance: 1e-12,
        charge_tolerance_c: 1e-19,
        energy_tolerance_j: 1e-15,
        max_newton_iterations: 64,
        max_backtracks: 40,
        max_operator_applications: 4096,
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 128,
            max_corrections: 2,
        },
    }
}
fn controls() -> HystereticPathControls {
    HystereticPathControls {
        field: native(),
        relative_tolerance: 0.0,
        absolute_polarization_tolerance_c_m2: 1e-8,
        absolute_heat_tolerance_j: 1e-9,
        absolute_work_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-13,
        maximum_field_solves: 8192,
        maximum_operator_applications: 400_000,
        maximum_refinement_depth: 20,
    }
}
fn field(
    mesh: TriComplex2,
    terminals: &[&[usize]],
    switching: bool,
) -> PlanarHystereticElectrostatics {
    let count = mesh.faces().len();
    PlanarHystereticElectrostatics::new(
        mesh,
        vec![source(switching); count],
        terminals,
        VectorDielectricHysteresisControls {
            relative_tolerance: 1e-12,
            absolute_field_tolerance_v_m: 1e-8,
            absolute_displacement_tolerance_c_m2: 1e-18,
            work_tolerance_j_m3: 1e-10,
            maximum_branch_iterations: 96,
            maximum_inverse_iterations: 96,
            maximum_backtracks: 32,
        },
        WeightedAssemblyLimits {
            max_cells: 32,
            max_dofs: 64,
            max_triplets: 288,
        },
        &mut || false,
    )
    .unwrap()
}
fn slab(switching: bool) -> PlanarHystereticElectrostatics {
    let mesh = TriComplex2::from_indexed_triangles(
        "affine-slab",
        vec![[0.0, 0.0], [GAP, 0.0], [GAP, WIDTH], [0.0, WIDTH]],
        vec![[0, 1, 2], [0, 2, 3]],
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    field(mesh, &[&[0, 3], &[1, 2]], switching)
}
fn hexagon() -> PlanarHystereticElectrostatics {
    let height = GAP * 0.8660254037844386;
    let mesh = TriComplex2::from_indexed_triangles(
        "changing-vector-path",
        vec![
            [GAP, 0.0],
            [GAP / 2.0, height],
            [-GAP / 2.0, height],
            [-GAP, 0.0],
            [-GAP / 2.0, -height],
            [GAP / 2.0, -height],
            [0.0, 0.0],
        ],
        vec![
            [6, 0, 1],
            [6, 1, 2],
            [6, 2, 3],
            [6, 3, 4],
            [6, 4, 5],
            [6, 5, 0],
        ],
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    field(mesh, &[&[0, 1], &[2, 3], &[4, 5]], true)
}
fn voltage(v: f64) -> [ElectrodeControl; 2] {
    [
        ElectrodeControl::Voltage(0.0),
        ElectrodeControl::Voltage(-v),
    ]
}
fn vector_voltage(a: f64, b: f64) -> [ElectrodeControl; 3] {
    [
        ElectrodeControl::Voltage(0.0),
        ElectrodeControl::Voltage(a),
        ElectrodeControl::Voltage(b),
    ]
}
fn memories(field: &PlanarHystereticElectrostatics) -> Vec<Vec<[f64; 3]>> {
    (0..field.problem().materials().len())
        .map(|i| field.cell_history(i).unwrap().polarizations_c_m2().to_vec())
        .collect()
}
fn memory_error(a: &PlanarHystereticElectrostatics, b: &PlanarHystereticElectrostatics) -> f64 {
    memories(a)
        .iter()
        .zip(memories(b))
        .flat_map(|(a, b)| {
            a.iter()
                .zip(b)
                .map(|(a, b)| (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]))
                .collect::<Vec<_>>()
        })
        .fold(0.0, f64::max)
}

#[test]
fn linear_source_work_converges_without_invented_heat_and_preserves_bias_and_charge_control() {
    let initial = slab(false);
    let v = 500.0;
    let capacitance = EPS * AREA / GAP;
    let stored = 0.5 * capacitance * v * v;
    let mut loose = controls();
    loose.absolute_work_tolerance_j = stored;
    let a = initial
        .solve_affine_path_candidate(&[0.0; 2], &voltage(v), loose, &mut || false)
        .unwrap();
    assert_eq!(a.report().summary.field_solves, 3);
    assert_eq!(a.report().summary.accepted_substeps, 2);
    near(
        a.report().account.numerical_dissipation_j,
        stored / 2.0,
        1e-18,
    );
    let mut tight = controls();
    tight.absolute_work_tolerance_j = stored / 100.0;
    let b = initial
        .solve_affine_path_candidate(&[0.0; 2], &voltage(v), tight, &mut || false)
        .unwrap();
    let report = b.report();
    near(report.account.stored_change_j, stored, 1e-18);
    assert_eq!(report.account.switching_heat_j, 0.0);
    assert!(report.account.numerical_dissipation_j < tight.absolute_work_tolerance_j);
    assert!(
        report.account.numerical_dissipation_j < a.report().account.numerical_dissipation_j / 32.0
    );
    near(
        report.account.endpoint_source_work_j,
        stored + report.account.numerical_dissipation_j,
        1e-18,
    );
    near(
        report.summary.estimated_source_work_error_j,
        report.account.numerical_dissipation_j,
        1e-18,
    );
    assert!(report.summary.maximum_accepted_defect_ratio <= 1.0);
    let biased = initial
        .solve_affine_path_candidate(
            &[0.0; 2],
            &[
                ElectrodeControl::Voltage(1e12),
                ElectrodeControl::Voltage(1e12 - v),
            ],
            tight,
            &mut || false,
        )
        .unwrap();
    assert_eq!(biased.field().reference_voltage_v(), 1e12);
    assert_eq!(
        biased.field().relative_voltage_v(),
        b.field().relative_voltage_v()
    );
    assert_eq!(
        biased.field().terminal_charge_c(),
        b.field().terminal_charge_c()
    );
    assert_eq!(biased.report().summary, b.report().summary);
    let charge = initial
        .solve_affine_path_candidate(
            &[0.0; 2],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(-capacitance * v),
            ],
            tight,
            &mut || false,
        )
        .unwrap();
    near(charge.field().terminal_voltage_v()[1], -v, 1e-8);
    near(
        charge.report().account.endpoint_source_work_j,
        report.account.endpoint_source_work_j,
        1e-16,
    );
    let expected_cells = b.report().account.cells.clone();
    let expected_summary = b.report().summary;
    let (native_candidate, summary) = b.into_electrostatic_candidate();
    assert_eq!(summary, expected_summary);
    assert_eq!(native_candidate.step().cells, expected_cells);
    assert_eq!(initial.stored_energy_j(), 0.0);
}

#[test]
fn major_loop_retains_exact_switching_heat_and_separate_refined_work() {
    let mut model = slab(true)
        .solve_candidate(&[0.0; 2], &voltage(500.0), native(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let initial_memory = memories(&model);
    let initial_energy = model.stored_energy_j();
    let expected_heat = 4.0 * EC * PS * VOLUME;
    let mut c = controls();
    c.absolute_work_tolerance_j = expected_heat / 200.0;
    let (mut heat, mut work, mut numerical, mut stored) = (0.0, 0.0, 0.0, 0.0);
    for v in [-500.0, 500.0] {
        let candidate = model
            .solve_affine_path_candidate(&[0.0; 2], &voltage(v), c, &mut || false)
            .unwrap();
        let report = candidate.report();
        heat += report.account.switching_heat_j;
        work += report.account.endpoint_source_work_j;
        numerical += report.account.numerical_dissipation_j;
        stored += report.account.stored_change_j;
        near(
            report
                .account
                .cells
                .iter()
                .map(|cell| cell.switching_heat_j)
                .sum(),
            report.account.switching_heat_j,
            1e-18,
        );
        model = candidate.into_parts().0;
    }
    near(heat, expected_heat, 1e-17);
    near(stored, 0.0, 1e-17);
    near(model.stored_energy_j(), initial_energy, 1e-17);
    near(work, expected_heat + numerical, 1e-17);
    assert!(numerical > 0.0 && numerical < expected_heat / 50.0);
    assert_eq!(memories(&model), initial_memory);
}

#[test]
fn changing_vector_history_converges_to_an_independently_subdivided_loading_path() {
    let model = hexagon()
        .solve_candidate(
            &[0.0; 6],
            &vector_voltage(400.0, -150.0),
            native(),
            &mut || false,
        )
        .unwrap()
        .into_parts()
        .0;
    let target = vector_voltage(-150.0, 400.0);
    let mut reference = model.clone();
    let (mut reference_heat, mut reference_work) = (0.0, 0.0);
    // Independently authored uniform partition, no adaptive controller or
    // full/half estimator is involved in this reference trajectory.
    for k in 1..=2048 {
        let s = k as f64 / 2048.0;
        let candidate = reference
            .solve_candidate(
                &[0.0; 6],
                &vector_voltage(400.0 - 550.0 * s, -150.0 + 550.0 * s),
                native(),
                &mut || false,
            )
            .unwrap();
        reference_heat += candidate.step().switching_heat_j;
        reference_work += candidate.step().endpoint_source_work_j;
        reference = candidate.into_parts().0;
    }
    let mut loose = controls();
    loose.absolute_polarization_tolerance_c_m2 = 5e-7;
    loose.absolute_heat_tolerance_j = 2e-9;
    loose.absolute_work_tolerance_j = 2e-9;
    let coarse = model
        .solve_affine_path_candidate(&[0.0; 6], &target, loose, &mut || false)
        .unwrap();
    let mut tight = loose;
    tight.absolute_polarization_tolerance_c_m2 /= 8.0;
    tight.absolute_heat_tolerance_j /= 8.0;
    tight.absolute_work_tolerance_j /= 8.0;
    let fine = model
        .solve_affine_path_candidate(&[0.0; 6], &target, tight, &mut || false)
        .unwrap();
    let old = memories(&model);
    let new = memories(fine.field());
    assert!(
        old.iter()
            .zip(&new)
            .any(|(a, b)| (a[0][0] * b[0][1] - a[0][1] * b[0][0]).abs() > 1e-11)
    );
    assert!(fine.report().summary.field_solves > coarse.report().summary.field_solves);
    assert!(memory_error(fine.field(), &reference) < memory_error(coarse.field(), &reference));
    assert!(
        (fine.report().account.switching_heat_j - reference_heat).abs()
            < (coarse.report().account.switching_heat_j - reference_heat).abs()
    );
    assert!(
        (fine.report().account.endpoint_source_work_j - reference_work).abs()
            < (coarse.report().account.endpoint_source_work_j - reference_work).abs()
    );
    assert!(
        fine.report().summary.estimated_polarization_error_c_m2
            <= tight.absolute_polarization_tolerance_c_m2
    );
    assert!(
        fine.report().summary.estimated_switching_heat_error_j <= tight.absolute_heat_tolerance_j
    );
    assert!(fine.report().summary.estimated_source_work_error_j <= tight.absolute_work_tolerance_j);
    // A whole-path operator cap includes coarse and discarded trials.
    let mut bounded = tight;
    bounded.maximum_operator_applications = fine.report().summary.operator_applications - 1;
    assert!(matches!(
        model.solve_affine_path_candidate(&[0.0; 6], &target, bounded, &mut || false),
        Err(HystereticPathError::OperatorBudget)
    ));
}

#[test]
fn budgets_refinement_and_late_cancellation_leave_the_original_history_replayable() {
    let model = slab(true)
        .solve_candidate(&[0.0; 2], &voltage(200.0), native(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let original_memory = memories(&model);
    let original_charges = model.terminal_charge_c().to_vec();
    let mut polls = 0;
    let complete = model
        .solve_affine_path_candidate(&[0.0; 2], &voltage(-500.0), controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 10);
    let mut count = 0;
    assert!(matches!(
        model.solve_affine_path_candidate(&[0.0; 2], &voltage(-500.0), controls(), &mut || {
            count += 1;
            count == polls
        }),
        Err(HystereticPathError::Cancelled)
    ));
    let mut work = controls();
    work.maximum_field_solves = complete.report().summary.field_solves - 1;
    assert!(matches!(
        model.solve_affine_path_candidate(&[0.0; 2], &voltage(-500.0), work, &mut || false),
        Err(HystereticPathError::FieldSolveBudget)
    ));
    let mut shallow = controls();
    shallow.maximum_refinement_depth = 0;
    assert!(matches!(
        model.solve_affine_path_candidate(&[0.0; 2], &voltage(-500.0), shallow, &mut || false),
        Err(HystereticPathError::RefinementLimit { .. })
    ));
    assert!(
        model
            .solve_affine_path_candidate(&[0.0; 2], &voltage(5000.0), controls(), &mut || false)
            .is_err()
    );
    assert_eq!(memories(&model), original_memory);
    assert_eq!(model.terminal_charge_c(), original_charges);
    let replay = model
        .solve_affine_path_candidate(&[0.0; 2], &voltage(-500.0), controls(), &mut || false)
        .unwrap();
    assert_eq!(replay.report().summary, complete.report().summary);
    assert_eq!(
        replay.report().account.cells,
        complete.report().account.cells
    );
    assert_eq!(memories(replay.field()), memories(complete.field()));
}
