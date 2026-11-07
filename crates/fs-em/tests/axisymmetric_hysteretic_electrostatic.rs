//! Original revolved storage, switching and source-work oracles.
//! Synthetic constant fits; these are numerical, not measured-material claims.
#![cfg(feature = "dielectric-hysteresis")]

use fs_em::electrostatic::nonlinear::{
    DielectricMaterial, ElectrodeControl, NonlinearElectrostaticControls,
    axisymmetric::{
        AxisymmetricNonlinearElectrostatics,
        hysteretic::{
            AxisymmetricHystereticElectrostatics, HystereticPathControls, HystereticPathError,
        },
    },
    hysteretic::{FrozenHystereticDielectric, PlanarHystereticElectrostatics},
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
use std::{f64::consts::TAU, sync::Arc};

const EPS: f64 = 1e-10;
const CHI: f64 = 2e-10;
const EC: f64 = 1e5;
const PS: f64 = 4e-5;
const HEIGHT: f64 = 1e-3;

fn near(actual: f64, expected: f64, scale: f64) {
    let tolerance = 2e-9 * scale.max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; allowance {tolerance:e}",
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
                    source: "authored cylindrical dielectric fixture".into(),
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
            process: "isotropic cylindrical play".into(),
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

fn local() -> VectorDielectricHysteresisControls {
    VectorDielectricHysteresisControls {
        relative_tolerance: 1e-12,
        absolute_field_tolerance_v_m: 1e-8,
        absolute_displacement_tolerance_c_m2: 1e-18,
        work_tolerance_j_m3: 1e-10,
        maximum_branch_iterations: 96,
        maximum_inverse_iterations: 96,
        maximum_backtracks: 32,
    }
}
fn controls() -> NonlinearElectrostaticControls {
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
fn path_controls() -> HystereticPathControls {
    HystereticPathControls {
        field: controls(),
        relative_tolerance: 0.0,
        absolute_polarization_tolerance_c_m2: 1e-8,
        absolute_heat_tolerance_j: 1e-10,
        absolute_work_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-13,
        maximum_field_solves: 8192,
        maximum_operator_applications: 400_000,
        maximum_refinement_depth: 20,
    }
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 32,
        max_dofs: 64,
        max_triplets: 288,
    }
}
fn annulus(inner: f64, outer: f64, metric: Metric2) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "axial-cylindrical-dielectric",
        vec![[inner, 0.0], [outer, 0.0], [outer, HEIGHT], [inner, HEIGHT]],
        vec![[0, 1, 2], [0, 2, 3]],
        metric,
    )
    .unwrap()
}
fn axial(
    inner: f64,
    outer: f64,
    angle: f64,
    switching: bool,
) -> AxisymmetricHystereticElectrostatics {
    AxisymmetricHystereticElectrostatics::new(
        annulus(inner, outer, Metric2::axisymmetric(angle).unwrap()),
        vec![source(switching); 2],
        &[&[0, 1], &[2, 3]],
        local(),
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn axial_major_and_minor_loops_use_original_revolved_volume_and_local_directions() {
    let (inner, outer) = (3e-3, 5e-3);
    for angle in [TAU, 0.7] {
        let area = 0.5 * angle * (outer * outer - inner * inner);
        let volume = area * HEIGHT;
        let mut field = axial(inner, outer, angle, true);
        let (mut old_p, mut old_d, mut old_w) = (0.0_f64, 0.0_f64, 0.0_f64);
        for voltage in [150.0, 500.0, 100.0, -350.0, 0.0] {
            let e = voltage / HEIGHT;
            // Independent collinear play oracle, including saturation and reversal.
            let p = old_p.max(CHI * (e - EC)).min(CHI * (e + EC)).clamp(-PS, PS);
            let d = EPS * e + p;
            let w = 0.5 * EPS * e * e + p * p / (2.0 * CHI);
            let heat = volume * EC * (p - old_p).abs();
            let work = volume * e * (d - old_d);
            let candidate = field
                .solve_candidate(
                    &[0.0; 2],
                    &[
                        ElectrodeControl::Voltage(voltage),
                        ElectrodeControl::Voltage(0.0),
                    ],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(candidate.field().terminal_charge_c()[0], area * d, 1e-9);
            near(candidate.field().terminal_charge_c()[1], -area * d, 1e-9);
            near(candidate.field().stored_energy_j(), volume * w, 1e-6);
            near(candidate.step().switching_heat_j, heat, 1e-6);
            near(candidate.step().stored_change_j, volume * (w - old_w), 1e-6);
            near(candidate.step().endpoint_source_work_j, work, 1e-6);
            near(
                candidate.step().numerical_dissipation_j,
                work - volume * (w - old_w) - heat,
                1e-6,
            );
            near(
                candidate
                    .step()
                    .cells
                    .iter()
                    .map(|cell| cell.volume_m3)
                    .sum(),
                volume,
                volume,
            );
            for cell in 0..2 {
                let history = candidate.field().cell_history(cell).unwrap();
                near(history.field_v_m()[0], 0.0, EC);
                near(history.field_v_m()[1], e, EC);
                assert_eq!(history.field_v_m()[2], 0.0);
                near(history.polarizations_c_m2()[0][0], 0.0, PS);
                near(history.polarizations_c_m2()[0][1], p, PS);
                assert_eq!(history.polarizations_c_m2()[0][2], 0.0);
            }
            field = candidate.into_parts().0;
            (old_p, old_d, old_w) = (p, d, w);
        }
    }
}

#[test]
fn changing_space_charge_pairs_exact_radial_moments_and_retained_field_work() {
    let (inner, outer, height, angle) = (2e-3, 6e-3, 2e-3, 0.8);
    let width = outer - inner;
    let center = inner + width / 4.0;
    let mesh = TriComplex2::from_indexed_triangles(
        "off-center-charged-cylindrical-fan",
        vec![
            [inner, 0.0],
            [outer, 0.0],
            [outer, height],
            [inner, height],
            [center, height / 2.0],
        ],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap();
    let mut field = AxisymmetricHystereticElectrostatics::new(
        mesh,
        vec![source(false); 4],
        &[&[0, 1, 2, 3]],
        local(),
        limits(),
        &mut || false,
    )
    .unwrap();
    // Exact one-free-node P1 oracle. The off-center radius makes the weak
    // source moment differ from total_volume/3, exposing an arithmetic-mean bug.
    let bottom_volume = angle * width * height * (inner + outer + center) / 12.0;
    let right_volume = angle * (outer - center) * height * (2.0 * outer + center) / 6.0;
    let left_volume = angle * (center - inner) * height * (2.0 * inner + center) / 6.0;
    let volume = angle * width * height * (inner + outer) / 2.0;
    let moment = angle * width * height * (1.5 * (inner + outer) + center) / 12.0;
    let stiffness = EPS
        * (8.0 * bottom_volume / (height * height)
            + right_volume / (outer - center).powi(2)
            + left_volume / (center - inner).powi(2));
    assert!((moment / volume - 1.0 / 3.0).abs() > 0.01);
    let (mut old_rho, mut old_peak) = (0.0, 0.0);
    for rho in [0.01, 0.02, -0.01] {
        let peak = rho * moment / stiffness;
        let energy = 0.5 * stiffness * peak * peak;
        let candidate = field
            .solve_candidate(
                &[rho; 4],
                &[ElectrodeControl::Voltage(0.0)],
                controls(),
                &mut || false,
            )
            .unwrap();
        near(candidate.field().relative_voltage_v()[4], peak, peak.abs());
        near(candidate.field().stored_energy_j(), energy, energy);
        near(
            candidate.field().terminal_charge_c()[0],
            -rho * volume,
            rho.abs() * volume,
        );
        near(
            candidate.step().endpoint_source_work_j,
            peak * (rho - old_rho) * moment,
            energy,
        );
        near(
            candidate.step().numerical_dissipation_j,
            0.5 * stiffness * (peak - old_peak).powi(2),
            energy,
        );
        assert_eq!(candidate.step().switching_heat_j, 0.0);
        near(
            candidate.field().cell_history(1).unwrap().field_v_m()[0],
            peak / (outer - center),
            EC,
        );
        near(
            candidate.field().cell_history(3).unwrap().field_v_m()[0],
            -peak / (center - inner),
            EC,
        );
        assert!(
            candidate.step().field_source_work_residual_j.abs() <= controls().energy_tolerance_j
        );
        field = candidate.into_parts().0;
        (old_rho, old_peak) = (rho, peak);
    }
}

#[test]
fn adaptive_cylindrical_paths_preserve_independent_gauges_and_cumulative_heat() {
    let (radius, angle) = (5e-3, 0.8);
    let area = 0.5 * angle * radius * radius;
    let volume = area * HEIGHT;
    let mesh = TriComplex2::from_indexed_triangles(
        "two-independent-axis-touching-cylinders",
        vec![
            [0.0, 0.0],
            [radius, 0.0],
            [radius, HEIGHT],
            [0.0, HEIGHT],
            [0.0, 2.0 * HEIGHT],
            [radius, 2.0 * HEIGHT],
            [radius, 3.0 * HEIGHT],
            [0.0, 3.0 * HEIGHT],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap();
    let initial = AxisymmetricHystereticElectrostatics::new(
        mesh,
        vec![source(true); 4],
        &[&[0, 1], &[2, 3], &[4, 5], &[6, 7]],
        local(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let e = 2.0 * EC;
    let p = CHI * (e - EC);
    let charge = area * (EPS * e + p);
    let boundary = |first, second| {
        [
            ElectrodeControl::Voltage(first),
            ElectrodeControl::Charge(-charge),
            ElectrodeControl::Voltage(second),
            ElectrodeControl::Charge(charge),
        ]
    };
    let mut budget = path_controls();
    budget.maximum_field_solves = 1;
    assert!(matches!(
        initial.solve_affine_path_candidate(&[0.0; 4], &boundary(0.0, 0.0), budget, &mut || false,),
        Err(HystereticPathError::FieldSolveBudget)
    ));
    assert!(matches!(
        initial.solve_affine_path_candidate(
            &[0.0; 4],
            &boundary(0.0, 0.0),
            path_controls(),
            &mut || true,
        ),
        Err(HystereticPathError::Cancelled)
    ));
    assert_eq!(initial.stored_energy_j(), 0.0);
    let ordinary = initial
        .solve_affine_path_candidate(&[0.0; 4], &boundary(0.0, 0.0), path_controls(), &mut || {
            false
        })
        .unwrap();
    let biased = initial
        .solve_affine_path_candidate(
            &[0.0; 4],
            &boundary(1e250, -1e250),
            path_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        ordinary.field().relative_voltage_v(),
        biased.field().relative_voltage_v()
    );
    assert_eq!(
        ordinary.field().terminal_charge_c(),
        biased.field().terminal_charge_c()
    );
    assert_eq!(ordinary.report().summary, biased.report().summary);
    assert_eq!(
        &biased.field().reference_voltage_by_vertex_v()[..4],
        &[1e250; 4]
    );
    assert_eq!(
        &biased.field().reference_voltage_by_vertex_v()[4..],
        &[-1e250; 4]
    );
    near(
        biased.field().stored_energy_j(),
        2.0 * volume * (0.5 * EPS * e * e + p * p / (2.0 * CHI)),
        1e-7,
    );
    near(
        biased.report().account.switching_heat_j,
        2.0 * volume * EC * p,
        1e-7,
    );
    assert!(biased.report().summary.accepted_substeps > 2);
    assert!(biased.report().summary.maximum_accepted_defect_ratio <= 1.0);
    near(
        biased.report().account.endpoint_source_work_j,
        biased.field().stored_energy_j()
            + biased.report().account.switching_heat_j
            + biased.report().account.numerical_dissipation_j,
        1e-7,
    );
    for cell in 0..4 {
        let sign = if cell < 2 { 1.0 } else { -1.0 };
        near(
            biased
                .field()
                .cell_history(cell)
                .unwrap()
                .polarizations_c_m2()[0][1],
            sign * p,
            PS,
        );
    }
    let (native, summary) = biased.into_electrostatic_candidate();
    assert_eq!(
        native.field().mesh().metric(),
        Metric2::axisymmetric(angle).unwrap()
    );
    assert!(summary.accepted_substeps > 2);
    // Switch the voltage reference independently in both components, retaining
    // exactly the same physical charge and polarization at the next endpoint.
    let (accepted, _) = native.into_parts();
    let continued = accepted
        .solve_candidate(
            &[0.0; 4],
            &[
                ElectrodeControl::Charge(charge),
                ElectrodeControl::Voltage(1e250),
                ElectrodeControl::Charge(-charge),
                ElectrodeControl::Voltage(-1e250),
            ],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(
        continued.field().stored_energy_j(),
        accepted.stored_energy_j(),
        1e-7,
    );
    near(continued.step().switching_heat_j, 0.0, 1e-7);
    assert_eq!(
        continued.field().reference_terminal_for_terminal(0),
        Some(1)
    );
    assert_eq!(
        continued.field().reference_terminal_for_terminal(2),
        Some(3)
    );
}

#[test]
fn retained_owners_preserve_geometry_guards_and_refused_history() {
    let radial = annulus(0.0, 5e-3, Metric2::axisymmetric(TAU).unwrap());
    let planar = annulus(0.0, 5e-3, Metric2::planar(1e-3).unwrap());
    assert!(
        PlanarHystereticElectrostatics::new(
            radial.clone(),
            vec![source(true); 2],
            &[&[0, 1], &[2, 3]],
            local(),
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        AxisymmetricHystereticElectrostatics::new(
            planar,
            vec![source(true); 2],
            &[&[0, 1], &[2, 3]],
            local(),
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        AxisymmetricHystereticElectrostatics::new(
            radial.clone(),
            vec![source(true); 2],
            &[&[0, 3], &[1, 2]],
            local(),
            limits(),
            &mut || false,
        )
        .is_err()
    );
    let frozen = FrozenHystereticDielectric::new(source(true), local(), &mut || false).unwrap();
    assert!(
        AxisymmetricNonlinearElectrostatics::new(
            radial,
            vec![DielectricMaterial::Hysteretic(frozen); 2],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    let initial = axial(0.0, 5e-3, TAU, true);
    assert!(
        initial
            .solve_candidate(
                &[0.0; 2],
                &[
                    ElectrodeControl::Voltage(1e9),
                    ElectrodeControl::Voltage(0.0)
                ],
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(initial.stored_energy_j(), 0.0);
    assert!(initial.relative_voltage_v().iter().all(|&v| v == 0.0));
    assert_eq!(
        initial.cell_history(0).unwrap().polarizations_c_m2(),
        &[[0.0; 3]]
    );
}
