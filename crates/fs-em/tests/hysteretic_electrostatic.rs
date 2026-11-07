//! Native Gauss/electrode solves with retained polarization and independent work.
//! Synthetic SI fixtures, not measured ferroelectric qualification.
#![cfg(feature = "dielectric-hysteresis")]

use fs_em::electrostatic::nonlinear::{
    ElectrodeControl, NonlinearElectrostaticControls, hysteretic::PlanarHystereticElectrostatics,
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
use fs_qty::{QuantitySpec, Temperature};
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
fn near(actual: f64, expected: f64, scale: f64) {
    let tolerance = 2e-9 * scale.max(expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance={tolerance:e}"
    );
}
fn source(epsilon: f64, switching: bool) -> ResolvedVectorDielectricHysteresis {
    source_at_temperature(epsilon, switching, None)
}
fn temperature_query(temperature_k: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature_k,
        )
        .unwrap()
}
fn source_at_temperature(
    epsilon: f64,
    switching: bool,
    temperature_k: Option<f64>,
) -> ResolvedVectorDielectricHysteresis {
    let query = temperature_k.map(temperature_query).unwrap_or_default();
    let validity = if temperature_k.is_some() {
        ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            280.0,
            360.0,
        )
    } else {
        ValidityDomain::unconstrained()
    };
    let branch = DielectricHysteresisBranchSources {
        susceptibility: PropertyKey::new("switching_slope", ELECTRIC_PERMITTIVITY_DIMS),
        pinning_field: PropertyKey::new("pinning", ELECTRIC_FIELD_DIMS),
        saturation_polarization: PropertyKey::new("polarization_limit", ELECTRIC_DISPLACEMENT_DIMS),
    };
    let roles = DielectricHysteresisSources {
        background_permittivity: PropertyKey::new("epsilon_background", ELECTRIC_PERMITTIVITY_DIMS),
        field_limit: PropertyKey::new("field_limit", ELECTRIC_FIELD_DIMS),
        branches: if switching {
            vec![branch.clone()]
        } else {
            vec![]
        },
    };
    let mut entries = vec![
        (roles.background_permittivity.clone(), epsilon),
        (roles.field_limit.clone(), 2e6),
    ];
    if switching {
        entries.extend([
            (branch.susceptibility, CHI),
            (branch.pinning_field, EC),
            (branch.saturation_polarization, PS),
        ]);
    }
    let mut claims = ClaimSet::new();
    for (key, value) in entries {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored spatial switching fixture".into(),
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
            process: "analytic spatial play fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedVectorDielectricHysteresis::from_isotropic_source(
        ResolvedDielectricHysteresis::resolve(
            Arc::new(card),
            &query,
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
        max_newton_iterations: 48,
        max_backtracks: 32,
        max_operator_applications: 4096,
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 128,
            max_corrections: 2,
        },
    }
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 32,
        max_dofs: 64,
        max_triplets: 288,
    }
}
fn mesh(rotation: [f64; 2]) -> TriComplex2 {
    let [c, s] = rotation;
    let points = [
        [0.0, 0.0],
        [GAP / 2.0, 0.0],
        [GAP, 0.0],
        [0.0, WIDTH],
        [GAP / 2.0, WIDTH],
        [GAP, WIDTH],
    ];
    TriComplex2::from_indexed_triangles(
        "retained-polarization-slab",
        points.map(|[x, y]| [c * x - s * y, s * x + c * y]).to_vec(),
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap()
}
fn slab(switching: bool, rotation: [f64; 2]) -> PlanarHystereticElectrostatics {
    PlanarHystereticElectrostatics::new(
        mesh(rotation),
        vec![source(EPS, switching); 4],
        &[&[0, 3], &[2, 5]],
        local(),
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn voltage(value: f64) -> [ElectrodeControl; 2] {
    [
        ElectrodeControl::Voltage(value),
        ElectrodeControl::Voltage(0.0),
    ]
}

#[test]
fn disconnected_histories_keep_independent_biases_when_reference_electrodes_change() {
    let single = mesh([1.0, 0.0]);
    let mut vertices = single.vertices().to_vec();
    vertices.extend(single.vertices().iter().map(|&[x, y]| [x, y + 2.0 * WIDTH]));
    let triangles = vec![
        [0, 1, 4],
        [0, 4, 3],
        [1, 2, 5],
        [1, 5, 4],
        [6, 7, 10],
        [6, 10, 9],
        [7, 8, 11],
        [7, 11, 10],
    ];
    let independent = TriComplex2::from_indexed_triangles(
        "independently-biased-switching-slabs",
        vertices,
        triangles,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    let mut ordinary = PlanarHystereticElectrostatics::new(
        independent,
        vec![source(EPS, true); 8],
        &[&[0, 3], &[2, 5], &[6, 9], &[8, 11]],
        local(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut biased = ordinary.clone();
    for (step, displacement) in [2.5e-5, -3.5e-5, 1.5e-5, -2.0e-5].into_iter().enumerate() {
        let q = AREA * displacement;
        let switched = step % 2 == 1;
        let boundary = |first, second| {
            if switched {
                [
                    ElectrodeControl::Charge(-q),
                    ElectrodeControl::Voltage(first),
                    ElectrodeControl::Charge(0.75 * q),
                    ElectrodeControl::Voltage(second),
                ]
            } else {
                [
                    ElectrodeControl::Voltage(first),
                    ElectrodeControl::Charge(q),
                    ElectrodeControl::Voltage(second),
                    ElectrodeControl::Charge(-0.75 * q),
                ]
            }
        };
        let a = ordinary
            .solve_candidate(&[0.0; 8], &boundary(0.0, 0.0), controls(), &mut || false)
            .unwrap();
        let b = biased
            .solve_candidate(&[0.0; 8], &boundary(1e250, -1e250), controls(), &mut || {
                false
            })
            .unwrap();
        assert_eq!(
            a.field().relative_voltage_v(),
            b.field().relative_voltage_v()
        );
        assert_eq!(a.field().terminal_charge_c(), b.field().terminal_charge_c());
        assert_eq!(a.field().stored_energy_j(), b.field().stored_energy_j());
        assert_eq!(a.step().switching_heat_j, b.step().switching_heat_j);
        assert_eq!(
            a.step().endpoint_source_work_j,
            b.step().endpoint_source_work_j
        );
        assert_eq!(
            a.step().numerical_dissipation_j,
            b.step().numerical_dissipation_j
        );
        assert!(b.step().maximum_component_charge_residual_c <= controls().charge_tolerance_c);
        let first_reference = if switched { 1 } else { 0 };
        assert_eq!(
            &b.field().reference_terminal_by_vertex()[..6],
            &[first_reference; 6]
        );
        assert_eq!(
            &b.field().reference_terminal_by_vertex()[6..],
            &[first_reference + 2; 6]
        );
        assert_eq!(&b.field().reference_voltage_by_vertex_v()[..6], &[1e250; 6]);
        assert_eq!(
            &b.field().reference_voltage_by_vertex_v()[6..],
            &[-1e250; 6]
        );
        for cell in 0..8 {
            assert_eq!(
                a.field().cell_history(cell).unwrap().polarizations_c_m2(),
                b.field().cell_history(cell).unwrap().polarizations_c_m2(),
            );
        }
        ordinary = a.into_parts().0;
        biased = b.into_parts().0;
    }
}

#[test]
fn linear_voltage_charge_and_common_bias_limits_close_independent_work() {
    let initial = slab(false, [1.0, 0.0]);
    let v = 400.0;
    let c = AREA * EPS / GAP;
    let candidate = initial
        .solve_candidate(&[0.0; 4], &voltage(v), controls(), &mut || false)
        .unwrap();
    let field = candidate.field();
    near(field.terminal_charge_c()[0], c * v, 1e-9);
    near(field.terminal_charge_c()[1], -c * v, 1e-9);
    near(field.stored_energy_j(), 0.5 * c * v * v, 1e-6);
    near(candidate.step().endpoint_source_work_j, c * v * v, 1e-6);
    near(
        candidate.step().numerical_dissipation_j,
        0.5 * c * v * v,
        1e-6,
    );
    assert_eq!(candidate.step().switching_heat_j, 0.0);
    assert!(candidate.step().nonlinear.newton_iterations > 0);
    let charged = initial
        .solve_candidate(
            &[0.0; 4],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(-c * v),
            ],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(charged.field().terminal_voltage_v()[1], -v, v);
    near(
        charged.field().stored_energy_j(),
        field.stored_energy_j(),
        1e-6,
    );
    let shifted = initial
        .solve_candidate(
            &[0.0; 4],
            &[
                ElectrodeControl::Voltage(1e12 + v),
                ElectrodeControl::Voltage(1e12),
            ],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        shifted.field().relative_voltage_v(),
        field.relative_voltage_v()
    );
    assert_eq!(
        shifted.field().terminal_charge_c(),
        field.terminal_charge_c()
    );
    let held = field
        .solve_candidate(&[0.0; 4], &voltage(v), controls(), &mut || false)
        .unwrap();
    near(held.step().stored_change_j, 0.0, 1e-6);
    near(held.step().endpoint_source_work_j, 0.0, 1e-6);
    assert_eq!(initial.stored_energy_j(), 0.0);
}

#[test]
fn spatial_major_and_minor_histories_match_independent_scalar_work_channels() {
    let fit = source(EPS, true);
    let scalar = fit.scalar_source().law();
    let mut memory = scalar.initial_state();
    let mut field = slab(true, [1.0, 0.0]);
    for e in [5e5, 0.0, 2e5, -5e5, -1e5, 0.0, 5e5] {
        let expected = scalar.candidate_field(&memory, e).unwrap();
        let old_d = memory.displacement_c_m2();
        let old = field.cell_history(0).unwrap().clone();
        let candidate = field
            .solve_candidate(&[0.0; 4], &voltage(e * GAP), controls(), &mut || false)
            .unwrap();
        assert_eq!(field.cell_history(0).unwrap(), &old);
        for cell in 0..4 {
            let h = candidate.field().cell_history(cell).unwrap();
            near(h.field_v_m()[0], e, 5e5);
            near(h.field_v_m()[1], 0.0, 5e5);
            near(
                h.displacement_c_m2()[0],
                expected.state.displacement_c_m2(),
                1e-4,
            );
            near(
                h.polarizations_c_m2()[0][0],
                expected.state.polarizations_c_m2()[0],
                PS,
            );
        }
        let step = candidate.step();
        near(
            step.switching_heat_j,
            VOLUME * expected.irreversible_heat_j_m3,
            1e-6,
        );
        near(
            step.stored_change_j,
            VOLUME * expected.energy_change_j_m3,
            1e-6,
        );
        let endpoint = VOLUME * e * (expected.state.displacement_c_m2() - old_d);
        near(step.endpoint_source_work_j, endpoint, 1e-6);
        near(
            step.numerical_dissipation_j,
            endpoint - VOLUME * (expected.energy_change_j_m3 + expected.irreversible_heat_j_m3),
            1e-6,
        );
        assert!(step.numerical_dissipation_j >= 0.0);
        assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
        memory = expected.state;
        field = candidate.into_parts().0;
    }
}

#[test]
fn floating_conductor_solves_two_switching_media_with_one_total_charge() {
    let split = TriComplex2::from_indexed_triangles(
        "switching-floating-electrode",
        vec![
            [0.0, 0.0],
            [GAP, 0.0],
            [GAP, WIDTH],
            [0.0, WIDTH],
            [GAP, 0.0],
            [2.0 * GAP, 0.0],
            [2.0 * GAP, WIDTH],
            [GAP, WIDTH],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    let model = PlanarHystereticElectrostatics::new(
        split,
        vec![
            source(EPS, true),
            source(EPS, true),
            source(EPS / 2.0, true),
            source(EPS / 2.0, true),
        ],
        &[&[0, 3], &[1, 2, 4, 7], &[5, 6]],
        local(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let v = 400.0;
    let result = model
        .solve_candidate(
            &[0.0; 4],
            &[
                ElectrodeControl::Voltage(v),
                ElectrodeControl::Charge(0.0),
                ElectrodeControl::Voltage(0.0),
            ],
            controls(),
            &mut || false,
        )
        .unwrap();
    // Both branches slide without saturation: D=(epsilon+chi)E-chi Ec.
    let k = (v / GAP) / (1.0 / (EPS + CHI) + 1.0 / (EPS / 2.0 + CHI));
    let e1 = k / (EPS + CHI);
    let e2 = k / (EPS / 2.0 + CHI);
    let d = k - CHI * EC;
    near(result.field().terminal_charge_c()[0], AREA * d, 1e-9);
    near(result.field().terminal_charge_c()[1], 0.0, 1e-9);
    near(result.field().terminal_charge_c()[2], -AREA * d, 1e-9);
    near(result.field().terminal_voltage_v()[1], -GAP * e1, v);
    let p1 = CHI * (e1 - EC);
    let p2 = CHI * (e2 - EC);
    let energy = VOLUME
        * (0.5 * EPS * e1 * e1
            + p1 * p1 / (2.0 * CHI)
            + 0.25 * EPS * e2 * e2
            + p2 * p2 / (2.0 * CHI));
    near(result.field().stored_energy_j(), energy, 1e-6);
    near(
        result.step().switching_heat_j,
        VOLUME * EC * (p1 + p2),
        1e-6,
    );
    near(result.step().endpoint_source_work_j, v * AREA * d, 1e-6);
    assert!(result.step().maximum_floating_charge_residual_c < controls().charge_tolerance_c);
}

#[test]
fn rotated_geometry_preserves_vector_memory_charge_and_losses() {
    let mut ordinary = slab(true, [1.0, 0.0]);
    let mut rotated = slab(true, [0.6, 0.8]);
    for v in [150.0, 500.0, 0.0, -300.0] {
        let a = ordinary
            .solve_candidate(&[0.0; 4], &voltage(v), controls(), &mut || false)
            .unwrap();
        let b = rotated
            .solve_candidate(&[0.0; 4], &voltage(v), controls(), &mut || false)
            .unwrap();
        near(a.step().switching_heat_j, b.step().switching_heat_j, 1e-6);
        near(
            a.step().numerical_dissipation_j,
            b.step().numerical_dissipation_j,
            1e-6,
        );
        near(
            a.field().terminal_charge_c()[0],
            b.field().terminal_charge_c()[0],
            1e-9,
        );
        for cell in 0..4 {
            let pa = a.field().cell_history(cell).unwrap().polarizations_c_m2()[0];
            let pb = b.field().cell_history(cell).unwrap().polarizations_c_m2()[0];
            near(pb[0], 0.6 * pa[0] - 0.8 * pa[1], PS);
            near(pb[1], 0.8 * pa[0] + 0.6 * pa[1], PS);
        }
        ordinary = a.into_parts().0;
        rotated = b.into_parts().0;
    }
}

#[test]
fn changing_volume_charge_supplies_its_own_work_without_invented_switching_heat() {
    let initial = slab(false, [1.0, 0.0]);
    let rho = 0.02;
    let expected_peak = rho * GAP * GAP / (8.0 * EPS);
    let energy = rho * rho * AREA * GAP * GAP * GAP / (32.0 * EPS);
    let first = initial
        .solve_candidate(&[rho; 4], &voltage(0.0), controls(), &mut || false)
        .unwrap();
    near(
        first.field().relative_voltage_v()[1],
        expected_peak,
        expected_peak,
    );
    near(
        first.field().relative_voltage_v()[4],
        expected_peak,
        expected_peak,
    );
    near(first.field().stored_energy_j(), energy, energy);
    near(first.step().endpoint_source_work_j, 2.0 * energy, energy);
    near(first.step().numerical_dissipation_j, energy, energy);
    assert_eq!(first.step().switching_heat_j, 0.0);
    let next = first
        .field()
        .solve_candidate(&[2.0 * rho; 4], &voltage(0.0), controls(), &mut || false)
        .unwrap();
    near(next.step().stored_change_j, 3.0 * energy, energy);
    near(next.step().endpoint_source_work_j, 4.0 * energy, energy);
    near(next.step().numerical_dissipation_j, energy, energy);
}

#[test]
fn cancellation_budgets_and_field_domain_leave_the_original_history_replayable() {
    let model = slab(true, [1.0, 0.0])
        .solve_candidate(&[0.0; 4], &voltage(200.0), controls(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let memory: Vec<_> = (0..4)
        .map(|i| model.cell_history(i).unwrap().clone())
        .collect();
    let charges = model.terminal_charge_c().to_vec();
    let mut polls = 0;
    assert!(
        model
            .solve_candidate(&[0.0; 4], &voltage(-300.0), controls(), &mut || {
                polls += 1;
                polls > 8
            })
            .is_err()
    );
    let mut limited = controls();
    limited.max_operator_applications = 0;
    assert!(
        model
            .solve_candidate(&[0.0; 4], &voltage(-300.0), limited, &mut || false)
            .is_err()
    );
    assert!(
        model
            .solve_candidate(&[0.0; 4], &voltage(5000.0), controls(), &mut || false)
            .is_err()
    );
    for (i, state) in memory.iter().enumerate() {
        assert_eq!(model.cell_history(i).unwrap(), state);
    }
    assert_eq!(model.terminal_charge_c(), charges);
    let a = model
        .solve_candidate(&[0.0; 4], &voltage(-300.0), controls(), &mut || false)
        .unwrap();
    let b = model
        .solve_candidate(&[0.0; 4], &voltage(-300.0), controls(), &mut || false)
        .unwrap();
    assert_eq!(a.field().terminal_charge_c(), b.field().terminal_charge_c());
    assert_eq!(
        a.field().relative_voltage_v(),
        b.field().relative_voltage_v()
    );
    assert_eq!(a.step().cells, b.step().cells);
}

fn thermal_slab() -> PlanarHystereticElectrostatics {
    PlanarHystereticElectrostatics::new(
        mesh([1.0, 0.0]),
        vec![source_at_temperature(EPS, true, Some(300.0)); 4],
        &[&[0, 3], &[2, 5]],
        local(),
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn source_temperature_refresh_preserves_poled_field_and_pending_increment() {
    let field = thermal_slab()
        .solve_candidate(&[0.0; 4], &voltage(500.0), controls(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let queries = [310.0, 315.0, 320.0, 325.0].map(temperature_query);
    let refreshed = field.at_source_queries(&queries, &mut || false).unwrap();
    assert!(std::ptr::eq(
        field.problem().mesh(),
        refreshed.problem().mesh()
    ));
    assert_eq!(field.relative_voltage_v(), refreshed.relative_voltage_v());
    assert_eq!(field.terminal_charge_c(), refreshed.terminal_charge_c());
    assert_eq!(
        field.stored_energy_j().to_bits(),
        refreshed.stored_energy_j().to_bits()
    );
    for (i, query) in queries.iter().enumerate() {
        assert_eq!(field.cell_history(i), refreshed.cell_history(i));
        assert_eq!(
            refreshed.problem().materials()[i].material().query(),
            *query
        );
        assert_eq!(
            field.problem().materials()[i].material().query(),
            temperature_query(300.0)
        );
    }
    // The same constant fit must yield the same continued reversal trajectory.
    let a = field
        .solve_candidate(&[0.0; 4], &voltage(-300.0), controls(), &mut || false)
        .unwrap();
    let b = refreshed
        .solve_candidate(&[0.0; 4], &voltage(-300.0), controls(), &mut || false)
        .unwrap();
    assert_eq!(a.step().cells, b.step().cells);
    assert_eq!(a.field().terminal_charge_c(), b.field().terminal_charge_c());
    let expected_memory: Vec<_> = (0..4)
        .map(|i| a.field().cell_history(i).unwrap().clone())
        .collect();
    let expected_work = a.step().endpoint_source_work_j.to_bits();
    let rebound = a.with_source_queries(&queries, &mut || false).unwrap();
    assert_eq!(
        rebound.step().endpoint_source_work_j.to_bits(),
        expected_work
    );
    for (i, memory) in expected_memory.iter().enumerate() {
        assert_eq!(rebound.field().cell_history(i).unwrap(), memory);
        assert_eq!(
            rebound.field().problem().materials()[i].material().query(),
            queries[i]
        );
    }
}

#[test]
fn incomplete_cancelled_and_out_of_domain_source_refreshes_are_atomic() {
    let field = thermal_slab()
        .solve_candidate(&[0.0; 4], &voltage(500.0), controls(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let before: Vec<_> = (0..4)
        .map(|i| field.cell_history(i).unwrap().clone())
        .collect();
    let queries = [310.0, 315.0, 320.0, 325.0].map(temperature_query);
    assert!(
        field
            .at_source_queries(&queries[..3], &mut || false)
            .is_err()
    );
    let mut invalid = queries.clone();
    invalid[3] = temperature_query(361.0);
    assert!(field.at_source_queries(&invalid, &mut || false).is_err());
    let mut polls = 0;
    assert!(
        field
            .at_source_queries(&queries, &mut || {
                polls += 1;
                polls > 3
            })
            .is_err()
    );
    assert!(
        field
            .solve_candidate(&[0.0; 4], &voltage(-300.0), controls(), &mut || false)
            .unwrap()
            .with_source_queries(&invalid, &mut || false)
            .is_err()
    );
    for (i, memory) in before.iter().enumerate() {
        assert_eq!(field.cell_history(i).unwrap(), memory);
        assert_eq!(
            field.problem().materials()[i].material().query(),
            temperature_query(300.0)
        );
    }
    let a = field.at_source_queries(&queries, &mut || false).unwrap();
    let b = field.at_source_queries(&queries, &mut || false).unwrap();
    assert_eq!(a.terminal_charge_c(), b.terminal_charge_c());
    for i in 0..4 {
        assert_eq!(a.cell_history(i), b.cell_history(i));
    }
}
