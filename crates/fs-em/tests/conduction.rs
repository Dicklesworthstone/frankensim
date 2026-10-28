//! G1/G3/G4: spatial source-defined Ohmic current, native circuit loss and heat.
//! These are synthetic physical-law checks, not measured material qualification.
mod common;
#[path = "conduction/response.rs"]
mod response;

use common::{config, ends, layered, limits, near};
use fs_em::conduction::{
    ConductionControls, ConductionCurrentTrace, ConductionError, PlanarConduction,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ElectricalPropertyBasis, ResolvedConductorState, resolve_conductor_state_with_basis,
    },
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_qty::Temperature;
use fs_rep_mesh::{Metric2, TriComplex2};

fn source(
    value: PropertyValue,
    validity: ValidityDomain,
    basis: ElectricalPropertyBasis,
) -> MaterialCard {
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(basis.property_name(), basis.dimensions()),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic spatial Ohmic conductor".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(
    card: &MaterialCard,
    point: &QueryPoint,
    basis: ElectricalPropertyBasis,
) -> ResolvedConductorState {
    resolve_conductor_state_with_basis(
        card,
        point,
        MaterialPropertySelection::SingleClaimOnly,
        basis,
    )
    .unwrap()
}
fn conductor(sigma: f64, basis: ElectricalPropertyBasis) -> ResolvedConductorState {
    let value = match basis {
        ElectricalPropertyBasis::Resistivity => sigma.recip(),
        ElectricalPropertyBasis::Conductivity => sigma,
    };
    resolve(
        &source(
            PropertyValue::Scalar {
                value,
                dims: basis.dimensions(),
            },
            ValidityDomain::unconstrained(),
            basis,
        ),
        &QueryPoint::new(),
        basis,
    )
}
fn controls() -> ConductionControls {
    ConductionControls {
        linear: config(),
        current_tolerance_a: 1e-10,
        power_tolerance_w: 1e-10,
    }
}
fn slab(
    left: ResolvedConductorState,
    right: ResolvedConductorState,
    depth: f64,
) -> PlanarConduction {
    PlanarConduction::new(
        layered(0.5, 1.0, depth),
        vec![left.clone(), left, right.clone(), right],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn prescribed_current_uses_actual_edge_area_and_material_field() {
    let low = conductor(2.0, ElectricalPropertyBasis::Resistivity);
    let high = conductor(4.0, ElectricalPropertyBasis::Conductivity);
    for (depth, offset, density) in [(0.2, 0.0, 3.0), (0.6, 8.0, -3.0)] {
        let problem = slab(high.clone(), low.clone(), depth);
        let current = [ConductionCurrentTrace {
            vertices: [5, 2],
            inward_current_density_a_m2: density,
        }];
        let field = problem
            .solve_with_current_traces(
                &[(0, offset), (3, offset)],
                &current,
                controls(),
                &mut || false,
            )
            .unwrap();
        // Independently integrated series slab: J is prescribed, while each
        // material determines its own voltage drop and irreversible J²/sigma.
        let drop = density * (0.5 / 4.0 + 0.5 / 2.0);
        near(field.voltage_v()[1], offset + density * 0.5 / 4.0);
        near(field.voltage_v()[2], offset + drop);
        near(field.voltage_v()[5], offset + drop);
        for i in 0..4 {
            let sigma = if i < 2 { 4.0 } else { 2.0 };
            near(field.electric_field_v_m()[i][0], -density / sigma);
            near(field.current_density_a_m2()[i][0], -density);
            near(field.current_density_a_m2()[i][1], 0.0);
            near(
                field.cell_power_w()[i],
                0.25 * depth * density.powi(2) / sigma,
            );
        }
        near(field.power_w(), depth * density * drop);
        for v in [2, 5] {
            near(
                field.nodal_applied_current_a()[v].unwrap(),
                0.5 * depth * density,
            );
            near(
                field.nodal_boundary_current_a()[v].unwrap(),
                0.5 * depth * density,
            );
        }
        near(
            field.nodal_boundary_current_a()[0].unwrap()
                + field.nodal_boundary_current_a()[3].unwrap(),
            -depth * density,
        );
        near(field.boundary_current_a(), 0.0);
        near(field.maximum_free_current_residual_a(), 0.0);
        near(field.power_balance_residual_w(), 0.0);
        let by_voltage = problem
            .solve(
                &[
                    (0, offset),
                    (3, offset),
                    (2, offset + drop),
                    (5, offset + drop),
                ],
                controls(),
                &mut || false,
            )
            .unwrap();
        for (a, b) in field.cell_power_w().iter().zip(by_voltage.cell_power_w()) {
            near(*a, *b);
        }
    }
}

#[test]
fn distributed_currents_and_voltage_corners_have_one_source_work_account() {
    let material = conductor(2.0, ElectricalPropertyBasis::Conductivity);
    let problem = slab(material.clone(), material, 0.2);
    // The two top edges meet at node 4. Node 3 also lies on the voltage trace;
    // its prescribed-current contribution must not become a second source.
    let currents = [
        ConductionCurrentTrace {
            vertices: [3, 4],
            inward_current_density_a_m2: 2.0,
        },
        ConductionCurrentTrace {
            vertices: [4, 5],
            inward_current_density_a_m2: 4.0,
        },
    ];
    let baseline = problem
        .solve_with_current_traces(&[(0, 0.0), (3, 0.0)], &currents, controls(), &mut || false)
        .unwrap();
    let reversed = [
        ConductionCurrentTrace {
            vertices: [5, 4],
            ..currents[1]
        },
        ConductionCurrentTrace {
            vertices: [4, 3],
            ..currents[0]
        },
    ];
    let repeat = problem
        .solve_with_current_traces(&[(0, 0.0), (3, 0.0)], &reversed, controls(), &mut || false)
        .unwrap();
    assert_eq!(baseline.voltage_v(), repeat.voltage_v());
    assert_eq!(baseline.cell_power_w(), repeat.cell_power_w());
    near(baseline.nodal_applied_current_a()[3].unwrap(), 0.1);
    near(baseline.nodal_applied_current_a()[4].unwrap(), 0.3);
    near(baseline.nodal_applied_current_a()[5].unwrap(), 0.2);
    let mut source_work = 0.0;
    for v in 0..6 {
        let applied = baseline.nodal_applied_current_a()[v].unwrap_or(0.0);
        source_work += baseline.voltage_v()[v] * applied;
        if [0, 3].contains(&v) {
            let voltage_source = baseline.nodal_boundary_current_a()[v].unwrap() - applied;
            source_work += baseline.voltage_v()[v] * voltage_source;
        }
    }
    near(source_work, baseline.power_w());
    near(baseline.boundary_current_a(), 0.0);
    let shifted = problem
        .solve_with_current_traces(&[(0, 7.0), (3, 7.0)], &currents, controls(), &mut || false)
        .unwrap();
    near(shifted.power_w(), baseline.power_w());
    for (a, b) in shifted.voltage_v().iter().zip(baseline.voltage_v()) {
        near(*a - 7.0, *b);
    }
}

#[test]
fn current_traces_refuse_interior_duplicate_overconstrained_and_cancelled_loads() {
    let material = conductor(2.0, ElectricalPropertyBasis::Conductivity);
    let problem = slab(material.clone(), material, 0.2);
    let current = ConductionCurrentTrace {
        vertices: [2, 5],
        inward_current_density_a_m2: 3.0,
    };
    for bad in [
        vec![ConductionCurrentTrace {
            vertices: [1, 4],
            ..current
        }],
        vec![ConductionCurrentTrace {
            vertices: [0, 3],
            ..current
        }],
        vec![ConductionCurrentTrace {
            vertices: [2, usize::MAX],
            ..current
        }],
        vec![ConductionCurrentTrace {
            inward_current_density_a_m2: f64::NAN,
            ..current
        }],
        vec![
            current,
            ConductionCurrentTrace {
                vertices: [5, 2],
                ..current
            },
        ],
    ] {
        assert!(matches!(
            problem
                .solve_with_current_traces(&[(0, 0.0), (3, 0.0)], &bad, controls(), &mut || false,),
            Err(ConductionError::Invalid(_))
        ));
    }
    assert!(
        problem
            .solve_with_current_traces(&[], &[current], controls(), &mut || false,)
            .is_err()
    );
    let mut polls = 0;
    let accepted = problem
        .solve_with_current_traces(&[(0, 0.0), (3, 0.0)], &[current], controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        problem.solve_with_current_traces(
            &[(0, 0.0), (3, 0.0)],
            &[current],
            controls(),
            &mut || {
                calls += 1;
                calls == polls
            },
        ),
        Err(ConductionError::Cancelled)
    ));
    let again = problem
        .solve_with_current_traces(&[(0, 0.0), (3, 0.0)], &[current], controls(), &mut || false)
        .unwrap();
    assert_eq!(accepted.voltage_v(), again.voltage_v());
}

#[test]
fn series_material_swap_redistributes_electric_field_and_physical_heat() {
    let low = conductor(2.0, ElectricalPropertyBasis::Resistivity);
    let high = conductor(4.0, ElectricalPropertyBasis::Conductivity);
    let uniform = slab(low.clone(), low.clone(), 0.2);
    let changed = slab(high.clone(), low.clone(), 0.2);
    let baseline = uniform
        .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    let reduced = changed
        .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    near(
        baseline.resistance_ohm(),
        low.clone()
            .into_uniform_conductor(1.0, 0.2)
            .unwrap()
            .resistance_ohm(),
    );
    // Independent piecewise-series result, including actual extrusion area.
    let resistance = (0.5 / 4.0 + 0.5 / 2.0) / 0.2;
    near(reduced.resistance_ohm(), resistance);
    assert!(reduced.resistance_ohm() < baseline.resistance_ohm());
    let field = reduced.unit_voltage_field();
    let current = resistance.recip();
    for cell in 0..4 {
        let expected_e = if cell < 2 { -2.0 / 3.0 } else { -4.0 / 3.0 };
        near(field.electric_field_v_m()[cell][0], expected_e);
        near(field.electric_field_v_m()[cell][1], 0.0);
        near(field.current_density_a_m2()[cell][0], -current / 0.2);
        near(field.current_density_a_m2()[cell][1], 0.0);
        near(
            field.cell_power_w()[cell],
            0.05 * expected_e * (-current / 0.2),
        );
        assert_eq!(
            &field.problem().conductors()[cell],
            if cell < 2 { &high } else { &low }
        );
    }
    // More resistive material heats MORE in this series path at equal current.
    near(field.cell_power_w()[2], 2.0 * field.cell_power_w()[0]);
    near(field.power_w(), current);
    let i = field.nodal_boundary_current_a();
    near(i[2].unwrap() + i[5].unwrap(), current);
    near(i[0].unwrap() + i[3].unwrap(), -current);
    assert!(i[1].is_none() && i[4].is_none());
    assert!(field.report().unwrap().converged_euclidean());
    near(field.boundary_current_a(), 0.0);
    near(field.maximum_free_current_residual_a(), 0.0);
    near(field.power_balance_residual_w(), 0.0);
}

#[test]
fn parallel_paths_thickness_voltage_reversal_and_gauge_preserve_power() {
    let left = conductor(4.0, ElectricalPropertyBasis::Conductivity);
    let right = conductor(2.0, ElectricalPropertyBasis::Resistivity);
    let mut power = Vec::new();
    for depth in [0.2, 0.6] {
        let problem = slab(left.clone(), right.clone(), depth);
        // Drive the OTHER pair of edges: these materials now conduct in parallel.
        let resistor = problem
            .resistance(&[3, 4, 5], &[0, 1, 2], controls(), &mut || false)
            .unwrap();
        let conductance = depth * (0.5 * 4.0 + 0.5 * 2.0);
        near(resistor.resistance_ohm(), 1.0 / conductance);
        let forward = resistor
            .field_at_voltage(3.0, 1e-10, 1e-10, &mut || false)
            .unwrap();
        let reverse = resistor
            .field_at_voltage(-3.0, 1e-10, 1e-10, &mut || false)
            .unwrap();
        let offset = problem
            .solve(
                &[
                    (0, 8.0),
                    (1, 8.0),
                    (2, 8.0),
                    (3, 11.0),
                    (4, 11.0),
                    (5, 11.0),
                ],
                controls(),
                &mut || false,
            )
            .unwrap();
        for cell in 0..4 {
            near(forward.electric_field_v_m()[cell][0], 0.0);
            near(forward.electric_field_v_m()[cell][1], -3.0);
            near(
                forward.current_density_a_m2()[cell][1],
                if cell < 2 { -12.0 } else { -6.0 },
            );
            for axis in 0..2 {
                near(
                    reverse.current_density_a_m2()[cell][axis],
                    -forward.current_density_a_m2()[cell][axis],
                );
                near(
                    offset.current_density_a_m2()[cell][axis],
                    forward.current_density_a_m2()[cell][axis],
                );
            }
            near(reverse.cell_power_w()[cell], forward.cell_power_w()[cell]);
        }
        near(forward.power_w(), 9.0 * conductance);
        near(offset.power_w(), forward.power_w());
        // In parallel the more conductive material has MORE heat at equal E.
        near(forward.cell_power_w()[0], 2.0 * forward.cell_power_w()[2]);
        let zero = resistor
            .field_at_voltage(0.0, 1e-10, 1e-10, &mut || false)
            .unwrap();
        assert_eq!(zero.power_w(), 0.0);
        assert_eq!(zero.current_density_a_m2(), [[0.0; 2]; 4]);
        power.push(forward.power_w());
    }
    near(power[1], 3.0 * power[0]);
}

#[test]
fn common_voltage_shift_preserves_a_free_interior_node_current_and_heat() {
    let mesh = TriComplex2::from_indexed_triangles(
        "conductor-free-interior-gauge",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.25, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(0.2).unwrap(),
    )
    .unwrap();
    let problem = PlanarConduction::new(
        mesh,
        vec![conductor(2.0, ElectricalPropertyBasis::Conductivity); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let boundary = [(0, 0.0), (1, 4.0), (2, 4.0), (3, 0.0)];
    let original = problem.solve(&boundary, controls(), &mut || false).unwrap();
    let reference = 2.0_f64.powi(40);
    let shifted = problem
        .solve(
            &boundary.map(|(node, voltage)| (node, reference + voltage)),
            controls(),
            &mut || false,
        )
        .unwrap();
    // Independent affine solution V=4*x: the unconstrained INTERIOR node is
    // at one volt, E=(-4,0), J=(-8,0), and sigma*|E|²*volume gives each heat.
    // Integer voltage differences remain exactly representable at this gauge;
    // the offset must not set the tolerance for the free-node solve.
    near(original.voltage_v()[4], 1.0);
    near(shifted.voltage_v()[4] - reference, 1.0);
    assert!(original.report().unwrap().converged_euclidean());
    assert!(shifted.report().unwrap().converged_euclidean());
    assert!(shifted.nodal_boundary_current_a()[4].is_none());
    for (cell, area) in [0.25, 0.375, 0.25, 0.125].into_iter().enumerate() {
        for field in [&original, &shifted] {
            near(field.electric_field_v_m()[cell][0], -4.0);
            near(field.electric_field_v_m()[cell][1], 0.0);
            near(field.current_density_a_m2()[cell][0], -8.0);
            near(field.current_density_a_m2()[cell][1], 0.0);
            near(field.cell_power_w()[cell], 6.4 * area);
        }
    }
    near(original.power_w(), 6.4);
    near(shifted.power_w(), original.power_w());
    near(shifted.maximum_free_current_residual_a(), 0.0);
    near(shifted.power_balance_residual_w(), 0.0);
}

#[test]
fn regional_conductor_changes_native_rc_discharge_and_partitions_its_joule_heat() {
    let mut endpoints = Vec::new();
    for sigma_left in [2.0, 4.0] {
        let problem = slab(
            conductor(sigma_left, ElectricalPropertyBasis::Conductivity),
            conductor(2.0, ElectricalPropertyBasis::Resistivity),
            0.2,
        );
        let resistor = problem
            .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
            .unwrap();
        let capacitance = 0.4;
        let dae = assemble_circuit(&CircuitGraph {
            node_count: 2,
            branches: vec![
                (
                    1,
                    0,
                    Branch::Capacitor {
                        farads: capacitance,
                    },
                ),
                (1, 0, resistor.circuit_branch()),
            ],
            transformers: vec![],
        })
        .unwrap();
        let mut initial = vec![0.0; dae.system.state_dim()];
        initial[dae.charge_index[0]] = capacitance;
        let mut state = dae.consistent_initial_state(&initial, &[]).unwrap();
        let dt = 0.1;
        let ratio = (1.0 - dt / (2.0 * resistor.resistance_ohm() * capacitance))
            / (1.0 + dt / (2.0 * resistor.resistance_ohm() * capacitance));
        let mut total_heat = 0.0;
        for _ in 0..3 {
            let before = state[dae.charge_index[0]] / capacitance;
            let (step, defect) = dae.step_audited(&state, &[], dt).unwrap();
            let after = step.x[dae.charge_index[0]] / capacitance;
            near(after, before * ratio);
            let field = resistor
                .field_at_voltage(0.5 * (before + after), 1e-10, 1e-10, &mut || false)
                .unwrap();
            let spatial_heat: f64 = field.cell_power_w().iter().map(|power| dt * power).sum();
            near(spatial_heat, step.dissipated);
            near(dt * field.power_w(), step.dissipated);
            near(
                0.5 * capacitance * (before * before - after * after),
                step.dissipated,
            );
            assert!(defect < 1e-10 && spatial_heat > 0.0);
            total_heat += spatial_heat;
            state = step.x;
        }
        let final_voltage = state[dae.charge_index[0]] / capacitance;
        near(
            total_heat + dae.system.hamiltonian(&state),
            0.5 * capacitance,
        );
        endpoints.push(final_voltage);
    }
    assert!(
        endpoints[1] < endpoints[0],
        "the material-only swap accelerates discharge"
    );
}

#[test]
fn temperature_and_independent_property_basis_refresh_the_same_physical_operator() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source(
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 0.5), (400.0, 1.0)],
            dims: basis.dimensions(),
        },
        ValidityDomain::unconstrained().with("T", 300.0, 400.0),
        basis,
    );
    let cold = resolve(&card, &QueryPoint::new().with("T", 300.0).unwrap(), basis);
    let warm_snapshot = cold
        .material()
        .at_query_point(&card, &QueryPoint::new().with("T", 400.0).unwrap())
        .unwrap();
    let warm = ResolvedConductorState::try_from_material(warm_snapshot).unwrap();
    let original = slab(cold.clone(), cold, 0.2);
    let changed = slab(warm.clone(), warm, 0.2);
    let same = slab(
        conductor(2.0, ElectricalPropertyBasis::Conductivity),
        conductor(2.0, ElectricalPropertyBasis::Conductivity),
        0.2,
    );
    let a = original
        .solve(&ends(1.0), controls(), &mut || false)
        .unwrap();
    let b = changed
        .solve(&ends(1.0), controls(), &mut || false)
        .unwrap();
    let c = same.solve(&ends(1.0), controls(), &mut || false).unwrap();
    near(b.power_w(), 0.5 * a.power_w());
    near(c.power_w(), a.power_w());
    assert_eq!(a.electric_field_v_m(), c.electric_field_v_m());
    for model in [&original, &changed] {
        assert_eq!(
            model.conductors()[0].material().card_identity(),
            card.content_hash()
        );
        assert_eq!(model.conductors()[0].source_basis(), basis);
    }
    assert_ne!(original.conductors()[0], changed.conductors()[0]);
    assert!(
        resolve_conductor_state_with_basis(
            &card,
            &QueryPoint::new().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            basis,
        )
        .is_err()
    );
    // Re-querying and refusing another state cannot change the accepted field.
    assert_eq!(
        a.power_w(),
        original
            .solve(&ends(1.0), controls(), &mut || false)
            .unwrap()
            .power_w()
    );
}

#[test]
fn source_model_electrode_work_and_late_cancellation_refuse_without_partial_field() {
    let basis = ElectricalPropertyBasis::Conductivity;
    let card = source(
        PropertyValue::Scalar {
            value: 2.0,
            dims: basis.dimensions(),
        },
        ValidityDomain::unconstrained().with("frequency", 100.0, 1000.0),
        basis,
    );
    let ac = resolve(
        &card,
        &QueryPoint::new().with("frequency", 500.0).unwrap(),
        basis,
    );
    assert!(matches!(
        PlanarConduction::new(layered(0.5, 1.0, 0.2), vec![ac; 4], limits(), &mut || false),
        Err(ConductionError::Invalid(_))
    ));
    let source = conductor(2.0, basis);
    assert!(
        PlanarConduction::new(
            layered(0.5, 1.0, 0.2),
            vec![source.clone(); 3],
            limits(),
            &mut || false
        )
        .is_err()
    );
    let mut limited = limits();
    limited.max_cells = 3;
    assert!(
        PlanarConduction::new(
            layered(0.5, 1.0, 0.2),
            vec![source.clone(); 4],
            limited,
            &mut || false
        )
        .is_err()
    );
    let problem = slab(source.clone(), source, 0.2);
    for (positive, negative) in [
        (vec![2], vec![0, 3]),
        (vec![2, 5], vec![2, 3]),
        (vec![2, 2, 5], vec![0, 3]),
    ] {
        assert!(
            problem
                .resistance(&positive, &negative, controls(), &mut || false)
                .is_err()
        );
    }
    let mut no_work = controls();
    no_work.linear.max_iterations = 0;
    assert!(matches!(
        problem.solve(&ends(1.0), no_work, &mut || false),
        Err(ConductionError::Nonconvergence(_))
    ));
    assert!(
        problem
            .solve(&[(0, 0.0)], controls(), &mut || false)
            .is_err()
    );
    let resistor = problem
        .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    let mut count = 0;
    let accepted = resistor
        .field_at_voltage(2.0, 1e-10, 1e-10, &mut || {
            count += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    let refused = resistor.field_at_voltage(2.0, 1e-10, 1e-10, &mut || {
        calls += 1;
        calls == count
    });
    assert!(matches!(refused, Err(ConductionError::Cancelled)));
    let again = resistor
        .field_at_voltage(2.0, 1e-10, 1e-10, &mut || false)
        .unwrap();
    assert_eq!(accepted.voltage_v(), again.voltage_v());
    assert_eq!(accepted.cell_power_w(), again.cell_power_w());
    assert!(
        resistor
            .field_at_voltage(f64::NAN, 1e-10, 1e-10, &mut || false)
            .is_err()
    );
    assert!(
        resistor
            .field_at_voltage(f64::MAX, 1e-10, 1e-10, &mut || false)
            .is_err()
    );
}
