//! G1/G3/G4: spatial source-defined Ohmic current, native circuit loss and heat.
//! These are synthetic physical-law checks, not measured material qualification.
mod common;

use common::{config, ends, layered, limits, near};
use fs_em::conduction::{ConductionControls, ConductionError, PlanarConduction};
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
