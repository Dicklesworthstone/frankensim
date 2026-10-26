//! Synthetic source-data regressions, not experimentally qualified dielectrics.
mod common;
#[path = "electrostatic/force.rs"]
mod force;
#[path = "electrostatic/network.rs"]
mod network;
use common::{config, layered, limits};
use fs_em::electrostatic::{ElectrostaticControls, ElectrostaticError, PlanarElectrostatics};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, ResolvedCapacitor,
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_qty::Temperature;
use fs_rep_mesh::{Metric2, TriComplex2};

fn source(value: PropertyValue, validity: ValidityDomain) -> MaterialCard {
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic spatial dielectric".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(card: &MaterialCard, point: &QueryPoint) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        point,
        &[ScalarPropertyRequirement::try_new(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn dielectric(epsilon: f64) -> ResolvedMaterialStatePoint {
    resolve(
        &source(
            PropertyValue::Scalar {
                value: epsilon,
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            },
            ValidityDomain::unconstrained(),
        ),
        &QueryPoint::new(),
    )
}
fn controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: config(),
        charge_tolerance_c: 1e-15,
        energy_tolerance_j: 1e-15,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-9 * b.abs().max(1e-20), "{a} != {b}");
}
fn slab(
    left: ResolvedMaterialStatePoint,
    right: ResolvedMaterialStatePoint,
) -> PlanarElectrostatics {
    PlanarElectrostatics::new(
        layered(0.5, 1.0, 0.2),
        vec![left.clone(), left, right.clone(), right],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn dielectric_swap_redistributes_field_and_matches_series_capacitance() {
    let low = dielectric(2e-11);
    let high = dielectric(4e-11);
    let uniform = slab(low.clone(), low.clone());
    let layered = slab(high.clone(), low.clone());
    let baseline = uniform
        .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    let changed = layered
        .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    let reference = ResolvedCapacitor::try_from_material(low.clone(), 0.2, 1.0).unwrap();
    near(baseline.capacitance_f(), reference.capacitance_f());
    // Independent two-slab series formula; geometry and solver are unchanged.
    let expected = 0.2 / (0.5 / 4e-11 + 0.5 / 2e-11);
    near(changed.capacitance_f(), expected);
    assert!(changed.capacitance_f() > baseline.capacitance_f());
    for cell in 0..4 {
        let field = changed.unit_voltage_field();
        near(
            baseline.unit_voltage_field().electric_field_v_m()[cell][0],
            -1.0,
        );
        near(
            field.electric_field_v_m()[cell][0],
            if cell < 2 { -2.0 / 3.0 } else { -4.0 / 3.0 },
        );
        assert!(field.electric_field_v_m()[cell][1].abs() < 1e-12);
        near(field.displacement_c_m2()[cell][0], -expected / 0.2);
        assert_eq!(
            &field.problem().materials()[cell],
            if cell < 2 { &high } else { &low }
        );
    }
    near(changed.unit_voltage_field().energy_j(), 0.5 * expected);
    let q = changed.unit_voltage_field().nodal_boundary_charge_c();
    near(q[2].unwrap() + q[5].unwrap(), expected);
    near(q[0].unwrap() + q[3].unwrap(), -expected);
    assert!(q[1].is_none() && q[4].is_none());
    assert!(
        changed
            .unit_voltage_field()
            .report()
            .unwrap()
            .converged_euclidean()
    );
    // Physical voltage scaling plus a common offset: the field and energy
    // depend on the voltage difference, never the arbitrary gauge.
    let driven = layered
        .solve(
            &[0.0; 4],
            &[(0, 2.0), (3, 2.0), (2, 8.0), (5, 8.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(driven.energy_j(), 0.5 * expected * 36.0);
    near(driven.electric_field_v_m()[0][0], -4.0);
    near(driven.electric_field_v_m()[2][0], -8.0);
    assert!(driven.charge_balance_residual_c().abs() < 1e-24);
    assert!(driven.energy_balance_residual_j().abs() < 1e-23);
}

#[test]
fn free_charge_and_electrode_reactions_close_gauss_and_energy() {
    let epsilon = 4e-11;
    let mesh = || {
        TriComplex2::from_indexed_triangles(
            "dielectric-volume-charge",
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
            vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
            Metric2::planar(0.25).unwrap(),
        )
        .unwrap()
    };
    let problem = |epsilon| {
        PlanarElectrostatics::new(mesh(), vec![dielectric(epsilon); 4], limits(), &mut || {
            false
        })
        .unwrap()
    };
    let original = problem(epsilon);
    let changed = problem(2.0 * epsilon);
    let boundary = [(0, 0.0), (1, 0.0), (2, 0.0), (3, 0.0)];
    let rho = 12.0 * epsilon;
    let a = original
        .solve(&[rho; 4], &boundary, controls(), &mut || false)
        .unwrap();
    let b = changed
        .solve(&[rho; 4], &boundary, controls(), &mut || false)
        .unwrap();
    // K_center=4*depth*epsilon, consistent load=rho*depth/3.
    near(a.voltage_v()[4], 1.0);
    near(b.voltage_v()[4], 0.5);
    near(a.volume_charge_c(), 0.25 * rho);
    near(a.boundary_charge_c(), -0.25 * rho);
    near(a.energy_j(), 0.5 * epsilon);
    near(b.energy_j(), 0.25 * epsilon);
    assert!(a.maximum_free_charge_residual_c() < 1e-24);
    assert!(a.charge_balance_residual_c().abs() < 1e-24);
    let cancelled = original.solve(&[rho; 4], &boundary, controls(), &mut || true);
    assert!(matches!(cancelled, Err(ElectrostaticError::Cancelled)));
    let again = original
        .solve(&[rho; 4], &boundary, controls(), &mut || false)
        .unwrap();
    assert_eq!(a.voltage_v(), again.voltage_v());
    assert_eq!(a.energy_j(), again.energy_j());
}

#[test]
fn spatial_dielectric_changes_native_circuit_charge_and_decay() {
    let mut voltages = Vec::new();
    for epsilon in [4e-6, 8e-6] {
        // Deliberately synthetic high epsilon keeps this tiny circuit fixture
        // well scaled. The independent analytical result is C=epsilon*A/d.
        let problem = slab(dielectric(epsilon), dielectric(epsilon));
        let capacitor = problem
            .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
            .unwrap();
        let c = capacitor.capacitance_f();
        near(c, 0.2 * epsilon);
        let circuit = assemble_circuit(&CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, capacitor.circuit_branch()),
                (1, 0, Branch::Resistor { ohms: 1000.0 }),
            ],
            transformers: vec![],
        })
        .unwrap();
        let mut initial = vec![0.0; circuit.system.state_dim()];
        initial[circuit.charge_index[0]] = c; // identical one-volt initial condition
        let initial = circuit.consistent_initial_state(&initial, &[]).unwrap();
        let (step, defect) = circuit.step_audited(&initial, &[], 1e-4).unwrap();
        let ratio = (1.0 - 1e-4 / (2.0 * 1000.0 * c)) / (1.0 + 1e-4 / (2.0 * 1000.0 * c));
        near(step.x[circuit.charge_index[0]], c * ratio);
        near(step.x[circuit.charge_index[0]] / c, ratio);
        near(step.dissipated, 0.5 * c * (1.0 - ratio * ratio));
        assert!(defect < 1e-15 && step.dissipated > 0.0);
        voltages.push(step.x[circuit.charge_index[0]] / c);
    }
    assert!(
        voltages[1] > voltages[0],
        "the data-only dielectric swap increases RC storage time"
    );
}

#[test]
fn fixed_source_temperature_refresh_changes_operator_without_erasing_identity() {
    let card = source(
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 2e-11), (400.0, 4e-11)],
            dims: ELECTRIC_PERMITTIVITY_DIMS,
        },
        ValidityDomain::unconstrained().with("T", 300.0, 400.0),
    );
    let cold = resolve(&card, &QueryPoint::new().with("T", 300.0).unwrap());
    let warm = resolve(&card, &QueryPoint::new().with("T", 400.0).unwrap());
    let original = slab(cold.clone(), cold);
    let refreshed = slab(warm.clone(), warm);
    let a = original
        .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    let b = refreshed
        .capacitance(&[2, 5], &[0, 3], controls(), &mut || false)
        .unwrap();
    near(b.capacitance_f(), 2.0 * a.capacitance_f());
    assert_eq!(original.materials()[0].card_identity(), card.content_hash());
    assert_eq!(
        refreshed.materials()[0].card_identity(),
        card.content_hash()
    );
    assert_ne!(original.materials()[0], refreshed.materials()[0]);
    assert!(
        ResolvedCapacitor::resolve(
            &card,
            &QueryPoint::new().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            0.2,
            1.0
        )
        .is_err()
    );
}

#[test]
fn unsupported_sources_geometry_work_and_electrodes_refuse() {
    let frequency = source(
        PropertyValue::Scalar {
            value: 2e-11,
            dims: ELECTRIC_PERMITTIVITY_DIMS,
        },
        ValidityDomain::unconstrained().with("frequency", 100.0, 1000.0),
    );
    let limited = resolve(
        &frequency,
        &QueryPoint::new().with("frequency", 500.0).unwrap(),
    );
    assert!(matches!(
        PlanarElectrostatics::new(
            layered(0.5, 1.0, 0.2),
            vec![limited; 4],
            limits(),
            &mut || false
        ),
        Err(ElectrostaticError::Material(_))
    ));
    assert!(matches!(
        PlanarElectrostatics::new(
            layered(0.5, 1.0, 0.2),
            vec![common::material(1e-6); 4],
            limits(),
            &mut || false
        ),
        Err(ElectrostaticError::Material(_))
    ));
    let small = fs_feec::weighted::WeightedAssemblyLimits {
        max_cells: 3,
        ..limits()
    };
    assert!(matches!(
        PlanarElectrostatics::new(
            layered(0.5, 1.0, 0.2),
            vec![dielectric(2e-11); 4],
            small,
            &mut || false
        ),
        Err(ElectrostaticError::Invalid(_))
    ));
    let problem = slab(dielectric(2e-11), dielectric(2e-11));
    for (positive, negative) in [
        (vec![2], vec![0, 3]),
        (vec![2, 5], vec![2, 3]),
        (vec![2, 2, 5], vec![0, 3]),
    ] {
        assert!(matches!(
            problem.capacitance(&positive, &negative, controls(), &mut || false),
            Err(ElectrostaticError::Invalid(_))
        ));
    }
    let mut no_work = controls();
    no_work.linear.max_iterations = 0;
    assert!(matches!(
        problem.capacitance(&[2, 5], &[0, 3], no_work, &mut || false),
        Err(ElectrostaticError::Nonconvergence(_))
    ));
    assert!(matches!(
        problem.solve(&[0.0; 4], &[(0, 0.0)], controls(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        problem.solve(&[0.0; 3], &[(0, 0.0), (3, 0.0)], controls(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        problem.capacitance(&[2, 5], &[0, 3], controls(), &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
}
