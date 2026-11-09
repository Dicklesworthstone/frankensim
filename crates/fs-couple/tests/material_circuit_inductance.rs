//! G1/G3 source-to-terminal magnetic work, declared gap/leakage and native RL.
//! Manufactured H(B) sources are implementation oracles, not measured devices.
#[path = "material_circuit_inductance/shared.rs"]
mod shared;

use fs_couple::material_circuit::inductance::{
    MaterialReversibleInductor, ReversibleInductorCircuit, ReversibleInductorGeometry,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
use fs_qty::QuantitySpec;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn material(points: &[(f64, f64)]) -> ReversibleMagneticMaterial {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(REVERSIBLE_FIELD_STRENGTH_PROPERTY, FIELD_STRENGTH_DIMS),
            value: PropertyValue::Curve {
                abscissa: "B".into(),
                abscissa_dims: FLUX_DENSITY_DIMS,
                knots: points.to_vec(),
                dims: FIELD_STRENGTH_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "B",
                QuantitySpec::dimensional(FLUX_DENSITY_DIMS),
                0.0,
                points.last().unwrap().0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured saturating magnetic path".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic nonlinear core".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap();
    ReversibleMagneticMaterial::resolve(
        &card,
        &point,
        "B",
        points.last().unwrap().0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn geometry() -> ReversibleInductorGeometry {
    ReversibleInductorGeometry {
        turns: 3.0,
        core_area_m2: 2.0,
        core_length_m: 4.0,
        gap_reluctance_per_h: 5.0,
        leakage_inductance_h: 0.2,
    }
}

#[test]
fn every_source_knot_maps_to_native_storage_and_core_gap_leakage_energy() {
    let law = material(&[(0.0, 0.0), (1.0, 2.0), (3.0, 14.0)]);
    let inductor = MaterialReversibleInductor::new(0, law.clone(), geometry()).unwrap();
    assert_eq!(inductor.material(), &law);
    assert_eq!(inductor.curve().knots().len(), 5);
    for (index, b) in [-3.0_f64, -1.0, 0.0, 1.0, 3.0].into_iter().enumerate() {
        let h = b.signum()
            * if b.abs() <= 1.0 {
                2.0 * b.abs()
            } else {
                6.0 * b.abs() - 4.0
            };
        let current = (4.0 * h + 10.0 * b) / 3.0;
        let linkage = 6.0 * b + 0.2 * current;
        near(inductor.curve().knots()[index].0, linkage, 2e-14);
        near(inductor.curve().knots()[index].1, current, 2e-14);
        assert_eq!(
            inductor
                .operating_point(inductor.curve().knots()[index].0)
                .unwrap()
                .core_flux_density_t,
            b
        );
    }
    for b in [-3.0_f64, -2.0, -0.25, 0.0, 0.25, 2.0, 3.0] {
        let h = b.signum()
            * if b.abs() <= 1.0 {
                2.0 * b.abs()
            } else {
                6.0 * b.abs() - 4.0
            };
        let density = if b.abs() <= 1.0 {
            b * b
        } else {
            3.0 * b * b - 4.0 * b.abs() + 2.0
        };
        let current = (4.0 * h + 10.0 * b) / 3.0;
        let linkage = 6.0 * b + 0.2 * current;
        let point = inductor.operating_point(linkage).unwrap();
        near(point.current_a, current, 2e-14);
        near(point.core_flux_density_t, b, 2e-14);
        near(point.core_field_strength_a_m, h, 3e-14);
        near(point.core_energy_j, 8.0 * density, 2e-12);
        near(point.gap_energy_j, 10.0 * b * b, 2e-12);
        near(point.leakage_energy_j, 0.1 * current * current, 2e-12);
        near(
            point.stored_energy_j,
            point.core_energy_j + point.gap_energy_j + point.leakage_energy_j,
            4e-12,
        );
        let slope = if b.abs() < 1.0 { 2.0 } else { 6.0 };
        near(
            point.differential_inductance_h,
            0.2 + 18.0 / (4.0 * slope + 10.0),
            2e-14,
        );
        near(inductor.flux_linkage_wb(current).unwrap(), linkage, 2e-14);
    }
    let point = inductor
        .operating_point(inductor.flux_linkage_wb(52.0 / 3.0).unwrap())
        .unwrap();
    assert!((point.stored_energy_j - 0.5 * point.flux_linkage_wb * point.current_a).abs() > 1.0);
}

#[test]
fn turns_scale_flux_current_and_incremental_inductance_without_changing_field_energy() {
    let law = material(&[(0.0, 0.0), (1.0, 2.0), (3.0, 14.0)]);
    let mut g = geometry();
    g.leakage_inductance_h = 0.0;
    let first = MaterialReversibleInductor::new(0, law.clone(), g).unwrap();
    g.turns *= 2.0;
    let second = MaterialReversibleInductor::new(0, law, g).unwrap();
    for b in [-2.5, -0.25, 0.25, 2.5] {
        let a = first.operating_point(6.0 * b).unwrap();
        let c = second.operating_point(12.0 * b).unwrap();
        near(c.current_a, 0.5 * a.current_a, 2e-14);
        near(c.stored_energy_j, a.stored_energy_j, 2e-12);
        near(
            c.differential_inductance_h,
            4.0 * a.differential_inductance_h,
            2e-14,
        );
        near(c.core_flux_density_t, a.core_flux_density_t, 2e-14);
    }
}

#[test]
fn leakage_dominated_linkage_keeps_resolvable_core_field_observations() {
    let g = ReversibleInductorGeometry {
        turns: 1.0,
        core_area_m2: 1.0,
        core_length_m: 1.0,
        gap_reluctance_per_h: 0.0,
        leakage_inductance_h: 1e16,
    };
    let inductor =
        MaterialReversibleInductor::new(0, material(&[(0.0, 0.0), (1.0, 2.0), (3.0, 14.0)]), g)
            .unwrap();
    let lambda = inductor.flux_linkage_wb(1.5).unwrap();
    let point = inductor.operating_point(lambda).unwrap();
    near(point.core_flux_density_t, 0.75, 2e-15);
    near(point.core_energy_j, 0.75 * 0.75, 3e-15);
    let cancellation_prone_field = lambda - g.leakage_inductance_h * point.current_a;
    assert!((cancellation_prone_field - point.core_flux_density_t).abs() > 0.25);
}

fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 1.0 }),
            (2, 0, Branch::Inductor { henries: 999.0 }),
        ],
        transformers: vec![],
    }
}
fn driven_source() -> MaterialReversibleInductor {
    MaterialReversibleInductor::new(
        2,
        material(&[(0.0, 0.0), (1.0, 1.0), (2.0, 6.0)]),
        ReversibleInductorGeometry {
            turns: 1.0,
            core_area_m2: 1.0,
            core_length_m: 1.0,
            gap_reluctance_per_h: 0.0,
            leakage_inductance_h: 0.0,
        },
    )
    .unwrap()
}

#[test]
fn sourced_saturation_drives_the_original_circuit_and_closes_work_and_resistor_heat() {
    let source = driven_source();
    let circuit = ReversibleInductorCircuit::new(
        graph(),
        vec![source.clone()],
        &[],
        &[],
        &[],
        &[],
        &mut || false,
    )
    .unwrap();
    assert_eq!(circuit.sources(), &[source]);
    assert_eq!(
        circuit.graph().branches[2].2,
        Branch::Inductor { henries: 1.0 }
    );
    let mut state = circuit
        .consistent_initial_state(&[0.0], &[], &[3.0])
        .unwrap();
    let (mut work, mut heat) = (0.0, 0.0);
    for _ in 0..80 {
        let (step, _) = circuit
            .dae()
            .step_audited(&state, &[3.0], 1.0 / 80.0)
            .unwrap();
        near(step.delta_h + step.dissipated, step.supplied, 2e-11);
        assert!(step.dissipated >= 0.0);
        work += step.supplied;
        heat += step.dissipated;
        state = step.x;
        let point = circuit.operating_points(&state).unwrap()[0];
        near(point.stored_energy_j + heat, work, 3e-10);
        near(
            point.current_a,
            circuit.dae().inductor_currents(&state).unwrap()[0],
            2e-13,
        );
    }
    let crossing = 1.5_f64.ln();
    let exact = 1.4 - 0.4 * (-5.0 * (1.0 - crossing)).exp();
    near(state[circuit.dae().flux_index[0]], exact, 5e-5);
    assert!(circuit.operating_points(&state).unwrap()[0].core_flux_density_t > 1.0);
}

#[test]
fn geometry_source_support_branch_ownership_and_cancellation_refuse() {
    let law = material(&[(0.0, 0.0), (1.0, 2.0), (3.0, 14.0)]);
    for bad in [
        ReversibleInductorGeometry {
            turns: 0.0,
            ..geometry()
        },
        ReversibleInductorGeometry {
            core_area_m2: -1.0,
            ..geometry()
        },
        ReversibleInductorGeometry {
            core_length_m: f64::NAN,
            ..geometry()
        },
        ReversibleInductorGeometry {
            gap_reluctance_per_h: -1.0,
            ..geometry()
        },
        ReversibleInductorGeometry {
            leakage_inductance_h: -1.0,
            ..geometry()
        },
    ] {
        assert!(MaterialReversibleInductor::new(0, law.clone(), bad).is_err());
    }
    let source = driven_source();
    assert!(source.operating_point(2.00001).is_err());
    assert!(source.current_a(f64::NAN).is_err());
    assert!(source.flux_linkage_wb(6.00001).is_err());
    assert!(
        ReversibleInductorCircuit::new(
            graph(),
            vec![source.clone()],
            &[],
            &[],
            &[],
            &[],
            &mut || true,
        )
        .is_err()
    );
    assert!(
        ReversibleInductorCircuit::new(
            graph(),
            vec![source.clone(), source.clone()],
            &[],
            &[],
            &[],
            &[],
            &mut || false,
        )
        .is_err()
    );
    assert!(
        ReversibleInductorCircuit::new(
            graph(),
            vec![source.clone()],
            &[],
            &[],
            &[source.binding()],
            &[],
            &mut || false,
        )
        .is_err()
    );
    let mut coupled = graph();
    coupled
        .branches
        .push((2, 0, Branch::Inductor { henries: 1.0 }));
    assert!(
        ReversibleInductorCircuit::new(
            coupled,
            vec![source],
            &[],
            &[],
            &[],
            &[MutualInductance {
                branches: (2, 3),
                henries: 0.1
            }],
            &mut || false,
        )
        .is_err()
    );
}
