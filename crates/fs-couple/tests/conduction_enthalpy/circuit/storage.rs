//! Source-only dielectric/core swaps drive native storage, actual fields and heat.
use super::*;
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
    MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
};
use fs_qty::Dims;

fn card(name: &str, dims: Dims, value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic source storage; not a measured dielectric or magnetic core"
                    .into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic field material".into(),
            phase: "solid".into(),
            process: "fixed-state numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub(super) fn capacitor(branch: usize, epsilon: f64) -> MaterialStorage {
    MaterialStorage {
        branch,
        element: MaterialStorageElement::Capacitor(
            ResolvedCapacitor::resolve(
                &card(
                    ELECTRIC_PERMITTIVITY_PROPERTY,
                    ELECTRIC_PERMITTIVITY_DIMS,
                    epsilon,
                ),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                2.0,
                0.5,
            )
            .unwrap(),
        ),
    }
}

fn inductor(branch: usize, mu: f64) -> MaterialStorage {
    MaterialStorage {
        branch,
        element: MaterialStorageElement::Inductor(
            ResolvedInductor::resolve(
                &card(
                    MAGNETIC_PERMEABILITY_PROPERTY,
                    MAGNETIC_PERMEABILITY_DIMS,
                    mu,
                ),
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                0.5,
                1.0,
                4,
            )
            .unwrap(),
        ),
    }
}

#[test]
fn g1_source_only_dielectric_and_core_swaps_change_actual_storage_fields_and_heat() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 999.0, None);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let layout = network(&bindings);
    let dt = 0.25;
    for electric in [true, false] {
        let mut outcomes = Vec::new();
        for value in [0.125, 0.25] {
            // Only the material property changes; geometry, initial charge/flux,
            // current field, thermal cells and integration controls are identical.
            let source = if electric {
                capacitor(1, value)
            } else {
                inductor(1, value)
            };
            let coefficient = if electric { 4.0 * value } else { 8.0 * value };
            let mut graph = rc(f64::NAN);
            if !electric {
                graph.branches[1] = (0, 1, Branch::Inductor { henries: f64::NAN });
            }
            let circuit = ConductionCircuit::new_with_material_storage(
                problem(
                    &cards,
                    [ElectricalPropertyBasis::Resistivity; 2],
                    300.0,
                    1.0,
                ),
                graph,
                terminal(),
                vec![],
                vec![source.clone()],
                0.0,
                if electric { vec![] } else { vec![4.0] },
                if electric { vec![4.0] } else { vec![] },
                &[],
                field_controls(),
                || false,
            )
            .unwrap();
            let next = layout
                .step_circuit_candidate(
                    &thermal,
                    dt,
                    &circuit,
                    &[],
                    &[],
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            // Independent midpoint RC / RL solution at the actual R=1+3=4 ohm.
            let ratio = if electric {
                (8.0 * coefficient - dt) / (8.0 * coefficient + dt)
            } else {
                (2.0 * coefficient - 4.0 * dt) / (2.0 * coefficient + 4.0 * dt)
            };
            let end = 4.0 * ratio;
            let current = if electric {
                (4.0 - end) / dt
            } else {
                (4.0 + end) / (2.0 * coefficient)
            };
            let stored_before = 8.0 / coefficient;
            let stored_after = end * end / (2.0 * coefficient);
            let heat = dt * 4.0 * current * current;
            let electrical = &next.report().electrical;
            near(
                if electric {
                    next.circuit().charges_c()[0]
                } else {
                    next.circuit().fluxes_wb()[0]
                },
                end,
            );
            near(electrical.spatial_current_a, current);
            near(electrical.spatial_voltage_v, 4.0 * current);
            near(electrical.stored_before_j, stored_before);
            near(electrical.stored_after_j, stored_after);
            near(electrical.heat_j, heat);
            near(stored_before - stored_after, heat);
            near(electrical.native.supplied, 0.0);
            for (i, &power) in electrical.spatial.cell_power_w.iter().enumerate() {
                near(
                    power,
                    0.5 * current * current * if i < 2 { 1.0 } else { 3.0 },
                );
                near(electrical.spatial.current_density_a_m2[i][0], -current);
            }
            // Independent two-cell backward-Euler heat/conduction balance.
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + dt * current * current * (2.0 + 4.0 * dt) / (2.0 + 3.0 * dt),
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + dt * current * current * (3.0 + 4.0 * dt) / (2.0 + 3.0 * dt),
            );
            near(
                next.thermal().energy_j() + stored_after,
                thermal.energy_j() + stored_before,
            );
            assert_eq!(next.circuit().storage(), &[source.clone()]);
            let second = layout
                .step_circuit_candidate(
                    next.thermal(),
                    2.0 * dt,
                    next.circuit(),
                    &[],
                    &[],
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            near(
                if electric {
                    second.circuit().charges_c()[0]
                } else {
                    second.circuit().fluxes_wb()[0]
                },
                end * ratio,
            );
            assert_eq!(second.circuit().storage(), &[source]);
            near(
                second.thermal().energy_j() + second.report().electrical.stored_after_j,
                thermal.energy_j() + stored_before,
            );
            outcomes.push((end, heat));
        }
        assert!(outcomes[1].0 > outcomes[0].0);
        assert!(outcomes[1].1 < outcomes[0].1);
    }
}

#[test]
fn g3_live_storage_reconfiguration_accounts_work_and_preserves_spatial_sources() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 0.0, None);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let bindings = cells(&cards, &curves, &context, &key);
    let mut graph = rc(f64::NAN);
    graph
        .branches
        .push((0, 1, Branch::Inductor { henries: f64::NAN }));
    let circuit = ConductionCircuit::new_with_material_storage(
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        ),
        graph,
        terminal(),
        vec![],
        vec![inductor(2, 0.25), capacitor(1, 0.125)],
        0.0,
        vec![2.0],
        vec![4.0],
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    let change = circuit
        .rebind_storage_candidate(vec![inductor(2, 0.5), capacitor(1, 0.25)], &[], || false)
        .unwrap();
    assert_eq!(change.circuit().charges_c(), circuit.charges_c());
    assert_eq!(change.circuit().fluxes_wb(), circuit.fluxes_wb());
    assert_eq!(change.circuit().time_s(), circuit.time_s());
    assert_eq!(
        change.circuit().problem().conductors(),
        circuit.problem().conductors()
    );
    assert_eq!(
        change.circuit().graph().branches[0],
        circuit.graph().branches[0]
    );
    near(change.parameter_work_j(), -8.5);
    assert_eq!(
        change.work().iter().map(|w| w.branch).collect::<Vec<_>>(),
        vec![1, 2]
    );
    for work in change.work() {
        near(work.stored_after_j, 0.5 * work.stored_before_j);
    }
    let reverse = change
        .circuit()
        .rebind_storage_candidate(circuit.storage().to_vec(), &[], || false)
        .unwrap();
    near(reverse.parameter_work_j(), 8.5);
    assert_eq!(reverse.circuit().storage(), circuit.storage());
    let next = network(&bindings)
        .step_circuit_candidate(
            &thermal,
            0.25,
            change.circuit(),
            &[],
            &[],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    near(next.report().electrical.stored_before_j, 8.5);
    near(
        next.thermal().energy_j() + next.report().electrical.stored_after_j,
        thermal.energy_j() + 17.0 + change.parameter_work_j(),
    );
    near(
        next.thermal().energy_j() - thermal.energy_j(),
        next.report().electrical.heat_j,
    );
    assert_eq!(next.circuit().storage(), change.circuit().storage());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.charges_c(), &[4.0]);
    assert_eq!(circuit.fluxes_wb(), &[2.0]);
}

#[test]
fn g4_storage_refusals_preserve_charge_flux_and_material_sources() {
    let (a, _) = material(1.0, 1.0, 0.0);
    let (b, _) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let circuit = ConductionCircuit::new_with_material_storage(
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        ),
        rc(f64::NAN),
        terminal(),
        vec![],
        vec![capacitor(1, 0.125)],
        0.0,
        vec![],
        vec![4.0],
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    let original = circuit.storage().to_vec();
    for bad in [
        vec![],
        vec![capacitor(1, 0.25), capacitor(1, 0.25)],
        vec![inductor(1, 0.25)],
        vec![capacitor(0, 0.25)],
    ] {
        assert!(
            circuit
                .rebind_storage_candidate(bad, &[], || false)
                .is_err()
        );
    }
    let mut checks = 0;
    circuit
        .rebind_storage_candidate(original.clone(), &[], || {
            checks += 1;
            false
        })
        .unwrap();
    let mut polled = 0;
    assert!(matches!(
        circuit.rebind_storage_candidate(vec![capacitor(1, 0.25)], &[], || {
            polled += 1;
            polled == checks
        }),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(circuit.storage(), &original);
    assert_eq!(circuit.charges_c(), &[4.0]);
    assert!(circuit.fluxes_wb().is_empty());
    assert_eq!(circuit.time_s(), 0.0);

    let mut graph = rc(f64::NAN);
    graph
        .branches
        .push((1, 0, Branch::VoltageSource { port: 0 }));
    let clamped = ConductionCircuit::new_with_material_storage(
        problem(
            &cards,
            [ElectricalPropertyBasis::Resistivity; 2],
            300.0,
            1.0,
        ),
        graph,
        terminal(),
        vec![],
        vec![capacitor(1, 0.125)],
        0.0,
        vec![],
        vec![1.0],
        &[2.0],
        field_controls(),
        || false,
    )
    .unwrap();
    assert!(
        clamped
            .rebind_storage_candidate(vec![capacitor(1, 0.25)], &[2.0], || false)
            .is_err()
    );
    assert_eq!(clamped.charges_c(), &[1.0]);
    assert_eq!(clamped.storage(), &[capacitor(1, 0.125)]);
}
