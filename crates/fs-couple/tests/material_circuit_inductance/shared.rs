//! Common material, all original windings, and native transformer transients.

use super::{material, near};
use fs_couple::material_circuit::inductance::{
    MaterialReversibleInductor, MaterialSharedCoreInductor, ReversibleCoreGeometry,
    ReversibleInductorCircuit, ReversibleInductorGeometry, SharedCoreWinding,
};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};

fn geometry() -> ReversibleCoreGeometry {
    ReversibleCoreGeometry {
        core_area_m2: 2.0,
        core_length_m: 4.0,
        gap_reluctance_per_h: 5.0,
    }
}

fn windings() -> Vec<SharedCoreWinding> {
    vec![
        SharedCoreWinding {
            branch: 2,
            turns: 2.0,
            leakage_inductance_h: 0.25,
        },
        SharedCoreWinding {
            branch: 4,
            turns: -1.0,
            leakage_inductance_h: 0.5,
        },
    ]
}

#[test]
fn one_source_core_is_shared_and_its_energy_is_partitioned_once() {
    let source = material(&[(0.0, 0.0), (1.0, 2.0), (3.0, 14.0)]);
    let core = MaterialSharedCoreInductor::new(windings(), source.clone(), geometry()).unwrap();
    assert_eq!(core.material(), &source);
    assert_eq!(core.geometry(), geometry());
    assert_eq!(core.windings(), windings());
    assert_eq!(
        core.law().core_flux_mmf_knots(),
        &[
            (-6.0, -86.0),
            (-2.0, -18.0),
            (0.0, 0.0),
            (2.0, 18.0),
            (6.0, 86.0)
        ]
    );

    // An additional (2 A, 4 A) contributes zero signed core MMF, but both
    // physical leakage stores remain. Binary-exact turns/leakage keep the
    // source-edge fixture representable without extrapolation by rounding.
    for b in [-3.0_f64, -2.0, -0.25, 0.0, 0.25, 2.0, 3.0] {
        let h = b.signum()
            * if b.abs() <= 1.0 {
                2.0 * b.abs()
            } else {
                6.0 * b.abs() - 4.0
            };
        let mmf = 4.0 * h + 10.0 * b;
        let currents = [mmf / 2.0 + 2.0, 4.0];
        let fluxes = core.flux_linkages_wb(&currents).unwrap();
        let point = core.operating_point(&fluxes).unwrap();
        near(fluxes[0], 4.0 * b + 0.25 * currents[0], 3e-13);
        near(fluxes[1], -2.0 * b + 0.5 * currents[1], 3e-13);
        near(point.common_flux_wb, 2.0 * b, 2e-13);
        near(point.core_flux_density_t, b, 2e-13);
        near(point.core_field_strength_a_m, h, 1e-12);
        for (&actual, &expected) in point.currents_a.iter().zip(&currents) {
            near(actual, expected, 2e-12);
        }
        let density = if b.abs() <= 1.0 {
            b * b
        } else {
            3.0 * b * b - 4.0 * b.abs() + 2.0
        };
        near(point.core_energy_j, 8.0 * density, 3e-11);
        near(point.gap_energy_j, 10.0 * b * b, 3e-11);
        near(
            point.leakage_energy_j[0],
            0.125 * currents[0] * currents[0],
            3e-11,
        );
        near(
            point.leakage_energy_j[1],
            0.25 * currents[1] * currents[1],
            3e-11,
        );
        near(
            point.stored_energy_j,
            point.core_energy_j + point.gap_energy_j + point.leakage_energy_j.iter().sum::<f64>(),
            5e-11,
        );
        let reluctance = if b.abs() < 1.0 { 9.0 } else { 17.0 };
        for (&actual, expected) in point.differential_inductance_h.iter().zip([
            0.25 + 4.0 / reluctance,
            -2.0 / reluctance,
            -2.0 / reluctance,
            0.5 + 1.0 / reluctance,
        ]) {
            near(actual, expected, 1e-13);
        }
    }
}

#[test]
fn one_winding_limit_matches_the_original_material_inductor() {
    let g = geometry();
    let source = material(&[(0.0, 0.0), (1.0, 2.0), (3.0, 14.0)]);
    let scalar = MaterialReversibleInductor::new(
        2,
        source.clone(),
        ReversibleInductorGeometry {
            turns: 2.0,
            core_area_m2: g.core_area_m2,
            core_length_m: g.core_length_m,
            gap_reluctance_per_h: g.gap_reluctance_per_h,
            leakage_inductance_h: 0.25,
        },
    )
    .unwrap();
    let core = MaterialSharedCoreInductor::new(vec![windings()[0]], source, g).unwrap();
    for current in [-25.0, -10.0, -1.0, 0.0, 1.0, 10.0, 25.0] {
        let flux = scalar.flux_linkage_wb(current).unwrap();
        near(core.flux_linkages_wb(&[current]).unwrap()[0], flux, 1e-13);
        let a = scalar.operating_point(flux).unwrap();
        let b = core.operating_point(&[flux]).unwrap();
        near(a.current_a, b.currents_a[0], 5e-13);
        near(a.core_flux_density_t, b.core_flux_density_t, 1e-13);
        near(a.stored_energy_j, b.stored_energy_j, 1e-11);
        near(
            a.differential_inductance_h,
            b.differential_inductance_h[0],
            1e-13,
        );
    }
}

fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 1.0 }),
            (2, 0, Branch::Inductor { henries: 999.0 }),
            (3, 0, Branch::Resistor { ohms: 2.0 }),
            (3, 0, Branch::Inductor { henries: 999.0 }),
        ],
        transformers: vec![],
    }
}

fn transformer_source() -> MaterialSharedCoreInductor {
    MaterialSharedCoreInductor::new(
        vec![
            SharedCoreWinding {
                branch: 2,
                turns: 1.0,
                leakage_inductance_h: 0.1,
            },
            SharedCoreWinding {
                branch: 4,
                turns: 1.0,
                leakage_inductance_h: 0.2,
            },
        ],
        material(&[(0.0, 0.0), (1.0, 1.0), (5.0, 21.0)]),
        ReversibleCoreGeometry {
            core_area_m2: 1.0,
            core_length_m: 1.0,
            gap_reluctance_per_h: 0.0,
        },
    )
    .unwrap()
}

#[test]
fn loaded_transformer_crosses_saturation_with_one_core_store_and_native_heat() {
    let source = transformer_source();
    let circuit = ReversibleInductorCircuit::with_shared_cores(
        graph(),
        vec![],
        vec![source.clone()],
        &[],
        &[],
        &[],
        &[],
        &mut || false,
    )
    .unwrap();
    assert_eq!(circuit.shared_cores(), &[source]);
    assert!(circuit.sources().is_empty());
    assert_eq!(
        circuit.graph().branches[2].2,
        Branch::Inductor { henries: 1.1 }
    );
    assert_eq!(
        circuit.graph().branches[4].2,
        Branch::Inductor { henries: 1.2 }
    );
    let mut state = circuit
        .consistent_initial_state(&[0.0, 0.0], &[], &[3.0])
        .unwrap();
    let (mut work, mut heat, mut secondary_peak) = (0.0, 0.0, 0.0_f64);
    let mut crossed_knee = false;
    for _ in 0..300 {
        let (step, _) = circuit.dae().step_audited(&state, &[3.0], 0.01).unwrap();
        near(step.delta_h + step.dissipated, step.supplied, 1e-10);
        assert!(step.dissipated >= 0.0);
        work += step.supplied;
        heat += step.dissipated;
        state = step.x;
        let points = circuit.shared_operating_points(&state).unwrap();
        let point = &points[0];
        near(point.stored_energy_j + heat, work, 1e-8);
        for (&native, &actual) in circuit
            .dae()
            .inductor_currents(&state)
            .unwrap()
            .iter()
            .zip(&point.currents_a)
        {
            near(native, actual, 2e-12);
        }
        secondary_peak = secondary_peak.max(point.currents_a[1].abs());
        crossed_knee |= point.core_flux_density_t > 1.0;
    }
    assert!(
        crossed_knee,
        "primary bias must reach the supplied common-core saturation knee"
    );
    assert!(
        secondary_peak > 0.1,
        "the loaded secondary must receive real induced current"
    );
    assert!(heat > 0.0);
}

#[test]
fn source_support_leakage_branch_ownership_and_cancellation_refuse() {
    let source = transformer_source();
    assert!(source.flux_linkages_wb(&[100.0, 0.0]).is_err());
    assert!(source.operating_point(&[1.0]).is_err());
    let construct = |windings| {
        MaterialSharedCoreInductor::new(
            windings,
            material(&[(0.0, 0.0), (1.0, 1.0), (5.0, 21.0)]),
            ReversibleCoreGeometry {
                core_area_m2: 1.0,
                core_length_m: 1.0,
                gap_reluctance_per_h: 0.0,
            },
        )
    };
    assert!(construct(vec![]).is_err());
    for invalid_leakage in [0.0, -0.1, f64::NAN] {
        let mut windings = source.windings().to_vec();
        windings[0].leakage_inductance_h = invalid_leakage;
        assert!(construct(windings).is_err());
    }
    let mut duplicate = source.windings().to_vec();
    duplicate[1].branch = duplicate[0].branch;
    assert!(construct(duplicate).is_err());
    assert!(
        ReversibleInductorCircuit::with_shared_cores(
            graph(),
            vec![],
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
        ReversibleInductorCircuit::with_shared_cores(
            graph(),
            vec![],
            vec![source.clone()],
            &[],
            &[],
            &[],
            &[MutualInductance {
                branches: (2, 4),
                henries: 0.1
            }],
            &mut || false,
        )
        .is_err()
    );
    let scalar = MaterialReversibleInductor::new(
        2,
        source.material().clone(),
        ReversibleInductorGeometry {
            turns: 1.0,
            core_area_m2: 1.0,
            core_length_m: 1.0,
            gap_reluctance_per_h: 0.0,
            leakage_inductance_h: 0.1,
        },
    )
    .unwrap();
    assert!(
        ReversibleInductorCircuit::with_shared_cores(
            graph(),
            vec![scalar],
            vec![source.clone()],
            &[],
            &[],
            &[],
            &[],
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        ReversibleInductorCircuit::with_shared_cores(
            graph(),
            vec![],
            vec![source],
            &[],
            &[],
            &[],
            &[],
            &mut || true,
        ),
        Err(fs_couple::material_circuit::MaterialCircuitError::Cancelled)
    ));
}

#[test]
fn material_saturation_drives_the_native_converter_and_its_finite_bus_work() {
    use fs_phs::circuit::converter::{
        AveragedBuckBoost, BridgeModulation, BuckBoostParameters, NonlinearConverterControls,
    };

    let source = MaterialReversibleInductor::new(
        0,
        material(&[(0.0, 0.0), (1.0, 1.0), (3.0, 9.0)]),
        ReversibleInductorGeometry {
            turns: 2.0,
            core_area_m2: 0.5,
            core_length_m: 2.0,
            gap_reluctance_per_h: 3.0,
            leakage_inductance_h: 0.1,
        },
    )
    .unwrap();
    let capacitance = [0.5, 0.25];
    let model = AveragedBuckBoost::new(
        BuckBoostParameters {
            // Consuming source binding must replace this descriptive linear
            // store completely, including physical current initialization.
            inductance_h: 999.0,
            bus_capacitance_f: capacitance,
            winding_resistance_ohm: 0.05,
            input_on_resistance_ohm: [0.02, 0.03],
            output_on_resistance_ohm: [0.04, 0.05],
        },
        BridgeModulation {
            input_high_fraction: 0.8,
            output_high_fraction: 0.6,
        },
    )
    .unwrap()
    .with_inductor_curve(source.curve().clone())
    .unwrap();
    let mut state = model.state(0.2, [5.0, 1.0]).unwrap();
    near(state[0], source.flux_linkage_wb(0.2).unwrap(), 1e-14);
    let initial_magnetic = source.operating_point(state[0]).unwrap().stored_energy_j;
    let initial_total = initial_magnetic + 0.5 * capacitance[0] * 25.0 + 0.5 * capacitance[1];
    let controls = NonlinearConverterControls {
        absolute_energy_tolerance_j: 2e-11,
        relative_balance_tolerance: 2e-10,
        ..Default::default()
    };
    let mut prepared = model.prepare_controlled(0.01, controls, || false).unwrap();
    let mut heat = [0.0; 5];
    let mut crossed_material_knee = false;
    for _ in 0..180 {
        let before = source.operating_point(state[0]).unwrap();
        let mut next = [0.0; 3];
        let mut voltage = [0.0; 2];
        let record = prepared
            .step_controlled(&state, &[0.0; 2], &mut next, &mut voltage, || false)
            .unwrap();
        let after = source.operating_point(next[0]).unwrap();
        near(
            record.magnetic_energy_change_j,
            after.stored_energy_j - before.stored_energy_j,
            2e-11,
        );
        let mut bus_energy = 0.0;
        for p in 0..2 {
            let old_voltage = state[p + 1] / capacitance[p];
            let new_voltage = next[p + 1] / capacitance[p];
            let old_energy = 0.5 * capacitance[p] * old_voltage * old_voltage;
            let new_energy = 0.5 * capacitance[p] * new_voltage * new_voltage;
            near(
                record.bus_energy_change_j[p],
                new_energy - old_energy,
                2e-11,
            );
            near(record.terminal_work_j[p], 0.0, 0.0);
            bus_energy += new_energy;
        }
        for (total, increment) in heat.iter_mut().zip(record.conduction.joules) {
            assert!(increment >= 0.0);
            *total += increment;
        }
        near(
            after.stored_energy_j + bus_energy + heat.iter().sum::<f64>(),
            initial_total,
            2e-9,
        );
        crossed_material_knee |= after.core_flux_density_t.abs() > 1.0;
        state = next;
    }
    assert!(crossed_material_knee);
    assert!(heat.iter().all(|value| *value > 0.0));
}
