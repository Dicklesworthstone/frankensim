//! Original sourced play coefficients, one common memory, and actual native ports.
use super::*;
use crate::material_circuit::hysteretic::shared::{
    SharedHystereticCircuit, SharedHystereticCircuitInitial,
};
use fs_material::field_storage::hysteresis::shared::{
    SharedHystereticCore, SharedHystereticWinding,
};
use fs_phs::circuit::MutualInductance;

const RATIOS: [f64; 2] = [2.0, -1.0];
const LEAKAGE: [f64; 2] = [0.3, 0.4];
const RESISTANCE: [f64; 2] = [0.7, 1.1];
const PORTS: [f64; 3] = [5.0, 17.0, -3.0]; // Preserve the authored gap at port 1.

fn shared_core(pinning: f64, reverse: bool) -> SharedHystereticCore {
    let (card, _) = source::core_material(pinning, 310.0);
    SharedHystereticCore::new(
        source::sourced_core(&card),
        vec![
            SharedHystereticWinding {
                turns_ratio: RATIOS[0],
                leakage_h: LEAKAGE[0],
            },
            SharedHystereticWinding {
                turns_ratio: if reverse { -RATIOS[1] } else { RATIOS[1] },
                leakage_h: LEAKAGE[1],
            },
        ],
    )
    .unwrap()
}
fn graph(reverse: bool) -> CircuitGraph {
    CircuitGraph {
        node_count: 5,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: f64::NAN }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
            (3, 0, Branch::VoltageSource { port: 2 }),
            (3, 4, Branch::Resistor { ohms: f64::NAN }),
            if reverse {
                (0, 4, Branch::Inductor { henries: f64::NAN })
            } else {
                (4, 0, Branch::Inductor { henries: f64::NAN })
            },
        ],
        transformers: vec![],
    }
}
fn resistors() -> Vec<MaterialResistor> {
    RESISTANCE
        .into_iter()
        .zip([1, 4])
        .map(|(r, branch)| resistor(branch, &material(r, 0.0).0))
        .collect()
}
fn initial(core: &SharedHystereticCore) -> SharedHystereticCircuitInitial {
    SharedHystereticCircuitInitial {
        magnetic: core.initial_state(),
        remaining_fluxes_wb: vec![],
        charges_c: vec![],
        external_inputs: vec![0.0; 3],
        time_s: 0.0,
    }
}
fn circuit(pinning: f64, reverse: bool) -> SharedHystereticCircuit {
    let core = shared_core(pinning, reverse);
    let start = initial(&core);
    SharedHystereticCircuit::new(graph(reverse), vec![2, 5], core, resistors(), start).unwrap()
}

#[test]
fn linear_shared_limit_matches_full_native_mutual_storage_with_unrelated_rlc() {
    let mut graph = graph(false);
    let mut sources = resistors();
    graph.node_count += 1;
    graph.branches.extend([
        (5, 0, Branch::Inductor { henries: 0.2 }),
        (1, 5, Branch::Resistor { ohms: 0.4 }),
        (5, 0, Branch::Capacitor { farads: 0.3 }),
    ]);
    sources.push(resistor(7, &material(0.4, 0.0).0));
    let original = shared_core(0.0, false);
    // Deliberately retain the caller's nonascending branch/winding order.
    let shared = SharedHystereticCore::new(
        original.core().clone(),
        vec![original.windings()[1], original.windings()[0]],
    )
    .unwrap();
    let mut start = initial(&shared);
    start.remaining_fluxes_wb = vec![0.0];
    start.charges_c = vec![0.0];
    let mut history =
        SharedHystereticCircuit::new(graph.clone(), vec![5, 2], shared, sources.clone(), start)
            .unwrap();
    let permeability = 0.01 + 0.2; // Zero pinning, unsaturated native source limit.
    for (winding, branch) in [2, 5].into_iter().enumerate() {
        graph.branches[branch].2 = Branch::Inductor {
            henries: LEAKAGE[winding] + RATIOS[winding].powi(2) * permeability,
        };
    }
    let mut native = MaterialCircuit::new_with_mutual_inductance(
        graph,
        sources.clone(),
        vec![MutualInductance {
            branches: (2, 5),
            henries: RATIOS[0] * RATIOS[1] * permeability,
        }],
        vec![0.0; 3],
        vec![0.0],
        &[0.0; 3],
    )
    .unwrap();
    let ports = [0.8, 17.0, -0.4];
    for tick in 1..=10 {
        let a = native
            .step_candidate(sources.clone(), &ports, 0.005, controls().circuit, || false)
            .unwrap();
        let b = history
            .step_candidate(
                f64::from(tick) * 0.005,
                sources.clone(),
                &ports,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(b.circuit().branch_indices(), &[5, 2]);
        let linkage = b.circuit().magnetic_state().flux_linkages_wb();
        near(linkage[0], a.circuit().fluxes_wb()[1], 2e-9);
        near(linkage[1], a.circuit().fluxes_wb()[0], 2e-9);
        near(
            b.circuit().native_circuit().fluxes_wb()[0],
            a.circuit().fluxes_wb()[2],
            2e-9,
        );
        near(
            b.circuit().native_circuit().charges_c()[0],
            a.circuit().charges_c()[0],
            2e-9,
        );
        near(b.step().supplied_work_j, a.step().native.supplied, 2e-9);
        near(b.step().heat_j, a.step().heat_j, 2e-9);
        near(b.step().magnetic.core.core_heat_j, 0.0, 0.0);
        near(b.step().stored_after_j, a.step().stored_after_j, 2e-9);
        assert_eq!(
            b.step()
                .circuit
                .resistors
                .iter()
                .map(|heat| heat.material.branch)
                .collect::<Vec<_>>(),
            vec![1, 4, 7]
        );
        native = a.into_parts().0;
        history = b.into_parts().0;
    }
}

// Independent first-monotone-segment play formula and analytic elimination of
// each native RL branch. Only a scalar bracketed root remains in this oracle.
fn first_step(pinning: f64, dt: f64) -> ([f64; 2], [f64; 2], f64, f64, f64) {
    let evaluate = |field: f64| {
        let polarization = (0.2 * (field - pinning)).clamp(0.0, 1.0);
        let flux = 0.01 * field + polarization;
        let energy = 0.005 * field * field + polarization * polarization / 0.4;
        let heat = pinning * polarization;
        let loading = RATIOS[0].powi(2) / LEAKAGE[0] + RATIOS[1].powi(2) / LEAKAGE[1];
        let drive = field + loading * flux;
        let mean_flux = if drive == 0.0 {
            0.0
        } else {
            0.5 * flux + (0.5 * field * flux - energy - heat) / drive
        };
        let linkages = std::array::from_fn(|i| {
            let voltage = PORTS[2 * i];
            (dt * voltage + dt * RESISTANCE[i] * RATIOS[i] * mean_flux / LEAKAGE[i])
                / (1.0 + dt * RESISTANCE[i] / (2.0 * LEAKAGE[i]))
        });
        let residual =
            RATIOS[0] * linkages[0] / LEAKAGE[0] + RATIOS[1] * linkages[1] / LEAKAGE[1] - drive;
        let means =
            std::array::from_fn(|i| (0.5 * linkages[i] - RATIOS[i] * mean_flux) / LEAKAGE[i]);
        (residual, linkages, means, flux, field, heat)
    };
    let (mut low, mut high) = (0.0, 20.0);
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if evaluate(middle).0 > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let (_, linkages, means, flux, field, heat) = evaluate(0.5 * (low + high));
    (linkages, means, flux, field, heat)
}

#[test]
fn two_driven_windings_cross_pinning_with_one_real_core_heat_and_exact_chord_work() {
    let mut results = Vec::new();
    for pinning in [1.0, 2.0] {
        let model = circuit(pinning, false);
        let saved = model.clone();
        let next = model
            .step_candidate(0.1, resistors(), &PORTS, controls(), || false)
            .unwrap();
        let (linkages, means, flux, field, core_heat) = first_step(pinning, 0.1);
        assert!(field > pinning && core_heat > 0.0);
        near(next.circuit().reference_flux_linkage_wb(), flux, 3e-9);
        near(
            next.circuit().magnetic_state().core_state().field_a_m(),
            field,
            3e-9,
        );
        near(next.step().magnetic.core.core_heat_j, core_heat, 3e-9);
        let mut winding_heat = 0.0;
        for i in 0..2 {
            near(
                next.circuit().magnetic_state().flux_linkages_wb()[i],
                linkages[i],
                3e-9,
            );
            near(
                next.circuit().magnetic_state().currents_a()[i],
                (linkages[i] - RATIOS[i] * flux) / LEAKAGE[i],
                3e-9,
            );
            near(next.step().magnetic.mean_currents_a[i], means[i], 3e-9);
            near(
                next.step().terminal_voltages_v[i],
                PORTS[2 * i] - RESISTANCE[i] * means[i],
                3e-9,
            );
            winding_heat += 0.1 * RESISTANCE[i] * means[i] * means[i];
        }
        let endpoint = next.circuit().magnetic_state().currents_a();
        near(
            RATIOS[0] * endpoint[0] + RATIOS[1] * endpoint[1],
            field,
            3e-9,
        );
        // An endpoint-average current is wrong across the first pinning knee.
        assert!((0.5 * endpoint[0] - means[0]).abs() > 1e-4);
        near(next.step().circuit.heat_j, winding_heat, 3e-9);
        near(next.step().heat_j, winding_heat + core_heat, 3e-9);
        near(next.step().energy_residual_j, 0.0, 3e-9);
        near(
            next.step().supplied_work_j,
            0.1 * (PORTS[0] * means[0] + PORTS[2] * means[1]),
            3e-9,
        );
        assert_eq!(model, saved);
        let warmed = next
            .circuit()
            .with_core_source_query(&point(305.0))
            .unwrap();
        assert_eq!(warmed.magnetic_state(), next.circuit().magnetic_state());
        assert_eq!(
            warmed.core().core().source().unwrap().material().query(),
            point(305.0)
        );
        assert!(
            next.circuit()
                .with_core_source_query(&point(311.0))
                .is_err()
        );
        results.push((field, core_heat, winding_heat));
    }
    assert!((results[0].0 - results[1].0).abs() > 0.1);
    assert!((results[0].1 - results[1].1).abs() > 1e-3);
    assert!((results[0].2 - results[1].2).abs() > 1e-3);
}

fn voltage_driven(core: SharedHystereticCore) -> SharedHystereticCircuit {
    let start = initial(&core);
    SharedHystereticCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: f64::NAN }),
                (2, 0, Branch::VoltageSource { port: 2 }),
                (2, 0, Branch::Inductor { henries: f64::NAN }),
            ],
            transformers: vec![],
        },
        vec![1, 3],
        core,
        vec![],
        start,
    )
    .unwrap()
}

#[test]
fn shared_major_loop_retains_remanence_and_deposits_pinning_area_only_once() {
    let mut model = voltage_driven(shared_core(2.0, false));
    let mut polarization = 0.0_f64;
    let (mut cycle_heat, mut cycle_work) = (0.0, 0.0);
    let mut saturated = None;
    for (tick, field) in [8.0_f64, 0.0, -8.0, 0.0, 8.0].into_iter().enumerate() {
        let next_z = polarization
            .clamp(0.2 * (field - 2.0), 0.2 * (field + 2.0))
            .clamp(-1.0, 1.0);
        let flux = 0.01 * field + next_z;
        // A physical endpoint with I0=H/2 and I1=0, satisfying sum(r I)=H.
        let target = [
            LEAKAGE[0] * field / 2.0 + RATIOS[0] * flux,
            RATIOS[1] * flux,
        ];
        let old = model.magnetic_state().flux_linkages_wb();
        let ports = [target[0] - old[0], 17.0, target[1] - old[1]];
        let next = model
            .step_candidate((tick + 1) as f64, vec![], &ports, controls(), || false)
            .unwrap();
        let material = next.circuit().magnetic_state().core_state();
        near(material.field_a_m(), field, 3e-8);
        near(material.flux_density_t(), flux, 3e-9);
        near(material.polarizations_t()[0], next_z, 3e-9);
        near(
            next.step().magnetic.core.core_heat_j,
            2.0 * (next_z - polarization).abs(),
            3e-9,
        );
        near(next.step().circuit.heat_j, 0.0, 0.0);
        near(next.step().numerical_dissipation_j, 0.0, 0.0);
        near(
            next.circuit().magnetic_state().currents_a()[0],
            field / 2.0,
            3e-8,
        );
        near(next.circuit().magnetic_state().currents_a()[1], 0.0, 3e-8);
        near(next.step().energy_residual_j, 0.0, 3e-8);
        if tick == 0 {
            saturated = Some(next.circuit().magnetic_state().clone());
        } else {
            cycle_heat += next.step().magnetic.core.core_heat_j;
            cycle_work += next.step().supplied_work_j;
        }
        if field == 0.0 {
            near(flux.abs(), 0.4, 1e-14);
            assert!(material.stored_energy_j_m3() > 0.0);
        }
        polarization = next_z;
        model = next.into_parts().0;
    }
    near(cycle_heat, 4.0 * 2.0 * 1.0, 3e-8);
    near(cycle_work, cycle_heat, 3e-8);
    near(
        model.magnetic_state().stored_energy_j(),
        saturated.unwrap().stored_energy_j(),
        3e-8,
    );
}

#[test]
fn winding_orientation_reverses_terminal_coordinates_without_changing_shared_history() {
    let forward = circuit(2.0, false);
    let reverse = circuit(2.0, true);
    let a = forward
        .step_candidate(0.1, resistors(), &PORTS, controls(), || false)
        .unwrap();
    let b = reverse
        .step_candidate(0.1, resistors(), &PORTS, controls(), || false)
        .unwrap();
    near(
        a.circuit().reference_flux_linkage_wb(),
        b.circuit().reference_flux_linkage_wb(),
        3e-9,
    );
    near(
        a.step().magnetic.core.core_heat_j,
        b.step().magnetic.core.core_heat_j,
        3e-9,
    );
    near(a.step().heat_j, b.step().heat_j, 3e-9);
    near(a.step().supplied_work_j, b.step().supplied_work_j, 3e-9);
    for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
        near(
            a.circuit().magnetic_state().flux_linkages_wb()[i],
            sign * b.circuit().magnetic_state().flux_linkages_wb()[i],
            3e-9,
        );
        near(
            a.step().terminal_voltages_v[i],
            sign * b.step().terminal_voltages_v[i],
            3e-9,
        );
        near(
            a.step().magnetic.mean_currents_a[i],
            sign * b.step().magnetic.mean_currents_a[i],
            3e-9,
        );
    }
}

#[test]
fn shared_turns_domains_budgets_and_late_cancellation_are_transactional() {
    let model = circuit(2.0, false);
    let saved = model.clone();
    let mut polls = 0;
    let first = model
        .step_candidate(0.1, resistors(), &PORTS, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let repeated = model
        .step_candidate(0.1, resistors(), &PORTS, controls(), || false)
        .unwrap();
    assert_eq!(first.circuit(), repeated.circuit());
    let mut calls = 0;
    assert!(matches!(
        model.step_candidate(0.1, resistors(), &PORTS, controls(), || {
            calls += 1;
            calls == polls
        }),
        Err(HystereticCircuitError::Cancelled)
    ));
    let mut limited = controls();
    limited.maximum_circuit_trials = 4;
    assert!(matches!(
        model.step_candidate(0.1, resistors(), &PORTS, limited, || false),
        Err(HystereticCircuitError::TrialBudget)
    ));
    assert_eq!(model, saved);

    let linear = circuit(0.0, false);
    let before = linear.clone();
    assert!(matches!(
        linear.step_candidate(10.0, resistors(), &[0.2, 17.0, -0.1], controls(), || false),
        Err(HystereticCircuitError::SharedMagneticReversal { .. })
    ));
    assert_eq!(linear, before);
    let unsupported = voltage_driven(shared_core(2.0, false));
    let before = unsupported.clone();
    assert!(
        unsupported
            .step_candidate(1.0, vec![], &[1000.0, 17.0, -1000.0], controls(), || false)
            .is_err()
    );
    assert_eq!(unsupported, before);

    let core = shared_core(2.0, false);
    assert!(
        SharedHystereticCircuit::new(
            graph(false),
            vec![2, 2],
            core.clone(),
            resistors(),
            initial(&core)
        )
        .is_err()
    );
    assert!(
        SharedHystereticCircuit::new(
            graph(false),
            vec![2, 4],
            core.clone(),
            resistors(),
            initial(&core)
        )
        .is_err()
    );
}
