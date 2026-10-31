//! G1/G3/G4: physical dielectric memory uses native circuit history and heat.
use fs_couple::material_circuit::{
    MaterialCircuitError, MaterialCircuitSolve, dielectric::DebyeCircuit,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor},
    },
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_qty::Time;

fn source(c_inf: f64, poles: &[(f64, f64)]) -> ResolvedDebyeCapacitor {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, value| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic Debye circuit fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let geometry = 10.0;
    add("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, c_inf / geometry);
    let mut names = Vec::new();
    for (index, &(capacitance, time)) in poles.iter().enumerate() {
        let permittivity_increment = format!("delta_{index}");
        let relaxation_time = format!("tau_{index}");
        add(
            &permittivity_increment,
            ELECTRIC_PERMITTIVITY_DIMS,
            capacitance / geometry,
        );
        add(&relaxation_time, Time::DIMS, time);
        names.push(DebyePoleProperties {
            permittivity_increment,
            relaxation_time,
        });
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedDebyeCapacitor::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names,
        0.01,
        0.001,
    )
    .unwrap()
}
fn graph(extra: &[Branch]) -> CircuitGraph {
    CircuitGraph {
        node_count: 2,
        branches: std::iter::once(Branch::Capacitor { farads: f64::NAN })
            .chain(extra.iter().copied())
            .map(|branch| (1, 0, branch))
            .collect(),
        transformers: vec![],
    }
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-9,
        power_tolerance_w: 1e-9,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}

#[test]
fn held_voltage_relaxes_multiple_source_poles_and_returns_physical_heat() {
    let v = 2.0;
    let poles = [(0.3, 0.1), (0.4, 0.7)];
    let material = source(0.2, &poles);
    let mut circuit = DebyeCircuit::new(
        graph(&[Branch::VoltageSource { port: 0 }]),
        0,
        material.clone(),
        vec![],
        vec![0.2 * v],
        vec![0.0; 2],
        &[v],
        || false,
    )
    .unwrap();
    let initial = circuit.dielectric_energy_j().unwrap();
    let mut total_heat = [0.0; 2];
    let mut total_work = 0.0;
    let dt = 0.0005;
    for _ in 0..400 {
        let proposal = circuit
            .step_candidate(&[v], dt, controls(), || false)
            .unwrap();
        assert_eq!(proposal.step().other_resistor_heat_j, 0.0);
        for (i, loss) in proposal.step().polarization.iter().enumerate() {
            assert_eq!(loss.pole, i);
            assert!(loss.heat_j >= 0.0);
            total_heat[i] += loss.heat_j;
        }
        near(
            proposal.step().native.dissipated,
            proposal.step().polarization_heat_j,
            1e-12,
        );
        total_work += proposal.step().dielectric_terminal_work_j;
        circuit = proposal.into_parts().0;
    }
    assert_eq!(circuit.source(), &material);
    near(circuit.terminal_voltage_v(), v, 1e-12);
    for (i, &(c, tau)) in poles.iter().enumerate() {
        let exponential = (-0.2_f64 / tau).exp();
        near(
            circuit.polarization_charges_c()[i],
            c * v * (1.0 - exponential),
            4e-7,
        );
        near(
            total_heat[i],
            0.5 * c * v * v * (1.0 - exponential * exponential),
            4e-7,
        );
    }
    near(
        circuit.dielectric_energy_j().unwrap() - initial + total_heat.iter().sum::<f64>(),
        total_work,
        1e-10,
    );
}

#[test]
fn open_circuit_absorption_recovers_voltage_preserving_total_charge() {
    let (c_inf, c_pole, tau, precharge) = (0.2, 0.3, 0.1, 2.0);
    let initial_polarization = c_pole * precharge;
    let mut circuit = DebyeCircuit::new(
        graph(&[]),
        0,
        source(c_inf, &[(c_pole, tau)]),
        vec![],
        vec![0.0],
        vec![initial_polarization],
        &[],
        || false,
    )
    .unwrap();
    let initial_energy = circuit.dielectric_energy_j().unwrap();
    let mut heat = 0.0;
    let duration = 0.08;
    let steps = 256;
    for _ in 0..steps {
        let proposal = circuit
            .step_candidate(&[], duration / f64::from(steps), controls(), || false)
            .unwrap();
        assert_eq!(proposal.step().native.supplied, 0.0);
        near(proposal.step().dielectric_terminal_work_j, 0.0, 1e-12);
        heat += proposal.step().polarization_heat_j;
        circuit = proposal.into_parts().0;
        near(
            circuit.instantaneous_charge_c() + circuit.polarization_charges_c()[0],
            initial_polarization,
            2e-11,
        );
    }
    // Independent two-capacitor charge conservation yields a DIFFERENT decay
    // time than a voltage-clamped Debye pole: tau_effective=tau*C_inf/(C_inf+C_k).
    let effective_tau = tau * c_inf / (c_inf + c_pole);
    let expected =
        initial_polarization / (c_inf + c_pole) * (1.0 - (-duration / effective_tau).exp());
    near(circuit.terminal_voltage_v(), expected, 2e-6);
    assert!(circuit.terminal_voltage_v() > 1.0 && heat > 0.0);
    near(
        initial_energy,
        circuit.dielectric_energy_j().unwrap() + heat,
        1e-10,
    );
}

#[test]
fn temporal_refinement_resolves_the_same_analytic_debye_step() {
    let run = |steps| {
        let mut circuit = DebyeCircuit::new(
            graph(&[Branch::VoltageSource { port: 0 }]),
            0,
            source(0.2, &[(0.3, 0.1)]),
            vec![],
            vec![0.2],
            vec![0.0],
            &[1.0],
            || false,
        )
        .unwrap();
        for _ in 0..steps {
            circuit = circuit
                .step_candidate(&[1.0], 0.1 / f64::from(steps), controls(), || false)
                .unwrap()
                .into_parts()
                .0;
        }
        (circuit.polarization_charges_c()[0] - 0.3 * (1.0 - (-1.0_f64).exp())).abs()
    };
    let coarse = run(32);
    let fine = run(64);
    assert!(
        coarse / fine > 3.8 && coarse / fine < 4.2,
        "coarse={coarse}, fine={fine}"
    );
}

#[test]
fn native_current_driven_ac_matches_material_admittance_and_loss() {
    let omega = 2.0;
    let amplitude = 0.7;
    let source = source(0.2, &[(0.3, 0.4), (0.1, 0.07)]);
    let y = source.admittance(omega).unwrap();
    let denominator = y.conductance_s * y.conductance_s + y.susceptance_s * y.susceptance_s;
    let vr = amplitude * y.conductance_s / denominator;
    let vi = -amplitude * y.susceptance_s / denominator;
    let charges: Vec<_> = source
        .poles()
        .iter()
        .map(|pole| {
            let x = omega * pole.relaxation_time_s;
            pole.capacitance_f * (vr + vi * x) / (1.0 + x * x)
        })
        .collect();
    let mut circuit = DebyeCircuit::new(
        graph(&[Branch::CurrentSource { port: 0 }]),
        0,
        source.clone(),
        vec![],
        vec![source.instantaneous_capacitance_f() * vr],
        charges,
        &[amplitude],
        || false,
    )
    .unwrap();
    let period = 2.0 * core::f64::consts::PI / omega;
    let steps = 512;
    let dt = period / f64::from(steps);
    let mut heat = 0.0;
    for i in 0..steps {
        let t_mid = (f64::from(i) + 0.5) * dt;
        let proposal = circuit
            .step_candidate(&[amplitude * (omega * t_mid).cos()], dt, controls(), || {
                false
            })
            .unwrap();
        heat += proposal.step().polarization_heat_j;
        circuit = proposal.into_parts().0;
        let phase = omega * f64::from(i + 1) * dt;
        near(
            circuit.terminal_voltage_v(),
            vr * phase.cos() - vi * phase.sin(),
            3e-5,
        );
    }
    let continuous_heat = 0.5 * period * y.conductance_s * amplitude * amplitude / denominator;
    near(heat, continuous_heat, 3e-5);
}

#[test]
fn lossless_limit_and_external_rlc_reuse_the_native_solver() {
    let mut original = graph(&[
        Branch::Resistor { ohms: 3.0 },
        Branch::Inductor { henries: 0.5 },
    ]);
    original.branches[0].2 = Branch::Capacitor { farads: 0.2 };
    let dae = assemble_circuit(&original).unwrap();
    let mut x = vec![0.0; dae.system.state_dim()];
    x[dae.flux_index[0]] = 0.4;
    x[dae.charge_index[0]] = 0.2;
    let x = dae.consistent_initial_state(&x, &[]).unwrap();
    let circuit = DebyeCircuit::new(
        original,
        0,
        source(0.2, &[]),
        vec![0.4],
        vec![0.2],
        vec![],
        &[],
        || false,
    )
    .unwrap();
    for (&a, &b) in circuit.state().iter().zip(&x) {
        near(a, b, 1e-11);
    }
    // The candidate first rechecks consistency at the held inputs. Native
    // initialization preserves charge/flux but reflects the algebraic endpoint
    // about its physical effort, so both native paths must start in that same
    // endpoint convention before comparing the COMPLETE descriptor state.
    let initial = dae.consistent_initial_state(&x, &[]).unwrap();
    for &slot in dae.charge_index.iter().chain(&dae.flux_index) {
        assert_eq!(initial[slot].to_bits(), x[slot].to_bits());
    }
    let expected = dae.step_audited(&initial, &[], 0.01).unwrap().0;
    let actual = circuit
        .step_candidate(&[], 0.01, controls(), || false)
        .unwrap();
    assert_eq!(actual.step().polarization_heat_j, 0.0);
    assert!(actual.step().other_resistor_heat_j > 0.0);
    for (&a, &b) in actual.step().native.x.iter().zip(&expected.x) {
        near(a, b, 1e-11);
    }
    near(
        actual.step().other_resistor_heat_j,
        expected.dissipated,
        1e-12,
    );
    near(actual.step().native.delta_h, expected.delta_h, 1e-12);
    near(actual.step().native.supplied, expected.supplied, 1e-12);
}

#[test]
fn polarization_and_external_load_heat_are_separate_positive_accounts() {
    let circuit = DebyeCircuit::new(
        graph(&[Branch::Resistor { ohms: 2.0 }]),
        0,
        source(0.2, &[(0.3, 0.1)]),
        vec![],
        vec![0.2],
        vec![0.0],
        &[],
        || false,
    )
    .unwrap();
    let proposal = circuit
        .step_candidate(&[], 0.01, controls(), || false)
        .unwrap();
    assert!(proposal.step().polarization_heat_j > 0.0);
    assert!(proposal.step().other_resistor_heat_j > 0.0);
    near(
        proposal.step().polarization_heat_j + proposal.step().other_resistor_heat_j,
        proposal.step().native.dissipated,
        1e-12,
    );
}

#[test]
fn invalid_history_ideal_impulses_and_cancellation_preserve_accepted_state() {
    assert!(
        DebyeCircuit::new(
            graph(&[]),
            0,
            source(0.2, &[(0.3, 0.1)]),
            vec![],
            vec![0.0],
            vec![],
            &[],
            || false
        )
        .is_err()
    );
    assert!(
        DebyeCircuit::new(
            graph(&[Branch::VoltageSource { port: 0 }]),
            0,
            source(0.2, &[(0.3, 0.1)]),
            vec![],
            vec![0.0],
            vec![0.0],
            &[1.0],
            || false
        )
        .is_err()
    );
    let circuit = DebyeCircuit::new(
        graph(&[Branch::VoltageSource { port: 0 }]),
        0,
        source(0.2, &[(0.3, 0.1)]),
        vec![],
        vec![0.2],
        vec![0.0],
        &[1.0],
        || false,
    )
    .unwrap();
    let before = circuit.state().to_vec();
    assert!(
        circuit
            .step_candidate(&[2.0], 0.1, controls(), || false)
            .is_err()
    );
    assert!(matches!(
        circuit.step_candidate(&[1.0], 0.1, controls(), || true),
        Err(MaterialCircuitError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        circuit.step_candidate(&[1.0], 0.1, controls(), || {
            polls += 1;
            polls == 4
        }),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert_eq!(circuit.state(), before);
    assert!(
        circuit
            .step_candidate(&[1.0], -0.1, controls(), || false)
            .is_err()
    );
}

#[test]
fn small_polarization_and_inductor_work_survive_large_dielectric_bias() {
    // The fixed instantaneous energy is 5e15 J, whose representable spacing is
    // larger than ALL work in this interval. Independent physical reservoirs
    // still exchange charge, flux and heat; endpoint subtraction loses it.
    let bias = 1e8;
    let material = source(1.0, &[(1e-16, 0.1)]);
    let polarization_capacitance = material.poles()[0].capacitance_f;
    let circuit = DebyeCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
                (2, 0, Branch::Inductor { henries: 1.0 }),
                (2, 0, Branch::Resistor { ohms: 1.0 }),
            ],
            transformers: vec![],
        },
        0,
        material,
        vec![1.0],
        vec![bias],
        vec![0.0],
        &[bias],
        || false,
    )
    .unwrap();
    let initial_h = circuit.dae().system.hamiltonian(circuit.state());
    let proposal = circuit
        .step_candidate(&[bias], 0.01, controls(), || false)
        .unwrap();
    let accepted = proposal.circuit();
    let step = proposal.step();
    assert_eq!(step.dielectric_energy_j[0], step.dielectric_energy_j[1]);
    assert_eq!(
        accepted.dae().system.hamiltonian(accepted.state()),
        initial_h
    );
    assert_eq!(
        accepted.instantaneous_charge_c(),
        circuit.instantaneous_charge_c()
    );

    // The initially empty polarization capacitor gives an independent small
    // endpoint-energy oracle, without cancellation from the held reservoir.
    let charge = accepted.polarization_charges_c()[0];
    let polarization_energy = 0.5 * charge * charge / polarization_capacitance;
    assert!(charge > 0.0 && polarization_energy > 0.0);
    near(step.dielectric_energy_change_j, polarization_energy, 1e-12);
    near(
        step.dielectric_terminal_work_j,
        polarization_energy + step.polarization_heat_j,
        1e-12,
    );
    assert!(step.polarization_heat_j > 0.0);

    // The separate RL loop starts with 0.5 J, so its ordinary endpoint values
    // resolve the magnetic release independently of the large electric bias.
    let flux = accepted.state()[accepted.dae().flux_index[0]];
    let magnetic_release = 0.5 - 0.5 * flux * flux;
    assert!(magnetic_release > 0.0);
    near(step.other_resistor_heat_j, magnetic_release, 1e-12);
    near(
        step.native.delta_h,
        polarization_energy - magnetic_release,
        1e-12,
    );
    assert!(step.dielectric_energy_residual_j.abs() < 1e-12);
    assert!(step.circuit_energy_residual_j.abs() < 1e-12);
}
