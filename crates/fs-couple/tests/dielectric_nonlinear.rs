//! Nonlinear instantaneous displacement with native Debye memory and heat.
//! Source curves/coefficients are authored numerical fixtures, not measurements.
use fs_couple::material_circuit::capacitance::dielectric::MaterialDielectricCapacitor;
use fs_couple::material_circuit::dielectric::DebyeCircuit;
use fs_couple::material_circuit::{MaterialCircuitError, MaterialCircuitSolve};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::field_storage::dielectric::reversible::{
    ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
    ReversibleDielectricMaterial,
};
use fs_material::field_storage::dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor};
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{QuantitySpec, Temperature, Time};

const AREA: f64 = 0.5;
const GAP: f64 = 0.25;
const POLES: [(f64, f64); 2] = [(1.0, 0.2), (0.5, 0.7)]; // chi [F/m], tau [s]

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        250.0,
        450.0,
    )
}
fn names(count: usize) -> Vec<DebyePoleProperties> {
    (0..count)
        .map(|i| DebyePoleProperties {
            permittivity_increment: format!("chi_{i}"),
            relaxation_time: format!("tau_{i}"),
        })
        .collect()
}
fn fixture(
    linear: bool,
    background: f64,
    poles: &[(f64, f64)],
    variable_increment: bool,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, value, interpolation, validity| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored nonlinear Debye fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(
        "epsilon_inf",
        ELECTRIC_PERMITTIVITY_DIMS,
        PropertyValue::Scalar {
            value: background,
            dims: ELECTRIC_PERMITTIVITY_DIMS,
        },
        InterpolationPolicy::ConstantWithinValidity,
        validity(),
    );
    let knots = if linear {
        vec![(0.0, 0.0), (4.0, 8.0)]
    } else {
        vec![(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)]
    };
    add(
        REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ELECTRIC_DISPLACEMENT_DIMS,
        PropertyValue::Curve {
            abscissa: "E".into(),
            abscissa_dims: ELECTRIC_FIELD_DIMS,
            knots,
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        InterpolationPolicy::LinearInside,
        validity().with_quantity(
            "E",
            QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
            0.0,
            4.0,
        ),
    );
    for (i, &(chi, tau)) in poles.iter().enumerate() {
        let value = if variable_increment {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.5 * chi), (450.0, 1.5 * chi)],
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            }
        } else {
            PropertyValue::Scalar {
                value: chi,
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            }
        };
        add(
            &format!("chi_{i}"),
            ELECTRIC_PERMITTIVITY_DIMS,
            value,
            if variable_increment {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            validity(),
        );
        add(
            &format!("tau_{i}"),
            Time::DIMS,
            PropertyValue::Scalar {
                value: tau,
                dims: Time::DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
            validity(),
        );
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "nonlinear instantaneous plus explicit passive poles".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolved(card: &MaterialCard, count: usize, temperature: f64) -> ResolvedDebyeCapacitor {
    ResolvedDebyeCapacitor::resolve(
        card,
        &point(temperature),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names(count),
        AREA,
        GAP,
    )
    .unwrap()
}
fn instantaneous(card: &MaterialCard, temperature: f64, area: f64) -> MaterialDielectricCapacitor {
    let point = point(temperature)
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    let material = ReversibleDielectricMaterial::resolve(
        card,
        &point,
        "E",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MaterialDielectricCapacitor::new(0, material, area, GAP).unwrap()
}
fn graph(extra: &[Branch]) -> CircuitGraph {
    CircuitGraph {
        node_count: 2,
        branches: std::iter::once(Branch::Capacitor { farads: f64::NAN })
            .chain(extra.iter().copied())
            .map(|b| (1, 0, b))
            .collect(),
        transformers: vec![],
    }
}
fn nonlinear(
    card: &MaterialCard,
    graph: CircuitGraph,
    charges: Vec<f64>,
    polarization: Vec<f64>,
    inputs: &[f64],
) -> DebyeCircuit {
    DebyeCircuit::new_with_reversible_instantaneous(
        graph,
        resolved(card, polarization.len(), 300.0),
        instantaneous(card, 300.0, AREA),
        vec![],
        charges,
        polarization,
        inputs,
        || false,
    )
    .unwrap()
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 2e-10,
        power_tolerance_w: 2e-10,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
// Independent integrals and inverse of the original source's D(E) segments.
fn density(d: f64) -> f64 {
    let d = d.abs();
    if d <= 2.0 {
        d * d / 4.0
    } else if d <= 3.0 {
        let x = d - 2.0;
        1.0 + x + x * x
    } else {
        let x = d - 3.0;
        3.0 + 3.0 * x + x * x / 8.0
    }
}
fn voltage(q: f64) -> f64 {
    let d = (q / AREA).abs();
    q.signum()
        * GAP
        * if d <= 2.0 {
            d / 2.0
        } else if d <= 3.0 {
            1.0 + 2.0 * (d - 2.0)
        } else {
            3.0 + (d - 3.0) / 4.0
        }
}
fn energy(q: f64, p: &[f64]) -> f64 {
    AREA * GAP * density(q / AREA)
        + p.iter()
            .zip(POLES)
            .map(|(&p, (chi, _))| 0.5 * p * p / (chi * AREA / GAP))
            .sum::<f64>()
}
fn rk4<const N: usize>(
    mut x: [f64; N],
    duration: f64,
    steps: u32,
    rhs: impl Fn([f64; N]) -> [f64; N],
) -> [f64; N] {
    let dt = duration / f64::from(steps);
    for _ in 0..steps {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = rhs(std::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = rhs(std::array::from_fn(|i| x[i] + dt * c[i]));
        x = std::array::from_fn(|i| x[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}

#[test]
fn linear_and_zero_pole_limits_recover_the_original_native_debye_owner() {
    for poles in [&POLES[..], &[][..]] {
        let card = fixture(true, 2.0, poles, false);
        let polarization = if poles.is_empty() {
            vec![]
        } else {
            vec![0.3, -0.1]
        };
        let extra = [Branch::Resistor { ohms: 0.7 }];
        let mut old = DebyeCircuit::new(
            graph(&extra),
            0,
            resolved(&card, poles.len(), 300.0),
            vec![],
            vec![0.8],
            polarization.clone(),
            &[],
            || false,
        )
        .unwrap();
        let mut new = nonlinear(&card, graph(&extra), vec![0.8], polarization, &[]);
        for _ in 0..80 {
            let a = old
                .step_candidate(&[], 0.002, controls(), || false)
                .unwrap();
            let b = new
                .step_candidate(&[], 0.002, controls(), || false)
                .unwrap();
            near(
                a.circuit().instantaneous_charge_c(),
                b.circuit().instantaneous_charge_c(),
                3e-12,
            );
            for (a, b) in a
                .circuit()
                .polarization_charges_c()
                .into_iter()
                .zip(b.circuit().polarization_charges_c())
            {
                near(a, b, 3e-12);
            }
            near(
                a.step().dielectric_energy_change_j,
                b.step().dielectric_energy_change_j,
                2e-13,
            );
            near(
                a.step().polarization_heat_j,
                b.step().polarization_heat_j,
                2e-13,
            );
            near(
                a.step().other_resistor_heat_j,
                b.step().other_resistor_heat_j,
                2e-13,
            );
            old = a.into_parts().0;
            new = b.into_parts().0;
        }
        assert!(old.reversible_instantaneous().is_none());
        assert_eq!(
            new.reversible_instantaneous()
                .unwrap()
                .material()
                .material()
                .card_identity(),
            card.content_hash()
        );
    }
}

#[test]
fn nonlinear_absorption_recovers_voltage_across_a_source_knot_with_positive_heat() {
    let card = fixture(false, 2.0, &POLES, false);
    let exact = rk4([2.0, 1.0], 1.0, 16000, |p| {
        let v = voltage(3.0 - p[0] - p[1]);
        std::array::from_fn(|i| (POLES[i].0 * AREA / GAP * v - p[i]) / POLES[i].1)
    });
    let run = |sign: f64, steps: u32| {
        let mut circuit = nonlinear(&card, graph(&[]), vec![0.0], vec![2.0 * sign, sign], &[]);
        let accepted = circuit.clone();
        let initial = energy(0.0, &[2.0, 1.0]);
        let mut heat = 0.0;
        for _ in 0..steps {
            let proposal = circuit
                .step_candidate(&[], 1.0 / f64::from(steps), controls(), || false)
                .unwrap();
            let next = proposal.circuit();
            let step = proposal.step();
            assert!(step.polarization.iter().all(|loss| loss.heat_j >= 0.0));
            near(step.dielectric_terminal_work_j, 0.0, 2e-11);
            near(next.terminal_charge_c().unwrap(), 3.0 * sign, 2e-11);
            heat += step.polarization_heat_j;
            near(
                next.dielectric_energy_j().unwrap(),
                energy(
                    next.instantaneous_charge_c(),
                    &next.polarization_charges_c(),
                ),
                2e-12,
            );
            near(next.dielectric_energy_j().unwrap() + heat, initial, 3e-10);
            near(
                step.dielectric_energy_residual_j,
                0.0,
                controls().energy_tolerance_j,
            );
            assert!(
                std::ptr::eq(next.dae(), circuit.dae()),
                "ordinary steps reuse the prepared native descriptor"
            );
            circuit = proposal.into_parts().0;
        }
        assert!(
            circuit.instantaneous_charge_c().abs() > 1.0,
            "instantaneous response must cross the original knot"
        );
        assert!(heat > 0.0);
        assert_eq!(accepted.instantaneous_charge_c(), 0.0);
        assert_eq!(accepted.polarization_charges_c(), vec![2.0 * sign, sign]);
        circuit
            .polarization_charges_c()
            .into_iter()
            .zip(exact)
            .map(|(a, b)| (a - sign * b).abs())
            .fold(0.0_f64, f64::max)
    };
    let coarse = run(1.0, 40);
    let fine = run(1.0, 80);
    assert!(fine < coarse / 2.8, "coarse {coarse:e}, fine {fine:e}");
    assert!(
        fine < 4e-5,
        "independent nonlinear absorption error {fine:e}"
    );
    near(run(-1.0, 80), fine, 1e-10);
}

#[test]
fn held_voltage_retains_nonlinear_instantaneous_energy_and_analytic_debye_heat() {
    let card = fixture(false, 2.0, &POLES, false);
    let v = 0.6;
    let instant = instantaneous(&card, 300.0, AREA);
    let q = instant.charge_c(v).unwrap();
    let mut circuit = nonlinear(
        &card,
        graph(&[Branch::VoltageSource { port: 0 }]),
        vec![q],
        vec![0.0; 2],
        &[v],
    );
    let h0 = AREA * GAP * density(q / AREA);
    near(circuit.dielectric_energy_j().unwrap(), h0, 2e-15);
    assert!((h0 - 0.5 * q * v).abs() > 0.1);
    let (duration, dt) = (0.1, 0.0005);
    let mut heat = [0.0; 2];
    let mut supplied = 0.0;
    for _ in 0..200 {
        let proposal = circuit
            .step_candidate(&[v], dt, controls(), || false)
            .unwrap();
        for (i, loss) in proposal.step().polarization.iter().enumerate() {
            heat[i] += loss.heat_j;
            near(
                loss.heat_j,
                dt * loss.current_a.powi(2) * circuit.source().poles()[i].resistance_ohm,
                2e-14,
            );
        }
        supplied += proposal.step().dielectric_terminal_work_j;
        near(proposal.circuit().instantaneous_charge_c(), q, 2e-12);
        near(proposal.circuit().terminal_voltage_v(), v, 2e-12);
        circuit = proposal.into_parts().0;
    }
    for (i, (chi, tau)) in POLES.into_iter().enumerate() {
        let c = chi * AREA / GAP;
        let decay = (-duration / tau).exp();
        near(
            circuit.polarization_charges_c()[i],
            c * v * (1.0 - decay),
            3e-7,
        );
        near(heat[i], 0.5 * c * v * v * (1.0 - decay * decay), 2e-7);
    }
    near(
        circuit.dielectric_energy_j().unwrap() - h0 + heat.iter().sum::<f64>(),
        supplied,
        2e-10,
    );
}

#[test]
fn driven_warm_memory_reports_bulk_and_polarization_heat_below_a_large_stationary_bias() {
    let card = fixture(false, 2.0, &POLES, false);
    let graph = CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: f64::NAN }),
            (2, 1, Branch::Resistor { ohms: 0.15 }),
            (2, 0, Branch::VoltageSource { port: 0 }),
            (3, 0, Branch::Capacitor { farads: 1.0 }),
        ],
        transformers: vec![],
    };
    let mut circuit = nonlinear(&card, graph, vec![0.2, 1e8], vec![0.3, -0.1], &[0.7]);
    let reference = rk4([0.2, 0.3, -0.1], 1.0, 16000, |x| {
        let v = voltage(x[0]);
        let p0 = (2.0 * v - x[1]) / 0.2;
        let p1 = (v - x[2]) / 0.7;
        [(0.7 - v) / 0.15 - p0 - p1, p0, p1]
    });
    let h0 = energy(0.2, &[0.3, -0.1]);
    let mut heat = 0.0;
    let mut work = 0.0;
    let bias = circuit.dae().charge_index[1];
    for _ in 0..200 {
        let proposal = circuit
            .step_candidate(&[0.7], 0.005, controls(), || false)
            .unwrap();
        let step = proposal.step();
        assert_eq!(step.other_resistors.len(), 1);
        assert_eq!(step.other_resistors[0].branch, 1);
        near(
            step.other_resistor_heat_j,
            0.005 * 0.15 * step.other_resistors[0].current_a.powi(2),
            2e-13,
        );
        assert!(step.polarization_heat_j >= 0.0 && step.other_resistor_heat_j > 0.0);
        near(step.native.delta_h, step.dielectric_energy_change_j, 2e-10);
        assert_eq!(
            proposal.circuit().state()[bias].to_bits(),
            1e8_f64.to_bits()
        );
        heat += step.polarization_heat_j + step.other_resistor_heat_j;
        work += step.native.supplied;
        near(
            proposal.circuit().dielectric_energy_j().unwrap() + heat,
            h0 + work,
            4e-9,
        );
        circuit = proposal.into_parts().0;
    }
    let actual = [
        circuit.instantaneous_charge_c(),
        circuit.polarization_charges_c()[0],
        circuit.polarization_charges_c()[1],
    ];
    for (a, b) in actual.into_iter().zip(reference) {
        near(a, b, 3e-5);
    }
    assert!(actual[0] > 1.0);
}

#[test]
fn unresolved_native_partition_work_refuses_instead_of_becoming_polarization_heat() {
    let card = fixture(false, 2.0, &POLES[..1], false);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: f64::NAN }),
            (2, 0, Branch::Capacitor { farads: 1.0 }),
            (2, 0, Branch::Resistor { ohms: 1.0 }),
        ],
        transformers: vec![],
    };
    let circuit = nonlinear(&card, graph, vec![0.9, 1.0], vec![2.0], &[]);
    let before = circuit.state().to_vec();
    let initial = circuit
        .dae()
        .consistent_initial_state(circuit.state(), &[])
        .unwrap();
    let (native, _) = circuit.dae().step_audited(&initial, &[], 0.05).unwrap();
    // The unrelated quadratic store has a known exact energy. A global
    // Gonzalez correction may distribute the nonlinear remainder to its effort;
    // the native whole-circuit balance alone cannot certify this local transfer.
    let q_other = circuit.dae().charge_index[1];
    let node = circuit.dae().node_potential_index[1];
    let midpoint = 0.5 * (initial[node] + native.x[node]);
    let local_defect =
        0.5 * (native.x[q_other] - initial[q_other]) * (native.x[q_other] + initial[q_other])
            + 0.05 * midpoint * midpoint;
    assert!(
        local_defect.abs() > 1e-8,
        "fixture must expose actual partition work: {local_defect:e}"
    );
    assert!(matches!(
        circuit.step_candidate(&[], 0.05, controls(), || false),
        Err(MaterialCircuitError::Balance(_))
    ));
    assert_eq!(circuit.state(), before);
    assert!(
        circuit
            .step_candidate(&[], 1e-4, controls(), || false)
            .is_ok()
    );
}

#[test]
fn source_background_context_refresh_and_clamp_refusals_preserve_all_memory() {
    let card = fixture(false, 2.0, &POLES, false);
    let circuit = nonlinear(&card, graph(&[]), vec![0.2], vec![0.3, -0.1], &[]);
    let before = circuit.state().to_vec();
    assert!(matches!(
        circuit.step_candidate(&[], 0.01, controls(), || true),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert!(
        circuit
            .rebind_source_candidate(circuit.source().clone(), &[], || false)
            .is_err()
    );
    assert!(
        circuit
            .rebind_source_and_resistors_candidate(circuit.source().clone(), vec![], &[], || false)
            .is_err()
    );
    assert_eq!(circuit.state(), before);
    let construct = |source, instant, charges, ports| {
        DebyeCircuit::new_with_reversible_instantaneous(
            graph(&[Branch::VoltageSource { port: 0 }]),
            source,
            instant,
            vec![],
            charges,
            vec![0.3, -0.1],
            ports,
            || false,
        )
    };
    assert!(
        construct(
            resolved(&card, 2, 300.0),
            instantaneous(&card, 350.0, AREA),
            vec![0.2],
            &[0.05]
        )
        .is_err()
    );
    assert!(
        construct(
            resolved(&card, 2, 300.0),
            instantaneous(&card, 300.0, 0.4),
            vec![0.2],
            &[0.05]
        )
        .is_err()
    );
    let other = fixture(false, 1.0, &POLES, false);
    assert!(
        construct(
            resolved(&other, 2, 300.0),
            instantaneous(&other, 300.0, AREA),
            vec![0.2],
            &[0.05]
        )
        .is_err(),
        "same-card nonlinear tangent may not contradict declared epsilon_inf"
    );
    assert!(
        construct(
            resolved(&card, 2, 300.0),
            instantaneous(&other, 300.0, AREA),
            vec![0.2],
            &[0.05]
        )
        .is_err()
    );
    let variable = fixture(false, 2.0, &POLES, true);
    assert!(
        construct(
            resolved(&variable, 2, 300.0),
            instantaneous(&variable, 300.0, AREA),
            vec![0.2],
            &[0.05]
        )
        .is_err()
    );
    assert!(
        construct(
            resolved(&card, 2, 300.0),
            instantaneous(&card, 300.0, AREA),
            vec![0.2],
            &[0.6]
        )
        .is_err(),
        "an ideal source cannot project the original instantaneous charge"
    );
    assert!(
        construct(
            resolved(&card, 2, 300.0),
            instantaneous(&card, 300.0, AREA),
            vec![3.6],
            &[1.1]
        )
        .is_err()
    );
}
