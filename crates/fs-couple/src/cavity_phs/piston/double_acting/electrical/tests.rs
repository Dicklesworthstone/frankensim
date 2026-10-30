//! Native shared gas/momentum and Kirchhoff owners, with synthetic material data.
use super::*;
use crate::cavity_phs::piston::{
    GasPiston,
    tests::{controls as gas_controls, gas},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor};
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::CircuitGraph;

fn resistor(r: f64) -> MaterialResistor {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: r,
                dims: ELECTRICAL_RESISTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic motion circuit resistivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    MaterialResistor {
        branch: 1,
        conductor: resolve_uniform_conductor(
            &card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn circuit(
    r: f64,
    charge: f64,
    extra: bool,
) -> MaterialCircuit {
    let mut branches = vec![
        (1, 0, Branch::Capacitor { farads: 0.01 }),
        (1, 2, Branch::Resistor { ohms: 999.0 }),
        (
            2,
            if extra { 3 } else { 0 },
            Branch::VoltageSource { port: 0 },
        ),
    ];
    if extra {
        branches.push((3, 0, Branch::VoltageSource { port: 1 }));
    }
    MaterialCircuit::new(
        CircuitGraph {
            node_count: if extra { 4 } else { 3 },
            branches,
            transformers: vec![],
        },
        vec![resistor(r)],
        vec![],
        vec![charge],
        if extra { &[0.0, 0.0] } else { &[0.0] },
    )
    .unwrap()
}
pub(in crate::cavity_phs::piston::double_acting) fn piston(velocity: f64) -> DoubleActingGasPiston {
    DoubleActingGasPiston::new(
        GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, velocity).unwrap(),
        gas(0.75, 0.001),
        0.01,
    )
    .unwrap()
}
pub(in crate::cavity_phs::piston::double_acting) fn input() -> DoublePistonInput {
    DoublePistonInput {
        duration_s: 0.00005,
        heat_j: [0.0; 2],
        external_force_n: 0.0,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn mechanical() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: gas_controls(),
        opposing_gas: gas_controls().gas,
        volume_tolerance_m3: 1e-15,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn solve() -> ElectricalMotionSolve {
    ElectricalMotionSolve {
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-8,
            power_tolerance_w: 1e-8,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
        maximum_iterations: 32,
        maximum_trials: 300,
        difference_velocity_m_s: 1e-4,
        velocity_tolerance_m_s: 1e-9,
        energy_tolerance_j: 1e-7,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn electrical<'a>(
    c: &'a MaterialCircuit,
    inputs: &'a [f64],
    factor: f64,
) -> ElectricalMotionInput<'a> {
    ElectricalMotionInput {
        circuit: c,
        resistors: c.resistors(),
        inputs,
        port: TranslationalCircuitPort {
            branch: 2,
            factor_n_a: factor,
        },
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-7, "{a} != {b}");
}

#[test]
fn finite_charge_drives_motion_and_back_emf_reduces_native_current() {
    let mut c = circuit(2.0, 0.2, false);
    let mut p = piston(0.0);
    let initial = p.clone();
    let saved = c.clone();
    let mut heat = 0.0;
    for index in 0..40 {
        let proposal = p
            .electrical_candidate(
                electrical(&c, &[0.0], 2.0),
                input(),
                mechanical(),
                solve(),
                p.primary().velocity_m_s(),
                || false,
            )
            .unwrap();
        let r = proposal.motion();
        assert!(r.current_a > 0.0 && r.force_n > 0.0);
        if index == 0 {
            assert!(r.back_emf_v > 0.0 && r.mechanical_work_j > 0.0 && r.electrical_work_j < 0.0);
        }
        near(r.force_n, 2.0 * r.current_a);
        near(
            r.back_emf_v,
            2.0 * f64::midpoint(
                proposal.piston_step().velocity_m_s[0],
                proposal.piston_step().velocity_m_s[1],
            ),
        );
        let locked = c
            .step_candidate(
                c.resistors().to_vec(),
                &[0.0],
                input().duration_s,
                solve().circuit,
                || false,
            )
            .unwrap();
        // The gas spring can reverse motion: then emf aids, rather than opposes, discharge.
        near(
            r.current_a - locked.step().resistors[0].current_a,
            -r.back_emf_v / (2.0 + input().duration_s / (2.0 * 0.01)),
        );
        near(proposal.total_energy_residual_j(), 0.0);
        near(r.transfer_residual_j, 0.0);
        let (next_c, next_p, e, _, _, _) = proposal.into_parts();
        heat += e.heat_j;
        c = next_c;
        p = next_p;
    }
    let du = p.primary().gas().mass_kg()
        * (p.primary().gas().state().specific_internal_energy_j_kg()
            - initial
                .primary()
                .gas()
                .state()
                .specific_internal_energy_j_kg())
        + p.opposing().mass_kg()
            * (p.opposing().state().specific_internal_energy_j_kg()
                - initial.opposing().state().specific_internal_energy_j_kg());
    near(
        0.5 * c.charges_c()[0].powi(2) / 0.01
            + heat
            + du
            + 0.5 * p.primary().mass_kg() * p.primary().velocity_m_s().powi(2),
        2.0,
    );
    assert_eq!(saved.charges_c(), &[0.2]);
}

#[test]
fn moving_gas_mechanism_generates_charge_with_a_reaction_force() {
    let c = circuit(2.0, 0.0, false);
    let p = piston(0.1);
    let proposal = p
        .electrical_candidate(
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            0.1,
            || false,
        )
        .unwrap();
    let r = proposal.motion();
    assert!(
        r.current_a < 0.0
            && r.force_n < 0.0
            && r.mechanical_work_j < 0.0
            && r.electrical_work_j > 0.0
    );
    assert!(proposal.circuit().circuit().charges_c()[0] > 0.0);
    assert!(proposal.circuit().step().stored_after_j > 0.0);
    let mut unloaded = p.clone();
    unloaded.step(input(), mechanical(), || false).unwrap();
    assert!(proposal.piston().primary().velocity_m_s() < unloaded.primary().velocity_m_s());
    near(proposal.total_energy_residual_j(), 0.0);
}

#[test]
fn resistance_data_and_winding_orientation_change_the_physical_response() {
    let p = piston(0.0);
    let low = circuit(2.0, 0.2, false);
    let high = circuit(8.0, 0.2, false);
    let run = |c: &MaterialCircuit, k| {
        p.electrical_candidate(
            electrical(c, &[0.0], k),
            input(),
            mechanical(),
            solve(),
            0.0,
            || false,
        )
        .unwrap()
    };
    let a = run(&low, 2.0);
    let b = run(&high, 2.0);
    let reverse = run(&low, -2.0);
    assert!(a.piston().primary().velocity_m_s() > b.piston().primary().velocity_m_s());
    assert!(
        a.piston().primary().gas().state().state().temperature
            < b.piston().primary().gas().state().state().temperature
    );
    assert!(
        a.piston().opposing().state().state().temperature
            > b.piston().opposing().state().state().temperature
    );
    near(
        a.piston().primary().velocity_m_s(),
        -reverse.piston().primary().velocity_m_s(),
    );
    near(
        a.circuit().circuit().charges_c()[0],
        reverse.circuit().circuit().charges_c()[0],
    );
    assert_eq!(
        a.circuit().step().resistors[0]
            .material
            .conductor
            .material(),
        low.resistors()[0].conductor.material()
    );
}

#[test]
fn explicit_zero_coupling_matches_both_original_owners() {
    let p = piston(0.0);
    let c = circuit(2.0, 0.2, false);
    let proposal = p
        .electrical_candidate(
            electrical(&c, &[0.0], 0.0),
            input(),
            mechanical(),
            solve(),
            0.0,
            || false,
        )
        .unwrap();
    let expected = c
        .step_candidate(
            c.resistors().to_vec(),
            &[0.0],
            input().duration_s,
            solve().circuit,
            || false,
        )
        .unwrap();
    assert_eq!(proposal.piston(), &p);
    assert_eq!(proposal.circuit().circuit(), expected.circuit());
    assert_eq!(proposal.motion().force_n, 0.0);
    assert_eq!(proposal.motion().mechanical_work_j, 0.0);
}

#[test]
fn other_voltage_sources_remain_external_not_internal_transfer() {
    let p = piston(0.0);
    let c = circuit(2.0, 0.0, true);
    let proposal = p
        .electrical_candidate(
            electrical(&c, &[0.0, 5.0], 2.0),
            input(),
            mechanical(),
            solve(),
            0.0,
            || false,
        )
        .unwrap();
    let r = proposal.motion();
    let n = &proposal.circuit().step().native;
    near(
        r.external_electrical_work_j,
        input().duration_s * 5.0 * n.y[1],
    );
    assert!(r.external_electrical_work_j > 0.0);
    near(
        r.external_electrical_work_j + r.electrical_work_j,
        n.supplied,
    );
    near(proposal.total_energy_residual_j(), 0.0);
}

#[test]
fn missing_materials_ambiguous_ports_work_limits_and_final_cancel_return_no_state() {
    let p = piston(0.0);
    let c = circuit(2.0, 0.2, false);
    let old_p = p.clone();
    let old_c = c.clone();
    let call = |e, s| p.electrical_candidate(e, input(), mechanical(), s, 0.0, || false);
    assert!(
        call(
            ElectricalMotionInput {
                resistors: &[],
                ..electrical(&c, &[0.0], 2.0)
            },
            solve()
        )
        .is_err()
    );
    assert!(call(electrical(&c, &[1.0], 2.0), solve()).is_err());
    assert!(
        call(
            ElectricalMotionInput {
                port: TranslationalCircuitPort {
                    branch: 1,
                    factor_n_a: 2.0
                },
                ..electrical(&c, &[0.0], 2.0)
            },
            solve()
        )
        .is_err()
    );
    assert!(matches!(
        call(
            electrical(&c, &[0.0], 2.0),
            ElectricalMotionSolve {
                maximum_trials: 1,
                ..solve()
            }
        ),
        Err(ElectricalMotionError::Budget)
    ));
    let mut polls = 0;
    p.electrical_candidate(
        electrical(&c, &[0.0], 2.0),
        input(),
        mechanical(),
        solve(),
        0.0,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        p.electrical_candidate(
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            0.0,
            || {
                count += 1;
                count == polls
            }
        ),
        Err(ElectricalMotionError::Cancelled)
    ));
    assert_eq!(p, old_p);
    assert_eq!(c, old_c);
}
