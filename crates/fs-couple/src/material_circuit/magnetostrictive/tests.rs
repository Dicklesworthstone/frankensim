//! Independent physical reductions with authored synthetic material claims.
use super::*;
use crate::magnetostrictive::{MagnetostrictiveLimits, MagnetostrictiveRod};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor};
use fs_material::field_storage::MAGNETIC_PERMEABILITY_DIMS;
use fs_material::magnetostrictive::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims, QuantitySpec};

fn card(roles: &[(&str, Dims, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, dims, value) in roles {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic reciprocal circuit regression".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored magnetomechanical fixture".into(),
            phase: "solid".into(),
            process: "fixed reversible tangent".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn rod(d: f64, density: f64, turns: u32) -> MagnetostrictiveRod {
    let card = card(&[
        (COMPLIANCE_33_H_PROPERTY, COMPLIANCE_DIMS, 0.1),
        (PIEZOMAGNETIC_33_PROPERTY, PIEZOMAGNETIC_DIMS, d),
        (PERMEABILITY_33_T_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, 2.0),
        ("density", Density::DIMS, density),
    ]);
    let point = QueryPoint::new()
        .with_quantity(
            "bias_H",
            QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS),
            10.0,
        )
        .unwrap();
    let law = LongitudinalMagnetostrictive::resolve(
        &card,
        &point,
        "bias_H",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    MagnetostrictiveRod::new(
        law,
        2.0,
        0.5,
        turns,
        MagnetostrictiveLimits {
            maximum_abs_strain: 0.05,
            maximum_abs_field_increment_a_m: 0.5,
        },
    )
    .unwrap()
}
fn resistor(branch: usize, resistance: f64) -> MaterialResistor {
    let card = card(&[(
        "electrical_resistivity",
        ELECTRICAL_RESISTIVITY_DIMS,
        resistance,
    )]);
    MaterialResistor {
        branch,
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
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-10,
        power_tolerance_w: 1e-11,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance {tolerance:.3e}"
    );
}
fn driven(
    rod: MagnetostrictiveRod,
    r: f64,
    x: f64,
    v: f64,
    current: f64,
    reverse: bool,
) -> MagnetostrictiveCircuit {
    let linkage = rod.clamped_inductance_h() * current + rod.coupling_wb_m() * x;
    MagnetostrictiveCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: 999.0 }),
                (
                    if reverse { 0 } else { 2 },
                    if reverse { 2 } else { 0 },
                    Branch::Inductor { henries: 999.0 },
                ),
            ],
            transformers: vec![],
        },
        vec![resistor(1, r)],
        vec![linkage],
        vec![],
        2,
        rod,
        x,
        v,
        &[0.0],
        || false,
    )
    .unwrap()
}
fn passive(
    rod: MagnetostrictiveRod,
    r: f64,
    x: f64,
    v: f64,
    current: f64,
) -> MagnetostrictiveCircuit {
    let linkage = rod.clamped_inductance_h() * current + rod.coupling_wb_m() * x;
    MagnetostrictiveCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Resistor { ohms: 999.0 }),
                (1, 0, Branch::Inductor { henries: 999.0 }),
            ],
            transformers: vec![],
        },
        vec![resistor(0, r)],
        vec![linkage],
        vec![],
        1,
        rod,
        x,
        v,
        &[],
        || false,
    )
    .unwrap()
}

#[test]
fn native_rl_rod_step_matches_independent_reciprocal_force_and_motional_voltage() {
    for (d, rho, n, r) in [
        (0.2, 3.0, 3, 0.8),
        (-0.2, 3.0, 3, 0.8),
        (0.0, 3.0, 3, 0.8),
        (0.2, 6.0, 3, 0.8),
        (0.2, 3.0, 6, 0.8),
        (0.2, 3.0, 3, 4.0),
    ] {
        let (x0, v0, i0, force, voltage, dt) = (0.002, -0.003, 0.004, 0.03, 0.05, 0.2);
        let model = driven(rod(d, rho, n), r, x0, v0, i0, false);
        let saved = model.clone();
        let candidate = model
            .step_candidate(
                model.circuit().resistors().to_vec(),
                &[voltage],
                force,
                dt,
                controls(),
                || false,
            )
            .unwrap();
        // Independent two-equation midpoint reduction of
        // m v'=F-k_H x+theta i; L_c i'=V-R i-theta v.
        let mass = rho / 3.0;
        let kh = 2.5;
        let lc = (2.0 - d * d / 0.1) * f64::from(n).powi(2) * 0.25;
        let theta = f64::from(n) * 0.25 * d / 0.1;
        let h = 0.5 * dt;
        let a = 1.0 + h * h * kh / mass;
        let b = lc + h * r;
        let rhs_v = v0 + h * (force - kh * x0) / mass;
        let rhs_i = lc * i0 + h * voltage;
        let determinant = a * b + h * h * theta * theta / mass;
        let midpoint_v = (rhs_v * b + h * theta * rhs_i / mass) / determinant;
        let midpoint_i = (a * rhs_i - h * theta * rhs_v) / determinant;
        let x1 = x0 + dt * midpoint_v;
        let v1 = 2.0 * midpoint_v - v0;
        let i1 = 2.0 * midpoint_i - i0;
        near(candidate.model().extension_m(), x1, 3e-12);
        near(candidate.model().velocity_m_s(), v1, 3e-12);
        near(candidate.model().current_increment_a().unwrap(), i1, 3e-12);
        near(
            candidate.model().flux_linkage_increment_wb(),
            lc * i1 + theta * x1,
            3e-12,
        );
        let report = candidate.step();
        near(
            report.incremental_electrical_supply_j,
            dt * voltage * midpoint_i,
            3e-12,
        );
        near(
            report.incremental_mechanical_supply_j,
            dt * force * midpoint_v,
            3e-12,
        );
        near(
            report.incremental_resistor_dissipation_j,
            dt * r * midpoint_i * midpoint_i,
            3e-12,
        );
        near(report.incremental_resistors[0].current_a, midpoint_i, 3e-12);
        near(report.incremental_energy_residual_j, 0.0, 3e-12);
        near(
            report.stored_after_j,
            0.5 * mass * v1 * v1 + 0.5 * kh * x1 * x1 + 0.5 * lc * i1 * i1,
            3e-12,
        );
        assert_eq!(model, saved);
        assert_eq!(candidate.model().rod().material(), model.rod().material());
    }
}

#[test]
fn zero_coupling_preserves_interleaved_native_rlc_branches_and_source_port_gaps() {
    let rod = rod(0.0, 3.0, 3);
    let mut graph = CircuitGraph {
        node_count: 4,
        branches: vec![
            (3, 0, Branch::Capacitor { farads: 0.4 }),
            (2, 3, Branch::Inductor { henries: 0.7 }),
            (1, 2, Branch::Resistor { ohms: 999.0 }),
            (1, 0, Branch::VoltageSource { port: 2 }),
            (2, 0, Branch::Capacitor { farads: 0.3 }),
            (3, 0, Branch::Inductor { henries: 999.0 }),
        ],
        transformers: vec![],
    };
    let mut coupled = MagnetostrictiveCircuit::new(
        graph.clone(),
        vec![resistor(2, 0.8)],
        vec![0.0; 2],
        vec![0.0; 2],
        5,
        rod.clone(),
        0.0,
        0.0,
        &[0.0; 3],
        || false,
    )
    .unwrap();
    graph.branches[5].2 = Branch::Inductor {
        henries: rod.clamped_inductance_h(),
    };
    let mut original = MaterialCircuit::new(
        graph,
        vec![resistor(2, 0.8)],
        vec![0.0; 2],
        vec![0.0; 2],
        &[0.0; 3],
    )
    .unwrap();
    for _ in 0..16 {
        let p = coupled
            .step_candidate(
                coupled.circuit().resistors().to_vec(),
                &[0.0, 0.0, 0.03],
                0.0,
                0.02,
                controls(),
                || false,
            )
            .unwrap();
        let q = original
            .step_candidate(
                original.resistors().to_vec(),
                &[0.0, 0.0, 0.03],
                0.02,
                controls(),
                || false,
            )
            .unwrap();
        for (&a, &b) in p
            .model()
            .circuit()
            .fluxes_wb()
            .iter()
            .zip(q.circuit().fluxes_wb())
        {
            near(a, b, 3e-12);
        }
        for (&a, &b) in p
            .model()
            .circuit()
            .charges_c()
            .iter()
            .zip(q.circuit().charges_c())
        {
            near(a, b, 3e-12);
        }
        near(
            p.step().incremental_resistor_dissipation_j,
            q.step().heat_j,
            3e-12,
        );
        near(
            p.step().incremental_electrical_supply_j,
            q.step().native.supplied,
            3e-12,
        );
        near(p.step().stored_after_j, q.step().stored_after_j, 3e-12);
        assert_eq!(p.model().extension_m(), 0.0);
        assert_eq!(p.model().velocity_m_s(), 0.0);
        assert_eq!(p.step().incremental_resistors[0].material.branch, 2);
        coupled = p.into_parts().0;
        original = q.into_parts().0;
    }
    assert!(coupled.circuit().fluxes_wb().iter().all(|v| v.abs() > 1e-8));
    assert!(coupled.circuit().charges_c().iter().all(|v| v.abs() > 1e-8));
}

#[test]
fn mechanical_vibration_drives_real_resistive_loading_without_an_extra_damping_fit() {
    let mut lossy = passive(rod(0.2, 3.0, 3), 1.0, 0.01, 0.0, 0.0);
    let mut uncoupled = passive(rod(0.0, 3.0, 3), 1.0, 0.01, 0.0, 0.0);
    let initial = lossy.rod_incremental_energy_j().unwrap();
    let mut dissipated = 0.0;
    let mut maximum_current = 0.0_f64;
    for _ in 0..100 {
        let p = lossy
            .step_candidate(
                lossy.circuit().resistors().to_vec(),
                &[],
                0.0,
                0.05,
                controls(),
                || false,
            )
            .unwrap();
        let s = p.step();
        assert_eq!(s.incremental_electrical_supply_j, 0.0);
        assert_eq!(s.incremental_mechanical_supply_j, 0.0);
        assert!(s.incremental_resistor_dissipation_j >= 0.0);
        near(
            s.incremental_resistor_dissipation_j,
            0.05 * s.incremental_resistors[0].current_a.powi(2),
            3e-12,
        );
        dissipated += s.incremental_resistor_dissipation_j;
        near(s.stored_after_j + dissipated, initial, 3e-12);
        maximum_current = maximum_current.max(p.model().current_increment_a().unwrap().abs());
        lossy = p.into_parts().0;
        let q = uncoupled
            .step_candidate(
                uncoupled.circuit().resistors().to_vec(),
                &[],
                0.0,
                0.05,
                controls(),
                || false,
            )
            .unwrap();
        near(q.step().incremental_resistor_dissipation_j, 0.0, 1e-22);
        uncoupled = q.into_parts().0;
    }
    assert!(maximum_current > 1e-4);
    assert!(dissipated > 0.01 * initial);
    assert!(lossy.rod_incremental_energy_j().unwrap() < 0.99 * initial);
    near(
        uncoupled.rod_incremental_energy_j().unwrap(),
        initial,
        3e-12,
    );
}

#[test]
fn coupled_damped_eigenmode_converges_to_independent_continuous_exponential() {
    // Pick a physical R that gives exactly s=-1 in
    // (m*s^2+k_H)*(L_c*s+R)+theta^2*s=0. This independent
    // eigenmode exercises reciprocal mechanics AND resistance simultaneously.
    let (mass, kh, lc, theta, s) = (1.0, 2.5, 3.6, 1.5, -1.0_f64);
    let resistance = -lc * s - theta * theta * s / (mass * s * s + kh);
    let x0 = 0.001;
    let v0 = s * x0;
    let i0 = (mass * s * s + kh) * x0 / theta;
    let lambda0 = lc * i0 + theta * x0;
    let mut previous = None;
    for count in [8, 16, 32] {
        let mut model = passive(rod(0.2, 3.0, 3), resistance, x0, v0, i0);
        let dt = 1.0 / f64::from(count);
        for _ in 0..count {
            model = model
                .step_candidate(
                    model.circuit().resistors().to_vec(),
                    &[],
                    0.0,
                    dt,
                    controls(),
                    || false,
                )
                .unwrap()
                .into_parts()
                .0;
        }
        let factor = s.exp();
        let errors = [
            (model.extension_m() - x0 * factor).abs(),
            (model.velocity_m_s() - v0 * factor).abs(),
            (model.current_increment_a().unwrap() - i0 * factor).abs(),
            (model.flux_linkage_increment_wb() - lambda0 * factor).abs(),
        ];
        if let Some(old) = previous {
            let old: [f64; 4] = old;
            for (coarse, fine) in old.into_iter().zip(errors) {
                assert!(
                    (3.9..4.1).contains(&(coarse / fine)),
                    "ratio {}",
                    coarse / fine
                );
            }
        }
        previous = Some(errors);
    }
    let errors = previous.unwrap();
    assert!(errors[0] < 3.1e-8 && errors[1] < 3.1e-8);
    assert!(errors[2] < 7.1e-8 && errors[3] < 3.1e-7);
}

#[test]
fn reversing_coil_terminals_reverses_actuation_and_sensing_but_preserves_loss() {
    let mut a = driven(rod(0.2, 3.0, 3), 0.8, 0.0, 0.0, 0.0, false);
    let mut b = driven(rod(0.2, 3.0, 3), 0.8, 0.0, 0.0, 0.0, true);
    for _ in 0..12 {
        let p = a
            .step_candidate(
                a.circuit().resistors().to_vec(),
                &[0.05],
                0.0,
                0.05,
                controls(),
                || false,
            )
            .unwrap();
        let q = b
            .step_candidate(
                b.circuit().resistors().to_vec(),
                &[0.05],
                0.0,
                0.05,
                controls(),
                || false,
            )
            .unwrap();
        near(p.model().extension_m(), -q.model().extension_m(), 3e-12);
        near(
            p.model().current_increment_a().unwrap(),
            -q.model().current_increment_a().unwrap(),
            3e-12,
        );
        near(
            p.model().flux_linkage_increment_wb(),
            -q.model().flux_linkage_increment_wb(),
            3e-12,
        );
        near(
            p.step().incremental_electrical_supply_j,
            q.step().incremental_electrical_supply_j,
            3e-12,
        );
        near(
            p.step().incremental_resistor_dissipation_j,
            q.step().incremental_resistor_dissipation_j,
            3e-12,
        );
        a = p.into_parts().0;
        b = q.into_parts().0;
    }
    assert!(a.extension_m() > 1e-5);
}

#[test]
fn true_coupled_current_constraints_and_refusals_preserve_original_history() {
    let rod = rod(0.2, 3.0, 3);
    let x0 = 0.01;
    let linkage = rod.coupling_wb_m() * x0;
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::CurrentSource { port: 0 }),
            (1, 0, Branch::Inductor { henries: 999.0 }),
        ],
        transformers: vec![],
    };
    let model = MagnetostrictiveCircuit::new(
        graph.clone(),
        vec![],
        vec![linkage],
        vec![],
        1,
        rod.clone(),
        x0,
        0.0,
        &[0.0],
        || false,
    )
    .unwrap();
    near(model.current_increment_a().unwrap(), 0.0, 1e-14);
    let dt = 0.1;
    let p = model
        .step_candidate(vec![], &[0.0], 0.0, dt, controls(), || false)
        .unwrap();
    let a = 2.5 * dt * dt / 4.0;
    near(p.model().extension_m(), x0 * (1.0 - a) / (1.0 + a), 3e-12);
    near(p.model().current_increment_a().unwrap(), 0.0, 3e-12);
    near(
        p.model().rod_incremental_energy_j().unwrap(),
        model.rod_incremental_energy_j().unwrap(),
        3e-12,
    );
    assert!(
        MagnetostrictiveCircuit::new(
            graph,
            vec![],
            vec![linkage + 0.1],
            vec![],
            1,
            rod,
            x0,
            0.0,
            &[0.0],
            || false
        )
        .is_err()
    );
    let saved = model.clone();
    assert!(
        model
            .step_candidate(vec![], &[0.1], 0.0, dt, controls(), || false)
            .is_err()
    );
    assert!(
        model
            .step_candidate(vec![], &[0.0], 1e9, dt, controls(), || false)
            .is_err()
    );
    assert!(
        model
            .step_candidate(vec![], &[f64::NAN], 0.0, dt, controls(), || false)
            .is_err()
    );
    assert_eq!(model, saved);
    let mut polls = 0;
    let accepted = model
        .step_candidate(vec![], &[0.0], 0.0, dt, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        model.step_candidate(vec![], &[0.0], 0.0, dt, controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(MagnetostrictiveCircuitError::Circuit(
            MaterialCircuitError::Cancelled
        ))
    ));
    assert_eq!(model, saved);
    let retry = model
        .step_candidate(vec![], &[0.0], 0.0, dt, controls(), || false)
        .unwrap();
    assert_eq!(retry.model(), accepted.model());
}
