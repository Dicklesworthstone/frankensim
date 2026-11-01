//! G1/G3/G4: original Maxwell history drives ABA motion and finite layer heat.
use super::*;
use fs_material::section::relaxation::LayerRelaxationProperties;
use fs_mbd_material::thermal::heating::articulated_motion::torsion::ShaftRelaxation;
use fs_qty::Time;

const DELTA: &str = "relaxing_shear_modulus";
const TAU: &str = "shear_relaxation_time";

fn memory_card(
    base: MaterialCard,
    delta: f64,
    tau: f64,
    slopes: [f64; 2],
    upper: f64,
) -> MaterialCard {
    let mut claims = base.claims().clone();
    for (name, dims, value, slope) in [
        (DELTA, Pressure::DIMS, delta, slopes[0]),
        (TAU, Time::DIMS, tau, slopes[1]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if slope == 0.0 {
                    PropertyValue::Scalar { value, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (upper, value + slope * (upper - 300.0))],
                        dims,
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, upper)
                    .with("omega", 0.0, 200.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if slope == 0.0 {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic Maxwell shaft, not measurements".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn definition<'a>(
    cards: &[&'a MaterialCard],
    twists: Vec<f64>,
) -> (TorsionalJointShaft<'a>, ShaftRelaxation) {
    let extra = [
        ScalarPropertyRequirement::try_new(DELTA, Pressure::DIMS, ScalarAdmissibility::NonNegative)
            .unwrap(),
        ScalarPropertyRequirement::try_new(TAU, Time::DIMS, ScalarAdmissibility::StrictlyPositive)
            .unwrap(),
    ];
    (
        shaft_with_requirements(cards, 300.0, &extra),
        ShaftRelaxation {
            properties: cards
                .iter()
                .map(|_| LayerRelaxationProperties {
                    equilibrium_modulus: G.into(),
                    branches: vec![(DELTA.into(), TAU.into())],
                })
                .collect(),
            initial_spring_twists_rad: twists,
            absolute_spring_twist_tolerance_rad: 1e-11,
        },
    )
}

/// Independent solution of I*q''+c*q'+k*q+a*z=0, z'=q'-z/tau.
/// The supplied positive decay rates are the three roots of its characteristic
/// polynomial. No integration or material-owner helper is used in this oracle.
fn modes(
    q0: f64,
    v0: f64,
    z0: f64,
    k_i: f64,
    c_i: f64,
    a_i: f64,
    rates: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let acceleration = -k_i * q0 - c_i * v0 - a_i * z0;
    let mut q = [0.0; 3];
    let mut z = [0.0; 3];
    for i in 0..3 {
        let j = (i + 1) % 3;
        let l = (i + 2) % 3;
        q[i] = (acceleration + (rates[j] + rates[l]) * v0 + rates[j] * rates[l] * q0)
            / ((rates[i] - rates[j]) * (rates[i] - rates[l]));
        z[i] = -(rates[i] * rates[i] - c_i * rates[i] + k_i) * q[i] / a_i;
    }
    (q, z)
}
fn sample(coefficients: &[f64; 3], rates: &[f64; 3], time: f64) -> f64 {
    coefficients
        .iter()
        .zip(rates)
        .map(|(&c, &r)| c * (-r * time).exp())
        .sum()
}
fn square_integral(coefficients: &[f64; 3], rates: &[f64; 3], time: f64) -> f64 {
    let mut value = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let rate = rates[i] + rates[j];
            value += coefficients[i] * coefficients[j] * -(-rate * time).exp_m1() / rate;
        }
    }
    value
}

#[test]
fn original_layer_spectra_drive_motion_memory_and_independent_physical_heat() {
    let mut outcomes = Vec::new();
    for (g, etas, deltas, tau, rates) in [
        (8.0, [3.0, 11.0], [0.2, 0.6], 2.0 / 3.0, [1.0, 2.0, 3.0]),
        (19.2, [7.0, 15.0], [0.15, 0.35], 0.4, [2.0, 3.0, 4.0]),
    ] {
        let cards = [0, 1].map(|i| {
            memory_card(
                card(g, etas[i], [0.0; 2], 500.0, 10.0 + i as f64 * 10.0, false),
                deltas[i],
                tau,
                [0.0; 2],
                500.0,
            )
        });
        let refs = cards.each_ref();
        let curves = cards.each_ref().map(caloric);
        let reference = mechanism(&refs);
        let (q0, v0, z0, dt) = (0.05, 0.01, 0.03, 0.4);
        let old = motion(&reference, &refs, &curves.each_ref(), &[0.0; 2], q0, v0)
            .with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![z0; 2])], || Ok(()))
            .unwrap();
        let polar = [1.0 / core::f64::consts::TAU, 3.0 / core::f64::consts::TAU];
        let inertia = 2.0 * (polar[0] + polar[1]);
        let k = g * (polar[0] + polar[1]);
        let c = etas[0] * polar[0] + etas[1] * polar[1];
        let a = deltas[0] * polar[0] + deltas[1] * polar[1];
        let (qm, zm) = modes(q0, v0, z0, k / inertia, c / inertia, a / inertia, rates);
        let vm = [0, 1, 2].map(|i| -rates[i] * qm[i]);
        let q = sample(&qm, &rates, dt);
        let v = sample(&vm, &rates, dt);
        let z = sample(&zm, &rates, dt);
        let integral_v2 = square_integral(&vm, &rates, dt);
        let integral_z2 = square_integral(&zm, &rates, dt);
        let (next, audit) = old.advance_to(dt, step(), || Ok(())).unwrap();
        near(next.joint_position()[0], q, 3e-9);
        near(next.joint_velocity()[0], v, 3e-9);
        let state = next.torsional_joint_shafts().next().unwrap();
        let memory = state.relaxation.as_ref().unwrap();
        for (i, &twist) in memory.spring_twists_rad.iter().enumerate() {
            near(twist, z, 3e-9);
            assert_eq!(memory.branches[i].layer_index, i);
            assert_eq!(memory.branches[i].relaxation_time_s, tau);
            let heat = polar[i] * (etas[i] * integral_v2 + deltas[i] / tau * integral_z2);
            let cell = &next.thermal().cells()[&CellId(i as u64 + 1)];
            near(cell.mass_kg * cell.specific_h_j_kg, heat, 3e-9);
            assert_eq!(
                state.section.shear().elastic().layers()[i]
                    .material
                    .card_identity(),
                cards[i].content_hash()
            );
        }
        near(
            audit.viscous_heat_j,
            c * integral_v2 + a / tau * integral_z2,
            3e-9,
        );
        near(
            state.elastic_energy_j,
            0.5 * k * q * q + 0.5 * a * z * z,
            3e-9,
        );
        near(state.generalized_torque_nm, -k * q - c * v - a * z, 3e-9);
        near(state.dissipated_power_w, c * v * v + a / tau * z * z, 3e-9);
        near(audit.parameter_work_j, 0.0, 1e-14);
        near(audit.combined_residual_j, 0.0, 3e-9);
        let (middle, _) = old.advance_to(dt / 2.0, step(), || Ok(())).unwrap();
        let (resumed, _) = middle.advance_to(dt, step(), || Ok(())).unwrap();
        near(resumed.joint_position()[0], q, 3e-9);
        near(
            resumed.thermal().energy_j(),
            next.thermal().energy_j(),
            3e-9,
        );
        near(
            resumed
                .torsional_joint_shafts()
                .next()
                .unwrap()
                .relaxation
                .as_ref()
                .unwrap()
                .spring_twists_rad[0],
            z,
            3e-9,
        );
        assert_eq!(
            old.torsional_joint_shafts()
                .next()
                .unwrap()
                .relaxation
                .as_ref()
                .unwrap()
                .spring_twists_rad,
            vec![z0; 2]
        );
        outcomes.push(q);
    }
    assert!((outcomes[0] - outcomes[1]).abs() > 1e-3);
}

#[test]
fn finite_heat_refreshes_original_relaxation_coefficients_without_resetting_memory() {
    let card = memory_card(
        card(8.0, 9.0, [0.0; 2], 500.0, 10.0, false),
        0.5,
        2.0 / 3.0,
        [0.2, 0.02],
        500.0,
    );
    let curve = caloric(&card);
    let refs = [&card];
    let reference = mechanism(&refs);
    let old = motion(&reference, &refs, &[&curve], &[10.0], 0.05, 0.01)
        .with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![0.03])], || Ok(()))
        .unwrap();
    let dt = 0.2;
    let (qm, zm) = modes(0.05, 0.01, 0.03, 4.0, 4.5, 0.25, [1.0, 2.0, 3.0]);
    let vm = [0, 1, 2].map(|i| -((i + 1) as f64) * qm[i]);
    let rates = [1.0, 2.0, 3.0];
    let polar = 1.0 / core::f64::consts::TAU;
    let heat =
        polar * (9.0 * square_integral(&vm, &rates, dt) + 0.75 * square_integral(&zm, &rates, dt));
    let z = sample(&zm, &rates, dt);
    let temperature = 300.0 + (10.0 * dt + heat) / 20.0;
    let (next, audit) = old.advance_to(dt, step(), || Ok(())).unwrap();
    let state = next.torsional_joint_shafts().next().unwrap();
    let memory = state.relaxation.as_ref().unwrap();
    near(state.phases[0].temperature_k(), temperature, 3e-9);
    near(memory.spring_twists_rad[0], z, 3e-9);
    near(
        memory.branches[0].relaxation_time_s,
        2.0 / 3.0 + 0.02 * (temperature - 300.0),
        3e-10,
    );
    near(
        memory.branches[0].relaxing_stiffness_n_m2,
        polar * (0.5 + 0.2 * (temperature - 300.0)),
        3e-10,
    );
    near(
        audit.parameter_work_j,
        0.5 * polar * 0.2 * (temperature - 300.0) * z * z,
        3e-10,
    );
    near(audit.viscous_heat_j, heat, 3e-9);
    near(audit.combined_residual_j, 0.0, 3e-9);
    let (continued, continued_audit) = next.advance_to(2.0 * dt, step(), || Ok(())).unwrap();
    let restart = motion_at_enthalpy(
        &reference,
        &refs,
        &[&curve],
        &[10.0],
        next.joint_position()[0],
        next.joint_velocity()[0],
        next.thermal().cells()[&CellId(1)].specific_h_j_kg,
    )
    .with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![z])], || Ok(()))
    .unwrap();
    let (replayed, _) = restart.advance_to(dt, step(), || Ok(())).unwrap();
    near(
        continued.joint_position()[0],
        replayed.joint_position()[0],
        3e-9,
    );
    near(
        continued.thermal().cells()[&CellId(1)].specific_h_j_kg,
        replayed.thermal().cells()[&CellId(1)].specific_h_j_kg,
        3e-9,
    );
    near(continued_audit.combined_residual_j, 0.0, 3e-9);
    let erased = motion_at_enthalpy(
        &reference,
        &refs,
        &[&curve],
        &[10.0],
        next.joint_position()[0],
        next.joint_velocity()[0],
        next.thermal().cells()[&CellId(1)].specific_h_j_kg,
    )
    .with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![0.0])], || Ok(()))
    .unwrap();
    let (erased, _) = erased.advance_to(dt, step(), || Ok(())).unwrap();
    assert!((continued.joint_velocity()[0] - erased.joint_velocity()[0]).abs() > 1e-5);
}

#[test]
fn explicitly_empty_spectrum_keeps_the_original_kelvin_shaft_path_exact() {
    let card = card(8.0, 0.4, [0.0; 2], 500.0, 10.0, false);
    let curve = caloric(&card);
    let refs = [&card];
    let reference = mechanism(&refs);
    let base = motion(&reference, &refs, &[&curve], &[0.0], 0.05, 0.01);
    let original = base
        .with_torsional_joint_shafts(vec![shaft(&refs)], || Ok(()))
        .unwrap();
    let memory = ShaftRelaxation {
        properties: vec![LayerRelaxationProperties {
            equilibrium_modulus: G.into(),
            branches: vec![],
        }],
        initial_spring_twists_rad: vec![],
        absolute_spring_twist_tolerance_rad: 1e-11,
    };
    let extended = base
        .with_relaxing_torsional_joint_shafts(vec![(shaft(&refs), memory)], || Ok(()))
        .unwrap();
    let (a, aa) = original.advance_to(0.4, step(), || Ok(())).unwrap();
    let (b, ab) = extended.advance_to(0.4, step(), || Ok(())).unwrap();
    assert_eq!(a.joint_position(), b.joint_position());
    assert_eq!(a.joint_velocity(), b.joint_velocity());
    assert_eq!(a.thermal().energy_j(), b.thermal().energy_j());
    assert_eq!(aa, ab);
}

#[test]
fn source_history_domain_and_late_cancellation_refusals_keep_the_accepted_state() {
    let card = memory_card(
        card(8.0, 9.0, [0.0; 2], 500.0, 10.0, false),
        0.5,
        2.0 / 3.0,
        [0.0; 2],
        500.0,
    );
    let curve = caloric(&card);
    let refs = [&card];
    let reference = mechanism(&refs);
    let base = motion(&reference, &refs, &[&curve], &[0.0], 0.05, 0.01);
    assert!(
        base.with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![])], || Ok(()))
            .is_err()
    );
    assert!(
        base.with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![1.0])], || Ok(()))
            .is_err()
    );
    let mut wrong = definition(&refs, vec![0.03]);
    wrong.1.properties[0].branches[0].1 = G.into();
    assert!(
        base.with_relaxing_torsional_joint_shafts(vec![wrong], || Ok(()))
            .is_err()
    );
    let old = base
        .with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![0.03])], || Ok(()))
        .unwrap();
    let mut checkpoints = 0;
    let (expected, _) = old
        .advance_to(0.4, step(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    assert!(
        old.advance_to(0.4, step(), || {
            seen += 1;
            if seen == checkpoints - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .is_err()
    );
    assert_eq!(old.thermal().time_s(), 0.0);
    assert_eq!(old.joint_position(), &[0.05]);
    assert_eq!(
        old.torsional_joint_shafts()
            .next()
            .unwrap()
            .relaxation
            .as_ref()
            .unwrap()
            .spring_twists_rad,
        vec![0.03]
    );
    let (retry, _) = old.advance_to(0.4, step(), || Ok(())).unwrap();
    assert_eq!(retry.joint_position(), expected.joint_position());
    assert_eq!(
        retry.torsional_joint_shafts().next(),
        expected.torsional_joint_shafts().next()
    );
    let narrow = memory_card(
        super::card(8.0, 9.0, [0.0; 2], 500.0, 10.0, false),
        0.5,
        2.0 / 3.0,
        [0.0; 2],
        300.001,
    );
    let curve = caloric(&narrow);
    let refs = [&narrow];
    let reference = mechanism(&refs);
    let unsupported = motion(&reference, &refs, &[&curve], &[100.0], 0.05, 0.01)
        .with_relaxing_torsional_joint_shafts(vec![definition(&refs, vec![0.03])], || Ok(()))
        .unwrap();
    assert!(unsupported.advance_to(0.2, step(), || Ok(())).is_err());
    assert_eq!(unsupported.thermal().time_s(), 0.0);
}
