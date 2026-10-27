//! G1/G3 real shared-step sensor/actuator and circuit tests on authored material data.
use fs_couple::piezoelectric::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::{PortHamiltonian, QuadraticStorage, interconnect, step};
use fs_qty::Density;

fn rod(d: f64, density: f64) -> PiezoelectricRod {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, 0.1),
        (STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS, d),
        (PERMITTIVITY_33_T_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, 2.0),
        ("density", Density::DIMS, density),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic reciprocal sensor/actuator fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored piezoelectric solid".into(),
            phase: "solid".into(),
            process: "synthetic longitudinal reduction".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let law = LongitudinalPiezoelectric::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PiezoelectricRod::new(
        law,
        2.0,
        0.5,
        SmallSignalLimits {
            max_abs_strain: 0.05,
            max_abs_electric_field_v_m: 10.0,
        },
    )
    .unwrap()
}
fn controls() -> PiezoelectricSolve {
    PiezoelectricSolve {
        energy_tolerance_j: 1e-11,
        voltage_tolerance_v: 1e-11,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance {tolerance}"
    );
}

#[test]
fn g1_material_data_alone_change_the_solved_sensor_motion_voltage_and_inertia() {
    for (d, density) in [(0.2, 3.0), (-0.2, 3.0), (0.4, 3.0), (0.2, 6.0), (0.0, 3.0)] {
        let rod = rod(d, density);
        let mass = density / 3.0; // rho A L/3 with A=0.5, L=2
        let epsilon = 2.0 - d * d / 0.1;
        let stiffness = 0.25 * (10.0 + (d / 0.1) * (d / 0.1) / epsilon);
        near(rod.modal_mass_kg(), mass, 1e-14);
        let body = PiezoelectricBody::new(rod, 0.0, 0.0, 0.0).unwrap();
        let before = body.clone();
        let (force, dt) = (0.03, 0.2);
        let candidate = body
            .step_candidate(force, ElectricalDrive::Current(0.0), dt, controls(), || {
                false
            })
            .unwrap();
        let next = candidate.body();
        // Independent linear midpoint oscillator, fixed charge: k_D replaces k_E.
        let extension = force * dt * dt / (2.0 * mass) / (1.0 + stiffness * dt * dt / (4.0 * mass));
        near(next.extension_m(), extension, 2e-13);
        near(next.velocity_m_s(), 2.0 * extension / dt, 2e-12);
        near(
            next.voltage_v().unwrap(),
            -(d / 0.1) * extension / epsilon,
            2e-12,
        );
        near(next.charge_c(), 0.0, 1e-14);
        near(next.stored_energy_j().unwrap(), force * extension, 1e-12);
        assert_eq!(candidate.step().electrical_supply_j, 0.0);
        assert_eq!(candidate.step().native.dissipated, 0.0);
        assert_eq!(body, before); // proposing never advances accepted history
        assert_eq!(next.rod().material(), body.rod().material());
    }
}

#[test]
fn g1_short_and_open_circuits_have_their_distinct_actual_mechanical_stiffness() {
    let rod = rod(0.2, 3.0);
    let (extension, dt) = (0.01, 0.2);
    for short in [false, true] {
        let charge = if short {
            rod.coupling_c_m() * extension
        } else {
            0.0
        };
        let body = PiezoelectricBody::new(rod.clone(), extension, 0.0, charge).unwrap();
        let drive = if short {
            ElectricalDrive::Voltage(0.0)
        } else {
            ElectricalDrive::Current(0.0)
        };
        let candidate = body
            .step_candidate(0.0, drive, dt, controls(), || false)
            .unwrap();
        let k = if short { 2.5 } else { 3.125 };
        let a = k * dt * dt / 4.0;
        let expected = extension * (1.0 - a) / (1.0 + a);
        near(candidate.body().extension_m(), expected, 2e-12);
        near(
            candidate.body().stored_energy_j().unwrap(),
            body.stored_energy_j().unwrap(),
            1e-12,
        );
        assert_eq!(candidate.step().electrical_supply_j, 0.0);
        if short {
            near(candidate.body().voltage_v().unwrap(), 0.0, 1e-12);
            near(
                candidate.body().charge_c(),
                rod.coupling_c_m() * expected,
                2e-12,
            );
        } else {
            near(candidate.body().charge_c(), 0.0, 1e-14);
        }
    }
}

#[test]
fn g1_voltage_and_current_drives_exchange_work_with_real_material_motion() {
    for d in [-0.2, 0.0, 0.2] {
        let rod = rod(d, 3.0);
        let (voltage, dt) = (0.2, 0.1);
        let charge = rod.axially_clamped_capacitance_f() * voltage;
        let body = PiezoelectricBody::new(rod.clone(), 0.0, 0.0, charge).unwrap();
        let candidate = body
            .step_candidate(
                0.0,
                ElectricalDrive::Voltage(voltage),
                dt,
                controls(),
                || false,
            )
            .unwrap();
        let expected = rod.coupling_c_m() * voltage * dt * dt / 2.0 / (1.0 + 2.5 * dt * dt / 4.0);
        near(candidate.body().extension_m(), expected, 2e-12);
        near(
            candidate.body().charge_c(),
            charge + rod.coupling_c_m() * expected,
            2e-12,
        );
        near(candidate.body().voltage_v().unwrap(), voltage, 1e-12);
        near(
            candidate.step().electrical_supply_j,
            voltage * rod.coupling_c_m() * expected,
            1e-12,
        );
        near(candidate.step().energy_residual_j, 0.0, 1e-12);
        // The direct/converse sign follows the same signed d in the material card.
        if d != 0.0 {
            assert_eq!(candidate.body().extension_m().signum(), d.signum());
        }
    }
    let rod = rod(0.2, 3.0);
    let body = PiezoelectricBody::new(rod.clone(), 0.0, 0.0, 0.0).unwrap();
    let (current, dt) = (0.05, 0.2);
    let candidate = body
        .step_candidate(
            0.0,
            ElectricalDrive::Current(current),
            dt,
            controls(),
            || false,
        )
        .unwrap();
    let midpoint_charge = current * dt / 2.0;
    let expected =
        (rod.coupling_c_m() / rod.axially_clamped_capacitance_f()) * midpoint_charge * dt * dt
            / 2.0
            / (1.0 + rod.open_circuit_stiffness_n_m() * dt * dt / 4.0);
    near(candidate.body().charge_c(), current * dt, 1e-13);
    near(candidate.body().extension_m(), expected, 2e-12);
    near(candidate.step().energy_residual_j, 0.0, 1e-12);
}

#[test]
fn g3_the_same_material_component_composes_with_native_rl_circuit_ports() {
    let rod = rod(0.2, 3.0);
    // A series R-L branch has flux as its native storage, voltage input and
    // current output. Interconnection closes Qdot=-I and phidot=V-RI.
    let branch = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.1],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![2.0], 1).unwrap()),
    )
    .unwrap();
    let system = interconnect(rod.port_hamiltonian().unwrap(), branch, &[(1, 0)]).unwrap();
    let mut state = vec![0.0, 0.0, 0.02, 0.0];
    let initial = system.hamiltonian(&state);
    let mut dissipated = 0.0;
    for _ in 0..20 {
        let record = step(&system, &state, &[0.0], 0.02).unwrap();
        dissipated += record.dissipated;
        near(record.balance_residual(), 0.0, 1e-12);
        state = record.x;
    }
    assert!(state[0].abs() > 1e-5); // actual converse mechanical response
    assert!(state[3].abs() > 1e-3); // actual electrical current/flux exchange
    assert!(dissipated > 0.0);
    near(system.hamiltonian(&state) + dissipated, initial, 1e-11);
}

#[test]
fn g4_incompatible_clamp_limits_and_cancellation_leave_the_whole_state_unchanged() {
    let body = PiezoelectricBody::new(rod(0.2, 3.0), 0.0, 0.0, 0.0).unwrap();
    let before = body.clone();
    assert!(
        body.step_candidate(0.0, ElectricalDrive::Voltage(0.2), 0.1, controls(), || {
            false
        })
        .is_err()
    );
    assert!(
        body.step_candidate(0.03, ElectricalDrive::Current(0.0), 0.1, controls(), || {
            true
        })
        .is_err()
    );
    assert!(
        body.step_candidate(1e9, ElectricalDrive::Current(0.0), 1.0, controls(), || {
            false
        })
        .is_err()
    );
    assert!(
        body.step_candidate(
            0.0,
            ElectricalDrive::Current(f64::NAN),
            0.1,
            controls(),
            || false
        )
        .is_err()
    );
    let mut calls = 0;
    assert!(
        body.step_candidate(0.03, ElectricalDrive::Current(0.0), 0.1, controls(), || {
            calls += 1;
            calls == 3
        })
        .is_err()
    );
    assert_eq!(body, before);
}
