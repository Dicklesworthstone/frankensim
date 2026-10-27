//! Native coil/core motion from synthetic reciprocal tangent data.
use fs_couple::magnetostrictive::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::MAGNETIC_PERMEABILITY_DIMS;
use fs_material::magnetostrictive::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::{PortHamiltonian, QuadraticStorage, interconnect, step};
use fs_qty::{Density, QuantitySpec};

fn rod(d: f64, density: f64, turns: u32) -> MagnetostrictiveRod {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_H_PROPERTY, COMPLIANCE_DIMS, 0.1),
        (PIEZOMAGNETIC_33_PROPERTY, PIEZOMAGNETIC_DIMS, d),
        (PERMEABILITY_33_T_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, 2.0),
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
                    source: "synthetic biased coil/core fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic magnetostrictive core".into(),
            phase: "solid".into(),
            process: "reversible axis-3 tangent".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
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
fn controls() -> MagnetostrictiveSolve {
    MagnetostrictiveSolve {
        energy_tolerance_j: 1e-11,
        current_tolerance_a: 1e-11,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance {tolerance}"
    );
}

#[test]
fn g1_material_data_and_turns_change_computed_motion_and_sensing() {
    for (d, density, turns) in [
        (0.2, 3.0, 3),
        (-0.2, 3.0, 3),
        (0.4, 3.0, 3),
        (0.2, 6.0, 3),
        (0.2, 3.0, 6),
        (0.0, 3.0, 3),
    ] {
        let rod = rod(d, density, turns);
        let mass = density / 3.0;
        let determinant = 0.1 * 2.0 - d * d;
        let stiffness = 0.25 * 2.0 / determinant;
        let body = MagnetostrictiveBody::new(rod.clone(), 0.0, 0.0, 0.0).unwrap();
        let (force, dt) = (0.03, 0.2);
        let candidate = body
            .step_candidate(
                force,
                CoilDrive::VoltageIncrement(0.0),
                dt,
                controls(),
                || false,
            )
            .unwrap();
        let extension = force * dt * dt / (2.0 * mass) / (1.0 + stiffness * dt * dt / (4.0 * mass));
        near(candidate.body().extension_m(), extension, 2e-12);
        near(candidate.body().velocity_m_s(), 2.0 * extension / dt, 2e-11);
        near(
            candidate.body().current_increment_a().unwrap(),
            -d * extension / (f64::from(turns) * determinant),
            2e-12,
        );
        near(candidate.body().flux_linkage_increment_wb(), 0.0, 1e-14);
        near(
            candidate.body().incremental_energy_j().unwrap(),
            force * extension,
            2e-12,
        );
        near(
            rod.clamped_inductance_h(),
            (2.0 - d * d / 0.1) * f64::from(turns).powi(2) * 0.25,
            1e-13,
        );
        assert_eq!(candidate.step().incremental_electrical_supply_j, 0.0);
        assert_eq!(body.state(), [0.0; 3]);
        assert_eq!(candidate.body().rod().material(), body.rod().material());
    }
}

#[test]
fn g1_fixed_current_and_fixed_flux_have_distinct_mechanical_frequencies() {
    let rod = rod(0.2, 3.0, 3);
    let (x0, dt) = (0.01, 0.2);
    for fixed_current in [false, true] {
        let lambda = if fixed_current {
            rod.coupling_wb_m() * x0
        } else {
            0.0
        };
        let body = MagnetostrictiveBody::new(rod.clone(), x0, 0.0, lambda).unwrap();
        let drive = if fixed_current {
            CoilDrive::CurrentIncrement(0.0)
        } else {
            CoilDrive::VoltageIncrement(0.0)
        };
        let candidate = body
            .step_candidate(0.0, drive, dt, controls(), || false)
            .unwrap();
        let k = if fixed_current { 2.5 } else { 3.125 };
        let a = k * dt * dt / 4.0;
        let expected = x0 * (1.0 - a) / (1.0 + a);
        near(candidate.body().extension_m(), expected, 2e-12);
        near(
            candidate.body().incremental_energy_j().unwrap(),
            body.incremental_energy_j().unwrap(),
            1e-12,
        );
        if fixed_current {
            near(candidate.body().current_increment_a().unwrap(), 0.0, 1e-12);
            near(
                candidate.body().flux_linkage_increment_wb(),
                rod.coupling_wb_m() * expected,
                1e-12,
            );
        } else {
            near(candidate.body().flux_linkage_increment_wb(), lambda, 1e-13);
        }
    }
}

#[test]
fn g1_held_current_and_voltage_generate_reciprocal_motion_and_work() {
    for d in [-0.2, 0.0, 0.2] {
        let rod = rod(d, 3.0, 3);
        let (current, dt) = (0.2, 0.1);
        let initial_flux = rod.clamped_inductance_h() * current;
        let body = MagnetostrictiveBody::new(rod.clone(), 0.0, 0.0, initial_flux).unwrap();
        let candidate = body
            .step_candidate(
                0.0,
                CoilDrive::CurrentIncrement(current),
                dt,
                controls(),
                || false,
            )
            .unwrap();
        let expected = rod.coupling_wb_m() * current * dt * dt / 2.0 / (1.0 + 2.5 * dt * dt / 4.0);
        near(candidate.body().extension_m(), expected, 2e-12);
        near(
            candidate.body().flux_linkage_increment_wb(),
            initial_flux + rod.coupling_wb_m() * expected,
            1e-12,
        );
        near(
            candidate.body().current_increment_a().unwrap(),
            current,
            1e-12,
        );
        near(
            candidate.step().incremental_electrical_supply_j,
            current * rod.coupling_wb_m() * expected,
            1e-12,
        );
        near(candidate.step().incremental_energy_residual_j, 0.0, 1e-12);
        if d != 0.0 {
            assert_eq!(candidate.body().extension_m().signum(), d.signum());
        }
    }
    let rod = rod(0.2, 3.0, 3);
    let body = MagnetostrictiveBody::new(rod.clone(), 0.0, 0.0, 0.0).unwrap();
    let (voltage, dt) = (0.05, 0.2);
    let candidate = body
        .step_candidate(
            0.0,
            CoilDrive::VoltageIncrement(voltage),
            dt,
            controls(),
            || false,
        )
        .unwrap();
    let midpoint_linkage = voltage * dt / 2.0;
    let expected = (rod.coupling_wb_m() / rod.clamped_inductance_h()) * midpoint_linkage * dt * dt
        / 2.0
        / (1.0 + rod.fixed_flux_stiffness_n_m() * dt * dt / 4.0);
    near(
        candidate.body().flux_linkage_increment_wb(),
        voltage * dt,
        1e-13,
    );
    near(candidate.body().extension_m(), expected, 2e-12);
    near(candidate.step().incremental_energy_residual_j, 0.0, 1e-12);
}

#[test]
fn g1_native_coil_capacitor_interconnection_recovers_both_coupled_modes() {
    let rod = rod(0.2, 3.0, 3);
    let (mass, kh, kb, lc, theta, c) = (1.0, 2.5, 3.125, 3.6, 1.5, 0.4);
    // Independent eigenvalues of the second-order x/lambda equations:
    // trace=k_B/m+1/(L_c*C), determinant=k_H/(m*L_c*C).
    let trace: f64 = kb / mass + 1.0 / (lc * c);
    let determinant = kh / (mass * lc * c);
    let split = (trace * trace - 4.0 * determinant).sqrt();
    for frequency_squared in [0.5 * (trace - split), 0.5 * (trace + split)] {
        let capacitor = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![0.0],
            vec![1.0],
            Box::new(QuadraticStorage::new(vec![1.0 / c], 1).unwrap()),
        )
        .unwrap();
        let system = interconnect(rod.port_hamiltonian().unwrap(), capacitor, &[(1, 0)]).unwrap();
        let x0 = 0.001;
        let lambda0 = (kb - mass * frequency_squared) * lc / theta * x0;
        let state = vec![x0, 0.0, lambda0, 0.0];
        let dt = 0.2;
        let candidate = step(&system, &state, &[0.0], dt).unwrap();
        let a = frequency_squared * dt * dt / 4.0;
        let factor = (1.0 - a) / (1.0 + a);
        near(candidate.x[0], factor * x0, 2e-12);
        near(
            candidate.x[1],
            -mass * frequency_squared * dt * x0 / (1.0 + a),
            2e-12,
        );
        near(candidate.x[2], factor * lambda0, 2e-12);
        near(
            system.hamiltonian(&candidate.x),
            system.hamiltonian(&state),
            2e-12,
        );
        assert!(candidate.x[3].abs() > 1e-7);
        near(candidate.balance_residual(), 0.0, 2e-12);
    }
}

#[test]
fn g4_clamp_amplitude_and_late_cancellation_preserve_joint_history() {
    let body = MagnetostrictiveBody::new(rod(0.2, 3.0, 3), 0.0, 0.0, 0.0).unwrap();
    let before = body.clone();
    assert!(
        body.step_candidate(
            0.0,
            CoilDrive::CurrentIncrement(0.2),
            0.1,
            controls(),
            || false
        )
        .is_err()
    );
    assert!(
        body.step_candidate(
            1e9,
            CoilDrive::VoltageIncrement(0.0),
            1.0,
            controls(),
            || false
        )
        .is_err()
    );
    assert!(
        body.step_candidate(
            0.0,
            CoilDrive::VoltageIncrement(f64::NAN),
            0.1,
            controls(),
            || false
        )
        .is_err()
    );
    let mut count = 0;
    let accepted = body
        .step_candidate(
            0.03,
            CoilDrive::VoltageIncrement(0.0),
            0.1,
            controls(),
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        body.step_candidate(
            0.03,
            CoilDrive::VoltageIncrement(0.0),
            0.1,
            controls(),
            || {
                seen += 1;
                seen == count
            }
        ),
        Err(MagnetostrictiveCouplingError::Cancelled)
    ));
    assert_eq!(body, before);
    let retry = body
        .step_candidate(
            0.03,
            CoilDrive::VoltageIncrement(0.0),
            0.1,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(retry.body(), accepted.body());
    let rod = body.rod().clone();
    let mut limits = rod.limits();
    limits.maximum_abs_field_increment_a_m = 10.0;
    assert!(MagnetostrictiveRod::new(rod.material().clone(), 2.0, 0.5, 3, limits).is_err());
    assert!(MagnetostrictiveRod::new(rod.material().clone(), 2.0, 0.5, 0, rod.limits()).is_err());
}
