#![cfg(feature = "enthalpy-coupling")]
//! Spatial latent heat, original Cp(T), and staged reservoir energy balances.

use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::ThermalLaw;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ConductivityTensor, EnergyConfig, FaceContactResistance, RotatedSolidConductivity,
    SolidMaterial, TENSOR_SDIRK2_GAMMA, TensorCaloricEndpoint, TensorCaloricLaw,
    TensorCaloricStorage, TensorNonlinearConfig, ThermalFace, ThermalSetup, Voxel, VoxelDomain,
    advance_caloric_tensor_conduction, advance_caloric_tensor_conduction_sdirk2,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const ROTATION: [[f64; 3]; 3] = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];

fn controls() -> TensorNonlinearConfig {
    TensorNonlinearConfig {
        linear: EnergyConfig {
            tolerance: 1e-10,
            ..EnergyConfig::default()
        },
        ..TensorNonlinearConfig::default()
    }
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-6, "{actual} != {expected}");
}

fn melting(density: f64) -> EquilibriumEnthalpyPhaseCurve {
    // Cp=2 J/(kg K) on both sides of a 100 J/kg isothermal latent plateau.
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xA1; 32]),
        [
            (-100.0, 250.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (300.0, 400.0, 1.0),
        ]
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: density,
        })
        .to_vec(),
    )
    .unwrap()
}

fn constant_cp(density: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        ContentHash([0xA2; 32]),
        SolidLiquidPhase::Solid,
        -100.0,
        &[250.0, 450.0].map(|t| HeatCapacityKnot {
            temperature_k: t,
            specific_heat_capacity_j_kg_k: 2.0,
            bulk_density_kg_m3: density,
        }),
    )
    .unwrap()
}

fn variable_cp(datum: f64, density: f64) -> EquilibriumEnthalpyPhaseCurve {
    // h(T)-datum=2*(T-300)+0.02*(T-300)^2, exactly, not a chord in T(h).
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        ContentHash([0xA3; 32]),
        SolidLiquidPhase::Solid,
        datum - 32.0,
        &[280.0, 400.0].map(|t| HeatCapacityKnot {
            temperature_k: t,
            specific_heat_capacity_j_kg_k: 2.0 + 0.04 * (t - 300.0),
            bulk_density_kg_m3: density,
        }),
    )
    .unwrap()
}

fn contact_setup(resistance: f64) -> ThermalSetup {
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: resistance,
    });
    setup
}

/// Two unit cells have half-cell resistances .5+.5 and contact resistance 2.
/// The receiver starts at h=20, T=300, mass=2; the donor at h=120, T=360, mass=3.
/// "carry" is the actual prior-stage energy delivered to the receiver, in J.
fn expected_contact_energy(a: f64, carry: f64) -> [f64; 2] {
    let hot_excess = (60.0 - carry / 6.0) / (1.0 + a / 18.0);
    let cold_h = 20.0 + (carry + a * hot_excess / 3.0) / 2.0;
    if cold_h <= 100.0 {
        [cold_h, 2.0 * hot_excess]
    } else {
        // After melting: 2*Tc_excess+3*Th_excess=100 from total energy.
        let cold_excess = (a * 100.0 / 9.0 + carry - 160.0) / (4.0 + 5.0 * a / 9.0);
        let hot_excess = (100.0 - 2.0 * cold_excess) / 3.0;
        [100.0 + 2.0 * cold_excess, 2.0 * hot_excess]
    }
}

fn check_contact_endpoint(
    endpoint: &TensorCaloricEndpoint,
    expected_h: [f64; 2],
    receiver: &EquilibriumEnthalpyPhaseCurve,
    donor: &EquilibriumEnthalpyPhaseCurve,
) {
    let heat = &endpoint.solution.solution;
    let h = &endpoint.specific_energy_j_kg;
    for cell in 0..2 {
        close(h[cell], expected_h[cell]);
    }
    let cold = receiver.values_at_specific_enthalpy(h[0]).unwrap();
    let hot = donor.values_at_specific_enthalpy(h[1]).unwrap();
    close(heat.temperature[0], cold.temperature_k);
    close(heat.temperature[1], hot.temperature_k);
    if h[0] < 100.0 {
        assert_eq!(heat.temperature[0], 300.0);
        close(cold.liquid_mass_fraction, h[0] / 100.0);
        close(receiver.temperature_slope_kg_k_j(h[0]).unwrap(), 0.0);
    } else {
        assert!(heat.temperature[0] > 300.0);
        assert_eq!(cold.liquid_mass_fraction, 1.0);
    }

    let q = (hot.temperature_k - cold.temperature_k) / 3.0;
    close(heat.outward_heat_w[0][1], -q);
    close(heat.outward_heat_w[1][0], q);
    close(heat.face_temperature_k[0][1], cold.temperature_k + 0.5 * q);
    close(heat.face_temperature_k[1][0], hot.temperature_k - 0.5 * q);
    close(
        heat.face_temperature_k[1][0] - heat.face_temperature_k[0][1],
        2.0 * q,
    );
    for cell in 0..2 {
        for face in 0..6 {
            if (cell, face) != (0, 1) && (cell, face) != (1, 0) {
                close(heat.outward_heat_w[cell][face], 0.0);
            }
        }
    }
    close(2.0 * (h[0] - 20.0) + 3.0 * (h[1] - 120.0), 0.0);

    // Independent phase-aware entropy: latent heat/Tm plus sensible Cp*ln(T).
    let cold_entropy = (h[0].min(100.0) - 20.0) / 300.0
        + if h[0] > 100.0 {
            2.0 * (cold.temperature_k / 300.0).ln()
        } else {
            0.0
        };
    let entropy = 2.0 * cold_entropy + 6.0 * (hot.temperature_k / 360.0).ln();
    let original_entropy = 2.0
        * ThermalLaw::specific_entropy_change_j_kg_k(receiver, 20.0, h[0])
            .unwrap()
            .unwrap()
        + 3.0
            * ThermalLaw::specific_entropy_change_j_kg_k(donor, 120.0, h[1])
                .unwrap()
                .unwrap();
    close(original_entropy, entropy);
    assert!(entropy > 0.0);
    assert!(q * (1.0 / cold.temperature_k - 1.0 / hot.temperature_k) > 0.0);
    assert!(q * (1.0 / heat.face_temperature_k[0][1] - 1.0 / heat.face_temperature_k[1][0]) > 0.0);
}

#[test]
fn contact_transfers_energy_through_a_latent_plateau_and_into_liquid() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("fixed unit conductivity", 1.0)];
    let conductivity = ConductivityTensor::new(IDENTITY).unwrap();
    let receiver = melting(2.0);
    let donor = constant_cp(3.0);
    let caloric: [&dyn TensorCaloricLaw; 2] = [&receiver, &donor];
    let mass = [2.0, 3.0];
    let previous = [20.0, 120.0];
    let setup = contact_setup(2.0);
    let saved = setup.clone();

    for dt in [3.0, 18.0] {
        for sdirk in [false, true] {
            let storage = TensorCaloricStorage {
                laws: &caloric,
                reference_mass_kg: &mass,
                previous_specific_energy_j_kg: &previous,
            };
            let (endpoint, record, expected) = if sdirk {
                let a = TENSOR_SDIRK2_GAMMA * dt;
                let b = (1.0 - TENSOR_SDIRK2_GAMMA) * dt;
                let step = advance_caloric_tensor_conduction_sdirk2(
                    &domain,
                    &solids,
                    [&setup, &setup],
                    &[&conductivity, &conductivity],
                    storage,
                    [&[], &[]],
                    4.0,
                    dt,
                    &controls(),
                    &CancelGate::new(),
                )
                .unwrap();
                let [first, last] = step.stages;
                let first_expected = expected_contact_energy(a, 0.0);
                check_contact_endpoint(&first, first_expected, &receiver, &donor);
                let q1 = first_expected[1] / 6.0;
                close(step.stage_time_s[0], 4.0 + a);
                close(step.stage_time_s[1], 4.0 + dt);
                (last, step.record, expected_contact_energy(a, b * q1))
            } else {
                let step = advance_caloric_tensor_conduction(
                    &domain,
                    &solids,
                    &setup,
                    &[&conductivity, &conductivity],
                    storage,
                    &[],
                    4.0,
                    dt,
                    &controls(),
                    &CancelGate::new(),
                )
                .unwrap();
                (step.endpoint, step.record, expected_contact_energy(dt, 0.0))
            };
            check_contact_endpoint(&endpoint, expected, &receiver, &donor);
            assert_eq!(endpoint.specific_energy_j_kg[0] < 100.0, dt == 3.0);
            close(record.stored_energy_change_j, 0.0);
            close(record.source_j, 0.0);
            close(record.boundary_outflow_j, 0.0);
            close(record.closure_j, 0.0);
            close(record.time_s, 4.0 + dt);
        }
    }
    assert_eq!(previous, [20.0, 120.0]);
    assert_eq!(setup, saved);
}

#[test]
fn variable_cp_and_rotated_original_k_integrals_close_every_face_and_energy_datum() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("rotated K(T)=T/300 * K0", 1.0)
        .with_orthotropic([4.0, 1.0, 2.0])
        .with_conductivity_table(&[(250.0, 250.0 / 300.0), (500.0, 500.0 / 300.0)])];
    let conductivity = RotatedSolidConductivity::new(&solids[0], ROTATION).unwrap();
    let target: f64 = 340.0;
    let gradient = [0.0, 10.0, -4.0];
    // In the Kirchhoff coordinate phi=T^2/600, traces have affine gradient
    // plus an isotropic -5 offset. Native HMM gives 140 W total outward heat.
    let faces = std::array::from_fn(|face| {
        let sign = if face & 1 == 0 { -1.0 } else { 1.0 };
        ThermalFace::Temperature(
            (target * target + 600.0 * (-5.0 + sign * gradient[face / 2] / 2.0)).sqrt(),
        )
    });
    let mut setup = ThermalSetup::new(faces);
    setup.power_w = vec![332.0]; // 1.5*(h(340)-h(320))/.5 + 140.
    let expected_face_heat = [43.6, 14.8, 41.6, 0.0, 12.0, 28.0];

    for datum in [0.0, 1e9] {
        let curve = variable_cp(datum, 1.5);
        let caloric: [&dyn TensorCaloricLaw; 1] = [&curve];
        let previous = [datum + 48.0];
        let step = advance_caloric_tensor_conduction(
            &domain,
            &solids,
            &setup,
            &[&conductivity],
            TensorCaloricStorage {
                laws: &caloric,
                reference_mass_kg: &[1.5],
                previous_specific_energy_j_kg: &previous,
            },
            &[],
            0.0,
            0.5,
            &controls(),
            &CancelGate::new(),
        )
        .unwrap();
        close(step.endpoint.specific_energy_j_kg[0] - datum, 112.0);
        close(step.endpoint.solution.solution.temperature[0], target);
        for face in 0..6 {
            close(
                step.endpoint.solution.solution.outward_heat_w[0][face],
                expected_face_heat[face],
            );
        }
        close(step.record.stored_energy_change_j, 96.0);
        close(step.record.source_j, 166.0);
        close(step.record.boundary_outflow_j, 70.0);
        close(step.record.closure_j, 0.0);
        assert_eq!(previous, [datum + 48.0]);
    }
}

#[test]
fn sdirk_integrates_variable_cp_in_energy_coordinates_at_second_order() {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    // k(T)=Cp(T); a unit face at distance .5 carries Q=2*h(T) to 300 K.
    // With fixed mass 2, the exact caloric ODE is dh/dt=-h.
    let solids = [SolidMaterial::new("k(T)=Cp(T)", 2.0)
        .with_conductivity_table(&[(280.0, 1.2), (400.0, 6.0)])];
    let conductivity = RotatedSolidConductivity::new(&solids[0], IDENTITY).unwrap();
    let curve = variable_cp(0.0, 2.0);
    let caloric: [&dyn TensorCaloricLaw; 1] = [&curve];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.faces[0] = ThermalFace::Temperature(300.0);
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let dt = 1.0 / f64::from(count);
        let a = TENSOR_SDIRK2_GAMMA * dt;
        let b = (1.0 - TENSOR_SDIRK2_GAMMA) * dt;
        let amplification = (1.0 - (b - a)) / (1.0 + a).powi(2);
        let mut h = [112.0];
        for index in 0..count {
            let step = advance_caloric_tensor_conduction_sdirk2(
                &domain,
                &solids,
                [&setup, &setup],
                &[&conductivity],
                TensorCaloricStorage {
                    laws: &caloric,
                    reference_mass_kg: &[2.0],
                    previous_specific_energy_j_kg: &h,
                },
                [&[], &[]],
                f64::from(index) * dt,
                dt,
                &controls(),
                &CancelGate::new(),
            )
            .unwrap();
            close(step.stages[0].specific_energy_j_kg[0], h[0] / (1.0 + a));
            close(step.stages[1].specific_energy_j_kg[0], h[0] * amplification);
            close(
                step.record.stored_energy_change_j,
                2.0 * (step.stages[1].specific_energy_j_kg[0] - h[0]),
            );
            close(step.record.source_j, 0.0);
            close(
                step.record.boundary_outflow_j,
                -step.record.stored_energy_change_j,
            );
            close(step.record.closure_j, 0.0);
            h[0] = step.stages[1].specific_energy_j_kg[0];
        }
        close(h[0], 112.0 * amplification.powi(count));
        errors.push((h[0] - 112.0 * (-1.0_f64).exp()).abs());
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn changing_prescribed_masks_keep_explicit_latent_energy_and_reservoir_carry() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("mask-switching finite contact", 1.0)];
    let conductivity = ConductivityTensor::new(IDENTITY).unwrap();
    let receiver = melting(2.0);
    let donor = melting(3.0);
    let caloric: [&dyn TensorCaloricLaw; 2] = [&receiver, &donor];
    let previous = [20.0, 120.0]; // T=[300,310]; both use their own mass.
    let mut first = contact_setup(0.25);
    first.fixed_temperature = vec![(0, 300.0)];
    first.power_w = vec![1000.0, 2.0]; // Fixed-cell source is owned by its thermostat.
    let mut second = first.clone();
    second.fixed_temperature = vec![(1, 330.0)];
    second.power_w = vec![3.0, 2000.0];
    let saved = [first.clone(), second.clone()];
    let a = TENSOR_SDIRK2_GAMMA;
    let b = 1.0 - TENSOR_SDIRK2_GAMMA;
    let step = advance_caloric_tensor_conduction_sdirk2(
        &domain,
        &solids,
        [&first, &second],
        &[&conductivity, &conductivity],
        TensorCaloricStorage {
            laws: &caloric,
            reference_mass_kg: &[2.0, 3.0],
            previous_specific_energy_j_kg: &previous,
        },
        [&[(0, 40.0)], &[(1, 160.0)]],
        0.0,
        1.0,
        &controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let h1 = &step.stages[0].specific_energy_j_kg;
    let h2 = &step.stages[1].specific_energy_j_kg;
    let heat1 = &step.stages[0].solution.solution;
    let heat2 = &step.stages[1].solution.solution;
    let q1 = 0.8 * (300.0 - (300.0 + (60.0 + 2.0 * a) / (6.0 + 0.8 * a)));
    let q2 = 0.8 * (300.0 - 330.0);
    let stage1_left_rate = 2.0 * (40.0 - previous[0]) / a;
    let expected_h2_left = previous[0] + (a * (3.0 - q2) + b * stage1_left_rate) / 2.0;
    close(h1[0], 40.0);
    close(h1[1], 120.0 + a * (2.0 + q1) / 3.0);
    close(h2[0], expected_h2_left);
    close(h2[1], 160.0);
    assert!(h2[0] > h1[0] && h2[0] < 100.0);
    assert_eq!(heat1.temperature[0], 300.0);
    assert_eq!(heat2.temperature[0], 300.0);
    close(heat2.temperature[1], 330.0);
    close(heat1.outward_heat_w[0][1], q1);
    close(heat2.outward_heat_w[0][1], q2);
    let reservoir1 = -q1 - stage1_left_rate;
    let reservoir2 = q2 - 3.0 * (160.0 - previous[1]) / a + b / a * (2.0 + q1);
    close(step.record.source_j, b * 2.0 + a * 3.0);
    close(
        step.record.boundary_outflow_j,
        b * reservoir1 + a * reservoir2,
    );
    close(
        step.record.stored_energy_change_j,
        2.0 * (h2[0] - previous[0]) + 3.0 * (h2[1] - previous[1]),
    );
    close(step.record.closure_j, 0.0);
    assert_eq!(previous, [20.0, 120.0]);
    assert_eq!([first, second], saved);
}

#[test]
fn ambiguous_or_inconsistent_prescribed_energy_and_cancellation_publish_no_state() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("explicit latent thermostat", 1.0)];
    let conductivity = ConductivityTensor::new(IDENTITY).unwrap();
    let receiver = melting(2.0);
    let donor = melting(3.0);
    let caloric: [&dyn TensorCaloricLaw; 2] = [&receiver, &donor];
    let previous = [20.0, 120.0];
    let mut setup = contact_setup(0.25);
    setup.fixed_temperature = vec![(0, 300.0)];
    let saved = setup.clone();
    let solve = |prescribed: &[(usize, f64)], gate: &CancelGate| {
        advance_caloric_tensor_conduction(
            &domain,
            &solids,
            &setup,
            &[&conductivity, &conductivity],
            TensorCaloricStorage {
                laws: &caloric,
                reference_mass_kg: &[2.0, 3.0],
                previous_specific_energy_j_kg: &previous,
            },
            prescribed,
            0.0,
            1.0,
            &controls(),
            gate,
        )
    };
    assert!(solve(&[(0, 40.0)], &CancelGate::new()).is_ok());
    for prescribed in [
        vec![],
        vec![(0, 160.0)],
        vec![(0, 40.0), (1, 120.0)],
        vec![(0, 40.0), (0, 40.0)],
    ] {
        assert!(matches!(
            solve(&prescribed, &CancelGate::new()),
            Err(ChtError::InvalidInput { .. })
        ));
    }
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(solve(&[(0, 40.0)], &gate).unwrap_err(), ChtError::Cancelled);

    let bad_second = [(0, 301.0)]; // Outside the original h span, even though T=300 is supported.
    assert!(matches!(
        advance_caloric_tensor_conduction_sdirk2(
            &domain,
            &solids,
            [&setup, &setup],
            &[&conductivity, &conductivity],
            TensorCaloricStorage {
                laws: &caloric,
                reference_mass_kg: &[2.0, 3.0],
                previous_specific_energy_j_kg: &previous,
            },
            [&[(0, 40.0)], &bad_second],
            0.0,
            1.0,
            &controls(),
            &CancelGate::new(),
        ),
        Err(ChtError::InvalidInput { .. })
    ));
    assert_eq!(previous, [20.0, 120.0]);
    assert_eq!(setup, saved);
}
