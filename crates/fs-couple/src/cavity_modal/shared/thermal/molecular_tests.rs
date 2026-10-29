//! Source binding and the original loaded SharedCavity runtime.
use super::*;
use fs_material::gas::GasSpec;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

fn source() -> HelmholtzModalCavity {
    HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: 1e-7,
            neck_radius_m: 1e-4,
            neck_length_m: 0.1,
            inner_end: fs_scenario::NeckEnd::Flanged,
            outer_end: fs_scenario::NeckEnd::Flanged,
        },
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
        CavityNeckLoss::Resistance(0.0),
    )
    .unwrap()
}

fn molecular(gas: &GasState) -> RelaxationState {
    // Deliberately synthetic slow rotation: an implementation oracle.
    let species = SpeciesId::new("A").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 101325.0 * 0.025,
                t_ref_k: 293.15,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 400.0,
            },
        }],
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas.clone(), vec![(species, 1.0)]).unwrap())
    .unwrap()
}

fn options() -> CylindricalCavityThermalOptions {
    CylindricalCavityThermalOptions {
        radius_m: 0.001,
        poles: 4,
        band_rad_s: [10.0, 1000.0],
        maximum_relative_error: 1e-3,
    }
}

fn coupled() -> SharedCavity {
    let source = source();
    let law = molecular(source.gas());
    let rate = 48000;
    let surface = CavitySurface {
        model: ModalAcousticTimeModel::try_new(
            rate,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 500.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap(),
        outward_volume_shapes: vec![1e-4],
    };
    let config = ModalCouplingConfig {
        max_modes: 8,
        max_connections: 6,
        max_setup_terms: 20000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-14,
        energy_relative_tolerance: 1e-8,
    };
    SharedCavity::new_with_cylindrical_molecular_wall(
        source,
        vec![surface],
        rate,
        ModalAcousticState {
            displacement_m_sqrt_kg: 2e-7,
            velocity_m_sqrt_kg_per_s: 0.0,
        },
        ModalAcousticTimeBudget::audible_reference(),
        config,
        None,
        options(),
        law,
        8,
        None,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn complex_close(a: C64, b: C64) {
    assert!((a - b).abs() < 2e-9 * b.abs().max(1e-25), "{a:?} != {b:?}");
}

#[test]
fn joint_compression_drives_loaded_motion_and_accounts_for_heat_once() {
    let gate = CancelGate::new_clock_free();
    let mut cavity = coupled();
    let source = cavity.source();
    let volume = source.geometry().volume_m3;
    let shapes = [1e-4, source.neck_volume_shape_m2_per_sqrt_kg()];
    let thermal = cavity.cavity_thermal_compression().unwrap();
    let law = thermal.molecular_relaxation().unwrap();
    let gas = law.gas();
    assert!(thermal.instantaneous_stiffness_n_m() > thermal.adiabatic_stiffness_n_m());
    assert_eq!(
        cavity.network().material_memory().len(),
        thermal.arms().len()
    );
    for omega in [50.0, 300.0, 500.0, 800.0] {
        // Direct constitutive sum and explicit two-coordinate Schur solve,
        // independent of the Maxwell eigenspectrum and coupled runtime.
        let cp = law.caloric_response(omega).unwrap().cp_j_kg_k;
        let heat = thermal
            .diffusion()
            .shunt_terms()
            .iter()
            .fold(C64::ZERO, |sum, term| {
                let weight = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
                sum + (cp + C64::new(0.0, gas.specific_heat_cp * term.rate_per_s / omega))
                    .recip()
                    .scale(weight)
            });
        let bulk = (C64::ONE - heat.scale(gas.specific_gas_constant))
            .recip()
            .scale(gas.pressure);
        let gain = bulk.scale(volume.recip());
        let a = C64::from_re(500.0_f64.powi(2) - omega * omega) + gain.scale(shapes[0].powi(2));
        let b = gain.scale(shapes[0] * shapes[1]);
        let d = C64::from_re(-omega * omega) + gain.scale(shapes[1].powi(2));
        let force = C64::new(1e-4, 2e-5);
        let determinant = a * d - b * b;
        let expected = [force * d / determinant, -force * b / determinant];
        let response = cavity
            .harmonic_response(omega, &[force], C64::ZERO, 1_000_000, &gate)
            .unwrap();
        for i in 0..2 {
            complex_close(response.network.displacement_m_sqrt_kg[i], expected[i]);
        }
        let delta_v = expected[0].scale(shapes[0]) + expected[1].scale(shapes[1]);
        complex_close(
            response.pressure_increment_pa,
            -bulk * delta_v.scale(volume.recip()),
        );
        let heat_w = -0.5 * omega * bulk.im / volume * delta_v.norm_sq();
        assert!(heat_w > 0.0);
        assert!((response.network.material_dissipation_w - heat_w).abs() < 2e-9 * heat_w);
        assert!(response.network.power_residual_w.abs() < 2e-9 * heat_w);
    }
    let initial = cavity.network().total_energy_j().unwrap();
    let (mut work, mut loss) = (0.0, 0.0);
    for sample in 0..128 {
        let external = if sample < 12 { 0.03 } else { 0.0 };
        let frame = cavity.step(&[0.0], external, &gate).unwrap();
        work += frame.external_work_j;
        loss += frame.total_dissipation_j();
        assert!(frame.material_dissipation_j >= 0.0);
    }
    let final_energy = cavity.network().total_energy_j().unwrap();
    assert!(loss > 0.0);
    assert!((final_energy + loss - initial - work).abs() <= 1e-19 + 2e-8 * (initial + work.abs()));
}

#[test]
fn joint_builder_rejects_stale_source_and_insufficient_or_cancelled_preparation() {
    let source = source();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 303.15, 101325.0).unwrap();
    let gate = CancelGate::new_clock_free();
    assert!(
        CylindricalThermalCompression::new_with_molecular_relaxation(
            &source,
            options(),
            molecular(&gas),
            8,
            &gate
        )
        .is_err()
    );
    assert!(
        CylindricalThermalCompression::new_with_molecular_relaxation(
            &source,
            options(),
            molecular(source.gas()),
            7,
            &gate
        )
        .is_err()
    );
    gate.request();
    assert!(matches!(
        CylindricalThermalCompression::new_with_molecular_relaxation(
            &source,
            options(),
            molecular(source.gas()),
            8,
            &gate
        ),
        Err(ModalCouplingError::Cancelled)
    ));
}
