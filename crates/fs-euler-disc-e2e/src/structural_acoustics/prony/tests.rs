use super::*;
use fs_couple::modal_acoustic_time::ModalAcousticMode;

fn budget() -> StructuralPronyAudioBudget {
    StructuralPronyAudioBudget {
        modal: ModalAcousticTimeBudget::audible_reference(),
        coupling: ModalCouplingConfig {
            max_modes: 64,
            max_connections: 0,
            max_setup_terms: 1_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1.0e6,
            maximum_abs_connection_force_n: 1.0e6,
            solve_relative_tolerance: 1.0e-11,
            energy_absolute_tolerance_j: 1.0e-12,
            energy_relative_tolerance: 1.0e-9,
        },
        prepared: PreparedModalBudget {
            max_modes: 64,
            max_dense_bytes: 536_870_912,
            max_multiply_terms: 100_000_000_000,
        },
        maximum_material_branches: 256,
    }
}

fn elastic(omegas: &[f64]) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        48_000,
        omegas
            .iter()
            .map(|&omega| ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::from_re(2.0),
            })
            .collect(),
        budget().modal,
    )
    .unwrap()
}

// Synthetic exact Prony coefficients isolate consumer projection and native
// memory dynamics. No measured fractional fit or material calibration is claimed.
fn material(terms: Vec<(f64, f64)>) -> LoweredModel {
    LoweredModel {
        model: GeneralizedMaxwell::new(2.0e9, terms).unwrap(),
        band: (10.0, 10_000.0),
        sup_rel_err: 0.0,
        sup_rel_storage_err: 0.0,
        sup_rel_loss_err: 0.0,
        verification_points: 2,
    }
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}",
    );
}

#[test]
fn g1_native_prony_response_contains_storage_dispersion_and_correct_heat() {
    let gate = CancelGate::new_clock_free();
    let law = material(vec![(1.6e9, 0.001), (0.4e9, 0.004)]);
    let omega0 = 1_000.0;
    let system = bind_prony_memory(elastic(&[omega0]), &law, budget(), &gate).unwrap();
    assert_eq!(system.material_relaxation_branches().len(), 2);
    for omega in [200.0, 700.0, 1_000.0, 1_400.0, 4_000.0] {
        let observed = system
            .harmonic_response(omega, &[C64::from_re(1.0)], 100_000, &gate)
            .unwrap();
        let (storage, loss) = law.modulus_checked(omega).unwrap();
        // The shared harmonic owner uses exp(-i omega t), hence -E''.
        let expected = C64::from_re(1.0)
            / C64::new(
                omega0 * omega0 * storage / law.model.e_inf - omega * omega,
                -omega0 * omega0 * loss / law.model.e_inf,
            );
        close(
            observed.displacement_m_sqrt_kg[0].re,
            expected.re,
            2.0e-10,
            1.0e-16,
        );
        close(
            observed.displacement_m_sqrt_kg[0].im,
            expected.im,
            2.0e-10,
            1.0e-16,
        );
        let q2 = expected.re * expected.re + expected.im * expected.im;
        let heat = 0.5 * omega * omega0 * omega0 * (loss / law.model.e_inf) * q2;
        close(observed.material_dissipation_w, heat, 2.0e-10, 1.0e-16);
        close(observed.external_power_w, heat, 2.0e-10, 1.0e-16);
        assert_eq!(observed.component_dissipation_w, 0.0);
        let mut stored = 0.25 * (omega0 * omega0 + omega * omega) * q2;
        for &(increment, tau) in &law.model.terms {
            let x = omega * tau;
            stored += 0.25
                * omega0
                * omega0
                * (increment / law.model.e_inf)
                * (x * x / (1.0 + x * x))
                * q2;
        }
        close(observed.mean_stored_energy_j, stored, 2.0e-10, 1.0e-16);
    }
    // Real dynamic stiffness at the old elastic frequency is positive. A
    // damping-only eta(omega0)/2 oscillator would incorrectly leave it zero.
    let at_old = system
        .harmonic_response(omega0, &[C64::from_re(1.0)], 100_000, &gate)
        .unwrap();
    assert!(at_old.displacement_m_sqrt_kg[0].re > 1.0e-7);
    let (mut low, mut high) = (omega0, omega0 * det::sqrt(2.0));
    for _ in 0..60 {
        let middle = low + 0.5 * (high - low);
        let storage = law.modulus_checked(middle).unwrap().0;
        if omega0 * omega0 * storage / law.model.e_inf > middle * middle {
            low = middle;
        } else {
            high = middle;
        }
    }
    let shifted = low + 0.5 * (high - low);
    assert!(shifted > 1.2 * omega0);
    let at_shifted = system
        .harmonic_response(shifted, &[C64::from_re(1.0)], 100_000, &gate)
        .unwrap();
    assert!(
        at_shifted.displacement_m_sqrt_kg[0].re.abs()
            <= 1.0e-12 * at_shifted.displacement_m_sqrt_kg[0].im.abs()
    );
}

#[test]
fn g1_native_prony_work_storage_and_heat_close_on_the_actual_trajectory() {
    let gate = CancelGate::new_clock_free();
    let law = material(vec![(1.6e9, 0.001)]);
    let mut system = bind_prony_memory(elastic(&[1_000.0]), &law, budget(), &gate).unwrap();
    let mut work = 0.0;
    let mut heat = 0.0;
    let mut maximum_material_storage: f64 = 0.0;
    for sample in 0..512 {
        let force = if sample < 64 { 0.2 } else { 0.0 };
        let frame = system.step_under_gate(&[force], &gate).unwrap().clone();
        work += frame.external_work_j;
        heat += frame.material_dissipation_j;
        maximum_material_storage = maximum_material_storage.max(frame.material_energy_j);
        assert!(frame.material_dissipation_j >= 0.0);
        assert_eq!(frame.component_dissipation_j, 0.0);
        assert_eq!(frame.connection_dissipation_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        close(frame.total_energy_j() + heat, work, 2.0e-9, 1.0e-13);
    }
    assert!(maximum_material_storage > 0.0);
    assert!(heat > 0.0);
    assert!(system.material_memory().iter().any(|memory| *memory != 0.0));
    let states = system.components()[0].states().to_vec();
    let memory = system.material_memory().to_vec();
    let samples = system.samples_rendered();
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(system.step_under_gate(&[0.1], &stopped).is_err());
    assert!(system.step_under_gate(&[f64::NAN], &gate).is_err());
    assert_eq!(system.components()[0].states(), states);
    assert_eq!(system.material_memory(), memory);
    assert_eq!(system.samples_rendered(), samples);
}

#[test]
fn g0_zero_prony_arms_recover_the_existing_elastic_pressure_trajectory() {
    let gate = CancelGate::new_clock_free();
    let law = material(Vec::new());
    let mut reference = elastic(&[700.0, 1_300.0]);
    let mut system = bind_prony_memory(elastic(&[700.0, 1_300.0]), &law, budget(), &gate).unwrap();
    for sample in 0..64 {
        let force = [if sample < 12 { 0.3 } else { 0.0 }, -0.1];
        let expected = reference.step(&force).unwrap();
        let actual = system.step_under_gate(&force, &gate).unwrap().clone();
        for (a, b) in system.components()[0]
            .states()
            .iter()
            .zip(reference.states())
        {
            close(
                a.displacement_m_sqrt_kg,
                b.displacement_m_sqrt_kg,
                2.0e-9,
                1.0e-15,
            );
            close(
                a.velocity_m_sqrt_kg_per_s,
                b.velocity_m_sqrt_kg_per_s,
                2.0e-9,
                1.0e-14,
            );
        }
        close(
            actual.observer_pressure_pa,
            expected.observer_pressure_pa,
            2.0e-9,
            1.0e-14,
        );
        close(
            actual.total_energy_j(),
            expected.total_modal_energy_j,
            2.0e-9,
            1.0e-16,
        );
        assert_eq!(actual.material_energy_j, 0.0);
        assert_eq!(actual.material_dissipation_j, 0.0);
        assert!(system.material_memory().is_empty());
    }
}

#[test]
fn g0_prony_admission_refuses_instantaneous_band_exit_and_branch_budget() {
    let gate = CancelGate::new_clock_free();
    let mut law = material(vec![(1.6e9, 0.001)]);
    law.band.1 = 1_100.0; // omega_eq fits, but the instantaneous modulus does not.
    assert!(bind_prony_memory(elastic(&[1_000.0]), &law, budget(), &gate).is_err());
    law.band.1 = 10_000.0;
    let mut limited = budget();
    limited.maximum_material_branches = 0;
    assert!(bind_prony_memory(elastic(&[1_000.0]), &law, limited, &gate).is_err());
    law.model.terms[0].1 = f64::NAN;
    assert!(bind_prony_memory(elastic(&[1_000.0]), &law, budget(), &gate).is_err());
}
