//! Original molecular law, full spatial modes, independent heat and pole oracles.
use fs_duct::RadiusProfile;
use fs_duct::modal::ModeCount;
use fs_duct::modal::time_domain::{
    ModalEnd, ModalMolecularBudget, assemble_lossless, assemble_lossless_regional,
};
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, PoleOptions};
use fs_thermochem::SpeciesId;

#[path = "modal_molecular_relaxation/field.rs"]
mod field;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn law(count: usize, same_rate: bool) -> RelaxationState {
    let gas = gas();
    let id = SpeciesId::new("synthetic").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes: (0..count)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + i as f64 * 600.0,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: gas.pressure
                        * (2e-5 + if same_rate { 0.0 } else { i as f64 * 3e-4 }),
                    t_ref_k: gas.temperature,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap()
}
fn budget() -> ModalMolecularBudget {
    ModalMolecularBudget {
        maximum_states: 124,
        maximum_payload_bytes: 32 << 20,
        maximum_preparation_work: 1_000_000_000,
    }
}
fn harmonic() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 124,
        max_dense_bytes: 32 << 20,
        max_dense_work: 40_000_000,
        relative_residual_tolerance: 2e-10,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance * a.abs().max(b.abs()).max(1e-30),
        "{a:?} != {b:?}"
    );
}

#[test]
fn g2_complete_modal_memory_matches_original_caloric_weak_equation_and_independent_heat() {
    let law = law(2, false);
    for flare in [0.0, 7.0] {
        let profile = RadiusProfile::exponential(0.05, flare, 0.035).unwrap();
        let image = assemble_lossless(
            &profile,
            &gas(),
            ModeCount(2),
            6,
            ModalEnd::Ports,
            124,
            || false,
        )
        .unwrap()
        .with_molecular_relaxation(law.clone(), budget(), || false)
        .unwrap();
        let source = image.source();
        let n = source.coordinates();
        assert_eq!(image.system().state_dim(), 5 * n);
        assert_eq!(source.molecular_relaxation(), Some(&law));
        assert_eq!(source.molecular_process_indices(), &[0, 1]);
        assert_eq!(source.molecular_state_range(), Some(2 * n..5 * n));
        let model = HarmonicModel::from_ode(image.system(), harmonic(), || false).unwrap();
        let frozen_bulk = gas().density * law.frozen_sound_speed_m_s().powi(2);
        for omega in [100.0, 3000.0, 19000.0, 40000.0] {
            let drive = [
                C64::new(1e-5, 2e-6),
                C64::new(-3e-6, 4e-6),
                C64::new(2e-6, -1e-6),
                C64::new(1e-6, 3e-6),
            ];
            let result = model
                .factor(omega, || false)
                .unwrap()
                .response(&drive, || false)
                .unwrap();
            let compliance = law
                .caloric_response(omega)
                .unwrap()
                .bulk_modulus_pa
                .recip()
                .scale(frozen_bulk);
            let mut matrix = vec![C64::ZERO; n * n];
            for i in 0..n {
                for j in 0..n {
                    matrix[i * n + j] = C64::from_re(source.stiffness()[i * n + j])
                        - compliance.scale(omega * omega * source.mass()[i * n + j]);
                }
            }
            let mut rhs = vec![C64::ZERO; n];
            rhs[..2].copy_from_slice(&drive[..2]);
            rhs[n - 2..].copy_from_slice(&drive[2..]);
            lu_complex(&matrix, n).unwrap().solve(&mut rhs);
            for (port, index) in [0, 1, n - 2, n - 1].into_iter().enumerate() {
                close(
                    result.output[port],
                    C64::new(0.0, -omega) * rhs[index],
                    2e-8,
                );
            }
            let pressures: Vec<_> = (0..n)
                .map(|i| {
                    (0..n).fold(C64::ZERO, |sum, j| {
                        sum + result.state[n + j].scale(source.inverse_mass()[i * n + j])
                    })
                })
                .collect();
            let mut norm = 0.0;
            for i in 0..n {
                for j in 0..n {
                    norm += (pressures[i].conj() * pressures[j]).re * source.mass()[i * n + j];
                }
            }
            let expected = 0.5 * omega * compliance.im * norm;
            assert!(expected > 0.0);
            let heat = source
                .molecular_mean_heat_per_process_w(&result.state)
                .unwrap();
            assert_eq!(heat.len(), 2);
            assert!(heat.iter().all(|h| *h > 0.0));
            close(
                C64::from_re(heat.iter().sum()),
                C64::from_re(expected),
                2e-8,
            );
            close(
                C64::from_re(result.dissipated_power_w),
                C64::from_re(expected),
                2e-8,
            );
            assert!(result.power_residual_w.abs() < 2e-9 * expected);
            let real: Vec<_> = result.state.iter().map(|z| z.re).collect();
            let imaginary: Vec<_> = result.state.iter().map(|z| z.im).collect();
            let a = source.molecular_heat_per_process_w(&real).unwrap();
            let b = source.molecular_heat_per_process_w(&imaginary).unwrap();
            for j in 0..2 {
                close(
                    C64::from_re(heat[j]),
                    C64::from_re(0.5 * (a[j] + b[j])),
                    2e-12,
                );
            }
        }
    }
}

// Independent molecular bulk response at complex s; not the assembled generator.
fn bulk(law: &RelaxationState, s: C64) -> (C64, C64) {
    let mut cv = C64::from_re(law.frozen_cv_j_kg_k());
    let mut derivative = C64::ZERO;
    for p in law.processes() {
        let denominator = C64::ONE + s.scale(p.relaxation_time_s);
        cv = cv + C64::from_re(p.heat_capacity_j_kg_k) / denominator;
        derivative = derivative
            - C64::from_re(p.heat_capacity_j_kg_k * p.relaxation_time_s)
                / (denominator * denominator);
    }
    let gas = law.gas();
    let k = (C64::ONE + C64::from_re(gas.specific_gas_constant) / cv).scale(gas.pressure);
    let dk = -(derivative / (cv * cv)).scale(gas.pressure * gas.specific_gas_constant);
    (k, dk)
}
#[test]
fn g1_molecular_loaded_resonance_converges_to_original_complex_caloric_root() {
    let law = law(2, false);
    let length = 0.25;
    let beta = core::f64::consts::PI / (2.0 * length);
    let mut exact = C64::new(-10.0, -gas().sound_speed * beta);
    for _ in 0..30 {
        let (k, dk) = bulk(&law, exact);
        let f = exact * exact + k.scale(beta * beta / gas().density);
        let df = exact.scale(2.0) + dk.scale(beta * beta / gas().density);
        exact = exact - f / df;
    }
    assert!(exact.re < -1.0);
    let mut errors = Vec::new();
    for elements in [4, 8, 16] {
        let profile = RadiusProfile::exponential(0.05, 0.0, length).unwrap();
        let image = assemble_lossless(
            &profile,
            &gas(),
            ModeCount(1),
            elements,
            ModalEnd::PressureRelease,
            124,
            || false,
        )
        .unwrap()
        .with_molecular_relaxation(law.clone(), budget(), || false)
        .unwrap();
        let model = HarmonicModel::from_ode(image.system(), harmonic(), || false).unwrap();
        let frequency = exact.im.abs() * 0.98;
        let state = model
            .factor(frequency, || false)
            .unwrap()
            .response(&[C64::ONE], || false)
            .unwrap()
            .state;
        let pole = model
            .resonance_pole(
                C64::new(exact.re, -frequency),
                &state,
                &PoleOptions {
                    spectral_scale: exact.abs(),
                    step_tolerance: 1e-11,
                    residual_tolerance: 1e-11,
                    max_iterations: 40,
                    max_dimension: 124,
                    max_dense_bytes: 32 << 20,
                    max_work: 100_000_000,
                },
                || false,
            )
            .unwrap();
        assert!(pole.laplace_pole_per_s.re < 0.0);
        assert!(pole.dissipated_power_w > 0.0);
        assert!(
            (2.0 * pole.laplace_pole_per_s.re * pole.stored_energy_j + pole.dissipated_power_w)
                .abs()
                < 1e-8 * pole.dissipated_power_w
        );
        errors.push((pole.laplace_pole_per_s - exact).abs() / exact.abs());
    }
    assert!(errors[0] > 1e-3, "{errors:?}");
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 5e-4, "{errors:?}");
}

#[test]
fn g3_empty_law_is_bitwise_and_equal_rate_original_reservoirs_stay_distinct() {
    let profile = RadiusProfile::exponential(0.03, 3.0, 0.04).unwrap();
    let base = assemble_lossless(
        &profile,
        &gas(),
        ModeCount(2),
        4,
        ModalEnd::Ports,
        124,
        || false,
    )
    .unwrap();
    let attached = assemble_lossless(
        &profile,
        &gas(),
        ModeCount(2),
        4,
        ModalEnd::Ports,
        124,
        || false,
    )
    .unwrap()
    .with_molecular_relaxation(law(0, false), budget(), || false)
    .unwrap();
    assert_eq!(base.system().structure(), attached.system().structure());
    assert_eq!(base.source().mass(), attached.source().mass());
    assert_eq!(
        base.source().inverse_mass(),
        attached.source().inverse_mass()
    );
    assert_eq!(
        base.system().storage().quadratic_part().unwrap().matrix,
        attached.system().storage().quadratic_part().unwrap().matrix
    );
    let original = law(2, true);
    let regional = assemble_lossless_regional(
        &profile,
        &[gas(); 4],
        ModeCount(2),
        ModalEnd::Ports,
        124,
        || false,
    )
    .unwrap()
    .with_molecular_relaxation(original.clone(), budget(), || false)
    .unwrap();
    assert_eq!(regional.source().molecular_process_indices(), &[0, 1]);
    assert_eq!(
        regional.source().state_dimension(),
        5 * regional.source().coordinates()
    );
    assert_eq!(regional.source().molecular_relaxation(), Some(&original));
    assert_eq!(
        regional.source().regional_gases(),
        Some([gas(); 4].as_slice())
    );
}

#[test]
fn g4_full_source_capacity_material_and_cancellation_refuse_without_truncation() {
    let profile = RadiusProfile::exponential(0.03, 0.0, 0.04).unwrap();
    let make = || {
        assemble_lossless(
            &profile,
            &gas(),
            ModeCount(2),
            4,
            ModalEnd::Ports,
            124,
            || false,
        )
        .unwrap()
    };
    for limits in [
        ModalMolecularBudget {
            maximum_states: 49,
            ..budget()
        },
        ModalMolecularBudget {
            maximum_payload_bytes: 1,
            ..budget()
        },
        ModalMolecularBudget {
            maximum_preparation_work: 1,
            ..budget()
        },
    ] {
        assert!(
            make()
                .with_molecular_relaxation(law(2, false), limits, || false)
                .is_err()
        );
    }
    let accepted = make()
        .with_molecular_relaxation(
            law(2, false),
            ModalMolecularBudget {
                maximum_states: 50,
                ..budget()
            },
            || false,
        )
        .unwrap();
    assert_eq!(accepted.system().state_dim(), 50);
    assert!(
        accepted
            .with_molecular_relaxation(law(2, false), budget(), || false)
            .is_err()
    );
    let warm = GasState::try_new(&GasSpec::dry_air_ussa1976(), 330.0, 101325.0).unwrap();
    let nonuniform = assemble_lossless_regional(
        &profile,
        &[gas(), warm],
        ModeCount(2),
        ModalEnd::Ports,
        124,
        || false,
    )
    .unwrap();
    assert!(
        nonuniform
            .with_molecular_relaxation(law(2, false), budget(), || false)
            .is_err()
    );
    assert!(
        make()
            .with_molecular_relaxation(law(2, false), budget(), || true)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        make()
            .with_molecular_relaxation(law(2, false), budget(), || {
                polls += 1;
                polls == 4
            })
            .is_err()
    );
}
