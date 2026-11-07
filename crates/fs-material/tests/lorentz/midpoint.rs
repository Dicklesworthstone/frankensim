//! Original source-bank midpoint closure; the numerical leaf has native oracles.
use super::{card, frequency};
use fs_matdb::QueryPoint;
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::lorentz::{LorentzPoleProperties, LorentzState, ResolvedLorentzMaterial},
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{Dims, QuantitySpec, semantic::FrequencyConvention};

fn material(count: usize) -> ResolvedLorentzMaterial {
    let epsilon = QuantitySpec::dimensional(ELECTRIC_PERMITTIVITY_DIMS);
    let angular = frequency(FrequencyConvention::Angular);
    let rate = QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0]));
    let source = card(
        &[
            ("epsilon_inf", epsilon, 2e-11),
            ("delta1", epsilon, 3e-11),
            ("omega1", angular, 1000.0),
            ("gamma1", rate, 25.0),
            ("delta2", epsilon, 5e-11),
            ("omega2", angular, 2500.0),
            ("gamma2", rate, 60.0),
        ],
        73,
    );
    let names = [
        LorentzPoleProperties {
            permittivity_increment: "delta1".into(),
            angular_frequency: "omega1".into(),
            damping_rate: "gamma1".into(),
        },
        LorentzPoleProperties {
            permittivity_increment: "delta2".into(),
            angular_frequency: "omega2".into(),
            damping_rate: "gamma2".into(),
        },
    ];
    ResolvedLorentzMaterial::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names[..count],
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale.max(1e-300),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn states() -> [LorentzState<3>; 2] {
    [
        LorentzState {
            polarization_c_m2: [1e-7, -2e-7, 3e-8],
            current_a_m2: [2e-4, -3e-4, 1e-4],
        },
        LorentzState {
            polarization_c_m2: [-3e-8, 4e-8, 1e-7],
            current_a_m2: [1e-4, 2e-4, -3e-4],
        },
    ]
}

#[test]
fn g1_source_bank_condensation_retains_original_states_work_and_material_identity() {
    let (dt, e0, e1) = (1e-4, [2e4, -3e4, 0.5e4], [-1e4, 2e4, 4e4]);
    for count in [1, 2] {
        let material = material(count);
        let original = material.clone();
        let all_states = states();
        let before = &all_states[..count];
        let prepared = material.prepare_midpoint(dt).unwrap();
        assert_eq!(prepared.duration_s(), dt);
        let epsilon = material.midpoint_permittivity_f_m(dt).unwrap();
        let history = material.midpoint_history_c_m2(dt, e0, before).unwrap();
        assert!(epsilon > material.high_frequency_permittivity_f_m());
        assert_eq!(epsilon, prepared.midpoint_permittivity_f_m());
        assert_eq!(history, prepared.midpoint_history_c_m2(e0, before).unwrap());
        let step = material.step_midpoint(dt, e0, e1, before).unwrap();
        assert_eq!(step, prepared.step_midpoint(e0, e1, before).unwrap());
        assert_eq!(material, original);
        for i in 0..3 {
            near(
                step.displacement_after_c_m2[i],
                epsilon * e1[i] + history[i],
                1e-5,
            );
        }
        let energy = |field: [f64; 3], states: &[LorentzState<3>]| {
            let mut value = 0.5
                * material.high_frequency_permittivity_f_m()
                * field.iter().map(|e| e * e).sum::<f64>();
            for (pole, state) in material.poles().iter().zip(states) {
                for i in 0..3 {
                    value += (state.polarization_c_m2[i].powi(2)
                        + (state.current_a_m2[i] / pole.angular_frequency_rad_s()).powi(2))
                        / (2.0 * pole.permittivity_increment_f_m());
                }
            }
            value
        };
        near(
            step.stored_before_j_m3,
            energy(e0, before),
            step.stored_before_j_m3,
        );
        near(
            step.stored_after_j_m3,
            energy(e1, &step.states_after),
            step.stored_after_j_m3,
        );
        near(
            material
                .stored_energy_density_j_m3(e1, &step.states_after)
                .unwrap(),
            step.stored_after_j_m3,
            step.stored_after_j_m3,
        );
        let work: f64 = (0..3)
            .map(|i| {
                0.5 * (e0[i] + e1[i])
                    * (step.displacement_after_c_m2[i] - step.displacement_before_c_m2[i])
            })
            .sum();
        let mut heat = 0.0;
        for ((pole, old), new) in material.poles().iter().zip(before).zip(&step.states_after) {
            for i in 0..3 {
                let jm = 0.5 * (old.current_a_m2[i] + new.current_a_m2[i]);
                heat += dt * pole.damping_rate_per_s() * jm * jm
                    / (pole.permittivity_increment_f_m() * pole.angular_frequency_rad_s().powi(2));
            }
        }
        assert!(heat > 0.0);
        let scale = step.stored_before_j_m3 + step.stored_after_j_m3 + work.abs() + heat;
        near(step.supplied_work_j_m3, work, scale);
        near(step.heat_j_m3, heat, scale);
        near(
            step.stored_change_j_m3,
            step.stored_after_j_m3 - step.stored_before_j_m3,
            scale,
        );
        near(
            step.stored_change_j_m3 + step.heat_j_m3,
            step.supplied_work_j_m3,
            scale,
        );
        near(step.energy_residual_j_m3, 0.0, scale);
    }
}

#[test]
fn g1_empty_source_bank_is_the_original_lossless_instantaneous_law() {
    let material = material(0);
    let (e0, e1) = ([2.0, -3.0, 4.0], [-1.0, 5.0, 2.0]);
    let prepared = material.prepare_midpoint(0.4).unwrap();
    assert_eq!(prepared.midpoint_permittivity_f_m(), 2e-11);
    assert_eq!(prepared.midpoint_history_c_m2(e0, &[]).unwrap(), [0.0; 3]);
    let step = prepared.step_midpoint(e0, e1, &[]).unwrap();
    assert_eq!(step, material.step_midpoint(0.4, e0, e1, &[]).unwrap());
    assert!(step.states_after.is_empty());
    assert_eq!(step.heat_j_m3, 0.0);
    assert_eq!(step.energy_residual_j_m3, 0.0);
    for i in 0..3 {
        assert_eq!(step.displacement_after_c_m2[i], 2e-11 * e1[i]);
    }
    near(
        step.stored_change_j_m3,
        step.stored_after_j_m3 - step.stored_before_j_m3,
        step.stored_after_j_m3,
    );
    assert_eq!(step.supplied_work_j_m3, step.stored_change_j_m3);
}

#[test]
fn g0_source_bank_refusals_preserve_input_and_cached_replay() {
    let material = material(2);
    let before = states();
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(material.prepare_midpoint(dt).is_err());
        assert!(material.midpoint_permittivity_f_m(dt).is_err());
        assert!(
            material
                .midpoint_history_c_m2(dt, [0.0; 3], &before)
                .is_err()
        );
        assert!(
            material
                .step_midpoint(dt, [0.0; 3], [0.0; 3], &before)
                .is_err()
        );
    }
    let prepared = material.prepare_midpoint(1e-4).unwrap();
    let reference = prepared.step_midpoint([0.0; 3], [1.0; 3], &before).unwrap();
    assert!(
        prepared
            .midpoint_history_c_m2([0.0; 3], &before[..1])
            .is_err()
    );
    assert!(material.displacement_c_m2([0.0; 3], &before[..1]).is_err());
    assert!(
        material
            .stored_energy_density_j_m3([0.0; 3], &before[..1])
            .is_err()
    );
    assert!(
        prepared
            .step_midpoint([0.0; 3], [0.0, 1e300, 0.0], &before)
            .is_err()
    );
    assert_eq!(before, states());
    assert_eq!(
        reference,
        prepared.step_midpoint([0.0; 3], [1.0; 3], &before).unwrap()
    );
    let empty = self::material(0);
    assert!(empty.displacement_c_m2::<0>([], &[]).is_err());
    assert!(empty.step_midpoint(0.1, [0.0; 4], [0.0; 4], &[]).is_err());
}
