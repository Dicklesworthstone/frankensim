//! Independent piecewise-path checks for one core shared by multiple windings.
//! Coefficients below are authored numerical fixtures, not material measurements.
use fs_material::field_storage::hysteresis::{
    HysteresisBranch, ScalarHysteresis, UniformHystereticCore,
    shared::{SharedHystereticCore, SharedHystereticWinding},
};

const RATIOS: [f64; 2] = [2.0, -1.0];
const LEAKAGE: [f64; 2] = [0.5, 1.0];
const K: f64 = 9.0;

fn core() -> UniformHystereticCore {
    UniformHystereticCore::new(
        ScalarHysteresis::new(
            0.1,
            30.0,
            vec![HysteresisBranch {
                susceptibility_h_m: 0.5,
                pinning_field_a_m: 0.2,
                saturation_polarization_t: 1.0,
            }],
        )
        .unwrap(),
        1.0,
        1.0,
        1.0,
    )
    .unwrap()
}
fn model(ratios: [f64; 2], leakage: [f64; 2]) -> SharedHystereticCore {
    SharedHystereticCore::new(
        core(),
        ratios
            .into_iter()
            .zip(leakage)
            .map(|(turns_ratio, leakage_h)| SharedHystereticWinding {
                turns_ratio,
                leakage_h,
            })
            .collect(),
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-10 * b.abs().max(1.0),
        "{a:.16e} != {b:.16e}"
    );
}
fn polarization(h: f64, old: f64) -> f64 {
    old.max(0.5 * (h - 0.2))
        .min(0.5 * (h + 0.2))
        .clamp(-1.0, 1.0)
}
fn drive(x: [f64; 2]) -> f64 {
    4.0 * x[0] - x[1]
}
fn at_field(h: f64, old: f64) -> f64 {
    h + K * (0.1 * h + polarization(h, old))
}
fn endpoint(x: [f64; 2], old: f64) -> ([f64; 2], f64, f64, f64) {
    let wanted = drive(x);
    let (mut lo, mut hi) = (-30.0, 30.0);
    assert!(at_field(lo, old) <= wanted && wanted <= at_field(hi, old));
    for _ in 0..100 {
        let h = f64::midpoint(lo, hi);
        if at_field(h, old) < wanted {
            lo = h;
        } else {
            hi = h;
        }
    }
    let h = f64::midpoint(lo, hi);
    let z = polarization(h, old);
    let b = 0.1 * h + z;
    let currents = std::array::from_fn(|i| (x[i] - RATIOS[i] * b) / LEAKAGE[i]);
    let energy =
        0.05 * h * h + z * z + 0.25 * currents[0] * currents[0] + 0.5 * currents[1] * currents[1];
    (currents, h, z, energy)
}
// Independently integrate linear currents on each exact source active-set span.
// The implementation under test uses integration by parts, not this quadrature.
fn chord_mean(a: [f64; 2], b: [f64; 2], old: f64) -> [f64; 2] {
    let (s0, ds) = (drive(a), drive(b) - drive(a));
    let mut cuts = vec![0.0, 1.0];
    if ds != 0.0 {
        for h in [2.0 * old - 0.2, 2.0 * old + 0.2, -2.2, -1.8, 1.8, 2.2] {
            let t = (at_field(h, old) - s0) / ds;
            if t > 0.0 && t < 1.0 {
                cuts.push(t);
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let mut mean = [0.0; 2];
    for cut in cuts.windows(2) {
        let x = |t: f64| std::array::from_fn(|i| a[i] + t * (b[i] - a[i]));
        let ia = endpoint(x(cut[0]), old).0;
        let ib = endpoint(x(cut[1]), old).0;
        for i in 0..2 {
            mean[i] += (cut[1] - cut[0]) * f64::midpoint(ia[i], ib[i]);
        }
    }
    mean
}

#[test]
fn one_core_history_closes_multiwinding_work_through_knees_and_reversal() {
    let model = model(RATIOS, LEAKAGE);
    let mut state = model.initial_state();
    for target in [[5.0, -2.0], [-4.0, 3.0], [1.0, 0.0], [0.0, 0.0]] {
        let saved = state.clone();
        let old_z = state.core_state().polarizations_t()[0];
        let previous: [f64; 2] = state.flux_linkages_wb().try_into().unwrap();
        let expected = endpoint(target, old_z);
        let mean = chord_mean(previous, target, old_z);
        let next = model.candidate_linkages(&state, &target).unwrap();
        assert_eq!(state, saved);
        near(next.state.core_state().field_a_m(), expected.1);
        near(next.state.core_state().polarizations_t()[0], expected.2);
        near(next.stored_energy_j, expected.3);
        near(next.core.core_heat_j, 0.2 * (expected.2 - old_z).abs());
        let mut work = 0.0;
        for i in 0..2 {
            near(next.state.currents_a()[i], expected.0[i]);
            near(next.mean_currents_a[i], mean[i]);
            work += mean[i] * (target[i] - previous[i]);
        }
        near(work, next.magnetic_work_j);
        near(work, next.stored_change_j + next.core.core_heat_j);
        near(
            next.stored_change_j,
            next.stored_energy_j - state.stored_energy_j(),
        );
        near(
            2.0 * next.state.currents_a()[0] - next.state.currents_a()[1],
            expected.1,
        );
        state = next.state;
    }
    assert_ne!(state.core_state().polarizations_t(), &[0.0]);
}

#[test]
fn differential_mode_keeps_core_memory_and_inductance_is_full_matrix() {
    let model = model(RATIOS, LEAKAGE);
    let initial = model.initial_state();
    let tangent = model.incremental_inductance_h(&initial).unwrap();
    for (actual, expected) in tangent.into_iter().zip([0.9, -0.2, -0.2, 1.1]) {
        near(actual, expected);
    }
    let before = model
        .candidate_linkages(&initial, &[0.8, 0.2])
        .unwrap()
        .state;
    // r^T D^-1 delta(lambda)=4*.25-1=0: only leakage changes.
    let target = [
        before.flux_linkages_wb()[0] + 0.25,
        before.flux_linkages_wb()[1] + 1.0,
    ];
    let next = model.candidate_linkages(&before, &target).unwrap();
    assert_eq!(next.state.core_state(), before.core_state());
    assert_eq!(next.core.core_heat_j, 0.0);
    for i in 0..2 {
        near(
            next.mean_currents_a[i],
            f64::midpoint(before.currents_a()[i], next.state.currents_a()[i]),
        );
    }
    near(next.magnetic_work_j, next.stored_change_j);
}

#[test]
fn polarity_permutation_and_explicit_prepared_remanence_preserve_physics() {
    let original = model(RATIOS, LEAKAGE);
    let reversed = model([-2.0, -1.0], LEAKAGE);
    let permuted = model([-1.0, 2.0], [1.0, 0.5]);
    let a = original
        .candidate_linkages(&original.initial_state(), &[5.0, -2.0])
        .unwrap();
    let b = reversed
        .candidate_linkages(&reversed.initial_state(), &[-5.0, -2.0])
        .unwrap();
    let c = permuted
        .candidate_linkages(&permuted.initial_state(), &[-2.0, 5.0])
        .unwrap();
    near(a.stored_energy_j, b.stored_energy_j);
    near(a.stored_energy_j, c.stored_energy_j);
    near(a.core.core_heat_j, b.core.core_heat_j);
    near(a.core.core_heat_j, c.core.core_heat_j);
    near(a.mean_currents_a[0], -b.mean_currents_a[0]);
    near(a.mean_currents_a[0], c.mean_currents_a[1]);
    let magnetized = original
        .core()
        .candidate_current(&original.core().initial_state(), 3.0)
        .unwrap();
    let remanent = original
        .core()
        .candidate_current(&magnetized.material.state, 0.0)
        .unwrap();
    let phi = remanent.flux_linkage_wb;
    let linkages = vec![0.25 + 2.0 * phi, 1.0 - phi]; // currents (.5,1), zero total MMF
    let prepared = original
        .state_from_core(remanent.material.state.clone(), linkages.clone())
        .unwrap();
    near(prepared.currents_a()[0], 0.5);
    near(prepared.currents_a()[1], 1.0);
    assert!(prepared.core_state().flux_density_t() > 0.0);
    assert_eq!(
        original.hold_candidate(&prepared).unwrap().core.core_heat_j,
        0.0
    );
    let mut wrong = linkages;
    wrong[0] += 0.5;
    assert!(
        original
            .state_from_core(remanent.material.state, wrong)
            .is_err()
    );
}

#[test]
fn source_and_geometry_refusals_leave_all_winding_and_core_history_unchanged() {
    let model = model(RATIOS, LEAKAGE);
    let state = model.initial_state();
    let saved = state.clone();
    for bad in [vec![1.0], vec![f64::NAN, 0.0], vec![1e5, 0.0]] {
        assert!(model.candidate_linkages(&state, &bad).is_err());
        assert_eq!(state, saved);
    }
    let another = model_with_other_geometry();
    assert!(another.candidate_linkages(&state, &[1.0, 0.0]).is_err());
    for winding in [
        SharedHystereticWinding {
            turns_ratio: 0.0,
            leakage_h: 1.0,
        },
        SharedHystereticWinding {
            turns_ratio: 1.0,
            leakage_h: 0.0,
        },
        SharedHystereticWinding {
            turns_ratio: f64::NAN,
            leakage_h: 1.0,
        },
    ] {
        assert!(SharedHystereticCore::new(core(), vec![winding]).is_err());
    }
    assert!(SharedHystereticCore::new(core(), vec![]).is_err());
}
fn model_with_other_geometry() -> SharedHystereticCore {
    model(RATIOS, [0.6, 1.0])
}
