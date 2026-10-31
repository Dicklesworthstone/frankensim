//! G0/G1 source-backed gap/work fixtures; authored numbers, not measurements.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY, ResolvedInductor,
        hysteresis::{
            HysteresisBranch, ScalarHysteresis, UniformHystereticCore,
            shared::{SharedHystereticCore, SharedHystereticWinding},
        },
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{QuantitySpec, Temperature};

const R: [f64; 2] = [2.0, -1.0];
const D: [f64; 2] = [0.3, 0.4];
const BRANCHES: [(f64, f64, f64); 2] = [(0.5, 0.2, 1.0), (0.13, 0.7, 0.35)];
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1.0),
        "{a:.16e} != {b:.16e}"
    );
}
fn medium(mu: f64, gap: f64, turns: u32, varying: bool) -> ResolvedInductor {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
            value: if varying {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(290.0, 0.9 * mu), (310.0, 1.1 * mu)],
                    dims: MAGNETIC_PERMEABILITY_DIMS,
                }
            } else {
                PropertyValue::Scalar {
                    value: mu,
                    dims: MAGNETIC_PERMEABILITY_DIMS,
                }
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: if varying {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored uniform series-gap test".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture gap".into(),
            phase: "fixture".into(),
            process: "numerical".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    ResolvedInductor::resolve(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        0.5,
        gap,
        turns,
    )
    .unwrap()
}
fn model(hysteretic: bool, reference: f64, mu_gap: f64) -> SharedHystereticCore {
    let branches = if hysteretic {
        BRANCHES
            .iter()
            .map(|&(chi, hc, zs)| HysteresisBranch {
                susceptibility_h_m: chi,
                pinning_field_a_m: hc,
                saturation_polarization_t: zs,
            })
            .collect()
    } else {
        vec![]
    };
    let core = UniformHystereticCore::new(
        ScalarHysteresis::new(0.1, 30.0, branches).unwrap(),
        reference,
        1.0,
        1.0,
    )
    .unwrap();
    SharedHystereticCore::new(
        core,
        R.into_iter()
            .zip(D)
            .map(|(r, leakage_h)| SharedHystereticWinding {
                turns_ratio: r / reference,
                leakage_h,
            })
            .collect(),
    )
    .unwrap()
    .with_gap(medium(mu_gap, 0.8, 1, false), [0.1, 3.0])
    .unwrap()
}
fn linear(x: [f64; 2], g: f64, c: f64) -> ([f64; 2], f64, f64) {
    let s = R[0] * x[0] / D[0] + R[1] * x[1] / D[1];
    let k = R[0] * R[0] / D[0] + R[1] * R[1] / D[1];
    let phi = s / (k + 10.0 + c * g);
    let currents = std::array::from_fn(|i| (x[i] - R[i] * phi) / D[i]);
    let energy = 0.5 * (10.0 + c * g) * phi * phi
        + 0.5 * D[0] * currents[0] * currents[0]
        + 0.5 * D[1] * currents[1] * currents[1];
    (currents, phi, energy)
}

#[test]
fn original_gap_source_changes_full_inductance_current_energy_and_force() {
    let x = [1.0, -0.4];
    let mut fluxes = Vec::new();
    for mu in [0.8, 1.6] {
        let model = model(false, 1.0, mu);
        let state = model.initial_state();
        let c = 1.0 / (0.5 * mu);
        let effective = 1.0 / (10.0 + c * 0.8);
        let matrix = model.incremental_inductance_h(&state).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                close(
                    matrix[2 * i + j],
                    R[i] * R[j] * effective + if i == j { D[i] } else { 0.0 },
                    1e-13,
                );
            }
        }
        let next = model.candidate_linkages(&state, &x).unwrap();
        let (currents, phi, energy) = linear(x, 0.8, c);
        close(next.core.flux_linkage_wb, phi, 1e-13);
        close(next.stored_energy_j, energy, 1e-13);
        close(next.magnetic_work_j, energy, 1e-13);
        for i in 0..2 {
            close(next.state.currents_a()[i], currents[i], 1e-13);
            close(next.mean_currents_a[i], 0.5 * currents[i], 1e-13);
        }
        let gap = next.gap.as_ref().unwrap();
        close(gap.force_n, -0.5 * c * phi * phi, 1e-13);
        close(gap.mean_force_n, -c * phi * phi / 6.0, 1e-13);
        assert_eq!(gap.mechanical_work_j, 0.0);
        assert_eq!(next.core.core_heat_j, 0.0);
        assert_eq!(next.state.gap_m(), Some(0.8));
        fluxes.push(phi);
    }
    assert!(fluxes[1] > fluxes[0]);
}

#[test]
fn held_linkages_exchange_gap_work_and_constant_flux_keeps_original_memory() {
    let model = model(false, 1.0, 0.8);
    let x = [1.0, -0.4];
    let before = model
        .candidate_linkages(&model.initial_state(), &x)
        .unwrap()
        .state;
    let next = model.candidate_linkages_and_gap(&before, &x, 1.4).unwrap();
    let (p0, p1) = (linear(x, 0.8, 2.5).1, linear(x, 1.4, 2.5).1);
    let gap = next.gap.as_ref().unwrap();
    close(gap.mean_force_n, -1.25 * p0 * p1, 2e-13);
    close(next.magnetic_work_j, 0.0, 2e-13);
    close(next.stored_change_j, -gap.mechanical_work_j, 2e-13);
    assert!(gap.mean_force_n < 0.0 && next.stored_change_j > 0.0);
    let epsilon = 1e-5;
    let minus = model
        .candidate_linkages_and_gap(&before, &x, 0.8 - epsilon)
        .unwrap();
    let plus = model
        .candidate_linkages_and_gap(&before, &x, 0.8 + epsilon)
        .unwrap();
    close(
        -(plus.stored_energy_j - minus.stored_energy_j) / (2.0 * epsilon),
        -1.25 * p0 * p0,
        2e-9,
    );

    // Exactly representable common-flux chord: delta(s)=delta(R)*Phi.
    let core = UniformHystereticCore::new(
        ScalarHysteresis::new(0.125, 30.0, vec![]).unwrap(),
        1.0,
        1.0,
        1.0,
    )
    .unwrap();
    let prepared = core.candidate_current(&core.initial_state(), 2.0).unwrap();
    let constant = SharedHystereticCore::new(
        core,
        vec![
            SharedHystereticWinding {
                turns_ratio: 1.0,
                leakage_h: 1.0,
            },
            SharedHystereticWinding {
                turns_ratio: -1.0,
                leakage_h: 1.0,
            },
        ],
    )
    .unwrap()
    .with_gap(medium(2.0, 1.0, 1, false), [0.1, 3.0])
    .unwrap();
    let before = constant
        .state_from_core(prepared.material.state, vec![2.5, -0.25])
        .unwrap();
    let after = constant
        .candidate_linkages_and_gap(&before, &[2.625, -0.25], 1.5)
        .unwrap();
    assert_eq!(after.state.core_state(), before.core_state());
    assert_eq!(after.core.core_heat_j, 0.0);
    close(after.gap.as_ref().unwrap().mean_force_n, -0.03125, 1e-15);
    close(after.magnetic_work_j, 2.3125 * 0.125, 1e-15);
}

// Independent monotone bisection of the authored play law at each quadrature
// sample. This does not use production events, slope trees or path moments.
fn point(x: [f64; 2], g: f64, old: &[f64]) -> ([f64; 2], f64, f64, Vec<f64>, f64) {
    let s = R[0] * x[0] / D[0] + R[1] * x[1] / D[1];
    let loading = R[0] * R[0] / D[0] + R[1] * R[1] / D[1] + 2.5 * g;
    let polar = |h: f64| {
        BRANCHES
            .iter()
            .zip(old)
            .map(|(&(chi, hc, zs), &z)| z.max(chi * (h - hc)).min(chi * (h + hc)).clamp(-zs, zs))
            .collect::<Vec<_>>()
    };
    let (mut lo, mut hi) = (-30.0, 30.0);
    for _ in 0..90 {
        let h = f64::midpoint(lo, hi);
        let b = 0.1 * h + polar(h).iter().sum::<f64>();
        if h + loading * b < s {
            lo = h;
        } else {
            hi = h;
        }
    }
    let h = f64::midpoint(lo, hi);
    let z = polar(h);
    let b = 0.1 * h + z.iter().sum::<f64>();
    let currents = std::array::from_fn(|i| (x[i] - R[i] * b) / D[i]);
    let energy = 0.05 * h * h
        + z.iter()
            .zip(BRANCHES)
            .map(|(&z, (chi, _, _))| z * z / (2.0 * chi))
            .sum::<f64>()
        + 1.25 * g * b * b
        + 0.5 * D[0] * currents[0] * currents[0]
        + 0.5 * D[1] * currents[1] * currents[1];
    (currents, -1.25 * b * b, h, z, energy)
}

#[test]
fn moving_gap_integrates_pin_saturation_and_reversal_from_one_original_history() {
    let reversed = model(true, -1.0, 0.8);
    let model = model(true, 1.0, 0.8);
    let mut state = model.initial_state();
    let mut reversed_state = reversed.initial_state();
    for (target, gap) in [([5.0, -2.0], 1.4), ([-4.0, 3.0], 0.35), ([0.8, -0.4], 2.2)] {
        let before = state.clone();
        let x0: [f64; 2] = state.flux_linkages_wb().try_into().unwrap();
        let g0 = state.gap_m().unwrap();
        let old = state.core_state().polarizations_t();
        let expected = point(target, gap, old);
        let next = model
            .candidate_linkages_and_gap(&state, &target, gap)
            .unwrap();
        let other = reversed
            .candidate_linkages_and_gap(&reversed_state, &target, gap)
            .unwrap();
        let mut means = [0.0; 3];
        let n = 8192;
        for j in 0..=n {
            let t = f64::from(j) / f64::from(n);
            let x = std::array::from_fn(|i| x0[i] + t * (target[i] - x0[i]));
            let sample = point(x, g0 + t * (gap - g0), old);
            let weight = if j == 0 || j == n {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            means[0] += weight * sample.0[0];
            means[1] += weight * sample.0[1];
            means[2] += weight * sample.1;
        }
        for value in &mut means {
            *value /= 3.0 * f64::from(n);
        }
        for i in 0..2 {
            close(next.mean_currents_a[i], means[i], 2e-7);
            close(next.mean_currents_a[i], other.mean_currents_a[i], 2e-12);
            close(next.state.currents_a()[i], expected.0[i], 2e-12);
        }
        close(next.gap.as_ref().unwrap().mean_force_n, means[2], 2e-7);
        close(next.state.core_state().field_a_m(), expected.2, 2e-12);
        close(next.stored_energy_j, expected.4, 2e-12);
        let heat = old
            .iter()
            .zip(&expected.3)
            .zip(BRANCHES)
            .map(|((&a, &b), (_, hc, _))| hc * (b - a).abs())
            .sum::<f64>();
        close(next.core.core_heat_j, heat, 2e-12);
        close(next.core.core_heat_j, other.core.core_heat_j, 2e-12);
        let electrical = next
            .mean_currents_a
            .iter()
            .enumerate()
            .map(|(i, &v)| v * (target[i] - x0[i]))
            .sum::<f64>();
        close(
            electrical,
            next.stored_change_j + heat + next.gap.as_ref().unwrap().mechanical_work_j,
            2e-11,
        );
        assert!(next.gap.as_ref().unwrap().mean_force_n <= 0.0);
        assert_eq!(state, before);
        state = next.state;
        reversed_state = other.state;
    }
}

#[test]
fn physical_source_geometry_and_history_cannot_be_replaced_on_refusal() {
    let core = UniformHystereticCore::new(
        ScalarHysteresis::new(0.1, 30.0, vec![]).unwrap(),
        1.0,
        1.0,
        1.0,
    )
    .unwrap();
    let plain = SharedHystereticCore::new(
        core,
        vec![SharedHystereticWinding {
            turns_ratio: 1.0,
            leakage_h: 0.3,
        }],
    )
    .unwrap();
    assert!(
        plain
            .clone()
            .with_gap(medium(0.8, 0.8, 2, false), [0.1, 3.0])
            .is_err()
    );
    assert!(
        plain
            .clone()
            .with_gap(medium(0.8, 0.8, 1, true), [0.1, 3.0])
            .is_err()
    );
    let bound = plain
        .clone()
        .with_gap(medium(0.8, 0.8, 1, false), [0.1, 3.0])
        .unwrap();
    assert!(
        bound
            .clone()
            .with_gap(medium(0.8, 0.8, 1, false), [0.1, 3.0])
            .is_err()
    );
    let state = bound.initial_state();
    let saved = state.clone();
    for gap in [0.1, 3.0, -1.0, f64::NAN] {
        assert!(
            bound
                .candidate_linkages_and_gap(&state, &[1.0], gap)
                .is_err()
        );
        assert_eq!(state, saved);
    }
    assert!(
        bound
            .candidate_linkages_and_gap(&state, &[1e8], 0.8)
            .is_err()
    );
    assert!(plain.candidate_linkages(&state, &[1.0]).is_err());
    assert!(
        plain
            .candidate_linkages_and_gap(&plain.initial_state(), &[1.0], 0.8)
            .is_err()
    );
    assert_eq!(state, saved);
}
