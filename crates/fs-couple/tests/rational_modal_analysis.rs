//! G0/G1/G3/G4: complete joint original-state descriptor and physical exchange.
use core::f64::consts::PI;
use fs_couple::mm_line::continuous::rational::{
    ModalRationalLoad, RationalModalLine, RationalModalLineBudget,
    analysis::{
        BeynOptions, PoleOptions, PreparedRationalModalAnalysis, RationalModalAnalysisBudget,
        RationalModalAnalysisError, RationalModalObservation, ResonanceRegion,
    },
};
use fs_duct::RadiusProfile;
use fs_exec::CancelGate;
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};
use fs_vfit::{
    mimo::{
        MimoModel,
        load::{MimoLoadLimits, PreparedMimoLoad},
    },
    model::PoleTerm,
};

const RADIUS: f64 = 0.08;
const LENGTH: f64 = 0.06;
fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn reference() -> f64 {
    air().density * air().sound_speed / (PI * RADIUS * RADIUS)
}
fn budget() -> RationalModalAnalysisBudget {
    RationalModalAnalysisBudget {
        maximum_dimension: 256,
        maximum_payload_bytes: 64 << 20,
        maximum_work: 1_000_000_000,
        reference_rate_per_s: 6000.0,
        relative_residual_tolerance: 1e-10,
    }
}
fn source(ports: usize, direct: bool) -> MimoModel {
    let z = reference();
    let mutual = if ports == 1 {
        vec![0.4]
    } else {
        vec![0.8, 0.18, 0.18, 0.6]
    };
    MimoModel {
        ports,
        poles: vec![PoleTerm::Real {
            pole: -400.0,
            residue: 0.0,
        }],
        residues: vec![mutual.iter().map(|r| C64::from_re(r * z * 400.0)).collect()],
        // For two ports this is genuinely mixed rank, with a null direction
        // different from the mutual memory matrix's principal directions.
        d: vec![if direct { 0.2 * z } else { 0.0 }; ports * ports],
        e: vec![0.0; ports * ports],
        reciprocal: true,
    }
}
fn load(model: &MimoModel, dt: f64, history: bool) -> ModalRationalLoad {
    let mut load = PreparedMimoLoad::new(
        model,
        dt,
        MimoLoadLimits {
            band_rad_s: [1.0, 4e6],
            sampled_passivity_points: 33,
            max_states: 32,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    if history {
        load.step_mean_into(&vec![1e-6; model.ports], &mut vec![0.0; model.ports])
            .unwrap();
    }
    ModalRationalLoad::new(load, RADIUS, air()).unwrap()
}
fn profile() -> RadiusProfile {
    RadiusProfile::exponential(RADIUS, 0.0, LENGTH).unwrap()
}
fn close_complex(a: C64, b: C64, tol: f64) {
    assert!(
        (a - b).abs() < tol,
        "{a:?} != {b:?}, error={:e}",
        (a - b).abs()
    );
}
fn check_power(r: &RationalModalObservation) {
    let scale = r.inlet_power_w.abs()
        + r.duct_mouth_power_w.abs()
        + r.load_power_w.abs()
        + r.incident_power_w
        + r.reflected_power_w
        + 1e-100;
    assert!(
        r.duct_balance_residual_w.abs() < 1e-8 * scale,
        "duct balance {r:?}"
    );
    assert!(r.interface_power_residual_w.abs() < 1e-10 * scale);
    assert!(r.wave_power_residual_w.abs() < 1e-10 * scale);
    assert!(r.relative_equation_residual <= budget().relative_residual_tolerance);
}
fn sinusoid(z: C64) -> (C64, C64) {
    (
        C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh()),
        C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh()),
    )
}
// Independent continuous cylindrical transfer equations, including evanescent
// modes, cross-mode mouth impedance and the original signed volume-flow basis.
fn exact_reflection(model: &MimoModel, gammas: &[f64], omega: f64, incident: &[C64]) -> Vec<C64> {
    let gas = air();
    let p = model.ports;
    let z0 = reference();
    let load = model.eval(C64::new(0.0, -omega));
    let mut lhs = vec![C64::ZERO; p * p];
    let mut pressure_map = vec![C64::ZERO; p * p];
    for k in 0..p {
        let beta =
            C64::from_re((omega / gas.sound_speed).powi(2) - (gammas[k] / RADIUS).powi(2)).sqrt();
        let characteristic = C64::from_re(gas.density * omega / (PI * RADIUS * RADIUS)) / beta;
        let (sin, cos) = sinusoid(beta.scale(LENGTH));
        for j in 0..p {
            let identity = if j == k { C64::ONE } else { C64::ZERO };
            let pressure =
                cos * load[k * p + j] - C64::new(0.0, 1.0) * characteristic * sin * identity;
            let flow =
                C64::new(0.0, -1.0) * sin / characteristic * load[k * p + j] + cos * identity;
            pressure_map[k * p + j] = pressure;
            lhs[k * p + j] = pressure + flow.scale(z0);
        }
    }
    let mut outward: Vec<_> = incident.iter().map(|a| a.scale(2.0)).collect();
    lu_complex(&lhs, p).unwrap().solve(&mut outward);
    (0..p)
        .map(|k| {
            (0..p).fold(C64::ZERO, |a, j| a + pressure_map[k * p + j] * outward[j]) - incident[k]
        })
        .collect()
}

#[test]
fn g1_zero_and_mixed_rank_mouths_converge_to_full_modal_transfer() {
    let gate = CancelGate::new();
    let omega = 6000.0;
    let incident = [C64::new(1.0, 0.2), C64::new(-0.4, 0.3)];
    for direct in [false, true] {
        let model = source(2, direct);
        let boundary = load(&model, 1e-5, false);
        let mut errors = Vec::new();
        for elements in [4, 8, 16] {
            let prepared = PreparedRationalModalAnalysis::new(
                &profile(),
                &air(),
                elements,
                &boundary,
                budget(),
                &gate,
            )
            .unwrap();
            assert_eq!(prepared.dimension(), 4 * (elements + 1) + 2 + 2);
            assert_eq!(prepared.load_dimension(), 2);
            let factor = prepared.factor(omega, &gate).unwrap();
            let response = factor.response(&incident, &gate).unwrap();
            let exact = exact_reflection(
                &model,
                prepared.source().radial_eigenvalues(),
                omega,
                &incident,
            );
            let error = response
                .reflected_pressure_pa
                .iter()
                .zip(&exact)
                .map(|(a, b)| (*a - *b).abs())
                .fold(0.0, f64::max);
            errors.push(error);
            check_power(&response);
            for (&duct, &load) in response
                .duct_mouth_pressure_pa
                .iter()
                .zip(&response.load_pressure_pa)
            {
                close_complex(duct, load, 1e-10);
            }
            // Same LU serves independent excitation columns; mutual response
            // must survive, and refusing a malformed drive does not consume it.
            assert!(factor.response(&incident[..1], &gate).is_err());
            let cross = factor.response(&[C64::ONE, C64::ZERO], &gate).unwrap();
            assert!(cross.reflected_pressure_pa[1].abs() > 1e-7);
            check_power(&cross);
        }
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 1e-3, "{errors:?}");
    }
}

fn poles() -> PoleOptions {
    PoleOptions {
        spectral_scale: 600.0,
        step_tolerance: 1e-12,
        residual_tolerance: 1e-13,
        max_iterations: 32,
        max_dimension: 256,
        max_dense_bytes: 32 << 20,
        max_work: 1_000_000_000,
    }
}
fn contour(nodes: usize) -> BeynOptions {
    BeynOptions {
        nodes,
        probe_columns: 2,
        rank_relative_tolerance: 1e-9,
        count_tolerance: 1e-6,
        quadrature_tolerance: 1e-6,
        contour_min_reciprocal_condition: 1e-12,
        linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-12,
        svd_max_sweeps: 128,
        eigen_max_sweeps: 256,
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_work: 1_000_000_000_000,
        refinement: poles(),
    }
}

#[test]
fn g1_joint_decay_pole_converges_to_loaded_continuous_reflection_pole() {
    let model = source(1, true);
    let boundary = load(&model, 1e-5, false);
    let gate = CancelGate::new();
    // Z(s)+Z0=0 gives the loaded reflection pole. It differs from the bare
    // material pole -400, and its decay is physical signed exterior exchange.
    let expected = -400.0 - model.residues[0][0].re / (model.d[0] + reference());
    let region = ResonanceRegion {
        center: C64::from_re(expected),
        real_radius: 60.0,
        imaginary_radius: 60.0,
    };
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let prepared = PreparedRationalModalAnalysis::new(
            &profile(),
            &air(),
            elements,
            &boundary,
            budget(),
            &gate,
        )
        .unwrap();
        let result = prepared
            .resonance_spectrum(region, &contour(64), &gate)
            .unwrap();
        assert_eq!(result.argument_count, 1);
        assert_eq!(result.resolved_rank, 1);
        let pole = &result.poles[0];
        assert!(pole.laplace_pole_per_s.im.abs() < 1e-7);
        assert!(pole.laplace_pole_per_s.re < 0.0);
        assert!(!pole.oscillation_within_source_band);
        errors.push((pole.laplace_pole_per_s.re - expected).abs());
        check_power(&pole.physical);
        assert!(pole.physical.state[prepared.duct_dimension()].abs() > 0.0);
        let refined = prepared
            .resonance_pole(
                pole.laplace_pole_per_s.scale(1.001),
                &pole.physical.state,
                &poles(),
                &gate,
            )
            .unwrap();
        close_complex(refined.laplace_pole_per_s, pole.laplace_pole_per_s, 1e-6);
    }
    assert!(
        errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 0.01, "{errors:?}");
}

#[test]
fn g0_bare_component_singularities_and_unobservable_memory_are_retained() {
    let mut model = source(1, true);
    model.residues[0][0] = C64::ZERO;
    let boundary = load(&model, 1e-5, false);
    let gate = CancelGate::new();
    let prepared =
        PreparedRationalModalAnalysis::new(&profile(), &air(), 1, &boundary, budget(), &gate)
            .unwrap();
    let mut seed = vec![C64::ZERO; prepared.dimension()];
    seed[prepared.duct_dimension()] = C64::ONE;
    // The exact original rational pole survives despite zero transfer residue.
    // Eliminating or inverting the bare impedance cannot represent this state.
    let pole = prepared
        .resonance_pole(C64::from_re(-399.0), &seed, &poles(), &gate)
        .unwrap();
    close_complex(pole.laplace_pole_per_s, C64::from_re(-400.0), 1e-9);
    assert_eq!(pole.physical.state.len(), 6);
    assert!(pole.physical.inlet_pressure_pa[0].abs() < 1e-12);
    let result = prepared
        .resonance_spectrum(
            ResonanceRegion {
                center: C64::from_re(-400.0),
                real_radius: 30.0,
                imaginary_radius: 30.0,
            },
            &contour(64),
            &gate,
        )
        .unwrap();
    assert_eq!(result.argument_count, 1);
    // Exact nonzero free-free frequency of the one-element consistent FEM.
    // The coupled factor remains regular at this bare duct resonance.
    let bare = 12.0_f64.sqrt() * air().sound_speed / LENGTH;
    let response = prepared
        .factor(bare, &gate)
        .unwrap()
        .response(&[C64::ONE], &gate)
        .unwrap();
    check_power(&response);
}

#[test]
fn g3_clock_history_invariance_and_original_runtime_adapter() {
    let model = source(1, true);
    let gate = CancelGate::new();
    let mut previous: Option<Vec<C64>> = None;
    for (dt, history) in [(1e-5, false), (2e-5, true), (5e-5, true)] {
        let boundary = load(&model, dt, history);
        let original = boundary.load().pole_states().to_vec();
        let prepared =
            PreparedRationalModalAnalysis::new(&profile(), &air(), 2, &boundary, budget(), &gate)
                .unwrap();
        let state = prepared
            .factor(2000.0, &gate)
            .unwrap()
            .response(&[C64::new(1.0, -0.3)], &gate)
            .unwrap()
            .state;
        if let Some(ref last) = previous {
            assert_eq!(last, &state);
        }
        previous = Some(state);
        assert_eq!(original, boundary.load().pole_states());
    }
    let mut runtime = RationalModalLine::new(
        &profile(),
        &air(),
        1,
        load(&model, 1e-5, false),
        RationalModalLineBudget {
            linear: LinearSplitBudget {
                max_modes: 32,
                max_omega_dt: 1.0,
            },
            exact: DenseExponentialBudget {
                max_dimension: 64,
                max_dense_bytes: 32 << 20,
                max_multiply_terms: 1_000_000_000,
            },
            relative_tolerance: 1e-8,
        },
        &gate,
    )
    .unwrap();
    let before = runtime
        .prepare_analysis(budget(), &gate)
        .unwrap()
        .factor(2000.0, &gate)
        .unwrap()
        .response(&[C64::ONE], &gate)
        .unwrap()
        .state;
    runtime.push(&[1.0], &gate).unwrap();
    let accepted = runtime.state().to_vec();
    let after = runtime
        .prepare_analysis(budget(), &gate)
        .unwrap()
        .factor(2000.0, &gate)
        .unwrap()
        .response(&[C64::ONE], &gate)
        .unwrap()
        .state;
    assert_eq!(before, after);
    assert_eq!(accepted, runtime.state());
}

#[test]
fn g4_budgets_band_and_cancellation_refuse_without_source_mutation() {
    let boundary = load(&source(2, true), 1e-5, true);
    let gate = CancelGate::new();
    let state = boundary.load().pole_states().to_vec();
    for b in [
        RationalModalAnalysisBudget {
            maximum_dimension: 1,
            ..budget()
        },
        RationalModalAnalysisBudget {
            maximum_payload_bytes: 1,
            ..budget()
        },
        RationalModalAnalysisBudget {
            maximum_work: 1,
            ..budget()
        },
    ] {
        assert!(
            PreparedRationalModalAnalysis::new(&profile(), &air(), 1, &boundary, b, &gate).is_err()
        );
    }
    let prepared =
        PreparedRationalModalAnalysis::new(&profile(), &air(), 1, &boundary, budget(), &gate)
            .unwrap();
    assert!(prepared.factor(0.5, &gate).is_err());
    assert!(prepared.factor(4e6 + 1.0, &gate).is_err());
    assert!(
        prepared
            .resonance_spectrum(
                ResonanceRegion {
                    center: C64::from_re(-500.0),
                    real_radius: 50.0,
                    imaginary_radius: 50.0
                },
                &BeynOptions {
                    max_work: 1,
                    ..contour(64)
                },
                &gate
            )
            .is_err()
    );
    gate.request();
    assert!(matches!(
        prepared.factor(1000.0, &gate),
        Err(RationalModalAnalysisError::Cancelled)
    ));
    assert!(matches!(
        PreparedRationalModalAnalysis::new(&profile(), &air(), 1, &boundary, budget(), &gate),
        Err(RationalModalAnalysisError::Cancelled)
    ));
    assert_eq!(state, boundary.load().pole_states());
}

#[path = "rational_modal_analysis/regional.rs"]
mod regional;
