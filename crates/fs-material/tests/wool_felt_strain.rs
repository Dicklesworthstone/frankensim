//! Native felt strain, compaction history and constitutive-work regressions.
use fs_material::fiber::{Uniaxial, WoolFelt, WoolFeltState};

fn felt() -> WoolFelt {
    WoolFelt::new(2.0e6, 0.2, 2.5, 3.0, 0.25, 0.7).unwrap()
}

#[test]
fn strain_limit_is_inclusive_and_refusal_preserves_compaction() {
    let law = felt();
    let state = law.response(0.4, &law.initial_state()).unwrap().state;
    let limit = law.response(law.eps_densify, &state).unwrap();
    assert_eq!(limit.state.eps_max, law.eps_densify);
    assert_eq!(limit.stress_pa, law.envelope(law.eps_densify).0);
    assert_eq!(state.eps_max, 0.4); // trial evaluation is not a commit
    for strain in [
        f64::from_bits(law.eps_densify.to_bits() + 1),
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ] {
        assert!(law.response(strain, &state).is_err());
        assert!(law.stress(strain, &state).is_nan());
        assert!(law.tangent(strain, &state).is_nan());
        assert_eq!(law.update_state(strain, &state), state);
    }
}

#[test]
fn compaction_history_must_belong_to_the_original_envelope() {
    let law = felt();
    let state = law.response(0.4, &law.initial_state()).unwrap().state;
    for bad in [
        WoolFeltState {
            eps_max: -0.1,
            sig_max: 0.0,
        },
        WoolFeltState {
            eps_max: 0.0,
            sig_max: 1.0,
        },
        WoolFeltState {
            eps_max: 0.8,
            sig_max: law.envelope(0.8).0,
        },
        WoolFeltState {
            eps_max: f64::NAN,
            ..state.clone()
        },
        WoolFeltState {
            sig_max: f64::INFINITY,
            ..state.clone()
        },
        WoolFeltState {
            sig_max: -1.0,
            ..state.clone()
        },
        WoolFeltState {
            sig_max: 0.9 * state.sig_max,
            ..state.clone()
        },
    ] {
        // Opening a gap must not hide an invalid stored maximum.
        assert!(law.response(-0.1, &bad).is_err());
    }
    let changed_envelope = WoolFelt {
        f_ref: 0.9 * law.f_ref,
        ..law.clone()
    };
    assert!(changed_envelope.response(0.2, &state).is_err());
    for bad in [
        WoolFelt {
            f_ref: f64::INFINITY,
            ..law.clone()
        },
        WoolFelt {
            eps_ref: f64::NAN,
            ..law.clone()
        },
        WoolFelt {
            p: f64::NAN,
            ..law.clone()
        },
        WoolFelt {
            q: law.p - 0.1,
            ..law.clone()
        },
        WoolFelt {
            crush_fraction: 1.0,
            ..law.clone()
        },
        WoolFelt {
            eps_densify: f64::INFINITY,
            ..law.clone()
        },
        WoolFelt {
            eps_densify: law.eps_ref,
            ..law.clone()
        },
    ] {
        assert!(
            WoolFelt::new(
                bad.f_ref,
                bad.eps_ref,
                bad.p,
                bad.q,
                bad.crush_fraction,
                bad.eps_densify
            )
            .is_err()
        );
        assert!(bad.response(0.0, &law.initial_state()).is_err());
    }
}

#[test]
fn residual_gap_reload_and_loading_tangents_match_physical_branches() {
    let law = felt();
    let state = law.response(0.4, &law.initial_state()).unwrap().state;
    let residual = law.eps_residual(&state);
    for strain in [-0.1, 0.0, residual] {
        let response = law.response(strain, &state).unwrap();
        assert_eq!((response.stress_pa, response.tangent_pa), (0.0, 0.0));
        assert_eq!(response.state, state);
    }
    for strain in [0.15, 0.3, 0.5] {
        let response = law.response(strain, &state).unwrap();
        let h = 1e-7;
        let fd = (law.response(strain + h, &state).unwrap().stress_pa
            - law.response(strain - h, &state).unwrap().stress_pa)
            / (2.0 * h);
        assert!((fd / response.tangent_pa - 1.0).abs() < 1e-9);
        assert!(response.stress_pa >= 0.0);
        assert!(response.tangent_pa > 0.0);
    }
    let halfway = law
        .response(0.5 * (residual + state.eps_max), &state)
        .unwrap();
    assert!((halfway.stress_pa / state.sig_max - 0.125).abs() < 1e-15);
}

#[test]
fn tiny_compression_retains_a_representable_tangent_after_stress_underflows() {
    let law = WoolFelt::new(1.0, 0.2, 2.0, 2.0, 0.0, 0.7).unwrap();
    let strain = 1e-200;
    let response = law.response(strain, &law.initial_state()).unwrap();
    assert_eq!(response.stress_pa, 0.0);
    let expected = 2.0 * strain / (law.eps_ref * law.eps_ref);
    assert!((response.tangent_pa / expected - 1.0).abs() < 1e-13);

    let history = law.response(0.4, &law.initial_state()).unwrap().state;
    let reloading = law.response(strain, &history).unwrap();
    assert_eq!(reloading.stress_pa, 0.0);
    assert!((reloading.tangent_pa / expected - 1.0).abs() < 1e-13);
}

#[test]
fn admitted_cycle_work_matches_the_dissipation_formula() {
    for law in [
        felt(),
        WoolFelt::new(2.0e6, 0.2, 2.5, 2.5, 0.0, 0.7).unwrap(),
    ] {
        let peak = 0.5;
        let state = law.response(peak, &law.initial_state()).unwrap().state;
        let steps = 4096;
        let mut loading = 0.0;
        let mut returning = 0.0;
        for i in 0..steps {
            let strain = peak * (f64::from(i) + 0.5) / f64::from(steps);
            loading += law
                .response(strain, &law.initial_state())
                .unwrap()
                .stress_pa;
            returning += law.response(strain, &state).unwrap().stress_pa;
        }
        let work = (loading - returning) * peak / f64::from(steps);
        let expected = state.sig_max
            * peak
            * (1.0 / (law.p + 1.0) - (1.0 - law.crush_fraction) / (law.q + 1.0));
        assert!((work - expected).abs() / (state.sig_max * peak) < 1e-7);
        assert!(work >= -1e-8);
    }
}
