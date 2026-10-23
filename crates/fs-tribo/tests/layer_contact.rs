//! G0/G1/G3 tests of the actual held-effort layer law; all inputs are manufactured.
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use fs_tribo::layer_contact::ThirdBodyLayerContact;

fn law(yield_force: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(1000.0, yield_force, 2.0, 300.0, 0.02).unwrap()
}
fn input(slip: [f64; 2], dt: f64, compliance: f64) -> TangentialStepInput {
    TangentialStepInput { normal_n: 1.0, slip_increment_m: slip, dt_s: dt,
        contact_temperature_k: Some(300.0), nominal_pressure_pa: None,
        compliance_m_per_n: compliance, solve: TangentialSolveConfig {
            max_iterations: 8, maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-12,
        } }
}
fn close(a: f64, b: f64, scale: f64) {
    assert!((a-b).abs() <= 2e-12*scale.max(a.abs()).max(b.abs()), "{a} != {b}");
}

#[test]
fn g1_elastic_loading_and_unloading_return_all_work_without_numerical_heat() {
    let law = law(10.0);
    let initial = law.initial_state();
    let applied = input([0.002, 0.001], 0.01, 0.0001);
    let loaded = law.advance(&initial, &applied).unwrap();
    assert_eq!(loaded.regime, FrictionRegime::Sticking);
    assert_ne!(loaded.slip_increment_m, [0.0; 2], "elastic sticking still stores shear");
    assert_eq!(loaded.dissipated_j, 0.0);
    assert_eq!(loaded.state.plastic_slip_m(), [0.0; 2]);
    for i in 0..2 {
        let force = 1000.0 * applied.slip_increment_m[i] / (2.0 + 1000.0*0.0001);
        close(loaded.traction_n[i], force, 1.0);
    }
    let e = loaded.state.elastic_slip_m();
    let unloaded = law.advance(&loaded.state, &input([-e[0], -e[1]], 0.01, 0.0)).unwrap();
    assert_eq!(unloaded.state.elastic_slip_m(), [0.0; 2]);
    assert_eq!(unloaded.dissipated_j, 0.0);
    assert!(unloaded.stored_change_j < 0.0);
    close(loaded.stored_change_j + unloaded.stored_change_j, 0.0, 0.001);
    close(loaded.closure_residual_j, 0.0, 0.001);
    close(unloaded.closure_residual_j, 0.0, 0.001);
}

#[test]
fn g1_yielded_radial_return_matches_independent_series_spring_bingham_equations() {
    let law = law(0.5);
    let initial = law.initial_state();
    let direction = [0.6, 0.8];
    let x = input([0.01*direction[0], 0.01*direction[1]], 0.01, 0.001);
    let step = law.advance(&initial, &x).unwrap();
    // Solve scalar series equations directly:
    // e = du_free - C F - dt(F-Y)/c and F = k e/2.
    let force = (1000.0*0.01/2.0 + 1000.0*0.01*0.5/(2.0*2.0))
        / (1.0 + 1000.0*0.001/2.0 + 1000.0*0.01/(2.0*2.0));
    let plastic = 0.01*(force-0.5)/2.0;
    assert_eq!(step.regime, FrictionRegime::Sliding);
    for i in 0..2 {
        close(step.traction_n[i], force*direction[i], 1.0);
        close(step.state.plastic_slip_m()[i], plastic*direction[i], 0.01);
        close(step.state.elastic_slip_m()[i], 2.0*force/1000.0*direction[i], 0.01);
    }
    let heat = 0.5*plastic + 2.0/0.01*plastic*plastic;
    close(step.dissipated_j, heat, 0.01);
    assert!(step.dissipated_j > 0.0);
    close(step.closure_residual_j, 0.0, 0.01);

    // The coupled owner checks the final force with its actual solved motion.
    let actual = law.advance(&initial, &input(step.slip_increment_m, x.dt_s, 0.0)).unwrap();
    for i in 0..2 {
        close(actual.traction_n[i], step.traction_n[i], 1.0);
        close(actual.state.elastic_slip_m()[i], step.state.elastic_slip_m()[i], 0.01);
        close(actual.state.plastic_slip_m()[i], step.state.plastic_slip_m()[i], 0.01);
    }
    close(actual.dissipated_j, step.dissipated_j, 0.01);
}

#[test]
fn g1_zero_yield_maxwell_limit_converges_quadratically_without_extra_dissipation() {
    let model = law(0.0);
    let time = 0.01_f64;
    let velocity = 0.001;
    let exact = 2.0/1000.0*velocity*(1.0-(-1000.0*time/2.0_f64).exp());
    let run = |steps: usize| {
        let dt = time/steps as f64;
        let mut state = model.initial_state();
        for _ in 0..steps {
            state = model.advance(&state, &input([velocity*dt, 0.0], dt, 0.0)).unwrap().state;
        }
        (state.elastic_slip_m()[0]-exact).abs()
    };
    let e1 = run(40);
    let e2 = run(80);
    let e3 = run(160);
    assert!(e1/e2 > 3.8 && e1/e2 < 4.2);
    assert!(e2/e3 > 3.8 && e2/e3 < 4.2);
}

#[test]
fn g0_reversing_two_axis_history_closes_physical_storage_and_plastic_work() {
    let model = law(0.1);
    let mut state = model.initial_state();
    let mut work = 0.0;
    let mut heat = 0.0;
    for i in 0..200 {
        let angle = i as f64*0.17;
        let x = input([0.0002*angle.cos(), 0.0002*angle.sin()], 0.0005, 0.0001);
        let old_plastic = state.plastic_slip_m();
        let step = model.advance(&state, &x).unwrap();
        let p = step.state.plastic_slip_m();
        let distance = (p[0]-old_plastic[0]).hypot(p[1]-old_plastic[1]);
        let physical_heat = 0.1*distance + 2.0/x.dt_s*distance*distance;
        close(step.dissipated_j, physical_heat, 1e-5);
        assert!(step.dissipated_j >= 0.0);
        work += step.traction_n[0]*step.slip_increment_m[0]
            + step.traction_n[1]*step.slip_increment_m[1];
        heat += step.dissipated_j;
        close(step.closure_residual_j, 0.0, 1e-5);
        state = step.state;
    }
    close(model.stored_energy_j(&state)+heat, work, work.abs());
    assert!(heat > 0.0);
}

#[test]
fn g0_refusal_preserves_strain_on_opening_temperature_change_or_limit_exhaustion() {
    let model = law(0.5);
    let state = model.advance(&model.initial_state(), &input([0.001, 0.0], 0.001, 0.0)).unwrap().state;
    let snapshot = state.clone();
    let base = input([0.001, 0.0], 0.001, 0.0);
    let mut cases = vec![
        TangentialStepInput { normal_n: 0.0, ..base },
        TangentialStepInput { contact_temperature_k: Some(301.0), ..base },
        TangentialStepInput { contact_temperature_k: None, ..base },
        TangentialStepInput { dt_s: 0.0, ..base },
        TangentialStepInput { compliance_m_per_n: -1.0, ..base },
        TangentialStepInput { slip_increment_m: [f64::NAN, 0.0], ..base },
        TangentialStepInput { slip_increment_m: [1.0, 0.0], ..base },
    ];
    let mut force_limited = base;
    force_limited.solve.maximum_force_n = 1e-6;
    cases.push(force_limited);
    for x in cases {
        assert!(model.advance(&state, &x).is_err());
        assert_eq!(state, snapshot);
    }
    let other = ThirdBodyLayerContact::new(2000.0, 0.5, 2.0, 300.0, 0.02).unwrap();
    assert!(other.advance(&state, &base).is_err());
    assert!(other.stored_energy_j(&state).is_nan());
    assert!(model.advance(&model.initial_state(), &TangentialStepInput { normal_n: 0.0, ..base }).is_err());
}
