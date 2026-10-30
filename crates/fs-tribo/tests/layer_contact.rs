//! G0/G1/G3 tests of the actual held-effort layer law; all inputs are manufactured.
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use fs_tribo::layer_contact::{ThirdBodyLayerContact, ThirdBodyLayerRelaxationArm};

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

#[test]
fn thermal_continuation_preserves_strain_and_releases_only_actual_elastic_energy() {
    let cold = law(0.5);
    let loaded = cold.advance(&cold.initial_state(), &input([0.01, 0.002], 0.01, 0.0)).unwrap();
    assert!(loaded.state.elastic_slip_m()[0] > 0.0 && loaded.state.plastic_slip_m()[0] > 0.0);
    let warm = ThirdBodyLayerContact::new(1000.0, 0.1, 0.8, 310.0, 0.02).unwrap();
    let continued = warm.continue_state(&loaded.state).unwrap();
    assert_eq!(continued.elastic_slip_m(), loaded.state.elastic_slip_m());
    assert_eq!(continued.plastic_slip_m(), loaded.state.plastic_slip_m());
    assert_eq!(warm.stored_energy_j(&continued), cold.stored_energy_j(&loaded.state));
    let mut held = input([0.0; 2], 0.01, 0.0);
    held.contact_temperature_k = Some(310.0);
    let relaxed = warm.advance(&continued, &held).unwrap();
    assert!(relaxed.dissipated_j > 0.0 && relaxed.stored_change_j < 0.0);
    close(relaxed.stored_change_j + relaxed.dissipated_j, 0.0, 0.01);
    assert!(relaxed.state.plastic_slip_m()[0] > continued.plastic_slip_m()[0]);

    // Temperature is no permission to alter the elastic potential or erase
    // stored work, and a changed law still cannot consume a foreign history.
    assert!(warm.advance(&loaded.state, &held).is_err());
    for altered in [
        ThirdBodyLayerContact::new(1001.0, 0.1, 0.8, 310.0, 0.02).unwrap(),
        ThirdBodyLayerContact::new(1000.0, 0.1, 0.8, 310.0, 0.01).unwrap(),
    ] {
        assert!(altered.continue_state(&loaded.state).is_err());
    }
    let saved = continued.clone();
    held.normal_n = 0.0;
    assert!(warm.advance(&continued, &held).is_err());
    assert_eq!(continued, saved);
}

#[test]
fn g1_affine_layer_thermodynamics_preserves_strain_and_distinguishes_free_from_internal_energy() {
    // k(T)=1600-2*T. This manufactured free-energy law is not a rosin fit.
    let cold = ThirdBodyLayerContact::new(1000.0, 100.0, 2.0, 300.0, 0.02).unwrap();
    let loaded = cold.advance(&cold.initial_state(), &input([0.002, 0.001], 0.01, 0.0)).unwrap();
    let warm = ThirdBodyLayerContact::new(960.0, 100.0, 2.0, 320.0, 0.02).unwrap();
    let transfer = warm.continue_state_with_work(&loaded.state).unwrap();
    assert_eq!(transfer.state.elastic_slip_m(), loaded.state.elastic_slip_m());
    assert_eq!(transfer.state.plastic_slip_m(), loaded.state.plastic_slip_m());
    let before = cold.thermodynamics(&loaded.state, -2.0).unwrap();
    let after = warm.thermodynamics(&transfer.state, -2.0).unwrap();
    close(before.free_energy_j, 0.0025, 0.001);
    close(before.entropy_j_k, 5e-6, 1e-6);
    close(before.internal_energy_j, 0.004, 0.001);
    close(after.internal_energy_j, before.internal_energy_j, 0.001);
    close(transfer.coefficient_free_energy_change_j, -0.0001, 0.001);
    close(after.free_energy_j-before.free_energy_j,
        transfer.coefficient_free_energy_change_j, 0.001);
    close(transfer.coefficient_free_energy_change_j
        +320.0*after.entropy_j_k-300.0*before.entropy_j_k, 0.0, 0.001);
    assert!(warm.continue_state(&loaded.state).is_err());
    assert!(warm.thermodynamics(&loaded.state, -2.0).is_err());
    assert!(warm.thermodynamics(&transfer.state, f64::NAN).is_err());
}

#[test]
fn g1_affine_layer_elastic_loading_cools_and_closes_first_law_without_plastic_heat() {
    let slope = -2.0;
    let capacity = 0.01; // finite, independent caloric capacity [J/K]
    let e2 = 0.002_f64.powi(2) + 0.001_f64.powi(2);
    let t1 = 300.0/(1.0-slope*e2/(2.0*capacity));
    let cold = ThirdBodyLayerContact::new(1000.0, 100.0, 2.0, 300.0, 0.02).unwrap();
    let warm = ThirdBodyLayerContact::new(1600.0+slope*t1, 100.0, 2.0, t1, 0.02).unwrap();
    let transfer = warm.continue_state_with_work(&cold.initial_state()).unwrap();
    let mut drive = input([0.002, 0.001], 0.01, 0.0);
    drive.contact_temperature_k = Some(t1);
    let step = warm.advance(&transfer.state, &drive).unwrap();
    let end = warm.thermodynamics(&step.state, slope).unwrap();
    let heat = -transfer.coefficient_free_energy_change_j - t1*end.entropy_j_k;
    let work = step.traction_n[0]*step.slip_increment_m[0]
        +step.traction_n[1]*step.slip_increment_m[1];
    assert!(t1 < 300.0 && heat < 0.0);
    assert_eq!(step.dissipated_j, 0.0);
    close(capacity*(t1-300.0), heat, 0.001);
    close(end.internal_energy_j + capacity*(t1-300.0), work, 0.001);
    // Backward endpoint temperature gives nonnegative O(dt^2) entropy error;
    // neither that error nor reversible cooling is material dissipation.
    let entropy_change = capacity*(t1/300.0).ln() + end.entropy_j_k;
    assert!(entropy_change >= -1e-15 && entropy_change < 2e-9);
}

fn hereditary(yield_force: f64) -> ThirdBodyLayerContact {
    law(yield_force).with_relaxation(vec![
        ThirdBodyLayerRelaxationArm { stiffness_n_m: 2000.0, relaxation_time_s: 0.01 },
        ThirdBodyLayerRelaxationArm { stiffness_n_m: 700.0, relaxation_time_s: 0.04 },
    ]).unwrap()
}

#[test]
fn g1_affine_maxwell_moduli_preserve_full_entropy_internal_energy_and_heat_work_identity() {
    // Independent complete free energy: k_inf=1600-2T,
    // k_0=2300-T, k_1=1300-2T. Internal-energy coefficients are their
    // intercepts, so a temperature change at fixed strains leaves U unchanged.
    let cold = hereditary(100.0);
    let loaded = cold.advance(&cold.initial_state(), &input([0.002, 0.001], 0.001, 0.0)).unwrap();
    let warm = ThirdBodyLayerContact::new(960.0, 100.0, 2.0, 320.0, 0.02).unwrap()
        .with_relaxation(vec![
            ThirdBodyLayerRelaxationArm { stiffness_n_m: 1980.0, relaxation_time_s: 0.01 },
            ThirdBodyLayerRelaxationArm { stiffness_n_m: 660.0, relaxation_time_s: 0.04 },
        ]).unwrap();
    let transfer = warm.continue_state_with_work(&loaded.state).unwrap();
    assert_eq!(transfer.state.relaxation_slip_m(), loaded.state.relaxation_slip_m());
    assert_eq!(transfer.state.elastic_slip_m(), loaded.state.elastic_slip_m());
    assert_eq!(transfer.state.plastic_slip_m(), loaded.state.plastic_slip_m());
    let norm2 = |x: [f64; 2]| x[0]*x[0] + x[1]*x[1];
    let e2 = norm2(loaded.state.elastic_slip_m());
    let r0 = norm2(loaded.state.relaxation_slip_m()[0]);
    let r1 = norm2(loaded.state.relaxation_slip_m()[1]);
    let entropy = 0.5*(2.0*e2+r0+2.0*r1);
    let internal = 0.5*(1600.0*e2+2300.0*r0+1300.0*r1);
    let before = cold.thermodynamics_with_relaxation(&loaded.state, -2.0, &[-1.0, -2.0]).unwrap();
    let continued = warm.thermodynamics_with_relaxation(&transfer.state, -2.0, &[-1.0, -2.0]).unwrap();
    close(before.entropy_j_k, entropy, entropy);
    close(before.internal_energy_j, internal, internal);
    close(continued.internal_energy_j, internal, internal);
    close(transfer.coefficient_free_energy_change_j, -20.0*entropy, 20.0*entropy);
    close(continued.free_energy_j-before.free_energy_j,
        transfer.coefficient_free_energy_change_j, internal);
    let step = warm.advance(&transfer.state, &TangentialStepInput {
        contact_temperature_k: Some(320.0), ..input([0.0002, -0.0001], 0.001, 0.0)
    }).unwrap();
    let after = warm.thermodynamics_with_relaxation(&step.state, -2.0, &[-1.0, -2.0]).unwrap();
    let heat = step.dissipated_j - transfer.coefficient_free_energy_change_j
        - (320.0*after.entropy_j_k - 300.0*before.entropy_j_k);
    let work = step.traction_n[0]*step.slip_increment_m[0]
        + step.traction_n[1]*step.slip_increment_m[1];
    close(after.internal_energy_j-before.internal_energy_j+heat, work, internal);
    assert!(step.dissipated_j > 0.0);
    assert_eq!(step.state.plastic_slip_m(), [0.0; 2]);
    assert!(cold.thermodynamics_with_relaxation(&loaded.state, -2.0, &[-1.0]).is_err());
    assert!(cold.thermodynamics_with_relaxation(&loaded.state, -2.0, &[-1.0, f64::NAN]).is_err());
    assert!(warm.thermodynamics_with_relaxation(&loaded.state, -2.0, &[-1.0, -2.0]).is_err());
    assert_eq!(cold.thermodynamics(&loaded.state, -2.0).unwrap(),
        cold.thermodynamics_with_relaxation(&loaded.state, -2.0, &[0.0, 0.0]).unwrap());
}

#[test]
fn g1_subyield_shear_relaxes_to_equilibrium_with_second_order_material_memory() {
    let model = hereditary(100.0);
    let loaded = model.advance(&model.initial_state(),
        &input([0.002, 0.001], 0.0002, 0.0)).unwrap();
    assert_eq!(loaded.regime, FrictionRegime::Sticking);
    assert!(loaded.dissipated_j > 0.0, "subyield viscoelastic loss is physical");
    assert_eq!(loaded.state.plastic_slip_m(), [0.0; 2]);
    let run = |steps: usize| {
        let mut state = loaded.state.clone();
        let mut heat = 0.0;
        let dt = 0.02/steps as f64;
        for _ in 0..steps {
            let step = model.advance(&state, &input([0.0; 2], dt, 0.0)).unwrap();
            let channels = model.step_dissipation(&state, &step.state, dt).unwrap();
            assert_eq!(channels.plastic_j, 0.0);
            close(channels.relaxation_j, step.dissipated_j, 1e-5);
            close(step.closure_residual_j, 0.0, 1e-5);
            heat += step.dissipated_j;
            state = step.state;
        }
        assert_eq!(state.elastic_slip_m(), loaded.state.elastic_slip_m());
        assert_eq!(state.plastic_slip_m(), [0.0; 2]);
        close(model.stored_energy_j(&state)+heat,
            model.stored_energy_j(&loaded.state), 0.01);
        // The independent continuous held-strain Maxwell solution is one
        // exponential per SOURCE arm. No force/spectrum fitting is involved.
        let mut error = 0.0_f64;
        for ((arm, initial), final_slip) in model.relaxation().iter()
            .zip(loaded.state.relaxation_slip_m()).zip(state.relaxation_slip_m()) {
            for axis in 0..2 {
                let exact = initial[axis]*(-0.02/arm.relaxation_time_s).exp();
                error = error.max((final_slip[axis]-exact).abs());
            }
        }
        error
    };
    let errors = [run(20), run(40), run(80)];
    for pair in errors.windows(2) { assert!((3.8..4.2).contains(&(pair[0]/pair[1]))); }
}

#[test]
fn g1_maxwell_and_plastic_losses_close_the_same_two_axis_contact_return() {
    let model = hereditary(0.1);
    let mut state = model.initial_state();
    let mut heat = 0.0;
    let mut work = 0.0;
    for i in 0..120 {
        let angle = i as f64*0.19;
        let drive = input([0.0002*angle.cos(), 0.0002*angle.sin()], 0.0005, 0.0001);
        let step = model.advance(&state, &drive).unwrap();
        let channels = model.step_dissipation(&state, &step.state, drive.dt_s).unwrap();
        close(channels.plastic_j+channels.relaxation_j, step.dissipated_j, 1e-5);
        assert!(channels.plastic_j >= 0.0 && channels.relaxation_j >= 0.0);
        let actual = model.advance(&state, &input(step.slip_increment_m, drive.dt_s, 0.0)).unwrap();
        for axis in 0..2 { close(actual.traction_n[axis], step.traction_n[axis], 1.0); }
        for ((arm, old), new) in model.relaxation().iter()
            .zip(state.relaxation_slip_m()).zip(step.state.relaxation_slip_m()) {
            for axis in 0..2 {
                let de = step.state.elastic_slip_m()[axis]-state.elastic_slip_m()[axis];
                close(new[axis]-old[axis], de-drive.dt_s/arm.relaxation_time_s
                    * (0.5*new[axis]+0.5*old[axis]), 0.001);
            }
        }
        work += step.traction_n[0]*step.slip_increment_m[0]
            +step.traction_n[1]*step.slip_increment_m[1];
        heat += step.dissipated_j;
        close(step.closure_residual_j, 0.0, 1e-5);
        state = step.state;
    }
    assert!(heat > 0.0);
    close(model.stored_energy_j(&state)+heat, work, work.abs());
}

#[test]
fn g0_material_memory_continues_thermal_history_and_refuses_unresolved_or_replaced_arms() {
    let cold = hereditary(100.0);
    let loaded = cold.advance(&cold.initial_state(), &input([0.002, 0.001], 0.001, 0.0)).unwrap();
    let warm = ThirdBodyLayerContact::new(1000.0, 100.0, 2.0, 310.0, 0.02).unwrap()
        .with_relaxation(cold.relaxation().iter().map(|arm| ThirdBodyLayerRelaxationArm {
            stiffness_n_m: arm.stiffness_n_m, relaxation_time_s: 0.5*arm.relaxation_time_s,
        }).collect()).unwrap();
    let continued = warm.continue_state(&loaded.state).unwrap();
    assert_eq!(continued.relaxation_slip_m(), loaded.state.relaxation_slip_m());
    assert_eq!(warm.stored_energy_j(&continued), cold.stored_energy_j(&loaded.state));
    let held = TangentialStepInput { contact_temperature_k: Some(310.0), ..input([0.0;2], 0.001, 0.0) };
    let relaxed = warm.advance(&continued, &held).unwrap();
    assert!(relaxed.dissipated_j > 0.0);
    close(relaxed.stored_change_j+relaxed.dissipated_j, 0.0, 0.01);
    let snapshot = continued.clone();
    assert!(warm.advance(&continued, &TangentialStepInput { dt_s: 0.02, ..held }).is_err());
    assert_eq!(continued, snapshot);
    assert!(law(100.0).continue_state_with_work(&loaded.state).is_err());
    for arm in [
        ThirdBodyLayerRelaxationArm { stiffness_n_m: 0.0, relaxation_time_s: 0.01 },
        ThirdBodyLayerRelaxationArm { stiffness_n_m: 1.0, relaxation_time_s: f64::NAN },
    ] { assert!(law(1.0).with_relaxation(vec![arm]).is_err()); }
    assert!(law(1.0).with_relaxation(vec![cold.relaxation()[0]; 65]).is_err());
}
