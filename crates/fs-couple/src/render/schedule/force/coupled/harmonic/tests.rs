use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode,ModalAcousticTimeBudget};
use crate::render::plate::impact::radiation::{Model,Pole};

fn config()->ModalCouplingConfig {
    ModalCouplingConfig {max_modes:16,max_connections:8,max_setup_terms:1_000_000,
        nyquist_guard_fraction:0.9,maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,
        maximum_abs_connection_force_n:1e6,solve_relative_tolerance:1e-11,
        energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-9}
}
fn mode(omega:f64,zeta:f64)->ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(48_000,vec![ModalAcousticMode {
        angular_frequency_rad_s:omega,damping_ratio:zeta,pressure_per_modal_velocity:C64::ZERO,
    }],ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn close(a:C64,b:C64) {
    assert!((a-b).abs()<=3e-12*a.abs().max(b.abs())+1e-18,"{a:?} != {b:?}");
}
fn power_closes(r:&ModalHarmonicResponse) {
    assert!(r.solve_relative_residual<=1e-11);
    assert!(r.mean_stored_energy_j>=0.0);
    assert!(r.power_residual_w.abs()<=2e-11*r.external_power_w.abs()+1e-18,
        "external={} dissipated={} residual={}",r.external_power_w,r.total_dissipation_w(),r.power_residual_w);
}

#[test]
fn harmonic_two_body_response_matches_independent_dynamic_stiffness_and_power() {
    // G1: direct 2x2 inversion of physical coupled equations, with signed
    // cross terms and complex forcing. Rest length changes only the DC state.
    let gate=CancelGate::new_clock_free();let omega=650.0;
    let mut system=CoupledModalSystem::new_with_columns(vec![mode(500.0,0.02),mode(900.0,0.03)],
        vec![ModalColumnConnection {column:vec![1.0,-0.4],stiffness_n_m:5e4,
            damping_n_s_m:12.0,rest_extension_m:1e-5}],config(),&gate).unwrap();
    system.step(&[0.02,-0.01]).unwrap();
    let before:Vec<_>=system.components().iter().map(|m|m.states().to_vec()).collect();
    let frame=system.last_frame().unwrap().clone();
    let force=[C64::new(0.3,0.1),C64::new(-0.2,0.05)];
    let r=system.harmonic_response(omega,&force,1000,&gate).unwrap();
    let a=C64::new(300_000.0-omega*omega,-32.0*omega);
    let b=C64::new(-20_000.0,4.8*omega);
    let d=C64::new(818_000.0-omega*omega,-55.92*omega);
    let determinant=a*d-b*b;
    close(r.displacement_m_sqrt_kg[0],(d*force[0]-b*force[1])/determinant);
    close(r.displacement_m_sqrt_kg[1],(a*force[1]-b*force[0])/determinant);
    let s=C64::new(0.0,-omega);
    let extension=r.displacement_m_sqrt_kg[0]-r.displacement_m_sqrt_kg[1].scale(0.4);
    close(r.connection_forces_n[0],-C64::new(5e4,-12.0*omega)*extension);
    for (&q,&v) in r.displacement_m_sqrt_kg.iter().zip(&r.velocity_m_sqrt_kg_per_s) { close(v,s*q); }
    assert!(r.component_dissipation_w>0.0 && r.connection_dissipation_w>0.0);
    power_closes(&r);
    assert_eq!(system.last_frame().unwrap(),&frame);
    for (model,expected) in system.components().iter().zip(before) {assert_eq!(model.states(),expected);}
}

fn material_and_air(zeta:f64)->CoupledModalSystem {
    let gate=CancelGate::new_clock_free();
    CoupledModalSystem::new(vec![mode(500.0,0.02)],vec![],config(),&gate).unwrap()
        .with_material_relaxation(vec![fs_phs::RelaxationBranch {
            projection:vec![1.0,0.0],stiffness:3e4,relaxation_time_s:0.002,
        }],Some(vec![0.007]),&gate).unwrap()
        .with_radiation_load(&Model {ports:1,poles:vec![Pole {
            omega:700.0,zeta,coupling:vec![100.0],
        }]},&[0],1,&gate).unwrap()
}

#[test]
fn harmonic_maxwell_and_air_match_analytic_impedances_without_changing_histories() {
    // G1/G3: eliminate the Maxwell ODE and air oscillator independently.
    // Both returned histories and their separate physical loss channels agree.
    let gate=CancelGate::new_clock_free();let mut system=material_and_air(0.1);
    system.step(&[0.03]).unwrap();
    let memory=system.material_memory().to_vec();let air=system.radiation_state().to_vec();
    let state=system.components()[0].states().to_vec();let frame=system.last_frame().unwrap().clone();
    for omega in [70.0,550.0,700.0,1500.0] {
        let s=C64::new(0.0,-omega);let force=C64::new(0.2,0.03);
        let r=system.harmonic_response(omega,&[force],1000,&gate).unwrap();
        let air_den=s*s+s.scale(140.0)+C64::from_re(490_000.0);
        let material=s.scale(60.0)/(C64::from_re(1.0)+s.scale(0.002));
        let dynamic=s*s+s.scale(20.0)+C64::from_re(250_000.0)+material+(s*s).scale(10_000.0)/air_den;
        let q=force/dynamic;let air_q=s.scale(100.0)*q/air_den;
        close(r.displacement_m_sqrt_kg[0],q);
        close(r.material_memory[0],q.scale(30_000.0_f64.sqrt())/(C64::from_re(1.0)+s.scale(0.002)));
        close(r.radiation_state[0],air_q);close(r.radiation_state[1],s*air_q);
        assert!(r.component_dissipation_w>0.0 && r.material_dissipation_w>0.0 && r.radiation_dissipation_w>0.0);
        power_closes(&r);
    }
    assert_eq!(system.material_memory(),memory);assert_eq!(system.radiation_state(),air);
    assert_eq!(system.components()[0].states(),state);assert_eq!(system.last_frame().unwrap(),&frame);
}

#[test]
fn harmonic_bare_lossless_air_pole_is_a_finite_joint_antiresonance() {
    // The isolated load impedance is singular at 700, but the JOINT pencil is
    // nonsingular: q_body=0 and g*p_air=F. Retaining air unknowns is essential.
    let gate=CancelGate::new_clock_free();let system=material_and_air(0.0);
    assert!(system.radiation_model().unwrap().impedance(700.0).is_err());
    let force=C64::new(0.2,0.03);
    let r=system.harmonic_response(700.0,&[force],1000,&gate).unwrap();
    assert!(r.displacement_m_sqrt_kg[0].abs()<1e-18);
    close(r.radiation_state[1],force.scale(0.01));
    close(r.radiation_state[0],force/C64::new(0.0,-70_000.0));
    assert_eq!(r.radiation_dissipation_w,0.0);
    power_closes(&r);
}

#[test]
fn harmonic_maxwell_response_is_the_refined_limit_of_the_original_time_owner() {
    // G1: start on the continuous periodic solution, then use the actual exact
    // held-force owner. Midpoint samples of a sinusoid converge at second order;
    // this independently checks the phasor convention and Maxwell state chart.
    let gate=CancelGate::new_clock_free();
    let arm=fs_phs::RelaxationBranch {projection:vec![1.0,0.0],stiffness:1e6,relaxation_time_s:0.0005};
    let reference=CoupledModalSystem::new(vec![mode(4000.0,0.03)],vec![],config(),&gate).unwrap()
        .with_material_relaxation(vec![arm.clone()],None,&gate).unwrap();
    let omega=2000.0*core::f64::consts::PI;
    let response=reference.harmonic_response(omega,&[C64::from_re(1.0)],1000,&gate).unwrap();
    let q=response.displacement_m_sqrt_kg[0];let v=response.velocity_m_sqrt_kg_per_s[0];
    let z=response.material_memory[0];
    let amplitude=det::sqrt(q.scale(4000.0).norm_sq()+v.norm_sq()+z.norm_sq());
    let mut errors=Vec::new();
    for rate in [8000,16000,32000] {
        let mut model=ModalAcousticTimeModel::try_new(rate,vec![ModalAcousticMode {
            angular_frequency_rad_s:4000.0,damping_ratio:0.03,pressure_per_modal_velocity:C64::ZERO,
        }],ModalAcousticTimeBudget::audible_reference()).unwrap();
        model.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg:q.re,velocity_m_sqrt_kg_per_s:v.re,
        }]).unwrap();
        let mut system=CoupledModalSystem::new(vec![model],vec![],config(),&gate).unwrap()
            .with_material_relaxation(vec![arm.clone()],Some(vec![z.re]),&gate).unwrap();
        for tick in 0..4*(rate/1000) {
            let phase=omega*(f64::from(tick)+0.5)/f64::from(rate);
            system.step(&[phase.cos()]).unwrap();
        }
        let endpoint=system.components()[0].states()[0];
        let dq=4000.0*(endpoint.displacement_m_sqrt_kg-q.re);
        let dv=endpoint.velocity_m_sqrt_kg_per_s-v.re;
        let dz=system.material_memory()[0]-z.re;
        errors.push(det::sqrt(dq*dq+dv*dv+dz*dz)/amplitude);
    }
    assert!(errors[2]<0.005,"{errors:?}");
    for pair in errors.windows(2) {assert!((3.2..4.8).contains(&(pair[0]/pair[1])),"{errors:?}");}
}

#[test]
fn harmonic_refuses_singular_and_unbounded_queries_without_mutation() {
    let gate=CancelGate::new_clock_free();let system=material_and_air(0.1);
    let before=system.material_memory().to_vec();
    assert!(system.harmonic_response(550.0,&[C64::from_re(1.0)],0,&gate).is_err());
    assert!(system.harmonic_response(0.0,&[C64::from_re(1.0)],1000,&gate).is_err());
    assert!(system.harmonic_response(550.0,&[],1000,&gate).is_err());
    assert!(system.harmonic_response(550.0,&[C64::new(f64::NAN,0.0)],1000,&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(system.harmonic_response(550.0,&[C64::ZERO],1000,&cancelled),Err(ModalCouplingError::Cancelled)));
    let singular=CoupledModalSystem::new(vec![mode(500.0,0.0)],vec![],config(),&gate).unwrap();
    assert!(singular.harmonic_response(500.0,&[C64::from_re(1.0)],1000,&gate).is_err());
    let zero=system.harmonic_response(550.0,&[C64::ZERO],1000,&gate).unwrap();
    assert_eq!(zero.mean_stored_energy_j,0.0);assert_eq!(zero.total_dissipation_w(),0.0);
    assert_eq!(system.material_memory(),before);assert_eq!(system.samples_rendered(),0);
}

fn prepared_budget(modes: usize) -> PreparedModalBudget {
    PreparedModalBudget {
        max_modes: modes,
        max_dense_bytes: 2_000_000,
        max_multiply_terms: 100_000_000,
    }
}

#[test]
fn harmonic_constant_loss_resonance_uses_the_coupled_mode_loss_factor() {
    // G1: one joint carries fraction s of the coupled mode's strain energy.
    // Its admitted damping is c_eff = eta*s*omega*G(omega), independent of
    // the old anchoring frequency. Keeping the original dashpot instead would
    // multiply that loss by omega/omega_ref and change the resonance height.
    let gate = CancelGate::new_clock_free();
    let (share, shape, eta, reference) = (0.6, 0.8, 0.04, 900.0);
    let force = C64::new(0.3, -0.1);
    for omega in [300.0, 1700.0, 6000.0] {
        let stiffness = share * omega * omega / (shape * shape);
        let damping = eta * stiffness / reference;
        let system = CoupledModalSystem::new_prepared_with_constant_loss(
            vec![mode(omega * (1.0_f64 - share).sqrt(), 0.0)],
            vec![ModalConnection {
                left: ModalAttachment {
                    component: 0,
                    shapes: vec![shape],
                },
                right: ModalAttachment {
                    component: 0,
                    shapes: vec![0.0],
                },
                stiffness_n_m: stiffness,
                damping_n_s_m: damping,
                rest_extension_m: 0.0,
            }],
            &[ConstantLossConnection {
                connection: 0,
                reference_omega_rad_s: reference,
            }],
            config(),
            prepared_budget(1),
            &gate,
        )
        .unwrap();
        let response = system
            .harmonic_response(omega, &[force], 1000, &gate)
            .unwrap();
        let spectrum = crate::shear_shell::constant_loss::ConstantLossSpectrum::fit().shape(omega);
        let effective_damping = eta * share * omega * spectrum;
        let expected = force / C64::new(0.0, -omega * effective_damping);
        close(response.displacement_m_sqrt_kg[0], expected);
        close(
            response.generalized_connection_forces_n_per_sqrt_kg[0],
            expected.scale(-share * omega * omega)
                - response.velocity_m_sqrt_kg_per_s[0].scale(effective_damping),
        );
        assert!(
            (response.external_power_w / (0.5 * force.norm_sq() / effective_damping) - 1.0).abs()
                < 1e-11
        );
        let old_dashpot_response = force / C64::new(0.0, -omega * damping * shape * shape);
        assert!(
            (response.displacement_m_sqrt_kg[0] - old_dashpot_response).abs()
                > 0.1 * expected.abs()
        );
        power_closes(&response);
        assert_eq!(system.samples_rendered(), 0);
        assert_eq!(
            system.components()[0].states()[0],
            ModalAcousticState::default()
        );
        // One additional rank is charged for the separate damping direction.
        assert!(system.harmonic_response(omega, &[force], 5, &gate).is_err());
        assert!(system.harmonic_response(omega, &[force], 6, &gate).is_ok());
    }
}

fn anchored_pair(rate: u32, initial: Option<&[ModalAcousticState]>) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut model = ModalAcousticTimeModel::try_new(
        rate,
        [(2000.0, 0.006), (9000.0, 0.009)]
            .iter()
            .map(|&(omega, zeta)| ModalAcousticMode {
                angular_frequency_rad_s: omega,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    if let Some(states) = initial {
        model.restore_states(states).unwrap();
    }
    CoupledModalSystem::new_prepared_with_constant_loss(
        vec![model],
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0, -0.7],
            },
            right: ModalAttachment {
                component: 0,
                shapes: vec![0.0, 0.0],
            },
            stiffness_n_m: 4e6,
            damping_n_s_m: 0.04 * 4e6 / 1200.0,
            rest_extension_m: 0.0,
        }],
        &[ConstantLossConnection {
            connection: 0,
            reference_omega_rad_s: 1200.0,
        }],
        config(),
        prepared_budget(2),
        &gate,
    )
    .unwrap()
}

#[test]
fn harmonic_constant_loss_matches_original_time_owner_with_distinct_reaction_directions() {
    // G1/G3: a nonproportional two-mode joint has d=g(K)b NOT parallel to b.
    // Initialize on the continuous periodic orbit, then advance the original
    // prepared time owner under midpoint-held sine forces. Both the endpoint
    // and integrated physical dissipation converge at second order. This
    // catches a harmonic-only copy of the wrong collocated damping law.
    let gate = CancelGate::new_clock_free();
    let reference = anchored_pair(48_000, None);
    let omega = core::f64::consts::TAU * 1000.0;
    let force = [C64::new(0.3, 0.1), C64::new(-0.2, 0.05)];
    let response = reference
        .harmonic_response(omega, &force, 1000, &gate)
        .unwrap();
    let b = [1.0, -0.7];
    let damping_reaction: Vec<_> = b.iter().enumerate().map(|(i, &b)|
        response.generalized_connection_forces_n_per_sqrt_kg[i]
            - response.connection_spring_forces_n[0].scale(b)).collect();
    assert!((damping_reaction[1].scale(b[0]) - damping_reaction[0].scale(b[1])).abs()
        > 0.05 * damping_reaction.iter().map(|f| f.abs()).fold(0.0, f64::max));
    let mut reaction_power = 0.0;
    for (i, (w, zeta)) in [(2000.0, 0.006), (9000.0, 0.009)].into_iter().enumerate() {
        // Original component equation, independently of either connection
        // column: D_i q_i = external_i + complete connection reaction_i.
        close(
            C64::new(w*w - omega*omega, -2.0*omega*zeta*w) * response.displacement_m_sqrt_kg[i],
            force[i] + response.generalized_connection_forces_n_per_sqrt_kg[i],
        );
        reaction_power += 0.5
            * (response.generalized_connection_forces_n_per_sqrt_kg[i]
                * response.velocity_m_sqrt_kg_per_s[i].conj())
            .re;
    }
    assert!(
        (reaction_power + response.connection_dissipation_w).abs()
            < 1e-11 * response.connection_dissipation_w
    );
    power_closes(&response);
    let initial: Vec<_> = response
        .displacement_m_sqrt_kg
        .iter()
        .zip(&response.velocity_m_sqrt_kg_per_s)
        .map(|(q, v)| ModalAcousticState {
            displacement_m_sqrt_kg: q.re,
            velocity_m_sqrt_kg_per_s: v.re,
        })
        .collect();
    let norm = response
        .displacement_m_sqrt_kg
        .iter()
        .zip(&response.velocity_m_sqrt_kg_per_s)
        .zip([2000.0, 9000.0])
        .map(|((q, v), w)| q.scale(w).norm_sq() + v.norm_sq())
        .sum::<f64>()
        .sqrt();
    let mut state_errors = Vec::new();
    let mut loss_errors = Vec::new();
    for rate in [24_000, 48_000, 96_000] {
        let mut system = anchored_pair(rate, Some(&initial));
        let mut dissipation = 0.0;
        let steps = 8 * (rate / 1000);
        for tick in 0..steps {
            let (sin, cos) = det::sin_cos(omega * (f64::from(tick) + 0.5) / f64::from(rate));
            let held = force.map(|f| f.re * cos + f.im * sin);
            dissipation += system.step(&held).unwrap().total_dissipation_j();
        }
        let error = system.components()[0]
            .states()
            .iter()
            .zip(&initial)
            .zip([2000.0, 9000.0])
            .map(|((actual, expected), w)| {
                let q = w * (actual.displacement_m_sqrt_kg - expected.displacement_m_sqrt_kg);
                let v = actual.velocity_m_sqrt_kg_per_s - expected.velocity_m_sqrt_kg_per_s;
                q * q + v * v
            })
            .sum::<f64>()
            .sqrt()
            / norm;
        state_errors.push(error);
        loss_errors.push((dissipation / (0.008 * response.total_dissipation_w()) - 1.0).abs());
    }
    assert!(state_errors[2] < 0.002, "{state_errors:?}");
    assert!(loss_errors[2] < 0.002, "{loss_errors:?}");
    for errors in [&state_errors, &loss_errors] {
        for pair in errors.windows(2) {
            assert!((3.2..4.8).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
    }
}
