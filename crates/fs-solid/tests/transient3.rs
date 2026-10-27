//! G1/G3/G4 inertial P1 response, physical work and transactional checkpoints.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::transient::{TransientConfig3, TransientError3};

#[path = "transient3/thermal_damping.rs"]
mod thermal_damping;

const NODES: [[f64; 3]; 4] = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]];
const TETS: [[usize; 4]; 1] = [[0,1,2,3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free(); if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 73, kernel_id: 73, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic)))
}
fn problem<'a>(material: &'a TetElasticMaterial, fixed: &'a [usize]) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem { nodes_m: &NODES, tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material), fixed_dofs: fixed, budget: TetAssemblyBudget::default() }
}
fn geometry(cx: &Cx<'_>) -> TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::default(), 4, cx).unwrap()
}
fn config() -> TransientConfig3 {
    TransientConfig3 { maximum_displacement_m: 1.0, ..Default::default() }
}

#[test]
fn g1_single_free_component_has_consistent_mass_frequency_and_second_order() {
    with_cx(false, |cx| {
        // E=60, nu=0, rho=6, V=1/6: K_33=10 N/m, M_33=0.1 kg.
        // This is a physical fully restrained tet with its x-node x DOF free.
        let material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material, &fixed); let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let mut initial = ZERO; initial[1][0] = 0.001;
        let mut errors = Vec::new();
        for count in [25,50,100] {
            let mut run = load.begin_transient(&problem,None,0.0,&initial,&ZERO,config(),cx).unwrap();
            assert!((run.state().acceleration[0] + 0.1).abs() < 1e-13);
            assert_eq!(run.free_dofs(), &[3]);
            let energy = run.state().total_energy_j;
            for _ in 0..count {
                let report = run.advance(&load,0.5/count as f64,cx).unwrap();
                assert!(report.balance_error_j.abs() < 1e-13);
                assert!((report.quadratic_energy_j-energy).abs() < 1e-13);
            }
            let full = run.nodal_motion(cx).unwrap();
            assert_eq!(full[0][0], [0.0;3]);
            errors.push((full[0][1][0] - 0.001*5.0_f64.cos()).abs());
            assert_eq!(full[0][1][0], run.state().displacement[0]);
        }
        for ratio in [errors[0]/errors[1],errors[1]/errors[2]] { assert!((3.7..4.3).contains(&ratio), "{errors:?}"); }
        assert!(errors[2] < 2e-6);
    });
}

#[test]
fn g1_free_body_accelerates_without_artificial_rigid_mode_clamps() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0,60.0,0.25,[73;32]).unwrap();
        let problem = problem(&material,&[]); let geometry = geometry(cx);
        let acceleration = [0.2,-0.4,0.6];
        let load = geometry.assemble(&[],VolumeLoad3::Uniform(acceleration.map(|v| 6.0*v)),&[],cx).unwrap();
        let mut run = load.begin_transient(&problem,None,0.0,&ZERO,&ZERO,config(),cx).unwrap();
        for _ in 0..50 { run.advance(&load,0.01,cx).unwrap(); }
        let motion = run.nodal_motion(cx).unwrap(); let t = run.state().time_s;
        for node in 0..4 { for c in 0..3 {
            assert!((motion[0][node][c]-0.5*acceleration[c]*t*t).abs() < 1e-10);
            assert!((motion[1][node][c]-acceleration[c]*t).abs() < 1e-10);
            assert!((motion[2][node][c]-acceleration[c]).abs() < 1e-9);
        } }
        // Reference volume=1/6 and density=6 give total physical mass one kg.
        let expected = 0.5*acceleration.iter().map(|a| (a*t).powi(2)).sum::<f64>();
        assert!((run.state().kinetic_energy_j-expected).abs() < 1e-11);
    });
}

#[test]
fn g1_rayleigh_damping_matches_the_independent_damped_oscillator() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material,&fixed); let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let mut initial = ZERO; initial[1][0] = 0.001;
        // C=.4 M + .006 K=.1 kg/s, decay rate=C/(2M)=.5/s.
        let cfg = TransientConfig3 { mass_damping_per_s:0.4, stiffness_damping_s:0.006, ..config() };
        let mut run = load.begin_transient(&problem,None,0.0,&initial,&ZERO,cfg,cx).unwrap();
        let initial_energy = run.state().total_energy_j; let mut dissipated = 0.0;
        for _ in 0..500 {
            let previous = run.state().total_energy_j;
            let report = run.advance(&load,0.001,cx).unwrap();
            assert!(report.damping_loss_j >= 0.0);
            assert!(report.quadratic_energy_j <= previous+1e-15);
            dissipated += report.damping_loss_j;
        }
        let t = run.state().time_s; let w = 99.75_f64.sqrt();
        let expected = 0.001*(-0.5*t).exp()*((w*t).cos()+0.5/w*(w*t).sin());
        assert!((run.state().displacement[0]-expected).abs() < 5e-8);
        assert!((run.state().total_energy_j+dissipated-initial_energy).abs() < 1e-13);
    });
}

#[test]
fn g3_frozen_thermal_load_and_external_load_have_separate_work_accounts() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material,&fixed); let geometry = geometry(cx);
        let empty = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let force = geometry.assemble(&[],VolumeLoad3::None,&[(3,0.01)],cx).unwrap();
        let thermal = TetThermalStrainState::try_new(301.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            material.material_state_identity,fs_blake3::ContentHash([74;32]),0.01).unwrap();
        let mut hot = empty.begin_transient(&problem,Some(TetThermalStrainField::Uniform(&thermal)),
            0.0,&ZERO,&ZERO,config(),cx).unwrap();
        let mut pulled = force.begin_transient(&problem,None,0.0,&ZERO,&ZERO,config(),cx).unwrap();
        for _ in 0..50 {
            let a = hot.advance(&empty,0.002,cx).unwrap(); let b = pulled.advance(&force,0.002,cx).unwrap();
            assert!((hot.state().displacement[0]-pulled.state().displacement[0]).abs() < 1e-13);
            assert_eq!(a.mechanical_work_j,0.0); assert_eq!(b.thermal_equivalent_work_j,0.0);
            assert!((a.thermal_equivalent_work_j-b.mechanical_work_j).abs() < 1e-14);
        }
    });
}

#[test]
fn g1_thermal_eigenstrain_ramp_drives_the_same_inertial_motion_and_refreshed_stress() {
    with_cx(false, |cx| {
        let original = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let refreshed = TetElasticMaterial::try_new(6.0,60.0,0.0,[75;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let initial_problem = problem(&original,&fixed);
        let next_problem = problem(&refreshed,&fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let frame = fs_blake3::ContentHash([74;32]);
        let mut errors = Vec::new();
        for count in [50,100] {
            let mut run = load.begin_transient(&initial_problem,None,0.0,&ZERO,&ZERO,config(),cx).unwrap();
            let mut material_work = 0.0;
            for i in 1..=count {
                let time = i as f64/count as f64;
                let strain = TetThermalStrainState::try_new(300.0+time,300.0,
                    [0.001*time,0.0,0.0,0.0,0.0,0.0],refreshed.material_state_identity,frame,0.01).unwrap();
                let (next, report) = run.advance_thermal_candidate(&next_problem,
                    TetThermalStrainField::Uniform(&strain),&load,1.0/count as f64,frame,0.01,cx).unwrap();
                assert!(report.physical_balance_error_j.abs() < 1e-13);
                assert!(report.native.balance_error_j.abs() < 1e-13);
                assert_eq!(report.native.stiffness_parameter_work_j,0.0);
                assert_eq!(report.stiffness_material_work_j,0.0);
                material_work += report.thermal_material_work_j;
                run = next;
            }
            // m=.1 kg, k=10 N/m, free length grows at s=.001 m/s:
            // q(t)=s*(t-sin(10*t)/10), with zero initial q and velocity.
            let expected = 0.001*(1.0-10.0_f64.sin()/10.0);
            errors.push((run.state().displacement[0]-expected).abs());
            let response = run.recover_response(frame,0.01,cx).unwrap();
            let stress = 60.0*(run.state().displacement[0]-0.001);
            assert!((response.stress.elements()[0].stress_mandel_pa()[0]-stress).abs() < 1e-12);
            assert!((response.kinetic_plus_elastic_energy_j-material_work).abs() < 1e-13);
            let TetMaterialField::Uniform(material) = run.material_problem().materials else { panic!("uniform source expected") };
            assert_eq!(material.material_state_identity,refreshed.material_state_identity);
        }
        assert!((3.5..4.5).contains(&(errors[0]/errors[1])), "{errors:?}");
    });
}

#[test]
fn g1_fully_clamped_heating_and_cooling_account_for_elastic_storage_without_motion() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let problem = problem(&material,&fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let frame = fs_blake3::ContentHash([74;32]);
        let cold = load.begin_transient(&problem,None,0.0,&ZERO,&ZERO,config(),cx).unwrap();
        let hot_strain = TetThermalStrainState::try_new(310.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            material.material_state_identity,frame,0.01).unwrap();
        let (hot, report) = cold.advance_thermal_candidate(&problem,
            TetThermalStrainField::Uniform(&hot_strain),&load,0.01,frame,0.01,cx).unwrap();
        let stored = 0.5*(1.0/6.0)*60.0*0.001_f64.powi(2);
        assert_eq!(report.native.thermal_equivalent_work_j,0.0);
        assert_eq!(report.native.quadratic_energy_j,0.0);
        assert!((report.thermal_material_work_j-stored).abs() < 1e-15);
        assert!((report.kinetic_plus_elastic_after_j-stored).abs() < 1e-15);
        assert!(report.physical_balance_error_j.abs() < 1e-15);
        assert_eq!(hot.nodal_motion(cx).unwrap(),[ZERO.to_vec(),ZERO.to_vec(),ZERO.to_vec()]);
        let response = hot.recover_response(frame,0.01,cx).unwrap();
        assert!((response.stress.elements()[0].stress_mandel_pa()[0]+0.06).abs() < 1e-14);
        let zero = TetThermalStrainState::try_new(300.0,300.0,[0.0;6],
            material.material_state_identity,frame,0.01).unwrap();
        let (cooled, cooling) = hot.advance_thermal_candidate(&problem,
            TetThermalStrainField::Uniform(&zero),&load,0.01,frame,0.01,cx).unwrap();
        assert!((cooling.thermal_material_work_j+stored).abs() < 1e-15);
        assert_eq!(cooling.kinetic_plus_elastic_after_j,0.0);
        assert!(cooled.recover_response(frame,0.01,cx).unwrap().stress.elements()[0]
            .stress_mandel_pa().iter().all(|s| *s == 0.0));
        assert_eq!(cold.state().time_s,0.0);
    });
}

#[test]
fn g1_temperature_dependent_stiffness_uses_endpoint_operator_and_signed_material_work() {
    with_cx(false, |cx| {
        // One free axial DOF: m=.1, k0=10, k1=6. The independently solved
        // Newmark scalar reference includes changing eigenstrain and force.
        let initial_material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let heated_material = TetElasticMaterial::try_new(6.0,36.0,0.0,[75;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let initial_problem = problem(&initial_material,&fixed);
        let heated_problem = problem(&heated_material,&fixed);
        let geometry = geometry(cx);
        let old_load = geometry.assemble(&[],VolumeLoad3::None,&[(3,0.0002)],cx).unwrap();
        let next_load = geometry.assemble(&[],VolumeLoad3::None,&[(3,-0.0001)],cx).unwrap();
        let frame = fs_blake3::ContentHash([74;32]);
        let old_strain = TetThermalStrainState::try_new(300.0,290.0,[0.0003,0.0,0.0,0.0,0.0,0.0],
            initial_material.material_state_identity,frame,0.01).unwrap();
        let next_strain = TetThermalStrainState::try_new(301.0,290.0,[0.0012,0.0,0.0,0.0,0.0,0.0],
            heated_material.material_state_identity,frame,0.01).unwrap();
        let mut u = ZERO; u[1][0] = 0.002;
        let mut v = ZERO; v[1][0] = -0.001;
        let cfg = TransientConfig3 { mass_damping_per_s:0.4, ..config() };
        let run = old_load.begin_transient(&initial_problem,Some(TetThermalStrainField::Uniform(&old_strain)),
            0.0,&u,&v,cfg,cx).unwrap();
        let saved = run.state().clone();
        let (heated, report) = run.advance_thermal_candidate(&heated_problem,
            TetThermalStrainField::Uniform(&next_strain),&next_load,0.02,frame,0.01,cx).unwrap();
        // Exact rational scalar solution, independently evaluated.
        let q1 = 0.001958712871287129;
        let v1 = -0.003128712871287129;
        assert!((heated.state().displacement[0]-q1).abs() < 1e-14);
        assert!((heated.state().velocity[0]-v1).abs() < 1e-13);
        assert!((report.native.stiffness_parameter_work_j+7.834851485148515e-6).abs() < 1e-14);
        assert!((report.stiffness_material_work_j+2.5796237623762375e-6).abs() < 1e-14);
        assert!((report.native.damping_loss_j-3.4092539947064014e-9).abs() < 1e-15);
        assert!(report.native.balance_error_j.abs() < 1e-13);
        assert!(report.physical_balance_error_j.abs() < 1e-13);
        let response = heated.recover_response(frame,0.01,cx).unwrap();
        assert!((response.stress.elements()[0].stress_mandel_pa()[0]-36.0*(q1-0.0012)).abs() < 1e-12);
        assert_eq!(run.state(),&saved);

        // Continue from K1 and retained q1/v1/a1 while the source stiffens.
        let final_strain = TetThermalStrainState::try_new(302.0,290.0,[0.0006,0.0,0.0,0.0,0.0,0.0],
            initial_material.material_state_identity,frame,0.01).unwrap();
        let (stiffened, second) = heated.advance_thermal_candidate(&initial_problem,
            TetThermalStrainField::Uniform(&final_strain),&next_load,0.02,frame,0.01,cx).unwrap();
        let qp = q1+0.02*v1+0.0001*heated.state().acceleration[0];
        let vp = v1+0.01*heated.state().acceleration[0];
        let a2 = (-0.0001+10.0*0.0006-10.0*qp-0.04*vp)/(0.1+0.0004+0.001);
        let q2 = qp+0.0001*a2;
        assert!((stiffened.state().displacement[0]-q2).abs() < 1e-14);
        assert!((second.native.stiffness_parameter_work_j-0.5*q1*4.0*q2).abs() < 1e-14);
        assert!(second.physical_balance_error_j.abs() < 1e-13);
        assert_eq!(stiffened.state().step,2);
    });
}

#[test]
fn g1_clamped_softening_and_expansion_keep_distinct_physical_work_accounts() {
    with_cx(false, |cx| {
        let cold_material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let hot_material = TetElasticMaterial::try_new(6.0,30.0,0.0,[75;32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let cold_problem = problem(&cold_material,&fixed);
        let hot_problem = problem(&hot_material,&fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let frame = fs_blake3::ContentHash([74;32]);
        let cold_strain = TetThermalStrainState::try_new(310.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            cold_material.material_state_identity,frame,0.01).unwrap();
        let hot_strain = TetThermalStrainState::try_new(320.0,300.0,[0.002,0.0,0.0,0.0,0.0,0.0],
            hot_material.material_state_identity,frame,0.01).unwrap();
        let cold = load.begin_transient(&cold_problem,Some(TetThermalStrainField::Uniform(&cold_strain)),
            0.0,&ZERO,&ZERO,config(),cx).unwrap();
        let (hot, report) = cold.advance_thermal_candidate(&hot_problem,
            TetThermalStrainField::Uniform(&hot_strain),&load,0.01,frame,0.01,cx).unwrap();
        assert_eq!(report.native.thermal_equivalent_work_j,0.0);
        assert_eq!(report.native.stiffness_parameter_work_j,0.0);
        assert!((report.kinetic_plus_elastic_before_j-5e-6).abs() < 1e-15);
        assert!((report.kinetic_plus_elastic_after_j-1e-5).abs() < 1e-15);
        assert!((report.thermal_material_work_j-1e-5).abs() < 1e-15);
        assert!((report.stiffness_material_work_j+5e-6).abs() < 1e-15);
        assert!(report.physical_balance_error_j.abs() < 1e-15);
        assert_eq!(hot.nodal_motion(cx).unwrap(),[ZERO.to_vec(),ZERO.to_vec(),ZERO.to_vec()]);
    });
}

#[test]
fn g4_thermal_refresh_refuses_changed_mass_and_late_strain_without_mutating_history() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let changed = TetElasticMaterial::try_new(7.0,60.0,0.0,[75;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let original = problem(&material,&fixed); let wrong = problem(&changed,&fixed);
        let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let run = load.begin_transient(&original,None,0.0,&ZERO,&ZERO,config(),cx).unwrap();
        let saved = run.state().clone();
        let frame = fs_blake3::ContentHash([74;32]);
        let strain = TetThermalStrainState::try_new(310.0,300.0,[0.01,0.0,0.0,0.0,0.0,0.0],
            material.material_state_identity,frame,0.1).unwrap();
        assert!(matches!(run.advance_thermal_candidate(&wrong,TetThermalStrainField::Uniform(&strain),
            &load,0.01,frame,0.1,cx),Err(TransientError3::Invalid(_))));
        let softened = TetElasticMaterial::try_new(6.0,30.0,0.0,[75;32]).unwrap();
        let softened_problem = problem(&softened,&fixed);
        let softened_strain = TetThermalStrainState::try_new(310.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            softened.material_state_identity,frame,0.1).unwrap();
        let damped = load.begin_transient(&original,None,0.0,&ZERO,&ZERO,
            TransientConfig3 { stiffness_damping_s:0.01, ..config() },cx).unwrap();
        let (damped_next, damped_report) = damped.advance_thermal_candidate(&softened_problem,
            TetThermalStrainField::Uniform(&softened_strain),&load,0.01,frame,0.1,cx).unwrap();
        assert!(damped_report.native.damping_loss_j >= 0.0);
        assert!(damped_report.physical_balance_error_j.abs() < 1e-12);
        assert_eq!(damped_next.state().time_s,0.01);
        assert_eq!(damped.state(),&saved);
        // The native inertial candidate is computed, then physical stress
        // recovery refuses an over-limit elastic strain. Accepted history stays.
        assert!(run.advance_thermal_candidate(&original,TetThermalStrainField::Uniform(&strain),
            &load,0.01,frame,0.001,cx).is_err());
        with_cx(true, |cancelled| {
            assert!(matches!(run.advance_thermal_candidate(&original,TetThermalStrainField::Uniform(&strain),
                &load,0.01,frame,0.1,cancelled),Err(TransientError3::Cancelled)));
        });
        assert_eq!(run.state(),&saved);
        assert!(run.thermal_strain_field().is_none());
    });
}

#[test]
fn g4_checkpoint_resume_limits_and_cancellation_never_publish_a_partial_step() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0,60.0,0.0,[73;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let problem = problem(&material,&fixed); let geometry = geometry(cx);
        let load = geometry.assemble(&[],VolumeLoad3::None,&[(3,0.01)],cx).unwrap();
        let cfg = TransientConfig3 { maximum_steps: 3, ..config() };
        let mut run = load.begin_transient(&problem,None,0.0,&ZERO,&ZERO,cfg,cx).unwrap();
        run.advance(&load,0.01,cx).unwrap(); let mut fork = run.clone();
        let before = run.state().clone();
        with_cx(true, |cancelled| {
            assert_eq!(run.advance(&load,0.01,cancelled),Err(TransientError3::Cancelled));
        });
        assert_eq!(run.state(),&before);
        for dt in [0.0,-0.01,f64::NAN,f64::INFINITY,1e-300] {
            assert!(run.advance(&load,dt,cx).is_err()); assert_eq!(run.state(),&before);
        }
        for _ in 0..2 {
            assert_eq!(run.advance(&load,0.01,cx).unwrap(),fork.advance(&load,0.01,cx).unwrap());
            assert_eq!(run.state(),fork.state());
        }
        let before = run.state().clone();
        assert_eq!(run.advance(&load,0.01,cx),Err(TransientError3::Budget));
        assert_eq!(run.state(),&before);
        let mut capped = load.begin_transient(&problem,None,0.0,&ZERO,&ZERO,
            TransientConfig3 { maximum_displacement_m:1e-6,..config() },cx).unwrap();
        let before = capped.state().clone();
        assert!(capped.advance(&load,0.1,cx).is_err()); assert_eq!(capped.state(),&before);
    });
}

#[test]
fn g0_missing_motion_bad_damping_foreign_loads_and_iteration_exhaustion_refuse() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0,60.0,0.25,[73;32]).unwrap();
        let fixed = [0,1,2,4,5,8]; let problem = problem(&material,&fixed); let geometry = geometry(cx);
        let empty = geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        assert!(empty.begin_transient(&problem,None,0.0,&[],&ZERO,config(),cx).is_err());
        let mut bad = ZERO; bad[0][0]=0.001;
        assert!(empty.begin_transient(&problem,None,0.0,&bad,&ZERO,config(),cx).is_err());
        for cfg in [TransientConfig3 { mass_damping_per_s:-1.0,..config() },
            TransientConfig3 { maximum_free_dofs:1,..config() },
            TransientConfig3 { maximum_element_entries:143,..config() }] {
            assert!(empty.begin_transient(&problem,None,0.0,&ZERO,&ZERO,cfg,cx).is_err());
        }
        let mut run = empty.begin_transient(&problem,None,0.0,&ZERO,&ZERO,
            TransientConfig3 { solve:TetStaticSolveConfig { maximum_iterations:1,..Default::default() },..config() },cx).unwrap();
        let load = geometry.assemble(&[],VolumeLoad3::None,&[(3,0.01),(6,-0.02),(10,0.03)],cx).unwrap();
        let before = run.state().clone();
        assert!(matches!(run.advance(&load,0.1,cx),Err(TransientError3::Solve { .. })));
        assert_eq!(run.state(),&before);
        let shifted = NODES.map(|p| p.map(|x| x+2.0));
        let other = TetLoadGeometry3::new(&shifted,&TETS,TetAssemblyBudget::default(),4,cx).unwrap();
        let foreign = other.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        assert!(run.advance(&foreign,0.01,cx).is_err()); assert_eq!(run.state(),&before);
    });
}
