use super::*;

struct Diagonal(Vec<f64>);
impl LinearOp for Diagonal {
    fn n(&self) -> usize { self.0.len() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for ((out, &k), &u) in y.iter_mut().zip(&self.0).zip(x) { *out = k*u; }
    }
}
fn law() -> BilinearCohesive { BilinearCohesive::new(1e9,1e6,2000.0,2e9,0.02).unwrap() }
fn patch<'a>(id: u64, area: f64, law: &'a BilinearCohesive) -> CohesivePatch<'a> {
    CohesivePatch::new(id,area,law,vec![JumpTerm { dof: 0, local: [1.0,0.0,0.0] }]).unwrap()
}
fn budget() -> Budget { Budget { max_dofs: 32, max_patches: 8, max_jump_terms: 64, max_krylov_basis_entries: 50_000 } }
fn config() -> Config {
    Config { newton: NewtonKrylovConfig { absolute_tolerance: 1e-9, relative_tolerance: 1e-13,
        linear_restart: 8, forcing_maximum: 0.05, ..NewtonKrylovConfig::default() },
        max_newton_iterations: 80, max_operator_evaluations: 2000,
        force_tolerance_n: 1e-7, local_energy_tolerance_j_m2: 1e-7,
        total_cohesive_energy_tolerance_j: 1e-9, maximum_displacement_m: 0.02 }
}
fn close(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() < tolerance, "{a:e} != {b:e}"); }

#[test]
fn compliant_pull_fixture_matches_closed_form_softening_and_releases_load() {
    let law = law();
    let stiffness = 5e6;
    let area = 1e-4;
    let elastic = Diagonal(vec![stiffness]);
    let system = CohesiveSystem::new(&elastic,vec![patch(1,area,&law)],budget()).unwrap();
    let mut state = system.bonded_state();
    let mut work = 0.0;
    for actuator in [0.0004,0.0012,0.0022,0.0032,0.0045] {
        let result = state.solve_controlled(&[stiffness*actuator],&[[0.0;3]],config(),&mut || false).unwrap();
        let elastic_jump = stiffness*actuator/(stiffness+area*1e9);
        let expected = if elastic_jump <= law.onset_separation_m() { elastic_jump }
            else if actuator >= law.failure_separation_m() { actuator }
            else { (stiffness*actuator-area*law.softening_stiffness_pa_m()*law.failure_separation_m())
                /(stiffness-area*law.softening_stiffness_pa_m()) };
        close(result.state.displacements_m()[0],expected,1e-12);
        let reaction = stiffness*(actuator-result.state.displacements_m()[0]);
        close(reaction,area*result.state.points()[0].traction_pa()[0],1e-7);
        assert!(result.state.points()[0].damage() >= state.points()[0].damage());
        assert!(result.maximum_force_residual_n < 1e-7);
        work += result.cohesive_energy_j.work;
        state = result.state;
    }
    assert!(state.points()[0].failed());
    close(work,area*law.fracture_energy_j_m2(),1e-9);
    close(state.points()[0].traction_pa()[0],0.0,1e-9);
    let reclosed = state.solve_controlled(&[-500.0],&[[0.0;3]],config(),&mut || false).unwrap();
    assert!(reclosed.state.points()[0].failed());
    assert!(reclosed.state.points()[0].traction_pa()[0] < 0.0);
    close(reclosed.cohesive_energy_j.dissipated,0.0,1e-12);
}

#[test]
fn prescribed_jump_moves_the_physical_interface_and_unloads_without_healing() {
    let law = law();
    let elastic = Diagonal(vec![2e6]);
    let system = CohesiveSystem::new(&elastic,vec![patch(10,1e-4,&law)],budget()).unwrap();
    let initial = system.bonded_state();
    let result = initial.solve_controlled(&[0.0],&[[0.002,0.0,0.0]],config(),&mut || false).unwrap();
    assert!(result.state.displacements_m()[0] < 0.0);
    close(result.state.points()[0].jump_m()[0],result.state.displacements_m()[0]+0.002,1e-14);
    assert!(result.state.points()[0].damage() > 0.0);
    let saved_damage = result.state.points()[0].damage();
    let unload = result.state.solve_controlled(&[0.0],&[[0.0005,0.0,0.0]],config(),&mut || false).unwrap();
    close(unload.state.points()[0].damage(),saved_damage,1e-14);
    close(unload.cohesive_energy_j.dissipated,0.0,1e-12);
}

#[test]
fn equal_area_partition_preserves_response_and_fracture_work() {
    let law = law();
    let elastic = Diagonal(vec![5e6]);
    let whole = CohesiveSystem::new(&elastic,vec![patch(1,1e-4,&law)],budget()).unwrap();
    let split = CohesiveSystem::new(&elastic,vec![patch(1,4e-5,&law),patch(2,6e-5,&law)],budget()).unwrap();
    let mut a = whole.bonded_state(); let mut b = split.bonded_state();
    for force in [5000.0,10000.0,3000.0,18000.0,22500.0] {
        let aa = a.solve_controlled(&[force],&[[0.0;3]],config(),&mut || false).unwrap();
        let bb = b.solve_controlled(&[force],&[[0.0;3];2],config(),&mut || false).unwrap();
        close(aa.state.displacements_m()[0],bb.state.displacements_m()[0],1e-12);
        close(aa.cohesive_energy_j.work,bb.cohesive_energy_j.work,1e-10);
        close(aa.cohesive_energy_j.dissipated,bb.cohesive_energy_j.dissipated,1e-10);
        a = aa.state; b = bb.state;
    }
}

#[test]
fn sparse_jump_scatter_is_work_conjugate_and_structural_tangent_is_consistent() {
    let law = law();
    let terms = vec![
        JumpTerm { dof: 0, local: [-0.6,0.8,0.0] },
        JumpTerm { dof: 1, local: [-0.8,-0.6,0.0] },
        JumpTerm { dof: 2, local: [0.0,0.0,-1.0] },
        JumpTerm { dof: 3, local: [0.6,-0.8,0.0] },
        JumpTerm { dof: 4, local: [0.8,0.6,0.0] },
        JumpTerm { dof: 5, local: [0.0,0.0,1.0] },
    ];
    let patch = CohesivePatch::new(7,1e-4,&law,terms).unwrap();
    let direction = [0.3,-0.4,0.2,0.1,0.5,-0.1];
    let traction = [3e5,-2e4,8e4];
    let mut forces = [0.0;6]; patch.scatter(traction,&mut forces).unwrap();
    close(direction.iter().zip(forces).map(|(u,f)| u*f).sum(),
        patch.area_m2()*dot(traction,patch.jump(&direction,[0.0;3]).unwrap()),1e-12);
    for i in 0..3 { close(forces[i]+forces[i+3],0.0,1e-12); }
    let elastic = Diagonal(vec![2e6;6]);
    let system = CohesiveSystem::new(&elastic,vec![patch],budget()).unwrap();
    let state = system.bonded_state();
    let mut cancel = || false;
    let problem = Equilibrium { state: &state, force: &[0.0;6], offsets: &[[0.0;3]], config: config(),
        cancelled: RefCell::new(&mut cancel), fault: RefCell::new(None), evaluations: Cell::new(0) };
    let x = [0.0,0.0,0.0,0.002,0.001,0.0004];
    let mut jac = [0.0;6]; problem.action(&x,Some(&direction),&mut jac).unwrap();
    let h = 1e-8;
    let plus = std::array::from_fn::<_,6,_>(|i| x[i]+h*direction[i]);
    let minus = std::array::from_fn::<_,6,_>(|i| x[i]-h*direction[i]);
    let mut a = [0.0;6]; let mut b = [0.0;6];
    problem.action(&plus,None,&mut a).unwrap(); problem.action(&minus,None,&mut b).unwrap();
    for i in 0..6 { close((a[i]-b[i])/(2.0*h),jac[i],0.02); }
    let mut load = [0.0;6]; problem.action(&x,None,&mut load).unwrap();
    let solved = state.solve_controlled(&load,&[[0.0;3]],config(),&mut || false).unwrap();
    for (actual, expected) in solved.state.displacements_m().iter().zip(x) { close(*actual,expected,1e-11); }
}

#[test]
fn failed_cancelled_and_budgeted_trials_never_commit_damage_and_replay_is_identical() {
    let law = law();
    let elastic = Diagonal(vec![5e6]);
    let system = CohesiveSystem::new(&elastic,vec![patch(1,1e-4,&law)],budget()).unwrap();
    let state = system.bonded_state().solve_controlled(&[10000.0],&[[0.0;3]],config(),&mut || false).unwrap().state;
    let saved = state.clone();
    let mut count = 0;
    let cancelled = state.solve_controlled(&[18000.0],&[[0.0;3]],config(),&mut || { count += 1; count >= 12 });
    assert!(matches!(cancelled,Err(Error::Cancelled)));
    let mut limited = config(); limited.max_operator_evaluations = 1;
    assert!(matches!(state.solve_controlled(&[18000.0],&[[0.0;3]],limited,&mut || false),Err(Error::Budget(_))));
    assert!(state.solve_controlled(&[18000.0],&[[0.03,0.0,0.0]],config(),&mut || false).is_err());
    assert_eq!(state.points(),saved.points()); assert_eq!(state.displacements_m(),saved.displacements_m());
    let a = state.solve_controlled(&[18000.0],&[[0.0;3]],config(),&mut || false).unwrap();
    let b = saved.solve_controlled(&[18000.0],&[[0.0;3]],config(),&mut || false).unwrap();
    assert_eq!(a.state.points(),b.state.points()); assert_eq!(a.state.displacements_m(),b.state.displacements_m());
    assert_eq!(a.cohesive_energy_j,b.cohesive_energy_j);
}

#[test]
fn invalid_maps_dimensions_inputs_and_krylov_capacity_refuse() {
    let law = law(); let elastic = Diagonal(vec![5e6]);
    assert!(CohesivePatch::new(1,0.0,&law,vec![JumpTerm { dof: 0, local: [1.0,0.0,0.0] }]).is_err());
    assert!(CohesivePatch::new(1,1e-4,&law,vec![JumpTerm { dof: 0, local: [1.0,0.0,0.0] };2]).is_err());
    assert!(CohesiveSystem::new(&elastic,vec![patch(1,1e-4,&law),patch(1,1e-4,&law)],budget()).is_err());
    let system = CohesiveSystem::new(&elastic,vec![patch(1,1e-4,&law)],budget()).unwrap();
    let state = system.bonded_state();
    assert!(state.solve_controlled(&[],&[[0.0;3]],config(),&mut || false).is_err());
    assert!(state.solve_controlled(&[f64::NAN],&[[0.0;3]],config(),&mut || false).is_err());
    let mut cfg = config(); cfg.newton.linear_restart = 1000;
    assert!(matches!(state.solve_controlled(&[1.0],&[[0.0;3]],cfg,&mut || false),Err(Error::Budget(_))));
    let mut cfg = config(); cfg.maximum_displacement_m = 1e-6;
    assert!(state.solve_controlled(&[10000.0],&[[0.0;3]],cfg,&mut || false).is_err());
}

#[test]
fn series_planes_distinguish_interface_from_declared_feature_failure() {
    let weak = BilinearCohesive::new(1e9,2e5,100.0,2e9,0.01).unwrap();
    let strong = law();
    let elastic = Diagonal(vec![0.0,1e6]);
    for (bond,feature,failed) in [(&weak,&strong,0),(&strong,&weak,1)] {
        let system = CohesiveSystem::new(&elastic,vec![patch(1,1e-4,bond),
            CohesivePatch::new(2,1e-4,feature,vec![
                JumpTerm { dof: 0, local: [-1.0,0.0,0.0] },
                JumpTerm { dof: 1, local: [1.0,0.0,0.0] },
            ]).unwrap()],budget()).unwrap();
        let mut state = system.bonded_state();
        let mut work = 0.0;
        for step in 0..=60 {
            let actuator = f64::from(step)*2e-5;
            let result = state.solve_controlled(&[0.0,1e6*actuator],&[[0.0;3];2],config(),&mut || false).unwrap();
            work += result.cohesive_energy_j.work;
            state = result.state;
            let reaction = 1e6*(actuator-state.displacements_m()[1]);
            for point in state.points() { close(1e-4*point.traction_pa()[0],reaction,1e-7); }
        }
        assert!(state.points()[failed].failed());
        assert_eq!(state.points()[1-failed].damage(),0.0);
        close(work,0.01,1e-9);
    }
}


#[test]
fn anisotropic_mixed_mode_law_uses_native_equilibrium_and_retains_onset_history() {
    use fs_material::cohesive::mixed_mode::MixedModeParameters;
    let law = BilinearCohesive::new_mixed_mode(
        MixedModeParameters {
            normal_stiffness_pa_m: 1e9,
            tangential_stiffness_pa_m: 4e8,
            normal_strength_pa: 2e5,
            tangential_strength_pa: 1.6e5,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 400.0,
            bk_exponent: 2.0,
        },
        2e9,
        0.02,
    )
    .unwrap();
    let elastic = Diagonal(vec![2e9, 2e9]);
    let patch = CohesivePatch::new(
        90,
        1.0,
        &law,
        vec![
            JumpTerm {
                dof: 0,
                local: [1.0, 0.0, 0.0],
            },
            JumpTerm {
                dof: 1,
                local: [0.0, 1.0, 0.0],
            },
        ],
    )
    .unwrap();
    let system = CohesiveSystem::new(&elastic, vec![patch], budget()).unwrap();
    let mut state = system.bonded_state();
    let mut expected = law.bonded_point();
    for target in [
        [0.00005, 0.00003, 0.0],
        [0.0004, 0.0007, 0.0],
        [-0.0002, 0.0015, 0.0],
        [0.0001, 0.0002, 0.0],
        [0.0, 0.012, 0.0],
        [-0.001, 0.001, 0.0],
    ] {
        let prediction = expected.trial(target, 1e-7).unwrap();
        let force = [
            2e9 * target[0] + prediction.traction_pa[0],
            2e9 * target[1] + prediction.traction_pa[1],
        ];
        let mut controls = config();
        controls.total_cohesive_energy_tolerance_j = 1e-7;
        let next = state
            .solve_controlled(&force, &[[0.0; 3]], controls, &mut || false)
            .unwrap();
        for (actual, wanted) in next.state.displacements_m().iter().zip(target) {
            close(*actual, wanted, 2e-12);
        }
        close(
            next.cohesive_energy_j.dissipated,
            prediction.energy.dissipated,
            1e-6,
        );
        if let Some(old) = state.points()[0].mixed_mode_history() {
            assert_eq!(next.state.points()[0].mixed_mode_history(), Some(old));
        }
        assert!(next.state.points()[0].damage() >= state.points()[0].damage());
        expected = prediction.point;
        state = next.state;
    }
    assert!(state.points()[0].failed());
    close(state.points()[0].traction_pa()[0], -2e6, 1e-4);
    assert_eq!(state.points()[0].traction_pa()[1], 0.0);
}
