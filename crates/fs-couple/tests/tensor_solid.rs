//! Synthetic solid/material oracles, not measured specimen validation.
use fs_couple::tensor_material::{
    InitialMemory, TensorMaterialBody, TensorMaterialBudget,
    solid::TensorSolidBody,
};
use fs_exec::{CancelGate, Cx};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::visco_tensor::TensorProjectionBudget;
use fs_solid::linear3::{
    TetAssemblyBudget,
    tet10::{TetQuadrature, promote_p1_to_p2},
    visco::{SolidMeshRef, SolidProjectionBudget, TetViscoField, TetViscoMaterial,
        ViscoElements, ViscoStateBudget},
};

const NODES: [[f64;3];4] = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const FIXED: [usize;11] = [0,1,2,3,4,5,6,7,8,9,10];
fn with_cx<R>(f:impl FnOnce(&Cx<'_>,&CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| f(&Cx::new(&gate,arena,
        fs_exec::StreamKey {seed:83,kernel_id:83,tile:0,iteration:0},
        fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn diagonal(values:[f64;6])->Mandel6 {
    core::array::from_fn(|i|core::array::from_fn(|j|if i==j {values[i]} else {0.0}))
}
fn material(viscous_time:f64,branch_fraction:f64)->TetViscoMaterial {
    let equilibrium=diagonal([120.0,90.0,120.0,60.0,50.0,40.0]);
    let branches=if branch_fraction==0.0 {vec![]} else {
        vec![PronyBranch::new(equilibrium.map(|r|r.map(|x|x*branch_fraction)),0.08).unwrap()]
    };
    TetViscoMaterial {
        law:TensorProny::new(equilibrium,
            Some(equilibrium.map(|r|r.map(|x|x*viscous_time))),branches,(0.0,1e5)).unwrap(),
        density_kg_m3:2.0,material_state_identity:[83;32],
    }
}
fn projection_budget()->SolidProjectionBudget {
    SolidProjectionBudget {max_coordinates:8,max_basis_scalars:4096,
        max_material_points:64,max_scalar_terms:1_000_000}
}
fn dynamics_budget()->TensorMaterialBudget {
    TensorMaterialBudget {max_state_coordinates:112,max_material_points:64}
}
fn system(mesh:SolidMeshRef<'_>,materials:TetViscoField<'_>,fixed:&[usize],cx:&Cx<'_>)->ViscoElements {
    ViscoElements::new(mesh,materials,fixed,TetAssemblyBudget::standard(),
        ViscoStateBudget {max_internal_variables:10_000},cx).unwrap()
}
fn close(a:f64,b:f64,tol:f64) {
    assert!((a-b).abs()<=tol*(1.0+b.abs()),"{a:e} != {b:e}");
}
fn dot(a:&[f64],b:&[f64])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}

#[test]
fn full_mass_drives_work_conjugate_force_and_inertia_without_lumping() {
    let gate=CancelGate::new_clock_free();
    let law=TensorProny::new_mandel(diagonal([1.0;6]),vec![],0.1).unwrap();
    // Two strain-free translation coordinates isolate consistent inertia.
    let point=law.project_strain(&[[0.0;6];2],1.0,
        TensorProjectionBudget {max_coordinates:2,max_scalar_terms:1000},||false).unwrap();
    let body=TensorMaterialBody::from_mass_matrix(2,vec![2.0,1.0,1.0,3.0],
        vec![point],dynamics_budget(),&gate).unwrap();
    assert_eq!(body.mass_matrix_kg(),vec![2.0,1.0,1.0,3.0]);
    let mut h=body.history(&[0.0,0.0],&[0.2,-0.1],InitialMemory::Relaxed,&gate).unwrap();
    close(h.state()[2],0.3,1e-13);close(h.state()[3],-0.1,1e-13);
    close(h.stored_energy_j(),0.035,1e-13);
    let step=h.advance(&[3.0,-1.0],0.4,&gate).unwrap();
    // Independent M^-1 f = [2,-1]; q=v0*t+a*t^2/2, p=p0+f*t.
    for (actual,expected) in h.state().iter().zip([0.24,-0.12,1.5,-0.5]) {
        close(*actual,expected,2e-12);
    }
    let velocity=h.velocity_m_s();close(velocity[0],1.0,2e-12);close(velocity[1],-0.5,2e-12);
    close(step.y[0],0.6,2e-12);close(step.y[1],-0.3,2e-12);
    close(step.supplied,0.84,2e-12);close(step.dissipated,0.0,1e-14);
    close(h.balance_residual_j(),0.0,2e-12);
}

#[test]
fn original_tet4_and_tet10_quadrature_preserves_shear_mass_and_virtual_work() {
    with_cx(|cx,gate| {
        let material=material(0.01,0.25);
        let quadrature=TetQuadrature {stiffness_points_per_dir:2,mass_points_per_dir:4};
        let (qnodes,qtets)=promote_p1_to_p2(&NODES,&TETS,None,cx).unwrap();
        let cases=[
            (NODES.as_slice(),SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS}),
            (qnodes.as_slice(),SolidMeshRef::Tet10 {nodes_m:&qnodes,elements:&qtets,quadrature}),
        ];
        let gradients=[
            [[1.0,0.0,0.0],[0.0,0.3,0.0],[0.0,0.0,0.2]],
            [[0.2,1.0,0.0],[0.0,0.0,0.5],[0.2,0.0,0.0]],
        ];
        let strains=[
            [1.0,0.3,0.2,0.0,0.0,0.0],
            [0.2,0.0,0.0,0.5,0.25,0.1],
        ];
        for (nodes,mesh) in cases {
            let source=system(mesh,TetViscoField::Uniform(&material),&[],cx);
            let basis:Vec<Vec<f64>>=gradients.iter().map(|gradient|
                nodes.iter().flat_map(|p| gradient.map(|row|dot(&row,p))).collect()).collect();
            let projected=source.project_basis(&basis,projection_budget(),cx).unwrap();
            let mut stiffness=vec![0.0;4];let mut viscosity=vec![0.0;4];let mut branch=vec![0.0;4];
            for point in projected.material_points() {
                for i in 0..4 {stiffness[i]+=point.equilibrium_stiffness()[i];viscosity[i]+=point.viscosity_matrix()[i];}
                for b in point.branches() {for row in b.rows() {
                    for i in 0..2 {for j in 0..2 {branch[2*i+j]+=row[i]*row[j];}}
                }}
            }
            let volume=1.0/6.0;
            let weights=[1.0,1.0,1.0,2.0_f64.sqrt(),2.0_f64.sqrt(),2.0_f64.sqrt()];
            for i in 0..2 {for j in 0..2 {
                // Analytic affine-strain energy; this catches Mandel/tensor
                // shear confusion independently of either element projection.
                let expected=volume*(0..6).map(|k|
                    strains[i][k]*weights[k]*material.law.c_eq()[k][k]*strains[j][k]*weights[k]).sum::<f64>();
                close(stiffness[2*i+j],expected,2e-11);
                close(viscosity[2*i+j],0.01*expected,2e-11);
                close(branch[2*i+j],0.25*expected,2e-11);
                // Unit-tetrahedron integrals: integral(x_i*x_i)=V/10,
                // integral(x_i*x_j)=V/20. All inertia cross terms survive.
                let mut expected_mass=0.0;
                for component in 0..3 {for a in 0..3 {for b in 0..3 {
                    expected_mass+=material.density_kg_m3*volume*gradients[i][component][a]
                        *gradients[j][component][b]*if a==b {0.1} else {0.05};
                }}}
                close(projected.mass_matrix_kg()[2*i+j],expected_mass,2e-12);
            }}
            assert!(projected.mass_matrix_kg()[1]>0.0);
            let force:Vec<_>=(0..nodes.len()).map(|i|[0.1+0.02*i as f64,-0.2,0.3]).collect();
            let generalized=projected.project_nodal_force(&force,gate).unwrap();
            let velocity=projected.expand_nodal(&[0.05,-0.03],gate).unwrap();
            let nodal_power=force.iter().zip(&velocity).map(|(f,v)|dot(f,v)).sum::<f64>();
            close(nodal_power,dot(&generalized,&[0.05,-0.03]),1e-13);
        }
    })
}

#[test]
fn actual_solid_motion_retains_maxwell_prehistory_and_physical_heat() {
    with_cx(|cx,gate| {
        let material=material(0.003,0.5);
        let source=system(SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS},
            TetViscoField::Uniform(&material),&FIXED,cx);
        let projection=source.project_basis(&[vec![1.0]],projection_budget(),cx).unwrap();
        let body=TensorSolidBody::new(projection,dynamics_budget(),gate).unwrap();
        assert_eq!(body.body().state_dimension(),8); // q,p and all six branch rows
        let mut relaxed=body.history(&[0.001],&[0.0],InitialMemory::Relaxed,gate).unwrap();
        let mut loaded=body.history(&[0.001],&[0.0],InitialMemory::Instantaneous,gate).unwrap();
        let initial=loaded.material_history().stored_energy_j();
        close(initial,0.5*(120.0+60.0)/6.0*0.001_f64.powi(2),1e-13);
        let force=[[0.0;3];4];
        relaxed.advance(&force,0.001,1e-11,gate).unwrap();
        loaded.advance(&force,0.001,1e-11,gate).unwrap();
        assert!(loaded.nodal_velocity_m_s(gate).unwrap()[3][2]
            < relaxed.nodal_velocity_m_s(gate).unwrap()[3][2]);
        let mut heat=loaded.material_history().dissipated_heat_j();
        for _ in 0..40 {
            let step=loaded.advance(&force,0.025,1e-11,gate).unwrap();
            assert!(step.element_heat_j[0]>=0.0);
            close(step.element_heat_j[0],step.dynamics.dissipated,1e-12);
            heat+=step.element_heat_j[0];
        }
        assert!(loaded.material_history().stored_energy_j()<0.5*initial);
        close(loaded.material_history().stored_energy_j()+heat,initial,2e-11);
        close(loaded.material_history().balance_residual_j(),0.0,2e-11);
    })
}

#[test]
fn data_only_regional_loss_swap_changes_motion_and_heat_at_the_correct_elements() {
    with_cx(|cx,gate| {
        let mut nodes=NODES.to_vec();
        nodes.extend(NODES.map(|mut p|{p[0]+=2.0;p}));
        let tets=[[0,1,2,3],[4,5,6,7]];
        let fixed:Vec<usize>=(0..24).filter(|&d|d!=11 && d!=23).collect();
        // First source dissipates through Kelvin viscosity, second through
        // hereditary Maxwell memory. Both must reach their own element heat.
        let kelvin=material(0.01,0.0);let maxwell=material(0.0,0.5);
        let fields=[kelvin.clone(),maxwell.clone()];
        let swapped=[maxwell,kelvin];
        let a=system(SolidMeshRef::Tet4 {nodes_m:&nodes,elements:&tets},
            TetViscoField::PerElement(&fields),&fixed,cx);
        let b=system(SolidMeshRef::Tet4 {nodes_m:&nodes,elements:&tets},
            TetViscoField::PerElement(&swapped),&fixed,cx);
        let basis=[vec![1.0,0.0],vec![0.0,1.0]];
        let ab=TensorSolidBody::new(a.project_basis(&basis,projection_budget(),cx).unwrap(),
            dynamics_budget(),gate).unwrap();
        let bb=TensorSolidBody::new(b.project_basis(&basis,projection_budget(),cx).unwrap(),
            dynamics_budget(),gate).unwrap();
        let mut ah=ab.history(&[0.001,0.001],&[0.1,0.1],InitialMemory::Relaxed,gate).unwrap();
        let mut bh=bb.history(&[0.001,0.001],&[0.1,0.1],InitialMemory::Relaxed,gate).unwrap();
        let force=vec![[0.0;3];8];
        let astep=ah.advance(&force,0.02,1e-10,gate).unwrap();
        let bstep=bh.advance(&force,0.02,1e-10,gate).unwrap();
        assert!(astep.element_heat_j.iter().all(|x|*x>0.0));
        assert!((astep.element_heat_j[0]-astep.element_heat_j[1]).abs()>1e-9);
        close(astep.element_heat_j[0],bstep.element_heat_j[1],1e-11);
        close(astep.element_heat_j[1],bstep.element_heat_j[0],1e-11);
        close(ah.material_history().displacement_m()[0],bh.material_history().displacement_m()[1],1e-11);
        close(astep.element_heat_j.iter().sum(),astep.dynamics.dissipated,1e-11);
        let state=ah.material_history().state().to_vec();
        let time=ah.material_history().time_s();let heat=ah.material_history().dissipated_heat_j();
        assert!(ah.advance(&force,0.01,0.0,gate).is_err());
        assert_eq!(ah.material_history().state(),state.as_slice());
        assert_eq!(ah.material_history().time_s(),time);
        assert_eq!(ah.material_history().dissipated_heat_j(),heat);
    })
}

#[test]
fn constrained_nonorthogonal_source_basis_preserves_nodal_port_work() {
    with_cx(|cx,gate| {
        let material=material(0.001,0.1);
        let fixed:Vec<usize>=(0..12).filter(|&d|d!=7 && d!=11).collect();
        let source=system(SolidMeshRef::Tet4 {nodes_m:&NODES,elements:&TETS},
            TetViscoField::Uniform(&material),&fixed,cx);
        let basis=[vec![1.0,0.5],vec![-0.4,1.2]];
        let projected=source.project_basis(&basis,projection_budget(),cx).unwrap();
        assert!(projected.mass_matrix_kg()[1].abs()>0.0);
        let body=TensorSolidBody::new(projected,dynamics_budget(),gate).unwrap();
        let mut history=body.history(&[0.0,0.0],&[0.0,0.0],InitialMemory::Relaxed,gate).unwrap();
        let mut force=[[0.0;3];4];force[0]=[100.0;3];force[2][1]=0.02;force[3][2]=-0.01;
        let step=history.advance(&force,0.01,1e-10,gate).unwrap();
        assert_eq!(step.mean_nodal_velocity_m_s[0],[0.0;3]);
        let work=0.01*force.iter().zip(&step.mean_nodal_velocity_m_s)
            .map(|(f,v)|dot(f,v)).sum::<f64>();
        close(work,step.dynamics.supplied,1e-13);
        close(history.material_history().stored_energy_j()+step.element_heat_j.iter().sum::<f64>(),
            work,1e-12);
        // Dependent columns are rejected by full SPD inertia, never repaired.
        let dependent=source.project_basis(&[vec![1.0,0.0],vec![1.0,0.0]],projection_budget(),cx).unwrap();
        assert!(TensorSolidBody::new(dependent,dynamics_budget(),gate).is_err());
    })
}
