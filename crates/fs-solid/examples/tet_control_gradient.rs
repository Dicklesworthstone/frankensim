//! Two traction-pulse amplitudes designed with the native trajectory adjoint.
//! cargo run -p fs-solid --example tet_control_gradient --release
//! Synthetic coarse P1 cantilever; no continuum beam-accuracy qualification.
//! The objective is 1/2*(tip-target)^2 + 1/2*(tau*tip_velocity)^2
//! + 1/2*regularization*|pressure_amplitudes|^2, in m^2. Fixed geometry,
//! material, damping, time grid and pulse shapes make it a convex quadratic.
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3,SurfaceLoad3,VolumeLoad3};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::MotionCotangent3;
use std::error::Error;

fn pulse(time_s:f64,index:usize)->f64 {
    (1.-((time_s-[0.01,0.03][index])/0.01).abs()).max(0.)
}
fn main()->Result<(),Box<dyn Error>> {
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena|->Result<(),Box<dyn Error>> {
        let cx=Cx::new(&gate,arena,StreamKey {seed:92,kernel_id:92,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let nodes=[[0.,0.,0.],[0.,0.,0.02],[0.,0.02,0.],[0.,0.02,0.02],
            [0.1,0.,0.],[0.1,0.,0.02],[0.1,0.02,0.],[0.1,0.02,0.02]];
        let tets=[[0,4,6,7],[0,6,2,7],[0,2,3,7],[0,3,1,7],[0,1,5,7],[0,5,4,7]];
        let material=TetElasticMaterial::try_new(1000.,1e6,0.25,[92;32])?;
        let fixed:Vec<_>=(0..12).collect();
        let problem=TetLinearElasticProblem {nodes_m:&nodes,tetrahedra:&tets,
            materials:TetMaterialField::Uniform(&material),fixed_dofs:&fixed,budget:Default::default()};
        let geometry=TetLoadGeometry3::new(&nodes,&tets,problem.budget,24,&cx)?;
        let faces:Vec<_>=geometry.faces().iter().enumerate()
            .filter(|(_,f)|f.nodes.iter().all(|&i|nodes[i][0]==0.1)).map(|(i,_)|i).collect();
        let unit=geometry.assemble(&[SurfaceLoad3::Traction {faces:&faces,pascals:[0.,0.,1.]}],VolumeLoad3::None,&[],&cx)?;
        let empty=geometry.assemble(&[],VolumeLoad3::None,&[],&cx)?;
        let cfg=TransientConfig3 {maximum_steps:200,mass_damping_per_s:0.2,
            solve:TetStaticSolveConfig {relative_residual_tolerance:1e-12,..Default::default()},..Default::default()};
        let initial=empty.begin_transient(&problem,None,0.,&[[0.;3];8],&[[0.;3];8],cfg,&cx)?;
        let initial=initial.record_controls(Default::default(),Default::default(),&cx)?;
        let target=1e-5;let tau=0.005;let regularization=1e-16;
        let evaluate=|parameters:[f64;2]|->Result<(f64,[f64;2],f64,f64),Box<dyn Error>> {
            let mut recording=initial.clone();
            for _ in 0..200 {
                let t=recording.run().state().time_s+0.0002;
                let p=parameters[0]*pulse(t,0)+parameters[1]*pulse(t,1);
                let load=geometry.assemble(&[SurfaceLoad3::Traction {faces:&faces,pascals:[0.,0.,p]}],VolumeLoad3::None,&[],&cx)?;
                recording.advance(&load,0.0002,&cx)?;
            }
            let n=recording.run().free_dofs().len();
            let tip=recording.run().free_dofs().binary_search(&23).map_err(|_|"missing free tip component")?;
            let q=recording.run().state().displacement[tip];let v=recording.run().state().velocity[tip];
            let mut seed=MotionCotangent3 {q:vec![0.;n],v:vec![0.;n],a:vec![0.;n]};
            seed.q[tip]=q-target;seed.v[tip]=tau*tau*v;
            let pullback=recording.reverse(&seed,&cx)?;
            let mut gradient=parameters.map(|p|regularization*p);
            for sample in &pullback.forces {
                let pressure=recording.adjoint().contract_load_direction(&sample.values,&unit,&cx)?;
                for j in 0..2 {gradient[j]+=pulse(sample.time_s,j)*pressure;}
            }
            let objective=0.5*(q-target).powi(2)+0.5*(tau*v).powi(2)
                +0.5*regularization*parameters.iter().map(|p|p*p).sum::<f64>();
            Ok((objective,gradient,q,v))
        };
        let zero=evaluate([0.,0.])?;
        let x=evaluate([1.,0.])?;let y=evaluate([0.,1.])?;
        // Gradient differences give the exact discrete quadratic Hessian.
        let h00=x.1[0]-zero.1[0];let h10=x.1[1]-zero.1[1];
        let h01=y.1[0]-zero.1[0];let h11=y.1[1]-zero.1[1];
        let det=h00*h11-h01*h10;
        if !det.is_finite() || det<=0. {return Err("control Hessian is not positive definite".into());}
        let optimum=[(-h11*zero.1[0]+h01*zero.1[1])/det,(h10*zero.1[0]-h00*zero.1[1])/det];
        let answer=evaluate(optimum)?;
        if answer.0>=zero.0 {return Err("control objective did not improve".into());}
        println!("pulse0_pa,pulse1_pa,initial_objective_m2,optimized_objective_m2,tip_z_m,tip_vz_m_s,gradient0,gradient1");
        println!("{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            optimum[0],optimum[1],zero.0,answer.0,answer.2,answer.3,answer.1[0],answer.1[1]);
        Ok(())
    })
}
