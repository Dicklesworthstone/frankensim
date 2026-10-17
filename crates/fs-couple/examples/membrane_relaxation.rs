//! Force-driven P1 membrane strip, using explicitly SYNTHETIC material data.
//! cargo run -p fs-couple --example membrane_relaxation -- [cells:1..4] [steps:1..4096]
//! CSV is mechanics/heat, not microphone pressure or a measured-specimen claim.
use fs_couple::tensor_material::{InitialMemory, TensorMaterialBudget};
use fs_couple::tensor_material::membrane::{MembraneEdgeTraction, MembraneFacetMaterial, TensorMembrane};
use fs_exec::CancelGate;
use fs_material::{IsotropicElastic, visco_tensor::{TensorProny, TensorPronyBranch, TensorProjectionBudget}};
use fs_plate::shell::ShellMesh;
use std::io::Write;

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()>2 {return Err("usage: membrane_relaxation [cells:1..4] [steps:1..4096]".into());}
    let cells=args.first().map(String::as_str).unwrap_or("2").parse::<usize>()?;
    let steps=args.get(1).map(String::as_str).unwrap_or("32").parse::<usize>()?;
    if !(1..=4).contains(&cells) || !(1..=4096).contains(&steps) {
        return Err("cells must be 1..4 and steps 1..4096; no silent mesh truncation".into());
    }
    let nodes:Vec<_>=(0..=cells).flat_map(|i|[[i as f64/cells as f64,0.0,0.0],
        [i as f64/cells as f64,0.1,0.0]]).collect();
    let tris=(0..cells).flat_map(|i|[[2*i,2*i+2,2*i+3],[2*i,2*i+3,2*i+1]]).collect();
    // Suppress transverse motion explicitly. This is an axial membrane strip,
    // not a bending plate or a prestressed transverse drumhead.
    let mut fixed:Vec<_>=(0..nodes.len()).flat_map(|i|[3*i+1,3*i+2]).collect();
    fixed.extend([0,3]);
    let c=IsotropicElastic::new(1000.0,0.0,0.1)?.stiffness();
    let law=TensorProny::from_tensor_shear(c,vec![TensorPronyBranch::from_stiffness_mandel(
        c.map(|r|r.map(|v|0.4*v)),0.02)?],0.1)?;
    let material=MembraneFacetMaterial {law,thickness_m:0.1,density_kg_m3:10.0};
    let gate=CancelGate::new_clock_free();
    let mesh=TensorMembrane::new(ShellMesh::new(nodes,tris)?,vec![material;2*cells],&fixed,
        TensorProjectionBudget {max_coordinates:128,max_scalar_terms:100_000_000},
        TensorMaterialBudget {max_state_coordinates:128,max_material_points:128},&gate)?;
    let zeros=vec![[0.0;3];mesh.mesh().nodes.len()];
    let mut history=mesh.history(&zeros,&zeros,InitialMemory::Relaxed,&gate)?;
    let unloaded_faces=vec![[[0.0;3];3];mesh.mesh().tris.len()];
    let edge=MembraneEdgeTraction {nodes:[2*cells,2*cells+1],values_n_m:[[0.02,0.0,0.0];2]};
    let pulse=mesh.distributed_nodal_forces(&unloaded_faces,&[edge],[0.0;3],&gate)?;
    eprintln!("synthetic membrane: cells={cells}, nodes={}, triangles={}, free_displacements={}, retained_state={}, mass_kg={}, dt_s=0.002; no bending/prestress/radiation",
        mesh.mesh().nodes.len(),mesh.mesh().tris.len(),mesh.body().coordinates(),
        mesh.body().state_dimension(),mesh.nodal_masses_kg().iter().sum::<f64>());
    let mut out=std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(out,"time_s,tip_x_m,stored_j,work_j,heat_j,balance_residual_j")?;
    for i in 0..steps {
        mesh.advance_nodal(&mut history,if i<4 {&pulse} else {&zeros},0.002,&gate)?;
        let u=mesh.expand_coordinates(history.displacement_m())?;
        let tip=0.5*(u[2*cells][0]+u[2*cells+1][0]);
        writeln!(out,"{:.17e},{tip:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",history.time_s(),
            history.stored_energy_j(),history.supplied_work_j(),history.dissipated_heat_j(),history.balance_residual_j())?;
    }
    out.flush()?;
    Ok(())
}
