//! A small, explicitly synthetic cantilever with a ramped transverse traction.
//! Run: cargo run -p fs-solid --example tet_transient --release > motion.csv
//! The coarse P1 mesh demonstrates integration; it does not qualify beam accuracy.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{FrameTree, WORLD};
use fs_scenario::payload::{OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::history::{SurfaceHistory3, TetLoadHistory3, VolumeHistory3, TRACTION_DIMS};
use fs_solid::loads3::transient::TransientConfig3;
use std::io::{self, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 73, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let nodes = [[0.,0.,0.],[0.,0.,0.02],[0.,0.02,0.],[0.,0.02,0.02],
            [0.1,0.,0.],[0.1,0.,0.02],[0.1,0.02,0.],[0.1,0.02,0.02]];
        let tets = [[0,4,6,7],[0,6,2,7],[0,2,3,7],[0,3,1,7],[0,1,5,7],[0,5,4,7]];
        let material = TetElasticMaterial::try_new(1000.0,1e6,0.25,[73;32])?;
        let fixed: Vec<_> = (0..12).collect();
        let problem = TetLinearElasticProblem { nodes_m:&nodes, tetrahedra:&tets,
            materials:TetMaterialField::Uniform(&material), fixed_dofs:&fixed, budget:TetAssemblyBudget::default() };
        let geometry = TetLoadGeometry3::new(&nodes,&tets,problem.budget,24,&cx)?;
        let tip_faces: Vec<_> = geometry.faces().iter().enumerate()
            .filter(|(_,f)| f.nodes.iter().all(|&i| nodes[i][0] == 0.1)).map(|(i,_)| i).collect();
        let meta = PayloadMeta::new(QuantityContract::Dimensions(TRACTION_DIMS),
            PayloadId::new("basis/world-cartesian")?,WORLD,OrientationParity::Even,ReferenceSemantics::Continuous)?;
        let times = [0.0,0.01,0.02].into_iter().map(|t| QtyAny::new(t,Dims([0,0,1,0,0,0]))).collect();
        let values = [0.0,100.0,0.0].into_iter()
            .map(|z| [0.0,0.0,z].map(|v| QtyAny::new(v,TRACTION_DIMS)).to_vec()).collect();
        let traction = VectorPayload::new(meta,SampleSource::table(times,values,
            TableInterpolation::Linear,OutsideDomainPolicy::Clamp)?)?;
        let surfaces = [SurfaceHistory3::Traction { faces:&tip_faces,value:&traction }];
        let frames = FrameTree::new();
        let history = TetLoadHistory3 { surfaces:&surfaces, volume:VolumeHistory3::Fixed(VolumeLoad3::None),
            nodal_forces_n:&[],frames:&frames };
        let mut run = history.begin_transient(&geometry,&problem,None,0.0,&[[0.0;3];8],&[[0.0;3];8],
            TransientConfig3 { maximum_steps:1000, mass_damping_per_s:0.2,
                maximum_displacement_m:0.01,..Default::default() },&cx)?;
        let stdout = io::stdout(); let mut output = io::BufWriter::new(stdout.lock());
        writeln!(output,"time_s,tip_z_m,tip_vz_m_s,tip_az_m_s2,quadratic_energy_j,balance_error_j,element0_sigma_xx_pa,support_z_n,elastic_energy_j,free_force_residual_n")?;
        for index in 0..1000 {
            let report = run.advance_history(&history,1e-4,&cx)?;
            if index % 10 == 0 || index == 999 {
                let motion = run.nodal_motion(&cx)?;
                let response = run.recover_response(fs_blake3::ContentHash([73;32]),0.01,&cx)?;
                let support_z: f64 = response.support_reactions_n.iter().map(|r| r[2]).sum();
                writeln!(output,"{:.8},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                    report.time_s,motion[0][7][2],motion[1][7][2],motion[2][7][2],
                    report.quadratic_energy_j,report.balance_error_j,
                    response.stress.elements()[0].stress_mandel_pa()[0],support_z,
                    response.stress.elastic_energy_j(),response.free_force_residual_n)?;
            }
        }
        output.flush()?; Ok(())
    })
}
