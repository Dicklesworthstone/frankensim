use super::*;
use fs_material::gas::{GasSpec, GasState};
fn loss(gas:GasState)->ContractionLoss {
    ContractionLoss::new(0.01,0.004,gas,ContractionSpec {area_ratio:0.6,
        maximum_contracted_mach:0.2,provenance:"authored port-equation fixture".into()}).unwrap()
}

#[test]
fn nonlinear_mouth_obeys_both_pressure_faces_and_work_including_extreme_reflection() {
    for t in [273.15,293.15,333.15] { for p in [50000.0,101325.0,180000.0] {
        let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),t,p).unwrap();
        let loss=loss(gas);let z=gas.density*gas.sound_speed/(core::f64::consts::PI*0.006*0.006);
        let (l,r,dt)=(100.0,10000.0,1.0/48000.0);
        for reflection in [-1.0,-0.999999999999,-0.5,0.0,0.7,0.999999999999,1.0] {
            for q in [-1e-4_f64,0.0,1e-4] {for old in [-5e-5,0.0,5e-5] {
                let outgoing=37.0;
                let incoming=outgoing-z*q;
                let history=incoming-reflection*outgoing;
                let dp=0.5*gas.density*q*q.abs()/(0.6*0.01*0.004_f64).powi(2);
                let source=l/dt*(q-old)+r*q+dp+outgoing+incoming;
                let (step,power)=solve(WindBoundary {direct_reflection:reflection,history_pressure_pa:history},
                    source,old,l,r,z,dt,&loss).unwrap();
                assert!((step.flow_m3_s-q).abs()<1e-15);
                assert_eq!(step.incoming_pa,reflection*step.outgoing_pa+history);
                let work=source*step.flow_m3_s*dt;
                let pipe=step.pressure_pa*step.flow_m3_s*dt;
                let linear=r*step.flow_m3_s.powi(2)*dt;
                let storage=0.5*l*(step.flow_m3_s.powi(2)-old*old);
                let numerical=0.5*l*(step.flow_m3_s-old).powi(2);
                let separation=power*dt;
                let scale=work.abs()+pipe.abs()+linear+storage.abs()+numerical+separation;
                assert!((work-pipe-linear-storage-numerical-separation).abs()<2e-12*scale.max(1e-30));
                assert!(separation>=0.0);
            }}
        }
    }}
}

#[test]
fn nonlinear_loss_can_resolve_pressure_release_without_inventing_linear_damping() {
    let g=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap();
    let law=loss(g);
    let (step,power)=solve(WindBoundary {direct_reflection:-1.0,history_pressure_pa:0.0},
        100.0,0.0,0.0,0.0,1e6,1.0/48000.0,&law).unwrap();
    assert_eq!(step.pressure_pa,0.0);
    assert!((power-100.0*step.flow_m3_s).abs()<1e-12*power);
    assert!(solve(WindBoundary {direct_reflection:1.01,history_pressure_pa:0.0},
        1.0,0.0,1.0,0.0,1e6,1.0/48000.0,&law).is_err());
}
