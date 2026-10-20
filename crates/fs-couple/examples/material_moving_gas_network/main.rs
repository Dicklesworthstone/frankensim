//! Two material-inertia walls with common finite inlet/outlet gas plumbing.
//! Synthetic material inputs; no experimental gas or structural fidelity claim.
#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

use std::{collections::BTreeMap, error::Error, io::Read};
use fs_couple::cavity_phs::piston::{PistonInput, PistonSolve};
use fs_couple::cavity_phs::piston::material::{MaterialGasPiston, graph::{
    step_candidate, MovingGasDrive, MovingGasSolve,
}};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_material::thermochemical_gas::energy::parcel::transfer::{snapshot::GasSnapshotSolve,
    nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink}};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{ResolvedBody, RegionBinding, VolumeMoments, resolve_body};
use fs_qty::{Density, Temperature};
type Result<T> = std::result::Result<T, Box<dyn Error>>;
type Walls = BTreeMap<u64,MaterialGasPiston>;
type Reservoirs = BTreeMap<u64,ClosedGasParcel>;
const EXTERIOR: f64 = 100000.0;

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn body(density: f64) -> Result<ResolvedBody> {
    if !density.is_finite() || density <= 0.0 { return Err("positive finite synthetic wall density required".into()); }
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new("density", Density::DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0,density),(400.0,density)], dims: Density::DIMS },
        validity: ValidityDomain::unconstrained().with("T",200.0,400.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "SYNTHETIC example density input".into(),
            license: "CC0-1.0".into(), artifact: None } })?;
    let material = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic wall".into(),
        phase: "solid".into(), process: "example".into(), revision: 0 },claims,vec![])?;
    let point = QueryPoint::new().with("T",300.0)?;
    Ok(resolve_body(&[RegionBinding { key: "wall",material: &material,point: &point,
        geometry: VolumeMoments::cuboid(Vec3::new(0.1,0.1,0.01))?,
        selection: MaterialPropertySelection::SingleClaimOnly }],1,||false)?)
}
fn initial(left: &str, right: &str, provenance: &str, density: f64) -> Result<(Walls,Reservoirs)> {
    // The paired files share one explicitly declared source set. Inconsistent
    // same-species source laws still fail in the original inventory owner.
    let mut supply = ClosedGasParcel::new(source::load(left,provenance)?,source::point(400.0,EXTERIOR)?,0.001)?;
    supply.add_heat(150.0,energy(),||false)?;
    let sink = ClosedGasParcel::new(source::load(right,provenance)?,source::point(400.0,EXTERIOR)?,0.001)?;
    let mut chamber = sink.clone(); chamber.add_heat(30.0,energy(),||false)?;
    Ok((BTreeMap::from([
        (2,MaterialGasPiston::new(chamber.clone(),body(density)?,0.01,0.02)?),
        (3,MaterialGasPiston::new(chamber,body(8000.0)?,0.01,0.02)?)]),
        BTreeMap::from([(1,supply),(4,sink)])))
}
fn links() -> Vec<SnapshotNozzleLink> {
    [(1,2),(1,3),(2,4),(3,4)].into_iter().enumerate().map(|(i,(a,b))| SnapshotNozzleLink {
        id: i as u64, parcels: [a,b], nozzle: GasNozzle { area_m2: 1e-5, discharge_coefficient: 0.8 },
        nozzle_solve: energy() }).collect()
}
fn advance(walls: &Walls,reservoirs: &Reservoirs,dt: f64,heat_w: f64) -> Result<(Walls,Reservoirs,f64)> {
    let controls = PistonSolve { volume_bracket_m3: [0.0007,0.0015], gas: energy(),
        momentum_tolerance_n_s: 2e-7*dt, energy_tolerance_j: 1e-7, maximum_iterations: 128 };
    let drives = [2,3].map(|id|(id,MovingGasDrive { input: PistonInput { duration_s: dt,
        heat_j: heat_w*dt, exterior_pressure_pa: EXTERIOR, external_force_n: 0.0 }, solve: controls })).into();
    let endpoints = (1..=4).map(|id|(id,energy())).collect();
    let solve = MovingGasSolve { exchange: GasSnapshotSolve { maximum_withdrawn_fraction: 0.1,
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }, total_energy_tolerance_j: 1e-7 };
    let candidate = step_candidate(walls,reservoirs,&links(),dt,&endpoints,&drives,solve,||false)?;
    let (walls,reservoirs,report) = candidate.into_parts();
    Ok((walls,reservoirs,report.total_energy_residual_j))
}
fn total_energy(walls: &Walls,reservoirs: &Reservoirs) -> f64 {
    reservoirs.values().map(|g|g.mass_kg()*g.state().specific_internal_energy_j_kg()).sum::<f64>()
        + walls.values().map(|w| w.gas().mass_kg()*w.gas().state().specific_internal_energy_j_kg()
            + 0.5*w.piston().mass_kg()*w.piston().velocity_m_s().powi(2)
            + EXTERIOR*w.gas().volume_m3()).sum::<f64>()
}
fn read_source(path: &str) -> Result<String> {
    let mut text=String::new();
    std::fs::File::open(path)?.take(65537).read_to_string(&mut text)?;
    if text.len()>65536 { return Err("synthetic gas source exceeds 64 KiB".into()); }
    Ok(text)
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(3..=5).contains(&args.len()) {
        return Err("usage: material_moving_gas_network SUPPLY.tsv CHAMBER.tsv WALL_DENSITY_kg_m3 [steps=40] [dt_s=0.00005]".into());
    }
    let density: f64=args[2].parse()?;
    let steps: usize=args.get(3).map_or(Ok(40),|s|s.parse())?;
    let dt: f64=args.get(4).map_or(Ok(0.00005),|s|s.parse())?;
    if !(1..=10000).contains(&steps) || !dt.is_finite() || dt<=0.0 || !(dt*steps as f64).is_finite() {
        return Err("bounded count and positive representable physical clock required".into());
    }
    let provenance=format!("paired synthetic gas source set: {} + {}",args[0],args[1]);
    let (mut walls,mut reservoirs)=initial(&read_source(&args[0])?,&read_source(&args[1])?,&provenance,density)?;
    let initial_energy=total_energy(&walls,&reservoirs); let mut supplied=0.0;
    println!("time_s,id,moving,gas_mass_kg,volume_m3,temperature_k,pressure_pa,sound_speed_m_s,viscosity_pa_s,conductivity_w_m_k,wall_mass_kg,wall_velocity_m_s,supplied_heat_j,network_energy_defect_j,step_energy_residual_j");
    for i in 0..=steps {
        let residual=if i==0 {0.0} else {
            let power=if i<=steps/2 {100.0} else {0.0};
            let (next_w,next_r,residual)=advance(&walls,&reservoirs,dt,power)?;
            walls=next_w;reservoirs=next_r;supplied+=2.0*power*dt;residual
        };
        let defect=total_energy(&walls,&reservoirs)-initial_energy-supplied;
        for id in 1..=4 {
            let (gas,moving,mass,velocity)=if let Some(wall)=walls.get(&id) {
                (wall.gas(),true,wall.piston().mass_kg(),wall.piston().velocity_m_s())
            } else {(&reservoirs[&id],false,0.0,0.0)};
            let state=gas.state().state();
            println!("{:.9},{},{},{:.12e},{:.12e},{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.9},{:.12e},{:.12e},{:.12e},{:.12e}",
                i as f64*dt,id,u8::from(moving),gas.mass_kg(),gas.volume_m3(),state.temperature,
                state.pressure,state.sound_speed,state.dynamic_viscosity,state.thermal_conductivity,
                mass,velocity,supplied,defect,residual);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const A:&str=include_str!("../../../fs-material/examples/material_gas_spring/mixture-a.tsv");
    const B:&str=include_str!("../../../fs-material/examples/material_gas_spring/mixture-b.tsv");
    fn run(left:&str,right:&str,density:f64)->(Walls,Reservoirs) {
        let (mut w,mut r)=initial(left,right,"paired synthetic test sources",density).unwrap();
        let initial=total_energy(&w,&r);let mut heat=0.0;
        for i in 0..40 {
            let power=if i<20 {100.0} else {0.0};
            let (nw,nr,residual)=advance(&w,&r,0.00005,power).unwrap();
            w=nw;r=nr;heat+=2.0*power*0.00005;
            assert!(residual.abs()<1e-7);
            assert!((total_energy(&w,&r)-initial-heat).abs()<1e-6);
        }
        (w,r)
    }
    #[test]
    fn original_gas_files_and_solid_density_change_real_network_trajectories() {
        let (light,lr)=run(A,B,1000.0);let (dense,dr)=run(A,B,8000.0);
        assert!(light[&2].piston().velocity_m_s()>dense[&2].piston().velocity_m_s());
        assert_ne!(lr[&1].mass_kg(),dr[&1].mass_kg());
        let (swapped,_)=run(B,A,1000.0);
        assert_ne!(light[&2].gas().state().state().sound_speed,swapped[&2].gas().state().state().sound_speed);
        assert_eq!(light[&2].piston().mass_kg(),light[&2].body().inertia().mass());
    }
    #[test]
    fn invalid_density_and_oversized_clock_do_not_reset_or_partially_advance_sources() {
        for rho in [0.0,-1.0,f64::NAN,f64::INFINITY] { assert!(body(rho).is_err()); }
        let (w,r)=initial(A,B,"paired synthetic test sources",1000.0).unwrap();
        let saved=(w.clone(),r.clone());
        assert!(advance(&w,&r,1.0,100.0).is_err());
        assert_eq!((w,r),saved);
    }
}
