//! Existing density resolver, existing multiport gas owner and actual motion.
use super::*;
use crate::cavity_phs::piston::tests::{controls, drive, gas};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_qty::{Density, Temperature};

fn body(density: f64) -> fs_mbd_material::ResolvedBody {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new("density", Density::DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0,density),(400.0,density)], dims: Density::DIMS },
        validity: ValidityDomain::unconstrained().with("T",200.0,400.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic network wall density".into(),
            license: "CC0-1.0".into(), artifact: None } }).unwrap();
    let material = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic wall".into(),
        phase: "solid".into(), process: "authored".into(), revision: 0 }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T",300.0).unwrap();
    resolve_body(&[RegionBinding { key: "wall", material: &material, point: &point,
        geometry: VolumeMoments::cuboid(Vec3::new(0.1,0.1,0.01)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly }],1,||false).unwrap()
}
fn solve() -> MovingGasSolve {
    MovingGasSolve { exchange: GasSnapshotSolve { maximum_withdrawn_fraction: 0.1,
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }, total_energy_tolerance_j: 1e-7 }
}
fn fixture(density: f64) -> (BTreeMap<u64,MaterialGasPiston>, BTreeMap<u64,ClosedGasParcel>, Vec<SnapshotNozzleLink>) {
    let mut supply = gas(0.25,0.001); supply.add_heat(150.0,controls().gas,||false).unwrap();
    let mut cavity = gas(0.75,0.001); cavity.add_heat(30.0,controls().gas,||false).unwrap();
    let walls = BTreeMap::from([
        (2,MaterialGasPiston::new(cavity.clone(),body(density),0.01,0.02).unwrap()),
        (3,MaterialGasPiston::new(cavity,body(8000.0),0.01,0.02).unwrap())]);
    let links = [(1,2),(1,3),(2,4),(3,4)].into_iter().enumerate().map(|(i,(a,b))| SnapshotNozzleLink {
        id: i as u64, parcels: [a,b], nozzle: GasNozzle { area_m2: 1e-5, discharge_coefficient: 0.8 },
        nozzle_solve: controls().gas }).collect();
    (walls,BTreeMap::from([(1,supply),(4,gas(0.75,0.001))]),links)
}
fn endpoint_controls() -> BTreeMap<u64,GasEnergySolve> { (1..=4).map(|id|(id,controls().gas)).collect() }
fn drives(dt: f64) -> BTreeMap<u64,MovingGasDrive> {
    [2,3].map(|id|(id,MovingGasDrive { input: super::super::super::PistonInput { duration_s: dt,
        heat_j: 100.0*dt, ..drive(0.0) }, solve: controls() })).into()
}

#[test]
fn source_bound_network_is_identical_to_shared_numeric_path_and_retains_each_body() {
    let (walls,reservoirs,links) = fixture(1000.0);
    let numeric = walls.iter().map(|(&id,w)|(id,w.piston().clone())).collect();
    let dt = drive(0.0).duration_s;
    let expected = graph::step_candidate(&numeric,&reservoirs,&links,dt,
        &endpoint_controls(),&drives(dt),solve(),||false).unwrap();
    let actual = step_candidate(&walls,&reservoirs,&links,dt,
        &endpoint_controls(),&drives(dt),solve(),||false).unwrap();
    assert_eq!(actual.reservoirs(),expected.reservoirs());
    assert_eq!(actual.step(),expected.step());
    for (&id,wall) in actual.pistons() {
        assert_eq!(wall.body(),walls[&id].body());
        assert_eq!(wall.piston(),&expected.pistons()[&id]);
        assert_eq!(wall.piston().ports,expected.pistons()[&id].ports);
        assert_eq!(wall.piston().mass_kg(),wall.body().inertia().mass());
        assert_eq!(wall.body().regions().len(),1);
    }
}

fn run(density: f64, count: usize) -> MovingGasCandidate<MaterialGasPiston> {
    let (mut walls,mut reservoirs,links) = fixture(density);
    let dt = 0.002/count as f64;
    for i in 0..count {
        let candidate = step_candidate(&walls,&reservoirs,&links,dt,
            &endpoint_controls(),&drives(dt),solve(),||false).unwrap();
        if i == count-1 { return candidate; }
        let (w,r,_) = candidate.into_parts(); walls=w;reservoirs=r;
    }
    unreachable!("test uses nonzero count")
}

#[test]
fn only_solid_density_changes_motion_and_then_flow_and_the_common_gas_medium() {
    let light = run(1000.0,40); let dense = run(8000.0,40);
    let (a,b) = (&light.pistons()[&2],&dense.pistons()[&2]);
    assert!(a.piston().velocity_m_s() > b.piston().velocity_m_s());
    assert!(a.gas().volume_m3() > b.gas().volume_m3());
    assert_ne!(a.gas().state().state().temperature,b.gas().state().state().temperature);
    assert_ne!(a.gas().state().state().sound_speed,b.gas().state().state().sound_speed);
    // Motion changes next-step pressure and therefore the shared finite supply.
    assert_ne!(light.reservoirs()[&1].mass_kg(),dense.reservoirs()[&1].mass_kg());
    assert_eq!(light.pistons()[&3].body(),dense.pistons()[&3].body());
    for (p,density) in [(a,1000.0),(b,8000.0)] {
        assert_eq!(p.body().regions()[0].material().property("density").unwrap().value_si(),density);
    }
}

#[test]
fn actual_material_motion_refines_with_the_declared_first_order_split() {
    let velocities: Vec<_> = [32,64,128].into_iter().map(|n| run(8000.0,n).pistons()[&2].piston().velocity_m_s()).collect();
    let ratio = ((velocities[0]-velocities[1])/(velocities[1]-velocities[2])).abs();
    assert!(ratio > 1.7 && ratio < 2.3,"{velocities:?}, ratio {ratio}");
}

#[test]
fn final_material_publication_cancellation_preserves_sources_and_prior_port_history() {
    let (walls,reservoirs,links) = fixture(1000.0);
    let saved = walls.clone(); let dt=drive(0.0).duration_s;
    let mut calls=0;
    step_candidate(&walls,&reservoirs,&links,dt,&endpoint_controls(),&drives(dt),solve(),
        || { calls+=1;false }).unwrap();
    let mut count=0;
    assert_eq!(step_candidate(&walls,&reservoirs,&links,dt,&endpoint_controls(),&drives(dt),solve(),
        || { count+=1;count==calls }),Err(MovingGasError::Cancelled));
    let mut bad = drives(dt); bad.get_mut(&3).unwrap().input.heat_j=1e12;
    assert!(step_candidate(&walls,&reservoirs,&links,dt,&endpoint_controls(),&bad,solve(),||false).is_err());
    assert_eq!(walls,saved);
    for id in walls.keys() { assert_eq!(walls[id].piston().ports,saved[id].piston().ports); }
}
