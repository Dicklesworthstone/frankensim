//! G1/G3: original heat transport changes the continuing Maxwell/fluid/circuit state.
//! Manufactured two-cell sources; parallel plates are an explicitly imposed limit.
use super::*;
use crate::circuit_enthalpy::contacts::{heat::*, material::ElectricalContactNetwork};
use crate::circuit_enthalpy::material::ConductionLink;
use fs_conduction::{GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity, ViewFactorMatrix,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermal::{interface::ThermalInterfacePoint, conduction::{LayerConductionConfig,
    network::{LayerLinkBinding, PreparedLayerConduction}, source::ResolvedConductionLayer}};

const LINK: u64 = 19;
const LINK_AREA: f64 = 0.0005;
const FACE_SIDE: f64 = 0.01;
const DT: f64 = 1.0/256.0;
const SIGMA: f64 = 5.670374419e-8;
#[derive(Clone, Copy)]
struct Spec { powered: bool, paths: bool, integrated: bool, emissivity: f64, resistance: f64, upper: f64, source_w: f64 }
fn spec() -> Spec { Spec { powered: true, paths: true, integrated: false,
    emissivity: 0.2, resistance: 1.5, upper: 340.0, source_w: 0.0 } }
fn solve() -> FrameThermalSolve {
    let mut c = controls(); c.maximum_material_trials = 256; c.maximum_newton_iterations = 20;
    c.difference_step_j = 1e-10; c.cell_energy_tolerance_j = 1e-12;
    c.thermal_energy_tolerance_j = 1e-12; c.newton.absolute_tolerance = 1e-13; c
}
fn input() -> RelaxationElectricalInput<'static> {
    RelaxationElectricalInput { ports: &[], solve: electrical_controls() }
}
fn source(layer: usize, s: Spec) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = powered_source(layer, 1.0, false);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() { claims.insert_claim(claim.clone()).unwrap(); }
    for (name, dims, value, slope, upper) in [
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, [0.04,0.08][layer],
            if s.integrated { 0.03 } else { 0.0 }, 340.0),
        (fs_conduction::SURFACE_EMISSIVITY_PROPERTY, fs_conduction::EMISSIVITY_DIMS,
            if layer == 0 { s.emissivity } else { 0.6 }, 0.0, if layer == 0 { s.upper } else { 340.0 }),
    ] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name,dims),
            value: if slope == 0.0 { PropertyValue::Scalar { value,dims } } else {
                PropertyValue::Curve { abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                    knots:vec![(290.0,value*(1.0-10.0*slope)),(340.0,value*(1.0+40.0*slope))],dims } },
            validity:ValidityDomain::unconstrained().with("T",290.0,upper),
            interpolation:if slope == 0.0 {InterpolationPolicy::ConstantWithinValidity} else {InterpolationPolicy::LinearInside},
            uncertainty:UncertaintyModel::Unstated,observations:vec![],provenance:Provenance {
                source:"manufactured coupled transport".into(),license:"test-only".into(),artifact:None},
        }).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(),claims,vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",lower:&point(290.0),upper:&point(340.0),
        selection:SelectionPolicy::SingleClaimOnly},SolidLiquidPhase::Solid,0.0,1e-6).unwrap();
    (card,curve)
}
fn fixture<R>(s: Spec, run: impl FnOnce(&MaterialFluidFrame<'_,FluidTransportState>,
    &MaterialNetwork<'_>, &MaterialNetwork<'_>, &[(usize,ElectricalPropertyBasis)],
    DepositionSimulation<'_>, &PreparedHeatNetwork<'_,'_>, &Cx<'_>) -> R) -> R {
    let (a,ac) = source(0,s); let (b,bc) = source(1,s); let cards = [a,b]; let curves = [ac,bc];
    let f = frame(&cards); let m = model(&f); let query = point(300.0);
    let mut events: Vec<_> = (0..2).map(|i| TimedEvent { time_s:0.0,event:Event::Birth(DepositedCell {
        id:CellId(i as u64),segment_id:0,body_id:0,material_id:MaterialId(i as u64),
        mass_kg:[3.0,5.0][i]*area(i),specific_h_j_kg:10.0*CP,
        geometry:CellGeometry {start_m:[0.0;3],end_m:[1.0,0.0,0.0],volume_m3:area(i),surface_area_m2:1.0},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:if i==0 {s.source_w} else {0.0}},
    })}).collect();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:LINK,a:CellId(0),b:CellId(1),area_m2:LINK_AREA,conductance_w_k:999.0})});
    let thermal = DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:128}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:DT,energy_tolerance_j:1e-11},||Ok(())).unwrap().0;
    let key = PropertyKey::new("thermal_conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells: [_;2] = core::array::from_fn(|i| MaterialCell {cell:CellId(i as u64),material_id:MaterialId(i as u64),
        curve:curves[i].curve(),point:ThermalMaterialPoint {material:&cards[i],context:&query,
            temperature_axis:"T",liquid_fraction_axis:None},conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let links = [ConductionLink {link:LINK,distance_a_m:0.05,distance_b_m:0.05}];
    let resistors: [_;2] = core::array::from_fn(|i| CellResistor {destination:ResistorCell {branch:i+1,cell:CellId(i as u64)},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:area(i)});
    let plain = MaterialNetwork {cells:&cells,links:&links,resistors:&[]};
    let powered = MaterialNetwork {cells:&cells,links:&links,resistors:&resistors};
    let bases = [(1,ElectricalPropertyBasis::Resistivity),(2,ElectricalPropertyBasis::Resistivity)];
    let electrical = ElectricalContactNetwork {bulk:if s.powered {&powered} else {&plain},contacts:&[],
        electrical_bases:if s.powered {&bases} else {&[]}};
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(RP,RD),value:PropertyValue::Scalar {value:s.resistance,dims:RD},
        validity:ValidityDomain::unconstrained().with("Ta",290.0,340.0).with("Tb",290.0,340.0),
        interpolation:InterpolationPolicy::ConstantWithinValidity,uncertainty:UncertaintyModel::Unstated,
        observations:vec![],provenance:Provenance {source:"manufactured thermal interface".into(),license:"test-only".into(),artifact:None}}).unwrap();
    let interface = InterfaceSystemCard::assemble(SurfaceSpec {material:cards[0].id().clone(),texture_frame:"A".into()},
        SurfaceSpec {material:cards[1].id().clone(),texture_frame:"B".into()},SystemContext {medium:"fixture".into(),
            third_body:None,environment:"fixed".into(),history:"synthetic".into()},claims,vec![]).unwrap();
    let context = QueryPoint::new().with("Ta",300.0).unwrap().with("Tb",300.0).unwrap();
    let contacts = [MaterialContactLink {link:LINK,point:ThermalInterfacePoint {interface:&interface,
        materials:[&cards[0],&cards[1]],context:&context,temperature_axes:["Ta","Tb"],
        liquid_fraction_axes:[None,None],mean_temperature_axis:None},
        sides:core::array::from_fn(|i| ContactCell {cell:cells[i].cell,material_id:cells[i].material_id,curve:cells[i].curve}),
        selection:ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)}];
    let radiating = [0,1].map(|i| RadiatingCell {surface:i,cell:cells[i].cell,material_id:cells[i].material_id,
        curve:cells[i].curve,material:cells[i].point,selection:MaterialPropertySelection::SingleClaimOnly});
    let (complex,positions) = fs_conduction::fixtures::box_grid([1,1,1],[FACE_SIDE;3]);
    let mesh = fs_conduction::mesh::ConductionMesh::new(complex,positions).unwrap();
    let surfaces = (0..2).map(|i| {
        let epsilon = SurfaceEmissivity::from_thermal_state("fixture",cells[i].point,
            cells[i].curve.state_at_specific_enthalpy(10.0*CP).unwrap(),MaterialPropertySelection::SingleClaimOnly).unwrap().0;
        RadiationSurface::new(&mesh,format!("side {i}"),|face|
            fs_conduction::fixtures::on_box_face(face.centroid[0],FACE_SIDE*i as f64),epsilon).unwrap()
    }).collect();
    let enclosure = GrayDiffuseEnclosure::new(surfaces,
        ViewFactorMatrix::infinite_parallel_plates(FACE_SIDE*FACE_SIDE).unwrap()).unwrap();
    let radiation = MaterialRadiation {enclosure:&enclosure,cells:&radiating,
        power_tolerance_w:1e-11,radiosity_tolerance_w_m2:1e-8};
    let assembly = MaterialHeatNetwork {electrical:&electrical,
        thermal_contacts:if s.paths && !s.integrated {&contacts} else {&[]},
        radiation:if s.paths {Some(&radiation)} else {None}};
    let layers = [0,1].map(|i| ResolvedConductionLayer::new(cells[i].curve,cells[i].point,&key,
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),SolidLiquidPhase::Solid,[290.0,340.0],||Ok(())).unwrap());
    let integrated = PreparedLayerConduction::new(&thermal,vec![LayerLinkBinding {link:LINK,
        layers:[&layers[0],&layers[1]],distances_m:[0.05,0.05],contact_area_resistance_m2_k_w:0.0}],1,||Ok(())).unwrap();
    let gate = CancelGate::new_clock_free(); let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate,arena,fs_exec::StreamKey {seed:73,kernel_id:52,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        let prepared = assembly.prepare(&cx,&thermal).unwrap();
        let prepared = if s.integrated {prepared.with_integrated_conduction(&cx,&thermal,&integrated,
            LayerConductionConfig {max_iterations:128,absolute_heat_tolerance_w:1e-13,relative_heat_tolerance:1e-11}).unwrap()} else {prepared};
        run(&m,&plain,&powered,&bases,thermal,&prepared,&cx)
    })
}
fn flux(s:Spec,t:[f64;2])->(f64,f64) {
    let conductive = if s.integrated {
        // Common affine shape: integrate k(T) over BOTH original series layers.
        LINK_AREA/(0.05/0.04+0.05/0.08)*(t[0]-t[1])*(1.0+0.015*(t[0]+t[1]-600.0))
    } else { LINK_AREA/(0.05/0.04+0.05/0.08+if s.paths {s.resistance} else {0.0})*(t[0]-t[1]) };
    let radiation = if s.paths {FACE_SIDE*FACE_SIDE*SIGMA*(t[0]-t[1])*(t[0]+t[1])
        *(t[0]*t[0]+t[1]*t[1])/(1.0/s.emissivity+1.0/0.6-1.0)} else {0.0};
    (conductive,radiation)
}

#[test]
fn source_heat_transport_changes_charge_rates_and_motion_with_local_and_global_balances() {
    let mut ends = Vec::new();
    for s in [spec(),Spec {emissivity:0.8,..spec()},Spec {resistance:9.0,..spec()},
        Spec {integrated:true,..spec()},Spec {powered:false,..spec()},Spec {source_w:0.00002,..spec()}] {
        fixture(s,|m,plain,powered,bases,mut heat,network,cx| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
            if s.powered {owner.bind_material_circuit(&heat,powered,circuit(powered,&heat,bases,None),bases,&gate).unwrap();}
            let initial=owner.mechanical().total_energy_j().unwrap()+caloric(&heat)
                +if s.powered {INITIAL_CHARGE.powi(2)/(2.0*CAPACITANCE)} else {0.0};
            let mut force=placement.project_nodal_force(&[0.0,0.0,0.0,0.0,0.002,0.0],&gate).unwrap();
            force.resize(owner.mechanical().mode_count(),0.0);
            let (mut outside,mut work)=(0.0,0.0);
            for n in 1..=12 {
                let cells=heat.cells().clone();let q0=owner.electrical().map_or(0.0,|e|e.charges_c()[0]);
                let report=owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network,cx},
                    s.powered.then_some(input()),solve(),&[0.0;2],&gate).unwrap();
                let t=[heat.temperature_k(CellId(0)).unwrap(),heat.temperature_k(CellId(1)).unwrap()];
                let (conduction,radiation)=flux(s,t);let mut joule=[0.0;2];let mut stored=0.0;
                if let Some(e)=&report.electrical {
                    let r=[0.2*(1.0+0.02*(t[0]-300.0))/area(0),0.6*(1.0+0.02*(t[1]-300.0))/area(1)];
                    let a=DT/(2.0*CAPACITANCE*(r[0]+r[1]));let q=q0*(1.0-a)/(1.0+a);
                    assert!((owner.electrical().unwrap().charges_c()[0]-q).abs()<1e-12);
                    let current=(q0+q)/(2.0*CAPACITANCE*(r[0]+r[1]));
                    for i in 0..2 {joule[i]=DT*current*current*r[i];assert!((e.resistors[i].heat_j-joule[i]).abs()<1e-12);}
                    stored=e.stored_after_j;work+=e.native.supplied;
                }
                for i in 0..2 {
                    let id=CellId(i as u64);let c=&heat.cells()[&id];
                    let loss=report.layers.iter().find(|l|l.destination.cell==id).unwrap().heat_j;
                    let exchange=DT*(conduction+radiation)*if i==0 {-1.0} else {1.0};
                    let authored=DT*if i==0 {s.source_w} else {0.0};
                    assert!((c.mass_kg*(c.specific_h_j_kg-cells[&id].specific_h_j_kg)-loss-joule[i]-exchange-authored).abs()<3e-12);
                }
                let transfer=report.heat_transfer.as_ref().unwrap();
                assert_eq!(transfer.integrated_conduction.len(),usize::from(s.integrated));
                assert_eq!(transfer.thermal_contacts.len(),usize::from(!s.integrated));
                if s.integrated {assert!((transfer.integrated_conduction[0].applied_power_a_to_b_w-conduction).abs()<1e-12);}
                let rad=transfer.radiation.as_ref().unwrap();
                assert!((rad.sources_w()[&CellId(0)]+radiation).abs()<1e-11);
                assert!((report.radiation_exchange_entropy_j_k-DT*radiation*(1.0/t[1]-1.0/t[0])).abs()<1e-13);
                assert_eq!(report.entropy.has_unclassified_sources,s.source_w!=0.0);
                assert!((report.entropy.unclassified_source_at_cell_j_k-DT*s.source_w/t[0]).abs()<1e-13);
                assert!(report.entropy.accounting_residual_j_k.abs()<1e-12);
                assert!((report.external_thermal_heat_j-DT*s.source_w).abs()<1e-12);
                outside+=report.external_mechanical_loss_j+report.mechanical.numerical_dissipation_j;
                work+=report.mechanical.external_work_j+DT*s.source_w;
                near(owner.mechanical().total_energy_j().unwrap()+caloric(&heat)+stored-initial+outside-work,0.0);
                assert_eq!(owner.mechanical().samples_rendered(),n);assert_eq!(heat.time_s(),n as f64*DT);
            }
            ends.push((placement.physical_state(owner.mechanical().components(),&gate).unwrap()[0],
                heat.temperature_k(CellId(0)).unwrap(),owner.electrical().map_or(0.0,|e|e.charges_c()[0])));
        });
    }
    for i in [1,2] { // ONLY emissivity or thermal contact data changed.
        assert!((ends[i].0-ends[0].0).abs()>1e-12);
        assert!((ends[i].1-ends[0].1).abs()>1e-5);
        assert!((ends[i].2-ends[0].2).abs()>1e-10);
    }
}

#[test]
fn empty_additional_transport_preserves_bulk_path_for_powered_and_unpowered_histories() {
    for powered_flag in [false,true] {
        fixture(Spec {powered:powered_flag,paths:false,..spec()},|m,plain,powered,bases,heat,network,cx| {
            let gate=CancelGate::new_clock_free();let mut states=Vec::new();
            for full in [false,true] {
                let mut heat=heat.clone();let (_,material,native)=native(m);
                let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
                if powered_flag {owner.bind_material_circuit(&heat,powered,circuit(powered,&heat,bases,None),bases,&gate).unwrap();}
                let force=vec![0.0;owner.mechanical().mode_count()];
                for _ in 0..3 {
                    if full {owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network,cx},
                        powered_flag.then_some(input()),solve(),&[0.0;2],&gate).unwrap();}
                    else if powered_flag {owner.step_electrical(&mut heat,&force,input(),solve(),&[0.0;2],&gate).unwrap();}
                    else {owner.step(&mut heat,&force,solve(),&[0.0;2],&gate).unwrap();}
                }
                states.push((heat.cells().clone(),owner.mechanical().material_memory().to_vec(),
                    owner.mechanical().components().iter().map(|c|c.states().to_vec()).collect::<Vec<_>>(),owner.electrical().cloned()));
            }
            assert_eq!(states[0],states[1]);
        });
    }
}

#[test]
fn late_transport_domain_budget_and_independent_cancellation_preserve_every_history() {
    for narrow in [false,true] {
        fixture(Spec {upper:if narrow {300.000001} else {340.0},..spec()},|m,plain,powered,bases,mut heat,network,cx| {
            let gate=CancelGate::new_clock_free();let (placement,material,native)=native(m);
            let mut owner=FrameRelaxationThermalSystem::new(material,native,&heat,plain,&DEST,budget(),&gate).unwrap();
            owner.bind_material_circuit(&heat,powered,circuit(powered,&heat,bases,None),bases,&gate).unwrap();
            let original=owner.electrical().unwrap().clone();let cells=heat.cells().clone();
            let memory=owner.mechanical().material_memory().to_vec();let material=owner.material().members()[0].source.clone();
            let physical=placement.physical_state(owner.mechanical().components(),&gate).unwrap();
            let force=vec![0.0;owner.mechanical().mode_count()];
            let mut exhausted=solve();exhausted.maximum_material_trials=1;
            assert!(owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network,cx},Some(input()),exhausted,&[0.0;2],&gate).is_err());
            if narrow {assert!(owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network,cx},Some(input()),solve(),&[0.0;2],&gate).is_err());}
            let cancelled=CancelGate::new_clock_free();cancelled.request();
            let native_gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
            pool.scope(|arena| {
                let cx=Cx::new(&cancelled,arena,fs_exec::StreamKey {seed:1,kernel_id:2,tile:0,iteration:0},fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
                assert!(matches!(owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network,cx:&cx},Some(input()),solve(),&[0.0;2],&native_gate),Err(CircuitEnthalpyError::Cancelled)));
            });
            assert!(native_gate.is_requested());assert_eq!(owner.electrical().unwrap(),&original);
            assert_eq!(heat.cells(),&cells);assert_eq!(heat.time_s(),0.0);assert_eq!(owner.mechanical().samples_rendered(),0);
            assert_eq!(owner.mechanical().material_memory(),memory);assert_eq!(owner.material().members()[0].source,material);
            assert_eq!(placement.physical_state(owner.mechanical().components(),&gate).unwrap(),physical);
            if !narrow {owner.step_heat_network(&mut heat,&force,RelaxationHeatInput {network,cx},Some(input()),solve(),&[0.0;2],&gate).unwrap();}
        });
    }
}
