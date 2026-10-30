//! The original frame, circuit, interface and enclosure owners on one clock.
//! These are manufactured material cards and an explicit parallel-plate limit.
use super::*;
use fs_conduction::{GrayDiffuseEnclosure, RadiationSurface, SurfaceEmissivity, ViewFactorMatrix,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_couple::circuit_enthalpy::{contacts::{heat::*, material::*}, field_storage::FieldThermalBinding};
use fs_couple::frame_vibration::thermoelastic::{FrameHeatInput, FrameHeatError};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_material::{field_storage::{ResolvedCapacitor, ResolvedInductor, thermal::ThermalFieldStorage,
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
    MAGNETIC_PERMEABILITY_PROPERTY}, electrical_interface::{ELECTRICAL_CONTACT_RESISTANCE_DIMS as ED,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as EP}, thermal::interface::ThermalInterfacePoint};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext, ClaimSelection};
use fs_exec::Cx;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const OTHER: CellId = CellId(19);
const FRACTION: f64 = 0.25;
const FACE_AREA: f64 = 0.01;
const SIGMA: f64 = 5.670374419e-8;
#[derive(Clone, Copy)]
struct HeatSpec { magnetic: bool, powered: bool, paths: bool, cached: bool,
    emissivity: f64, rt: f64, upper: f64, fraction: f64 }
fn spec() -> HeatSpec { HeatSpec { magnetic: false, powered: true, paths: true, cached: true,
    emissivity: 0.2, rt: 0.3, upper: 600.0, fraction: FRACTION } }
fn mechanical_spec() -> Spec { Spec { alpha: 0.00002, e_slope: -1.0, ..super::spec() } }
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static)), t).unwrap().with("omega", 2048.0).unwrap()
}
fn source(i: usize, s: HeatSpec) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let (field, dims) = if s.magnetic { (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS) }
        else { (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS) };
    for (name, dims, value, slope) in [
        ("density", Density::DIMS, RHO, 0.0), ("young_modulus", Pressure::DIMS, E0, -1.0),
        ("eta", DynViscosity::DIMS, ETA, 0.0),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, 0.00002, 0.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, (i+1) as f64, 0.0),
        ("thermal_conductivity", fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS, 0.02, 0.0),
        ("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 5.0*area(), 0.01*area()),
        (field, dims, CAPACITANCE/area(), 0.001*CAPACITANCE/area()),
        (fs_conduction::SURFACE_EMISSIVITY_PROPERTY, fs_conduction::EMISSIVITY_DIMS,
            if i==0 { s.emissivity } else { 0.6 }, 0.0),
    ] {
        let upper = if i==0 && name==fs_conduction::SURFACE_EMISSIVITY_PROPERTY { s.upper } else { 600.0 };
        let mut validity = ValidityDomain::unconstrained().with_quantity("T",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static)), 200.0, upper);
        if name=="eta" { validity=validity.with("omega",0.0,4096.0); }
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(name,dims),
            value:if slope==0.0 { PropertyValue::Scalar { value,dims } } else {
                PropertyValue::Curve { abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                    knots:vec![(200.0,value-100.0*slope),(600.0,value+300.0*slope)],dims } },
            validity,uncertainty:UncertaintyModel::Unstated,
            interpolation:if slope==0.0 {InterpolationPolicy::ConstantWithinValidity} else {InterpolationPolicy::LinearInside},
            observations:vec![],provenance:Provenance {source:"manufactured full-heat frame".into(),license:"test-only".into(),artifact:None},
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:format!("body {i}"),phase:"solid".into(),
        process:"manufactured".into(),revision:0},claims,vec![]).unwrap()
}
fn solve() -> FrameThermalSolve {
    let mut c=controls();c.maximum_cells=2;c.maximum_material_trials=256;c.newton.linear_restart=2;c
}
fn fixture<R>(s:HeatSpec,run:impl FnOnce(&Cx<'_>,&CancelGate,FrameThermoelasticSystem<'_,'_>,
    DepositionSimulation<'_>,&PreparedHeatNetwork<'_,'_>)->R)->R {
    let cards=[source(0,s),source(1,s)];let context=point(300.0);
    let curves=cards.each_ref().map(|card|ResolvedHeatCapacityCurve::try_new(card,HeatCapacitySourceQuery {
        heat_capacity:&PropertyKey::new("specific_heat_capacity",SPECIFIC_HEAT_CAPACITY_DIMS),
        density:&PropertyKey::new("density",Density::DIMS),temperature_axis:"T",lower:&point(200.0),upper:&point(600.0),
        selection:SelectionPolicy::SingleClaimOnly},SolidLiquidPhase::Solid,0.0,1e-9).unwrap());
    let mut events:Vec<_>=(0..2).map(|i|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:[CELL,OTHER][i],segment_id:i as u64,body_id:i as u64,material_id:MaterialId(i as u64),mass_kg:RHO*area(),
        specific_h_j_kg:(if i==1 && !s.powered {130.0} else {100.0})*(i+1) as f64,
        geometry:CellGeometry {start_m:[0.0,i as f64,0.0],end_m:[1.0,i as f64,0.0],volume_m3:area(),
            surface_area_m2:2.0*area()+2.0*core::f64::consts::PI*RADIUS},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})}).collect();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {id:3,a:CELL,b:OTHER,area_m2:area(),conductance_w_k:999.0})});
    let thermal=DepositionSimulation::new(events,&[(MaterialId(0),curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1),curves[1].curve() as &dyn ThermalLaw)],DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:128}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:DT,energy_tolerance_j:1e-10},||Ok(())).unwrap().0;
    let key=PropertyKey::new("thermal_conductivity",fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS);
    let cells=[0,1].map(|i|MaterialCell {cell:[CELL,OTHER][i],material_id:MaterialId(i as u64),curve:curves[i].curve(),
        point:ThermalMaterialPoint {material:&cards[i],context:&context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let links=[ConductionLink {link:3,distance_a_m:0.1,distance_b_m:0.1}];
    let thermal_source=MaterialNetwork {cells:&cells,resistors:&[],links:&links};
    let resistors=[CellResistor {destination:ResistorCell {branch:0,cell:CELL},selection:MaterialPropertySelection::SingleClaimOnly,length_m:1.0,area_m2:area()}];
    let circuit_source=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let mut claims=ClaimSet::new();
    for (name,dims,value,slope) in [(RP,RD,s.rt,0.0),(EP,ED,2.0*area(),0.004*area())] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:if slope==0.0 {PropertyValue::Scalar {value,dims}} else {PropertyValue::Curve {abscissa:"Ta".into(),
                abscissa_dims:Temperature::DIMS,knots:vec![(200.0,value-100.0*slope),(600.0,value+300.0*slope)],dims}},
            validity:ValidityDomain::unconstrained().with("Ta",200.0,600.0).with("Tb",200.0,600.0),uncertainty:UncertaintyModel::Unstated,
            interpolation:if slope==0.0 {InterpolationPolicy::ConstantWithinValidity} else {InterpolationPolicy::LinearInside},observations:vec![],
            provenance:Provenance {source:"manufactured ordered heat interface".into(),license:"test-only".into(),artifact:None}}).unwrap();
    }
    let interface=InterfaceSystemCard::assemble(SurfaceSpec {material:cards[0].id().clone(),texture_frame:"A".into()},
        SurfaceSpec {material:cards[1].id().clone(),texture_frame:"B".into()},SystemContext {medium:"test".into(),third_body:None,
            environment:"test".into(),history:"fixed".into()},claims,vec![]).unwrap();
    let interface_context=QueryPoint::new().with("Ta",300.0).unwrap().with("Tb",300.0).unwrap();
    let interface_point=ThermalInterfacePoint {interface:&interface,materials:[&cards[0],&cards[1]],context:&interface_context,
        temperature_axes:["Ta","Tb"],liquid_fraction_axes:[None,None],mean_temperature_axis:None};
    let contacts=[MaterialElectricalContact {link:3,branch:0,point:interface_point,selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:FRACTION}];
    let heat_contacts=[MaterialElectricalContact {link:3,branch:0,point:interface_point,
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:s.fraction}];
    let bases=[(0,ElectricalPropertyBasis::Resistivity)];
    let electrical=ElectricalContactNetwork {bulk:if s.powered {&circuit_source} else {&thermal_source},
        contacts:if s.powered {&contacts} else {&[]},electrical_bases:if s.powered {&bases} else {&[]}};
    let heat_electrical=ElectricalContactNetwork {bulk:electrical.bulk,
        contacts:if s.powered {&heat_contacts} else {&[]},electrical_bases:electrical.electrical_bases};
    let thermal_contacts=[MaterialContactLink {link:3,point:interface_point,
        sides:[0,1].map(|i|ContactCell {cell:cells[i].cell,material_id:cells[i].material_id,curve:cells[i].curve}),
        selection:ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)}];
    let radiating=[0,1].map(|i|RadiatingCell {surface:i,cell:cells[i].cell,material_id:cells[i].material_id,
        curve:cells[i].curve,material:cells[i].point,selection:MaterialPropertySelection::SingleClaimOnly});
    let (complex,positions)=fs_conduction::fixtures::box_grid([1,1,1],[0.1;3]);
    let mesh=fs_conduction::mesh::ConductionMesh::new(complex,positions).unwrap();
    let surfaces=(0..2).map(|i| {
        let emissivity=SurfaceEmissivity::from_thermal_state("test",cells[i].point,
            cells[i].curve.state_at_specific_enthalpy(thermal.cells()[&cells[i].cell].specific_h_j_kg).unwrap(),
            MaterialPropertySelection::SingleClaimOnly).unwrap().0;
        RadiationSurface::new(&mesh,format!("face {i}"),|face|fs_conduction::fixtures::on_box_face(face.centroid[0],0.1*i as f64),emissivity).unwrap()
    }).collect();
    // Explicit infinite-plate limiting view factors; cube faces supply area only.
    let enclosure=GrayDiffuseEnclosure::new(surfaces,ViewFactorMatrix::infinite_parallel_plates(FACE_AREA).unwrap()).unwrap();
    let radiation=MaterialRadiation {enclosure:&enclosure,cells:&radiating,power_tolerance_w:1e-10,radiosity_tolerance_w_m2:1e-8};
    let assembly=MaterialHeatNetwork {electrical:&heat_electrical,thermal_contacts:if s.paths {&thermal_contacts} else {&[]},
        radiation:if s.paths {Some(&radiation)} else {None}};
    let gate=CancelGate::new_clock_free();
    let component=frame_component(&cards[0],&context,0.01);
    let mut owner=component.into_thermoelastic_system(&thermal,&thermal_source,&[FrameLayerCell {element:0,layer:0,cell:CELL}],0.05,
        ModalCouplingConfig {max_modes:1,max_connections:6,max_setup_terms:4096,nyquist_guard_fraction:0.9,
            maximum_total_energy_j:10.0,maximum_abs_pressure_pa:100.0,maximum_abs_connection_force_n:1000.0,
            solve_relative_tolerance:1e-11,energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-10},
        s.cached.then_some(PreparedModalBudget {max_modes:1,max_dense_bytes:8*1024*1024,max_multiply_terms:1_000_000_000}),&gate).unwrap();
    if s.powered {
        let prepared=electrical.prepare(&thermal,||false).unwrap();let initial=prepared.resolve(&thermal,||false).unwrap();
        let selection=MaterialPropertySelection::SingleClaimOnly;
        let (branch,element,law)=if s.magnetic {
            let value=ResolvedInductor::resolve(&cards[0],&context,selection,area(),1.0,1).unwrap();
            let law=ThermalFieldStorage::inductor(&cards[0],&value,"T",[200.0,600.0]).unwrap();
            (Branch::Inductor {henries:f64::NAN},MaterialStorageElement::Inductor(value),law)
        } else {
            let value=ResolvedCapacitor::resolve(&cards[0],&context,selection,area(),1.0).unwrap();
            let law=ThermalFieldStorage::capacitor(&cards[0],&value,"T",[200.0,600.0]).unwrap();
            (Branch::Capacitor {farads:f64::NAN},MaterialStorageElement::Capacitor(value),law)
        };
        let circuit=MaterialCircuit::new_with_storage(CircuitGraph {node_count:2,branches:vec![(1,0,Branch::Resistor {ohms:f64::NAN}),(1,0,branch)],
            transformers:vec![]},initial.bulk.resistors,vec![MaterialStorage {branch:1,element}],
            if s.magnetic {vec![CHARGE]} else {vec![]},if s.magnetic {vec![]} else {vec![CHARGE]},&[]).unwrap()
            .rebind_contacts_candidate(initial.contacts,&[]).unwrap();
        owner=owner.with_material_contact_field_circuit(&thermal,&electrical,circuit,
            vec![FieldThermalBinding {branch:1,cell:CELL,law,maximum_abs_coordinate:CHARGE}],&gate).unwrap();
    }
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:73,kernel_id:52,tile:0,iteration:0},fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        let prepared=assembly.prepare(&cx,&thermal).unwrap();run(&cx,&gate,owner,thermal,&prepared)
    })
}
fn coordinate(owner:&FrameThermoelasticSystem<'_,'_>,s:HeatSpec)->f64 {
    owner.electrical().map_or(0.0,|c|if s.magnetic {c.fluxes_wb()[0]} else {c.charges_c()[0]})
}
fn field_u(t:f64,z:f64)->f64 {let f=1.0+0.001*(t-300.0);0.5*z*z/(CAPACITANCE*f)*(1.0+0.001*t/f)}
fn radiative(ta:f64,tb:f64,s:HeatSpec)->f64 {
    if s.paths {FACE_AREA*SIGMA*(ta-tb)*(ta+tb)*(ta*ta+tb*tb)/(1.0/s.emissivity+1.0/0.6-1.0)} else {0.0}
}
fn bisect(mut a:f64,mut b:f64,f:impl Fn(f64)->f64)->f64 {
    assert!(f(a)<0.0 && f(b)>0.0);
    for _ in 0..80 {let m=0.5*(a+b);if f(m)>0.0 {b=m;} else {a=m;}}
    0.5*(a+b)
}
// Independent exact frozen-T oscillator, RC/RL midpoint and two material heat
// equations. No native report supplies the oracle's motion, loss or temperature.
fn reference(s:HeatSpec,t0:[f64;2],z0:f64,x0:[f64;2])->([f64;2],f64,[f64;2],[f64;3]) {
    let (x,loss)=motion(mechanical_spec(),t0[0],x0);
    let g=area()/(10.0+if s.paths {s.rt} else {0.0});
    let evaluate=|ta:f64| {
        let k=CAPACITANCE*(1.0+0.001*(ta-300.0));let r=7.0*(1.0+0.002*(ta-300.0));
        let ratio=if s.magnetic {DT*r/(2.0*k)} else {DT/(2.0*k*r)};
        let z=if s.powered {z0*(1.0-ratio)/(1.0+ratio)} else {0.0};
        let i=if s.powered {0.5*(z0+z)/k/if s.magnetic {1.0} else {r}} else {0.0};
        let bulk=DT*i*i*5.0*(1.0+0.002*(ta-300.0));let contact=0.4*bulk;
        let tb=bisect(200.0,600.0,|tb|2.0*capacity()*(tb-t0[1])-(1.0-FRACTION)*contact-DT*(g*(ta-tb)+radiative(ta,tb,s)));
        let balance=capacity()*(ta-t0[0])+2.0*capacity()*(tb-t0[1])+0.5*mass()*(x[1]*x[1]-x0[1]*x0[1])
            +potential(mechanical_spec(),ta,x[0]).1-potential(mechanical_spec(),t0[0],x0[0]).1+field_u(ta,z)-field_u(t0[0],z0);
        (balance,[ta,tb],z,[loss,bulk,contact])
    };
    let ta=bisect(200.0,550.0,|ta|evaluate(ta).0);let (_,t,z,h)=evaluate(ta);(t,z,x,h)
}

#[test]
fn full_heat_sources_drive_native_frames_and_field_storage_on_one_endpoint() {
    for magnetic in [false,true] {for cached in [false,true] {for powered in [false,true] {
        let s=HeatSpec {magnetic,cached,powered,..spec()};
        fixture(s,|cx,gate,mut owner,mut thermal,heat| {
            let initial=thermal.energy_j()+owner.mechanical_internal_energy_j().unwrap()+field_u(300.0,coordinate(&owner,s));
            for n in 1..=12 {
                let expected=reference(s,[thermal.temperature_k(CELL).unwrap(),thermal.temperature_k(OTHER).unwrap()],coordinate(&owner,s),state(&owner));
                let report=owner.step_heat_network(&mut thermal,&[0.0],FrameHeatInput {network:heat,cx},
                    powered.then_some(input(&[])),solve(),&[0.0;2],gate).unwrap();
                near(thermal.temperature_k(CELL).unwrap(),expected.0[0],2e-7);near(thermal.temperature_k(OTHER).unwrap(),expected.0[1],2e-7);
                near(coordinate(&owner,s),expected.1,2e-9);near(state(&owner)[0],expected.2[0],2e-10);near(state(&owner)[1],expected.2[1],2e-9);
                near(report.deposited_material_heat_j,expected.3[0],2e-10);near(report.deposited_joule_heat_j,expected.3[1]+expected.3[2],2e-9);
                let transfer=report.heat_transfer.as_ref().unwrap();let rad=transfer.radiation.as_ref().unwrap();
                near(transfer.thermal_contacts[0].conductance_w_k,area()/(10.0+s.rt),2e-12);
                near(rad.sources_w()[&CELL],-radiative(expected.0[0],expected.0[1],s),2e-8);
                near(report.entropy.radiation_exchange_entropy_j_k,
                    DT*radiative(expected.0[0],expected.0[1],s)*(1.0/expected.0[1]-1.0/expected.0[0]),2e-12);
                near(report.external_thermal_heat_j,0.0,2e-10);near(report.total_energy_residual_j,0.0,3e-10);
                near(report.entropy.accounting_residual_j_k,0.0,2e-11);assert!(!report.entropy.has_unclassified_sources);
                near(thermal.energy_j()+owner.mechanical_internal_energy_j().unwrap()+field_u(expected.0[0],expected.1),initial,4e-9);
                assert_eq!(owner.mechanical().samples_rendered(),n);near(thermal.time_s(),n as f64*DT,0.0);
            }
        });
    }}}
}

#[test]
fn changing_only_emissivity_or_thermal_contact_changes_later_motion_and_charge() {
    let mut endpoints=Vec::new();
    for (emissivity,rt) in [(0.2,0.3),(0.8,0.3),(0.2,30.0)] {
        let s=HeatSpec {emissivity,rt,..spec()};fixture(s,|cx,gate,mut owner,mut thermal,heat| {
            for _ in 0..32 {owner.step_heat_network(&mut thermal,&[0.0],FrameHeatInput {network:heat,cx},Some(input(&[])),solve(),&[0.0;2],gate).unwrap();}
            endpoints.push([thermal.temperature_k(CELL).unwrap(),thermal.temperature_k(OTHER).unwrap(),state(&owner)[0],coordinate(&owner,s)]);
        });
    }
    for changed in &endpoints[1..] {
        for (i,(a,b)) in endpoints[0].iter().zip(changed).enumerate() {
            assert!((a-b).abs()>if i==2 {1e-9} else {1e-7});
        }
    }
}

#[test]
fn absent_heat_paths_match_existing_arithmetic_and_refused_warm_steps_do_not_publish() {
    let run=|complete:bool,refuse:bool|fixture(HeatSpec {paths:false,..spec()},|cx,gate,mut owner,mut thermal,heat| {
        for _ in 0..3 {
            if complete {owner.step_heat_network(&mut thermal,&[0.0],FrameHeatInput {network:heat,cx},Some(input(&[])),solve(),&[0.0;2],gate).unwrap();}
            else {owner.step_electrical(&mut thermal,&[0.0],input(&[]),solve(),&[0.0;2],gate).unwrap();}
        }
        if refuse {
            let accepted=owner.electrical().unwrap().clone();let q=state(&owner);let cells=thermal.cells().clone();
            let last=owner.mechanical().last_frame().unwrap().clone();let mut budget=solve();budget.maximum_material_trials=1;
            assert!(owner.step_heat_network(&mut thermal,&[0.0],FrameHeatInput {network:heat,cx},Some(input(&[])),budget,&[0.0;2],gate).is_err());
            assert_eq!(owner.electrical().unwrap(),&accepted);assert_eq!(state(&owner),q);assert_eq!(thermal.cells(),&cells);
            assert_eq!(owner.mechanical().samples_rendered(),3);assert_eq!(owner.mechanical().last_frame(),Some(&last));
        }
        let report=if complete {owner.step_heat_network(&mut thermal,&[0.0],FrameHeatInput {network:heat,cx},Some(input(&[])),solve(),&[0.0;2],gate).unwrap()}
            else {owner.step_electrical(&mut thermal,&[0.0],input(&[]),solve(),&[0.0;2],gate).unwrap()};
        (state(&owner),owner.electrical().unwrap().clone(),thermal.cells().clone(),report.field_storage,report.entropy)
    });
    assert_eq!(run(true,false),run(false,false));assert_eq!(run(true,true),run(true,false));
}

#[test]
fn source_failure_heat_partition_mismatch_and_cancellation_keep_original_histories() {
    for (upper,fraction) in [(300.0001,FRACTION),(600.0,0.75)] {
        fixture(HeatSpec {upper,fraction,..spec()},|cx,gate,mut owner,mut thermal,heat| {
            let cells=thermal.cells().clone();let electrical=owner.electrical().unwrap().clone();let q=state(&owner);
            let result=owner.step_heat_network(&mut thermal,&[0.0],FrameHeatInput {network:heat,cx},Some(input(&[])),solve(),&[0.0;2],gate);
            if fraction==FRACTION {assert!(matches!(result,Err(FrameHeatError::Heat(HeatNetworkError::Radiation(_)))));}
            else {assert!(matches!(result,Err(FrameHeatError::Frame(CircuitEnthalpyError::Invalid("complete heat network changed original contact heat partitions")))));}
            assert_eq!(thermal.cells(),&cells);assert_eq!(owner.electrical().unwrap(),&electrical);assert_eq!(state(&owner),q);
            assert_eq!(owner.mechanical().samples_rendered(),0);assert_eq!(thermal.time_s(),0.0);
        });
    }
    fixture(spec(),|cx,gate,mut owner,mut thermal,heat| {
        let cells=thermal.cells().clone();let q=state(&owner);gate.request();
        assert!(matches!(owner.step_heat_network(&mut thermal,&[0.0],FrameHeatInput {network:heat,cx},Some(input(&[])),solve(),&[0.0;2],gate),
            Err(FrameHeatError::Frame(CircuitEnthalpyError::Cancelled))));
        assert_eq!(thermal.cells(),&cells);assert_eq!(state(&owner),q);assert_eq!(owner.mechanical().samples_rendered(),0);
    });
}
