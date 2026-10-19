//! Sourced liquid viscosity -> joint mechanical force -> finite material heat.
//! Synthetic constitutive data test wiring and balances, not measured fluids.
#[path = "material_film_heat/feedback.rs"]
mod feedback;
use fs_blake3::ContentHash;
use fs_couple::render::plate::impact::{ImpactBody, ImpactConfig, ImpactSystem};
use fs_couple::render::plate::impact::material::heating::{
    CaloricReservoir, HeatedImpactSystem, SqueezeFilmHeatRoute,
};
use fs_couple::render::plate::impact::material::heating::network::ReservoirCell;
use fs_couple::render::plate::impact::squeeze::{FilmCell, FilmChannel, FilmLimits, GapPort};
use fs_couple::render::plate::impact::squeeze::material::MaterialSqueezeFilm;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::{GasSpec, GasState};
use fs_material::liquid::{ResolvedLiquidState, resolve_liquid_state};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EnthalpyPhaseKnot, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, QuantitySpec};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};

const AREA: f64 = 1e-4;
const GAP: f64 = 2e-4;
const WIDTH: f64 = 0.01;
const LENGTH: f64 = 0.02;
const DT: f64 = 1e-6;
const MASS: f64 = 0.001;

fn liquid(mu: f64, cp: f64) -> ResolvedLiquidState {
    let mut claims=ClaimSet::new();
    for (name,value,dims,quantity) in [
        ("density",1000.0,Dims([-3,1,0,0,0,0]),QuantitySpec::dimensional(Dims([-3,1,0,0,0,0]))),
        ("specific_heat_capacity",cp,Dims([2,0,-2,-1,0,0]),QuantitySpec::dimensional(Dims([2,0,-2,-1,0,0]))),
        ("dynamic_viscosity",mu,Dims([-1,1,-1,0,0,0]),QuantitySpec::dimensional(Dims([-1,1,-1,0,0,0]))),
        ("thermal_conductivity",0.6,Dims([1,1,-3,-1,0,0]),QuantitySpec::semantic(
            SemanticType::new(QuantityKind::ThermalConductivity,ValueForm::Static))),
    ] {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::with_quantity(name,quantity),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",290.0,310.0).with("p",100000.0,100000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic film heat fixture".into(),
                license:"test-only".into(),artifact:None},
        }).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic liquid".into(),
        phase:"liquid".into(),process:"numerical fixture".into(),revision:0},claims,vec![]).unwrap();
    resolve_liquid_state(&card,&QueryPoint::new().with("T",300.0).unwrap().with("p",100000.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn film<F: FluidTransport + ?Sized>(source:&F) -> Result<MaterialSqueezeFilm<'_,F>,
    fs_couple::render::plate::impact::ImpactError> {
    let gap=GapPort {reference_m:GAP,closure:vec![1.0/0.2_f64.sqrt(),-1.0/0.4_f64.sqrt()]};
    MaterialSqueezeFilm::new(source,vec![FilmCell {area_m2:AREA,gap:gap.clone()}],
        vec![FilmChannel {from:0,to:None,width_m:WIDTH,length_m:LENGTH,gap}],2,
        FilmLimits {minimum_cell_gap_m:1e-5,maximum_gap_m:5e-4,maximum_pressure_pa:10000.0})
}
fn mechanics(source:&impl FluidTransport, extra_damping:f64) -> ImpactSystem {
    let (mut a,_)=ImpactBody::free_mass(0.2,0.0,0.001).unwrap();
    let (b,_)=ImpactBody::free_mass(0.4,0.0,0.0).unwrap();
    a.damping_per_s[0]=extra_damping;
    let system=ImpactSystem::new(vec![a,b],vec![],vec![],vec![],ImpactConfig {
        dt_s:DT,max_steps:512,maximum_energy_j:1.0,energy_absolute_tolerance_j:1e-14,
        energy_relative_tolerance:1e-7,maximum_generalized_force:1.0,
    }).unwrap();
    let (system,index)=film(source).unwrap().attach(system).unwrap();
    assert_eq!(index,0);system
}
fn curve(id:ContentHash,cp:f64,upper:f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(id,SolidLiquidPhase::Liquid,vec![
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:300.0,
            liquid_mass_fraction:1.0,bulk_density_kg_m3:1000.0},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:upper,temperature_k:300.0+upper/cp,
            liquid_mass_fraction:1.0,bulk_density_kg_m3:1000.0},
    ]).unwrap()
}
fn stores(c:&EquilibriumEnthalpyPhaseCurve) -> Vec<CaloricReservoir<'_>> {
    ["fluid","sink"].into_iter().map(|key|CaloricReservoir::new(key.into(),c,
        c.state_at_specific_enthalpy(0.0).unwrap(),MASS,1e-6,0.0006).unwrap()).collect()
}
fn route() -> SqueezeFilmHeatRoute {
    SqueezeFilmHeatRoute {film:0,to_a:"fluid".into(),to_b:"sink".into(),fraction_to_a:1.0}
}
fn heated<'a>(source:&impl FluidTransport,c:&'a EquilibriumEnthalpyPhaseCurve,
    extra_damping:f64) -> HeatedImpactSystem<'a> {
    HeatedImpactSystem::new_with_films(mechanics(source,extra_damping).prepare_analytic().unwrap(),
        stores(c),vec![],vec![route()],1e-13).unwrap()
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a} != {b}");}

#[test]
fn source_card_viscosity_drives_analytic_pressure_force_work_and_retains_the_query() {
    let a=liquid(0.01,0.2);let b=liquid(0.02,0.2);
    let velocity=[0.001*0.2_f64.sqrt(),0.0];
    for source in [&a,&b] {
        let binding=film(source).unwrap();
        assert!(std::ptr::eq(binding.source(),source));
        assert_eq!(binding.source().material().query_point(),source.material().query_point());
        assert_eq!(binding.transport().dynamic_viscosity_pa_s(),source.dynamic_viscosity_pa_s());
        let mut force=[0.0;2];let mut pressure=[0.0;1];
        let report=binding.film().evaluate_into(&[0.0;2],&velocity,&mut force,&mut pressure).unwrap();
        // Independent single-channel Poiseuille balance: Q=A*u, p=Q/G,
        // physical force=A*p; virtual work uses the mass-normalized port.
        let conductance=WIDTH*GAP.powi(3)/(12.0*source.dynamic_viscosity_pa_s()*LENGTH);
        let exact_pressure=AREA*0.001/conductance;
        near(pressure[0],exact_pressure,1e-10);
        near(force[0]*0.2_f64.sqrt(),AREA*exact_pressure,1e-14);
        near(force[1]*0.4_f64.sqrt(),-AREA*exact_pressure,1e-14);
        near(report.dissipated_power_w,AREA*0.001*exact_pressure,1e-16);
    }
    let mut sa=mechanics(&a,0.0).prepare_analytic().unwrap();
    let mut sb=mechanics(&b,0.0).prepare_analytic().unwrap();
    for _ in 0..64 {
        sa.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
        sb.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
    }
    assert!(sb.state()[3]>sa.state()[3]);
    assert!(sb.stored_energy_j()<sa.stored_energy_j());
}

#[test]
fn actual_film_work_heats_material_once_and_preserves_motion_in_both_newton_images() {
    let source=liquid(0.01,0.2);
    let c=curve(source.material().card_identity(),source.specific_heat_j_kg_k(),1.0);
    let gate=CancelGate::new_clock_free();
    for analytic in [false,true] {
        let mut plain=mechanics(&source,0.0).prepare().unwrap();plain.set_analytic_newton(analytic);
        let mut coupled=heated(&source,&c,0.0);coupled.set_analytic_newton(analytic);
        let initial=plain.stored_energy_j();let mut total=0.0;
        for _ in 0..128 {
            let before=plain.state().to_vec();let expected=plain.step(&[0.0;2],&gate).unwrap();
            let frame=coupled.step(&[0.0;2],&gate).unwrap();
            assert_eq!(plain.state(),coupled.mechanics().state());
            assert_eq!(expected.dissipated_energy_j,frame.mechanical.dissipated_energy_j);
            near(frame.deposited_heat_j,frame.mechanical.dissipated_energy_j,1e-20);
            // These free masses have quadratic storage, so their exact effort
            // is midpoint velocity. Verify heat independently from the graph.
            let after=plain.state();
            let closure=0.5*((before[0]+after[0])/0.2_f64.sqrt()
                -(before[2]+after[2])/0.4_f64.sqrt());
            let relative=0.5*((before[1]+after[1])/0.2_f64.sqrt()
                -(before[3]+after[3])/0.4_f64.sqrt());
            let conductance=WIDTH*(GAP-closure).powi(3)/(12.0*source.dynamic_viscosity_pa_s()*LENGTH);
            let independent_heat=DT*(AREA*relative).powi(2)/conductance;
            near(frame.deposited_heat_j,independent_heat,1e-20);
            total+=frame.deposited_heat_j;
        }
        assert!(total>1e-9);
        let fluid=&coupled.reservoirs()[0];
        near(fluid.received_heat_j(),total,1e-20);
        near(fluid.state().temperature_k()-300.0,total/(MASS*source.specific_heat_j_kg_k()),1e-12);
        assert_eq!(coupled.reservoirs()[1].received_heat_j(),0.0);
        near(coupled.mechanics().stored_energy_j()+fluid.storage_change_j(),initial,1e-13);
    }
}

#[test]
fn unrelated_damping_stays_external_and_heat_capacity_changes_temperature_only() {
    let a=liquid(0.01,0.2);let b=liquid(0.01,0.4);
    let ca=curve(a.material().card_identity(),a.specific_heat_j_kg_k(),1.0);
    let cb=curve(b.material().card_identity(),b.specific_heat_j_kg_k(),1.0);
    let mut sa=heated(&a,&ca,100.0);let mut sb=heated(&b,&cb,100.0);
    let initial=sa.mechanics().stored_energy_j();let mut outside=0.0;
    for _ in 0..128 {
        let fa=sa.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
        let fb=sb.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
        assert_eq!(fa.deposited_heat_j,fb.deposited_heat_j);
        assert_eq!(sa.mechanics().state(),sb.mechanics().state());
        outside+=fa.remaining_mechanical_loss_j;
    }
    assert!(outside>0.0);
    near(sa.reservoirs()[0].state().temperature_k()-300.0,
        2.0*(sb.reservoirs()[0].state().temperature_k()-300.0),1e-12);
    near(sa.mechanics().stored_energy_j()+sa.reservoirs()[0].storage_change_j()+outside,initial,1e-13);
}

fn network(c:&EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events=Vec::new();
    for id in [1,2] {
        events.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id:CellId(id),segment_id:id,body_id:id,material_id:MaterialId(1),mass_kg:MASS,
            specific_h_j_kg:0.0,geometry:CellGeometry {start_m:[id as f64,0.0,0.0],
                end_m:[id as f64+0.01,0.0,0.0],volume_m3:1e-6,surface_area_m2:0.0006},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0},
        })});
    }
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:7,a:CellId(1),b:CellId(2),area_m2:1e-6,conductance_w_k:1.0,
    })});
    DepositionSimulation::new(events,&[(MaterialId(1),c as &dyn ThermalLaw)],
        DepositionBudget {max_cells:2,max_links:1,max_events:8,max_substeps:16}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:DT/4.0,energy_tolerance_j:1e-13},||Ok(())).unwrap().0
}

#[test]
fn deposited_film_loss_and_conduction_share_the_existing_network_substeps() {
    let source=liquid(0.01,0.2);
    let c=curve(source.material().card_identity(),source.specific_heat_j_kg_k(),1.0);
    let mut system=heated(&source,&c,0.0);let mut thermal=network(&c);let mut observer=();
    let initial=system.mechanics().stored_energy_j();let mut substeps=0;
    let bindings=[ReservoirCell {reservoir:"fluid".into(),cell:CellId(1)},
        ReservoirCell {reservoir:"sink".into(),cell:CellId(2)}];
    for _ in 0..128 {
        let frame=system.step_with_thermal_network(&mut thermal,&mut observer,&bindings,
            StepConfig {max_dt_s:DT/4.0,energy_tolerance_j:1e-13},&[0.0;2],
            &CancelGate::new_clock_free(),|_,_|Ok(None)).unwrap();
        substeps+=frame.network.substeps;
        near(frame.external_heat_j,0.0,1e-19);
        near(frame.impact.thermal.transfers.iter().map(|x|x.report.deposited_heat_j).sum(),
            frame.impact.thermal.impact.deposited_heat_j,1e-19);
        for (i,binding) in bindings.iter().enumerate() {
            assert_eq!(system.reservoirs()[i].state().specific_enthalpy_j_kg(),
                thermal.cells()[&binding.cell].specific_h_j_kg);
        }
        assert!(thermal.cells().values().all(|cell|cell.boundary.source_w==0.0));
    }
    assert!(substeps>=512);
    assert!(system.reservoirs()[1].state().temperature_k()>300.0);
    assert_eq!(system.reservoirs()[1].received_heat_j(),0.0);
    near(system.mechanics().stored_energy_j()+system.reservoirs().iter()
        .map(|store|store.storage_change_j()).sum::<f64>(),initial,1e-13);
}

#[test]
fn thermal_domain_cancellation_invalid_routes_and_incoherent_gas_do_not_publish() {
    let source=liquid(0.01,0.2);
    let small=curve(source.material().card_identity(),0.2,1e-12);
    let mut system=heated(&source,&small,0.0);
    let motion=system.mechanics().state().to_vec();let stores=system.reservoirs().to_vec();
    assert!(system.step(&[0.0;2],&CancelGate::new_clock_free()).is_err());
    assert_eq!(system.mechanics().state(),motion);assert_eq!(system.reservoirs(),stores);
    assert_eq!(system.mechanics().samples(),0);
    let stopped=CancelGate::new_clock_free();stopped.request();
    assert!(system.step(&[0.0;2],&stopped).is_err());
    assert_eq!(system.mechanics().state(),motion);assert_eq!(system.reservoirs(),stores);
    for routes in [vec![route(),route()],vec![SqueezeFilmHeatRoute {film:1,..route()}],
        vec![SqueezeFilmHeatRoute {fraction_to_a:f64::NAN,..route()}]] {
        assert!(HeatedImpactSystem::new_with_films(mechanics(&source,0.0).prepare().unwrap(),
            stores.clone(),vec![],routes,1e-13).is_err());
    }
    let mut gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),300.0,100000.0).unwrap();
    assert!(film(&gas).is_ok());
    gas.prandtl*=2.0;
    assert!(film(&gas).is_err());
}
