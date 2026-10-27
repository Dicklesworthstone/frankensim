//! Complete material transport into distinct native acoustic chambers.
use super::*;
use std::collections::BTreeMap;
use fs_couple::bernoulli_aperture::wall::material::{
    FiniteGasDgBinding, FiniteGasDgNode, FiniteGasDgTransportInput,
    FiniteGasDgTransportSolve, FiniteGasDgTransportError, advance_finite_gas_dg_transport,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::{SnapshotNozzleLink, SnapshotNozzleError},
    species::{
        diffusion::{BinaryDiffusionMaterial, GasDiffusionGeometry,
            BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS,
            network::{DiffusionNetworkLink, DiffusionFaceSolve}},
        network::{SpeciesNetworkSolve, transport::{GasTransportSolve, GasTransportError}},
    },
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

fn field_for(gas: &ClosedGasParcel) -> AxialGasField {
    AxialGasField::new(vec![0.0,LENGTH], gas.state().state().pressure,
        vec![PrimitiveGasState {
            temperature_k: gas.state().state().temperature,
            composition: gas.model().thermodynamics().mole_composition().clone(),
        };2], gas.model().clone(),
        QueryPoint::new().with_quantity("material-lot", QuantitySpec::dimensional(Dims([0;6])),42.0).unwrap(),
        ThermochemicalPhase::FrozenGas,
    ).unwrap()
}
fn fixture() -> (BTreeMap<u64,ClosedGasParcel>,BTreeMap<u64,FiniteGasDgNode>) {
    let gases=BTreeMap::from([
        (0,parcel(0.8,400.0,101325.0)), (1,parcel(0.2,400.0,101325.0)),
        (2,parcel(0.2,400.0,70000.0)),
    ]);
    let mut nodes=BTreeMap::new();
    for id in [0,1] {
        let duct=acoustic(&field_for(&gases[&id]),true);
        let binding=if id==0 {
            FiniteGasDgBinding::Uniform(kelvin_wall(1400.0,2e8,2.0,2e-5).bind_dg(duct).unwrap())
        } else {
            let source=source_wall(1800.0,4e8);
            let right=LocalMaterialWall::new(source.material().clone(),2e-5,
                source.axial_constraint(),source.loss().clone(),source.band_hz()).unwrap();
            FiniteGasDgBinding::Regional(RegionalMaterialWallDg::new(
                vec![kelvin_wall(1000.0,1e8,2.0,2e-5),right],duct,
            ).unwrap())
        };
        let propagator=binding.duct().propagator(1e-6,1e-13).unwrap();
        let mut state=binding.duct().project_initial_state(|x|1000.0+300.0*x,|_|2e-6);
        for _ in 0..8 { propagator.step(&mut state,&[]).unwrap(); }
        nodes.insert(id,FiniteGasDgNode {binding,state,propagator});
    }
    (gases,nodes)
}
fn inputs() -> BTreeMap<u64,FiniteGasDgTransportInput> {
    [0,1,2].into_iter().map(|id|(id,FiniteGasDgTransportInput {
        gas:energy_solve(), supplied_heat_j:if id==0 {1e-5} else if id==2 {-2e-6} else {0.0},
        acoustic:(id!=2).then_some(FiniteGasDgSolve {
            work_bracket_j:[-1e-7,1e-7],..controls(true)
        }),
        ports:if id==0 {vec![PortInput {id:PortId(7),effort:1e-8}]} else {vec![]},
    })).collect()
}
fn solve() -> FiniteGasDgTransportSolve {
    FiniteGasDgTransportSolve {
        transport:GasTransportSolve {
            exchange:SpeciesNetworkSolve {maximum_species_fraction:0.05,
                mass_tolerance_kg:1e-16,energy_tolerance_j:1e-11},
            entropy_tolerance_j_k:1e-12,
        },
        energy_tolerance_j:5e-11,
    }
}
fn nozzle() -> SnapshotNozzleLink {
    SnapshotNozzleLink {id:9,parcels:[0,2],
        nozzle:GasNozzle {area_m2:1e-6,discharge_coefficient:0.8},
        nozzle_solve:energy_solve()}
}
fn diffusivity(value: f64) -> MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY,BINARY_DIFFUSIVITY_DIMS),
        value:PropertyValue::Scalar {value,dims:BINARY_DIFFUSIVITY_DIMS},
        validity:ValidityDomain::unconstrained()
            .with_quantity("temperature",QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,ValueForm::Static)),200.0,900.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance {
            source:"synthetic binary diffusion coefficient; not measured gas validation".into(),
            license:"CC0-1.0".into(),artifact:None,
        },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:"A/B pair".into(),phase:"gas".into(),
        process:"manufactured network fixture".into(),revision:0},claims,vec![]).unwrap()
}
fn pair<'a>(card:&'a MaterialCard,point:&'a QueryPoint) -> [BinaryDiffusionMaterial<'a>;1] {
    [BinaryDiffusionMaterial {species:[SpeciesId::new("A").unwrap(),SpeciesId::new("B").unwrap()],
        material:card,context:point,selection:MaterialPropertySelection::SingleClaimOnly}]
}
fn face<'a>(pairs:&'a [BinaryDiffusionMaterial<'a>]) -> DiffusionNetworkLink<'a> {
    DiffusionNetworkLink {id:4,parcels:[0,1],pairs,
        geometry:GasDiffusionGeometry {area_m2:2e-6,distance_m:1e-5},
        solve:DiffusionFaceSolve {maximum_temperature_difference_k:1.0,
            maximum_pressure_difference_pa:200.0,gradient_tolerance_per_m:1e-7,
            molar_flux_tolerance_mol_m2_s:1e-10}}
}
fn species_mass(gases:&BTreeMap<u64,ClosedGasParcel>,name:&str) -> f64 {
    let id=SpeciesId::new(name).unwrap();
    gases.values().map(|gas| {
        let law=gas.model().thermodynamics();
        let masses:Vec<_>=law.components().iter().map(|c|c.molar_mass()).collect();
        let fractions=law.mole_composition().to_mass_fractions(&masses).unwrap();
        law.components().iter().zip(fractions.fractions())
            .filter(|(c,_)|c.species()==&id).map(|(_,f)|gas.mass_kg()*f).sum::<f64>()
    }).sum()
}
fn unchanged(
    nodes:&BTreeMap<u64,FiniteGasDgNode>, old:&BTreeMap<u64,FiniteGasDgNode>,
) {
    for (&id,node) in nodes {
        assert_eq!(node.state,old[&id].state);
        assert_eq!(node.binding.duct().gas_field(),old[&id].binding.duct().gas_field());
        assert_eq!(node.binding.duct().operator_identity(),old[&id].binding.duct().operator_identity());
        let (mut a,mut b)=(node.state.clone(),old[&id].state.clone());
        assert_eq!(node.propagator.step(&mut a,&[]).unwrap(),old[&id].propagator.step(&mut b,&[]).unwrap());
    }
}

#[test]
fn source_diffusion_and_nozzle_share_inventory_then_drive_two_actual_continuing_ducts() {
    let mut histories=Vec::new();
    for coefficient in [2e-5,4e-5] {
        let (mut gases,mut nodes)=fixture();
        let (old_gas,old_nodes)=(gases.clone(),nodes.clone());
        let card=diffusivity(coefficient);
        let context=QueryPoint::new();
        let pairs=pair(&card,&context);
        let face=face(&pairs);
        let input=inputs();
        let original_nozzle=gases[&0].nozzle_flow_to(
            gases[&2].state().state().pressure,nozzle().nozzle,energy_solve(),||false,
        ).unwrap();
        let report=advance_finite_gas_dg_transport(
            &mut gases,&mut nodes,&input,&[nozzle()],&[face.clone()],1e-6,solve(),&mut ||false,
        ).unwrap();
        assert_eq!(report.transport.nozzles[0].flow,original_nozzle);
        assert_eq!(report.transport.nozzles[0].transferred_mass_kg,original_nozzle.mass_flow_kg_s*1e-6);
        // Independent binary Maxwell-Stefan limit, in the original molar frame.
        let expected=101325.0/(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0)
            *coefficient*(0.8-0.2)/face.geometry.distance_m;
        let flux=&report.transport.diffusion[0].flux;
        close(flux.species[0].molar_flux_mol_m2_s,expected);
        close(flux.species[1].molar_flux_mol_m2_s,-expected);
        assert_eq!(flux.pairs[0].diffusivity_m2_s,coefficient);
        assert_eq!(flux.pairs[0].material.card_identity(),card.content_hash());
        for name in ["A","B"] {
            assert!((species_mass(&gases,name)-species_mass(&old_gas,name)).abs()<1e-16);
        }
        let gas_change:f64=gases.iter().map(|(id,g)|energy_change(&old_gas[id],g)).sum();
        let wave_change:f64=nodes.iter().map(|(id,n)|n.state.stored_energy_j()-old_nodes[id].state.stored_energy_j()).sum();
        let total_heat:f64=input.values().map(|i|i.supplied_heat_j).sum();
        let residual=gas_change+wave_change+report.acoustic_dissipation_j-total_heat-report.port_work_j;
        assert!(residual.abs()<=solve().energy_tolerance_j,"{residual}");
        assert_eq!(report.supplied_heat_j,total_heat);
        assert!(report.combined_energy_residual_j.abs()<=solve().energy_tolerance_j);
        assert_eq!(report.reservoirs[&2].supplied_heat_j,-2e-6);
        for (&id,node) in &nodes {
            let mut marched=old_nodes[&id].state.clone();
            old_nodes[&id].propagator.step(&mut marched,&input[&id].ports).unwrap();
            let previous=old_nodes[&id].binding.duct().compliant_wall_state(&marched).unwrap();
            same_motion(&node.binding.duct().compliant_wall_state(&node.state).unwrap(),&previous);
            if id==1 {assert!(previous.relaxation_stored_energy_j>0.0);}
            assert_eq!(node.state.step_count(),9);
            assert_eq!(node.state.initial_energy_j(),old_nodes[&id].state.initial_energy_j());
            assert_eq!(node.binding.duct().gas_field().unwrap().uniform_state(),Some(gases[&id].state()));
            assert_eq!(gases[&id].state().source_packs().len(),2);
            assert_ne!(gases[&id].state().state().sound_speed,old_gas[&id].state().state().sound_speed);
            assert_ne!(node.binding.duct().viscothermal_model().unwrap().series_terms()[0].rate_per_s,
                old_nodes[&id].binding.duct().viscothermal_model().unwrap().series_terms()[0].rate_per_s);
        }
        // A second interval samples the new inventories once and continues both clocks.
        advance_finite_gas_dg_transport(
            &mut gases,&mut nodes,&input,&[nozzle()],&[face],1e-6,solve(),&mut ||false,
        ).unwrap();
        let mut trace=Vec::new();
        for node in nodes.values_mut() {
            assert_eq!(node.state.step_count(),10);
            let parameter_work=node.state.cumulative_parameter_work_j();
            for _ in 0..8 {
                let frame=node.propagator.step(&mut node.state,&[]).unwrap();
                node.propagator.sync_nodal(&mut node.state).unwrap();
                trace.extend_from_slice(node.state.pressure());
                assert_eq!(frame.ledger.cumulative_parameter_work_j,parameter_work);
            }
        }
        histories.push(trace);
    }
    let difference:f64=histories[0].iter().zip(&histories[1]).map(|(a,b)|(a-b).powi(2)).sum();
    let signal:f64=histories[0].iter().map(|x|x*x).sum();
    assert!(difference>1e-18*signal);
}

#[test]
fn final_reservoir_heating_reversal_refuses_after_local_acoustic_work_closes() {
    let (mut gases,mut nodes)=fixture();
    gases.insert(2,parcel(0.2,400.0,100000.0));
    let (old_gas,old_nodes)=(gases.clone(),nodes.clone());
    let mut input=inputs();
    input.get_mut(&2).unwrap().supplied_heat_j=0.02;
    let error=advance_finite_gas_dg_transport(
        &mut gases,&mut nodes,&input,&[nozzle()],&[],1e-6,solve(),&mut ||false,
    ).unwrap_err();
    assert!(matches!(error,FiniteGasDgTransportError::Transport(
        GasTransportError::Nozzle(SnapshotNozzleError::EquilibriumCrossed {id:9})
    )),"{error:?}");
    assert_eq!(gases,old_gas);
    unchanged(&nodes,&old_nodes);
}

#[test]
fn second_duct_refusal_and_last_poll_cancellation_preserve_every_native_owner() {
    let (mut gases,mut nodes)=fixture();
    let (old_gas,old_nodes)=(gases.clone(),nodes.clone());
    let card=diffusivity(2e-5);
    let context=QueryPoint::new();
    let pairs=pair(&card,&context);
    let face=face(&pairs);
    let mut input=inputs();
    input.get_mut(&1).unwrap().acoustic.as_mut().unwrap().work_bracket_j=[0.0,0.0];
    assert!(matches!(advance_finite_gas_dg_transport(
        &mut gases,&mut nodes,&input,&[nozzle()],&[face.clone()],1e-6,solve(),&mut ||false,
    ),Err(FiniteGasDgTransportError::Acoustic {parcel:1,..})));
    assert_eq!(gases,old_gas);
    unchanged(&nodes,&old_nodes);
    let input=inputs();
    let mut polls=0;
    advance_finite_gas_dg_transport(
        &mut old_gas.clone(),&mut old_nodes.clone(),&input,&[nozzle()],&[face.clone()],1e-6,
        solve(),&mut ||{polls+=1;false},
    ).unwrap();
    let mut count=0;
    assert!(matches!(advance_finite_gas_dg_transport(
        &mut gases,&mut nodes,&input,&[nozzle()],&[face],1e-6,solve(),
        &mut ||{count+=1;count==polls},
    ),Err(FiniteGasDgTransportError::Cancelled)));
    assert_eq!(gases,old_gas);
    unchanged(&nodes,&old_nodes);
}

