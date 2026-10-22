//! Same source cards drive current, fluid momentum and total-h transport.
use super::*;
use std::cell::{Cell, RefCell};
use std::rc::Rc;
use fs_lbm::conjugate::enthalpy::advection::feedback::electrical::{
    MaterialFlowElectricalError, advance_joule_enthalpy_with_material_flow,
};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::{ChtError, ElectricalContact, ElectricSetup, Electrode, ElectrodeDrive};
use fs_material::MaterialError;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_material::thermal::source::ThermalMaterialPoint;

fn source_card(mu: f64, rho_e: f64, k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, variable) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0, false),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 2.0, false),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu, true),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, k, true),
        (ElectricalPropertyBasis::Resistivity.property_name(),
         ElectricalPropertyBasis::Resistivity.dimensions(), rho_e, true),
    ] {
        let mut claim = scalar(name, dims, value);
        claim.validity = claim.validity.with("liquid", 0.0, 1.0);
        if variable {
            claim.value = PropertyValue::Curve {
                abscissa: "liquid".into(), abscissa_dims: Dims([0; 6]),
                knots: vec![(0.0, value), (1.0, 2.0*value)], dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

#[derive(Clone, Debug, PartialEq)]
struct MeltSource {
    phase: EquilibriumPhaseState,
    conductor: ResolvedConductorState,
    transport: FluidTransportState,
    query: usize,
}
impl FluidTransport for MeltSource {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> { Ok(self.transport) }
}
fn source_error(error: impl core::fmt::Display) -> MaterialError {
    MaterialError::State { what: error.to_string() }
}
fn resolve(
    card: &MaterialCard, curve: &EquilibriumEnthalpyPhaseCurve,
    cell: &DepositedCell, temperature: f64, query: usize,
) -> Result<MeltSource, MaterialError> {
    let phase = curve.state_at_specific_enthalpy(cell.specific_h_j_kg).map_err(source_error)?;
    assert_eq!(phase.temperature_k().to_bits(), temperature.to_bits());
    let context = point(500.0, 100000.0).with("liquid", 0.0).unwrap();
    let binding = ThermalMaterialPoint {
        material: card, context: &context, temperature_axis: "T", liquid_fraction_axis: Some("liquid"),
    };
    let requirements = [
        ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_new(DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_with_quantity(THERMAL_CONDUCTIVITY_PROPERTY,
            QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static)),
            ScalarAdmissibility::StrictlyPositive).unwrap(),
        ElectricalPropertyBasis::Resistivity.requirement().unwrap(),
    ];
    let material = resolve_material_state_point(card, &binding.at(phase).map_err(source_error)?,
        &requirements, MaterialPropertySelection::SingleClaimOnly).map_err(source_error)?;
    for property in material.properties() { binding.admit_static(property).map_err(source_error)?; }
    let value = |name| material.property(name).unwrap().value_si();
    let transport = FluidTransportState::try_new(value(DENSITY_PROPERTY),
        value(SPECIFIC_HEAT_CAPACITY_PROPERTY), value(THERMAL_CONDUCTIVITY_PROPERTY),
        value(DYNAMIC_VISCOSITY_PROPERTY))?;
    let conductor = ResolvedConductorState::try_from_material_with_basis(material,
        ElectricalPropertyBasis::Resistivity).map_err(source_error)?;
    Ok(MeltSource { phase, conductor, transport, query })
}
fn setup(voltage: f64) -> ElectricSetup {
    ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo: [0,0,0], hi: [1,1,1], drive: ElectrodeDrive::Voltage(voltage) },
        Electrode { lo: [3,0,0], hi: [4,1,1], drive: ElectrodeDrive::Voltage(0.0) },
    ] }
}
fn electric_input<'a>(
    domain: &'a VoxelDomain, setup: &'a ElectricSetup, ids: &'a [CellId],
    contacts: &'a [ElectricalContact],
) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput { domain, setup, contacts, cell_ids: ids, max_cells: 4,
        max_current_solves: 100, power_tolerance_w: 1e-8, energy_tolerance_j: 1e-8 }
}
fn resistivity(sample: &fs_lbm::conjugate::enthalpy::advection::feedback::FlowMaterialSample<MeltSource>)
    -> Result<Option<f64>, ChtError> {
    Ok(Some(sample.material.conductor.resistivity_ohm_m()))
}
fn linked_initial<'a>(curve: &'a EquilibriumEnthalpyPhaseCurve, h: [f64; 4])
    -> DepositionSimulation<'a> {
    let domain = VoxelDomain::new(4,1,1,1.0).unwrap();
    let base = initial(&domain, curve, 0.0);
    let mut events: Vec<_> = base.cells().values().cloned().map(|mut cell| {
        cell.specific_h_j_kg = h[cell.id.0 as usize];
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    for c in 0..3 { events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: c, a: CellId(c), b: CellId(c+1), area_m2: 1.0, conductance_w_k: 0.0,
    }) }); }
    DepositionSimulation::new(events, &[(MAT, curve as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 4, max_links: 64, max_events: 7, max_substeps: 1000 })
        .unwrap().advance_to(0.0, config(), || Ok(())).unwrap().0
}
fn step_config() -> StepConfig { StepConfig { max_dt_s: 0.0125, energy_tolerance_j: 1e-8 } }

#[test]
fn one_phase_material_epoch_closes_current_mushy_flow_and_total_enthalpy() {
    let domain = VoxelDomain::new(4,1,1,1.0).unwrap();
    let gate = CancelGate::new();
    let ids = [CellId(0),CellId(1),CellId(2),CellId(3)];
    let flow_ids = ids.map(Some);
    let config_flow = momentum();
    let setup = setup(10.0);
    let contacts = [ElectricalContact { cells: [1,2], resistance_ohm_m2: 0.5,
        heat_fraction_to_first: 0.25 }];
    let mut outputs = Vec::new();
    for scale in [1.0, 2.0] {
        let card = source_card(scale, scale, 0.01);
        let curve = phase_curve(&card, [1.0;2]);
        let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
        let old = linked_initial(&curve, [50.0;4]);
        let queries = Rc::new(Cell::new(0));
        let reads = RefCell::new(Vec::new());
        let observer = MaterialFlowFeedback::try_new(&old, &domain, input(&flow_ids,&curve),
            &config_flow, Heater::default(), &gate, 100, |_, cell: &DepositedCell, t| {
                queries.set(queries.get()+1);
                resolve(&card,&curve,cell,t,queries.get())
            }).unwrap().with_mushy_mobility(&mobility).unwrap()
            .with_fluid_conduction(&old,&[0,1,2]).unwrap();
        let (next, accepted, report) = advance_joule_enthalpy_with_material_flow(
            &old,0.1,step_config(),electric_input(&domain,&setup,&ids,&contacts),&observer,
            |sample| { reads.borrow_mut().push((sample.voxel,sample.material.query)); resistivity(sample) },
        ).unwrap();
        let mut h = [50.0_f64;4];
        let mut boundary_j = 0.0; let mut electrode_j = 0.0;
        let mut last_flow = 0.0; let mut last_current = 0.0;
        for &dt in &accepted.inner().durations {
            let phase = h.map(|h| h/100.0);
            let mu = phase.map(|f| scale*(1.0+f));
            let rho_e = phase.map(|f| scale*(1.0+f));
            let drag = std::array::from_fn::<_,4,_>(|c|
                mu[c]*(1.0-phase[c]).powi(2)/phase[c].powi(3));
            last_flow = 4.0/drag.iter().sum::<f64>();
            last_current = 10.0/(rho_e[1]+rho_e[2]+0.5);
            let power = [0.0,last_current.powi(2)*(rho_e[1]+0.125),
                last_current.powi(2)*(rho_e[2]+0.375),0.0];
            boundary_j += dt*last_flow*(50.0-h[3]);
            electrode_j += dt*10.0*last_current;
            let previous = h;
            for c in 0..4 {
                let upstream = if c==0 {50.0} else {previous[c-1]};
                h[c] += dt*((if c==1 {100.0} else {0.0})+power[c]
                    +last_flow*(upstream-previous[c]));
            }
        }
        let epoch = accepted.applied().unwrap();
        near(epoch.flow.field.flux_x(2,0,0),last_flow);
        near(report.last_applied.as_ref().unwrap().electric.electrodes[0].current_a,last_current);
        for sample in accepted.endpoint() {
            near(next.cells()[&sample.cell].specific_h_j_kg,h[sample.voxel]);
            assert_eq!(sample.temperature_k,500.0);
            assert_eq!(sample.material.conductor.material().card_identity(),card.content_hash());
            near(sample.transport.dynamic_viscosity_pa_s(),scale*(1.0+h[sample.voxel]/100.0));
            near(sample.transport.thermal_conductivity_w_m_k(),0.01*(1.0+h[sample.voxel]/100.0));
        }
        near(report.electrode_energy_j,electrode_j);
        near(report.deposited_joule_j,electrode_j);
        near(accepted.audit().external_enthalpy_j,boundary_j);
        near(next.energy_j()-old.energy_j(),10.0+electrode_j+boundary_j);
        near(report.thermal.external_heat_j,10.0+electrode_j+boundary_j);
        assert_eq!(next.mass_kg(),old.mass_kg());
        let steps = accepted.momentum_solves();
        assert!(steps>2);
        assert_eq!(report.current_solves,steps+1);
        assert_eq!(queries.get(),4*(1+2*steps)); // joining adds ZERO material queries
        let reads = reads.borrow();
        assert_eq!(reads.len(),2*(steps+1)); // ideal electrodes use no bulk getter
        assert_eq!(&reads[reads.len()-2..],&[
            (1,accepted.endpoint()[1].material.query),(2,accepted.endpoint()[2].material.query)]);
        assert_eq!(&reads[reads.len()-4..reads.len()-2],&[
            (1,epoch.samples[1].material.query),(2,epoch.samples[2].material.query)]);
        assert_eq!(observer.momentum_solves(),0);
        assert!(observer.inner().durations.is_empty());
        outputs.push((last_current,last_flow,h[1]));
    }
    assert!(outputs[1].0<outputs[0].0);
    assert!(outputs[1].1<outputs[0].1);
    assert!((outputs[1].2-outputs[0].2).abs()>0.1);
}

#[test]
fn the_same_samples_conduct_heat_and_zero_voltage_recovers_the_original_flow_owner() {
    let domain = VoxelDomain::new(4,1,1,1.0).unwrap();
    let gate = CancelGate::new();
    let ids = [CellId(0),CellId(1),CellId(2),CellId(3)];
    let flow_ids = ids.map(Some);
    let config_flow = throughflow();
    let mut outputs = Vec::new();
    for conductivity in [0.01,5.0] {
        let card = source_card(1.0,1.0,conductivity);
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),[
            EnthalpyPhaseKnot { specific_enthalpy_j_kg:0.0,temperature_k:300.0,
                liquid_mass_fraction:1.0,bulk_density_kg_m3:1.0 },
            EnthalpyPhaseKnot { specific_enthalpy_j_kg:200.0,temperature_k:400.0,
                liquid_mass_fraction:1.0,bulk_density_kg_m3:1.0 },
        ].to_vec()).unwrap();
        let old = linked_initial(&curve,[0.0,20.0,40.0,60.0]);
        let mut flow_input = input(&flow_ids,&curve);
        flow_input.inlet_specific_h_j_kg = [Some(0.0);6];
        let observer = MaterialFlowFeedback::try_new(&old,&domain,flow_input,&config_flow,
            Heater::default(),&gate,100,|_,cell:&DepositedCell,t|resolve(&card,&curve,cell,t,0))
            .unwrap().with_fluid_conduction(&old,&[0,1,2]).unwrap();
        let voltage = setup(10.0);
        let (next,accepted,report) = advance_joule_enthalpy_with_material_flow(
            &old,0.1,step_config(),electric_input(&domain,&voltage,&ids,&[]),&observer,resistivity).unwrap();
        let mut h = [0.0,20.0,40.0,60.0];
        for &dt in &accepted.inner().durations {
            let previous = h;
            for c in 0..4 {
                let inlet = if c==0 {0.0} else {previous[c-1]};
                let mut watts = 0.25*(inlet-previous[c]);
                if c>0 {watts += conductivity*(previous[c-1]-previous[c]);}
                if c<3 {watts += conductivity*(previous[c+1]-previous[c]);}
                if c==1 {watts +=100.0;}
                if c==1 || c==2 {watts +=12.5;} // I=10/(2+2), I²*rho_e=12.5 W
                h[c] += dt*watts;
            }
        }
        for (c,expected) in h.into_iter().enumerate() {
            near(next.cells()[&ids[c]].specific_h_j_kg,expected);
        }
        near(report.deposited_joule_j,2.5);
        let zero = setup(0.0);
        let (off,_,electrical) = advance_joule_enthalpy_with_material_flow(
            &old,0.1,step_config(),electric_input(&domain,&zero,&ids,&[]),&observer,resistivity).unwrap();
        let (baseline,_,_) = old.advance_observed(&observer,0.1,step_config(),||Ok(())).unwrap();
        assert_eq!(off.cells(),baseline.cells());
        assert_eq!(electrical.deposited_joule_j,0.0);
        outputs.push(next.cells()[&ids[0]].specific_h_j_kg);
    }
    assert!(outputs[1]>outputs[0]);
}

#[test]
fn mismatched_maps_endpoint_budget_and_source_refusal_leave_both_histories_unchanged() {
    let card = source_card(1.0,1.0,0.01);
    let curve = phase_curve(&card,[1.0;2]);
    let domain = VoxelDomain::new(4,1,1,1.0).unwrap();
    let copy = VoxelDomain::new(4,1,1,1.0).unwrap();
    let gate = CancelGate::new();
    let ids = [CellId(0),CellId(1),CellId(2),CellId(3)];
    let flow_ids = ids.map(Some);
    let config_flow = throughflow();
    let old = linked_initial(&curve,[50.0;4]);
    let observer = MaterialFlowFeedback::try_new(&old,&domain,input(&flow_ids,&curve),
        &config_flow,Heater::default(),&gate,100,|_,cell:&DepositedCell,t|
            resolve(&card,&curve,cell,t,0)).unwrap();
    let setup = setup(10.0);
    assert!(advance_joule_enthalpy_with_material_flow(&old,0.01,step_config(),
        electric_input(&copy,&setup,&ids,&[]),&observer,resistivity).is_err());
    let wrong = [ids[1],ids[0],ids[2],ids[3]];
    assert!(advance_joule_enthalpy_with_material_flow(&old,0.01,step_config(),
        electric_input(&domain,&setup,&wrong,&[]),&observer,resistivity).is_err());
    let mut limited = electric_input(&domain,&setup,&ids,&[]);
    limited.max_current_solves = 1; // the applied solve succeeds; final support refuses
    assert!(matches!(advance_joule_enthalpy_with_material_flow(&old,0.01,step_config(),
        limited,&observer,resistivity),Err(MaterialFlowElectricalError::Electrical(JouleEnthalpyError::CurrentBudget))));
    assert!(advance_joule_enthalpy_with_material_flow(&old,0.01,step_config(),
        electric_input(&domain,&setup,&ids,&[]),&observer,|sample| {
            if sample.specific_enthalpy_j_kg>50.0 {
                Err(ChtError::InvalidInput {field:"fixture.source",reason:"endpoint source refused".into()})
            } else {resistivity(sample)}
        }).is_err());
    assert_eq!(old.time_s(),0.0);
    assert!(old.cells().values().all(|cell|cell.specific_h_j_kg==50.0));
    assert_eq!(observer.momentum_solves(),0);
    assert!(observer.applied().is_none());
    assert!(observer.inner().durations.is_empty());
    let (a,fa,ra) = advance_joule_enthalpy_with_material_flow(&old,0.01,step_config(),
        electric_input(&domain,&setup,&ids,&[]),&observer,resistivity).unwrap();
    let (b,fb,rb) = advance_joule_enthalpy_with_material_flow(&old,0.01,step_config(),
        electric_input(&domain,&setup,&ids,&[]),&observer,resistivity).unwrap();
    assert_eq!(a.cells(),b.cells());
    assert_eq!(fa.endpoint(),fb.endpoint());
    assert_eq!(ra,rb);
    // Reusing the old observer for a changed zero-duration state must not use stale samples.
    assert!(advance_joule_enthalpy_with_material_flow(&a,a.time_s(),step_config(),
        electric_input(&domain,&setup,&ids,&[]),&observer,resistivity).is_err());
}

