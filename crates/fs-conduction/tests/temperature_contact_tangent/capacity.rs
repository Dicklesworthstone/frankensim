//! Fresh coupled trajectory probes: never differentiate a frozen R(T) endpoint.
use super::*;
use fs_conduction::{ConductionProblem, ConductivityModel, ConductivityTable, ScalarField, SolveConfig};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::interface::temperature::ContactIterationConfig;

const INITIAL: [f64;8]=[300.,300.,300.,300.,360.,360.,360.,360.];
const SOURCE: [f64;8]=[0.,0.,0.,0.,30.,30.,30.,30.];
const D_INITIAL: [f64;8]=[2.2,1.9,2.7,2.5,1.7,2.9,2.4,1.4];
const D_SOURCE: [f64;8]=[2.,3.,4.,1.,-1.,2.,-3.,4.];
const DT: [f64;3]=[0.15,0.35,0.5];
fn config() -> StepConfig {
    let mut linear=SolveConfig::default().linear;
    linear.tolerance=1e-12;linear.max_iterations=2048;linear.restart=8;
    StepConfig {linear,energy_tolerance_j:1e-9}
}
fn coupling() -> ContactIterationConfig {
    ContactIterationConfig {max_iterations:96,relaxation:0.5,residual_tolerance_j:1e-10}
}
fn material(variable: bool) -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(if variable {
        vec![(280.,1.6),(380.,3.6)]
    } else {vec![(280.,2.),(380.,2.)]}).unwrap())
}
fn nonlinear() -> NonlinearStepConfig {
    NonlinearStepConfig {residual_rtol:1e-13,residual_atol_j:1e-12,..Default::default()}
}

fn trajectory(cx: &Cx<'_>, variable: bool, delta: f64, differentiate: bool) -> (Vec<f64>,f64,f64) {
    let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let card=card(1.+0.1*delta);let set=contacts(cx,&mesh,&bc,&card,false);
    let capacities=[VolumetricHeatCapacity::declared(10.*(1.+0.04*delta)).unwrap(),
        VolumetricHeatCapacity::declared(20.*(1.+0.04*delta)).unwrap()];
    let engine=BackwardEuler::per_element(cx,&mesh,&capacities).unwrap();
    let material=material(variable);
    let source=ScalarField::Nodal(SOURCE.iter().zip(D_SOURCE).map(|(q,d)|q+delta*d).collect());
    let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
        element_materials:None,source:&source};
    let mut states=vec![INITIAL.iter().zip(D_INITIAL).map(|(t,d)|t+delta*d).collect::<Vec<_>>()];
    for dt in DT {
        let old=states.last().unwrap();let before=old.clone();
        let next=set.advance_capacity(cx,&engine,problem,old,dt,config(),Some(nonlinear()),coupling()).unwrap();
        assert_eq!(&before,old);
        states.push(next.conduction.temperature);
    }
    if !differentiate {return (states.pop().unwrap(),0.,0.);}
    assert_eq!(delta,0.);
    let mut carry=vec![0.;8];carry[0]=1.;
    let mut frozen_carry=carry.clone();let mut derivative=0.;let mut frozen_derivative=0.;
    for index in (0..DT.len()).rev() {
        let linear=engine.linearize_temperature_contact_endpoint(
            cx,problem,&set,&states[index],&states[index+1],DT[index],config(),
            coupling().residual_tolerance_j,&[],1,
        ).unwrap();
        assert_eq!(linear.primal().temperature,states[index+1]);
        assert_eq!(linear.primal().krylov_iterations,0,"no replacement primal solve");
        assert!(linear.uses_nonlinear_jacobian(),"R-prime needs nonsymmetric solves even for constant k");
        let gradient=linear.pullback_with_resistance_controls(cx,&carry,&[],&[],1).unwrap();
        assert_eq!(gradient.resistance_multipliers.len(),1);
        assert_eq!(gradient.resistance_multipliers[0].interface,"joint");
        assert_eq!(gradient.resistance_multipliers[0].card_identity,card.content_hash());
        let dual=&gradient.transport.nodal_load;
        assert!(gradient.transport.relative_residual<config().linear.tolerance);
        derivative+=dot(&linear.source_density_pullback(cx,dual).unwrap(),&D_SOURCE)
            +0.04*linear.capacity_multiplier_pullback(cx,dual).unwrap()
            +0.1*gradient.resistance_multipliers[0].derivative;
        carry=linear.previous_temperature_pullback(cx,dual).unwrap();
        // Explicitly wrong comparison: same accepted T and actual R values,
        // but omitting the temperature dependence of that source law.
        let physical=linear.contact().physical_contacts().frozen_interfaces();
        let frozen=engine.linearize_prescribed_endpoint(cx,problem,Some(physical),
            &states[index],&states[index+1],DT[index],config(),coupling().residual_tolerance_j,&[]).unwrap();
        let fg=frozen.pullback(cx,&frozen_carry,&[],&[]).unwrap();
        frozen_derivative+=dot(&frozen.source_density_pullback(cx,&fg.nodal_load).unwrap(),&D_SOURCE)
            +0.04*frozen.capacity_multiplier_pullback(cx,&fg.nodal_load).unwrap()
            +0.1*physical.matching_resistance_scale_pullback(cx,&states[index+1],&fg.nodal_load,1).unwrap()[0].derivative;
        frozen_carry=frozen.previous_temperature_pullback(cx,&fg.nodal_load).unwrap();
    }
    derivative+=dot(&carry,&D_INITIAL);frozen_derivative+=dot(&frozen_carry,&D_INITIAL);
    (states.pop().unwrap(),derivative,frozen_derivative)
}

#[test]
fn g3_capacity_contact_adjoint_matches_fresh_histories_sources_capacities_and_source_curves() {
    context(&CancelGate::new_clock_free(),|cx| {
        for (variable,expected_t,expected_d) in [
            (false,323.48193256858633,0.6140567654970317),
            (true,323.3907563852738,0.5977465447399422),
        ] {
            let (field,derivative,frozen)=trajectory(cx,variable,0.,true);
            near(field[0],expected_t,2e-6);near(derivative,expected_d,2e-6);
            assert!((derivative-frozen).abs()>0.2,"fixture must expose the missing contact-state feedback");
            for epsilon in [0.001,0.0005] {
                let plus=trajectory(cx,variable,epsilon,false).0;
                let minus=trajectory(cx,variable,-epsilon,false).0;
                near(derivative,(plus[0]-minus[0])/(2.*epsilon),5e-5);
            }
            let (again,repeat,_)=trajectory(cx,variable,0.,true);
            assert_eq!(field,again);assert_eq!(derivative,repeat);
        }
    });
}

#[test]
fn g0_capacity_contact_linearization_preserves_physical_and_work_gates() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card=card(1.);let set=contacts(cx,&mesh,&bc,&card,false);
        let material=material(false);let source=ScalarField::Uniform(0.);
        let problem=ConductionProblem {mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let engine=BackwardEuler::uniform(cx,&mesh,VolumetricHeatCapacity::declared(10.).unwrap()).unwrap();
        let old=[330.;8];
        let call=|candidate:&[f64],limit,config| engine.linearize_temperature_contact_endpoint(
            cx,problem,&set,&old,candidate,0.35,config,1e-10,&[],limit,
        );
        assert!(call(&[331.;8],1,config()).is_err(),"unchanged target rejects an unbalanced candidate");
        assert!(call(&old,0,config()).is_err());
        let mut bad=config();bad.linear.restart=0;
        assert!(call(&old,1,bad).is_err(),"no silent symmetric fallback");
        let linear=call(&old,1,config()).unwrap();
        assert!(linear.pullback_with_resistance_controls(cx,&[0.;8],&[],&[],0).is_err());
        let zero=linear.pullback_with_resistance_controls(cx,&[0.;8],&[],&[],1).unwrap();
        assert_eq!(zero.transport.iterations,0);assert_eq!(zero.transport.nodal_load,vec![0.;8]);
        assert_eq!(zero.resistance_multipliers[0].derivative,0.);
        let gate=CancelGate::new_clock_free();gate.request();
        context(&gate,|stopped| {
            assert!(matches!(linear.pullback_with_resistance_controls(stopped,&[0.;8],&[],&[],1),
                Err(ConductionError::Cancelled {..})));
        });
    });
}
