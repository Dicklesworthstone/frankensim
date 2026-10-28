//! G0/G1/G3: independent slab formulas, frames, source domains and conservation.
use super::*;
use crate::gas::GasSpec;
use crate::thermochemical_gas::tests::{model, point};
use crate::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn support() -> PorousGasSupport {
    PorousGasSupport { temperature_k: (250.0, 600.0), pressure_pa: (10000.0, 1e6),
        maximum_superficial_speed_m_s: 20.0 }
}
fn gas() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 150000.0).unwrap() }
fn layer(k: f64, thickness: f64) -> PorousGasLayer {
    PorousGasLayer::try_new(IntrinsicPermeability::isotropic(k).unwrap(), thickness, 0.0,
        "synthetic Darcy slab, not a measured material").unwrap()
}
fn barrier(k: f64, thickness: f64) -> DarcyGasBarrier {
    DarcyGasBarrier::try_new(vec![layer(k, thickness)], [1.0,0.0,0.0], support()).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 3e-12*a.abs().max(b.abs()).max(1e-300), "{a} != {b}"); }

#[test]
fn compressible_slab_conserves_mass_not_equal_volume_flow() {
    let g = gas(); let b = barrier(1e-12, 0.01);
    let q = b.mass_flow(&g, 200000.0, 100000.0, 0.02).unwrap();
    let exact = 0.02*1e-12*(200000.0_f64.powi(2)-100000.0_f64.powi(2))
        /(2.0*g.dynamic_viscosity*g.specific_gas_constant*g.temperature*0.01);
    near(q.mass_rate_left_to_right_kg_s, exact);
    near(q.right_volume_rate_m3_s, 2.0*q.left_volume_rate_m3_s);
    let rt = g.specific_gas_constant*g.temperature;
    near(q.left_volume_rate_m3_s*200000.0/rt, exact);
    near(q.right_volume_rate_m3_s*100000.0/rt, exact);
    near(q.drag_dissipation_w, exact*rt*2.0_f64.ln());
    let reverse = b.mass_flow(&g, 100000.0, 200000.0, 0.02).unwrap();
    assert_eq!(reverse.mass_rate_left_to_right_kg_s, -q.mass_rate_left_to_right_kg_s);
    assert_eq!(reverse.drag_dissipation_w, q.drag_dissipation_w);
    assert_eq!(reverse.left_volume_rate_m3_s, -q.right_volume_rate_m3_s);
}

#[test]
fn layered_harmonic_resistance_and_exact_sealed_limits() {
    let g = gas();
    let layers = vec![layer(1e-12, 0.01), layer(4e-12, 0.02)];
    let b = DarcyGasBarrier::try_new(layers, [0.0,1.0,0.0], support()).unwrap();
    assert_eq!(b.geometric_resistance_per_m(), Some(1.5e10));
    let equivalent = barrier(2e-12, 0.03);
    near(b.mass_flow(&g, 200000.0, 100000.0, 0.02).unwrap().mass_rate_left_to_right_kg_s,
        equivalent.mass_flow(&g, 200000.0, 100000.0, 0.02).unwrap().mass_rate_left_to_right_kg_s);
    let sealed = DarcyGasBarrier::try_new(vec![layer(1e-12, 0.01), layer(0.0, 0.001)],
        [1.0,0.0,0.0], support()).unwrap();
    assert_eq!(sealed.geometric_resistance_per_m(), None);
    assert_eq!(sealed.mass_flow(&g, 200000.0, 100000.0, 0.02).unwrap().mass_rate_left_to_right_kg_s, 0.0);
    assert_eq!(sealed.linear_volume_resistance(&g, 0.02).unwrap(), None);
    assert_eq!(b.mass_flow(&g, 200000.0, 100000.0, 0.0).unwrap().mass_rate_left_to_right_kg_s, 0.0);
    assert_eq!(b.mass_flow(&g, 150000.0, 150000.0, 0.02).unwrap().drag_dissipation_w, 0.0);
}

#[test]
fn orthotropic_normal_and_rotated_frame_obey_tensor_law() {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let k = IntrinsicPermeability::orthotropic([1e-12, 9e-12, 0.0],
        [[s,s,0.0],[-s,s,0.0],[0.0,0.0,1.0]]).unwrap();
    near(k.normal_permeability([1.0,0.0,0.0]).unwrap(), 5e-12);
    near(k.normal_permeability([s,s,0.0]).unwrap(), 1e-12);
    assert_eq!(k.normal_permeability([0.0,0.0,1.0]).unwrap(), 0.0);
    assert!(k.normal_permeability([2.0,0.0,0.0]).is_err());
    assert!(IntrinsicPermeability::orthotropic([1e-12;3], [[1.,0.,0.];3]).is_err());
    let world = IntrinsicPermeability::orthotropic([1e-12,9e-12,0.0],
        [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).unwrap();
    near(k.normal_permeability([1.,0.,0.]).unwrap(), world.normal_permeability([s,-s,0.]).unwrap());
}

#[test]
fn area_thickness_permeability_and_gas_substitution_scale_physically() {
    let g = gas(); let b = barrier(1e-12, 0.01);
    let q = b.mass_flow(&g, 160000.0, 140000.0, 0.01).unwrap().mass_rate_left_to_right_kg_s;
    near(b.mass_flow(&g, 160000.0, 140000.0, 0.02).unwrap().mass_rate_left_to_right_kg_s, 2.0*q);
    near(barrier(2e-12,0.01).mass_flow(&g,160000.0,140000.0,0.01).unwrap().mass_rate_left_to_right_kg_s, 2.0*q);
    near(barrier(1e-12,0.02).mass_flow(&g,160000.0,140000.0,0.01).unwrap().mass_rate_left_to_right_kg_s, 0.5*q);
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 500.0, 150000.0).unwrap();
    let q_hot = b.mass_flow(&hot,160000.0,140000.0,0.01).unwrap().mass_rate_left_to_right_kg_s;
    near(q_hot/q, g.dynamic_viscosity*g.temperature/(hot.dynamic_viscosity*hot.temperature));
    near(b.linear_volume_resistance(&g,0.01).unwrap().unwrap(), g.dynamic_viscosity*0.01/(1e-12*0.01));
}

#[test]
fn tiny_pressure_drop_keeps_flux_and_positive_drag() {
    let g=gas(); let b=barrier(1e-12,0.01); let p=150000.0_f64;
    let next=f64::from_bits(p.to_bits()+1);
    let q=b.mass_flow(&g,next,p,0.01).unwrap();
    assert!(q.mass_rate_left_to_right_kg_s>0.0 && q.drag_dissipation_w>0.0);
    near(q.mass_rate_left_to_right_kg_s,
        (next-p)*0.01*1e-12/(g.dynamic_viscosity*0.01)*((0.5*next+0.5*p)/(g.specific_gas_constant*g.temperature)));
}

#[test]
fn klinkenberg_is_exact_and_bound_to_calibrated_gas_composition() {
    let m=model(vec![0.75,0.25],0.0); let state=m.evaluate(&point(300.0,150000.0)).unwrap();
    let mut l=layer(1e-12,0.01); l.klinkenberg_slip_pa=50000.0;
    assert!(DarcyGasBarrier::try_new(vec![l.clone()],[1.,0.,0.],support()).is_err());
    let b=DarcyGasBarrier::try_new_for_gas_model(vec![l.clone()],[1.,0.,0.],support(),&m).unwrap();
    let q=b.mass_flow_for_model(state.state(),&m,200000.0,100000.0,0.01).unwrap();
    let base=barrier(1e-12,0.01).mass_flow(state.state(),200000.0,100000.0,0.01).unwrap();
    near(q.mass_rate_left_to_right_kg_s/base.mass_rate_left_to_right_kg_s,4.0/3.0);
    assert!(b.mass_flow(state.state(),200000.0,100000.0,0.01).is_err());
    assert!(b.mass_flow_for_model(state.state(),&model(vec![0.25,0.75],0.0),200000.0,100000.0,0.01).is_err());
    assert!(b.linear_volume_resistance(state.state(),0.01).is_err());
    assert_eq!(b.klinkenberg_slip_pa(),Some(50000.0));
    let heterogeneous=DarcyGasBarrier::try_new_for_gas_model(vec![l,layer(1e-12,0.01)],
        [1.,0.,0.],support(),&m).unwrap();
    assert_eq!(heterogeneous.klinkenberg_slip_pa(),None);
}

fn slip_layer(k: f64, thickness: f64, slip: f64) -> PorousGasLayer {
    PorousGasLayer::try_new(IntrinsicPermeability::isotropic(k).unwrap(),thickness,slip,
        "synthetic layer-specific gas calibration, not experimental evidence").unwrap()
}

fn manufactured_laminate(m: &ThermochemicalGasModel) -> DarcyGasBarrier {
    // Each r_i equals its independent integrated potential for C=j*mu*R*T=1:
    // 90000*(195000+10000) and 70000*(115000+100000), respectively.
    DarcyGasBarrier::try_new_for_gas_model(vec![
        slip_layer(1e-12,0.01845,10000.0),slip_layer(1e-12,0.01505,100000.0),
    ],[1.,0.,0.],support(),m).unwrap()
}

#[test]
fn heterogeneous_slip_recovers_manufactured_interface_and_each_layer_mass_balance() {
    let m=model(vec![0.75,0.25],0.0); let state=m.evaluate(&point(300.0,150000.0)).unwrap();
    let g=state.state(); let b=manufactured_laminate(&m);
    let profile=b.mass_flow_profile_for_model(g,&m,240000.0,80000.0,0.02).unwrap();
    let pressures=profile.interface_pressures_pa.as_ref().unwrap();
    for (actual,expected) in pressures.iter().zip([240000.0,150000.0,80000.0]) { near(*actual,expected); }
    let scaled_flux=profile.flow.mass_flux_left_to_right_kg_m2_s*g.dynamic_viscosity
        *g.specific_gas_constant*g.temperature;
    near(scaled_flux,1.0);
    for (i,layer) in b.layers().iter().enumerate() {
        let r=layer.thickness_m()/layer.permeability().normal_permeability(b.normal()).unwrap();
        let independent=(pressures[i]-pressures[i+1])
            *(0.5*pressures[i]+0.5*pressures[i+1]+layer.klinkenberg_slip_pa())/r;
        near(independent,scaled_flux);
    }
    assert!(profile.iterations>0 && profile.iterations<=MAX_LAYERED_DARCY_ITERATIONS);
    assert!(profile.normalized_pressure_residual.abs()<=b.pressure_residual_tolerance());
    assert_eq!(profile.flow,b.mass_flow_for_model(g,&m,240000.0,80000.0,0.02).unwrap());
    let averaged=(240000.0-80000.0)*(160000.0+55000.0)/b.geometric_resistance_per_m().unwrap();
    assert!((averaged-scaled_flux).abs()>0.01*scaled_flux);
}

#[test]
fn heterogeneous_orientation_reversal_and_layer_subdivision_preserve_the_physical_stack() {
    let m=model(vec![0.75,0.25],0.0); let state=m.evaluate(&point(300.0,150000.0)).unwrap();
    let g=state.state(); let b=manufactured_laminate(&m);
    let original=b.mass_flow_profile_for_model(g,&m,240000.0,80000.0,0.02).unwrap();
    let mut reversed=b.layers().to_vec(); reversed.reverse();
    let reversed=DarcyGasBarrier::try_new_for_gas_model(reversed,[1.,0.,0.],support(),&m).unwrap();
    let flipped=reversed.mass_flow_profile_for_model(g,&m,80000.0,240000.0,0.02).unwrap();
    assert_eq!(flipped.flow.mass_rate_left_to_right_kg_s,-original.flow.mass_rate_left_to_right_kg_s);
    let mut pressures=original.interface_pressures_pa.clone().unwrap(); pressures.reverse();
    assert_eq!(flipped.interface_pressures_pa.unwrap(),pressures);
    // Swapping only reservoirs is physically different for unequal b_i.
    let opposed=b.mass_flow_for_model(g,&m,80000.0,240000.0,0.02).unwrap();
    assert!((opposed.mass_rate_left_to_right_kg_s+original.flow.mass_rate_left_to_right_kg_s).abs()
        >1e-3*original.flow.mass_rate_left_to_right_kg_s);
    let split=b.layers().iter().flat_map(|layer| {
        let mut half=layer.clone(); half.thickness_m*=0.5; [half.clone(),half]
    }).collect();
    let split=DarcyGasBarrier::try_new_for_gas_model(split,[1.,0.,0.],support(),&m).unwrap();
    let refined=split.mass_flow_profile_for_model(g,&m,240000.0,80000.0,0.02).unwrap();
    near(refined.flow.mass_rate_left_to_right_kg_s,original.flow.mass_rate_left_to_right_kg_s);
    near(refined.interface_pressures_pa.unwrap()[2],150000.0);
    let scaled=b.layers().iter().map(|layer| {
        let mut layer=layer.clone(); layer.thickness_m*=2.0; layer
    }).collect();
    let scaled=DarcyGasBarrier::try_new_for_gas_model(scaled,[1.,0.,0.],support(),&m).unwrap();
    let changed=scaled.mass_flow_profile_for_model(g,&m,240000.0,80000.0,0.04).unwrap();
    near(changed.flow.mass_rate_left_to_right_kg_s,original.flow.mass_rate_left_to_right_kg_s);
    assert_eq!(changed.interface_pressures_pa,original.interface_pressures_pa);
}

#[test]
fn heterogeneous_equilibrium_and_adjacent_float_drop_use_exact_small_signal_series_limit() {
    let m=model(vec![0.75,0.25],0.0); let state=m.evaluate(&point(300.0,150000.0)).unwrap();
    let g=state.state(); let b=manufactured_laminate(&m); let p=150000.0_f64;
    let resistance:f64=b.layers().iter().map(|layer| layer.thickness_m()
        /layer.permeability().normal_permeability(b.normal()).unwrap()
        /(1.0+layer.klinkenberg_slip_pa()/p)).sum();
    let conductance=0.02/(g.dynamic_viscosity*resistance);
    let equilibrium=b.mass_flow_profile_for_model(g,&m,p,p,0.02).unwrap();
    near(equilibrium.flow.mean_volume_conductance_m3_pa_s,conductance);
    assert_eq!(equilibrium.flow.mass_rate_left_to_right_kg_s,0.0);
    assert_eq!(equilibrium.interface_pressures_pa.unwrap(),vec![p;3]);
    assert_eq!(equilibrium.iterations,0);
    let next=f64::from_bits(p.to_bits()+1);
    let tiny=b.mass_flow_profile_for_model(g,&m,next,p,0.02).unwrap();
    near(tiny.flow.mass_rate_left_to_right_kg_s,
        conductance*(next-p)*p/(g.specific_gas_constant*g.temperature));
    assert!(tiny.flow.drag_dissipation_w>0.0);
    assert!(tiny.normalized_pressure_residual.abs()<=b.pressure_residual_tolerance());
}

#[test]
fn common_factor_profile_keeps_analytic_fast_path_and_sealed_interiors_are_not_invented() {
    let m=model(vec![0.75,0.25],0.0); let state=m.evaluate(&point(300.0,150000.0)).unwrap();
    let g=state.state();
    for slip in [0.0,50000.0] {
        let b=DarcyGasBarrier::try_new_for_gas_model(vec![slip_layer(1e-12,0.01,slip);2],
            [1.,0.,0.],support(),&m).unwrap();
        let profile=b.mass_flow_profile_for_model(g,&m,200000.0,100000.0,0.01).unwrap();
        assert_eq!(profile.iterations,0);
        assert_eq!(b.klinkenberg_slip_pa(),Some(slip));
        let expected=(0.5*(200000.0+slip).powi(2)+0.5*(100000.0+slip).powi(2)).sqrt()-slip;
        near(profile.interface_pressures_pa.unwrap()[1],expected);
        near(profile.flow.mass_flux_left_to_right_kg_m2_s,
            100000.0*(150000.0+slip)/(g.dynamic_viscosity*g.specific_gas_constant*g.temperature*2e10));
    }
    let sealed=DarcyGasBarrier::try_new_for_gas_model(vec![layer(0.0,0.01),slip_layer(1e-12,0.01,50000.0)],
        [1.,0.,0.],support(),&m).unwrap();
    for (b,area) in [(sealed,0.01),(manufactured_laminate(&m),0.0)] {
        let profile=b.mass_flow_profile_for_model(g,&m,200000.0,100000.0,area).unwrap();
        assert_eq!(profile.flow.mass_rate_left_to_right_kg_s,0.0);
        assert_eq!(profile.interface_pressures_pa,None);
        assert_eq!(profile.iterations,0);
    }
}

#[test]
fn heterogeneous_slip_preserves_calibration_source_support_and_bounded_refusals() {
    let m=model(vec![0.75,0.25],0.0); let state=m.evaluate(&point(300.0,150000.0)).unwrap();
    let g=state.state();
    let source=PorousGasLayer::from_material_isotropic(material(INTRINSIC_PERMEABILITY_DIMS),
        "intrinsic_permeability",0.01,None).unwrap();
    let layers=vec![source,slip_layer(1e-12,0.01,50000.0)];
    assert!(DarcyGasBarrier::try_new(layers.clone(),[1.,0.,0.],support()).is_err());
    let b=DarcyGasBarrier::try_new_for_gas_model(layers,[1.,0.,0.],support(),&m).unwrap();
    assert!(b.mass_flow(g,200000.0,100000.0,0.01).is_err());
    assert!(b.mass_flow_profile(g,200000.0,100000.0,0.01).is_err());
    assert!(b.linear_volume_resistance(g,0.01).is_err());
    assert!(b.mass_flow_profile_for_model(g,&model(vec![0.25,0.75],0.0),200000.0,100000.0,0.01).is_err());
    assert!(b.mass_flow_profile_for_model(g,&m,200000.0,99999.0,0.01).is_err());
    let profile=b.mass_flow_profile_for_model(g,&m,200000.0,100000.0,0.01).unwrap();
    assert!(profile.interface_pressures_pa.unwrap().iter().all(|p| (100000.0..=200000.0).contains(p)));
    let query=LayeredQuery::new(200000.0,100000.0,b.maximum_slip_pa).unwrap();
    assert!(b.solve_layered(query,0).is_err());
    let nontrivial=manufactured_laminate(&m);
    let query=LayeredQuery::new(240000.0,80000.0,nontrivial.maximum_slip_pa).unwrap();
    assert!(nontrivial.solve_layered(query,1).is_err());
}

#[test]
fn extreme_slip_contrast_and_maximum_laminate_size_remain_bounded_and_repeatable() {
    let m=model(vec![0.75,0.25],0.0); let state=m.evaluate(&point(300.0,150000.0)).unwrap();
    let g=state.state();
    let b=DarcyGasBarrier::try_new_for_gas_model(vec![layer(1e-12,0.01),slip_layer(1e-12,0.01,1e100)],
        [1.,0.,0.],support(),&m).unwrap();
    let result=b.mass_flow_profile_for_model(g,&m,200000.0,100000.0,0.01).unwrap();
    let limit=barrier(1e-12,0.01).mass_flow(g,200000.0,100000.0,0.01).unwrap();
    near(result.flow.mass_rate_left_to_right_kg_s,limit.mass_rate_left_to_right_kg_s);
    assert!(result.iterations<=MAX_LAYERED_DARCY_ITERATIONS);
    let layers=(0..MAX_POROUS_LAYERS).map(|i| slip_layer(1e-12,0.01,
        if i%2==0 {0.0} else {50000.0})).collect();
    let b=DarcyGasBarrier::try_new_for_gas_model(layers,[1.,0.,0.],support(),&m).unwrap();
    let first=b.mass_flow_profile_for_model(g,&m,200000.0,100000.0,0.01).unwrap();
    assert_eq!(first,b.mass_flow_profile_for_model(g,&m,200000.0,100000.0,0.01).unwrap());
    let pressures=first.interface_pressures_pa.unwrap();
    assert_eq!(pressures.len(),MAX_POROUS_LAYERS+1);
    assert!(pressures.windows(2).all(|p| p[0]>=p[1] && p[1]>0.0));
    assert!(first.normalized_pressure_residual.abs()<=b.pressure_residual_tolerance());
}

fn material(dims: Dims) -> ResolvedMaterialStatePoint {
    material_axes(dims,0,"temperature","pressure")
}
fn material_axes(dims: Dims, convention: u8, temperature_axis: &str,
    pressure_axis: &str) -> ResolvedMaterialStatePoint {
    let mut claims=ClaimSet::new();
    let mut validity=ValidityDomain::unconstrained();
    let mut query=fs_matdb::QueryPoint::new();
    for (axis,kind,dimension,lo,hi,value) in [
        (temperature_axis,QuantityKind::AbsoluteTemperature,fs_qty::Temperature::DIMS,280.0,350.0,300.0),
        (pressure_axis,QuantityKind::Pressure,fs_qty::Pressure::DIMS,100000.0,200000.0,150000.0)] {
        if convention==2 {
            validity=validity.with(axis,lo,hi); query=query.with(axis,value).unwrap();
        } else {
            let quantity=if convention==0 { QuantitySpec::semantic(SemanticType::new(kind,ValueForm::Static)) }
                else { QuantitySpec::dimensional(dimension) };
            validity=validity.with_quantity(axis,quantity,lo,hi);
            query=query.with_quantity(axis,quantity,value).unwrap();
        }
    }
    claims.insert_claim(PropertyClaim { key:PropertyKey::new("intrinsic_permeability",dims),
        value:PropertyValue::Scalar {value:1e-12,dims},
        validity,
        uncertainty:UncertaintyModel::Unstated, interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![], provenance:Provenance {source:"synthetic source-bound Darcy slab".into(),
            license:"CC0-1.0".into(),artifact:None} }).unwrap();
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic porous solid".into(),
        phase:"solid".into(),process:"analytical test".into(),revision:0},claims,vec![]).unwrap();
    let requirement=ScalarPropertyRequirement::try_new("intrinsic_permeability",dims,
        ScalarAdmissibility::NonNegative).unwrap();
    resolve_material_state_point(&card,&query,&[requirement],
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}

#[test]
fn source_bounds_are_intersected_and_wrong_quantities_refuse() {
    let source=material(INTRINSIC_PERMEABILITY_DIMS);
    let l=PorousGasLayer::from_material_isotropic(source.clone(),"intrinsic_permeability",0.01,None).unwrap();
    assert_eq!(l.source(),Some(&source));
    let b=DarcyGasBarrier::try_new(vec![l],[1.,0.,0.],support()).unwrap();
    assert_eq!(b.support().temperature_k,(280.0,350.0));
    assert_eq!(b.support().pressure_pa,(100000.0,200000.0));
    assert!(b.mass_flow(&gas(),200000.0,99000.0,0.01).is_err());
    let hot=GasState::try_new(&GasSpec::dry_air_ussa1976(),360.0,150000.0).unwrap();
    assert!(b.mass_flow(&hot,160000.0,140000.0,0.01).is_err());
    assert!(PorousGasLayer::from_material_isotropic(material(Dims([1,0,-1,0,0,0])),
        "intrinsic_permeability",0.01,None).is_err());
}

#[test]
fn legacy_and_dimension_only_domains_still_bind_actual_gas_and_aliases_refuse() {
    for convention in [1,2] {
        for temperature_axis in ["T","temperature"] {
            let source=material_axes(INTRINSIC_PERMEABILITY_DIMS,convention,temperature_axis,"pressure");
            let layer=PorousGasLayer::from_material_isotropic(source,"intrinsic_permeability",0.01,None).unwrap();
            let b=DarcyGasBarrier::try_new(vec![layer],[1.,0.,0.],support()).unwrap();
            assert_eq!(b.support().temperature_k,(280.0,350.0));
            assert_eq!(b.support().pressure_pa,(100000.0,200000.0));
            assert!(b.mass_flow(&gas(),210000.0,150000.0,0.01).is_err());
            let hot=GasState::try_new(&GasSpec::dry_air_ussa1976(),360.0,150000.0).unwrap();
            assert!(b.mass_flow(&hot,160000.0,140000.0,0.01).is_err());
        }
    }
    for (temperature,pressure) in [("reference_temperature","pressure"),("temperature","partial_pressure")] {
        let source=material_axes(INTRINSIC_PERMEABILITY_DIMS,0,temperature,pressure);
        assert!(PorousGasLayer::from_material_isotropic(source,"intrinsic_permeability",0.01,None).is_err());
    }
}

#[test]
fn invalid_inputs_and_speed_outside_linear_regime_refuse() {
    for k in [-1.0,f64::NAN,f64::INFINITY] { assert!(IntrinsicPermeability::isotropic(k).is_err()); }
    assert!(DarcyGasBarrier::try_new(vec![],[1.,0.,0.],support()).is_err());
    assert!(DarcyGasBarrier::try_new(vec![layer(1e-12,0.01);65],[1.,0.,0.],support()).is_err());
    assert!(DarcyGasBarrier::try_new(vec![layer(1e-12,0.01);64],[1.,0.,0.],support()).is_ok());
    let g=gas(); let b=barrier(1e-12,0.01);
    for pressure in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(b.mass_flow(&g,pressure,150000.0,0.01).is_err());
    }
    let mut incoherent=g; incoherent.density*=2.0;
    assert!(b.mass_flow(&incoherent,160000.0,140000.0,0.01).is_err());
    let limited=DarcyGasBarrier::try_new(vec![layer(1e-12,0.01)],[1.,0.,0.],
        PorousGasSupport {maximum_superficial_speed_m_s:1e-5,..support()}).unwrap();
    assert!(limited.mass_flow(&g,200000.0,100000.0,0.01).is_err());
}

#[test]
fn acoustic_and_darcy_owners_share_intrinsic_geometry() {
    let acoustic=super::super::gas::PorousGasLaw::try_new(
        super::super::PorousCard::melamine_uf_nguyen2024(),1.81e-5).unwrap();
    near(acoustic.intrinsic_permeability().normal_permeability([1.,0.,0.]).unwrap(),
        acoustic.viscous_permeability_m2());
    let l=PorousGasLayer::try_new(acoustic.intrinsic_permeability(),0.1,0.0,
        acoustic.reference_card().provenance).unwrap();
    let b=DarcyGasBarrier::try_new(vec![l],[1.,0.,0.],support()).unwrap();
    near(b.linear_volume_resistance(&gas(),1.0).unwrap().unwrap(),
        gas().dynamic_viscosity*0.1/acoustic.viscous_permeability_m2());
}
