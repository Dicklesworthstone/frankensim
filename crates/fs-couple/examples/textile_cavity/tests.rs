use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,PropertyClaim,PropertyKey,PropertyValue,Provenance,UncertaintyModel};
use fs_qty::Dims;

fn specimen(mass:f64,stiffness:f64)->MaterialCard {
    let mut claims=ClaimSet::new();
    for (key,value,dims) in [
        ("areal_density",mass,Dims([-2,1,0,0,0,0])),("thickness",0.001,Dims([1,0,0,0,0,0])),
        ("warp_stiffness",stiffness,Dims([0,1,-2,0,0,0])),("weft_stiffness",stiffness,Dims([0,1,-2,0,0,0])),
        ("poisson_warp_weft",0.0,Dims::NONE),("shear_stiffness",stiffness/2.0,Dims([0,1,-2,0,0,0])),
    ] {
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),value:PropertyValue::Scalar {value,dims},
            interpolation:InterpolationPolicy::ConstantWithinValidity,
            validity:ValidityDomain::unconstrained().with("T",290.0,310.0).with("ftx_characterization",0.0,0.0),
            uncertainty:UncertaintyModel::Unstated,observations:vec![],provenance:Provenance {
                source:"synthetic sheet for coupled-cavity regression; not measured fabric".into(),license:"test".into(),artifact:None}}).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic sheet".into(),phase:"solid".into(),
        process:"fixed synthetic specimen".into(),revision:0},claims,vec![]).unwrap()
}
fn input()->Input {Input::parse(include_str!("air.txt")).unwrap()}

#[test]
fn material_only_mass_substitution_propagates_through_actual_gas_feedback() {
    let a=specimen(0.2,8000.0);let b=specimen(0.2,8000.0);let heavy=specimen(0.6,8000.0);let i=input();
    let baseline=simulate(&a,&b,&i).unwrap();let changed=simulate(&a,&heavy,&i).unwrap();
    assert_eq!(baseline.frequencies[0],changed.frequencies[0]);
    for (x,y) in baseline.frequencies[1].iter().zip(&changed.frequencies[1]) {
        assert!((x/y-3.0_f64.sqrt()).abs()<1e-8);
    }
    assert!((changed.masses[1]/baseline.masses[1]-3.0).abs()<1e-12);
    assert!(baseline.rows.iter().any(|r|r.b.abs()>1e-12)); // undriven wall moves
    assert!(baseline.rows.iter().zip(&changed.rows).any(|(a,b)|(a.a-b.a).abs()>1e-12)); // backreaction
    for r in baseline.rows.iter().chain(&changed.rows) {
        assert!(r.residual.abs()<1e-10);
        assert!(r.volume_fraction.abs()<=i.max_volume_fraction*(1.0+1e-12));
        assert!(r.a.abs()<=i.max_displacement && r.b.abs()<=i.max_displacement);
    }
    assert_eq!(baseline.rows,simulate(&a,&b,&i).unwrap().rows);
}

#[test]
fn gas_mixture_only_change_alters_pressure_and_both_source_wall_histories() {
    let a=specimen(0.2,8000.0);let b=specimen(0.3,6000.0);
    let air=input();let co2=Input::parse(include_str!("co2.txt")).unwrap();
    let before=simulate(&a,&b,&air).unwrap();let after=simulate(&a,&b,&co2).unwrap();
    assert_eq!(before.masses,after.masses);assert_eq!(before.frequencies,after.frequencies);
    assert_ne!(air.gas.density,co2.gas.density);assert_ne!(air.gas.sound_speed,co2.gas.sound_speed);
    assert!(before.rows.iter().zip(&after.rows).any(|(a,b)|(a.pressure-b.pressure).abs()>1e-6));
    assert!(before.rows.iter().zip(&after.rows).any(|(a,b)|(a.a-b.a).abs()>1e-12 && (a.b-b.b).abs()>1e-12));
}

#[test]
fn shared_gas_parser_and_physical_limits_refuse_without_partial_output() {
    let text=include_str!("air.txt");
    for bad in [text.replace("O2 0.21","O2 21"),text.replace("2400 240","2400 2401"),
        text.replace("0.0001 0.01","0.0001 0.5"),format!("{text} extra"),text.replace("N2","unknown")] {
        assert!(Input::parse(&bad).is_err());
    }
    let legacy=text.replace("mixture 293.15 101325 mole 2 N2 0.79 O2 0.21","293.15 101325 0");
    assert_eq!(Input::parse(&legacy).unwrap().gas,GasState::try_new_moist_air(293.15,101325.0,0.0).unwrap());
    let a=specimen(0.2,8000.0);let mut i=input();i.temperature=500.0;
    assert!(simulate(&a,&a,&i).is_err());
    let mut i=input();i.force=-1e10;assert!(simulate(&a,&a,&i).is_err());
}
