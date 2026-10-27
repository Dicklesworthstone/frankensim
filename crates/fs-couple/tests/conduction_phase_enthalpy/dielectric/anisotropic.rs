//! G1/G3/G4: the original spatial charge/enthalpy owner uses sourced epsilon axes.
use super::*;
use fs_material::field_storage::planar::{
    PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, PlanarDielectricLaw,
};

fn material(varying_second: bool) -> MaterialCard {
    let base = super::material(3., 1., false);
    let mut claims = base.claims().clone();
    for (i, name) in PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES.into_iter().enumerate() {
        let value = [2., 4.][i];
        let varying = i == 1 && varying_second;
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, ELECTRIC_PERMITTIVITY_DIMS),
            value: if varying {
                PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(300., value), (350., value+1.)],
                    dims: ELECTRIC_PERMITTIVITY_DIMS,
                }
            } else {
                PropertyValue::Scalar { value, dims: ELECTRIC_PERMITTIVITY_DIMS }
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T", QuantitySpec::dimensional(Temperature::DIMS), 300., 350.,
            ),
            interpolation: if varying { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored constant principal dielectric fixture".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn conductor(card: &MaterialCard, context: &QueryPoint, tensor: bool) -> ResolvedPlanarConductor {
    let original = super::conductor(card, context, tensor);
    let mut requirements: Vec<_> = original.material().properties().iter()
        .map(|property| property.requirement().clone()).collect();
    requirements.extend(PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES.map(|name| {
        ScalarPropertyRequirement::try_new(
            name, ELECTRIC_PERMITTIVITY_DIMS, ScalarAdmissibility::StrictlyPositive,
        ).unwrap()
    }));
    let bundle = resolve_material_state_point(
        card, context, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    original.with_material(bundle).unwrap()
}
fn law(transverse: bool, geometry: u8) -> PlanarDielectricLaw {
    PlanarDielectricLaw::Orthotropic(
        PlanarElectricalFrame::new(
            ContentHash([9;32]), ContentHash([geometry;32]),
            core::f64::consts::FRAC_PI_4
                + if transverse { core::f64::consts::FRAC_PI_2 } else { 0. },
        ).unwrap(),
    )
}
fn field(card: &MaterialCard, context: &QueryPoint, tensor: bool, transverse: bool)
    -> PlanarLeakyDielectric<ResolvedPlanarConductor>
{
    let conduction = PlanarConduction::new(
        mesh(), vec![conductor(card,context,tensor);4], limits(), &mut || false,
    ).unwrap();
    PlanarLeakyDielectric::new_with_dielectric_laws(
        conduction, vec![law(transverse,2);4], vec![0.,8.,8.,0.,4.], &[0,3], 0.,
        limits(), electric_controls(), &mut || false,
    ).unwrap()
}

// Independent scalar elimination on the latent plateau. Physical R(f)=1+3f,
// Q=epsilon_parallel*(V0²-V1²)/2 and the midpoint charge equation determine h1.
// Bisection stays inside the authored f in [0,1] source interval.
fn reference(v0: f64, h0: f64, epsilon: f64) -> (f64, f64) {
    let voltage = |h: f64| {
        let resistance = 1.+3.*(h-20.)/100.;
        v0*(2.*epsilon*resistance-1.)/(2.*epsilon*resistance+1.)
    };
    let residual = |h: f64| h-h0-0.5*epsilon*(v0*v0-voltage(h).powi(2));
    let (mut low,mut high)=(h0,120_f64.min(h0+0.5*epsilon*v0*v0));
    assert!(residual(low)<0. && residual(high)>0.);
    for _ in 0..100 {
        let middle=0.5*(low+high);
        if residual(middle)<0. { low=middle; } else { high=middle; }
    }
    let h=0.5*(low+high);
    (voltage(h),h)
}

#[test]
fn dielectric_axes_change_charge_relaxation_and_native_phase_heat_feedback() {
    let card=material(false);
    let caloric=curve(&card,1.);
    let context=point(310.,0.);
    let cells=[binding(&card,&caloric,&context)];
    let layout=network(&cells);
    let mut retained=Vec::new();
    for tensor in [false,true] {
        for transverse in [false,true] {
            let epsilon=if transverse {4.} else {2.};
            let thermal=initial(&caloric,20.);
            let old=field(&card,&context,tensor,transverse);
            let first=layout.step_dielectric_candidate(
                &thermal,1.,&old,&[(0,0.),(3,0.)],limits(),electric_controls(),
                config(),&[25.],||false,
            ).unwrap();
            let (v1,h1)=reference(8.,20.,epsilon);
            near(first.field().voltage_v()[1],v1);
            near(first.thermal().energy_j(),h1);
            near(first.thermal().temperature_k(CellId(7)).unwrap(),310.);
            near(first.report().deposited_joule_j,h1-20.);
            near(first.report().electrical.supplied_work_j,0.);
            near(first.field().stored_energy_j(),0.5*epsilon*v1*v1);
            near(first.thermal().energy_j()+first.field().stored_energy_j(),20.+32.*epsilon);
            check_sources(old.conduction(),first.field().conduction(),310.,(h1-20.)/100.,3.);
            let next=layout.step_dielectric_candidate(
                first.thermal(),2.,first.field(),&[(0,0.),(3,0.)],
                limits(),electric_controls(),config(),&[10.],||false,
            ).unwrap();
            let (v2,h2)=reference(v1,h1,epsilon);
            near(next.field().voltage_v()[1],v2);
            near(next.thermal().energy_j(),h2);
            near(next.thermal().energy_j()+next.field().stored_energy_j(),20.+32.*epsilon);
            assert_eq!(next.field().time_s(),next.thermal().time_s());
            assert_eq!(next.field().dielectric().permittivity_laws(),&[law(transverse,2);4]);
            for (before,after) in old.dielectric().materials().iter().zip(next.field().dielectric().materials()) {
                assert_eq!(before.card_identity(),after.card_identity());
                near(after.query().axes()[FRACTION],(h2-20.)/100.);
                for (name,value) in PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES.into_iter().zip([2.,4.]) {
                    let a=before.property(name).unwrap();
                    let b=after.property(name).unwrap();
                    near(b.value_si(),value);
                    assert_eq!(a.answer().receipt.selected,b.answer().receipt.selected);
                    assert_eq!(a.source_axes(),b.source_axes());
                }
            }
            assert_eq!(old.time_s(),0.);
            near(old.stored_energy_j(),32.*epsilon);
            near(thermal.energy_j(),20.);
            retained.push(v2);
        }
    }
    assert!(retained[1]>retained[0]);
    near(retained[0],retained[2]);
    near(retained[1],retained[3]);
}

#[test]
fn independent_dielectric_frame_and_source_refusals_leave_charge_unchanged() {
    let context=point(310.,0.);
    let card=material(false);
    let build=|source:ResolvedPlanarConductor,law:PlanarDielectricLaw| {
        let conduction=PlanarConduction::new(mesh(),vec![source;4],limits(),&mut||false).unwrap();
        PlanarLeakyDielectric::new_with_dielectric_laws(
            conduction,vec![law;4],vec![0.,8.,8.,0.,4.],&[0,3],0.,
            limits(),electric_controls(),&mut||false,
        )
    };
    // Independent material frames are admitted; incompatible mesh frames refuse.
    assert!(build(conductor(&card,&context,true),law(false,3)).is_err());
    let varying=material(true);
    assert!(build(conductor(&varying,&context,true),law(false,2)).is_err());
    // A scalar permittivity in the same bundle cannot fill an absent principal axis.
    assert!(build(super::conductor(&card,&context,true),law(false,2)).is_err());
    let old=field(&card,&context,true,false);
    let charges=old.nodal_charge_c().to_vec();
    let coefficients=old.dielectric().permittivity_coefficients().to_vec();
    let mut polls=0;
    old.step_to_candidate(
        1.,old.conduction().conductors().to_vec(),&[(0,0.),(3,0.)],
        limits(),electric_controls(),&mut||{polls+=1;false},
    ).unwrap();
    let mut calls=0;
    assert!(matches!(
        old.step_to_candidate(
            1.,old.conduction().conductors().to_vec(),&[(0,0.),(3,0.)],
            limits(),electric_controls(),&mut||{calls+=1;calls>=polls},
        ),
        Err(fs_em::electrostatic::relaxation::DielectricRelaxationError::Cancelled)
    ));
    assert_eq!(old.nodal_charge_c(),charges.as_slice());
    assert_eq!(old.dielectric().permittivity_coefficients(),coefficients.as_slice());
    assert_eq!(old.time_s(),0.);
    near(old.stored_energy_j(),64.);
}
