//! Unequal physical reservoir temperatures reach subsequent source-bound impacts.
use super::*;
use fs_couple::render::plate::impact::material::heating::feedback::resolve_contact_at_temperatures;
use fs_conduction::interface::{InterfaceResistance,AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS};
use fs_conduction::interface::lumped_pair::{advance_contact_pair,ContactPairReport};
use fs_conduction::lumped::LumpedThermalTransport;
use fs_matdb::{SelectionPolicy,SurfaceSpec,SystemContext};
use fs_qty::QuantitySpec;

fn pair_cards<'a>(a:&'a MaterialCard,b:&'a MaterialCard,i:&'a InterfaceSystemCard)->PointContactMaterialCards<'a> {
    PointContactMaterialCards {surface_a:a,surface_b:b,interface:i,
        bulk_selection:[MaterialPropertySelection::SingleClaimOnly,MaterialPropertySelection::SingleClaimOnly],
        interface_selection:MaterialPropertySelection::SingleClaimOnly,
        model_selection:NormalContactModelSelection::SingleSupported,authority:InputAuthority::SyntheticFixture}
}
fn pair_interface(a:&MaterialCard,b:&MaterialCard)->InterfaceSystemCard {
    crate::common::bind(a,b,Some(0.4),0.0).normal_state().elastic().interface().card().clone()
}
fn pair<'a>(cards:&PointContactMaterialCards<'_>,curves:[&'a EquilibriumEnthalpyPhaseCurve;2],
    enthalpies:[f64;2])->HeatedImpactSystem<'a> {
    let (mechanics,_)=heated(cards,curves[0],0.0).into_parts();
    let stores=["left","right"].into_iter().enumerate().map(|(i,key)|CaloricReservoir::new(
        key.into(),curves[i],curves[i].state_at_specific_enthalpy(enthalpies[i]).unwrap(),
        0.001,1e-6,0.0006).unwrap()).collect();
    let mut split=route("touch");split.fraction_to_a=0.7;
    HeatedImpactSystem::new(mechanics,stores,vec![split],1e-11).unwrap()
}
fn refresh(c:&MaterialContact,stores:&[CaloricReservoir<'_>],cards:&PointContactMaterialCards<'_>)
    ->Result<Option<SourcePointContact>,ImpactError> {
    // Explicit fixture assumption: side A supplies the interface-law trace.
    // The bulk on side B must still use ITS OWN accepted temperature.
    let source=resolve_contact_at_temperatures(c,cards,&point(),[&stores[0],&stores[1]],
        stores[0].state().temperature_k())?;
    let elastic=source.model().normal_state().elastic();
    for (s,bulk) in stores.iter().zip([elastic.surface_a_state(),elastic.surface_b_state()]) {
        assert_eq!(bulk.resolved().query_point(),&[("T".into(),s.state().temperature_k())]);
    }
    Ok(Some(source))
}

#[test]
fn changing_only_one_caloric_card_changes_a_later_impact_without_averaging_the_pair() {
    let left_fast=material(1e-6);let left_slow=material(2e-6);let right=material(2e-6);
    let ia=pair_interface(&left_fast,&right);let ib=pair_interface(&left_slow,&right);
    let a=pair_cards(&left_fast,&right,&ia);let b=pair_cards(&left_slow,&right,&ib);
    let ca=caloric(&left_fast,1e-6);let cb=caloric(&left_slow,2e-6);let cr=caloric(&right,2e-6);
    let mut fast=pair(&a,[&ca,&cr],[0.0;2]);let mut slow=pair(&b,[&cb,&cr],[0.0;2]);
    let gate=CancelGate::new_clock_free();let initial=fast.mechanics().stored_energy_j();
    let mut work=0.0;let mut peaks=[[0.0_f64;2];2];let mut refreshed=0;let mut held=0;
    for sample in 0..6000 {
        let forces=if (2000..2500).contains(&sample){[1.0,-1.0]}else{[0.0;2]};
        let f=fast.step_with_material_feedback(&forces,&gate,|c,s|refresh(c,s,&a),insulated).unwrap();
        slow.step_with_material_feedback(&forces,&gate,|c,s|refresh(c,s,&b),insulated).unwrap();
        refreshed+=f.refreshed_contacts;held+=f.held_loaded_contacts;
        work+=f.thermal.impact.mechanical.supplied_work_j;
        let episode=usize::from(sample>=2000);
        peaks[0][episode]=peaks[0][episode].max(fast.mechanics().material_contact_observation(0).unwrap().normal_force_n);
        peaks[1][episode]=peaks[1][episode].max(slow.mechanics().material_contact_observation(0).unwrap().normal_force_n);
        if sample<2000 {assert_eq!(fast.mechanics().state(),slow.mechanics().state());}
        assert!(f.thermal.impact.combined_energy_residual_j.abs()<1e-9);
    }
    assert_eq!(peaks[0][0],peaks[1][0]);assert!(peaks[0][1]>1.002*peaks[1][1]);
    assert!(refreshed>=2 && held>0);
    assert!(fast.reservoirs()[0].state().temperature_k()>fast.reservoirs()[1].state().temperature_k());
    assert!(fast.reservoirs()[0].state().temperature_k()>slow.reservoirs()[0].state().temperature_k());
    near(fast.mechanics().stored_energy_j()+fast.reservoirs().iter().map(|r|r.storage_change_j()).sum::<f64>()-work,initial,1e-9);
    assert!(resolve_isothermal_contact(&fast.mechanics().material_contacts()[0],&a,&point(),
        [&fast.reservoirs()[0],&fast.reservoirs()[1]]).is_err());
}

fn resistance(cards:&PointContactMaterialCards<'_>)->InterfaceResistance {
    let dims=AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,dims),
        value:PropertyValue::Scalar {value:1e-4,dims},validity:ValidityDomain::unconstrained(),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance {source:"synthetic finite thermal interface".into(),
            license:"CC0-1.0".into(),artifact:None}}).unwrap();
    let interface=InterfaceSystemCard::assemble(
        SurfaceSpec {material:cards.surface_a.id().clone(),texture_frame:"surface-a".into()},
        SurfaceSpec {material:cards.surface_b.id().clone(),texture_frame:"surface-b".into()},
        SystemContext {medium:"dry".into(),third_body:None,environment:"test".into(),history:"initial".into()},
        claims,vec![]).unwrap();
    InterfaceResistance::from_card("touch",&interface,&point(),SelectionPolicy::SingleClaimOnly).unwrap()
}

#[test]
fn unequal_temperature_feedback_composes_with_actual_contact_area_and_finite_heat_exchange() {
    let a=material(0.02);let b=material(0.04);let i=pair_interface(&a,&b);let cards=pair_cards(&a,&b,&i);
    let ca=caloric(&a,0.02);let cb=caloric(&b,0.04);let interface=resistance(&cards);
    let transport=LumpedThermalTransport::try_declared(100.0,0.0).unwrap();
    context(|gate,cx| {
        let mut s=pair(&cards,[&ca,&cb],[0.2,0.0]); // 310 K against 300 K.
        let initial=s.mechanics().stored_energy_j();let mut transferred=0.0;let mut active=0;let mut dry=0;
        for _ in 0..1800 {
            let f=s.step_with_material_feedback(&[0.0;2],gate,|c,r|refresh(c,r,&cards),|trial| {
                let r=trial.reservoirs();let area=trial.contact_patch_area_m2("touch")?;
                let result=advance_contact_pair(cx,[r[0].thermal,r[1].thermal],
                    [&transport,&transport],&interface,area,1e-13)
                    .map_err(|e|ImpactError::Owner(e.to_string()))?;
                Ok(CoupledThermalStep {transfers:Vec::from(result),boundary_heat_j:0.0})
            }).unwrap();
            let report:&ContactPairReport=&f.thermal.transfers[0].report;
            assert_eq!(report.interface,interface);
            assert_eq!(f.thermal.boundary_heat_j,0.0);
            assert_eq!(f.thermal.transfers[0].external_heat_j,-f.thermal.transfers[1].external_heat_j);
            if report.area_m2>0.0 {active+=1;transferred+=f.thermal.transfers[1].external_heat_j;}
            else {dry+=1;assert_eq!(f.thermal.transfers[0].external_heat_j,0.0);}
            assert!(f.thermal.impact.combined_energy_residual_j.abs()<1e-9);
        }
        assert!(active>0 && dry>0 && transferred>1e-7);
        assert!(s.reservoirs()[0].state().temperature_k()<310.0);
        assert!(s.reservoirs()[1].state().temperature_k()>300.0);
        near(s.mechanics().stored_energy_j()+s.reservoirs().iter().map(|r|r.storage_change_j()).sum::<f64>(),initial,1e-9);
        near(s.reservoirs().iter().map(|r|r.boundary_heat_j()).sum(),0.0,1e-15);
    });
}

#[test]
fn unequal_source_refresh_rolls_back_with_transport_solver_or_late_cancellation() {
    let a=material(0.02);let b=material(0.04);let i=pair_interface(&a,&b);let cards=pair_cards(&a,&b,&i);
    let ca=caloric(&a,0.02);let cb=caloric(&b,0.04);let mut s=pair(&cards,[&ca,&cb],[0.2,0.0]);
    let gate=CancelGate::new_clock_free();let sources=s.mechanics().material_contacts().to_vec();
    let motion=s.mechanics().state().to_vec();let stores=s.reservoirs().to_vec();
    assert!(s.step_with_material_feedback(&[0.0;2],&gate,|c,r|refresh(c,r,&cards),
        |_|Err::<CoupledThermalStep<()>,_>(ImpactError::Invalid("finite thermal owner refused"))).is_err());
    assert_eq!(s.mechanics().material_contacts(),sources);
    s.set_iteration_limit(0).unwrap();
    assert!(s.step_with_material_feedback(&[0.0;2],&gate,|c,r|refresh(c,r,&cards),insulated).is_err());
    s.set_iteration_limit(50).unwrap();let late=CancelGate::new_clock_free();
    assert!(s.step_with_material_feedback(&[0.0;2],&late,|c,r|refresh(c,r,&cards),|trial| {
        let result=insulated(trial);late.request();result
    }).is_err());
    assert_eq!(s.mechanics().material_contacts(),sources);assert_eq!(s.mechanics().state(),motion);
    assert_eq!(s.reservoirs(),stores);assert_eq!(s.mechanics().samples(),0);
    let f=s.step_with_material_feedback(&[0.0;2],&gate,|c,r|refresh(c,r,&cards),insulated).unwrap();
    assert_eq!(f.refreshed_contacts,1);assert_eq!(s.mechanics().samples(),1);
    let elastic=s.mechanics().material_contacts()[0].source().model().normal_state().elastic();
    assert_eq!(elastic.surface_a_state().resolved().query_point(),&[("T".into(),310.0)]);
    assert_eq!(elastic.surface_b_state().resolved().query_point(),&[("T".into(),300.0)]);
}

#[test]
fn independent_temperature_binding_retains_strict_source_phase_density_and_context_admission() {
    let a=material(0.02);let b=material(0.04);let i=pair_interface(&a,&b);let cards=pair_cards(&a,&b,&i);
    let ca=caloric(&a,0.02);let cb=caloric(&b,0.04);let s=pair(&cards,[&ca,&cb],[0.2,0.0]);
    let contact=&s.mechanics().material_contacts()[0];let r=s.reservoirs();
    let resolve=|stores,reference:&QueryPoint,t|resolve_contact_at_temperatures(contact,&cards,reference,stores,t);
    let law=resolve([&r[0],&r[1]],&point(),305.0).unwrap();
    let e=law.model().normal_state().elastic();
    assert_eq!(e.state_coordinate("T"),Some(305.0));
    near(e.reduced_modulus_pa(),1.0/(0.9375/3.25e7+0.9375/2.5e7),1e-7);
    assert!(resolve([&r[1],&r[0]],&point(),305.0).is_err());
    for t in [0.0,-1.0,f64::NAN,f64::INFINITY,330.0] {assert!(resolve([&r[0],&r[1]],&point(),t).is_err());}
    let extra=point().with("p",1e5).unwrap();assert!(resolve([&r[0],&r[1]],&extra,305.0).is_err());
    let typed=point().with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),300.0).unwrap();
    assert!(resolve([&r[0],&r[1]],&typed,305.0).is_err());
    let cold=CaloricReservoir::new("left".into(),&ca,ca.state_at_specific_enthalpy(0.0).unwrap(),0.001,1e-6,0.0006).unwrap();
    assert_eq!(resolve([&cold,&r[1]],&point(),300.0).unwrap(),
        resolve_isothermal_contact(contact,&cards,&point(),[&cold,&r[1]]).unwrap());
    for density in [900.0,1000.0] {
        let curve=EquilibriumEnthalpyPhaseCurve::try_new(a.content_hash(),vec![
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:310.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:density},
            EnthalpyPhaseKnot {specific_enthalpy_j_kg:1.0,temperature_k:310.0,liquid_mass_fraction:1.0,bulk_density_kg_m3:density},
        ]).unwrap();
        // Separately reject inconsistent solid density and a liquid fraction.
        let h=if density==900.0 {0.0}else{0.5};
        let bad=CaloricReservoir::new("bad".into(),&curve,curve.state_at_specific_enthalpy(h).unwrap(),0.001,1e-6,0.0006).unwrap();
        // Called directly: the `resolve` closure fixes one borrow lifetime for its stores.
        assert!(resolve_contact_at_temperatures(contact,&cards,&point(),[&bad,&r[1]],305.0).is_err());
    }
}

#[test]
fn source_dependent_contact_heat_reaches_both_native_reservoirs_with_closed_feedback() {
    use fs_conduction::interface::lumped_pair::MaterialContactLaw;
    use fs_conduction::interface::temperature::ContactIterationConfig;
    use fs_material::thermal::interface::ThermalInterfacePoint;

    let a = material(0.02);
    let b = material(0.04);
    let mechanical_interface = pair_interface(&a, &b);
    let cards = pair_cards(&a, &b, &mechanical_interface);
    let ca = caloric(&a, 0.02);
    let cb = caloric(&b, 0.04);
    let kelvin = QuantitySpec::dimensional(Temperature::DIMS);
    let dims = AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
            value: PropertyValue::Curve {
                abscissa: "meanT".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 1e-4), (320.0, 3e-4)],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("aT", kelvin, 300.0, 320.0)
                .with_quantity("bT", kelvin, 300.0, 320.0)
                .with_quantity("meanT", kelvin, 300.0, 320.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic temperature-dependent contact resistance, not measured material"
                    .into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let thermal_interface = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.id().clone(),
            texture_frame: "surface-a".into(),
        },
        SurfaceSpec {
            material: b.id().clone(),
            texture_frame: "surface-b".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "test".into(),
            history: "initial".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let conditions = QueryPoint::new()
        .with_quantity("aT", kelvin, 310.0)
        .unwrap()
        .with_quantity("bT", kelvin, 300.0)
        .unwrap()
        .with_quantity("meanT", kelvin, 305.0)
        .unwrap();
    let binding = ThermalInterfacePoint {
        interface: &thermal_interface,
        materials: [&a, &b],
        context: &conditions,
        temperature_axes: ["aT", "bT"],
        liquid_fraction_axes: [None, None],
        mean_temperature_axis: Some("meanT"),
    };
    let transport = LumpedThermalTransport::try_declared(100.0, 0.0).unwrap();
    let iteration = ContactIterationConfig {
        max_iterations: 64,
        relaxation: 0.7,
        residual_tolerance_j: 1e-13,
    };

    context(|gate, cx| {
        let mut system = pair(&cards, [&ca, &cb], [0.2, 0.0]);
        let law = MaterialContactLaw::new(
            binding,
            [
                system.reservoirs()[0].state(),
                system.reservoirs()[1].state(),
            ],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let reference = law
            .reference()
            .property(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY)
            .unwrap();
        let selected = reference.answer().receipt.selected;
        let initial_resistance = reference.value_si();
        near(initial_resistance, 1.5e-4, 1e-18);
        let initial_energy = system.mechanics().stored_energy_j();
        let mut transferred = 0.0;
        let mut active = 0;
        let mut changed_source_value = false;

        for _ in 0..1800 {
            let frame = system
                .step_with_material_feedback(
                    &[0.0; 2],
                    gate,
                    |contact, reservoirs| refresh(contact, reservoirs, &cards),
                    |trial| {
                        let reservoirs = trial.reservoirs();
                        let area = trial.contact_patch_area_m2("touch")?;
                        let transfers = law
                            .advance(
                                cx,
                                [reservoirs[0].thermal, reservoirs[1].thermal],
                                [&transport, &transport],
                                area,
                                1e-14,
                                iteration,
                            )
                            .map_err(|error| ImpactError::Owner(error.to_string()))?;
                        Ok(CoupledThermalStep {
                            transfers: Vec::from(transfers),
                            boundary_heat_j: 0.0,
                        })
                    },
                )
                .unwrap();

            let transfers = &frame.thermal.transfers;
            let report = &transfers[0].report;
            assert_eq!(frame.thermal.boundary_heat_j, 0.0);
            assert_eq!(transfers[0].external_heat_j, -transfers[1].external_heat_j);
            assert!(frame.thermal.impact.combined_energy_residual_j.abs() < 1e-9);
            assert!(report.physical_contact_residual_j.abs() <= iteration.residual_tolerance_j);
            assert!(report.endpoint_maximum_biot.is_finite());
            assert!(report.endpoint_maximum_biot <= 0.1);
            assert!(report.iterations <= iteration.max_iterations);

            let temperatures = [
                system.reservoirs()[0].state().temperature_k(),
                system.reservoirs()[1].state().temperature_k(),
            ];
            assert_eq!(
                report.endpoint.state_coordinate("aT"),
                Some(temperatures[0])
            );
            assert_eq!(
                report.endpoint.state_coordinate("bT"),
                Some(temperatures[1])
            );
            let mean = 0.5 * temperatures[0] + 0.5 * temperatures[1];
            assert_eq!(report.endpoint.state_coordinate("meanT"), Some(mean));
            let endpoint = report
                .endpoint
                .property(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY)
                .unwrap();
            assert_eq!(endpoint.answer().receipt.selected, selected);
            near(
                endpoint.value_si(),
                1e-4 * (1.0 + 0.1 * (mean - 300.0)),
                1e-16,
            );
            changed_source_value |= (endpoint.value_si() - initial_resistance).abs() > 1e-12;
            if report.applied.area_m2 > 0.0 {
                active += 1;
                transferred += transfers[1].external_heat_j;
            } else {
                assert_eq!(transfers[0].external_heat_j, 0.0);
            }
        }

        assert!(active > 0 && transferred > 1e-7);
        assert!(changed_source_value);
        assert!(system.reservoirs()[0].state().temperature_k() < 310.0);
        assert!(system.reservoirs()[1].state().temperature_k() > 300.0);
        near(
            system.mechanics().stored_energy_j()
                + system
                    .reservoirs()
                    .iter()
                    .map(|r| r.storage_change_j())
                    .sum::<f64>(),
            initial_energy,
            1e-9,
        );
        near(
            system
                .reservoirs()
                .iter()
                .map(|r| r.boundary_heat_j())
                .sum(),
            0.0,
            1e-15,
        );

        // Refusal after a real mechanical trial must leave both accepted histories intact.
        let motion = system.mechanics().state().to_vec();
        let samples = system.mechanics().samples();
        let sources = system.mechanics().material_contacts().to_vec();
        let stores = system.reservoirs().to_vec();
        let mut thermal_called = false;
        let refused = system.step_with_coupled_transport(&[0.0; 2], gate, |trial| {
            thermal_called = true;
            let reservoirs = trial.reservoirs();
            let area = trial.contact_patch_area_m2("touch")?;
            let transfers = law
                .advance(
                    cx,
                    [reservoirs[0].thermal, reservoirs[1].thermal],
                    [&transport, &transport],
                    area,
                    1e-14,
                    ContactIterationConfig {
                        max_iterations: 0,
                        ..iteration
                    },
                )
                .map_err(|error| ImpactError::Owner(error.to_string()))?;
            Ok(CoupledThermalStep {
                transfers: Vec::from(transfers),
                boundary_heat_j: 0.0,
            })
        });
        assert!(thermal_called && refused.is_err());
        assert_eq!(system.mechanics().state(), motion);
        assert_eq!(system.mechanics().samples(), samples);
        assert_eq!(system.mechanics().material_contacts(), sources);
        assert_eq!(system.reservoirs(), stores);
    });
}
