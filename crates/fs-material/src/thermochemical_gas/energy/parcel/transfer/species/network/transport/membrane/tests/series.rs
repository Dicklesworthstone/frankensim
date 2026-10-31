//! G0/G1/G3: resistance-in-series oracles and shared conservative transport.
//! Synthetic coefficients establish numerical behavior, not polymer validation.
use super::*;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};

fn layer<'a>(thickness_m: f64, materials: &'a [MembranePermeabilityMaterial<'a>])
    -> MembraneSeriesElement<'a>
{
    MembraneSeriesElement::Layer { thickness_m, materials }
}
fn path<'a>(elements: &'a [MembraneSeriesElement<'a>]) -> MembranePath<'a> {
    MembranePath::Series { area_m2: geometry().area_m2, elements }
}
fn bounded_card(value: f64, pressures: [f64; 2]) -> MaterialCard {
    source(PropertyValue::Scalar { value, dims: MOLAR_PERMEABILITY_DIMS },
        MOLAR_PERMEABILITY_DIMS, InterpolationPolicy::ConstantWithinValidity, pressures)
}
fn named_card(name: &str, value: f64) -> MaterialCard {
    let original = card(value);
    let mut id = original.id().clone();
    id.chemistry = name.into();
    MaterialCard::assemble(id, original.claims().clone(), vec![]).unwrap()
}
fn interface(value: f64, dims: Dims, left: &MaterialCard, right: &MaterialCard)
    -> InterfaceSystemCard
{
    let mut claim = card(value).claims_for(MOLAR_PERMEABILITY_PROPERTY)[0].1.clone();
    claim.key = PropertyKey::new(MOLAR_INTERFACE_PERMEANCE_PROPERTY, dims);
    claim.value = PropertyValue::Scalar { value, dims };
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec { material: left.id().clone(), texture_frame: "left test surface".into() },
        SurfaceSpec { material: right.id().clone(), texture_frame: "right test surface".into() },
        SystemContext { medium: "direct dense-layer contact".into(), third_body: None,
            environment: "synthetic ideal dilute permeants".into(), history: "unchanged".into() },
        claims, vec![]).unwrap()
}
fn interface_bindings<'a>(cards: [&'a InterfaceSystemCard; 2], q: &'a QueryPoint)
    -> [MembraneInterfaceMaterial<'a>; 2]
{
    ["A", "B"].into_iter().zip(cards).map(|(name, interface)| MembraneInterfaceMaterial {
        species: SpeciesId::new(name).unwrap(), interface, context: q,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }).collect::<Vec<_>>().try_into().unwrap()
}

#[test]
fn heterogeneous_layers_follow_harmonic_resistance_and_species_counterflow() {
    let (left, right) = (gas(0.8, 100000.0), gas(0.2, 100000.0));
    let (a, b, c, d, q) = (card(2e-14), card(5e-14), card(6e-14), card(1e-14), point(400.0,100000.0));
    let first = bindings([&a, &b], &q);
    let second = bindings([&c, &d], &q);
    let layers = [layer(1e-4, &first), layer(2e-4, &second)];
    let flux = left.membrane_path_flux_to(&right, path(&layers), 0.0, || false).unwrap();
    let expected = [1.0/(1e-4/2e-14 + 2e-4/6e-14), 1.0/(1e-4/5e-14 + 2e-4/1e-14)];
    for (i, row) in flux.species.iter().enumerate() {
        near(row.permeance.value_si(), expected[i]);
        let signed_difference = if i == 0 { 60000.0 } else { -60000.0 };
        near(row.molar_flux_mol_m2_s, expected[i]*signed_difference);
        for part in &row.elements {
            near(part.permeance.value_si()*(part.partial_pressures_pa[0]-part.partial_pressures_pa[1]),
                row.molar_flux_mol_m2_s);
        }
        assert_eq!(row.elements[0].partial_pressures_pa[1], row.elements[1].partial_pressures_pa[0]);
    }
    let reversed_layers = [layers[1], layers[0]];
    let reordered = left.membrane_path_flux_to(&right, path(&reversed_layers), 0.0, || false).unwrap();
    let reverse = right.membrane_path_flux_to(&left, path(&layers), 0.0, || false).unwrap();
    for i in 0..2 {
        near(reordered.species[i].molar_flux_mol_m2_s, flux.species[i].molar_flux_mol_m2_s);
        near(reverse.species[i].molar_flux_mol_m2_s, -flux.species[i].molar_flux_mol_m2_s);
    }
    assert_ne!(flux.interface_total_pressures_pa[1], flux.face_pressure_pa);
}

#[test]
fn an_explicit_interface_adds_resistance_and_retains_its_ordered_source() {
    let (left, right) = (gas(0.8,100000.0), gas(0.2,100000.0));
    let (a,b,c,d,q) = (named_card("layer one/A",2e-14), named_card("layer one/B",5e-14),
        named_card("layer two/A",6e-14), named_card("layer two/B",1e-14), point(400.0,100000.0));
    let (ia,ib) = (interface(1e-11,MOLAR_PERMEANCE_DIMS,&a,&c),
        interface(3e-11,MOLAR_PERMEANCE_DIMS,&b,&d));
    let first = bindings([&a,&b],&q);
    let second = bindings([&c,&d],&q);
    let interfaces = interface_bindings([&ia,&ib],&q);
    let elements = [layer(1e-4,&first), MembraneSeriesElement::Interface {materials:&interfaces},
        layer(2e-4,&second)];
    let flux = left.membrane_path_flux_to(&right,path(&elements),0.0,||false).unwrap();
    near(flux.species[0].permeance.value_si(),1.0/(1e-4/2e-14+1.0/1e-11+2e-4/6e-14));
    near(flux.species[1].permeance.value_si(),1.0/(1e-4/5e-14+1.0/3e-11+2e-4/1e-14));
    for (row, original) in flux.species.iter().zip([&ia,&ib]) {
        let MembraneElementLaw::Interface {material} = &row.elements[1].law else {panic!("interface source missing")};
        assert_eq!(material.card_identity(),original.content_hash());
        assert_eq!(material.surface_materials(),&[original.surface_a().material.clone(),original.surface_b().material.clone()]);
        assert_eq!(material.property(MOLAR_INTERFACE_PERMEANCE_PROPERTY).unwrap().requirement().quantity(),
            QuantitySpec::dimensional(MOLAR_PERMEANCE_DIMS));
        for element in &row.elements {
            near(element.permeance.value_si()*(element.partial_pressures_pa[0]-element.partial_pressures_pa[1]),
                row.molar_flux_mol_m2_s);
        }
    }
    let wrong_order = interface(1e-11,MOLAR_PERMEANCE_DIMS,&c,&a);
    let wrong_bindings = interface_bindings([&wrong_order,&ib],&q);
    let wrong = [elements[0],MembraneSeriesElement::Interface {materials:&wrong_bindings},elements[2]];
    assert!(left.membrane_path_flux_to(&right,path(&wrong),0.0,||false).is_err());
    let wrong_units = interface(1e-11,MOLAR_PERMEABILITY_DIMS,&a,&c);
    let wrong_bindings = interface_bindings([&wrong_units,&ib],&q);
    let wrong = [elements[0],MembraneSeriesElement::Interface {materials:&wrong_bindings},elements[2]];
    assert!(left.membrane_path_flux_to(&right,path(&wrong),0.0,||false).is_err());
}

#[test]
fn local_material_domains_follow_actual_interfaces_not_both_reservoirs() {
    let pure = pure("A",0.028,3.5);
    let left = ClosedGasParcel::new(pure.model().clone(),point(400.0,200000.0),0.001).unwrap();
    let right = pure;
    // Resistances 1e9 and 9e9 give p_interface = 190 kPa. Neither material
    // covers the whole reservoir range, and the first excludes mean p=150 kPa.
    let a = bounded_card(1e-13,[180000.0,210000.0]);
    let b = bounded_card(1e-4/9e9,[90000.0,195000.0]);
    let (qa,qb) = (point(400.0,195000.0),point(400.0,145000.0));
    let binding = |material,context| MembranePermeabilityMaterial {species:SpeciesId::new("A").unwrap(),
        material,context,selection:MaterialPropertySelection::SingleClaimOnly};
    let (first,second) = ([binding(&a,&qa)],[binding(&b,&qb)]);
    let elements = [layer(1e-4,&first),layer(1e-4,&second)];
    let flux = left.membrane_path_flux_to(&right,path(&elements),0.0,||false).unwrap();
    near(flux.interface_total_pressures_pa[1],190000.0);
    near(flux.species[0].molar_flux_mol_m2_s,1e-5);
    for (sample, original, pressure) in [( &flux.species[0].elements[0],&a,195000.0),
        (&flux.species[0].elements[1],&b,145000.0)] {
        let MembraneElementLaw::Layer {material,..} = &sample.law else {panic!("layer source missing")};
        assert_eq!(material.card_identity(),original.content_hash());
        near(material.query_point().iter().find(|(name,_)|name=="pressure").unwrap().1,pressure);
    }
}

#[test]
fn counterpermeation_checks_internal_pressure_outside_reservoir_envelope() {
    let (left,right,q) = (gas(0.8,100000.0),gas(0.2,100000.0),point(400.0,100000.0));
    let (a,b,c,d) = (card(1e-12),card(1e-16),card(1e-16),card(1e-12));
    let (first,second) = (bindings([&a,&b],&q),bindings([&c,&d],&q));
    let elements = [layer(1e-4,&first),layer(1e-4,&second)];
    let flux = left.membrane_path_flux_to(&right,path(&elements),0.0,||false).unwrap();
    assert!(flux.interface_total_pressures_pa[1]>159000.0);
    // Both reservoir checks and the initial mean-state query would pass.
    // The actual inner interface invalidates this source.
    let narrow = bounded_card(1e-12,[90000.0,110000.0]);
    let first = bindings([&narrow,&b],&q);
    let elements = [layer(1e-4,&first),layer(1e-4,&second)];
    assert!(left.membrane_path_flux_to(&right,path(&elements),0.0,||false).is_err());
}

#[test]
fn exact_layer_and_interface_seals_do_not_become_small_leaks() {
    let (left,right,q) = (gas(0.8,100000.0),gas(0.2,100000.0),point(400.0,100000.0));
    let (zero,open) = (card(0.0),card(2e-14));
    let (first,second) = (bindings([&zero,&open],&q),bindings([&open,&open],&q));
    let layers = [layer(1e-4,&first),layer(1e-4,&second)];
    let flux = left.membrane_path_flux_to(&right,path(&layers),0.0,||false).unwrap();
    assert_eq!(flux.species[0].molar_flux_mol_m2_s,0.0);
    assert_eq!(flux.species[0].elements[1].partial_pressures_pa,[20000.0;2]);
    assert!(flux.species[1].molar_flux_mol_m2_s<0.0);
    let (ia,ib) = (interface(0.0,MOLAR_PERMEANCE_DIMS,&open,&open),
        interface(1e-10,MOLAR_PERMEANCE_DIMS,&open,&open));
    let interfaces = interface_bindings([&ia,&ib],&q);
    let elements = [layer(1e-4,&second),MembraneSeriesElement::Interface {materials:&interfaces},
        layer(1e-4,&second)];
    let flux = left.membrane_path_flux_to(&right,path(&elements),0.0,||false).unwrap();
    assert_eq!(flux.species[0].permeance.value_si(),0.0);
    assert_eq!(flux.species[0].elements[0].partial_pressures_pa,[80000.0;2]);
    assert_eq!(flux.species[0].elements[2].partial_pressures_pa,[20000.0;2]);
    let isolated = [layers[0],layers[0]];
    assert!(matches!(left.membrane_path_flux_to(&right,path(&isolated),0.0,||false),
        Err(GasEnergyError::Invalid(message)) if message.contains("undetermined")));
}

#[test]
fn very_resistive_series_avoids_overflowing_reciprocals() {
    let (left,right,q,small) = (gas(0.8,100000.0),gas(0.2,100000.0),point(400.0,100000.0),card(1e-300));
    let materials = bindings([&small,&small],&q);
    let layers = [layer(1e8,&materials);2];
    let flux = left.membrane_path_flux_to(&right,path(&layers),0.0,||false).unwrap();
    assert!((1.0/1e-308_f64).is_finite());
    assert!((2.0/1e-308_f64).is_infinite());
    for row in &flux.species {
        assert!((row.permeance.value_si()/5e-309-1.0).abs()<1e-14);
        assert!((row.molar_flux_mol_m2_s.abs()/3e-304-1.0).abs()<1e-14);
        assert_ne!(row.mass_rate_left_to_right_kg_s,0.0);
    }
}

#[test]
fn laminate_network_spends_each_species_and_partial_enthalpy_once() {
    let parcels = BTreeMap::from([(0,pure("A",0.028,3.5)),(1,pure("B",0.004,2.5))]);
    let saved = parcels.clone();
    let endpoints = parcels.keys().map(|&id|(id,energy())).collect();
    let (a,b,q) = (card(2e-14),card(5e-14),point(400.0,100000.0));
    let materials = bindings([&a,&b],&q);
    let layers = [layer(1e-4,&materials),layer(2e-4,&materials)];
    let links = [MembraneNetworkLink {id:9,parcels:[0,1],path:path(&layers),maximum_temperature_difference_k:0.0}];
    let candidate = step_candidate_with_paths(&parcels,GasTransportPaths {membranes:&links,..Default::default()},
        0.01,&endpoints,solve(),||false).unwrap();
    let step = candidate.step();
    assert_eq!(step.exchange.fluxes.len(),2);
    assert!(step.exchange.total_mass_residual_kg.abs()<1e-13);
    assert!(step.exchange.maximum_species_mass_residual_kg<1e-13);
    assert!(step.exchange.total_energy_residual_j.abs()<1e-8);
    assert!(step.entropy_change_j_k>=-1e-10);
    for row in &step.exchange.fluxes {
        let rate = step.membranes[0].flux.species.iter().find(|s|s.species==row.transfer.species).unwrap();
        near(row.transfer.mass_kg,rate.mass_rate_left_to_right_kg_s.abs()*0.01);
        let species = &parcels[&row.transfer.donor].model().thermodynamics().components()[0];
        let h = species.evaluate(fs_qty::Temperature::new(400.0)).unwrap().properties().h().value()
            /species.molar_mass().value();
        near(row.carried_enthalpy_j,row.transfer.mass_kg*h);
    }
    assert_eq!(parcels,saved);
    assert!(matches!(step_candidate_with_paths(&parcels,GasTransportPaths {membranes:&links,..Default::default()},
        0.01,&endpoints,solve(),||true),Err(GasTransportError::Cancelled)));
    assert_eq!(parcels,saved);
}

#[test]
fn series_geometry_species_work_bounds_and_mid_sample_cancellation_are_explicit() {
    let (left,right,q,a) = (gas(0.8,100000.0),gas(0.2,100000.0),point(400.0,100000.0),card(2e-14));
    let materials = bindings([&a,&a],&q);
    assert!(left.membrane_path_flux_to(&right,path(&[]),0.0,||false).is_err());
    for thickness in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(left.membrane_path_flux_to(&right,path(&[layer(thickness,&materials)]),0.0,||false).is_err());
    }
    let maximum = [layer(1e-4,&materials);MAX_MEMBRANE_ELEMENTS];
    left.membrane_path_flux_to(&right,path(&maximum),0.0,||false).unwrap();
    let too_many = [layer(1e-4,&materials);MAX_MEMBRANE_ELEMENTS+1];
    assert!(left.membrane_path_flux_to(&right,path(&too_many),0.0,||false).is_err());
    let missing = [layer(1e-4,&materials[..1])];
    assert!(left.membrane_path_flux_to(&right,path(&missing),0.0,||false).is_err());
    let inactive = [layer(1e-4,&[])];
    assert!(left.membrane_path_flux_to(&right,MembranePath::Series {area_m2:0.0,elements:&inactive},
        0.0,||false).unwrap().species.is_empty());
    let mut polls = 0;
    assert_eq!(left.membrane_path_flux_to(&right,path(&maximum),0.0,||{polls+=1;polls==100}),
        Err(GasEnergyError::Cancelled));
}
