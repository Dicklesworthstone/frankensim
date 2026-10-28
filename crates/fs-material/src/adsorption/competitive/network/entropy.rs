//! Closed gas+adsorbate entropy with conserved species references cancelled.
use super::*;

pub(super) fn admit_catalog(parcels:&BTreeMap<u64,ClosedGasParcel>,catalog:&ThermochemicalGasModel,
    cancelled:&mut impl FnMut()->bool)->Result<()> {
    for gas in parcels.values() {
        for source in gas.model().thermodynamics().components() {
            poll(cancelled)?;
            if !catalog.thermodynamics().components().iter().any(|s|s==source) {
                return Err(invalid("joint entropy catalog must retain every original gas thermochemical source"));
            }
        }
    }
    Ok(())
}

pub(super) fn mass_balance(old:&BTreeMap<u64,ClosedGasParcel>,new:&BTreeMap<u64,ClosedGasParcel>,
    solids:&[(&AdsorptionNetworkBinding<'_>,PreparedAdsorption)],catalog:&ThermochemicalGasModel,
    cancelled:&mut impl FnMut()->bool)->Result<(f64,f64)> {
    let mut changes=Vec::with_capacity(catalog.thermodynamics().components().len());
    let mut maximum:f64=0.0;
    for source in catalog.thermodynamics().components() {
        let mut terms=Vec::with_capacity(old.len()+solids.len());
        for (&id,gas) in old {
            poll(cancelled)?;
            terms.push(new[&id].species_mass_kg(source.species())?-gas.species_mass_kg(source.species())?);
        }
        for (binding,proposal) in solids {
            poll(cancelled)?;
            if let Ok(i)=binding.solid.model.channels.binary_search_by(|law|law.species.cmp(source.species())) {
                terms.push(proposal.next.moles[i]*source.molar_mass().value()
                    -binding.solid.moles[i]*source.molar_mass().value());
            }
        }
        let change=sum(terms)?; maximum=maximum.max(change.abs()); changes.push(change);
    }
    Ok((sum(changes)?,maximum))
}

pub(super) fn change(old:&BTreeMap<u64,ClosedGasParcel>,new:&BTreeMap<u64,ClosedGasParcel>,
    solids:&[(&AdsorptionNetworkBinding<'_>,PreparedAdsorption)],catalog:&ThermochemicalGasModel,
    exchange:&SpeciesNetworkStep,cancelled:&mut impl FnMut()->bool)->Result<f64> {
    let mut references=BTreeMap::new();
    for (binding,_) in solids {
        for law in &binding.solid.model.channels {
            references.entry(law.species.clone()).or_insert(law.source.reference_temperature_k);
        }
    }
    for gas in old.values() {
        for source in gas.model().thermodynamics().components() {
            references.entry(source.species().clone()).or_insert(gas.state().state().temperature);
        }
    }
    // Use intended stream increments, not rounded inventories multiplied by a
    // standard-state entropy constant. Their actual mass defects are separate.
    let mut changes:BTreeMap<(u64,SpeciesId),Vec<f64>>=BTreeMap::new();
    for flux in &exchange.fluxes {
        let row=&flux.transfer;
        changes.entry((row.donor,row.species.clone())).or_default().push(-row.mass_kg);
        changes.entry((row.receiver,row.species.clone())).or_default().push(row.mass_kg);
    }
    for row in &exchange.boundaries {
        changes.entry((row.parcel,row.transfer.species.clone())).or_default().push(row.transfer.net_mass_into_gas_kg);
    }
    let mut terms=Vec::new();
    for (binding,proposal) in solids {
        poll(cancelled)?;
        let old=binding.solid; let next=&proposal.next;
        terms.push(old.dry_mass_kg*old.model.dry().specific_entropy_change(old.temperature_k,next.temperature_k)?);
        terms.push(site_mixing_change(old,next)?);
        for (i,law) in old.model.channels.iter().enumerate() {
            poll(cancelled)?;
            let dn=next.moles[i]-old.moles[i]; let reference=law.source.reference_temperature_k;
            let common=references[&law.species];
            let standard=-law.reference_heat_j_mol/reference+R*det::ln(law.reference_affinity_pa_inverse)
                +law.gas_ds(common,reference)?;
            terms.push(next.moles[i]*law.bound_ds(old.temperature_k,next.temperature_k)?
                +dn*law.bound_ds(reference,old.temperature_k)?+dn*standard);
        }
    }
    for (&id,gas) in old {
        let next=&new[&id]; let (t0,t1)=(gas.state().state().temperature,next.state().state().temperature);
        for source in catalog.thermodynamics().components() {
            poll(cancelled)?;
            let mass=source.molar_mass().value();
            let n0=gas.species_mass_kg(source.species())?/mass;
            let n1=next.species_mass_kg(source.species())?/mass;
            if n0==0.0 && n1==0.0 { continue; }
            let reference=*references.get(source.species()).ok_or(invalid("joint entropy lost a species reference"))?;
            let delta=sum(changes.get(&(id,source.species().clone())).into_iter().flatten().copied())?/mass;
            let ds=|a,b|source.molar_entropy_change(Temperature::new(a),Temperature::new(b))
                .map(|s|s.value()).map_err(gas_failure);
            let log_p=|n:f64,t:f64|det::ln(n)+det::ln(R)+det::ln(t)-det::ln(gas.volume_m3());
            let pressure=if n0==0.0 { n1*log_p(n1,t1) } else if n1==0.0 { -n0*log_p(n0,t0) }
                else { delta*log_p(n0,t0)+n1*(log_ratio(n1,n0)+log_ratio(t1,t0)) };
            terms.push(n1*ds(t0,t1)?+delta*ds(reference,t0)?-R*pressure);
        }
    }
    sum(terms)
}
