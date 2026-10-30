//! Competitive Langmuir exchange with ONE finite site population and heat store.
//!
//! Each retained single-species source supplies its original affinity, kinetics,
//! bound calorics and thermochemical reference. All channels must explicitly
//! share one capacity, dry-skeleton chart and reference temperature. A molecule
//! occupies one site; independent pure-species capacities are NOT added.
//!
//! Frozen rates obey dn_i/dt = a_i*(C-sum(n_j))-d_i*n_i, a_i=d_i*b_i*p_i.
//! The diagonal-plus-rank-one backward-Euler solve is O(species), preserves
//! positive site inventories without clipping, and reports both gross streams.
//! Gas withdrawals share the OLD gas inventory; one original gas reconstruction
//! updates every species, temperature, pressure and transport property together.
//! One common solid temperature includes dry storage once and every bound phase.
//!
//! Shared-site mixing entropy includes exactly one vacancy species. Conserved
//! species reference entropies cancel before rounding. Finite-bath entropy and
//! inventory gates remain active. A species can cross its eventual equilibrium
//! loading during physical competitive displacement, so no independent scalar
//! equilibrium-crossing restriction is imposed. This is first-order frozen-
//! forcing integration, not an exact kinetic solution or unconditional error
//! guarantee. Near-isothermal applicability, source domains and refinement
//! remain explicit. No multilayer sites, multisite molecules, swelling or
//! lateral adsorbate interactions are inferred.
//!
//! Site balance: COMSOL, Governing Equations for the Surface Concentrations,
//! https://doc.comsol.com/6.3/doc/com.comsol.help.echem/echem_ug_chemsptrans.07.112.html
//! Standard-state thermodynamics remain those of the parent Langmuir owner.

use super::*;
use std::collections::BTreeMap;

pub mod sites;
/// Source-bound isothermal diffusion of bound species through shared vacancies.
pub mod surface;
/// Shared-snapshot gas transport and finite competitive adsorption.
pub mod network;
pub use sites::{CompetitiveSiteStep, backward_euler_sites};
pub const MAX_COMPETITIVE_SPECIES: usize = 32;

fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let (mut value, mut correction) = (0.0, 0.0);
    for term in values {
        let adjusted = finite(term - correction)?;
        let next = finite(value + adjusted)?;
        correction = finite((next - value) - adjusted)?;
        value = next;
    }
    Ok(value)
}

/// An explicitly declared common lattice of singly occupied adsorption sites.
/// Sorting makes runtime state, sums and source reports canonical by SpeciesId.
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveAdsorbent {
    channels: Vec<LangmuirAdsorbent>,
    temperature_range_k: [f64; 2],
}
impl CompetitiveAdsorbent {
    pub fn try_new(mut channels: Vec<LangmuirAdsorbent>) -> Result<Self> {
        if channels.is_empty() || channels.len() > MAX_COMPETITIVE_SPECIES {
            return Err(invalid("a competitive site population needs 1..32 explicit adsorbate laws"));
        }
        channels.sort_by(|a,b| a.species.cmp(&b.species));
        let first = &channels[0];
        let mut range = first.temperature_range_k;
        for (i, law) in channels.iter().enumerate() {
            if (i > 0 && law.species == channels[i-1].species)
                || law.capacity_mol_kg_dry != first.capacity_mol_kg_dry
                || law.calorics.dry_skeleton != first.calorics.dry_skeleton
                || law.source.reference_temperature_k != first.source.reference_temperature_k {
                return Err(invalid("competitive species must be unique and share the exact capacity, dry chart and reference temperature"));
            }
            range[0] = range[0].max(law.temperature_range_k[0]);
            range[1] = range[1].min(law.temperature_range_k[1]);
        }
        if range[0] >= range[1] { return Err(invalid("competitive caloric source ranges do not overlap")); }
        Ok(Self { channels, temperature_range_k: range })
    }
    pub fn channels(&self) -> &[LangmuirAdsorbent] { &self.channels }
    pub fn capacity_mol_kg_dry(&self) -> f64 { self.channels[0].capacity_mol_kg_dry }
    pub const fn temperature_range_k(&self) -> [f64; 2] { self.temperature_range_k }
    /// Original finite dry-skeleton chart, shared by all retained site channels.
    /// Bind source-qualified heat paths to this exact borrowed caloric object.
    pub fn dry_skeleton_calorics(&self) -> &ResolvedHeatCapacityCurve { self.dry() }
    fn dry(&self) -> &ResolvedHeatCapacityCurve { &self.channels[0].calorics.dry_skeleton }
}

/// Dry mass, per-species bound inventories and one finite solid temperature.
/// Immutable physical source records are shared across all candidate checkpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveSolid {
    model: Arc<CompetitiveAdsorbent>,
    dry_mass_kg: f64,
    moles: Vec<f64>,
    temperature_k: f64,
}
impl CompetitiveSolid {
    /// Exactly one loading (including explicit zeros) per source species.
    pub fn try_new(model: Arc<CompetitiveAdsorbent>, dry_mass_kg: f64,
        loadings: &[(SpeciesId, f64)], temperature_k: f64) -> Result<Self> {
        if loadings.len() != model.channels.len() { return Err(invalid("every competitive species needs one explicit loading")); }
        let mut ordered = BTreeMap::new();
        for (id, value) in loadings {
            if ordered.insert(id, *value).is_some() { return Err(invalid("duplicate competitive loading species")); }
        }
        let moles = model.channels.iter().map(|law| ordered.get(&law.species).copied()
            .ok_or(invalid("competitive loading names a different species roster"))).collect::<Result<Vec<_>>>()?;
        Self::ordered(model, dry_mass_kg, moles, temperature_k)
    }
    fn ordered(model: Arc<CompetitiveAdsorbent>, dry_mass_kg: f64,
        moles: Vec<f64>, temperature_k: f64) -> Result<Self> {
        let capacity = dry_mass_kg * model.capacity_mol_kg_dry();
        if !positive(dry_mass_kg) || !positive(capacity) || !positive(temperature_k)
            || moles.len() != model.channels.len() || moles.iter().any(|&n| !nonnegative(n))
            || sum(moles.iter().copied())? > capacity
            || temperature_k < model.temperature_range_k[0] || temperature_k > model.temperature_range_k[1] {
            return Err(invalid("competitive state needs supported temperature and nonnegative loadings within ONE finite capacity"));
        }
        model.dry().state_at_temperature(temperature_k)?;
        for (law,&n) in model.channels.iter().zip(&moles) {
            if n > 0.0 && (!positive(n * law.molar_mass()) || n / capacity == 0.0) {
                return Err(invalid("positive competitive trace inventory is unrepresentable"));
            }
            law.calorics.bound_adsorbate.state_at_temperature(temperature_k)?;
        }
        Ok(Self { model, dry_mass_kg, moles, temperature_k })
    }
    pub fn model(&self) -> &CompetitiveAdsorbent { &self.model }
    pub const fn dry_mass_kg(&self) -> f64 { self.dry_mass_kg }
    pub const fn temperature_k(&self) -> f64 { self.temperature_k }
    pub fn capacity_moles(&self) -> f64 { self.dry_mass_kg * self.model.capacity_mol_kg_dry() }
    pub fn adsorbed_moles(&self, id: &SpeciesId) -> Option<f64> {
        self.model.channels.binary_search_by(|law| law.species.cmp(id)).ok().map(|i| self.moles[i])
    }
    pub fn loadings(&self) -> impl Iterator<Item = (&SpeciesId, f64)> {
        self.model.channels.iter().zip(&self.moles).map(|(law,&n)| (&law.species,n))
    }
    pub fn vacancy_moles(&self) -> Result<f64> { finite(self.capacity_moles()-sum(self.moles.iter().copied())?) }
}

/// Physical approximation controls are separate from absolute balance budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompetitiveSolve {
    pub gas: GasSpeciesBoundarySolve,
    pub maximum_gas_gross_fraction: f64,
    pub maximum_site_turnover: f64,
    pub maximum_gas_surface_temperature_difference_k: f64,
    pub maximum_temperature_change_k: f64,
    pub entropy_tolerance_j_k: f64,
    pub maximum_iterations: usize,
}
impl CompetitiveSolve {
    fn validate(self, duration: f64) -> Result<()> {
        AdsorptionSolve { gas: self.gas, maximum_gas_gross_fraction: self.maximum_gas_gross_fraction,
            maximum_site_turnover: self.maximum_site_turnover,
            maximum_gas_surface_temperature_difference_k: self.maximum_gas_surface_temperature_difference_k,
            maximum_temperature_change_k: self.maximum_temperature_change_k,
            equilibrium_coverage_tolerance: 0.0, entropy_tolerance_j_k: self.entropy_tolerance_j_k,
            maximum_iterations: self.maximum_iterations }.validate(duration)
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveSpeciesStep {
    pub species: SpeciesId,
    /// Initial and final shared-site equilibrium coverages, not scalar isotherms.
    pub equilibrium_coverage: [f64; 2],
    pub source: [LangmuirProperties; 2],
    pub energy_into_gas_j: f64,
    pub mass_residual_kg: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveStep {
    pub duration_s: f64,
    /// Arrays have the same canonical SpeciesId order as species below.
    pub sites: CompetitiveSiteStep,
    pub species: Vec<CompetitiveSpeciesStep>,
    pub gas_temperature_k: [f64; 2],
    pub solid_temperature_k: [f64; 2],
    pub dry_sensible_energy_change_j: f64,
    pub bound_sensible_energy_change_j: f64,
    pub solid_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_k: f64,
    pub thermal_iterations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveCandidate<'catalog> {
    solid: CompetitiveSolid,
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    step: CompetitiveStep,
}
impl<'catalog> CompetitiveCandidate<'catalog> {
    pub const fn solid(&self) -> &CompetitiveSolid { &self.solid }
    pub fn gas(&self) -> &ClosedGasParcel { self.gas.gas() }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> { &self.gas }
    pub const fn step(&self) -> &CompetitiveStep { &self.step }
    pub fn into_parts(self) -> (CompetitiveSolid, ClosedGasParcel, &'catalog ThermochemicalGasModel, CompetitiveStep) {
        let (gas,catalog,_) = self.gas.into_parts(); (self.solid,gas,catalog,self.step)
    }
}

fn equilibrium_coverages(properties: &[LangmuirProperties], model: &CompetitiveAdsorbent,
    gas: &ClosedGasParcel) -> Result<Vec<f64>> {
    let mut activity = Vec::with_capacity(properties.len());
    for (p,law) in properties.iter().zip(&model.channels) {
        let pressure = gas.species_partial_pressure_pa(&law.species)?;
        let a = finite(p.affinity_pa_inverse * pressure)?;
        if pressure > 0.0 && a == 0.0 { return Err(invalid("competitive gas activity underflow")); }
        activity.push(a);
    }
    let denominator = finite(1.0 + sum(activity.iter().copied())?)?;
    let vacancy = 1.0 / denominator;
    if !positive(vacancy) { return Err(invalid("finite competitive equilibrium lost vacancy")); }
    for a in &mut activity {
        let theta = *a / denominator;
        if *a > 0.0 && !positive(theta) { return Err(invalid("competitive equilibrium lost a trace species")); }
        *a = theta;
    }
    Ok(activity)
}

fn sensible_change(solid: &CompetitiveSolid, moles: &[f64], temperature: f64) -> Result<[f64;2]> {
    let dry = finite(solid.dry_mass_kg * solid.model.dry().specific_enthalpy_change(solid.temperature_k,temperature)?)?;
    let mut parts = [0.0; MAX_COMPETITIVE_SPECIES];
    for (i,(law,&n)) in solid.model.channels.iter().zip(moles).enumerate() {
        parts[i] = finite((n*law.molar_mass()) * law.calorics.bound_adsorbate
            .specific_enthalpy_change(solid.temperature_k,temperature)?)?;
    }
    Ok([dry,sum(parts[..moles.len()].iter().copied())?])
}
fn solve_temperature(solid: &CompetitiveSolid, moles: &[f64], heat: f64,
    tolerance: f64, iterations: usize, cancelled: &mut impl FnMut()->bool) -> Result<(f64,[f64;2],usize)> {
    if heat == 0.0 { return Ok((solid.temperature_k,[0.0;2],0)); }
    let [mut lo,mut hi] = solid.model.temperature_range_k;
    let a=sensible_change(solid,moles,lo)?; let b=sensible_change(solid,moles,hi)?;
    if heat < sum(a)? || heat > sum(b)? { return Err(AdsorptionError::ThermalRange); }
    let mut t=solid.temperature_k; let mut residual=f64::INFINITY;
    for iteration in 1..=iterations {
        poll(cancelled)?;
        let change=sensible_change(solid,moles,t)?;
        residual=finite(sum(change)?-heat)?;
        if residual.abs() <= tolerance && t != solid.temperature_k { return Ok((t,change,iteration)); }
        if residual>0.0 { hi=t; } else { lo=t; }
        let mut capacity=solid.dry_mass_kg*solid.model.dry().secant_heat_capacity(t,t)?;
        for (law,&n) in solid.model.channels.iter().zip(moles) {
            capacity=finite(capacity+n*law.molar_mass()*law.calorics.bound_adsorbate.secant_heat_capacity(t,t)?)?;
        }
        if !positive(capacity) { return Err(invalid("competitive solid lost positive caloric capacity")); }
        let newton=t-residual/capacity;
        let next=if newton.is_finite() && newton>=lo && newton<=hi && newton!=t { newton } else { 0.5*lo+0.5*hi };
        if next==t { break; } t=next;
    }
    Err(AdsorptionError::NotConverged { residual_j:residual })
}

fn site_mixing_change(old: &CompetitiveSolid, new: &CompetitiveSolid) -> Result<f64> {
    let capacity=old.capacity_moles();
    let change=|a:f64,b:f64| {
        if a==b { 0.0 }
        else if a>0.0 && b>0.0 { (b-a)*(log_ratio(a,capacity)+1.0)+mixing_remainder(a,b) }
        else if a==0.0 { b*log_ratio(b,capacity) }
        else { -a*log_ratio(a,capacity) }
    };
    let components=old.moles.iter().zip(&new.moles).map(|(&a,&b)| change(a,b));
    finite(-R*sum(components.chain(std::iter::once(change(old.vacancy_moles()?,new.vacancy_moles()?))))?)
}
fn combined_entropy_change(old: &CompetitiveSolid, new: &CompetitiveSolid,
    gas: &ClosedGasParcel, next_gas: &ClosedGasParcel, catalog: &ThermochemicalGasModel,
    cancelled: &mut impl FnMut()->bool) -> Result<f64> {
    let reference=old.model.channels[0].source.reference_temperature_k;
    let mut terms=Vec::with_capacity(2+old.moles.len()+catalog.thermodynamics().components().len());
    terms.push(old.dry_mass_kg*old.model.dry().specific_entropy_change(old.temperature_k,new.temperature_k)?);
    terms.push(site_mixing_change(old,new)?);
    for (i,law) in old.model.channels.iter().enumerate() {
        poll(cancelled)?;
        let dn=new.moles[i]-old.moles[i];
        terms.push(new.moles[i]*law.bound_ds(old.temperature_k,new.temperature_k)?
            +dn*law.bound_ds(reference,old.temperature_k)?
            +dn*(-law.reference_heat_j_mol/reference+R*det::ln(law.reference_affinity_pa_inverse)));
    }
    let (t0,t1)=(gas.state().state().temperature,next_gas.state().state().temperature);
    for source in catalog.thermodynamics().components() {
        poll(cancelled)?;
        let n0=gas.species_mass_kg(source.species())?/source.molar_mass().value();
        let n1=next_gas.species_mass_kg(source.species())?/source.molar_mass().value();
        if n0==0.0 && n1==0.0 { continue; }
        let delta=old.model.channels.binary_search_by(|law| law.species.cmp(source.species()))
            .ok().map_or(0.0,|i|old.moles[i]-new.moles[i]);
        let ds=|a,b| source.molar_entropy_change(Temperature::new(a),Temperature::new(b))
            .map(|s|s.value()).map_err(gas_failure);
        let log_p=|n:f64,t:f64|det::ln(n)+det::ln(R)+det::ln(t)-det::ln(gas.volume_m3());
        let pressure=if n0==0.0 { n1*log_p(n1,t1) } else if n1==0.0 { -n0*log_p(n0,t0) }
            else { delta*log_p(n0,t0)+n1*(log_ratio(n1,n0)+log_ratio(t1,t0)) };
        terms.push(n1*ds(t0,t1)?+delta*ds(reference,t0)?-R*pressure);
    }
    sum(terms)
}

// Shared by pair and network coupling. Sampling changes no gas inventory and
// returns the original finite solid proposal, not an independently evolved bath.
struct PreparedAdsorption {
    next: CompetitiveSolid,
    sites: CompetitiveSiteStep,
    before: Vec<LangmuirProperties>,
    transfers: Vec<SpeciesBoundaryTransfer>,
    boundary: Vec<crate::thermochemical_gas::energy::parcel::transfer::species::network::boundary::BoundarySpeciesTransfer>,
    sensible: [f64; 2],
    heat: f64,
    binding: Vec<f64>,
    iterations: usize,
}
fn prepare_adsorption(solid: &CompetitiveSolid, gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel, duration_s: f64, controls: CompetitiveSolve,
    conductive_heat_j: f64,
    mut cancelled: impl FnMut()->bool) -> Result<PreparedAdsorption> {
    use crate::thermochemical_gas::energy::parcel::transfer::species::network::boundary::BoundarySpeciesTransfer;
    controls.validate(duration_s)?; poll(&mut cancelled)?;
    let (ts,tg)=(solid.temperature_k,gas.state().state().temperature);
    if (ts-tg).abs()>controls.maximum_gas_surface_temperature_difference_k {
        return Err(invalid("competitive rates exceed declared near-isothermal support"));
    }
    let mut before=Vec::with_capacity(solid.moles.len());
    let mut adsorption=Vec::with_capacity(solid.moles.len()); let mut desorption=Vec::with_capacity(solid.moles.len());
    for law in &solid.model.channels {
        poll(&mut cancelled)?; law.admit_catalog(catalog)?;
        let p=law.properties(ts,gas.state().state().pressure)?;
        let partial=gas.species_partial_pressure_pa(&law.species)?;
        let activity=finite(p.affinity_pa_inverse*partial)?;
        let a=finite(p.desorption_rate_s_inverse*activity)?;
        if partial>0.0 && (activity==0.0 || (p.desorption_rate_s_inverse>0.0 && a==0.0)) {
            return Err(invalid("competitive adsorption forcing underflow"));
        }
        adsorption.push(a); desorption.push(p.desorption_rate_s_inverse); before.push(p);
    }
    let sites=backward_euler_sites(&solid.moles,solid.capacity_moles(),&adsorption,&desorption,duration_s,&mut cancelled)?;
    let turnover=sum(sites.gross_adsorption_mol.iter().chain(&sites.gross_desorption_mol).copied())?/solid.capacity_moles();
    let mut maximum_gross: f64=0.0;
    for (i,law) in solid.model.channels.iter().enumerate() {
        let amount=finite(sites.gross_adsorption_mol[i]*law.molar_mass())?;
        if sites.gross_adsorption_mol[i]>0.0 && amount==0.0 { return Err(invalid("competitive gross gas withdrawal underflow")); }
        let fraction=if amount==0.0 { 0.0 } else { amount/gas.species_mass_kg(&law.species)? };
        if !nonnegative(fraction) { return Err(AdsorptionError::InventoryLimit { gross_gas_fraction:fraction,site_turnover:turnover }); }
        maximum_gross=maximum_gross.max(fraction);
    }
    if maximum_gross>controls.maximum_gas_gross_fraction || !nonnegative(turnover) || turnover>controls.maximum_site_turnover {
        return Err(AdsorptionError::InventoryLimit { gross_gas_fraction:maximum_gross,site_turnover:turnover });
    }
    let mut transfers=Vec::with_capacity(solid.moles.len());
    let mut boundary=Vec::with_capacity(solid.moles.len());
    let mut heats=Vec::with_capacity(solid.moles.len()); let mut binding=Vec::with_capacity(solid.moles.len());
    for (i,law) in solid.model.channels.iter().enumerate() {
        poll(&mut cancelled)?;
        let dn=sites.final_moles[i]-solid.moles[i]; let dm=finite(dn*law.molar_mass())?;
        if dn!=0.0 && dm==0.0 { return Err(invalid("competitive mass increment underflow")); }
        let site_residual=finite(((sites.gross_adsorption_mol[i]-sites.gross_desorption_mol[i])-dn)*law.molar_mass())?;
        if site_residual.abs()>controls.gas.mass_tolerance_kg { return Err(invalid("competitive gross streams failed species balance")); }
        let incoming=catalog.species_specific_enthalpy_j_kg(&law.species,tg)?*law.molar_mass();
        let dh=law.gas_dh(ts,tg)?;
        let energy=finite(-dn*incoming-sites.gross_desorption_mol[i]*dh)?;
        heats.push(finite(dn*(before[i].adsorption_heat_j_mol+dh)+sites.gross_desorption_mol[i]*dh)?);
        let hb=law.calorics.bound_adsorbate.state_at_temperature(ts)?.specific_enthalpy_j_kg()
            +law.calorics.bound_enthalpy_offset_j_kg;
        binding.push(finite(dm*hb)?);
        let mass_in=finite(sites.gross_desorption_mol[i]*law.molar_mass())?;
        let mass_out=finite(sites.gross_adsorption_mol[i]*law.molar_mass())?;
        if sites.gross_desorption_mol[i]>0.0 && mass_in==0.0 { return Err(invalid("competitive gross desorption underflow")); }
        boundary.push(BoundarySpeciesTransfer { species:law.species.clone(), incoming_mass_kg:mass_in,
            outgoing_mass_kg:mass_out, net_mass_into_gas_kg:-dm,
            relative_energy_into_gas_j:finite(-sites.gross_desorption_mol[i]*dh)? });
        transfers.push(SpeciesBoundaryTransfer { species:law.species.clone(),mass_into_gas_kg:-dm,energy_into_gas_j:energy });
    }
    // The original finite dry+bound caloric inversion receives all physical
    // heat once. Internal conduction is debited by the shared network owner;
    // it neither changes the kinetic snapshot nor supplies a second h chart.
    let heat=if conductive_heat_j==0.0 { sum(heats)? }
        else { sum(heats.into_iter().chain(std::iter::once(conductive_heat_j)))? };
    let (temperature,sensible,iterations)=solve_temperature(solid,&sites.final_moles,heat,
        0.5*controls.gas.energy_tolerance_j,controls.maximum_iterations,&mut cancelled)?;
    if (temperature-ts).abs()>controls.maximum_temperature_change_k { return Err(invalid("competitive solid heating exceeds frozen-rate allowance")); }
    let next=CompetitiveSolid::ordered(Arc::clone(&solid.model),solid.dry_mass_kg,sites.final_moles.clone(),temperature)?;
    if sites.vacancy_moles[1]>0.0 && next.vacancy_moles()?<=0.0 { return Err(invalid("competitive endpoint lost a positive vacancy inventory")); }
    poll(&mut cancelled)?;
    Ok(PreparedAdsorption { next,sites,before,transfers,boundary,sensible,heat,binding,iterations })
}

/// Simultaneous competitive uptake/desorption, finite solid heat and one gas
/// reconstruction. Caller accepts the whole immutable candidate or none of it.
pub fn advance_candidate<'catalog>(solid: &CompetitiveSolid, gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel, duration_s: f64, controls: CompetitiveSolve,
    mut cancelled: impl FnMut()->bool) -> Result<CompetitiveCandidate<'catalog>> {
    use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::SpeciesBoundaryRelativeTransfer;
    let PreparedAdsorption { next,sites,before,transfers,boundary,sensible,heat,binding,iterations } =
        prepare_adsorption(solid,gas,catalog,duration_s,controls,0.0,&mut cancelled)?;
    let (ts,tg)=(solid.temperature_k,gas.state().state().temperature);
    let temperature=next.temperature_k;
    // The source already separates each stream's old-gas partial enthalpy from
    // its desorption-temperature correction. Preserve that small physical heat
    // instead of recovering it from the rounded absolute-energy report rows.
    let relative:Vec<_>=boundary.into_iter().map(|row|SpeciesBoundaryRelativeTransfer {
        species:row.species,mass_into_gas_kg:row.net_mass_into_gas_kg,
        relative_energy_into_gas_j:row.relative_energy_into_gas_j,
    }).collect();
    let candidate=gas.exchange_species_boundary_relative(catalog,&relative,GasSpeciesBoundarySolve {
        energy_tolerance_j:0.5*controls.gas.energy_tolerance_j,..controls.gas },&mut cancelled)?;
    let gas_next=candidate.gas(); let tg_next=gas_next.state().state().temperature;
    if (tg_next-tg).abs()>controls.maximum_temperature_change_k
        || (tg_next-temperature).abs()>controls.maximum_gas_surface_temperature_difference_k {
        return Err(invalid("competitive endpoint temperatures exceed frozen-rate support"));
    }
    let mut after=Vec::with_capacity(solid.moles.len());
    for law in &solid.model.channels { poll(&mut cancelled)?; after.push(law.properties(temperature,gas_next.state().state().pressure)?); }
    let targets0=equilibrium_coverages(&before,&solid.model,gas)?;
    let targets1=equilibrium_coverages(&after,&solid.model,gas_next)?;
    let mut rows=Vec::with_capacity(solid.moles.len()); let mut residuals=Vec::with_capacity(solid.moles.len());
    let mut maximum=candidate.step().maximum_species_mass_residual_kg;
    for (i,(a,b)) in before.into_iter().zip(after).enumerate() {
        let law=&solid.model.channels[i];
        let residual=finite(gas_next.species_mass_kg(&law.species)?-gas.species_mass_kg(&law.species)?
            +(next.moles[i]-solid.moles[i])*law.molar_mass())?;
        maximum=maximum.max(residual.abs()); residuals.push(residual);
        rows.push(CompetitiveSpeciesStep { species:law.species.clone(),equilibrium_coverage:[targets0[i],targets1[i]],
            source:[a,b],energy_into_gas_j:transfers[i].energy_into_gas_j,mass_residual_kg:residual });
    }
    // Include all carrier mass rounding via the original gas owner's total account.
    let mass_residual=finite(candidate.step().mass_residual_kg)?;
    if maximum>controls.gas.mass_tolerance_kg || mass_residual.abs()>controls.gas.mass_tolerance_kg
        || sum(residuals)?.abs()>controls.gas.mass_tolerance_kg { return Err(invalid("competitive gas/site inventory balance failed")); }
    let residual=finite((sum(sensible)?-heat)+candidate.step().energy_residual_j)?;
    if residual.abs()>controls.gas.energy_tolerance_j { return Err(invalid("competitive joint energy balance failed")); }
    let entropy=combined_entropy_change(solid,&next,gas,gas_next,catalog,&mut cancelled)?;
    if entropy < -controls.entropy_tolerance_j_k { return Err(AdsorptionError::EntropyDecreased { change_j_k:entropy }); }
    let step=CompetitiveStep { duration_s,sites,species:rows,gas_temperature_k:[tg,tg_next],solid_temperature_k:[ts,temperature],
        dry_sensible_energy_change_j:sensible[0],bound_sensible_energy_change_j:sensible[1],
        solid_energy_change_j:finite(sum(sensible)?+sum(binding)?)?,gas_energy_change_j:candidate.step().stored_energy_change_j,
        mass_residual_kg:mass_residual,maximum_species_mass_residual_kg:maximum,energy_residual_j:residual,
        entropy_change_j_k:entropy,thermal_iterations:iterations };
    poll(&mut cancelled)?;
    Ok(CompetitiveCandidate { solid:next,gas:candidate,step })
}

#[cfg(test)]
mod tests;
