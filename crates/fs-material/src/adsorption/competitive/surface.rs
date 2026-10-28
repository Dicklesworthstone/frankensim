//! Isothermal exclusion hopping between two finite populations of common sites.
//!
//! A source-defined vacancy-hopping diffusivity D_i [m2/s], site density Gamma
//! [mol/m2], connecting width w and distance L give G_i=D_i*Gamma*w/L [mol/s].
//! Flux is G_i*(theta_Li*theta_Rv-theta_Ri*theta_Lv). Vacancies are shared among
//! all species, not independent binary vacancies. A finite population's exposed
//! site area is C/Gamma; width is NOT a bulk cross-sectional area.
//!
//! Conservation eliminates right inventories. With N_i=n_Li+n_Ri and total
//! vacancies V=v_L+v_R, dn_Li/dt=(G_i/(C_L*C_R))*(N_i*v_L-V*n_Li).
//! The original competitive O(n) backward-Euler solve therefore applies exactly
//! to these discretized equations; no nonlinear iteration or dense inverse is
//! required. Directional gross streams are endpoint-rate quadratures, not exact
//! continuous-time hop counts. They may include recrossing within the interval.
//!
//! Both solids must share the SAME original material model and temperature.
//! Bound-species internal energy travels with the actual inventory. Ideal-site
//! mixing contributes entropy but no extra heat. Neither gas nor a thermostat
//! participates. Unequal temperatures, unlike surfaces, changing geometry,
//! coverage-dependent mobility, nonideal lateral interactions and multilayer
//! adsorption require other constitutive owners and are refused here.

use super::*;

pub const SURFACE_HOPPING_DIFFUSIVITY_PROPERTY: &str = "adsorbed_vacancy_hopping_diffusivity";
pub const SURFACE_HOPPING_DIFFUSIVITY_DIMS: Dims = Dims([2,0,-1,0,0,0]);

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDiffusionGeometry {
    pub connecting_width_m: f64,
    pub distance_m: f64,
    pub site_density_mol_m2: f64,
}
/// Select one mobility claim from each original adsorption channel's material.
/// Missing mobility is not replaced by gas diffusivity or a presumed zero.
#[derive(Clone, Debug)]
pub struct SurfaceDiffusionSelection {
    pub species: SpeciesId,
    pub selection: MaterialPropertySelection,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceDiffusionCoefficient {
    pub species: SpeciesId,
    pub diffusivity_m2_s: f64,
    pub conductance_mol_s: f64,
    pub source: ResolvedMaterialStatePoint,
}
#[derive(Clone, Debug, PartialEq)]
struct PreparedSurface {
    model: Arc<CompetitiveAdsorbent>,
    temperature_k: f64,
    geometry: SurfaceDiffusionGeometry,
    coefficients: Vec<SurfaceDiffusionCoefficient>,
}
/// Sources and fixed geometry are prepared once and cheaply shared by reports.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceDiffusionLink { prepared: Arc<PreparedSurface> }
impl SurfaceDiffusionLink {
    pub fn try_new(model: Arc<CompetitiveAdsorbent>, temperature_k: f64,
        geometry: SurfaceDiffusionGeometry, selections: &[SurfaceDiffusionSelection],
        mut cancelled: impl FnMut()->bool) -> Result<Self> {
        poll(&mut cancelled)?;
        if !nonnegative(geometry.connecting_width_m) || !positive(geometry.distance_m)
            || !positive(geometry.site_density_mol_m2) || !positive(temperature_k)
            || temperature_k<model.temperature_range_k[0] || temperature_k>model.temperature_range_k[1]
            || selections.len()!=model.channels.len() {
            return Err(invalid("surface diffusion needs supported temperature, explicit geometry and one mobility selection per species"));
        }
        let mut bindings=BTreeMap::new();
        for selection in selections {
            if bindings.insert(&selection.species,&selection.selection).is_some() {
                return Err(invalid("duplicate surface diffusivity species"));
            }
        }
        let requirement=ScalarPropertyRequirement::try_new(SURFACE_HOPPING_DIFFUSIVITY_PROPERTY,
            SURFACE_HOPPING_DIFFUSIVITY_DIMS,ScalarAdmissibility::NonNegative)?;
        let mut coefficients=Vec::with_capacity(model.channels.len());
        for law in &model.channels {
            poll(&mut cancelled)?;
            let selection=bindings.get(&law.species).ok_or(invalid("surface mobility species roster differs from adsorption"))?;
            let pressure=law.source.point.axes().get("pressure").copied().ok_or(invalid("original adsorption context lost pressure"))?;
            law.properties(temperature_k,pressure)?;
            let point=source_point(&law.source.point,temperature_k,pressure)?;
            let source=resolve_material_state_point(&law.source.material,&point,
                std::slice::from_ref(&requirement),(**selection).clone())?;
            let property=source.property(SURFACE_HOPPING_DIFFUSIVITY_PROPERTY).ok_or(invalid("selected surface mobility is missing"))?;
            if !property.is_constant_at_fixed_temperature()
                || !property.has_only_owned_and_fixed_external_conditions("temperature",temp_quantity())
                || property.source_axes().iter().any(|(axis,_)|axis!="temperature") {
                return Err(invalid("surface hopping mobility must be constant at the actual fixed temperature with no unowned coverage or dynamic dependence"));
            }
            let diffusivity=property.value_si();
            let conductance=finite((diffusivity*geometry.site_density_mol_m2)
                *(geometry.connecting_width_m/geometry.distance_m))?;
            if !nonnegative(conductance) || (diffusivity>0.0 && geometry.connecting_width_m>0.0 && conductance==0.0) {
                return Err(invalid("nonzero surface hopping conductance is unrepresentable"));
            }
            coefficients.push(SurfaceDiffusionCoefficient { species:law.species.clone(),diffusivity_m2_s:diffusivity,
                conductance_mol_s:conductance,source });
        }
        poll(&mut cancelled)?;
        Ok(Self { prepared:Arc::new(PreparedSurface { model,temperature_k,geometry,coefficients }) })
    }
    pub fn coefficients(&self)->&[SurfaceDiffusionCoefficient] { &self.prepared.coefficients }
    pub fn geometry(&self)->SurfaceDiffusionGeometry { self.prepared.geometry }
    pub fn temperature_k(&self)->f64 { self.prepared.temperature_k }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceDiffusionSolve {
    /// Net removal divided by OLD donor species inventory, in (0,1].
    pub maximum_net_species_fraction: f64,
    /// Total forward+reverse hopping / smaller site's total capacity.
    pub maximum_site_turnover: f64,
    pub mass_tolerance_kg: f64,
    pub energy_tolerance_j: f64,
    pub entropy_tolerance_j_k: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceSpeciesTransfer {
    pub species: SpeciesId,
    pub moles_left_to_right: f64,
    pub gross_left_to_right_mol: f64,
    pub gross_right_to_left_mol: f64,
    pub bound_energy_left_to_right_j: f64,
    pub species_mass_residual_kg: f64,
    pub kinetic_mass_residual_kg: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceDiffusionStep {
    pub duration_s: f64,
    pub link: SurfaceDiffusionLink,
    pub species: Vec<SurfaceSpeciesTransfer>,
    pub site_turnover: f64,
    pub maximum_net_species_fraction: f64,
    pub mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_k: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceDiffusionCandidate {
    left: CompetitiveSolid,
    right: CompetitiveSolid,
    step: SurfaceDiffusionStep,
}
impl SurfaceDiffusionCandidate {
    pub const fn left(&self)->&CompetitiveSolid { &self.left }
    pub const fn right(&self)->&CompetitiveSolid { &self.right }
    pub const fn step(&self)->&SurfaceDiffusionStep { &self.step }
    pub fn into_parts(self)->(CompetitiveSolid,CompetitiveSolid,SurfaceDiffusionStep) { (self.left,self.right,self.step) }
}

impl SurfaceDiffusionLink {
    /// Advance both bound phases together without spending any gas inventory.
    /// Finite-step budgets, original equations and total mixing entropy are
    /// checked before either immutable endpoint is returned.
    pub fn advance_candidate(&self,left:&CompetitiveSolid,right:&CompetitiveSolid,
        duration_s:f64,controls:SurfaceDiffusionSolve,mut cancelled:impl FnMut()->bool)
        ->Result<SurfaceDiffusionCandidate> {
        poll(&mut cancelled)?;
        if !nonnegative(duration_s) || !positive(controls.maximum_net_species_fraction)
            || controls.maximum_net_species_fraction>1.0 || !positive(controls.maximum_site_turnover)
            || !positive(controls.mass_tolerance_kg) || !positive(controls.energy_tolerance_j)
            || !positive(controls.entropy_tolerance_j_k) || std::ptr::eq(left,right) {
            return Err(invalid("surface diffusion requires distinct finite owners and positive physical/numerical budgets"));
        }
        for solid in [left,right] {
            if solid.temperature_k!=self.prepared.temperature_k
                || (!Arc::ptr_eq(&solid.model,&self.prepared.model) && solid.model!=self.prepared.model) {
                return Err(invalid("surface link requires its original common material and exactly equal fixed temperatures"));
            }
        }
        let (cl,cr)=(left.capacity_moles(),right.capacity_moles());
        let vacancies=sum([left.vacancy_moles()?,right.vacancy_moles()?])?;
        let n=left.moles.len(); let mut a=[0.0;MAX_COMPETITIVE_SPECIES]; let mut d=a;
        for i in 0..n {
            poll(&mut cancelled)?;
            let total=finite(left.moles[i]+right.moles[i])?;
            let rate=self.prepared.coefficients[i].conductance_mol_s/cl;
            a[i]=finite(rate*(total/cr))?; d[i]=finite(rate*(vacancies/cr))?;
            if self.prepared.coefficients[i].conductance_mol_s>0.0
                && ((total>0.0 && a[i]==0.0) || (vacancies>0.0 && d[i]==0.0)) {
                return Err(invalid("surface inventory-scaled kinetic coefficient underflow"));
            }
        }
        let site=backward_euler_sites(&left.moles,cl,&a[..n],&d[..n],duration_s,&mut cancelled)?;
        let next_left=CompetitiveSolid::ordered(Arc::clone(&left.model),left.dry_mass_kg,site.final_moles,left.temperature_k)?;
        let mut right_moles=Vec::with_capacity(n);
        for i in 0..n { right_moles.push(finite(right.moles[i]+(left.moles[i]-next_left.moles[i]))?); }
        let next_right=CompetitiveSolid::ordered(Arc::clone(&right.model),right.dry_mass_kg,right_moles,right.temperature_k)?;
        let (vl,vr)=(next_left.vacancy_moles()?,next_right.vacancy_moles()?);
        if site.vacancy_moles[1]>0.0 && vl<=0.0 { return Err(invalid("surface endpoint lost positive vacant sites")); }
        let mut rows=Vec::with_capacity(n); let mut gross=Vec::with_capacity(n);
        let mut masses=Vec::with_capacity(n); let mut energies=Vec::with_capacity(n);
        let mut largest: f64=0.0;
        for i in 0..n {
            poll(&mut cancelled)?;
            let law=&left.model.channels[i]; let g=self.prepared.coefficients[i].conductance_mol_s;
            let forward=finite((duration_s*g)*((next_left.moles[i]/cl)*(vr/cr)))?;
            let reverse=finite((duration_s*g)*((next_right.moles[i]/cr)*(vl/cl)))?;
            if duration_s>0.0 && g>0.0 && ((next_left.moles[i]>0.0 && vr>0.0 && forward==0.0)
                || (next_right.moles[i]>0.0 && vl>0.0 && reverse==0.0)) {
                return Err(invalid("nonzero surface hopping stream underflow"));
            }
            let transfer=left.moles[i]-next_left.moles[i];
            let fraction=if transfer>0.0 { transfer/left.moles[i] } else if transfer<0.0 { -transfer/right.moles[i] } else { 0.0 };
            if !nonnegative(fraction) { return Err(invalid("surface diffusion tried to withdraw an absent donor species")); }
            largest=largest.max(fraction);
            let molecular_residual=finite((next_left.moles[i]-left.moles[i])+(next_right.moles[i]-right.moles[i]))?;
            let mass_residual=finite(molecular_residual*law.molar_mass())?;
            let kinetic_residual=finite((transfer-(forward-reverse))*law.molar_mass())?;
            if mass_residual.abs()>controls.mass_tolerance_kg || kinetic_residual.abs()>controls.mass_tolerance_kg {
                return Err(invalid("original surface hopping equations or per-species inventory failed"));
            }
            let hb=finite((law.calorics.bound_adsorbate.state_at_temperature(left.temperature_k)?.specific_enthalpy_j_kg()
                +law.calorics.bound_enthalpy_offset_j_kg)*law.molar_mass())?;
            masses.push(mass_residual); energies.push(finite(molecular_residual*hb)?);
            gross.push(finite(forward+reverse)?);
            rows.push(SurfaceSpeciesTransfer { species:law.species.clone(),moles_left_to_right:transfer,
                gross_left_to_right_mol:forward,gross_right_to_left_mol:reverse,
                bound_energy_left_to_right_j:finite(transfer*hb)?,species_mass_residual_kg:mass_residual,
                kinetic_mass_residual_kg:kinetic_residual });
        }
        let turnover=sum(gross)?/cl.min(cr);
        if !nonnegative(turnover) || turnover>controls.maximum_site_turnover || largest>controls.maximum_net_species_fraction {
            return Err(invalid("surface hopping exceeds its gross-turnover or net-species withdrawal allowance"));
        }
        let mass_residual=sum(masses)?; let energy_residual=sum(energies)?;
        if mass_residual.abs()>controls.mass_tolerance_kg || energy_residual.abs()>controls.energy_tolerance_j {
            return Err(invalid("combined surface species mass or binding-energy balance failed"));
        }
        // Equal temperatures and exact common species laws cancel every bound
        // reference entropy. The ONE vacancy per solid supplies all mixing.
        let entropy=sum([site_mixing_change(left,&next_left)?,site_mixing_change(right,&next_right)?])?;
        if entropy < -controls.entropy_tolerance_j_k { return Err(AdsorptionError::EntropyDecreased { change_j_k:entropy }); }
        let step=SurfaceDiffusionStep { duration_s,link:self.clone(),species:rows,site_turnover:turnover,
            maximum_net_species_fraction:largest,mass_residual_kg:mass_residual,energy_residual_j:energy_residual,
            entropy_change_j_k:entropy };
        poll(&mut cancelled)?;
        Ok(SurfaceDiffusionCandidate { left:next_left,right:next_right,step })
    }
}

#[cfg(test)]
mod tests;
