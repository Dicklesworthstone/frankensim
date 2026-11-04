//! Ordered-layer dusty-gas resistance with the native Darcy laminate.
//!
//! Each rigid layer contributes its full Dusty-Gas friction matrix integrated
//! across physical thickness. Molecular and Knudsen resistances are added as
//! matrices, not by averaging pore radius, porosity, tortuosity or diffusivity.
//! The existing layered Darcy owner supplies the common viscous molar flux.
//!
//! Klinkenberg slip and Forchheimer inertia are refused because the Dusty-Gas
//! wall term already owns rarefaction and this rung does not derive an inertial
//! multicomponent closure. A fully gas-blocking serial layer closes the path.

use super::*;
use crate::porous::permeability::MAX_POROUS_LAYERS;

#[derive(Clone, Debug, PartialEq)]
pub struct LayeredDustyGasMedium {
    barrier: Arc<DarcyGasBarrier>,
    pores: Arc<[DustyGasPores]>,
}
impl LayeredDustyGasMedium {
    pub fn try_new(barrier: DarcyGasBarrier, pores: Vec<DustyGasPores>)
        -> Result<Self, GasEnergyError>
    {
        if barrier.layers().is_empty() || barrier.layers().len() > MAX_POROUS_LAYERS
            || barrier.layers().len() != pores.len() {
            return Err(invalid("layered dusty gas needs one explicit pore geometry per bounded Darcy layer"));
        }
        if barrier.inertial_geometric_resistance() != 0.0 {
            return Err(invalid("layered dusty gas does not qualify an added Forchheimer inertial law"));
        }
        if barrier.layers().iter().any(|layer| layer.klinkenberg_slip_pa() != 0.0) {
            return Err(invalid("layered dusty gas requires no-slip intrinsic permeability; Klinkenberg rarefaction would double count wall transport"));
        }
        Ok(Self { barrier: Arc::new(barrier), pores: Arc::from(pores) })
    }
    pub fn barrier(&self) -> &DarcyGasBarrier { &self.barrier }
    pub fn pore_layers(&self) -> &[DustyGasPores] { &self.pores }
    pub fn thickness_m(&self) -> f64 {
        self.barrier.layers().iter().map(|layer| layer.thickness_m()).sum()
    }
    fn layer_factor(&self, layer: usize) -> f64 {
        self.barrier.layers()[layer].wet_state().map_or(
            self.pores[layer].effective_diffusivity_factor,
            |state| state.effective_diffusivity_factor(),
        )
    }
    fn layer_radius_m(&self, layer: usize) -> f64 {
        self.barrier.layers()[layer].wet_state().map_or(
            self.pores[layer].mean_pore_radius_m,
            |state| state.mean_gas_pore_radius_m(),
        )
    }
    fn gas_path_open(&self) -> bool {
        self.barrier.layers().iter().enumerate()
            .all(|(i, _)| self.layer_factor(i) > 0.0 && self.layer_radius_m(i) > 0.0)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LayeredDustyGasSpeciesFlux {
    pub species: SpeciesId,
    pub face_mole_fraction: f64,
    pub effective_knudsen_diffusivity_m2_s: f64,
    pub knudsen_resistance_s_m: f64,
    pub diffusive_molar_flux_mol_m2_s: f64,
    pub viscous_molar_flux_mol_m2_s: f64,
    pub molar_flux_mol_m2_s: f64,
    pub mass_rate_left_to_right_kg_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LayeredDustyGasFlux {
    pub medium: LayeredDustyGasMedium,
    pub species: Vec<LayeredDustyGasSpeciesFlux>,
    pub pairs: Vec<BinaryDiffusionSample>,
    pub face_state: Option<ResolvedThermochemicalGasState>,
    pub darcy: Option<DarcyGasFlow>,
    pub maximum_integrated_composition_residual: f64,
    pub entropy_production_w_m2_k: f64,
    pub net_molar_flux_mol_m2_s: f64,
}

impl ClosedGasParcel {
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn layered_dusty_gas_flux_to(
        &self, right: &Self, medium: &LayeredDustyGasMedium, area_m2: f64,
        pairs: &[BinaryDiffusionMaterial<'_>], gas_context: &QueryPoint,
        controls: DustyGasFaceSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayeredDustyGasFlux, GasEnergyError> {
        poll(&mut cancelled)?; controls.validate()?;
        if !area_m2.is_finite() || area_m2 < 0.0
            || pairs.len() > MAX_DIFFUSION_SPECIES*(MAX_DIFFUSION_SPECIES-1)/2 {
            return Err(invalid("nonnegative area and bounded layered dusty-gas pair list required"));
        }
        let mut report = LayeredDustyGasFlux {
            medium: medium.clone(), species: Vec::new(), pairs: Vec::new(),
            face_state: None, darcy: None, maximum_integrated_composition_residual: 0.0,
            entropy_production_w_m2_k: 0.0, net_molar_flux_mol_m2_s: 0.0,
        };
        if area_m2 == 0.0 || !medium.gas_path_open() {
            poll(&mut cancelled)?; return Ok(report);
        }
        let distance = medium.thickness_m();
        if !positive(distance) {
            return Err(invalid("layered dusty-gas total thickness is unrepresentable"));
        }

        let (left, right_state) = (self.state().state(), right.state().state());
        let temperature = 0.5*left.temperature + 0.5*right_state.temperature;
        let pressure = 0.5*left.pressure + 0.5*right_state.pressure;
        let pressure_difference = right_state.pressure-left.pressure;
        if (left.temperature-right_state.temperature).abs() > controls.maximum_temperature_difference_k
            || pressure_difference.abs()/pressure > controls.maximum_relative_pressure_difference {
            return Err(invalid("gas endpoints exceed the declared layered dusty-gas face approximation"));
        }
        let support = medium.barrier.support();
        if [left.temperature,right_state.temperature].iter().any(|&t|
            t < support.temperature_k.0 || t > support.temperature_k.1) {
            return Err(invalid("layered dusty-gas endpoint temperature is outside porous source support"));
        }

        let joined = catalog(self,right)?;
        let components = joined.thermodynamics().components();
        let n = components.len();
        if n == 0 || n > MAX_DIFFUSION_SPECIES || pairs.len() != n*(n-1)/2 {
            return Err(invalid("layered dusty gas needs one selected binary material per active pair and at most 32 species"));
        }
        let ids:BTreeMap<_,_> = components.iter().enumerate()
            .map(|(i,s)|(s.species().clone(),i)).collect();
        let mut ordered=BTreeMap::new();
        for pair in pairs {
            poll(&mut cancelled)?;
            let i=*ids.get(&pair.species[0]).ok_or(invalid("unknown layered dusty-gas pair species"))?;
            let j=*ids.get(&pair.species[1]).ok_or(invalid("unknown layered dusty-gas pair species"))?;
            if i==j || ordered.insert((i.min(j),i.max(j)),pair).is_some() {
                return Err(invalid("duplicate or self layered dusty-gas binary pair"));
            }
        }

        let (left_x,right_x)=(fractions(self),fractions(right));
        let concentration=pressure/(R*temperature);
        if !positive(concentration) { return Err(invalid("layered dusty-gas concentration is unrepresentable")); }
        let mut x=Vec::with_capacity(n);
        let mut delta=Vec::with_capacity(n);
        let mut knudsen_resistance=Vec::with_capacity(n);
        let mut effective_knudsen=Vec::with_capacity(n);
        let mut minimum_effective=f64::INFINITY;
        for species in components {
            poll(&mut cancelled)?;
            let lx=left_x.get(species.species()).copied().unwrap_or(0.0);
            let rx=right_x.get(species.species()).copied().unwrap_or(0.0);
            let a=left.pressure*lx; let b=right_state.pressure*rx;
            let fraction=(0.5*a+0.5*b)/pressure;
            let drive=(b-a)/pressure;
            let molecular_speed=det::sqrt((8.0/core::f64::consts::PI)
                *(R*temperature/species.molar_mass().value()));
            let mut resistance=Sum::default();
            for (layer,porous) in medium.barrier.layers().iter().enumerate() {
                let dk=((2.0/3.0)*medium.layer_radius_m(layer))
                    *medium.layer_factor(layer)*molecular_speed;
                if !positive(dk) { return Err(invalid("layered dusty-gas Knudsen coefficient is unrepresentable")); }
                resistance.add(porous.thickness_m()/dk)?;
            }
            let effective=distance/resistance.value;
            if !positive(fraction) || !drive.is_finite() || !positive(resistance.value)
                || !positive(effective) || (lx!=0.0 && a==0.0) || (rx!=0.0 && b==0.0)
                || (a!=b && drive==0.0) {
                return Err(invalid("layered dusty-gas face species or integrated Knudsen resistance is unrepresentable"));
            }
            x.push(fraction); delta.push(drive); knudsen_resistance.push(resistance.value);
            effective_knudsen.push(effective); minimum_effective=minimum_effective.min(effective);
        }
        fs_thermal::enthalpy::advection::species::close_generated_simplex(&mut x).map_err(source_error)?;
        let composition=Composition::new(CompositionBasis::MoleFraction,x.clone()).map_err(source_error)?;
        let face_model=joined.with_composition(composition)?;
        let point=updated_point(gas_context,temperature,pressure)?;
        let face=face_model.evaluate(&point)?;
        let gas=face.state();
        let darcy=medium.barrier.mass_flow_for_model(gas,&face_model,
            left.pressure,right_state.pressure,area_m2)?;
        let viscous_molar_flux=darcy.mass_flux_left_to_right_kg_m2_s
            /(R/gas.specific_gas_constant);
        if !viscous_molar_flux.is_finite() {
            return Err(invalid("layered dusty-gas viscous flux is unrepresentable"));
        }

        let requirement=ScalarPropertyRequirement::try_new(BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,ScalarAdmissibility::StrictlyPositive).map_err(source_error)?;
        let mut pair_resistance=vec![0.0;n*n];
        for ((i,j),binding) in ordered {
            poll(&mut cancelled)?;
            let point=updated_point(binding.context,temperature,pressure)?;
            let material:ResolvedMaterialStatePoint=resolve_material_state_point(
                binding.material,&point,std::slice::from_ref(&requirement),binding.selection.clone())
                .map_err(source_error)?;
            let d=material.property(BINARY_DIFFUSIVITY_PROPERTY)
                .ok_or(invalid("resolved layered dusty-gas binary diffusivity is missing"))?.value_si();
            for state in [left,right_state] {
                poll(&mut cancelled)?;
                let boundary=updated_point(binding.context,state.temperature,state.pressure)?;
                material.at_query_point(binding.material,&boundary).map_err(source_error)?;
            }
            let mut resistance=Sum::default();
            for (layer,porous) in medium.barrier.layers().iter().enumerate() {
                let effective=medium.layer_factor(layer)*d;
                if !positive(effective) {
                    return Err(invalid("layered effective binary diffusivity is unrepresentable"));
                }
                resistance.add(porous.thickness_m()/effective)?;
            }
            let effective=distance/resistance.value;
            if !positive(resistance.value)||!positive(effective) {
                return Err(invalid("layered binary diffusion resistance is unrepresentable"));
            }
            pair_resistance[i*n+j]=resistance.value; pair_resistance[j*n+i]=resistance.value;
            minimum_effective=minimum_effective.min(effective);
            report.pairs.push(BinaryDiffusionSample {
                species:[components[i].species().clone(),components[j].species().clone()],
                diffusivity_m2_s:d,material,
            });
        }

        let scale=minimum_effective/distance;
        if !positive(scale) { return Err(invalid("layered dusty-gas resistance scale is unrepresentable")); }
        let root_x:Vec<_>=x.iter().map(|&value|det::sqrt(value)).collect();
        let mut matrix=vec![0.0;n*n]; let mut weighted=vec![0.0;n];
        for i in 0..n {
            poll(&mut cancelled)?;
            let wall=scale*knudsen_resistance[i];
            if !positive(wall) { return Err(invalid("layered dusty-gas wall friction scaling underflow")); }
            matrix[i*n+i]=wall;
            weighted[i]=-(concentration*scale)*delta[i]/root_x[i];
            if delta[i]!=0.0 && weighted[i]==0.0 {
                return Err(invalid("layered dusty-gas driving force underflowed during matrix scaling"));
            }
            for j in 0..i {
                let ratio=scale*pair_resistance[i*n+j];
                if !positive(ratio) { return Err(invalid("layered dusty-gas molecular friction scaling underflow")); }
                matrix[i*n+i]+=x[j]*ratio; matrix[j*n+j]+=x[i]*ratio;
                let off=-(root_x[i]*root_x[j])*ratio;
                matrix[i*n+j]=off; matrix[j*n+i]=off;
            }
        }
        if matrix.iter().chain(&weighted).any(|v|!v.is_finite()) {
            return Err(invalid("layered dusty-gas friction system is unrepresentable"));
        }
        poll(&mut cancelled)?;
        let factor=fs_la::factor::cholesky(&matrix,n).map_err(source_error)?;
        factor.solve(&mut weighted);
        poll(&mut cancelled)?;

        let mut reduced=Vec::with_capacity(n); let mut net=Sum::default(); let mut gross=Sum::default();
        for (i,species) in components.iter().enumerate() {
            poll(&mut cancelled)?;
            let diffusive=root_x[i]*weighted[i]; let viscous=x[i]*viscous_molar_flux;
            let total=diffusive+viscous;
            let rate=(total*species.molar_mass().value())*area_m2;
            if !total.is_finite() || !rate.is_finite() || (total!=0.0&&rate==0.0)
                || (weighted[i]!=0.0&&diffusive==0.0) || (viscous_molar_flux!=0.0&&viscous==0.0) {
                return Err(invalid("layered dusty-gas species flux is unrepresentable"));
            }
            reduced.push(total-viscous); net.add(total)?; gross.add(total.abs())?;
            report.species.push(LayeredDustyGasSpeciesFlux {
                species:species.species().clone(),face_mole_fraction:x[i],
                effective_knudsen_diffusivity_m2_s:effective_knudsen[i],
                knudsen_resistance_s_m:knudsen_resistance[i],
                diffusive_molar_flux_mol_m2_s:diffusive,viscous_molar_flux_mol_m2_s:viscous,
                molar_flux_mol_m2_s:total,mass_rate_left_to_right_kg_s:rate,
            });
        }
        let minimum_concentration=left.pressure.min(right_state.pressure)/(R*temperature);
        if gross.value/minimum_concentration > support.maximum_superficial_speed_m_s {
            return Err(invalid("gross layered dusty-gas superficial transport exceeds porous source support"));
        }

        let mut dissipation=Sum::default();
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut residual=Sum::default(); residual.add(delta[i])?;
            residual.add(reduced[i]*knudsen_resistance[i]/concentration)?;
            let yi=reduced[i]/root_x[i];
            dissipation.add(R*yi*yi*knudsen_resistance[i]/concentration)?;
            for j in 0..n {
                if i==j { continue; }
                let resistance=pair_resistance[i*n+j];
                residual.add((x[j]*reduced[i]-x[i]*reduced[j])*resistance/concentration)?;
                if j<i {
                    let yj=reduced[j]/root_x[j];
                    let drag=root_x[j]*yi-root_x[i]*yj;
                    dissipation.add(R*drag*drag*resistance/concentration)?;
                }
            }
            report.maximum_integrated_composition_residual=
                report.maximum_integrated_composition_residual.max(residual.value.abs());
        }
        let viscous_entropy=-viscous_molar_flux*R*(pressure_difference/pressure);
        if !viscous_entropy.is_finite() || viscous_entropy<0.0 {
            return Err(invalid("layered dusty-gas viscous dissipation is unrepresentable"));
        }
        dissipation.add(viscous_entropy)?;
        let integrated_tolerance=controls.gradient_tolerance_per_m*distance;
        if !positive(integrated_tolerance)
            || report.maximum_integrated_composition_residual>integrated_tolerance {
            return Err(invalid("integrated layered dusty-gas equations did not close within the gradient tolerance"));
        }
        if dissipation.value<0.0 || !dissipation.value.is_finite() {
            return Err(invalid("layered dusty-gas entropy production is unrepresentable"));
        }
        report.entropy_production_w_m2_k=dissipation.value;
        report.net_molar_flux_mol_m2_s=net.value;
        report.face_state=Some(face); report.darcy=Some(darcy);
        poll(&mut cancelled)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
