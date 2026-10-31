//! Isothermal ideal-gas dusty-gas transport in a homogeneous porous slab.
//!
//! Fluxes are relative to the STATIONARY SOLID and per bulk cross-sectional
//! area. Molecular collisions, Knudsen wall collisions and viscous Darcy flow
//! share one physical path. No zero-molar or zero-mass-flux constraint applies.
//! With A_ii=1/DK_i+sum_j(x_j/Dij), A_ij=-x_i/Dij, solve A Nd=-grad(c),
//! then N_i=Nd_i+x_i Nv. A*x=x_i/DK_i makes this decomposition exact.
//! The existing no-slip Darcy owner supplies Nv; Klinkenberg slip would count
//! overlapping rarefaction physics twice and is refused here.
//!
//! D_ij_eff=(porosity/tortuosity)*D_ij and
//! DK_i=(porosity/tortuosity)*(2*r/3)*sqrt(8*R*T/(pi*M_i)). Tortuosity is the
//! transport factor in this equation, not an implicitly squared path ratio or
//! the JCA acoustic high-frequency tortuosity. Primary formulation:
//! https://cantera.org/3.1/cxx/dc/de6/classCantera_1_1DustyGasTransport.html
//!
//! The face uses concentration-weighted composition at one declared mean T;
//! the admitted small temperature difference is ignored in the driving force.
//! No thermal transpiration, Soret effect, adsorption, surface diffusion,
//! pore storage, pore-size distribution or Forchheimer
//! inertia is inferred. Central face sampling requires spatial/time refinement;
//! it is not a resolved steady nonlinear composition profile through the slab.

use super::{BinaryDiffusionMaterial, BinaryDiffusionSample, BINARY_DIFFUSIVITY_DIMS,
    BINARY_DIFFUSIVITY_PROPERTY, MAX_DIFFUSION_SPECIES, fractions};
use super::super::{ClosedGasParcel, GasEnergyError, SpeciesId, Sum, catalog,
    invalid, poll, positive, source_error};
use crate::porous::permeability::{DarcyGasBarrier, DarcyGasFlow};
use crate::state_point::{ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use crate::thermochemical_gas::{ResolvedThermochemicalGasState,
    energy::updated_point};
use fs_matdb::QueryPoint;
use fs_math::det;
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::{collections::BTreeMap, sync::Arc};

/// Authored pore geometry. These inputs do not acquire measured authority from
/// a citation; pore-size/tortuosity applicability must fit the physical solid.
#[derive(Clone, Debug, PartialEq)]
pub struct DustyGasPores {
    mean_pore_radius_m: f64,
    gas_filled_porosity: f64,
    transport_tortuosity_factor: f64,
    effective_diffusivity_factor: f64,
    provenance: String,
}
impl DustyGasPores {
    pub fn try_new(mean_pore_radius_m: f64, gas_filled_porosity: f64,
        transport_tortuosity_factor: f64, provenance: impl Into<String>)
        -> Result<Self, GasEnergyError>
    {
        let provenance = provenance.into();
        let factor = gas_filled_porosity / transport_tortuosity_factor;
        if !positive(mean_pore_radius_m) || !positive(gas_filled_porosity)
            || gas_filled_porosity > 1.0 || !transport_tortuosity_factor.is_finite()
            || transport_tortuosity_factor < 1.0 || !positive(factor)
            || provenance.trim().is_empty() || provenance.len() > 4096 {
            return Err(invalid("dusty gas needs positive pore radius, porosity in (0,1], transport tortuosity >=1 and bounded provenance"));
        }
        Ok(Self { mean_pore_radius_m, gas_filled_porosity, transport_tortuosity_factor,
            effective_diffusivity_factor: factor, provenance })
    }
    pub const fn mean_pore_radius_m(&self) -> f64 { self.mean_pore_radius_m }
    pub const fn gas_filled_porosity(&self) -> f64 { self.gas_filled_porosity }
    pub const fn transport_tortuosity_factor(&self) -> f64 { self.transport_tortuosity_factor }
    pub const fn effective_diffusivity_factor(&self) -> f64 { self.effective_diffusivity_factor }
    pub fn provenance(&self) -> &str { &self.provenance }
}

/// One homogeneous reference pore geometry and the original permeability law.
/// The Arc shares immutable source records across reports. Zero intrinsic k
/// disables viscous flow only; explicit pore diffusivity can still conduct gas.
#[derive(Clone, Debug, PartialEq)]
pub struct DustyGasMedium {
    barrier: Arc<DarcyGasBarrier>,
    pores: DustyGasPores,
}
impl DustyGasMedium {
    pub fn try_new(barrier: DarcyGasBarrier, pores: DustyGasPores)
        -> Result<Self, GasEnergyError>
    {
        if barrier.inertial_geometric_resistance() != 0.0 {
            return Err(invalid("dusty-gas molecular/Knudsen transport does not qualify an added Forchheimer inertial law"));
        }
        if barrier.layers().len() != 1 || barrier.layers().iter()
            .any(|layer| layer.klinkenberg_slip_pa() != 0.0) {
            return Err(invalid("dusty gas needs one homogeneous no-slip layer; layered pore compositions and added Klinkenberg slip are unsupported"));
        }
        Ok(Self { barrier: Arc::new(barrier), pores })
    }
    pub fn barrier(&self) -> &DarcyGasBarrier { &self.barrier }
    /// Authored reference geometry. The effective getters below include the
    /// current wet response and are the geometry consumed by transport/storage.
    pub const fn reference_pores(&self) -> &DustyGasPores { &self.pores }
    pub fn thickness_m(&self) -> f64 { self.barrier.layers()[0].thickness_m() }
    pub fn wet_state(&self) -> Option<&crate::porous::wet::WetPoreGasState> {
        self.barrier.layers()[0].wet_state()
    }
    pub fn gas_filled_porosity(&self) -> f64 {
        self.wet_state().map_or(self.pores.gas_filled_porosity, |state| state.gas_filled_porosity())
    }
    pub fn mean_gas_pore_radius_m(&self) -> f64 {
        self.wet_state().map_or(self.pores.mean_pore_radius_m, |state| state.mean_gas_pore_radius_m())
    }
    pub fn effective_diffusivity_factor(&self) -> f64 {
        self.wet_state().map_or(self.pores.effective_diffusivity_factor,
            |state| state.effective_diffusivity_factor())
    }
    /// Fully liquid-filled pores have no finite effective gas tortuosity.
    pub fn transport_tortuosity_factor(&self) -> Option<f64> {
        match self.wet_state() {
            Some(state) => state.transport_tortuosity_factor(),
            None => Some(self.pores.transport_tortuosity_factor),
        }
    }
    /// Snapshot the prescribed saturation and pore geometry. Darcy mobility
    /// uses krg once; molecular and Knudsen diffusion use the independently
    /// calibrated absolute diffusion factor once. Rewetting never compounds
    /// an earlier factor. Liquid motion, displaced-gas work and pore-body volume
    /// changes remain the enclosing coupled owner's responsibility.
    pub fn with_wet_state(&self, state: crate::porous::wet::WetPoreGasState)
        -> Result<Self, GasEnergyError>
    {
        Ok(Self { barrier: Arc::new(self.barrier.with_wet_states(&[state])?),
            pores: self.pores.clone() })
    }
}

/// Explicit face approximation and original-equation numerical admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DustyGasFaceSolve {
    pub maximum_temperature_difference_k: f64,
    /// |p_right-p_left| / p_mean. Zero admits only equal total pressure.
    pub maximum_relative_pressure_difference: f64,
    /// Maximum original dusty-gas row residual divided by face c [1/m].
    pub gradient_tolerance_per_m: f64,
}
impl DustyGasFaceSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        if !self.maximum_temperature_difference_k.is_finite()
            || self.maximum_temperature_difference_k < 0.0
            || !self.maximum_relative_pressure_difference.is_finite()
            || !(0.0..2.0).contains(&self.maximum_relative_pressure_difference)
            || !positive(self.gradient_tolerance_per_m) {
            return Err(invalid("dusty-gas face needs nonnegative temperature allowance, relative pressure allowance in [0,2), and positive residual tolerance"));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DustyGasSpeciesFlux {
    pub species: SpeciesId,
    pub face_mole_fraction: f64,
    pub knudsen_diffusivity_m2_s: f64,
    /// Molecular/Knudsen contribution, before the common viscous velocity.
    pub diffusive_molar_flux_mol_m2_s: f64,
    pub viscous_molar_flux_mol_m2_s: f64,
    pub molar_flux_mol_m2_s: f64,
    pub mass_rate_left_to_right_kg_s: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DustyGasFlux {
    pub medium: DustyGasMedium,
    pub species: Vec<DustyGasSpeciesFlux>,
    /// Original binary diffusivities before the explicit porosity/tortuosity factor.
    pub pairs: Vec<BinaryDiffusionSample>,
    /// None for zero area or fully flooded pores; inactive gas sources are not queried.
    pub face_state: Option<ResolvedThermochemicalGasState>,
    pub darcy: Option<DarcyGasFlow>,
    pub maximum_gradient_residual_per_m: f64,
    /// Local positive friction + Darcy dissipation / T, not extra gas heat or
    /// a certificate for a finite update. The exchange owner checks endpoint S.
    pub local_entropy_production_w_m3_k: f64,
    /// A physical result, never a zero-reference-frame constraint.
    pub net_molar_flux_mol_m2_s: f64,
}

impl ClosedGasParcel {
    /// Sample the stationary-solid dusty-gas law without spending any inventory.
    /// Every active binary pair requires its original selected material source.
    /// `gas_context` explicitly owns other face gas-source coordinates; neither
    /// endpoint's unrelated environment is silently chosen as a face condition.
    /// All species rows, including trace and absent-side species, remain intact.
    #[allow(clippy::too_many_arguments)]
    pub fn dusty_gas_flux_to(
        &self, right: &Self, medium: &DustyGasMedium, area_m2: f64,
        pairs: &[BinaryDiffusionMaterial<'_>], gas_context: &QueryPoint,
        controls: DustyGasFaceSolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<DustyGasFlux, GasEnergyError> {
        poll(&mut cancelled)?; controls.validate()?;
        if !area_m2.is_finite() || area_m2 < 0.0
            || pairs.len() > MAX_DIFFUSION_SPECIES*(MAX_DIFFUSION_SPECIES-1)/2 {
            return Err(invalid("nonnegative area and bounded dusty-gas pair list required"));
        }
        let mut report = DustyGasFlux { medium: medium.clone(), species: Vec::new(),
            pairs: Vec::new(), face_state: None, darcy: None,
            maximum_gradient_residual_per_m: 0.0,
            local_entropy_production_w_m3_k: 0.0, net_molar_flux_mol_m2_s: 0.0 };
        if area_m2 == 0.0 || medium.gas_filled_porosity() == 0.0 {
            poll(&mut cancelled)?; return Ok(report);
        }
        let (left, right_state) = (self.state().state(), right.state().state());
        let temperature = 0.5*left.temperature + 0.5*right_state.temperature;
        let pressure = 0.5*left.pressure + 0.5*right_state.pressure;
        let pressure_difference = right_state.pressure-left.pressure;
        if (left.temperature-right_state.temperature).abs() > controls.maximum_temperature_difference_k
            || pressure_difference.abs()/pressure > controls.maximum_relative_pressure_difference {
            return Err(invalid("gas endpoints exceed the declared dusty-gas face approximation"));
        }
        let support = medium.barrier.support();
        if [left.temperature, right_state.temperature].iter().any(|&t|
            t < support.temperature_k.0 || t > support.temperature_k.1) {
            return Err(invalid("dusty-gas endpoint temperature is outside porous source support"));
        }
        let joined = catalog(self, right)?;
        let components = joined.thermodynamics().components();
        let n = components.len();
        if n == 0 || n > MAX_DIFFUSION_SPECIES || pairs.len() != n*(n-1)/2 {
            return Err(invalid("dusty gas needs one selected binary material per active pair and at most 32 species"));
        }
        let ids: BTreeMap<_, _> = components.iter().enumerate()
            .map(|(i, s)| (s.species().clone(), i)).collect();
        let mut ordered = BTreeMap::new();
        for pair in pairs {
            poll(&mut cancelled)?;
            let i = *ids.get(&pair.species[0]).ok_or(invalid("unknown dusty-gas pair species"))?;
            let j = *ids.get(&pair.species[1]).ok_or(invalid("unknown dusty-gas pair species"))?;
            if i == j || ordered.insert((i.min(j), i.max(j)), pair).is_some() {
                return Err(invalid("duplicate or self dusty-gas binary pair"));
            }
        }
        let (left_x, right_x) = (fractions(self), fractions(right));
        let distance = medium.thickness_m();
        let concentration = pressure/(R*temperature);
        if !positive(concentration) { return Err(invalid("dusty-gas concentration is unrepresentable")); }
        let mut x = Vec::with_capacity(n);
        let mut gradient = Vec::with_capacity(n);
        let mut knudsen = Vec::with_capacity(n);
        let mut scale = f64::INFINITY;
        for species in components {
            poll(&mut cancelled)?;
            let lx = left_x.get(species.species()).copied().unwrap_or(0.0);
            let rx = right_x.get(species.species()).copied().unwrap_or(0.0);
            let a = left.pressure*lx;
            let b = right_state.pressure*rx;
            // Concentration-weighted, not arithmetic-average mole fractions.
            let fraction = (0.5*a+0.5*b)/pressure;
            let row_gradient = ((b-a)/pressure)/distance;
            let molecular_speed = det::sqrt((8.0/core::f64::consts::PI)
                *(R*temperature/species.molar_mass().value()));
            let dk = ((2.0/3.0)*medium.mean_gas_pore_radius_m())
                *medium.effective_diffusivity_factor()*molecular_speed;
            if !positive(fraction) || !row_gradient.is_finite() || !positive(dk)
                || (lx != 0.0 && a == 0.0) || (rx != 0.0 && b == 0.0)
                || (a != b && row_gradient == 0.0) {
                return Err(invalid("dusty-gas face species or Knudsen coefficient is unrepresentable"));
            }
            x.push(fraction); gradient.push(row_gradient); knudsen.push(dk); scale = scale.min(dk);
        }
        let composition = Composition::new(CompositionBasis::MoleFraction, x.clone()).map_err(source_error)?;
        let face_model = joined.with_composition(composition)?;
        let point = updated_point(gas_context, temperature, pressure)?;
        let face = face_model.evaluate(&point)?;
        let gas = face.state();
        let darcy = medium.barrier.mass_flow_for_model(gas, &face_model,
            left.pressure, right_state.pressure, area_m2)?;
        let viscous_molar_flux = darcy.mass_flux_left_to_right_kg_m2_s
            /(R/gas.specific_gas_constant);
        if !viscous_molar_flux.is_finite() { return Err(invalid("dusty-gas viscous flux is unrepresentable")); }
        let requirement = ScalarPropertyRequirement::try_new(BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS, ScalarAdmissibility::StrictlyPositive).map_err(source_error)?;
        let mut diffusion = vec![0.0; n*n];
        for ((i, j), binding) in ordered {
            poll(&mut cancelled)?;
            let point = updated_point(binding.context, temperature, pressure)?;
            let material: ResolvedMaterialStatePoint = resolve_material_state_point(
                binding.material, &point, std::slice::from_ref(&requirement), binding.selection.clone())
                .map_err(source_error)?;
            let d = material.property(BINARY_DIFFUSIVITY_PROPERTY)
                .ok_or(invalid("resolved dusty-gas binary diffusivity is missing"))?.value_si();
            for state in [left, right_state] {
                poll(&mut cancelled)?;
                let boundary = updated_point(binding.context, state.temperature, state.pressure)?;
                material.at_query_point(binding.material, &boundary).map_err(source_error)?;
            }
            let effective = medium.effective_diffusivity_factor()*d;
            if !positive(effective) { return Err(invalid("effective binary diffusivity is unrepresentable")); }
            diffusion[i*n+j] = effective; diffusion[j*n+i] = effective; scale = scale.min(effective);
            report.pairs.push(BinaryDiffusionSample { species: [components[i].species().clone(),
                components[j].species().clone()], diffusivity_m2_s: d, material });
        }
        // Similarity by sqrt(x) gives a symmetric positive definite friction
        // matrix. Cholesky uses all rows, with no nullspace/gauge replacement.
        let root_x: Vec<_> = x.iter().map(|&value| det::sqrt(value)).collect();
        let mut matrix = vec![0.0; n*n];
        let mut weighted = vec![0.0; n];
        for i in 0..n {
            poll(&mut cancelled)?;
            let wall = scale/knudsen[i];
            if !positive(wall) { return Err(invalid("dusty-gas wall friction scaling underflow")); }
            matrix[i*n+i] = wall;
            weighted[i] = -(concentration*scale)*gradient[i]/root_x[i];
            if gradient[i] != 0.0 && weighted[i] == 0.0 {
                return Err(invalid("dusty-gas driving force underflowed during matrix scaling"));
            }
            for j in 0..i {
                let ratio = scale/diffusion[i*n+j];
                if !positive(ratio) { return Err(invalid("dusty-gas molecular friction scaling underflow")); }
                matrix[i*n+i] += x[j]*ratio;
                matrix[j*n+j] += x[i]*ratio;
                let off = -(root_x[i]*root_x[j])*ratio;
                matrix[i*n+j] = off; matrix[j*n+i] = off;
            }
        }
        if matrix.iter().chain(&weighted).any(|value| !value.is_finite()) {
            return Err(invalid("dusty-gas friction system is unrepresentable"));
        }
        poll(&mut cancelled)?;
        let factor = fs_la::factor::cholesky(&matrix, n).map_err(source_error)?;
        factor.solve(&mut weighted);
        poll(&mut cancelled)?;
        let mut reduced = Vec::with_capacity(n);
        let mut net = Sum::default(); let mut gross = Sum::default();
        for (i, species) in components.iter().enumerate() {
            poll(&mut cancelled)?;
            let diffusive = root_x[i]*weighted[i];
            let viscous = x[i]*viscous_molar_flux;
            let total = diffusive+viscous;
            let rate = (total*species.molar_mass().value())*area_m2;
            if !total.is_finite() || !rate.is_finite() || (total != 0.0 && rate == 0.0)
                || (weighted[i] != 0.0 && diffusive == 0.0)
                || (viscous_molar_flux != 0.0 && viscous == 0.0) {
                return Err(invalid("dusty-gas species flux is unrepresentable"));
            }
            // Check the actual assembled total: adding a dominant viscous term
            // must not silently discard a resolved diffusive contribution.
            reduced.push(total-viscous);
            net.add(total)?; gross.add(total.abs())?;
            report.species.push(DustyGasSpeciesFlux { species: species.species().clone(),
                face_mole_fraction: x[i], knudsen_diffusivity_m2_s: knudsen[i],
                diffusive_molar_flux_mol_m2_s: diffusive, viscous_molar_flux_mol_m2_s: viscous,
                molar_flux_mol_m2_s: total, mass_rate_left_to_right_kg_s: rate });
        }
        // A net counterflow can cancel; use gross molar transport for conservative
        // superficial-speed admission at the lower-concentration face.
        let minimum_concentration = left.pressure.min(right_state.pressure)/(R*temperature);
        if gross.value/minimum_concentration > support.maximum_superficial_speed_m_s {
            return Err(invalid("gross dusty-gas superficial transport exceeds porous source support"));
        }
        let mut dissipation = Sum::default();
        for i in 0..n {
            poll(&mut cancelled)?;
            let wall_denominator = concentration*knudsen[i];
            if !positive(wall_denominator) { return Err(invalid("dusty-gas residual scaling is unrepresentable")); }
            let mut residual = Sum::default(); residual.add(gradient[i])?;
            residual.add(reduced[i]/wall_denominator)?;
            let yi = reduced[i]/root_x[i];
            let wall = yi/det::sqrt(wall_denominator);
            dissipation.add(R*wall*wall)?;
            for j in 0..n {
                if i == j { continue; }
                let denominator = concentration*diffusion[i*n+j];
                if !positive(denominator) { return Err(invalid("dusty-gas binary residual scaling is unrepresentable")); }
                residual.add((x[j]*reduced[i]-x[i]*reduced[j])/denominator)?;
                if j < i {
                    let yj = reduced[j]/root_x[j];
                    let drag = (root_x[j]*yi-root_x[i]*yj)/det::sqrt(denominator);
                    dissipation.add(R*drag*drag)?;
                }
            }
            report.maximum_gradient_residual_per_m = report.maximum_gradient_residual_per_m.max(residual.value.abs());
        }
        let viscous_entropy = -viscous_molar_flux*R*(pressure_difference/pressure)/distance;
        if !viscous_entropy.is_finite() || viscous_entropy < 0.0 {
            return Err(invalid("dusty-gas viscous dissipation is unrepresentable"));
        }
        dissipation.add(viscous_entropy)?;
        if report.maximum_gradient_residual_per_m > controls.gradient_tolerance_per_m {
            return Err(invalid("original dusty-gas equations did not close within the gradient tolerance"));
        }
        report.local_entropy_production_w_m3_k = dissipation.value;
        report.net_molar_flux_mol_m2_s = net.value;
        report.face_state = Some(face); report.darcy = Some(darcy);
        poll(&mut cancelled)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
