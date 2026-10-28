//! Source-backed isothermal/isobaric ideal-gas Maxwell-Stefan counterdiffusion.
//!
//! At the central face, solve sum_j (x_j J_i - x_i J_j)/D_ij = -c grad(x_i)
//! with sum_i J_i=0 (MOLAR-average frame, NOT zero mass flux). Use the existing
//! dense LU, then check EVERY original equation and the molar-frame constraint.
//! A finite-volume step delegates species/partial-enthalpy transport and full
//! endpoint reconstruction to exchange_species_with. Diffusivity is explicit
//! pair material DATA; viscosity, species labels and air defaults never infer it.
//!
//! This local closure neglects thermal/pressure diffusion and external species
//! forces. Caller-declared T/p difference limits admit that approximation; they
//! are not an error estimate. No bulk flow, Knudsen law or membrane selectivity.
//! Primary formulation: COMSOL, Multicomponent Gas Diffusion: Maxwell-Stefan
//! Description (equations 3-75--3-78); Cantera, MultiTransport reference.

/// All diffusion faces consume one unchanged material-network snapshot.
pub mod network;
/// Molecular friction, Knudsen wall transport and Darcy flow in one porous path.
pub mod dusty_gas;

use super::{ClosedGasParcel, GasEnergyError, GasSpeciesExchangeCandidate,
    GasSpeciesExchangeSolve, SpeciesId, SpeciesMassTransfer, Sum, catalog, invalid, poll, positive, source_error};
use crate::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::Dims;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use std::collections::BTreeMap;

/// Hard bound on the indivisible local dense factorization, not a mesh limit.
pub const MAX_DIFFUSION_SPECIES: usize = 32;
/// Molecular diffusivity has units m2/s, distinct from conductivity or viscosity.
pub const BINARY_DIFFUSIVITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);
pub const BINARY_DIFFUSIVITY_PROPERTY: &str = "binary_diffusivity";

/// One unordered physical species pair and its explicitly selected source law.
#[derive(Clone, Debug)]
pub struct BinaryDiffusionMaterial<'a> {
    pub species: [SpeciesId; 2],
    pub material: &'a MaterialCard,
    /// Other source axes are preserved; T/p are replaced by actual face values.
    pub context: &'a QueryPoint,
    pub selection: MaterialPropertySelection,
}
/// One-dimensional central face, with positive left -> right orientation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasDiffusionGeometry {
    pub area_m2: f64,
    pub distance_m: f64,
}
/// Physical approximation and numerical residual admission, all explicit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasDiffusionSolve {
    pub maximum_temperature_difference_k: f64,
    pub maximum_pressure_difference_pa: f64,
    /// Maximum absolute residual of each original Maxwell-Stefan row [1/m].
    pub gradient_tolerance_per_m: f64,
    /// Maximum absolute sum of molar fluxes [mol/(m2 s)].
    pub molar_flux_tolerance_mol_m2_s: f64,
    /// Maximum tolerated rounding decrease in actual pair entropy [J/K].
    pub entropy_tolerance_j_k: f64,
    pub exchange: GasSpeciesExchangeSolve,
}
impl GasDiffusionSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        self.exchange.validate()?;
        if ![self.maximum_temperature_difference_k, self.maximum_pressure_difference_pa,
            self.gradient_tolerance_per_m, self.molar_flux_tolerance_mol_m2_s,
            self.entropy_tolerance_j_k].iter().all(|v| positive(*v)) {
            return Err(invalid("positive finite diffusion support and residual tolerances required"));
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryDiffusionSample {
    /// Canonical unordered species pair.
    pub species: [SpeciesId; 2],
    pub diffusivity_m2_s: f64,
    pub material: ResolvedMaterialStatePoint,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesDiffusionFlux {
    pub species: SpeciesId,
    pub molar_flux_mol_m2_s: f64,
    /// Area-integrated signed rate, not kg/(m2 s).
    pub mass_rate_left_to_right_kg_s: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct GasDiffusionFlux {
    pub species: Vec<SpeciesDiffusionFlux>,
    pub pairs: Vec<BinaryDiffusionSample>,
    pub face_temperature_k: f64,
    pub face_pressure_pa: f64,
    pub maximum_gradient_residual_per_m: f64,
    pub molar_frame_residual_mol_m2_s: f64,
}
/// Original flux plus a checked, complete two-volume material proposal.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDiffusionCandidate {
    exchange: GasSpeciesExchangeCandidate,
    flux: GasDiffusionFlux,
    entropy_change_j_k: f64,
}
impl GasDiffusionCandidate {
    pub const fn exchange(&self) -> &GasSpeciesExchangeCandidate { &self.exchange }
    pub const fn flux(&self) -> &GasDiffusionFlux { &self.flux }
    pub const fn entropy_change_j_k(&self) -> f64 { self.entropy_change_j_k }
    pub fn into_parts(self) -> (GasSpeciesExchangeCandidate, GasDiffusionFlux, f64) {
        (self.exchange, self.flux, self.entropy_change_j_k)
    }
}
fn fractions(parcel: &ClosedGasParcel) -> BTreeMap<SpeciesId, f64> {
    parcel.model().thermodynamics().components().iter()
        .zip(parcel.model().thermodynamics().mole_composition().fractions())
        .map(|(s, &x)| (s.species().clone(), x)).collect()
}

impl ClosedGasParcel {
    /// Resolve the whole multicomponent diffusion system at unchanged gas states.
    /// Pair order and orientation do not affect the result. Exactly one source
    /// per active unordered pair is required (zero area is explicitly closed).
    /// The original source join refuses incompatible caloric/transport identities.
    /// LU is bounded to 32 species and polled immediately before/after; all
    /// original unscaled equations are checked, including the replaced gauge row.
    pub fn diffusion_flux_to(
        &self, right: &Self, geometry: GasDiffusionGeometry,
        pairs: &[BinaryDiffusionMaterial<'_>], controls: GasDiffusionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasDiffusionFlux, GasEnergyError> {
        controls.validate()?; poll(&mut cancelled)?;
        if !geometry.area_m2.is_finite() || geometry.area_m2 < 0.0 || !positive(geometry.distance_m)
            || pairs.len() > MAX_DIFFUSION_SPECIES*(MAX_DIFFUSION_SPECIES-1)/2 {
            return Err(invalid("bounded pair list, nonnegative diffusion area and positive distance required"));
        }
        let (a, b) = (self.state().state(), right.state().state());
        let temperature = 0.5*a.temperature + 0.5*b.temperature;
        let pressure = 0.5*a.pressure + 0.5*b.pressure;
        let mut report = GasDiffusionFlux { species: Vec::new(), pairs: Vec::new(),
            face_temperature_k: temperature, face_pressure_pa: pressure,
            maximum_gradient_residual_per_m: 0.0, molar_frame_residual_mol_m2_s: 0.0 };
        if geometry.area_m2 == 0.0 { poll(&mut cancelled)?; return Ok(report); }
        if (a.temperature-b.temperature).abs() > controls.maximum_temperature_difference_k
            || (a.pressure-b.pressure).abs() > controls.maximum_pressure_difference_pa {
            return Err(invalid("gas states exceed declared isothermal/isobaric diffusion support"));
        }
        let catalog = catalog(self, right)?;
        let components = catalog.thermodynamics().components();
        let n = components.len();
        if n > MAX_DIFFUSION_SPECIES || pairs.len() != n*(n-1)/2 {
            return Err(invalid("bounded diffusion requires exactly one material for each active pair"));
        }
        let ids: BTreeMap<_, _> = components.iter().enumerate().map(|(i, s)| (s.species().clone(), i)).collect();
        let mut ordered = BTreeMap::new();
        for pair in pairs {
            poll(&mut cancelled)?;
            let i = *ids.get(&pair.species[0]).ok_or(invalid("unknown diffusion pair species"))?;
            let j = *ids.get(&pair.species[1]).ok_or(invalid("unknown diffusion pair species"))?;
            if i == j || ordered.insert((i.min(j), i.max(j)), pair).is_some() {
                return Err(invalid("duplicate or self diffusion pair"));
            }
        }
        let (left_x, right_x) = (fractions(self), fractions(right));
        let mut x = Vec::with_capacity(n); let mut gradient = Vec::with_capacity(n);
        for species in components {
            let left = left_x.get(species.species()).copied().unwrap_or(0.0);
            let right = right_x.get(species.species()).copied().unwrap_or(0.0);
            x.push(0.5*left + 0.5*right);
            gradient.push((right-left)/geometry.distance_m);
        }
        let concentration = pressure/(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*temperature);
        if !positive(concentration) || gradient.iter().any(|g| !g.is_finite()) || x.iter().any(|v| !positive(*v)) {
            return Err(invalid("diffusion face concentration or gradient is unrepresentable"));
        }
        let requirement = ScalarPropertyRequirement::try_new(BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS, ScalarAdmissibility::StrictlyPositive).map_err(source_error)?;
        let mut diffusion = vec![0.0; n*n];
        let mut scale = f64::INFINITY;
        for ((i, j), binding) in ordered {
            poll(&mut cancelled)?;
            let point = super::super::updated_point(binding.context, temperature, pressure)?;
            let material = resolve_material_state_point(binding.material, &point,
                std::slice::from_ref(&requirement), binding.selection.clone()).map_err(source_error)?;
            let d = material.property(BINARY_DIFFUSIVITY_PROPERTY)
                .ok_or(invalid("resolved binary diffusivity missing"))?.value_si();
            diffusion[i*n+j] = d; diffusion[j*n+i] = d; scale = scale.min(d);
            report.pairs.push(BinaryDiffusionSample { species: [components[i].species().clone(),
                components[j].species().clone()], diffusivity_m2_s: d, material });
        }
        let mut flux = vec![0.0; n];
        if n > 1 {
            // D_min/D_ij scales the friction matrix without changing its solution.
            let mut matrix = vec![0.0; n*n];
            for i in 0..n {
                for j in 0..n {
                    if i != j {
                        let ratio = scale/diffusion[i*n+j];
                        if !positive(ratio) { return Err(invalid("diffusion coefficient ratio underflow")); }
                        matrix[i*n+i] += x[j]*ratio;
                        matrix[i*n+j] -= x[i]*ratio;
                    }
                }
                flux[i] = -(concentration*scale)*gradient[i];
            }
            // This constraint defines the molar-average frame, not a mass frame.
            matrix[(n-1)*n..].fill(1.0); flux[n-1] = 0.0;
            if matrix.iter().chain(&flux).any(|v| !v.is_finite()) {
                return Err(invalid("diffusion linear system is unrepresentable"));
            }
            poll(&mut cancelled)?;
            let factor = fs_la::factor::lu(&matrix, n).map_err(source_error)?;
            factor.solve(&mut flux);
            poll(&mut cancelled)?;
            if flux.iter().any(|v| !v.is_finite()) { return Err(invalid("diffusion flux solve overflow")); }
            for i in 0..n {
                let mut residual = Sum::default(); residual.add(gradient[i])?;
                for j in 0..n {
                    if i != j {
                        let denominator = concentration*diffusion[i*n+j];
                        if !positive(denominator) { return Err(invalid("diffusion residual scaling is unrepresentable")); }
                        residual.add((x[j]*flux[i]-x[i]*flux[j])/denominator)?;
                    }
                }
                report.maximum_gradient_residual_per_m = report.maximum_gradient_residual_per_m.max(residual.value.abs());
            }
        }
        let mut net = Sum::default();
        for (species, &j) in components.iter().zip(&flux) {
            net.add(j)?;
            let rate = (j*species.molar_mass().value())*geometry.area_m2;
            if !rate.is_finite() || (j != 0.0 && rate == 0.0) {
                return Err(invalid("diffusion species mass rate is unrepresentable"));
            }
            report.species.push(SpeciesDiffusionFlux { species: species.species().clone(),
                molar_flux_mol_m2_s: j, mass_rate_left_to_right_kg_s: rate });
        }
        report.molar_frame_residual_mol_m2_s = net.value;
        if report.maximum_gradient_residual_per_m > controls.gradient_tolerance_per_m
            || net.value.abs() > controls.molar_flux_tolerance_mol_m2_s {
            return Err(invalid("diffusion original equations or molar-frame constraint did not close"));
        }
        poll(&mut cancelled)?;
        Ok(report)
    }

    /// Explicit first-order diffusion from one original pair snapshot. Partial
    /// species enthalpy crosses once; no mixture-enthalpy surrogate or extra heat.
    /// Returns no endpoint on source/budget/consumption/energy/entropy refusal.
    /// The whole-pair entropy check cancels conserved species reference datums;
    /// numerical species-mass residuals remain in the separate mass accounts.
    /// This is not an unconditional stability theorem. Shorten the physical
    /// interval rather than clipping trace species.
    pub fn diffuse_with(
        &self, right: &Self, geometry: GasDiffusionGeometry, duration_s: f64,
        pairs: &[BinaryDiffusionMaterial<'_>], controls: GasDiffusionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasDiffusionCandidate, GasEnergyError> {
        if !positive(duration_s) { return Err(invalid("positive finite diffusion duration required")); }
        let flux = self.diffusion_flux_to(right, geometry, pairs, controls, &mut cancelled)?;
        let mut rows = Vec::with_capacity(flux.species.len());
        for species in &flux.species {
            let mass = species.mass_rate_left_to_right_kg_s*duration_s;
            if !mass.is_finite() || (species.mass_rate_left_to_right_kg_s != 0.0 && mass == 0.0) {
                return Err(invalid("integrated diffusion mass vanished or overflowed"));
            }
            rows.push(SpeciesMassTransfer { species: species.species.clone(), mass_left_to_right_kg: mass });
        }
        let exchange = self.exchange_species_with(right, &rows, controls.exchange, &mut cancelled)?;
        let entropy = super::entropy::closed_entropy_change(
            [(0, self, exchange.left()), (1, right, exchange.right())],
            controls.exchange.mass_tolerance_kg, &mut cancelled,
        ).map_err(|error| match error {
            super::entropy::Error::Invalid(reason) => invalid(reason),
            super::entropy::Error::Parcel { error, .. } => error,
            super::entropy::Error::Cancelled => GasEnergyError::Cancelled,
        })?;
        if entropy < -controls.entropy_tolerance_j_k {
            return Err(invalid("finite diffusion step decreases pair entropy; reduce the interval"));
        }
        poll(&mut cancelled)?;
        Ok(GasDiffusionCandidate { exchange, flux, entropy_change_j_k: entropy })
    }
}

#[cfg(test)]
mod tests;
