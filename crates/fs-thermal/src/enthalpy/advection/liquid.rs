//! Liquid enthalpy transport with a stationary pure-substance solid skeleton.
//!
//! Storage remains the original mixture `mass*h`. The moving liquid carries
//! the liquid-side plateau enthalpy in coexistence, and the original h when
//! fully liquid. Thus balanced isothermal liquid circulation cannot transport
//! a stationary solid's latent-energy fraction. This supplies no momentum,
//! phase-volume change, species segregation or second heat integrator.

use super::{Advected, MassBoundary, finite, invalid};
use crate::enthalpy::{
    DepositedCell, DepositionError, DepositionSimulation, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};

/// One original source-chart knot, including its phase and density data.
/// The caller supplies the complete immutable chart, not selected samples or
/// inferred phase fractions. Its caloric identity remains the binding authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidCarrierKnot {
    pub specific_enthalpy_j_kg: f64,
    pub temperature_k: f64,
    pub liquid_mass_fraction: f64,
    pub bulk_density_kg_m3: f64,
}

/// Validated isothermal phase chart for liquid-only mass carriage.
///
/// One fixed density applies throughout the original chart. Source receipts,
/// pressure/composition context and the validity of the stationary-skeleton
/// approximation belong to its caller; this descriptor invents none of them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidEnthalpyCarrier {
    material: MaterialId,
    caloric_identity: [u8; 32],
    enthalpy_bounds_j_kg: [f64; 2],
    plateau_enthalpy_bounds_j_kg: [f64; 2],
    melting_temperature_k: f64,
    density_kg_m3: f64,
}

impl LiquidEnthalpyCarrier {
    /// Bind the complete original chart to its native enthalpy law.
    ///
    /// Requires an isothermal solid/liquid plateau, affine liquid mass fraction
    /// on that plateau, strictly sensible pure branches, and equal constant
    /// phase densities. Each supplied knot's T(h) is checked against the law;
    /// its retained identity is checked again at every source/residual epoch.
    /// Passing a subset or an independently authored phase chart is not a
    /// source binding. Material adapters should pass their original knots.
    pub fn from_source_chart(
        material: MaterialId,
        law: &dyn ThermalLaw,
        knots: &[LiquidCarrierKnot],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if !(2..=65_536).contains(&knots.len())
            || law.storage_basis() != ThermalStorageBasis::Enthalpy
        {
            return Err(invalid("liquid carriage needs a bounded complete enthalpy phase chart"));
        }
        let caloric_identity = law.caloric_identity()
            .ok_or(invalid("liquid carriage requires the original caloric identity"))?;
        let rho = knots[0].bulk_density_kg_m3;
        for (i, knot) in knots.iter().enumerate() {
            if i % 256 == 0 { checkpoint()?; }
            if !knot.specific_enthalpy_j_kg.is_finite()
                || !knot.temperature_k.is_finite() || knot.temperature_k <= 0.0
                || !knot.liquid_mass_fraction.is_finite()
                || !(0.0..=1.0).contains(&knot.liquid_mass_fraction)
                || !knot.bulk_density_kg_m3.is_finite() || knot.bulk_density_kg_m3 <= 0.0
                || knot.bulk_density_kg_m3.to_bits() != rho.to_bits()
                || law.temperature_k(knot.specific_enthalpy_j_kg)? != knot.temperature_k
            {
                return Err(invalid("liquid carriage requires original phase knots and equal constant density"));
            }
            if i != 0 {
                let previous = knots[i - 1];
                if knot.specific_enthalpy_j_kg <= previous.specific_enthalpy_j_kg
                    || knot.temperature_k < previous.temperature_k
                    || knot.liquid_mass_fraction < previous.liquid_mass_fraction
                    || (knot.liquid_mass_fraction == previous.liquid_mass_fraction
                        && knot.temperature_k <= previous.temperature_k)
                {
                    return Err(invalid("liquid carriage requires ordered sensible branches and one latent plateau"));
                }
            }
        }
        if knots[0].liquid_mass_fraction != 0.0
            || knots[knots.len() - 1].liquid_mass_fraction != 1.0
        {
            return Err(invalid("liquid carriage requires complete solid and liquid phase endpoints"));
        }
        let first_mobile = knots.iter().position(|k| k.liquid_mass_fraction > 0.0)
            .ok_or(invalid("liquid carriage lacks a liquid phase"))?;
        let solid = knots[first_mobile - 1];
        let liquid_index = knots.iter().position(|k| k.liquid_mass_fraction == 1.0)
            .ok_or(invalid("liquid carriage lacks a pure liquid endpoint"))?;
        let liquid = knots[liquid_index];
        let latent = finite(liquid.specific_enthalpy_j_kg - solid.specific_enthalpy_j_kg)?;
        if latent <= 0.0 {
            return Err(invalid("liquid carriage requires positive representable latent enthalpy"));
        }
        for (i, knot) in knots[first_mobile - 1..=liquid_index].iter().enumerate() {
            if i % 256 == 0 { checkpoint()?; }
            let fraction = finite(knot.specific_enthalpy_j_kg - solid.specific_enthalpy_j_kg)? / latent;
            if knot.temperature_k != solid.temperature_k
                || (fraction - knot.liquid_mass_fraction).abs() > 64.0 * f64::EPSILON
            {
                return Err(invalid("liquid carriage requires a pure-substance isothermal affine latent plateau"));
            }
        }
        checkpoint()?;
        Ok(Self {
            material,
            caloric_identity,
            enthalpy_bounds_j_kg: [knots[0].specific_enthalpy_j_kg, knots[knots.len() - 1].specific_enthalpy_j_kg],
            plateau_enthalpy_bounds_j_kg: [solid.specific_enthalpy_j_kg, liquid.specific_enthalpy_j_kg],
            melting_temperature_k: solid.temperature_k,
            density_kg_m3: rho,
        })
    }

    pub const fn material_id(&self) -> MaterialId { self.material }
    pub const fn caloric_identity(&self) -> [u8; 32] { self.caloric_identity }
    pub const fn plateau_enthalpy_bounds_j_kg(&self) -> [f64; 2] { self.plateau_enthalpy_bounds_j_kg }
    pub const fn melting_temperature_k(&self) -> f64 { self.melting_temperature_k }
    pub const fn density_kg_m3(&self) -> f64 { self.density_kg_m3 }

    pub(super) fn admit_law(&self, thermal: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        if thermal.material_storage_basis(self.material)? != ThermalStorageBasis::Enthalpy
            || thermal.material_caloric_identity(self.material)? != Some(self.caloric_identity)
        {
            return Err(invalid("liquid-carrier caloric identity or storage basis changed"));
        }
        Ok(())
    }

    pub(super) fn admit_boundary(&self, boundary: &MassBoundary) -> Result<(), DepositionError> {
        if let MassBoundary::Inlet { material_id, specific_h_j_kg: h, .. } = *boundary {
            if material_id != self.material {
                return Err(invalid("liquid-carrier inlet must use the original material and enthalpy datum"));
            }
            if !h.is_finite() || h < self.enthalpy_bounds_j_kg[0] || h > self.enthalpy_bounds_j_kg[1] {
                return Err(DepositionError::MaterialDomain);
            }
            if h < self.plateau_enthalpy_bounds_j_kg[1] {
                return Err(invalid("liquid-carrier inlet must be entirely liquid within its original phase chart"));
            }
        }
        Ok(())
    }

    pub(super) fn carried_h(&self, cell: &DepositedCell) -> Result<f64, DepositionError> {
        let expected_mass = finite(self.density_kg_m3 * cell.geometry.volume_m3)?;
        if cell.material_id != self.material || expected_mass <= 0.0
            || (cell.mass_kg - expected_mass).abs()
                > 64.0 * f64::EPSILON * cell.mass_kg.max(expected_mass)
        {
            return Err(invalid("liquid carriage changed original material or fixed-density inventory"));
        }
        let h = cell.specific_h_j_kg;
        if !h.is_finite() || h < self.enthalpy_bounds_j_kg[0] || h > self.enthalpy_bounds_j_kg[1] {
            return Err(DepositionError::MaterialDomain);
        }
        if h <= self.plateau_enthalpy_bounds_j_kg[0] {
            return Err(invalid("nonzero liquid mass flow crosses a fully solid cell"));
        }
        Ok(h.max(self.plateau_enthalpy_bounds_j_kg[1]))
    }
}

impl<O> Advected<O> {
    /// Transport the moving liquid's h while retaining stationary-solid energy
    /// in the original bulk cell. Every cell balances liquid mass, including
    /// any explicit inlets/outlets. The resolved mass rate is liquid mass per
    /// second, not total mixture motion or intrinsic liquid velocity; its owner
    /// supplies momentum. The solid skeleton carries no mass through a boundary.
    ///
    /// An inlet must be pure liquid at or above the original liquidus enthalpy,
    /// within that same material chart and enthalpy datum. An outlet carries the
    /// current liquid enthalpy, including in a partially solid cell. Boundary
    /// cells must contain liquid. The original external-enthalpy account owns
    /// their signed heat input; no separate latent, pressure or kinetic source
    /// is inferred. Source pressure/work conventions remain the flow owner's
    /// responsibility, just as for ordinary single-phase stream enthalpy.
    ///
    /// Existing conservative transfers, component-reference subtraction and
    /// physical energy gates remain authoritative. The carried-h slope is in
    /// [0,1], so the original outgoing-mass CFL remains conservative. A bulk-h
    /// entropy audit, species mixing or nonempty bulk-h reconstruction refuses.
    pub fn with_liquid_carrier(mut self, carrier: LiquidEnthalpyCarrier) -> Result<Self, DepositionError> {
        if self.stationary_species.is_some() || !self.selective_outlets.is_empty()
            || self.entropy.is_some()
        {
            return Err(invalid("liquid-carrier transport excludes species and bulk-entropy mixing"));
        }
        for boundary in &self.boundaries { carrier.admit_boundary(boundary)?; }
        self.liquid_carrier = Some(carrier);
        self.pending = None;
        self.boundary_pending = None;
        Ok(self)
    }

    /// Original phase descriptor, absent for homogeneous mixture transport.
    pub fn liquid_carrier(&self) -> Option<&LiquidEnthalpyCarrier> {
        self.liquid_carrier.as_ref()
    }
}

#[cfg(test)]
mod tests;
