//! Original mixture bulk viscosity at the actual compressing gas state.
//!
//! This is the Newtonian coefficient zeta in tau = 2 mu dev(D) + zeta tr(D) I.
//! The original molecular face owns both stress and its conservative energy
//! flux; adding a separate viscous heat source would count that work twice.
//! A source must support a local, instantaneous coefficient. A frequency-
//! dependent apparent viscosity is not an admissible replacement for molecular
//! relaxation memory, and no mixture rule is inferred from pure-species data.
//! Callers composing a separate relaxation channel must ensure the supplied
//! coefficient does not include that same loss a second time.

use super::thermal_cross::{face_query, reference_state};
use super::{material_input, poll};
use crate::multicomponent::dissipation::{
    DissipationError, TransportCoefficients, TransportLaw, validate_binary,
};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;

/// Explicit assignment of a mixture's bulk-viscosity source in Pa s.
///
/// This coefficient is distinct from shear viscosity despite equal dimensions.
/// The caller supplies the physical assignment; neither a property name nor a
/// caloric/collision source implies the constitutive law or its calibration.
#[derive(Clone, Debug, PartialEq)]
pub struct BulkViscosityMaterial {
    pub material: MaterialCard,
    pub property: PropertyKey,
    /// Applied at the reference state, then the exact selected claim is pinned.
    pub selection: MaterialPropertySelection,
}

/// One complete successful original-source query and its composition.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBulkViscosity {
    temperature: Temperature,
    pressure: Pressure,
    requested_fractions: Vec<f64>,
    mole_fractions: Vec<f64>,
    bulk_viscosity: DynViscosity,
    material: ResolvedMaterialStatePoint,
}
impl ResolvedBulkViscosity {
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Canonical fractions actually supplied to every constitutive owner.
    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }
    /// Nonnegative source coefficient, including an explicitly supported zero.
    pub const fn bulk_viscosity(&self) -> DynViscosity {
        self.bulk_viscosity
    }
    /// Original card, pinned claim, full query and usage receipts.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    fn matches(&self, t: Temperature, p: Pressure, x: &[f64]) -> bool {
        self.temperature == t && self.pressure == p && self.requested_fractions == x
    }
}

/// Material bulk viscosity composed with any existing gas transport law.
///
/// Only the wrapped law's bulk-viscosity output is replaced. Shear viscosity,
/// conductivity, diffusion and reciprocal thermal cross coefficients retain
/// their original owners. Thus the same material source reaches the existing
/// axial, planar, spatial, moving-wall and reacting compression paths.
///
/// Every query uses actual T,p and the complete canonical mole composition,
/// with unchanged external coordinates. Source boundaries refuse instead of
/// selecting another claim or keeping an old coefficient. Repeated queries at
/// exactly the same state reuse one successful source resolution. Material
/// queries may allocate; the adapter does not claim allocation-free execution.
/// Failure/cancellation preserves caller outputs and the last complete bulk
/// receipt. Wrapped workspaces/receipts are not accepted evolving gas state.
pub struct MaterialBulkViscosityTransport<L: TransportLaw> {
    law: L,
    bulk_source: BulkViscosityMaterial,
    composition_axes: Vec<(SpeciesId, String)>,
    reference_query: QueryPoint,
    resolved: ResolvedBulkViscosity,
    binary_workspace: Vec<f64>,
    heat_workspace: Vec<f64>,
}
impl<L: TransportLaw> MaterialBulkViscosityTransport<L> {
    /// Bind a source and one typed mole-fraction coordinate for every species.
    ///
    /// Reference T,p must be absolute K and Pa under `temperature`/`pressure`.
    /// Axis assignment order is immaterial. Aliases of evolving coordinates,
    /// incomplete species axes, nonstatic source support and extrapolation
    /// refuse through the original material/transport admission rules.
    pub fn new(
        law: L,
        bulk_source: BulkViscosityMaterial,
        mut composition_axes: Vec<(SpeciesId, String)>,
        reference_query: QueryPoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        let (temperature, pressure, fractions) =
            reference_state(law.source(), &mut composition_axes, &reference_query)?;
        let (query, canonical) = face_query(
            &reference_query,
            &composition_axes,
            temperature,
            pressure,
            &fractions,
        )?;
        if bulk_source.property.quantity() != QuantitySpec::dimensional(DynViscosity::DIMS) {
            return Err(DissipationError::Invalid(
                "bulk viscosity requires Pa s quantity",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_with_key(
            &bulk_source.property,
            ScalarAdmissibility::NonNegative,
        )
        .map_err(|e| material_input("bulk viscosity requirement", e))?;
        let material = resolve_material_state_point(
            &bulk_source.material,
            &query,
            &[requirement],
            bulk_source.selection.clone(),
        )
        .map_err(|e| material_input("bulk viscosity source", e))?;
        let bulk_viscosity = coefficient(&bulk_source, &material)?;
        let n = canonical.len();
        let count = n
            .checked_mul(n)
            .ok_or(DissipationError::Invalid("binary diffusion size overflow"))?;
        let allocate = |len| {
            let mut workspace = Vec::new();
            workspace
                .try_reserve_exact(len)
                .map_err(|_| DissipationError::Flow(crate::multicomponent::Error::Allocation))?;
            workspace.resize(len, 0.0);
            Ok::<_, DissipationError>(workspace)
        };
        let binary_workspace = allocate(count)?;
        let heat_workspace = allocate(n)?;
        poll(&mut cancelled)?;
        Ok(Self {
            law,
            bulk_source,
            composition_axes,
            reference_query,
            resolved: ResolvedBulkViscosity {
                temperature,
                pressure,
                requested_fractions: fractions,
                mole_fractions: canonical,
                bulk_viscosity,
                material,
            },
            binary_workspace,
            heat_workspace,
        })
    }

    /// Original law and all of its source assignments.
    pub const fn law(&self) -> &L {
        &self.law
    }
    pub const fn bulk_source(&self) -> &BulkViscosityMaterial {
        &self.bulk_source
    }
    pub fn composition_axes(&self) -> &[(SpeciesId, String)] {
        &self.composition_axes
    }
    pub const fn reference_query(&self) -> &QueryPoint {
        &self.reference_query
    }
    pub const fn last_resolved(&self) -> &ResolvedBulkViscosity {
        &self.resolved
    }

    fn candidate(
        &self,
        t: Temperature,
        p: Pressure,
        fractions: &[f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Option<ResolvedBulkViscosity>, DissipationError> {
        poll(cancelled)?;
        if self.resolved.matches(t, p, fractions) {
            return Ok(None);
        }
        let (query, canonical) = face_query(
            &self.reference_query,
            &self.composition_axes,
            t,
            p,
            fractions,
        )?;
        let material = self
            .resolved
            .material
            .at_query_point(&self.bulk_source.material, &query)
            .map_err(|e| material_input("bulk viscosity refresh", e))?;
        let bulk_viscosity = coefficient(&self.bulk_source, &material)?;
        poll(cancelled)?;
        Ok(Some(ResolvedBulkViscosity {
            temperature: t,
            pressure: p,
            requested_fractions: fractions.to_vec(),
            mole_fractions: canonical,
            bulk_viscosity,
            material,
        }))
    }
}

impl<L: TransportLaw> TransportLaw for MaterialBulkViscosityTransport<L> {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.law.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        output: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        poll(cancelled)?;
        if output.len() != self.binary_workspace.len() {
            return Err(DissipationError::Invalid(
                "bulk viscosity binary output axis",
            ));
        }
        let candidate = self.candidate(temperature, pressure, fractions, cancelled)?;
        let resolved = candidate.as_ref().unwrap_or(&self.resolved);
        let mut coefficients = self.law.evaluate(
            temperature,
            pressure,
            &resolved.mole_fractions,
            &mut self.binary_workspace,
            cancelled,
        )?;
        validate_binary(&self.binary_workspace, self.heat_workspace.len())?;
        coefficients.bulk_viscosity = resolved.bulk_viscosity;
        poll(cancelled)?;
        output.copy_from_slice(&self.binary_workspace);
        if let Some(candidate) = candidate {
            self.resolved = candidate;
        }
        Ok(coefficients)
    }

    fn heats_of_transport(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        output: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        poll(cancelled)?;
        if output.len() != self.heat_workspace.len() {
            return Err(DissipationError::Invalid("bulk viscosity heat output axis"));
        }
        let candidate = self.candidate(temperature, pressure, fractions, cancelled)?;
        let resolved = candidate.as_ref().unwrap_or(&self.resolved);
        self.heat_workspace.fill(0.0);
        let present = self.law.heats_of_transport(
            temperature,
            pressure,
            &resolved.mole_fractions,
            &mut self.heat_workspace,
            cancelled,
        )?;
        if present && self.heat_workspace.iter().any(|v| !v.is_finite()) {
            return Err(DissipationError::Invalid(
                "nonfinite delegated heat of transport",
            ));
        }
        if !present {
            self.heat_workspace.fill(0.0);
        }
        poll(cancelled)?;
        output.copy_from_slice(&self.heat_workspace);
        if let Some(candidate) = candidate {
            self.resolved = candidate;
        }
        Ok(present)
    }
}

fn coefficient(
    source: &BulkViscosityMaterial,
    resolved: &ResolvedMaterialStatePoint,
) -> Result<DynViscosity, DissipationError> {
    let property = resolved
        .property_by_key(&source.property)
        .ok_or(DissipationError::Invalid("missing resolved bulk viscosity"))?;
    if !property.has_only_fixed_external_conditions(Some("temperature"), false) {
        return Err(DissipationError::Invalid(
            "bulk viscosity requires static local T,p,composition and external-condition support",
        ));
    }
    Ok(DynViscosity::new(property.value_si()))
}
