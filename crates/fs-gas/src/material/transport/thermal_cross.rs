//! Original material claims for the shared reciprocal thermal-diffusion law.

use super::{admit_state, material_input, poll, typed};
use crate::multicomponent::dissipation::{
    DissipationError, TransportCoefficients, TransportLaw, validate_binary,
};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{MassSpecificEnergyQuantityV1, SpeciesId};

/// Explicit convention of the wrapped law's thermal conductivity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalConductivityConvention {
    /// Conductivity at zero diffusive species flux, excluding reciprocal heat.
    /// A measured conductivity at another flux constraint cannot be used here
    /// without a separately justified conversion by its source owner.
    ZeroDiffusiveSpeciesFlux,
}

/// One explicitly assigned kinetic heat-of-transport source, not a caloric law.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesHeatOfTransportMaterial {
    /// Original canonical species; neither material nor property names infer it.
    pub species: SpeciesId,
    /// Immutable source card supporting this species' kinetic cross coefficient.
    pub material: MaterialCard,
    /// Explicit finite signed mass-specific energy [J/kg] property assignment.
    /// A partial enthalpy is NOT this coefficient, despite identical dimensions.
    pub property: PropertyKey,
    /// Source selection applied once at the reference state, then pinned.
    pub selection: MaterialPropertySelection,
}

/// One complete successful thermal-cross query in the original species order.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMaterialHeatsOfTransport {
    temperature: Temperature,
    pressure: Pressure,
    requested_fractions: Vec<f64>,
    mole_fractions: Vec<f64>,
    heats: Vec<MassSpecificEnergyQuantityV1>,
    components: Vec<ResolvedMaterialStatePoint>,
}
impl ResolvedMaterialHeatsOfTransport {
    /// Actual common face temperature.
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Actual common face pressure.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Same canonical composition supplied to the wrapped molecular law.
    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }
    /// Signed kinetic cross coefficients, not partial thermodynamic enthalpies.
    pub fn heats(&self) -> &[MassSpecificEnergyQuantityV1] {
        &self.heats
    }
    /// Complete selected claims, typed coordinates, cards and usage receipts.
    pub fn components(&self) -> &[ResolvedMaterialStatePoint] {
        &self.components
    }

    fn matches(&self, t: Temperature, p: Pressure, x: &[f64]) -> bool {
        self.temperature == t && self.pressure == p && self.requested_fractions == x
    }
}

/// Material-sourced Soret/Dufour coefficients for every shared gas-face owner.
///
/// Wrap `MaterialTransport` to combine source-card viscosity/conductivity with
/// source-card cross coefficients. Other `TransportLaw` implementations retain
/// their original molecular ownership. The explicitly declared conductivity
/// convention applies to the wrapped law's actual output, not to unused fields
/// in an inner binary-diffusivity provider.
///
/// The reference query selects exact source claims. Every subsequent query
/// refreshes those same claims at actual T,p,composition; it cannot switch to a
/// different claim when a source boundary is crossed. Scalar, curve and surface
/// interpolation remain with the original material resolver. Every species,
/// including exact zeros, has a declared coefficient and query coordinate.
/// No coefficient is inferred from caloric data, names or collision parameters.
///
/// Runtime queries allocate through the original material resolver. Repeated
/// `evaluate`/`heats_of_transport` calls for the same face reuse one successful
/// cross query. Failure/cancellation preserves outputs and the last successful
/// cross receipts; private workspaces are not accepted fluid state.
pub struct MaterialSoretDufourTransport<L: TransportLaw> {
    law: L,
    convention: ThermalConductivityConvention,
    sources: Vec<SpeciesHeatOfTransportMaterial>,
    composition_axes: Vec<(SpeciesId, String)>,
    reference_query: QueryPoint,
    resolved: ResolvedMaterialHeatsOfTransport,
    binary_workspace: Vec<f64>,
}
impl<L: TransportLaw> MaterialSoretDufourTransport<L> {
    /// Bind complete canonical species assignments and the original typed query.
    ///
    /// `reference_query` must contain `temperature` [absolute K], `pressure`
    /// [Pa], and one explicitly typed mole fraction for every species/axis pair.
    /// Input assignment order is immaterial. Additional external coordinates
    /// are retained; aliases for the evolving T,p,composition coordinates refuse
    /// rather than silently freezing a source coordinate at its initial value.
    pub fn new(
        law: L,
        convention: ThermalConductivityConvention,
        mut sources: Vec<SpeciesHeatOfTransportMaterial>,
        mut composition_axes: Vec<(SpeciesId, String)>,
        reference_query: QueryPoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        let models = law.source().models();
        let n = models.len();
        sources.sort_by(|a, b| a.species.cmp(&b.species));
        composition_axes.sort_by(|a, b| a.0.cmp(&b.0));
        if sources.len() != n
            || composition_axes.len() != n
            || sources
                .iter()
                .zip(models)
                .any(|(s, m)| &s.species != m.species())
            || composition_axes
                .iter()
                .zip(models)
                .any(|((s, _), m)| s != m.species())
            || composition_axes.iter().enumerate().any(|(i, (_, axis))| {
                axis.trim().is_empty()
                    || axis.len() > 128
                    || axis == "temperature"
                    || axis == "pressure"
                    || composition_axes[..i].iter().any(|(_, prior)| prior == axis)
            })
        {
            return Err(DissipationError::Invalid(
                "thermal cross source species or coordinate axis",
            ));
        }
        let temperature = Temperature::new(coordinate(
            &reference_query,
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
        )?);
        let pressure = Pressure::new(coordinate(
            &reference_query,
            "pressure",
            typed(QuantityKind::Pressure),
        )?);
        let mole_quantity = typed(QuantityKind::Composition(CompositionBasis::MoleFraction));
        let fractions = composition_axes
            .iter()
            .map(|(_, axis)| coordinate(&reference_query, axis, mole_quantity))
            .collect::<Result<Vec<_>, _>>()?;
        for (axis, quantity) in reference_query.axis_quantities() {
            if (quantity.dims() == Temperature::DIMS && axis != "temperature")
                || (quantity.dims() == Pressure::DIMS && axis != "pressure")
                || (quantity
                    .semantic_type()
                    .is_some_and(|q| matches!(q.kind(), QuantityKind::Composition(_)))
                    && !composition_axes.iter().any(|(_, own)| own == axis))
            {
                return Err(DissipationError::Invalid(
                    "unowned evolving thermal cross coordinate",
                ));
            }
        }
        let (query, canonical) = face_query(
            &reference_query,
            &composition_axes,
            temperature,
            pressure,
            &fractions,
        )?;
        let mut components = Vec::with_capacity(n);
        let mut heats = Vec::with_capacity(n);
        for source in &sources {
            poll(&mut cancelled)?;
            if source.property.quantity()
                != QuantitySpec::dimensional(MassSpecificEnergyQuantityV1::DIMS)
            {
                return Err(DissipationError::Invalid(
                    "heat of transport requires signed J/kg quantity",
                ));
            }
            let requirement = ScalarPropertyRequirement::try_with_key(
                &source.property,
                ScalarAdmissibility::Finite,
            )
            .map_err(|e| material_input("heat of transport requirement", e))?;
            let resolved = resolve_material_state_point(
                &source.material,
                &query,
                &[requirement],
                source.selection.clone(),
            )
            .map_err(|e| material_input("heat of transport source", e))?;
            heats.push(heat(source, &resolved)?);
            components.push(resolved);
        }
        let count = n
            .checked_mul(n)
            .ok_or(DissipationError::Invalid("binary diffusion size overflow"))?;
        let mut binary_workspace = Vec::new();
        binary_workspace
            .try_reserve_exact(count)
            .map_err(|_| DissipationError::Flow(crate::multicomponent::Error::Allocation))?;
        binary_workspace.resize(count, 0.0);
        poll(&mut cancelled)?;
        Ok(Self {
            law,
            convention,
            sources,
            composition_axes,
            reference_query,
            resolved: ResolvedMaterialHeatsOfTransport {
                temperature,
                pressure,
                requested_fractions: fractions,
                mole_fractions: canonical,
                heats,
                components,
            },
            binary_workspace,
        })
    }

    /// Original molecular law, including its own material/source assignments.
    pub const fn law(&self) -> &L {
        &self.law
    }
    /// Explicit conductivity convention used by the reciprocal closure.
    pub const fn conductivity_convention(&self) -> ThermalConductivityConvention {
        self.convention
    }
    /// Original complete kinetic source assignments in canonical species order.
    pub fn sources(&self) -> &[SpeciesHeatOfTransportMaterial] {
        &self.sources
    }
    /// Canonical species-to-query-axis mapping; no species name is interpreted.
    pub fn composition_axes(&self) -> &[(SpeciesId, String)] {
        &self.composition_axes
    }
    /// Authored reference query, including unchanged external conditions.
    pub const fn reference_query(&self) -> &QueryPoint {
        &self.reference_query
    }
    /// Last complete successful query, initially the source-binding reference.
    pub const fn last_resolved(&self) -> &ResolvedMaterialHeatsOfTransport {
        &self.resolved
    }

    fn resolve(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<ResolvedMaterialHeatsOfTransport, DissipationError> {
        poll(cancelled)?;
        let (query, canonical) = face_query(
            &self.reference_query,
            &self.composition_axes,
            temperature,
            pressure,
            fractions,
        )?;
        let mut components = Vec::with_capacity(self.sources.len());
        let mut heats = Vec::with_capacity(self.sources.len());
        for (source, original) in self.sources.iter().zip(&self.resolved.components) {
            poll(cancelled)?;
            let resolved = original
                .at_query_point(&source.material, &query)
                .map_err(|e| material_input("heat of transport refresh", e))?;
            heats.push(heat(source, &resolved)?);
            components.push(resolved);
        }
        poll(cancelled)?;
        Ok(ResolvedMaterialHeatsOfTransport {
            temperature,
            pressure,
            requested_fractions: fractions.to_vec(),
            mole_fractions: canonical,
            heats,
            components,
        })
    }
}
impl<L: TransportLaw> TransportLaw for MaterialSoretDufourTransport<L> {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.law.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary_diffusivities: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        poll(cancelled)?;
        if binary_diffusivities.len() != self.binary_workspace.len() {
            return Err(DissipationError::Invalid(
                "thermal cross binary output axis",
            ));
        }
        let candidate = if self.resolved.matches(temperature, pressure, fractions) {
            None
        } else {
            Some(self.resolve(temperature, pressure, fractions, cancelled)?)
        };
        let resolved = candidate.as_ref().unwrap_or(&self.resolved);
        let coefficients = self.law.evaluate(
            temperature,
            pressure,
            &resolved.mole_fractions,
            &mut self.binary_workspace,
            cancelled,
        )?;
        validate_binary(&self.binary_workspace, self.sources.len())?;
        poll(cancelled)?;
        binary_diffusivities.copy_from_slice(&self.binary_workspace);
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
        output_j_kg: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        poll(cancelled)?;
        if output_j_kg.len() != self.sources.len() {
            return Err(DissipationError::Invalid("thermal cross heat output axis"));
        }
        let candidate = if self.resolved.matches(temperature, pressure, fractions) {
            None
        } else {
            Some(self.resolve(temperature, pressure, fractions, cancelled)?)
        };
        let resolved = candidate.as_ref().unwrap_or(&self.resolved);
        poll(cancelled)?;
        for (output, heat) in output_j_kg.iter_mut().zip(&resolved.heats) {
            *output = heat.value();
        }
        if let Some(candidate) = candidate {
            self.resolved = candidate;
        }
        Ok(true)
    }
}

fn face_query(
    reference: &QueryPoint,
    axes: &[(SpeciesId, String)],
    t: Temperature,
    p: Pressure,
    x: &[f64],
) -> Result<(QueryPoint, Vec<f64>), DissipationError> {
    admit_state(axes.len(), t, p, x)?;
    let mut canonical = x.to_vec();
    fs_thermal::enthalpy::advection::species::close_generated_simplex(&mut canonical)
        .map_err(|e| material_input("heat of transport composition", e))?;
    let mut query = reference
        .clone()
        .with("temperature", t.value())
        .and_then(|q| q.with("pressure", p.value()))
        .map_err(|e| material_input("heat of transport coordinates", e))?;
    for ((_, axis), fraction) in axes.iter().zip(&canonical) {
        query = query
            .with(axis, *fraction)
            .map_err(|e| material_input("heat of transport composition coordinate", e))?;
    }
    Ok((query, canonical))
}
fn coordinate(
    query: &QueryPoint,
    axis: &str,
    quantity: QuantitySpec,
) -> Result<f64, DissipationError> {
    if query.axis_quantities().get(axis) != Some(&quantity) {
        return Err(DissipationError::Invalid(
            "thermal cross coordinate quantity",
        ));
    }
    query
        .axes()
        .get(axis)
        .copied()
        .ok_or(DissipationError::Invalid(
            "missing thermal cross coordinate",
        ))
}
fn heat(
    source: &SpeciesHeatOfTransportMaterial,
    resolved: &ResolvedMaterialStatePoint,
) -> Result<MassSpecificEnergyQuantityV1, DissipationError> {
    let property = resolved
        .property_by_key(&source.property)
        .ok_or(DissipationError::Invalid(
            "missing resolved heat of transport",
        ))?;
    if !property.has_only_fixed_external_conditions(Some("temperature"), false) {
        return Err(DissipationError::Invalid(
            "heat of transport requires static local T,p,composition and external-condition support",
        ));
    }
    Ok(MassSpecificEnergyQuantityV1::new(property.value_si()))
}
