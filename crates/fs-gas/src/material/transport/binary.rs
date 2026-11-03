//! Original binary-diffusivity material cards in every continuum gas consumer.

use super::{MaterialTransport, admit_state, material_input, poll, typed};
use crate::multicomponent::dissipation::{
    DissipationError, TransportCoefficients, TransportLaw, validate_binary,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
    BinaryDiffusionSample,
};
use fs_qty::semantic::QuantityKind;
use fs_qty::{DynViscosity, Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

/// Complete material-card molecular law for the existing gas-face owners.
///
/// The original material gas supplies Wilke/WMS viscosity and conductivity;
/// explicitly assigned pair cards supply the canonical symmetric `D_ij` matrix
/// in m²/s; bulk viscosity remains an explicit nonnegative input. This consumes
/// the same [`BinaryDiffusionMaterial`] assignments as material parcel diffusion
/// without retaining its isothermal/isobaric or molar-frame restrictions: the
/// consuming continuum face owns pressure drive, the mass frame and species
/// enthalpy transport. [`super::MaterialSoretDufourTransport`] can add its original
/// independently supplied reciprocal thermal cross coefficients unchanged.
///
/// Each pair's context is its reference query and retains external coordinates.
/// Its selected claim is pinned there and refreshed at actual face T,p; source
/// boundaries never trigger an implicit switch to another claim. Pair cards
/// describe binary ideal-gas data, not a mixture-composition-dependent law.
/// Composition coordinates or aliases for evolving T,p therefore refuse.
/// Source axes and applicability bounds must carry supported static external
/// meanings; untyped applicability restrictions and frequency/field/rate
/// dependence refuse. Unused context labels remain in the original receipts.
/// All pairs are required, including species currently at exact zero fraction;
/// the original continuum source bounds the bank to 128 species.
/// Source queries allocate; this is not the allocation-free kinetic-data path.
pub struct MaterialBinaryTransport {
    inner: MaterialTransport<BinaryCardLaw>,
}

impl MaterialBinaryTransport {
    /// Bind one original source bank and exactly one card per unordered pair.
    /// Source models must match the material bank in full, not just by label.
    /// Pair ordering and orientation are immaterial. Each pair reference query
    /// must declare absolute `temperature` and `pressure` in coherent SI.
    pub fn new(
        material: ThermochemicalGasModel,
        context: QueryPoint,
        source: IdealGasEquilibriumV1,
        pairs: &[BinaryDiffusionMaterial<'_>],
        bulk_viscosity: DynViscosity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if material.thermodynamics().components() != source.models()
            || !bulk_viscosity.value().is_finite()
            || bulk_viscosity.value() < 0.0
        {
            return Err(DissipationError::Invalid(
                "binary material transport caloric source or bulk viscosity",
            ));
        }
        let diffusion = BinaryCardLaw::new(source, pairs, &mut cancelled)?;
        let inner =
            MaterialTransport::new(material, context, diffusion, bulk_viscosity, &mut cancelled)?;
        Ok(Self { inner })
    }

    /// Original species transport cards and caloric/material identity.
    pub const fn material(&self) -> &ThermochemicalGasModel {
        self.inner.material()
    }

    /// Original external context for the gas's scalar transport properties.
    pub const fn context(&self) -> &QueryPoint {
        self.inner.context()
    }

    /// Explicit bulk viscosity, never inferred from binary diffusion data.
    pub const fn bulk_viscosity(&self) -> DynViscosity {
        self.inner.bulk_viscosity()
    }

    /// Canonical original assignments, still using the shared material type.
    pub fn pair_sources(&self) -> impl ExactSizeIterator<Item = BinaryDiffusionMaterial<'_>> {
        self.inner
            .diffusion
            .pairs
            .iter()
            .map(|pair| BinaryDiffusionMaterial {
                species: pair.species.clone(),
                material: &pair.material,
                context: &pair.context,
                selection: pair.selection.clone(),
            })
    }

    /// Last complete successful pair query in canonical unordered-pair order.
    /// Initially these are the binding reference queries. Failed/cancelled
    /// complete evaluations preserve them and the caller's binary output.
    /// A query receipt is not an accepted evolving gas state.
    pub fn last_resolved(&self) -> &[BinaryDiffusionSample] {
        &self.inner.diffusion.resolved
    }
}

impl TransportLaw for MaterialBinaryTransport {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.inner.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        binary_diffusivities: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        self.inner.diffusion.staged = None;
        let coefficients = self.inner.evaluate(
            temperature,
            pressure,
            mole_fractions,
            binary_diffusivities,
            cancelled,
        )?;
        if let Some(resolved) = self.inner.diffusion.staged.take() {
            self.inner.diffusion.resolved = resolved;
        }
        Ok(coefficients)
    }
}

struct PairCard {
    species: [SpeciesId; 2],
    indices: [usize; 2],
    material: MaterialCard,
    context: QueryPoint,
    selection: MaterialPropertySelection,
}

// This private provider is used only by MaterialTransport, which deliberately
// discards its scalar fields and evaluates the original gas's scalar cards.
// The public complete law can therefore never publish these placeholders.
struct BinaryCardLaw {
    source: IdealGasEquilibriumV1,
    pairs: Vec<PairCard>,
    resolved: Vec<BinaryDiffusionSample>,
    staged: Option<Vec<BinaryDiffusionSample>>,
}

impl BinaryCardLaw {
    fn new(
        source: IdealGasEquilibriumV1,
        bindings: &[BinaryDiffusionMaterial<'_>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        let n = source.models().len();
        if !(1..=crate::multicomponent::MAX_SPECIES).contains(&n)
            || bindings.len() != n * (n - 1) / 2
        {
            return Err(DissipationError::Invalid(
                "one binary material per species pair required",
            ));
        }
        let mut pairs = Vec::with_capacity(bindings.len());
        for binding in bindings {
            poll(cancelled)?;
            let mut indices = [0; 2];
            for (index, id) in indices.iter_mut().zip(&binding.species) {
                *index = source
                    .models()
                    .iter()
                    .position(|m| m.species() == id)
                    .ok_or(DissipationError::Invalid("unknown binary material species"))?;
            }
            indices.sort_unstable();
            if indices[0] == indices[1] {
                return Err(DissipationError::Invalid("self binary material pair"));
            }
            admit_context(binding.context)?;
            pairs.push(PairCard {
                species: indices.map(|i| source.models()[i].species().clone()),
                indices,
                material: binding.material.clone(),
                context: binding.context.clone(),
                selection: binding.selection.clone(),
            });
        }
        pairs.sort_by_key(|pair| pair.indices);
        if pairs.windows(2).any(|p| p[0].indices == p[1].indices) {
            return Err(DissipationError::Invalid("duplicate binary material pair"));
        }
        let requirement = ScalarPropertyRequirement::try_new(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(|e| material_input("binary diffusion requirement", e))?;
        let mut resolved = Vec::with_capacity(pairs.len());
        for pair in &pairs {
            poll(cancelled)?;
            let material = resolve_material_state_point(
                &pair.material,
                &pair.context,
                std::slice::from_ref(&requirement),
                pair.selection.clone(),
            )
            .map_err(|e| material_input("binary diffusion source", e))?;
            resolved.push(sample(pair, material)?);
        }
        poll(cancelled)?;
        Ok(Self {
            source,
            pairs,
            resolved,
            staged: None,
        })
    }
}

impl TransportLaw for BinaryCardLaw {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        &self.source
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
        let n = self.source.models().len();
        admit_state(n, temperature, pressure, fractions)?;
        if output.len() != n * n {
            return Err(DissipationError::Invalid("binary material output axis"));
        }
        let mut resolved = Vec::with_capacity(self.pairs.len());
        for (pair, pinned) in self.pairs.iter().zip(&self.resolved) {
            poll(cancelled)?;
            let point = pair
                .context
                .clone()
                .with("temperature", temperature.value())
                .and_then(|q| q.with("pressure", pressure.value()))
                .map_err(|e| material_input("binary diffusion coordinates", e))?;
            let material = pinned
                .material
                .at_query_point(&pair.material, &point)
                .map_err(|e| material_input("binary diffusion refresh", e))?;
            resolved.push(sample(pair, material)?);
        }
        poll(cancelled)?;
        output.fill(0.0);
        for (pair, sample) in self.pairs.iter().zip(&resolved) {
            let [i, j] = pair.indices;
            output[i * n + j] = sample.diffusivity_m2_s;
            output[j * n + i] = sample.diffusivity_m2_s;
        }
        validate_binary(output, n)?;
        self.staged = Some(resolved);
        Ok(TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.0),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.0),
        })
    }
}

fn sample(
    pair: &PairCard,
    material: fs_material::state_point::ResolvedMaterialStatePoint,
) -> Result<BinaryDiffusionSample, DissipationError> {
    let property =
        material
            .property(BINARY_DIFFUSIVITY_PROPERTY)
            .ok_or(DissipationError::Invalid(
                "missing resolved binary diffusivity",
            ))?;
    if !property.has_only_fixed_external_conditions(Some("temperature"), false) {
        return Err(DissipationError::Invalid(
            "binary diffusion requires static local source support",
        ));
    }
    let diffusivity = property.value_si();
    if !diffusivity.is_finite() || diffusivity <= 0.0 {
        return Err(DissipationError::Invalid(
            "binary diffusivity must be positive finite m2/s",
        ));
    }
    Ok(BinaryDiffusionSample {
        species: pair.species.clone(),
        diffusivity_m2_s: diffusivity,
        material,
    })
}

fn admit_context(context: &QueryPoint) -> Result<(), DissipationError> {
    for (axis, kind) in [
        ("temperature", QuantityKind::AbsoluteTemperature),
        ("pressure", QuantityKind::Pressure),
    ] {
        if context.axis_quantities().get(axis) != Some(&typed(kind))
            || !context
                .axes()
                .get(axis)
                .is_some_and(|v| v.is_finite() && *v > 0.0)
        {
            return Err(DissipationError::Invalid(
                "binary source needs typed reference T,p",
            ));
        }
    }
    if context.axis_quantities().iter().any(|(axis, quantity)| {
        (quantity.dims() == Temperature::DIMS && axis != "temperature")
            || (quantity.dims() == Pressure::DIMS && axis != "pressure")
            || quantity
                .semantic_type()
                .is_some_and(|q| matches!(q.kind(), QuantityKind::Composition(_)))
    }) {
        return Err(DissipationError::Invalid(
            "unowned evolving binary diffusion coordinate",
        ));
    }
    Ok(())
}
