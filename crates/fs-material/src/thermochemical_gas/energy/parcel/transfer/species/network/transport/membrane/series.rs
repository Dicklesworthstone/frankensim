//! Heterogeneous dense layers and explicit interlayer resistance in series.
//!
//! At fixed temperature, 1/G_i = sum_j(L_j/P_ji) + sum_j(1/H_ji),
//! with J_i = G_i (p_i,left - p_i,right). Layer/interface conductances are
//! scaled by their minimum before summing, avoiding reciprocals of tiny
//! permeances. The same J_i crosses every element and is spent only once.
//! Lillepärg et al., Membranes 9(2):22 (2019), Eqs. 14–15:
//! <https://doi.org/10.3390/membranes9020022>.
//!
//! The interior partial pressures are gas-equivalent Henry potentials, not
//! actual gas pockets, polymer stress or an additional caloric inventory.
//! All layers share one admitted mean temperature. Explicit interface cards
//! describe reciprocal linear transfer for the ordered adjacent materials;
//! absent interfaces mean local equilibrium. No porous support, bypass defect,
//! nonlinear sorption or unresolved thickness effect is inferred.

use super::*;
use crate::state_point::{ResolvedInterfaceStatePoint, ResolvedScalarProperty,
    resolve_interface_state_point};
use fs_matdb::InterfaceSystemCard;

pub const MAX_MEMBRANE_ELEMENTS: usize = 64;
pub const MOLAR_INTERFACE_PERMEANCE_PROPERTY: &str = "gas_molar_interface_permeance";

/// Per-species interlayer transfer coefficient on the actual ordered system.
#[derive(Clone, Debug)]
pub struct MembraneInterfaceMaterial<'a> {
    pub species: SpeciesId,
    pub interface: &'a InterfaceSystemCard,
    /// Fixed process/mixture/humidity context. An optional canonical pressure
    /// selects the original pressure-independent claim before interfaces are
    /// solved; its actual support is checked again at both local boundaries.
    pub context: &'a QueryPoint,
    pub selection: MaterialPropertySelection,
}

/// One ordered piece of the physical laminate. Each element explicitly binds
/// EVERY catalog species, including an authored zero for an impermeable one.
#[derive(Clone, Copy, Debug)]
pub enum MembraneSeriesElement<'a> {
    Layer {
        thickness_m: f64,
        materials: &'a [MembranePermeabilityMaterial<'a>],
    },
    /// An interlayer interface must have a dense layer immediately on each
    /// side. Its card must identify those original ordered material states.
    Interface { materials: &'a [MembraneInterfaceMaterial<'a>] },
}

#[derive(Clone, Copy, Debug)]
pub enum MembranePath<'a> {
    Film {
        geometry: MembraneGeometry,
        materials: &'a [MembranePermeabilityMaterial<'a>],
    },
    /// A single shared area; a parallel leak is a separate network path.
    Series { area_m2: f64, elements: &'a [MembraneSeriesElement<'a>] },
}

/// Original physical coefficients; an effective laminate permeability is
/// deliberately not substituted for these source-bearing component laws.
#[derive(Clone, Debug, PartialEq)]
pub enum MembraneElementLaw {
    Layer {
        thickness_m: f64,
        permeability: MolarPermeability,
        material: ResolvedMaterialStatePoint,
    },
    Interface { material: ResolvedInterfaceStatePoint },
}

#[derive(Clone, Debug, PartialEq)]
pub struct MembraneElementSample {
    pub law: MembraneElementLaw,
    /// Per-area pressure-normalized molar flux coefficient of THIS element.
    pub permeance: MolarPermeance,
    pub partial_pressures_pa: [f64; 2],
}

#[derive(Clone, Copy)]
enum Binding<'a> {
    Layer(f64, &'a MembranePermeabilityMaterial<'a>),
    Interface(&'a MembraneInterfaceMaterial<'a>),
}

fn invalid(message: &'static str) -> GasEnergyError { GasEnergyError::Invalid(message) }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasEnergyError> {
    if cancelled() { Err(GasEnergyError::Cancelled) } else { Ok(()) }
}

fn admit_henry(property: &ResolvedScalarProperty) -> Result<(), GasEnergyError> {
    if !property.is_constant_at_fixed_temperature()
        || property.source_axes().iter().any(|(name, _)| name != "temperature")
        || !property.is_static_at_fixed_external_conditions()
        || property.answer().evidence.model.validity.bounds().keys().any(|name| {
            name == "T" || property.answer().receipt.axis_quantities.get(name).is_some_and(|q|
                (q.dims() == fs_qty::Temperature::DIMS && name != "temperature")
                || (q.dims() == fs_qty::Pressure::DIMS && name != "pressure"))
        }) {
        return Err(invalid("linear solution-diffusion needs a constant-at-fixed-temperature source with canonical T/p axes"));
    }
    Ok(())
}

impl Binding<'_> {
    fn sample(self, temperature: f64, default_pressure: f64)
        -> Result<MembraneElementSample, GasEnergyError>
    {
        let context = match self { Self::Layer(_, b) => b.context, Self::Interface(b) => b.context };
        let pressure = context.axes().get("pressure").copied().unwrap_or(default_pressure);
        if !positive(pressure) { return Err(invalid("membrane source selection pressure must be positive")); }
        let point = updated_point(context, temperature, pressure)?;
        let (law, permeance) = match self {
            Self::Layer(thickness_m, binding) => {
                let requirement = ScalarPropertyRequirement::try_new(MOLAR_PERMEABILITY_PROPERTY,
                    MOLAR_PERMEABILITY_DIMS, ScalarAdmissibility::NonNegative).map_err(source_error)?;
                let material = resolve_material_state_point(binding.material, &point,
                    &[requirement], binding.selection.clone()).map_err(source_error)?;
                let property = material.property(MOLAR_PERMEABILITY_PROPERTY)
                    .ok_or(invalid("resolved membrane permeability is missing"))?;
                admit_henry(property)?;
                let permeability = MolarPermeability::try_new(property.value_si())?;
                let permeance = permeability.permeance(thickness_m)?;
                (MembraneElementLaw::Layer { thickness_m, permeability, material }, permeance)
            }
            Self::Interface(binding) => {
                let requirement = ScalarPropertyRequirement::try_new(MOLAR_INTERFACE_PERMEANCE_PROPERTY,
                    MOLAR_PERMEANCE_DIMS, ScalarAdmissibility::NonNegative).map_err(source_error)?;
                let material = resolve_interface_state_point(binding.interface, &point,
                    &[requirement], binding.selection.clone()).map_err(source_error)?;
                let property = material.property(MOLAR_INTERFACE_PERMEANCE_PROPERTY)
                    .ok_or(invalid("resolved interlayer molar permeance is missing"))?;
                admit_henry(property)?;
                let permeance = MolarPermeance::try_new(property.value_si())?;
                (MembraneElementLaw::Interface { material }, permeance)
            }
        };
        Ok(MembraneElementSample { law, permeance, partial_pressures_pa: [0.0; 2] })
    }

    /// Preserve exact claims while checking ACTUAL interface conditions. Source
    /// selection pressure is not permission to extend a material's domain.
    fn admit(self, sample: &mut MembraneElementSample, temperatures: [f64; 2],
        pressures: [f64; 2], face_temperature: f64, cancelled: &mut impl FnMut() -> bool)
        -> Result<(), GasEnergyError>
    {
        let pressure = 0.5 * pressures[0] + 0.5 * pressures[1];
        match (self, &mut sample.law) {
            (Self::Layer(_, binding), MembraneElementLaw::Layer { material, .. }) => {
                for i in 0..2 {
                    poll(cancelled)?;
                    material.at_query_point(binding.material,
                        &updated_point(binding.context, temperatures[i], pressures[i])?).map_err(source_error)?;
                }
                *material = material.at_query_point(binding.material,
                    &updated_point(binding.context, face_temperature, pressure)?).map_err(source_error)?;
            }
            (Self::Interface(binding), MembraneElementLaw::Interface { material }) => {
                for i in 0..2 {
                    poll(cancelled)?;
                    material.at_query_point(binding.interface,
                        &updated_point(binding.context, temperatures[i], pressures[i])?).map_err(source_error)?;
                }
                *material = material.at_query_point(binding.interface,
                    &updated_point(binding.context, face_temperature, pressure)?).map_err(source_error)?;
            }
            _ => return Err(invalid("membrane source/element kinds differ")),
        }
        Ok(())
    }
}

fn series_permeance(samples: &mut [MembraneElementSample], pressures: [f64; 2])
    -> Result<MolarPermeance, GasEnergyError>
{
    let mut seal = None;
    let mut scale = f64::INFINITY;
    for (i, sample) in samples.iter().enumerate() {
        if sample.permeance.0 == 0.0 {
            if seal.replace(i).is_some() {
                return Err(invalid("multiple impermeable elements leave an interior species potential undetermined; a stored-state model is required"));
            }
        } else { scale = scale.min(sample.permeance.0); }
    }
    if let Some(seal) = seal {
        for (i, sample) in samples.iter_mut().enumerate() {
            sample.partial_pressures_pa = if i < seal { [pressures[0]; 2] }
                else if i > seal { [pressures[1]; 2] } else { pressures };
        }
        return Ok(MolarPermeance(0.0));
    }
    let mut total = Sum::default();
    for sample in samples.iter() {
        let weight = scale / sample.permeance.0;
        if !positive(weight) { return Err(invalid("series membrane resistance ratio is unrepresentable")); }
        total.add(weight)?;
    }
    let permeance = scale / total.value;
    if !positive(permeance) { return Err(invalid("series membrane permeance is unrepresentable")); }
    let mut prefix = Sum::default();
    let mut left = pressures[0];
    let last = samples.len() - 1;
    for (i, sample) in samples.iter_mut().enumerate() {
        prefix.add(scale / sample.permeance.0)?;
        let fraction = (prefix.value / total.value).clamp(0.0, 1.0);
        let right = if i == last { pressures[1] }
            else { (1.0-fraction)*pressures[0] + fraction*pressures[1] };
        if !right.is_finite() || right < 0.0 {
            return Err(invalid("membrane interface partial pressure is unrepresentable"));
        }
        sample.partial_pressures_pa = [left, right];
        left = right;
    }
    Ok(MolarPermeance(permeance))
}

impl ClosedGasParcel {
    /// Resolve one film or laminate on the unchanged endpoint inventories.
    ///
    /// A layer's canonical `pressure` context, when present, only anchors its
    /// original pressure-independent claim. After the flux solve, the SAME
    /// claim must support both actual local equivalent-gas total pressures.
    /// Reservoir-only checks are insufficient: selective counterpermeation can
    /// put an interior total pressure outside the reservoir range. All fixed
    /// process/mixture/humidity context remains caller-qualified.
    ///
    /// One exact seal per species is supported. Multiple seals isolate an
    /// undetermined interior species potential and refuse this no-storage rung.
    pub fn membrane_path_flux_to(&self, right: &Self, path: MembranePath<'_>,
        maximum_temperature_difference_k: f64, mut cancelled: impl FnMut() -> bool)
        -> Result<MembraneFlux, GasEnergyError>
    {
        match path {
            MembranePath::Film { geometry, materials } => sample_path(self, right,
                geometry.area_m2, &[MembraneSeriesElement::Layer {
                    thickness_m: geometry.thickness_m, materials }],
                maximum_temperature_difference_k, &mut cancelled),
            MembranePath::Series { area_m2, elements } => sample_path(self, right,
                area_m2, elements, maximum_temperature_difference_k, &mut cancelled),
        }
    }
}

fn sample_path(left: &ClosedGasParcel, right: &ClosedGasParcel, area_m2: f64,
    elements: &[MembraneSeriesElement<'_>], maximum_temperature_difference_k: f64,
    cancelled: &mut impl FnMut() -> bool) -> Result<MembraneFlux, GasEnergyError>
{
    poll(cancelled)?;
    if !area_m2.is_finite() || area_m2 < 0.0 || elements.is_empty()
        || elements.len() > MAX_MEMBRANE_ELEMENTS
        || !maximum_temperature_difference_k.is_finite() || maximum_temperature_difference_k < 0.0 {
        return Err(invalid("membrane needs bounded nonempty elements, nonnegative area and temperature allowance"));
    }
    for (i, element) in elements.iter().enumerate() {
        poll(cancelled)?;
        match element {
            MembraneSeriesElement::Layer { thickness_m, materials }
                if positive(*thickness_m) && materials.len() <= MAX_MEMBRANE_SPECIES => {},
            MembraneSeriesElement::Interface { materials }
                if materials.len() <= MAX_MEMBRANE_SPECIES && i > 0 && i+1 < elements.len()
                    && matches!(elements[i-1], MembraneSeriesElement::Layer { .. })
                    && matches!(elements[i+1], MembraneSeriesElement::Layer { .. }) => {},
            _ => return Err(invalid("membrane layers need positive thickness and bounded species; interfaces require adjacent dense layers")),
        }
    }
    let (a, b) = (left.state().state(), right.state().state());
    let mut flux = MembraneFlux { face_temperature_k: 0.5*a.temperature+0.5*b.temperature,
        face_pressure_pa: 0.5*a.pressure+0.5*b.pressure, species: Vec::new(),
        interface_total_pressures_pa: Vec::new() };
    if area_m2 == 0.0 { poll(cancelled)?; return Ok(flux); }
    if (a.temperature-b.temperature).abs() > maximum_temperature_difference_k {
        return Err(invalid("gas states exceed the declared isothermal membrane support"));
    }
    let joined = catalog(left, right)?;
    let components = joined.thermodynamics().components();
    if components.is_empty() || components.len() > MAX_MEMBRANE_SPECIES {
        return Err(invalid("membrane requires a bounded nonempty species catalog"));
    }
    let mut ordered = Vec::with_capacity(elements.len());
    for element in elements {
        let mut map = BTreeMap::new();
        match element {
            MembraneSeriesElement::Layer { thickness_m, materials } => {
                for binding in *materials {
                    poll(cancelled)?;
                    if map.insert(&binding.species, Binding::Layer(*thickness_m, binding)).is_some() {
                        return Err(invalid("duplicate membrane layer species permeability"));
                    }
                }
            }
            MembraneSeriesElement::Interface { materials } => {
                for binding in *materials {
                    poll(cancelled)?;
                    if map.insert(&binding.species, Binding::Interface(binding)).is_some() {
                        return Err(invalid("duplicate membrane interface species permeance"));
                    }
                }
            }
        }
        if map.len() != components.len() || components.iter().any(|s| !map.contains_key(s.species())) {
            return Err(invalid("each membrane element needs one explicit coefficient per catalog species"));
        }
        ordered.push(map);
    }
    for species in components {
        let id = species.species();
        let mut samples = Vec::with_capacity(elements.len());
        for (i, bindings) in ordered.iter().enumerate() {
            poll(cancelled)?;
            let binding = bindings[id];
            if let Binding::Interface(interface) = binding {
                let (Binding::Layer(_, before), Binding::Layer(_, after)) =
                    (ordered[i-1][id], ordered[i+1][id]) else {
                        return Err(invalid("membrane interface lacks its adjacent material bindings"));
                    };
                if interface.interface.surface_a().material != *before.material.id()
                    || interface.interface.surface_b().material != *after.material.id() {
                    return Err(invalid("membrane interface source does not identify its ordered adjacent materials"));
                }
            }
            samples.push(binding.sample(flux.face_temperature_k, flux.face_pressure_pa)?);
        }
        let pressures = [partial_pressure(left, id)?, partial_pressure(right, id)?];
        let permeance = series_permeance(&mut samples, pressures)?;
        let difference = pressures[0]-pressures[1];
        let molar_flux = permeance.0*difference;
        let rate = (molar_flux*species.molar_mass().value())*area_m2;
        if !molar_flux.is_finite() || !rate.is_finite()
            || (permeance.0 != 0.0 && difference != 0.0 && rate == 0.0) {
            return Err(invalid("membrane species flux is unrepresentable"));
        }
        flux.species.push(MembraneSpeciesFlux { species: id.clone(), partial_pressures_pa: pressures,
            permeance, molar_flux_mol_m2_s: molar_flux, mass_rate_left_to_right_kg_s: rate,
            elements: samples });
    }
    flux.species.sort_by(|a, b| a.species.cmp(&b.species));
    for i in 0..=elements.len() {
        poll(cancelled)?;
        let mut total = Sum::default();
        for row in &flux.species {
            total.add(if i == elements.len() { row.elements[i-1].partial_pressures_pa[1] }
                else { row.elements[i].partial_pressures_pa[0] })?;
        }
        flux.interface_total_pressures_pa.push(total.value);
    }
    // Reservoir totals come directly from their coherent EOS, avoiding a
    // one-ULP composition sum accidentally extending/narrowing source support.
    flux.interface_total_pressures_pa[0] = a.pressure;
    flux.interface_total_pressures_pa[elements.len()] = b.pressure;
    for row in &mut flux.species {
        for (i, sample) in row.elements.iter_mut().enumerate() {
            poll(cancelled)?;
            let temperatures = [if i == 0 { a.temperature } else { flux.face_temperature_k },
                if i+1 == elements.len() { b.temperature } else { flux.face_temperature_k }];
            let pressures = [flux.interface_total_pressures_pa[i], flux.interface_total_pressures_pa[i+1]];
            ordered[i][&row.species].admit(sample, temperatures, pressures, flux.face_temperature_k, cancelled)?;
        }
    }
    poll(cancelled)?;
    Ok(flux)
}
