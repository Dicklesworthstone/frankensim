//! Reciprocal material convection between two original finite thermal cells.
//! Existing thermal links own area, both energy debits, CFL and integration.

use super::super::{CellConvectionEvaluation, evaluate_convection_state};
use super::*;

/// An existing thermal contact whose total conductance comes from convection.
/// `fluid_cell` names one endpoint; the other is the wall. The original link's
/// area and per-side coverage remain authoritative. Zero authored conductance
/// is required so a separate Fourier path cannot silently disappear.
#[derive(Clone, Copy, Debug)]
pub struct SourceConvectiveLink {
    pub link: u64,
    pub fluid_cell: CellId,
    pub characteristic_length: Length,
    /// Reservoir means the actual finite FLUID cell's bulk temperature.
    pub property_temperature: ConvectionPropertyTemperature,
    pub law: SourceCellConvectionLaw,
}

/// Complete material source and interface evaluation at one proposed endpoint.
#[derive(Debug)]
pub struct ConvectiveLinkSample<F> {
    pub link: u64,
    pub wall_cell: CellId,
    pub fluid_cell: CellId,
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    pub property_temperature_k: f64,
    pub area_m2: f64,
    pub conductance_w_k: f64,
    /// Signed wall-to-fluid power. The native link applies its exact opposite
    /// to the wall; this is internal exchange, never external source power.
    pub heat_into_fluid_w: f64,
    pub material: F,
    pub thermal: ThermalFluidTransportState,
    /// None only at exact natural-convection thermal equilibrium.
    pub evaluation: Option<CellConvectionEvaluation>,
}

/// Step-local native link conductances and their retained source owners.
#[derive(Debug)]
pub struct ResolvedLinkConvection<F> {
    conductances_w_k: BTreeMap<u64, f64>,
    samples: Vec<ConvectiveLinkSample<F>>,
}
impl<F> ResolvedLinkConvection<F> {
    pub fn conductances_w_k(&self) -> &BTreeMap<u64, f64> {
        &self.conductances_w_k
    }
    pub fn samples(&self) -> &[ConvectiveLinkSample<F>] {
        &self.samples
    }

    /// Insert into an existing implicit trial's conductance map. Each link
    /// must have one owner; even an identical competing override refuses.
    /// Heat belongs in this map, not a second pair of generated source rates.
    pub fn apply_to(
        &self,
        destination: &mut BTreeMap<u64, f64>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ConvectionSourceError> {
        check(&mut cancelled)?;
        for link in self.conductances_w_k.keys() {
            check(&mut cancelled)?;
            if destination.contains_key(link) {
                return Err(ConvectionSourceError::Invalid(
                    "convective link already has a conductance owner",
                ));
            }
        }
        check(&mut cancelled)?;
        destination.extend(self.conductances_w_k.iter().map(|(&id, &g)| (id, g)));
        Ok(())
    }
}

/// Resolve reciprocal convection using both ACTUAL caloric temperatures.
///
/// The source resolver sees the original thermal state, fluid identity and
/// declared bulk/film property temperature. Bind it to the same fluid material
/// and fixed conditions as the native caloric law; a generic transport trait
/// cannot prove this identity or invent pressure/composition evolution. Retain
/// prepared material laws across calls so nonlinear trials do not select new
/// sources. Every returned material owner stays with its heat-transfer result.
///
/// This returns total conductances for existing links. Use them directly in
/// `BackwardEulerStep` coefficients or the native observer's conductance hook:
/// the original thermal owner then debits and credits the SAME heat and applies
/// its combined outgoing-conductance timestep bound. It neither changes cell
/// mass nor creates reservoir heat, a second fluid temperature or a solver.
///
/// Links must be explicit non-overlapping interface patches in the admitted
/// geometry. Their original coverage excludes ambient exchange on those faces.
/// Several walls may exchange with one finite fluid cell through distinct links.
/// Uniform-cell validity, correlation support, prescribed relative flow and
/// source time/space refinement remain the caller's physical model obligations.
/// No momentum, mass transport, phase motion or pressure work is added here.
pub fn evaluate_source_link_convection<F: ThermalFluidTransport, E>(
    thermal: &DepositionSimulation<'_>,
    bindings: &[SourceConvectiveLink],
    maximum_links: usize,
    mut resolve: impl FnMut(
        &SourceConvectiveLink,
        &DepositionSimulation<'_>,
        Temperature,
    ) -> Result<F, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ResolvedLinkConvection<F>, SourceCellConvectionError<E>> {
    check(&mut cancelled)?;
    if bindings.len() > maximum_links || bindings.len() > thermal.links().len() {
        return Err(ConvectionSourceError::Invalid(
            "convective links exceed their explicit work bound",
        )
        .into());
    }
    let mut ids = BTreeSet::new();
    let mut ordered = Vec::with_capacity(bindings.len());
    for binding in bindings {
        check(&mut cancelled)?;
        let link = thermal
            .links()
            .get(&binding.link)
            .ok_or(ConvectionSourceError::Invalid("inactive convective link"))?;
        if !ids.insert(binding.link) {
            return Err(ConvectionSourceError::Invalid("duplicate convective link").into());
        }
        if binding.fluid_cell != link.a && binding.fluid_cell != link.b {
            return Err(ConvectionSourceError::Invalid(
                "convective fluid must be an original link endpoint",
            )
            .into());
        }
        if link.conductance_w_k != 0.0 {
            return Err(ConvectionSourceError::Invalid(
                "convective links require zero authored conductance",
            )
            .into());
        }
        if !binding.characteristic_length.value().is_finite()
            || binding.characteristic_length.value() <= 0.0
        {
            return Err(ConvectionSourceError::Invalid(
                "convective links require positive finite characteristic length",
            )
            .into());
        }
        ordered.push(binding);
    }
    ordered.sort_by_key(|binding| binding.link);
    let mut conductances_w_k = BTreeMap::new();
    let mut samples = Vec::with_capacity(ordered.len());
    for binding in ordered {
        check(&mut cancelled)?;
        let link = &thermal.links()[&binding.link];
        let fluid_cell = binding.fluid_cell;
        let wall_cell = if fluid_cell == link.a { link.b } else { link.a };
        let wall_temperature_k = thermal.temperature_k(wall_cell)?;
        let fluid_temperature_k = thermal.temperature_k(fluid_cell)?;
        let property_temperature_k = match binding.property_temperature {
            ConvectionPropertyTemperature::Reservoir => fluid_temperature_k,
            ConvectionPropertyTemperature::FilmMean => {
                0.5 * wall_temperature_k + 0.5 * fluid_temperature_k
            }
        };
        let material = resolve(binding, thermal, Temperature::new(property_temperature_k))
            .map_err(|source| SourceCellConvectionError::Source {
                cell: fluid_cell,
                property_temperature_k,
                source,
            })?;
        let snapshot = material
            .thermal_transport_state()
            .map_err(|error| ConvectionSourceError::Material(error.into()))?;
        let transport = material
            .transport_state()
            .map_err(|error| ConvectionSourceError::Material(error.into()))?;
        if snapshot.temperature_k() != property_temperature_k || snapshot.transport() != transport {
            return Err(ConvectionSourceError::Invalid(
                "fluid source does not match the requested temperature and transport state",
            )
            .into());
        }
        let law = match binding.law {
            SourceCellConvectionLaw::Forced {
                correlation,
                speed_m_s,
                settings,
            } => CellConvectionLaw::Forced {
                correlation,
                speed_m_s,
                settings,
            },
            SourceCellConvectionLaw::VerticalPlate { gravity_m_s2 } => {
                CellConvectionLaw::VerticalPlate {
                    gravity_m_s2,
                    expansion_per_k: snapshot.volumetric_expansion_per_k().ok_or(
                        ConvectionSourceError::Invalid(
                            "natural convection requires sourced isobaric volumetric expansion",
                        ),
                    )?,
                }
            }
        };
        let evaluation = evaluate_convection_state(
            &transport,
            wall_temperature_k,
            Temperature::new(fluid_temperature_k),
            binding.characteristic_length,
            law,
        )
        .map_err(ConvectionSourceError::Material)?;
        let conductance_w_k = match &evaluation {
            Some(sample) => finite(link.area_m2 * sample.coefficient.value.value())?,
            None => 0.0,
        };
        let difference = finite(wall_temperature_k - fluid_temperature_k)?;
        let heat_into_fluid_w = finite(conductance_w_k * difference)?;
        if (evaluation.is_some() && link.area_m2 > 0.0 && conductance_w_k == 0.0)
            || (conductance_w_k != 0.0 && difference != 0.0 && heat_into_fluid_w == 0.0)
        {
            return Err(ConvectionSourceError::Invalid(
                "convective link transfer underflowed to zero",
            )
            .into());
        }
        conductances_w_k.insert(binding.link, conductance_w_k);
        samples.push(ConvectiveLinkSample {
            link: binding.link,
            wall_cell,
            fluid_cell,
            wall_temperature_k,
            fluid_temperature_k,
            property_temperature_k,
            area_m2: link.area_m2,
            conductance_w_k,
            heat_into_fluid_w,
            material,
            thermal: snapshot,
            evaluation,
        });
    }
    check(&mut cancelled)?;
    Ok(ResolvedLinkConvection {
        conductances_w_k,
        samples,
    })
}
