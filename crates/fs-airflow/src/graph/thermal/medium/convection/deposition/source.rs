//! Re-evaluate the original fluid source at each finite-body thermal trial.
//!
//! A correlation's property temperature and the physical reservoir temperature
//! are different quantities. In particular, film properties use the mean wall
//! and reservoir temperature while heat and Rayleigh forcing retain the full
//! wall/reservoir difference. No additional thermal store or solver is created.

/// Material fluid cooling in the complete shared heat/circuit solve.
pub mod heat;
/// Source-backed convection on the existing explicit thermal timeline.
pub mod observer;
/// Reciprocal convection on the original finite-cell thermal links.
pub mod links;
/// Quasi-steady enclosed gas layers between two finite isothermal walls.
pub mod enclosure;
/// Source-owned axial storage, conduction and finite-fluid cooling of fins.
pub mod fin;

use super::*;
use fs_material::fluid_transport::{ThermalFluidTransport, ThermalFluidTransportState};
use fs_material::liquid::ResolvedLiquidState;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Authored property-evaluation convention required by the selected correlation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConvectionPropertyTemperature {
    Reservoir,
    FilmMean,
}

/// Expansion is supplied by the same resolved material as the transport tuple.
#[derive(Clone, Copy, Debug)]
pub enum SourceCellConvectionLaw {
    /// Original forced relation. Sieder-Tate additionally resolves this same
    /// fluid at the wall for viscosity; it requires Reservoir bulk properties
    /// and an explicit circular tube heated length/diameter in `settings`.
    /// For SiederTateLaminar, the reference temperature is the authored
    /// arithmetic mean bulk temperature over the WHOLE heated tube. It is
    /// neither a local downstream fluid state nor a log-mean temperature;
    /// the caller retains that mean convention and uniform wall temperature.
    /// GnielinskiDeveloping also requires bulk properties and physical heated
    /// L/D, with its explicitly constant-property K=1 approximation.
    Forced {
        correlation: CorrelationId,
        speed_m_s: f64,
        settings: ForcedConvectionSettings,
    },
    /// External flat-plate strip using the actual mean-film source state.
    /// Characteristic length is the downstream distance from the physical
    /// leading edge; `start` is measured from that same edge. Only the original
    /// laminar and turbulent mean cards are admitted. Rebind after reversing
    /// flow; local source updates retain the held-property approximation.
    ForcedPlateSection {
        correlation: CorrelationId,
        speed_m_s: f64,
        start: Length,
    },
    /// Developing turbulent circular-tube interval. Characteristic length is
    /// diameter; start and end retain the physical heated entrance. Reservoir
    /// bulk properties are mandatory for the explicit K=1 model. Re-evaluating
    /// a cell does not restart its inlet or introduce a wall-property correction.
    ForcedTubeSection {
        speed_m_s: f64,
        start: Length,
        end: Length,
    },
    /// Gnielinski's liquid wall-Pr correction on the same turbulent interval.
    /// Primary properties must be actual bulk properties. The SAME resolver
    /// is queried at the wall temperature and must retain original liquid
    /// owners with the same card, selected claims and non-temperature context.
    /// Numeric-only and gas sources refuse rather than imply liquid identity.
    /// The empirical correction and finite-section model retain their original
    /// validity limits; no single-phase or experimental-error bound is inferred.
    ForcedTubeSectionLiquidWall {
        speed_m_s: f64,
        start: Length,
        end: Length,
    },
    /// Hausen laminar circular-tube interval from one physical heated entrance.
    /// Characteristic length is diameter. Resolve the SAME finite fluid at its
    /// actual bulk temperature each trial; no surface correction is appended.
    ForcedLaminarTubeSection {
        speed_m_s: f64,
        start: Length,
        end: Length,
    },
    /// Developing laminar rectangular channel with four isothermal walls.
    /// Characteristic length is hydraulic diameter; active area must cover the
    /// whole perimeter of this interval. The original table still admits only
    /// its sourced aspect ratio and Pr, using actual bulk transport properties.
    ForcedRectangularDuctSection {
        speed_m_s: f64,
        aspect_ratio: f64,
        start: Length,
        end: Length,
    },
    /// Isolated isothermal sphere with Whitaker's original bulk/surface
    /// property convention. Select `Reservoir` for the primary state; the
    /// SAME resolver is additionally queried at the actual wall temperature.
    /// Characteristic length is diameter. No film substitution is admitted.
    ForcedSphere {
        speed_m_s: f64,
    },
    VerticalPlate {
        gravity_m_s2: f64,
    },
    /// Geometry selects the shared natural-convection card. Horizontal plate,
    /// cylinder and sphere cards require `FilmMean` properties; Rayleigh and
    /// heat use the full wall/bulk temperature difference. Expansion comes from the
    /// same original material source, including its sign.
    Natural {
        geometry: NaturalConvectionGeometry,
        gravity_m_s2: f64,
    },
    /// Laminar aiding vertical plate; expansion and transport share the same
    /// source. This adapter declares a mean-film frozen-property approximation.
    AidingVerticalPlate {
        settings: AidingVerticalPlateSettings,
    },
}

impl SourceCellConvectionLaw {
    /// These original relations require a second query of the SAME fluid at
    /// its actual wall temperature. Only that state's viscosity is consumed.
    pub(crate) fn needs_surface_viscosity(self) -> bool {
        matches!(self, Self::ForcedSphere { .. }
            | Self::Forced { correlation: CorrelationId::SiederTateTurbulent
                | CorrelationId::SiederTateLaminar, .. })
    }

    pub(crate) fn needs_surface_transport(self) -> bool {
        self.needs_surface_viscosity()
            || matches!(self, Self::ForcedTubeSectionLiquidWall { .. })
    }

    // Admit prescribed strip geometry before invoking an external resolver.
    // The same native law also checks it when consuming a borrowed snapshot.
    pub(crate) fn admit_geometry(self, length: Length) -> Result<(), ConvectionSourceError> {
        if let Self::ForcedPlateSection { correlation, speed_m_s, start } = self {
            CellConvectionLaw::ForcedPlateSection { correlation, speed_m_s, start }
                .admit_area(0.0, length)?;
        }
        if let Self::ForcedTubeSection { speed_m_s, start, end }
            | Self::ForcedTubeSectionLiquidWall { speed_m_s, start, end } = self
        {
            CellConvectionLaw::ForcedTubeSection { speed_m_s, start, end }
                .admit_area(0.0, length)?;
        }
        if let Self::ForcedLaminarTubeSection { speed_m_s, start, end } = self {
            CellConvectionLaw::ForcedLaminarTubeSection { speed_m_s, start, end }
                .admit_area(0.0, length)?;
        }
        if let Self::ForcedRectangularDuctSection { speed_m_s, aspect_ratio, start, end } = self {
            CellConvectionLaw::ForcedRectangularDuctSection { speed_m_s, aspect_ratio, start, end }
                .admit_area(0.0, length)?;
        }
        if let Self::Forced {
            correlation: CorrelationId::SiederTateTurbulent | CorrelationId::SiederTateLaminar
                | CorrelationId::GnielinskiDeveloping,
            settings, ..
        } = self {
            if settings.aspect_ratio.is_some()
                || settings.length_over_hydraulic_diameter
                    .is_none_or(|ratio| !ratio.is_finite() || ratio <= 0.0)
            {
                return Err(ConvectionSourceError::Invalid(
                    "tube convection requires a circular tube and explicit positive heated length/diameter",
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn admit_property_temperature(
        self,
        convention: ConvectionPropertyTemperature,
    ) -> Result<(), ConvectionSourceError> {
        if self.needs_surface_viscosity()
            && convention != ConvectionPropertyTemperature::Reservoir
        {
            return Err(ConvectionSourceError::Invalid(
                "viscosity-corrected convection requires bulk properties plus separately resolved wall viscosity",
            ));
        }
        if matches!(self, Self::Forced { correlation: CorrelationId::GnielinskiDeveloping, .. }
            | Self::ForcedTubeSection { .. }
            | Self::ForcedTubeSectionLiquidWall { .. }
            | Self::ForcedLaminarTubeSection { .. }
            | Self::ForcedRectangularDuctSection { .. })
            && convention != ConvectionPropertyTemperature::Reservoir
        {
            return Err(ConvectionSourceError::Invalid(
                "developing tube and rectangular-duct sections require actual bulk properties with their original held-property approximation",
            ));
        }
        if matches!(
            self,
            Self::Natural {
                geometry: NaturalConvectionGeometry::HorizontalCylinder
                    | NaturalConvectionGeometry::HorizontalPlate(_)
                    | NaturalConvectionGeometry::Sphere,
                ..
            }
                | Self::AidingVerticalPlate { .. }
                | Self::ForcedPlateSection { .. }
        ) && convention != ConvectionPropertyTemperature::FilmMean
        {
            return Err(ConvectionSourceError::Invalid(
                "plate sections, horizontal or spherical natural convection and aiding mixed convection require mean wall/bulk film properties",
            ));
        }
        Ok(())
    }

    // Liquid-only corrections require an actual original source, not phase
    // inference from a numeric tuple or an authored material name.
    fn admit_liquid_sources<F: ThermalFluidTransport>(
        self,
        primary: &SurfaceConvectionMaterialSample<F>,
        surface: Option<&SurfaceConvectionMaterialSample<F>>,
    ) -> Result<(), ConvectionSourceError> {
        if !matches!(self, Self::ForcedTubeSectionLiquidWall { .. }) {
            return Ok(());
        }
        let bulk = original_liquid(primary)?;
        let surface = surface.ok_or(ConvectionSourceError::Invalid(
            "liquid wall-Pr correction requires its original wall-temperature source",
        ))?;
        let wall = original_liquid(surface)?;
        let bulk_material = bulk.material();
        let wall_material = wall.material();
        let bulk_query = bulk_material.query();
        let wall_query = wall_material.query();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        let mut temperatures = bulk_query.axis_quantities().iter()
            .filter(|(_, quantity)| **quantity == absolute);
        let (temperature_axis, _) = temperatures.next().ok_or(
            ConvectionSourceError::Invalid(
                "liquid wall-Pr correction requires one explicitly typed absolute-temperature coordinate",
            ),
        )?;
        if temperatures.next().is_some()
            || bulk_query.axes().get(temperature_axis) != Some(&primary.thermal.temperature_k())
            || wall_query.axes().get(temperature_axis) != Some(&surface.thermal.temperature_k())
        {
            return Err(ConvectionSourceError::Invalid(
                "liquid wall-Pr source temperature must be unique and match each sampled state",
            ));
        }
        let bulk_properties = bulk_material.properties();
        let wall_properties = wall_material.properties();
        if bulk_material.card_identity() != wall_material.card_identity()
            || bulk.requirements() != wall.requirements()
            || bulk_properties.len() != wall_properties.len()
            || bulk_properties.iter().zip(wall_properties).any(|(a, b)| {
                a.requirement() != b.requirement()
                    || a.answer().receipt.selected != b.answer().receipt.selected
            })
            || bulk_query.axis_quantities() != wall_query.axis_quantities()
            || bulk_query.axes().len() != wall_query.axes().len()
            || bulk_query.axes().iter().any(|(axis, value)| {
                axis != temperature_axis && wall_query.axes().get(axis) != Some(value)
            })
        {
            return Err(ConvectionSourceError::Invalid(
                "liquid wall-Pr correction must retain one card, selected source bundle and non-temperature context",
            ));
        }
        Ok(())
    }

    // Reservoir boundaries and reciprocal finite links retain one law mapping
    // and one source-expansion requirement.
    pub(crate) fn resolve(
        self,
        source: &ThermalFluidTransportState,
        surface: Option<&ThermalFluidTransportState>,
    ) -> Result<CellConvectionLaw, ConvectionSourceError> {
        let expansion = || {
            source
                .volumetric_expansion_per_k()
                .ok_or(ConvectionSourceError::Invalid(
                    "natural convection requires sourced isobaric volumetric expansion",
                ))
        };
        Ok(match self {
            Self::Forced {
                correlation: CorrelationId::SiederTateTurbulent,
                speed_m_s,
                settings,
            } => CellConvectionLaw::ForcedTube {
                speed_m_s,
                length_over_diameter: settings.length_over_hydraulic_diameter.ok_or(
                    ConvectionSourceError::Invalid("Sieder-Tate requires an explicit heated length/diameter"),
                )?,
                surface_transport: surface.ok_or(ConvectionSourceError::Invalid(
                    "Sieder-Tate is missing its original fluid evaluated at wall temperature",
                ))?.transport(),
            },
            Self::Forced {
                correlation: CorrelationId::SiederTateLaminar,
                speed_m_s,
                settings,
            } => CellConvectionLaw::ForcedLaminarTube {
                speed_m_s,
                length_over_diameter: settings.length_over_hydraulic_diameter.ok_or(
                    ConvectionSourceError::Invalid("Sieder-Tate requires an explicit heated length/diameter"),
                )?,
                surface_transport: surface.ok_or(ConvectionSourceError::Invalid(
                    "Sieder-Tate is missing its original fluid evaluated at wall temperature",
                ))?.transport(),
            },
            Self::Forced {
                correlation,
                speed_m_s,
                settings,
            } => CellConvectionLaw::Forced {
                correlation,
                speed_m_s,
                settings,
            },
            Self::ForcedPlateSection { correlation, speed_m_s, start } => {
                CellConvectionLaw::ForcedPlateSection { correlation, speed_m_s, start }
            }
            Self::ForcedTubeSection { speed_m_s, start, end } => {
                CellConvectionLaw::ForcedTubeSection { speed_m_s, start, end }
            }
            Self::ForcedTubeSectionLiquidWall { speed_m_s, start, end } => {
                CellConvectionLaw::ForcedTubeSectionLiquidWall {
                    speed_m_s, start, end,
                    surface_transport: surface.ok_or(ConvectionSourceError::Invalid(
                        "liquid wall-Pr correction is missing the actual wall-temperature source",
                    ))?.transport(),
                }
            }
            Self::ForcedLaminarTubeSection { speed_m_s, start, end } => {
                CellConvectionLaw::ForcedLaminarTubeSection { speed_m_s, start, end }
            }
            Self::ForcedRectangularDuctSection { speed_m_s, aspect_ratio, start, end } => {
                CellConvectionLaw::ForcedRectangularDuctSection { speed_m_s, aspect_ratio, start, end }
            }
            Self::ForcedSphere { speed_m_s } => CellConvectionLaw::ForcedSphere {
                speed_m_s,
                surface_transport: surface.ok_or(ConvectionSourceError::Invalid(
                    "forced sphere is missing its original fluid evaluated at wall temperature",
                ))?.transport(),
            },
            Self::VerticalPlate { gravity_m_s2 } => CellConvectionLaw::VerticalPlate {
                gravity_m_s2,
                expansion_per_k: expansion()?,
            },
            Self::Natural {
                geometry,
                gravity_m_s2,
            } => CellConvectionLaw::Natural {
                geometry,
                gravity_m_s2,
                expansion_per_k: expansion()?,
            },
            Self::AidingVerticalPlate { settings } => CellConvectionLaw::AidingVerticalPlate {
                settings,
                expansion_per_k: expansion()?,
            },
        })
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SourceConvectiveCell {
    pub cell: CellId,
    pub reference_temperature: Temperature,
    pub characteristic_length: Length,
    pub property_temperature: ConvectionPropertyTemperature,
    pub law: SourceCellConvectionLaw,
}

impl SourceConvectiveCell {
    // Admit authored inputs without querying fluid properties or a not-yet-born
    // cell. The observer uses the same checks for its complete future roster.
    fn admit(&self) -> Result<(), ConvectionSourceError> {
        self.law.admit_property_temperature(self.property_temperature)?;
        self.law.admit_geometry(self.characteristic_length)?;
        if !self.reference_temperature.value().is_finite()
            || self.reference_temperature.value() <= 0.0
            || !self.characteristic_length.value().is_finite()
            || self.characteristic_length.value() <= 0.0
        {
            return Err(ConvectionSourceError::Invalid(
                "source convection needs positive finite reservoir temperature and length",
            ));
        }
        Ok(())
    }
}

/// Concrete material owner and its checked numerical view at this trial.
#[derive(Debug)]
pub struct CellConvectionMaterialSample<F> {
    pub cell: CellId,
    pub property_temperature_k: f64,
    pub material: F,
    pub thermal: ThermalFluidTransportState,
    /// Present for a selected wall-property relation: the SAME fluid at the
    /// actual wall temperature, with no separate thermal storage or mass.
    pub surface: Option<SurfaceConvectionMaterialSample<F>>,
}

/// Original fluid property owner evaluated at the surface temperature.
#[derive(Debug)]
pub struct SurfaceConvectionMaterialSample<F> {
    pub material: F,
    pub thermal: ThermalFluidTransportState,
}

// A source witness is useful only if it describes this SAME sampled state.
fn original_liquid<F: ThermalFluidTransport>(
    sample: &SurfaceConvectionMaterialSample<F>,
) -> Result<&ResolvedLiquidState, ConvectionSourceError> {
    let liquid = sample.material.resolved_liquid_source().ok_or(
        ConvectionSourceError::Invalid(
            "liquid wall-Pr correction requires an original resolved liquid source",
        ),
    )?;
    let thermal = liquid.thermal_transport_state()
        .map_err(|error| ConvectionSourceError::Material(error.into()))?;
    if thermal != sample.thermal {
        return Err(ConvectionSourceError::Invalid(
            "liquid source witness does not match the retained thermal transport state",
        ));
    }
    Ok(liquid)
}

// One source/temperature admission for reservoir, film and surface queries.
// Both finite-cell boundaries and reciprocal links retain the returned owner.
pub(crate) fn material_at<F: ThermalFluidTransport, E>(
    cell: CellId,
    temperature_k: f64,
    resolve: impl FnOnce() -> Result<F, E>,
) -> Result<SurfaceConvectionMaterialSample<F>, SourceCellConvectionError<E>> {
    let material = resolve().map_err(|source| SourceCellConvectionError::Source {
        cell, property_temperature_k: temperature_k, source,
    })?;
    let thermal = material.thermal_transport_state()
        .map_err(|error| ConvectionSourceError::Material(error.into()))?;
    let transport = material.transport_state()
        .map_err(|error| ConvectionSourceError::Material(error.into()))?;
    if thermal.temperature_k() != temperature_k || thermal.transport() != transport {
        return Err(ConvectionSourceError::Invalid(
            "fluid source does not match the requested temperature and transport state",
        ).into());
    }
    Ok(SurfaceConvectionMaterialSample { material, thermal })
}

// Resolve one physical wall state, independently of contact/event ordering.
// Birth admission uses this same material path but never applies a heat rate
// at the provisional surface area before same-time contacts have been applied.
fn resolve_cell_material<F: ThermalFluidTransport, E>(
    binding: &SourceConvectiveCell,
    wall: f64,
    resolve: &mut impl FnMut(&SourceConvectiveCell, Temperature) -> Result<F, E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<CellConvectionMaterialSample<F>, SourceCellConvectionError<E>> {
    check(cancelled)?;
    let property_temperature_k = match binding.property_temperature {
        ConvectionPropertyTemperature::Reservoir => binding.reference_temperature.value(),
        ConvectionPropertyTemperature::FilmMean => {
            0.5 * wall + 0.5 * binding.reference_temperature.value()
        }
    };
    let primary = material_at(binding.cell, property_temperature_k, || {
        resolve(binding, Temperature::new(property_temperature_k))
    })?;
    let surface = if binding.law.needs_surface_transport() {
        check(cancelled)?;
        Some(material_at(binding.cell, wall, || {
            resolve(binding, Temperature::new(wall))
        })?)
    } else { None };
    binding.law.admit_liquid_sources(&primary, surface.as_ref())?;
    Ok(CellConvectionMaterialSample {
        cell: binding.cell,
        property_temperature_k,
        material: primary.material,
        thermal: primary.thermal,
        surface,
    })
}

/// Original heat-source result together with every retained material source.
#[derive(Debug)]
pub struct ResolvedCellConvectionSources<F> {
    convection: CellConvectionSources,
    materials: Vec<CellConvectionMaterialSample<F>>,
}
impl<F> ResolvedCellConvectionSources<F> {
    pub const fn convection(&self) -> &CellConvectionSources {
        &self.convection
    }
    pub fn materials(&self) -> &[CellConvectionMaterialSample<F>] {
        &self.materials
    }
}

#[derive(Debug)]
pub enum SourceCellConvectionError<E> {
    Source {
        cell: CellId,
        property_temperature_k: f64,
        source: E,
    },
    Convection(ConvectionSourceError),
}
impl<E> From<ConvectionSourceError> for SourceCellConvectionError<E> {
    fn from(error: ConvectionSourceError) -> Self {
        Self::Convection(error)
    }
}
impl<E> From<DepositionError> for SourceCellConvectionError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Convection(ConvectionSourceError::from(error))
    }
}
impl<E: fmt::Debug> fmt::Display for SourceCellConvectionError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "source-resolved cell convection: {self:?}")
    }
}
impl<E: fmt::Debug> std::error::Error for SourceCellConvectionError<E> {}

/// Resolve source-backed transport and expansion inside the original residual.
///
/// Call with the proposed enthalpy endpoint, not only the interval's initial
/// state. The resolver receives the declared property temperature exactly once
/// per primary state, in CellId order. Spheres and Sieder-Tate tubes also query the
/// same resolver at the actual wall temperature and retain that second owner.
/// It must retain one material law across temperatures and trials
/// (for example a prepared `LiquidTemperatureLaw` or a fixed gas-mixture model)
/// and preserve its pressure, composition, selected claims and source support.
/// This generic transport interface cannot itself establish source identity.
///
/// Natural convection requires known expansion, including its sign; forced
/// convection may use a source whose expansion is unknown. Both material views
/// must describe the same transport tuple and requested temperature. Source
/// refusal, invalid bindings and cancellation return no partial result or
/// mutation of the input thermal history.
///
/// Original exposed-area accounting, correlation support, equilibrium handling
/// and boundary ownership are delegated to `evaluate_material_cell_convection`.
/// The reservoir remains prescribed; geometry is fixed, the correlation must
/// support the chosen property convention, and this is not a resolved flow or
/// evolving reservoir-energy model.
pub fn evaluate_source_cell_convection<F: ThermalFluidTransport, E>(
    thermal: &DepositionSimulation<'_>,
    bindings: &[SourceConvectiveCell],
    maximum_cells: usize,
    mut resolve: impl FnMut(&SourceConvectiveCell, Temperature) -> Result<F, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ResolvedCellConvectionSources<F>, SourceCellConvectionError<E>> {
    check(&mut cancelled)?;
    if bindings.len() > maximum_cells || bindings.len() > thermal.cells().len() {
        return Err(ConvectionSourceError::Invalid(
            "convective cells exceed their explicit work bound",
        )
        .into());
    }
    // Reject an invalid assembly before invoking any external source resolver.
    let mut ids = BTreeSet::new();
    let mut ordered = Vec::with_capacity(bindings.len());
    for binding in bindings {
        check(&mut cancelled)?;
        binding.admit()?;
        if !ids.insert(binding.cell) {
            return Err(ConvectionSourceError::Invalid("duplicate convective cell").into());
        }
        let cell = thermal
            .cells()
            .get(&binding.cell)
            .ok_or(ConvectionSourceError::Invalid("inactive convective cell"))?;
        if cell.boundary.htc_w_m2_k != 0.0 {
            return Err(ConvectionSourceError::Invalid(
                "material convection would duplicate the authored ambient coefficient",
            )
            .into());
        }
        ordered.push(binding);
    }
    ordered.sort_by_key(|binding| binding.cell);
    let mut materials = Vec::with_capacity(ordered.len());
    for binding in &ordered {
        check(&mut cancelled)?;
        let wall = thermal
            .temperature_k(binding.cell)
            .map_err(ConvectionSourceError::from)?;
        materials.push(resolve_cell_material(binding, wall, &mut resolve, &mut cancelled)?);
    }
    let mut resolved = Vec::with_capacity(ordered.len());
    for (binding, material) in ordered.iter().zip(&materials) {
        check(&mut cancelled)?;
        let law = binding.law.resolve(&material.thermal,
            material.surface.as_ref().map(|sample| &sample.thermal))?;
        resolved.push(MaterialConvectiveCell {
            cell: binding.cell,
            fluid: &material.thermal,
            reference_temperature: binding.reference_temperature,
            characteristic_length: binding.characteristic_length,
            law,
        });
    }
    let convection =
        evaluate_material_cell_convection(thermal, &resolved, maximum_cells, &mut cancelled)?;
    check(&mut cancelled)?;
    Ok(ResolvedCellConvectionSources {
        convection,
        materials,
    })
}
