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

use super::*;
use fs_material::fluid_transport::{ThermalFluidTransport, ThermalFluidTransportState};

/// Authored property-evaluation convention required by the selected correlation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConvectionPropertyTemperature {
    Reservoir,
    FilmMean,
}

/// Expansion is supplied by the same resolved material as the transport tuple.
#[derive(Clone, Copy, Debug)]
pub enum SourceCellConvectionLaw {
    Forced {
        correlation: CorrelationId,
        speed_m_s: f64,
        settings: ForcedConvectionSettings,
    },
    VerticalPlate {
        gravity_m_s2: f64,
    },
    /// Geometry selects the shared natural-convection card. Horizontal plate
    /// and cylinder cards require `FilmMean` properties; Rayleigh and heat use
    /// the full wall/bulk temperature difference. Expansion comes from the
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
    fn admit_property_temperature(
        self,
        convention: ConvectionPropertyTemperature,
    ) -> Result<(), ConvectionSourceError> {
        if matches!(
            self,
            Self::Natural {
                geometry: NaturalConvectionGeometry::HorizontalCylinder
                    | NaturalConvectionGeometry::HorizontalPlate(_),
                ..
            }
                | Self::AidingVerticalPlate { .. }
        ) && convention != ConvectionPropertyTemperature::FilmMean
        {
            return Err(ConvectionSourceError::Invalid(
                "horizontal natural and aiding mixed convection require mean wall/bulk film properties",
            ));
        }
        Ok(())
    }

    // Reservoir boundaries and reciprocal finite links retain one law mapping
    // and one source-expansion requirement.
    fn resolve(
        self,
        source: &ThermalFluidTransportState,
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
                correlation,
                speed_m_s,
                settings,
            } => CellConvectionLaw::Forced {
                correlation,
                speed_m_s,
                settings,
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

/// Concrete material owner and its checked numerical view at this trial.
#[derive(Debug)]
pub struct CellConvectionMaterialSample<F> {
    pub cell: CellId,
    pub property_temperature_k: f64,
    pub material: F,
    pub thermal: ThermalFluidTransportState,
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
/// per cell, in CellId order. It must retain one material law across trials
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
        binding
            .law
            .admit_property_temperature(binding.property_temperature)?;
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
        if !binding.reference_temperature.value().is_finite()
            || binding.reference_temperature.value() <= 0.0
            || !binding.characteristic_length.value().is_finite()
            || binding.characteristic_length.value() <= 0.0
        {
            return Err(ConvectionSourceError::Invalid(
                "source convection needs positive finite reservoir temperature and length",
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
        let property_temperature_k = match binding.property_temperature {
            ConvectionPropertyTemperature::Reservoir => binding.reference_temperature.value(),
            ConvectionPropertyTemperature::FilmMean => {
                0.5 * wall + 0.5 * binding.reference_temperature.value()
            }
        };
        let material =
            resolve(binding, Temperature::new(property_temperature_k)).map_err(|source| {
                SourceCellConvectionError::Source {
                    cell: binding.cell,
                    property_temperature_k,
                    source,
                }
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
        materials.push(CellConvectionMaterialSample {
            cell: binding.cell,
            property_temperature_k,
            material,
            thermal: snapshot,
        });
    }
    let mut resolved = Vec::with_capacity(ordered.len());
    for (binding, material) in ordered.iter().zip(&materials) {
        check(&mut cancelled)?;
        let law = binding.law.resolve(&material.thermal)?;
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
