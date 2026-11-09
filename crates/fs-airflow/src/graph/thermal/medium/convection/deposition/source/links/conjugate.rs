//! Source-integrated wall conduction and convection on original finite links.
//!
//! The existing material layer owner eliminates the fluid-facing temperature.
//! Existing thermal links apply the resulting combined conductance once.

use super::super::super::ConvectionSourceError;
use super::{
    ConvectiveLinkSample, SourceCellConvectionError, SourceConvectiveLink, admit_link,
    evaluate_source_link_at_surface,
};
use fs_matdb::PropertyLineSlice;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_material::thermal::conduction::{
    ConvectiveLayerError, ConvectiveLayerResponse, LayerConductionConfig, ReservoirLayerGeometry,
    source::{LayerSourceError, ResolvedConductionLayer},
};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{
    CellGeometry, CellId, ContactPatch, DepositionError, DepositionSimulation, MaterialId,
    ThermalLaw,
};
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// One existing contact and its source-qualified wall-to-surface material path.
/// The caloric object must be the exact original wall-cell law. Geometry spans
/// that cell's temperature trace to the fluid-facing surface and owns the exact
/// native contact area. This adds neither wall nor fluid storage.
#[derive(Clone, Copy)]
pub struct SourceConjugateLinkBinding<'a, 'm> {
    pub convection: SourceConvectiveLink,
    pub wall_layer: &'a ResolvedConductionLayer<'m>,
    pub geometry: ReservoirLayerGeometry,
}

/// Preserve the original source/domain failure or unclosed surface heat.
#[derive(Debug)]
pub enum SourceConjugateError<E = core::convert::Infallible> {
    Invalid(&'static str),
    Convection(SourceCellConvectionError<E>),
    Layer(LayerSourceError),
    Thermal(DepositionError),
    SurfaceBalance {
        link: u64,
        residual_w: f64,
        iterations: usize,
    },
}
impl<E> From<SourceCellConvectionError<E>> for SourceConjugateError<E> {
    fn from(error: SourceCellConvectionError<E>) -> Self {
        Self::Convection(error)
    }
}
impl<E> From<ConvectionSourceError> for SourceConjugateError<E> {
    fn from(error: ConvectionSourceError) -> Self {
        Self::Convection(error.into())
    }
}
impl<E> From<LayerSourceError> for SourceConjugateError<E> {
    fn from(error: LayerSourceError) -> Self {
        Self::Layer(error)
    }
}
impl<E> From<DepositionError> for SourceConjugateError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SourceConjugateError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source conjugate convection: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SourceConjugateError<E> {}

struct PreparedLink<'a, 'm> {
    binding: SourceConjugateLinkBinding<'a, 'm>,
    patch: ContactPatch,
    wall_cell: CellId,
    wall_geometry: CellGeometry,
    wall_material: MaterialId,
}

/// Original contact, geometry and source pins, without a frozen coefficient.
/// Retain this owner and the same fluid resolver across nonlinear trials and
/// accepted steps. The fluid callback must represent the native fluid cell's
/// original material, pressure and composition; generic transport alone does
/// not establish its identity with a caloric chart.
pub struct PreparedSourceConjugateLinks<'a, 'm> {
    links: Vec<PreparedLink<'a, 'm>>,
    convection: Vec<SourceConvectiveLink>,
}

/// One actual material/film balance; bare film hA is diagnostic, not applicable
/// as another native link. The stored wall trace and solved surface remain
/// distinct, and the original conductivity line and fluid sources are retained.
#[derive(Debug)]
pub struct SourceConjugateLinkSample<F> {
    pub link: u64,
    pub wall_cell: CellId,
    pub fluid_cell: CellId,
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    pub geometry: ReservoirLayerGeometry,
    pub source_line: Arc<PropertyLineSlice>,
    pub layer: ConvectiveLayerResponse<ConvectiveLinkSample<F>>,
    pub applied_heat_into_fluid_w: f64,
    pub representation_error_w: f64,
    pub representation_allowance_w: f64,
}

/// Only combined wall/film conductances can be applied to native reciprocal
/// links. No separate bare-film map or pair of generated source powers exists.
#[derive(Debug)]
pub struct ResolvedSourceConjugateHeat<F> {
    samples: Vec<SourceConjugateLinkSample<F>>,
    conductances_w_k: BTreeMap<u64, f64>,
}
impl<F> ResolvedSourceConjugateHeat<F> {
    pub fn samples(&self) -> &[SourceConjugateLinkSample<F>] {
        &self.samples
    }
    pub fn conductances_w_k(&self) -> &BTreeMap<u64, f64> {
        &self.conductances_w_k
    }

    /// Add the combined coefficient exactly once; conflicts and cancellation
    /// refuse before modifying the destination.
    pub fn apply_to(
        &self,
        destination: &mut BTreeMap<u64, f64>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SourceConjugateError> {
        checkpoint()?;
        for link in self.conductances_w_k.keys() {
            checkpoint()?;
            if destination.contains_key(link) {
                return Err(SourceConjugateError::Invalid(
                    "conjugate link already has a conductance owner",
                ));
            }
        }
        checkpoint()?;
        destination.extend(
            self.conductances_w_k
                .iter()
                .map(|(&id, &value)| (id, value)),
        );
        Ok(())
    }
}

impl<'a, 'm> PreparedSourceConjugateLinks<'a, 'm> {
    /// Bind every original area/contact/caloric identity before fluid queries.
    /// The layer remains homogeneous, source-free, pure-phase and quasi-steady;
    /// no additional contact resistance or lateral spreading is inferred.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        bindings: &[SourceConjugateLinkBinding<'a, 'm>],
        maximum_links: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SourceConjugateError> {
        checkpoint()?;
        if bindings.len() > maximum_links || bindings.len() > thermal.links().len() {
            return Err(SourceConjugateError::Invalid(
                "conjugate links exceed their explicit work bound",
            ));
        }
        let mut ids = BTreeSet::new();
        let mut links = Vec::with_capacity(bindings.len());
        for &binding in bindings {
            checkpoint()?;
            admit_link(thermal, &binding.convection)?;
            if !ids.insert(binding.convection.link) {
                return Err(SourceConjugateError::Invalid("duplicate conjugate link"));
            }
            let patch = thermal
                .contact_patch(binding.convection.link)
                .ok_or(SourceConjugateError::Invalid("inactive conjugate contact"))?;
            if patch.link.area_m2 != binding.geometry.surface_area_m2() {
                return Err(SourceConjugateError::Invalid(
                    "conjugate geometry must use the exact original contact area",
                ));
            }
            let wall_cell = if binding.convection.fluid_cell == patch.link.a {
                patch.link.b
            } else {
                patch.link.a
            };
            let wall = thermal
                .cells()
                .get(&wall_cell)
                .ok_or(SourceConjugateError::Invalid("inactive conjugate wall"))?;
            let prepared = PreparedLink {
                binding,
                patch,
                wall_cell,
                wall_geometry: wall.geometry,
                wall_material: wall.material_id,
            };
            prepared.temperatures::<core::convert::Infallible>(thermal)?;
            links.push(prepared);
        }
        links.sort_by_key(|link| link.binding.convection.link);
        let convection = links.iter().map(|link| link.binding.convection).collect();
        checkpoint()?;
        Ok(Self { links, convection })
    }

    /// Immutable original film declarations, in deterministic native link order.
    pub fn convection_bindings(&self) -> &[SourceConvectiveLink] {
        &self.convection
    }

    /// Re-solve each actual surface at this same proposed caloric endpoint.
    ///
    /// The existing material integral-coordinate solver owns k(T), its bracket,
    /// power residual and work bound. Every film trial uses the existing source
    /// adapter at its actual surface temperature, including film properties,
    /// wall viscosity and admitted liquid wall Pr. No nested scalar solver,
    /// extrapolation, second energy debit or extra fluid state is introduced.
    pub fn resolve<F: ThermalFluidTransport, E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        config: LayerConductionConfig,
        mut resolve_fluid: impl FnMut(
            &SourceConvectiveLink,
            &DepositionSimulation<'_>,
            Temperature,
        ) -> Result<F, E>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedSourceConjugateHeat<F>, SourceConjugateError<E>> {
        checkpoint()?;
        // Admit every saved physical identity before invoking a fluid resolver.
        let mut endpoints = Vec::with_capacity(self.links.len());
        for link in &self.links {
            checkpoint()?;
            endpoints.push(link.temperatures::<E>(thermal)?);
        }
        let checkpoint = RefCell::new(checkpoint);
        let mut conductances_w_k = BTreeMap::new();
        let mut samples = Vec::with_capacity(self.links.len());
        for (prepared, [wall, bulk]) in self.links.iter().zip(endpoints) {
            checkpoint.borrow_mut()()?;
            let binding = &prepared.binding;
            let layer = binding
                .geometry
                .evaluate_convective_material_to_bath(
                    binding.wall_layer.conductivity_curve(),
                    [wall, bulk],
                    config,
                    |surface| {
                        let film = evaluate_source_link_at_surface(
                            thermal,
                            &binding.convection,
                            surface,
                            &mut resolve_fluid,
                            || checkpoint.borrow_mut()().is_err(),
                        )?;
                        Ok::<_, SourceCellConvectionError<E>>((
                            film.heat_into_fluid_w,
                            film.conductance_w_k,
                            film,
                        ))
                    },
                    || checkpoint.borrow_mut()(),
                )
                .map_err(|error| surface_error(binding.convection.link, error))?;
            let difference = wall - bulk;
            let conductance = layer.conductance_w_k;
            let applied = conductance * difference;
            let error = applied - layer.power_to_bath_w;
            let allowance = 8.0 * f64::EPSILON * applied.abs().max(layer.power_to_bath_w.abs());
            if !conductance.is_finite()
                || conductance < 0.0
                || !applied.is_finite()
                || !error.is_finite()
                || error.abs() > allowance
                || (difference != 0.0 && applied == 0.0)
                || (difference == 0.0 && layer.power_to_bath_w != 0.0)
            {
                return Err(DepositionError::Arithmetic.into());
            }
            conductances_w_k.insert(binding.convection.link, conductance);
            samples.push(SourceConjugateLinkSample {
                link: binding.convection.link,
                wall_cell: prepared.wall_cell,
                fluid_cell: binding.convection.fluid_cell,
                wall_temperature_k: wall,
                fluid_temperature_k: bulk,
                geometry: binding.geometry,
                source_line: Arc::clone(binding.wall_layer.source_line()),
                layer,
                applied_heat_into_fluid_w: applied,
                representation_error_w: error,
                representation_allowance_w: allowance,
            });
        }
        checkpoint.borrow_mut()()?;
        Ok(ResolvedSourceConjugateHeat {
            samples,
            conductances_w_k,
        })
    }
}

impl PreparedLink<'_, '_> {
    fn temperatures<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<[f64; 2], SourceConjugateError<E>> {
        if thermal.contact_patch(self.binding.convection.link) != Some(self.patch) {
            return Err(SourceConjugateError::Invalid(
                "conjugate contact changed; reprepare its geometry",
            ));
        }
        let wall = thermal
            .cells()
            .get(&self.wall_cell)
            .ok_or(SourceConjugateError::Invalid("inactive conjugate wall"))?;
        if wall.geometry != self.wall_geometry
            || wall.material_id != self.wall_material
            || !thermal.uses_material_law(
                wall.material_id,
                self.binding.wall_layer.caloric_curve() as &dyn ThermalLaw,
            )
        {
            return Err(SourceConjugateError::Invalid(
                "conjugate wall requires its original geometry and caloric object",
            ));
        }
        // Retain the original pure-phase/enthalpy admission. Its restricted
        // source endpoints may snap an inverse round trip to a source bound;
        // native heat instead uses the actual original caloric temperature.
        let supported = self
            .binding
            .wall_layer
            .temperature_at_enthalpy(wall.specific_h_j_kg)?;
        let actual = thermal.temperature_k(self.wall_cell)?;
        if actual != supported {
            // Admit the UNSNAPPED native trace or refuse. Never extend source
            // support or change the temperature used by the native G*DeltaT.
            self.binding
                .wall_layer
                .conductivity_curve()
                .conductivity_w_m_k(actual)?;
        }
        Ok([
            actual,
            thermal.temperature_k(self.binding.convection.fluid_cell)?,
        ])
    }
}

fn surface_error<E>(
    link: u64,
    error: ConvectiveLayerError<SourceCellConvectionError<E>>,
) -> SourceConjugateError<E> {
    match error {
        ConvectiveLayerError::Thermal(error) => SourceConjugateError::Thermal(error),
        ConvectiveLayerError::Film(error) => SourceConjugateError::Convection(error),
        ConvectiveLayerError::SurfaceBalance {
            residual_w,
            iterations,
        } => SourceConjugateError::SurfaceBalance {
            link,
            residual_w,
            iterations,
        },
    }
}
