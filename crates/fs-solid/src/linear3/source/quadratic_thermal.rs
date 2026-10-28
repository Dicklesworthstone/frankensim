//! Nodal temperatures lowered through the original source expansion integrals.
//!
//! Interpolate temperature first, then integrate the selected expansion law at
//! each native stiffness point. Interpolating nodal eigenstrains instead gives
//! a different constitutive model whenever alpha depends on temperature.

use super::{SolidMaterialError, ThermalExpansionModel, ThermalExpansionSource};
use crate::linear3::{
    QuadraticThermalStrainField, TetElasticError, TetMaterialField, TetThermalStrainState,
};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};

/// Absolute nodal temperatures [K], in the structural mesh's node order.
///
/// `reference_k` is the original stress-free field, not the previous time step.
/// Both fields use the native quadratic interpolation, including midside nodes.
/// A different thermal discretization must first supply its explicit transfer
/// onto these nodes; neither enthalpy nor a conservative transfer is inferred.
#[derive(Clone, Copy, Debug)]
pub struct NodalThermalTemperatures<'a> {
    pub reference_k: &'a [f64],
    pub current_k: &'a [f64],
}

/// Original expansion source and its fixed material axes for one element.
///
/// The reference must come from this card and contain the expansion schema of
/// `model`. Its selected claims stay pinned at every integration point. Its
/// nonthermal coordinates and quantity conventions remain unchanged; only its
/// `T` coordinate is replaced by the interpolated reference/current fields.
/// Elasticity remains the problem's explicitly supplied per-element tensor.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticThermalSource<'a> {
    pub card: &'a MaterialCard,
    pub reference_expansion: &'a ResolvedMaterialStatePoint,
    pub model: ThermalExpansionModel,
}

/// Source assignments in element order, independent of nodal temperature order.
#[derive(Clone, Copy, Debug)]
pub enum QuadraticThermalSources<'a> {
    Uniform(&'a QuadraticThermalSource<'a>),
    PerElement(&'a [QuadraticThermalSource<'a>]),
}

impl QuadraticThermalSources<'_> {
    fn at(&self, element: usize) -> &QuadraticThermalSource<'_> {
        match self {
            Self::Uniform(source) => source,
            Self::PerElement(sources) => &sources[element],
        }
    }
}

/// Owned source receipts and native eigenstrains in stiffness-quadrature order.
///
/// A borrow converts directly into the existing load, equilibrium and vibration
/// tangent APIs. Re-resolve after changing temperature, quadrature, element
/// order, geometry, material identity or material axes. Reference mass and the
/// supplied elastic tensors are untouched. This is prescribed-temperature
/// loading, not heat transport or reciprocal thermoelastic heat feedback.
#[derive(Clone, Debug)]
pub struct ResolvedQuadraticThermalStrains {
    sources: Vec<ThermalExpansionSource>,
    strains: Vec<TetThermalStrainState>,
    points_per_element: usize,
}

impl ResolvedQuadraticThermalStrains {
    /// Original integrated sources, element-major and then quadrature-major.
    #[must_use]
    pub fn sources(&self) -> &[ThermalExpansionSource] {
        &self.sources
    }

    #[must_use]
    pub fn strains(&self) -> &[TetThermalStrainState] {
        &self.strains
    }

    #[must_use]
    pub const fn points_per_element(&self) -> usize {
        self.points_per_element
    }

    #[must_use]
    pub fn field(&self) -> QuadraticThermalStrainField<'_> {
        QuadraticThermalStrainField::AtQuadrature {
            states: &self.strains,
            points_per_element: self.points_per_element,
        }
    }
}

impl<'a> From<&'a ResolvedQuadraticThermalStrains> for QuadraticThermalStrainField<'a> {
    fn from(resolved: &'a ResolvedQuadraticThermalStrains) -> Self {
        resolved.field()
    }
}

fn invalid(what: &'static str) -> SolidMaterialError {
    TetElasticError::InvalidThermalStrain { what }.into()
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidMaterialError> {
    cx.checkpoint()
        .map_err(|_| TetElasticError::Cancelled.into())
}

fn temperature_point(
    template: &QueryPoint,
    temperature_k: f64,
) -> Result<QueryPoint, SolidMaterialError> {
    if !temperature_k.is_finite() || temperature_k <= 0.0 || !template.axes().contains_key("T") {
        return Err(invalid(
            "positive finite interpolated absolute T and the original T coordinate required",
        ));
    }
    template
        .clone()
        .with("T", temperature_k)
        .map_err(|_| invalid("could not retain thermal source coordinates"))
}

// Geometry/topology/mesh budgets have already passed the native problem's
// admission. The callback is that family's existing point evaluator; no second
// shape, Jacobian, quadrature or material-integration implementation lives here.
#[allow(clippy::too_many_arguments)]
pub(crate) fn resolve<const N: usize, const D: usize>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    temperatures: NodalThermalTemperatures<'_>,
    bindings: QuadraticThermalSources<'_>,
    quadrature: &[([f64; D], f64)],
    mesh_frame: ContentHash,
    strain_limit: f64,
    point: impl Fn([[f64; 3]; N], [f64; D], usize, usize) -> Result<[f64; N], TetElasticError>,
    cx: &Cx<'_>,
) -> Result<ResolvedQuadraticThermalStrains, SolidMaterialError> {
    checkpoint(cx)?;
    if temperatures.reference_k.len() != nodes_m.len()
        || temperatures.current_k.len() != nodes_m.len()
    {
        return Err(invalid(
            "one reference and current temperature per structural node required",
        ));
    }
    if let QuadraticThermalSources::PerElement(sources) = bindings
        && sources.len() != elements.len()
    {
        return Err(invalid("one expansion source per element required"));
    }
    if mesh_frame == ContentHash([0; 32])
        || !strain_limit.is_finite()
        || strain_limit <= 0.0
        || strain_limit >= 1.0
    {
        return Err(invalid(
            "nonzero mesh frame and a finite small-strain limit in (0,1) required",
        ));
    }
    for (node, (&reference, &current)) in temperatures
        .reference_k
        .iter()
        .zip(temperatures.current_k)
        .enumerate()
    {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        if !reference.is_finite() || reference <= 0.0 || !current.is_finite() || current <= 0.0 {
            return Err(invalid(
                "nodal absolute temperatures must be positive and finite",
            ));
        }
    }
    let count = elements
        .len()
        .checked_mul(quadrature.len())
        .ok_or(TetElasticError::SizeOverflow)?;
    let mut sources: Vec<ThermalExpansionSource> = Vec::with_capacity(count);
    let mut strains = Vec::with_capacity(count);
    for (element, connectivity) in elements.iter().enumerate() {
        checkpoint(cx)?;
        let binding = bindings.at(element);
        if binding.reference_expansion.card_identity() != binding.card.content_hash() {
            return Err(invalid(
                "expansion reference must belong to the original material card",
            ));
        }
        // Prepare the exact source selection once per element. The existing
        // source owner performs the actual supported-path integral at each T.
        let selection = MaterialPropertySelection::PinnedByProperty(
            binding
                .reference_expansion
                .properties()
                .iter()
                .map(|property| {
                    (
                        property.requirement().name().to_owned(),
                        property.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let template = binding.reference_expansion.query();
        let reference_anchor = temperatures.reference_k[connectivity[0]];
        let current_anchor = temperatures.current_k[connectivity[0]];
        let uniform = connectivity.iter().all(|&node| {
            temperatures.reference_k[node].to_bits() == reference_anchor.to_bits()
                && temperatures.current_k[node].to_bits() == current_anchor.to_bits()
        });
        let start = sources.len();
        let nodes = connectivity.map(|node| nodes_m[node]);
        for (index, &(coordinate, _)) in quadrature.iter().enumerate() {
            checkpoint(cx)?;
            let shapes = point(nodes, coordinate, element, index)?;
            let source = if uniform && index != 0 {
                // The common constant-field path resolves its law only once.
                sources[start].clone()
            } else {
                // Anchoring removes a rounded partition-of-unity residual
                // times the absolute temperature and preserves uniform bits.
                let mut reference_k = reference_anchor;
                let mut current_k = current_anchor;
                for local in 0..N {
                    reference_k += shapes[local]
                        * (temperatures.reference_k[connectivity[local]] - reference_anchor);
                    current_k += shapes[local]
                        * (temperatures.current_k[connectivity[local]] - current_anchor);
                }
                // Higher-order shape functions need not be positive. Reject
                // nonphysical/source-domain overshoot; never clamp or average it.
                let reference = temperature_point(&template, reference_k)?;
                let current = temperature_point(&template, current_k)?;
                binding.model.resolve(
                    binding.card,
                    &reference,
                    &current,
                    selection.clone(),
                    materials.at(element),
                    mesh_frame,
                    strain_limit,
                )?
            };
            strains.push(*source.thermal_strain());
            sources.push(source);
        }
    }
    checkpoint(cx)?;
    Ok(ResolvedQuadraticThermalStrains {
        sources,
        strains,
        points_per_element: quadrature.len(),
    })
}
