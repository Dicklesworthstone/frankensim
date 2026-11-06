//! Finite-strain membranes from the existing three-dimensional energy owner.
//!
//! Material coordinates are warp, weft and reference normal. A positive
//! Cholesky representative of `C_in_plane = I + 2 E` removes rigid rotation;
//! native plane stress either relaxes the entire transverse column or fixes
//! elastic `det(Fe)=1` while relaxing both transverse shears. With a selected
//! natural reference G, that constraint gives total `det(F)=det(G)`. Multiplying
//! the reduced potential by the original thickness gives work per reference
//! area. Its exact chain rule supplies second-Piola resultants and their
//! Green-strain tangent, including the stress times second-kinematic-derivative
//! term. The physical membrane element retains its own geometric stiffness.

use super::{
    AREAL_DENSITY, MEMBRANE_PROPERTIES, ResolvedSheet, SheetError, SheetLaw, SheetResponse,
    THICKNESS,
};
use crate::hyper::cards::ResolvedHyperelastic;
use crate::hyper::mixed::VolumetricTreatment;
use crate::hyper::plane_stress::{PlaneStressResponse, PlaneStressSettings};
use crate::hyper::{FibreInvariants, HyperelasticModel};
use crate::state_point::{
    MeanThermalExpansionStrain, ResolvedMaterialStatePoint, is_fixed_external_coordinate,
};
use fs_ad::{Dual, Dual64, Real};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

/// Immutable, source-bound hyperelastic membrane at fixed external conditions.
///
/// Reference mass is independent of the bulk law and thickness. The original
/// source bundle and exact model card survive extraction into `MembraneLaw`.
/// This is the native locally stable plane-stress branch, with the original
/// finite compliance or an explicitly imposed exact volume constraint. The
/// constraint is a modeling choice, not additional source evidence. No wrinkling,
/// bending or residual stress is inferred. An explicit natural reference may
/// map the mesh reference to the model's stress-free material frame; arbitrary
/// fibre directions remain in that frame.
#[derive(Clone, Debug)]
pub struct HyperelasticMembrane {
    material: ResolvedMaterialStatePoint,
    model: ResolvedHyperelastic,
    reference_thickness_m: f64,
    settings: PlaneStressSettings,
    volumetric: VolumetricTreatment,
    natural_reference: Option<NaturalReference>,
}

#[derive(Clone, Debug, PartialEq)]
struct NaturalReference {
    stretches: [f64; 3],
    inverse: [f64; 3],
    volume_ratio: f64,
    thermal: Option<Box<MeanThermalSource>>,
}

#[derive(Clone, Debug, PartialEq)]
struct MeanThermalSource {
    expansion: [MeanThermalExpansionStrain; 3],
    temperature_axes: [String; 2],
}

impl PartialEq for HyperelasticMembrane {
    fn eq(&self, other: &Self) -> bool {
        self.material == other.material
            && self.model == other.model
            && self.reference_thickness_m == other.reference_thickness_m
            && self.volumetric == other.volumetric
            && self.natural_reference == other.natural_reference
            && self.settings.relative_tolerance == other.settings.relative_tolerance
            && self.settings.max_iterations == other.settings.max_iterations
            && self.settings.max_backtracks == other.settings.max_backtracks
    }
}

impl HyperelasticMembrane {
    pub(super) fn from_material_source(
        sheet: &ResolvedSheet,
        card: &MaterialCard,
        model_index: usize,
        settings: PlaneStressSettings,
        volumetric: VolumetricTreatment,
    ) -> Result<Self, SheetError> {
        if sheet.material.card_identity() != card.content_hash() {
            return Err(SheetError::Invalid(
                "sheet and hyperelastic source card differ",
            ));
        }
        if sheet.hyperelastic.is_some()
            || sheet.tabulated.is_some()
            || MEMBRANE_PROPERTIES
                .iter()
                .any(|k| sheet.material.property(k).is_some())
            || super::tabulated::CURVE_PROPERTIES
                .iter()
                .chain([&super::tabulated::COUPLING_PROPERTY])
                .any(|k| sheet.material.property(k).is_some())
        {
            return Err(SheetError::Invalid(
                "hyperelastic selection cannot replace or add to an existing membrane capability",
            ));
        }
        for key in [AREAL_DENSITY, THICKNESS] {
            if !sheet
                .material
                .property(key)
                .is_some_and(|p| p.is_static_at_fixed_external_conditions())
            {
                return Err(SheetError::Invalid(
                    "hyperelastic reference mass and thickness require fixed external source conditions",
                ));
            }
        }
        let card = card.models().get(model_index).ok_or(SheetError::Invalid(
            "hyperelastic model index is absent from the sheet source card",
        ))?;
        let model = ResolvedHyperelastic::try_from_model_card_at(card, &sheet.material.query())?;
        let law = Self {
            material: sheet.material.clone(),
            model,
            reference_thickness_m: sheet.thickness,
            settings,
            volumetric,
            natural_reference: None,
        };
        // Validate local controls and the supplied model at its reference.
        // HGO may have a one-sided tangent here, so this is not a vibration
        // admission or a certificate of stability at every later state.
        law.response([0.0; 3])?;
        Ok(law)
    }

    /// Complete original property bundle, including independent mass receipts.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Exact selected model card and held typed context. A thermal reference
    /// retains the original mass/thickness context in `material()` and binds
    /// this model to its explicitly sourced endpoint temperature.
    pub const fn model(&self) -> &ResolvedHyperelastic {
        &self.model
    }

    /// Original reference thickness [m], never the locally relaxed thickness.
    pub const fn reference_thickness_m(&self) -> f64 {
        self.reference_thickness_m
    }

    /// Bounded controls used by the original local plane-stress owner.
    pub const fn settings(&self) -> PlaneStressSettings {
        self.settings
    }

    /// Original finite-compliance law or explicitly imposed elastic `det(Fe)=1`.
    /// A selected natural reference G then gives total `det(F)=det(G)`.
    /// Exact incompressibility retains the source model's native remainder;
    /// it never evaluates or substitutes its finite bulk penalty.
    pub const fn volumetric_treatment(&self) -> VolumetricTreatment {
        self.volumetric
    }

    /// Stress-free lengths per original mesh-reference length, in the fixed
    /// warp, weft and normal axes. The unmodified reference returns [1; 3].
    /// Original thickness and areal mass remain measured in the mesh reference.
    pub fn reference_stretches(&self) -> [f64; 3] {
        self.natural_reference
            .as_ref()
            .map_or([1.0; 3], |r| r.stretches)
    }

    /// Original finite endpoint expansion sources in warp/weft/normal order.
    pub fn mean_thermal_expansion(&self) -> Option<&[MeanThermalExpansionStrain; 3]> {
        self.natural_reference
            .as_ref()?
            .thermal
            .as_ref()
            .map(|t| &t.expansion)
    }

    /// Exact reference/endpoint axis roles of the retained mean sources.
    pub fn mean_thermal_temperature_axes(&self) -> Option<[&str; 2]> {
        self.natural_reference
            .as_ref()?
            .thermal
            .as_ref()
            .map(|t| t.temperature_axes.each_ref().map(String::as_str))
    }

    pub(super) fn with_reference_stretches(
        mut self,
        stretches: [f64; 3],
    ) -> Result<Self, SheetError> {
        let inverse = stretches.map(f64::recip);
        let volume_ratio = stretches[0] * stretches[1] * stretches[2];
        if self.natural_reference.is_some() {
            return Err(SheetError::Invalid(
                "a membrane natural reference may be selected only once",
            ));
        }
        if stretches
            .iter()
            .chain(&inverse)
            .chain([&volume_ratio])
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(SheetError::Invalid(
                "finite positive natural lengths, inverse and volume required",
            ));
        }
        self.natural_reference = Some(NaturalReference {
            stretches,
            inverse,
            volume_ratio,
            thermal: None,
        });
        Ok(self)
    }

    pub(super) fn with_mean_thermal_expansion(
        self,
        card: &MaterialCard,
        expansion: [MeanThermalExpansionStrain; 3],
        reference_axis: &str,
        endpoint_axis: &str,
    ) -> Result<Self, SheetError> {
        let point = expansion[0].resolved().query();
        let reference_temperature = expansion[0].reference_temperature_k();
        let endpoint_temperature = expansion[0].endpoint_temperature_k();
        if card.content_hash() != self.material.card_identity()
            || reference_axis == endpoint_axis
            || point.axes().contains_key("T") && reference_axis != "T" && endpoint_axis != "T"
        {
            return Err(SheetError::Invalid(
                "same material and unambiguous thermal temperature roles required",
            ));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        for (axis, value) in [
            (reference_axis, reference_temperature),
            (endpoint_axis, endpoint_temperature),
        ] {
            if point.axes().get(axis) != Some(&value)
                || !point.axis_quantities().get(axis).is_some_and(|q| {
                    *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS)
                })
            {
                return Err(SheetError::Invalid(
                    "thermal source roles must retain their typed endpoint temperatures",
                ));
            }
        }
        for source in &expansion {
            if source.resolved().card_identity() != card.content_hash()
                || source.resolved().query() != point
                || source.reference_temperature_k() != reference_temperature
                || source.endpoint_temperature_k() != endpoint_temperature
            {
                return Err(SheetError::Invalid(
                    "all three thermal sources must share the exact card and context",
                ));
            }
            let property = &source.resolved().properties()[0];
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(SheetError::Invalid(
                    "thermal reference lost its original selected claim",
                ))?;
            for axis in claim
                .validity
                .bounds()
                .keys()
                .chain(property.source_axes().iter().map(|(axis, _)| axis))
            {
                if axis != reference_axis
                    && axis != endpoint_axis
                    && !is_fixed_external_coordinate(
                        axis,
                        point.axis_quantities().get(axis).copied(),
                        point.axes().get(axis).copied(),
                    )
                {
                    return Err(SheetError::Invalid(
                        "thermal reference requires fixed external source coordinates",
                    ));
                }
            }
        }
        let context = |temperature_axis: &str| -> Result<QueryPoint, SheetError> {
            let mut context = QueryPoint::new();
            for (axis, value) in point.axes() {
                if axis == reference_axis || axis == endpoint_axis {
                    continue;
                }
                context = match point.axis_quantities().get(axis) {
                    Some(quantity) => context.with_quantity(axis, *quantity, *value),
                    None => context.with(axis, *value),
                }
                .map_err(|_| SheetError::Invalid("unrepresentable fixed thermal context"))?;
            }
            context
                .with_quantity(
                    "T",
                    point.axis_quantities()[temperature_axis],
                    point.axes()[temperature_axis],
                )
                .map_err(|_| SheetError::Invalid("unrepresentable thermal endpoint context"))
        };
        let reference_context = context(reference_axis)?;
        if self.material.query() != reference_context || *self.model.context() != reference_context
        {
            return Err(SheetError::Invalid(
                "sheet mass, thickness and initial model must belong to the thermal reference context",
            ));
        }
        let model = ResolvedHyperelastic::try_from_model_card_at(
            self.model.model_card(),
            &context(endpoint_axis)?,
        )?;
        let stretches = expansion
            .each_ref()
            .map(|source| 1.0 + source.engineering_strain());
        let mut result = self.with_reference_stretches(stretches)?;
        result.model = model;
        result
            .natural_reference
            .as_mut()
            .expect("reference just selected")
            .thermal = Some(Box::new(MeanThermalSource {
            expansion,
            temperature_axes: [reference_axis.into(), endpoint_axis.into()],
        }));
        Ok(result)
    }

    fn reduced_response(&self, strain: [f64; 3]) -> Result<PlaneStressResponse, SheetError> {
        let mut f = admitted_representative(strain)?;
        // G maps the mesh reference to the native model's natural frame.
        // Preserve the original arithmetic when no reference is selected.
        if let Some(reference) = &self.natural_reference {
            if reference.stretches != [1.0; 3] {
                for (a, value) in f.iter_mut().enumerate() {
                    *value *= reference.inverse[a % 2];
                }
            }
        }
        let mut response = match self.volumetric {
            VolumetricTreatment::OriginalMaterial => {
                self.model.law().plane_stress(&f, self.settings)?
            }
            VolumetricTreatment::ExactIncompressible => {
                let response = self
                    .model
                    .law()
                    .incompressible_plane_stress(&f, self.settings)?;
                PlaneStressResponse {
                    deformation_gradient: response.deformation_gradient,
                    energy_density_pa: response.energy_density_pa,
                    piola: response.piola,
                    tangent: response.tangent,
                    relative_transverse_residual: response.relative_transverse_residual,
                    iterations: response.iterations,
                }
            }
        };
        if let Some(reference) = &self.natural_reference {
            if reference.stretches != [1.0; 3] {
                // Native transverse relaxation commutes with this diagonal
                // pullback. Exact incompressibility constrains det(Fe), not
                // total volume. Keep Fe for elastic stretch/fibre admission.
                response.energy_density_pa *= reference.volume_ratio;
                for a in 0..4 {
                    let scale = reference.volume_ratio * reference.inverse[a % 2];
                    response.piola[a] *= scale;
                    for b in 0..4 {
                        response.tangent[a][b] *= scale * reference.inverse[b % 2];
                    }
                }
                if !response.energy_density_pa.is_finite()
                    || response
                        .piola
                        .iter()
                        .chain(response.tangent.iter().flatten())
                        .any(|value| !value.is_finite())
                {
                    return Err(SheetError::Invalid(
                        "unrepresentable natural-reference membrane pullback",
                    ));
                }
            }
        }
        Ok(response)
    }

    /// Require a two-sided material tangent for small-vibration extraction.
    ///
    /// Native nonlinear equilibrium retains HGO's one-sided fibre activation
    /// convention. Here a conservative roundoff neighborhood of each activation
    /// surface refuses; no artificial tensile stiffness or smoothing is added.
    pub fn admit_linearization(&self, strain: [f64; 3]) -> Result<(), SheetError> {
        let response = self.reduced_response(strain)?;
        let law = self.model.law();
        let stretches = law.principal_stretches(&response.deformation_gradient)?;
        let margin = 64.0 * f64::EPSILON;
        if stretches[0] <= (1.0 + margin) / law.stretch_limit
            || stretches[2] >= (1.0 - margin) * law.stretch_limit
        {
            return Err(SheetError::Invalid(
                "small vibration requires two-sided support inside the hyperelastic stretch domain",
            ));
        }
        if let HyperelasticModel::HolzapfelGasserOgden {
            k1,
            dispersion,
            fibres,
            invariants,
            ..
        } = &law.model
        {
            if *k1 > 0.0 {
                let f = response.deformation_gradient;
                let i1: f64 = f.iter().map(|v| v * v).sum();
                // The accepted representative is upper triangular.
                let j = f[0] * f[4] * f[8];
                let scale = match invariants {
                    FibreInvariants::Full => 1.0,
                    FibreInvariants::Isochoric => {
                        fs_math::det::exp((-2.0 / 3.0) * fs_math::det::ln(j))
                    }
                };
                for a in fibres {
                    let i4: f64 = (0..3)
                        .map(|row| {
                            let fa: f64 = (0..3).map(|col| f[3 * row + col] * a[col]).sum();
                            fa * fa
                        })
                        .sum();
                    let e = dispersion * (scale * i1 - 3.0)
                        + (1.0 - 3.0 * dispersion) * (scale * i4 - 1.0);
                    let guard = 64.0
                        * f64::EPSILON
                        * (1.0
                            + dispersion * (scale * i1).abs()
                            + (1.0 - 3.0 * dispersion) * (scale * i4).abs());
                    if !e.is_finite() || e.abs() <= guard {
                        return Err(SheetError::Invalid(
                            "no resolved two-sided membrane tangent at HGO fibre activation",
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

// F2^T F2 = [[1+2 E11, 2 E12], [2 E12, 1+2 E22]].
// The third engineering coordinate is already 2 E12, without another factor.
fn representative<T: Real>(strain: [T; 3]) -> [T; 4] {
    let two = T::from_f64(2.0);
    let r = (T::one() + two * strain[0]).sqrt();
    let s = strain[2] / r;
    let t = (T::one() + two * strain[1] - s * s).sqrt();
    [r, s, T::zero(), t]
}

pub(crate) fn admitted_representative(strain: [f64; 3]) -> Result<[f64; 4], SheetError> {
    if strain.iter().any(|v| !v.is_finite()) {
        return Err(SheetError::Invalid("finite Green strain required"));
    }
    let f = representative(strain);
    if f.iter().any(|v| !v.is_finite()) || f[0] <= 0.0 || f[3] <= 0.0 {
        return Err(SheetError::Invalid(
            "membrane Green strain must define a positive-definite in-plane metric",
        ));
    }
    Ok(f)
}

impl SheetLaw for HyperelasticMembrane {
    fn response(&self, strain: [f64; 3]) -> Result<SheetResponse, SheetError> {
        let response = self.reduced_response(strain)?;
        response_from_plane_stress(strain, response, self.reference_thickness_m)
    }
}

// Shared by elastic and frozen-history membrane potentials. The admitted
// metric and relaxed three-dimensional response must refer to the SAME strain.
// This changes kinematic coordinates only; it cannot change the native energy,
// material history, transverse branch, or physical meaning of the potential.
pub(crate) fn response_from_plane_stress(
    strain: [f64; 3],
    response: PlaneStressResponse,
    reference_thickness_m: f64,
) -> Result<SheetResponse, SheetError> {
    let mut jacobian = [[0.0; 3]; 4];
    let mut hessian = [[[0.0; 3]; 3]; 4];
    // Three small nested-AD sweeps differentiate only the kinematic map.
    // The constitutive update and transverse Schur complement are shared
    // with all native plane-stress consumers, without a second Newton solve.
    for col in 0..3 {
        let vars = core::array::from_fn(|row| Dual::<Dual64<3>, 1> {
            re: Dual64::variable(strain[row], row),
            eps: [Dual64::constant(if row == col { 1.0 } else { 0.0 })],
        });
        let f = representative(vars);
        for a in 0..4 {
            if col == 0 {
                jacobian[a] = f[a].re.eps;
            }
            for row in 0..3 {
                hessian[a][row][col] = f[a].eps[0].eps[row];
            }
        }
    }
    let mut resultants_n_m = [0.0; 3];
    let mut tangent_n_m = [0.0; 9];
    for i in 0..3 {
        for a in 0..4 {
            resultants_n_m[i] += response.piola[a] * jacobian[a][i];
            for j in 0..3 {
                tangent_n_m[3 * i + j] += response.piola[a] * hessian[a][i][j];
                for b in 0..4 {
                    tangent_n_m[3 * i + j] +=
                        jacobian[a][i] * response.tangent[a][b] * jacobian[b][j];
                }
            }
        }
        resultants_n_m[i] *= reference_thickness_m;
        for j in 0..3 {
            tangent_n_m[3 * i + j] *= reference_thickness_m;
        }
    }
    let energy_j_m2 = reference_thickness_m * response.energy_density_pa;
    if !energy_j_m2.is_finite()
        || resultants_n_m
            .iter()
            .chain(&tangent_n_m)
            .any(|v| !v.is_finite())
    {
        return Err(SheetError::Invalid(
            "unrepresentable finite-strain membrane response",
        ));
    }
    Ok(SheetResponse {
        energy_j_m2,
        resultants_n_m,
        tangent_n_m,
    })
}
