//! Oriented thermoelasticity from the original scalar material claims.
//!
//! Affine E_i, nu_ij and G_ij do NOT make the stiffness affine. The native
//! orthotropic compliance is inverted by `OrthotropicElastic`; its thermal
//! derivatives use C'=-C S' C and C''=2 C S' C S' C-C S'' C. Independently
//! oriented instantaneous expansion is integrated by the original source
//! owner. All derivatives are at fixed TOTAL Mandel strain in the target frame.
//!
//! Whole-interval admission uses outward-rounded coefficient bounds, a
//! sufficient positive-compliance determinant bound and an operator-norm bound
//! on the elastic heat capacity. It may refuse a valid but poorly bounded
//! source interval; endpoint samples alone never establish that admission.

use super::*;
use crate::Tangent6;
use crate::elastic::OrthotropicElastic;
use crate::state_point::{
    ORTHOTROPIC_POISSON_RATIO_PROPERTIES, ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES, integrate_orthotropic_thermal_expansion,
};
use crate::tensor::{mandel_rotation, rotate_mandel};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Explicit source schema and fixed material axes. Elastic and thermal
/// symmetries are independent; neither is inferred from a material name.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThermoelasticSymmetry {
    Isotropic,
    /// Columns are the original principal directions in the target frame.
    OrientedOrthotropic { principal_to_target: [[f64; 3]; 3] },
}

impl ThermoelasticSymmetry {
    fn rotation(self) -> [[f64; 3]; 3] {
        match self {
            Self::Isotropic => IDENTITY,
            Self::OrientedOrthotropic { principal_to_target } => principal_to_target,
        }
    }
    fn admit(self) -> Result<(), ThermoelasticError> {
        let q = self.rotation();
        if q.iter().flatten().any(|x| !x.is_finite()) {
            return Err(ThermoelasticError::Invalid("finite thermoelastic material axes required"));
        }
        for i in 0..3 {
            for j in 0..3 {
                let dot = (0..3).map(|k| q[i][k] * q[j][k]).sum::<f64>();
                if (dot - if i == j { 1.0 } else { 0.0 }).abs() > 1.0e-12 {
                    return Err(ThermoelasticError::Invalid("orthonormal thermoelastic material axes required"));
                }
            }
        }
        let det = q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1])
            - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0])
            + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]);
        if !det.is_finite() || (det - 1.0).abs() > 1.0e-12 {
            return Err(ThermoelasticError::Invalid("proper thermoelastic material rotation required"));
        }
        Ok(())
    }
}

/// One source-pinned thermal potential for the existing mechanical and heat
/// owners. The original isotropic implementation and arithmetic are retained.
#[derive(Clone, Debug, PartialEq)]
pub enum AffineThermoelastic {
    Isotropic(AffineIsotropicThermoelastic),
    Tensor(AffineTensorThermoelastic),
}

impl AffineThermoelastic {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        card: &MaterialCard, elastic: &ResolvedMaterialStatePoint,
        expansion_reference: &ResolvedMaterialStatePoint,
        stiffness: ThermoelasticSymmetry, expansion: ThermoelasticSymmetry,
        temperature_range_k: [f64; 2], maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, ThermoelasticError> {
        if stiffness == ThermoelasticSymmetry::Isotropic && expansion == ThermoelasticSymmetry::Isotropic {
            Ok(Self::Isotropic(AffineIsotropicThermoelastic::new(card, elastic,
                expansion_reference, temperature_range_k, maximum_absolute_mandel_strain)?))
        } else {
            Ok(Self::Tensor(AffineTensorThermoelastic::new(card, elastic,
                expansion_reference, stiffness, expansion, temperature_range_k,
                maximum_absolute_mandel_strain)?))
        }
    }
    pub fn temperature_range_k(&self) -> [f64; 2] {
        match self { Self::Isotropic(x) => x.temperature_range_k(), Self::Tensor(x) => x.temperature_range_k }
    }
    pub fn maximum_absolute_mandel_strain(&self) -> f64 {
        match self { Self::Isotropic(x) => x.maximum_absolute_mandel_strain(), Self::Tensor(x) => x.maximum_absolute_mandel_strain }
    }
    pub fn at_temperature(&self, t: f64) -> Result<ResolvedThermoelastic, ThermoelasticError> {
        match self { Self::Isotropic(x) => x.at_temperature(t), Self::Tensor(x) => x.at_temperature(t) }
    }
    pub fn at_materials(&self, elastic: &ResolvedMaterialStatePoint, expansion: &ResolvedMaterialStatePoint)
        -> Result<ResolvedThermoelastic, ThermoelasticError> {
        match self { Self::Isotropic(x) => x.at_materials(elastic, expansion), Self::Tensor(x) => x.at_materials(elastic, expansion) }
    }
    pub fn minimum_heat_capacity_j_k(&self, volume_m3: f64) -> Result<f64, ThermoelasticError> {
        match self { Self::Isotropic(x) => x.minimum_heat_capacity_j_k(volume_m3), Self::Tensor(x) => x.minimum_heat_capacity_j_k(volume_m3) }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Coefficient {
    name: &'static str, dims: Dims, reference: f64, slope: f64, endpoints: [f64; 2],
}
impl Coefficient {
    fn new(card: &MaterialCard, source: &ResolvedMaterialStatePoint, name: &'static str,
        dims: Dims, range: [f64; 2], positive: bool) -> Result<Self, ThermoelasticError> {
        let reference = value(source, name, dims)?;
        let (_, slope, endpoints) = affine_source_path(card, source, name, "T", range, positive).map_err(field)?;
        if slope == 0.0 && endpoints[0] != endpoints[1] {
            return Err(ThermoelasticError::Invalid("nonconstant tensor thermal derivative underflows"));
        }
        Ok(Self { name, dims, reference, slope, endpoints })
    }
    fn maximum_abs(&self) -> f64 { self.endpoints[0].abs().max(self.endpoints[1].abs()) }
    fn minimum(&self) -> f64 { self.endpoints[0].min(self.endpoints[1]) }
    fn maximum(&self) -> f64 { self.endpoints[0].max(self.endpoints[1]) }
    fn at(&self, source: &ResolvedMaterialStatePoint, delta: f64) -> Result<f64, ThermoelasticError> {
        let actual = value(source, self.name, self.dims)?;
        let expected = finite(self.slope.mul_add(delta, self.reference))?;
        if (actual - expected).abs() > 512.0 * f64::EPSILON * self.maximum_abs().max(actual.abs()).max(expected.abs()) {
            return Err(ThermoelasticError::Invalid("tensor thermoelastic source left its admitted affine path"));
        }
        Ok(actual)
    }
}

/// A bounded smooth tensor source path. Stiffness is the native compliance
/// inverse, with fully independent affine engineering constants; expansion can
/// have different fixed axes. Current density is never queried or synthesized.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineTensorThermoelastic {
    card: MaterialCard,
    elastic_reference: ResolvedMaterialStatePoint,
    expansion_reference: ResolvedMaterialStatePoint,
    stiffness: ThermoelasticSymmetry,
    expansion: ThermoelasticSymmetry,
    temperature_range_k: [f64; 2],
    reference_temperature_k: f64,
    maximum_absolute_mandel_strain: f64,
    elastic_coefficients: Vec<Coefficient>,
    expansion_coefficients: Vec<Coefficient>,
    curvature_bound_j_m3_k2: f64,
}

impl AffineTensorThermoelastic {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        card: &MaterialCard, elastic: &ResolvedMaterialStatePoint,
        expansion_reference: &ResolvedMaterialStatePoint,
        stiffness: ThermoelasticSymmetry, expansion: ThermoelasticSymmetry,
        temperature_range_k: [f64; 2], maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, ThermoelasticError> {
        if !maximum_absolute_mandel_strain.is_finite() || maximum_absolute_mandel_strain <= 0.0
            || maximum_absolute_mandel_strain >= 1.0 {
            return Err(ThermoelasticError::Invalid("a finite small-strain component bound in (0,1) is required"));
        }
        stiffness.admit()?; expansion.admit()?;
        let reference = elastic.at_query_point(card, &expansion_reference.query())?;
        admit_source_temperature_refresh(elastic, &reference).map_err(field)?;
        let reference_temperature_k = temperature(expansion_reference)?;
        let elastic_schema: Vec<_> = match stiffness {
            ThermoelasticSymmetry::Isotropic => vec![
                (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, true), (POISSON_RATIO_PROPERTY, Dims::NONE, false)],
            ThermoelasticSymmetry::OrientedOrthotropic { .. } => ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES
                .map(|name| (name, Pressure::DIMS, true)).into_iter()
                .chain(ORTHOTROPIC_POISSON_RATIO_PROPERTIES.map(|name| (name, Dims::NONE, false)))
                .chain(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES.map(|name| (name, Pressure::DIMS, true))).collect(),
        };
        let expansion_names: Vec<_> = match expansion {
            ThermoelasticSymmetry::Isotropic => vec![LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY],
            ThermoelasticSymmetry::OrientedOrthotropic { .. } => ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES.to_vec(),
        };
        let elastic_coefficients = elastic_schema.into_iter().map(|(name, dims, positive)|
            Coefficient::new(card, &reference, name, dims, temperature_range_k, positive)).collect::<Result<Vec<_>, _>>()?;
        let expansion_coefficients = expansion_names.into_iter().map(|name|
            Coefficient::new(card, expansion_reference, name, Dims([0,0,0,-1,0,0]), temperature_range_k, false))
            .collect::<Result<Vec<_>, _>>()?;
        let mut law = Self { card: card.clone(), elastic_reference: reference,
            expansion_reference: expansion_reference.clone(), stiffness, expansion,
            temperature_range_k, reference_temperature_k, maximum_absolute_mandel_strain,
            elastic_coefficients, expansion_coefficients, curvature_bound_j_m3_k2: 0.0 };
        law.curvature_bound_j_m3_k2 = law.curvature_bound()?;
        for t in temperature_range_k { law.at_temperature(t)?; }
        Ok(law)
    }

    pub fn at_temperature(&self, t: f64) -> Result<ResolvedThermoelastic, ThermoelasticError> {
        if !t.is_finite() || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&t) {
            return Err(ThermoelasticError::Invalid("temperature leaves the admitted tensor thermoelastic source interval"));
        }
        let point = self.elastic_reference.query().with("T", t)
            .map_err(|_| ThermoelasticError::Invalid("tensor thermal query is unrepresentable"))?;
        let elastic = self.elastic_reference.at_query_point(&self.card, &point)?;
        let expansion = self.expansion_reference.at_query_point(&self.card, &point)?;
        let delta = t - self.reference_temperature_k;
        let coefficients = self.elastic_coefficients.iter().map(|x| x.at(&elastic, delta)).collect::<Result<Vec<_>, _>>()?;
        let slopes: Vec<_> = self.elastic_coefficients.iter().map(|x| x.slope).collect();
        let (mut stiffness, mut first, mut second) = principal_stiffness(self.stiffness,
            &coefficients, &slopes, self.maximum_absolute_mandel_strain)?;
        let q = self.stiffness.rotation();
        // The original solid owner also applies its shared Bond rotation for
        // an identity-oriented orthotropic source. Preserve its exact tensor.
        if matches!(self.stiffness, ThermoelasticSymmetry::OrientedOrthotropic { .. }) {
            stiffness = rotate_mandel(stiffness, q);
            first = rotate_mandel(first, q);
            second = rotate_mandel(second, q);
        }
        let alpha = self.expansion_coefficients.iter().map(|x| x.at(&expansion, delta)).collect::<Result<Vec<_>, _>>()?;
        let selected = MaterialPropertySelection::PinnedByProperty(self.expansion_coefficients.iter().map(|x| {
            let claim = self.expansion_reference.property(x.name).expect("complete expansion schema admitted").answer().receipt.selected;
            (x.name.to_owned(), claim)
        }).collect());
        let (free, instantaneous, derivative) = match self.expansion {
            ThermoelasticSymmetry::Isotropic => {
                let path = integrate_isotropic_thermal_expansion(&self.card, &self.expansion_reference.query(), &point, selected)?;
                ([path.free_linear_strain(); 3], [alpha[0]; 3], [self.expansion_coefficients[0].slope; 3])
            }
            ThermoelasticSymmetry::OrientedOrthotropic { .. } => {
                let path = integrate_orthotropic_thermal_expansion(&self.card, &self.expansion_reference.query(), &point, selected)?;
                (*path.free_principal_strains(), [alpha[0], alpha[1], alpha[2]],
                    core::array::from_fn(|i| self.expansion_coefficients[i].slope))
            }
        };
        let response = TensorThermalResponse {
            stiffness_mandel_pa: stiffness, stiffness_derivative_mandel_pa_k: first,
            stiffness_second_derivative_mandel_pa_k2: second,
            free_strain_mandel: rotate_diagonal(free, self.expansion.rotation()),
            expansion_mandel_per_k: rotate_diagonal(instantaneous, self.expansion.rotation()),
            expansion_derivative_mandel_per_k2: rotate_diagonal(derivative, self.expansion.rotation()),
        };
        if response.values().any(|x| !x.is_finite()) {
            return Err(ThermoelasticError::Invalid("oriented tensor thermal response is unrepresentable"));
        }
        Ok(ResolvedThermoelastic { elastic, expansion, temperature_k: t,
            maximum_absolute_mandel_strain: self.maximum_absolute_mandel_strain,
            response: ResolvedThermoelasticLaw::Tensor(response) })
    }

    pub fn at_materials(&self, elastic: &ResolvedMaterialStatePoint, expansion: &ResolvedMaterialStatePoint)
        -> Result<ResolvedThermoelastic, ThermoelasticError> {
        admit_source_temperature_refresh(&self.elastic_reference, elastic).map_err(field)?;
        admit_source_temperature_refresh(&self.expansion_reference, expansion).map_err(field)?;
        let t = temperature(elastic)?;
        if t.to_bits() != temperature(expansion)?.to_bits() || elastic.query() != expansion.query() {
            return Err(ThermoelasticError::Invalid("tensor elasticity and expansion must share the actual physical query"));
        }
        let at = self.at_temperature(t)?;
        for (expected, actual) in [(&at.elastic, elastic), (&at.expansion, expansion)] {
            if expected.properties().len() != actual.properties().len() || expected.properties().iter().zip(actual.properties())
                .any(|(a,b)| a.value_si().to_bits() != b.value_si().to_bits()) {
                return Err(ThermoelasticError::Invalid("native tensor mechanics and caloric potential selected different source values"));
            }
        }
        Ok(at)
    }

    pub fn minimum_heat_capacity_j_k(&self, volume_m3: f64) -> Result<f64, ThermoelasticError> {
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(ThermoelasticError::Invalid("positive finite reference volume is required"));
        }
        Ok(-upper_mul(upper_mul(self.temperature_range_k[1], volume_m3)?, self.curvature_bound_j_m3_k2)?)
    }

    fn curvature_bound(&self) -> Result<f64, ThermoelasticError> {
        let (c, first, second) = stiffness_bounds(self.stiffness, &self.elastic_coefficients)?;
        let stiffness_rotation = rotation_norm_bound(self.stiffness.rotation())?;
        let stiffness_scale = upper_mul(stiffness_rotation, stiffness_rotation)?;
        let (c, first, second) = (upper_mul(c, stiffness_scale)?,
            upper_mul(first, stiffness_scale)?, upper_mul(second, stiffness_scale)?);
        let mut alpha = 0.0; let mut alpha_first = 0.0;
        for i in 0..3 {
            let a = &self.expansion_coefficients[if self.expansion_coefficients.len() == 1 { 0 } else { i }];
            alpha = upper_add(alpha, a.maximum_abs())?;
            alpha_first = upper_add(alpha_first, a.slope.abs())?;
        }
        let expansion_rotation = rotation_norm_bound(self.expansion.rotation())?;
        alpha = upper_mul(alpha, expansion_rotation)?;
        alpha_first = upper_mul(alpha_first, expansion_rotation)?;
        // ||e||_2 <= sqrt(6)*component_limit <= 3*component_limit. The
        // response checks BOTH total and elastic strain. Actual finite Bond
        // transforms are bounded too; approximate orthogonality is not used
        // as an exact norm identity in this admission.
        let strain = upper_mul(3.0, self.maximum_absolute_mandel_strain)?;
        let curvature = upper_div(upper_mul(second, upper_mul(strain, strain)?)?, 2.0)?;
        let mixed = upper_mul(2.0, upper_mul(alpha, upper_mul(first, strain)?)?)?;
        let expansion = upper_mul(alpha, upper_mul(c, alpha)?)?;
        let acceleration = upper_mul(alpha_first, upper_mul(c, strain)?)?;
        upper_add(upper_add(curvature, mixed)?, upper_add(expansion, acceleration)?)
    }
}

fn rotate_diagonal(principal: [f64; 3], q: [[f64; 3]; 3]) -> [f64; 6] {
    let source = [principal[0], principal[1], principal[2], 0.0, 0.0, 0.0];
    if q == IDENTITY { return source; }
    let transform = mandel_rotation(q);
    core::array::from_fn(|i| (0..6).fold(0.0, |sum, j| transform[i][j].mul_add(source[j], sum)))
}

fn rotation_norm_bound(q: [[f64; 3]; 3]) -> Result<f64, ThermoelasticError> {
    let transform = mandel_rotation(q);
    let mut bound: f64 = 1.0;
    for i in 0..6 {
        let mut row = 0.0; let mut column = 0.0;
        for j in 0..6 {
            row = upper_add(row, transform[i][j].abs())?;
            column = upper_add(column, transform[j][i].abs())?;
        }
        bound = bound.max(row).max(column);
    }
    Ok(bound)
}

/// Exact first/second temperature derivatives of the native compliance inverse.
fn principal_stiffness(symmetry: ThermoelasticSymmetry, values: &[f64], slopes: &[f64], limit: f64)
    -> Result<(Tangent6, Tangent6, Tangent6), ThermoelasticError> {
    let (e, nu, et, nut, c) = match symmetry {
        ThermoelasticSymmetry::Isotropic => {
            let law = IsotropicElastic::new(values[0], values[1], limit)?;
            ([values[0]; 3], [values[1]; 3], [slopes[0]; 3], [slopes[1]; 3], law.stiffness())
        }
        ThermoelasticSymmetry::OrientedOrthotropic { .. } => {
            let e = [values[0], values[1], values[2]];
            let nu = [values[3], values[4], values[5]];
            let law = OrthotropicElastic::new(e, nu, [values[6], values[7], values[8]], limit)?;
            (e, nu, [slopes[0], slopes[1], slopes[2]], [slopes[3], slopes[4], slopes[5]], law.stiffness())
        }
    };
    let mut st = [[0.0; 3]; 3]; let mut stt = [[0.0; 3]; 3];
    for i in 0..3 {
        st[i][i] = finite(-(et[i] / e[i]) / e[i])?;
        stt[i][i] = finite(-2.0 * (et[i] / e[i]) * st[i][i])?;
    }
    for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
        let slope = finite((-nut[k] + nu[k] * (et[i] / e[i])) / e[i])?;
        let curvature = finite(-2.0 * (et[i] / e[i]) * slope)?;
        st[i][j] = slope; st[j][i] = slope;
        stt[i][j] = curvature; stt[j][i] = curvature;
    }
    let mut first = [[0.0; 6]; 6]; let mut second = [[0.0; 6]; 6];
    for i in 0..3 { for j in i..3 {
        first[i][j] = finite(-(0..3).map(|a| (0..3).map(|b| c[i][a] * st[a][b] * c[b][j]).sum::<f64>()).sum::<f64>())?;
        first[j][i] = first[i][j];
    }}
    for i in 0..3 { for j in i..3 {
        second[i][j] = finite(-(0..3).map(|a| (0..3).map(|b|
            first[i][a] * st[a][b] * c[b][j] + c[i][a] * stt[a][b] * c[b][j]
                + c[i][a] * st[a][b] * first[b][j]).sum::<f64>()).sum::<f64>())?;
        second[j][i] = second[i][j];
    }}
    for i in 3..6 {
        match symmetry {
            ThermoelasticSymmetry::Isotropic => {
                let q = 1.0 + values[1];
                first[i][i] = finite((slopes[0] - values[0] * slopes[1] / q) / q)?;
                second[i][i] = finite(-2.0 * slopes[1] * first[i][i] / q)?;
            }
            ThermoelasticSymmetry::OrientedOrthotropic { .. } => { first[i][i] = finite(2.0 * slopes[i + 3])?; }
        }
    }
    if c.iter().flatten().any(|x| !x.is_finite()) {
        return Err(ThermoelasticError::Invalid("native orthotropic stiffness is unrepresentable"));
    }
    Ok((c, first, second))
}

fn positive_lower(x: f64) -> Result<f64, ThermoelasticError> {
    let value = x.next_down();
    if !value.is_finite() || value <= 0.0 {
        return Err(ThermoelasticError::Invalid("whole thermal interval cannot establish positive finite orthotropic compliance"));
    }
    Ok(value)
}
fn row_norm(a: &[[f64; 3]; 3]) -> Result<f64, ThermoelasticError> {
    let mut norm: f64 = 0.0;
    for row in a { let mut sum = 0.0; for &v in row { sum = upper_add(sum, v)?; } norm = norm.max(sum); }
    Ok(norm)
}

// All entries below are nonnegative outward upper bounds. Let R be the
// compliance normalized to unit diagonal. Its determinant is
// 1-nu12² E2/E1-nu13² E3/E1-nu23² E3/E2-2 nu12 nu13 nu23 E3/E1.
// Bounding every subtracted term by magnitude gives a sufficient Sylvester
// admission over the ENTIRE affine source interval, including signed nu.
fn stiffness_bounds(symmetry: ThermoelasticSymmetry, coefficients: &[Coefficient])
    -> Result<(f64, f64, f64), ThermoelasticError> {
    let isotropic = symmetry == ThermoelasticSymmetry::Isotropic;
    let e: [&Coefficient; 3] = core::array::from_fn(|i| &coefficients[if isotropic { 0 } else { i }]);
    let nu: [&Coefficient; 3] = core::array::from_fn(|i| &coefficients[if isotropic { 1 } else { i + 3 }]);
    let emin = e.map(|x| positive_lower(x.minimum())).into_iter().collect::<Result<Vec<_>, _>>()?;
    let emax = e.map(|x| x.maximum().next_up());
    let nmax = nu.map(Coefficient::maximum_abs);
    let mut defect = 0.0;
    for (k, (i, j)) in [(0,1), (0,2), (1,2)].into_iter().enumerate() {
        defect = upper_add(defect, upper_mul(upper_mul(nmax[k], nmax[k])?, upper_div(emax[j], emin[i])?)?)?;
    }
    let product = upper_mul(nmax[0], upper_mul(nmax[1], nmax[2])?)?;
    defect = upper_add(defect, upper_mul(2.0, upper_mul(product, upper_div(emax[2], emin[0])?)?)?)?;
    let normalized_det = positive_lower(1.0 - defect)?;
    let mut st = [[0.0; 3]; 3]; let mut stt = st;
    for i in 0..3 {
        let relative = upper_div(e[i].slope.abs(), emin[i])?;
        st[i][i] = upper_div(relative, emin[i])?;
        stt[i][i] = upper_mul(2.0, upper_mul(relative, st[i][i])?)?;
    }
    for (k, (i, j)) in [(0,1), (0,2), (1,2)].into_iter().enumerate() {
        let relative = upper_div(e[i].slope.abs(), emin[i])?;
        st[i][j] = upper_div(upper_add(nu[k].slope.abs(), upper_mul(nmax[k], relative)?)?, emin[i])?;
        st[j][i] = st[i][j];
        stt[i][j] = upper_mul(2.0, upper_mul(relative, st[i][j])?)?; stt[j][i] = stt[i][j];
    }
    // Algebraically scaled cofactors avoid cubed dimensional moduli and keep
    // strong principal anisotropy out of unrelated directions' bounds.
    let mut inverse = [[0.0; 3]; 3];
    for i in 0..3 { inverse[i][i] = upper_div(emax[i], normalized_det)?; }
    inverse[0][1] = upper_div(upper_add(upper_mul(nmax[0],emax[1])?,
        upper_mul(upper_mul(nmax[1],nmax[2])?,emax[2])?)?,normalized_det)?;
    inverse[0][2] = upper_div(upper_mul(emax[2],upper_add(nmax[1],upper_mul(nmax[0],nmax[2])?)?)?,normalized_det)?;
    inverse[1][2] = upper_div(upper_mul(emax[2],upper_add(nmax[2],
        upper_mul(upper_mul(nmax[0],nmax[1])?,upper_div(emax[1],emin[0])?)?)?)?,normalized_det)?;
    for i in 0..3 { for j in 0..i { inverse[i][j] = inverse[j][i]; } }
    let first_matrix = bound_product(&bound_product(&inverse,&st)?,&inverse)?;
    let left = bound_product(&bound_product(&first_matrix,&st)?,&inverse)?;
    let middle = bound_product(&bound_product(&inverse,&stt)?,&inverse)?;
    let right = bound_product(&bound_product(&inverse,&st)?,&first_matrix)?;
    let mut second_matrix = [[0.0;3];3];
    for i in 0..3 {for j in 0..3 {
        second_matrix[i][j]=upper_add(upper_add(left[i][j],middle[i][j])?,right[i][j])?;
    }}
    let mut c = row_norm(&inverse)?;
    let mut first = row_norm(&first_matrix)?;
    let mut second = row_norm(&second_matrix)?;
    if isotropic {
        let q = positive_lower(1.0 + coefficients[1].minimum())?;
        positive_lower(1.0 - 2.0 * coefficients[1].maximum())?;
        let shear = upper_div(emax[0], q)?;
        let shear_first = upper_div(upper_add(coefficients[0].slope.abs(),
            upper_mul(shear, coefficients[1].slope.abs())?)?, q)?;
        let shear_second = upper_div(upper_mul(2.0, upper_mul(coefficients[1].slope.abs(), shear_first)?)?, q)?;
        c = c.max(shear); first = first.max(shear_first); second = second.max(shear_second);
    } else {
        for shear in &coefficients[6..9] {
            c = c.max(upper_mul(2.0, shear.maximum())?);
            first = first.max(upper_mul(2.0, shear.slope.abs())?);
        }
    }
    Ok((c,first,second))
}

fn bound_product(a:&[[f64;3];3],b:&[[f64;3];3])->Result<[[f64;3];3],ThermoelasticError> {
    let mut result=[[0.0;3];3];
    for i in 0..3 {for j in 0..3 {for k in 0..3 {
        result[i][j]=upper_add(result[i][j],upper_mul(a[i][k],b[k][j])?)?;
    }}}
    Ok(result)
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct TensorThermalResponse {
    pub(super) stiffness_mandel_pa: Tangent6,
    pub(super) stiffness_derivative_mandel_pa_k: Tangent6,
    pub(super) stiffness_second_derivative_mandel_pa_k2: Tangent6,
    pub(super) free_strain_mandel: [f64; 6],
    pub(super) expansion_mandel_per_k: [f64; 6],
    pub(super) expansion_derivative_mandel_per_k2: [f64; 6],
}

impl TensorThermalResponse {
    fn values(&self) -> impl Iterator<Item=f64> {
        self.stiffness_mandel_pa.into_iter().flatten()
            .chain(self.stiffness_derivative_mandel_pa_k.into_iter().flatten())
            .chain(self.stiffness_second_derivative_mandel_pa_k2.into_iter().flatten())
            .chain(self.free_strain_mandel).chain(self.expansion_mandel_per_k)
            .chain(self.expansion_derivative_mandel_per_k2)
    }
    pub(super) fn state(&self, t: f64, total: [f64; 6], volume: f64, limit: f64)
        -> Result<ThermoelasticState, ThermoelasticError> {
        if !volume.is_finite() || volume <= 0.0 || total.iter().any(|x| !x.is_finite() || x.abs() > limit) {
            return Err(ThermoelasticError::Invalid("strain/volume leaves the admitted tensor thermoelastic envelope"));
        }
        let e: [f64; 6] = core::array::from_fn(|i| total[i] - self.free_strain_mandel[i]);
        if e.iter().any(|x| !x.is_finite() || x.abs() > limit) {
            return Err(ThermoelasticError::Invalid("elastic strain leaves the admitted tensor small-strain envelope"));
        }
        let apply = |c: &Tangent6, x: [f64; 6]| -> [f64; 6] {
            core::array::from_fn(|i| (0..6).fold(0.0, |sum,j| c[i][j].mul_add(x[j],sum)))
        };
        let dot = |a: [f64;6], b: [f64;6]| -> f64 { a.into_iter().zip(b).map(|(x,y)| x*y).sum() };
        let stress = apply(&self.stiffness_mandel_pa,e);
        let first_stress = apply(&self.stiffness_derivative_mandel_pa_k,e);
        let thermal_stress = apply(&self.stiffness_mandel_pa,self.expansion_mandel_per_k);
        let beta = core::array::from_fn(|i| first_stress[i] - thermal_stress[i]);
        let free_energy_j = finite(0.5 * volume * dot(e,stress))?;
        let entropy_j_k = finite(volume * (dot(self.expansion_mandel_per_k,stress) - 0.5*dot(e,first_stress)))?;
        let internal_energy_j = finite(t.mul_add(entropy_j_k,free_energy_j))?;
        let curvature = 0.5*dot(e,apply(&self.stiffness_second_derivative_mandel_pa_k2,e))
            - 2.0*dot(self.expansion_mandel_per_k,first_stress)
            + dot(self.expansion_mandel_per_k,thermal_stress)
            - dot(self.expansion_derivative_mandel_per_k2,stress);
        let heat_capacity_j_k = finite(-t * volume * curvature)?;
        if free_energy_j < 0.0 || stress.iter().chain(beta.iter()).any(|x| !x.is_finite()) {
            return Err(ThermoelasticError::Invalid("tensor thermoelastic response is unrepresentable"));
        }
        Ok(ThermoelasticState { elastic_strain_mandel: e, stress_mandel_pa: stress,
            stress_temperature_derivative_mandel_pa_k: beta, free_energy_j, entropy_j_k,
            internal_energy_j, heat_capacity_j_k })
    }
}

#[cfg(test)]
mod tests;
