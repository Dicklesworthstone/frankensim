//! Thermal derivatives of the SAME longitudinal strain/charge potential.
//!
//! A(T)=[[s_E,d],[d,epsilon_T]] maps [stress,E] to x=[strain,D]. With fixed
//! axes and reference volume, psi=V*x^T A^-1*x/2. For an affine, positive
//! definite A, S=V*y^T A'*y/2 and C_x=-T*V*(A'*y)^T A^-1*(A'*y), where
//! y=A^-1*x. Signed d and d' are retained; diagonal-only derivative formulas
//! would discard real electromechanical thermal exchange.
//!
//! This uses the original small-signal constitutive response and the common
//! affine-source/capacity admission. It adds no pyroelectric offset, thermal
//! expansion, depoling, hysteresis, phase transition or mechanical integrator.

use super::*;
use crate::field_storage::{
    FieldStorageError, admit_source_temperature_refresh,
    thermal::{FieldThermalState, affine_source_path, capacity_bound},
};

#[derive(Clone, Debug, PartialEq)]
pub struct ThermalLongitudinalPiezoelectric {
    card: MaterialCard,
    original: LongitudinalPiezoelectric,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    reference_temperature_k: f64,
    slopes: [f64; 3],
    // A=S G S with S=diag(sqrt(s_ref),sqrt(epsilon_ref)). G is dimensionless;
    // its endpoint eigenvalue lower bounds every affine interior G by convexity.
    scales: [f64; 2],
    minimum_eigenvalue: f64,
    derivative_norm_bound_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalPiezoelectric {
    material: LongitudinalPiezoelectric,
    temperature_k: f64,
    slopes: [f64; 3],
}

fn field_error(error: FieldStorageError) -> PiezoelectricError {
    match error {
        FieldStorageError::Material(error) => PiezoelectricError::Material(error),
        FieldStorageError::Invalid(message) => PiezoelectricError::Invalid(message),
    }
}
fn finite(value: f64) -> Result<f64, PiezoelectricError> {
    if value.is_finite() { Ok(value) }
    else { Err(PiezoelectricError::Invalid("piezoelectric thermal arithmetic is unrepresentable")) }
}
fn coefficients(material: &LongitudinalPiezoelectric) -> [f64; 3] {
    [material.compliance_at_fixed_field_pa_inv(), material.strain_charge_coefficient_m_v(),
        material.permittivity_at_fixed_stress_f_m()]
}

impl ThermalLongitudinalPiezoelectric {
    /// Pin the complete original source and affine s_E(T), signed d(T), and
    /// epsilon_T(T). Positive coupled stability at both endpoints admits the
    /// whole affine interval. Reference density must be constant throughout.
    pub fn new(
        card: &MaterialCard,
        original: &LongitudinalPiezoelectric,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, PiezoelectricError> {
        let mut slopes = [0.0; 3];
        let mut reference_temperature_k = 0.0;
        for (i, name) in [COMPLIANCE_33_E_PROPERTY, STRAIN_CHARGE_33_PROPERTY,
            PERMITTIVITY_33_T_PROPERTY].into_iter().enumerate() {
            let (reference, slope, _) = affine_source_path(card, original.material(), name,
                temperature_axis, temperature_range_k, i != 1).map_err(field_error)?;
            reference_temperature_k = reference;
            slopes[i] = slope;
        }
        let (_, density_slope, _) = affine_source_path(card, original.material(), DENSITY_PROPERTY,
            temperature_axis, temperature_range_k, true).map_err(field_error)?;
        if density_slope != 0.0 {
            return Err(PiezoelectricError::Invalid("fixed-reference piezoelectric thermal motion requires constant source density"));
        }
        let scales = [fs_math::det::sqrt(original.compliance_at_fixed_field_pa_inv()),
            fs_math::det::sqrt(original.permittivity_at_fixed_stress_f_m())];
        let cross = finite(slopes[1] / scales[0] / scales[1])?;
        let derivative_norm_bound_per_k = finite(
            (slopes[0] / scales[0] / scales[0]).abs() + cross.abs(),
        )?.max(finite(cross.abs() + (slopes[2] / scales[1] / scales[1]).abs())?);
        if slopes.iter().any(|&s| s != 0.0) && derivative_norm_bound_per_k == 0.0 {
            return Err(PiezoelectricError::Invalid("coupled thermal derivative bound underflows"));
        }
        let mut law = Self { card: card.clone(), original: original.clone(),
            temperature_axis: temperature_axis.to_owned(), temperature_range_k,
            reference_temperature_k, slopes, scales, minimum_eigenvalue: f64::INFINITY,
            derivative_norm_bound_per_k };
        for temperature in temperature_range_k {
            let material = law.at_temperature(temperature)?.material;
            let [s,d,e] = coefficients(&material);
            let (a,b,c) = (s/scales[0]/scales[0], d/scales[0]/scales[1], e/scales[1]/scales[1]);
            let determinant = finite(a.mul_add(c,-b*b))?;
            let largest = finite(0.5*(a+c+(a-c).hypot(2.0*b)))?;
            let smallest = finite(determinant/largest)? * (1.0-4096.0*f64::EPSILON);
            if smallest <= 0.0 {
                return Err(PiezoelectricError::Invalid("complete affine piezoelectric block must remain representably positive definite"));
            }
            law.minimum_eigenvalue = law.minimum_eigenvalue.min(smallest);
        }
        Ok(law)
    }
    pub fn card(&self) -> &MaterialCard { &self.card }
    pub fn original(&self) -> &LongitudinalPiezoelectric { &self.original }
    pub fn temperature_axis(&self) -> &str { &self.temperature_axis }
    pub fn temperature_range_k(&self) -> [f64; 2] { self.temperature_range_k }

    /// Re-query the original full bundle, without dropping any other source role.
    pub fn at_temperature(&self, temperature_k: f64) -> Result<ResolvedThermalPiezoelectric, PiezoelectricError> {
        if !temperature_k.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&temperature_k) {
            return Err(PiezoelectricError::Invalid("piezoelectric temperature is outside its supported source interval"));
        }
        let point = self.original.material().query().with(&self.temperature_axis,temperature_k)
            .map_err(|_| PiezoelectricError::Invalid("piezoelectric thermal query is unrepresentable"))?;
        let material = LongitudinalPiezoelectric::from_material_state(
            self.original.material().at_query_point(&self.card,&point)?,
        )?;
        for ((original, actual), slope) in coefficients(&self.original).into_iter()
            .zip(coefficients(&material)).zip(self.slopes) {
            let expected = finite(slope.mul_add(temperature_k-self.reference_temperature_k,original))?;
            // A signed cross coefficient may pass through zero. Compare with
            // its admitted path scale, without clamping the actual source value.
            let span = (self.temperature_range_k[0]-self.reference_temperature_k).abs()
                .max((self.temperature_range_k[1]-self.reference_temperature_k).abs());
            let scale = original.abs().max(finite(slope.abs()*span)?)
                .max(actual.abs()).max(expected.abs());
            if (actual-expected).abs() > 512.0*f64::EPSILON*scale {
                return Err(PiezoelectricError::Invalid("piezoelectric source left its admitted affine path"));
            }
        }
        if material.density_kg_m3().to_bits() != self.original.density_kg_m3().to_bits() {
            return Err(PiezoelectricError::Invalid("piezoelectric reference density changed"));
        }
        Ok(ResolvedThermalPiezoelectric { material, temperature_k, slopes: self.slopes })
    }

    pub fn at_material(&self, material: &LongitudinalPiezoelectric)
        -> Result<ResolvedThermalPiezoelectric, PiezoelectricError> {
        admit_source_temperature_refresh(self.original.material(),material.material()).map_err(field_error)?;
        let temperature = material.material().query_point().iter()
            .find_map(|(name,value)| (name == &self.temperature_axis).then_some(*value))
            .ok_or(PiezoelectricError::Invalid("piezoelectric temperature coordinate is missing"))?;
        let at = self.at_temperature(temperature)?;
        if at.material.material().properties().iter().zip(material.material().properties())
            .any(|(a,b)| a.value_si().to_bits() != b.value_si().to_bits()) {
            return Err(PiezoelectricError::Invalid("piezoelectric dynamics and caloric law selected different source values"));
        }
        Ok(at)
    }

    /// Conservative full-interval lower bound on fixed-(strain,D) field heat
    /// capacity. Caller adds the actual whole-rod caloric capacity and enforces
    /// both coordinate limits. This is a model bound, not a numerical certificate.
    pub fn minimum_heat_capacity_j_k(&self, maximum_abs_strain: f64,
        maximum_abs_displacement_c_m2: f64, volume_m3: f64) -> Result<f64, PiezoelectricError> {
        if [maximum_abs_strain,maximum_abs_displacement_c_m2].iter()
            .any(|v| !v.is_finite() || *v < 0.0) || !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(PiezoelectricError::Invalid("finite nonnegative piezoelectric coordinate bounds and positive volume required"));
        }
        let norm = finite((maximum_abs_strain/self.scales[0])
            .hypot(maximum_abs_displacement_c_m2/self.scales[1]))?;
        let worst = finite(0.5*volume_m3*norm*(norm/self.minimum_eigenvalue))?;
        let ratio = finite(self.derivative_norm_bound_per_k/self.minimum_eigenvalue)?;
        if (maximum_abs_strain != 0.0 || maximum_abs_displacement_c_m2 != 0.0) && worst == 0.0
            || self.derivative_norm_bound_per_k != 0.0 && ratio == 0.0 {
            return Err(PiezoelectricError::Invalid("piezoelectric thermal capacity bound underflows"));
        }
        capacity_bound(worst,self.temperature_range_k[1],ratio).map_err(field_error)
    }
}

impl ResolvedThermalPiezoelectric {
    pub fn material(&self) -> &LongitudinalPiezoelectric { &self.material }
    pub fn temperature_k(&self) -> f64 { self.temperature_k }

    /// Exact integrated potentials of the original coupled compliance law.
    pub fn state(&self, strain: f64, displacement_c_m2: f64, volume_m3: f64)
        -> Result<FieldThermalState, PiezoelectricError> {
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(PiezoelectricError::Invalid("piezoelectric thermal volume must be finite positive"));
        }
        let y = self.material.response(strain,displacement_c_m2)?;
        let derivative = [finite(self.slopes[0].mul_add(y.stress_pa,self.slopes[1]*y.electric_field_v_m))?,
            finite(self.slopes[1].mul_add(y.stress_pa,self.slopes[2]*y.electric_field_v_m))?];
        let free_energy_j = finite(volume_m3*y.stored_energy_j_m3)?;
        let entropy_j_k = finite(0.5*volume_m3*y.stress_pa.mul_add(derivative[0],y.electric_field_v_m*derivative[1]))?;
        // The same native inverse block evaluates (A'y)^T A^-1(A'y)/2.
        let quadratic = self.material.response(derivative[0],derivative[1])?.stored_energy_j_m3;
        let heat_capacity_j_k = finite(-2.0*self.temperature_k*volume_m3*quadratic)?;
        let internal_energy_j = finite(self.temperature_k.mul_add(entropy_j_k,free_energy_j))?;
        if y.stored_energy_j_m3 != 0.0 && free_energy_j == 0.0 || quadratic != 0.0 && heat_capacity_j_k == 0.0 {
            return Err(PiezoelectricError::Invalid("piezoelectric thermal contribution underflows"));
        }
        Ok(FieldThermalState { free_energy_j,entropy_j_k,internal_energy_j,heat_capacity_j_k })
    }
}
