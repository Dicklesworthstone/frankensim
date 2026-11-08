//! Full source-bound three-dimensional piezoelectric thermal storage.
//!
//! At fixed Mandel strain S and Cartesian electric displacement D, the optional
//! sourced thermal expansion alpha and CLAMPED pyroelectric coefficient p give
//! A=alpha*(T-Tref), P=p*(T-Tref), x=S-A, and
//! psi/V = x:C:x/2 + (D-eS-P).epsilonS^-1.(D-eS-P)/2.
//! With affine ORIGINAL stress-charge C(T), e(T), epsilonS(T), constant alpha,p,
//! and w=e'S+p+epsilonS'E, entropy/V is
//! alpha.C.x-x.C'.x/2+E.(e'S+p)+E.epsilonS'.E/2, while capacity/V is
//! -T*(alpha.C.alpha-2*alpha.C'.x+w.epsilonS^-1.w).
//! Without explicit bias sources alpha=p=0, recovering the original law.
//! This is the complete local strain/charge potential, without longitudinal
//! projection, relaxed coordinates or an electrode assumption.
//!
//! The native tensor owner retains the source/frame and supplies both field
//! recovery and inverse dielectric factor actions. A response uses fixed-size
//! vectors, without source queries, heap buffers or fresh factorizations.
//! Temperature refresh re-queries the complete pinned source bundle.
//!
//! The capacity bound is analytic model admission with roundoff padding, not
//! interval arithmetic. The consuming thermal owner adds its actual caloric
//! capacity and enforces the supplied coordinate bounds. Field entropy/capacity
//! describe reversible exchange; they are not dissipative Joule heat.
//! The caloric chart is the stress-free, zero-field reference: psi=0 at S=A,E=0.
//! A measured clamped heat capacity requires conversion before use as that chart.
//!
//! Stress-charge conventions:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.043.html>.
//! Helmholtz entropy/internal-energy relations:
//! <https://web.mit.edu/beretta/www/AdvancedThermodynamics_2_43_Lecture6.htm>.

use super::{field_error, finite};
use crate::field_storage::{
    admit_source_temperature_refresh,
    thermal::{FieldThermalState, affine_source_path},
};
use crate::piezoelectric::{PiezoelectricError, tensor::TensorPiezoelectric};
use crate::state_point::{DENSITY_PROPERTY, ResolvedMaterialStatePoint};
use crate::tensor::{mandel_rotation, rotate_mandel};
use fs_matdb::MaterialCard;

mod bias;
pub use bias::PYROELECTRIC_CLAMPED_DIMS;
use bias::{ResolvedThermalBias, ThermalBiasSource};

type Result<T> = core::result::Result<T, PiezoelectricError>;
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

#[derive(Clone, Copy, Debug)]
struct TensorTemperatureSlopes {
    stiffness: [[f64; 6]; 6],
    coupling: [[f64; 6]; 3],
    permittivity: [[f64; 3]; 3],
}

/// Full affine stress-charge source trajectory in the ORIGINAL tensor's fixed
/// current frame. Density, reference volume and all nonthermal coordinates are
/// held. Explicit thermal expansion and clamped pyroelectric sources can be
/// attached with [`Self::with_thermal_bias`]. Hysteresis, depoling, phase change,
/// strain/field amplitude envelopes and measured-material accuracy are not
/// inferred. The spatial consumer owns fields, dynamics and caloric storage.
#[derive(Clone, Debug)]
pub struct ThermalTensorPiezoelectric {
    card: MaterialCard,
    original: TensorPiezoelectric,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    slopes: TensorTemperatureSlopes,
    minimum_permittivity_f_m: f64,
    maximum_coupling_norm_c_m2: f64,
    coupling_slope_norm: f64,
    permittivity_slope_norm: f64,
    bias: Option<ThermalBiasSource>,
}

/// Native full tensor and exact thermal derivatives at one pinned temperature.
#[derive(Clone, Debug)]
pub struct ResolvedThermalTensorPiezoelectric {
    material: TensorPiezoelectric,
    temperature_k: f64,
    slopes: TensorTemperatureSlopes,
    bias: Option<ResolvedThermalBias>,
}

fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (a, b)| a.mul_add(*b, sum))
}
fn apply<const R: usize, const C: usize>(matrix: &[[f64; C]; R], x: &[f64; C]) -> [f64; R] {
    core::array::from_fn(|i| dot(&matrix[i], x))
}
fn frobenius_norm<const R: usize, const C: usize>(matrix: &[[f64; C]; R]) -> Result<f64> {
    finite(matrix.iter().flatten().fold(0.0_f64, |value, x| value.hypot(*x)))
}

impl TensorTemperatureSlopes {
    fn rotate(&mut self, q: [[f64; 3]; 3]) -> Result<()> {
        // Same tensor ranks and exact symmetric triangle convention as the
        // native material. No diagonal truncation or shared eigenbasis.
        self.stiffness = rotate_mandel(self.stiffness, q);
        let m = mandel_rotation(q);
        let mut coupling = [[0.0; 6]; 3];
        for (i, row) in coupling.iter_mut().enumerate() {
            for (a, value) in row.iter_mut().enumerate() {
                *value = (0..3).fold(0.0, |sum, j| {
                    (0..6).fold(sum, |sum, b| {
                        (q[i][j] * self.coupling[j][b]).mul_add(m[a][b], sum)
                    })
                });
            }
        }
        let mut permittivity = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let value = (0..3).fold(0.0, |sum, a| {
                    (0..3).fold(sum, |sum, b| {
                        (q[i][a] * self.permittivity[a][b]).mul_add(q[j][b], sum)
                    })
                });
                permittivity[i][j] = value;
                permittivity[j][i] = value;
            }
        }
        if self.stiffness.iter().flatten().chain(coupling.iter().flatten())
            .chain(permittivity.iter().flatten()).any(|value| !value.is_finite()) {
            return Err(PiezoelectricError::Invalid("rotated full piezoelectric thermal derivative is unrepresentable"));
        }
        self.coupling = coupling;
        self.permittivity = permittivity;
        Ok(())
    }
}

fn source_tensor(
    original: &TensorPiezoelectric,
    source: ResolvedMaterialStatePoint,
) -> Result<TensorPiezoelectric> {
    let tensor = TensorPiezoelectric::from_material_state(source, original.source_frame_identity())?;
    if original.frame_identity() == original.source_frame_identity()
        && *original.source_to_frame_rotation() == IDENTITY {
        Ok(tensor)
    } else {
        tensor.rotated(original.frame_identity(), *original.source_to_frame_rotation())
    }
}

fn admit_same_coefficients<const R: usize, const C: usize>(
    original: &[[f64; C]; R],
    refreshed: &[[f64; C]; R],
) -> Result<()> {
    let scale = original.iter().flatten().chain(refreshed.iter().flatten())
        .fold(0.0_f64, |scale, value| scale.max(value.abs()));
    if original.iter().flatten().zip(refreshed.iter().flatten())
        .any(|(a, b)| (a-b).abs() > 1024.0*f64::EPSILON*scale) {
        return Err(PiezoelectricError::Invalid(
            "full thermal source must reproduce the original complete stress-charge tensor",
        ));
    }
    Ok(())
}

impl ThermalTensorPiezoelectric {
    /// Bind all original affine C^E, signed e, epsilon^S and constant-density
    /// source roles over one explicit positive absolute-temperature interval.
    ///
    /// Full native C^E and epsilon^S endpoint SPD establishes stability of the
    /// affine source throughout the interval, including directions not excited
    /// by a particular field. An affine strain-charge source is a different
    /// trajectory and is not silently reinterpreted as affine stress-charge.
    pub fn new(
        card: &MaterialCard,
        original: &TensorPiezoelectric,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self> {
        admit_source_temperature_refresh(original.material(), original.material())
            .map_err(field_error)?;
        let mut slopes = TensorTemperatureSlopes {
            stiffness: [[0.0; 6]; 6],
            coupling: [[0.0; 6]; 3],
            permittivity: [[0.0; 3]; 3],
        };
        for (prefix, rows, columns, block) in [
            ("piezoelectric_stiffness_e_mandel", 6, 6, 0),
            ("piezoelectric_stress_charge_mandel", 3, 6, 1),
            ("piezoelectric_permittivity_s", 3, 3, 2),
        ] {
            for i in 0..rows {
                for j in 0..columns {
                    let name = format!("{prefix}_{}{}", i+1, j+1);
                    let (_, slope, _) = affine_source_path(
                        card, original.material(), &name, temperature_axis,
                        temperature_range_k, block != 1 && i == j,
                    ).map_err(field_error)?;
                    match block {
                        0 => slopes.stiffness[i][j] = slope,
                        1 => slopes.coupling[i][j] = slope,
                        _ => slopes.permittivity[i][j] = slope,
                    }
                }
            }
        }
        let (_, density_slope, _) = affine_source_path(
            card, original.material(), DENSITY_PROPERTY, temperature_axis,
            temperature_range_k, true,
        ).map_err(field_error)?;
        if density_slope != 0.0 {
            return Err(PiezoelectricError::Invalid(
                "fixed-reference tensor piezoelectric thermal storage requires constant source density",
            ));
        }
        if *original.source_to_frame_rotation() != IDENTITY {
            slopes.rotate(*original.source_to_frame_rotation())?;
        }
        // A bundle may retain more than one tensor schema. Select the actual
        // requested source law, not inconsistent coefficients from another role.
        let refreshed = source_tensor(original, original.material().clone())?;
        admit_same_coefficients(original.stiffness_mandel_pa(), refreshed.stiffness_mandel_pa())?;
        admit_same_coefficients(original.coupling_c_m2(), refreshed.coupling_c_m2())?;
        admit_same_coefficients(original.permittivity_f_m(), refreshed.permittivity_f_m())?;
        let up = 1.0 + 4096.0*f64::EPSILON;
        let coupling_slope_norm = finite(frobenius_norm(&slopes.coupling)?*up)?;
        let permittivity_slope_norm = finite(frobenius_norm(&slopes.permittivity)?*up)?;
        let mut law = Self {
            card: card.clone(), original: original.clone(),
            temperature_axis: temperature_axis.to_owned(), temperature_range_k,
            slopes, minimum_permittivity_f_m: f64::INFINITY,
            maximum_coupling_norm_c_m2: 0.0,
            coupling_slope_norm, permittivity_slope_norm, bias: None,
        };
        for temperature in temperature_range_k {
            let at = law.at_temperature(temperature)?;
            let material = at.material();
            let mut inverse = [[0.0; 3]; 3];
            for j in 0..3 {
                let column = material.solve_dielectric(core::array::from_fn(|i| {
                    if i == j { 1.0 } else { 0.0 }
                }))?;
                for i in 0..3 { inverse[i][j] = column[i]; }
            }
            // lambda_min(epsilon) >= 1/||epsilon^-1||_F. Endpoint
            // bounds extend through an affine SPD interval by convexity.
            let lower = finite(1.0/frobenius_norm(&inverse)?)?;
            if lower <= 0.0 {
                return Err(PiezoelectricError::Invalid("tensor piezoelectric dielectric lower bound underflows"));
            }
            law.minimum_permittivity_f_m = law.minimum_permittivity_f_m.min(lower);
            law.maximum_coupling_norm_c_m2 = law.maximum_coupling_norm_c_m2
                .max(frobenius_norm(material.coupling_c_m2())?);
        }
        law.minimum_permittivity_f_m = finite(
            law.minimum_permittivity_f_m*(1.0-4096.0*f64::EPSILON),
        )?;
        law.maximum_coupling_norm_c_m2 = finite(law.maximum_coupling_norm_c_m2*up)?;
        if law.minimum_permittivity_f_m <= 0.0 {
            return Err(PiezoelectricError::Invalid("tensor piezoelectric dielectric interval bound underflows"));
        }
        Ok(law)
    }

    pub const fn card(&self) -> &MaterialCard { &self.card }
    pub const fn original(&self) -> &TensorPiezoelectric { &self.original }
    pub fn temperature_axis(&self) -> &str { &self.temperature_axis }
    pub const fn temperature_range_k(&self) -> [f64; 2] { self.temperature_range_k }

    /// Re-query every original source role at one temperature, retaining the
    /// exact source selections, fixed external coordinates and frame map.
    pub fn at_temperature(&self, temperature_k: f64) -> Result<ResolvedThermalTensorPiezoelectric> {
        if !temperature_k.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&temperature_k) {
            return Err(PiezoelectricError::Invalid("full piezoelectric temperature is outside its admitted source interval"));
        }
        let original_point = self.original.material().query();
        let point = original_point.clone().with(&self.temperature_axis, temperature_k)
            .map_err(|_| PiezoelectricError::Invalid("full piezoelectric thermal query is unrepresentable"))?;
        let source = self.original.material().at_query_point(&self.card, &point)?;
        // Refresh validates every source domain even at the reference point.
        // Pinning changes receipt-policy hashes, not the physical query or
        // source values; preserve the complete original snapshot at that point.
        let material = if point == original_point {
            self.original.clone()
        } else {
            source_tensor(&self.original, source)?
        };
        if material.density_kg_m3().to_bits() != self.original.density_kg_m3().to_bits() {
            return Err(PiezoelectricError::Invalid("full piezoelectric reference density changed"));
        }
        let bias = self.bias.as_ref().map(|bias| {
            bias.at_temperature(&self.card, &material, temperature_k)
        }).transpose()?;
        Ok(ResolvedThermalTensorPiezoelectric {
            material, temperature_k, slopes: self.slopes, bias,
        })
    }

    /// Admit the same actual tensor used by a spatial dynamics owner. Source
    /// schema, selected claims, coefficient values and named frame must agree.
    pub fn at_material(&self, material: &TensorPiezoelectric) -> Result<ResolvedThermalTensorPiezoelectric> {
        admit_source_temperature_refresh(self.original.material(), material.material())
            .map_err(field_error)?;
        let temperature = material.material().query_point().iter()
            .find_map(|(name, value)| (name == &self.temperature_axis).then_some(*value))
            .ok_or(PiezoelectricError::Invalid("full piezoelectric temperature coordinate is missing"))?;
        let mut at = self.at_temperature(temperature)?;
        // Selection policy is evidence metadata, not a constitutive coordinate.
        // Compare the actual complete law and frame after source/context
        // admission instead of requiring policy-sensitive identity hashes.
        if at.material.source_frame_identity() != material.source_frame_identity()
            || at.material.frame_identity() != material.frame_identity()
            || at.material.source_to_frame_rotation() != material.source_to_frame_rotation()
            || at.material.density_kg_m3().to_bits() != material.density_kg_m3().to_bits()
            || at.material.stiffness_mandel_pa() != material.stiffness_mandel_pa()
            || at.material.coupling_c_m2() != material.coupling_c_m2()
            || at.material.permittivity_f_m() != material.permittivity_f_m()
            || at.material.material().properties().iter().zip(material.material().properties())
                .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits()) {
            return Err(PiezoelectricError::Invalid(
                "full piezoelectric dynamics and thermal law selected different tensors or source values",
            ));
        }
        // Retain the caller's admitted actual tensor, factors and source receipt.
        at.material = material.clone();
        Ok(at)
    }

    /// Conservative whole-interval lower bound on this field's fixed-(S,D)
    /// heat capacity [J/K] for ||S_mandel||_2 <= maximum_strain_mandel_norm and
    /// ||D||_2 <= maximum_displacement_norm_c_m2. Mandel strain norm equals
    /// the physical symmetric tensor's Frobenius norm, including shear.
    ///
    /// With epsilon >= epsilon_min I and ||e|| <= e_max throughout,
    /// ||E|| <= (D_max+e_max S_max+||p|| max|T-Tref|)/epsilon_min.
    /// Thermal bias adds the directional bounds max(alpha.C.alpha),
    /// 2||C'alpha|| (S_max+||alpha|| max|T-Tref|), and the inverse-dielectric
    /// norm of p. These endpoint directional bounds retain anisotropy without
    /// replacing every thermal direction by the stiffest material direction.
    /// The caller enforces both coordinate limits and adds the actual
    /// stress-free, zero-field caloric heat capacity.
    /// This bound alone neither advances heat nor proves a relaxed field solve.
    pub fn minimum_heat_capacity_j_k(
        &self,
        maximum_strain_mandel_norm: f64,
        maximum_displacement_norm_c_m2: f64,
        volume_m3: f64,
    ) -> Result<f64> {
        if [maximum_strain_mandel_norm, maximum_displacement_norm_c_m2].iter()
            .any(|value| !value.is_finite() || *value < 0.0)
            || !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(PiezoelectricError::Invalid(
                "full piezoelectric capacity needs finite nonnegative coordinate bounds and positive volume",
            ));
        }
        let mechanical = finite(self.coupling_slope_norm*maximum_strain_mandel_norm)?;
        if self.coupling_slope_norm != 0.0 && maximum_strain_mandel_norm != 0.0 && mechanical == 0.0 {
            return Err(PiezoelectricError::Invalid("full piezoelectric coupling derivative bound underflows"));
        }
        let dielectric = if self.permittivity_slope_norm == 0.0 {
            0.0
        } else {
            let polarization = finite(self.maximum_coupling_norm_c_m2*maximum_strain_mandel_norm)?;
            if self.maximum_coupling_norm_c_m2 != 0.0 && maximum_strain_mandel_norm != 0.0 && polarization == 0.0 {
                return Err(PiezoelectricError::Invalid("full piezoelectric polarization bound underflows"));
            }
            let pyroelectric = self.bias.as_ref().map_or(0.0, |bias| {
                bias.polarization_norm * bias.maximum_temperature_change_k
            });
            let charge = finite(maximum_displacement_norm_c_m2+polarization+pyroelectric)?;
            let field = finite(charge/self.minimum_permittivity_f_m)?;
            let derivative = finite(self.permittivity_slope_norm*field)?;
            if charge != 0.0 && (field == 0.0 || derivative == 0.0) {
                return Err(PiezoelectricError::Invalid("full piezoelectric dielectric derivative bound underflows"));
            }
            derivative
        };
        let derivative = finite(mechanical+dielectric)?;
        let action = finite(derivative/fs_math::det::sqrt(self.minimum_permittivity_f_m))?;
        if derivative != 0.0 && action == 0.0 {
            return Err(PiezoelectricError::Invalid("full piezoelectric capacity bound underflows"));
        }
        let (action, elastic_curvature) = if let Some(bias) = &self.bias {
            let strain = finite(maximum_strain_mandel_norm
                + bias.expansion_norm*bias.maximum_temperature_change_k)?;
            (finite(action+bias.maximum_pyroelectric_inverse_norm)?,
                finite(bias.maximum_expansion_curvature
                    + 2.0*bias.expansion_slope_norm*strain)?)
        } else { (action, 0.0) };
        let curvature = finite(action*action+elastic_curvature)?;
        if curvature == 0.0 {
            return if action == 0.0 { Ok(0.0) }
            else { Err(PiezoelectricError::Invalid("full piezoelectric capacity bound underflows")) };
        }
        let magnitude = finite(self.temperature_range_k[1]*volume_m3*curvature)?;
        if magnitude == 0.0 {
            return Err(PiezoelectricError::Invalid("full piezoelectric capacity bound underflows"));
        }
        finite(-magnitude*(1.0+4096.0*f64::EPSILON))
    }
}

impl ResolvedThermalTensorPiezoelectric {
    pub const fn material(&self) -> &TensorPiezoelectric { &self.material }
    pub const fn temperature_k(&self) -> f64 { self.temperature_k }

    /// Integrated Helmholtz energy, entropy, internal energy and heat capacity
    /// at fixed complete strain/displacement and fixed positive reference volume.
    /// The native material factors supply psi and E with the sourced offsets.
    /// With changing C(T) and thermal expansion the field capacity can have
    /// either sign. Finite caloric storage uses the stress-free zero-field chart;
    /// reversible field entropy is never counted as dissipative ceramic heat.
    pub fn state(
        &self,
        strain_mandel: [f64; 6],
        displacement_c_m2: [f64; 3],
        volume_m3: f64,
    ) -> Result<FieldThermalState> {
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(PiezoelectricError::Invalid("full piezoelectric thermal volume must be finite positive"));
        }
        let response = self.charge_response(strain_mandel, displacement_c_m2)?;
        let field = response.electric_field_v_m;
        let thermal_strain = self.thermal_strain_mandel();
        let elastic_strain = core::array::from_fn(|i| strain_mandel[i]-thermal_strain[i]);
        let elastic_slope = apply(&self.slopes.stiffness, &elastic_strain);
        let mut polarization_slope = apply(&self.slopes.coupling, &strain_mandel);
        let (expansion_entropy, expansion_curvature) = if let Some(bias) = &self.bias {
            for (value, p) in polarization_slope.iter_mut().zip(bias.pyroelectric) { *value += p; }
            (dot(&bias.expansion, &apply(self.material.stiffness_mandel_pa(), &elastic_strain)),
                bias.expansion_curvature-2.0*dot(&bias.expansion, &elastic_slope))
        } else { (0.0, 0.0) };
        let dielectric_slope = apply(&self.slopes.permittivity, &field);
        let entropy_density = finite(
            (-0.5_f64).mul_add(dot(&elastic_strain, &elastic_slope),
                0.5_f64.mul_add(dot(&field, &dielectric_slope),
                    dot(&field, &polarization_slope)+expansion_entropy)),
        )?;
        let drive = core::array::from_fn(|i| polarization_slope[i]+dielectric_slope[i]);
        let action = self.material.dielectric_inverse_factor_action(drive)?;
        // The dielectric curvature remains a nonnegative squared native factor
        // action; the thermoelastic cross term can change the complete sign.
        let dielectric_curvature = finite(dot(&action, &action))?;
        let curvature = finite(dielectric_curvature+expansion_curvature)?;
        let free_energy_j = finite(volume_m3*response.stored_energy_j_m3)?;
        let entropy_j_k = finite(volume_m3*entropy_density)?;
        let heat_capacity_j_k = finite(-self.temperature_k*volume_m3*curvature)?;
        let internal_energy_j = finite(self.temperature_k.mul_add(entropy_j_k, free_energy_j))?;
        if response.stored_energy_j_m3 != 0.0 && free_energy_j == 0.0
            || entropy_density != 0.0 && entropy_j_k == 0.0
            || action.iter().any(|&value| value != 0.0) && dielectric_curvature == 0.0
            || curvature != 0.0 && heat_capacity_j_k == 0.0 {
            return Err(PiezoelectricError::Invalid("full piezoelectric thermal contribution underflows"));
        }
        Ok(FieldThermalState { free_energy_j, entropy_j_k, internal_energy_j, heat_capacity_j_k })
    }
}
