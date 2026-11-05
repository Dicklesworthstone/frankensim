//! Passive material memory around the original reciprocal piezoelectric tensor.
//!
//! The base stiffness is the RELAXED constant-field stiffness C0; its clamped
//! permittivity is the INSTANTANEOUS epsilon_inf. These roles must be declared
//! by the source model: adding poles to a measured static permittivity would
//! count their relaxed polarization twice. The original real coupling e stays
//! reciprocal. No empirical quality factor or complex coupling is inferred.
//!
//! Native Maxwell branches C_j=B_j^T B_j have qdot=B_j Sdot-q/tau. A directed
//! Debye pole has pdot=(delta_epsilon n.E-p)/tau. Together their storage is
//!   psi = S.C0.S/2 + (D-eS-sum(n p)).epsilon_inf^-1.(D-eS-sum(n p))/2
//!       + sum(|q|^2)/2 + sum(p^2/delta_epsilon)/2.
//! Thus sigma=C0 S-e^T E+sum(B^T q), D=eS+epsilon_inf E+sum(n p), and heat is
//! sum(|q|^2/tau)+sum(tau pdot^2/delta_epsilon), always nonnegative. Mechanical
//! branches use the original tensor-Prony PSD factors and Mandel rotation.
//! Directed positive dielectric poles admit rank-deficient anisotropy without
//! diagonalizing the piezoelectric tensor or assuming common principal axes.
//!
//! Peak phasors use exp(-i omega t). Frequency coefficients and field recovery
//! allocate nothing; source queries and frame preparation occur only at binding.
//! This supplies a fixed-temperature harmonic constitutive response, not a
//! transient integrator, thermal update, hysteresis, depoling or measured fit.

use super::{Result, TensorPiezoelectric, invalid};
use crate::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::tensor::mandel_rotation;
use crate::visco::tensor::{PronyBranch, TensorProny};
use crate::visco_tensor::source::TensorBranchProperties;
use fs_blake3::ContentHash;
use fs_matdb::MaterialCard;
use fs_math::{c64::C64, det};
use fs_qty::{Dims, Pressure, QuantitySpec, Time};
use std::sync::Arc;

/// Original full Mandel stiffness increments; zero/rank-deficient branches are valid.
pub const MAX_MECHANICAL_BRANCHES: usize = 64;
/// Original directed polarization processes, independent of the mesh size.
pub const MAX_DIELECTRIC_POLES: usize = 16;
/// Existing material-query bound on distinct additional property names. Reused
/// same-unit names resolve once; the complete base tensor has its own bundle.
pub const MAX_RELAXATION_SOURCE_PROPERTIES: usize =
    crate::state_point::MAX_MATERIAL_STATE_PROPERTIES;

/// One positive dielectric increment and its direction in the BASE SOURCE frame.
/// All three dimensionless direction components are explicitly sourced, including
/// zeros. Their norm must already be one; binding never normalizes source data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DirectedDebyeProperties {
    /// Positive absolute permittivity increment [F/m].
    pub permittivity_increment: String,
    /// Positive polarization relaxation time [s].
    pub relaxation_time: String,
    /// Source-property names for the oriented unit vector [nx,ny,nz].
    pub direction: [String; 3],
}

/// Source-derived polarization process, already expressed in the current frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DirectedDebyePole {
    pub permittivity_increment_f_m: f64,
    pub relaxation_time_s: f64,
    pub direction: [f64; 3],
}

/// Original base tensor plus source-resolved passive material memory.
/// The complete base and auxiliary source bundles remain separately inspectable.
#[derive(Clone, Debug)]
pub struct RelaxingTensorPiezoelectric {
    base: TensorPiezoelectric,
    material: Arc<ResolvedMaterialStatePoint>,
    mechanical_properties: Arc<[TensorBranchProperties]>,
    dielectric_properties: Arc<[DirectedDebyeProperties]>,
    mechanical: TensorProny,
    dielectric: Arc<[DirectedDebyePole]>,
}

/// Peak complex constitutive coefficients for exp(-i omega t).
/// Increments are accumulated DIRECTLY, preserving weak loss/storage next to a
/// large base coefficient. Consumers must not recover them by subtraction.
#[derive(Clone, Copy, Debug)]
pub struct PiezoelectricFrequencyCoefficients {
    pub stiffness_mandel_pa: [[C64; 6]; 6],
    pub stiffness_increment_mandel_pa: [[C64; 6]; 6],
    pub permittivity_f_m: [[C64; 3]; 3],
    pub permittivity_increment_f_m: [[C64; 3]; 3],
}

/// Allocation-free harmonic fields and their independent physical storage/loss.
/// Memory energies are subaccounts of total stored energy, not additional heat.
#[derive(Clone, Copy, Debug)]
pub struct RelaxingPiezoelectricHarmonicResponse {
    pub stress_mandel_pa: [C64; 6],
    pub displacement_c_m2: [C64; 3],
    /// Stress supplied by the retained Maxwell memories only [Pa].
    pub additional_stress_mandel_pa: [C64; 6],
    /// Displacement supplied by the retained Debye polarizations only [C/m²].
    pub polarization_c_m2: [C64; 3],
    pub average_stored_energy_j_m3: f64,
    pub average_mechanical_memory_energy_j_m3: f64,
    pub average_dielectric_memory_energy_j_m3: f64,
    pub average_mechanical_heat_w_m3: f64,
    pub average_dielectric_heat_w_m3: f64,
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn scalar(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "piezoelectric relaxation response is not representable",
        ))
    }
}

/// Maxwell's (-i x)/(1-i x) and Debye's 1/(1-i x), without subtracting
/// nearly equal complex numbers or overflowing x squared. The reciprocal path
/// also retains a subnormal 1/x when the intermediate omega*tau overflows.
fn relaxation_factors(omega: f64, tau: f64) -> Result<(C64, C64)> {
    if omega == 0.0 {
        return Ok((C64::ZERO, C64::ONE));
    }
    let x = omega * tau;
    if x == 0.0 {
        return Err(invalid("piezoelectric relaxation rate ratio underflows"));
    }
    let (maxwell, debye) = if x <= 1.0 {
        let denominator = 1.0 + x * x;
        (
            C64::new(x * x / denominator, -x / denominator),
            C64::new(1.0 / denominator, x / denominator),
        )
    } else {
        let inverse = (1.0 / omega) / tau;
        if inverse == 0.0 {
            return Err(invalid(
                "piezoelectric relaxation inverse rate ratio underflows",
            ));
        }
        let denominator = 1.0 + inverse * inverse;
        (
            C64::new(1.0 / denominator, -inverse / denominator),
            C64::new(inverse * inverse / denominator, inverse / denominator),
        )
    };
    Ok((maxwell, debye))
}

fn requirements(
    mechanical: &[TensorBranchProperties],
    dielectric: &[DirectedDebyeProperties],
) -> Result<Vec<ScalarPropertyRequirement>> {
    if mechanical.len() > MAX_MECHANICAL_BRANCHES || dielectric.len() > MAX_DIELECTRIC_POLES {
        return Err(invalid(
            "piezoelectric relaxation exceeds its bounded source count",
        ));
    }
    let mut result: Vec<ScalarPropertyRequirement> = Vec::new();
    let mut add = |name: &str, dims, admission| -> Result<()> {
        if name.is_empty() || name.len() > 128 {
            return Err(invalid(
                "relaxation properties require names of 1..=128 bytes",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_new(name, dims, admission)?;
        if let Some(previous) = result.iter().find(|r| r.name() == name) {
            if previous.quantity() != requirement.quantity() {
                return Err(invalid(
                    "one relaxation source name has incompatible quantity roles",
                ));
            }
        } else {
            result.push(requirement);
            if result.len() > MAX_RELAXATION_SOURCE_PROPERTIES {
                return Err(invalid(
                    "relaxation exceeds the native 64-property source bundle limit",
                ));
            }
        }
        Ok(())
    };
    for branch in mechanical {
        for name in &branch.stiffness_upper {
            add(name, Pressure::DIMS, ScalarAdmissibility::Finite)?;
        }
        add(
            &branch.relaxation_time,
            Time::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
    }
    for pole in dielectric {
        add(
            &pole.permittivity_increment,
            ELECTRIC_PERMITTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        add(
            &pole.relaxation_time,
            Time::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        for name in &pole.direction {
            add(name, Dims::NONE, ScalarAdmissibility::Finite)?;
        }
    }
    Ok(result)
}

fn rotate_direction(rotation: [[f64; 3]; 3], direction: [f64; 3]) -> Result<[f64; 3]> {
    let result = rotation.map(|row| {
        row.into_iter()
            .zip(direction)
            .fold(0.0, |sum, (a, b)| a.mul_add(b, sum))
    });
    if result.iter().any(|x| !x.is_finite()) {
        return Err(invalid(
            "relaxation direction rotation is not representable",
        ));
    }
    Ok(result)
}

impl RelaxingTensorPiezoelectric {
    /// Bind material memory to the SAME immutable card and complete physical
    /// query as the original tensor. The base represents C0 and epsilon_inf;
    /// original source claims must support those respective limiting roles.
    /// Empty lists give the original reversible material exactly.
    ///
    /// Every branch tensor is in the original base source frame. Native PSD
    /// factors preserve signed off-diagonal entries and rank-deficient loss.
    /// Source parameters are constant at the held external environment; an
    /// applicability-only frequency band is intersected across ALL selections.
    /// The original query owner admits at most 64 distinct additional names.
    pub fn resolve(
        base: TensorPiezoelectric,
        card: &MaterialCard,
        mechanical: Vec<TensorBranchProperties>,
        dielectric: Vec<DirectedDebyeProperties>,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        if card.content_hash() != base.material().card_identity() {
            return Err(invalid(
                "relaxation requires the original piezoelectric material card",
            ));
        }
        let requirements = requirements(&mechanical, &dielectric)?;
        let material = if requirements.is_empty() {
            base.material().clone()
        } else {
            resolve_material_state_point(card, &base.material().query(), &requirements, selection)?
        };
        if material.query_point() != base.material().query_point()
            || material.axis_quantities() != base.material().axis_quantities()
        {
            return Err(invalid(
                "piezoelectric and relaxation sources require the same physical query",
            ));
        }
        let mut band = (0.0_f64, f64::MAX);
        for property in material.properties() {
            if !property.value_si().is_finite() || !property.is_constant_at_fixed_environment() {
                return Err(invalid(
                    "relaxation parameters must be constant at the fixed environment",
                ));
            }
            if let Some((low, high)) = property.angular_frequency_band_rad_s()? {
                band = (band.0.max(low), band.1.min(high));
            }
        }
        for requirement in &requirements {
            let property = material
                .property(requirement.name())
                .ok_or_else(|| invalid("original relaxation property is missing"))?;
            if property.requirement().quantity() != QuantitySpec::dimensional(requirement.dims()) {
                return Err(invalid(
                    "relaxation coefficients require their exact declared SI units",
                ));
            }
            if let Some(original) = base.material().property(requirement.name()) {
                if original.requirement().quantity() != property.requirement().quantity()
                    || original.answer().receipt.selected != property.answer().receipt.selected
                    || original.value_si().to_bits() != property.value_si().to_bits()
                {
                    return Err(invalid(
                        "relaxation cannot replace an already selected base property",
                    ));
                }
            }
        }
        let value = |name: &str| -> f64 {
            material
                .property(name)
                .expect("complete admitted relaxation requirements")
                .value_si()
        };
        let source_to_frame = *base.source_to_frame_rotation();
        let rotation = mandel_rotation(source_to_frame);
        let mut branches = Vec::with_capacity(mechanical.len());
        for properties in &mechanical {
            let mut stiffness = [[0.0; 6]; 6];
            let mut entry = 0;
            for row in 0..6 {
                for column in row..6 {
                    stiffness[row][column] = value(&properties.stiffness_upper[entry]);
                    stiffness[column][row] = stiffness[row][column];
                    entry += 1;
                }
            }
            let source =
                PronyBranch::new(stiffness, value(&properties.relaxation_time)).map_err(|_| {
                    invalid("mechanical relaxation requires a PSD Mandel tensor and positive time")
                })?;
            let gram = source
                .gram_rows
                .as_ref()
                .ok_or_else(|| invalid("native Prony branch did not retain its admitted factor"))?;
            // SAME native Mandel map as the complete base. Transform B, not a
            // refactored/diagonalized loss tensor, to retain null directions.
            let rows = core::array::from_fn(|i| {
                core::array::from_fn(|j| {
                    (0..6).fold(0.0, |sum, k| gram[i][k].mul_add(rotation[j][k], sum))
                })
            });
            branches.push(
                PronyBranch::from_gram_rows(rows, source.tau_s).map_err(|_| {
                    invalid("rotated mechanical relaxation factor is not representable")
                })?,
            );
        }
        // This native constructor preserves the base C0 coefficients exactly.
        let mechanical_law = TensorProny::new(*base.stiffness_mandel_pa(), None, branches, band)
            .map_err(|_| {
                invalid("complete mechanical relaxation or frequency band is inadmissible")
            })?;
        let mut poles = Vec::with_capacity(dielectric.len());
        for properties in &dielectric {
            let direction = properties.direction.each_ref().map(|name| value(name));
            let norm = direction.into_iter().fold(0.0, |sum, x| x.mul_add(x, sum));
            let delta = value(&properties.permittivity_increment);
            let tau = value(&properties.relaxation_time);
            if !norm.is_finite()
                || (norm - 1.0).abs() > 1.0e-12
                || !delta.is_finite()
                || delta <= 0.0
                || !tau.is_finite()
                || tau <= 0.0
            {
                return Err(invalid(
                    "Debye poles require positive increment/time and an explicitly unit direction",
                ));
            }
            poles.push(DirectedDebyePole {
                permittivity_increment_f_m: delta,
                relaxation_time_s: tau,
                direction: rotate_direction(source_to_frame, direction)?,
            });
        }
        let model = Self {
            base,
            material: Arc::new(material),
            mechanical_properties: mechanical.into(),
            dielectric_properties: dielectric.into(),
            mechanical: mechanical_law,
            dielectric: poles.into(),
        };
        // Positive increments must still have a representable relaxed total.
        let mut relaxed = *model.base.permittivity_f_m();
        for pole in model.dielectric.iter() {
            for row in 0..3 {
                for column in 0..3 {
                    relaxed[row][column] = scalar(
                        relaxed[row][column]
                            + pole.permittivity_increment_f_m
                                * pole.direction[row]
                                * pole.direction[column],
                    )?;
                }
            }
        }
        Ok(model)
    }

    /// Original reciprocal tensor in the current frame, including its identity.
    pub const fn base(&self) -> &TensorPiezoelectric {
        &self.base
    }
    /// Original additional source selections at the exact base physical query.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Native Maxwell law, retaining original C0, PSD factors and source band.
    pub const fn mechanical_law(&self) -> &TensorProny {
        &self.mechanical
    }
    pub fn mechanical_properties(&self) -> &[TensorBranchProperties] {
        &self.mechanical_properties
    }
    pub fn dielectric_properties(&self) -> &[DirectedDebyeProperties] {
        &self.dielectric_properties
    }
    pub fn dielectric_poles(&self) -> &[DirectedDebyePole] {
        &self.dielectric
    }
    pub fn frequency_band_rad_s(&self) -> [f64; 2] {
        let (low, high) = self.mechanical.band_rad_s();
        [low, high]
    }

    /// Co-rotate the complete reciprocal base and every memory direction. All
    /// original source selections and relaxation times remain unchanged.
    pub fn rotated(&self, target_frame: ContentHash, rotation: [[f64; 3]; 3]) -> Result<Self> {
        let base = self.base.rotated(target_frame, rotation)?;
        let mechanical = self
            .mechanical
            .rotated(rotation)
            .map_err(|_| invalid("native mechanical relaxation frame rotation refused"))?;
        let dielectric = self
            .dielectric
            .iter()
            .map(|pole| {
                Ok(DirectedDebyePole {
                    direction: rotate_direction(rotation, pole.direction)?,
                    ..*pole
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            base,
            mechanical,
            dielectric: dielectric.into(),
            material: self.material.clone(),
            mechanical_properties: self.mechanical_properties.clone(),
            dielectric_properties: self.dielectric_properties.clone(),
        })
    }

    fn admit_frequency(&self, omega: f64) -> Result<()> {
        let [low, high] = self.frequency_band_rad_s();
        if !omega.is_finite() || omega < low || omega > high {
            return Err(invalid(
                "frequency leaves the original piezoelectric relaxation source band",
            ));
        }
        Ok(())
    }

    /// Frequency-domain coefficients from the declared causal memories. No
    /// large-base subtraction, source query, allocation or matrix inversion.
    pub fn frequency_coefficients(&self, omega: f64) -> Result<PiezoelectricFrequencyCoefficients> {
        self.admit_frequency(omega)?;
        let mut stiffness_increment = [[C64::ZERO; 6]; 6];
        let mut permittivity_increment = [[C64::ZERO; 3]; 3];
        for branch in self.mechanical.branches() {
            let (factor, _) = relaxation_factors(omega, branch.tau_s)?;
            for row in 0..6 {
                for column in 0..6 {
                    stiffness_increment[row][column] =
                        stiffness_increment[row][column] + factor.scale(branch.c_pa[row][column]);
                }
            }
        }
        for pole in self.dielectric.iter() {
            let (_, factor) = relaxation_factors(omega, pole.relaxation_time_s)?;
            for row in 0..3 {
                for column in 0..3 {
                    permittivity_increment[row][column] = permittivity_increment[row][column]
                        + factor.scale(
                            pole.permittivity_increment_f_m
                                * pole.direction[row]
                                * pole.direction[column],
                        );
                }
            }
        }
        let stiffness = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                stiffness_increment[i][j] + C64::from_re(self.base.stiffness_mandel_pa()[i][j])
            })
        });
        let permittivity = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                permittivity_increment[i][j] + C64::from_re(self.base.permittivity_f_m()[i][j])
            })
        });
        if stiffness
            .iter()
            .flatten()
            .chain(stiffness_increment.iter().flatten())
            .chain(permittivity.iter().flatten())
            .chain(permittivity_increment.iter().flatten())
            .any(|&z| !finite(z))
        {
            return Err(invalid(
                "complex piezoelectric relaxation coefficients are not representable",
            ));
        }
        Ok(PiezoelectricFrequencyCoefficients {
            stiffness_mandel_pa: stiffness,
            stiffness_increment_mandel_pa: stiffness_increment,
            permittivity_f_m: permittivity,
            permittivity_increment_f_m: permittivity_increment,
        })
    }

    /// Recover actual phasor fields and independent memory storage/heat. Each
    /// mechanical contribution is a native Gram norm; each dielectric one is
    /// a positive polarization/friction norm. These are not inferred from a
    /// difference of total input power or imaginary diagonal coefficients.
    /// Field amplitudes remain subject to the consuming spatial owner's bounds.
    /// At zero frequency these means denote the continuous sinusoidal limit,
    /// not the stored energy of a static field with the same real amplitude.
    pub fn harmonic_response(
        &self,
        omega: f64,
        strain_mandel: [C64; 6],
        electric_field_v_m: [C64; 3],
    ) -> Result<RelaxingPiezoelectricHarmonicResponse> {
        self.admit_frequency(omega)?;
        if strain_mandel
            .iter()
            .chain(&electric_field_v_m)
            .any(|&z| !finite(z))
        {
            return Err(invalid(
                "piezoelectric relaxation requires finite peak strain and field",
            ));
        }
        let re = self.base.field_response(
            strain_mandel.map(|z| z.re),
            electric_field_v_m.map(|z| z.re),
        )?;
        let im = self.base.field_response(
            strain_mandel.map(|z| z.im),
            electric_field_v_m.map(|z| z.im),
        )?;
        let mut additional_stress = [C64::ZERO; 6];
        let mut polarization = [C64::ZERO; 3];
        let mut mechanical_energy = 0.0;
        let mut dielectric_energy = 0.0;
        let mut mechanical_heat = 0.0;
        let mut dielectric_heat = 0.0;
        for (branch, factor) in self
            .mechanical
            .branches()
            .iter()
            .zip(self.mechanical.factors())
        {
            let (weight, _) = relaxation_factors(omega, branch.tau_s)?;
            let gram = branch.resolved_gram_rows(factor);
            let heat_scale = core::f64::consts::FRAC_1_SQRT_2 / det::sqrt(branch.tau_s);
            for row in gram {
                let amplitude = row
                    .iter()
                    .zip(&strain_mandel)
                    .fold(C64::ZERO, |sum, (&a, &s)| sum + s.scale(a));
                let memory = weight * amplitude;
                mechanical_energy = scalar(mechanical_energy + memory.scale(0.5).norm_sq())?;
                mechanical_heat = scalar(mechanical_heat + memory.scale(heat_scale).norm_sq())?;
                for i in 0..6 {
                    additional_stress[i] = additional_stress[i] + memory.scale(row[i]);
                }
            }
        }
        for pole in self.dielectric.iter() {
            let (friction_weight, polarization_weight) =
                relaxation_factors(omega, pole.relaxation_time_s)?;
            let root = det::sqrt(pole.permittivity_increment_f_m);
            let along = pole
                .direction
                .iter()
                .zip(&electric_field_v_m)
                .fold(C64::ZERO, |sum, (&n, &e)| sum + e.scale(n));
            let amplitude = along.scale(root);
            let memory = polarization_weight * amplitude;
            let friction = friction_weight * amplitude;
            let heat_scale = core::f64::consts::FRAC_1_SQRT_2 / det::sqrt(pole.relaxation_time_s);
            dielectric_energy = scalar(dielectric_energy + memory.scale(0.5).norm_sq())?;
            dielectric_heat = scalar(dielectric_heat + friction.scale(heat_scale).norm_sq())?;
            for i in 0..3 {
                polarization[i] = polarization[i] + memory.scale(root * pole.direction[i]);
            }
        }
        let stress = core::array::from_fn(|i| {
            C64::new(re.stress_mandel_pa[i], im.stress_mandel_pa[i]) + additional_stress[i]
        });
        let displacement = core::array::from_fn(|i| {
            C64::new(re.displacement_c_m2[i], im.displacement_c_m2[i]) + polarization[i]
        });
        if stress
            .iter()
            .chain(&additional_stress)
            .chain(&displacement)
            .chain(&polarization)
            .any(|&z| !finite(z))
        {
            return Err(invalid(
                "piezoelectric relaxation fields are not representable",
            ));
        }
        let energy = scalar(
            0.5 * re.stored_energy_j_m3
                + 0.5 * im.stored_energy_j_m3
                + mechanical_energy
                + dielectric_energy,
        )?;
        Ok(RelaxingPiezoelectricHarmonicResponse {
            stress_mandel_pa: stress,
            displacement_c_m2: displacement,
            additional_stress_mandel_pa: additional_stress,
            polarization_c_m2: polarization,
            average_stored_energy_j_m3: energy,
            average_mechanical_memory_energy_j_m3: mechanical_energy,
            average_dielectric_memory_energy_j_m3: dielectric_energy,
            average_mechanical_heat_w_m3: mechanical_heat,
            average_dielectric_heat_w_m3: dielectric_heat,
        })
    }
}
