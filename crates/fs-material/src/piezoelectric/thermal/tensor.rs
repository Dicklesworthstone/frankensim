//! Thermal reduction of the complete original affine strain-charge tensor.
//!
//! At fixed direction n and r=mandel(n⊗n), the reciprocal scalar compliance
//! block is the projection P^T A(T) P of A=[[s^E,d^T],[d,epsilon^T]]. A is
//! affine when every original source coefficient is affine. Positive complete
//! endpoint blocks then guarantee positivity throughout the interval. This
//! retains transverse relaxation and every off-diagonal thermal derivative.
//!
//! The source convention is the strain-charge form in COMSOL's Piezoelectricity
//! theory; the parent module supplies the same Helmholtz derivatives:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.043.html>.
//! <https://web.mit.edu/beretta/www/AdvancedThermodynamics_2_43_Lecture6.htm>.

use super::*;
use crate::piezoelectric::tensor::TensorPiezoelectric;
use crate::tensor::mandel_basis;

mod stress_charge;

impl TensorThermalReduction {
    fn tensor(
        &self,
        source: ResolvedMaterialStatePoint,
    ) -> Result<TensorPiezoelectric, PiezoelectricError> {
        let tensor = if self.stress_charge.is_some() {
            TensorPiezoelectric::from_material_state(source, self.source_frame)?
        } else {
            TensorPiezoelectric::from_strain_charge_material_state(source, self.source_frame)?
        };
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let tensor = if self.frame == self.source_frame && self.source_to_frame == identity {
            tensor
        } else {
            tensor.rotated(self.frame, self.source_to_frame)?
        };
        Ok(tensor)
    }

    pub(super) fn material(
        &self,
        source: ResolvedMaterialStatePoint,
    ) -> Result<(LongitudinalPiezoelectric, Option<TensorThermalDerivatives>), PiezoelectricError> {
        let tensor = self.tensor(source)?;
        let material = tensor.longitudinal(self.direction)?;
        let direction = material.tensor_reduction().expect("native reduction receipt").direction;
        let derivatives = self.stress_charge.as_ref()
            .map(|source| source.derivatives(&tensor, direction))
            .transpose()?;
        Ok((material, derivatives))
    }
}

impl ThermalLongitudinalPiezoelectric {
    /// Bind a full anisotropic strain-charge source to the existing longitudinal
    /// rod/circuit thermal law. The original source tensor must declare all
    /// s^E, signed d, epsilon^T and density roles. For affine C^E, e and
    /// epsilon^S input, use `from_stress_charge_tensor`, which retains the
    /// non-affine inverse derivatives.
    ///
    /// The fixed physical projection imposes T=t n⊗n and E=E_n n. Transverse
    /// strain and electric displacement relax exactly as in `longitudinal`.
    /// Original frame transformations, all source selections and all external
    /// coordinates survive every temperature trial. Both COMPLETE 9x9 coupled
    /// endpoint blocks must be positive, including directions the rod does not
    /// excite. Binding and temperature refresh use the original tensor owner;
    /// repeated storage/entropy/heat-capacity responses remain scalar.
    ///
    /// This is a small-signal homogeneous constraint, not an insulating lateral
    /// boundary solve, thermal-expansion or pyroelectric offset model. The
    /// synthetic G0/G1/G3 tests are implementation references, not measured
    /// ceramic validation, a runtime qualification or a numerical certificate.
    pub fn from_tensor(
        card: &MaterialCard,
        original: &TensorPiezoelectric,
        direction: [f64; 3],
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, PiezoelectricError> {
        let longitudinal = original.longitudinal(direction)?;
        let direction = longitudinal
            .tensor_reduction()
            .expect("native reduction receipt")
            .direction;
        let mut reduction = TensorThermalReduction {
            source_frame: original.source_frame_identity(),
            frame: original.frame_identity(),
            source_to_frame: *original.source_to_frame_rotation(),
            direction,
            coefficient_roundoff_scales: [0.0; 3],
            stress_charge: None,
        };
        // Original source coordinates, rather than scalar surrogate claims,
        // admit all 64 roles and retain their exact source-domain boundaries.
        // The basis change is fixed throughout this thermal trajectory.
        let n: [f64; 3] = core::array::from_fn(|i| {
            (0..3).fold(0.0, |sum, j| {
                reduction.source_to_frame[j][i].mul_add(direction[j], sum)
            })
        });
        let basis = mandel_basis();
        let r: [f64; 6] = core::array::from_fn(|a| {
            (0..3).fold(0.0, |sum, i| {
                (0..3).fold(sum, |sum, j| {
                    (basis[a][i][j] * n[i]).mul_add(n[j], sum)
                })
            })
        });
        let mut slopes: [f64; 3] = [0.0; 3];
        let mut reference_temperature_k = 0.0;
        for (prefix, rows, columns, output) in [
            ("piezoelectric_compliance_e_mandel", 6, 6, 0),
            ("piezoelectric_strain_charge_mandel", 3, 6, 1),
            ("piezoelectric_permittivity_t", 3, 3, 2),
        ] {
            for i in 0..rows {
                for j in 0..columns {
                    let name = format!("{prefix}_{}{}", i + 1, j + 1);
                    let (reference, slope, endpoints) = affine_source_path(
                        card,
                        original.material(),
                        &name,
                        temperature_axis,
                        temperature_range_k,
                        output != 1 && i == j,
                    )
                    .map_err(field_error)?;
                    reference_temperature_k = reference;
                    let weight = match output {
                        0 => r[i] * r[j],
                        1 => n[i] * r[j],
                        _ => n[i] * n[j],
                    };
                    slopes[output] = finite(slope.mul_add(weight, slopes[output]))?;
                    let scale = endpoints[0].abs().max(endpoints[1].abs());
                    reduction.coefficient_roundoff_scales[output] = finite(
                        weight.abs().mul_add(
                            scale,
                            reduction.coefficient_roundoff_scales[output],
                        ),
                    )?;
                }
            }
        }
        Self::prepare(
            card,
            &longitudinal,
            temperature_axis,
            temperature_range_k,
            reference_temperature_k,
            slopes,
            Some(reduction),
        )
    }
}
