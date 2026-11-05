//! Local stresses from the same two-probe equilibrium field that defines H.
//!
//! A first-order shell's section shear strain is an energy coordinate. Applying
//! each ply's G directly to that constant coordinate does not recover its
//! physical transverse stress. Here the original quadratic traction profile
//! supplies continuous interface traction and zero exterior-face traction;
//! the ply's uncorrected compliance then supplies its local engineering strain.
//! This remains the declared cylindrical-bending reduction, not a free-edge,
//! delamination, transverse-normal-stress or general 3-D warping solution.

use super::{equilibrated_with_profiles, invalid, original_layers};
use crate::laminate::{LaminateError, LaminateSection, LaminateShearModel};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
struct LayerProfile {
    // Dimensionless polynomial in fraction through this layer. Dividing the
    // resultant by total thickness gives physical stress in Pa.
    traction: [[f64; 4]; 3],
    // Cholesky factor of the ORIGINAL rotated, uncorrected physical G [Pa].
    lower: [f64; 3],
    sin: f64,
    cos: f64,
}

/// Immutable O(layers) preparation retaining the exact source section.
/// Local evaluation is allocation-free and constant work, independent of the
/// number of layers. A single preparation may be shared by many shell facets.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateShearRecovery {
    section: Arc<LaminateSection>,
    profiles: Vec<LayerProfile>,
}

/// Local transverse shear at one physical point within an original ply.
/// Components are physical xz/yz stresses and engineering xz/yz strains;
/// material components instead use the original ply's 13/23 directions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerShearResponse {
    /// Coordinate relative to the original reference surface [m].
    pub z_m: f64,
    /// Undoubled physical [tau_xz, tau_yz] in section axes [Pa].
    pub stress_pa: [f64; 2],
    /// Physical engineering [gamma_xz, gamma_yz] in section axes.
    pub engineering_strain: [f64; 2],
    /// Undoubled physical [tau_13, tau_23] in original ply axes [Pa].
    pub material_stress_pa: [f64; 2],
    /// Physical engineering [gamma_13, gamma_23] in original ply axes.
    pub material_engineering_strain: [f64; 2],
    /// Local recoverable shear energy density, tau dot gamma / 2 [J/m^3].
    pub energy_density_j_m3: f64,
}

impl LaminateShearRecovery {
    /// Retain the original equilibrated elastic section and prepare its
    /// continuous traction field. Directional Whitney and Prony reductions
    /// refuse rather than receiving a different implied recovery model.
    /// No constitutive fit, correction factor or material query is introduced.
    pub fn try_new(section: Arc<LaminateSection>) -> Result<Self, LaminateError> {
        let LaminateShearModel::EquilibratedCylindrical { bending_axis_rad } =
            section.shear_model()
        else {
            return Err(invalid(
                "shear recovery requires the original equilibrated cylindrical model",
            ));
        };
        let layers = original_layers(&section)?;
        let mut polynomials = Vec::with_capacity(layers.len());
        let (stiffness, _) = equilibrated_with_profiles(
            &layers,
            bending_axis_rad,
            invalid,
            &mut || Ok(()),
            &mut |polynomial| {
                if polynomial.iter().flatten().any(|x| !x.is_finite()) {
                    return Err(invalid("unrepresentable equilibrated shear profile"));
                }
                polynomials.push(polynomial);
                Ok(())
            },
        )?;
        // Both paths call the SAME numerical kernel, not independently fitted
        // or recomputed stiffness formulas. No tolerance can mask stale H.
        if stiffness != *section.shear_pa_m() {
            return Err(invalid(
                "shear recovery disagrees with its original section stiffness",
            ));
        }
        let mut profiles = Vec::with_capacity(layers.len());
        for ((layer, original), traction) in layers.iter().zip(section.layers()).zip(polynomials) {
            let l00 = layer.g[0].sqrt();
            let l10 = layer.g[2] / l00;
            let l11 = (layer.g[3] - l10 * l10).sqrt();
            if !l00.is_finite() || l00 <= 0.0 || !l10.is_finite() || !l11.is_finite() || l11 <= 0.0
            {
                return Err(invalid("unrepresentable original ply shear compliance"));
            }
            // Use the original laminate's rotation convention.
            let cos = fs_math::det::cos(original.orientation_rad);
            let sin = fs_math::det::sin(original.orientation_rad);
            profiles.push(LayerProfile {
                traction,
                lower: [l00, l10, l11],
                sin,
                cos,
            });
        }
        Ok(Self { section, profiles })
    }

    /// Exact retained section, including original spatial layer/source order.
    pub fn section(&self) -> &LaminateSection {
        &self.section
    }

    /// Recover one ply point from the actual section resultants [Qx,Qy] [N/m].
    /// `fraction` ranges inclusively from that ply's bottom (0) to top (1).
    /// At a bonded interface both plies have the same stress to roundoff,
    /// while their strains may jump according to their different compliance.
    pub fn evaluate(
        &self,
        layer: usize,
        fraction: f64,
        resultant_n_m: [f64; 2],
    ) -> Result<LayerShearResponse, LaminateError> {
        let profile = self
            .profiles
            .get(layer)
            .ok_or(invalid("shear recovery layer out of range"))?;
        if !fraction.is_finite()
            || !(0.0..=1.0).contains(&fraction)
            || resultant_n_m.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "finite resultants and a point inside the original ply required",
            ));
        }
        let matrix: [f64; 4] = core::array::from_fn(|i| {
            (profile.traction[2][i] * fraction + profile.traction[1][i]) * fraction
                + profile.traction[0][i]
        });
        let q = resultant_n_m.map(|x| x / self.section.total_thickness_m());
        let stress_pa = [
            matrix[0] * q[0] + matrix[1] * q[1],
            matrix[2] * q[0] + matrix[3] * q[1],
        ];
        let [l00, l10, l11] = profile.lower;
        let w0 = stress_pa[0] / l00;
        let w1 = (stress_pa[1] - l10 * w0) / l11;
        let gamma1 = w1 / l11;
        let engineering_strain = [(w0 - l10 * gamma1) / l00, gamma1];
        // Sum of squares avoids cancellation in strongly coupled shear work.
        let energy_density_j_m3 = 0.5 * w0 * w0 + 0.5 * w1 * w1;
        let rotate = |v: [f64; 2]| {
            [
                profile.cos * v[0] + profile.sin * v[1],
                -profile.sin * v[0] + profile.cos * v[1],
            ]
        };
        let material_stress_pa = rotate(stress_pa);
        let material_engineering_strain = rotate(engineering_strain);
        let geometry = &self.section.layer_integrals()[layer];
        let z_m = geometry.z_bottom_m + fraction * self.section.layers()[layer].thickness_m;
        if stress_pa
            .iter()
            .chain(&engineering_strain)
            .chain(&material_stress_pa)
            .chain(&material_engineering_strain)
            .any(|x| !x.is_finite())
            || !energy_density_j_m3.is_finite()
            || !z_m.is_finite()
        {
            return Err(invalid("unrepresentable recovered ply shear response"));
        }
        Ok(LayerShearResponse {
            z_m,
            stress_pa,
            engineering_strain,
            material_stress_pa,
            material_engineering_strain,
            energy_density_j_m3,
        })
    }

    /// Recover from the original energy-conjugate section shear coordinate.
    /// It first forms Q = H gamma with the SAME cached coupled stiffness.
    pub fn from_section_strain(
        &self,
        layer: usize,
        fraction: f64,
        engineering_shear: [f64; 2],
    ) -> Result<LayerShearResponse, LaminateError> {
        if engineering_shear.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite section engineering shear required"));
        }
        let h = self.section.shear_pa_m();
        self.evaluate(
            layer,
            fraction,
            [
                h[0] * engineering_shear[0] + h[1] * engineering_shear[1],
                h[2] * engineering_shear[0] + h[3] * engineering_shear[1],
            ],
        )
    }
}
