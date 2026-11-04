//! Finite-strain membranes from the existing three-dimensional energy owner.
//!
//! Material coordinates are warp, weft and reference normal. A positive
//! Cholesky representative of `C_in_plane = I + 2 E` removes rigid rotation;
//! native plane stress relaxes the entire transverse column of F. Multiplying
//! the reduced potential by the original thickness gives work per reference
//! area. Its exact chain rule supplies second-Piola resultants and their
//! Green-strain tangent, including the stress times second-kinematic-derivative
//! term. The physical membrane element retains its own geometric stiffness.

use super::{
    AREAL_DENSITY, MEMBRANE_PROPERTIES, ResolvedSheet, SheetError, SheetLaw, SheetResponse,
    THICKNESS,
};
use crate::hyper::cards::ResolvedHyperelastic;
use crate::hyper::plane_stress::{PlaneStressResponse, PlaneStressSettings};
use crate::hyper::{FibreInvariants, HyperelasticModel};
use crate::state_point::ResolvedMaterialStatePoint;
use fs_ad::{Dual, Dual64, Real};
use fs_matdb::MaterialCard;

/// Immutable, source-bound hyperelastic membrane at fixed external conditions.
///
/// Reference mass is independent of the bulk law and thickness. The original
/// source bundle and exact model card survive extraction into `MembraneLaw`.
/// This is the native locally stable plane-stress branch, with no wrinkling,
/// bending, residual stress, thermal reference change or incompressibility
/// inferred. Arbitrary fibre directions remain in the reference material frame.
#[derive(Clone, Debug)]
pub struct HyperelasticMembrane {
    material: ResolvedMaterialStatePoint,
    model: ResolvedHyperelastic,
    reference_thickness_m: f64,
    settings: PlaneStressSettings,
}

impl PartialEq for HyperelasticMembrane {
    fn eq(&self, other: &Self) -> bool {
        self.material == other.material
            && self.model == other.model
            && self.reference_thickness_m == other.reference_thickness_m
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

    /// Exact selected model card and the same fixed, typed external context.
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

    fn reduced_response(&self, strain: [f64; 3]) -> Result<PlaneStressResponse, SheetError> {
        let f = admitted_representative(strain)?;
        Ok(self.model.law().plane_stress(&f, self.settings)?)
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

fn admitted_representative(strain: [f64; 3]) -> Result<[f64; 4], SheetError> {
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
            resultants_n_m[i] *= self.reference_thickness_m;
            for j in 0..3 {
                tangent_n_m[3 * i + j] *= self.reference_thickness_m;
            }
        }
        let energy_j_m2 = self.reference_thickness_m * response.energy_density_pa;
        if !energy_j_m2.is_finite()
            || resultants_n_m
                .iter()
                .chain(&tangent_n_m)
                .any(|v| !v.is_finite())
        {
            return Err(SheetError::Invalid(
                "unrepresentable hyperelastic membrane response",
            ));
        }
        Ok(SheetResponse {
            energy_j_m2,
            resultants_n_m,
            tangent_n_m,
        })
    }
}
