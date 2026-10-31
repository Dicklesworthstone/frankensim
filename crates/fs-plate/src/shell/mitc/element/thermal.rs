//! Thermal force and stress recovery in the original enriched element space.

use super::{PreparedStrainField, stiffness_rule};
use crate::shell::mitc::thermal::{FacetThermalStrain, LaminateThermalField, ThermalFacetResponse};
use crate::shell::mitc::{PlateError, ShearShellSection};

type Matrices = ([f64; 9], [f64; 9], [f64; 9], [[f64; 2]; 2]);

pub(in crate::shell::mitc) enum ThermalLoading {
    Strain(FacetThermalStrain),
    Laminate(LaminateThermalField),
}

pub(in crate::shell::mitc) struct PreparedThermalFacet {
    field: PreparedStrainField,
    matrices: Matrices,
    loading: ThermalLoading,
    bubble_offset: [f64; 2],
    pub(in crate::shell::mitc) load: [f64; 18],
}

fn resultants(
    (a, b, d, shear): &Matrices,
    membrane: [f64; 3],
    gradient: [f64; 3],
    gamma: [f64; 2],
) -> ([f64; 3], [f64; 3], [f64; 2]) {
    let n = core::array::from_fn(|i| {
        (0..3)
            .map(|j| a[3 * i + j] * membrane[j] + b[3 * i + j] * gradient[j])
            .sum()
    });
    let m = core::array::from_fn(|i| {
        (0..3)
            .map(|j| b[3 * i + j] * membrane[j] + d[3 * i + j] * gradient[j])
            .sum()
    });
    let q = core::array::from_fn(|i| (0..2).map(|j| shear[i][j] * gamma[j]).sum());
    (n, m, q)
}

impl PreparedThermalFacet {
    pub(in crate::shell::mitc) fn strain_field(
        &self,
        nodal: &[f64; 18],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<crate::shell::mitc::strain::RealElementStrain, PlateError> {
        let mut q = self.field.expand(nodal);
        for i in 0..2 {
            q[18 + i] += self.bubble_offset[i];
        }
        crate::shell::mitc::strain::real_field_from_expanded(&self.field, &q, cancelled)
    }

    pub(in crate::shell::mitc) fn new(
        field: PreparedStrainField,
        section: &ShearShellSection,
        loading: ThermalLoading,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        let matrices = section.coupled_matrices();
        let mut full = [0.0; 20];
        for (lambda, weight) in stiffness_rule() {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            let (n, m) = match &loading {
                ThermalLoading::Strain(strain) => {
                    let point = strain.at(lambda);
                    let (n, m, _) = resultants(
                        &matrices,
                        point.membrane(),
                        point.gradient_per_m(),
                        [0.0; 2],
                    );
                    (n, m)
                }
                ThermalLoading::Laminate(laminate) => {
                    (laminate.membrane_force_n_m(), laminate.bending_moment_n())
                }
            };
            for (c, force) in full.iter_mut().enumerate() {
                let (bm, bb, _) = field.kin.basis(&field.g, lambda[1], lambda[2], c);
                *force += field.g.area_m2
                    * weight
                    * (0..3).map(|i| bm[i] * n[i] + bb[i] * m[i]).sum::<f64>();
            }
        }
        // q_b = -K_bb^-1 K_bn q_n + K_bb^-1 f_th,b. The same
        // homogeneous map condenses load by T^T; the AFFINE offset must also
        // be retained in stress recovery, including when all nodal q are zero.
        let bubble_offset = core::array::from_fn(|i| {
            field.bubble_compliance[i][0] * full[18] + field.bubble_compliance[i][1] * full[19]
        });
        let load = core::array::from_fn(|i| {
            full[i] + field.bubble[0][i] * full[18] + field.bubble[1][i] * full[19]
        });
        if full
            .iter()
            .chain(&load)
            .chain(&bubble_offset)
            .any(|v| !v.is_finite())
        {
            return Err(PlateError::BadSection {
                what: "nonfinite condensed thermal shell load",
            });
        }
        Ok(Self {
            field,
            matrices,
            loading,
            bubble_offset,
            load,
        })
    }

    pub(in crate::shell::mitc) fn recover(
        &self,
        nodal: &[f64; 18],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermalFacetResponse, PlateError> {
        let mut q = self.field.expand(nodal);
        for i in 0..2 {
            q[18 + i] += self.bubble_offset[i];
        }
        let mut out = ThermalFacetResponse {
            membrane_force_n_m: [0.0; 3],
            bending_moment_n: [0.0; 3],
            transverse_shear_n_m: [0.0; 2],
            strain_energy_j: 0.0,
            material_nodal_force: [0.0; 18],
        };
        for (lambda, weight) in stiffness_rule() {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            let (mut eps, mut curvature, gamma) = self.field.eval(&q, lambda[1], lambda[2]);
            let (n, m, shear, energy) = match &self.loading {
                ThermalLoading::Strain(strain) => {
                    let point = strain.at(lambda);
                    for i in 0..3 {
                        eps[i] -= point.membrane()[i];
                        curvature[i] -= point.gradient_per_m()[i];
                    }
                    let (n, m, shear) = resultants(&self.matrices, eps, curvature, gamma);
                    let energy = 0.5
                        * ((0..3)
                            .map(|i| eps[i] * n[i] + curvature[i] * m[i])
                            .sum::<f64>()
                            + (0..2).map(|i| gamma[i] * shear[i]).sum::<f64>());
                    (n, m, shear, energy)
                }
                ThermalLoading::Laminate(laminate) => {
                    let response = laminate.evaluate(eps, curvature, cancelled)?;
                    let (_, _, shear) = resultants(&self.matrices, [0.0; 3], [0.0; 3], gamma);
                    let energy = response.strain_energy_j_m2
                        + 0.5 * (0..2).map(|i| gamma[i] * shear[i]).sum::<f64>();
                    (
                        response.membrane_force_n_m,
                        response.bending_moment_n,
                        shear,
                        energy,
                    )
                }
            };
            for i in 0..3 {
                out.membrane_force_n_m[i] += weight * n[i];
                out.bending_moment_n[i] += weight * m[i];
            }
            for i in 0..2 {
                out.transverse_shear_n_m[i] += weight * shear[i];
            }
            out.strain_energy_j += weight * self.field.g.area_m2 * energy;
            for c in 0..18 {
                let (bm, bb, _) = self.field.kin.basis(&self.field.g, lambda[1], lambda[2], c);
                let cov = [
                    self.field.constant[0][c] + self.field.twisting[c] * (lambda[2] - 1.0 / 3.0),
                    self.field.constant[1][c] + self.field.twisting[c] * (1.0 / 3.0 - lambda[1]),
                ];
                let work = (0..3).map(|i| bm[i] * n[i] + bb[i] * m[i]).sum::<f64>()
                    + (0..2)
                        .map(|i| {
                            (self.field.g.gradient[1][i] * cov[0]
                                + self.field.g.gradient[2][i] * cov[1])
                                * shear[i]
                        })
                        .sum::<f64>();
                out.material_nodal_force[c] += weight * self.field.g.area_m2 * work;
            }
        }
        if out
            .membrane_force_n_m
            .iter()
            .chain(&out.bending_moment_n)
            .chain(&out.transverse_shear_n_m)
            .chain(&out.material_nodal_force)
            .any(|v| !v.is_finite())
            || !out.strain_energy_j.is_finite()
        {
            return Err(PlateError::BadSection {
                what: "nonfinite thermal shell response",
            });
        }
        Ok(out)
    }
}
