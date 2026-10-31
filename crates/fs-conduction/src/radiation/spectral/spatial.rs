//! Positive adaptive surface quadrature for nonpolynomial Planck-band loads.

use super::super::pointwise::TRIANGLE_DEGREE_FIVE;
use super::{BlackbodyBand, SpectralPointwiseControls, SpectralSpatialReport, finite, invalid};
use crate::ConductionError;
use fs_exec::Cx;

#[derive(Clone, Copy, Default)]
pub(in crate::radiation) struct FaceIntegral {
    /// Integral N_i [E_band(T_h)-E_band(T_reference)] dA.
    pub load: [f64; 3],
    /// Integral N_i N_j E_band'(T_h) dA, before surface emissivity.
    pub tangent: [[f64; 3]; 3],
}

impl FaceIntegral {
    fn add(&mut self, other: &Self) -> Result<(), ConductionError> {
        for i in 0..3 {
            self.load[i] = finite(self.load[i] + other.load[i])?;
            for j in 0..3 {
                self.tangent[i][j] = finite(self.tangent[i][j] + other.tangent[i][j])?;
            }
        }
        Ok(())
    }
}

pub(in crate::radiation) struct Work {
    pub controls: SpectralPointwiseControls,
    pub report: SpectralSpatialReport,
}

impl Work {
    pub fn new(controls: SpectralPointwiseControls) -> Self {
        Self {
            controls,
            report: SpectralSpatialReport::default(),
        }
    }
}

struct Face<'a, 'scope> {
    cx: &'a Cx<'scope>,
    band: BlackbodyBand,
    temperatures: [f64; 3],
    reference: f64,
}

pub(in crate::radiation) fn integrate(
    cx: &Cx<'_>,
    band: BlackbodyBand,
    temperatures: [f64; 3],
    reference: f64,
    area: f64,
    budget_fraction: f64,
    work: &mut Work,
) -> Result<FaceIntegral, ConductionError> {
    let face = Face {
        cx,
        band,
        temperatures,
        reference,
    };
    let triangle = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let coarse = face.rule(triangle, area, work)?;
    face.refine(triangle, area, budget_fraction, 0, coarse, work)
}

impl Face<'_, '_> {
    fn rule(
        &self,
        triangle: [[f64; 3]; 3],
        area: f64,
        work: &mut Work,
    ) -> Result<FaceIntegral, ConductionError> {
        self.cx
            .checkpoint()
            .map_err(|_| ConductionError::Cancelled {
                stage: "spectral-surface-quadrature",
                at: work.report.planck_evaluations,
            })?;
        let mut out = FaceIntegral::default();
        for (local, weight) in TRIANGLE_DEGREE_FIVE {
            if work.report.planck_evaluations >= work.controls.maximum_spatial_evaluations {
                return Err(invalid("spectral spatial evaluation budget exhausted"));
            }
            work.report.planck_evaluations += 1;
            let bary = std::array::from_fn::<_, 3, _>(|i| {
                local[0].mul_add(
                    triangle[0][i],
                    local[1].mul_add(triangle[1][i], local[2] * triangle[2][i]),
                )
            });
            let t = finite(
                self.temperatures[0]
                    + (self.temperatures[1] - self.temperatures[0]).mul_add(
                        bary[1],
                        (self.temperatures[2] - self.temperatures[0]) * bary[2],
                    ),
            )?;
            let measure = finite(area * weight)?;
            let emission = finite(measure * self.band.difference(t, self.reference)?)?;
            let slope = finite(measure * self.band.evaluate(t)?.temperature_derivative_w_m2_k)?;
            for i in 0..3 {
                out.load[i] = finite(bary[i].mul_add(emission, out.load[i]))?;
                for j in 0..=i {
                    out.tangent[i][j] =
                        finite((bary[i] * bary[j]).mul_add(slope, out.tangent[i][j]))?;
                }
            }
        }
        for i in 0..3 {
            for j in 0..i {
                out.tangent[j][i] = out.tangent[i][j];
            }
        }
        Ok(out)
    }

    fn refine(
        &self,
        triangle: [[f64; 3]; 3],
        area: f64,
        budget_fraction: f64,
        depth: usize,
        coarse: FaceIntegral,
        work: &mut Work,
    ) -> Result<FaceIntegral, ConductionError> {
        let midpoint = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|i| 0.5 * a[i] + 0.5 * b[i]);
        let [a, b, c] = triangle;
        let (ab, bc, ca) = (midpoint(a, b), midpoint(b, c), midpoint(c, a));
        let children = [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]];
        let mut parts = [FaceIntegral::default(); 4];
        let mut fine = FaceIntegral::default();
        for (index, &child) in children.iter().enumerate() {
            parts[index] = self.rule(child, 0.25 * area, work)?;
            fine.add(&parts[index])?;
        }
        let mut power_error = 0.0;
        let mut power_scale = 0.0;
        let mut tangent_error = 0.0;
        let mut tangent_scale = 0.0;
        for i in 0..3 {
            power_error = finite(power_error + (fine.load[i] - coarse.load[i]).abs())?;
            power_scale = finite(power_scale + fine.load[i].abs())?;
            for j in 0..3 {
                tangent_error =
                    finite(tangent_error + (fine.tangent[i][j] - coarse.tangent[i][j]).abs())?;
                tangent_scale = finite(tangent_scale + fine.tangent[i][j].abs())?;
            }
        }
        let power_gate = finite(
            work.controls.absolute_spatial_power_w * budget_fraction
                + work.controls.relative_spatial_tolerance * power_scale,
        )?;
        let tangent_gate = finite(
            work.controls.absolute_spatial_tangent_w_k * budget_fraction
                + work.controls.relative_spatial_tolerance * tangent_scale,
        )?;
        if power_error <= power_gate && tangent_error <= tangent_gate {
            work.report.accepted_subtriangles += 4;
            work.report.estimated_emission_error_w =
                finite(work.report.estimated_emission_error_w + power_error)?;
            work.report.estimated_tangent_error_w_k =
                finite(work.report.estimated_tangent_error_w_k + tangent_error)?;
            return Ok(fine);
        }
        if depth + 1 >= work.controls.maximum_spatial_depth {
            return Err(invalid(
                "spectral spatial quadrature did not meet its power/tangent estimates before the depth limit",
            ));
        }
        let mut result = FaceIntegral::default();
        for (index, child) in children.into_iter().enumerate() {
            result.add(&self.refine(
                child,
                0.25 * area,
                0.25 * budget_fraction,
                depth + 1,
                parts[index],
                work,
            )?)?;
        }
        Ok(result)
    }
}
