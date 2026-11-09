//! Direct native scattering constraints; no admittance or transmission inverse.
use super::super::physics::{CavityOpeningResponse, field_response};
use super::{C64, TwoPort, VibroError, invalid};
use fs_math::dd::Dd;

pub(in crate::vibroacoustic) struct DistributedOpeningFrequency<'a> {
    pub(in crate::vibroacoustic) structural_modes: usize,
    pub(in crate::vibroacoustic) pressure_projection: &'a [Vec<f64>],
    pub(in crate::vibroacoustic) exterior: &'a [C64],
    pub(in crate::vibroacoustic) ducts: &'a [TwoPort],
    pub(in crate::vibroacoustic) tolerance: f64,
}

impl DistributedOpeningFrequency<'_> {
    pub(in crate::vibroacoustic) fn unknown_count(&self) -> usize {
        2 * self.ducts.len()
    }

    /// The two scattering rows are b-S*a=0, where
    /// a=[(p_inner+R Q_inner)/2,(p_outer-R Q_outer)/2] and
    /// b=[(p_inner-R Q_inner)/2,(p_outer+R Q_outer)/2].
    /// p_inner=H^T a_cavity; p_outer=Z_Qs v+Z_QQ Q_outer.
    /// Only Q_inner enters the cavity weak equation. Every exterior port and
    /// both end flows remain present at pressure nodes and opaque transmission.
    pub(in crate::vibroacoustic) fn add_operator(
        &self,
        omega: f64,
        bulk: f64,
        lambdas: &[f64],
        matrix: &mut [C64],
        keep_s: usize,
        keep_a: usize,
    ) -> Result<(), VibroError> {
        let holes = self.ducts.len();
        let ports = self.structural_modes + holes;
        let inner = keep_s + keep_a;
        let outer = inner + holes;
        let n = outer + holes;
        let rate = C64::new(0.0, -omega);
        for r in 0..keep_s {
            for s in 0..keep_s {
                matrix[r * n + s] = matrix[r * n + s] + rate * self.exterior[r * ports + s];
            }
            for j in 0..holes {
                matrix[r * n + outer + j] = self.exterior[r * ports + self.structural_modes + j];
            }
        }
        for q in 0..keep_a {
            let scale = rate.scale(bulk / lambdas[q]);
            for j in 0..holes {
                matrix[(keep_s + q) * n + inner + j] = scale.scale(self.pressure_projection[q][j]);
            }
        }
        for (j, duct) in self.ducts.iter().enumerate() {
            let s = duct.scattering;
            let r = duct.reference_impedance;
            // Coefficients multiply [p_inner,p_outer,Q_inner,Q_outer].
            let rows = [
                [
                    C64::ONE - s[0],
                    -s[1],
                    -(C64::ONE + s[0]).scale(r),
                    s[1].scale(r),
                ],
                [
                    -s[2],
                    C64::ONE - s[3],
                    -s[2].scale(r),
                    (C64::ONE + s[3]).scale(r),
                ],
            ];
            let source = self.structural_modes + j;
            for (end, c) in rows.iter().enumerate() {
                let row = inner + end * holes + j;
                for mode in 0..keep_s {
                    matrix[row * n + mode] = c[1] * rate * self.exterior[source * ports + mode];
                }
                for q in 0..keep_a {
                    matrix[row * n + keep_s + q] = c[0].scale(self.pressure_projection[q][j]);
                }
                for k in 0..holes {
                    matrix[row * n + outer + k] =
                        c[1] * self.exterior[source * ports + self.structural_modes + k];
                }
                matrix[row * n + inner + j] = c[2];
                matrix[row * n + outer + j] = matrix[row * n + outer + j] + c[3];
            }
        }
        if matrix
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid(
                "distributed opening scattering operator is unrepresentable",
            ));
        }
        Ok(())
    }

    pub(in crate::vibroacoustic) fn response(
        &self,
        omega: f64,
        displacement: &[C64],
        pressure: &[C64],
        flow: &[C64],
    ) -> Result<CavityOpeningResponse, VibroError> {
        let holes = self.ducts.len();
        let (inner, outer) = flow.split_at(holes);
        let mut result = field_response(
            self.structural_modes,
            self.pressure_projection,
            self.exterior,
            omega,
            displacement,
            pressure,
            inner,
            outer,
        )?;
        let mut heat = Dd::ZERO;
        let dd = Dd::from_f64;
        for (j, duct) in self.ducts.iter().enumerate() {
            let p = [result.cavity_pressure_pa[j], result.exterior_pressure_pa[j]];
            let q = [inner[j], -outer[j]]; // Native TwoPort flows point into BOTH faces.
            let r = duct.reference_impedance;
            let incoming = core::array::from_fn(|end| {
                let rq = q[end].scale(r);
                C64::new(
                    f64::midpoint(p[end].re, rq.re),
                    f64::midpoint(p[end].im, rq.im),
                )
            });
            let original = duct
                .power_for_incident_pressure(incoming)
                .map_err(|_| invalid("native two-sided opening observation refused"))?;
            let mut scale = 0.0_f64;
            let mut error = 0.0_f64;
            let mut face_work = Dd::ZERO;
            for end in 0..2 {
                scale = scale
                    .max(p[end].abs())
                    .max(original.pressure_pa[end].abs())
                    .max(q[end].abs() * r)
                    .max(original.inward_volume_flow_m3_s[end].abs() * r);
                error = error
                    .max((p[end] - original.pressure_pa[end]).abs())
                    .max((q[end] - original.inward_volume_flow_m3_s[end]).abs() * r);
                face_work = face_work
                    + dd(0.5) * (dd(p[end].re) * dd(q[end].re) + dd(p[end].im) * dd(q[end].im));
            }
            if !scale.is_finite() || !error.is_finite() || error > self.tolerance * scale {
                return Err(invalid(
                    "distributed opening original scattering faces failed residual admission",
                ));
            }
            let physical_work = face_work.to_f64();
            let work_scale = original
                .incident_power_w
                .max(original.outgoing_power_w)
                .max(physical_work.abs());
            if !physical_work.is_finite()
                || (physical_work - original.absorbed_power_w).abs() > self.tolerance * work_scale
                || original.absorbed_power_w < -self.tolerance * work_scale
            {
                return Err(invalid(
                    "distributed opening work is inconsistent or materially active",
                ));
            }
            heat = heat + dd(original.absorbed_power_w);
            result.distributed_duct_power.push(original);
        }
        result.neck_dissipation_w = heat.to_f64();
        if !result.neck_dissipation_w.is_finite() {
            return Err(invalid("distributed opening absorption is unrepresentable"));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
