//! Reciprocal pressure/flow constraints for finite ducts internal to the cavity bank.
use super::{C64, CavityDuctConnection, TwoPort, VibroError, invalid};
use crate::vibroacoustic::openings::observe_duct;
use fs_math::dd::Dd;

/// Native two-ended duct fields and physical work in original connection order.
#[derive(Debug, Clone)]
pub struct InternalCavityDuctResponse {
    /// Positive from the left compartment into the duct [m^3/s].
    pub left_volume_flow_m3_s: Vec<C64>,
    /// Positive from the duct into the right compartment [m^3/s].
    pub right_volume_flow_m3_s: Vec<C64>,
    pub left_pressure_pa: Vec<C64>,
    pub right_pressure_pa: Vec<C64>,
    /// Original native face fields; right-face inward flow has opposite sign.
    pub duct_power: Vec<fs_duct::compose::TwoPortPower>,
    /// Signed work leaving each physical compartment through INTERNAL ducts.
    /// Transfer to another compartment is not another exterior/heat debit.
    pub compartment_outward_power_w: Vec<f64>,
    /// Sum of physical compartment work, checked through original face fields.
    pub net_cavity_work_w: f64,
    /// Native internal acoustic absorption, debited once in the full balance.
    /// Signed roundoff is retained, as on the exterior opening-duct route.
    pub dissipation_w: f64,
}
pub(in crate::vibroacoustic) struct InternalDuctFrequency<'a> {
    pub(super) left_projection: &'a [Vec<f64>],
    pub(super) right_projection: &'a [Vec<f64>],
    pub(super) connections: &'a [CavityDuctConnection],
    pub(super) compartments: usize,
    pub(super) ducts: &'a [TwoPort],
    pub(super) tolerance: f64,
}
impl InternalDuctFrequency<'_> {
    pub(in crate::vibroacoustic) fn unknown_count(&self) -> usize {
        2 * self.ducts.len()
    }
    #[allow(clippy::too_many_arguments)]
    pub(in crate::vibroacoustic) fn add_operator(
        &self,
        omega: f64,
        bulk: f64,
        lambdas: &[f64],
        matrix: &mut [C64],
        keep_s: usize,
        keep_a: usize,
        offset: usize,
        n: usize,
    ) -> Result<(), VibroError> {
        let holes = self.ducts.len();
        if n.checked_mul(n) != Some(matrix.len())
            || offset < keep_s + keep_a
            || offset + 2 * holes > n
        {
            return Err(invalid(
                "internal duct constraints require the complete simultaneous layout",
            ));
        }
        let right = offset + holes;
        let rate = C64::new(0.0, -omega);
        for q in 0..keep_a {
            let scale = rate.scale(bulk / lambdas[q]);
            for j in 0..holes {
                // Outward cavity flow: +Q_left at left, -Q_right at right.
                let at = (keep_s + q) * n;
                matrix[at + offset + j] =
                    matrix[at + offset + j] + scale.scale(self.left_projection[q][j]);
                matrix[at + right + j] =
                    matrix[at + right + j] - scale.scale(self.right_projection[q][j]);
            }
        }
        for (j, duct) in self.ducts.iter().enumerate() {
            let s = duct.scattering;
            let r = duct.reference_impedance;
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
            for (end, c) in rows.iter().enumerate() {
                let row = offset + end * holes + j;
                for q in 0..keep_a {
                    matrix[row * n + keep_s + q] = c[0].scale(self.left_projection[q][j])
                        + c[1].scale(self.right_projection[q][j]);
                }
                matrix[row * n + offset + j] = c[2];
                matrix[row * n + right + j] = c[3];
            }
        }
        if matrix
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid("internal duct operator is unrepresentable"));
        }
        Ok(())
    }
    pub(in crate::vibroacoustic) fn response(
        &self,
        pressure: &[C64],
        flow: &[C64],
    ) -> Result<InternalCavityDuctResponse, VibroError> {
        let (left, right) = flow.split_at(self.ducts.len());
        let dd = Dd::from_f64;
        let mut compartment = vec![Dd::ZERO; self.compartments];
        let mut heat = Dd::ZERO;
        let mut lp = Vec::with_capacity(self.ducts.len());
        let mut rp = Vec::with_capacity(self.ducts.len());
        let mut observations = Vec::with_capacity(self.ducts.len());
        for (j, (duct, link)) in self.ducts.iter().zip(self.connections).enumerate() {
            let a = project(self.left_projection, j, pressure);
            let b = project(self.right_projection, j, pressure);
            let original = observe_duct(duct, [a, b], [left[j], -right[j]], self.tolerance)?;
            compartment[link.left_compartment] = compartment[link.left_compartment]
                + dd(0.5) * (dd(a.re) * dd(left[j].re) + dd(a.im) * dd(left[j].im));
            compartment[link.right_compartment] = compartment[link.right_compartment]
                - dd(0.5) * (dd(b.re) * dd(right[j].re) + dd(b.im) * dd(right[j].im));
            heat = heat + dd(original.absorbed_power_w);
            lp.push(a);
            rp.push(b);
            observations.push(original);
        }
        let net = compartment
            .iter()
            .fold(Dd::ZERO, |sum, &x| sum + x)
            .to_f64();
        let work: Vec<_> = compartment.into_iter().map(Dd::to_f64).collect();
        let heat = heat.to_f64();
        if !net.is_finite() || !heat.is_finite() || work.iter().any(|x| !x.is_finite()) {
            return Err(invalid("internal cavity work is unrepresentable"));
        }
        Ok(InternalCavityDuctResponse {
            left_volume_flow_m3_s: left.to_vec(),
            right_volume_flow_m3_s: right.to_vec(),
            left_pressure_pa: lp,
            right_pressure_pa: rp,
            duct_power: observations,
            compartment_outward_power_w: work,
            net_cavity_work_w: net,
            dissipation_w: heat,
        })
    }
}
fn project(map: &[Vec<f64>], j: usize, pressure: &[C64]) -> C64 {
    let dd = Dd::from_f64;
    let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
    for (q, p) in pressure.iter().enumerate() {
        re = re + dd(p.re) * dd(map[q][j]);
        im = im + dd(p.im) * dd(map[q][j]);
    }
    C64::new(re.to_f64(), im.to_f64())
}

pub(in crate::vibroacoustic) fn flow_delta(
    full: Option<&InternalCavityDuctResponse>,
    reduced: Option<&InternalCavityDuctResponse>,
) -> Result<f64, VibroError> {
    let Some(full) = full else {
        return Ok(0.0);
    };
    let reduced =
        reduced.ok_or_else(|| invalid("internal duct convergence lost its original graph"))?;
    let left = super::super::openings::relative_flow_delta(
        &full.left_volume_flow_m3_s,
        &reduced.left_volume_flow_m3_s,
    )?;
    let right = super::super::openings::relative_flow_delta(
        &full.right_volume_flow_m3_s,
        &reduced.right_volume_flow_m3_s,
    )?;
    Ok(left.max(right))
}

#[cfg(test)]
mod tests;
