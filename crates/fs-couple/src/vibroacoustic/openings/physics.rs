//! Original mixed pressure/displacement/flow equations and independent work.
use super::{C64, VibroError, invalid};
use fs_math::dd::Dd;

/// Solved physical opening flows and the matching complete exterior excitation.
#[derive(Clone, Debug)]
pub struct CavityOpeningResponse {
    /// Positive OUT into the exterior, one complex peak mouth flow [m^3/s].
    /// Equal to cavity flow for the compact owner; generally different for a
    /// finite compressible duct.
    pub volume_flow_m3_s: Vec<C64>,
    /// Positive FROM the cavity INTO each neck/duct [m^3/s]. This is the
    /// flow entering the original cavity weak equations and source-work audit.
    pub cavity_volume_flow_m3_s: Vec<C64>,
    /// Native two-sided scattering observations for distributed ducts, in
    /// opening order. Empty for compact one-flow necks. Right-port inward
    /// flow is the NEGATIVE of the exterior mouth flow above.
    pub distributed_duct_power: Vec<fs_duct::compose::TwoPortPower>,
    /// Work-conjugate inner pressures, H^T a [Pa].
    pub cavity_pressure_pa: Vec<C64>,
    /// Outer pressures from the complete exterior impedance [Pa].
    pub exterior_pressure_pa: Vec<C64>,
    /// [Full original structural modal velocities, outward opening flows].
    /// Truncated structural coordinates are zero at their ORIGINAL indices.
    /// These amplitudes can drive the same native BEM velocity-port solutions
    /// for coherent pressure/particle-velocity observation; no source re-solve
    /// or replacement radiation model is required.
    pub radiation_port_velocities: Vec<C64>,
    /// Signed pressure work leaving the cavity through these openings [W].
    pub cavity_outward_power_w: f64,
    /// Internal acoustic absorption [W], excluding the exterior load.
    /// Compact necks measure viscous channel heat. Distributed ducts retain
    /// native two-sided viscous/thermal/molecular absorption, with signed
    /// floating-point roundoff and no repair or additional power debit.
    pub neck_dissipation_w: f64,
    /// Complete joint exterior work, including wall/mouth interference [W].
    /// This is the same quantity reported as PowerBreakdown::radiated,
    /// not an additional debit or a sum of independent source powers.
    pub exterior_power_w: f64,
}

pub(in crate::vibroacoustic) struct OpeningFrequency<'a> {
    pub(in crate::vibroacoustic) structural_modes: usize,
    pub(in crate::vibroacoustic) pressure_projection: &'a [Vec<f64>],
    pub(in crate::vibroacoustic) exterior: &'a [C64],
    pub(in crate::vibroacoustic) neck: &'a [C64],
}

impl OpeningFrequency<'_> {
    pub(in crate::vibroacoustic) fn count(&self) -> usize {
        self.neck.len()
    }

    /// Retain Q as an unknown even at a neck or cavity resonance.
    /// Rows: D_s b-C a+Z_ss v+Z_sQ Q=F;
    /// D_a a-S omega^2 C^T b-i omega S H Q=-i omega S Q_in;
    /// Z_Qs v+(Z_QQ+Z_neck)Q-H^T a=0, with v=-i omega b.
    pub(in crate::vibroacoustic) fn add_operator(
        &self,
        omega: f64,
        bulk: f64,
        lambdas: &[f64],
        matrix: &mut [C64],
        keep_s: usize,
        keep_a: usize,
    ) -> Result<(), VibroError> {
        let holes = self.count();
        let ports = self.structural_modes + holes;
        let offset = keep_s + keep_a;
        let n = offset + holes;
        let rate = C64::new(0.0, -omega);
        for r in 0..keep_s {
            for s in 0..keep_s {
                matrix[r * n + s] = matrix[r * n + s] + rate * self.exterior[r * ports + s];
            }
            for j in 0..holes {
                matrix[r * n + offset + j] = self.exterior[r * ports + self.structural_modes + j];
            }
        }
        for q in 0..keep_a {
            let scale = rate.scale(bulk / lambdas[q]);
            for j in 0..holes {
                matrix[(keep_s + q) * n + offset + j] = scale.scale(self.pressure_projection[q][j]);
            }
        }
        for j in 0..holes {
            let row = offset + j;
            let source = self.structural_modes + j;
            for s in 0..keep_s {
                matrix[row * n + s] = rate * self.exterior[source * ports + s];
            }
            for q in 0..keep_a {
                matrix[row * n + keep_s + q] = C64::from_re(-self.pressure_projection[q][j]);
            }
            for k in 0..holes {
                matrix[row * n + offset + k] =
                    self.exterior[source * ports + self.structural_modes + k];
            }
            matrix[row * n + row] = matrix[row * n + row] + self.neck[j];
        }
        if matrix
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid("joint opening operator is unrepresentable"));
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
        let dd = Dd::from_f64;
        let mut heat = Dd::ZERO;
        for (j, &q) in flow.iter().enumerate() {
            heat =
                heat + dd(0.5) * dd(self.neck[j].re) * (dd(q.re) * dd(q.re) + dd(q.im) * dd(q.im));
        }
        let heat = heat.to_f64();
        if !heat.is_finite() {
            return Err(invalid("compact opening heat is unrepresentable"));
        }
        let mut result = field_response(
            self.structural_modes,
            self.pressure_projection,
            self.exterior,
            omega,
            displacement,
            pressure,
            flow,
            flow,
        )?;
        result.neck_dissipation_w = heat;
        Ok(result)
    }
}

/// One physical pressure/exterior-work owner for compact and distributed flow.
#[allow(clippy::too_many_arguments)]
pub(super) fn field_response(
    structural_modes: usize,
    pressure_projection: &[Vec<f64>],
    exterior_impedance: &[C64],
    omega: f64,
    displacement: &[C64],
    pressure: &[C64],
    inner_flow: &[C64],
    outer_flow: &[C64],
) -> Result<CavityOpeningResponse, VibroError> {
    let holes = outer_flow.len();
    let ports = structural_modes + holes;
    let mut velocity = vec![C64::ZERO; ports];
    for (r, &b) in displacement.iter().enumerate() {
        velocity[r] = C64::new(omega * b.im, -omega * b.re);
    }
    velocity[structural_modes..].copy_from_slice(outer_flow);
    let mut inner = Vec::with_capacity(holes);
    let mut outer = Vec::with_capacity(holes);
    let mut outward_work = Dd::ZERO;
    let dd = Dd::from_f64;
    for (j, &q) in inner_flow.iter().enumerate() {
        let p = products(
            pressure
                .iter()
                .enumerate()
                .map(|(mode, &a)| (a, C64::from_re(pressure_projection[mode][j]))),
        );
        let row = (structural_modes + j) * ports;
        let exterior = products(
            velocity
                .iter()
                .enumerate()
                .map(|(r, &v)| (exterior_impedance[row + r], v)),
        );
        outward_work = outward_work + dd(0.5) * (dd(p.re) * dd(q.re) + dd(p.im) * dd(q.im));
        inner.push(p);
        outer.push(exterior);
    }
    let cavity_outward_power_w = outward_work.to_f64();
    let exterior_power_w =
        super::super::radiation::flow_power(exterior_impedance, ports, &velocity);
    if velocity
        .iter()
        .chain(&inner)
        .chain(&outer)
        .chain(inner_flow)
        .any(|z| !z.re.is_finite() || !z.im.is_finite())
        || !cavity_outward_power_w.is_finite()
        || !exterior_power_w.is_finite()
    {
        return Err(invalid(
            "joint opening field or physical work is unrepresentable",
        ));
    }
    Ok(CavityOpeningResponse {
        volume_flow_m3_s: outer_flow.to_vec(),
        cavity_volume_flow_m3_s: inner_flow.to_vec(),
        distributed_duct_power: Vec::new(),
        cavity_pressure_pa: inner,
        exterior_pressure_pa: outer,
        radiation_port_velocities: velocity,
        cavity_outward_power_w,
        neck_dissipation_w: 0.0,
        exterior_power_w,
    })
}

fn products(values: impl Iterator<Item = (C64, C64)>) -> C64 {
    let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
    let dd = Dd::from_f64;
    for (a, b) in values {
        re = re + dd(a.re) * dd(b.re) - dd(a.im) * dd(b.im);
        im = im + dd(a.re) * dd(b.im) + dd(a.im) * dd(b.re);
    }
    C64::new(re.to_f64(), im.to_f64())
}

/// Observe the new physical unknowns separately from the reference modal
/// energy metric. Full and reduced models retain the SAME opening basis.
pub(in crate::vibroacoustic) fn flow_delta(
    full: Option<&CavityOpeningResponse>,
    reduced: Option<&CavityOpeningResponse>,
) -> Result<f64, VibroError> {
    let Some(full) = full else {
        return Ok(0.0);
    };
    let reduced =
        reduced.ok_or_else(|| invalid("opening convergence lost its original flow basis"))?;
    let mut worst = 0.0_f64;
    for (full, reduced) in [
        (&full.volume_flow_m3_s, &reduced.volume_flow_m3_s),
        (
            &full.cavity_volume_flow_m3_s,
            &reduced.cavity_volume_flow_m3_s,
        ),
    ] {
        if full.len() != reduced.len() {
            return Err(invalid(
                "opening convergence changed its original flow basis",
            ));
        }
        let (mut norm, mut error) = (0.0, 0.0);
        for (&a, &b) in full.iter().zip(reduced) {
            norm = fs_math::det::hypot(norm, a.abs());
            error = fs_math::det::hypot(error, (a - b).abs());
        }
        if !norm.is_finite() || !error.is_finite() {
            return Err(invalid("opening convergence flow norm is unrepresentable"));
        }
        let delta = if norm == 0.0 {
            if error == 0.0 { 0.0 } else { f64::INFINITY }
        } else {
            error / norm
        };
        worst = worst.max(delta);
    }
    Ok(worst)
}

#[cfg(test)]
mod tests;
