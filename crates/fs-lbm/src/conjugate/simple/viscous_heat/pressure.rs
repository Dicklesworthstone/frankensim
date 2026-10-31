//! Work-consistent reservoir mouths on the existing native Stokes operator.
//!
//! The reservoir prescribes NORMAL TRACTION -p_res*n and zero TANGENTIAL
//! VELOCITY. Reaction shear performs no boundary work. This is not the legacy
//! static-pressure/zero-gradient opening, nor a zero-shear-traction boundary.
//! A mouth shear half-edge uses gamma=du_n/dt+du_t/dn, including the zero
//! tangential trace at distance h_n/2. Both momentum rows and local heat use
//! the same two-cell harmonic viscosity and physical half-edge volume.
use super::*;

/// Quasisteady source-viscosity flow through opposed reservoir-traction mouths.
///
/// Exactly one opposed pair of Pressure/Outlet planes is required; remaining
/// planes are walls/symmetry. At each mouth the configured pressure denotes
/// normal mechanical traction, while the tangential velocity trace is zero.
/// Pressure work is independently -sum((p_res-p_ref)*Q_out), never a residual
/// converted to heat. This shares native stress and physical dissipation with
/// the closed Stokes owner; it has no transient/open kinetic-energy claim.
pub fn creeping_flow_with_pressure_work_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    viscosity_pa_s: &[f64],
    config: &SimpleConfig,
    controls: ViscousHeatControls,
    gate: &CancelGate,
) -> Result<(FvFlow, ViscousPower), ChtError> {
    creeping_material_heat_with_darcy_tensor(
        domain,
        fluid,
        viscosity_pa_s,
        None,
        None,
        None,
        None,
        config,
        controls,
        true,
        gate,
    )
}

pub(in crate::conjugate) fn admit_pressure_work_heating(
    domain: &VoxelDomain,
    config: &SimpleConfig,
    controls: ViscousHeatControls,
) -> Result<usize, ChtError> {
    let mut axis = None;
    let mut count = 0;
    for (side, face) in config.faces.iter().enumerate() {
        match face {
            FvBoundary::Pressure { .. } | FvBoundary::Outlet => {
                if axis.is_some_and(|a| a != side / 2) {
                    return Err(invalid_power(
                        "pressure-work heat needs one opposed reservoir-mouth pair",
                    ));
                }
                axis = Some(side / 2);
                count += 1;
                finite_value(
                    "pressure_heat.reservoir_pressure",
                    face.prescribed_pressure().expect("pressure face"),
                )?;
            }
            FvBoundary::Wall { .. } | FvBoundary::Symmetry => {}
            _ => {
                return Err(invalid_power(
                    "pressure-work heat admits reservoir traction mouths, walls and symmetry only",
                ));
            }
        }
    }
    if count != 2 {
        return Err(invalid_power(
            "pressure-work heat needs one opposed reservoir-mouth pair",
        ));
    }
    let axis = axis.expect("two pressure mouths");
    // Reuse all source, drive, wall, turbulence and numerical admission. Only
    // the two separately derived mouth conditions differ from the closed mode.
    let mut closed = config.clone();
    closed.faces[2 * axis] = FvBoundary::Symmetry;
    closed.faces[2 * axis + 1] = FvBoundary::Symmetry;
    admit_heating(domain, &closed, controls)?;
    Ok(axis)
}

impl Solver<'_> {
    /// Native assembly mirrors a mouth's interior cell. Its actual physical
    /// volume is one HALF cell, so work traversals use weight 1/2 on that row.
    /// Interior faces retain their original two-cell arithmetic exactly.
    pub(in crate::conjugate::simple) fn power_face_cells(
        &self,
        a: usize,
        index: usize,
    ) -> ([usize; 3], [usize; 3], f64) {
        let f = self.comps[a].coords(index);
        let mut minus = f;
        if f[a] > 0 {
            minus[a] -= 1;
        }
        let plus = if f[a] == self.n[a] { minus } else { f };
        let weight = if f[a] == 0 || f[a] == self.n[a] {
            0.5
        } else {
            1.0
        };
        (minus, plus, weight)
    }

    /// The SAME coefficient for both mouth-edge momentum derivatives and heat.
    pub(in crate::conjugate::simple) fn pressure_mouth_mu(
        &self,
        tangent: usize,
        minus: [usize; 3],
        plus: [usize; 3],
    ) -> f64 {
        let wm = self.domain.width(tangent, minus[tangent]);
        let wp = self.domain.width(tangent, plus[tangent]);
        let mm = self.mu_eff(minus);
        let mp = self.mu_eff(plus);
        let minimum = mm.min(mp);
        minimum / ((wm / (wm + wp)) * (minimum / mm) + (wp / (wm + wp)) * (minimum / mp))
    }

    pub(super) fn pressure_mouth_heat(
        &self,
        heat: &mut [f64],
        quadrature: &mut f64,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        let Some(a) = self.pressure_work_axis else {
            return Ok(());
        };
        let domain = self.domain;
        for side in 0..2 {
            let plane = if side == 0 { 0 } else { self.n[a] - 1 };
            let inward = if side == 0 { 1.0 } else { -1.0 };
            for c in 0..domain.cell_count() {
                if c % 256 == 0 {
                    poll(gate)?;
                }
                let plus = domain.coords(c);
                if plus[a] != plane || !domain.is_fluid(c) {
                    continue;
                }
                for d in 0..3 {
                    if d == a || plus[d] == 0 {
                        continue;
                    }
                    let mut minus = plus;
                    minus[d] -= 1;
                    if !domain.is_fluid(cell_of(self.n, minus)) {
                        continue;
                    }
                    let wm = domain.width(d, minus[d]);
                    let wp = domain.width(d, plus[d]);
                    let ld = 0.5 * (wm + wp);
                    let la = 0.5 * domain.width(a, plus[a]);
                    let normal = (self.vel[a][self.cell_face(a, plus, side == 1)]
                        - self.vel[a][self.cell_face(a, minus, side == 1)])
                        / ld;
                    let tangent = self.vel[d][self.cell_face(d, minus, true)];
                    let gamma = normal + inward * tangent / la;
                    let mu = self.pressure_mouth_mu(d, minus, plus);
                    let power = finite_value(
                        "pressure_heat.mouth_shear",
                        mu * la * ld * domain.width(3 - a - d, plus[3 - a - d]) * gamma * gamma,
                    )?;
                    *quadrature = finite_value("pressure_heat.shear_sum", *quadrature + power)?;
                    let minimum = self.mu_eff(minus).min(self.mu_eff(plus));
                    let weights = [
                        wm / (wm + wp) * minimum / self.mu_eff(minus),
                        wp / (wm + wp) * minimum / self.mu_eff(plus),
                    ];
                    let total = weights[0] + weights[1];
                    for (at, weight) in [(minus, weights[0]), (plus, weights[1])] {
                        let cell = cell_of(self.n, at);
                        heat[cell] = finite_value(
                            "pressure_heat.local_shear",
                            heat[cell] + power * (weight / total),
                        )?;
                    }
                }
            }
        }
        Ok(())
    }

    pub(super) fn pressure_work(
        &self,
        field: &FlowField,
        gate: &CancelGate,
    ) -> Result<([f64; 6], f64), ChtError> {
        let mut power = [0.0; 6];
        let Some(a) = self.pressure_work_axis else {
            return Ok((power, 0.0));
        };
        for side in 0..2 {
            let face = 2 * a + side;
            let pressure = self.faces[face]
                .prescribed_pressure()
                .expect("admitted pressure mouth")
                - self.pressure_reference_pa;
            for c in 0..self.domain.cell_count() {
                if c % 256 == 0 {
                    poll(gate)?;
                }
                if self.domain.is_fluid(c) && self.domain.neighbor(c, face).is_none() {
                    power[face] = finite_value(
                        "pressure_heat.reservoir_work",
                        power[face] - pressure * field.outward(self.domain, c, face),
                    )?;
                }
            }
        }
        let total = finite_value("pressure_heat.supply", power.iter().sum())?;
        Ok((power, total))
    }
}
