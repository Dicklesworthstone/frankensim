//! Lateral moments of the accepted distributed aerodynamic solution.
//!
//! Reuse the circulation and washed incident velocities from the actual
//! propeller/wing solve. There is no second aerodynamic solve and no gust
//! coefficient or random kick. The existing reference-flow induced/parasite
//! drag acts on the centerline and supplies no differential yaw at this tier.

use crate::Refusal;
use crate::aircraft::{ForceBuildup, OpenLoopDesign};

/// Body-FRD roll and yaw moments about the aircraft CG, in N m.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LateralMomentsNm {
    /// Positive about +x: right wing down.
    pub roll_nm: f64,
    /// Positive about +z: nose right.
    pub yaw_nm: f64,
}

impl LateralMomentsNm {
    fn add_force(&mut self, r: [f64; 3], f: [f64; 3]) {
        let roll_a = r[1] * f[2];
        let roll_b = r[2] * f[1];
        let yaw_a = r[0] * f[1];
        let yaw_b = r[1] * f[0];
        self.roll_nm += roll_a - roll_b;
        self.yaw_nm += yaw_a - yaw_b;
    }

    pub(super) fn is_finite(self) -> bool {
        self.roll_nm.is_finite() && self.yaw_nm.is_finite()
    }
}

/// Physical lateral moment contributions from a distributed load solve.
/// The established reduced control/damping moments are NOT included here.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DistributedLateralLoads {
    /// Kutta-Joukowsky panel forces at the bound-vortex midpoints.
    pub circulation: LateralMomentsNm,
    /// Local chord-plane profile drag at each strip reference line.
    pub profile_drag: LateralMomentsNm,
    /// Separate axial propeller thrusts at their original disk centers.
    pub propellers: LateralMomentsNm,
}

impl DistributedLateralLoads {
    /// Net moments, summed in the same declared contribution order.
    #[must_use]
    pub fn net(self) -> LateralMomentsNm {
        LateralMomentsNm {
            roll_nm: (self.circulation.roll_nm + self.profile_drag.roll_nm)
                + self.propellers.roll_nm,
            yaw_nm: (self.circulation.yaw_nm + self.profile_drag.yaw_nm)
                + self.propellers.yaw_nm,
        }
    }
}

impl OpenLoopDesign {
    /// Assemble roll/yaw moments from an already accepted distributed solve.
    ///
    /// Supply the same design, effective canard setting and density as to
    /// `force_buildup_with_incident`. Accepted local velocities already contain
    /// prop wash and rigid pitch motion; neither is added again. Layout work
    /// is linear in the panel count, without another influence matrix or
    /// nonlinear iteration. No section gust memory, differential induced-drag
    /// model, rotor reaction torque or crossflow force is invented here.
    ///
    /// # Errors
    /// `distributed-loads-invalid` on missing/inconsistent/nonfinite accepted
    /// fields, or a non-representable resultant. Layout refusals pass through.
    pub fn distributed_lateral_loads(
        &self,
        loads: &ForceBuildup,
        delta_canard_rad: f64,
        rho_kg_m3: f64,
    ) -> Result<DistributedLateralLoads, Refusal> {
        if !rho_kg_m3.is_finite() || rho_kg_m3 <= 0.0
            || !delta_canard_rad.is_finite()
            || !self.cd0_wing.is_finite() || self.cd0_wing < 0.0
            || !self.cg_m.iter().all(|v| v.is_finite())
            || !loads.thrust_n.iter().all(|v| v.is_finite())
        {
            return Err(invalid("invalid density, geometry, drag or thrust"));
        }
        let local = loads.coupled.strip_velocity_mps.as_ref()
            .ok_or_else(|| invalid("distributed moments require the accepted local force velocities"))?;
        let (panels, strips) = self.layout(delta_canard_rad)?;
        if loads.coupled.gamma.len() != panels.len() || local.len() != strips.len()
            || !loads.coupled.gamma.iter().all(|v| v.is_finite())
        {
            return Err(invalid("circulation or strip velocities do not match the aircraft layout"));
        }
        let mut result = DistributedLateralLoads::default();
        for (strip, &velocity) in strips.iter().zip(local) {
            let speed2 = velocity[0] * velocity[0] + velocity[2] * velocity[2];
            if !velocity.iter().all(|v| v.is_finite()) || velocity[0] <= 0.0
                || !speed2.is_finite() || speed2 <= 0.0
            {
                return Err(invalid("invalid accepted strip velocity"));
            }
            for &j in &strip.panel_indices {
                let panel = &panels[j];
                let seg = std::array::from_fn::<_, 3, _>(|a| panel.b[a] - panel.a[a]);
                let arm = std::array::from_fn(|a| 0.5 * (panel.a[a] + panel.b[a]) - self.cg_m[a]);
                let fx_a = velocity[1] * seg[2];
                let fx_b = velocity[2] * seg[1];
                let fy_a = velocity[2] * seg[0];
                let fy_b = velocity[0] * seg[2];
                let fz_a = velocity[0] * seg[1];
                let fz_b = velocity[1] * seg[0];
                let scale = rho_kg_m3 * loads.coupled.gamma[j];
                result.circulation.add_force(arm, [
                    scale * (fx_a - fx_b), scale * (fy_a - fy_b), scale * (fz_a - fz_b),
                ]);
            }
            let front = &panels[strip.panel_indices[0]];
            let arm = std::array::from_fn(|a| 0.5 * (front.a[a] + front.b[a]) - self.cg_m[a]);
            let width = (front.b[1] - front.a[1]).abs();
            let drag = 0.5 * rho_kg_m3 * speed2 * strip.chord_m * width * self.cd0_wing;
            let speed = speed2.sqrt();
            result.profile_drag.add_force(arm, [
                -drag * (velocity[0] / speed), 0.0, -drag * (velocity[2] / speed),
            ]);
        }
        for (disk, &thrust) in self.disks.iter().zip(&loads.thrust_n) {
            let arm = std::array::from_fn(|a| disk.center_m[a] - self.cg_m[a]);
            result.propellers.add_force(arm, [thrust, 0.0, 0.0]);
        }
        if ![result.circulation, result.profile_drag, result.propellers, result.net()]
            .iter().all(|v| v.is_finite())
        {
            return Err(invalid("nonfinite distributed lateral moment"));
        }
        Ok(result)
    }
}

fn invalid(message: &str) -> Refusal {
    Refusal {
        code: "distributed-loads-invalid",
        message: message.into(),
        ranked_repairs: vec!["use the original distributed solve, layout, density and accepted force velocities".into()],
    }
}
