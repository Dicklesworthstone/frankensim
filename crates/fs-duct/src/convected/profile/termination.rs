//! Original conjugate end loads on a finite-length variable-flow profile.
//! `H' = Z_Hm m_out' + H_drive'` uses the normal out of the duct at BOTH
//! ends. Local c/(rho A) references differ; power-normalized characteristics
//! keep the feedback and the independently observed physical work consistent.

use super::{IsentropicProfileExcitation, PreparedIsentropicProfile};
use crate::DuctError;
use crate::convected::termination::{FlowImpedance, FlowTerminalObservation, FlowTerminationPower};
use crate::convected::{bad, finite, positive};
use fs_math::{c64::C64, det};

/// Prepared fixed-size, coherent load/source solve on the original profile.
#[derive(Debug)]
pub struct PreparedProfileTerminations<'prepared, 'flow, 'source> {
    profile: &'prepared PreparedIsentropicProfile<'flow, 'source>,
    impedances: [FlowImpedance; 2],
    source_to_wave: [C64; 4],
    condition: f64,
}

/// Driven field with independent physical terminal observations and power.
#[derive(Debug)]
pub struct TerminatedIsentropicProfile<'prepared, 'flow, 'source> {
    field: IsentropicProfileExcitation<'prepared, 'flow, 'source>,
    terminals: [FlowTerminalObservation; 2],
    power: FlowTerminationPower,
}

impl<'flow, 'source> PreparedIsentropicProfile<'flow, 'source> {
    /// Prepare authored passive H'/m_out' loads [left,right].
    ///
    /// `maximum_condition >= 1` bounds the infinity-norm condition estimate of
    /// the complete 2x2 matrix `I - diag(reflection) S`. Undamped resonances
    /// and excessive cancellation refuse; no damping floor is introduced.
    pub fn prepare_terminations(
        &self,
        impedances: [FlowImpedance; 2],
        maximum_condition: f64,
    ) -> Result<PreparedProfileTerminations<'_, 'flow, 'source>, DuctError> {
        if !maximum_condition.is_finite() || maximum_condition < 1.0 {
            return Err(bad(
                "profile end condition bound must be finite and at least one",
            ));
        }
        let mut reflection = [C64::ZERO; 2];
        let mut gain = [C64::ZERO; 2];
        for i in 0..2 {
            let z = impedances[i].value_j_s_kg2();
            let r = self.reference_impedances[i];
            let scale = z.re.max(z.im.abs()).max(r);
            let scaled_r = r / scale;
            let scaled_z = C64::new(z.re / scale, z.im / scale);
            let denominator = scaled_z + C64::from_re(scaled_r);
            reflection[i] = (scaled_z - C64::from_re(scaled_r)) / denominator;
            gain[i] = (C64::from_re(scaled_r) / denominator)
                .scale(core::f64::consts::FRAC_1_SQRT_2 / det::sqrt(r));
            if !positive(scaled_r)
                || !finite(reflection[i])
                || !finite(gain[i])
                || gain[i].abs() == 0.0
            {
                return Err(bad("profile end reflection/source gain is unrepresentable"));
            }
        }
        let s = self.scattering;
        let a = C64::ONE - reflection[0] * s[0];
        let b = -reflection[0] * s[1];
        let c = -reflection[1] * s[2];
        let d = C64::ONE - reflection[1] * s[3];
        let determinant = a * d - b * c;
        let norm = (a.abs() + b.abs()).max(c.abs() + d.abs());
        let adjugate_norm = (d.abs() + b.abs()).max(c.abs() + a.abs());
        let condition = norm * adjugate_norm / determinant.abs();
        if !condition.is_finite() || condition > maximum_condition || determinant.abs() == 0.0 {
            return Err(DuctError::Singular);
        }
        let source_to_wave = [
            d * gain[0] / determinant,
            -b * gain[1] / determinant,
            -c * gain[0] / determinant,
            a * gain[1] / determinant,
        ];
        if !source_to_wave.into_iter().all(finite) {
            return Err(bad("profile end driven response is unrepresentable"));
        }
        Ok(PreparedProfileTerminations {
            profile: self,
            impedances,
            source_to_wave,
            condition,
        })
    }
}

impl<'prepared, 'flow, 'source> PreparedProfileTerminations<'prepared, 'flow, 'source> {
    /// Original prepared finite-length flow profile.
    pub const fn profile(&self) -> &'prepared PreparedIsentropicProfile<'flow, 'source> {
        self.profile
    }
    /// Original authored loads [left,right].
    pub const fn impedances(&self) -> [FlowImpedance; 2] {
        self.impedances
    }
    /// Infinity-norm condition estimate of the complete boundary solve.
    pub const fn condition_estimate(&self) -> f64 {
        self.condition
    }

    /// Solve coherent in-series source enthalpies [J/kg] without allocation,
    /// spatial integration or mean-state re-evaluation.
    pub fn solve(
        &self,
        source_enthalpy_j_kg: [C64; 2],
    ) -> Result<TerminatedIsentropicProfile<'prepared, 'flow, 'source>, DuctError> {
        let h = source_enthalpy_j_kg;
        if !h.into_iter().all(finite) {
            return Err(bad("profile end drives must be finite"));
        }
        let t = self.source_to_wave;
        let field = self
            .profile
            .excite_power_waves([t[0] * h[0] + t[1] * h[1], t[2] * h[0] + t[3] * h[1]])?;
        let sections = [
            field.at_node(0)?,
            field.at_node(self.profile.nodes.len() - 1)?,
        ];
        let terminals = core::array::from_fn(|i| {
            let enthalpy = sections[i].stagnation_enthalpy_perturbation_j_kg;
            let mass =
                sections[i]
                    .mass_flow_perturbation_kg_s
                    .scale(if i == 0 { -1.0 } else { 1.0 });
            let z = self.impedances[i].value_j_s_kg2();
            FlowTerminalObservation {
                enthalpy_j_kg: enthalpy,
                outward_mass_flow_kg_s: mass,
                source_enthalpy_j_kg: h[i],
                boundary_residual_j_kg: enthalpy - z * mass - h[i],
                absorbed_disturbance_power_w: mass
                    .scale(det::sqrt(z.re) * core::f64::consts::FRAC_1_SQRT_2)
                    .norm_sq(),
                source_work_into_duct_w: -0.5 * (h[i] * mass.conj()).re,
                outward_disturbance_power_w: 0.5 * (enthalpy * mass.conj()).re,
            }
        });
        let source_work =
            terminals[0].source_work_into_duct_w + terminals[1].source_work_into_duct_w;
        let absorbed =
            terminals[0].absorbed_disturbance_power_w + terminals[1].absorbed_disturbance_power_w;
        let power = FlowTerminationPower {
            source_work_into_duct_w: source_work,
            absorbed_disturbance_power_w: absorbed,
            balance_residual_w: source_work - absorbed,
            interior_flux_residual_w: terminals[0].outward_disturbance_power_w
                + terminals[1].outward_disturbance_power_w,
        };
        if terminals.iter().any(|end| {
            !finite(end.boundary_residual_j_kg)
                || ![
                    end.absorbed_disturbance_power_w,
                    end.source_work_into_duct_w,
                    end.outward_disturbance_power_w,
                ]
                .into_iter()
                .all(f64::is_finite)
        }) || ![
            source_work,
            absorbed,
            power.balance_residual_w,
            power.interior_flux_residual_w,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(bad("profile end work is unrepresentable"));
        }
        Ok(TerminatedIsentropicProfile {
            field,
            terminals,
            power,
        })
    }
}

impl<'prepared, 'flow, 'source> TerminatedIsentropicProfile<'prepared, 'flow, 'source> {
    /// Full coherent acoustic field on all retained profile sections.
    pub const fn field(&self) -> &IsentropicProfileExcitation<'prepared, 'flow, 'source> {
        &self.field
    }
    /// Physical outward terminal variables, boundary residuals and work.
    pub const fn terminals(&self) -> &[FlowTerminalObservation; 2] {
        &self.terminals
    }
    /// Independent source/load and interior-flux balances.
    pub const fn power(&self) -> &FlowTerminationPower {
        &self.power
    }
}
