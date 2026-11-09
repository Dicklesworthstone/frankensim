//! Passive external RCR storage and optional series inertance, per unit depth.

use crate::ns::boundary::bdf::StokesBdfStep;
use fs_math::dd::Dd;

/// Linear RCR/Windkessel load, optionally with series fluid inertance.
/// All pressures are increments to the original port traction, normalized by
/// ONE fixed reference density. Parameters and compliance reference stay
/// fixed through the supplied histories; the reservoir is held per interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesWindkesselPort {
    /// Endpoint reservoir pressure increment [m²/s²].
    pub reservoir_pressure_m2_s2: f64,
    /// Nonnegative proximal resistance [1/s]. Zero is admitted.
    pub proximal_resistance_s_inv: f64,
    /// Strictly positive distal resistance to the same reservoir [1/s].
    pub distal_resistance_s_inv: f64,
    /// Strictly positive volume-per-depth compliance times density [s²].
    pub compliance_s2: f64,
    /// Nonnegative series inertance divided by density [dimensionless].
    /// Zero gives the ordinary three-element RCR load.
    pub inertance: f64,
}

/// External storage, physical loss and time-discretization accounting for one
/// accepted port. Multiply energies/powers by reference density for J/m, W/m.
/// These terms are outside the resolved fluid and are never added to its heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesWindkesselBalance {
    /// Compliance pressure relative to its reservoir, x [m²/s²].
    pub compliance_pressure_m2_s2: f64,
    /// Stored compliance energy C*x²/2 [m⁴/s²].
    pub compliance_energy_m4_s2: f64,
    /// Stored series-inertance energy L*Q²/2 [m⁴/s²].
    pub inertance_energy_m4_s2: f64,
    /// Proximal resistor loss Rp*Q² [m⁴/s³], nonnegative.
    pub proximal_dissipation_m4_s3: f64,
    /// Distal resistor loss x²/Rd [m⁴/s³], nonnegative.
    pub distal_dissipation_m4_s3: f64,
    /// Reservoir pressure increment times outward Q [m⁴/s³], signed.
    pub reservoir_power_m4_s3: f64,
    /// Added pressure times outward Q [m⁴/s³], entering the external load.
    /// The original traction's work is not included a second time.
    pub port_power_m4_s3: f64,
    /// (E_new - E_previous)/actual_timestep [m⁴/s³].
    pub stored_energy_change_rate_m4_s3: f64,
    /// C*x*D_t(x) + L*Q*D_t(Q) minus the physical energy-change rate
    /// [m⁴/s³]. This is numerical, not resistor heat; it can be negative for
    /// BDF2. No single-state energy monotonicity is claimed for variable steps.
    pub bdf_energy_defect_m4_s3: f64,
    /// Port power minus reservoir power, both losses, energy change and BDF
    /// defect [m⁴/s³]. Computed independently from the accepted endpoint laws.
    pub power_balance_residual_m4_s3: f64,
    /// Signed p - p_res - Rp*Q - L*D_t(Q) - x [m²/s²].
    pub pressure_law_residual_m2_s2: f64,
    /// Signed C*D_t(x) - Q + x/Rd [m²/s].
    pub compliance_flow_residual_m2_s: f64,
}

pub(in crate::ns::boundary) struct ReducedWindkesselPort {
    pub(super) pressure_offset: f64,
    pub(super) resistance: f64,
    load: StokesWindkesselPort,
    step: StokesBdfStep,
    previous_q: f64,
    previous_x: f64,
    history_q: f64,
    history_x: f64,
    capacitor_response: f64,
    capacitor_offset: f64,
}

impl ReducedWindkesselPort {
    /// Q histories come from the original native velocity edge moments.
    pub(super) fn new(
        load: StokesWindkesselPort,
        step: StokesBdfStep,
        previous_q: f64,
        previous_x: f64,
        older: Option<(f64, f64)>,
    ) -> Result<Self, &'static str> {
        if ![
            load.reservoir_pressure_m2_s2,
            load.proximal_resistance_s_inv,
            load.distal_resistance_s_inv,
            load.compliance_s2,
            load.inertance,
            previous_q,
            previous_x,
        ]
        .iter()
        .all(|v| v.is_finite())
            || load.proximal_resistance_s_inv < 0.0
            || load.distal_resistance_s_inv <= 0.0
            || load.compliance_s2 <= 0.0
            || load.inertance < 0.0
            || older.is_some_and(|(q, x)| !q.is_finite() || !x.is_finite())
        {
            return Err("finite histories and passive RCR/inertance parameters required");
        }
        let (old_q, old_x) = older.unwrap_or((previous_q, previous_x));
        let history_q = step.history(previous_q, old_q);
        let history_x = step.history(previous_x, old_x);
        let mass = load.compliance_s2 / step.effective_timestep_s;
        let conductance = load.distal_resistance_s_inv.recip();
        let denominator = mass + conductance;
        let capacitor_response = denominator.recip();
        let inertia = load.inertance / step.effective_timestep_s;
        if ![mass, conductance, denominator, capacitor_response]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
            || ![history_q, history_x, inertia]
                .iter()
                .all(|v| v.is_finite())
            || (load.inertance > 0.0 && inertia == 0.0)
        {
            return Err("unrepresentable RCR/inertance history or endpoint coefficient");
        }
        let d = Dd::from_f64;
        let capacitor_offset = (d(mass / denominator) * d(history_x)).to_f64();
        let pressure_offset = (d(load.reservoir_pressure_m2_s2) + d(capacitor_offset)
            - d(inertia) * d(history_q))
        .to_f64();
        let resistance =
            (d(load.proximal_resistance_s_inv) + d(inertia) + d(capacitor_response)).to_f64();
        if ![capacitor_offset, pressure_offset, resistance]
            .iter()
            .all(|v| v.is_finite())
            || resistance <= 0.0
        {
            return Err("unrepresentable effective hydraulic resistance or pressure");
        }
        Ok(Self {
            pressure_offset,
            resistance,
            load,
            step,
            previous_q,
            previous_x,
            history_q,
            history_x,
            capacitor_response,
            capacitor_offset,
        })
    }

    /// Reconstruct storage from the ACTUAL integrated endpoint Q and evaluate
    /// the original two laws, not just the eliminated effective resistance.
    pub(super) fn report(
        &self,
        q: f64,
        pressure: f64,
    ) -> Result<(StokesWindkesselBalance, f64, f64), &'static str> {
        let d = Dd::from_f64;
        let l = self.load;
        let x = (d(self.capacitor_offset) + d(self.capacitor_response) * d(q)).to_f64();
        let dq = (d(q) - d(self.history_q)) / d(self.step.effective_timestep_s);
        let dx = (d(x) - d(self.history_x)) / d(self.step.effective_timestep_s);
        let inertial_pressure = d(l.inertance) * dq;
        let proximal_pressure = d(l.proximal_resistance_s_inv) * d(q);
        let capacitor_flow = d(l.compliance_s2) * dx;
        let distal_flow = d(x) / d(l.distal_resistance_s_inv);
        let ex = d(0.5) * d(l.compliance_s2) * d(x) * d(x);
        let eq = d(0.5) * d(l.inertance) * d(q) * d(q);
        let previous_energy = d(0.5)
            * (d(l.compliance_s2) * d(self.previous_x) * d(self.previous_x)
                + d(l.inertance) * d(self.previous_q) * d(self.previous_q));
        let energy_change = (ex + eq - previous_energy) / d(self.step.timestep_s);
        let numerical_defect = d(x) * capacitor_flow + d(q) * inertial_pressure - energy_change;
        let proximal_loss = proximal_pressure * d(q);
        let distal_loss = distal_flow * d(x);
        let reservoir_power = d(l.reservoir_pressure_m2_s2) * d(q);
        let port_power = d(pressure) * d(q);
        let balance = StokesWindkesselBalance {
            compliance_pressure_m2_s2: x,
            compliance_energy_m4_s2: ex.to_f64(),
            inertance_energy_m4_s2: eq.to_f64(),
            proximal_dissipation_m4_s3: proximal_loss.to_f64(),
            distal_dissipation_m4_s3: distal_loss.to_f64(),
            reservoir_power_m4_s3: reservoir_power.to_f64(),
            port_power_m4_s3: port_power.to_f64(),
            stored_energy_change_rate_m4_s3: energy_change.to_f64(),
            bdf_energy_defect_m4_s3: numerical_defect.to_f64(),
            power_balance_residual_m4_s3: (port_power
                - reservoir_power
                - proximal_loss
                - distal_loss
                - energy_change
                - numerical_defect)
                .to_f64(),
            pressure_law_residual_m2_s2: (d(pressure)
                - d(l.reservoir_pressure_m2_s2)
                - proximal_pressure
                - inertial_pressure
                - d(x))
            .to_f64(),
            compliance_flow_residual_m2_s: (capacitor_flow - d(q) + distal_flow).to_f64(),
        };
        let pressure_scale = pressure
            .abs()
            .max(l.reservoir_pressure_m2_s2.abs())
            .max(proximal_pressure.to_f64().abs())
            .max(inertial_pressure.to_f64().abs())
            .max(x.abs());
        let flow_scale = q
            .abs()
            .max(capacitor_flow.to_f64().abs())
            .max(distal_flow.to_f64().abs());
        if ![
            q,
            pressure,
            x,
            pressure_scale,
            flow_scale,
            balance.compliance_energy_m4_s2,
            balance.inertance_energy_m4_s2,
            balance.proximal_dissipation_m4_s3,
            balance.distal_dissipation_m4_s3,
            balance.reservoir_power_m4_s3,
            balance.port_power_m4_s3,
            balance.stored_energy_change_rate_m4_s3,
            balance.bdf_energy_defect_m4_s3,
            balance.power_balance_residual_m4_s3,
            balance.pressure_law_residual_m2_s2,
            balance.compliance_flow_residual_m2_s,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err("unrepresentable hydraulic storage, power or original dynamic law");
        }
        Ok((balance, pressure_scale, flow_scale))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load() -> StokesWindkesselPort {
        StokesWindkesselPort {
            reservoir_pressure_m2_s2: 1.7,
            proximal_resistance_s_inv: 0.3,
            distal_resistance_s_inv: 2.1,
            compliance_s2: 0.6,
            inertance: 0.4,
        }
    }

    /// G1: independently solve the original two-state equations for a prescribed
    /// port pressure. This checks storage feedback, reversal and the full BDF
    /// derivative, rather than deriving the expected value from the reduction.
    #[test]
    fn endpoint_matches_independent_two_by_two_system_and_power() {
        let l = load();
        for old_dt in [None, Some(0.1), Some(0.07), Some(0.2)] {
            let dt = 0.1;
            let (qp, xp, qo, xo) = (0.2, -0.3, -0.1, 0.4);
            let (a0, a1, a2) = old_dt.map_or((1.0, -1.0, 0.0), |h| {
                let r = dt / h;
                ((1.0 + 2.0 * r) / (1.0 + r), -(1.0 + r), r * r / (1.0 + r))
            });
            let step = StokesBdfStep::new(dt, old_dt).unwrap();
            let reduced =
                ReducedWindkesselPort::new(l, step, qp, xp, old_dt.map(|_| (qo, xo))).unwrap();
            for p in [-3.0, 0.0, 3.0] {
                let a = l.proximal_resistance_s_inv + l.inertance * a0 / dt;
                let b = l.compliance_s2 * a0 / dt + 1.0 / l.distal_resistance_s_inv;
                let rhs_p = p - l.reservoir_pressure_m2_s2 - l.inertance * (a1 * qp + a2 * qo) / dt;
                let rhs_x = -l.compliance_s2 * (a1 * xp + a2 * xo) / dt;
                let expected_q = (rhs_p * b - rhs_x) / (a * b + 1.0);
                let expected_x = (rhs_p + a * rhs_x) / (a * b + 1.0);
                let q = (p - reduced.pressure_offset) / reduced.resistance;
                let (balance, _, _) = reduced.report(q, p).unwrap();
                assert!((q - expected_q).abs() < 2e-15);
                assert!((balance.compliance_pressure_m2_s2 - expected_x).abs() < 2e-15);
                assert!(balance.pressure_law_residual_m2_s2.abs() < 5e-15);
                assert!(balance.compliance_flow_residual_m2_s.abs() < 5e-15);
                assert!(balance.power_balance_residual_m4_s3.abs() < 5e-15);
                assert!(balance.proximal_dissipation_m4_s3 >= 0.0);
                assert!(balance.distal_dissipation_m4_s3 >= 0.0);
                if old_dt.is_none() {
                    let numerical_loss = (l.compliance_s2 * (expected_x - xp).powi(2)
                        + l.inertance * (expected_q - qp).powi(2))
                        / (2.0 * dt);
                    assert!((balance.bdf_energy_defect_m4_s3 - numerical_loss).abs() < 2e-15);
                }
            }
        }
    }

    /// G1: a closed port relaxes C*x'=-x/Rd; a precharged compliance can return
    /// energy. Physical resistor loss must not include the numerical BDF term.
    #[test]
    fn zero_inertance_relaxation_and_steady_series_resistance() {
        let mut l = load();
        l.inertance = 0.0;
        l.proximal_resistance_s_inv = 0.0;
        let dt = 0.2;
        let step = StokesBdfStep::new(dt, None).unwrap();
        let reduced = ReducedWindkesselPort::new(l, step, 0.0, 2.0, None).unwrap();
        let (b, _, _) = reduced.report(0.0, reduced.pressure_offset).unwrap();
        let expected = 2.0 / (1.0 + dt / (l.compliance_s2 * l.distal_resistance_s_inv));
        assert!((b.compliance_pressure_m2_s2 - expected).abs() < 1e-15);
        assert!(b.stored_energy_change_rate_m4_s3 < 0.0);
        assert!(
            (b.stored_energy_change_rate_m4_s3
                + b.distal_dissipation_m4_s3
                + b.bdf_energy_defect_m4_s3)
                .abs()
                < 2e-15
        );
        for q in [-0.7, 0.0, 0.7] {
            let x = l.distal_resistance_s_inv * q;
            let p = l.reservoir_pressure_m2_s2
                + (l.proximal_resistance_s_inv + l.distal_resistance_s_inv) * q;
            let reduced = ReducedWindkesselPort::new(l, step, q, x, None).unwrap();
            let (b, _, _) = reduced.report(q, p).unwrap();
            assert!(b.pressure_law_residual_m2_s2.abs() < 1e-15);
            assert!(b.compliance_flow_residual_m2_s.abs() < 1e-15);
            assert!(b.stored_energy_change_rate_m4_s3.abs() < 2e-15);
        }
    }

    /// G1: BDF1 startup plus the actual alternating-step BDF2 storage update
    /// converges quadratically to the closed-port exponential relaxation.
    #[test]
    fn unequal_step_relaxation_is_second_order_and_defect_is_signed() {
        let mut l = load();
        l.inertance = 0.0;
        let error = |count: usize| {
            let base_dt = 1.0 / count as f64;
            let mut recent = 2.0;
            let mut older = 2.0;
            let mut old_dt = None;
            for i in 0..count {
                let dt = base_dt * if i % 2 == 0 { 0.75 } else { 1.25 };
                let step = StokesBdfStep::new(dt, old_dt).unwrap();
                let reduced =
                    ReducedWindkesselPort::new(l, step, 0.0, recent, old_dt.map(|_| (0.0, older)))
                        .unwrap();
                let (balance, _, _) = reduced.report(0.0, reduced.pressure_offset).unwrap();
                older = recent;
                recent = balance.compliance_pressure_m2_s2;
                old_dt = Some(dt);
            }
            (recent - 2.0 * (-1.0 / (l.compliance_s2 * l.distal_resistance_s_inv)).exp()).abs()
        };
        let errors = [error(20), error(40), error(80)];
        assert!(errors[0] / errors[1] > 3.7, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.7, "{errors:?}");
        let step = StokesBdfStep::new(0.1, Some(0.1)).unwrap();
        let reduced = ReducedWindkesselPort::new(l, step, 0.0, 0.0, Some((0.0, 1.0))).unwrap();
        let (balance, _, _) = reduced.report(0.0, reduced.pressure_offset).unwrap();
        assert!(balance.bdf_energy_defect_m4_s3 < 0.0);
        assert!(balance.distal_dissipation_m4_s3 > 0.0);
        assert!(balance.power_balance_residual_m4_s3.abs() < 1e-15);
    }

    #[test]
    fn nonpassive_nonfinite_and_unrepresentable_loads_are_refused() {
        let step = StokesBdfStep::new(0.1, None).unwrap();
        let l = load();
        for invalid in [
            StokesWindkesselPort {
                compliance_s2: 0.0,
                ..l
            },
            StokesWindkesselPort {
                distal_resistance_s_inv: 0.0,
                ..l
            },
            StokesWindkesselPort {
                proximal_resistance_s_inv: -1.0,
                ..l
            },
            StokesWindkesselPort {
                inertance: -1.0,
                ..l
            },
            StokesWindkesselPort {
                reservoir_pressure_m2_s2: f64::NAN,
                ..l
            },
            StokesWindkesselPort {
                compliance_s2: f64::MAX,
                ..l
            },
        ] {
            assert!(ReducedWindkesselPort::new(invalid, step, 0.0, 0.0, None).is_err());
        }
        assert!(ReducedWindkesselPort::new(l, step, 0.0, f64::INFINITY, None).is_err());
        assert!(ReducedWindkesselPort::new(l, step, 0.0, 0.0, Some((f64::NAN, 0.0))).is_err());
    }
}
