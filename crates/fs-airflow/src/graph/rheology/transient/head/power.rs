//! Original pump/drive power values, independent of a thermal destination owner.
//!
//! The equipment evaluator constructs these same values; its historical public
//! paths re-export them. Native head evaluations retain that exact output and
//! check its balance against their actual hydraulic power. No source query,
//! efficiency estimate, shaft dynamics or thermal integration belongs here.

/// Authored fractions of the two distinct losses delivered to the finite
/// coolant inventory. The complements leave to unmodelled surroundings.
/// Neither fraction includes hydraulic work, whose eventual dissipation belongs
/// to the original flow network. Zero is external cooling; one is complete
/// coupling of that loss to explicitly identified thermal cells.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PumpLossRouting {
    /// Fraction of shaft input minus hydraulic output heating the inventory.
    pub pump_to_reservoir: f64,
    /// Fraction of electrical input minus shaft input heating the inventory.
    pub drive_to_reservoir: f64,
}

/// Actual equipment consumption and separated losses at one admitted state.
/// Electrical input = hydraulic output + reservoir heat + surroundings heat.
/// The existing network owns hydraulic output and its eventual dissipation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PumpPowerState {
    /// Mechanical input delivered to the pump shaft.
    pub shaft_input_w: f64,
    /// Total electrical input consumed by the drive.
    pub electrical_input_w: f64,
    /// Shaft input minus hydraulic output, never negative.
    pub pump_loss_w: f64,
    /// Electrical input minus shaft input, never negative.
    pub drive_loss_w: f64,
    /// Explicitly routed pump and drive heat entering the finite inventory.
    pub reservoir_heat_w: f64,
    /// Remaining pump and drive heat delivered to unmodelled surroundings.
    pub surroundings_heat_w: f64,
}

/// Exact original evaluator output together with its original loss partition.
/// The accepted head also retains the source object and original flow sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowHeadPowerReceipt {
    pub state: PumpPowerState,
    pub routing: PumpLossRouting,
}

impl FlowHeadPowerReceipt {
    /// Validate the original nonnegative pump/drive model without changing any
    /// value. Signed recovering head sources remain usable without this optional
    /// pump-power model. Zero hydraulic flow may have nonzero input/loss heat.
    pub fn validate(self, hydraulic_power_w: f64) -> Result<(), &'static str> {
        let s = self.state;
        let r = self.routing;
        if !hydraulic_power_w.is_finite()
            || hydraulic_power_w < 0.0
            || [
                s.shaft_input_w,
                s.electrical_input_w,
                s.pump_loss_w,
                s.drive_loss_w,
                s.reservoir_heat_w,
                s.surroundings_heat_w,
            ]
            .iter()
            .any(|p| !p.is_finite() || *p < 0.0)
            || [r.pump_to_reservoir, r.drive_to_reservoir]
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || s.shaft_input_w < hydraulic_power_w
            || s.electrical_input_w < s.shaft_input_w
        {
            return Err("original pump input, loss and routing support required");
        }
        let pump = s.shaft_input_w - hydraulic_power_w;
        let drive = s.electrical_input_w - s.shaft_input_w;
        let reservoir = r.pump_to_reservoir * s.pump_loss_w + r.drive_to_reservoir * s.drive_loss_w;
        let surroundings = (1.0 - r.pump_to_reservoir) * s.pump_loss_w
            + (1.0 - r.drive_to_reservoir) * s.drive_loss_w;
        // Each loss/route is its own physical channel. A large shaft input or
        // unrelated loss must not hide invented small heat in another channel.
        for (actual, expected) in [
            (s.pump_loss_w, pump),
            (s.drive_loss_w, drive),
            (s.reservoir_heat_w, reservoir),
            (s.surroundings_heat_w, surroundings),
        ] {
            if !expected.is_finite()
                || (actual - expected).abs()
                    > 64.0 * f64::EPSILON * actual.abs().max(expected.abs())
            {
                return Err(
                    "original pump power and loss partition disagree with actual hydraulic output",
                );
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_heat_cannot_hide_behind_large_power_or_another_loss() {
        let original = FlowHeadPowerReceipt {
            state: PumpPowerState {
                shaft_input_w: 1e16,
                electrical_input_w: 1e16,
                pump_loss_w: 0.0,
                drive_loss_w: 0.0,
                reservoir_heat_w: 0.0,
                surroundings_heat_w: 0.0,
            },
            routing: PumpLossRouting {
                pump_to_reservoir: 1.0,
                drive_to_reservoir: 0.0,
            },
        };
        assert!(original.validate(1e16).is_ok());
        let mut invented = original;
        invented.state.pump_loss_w = 1.0;
        invented.state.reservoir_heat_w = 1.0;
        assert!(invented.validate(1e16).is_err());

        let mut misrouted = FlowHeadPowerReceipt {
            state: PumpPowerState {
                shaft_input_w: 2e16,
                electrical_input_w: 2e16,
                pump_loss_w: 1e16,
                drive_loss_w: 0.0,
                reservoir_heat_w: 0.0,
                surroundings_heat_w: 1e16,
            },
            routing: PumpLossRouting {
                pump_to_reservoir: 0.0,
                drive_to_reservoir: 1.0,
            },
        };
        assert!(misrouted.validate(1e16).is_ok());
        // Large pump waste all leaves outside; it cannot fund one watt falsely
        // credited to the independently routed, lossless drive.
        misrouted.state.reservoir_heat_w = 1.0;
        assert!(misrouted.validate(1e16).is_err());
    }
}
