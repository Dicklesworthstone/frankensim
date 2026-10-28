//! Coordinate semantics for one shared search, ledger and export pipeline.
//! Keeping units and identity together prevents amplitude from being labeled W.
use super::{DesignControl, DesignPolicy, FAN_SCOPE, SCOPE};

pub(super) struct Presentation {
    pub coordinate: &'static str,
    pub objective: &'static str,
    pub entity: &'static str,
    pub tolerance: &'static str,
    pub bracket_width: &'static str,
    pub schema: &'static str,
    pub case_schema: &'static str,
    pub label: &'static str,
    pub unit: &'static str,
    pub title: &'static str,
    pub column: &'static str,
    pub scope: &'static str,
    pub capability: &'static str,
    pub baseline_refusal: &'static str,
    pub resolution_refusal: &'static str,
}

pub(super) const fn presentation(policy: DesignPolicy) -> Presentation {
    match policy.control {
        DesignControl::DeliveredPower { .. } => Presentation {
            coordinate: "power_w", objective: "maximize-delivered-regional-power",
            entity: "power_region", tolerance: "power_tolerance_w", bracket_width: "bracket_width_w",
            schema: "frankensim.cli.native-power-design-result.v1",
            case_schema: "frankensim.cli.native-power-case-design-result.v1",
            label: "delivered power", unit: " W", title: "Native power design", column: "Delivered W",
            scope: SCOPE,
            capability: "{\"ops\":[\"native-delivered-power-design\"],\"physical_ops\":\"inherited-from-retained-base-project\"}",
            baseline_refusal: "minimum-power-infeasible", resolution_refusal: "power-resolution-exhausted",
        },
        DesignControl::FanSpeedRatio { .. } => Presentation {
            coordinate: "fan_speed_ratio", objective: "minimize-fan-speed-ratio",
            entity: "fan_bank", tolerance: "fan_speed_ratio_tolerance", bracket_width: "bracket_width_ratio",
            schema: "frankensim.cli.native-fan-design-result.v1",
            case_schema: "frankensim.cli.native-fan-case-design-result.v1",
            label: "fan speed ratio", unit: " (absolute ratio)", title: "Native fan design", column: "Fan speed ratio",
            scope: FAN_SCOPE,
            capability: "{\"ops\":[\"native-fan-speed-design\"],\"physical_ops\":\"inherited-from-retained-base-project\"}",
            baseline_refusal: "maximum-fan-speed-infeasible", resolution_refusal: "fan-speed-resolution-exhausted",
        },
        DesignControl::PowerScheduleScale { .. } => Presentation {
            coordinate: "power_schedule_scale", objective: "maximize-regional-power-schedule-scale",
            entity: "schedule_region", tolerance: "schedule_scale_tolerance", bracket_width: "bracket_width_scale",
            schema: "frankensim.cli.native-schedule-design-result.v1",
            case_schema: "frankensim.cli.native-schedule-case-design-result.v1",
            label: "schedule amplitude", unit: " (original schedule multiplier)",
            title: "Native pulse amplitude design", column: "Schedule amplitude",
            scope: "Maximize one original regional delivered-watt schedule's dimensionless amplitude under the explicit nondecreasing-temperature-with-schedule-scale assumption. Each complete native trajectory starts from the original initial state; switching times, off intervals, static duty and unrelated schedules stay unchanged. Feasibility uses the fine-grid sampled peak including the initial state against the unchanged limit minus margin, never only the cooled final state. The incumbent is the highest-amplitude observed feasible project. No continuous-time or continuum temperature bound, global monotonicity proof, probability guarantee, physical validation or compliance signoff is supplied.",
            capability: "{\"ops\":[\"native-power-schedule-design\"],\"physical_ops\":\"inherited-from-retained-base-project\"}",
            baseline_refusal: "minimum-schedule-amplitude-infeasible",
            resolution_refusal: "schedule-amplitude-resolution-exhausted",
        },
    }
}
