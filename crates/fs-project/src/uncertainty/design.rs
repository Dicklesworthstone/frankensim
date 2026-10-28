//! Deterministic power and fan sizing of the existing native cooling model.
//! This reuses coordinate/asset binding, never a probability law or sampler.

use super::*;

mod cases;
pub use cases::DesignOperatingCase;

/// The physical coordinate and direction of a native cooling design search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DesignControl {
    /// Maximize delivered regional watts under a nondecreasing temperature law.
    DeliveredPower {
        /// Maximum passing/failing delivered-power bracket width, W.
        tolerance_w: f64,
    },
    /// Maximize one region's original delivered-watt schedule amplitude.
    /// Event times and other scheduled regions remain fixed.
    PowerScheduleScale {
        /// Maximum passing/failing dimensionless amplitude bracket width.
        tolerance: f64,
    },
    /// Minimize one bank's absolute source-curve speed ratio under an explicitly
    /// assumed nonincreasing temperature law. Other banks remain unchanged.
    FanSpeedRatio {
        /// Maximum passing/failing absolute speed-ratio bracket width.
        tolerance: f64,
    },
}

impl DesignControl {
    /// Bracket tolerance in the selected physical coordinate's coherent units.
    #[must_use]
    pub const fn tolerance(self) -> f64 {
        match self {
            Self::DeliveredPower { tolerance_w } => tolerance_w,
            Self::FanSpeedRatio { tolerance } | Self::PowerScheduleScale { tolerance } => tolerance,
        }
    }
}

/// Bracket tolerances for one explicitly monotone physical coordinate.
/// Feasibility always comes from an evaluated native temperature, not a model
/// prediction. These tolerances provide no continuous-time or continuum bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignPolicy {
    /// Typed design coordinate, optimization direction and coordinate tolerance.
    pub control: DesignControl,
    /// Maximum temperature slack at a bracketed passing endpoint, K.
    pub temperature_tolerance_k: f64,
}

pub(super) fn parse(root: &Node) -> Result<UncertaintyStudy> {
    let nodes = list(root)?;
    let f = fields_with_optional(&nodes[1..], &["version", "project", "evaluations", "seed", "wall-time",
        "qoi", "geometry", "materials", "interfaces", "control", "temperature-tolerance",
        "assumption"], &["operating-cases"])?;
    if integer(f["version"])? != 1 {
        return Err(error("native cooling design requires version 1"));
    }
    let samples = usize::try_from(integer(f["evaluations"])? )
        .map_err(|_| error("design evaluation count overflow"))?;
    if !(2..=MAX_SAMPLES).contains(&samples) {
        return Err(error("native design evaluations must be in 2..=256"));
    }
    let wall_seconds = quantity(f["wall-time"], crate::spec::dims::TIME)?;
    if !(wall_seconds > 0.0 && wall_seconds <= 86_400.0) {
        return Err(error("design wall-time must be in (0, 86400] seconds"));
    }
    let qoi = text(f["qoi"])?;
    if qoi != "temperature-max" {
        return Err(error("native cooling design requires temperature-max"));
    }
    let control = list(f["control"])?;
    let kind = control.first().ok_or_else(|| error("empty design control"))?;
    let target = match &kind.kind {
        NodeKind::Symbol(name) if name == "regional-power" => Target::Power,
        NodeKind::Symbol(name) if name == "fan-speed-ratio" => Target::FanSpeedRatio,
        NodeKind::Symbol(name) if name == "power-schedule-scale" => Target::PowerScheduleScale,
        _ => return Err(error("design control must be regional-power, power-schedule-scale or fan-speed-ratio")),
    };
    symbol(f["assumption"], match target {
        Target::Power => "nondecreasing-temperature-with-power",
        Target::PowerScheduleScale => "nondecreasing-temperature-with-schedule-scale",
        _ => "nonincreasing-temperature-with-fan-speed",
    })?;
    let c = fields(&control[1..], &["entity", "minimum", "maximum", "tolerance"])?;
    let low = quantity(c["minimum"], target.dims())?;
    let high = quantity(c["maximum"], target.dims())?;
    let tolerance = quantity(c["tolerance"], target.dims())?;
    let temperature_tolerance_k = quantity(f["temperature-tolerance"], crate::spec::dims::TEMPERATURE)?;
    if !(low >= 0.0 && high > low && tolerance > 0.0
        && temperature_tolerance_k > 0.0 && (high - low).is_finite()) {
        return Err(error("design requires 0 <= minimum < maximum and positive coordinate/temperature tolerances"));
    }
    if target == Target::FanSpeedRatio && low == 0.0 {
        return Err(error("fan design requires a strictly positive minimum speed ratio"));
    }
    let control = match target {
        Target::Power => DesignControl::DeliveredPower { tolerance_w: tolerance },
        Target::PowerScheduleScale => DesignControl::PowerScheduleScale { tolerance },
        _ => DesignControl::FanSpeedRatio { tolerance },
    };
    let geometry_nodes = list(f["geometry"])?;
    if geometry_nodes.is_empty() || geometry_nodes.len() > 32 {
        return Err(error("declare 1..=32 mesh sources"));
    }
    let mut geometry = Vec::new();
    let mut roles = BTreeSet::new();
    for node in geometry_nodes {
        let n = list(node)?;
        symbol(n.first().ok_or_else(|| error("empty mesh source"))?, "mesh")?;
        let m = fields(&n[1..], &["role", "path", "unit", "max-hole-edges"])?;
        let role = text(m["role"])?;
        if !roles.insert(role.clone()) {
            return Err(error("duplicate geometry role"));
        }
        let max_hole_edges = usize::try_from(integer(m["max-hole-edges"])? )
            .map_err(|_| error("repair count overflow"))?;
        if max_hole_edges > 10_000 {
            return Err(error("repair edge budget exceeds 10000"));
        }
        geometry.push(MeshSource { role, path: text(m["path"])?, unit: text(m["unit"])?, max_hole_edges });
    }
    let coordinate = UniformParameter {
        name: match target {
            Target::Power => "delivered-power",
            Target::PowerScheduleScale => "schedule-amplitude",
            _ => "fan-speed",
        }.into(),
        target, entity: text(c["entity"])?, low, high,
    };
    let (mut operating_parameters, operating_cases) = f.get("operating-cases")
        .map(|node| cases::parse(node, &coordinate)).transpose()?.unwrap_or_default();
    let mut parameters = vec![coordinate];
    parameters.append(&mut operating_parameters);
    Ok(UncertaintyStudy {
        canonical: fs_ir::sexpr::print(root).map_err(|e| error(e.to_string()))?,
        project: text(f["project"])?, samples, seed: integer(f["seed"])?, wall_seconds, qoi,
        geometry, materials: paths(f["materials"])?, interfaces: paths(f["interfaces"])? ,
        // The existing physical-coordinate binder owns every mutation. Search
        // bounds and finite cases never acquire a probability interpretation.
        parameters, operating_cases,
        latent_correlation: None, compliance: None, qmc: None, lhs: None,
        sobol_sensitivity: false, mean_control: None,
        design: Some(DesignPolicy { control, temperature_tolerance_k }),
    })
}

/// A schedule must never silently override the coordinate being optimized.
pub(super) fn admit(study: &UncertaintyStudy, base: &ProjectSpec) -> Result<()> {
    let Some(policy) = study.design else { return Ok(()); };
    let setup = base.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| error("native cooling design requires declared conduction"))?;
    if matches!(policy.control, DesignControl::FanSpeedRatio { .. }) {
        // A fan that never participates in the physical thermal model cannot
        // be sized by a declared convection coefficient unrelated to its flow.
        // Exact bank/domain admission remains owned by the existing binder.
        if !setup.boundaries.iter().any(|boundary|
            matches!(boundary.condition, ThermalBoundaryCondition::AirflowConvection { .. })) {
            return Err(error("fan-speed design requires an airflow-convection boundary so the native thermal solve recomputes fan-driven cooling"));
        }
    }
    if matches!(policy.control, DesignControl::PowerScheduleScale { .. }) {
        let region = &study.parameters[0].entity;
        if !setup.transient.as_ref().and_then(|time|
            time.power_schedules.iter().find(|schedule| &schedule.region == region))
            .is_some_and(|schedule| schedule.steps.iter().any(|step| step.watts.value > 0.0))
        {
            return Err(error("schedule-amplitude design requires an original schedule with nonzero delivered heating"));
        }
    }
    // Operating-case power is delivered watts too. A time schedule must not
    // silently replace either the search control or an authored operating case.
    for parameter in study.parameters.iter().filter(|p| p.target == Target::Power) {
        if !setup.regions.iter().any(|region| region.region == parameter.entity) {
            return Err(error("regional-power design must name a volumetric conduction region"));
        }
        let power: Vec<_> = base.power.iter().flatten()
            .filter(|row| row.region == parameter.entity).collect();
        if power.len() != 1 || !(power[0].duty.is_finite() && power[0].duty > 0.0) {
            return Err(error("delivered-power design needs one matching power row with positive unchanged duty"));
        }
        if setup.transient.as_ref().is_some_and(|time|
            time.power_schedules.iter().any(|schedule| schedule.region == parameter.entity)) {
            return Err(error("a transient power schedule overrides this static regional-power control; select an unscheduled region or explicitly use power-schedule-scale"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SOURCE: &str = r#"(fsim-native-design-study
        :version 1 :project "cooling-reference.fsim" :evaluations 24 :seed 29
        :wall-time 120s :qoi "temperature-max"
        :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
        :materials ("aa6061.fsmcdpk") :interfaces ()
        :control (regional-power :entity "air" :minimum 0W :maximum 20W :tolerance 0.01W)
        :temperature-tolerance 0.01K :assumption nondecreasing-temperature-with-power)"#;

    fn base() -> ProjectSpec {
        crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-reference.fsim"))).unwrap().decoded.spec
    }

    #[test]
    fn g0_design_is_canonical_without_a_sampling_or_probability_policy() {
        let study = UncertaintyStudy::parse(SOURCE).unwrap();
        assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
        assert_eq!(study.samples(), 24);
        assert_eq!(study.design().unwrap().control, DesignControl::DeliveredPower { tolerance_w: 0.01 });
        assert!(study.compliance().is_none());
        assert!(study.mean_control().is_none());
        for (from, to) in [
            (":evaluations 24", ":evaluations 1"),
            (":maximum 20W", ":maximum 0W"),
            (":minimum 0W", ":minimum -1W"),
            (":tolerance 0.01W", ":tolerance 0.01K"),
            (":temperature-tolerance 0.01K", ":temperature-tolerance 0K"),
            ("nondecreasing-temperature-with-power", "assumed-safe"),
            (":version 1", ":version 2"),
            (":qoi \"temperature-max\"", ":qoi \"temperature-max\" :method monte-carlo"),
        ] {
            assert!(UncertaintyStudy::parse(&SOURCE.replace(from, to)).is_err(), "{to}");
        }
    }

    #[test]
    fn g0_scheduled_replacement_cannot_turn_the_design_coordinate_into_a_no_op() {
        use crate::{ConductionTransient, TransientInitialTemperature, TransientRegionCapacity,
            TransientRegionPower, TransientPowerStep, spec::dims};
        let mut base = base();
        let setup = base.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.transient = Some(ConductionTransient {
            initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(293.15, dims::TEMPERATURE)),
            horizon: QtyAny::new(1.0, dims::TIME),
            max_step: QtyAny::new(1.0, dims::TIME), max_steps: 3,
            energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
            capacities: setup.regions.iter().map(|row| TransientRegionCapacity {
                region: row.region.clone(),
                volumetric_heat_capacity: QtyAny::new(100.0, dims::VOLUMETRIC_HEAT_CAPACITY),
                source: "synthetic storage".into(),
            }).collect(),
            enthalpy: None,
            power_schedules: vec![TransientRegionPower { region: "air".into(),
                source: "authored delivered-power replacement".into(),
                steps: vec![TransientPowerStep { until: QtyAny::new(1.0, dims::TIME),
                    watts: QtyAny::new(5.0, dims::POWER) }] }],
        });
        let error = UncertaintyStudy::parse(SOURCE).unwrap().bind(&base).unwrap_err();
        assert!(error.detail.contains("overrides"), "{error:?}");
    }

    #[test]
    fn g1_delivered_power_preserves_duty_and_all_other_native_fields() {
        let mut base = base();
        base.power.as_mut().unwrap()[0].duty = 0.25;
        let bound = UncertaintyStudy::parse(SOURCE).unwrap().bind(&base).unwrap();
        let candidate = bound.sample_project(&[7.0]).unwrap();
        let power = &candidate.power.as_ref().unwrap()[0];
        assert_eq!(power.watts.value, 28.0);
        assert_eq!(power.duty, 0.25);
        let mut expected = base.clone();
        expected.power.as_mut().unwrap()[0].watts.value = 28.0;
        assert_eq!(candidate, expected);
        assert_eq!(bound.base(), &base);
        base.power.as_mut().unwrap()[0].duty = 0.0;
        assert!(UncertaintyStudy::parse(SOURCE).unwrap().bind(&base).is_err());
        assert!(UncertaintyStudy::parse(&SOURCE.replace(":entity \"air\"", ":entity \"missing\""))
            .unwrap().bind(bound.base()).is_err());
    }

    fn fan_source() -> String {
        SOURCE.replace("regional-power :entity \"air\" :minimum 0W :maximum 20W :tolerance 0.01W",
            "fan-speed-ratio :entity \"fixture-bank\" :minimum 0.6 :maximum 1.4 :tolerance 0.001")
            .replace("nondecreasing-temperature-with-power", "nonincreasing-temperature-with-fan-speed")
    }

    fn fan_base() -> ProjectSpec {
        let mut project = base();
        let cooling = project.cooling.as_mut().unwrap();
        cooling.fan_system.as_mut().unwrap().banks[0].speed_ratio = 0.8;
        cooling.conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::AirflowConvection {
                branch: "air".into(), order: 0,
                inlet_temperature: QtyAny::new(295.0, crate::spec::dims::TEMPERATURE),
                hydraulic_diameter: QtyAny::new(0.02, crate::spec::dims::LENGTH),
                flow_area: QtyAny::new(0.004, crate::spec::dims::AREA),
                channel_length: QtyAny::new(0.3, crate::spec::dims::LENGTH),
                correlation: "convection.gnielinski".into(),
            };
        project
    }

    #[test]
    fn g0_fan_design_requires_dimensionless_speed_and_its_own_monotonicity_assumption() {
        let source = fan_source();
        let study = UncertaintyStudy::parse(&source).unwrap();
        assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
        assert_eq!(study.parameters()[0].target, Target::FanSpeedRatio);
        assert_eq!(study.design().unwrap().control, DesignControl::FanSpeedRatio { tolerance: 0.001 });
        for (from, to) in [
            (":minimum 0.6", ":minimum 0"),
            (":minimum 0.6", ":minimum 0.6W"),
            (":maximum 1.4", ":maximum 0.6"),
            (":tolerance 0.001", ":tolerance 0.001W"),
            (":tolerance 0.001", ":tolerance 0"),
            ("nonincreasing-temperature-with-fan-speed", "nondecreasing-temperature-with-power"),
            ("fan-speed-ratio", "fan-electrical-power"),
        ] {
            assert!(UncertaintyStudy::parse(&source.replace(from, to)).is_err(), "{to}");
        }
    }

    #[test]
    fn g1_fan_design_changes_only_the_absolute_named_bank_ratio() {
        let mut base = fan_base();
        base.power.as_mut().unwrap()[0].duty = 0.25;
        let system = base.cooling.as_mut().unwrap().fan_system.as_mut().unwrap();
        let mut other = system.banks[0].clone();
        other.bank_id = "fixed-bank".into();
        other.speed_ratio = 1.2;
        system.banks.push(other);
        system.topology = crate::fansystem::FanSystemTopology::Series(
            vec!["fixture-bank".into(), "fixed-bank".into()]);
        let bound = UncertaintyStudy::parse(&fan_source()).unwrap().bind(&base).unwrap();
        for ratio in [0.6, 1.0, 1.4, 0.6] {
            let candidate = bound.sample_project(&[ratio]).unwrap();
            let mut expected = base.clone();
            expected.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio = ratio;
            assert_eq!(candidate, expected, "power, duty, curves, other banks and geometry stay unchanged");
        }
        assert_eq!(bound.base(), &base);
    }

    #[test]
    fn g0_fan_design_cannot_expand_domains_or_size_disconnected_cooling() {
        let source = fan_source();
        assert!(UncertaintyStudy::parse(&source).unwrap().bind(&base()).is_err(),
            "a scalar convection coefficient does not connect a fan to thermal physics");
        for (from, to) in [
            (":minimum 0.6", ":minimum 0.49"),
            (":maximum 1.4", ":maximum 2.01"),
            (":entity \"fixture-bank\"", ":entity \"missing\""),
        ] {
            assert!(UncertaintyStudy::parse(&source.replace(from, to)).unwrap().bind(&fan_base()).is_err(), "{to}");
        }
        let mut base = fan_base();
        base.cooling.as_mut().unwrap().fan_system = None;
        assert!(UncertaintyStudy::parse(&source).unwrap().bind(&base).is_err());
    }
}
