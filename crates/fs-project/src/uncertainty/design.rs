//! Deterministic delivered-power sizing of the existing native cooling model.
//! This reuses coordinate/asset binding, never a probability law or sampler.

use super::*;

/// Bracket tolerances for one explicitly monotone delivered-power coordinate.
/// Feasibility always comes from an evaluated native temperature, not a model
/// prediction. These tolerances provide no continuous-time or continuum bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignPolicy {
    /// Maximum passing/failing delivered-power bracket width, W.
    pub power_tolerance_w: f64,
    /// Maximum temperature slack at a bracketed passing endpoint, K.
    pub temperature_tolerance_k: f64,
}

pub(super) fn parse(root: &Node) -> Result<UncertaintyStudy> {
    let nodes = list(root)?;
    let f = fields(&nodes[1..], &["version", "project", "evaluations", "seed", "wall-time",
        "qoi", "geometry", "materials", "interfaces", "control", "temperature-tolerance",
        "assumption"])?;
    if integer(f["version"])? != 1 {
        return Err(error("native delivered-power design requires version 1"));
    }
    symbol(f["assumption"], "nondecreasing-temperature-with-power")?;
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
        return Err(error("native delivered-power design requires temperature-max"));
    }
    let control = list(f["control"])?;
    symbol(control.first().ok_or_else(|| error("empty design control"))?, "regional-power")?;
    let c = fields(&control[1..], &["entity", "minimum", "maximum", "tolerance"])?;
    let low = quantity(c["minimum"], crate::spec::dims::POWER)?;
    let high = quantity(c["maximum"], crate::spec::dims::POWER)?;
    let power_tolerance_w = quantity(c["tolerance"], crate::spec::dims::POWER)?;
    let temperature_tolerance_k = quantity(f["temperature-tolerance"], crate::spec::dims::TEMPERATURE)?;
    if !(low >= 0.0 && high > low && power_tolerance_w > 0.0
        && temperature_tolerance_k > 0.0 && (high - low).is_finite()) {
        return Err(error("design requires 0 <= minimum < maximum and positive power/temperature tolerances"));
    }
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
    Ok(UncertaintyStudy {
        canonical: fs_ir::sexpr::print(root).map_err(|e| error(e.to_string()))?,
        project: text(f["project"])?, samples, seed: integer(f["seed"])?, wall_seconds, qoi,
        geometry, materials: paths(f["materials"])?, interfaces: paths(f["interfaces"])? ,
        // The existing physical-coordinate binder owns field mutation. On this
        // deterministic path these are search bounds, not a uniform input law.
        parameters: vec![UniformParameter { name: "delivered-power".into(),
            target: Target::Power, entity: text(c["entity"])?, low, high }],
        latent_correlation: None, compliance: None, qmc: None, lhs: None,
        sobol_sensitivity: false, mean_control: None,
        design: Some(DesignPolicy { power_tolerance_w, temperature_tolerance_k }),
    })
}

/// A schedule must never silently override the coordinate being optimized.
pub(super) fn admit(study: &UncertaintyStudy, base: &ProjectSpec) -> Result<()> {
    if study.design.is_none() { return Ok(()); }
    let control = &study.parameters[0];
    let setup = base.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| error("native power design requires declared conduction"))?;
    if !setup.regions.iter().any(|region| region.region == control.entity) {
        return Err(error("regional-power design must name a volumetric conduction region"));
    }
    let power: Vec<_> = base.power.iter().flatten().filter(|row| row.region == control.entity).collect();
    if power.len() != 1 || !(power[0].duty.is_finite() && power[0].duty > 0.0) {
        return Err(error("delivered-power design needs one matching power row with positive unchanged duty"));
    }
    if setup.transient.as_ref().is_some_and(|time|
        time.power_schedules.iter().any(|schedule| schedule.region == control.entity)) {
        return Err(error("a transient power schedule overrides this static regional-power control; select an unscheduled region"));
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
        assert_eq!(study.design().unwrap().power_tolerance_w, 0.01);
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
}
