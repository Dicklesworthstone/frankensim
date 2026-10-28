//! One amplitude scales the original delivered-watt history, not static power.
//! Switching times, duty, other regions and the physical initial state stay fixed.

use super::{ProjectSpec, Result, Target, UncertaintyStudy, UniformParameter, error};

pub(super) fn admit(study: &UncertaintyStudy, base: &ProjectSpec) -> Result<()> {
    for parameter in study.parameters.iter().filter(|p| p.target == Target::PowerScheduleScale) {
        let setup = base.cooling.as_ref().and_then(|c| c.conduction.as_ref())
            .ok_or_else(|| error("power-schedule-scale requires native conduction"))?;
        if !setup.regions.iter().any(|region| region.region == parameter.entity) {
            return Err(error("power-schedule-scale must name a volumetric conduction region"));
        }
        let time = setup.transient.as_ref()
            .ok_or_else(|| error("power-schedule-scale requires a declared transient"))?;
        if time.power_schedules.iter().filter(|s| s.region == parameter.entity).count() != 1 {
            return Err(error("power-schedule-scale requires exactly one original regional power schedule"));
        }
        if study.parameters.iter().any(|p| p.target == Target::Power && p.entity == parameter.entity) {
            return Err(error("a schedule amplitude and its overridden static power cannot be separate study coordinates"));
        }
    }
    Ok(())
}

pub(super) fn apply(project: &mut ProjectSpec, parameter: &UniformParameter, scale: f64) -> Result<()> {
    if !scale.is_finite() || scale < 0.0 {
        return Err(error("power-schedule-scale requires a finite nonnegative amplitude"));
    }
    let time = project.cooling.as_mut().and_then(|c| c.conduction.as_mut())
        .and_then(|c| c.transient.as_mut())
        .ok_or_else(|| error("power-schedule-scale requires a declared transient"))?;
    let mut matches = time.power_schedules.iter_mut().filter(|s| s.region == parameter.entity);
    let schedule = matches.next()
        .ok_or_else(|| error("power-schedule-scale has no matching regional schedule"))?;
    if matches.next().is_some() {
        return Err(error("power-schedule-scale has an ambiguous regional schedule"));
    }
    // Preflight all products before modifying even this disposable candidate.
    // A representable positive pulse must not silently become an off interval.
    for step in &schedule.steps {
        scaled_watts(step.watts.value, scale)?;
    }
    for step in &mut schedule.steps {
        step.watts.value *= scale;
    }
    Ok(())
}

fn scaled_watts(watts: f64, scale: f64) -> Result<f64> {
    let result = watts * scale;
    if !watts.is_finite() || watts < 0.0 || !result.is_finite()
        || (watts > 0.0 && scale > 0.0 && result == 0.0)
    {
        Err(error("scaled delivered power is invalid or not representable; no pulse is clipped or erased"))
    } else {
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::QtyAny;
    use crate::{ConductionTransient, TransientInitialTemperature, TransientRegionCapacity,
        TransientRegionPower, TransientPowerStep, spec::dims};

    const SOURCE: &str = r#"(fsim-uncertainty-study
        :version 1 :project "cooling-reference.fsim" :samples 4 :seed 29 :wall-time 120s
        :method monte-carlo :correlation independent :qoi "temperature-peak"
        :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
        :materials ("aa6061.fsmcdpk") :interfaces ()
        :parameters ((uniform :name "pulse" :target power-schedule-scale :entity "air" :low 0 :high 2)))"#;

    fn base() -> ProjectSpec {
        let mut base = crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-reference.fsim"))).unwrap().decoded.spec;
        base.power.as_mut().unwrap()[0].duty = 0.25;
        let setup = base.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        setup.transient = Some(ConductionTransient {
            initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(293.15, dims::TEMPERATURE)),
            horizon: QtyAny::new(2.0, dims::TIME), max_step: QtyAny::new(0.5, dims::TIME),
            max_steps: 24, energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
            capacities: setup.regions.iter().map(|row| TransientRegionCapacity {
                region: row.region.clone(), volumetric_heat_capacity: QtyAny::new(100.0, dims::VOLUMETRIC_HEAT_CAPACITY),
                source: "synthetic capacity".into(),
            }).collect(), enthalpy: None,
            power_schedules: vec![TransientRegionPower {
                region: "air".into(), source: "synthetic delivered pulse".into(),
                steps: [(0.25, 8.0), (0.75, 0.0), (2.0, 2.0)].into_iter()
                    .map(|(until, watts)| TransientPowerStep {
                        until: QtyAny::new(until, dims::TIME), watts: QtyAny::new(watts, dims::POWER),
                    }).collect(),
            }],
        });
        base
    }

    fn schedule(project: &ProjectSpec) -> &TransientRegionPower {
        &project.cooling.as_ref().unwrap().conduction.as_ref().unwrap()
            .transient.as_ref().unwrap().power_schedules[0]
    }

    #[test]
    fn g0_schedule_scale_is_explicit_dimensionless_nonnegative_and_canonical() {
        let parsed = UncertaintyStudy::parse(SOURCE).unwrap();
        assert_eq!(UncertaintyStudy::parse(parsed.canonical()).unwrap(), parsed);
        assert_eq!(parsed.parameters()[0].target.unit(), "1");
        for (from, to) in [(":low 0", ":low -1"), (":high 2", ":high 2W"),
            ("power-schedule-scale", "power-schedule"), (":high 2", ":high -1")] {
            assert!(UncertaintyStudy::parse(&SOURCE.replace(from, to)).is_err(), "{to}");
        }
    }

    #[test]
    fn g1_g3_scaling_preserves_events_off_periods_duty_and_original_history() {
        let original = base();
        let bound = UncertaintyStudy::parse(SOURCE).unwrap().bind(&original).unwrap();
        for scale in [2.0, 0.0, 0.5, 1.0, 2.0] {
            let candidate = bound.sample_project(&[scale]).unwrap();
            let mut expected = original.clone();
            for step in &mut expected.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .transient.as_mut().unwrap().power_schedules[0].steps {
                step.watts.value *= scale;
            }
            assert_eq!(candidate, expected);
            assert_eq!(bound.base(), &original);
            let mut start = 0.0;
            let energy: f64 = schedule(&candidate).steps.iter().map(|s| {
                let heat = (s.until.value - start) * s.watts.value;
                start = s.until.value;
                heat
            }).sum();
            assert_eq!(energy, 4.5 * scale, "delivered watts must not receive duty twice");
            assert_eq!(candidate.cooling.as_ref().unwrap().conduction.as_ref().unwrap().transient.as_ref().unwrap()
                .coarse_step_ends_s().unwrap(), original.cooling.as_ref().unwrap().conduction.as_ref().unwrap()
                .transient.as_ref().unwrap().coarse_step_ends_s().unwrap());
        }
    }

    #[test]
    fn g0_missing_ambiguous_and_overridden_coordinates_refuse_without_mutation() {
        let original = base();
        for mode in 0..3 {
            let mut invalid = original.clone();
            let setup = invalid.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
            match mode {
                0 => setup.transient = None,
                1 => setup.transient.as_mut().unwrap().power_schedules.clear(),
                _ => {
                    let time = setup.transient.as_mut().unwrap();
                    time.power_schedules.push(time.power_schedules[0].clone());
                }
            }
            let before = invalid.clone();
            assert!(UncertaintyStudy::parse(SOURCE).unwrap().bind(&invalid).is_err());
            assert_eq!(invalid, before);
        }
        let conflict = SOURCE.replace(":high 2)))", ":high 2) (uniform :name \"static\" :target power :entity \"air\" :low 1W :high 2W)))");
        assert!(UncertaintyStudy::parse(&conflict).unwrap().bind(&original).is_err());
        let wrong = SOURCE.replace(":entity \"air\"", ":entity \"missing\"");
        assert!(UncertaintyStudy::parse(&wrong).unwrap().bind(&original).is_err());
    }

    #[test]
    fn g0_unrepresentable_scaling_is_atomic_and_zero_is_not_a_source_domain_escape() {
        let parameter = UncertaintyStudy::parse(SOURCE).unwrap().parameters()[0].clone();
        for (watts, scale) in [(f64::MAX, 2.0), (f64::from_bits(1), 0.25)] {
            let mut project = base();
            project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap()
                .power_schedules[0].steps[2].watts.value = watts;
            let before = project.clone();
            assert!(apply(&mut project, &parameter, scale).is_err());
            assert_eq!(project, before);
        }
        assert_eq!(scaled_watts(f64::MAX, 0.0).unwrap(), 0.0);
        assert!(scaled_watts(f64::NAN, 0.0).is_err());
        let mut original = base();
        original.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap()
            .power_schedules[0].steps[0].watts.value = -1.0;
        assert!(UncertaintyStudy::parse(&SOURCE.replace(":high 2", ":high 0")).unwrap().bind(&original).is_err());
    }
    fn design_source() -> &'static str {
        r#"(fsim-native-design-study :version 1 :project "cooling-reference.fsim"
          :evaluations 32 :seed 29 :wall-time 120s :qoi "temperature-max"
          :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
          :materials ("aa6061.fsmcdpk") :interfaces ()
          :control (power-schedule-scale :entity "air" :minimum 0 :maximum 2 :tolerance 0.001)
          :temperature-tolerance 0.01K :assumption nondecreasing-temperature-with-schedule-scale)"#
    }

    #[test]
    fn g0_design_amplitude_has_explicit_units_assumption_and_nonzero_original_workload() {
        use crate::uncertainty::DesignControl;
        let parsed = UncertaintyStudy::parse(design_source()).unwrap();
        assert_eq!(parsed.design().unwrap().control, DesignControl::PowerScheduleScale { tolerance: 0.001 });
        assert_eq!(UncertaintyStudy::parse(parsed.canonical()).unwrap(), parsed);
        for (from, to) in [(":minimum 0", ":minimum -1"), (":maximum 2", ":maximum 2W"),
            (":tolerance 0.001", ":tolerance 0.001W"),
            ("nondecreasing-temperature-with-schedule-scale", "nondecreasing-temperature-with-power")] {
            assert!(UncertaintyStudy::parse(&design_source().replace(from, to)).is_err());
        }
        let mut zero = base();
        for step in &mut zero.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .transient.as_mut().unwrap().power_schedules[0].steps { step.watts.value = 0.0; }
        assert!(parsed.bind(&zero).is_err(), "an all-zero schedule is not a workload sizing coordinate");
    }

    #[test]
    fn g1_design_cases_share_the_original_schedule_and_do_not_reapply_duty() {
        let original = base();
        let source = design_source().trim_end_matches(')').to_owned() + r#"
          :operating-cases (
            (case :name "cool" :inputs ((convection-temperature :entity "air" :value 295K)))
            (case :name "hot" :inputs ((convection-temperature :entity "air" :value 303K)))))"#;
        let bound = UncertaintyStudy::parse(&source).unwrap().bind(&original).unwrap();
        for (amplitude, temperature) in [(2.0, 303.0), (0.0, 295.0), (1.0, 303.0), (2.0, 303.0)] {
            let candidate = bound.sample_project(&[amplitude, temperature]).unwrap();
            assert_eq!(candidate.power, original.power);
            assert_eq!(schedule(&candidate).steps[0].watts.value, 8.0 * amplitude);
            assert_eq!(schedule(&candidate).steps[1].watts.value, 0.0);
            assert_eq!(schedule(&candidate).steps[2].watts.value, 2.0 * amplitude);
            assert_eq!(bound.base(), &original);
        }
        assert!(bound.sample_project(&[1.0, 299.0]).is_err(), "no unauthored intermediate case");
    }

    #[test]
    fn g0_operating_case_can_disable_a_schedule_without_overwriting_the_search_control() {
        let source = design_source().replace(
            "power-schedule-scale :entity \"air\" :minimum 0 :maximum 2 :tolerance 0.001",
            "fan-speed-ratio :entity \"fixture-bank\" :minimum 0.6 :maximum 1.4 :tolerance 0.001")
            .replace("nondecreasing-temperature-with-schedule-scale", "nonincreasing-temperature-with-fan-speed");
        let source = source.trim_end_matches(')').to_owned() + r#"
          :operating-cases (
            (case :name "off" :inputs ((power-schedule-scale :entity "air" :value 0)))
            (case :name "on" :inputs ((power-schedule-scale :entity "air" :value 2)))))"#;
        let parsed = UncertaintyStudy::parse(&source).unwrap();
        assert_eq!(parsed.operating_cases()[0].values, [0.0]);
        assert_eq!(parsed.operating_cases()[1].values, [2.0]);
        let bad = source.replace(":value 0", ":value -1");
        assert!(UncertaintyStudy::parse(&bad).is_err());
    }

}
