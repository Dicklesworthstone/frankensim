//! Finite authored operating cases for native design. These values are not
//! sampled, independent marginals, probability laws, or an uncertainty box.
use super::*;

/// One complete operating condition. Values follow study.parameters()[1..];
/// the first parameter is the independently searched design coordinate.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignOperatingCase {
    pub name: String,
    pub values: Vec<f64>,
}

pub(super) fn parse(node: &Node, control: &UniformParameter)
    -> Result<(Vec<UniformParameter>, Vec<DesignOperatingCase>)>
{
    let declarations = list(node)?;
    if declarations.is_empty() || declarations.len() > 16 {
        return Err(error("declare 1..=16 finite operating cases"));
    }
    let mut names = BTreeSet::new();
    let mut parameters = Vec::<UniformParameter>::new();
    let mut cases = Vec::with_capacity(declarations.len());
    for declaration in declarations {
        let n = list(declaration)?;
        symbol(n.first().ok_or_else(|| error("empty operating case"))?, "case")?;
        let f = fields(&n[1..], &["name", "inputs"])?;
        let name = text(f["name"])?;
        if name.is_empty() || !names.insert(name.clone()) {
            return Err(error("operating case names must be nonempty and distinct"));
        }
        let inputs = list(f["inputs"])?;
        if inputs.is_empty() || inputs.len() >= 32 {
            return Err(error("each operating case needs 1..=31 physical coordinates"));
        }
        let mut values = BTreeMap::new();
        for input in inputs {
            let n = list(input)?;
            let kind = n.first().ok_or_else(|| error("empty operating-case input"))?;
            let target = if matches!(&kind.kind, NodeKind::Symbol(name) if name == "regional-power") {
                Target::Power
            } else { parse_target(kind)? };
            let f = fields(&n[1..], &["entity", "value"])?;
            let entity = text(f["entity"])?;
            if target == control.target && entity == control.entity {
                return Err(error("an operating case cannot overwrite the design control"));
            }
            let value = quantity(f["value"], target.dims())?;
            if (target == Target::Power && value < 0.0)
                || (target != Target::Power && target != Target::HeatFlux && value <= 0.0)
            {
                return Err(error("invalid operating-case physical value"));
            }
            if values.insert((target, entity.clone()), value).is_some() {
                return Err(error("duplicate operating-case physical target"));
            }
            if cases.is_empty() {
                parameters.push(UniformParameter {
                    name: format!("operating-coordinate-{}", parameters.len() + 1),
                    target, entity, low: value, high: value,
                });
            }
        }
        if values.len() != parameters.len()
            || parameters.iter().any(|p| !values.contains_key(&(p.target, p.entity.clone())))
        {
            return Err(error("every operating case must supply the same physical coordinates"));
        }
        let values = parameters.iter_mut().map(|parameter| {
            let value = values[&(parameter.target, parameter.entity.clone())];
            parameter.low = parameter.low.min(value);
            parameter.high = parameter.high.max(value);
            value
        }).collect();
        cases.push(DesignOperatingCase { name, values });
    }
    Ok((parameters, cases))
}

#[cfg(test)]
mod tests {
    use super::*;
    const SOURCE: &str = r#"(fsim-native-design-study
      :version 1 :project "cooling-reference.fsim" :evaluations 32 :seed 29
      :wall-time 120s :qoi "temperature-max"
      :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
      :materials ("aa6061.fsmcdpk") :interfaces ()
      :control (regional-power :entity "air" :minimum 0W :maximum 20W :tolerance 0.01W)
      :temperature-tolerance 0.01K :assumption nondecreasing-temperature-with-power
      :operating-cases (
        (case :name "cool" :inputs (
          (convection-temperature :entity "air" :value 295K)
          (convection-coefficient :entity "air" :value 12kg·s^-3·K^-1)))
        (case :name "hot" :inputs (
          (convection-coefficient :entity "air" :value 8kg·s^-3·K^-1)
          (convection-temperature :entity "air" :value 303K)))))"#;

    #[test]
    fn cases_preserve_complete_authored_combinations_and_physical_coordinates() {
        let study = UncertaintyStudy::parse(SOURCE).unwrap();
        assert_eq!(study.parameters().len(), 3);
        assert_eq!(study.operating_cases()[0].values, [295.0, 12.0]);
        assert_eq!(study.operating_cases()[1].values, [303.0, 8.0],
            "input order may differ without changing the case");
        assert_eq!(study.parameters()[1].low, 295.0);
        assert_eq!(study.parameters()[1].high, 303.0);
        assert!(study.compliance().is_none());
        assert!(study.latent_correlation().is_none());
        let mut base = crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-reference.fsim"))).unwrap().decoded.spec;
        base.power.as_mut().unwrap()[0].duty = 0.5;
        let bound = study.bind(&base).unwrap();
        assert!(bound.sample_project(&[7.0, 303.0, 12.0]).is_err(),
            "marginal extrema do not create an unauthored operating case");
        let changed = bound.sample_project(&[7.0, 303.0, 8.0]).unwrap();
        assert_eq!(changed.power.as_ref().unwrap()[0].watts.value, 14.0);
        assert_eq!(changed.power.as_ref().unwrap()[0].duty, 0.5);
        let mut expected = base.clone();
        expected.power.as_mut().unwrap()[0].watts.value = 14.0;
        expected.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::Convection {
                coefficient: QtyAny::new(8.0, Target::ConvectionCoefficient.dims()),
                reference_temperature: QtyAny::new(303.0, Target::ConvectionTemperature.dims()),
            };
        assert_eq!(changed, expected);
        assert_eq!(bound.base(), &base);
    }

    #[test]
    fn cases_refuse_aliases_missing_values_units_and_probability_extensions() {
        for source in [
            SOURCE.replace(":name \"hot\"", ":name \"cool\""),
            SOURCE.replace("(convection-temperature :entity \"air\" :value 303K)", ""),
            SOURCE.replace("(convection-temperature :entity \"air\" :value 303K)",
                "(regional-power :entity \"air\" :value 3W)"),
            SOURCE.replace(":value 303K", ":value 303W"),
            SOURCE.replace(":value 8kg·s^-3·K^-1", ":value 0kg·s^-3·K^-1"),
            SOURCE.replace(":wall-time 120s", ":wall-time 120s :correlation independent"),
            SOURCE.replace("(convection-temperature :entity \"air\" :value 303K)",
                "(convection-coefficient :entity \"air\" :value 10kg·s^-3·K^-1)"),
        ] {
            assert!(UncertaintyStudy::parse(&source).is_err(), "{source}");
        }
        let base = crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-reference.fsim"))).unwrap().decoded.spec;
        assert!(UncertaintyStudy::parse(&SOURCE.replace(":value 303K", ":value 330K"))
            .unwrap().bind(&base).is_err(), "cases cannot expand the original envelope");
    }

    #[test]
    fn fan_cases_keep_delivered_power_under_unchanged_duty() {
        let source = SOURCE.replace(
            "regional-power :entity \"air\" :minimum 0W :maximum 20W :tolerance 0.01W",
            "fan-speed-ratio :entity \"fixture-bank\" :minimum 0.8 :maximum 1.2 :tolerance 0.01")
            .replace("nondecreasing-temperature-with-power", "nonincreasing-temperature-with-fan-speed")
            .replace("convection-temperature", "air-inlet-temperature")
            .replace("convection-coefficient", "regional-power")
            .replace("12kg·s^-3·K^-1", "12W").replace("8kg·s^-3·K^-1", "8W");
        let mut base = crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../data/reference-project/cooling-reference.fsim"))).unwrap().decoded.spec;
        base.power.as_mut().unwrap()[0].duty = 0.5;
        base.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::AirflowConvection {
                branch: "air".into(), order: 0,
                inlet_temperature: QtyAny::new(295.0, crate::spec::dims::TEMPERATURE),
                hydraulic_diameter: QtyAny::new(0.02, crate::spec::dims::LENGTH),
                flow_area: QtyAny::new(0.004, crate::spec::dims::AREA),
                channel_length: QtyAny::new(0.3, crate::spec::dims::LENGTH),
                correlation: "convection.gnielinski".into(),
            };
        let bound = UncertaintyStudy::parse(&source).unwrap().bind(&base).unwrap();
        let candidate = bound.sample_project(&[1.1, 303.0, 8.0]).unwrap();
        assert_eq!(candidate.power.as_ref().unwrap()[0].watts.value, 16.0);
        assert_eq!(candidate.cooling.as_ref().unwrap().fan_system.as_ref().unwrap()
            .banks[0].speed_ratio, 1.1);
    }
}
