//! Choose the existing native derivative producer before calibration/recovery.
use super::{OUTPUT, Result, Target, UniformParameter, invalid};
use fs_project::{ProjectSpec, spec::OutputRequest};

const BOUNDARY: &str = "temperature-max-boundary-adjoint";
const CONTACT: &str = "temperature-max-contact-adjoint";

/// Called only on the calibration clone. Probability children and original
/// output declarations remain unchanged. Repeating the choice is idempotent.
pub(super) fn configure(spec: &mut ProjectSpec, parameters: &[UniformParameter]) -> Result<()> {
    let boundary = parameters.iter().any(|p| p.target == Target::FixedTemperature && p.low != p.high);
    let outputs = spec.outputs.get_or_insert_with(Vec::new);
    let requests: Vec<_> = outputs.iter().enumerate()
        .filter(|(_, row)| matches!(row.name.as_str(), OUTPUT | BOUNDARY | CONTACT))
        .map(|(i, _)| i).collect();
    if requests.len() > 1 || requests.first().is_some_and(|&i| outputs[i].kind != "report") {
        return Err(invalid("nominal calibration requires one unambiguous native adjoint report"));
    }
    if let Some(&i) = requests.first() {
        if boundary && outputs[i].name == CONTACT {
            return Err(invalid("fixed-temperature calibration needs temperature-max-boundary-adjoint, not the mutually exclusive contact report"));
        }
        if boundary { outputs[i].name = BOUNDARY.into(); }
    } else {
        outputs.push(OutputRequest { name: if boundary { BOUNDARY } else { OUTPUT }.into(),
            kind: "report".into(), region: None });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn base() -> ProjectSpec {
        fs_project::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../examples/contact-pair/contact-pair.fsim"))).unwrap().decoded.spec
    }
    fn parameter(constant: bool) -> UniformParameter {
        UniformParameter { name: "wall".into(), target: Target::FixedTemperature,
            entity: "cold".into(), low: 280.0, high: if constant {280.0} else {290.0} }
    }
    fn request(name: &str) -> OutputRequest {
        OutputRequest { name: name.into(), kind: "report".into(), region: None }
    }
    #[test]
    fn prescribed_calibration_changes_only_the_report_and_repeats_exactly() {
        for existing in [None, Some(OUTPUT), Some(BOUNDARY)] {
            let mut project = base();
            if let Some(name) = existing { project.outputs.as_mut().unwrap().push(request(name)); }
            let mut expected = project.clone();
            if existing.is_some() { expected.outputs.as_mut().unwrap().last_mut().unwrap().name = BOUNDARY.into(); }
            else { expected.outputs.as_mut().unwrap().push(request(BOUNDARY)); }
            configure(&mut project, &[parameter(false)]).unwrap();
            assert_eq!(project, expected, "physical inputs and unrelated outputs must not change");
            configure(&mut project, &[parameter(false)]).unwrap();
            assert_eq!(project, expected, "recovery must choose the same canonical project");
        }
    }
    #[test]
    fn existing_extensions_are_reused_and_constant_fixed_inputs_need_no_lift() {
        for name in [OUTPUT, BOUNDARY, CONTACT] {
            let mut project = base(); project.outputs.as_mut().unwrap().push(request(name));
            let before = project.clone();
            configure(&mut project, &[parameter(true)]).unwrap();
            assert_eq!(project, before, "never append a conflicting standard report");
        }
        let mut project = base();
        configure(&mut project, &[parameter(true)]).unwrap();
        assert_eq!(project.outputs.unwrap().last().unwrap().name, OUTPUT);
    }
    #[test]
    fn ambiguous_or_incompatible_report_intent_refuses_without_replacing_it() {
        for names in [vec![CONTACT], vec![OUTPUT, BOUNDARY], vec![BOUNDARY, BOUNDARY]] {
            let mut project = base();
            project.outputs.as_mut().unwrap().extend(names.into_iter().map(request));
            let before = project.clone();
            assert!(configure(&mut project, &[parameter(false)]).is_err());
            assert_eq!(project, before);
        }
        let mut project = base();
        project.outputs.as_mut().unwrap().push(OutputRequest { kind: "scalar".into(), ..request(BOUNDARY) });
        assert!(configure(&mut project, &[]).is_err());
    }
}
