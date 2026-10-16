//! Choose the existing native derivative producer before calibration/recovery.
use super::{OUTPUT, Result, Target, UniformParameter, invalid};
use fs_project::{ProjectSpec, spec::OutputRequest};

const BOUNDARY: &str = "temperature-max-boundary-adjoint";
pub(super) const CONTACT: &str = "temperature-max-contact-adjoint";
pub(super) const COMBINED: &str = "temperature-max-contact-boundary-adjoint";

/// Called only on the calibration clone. Probability children and original
/// output declarations remain unchanged. Repeating the choice is idempotent.
pub(super) fn configure(spec: &mut ProjectSpec, parameters: &[UniformParameter]) -> Result<()> {
    let boundary = parameters.iter().any(|p| p.target == Target::FixedTemperature && p.low != p.high);
    let contact = parameters.iter().any(|p| p.target.contact_axis().is_some() && p.low != p.high);
    let outputs = spec.outputs.get_or_insert_with(Vec::new);
    let requests: Vec<_> = outputs.iter().enumerate()
        .filter(|(_, row)| matches!(row.name.as_str(), OUTPUT | BOUNDARY | CONTACT | COMBINED))
        .map(|(i, _)| i).collect();
    if requests.len() > 1 || requests.first().is_some_and(|&i| outputs[i].kind != "report") {
        return Err(invalid("nominal calibration requires one unambiguous native adjoint report"));
    }
    // Keep explicitly requested control families when adding those required by
    // calibration. Combining their contractions does not add a second solve.
    let existing = requests.first().map(|&i| outputs[i].name.as_str());
    let contact = contact || matches!(existing, Some(CONTACT | COMBINED));
    let boundary = boundary || matches!(existing, Some(BOUNDARY | COMBINED));
    let desired = match (contact, boundary) {
        (true, true) => COMBINED,
        (true, false) => CONTACT,
        (false, true) => BOUNDARY,
        (false, false) => OUTPUT,
    };
    if let Some(&i) = requests.first() {
        outputs[i].name = desired.into();
    } else {
        outputs.push(OutputRequest { name: desired.into(),
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
        for name in [OUTPUT, BOUNDARY, CONTACT, COMBINED] {
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
        for names in [vec![OUTPUT, BOUNDARY], vec![BOUNDARY, BOUNDARY],
            vec![CONTACT, COMBINED], vec![COMBINED, COMBINED]] {
            let mut project = base();
            project.outputs.as_mut().unwrap().extend(names.into_iter().map(request));
            let before = project.clone();
            assert!(configure(&mut project, &[parameter(false)]).is_err());
            assert_eq!(project, before);
        }
        for name in [OUTPUT, CONTACT, BOUNDARY, COMBINED] {
            let mut project = base();
            project.outputs.as_mut().unwrap().push(OutputRequest { kind: "scalar".into(), ..request(name) });
            let before = project.clone();
            assert!(configure(&mut project, &[]).is_err());
            assert_eq!(project, before);
        }
    }
    #[test]
    fn joint_calibration_selects_contact_only_and_preserves_physical_inputs() {
        let contact = UniformParameter { name:"clamp".into(),entity:"cold-hot-joint".into(),
            target:Target::ContactPressure,low:500000.0,high:1500000.0 };
        for existing in [None,Some(OUTPUT),Some(CONTACT)] {
            let mut project=base();
            if let Some(name)=existing {project.outputs.as_mut().unwrap().push(request(name));}
            let mut expected=project.clone();
            if existing.is_some() {expected.outputs.as_mut().unwrap().last_mut().unwrap().name=CONTACT.into();}
            else {expected.outputs.as_mut().unwrap().push(request(CONTACT));}
            configure(&mut project,std::slice::from_ref(&contact)).unwrap();
            assert_eq!(project,expected);
            configure(&mut project,std::slice::from_ref(&contact)).unwrap();
            assert_eq!(project,expected);
        }
        let mut project=base();
        let mut constant=contact;constant.high=constant.low;
        configure(&mut project,&[constant]).unwrap();
        assert_eq!(project.outputs.unwrap().last().unwrap().name,OUTPUT);
    }
    #[test]
    fn joint_and_prescribed_controls_select_one_combined_report_without_erasing_intent() {
        for existing in [None, Some(OUTPUT), Some(CONTACT), Some(BOUNDARY), Some(COMBINED)] {
            for (vary_contact, vary_wall) in [(true, true), (true, false), (false, true), (false, false)] {
                let contact = UniformParameter { name:"clamp".into(),entity:"cold-hot-joint".into(),
                    target:Target::ContactPressure,low:500000.0,
                    high:if vary_contact {1500000.0} else {500000.0} };
                let parameters = [contact, parameter(!vary_wall)];
                let mut project = base();
                if let Some(name) = existing { project.outputs.as_mut().unwrap().push(request(name)); }
                let want_contact = vary_contact || matches!(existing, Some(CONTACT | COMBINED));
                let want_wall = vary_wall || matches!(existing, Some(BOUNDARY | COMBINED));
                let name = match (want_contact, want_wall) {
                    (true, true) => COMBINED, (true, false) => CONTACT,
                    (false, true) => BOUNDARY, (false, false) => OUTPUT,
                };
                let mut expected = project.clone();
                if existing.is_some() { expected.outputs.as_mut().unwrap().last_mut().unwrap().name = name.into(); }
                else { expected.outputs.as_mut().unwrap().push(request(name)); }
                configure(&mut project, &parameters).unwrap();
                assert_eq!(project, expected, "only the necessary report may change");
                configure(&mut project, &parameters).unwrap();
                assert_eq!(project, expected, "recovery must select identical calibration source");
            }
        }
    }
}
