//! Sample only the declared nominal manufactured coordinate. Manufacturing
//! half-widths remain source-support requirements, never inferred distributions.
use super::{ProjectSpec, Result, Target, UniformParameter, error};
use crate::InterfaceState;

pub(super) fn apply(project: &mut ProjectSpec, parameter: &UniformParameter, value: f64) -> Result<()> {
    let temperature = project.envelope.as_ref().ok_or_else(|| error("missing operating envelope"))?.ambient_lo.value;
    let bindings = project.interface_cards.as_mut().ok_or_else(|| error("contact input requires interface-card bindings"))?;
    let mut matches = bindings.iter_mut().filter(|b| b.interface == parameter.entity);
    let binding = matches.next().ok_or_else(|| error("contact input names no interface-card binding"))?;
    if matches.next().is_some() { return Err(error("contact input names an ambiguous interface")); }
    let coordinate = match (&mut binding.state, parameter.target) {
        (InterfaceState::DryContact { pressure, .. }, Target::ContactPressure) => pressure,
        (InterfaceState::Tim { thickness, .. } | InterfaceState::Adhesive { thickness, .. }, Target::ContactThickness) => thickness,
        (InterfaceState::GapWithFluid { gap, .. }, Target::ContactGap) => gap,
        (InterfaceState::BoltedWithPattern { torque, .. }, Target::ContactTorque) => torque,
        _ => return Err(error("contact probability target does not match the declared manufactured joint class")),
    };
    coordinate.value = value;
    // Preserve the original uncertainty half-width; it must still form a
    // physical support band around each newly sampled nominal coordinate.
    crate::interface_state::query_point(&binding.state, "T", temperature,
        crate::interface_state::StatePoint::Nominal)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::UncertaintyStudy;
    use fs_qty::QtyAny;
    use crate::spec::dims;
    const STUDY: &str = r#"(fsim-uncertainty-study :version 3 :project "contact.fsim"
      :samples 4 :seed 29 :wall-time 120s :method monte-carlo :correlation independent
      :qoi "temperature-max" :geometry (
       (mesh :role "cold-body" :path "cold-body.stl" :unit "m" :max-hole-edges 0)
       (mesh :role "hot-body" :path "hot-body.stl" :unit "m" :max-hole-edges 0))
      :materials () :interfaces () :mean-control (nominal-adjoint :max-solves 1)
      :parameters ((uniform :name "joint" :target contact-pressure :entity "cold-hot-joint" :low 500000Pa :high 1500000Pa)))"#;
    fn base() -> ProjectSpec {
        crate::parse_sexpr_migrating(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../examples/contact-pair/contact-pair.fsim"))).unwrap().decoded.spec
    }
    #[test]
    fn contact_sampling_preserves_cards_half_widths_and_all_unrelated_inputs() {
        let length = |v| QtyAny::new(v, dims::LENGTH);
        let torque = |v| QtyAny::new(v, fs_qty::Dims([2,1,-2,0,0,0]));
        for (state, target, unit, low, high) in [
            (base().interface_cards.unwrap()[0].state.clone(), "contact-pressure", "Pa", 500000.0_f64, 1500000.0_f64),
            (InterfaceState::Tim { thickness:length(0.003),thickness_half_width:length(0.0001) }, "contact-thickness", "m", 0.002, 0.004),
            (InterfaceState::Adhesive { thickness:length(0.003),thickness_half_width:length(0.0001) }, "contact-thickness", "m", 0.002, 0.004),
            (InterfaceState::GapWithFluid { gap:length(0.003),gap_half_width:length(0.0001),fluid:"air".into() }, "contact-gap", "m", 0.002, 0.004),
            (InterfaceState::BoltedWithPattern { torque:torque(4.0),torque_half_width:torque(0.1),bolt_count:4,pattern:"corners".into() }, "contact-torque", "N*m", 3.0, 5.0),
        ] {
            let text = STUDY.replace("contact-pressure",target)
                .replace("500000Pa",&format!("{low}{unit}"))
                .replace("1500000Pa",&format!("{high}{unit}"));
            let study = UncertaintyStudy::parse(&text).unwrap();
            assert_eq!(study, UncertaintyStudy::parse(study.canonical()).unwrap());
            let mut original = base(); original.interface_cards.as_mut().unwrap()[0].state = state;
            let bound = study.bind(&original).unwrap();
            assert_eq!(bound.study().parameters()[0].target.unit(),unit);
            assert_eq!(bound.study().mean_control().unwrap().probe_count(bound.study().parameters()),1);
            for value in [low, low.midpoint(high), high, low] {
                let sample = bound.sample_project(&[value]).unwrap();
                let mut expected = original.clone();
                match &mut expected.interface_cards.as_mut().unwrap()[0].state {
                    InterfaceState::DryContact { pressure,.. } => pressure.value=value,
                    InterfaceState::Tim { thickness,.. } | InterfaceState::Adhesive { thickness,.. } => thickness.value=value,
                    InterfaceState::GapWithFluid { gap,.. } => gap.value=value,
                    InterfaceState::BoltedWithPattern { torque,.. } => torque.value=value,
                }
                assert_eq!(sample, expected);
            }
            assert_eq!(bound.base(),&original);
            assert!(bound.sample_project(&[f64::NAN]).is_err());
            assert!(bound.sample_project(&[high * 2.0]).is_err());
        }
    }
    #[test]
    fn invalid_units_ownership_class_and_nonpositive_manufacturing_bands_refuse() {
        for text in [STUDY.replace("500000Pa","500000K"),STUDY.replace("500000Pa","0Pa")] {
            assert!(UncertaintyStudy::parse(&text).is_err());
        }
        let study = || UncertaintyStudy::parse(STUDY).unwrap();
        assert!(UncertaintyStudy::parse(&STUDY.replace("cold-hot-joint","missing")).unwrap().bind(&base()).is_err());
        let mut wrong=base();
        wrong.interface_cards.as_mut().unwrap()[0].state=InterfaceState::Tim {
            thickness:QtyAny::new(0.003,dims::LENGTH),thickness_half_width:QtyAny::new(0.0001,dims::LENGTH) };
        assert!(study().bind(&wrong).is_err());
        // Positive random support can still be invalid once the UNCHANGED
        // manufacturing half-width is applied around a sample.
        assert!(UncertaintyStudy::parse(&STUDY.replace("500000Pa","50000Pa")).unwrap().bind(&base()).is_err());
    }
}
