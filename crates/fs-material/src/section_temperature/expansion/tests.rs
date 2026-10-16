use super::*;
use crate::section::{ConcentricMaterialLayer, thermal::AxialConstraint};
use crate::state_point::{INVERSE_TEMPERATURE_DIMS, ScalarAdmissibility, ScalarPropertyRequirement};
use fs_evidence::ValidityDomain;
use fs_matdb::{MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    InterpolationPolicy, Provenance, UncertaintyModel};

pub(crate) fn card(alpha: f64) -> (MaterialCard, Vec<ScalarPropertyRequirement>) {
    let (base, mut requirements) = crate::section_temperature::tests::card(1.0);
    let mut claims = base.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ALPHA, INVERSE_TEMPERATURE_DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, alpha), (400.0, 2.0*alpha)], dims: INVERSE_TEMPERATURE_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0)
            .with("pressure", 100000.0, 100000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "synthetic expansion, not a measured material".into(),
            license: "test-only".into(), artifact: None },
    }).unwrap();
    requirements.push(ScalarPropertyRequirement::try_new(ALPHA, INVERSE_TEMPERATURE_DIMS,
        ScalarAdmissibility::Finite).unwrap());
    (MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "expansion test".into(), revision: 0 }, claims, vec![]).unwrap(), requirements)
}
fn close(a: f64, b: f64) { assert!((a-b).abs() <= 1e-11*a.abs().max(b.abs()).max(1.0), "{a} != {b}"); }
fn section(cards: &[&MaterialCard], reqs: &[Vec<ScalarPropertyRequirement>]) -> ConcentricElasticSection {
    let layers = cards.iter().zip(reqs).enumerate().map(|(i,(card,req))| {
        let resolved = crate::section_temperature::tests::section(card, req);
        ConcentricMaterialLayer { outer_radius_m: (i+1) as f64*0.001,
            material: resolved.layers()[0].material.clone() }
    }).collect();
    ConcentricElasticSection::try_new(0.0, layers).unwrap()
}

#[test]
fn integrates_alpha_instead_of_multiplying_its_endpoint_value() {
    let (a,ra) = card(1e-5); let (b,rb) = card(2e-5);
    let source = section(&[&a,&b], &[ra,rb]);
    let law = ConcentricSectionTemperatureLaw::try_new(source.clone(), &[&a,&b], &["T","T"]).unwrap();
    let hot = law.thermoelastic_at_temperatures_k(&[350.0,400.0]).unwrap();
    let strains = [0.000625, 0.003];
    for (path, expected) in hot.expansion_paths().iter().zip(strains) {
        close(path.free_linear_strain(), expected);
    }
    let held = hot.response_with_axial_constraint(AxialConstraint::FixedStrain(0.002),0.01).unwrap();
    let expected: f64 = hot.section().layer_integrals().iter().zip(strains)
        .map(|(layer,eth)| layer.axial_stiffness_n()*(0.002-eth)).sum();
    close(held.force_n(), expected);
    assert!(held.layers()[0].force_n()>0.0 && held.layers()[1].force_n()<0.0);
    assert_eq!(law.reference(), &source);
    assert_eq!(hot.section().layers()[0].material.property(ALPHA).unwrap().answer().receipt.selected,
        source.layers()[0].material.property(ALPHA).unwrap().answer().receipt.selected);
}

#[test]
fn free_force_and_elastic_supports_share_the_same_layer_law() {
    let (a,ra) = card(1e-5); let (b,rb) = card(2e-5);
    let law = ConcentricSectionTemperatureLaw::try_new(section(&[&a,&b], &[ra,rb]), &[&a,&b], &["T","T"]).unwrap();
    let hot = law.thermoelastic_at_temperatures_k(&[400.0,400.0]).unwrap();
    let free = hot.response_with_axial_constraint(AxialConstraint::FixedForce(0.0),0.02).unwrap();
    close(free.force_n(),0.0);
    assert!(free.layers().iter().all(|l| l.stress_pa()!=0.0));
    assert!(free.strain_energy_per_length_j_m()>0.0);
    let force = hot.response_with_axial_constraint(AxialConstraint::FixedForce(20.0),0.02).unwrap();
    close(force.force_n(),20.0);
    for stiffness in [0.0, 1e3, 1e8, 1e12] {
        let spring = hot.response_with_axial_constraint(AxialConstraint::LinearSpring {
            stiffness_n_m: stiffness, reference_length_m: 0.5, anchor_strain: 0.001, preload_n: 20.0,
        },0.02).unwrap();
        // A very stiff support subtracts nearly equal strains; compare on the
        // force scale of that cancellation, not an invented exact equality.
        let expected = 20.0 + stiffness*0.5*(0.001-spring.total_strain());
        assert!((spring.force_n()-expected).abs()<1e-7*spring.force_n().abs().max(1.0));
        if stiffness==0.0 { close(spring.total_strain(),force.total_strain()); }
        assert!(spring.total_strain()>=0.001 && spring.total_strain()<=force.total_strain());
    }
}

#[test]
fn contraction_cooling_and_reference_replay_do_not_inherit_a_positive_alpha() {
    let (a,req) = card(-1e-5);
    let reference = crate::section_temperature::tests::section(&a,&req);
    let law = ConcentricSectionTemperatureLaw::try_new(reference, &[&a,&a], &["T","T"]).unwrap();
    let hot = law.thermoelastic_at_temperatures_k(&[400.0,400.0]).unwrap();
    assert!(hot.response_with_axial_constraint(AxialConstraint::FixedStrain(0.0),0.01).unwrap().force_n()>0.0);
    let warm = law.thermoelastic_at_temperatures_k(&[350.0,350.0]).unwrap();
    close(warm.expansion_paths()[0].free_linear_strain(),-0.000625);
    let cold = law.thermoelastic_at_temperatures_k(&[300.0,300.0]).unwrap();
    assert_eq!(cold.response_with_axial_constraint(AxialConstraint::FixedStrain(0.0),0.01).unwrap().force_n(),0.0);
}

#[test]
fn missing_expansion_invalid_supports_and_excessive_layer_strain_refuse() {
    let (plain,req) = crate::section_temperature::tests::card(1.0);
    let law = ConcentricSectionTemperatureLaw::try_new(crate::section_temperature::tests::section(&plain,&req),
        &[&plain,&plain], &["T","T"]).unwrap();
    assert!(law.thermoelastic_at_temperatures_k(&[300.0,300.0]).is_err());
    let (a,req) = card(1e-5);
    let law = ConcentricSectionTemperatureLaw::try_new(crate::section_temperature::tests::section(&a,&req),
        &[&a,&a], &["T","T"]).unwrap();
    let hot = law.thermoelastic_at_temperatures_k(&[400.0,400.0]).unwrap();
    for constraint in [AxialConstraint::FixedForce(f64::NAN), AxialConstraint::FixedStrain(f64::INFINITY),
        AxialConstraint::LinearSpring { stiffness_n_m:-1.0, reference_length_m:1.0,anchor_strain:0.0,preload_n:0.0 },
        AxialConstraint::LinearSpring { stiffness_n_m:1.0, reference_length_m:0.0,anchor_strain:0.0,preload_n:0.0 }] {
        assert!(hot.response_with_axial_constraint(constraint,0.01).is_err());
    }
    // Small net force cannot waive each layer's thermal-strain limit.
    assert!(hot.response_with_axial_constraint(AxialConstraint::FixedForce(0.0),1e-4).is_err());
    assert!(law.thermoelastic_at_temperatures_k(&[300.0,401.0]).is_err());
}
