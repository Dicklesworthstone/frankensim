//! G0/G1/G3: one material card drives the existing Maxwell and heat-feedback laws.
//! All source data are synthetic; these tests do not qualify a real specimen.
use core::convert::Infallible;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStep,
};
use fs_material::thermovisco::source::{
    ResolvedThermoMaxwell, SourceError, ThermoMaxwellSourceQuery,
};
use fs_material::thermovisco::thermal::{FeedbackConfig, MechanicalControl};
use fs_material::thermovisco::{HistoryIntegration, Limits, ThermalKnot, ThermoMaxwell};
use fs_qty::semantic::{
    FrequencyConvention, QuantityKind, SemanticType, StrainBasis, StrainComponent, ValueForm,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature, Time};

const ONE: Dims = Dims([0; 6]);
fn limits() -> Limits {
    Limits {
        max_strain: 0.04,
        max_strain_rate_s: 10.0,
    }
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
fn add(claims: &mut ClaimSet, name: &str, dims: Dims, knots: Vec<(f64, f64)>) {
    let lo = knots[0].0;
    let hi = knots.last().unwrap().0;
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims,
            },
            validity: ValidityDomain::unconstrained().with("T", lo, hi),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic Maxwell source".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
}
fn claims(cp: f64) -> ClaimSet {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("E_inf", Pressure::DIMS, 1e8),
        ("E0", Pressure::DIMS, 3e8),
        ("tau0", Time::DIMS, 1.0),
        ("scale", ONE, 1.0),
        ("free", ONE, 0.0),
        ("density", Density::DIMS, 1000.0),
        ("cp", Dims([2, 0, -2, -1, 0, 0]), cp),
    ] {
        add(
            &mut claims,
            name,
            dims,
            vec![(280.0, value), (320.0, value)],
        );
    }
    add(
        &mut claims,
        "shift",
        ONE,
        vec![(280.0, 1.0), (300.0, 1.0), (310.0, 0.4), (320.0, 0.05)],
    );
    claims
}
fn card(claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "source coupling test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(
    card: &MaterialCard,
    lower: &QueryPoint,
    upper: &QueryPoint,
) -> Result<ResolvedThermoMaxwell, SourceError> {
    ResolvedThermoMaxwell::try_new(
        card,
        ThermoMaxwellSourceQuery {
            equilibrium_modulus: &PropertyKey::new("E_inf", Pressure::DIMS),
            branches: &[(
                PropertyKey::new("E0", Pressure::DIMS),
                PropertyKey::new("tau0", Time::DIMS),
            )],
            modulus_scale: &PropertyKey::new("scale", ONE),
            relaxation_scale: &PropertyKey::new("shift", ONE),
            free_strain: &PropertyKey::new("free", ONE),
            temperature_axis: "T",
            lower,
            upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        limits(),
        || false,
    )
}
fn caloric(
    card: &MaterialCard,
    cp: f64,
    lower: &QueryPoint,
    upper: &QueryPoint,
) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("cp", Dims([2, 0, -2, -1, 0, 0])),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower,
            upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * cp,
        1e-7,
    )
    .unwrap()
}

#[test]
fn all_source_breakpoints_drive_the_original_history_update() {
    let card = card(claims(5.0));
    let source = resolve(&card, &point(280.0), &point(320.0)).unwrap();
    assert_eq!(source.samples().len(), 4);
    let direct = ThermoMaxwell::new(
        1e8,
        vec![(3e8, 1.0)],
        vec![
            ThermalKnot {
                temperature_k: 280.0,
                modulus_scale: 1.0,
                relaxation_scale: 1.0,
                free_strain: 0.0,
            },
            ThermalKnot {
                temperature_k: 300.0,
                modulus_scale: 1.0,
                relaxation_scale: 1.0,
                free_strain: 0.0,
            },
            ThermalKnot {
                temperature_k: 310.0,
                modulus_scale: 1.0,
                relaxation_scale: 0.4,
                free_strain: 0.0,
            },
            ThermalKnot {
                temperature_k: 320.0,
                modulus_scale: 1.0,
                relaxation_scale: 0.05,
                free_strain: 0.0,
            },
        ],
        limits(),
    )
    .unwrap();
    assert_eq!(source.law(), &direct);
    let mut a = source.law().birth(300.0, 0.0).unwrap();
    let mut b = direct.birth(300.0, 0.0).unwrap();
    for (t, strain) in [(307.0, 0.01), (313.0, 0.02), (303.0, 0.01), (300.0, 0.0)] {
        let sa = a.advance(t, strain, 0.2, 1e-6).unwrap();
        let sb = b.advance(t, strain, 0.2, 1e-6).unwrap();
        assert_eq!(sa.point, sb.point);
        assert_eq!(sa.stress_pa, sb.stress_pa);
        a = sa.point;
        b = sb.point;
    }
    for sample in source.samples() {
        assert_eq!(sample.card_identity(), card.content_hash());
    }
}

#[test]
fn replacing_only_caloric_source_data_changes_actual_stress_and_self_heating() {
    let mut output = Vec::new();
    for cp in [5.0, 10.0] {
        let card = card(claims(cp));
        let source = resolve(&card, &point(280.0), &point(320.0)).unwrap();
        let thermal = caloric(&card, cp, &point(280.0), &point(320.0));
        let old = source
            .caloric_point(&thermal, 300.0, 0.0, 1e-6, 0.001)
            .unwrap();
        let step = old
            .advance(
                0.02,
                0.2,
                FeedbackConfig {
                    max_iterations: 100,
                    temperature_tolerance_k: 1e-10,
                    energy_tolerance_j: 1e-9,
                    relaxation: 1.0,
                },
                |input| {
                    let h = input.initial.specific_enthalpy_j_kg()
                        + input.internal_heat_j / input.mass_kg;
                    Ok::<_, Infallible>(UniformEnthalpyStep {
                        state: input.curve.state_at_specific_enthalpy(h).unwrap(),
                        external_heat_j: 0.0,
                        energy_residual_tolerance_j: 1e-10,
                        report: (),
                    })
                },
                || false,
            )
            .unwrap();
        output.push((
            step.next.caloric_state().temperature_k(),
            step.next.mechanical().stress_pa(),
        ));
        assert!(step.combined_residual_j.abs() < 1e-9);
        assert_eq!(old.time_s(), 0.0);
    }
    assert!(output[0].0 > output[1].0);
    assert!(output[0].1 < output[1].1);
}

#[test]
fn source_selected_exact_ramps_retain_receipts_and_drive_caloric_feedback() {
    let card = card(claims(100.0));
    let source = resolve(&card, &point(280.0), &point(320.0)).unwrap();
    let keys = source.keys().to_vec();
    let samples = source.samples().len();
    let source = source.with_history_integration(HistoryIntegration::ExactStrainRamp);
    assert_eq!(source.law().history_integration(), HistoryIntegration::ExactStrainRamp);
    assert_eq!(source.keys(), keys.as_slice());
    assert_eq!(source.samples().len(), samples);
    for sample in source.samples() {
        assert_eq!(sample.card_identity(), card.content_hash());
    }
    let thermal = caloric(&card, 100.0, &point(280.0), &point(320.0));
    let old = source.caloric_point(&thermal, 300.0, 0.0, 1e-6, 0.001).unwrap();
    let step = old.advance(0.02, 0.2, FeedbackConfig {
        max_iterations: 100, temperature_tolerance_k: 1e-10,
        energy_tolerance_j: 1e-9, relaxation: 1.0,
    }, |input| {
        Ok::<_, Infallible>(UniformEnthalpyStep {
            state: input.curve.state_at_specific_enthalpy(
                input.initial.specific_enthalpy_j_kg()+input.internal_heat_j/input.mass_kg).unwrap(),
            external_heat_j: 0.0, energy_residual_tolerance_j: 1e-10, report: (),
        })
    }, || false).unwrap();
    let native = old.mechanical().advance(step.next.caloric_state().temperature_k(),
        0.02, 0.2, 1e-3).unwrap();
    assert_eq!(step.next.mechanical(), &native.point);
    assert_eq!(step.mechanical_energy.numerical_dissipation, 0.0);
    assert!(step.applied_heat_j > 0.0 && step.combined_residual_j.abs() < 1e-9);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn original_source_knots_admit_affine_spring_heat_and_refuse_a_hidden_entropy_jump() {
    for (middle, admissible) in [(1.0, true), (1.1, false)] {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("E_inf", Pressure::DIMS, 1e8), ("E0", Pressure::DIMS, 3e8),
            ("tau0", Time::DIMS, 1.0), ("shift", ONE, 1.0), ("free", ONE, 0.0),
            ("density", Density::DIMS, 1000.0),
            ("cp", Dims([2, 0, -2, -1, 0, 0]), 100.0),
        ] {
            add(&mut claims, name, dims, vec![(280.0, value), (320.0, value)]);
        }
        add(&mut claims, "scale", ONE,
            vec![(280.0, 1.1), (300.0, middle), (320.0, 0.9)]);
        let card = card(claims);
        let source = resolve(&card, &point(280.0), &point(320.0)).unwrap();
        let heat = caloric(&card, 100.0, &point(280.0), &point(320.0));
        // The constitutive law retains the actual middle source knot. The
        // thermal owner must not infer a slope using only endpoint values that
        // are shared by the affine and kinked source variants.
        assert_eq!(source.samples().len(), 3);
        let old = source.caloric_point(&heat, 300.0, 0.0, 1e-6, 0.001);
        if !admissible {
            assert!(matches!(old, Err(SourceError::Feedback(_))));
            continue;
        }
        let old = old.unwrap();
        let step = old.advance(0.02, 0.2, FeedbackConfig {
            max_iterations: 100, temperature_tolerance_k: 1e-10,
            energy_tolerance_j: 1e-9, relaxation: 1.0,
        }, |input| {
            Ok::<_, Infallible>(UniformEnthalpyStep {
                state: input.curve.state_at_specific_enthalpy(
                    input.initial.specific_enthalpy_j_kg()+input.internal_heat_j/input.mass_kg).unwrap(),
                external_heat_j: 0.0, energy_residual_tolerance_j: 1e-10, report: (),
            })
        }, || false).unwrap();
        let q = 0.02_f64/1.2;
        let xi = 0.5e-6*(1e8*0.02_f64.powi(2)+3e8*q*q);
        let base_loss = 1e-6*3e8*0.2*q*q;
        let expected_t = (0.1*300.0+2.5*base_loss)/(0.1+0.005*(base_loss+xi));
        assert!((step.next.caloric_state().temperature_k()-expected_t).abs() < 1e-9);
        assert!((step.next.mechanical_internal_energy_j()-2.5*xi).abs() < 1e-14);
        assert!(step.mechanical_energy.dissipated > 0.0);
        assert!(step.reversible_heat_j < 0.0 && step.applied_heat_j < 0.0);
        assert_eq!(step.next.caloric_curve().material_card_identity(), card.content_hash());
        assert!(core::ptr::eq(step.next.mechanical().law(), source.law()));
        assert_eq!(old.time_s(), 0.0);
    }
}

#[test]
fn interior_conflicts_and_hidden_reference_coefficient_changes_are_not_averaged() {
    let mut conflict = claims(5.0);
    add(
        &mut conflict,
        "E0",
        Pressure::DIMS,
        vec![(305.0, 2e8), (315.0, 2e8)],
    );
    assert!(resolve(&card(conflict), &point(280.0), &point(320.0)).is_err());
    // The envelope endpoints alone would falsely appear constant.
    let mut hidden = ClaimSet::new();
    for (name, dims, value) in [
        ("E_inf", Pressure::DIMS, 1e8),
        ("tau0", Time::DIMS, 1.0),
        ("scale", ONE, 1.0),
        ("shift", ONE, 1.0),
        ("free", ONE, 0.0),
    ] {
        add(
            &mut hidden,
            name,
            dims,
            vec![(280.0, value), (320.0, value)],
        );
    }
    add(
        &mut hidden,
        "E0",
        Pressure::DIMS,
        vec![(280.0, 3e8), (309.0, 2e8), (320.0, 3e8)],
    );
    assert!(matches!(
        resolve(&card(hidden), &point(280.0), &point(320.0)),
        Err(SourceError::Invalid(_))
    ));
}

#[test]
fn foreign_caloric_cards_or_fixed_contexts_cannot_be_bound() {
    let a = card(claims(5.0));
    let b = card(claims(10.0));
    let lower = point(280.0).with("p", 100000.0).unwrap();
    let upper = point(320.0).with("p", 100000.0).unwrap();
    let source = resolve(&a, &lower, &upper).unwrap();
    let same = caloric(&a, 5.0, &lower, &upper);
    assert_eq!(source.reference_density_kg_m3(&same).unwrap(), 1000.0);
    assert!(source.caloric_point(&same, 300.0, 0.0, 1e-6, 0.001).is_ok());
    let other = caloric(&b, 10.0, &lower, &upper);
    assert!(source.reference_density_kg_m3(&other).is_err());
    assert!(
        source
            .caloric_point(&other, 300.0, 0.0, 1e-6, 0.001)
            .is_err()
    );
    let wrong = caloric(
        &a,
        5.0,
        &point(280.0).with("p", 200000.0).unwrap(),
        &point(320.0).with("p", 200000.0).unwrap(),
    );
    assert!(
        source
            .caloric_point(&wrong, 300.0, 0.0, 1e-6, 0.001)
            .is_err()
    );
    assert!(resolve(&a, &lower, &point(320.0).with("p", 200000.0).unwrap()).is_err());
}

#[test]
fn source_binding_does_not_relax_caloric_physics_or_source_units() {
    let mut c = ClaimSet::new();
    for (name, dims, value) in [
        ("E_inf", Pressure::DIMS, 1e8),
        ("E0", Pressure::DIMS, 3e8),
        ("tau0", Time::DIMS, 1.0),
        ("scale", ONE, 1.0),
        ("shift", ONE, 1.0),
        ("density", Density::DIMS, 1000.0),
        ("cp", Dims([2, 0, -2, -1, 0, 0]), 5.0),
    ] {
        add(&mut c, name, dims, vec![(280.0, value), (320.0, value)]);
    }
    // This affine expansion is mechanically admitted but exceeds the complete
    // caloric capacity bound; source binding must retain that physical refusal.
    add(&mut c, "free", ONE, vec![(280.0, 0.0), (320.0, 0.04)]);
    let c = card(c);
    let source = resolve(&c, &point(280.0), &point(320.0)).unwrap();
    let thermal = caloric(&c, 5.0, &point(280.0), &point(320.0));
    assert!(
        source
            .caloric_point(&thermal, 300.0, 0.0, 1e-6, 0.001)
            .is_err()
    );
    assert!(
        source
            .law()
            .birth(280.0, 0.0)
            .unwrap()
            .advance(320.0, 0.0, 1.0, 1e-6)
            .unwrap()
            .stress_pa
            < 0.0
    );
    let wrong = PropertyKey::new("E0", Time::DIMS);
    assert!(
        ResolvedThermoMaxwell::try_new(
            &c,
            ThermoMaxwellSourceQuery {
                equilibrium_modulus: &wrong,
                branches: &[],
                modulus_scale: &PropertyKey::new("scale", ONE),
                relaxation_scale: &PropertyKey::new("shift", ONE),
                free_strain: &PropertyKey::new("free", ONE),
                temperature_axis: "T",
                lower: &point(280.0),
                upper: &point(320.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            limits(),
            || false
        )
        .is_err()
    );
}

#[test]
fn cancellation_refuses_before_publication_and_zero_branches_still_need_declared_laws() {
    let c = card(claims(5.0));
    let mut calls = 0;
    let query = ThermoMaxwellSourceQuery {
        equilibrium_modulus: &PropertyKey::new("E_inf", Pressure::DIMS),
        branches: &[],
        modulus_scale: &PropertyKey::new("scale", ONE),
        relaxation_scale: &PropertyKey::new("shift", ONE),
        free_strain: &PropertyKey::new("free", ONE),
        temperature_axis: "T",
        lower: &point(280.0),
        upper: &point(320.0),
        selection: SelectionPolicy::SingleClaimOnly,
    };
    assert!(matches!(
        ResolvedThermoMaxwell::try_new(&c, query, limits(), || {
            calls += 1;
            calls == 5
        }),
        Err(SourceError::Cancelled)
    ));
    let source = ResolvedThermoMaxwell::try_new(&c, query, limits(), || false).unwrap();
    assert_eq!(source.law().branch_count(), 0);
    let thermal = caloric(&c, 5.0, &point(280.0), &point(320.0));
    let state = source
        .caloric_point(&thermal, 300.0, 0.0, 1e-6, 0.001)
        .unwrap();
    let control = MechanicalControl::Stress {
        stress_pa: 2e6,
        tolerance_pa: 1e-6,
    };
    let step = state
        .advance_controlled(
            control,
            0.2,
            FeedbackConfig {
                max_iterations: 5,
                temperature_tolerance_k: 1e-8,
                energy_tolerance_j: 1e-9,
                relaxation: 1.0,
            },
            |input| {
                assert_eq!(input.internal_heat_j, 0.0);
                Ok::<_, Infallible>(UniformEnthalpyStep {
                    state: input.initial,
                    external_heat_j: 0.0,
                    energy_residual_tolerance_j: 1e-10,
                    report: (),
                })
            },
            || false,
        )
        .unwrap();
    assert_eq!(step.next.caloric_state().temperature_k(), 300.0);
}

#[test]
fn original_free_strain_claims_drive_signed_expansion_and_refuse_kinks_or_inadequate_capacity() {
    for (alpha,middle_offset,cp,admissible) in [
        (1e-4,0.0,5.0,true),(-1e-4,0.0,5.0,true),
        (1e-4,1e-4,5.0,false),(1e-4,0.0,0.5,false),
    ] {
        let mut claims=ClaimSet::new();
        for (name,dims,value) in [
            ("E_inf",Pressure::DIMS,1e8),("E0",Pressure::DIMS,3e8),
            ("tau0",Time::DIMS,1.0),("scale",ONE,1.0),("shift",ONE,1.0),
            ("density",Density::DIMS,1000.0),("cp",Dims([2,0,-2,-1,0,0]),cp),
        ] { add(&mut claims,name,dims,vec![(280.0,value),(320.0,value)]); }
        add(&mut claims,"free",ONE,
            vec![(280.0,-20.0*alpha),(300.0,middle_offset),(320.0,20.0*alpha)]);
        let card=card(claims);
        let source=resolve(&card,&point(280.0),&point(320.0)).unwrap();
        let heat=caloric(&card,cp,&point(280.0),&point(320.0));
        assert_eq!(source.samples().len(),3);
        let old=source.caloric_point(&heat,300.0,0.002,1e-6,0.001);
        if !admissible {
            assert!(matches!(old,Err(SourceError::Feedback(_))));
            continue;
        }
        let old=old.unwrap();
        let step=old.advance_controlled(MechanicalControl::Stress {
            stress_pa:0.0,tolerance_pa:1e-6 },0.2,FeedbackConfig {
                max_iterations:10,temperature_tolerance_k:1e-10,
                energy_tolerance_j:1e-9,relaxation:1.0 },|input| {
            let target_h=15.0*cp;
            Ok::<_,Infallible>(UniformEnthalpyStep {
                state:input.curve.state_at_specific_enthalpy(target_h).unwrap(),
                external_heat_j:input.mass_kg*(target_h-input.initial.specific_enthalpy_j_kg())
                    -input.internal_heat_j,
                energy_residual_tolerance_j:1e-10,report:(),
            })
        },||false).unwrap();
        assert!((step.next.mechanical().total_strain()-(0.002+15.0*alpha)).abs()<1e-15);
        assert!(step.next.mechanical().stress_pa().abs()<1e-7);
        assert!(step.next.mechanical().branch_strains()[0].abs()<1e-15);
        assert!(step.next.mechanical_internal_energy_j().abs()<1e-12);
        assert!(step.reversible_heat_j.abs()<1e-12);
        assert!(step.combined_residual_j.abs()<1e-10);
        assert_eq!(step.next.caloric_curve().material_card_identity(),card.content_hash());
        assert!(core::ptr::eq(step.next.mechanical().law(),source.law()));
        assert_eq!(old.mechanical().total_strain(),0.002);
        assert_eq!(old.time_s(),0.0);
    }
}

fn static_quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn pressure_point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", static_quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("p", static_quantity(QuantityKind::Pressure), p)
        .unwrap()
}

fn surface_card(mut change: impl FnMut(&mut PropertyClaim)) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("E_inf", Pressure::DIMS, [1e8; 4]),
        ("E0", Pressure::DIMS, [3e8; 4]),
        ("tau0", Time::DIMS, [1.0; 4]),
        ("scale", ONE, [1.0; 4]),
        ("shift", ONE, [1.0, 1.0, 1.0, 0.1]),
        ("free", ONE, [0.0; 4]),
        ("density", Density::DIMS, [1000.0; 4]),
        ("cp", Dims([2, 0, -2, -1, 0, 0]), [5.0; 4]),
    ] {
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "p".into()],
                axis_dims: [Temperature::DIMS, Pressure::DIMS],
                samples: [[280.0, 1e5], [320.0, 1e5], [280.0, 3e5], [320.0, 3e5]]
                    .into_iter()
                    .zip(values)
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    static_quantity(QuantityKind::AbsoluteTemperature),
                    280.0,
                    320.0,
                )
                .with_quantity("p", static_quantity(QuantityKind::Pressure), 1e5, 3e5),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic pressure-conditioned Maxwell surface".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        };
        change(&mut claim);
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

#[test]
fn fixed_pressure_surface_edges_drive_the_original_maxwell_history() {
    let card = surface_card(|_| {});
    let lower = pressure_point(280.0, 2e5);
    let upper = pressure_point(320.0, 2e5);
    let source = resolve(&card, &lower, &upper).unwrap();
    let temperatures: Vec<_> = source
        .samples()
        .iter()
        .map(|sample| {
            sample
                .query_point()
                .iter()
                .find(|(axis, _)| axis == "T")
                .unwrap()
                .1
        })
        .collect();
    // The diagonal cell edge at 300 K is essential: the shift is constant
    // before it and changes slope afterward, despite having no sample here.
    assert_eq!(temperatures, vec![280.0, 300.0, 320.0]);
    let direct = ThermoMaxwell::new(
        1e8,
        vec![(3e8, 1.0)],
        vec![
            ThermalKnot {
                temperature_k: 280.0,
                modulus_scale: 1.0,
                relaxation_scale: 1.0,
                free_strain: 0.0,
            },
            ThermalKnot {
                temperature_k: 300.0,
                modulus_scale: 1.0,
                relaxation_scale: 1.0,
                free_strain: 0.0,
            },
            ThermalKnot {
                temperature_k: 320.0,
                modulus_scale: 1.0,
                relaxation_scale: 0.55,
                free_strain: 0.0,
            },
        ],
        limits(),
    )
    .unwrap();
    let mut a = source.law().birth(280.0, 0.0).unwrap();
    let mut b = direct.birth(280.0, 0.0).unwrap();
    for (temperature, strain) in [(290.0, 0.01), (307.0, 0.02), (315.0, 0.01), (298.0, 0.0)] {
        let actual = a.advance(temperature, strain, 0.2, 1e-6).unwrap();
        let expected = b.advance(temperature, strain, 0.2, 1e-6).unwrap();
        assert!((actual.stress_pa - expected.stress_pa).abs() < 1e-6);
        assert!(
            (actual.point.stored_energy_j_m3() - expected.point.stored_energy_j_m3()).abs() < 1e-6
        );
        a = actual.point;
        b = expected.point;
    }
    for sample in source.samples() {
        assert_eq!(sample.card_identity(), card.content_hash());
        assert!(sample.query_point().contains(&("p".into(), 2e5)));
        for key in source.keys() {
            let receipt = &sample.property(key.name()).unwrap().answer().receipt;
            card.claims().verify_receipt(receipt).unwrap();
            assert_eq!(
                receipt.selected,
                source.samples()[0]
                    .property(key.name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    let thermal = caloric(&card, 5.0, &lower, &upper);
    assert!(
        source
            .caloric_point(&thermal, 300.0, 0.0, 1e-6, 0.001)
            .is_ok()
    );
}

#[test]
fn surface_paths_refuse_holes_discrete_spans_and_hidden_branch_changes() {
    let gap = surface_card(|claim| {
        if claim.key.name() != "shift" {
            return;
        }
        let right = f64::from_bits(300_f64.to_bits() + 1);
        claim.value = PropertyValue::Surface2D {
            axes: ["T".into(), "p".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples: [
                [280.0, 1e5],
                [300.0, 1e5],
                [280.0, 3e5],
                [300.0, 3e5],
                [right, 1e5],
                [320.0, 1e5],
                [right, 3e5],
                [320.0, 3e5],
            ]
            .map(|at| (at, 1.0))
            .to_vec(),
            triangles: vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
            dims: ONE,
        };
    });
    assert!(
        resolve(
            &gap,
            &pressure_point(280.0, 2e5),
            &pressure_point(320.0, 2e5)
        )
        .is_err()
    );
    let discrete = surface_card(|claim| {
        if claim.key.name() == "shift" {
            claim.interpolation = InterpolationPolicy::TabulatedOnly;
        }
    });
    // Exact endpoint samples cannot authorize the unsupported interval.
    assert!(
        resolve(
            &discrete,
            &pressure_point(280.0, 1e5),
            &pressure_point(320.0, 1e5)
        )
        .is_err()
    );
    let hidden = surface_card(|claim| {
        if claim.key.name() == "E0" {
            let PropertyValue::Surface2D { samples, .. } = &mut claim.value else {
                unreachable!()
            };
            for ((_, value), replacement) in samples.iter_mut().zip([3e8, 2e8, 2e8, 3e8]) {
                *value = replacement;
            }
        }
    });
    // At p=2e5 the two endpoints agree, but the crossed diagonal has a
    // different reference modulus. A single shared multiplier cannot model it.
    assert!(matches!(
        resolve(
            &hidden,
            &pressure_point(280.0, 2e5),
            &pressure_point(320.0, 2e5)
        ),
        Err(SourceError::Invalid(_))
    ));
}

#[test]
fn surface_context_preserves_semantic_kinds_and_temperature_axis_dimensions() {
    let source = surface_card(|_| {});
    let stress_point = |t| {
        QueryPoint::new()
            .with_quantity("T", static_quantity(QuantityKind::AbsoluteTemperature), t)
            .unwrap()
            .with_quantity("p", static_quantity(QuantityKind::Stress), 2e5)
            .unwrap()
    };
    assert!(resolve(&source, &stress_point(280.0), &stress_point(320.0)).is_err());
    let wrong = surface_card(|claim| {
        claim.validity = ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0)
            .with_quantity("p", static_quantity(QuantityKind::Pressure), 1e5, 3e5);
        if claim.key.name() == "shift" {
            let PropertyValue::Surface2D { axis_dims, .. } = &mut claim.value else {
                unreachable!()
            };
            axis_dims[0] = Pressure::DIMS;
        }
    });
    let untyped = |t| {
        point(t)
            .with_quantity("p", static_quantity(QuantityKind::Pressure), 2e5)
            .unwrap()
    };
    // Untyped legacy T can be queried by matdb, but it cannot make a
    // pressure-dimensional source axis into an absolute-temperature law.
    wrong
        .claims()
        .query_typed(
            &PropertyKey::new("shift", ONE),
            &untyped(300.0),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!(matches!(
        resolve(&wrong, &untyped(280.0), &untyped(320.0)),
        Err(SourceError::Invalid(_))
    ));
}

#[test]
fn mechanical_and_second_temperature_conditions_cannot_be_frozen_by_source_slicing() {
    for (quantity, lo, hi, at) in [
        (
            static_quantity(QuantityKind::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Normal,
            }),
            0.0,
            0.02,
            0.01,
        ),
        (static_quantity(QuantityKind::Stress), 1e5, 3e5, 2e5),
        (
            QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
            0.0,
            100.0,
            50.0,
        ),
        (
            static_quantity(QuantityKind::AbsoluteTemperature),
            280.0,
            320.0,
            300.0,
        ),
    ] {
        for value_depends in [false, true] {
            let card = surface_card(|claim| {
                if claim.key.name() != "shift" {
                    return;
                }
                claim.validity =
                    claim
                        .validity
                        .clone()
                        .with_quantity("condition", quantity, lo, hi);
                if value_depends {
                    let PropertyValue::Surface2D {
                        axes,
                        axis_dims,
                        samples,
                        ..
                    } = &mut claim.value
                    else {
                        unreachable!()
                    };
                    axes[1] = "condition".into();
                    axis_dims[1] = quantity.dims();
                    for (at, _) in samples {
                        at[1] = if at[1] == 1e5 { lo } else { hi };
                    }
                } else {
                    claim.value = PropertyValue::Scalar {
                        value: 1.0,
                        dims: ONE,
                    };
                    claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
                }
            });
            let context = |t| {
                pressure_point(t, 2e5)
                    .with_quantity("condition", quantity, at)
                    .unwrap()
            };
            assert!(
                resolve(&card, &context(280.0), &context(320.0)).is_err(),
                "cannot freeze {quantity:?}, including applicability-only bounds"
            );
        }
    }
    let legacy = surface_card(|claim| {
        if claim.key.name() == "shift" {
            claim.validity = claim.validity.clone().with("strain_rate", 0.0, 100.0);
        }
    });
    let context = |t| pressure_point(t, 2e5).with("strain_rate", 50.0).unwrap();
    assert!(resolve(&legacy, &context(280.0), &context(320.0)).is_err());
}

#[test]
fn frequency_applicability_retains_static_coverage_without_freezing_frequency_dependence() {
    let frequency = static_quantity(QuantityKind::Frequency(FrequencyConvention::Angular));
    for typed in [false, true] {
        for low in [0.0, 1.0] {
            let card = surface_card(|claim| {
                claim.validity = if typed {
                    claim
                        .validity
                        .clone()
                        .with_quantity("omega", frequency, low, 100.0)
                } else {
                    claim.validity.clone().with("omega", low, 100.0)
                };
            });
            let context = |t| {
                let point = pressure_point(t, 2e5);
                if typed {
                    point.with_quantity("omega", frequency, 50.0).unwrap()
                } else {
                    point.with("omega", 50.0).unwrap()
                }
            };
            assert_eq!(
                resolve(&card, &context(280.0), &context(320.0)).is_ok(),
                low == 0.0
            );
        }
    }
    let dependent = surface_card(|claim| {
        if claim.key.name() != "shift" {
            return;
        }
        claim.validity = claim
            .validity
            .clone()
            .with_quantity("omega", frequency, 0.0, 100.0);
        let PropertyValue::Surface2D {
            axes,
            axis_dims,
            samples,
            ..
        } = &mut claim.value
        else {
            unreachable!()
        };
        axes[1] = "omega".into();
        axis_dims[1] = fs_qty::Frequency::DIMS;
        for (at, _) in samples {
            at[1] = if at[1] == 1e5 { 0.0 } else { 100.0 };
        }
    });
    // Even zero is not an implicit authority to turn frequency-dependent data
    // into the time-domain spring/dashpot coefficients of this owner.
    for omega in [0.0, 50.0] {
        let context = |t| {
            pressure_point(t, 2e5)
                .with_quantity("omega", frequency, omega)
                .unwrap()
        };
        assert!(resolve(&dependent, &context(280.0), &context(320.0)).is_err());
    }
}
