//! Original hardening and independent kinetics drive native timed deformation.
use super::*;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::tabulated::rate::axial::UniaxialTabulatedOverstress;
use fs_material::plastic::tabulated::rate::source::{
    ResolvedTabulatedOverstressJ2, TabulatedOverstressKinetics, TabulatedOverstressSourceRequest,
    resolve_tabulated_overstress_j2, tabulated_overstress_requirements,
};

const ETA: f64 = 40.0;
const SCALE: f64 = 2.0;
const EXPONENT: f64 = 3.2;
const DT: f64 = 0.2;
const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-9,
    yield_pa: 1e-9,
};

fn step_request(kinetics: TabulatedOverstressKinetics) -> TabulatedOverstressSourceRequest {
    TabulatedOverstressSourceRequest {
        plasticity: TabulatedJ2SourceRequest {
            kinematic_hardening: TabulatedKinematicHardening::Prager,
            strain_limit: 0.05,
            maximum_points: KNOTS.len(),
        },
        kinetics,
        duration_s: DT,
    }
}
fn timed_properties(power: bool) -> Vec<PropertyClaim> {
    let mut result = properties();
    result.push(scalar(
        OVERSTRESS_VISCOSITY_PROPERTY,
        OVERSTRESS_VISCOSITY_DIMS,
        ETA,
    ));
    if power {
        result.push(scalar(
            OVERSTRESS_STRESS_SCALE_PROPERTY,
            Pressure::DIMS,
            SCALE,
        ));
        result.push(scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, EXPONENT));
    }
    result
}
fn resolved(
    card: &MaterialCard,
    kinetics: TabulatedOverstressKinetics,
) -> ResolvedTabulatedOverstressJ2 {
    resolve_tabulated_overstress_j2(
        card,
        &context(),
        &step_request(kinetics),
        MaterialPropertySelection::SingleClaimOnly,
        SINGLE,
    )
    .unwrap()
}
fn tensor(axial: f64) -> [f64; 6] {
    [axial, -0.5 * axial, -0.5 * axial, 0.0, 0.0, 0.0]
}

#[test]
fn source_rate_retains_density_and_curve_evidence_and_native_physical_energy() {
    let card = card(timed_properties(true));
    let mut requirements = tabulated_overstress_requirements(
        TabulatedKinematicHardening::Prager,
        TabulatedOverstressKinetics::PowerLaw,
    )
    .unwrap();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        &card,
        &context(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let curve = resolve_tabulated_hardening(
        &card,
        &context(),
        &request(YieldStressMeasure::Cauchy),
        SINGLE,
    )
    .unwrap();
    let plastic = ResolvedTabulatedJ2::try_from_material(
        bundle.clone(),
        curve.clone(),
        TabulatedKinematicHardening::Prager,
        0.05,
    )
    .unwrap();
    let law = ResolvedTabulatedOverstressJ2::try_from_source(
        plastic,
        TabulatedOverstressKinetics::PowerLaw,
        DT,
    )
    .unwrap();
    assert_eq!(law.material(), &bundle);
    assert_eq!(law.curve(), &curve);
    assert_eq!(
        law.material()
            .property(DENSITY_PROPERTY)
            .unwrap()
            .value_si(),
        7800.0
    );
    assert_eq!(law.physical_duration_s(), Some(DT));
    for property in law.material().properties() {
        card.claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let alpha = 0.009;
    let viscous = SCALE * (ETA * alpha / (DT * SCALE)).powf(1.0 / EXPONENT);
    let strain = tensor(((1500.0 + 160.0) * alpha + 5.5 + viscous) / 1500.0);
    let old = law.initial_state();
    let response = law.response(&strain, &old).unwrap();
    let increment = law
        .advance_with_energy(&[0.0; 6], &old, &strain, TOL)
        .unwrap();
    near(response.state.alpha, alpha, 2e-14);
    near(
        response.stress_pa[0],
        (2.0 / 3.0) * (5.5 + 160.0 * alpha + viscous),
        2e-11,
    );
    near(
        increment.energy.physical_dissipation_j_m3,
        (3.0 + viscous) * alpha,
        2e-12,
    );
    near(
        increment.energy.final_state.hardening_j_m3,
        curve.hardening().hardening_energy_j_m3(alpha).unwrap(),
        2e-13,
    );
    near(increment.energy.residual_j_m3, 0.0, 2e-12);
    let (state, stress, account) = law
        .native_energy_increment(&[0.0; 6], &old, &strain, TOL)
        .unwrap();
    assert_eq!(state, increment.state);
    assert_eq!(stress, increment.stress_pa);
    assert_eq!(account, increment.energy);
    assert_eq!(old.alpha, 0.0);
}

#[test]
fn duration_rebinding_shares_calibration_and_relaxes_original_overstressed_history() {
    let card = card(timed_properties(true));
    let law = resolved(&card, TabulatedOverstressKinetics::PowerLaw);
    let strain = tensor(0.019);
    let first = law
        .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, TOL)
        .unwrap();
    assert!(first.state.alpha > 0.0);
    assert!(
        law.source()
            .law()
            .yield_function_with_state(&first.stress_pa, &first.state)
            > 0.0
    );
    let quick = law.for_duration(0.01).unwrap();
    let slow = law.for_duration(1.0).unwrap();
    for next in [&quick, &slow] {
        assert!(std::ptr::eq(law.source(), next.source()));
        assert!(std::ptr::eq(
            law.law().backbone().hardening().points(),
            next.law().backbone().hardening().points()
        ));
        assert_eq!(law.curve(), next.curve());
        next.validate_history(&strain, &first.state).unwrap();
    }
    let short = quick
        .advance_with_energy(&strain, &first.state, &strain, TOL)
        .unwrap();
    let long = slow
        .advance_with_energy(&strain, &first.state, &strain, TOL)
        .unwrap();
    assert!(short.state.alpha > first.state.alpha);
    assert!(long.state.alpha > short.state.alpha);
    assert!(long.stress_pa[0] < short.stress_pa[0]);
    assert!(short.energy.physical_dissipation_j_m3 > 0.0);
    assert!(long.energy.physical_dissipation_j_m3 > 0.0);
    near(short.energy.mechanical_work_j_m3, 0.0, 1e-15);
    near(long.energy.mechanical_work_j_m3, 0.0, 1e-15);
    near(short.energy.residual_j_m3, 0.0, TOL.energy_j_m3);
    near(long.energy.residual_j_m3, 0.0, TOL.energy_j_m3);
    assert_eq!(law.duration_s(), DT);
    assert_eq!(quick.duration_s(), 0.01);
    assert_eq!(slow.duration_s(), 1.0);
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(law.for_duration(invalid).is_err());
    }
    assert_eq!(
        first.state,
        law.response(&strain, &law.initial_state()).unwrap().state
    );
}

#[test]
fn linear_kinetics_require_only_viscosity_and_never_infer_a_power_law() {
    let card = card(timed_properties(false));
    let linear = resolved(&card, TabulatedOverstressKinetics::Linear);
    assert_eq!(linear.law().stress_scale_pa(), None);
    assert_eq!(linear.law().exponent(), 1.0);
    assert!(
        resolve_tabulated_overstress_j2(
            &card,
            &context(),
            &step_request(TabulatedOverstressKinetics::PowerLaw),
            MaterialPropertySelection::SingleClaimOnly,
            SINGLE
        )
        .is_err()
    );
    let alpha = 0.003;
    let strain = tensor(((1500.0 + 160.0) * alpha + 4.0 + ETA * alpha / DT) / 1500.0);
    let response = linear.response(&strain, &linear.initial_state()).unwrap();
    near(response.state.alpha, alpha, 2e-15);
    let longer = linear.for_duration(2.0).unwrap();
    assert_eq!(longer.kinetics(), TabulatedOverstressKinetics::Linear);
    assert_eq!(longer.law().stress_scale_pa(), None);
    assert!(std::ptr::eq(linear.source(), longer.source()));
}

#[test]
fn source_kinetic_units_bounds_and_complete_bundle_are_required() {
    for variant in 0..4 {
        let mut props = timed_properties(true);
        let name = if variant == 2 {
            OVERSTRESS_EXPONENT_PROPERTY
        } else {
            OVERSTRESS_VISCOSITY_PROPERTY
        };
        let property = props.iter_mut().find(|p| p.key.name() == name).unwrap();
        match variant {
            0 => {
                property.key = PropertyKey::new(name, Pressure::DIMS);
                property.value = PropertyValue::Scalar {
                    value: ETA,
                    dims: Pressure::DIMS,
                };
            }
            1 => {
                property.value = PropertyValue::Scalar {
                    value: 0.0,
                    dims: OVERSTRESS_VISCOSITY_DIMS,
                };
            }
            2 => {
                property.value = PropertyValue::Scalar {
                    value: 0.5,
                    dims: Dims::NONE,
                };
            }
            _ => {
                property.key =
                    PropertyKey::new("fluid_dynamic_viscosity", OVERSTRESS_VISCOSITY_DIMS);
            }
        }
        assert!(
            resolve_tabulated_overstress_j2(
                &card(props),
                &context(),
                &step_request(TabulatedOverstressKinetics::PowerLaw),
                MaterialPropertySelection::SingleClaimOnly,
                SINGLE
            )
            .is_err()
        );
    }
    let card = card(timed_properties(true));
    let plastic = resolve_tabulated_j2(
        &card,
        &context(),
        &step_request(TabulatedOverstressKinetics::Linear).plasticity,
        MaterialPropertySelection::SingleClaimOnly,
        SINGLE,
    )
    .unwrap();
    assert!(
        ResolvedTabulatedOverstressJ2::try_from_source(
            plastic,
            TabulatedOverstressKinetics::Linear,
            DT
        )
        .is_err(),
        "adopting a bundle must not query omitted kinetic roles implicitly"
    );
}

#[test]
fn sourced_axial_flow_preserves_original_curve_and_kinetics_across_physical_steps() {
    let card = card(timed_properties(true));
    let source = resolved(&card, TabulatedOverstressKinetics::PowerLaw);
    let law = UniaxialTabulatedOverstress::from_source(source.clone()).unwrap();
    let old = law.initial_state();
    let alpha = 0.009;
    let viscous = SCALE * (ETA * alpha / (DT * SCALE)).powf(1.0 / EXPONENT);
    let strain = ((1200.0 + 160.0) * alpha + 5.5 + viscous) / 1200.0;
    let next = law.response(strain, &old).unwrap();
    near(next.state.equivalent_plastic_strain(), alpha, 2e-14);
    near(next.state.plastic_strain(), alpha, 2e-14);
    near(next.stress_pa, 5.5 + 160.0 * alpha + viscous, 2e-11);
    assert_eq!(law.material(), Some(source.material()));
    assert!(std::ptr::eq(
        law.source().unwrap().source(),
        source.source()
    ));
    let rebound = law.for_duration(0.7).unwrap();
    assert!(std::ptr::eq(
        rebound.source().unwrap().source(),
        source.source()
    ));
    assert!(std::ptr::eq(
        rebound.backbone().backbone().hardening().points(),
        law.backbone().backbone().hardening().points()
    ));
    assert_eq!(rebound.source().unwrap().duration_s(), 0.7);
    assert_eq!(law.source().unwrap().duration_s(), DT);
    let relaxed = rebound.response(strain, &next.state).unwrap();
    assert!(relaxed.state.equivalent_plastic_strain() > next.state.equivalent_plastic_strain());
    assert!(relaxed.stress_pa < next.stress_pa);
    let increment = rebound
        .advance_with_energy(strain, &next.state, strain, TOL)
        .unwrap();
    near(increment.energy.mechanical_work_j_m3, 0.0, 2e-12);
    near(increment.energy.residual_j_m3, 0.0, 2e-12);
    for stress in &increment.stress_pa[1..] {
        near(*stress, 0.0, TOL.yield_pa);
    }
    assert_eq!(old.plastic_strain(), 0.0);
    assert_eq!(old.equivalent_plastic_strain(), 0.0);
}
