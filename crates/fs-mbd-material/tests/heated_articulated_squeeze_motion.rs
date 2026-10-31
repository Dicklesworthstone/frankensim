//! Native Reynolds pressure -> articulated effort -> finite heat -> original liquid.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, LiquidPropertyRequirements,
};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_PROPERTY, resolve_material_state_point,
};
use fs_mbd::articulated::JointModel;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::articulated_motion::squeeze::{SqueezeFilm, SqueezeLimits};
use fs_mbd_material::thermal::heating::articulated_motion::{Controls, HeatedArticulatedMotion};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;
use fs_tribo::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm};

const AREA: f64 = 1e-4;
const SPEED: f64 = 0.04;
const SATURATION: &str = "saturation_pressure";

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn lubricant(mu: [f64; 2], saturation: [f64; 2], pressure: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (DENSITY_PROPERTY, Density::DIMS, [2., 2.]),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, [1., 1.]),
        (THERMAL_CONDUCTIVITY_PROPERTY, Dims([1, 1, -3, -1, 0, 0]), [0.5, 0.5]),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu),
        (SATURATION, Pressure::DIMS, saturation),
    ] {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
        } else { QuantitySpec::dimensional(dims) };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name, quantity),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300., values[0]), (400., values[1])], dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 300., 400.)
                .with_quantity("pressure", pressure_quantity(), pressure[0], pressure[1]),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance {
                source: "synthetic squeeze fixture, not measured oil".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic lubricant".into(), phase: "liquid".into(),
        process: "fixture".into(), revision: 0,
    }, claims, vec![]).unwrap()
}

fn mechanism(solid: &MaterialCard, two: bool) -> ResolvedArticulated {
    let p = point(300., 1e5);
    let regions = [RegionBinding {
        key: "a", material: solid, point: &p,
        geometry: VolumeMoments::cuboid(Vec3::new(1., 1., 1.)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    let slide = JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap();
    let mut links = vec![
        MaterialLink { name: "mount", parent: None, parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED, regions: &regions },
        MaterialLink { name: "piston", parent: Some(0), parent_from_child_zero: Se3::identity(),
            joint: slide, regions: &regions },
    ];
    if two {
        links.push(MaterialLink { name: "second", parent: Some(0),
            parent_from_child_zero: Se3::identity(), joint: slide, regions: &regions });
    }
    resolve_articulated(&links, 3, || false).unwrap()
}

fn moving<'a>(reference: &'a ResolvedArticulated, solid: &'a MaterialCard,
    heat: &'a ResolvedHeatCapacityCurve, speed: Vec<f64>, dt: f64,
) -> HeatedArticulatedMotion<'a> {
    let bindings = reference.model().links().iter().enumerate().map(|(i, link)| LinkHeating {
        link: link.name(), regions: vec![region(solid, heat, "a", i as u64, 0.)],
    }).collect();
    let heated = HeatedArticulated::new(reference, bindings, vec![], budget(), config(), || Ok(())).unwrap();
    let n = speed.len();
    HeatedArticulatedMotion::new(heated, Se3::identity(), vec![0.; n], speed,
        vec![0.; n], Vec3::default(), Controls {
            maximum_step_s: dt, relative_state_tolerance: 1e-11,
            absolute_position_tolerance: vec![1e-12; n],
            absolute_velocity_tolerance: vec![1e-12; n],
            maximum_rk_attempts_per_substep: 256, relative_momentum_tolerance: 1e-11,
            absolute_momentum_tolerance: vec![1e-12; n], maximum_parameter_work_j: 0.,
            maximum_numerical_energy_change_j: 1e-9,
        }, || Ok(())).unwrap()
}

fn film(oil: &MaterialCard, closure: Vec<f64>) -> SqueezeFilm<'_> {
    let n = closure.len();
    let point = QueryPoint::new().with("T", 300.).unwrap()
        .with_quantity("pressure", pressure_quantity(), 1e5).unwrap();
    let mut requirements = LiquidPropertyRequirements::canonical(LiquidViscosityBasis::Dynamic)
        .unwrap().properties().to_vec();
    requirements.push(ScalarPropertyRequirement::try_new(
        SATURATION, Pressure::DIMS, ScalarAdmissibility::NonNegative).unwrap());
    let material = resolve_material_state_point(oil, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    SqueezeFilm {
        heat_cell: CellId(1),
        // A chamber empties through a fixed thin slit. Its volume changes, but
        // the outlet stays fixed, giving c = 12 mu L A²/(w h³) = 120 mu.
        geometry: ResistiveFilm::new(vec![FilmCell {
            area_m2: AREA, gap: GapPort { reference_m: 0.1, closure },
        }], vec![FilmChannel {
            from: 0, to: None, width_m: 0.1, length_m: 0.1,
            gap: GapPort { reference_m: 0.001, closure: vec![0.; n] },
        }], n, 123., FilmLimits {
            minimum_cell_gap_m: 0.01, maximum_gap_m: 0.2, maximum_pressure_pa: 1000.,
        }).unwrap(),
        liquid: LiquidTemperatureLaw::try_from_material(oil, material, "T", LiquidViscosityBasis::Dynamic).unwrap(),
        pressure_axis: "pressure", saturation_pressure_property: SATURATION,
        limits: SqueezeLimits {
            maximum_relative_gauge_pressure: 0.01,
            maximum_relative_viscosity_change: 0.5, absolute_loss_tolerance_j: 1e-12,
        },
    }
}

#[test]
fn source_viscosity_swap_changes_signed_motion_pressure_and_heat_once() {
    let solid = card(10., 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid, false);
    let mut final_speeds = Vec::new();
    for mu in [0.01, 0.02] {
        let oil = lubricant([mu, mu], [2000., 2000.], [99000., 101000.]);
        for sign in [-1., 1.] {
            let old = moving(&reference, &solid, &heat, vec![sign * SPEED], 0.0625)
                .with_squeeze_films(vec![film(&oil, vec![1.])], || Ok(())).unwrap();
            let c = 12. * mu * 0.1 * AREA.powi(2) / (0.1 * 0.001_f64.powi(3));
            let decay = (-c / 2.).exp(); // Native density × unit volume gives m=2 kg.
            let (next, audit) = old.advance_to(1., config(), || Ok(())).unwrap();
            let state = next.squeeze_films().next().unwrap();
            let loss = SPEED.powi(2) * (1. - decay * decay);
            close(next.joint_velocity()[0], sign * SPEED * decay);
            close(next.joint_position()[0], sign * 2. * SPEED * (1. - decay) / c);
            close(state.generalized_resistance[0], c * next.joint_velocity()[0]);
            close(state.gauge_pressure_pa[0], state.generalized_resistance[0] / AREA);
            close(state.minimum_gap_m, 0.1 - next.joint_position()[0]);
            close(state.dissipated_power_w, c * next.joint_velocity()[0].powi(2));
            close(state.dissipated_heat_j, loss);
            close(audit.viscous_heat_j, loss);
            close(audit.thermal.external_heat_j, loss);
            close(next.thermal().temperature_k(CellId(1)).unwrap(), 300. + loss / 20.);
            close(next.thermal().temperature_k(CellId(0)).unwrap(), 300.);
            close(audit.mechanical_after.kinetic_j + loss, audit.mechanical_before.kinetic_j);
            close(audit.combined_residual_j, 0.);
            close(state.ambient_pressure_pa, 1e5);
            close(state.saturation_pressure_pa, 2000.);
            assert_eq!(state.liquid.material().properties().len(), 5);
            assert_eq!(old.thermal().time_s(), 0.);
            if sign > 0. { final_speeds.push(next.joint_velocity()[0]); }
        }
    }
    assert!(final_speeds[1] < final_speeds[0]);
}

#[test]
fn one_pressure_graph_couples_sibling_links_and_conserves_total_momentum() {
    let solid = card(10., 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid, true);
    let oil = lubricant([0.01, 0.01], [2000., 2000.], [99000., 101000.]);
    let old = moving(&reference, &solid, &heat, vec![0.04, -0.02], 0.0625)
        .with_squeeze_films(vec![film(&oil, vec![1., -1.])], || Ok(())).unwrap();
    let (next, audit) = old.advance_to(1., config(), || Ok(())).unwrap();
    // Each sibling has mass 2: the relative mode has reduced mass 1 kg.
    let relative = 0.06 * (-1.2_f64).exp();
    close(next.joint_velocity()[0], 0.01 + 0.5 * relative);
    close(next.joint_velocity()[1], 0.01 - 0.5 * relative);
    close(next.joint_momentum().unwrap().iter().sum(), old.joint_momentum().unwrap().iter().sum());
    let state = next.squeeze_films().next().unwrap();
    close(state.generalized_resistance[0], 1.2 * relative);
    close(state.generalized_resistance[1], -1.2 * relative);
    close(state.dissipated_power_w, 1.2 * relative.powi(2));
    close(audit.viscous_heat_j, 0.5 * (0.06_f64.powi(2) - relative.powi(2)));
    close(next.thermal().temperature_k(CellId(2)).unwrap(), 300.);
    close(next.thermal().temperature_k(CellId(1)).unwrap(), 300. + audit.viscous_heat_j / 20.);
    close(audit.combined_residual_j, 0.);
}

#[test]
fn original_source_self_heating_feedback_converges_and_restarts_without_duplicate_heat() {
    let solid = card(1e-5, 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid, false);
    let oil = lubricant([0.01, 0.001], [2000., 3000.], [99000., 101000.]);
    let mut errors = Vec::new();
    for dt in [0.0625, 0.03125, 0.015625] {
        let old = moving(&reference, &solid, &heat, vec![SPEED], dt)
            .with_squeeze_films(vec![film(&oil, vec![1.])], || Ok(())).unwrap();
        // E+C*T is constant, C=2e-5 J/K. Thus mu/mu0=a+b*v² and
        // d(v²)/dt=-2*(c0/m)*(a+b*v²)*v² has this Bernoulli solution.
        let b = 0.009 / 2e-5;
        let a = 1. - b * SPEED.powi(2);
        let e = (-1.2 * a).exp();
        let expected = (a * SPEED.powi(2) * e / (a + b * SPEED.powi(2) * (1. - e))).sqrt();
        let (next, audit) = old.advance_to(1., config(), || Ok(())).unwrap();
        let temperature = next.thermal().temperature_k(CellId(1)).unwrap();
        let state = next.squeeze_films().next().unwrap();
        let mu = 0.01 - 0.00009 * (temperature - 300.);
        close(state.liquid.dynamic_viscosity_pa_s(), mu);
        close(state.generalized_resistance[0], 120. * mu * next.joint_velocity()[0]);
        close(state.saturation_pressure_pa, 2000. + 10. * (temperature - 300.));
        close(state.liquid.material().query().axes()["T"], temperature);
        close(audit.viscous_heat_j, 2e-5 * (temperature - 300.));
        close(audit.combined_residual_j, 0.);
        let original = old.squeeze_films().next().unwrap();
        for (before, after) in original.liquid.material().properties().iter()
            .zip(state.liquid.material().properties()) {
            assert_eq!(before.answer().receipt.selected, after.answer().receipt.selected);
        }
        assert_ne!(original.liquid.material().identity(), state.liquid.material().identity());
        assert!(next.joint_velocity()[0] > SPEED * (-0.6_f64).exp());
        errors.push((next.joint_velocity()[0] - expected).abs());
        let (half, first) = old.advance_to(0.5, config(), || Ok(())).unwrap();
        let (split, second) = half.advance_to(1., config(), || Ok(())).unwrap();
        assert_eq!(next.joint_position(), split.joint_position());
        assert_eq!(next.joint_velocity(), split.joint_velocity());
        assert_eq!(next.thermal().cells(), split.thermal().cells());
        assert_eq!(next.squeeze_films().next(), split.squeeze_films().next());
        close(audit.viscous_heat_j, first.viscous_heat_j + second.viscous_heat_j);
        close(state.dissipated_heat_j, audit.viscous_heat_j);
    }
    assert!(errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1], "{errors:?}");
}

#[test]
fn pressure_support_saturation_and_late_cancel_refuse_without_publishing_a_candidate() {
    let solid = card(1e-5, 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid, false);
    let moving = moving(&reference, &solid, &heat, vec![SPEED], 0.0625);
    let narrow = lubricant([0.01, 0.01], [2000., 2000.], [99900., 100100.]);
    assert!(moving.with_squeeze_films(vec![film(&narrow, vec![1.])], || Ok(())).is_err());
    // Suction remains positive absolute pressure, but crosses the source's
    // saturation pressure. It must refuse instead of clamping negative gauge.
    let boiling = lubricant([0.01, 0.01], [99900., 99900.], [99000., 101000.]);
    assert!(moving.with_squeeze_films(vec![film(&boiling, vec![-1.])], || Ok(())).is_err());
    let rising = lubricant([0.01, 0.01], [99000., 105000.], [99000., 101000.]);
    let old = moving.with_squeeze_films(vec![film(&rising, vec![1.])], || Ok(())).unwrap();
    let snapshot = old.squeeze_films().next().unwrap().clone();
    assert!(old.advance_to(1., config(), || Ok(())).is_err());
    assert_eq!(old.squeeze_films().next().unwrap(), &snapshot);
    assert_eq!(old.joint_position(), &[0.]);
    assert_eq!(old.joint_velocity(), &[SPEED]);
    assert_eq!(old.thermal().time_s(), 0.);
    close(old.thermal().temperature_k(CellId(1)).unwrap(), 300.);

    let oil = lubricant([0.01, 0.001], [2000., 3000.], [99000., 101000.]);
    let old = moving.with_squeeze_films(vec![film(&oil, vec![1.])], || Ok(())).unwrap();
    let snapshot = old.squeeze_films().next().unwrap().clone();
    let mut calls = 0;
    old.advance_to(1., config(), || { calls += 1; Ok(()) }).unwrap();
    let mut count = 0;
    assert!(old.advance_to(1., config(), || {
        count += 1;
        if count == calls - 1 { Err(DepositionError::Cancelled) } else { Ok(()) }
    }).is_err());
    assert_eq!(old.squeeze_films().next().unwrap(), &snapshot);
    assert_eq!(old.joint_position(), &[0.]);
    assert_eq!(old.joint_velocity(), &[SPEED]);
    assert_eq!(old.thermal().time_s(), 0.);
    close(old.thermal().temperature_k(CellId(1)).unwrap(), 300.);
}
