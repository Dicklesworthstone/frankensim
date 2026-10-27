//! Actual joint loss -> local enthalpy -> source viscosity -> subsequent motion.
#[path = "heating/common.rs"]
mod common;
use common::*;
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_mbd::articulated::JointModel;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::articulated_motion::viscous::{
    JournalBearing, JournalGeometry, JournalLimits,
};
use fs_mbd_material::thermal::heating::articulated_motion::{
    Controls, Error, HeatedArticulatedMotion,
};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;

fn lubricant(mu: [f64; 2], upper: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (DENSITY_PROPERTY, Density::DIMS, [2., 2.]),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [1., 1.],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            Dims([1, 1, -3, -1, 0, 0]),
            [0.5, 0.5],
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu),
    ] {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300., values[0]), (upper, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300., upper),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic journal fixture, not measured oil".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic lubricant".into(),
            phase: "liquid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn mechanism(card: &MaterialCard, two: bool) -> ResolvedArticulated {
    let p = point(300., 1e5);
    let regions = [RegionBinding {
        key: "a",
        material: card,
        point: &p,
        geometry: VolumeMoments::cuboid(Vec3::new(1., 1., 1.)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    let hinge = JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap();
    let mut links = vec![
        MaterialLink {
            name: "mount",
            parent: None,
            parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED,
            regions: &regions,
        },
        MaterialLink {
            name: "rotor",
            parent: Some(0),
            parent_from_child_zero: Se3::identity(),
            joint: hinge,
            regions: &regions,
        },
    ];
    if two {
        links.push(MaterialLink {
            name: "second",
            parent: Some(1),
            parent_from_child_zero: Se3::identity(),
            joint: hinge,
            regions: &regions,
        });
    }
    resolve_articulated(&links, 3, || false).unwrap()
}

fn moving<'a>(
    reference: &'a ResolvedArticulated,
    solid: &'a MaterialCard,
    heat: &'a ResolvedHeatCapacityCurve,
    speed: Vec<f64>,
    dt: f64,
) -> HeatedArticulatedMotion<'a> {
    let bindings = reference
        .model()
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| LinkHeating {
            link: link.name(),
            regions: vec![region(solid, heat, "a", i as u64, 0.)],
        })
        .collect();
    let heated =
        HeatedArticulated::new(reference, bindings, vec![], budget(), config(), || Ok(())).unwrap();
    let n = speed.len();
    HeatedArticulatedMotion::new(
        heated,
        Se3::identity(),
        vec![0.; n],
        speed,
        vec![0.; n],
        Vec3::default(),
        Controls {
            maximum_step_s: dt,
            relative_state_tolerance: 1e-11,
            absolute_position_tolerance: vec![1e-12; n],
            absolute_velocity_tolerance: vec![1e-12; n],
            maximum_rk_attempts_per_substep: 256,
            relative_momentum_tolerance: 1e-11,
            absolute_momentum_tolerance: vec![1e-12; n],
            maximum_parameter_work_j: 0.,
            maximum_numerical_energy_change_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
}

fn bearing(oil: &MaterialCard, link: usize, cell: u64) -> JournalBearing<'_> {
    JournalBearing {
        link,
        heat_cell: CellId(cell),
        geometry: JournalGeometry {
            inner_radius_m: 0.1,
            outer_radius_m: 0.11,
            length_m: 0.2,
        },
        limits: JournalLimits {
            maximum_gap_reynolds: 1.,
            maximum_relative_viscosity_change: 0.5,
            absolute_loss_tolerance_j: 1e-12,
        },
        liquid: LiquidTemperatureLaw::try_new(
            oil,
            &QueryPoint::new().with("T", 300.).unwrap(),
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    }
}

#[test]
fn source_viscosity_swap_sets_native_spin_down_and_exact_local_heat() {
    let solid = card(10., 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid, false);
    let mut speeds = Vec::new();
    for mu in [1., 2.] {
        let oil = lubricant([mu, mu], 500.);
        for sign in [-1., 1.] {
            let old = moving(&reference, &solid, &heat, vec![sign * 3.], 0.0625)
                .with_journal_bearings(vec![bearing(&oil, 1, 1)], || Ok(()))
                .unwrap();
            let resistance = old.journal_bearings().next().unwrap().resistance_nm_s;
            let decay = (-3. * resistance).exp(); // Unit cube: Izz=1/3 kg m².
            let (next, audit) = old.advance_to(1., config(), || Ok(())).unwrap();
            close(next.joint_velocity()[0], sign * 3. * decay);
            close(next.joint_position()[0], sign * (1. - decay) / resistance);
            let loss = 1.5 * (1. - decay * decay);
            close(audit.viscous_heat_j, loss);
            close(audit.thermal.external_heat_j, loss);
            close(
                next.thermal().temperature_k(CellId(1)).unwrap(),
                300. + loss / 20.,
            );
            close(next.thermal().temperature_k(CellId(0)).unwrap(), 300.);
            close(
                audit.mechanical_after.kinetic_j + loss,
                audit.mechanical_before.kinetic_j,
            );
            close(audit.combined_residual_j, 0.);
            let state = next.journal_bearings().next().unwrap();
            close(
                state.resisting_torque_nm,
                -resistance * next.joint_velocity()[0],
            );
            close(
                state.dissipated_power_w,
                resistance * next.joint_velocity()[0].powi(2),
            );
            close(
                state.liquid.material().query().axes()["T"],
                next.thermal().temperature_k(CellId(1)).unwrap(),
            );
            assert_eq!(old.thermal().time_s(), 0.);
            if sign > 0. {
                speeds.push(next.joint_velocity()[0]);
            }
        }
    }
    assert!(speeds[1] < speeds[0]);
}

#[test]
fn relative_joint_shear_decelerates_both_native_articulated_links() {
    let solid = card(10., 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let oil = lubricant([1., 1.], 500.);
    let reference = mechanism(&solid, true);
    let old = moving(&reference, &solid, &heat, vec![3., -1.], 0.0625)
        .with_journal_bearings(vec![bearing(&oil, 2, 2)], || Ok(()))
        .unwrap();
    let c = old.journal_bearings().next().unwrap().resistance_nm_s;
    let (next, audit) = old.advance_to(1., config(), || Ok(())).unwrap();
    let decay = (-6. * c).exp();
    // M=[[2/3,1/3],[1/3,1/3]]. A torque only on relative joint 2
    // accelerates BOTH coordinates; absolute system angular momentum is fixed.
    close(next.joint_velocity()[1], -decay);
    close(next.joint_velocity()[0], 2.5 + 0.5 * decay);
    close(
        next.joint_momentum().unwrap()[0],
        old.joint_momentum().unwrap()[0],
    );
    close(audit.viscous_heat_j, (1. - decay * decay) / 12.);
    close(next.thermal().temperature_k(CellId(1)).unwrap(), 300.);
    close(
        next.thermal().temperature_k(CellId(2)).unwrap(),
        300. + audit.viscous_heat_j / 20.,
    );
    close(audit.combined_residual_j, 0.);
}

#[test]
fn generated_heat_refreshes_pinned_viscosity_and_converges_to_coupled_solution() {
    let solid = card(0.01, 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid, false);
    let oil = lubricant([1., 0.1], 500.);
    let mut errors = Vec::new();
    let mut velocities = Vec::new();
    for dt in [0.0625, 0.03125, 0.015625] {
        let old = moving(&reference, &solid, &heat, vec![3.], dt)
            .with_journal_bearings(vec![bearing(&oil, 1, 1)], || Ok(()))
            .unwrap();
        let original = old.journal_bearings().next().unwrap();
        let factor = original.resistance_nm_s;
        // E+CT conserved gives mu=a+b*omega². The continuous coupled ODE
        // has this independently integrated Bernoulli solution for omega².
        let a = 1. - 0.0045 * 75.;
        let b = 0.0045 / 0.12;
        let e = (-6. * factor * a).exp();
        let expected = (a * 9. * e / (a + b * 9. * (1. - e))).sqrt();
        let (next, audit) = old.advance_to(1., config(), || Ok(())).unwrap();
        let temperature = next.thermal().temperature_k(CellId(1)).unwrap();
        let state = next.journal_bearings().next().unwrap();
        let viscosity = 1. - 0.0045 * (temperature - 300.);
        close(state.liquid.dynamic_viscosity_pa_s(), viscosity);
        close(state.resistance_nm_s, factor * viscosity);
        close(audit.viscous_heat_j, 0.02 * (temperature - 300.));
        close(audit.combined_residual_j, 0.);
        for (initial, current) in original
            .liquid
            .material()
            .properties()
            .iter()
            .zip(state.liquid.material().properties())
        {
            assert_eq!(
                initial.answer().receipt.selected,
                current.answer().receipt.selected
            );
        }
        assert_ne!(
            original.liquid.material().identity(),
            state.liquid.material().identity()
        );
        assert!(next.joint_velocity()[0] > 3. * (-3. * factor).exp());
        errors.push((next.joint_velocity()[0] - expected).abs());
        velocities.push(next.joint_velocity()[0]);
        let (half, first) = old.advance_to(0.5, config(), || Ok(())).unwrap();
        let (split, second) = half.advance_to(1., config(), || Ok(())).unwrap();
        assert_eq!(next.joint_position(), split.joint_position());
        assert_eq!(next.joint_velocity(), split.joint_velocity());
        assert_eq!(next.thermal().cells(), split.thermal().cells());
        close(
            audit.viscous_heat_j,
            first.viscous_heat_j + second.viscous_heat_j,
        );
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
    assert!(velocities[0] < velocities[1] && velocities[1] < velocities[2]);
}

#[test]
fn source_domain_and_late_cancellation_leave_mechanics_heat_and_source_state_unchanged() {
    let solid = card(0.01, 0., 500.);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid, false);
    let moving = moving(&reference, &solid, &heat, vec![3.], 0.0625);
    let limited = lubricant([1., 1.], 301.);
    let old = moving
        .with_journal_bearings(vec![bearing(&limited, 1, 1)], || Ok(()))
        .unwrap();
    assert!(matches!(
        old.advance_to(1., config(), || Ok(())),
        Err(Error::Material(_))
    ));
    assert_eq!(old.thermal().time_s(), 0.);
    assert_eq!(old.joint_velocity(), &[3.]);
    close(old.thermal().temperature_k(CellId(1)).unwrap(), 300.);
    close(
        old.journal_bearings()
            .next()
            .unwrap()
            .liquid
            .dynamic_viscosity_pa_s(),
        1.,
    );
    let oil = lubricant([1., 0.1], 500.);
    let old = moving
        .with_journal_bearings(vec![bearing(&oil, 1, 1)], || Ok(()))
        .unwrap();
    let mut calls = 0;
    old.advance_to(1., config(), || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    let snapshot = old.journal_bearings().next().unwrap().clone();
    let mut count = 0;
    assert!(
        old.advance_to(1., config(), || {
            count += 1;
            if count == calls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .is_err()
    );
    assert_eq!(old.joint_position(), &[0.]);
    assert_eq!(old.joint_velocity(), &[3.]);
    assert_eq!(old.thermal().time_s(), 0.);
    assert_eq!(old.journal_bearings().next().unwrap(), &snapshot);
    assert!(
        moving
            .with_journal_bearings(vec![bearing(&oil, 0, 1)], || Ok(()))
            .is_err()
    );
    assert!(
        moving
            .with_journal_bearings(vec![bearing(&oil, 1, 99)], || Ok(()))
            .is_err()
    );
    let mut too_fast = bearing(&oil, 1, 1);
    too_fast.limits.maximum_gap_reynolds = 1e-6;
    assert!(
        moving
            .with_journal_bearings(vec![too_fast], || Ok(()))
            .is_err()
    );
    let mut refine = bearing(&oil, 1, 1);
    refine.limits.maximum_relative_viscosity_change = 1e-6;
    let bounded = moving
        .with_journal_bearings(vec![refine], || Ok(()))
        .unwrap();
    assert!(matches!(
        bounded.advance_to(1., config(), || Ok(())),
        Err(Error::Budget(_))
    ));
}
