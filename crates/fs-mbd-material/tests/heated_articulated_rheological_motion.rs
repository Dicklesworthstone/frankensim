//! Original T/shear viscosity -> finite-gap journal torque -> native motion and heat.
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
    JournalBearing, JournalGeometry, JournalLimits, rheology::RheologicalJournalBearing,
};
use fs_mbd_material::thermal::heating::articulated_motion::{Controls, HeatedArticulatedMotion};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Frequency, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;
use fs_tribo::journal::JournalShearSolve;

const INERTIA: f64 = 1.0 / 3.0;
const SPEED: f64 = 2.0;
const CUTOFF: f64 = 5.0;

fn geometry() -> JournalGeometry {
    JournalGeometry {
        inner_radius_m: 0.1,
        outer_radius_m: 0.11,
        length_m: 0.2,
    }
}

fn lubricant(k: f64, thermal_slope: f64, shear_dependent: bool, upper: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 2.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            1.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            Dims([1, 1, -3, -1, 0, 0]),
            0.5,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, k),
    ] {
        let viscosity = name == DYNAMIC_VISCOSITY_PROPERTY;
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        let response = if viscosity && shear_dependent {
            // Each strip is affine, so source triangulation represents this
            // law exactly: mu = K max(gamma, cutoff) + s*(T-300).
            PropertyValue::Surface2D {
                axes: ["T".into(), "rate".into()],
                axis_dims: [Temperature::DIMS, Frequency::DIMS],
                samples: [
                    [300.0, 0.0],
                    [upper, 0.0],
                    [300.0, CUTOFF],
                    [upper, CUTOFF],
                    [300.0, 100.0],
                    [upper, 100.0],
                ]
                .into_iter()
                .map(|p| (p, k * p[1].max(CUTOFF) + thermal_slope * (p[0] - 300.0)))
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2], [2, 3, 4], [3, 5, 4]],
                dims,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, value),
                    (
                        upper,
                        value
                            + if viscosity {
                                thermal_slope * (upper - 300.0)
                            } else {
                                0.0
                            },
                    ),
                ],
                dims,
            }
        };
        let validity = if viscosity && shear_dependent {
            ValidityDomain::unconstrained()
                .with("T", 300.0, upper)
                .with_quantity(
                    "rate",
                    QuantitySpec::dimensional(Frequency::DIMS),
                    0.0,
                    100.0,
                )
        } else {
            ValidityDomain::unconstrained().with("T", 300.0, upper)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: response,
                validity,
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite-gap rheology fixture, not measured lubricant".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic rheological liquid".into(),
            phase: "liquid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn mechanism(solid: &MaterialCard) -> ResolvedArticulated {
    let p = point(300.0, 1e5);
    let regions = [RegionBinding {
        key: "a",
        material: solid,
        point: &p,
        geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    resolve_articulated(
        &[
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
                joint: JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
                regions: &regions,
            },
        ],
        2,
        || false,
    )
    .unwrap()
}

fn moving<'a>(
    reference: &'a ResolvedArticulated,
    solid: &'a MaterialCard,
    heat: &'a ResolvedHeatCapacityCurve,
    speed: f64,
    dt: f64,
) -> HeatedArticulatedMotion<'a> {
    let bindings = reference
        .model()
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| LinkHeating {
            link: link.name(),
            regions: vec![region(solid, heat, "a", i as u64, 0.0)],
        })
        .collect();
    let heated =
        HeatedArticulated::new(reference, bindings, vec![], budget(), config(), || Ok(())).unwrap();
    HeatedArticulatedMotion::new(
        heated,
        Se3::identity(),
        vec![0.0],
        vec![speed],
        vec![0.0],
        Vec3::default(),
        Controls {
            maximum_step_s: dt,
            relative_state_tolerance: 1e-11,
            absolute_position_tolerance: vec![1e-12],
            absolute_velocity_tolerance: vec![1e-12],
            maximum_rk_attempts_per_substep: 256,
            relative_momentum_tolerance: 1e-11,
            absolute_momentum_tolerance: vec![1e-12],
            maximum_parameter_work_j: 0.0,
            maximum_numerical_energy_change_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
}

fn bearing(oil: &MaterialCard) -> JournalBearing<'_> {
    JournalBearing {
        link: 1,
        heat_cell: CellId(1),
        geometry: geometry(),
        limits: JournalLimits {
            maximum_gap_reynolds: 1.0,
            maximum_relative_viscosity_change: 0.5,
            absolute_loss_tolerance_j: 1e-12,
        },
        liquid: LiquidTemperatureLaw::try_new(
            oil,
            &QueryPoint::new()
                .with("T", 300.0)
                .unwrap()
                .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
                .unwrap(),
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    }
}

fn rheology(oil: &MaterialCard) -> RheologicalJournalBearing<'_> {
    RheologicalJournalBearing {
        bearing: bearing(oil),
        shear_rate_axis: "rate",
        maximum_shear_rate_per_s: 100.0,
        solve: JournalShearSolve {
            maximum_iterations: 128,
            absolute_angular_rate_tolerance_per_s: 1e-13,
            relative_angular_rate_tolerance: 1e-13,
        },
        expansion: None,
    }
}

#[test]
fn shear_data_alone_changes_signed_native_motion_and_deposits_analytic_work_once() {
    let solid = card(10.0, 0.0, 500.0);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid);
    let mut speeds = Vec::new();
    for k in [0.01, 0.03] {
        let oil = lubricant(k, 0.0, true, 500.0);
        for sign in [-1.0, 1.0] {
            let old = moving(&reference, &solid, &heat, sign * SPEED, 0.0625)
                .with_rheological_journal_bearings(vec![rheology(&oil)], || Ok(()))
                .unwrap();
            // Independent exact n=2 Couette solution, integrated through I*w'=-c*w*|w|.
            // Both wall rates remain above cutoff, so the low-rate plateau is inactive.
            let g = geometry();
            let c = 2.0 * core::f64::consts::PI * g.length_m * k
                / (g.inner_radius_m.recip() - g.outer_radius_m.recip()).powi(2);
            let expected = SPEED / (1.0 + c * SPEED / INERTIA);
            let loss = 0.5 * INERTIA * (SPEED.powi(2) - expected.powi(2));
            let (next, audit) = old.advance_to(1.0, config(), || Ok(())).unwrap();
            close(next.joint_velocity()[0], sign * expected);
            close(
                next.joint_position()[0],
                sign * INERTIA / c * (c * SPEED / INERTIA).ln_1p(),
            );
            let state = next.journal_bearings().next().unwrap();
            close(state.resisting_torque_nm, -sign * c * expected.powi(2));
            close(state.dissipated_power_w, c * expected.powi(3));
            close(state.resistance_nm_s, c * expected);
            let shear = state.rheology.as_ref().unwrap();
            let wall_rates = [
                expected * g.outer_radius_m / (g.outer_radius_m - g.inner_radius_m),
                expected * g.inner_radius_m / (g.outer_radius_m - g.inner_radius_m),
            ];
            for (wall, (rate, expected_rate)) in [&shear.inner_liquid, &shear.outer_liquid]
                .into_iter()
                .zip(shear.shear_rate_per_s.into_iter().zip(wall_rates))
            {
                assert!(rate > CUTOFF);
                close(rate, expected_rate);
                close(wall.dynamic_viscosity_pa_s(), k * expected_rate);
                close(wall.material().query().axes()["rate"], expected_rate);
                assert_eq!(wall.material().card_identity(), oil.content_hash());
            }
            assert!(shear.angular_rate_defect_per_s.abs() < 4e-13);
            close(audit.viscous_heat_j, loss);
            close(audit.thermal.external_heat_j, loss);
            close(
                next.thermal().temperature_k(CellId(1)).unwrap(),
                300.0 + loss / 20.0,
            );
            close(next.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
            close(
                audit.mechanical_after.kinetic_j + loss,
                audit.mechanical_before.kinetic_j,
            );
            close(audit.combined_residual_j, 0.0);
            assert!(state.rheology.is_some());
            assert_eq!(old.joint_velocity(), &[sign * SPEED]);
            assert_eq!(old.thermal().time_s(), 0.0);
            if sign > 0.0 {
                speeds.push(next.joint_velocity()[0]);
            }
        }
    }
    assert!(speeds[1] < speeds[0]);
}

#[test]
fn a_constant_source_recovers_existing_newtonian_motion_and_heat() {
    let solid = card(10.0, 0.0, 500.0);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid);
    let oil = lubricant(0.2, 0.0, false, 500.0);
    let initial = moving(&reference, &solid, &heat, SPEED, 0.0625);
    let old = initial
        .with_journal_bearings(vec![bearing(&oil)], || Ok(()))
        .unwrap();
    let curved = initial
        .with_rheological_journal_bearings(vec![rheology(&oil)], || Ok(()))
        .unwrap();
    let (old, a) = old.advance_to(1.0, config(), || Ok(())).unwrap();
    let (curved, b) = curved.advance_to(1.0, config(), || Ok(())).unwrap();
    close(old.joint_position()[0], curved.joint_position()[0]);
    close(old.joint_velocity()[0], curved.joint_velocity()[0]);
    close(a.viscous_heat_j, b.viscous_heat_j);
    close(
        old.thermal().temperature_k(CellId(1)).unwrap(),
        curved.thermal().temperature_k(CellId(1)).unwrap(),
    );
}

#[test]
fn generated_heat_refreshes_original_shear_sources_and_restart_preserves_the_trajectory() {
    let solid = card(0.002, 0.0, 500.0);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid);
    let cold_oil = lubricant(0.01, 0.0, true, 500.0);
    let oil = lubricant(0.01, 0.001, true, 500.0);
    let initial = moving(&reference, &solid, &heat, SPEED, 0.0625);
    let cold = initial
        .with_rheological_journal_bearings(vec![rheology(&cold_oil)], || Ok(()))
        .unwrap();
    let old = initial
        .with_rheological_journal_bearings(vec![rheology(&oil)], || Ok(()))
        .unwrap();
    let (cold, _) = cold.advance_to(1.0, config(), || Ok(())).unwrap();
    let (next, audit) = old.advance_to(1.0, config(), || Ok(())).unwrap();
    let t = next.thermal().temperature_k(CellId(1)).unwrap();
    assert!(t > 310.0);
    assert!(next.joint_velocity()[0] < cold.joint_velocity()[0]);
    close(audit.viscous_heat_j, 0.004 * (t - 300.0));
    close(audit.combined_residual_j, 0.0);
    let state = next.journal_bearings().next().unwrap();
    close(
        state.liquid.dynamic_viscosity_pa_s(),
        0.05 + 0.001 * (t - 300.0),
    );
    close(state.liquid.material().query().axes()["T"], t);
    let shear = state.rheology.as_ref().unwrap();
    for (wall, rate) in [&shear.inner_liquid, &shear.outer_liquid]
        .into_iter()
        .zip(shear.shear_rate_per_s)
    {
        close(wall.material().query().axes()["T"], t);
        close(wall.material().query().axes()["rate"], rate);
        close(
            wall.dynamic_viscosity_pa_s(),
            0.01 * rate.max(CUTOFF) + 0.001 * (t - 300.0),
        );
        for (reference, current) in state
            .liquid
            .material()
            .properties()
            .iter()
            .zip(wall.material().properties())
        {
            assert_eq!(
                reference.answer().receipt.selected,
                current.answer().receipt.selected
            );
        }
    }
    for (before, after) in old
        .journal_bearings()
        .next()
        .unwrap()
        .liquid
        .material()
        .properties()
        .iter()
        .zip(state.liquid.material().properties())
    {
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
    }
    let (half, first) = old.advance_to(0.5, config(), || Ok(())).unwrap();
    let (split, second) = half.advance_to(1.0, config(), || Ok(())).unwrap();
    assert_eq!(next.joint_position(), split.joint_position());
    assert_eq!(next.joint_velocity(), split.joint_velocity());
    assert_eq!(next.thermal().cells(), split.thermal().cells());
    assert_eq!(
        next.journal_bearings().next(),
        split.journal_bearings().next()
    );
    close(
        audit.viscous_heat_j,
        first.viscous_heat_j + second.viscous_heat_j,
    );
}

#[test]
fn rate_support_temperature_support_and_late_cancellation_keep_all_accepted_history() {
    let solid = card(0.002, 0.0, 500.0);
    let heat = caloric(&solid, 1e5, SolidLiquidPhase::Solid);
    let reference = mechanism(&solid);
    let initial = moving(&reference, &solid, &heat, SPEED, 0.0625);
    let oil = lubricant(0.01, 0.001, true, 500.0);
    let mut narrow = rheology(&oil);
    narrow.maximum_shear_rate_per_s = 10.0;
    assert!(
        initial
            .with_rheological_journal_bearings(vec![narrow], || Ok(()))
            .is_err()
    );
    let narrow_temperature = lubricant(0.01, 0.001, true, 302.0);
    let old = initial
        .with_rheological_journal_bearings(vec![rheology(&narrow_temperature)], || Ok(()))
        .unwrap();
    let snapshot = old.journal_bearings().next().unwrap().clone();
    assert!(old.advance_to(1.0, config(), || Ok(())).is_err());
    assert_eq!(old.journal_bearings().next().unwrap(), &snapshot);
    assert_eq!(old.joint_position(), &[0.0]);
    assert_eq!(old.joint_velocity(), &[SPEED]);
    assert_eq!(old.thermal().time_s(), 0.0);
    close(old.thermal().temperature_k(CellId(1)).unwrap(), 300.0);

    let old = initial
        .with_rheological_journal_bearings(vec![rheology(&oil)], || Ok(()))
        .unwrap();
    let snapshot = old.journal_bearings().next().unwrap().clone();
    let mut calls = 0;
    old.advance_to(1.0, config(), || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    let mut count = 0;
    assert!(
        old.advance_to(1.0, config(), || {
            count += 1;
            if count == calls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .is_err()
    );
    assert_eq!(old.journal_bearings().next().unwrap(), &snapshot);
    assert_eq!(old.joint_position(), &[0.0]);
    assert_eq!(old.joint_velocity(), &[SPEED]);
    assert_eq!(old.thermal().time_s(), 0.0);
    close(old.thermal().temperature_k(CellId(1)).unwrap(), 300.0);
}
