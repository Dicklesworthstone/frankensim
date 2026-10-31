//! G1/G3/G4: original shaft materials drive motion, layer heat and thermal feedback.
#[path = "heating/common.rs"]
mod common;
use common::{budget, region};
use fs_evidence::ValidityDomain;
use fs_ga::{Mat3, Se3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section_temperature::ConcentricSectionTemperatureLaw;
use fs_material::state_point::*;
use fs_mbd::articulated::JointModel;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::{
    HeatedArticulated, LinkHeating,
    articulated_motion::{Controls, HeatedArticulatedMotion, torsion::TorsionalJointShaft},
};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;

const G: &str = "shear_modulus";
const ETA: &str = "solid_shear_viscosity";
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
        )
        .unwrap()
        .with("omega", 0.0)
        .unwrap()
}
fn card(
    g: f64,
    eta: f64,
    slopes: [f64; 2],
    upper: f64,
    cp: f64,
    frequency_sample: bool,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (DENSITY_PROPERTY, Density::DIMS, 2.0, 0.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 10.0, 0.0),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            0.0,
            0.0,
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            cp,
            0.0,
        ),
        (G, Pressure::DIMS, g, slopes[0]),
        (ETA, DynViscosity::DIMS, eta, slopes[1]),
    ] {
        let shear = name == G || name == ETA;
        let end = if shear { upper } else { 500.0 };
        let sampled = frequency_sample && name == ETA;
        let value = if sampled {
            PropertyValue::Curve {
                abscissa: "omega".into(),
                abscissa_dims: fs_qty::Frequency::DIMS,
                knots: vec![(0.0, value), (200.0, 2.0 * value)],
                dims,
            }
        } else if slope != 0.0 {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, value), (end, value + slope * (end - 300.0))],
                dims,
            }
        } else {
            PropertyValue::Scalar { value, dims }
        };
        let mut validity = ValidityDomain::unconstrained()
            .with("T", 300.0, end)
            .with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                1e5,
                2e5,
            );
        if shear {
            validity = validity.with("omega", 0.0, 200.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation: if sampled || slope != 0.0 {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic shaft, not measurements".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "layered circular shaft".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(300.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap()
}
fn mechanism(cards: &[&MaterialCard]) -> ResolvedArticulated {
    mechanism_at_temperature(cards, 300.0)
}
fn mechanism_at_temperature(cards: &[&MaterialCard], temperature_k: f64) -> ResolvedArticulated {
    let p = point(temperature_k);
    let make = |i: usize, card| RegionBinding {
        key: ["a", "b"][i],
        material: card,
        point: &p,
        geometry: {
            // Each concentric annulus has length and volume one. About its x
            // axis, J_i=(2i+1)/(2pi); the original density supplies inertia.
            let radii_squared = (2 * i + 1) as f64 / core::f64::consts::PI;
            let transverse = (3.0 * radii_squared + 1.0) / 12.0;
            VolumeMoments::new(
                1.0,
                Vec3::default(),
                Mat3 {
                    m: [
                        0.5 * radii_squared,
                        0.0,
                        0.0,
                        0.0,
                        transverse,
                        0.0,
                        0.0,
                        0.0,
                        transverse,
                    ],
                },
            )
            .unwrap()
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let mount = [make(0, cards[0])];
    let shaft = cards
        .iter()
        .enumerate()
        .map(|(i, &c)| make(i, c))
        .collect::<Vec<_>>();
    resolve_articulated(
        &[
            MaterialLink {
                name: "mount",
                parent: None,
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::FIXED,
                regions: &mount,
            },
            MaterialLink {
                name: "shaft",
                parent: Some(0),
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
                regions: &shaft,
            },
        ],
        3,
        || false,
    )
    .unwrap()
}
fn shaft<'a>(cards: &[&'a MaterialCard]) -> TorsionalJointShaft<'a> {
    shaft_at_temperature(cards, 300.0)
}
fn shaft_at_temperature<'a>(
    cards: &[&'a MaterialCard],
    temperature_k: f64,
) -> TorsionalJointShaft<'a> {
    let requirements = [
        (
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            ScalarAdmissibility::Finite,
        ),
        (G, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive),
        (ETA, DynViscosity::DIMS, ScalarAdmissibility::NonNegative),
    ]
    .map(|(name, dims, range)| ScalarPropertyRequirement::try_new(name, dims, range).unwrap());
    let layers = cards
        .iter()
        .enumerate()
        .map(|(i, card)| ConcentricMaterialLayer {
            outer_radius_m: ((i + 1) as f64 / core::f64::consts::PI).sqrt(),
            material: resolve_material_state_point(
                card,
                &point(temperature_k),
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        })
        .collect();
    let section = ConcentricElasticSection::try_new(0.0, layers).unwrap();
    TorsionalJointShaft {
        link: 1,
        cells: (1..=cards.len()).map(|i| CellId(i as u64)).collect(),
        law: ConcentricSectionTemperatureLaw::try_new(section, cards, &vec!["T"; cards.len()])
            .unwrap(),
        shear_modulus_property: G.into(),
        shear_viscosity_property: ETA.into(),
        reference_length_m: 1.0,
        reference_joint_angle_rad: 0.0,
        maximum_absolute_shear_strain: 0.1,
        maximum_angular_frequency_rad_s: 100.0,
        absolute_loss_tolerance_j: 1e-12,
    }
}
fn controls() -> Controls {
    Controls {
        maximum_step_s: 1.0,
        relative_state_tolerance: 1e-11,
        absolute_position_tolerance: vec![1e-10],
        absolute_velocity_tolerance: vec![1e-10],
        maximum_rk_attempts_per_substep: 2000,
        relative_momentum_tolerance: 1e-10,
        absolute_momentum_tolerance: vec![1e-10],
        maximum_parameter_work_j: 1.0,
        maximum_numerical_energy_change_j: 1e-8,
    }
}
fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-9,
    }
}
fn motion<'a>(
    reference: &'a ResolvedArticulated,
    cards: &[&'a MaterialCard],
    curves: &[&'a ResolvedHeatCapacityCurve],
    watts: &[f64],
    q: f64,
    v: f64,
) -> HeatedArticulatedMotion<'a> {
    motion_at_enthalpy(reference, cards, curves, watts, q, v, 0.0)
}
fn motion_at_enthalpy<'a>(
    reference: &'a ResolvedArticulated,
    cards: &[&'a MaterialCard],
    curves: &[&'a ResolvedHeatCapacityCurve],
    watts: &[f64],
    q: f64,
    v: f64,
    initial_specific_h_j_kg: f64,
) -> HeatedArticulatedMotion<'a> {
    let at_enthalpy = |mut region: fs_mbd_material::thermal::heating::ThermalRegion<'a>| {
        region.initial_specific_enthalpy_j_kg = initial_specific_h_j_kg;
        region
    };
    let child = cards
        .iter()
        .zip(curves)
        .enumerate()
        .map(|(i, (&c, &h))| at_enthalpy(region(c, h, ["a", "b"][i], i as u64 + 1, watts[i])))
        .collect();
    let heated = HeatedArticulated::new(
        reference,
        vec![
            LinkHeating {
                link: "mount",
                regions: vec![at_enthalpy(region(cards[0], curves[0], "a", 0, 0.0))],
            },
            LinkHeating {
                link: "shaft",
                regions: child,
            },
        ],
        vec![],
        budget(),
        step(),
        || Ok(()),
    )
    .unwrap();
    HeatedArticulatedMotion::new(
        heated,
        Se3::identity(),
        vec![q],
        vec![v],
        vec![0.0],
        Vec3::default(),
        controls(),
        || Ok(()),
    )
    .unwrap()
}
fn oscillator(q0: f64, v0: f64, inertia: f64, k: f64, c: f64, dt: f64) -> (f64, f64) {
    let beta = c / (2.0 * inertia);
    let omega = (k / inertia - beta * beta).sqrt();
    let b = (v0 + beta * q0) / omega;
    let e = (-beta * dt).exp();
    let (s, c) = (omega * dt).sin_cos();
    (
        e * (q0 * c + b * s),
        e * (v0 * c - (beta * b + omega * q0) * s),
    )
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

#[test]
fn material_shaft_storage_and_loss_drive_native_motion_and_each_finite_layer() {
    let mut responses = Vec::new();
    for (moduli, viscosities) in [
        ([8.0, 8.0], [0.0, 0.0]),
        ([8.0, 8.0], [0.4, 0.4]),
        ([8.0, 16.0], [0.4, 1.2]),
    ] {
        let cards = [
            card(moduli[0], viscosities[0], [0.0; 2], 500.0, 10.0, false),
            card(moduli[1], viscosities[1], [0.0; 2], 500.0, 20.0, false),
        ];
        let curves = cards.each_ref().map(caloric);
        let references = cards.each_ref();
        let reference = mechanism(&references);
        let q0 = 0.05;
        let v0 = 0.01;
        let old = motion(
            &reference,
            &references,
            &curves.each_ref(),
            &[0.0; 2],
            q0,
            v0,
        )
        .with_torsional_joint_shafts(vec![shaft(&references)], || Ok(()))
        .unwrap();
        let polar = [1.0 / core::f64::consts::TAU, 3.0 / core::f64::consts::TAU];
        let inertia = 2.0 * (polar[0] + polar[1]);
        let k = moduli[0] * polar[0] + moduli[1] * polar[1];
        let c = viscosities[0] * polar[0] + viscosities[1] * polar[1];
        let dt = 0.4;
        let (q, v) = oscillator(q0, v0, inertia, k, c, dt);
        let loss = 0.5 * inertia * (v0 * v0 - v * v) + 0.5 * k * (q0 * q0 - q * q);
        let (next, audit) = old.advance_to(dt, step(), || Ok(())).unwrap();
        near(next.joint_position()[0], q, 2e-9);
        near(next.joint_velocity()[0], v, 2e-9);
        near(audit.viscous_heat_j, loss, 2e-9);
        near(audit.parameter_work_j, 0.0, 1e-14);
        near(
            next.thermal().energy_j() - old.thermal().energy_j(),
            loss,
            2e-9,
        );
        near(audit.combined_residual_j, 0.0, 2e-9);
        for i in 0..2 {
            let expected_heat = if c == 0.0 {
                0.0
            } else {
                loss * viscosities[i] * polar[i] / c
            };
            let cell = &next.thermal().cells()[&CellId(i as u64 + 1)];
            near(cell.mass_kg * cell.specific_h_j_kg, expected_heat, 2e-9);
        }
        let state = next.torsional_joint_shafts().next().unwrap();
        near(state.elastic_energy_j, 0.5 * k * q * q, 2e-9);
        near(state.dissipated_power_w, c * v * v, 2e-9);
        for (i, layer) in state.section.shear().elastic().layers().iter().enumerate() {
            assert_eq!(layer.material.card_identity(), cards[i].content_hash());
        }
        assert_eq!(old.joint_position(), &[q0]);
        responses.push(next.joint_position()[0]);
    }
    assert!(responses[0] != responses[1] && responses[1] != responses[2]);
}

#[test]
fn finite_heat_refreshes_original_stiffness_and_shear_loss_before_continuing() {
    let card = card(8.0, 0.4, [-2.0, 2.0], 303.0, 10.0, false);
    let curve = caloric(&card);
    let cards = [&card];
    let reference = mechanism(&cards);
    let old = motion(&reference, &cards, &[&curve], &[20.0], 0.05, 0.0)
        .with_torsional_joint_shafts(vec![shaft(&cards)], || Ok(()))
        .unwrap();
    let polar = 1.0 / core::f64::consts::TAU;
    let inertia = 2.0 * polar;
    let dt = 0.2;
    let (q1, v1) = oscillator(0.05, 0.0, inertia, 8.0 * polar, 0.4 * polar, dt);
    let loss1 = 0.5 * 8.0 * polar * (0.05_f64.powi(2) - q1 * q1) - 0.5 * inertia * v1 * v1;
    let t1 = 300.0 + (20.0 * dt + loss1) / 20.0;
    let k1 = (8.0 - 2.0 * (t1 - 300.0)) * polar;
    let c1 = (0.4 + 2.0 * (t1 - 300.0)) * polar;
    let (one, first) = old.advance_to(dt, step(), || Ok(())).unwrap();
    let current = one.torsional_joint_shafts().next().unwrap();
    near(current.phases[0].temperature_k(), t1, 2e-10);
    near(current.stiffness_nm_rad, k1, 2e-10);
    near(current.viscous_resistance_nm_s, c1, 2e-10);
    near(
        first.parameter_work_j,
        0.5 * (k1 - 8.0 * polar) * q1 * q1,
        2e-9,
    );
    let (q2, v2) = oscillator(q1, v1, inertia, k1, c1, dt);
    let loss2 = 0.5 * inertia * (v1 * v1 - v2 * v2) + 0.5 * k1 * (q1 * q1 - q2 * q2);
    let (two, second) = one.advance_to(2.0 * dt, step(), || Ok(())).unwrap();
    near(two.joint_position()[0], q2, 2e-9);
    near(two.joint_velocity()[0], v2, 2e-9);
    near(second.viscous_heat_j, loss2, 2e-9);
    near(
        two.thermal().energy_j() - one.thermal().energy_j(),
        20.0 * dt + loss2,
        2e-9,
    );
    near(second.combined_residual_j, 0.0, 2e-9);
}

#[test]
fn cooling_below_reference_preserves_supported_shear_sources_and_parameter_work() {
    let card = card(8.0, 0.4, [-2.0, 2.0], 303.0, 10.0, false);
    let curve = caloric(&card);
    let cards = [&card];
    let reference_temperature = 301.0;
    let reference = mechanism_at_temperature(&cards, reference_temperature);
    let old = motion_at_enthalpy(&reference, &cards, &[&curve], &[-20.0], 0.05, 0.0, 10.0)
        .with_torsional_joint_shafts(
            vec![shaft_at_temperature(&cards, reference_temperature)],
            || Ok(()),
        )
        .unwrap();
    let polar = 1.0 / core::f64::consts::TAU;
    let inertia = 2.0 * polar;
    let k = 6.0 * polar;
    let c = 2.4 * polar;
    let dt = 0.2;
    let (q, v) = oscillator(0.05, 0.0, inertia, k, c, dt);
    let loss = 0.5 * k * (0.05_f64.powi(2) - q * q) - 0.5 * inertia * v * v;
    let temperature = reference_temperature + (-20.0 * dt + loss) / 20.0;
    assert!(temperature > 300.0 && temperature < reference_temperature);
    let (next, audit) = old.advance_to(dt, step(), || Ok(())).unwrap();
    let current = next.torsional_joint_shafts().next().unwrap();
    let refreshed_k = (8.0 - 2.0 * (temperature - 300.0)) * polar;
    near(next.joint_position()[0], q, 2e-9);
    near(next.joint_velocity()[0], v, 2e-9);
    near(current.phases[0].temperature_k(), temperature, 2e-10);
    near(current.stiffness_nm_rad, refreshed_k, 2e-10);
    near(
        current.viscous_resistance_nm_s,
        (0.4 + 2.0 * (temperature - 300.0)) * polar,
        2e-10,
    );
    near(audit.viscous_heat_j, loss, 2e-9);
    near(
        audit.parameter_work_j,
        0.5 * (refreshed_k - k) * q * q,
        2e-9,
    );
    near(
        next.thermal().energy_j() - old.thermal().energy_j(),
        -20.0 * dt + loss,
        2e-9,
    );
    near(audit.combined_residual_j, 0.0, 2e-9);
    let reference_query = point(reference_temperature);
    let material = &current.section.shear().elastic().layers()[0].material;
    assert_eq!(material.card_identity(), card.content_hash());
    assert_eq!(
        material.query().axis_quantities(),
        reference_query.axis_quantities()
    );
    for name in ["pressure", "omega"] {
        assert_eq!(material.query().axes()[name], reference_query.axes()[name]);
    }
}

#[test]
fn sampled_frequency_laws_incomplete_bands_and_duplicate_material_inventory_refuse() {
    for sampled in [false, true] {
        let card = card(8.0, 0.4, [0.0; 2], 500.0, 10.0, sampled);
        let curve = caloric(&card);
        let cards = [&card];
        let reference = mechanism(&cards);
        let old = motion(&reference, &cards, &[&curve], &[0.0], 0.05, 0.0);
        if sampled {
            assert!(
                old.with_torsional_joint_shafts(vec![shaft(&cards)], || Ok(()))
                    .is_err()
            );
        } else {
            let mut bad = shaft(&cards);
            bad.maximum_angular_frequency_rad_s = 300.0;
            assert!(
                old.with_torsional_joint_shafts(vec![bad], || Ok(()))
                    .is_err()
            );
            let mut bad = shaft(&cards);
            bad.reference_length_m = 2.0;
            assert!(
                old.with_torsional_joint_shafts(vec![bad], || Ok(()))
                    .is_err()
            );
            assert!(
                old.with_torsional_joint_shafts(vec![shaft(&cards), shaft(&cards)], || Ok(()))
                    .is_err()
            );
            let mut bad = shaft(&cards);
            bad.maximum_absolute_shear_strain = 0.001;
            assert!(
                old.with_torsional_joint_shafts(vec![bad], || Ok(()))
                    .is_err()
            );
        }
        assert_eq!(old.joint_position(), &[0.05]);
        assert_eq!(old.thermal().time_s(), 0.0);
    }
}

#[test]
fn thermal_source_exit_and_late_cancellation_keep_motion_heat_and_source_history_atomic() {
    let card = card(8.0, 0.4, [0.0; 2], 300.1, 10.0, false);
    let curve = caloric(&card);
    let cards = [&card];
    let reference = mechanism(&cards);
    let old = motion(&reference, &cards, &[&curve], &[20.0], 0.05, 0.0)
        .with_torsional_joint_shafts(vec![shaft(&cards)], || Ok(()))
        .unwrap();
    assert!(old.advance_to(0.2, step(), || Ok(())).is_err());
    let mut polls = 0;
    let (short, _) = old
        .advance_to(0.01, step(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut cancelled = 0;
    assert!(
        old.advance_to(0.01, step(), || {
            cancelled += 1;
            if cancelled == polls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .is_err()
    );
    let (retry, _) = old.advance_to(0.01, step(), || Ok(())).unwrap();
    assert_eq!(short.thermal().cells(), retry.thermal().cells());
    assert_eq!(short.joint_position(), retry.joint_position());
    assert_eq!(short.joint_velocity(), retry.joint_velocity());
    assert_eq!(old.thermal().energy_j(), 0.0);
    assert_eq!(old.joint_position(), &[0.05]);
    assert_eq!(
        old.torsional_joint_shafts().next().unwrap().phases[0].temperature_k(),
        300.0
    );
}
