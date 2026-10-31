//! G1/G3/G4: original interface resistance changes finite-body heat and frictional motion.
#[path = "../../fs-mbd-material/tests/heating/common.rs"]
mod common;
use common::*;
use fs_conduction::transient::deposition::contact::{
    ContactCell, MaterialContactLink, MaterialContactObserver, MaterialContactObserverError,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
};
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, Vec3};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_mbd::articulated::JointModel;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::thermal::heating::articulated_motion::{
    Controls, HeatedArticulatedMotion,
    dry::{
        DryBrakeLimits, DryJointBrake, JointSlipGeometry, KINETIC_FRICTION_PROPERTY,
        RateDependentDryJointBrake, STATIC_FRICTION_PROPERTY,
    },
};
use fs_mbd_material::thermal::heating::{HeatedArticulated, LinkHeating};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature, Velocity};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::enthalpy::*;
use fs_thermal::observers::ObserverCouplingError;
use fs_tribo::{HeatPartition, InputAuthority};

const LINK: u64 = 77;

// Explicit constant-conductivity fixture, consumed by the native bulk owner.
struct Conductivity(f64);
impl ConductivityLaw for Conductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.0)
    }
}

fn source(cards: [&MaterialCard; 2], resistance: f64, parent_upper: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for name in [STATIC_FRICTION_PROPERTY, KINETIC_FRICTION_PROPERTY, RP] {
        let thermal = name == RP;
        let dims = if thermal { RD } else { Dims::NONE };
        let value = if thermal {
            PropertyValue::Scalar {
                value: resistance,
                dims,
            }
        } else {
            PropertyValue::Surface2D {
                axes: ["TA".into(), "speed".into()],
                axis_dims: [Temperature::DIMS, Velocity::DIMS],
                samples: [[300.0, 0.0], [500.0, 0.0], [300.0, 2.0], [500.0, 2.0]]
                    .into_iter()
                    .map(|p| {
                        (
                            p,
                            if name == STATIC_FRICTION_PROPERTY {
                                3.0
                            } else {
                                0.2 + 0.01 * (p[0] - 300.0) + 0.1 * p[1]
                            },
                        )
                    })
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            }
        };
        let validity = ValidityDomain::unconstrained()
            .with("TA", 300.0, 500.0)
            .with("TB", 300.0, if thermal { parent_upper } else { 500.0 })
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0,
                2.0,
            );
        let validity = if thermal {
            validity
        } else {
            validity.with_quantity("speed", QuantitySpec::dimensional(Velocity::DIMS), 0.0, 2.0)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation: if thermal {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic shared thermal/friction interface, not measured data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "slider".into(),
        },
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "mount".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "fixed surface state".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("TA", QuantitySpec::dimensional(Temperature::DIMS), 310.0)
        .unwrap()
        .with_quantity("TB", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 2.0)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(Velocity::DIMS), 0.0)
        .unwrap()
}

fn interface_point<'a>(
    source: &'a InterfaceSystemCard,
    cards: [&'a MaterialCard; 2],
    context: &'a QueryPoint,
) -> ThermalInterfacePoint<'a> {
    ThermalInterfacePoint {
        interface: source,
        materials: [cards[1], cards[0]],
        context,
        temperature_axes: ["TA", "TB"],
        liquid_fraction_axes: [None, None],
        mean_temperature_axis: None,
    }
}

fn mechanism(cards: [&MaterialCard; 2]) -> ResolvedArticulated {
    let points = [point(300.0, 1e5), point(310.0, 1e5)];
    let regions: Vec<_> = cards
        .iter()
        .enumerate()
        .map(|(i, card)| {
            vec![RegionBinding {
                key: "a",
                material: card,
                point: &points[i],
                geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
                selection: MaterialPropertySelection::SingleClaimOnly,
            }]
        })
        .collect();
    resolve_articulated(
        &[
            MaterialLink {
                name: "mount",
                parent: None,
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::FIXED,
                regions: &regions[0],
            },
            MaterialLink {
                name: "slider",
                parent: Some(0),
                parent_from_child_zero: Se3::identity(),
                joint: JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
                regions: &regions[1],
            },
        ],
        2,
        || false,
    )
    .unwrap()
}

fn moving<'a>(
    reference: &'a ResolvedArticulated,
    cards: [&'a MaterialCard; 2],
    heat: [&'a ResolvedHeatCapacityCurve; 2],
    speed: f64,
    dt: f64,
) -> HeatedArticulatedMotion<'a> {
    let bindings = reference
        .model()
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| {
            let mut region = region(cards[i], heat[i], "a", i as u64, 0.0);
            if i == 1 {
                region.initial_specific_enthalpy_j_kg = 1.0;
            } // cp=.1, T=310 K.
            LinkHeating {
                link: link.name(),
                regions: vec![region],
            }
        })
        .collect();
    let contact = TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CellId(1),
            b: CellId(0),
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    };
    let heated = HeatedArticulated::new(
        reference,
        bindings,
        vec![contact],
        budget(),
        config(),
        || Ok(()),
    )
    .unwrap();
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

fn brake<'a>(
    source: &'a InterfaceSystemCard,
    cards: [&'a MaterialCard; 2],
    context: &'a QueryPoint,
) -> RateDependentDryJointBrake<'a> {
    RateDependentDryJointBrake {
        brake: DryJointBrake {
            link: 1,
            cells: [CellId(1), CellId(0)],
            geometry: JointSlipGeometry::Prismatic,
            normal_force_n: 2.0,
            nominal_area_m2: 1.0,
            pressure_axis: "pressure",
            point: interface_point(source, cards, context),
            selection: MaterialPropertySelection::SingleClaimOnly,
            authority: InputAuthority::SyntheticFixture,
            partition: HeatPartition::new(0.75, 0.25, 0.0).unwrap(),
            limits: DryBrakeLimits {
                regularization_m_s: 0.01,
                maximum_coefficient_change: 0.2,
                absolute_loss_tolerance_j: 1e-12,
            },
        },
        slip_speed_axis: "speed",
        maximum_slip_speed_m_s: 2.0,
    }
}

fn contacts<'a>(
    thermal: &DepositionSimulation<'_>,
    source: &'a InterfaceSystemCard,
    cards: [&'a MaterialCard; 2],
    heat: [&'a ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    conductivity: &'a Conductivity,
) -> MaterialContactObserver<'a, ConductionObserver<'a>> {
    let bulk = ConductionObserver::new(
        thermal,
        &[
            (MaterialId(10), conductivity),
            (MaterialId(11), conductivity),
        ],
        &[ConductiveLink {
            link_id: LINK,
            half_distances_m: [0.25; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    MaterialContactObserver::new(
        bulk,
        thermal,
        vec![MaterialContactLink {
            link: LINK,
            sides: [
                ContactCell {
                    cell: CellId(1),
                    material_id: MaterialId(11),
                    curve: heat[1].curve(),
                },
                ContactCell {
                    cell: CellId(0),
                    material_id: MaterialId(10),
                    curve: heat[0].curve(),
                },
            ],
            point: interface_point(source, cards, context),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }],
        || Ok(()),
    )
    .unwrap()
}

fn step(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}

#[test]
fn native_bulk_and_original_contact_resistance_converge_to_finite_body_exchange() {
    let parent = card(0.2, 0.0, 500.0);
    let child = card(0.1, 0.0, 500.0);
    let cards = [&parent, &child];
    let hp = caloric(&parent, 1e5, SolidLiquidPhase::Solid);
    let hc = caloric(&child, 1e5, SolidLiquidPhase::Solid);
    let heat = [&hp, &hc];
    let reference = mechanism(cards);
    let point = context();
    let source = source(cards, 1.0, 500.0);
    let conductivity = Conductivity(1.0);
    // G=1/(.25/k + R''/A + .25/k)=2/3 W/K; C_A=.2,C_B=.4.
    let g: f64 = 2.0 / 3.0;
    let equilibrium = (0.2 * 310.0 + 0.4 * 300.0) / 0.6;
    let difference = 10.0 * (-g * (1.0 / 0.2 + 1.0 / 0.4) * 0.5).exp();
    let exact = equilibrium + (2.0 / 3.0) * difference;
    let mut errors = Vec::new();
    for dt in [0.03125, 0.015625, 0.0078125] {
        let initial = moving(&reference, cards, heat, 0.0, dt);
        let observer = contacts(
            initial.thermal(),
            &source,
            cards,
            heat,
            &point,
            &conductivity,
        );
        let (next, accepted, audit) = initial
            .thermal()
            .advance_observed(&observer, 0.5, step(dt), || Ok(()))
            .unwrap();
        let ta = next.temperature_k(CellId(1)).unwrap();
        let tb = next.temperature_k(CellId(0)).unwrap();
        errors.push((ta - exact).abs());
        close(0.2 * ta + 0.4 * tb, 0.2 * 310.0 + 0.4 * 300.0);
        close(next.energy_j(), initial.thermal().energy_j());
        close(audit.external_heat_j, 0.0);
        let sample = &accepted.last_samples()[0];
        close(sample.bulk_conductance_w_k, 2.0);
        close(sample.conductance_w_k, g);
        assert_eq!(sample.series.as_ref().unwrap().terms().len(), 2);
        close(accepted.endpoint_sources()[0].phases[0].temperature_k(), ta);
        close(accepted.endpoint_sources()[0].phases[1].temperature_k(), tb);
        assert_eq!(
            accepted.endpoint_sources()[0].interface.receipt().selected,
            observer.endpoint_sources()[0].interface.receipt().selected
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.7 && pair[0] / pair[1] < 2.3,
            "{errors:?}"
        );
    }
}

#[test]
fn resistance_only_data_swap_recouples_cooling_friction_motion_and_restart() {
    let parent = card(0.2, 0.0, 500.0);
    let child = card(0.1, 0.0, 500.0);
    let cards = [&parent, &child];
    let hp = caloric(&parent, 1e5, SolidLiquidPhase::Solid);
    let hc = caloric(&child, 1e5, SolidLiquidPhase::Solid);
    let heat = [&hp, &hc];
    let reference = mechanism(cards);
    let point = context();
    let conductivity = Conductivity(1.0);
    let mut outputs = Vec::new();
    let mut friction_claims = Vec::new();
    for resistance in [1.0, 4.0] {
        let source = source(cards, resistance, 500.0);
        let initial = moving(&reference, cards, heat, 1.5, 0.015625)
            .with_rate_dependent_dry_joint_brakes(vec![brake(&source, cards, &point)], || Ok(()))
            .unwrap();
        let observer = contacts(
            initial.thermal(),
            &source,
            cards,
            heat,
            &point,
            &conductivity,
        );
        let (next, accepted, audit) = initial
            .advance_with(&observer, 0.5, step(0.015625), || Ok(()))
            .unwrap();
        let ta = next.thermal().temperature_k(CellId(1)).unwrap();
        let tb = next.thermal().temperature_k(CellId(0)).unwrap();
        let speed = next.joint_velocity()[0];
        let loss = 2.25 - speed * speed; // independently m/2=1 kg.
        close(audit.friction_heat_j, loss);
        close(
            next.thermal().energy_j() - initial.thermal().energy_j(),
            loss,
        );
        close(audit.thermal.external_heat_j, loss);
        close(audit.combined_residual_j, 0.0);
        let brake = next.dry_joint_brakes().next().unwrap();
        let rate = brake.rate_response.as_ref().unwrap();
        close(rate.resolved.state_coordinate("TA").unwrap(), ta);
        close(rate.resolved.state_coordinate("TB").unwrap(), tb);
        close(rate.resolved.state_coordinate("speed").unwrap(), speed);
        close(
            rate.law.zero_speed_kinetic_coefficient().unwrap(),
            0.2 + 0.01 * (ta - 300.0) + 0.1 * speed,
        );
        close(brake.work.surface_a_heat_j(), 0.75 * loss);
        close(brake.work.surface_b_heat_j(), 0.25 * loss);
        friction_claims.push(
            rate.resolved
                .properties()
                .iter()
                .map(|p| p.answer().receipt.selected)
                .collect::<Vec<_>>(),
        );
        let contact = &accepted.endpoint_sources()[0];
        assert_eq!(contact.interface.card_identity(), source.content_hash());
        assert_eq!(
            contact.interface.receipt().selected,
            observer.endpoint_sources()[0].interface.receipt().selected
        );
        close(contact.phases[0].temperature_k(), ta);
        close(contact.phases[1].temperature_k(), tb);
        close(accepted.last_samples()[0].bulk_conductance_w_k, 2.0);
        close(
            accepted.last_samples()[0].conductance_w_k,
            1.0 / (0.5 + resistance),
        );
        let (prefix, prefix_contact, _) = initial
            .advance_with(&observer, 0.25, step(0.015625), || Ok(()))
            .unwrap();
        let (resumed, resumed_contact, _) = prefix
            .advance_with(&prefix_contact, 0.5, step(0.015625), || Ok(()))
            .unwrap();
        assert_eq!(resumed.joint_velocity(), next.joint_velocity());
        assert_eq!(resumed.joint_position(), next.joint_position());
        assert_eq!(resumed.thermal().cells(), next.thermal().cells());
        assert_eq!(
            resumed_contact.endpoint_sources(),
            accepted.endpoint_sources()
        );
        assert_eq!(resumed_contact.last_samples(), accepted.last_samples());
        outputs.push((ta, speed));
    }
    assert_eq!(
        friction_claims[0], friction_claims[1],
        "only thermal-resistance data changed"
    );
    assert!(
        outputs[0].0 < outputs[1].0,
        "lower contact R cools the slider"
    );
    assert!(
        outputs[0].1 > outputs[1].1,
        "the same friction law senses that cooling"
    );
}

#[test]
fn final_contact_source_exit_and_late_cancellation_preserve_both_owners() {
    let parent = card(0.2, 0.0, 500.0);
    let child = card(0.1, 0.0, 500.0);
    let cards = [&parent, &child];
    let hp = caloric(&parent, 1e5, SolidLiquidPhase::Solid);
    let hc = caloric(&child, 1e5, SolidLiquidPhase::Solid);
    let heat = [&hp, &hc];
    let reference = mechanism(cards);
    let point = context();
    let conductivity = Conductivity(1.0);
    for upper in [300.1, 500.0] {
        let source = source(cards, 1.0, upper);
        let initial = moving(&reference, cards, heat, 1.5, 0.03125)
            .with_rate_dependent_dry_joint_brakes(vec![brake(&source, cards, &point)], || Ok(()))
            .unwrap();
        let observer = contacts(
            initial.thermal(),
            &source,
            cards,
            heat,
            &point,
            &conductivity,
        );
        let thermal = initial.thermal().cells().clone();
        let brake = initial.dry_joint_brakes().next().unwrap().clone();
        let samples = observer.last_samples().to_vec();
        let endpoints = observer.endpoint_sources().to_vec();
        let bulk = observer.inner().last_conductivities_w_m_k().clone();
        if upper < 301.0 {
            // One substep starts inside the source, then exits at its final
            // endpoint; a start-only resistance query would miss this failure.
            assert!(matches!(
                initial.advance_with(&observer, 0.03125, step(0.03125), || Ok(())),
                Err(ObserverCouplingError::First(
                    MaterialContactObserverError::Contact(_)
                ))
            ));
            let (retry, _, _) = initial
                .advance_with(&observer, 1.0 / 4096.0, step(0.03125), || Ok(()))
                .unwrap();
            assert!(retry.thermal().temperature_k(CellId(0)).unwrap() > 300.0);
        } else {
            let mut calls = 0;
            initial
                .advance_with(&observer, 0.5, step(0.03125), || {
                    calls += 1;
                    Ok(())
                })
                .unwrap();
            assert!(calls > 4);
            let mut attempted = 0;
            assert!(
                initial
                    .advance_with(&observer, 0.5, step(0.03125), || {
                        attempted += 1;
                        if attempted == calls - 1 {
                            Err(DepositionError::Cancelled)
                        } else {
                            Ok(())
                        }
                    })
                    .is_err()
            );
        }
        assert_eq!(initial.thermal().cells(), &thermal);
        assert_eq!(initial.joint_velocity(), &[1.5]);
        assert_eq!(initial.joint_position(), &[0.0]);
        assert_eq!(initial.dry_joint_brakes().next().unwrap(), &brake);
        assert_eq!(observer.last_samples(), samples);
        assert_eq!(observer.endpoint_sources(), endpoints);
        assert_eq!(observer.inner().last_conductivities_w_m_k(), &bulk);
        close(initial.thermal().time_s(), 0.0);
    }
}
