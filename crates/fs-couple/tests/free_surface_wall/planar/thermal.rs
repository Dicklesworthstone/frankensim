//! Finite recoil, actual thermal trajectories, and one-time source accounting.

use super::*;
use fs_couple::free_surface_wall::planar::thermal::FreeSurfaceThermalPlanarWall;
use fs_evidence::ValidityDomain;
use fs_lbm::material::free_surface::MaterialFreeSurfaceDomain3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat, FreeSurfaceEnthalpyControls, MaterialFreeSurfaceEnthalpy3,
};
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

fn card(hot_viscosity: f64, hot_density: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), [1000.0, hot_density]),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [2000.0; 2],
        ),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), [0.5; 2]),
        (
            "dynamic_viscosity",
            Dims([-1, 1, -1, 0, 0, 0]),
            [0.1, hot_viscosity],
        ),
    ] {
        let quantity = if name == "thermal_conductivity" {
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
                    knots: vec![(300.0, values[0]), (500.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 500.0)
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic thermal finite-wall fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "thermal finite-wall fixture".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn with_thermal<R>(
    wall: PlanarWall3,
    hot_viscosity: f64,
    buoyant: bool,
    exercise: impl FnOnce(MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = card(hot_viscosity, if buoyant { 980.0 } else { 1000.0 });
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap();
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        [
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 600000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ]
        .to_vec(),
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    let gravity = if buoyant {
        [0.01, 0.0, -0.02]
    } else {
        [0.0; 3]
    };
    let flow = plan_free_surface_from_material(&source, DX, DT, 0.1)
        .unwrap()
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                gravity_m_s2: gravity,
                ..domain()
            },
            |x, y, z| {
                let coordinate = [x, y, z][axis(wall)];
                if upper(wall) {
                    coordinate >= 4
                } else {
                    coordinate < 4
                }
            },
        )
        .unwrap();
    let h = |p: [usize; 3]| 100000.0 + 5000.0 * p[0] as f64 + 3000.0 * p[1] as f64;
    let flow = if buoyant {
        flow.with_boussinesq_enthalpy(
            &law,
            &caloric,
            FreeSurfaceEnthalpyControls::default(),
            BoussinesqConfig {
                gravity_m_s2: gravity,
                max_relative_density_departure: 0.03,
                max_relative_specific_heat_departure: 0.0,
            },
            h,
            &gate,
        )
    } else {
        flow.with_enthalpy(
            &law,
            &caloric,
            FreeSurfaceEnthalpyControls::default(),
            h,
            &gate,
        )
    }
    .unwrap()
    .with_time_origin(37.0)
    .unwrap();
    exercise(flow, &gate)
}

fn body_parameters(wall: PlanarWall3, ratio: f64) -> (MassProperties, RigidBodyState) {
    let normal = axis(wall);
    let mass = ratio * IMPEDANCE;
    let spin = 0.6 * ratio * SPIN_IMPEDANCE;
    let mut center = [0.0; 3];
    center[(normal + 1) % 3] = 1.5 * DX;
    center[(normal + 2) % 3] = 2.3 * DX;
    center[normal] = if upper(wall) { 10.0 * DX } else { -2.0 * DX };
    (
        inertia(mass, spin, normal),
        state(mass, spin, [0.015, -0.007, 1.3], vector(center), normal),
    )
}

#[test]
fn g1_thermal_constant_properties_reproduce_original_finite_wall_on_all_planes() {
    for wall in [
        PlanarWall3::XMin,
        PlanarWall3::XMax,
        PlanarWall3::YMin,
        PlanarWall3::YMax,
        PlanarWall3::ZMin,
        PlanarWall3::ZMax,
    ] {
        with_thermal(wall, 0.1, false, |liquid, gate| {
            let (properties, body) = body_parameters(wall, 3.0);
            let mut plain = FreeSurfacePlanarWall::new(
                liquid.flow().clone(),
                properties,
                body,
                wall,
                CELLS,
                gate,
            )
            .unwrap();
            let h0 = liquid.total_enthalpy_j().unwrap();
            let mut thermal =
                FreeSurfaceThermalPlanarWall::new(liquid, properties, body, wall, CELLS, gate)
                    .unwrap();
            let force = tangent([0.0003, -0.0002, 0.0], axis(wall));
            for index in 1..=2 {
                let reference = plain.step_candidate(force, 1e-7, controls(), gate).unwrap();
                let actual = thermal
                    .step_candidate(force, 1e-7, &[], controls(), gate)
                    .unwrap();
                assert_eq!(actual.motion.body(), reference.motion.body());
                assert_eq!(actual.wall, reference.wall);
                assert_eq!(
                    actual.motion.liquid().flow().lattice().grid().state_bits(),
                    reference.motion.liquid().lattice().grid().state_bits()
                );
                assert_eq!(actual.native_trials, 2);
                assert_eq!(actual.motion.accepted_steps(), index);
                assert_eq!(actual.motion.liquid().accepted_steps(), index);
                near(
                    actual.motion.liquid().physical_time_s(),
                    37.0 + index as f64 * DT,
                    1e-14,
                );
                near(
                    actual.motion.liquid().total_enthalpy_j().unwrap(),
                    h0,
                    2e-10,
                );
                assert!(actual.liquid_report.molecular_heat.is_none());
                assert_eq!(actual.liquid_report.enthalpy.external_heat_j, 0.0);
                assert_eq!(actual.liquid_report.enthalpy.viscous_heat_j, 0.0);
                plain = reference.motion;
                thermal = actual.motion;
            }
        });
    }
}

#[test]
fn g1_finite_recoil_uses_actual_heated_buoyant_flow_and_separate_work_accounts() {
    let wall = PlanarWall3::ZMin;
    for buoyant in [false, true] {
        with_thermal(wall, 0.02, buoyant, |liquid, gate| {
            let (properties, initial) = body_parameters(wall, 0.7);
            let mut motion =
                FreeSurfaceThermalPlanarWall::new(liquid, properties, initial, wall, CELLS, gate)
                    .unwrap();
            let force = Vec3::new(0.0003, -0.0002, 0.0);
            let heat = [FreeSurfaceCellHeat {
                coordinates: [3, 3, 1],
                energy_j: 0.01,
            }];
            for index in 1..=3 {
                let before = motion.body();
                let before_h = motion.liquid().total_enthalpy_j().unwrap();
                let before_mass = motion.liquid().flow().ledger_mass_kg().unwrap();
                let step = motion
                    .step_candidate(force, 1e-7, &heat, controls(), gate)
                    .unwrap();
                let after = step.motion.body();
                let linear = vector(step.wall.wall_impulse_kg_m_s)
                    .add(step.raw_guide_impulse_world)
                    .add(step.external_linear_impulse_world);
                let couple = vector(step.wall.wall_angular_impulse_kg_m2_s)
                    .add(step.raw_guide_angular_impulse_world)
                    .add(step.external_couple_impulse_world);
                close_vector(
                    after
                        .linear_momentum_world()
                        .sub(before.linear_momentum_world()),
                    linear,
                    2e-17,
                );
                close_vector(
                    angular(after).sub(angular(before)),
                    step.coupling_origin_world_m.cross(linear).add(couple),
                    2e-20,
                );
                let dk = after.kinetic_energy(properties).unwrap()
                    - before.kinetic_energy(properties).unwrap();
                near(dk, step.kinetic_energy_change_j, 2e-18);
                near(
                    dk + step.wall.fluid_work_j - step.external_work_j - step.raw_guide_work_j,
                    0.0,
                    2e-18,
                );
                assert!(step.wall.fluid_work_j.abs() > 0.0);
                near(
                    step.motion.liquid().total_enthalpy_j().unwrap() - before_h,
                    0.01,
                    2e-10,
                );
                near(
                    step.motion.liquid().flow().ledger_mass_kg().unwrap(),
                    before_mass,
                    2e-16,
                );
                assert_eq!(step.liquid_report.enthalpy.external_heat_j, 0.01);
                assert_eq!(step.liquid_report.enthalpy.viscous_heat_j, 0.0);
                assert!(step.liquid_report.enthalpy.mass_routes > 0);
                assert_eq!(step.motion.liquid().accepted_steps(), index);
                let h = step
                    .motion
                    .liquid()
                    .specific_enthalpy_at([3, 3, 1])
                    .unwrap();
                let mu = step
                    .motion
                    .liquid()
                    .material_at([3, 3, 1])
                    .unwrap()
                    .dynamic_viscosity_pa_s();
                near(mu, 0.1 - 0.08 * h / 400000.0, 1e-14);

                // Actual native readback at the solved motion must reproduce
                // the receiving thermal trajectory, not just a predicted load.
                let (replayed, replay) = motion
                    .liquid()
                    .planar_rigid_wall_motion_candidate_si(
                        step.wall.motion,
                        controls().maximum_relative_mass_defect,
                        &heat,
                        gate,
                    )
                    .unwrap();
                assert_eq!(replay.wall, step.wall);
                assert_eq!(replay.flow, step.liquid_report.enthalpy);
                assert_eq!(
                    replayed.flow().lattice().grid().state_bits(),
                    step.motion.liquid().flow().lattice().grid().state_bits()
                );
                assert_eq!(
                    replayed.total_enthalpy_j().unwrap(),
                    step.motion.liquid().total_enthalpy_j().unwrap()
                );
                assert_eq!(
                    replayed.material_at([3, 3, 1]),
                    step.motion.liquid().material_at([3, 3, 1])
                );
                motion = step.motion;
            }
        });
    }
}

#[test]
fn g3_molecular_source_is_original_and_applied_once_despite_the_pilot() {
    with_thermal(PlanarWall3::ZMin, 0.02, false, |liquid, gate| {
        let (properties, body) = body_parameters(PlanarWall3::ZMin, 3.0);
        let motion = FreeSurfaceThermalPlanarWall::new(
            liquid,
            properties,
            body,
            PlanarWall3::ZMin,
            CELLS,
            gate,
        )
        .unwrap()
        .step_candidate(Vec3::ZERO, 0.0, &[], controls(), gate)
        .unwrap()
        .motion;
        let mut original_heat_owner = motion.liquid().clone();
        let reference = original_heat_owner
            .step_with_viscous_heating(&[], gate)
            .unwrap();
        let molecular: f64 = reference.applied.iter().map(|cell| cell.heat_j).sum();
        assert!(molecular > 0.0);
        let authored: Vec<_> = reference
            .applied
            .iter()
            .map(|cell| FreeSurfaceCellHeat {
                coordinates: cell.coordinates,
                energy_j: cell.heat_j,
            })
            .collect();
        let implicit = motion
            .step_candidate_with_viscous_heating(Vec3::ZERO, 0.0, &[], controls(), gate)
            .unwrap();
        let explicit = motion
            .step_candidate(Vec3::ZERO, 0.0, &authored, controls(), gate)
            .unwrap();
        assert_eq!(
            implicit.liquid_report.molecular_heat.as_ref().unwrap(),
            &reference.applied
        );
        assert_eq!(implicit.liquid_report.enthalpy.viscous_heat_j, molecular);
        assert_eq!(implicit.liquid_report.enthalpy.external_heat_j, 0.0);
        near(
            implicit.motion.liquid().total_enthalpy_j().unwrap()
                - motion.liquid().total_enthalpy_j().unwrap(),
            molecular,
            2e-10,
        );
        assert_eq!(implicit.motion.body(), explicit.motion.body());
        assert_eq!(implicit.wall, explicit.wall);
        assert_eq!(
            implicit
                .motion
                .liquid()
                .flow()
                .lattice()
                .grid()
                .state_bits(),
            explicit
                .motion
                .liquid()
                .flow()
                .lattice()
                .grid()
                .state_bits()
        );
        assert_eq!(
            implicit.motion.liquid().total_enthalpy_j().unwrap(),
            explicit.motion.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(
            implicit.motion.liquid().accepted_steps(),
            motion.liquid().accepted_steps() + 1
        );
        assert_eq!(
            implicit.motion.accepted_steps(),
            motion.accepted_steps() + 1
        );
    });
}

#[test]
fn g4_late_mechanical_or_source_refusal_preserves_both_owners_and_retry() {
    with_thermal(PlanarWall3::ZMin, 0.02, true, |liquid, gate| {
        let (properties, body) = body_parameters(PlanarWall3::ZMin, 3.0);
        let motion = FreeSurfaceThermalPlanarWall::new(
            liquid,
            properties,
            body,
            PlanarWall3::ZMin,
            CELLS,
            gate,
        )
        .unwrap();
        let reference = motion.clone();
        let original = (
            motion.body(),
            motion.liquid().flow().lattice().grid().state_bits(),
            motion.liquid().total_enthalpy_j().unwrap(),
            motion.liquid().material_at([3, 3, 1]).cloned(),
            motion.liquid().physical_time_s(),
            motion.liquid().accepted_steps(),
            motion.elapsed_time_s(),
            motion.accepted_steps(),
        );
        let heat = [FreeSurfaceCellHeat {
            coordinates: [3, 3, 1],
            energy_j: 0.01,
        }];
        let mut limited = controls();
        limited.maximum_trial_cells -= 1;
        assert!(matches!(
            motion.step_candidate(Vec3::ZERO, 0.0, &heat, limited, gate),
            Err(PlanarWallError::WorkBudget { .. })
        ));
        limited = controls();
        limited.maximum_absolute_wall_work_j = 0.0;
        // This gate follows both complete thermal trials, including heat,
        // transported source/rheology state and the final native wall receipt.
        assert!(matches!(
            motion.step_candidate(Vec3::ZERO, 0.0, &heat, limited, gate),
            Err(PlanarWallError::Invalid(_))
        ));
        assert!(matches!(
            motion.step_candidate(
                Vec3::ZERO,
                0.0,
                &[FreeSurfaceCellHeat {
                    coordinates: [3, 3, 1],
                    energy_j: 10.0
                }],
                controls(),
                gate
            ),
            Err(PlanarWallError::Thermal(_))
        ));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            motion.step_candidate(Vec3::ZERO, 0.0, &heat, controls(), &cancelled),
            Err(PlanarWallError::Cancelled)
        ));
        assert_eq!(
            (
                motion.body(),
                motion.liquid().flow().lattice().grid().state_bits(),
                motion.liquid().total_enthalpy_j().unwrap(),
                motion.liquid().material_at([3, 3, 1]).cloned(),
                motion.liquid().physical_time_s(),
                motion.liquid().accepted_steps(),
                motion.elapsed_time_s(),
                motion.accepted_steps(),
            ),
            original
        );
        let retry = motion
            .step_candidate(Vec3::ZERO, 0.0, &heat, controls(), gate)
            .unwrap();
        let fresh = reference
            .step_candidate(Vec3::ZERO, 0.0, &heat, controls(), gate)
            .unwrap();
        assert_eq!(retry.motion.body(), fresh.motion.body());
        assert_eq!(retry.wall, fresh.wall);
        assert_eq!(retry.liquid_report, fresh.liquid_report);
        assert_eq!(
            retry.motion.liquid().flow().lattice().grid().state_bits(),
            fresh.motion.liquid().flow().lattice().grid().state_bits()
        );
    });
}
