//! Reciprocal drag stays on the original thermal and constitutive trajectories.
//! Synthetic source data exercise ownership and mechanics, not measured validity.

use super::*;
use crate::liquid_momentum::{a, v};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, Cx, ExecMode, StreamKey};
use fs_lbm::material::free_surface::MaterialFreeSurfaceInitialCell3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceEnthalpyControls, MaterialFreeSurfaceEnthalpy3,
};
use fs_lbm::material::free_surface::rheology::{
    IsothermalRheologyLaw, MaterialFreeSurfaceRheology3, RheologySolve,
};
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_qty::{Dims, Frequency, Temperature};

const CELL: [usize; 3] = [2, 1, 1];
const INDEX: usize = 22;
const DX: f64 = 0.001;
const DT: f64 = 0.0002;

fn evolving_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, scalar) in [
        ("density", Density::DIMS, 1000.0),
        ("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0]), 2000.0),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), 0.5),
        ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS, 0.1),
    ] {
        let viscosity = name == "dynamic_viscosity";
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        let value = if viscosity {
            PropertyValue::Surface2D {
                axes: ["T".into(), "rate".into()],
                axis_dims: [Temperature::DIMS, Frequency::DIMS],
                samples: [[300.0, 0.0], [500.0, 0.0], [300.0, 1000.0], [500.0, 1000.0]]
                    .into_iter()
                    .map(|p| (p, 0.1 + 0.0002 * (p[0] - 400.0) + 0.003 * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            }
        } else {
            PropertyValue::Scalar {
                value: scalar,
                dims,
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value,
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 500.0)
                    .with("pressure", 100_000.0, 100_000.0)
                    .with_quantity(
                        "rate",
                        QuantitySpec::dimensional(Frequency::DIMS),
                        0.0,
                        1000.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if viscosity {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic evolving-liquid reciprocal momentum fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "reciprocal source fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn with_liquids(
    exercise: impl FnOnce(
        &mut MaterialFreeSurfaceRheology3<'_, '_>,
        &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
        &CancelGate,
    ),
) {
    let card = evolving_card();
    let point = QueryPoint::new()
        .with("T", 400.0)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let isothermal = IsothermalRheologyLaw::try_new(
        &card,
        law.reference_state().material().clone(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        "T",
        "rate",
    )
    .unwrap();
    let source = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(),
        LiquidViscosityBasis::Dynamic,
        "T",
        "rate",
    )
    .unwrap();
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
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
                specific_enthalpy_j_kg: 400_000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ]
        .to_vec(),
    )
    .unwrap();
    let domain = MaterialFreeSurfaceDomain3 {
        cells: [4; 3],
        gravity_m_s2: [0.03, -0.02, 0.01],
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 64,
    };
    let mut cells = Vec::new();
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..3 {
                cells.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: if x == 2 { 0.5 } else { 1.0 },
                    velocity_m_s: [0.001 * z as f64, -0.001 * y as f64, 0.0],
                });
            }
        }
    }
    let solve = RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    };
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let mut rheological = MaterialFreeSurfaceRheology3::from_cells(
            &isothermal,
            DX,
            DT,
            0.1,
            domain,
            &cells,
            solve,
            &cx,
            &gate,
        )
        .unwrap();
        let mut thermal = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_rheology(
                domain,
                &cells,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                solve,
                |p| 200_000.0 + 1000.0 * p[0] as f64,
                &cx,
                &gate,
            )
            .unwrap()
            .with_time_origin(7.0)
            .unwrap();
        // Exercise accepted constitutive/thermal history, not just fresh wrappers.
        rheological.step_with_gate(&gate).unwrap();
        thermal.step(&gate).unwrap();
        exercise(&mut rheological, &mut thermal, &gate);
    });
}

fn body_at_cell() -> (RigidBodyState, MassProperties) {
    let mass = 3e-6;
    let properties = MassProperties::new(mass, Vec3::ZERO, Vec3::new(3e-12, 4e-12, 5e-12)).unwrap();
    let position = Vec3::new(
        (CELL[0] as f64 + 0.5) * DX,
        (CELL[1] as f64 + 0.5) * DX,
        (CELL[2] as f64 + 0.5) * DX,
    );
    let body = RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.04, 0.01, -0.02).scale(mass),
        Vec3::ZERO,
    )
    .unwrap();
    (body, properties)
}

fn drag<L: MaterialFreeSurfaceMomentum3 + ?Sized>(
    liquid: &mut L,
    rigid: &mut RigidBodyState,
    properties: MassProperties,
    directional: bool,
    gate: &CancelGate,
) -> (Vec3, f64) {
    if directional {
        let report = apply_heterogeneous_rigid_drag(
            liquid,
            rigid,
            properties,
            &[INDEX],
            &[heterogeneous::PrincipalLiquidDrag3 {
                rates_per_s: [3.0, 0.0, 0.0],
                axes_to_world: UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.4)
                    .unwrap(),
            }],
            Vec3::ZERO,
            0.2,
            heterogeneous::HeterogeneousRigidDragControls {
                exchange: controls(),
                substeps: 1,
                maximum_pair_exchanges: 1,
            },
            gate,
        )
        .unwrap();
        (
            report.exchange.liquid_impulse_kg_m_s,
            report.exchange.drag_dissipation_j,
        )
    } else {
        let report = apply_rigid_drag(
            liquid,
            rigid,
            properties,
            &[INDEX],
            Vec3::ZERO,
            3.0,
            0.2,
            controls(),
            gate,
        )
        .unwrap();
        (report.liquid_impulse_kg_m_s, report.drag_dissipation_j)
    }
}

#[allow(clippy::too_many_arguments)]
fn check_recoil(
    before: RigidBodyState,
    after: RigidBodyState,
    properties: MassProperties,
    mass: f64,
    initial: [f64; 3],
    final_velocity: [f64; 3],
    impulse: Vec3,
    loss: f64,
    directional: bool,
) {
    let velocity = before.center_of_mass_velocity_world(properties).unwrap();
    let mut slip = velocity.sub(v(initial));
    if directional {
        let axis = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.4)
            .unwrap()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        slip = axis.scale(axis.dot(slip));
    }
    let reduced = mass * properties.mass() / (mass + properties.mass());
    let exposure = 0.6 * (1.0 + mass / properties.mass());
    let expected = slip.scale(reduced * -(-exposure).exp_m1());
    let heat = 0.5 * reduced * slip.dot(slip) * -(-2.0 * exposure).exp_m1();
    for axis in 0..3 {
        assert!((a(impulse)[axis] - a(expected)[axis]).abs() < 4e-21);
        assert!((mass * (final_velocity[axis] - initial[axis]) - a(impulse)[axis]).abs() < 2e-20);
    }
    let body_change = after
        .linear_momentum_world()
        .sub(before.linear_momentum_world());
    assert!(body_change.add(impulse).norm_squared().sqrt() < 2e-20);
    assert!((loss - heat).abs() < 2e-21);
    let kinetic_change = after.kinetic_energy(properties).unwrap()
        - before.kinetic_energy(properties).unwrap()
        + 0.5 * mass * (v(final_velocity).dot(v(final_velocity)) - v(initial).dot(v(initial)));
    assert!((kinetic_change + loss).abs() < 2e-20);
    assert_eq!(after.pose(), before.pose());
}

#[test]
fn evolved_rheological_liquid_recoils_without_detaching_source_history() {
    with_liquids(|flow, _, gate| {
        for directional in [false, true] {
            let mut liquid = flow.clone();
            let before = liquid.clone();
            let (mut rigid, properties) = body_at_cell();
            let original = rigid;
            let [x, y, z] = CELL;
            let sample = liquid.sample_at_si(x, y, z).unwrap().unwrap();
            let units = liquid.momentum_scaling();
            let mass =
                liquid.lattice().region_mass(|x, y, z| [x, y, z] == CELL) * units.mass_unit_kg();
            let unit_bits = [
                units.cell_size_m(),
                units.mass_unit_kg(),
                units.velocity_si(1.0).unwrap(),
            ]
            .map(f64::to_bits);
            assert!(mass > 0.0);
            let (impulse, loss) = drag(&mut liquid, &mut rigid, properties, directional, gate);
            let after = liquid.sample_at_si(x, y, z).unwrap().unwrap();
            check_recoil(
                original,
                rigid,
                properties,
                mass,
                sample.velocity_m_s,
                after.velocity_m_s,
                impulse,
                loss,
                directional,
            );
            let scaling = liquid.momentum_scaling();
            assert_eq!(
                unit_bits,
                [
                    scaling.cell_size_m(),
                    scaling.mass_unit_kg(),
                    scaling.velocity_si(1.0).unwrap()
                ]
                .map(f64::to_bits)
            );
            assert!(core::ptr::eq(liquid.law(), before.law()));
            assert_eq!(
                liquid.lattice().grid().steps(),
                before.lattice().grid().steps()
            );
            assert_eq!(
                liquid.lattice().geometry_epoch(),
                before.lattice().geometry_epoch()
            );
            assert_eq!(
                liquid.ledger_mass_kg().unwrap().to_bits(),
                before.ledger_mass_kg().unwrap().to_bits()
            );
            for z in 0..4 {
                for y in 0..4 {
                    for x in 0..4 {
                        assert_eq!(liquid.material_at([x, y, z]), before.material_at([x, y, z]));
                        assert_eq!(liquid.rheology_at([x, y, z]), before.rheology_at([x, y, z]));
                    }
                }
            }
            liquid.step_with_gate(gate).unwrap();
            assert_eq!(
                liquid.lattice().grid().steps(),
                before.lattice().grid().steps() + 1
            );
        }
    });
}

#[test]
fn thermal_rheology_recoils_with_unchanged_enthalpy_and_accepted_time() {
    with_liquids(|_, flow, gate| {
        for directional in [false, true] {
            let mut liquid = flow.clone();
            let before = liquid.clone();
            let (mut rigid, properties) = body_at_cell();
            let original = rigid;
            let [x, y, z] = CELL;
            let sample = liquid.flow().sample_at_si(x, y, z).unwrap().unwrap();
            let mass = liquid
                .flow()
                .lattice()
                .region_mass(|x, y, z| [x, y, z] == CELL)
                * liquid.momentum_scaling().mass_unit_kg();
            let (impulse, loss) = drag(&mut liquid, &mut rigid, properties, directional, gate);
            let after = liquid.flow().sample_at_si(x, y, z).unwrap().unwrap();
            check_recoil(
                original,
                rigid,
                properties,
                mass,
                sample.velocity_m_s,
                after.velocity_m_s,
                impulse,
                loss,
                directional,
            );
            assert_eq!(liquid.accepted_steps(), before.accepted_steps());
            assert_eq!(
                liquid.physical_time_s().to_bits(),
                before.physical_time_s().to_bits()
            );
            assert_eq!(
                liquid.total_enthalpy_j().unwrap().to_bits(),
                before.total_enthalpy_j().unwrap().to_bits()
            );
            assert_eq!(
                liquid.flow().ledger_mass_kg().unwrap().to_bits(),
                before.flow().ledger_mass_kg().unwrap().to_bits()
            );
            assert_eq!(
                liquid.flow().lattice().geometry_epoch(),
                before.flow().lattice().geometry_epoch()
            );
            for z in 0..4 {
                for y in 0..4 {
                    for x in 0..4 {
                        let p = [x, y, z];
                        assert_eq!(
                            liquid.specific_enthalpy_at(p),
                            before.specific_enthalpy_at(p)
                        );
                        assert_eq!(liquid.material_at(p), before.material_at(p));
                        assert_eq!(liquid.rheology_at(p), before.rheology_at(p));
                    }
                }
            }
            liquid.step(gate).unwrap();
            assert_eq!(liquid.accepted_steps(), before.accepted_steps() + 1);
            assert!((liquid.physical_time_s() - before.physical_time_s() - DT).abs() < 1e-15);
        }
    });
}

#[test]
fn thermal_native_endpoint_speed_refusal_preserves_both_owners_and_retries() {
    with_liquids(|_, liquid, gate| {
        let [x, y, z] = CELL;
        let sample = liquid.flow().sample_at_si(x, y, z).unwrap().unwrap();
        let target = [0.0, 0.04, 0.0];
        let mut field = vec![[0.0; 3]; 64];
        field[INDEX] = core::array::from_fn(|a| {
            sample.density_kg_m3 * sample.fill_fraction * (target[a] - sample.velocity_m_s[a])
        });
        liquid
            .apply_impulse_density(
                &field,
                MaterialFreeSurfaceImpulseControls3 {
                    max_cells: 64,
                    maximum_speed_m_s: 0.08,
                    maximum_momentum_defect_kg_m_s: 1e-18,
                    maximum_work_defect_j: 1e-19,
                },
                gate,
            )
            .unwrap();
        let (_, properties) = body_at_cell();
        let mut rigid = RigidBodyState::new(
            body_at_cell().0.pose(),
            Vec3::new(0.04 * properties.mass(), 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap();
        let original = rigid;
        let before = liquid.clone();
        let law = heterogeneous::PrincipalLiquidDrag3 {
            rates_per_s: [1e6, 0.0, 0.0],
            axes_to_world: UnitQuaternion::IDENTITY,
        };
        let mut limits = heterogeneous::HeterogeneousRigidDragControls {
            exchange: LiquidRigidDragControls {
                maximum_speed_m_s: 0.05,
                ..controls()
            },
            substeps: 1,
            maximum_pair_exchanges: 1,
        };
        // Both initial speeds are 0.04 m/s. Directional exchange retains the
        // liquid's y speed while adding x momentum, exceeding 0.05 only at the
        // staged endpoint. The body candidate has already been formed.
        assert!(
            apply_heterogeneous_rigid_drag(
                liquid,
                &mut rigid,
                properties,
                &[INDEX],
                &[law],
                Vec3::ZERO,
                0.2,
                limits,
                gate,
            )
            .is_err()
        );
        assert_eq!(rigid, original);
        assert_eq!(
            liquid.flow().lattice().grid().state_bits(),
            before.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            liquid.total_enthalpy_j().unwrap().to_bits(),
            before.total_enthalpy_j().unwrap().to_bits()
        );
        assert_eq!(liquid.accepted_steps(), before.accepted_steps());
        assert_eq!(liquid.material_at(CELL), before.material_at(CELL));
        assert_eq!(liquid.rheology_at(CELL), before.rheology_at(CELL));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        limits.exchange.maximum_speed_m_s = 0.08;
        assert!(
            apply_heterogeneous_rigid_drag(
                liquid,
                &mut rigid,
                properties,
                &[INDEX],
                &[law],
                Vec3::ZERO,
                0.2,
                limits,
                &cancelled,
            )
            .is_err()
        );
        assert_eq!(rigid, original);
        assert_eq!(
            liquid.flow().lattice().grid().state_bits(),
            before.flow().lattice().grid().state_bits()
        );
        let mut direct = before.clone();
        let mut direct_body = original;
        let retry = apply_heterogeneous_rigid_drag(
            liquid,
            &mut rigid,
            properties,
            &[INDEX],
            &[law],
            Vec3::ZERO,
            0.2,
            limits,
            gate,
        )
        .unwrap();
        let expected = apply_heterogeneous_rigid_drag(
            &mut direct,
            &mut direct_body,
            properties,
            &[INDEX],
            &[law],
            Vec3::ZERO,
            0.2,
            limits,
            gate,
        )
        .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(rigid, direct_body);
        assert_eq!(
            liquid.flow().lattice().grid().state_bits(),
            direct.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            liquid.total_enthalpy_j().unwrap().to_bits(),
            before.total_enthalpy_j().unwrap().to_bits()
        );
    });
}
