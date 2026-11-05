//! Thermal fracture-source changes preserve an accepted body's physical motion.
//! Synthetic source curves exercise the native model; no measured calibration
//! or physical heat/expansion law is claimed by this fixture.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::dynamics::resolve_brittle_dynamic_material;
use fs_material::brittle::source::{
    BRITTLE_BK_EXPONENT_PROPERTY, BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
    BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY, BRITTLE_SHEAR_STRENGTH_PROPERTY,
    BRITTLE_TENSILE_STRENGTH_PROPERTY, ResolvedBrittleCrackBand,
};
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::TetAssemblyBudget;
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::dynamics::temperature::BrittleTetDynamicTemperatureError;
use fs_solid::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use fs_solid::brittle::solid::temperature::BrittleTetTemperatureError;
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetError, BrittleTetSolid};
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::IncrementSettings;

const YOUNG: f64 = 30e9;
const STRENGTH: f64 = 3e6;
const DENSITY: f64 = 2500.0;
const NODES: [[f64; 3]; 4] = [
    [0.0; 3],
    [0.03, 0.0, 0.0],
    [0.0, 0.02, 0.0],
    [0.0, 0.0, 0.01],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx<R>(cancelled: bool, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 918,
                kernel_id: 918,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn source_card(tag: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, hot_factor) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, YOUNG, 0.5),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0, 1.0),
        (
            BRITTLE_TENSILE_STRENGTH_PROPERTY,
            Pressure::DIMS,
            STRENGTH,
            0.25,
        ),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 4e6, 1.0),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
            1.0,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            200.0,
            1.0,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5, 1.0),
        // A hot density query changes, but the accepted reference inventory must not.
        (DENSITY_PROPERTY, Density::DIMS, DENSITY, 0.5),
    ] {
        let varied = hot_factor != 1.0;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varied {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (600.0, hot_factor * value)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: if varied {
                    ValidityDomain::unconstrained().with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.0,
                        600.0,
                    )
                } else {
                    ValidityDomain::unconstrained()
                },
                interpolation: if varied {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: format!("synthetic brittle thermal motion {tag}"),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic brittle thermal motion".into(),
            phase: "solid".into(),
            process: tag.into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn material(source: &MaterialCard, geometry: &TetLoadGeometry3<'_>) -> ResolvedBrittleCrackBand {
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    resolve_brittle_dynamic_material(
        source,
        &QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.03, 1e-7),
    )
    .unwrap()
}

fn work_controls() -> BrittleTetControls {
    BrittleTetControls {
        material_energy_tolerance_j: 1e-9,
        equilibrium_work_tolerance_j: 1e-9,
        maximum_endpoint_work_defect_j: 1.0,
    }
}

fn midpoint_controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-9,
        equilibrium_work_tolerance_j: 1e-9,
        kinetic_work_tolerance_j: 1e-9,
        total_energy_tolerance_j: 1e-9,
        difference_step_strain: Some(1e-8),
    }
}

fn preloaded<'a>(
    geometry: &'a TetLoadGeometry3<'a>,
    materials: &'a [ResolvedBrittleCrackBand],
    strain: f64,
    cx: &Cx<'_>,
) -> BrittleTetSolid<'a> {
    let prescribed: Vec<_> = (0..12)
        .filter(|&dof| dof != 7 && dof != 11)
        .map(|dof| (dof, if dof == 3 { 0.03 * strain } else { 0.0 }))
        .collect();
    BrittleTetSolid::try_new(geometry, materials, TetAssemblyBudget::standard(), cx)
        .unwrap()
        .advance(
            &prescribed,
            &[],
            IncrementSettings::default(),
            work_controls(),
            cx,
        )
        .unwrap()
        .solid
}

#[test]
fn thermal_midpoint_balances_parameter_fracture_and_motion_through_heating_and_cooling() {
    with_cx(false, |cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = source_card("history");
        let cold = [material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let solid = preloaded(&geometry, &cold, 0.75 * STRENGTH / YOUNG, cx);
        let start =
            BrittleTetMotion::try_from_solid(&solid, &[[0.2, -0.1, 0.05]; 4], 1.25, 4, cx).unwrap();
        let motion = start
            .advance_midpoint(
                1e-8,
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let snapshot = motion.solid().solution().clone();
        let old_loss = motion.solid().fracture_dissipation_j();
        let old_damage = snapshot.element_states[0].damage();
        assert_eq!(old_damage, 0.0);
        let dt = 1e-7;
        let heated = motion
            .advance_midpoint_at_temperature(
                dt,
                &hot,
                "T",
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        let next = &heated.mechanical.motion;
        assert_eq!(next.time_s(), motion.time_s() + dt);
        assert_eq!(next.accepted_steps(), 2);
        assert_eq!(next.remaining_steps(), 2);
        close(
            heated.mechanical.energy.inertia.kinetic_before_j,
            motion.kinetic_energy_j(),
            1e-15,
        );
        for node in 0..4 {
            for axis in 0..3 {
                close(
                    next.solid().solution().displacement_m[node][axis]
                        - snapshot.displacement_m[node][axis],
                    0.5 * dt
                        * (next.velocity_mps()[node][axis] + motion.velocity_mps()[node][axis]),
                    1e-16,
                );
            }
        }
        assert_eq!(next.reference_mass_kg(), motion.reference_mass_kg());
        assert_eq!(next.reference_density_kg_m3(), &[DENSITY]);
        assert!(std::ptr::eq(
            next.reference_density_kg_m3().as_ptr(),
            motion.reference_density_kg_m3().as_ptr()
        ));
        assert_eq!(hot[0].reference_density_kg_m3().unwrap(), 0.5 * DENSITY);
        assert!(std::ptr::eq(
            next.solid().materials()[0].law(),
            hot[0].law()
        ));
        let transferred = cold[0]
            .transfer_temperature(&snapshot.element_states[0], &hot[0], "T", 1e-7)
            .unwrap();
        let expected = hot[0]
            .law()
            .trial_with_mean_stress(
                &next.solid().solution().element_strain[0],
                &transferred.state,
            )
            .unwrap();
        assert_eq!(
            next.solid().solution().element_states[0],
            expected.endpoint.state
        );
        assert_eq!(
            heated.mechanical.element_energy[0],
            expected.endpoint.energy
        );
        assert!(transferred.state.damage() > old_damage);
        assert!(heated.transfer.fracture_dissipation_j > 0.0);
        assert_eq!(heated.mechanical.energy.inertia.numerical_loss_j, 0.0);
        assert_eq!(
            heated.mechanical.energy.stored_before_j,
            heated.transfer.stored_after_j
        );
        close(
            heated.transfer.parameter_work_j + heated.mechanical.energy.external_mean_work_j(),
            next.solid().stored_energy_j() - motion.solid().stored_energy_j()
                + next.solid().fracture_dissipation_j()
                - old_loss
                + heated.mechanical.energy.inertia.kinetic_change_j,
            1e-9,
        );
        close(heated.transfer.balance_error_j, 0.0, 1e-9);
        close(heated.total_balance_error_j, 0.0, 1e-9);
        let cooled = next
            .advance_midpoint_at_temperature(
                dt,
                &cold,
                "T",
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        let final_motion = &cooled.mechanical.motion;
        assert!(
            final_motion.solid().solution().element_states[0].damage()
                >= next.solid().solution().element_states[0].damage()
        );
        assert!(
            final_motion.solid().fracture_dissipation_j() >= next.solid().fracture_dissipation_j()
        );
        assert_eq!(final_motion.time_s(), next.time_s() + dt);
        assert_eq!(final_motion.accepted_steps(), 3);
        assert_eq!(final_motion.remaining_steps(), 1);
        assert_eq!(final_motion.reference_mass_kg(), motion.reference_mass_kg());
        close(
            heated.transfer.parameter_work_j
                + cooled.transfer.parameter_work_j
                + heated.mechanical.energy.external_mean_work_j()
                + cooled.mechanical.energy.external_mean_work_j(),
            final_motion.solid().stored_energy_j() - motion.solid().stored_energy_j()
                + final_motion.solid().fracture_dissipation_j()
                - old_loss
                + heated.mechanical.energy.inertia.kinetic_change_j
                + cooled.mechanical.energy.inertia.kinetic_change_j,
            2e-9,
        );
        assert_eq!(motion.solid().solution(), &snapshot);
        assert_eq!(motion.solid().fracture_dissipation_j(), old_loss);
    });
}

#[test]
fn hot_midpoint_uses_changed_tractions_and_frozen_mass_without_resetting_step_budget() {
    with_cx(false, |cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = source_card("acceleration");
        let cold = [material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let solid = preloaded(&geometry, &cold, 0.25 * STRENGTH / YOUNG, cx);
        let motion = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 1, cx).unwrap();
        let dt = 1e-7;
        let cold_step = motion
            .advance_midpoint(
                dt,
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        let heated = motion
            .advance_midpoint_at_temperature(
                dt,
                &hot,
                "T",
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        let hot_step = &heated.mechanical;
        assert!(hot_step.motion.velocity_mps()[1][0] < 0.0);
        assert!(
            hot_step.motion.velocity_mps()[1][0].abs()
                < cold_step.motion.velocity_mps()[1][0].abs()
        );
        close(hot_step.energy.inertia.kinetic_before_j, 0.0, 1e-18);
        assert!(hot_step.motion.kinetic_energy_j() > 0.0);
        assert_eq!(hot_step.energy.inertia.numerical_loss_j, 0.0);
        assert_eq!(hot[0].reference_density_kg_m3().unwrap(), 0.5 * DENSITY);
        assert!(heated.transfer.parameter_work_j < 0.0);
        close(
            hot_step.motion.kinetic_energy_j() + hot_step.motion.solid().stored_energy_j(),
            motion.solid().stored_energy_j() + heated.transfer.parameter_work_j,
            1e-9,
        );
        for accepted in [&hot_step.motion, &cold_step.motion] {
            assert_eq!(accepted.reference_mass_kg(), motion.reference_mass_kg());
            assert_eq!(accepted.reference_density_kg_m3(), &[DENSITY]);
            assert_eq!(accepted.time_s(), dt);
            assert_eq!(accepted.accepted_steps(), 1);
            assert_eq!(accepted.remaining_steps(), 0);
            assert!(matches!(
                accepted.advance_midpoint_at_temperature(
                    dt,
                    &hot,
                    "T",
                    &[],
                    &[],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    cx
                ),
                Err(BrittleTetDynamicTemperatureError::Dynamic(
                    BrittleTetDynamicError::StepBudget {
                        accepted_steps: 1,
                        maximum_steps: 1
                    }
                ))
            ));
        }
    });
}

#[test]
fn thermal_midpoint_source_domain_solve_and_cancellation_refusals_preserve_retryability() {
    with_cx(false, |cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = source_card("rollback");
        let cold = [material(&source, &geometry)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let other = source_card("foreign");
        let foreign = [material(&other, &geometry)
            .at_temperature(&other, "T", 600.0)
            .unwrap()];
        let solid = preloaded(&geometry, &cold, 0.002, cx);
        let motion =
            BrittleTetMotion::try_from_solid(&solid, &[[0.1, 0.0, 0.0]; 4], 2.0, 3, cx).unwrap();
        let snapshot = motion.solid().solution().clone();
        let storage = motion.solid().stored_energy_j();
        let loss = motion.solid().fracture_dissipation_j();
        assert!(matches!(
            motion.advance_midpoint_at_temperature(
                1e-7,
                &foreign,
                "T",
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx
            ),
            Err(BrittleTetDynamicTemperatureError::Transfer(
                BrittleTetTemperatureError::Source { element: 0, .. }
            ))
        ));
        assert!(matches!(
            motion.advance_midpoint_at_temperature(
                1e-7,
                &[],
                "T",
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx
            ),
            Err(BrittleTetDynamicTemperatureError::Transfer(
                BrittleTetTemperatureError::Solid(BrittleTetError::MaterialCount { .. })
            ))
        ));
        for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
            let controls = BrittleTetMidpointControls {
                total_energy_tolerance_j: tolerance,
                ..midpoint_controls()
            };
            assert!(
                motion
                    .advance_midpoint_at_temperature(
                        1e-7,
                        &hot,
                        "T",
                        &[],
                        &[],
                        IncrementSettings::default(),
                        controls,
                        cx
                    )
                    .is_err()
            );
        }
        assert!(
            motion
                .advance_midpoint_at_temperature(
                    1e-7,
                    &hot,
                    "not-temperature",
                    &[],
                    &[],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    cx
                )
                .is_err()
        );
        // A valid thermal transfer followed by an inadmissible mechanical
        // endpoint must not publish the private heated history or source.
        assert!(
            motion
                .advance_midpoint_at_temperature(
                    1e-7,
                    &hot,
                    "T",
                    &[(3, 0.1)],
                    &[],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    cx
                )
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(matches!(
                motion.advance_midpoint_at_temperature(
                    1e-7,
                    &hot,
                    "T",
                    &[],
                    &[],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    cancelled
                ),
                Err(BrittleTetDynamicTemperatureError::Dynamic(
                    BrittleTetDynamicError::Solid(BrittleTetError::Cancelled)
                ))
            ));
        });
        assert_eq!(motion.solid().solution(), &snapshot);
        assert_eq!(motion.solid().stored_energy_j(), storage);
        assert_eq!(motion.solid().fracture_dissipation_j(), loss);
        assert_eq!(motion.velocity_mps(), &[[0.1, 0.0, 0.0]; 4]);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
        assert_eq!(motion.remaining_steps(), 3);
        let resumed = motion
            .advance_midpoint_at_temperature(
                1e-7,
                &hot,
                "T",
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        let repeated = motion
            .advance_midpoint_at_temperature(
                1e-7,
                &hot,
                "T",
                &[],
                &[],
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            resumed.mechanical.motion.solid().solution(),
            repeated.mechanical.motion.solid().solution()
        );
        assert_eq!(resumed.element_transfer, repeated.element_transfer);
        assert_eq!(resumed.transfer, repeated.transfer);
        assert_eq!(resumed.mechanical.energy, repeated.mechanical.energy);
        assert_eq!(
            resumed.total_balance_error_j,
            repeated.total_balance_error_j
        );
    });
}
