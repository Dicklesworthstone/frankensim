//! Source-bound brittle equilibrium and energy through the original P1 solver.
//! Synthetic coefficients verify the stated model; they do not validate a real
//! material, an automatically selected crack plane or a converged crack path.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::{
    BRITTLE_BK_EXPONENT_PROPERTY, BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
    BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY, BRITTLE_SHEAR_STRENGTH_PROPERTY,
    BRITTLE_TENSILE_STRENGTH_PROPERTY, ResolvedBrittleCrackBand, resolve_brittle_crack_band,
};
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::TetAssemblyBudget;
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::temperature::BrittleTetTemperatureError;
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetError, BrittleTetSolid};
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

const YOUNG: f64 = 30.0e9;
const POISSON: f64 = 0.22;
const STRENGTH: f64 = 3.0e6;
const FRACTURE_WORK: f64 = 100.0;
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

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()).max(1.0e-14),
        "{actual:e} != {expected:e}"
    );
}

fn points(length: f64) -> [[f64; 3]; 4] {
    [
        [0.0; 3],
        [length, 0.0, 0.0],
        [0.0, 0.02, 0.0],
        [0.0, 0.0, 0.01],
    ]
}

fn card(strength: f64, tag: &str) -> MaterialCard {
    source_card(strength, tag, false)
}

fn source_card(strength: f64, tag: &str, thermal: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, YOUNG),
        (
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            if thermal { 0.0 } else { POISSON },
        ),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, strength),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 4.0e6),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            FRACTURE_WORK,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            200.0,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5),
    ] {
        let varied = thermal
            && (name == YOUNG_MODULUS_PROPERTY || name == BRITTLE_TENSILE_STRENGTH_PROPERTY);
        let hot_factor = if name == YOUNG_MODULUS_PROPERTY {
            0.5
        } else {
            0.25
        };
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
                // Typed like every query of this card: an untyped "T" axis
                // cannot answer a query that types it in kelvin.
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
                    source: format!("synthetic brittle coupon {tag}"),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle solid".into(),
            phase: "solid".into(),
            process: tag.into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn material(
    source: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
    element: usize,
) -> ResolvedBrittleCrackBand {
    let band = TetCrackBandGeometry::from_element(geometry, element, [1.0, 0.0, 0.0]).unwrap();
    resolve_brittle_crack_band(
        source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.03, 1e-7),
    )
    .unwrap()
}

fn controls() -> BrittleTetControls {
    BrittleTetControls {
        material_energy_tolerance_j: 1e-9,
        equilibrium_work_tolerance_j: 1e-9,
        maximum_endpoint_work_defect_j: 1.0,
    }
}

/// Two transverse strain coordinates remain genuinely free and find equilibrium.
fn prescribed(length: f64, strain: f64) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 7 && dof != 11)
        .map(|dof| (dof, if dof == 3 { length * strain } else { 0.0 }))
        .collect()
}

#[test]
fn source_bound_solid_softens_unloads_fails_and_recloses_with_physical_energy() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "opening");
        let materials = [material(&source, &geometry, 0)];
        let mut solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        assert!(solid.solution().increments.is_empty());
        assert!(std::ptr::eq(solid.materials()[0].law(), materials[0].law()));
        let band = solid.band_geometry()[0].clone();
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * band.band_width_m());
        let peak = 0.002;
        let peak_stress = STRENGTH * (failure - peak) / (failure - onset);
        let path = [
            0.5 * onset,
            peak,
            0.25 * peak,
            peak,
            1.2 * failure,
            -0.001,
            0.001,
        ];
        let expected = [
            0.5 * STRENGTH,
            peak_stress,
            0.25 * peak_stress,
            peak_stress,
            0.0,
            -0.001 * YOUNG,
            0.0,
        ];
        let mut total_work = 0.0;
        let mut previous_damage = 0.0;
        let mut previous_loss = 0.0;
        let mut peak_damage = 0.0;
        let mut peak_loss = 0.0;
        for (index, (&strain, &stress)) in path.iter().zip(&expected).enumerate() {
            let targets = prescribed(0.03, strain);
            let step = solid
                .advance(&targets, &[], IncrementSettings::default(), controls(), cx)
                .unwrap();
            let solution = step.solid.solution();
            let state = &solution.element_states[0];
            assert!((solution.element_stress[0][0] - stress).abs() < 1e-7 * STRENGTH);
            for value in &solution.element_stress[0][1..] {
                assert!(value.abs() < 1e-7 * STRENGTH);
            }
            close(solution.element_strain[0][0], strain, 2e-12);
            assert!((solution.element_strain[0][1] + POISSON * stress / YOUNG).abs() < 2e-10);
            assert!((solution.element_strain[0][2] + POISSON * stress / YOUNG).abs() < 2e-10);
            let axial = targets.iter().position(|(dof, _)| *dof == 3).unwrap();
            assert!(
                (solution.increments[0].reactions_n[axial]
                    - stress * band.reference_volume_m3() / 0.03)
                    .abs()
                    < 1e-5
            );
            assert!(state.damage() >= previous_damage);
            assert!(step.solid.fracture_dissipation_j() >= previous_loss);
            assert!(step.energy.fracture_dissipation_j >= 0.0);
            assert!(step.energy.material_balance_error_j.abs() < 1e-9);
            assert!(step.energy.equilibrium_work_residual_j.abs() < 1e-9);
            close(
                step.energy.fracture_dissipation_j,
                band.reference_volume_m3() * step.element_energy[0].dissipated,
                2e-12,
            );
            total_work += step.energy.material_work_j;
            assert!(
                (total_work - step.solid.stored_energy_j() - step.solid.fracture_dissipation_j())
                    .abs()
                    < 1e-9
            );
            if index == 1 {
                peak_damage = state.damage();
                peak_loss = step.solid.fracture_dissipation_j();
                assert!(peak_damage > 0.0 && peak_damage < 1.0);
            }
            if index == 2 || index == 3 {
                close(state.damage(), peak_damage, 2e-12);
                close(step.solid.fracture_dissipation_j(), peak_loss, 2e-12);
            }
            if index >= 4 {
                assert!(state.failed());
                close(
                    step.solid.fracture_dissipation_j(),
                    FRACTURE_WORK * band.projected_area_m2(),
                    2e-9,
                );
            }
            if index == 4 {
                // A failed endpoint has zero tensile force but positive exact
                // path work. This negative quadrature defect is not healing.
                assert!(step.energy.endpoint_work_defect_j < 0.0);
                assert!(step.energy.fracture_dissipation_j > 0.0);
            }
            assert_eq!(
                step.solid.materials()[0].material().card_identity(),
                source.content_hash()
            );
            previous_damage = state.damage();
            previous_loss = step.solid.fracture_dissipation_j();
            solid = step.solid;
        }
    });
}

#[test]
fn independent_sourced_strengths_and_widths_retain_local_history_and_surface_work() {
    with_cx(false, |cx| {
        let nodes: Vec<_> = points(0.03)
            .into_iter()
            .chain(points(0.015).map(|p| [p[0] + 0.1, p[1], p[2]]))
            .collect();
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let geometry =
            TetLoadGeometry3::new(&nodes, &tets, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let sources = [card(STRENGTH, "weaker"), card(2.0 * STRENGTH, "stronger")];
        let materials = [
            material(&sources[0], &geometry, 0),
            material(&sources[1], &geometry, 1),
        ];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let targets = |strain| {
            prescribed(0.03, strain)
                .into_iter()
                .chain(
                    prescribed(0.015, strain)
                        .into_iter()
                        .map(|(d, v)| (d + 12, v)),
                )
                .collect::<Vec<_>>()
        };
        let partial = solid
            .advance(
                &targets(0.00015),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert!(partial.solid.solution().element_states[0].damage() > 0.0);
        assert_eq!(partial.solid.solution().element_states[1].damage(), 0.0);
        assert!(
            partial.solid.solution().element_stress[1][0]
                > partial.solid.solution().element_stress[0][0]
        );
        let failed = partial
            .solid
            .advance(
                &targets(0.009),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let mut expected_loss = 0.0;
        for element in 0..2 {
            let band = &failed.solid.band_geometry()[element];
            assert!(failed.solid.solution().element_states[element].failed());
            close(
                failed.solid.solution().element_states[element].dissipated_j_m3()
                    * band.reference_volume_m3(),
                FRACTURE_WORK * band.projected_area_m2(),
                2e-9,
            );
            assert_eq!(
                failed.solid.materials()[element].material().card_identity(),
                sources[element].content_hash()
            );
            expected_loss += FRACTURE_WORK * band.projected_area_m2();
        }
        close(failed.solid.fracture_dissipation_j(), expected_loss, 2e-9);
        close(
            failed.solid.band_geometry()[0].band_width_m(),
            2.0 * failed.solid.band_geometry()[1].band_width_m(),
            2e-12,
        );
        let swapped = [materials[1].clone(), materials[0].clone()];
        assert!(matches!(
            BrittleTetSolid::try_new(&geometry, &swapped, TetAssemblyBudget::standard(), cx),
            Err(BrittleTetError::WidthMismatch { element: 0, .. })
        ));
    });
}

#[test]
fn force_control_and_applied_loads_on_moving_supports_keep_the_work_partition() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "work partition");
        let materials = [material(&source, &geometry, 0)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let strain = 0.5 * STRENGTH / YOUNG;
        let volume = solid.band_geometry()[0].reference_volume_m3();
        let force = YOUNG * strain * volume / 0.03;
        let fixed: Vec<_> = prescribed(0.03, 0.0)
            .into_iter()
            .filter(|(d, _)| *d != 3)
            .collect();
        let loaded = solid
            .advance(
                &fixed,
                &[(3, force)],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        close(
            loaded.solid.solution().displacement_m[1][0],
            0.03 * strain,
            2e-10,
        );
        close(
            loaded.solid.solution().element_stress[0][0],
            YOUNG * strain,
            2e-10,
        );
        let storage = 0.5 * YOUNG * strain * strain * volume;
        close(loaded.solid.stored_energy_j(), storage, 2e-10);
        assert_eq!(loaded.solid.fracture_dissipation_j(), 0.0);
        close(loaded.energy.material_work_j, storage, 2e-10);
        close(
            loaded.energy.dead_load_endpoint_work_j,
            2.0 * storage,
            2e-10,
        );
        assert_eq!(loaded.energy.support_endpoint_work_j, 0.0);
        close(loaded.energy.endpoint_work_defect_j, storage, 2e-10);
        let imposed = solid
            .advance(
                &prescribed(0.03, strain),
                &[(3, 0.4 * force)],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        close(
            imposed.energy.dead_load_endpoint_work_j,
            0.8 * storage,
            2e-10,
        );
        close(imposed.energy.support_endpoint_work_j, 1.2 * storage, 2e-10);
        close(
            imposed.energy.external_endpoint_work_defect_j(),
            storage,
            2e-10,
        );
        close(
            imposed.solid.stored_energy_j(),
            loaded.solid.stored_energy_j(),
            2e-10,
        );
    });
}

#[test]
fn geometry_domain_newton_energy_and_cancellation_refusals_preserve_accepted_solid() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "rollback");
        let materials = [material(&source, &geometry, 0)];
        assert!(matches!(
            BrittleTetSolid::try_new(&geometry, &[], TetAssemblyBudget::standard(), cx),
            Err(BrittleTetError::MaterialCount { .. })
        ));
        let mut wrong = materials[0].law().settings();
        wrong.normal = [0.0, 1.0, 0.0];
        let wrong_plane =
            [
                ResolvedBrittleCrackBand::try_from_material(materials[0].material().clone(), wrong)
                    .unwrap(),
            ];
        assert!(matches!(
            BrittleTetSolid::try_new(&geometry, &wrong_plane, TetAssemblyBudget::standard(), cx),
            Err(BrittleTetError::WidthMismatch { .. })
        ));
        let virgin =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let accepted = virgin
            .advance(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let snapshot = (*accepted.solution()).clone();
        let stored = accepted.stored_energy_j();
        let loss = accepted.fracture_dissipation_j();
        assert!(
            accepted
                .advance(
                    &prescribed(0.03, 0.08),
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        let targets = prescribed(0.03, 0.001);
        let stalled = IncrementSettings {
            max_newton_iterations: 1,
            ..IncrementSettings::default()
        };
        assert!(matches!(
            accepted.advance(&targets, &[], stalled, controls(), cx),
            Err(BrittleTetError::Solve(
                SmallStrainTetError::NewtonStalled { .. }
            ))
        ));
        let exact_endpoint = BrittleTetControls {
            maximum_endpoint_work_defect_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            accepted.advance(
                &targets,
                &[],
                IncrementSettings::default(),
                exact_endpoint,
                cx
            ),
            Err(BrittleTetError::EnergyGate { .. })
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                accepted.advance(
                    &targets,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    cancelled
                ),
                Err(BrittleTetError::Cancelled)
            ));
        });
        assert_eq!(accepted.solution(), &snapshot);
        assert_eq!(accepted.stored_energy_j(), stored);
        assert_eq!(accepted.fracture_dissipation_j(), loss);
        let resumed = accepted
            .advance(&targets, &[], IncrementSettings::default(), controls(), cx)
            .unwrap();
        let repeated = accepted
            .advance(&targets, &[], IncrementSettings::default(), controls(), cx)
            .unwrap();
        assert_eq!(resumed.solid.solution(), repeated.solid.solution());
        assert_eq!(resumed.energy, repeated.energy);
        close(
            resumed.solid.solution().element_stress[0][0],
            0.5 * snapshot.element_stress[0][0],
            2e-9,
        );
        close(resumed.solid.fracture_dissipation_j(), loss, 2e-12);
    });
}

fn thermal_material(
    source: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
    temperature: f64,
) -> ResolvedBrittleCrackBand {
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    resolve_brittle_crack_band(
        source,
        &QueryPoint::new()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                temperature,
            )
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.03, 1e-7),
    )
    .unwrap()
}

#[test]
fn temperature_transfer_re_equilibrates_a_supported_load_and_can_initiate_held_strain_damage() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = source_card(STRENGTH, "thermal loading", true);
        let cold = [thermal_material(&source, &geometry, 300.0)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx).unwrap();
        let damaged = virgin
            .advance(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let unloaded = damaged
            .advance(
                &prescribed(0.03, 0.0002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let old_stress = unloaded.solution().element_stress[0][0];
        let old_damage = unloaded.solution().element_states[0].damage();
        let volume = unloaded.band_geometry()[0].reference_volume_m3();
        let force = old_stress * volume / 0.03;
        let fixed: Vec<_> = prescribed(0.03, 0.0)
            .into_iter()
            .filter(|(d, _)| *d != 3)
            .collect();
        let heated = unloaded
            .advance_at_temperature(
                &hot,
                "T",
                &fixed,
                &[(3, force)],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let final_solid = &heated.mechanical.solid;
        close(
            final_solid.solution().element_stress[0][0],
            old_stress,
            2e-9,
        );
        // At retained damage, halving E doubles the unloading compliance.
        // This load remains below the evolved hot envelope.
        close(final_solid.solution().element_strain[0][0], 0.0004, 2e-8);
        assert!(final_solid.solution().element_states[0].damage() >= old_damage);
        close(
            final_solid.solution().element_states[0].damage(),
            old_damage,
            2e-10,
        );
        assert!(final_solid.fracture_dissipation_j() >= unloaded.fracture_dissipation_j());
        assert!(heated.transfer.parameter_work_j < 0.0);
        assert!(heated.transfer.fracture_dissipation_j >= 0.0);
        assert!(heated.mechanical.energy.material_work_j > 0.0);
        assert!(heated.transfer.balance_error_j.abs() < 1e-9);
        assert!(heated.total_balance_error_j.abs() < 1e-9);
        close(
            final_solid.stored_energy_j(),
            2.0 * unloaded.stored_energy_j(),
            2e-8,
        );
        assert_eq!(
            final_solid.materials()[0].material().state_coordinate("T"),
            Some(600.0)
        );
        assert_eq!(
            unloaded.materials()[0].material().state_coordinate("T"),
            Some(300.0)
        );
        assert_eq!(
            final_solid.materials()[0].material().card_identity(),
            source.content_hash()
        );
        close(
            heated.transfer.parameter_work_j,
            volume * heated.element_transfer[0].parameter_work_j_m3,
            2e-12,
        );

        // An intact point below its original strength crosses the hot
        // strength at the SAME held strain. No imposed thermal strain is used.
        let held_strain = 0.75 * STRENGTH / YOUNG;
        let held = virgin
            .advance(
                &prescribed(0.03, held_strain),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        assert_eq!(held.solution().element_states[0].damage(), 0.0);
        let heated_held = held
            .advance_at_temperature(
                &hot,
                "T",
                &prescribed(0.03, held_strain),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let hot_strength = 0.25 * STRENGTH;
        let hot_young = 0.5 * YOUNG;
        let failure = 2.0 * FRACTURE_WORK
            / (hot_strength * heated_held.mechanical.solid.band_geometry()[0].band_width_m());
        let expected_stress =
            hot_strength * (failure - held_strain) / (failure - hot_strength / hot_young);
        assert!(heated_held.mechanical.solid.solution().element_states[0].damage() > 0.0);
        assert!(heated_held.transfer.fracture_dissipation_j > 0.0);
        close(
            heated_held.mechanical.solid.solution().element_stress[0][0],
            expected_stress,
            2e-10,
        );
        assert_eq!(
            heated_held.mechanical.solid.solution().element_strain,
            held.solution().element_strain
        );
        assert_eq!(heated_held.mechanical.energy.material_work_j, 0.0);
        assert_eq!(heated_held.mechanical.energy.fracture_dissipation_j, 0.0);
        assert!(heated_held.total_balance_error_j.abs() < 1e-9);
    });
}

#[test]
fn temperature_source_overload_and_late_work_refusals_keep_the_cold_solid_retryable() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = source_card(STRENGTH, "thermal rollback", true);
        let cold = [thermal_material(&source, &geometry, 300.0)];
        let hot = [cold[0].at_temperature(&source, "T", 600.0).unwrap()];
        let virgin =
            BrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx).unwrap();
        let damaged = virgin
            .advance(
                &prescribed(0.03, 0.002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let accepted = damaged
            .advance(
                &prescribed(0.03, 0.0002),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let snapshot = (*accepted.solution()).clone();
        let old_storage = accepted.stored_energy_j();
        let old_loss = accepted.fracture_dissipation_j();
        let volume = accepted.band_geometry()[0].reference_volume_m3();
        let force = accepted.solution().element_stress[0][0] * volume / 0.03;
        let fixed: Vec<_> = prescribed(0.03, 0.0)
            .into_iter()
            .filter(|(d, _)| *d != 3)
            .collect();
        let strict = BrittleTetControls {
            maximum_endpoint_work_defect_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            accepted.advance_at_temperature(
                &hot,
                "T",
                &fixed,
                &[(3, force)],
                IncrementSettings::default(),
                strict,
                cx
            ),
            Err(BrittleTetTemperatureError::Solid(
                BrittleTetError::EnergyGate {
                    account: "signed endpoint integration defect",
                    ..
                }
            ))
        ));
        // Positive axial stress cannot exceed the hot source peak; the native
        // static solver must refuse this unsupported load, not reset damage.
        let overload = 1.2 * (0.25 * STRENGTH) * volume / 0.03;
        let bounded = IncrementSettings {
            max_newton_iterations: 8,
            ..IncrementSettings::default()
        };
        assert!(
            accepted
                .advance_at_temperature(
                    &hot,
                    "T",
                    &fixed,
                    &[(3, overload)],
                    bounded,
                    controls(),
                    cx
                )
                .is_err()
        );
        let other = source_card(STRENGTH, "different source", true);
        let wrong = [thermal_material(&other, &geometry, 600.0)];
        assert!(matches!(
            accepted.advance_at_temperature(
                &wrong,
                "T",
                &fixed,
                &[(3, force)],
                IncrementSettings::default(),
                controls(),
                cx
            ),
            Err(BrittleTetTemperatureError::Source { element: 0, .. })
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                accepted.advance_at_temperature(
                    &hot,
                    "T",
                    &fixed,
                    &[(3, force)],
                    IncrementSettings::default(),
                    controls(),
                    cancelled
                ),
                Err(BrittleTetTemperatureError::Solid(
                    BrittleTetError::Cancelled
                ))
            ));
        });
        assert_eq!(accepted.solution(), &snapshot);
        assert_eq!(accepted.stored_energy_j(), old_storage);
        assert_eq!(accepted.fracture_dissipation_j(), old_loss);
        assert_eq!(
            accepted.materials()[0].material().state_coordinate("T"),
            Some(300.0)
        );
        let resumed = accepted
            .advance_at_temperature(
                &hot,
                "T",
                &fixed,
                &[(3, force)],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let repeated = accepted
            .advance_at_temperature(
                &hot,
                "T",
                &fixed,
                &[(3, force)],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            resumed.mechanical.solid.solution(),
            repeated.mechanical.solid.solution()
        );
        assert_eq!(resumed.transfer, repeated.transfer);
        assert_eq!(resumed.mechanical.energy, repeated.mechanical.energy);
    });
}

#[path = "brittle_tet_solid/path.rs"]
mod path;

#[path = "brittle_tet_solid/initiation.rs"]
mod initiation;

#[path = "brittle_tet_solid/dynamics.rs"]
mod dynamics;

#[path = "brittle_tet_solid/expansion.rs"]
mod expansion;

#[path = "brittle_tet_solid/continuation.rs"]
mod continuation;
