//! Native thermal strain must change component equilibrium and tensile damage.
//! Synthetic coefficients test the admitted fixed-plane, small-strain model;
//! transverse stresses here do not initiate additional unsupplied crack planes.

use super::*;
use fs_material::brittle::source::expansion::{ResolvedThermalBrittle, resolve_thermal_brittle};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
};
use fs_solid::brittle::solid::expansion::{
    ThermalBrittleTetEnergy, ThermalBrittleTetError, ThermalBrittleTetSolid,
};

#[path = "expansion/initiation.rs"]
mod initiation;

const RIGID_SUPPORTS: [(usize, f64); 6] =
    [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];

/// E, nu and the fracture envelope stay fixed so eigenstrain is independently
/// observable. The selected instantaneous alpha curve has a real knot at 300 K.
fn expansion_card(tag: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, YOUNG),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, STRENGTH),
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
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: format!("synthetic thermal brittle solid {tag}"),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, 8.0e-6), (300.0, 1.0e-5), (400.0, 2.0e-5)],
                dims: INVERSE_TEMPERATURE_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                200.0,
                400.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: format!("synthetic instantaneous expansion {tag}"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured thermal brittle solid".into(),
            phase: "solid".into(),
            process: tag.into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn temperature(value: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), value)
        .unwrap()
}

fn resolved_expansion(
    source: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
    reference_k: f64,
    current_k: f64,
) -> ResolvedThermalBrittle {
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    resolve_thermal_brittle(
        source,
        &temperature(reference_k),
        &temperature(current_k),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.03, 1e-7),
    )
    .unwrap()
}

fn fixed() -> Vec<(usize, f64)> {
    (0..12).map(|dof| (dof, 0.0)).collect()
}

fn affine_targets(nodes: &[[f64; 3]], diagonal: [f64; 3]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, xyz)| {
            (0..3).map(move |axis| (3 * node + axis, diagonal[axis] * xyz[axis]))
        })
        .collect()
}

fn check_balances(energy: &ThermalBrittleTetEnergy) {
    for error in [
        energy.strain_split_error_j,
        energy.material_balance_error_j,
        energy.equilibrium_work_residual_j,
        energy.total_balance_error_j,
    ] {
        assert!(error.abs() < 1e-9, "{energy:?}");
    }
}

#[test]
fn free_expansion_crosses_a_source_knot_without_artificial_tensile_damage() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = expansion_card("free");
        let reference = resolved_expansion(&source, &geometry, 300.0, 300.0);
        let temperatures = [
            300.0, 304.0, 308.0, 312.0, 308.0, 304.0, 300.0, 296.0, 292.0, 288.0,
        ];
        let materials: Vec<_> = temperatures
            .iter()
            .map(|&t| [reference.at_temperature(&source, t).unwrap()])
            .collect();
        let mut solid = ThermalBrittleTetSolid::try_new(
            &geometry,
            &materials[0],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let mut previous_free = 0.0;
        for (index, next) in materials.iter().enumerate().skip(1) {
            let delta = temperatures[index] - 300.0;
            // Independent exact integrals of the two linear alpha segments.
            let slope = if delta >= 0.0 { 1e-7 } else { 2e-8 };
            let free = 1e-5 * delta + 0.5 * slope * delta * delta;
            close(next[0].free_linear_strain(), free, 2e-12);
            assert!((free - previous_free).abs() < STRENGTH / YOUNG);
            let step = solid
                .advance(
                    next,
                    &RIGID_SUPPORTS,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            for (node, xyz) in nodes.iter().enumerate() {
                for axis in 0..3 {
                    assert!(
                        (step.solid.solution().displacement_m[node][axis] - free * xyz[axis]).abs()
                            < 1e-13
                    );
                }
            }
            for component in 0..6 {
                let expected = if component < 3 { free } else { 0.0 };
                assert!(
                    (step.solid.solution().element_strain[0][component] - expected).abs() < 5e-14
                );
                assert!(step.solid.solution().element_states[0].strain()[component].abs() < 5e-14);
                assert!(step.solid.solution().element_stress[0][component].abs() < 1e-4);
            }
            assert_eq!(step.solid.solution().element_states[0].damage(), 0.0);
            assert_eq!(step.solid.fracture_dissipation_j(), 0.0);
            assert!(step.solid.stored_energy_j().abs() < 1e-12);
            assert!(step.energy.mechanical_work_j.abs() < 1e-10);
            assert!(step.energy.eigenstrain_work_j.abs() < 1e-10);
            assert!(step.energy.constitutive_strain_work_j.abs() < 1e-10);
            check_balances(&step.energy);
            previous_free = free;
            solid = step.solid;
        }
        // Accumulated free strain exceeds tensile onset; only the mechanical
        // strain is an admissible damage driver.
        assert!(previous_free.abs() > STRENGTH / YOUNG);
        assert!(std::ptr::eq(
            solid.materials()[0].band().law(),
            materials.last().unwrap()[0].band().law(),
        ));
    });
}

#[test]
fn restrained_cooling_drives_native_tensile_softening_and_full_thermal_work() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = expansion_card("restrained");
        let reference = [resolved_expansion(&source, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&source, 300.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let step = solid
            .advance(
                &cold,
                &fixed(),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();

        let strain = 0.5 * (1e-5 + 2e-5) * 100.0;
        let volume = 0.03 * 0.02 * 0.01 / 6.0;
        let width = 0.03 / 3.0;
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * width);
        let stress = STRENGTH * (failure - strain) / (failure - onset);
        let opening_work = 0.5 * STRENGTH * onset + 0.5 * (STRENGTH + stress) * (strain - onset);
        let opening_storage = 0.5 * stress * strain;
        // The supplied x plane damages. The two intact in-plane directions
        // retain their complete elastic storage under isotropic restraint.
        let transverse_storage = YOUNG * strain * strain;
        let expected_work = volume * (opening_work + transverse_storage);
        let expected_storage = volume * (opening_storage + transverse_storage);
        let expected_loss = volume * (opening_work - opening_storage);

        close(cold[0].free_linear_strain(), -strain, 2e-12);
        close(
            step.solid.band_geometry()[0].reference_volume_m3(),
            volume,
            2e-12,
        );
        close(step.solid.band_geometry()[0].band_width_m(), width, 2e-12);
        let state = &step.solid.solution().element_states[0];
        assert!(state.damage() > 0.0 && state.damage() < 1.0);
        close(state.damage(), 1.0 - stress / (YOUNG * strain), 2e-11);
        for component in 0..3 {
            assert_eq!(step.solid.solution().element_strain[0][component], 0.0);
            close(state.strain()[component], strain, 2e-12);
            close(
                step.solid.solution().element_stress[0][component],
                if component == 0 {
                    stress
                } else {
                    YOUNG * strain
                },
                2e-11,
            );
        }
        close(step.solid.stored_energy_j(), expected_storage, 2e-11);
        close(step.solid.fracture_dissipation_j(), expected_loss, 2e-11);
        close(step.energy.fracture_dissipation_j, expected_loss, 2e-11);
        close(step.energy.eigenstrain_work_j, expected_work, 2e-11);
        close(step.energy.constitutive_strain_work_j, expected_work, 2e-11);
        assert_eq!(step.energy.mechanical_work_j, 0.0);
        assert_eq!(step.energy.support_endpoint_work_j, 0.0);
        assert_eq!(step.energy.dead_load_endpoint_work_j, 0.0);
        assert_eq!(step.energy.endpoint_work_defect_j, 0.0);
        assert!(step.energy.parameter_work_j.abs() < 1e-12);
        assert!(step.energy.parameter_fracture_dissipation_j.abs() < 1e-12);
        let reactions = &step.solid.solution().increments[0].reactions_n;
        close(reactions[3], stress * volume / 0.03, 2e-11);
        close(reactions[7], YOUNG * strain * volume / 0.02, 2e-11);
        close(reactions[11], YOUNG * strain * volume / 0.01, 2e-11);
        check_balances(&step.energy);
    });
}

#[test]
fn reheating_and_mechanical_unloading_release_storage_without_healing() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = expansion_card("unloading");
        let reference = [resolved_expansion(&source, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&source, 300.0).unwrap()];
        let warm = [reference[0].at_temperature(&source, 350.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let cooled = solid
            .advance(
                &cold,
                &fixed(),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let storage = cooled.stored_energy_j();
        let loss = cooled.fracture_dissipation_j();
        let damage = cooled.solution().element_states[0].damage();
        let ratio: f64 = (0.5 * (1.5e-5 + 2e-5) * 50.0) / (0.5 * (1e-5 + 2e-5) * 100.0);
        let heated = cooled
            .advance(
                &warm,
                &fixed(),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        close(
            heated.solid.stored_energy_j(),
            storage * ratio * ratio,
            2e-11,
        );
        close(
            heated.energy.eigenstrain_work_j,
            storage * (ratio * ratio - 1.0),
            2e-11,
        );
        assert!(heated.energy.eigenstrain_work_j < 0.0);
        assert_eq!(heated.energy.mechanical_work_j, 0.0);
        assert!(heated.energy.parameter_work_j.abs() < 1e-12);
        assert!(heated.energy.fracture_dissipation_j.abs() < 1e-10);
        close(heated.solid.fracture_dissipation_j(), loss, 2e-11);
        close(
            heated.solid.solution().element_states[0].damage(),
            damage,
            2e-11,
        );
        for component in 0..3 {
            close(
                heated.solid.solution().element_stress[0][component],
                ratio * cooled.solution().element_stress[0][component],
                2e-11,
            );
        }
        check_balances(&heated.energy);
        let restored = heated
            .solid
            .advance(
                &reference,
                &fixed(),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert!(restored.solid.stored_energy_j().abs() < 1e-12);
        assert!(
            restored.solid.solution().element_stress[0]
                .iter()
                .all(|x| x.abs() < 1e-5)
        );
        close(restored.solid.fracture_dissipation_j(), loss, 2e-11);
        close(
            restored.solid.solution().element_states[0].damage(),
            damage,
            2e-11,
        );
        check_balances(&restored.energy);
        let reloaded = restored
            .solid
            .advance(
                &cold,
                &fixed(),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        close(reloaded.solid.stored_energy_j(), storage, 2e-11);
        close(reloaded.solid.fracture_dissipation_j(), loss, 2e-11);
        assert!(reloaded.energy.fracture_dissipation_j.abs() < 1e-10);

        // A distinct branch releases the SAME stored energy by allowing the
        // cold supports to contract, so work belongs to mechanical motion.
        let relaxed = cooled
            .advance(
                &cold,
                &affine_targets(&nodes, [cold[0].free_linear_strain(); 3]),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        close(relaxed.energy.mechanical_work_j, -storage, 2e-11);
        assert_eq!(relaxed.energy.eigenstrain_work_j, 0.0);
        assert!(relaxed.solid.stored_energy_j().abs() < 1e-12);
        close(relaxed.solid.fracture_dissipation_j(), loss, 2e-11);
        close(
            relaxed.solid.solution().element_states[0].damage(),
            damage,
            2e-11,
        );
        assert!(relaxed.energy.fracture_dissipation_j.abs() < 1e-10);
        // Endpoint reactions vanish at full unloading; their quadrature defect
        // is identified explicitly instead of being labeled fracture or heat.
        close(relaxed.energy.endpoint_work_defect_j, storage, 2e-11);
        check_balances(&relaxed.energy);
    });
}

#[test]
fn thermal_source_reference_and_late_work_refusals_preserve_the_complete_state() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = expansion_card("rollback");
        let reference = [resolved_expansion(&source, &geometry, 400.0, 400.0)];
        let cold = [reference[0].at_temperature(&source, 300.0).unwrap()];
        assert!(matches!(
            ThermalBrittleTetSolid::try_new(&geometry, &cold, TetAssemblyBudget::standard(), cx),
            Err(ThermalBrittleTetError::Solid(BrittleTetError::Invalid(_)))
        ));
        let virgin = ThermalBrittleTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let accepted = virgin
            .advance(
                &cold,
                &fixed(),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let snapshot = (*accepted.solution()).clone();
        let storage = accepted.stored_energy_j();
        let loss = accepted.fracture_dissipation_j();
        let shifted_reference = [resolved_expansion(&source, &geometry, 350.0, 300.0)];
        let foreign_card = expansion_card("different-original-source");
        let foreign = [resolved_expansion(&foreign_card, &geometry, 400.0, 300.0)];
        for replacement in [&shifted_reference, &foreign] {
            assert!(matches!(
                accepted.advance(
                    replacement,
                    &fixed(),
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    cx,
                ),
                Err(ThermalBrittleTetError::Source { .. })
            ));
        }
        let moved = affine_targets(&nodes, [0.0001, 0.0, 0.0]);
        let strict = BrittleTetControls {
            maximum_endpoint_work_defect_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            accepted.advance(&cold, &moved, &[], IncrementSettings::default(), strict, cx),
            Err(ThermalBrittleTetError::Solid(BrittleTetError::EnergyGate {
                account: "thermal mechanical endpoint integration defect", defect_j, ..
            })) if defect_j.abs() > 0.0
        ));
        let stalled = IncrementSettings {
            max_newton_iterations: 1,
            ..IncrementSettings::default()
        };
        assert!(matches!(
            accepted.advance(&cold, &moved, &[], stalled, controls(), cx),
            Err(ThermalBrittleTetError::Solid(BrittleTetError::Solve(
                SmallStrainTetError::NewtonStalled { .. }
            )))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                accepted.advance(
                    &cold,
                    &moved,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    cancelled
                ),
                Err(ThermalBrittleTetError::Solid(BrittleTetError::Cancelled))
            ));
        });
        assert_eq!(accepted.solution(), &snapshot);
        assert_eq!(accepted.stored_energy_j(), storage);
        assert_eq!(accepted.fracture_dissipation_j(), loss);
        assert_eq!(
            accepted.materials()[0]
                .band()
                .material()
                .state_coordinate("T"),
            Some(300.0)
        );
        let resumed = accepted
            .advance(
                &cold,
                &moved,
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let repeated = accepted
            .advance(
                &cold,
                &moved,
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(resumed.solid.solution(), repeated.solid.solution());
        assert_eq!(resumed.energy, repeated.energy);
        assert!(resumed.solid.fracture_dissipation_j() > loss);
        check_balances(&resumed.energy);
    });
}
