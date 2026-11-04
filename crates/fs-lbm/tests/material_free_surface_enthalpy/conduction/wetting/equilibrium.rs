//! Thermal adhesion uses the existing accepted heat and native flow history.

use super::*;
use fs_lbm::material::free_surface::wetting::EquilibriumAngleOrigin;
use fs_material::liquid::surface::EQUILIBRIUM_CONTACT_ANGLE_PROPERTY;
use fs_material::liquid::surface::adhesion::EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY;

fn equilibrium_interface(
    liquid: &MaterialCard,
    solid: &MaterialCard,
    origin: EquilibriumAngleOrigin,
    slope: f64,
    upper_temperature: f64,
) -> InterfaceSystemCard {
    let original = interface(liquid, solid, 0.0, 0.0, 500.0, 100.0);
    let (name, dims, reference) = match origin {
        EquilibriumAngleOrigin::SourcedAngle => (
            EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
            Dims::NONE,
            core::f64::consts::FRAC_PI_2,
        ),
        EquilibriumAngleOrigin::YoungDupre => (
            EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
            SURFACE_TENSION_DIMS,
            0.004,
        ),
    };
    let value = if slope == 0.0 {
        PropertyValue::Scalar {
            value: reference,
            dims,
        }
    } else {
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (300.0, (reference - 100.0 * slope).max(0.25 * reference)),
                (400.0, reference),
                (500.0, reference + 100.0 * slope),
            ],
            dims,
        }
    };
    let claim = PropertyClaim {
        key: PropertyKey::new(name, dims),
        value,
        validity: ValidityDomain::unconstrained().with("T", 300.0, upper_temperature),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        observations: vec![],
        provenance: Provenance {
            source: "synthetic reversible wetting thermodynamics; not measured".into(),
            license: "test-only".into(),
            artifact: None,
        },
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim.clone()).unwrap();
    if upper_temperature < 500.0 {
        // A different supported source exists beyond the first source's range.
        // A thermal trajectory must refuse, rather than silently switching.
        let mut fallback = claim;
        fallback.validity =
            ValidityDomain::unconstrained().with("T", upper_temperature + 0.000001, 500.0);
        claims.insert_claim(fallback).unwrap();
    }
    InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

fn with_equilibrium<R>(
    card: &MaterialCard,
    interface: &InterfaceSystemCard,
    origin: EquilibriumAngleOrigin,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(card, false);
    let law = LiquidTemperatureLaw::try_from_material(
        card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &cells(), cx)
            .unwrap();
        let thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        let mut thermal = thermal
            .with_equilibrium_planar_wetting(
                PlanarWall3::XMin,
                interface,
                origin,
                MaterialPropertySelection::SingleClaimOnly,
                64,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

fn check_equilibrium(thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>, slope: f64) {
    let original = thermal.wetting_source().unwrap();
    let selected = original
        .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
        .unwrap()
        .answer()
        .receipt
        .selected;
    for (&position, sample) in thermal.wetting_samples().unwrap() {
        let t = sample
            .material
            .query_point()
            .iter()
            .find(|(axis, _)| axis == "T")
            .unwrap()
            .1;
        if let Some(material) = thermal.material_at(position) {
            near(t, temperature(material), 2e-12);
        } else {
            assert_eq!(
                t, 400.0,
                "dry wall ghosts retain the declared reference context"
            );
        }
        let expected_work = 0.004 + slope * (t - 400.0);
        let work = sample
            .material
            .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
            .unwrap();
        near(work.value_si(), expected_work, 2e-12);
        assert_eq!(work.answer().receipt.selected, selected);
        near(
            sample.applied_angle_radians,
            (expected_work / 0.004 - 1.0).acos(),
            2e-12,
        );
        assert_eq!(sample.equilibrium, Some(EquilibriumAngleOrigin::YoungDupre));
        assert_eq!(sample.regime, ContactLineRegime2::Equilibrium);
        assert_eq!(sample.signed_speed_m_s, 0.0);
        assert!(sample.hydrodynamic.is_none());
        assert!(
            sample
                .material
                .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
                .is_none()
        );
        assert!(
            sample
                .material
                .property(ADVANCING_CONTACT_ANGLE_PROPERTY)
                .is_none()
        );
        assert!(
            !sample
                .material
                .query_point()
                .iter()
                .any(|(axis, _)| axis == DYNAMIC_WETTING_SPEED_AXIS)
        );
    }
}

#[test]
fn g1_constant_adhesion_and_measured_angle_keep_the_same_native_thermal_trajectory() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let run = |origin| {
        let interface = equilibrium_interface(&liquid, &solid, origin, 0.0, 500.0);
        with_equilibrium(&liquid, &interface, origin, |thermal, gate| {
            let original = thermal.wetting_source().unwrap().clone();
            let initial = thermal.wetting_samples().unwrap().clone();
            for sample in initial.values() {
                assert_eq!(sample.material, original);
                assert_eq!(sample.applied_angle_radians, core::f64::consts::FRAC_PI_2);
            }
            let energy = thermal.total_enthalpy_j().unwrap();
            let mut result = Vec::new();
            for _ in 0..3 {
                let report = thermal.step(gate).unwrap();
                near(thermal.total_enthalpy_j().unwrap(), energy, 2e-11);
                assert_eq!(report.external_heat_j, 0.0);
                for sample in thermal.wetting_samples().unwrap().values() {
                    assert_eq!(sample.material, original);
                }
                result.push(thermal.flow().lattice().grid().state_bits());
            }
            result
        })
    };
    assert_eq!(
        run(EquilibriumAngleOrigin::YoungDupre),
        run(EquilibriumAngleOrigin::SourcedAngle)
    );
}

#[test]
fn g1_g3_local_thermal_adhesion_drives_native_momentum_and_retains_original_sources() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let run = |slope| {
        let interface = equilibrium_interface(
            &liquid,
            &solid,
            EquilibriumAngleOrigin::YoungDupre,
            slope,
            500.0,
        );
        with_equilibrium(
            &liquid,
            &interface,
            EquilibriumAngleOrigin::YoungDupre,
            |thermal, gate| {
                let initial_energy = thermal.total_enthalpy_j().unwrap();
                let heat = [FreeSurfaceCellHeat {
                    coordinates: WALL_CELL,
                    energy_j: 0.04,
                }];
                let report = thermal.step_with_heat(&heat, gate).unwrap();
                near(
                    thermal.total_enthalpy_j().unwrap() - initial_energy,
                    0.04,
                    2e-11,
                );
                near(report.external_heat_j, 0.04, 2e-13);
                near(report.mass_before_kg, report.mass_after_kg, 2e-13);
                assert_eq!(report.viscous_heat_j, 0.0);
                check_equilibrium(thermal, slope);
                let boundary = &thermal.wetting_samples().unwrap()[&WALL_CELL];
                assert!(
                    boundary
                        .material
                        .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                        .unwrap()
                        .value_si()
                        >= 0.004
                );
                thermal.flow().lattice().grid().state_bits()
            },
        )
    };
    let fixed = run(0.0);
    let responsive = run(0.00002);
    assert_eq!(fixed.len(), responsive.len());
    let change = fixed
        .iter()
        .zip(&responsive)
        .map(|(&a, &b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        change > 1e-12,
        "adhesion at the heated wall must reach the same native momentum step"
    );
    assert_eq!(
        responsive,
        run(0.00002),
        "accepted thermal adhesion is deterministic"
    );
}

#[test]
fn g3_g4_work_support_complete_wetting_and_tension_support_refuse_atomically() {
    let solid = body::wall_card(0.5, 500.0);
    for (work_limit, slope, sigma_limit) in [
        (410.0, 0.00002, 500.0),
        (500.0, 0.0003, 500.0),
        (500.0, 0.00002, 410.0),
    ] {
        let liquid = source_card(0.004, false, sigma_limit);
        let interface = equilibrium_interface(
            &liquid,
            &solid,
            EquilibriumAngleOrigin::YoungDupre,
            slope,
            work_limit,
        );
        with_equilibrium(
            &liquid,
            &interface,
            EquilibriumAngleOrigin::YoungDupre,
            |thermal, gate| {
                let populations = thermal.flow().lattice().grid().state_bits();
                let samples = thermal.wetting_samples().unwrap().clone();
                let before_h = cells()
                    .iter()
                    .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                    .collect::<Vec<_>>();
                let energy = thermal.total_enthalpy_j().unwrap();
                let heat = [FreeSurfaceCellHeat {
                    coordinates: WALL_CELL,
                    energy_j: 0.04,
                }];
                assert!(thermal.step_with_heat(&heat, gate).is_err());
                assert_eq!(thermal.flow().lattice().grid().state_bits(), populations);
                assert_eq!(thermal.wetting_samples().unwrap(), &samples);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
                assert_eq!(
                    cells()
                        .iter()
                        .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                        .collect::<Vec<_>>(),
                    before_h
                );
                thermal.step(gate).unwrap();
                check_equilibrium(thermal, slope);
            },
        );
    }
}

#[test]
fn g3_equilibrium_adhesion_cannot_silently_acquire_cox_dynamics() {
    use fs_lbm::material::free_surface::wetting::{CoxVoinovControls, CoxVoinovLaw};
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let interface = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.00002,
        500.0,
    );
    with_equilibrium(
        &liquid,
        &interface,
        EquilibriumAngleOrigin::YoungDupre,
        |thermal, gate| {
            let law = CoxVoinovLaw::try_new(
                CoxVoinovControls {
                    microscopic_cutoff_m: 1e-8,
                    observation_length_m: 0.001,
                    maximum_capillary_number: 0.01,
                    minimum_angle_radians: 0.1,
                    maximum_angle_radians: 3.0,
                },
                gate,
            )
            .unwrap();
            assert!(thermal.clone().with_cox_voinov_wetting(law, gate).is_err());
            thermal.step(gate).unwrap();
            check_equilibrium(thermal, 0.00002);
        },
    );
}
