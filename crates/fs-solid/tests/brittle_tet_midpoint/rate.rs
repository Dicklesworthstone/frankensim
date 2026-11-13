//! Native fracture viscosity: independent scalar mechanics and energy oracles.

use super::*;
use fs_material::brittle::source::rate::{
    BRITTLE_NORMAL_VISCOSITY_PROPERTY, BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    BRITTLE_VISCOSITY_DIMS, resolve_viscous_brittle_dynamic_material,
};
use fs_qty::{QuantitySpec, Temperature};
use fs_solid::brittle::solid::BrittleTetControls;
use fs_solid::brittle::solid::dynamics::BrittleTetDynamicControls;
use fs_solid::brittle::solid::dynamics::time::BrittleTetTimeControls;

const ETA_NORMAL: f64 = 60.0;
const ETA_TANGENTIAL: f64 = 90.0;
const WIDTH: f64 = 1.0 / 3.0;

fn viscous_card() -> MaterialCard {
    let base = material_card();
    let mut claims = base.claims().clone();
    for (name, value) in [
        (BRITTLE_NORMAL_VISCOSITY_PROPERTY, ETA_NORMAL),
        (BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY, ETA_TANGENTIAL),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, BRITTLE_VISCOSITY_DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: BRITTLE_VISCOSITY_DIMS,
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured fracture-plane viscosity; not measured data".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn viscous_source(geometry: &TetLoadGeometry3<'_>, maximum_rate: f64) -> ResolvedBrittleCrackBand {
    let card = viscous_card();
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1., 0., 0.]).unwrap();
    near(band.band_width_m(), WIDTH, 1e-15);
    resolve_viscous_brittle_dynamic_material(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.02, 1e-9),
        maximum_rate,
    )
    .unwrap()
}

fn endpoint_controls() -> BrittleTetDynamicControls {
    BrittleTetDynamicControls {
        work: BrittleTetControls {
            material_energy_tolerance_j: 1e-10,
            equilibrium_work_tolerance_j: 1e-10,
            maximum_endpoint_work_defect_j: 1.0,
        },
        kinetic_work_tolerance_j: 1e-10,
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

#[test]
fn native_normal_and_shear_dynamics_use_physical_viscosity_and_actual_endpoint_rate() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [viscous_source(&geometry, 1.0)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let dt = 0.01;
        let initial_speed = 0.003;
        // Both cases have one free coordinate with consistent mass rho*V/10.
        // Shear uses engineering gamma=u_y/L, so its modulus is mu, not 2mu.
        for (component, stress_slot, modulus, viscosity) in [
            (0, 0, YOUNG, ETA_NORMAL),
            (1, 3, YOUNG / 2.0, ETA_TANGENTIAL),
        ] {
            let dof = 3 + component;
            let supports: Vec<_> = (0..12)
                .filter(|&index| index != dof)
                .map(|index| (index, 0.0))
                .collect();
            let mut velocity = [[0.0; 3]; 4];
            velocity[1][component] = initial_speed;
            let motion = BrittleTetMotion::try_from_solid(&solid, &velocity, 0.0, 4, cx).unwrap();
            let eta_bulk = WIDTH * viscosity;
            let mass = DENSITY * VOLUME / 10.0;
            let stiffness = VOLUME * modulus;
            let damping = VOLUME * eta_bulk;
            near(
                motion.solid().solution().element_stress[0][stress_slot],
                eta_bulk * initial_speed,
                1e-14,
            );
            assert_eq!(motion.solid().solution().element_states[0].damage(), 0.0);

            let backward = motion
                .advance(dt, &supports, &[], settings(), endpoint_controls(), cx)
                .unwrap();
            let u_be = (mass * initial_speed / dt) / (mass / (dt * dt) + damping / dt + stiffness);
            let v_be = u_be / dt;
            near(
                backward.motion.solid().solution().displacement_m[1][component],
                u_be,
                2e-13,
            );
            near(backward.motion.velocity_mps()[1][component], v_be, 2e-12);
            near(
                backward.motion.solid().solution().element_stress[0][stress_slot],
                modulus * u_be + eta_bulk * v_be,
                2e-11,
            );
            near(
                backward.energy.viscous_dissipation_j,
                damping * v_be * v_be * dt,
                2e-14,
            );
            near(
                backward.energy.inertia.numerical_loss_j,
                0.5 * mass * (v_be - initial_speed).powi(2),
                2e-14,
            );
            near(backward.energy.total_balance_error_j, 0.0, 1e-12);

            let midpoint = motion
                .advance_midpoint(dt, &supports, &[], settings(), controls(), cx)
                .unwrap();
            let u_mid = (2.0 * mass * initial_speed / dt)
                / (2.0 * mass / (dt * dt) + damping / dt + stiffness / 2.0);
            let mean_speed = u_mid / dt;
            let v_mid = 2.0 * mean_speed - initial_speed;
            near(
                midpoint.motion.solid().solution().displacement_m[1][component],
                u_mid,
                2e-13,
            );
            near(midpoint.motion.velocity_mps()[1][component], v_mid, 2e-12);
            near(
                midpoint.mean_element_stress_pa[0][stress_slot],
                modulus * u_mid / 2.0 + eta_bulk * mean_speed,
                2e-11,
            );
            near(
                midpoint.motion.solid().solution().element_stress[0][stress_slot],
                modulus * u_mid + eta_bulk * v_mid,
                2e-11,
            );
            assert!((v_mid - mean_speed).abs() > 1e-4);
            near(
                midpoint.energy.viscous_dissipation_j,
                damping * mean_speed * mean_speed * dt,
                2e-14,
            );
            assert_eq!(midpoint.energy.fracture_dissipation_j, 0.0);
            assert_eq!(midpoint.energy.inertia.numerical_loss_j, 0.0);
            near(
                midpoint.motion.kinetic_energy_j()
                    + midpoint.motion.solid().stored_energy_j()
                    + midpoint.energy.viscous_dissipation_j,
                motion.kinetic_energy_j(),
                2e-12,
            );
            assert_eq!(motion.solid().solution().displacement_m, [[0.0; 3]; 4]);
        }
    });
}

#[test]
fn rupture_integrates_surviving_bond_loss_and_failed_reclosure_has_no_dashpot() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [viscous_source(&geometry, 0.005)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let motion = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let strain = 0.007;
        let dt = 1.0;
        let failed = motion
            .advance_midpoint(dt, &prescribed(strain), &[], settings(), controls(), cx)
            .unwrap();
        // Independent integral of bilinear integrity: before onset it is one;
        // afterwards delta0*(delta_f-delta)/((delta_f-delta0)*delta).
        let delta0 = WIDTH * STRENGTH / YOUNG;
        let delta_f = 2.0 * GC / STRENGTH;
        let intact_path = delta0 * delta_f / (delta_f - delta0) * (delta_f / delta0).ln();
        let expected_viscous = (VOLUME / WIDTH) * ETA_NORMAL * (WIDTH * strain / dt) * intact_path;
        near(failed.energy.viscous_dissipation_j, expected_viscous, 2e-12);
        near(
            failed.energy.fracture_dissipation_j,
            (VOLUME / WIDTH) * GC,
            2e-12,
        );
        assert_eq!(
            failed.motion.solid().solution().element_states[0].damage(),
            1.0
        );
        near(
            failed.motion.solid().solution().element_stress[0][0],
            0.0,
            1e-12,
        );
        near(
            failed.energy.material_work_j,
            expected_viscous + (VOLUME / WIDTH) * GC,
            2e-12,
        );
        near(failed.energy.total_balance_error_j, 0.0, 2e-12);

        // Previously failed bonds no longer impose their viscosity rate bound.
        // Original normal compression remains, with no healed fracture or eta.
        let closed = failed
            .motion
            .advance_midpoint(0.01, &prescribed(-0.0002), &[], settings(), controls(), cx)
            .unwrap();
        assert_eq!(closed.energy.viscous_dissipation_j, 0.0);
        assert_eq!(closed.energy.fracture_dissipation_j, 0.0);
        assert_eq!(
            closed.motion.solid().solution().element_states[0].damage(),
            1.0
        );
        near(
            closed.motion.solid().solution().element_stress[0][0],
            -YOUNG * 0.0002,
            1e-10,
        );
        near(
            closed.motion.solid().fracture_dissipation_j(),
            (VOLUME / WIDTH) * GC,
            2e-12,
        );
    });
}

#[test]
fn adaptive_paths_sum_physical_viscosity_separately_from_fracture_and_numerics() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [viscous_source(&geometry, 1.0)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let motion = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 64, cx).unwrap();
        let duration = 0.1;
        let strain = 0.0002;
        let time_controls = BrittleTetTimeControls {
            initial_interval_s: duration,
            minimum_interval_s: duration / 1024.0,
            maximum_interval_s: duration,
            maximum_intervals: 16,
            maximum_trials: 128,
            displacement_tolerance_m: 0.01,
            velocity_tolerance_mps: 0.1,
            damage_tolerance: 1e-4,
            fracture_dissipation_tolerance_j: 1e-10,
            relative_tolerance: 0.0,
            maximum_damage_increment: 0.1,
            work: endpoint_controls(),
        };
        let backward = motion
            .advance_time(
                duration,
                &prescribed(strain),
                &[],
                settings(),
                time_controls,
                cx,
            )
            .unwrap();
        let midpoint = motion
            .advance_midpoint_time(
                duration,
                &prescribed(strain),
                &[],
                settings(),
                BrittleTetTimeControls {
                    initial_interval_s: time_controls.initial_interval_s,
                    minimum_interval_s: time_controls.minimum_interval_s,
                    maximum_interval_s: time_controls.maximum_interval_s,
                    maximum_intervals: time_controls.maximum_intervals,
                    maximum_trials: time_controls.maximum_trials,
                    displacement_tolerance_m: time_controls.displacement_tolerance_m,
                    velocity_tolerance_mps: time_controls.velocity_tolerance_mps,
                    damage_tolerance: time_controls.damage_tolerance,
                    fracture_dissipation_tolerance_j: time_controls
                        .fracture_dissipation_tolerance_j,
                    relative_tolerance: time_controls.relative_tolerance,
                    maximum_damage_increment: time_controls.maximum_damage_increment,
                    work: controls(),
                },
                cx,
            )
            .unwrap();
        let expected_viscous = VOLUME * WIDTH * ETA_NORMAL * strain * strain / duration;
        for (loss, material, fracture, defect) in [
            (
                backward.energy.viscous_dissipation_j,
                backward.energy.material_work_j,
                backward.energy.fracture_dissipation_j,
                backward.energy.total_balance_error_j,
            ),
            (
                midpoint.energy.viscous_dissipation_j,
                midpoint.energy.material_work_j,
                midpoint.energy.fracture_dissipation_j,
                midpoint.energy.total_balance_error_j,
            ),
        ] {
            near(loss, expected_viscous, 1e-13);
            near(
                material,
                0.5 * VOLUME * YOUNG * strain * strain + expected_viscous,
                1e-13,
            );
            assert_eq!(fracture, 0.0);
            near(defect, 0.0, 1e-12);
        }
        near(
            backward
                .records
                .iter()
                .flat_map(|record| record.energy)
                .map(|energy| energy.viscous_dissipation_j)
                .sum(),
            expected_viscous,
            1e-13,
        );
        near(
            midpoint
                .records
                .iter()
                .flat_map(|record| record.energy)
                .map(|energy| energy.viscous_dissipation_j)
                .sum(),
            expected_viscous,
            1e-13,
        );
        assert_eq!(midpoint.energy.inertia.numerical_loss_j, 0.0);
    });
}

#[test]
fn endpoint_rate_and_cancellation_refusals_preserve_the_entire_accepted_motion() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [viscous_source(&geometry, 0.001)];
        let solid = BrittleTetSolid::try_new(
            &geometry,
            &sources,
            TetAssemblyBudget {
                maximum_free_dofs: 0,
                ..budget()
            },
            cx,
        )
        .unwrap();
        let mut fast = [[0.0; 3]; 4];
        fast[1][0] = 0.004;
        assert!(BrittleTetMotion::try_from_solid(&solid, &fast, 0.0, 4, cx).is_err());
        let motion = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let saved = motion.solid().solution().clone();
        // The secant jump rate is admissible; the actual centered endpoint
        // velocity is twice as large and must refuse before publication.
        assert!(
            motion
                .advance_midpoint(0.1, &prescribed(0.0002), &[], settings(), controls(), cx)
                .is_err()
        );
        let accepted = motion
            .advance(
                0.1,
                &prescribed(0.0002),
                &[],
                settings(),
                endpoint_controls(),
                cx,
            )
            .unwrap();
        assert!(accepted.energy.viscous_dissipation_j > 0.0);
        assert_eq!(motion.solid().solution(), &saved);
        assert_eq!(motion.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(motion.time_s(), 0.0);
        assert_eq!(motion.accepted_steps(), 0);
        gate.request();
        assert!(matches!(
            motion.advance(
                0.1,
                &prescribed(0.0001),
                &[],
                settings(),
                endpoint_controls(),
                cx
            ),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::Cancelled))
        ));
        assert_eq!(motion.solid().solution(), &saved);
    });
}

#[test]
fn temperature_continuation_keeps_rate_source_work_and_original_reference_mass() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let base = viscous_card();
        let mut claims = ClaimSet::new();
        for (_, original) in base.claims().claims_ordered() {
            let mut claim = original.clone();
            let factor = match claim.key.name() {
                YOUNG_MODULUS_PROPERTY | DENSITY_PROPERTY => Some(0.5),
                BRITTLE_NORMAL_VISCOSITY_PROPERTY | BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY => {
                    Some(2.0)
                }
                _ => None,
            };
            if let Some(factor) = factor {
                let PropertyValue::Scalar { value, dims } = claim.value else {
                    panic!("scalar fixture required");
                };
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, value), (600.0, factor * value)],
                    dims,
                };
                // Typed like every query of this card: an untyped "T" axis
                // cannot answer a query that types it in kelvin.
                claim.validity = ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    600.0,
                );
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
            claims.insert_claim(claim).unwrap();
        }
        let card =
            MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap();
        let band = TetCrackBandGeometry::from_element(&geometry, 0, [1., 0., 0.]).unwrap();
        let query = QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap();
        let cold = [resolve_viscous_brittle_dynamic_material(
            &card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            band.material_settings(0.02, 1e-9),
            1.0,
        )
        .unwrap()];
        let hot = [cold[0].at_temperature(&card, "T", 600.0).unwrap()];
        let solid = BrittleTetSolid::try_new(&geometry, &cold, budget(), cx).unwrap();
        let virgin = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let initial_strain = 0.0002;
        let old = virgin
            .advance(
                1.0,
                &prescribed(initial_strain),
                &[],
                settings(),
                endpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let saved = old.solid().solution().clone();
        let strain = 0.0003;
        let dt = 0.1;
        let mean_rate = (strain - initial_strain) / dt;
        let hot_young = 0.5 * YOUNG;
        let hot_eta = 2.0 * ETA_NORMAL * WIDTH;
        let parameter_work = 0.5 * VOLUME * (hot_young - YOUNG) * initial_strain * initial_strain;
        let viscous_loss = VOLUME * hot_eta * mean_rate * mean_rate * dt;
        let backward = old
            .advance_at_temperature(
                dt,
                &hot,
                "T",
                &prescribed(strain),
                &[],
                settings(),
                endpoint_controls(),
                cx,
            )
            .unwrap();
        let midpoint = old
            .advance_midpoint_at_temperature(
                dt,
                &hot,
                "T",
                &prescribed(strain),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        for (next, transfer_work, loss, defect) in [
            (
                &backward.mechanical.motion,
                backward.transfer.parameter_work_j,
                backward.mechanical.energy.viscous_dissipation_j,
                backward.total_balance_error_j,
            ),
            (
                &midpoint.mechanical.motion,
                midpoint.transfer.parameter_work_j,
                midpoint.mechanical.energy.viscous_dissipation_j,
                midpoint.total_balance_error_j,
            ),
        ] {
            near(transfer_work, parameter_work, 2e-13);
            near(loss, viscous_loss, 2e-13);
            near(defect, 0.0, 2e-12);
            assert_eq!(next.reference_mass_kg(), old.reference_mass_kg());
            assert_eq!(next.reference_density_kg_m3(), &[DENSITY]);
            near(
                next.solid().materials()[0]
                    .reference_density_kg_m3()
                    .unwrap(),
                0.5 * DENSITY,
                1e-14,
            );
            assert!(std::ptr::eq(
                next.solid().materials()[0].law(),
                hot[0].law()
            ));
            assert_eq!(next.solid().fracture_dissipation_j(), 0.0);
        }
        near(
            backward.mechanical.motion.solid().solution().element_stress[0][0],
            hot_young * strain + hot_eta * mean_rate,
            2e-11,
        );
        let endpoint_rate = 2.0 * mean_rate - old.velocity_mps()[1][0];
        near(
            midpoint.mechanical.motion.solid().solution().element_stress[0][0],
            hot_young * strain + hot_eta * endpoint_rate,
            2e-11,
        );
        assert_eq!(old.solid().solution(), &saved);
    });
}
