//! Timed thermal eigenstrain, physical bond viscosity and native endpoint rates.
//! Manufactured coefficients provide independent mechanics and work oracles.

use super::*;
use fs_material::brittle::source::expansion::rate::resolve_thermal_viscous_brittle_dynamic_material;
use fs_material::brittle::source::rate::{
    BRITTLE_NORMAL_VISCOSITY_PROPERTY, BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    BRITTLE_VISCOSITY_DIMS,
};
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::expansion::initiation::ThermalTensileElasticTetSolid;

const ETA_NORMAL: f64 = 60.0;
const ETA_TANGENTIAL: f64 = 90.0;
const WIDTH: f64 = 1.0 / 3.0;
const BULK_ETA: f64 = WIDTH * ETA_NORMAL;

fn viscous_card(tag: &str, alpha: f64, varying_density: bool) -> MaterialCard {
    let base = card(tag, alpha, varying_density);
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

fn viscous_source(
    card: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
    reference_k: f64,
    maximum_rate_m_s: f64,
) -> ResolvedThermalBrittle {
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    near(band.band_width_m(), WIDTH, 1e-15);
    resolve_thermal_viscous_brittle_dynamic_material(
        card,
        &temperature(reference_k),
        &temperature(reference_k),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.02, 1e-9),
        maximum_rate_m_s,
    )
    .unwrap()
}

fn midpoint_controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        kinetic_work_tolerance_j: 1e-10,
        total_energy_tolerance_j: 1e-10,
        difference_step_strain: Some(1e-7),
    }
}

#[test]
fn free_expansion_cancels_mechanical_rate_but_centered_endpoint_is_not_the_chord_rate() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = viscous_card("matched expansion rate", 2e-6, true);
        let reference = [viscous_source(&card, &geometry, 300.0, 1.0)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let dt = 0.1;
        let free = hot[0].free_linear_strain();
        let rate = free / dt;
        near(free, 0.0002, 1e-16);
        near(
            reference[0]
                .free_strain_rate_from_temperature_rate(100.0 / dt)
                .unwrap(),
            rate,
            1e-15,
        );
        let velocity = NODES.map(|point| point.map(|coordinate| rate * coordinate));
        let old = ThermalBrittleTetMotion::try_from_solid_with_free_strain_rates(
            &solid,
            &velocity,
            &[rate],
            0.0,
            4,
            cx,
        )
        .unwrap();
        for &stress in &old.solid().solution().element_stress[0] {
            near(stress, 0.0, 1e-13);
        }
        let backward = old
            .advance(
                &hot,
                dt,
                &affine_targets(free),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let midpoint = old
            .advance_midpoint(
                &hot,
                dt,
                &affine_targets(free),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        for (next, loss, strain_work, thermal_work, defect) in [
            (
                &backward.motion,
                backward.energy.viscous_dissipation_j,
                backward.energy.constitutive_strain_work_j,
                backward.energy.eigenstrain_work_j,
                backward.energy.total_balance_error_j,
            ),
            (
                &midpoint.motion,
                midpoint.energy.viscous_dissipation_j,
                midpoint.energy.constitutive_strain_work_j,
                midpoint.energy.eigenstrain_work_j,
                midpoint.energy.total_balance_error_j,
            ),
        ] {
            near(loss, 0.0, 1e-16);
            near(strain_work, 0.0, 1e-16);
            near(thermal_work, 0.0, 1e-16);
            near(defect, 0.0, 1e-12);
            near(next.free_linear_strain_rates_per_s()[0], rate, 1e-15);
            for &stress in &next.solid().solution().element_stress[0] {
                near(stress, 0.0, 1e-12);
            }
            assert_eq!(next.solid().solution().element_states[0].damage(), 0.0);
            near(next.kinetic_energy_j(), old.kinetic_energy_j(), 1e-14);
            assert_eq!(next.reference_mass_kg(), old.reference_mass_kg());
            assert_eq!(next.reference_density_kg_m3(), &[DENSITY]);
            near(
                next.solid().materials()[0]
                    .band()
                    .reference_density_kg_m3()
                    .unwrap(),
                DENSITY / 2.0,
                1e-14,
            );
        }

        // The same zero mechanical chord from rest gives v1=2*free/dt.
        // Its endpoint mechanical rate is free/dt, although path viscous work
        // is zero. Publishing either mean or secant endpoint stress is wrong.
        let rest =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let accelerated = rest
            .advance_midpoint(
                &hot,
                dt,
                &affine_targets(free),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        near(accelerated.energy.viscous_dissipation_j, 0.0, 1e-16);
        near(accelerated.motion.velocity_mps()[1][0], 2.0 * rate, 1e-14);
        near(
            accelerated.motion.solid().solution().element_stress[0][0],
            BULK_ETA * rate,
            1e-12,
        );
        for axis in [1, 2] {
            near(
                accelerated.motion.solid().solution().element_stress[0][axis],
                0.0,
                1e-12,
            );
        }
    });
}

#[test]
fn restrained_heating_and_cooling_drive_physical_viscosity_without_nodal_motion() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        for (reference_k, next_k) in [(300.0, 400.0), (400.0, 300.0)] {
            let card = viscous_card("restrained rate work", 2e-6, false);
            let reference = [viscous_source(&card, &geometry, reference_k, 1.0)];
            let next = [reference[0].at_temperature(&card, next_k).unwrap()];
            let solid =
                ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
            let old = ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx)
                .unwrap();
            let dt = 0.1;
            let free = next[0].free_linear_strain();
            let loss = VOLUME * BULK_ETA * free * free / dt;
            let storage = 1.5 * VOLUME * YOUNG * free * free;
            let backward = old
                .advance(
                    &next,
                    dt,
                    &affine_targets(0.0),
                    &[],
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            let midpoint = old
                .advance_midpoint(
                    &next,
                    dt,
                    &affine_targets(0.0),
                    &[],
                    settings(),
                    midpoint_controls(),
                    cx,
                )
                .unwrap();
            for (motion, viscosity, mechanical, thermal, constitutive, defect) in [
                (
                    &backward.motion,
                    backward.energy.viscous_dissipation_j,
                    backward.energy.mechanical_work_j,
                    backward.energy.eigenstrain_work_j,
                    backward.energy.constitutive_strain_work_j,
                    backward.energy.total_balance_error_j,
                ),
                (
                    &midpoint.motion,
                    midpoint.energy.viscous_dissipation_j,
                    midpoint.energy.mechanical_work_j,
                    midpoint.energy.eigenstrain_work_j,
                    midpoint.energy.constitutive_strain_work_j,
                    midpoint.energy.total_balance_error_j,
                ),
            ] {
                near(viscosity, loss, 1e-13);
                near(mechanical, 0.0, 1e-16);
                near(thermal, storage + loss, 1e-13);
                near(constitutive, storage + loss, 1e-13);
                near(defect, 0.0, 1e-12);
                near(motion.solid().stored_energy_j(), storage, 1e-13);
                near(
                    motion.solid().solution().element_stress[0][0],
                    -YOUNG * free - BULK_ETA * free / dt,
                    1e-12,
                );
                for axis in [1, 2] {
                    near(
                        motion.solid().solution().element_stress[0][axis],
                        -YOUNG * free,
                        1e-12,
                    );
                }
                assert_eq!(motion.velocity_mps(), &[[0.0; 3]; 4]);
                assert_eq!(motion.solid().fracture_dissipation_j(), 0.0);
                near(motion.free_linear_strain_rates_per_s()[0], free / dt, 1e-15);
            }
        }
    });
}

#[test]
fn coupled_cooling_uses_the_mechanical_rate_in_both_native_scalar_momentum_solutions() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = viscous_card("cooling scalar dynamics", 2e-6, false);
        let reference = [viscous_source(&card, &geometry, 400.0, 1.0)];
        let cold = [reference[0].at_temperature(&card, 300.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let old =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let dt = 0.01;
        let free = cold[0].free_linear_strain();
        let mass = DENSITY * VOLUME / 10.0;
        let stiffness = YOUNG * VOLUME;
        let damping = BULK_ETA * VOLUME;
        let backward = old
            .advance(
                &cold,
                dt,
                &axial_supports(),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let midpoint = old
            .advance_midpoint(
                &cold,
                dt,
                &axial_supports(),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        for (motion, viscosity, mechanical_work, thermal_work, defect, factor) in [
            (
                &backward.motion,
                backward.energy.viscous_dissipation_j,
                backward.energy.mechanical_work_j,
                backward.energy.eigenstrain_work_j,
                backward.energy.total_balance_error_j,
                1.0,
            ),
            (
                &midpoint.motion,
                midpoint.energy.viscous_dissipation_j,
                midpoint.energy.mechanical_work_j,
                midpoint.energy.eigenstrain_work_j,
                midpoint.energy.total_balance_error_j,
                2.0,
            ),
        ] {
            // Native consistent mass, not row-sum lumping. Force uses endpoint
            // elasticity for BE, half endpoint elasticity for the path mean;
            // both use the constant mechanical path rate (u-free)/dt.
            let response_stiffness = stiffness / factor + damping / dt;
            let u = response_stiffness * free / (factor * mass / (dt * dt) + response_stiffness);
            let strain = u - free;
            let velocity = factor * u / dt;
            let endpoint_rate = velocity - free / dt;
            let mean_stress_x = YOUNG * strain / 2.0 + BULK_ETA * strain / dt;
            let expected_loss = VOLUME * BULK_ETA * strain * strain / dt;
            let expected_storage = 0.5 * VOLUME * YOUNG * (strain * strain + 2.0 * free * free);
            assert!(strain > 0.0 && strain < STRENGTH / YOUNG);
            near(motion.solid().solution().displacement_m[1][0], u, 2e-12);
            near(motion.velocity_mps()[1][0], velocity, 2e-10);
            near(
                motion.solid().solution().element_states[0].strain()[0],
                strain,
                2e-12,
            );
            near(
                motion.solid().solution().element_stress[0][0],
                YOUNG * strain + BULK_ETA * endpoint_rate,
                2e-10,
            );
            near(viscosity, expected_loss, 2e-12);
            near(mechanical_work, VOLUME * mean_stress_x * u, 2e-12);
            near(
                thermal_work,
                -VOLUME * free * (mean_stress_x - YOUNG * free),
                2e-12,
            );
            near(motion.solid().stored_energy_j(), expected_storage, 2e-12);
            near(defect, 0.0, 2e-12);
            assert_eq!(motion.solid().fracture_dissipation_j(), 0.0);
        }
        assert_eq!(midpoint.energy.inertia.numerical_loss_j, 0.0);
    });
}

#[test]
fn thermal_rupture_integrates_bond_loss_and_failed_reclosure_has_no_viscous_rate_bound() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = viscous_card("thermal viscous rupture", 7e-5, false);
        let reference = [viscous_source(&card, &geometry, 400.0, 0.005)];
        let cold = [reference[0].at_temperature(&card, 300.0).unwrap()];
        let warm = [reference[0].at_temperature(&card, 401.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let old =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let dt = 1.0;
        let free = cold[0].free_linear_strain();
        let failed = old
            .advance_midpoint(
                &cold,
                dt,
                &affine_targets(0.0),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        let delta0 = WIDTH * STRENGTH / YOUNG;
        let delta_f = 2.0 * GC / STRENGTH;
        let integrity_integral = delta0 * delta_f / (delta_f - delta0) * (delta_f / delta0).ln();
        let loss = (VOLUME / WIDTH) * ETA_NORMAL * (-WIDTH * free / dt) * integrity_integral;
        let fracture = (VOLUME / WIDTH) * GC;
        let parallel_storage = VOLUME * YOUNG * free * free;
        near(failed.energy.viscous_dissipation_j, loss, 2e-12);
        near(failed.energy.fracture_dissipation_j, fracture, 2e-12);
        near(
            failed.energy.eigenstrain_work_j,
            parallel_storage + fracture + loss,
            2e-12,
        );
        near(failed.energy.mechanical_work_j, 0.0, 1e-16);
        assert!(failed.motion.solid().solution().element_states[0].failed());
        near(
            failed.motion.solid().solution().element_stress[0][0],
            0.0,
            1e-12,
        );
        near(
            failed.motion.free_linear_strain_rates_per_s()[0],
            free / dt,
            1e-15,
        );

        // This physical thermal rate is far outside the former surviving-bond
        // ceiling. Failed bonds have no eta; original normal compression and
        // in-plane storage survive without resetting the fracture account.
        let closed = failed
            .motion
            .advance_midpoint(
                &warm,
                0.001,
                &affine_targets(0.0),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(closed.energy.viscous_dissipation_j, 0.0);
        assert_eq!(closed.energy.fracture_dissipation_j, 0.0);
        assert!(closed.motion.free_linear_strain_rates_per_s()[0] > 1.0);
        assert!(closed.motion.solid().solution().element_states[0].failed());
        near(
            closed.motion.solid().fracture_dissipation_j(),
            fracture,
            2e-12,
        );
        for axis in 0..3 {
            near(
                closed.motion.solid().solution().element_stress[0][axis],
                -YOUNG * warm[0].free_linear_strain(),
                2e-11,
            );
        }
        near(closed.energy.total_balance_error_j, 0.0, 2e-12);
    });
}

#[test]
fn source_rate_static_and_cancellation_refusals_preserve_the_complete_thermal_motion() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let card = viscous_card("thermal rate refusal", 4e-6, false);
        let reference = [viscous_source(&card, &geometry, 300.0, 0.001)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let solid = ThermalBrittleTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget {
                maximum_free_dofs: 0,
                ..budget()
            },
            cx,
        )
        .unwrap();
        assert!(
            solid
                .advance(
                    &hot,
                    &affine_targets(0.0),
                    &[],
                    settings(),
                    controls().work,
                    cx
                )
                .is_err()
        );
        assert!(
            ThermalTensileElasticTetSolid::try_new(&geometry, &reference, budget(), cx).is_err()
        );
        for rates in [&[][..], &[f64::NAN][..], &[0.004][..]] {
            assert!(
                ThermalBrittleTetMotion::try_from_solid_with_free_strain_rates(
                    &solid,
                    &[[0.0; 3]; 4],
                    rates,
                    0.0,
                    4,
                    cx
                )
                .is_err()
            );
        }
        // The source ceiling bounds the combined mechanical rate. Large
        // matching total/free rates cancel physically and must remain admissible.
        let matching_velocity = NODES.map(|point| point.map(|coordinate| 0.004 * coordinate));
        let matched = ThermalBrittleTetMotion::try_from_solid_with_free_strain_rates(
            &solid,
            &matching_velocity,
            &[0.004],
            0.0,
            4,
            cx,
        )
        .unwrap();
        for &stress in &matched.solid().solution().element_stress[0] {
            near(stress, 0.0, 1e-13);
        }
        let old =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let saved = old.solid().solution().clone();
        let free = hot[0].free_linear_strain();
        // Both paths have zero mechanical chord. Only midpoint's actual v1
        // creates an out-of-domain endpoint rate, so that late gate must refuse.
        assert!(
            old.advance_midpoint(
                &hot,
                0.1,
                &affine_targets(free),
                &[],
                settings(),
                midpoint_controls(),
                cx
            )
            .is_err()
        );
        let accepted = old
            .advance(
                &hot,
                0.1,
                &affine_targets(free),
                &[],
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        near(accepted.energy.viscous_dissipation_j, 0.0, 1e-16);
        let other_card = viscous_card("foreign identical coefficients", 4e-6, false);
        let other = [viscous_source(&other_card, &geometry, 300.0, 0.001)];
        assert!(
            old.advance(
                &other,
                0.1,
                &affine_targets(0.0),
                &[],
                settings(),
                controls(),
                cx
            )
            .is_err()
        );
        assert_eq!(old.solid().solution(), &saved);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(old.free_linear_strain_rates_per_s(), &[0.0]);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.accepted_steps(), 0);

        // Prescribed thermal kinematics are meaningful even without eta.
        let elastic_card = super::card("rate-independent thermal input", 4e-6, false);
        let elastic_sources = [source(&elastic_card, &geometry, 300.0, 300.0)];
        let elastic =
            ThermalBrittleTetSolid::try_new(&geometry, &elastic_sources, budget(), cx).unwrap();
        let admitted = ThermalBrittleTetMotion::try_from_solid_with_free_strain_rates(
            &elastic,
            &[[0.0; 3]; 4],
            &[0.123],
            0.0,
            4,
            cx,
        )
        .unwrap();
        assert_eq!(admitted.free_linear_strain_rates_per_s(), &[0.123]);

        gate.request();
        assert!(matches!(
            old.advance(
                &hot,
                0.1,
                &affine_targets(free),
                &[],
                settings(),
                controls(),
                cx
            ),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::Cancelled))
        ));
        assert_eq!(old.solid().solution(), &saved);
        assert_eq!(old.free_linear_strain_rates_per_s(), &[0.0]);
    });
}
