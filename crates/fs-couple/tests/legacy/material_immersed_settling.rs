//! Analytic inertia/Archimedes identities and real native fluid/body exchanges.
//! Synthetic sources test implementation, not calibrated sedimentation accuracy.
use fs_couple::material_immersed::settling::{
    SubmergedEllipsoid, SubmergedMaterialEllipsoids, SubmergedMaterialSpheres, SubmergedSphere,
};
use fs_couple::material_immersed::{BodyStateSi, FilmGeometry, ImmersedDomain};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::{Se3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{E3, Q3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

const DX: f64 = 0.01;
const DT: f64 = 0.1;
const RADIUS: f64 = 0.012;
const G: [f64; 3] = [0.0, -0.0001, 0.0];

fn card(properties: &[(&str, Dims, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, dims, value) in properties {
        let key = if name == "thermal_conductivity" {
            PropertyKey::with_quantity(
                name,
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                )),
            )
        } else {
            PropertyKey::new(name, dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic submerged-solid fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic material".into(),
            phase: "fixture".into(),
            process: "numerical test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
}
fn fluid(rho: f64, mu: f64) -> ResolvedLiquidState {
    resolve_liquid_state(
        &card(&[
            ("density", Density::DIMS, rho),
            ("dynamic_viscosity", DYNAMIC_VISCOSITY_DIMS, mu),
            (
                "specific_heat_capacity",
                SPECIFIC_HEAT_CAPACITY_DIMS,
                1000.0,
            ),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        ]),
        &point(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn body(rho: f64, geometry: VolumeMoments) -> ResolvedBody {
    resolve_body(
        &[RegionBinding {
            key: "ordinary solid",
            geometry,
            material: &card(&[("density", Density::DIMS, rho)]),
            point: &point(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}
fn domain() -> ImmersedDomain {
    ImmersedDomain {
        cells: [8; 3],
        periodic: [true; 3],
        cell_size_m: DX,
        timestep_s: DT,
        maximum_speed_m_s: 0.001,
        maximum_cells: 512,
        samples_per_axis: 4,
    }
}
fn film() -> FilmGeometry {
    FilmGeometry {
        cutoff_gap_m: 0.003,
        roughness_gap_m: 1e-5,
        tangential_gaps_m: None,
        maximum_pairs: 2,
        maximum_resistance_modes: 8,
    }
}
fn rest() -> BodyStateSi {
    BodyStateSi {
        center_m: [0.04; 3],
        orientation: So3::identity(),
        velocity_m_s: [0.0; 3],
        angular_momentum_body_kg_m2_s: [0.0; 3],
    }
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!((a - b).abs() <= absolute, "{a:.17e} != {b:.17e}");
}
type SphereModel<'a, const N: usize> = SubmergedMaterialSpheres<'a, ResolvedLiquidState, N>;
fn populations<const N: usize>(model: &SphereModel<'_, N>) -> Vec<[f64; Q3]> {
    (0..8)
        .flat_map(|z| {
            (0..8).flat_map(move |y| (0..8).map(move |x| model.grid().populations(x, y, z)))
        })
        .collect()
}
fn represented_momentum<const N: usize>(model: &SphereModel<'_, N>) -> [f64; 3] {
    let mut total = [0.0; 3];
    let unit = model.scaling().transport().density_kg_m3() * DX.powi(4) / DT;
    for values in populations(model) {
        for (q, value) in values.into_iter().enumerate() {
            let e = [E3[q].0, E3[q].1, E3[q].2];
            for i in 0..3 {
                total[i] += value * f64::from(e[i]) * unit;
            }
        }
    }
    for (carrier, split) in model.excess_carriers_lattice().iter().zip(model.inertias()) {
        for i in 0..3 {
            total[i] += split.excess.mass() * carrier.velocity[i] * DX / DT;
        }
    }
    total
}

#[test]
fn exact_source_inertia_subtraction_drives_real_fluid_and_reduced_gravity_once() {
    let gate = CancelGate::new();
    let sphere = VolumeMoments::sphere(RADIUS).unwrap();
    let mut terminal_samples = Vec::new();
    for (rho_s, rho_f) in [(3000.0, 1000.0), (1500.0, 1000.0), (3000.0, 1200.0)] {
        let liquid = fluid(rho_f, 0.11);
        let solid = body(rho_s, sphere);
        let source = SubmergedSphere::new(&solid, RADIUS).unwrap();
        let mut model =
            SphereModel::new(&liquid, domain(), [source], [rest()], G, film(), &gate).unwrap();
        let split = model.inertias()[0];
        assert_eq!(split.physical, solid.inertia());
        assert!(core::ptr::eq(model.sources()[0].material(), &solid));
        assert!(core::ptr::eq(model.scaling().material(), &liquid));
        close(split.displaced.mass(), rho_f * sphere.volume_m3(), 1e-16);
        close(
            split.excess.mass(),
            (rho_s - rho_f) * sphere.volume_m3(),
            1e-16,
        );
        for i in 0..9 {
            close(
                split.excess.inertia_com().m[i],
                (rho_s - rho_f) * sphere.inertia_m5().m[i],
                1e-20,
            );
        }
        let initial = populations(&model);
        for _ in 0..5 {
            let before = represented_momentum(&model);
            let receipt = model.step_forces_si([[0.0; 6]], [0.0; 3], &gate).unwrap();
            let after = represented_momentum(&model);
            let account = receipt.body_exchange_si(0).unwrap();
            assert_eq!(receipt.native.coupling.forcing.impulse, [0.0; 3]);
            for i in 0..3 {
                close(after[i] - before[i], split.excess.mass() * G[i] * DT, 2e-16);
                close(
                    account.buoyancy_impulse_ns[i],
                    -rho_f * sphere.volume_m3() * G[i] * DT,
                    1e-20,
                );
                close(account.dry_momentum_residual[i], 0.0, 2e-16);
                close(
                    receipt.native.coupling.bodies[0].hydrodynamic.body_impulse[i]
                        + receipt.native.coupling.bodies[0]
                            .hydrodynamic
                            .fluid_population_impulse[i],
                    0.0,
                    1e-14,
                );
            }
            close(account.represented_work_residual_j, 0.0, 1e-16);
            assert!(account.reduced_gravity_potential_change_j <= 0.0);
        }
        assert_ne!(populations(&model), initial);
        assert_eq!(model.accepted_steps(), 5);
        let state = model.state_si(0).unwrap();
        assert!(state.velocity_m_s[1] < 0.0);
        assert!(state.center_m[1] < rest().center_m[1]);
        let before = populations(&model);
        let carriers = *model.excess_carriers_lattice();
        let diagnosis = model.interior_diagnostics(&gate).unwrap()[0];
        assert!(diagnosis.rms_slip_m_s.is_finite());
        assert!(diagnosis.maximum_slip_m_s >= diagnosis.rms_slip_m_s - 1e-18);
        assert!(diagnosis.relative_volume_error < 0.15);
        assert_eq!(before, populations(&model));
        assert_eq!(carriers, *model.excess_carriers_lattice());
        terminal_samples.push(state.velocity_m_s[1]);
    }
    assert!(terminal_samples[0] < terminal_samples[1]); // More density contrast settles faster here.
    assert_ne!(terminal_samples[0], terminal_samples[2]); // Original fluid density reaches dynamics.
}

#[test]
fn rotated_ellipsoid_preserves_full_source_tensor_and_reconstructs_dry_spin() {
    let gate = CancelGate::new();
    let liquid = fluid(1000.0, 0.11);
    let axes = [0.014, 0.009, 0.006];
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap();
    let geometry = VolumeMoments::ellipsoid(Vec3::new(axes[0], axes[1], axes[2]))
        .unwrap()
        .transformed(Se3::from_parts(rotation, Vec3::new(0.0, 0.0, 0.0)).unwrap())
        .unwrap();
    let solid = body(2500.0, geometry);
    let source = SubmergedEllipsoid::new(&solid, axes, rotation).unwrap();
    let mut model =
        SubmergedMaterialEllipsoids::new(&liquid, domain(), [source], [rest()], G, &gate).unwrap();
    let split = model.inertias()[0];
    assert!(split.physical.inertia_com().m[1].abs() > 1e-10);
    for i in 0..9 {
        close(
            split.excess.inertia_com().m[i],
            solid.inertia().inertia_com().m[i] - 1000.0 * geometry.inertia_m5().m[i],
            1e-20,
        );
    }
    for _ in 0..3 {
        let step = model
            .step_forces_si([[0.0, 0.0, 0.0, 1e-10, -2e-10, 3e-10]], [0.0; 3], &gate)
            .unwrap();
        let account = step.body_exchange_si(0).unwrap();
        for value in account.dry_momentum_residual {
            close(value, 0.0, 1e-16);
        }
        close(account.represented_work_residual_j, 0.0, 1e-16);
        let state = model.state_si(0).unwrap();
        let carrier = model.excess_carriers_lattice()[0];
        let omega = carrier
            .inertia
            .velocity(carrier.angular_momentum_body)
            .unwrap()
            .scale(1.0 / DT);
        let expected = solid.inertia().inertia_com().apply(omega);
        for (a, b) in state
            .angular_momentum_body_kg_m2_s
            .into_iter()
            .zip([expected.x, expected.y, expected.z])
        {
            close(a, b, 1e-20);
        }
    }
    assert_eq!(model.accepted_steps(), 3);
    assert!(model.interior_diagnostics(&gate).unwrap()[0].sampled_volume_m3 > 0.0);
}

#[test]
fn simultaneous_settling_keeps_existing_film_reciprocity_and_source_refresh() {
    let gate = CancelGate::new();
    let liquid = fluid(1000.0, 0.11);
    let more_viscous = fluid(1000.0, 0.22);
    let light = body(1500.0, VolumeMoments::sphere(RADIUS).unwrap());
    let heavy = body(3000.0, VolumeMoments::sphere(RADIUS).unwrap());
    let mut initial = [rest(); 2];
    initial[0].center_m[1] = 0.0275;
    initial[1].center_m[1] = 0.0525;
    let mut model = SphereModel::new(
        &liquid,
        domain(),
        [
            SubmergedSphere::new(&light, RADIUS).unwrap(),
            SubmergedSphere::new(&heavy, RADIUS).unwrap(),
        ],
        initial,
        G,
        film(),
        &gate,
    )
    .unwrap();
    let receipt = model
        .step_forces_si([[0.0; 6]; 2], [0.0; 3], &gate)
        .unwrap();
    assert_eq!(receipt.native.coupling.lubrication.active_pairs, 1);
    assert!(receipt.native.coupling.lubrication.dissipated_energy >= 0.0);
    for i in 0..3 {
        close(
            receipt.native.coupling.lubrication.body_impulses[0][i]
                + receipt.native.coupling.lubrication.body_impulses[1][i],
            0.0,
            1e-14,
        );
    }
    for i in 0..2 {
        close(
            receipt
                .body_exchange_si(i)
                .unwrap()
                .represented_work_residual_j,
            0.0,
            1e-16,
        );
    }
    let before = populations(&model);
    let carriers = *model.excess_carriers_lattice();
    let splits = *model.inertias();
    model.rebind_fluid(&more_viscous, &gate).unwrap();
    assert_eq!(before, populations(&model));
    assert_eq!(carriers, *model.excess_carriers_lattice());
    assert_eq!(splits, *model.inertias());
    close(
        model.grid().viscosity(),
        0.22 / 1000.0 * DT / DX.powi(2),
        1e-16,
    );
    assert!(core::ptr::eq(model.scaling().material(), &more_viscous));
}

#[test]
fn light_neutral_and_moving_initial_states_refuse_and_native_failure_is_atomic() {
    let gate = CancelGate::new();
    let liquid = fluid(1000.0, 0.11);
    for density in [900.0, 1000.0] {
        let material = body(density, VolumeMoments::sphere(RADIUS).unwrap());
        assert!(
            SphereModel::new(
                &liquid,
                domain(),
                [SubmergedSphere::new(&material, RADIUS).unwrap()],
                [rest()],
                G,
                film(),
                &gate
            )
            .is_err()
        );
    }
    let material = body(3000.0, VolumeMoments::sphere(RADIUS).unwrap());
    let source = SubmergedSphere::new(&material, RADIUS).unwrap();
    let mut moving = rest();
    moving.velocity_m_s[0] = 1e-8;
    assert!(SphereModel::new(&liquid, domain(), [source], [moving], G, film(), &gate).is_err());
    let mut spinning = rest();
    spinning.angular_momentum_body_kg_m2_s[0] = 1e-12;
    assert!(SphereModel::new(&liquid, domain(), [source], [spinning], G, film(), &gate).is_err());
    let mut model =
        SphereModel::new(&liquid, domain(), [source], [rest()], G, film(), &gate).unwrap();
    let mut replay =
        SphereModel::new(&liquid, domain(), [source], [rest()], G, film(), &gate).unwrap();
    let before = populations(&model);
    let carriers = *model.excess_carriers_lattice();
    assert!(
        model
            .step_forces_si([[1e100, 0.0, 0.0, 0.0, 0.0, 0.0]], [0.0; 3], &gate)
            .is_err()
    );
    assert_eq!(before, populations(&model));
    assert_eq!(carriers, *model.excess_carriers_lattice());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        model
            .step_forces_si([[0.0; 6]], [0.0; 3], &cancelled)
            .is_err()
    );
    assert_eq!(model.accepted_steps(), 0);
    assert_eq!(before, populations(&model));
    let different_density = fluid(1100.0, 0.11);
    assert!(model.rebind_fluid(&different_density, &gate).is_err());
    assert_eq!(before, populations(&model));
    assert_eq!(carriers, *model.excess_carriers_lattice());
    let a = model.step_forces_si([[0.0; 6]], [0.0; 3], &gate).unwrap();
    let b = replay.step_forces_si([[0.0; 6]], [0.0; 3], &gate).unwrap();
    assert_eq!(a.native, b.native);
    assert_eq!(populations(&model), populations(&replay));
}
