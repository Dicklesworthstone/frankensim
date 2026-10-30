//! G1/G3/G4: native material-temperature updates on an accepted rod trajectory.
//! Manufactured source curves test the physical equations, not measured alloys.
use super::*;
use fs_couple::rod_motion::SpatialRodMotion;
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::INVERSE_TEMPERATURE_DIMS;
use fs_qty::Temperature;
use fs_solid::rod::material_field::gravity::ReferenceDisplacementVolume;
use fs_solid::rod::material_field::reference::RodReferenceGeometry;
use fs_solid::rod::material_field::thermal::source::{
    RodLayerSource, RodMaterialSource, RodSegmentSource,
};

const LENGTH: f64 = 0.8;
const RADIUS: f64 = 0.05;

fn card(name: &str, elastic_scale: f64, expansion_scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, lo, hi, dims) in [
        ("density", 1000.0, 500.0, Density::DIMS),
        (
            "young_modulus",
            1e4 * elastic_scale,
            5e3 * elastic_scale,
            Pressure::DIMS,
        ),
        (
            "shear_modulus",
            4e3 * elastic_scale,
            1e3 * elastic_scale,
            Pressure::DIMS,
        ),
        (
            "linear_thermal_expansion_coefficient",
            1e-5 * expansion_scale,
            3e-5 * expansion_scale,
            INVERSE_TEMPERATURE_DIMS,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, lo), (400.0, hi)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured continuing thermoelastic rod".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn source<'a>(cards: &[&'a MaterialCard], segments: usize, limit: f64) -> RodMaterialSource<'a> {
    RodMaterialSource::try_new(
        LENGTH,
        (0..segments)
            .map(|_| RodSegmentSource {
                inner_radius_m: 0.0,
                shear_correction: 0.8,
                layers: cards
                    .iter()
                    .enumerate()
                    .map(|(i, card)| RodLayerSource {
                        outer_radius_m: RADIUS * (i + 1) as f64 / cards.len() as f64,
                        card,
                        reference_point: QueryPoint::new().with("T", 300.0).unwrap(),
                        mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
                        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
                    })
                    .collect(),
            })
            .collect(),
        "shear_modulus",
        limit,
        &CancelGate::new(),
    )
    .unwrap()
}

fn axial(source: &RodMaterialSource<'_>, strain: f64) -> RodMotion {
    let gate = CancelGate::new();
    let specimen = source
        .resolve_at(source.reference_temperatures_k(), &gate)
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [None, Some(0.0), Some(0.0)],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
    let model = Arc::new(
        specimen
            .into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] = LENGTH * (1.0 + strain);
    RodMotion::new(model, pose, &[[0.0; 3]; 2], &[[0.0; 3]; 2], 0.0, &gate).unwrap()
}

fn alpha_integral(t: f64, scale: f64) -> f64 {
    let d = t - 300.0;
    scale * (1e-5 * d + 1e-7 * d * d)
}
fn young(t: f64, scale: f64) -> f64 {
    scale * (1e4 - 50.0 * (t - 300.0))
}
fn shear(t: f64, scale: f64) -> f64 {
    scale * (4e3 - 30.0 * (t - 300.0))
}

#[test]
fn g1_source_temperature_changes_transient_extension_and_frequency_without_reset() {
    let gate = CancelGate::new();
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    for (elastic_scale, expansion_scale) in [(1.0, 1.0), (2.0, 1.0), (1.0, 2.0)] {
        let card = card(
            "data-only transient material",
            elastic_scale,
            expansion_scale,
        );
        let source = source(&[&card], 1, 0.1);
        let initial = evolve(axial(&source, 0.01), 0.002, 30);
        let before = initial.observe([0.0; 3]).unwrap();
        let (hot, refresh) = initial
            .rebind_material_temperatures_candidate(&source, &[400.0], &gate)
            .unwrap();
        assert_eq!(hot.pose(), initial.pose());
        assert_eq!(hot.time_s(), initial.time_s());
        assert_eq!(hot.accepted_steps(), initial.accepted_steps());
        assert_eq!(hot.linear_momenta_kg_m_s(), initial.linear_momenta_kg_m_s());
        assert_eq!(hot.body_spins_kg_m2_s(), initial.body_spins_kg_m2_s());
        assert_eq!(hot.model().inertia(), initial.model().inertia());
        assert_eq!(refresh.temperatures_before_k, [300.0]);
        assert_eq!(refresh.temperatures_after_k, [400.0]);
        assert_eq!(
            hot.observe([0.0; 3]).unwrap().kinetic_energy_j,
            before.kinetic_energy_j
        );
        let current = &hot.model().source().sections()[0].elastic().layers()[0].material;
        near(current.property("density").unwrap().value_si(), 500.0, 0.0);
        near(
            hot.model().inertia()[1].mass_kg(),
            0.5 * 1000.0 * area * LENGTH,
            1e-14,
        );

        let x = initial.pose().positions_m[1][0];
        let strain = x / LENGTH - 1.0;
        let free = alpha_integral(400.0, expansion_scale);
        let old_storage = 0.5 * LENGTH * young(300.0, elastic_scale) * area * strain * strain;
        let new_storage =
            0.5 * LENGTH * young(400.0, elastic_scale) * area * (strain - free).powi(2);
        near(refresh.strain_energy_before_j, old_storage, 2e-14);
        near(refresh.strain_energy_after_j, new_storage, 2e-14);
        near(
            refresh.thermal_parameter_work_j,
            new_storage - old_storage,
            3e-14,
        );

        // Original reference density sets inertia; current rho(T) must not double
        // omega by silently reducing mass on unchanged section geometry.
        let omega =
            (young(400.0, elastic_scale) * area / LENGTH / hot.model().inertia()[1].mass_kg())
                .sqrt();
        let equilibrium = LENGTH * (1.0 + free);
        let velocity = hot.velocities_m_s().unwrap()[1][0];
        let duration = 0.2;
        let exact = equilibrium
            + (x - equilibrium) * fs_math::det::cos(omega * duration)
            + velocity / omega * fs_math::det::sin(omega * duration);
        let mut errors = Vec::new();
        for (dt, n) in [(0.004, 50), (0.002, 100), (0.001, 200)] {
            let after = evolve(hot.clone(), dt, n);
            errors.push((after.pose().positions_m[1][0] - exact).abs());
        }
        for i in 0..2 {
            assert!(
                (3.8..4.2).contains(&(errors[i] / errors[i + 1])),
                "{errors:?}"
            );
        }
        // The same accepted pose can cool back to its original coefficients.
        let (cold, back) = hot
            .rebind_material_temperatures_candidate(&source, &[300.0], &gate)
            .unwrap();
        assert_eq!(
            cold.model().source().sections(),
            initial.model().source().sections()
        );
        assert_eq!(cold.observe([0.0; 3]).unwrap(), before);
        near(
            back.thermal_parameter_work_j,
            -refresh.thermal_parameter_work_j,
            1e-15,
        );
    }
}

#[test]
fn g1_nonuniform_layer_temperatures_drive_native_support_forces_and_torsion() {
    let gate = CancelGate::new();
    let core = card("core", 1.0, 1.0);
    let wrap = card("wrap", 2.0, 2.0);
    let source = source(&[&core, &wrap], 2, 0.1);
    let angle = 0.01;
    let last_frame = So3::exp(So3Tangent::new(Vec3::new(angle, 0.0, 0.0))).unwrap();
    let specimen = source
        .resolve_at(&[300.0; 4], &gate)
        .unwrap()
        .with_supports(
            &[
                RodSupport {
                    node: 1,
                    position_m: [Some(0.4), Some(0.0), Some(0.0)],
                    frame: Some(So3::identity()),
                },
                RodSupport {
                    node: 2,
                    position_m: [Some(0.8), Some(0.0), Some(0.0)],
                    frame: Some(last_frame),
                },
            ],
            &gate,
        )
        .unwrap();
    let model = Arc::new(
        specimen
            .into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let cold = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3]; 3],
        &[[0.0; 3]; 3],
        0.7,
        &gate,
    )
    .unwrap();
    let temperatures = [350.0, 400.0, 325.0, 375.0];
    let (hot, refresh) = cold
        .rebind_material_temperatures_candidate(&source, &temperatures, &gate)
        .unwrap();
    let mut axial = [0.0; 2];
    let mut twist_stiffness = [0.0; 2];
    let mut stored = 0.0;
    for segment in 0..2 {
        for (layer, scale) in [1.0, 2.0].into_iter().enumerate() {
            let inner = RADIUS * layer as f64 / 2.0;
            let outer = RADIUS * (layer + 1) as f64 / 2.0;
            let area = core::f64::consts::PI * (outer * outer - inner * inner);
            let polar = 0.5 * core::f64::consts::PI * (outer.powi(4) - inner.powi(4));
            let t = temperatures[2 * segment + layer];
            let expansion = alpha_integral(t, scale);
            axial[segment] -= young(t, scale) * area * expansion;
            twist_stiffness[segment] += shear(t, scale) * polar;
            stored += 0.5 * 0.4 * young(t, scale) * area * expansion * expansion;
        }
    }
    let torque = twist_stiffness[1] * angle / 0.4;
    stored += 0.5 * twist_stiffness[1] * angle * angle / 0.4;
    near(refresh.strain_energy_after_j, stored, 2e-14);
    let response = hot.model().evaluate(hot.pose(), &gate).unwrap();
    let forces = [axial[0], axial[1] - axial[0], -axial[1]];
    let torques = [0.0, torque, -torque];
    for node in 0..3 {
        near(response.force_n[node][0], forces[node], 2e-12);
        near(response.torque_body_nm[node][0], torques[node], 2e-12);
    }
    let dt = 0.001;
    let (after, step) = hot
        .step_candidate(dt, &[RodNodalWrench::default(); 3], controls(dt), &gate)
        .unwrap();
    assert_eq!(after.pose(), hot.pose());
    for node in 0..3 {
        near(
            step.support_impulses[node].linear_n_s[0],
            -dt * forces[node],
            2e-14,
        );
        near(
            step.support_impulses[node].angular_nm_s[0],
            -dt * torques[node],
            2e-14,
        );
    }
    near(step.integration_defect_j, 0.0, 1e-15);
    near(step.source_work_quadrature_j, 0.0, 0.0);
}

#[test]
fn g3_curved_eccentric_inventory_and_body_loads_survive_temperature_continuation() {
    let gate = CancelGate::new();
    let card = card("curved core", 1.0, 1.0);
    let geometry = RodReferenceGeometry::try_new(
        LENGTH,
        vec![[0.0; 3], [0.39, 0.02, 0.0], [0.77, 0.06, 0.01]],
        [0.0, 0.04, 0.08]
            .into_iter()
            .map(|a| So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, a))).unwrap())
            .collect(),
    )
    .unwrap();
    let source = source(&[&card], 2, 0.1)
        .with_reference_geometry(geometry, &gate)
        .unwrap();
    let rod = source.resolve_at(&[300.0; 2], &gate).unwrap();
    // An axial component gives this disk-plus-point nodal inventory genuine
    // three-dimensional support, away from the planar inertia triangle limit.
    // Retain off-diagonal products and the original strict principal admission.
    let link = MaterialTensionLink::try_new(
        rod.sections()[0].elastic().clone(),
        [
            TensionLinkAttachment {
                node: 0,
                offset_body_m: [0.006, 0.01, 0.005],
            },
            TensionLinkAttachment {
                node: 2,
                offset_body_m: [0.006, 0.01, 0.005],
            },
        ],
        0.78,
        1,
        0.2,
    )
    .unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let rod = rod
        .with_tension_links(vec![link], &gate)
        .unwrap()
        .with_gravity([0.0, 0.0, -9.0], &gate)
        .unwrap()
        .with_submerged_fluid(&gas, ReferenceDisplacementVolume::FloodedVoids, &gate)
        .unwrap();
    let model = Arc::new(
        rod.into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let initial = SpatialRodMotion::new_spatial(
        model.clone(),
        model.initial_pose(),
        &[[0.1, 0.02, 0.0]; 3],
        &[[0.1, 0.02, 0.01]; 3],
        0.0,
        1e-12,
        &gate,
    )
    .unwrap();
    let (initial, _) = initial
        .step_candidate(1e-4, &[RodNodalWrench::default(); 3], controls(1e-4), &gate)
        .unwrap();
    let before = initial.observe([0.0; 3]).unwrap();
    let loads = initial.model().body_loads(initial.pose(), &gate).unwrap();
    let (hot, _) = initial
        .rebind_material_temperatures_candidate(&source, &[350.0, 375.0], &gate)
        .unwrap();
    assert_eq!(hot.pose(), initial.pose());
    assert_eq!(hot.linear_momenta_kg_m_s(), initial.linear_momenta_kg_m_s());
    assert_eq!(hot.body_spins_kg_m2_s(), initial.body_spins_kg_m2_s());
    assert_eq!(hot.model().inertia(), initial.model().inertia());
    assert_eq!(
        hot.model().source().reference_mass(),
        initial.model().source().reference_mass()
    );
    assert_eq!(
        hot.model().source().tension_links(),
        initial.model().source().tension_links()
    );
    assert_eq!(hot.model().body_loads(hot.pose(), &gate).unwrap(), loads);
    assert_eq!(hot.time_s(), initial.time_s());
    assert_eq!(hot.accepted_steps(), initial.accepted_steps());
    let observed = hot.observe([0.0; 3]).unwrap();
    assert_eq!(observed.kinetic_energy_j, before.kinetic_energy_j);
    assert_eq!(
        observed.linear_momentum_kg_m_s,
        before.linear_momentum_kg_m_s
    );
    assert_eq!(
        observed.angular_momentum_kg_m2_s,
        before.angular_momentum_kg_m2_s
    );
    let (_, step) = hot
        .step_candidate(1e-4, &[RodNodalWrench::default(); 3], controls(1e-4), &gate)
        .unwrap();
    for c in 0..3 {
        near(step.linear_balance_residual_n_s[c], 0.0, 2e-12);
        near(step.angular_balance_residual_nm_s[c], 0.0, 2e-12);
    }
}

#[test]
fn g4_foreign_sources_references_actual_strain_and_cancellation_preserve_history() {
    let gate = CancelGate::new();
    let card = card("original", 1.0, 1.0);
    let source = source(&[&card], 1, 0.015);
    let state = axial(&source, -0.014);
    let before = state.observe([0.0; 3]).unwrap();
    // New thermal strain is 0.002; actual elastic compression would be -0.016.
    // The source's reference pose is admissible, but THIS accepted pose is not.
    assert!(source.resolve_at(&[400.0], &gate).is_ok());
    assert!(
        state
            .rebind_material_temperatures_candidate(&source, &[400.0], &gate)
            .is_err()
    );
    for temperatures in [&[][..], &[f64::NAN][..], &[401.0][..]] {
        assert!(
            state
                .rebind_material_temperatures_candidate(&source, temperatures, &gate)
                .is_err()
        );
    }
    let foreign = self::card("other source", 1.0, 1.0);
    let foreign = self::source(&[&foreign], 1, 0.015);
    assert!(
        state
            .rebind_material_temperatures_candidate(&foreign, &[305.0], &gate)
            .is_err()
    );
    let changed_limit = self::source(&[&card], 1, 0.1);
    assert!(
        state
            .rebind_material_temperatures_candidate(&changed_limit, &[305.0], &gate)
            .is_err()
    );
    let explicit_geometry = RodReferenceGeometry::try_new(
        LENGTH,
        vec![[0.0; 3], [LENGTH, 0.0, 0.0]],
        vec![So3::identity(); 2],
    )
    .unwrap();
    let changed_reference = source
        .with_reference_geometry(explicit_geometry, &gate)
        .unwrap();
    assert!(
        state
            .rebind_material_temperatures_candidate(&changed_reference, &[305.0], &gate)
            .is_err()
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        state
            .rebind_material_temperatures_candidate(&source, &[305.0], &cancelled)
            .is_err()
    );
    assert_eq!(state.observe([0.0; 3]).unwrap(), before);
    assert_eq!(state.model().section_temperatures_k().unwrap(), [300.0]);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.accepted_steps(), 0);

    // An unheated reference can explicitly acquire its original thermal bind,
    // even though reference and resolved receipts use distinct selection modes.
    let reference = Arc::new(
        source
            .reference_rod()
            .clone()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let unbound = RodMotion::new(
        reference.clone(),
        reference.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let (bound, _) = unbound
        .rebind_material_temperatures_candidate(&source, &[350.0], &gate)
        .unwrap();
    assert!(unbound.model().source().thermal_sections().is_none());
    assert!(bound.model().source().thermal_sections().is_some());
    assert_eq!(bound.pose(), unbound.pose());
}
