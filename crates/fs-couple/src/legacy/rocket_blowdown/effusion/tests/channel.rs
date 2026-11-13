//! G1/G3/G4: resolved channel recoil, actual finite wall storage and atomic refusal.
mod conduction;
mod trajectory;

use super::super::channel::*;
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    thermal::transport::volumetric::{
        PreparedVolumetricThermalTransport, VOLUMETRIC_HEAT_CAPACITY_DIMS,
        VolumetricThermalPathQuery, timeline::FixedVolumeThermalLaw,
    },
    thermochemical_gas::energy::parcel::transfer::species::{
        effusion::channel::{
            ClausingControls, CylindricalMolecularGeometry, DiffuseMolecularChannel,
            momentum::DiffuseChannelMomentumInterface,
        },
        network::transport::effusion::channel::thermal::FixedVolumeChannelWall,
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

// The actual material owner consumes an authored Cv per volume, A+B*T.
// This is a reversible fixed-volume synthetic source, not an inferred solid Cp.
fn wall_material() -> PreparedVolumetricThermalTransport {
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = |t| {
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                t,
            )
            .unwrap()
    };
    let capacity = PropertyKey::new(
        "isochoric-volumetric-capacity",
        VOLUMETRIC_HEAT_CAPACITY_DIMS,
    );
    let conductivity = PropertyKey::with_quantity(
        "thermal-conductivity",
        quantity(QuantityKind::ThermalConductivity),
    );
    let span = [200.0, 2000.0];
    let mut claims = ClaimSet::new();
    for (key, values, dims) in [
        (
            capacity.clone(),
            span.map(|t| (t, 1e6 + 1000.0 * t)),
            VOLUMETRIC_HEAT_CAPACITY_DIMS,
        ),
        (
            conductivity.clone(),
            span.map(|t| (t, 1.0)),
            Dims([1, 1, -3, -1, 0, 0]),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: values.to_vec(),
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    span[0],
                    span[1],
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic reversible isochoric channel wall".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored finite channel wall".into(),
            phase: "solid".into(),
            process: "fixed volume and composition".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    PreparedVolumetricThermalTransport::try_new(
        &card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &capacity,
            conductivity: &conductivity,
            temperature_axis: "temperature",
            lower: &point(span[0]),
            upper: &point(span[1]),
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 2],
        },
    )
    .unwrap()
}

fn channel(before: &GasRocketState, aspect: f64, alpha: f64) -> DiffuseChannelMomentumInterface {
    let path = DiffuseMolecularChannel::try_new(
        interface(before.reservoir().model()).interface().clone(),
        CylindricalMolecularGeometry::try_new(
            1e-9,
            aspect * 1e-9,
            1e-7,
            "synthetic diffuse circular recoil passage",
        )
        .unwrap(),
    )
    .unwrap()
    .with_thermal_accommodation(
        alpha,
        ClausingControls::default(),
        "explicit synthetic molecular reset law",
        || false,
    )
    .unwrap();
    DiffuseChannelMomentumInterface::try_new(
        path,
        ClausingControls::default(),
        "diffuse Maxwellian reset or elastic speed law",
        || false,
    )
    .unwrap()
}

fn channel_outlet(
    interface: &DiffuseChannelMomentumInterface,
    normal: Vec3,
) -> EffusionChannelOutlet<'_> {
    EffusionChannelOutlet {
        id: 4,
        interface,
        area_m2: AREA,
        wall: 9,
        outward_normal_world: normal,
    }
}

fn environment() -> EffusionRocketEnvironment {
    EffusionRocketEnvironment {
        inertial_frame: FRAME,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
    }
}

fn check_channel_accounts(before: &GasRocketState, candidate: &EffusionChannelCandidate<'_, '_>) {
    let after = candidate.endpoint();
    let step = candidate.step();
    let r = &step.recoil;
    vector_near(
        after
            .vehicle()
            .velocity_m_s
            .scale(mass(after))
            .sub(before.vehicle().velocity_m_s.scale(mass(before)))
            .add(r.exhaust_momentum_world_kg_m_s)
            .sub(r.exhaust.external_impulse_n_s),
        Vec3::ZERO,
        2e-12,
    );
    near(
        internal(after) - internal(before) + step.thermal.wall_energy_change_j + kinetic(after)
            - kinetic(before)
            + r.exhaust_energy_world_j
            - r.exhaust.external_work_j,
        0.0,
        3e-8,
    );
    assert!(
        r.absolute_energy_residual_sum_j + r.energy_uncertainty_bound_j
            <= controls().transport.exchange.energy_tolerance_j
    );
    vector_near(
        r.momentum_residual_kg_m_s,
        Vec3::ZERO,
        controls().momentum_tolerance_n_s,
    );
    assert_eq!(
        after.vehicle().propellant_mass_kg,
        after.reservoir().mass_kg()
    );
    assert_eq!(after.vehicle().dry_mass_kg, before.vehicle().dry_mass_kg);
    assert_eq!(after.reservoir().volume_m3(), VOLUME);
}

#[test]
fn hot_and_cold_channels_use_collimated_escape_and_finite_nonlinear_wall_energy_once() {
    let before = initial(PURE, Vec3::ZERO);
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let channel = channel(&before, 2.0, 1.0);
    let normal = Vec3::new(0.6, 0.8, 0.0);
    let geometry = channel.channel().geometry();
    let tau = geometry.transmission_probability();
    let direct = geometry.direct_probability();
    let n = escaped_moles(0.028, 1.0);
    let n0 = before.reservoir().mass_kg() / 0.028;
    let incoming_speed = (core::f64::consts::PI * R * T0 / (2.0 * 0.028)).sqrt();
    for tw in [250.0, 900.0] {
        let walls = BTreeMap::from([(
            9,
            FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, tw).unwrap(),
        )]);
        let result = before
            .advance_effusion_channels(
                environment(),
                &walls,
                &[channel_outlet(&channel, normal)],
                DT,
                endpoint(),
                controls(),
                || false,
            )
            .unwrap();
        let wall_speed = incoming_speed * (tw / T0).sqrt();
        let impulse = 1.5
            * n
            * 0.028
            * (incoming_speed * channel.moments().direct_axial_cosine()
                + wall_speed * channel.moments().wall_transmitted_axial_cosine());
        let qwall = -2.0 * R * n * (1.0 - direct) * (tw - T0);
        let escape = 2.0 * R * n * (direct * T0 + (tau - direct) * tw);
        let gas_energy = -2.0 * R * n * tau * T0 + 2.0 * R * n * (1.0 - tau) * (tw - T0);
        let gas_t = (1.5 * R * n0 * T0 + gas_energy) / (1.5 * R * (n0 - tau * n));
        let c0 = 1e6 + 1000.0 * tw;
        let tw1 = tw + (2.0 * qwall / 1e-6) / ((c0 * c0 + 2.0 * 1000.0 * qwall / 1e-6).sqrt() + c0);
        near(
            result.step().recoil.expelled_mass_kg,
            tau * n * 0.028,
            1e-18,
        );
        near(result.step().recoil.relative_energy_j, escape, 2e-10);
        near(result.step().thermal.wall_energy_change_j, qwall, 1e-9);
        near(result.walls()[&9].temperature_k(), tw1, 2e-10);
        near(
            result.endpoint().reservoir().state().state().temperature,
            gas_t,
            3e-8,
        );
        vector_near(
            result.step().recoil.relative_momentum_kg_m_s,
            normal.scale(impulse),
            2e-16,
        );
        vector_near(
            result.step().recoil.exhaust.mechanics.thrust_impulse_n_s,
            normal.scale(-impulse),
            2e-16,
        );
        let native = &result.step().momentum[0].channel;
        assert!(native.wall_impulse_n_s.abs() > 1e-6);
        near(
            native.exterior_impulses_n_s[0] + native.wall_impulse_n_s,
            -impulse,
            2e-16,
        );
        let speed = impulse / (tau * n * 0.028);
        vector_near(
            result.endpoint().vehicle().velocity_m_s,
            normal.scale(-speed * (mass(&before) / mass(result.endpoint())).ln()),
            1e-11,
        );
        if tw > T0 {
            assert!(
                impulse > tau * AREA * P0 * DT / 2.0,
                "conductance times thin-aperture pressure misses the channel exit momentum"
            );
        }
        check_channel_accounts(&before, &result);
        let next = result
            .endpoint()
            .advance_effusion_channels(
                environment(),
                result.walls(),
                &[channel_outlet(&channel, normal)],
                DT,
                endpoint(),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(
            next.step().thermal.gas.outlets[0].wall_temperature_k,
            result.walls()[&9].temperature_k()
        );
        check_channel_accounts(result.endpoint(), &next);
        assert_eq!(walls[&9].mass_kg(), result.walls()[&9].mass_kg());
        assert_eq!(result.walls()[&9].law().identity(), law.identity());
    }
}

#[test]
fn opposed_channels_share_one_wall_and_keep_species_heat_without_net_thrust() {
    let before = initial(MIXTURE, Vec3::new(20.0, -10.0, 3.0));
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let channel = channel(&before, 1.0, 1.0);
    let a = channel_outlet(&channel, Vec3::new(1.0, 0.0, 0.0));
    let b = EffusionChannelOutlet {
        id: 1,
        outward_normal_world: Vec3::new(-1.0, 0.0, 0.0),
        ..a
    };
    let result = before
        .advance_effusion_channels(
            environment(),
            &walls,
            &[a, b],
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    let reordered = before
        .advance_effusion_channels(
            environment(),
            &walls,
            &[b, a],
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result.endpoint(), reordered.endpoint());
    assert_eq!(result.step(), reordered.step());
    assert_eq!(result.step().thermal.walls.len(), 1);
    assert_eq!(result.step().momentum[0].channel.id, 1);
    assert_eq!(result.step().recoil.relative_momentum_kg_m_s, Vec3::ZERO);
    assert_eq!(
        result.endpoint().vehicle().velocity_m_s,
        before.vehicle().velocity_m_s
    );
    assert!(result.walls()[&9].temperature_k() < walls[&9].temperature_k());
    let species = &result.step().thermal.gas.outlets[0].escaped_species;
    near(
        (species[1].mass_kg / 0.004) / (species[0].mass_kg / 0.028),
        7.0_f64.sqrt(),
        2e-14,
    );
    assert!(
        result.step().thermal.gas.outlets[0]
            .wall
            .reflected_species
            .iter()
            .all(|s| s.mass_kg > 0.0)
    );
    check_channel_accounts(&before, &result);
}

#[test]
fn loaded_channel_recoil_transforms_complete_escape_energy_under_a_galilean_boost() {
    let before = initial(PURE, Vec3::new(10.0, 20.0, -5.0));
    let boost = Vec3::new(200.0, -30.0, 10.0);
    let boosted = initial(PURE, before.vehicle().velocity_m_s.add(boost));
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let channel = channel(&before, 2.0, 1.0);
    let ports = [channel_outlet(&channel, Vec3::new(0.0, 1.0, 0.0))];
    let a = before
        .advance_effusion_channels(
            loaded_environment(),
            &walls,
            &ports,
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    let b = boosted
        .advance_effusion_channels(
            loaded_environment(),
            &walls,
            &ports,
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(a.endpoint().reservoir(), b.endpoint().reservoir());
    assert_eq!(a.step().thermal, b.step().thermal);
    let (ar, br) = (&a.step().recoil, &b.step().recoil);
    vector_near(
        b.endpoint()
            .vehicle()
            .velocity_m_s
            .sub(a.endpoint().vehicle().velocity_m_s),
        boost,
        1e-12,
    );
    vector_near(
        br.exhaust_momentum_world_kg_m_s
            .sub(ar.exhaust_momentum_world_kg_m_s),
        boost.scale(ar.expelled_mass_kg),
        2e-15,
    );
    near(
        br.exhaust_energy_world_j - ar.exhaust_energy_world_j,
        boost.dot(ar.exhaust_momentum_world_kg_m_s)
            + 0.5 * ar.expelled_mass_kg * boost.norm_squared(),
        1e-10,
    );
    assert!(ar.gravity_work_uncertainty_j >= 0.0);
    check_channel_accounts(&before, &a);
    check_channel_accounts(&boosted, &b);
}

#[test]
fn closed_elastic_and_zero_length_limits_preserve_original_wall_and_aperture_physics() {
    let before = initial(PURE, Vec3::ZERO);
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1900.0).unwrap(),
    )]);
    let normal = Vec3::new(1.0, 0.0, 0.0);
    for (h, alpha, area) in [(2.0, 1.0, 0.0), (2.0, 0.0, AREA), (0.0, 1.0, AREA)] {
        let path = channel(&before, h, alpha);
        let port = EffusionChannelOutlet {
            area_m2: area,
            ..channel_outlet(&path, normal)
        };
        let result = before
            .advance_effusion_channels(
                environment(),
                &walls,
                &[port],
                DT,
                endpoint(),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(
            result.walls()[&9].specific_energy_j_kg(),
            walls[&9].specific_energy_j_kg()
        );
        assert_eq!(result.step().thermal.wall_energy_change_j, 0.0);
        if area == 0.0 {
            assert_eq!(result.endpoint(), &before);
            assert_eq!(result.step().recoil.expelled_mass_kg, 0.0);
        } else if h == 0.0 {
            let thin = interface(before.reservoir().model());
            let aperture = before
                .advance_effusion(
                    FRAME,
                    &[outlet(&thin, normal)],
                    DT,
                    endpoint(),
                    controls(),
                    || false,
                )
                .unwrap();
            near(
                result.endpoint().reservoir().mass_kg(),
                aperture.endpoint().reservoir().mass_kg(),
                1e-18,
            );
            near(
                result.endpoint().reservoir().state().state().temperature,
                aperture.endpoint().reservoir().state().state().temperature,
                3e-8,
            );
            vector_near(
                result.endpoint().vehicle().velocity_m_s,
                aperture.endpoint().vehicle().velocity_m_s,
                1e-12,
            );
            near(
                result.step().recoil.relative_energy_j,
                aperture.step().relative_energy_j,
                1e-10,
            );
        }
        check_channel_accounts(&before, &result);
    }
}

#[test]
fn wall_mass_frames_gross_withdrawal_and_late_cancellation_refuse_atomically() {
    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let original_wall = walls[&9].specific_energy_j_kg();
    let path = channel(&before, 2.0, 1.0);
    let port = channel_outlet(&path, Vec3::new(1.0, 0.0, 0.0));
    let mut polls = 0;
    before
        .advance_effusion_channels(
            environment(),
            &walls,
            &[port],
            DT,
            endpoint(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut again = 0;
    assert!(matches!(
        before.advance_effusion_channels(
            environment(),
            &walls,
            &[port],
            DT,
            endpoint(),
            controls(),
            || {
                again += 1;
                again == polls
            }
        ),
        Err(EffusionChannelError::Recoil(EffusionRocketError::Cancelled))
    ));
    let oversized = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 2e-5, 900.0).unwrap(),
    )]);
    assert!(matches!(
        before.advance_effusion_channels(
            environment(),
            &oversized,
            &[port],
            DT,
            endpoint(),
            controls(),
            || false
        ),
        Err(EffusionChannelError::Invalid(_))
    ));
    let mismatched = EffusionRocketEnvironment {
        inertial_frame: FRAME + 1,
        ..environment()
    };
    assert!(matches!(
        before.advance_effusion_channels(
            mismatched,
            &walls,
            &[port],
            DT,
            endpoint(),
            controls(),
            || false
        ),
        Err(EffusionChannelError::Invalid(_))
    ));
    let incident_fraction = escaped_moles(0.028, 1.0) * 0.028 / before.reservoir().mass_kg();
    let transmitted_fraction =
        incident_fraction * path.channel().geometry().transmission_probability();
    let mut small = controls();
    small.transport.exchange.maximum_species_fraction =
        0.5 * (incident_fraction + transmitted_fraction);
    assert!(
        before
            .advance_effusion_channels(
                environment(),
                &walls,
                &[port],
                DT,
                endpoint(),
                small,
                || false
            )
            .is_err(),
        "reflected gross withdrawals cannot disappear behind net transmitted mass"
    );
    assert_eq!(before, saved);
    assert_eq!(walls[&9].specific_energy_j_kg(), original_wall);
}
