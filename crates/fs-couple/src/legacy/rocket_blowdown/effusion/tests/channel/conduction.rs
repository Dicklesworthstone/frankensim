//! G1/G3/G4: original material conduction changes later molecular recoil.
use super::*;
use fs_material::thermal::transport::volumetric::timeline::FixedVolumeConductionLink;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::effusion::channel::thermal::conduction::ChannelWallConduction;
use fs_thermal::enthalpy::{CellId, ThermalLink};

fn connection<'l, 'p>(law: &'l FixedVolumeThermalLaw<'p>) -> ChannelWallConduction<'l, 'p> {
    ChannelWallConduction::Cellwise(
        FixedVolumeConductionLink::try_aligned(
            &ThermalLink {
                id: 17,
                a: CellId(9),
                b: CellId(12),
                area_m2: 1e-4,
                // The prepared source, distances and contact law supply G.
                conductance_w_k: 0.0,
            },
            [law, law],
            [0.005; 2],
            0.0,
        )
        .unwrap(),
    )
}

fn heated_temperature(t: f64, heat: f64) -> f64 {
    let capacity = 1e6 + 1000.0 * t;
    t + (2.0 * heat / 1e-6)
        / ((capacity * capacity + 2.0 * 1000.0 * heat / 1e-6).sqrt() + capacity)
}

#[test]
fn internal_conduction_heats_unvented_walls_without_external_force_or_work() {
    let before = initial(PURE, Vec3::new(23.0, -7.0, 2.0));
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 400.0).unwrap()),
        (12, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1000.0).unwrap()),
    ]);
    let conduction = [connection(&law)];
    let result = before
        .advance_effusion_channels_with_conduction(
            environment(), &walls, &[], &conduction, DT, endpoint(), controls(), || false,
        )
        .unwrap();
    let heat = 1e-4 / (0.005 + 0.005) * (1000.0 - 400.0) * DT;
    near(result.walls()[&9].temperature_k(), heated_temperature(400.0, heat), 2e-10);
    near(result.walls()[&12].temperature_k(), heated_temperature(1000.0, -heat), 2e-10);
    near(result.step().thermal.wall_energy_change_j, 0.0, 1e-9);
    assert!(result.step().thermal.wall_entropy_change_j_k > 0.0);
    assert_eq!(result.step().thermal.walls.len(), 2);
    assert_eq!(result.step().recoil.relative_momentum_kg_m_s, Vec3::ZERO);
    assert_eq!(result.step().recoil.expelled_mass_kg, 0.0);
    assert_eq!(result.step().recoil.exhaust.external_work_j, 0.0);
    vector_near(result.endpoint().vehicle().velocity_m_s, before.vehicle().velocity_m_s, 1e-13);
    near(result.endpoint().reservoir().mass_kg(), before.reservoir().mass_kg(), 1e-18);
    near(result.endpoint().reservoir().state().state().temperature, T0, 1e-12);
    check_channel_accounts(&before, &result);
    // Publication is a candidate: the caller's heat inventories are untouched.
    near(walls[&9].temperature_k(), 400.0, 1e-12);
    near(walls[&12].temperature_k(), 1000.0, 1e-12);
}

#[test]
fn conduction_and_molecular_heat_share_one_wall_update_and_change_later_escape() {
    let before = initial(PURE, Vec3::new(17.0, -3.0, 0.0));
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 500.0).unwrap()),
        (12, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1000.0).unwrap()),
    ]);
    let interface = channel(&before, 2.0, 0.35);
    let outlets = [channel_outlet(&interface, Vec3::new(1.0, 0.0, 0.0))];
    let conduction = [connection(&law)];
    let unconnected = before
        .advance_effusion_channels(environment(), &walls, &outlets, DT, endpoint(), controls(), || false)
        .unwrap();
    let connected = before
        .advance_effusion_channels_with_conduction(
            environment(), &walls, &outlets, &conduction, DT, endpoint(), controls(), || false,
        )
        .unwrap();
    // Both molecular stages sample the SAME original gas and wall temperatures.
    assert_eq!(connected.step().momentum, unconnected.step().momentum);
    assert_eq!(connected.step().thermal.gas, unconnected.step().thermal.gas);
    let molecular_heat = connected.step().thermal.gas.outlets[0].wall.heat_into_wall_j;
    let conductive_heat = 1e-4 / 0.01 * (1000.0 - 500.0) * DT;
    near(
        connected.walls()[&9].temperature_k(),
        heated_temperature(500.0, molecular_heat + conductive_heat),
        2e-10,
    );
    near(connected.walls()[&12].temperature_k(), heated_temperature(1000.0, -conductive_heat), 2e-10);
    check_channel_accounts(&before, &connected);
    let next = connected.endpoint()
        .advance_effusion_channels_with_conduction(
            environment(), connected.walls(), &outlets, &conduction, DT, endpoint(), controls(), || false,
        )
        .unwrap();
    let unconnected_next = unconnected.endpoint()
        .advance_effusion_channels(environment(), unconnected.walls(), &outlets, DT, endpoint(), controls(), || false)
        .unwrap();
    assert!(next.step().recoil.relative_momentum_kg_m_s.x > unconnected_next.step().recoil.relative_momentum_kg_m_s.x);
    assert!(next.step().recoil.relative_energy_j > unconnected_next.step().recoil.relative_energy_j);
    assert!(next.endpoint().vehicle().velocity_m_s.x < unconnected_next.endpoint().vehicle().velocity_m_s.x);
    check_channel_accounts(connected.endpoint(), &next);
}

#[test]
fn conduction_source_mismatch_or_late_cancellation_publishes_no_vessel_or_wall() {
    let before = initial(PURE, Vec3::ZERO);
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let changed_basis = FixedVolumeThermalLaw::try_new(&material, 2000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 400.0).unwrap()),
        (12, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1000.0).unwrap()),
    ]);
    let saved = before.clone();
    let stored = walls.iter().map(|(&id, wall)| (id, wall.specific_energy_j_kg())).collect::<Vec<_>>();
    assert!(before.advance_effusion_channels_with_conduction(
        environment(), &walls, &[], &[connection(&changed_basis)], DT, endpoint(), controls(), || false,
    ).is_err());
    let conduction = [connection(&law)];
    let mut checks = 0;
    let first = before.advance_effusion_channels_with_conduction(
        environment(), &walls, &[], &conduction, DT, endpoint(), controls(), || { checks += 1; false },
    ).unwrap();
    let mut replay_checks = 0;
    assert!(before.advance_effusion_channels_with_conduction(
        environment(), &walls, &[], &conduction, DT, endpoint(), controls(), || {
            replay_checks += 1;
            replay_checks == checks
        },
    ).is_err());
    assert_eq!(before, saved);
    assert_eq!(walls.iter().map(|(&id, wall)| (id, wall.specific_energy_j_kg())).collect::<Vec<_>>(), stored);
    let replay = before.advance_effusion_channels_with_conduction(
        environment(), &walls, &[], &conduction, DT, endpoint(), controls(), || false,
    ).unwrap();
    assert_eq!(first.endpoint(), replay.endpoint());
    assert_eq!(first.step(), replay.step());
}
