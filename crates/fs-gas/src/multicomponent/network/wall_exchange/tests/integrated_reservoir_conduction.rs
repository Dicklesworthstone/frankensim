//! Original source/film heat from a prescribed reservoir, without wall storage.

use super::conduction::source;
use super::*;
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, WallConductionEvaluation,
};
use fs_material::thermal::conduction::ReservoirLayerGeometry;
use fs_thermal::enthalpy::DepositionError;

const RESERVOIR: f64 = 5.0;
const GAS_CAPACITY: f64 = 0.5;

fn reservoir(
    gas: &MixtureGas,
    geometry: ReservoirLayerGeometry,
    film: f64,
    span: [f64; 2],
    integrated: bool,
) -> (MixtureNetwork1d, NetworkWallOperator<'static>) {
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(0.5),
        0.4,
        vec![gas.conserved(&[1.0], 0.0, 2.0).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    // Every physical gas volume has Cgas=0.5, despite different rho and volume.
    let plenum = Plenum::new(0.1, gas.conserved(&[2.0], 0.0, 4.0).unwrap()).unwrap();
    let network = MixtureNetwork1d::new(vec![branch], vec![plenum], vec![], || false).unwrap();
    let (transport, _) = source([0.0, 1.0], [1.0, 0.2], span);
    let layer = MaterialWallConduction::from_geometry(transport, geometry).unwrap();
    let layer = if integrated {
        layer.with_integrated_profile().unwrap()
    } else {
        layer
    };
    // No MaterialWallCaloric binding: prescribed temperature is an external
    // reservoir, not a finite solid amount inferred from transport density/Cv.
    let wall = uniform_operator(
        &network,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Reservoir {
                temperature_k: RESERVOIR,
                conductance_w_per_k: film,
            },
        },
    )
    .with_material_conduction(
        vec![vec![Some(layer.clone()); 2]],
        vec![Some(layer)],
        || false,
    )
    .unwrap();
    (network, wall)
}

fn original_heat(temperature: f64, multiplier: f64, film: f64, dt: f64) -> f64 {
    // Eliminate the original gas and film equations:
    // Tg=Told+Q/Cgas; Ts=Told+B*Q; Q=dt*C*(Tw^2-Ts^2)/2.
    // The positive quadratic root uses only original physical input data.
    let b_face = 1.0 / GAS_CAPACITY + 1.0 / (dt * film);
    let a = 0.5 * dt * multiplier * b_face * b_face;
    let b = 1.0 + dt * multiplier * temperature * b_face;
    let c = 0.5 * dt * multiplier * (temperature * temperature - RESERVOIR * RESERVOIR);
    -2.0 * c / (b + (b * b - 4.0 * a * c).sqrt())
}

#[test]
fn prescribed_profiles_close_original_films_and_external_heat_over_repeated_steps() {
    let g = gas(1, false, 0.0);
    let geometries = [
        (ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(), 0.1),
        (
            ReservoirLayerGeometry::cylindrical_face(0.2, 0.1, 0.01).unwrap(),
            std::f64::consts::TAU * 0.01 / 2.0_f64.ln(),
        ),
    ];
    for (geometry, multiplier) in geometries {
        // Manufactured first endpoint: Tw=5, Ts=4, Tg=3, Q=0.5 J.
        // C is independently A/L or 2*pi*length/ln(r_outer/r_inner).
        let power = 0.5 * multiplier * (25.0 - 16.0);
        let film = power / (4.0 - 3.0);
        let dt = 0.5 / power;
        let (mut n, mut wall) = reservoir(&g, geometry, film, [3.5, 10.0], true);
        let layer = wall.branch_conduction(0, 0).unwrap().clone();
        assert!(layer.has_integrated_profile());
        assert!(layer.source().unwrap().at_temperature(3.0).is_err());
        assert!(layer.source().unwrap().at_temperature(4.0).is_ok());
        let mut temperature = 2.0;
        for step in 0..3 {
            let heat = original_heat(temperature, multiplier, film, dt);
            temperature += heat / GAS_CAPACITY;
            let face = temperature + heat / (dt * film);
            if step == 0 {
                close(heat, 0.5, 3e-13);
                close(temperature, 3.0, 3e-13);
                close(face, 4.0, 3e-13);
            }
            let before = fields(&n);
            let initial_energy = n.totals().absolute_energy(&g).unwrap();
            let report = wall
                .advance_for(&mut n, Time::new(dt), controls(), || false)
                .unwrap();
            let after = fields(&n);
            for i in 0..3 {
                let (old, new, volume, contact) = if i < 2 {
                    (
                        &before.0[0][i],
                        &after.0[0][i],
                        0.2,
                        wall.branch_contacts(0).unwrap()[i],
                    )
                } else {
                    (&before.1[0], &after.1[0], 0.1, wall.plenum_contacts()[0])
                };
                let tg = g.primitive(new.as_cell(), || false).unwrap().temperature;
                let q = volume * (new.energy - old.energy);
                let actual_face = tg + q / (dt * film);
                close(tg, temperature, 3e-10);
                close(q, heat, 3e-10);
                close(actual_face, face, 3e-10);
                close(
                    q,
                    dt * 0.5 * multiplier * (25.0 - actual_face * actual_face),
                    3e-10,
                );
                close(q, dt * film * (face - tg), 3e-10);
                assert_eq!(new.partial_densities, old.partial_densities);
                assert_eq!(new.momentum, old.momentum);
                let WallThermal::Reservoir {
                    temperature_k,
                    conductance_w_per_k,
                } = contact.thermal
                else {
                    panic!("conduction must preserve prescribed reservoir ownership");
                };
                assert_eq!(temperature_k, RESERVOIR);
                close(
                    conductance_w_per_k,
                    heat / (dt * (RESERVOIR - temperature)),
                    3e-10,
                );
            }
            assert!(
                wall.branch_heat_storage(0)
                    .unwrap()
                    .iter()
                    .all(Option::is_none)
            );
            assert!(wall.plenum_heat_storage().iter().all(Option::is_none));
            assert!(wall.shared_wall(0).is_none());
            assert_eq!(wall.branch_conduction(0, 0), Some(&layer));
            assert_eq!(wall.plenum_conduction(0), Some(&layer));
            close(report.heat_to_gas_j, 3.0 * heat, 3e-10);
            close(report.external_heat_to_gas_j, 3.0 * heat, 3e-10);
            assert_eq!(report.wall_energy_change_j, 0.0);
            close(
                n.totals().absolute_energy(&g).unwrap() - initial_energy,
                report.external_heat_to_gas_j,
                3e-10,
            );
            close(report.energy_residual_j, 0.0, 3e-10);
            close(report.maximum_thermal_residual_j, 0.0, 3e-10);
            assert_eq!(n.time(), Time::new(0.0));
        }
    }
}

#[test]
fn prescribed_profile_refusal_and_late_cancellation_preserve_accepted_state() {
    let g = gas(1, false, 0.0);
    let geometry = ReservoirLayerGeometry::planar(0.1, 1.0).unwrap();
    let dt = Time::new(10.0 / 9.0);
    let (mut denied, mut denied_wall) = reservoir(&g, geometry, 0.45, [4.1, 10.0], true);
    let before = fields(&denied);
    let contacts = denied_wall.branch_contacts(0).unwrap().to_vec();
    let plenums = denied_wall.plenum_contacts().to_vec();
    // Prescribed Tw=5 is supported, but the required material face Ts=4 is not.
    assert!(matches!(
        denied_wall.advance_for(&mut denied, dt, controls(), || false),
        Err(WallExchangeError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(fields(&denied), before);
    assert_eq!(denied_wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(denied_wall.plenum_contacts(), plenums);
    assert!(denied_wall.last_report().is_none());
    assert_eq!(denied.time(), Time::new(0.0));

    let (mut reference, mut reference_wall) = reservoir(&g, geometry, 0.45, [3.5, 10.0], true);
    let first = reference_wall
        .advance_for(&mut reference, dt, controls(), || false)
        .unwrap();
    let (mut retry, mut retry_wall) = reservoir(&g, geometry, 0.45, [3.5, 10.0], true);
    assert_eq!(
        retry_wall
            .advance_for(&mut retry, dt, controls(), || false)
            .unwrap(),
        first
    );
    let accepted = fields(&retry);
    let contacts = retry_wall.branch_contacts(0).unwrap().to_vec();
    let plenums = retry_wall.plenum_contacts().to_vec();
    let mut polls = 0;
    let next = reference_wall
        .advance_for(&mut reference, dt, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        retry_wall
            .advance_for(&mut retry, dt, controls(), || {
                seen += 1;
                seen >= polls
            })
            .unwrap_err(),
        WallExchangeError::Cancelled
    );
    assert_eq!(fields(&retry), accepted);
    assert_eq!(retry_wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(retry_wall.plenum_contacts(), plenums);
    assert_eq!(retry_wall.last_report(), Some(first));
    assert_eq!(retry.time(), Time::new(0.0));
    assert_eq!(
        retry_wall
            .advance_for(&mut retry, dt, controls(), || false)
            .unwrap(),
        next
    );
    assert_eq!(fields(&retry), fields(&reference));
    assert_eq!(
        retry_wall.branch_contacts(0),
        reference_wall.branch_contacts(0)
    );
    assert_eq!(
        retry_wall.plenum_contacts(),
        reference_wall.plenum_contacts()
    );
}

#[test]
fn original_pointwise_reservoir_path_keeps_its_sampling_and_external_balance() {
    let g = gas(1, false, 0.0);
    let (mut n, mut wall) = reservoir(
        &g,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
        0.45,
        [3.5, 10.0],
        false,
    );
    let layer = wall.branch_conduction(0, 0).unwrap();
    assert!(!layer.has_integrated_profile());
    assert_eq!(layer.evaluation(), WallConductionEvaluation::InitialState);
    // Original point sampling uses k(5)=5, solid C*k=0.5 W/K and film=0.45.
    let conductance = 1.0 / (1.0 / 0.45 + 1.0 / 0.5);
    let dt = 10.0 / 9.0;
    let mut temperature = 2.0;
    for step in 0..3 {
        let heat =
            dt * conductance * (RESERVOIR - temperature) / (1.0 + dt * conductance / GAS_CAPACITY);
        temperature += heat / GAS_CAPACITY;
        if step == 0 {
            assert!(
                (heat - 0.5).abs() > 1e-3,
                "opting in changes the physical law"
            );
        }
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap();
        for i in 0..2 {
            close(
                g.primitive(n.branches()[0].cell(i).unwrap(), || false)
                    .unwrap()
                    .temperature,
                temperature,
                3e-10,
            );
        }
        close(
            g.primitive(n.plenums()[0].state(), || false)
                .unwrap()
                .temperature,
            temperature,
            3e-10,
        );
        close(report.heat_to_gas_j, 3.0 * heat, 3e-10);
        close(report.external_heat_to_gas_j, 3.0 * heat, 3e-10);
        assert_eq!(report.wall_energy_change_j, 0.0);
        close(report.energy_residual_j, 0.0, 3e-10);
    }
}

#[test]
fn prescribed_wall_above_gas_temperature_support_keeps_a_supported_physical_endpoint() {
    use crate::dynamics::CaloricDomain;
    use fs_qty::Temperature;

    let original = gas(1, false, 0.0);
    let g = MixtureGas::new(
        original.source().clone(),
        CaloricDomain::new(Temperature::new(1.0), Temperature::new(10.0)),
        || false,
    )
    .unwrap();
    // Original Fourier/film/caloric equations give Tw=12, Ts=10, Tg=8,
    // P=.1*(12^2-10^2)/2=2.2 W, and Q=.2*P=.44 J for each Cgas=.5.
    // The wall is outside the GAS chart; no gas query at Tw is permissible.
    let initial_temperature = 7.12;
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.4,
        vec![g.conserved(&[1.0], 0.0, initial_temperature).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let plenum = Plenum::new(
        0.1,
        g.conserved(&[2.0], 0.0, 2.0 * initial_temperature).unwrap(),
    )
    .unwrap();
    let mut n = MixtureNetwork1d::new(vec![branch], vec![plenum], vec![], || false).unwrap();
    let (transport, _) = source([0.0, 1.0], [1.0, 0.2], [1.0, 20.0]);
    let layer = MaterialWallConduction::from_geometry(
        transport,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
    )
    .unwrap()
    .with_integrated_profile()
    .unwrap();
    let mut wall = uniform_operator(
        &n,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Reservoir {
                temperature_k: 12.0,
                conductance_w_per_k: 1.1,
            },
        },
    )
    .with_material_conduction(
        vec![vec![Some(layer.clone()); 2]],
        vec![Some(layer)],
        || false,
    )
    .unwrap();
    let before = n.totals().absolute_energy(&g).unwrap();
    let report = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || false)
        .unwrap();
    for cell in [
        n.branches()[0].cell(0).unwrap(),
        n.branches()[0].cell(1).unwrap(),
        n.plenums()[0].state(),
    ] {
        close(g.primitive(cell, || false).unwrap().temperature, 8.0, 3e-10);
    }
    close(report.heat_to_gas_j, 1.32, 3e-10);
    close(report.external_heat_to_gas_j, 1.32, 3e-10);
    assert_eq!(report.wall_energy_change_j, 0.0);
    close(
        n.totals().absolute_energy(&g).unwrap() - before,
        1.32,
        3e-10,
    );
    close(report.energy_residual_j, 0.0, 3e-10);
    close(report.maximum_thermal_residual_j, 0.0, 3e-10);
    for contact in wall
        .branch_contacts(0)
        .unwrap()
        .iter()
        .chain(wall.plenum_contacts())
    {
        let WallThermal::Reservoir {
            temperature_k,
            conductance_w_per_k,
        } = contact.thermal
        else {
            panic!("prescribed temperature must retain its external energy owner");
        };
        assert_eq!(temperature_k, 12.0);
        close(conductance_w_per_k, 0.55, 3e-10);
    }
}
