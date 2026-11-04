use super::*;
use crate::resistive_film::isothermal::AmbientGas;
use crate::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm};

const P0: f64 = 100_000.0;
const RT: f64 = 293.15 * 287.05;
const MU: f64 = 1.8e-5;

fn gap(h: f64, b: &[f64]) -> GapPort {
    GapPort {
        reference_m: h,
        closure: b.to_vec(),
    }
}
fn film() -> IsothermalFilm {
    let geometry = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: 0.002,
                gap: gap(0.001, &[1.0, -1.0]),
            },
            FilmCell {
                area_m2: 0.003,
                gap: gap(0.0012, &[0.7, -0.7]),
            },
        ],
        vec![
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: 0.04,
                length_m: 0.02,
                gap: gap(0.0011, &[0.85, -0.85]),
            },
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.03,
                length_m: 0.01,
                gap: gap(0.001, &[1.0, -1.0]),
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: 0.02,
                length_m: 0.015,
                gap: gap(0.0012, &[0.7, -0.7]),
            },
        ],
        2,
        MU,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 0.01,
            maximum_pressure_pa: 300_000.0,
        },
    )
    .unwrap();
    IsothermalFilm::at_ambient(
        geometry,
        &[0.0; 2],
        AmbientGas {
            pressure_pa: P0,
            temperature_k: 293.15,
            specific_gas_constant_j_kg_k: 287.05,
        },
    )
    .unwrap()
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-14),
        "{a} != {b}"
    );
}
fn observation_eq(a: GasObservation, b: GasObservation) {
    assert_eq!(a.mass_kg, b.mass_kg);
    assert_eq!(a.free_energy_j, b.free_energy_j);
    assert_eq!(a.minimum_pressure_pa, b.minimum_pressure_pa);
    assert_eq!(a.maximum_pressure_pa, b.maximum_pressure_pa);
}

#[test]
fn independent_supplies_recover_poiseuille_mass_and_complete_free_energy_work() {
    for (pressures, supplies) in [
        ([130_000.0, 90_000.0], [250_000.0, 50_000.0]),
        ([90_000.0, 130_000.0], [50_000.0, 250_000.0]),
    ] {
        let original = film();
        let z = [
            original.initial_state()[0] * pressures[0] / P0,
            original.initial_state()[1] * pressures[1] / P0,
        ];
        let driven = PressureDrivenGasFilm::new(
            original.clone(),
            &[None, Some(supplies[0]), Some(supplies[1])],
        )
        .unwrap();
        let (mut gq, mut e, mut flow) = ([0.0; 2], [0.0; 2], [0.0; 2]);
        driven
            .gradient_into(&[0.0; 2], &z, &mut gq, &mut e)
            .unwrap();
        let mut channels = [GasChannelFlow::default(); 3];
        let report = driven
            .flow_into(&[0.0; 2], &z, &e, &mut flow, &mut channels)
            .unwrap();
        let conductance = [
            0.04 * 0.0011_f64.powi(3) / (12.0 * MU * 0.02),
            0.03 * 0.001_f64.powi(3) / (12.0 * MU * 0.01),
            0.02 * 0.0012_f64.powi(3) / (12.0 * MU * 0.015),
        ];
        let ends = [
            (pressures[0], pressures[1]),
            (pressures[0], supplies[0]),
            (pressures[1], supplies[1]),
        ];
        let mut expected_mass = [0.0; 3];
        let mut expected_loss = 0.0;
        let mut expected_supply = 0.0;
        for k in 0..3 {
            let (pi, pj) = ends[k];
            expected_mass[k] = conductance[k] * (pi * pi - pj * pj) / (2.0 * RT);
            let loss = RT * expected_mass[k] * (pi / pj).ln();
            let supply = if k == 0 {
                0.0
            } else {
                -RT * expected_mass[k] * (pj / P0).ln()
            };
            close(channels[k].mass_flow_kg_s, expected_mass[k], 3e-14);
            close(channels[k].dissipated_power_w, loss, 3e-14);
            close(channels[k].reservoir_power_w, supply, 3e-14);
            assert!(channels[k].dissipated_power_w >= 0.0);
            expected_loss += loss;
            expected_supply += supply;
        }
        close(
            flow[0] * original.scales[0] / RT,
            expected_mass[0] + expected_mass[1],
            3e-14,
        );
        close(
            flow[1] * original.scales[1] / RT,
            -expected_mass[0] + expected_mass[2],
            3e-14,
        );
        close(
            report.reservoir_mass_rate_kg_s,
            -expected_mass[1] - expected_mass[2],
            3e-14,
        );
        close(report.dissipated_power_w, expected_loss, 3e-14);
        close(report.reservoir_power_w, expected_supply, 3e-14);
        close(
            report.mass_storage_power_w + expected_loss,
            expected_supply,
            3e-14,
        );
        observation_eq(
            driven.observe(&[0.0; 2], &z).unwrap(),
            original.observe(&[0.0; 2], &z).unwrap(),
        );
    }
}

#[test]
fn arbitrary_joint_effort_retains_signed_supply_and_mass_balance() {
    let driven =
        PressureDrivenGasFilm::new(film(), &[None, Some(190_000.0), Some(65_000.0)]).unwrap();
    let q = [0.0001, -0.00002];
    let z = driven.initial_state();
    for e in [[0.0; 2], [0.1, -0.08], [-0.5, 0.6], [1.0, 1.1]] {
        let mut flow = [0.0; 2];
        let mut channels = [GasChannelFlow::default(); 3];
        let report = driven
            .flow_into(&q, z, &e, &mut flow, &mut channels)
            .unwrap();
        let mass_rate = -(flow[0] * driven.film.scales[0] + flow[1] * driven.film.scales[1]) / RT;
        close(mass_rate, report.reservoir_mass_rate_kg_s, 5e-14);
        assert!(report.dissipated_power_w >= 0.0);
        close(
            report.mass_storage_power_w,
            -e[0] * flow[0] - e[1] * flow[1],
            2e-14,
        );
        close(
            report.mass_storage_power_w + report.dissipated_power_w,
            report.reservoir_power_w,
            5e-14,
        );
        if e == [0.0; 2] {
            assert!(report.reservoir_power_w > 0.0);
            assert_eq!(report.mass_storage_power_w, 0.0);
        }
    }
}

#[test]
fn pressure_supplies_support_a_static_load_with_nonzero_throughflow() {
    let mut base = film();
    // One original cell, two exterior passages; derive its steady p² directly.
    let g = &mut base.geometry;
    let cells = vec![g.cells[0].clone()];
    let mut outlet = g.channels[2].clone();
    outlet.from = 0;
    *g = ResistiveFilm::new(cells, vec![g.channels[1].clone(), outlet], 2, MU, g.limits).unwrap();
    base.initial_gap.truncate(1);
    base.scales.truncate(1);
    let a = 0.03 * 0.001_f64.powi(3) / 0.01;
    let b = 0.02 * 0.0012_f64.powi(3) / 0.015;
    let supply = [180_000.0_f64, 80_000.0_f64];
    let pressure = ((a * supply[0].powi(2) + b * supply[1].powi(2)) / (a + b)).sqrt();
    let driven = PressureDrivenGasFilm::new(base, &[Some(supply[0]), Some(supply[1])]).unwrap();
    let z = [driven.initial_state()[0] * pressure / P0];
    let (mut gq, mut effort, mut flow) = ([0.0; 2], [0.0], [0.0]);
    driven
        .gradient_into(&[0.0; 2], &z, &mut gq, &mut effort)
        .unwrap();
    let mut channels = [GasChannelFlow::default(); 2];
    let report = driven
        .flow_into(&[0.0; 2], &z, &effort, &mut flow, &mut channels)
        .unwrap();
    close(gq[0], 0.002 * (pressure - P0), 2e-14);
    assert_eq!(gq[0], -gq[1]);
    assert!(channels[0].mass_flow_kg_s < 0.0 && channels[1].mass_flow_kg_s > 0.0);
    close(
        channels[0].mass_flow_kg_s,
        -channels[1].mass_flow_kg_s,
        5e-14,
    );
    assert!(report.reservoir_mass_rate_kg_s.abs() < 1e-14);
    assert!(report.mass_storage_power_w.abs() < 1e-9);
    close(report.reservoir_power_w, report.dissipated_power_w, 5e-14);
}

#[test]
fn uniform_ambient_is_bitwise_legacy_and_viscosity_rebinding_preserves_inventory() {
    let original = film();
    let ambient =
        PressureDrivenGasFilm::new(original.clone(), &[None, Some(P0), Some(P0)]).unwrap();
    let q = [0.0001, -0.00003];
    let z = [original.scales[0] * 1.07, original.scales[1] * 0.94];
    let e = [0.03, -0.07];
    let dq = [0.0001, 0.00004];
    let dz = [0.03, -0.04];
    let de = [0.01, -0.03];
    let (mut a, mut b) = ([0.0; 2], [0.0; 2]);
    let loss = original.flow_into(&q, &z, &e, &mut a).unwrap();
    let mut channels = [GasChannelFlow::default(); 3];
    let report = ambient
        .flow_into(&q, &z, &e, &mut b, &mut channels)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(loss, report.dissipated_power_w);
    assert_eq!(report.reservoir_power_w, 0.0);
    original
        .flow_tangent_into(&q, &z, &e, &dq, &dz, &de, &mut a)
        .unwrap();
    ambient
        .flow_tangent_into(
            &q,
            &z,
            &e,
            &dq,
            &dz,
            &de,
            &[None, Some(0.0), Some(0.0)],
            &mut b,
        )
        .unwrap();
    assert_eq!(a, b);

    let driven = ambient
        .with_reservoir_pressures(&[None, Some(160_000.0), Some(70_000.0)])
        .unwrap();
    let current = driven.with_viscosity(3.0 * MU).unwrap();
    let spatial = driven.with_channel_viscosities(&[3.0 * MU; 3]).unwrap();
    for changed in [&current, &spatial] {
        assert_eq!(changed.initial_state(), driven.initial_state());
        assert_eq!(
            changed.reservoir_pressures_pa().collect::<Vec<_>>(),
            driven.reservoir_pressures_pa().collect::<Vec<_>>()
        );
        observation_eq(
            changed.observe(&q, &z).unwrap(),
            driven.observe(&q, &z).unwrap(),
        );
        let ra = driven.flow_into(&q, &z, &e, &mut a, &mut channels).unwrap();
        let rb = changed
            .flow_into(&q, &z, &e, &mut b, &mut channels)
            .unwrap();
        for i in 0..2 {
            close(b[i], a[i] / 3.0, 2e-14);
        }
        close(rb.dissipated_power_w, ra.dissipated_power_w / 3.0, 2e-14);
        close(rb.reservoir_power_w, ra.reservoir_power_w / 3.0, 2e-14);
    }
}

#[test]
fn exact_flow_tangent_includes_reservoir_mobility_and_chemical_effort() {
    for (q, supplies) in [
        ([0.0; 2], [P0, P0]),
        ([0.0001, -0.00003], [210_000.0, 60_000.0]),
    ] {
        let driven =
            PressureDrivenGasFilm::new(film(), &[None, Some(supplies[0]), Some(supplies[1])])
                .unwrap();
        let z = driven.initial_state();
        let e = [0.03, -0.06];
        let dq = [0.0001, 0.00004];
        let dz = [0.02, -0.04];
        let de = [0.03, 0.02];
        let db = [20_000.0, -17_000.0];
        let mut tangent = [0.0; 2];
        driven
            .flow_tangent_into(
                &q,
                z,
                &e,
                &dq,
                &dz,
                &de,
                &[None, Some(db[0]), Some(db[1])],
                &mut tangent,
            )
            .unwrap();
        let eps = 1e-5;
        let mut samples = [[0.0; 2]; 2];
        for (k, sign) in [1.0, -1.0].into_iter().enumerate() {
            let model = driven
                .with_reservoir_pressures(&[
                    None,
                    Some(supplies[0] + sign * eps * db[0]),
                    Some(supplies[1] + sign * eps * db[1]),
                ])
                .unwrap();
            let qa = core::array::from_fn::<_, 2, _>(|i| q[i] + sign * eps * dq[i]);
            let za = core::array::from_fn::<_, 2, _>(|i| z[i] + sign * eps * dz[i]);
            let ea = core::array::from_fn::<_, 2, _>(|i| e[i] + sign * eps * de[i]);
            model
                .flow_into(
                    &qa,
                    &za,
                    &ea,
                    &mut samples[k],
                    &mut [GasChannelFlow::default(); 3],
                )
                .unwrap();
        }
        for i in 0..2 {
            close(
                tangent[i],
                (samples[0][i] - samples[1][i]) / (2.0 * eps),
                2e-8,
            );
        }
    }
}

#[test]
fn sealed_supplies_can_rebind_and_reopen_without_resetting_mass() {
    let mut original = film();
    let passages = Arc::make_mut(&mut original.geometry.channels);
    for passage in passages {
        passage.gap = gap(0.0001, &[1.0, -1.0]);
    }
    let driven =
        PressureDrivenGasFilm::new(original, &[None, Some(200_000.0), Some(60_000.0)]).unwrap();
    let z = [
        driven.initial_state()[0] * 1.1,
        driven.initial_state()[1] * 0.9,
    ];
    let q = [0.0002, 0.0];
    let before = driven.observe(&q, &z).unwrap();
    let changed = driven
        .with_reservoir_pressures(&[None, Some(60_000.0), Some(200_000.0)])
        .unwrap();
    observation_eq(before, changed.observe(&q, &z).unwrap());
    let mut flow = [17.0; 2];
    let mut channels = [GasChannelFlow::default(); 3];
    let report = changed
        .flow_into(&q, &z, &[0.4, -0.3], &mut flow, &mut channels)
        .unwrap();
    assert_eq!(flow, [0.0; 2]);
    assert_eq!(report.dissipated_power_w, 0.0);
    assert_eq!(report.reservoir_power_w, 0.0);
    assert!(channels.iter().all(|c| c.mass_flow_kg_s == 0.0));
    changed
        .flow_into(&[0.0; 2], &z, &[0.4, -0.3], &mut flow, &mut channels)
        .unwrap();
    assert!(flow.iter().any(|x| *x != 0.0));
    assert_eq!(
        changed.observe(&[0.0; 2], &z).unwrap().mass_kg,
        before.mass_kg
    );
}

#[test]
fn boundary_domain_and_output_refusals_are_atomic() {
    let original = film();
    for pressures in [
        vec![],
        vec![None, Some(P0)],
        vec![Some(P0); 3],
        vec![None, None, Some(P0)],
        vec![None, Some(0.0), Some(P0)],
        vec![None, Some(-1.0), Some(P0)],
        vec![None, Some(f64::NAN), Some(P0)],
        vec![None, Some(f64::INFINITY), Some(P0)],
        vec![None, Some(400_001.0), Some(P0)],
    ] {
        assert!(PressureDrivenGasFilm::new(original.clone(), &pressures).is_err());
    }
    let driven =
        PressureDrivenGasFilm::new(original, &[None, Some(180_000.0), Some(70_000.0)]).unwrap();
    let z = driven.initial_state();
    let mut out = [17.0; 2];
    let sentinel = GasChannelFlow {
        mass_flow_kg_s: 19.0,
        dissipated_power_w: 23.0,
        reservoir_power_w: 29.0,
    };
    let mut channels = [sentinel; 3];
    for (q, mass, effort) in [
        ([0.002, 0.0], z.to_vec(), [0.0; 2]),
        ([0.0; 2], vec![0.0, z[1]], [0.0; 2]),
        ([0.0; 2], z.to_vec(), [f64::NAN, 0.0]),
        ([0.0; 2], z.to_vec(), [f64::MAX; 2]),
    ] {
        assert!(
            driven
                .flow_into(&q, &mass, &effort, &mut out, &mut channels)
                .is_err()
        );
        assert_eq!(out, [17.0; 2]);
        assert!(channels.iter().all(|c| c.mass_flow_kg_s == 19.0
            && c.dissipated_power_w == 23.0
            && c.reservoir_power_w == 29.0));
    }
    assert!(
        driven
            .flow_into(&[0.0; 2], z, &[0.0; 2], &mut out, &mut channels[..2])
            .is_err()
    );
    assert_eq!(out, [17.0; 2]);
    for directions in [
        vec![],
        vec![Some(0.0); 3],
        vec![None, None, Some(0.0)],
        vec![None, Some(f64::NAN), Some(0.0)],
    ] {
        assert!(
            driven
                .flow_tangent_into(
                    &[0.0; 2],
                    z,
                    &[0.0; 2],
                    &[0.0; 2],
                    &[0.0; 2],
                    &[0.0; 2],
                    &directions,
                    &mut out
                )
                .is_err()
        );
        assert_eq!(out, [17.0; 2]);
    }
}
