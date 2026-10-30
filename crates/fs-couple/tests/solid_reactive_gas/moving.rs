//! Original native owners, independently checked heat, impulse and pressure work.
use super::*;
use fs_couple::cavity_phs::piston::PistonSolve;
use fs_couple::cavity_phs::piston::material::{MaterialGasPiston, solid::*};
use fs_ga::Vec3;
use fs_matdb::{ClaimSet, MaterialCard, PropertyValue};
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_ports::{ModuleWindowAudit, PortedModule};
use std::collections::BTreeMap;

const DT: f64 = 0.0005;
const AREA: f64 = 0.01;
const FORCE: f64 = -25.0;

fn piston(gas: IsochoricGasParcel, density: f64, fixed: bool) -> MaterialGasPiston {
    let (source, _) = wall(1.0, false);
    let mut claims = ClaimSet::new();
    for (_, claim) in source.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Scalar {
                value: density,
                dims: fs_qty::Density::DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    let source = MaterialCard::assemble(source.id().clone(), claims, vec![]).unwrap();
    let body = resolve_body(
        &[RegionBinding {
            key: "piston",
            material: &source,
            point: &point(300.0),
            geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.01)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap();
    let piston = MaterialGasPiston::new(
        ClosedGasParcel::from_isochoric(gas),
        body,
        AREA,
        if fixed { 0.0 } else { 0.02 },
    )
    .unwrap();
    if fixed {
        piston.with_fixed_support().unwrap()
    } else {
        piston
    }
}

fn inputs<'a>(
    pistons: &'a BTreeMap<u64, MaterialGasPiston>,
    patches: &'a [SolidGasContact],
) -> [SolidPistonCompartment<'a>; 2] {
    [20, 10].map(|id| {
        let piston = &pistons[&id];
        let v = piston.gas().volume_m3();
        SolidPistonCompartment {
            id,
            piston,
            external_heat_j: 0.0,
            contacts: if id == 20 {
                &patches[..1]
            } else {
                &patches[1..]
            },
            reaction: None,
            exterior_pressure_pa: 100000.0,
            external_force_n: if id == 20 { FORCE } else { 0.0 },
            motion: PistonSolve {
                volume_bracket_m3: if id == 20 {
                    [v * 0.99, v * 1.01]
                } else {
                    [v, v]
                },
                gas: solve().gas,
                momentum_tolerance_n_s: 1e-10,
                energy_tolerance_j: 1e-7,
                maximum_iterations: 128,
            },
        }
    })
}

// Eliminate the two finite-gas BE equations, then solve the two wall equations.
fn heat_reference(t: [f64; 4], mass: [f64; 2], b: [f64; 2], k: f64) -> [f64; 4] {
    let capacity = mass.map(|m| m * CV);
    let g = [0.1, 0.2]
        .into_iter()
        .zip(b)
        .map(|(area, b)| DT * area / (0.001 / (0.03 * (1.0 + 3.0 * b)) + 0.02 / k))
        .collect::<Vec<_>>();
    let a = [0, 1].map(|i| g[i] * capacity[i] / (capacity[i] + g[i]));
    let h = DT * 0.5;
    let d = a.map(|a| 10.0 + h + a);
    let rhs = [10.0 * t[0] + a[0] * t[2], 10.0 * t[1] + a[1] * t[3]];
    let det = d[0] * d[1] - h * h;
    let wall = [
        (rhs[0] * d[1] + h * rhs[1]) / det,
        (rhs[1] * d[0] + h * rhs[0]) / det,
    ];
    [
        wall[0],
        wall[1],
        (capacity[0] * t[2] + g[0] * wall[0]) / (capacity[0] + g[0]),
        (capacity[1] * t[3] + g[1] * wall[1]) / (capacity[1] + g[1]),
    ]
}

// Ideal-gas constant-Cv midpoint, solved independently by plain bisection.
// Cv(Tm-T0) + R Tm (V1-V0)/(2 Vm) = 0; T1=2 Tm-T0.
fn motion_reference(
    t: f64,
    gas_mass: f64,
    volume: f64,
    velocity: f64,
    body_mass: f64,
) -> (f64, f64, f64) {
    let at = |v1: f64| {
        let dv = v1 - volume;
        let vm = (v1 + volume) / 2.0;
        let tm = t / (1.0 + RS * dv / (2.0 * CV * vm));
        let next_velocity = 2.0 * dv / (AREA * DT) - velocity;
        let residual = body_mass * (next_velocity - velocity)
            - DT * (AREA * (gas_mass * RS * tm / vm - 100000.0) + FORCE);
        (residual, 2.0 * tm - t, next_velocity)
    };
    let (mut lo, mut hi) = (volume * 0.99, volume * 1.01);
    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        if at(mid).0 < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let v = (lo + hi) / 2.0;
    let (_, t, u) = at(v);
    (t, v, u)
}

fn stored(thermal: &DepositionSimulation<'_>, pistons: &BTreeMap<u64, MaterialGasPiston>) -> f64 {
    thermal.energy_j()
        + pistons
            .values()
            .map(|p| {
                p.gas().mass_kg() * p.gas().state().specific_internal_energy_j_kg()
                    + 0.5 * p.piston().mass_kg() * p.piston().velocity_m_s().powi(2)
            })
            .sum::<f64>()
}

#[test]
fn wall_and_inertia_data_drive_native_motion_and_next_flow_with_one_continuing_clock() {
    with_cx(|cx, _| {
        let law = reaction(1.0, 0.0, 4.0);
        let patches = contacts();
        let mut observations = Vec::new();
        for (k, density) in [(0.1, 1000.0), (10.0, 1000.0), (10.0, 8000.0)] {
            let (card, caloric) = wall(k, false);
            let mut thermal = initial(caloric.curve());
            let context = point(300.0);
            let model = solid(&thermal, &card, caloric.curve(), &context, cx);
            let mut pistons = BTreeMap::from([
                (
                    20,
                    piston(gas_at(&law, 0.8, 400.0, 300000.0), density, false),
                ),
                (
                    10,
                    piston(gas_at(&law, 0.1, 300.0, 100000.0), density, true),
                ),
            ]);
            let bodies = pistons.clone();
            let initial_energy = stored(&thermal, &pistons);
            let (mut external_work, mut last_flow) = (0.0, 0.0);
            let mut t = [300.0, 300.0, 400.0, 300.0];
            let mut m = [pistons[&20].gas().mass_kg(), pistons[&10].gas().mass_kg()];
            let mut b = [0.8, 0.1];
            let (mut volume, mut velocity) = (VOLUME, 0.02);
            for tick in 1..=3 {
                let before = pistons.clone();
                let candidate = model
                    .step_with_material_pistons(
                        &thermal,
                        &inputs(&pistons, &patches),
                        &[link(1e-5)],
                        tick as f64 * DT,
                        None,
                        controls(),
                        &[0.0; 4],
                        cx,
                    )
                    .unwrap();
                let dm = DT * mass_flow(m[0] * RS * t[2] / volume, t[2], 1e-5);
                let after_flow = [
                    t[0],
                    t[1],
                    t[2] * (m[0] - 1.4 * dm) / (m[0] - dm),
                    (m[1] * t[3] + 1.4 * dm * t[2]) / (m[1] + dm),
                ];
                b[1] = (m[1] * b[1] + dm * b[0]) / (m[1] + dm);
                m = [m[0] - dm, m[1] + dm];
                t = heat_reference(after_flow, m, b, k);
                (t[2], volume, velocity) =
                    motion_reference(t[2], m[0], volume, velocity, bodies[&20].piston().mass_kg());
                near(
                    candidate.piston(20).unwrap().gas().volume_m3(),
                    volume,
                    1e-13,
                );
                near(
                    candidate.piston(20).unwrap().piston().velocity_m_s(),
                    velocity,
                    3e-7,
                );
                near(
                    candidate.thermal().temperature_k(CellId(30)).unwrap(),
                    t[0],
                    3e-7,
                );
                for (index, id) in [20, 10].into_iter().enumerate() {
                    let next = candidate.piston(id).unwrap();
                    near(next.gas().mass_kg(), m[index], 1e-12);
                    near(next.gas().state().state().temperature, t[index + 2], 3e-7);
                    near(
                        next.gas()
                            .model()
                            .thermodynamics()
                            .mole_composition()
                            .fractions()[1],
                        b[index],
                        1e-11,
                    );
                    assert_eq!(next.body(), bodies[&id].body());
                    assert_eq!(next.piston().mass_kg(), next.body().inertia().mass());
                    near(next.piston().last_window_duration_s().unwrap(), DT, 1e-18);
                    assert!(matches!(
                        next.piston().audit_last_window().unwrap(),
                        ModuleWindowAudit::Unavailable { .. }
                    ));
                }
                near(
                    candidate.report().heat.flow.links[0].transferred_mass_kg,
                    dm,
                    1e-12,
                );
                last_flow = candidate.report().heat.flow.links[0].flow.mass_flow_kg_s;
                for account in &candidate.report().pistons {
                    assert_eq!(
                        account.motion.velocity_m_s[0],
                        before[&account.id].piston().velocity_m_s()
                    );
                    external_work += account.motion.external_force_work_j
                        - account.motion.exterior_pressure_work_j;
                }
                assert!(
                    candidate.report().total_energy_residual_j.abs()
                        <= controls().total_energy_tolerance_j
                );
                near(
                    stored(candidate.thermal(), candidate.pistons())
                        - initial_energy
                        - external_work,
                    0.0,
                    2e-6,
                );
                assert_eq!(pistons, before);
                let (next_heat, next_pistons, _, _) = candidate.into_parts();
                thermal = next_heat;
                pistons = next_pistons;
            }
            observations.push((pistons[&20].piston().velocity_m_s(), last_flow, t[0]));
            // Closing the nozzle exposes the existing closed two-port audit.
            // Its fourth tick proves preceding flow/heat/motion never reset or
            // advanced the original mechanical port clock twice.
            let candidate = model
                .step_with_material_pistons(
                    &thermal,
                    &inputs(&pistons, &patches),
                    &[],
                    4.0 * DT,
                    None,
                    controls(),
                    &[0.0; 4],
                    cx,
                )
                .unwrap();
            let ModuleWindowAudit::Available(audit) = candidate
                .piston(20)
                .unwrap()
                .piston()
                .audit_last_window()
                .unwrap()
            else {
                panic!("closed chamber account missing")
            };
            assert_eq!(audit.spec().interval().end().tick(), 4);
            near(audit.residuals().first_law().value(), 0.0, 1e-6);
        }
        assert!(observations[1].2 > observations[0].2);
        assert!((observations[1].0 - observations[0].0).abs() > 1e-5);
        assert!((observations[1].1 - observations[0].1).abs() > 1e-9);
        assert!(observations[1].0 > observations[2].0);
        assert!((observations[1].1 - observations[2].1).abs() > 1e-7);
    });
}

#[test]
fn reacting_full_heat_network_retains_formation_energy_and_late_refusal_preserves_all_owners() {
    with_cx(|cx, gate| {
        let law = super::super::coupled::implicit(reaction(1.0, 20000.0, 4.0));
        let (card, caloric) = wall(1.0, false);
        let thermal = initial(caloric.curve());
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let patches = contacts();
        let pistons = BTreeMap::from([
            (
                20,
                piston(gas_at(&law, 0.8, 400.0, 300000.0), 8000.0, false),
            ),
            (10, piston(gas_at(&law, 0.1, 300.0, 100000.0), 8000.0, true)),
        ]);
        let before = pistons.clone();
        let saved_heat = thermal.cells().clone();
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = [30, 10].map(|id| MaterialCell {
            cell: CellId(id),
            material_id: MaterialId(0),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
        });
        let bulk_links = [ConductionLink {
            link: 7,
            distance_a_m: 0.3,
            distance_b_m: 0.7,
        }];
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &bulk_links,
        };
        let electrical =
            fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
        let network = fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let prepared = network.prepare(cx, &thermal).unwrap();
        let mut drives = inputs(&pistons, &patches);
        drives[1].reaction = Some((&law, chemical()));
        let run = |drives: &[SolidPistonCompartment<'_>]| {
            model.step_with_material_pistons_heat_network(
                &thermal,
                drives,
                &[link(1e-5)],
                DT,
                &prepared,
                controls(),
                &[0.0; 4],
                cx,
            )
        };
        let candidate = run(&drives).unwrap();
        let receiver = candidate.piston(10).unwrap().gas();
        let chemistry = candidate.report().heat.thermal.gases[0]
            .endpoint
            .as_ref()
            .unwrap();
        let dm = DT * mass_flow(300000.0, 400.0, 1e-5);
        let product = (pistons[&10].gas().mass_kg() * 0.1 + dm * 0.8) / MOLAR_MASS;
        let amount = receiver.mass_kg() / MOLAR_MASS;
        let t = receiver.state().state().temperature;
        let q = DT * R * t * t / (100000.0 * 300.0);
        let reverse = (-20000.0 / (R * t)).exp();
        let extent = q * (amount - (1.0 + reverse) * product) / (1.0 + q * (1.0 + reverse));
        near(chemistry.reaction.extent_mol, extent, 1e-9);
        near(chemistry.reaction.species[1].initial_mol, product, 1e-10);
        near(
            receiver
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions()[1],
            (product + extent) / amount,
            1e-10,
        );
        let external = candidate
            .report()
            .pistons
            .iter()
            .map(|p| p.motion.external_force_work_j - p.motion.exterior_pressure_work_j)
            .sum::<f64>();
        near(
            stored(candidate.thermal(), candidate.pistons())
                - stored(&thermal, &pistons)
                - external,
            0.0,
            1e-6,
        );
        assert!(candidate.report().heat.thermal.transport.is_some());
        assert!(candidate.piston(20).unwrap().gas().volume_m3() > VOLUME);
        assert_eq!(candidate.piston(10).unwrap().gas().volume_m3(), VOLUME);
        // Canonical id 10 completes before this id-20 mechanical refusal.
        // The complete heat/reaction candidate has already been solved as well.
        let mut bad = drives;
        bad[0].motion.volume_bracket_m3 = [VOLUME, VOLUME];
        assert!(matches!(
            run(&bad),
            Err(MovingSolidGasError::Piston { id: 20, .. })
        ));
        assert_eq!(pistons, before);
        assert_eq!(thermal.cells(), &saved_heat);
        assert_eq!(thermal.time_s(), 0.0);
        let retry = run(&drives).unwrap();
        assert_eq!(retry.pistons(), candidate.pistons());
        assert_eq!(retry.thermal().cells(), candidate.thermal().cells());
        for p in pistons.values() {
            assert_eq!(
                p.piston().audit_last_window().unwrap(),
                ModuleWindowAudit::NotStepped
            );
        }
        gate.request();
        assert!(matches!(run(&drives), Err(MovingSolidGasError::Cancelled)));
        assert_eq!(pistons, before);
        assert_eq!(thermal.cells(), &saved_heat);
    });
}
