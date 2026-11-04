//! G1/G3/G4: real acceleration, sourced heat feedback and one native mass owner.

use super::*;
use fs_lbm::d3q19::freesurface3::{
    ConversionRedistribution3, ConversionStats3, FreeSurfaceSample3,
};
use fs_lbm::d3q19::{E3, Q3, W3};
use fs_lbm::material::free_surface::MaterialRigidFrame3;
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_material::liquid::ResolvedLiquidState;

const DX: f64 = 0.001;
const DT: f64 = 0.001;
const CELL: [usize; 3] = [3; 3];

#[derive(Clone, Copy)]
struct Case {
    full: bool,
    buoyant: bool,
    gravity: [f64; 3],
    velocity: [f64; 3],
    controls: FreeSurfaceEnthalpyControls,
}

impl Default for Case {
    fn default() -> Self {
        Self {
            full: true,
            buoyant: false,
            gravity: [0.0; 3],
            velocity: [0.02, 0.0, 0.0],
            controls: FreeSurfaceEnthalpyControls::default(),
        }
    }
}

fn with_liquid(
    card: &MaterialCard,
    case: Case,
    h: impl Fn([usize; 3]) -> f64,
    exercise: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate),
) {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut cells = initial(case.full);
    for cell in &mut cells {
        cell.velocity_m_s = case.velocity;
    }
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells(
                MaterialFreeSurfaceDomain3 {
                    gravity_m_s2: case.gravity,
                    ..domain()
                },
                &cells,
                cx,
            )
            .unwrap()
            .with_conversion_redistribution(ConversionRedistribution3::ConservativeLocal);
        let mut liquid = if case.buoyant {
            flow.with_boussinesq_enthalpy(
                &law,
                &chart,
                case.controls,
                BoussinesqConfig {
                    gravity_m_s2: case.gravity,
                    max_relative_density_departure: 0.03,
                    max_relative_specific_heat_departure: 0.0,
                },
                h,
                &gate,
            )
        } else {
            flow.with_enthalpy(&law, &chart, case.controls, h, &gate)
        }
        .unwrap();
        exercise(&mut liquid, &gate);
    });
}

fn full_cell_heat(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    specific: f64,
) -> Vec<FreeSurfaceCellHeat> {
    initial(true)
        .iter()
        .map(|cell| {
            let [x, y, z] = cell.coordinates;
            let sample = liquid.flow().sample_at_si(x, y, z).unwrap().unwrap();
            FreeSurfaceCellHeat {
                coordinates: cell.coordinates,
                energy_j: sample.density_kg_m3 * DX.powi(3) * specific,
            }
        })
        .collect()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Solve the dimensional midpoint equation by elimination, independently of
// the production closed-form inverse and all production SI conversion helpers.
fn acceleration(
    frame: MaterialRigidFrame3,
    body: [f64; 3],
    p: [usize; 3],
    u: [f64; 3],
) -> [f64; 3] {
    let r = core::array::from_fn(|d| (p[d] as f64 + 0.5) * DX - frame.origin_m[d]);
    let euler = cross(frame.angular_acceleration_rad_s2, r);
    let centrifugal = cross(
        frame.angular_velocity_rad_s,
        cross(frame.angular_velocity_rad_s, r),
    );
    let coriolis = cross(frame.angular_velocity_rad_s, u);
    let [wx, wy, wz] = frame.angular_velocity_rad_s.map(|w| w * DT);
    let mut matrix = [[1.0, -wz, wy], [wz, 1.0, -wx], [-wy, wx, 1.0]];
    let mut rhs: [f64; 3] = core::array::from_fn(|d| {
        body[d] - frame.origin_acceleration_m_s2[d] - euler[d] - centrifugal[d] - 2.0 * coriolis[d]
    });
    for pivot in 0..3 {
        let scale = matrix[pivot][pivot];
        for column in pivot..3 {
            matrix[pivot][column] /= scale;
        }
        rhs[pivot] /= scale;
        for row in 0..3 {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            for column in pivot..3 {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    rhs
}

fn projection(q: usize, u: [f64; 3], force: [f64; 3]) -> f64 {
    let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
    let eu: f64 = e.iter().zip(u).map(|(e, u)| e * u).sum();
    W3[q]
        * (0..3)
            .map(|d| (3.0 * (e[d] - u[d]) + 9.0 * eu * e[d]) * force[d])
            .sum::<f64>()
}

fn population(bits: &[u64], p: [usize; 3], q: usize) -> f64 {
    // This full 8^3 fixture contains all eight Morton-ordered 4^3 tiles.
    let [x, y, z] = p;
    let slot = (x / 4) | ((y / 4) << 1) | ((z / 4) << 2);
    let lane = x % 4 + 4 * (y % 4 + 4 * (z % 4));
    f64::from_bits(bits[(slot * Q3 + q) * 64 + lane])
}

#[test]
fn g1_heated_rotation_matches_independent_si_midpoint_and_incoming_bgk_guo_populations() {
    let source = card(0.02, 1000.0);
    with_liquid(
        &source,
        Case::default(),
        |_| 200000.0,
        |liquid, gate| {
            let frame = MaterialRigidFrame3 {
                origin_m: [2.5 * DX, 4.5 * DX, 1.5 * DX],
                angular_velocity_rad_s: [1.5, -2.0, 4.0],
                angular_acceleration_rad_s2: [0.05, -0.04, 0.2],
                origin_acceleration_m_s2: [0.1, -0.08, 0.03],
            };
            let body = [0.2, -0.1, 0.05];
            liquid
                .step_in_rigid_frame_si(frame, body, &[], gate)
                .unwrap();
            assert_eq!(
                liquid.flow().lattice().force_density_at(3, 3, 3),
                Some([0.0; 3])
            );
            let bits = liquid.flow().lattice().grid().state_bits();
            let force_at = |p: [usize; 3]| {
                let sample = liquid.flow().lattice().sample_at(p[0], p[1], p[2]).unwrap();
                acceleration(frame, body, p, sample.velocity.map(|u| u * DX / DT))
                    .map(|a| sample.density * a * DT * DT / DX)
            };
            let center_force = force_at(CELL);
            let incoming = |tau: f64| -> [f64; Q3] {
                core::array::from_fn(|q| {
                    let e = [E3[q].0, E3[q].1, E3[q].2];
                    let p = core::array::from_fn(|d| (CELL[d] as i32 - e[d]) as usize);
                    let sample = liquid.flow().lattice().sample_at(p[0], p[1], p[2]).unwrap();
                    let u = sample.velocity;
                    let eu: f64 = (0..3).map(|d| f64::from(e[d]) * u[d]).sum();
                    let eq = W3[q]
                        * sample.density
                        * (1.0 + 3.0 * eu + 4.5 * eu * eu
                            - 1.5 * u.iter().map(|v| v * v).sum::<f64>());
                    let old = population(&bits, p, q);
                    // Initial endpoint force is zero. Reframing plus one BGK/Guo
                    // collision leaves half of the new Guo source in the outgoing PDF.
                    old - (old - eq) / tau + 0.5 * projection(q, u, force_at(p))
                })
            };
            // Uniform +100 kJ/kg takes T=400 to 450 K: mu=.04, tau=.62.
            let streamed = incoming(0.62);
            assert!(
                streamed
                    .iter()
                    .zip(incoming(0.68))
                    .any(|(hot, cold)| (hot - cold).abs() > 1e-12)
            );
            let rho: f64 = streamed.iter().sum();
            let u: [f64; 3] = core::array::from_fn(|d| {
                (streamed
                    .iter()
                    .enumerate()
                    .map(|(q, f)| f * f64::from([E3[q].0, E3[q].1, E3[q].2][d]))
                    .sum::<f64>()
                    + 0.5 * center_force[d])
                    / rho
            });
            let expected: [f64; Q3] =
                core::array::from_fn(|q| streamed[q] + 0.5 * projection(q, u, center_force));
            let heat = full_cell_heat(liquid, 100000.0);
            let report = liquid
                .step_in_rigid_frame_si(frame, body, &heat, gate)
                .unwrap();
            let actual = liquid.flow().lattice().sample_at(3, 3, 3).unwrap();
            assert_eq!(
                liquid.flow().lattice().force_density_at(3, 3, 3),
                Some([0.0; 3])
            );
            near(actual.density, rho, 3e-13);
            for d in 0..3 {
                assert!((actual.velocity[d] - u[d]).abs() < 3e-13);
            }
            let after = liquid.flow().lattice().grid().state_bits();
            for q in 0..Q3 {
                assert!((population(&after, CELL, q) - expected[q]).abs() < 3e-13);
            }
            let material = liquid.material_at(CELL).unwrap();
            near(material.dynamic_viscosity_pa_s(), 0.04, 2e-11);
            assert_eq!(material.material().card_identity(), source.content_hash());
            assert!(
                material
                    .material()
                    .query_point()
                    .contains(&("pressure".into(), 100000.0))
            );
            near(
                report.enthalpy_after_j - report.enthalpy_before_j,
                heat.iter().map(|q| q.energy_j).sum(),
                2e-11,
            );
            near(report.mass_after_kg, report.mass_before_kg, 2e-11);
            assert_eq!(report.viscous_heat_j, 0.0);
            assert_eq!(liquid.accepted_steps(), 2);
            assert_eq!(liquid.physical_time_s(), 2.0 * DT);
        },
    );
}

#[test]
fn g3_source_density_uses_effective_gravity_and_free_fall_creates_no_molecular_heat() {
    let source = card(0.02, 980.0);
    let gravity = [0.2, -0.1, -1.0];
    with_liquid(
        &source,
        Case {
            buoyant: true,
            gravity,
            velocity: [0.0; 3],
            ..Case::default()
        },
        |_| 200000.0,
        |liquid, gate| {
            let freefall = MaterialRigidFrame3 {
                origin_acceleration_m_s2: gravity,
                ..Default::default()
            };
            for specific_heat in [0.0, 100000.0] {
                let heat = full_cell_heat(liquid, specific_heat);
                let report = liquid
                    .step_in_rigid_frame_with_viscous_heating_si(freefall, gravity, &heat, gate)
                    .unwrap();
                assert!(report.flow.viscous_heat_j.abs() < 1e-20);
                assert!(
                    report
                        .applied
                        .iter()
                        .all(|cell| cell.shear_rate_per_s.abs() < 1e-8)
                );
                let sample = liquid.flow().sample_at_si(3, 3, 3).unwrap().unwrap();
                assert!(sample.velocity_m_s.iter().all(|u| u.abs() < 2e-14));
                let departure =
                    (liquid.material_at(CELL).unwrap().density_kg_m3() - 1000.0) / 1000.0;
                let original_force =
                    gravity.map(|g| g * DT / (DX / DT) + departure * g * DT / (DX / DT));
                assert_eq!(
                    liquid.flow().lattice().force_density_at(3, 3, 3),
                    Some(original_force)
                );
                near(
                    report.flow.enthalpy_after_j,
                    report.flow.enthalpy_before_j + report.flow.external_heat_j,
                    2e-11,
                );
            }
            near(
                liquid.material_at(CELL).unwrap().density_kg_m3(),
                985.0,
                2e-11,
            );
            near(
                liquid.material_at(CELL).unwrap().dynamic_viscosity_pa_s(),
                0.04,
                2e-11,
            );
            let before = liquid
                .flow()
                .sample_at_si(3, 3, 3)
                .unwrap()
                .unwrap()
                .velocity_m_s;
            let body = [0.35, -0.15, 0.08];
            let origin = [0.12, 0.2, -0.05];
            let report = liquid
                .step_in_rigid_frame_si(
                    MaterialRigidFrame3 {
                        origin_acceleration_m_s2: origin,
                        ..Default::default()
                    },
                    body,
                    &[],
                    gate,
                )
                .unwrap();
            let after = liquid
                .flow()
                .sample_at_si(3, 3, 3)
                .unwrap()
                .unwrap()
                .velocity_m_s;
            for d in 0..3 {
                assert!(
                    (after[d] - before[d] - DT * (985.0 / 1000.0) * (body[d] - origin[d])).abs()
                        < 3e-14
                );
            }
            // The accepted endpoint retains the original thermal force datum,
            // with no additional physical kick. A subsequent frame must be explicit.
            let force = liquid.flow().lattice().force_density_at(3, 3, 3).unwrap();
            for d in 0..3 {
                near(force[d], gravity[d] * DT * DT / DX * 0.985, 2e-11);
            }
            near(report.mass_after_kg, report.mass_before_kg, 2e-11);
            assert_eq!(report.external_heat_j, 0.0);
            assert_eq!(report.viscous_heat_j, 0.0);
        },
    );
}

#[test]
fn g1_frame_driven_interface_conversion_carries_actual_source_enthalpy_once() {
    let source = card(0.02, 1000.0);
    with_liquid(
        &source,
        Case {
            full: false,
            ..Case::default()
        },
        |p| 100000.0 + 50000.0 * p[0] as f64,
        |liquid, gate| {
            let frame = MaterialRigidFrame3 {
                angular_velocity_rad_s: [1.0, -0.5, 2.0],
                origin_acceleration_m_s2: [-0.1, 0.02, 0.0],
                ..Default::default()
            };
            let heat = [FreeSurfaceCellHeat {
                coordinates: [2, 3, 3],
                energy_j: 0.02,
            }];
            let report = liquid
                .step_in_rigid_frame_si(frame, [0.0; 3], &heat, gate)
                .unwrap();
            let conversions = liquid.flow().lattice().conversions();
            assert!(conversions.to_fluid > 0 && conversions.gas_to_interface > 0);
            assert!(report.mass_routes > 0);
            near(report.mass_after_kg, report.mass_before_kg, 2e-11);
            near(
                report.enthalpy_after_j,
                report.enthalpy_before_j + 0.02,
                2e-11,
            );
            assert_eq!(report.viscous_heat_j, 0.0);
            let h = liquid.specific_enthalpy_at([4, 3, 3]).unwrap();
            let material = liquid.material_at([4, 3, 3]).unwrap();
            assert!((150000.0..=270000.0).contains(&h));
            near(
                material.dynamic_viscosity_pa_s(),
                0.1 - 0.08 * h / 400000.0,
                2e-11,
            );
            assert_eq!(material.material().card_identity(), source.content_hash());
        },
    );
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    populations: Vec<u64>,
    cells: Vec<(
        Option<FreeSurfaceSample3>,
        Option<[f64; 3]>,
        Option<f64>,
        Option<ResolvedLiquidState>,
    )>,
    conversions: ConversionStats3,
    mass: f64,
    enthalpy: f64,
    steps: u64,
    time: f64,
}

fn snapshot(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> Snapshot {
    let lattice = liquid.flow().lattice();
    Snapshot {
        populations: lattice.grid().state_bits(),
        cells: initial(true)
            .iter()
            .map(|cell| {
                let [x, y, z] = cell.coordinates;
                (
                    lattice.sample_at(x, y, z),
                    lattice.force_density_at(x, y, z),
                    liquid.specific_enthalpy_at(cell.coordinates),
                    liquid.material_at(cell.coordinates).cloned(),
                )
            })
            .collect(),
        conversions: lattice.conversions(),
        mass: liquid.flow().ledger_mass_kg().unwrap(),
        enthalpy: liquid.total_enthalpy_j().unwrap(),
        steps: liquid.accepted_steps(),
        time: liquid.physical_time_s(),
    }
}

#[test]
fn g4_frame_source_rotation_cancellation_and_late_transport_refusals_are_atomic() {
    let source = card(0.02, 980.0);
    with_liquid(
        &source,
        Case {
            buoyant: true,
            ..Case::default()
        },
        |_| 200000.0,
        |liquid, gate| {
            let before = snapshot(liquid);
            for frame in [
                MaterialRigidFrame3 {
                    angular_velocity_rad_s: [0.0, 1.0, 0.0],
                    ..Default::default()
                },
                MaterialRigidFrame3 {
                    angular_acceleration_rad_s2: [1.0, 0.0, 0.0],
                    ..Default::default()
                },
            ] {
                assert!(matches!(
                    liquid.step_in_rigid_frame_si(frame, [0.0; 3], &[], gate),
                    Err(FreeSurfaceEnthalpyError::Invalid(
                        "rotating thermal buoyancy requires density-dependent inertial coupling"
                    ))
                ));
                assert_eq!(snapshot(liquid), before);
            }
            let frame = MaterialRigidFrame3::default();
            assert!(
                liquid
                    .step_in_rigid_frame_si(frame, [f64::NAN, 0.0, 0.0], &[], gate)
                    .is_err()
            );
            let heat = [FreeSurfaceCellHeat {
                coordinates: CELL,
                energy_j: 0.3,
            }];
            assert!(matches!(
                liquid.step_in_rigid_frame_si(frame, [0.0; 3], &heat, gate),
                Err(FreeSurfaceEnthalpyError::Material(_))
            ));
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(
                liquid
                    .step_in_rigid_frame_si(frame, [0.0; 3], &[], &cancelled)
                    .is_err()
            );
            assert_eq!(snapshot(liquid), before);
        },
    );
    let fixed = card(0.02, 1000.0);
    with_liquid(
        &fixed,
        Case {
            full: false,
            controls: FreeSurfaceEnthalpyControls {
                max_mass_routes: 1,
                ..FreeSurfaceEnthalpyControls::default()
            },
            ..Case::default()
        },
        |_| 200000.0,
        |liquid, gate| {
            let before = snapshot(liquid);
            let frame = MaterialRigidFrame3 {
                origin_acceleration_m_s2: [0.1, 0.0, 0.0],
                ..Default::default()
            };
            let error = liquid
                .step_in_rigid_frame_si(frame, [0.0; 3], &[], gate)
                .unwrap_err();
            assert!(
                format!("{error:?}").contains("route"),
                "expected native route-budget refusal: {error:?}"
            );
            assert_eq!(snapshot(liquid), before);
        },
    );
}
