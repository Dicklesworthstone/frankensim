use super::*;
use crate::acoustics::helmholtz::QuadraticGeometryControls;
use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
};
use crate::acoustics::helmholtz::tri_quadrature;
use crate::acoustics::thermoviscous::time::{
    ThermoviscousCaloricReservoir, ThermoviscousTimeLimits,
};
use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const LENGTH: f64 = 0.002;
const BOW: f64 = 0.35;
const REFERENCE_POINTS: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [LENGTH, 0.0, 0.0],
    [0.0, LENGTH, 0.0],
    [0.0, 0.0, LENGTH],
    [0.25 * LENGTH, 0.25 * LENGTH, 0.25 * LENGTH],
];

fn subdivided_tet(reverse: bool) -> TetComplex {
    let mut tets = vec![[0, 1, 2, 4], [0, 3, 1, 4], [0, 2, 3, 4], [1, 3, 2, 4]];
    if reverse {
        for tet in &mut tets {
            tet.swap(1, 2);
        }
    }
    TetComplex::from_tets(REFERENCE_POINTS.len(), tets)
}

fn with_context(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x544845524d414c,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

// Exact polynomial map (X,Y,Z) -> L*(X,Y,Z*(1+a*X)); det J=L^3*(1+a*X).
// Its volume and source moments below are analytic simplex integrals, not
// values recovered from the FE matrices under test.
fn map(mut p: [f64; 3]) -> [f64; 3] {
    p[2] *= 1.0 + BOW * p[0] / LENGTH;
    p
}

fn mids(mesh: &TetComplex) -> Vec<[f64; 3]> {
    mesh.edges
        .iter()
        .map(|e| {
            map(std::array::from_fn(|axis| {
                f64::midpoint(
                    REFERENCE_POINTS[e[0] as usize][axis],
                    REFERENCE_POINTS[e[1] as usize][axis],
                )
            }))
        })
        .collect()
}

fn curved<'c>(
    mesh: &'c TetComplex,
    order: usize,
    isothermal: bool,
    cx: &Cx<'_>,
) -> ThermoviscousFe<'c> {
    let gas = air();
    let medium = Medium {
        density_kg_m3: gas.density_kg_m3,
        sound_speed_m_s: ThermoviscousThermodynamics::ideal_gas(gas)
            .unwrap()
            .sound_speed_m_s(gas)
            .unwrap(),
    };
    let geometry = HelmholtzFe::new_quadratic(
        mesh,
        &REFERENCE_POINTS.map(map),
        &mids(mesh),
        2,
        order,
        medium,
        QuadraticGeometryControls::default(),
        cx,
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    ThermoviscousFe::from_geometry(
        geometry,
        gas,
        &ThermoviscousBoundary {
            no_slip_faces: faces.clone(),
            isothermal_faces: if isothermal { faces } else { vec![] },
            ports: vec![],
        },
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn curved_heat_release_uses_physical_positions_variable_volume_and_declared_order() {
    with_context(|cx| {
        let mesh = subdivided_tet(false);
        let model = curved(&mesh, 6, false, cx);
        let constant = model
            .volume_heat_source(6, &mut |_, _| Ok(LENGTH.powi(-3)), &mut || false)
            .unwrap();
        close(constant.integrated_shape(), 1.0 / 6.0 + BOW / 24.0, 3e-13);
        let mut calls = 0;
        let source = model
            .volume_heat_source(
                6,
                &mut |element, point| {
                    assert!(element < mesh.tets.len());
                    calls += 1;
                    Ok(point[2].powi(2) / LENGTH.powi(5))
                },
                &mut || false,
            )
            .unwrap();
        assert_eq!(calls, mesh.tets.len() * 6_usize.pow(3));
        let moment = 1.0 / 60.0 + BOW / 120.0 + BOW.powi(2) / 420.0 + BOW.powi(3) / 3360.0;
        close(source.integrated_shape(), moment, 3e-13);
        let vertex_load = (0..mesh.vertex_count)
            .map(|vertex| source.load[model.temperature_map[vertex].unwrap()])
            .sum::<f64>();
        close(vertex_load, moment, 3e-13);
        let mut rejected_calls = 0;
        assert!(
            model
                .volume_heat_source(
                    5,
                    &mut |_, _| {
                        rejected_calls += 1;
                        Ok(1.0)
                    },
                    &mut || false,
                )
                .is_err()
        );
        assert_eq!(rejected_calls, 0);
        assert!(matches!(
            model.volume_heat_source(6, &mut |_, _| Ok(1.0), &mut || true),
            Err(ThermoviscousError::Cancelled)
        ));
    });
}

#[test]
fn curved_original_molecular_field_has_analytic_storage_and_exact_exchange_decay() {
    with_context(|cx| {
        let mesh = subdivided_tet(true);
        let model = curved(&mesh, 6, false, cx);
        let reservoir = ThermoviscousCaloricReservoir {
            heat_capacity_j_kg_k: 100.0,
            relaxation_time_s: 0.001,
        };
        let dt = 0.3 * reservoir.relaxation_time_s;
        let prepared = model
            .prepare_time_with_caloric_reservoirs(
                dt,
                air().cp_j_kg_k - reservoir.heat_capacity_j_kg_k,
                &[reservoir],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut initial = prepared.zero_state();
        let temperature = 0.7;
        initial.molecular_temperature_k[0][..mesh.vertex_count].fill(temperature);
        let volume = LENGTH.powi(3) * (1.0 / 6.0 + BOW / 24.0);
        let energy = air().density_kg_m3 * reservoir.heat_capacity_j_kg_k
            / (2.0 * air().temperature_k)
            * volume
            * temperature
            * temperature;
        close(prepared.stored_energy_j(&initial).unwrap(), energy, 2e-12);
        let snapshot = initial.clone();
        let step = prepared.step(&initial, &[], &[], &mut || false).unwrap();
        assert_eq!(initial, snapshot);
        // Uniform insulated exchange preserves density. The ORIGINAL frozen
        // Cv and molecular capacity exchange heat, so both temperatures move;
        // this is not a clamped thermal field or a fitted decay time.
        let capacity = reservoir.heat_capacity_j_kg_k;
        let frozen_cv = air().cp_j_kg_k - capacity - air().gas_constant_j_kg_k();
        let decay = det::exp(-(1.0 + capacity / frozen_cv) * dt / reservoir.relaxation_time_s);
        let fluid_temperature = capacity / (frozen_cv + capacity) * temperature * (1.0 - decay);
        let molecular_temperature =
            (capacity + frozen_cv * decay) / (frozen_cv + capacity) * temperature;
        for (i, (&fluid, &molecular)) in step
            .state
            .temperature_k
            .iter()
            .zip(&step.state.molecular_temperature_k[0])
            .enumerate()
        {
            let original_vertex = i < mesh.vertex_count;
            let expected_fluid = if original_vertex {
                fluid_temperature
            } else {
                0.0
            };
            let expected_molecular = if original_vertex {
                molecular_temperature
            } else {
                0.0
            };
            assert!((fluid - expected_fluid).abs() < 3e-11);
            assert!((molecular - expected_molecular).abs() < 3e-11);
        }
        for &pressure in &step.state.pressure_pa {
            close(
                pressure,
                air().pressure_pa / air().temperature_k * fluid_temperature,
                3e-10,
            );
        }
        let final_energy = air().density_kg_m3 * volume / (2.0 * air().temperature_k)
            * (frozen_cv * fluid_temperature.powi(2) + capacity * molecular_temperature.powi(2));
        close(step.work.final_energy_j, final_energy, 3e-10);
        close(
            step.work.molecular_relaxation_heat_j,
            energy - final_energy,
            3e-10,
        );
        assert!(step.work.balance_residual_j().abs() < 3e-10 * energy);
    });
}

fn solid_region() -> ConductiveWallRegion {
    ConductiveWallRegion {
        faces: vec![[1, 2, 3]],
        wall: ConductiveThermalWall {
            transport: WallThermalTransport {
                conductivity_w_m_k: 0.003,
                volumetric_heat_capacity_j_m3_k: 2000.0,
            },
            depth: ThermalWallDepth::Slab {
                thickness_m: 0.00002,
                backing: ThermalBacking::Adiabatic,
            },
            minimum_tangential_scale_m: LENGTH,
        },
        geometry: BoundaryLayerGeometry {
            minimum_curvature_radius_m: LENGTH / (2.0 * BOW),
            minimum_half_gap_m: 0.0001,
            maximum_layer_ratio: 0.05,
        },
        temperature_k: air().temperature_k,
    }
}

fn temperature(point: [f64; 3]) -> f64 {
    0.3 * (1.0 + 0.1 * point[0] / LENGTH + 0.2 * point[2] / LENGTH)
}

#[test]
fn curved_wall_pointwise_heat_matches_original_form_and_physical_surface_ladder() {
    with_context(|cx| {
        // Independent graph-area formula for z/L=(1-X-Y)*(1+a*X).
        let reference = tri_quadrature(16)
            .iter()
            .map(|&(b, weight)| {
                let [x, _, z] = b;
                let jacobian = LENGTH
                    * LENGTH
                    * (1.0 + (BOW * z - 1.0 - BOW * x).powi(2) + (1.0 + BOW * x).powi(2)).sqrt();
                let temp = 0.3 * (1.0 + 0.1 * x + 0.2 * z * (1.0 + BOW * x));
                weight * jacobian * temp * temp
            })
            .sum::<f64>();
        for reverse in [false, true] {
            let mesh = subdivided_tet(reverse);
            let positions = REFERENCE_POINTS.map(map);
            let midpoints = mids(&mesh);
            let mut errors = Vec::new();
            for order in [4, 6, 8] {
                let model = curved(&mesh, order, false, cx)
                    .with_conductive_walls(&[solid_region()], &mut || false)
                    .unwrap();
                let mut temperatures = vec![C64::ZERO; model.space.ndof];
                for (i, &p) in positions.iter().enumerate() {
                    temperatures[i] = C64::from_re(temperature(p));
                }
                for (i, edge) in mesh.edges.iter().enumerate() {
                    temperatures[model.space.edge_off + i] = C64::from_re(
                        4.0 * temperature(midpoints[i])
                            - 2.0
                                * (temperature(positions[edge[0] as usize])
                                    + temperature(positions[edge[1] as usize])),
                    );
                }
                let admittances = model
                    .conductive
                    .admittances(model.fe.medium(), 1700.0, &mut || false)
                    .unwrap();
                let actual = model
                    .conductive
                    .dissipated_w(
                        &temperatures,
                        air().temperature_k,
                        &admittances,
                        &mut || false,
                    )
                    .unwrap();
                let mut mixed = vec![C64::ZERO; model.n];
                for (original, &offset) in model.temperature_map.iter().enumerate() {
                    if let Some(offset) = offset {
                        mixed[model.temperature_start + offset] = temperatures[original];
                    }
                }
                let form = model
                    .conductive
                    .entries
                    .iter()
                    .map(|entry| {
                        0.5 / air().temperature_k
                            * admittances[entry.region].re
                            * entry.mass
                            * (mixed[entry.row].conj() * mixed[entry.col]).re
                    })
                    .sum::<f64>();
                close(actual, form, 4e-12);
                let expected = admittances[0].re / (2.0 * air().temperature_k) * reference;
                let error = (actual - expected).abs() / expected.abs();
                assert!(error < 1e-8, "order={order}, error={error:e}");
                errors.push(error);
                let source = model
                    .volume_heat_source(order, &mut |_, _| Ok(LENGTH.powi(-3)), &mut || false)
                    .unwrap();
                let field = model
                    .prepare(1700.0, &mut || false)
                    .unwrap()
                    .response_with_heat_sources(
                        &[],
                        &[],
                        &[(&source, C64::from_re(1e-5))],
                        &mut || false,
                    )
                    .unwrap();
                assert_work(&field);
                assert!(field.power.solid_conductive_w > 0.0);
            }
            assert!(errors[1] < errors[0] + 1e-13, "{errors:?}");
            assert!(errors[2] < 2e-12, "{errors:?}");
        }
    });
}
