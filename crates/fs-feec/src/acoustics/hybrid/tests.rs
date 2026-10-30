use super::*;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
};
use crate::acoustics::helmholtz::boundary_layer::{BoundaryLayerGeometry, BoundaryLayerTransport};
use crate::acoustics::helmholtz::multiport::BoundaryMode;
use crate::acoustics::helmholtz::{HelmholtzFe, HelmholtzModes};
use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;
use crate::acoustics::thermoviscous::{ThermoviscousBoundary, ThermoviscousLimits};
use fs_rep_mesh::TetComplex;

fn gas() -> ThermoviscousMedium {
    ThermoviscousMedium {
        density_kg_m3: 101325.0 / (287.05 * 293.15),
        pressure_pa: 101325.0,
        temperature_k: 293.15,
        cp_j_kg_k: 1005.0,
        dynamic_viscosity_pa_s: 1.81e-5,
        thermal_conductivity_w_m_k: 0.0257,
        bulk_viscosity_pa_s: 1.1e-5,
    }
}
fn medium() -> Medium {
    let g = gas();
    Medium {
        density_kg_m3: g.density_kg_m3,
        sound_speed_m_s: (g.cp_j_kg_k / (g.cp_j_kg_k - g.gas_constant_j_kg_k()) * g.pressure_pa
            / g.density_kg_m3)
            .sqrt(),
    }
}
fn box_mesh(cells: [usize; 3], size: [f64; 3], offset: f64) -> (TetComplex, Vec<[f64; 3]>) {
    let index = |i, j, k| ((i * (cells[1] + 1) + j) * (cells[2] + 1) + k) as u32;
    let mut positions = Vec::new();
    for i in 0..=cells[0] {
        for j in 0..=cells[1] {
            for k in 0..=cells[2] {
                positions.push([
                    offset + size[0] * i as f64 / cells[0] as f64,
                    size[1] * j as f64 / cells[1] as f64,
                    size[2] * k as f64 / cells[2] as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..cells[0] {
        for j in 0..cells[1] {
            for k in 0..cells[2] {
                for order in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    let mut corner = [i, j, k];
                    let mut tet = [index(i, j, k); 4];
                    for (step, axis) in order.into_iter().enumerate() {
                        corner[axis] += 1;
                        tet[step + 1] = index(corner[0], corner[1], corner[2]);
                    }
                    tets.push(tet);
                }
            }
        }
    }
    (TetComplex::from_tets(positions.len(), tets), positions)
}
fn complete(fe: &HelmholtzFe<'_>) -> HelmholtzModes {
    let n = fe.ndof();
    HelmholtzModes {
        lambda: vec![0.0; n],
        frequency_hz: vec![0.0; n],
        residual: vec![0.0; n],
        shapes: (0..n)
            .map(|i| {
                let mut x = vec![0.0; n];
                x[i] = 1.0;
                x
            })
            .collect(),
        below_low: 0,
        below_high: n,
    }
}
fn ends(fe: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|f| {
            f.iter()
                .all(|&v| (fe.positions()[v as usize][0] - x).abs() < 1e-14)
        })
        .collect()
}
fn shapes(
    fe: &HelmholtzFe<'_>,
    start: f64,
    length: f64,
    width: f64,
    nonuniform: bool,
    lns: bool,
) -> Vec<BoundaryMode> {
    let left = ends(fe, start);
    let right = ends(fe, start + length);
    let mut modes = vec![BoundaryMode::uniform(&left), BoundaryMode::uniform(&right)];
    if nonuniform {
        let faces = if lns { left } else { right };
        modes.push(BoundaryMode {
            traces: faces
                .iter()
                .map(|&face| {
                    (
                        face,
                        face.map(|v| {
                            let y = fe.positions()[v as usize][1] / width;
                            if lns { 1.0 + 2.0 * y } else { y }
                        }),
                    )
                })
                .collect(),
        });
    }
    modes
}
fn boundary(
    fe: &HelmholtzFe<'_>,
    modes: &[BoundaryMode],
    isothermal: bool,
) -> ThermoviscousBoundary {
    let ends = modes
        .iter()
        .flat_map(|m| m.traces.iter().map(|(f, _)| *f))
        .collect::<BTreeSet<_>>();
    let walls = fe
        .boundary_faces()
        .into_iter()
        .filter(|f| !ends.contains(f))
        .collect::<Vec<_>>();
    ThermoviscousBoundary {
        no_slip_faces: walls.clone(),
        isothermal_faces: if isothermal { walls } else { vec![] },
        ports: modes.to_vec(),
    }
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs().max(1e-30)
}

// Independent mixed impedance/admittance closure, without power waves at the
// interface, the transformer S, or the hybrid condensation implementation.
fn independent(
    z: &[C64],
    y: &[C64],
    map: &[f64],
    h_ref: f64,
    l_ref: f64,
    incident: [C64; 2],
) -> Vec<C64> {
    let m = (map.len() as f64).sqrt() as usize;
    let n = m + 1;
    let hports = (1..n).collect::<Vec<_>>();
    let lports = (0..n).filter(|&p| p != 1).collect::<Vec<_>>();
    let mut a = vec![C64::ZERO; 4 * n * n];
    let mut rhs = vec![C64::ZERO; 2 * n];
    a[..n].copy_from_slice(&z[..n]);
    a[0] = a[0] + C64::from_re(h_ref);
    rhs[0] = incident[0].scale(2.0 * h_ref.sqrt());
    for i in 0..m {
        let row = 1 + i;
        a[row * 2 * n..row * 2 * n + n].copy_from_slice(&z[hports[i] * n..(hports[i] + 1) * n]);
        for j in 0..m {
            a[row * 2 * n + n + lports[j]] = C64::from_re(-map[j * m + i]);
        }
        let row = n + i;
        a[row * 2 * n + n..(row + 1) * 2 * n]
            .copy_from_slice(&y[lports[i] * n..(lports[i] + 1) * n]);
        for j in 0..m {
            a[row * 2 * n + hports[j]] = C64::from_re(map[i * m + j]);
        }
    }
    for j in 0..n {
        a[(2 * n - 1) * 2 * n + n + j] = y[n + j].scale(l_ref);
    }
    a[(2 * n - 1) * 2 * n + n + 1] = a[(2 * n - 1) * 2 * n + n + 1] + C64::ONE;
    rhs[2 * n - 1] = incident[1].scale(2.0 * l_ref.sqrt());
    lu_complex_with_cancel(&a, 2 * n, || false)
        .unwrap()
        .solve_with_cancel(&mut rhs, || false)
        .unwrap();
    rhs
}

#[test]
fn nonorthogonal_source_traces_close_original_fields_and_all_physical_work() {
    let size = [0.04, 0.02, 0.02];
    let (hm, hp) = box_mesh([2, 2, 2], size, 0.0);
    let (lm, lp) = box_mesh([2, 2, 2], [0.01, size[1], size[2]], size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = shapes(&hf, 0.0, size[0], size[1], true, false);
    let ls = shapes(&lf, size[0], 0.01, size[1], true, true);
    let hb = boundary(&hf, &hs, true);
    let lb = boundary(&lf, &ls, false);
    let g = gas();
    let geometry = BoundaryLayerGeometry {
        minimum_curvature_radius_m: f64::INFINITY,
        minimum_half_gap_m: 0.01,
        maximum_layer_ratio: 0.05,
    };
    let h = InteriorMultiport::new_with_boundary_layers(
        &hf,
        &complete(&hf),
        &hs,
        false,
        &hb.no_slip_faces,
        BoundaryLayerTransport {
            dynamic_viscosity_pa_s: g.dynamic_viscosity_pa_s,
            thermal_conductivity_w_m_k: g.thermal_conductivity_w_m_k,
            cp_j_kg_k: g.cp_j_kg_k,
            heat_capacity_ratio: g.cp_j_kg_k / (g.cp_j_kg_k - g.gas_constant_j_kg_k()),
        },
        geometry,
        &|| false,
    )
    .unwrap();
    let l = ThermoviscousFe::new(
        &lm,
        &lp,
        g,
        &lb,
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap()
    .with_conductive_walls(
        &[ConductiveWallRegion {
            faces: lb.no_slip_faces.clone(),
            temperature_k: g.temperature_k,
            geometry,
            wall: ConductiveThermalWall {
                transport: WallThermalTransport {
                    conductivity_w_m_k: 0.003,
                    volumetric_heat_capacity_j_m3_k: 2000.0,
                },
                depth: ThermalWallDepth::Slab {
                    thickness_m: 0.00002,
                    backing: ThermalBacking::Adiabatic,
                },
                minimum_tangential_scale_m: 0.004,
            },
        }],
        &mut || false,
    )
    .unwrap();
    let hybrid = HybridInterior::new(&h, &l, &[1, 2], &[0, 2], &|| false).unwrap();
    let expected = [1.0, 0.5, 2.0, 7.0 / 6.0];
    for (&actual, exact) in hybrid.interface().flow_map().iter().zip(expected) {
        assert!((actual - exact).abs() < 2e-14, "{actual} {exact}");
    }
    let f = 1700.0;
    let cp = C64::new(g.cp_j_kg_k, -25.0);
    let hr = [1.1e6, 0.8e6, 2.5e6];
    let lr = [1.7e6, 0.9e6, 3.1e6];
    let prepared = hybrid.prepare(f, cp, &hr, &lr, &|| false).unwrap();
    let incident = [C64::new(0.4, 0.1), C64::new(-0.3, 0.7)];
    let response = prepared.response(&incident, &|| false).unwrap();
    let power = response.power;
    assert!(power.helmholtz_wall.viscous_w > 0.0 && power.helmholtz_wall.thermal_w > 0.0);
    let lp = power.thermoviscous;
    assert!(
        lp.viscous_w > 0.0
            && lp.bulk_w > 0.0
            && lp.conductive_w > 0.0
            && lp.solid_conductive_w > 0.0
            && lp.caloric_w > 0.0
    );
    assert!(
        power.interface_residual_w().abs() < 1e-10 * power.external_input_w,
        "{power:?}"
    );
    assert!(
        power.balance_residual_w().abs() < 2e-8 * power.external_input_w,
        "{power:?}"
    );
    let pscale = response
        .helmholtz
        .pressure_pa
        .iter()
        .chain(&response.thermoviscous.pressure_drive_pa)
        .map(|z| z.abs())
        .fold(0.0_f64, f64::max);
    let qscale = response
        .helmholtz
        .inward_flow_m3_s
        .iter()
        .chain(&response.thermoviscous.inward_flow_m3_s)
        .map(|z| z.abs())
        .fold(0.0_f64, f64::max);
    assert!(
        response
            .interface_pressure_residual_pa
            .iter()
            .all(|z| z.abs() < 1e-10 * pscale)
    );
    assert!(
        response
            .interface_flow_residual_m3_s
            .iter()
            .all(|z| z.abs() < 1e-10 * qscale)
    );
    let z = h.impedance(f, &|| false).unwrap();
    let y = l
        .admittance_with_heat_capacity(f, cp, &mut || false)
        .unwrap();
    let oracle = independent(&z, &y, &expected, hr[0], lr[1], incident);
    for i in 0..3 {
        assert!(relative(response.helmholtz.inward_flow_m3_s[i], oracle[i]) < 2e-9);
        assert!(relative(response.thermoviscous.pressure_drive_pa[i], oracle[3 + i]) < 2e-9);
    }
    let s = prepared.scattering(&|| false).unwrap();
    assert!((s[1] - s[2]).abs() < 2e-10);
    for row in 0..2 {
        assert!(
            relative(
                response.outgoing[row],
                s[row * 2] * incident[0] + s[row * 2 + 1] * incident[1]
            ) < 2e-10
        );
    }
    let wave_work = 0.5
        * (incident.iter().map(|z| z.norm_sq()).sum::<f64>()
            - response.outgoing.iter().map(|z| z.norm_sq()).sum::<f64>());
    assert!((wave_work - power.dissipated_w()).abs() < 2e-8 * wave_work);
    let repeat = prepared.response(&incident, &|| false).unwrap();
    assert_eq!(
        response.thermoviscous.temperature_k,
        repeat.thermoviscous.temperature_k
    );
    assert_eq!(response.helmholtz.amplitudes, repeat.helmholtz.amplitudes);
    assert_eq!(s, prepared.scattering(&|| false).unwrap());
    // Internal reference resistance is only a coordinate choice, even with
    // unequal, nonorthogonal spatial bases and off-diagonal Gram entries.
    let other = hybrid
        .prepare(
            f,
            cp,
            &[hr[0], 3.3e6, 0.4e6],
            &[0.6e6, lr[1], 1.9e6],
            &|| false,
        )
        .unwrap()
        .scattering(&|| false)
        .unwrap();
    for (&a, b) in s.iter().zip(other) {
        assert!((a - b).abs() < 2e-10);
    }
}

#[test]
fn actual_lns_load_with_helmholtz_mesh_ladder_matches_independent_duct_propagation() {
    let length = 0.3;
    let size = [0.01, 0.03, 0.04];
    let (lm, lp) = box_mesh([2, 1, 1], size, length);
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let ls = shapes(&lf, length, size[0], size[1], false, true);
    let l = ThermoviscousFe::new(
        &lm,
        &lp,
        gas(),
        &boundary(&lf, &ls, true),
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    let f = 420.0;
    let cp = C64::from_re(gas().cp_j_kg_k);
    let r = medium().density_kg_m3 * medium().sound_speed_m_s / (size[1] * size[2]);
    let kl = std::f64::consts::TAU * f * length / medium().sound_speed_m_s;
    let zh = [
        C64::new(0.0, -r / kl.tan()),
        C64::new(0.0, -r / kl.sin()),
        C64::new(0.0, -r / kl.sin()),
        C64::new(0.0, -r / kl.tan()),
    ];
    let y = l
        .admittance_with_heat_capacity(f, cp, &mut || false)
        .unwrap();
    let incident = [C64::new(0.8, 0.1), C64::new(-0.2, 0.3)];
    let oracle = independent(&zh, &y, &[1.0], r, r, incident);
    let exact_pressure = zh[0] * oracle[0] + zh[1] * oracle[1];
    let mut errors = Vec::new();
    for nx in [2, 4, 8] {
        let (hm, hp) = box_mesh([nx, 1, 1], [length, size[1], size[2]], 0.0);
        let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
        let hs = shapes(&hf, 0.0, length, size[1], false, false);
        let h = InteriorMultiport::new(&hf, &complete(&hf), &hs, false, &|| false).unwrap();
        let hybrid = HybridInterior::new(&h, &l, &[1], &[0], &|| false).unwrap();
        let response = hybrid
            .prepare(f, cp, &[r, r], &[r, r], &|| false)
            .unwrap()
            .response(&incident, &|| false)
            .unwrap();
        let error = relative(response.pressure_pa[0], exact_pressure)
            .max(relative(response.inward_flow_m3_s[0], oracle[0]))
            .max(relative(response.pressure_pa[1], oracle[3]));
        eprintln!(
            "hybrid Helmholtz axial cells={nx} dofs={} independent load/propagation error={error:.15e}",
            h.dim()
        );
        errors.push(error);
    }
    assert!(
        errors[1] < errors[0] / 5.0 && errors[2] < errors[1] / 5.0,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-4, "{errors:?}");
}

#[test]
fn source_trace_geometry_rank_coordinate_and_cancellation_refusals() {
    let size = [0.02, 0.02, 0.02];
    let (hm, hp) = box_mesh([1, 1, 1], size, 0.0);
    let (lm, lp) = box_mesh([1, 1, 1], size, size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = shapes(&hf, 0.0, size[0], size[1], false, false);
    let ls = shapes(&lf, size[0], size[0], size[1], false, true);
    let h = InteriorMultiport::new(&hf, &complete(&hf), &hs, false, &|| false).unwrap();
    let l = ThermoviscousFe::new(
        &lm,
        &lp,
        gas(),
        &boundary(&lf, &ls, false),
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    for (hi, li) in [
        (vec![], vec![]),
        (vec![1], vec![0, 1]),
        (vec![2], vec![0]),
        (vec![1, 1], vec![0, 1]),
        (vec![0], vec![0]),
        (vec![0, 1], vec![0, 1]),
    ] {
        assert!(HybridInterior::new(&h, &l, &hi, &li, &|| false).is_err());
    }
    assert!(matches!(
        HybridInterior::new(&h, &l, &[1], &[0], &|| true),
        Err(HybridError::Cancelled)
    ));
    let duplicate_h = [hs[0].clone(), hs[1].clone(), hs[1].clone()];
    let duplicate_l = [ls[0].clone(), ls[1].clone(), ls[0].clone()];
    let dh = InteriorMultiport::new(&hf, &complete(&hf), &duplicate_h, false, &|| false).unwrap();
    let dl = ThermoviscousFe::new(
        &lm,
        &lp,
        gas(),
        &boundary(&lf, &duplicate_l, false),
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    assert!(HybridInterior::new(&dh, &dl, &[1, 2], &[0, 2], &|| false).is_err());
    // Reusing coincident geometry with the same outward side must refuse.
    let same_l = ThermoviscousFe::new(
        &hm,
        &hp,
        gas(),
        &boundary(&hf, &hs, false),
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    assert!(HybridInterior::new(&h, &same_l, &[1], &[1], &|| false).is_err());
    let hybrid = HybridInterior::new(&h, &l, &[1], &[0], &|| false).unwrap();
    assert_eq!(hybrid.external_port_count(), 2);
    let cp = C64::from_re(gas().cp_j_kg_k);
    assert!(
        hybrid
            .prepare(0.0, cp, &[1e6; 2], &[1e6; 2], &|| false)
            .is_err()
    );
    assert!(
        hybrid
            .prepare(400.0, cp, &[1e6], &[1e6; 2], &|| false)
            .is_err()
    );
    assert!(
        hybrid
            .prepare(400.0, cp, &[1e6; 2], &[-1.0, 1e6], &|| false)
            .is_err()
    );
    assert!(
        hybrid
            .prepare(400.0, C64::new(1005.0, 1.0), &[1e6; 2], &[1e6; 2], &|| {
                false
            })
            .is_err()
    );
    let prepared = hybrid
        .prepare(400.0, cp, &[1e6; 2], &[1e6; 2], &|| false)
        .unwrap();
    assert!(prepared.response(&[C64::ONE], &|| false).is_err());
    assert!(
        prepared
            .response(&[C64::new(f64::NAN, 0.0), C64::ONE], &|| false)
            .is_err()
    );
    assert!(matches!(
        prepared.response(&[C64::ONE; 2], &|| true),
        Err(HybridError::Cancelled)
    ));
    assert!(matches!(
        prepared.scattering(&|| true),
        Err(HybridError::Cancelled)
    ));
}
