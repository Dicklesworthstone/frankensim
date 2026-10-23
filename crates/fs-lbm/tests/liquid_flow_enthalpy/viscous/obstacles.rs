//! Fixed no-slip solids use the original material stress and enthalpy owners.
use super::*;
use fs_lbm::conjugate::{Voxel, creeping_flow_with_body_force_heat};

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 100.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 100.0,
    }
}

#[test]
fn material_wall_work_above_a_solid_bed_reaches_only_the_original_fluid_cells() {
    let mut domain = VoxelDomain::new(2, 3, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    domain.set(1, 0, 0, Voxel::Solid(0));
    let config = driven_box();
    let gate = CancelGate::new();
    let ids = [
        None,
        None,
        Some(CellId(2)),
        Some(CellId(3)),
        Some(CellId(4)),
        Some(CellId(5)),
    ];
    for mu in [100.0, 200.0] {
        let card = source(mu, false);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        // Retained, unmapped thermal cells have no fluid source or transfer.
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_viscous_heating(controls())
            .unwrap();
        let (next, accepted, thermal) = old
            .advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
            .unwrap();
        let epoch = accepted.applied().unwrap();
        let power = epoch.viscous_power.as_ref().unwrap();
        // The two fluid rows above the impermeable bed have one circulation q.
        // Independently minimizing D=mu*(24q²-4Uq+2U²) gives q=U/12.
        // The moving lid supplies 2mu*U(U-q)=11mu*U²/6.
        near(epoch.flow.field.flux_x(1, 1, 0), -1.0 / 12.0);
        near(epoch.flow.field.flux_x(1, 2, 0), 1.0 / 12.0);
        near(epoch.flow.field.flux_y(0, 2, 0), 1.0 / 12.0);
        near(epoch.flow.field.flux_y(1, 2, 0), -1.0 / 12.0);
        for x in 0..2 {
            assert_eq!(epoch.flow.field.flux_y(x, 1, 0), 0.0);
        }
        assert_eq!(epoch.flow.field.flux_x(1, 0, 0), 0.0);
        near(power.heat_w, 11.0 * mu / 6.0);
        near(power.wall_power_w[3], power.heat_w);
        near(power.supplied_wall_power_w, power.heat_w);
        near(power.pressure_constraint_power_w, 0.0);
        near(power.residual_w, 0.0);
        assert_eq!(epoch.samples.len(), 4);
        for c in 0..6 {
            let expected = match c {
                0 | 1 => 0.0,
                2 | 3 => mu / 24.0,
                _ => 7.0 * mu / 8.0,
            };
            near(power.cell_heat_w[c], expected);
            near(
                next.cells()[&CellId(c as u64)].specific_h_j_kg,
                0.5 * expected,
            );
        }
        let work = accepted.viscous_heating_audit().unwrap();
        assert_eq!(work.substeps, 1);
        near(work.generated_heat_j, 0.5 * power.heat_w);
        near(work.generated_heat_j, work.supplied_wall_work_j);
        near(next.energy_j() - old.energy_j(), work.supplied_wall_work_j);
        near(thermal.external_heat_j, work.generated_heat_j);
        near(thermal.residual_j, 0.0);
        near(next.mass_kg(), old.mass_kg());
        assert!(observer.applied().is_none());
    }
}

#[test]
fn layered_graded_fluid_uses_its_own_half_cell_resistances_at_the_solid_wall() {
    let x = [1.0, 2.0];
    let y = [0.5, 1.0];
    let depth = 0.75;
    let mut domain = VoxelDomain::graded([x.to_vec(), vec![0.3, y[0], y[1]], vec![depth]]).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(7));
    domain.set(1, 0, 0, Voxel::Solid(7));
    let mu = [100.0, 200.0, 300.0, 400.0];
    // The native geometry ignores these solid transport slots. In particular,
    // they must not set the source-fluid Reynolds gate or a wall resistance.
    let viscosity = [1e-20, 1e20, mu[0], mu[1], mu[2], mu[3]];
    let (flow, power) = creeping_flow_with_viscous_heat(
        &domain,
        &fluid(),
        &viscosity,
        &driven_box(),
        controls(),
        &CancelGate::new(),
    )
    .unwrap();

    // Independent one-coordinate Stokes problem for the four fluid cells.
    // D=Aq²-2Bq+C; the half-cell wall paths are parallel, while the reciprocal
    // interior shear uses the four area-weighted series resistances.
    let lx = 0.5 * (x[0] + x[1]);
    let ly = 0.5 * (y[0] + y[1]);
    let resistance: f64 = (0..4).map(|c| x[c % 2] * y[c / 2] / mu[c]).sum();
    let edge_mu = (x[0] + x[1]) * (y[0] + y[1]) / resistance;
    let gamma_per_q = ((1.0 / y[0] + 1.0 / y[1]) / ly - (1.0 / x[0] + 1.0 / x[1]) / lx) / depth;
    let shear_a = edge_mu * lx * ly * depth * gamma_per_q * gamma_per_q;
    let normal_a: f64 = (0..4)
        .map(|c| 4.0 * mu[c] / (x[c % 2] * y[c / 2] * depth))
        .sum();
    let wall_x = [
        depth * (x[0] * mu[0] + x[1] * mu[1]) / y[0],
        depth * (x[0] * mu[2] + x[1] * mu[3]) / y[1],
    ];
    let wall_y = [
        depth * (y[0] * mu[0] + y[1] * mu[2]) / x[0],
        depth * (y[0] * mu[1] + y[1] * mu[3]) / x[1],
    ];
    let a = normal_a
        + shear_a
        + (0..2)
            .map(|i| wall_x[i] / (y[i] * depth).powi(2) + wall_y[i] / (x[i] * depth).powi(2))
            .sum::<f64>();
    let b = wall_x[1] / (y[1] * depth);
    let q = b / a;
    near(flow.field.flux_x(1, 1, 0), -q);
    near(flow.field.flux_x(1, 2, 0), q);
    near(flow.field.flux_y(0, 2, 0), q);
    near(flow.field.flux_y(1, 2, 0), -q);
    near(power.heat_w, wall_x[1] - b * q);
    near(power.supplied_wall_power_w, power.heat_w);
    assert_eq!(power.cell_heat_w[0], 0.0);
    assert_eq!(power.cell_heat_w[1], 0.0);
    for c in 0..4 {
        let i = c % 2;
        let j = c / 2;
        let normal = 4.0 * mu[c] * q * q / (x[i] * y[j] * depth);
        let shear = shear_a * q * q * (x[i] * y[j] / mu[c]) / resistance;
        let u = if j == 0 {
            -q / (y[j] * depth)
        } else {
            q / (y[j] * depth)
        };
        let wall = if j == 0 { 0.0 } else { 1.0 };
        let along_x = depth * x[i] * mu[c] / y[j] * (wall - u).powi(2);
        let along_y = depth * y[j] * mu[c] / x[i] * (q / (x[i] * depth)).powi(2);
        near(power.cell_heat_w[c + 2], normal + shear + along_x + along_y);
    }
    near(power.pressure_constraint_power_w, 0.0);
    near(power.residual_w, 0.0);
}

fn ring() -> VoxelDomain {
    let mut domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    domain.set(1, 1, 0, Voxel::Solid(3));
    domain
}

#[test]
fn circulation_around_a_solid_has_exact_three_fluid_corner_power_and_no_solid_heat() {
    let domain = ring();
    let mu = 100.0;
    let mut viscosity = vec![mu; 9];
    viscosity[4] = 1e-20;
    let (flow, power) = creeping_flow_with_viscous_heat(
        &domain,
        &fluid(),
        &viscosity,
        &driven_box(),
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    // Eight fluid cells around one fixed obstacle have one circulating flux q.
    // There are no four-fluid shear edges. The four corner normal strains give
    // 16mu*q²; stationary exterior walls give 12mu*q², internal walls 16mu*q²,
    // and the moving lid 4mu*(U-q)². Hence D=mu*(48q²-8Uq+4U²), q=U/12.
    let q: f64 = 1.0 / 12.0;
    for x in 1..3 {
        near(flow.field.flux_x(x, 0, 0), -q);
        near(flow.field.flux_x(x, 2, 0), q);
        assert_eq!(flow.field.flux_x(x, 1, 0), 0.0);
    }
    for y in 1..3 {
        near(flow.field.flux_y(0, y, 0), q);
        near(flow.field.flux_y(2, y, 0), -q);
        assert_eq!(flow.field.flux_y(1, y, 0), 0.0);
    }
    let expected = [
        8.0 * q * q,
        4.0 * q * q,
        8.0 * q * q,
        4.0 * q * q,
        0.0,
        4.0 * q * q,
        7.0 * q * q + (1.0 - q).powi(2),
        2.0 * q * q + 2.0 * (1.0 - q).powi(2),
        7.0 * q * q + (1.0 - q).powi(2),
    ];
    for (actual, expected) in power.cell_heat_w.iter().zip(expected) {
        near(*actual, mu * expected);
    }
    assert_eq!(power.cell_heat_w[4], 0.0);
    near(power.heat_w, 11.0 * mu / 3.0);
    near(power.wall_tangential_force_n[3][0], 11.0 * mu / 3.0);
    near(power.supplied_wall_power_w, power.heat_w);
    near(power.wall_power_w[3], power.heat_w);
    near(power.pressure_constraint_power_w, 0.0);
    near(power.residual_w, 0.0);
}

#[test]
fn fixed_obstacle_at_rest_has_no_heat_and_unrepresented_wall_motion_refuses() {
    let domain = ring();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let (flow, power) =
        creeping_flow_with_viscous_heat(&domain, &fluid(), &[100.0; 9], &config, controls(), &gate)
            .unwrap();
    for velocity in flow.velocity_m_s.iter().flatten() {
        assert_eq!(*velocity, 0.0);
    }
    assert_eq!(power.cell_heat_w, vec![0.0; 9]);
    assert_eq!(power.heat_w, 0.0);
    assert_eq!(power.supplied_wall_power_w, 0.0);
    assert_eq!(power.residual_w, 0.0);

    let unresolved = VoxelDomain::from_fn(3, 2, 1, 1.0, |p| {
        if p[0] > 1.0 && p[0] < 2.0 {
            Voxel::Fluid
        } else {
            Voxel::Solid(0)
        }
    })
    .unwrap();
    assert!(matches!(
        creeping_flow_with_viscous_heat(
            &unresolved,
            &fluid(),
            &[100.0; 6],
            &driven_box(),
            controls(),
            &gate,
        ),
        Err(ChtError::InvalidInput {
            field: "simple.viscous_heat",
            ..
        })
    ));
}

#[test]
fn body_and_darcy_work_share_the_obstacle_flow_without_heating_the_solid() {
    let domain = ring();
    let mu = 100.0;
    let beta = 200.0;
    let f = 100.0;
    let force = [
        [-f, f, 0.0],
        [-f, 0.0, 0.0],
        [-f, -f, 0.0],
        [0.0, f, 0.0],
        [0.0; 3],
        [0.0, -f, 0.0],
        [f, f, 0.0],
        [f, 0.0, 0.0],
        [f, -f, 0.0],
    ];
    let mut drag = vec![beta; 9];
    drag[4] = 0.0;
    let (flow, power) = creeping_flow_with_body_force_heat(
        &domain,
        &fluid(),
        &[mu; 9],
        &force,
        Some(&drag),
        &driven_box(),
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    // Eight driven MAC faces: body supply=8Fq and Darcy heat=8 beta q².
    // They join the same independently reduced molecular quadratic above.
    let q = (4.0 * mu + 8.0 * f) / (48.0 * mu + 8.0 * beta);
    near(flow.field.flux_x(1, 2, 0), q);
    near(power.supplied_body_power_w, 8.0 * f * q);
    near(power.supplied_wall_power_w, 4.0 * mu * (1.0 - q));
    near(power.molecular_heat_w, mu * (48.0 * q * q - 8.0 * q + 4.0));
    near(power.darcy_heat_w, 8.0 * beta * q * q);
    near(
        power.heat_w,
        power.supplied_body_power_w + power.supplied_wall_power_w,
    );
    near(power.cell_heat_w.iter().sum(), power.heat_w);
    assert_eq!(power.cell_heat_w[4], 0.0);
    assert_eq!(power.cell_body_power_w.as_ref().unwrap()[4], 0.0);
    assert_eq!(power.cell_darcy_heat_w.as_ref().unwrap()[4], 0.0);
    near(power.pressure_constraint_power_w, 0.0);
    near(power.residual_w, 0.0);
}
