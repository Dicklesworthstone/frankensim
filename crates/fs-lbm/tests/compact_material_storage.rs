//! Physical storage with resolved material fields, reservoirs and compact junctions.

use fs_exec::CancelGate;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    AdvectiveEnthalpyField, AffineSpecificEnthalpy, CellSink, ChtError, CompactComponent,
    ConvectionScheme, EnergyConfig, FlowFace, FlowField, FluidProperties, SolidMaterial,
    ThermalFace, ThermalSetup, Voxel, VoxelDomain,
};

fn fluid(reference_cp: f64) -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: reference_cp,
        conductivity_w_m_k: 17.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

fn component(power_w: f64) -> CompactComponent {
    CompactComponent {
        lo: [0, 0, 1],
        hi: [1, 1, 2],
        board_face: Face3::ZMin,
        power_w,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: Some(6.0),
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-8,
        "actual={actual:.15}, expected={expected:.15}"
    );
}

#[test]
fn resolved_compact_storage_uses_supplied_k_and_rejects_duplicate_capacity() {
    let domain = VoxelDomain::from_fn(1, 1, 3, 1.0, |_| Voxel::Solid(0)).unwrap();
    // Reference conductivity and absent volumetric capacity must not replace
    // the explicitly supplied field and per-cell capacities.
    let solids = [SolidMaterial::new("reference only", 123.0)];
    let flow = FlowField::quiescent(&domain);
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.compact_components.push(component(0.0));
    let conductivity = [[2.0; 3]; 3];
    let previous = [300.0, 400.0, 300.0];
    let config = EnergyConfig::default();

    let (solution, record) = config
        .step_with_conductivity(
            &domain,
            &fluid(1.0),
            &solids,
            &flow,
            &setup,
            &conductivity,
            &[2.0, 0.0, 2.0],
            &previous,
            3.0,
            &CancelGate::new(),
        )
        .unwrap();

    // Each branch: G=1/(1+0.5/2)=4/5 W/K. The two outside cells
    // reduce to C=4 J/K, coupled to Cj=6. Weighted equilibrium is 360 K.
    // dt*Gtotal*(1/6+1/4)=2, so the temperature difference falls by 1/3.
    close(solution.temperature[0], 340.0);
    close(solution.temperature[1], 1120.0 / 3.0);
    close(solution.temperature[2], 340.0);
    close(solution.junctions[0].temperature_k, 1120.0 / 3.0);
    close(solution.junctions[0].board_w, 80.0 / 3.0);
    close(solution.junctions[0].case_w, 80.0 / 3.0);
    close(
        2.0 * solution.temperature[0]
            + 6.0 * solution.temperature[1]
            + 2.0 * solution.temperature[2],
        3600.0,
    );
    close(record.stored_energy_change_j, 0.0);
    close(record.source_j, 0.0);
    close(record.boundary_outflow_j, 0.0);
    close(record.closure_j, 0.0);
    close(solution.report.balance.sink_outflow_w, 0.0);

    // The supplied cell slice cannot add any second capacity to an alias,
    // even one small enough to disappear in ordinary floating-point sums.
    for duplicate in [6.0, f64::MIN_POSITIVE] {
        let result = config.step_with_conductivity(
            &domain,
            &fluid(1.0),
            &solids,
            &flow,
            &setup,
            &conductivity,
            &[2.0, duplicate, 2.0],
            &previous,
            3.0,
            &CancelGate::new(),
        );
        assert!(
            matches!(result, Err(ChtError::InvalidInput { .. })),
            "duplicate alias capacity={duplicate}"
        );
    }
}

#[test]
fn resolved_fixed_reservoir_jump_separates_storage_from_physical_sinks() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("reference only", 123.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature.push((0, 320.0));
    setup.power_w = vec![0.0, 8.0];
    setup.cell_sinks = vec![
        CellSink {
            cell: 0,
            conductance_w_k: 0.5,
            temperature_k: 300.0,
        },
        CellSink {
            cell: 1,
            conductance_w_k: 1.0,
            temperature_k: 300.0,
        },
    ];

    let (solution, record) = EnergyConfig::default()
        .step_with_conductivity(
            &domain,
            &fluid(1.0),
            &solids,
            &FlowField::quiescent(&domain),
            &setup,
            &[[1.0; 3]; 2],
            &[3.0, 2.0],
            &[300.0, 300.0],
            1.0,
            &CancelGate::new(),
        )
        .unwrap();

    // Resolved k=1 gives G=1. For theta=Tfree-300:
    // 2*theta + (theta-20) + theta = 8, hence theta=7.
    // The fixed reservoir supplies its 60 J jump plus 23 J reaction.
    // Physical sinks remove 10+7=17 J; net external outflow is 17-83.
    close(solution.temperature[0], 320.0);
    close(solution.temperature[1], 307.0);
    let balance = solution.report.balance;
    close(balance.source_w, 8.0);
    close(balance.fixed_cell_injection_w, 23.0);
    close(balance.boundary_outflow_w, 0.0);
    close(balance.sink_outflow_w, 17.0);
    // Storage is excluded from the physical steady-form balance. Its
    // residual is the free cell's 2*(307-300)/dt storage rate.
    close(balance.residual_w, 14.0);

    close(record.stored_energy_change_j, 74.0);
    close(record.source_j, 8.0);
    close(record.boundary_outflow_j, -66.0);
    close(record.closure_j, 0.0);
}

#[test]
fn resolved_enthalpy_advection_and_compact_storage_share_one_heat_balance() {
    // Two independent one-cell fluid channels run in x on either side of
    // the middle z-layer junction. No fluid passes through the component.
    let domain = VoxelDomain::from_fn(1, 1, 3, 1.0, |p| {
        if p[2] > 1.0 && p[2] < 2.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let flow = FlowField::from_face_velocity(
        &domain,
        [
            FlowFace::Fixed,
            FlowFace::Fixed,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
        ],
        |axis, _| if axis == 0 { 1.0 } else { 0.0 },
    )
    .unwrap();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    let mut setup = ThermalSetup::new(faces);
    setup.compact_components.push(component(36.0));
    let solids = [SolidMaterial::new("collapsed package", 123.0)];
    let config = EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        ..EnergyConfig::default()
    };

    for datum in [0.0, 1024.0] {
        let caloric = AffineSpecificEnthalpy {
            temperature_k: 300.0,
            enthalpy_j_kg: datum,
            heat_capacity_j_kg_k: 2.0,
        };
        let cells = [Some(caloric), None, Some(caloric)];
        let enthalpy = AdvectiveEnthalpyField {
            cells: &cells,
            boundary_j_kg: [Some(datum), Some(datum), None, None, None, None],
        };
        let (solution, record) = config
            .step_with_enthalpy(
                &domain,
                &fluid(999.0),
                &solids,
                &flow,
                &setup,
                &[[1.0; 3]; 3],
                enthalpy,
                &[2.0, 0.0, 2.0],
                &[300.0; 3],
                1.0,
                &CancelGate::new(),
            )
            .unwrap();

        // Each channel has storage 2, inlet conduction 2, advective Cp=2,
        // and junction G=2/3. In temperature rises:
        // (2+2+2+2/3)*theta_f=(2/3)*theta_j, so theta_f=theta_j/10.
        // The junction equation gives (6+4/3)*theta_j-(4/3)*theta_f=36.
        close(solution.temperature[0], 300.5);
        close(solution.temperature[1], 305.0);
        close(solution.temperature[2], 300.5);
        close(solution.junctions[0].board_w, 3.0);
        close(solution.junctions[0].case_w, 3.0);
        close(solution.solid_to_fluid_heat_w(&domain), 6.0);

        let balance = solution.report.balance;
        close(balance.source_w, 36.0);
        close(balance.advective_outflow_w, 2.0);
        close(balance.boundary_outflow_w, 4.0);
        close(balance.sink_outflow_w, 0.0);
        close(balance.residual_w, 32.0);
        close(record.stored_energy_change_j, 32.0);
        close(record.source_j, 36.0);
        close(record.boundary_outflow_j, 4.0);
        close(record.closure_j, 0.0);
    }
}
