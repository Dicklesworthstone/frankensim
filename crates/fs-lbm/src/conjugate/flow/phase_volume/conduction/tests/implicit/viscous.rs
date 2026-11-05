//! Molecular heat must enter the implicit endpoint before its one native kick.

use super::*;
use crate::conjugate::FvBoundary;
use crate::conjugate::flow::{
    MaterialMomentumViscosityInput, MaterialViscosityControls, StaggeredBoundaryPressure,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::Cell;

struct CountedTransport {
    state: FluidTransportState,
    calls: Cell<usize>,
}

impl FluidTransport for CountedTransport {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.calls.set(self.calls.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn molecular_heat_enters_the_asymmetric_bath_endpoint_before_native_momentum() {
    let material = card("T", &[(300.0, 2.0), (304.0, 2.0)]);
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 2.0), (4.0, 304.0, 1.0, 2.0)]);
    let conductivity = source(&material, &caloric, [0.0, 4.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let states = [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let momentum = [[2.0, 0.0, 0.0]; 2];
    let force = [[0.0; 3]; 2];
    let molecular = CountedTransport {
        state: FluidTransportState::try_new(2.0, 1.0, 2.0, 1.0).unwrap(),
        calls: Cell::new(0),
    };
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&molecular); 2];
    let viscosity = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: [
            FvBoundary::Outlet,
            FvBoundary::Outlet,
            FvBoundary::wall(),
            FvBoundary::wall(),
            FvBoundary::Symmetry,
            FvBoundary::Symmetry,
        ],
    };
    let viscosity_controls = MaterialViscosityControls {
        max_cells: 100,
        max_sweeps: 200,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-8,
        maximum_angular_impulse_defect_kg_m2_s: 1e-8,
        maximum_energy_defect_j: 1e-8,
    };
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let profile = StaggeredBoundaryPressure {
        x_pa: &[Some(0.0), None, Some(0.0)],
        y_pa: &[None; 4],
        z_pa: &[None; 4],
    };
    // Side-wall molecular stress gives u_mu=4/5 and Q_mu=8/25 J/cell.
    // At dt=1/8, G_internal=G_bath=2, Q_authored=[1/8,-1/8],
    // the actual BE equations are 10*h0-h1=589/50 and -h0+9*h1=1239/50.
    // Their solution includes molecular heat BEFORE evaluating the one-sided
    // 302 K bath: h*=[654/445,12979/4450], Q_bath=59/445 J.
    let h_star = [654.0 / 445.0, 12979.0 / 4450.0];
    let mut volume = input(&mass, &states, &[0.125, -0.125]);
    volume.boundary_inflow[0] = Some(caloric.state_at_specific_enthalpy(h_star[0]).unwrap());
    let mut incoming = [None; 6];
    incoming[0] = Some([0.8, 0.0, 0.0]);
    let data = ThermalVolumeMomentumInput {
        volume,
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: incoming,
    };
    let controls = super::super::momentum::coupled_controls(0.125);
    let gate = CancelGate::new_clock_free();
    let result = field
        .advance_isotropic_implicit_conductive_thermal_volume_momentum_with_viscosity_and_baths(
            &domain,
            Some((faces, profile)),
            &conductivity,
            &[bath(0, Face3::XMin, 302.0, 0.25)],
            data,
            viscosity,
            controls,
            viscosity_controls,
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(),
            &gate,
        )
        .unwrap();
    assert_eq!(molecular.calls.get(), 2);
    let expected_heat = [8799.0 / 17800.0, -6439.0 / 17800.0];
    for cell in 0..2 {
        close(
            result.conductive_enthalpy_increment_j[cell],
            expected_heat[cell],
        );
        close(
            result.evolution.viscosity.cell_molecular_heat_j[cell],
            8.0 / 25.0,
        );
    }
    close(result.report.boundary_heat_j[0], 59.0 / 445.0);
    close(
        result.report.solve.energy.external_heat_j,
        16.0 / 25.0 + 59.0 / 445.0,
    );
    close(result.evolution.viscosity.molecular_heat_j, 16.0 / 25.0);
    close(result.evolution.viscosity.backward_euler_loss_j, 2.0 / 25.0);
    close(
        result
            .evolution
            .evolution
            .thermal
            .volume
            .report
            .requested_enthalpy_increment_j,
        16.0 / 25.0 + 59.0 / 445.0,
    );

    // Independent reference: the original viscous owner receives only authored
    // plus analytically computed conductive heat, then supplies its own Q_mu.
    // It must reproduce the same physical impulse and one mass/momentum route.
    let external = [0.125 + expected_heat[0], -0.125 + expected_heat[1]];
    let reference = field
        .advance_piecewise_thermal_volume_momentum_with_viscosity_and_pressure_profile(
            &domain,
            faces,
            profile,
            &caloric,
            ThermalVolumeMomentumInput {
                volume: CoexistenceVolumeInput {
                    net_enthalpy_increment_j: &external,
                    ..volume
                },
                ..data
            },
            viscosity,
            controls.evolution,
            viscosity_controls,
            &gate,
        )
        .unwrap();
    assert_eq!(molecular.calls.get(), 4);
    assert_eq!(result.evolution.viscosity, reference.viscosity);
    assert_eq!(
        result.evolution.viscosity_reconstruction,
        reference.viscosity_reconstruction
    );
    let transport = &result.evolution.evolution.thermal.volume.transport;
    // Inflow h equals h0*. Each accepted route is 1/5 kg; the native implicit
    // donor update gives h0=h0*, h1=(10*h1*+h0*)/11 after the heat interval.
    close(transport.specific_enthalpy_j_kg[0], h_star[0]);
    close(transport.specific_enthalpy_j_kg[1], 13633.0 / 4895.0);
    for cell in 0..2 {
        close(transport.mass_kg[cell], 2.0);
        close(
            transport.specific_enthalpy_j_kg[cell],
            reference
                .evolution
                .thermal
                .volume
                .transport
                .specific_enthalpy_j_kg[cell],
        );
        for axis in 0..3 {
            close(
                result.evolution.evolution.momentum.momentum_kg_m_s[cell][axis],
                reference.evolution.momentum.momentum_kg_m_s[cell][axis],
            );
        }
        close(
            result.evolution.evolution.momentum.velocity_m_s[cell][0],
            0.8,
        );
    }
    assert!(
        result.report.absolute_enthalpy_defect_j
            <= controls
                .evolution
                .volume
                .thermal
                .volume
                .transport
                .maximum_enthalpy_defect_j
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]; 2]);
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 1.0);
}
